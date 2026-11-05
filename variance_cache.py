"""0 대칭·끝점 포함 노드에서 내부 F_i 최댓값 분산을 n별 JSON으로 보관한다.

노드 형태 (m = (n - 2) // 2개의 매개변수 0 < p_1 < ... < p_m < 1):
    n 짝수: [-1, -p_m, ..., -p_1, p_1, ..., p_m, 1]
    n 홀수: [-1, -p_m, ..., -p_1, 0, p_1, ..., p_m, 1]
예) n=6 -> [-1, -b, -a, a, b, 1] (params = (a, b)), n=5 -> [-1, -a, 0, a, 1].

캐시 파일은 variance_cache/n{n}/variance.json이며 매개변수의 정확한 유리수 값을 키로 쓴다.
계산 함수의 소스가 바뀌면 identity가 달라져 이전 캐시는 무시된다.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import inspect
import json
import math
import os
from pathlib import Path
import tempfile
import time

import sympy as sp

import basis_cache
import generate_basis
import maximum

CACHE_DIR = Path(__file__).resolve().parent / 'variance_cache'
_FORMAT_VERSION = 1


def parameter_count(n):
    if n < 3:
        raise ValueError("At least three nodes are required")
    return (n - 2) // 2


def symmetric_nodes(n, params):
    """매개변수 (p_1 < ... < p_m)로 0 대칭·끝점 포함 노드를 만든다."""
    params = [sp.Rational(str(p)) if isinstance(p, float) else sp.sympify(p) for p in params]
    if len(params) != parameter_count(n):
        raise ValueError(f"n={n} needs {parameter_count(n)} parameters")
    middle = [sp.Integer(0)] if n % 2 else []
    return [sp.Integer(-1), *(-p for p in reversed(params)), *middle, *params, sp.Integer(1)]


def inner_variance(node):
    """양 끝 기저 F_0, F_n을 제외한 F_1, ..., F_{n-1} 최댓값의 모분산."""
    inner = maximum.maximums_of_F(node)[1:-1]
    mean = sum(inner) / len(inner)
    return sp.simplify(sum((value - mean) ** 2 for value in inner) / len(inner))


def _identity(n):
    functions = (maximum.maximums_of_F, inner_variance, symmetric_nodes,
                 generate_basis.lagrange_poly, generate_basis.lagrange_basis,
                 generate_basis.weights, generate_basis.Normalized_Lagrange,
                 generate_basis.Cumulative_Lagrange, generate_basis.New_Basis)
    return {'format_version': _FORMAT_VERSION, 'sympy_version': sp.__version__, 'n': n,
            'algorithm': basis_cache._digest([inspect.getsource(f) for f in functions])}


def _key(params):
    return ','.join(str(p) for p in params)


def _path(n):
    return CACHE_DIR / f'n{n}' / 'variance.json'


def load(n):
    """{매개변수 키: {'value': float | None, 'exact': str | None}}. 없거나 무효면 {}."""
    try:
        data = json.loads(_path(n).read_text(encoding='utf-8'))
        if data['identity'] != _identity(n) or data['checksum'] != basis_cache._digest(data['points']):
            return {}
        return data['points']
    except (OSError, ValueError, TypeError, KeyError, RecursionError):
        return {}


def save(n, points):
    """원자적으로 기록하며, 실패해도 계산 결과 반환은 막지 않는다."""
    temp_path = None
    try:
        path = _path(n)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {'identity': _identity(n), 'points': points,
                'checksum': basis_cache._digest(points)}
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix='.variance-', suffix='.tmp', delete=False) as out:
            temp_path = Path(out.name)
            json.dump(data, out, indent=1, allow_nan=False)
        os.replace(temp_path, path)
    except (OSError, ValueError, TypeError, RecursionError):
        pass
    finally:
        if temp_path is not None:
            try:
                temp_path.unlink(missing_ok=True)
            except OSError:
                pass


def _compute(n, params):
    """한 점의 분산. 가중치가 0인 등 정의되지 않으면 None."""
    try:
        exact = inner_variance(symmetric_nodes(n, params))
    except ValueError:
        return {'value': None, 'exact': None}
    return {'value': float(exact), 'exact': str(exact)}


def _compute_task(n, key):
    return key, _compute(n, [sp.Rational(p) for p in key.split(',')])


def get_variances(n, param_points, *, workers=None, save_every=200, progress=True):
    """매개변수 점들의 분산을 캐시에서 읽고, 없는 점만 계산해 캐시에 추가한다.

    param_points: 각 원소가 (p_1, ..., p_m)인 iterable. 0 < p_1 < ... < p_m < 1이
    아닌 점은 계산하지 않고 NaN을 돌려준다.
    workers: 병렬 프로세스 수 (None이면 CPU 수, 1이면 현재 프로세스에서 계산).
    반환값은 param_points 순서의 float 리스트 (정의되지 않으면 NaN).
    """
    m = parameter_count(n)
    params_list = []
    for params in param_points:
        params = tuple(sp.Rational(str(p)) if isinstance(p, float) else sp.Rational(p)
                       for p in params)
        if len(params) != m:
            raise ValueError(f"n={n} needs {m} parameters per point")
        valid = all(0 < p < 1 for p in params) and all(l < r for l, r in zip(params, params[1:]))
        params_list.append(params if valid else None)

    points = load(n)
    requested = {_key(p) for p in params_list if p is not None}
    missing = sorted(requested - points.keys())
    if missing:
        start, done = time.time(), 0

        def record(key, result):
            nonlocal done
            points[key] = result
            done += 1
            if done % save_every == 0 or done == len(missing):
                save(n, points)
                if progress:
                    elapsed = time.time() - start
                    eta = elapsed / done * (len(missing) - done)
                    print(f"[variance n={n}] {done}/{len(missing)} computed, "
                          f"{elapsed:.0f}s elapsed, ~{eta:.0f}s left", flush=True)

        if progress:
            print(f"[variance n={n}] {len(requested) - len(missing)} points from cache, "
                  f"computing {len(missing)}", flush=True)
        if workers is None:
            workers = os.cpu_count() or 1
        if workers <= 1 or len(missing) < 2:
            for key in missing:
                record(*_compute_task(n, key))
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_compute_task, n, key) for key in missing]
                for future in as_completed(futures):
                    record(*future.result())

    result = []
    for params in params_list:
        value = None if params is None else points[_key(params)]['value']
        result.append(math.nan if value is None else value)
    return result
