"""0 대칭 노드에서 매개변수에 따른 내부 F_i 최댓값 분산을 3차원으로 그린다.

양 끝 노드는 ±a(endpoint_variable)이고, 격자 매개변수(0 < p < 1)는 a에 대한 비율이다. 실제 노드는 a·p이다.
a(endpoint_variable)와 결과를 저장할 폴더(cache_dir)는 호출하는 쪽이 인자로 넘긴다.
관찰 폴더 endpoint_variable_{a}의 스크립트는 cache_dir로 그 폴더의 variance_cache/를 넘긴다.
"""
from fractions import Fraction
import gzip
import hashlib
import inspect
import json
import multiprocessing
import os
from pathlib import Path
import time

import numpy as np
import sympy as sp

import calculate_extremum

import calculate_variance
import s_compass_search as compass_search


def symmetric_nodes(n, params, endpoint_variable):
    """끝 노드 ±a 사이의 대칭 노드. params는 a에 대한 비율이므로 a = 1이면 params가 곧 노드이다."""
    nodes = [-1, *(-p for p in reversed(params)), *([0] if n % 2 else []), *params, 1]
    a = Fraction(endpoint_variable)
    return nodes if a == 1 else [a * v for v in nodes]


def plane_nodes(n, a, b, endpoint_variable):
    """(a, b) 평면 격자점의 대칭 노드.

    endpoint_variable이 있으면 양 끝을 ±endpoint_variable로 고정한 n=6, 7 노드
    [-1, -b, -a, (0), a, b, 1](a배)이고, None이면 b가 끝 노드인 n=4, 5 노드 [-b, -a, (0), a, b]이다.
    """
    if endpoint_variable is None:
        return [-b, -a, *([0] if n % 2 else []), a, b]
    return symmetric_nodes(n, (a, b), endpoint_variable)


def plane_last_column(N, endpoint_variable):
    """b = j/N의 마지막 j. 끝 노드를 고정하면 b < 1이고, b가 끝 노드이면 b = 1까지 계산한다."""
    return N if endpoint_variable is None else N - 1


def volume_nodes(n, params, endpoint_variable):
    """(a, b, c) 격자점의 대칭 노드.

    endpoint_variable이 있으면 양 끝을 ±endpoint_variable로 고정한 n=8, 9 노드
    [-1, -c, -b, -a, (0), a, b, c, 1](a배)이고, None이면 c가 끝 노드인 n=6, 7 노드 [-c, -b, -a, (0), a, b, c]이다.
    """
    if endpoint_variable is None:
        a, b, c = params
        return [-c, -b, -a, *([0] if n % 2 else []), a, b, c]
    return symmetric_nodes(n, params, endpoint_variable)


def volume_slab_count(last, i):
    """a = i/N 슬랩의 (b, c) 개수: i < j < k <= last."""
    return (last - i) * (last - i - 1) // 2


def volume_last(N, n, cache_dir):
    """저장된 입체 격자의 마지막 c 인덱스. 첫 슬랩의 길이로 c = 1까지 계산했는지 알아낸다."""
    with np.load(volume_slab_path(cache_dir, n, N, 1)) as data:
        count = len(data['V'])
    for last in (N, N - 1):
        if count == volume_slab_count(last, 1):
            return last
    raise ValueError(f"Unexpected slab size {count} for N={N}")


def plane_rows(N, n, cache_dir, kind="max"):
    """저장된 평면 격자의 행 i(a = i/N) 범위. 첫 행의 길이로 b = 1까지 계산했는지 알아낸다."""
    with np.load(volume_slab_path(cache_dir, n, N, 1, kind)) as data:
        last = 1 + len(data['V'])
    return range(1, last)


def to_fraction(value):
    """SymPy 유리수는 문자열을 거쳐, 그 밖의 수(float, int)는 그 값 그대로 Fraction으로 바꾼다."""
    return Fraction(str(value)) if isinstance(value, sp.Basic) else Fraction(value)


def has_negative_weight(n, params, endpoint_variable):
    """가중치 중 0 이하가 있는지. 가중치가 정확히 0이면 정의되지 않으므로 True로 본다."""
    try:
        _, weights = calculate_extremum.precise_extrema_and_weights(symmetric_nodes(n, [float(p) for p in params], endpoint_variable), kind="max")
    except ValueError:
        return True
    return min(weights) < 0


def load_points(n, cache_dir):
    path = Path(cache_dir) / f'n{n}' / 'variance.json'
    return json.loads(path.read_text(encoding='utf-8'))['points'] if path.exists() else {}


def save_points(n, points, cache_dir):
    path = Path(cache_dir) / f'n{n}' / 'variance.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'points': points}, indent=1), encoding='utf-8')


def cached_variance(n, params, points, endpoint_variable):
    """params의 분산을 points에서 읽고, 없으면 계산해 추가한다. 정의되지 않으면 NaN."""
    key = ",".join(str(p) for p in params)
    if key not in points:
        try:
            M, _ = calculate_extremum.exact_extrema_and_weights(symmetric_nodes(n, [to_fraction(p) for p in params], endpoint_variable), kind="max")
            points[key] = {'value': float(calculate_variance.exact_variance(M[1:-1])), 'exact': None}
        except ValueError:  # 가중치가 0인 점
            points[key] = {'value': None, 'exact': None}
    value = points[key]['value']
    return np.nan if value is None else value


def variance_curve(a_values, *, n=5, endpoint_variable, cache_dir):
    """노드 [-1, -a, 0, a, 1](n=5) 또는 [-1, -a, a, 1](n=4)에서 a에 따른 분산을 구한다.

    0 < a < 1이 아니거나 가중치가 0이 되어 정의되지 않는 점은 NaN으로 둔다.
    """
    if n not in (4, 5):
        raise ValueError("3D curve needs exactly one parameter (n = 4 or 5)")
    points = load_points(n, cache_dir)
    computed = len(points)

    V = np.full(len(a_values), np.nan)
    for p, a in enumerate(a_values):
        a = sp.Rational(str(a))
        if 0 < a < 1:
            V[p] = cached_variance(n, (a,), points, endpoint_variable)

    if len(points) != computed:
        save_points(n, points, cache_dir)
    return V


def variance_grid(a_values, b_values, *, n=6, endpoint_variable, cache_dir):
    """(a, b) 격자에서 양 끝을 제외한 분산을 variance_cache/n{n}/variance.json을 거쳐 구한다.

    n=6이면 노드 [-1, -b, -a, a, b, 1], n=7이면 [-1, -b, -a, 0, a, b, 1].
    a < b가 아니거나 가중치가 0이 되어 정의되지 않는 점은 NaN으로 둔다.
    """
    if n not in (6, 7):
        raise ValueError("3D surface needs exactly two parameters (n = 6 or 7)")
    points = load_points(n, cache_dir)
    computed = len(points)

    Z = np.full((len(a_values), len(b_values)), np.nan)
    for p, a in enumerate(a_values):
        for q, b in enumerate(b_values):
            a, b = sp.Rational(str(a)), sp.Rational(str(b))
            if 0 < a < b < 1:
                Z[p, q] = cached_variance(n, (a, b), points, endpoint_variable)

    if len(points) != computed:
        save_points(n, points, cache_dir)
    return Z


def grid_dir(cache_dir, n, N, kind="max"):
    """격자 결과 폴더. kind="max"(F_i 최댓값 M_i)는 n{n}/N{N}/, kind="min"(최솟값 m_i)은 n{n}/N{N}_min/이다."""
    if kind not in ("max", "min"):
        # ValueError는 격자 계산에서 정의되지 않는 점(NaN)으로 처리되므로 TypeError를 쓴다.
        raise TypeError(f'kind must be "max" or "min", not {kind!r}')
    return Path(cache_dir) / f'n{n}' / (f'N{N}' if kind == "max" else f'N{N}_{kind}')


def volume_slab_path(cache_dir, n, N, i, kind="max"):
    return grid_dir(cache_dir, n, N, kind) / f'a{i:04d}.npz'


def compute_volume_slab(task):
    """a=i/N로 고정한 모든 (b, c)의 결과를 Fraction으로 정확히 계산해 저장한다.

    끝 노드를 고정하면 c < 1, endpoint_variable=None이면 c가 끝 노드이므로 c = 1까지 계산한다.

    a{i}.npz: V(내부 분산), V_all(끝점 포함 분산), M(F_0~F_n의 최댓값), weights, positive.
    실수 값은 정확한 유리수를 가장 가까운 배정밀도로 반올림한 것이다.
    a{i}.exact.txt.gz: 내부 분산의 정확한 값 "분자/분모"를 한 줄에 하나씩 저장한다.
    (b, c)는 저장하지 않으며 np.triu_indices의 순서로 복원한다. 정의되지 않은 점은 NaN이다.
    """
    n, N, i, endpoint_variable, cache_dir = task
    path = volume_slab_path(cache_dir, n, N, i)
    last = plane_last_column(N, endpoint_variable)
    count = volume_slab_count(last, i)
    V = np.full(count, np.nan)
    V_all = np.full(count, np.nan)
    M_values = np.full((count, n + 1), np.nan)
    weight_values = np.full((count, n), np.nan)
    positive = np.zeros(count, dtype=bool)

    path.parent.mkdir(parents=True, exist_ok=True)
    exact_path = path.with_name(path.stem + '.exact.txt.gz')
    exact_temp = path.with_name(path.stem + '.exact.tmp')
    a = Fraction(i, N)
    position = 0
    with gzip.open(exact_temp, 'wt', encoding='ascii', newline='\n') as exact_file:
        for j in range(i + 1, last):
            for k in range(j + 1, last + 1):
                try:
                    M, weights = calculate_extremum.exact_extrema_and_weights(
                        volume_nodes(n, (a, Fraction(j, N), Fraction(k, N)), endpoint_variable), kind="max")
                except ValueError:
                    exact_file.write('nan\n')
                else:
                    inner = calculate_variance.exact_variance(M[1:-1])
                    V[position] = float(inner)
                    V_all[position] = float(calculate_variance.exact_variance(M))
                    M_values[position] = [float(value) for value in M]
                    weight_values[position] = [float(value) for value in weights]
                    positive[position] = min(weights) > 0
                    exact_file.write(f'{inner.numerator}/{inner.denominator}\n')
                position += 1

    # npz가 완료 표시이므로 마지막에 만든다. 중단되어도 깨진 파일이 완료로 남지 않는다.
    os.replace(exact_temp, exact_path)
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, V=V, V_all=V_all, M=M_values, weights=weight_values, positive=positive)
    os.replace(temp_path, path)
    return count


def variance_volume(N, *, n=8, endpoint_variable, cache_dir, workers=None):
    """0 < a < b < c < 1인 간격 1/N 격자의 분산을 병렬로 계산한다.

    n=8이면 노드 [-1, -c, -b, -a, a, b, c, 1], n=9이면 가운데에 0이 추가된다.
    endpoint_variable=None이면 c가 끝 노드인 0 < a < b < c <= 1,
    n=6이면 노드 [-c, -b, -a, a, b, c], n=7이면 가운데에 0이 추가된다.
    a 값 하나가 끝날 때마다 variance_cache/n{n}/N{N}/에 저장하므로,
    중단 후 다시 실행하면 끝난 a는 건너뛰고 이어서 계산한다.
    """
    if n not in ((6, 7) if endpoint_variable is None else (8, 9)):
        raise ValueError("Volume needs exactly three parameters (n = 8, 9, or n = 6, 7 with endpoint_variable=None)")
    last = plane_last_column(N, endpoint_variable)
    tasks = [(n, N, i, endpoint_variable, cache_dir) for i in range(1, last - 1)
             if not volume_slab_path(cache_dir, n, N, i).exists()]
    if not tasks:
        return

    total = sum(volume_slab_count(last, i) for _, _, i, _, _ in tasks)
    print(f"n={n}, 1/N={1 / N:g}: {len(tasks)} slabs, {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_volume_slab, tasks):
            done += count
            elapsed = time.perf_counter() - start
            remaining = elapsed * (total - done) / done if done else 0
            print(f"{done:,}/{total:,} points, elapsed {elapsed / 3600:.2f} h, "
                  f"remaining about {remaining / 3600:.2f} h", flush=True)


def load_variance_volume(N, *, n=8, cache_dir):
    """양수 가중치 격자점들의 (a, b, c 인덱스, 분산)과 격자 위의 국소 최솟값을 반환한다.

    국소 최솟값은 주변 26개 격자점 중 양수 가중치인 점들보다 분산이 작은 점이다.
    interior가 False이면 주변에 양수 가중치가 아니거나 a<b<c를 벗어나는 점이 있다는 뜻이다.
    """
    last = volume_last(N, n, cache_dir)

    def layer(i):
        # 인덱스 0과 N + 1은 항상 NaN이다 (c = 1까지 계산해도 N + 1은 비어 있다).
        grid = np.full((N + 2, N + 2), np.nan)
        if 1 <= i <= last - 2:
            with np.load(volume_slab_path(cache_dir, n, N, i)) as data:
                V, positive = data['V'], data['positive']
            j, k = np.triu_indices(last - i, 1)
            keep = positive & (V > 0)
            grid[j[keep] + i + 1, k[keep] + i + 1] = V[keep]
        return grid

    offsets = [(dj, dk) for dj in (-1, 0, 1) for dk in (-1, 0, 1)]
    points = {'a': [], 'b': [], 'c': [], 'V': []}
    minima = {'a': [], 'b': [], 'c': [], 'V': [], 'interior': []}
    previous, current, following = layer(0), layer(1), layer(2)
    for i in range(1, last - 1):
        finite = np.isfinite(current)
        j, k = np.nonzero(finite)
        points['a'].append(np.full(len(j), i, dtype=np.uint16))
        points['b'].append(j.astype(np.uint16))
        points['c'].append(k.astype(np.uint16))
        points['V'].append(current[finite])

        # 인덱스 0과 N + 1은 항상 NaN이므로 roll로 생기는 경계 값은 비교에 영향이 없다.
        is_minimum = finite.copy()
        interior = finite.copy()
        for depth, other in enumerate((previous, current, following)):
            for dj, dk in offsets:
                if depth == 1 and dj == 0 and dk == 0:
                    continue
                neighbor = np.roll(other, (dj, dk), axis=(0, 1))
                is_minimum &= ~(neighbor <= current)
                interior &= np.isfinite(neighbor)
        j, k = np.nonzero(is_minimum)
        minima['a'].append(np.full(len(j), i, dtype=np.uint16))
        minima['b'].append(j.astype(np.uint16))
        minima['c'].append(k.astype(np.uint16))
        minima['V'].append(current[is_minimum])
        minima['interior'].append(interior[is_minimum])

        previous, current, following = current, following, layer(i + 2)

    return ({name: np.concatenate(parts) for name, parts in points.items()},
            {name: np.concatenate(parts) for name, parts in minima.items()})


def compute_plane_row(task):
    """a=i/N로 고정한 모든 b의 결과를 Fraction으로 정확히 계산해 저장한다. 노드는 plane_nodes를 따른다.

    저장 형식은 compute_volume_slab과 같고, j번째 값은 b=(i+1+j)/N이다.
    kind="min"이면 M 자리에 최솟값 m_i를 저장하고 V, V_all은 최솟값의 분산이다.
    """
    n, N, i, endpoint_variable, cache_dir, kind = task
    path = volume_slab_path(cache_dir, n, N, i, kind)
    last = plane_last_column(N, endpoint_variable)
    count = last - i
    V = np.full(count, np.nan)
    V_all = np.full(count, np.nan)
    M_values = np.full((count, n + 1), np.nan)
    weight_values = np.full((count, n), np.nan)
    positive = np.zeros(count, dtype=bool)

    path.parent.mkdir(parents=True, exist_ok=True)
    exact_path = path.with_name(path.stem + '.exact.txt.gz')
    exact_temp = path.with_name(path.stem + '.exact.tmp')
    a = Fraction(i, N)
    with gzip.open(exact_temp, 'wt', encoding='ascii', newline='\n') as exact_file:
        for position, j in enumerate(range(i + 1, last + 1)):
            try:
                M, weights = calculate_extremum.exact_extrema_and_weights(
                    plane_nodes(n, a, Fraction(j, N), endpoint_variable), kind=kind)
            except ValueError:
                exact_file.write('nan\n')
                continue
            inner = calculate_variance.exact_variance(M[1:-1])
            V[position] = float(inner)
            V_all[position] = float(calculate_variance.exact_variance(M))
            M_values[position] = [float(value) for value in M]
            weight_values[position] = [float(value) for value in weights]
            positive[position] = min(weights) > 0
            exact_file.write(f'{inner.numerator}/{inner.denominator}\n')

    os.replace(exact_temp, exact_path)
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, V=V, V_all=V_all, M=M_values, weights=weight_values, positive=positive)
    os.replace(temp_path, path)
    return count


def variance_plane(N, *, n=6, endpoint_variable, cache_dir, workers=None, kind="max"):
    """간격 1/N 격자의 (a, b) 평면에서 분산을 Fraction으로 정확히 계산한다.

    endpoint_variable이 있으면 0 < a < b < 1, n=6이면 노드 [-1, -b, -a, a, b, 1](양 끝 ±endpoint_variable),
    n=7이면 가운데에 0이 추가된다.
    endpoint_variable=None이면 b가 끝 노드인 0 < a < b <= 1, n=4이면 노드 [-b, -a, a, b], n=5이면 가운데에 0.
    a 값 하나가 끝날 때마다 cache_dir/n{n}/N{N}/에 저장하므로 중단해도 이어서 계산한다.
    kind="min"이면 F_i 최솟값 m_i의 분산을 계산해 cache_dir/n{n}/N{N}_min/에 저장한다.
    """
    if n not in ((4, 5) if endpoint_variable is None else (6, 7)):
        raise ValueError("Plane needs exactly two parameters (n = 6, 7, or n = 4, 5 with endpoint_variable=None)")
    last = plane_last_column(N, endpoint_variable)
    tasks = [(n, N, i, endpoint_variable, cache_dir, kind) for i in range(1, last)
             if not volume_slab_path(cache_dir, n, N, i, kind).exists()]
    if not tasks:
        return

    total = sum(last - task[2] for task in tasks)
    print(f"n={n}, 1/N={1 / N:g}, {kind}: {len(tasks)} rows, {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_plane_row, tasks, chunksize=4):
            done += count
            elapsed = time.perf_counter() - start
            remaining = elapsed * (total - done) / done if done else 0
            print(f"\r{done:,}/{total:,} points, elapsed {elapsed / 60:.1f} min, "
                  f"remaining about {remaining / 60:.1f} min", end="", flush=True)
    print()


def active_row_path(cache_dir, n, N, i, kind="max"):
    return volume_slab_path(cache_dir, n, N, i, kind).with_suffix('.active.npz')


def compute_active_row(task):
    """a=i/N로 고정한 모든 b에서 내부 F_1~F_{n-1}의 최댓값(kind="min"이면 최솟값)을 주는 후보 이름을 정확히 구해 저장한다.

    label은 F_i별 후보 이름을 ","로 이은 문자열이며, 가중치가 0인 점은 빈 문자열이다.
    """
    n, N, i, endpoint_variable, cache_dir, kind = task
    path = active_row_path(cache_dir, n, N, i, kind)
    labels = []
    a = Fraction(i, N)
    for j in range(i + 1, plane_last_column(N, endpoint_variable) + 1):
        try:
            _, _, active = calculate_extremum.exact_extrema_and_weights(
                plane_nodes(n, a, Fraction(j, N), endpoint_variable), active=True, kind=kind)
        except ValueError:
            labels.append('')
        else:
            labels.append(",".join(active[1:-1]))
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, labels=np.array(labels))
    os.replace(temp_path, path)
    return len(labels)


def active_plane(N, *, n=6, endpoint_variable, cache_dir, workers=None, kind="max"):
    """variance_plane과 같은 격자에서 최댓값(kind="min"이면 최솟값) 후보 라벨을 계산한다. 중단해도 이어서 계산한다."""
    last = plane_last_column(N, endpoint_variable)
    tasks = [(n, N, i, endpoint_variable, cache_dir, kind) for i in range(1, last)
             if not active_row_path(cache_dir, n, N, i, kind).exists()]
    if not tasks:
        return
    total = sum(last - task[2] for task in tasks)
    print(f"active candidates n={n}, 1/N={1 / N:g}, {kind}: {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_active_row, tasks, chunksize=4):
            done += count
            elapsed = time.perf_counter() - start
            print(f"\r{done:,}/{total:,} points, elapsed {elapsed / 60:.1f} min, remaining about "
                  f"{elapsed * (total - done) / done / 60:.1f} min", end="", flush=True)
    print()


def load_active_plane(N, *, n=6, cache_dir, kind="max"):
    """(N+1)x(N+1) 정수 라벨 격자(-1은 정의되지 않음)와 라벨 이름 목록을 반환한다."""
    ids = np.full((N + 1, N + 1), -1, dtype=np.int32)
    names = {}
    for i in plane_rows(N, n, cache_dir, kind):
        with np.load(active_row_path(cache_dir, n, N, i, kind)) as data:
            row = data['labels']
        ids[i, i + 1:i + 1 + len(row)] = [names.setdefault(label, len(names)) if label else -1 for label in row]
    return ids, list(names)


def load_variance_plane(N, *, n=6, cache_dir, kind="max"):
    """(N+1)x(N+1) 격자 V[i, j](a=i/N, b=j/N)와 positive를 반환한다. 정의되지 않으면 NaN."""
    V = np.full((N + 1, N + 1), np.nan)
    positive = np.zeros((N + 1, N + 1), dtype=bool)
    for i in plane_rows(N, n, cache_dir, kind):
        with np.load(volume_slab_path(cache_dir, n, N, i, kind)) as data:
            end = i + 1 + len(data['V'])
            V[i, i + 1:end] = data['V']
            positive[i, i + 1:end] = data['positive']
    return V, positive


def variance_colormap(low, high, color_max=0.0, red_at=4.0):
    """log10(분산)용 고정 색: color_max까지 viridis, 그 위는 주황, red_at에서 빨강, 최댓값은 검붉은색."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    split = max(color_max, low + 1.0)
    top = max(high, split + 1.0)
    red = min((red_at - split) / (top - split), 0.9)
    viridis = plt.get_cmap("viridis")
    cmap = LinearSegmentedColormap.from_list("variance", [
        *((0.5 * t, viridis(t)) for t in np.linspace(0, 1, 32)),
        (0.5 + 0.25 * red, "darkorange"),
        (0.5 + 0.5 * red, "red"),
        (1.0, "#3d0000"),
    ], N=1024)
    return cmap, TwoSlopeNorm(vmin=low, vcenter=split, vmax=top)


def cached_compass_search(n, h, h_limit, initial_node, *, cache_dir):
    """compass_search.compass_search_ver3의 결과를 variance_cache/n{n}/compass.json을 거쳐 구한다.

    h, h_limit, 시작 노드와 탐색 코드가 모두 같을 때만 저장된 결과를 재사용한다.
    """
    source = "".join(inspect.getsource(function) for function in (
        compass_search.compass_search_ver3, calculate_variance.precise_variance_of_maximums_of_F,
        calculate_extremum.precise_extrema_and_weights))
    code = hashlib.sha256(source.encode('utf-8')).hexdigest()[:16]
    key = "|".join([code, repr(h), repr(h_limit), *(repr(float(v)) for v in initial_node)])

    path = Path(cache_dir) / f'n{n}' / 'compass.json'
    runs = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    if key in runs:
        return runs[key]

    node_history, min_history, h_history = compass_search.compass_search_ver3(n, h, h_limit, initial_node)
    runs[key] = [[[float(v) for v in node] for node in node_history],
                 [float(v) for v in min_history], [float(v) for v in h_history]]
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + '.tmp')
    temp_path.write_text(json.dumps(runs), encoding='utf-8')
    os.replace(temp_path, path)
    return runs[key]


def add_toggles(fig, rect, toggles):
    """toggles = {이름: (처음 상태, 상태를 받는 함수)}로 체크박스를 만든다. 비어 있으면 None.

    체크박스가 가비지 컬렉션되면 동작하지 않으므로 호출자가 반환값을 보관해야 한다.
    """
    from matplotlib.widgets import CheckButtons

    if not toggles:
        return None
    names = list(toggles)
    check = CheckButtons(fig.add_axes(rect), names, [toggles[name][0] for name in names])
    state = {name: toggles[name][0] for name in names}

    def on_clicked(name):
        state[name] = not state[name]
        toggles[name][1](state[name])
        fig.canvas.draw_idle()

    check.on_clicked(on_clicked)
    return check


def set_visible(artists, visible):
    for artist in artists:
        artist.set_visible(visible)


def add_view_sliders(fig, ax, elev, azim):
    from matplotlib.widgets import Slider

    # 마우스로 드래그해도 회전하며, 슬라이더로 정확한 각도를 지정할 수도 있다.
    elev_slider = Slider(fig.add_axes([0.15, 0.08, 0.7, 0.03]), "elev", -90, 90, valinit=elev)
    azim_slider = Slider(fig.add_axes([0.15, 0.03, 0.7, 0.03]), "azim", -180, 180, valinit=azim)

    def update(_):
        ax.view_init(elev=elev_slider.val, azim=azim_slider.val)
        fig.canvas.draw_idle()

    elev_slider.on_changed(update)
    azim_slider.on_changed(update)

    def sync_sliders(_):
        # 마우스로 돌린 뒤 슬라이더 값도 현재 각도에 맞춘다.
        for slider, value in ((elev_slider, ax.elev), (azim_slider, ax.azim)):
            if abs(slider.val - value) > 1e-9:
                slider.eventson = False
                slider.set_val(((value + 180) % 360) - 180 if slider is azim_slider else value)
                slider.eventson = True

    fig.canvas.mpl_connect("button_release_event", sync_sliders)
    return elev_slider, azim_slider


def draw_variance_curve(a_values, V, *, n=5, endpoint_variable, histories=None):
    """매개변수가 a 하나일 때 x=a, y=log10(분산)인 2차원 곡선을 그린다.

    음수 가중치가 있는 구간은 점선으로 표시한다.
    histories는 draw_variance_surface와 같고, 반환값은 (figure, axes)이다.
    """
    import matplotlib.pyplot as plt

    A = np.array([float(v) for v in a_values])
    V = np.asarray(V, dtype=float)
    mask = np.isfinite(V) & (V > 0)
    if not mask.any():
        raise ValueError("At least one positive finite variance is required")
    L = np.full_like(V, np.nan)
    L[mask] = np.log10(V[mask])

    negative = np.array([mask[p] and has_negative_weight(n, (A[p],), endpoint_variable) for p in range(len(A))])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(A, np.where(negative, np.nan, L), color="black", linewidth=1, label="variance")
    if negative.any():
        ax.plot(A, np.where(negative, L, np.nan), color="black", linewidth=1, linestyle="--",
                label="negative weight")
    for index, (name, history) in enumerate((histories or {}).items()):
        node_history, min_history, _ = history
        history_a = [node[-2] for node in node_history]
        ax.scatter(history_a, np.log10(np.array(min_history, dtype=float)),
                   color=f"C{index % 10}", s=4, linewidths=0, label=name)
    ax.set_xlabel("a")
    ax.set_ylabel("log10(variance)")
    ax.set_title(f"Variance of max F_i, n={n}")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", markerscale=2)
    return fig, ax


def draw_variance_surface(a_values, b_values, Z, *, n=6, endpoint_variable, histories=None, elev=25, azim=-130,
                          max_points=400, max_mesh=60):
    """마우스 드래그와 elev/azim 슬라이더로 돌려볼 수 있는 3차원 그래프.

    음수 가중치가 있는 격자점은 작은 검은 점으로 표시한다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 작은 점으로 표시한다.
    반환값은 (figure, axes, sliders)이다. 슬라이더가 가비지 컬렉션되면
    동작하지 않으므로 호출한 쪽에서 sliders를 변수로 들고 있어야 한다.
    """
    import matplotlib.pyplot as plt

    A, B = np.meshgrid([float(v) for v in a_values], [float(v) for v in b_values], indexing="ij")
    Z = np.asarray(Z, dtype=float)
    mask = np.isfinite(Z) & (Z > 0)
    if not mask.any():
        raise ValueError("At least one positive finite variance is required")
    # 분산이 0.003 ~ 50 정도로 퍼져 있어 log10 축으로 그린다.
    L = np.full_like(Z, np.nan)
    L[mask] = np.log10(Z[mask])

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_axes([0.05, 0.18, 0.9, 0.77], projection="3d", computed_zorder=False)
    # 드래그/슬라이더 회전 시 매 프레임 전체를 다시 그리므로, 촘촘한 격자에서는
    # (0.01 간격 약 4,700점에서 프레임당 0.75초) 화면이 멈춘 것처럼 보인다.
    # 점이 적으면 삼각분할 곡면과 격자점을, 많으면 최대 max_mesh x max_mesh로 줄인
    # 격자 곡면만 그려 회전을 가볍게 한다.
    if mask.sum() <= max_points:
        if mask.sum() >= 3:
            ax.plot_trisurf(A[mask], B[mask], L[mask], cmap="viridis", alpha=0.85, edgecolor="none")
        ax.scatter(A[mask], B[mask], L[mask], color="gray", s=10)
    else:
        # 정의되지 않은 점(a >= b)은 마스킹하고, 불투명 격자 곡면만 그린다.
        ax.plot_surface(A, B, np.ma.masked_invalid(L), cmap="viridis", edgecolor="none",
                        rcount=min(max_mesh, A.shape[0]), ccount=min(max_mesh, A.shape[1]),
                        antialiased=False)
    # 점을 격자 곡면과 같은 간격으로 솎아 회전이 느려지지 않게 한다.
    stride = max(1, int(np.ceil(max(A.shape) / max_mesh)))
    negative = np.zeros_like(mask)
    for p in range(0, A.shape[0], stride):
        for q in range(0, A.shape[1], stride):
            negative[p, q] = mask[p, q] and has_negative_weight(n, (A[p, q], B[p, q]), endpoint_variable)
    if negative.any():
        ax.scatter(A[negative], B[negative], L[negative], color="black", s=1, linewidths=0,
                   depthshade=False, label="negative weight")
    for index, (name, history) in enumerate((histories or {}).items()):
        node_history, min_history, _ = history
        history_a = [node[-3] for node in node_history]
        history_b = [node[-2] for node in node_history]
        ax.scatter(history_a, history_b, np.log10(np.array(min_history, dtype=float)),
                   color=f"C{index % 10}", s=1, linewidths=0, depthshade=False, label=name)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("log10(variance)")
    ax.set_title(f"Variance of max F_i, n={n}, a < b")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper left", markerscale=4)
    ax.view_init(elev=elev, azim=azim)

    return fig, ax, add_view_sliders(fig, ax, elev, azim)


def draw_variance_volume(N, *, n=8, cache_dir, max_points=3000, points_limit=30000, max_size=10,
                         color_max=0.0, red_at=4.0, zoom_limit=50, elev=25, azim=-130,
                         histories=None, show_histories=True, curve=None, show_curve=True):
    """(a, b, c) 위치에 양수 가중치 격자점을 찍고, 분산은 점의 크기와 색으로 나타낸다.

    log10(분산)이 작을수록 점이 크다. 슬라이더 값 이하인 점만 표시하며,
    회전이 느려지지 않도록 표시할 점은 max_points개 이내로 고르게 솎는다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 가는 선으로 표시한다.
    show_histories는 처음 상태이고, histories가 있으면 "history" 체크박스로 켜고 끌 수 있다.
    curve는 zero_variance_curve.zero_variance_curve_volume의 점 목록(분산 = 0 곡선)이며 검은 점으로 표시한다.
    show_curve는 처음 상태이고, curve가 있으면 "curve" 체크박스로 켜고 끌 수 있다.
    반환값은 (figure, axes, widgets)이다. widgets(슬라이더, 체크박스)가 가비지 컬렉션되면
    동작하지 않으므로 호출한 쪽에서 변수로 들고 있어야 한다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.widgets import Slider

    points, minima = load_variance_volume(N, n=n, cache_dir=cache_dir)
    if len(points['V']) == 0:
        raise ValueError("No grid point has all positive weights")
    # 점이 1억 개를 넘을 수 있으므로(1/N = 0.001, c <= 1) 그림용 배열은 float32로 줄인다.
    L = np.log10(points['V']).astype(np.float32)
    del points['V']
    low, high = float(L.min()), float(L.max())
    span = (high - low) or 1.0
    # 색은 슬라이더와 무관하게 log10(분산) 값으로만 정해진다.
    # 빨강 이상은 발산에 가까운 영역으로 본다.
    cmap, norm = variance_colormap(low, high, color_max, red_at)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.02, 0.39, 0.82, 0.57], projection="3d", computed_zorder=False)
    ax.set_proj_type("ortho")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("c")
    ax.view_init(elev=elev, azim=azim)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.88, 0.44, 0.02, 0.48]),
                 label="log10(variance)")

    for interior, marker, color, size, label in ((True, "*", "magenta", 25, "local min"),
                                                 (False, "x", "gray", 8, "local min on boundary")):
        chosen = minima['interior'] == interior
        if chosen.any():
            ax.scatter(minima['a'][chosen] / N, minima['b'][chosen] / N, minima['c'][chosen] / N,
                       marker=marker, color=color, s=size, depthshade=False, zorder=3, axlim_clip=True,
                       label=f"{label} ({int(chosen.sum())})")
    history_lines = []
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        history_lines += ax.plot(path[:, -4], path[:, -3], path[:, -2], color=plt.get_cmap("tab20")(index % 20),
                                 linewidth=0.7, marker="o", markersize=1, zorder=4, axlim_clip=True, label=name)
    curve_dots = []
    if curve:
        xyz = np.array([[float(Fraction(p['a'])), float(p['b_decimal']), float(p['c_decimal'])] for p in curve])
        curve_dots.append(ax.scatter(*xyz.T, s=3, color="black", depthshade=False, zorder=4, axlim_clip=True,
                                     label=f"variance = 0 curve ({len(curve)})"))
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper left", fontsize=7, ncols=2)
    set_visible(history_lines, show_histories)
    set_visible(curve_dots, show_curve)

    A, B, C = (points[name].astype(np.float32) / np.float32(N) for name in ('a', 'b', 'c'))
    del points
    # 줌 중심은 처음에 분산이 가장 작은 격자점이며, 슬라이더로 옮길 수 있다.
    best = int(np.argmin(L))
    scatter = [None]

    def update(_=None):
        if scatter[0] is not None:
            scatter[0].remove()
        width = 1 / zoom_slider.val
        center = np.array([slider.val for slider in center_sliders])
        lower = np.clip(center - width / 2, 0, 1 - width)
        upper = lower + width
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
        ax.set_zlim(lower[2], upper[2])

        # 줌을 해도 점이 드물어지지 않도록 화면 안의 점들 중에서 고른다.
        threshold = threshold_slider.val
        selected = np.nonzero((L <= threshold)
                              & (lower[0] <= A) & (A <= upper[0])
                              & (lower[1] <= B) & (B <= upper[1])
                              & (lower[2] <= C) & (C <= upper[2]))[0]
        total = len(selected)
        shown = int(points_slider.val)
        if total > shown:
            selected = selected[np.linspace(0, total - 1, shown, dtype=int)]
        scatter[0] = ax.scatter(A[selected], B[selected], C[selected], c=L[selected], cmap=cmap,
                                norm=norm, s=0.5 + (max_size - 0.5) * (high - L[selected]) / span,
                                linewidths=0, depthshade=False, zorder=2, axlim_clip=True)
        ax.set_title(f"n={n}, grid 1/{N}: {len(selected):,} of {total:,} points, "
                     f"log10(var) <= {threshold:.2f}, zoom x{zoom_slider.val:.1f}")
        fig.canvas.draw_idle()

    center_sliders = [Slider(fig.add_axes([0.15, y, 0.7, 0.025]), f"{name} center", 0, 1,
                             valinit=float(value[best]))
                      for name, y, value in (("a", 0.35, A), ("b", 0.305, B), ("c", 0.26, C))]
    zoom_slider = Slider(fig.add_axes([0.15, 0.215, 0.7, 0.025]), "zoom", 1, zoom_limit, valinit=1)
    points_slider = Slider(fig.add_axes([0.15, 0.17, 0.7, 0.025]), "points", 100,
                           max(points_limit, max_points), valinit=max_points, valstep=100)
    threshold_slider = Slider(fig.add_axes([0.15, 0.125, 0.7, 0.025]), "log10(var) <=", low, high,
                              valinit=float(np.percentile(L, 5)))
    for slider in (*center_sliders, zoom_slider, points_slider, threshold_slider):
        slider.on_changed(update)

    def on_scroll(event):
        if event.inaxes is ax:
            zoom_slider.set_val(float(np.clip(zoom_slider.val * 1.25 ** event.step, 1, zoom_limit)))

    fig.canvas.mpl_connect("scroll_event", on_scroll)
    update()

    # 탐색 경로와 분산 0 곡선 켜고 끄기 (있을 때만)
    toggles = {}
    if history_lines:
        toggles["history"] = (show_histories, lambda on: set_visible(history_lines, on))
    if curve_dots:
        toggles["curve"] = (show_curve, lambda on: set_visible(curve_dots, on))
    check = add_toggles(fig, [0.88, 0.385, 0.11, 0.045], toggles)

    return fig, ax, (*center_sliders, zoom_slider, points_slider, threshold_slider,
                     *add_view_sliders(fig, ax, elev, azim), *([check] if check else []))


def label_text(label):
    return "(" + ", ".join(part.replace("tau", "τ") for part in label.split(",")) + ")"


def draw_active_regions(N, *, n=6, cache_dir, histories=None, show_histories=True, show_boundaries=False,
                        near_zero=None, curve=None, kind="max"):
    """(a, b) 평면을 내부 F_1~F_{n-1}의 최댓값 후보 조합별로 색칠한 지도.

    색이 바뀌는 곳이 분산이 미분 불가능할 수 있는 경계이다. 검은 실선/점선/일점쇄선은
    가중치 ω_1, ω_2, ω_3 = 0인 곡선(분산이 발산하는 곳)이고, 음수 가중치 영역은 흐리게 칠한다.
    후보가 정확히 같은 값으로 겹친 점(라벨에 "/")은 빨간 원, 분산 최소 격자점은 검은 테두리의 자홍색 점이다.
    histories(탐색 경로)는 검은 선이며, show_histories는 처음 상태이고 "history" 체크박스로 켜고 끈다.
    show_boundaries=True이면 최댓값 후보가 바뀌는 경계(분산이 미분 불가능할 수 있는 곳)를 흰 선으로 그린다.
    near_zero가 있으면 양수 가중치이고 분산이 near_zero보다 작은 격자점을 자홍색 점으로 표시한다.
    curve가 있으면 (a, b) 점 목록(예: zero_variance_curve의 분산 0 곡선)을 검은 선으로 겹쳐 그린다.
    kind="min"이면 최솟값 후보 조합과 최솟값 분산(active_plane, variance_plane의 kind="min" 결과)으로 그린다.
    반환값은 (figure, axes)이고, 체크박스는 가비지 컬렉션되지 않도록 figure에 붙여 둔다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    ids, names = load_active_plane(N, n=n, cache_dir=cache_dir, kind=kind)
    V, positive = load_variance_plane(N, n=n, cache_dir=cache_dir, kind=kind)
    weights = np.full((N + 1, N + 1, n), np.nan)
    for i in plane_rows(N, n, cache_dir, kind):
        with np.load(volume_slab_path(cache_dir, n, N, i, kind)) as data:
            weights[i, i + 1:i + 1 + len(data['weights'])] = data['weights']

    order = sorted(range(len(names)), key=lambda k: -(ids == k).sum())
    rank = np.full(len(names), -1)
    rank[order] = np.arange(len(order))
    # 영역이 10개를 넘으면 tab10 색이 겹치므로 tab20을 쓴다.
    palette = "tab10" if len(order) <= 10 else "tab20"
    base = plt.get_cmap(palette)(np.arange(len(order)) % plt.get_cmap(palette).N)
    cmap = ListedColormap(base)
    cmap.set_bad("white")
    shown = np.ma.masked_less(np.where(ids >= 0, rank[ids], -1), 0).astype(float)

    fig, ax = plt.subplots(figsize=(11, 8))
    fig.subplots_adjust(left=0.07, right=0.62, top=0.93, bottom=0.08)
    extent = (-0.5 / N, 1 + 0.5 / N, -0.5 / N, 1 + 0.5 / N)
    ax.imshow(shown.T, origin="lower", extent=extent, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5,
              interpolation="nearest")
    # 음수 가중치 영역은 흰색을 덮어 흐리게 보이게 한다.
    fade = np.zeros((N + 1, N + 1, 4))
    fade[(ids >= 0) & ~positive] = (1, 1, 1, 0.55)
    ax.imshow(fade.transpose(1, 0, 2), origin="lower", extent=extent, interpolation="nearest")

    grid = np.arange(N + 1) / N
    for k, style in zip(range((n + 1) // 2), ("-", "--", "-.", ":")):
        W = np.ma.masked_invalid(weights[:, :, k])
        if np.nanmin(W) < 0 < np.nanmax(W):
            ax.contour(grid, grid, W.T, levels=[0], colors="black", linewidths=1, linestyles=style)
    ties = np.array([("/" in name) for name in names])
    tie_i, tie_j = np.nonzero((ids >= 0) & ties[np.maximum(ids, 0)])
    if len(tie_i):
        ax.scatter(tie_i / N, tie_j / N, s=60, facecolors="none", edgecolors="red", linewidths=1.2, zorder=5)
    if show_boundaries:
        segments = boundary_segments(ids, N)
        if segments:
            from matplotlib.collections import LineCollection
            ax.add_collection(LineCollection(segments, colors="white", linewidths=0.9, zorder=4))
    near = np.zeros_like(positive)
    if near_zero is not None:
        near = positive & np.isfinite(V) & (V < near_zero)
        near_i, near_j = np.nonzero(near)
        ax.scatter(near_i / N, near_j / N, s=4, color="magenta", linewidths=0, zorder=5)
    if curve is not None and len(curve):
        curve = np.asarray(curve, dtype=float)
        # 가지가 여러 개인 곡선(같은 a에 b가 둘)도 이어 그리지 않도록 점으로 그린다.
        ax.scatter(curve[:, 0], curve[:, 1], s=2, color="black", zorder=6)
    P = np.where(positive & np.isfinite(V) & (V > 0), V, np.inf)
    best = np.unravel_index(np.argmin(P), P.shape)
    ax.scatter([best[0] / N], [best[1] / N], marker="o", s=12, color="magenta", edgecolors="black",
               linewidths=0.5, zorder=6)
    history_lines = []
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        history_lines += ax.plot(path[:, -3], path[:, -2], color="black", linewidth=0.8, zorder=4)
    set_visible(history_lines, show_histories)
    fig.history_toggle = add_toggles(fig, [0.64, 0.02, 0.14, 0.05],
                                     {"history": (show_histories, lambda on: set_visible(history_lines, on))}
                                     if history_lines else {})

    handles = [Patch(color=base[r], label=f"{label_text(names[k])}  {int((ids == k).sum()):,}")
               for r, k in enumerate(order)]
    handles += [plt.Line2D([], [], color="black", linestyle=style, label=f"ω_{k + 1} = 0" if 2 * k + 1 == n else f"ω_{k + 1} = ω_{n - k} = 0")
                for k, style in zip(range((n + 1) // 2), ("-", "--", "-.", ":"))]
    handles += [plt.Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor="red", linestyle="",
                           label=f"exact tie ({len(tie_i)})"),
                plt.Line2D([], [], marker="o", color="magenta", markeredgecolor="black", linestyle="", markersize=4,
                           label="min variance")]
    if show_boundaries:
        from matplotlib import patheffects
        handles.append(plt.Line2D([], [], color="white", linewidth=1.5, label=f"arg{kind} changes (maybe not differentiable)",
                                  path_effects=[patheffects.withStroke(linewidth=3, foreground="0.4")]))
    if near_zero is not None:
        handles.append(plt.Line2D([], [], marker="o", color="magenta", linestyle="", markersize=3,
                                  label=f"variance < {near_zero:g} ({int(near.sum()):,})"))
    if curve is not None and len(curve):
        handles.append(plt.Line2D([], [], color="black", linewidth=1.0, label=f"variance = 0 curve ({len(curve)})"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8,
              title=f"arg{kind} of (F_1, ..., F_{n - 1})", title_fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_title(f"n={n}, grid 1/{N}: active candidates of the {kind}imum (faded = negative weight)")
    if show_boundaries:
        ax.set_facecolor("0.85")  # 흰 경계선이 흰 배경(정의되지 않은 영역)과 구별되도록
    return fig, ax


def boundary_segments(ids, N):
    """이웃한 두 격자점의 라벨(최댓값 후보)이 다르면 그 사이의 칸 경계를 선분으로 돌려준다."""
    half = 0.5 / N
    segments = []
    for axis in (0, 1):
        first = ids[:-1, :] if axis == 0 else ids[:, :-1]
        second = ids[1:, :] if axis == 0 else ids[:, 1:]
        for p, q in zip(*np.nonzero((first != second) & (first >= 0) & (second >= 0))):
            a0, b0 = p / N, q / N
            if axis == 0:   # (p, q)와 (p+1, q) 사이: a = (p + 1/2)/N에 세로 선분
                segments.append([(a0 + half, b0 - half), (a0 + half, b0 + half)])
            else:           # (p, q)와 (p, q+1) 사이: b = (q + 1/2)/N에 가로 선분
                segments.append([(a0 - half, b0 + half), (a0 + half, b0 + half)])
    return segments


def draw_variance_plane(N, *, n=6, cache_dir, histories=None, max_mesh=100, zoom_limit=50,
                        color_max=0.0, red_at=4.0, elev=25, azim=-130, show_boundaries=True,
                        show_histories=True, negative_color="black", near_zero=None, kind="max"):
    """간격 1/N 격자의 log10(분산) 곡면을 그리고, zoom 슬라이더로 최솟값 근처를 확대한다.

    회전이 느려지지 않도록 화면 안의 격자를 최대 max_mesh x max_mesh로 솎아 그리며,
    확대할수록 솎는 간격이 줄어 원래 1/N 격자가 그대로 보인다.
    색은 log10(분산) 값으로만 정해진다. 음수 가중치 점은 작은 점(색은 negative_color),
    near_zero가 있으면 양수 가중치이고 분산이 near_zero보다 작은 점은 자홍색 점,
    양수 가중치 점들 사이의 격자 국소 최솟값은 작은 자홍색 점으로 표시한다.
    active_plane 결과가 있으면, 최댓값 후보가 바뀌는 경계(분산이 미분 불가능할 수 있는 곳)를 흰 선으로 그린다.
    kind="min"이면 최솟값 m_i의 분산과 최솟값 후보 경계(variance_plane, active_plane의 kind="min" 결과)를 그린다.
    show_boundaries는 처음 상태이고, 그림 오른쪽 위의 "boundaries" 체크박스로 실행 중에 켜고 끌 수 있다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 가는 선으로 표시한다.
    show_histories는 처음 상태이고, histories가 있으면 "history" 체크박스로 켜고 끌 수 있다.
    반환값은 (figure, axes, widgets)이다. widgets(슬라이더, 체크박스)가 가비지 컬렉션되면
    동작하지 않으므로 호출자가 보관해야 한다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.widgets import Slider

    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    V, positive = load_variance_plane(N, n=n, cache_dir=cache_dir, kind=kind)
    ids = (load_active_plane(N, n=n, cache_dir=cache_dir, kind=kind)[0]
           if active_row_path(cache_dir, n, N, 1, kind).exists() else None)
    defined = np.isfinite(V) & (V > 0)
    L = np.full_like(V, np.nan)
    L[defined] = np.log10(V[defined])
    good = defined & positive
    if not good.any():
        raise ValueError("No grid point has all positive weights")
    low, high = float(np.nanmin(L)), float(np.nanmax(L))
    cmap, norm = variance_colormap(low, high, color_max, red_at)

    # 양수 가중치 점들 사이의 격자 국소 최솟값 (주변 8점과 비교)
    P = np.where(good, L, np.nan)
    is_minimum = good.copy()
    interior = good.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                neighbor = np.roll(P, (di, dj), axis=(0, 1))
                is_minimum &= ~(neighbor <= P)
                interior &= np.isfinite(neighbor)
    is_minimum &= interior
    best = np.unravel_index(np.nanargmin(P), P.shape)

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_axes([0.02, 0.28, 0.82, 0.63], projection="3d", computed_zorder=False)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("log10(variance)")
    ax.view_init(elev=elev, azim=azim)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.88, 0.33, 0.02, 0.52]),
                 label="log10(variance)")

    i, j = np.nonzero(is_minimum)
    ax.scatter(i / N, j / N, L[i, j], marker="o", color="magenta", s=4, linewidths=0, depthshade=False,
               zorder=5, axlim_clip=True, label=f"grid local min ({len(i)})")
    history_lines = []
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        history_lines += ax.plot(path[:, -3], path[:, -2], np.log10(np.array(history[1], dtype=float)),
                                 color=plt.get_cmap("tab10")(index % 10), linewidth=0.8, marker="o",
                                 markersize=1, zorder=4, axlim_clip=True, label=name)
    ax.legend(loc="upper left", fontsize=7)
    set_visible(history_lines, show_histories)

    drawn = []
    boundaries = {'on': show_boundaries and ids is not None}

    def update(_=None):
        while drawn:
            drawn.pop().remove()
        width = 1 / zoom_slider.val
        center = np.array([slider.val for slider in center_sliders])
        lower = np.clip(center - width / 2, 0, 1 - width)
        first = np.floor(lower * N).astype(int)
        last = np.minimum(np.ceil((lower + width) * N).astype(int), N)
        stride = max(1, int(np.ceil((last - first).max() / max_mesh)))
        rows = np.arange(first[0], last[0] + 1, stride)
        cols = np.arange(first[1], last[1] + 1, stride)
        window = L[np.ix_(rows, cols)]
        A, B = np.meshgrid(rows / N, cols / N, indexing="ij")
        drawn.append(ax.plot_surface(A, B, np.ma.masked_invalid(window), cmap=cmap, norm=norm,
                                     rstride=1, cstride=1, edgecolor="none", antialiased=False,
                                     zorder=1))
        negative = np.isfinite(window) & ~positive[np.ix_(rows, cols)]
        dots = max(1, int(np.ceil(max(window.shape) / 60)))
        thin = np.zeros_like(negative)
        thin[::dots, ::dots] = True
        negative &= thin
        if negative.any():
            drawn.append(ax.scatter(A[negative], B[negative], window[negative], color=negative_color,
                                    s=1 if negative_color == "black" else 3,
                                    linewidths=0, depthshade=False, zorder=2, axlim_clip=True))
        if near_zero is not None:
            near = (np.isfinite(window) & positive[np.ix_(rows, cols)]
                    & (window < np.log10(near_zero)))
            if near.any():
                drawn.append(ax.scatter(A[near], B[near], window[near], color="magenta", s=4,
                                        linewidths=0, depthshade=False, zorder=4, axlim_clip=True))
        if boundaries['on']:
            # 이웃한 두 격자점의 최댓값 후보가 다르면 그 사이에 짧은 선분을 곡면 높이로 그린다.
            label = ids[np.ix_(rows, cols)]
            half = stride / (2 * N)
            segments = []
            for axis in (0, 1):
                first_part = (slice(None, -1), slice(None)) if axis == 0 else (slice(None), slice(None, -1))
                second_part = (slice(1, None), slice(None)) if axis == 0 else (slice(None), slice(1, None))
                change = ((label[first_part] != label[second_part]) & (label[first_part] >= 0)
                          & (label[second_part] >= 0)
                          & np.isfinite(window[first_part]) & np.isfinite(window[second_part]))
                for p, q in zip(*np.nonzero(change)):
                    a0, b0 = A[p, q], B[p, q]
                    z = (window[p, q] + window[p + (axis == 0), q + (axis == 1)]) / 2
                    if axis == 0:
                        segments.append([(a0 + half, b0 - half, z), (a0 + half, b0 + half, z)])
                    else:
                        segments.append([(a0 - half, b0 + half, z), (a0 + half, b0 + half, z)])
            if segments:
                lines = Line3DCollection(segments, colors="white", linewidths=1.2, zorder=3, axlim_clip=True)
                ax.add_collection3d(lines)
                drawn.append(lines)
        ax.set_xlim(rows[0] / N, rows[-1] / N)
        ax.set_ylim(cols[0] / N, cols[-1] / N)
        if np.isfinite(window).any():
            ax.set_zlim(np.nanmin(window), np.nanmax(window))
        ax.set_title(f"n={n}, {kind}imum variance, grid 1/{N}, zoom x{zoom_slider.val:.1f}: "
                     f"{len(rows)}x{len(cols)} mesh, step {stride}/{N}\n"
                     f"min log10(var) = {L[best]:.3f} at a={best[0] / N:g}, b={best[1] / N:g}")
        fig.canvas.draw_idle()

    # 줌 중심은 처음에 분산이 가장 작은 격자점이며, 슬라이더로 옮길 수 있다.
    center_sliders = [Slider(fig.add_axes([0.15, y, 0.7, 0.03]), f"{name} center", 0, 1,
                             valinit=best[axis] / N)
                      for axis, (name, y) in enumerate((("a", 0.23), ("b", 0.18)))]
    zoom_slider = Slider(fig.add_axes([0.15, 0.13, 0.7, 0.03]), "zoom", 1, zoom_limit, valinit=1)
    for slider in (*center_sliders, zoom_slider):
        slider.on_changed(update)

    def on_scroll(event):
        if event.inaxes is ax:
            zoom_slider.set_val(float(np.clip(zoom_slider.val * 1.25 ** event.step, 1, zoom_limit)))

    fig.canvas.mpl_connect("scroll_event", on_scroll)

    # 후보 경계(흰 선)와 탐색 경로 켜고 끄기. 그릴 것이 없는 항목은 체크박스에 넣지 않는다.
    def set_boundaries(on):
        boundaries['on'] = on
        update()

    toggles = {}
    if ids is not None:
        toggles["boundaries"] = (boundaries['on'], set_boundaries)
    if history_lines:
        toggles["history"] = (show_histories, lambda on: set_visible(history_lines, on))
    check = add_toggles(fig, [0.84, 0.95 - 0.045 * len(toggles), 0.14, 0.045 * len(toggles)], toggles)
    update()

    return fig, ax, (*center_sliders, zoom_slider, *add_view_sliders(fig, ax, elev, azim),
                     *([check] if check else []))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    a_values = [sp.Rational(k, 10) for k in range(1, 9)]   # 0.1 ~ 0.8
    b_values = [sp.Rational(k, 10) for k in range(2, 10)]  # 0.2 ~ 0.9
    Z = variance_grid(a_values, b_values)

    print("a\\b " + " ".join(f"{float(b):>9.1f}" for b in b_values))
    for a, row in zip(a_values, Z):
        print(f"{float(a):.1f} " + " ".join("        -" if np.isnan(v) else f"{v:9.5f}" for v in row))

    fig, ax, sliders = draw_variance_surface(a_values, b_values, Z)
    fig.savefig("variance_surface.png", dpi=130)
    plt.show()
