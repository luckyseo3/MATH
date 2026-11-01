import inspect
import re
import unicodedata
from pathlib import Path
from functools import lru_cache

import numpy as np
import sympy as sp

def get_basis_maxima(BasisArr):
    """각 기저 집합의 [-1, 1] 구간 최댓값 배열을 리스트로 반환한다.

    BasisArr는 [F1, F2, ...]이며 각 F는 x에 대한 SymPy 식의 배열이다.
    다른 매개변수는 미리 대입해야 한다. 반환값의 i번째 NumPy 배열은
    BasisArr[i]에 대응하며, 기저 집합마다 함수 개수가 달라도 된다.
    """
    return [_get_single_basis_maxima(F) for F in BasisArr]


def get_basis_minima(BasisArr):
    """각 기저 집합의 [-1, 1] 구간 최솟값 배열을 리스트로 반환한다."""
    return [_get_single_basis_minima(F) for F in BasisArr]


def _get_single_basis_maxima(F):
    """하나의 기저 집합에 대한 최댓값을 계산한다."""
    return np.array([_get_expression_maximum(sp.sympify(F_i)) for F_i in F], dtype=float)
    

def _get_single_basis_minima(F):
    """하나의 기저 집합에 대한 최솟값을 계산한다."""
    return np.array([_get_expression_minimum(sp.sympify(F_i)) for F_i in F], dtype=float)


@lru_cache(maxsize=512)
def _get_expression_maximum(expression):
    """불변 SymPy 식별로 최댓값을 캐시해 보고서 간 중복 계산을 피한다.

    식이 바뀌면 새로 계산하며, 최대 512개를 실행 중 메모리에만 보관한다.
    """
    x = sp.symbols('x')
    if expression.free_symbols - {x}:
        raise ValueError("Substitute all parameters other than x before evaluation.")
    return float(sp.maximum(expression, x, sp.Interval(-1, 1)))


@lru_cache(maxsize=512)
def _get_expression_minimum(expression):
    """불변 SymPy 식별로 최솟값을 캐시해 보고서 간 중복 계산을 피한다."""
    x = sp.symbols('x')
    if expression.free_symbols - {x}:
        raise ValueError("Substitute all parameters other than x before evaluation.")
    return float(sp.minimum(expression, x, sp.Interval(-1, 1)))


def _validate_maxima(maxima):
    """지표 계산에 사용하는 비어 있지 않은 유한한 1차원 실수 배열을 반환한다."""
    values = np.asarray(maxima, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("Maxima must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(values)):
        raise ValueError("Maxima must contain only finite values.")
    return values


def get_maxima_range(maxima):
    """최댓값 배열의 범위(max - min)를 반환한다.

    각주: 가장 높은 기저함수 봉우리와 가장 낮은 봉우리의 높이 차이다.
    0이면 모든 봉우리의 높이가 같으며, 원래 기저함수 값과 단위가 같다.
    """
    values = _validate_maxima(maxima)
    return float(np.ptp(values))


def get_maxima_ratio(maxima):
    """최댓값 배열의 비율(max / min)을 반환한다. 분모가 0이면 예외를 발생시킨다.

    각주: 봉우리 높이가 모두 양수일 때 가장 높은 봉우리가 가장 낮은
    봉우리의 몇 배인지 나타낸다. 1이면 높이가 같고, 단위에 무관하다.
    음수 높이가 포함되면 이러한 배수 해석은 적용되지 않는다.
    """
    values = _validate_maxima(maxima)
    minimum = float(np.min(values))
    if minimum == 0:
        raise ValueError("Maxima ratio is undefined when the minimum is zero.")
    return float(np.max(values) / minimum)


def get_maxima_mean(maxima):
    """최댓값 배열의 산술평균을 반환한다.

    각주: 각 기저함수 봉우리의 평균 높이로, 기저 전체의 대표적인
    최고 높이를 나타낸다. 곡선 전체의 평균값이나 면적을 뜻하지 않는다.
    """
    return float(np.mean(_validate_maxima(maxima)))


def get_maxima_std(maxima):
    """최댓값 배열의 모집단 표준편차(ddof=0)를 반환한다.

    각주: 봉우리 높이가 평균 높이에서 얼마나 흩어져 있는지 나타낸다.
    0이면 높이가 모두 같으며, 원래 기저함수 값과 단위가 같다.
    """
    return float(np.std(_validate_maxima(maxima), ddof=0))


def get_maxima_cv(maxima):
    """변동계수(모집단 표준편차 / 평균)를 반환한다. 백분율이 아닌 비율이다.

    각주: 평균 높이가 양수일 때 봉우리 높이의 상대적인 불균일성을
    나타내며, 모든 기저를 같은 양수 배율로 확대해도 변하지 않는다.
    평균이 0이면 정의되지 않고, 음수이면 통상적인 상대 산포 해석에 부적합하다.
    """
    mean = get_maxima_mean(maxima)
    if mean == 0:
        raise ValueError("Coefficient of variation is undefined when the mean is zero.")
    return get_maxima_std(maxima) / mean


def _resolve_basis_names(BasisArr, basis_names):
    """보고서 호출자 영역에서 이름을 찾거나 명시한 이름을 검증한다."""
    if basis_names is None:
        frame = inspect.currentframe()
        caller = None
        try:
            caller = frame.f_back.f_back if frame is not None else None
            namespace = {**caller.f_globals, **caller.f_locals} if caller else {}
            basis_names = []
            for index, F in enumerate(BasisArr):
                matches = [name for name, value in namespace.items() if value is F]
                basis_names.append(matches[0] if len(matches) == 1 else f"basis_{index + 1}")
        finally:
            del caller
            del frame
    else:
        if isinstance(basis_names, str):
            raise ValueError("Basis names must be a sequence of names, not a string.")
        basis_names = list(basis_names)
    if len(basis_names) != len(BasisArr):
        raise ValueError("Provide one name per basis set.")
    if any(not isinstance(name, str) or not name.strip()
           or any(char in name for char in "\r\n\t") for name in basis_names):
        raise ValueError("Basis names must be non-empty single-line strings without tabs.")
    if len(set(basis_names)) != len(basis_names):
        raise ValueError("Basis names must be unique; specify basis_names explicitly.")
    return basis_names


_METRICS = (
    ("range", get_maxima_range),
    ("ratio", get_maxima_ratio),
    ("mean", get_maxima_mean),
    ("std", get_maxima_std),
    ("cv", get_maxima_cv),
)


def _evaluate_maxima(maxima):
    """원래 정밀도로 지표를 계산하고 출력할 때만 소수점 네 자리로 표시한다."""
    return {name: f"{evaluate(maxima):.4f}" for name, evaluate in _METRICS}


def save_basis_report(BasisArr, file_path="basis_report.txt", basis_names=None, *, overwrite_existing=False):
    """여러 기저의 상세 보고서를 UTF-8 TXT에 추가하고 절대 경로를 반환한다.

    BasisArr는 [F1, F2, ...]이다. 이름은 각 원소의 변수명에서 찾으며,
    찾을 수 없으면 basis_1, basis_2 등을 사용한다. basis_names로 명시할 수 있다.
    같은 이름이 이미 기록되어 있으면 그 원소의 계산과 저장을 생략한다.
    overwrite_existing=True이면 같은 이름의 기록을 현재 계산 결과로 교체한다.
    입력에 없는 다른 이름의 기록은 유지한다.
    파일이 없으면 생성한다. 지표 계산 실패 시 파일은 변경하지 않는다.
    """
    BasisArr = list(BasisArr)
    names = _resolve_basis_names(BasisArr, basis_names)
    path = Path(file_path).resolve()
    existing_text = path.read_text(encoding="utf-8-sig") if path.exists() else ""
    existing_lines = set(existing_text.splitlines())
    pending = [(name, F) for name, F in zip(names, BasisArr)
               if overwrite_existing or f"basis의 이름: {name}" not in existing_lines]
    if not pending:
        return path

    pending_basis = [F for _, F in pending]
    maxima_arrays = get_basis_maxima(pending_basis)
    minima_arrays = get_basis_minima(pending_basis)
    blocks = []
    for (name, _), maxima, minima in zip(pending, maxima_arrays, minima_arrays):
        metrics = _evaluate_maxima(maxima)
        blocks.append("\n".join([
            f"basis의 이름: {name}",
            "basis 별 최댓값: [" + ", ".join(f"{value:.4f}" for value in maxima) + "]",
            "basis 별 최솟값: [" + ", ".join(f"{value:.4f}" for value in minima) + "]",
            "",
            *(f"{metric}: {value}" for metric, value in metrics.items()),
        ]))
    if overwrite_existing:
        updated_text = existing_text
        for (name, _), block in zip(pending, blocks):
            pattern = re.compile(
                r"^" + re.escape(f"basis의 이름: {name}")
                + r"\n.*?(?=^basis의 이름: |\Z)", re.MULTILINE | re.DOTALL
            )
            if pattern.search(updated_text):
                updated_text = pattern.sub(lambda _: block + "\n\n", updated_text)
            else:
                updated_text = updated_text.rstrip("\n") + ("\n\n" if updated_text else "") + block + "\n"
        path.write_text(updated_text.rstrip("\n") + "\n", encoding="utf-8")
        return path

    trailing_newlines = len(existing_text) - len(existing_text.rstrip("\n"))
    separator = "\n" * max(0, 2 - trailing_newlines) if existing_text else ""
    with path.open("a", encoding="utf-8") as report:
        report.write(separator + "\n\n".join(blocks) + "\n")
    return path


def save_basis_report_without_endpoint(
    BasisArr, file_path="basis_report_without_endpoint.txt", basis_names=None,
    *, overwrite_existing=False,
):
    """끝점 조건을 만족하는 경계 기저를 제외하고 기존 상세 형식으로 저장한다.

    각 집합에서 첫 함수의 F(-1)=1, 마지막 함수의 F(1)=1 조건을 각각
    검사해 해당 함수 전체를 제외한다(절대 허용오차 1e-10).
    내부 함수는 최댓값이 1이어도 유지한다. 경계 함수가 내부에서 1보다
    커져도 끝점 조건을 만족하면 제외한다. 양 끝의 x값만 제외하는 것은 아니다.
    표시되는 최댓값 배열과 다섯 지표 모두 제외 후의 함수들에 해당한다.
    각주: 고정된 끝점 조건을 갖는 경계 기저를 제외한 봉우리 높이의 균일성이다.
    남은 함수가 없으면 ValueError를 발생시키며 파일을 변경하지 않는다.
    이름 추론, 중복 생략, 명시적 갱신, 네 자리 출력은 기존 함수를 재사용한다.
    """
    BasisArr = list(BasisArr)
    names = _resolve_basis_names(BasisArr, basis_names)
    return save_basis_report(
        _exclude_endpoint_basis(BasisArr, names), file_path, names,
        overwrite_existing=overwrite_existing,
    )


def _exclude_endpoint_basis(BasisArr, names):
    """상세 보고서와 비교표가 동일한 끝점 제외 규칙을 공유한다."""
    x = sp.symbols('x')
    filtered = []
    for name, F in zip(names, BasisArr):
        expressions = [sp.sympify(F_i) for F_i in F]
        if not expressions:
            raise ValueError(f"Basis {name} is empty.")
        if any(expression.free_symbols - {x} for expression in expressions):
            raise ValueError("Substitute all parameters other than x before evaluation.")
        excluded = set()
        for index, endpoint in ((0, -1), (len(expressions) - 1, 1)):
            value = float(expressions[index].subs(x, endpoint))
            if np.isclose(value, 1.0, rtol=0, atol=1e-10):
                excluded.add(index)
        remaining = [expression for index, expression in enumerate(expressions) if index not in excluded]
        if not remaining:
            raise ValueError(f"Basis {name} has no functions left after endpoint exclusion.")
        filtered.append(remaining)
    return filtered


def save_basis_comparison_without_endpoint(
    BasisArr, file_path="basis_comparison_without_endpoint.txt", basis_names=None,
):
    """끝점 기저 제외 후 열은 기저, 행은 다섯 지표인 비교표를 저장한다.

    제외 조건은 save_basis_report_without_endpoint와 같으며, 최댓값 캐시와
    기존 비교표 저장 함수를 재사용한다. 소수점 네 자리로 표시하고,
    같은 파일은 현재 결과로 덮어쓴다. 저장한 절대 경로를 반환한다.
    """
    BasisArr = list(BasisArr)
    names = _resolve_basis_names(BasisArr, basis_names)
    return save_basis_comparison(_exclude_endpoint_basis(BasisArr, names), file_path, names)


def _text_width(text):
    """고정폭 글꼴에서 한글 등 전각 문자를 두 칸으로 계산한다."""
    return sum(0 if unicodedata.combining(char) else
               2 if unicodedata.east_asian_width(char) in "WF" else 1 for char in text)


def save_basis_comparison(BasisArr, file_path="basis_comparison.txt", basis_names=None):
    """열은 기저 이름, 행은 다섯 지표인 비교표를 UTF-8 TXT로 저장한다.

    최댓값 배열은 표시하지 않으며 측정값은 소수점 네 자리로 표시한다.
    basis_names의 동작은 상세 보고서와 같다. 같은 파일은 현재 비교표로 덮어쓴다.
    빈 BasisArr는 허용하지 않는다. 저장한 파일의 절대 경로를 반환한다.
    """
    BasisArr = list(BasisArr)
    if not BasisArr:
        raise ValueError("Provide at least one basis set for comparison.")
    names = _resolve_basis_names(BasisArr, basis_names)
    results = [_evaluate_maxima(maxima) for maxima in get_basis_maxima(BasisArr)]
    rows = [["metric", *names]]
    rows.extend([metric, *(result[metric] for result in results)] for metric, _ in _METRICS)
    widths = [max(_text_width(cell) for cell in column) for column in zip(*rows)]
    lines = [" | ".join(cell + " " * (width - _text_width(cell))
                        for cell, width in zip(row, widths)) for row in rows]
    lines.insert(1, "-+-".join("-" * width for width in widths))
    path = Path(file_path).resolve()
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def save_basis_degree_comparison(
    node_function, n_initial, n_min, n_max,
    file_path="basis_degree_comparison.txt", *, without_endpoint=False,
):
    """한 노드 배치 방식의 차수별 지표 비교표를 UTF-8 TXT로 저장한다.

    node_function은 generate_nodes.py의 함수처럼 n개 노드를 반환해야 한다.
    고정된 기저 배열 대신 생성 함수를 받는 이유는 차수마다 기저를 새로 만들기 때문이다.
    n_initial부터 기저를 생성하되 표의 열은 n_min부터 n_max까지 오름차순이다.
    행은 range, ratio, mean, std, cv이며 소수점 네 자리로 표시한다.
    기본값은 끝점 기저 포함, without_endpoint=True이면 기존 끝점 제외 규칙을 쓴다.
    모든 차수를 계산하므로 범위가 크면 시간이 걸린다. 기존 최댓값 캐시를 공유한다.
    같은 파일은 덮어쓰며 계산 실패 시 파일을 변경하지 않는다. 절대 경로를 반환한다.
    """
    import generate_basis as basis_generator

    if not callable(node_function):
        raise TypeError("node_function must be callable.")
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
           for value in (n_initial, n_min, n_max)):
        raise ValueError("Degree values must be integers.")
    if not 1 <= n_min <= n_initial <= n_max:
        raise ValueError("Require 1 <= n_min <= n_initial <= n_max.")

    degrees = list(range(n_min, n_max + 1))
    bases = {}
    for degree in [n_initial, *(n for n in degrees if n != n_initial)]:
        try:
            nodes = np.asarray(node_function(degree))
            if np.iscomplexobj(nodes):
                if np.any(np.abs(nodes.imag) > 1e-10):
                    raise ValueError("Nodes must be real.")
                nodes = nodes.real
            nodes = np.asarray(nodes, dtype=float)
            if (nodes.ndim != 1 or len(nodes) != degree
                    or not np.all(np.isfinite(nodes)) or len(np.unique(nodes)) != degree):
                raise ValueError("Expected n distinct finite nodes.")
            bases[degree] = basis_generator.get_New_Basis(nodes)
        except Exception as error:
            raise ValueError(f"Cannot generate basis for n={degree}: {error}") from error

    save = save_basis_comparison_without_endpoint if without_endpoint else save_basis_comparison
    return save([bases[n] for n in degrees], file_path, [f"n={n}" for n in degrees])
