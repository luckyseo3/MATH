"""끝 노드를 고정하지 않은 0-symmetric 노드에서 분산이 0인 곡선을 구한다.

n=4 노드 (-b, -a, a, b), n=5 노드 (-b, -a, 0, a, b)에서 대칭이므로 서로 다른 내부 최댓값은 M_1, M_2 두 개이다.
따라서 분산 = 0  <=>  M_1 = M_2 이고, 식 하나에 변수 (a, b) 두 개이므로 해는 곡선을 이룬다.

방법: draw_variance.variance_plane(endpoint_variable=None)이 저장한 1/N 격자에서, 행 a = i/N마다
g(b) = M_1 - M_2의 부호가 바뀌는 이웃 격자점 쌍(양 끝 모두 양수 가중치)을 찾아 그 사이에서 b를 이분법으로 푼다.
g는 calculate_extremum.precise_extrema_and_weights(mpmath)로 계산하고,
마지막 b는 digits자리 십진수로 반올림한 Fraction에서 calculate_extremum / calculate_variance로 분산을 정확히 다시 구한다.
결과는 cache_dir/n{n}/N{N}/zero_curve.json에 저장하며, 설정이 같으면 다시 계산하지 않는다.

kind="min"이면 같은 방법으로 최솟값 m_1 = m_2인 곡선을 구한다 (cache_dir/n{n}/N{N}_min/zero_curve.json).
common_zero_points는 두 곡선의 교점, 즉 최댓값과 최솟값의 분산이 동시에 0인 노드를 구한다.
"""
import json
import multiprocessing
import os
from fractions import Fraction
from pathlib import Path

import numpy as np

import calculate_extremum
import calculate_variance
import draw_variance


SYMBOL = {"max": "M", "min": "m"}  # 최댓값 M_i, 최솟값 m_i


def zero_curve_path(cache_dir, n, N, kind="max"):
    return draw_variance.grid_dir(cache_dir, n, N, kind) / 'zero_curve.json'


def gap(n, a, b, kind, dps):
    """g = E_1 - E_2 (E = M 또는 m)를 mpmath dps자리로 계산한다. 노드는 plane_nodes(n, a, b, None)."""
    E, _ = calculate_extremum.precise_extrema_and_weights(draw_variance.plane_nodes(n, a, b, None), dps, kind=kind)
    return E[1] - E[2]


def bisect(g, low, high, digits):
    """g(low)와 g(high)의 부호가 다른 [low, high]에서 g = 0인 점을 이분법으로 구해 digits자리 Fraction으로 반올림한다."""
    low, high = Fraction(low), Fraction(high)
    g_low = g(low)
    width = Fraction(1, 10 ** (digits + 2))
    while high - low > width:
        middle = (low + high) / 2
        g_middle = g(middle)
        if g_middle == 0:
            low = high = middle
            break
        if (g_middle > 0) == (g_low > 0):
            low, g_low = middle, g_middle
        else:
            high = middle
    return Fraction(round((low + high) / 2 * 10 ** digits), 10 ** digits)


def load_extremum_plane(N, *, n, cache_dir, kind="max"):
    """(N+1)x(N+1)x(n+1) 격자 E[i, j](a=i/N, b=j/N에서 F_0~F_n의 최댓값 또는 최솟값)와 positive."""
    E = np.full((N + 1, N + 1, n + 1), np.nan)
    positive = np.zeros((N + 1, N + 1), dtype=bool)
    for i in draw_variance.plane_rows(N, n, cache_dir, kind):
        with np.load(draw_variance.volume_slab_path(cache_dir, n, N, i, kind)) as data:
            end = i + 1 + len(data['M'])
            E[i, i + 1:end] = data['M']
            positive[i, i + 1:end] = data['positive']
    return E, positive


def sign_change_brackets(N, n, cache_dir, kind="max"):
    """행 i마다 g = E_1 - E_2의 부호가 바뀌는 (j, j+1) 쌍. 두 점 모두 양수 가중치여야 한다."""
    brackets = {}
    for i in draw_variance.plane_rows(N, n, cache_dir, kind):
        with np.load(draw_variance.volume_slab_path(cache_dir, n, N, i, kind)) as data:
            M, positive = data['M'], data['positive']
        g = M[:, 1] - M[:, 2]
        ok = positive & np.isfinite(g)
        change = ok[:-1] & ok[1:] & (np.sign(g[:-1]) != np.sign(g[1:]))
        pairs = [(i + 1 + int(k), i + 2 + int(k)) for k in np.nonzero(change)[0]]
        if pairs:
            brackets[i] = pairs
    return brackets


def solve_row(task):
    """a = i/N에서 각 구간 [j/N, (j+1)/N] 안의 E_1 = E_2 근 b를 이분법으로 digits자리까지 구한다."""
    from mpmath import mp
    n, N, i, pairs, digits, kind = task
    dps = digits + 20
    a = Fraction(i, N)

    points = []
    with mp.workdps(dps):
        for j_low, j_high in pairs:
            b = bisect(lambda b: gap(n, a, b, kind, dps), Fraction(j_low, N), Fraction(j_high, N), digits)

            # digits자리 b에서 Fraction으로 정확히 다시 계산한다.
            node = draw_variance.plane_nodes(n, a, b, None)
            E, weights, labels = calculate_extremum.exact_extrema_and_weights(node, active=True, kind=kind)
            variance = calculate_variance.exact_variance(E[1:-1])
            points.append({'a': str(a), 'b': str(b), 'b_decimal': f'{float(b):.17g}',
                           f'common_{kind}': float(E[1]), 'variance': float(variance),
                           'min_weight': float(min(weights)), 'labels': ",".join(labels[1:-1])})
    return points


def volume_gaps(n, a, b, c, dps):
    """(g_1, g_2) = (M_1 - M_2, M_2 - M_3)를 mpmath dps자리로 계산한다. 노드는 volume_nodes(n, (a, b, c), None)."""
    M, _ = calculate_extremum.precise_extrema_and_weights(draw_variance.volume_nodes(n, (a, b, c), None), dps,
                                                          kind="max")
    return [M[1] - M[2], M[2] - M[3]]


def volume_cells(N, n, cache_dir, i, last):
    """a = i/N 슬랩의 (b, c) 격자 칸 중 네 꼭짓점이 모두 양수 가중치이고 g_1, g_2의 부호가 모두 바뀌는 칸 (j, k)."""
    with np.load(draw_variance.volume_slab_path(cache_dir, n, N, i)) as data:
        M, positive = data['M'], data['positive']
    j, k = np.triu_indices(last - i, 1)
    cell_ok = np.ones((N + 1, N + 1), dtype=bool)
    for low in (1, 2):
        g = np.full((N + 2, N + 2), np.nan)
        g[j + i + 1, k + i + 1] = np.where(positive, M[:, low] - M[:, low + 1], np.nan)
        corners = np.stack([g[:-1, :-1], g[1:, :-1], g[:-1, 1:], g[1:, 1:]])
        finite = np.isfinite(corners).all(axis=0)
        signs = np.sign(np.where(finite, corners, 0))
        cell_ok &= finite & (signs.min(axis=0) != signs.max(axis=0))
    return [(int(j), int(k)) for j, k in zip(*np.nonzero(cell_ok))]


def solve_volume_row(task):
    """a = i/N에서 각 칸 중심부터 (g_1, g_2) = 0을 (b, c)에 대해 findroot(Newton)로 풀고 Fraction으로 정확히 확인한다."""
    from mpmath import mp
    n, N, i, cells, digits = task
    dps = digits + 20
    a = Fraction(i, N)
    points, messages = [], []
    with mp.workdps(dps):
        a_mp = mp.mpf(i) / N
        for j, k in cells:
            start = [mp.mpf(2 * j + 1) / (2 * N), mp.mpf(2 * k + 1) / (2 * N)]
            try:
                b, c = mp.findroot(lambda b, c: volume_gaps(n, a_mp, b, c, dps), start,
                                   tol=mp.mpf(10) ** (-2 * digits))
            except (ValueError, ZeroDivisionError) as error:
                messages.append(f"a={i}/{N}, cell ({j}, {k}): findroot failed ({str(error).splitlines()[0]})")
                continue
            # round(mpf)는 float를 거치므로 mp.nint로 반올림한다.
            b, c = (Fraction(int(mp.nint(value * 10 ** digits)), 10 ** digits) for value in (b, c))
            if not (a < b < c <= 1 and abs(b - Fraction(2 * j + 1, 2 * N)) <= Fraction(3, 2 * N)
                    and abs(c - Fraction(2 * k + 1, 2 * N)) <= Fraction(3, 2 * N)):
                messages.append(f"a={i}/{N}, cell ({j}, {k}): left the cell or a < b < c <= 1 "
                                f"(b={float(b):.6f}, c={float(c):.6f})")
                continue
            if any(abs(Fraction(p['b']) - b) + abs(Fraction(p['c']) - c) < Fraction(1, 10 ** (digits // 2))
                   for p in points):
                continue  # 이웃 칸에서 같은 해로 수렴했다.

            node = draw_variance.volume_nodes(n, (a, b, c), None)
            M, weights, labels = calculate_extremum.exact_extrema_and_weights(node, active=True, kind="max")
            points.append({'a': str(a), 'b': str(b), 'c': str(c),
                           'b_decimal': f'{float(b):.17g}', 'c_decimal': f'{float(c):.17g}',
                           'common_max': float(M[1]), 'variance': float(calculate_variance.exact_variance(M[1:-1])),
                           'min_weight': float(min(weights)), 'labels': ",".join(labels[1:-1]), 'cell': [j, k]})
    return points, messages


def zero_variance_curve_volume(N, *, n=6, cache_dir, digits=30, workers=None):
    """n=6, 7 (끝 노드 c가 변수) 노드 (-c, -b, -a, (0), a, b, c)에서 최댓값 분산이 0인 곡선 위의 점을 구한다.

    대칭이므로 서로 다른 내부 최댓값은 M_1, M_2, M_3이고, 분산 = 0 <=> M_1 = M_2 = M_3 (식 두 개, 변수 세 개)이다.
    draw_variance.variance_volume(N)이 저장한 격자에서 a = i/N마다, 네 꼭짓점이 모두 양수 가중치이고
    g_1 = M_1 - M_2, g_2 = M_2 - M_3의 부호가 모두 바뀌는 (b, c) 칸을 고른다.
    칸 중심에서 (g_1, g_2) = 0을 mpmath findroot로 풀고, 해가 그 칸(이웃 칸까지)과 a < b < c <= 1 안에 있을 때만 남긴다.
    digits자리로 반올림한 (b, c)에서 분산과 가중치를 Fraction으로 정확히 다시 계산한다.
    같은 a에 해가 여러 개일 수 있다(곡선의 여러 가지). 결과는 cache_dir/n{n}/N{N}/zero_curve.json에 저장한다.
    """
    if n not in (6, 7):
        raise ValueError("M_1 = M_2 = M_3 curve is for n = 6, 7 with free endpoints")
    path = zero_curve_path(cache_dir, n, N)
    settings = {'n': n, 'N': N, 'digits': digits,
                'method': 'cells with g_1, g_2 sign changes per a, findroot on (b, c), exact check'}
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['settings'] == settings:
            return saved['points']

    last = draw_variance.volume_last(N, n, cache_dir)
    tasks = [(n, N, i, cells, digits) for i in range(1, last - 1)
             if (cells := volume_cells(N, n, cache_dir, i, last))]
    print(f"n={n}: solving M_1 = M_2 = M_3 from {sum(len(t[3]) for t in tasks)} cells in {len(tasks)} rows",
          flush=True)
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        results = pool.map(solve_volume_row, tasks)
    points = [point for row, _ in results for point in row]
    messages = [message for _, row in results for message in row]
    print(f"n={n}: {len(points)} points, {len(messages)} cells rejected", flush=True)

    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps({'settings': settings, 'points': points, 'rejected': messages}, indent=1),
                    encoding='utf-8')
    os.replace(temp, path)
    return points


def draw_zero_curve_volume(points, *, n=6):
    """M_1 = M_2 = M_3 곡선을 (a, b), (a, c) 투영과 a에 따른 공통 최댓값으로 그린다. 점 색은 최댓값 후보 조합이다."""
    import matplotlib.pyplot as plt

    a = np.array([float(Fraction(p['a'])) for p in points])
    b = np.array([float(p['b_decimal']) for p in points])
    c = np.array([float(p['c_decimal']) for p in points])
    common = np.array([p['common_max'] for p in points])
    labels = [p['labels'] for p in points]
    names = list(dict.fromkeys(labels))
    colors = {name: plt.get_cmap("tab10")(k % 10) for k, name in enumerate(names)}
    best = int(np.argmin(common))

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    for name in names:
        mask = np.array([label == name for label in labels])
        text = draw_variance.label_text(name)
        for ax, y in zip(axes, (b, c, common)):
            ax.scatter(a[mask], y[mask], s=6, color=colors[name], label=f"{text}  {int(mask.sum())}")
    for ax, y in zip(axes, (b, c, common)):
        ax.scatter([a[best]], [y[best]], marker="o", s=25, color="magenta", edgecolors="black", zorder=5)
        ax.set_xlabel("a")
        ax.grid(True, alpha=0.3)
    axes[0].set(ylabel="b", title=f"n={n}: M_1 = M_2 = M_3 curve, {len(points)} points (a, b)")
    axes[1].set(ylabel="c", title="(a, c)")
    axes[2].set(ylabel="common maximum M_1 = M_2 = M_3",
                title=f"min common maximum {common[best]:.10f}\n"
                      f"at a={a[best]:.3f}, b={b[best]:.10f}, c={c[best]:.10f}")
    axes[0].legend(fontsize=7, title="argmax of (F_1, ..., F_{n-1})", title_fontsize=7)
    fig.tight_layout()
    return fig, axes


def zero_variance_curve(N, *, n=4, cache_dir, digits=30, workers=None, kind="max"):
    """분산이 0인 곡선 위의 점 목록을 구한다 (a = i/N마다, a 오름차순). 저장된 결과가 있으면 읽는다.

    kind="max"이면 최댓값 M_1 = M_2, kind="min"이면 최솟값 m_1 = m_2인 곡선이다
    (draw_variance.variance_plane(N, kind=kind)의 격자가 먼저 있어야 한다).
    각 점은 a, b(Fraction 문자열), common_max 또는 common_min(= E_1 = E_2), variance(정확한 분산의 float),
    min_weight, labels(F_1~F_{n-1}의 최댓값 또는 최솟값 후보)를 가진다.
    """
    if n not in (4, 5):
        raise ValueError("E_1 = E_2 curve is for n = 4, 5 with free endpoints")
    path = zero_curve_path(cache_dir, n, N, kind)
    settings = {'n': n, 'N': N, 'digits': digits, 'method': 'bisection on M_1 - M_2 (mpmath), exact check'}
    if kind != "max":  # 기존 최댓값 결과의 설정은 그대로 두어 다시 계산하지 않는다.
        settings |= {'kind': kind, 'method': 'bisection on m_1 - m_2 (mpmath), exact check'}
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['settings'] == settings:
            return saved['points']

    brackets = sign_change_brackets(N, n, cache_dir, kind)
    tasks = [(n, N, i, pairs, digits, kind) for i, pairs in sorted(brackets.items())]
    print(f"n={n}: solving {SYMBOL[kind]}_1 = {SYMBOL[kind]}_2 on {sum(len(t[3]) for t in tasks)} brackets "
          f"in {len(tasks)} rows", flush=True)
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        points = [point for row in pool.map(solve_row, tasks) for point in row]

    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps({'settings': settings, 'points': points}, indent=1), encoding='utf-8')
    os.replace(temp, path)
    return points


def draw_zero_curve(points, *, n=4, kind="max"):
    """분산이 0인 곡선 b(a)와, 곡선을 따라간 공통값 E_1 = E_2를 그린다. 점 색은 최댓값(최솟값) 후보 조합이다.

    kind="max"이면 공통 최댓값이 가장 작은 점, kind="min"이면 공통 최솟값이 가장 큰 점(|m|이 가장 작은 점)을 표시한다.
    """
    import matplotlib.pyplot as plt

    E = SYMBOL[kind]
    a = np.array([float(Fraction(p['a'])) for p in points])
    b = np.array([float(Fraction(p['b'])) for p in points])
    common = np.array([p[f'common_{kind}'] for p in points])
    labels = [p['labels'] for p in points]
    names = list(dict.fromkeys(labels))
    colors = {name: plt.get_cmap("tab10")(k % 10) for k, name in enumerate(names)}

    fig, (ax_curve, ax_max) = plt.subplots(1, 2, figsize=(13, 5.5))
    for name in names:
        mask = np.array([label == name for label in labels])
        text = draw_variance.label_text(name)
        ax_curve.scatter(a[mask], b[mask], s=6, color=colors[name], label=f"{text}  {int(mask.sum())}")
        ax_max.scatter(a[mask], np.abs(common[mask]), s=6, color=colors[name], label=text)
    best = int(np.argmin(np.abs(common)))
    for ax, y in ((ax_curve, b[best]), (ax_max, abs(common[best]))):
        ax.scatter([a[best]], [y], marker="o", s=20, color="magenta", edgecolors="black", zorder=5)
    ax_curve.set(xlabel="a", ylabel="b", xlim=(0, 1), ylim=(0, 1),
                 title=f"n={n}: {kind}imum variance = 0 curve ({E}_1 = {E}_2), {len(points)} points")
    ax_curve.plot([0, 1], [0, 1], color="0.7", linewidth=0.8)
    ax_curve.legend(fontsize=8, title=f"arg{kind} of (F_1, ..., F_{{n-1}})", title_fontsize=8)
    # 공통값의 크기가 크게 변하므로 로그 축으로 그린다. 내부 F_i는 ±1에서 0이므로 최솟값은 0 이하이고 -m_1을 그린다.
    ax_max.set_yscale("log")
    if kind == "max":
        ax_max.set(ylabel="common maximum M_1 = M_2 (log scale)",
                   title=f"min common maximum {common[best]:.10f}\nat a={a[best]:.3f}, b={b[best]:.10f}")
    else:
        ax_max.set(ylabel="-(common minimum) = -m_1 = -m_2 (log scale)",
                   title=f"max common minimum {common[best]:.10f}\nat a={a[best]:.3f}, b={b[best]:.10f}")
    ax_max.set(xlabel="a")
    ax_max.grid(True)
    fig.tight_layout()
    return fig, (ax_curve, ax_max)


def common_zero_path(cache_dir, n):
    return Path(cache_dir) / f'n{n}' / 'common_zero.json'


def common_zero_points(n, *, cache_dir, grids=None, digits=30):
    """최댓값의 분산과 최솟값의 분산이 동시에 0인 점 (M_1 = M_2이고 m_1 = m_2)을 모두 구한다.

    grids = {kind: N}은 variance_plane(N, kind=kind)로 저장한 격자이다 (기본값 max 1/1000, min 1/100).
    두 격자를 더 성긴 간격 1/N_c로 맞추고, 네 꼭짓점이 모두 양수 가중치이며 g_M = M_1 - M_2와
    g_m = m_1 - m_2의 부호가 모두 바뀌는 격자 칸마다 칸 중심에서 (g_M, g_m) = 0을 mpmath findroot(Newton)로 푼다.
    해는 digits자리로 반올림한 Fraction에서 두 분산과 가중치를 정확히 다시 계산한다.
    결과는 cache_dir/n{n}/common_zero.json에 저장하며, 설정이 같으면 다시 계산하지 않는다.
    반환값은 {'cells': 후보 칸 수, 'points': [해 ...]}이다.
    """
    from mpmath import mp

    grids = grids or {"max": 1000, "min": 100}
    path = common_zero_path(cache_dir, n)
    settings = {'n': n, 'grids': grids, 'digits': digits, 'method': 'cells with both sign changes, findroot, exact check'}
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['settings'] == settings:
            return saved

    N = min(grids.values())
    g, ok = {}, np.ones((N + 1, N + 1), dtype=bool)
    for kind, N_kind in grids.items():
        if N_kind % N:
            raise ValueError(f"grid 1/{N_kind} is not a refinement of 1/{N}")
        E, positive = load_extremum_plane(N_kind, n=n, cache_dir=cache_dir, kind=kind)
        step = N_kind // N
        g[kind] = (E[:, :, 1] - E[:, :, 2])[::step, ::step]
        ok &= positive[::step, ::step] & np.isfinite(g[kind])

    def corners(array):
        return np.stack([array[:-1, :-1], array[1:, :-1], array[:-1, 1:], array[1:, 1:]])

    cell_ok = corners(ok).all(axis=0)
    for kind in grids:
        signs = np.sign(corners(g[kind]))
        cell_ok &= signs.min(axis=0) != signs.max(axis=0)
    cells = [(int(i), int(j)) for i, j in zip(*np.nonzero(cell_ok))]

    dps = digits + 20
    points = []
    with mp.workdps(dps):
        for i, j in cells:
            start = [mp.mpf(2 * i + 1) / (2 * N), mp.mpf(2 * j + 1) / (2 * N)]
            try:
                a, b = mp.findroot(lambda a, b: (gap(n, a, b, "max", dps), gap(n, a, b, "min", dps)),
                                   start, tol=mp.mpf(10) ** (-2 * digits))
            except (ValueError, ZeroDivisionError) as error:  # 수렴하지 않으면 칸만 기록한다.
                print(f"cell a={i}/{N}, b={j}/{N}: findroot failed ({str(error).splitlines()[0]})", flush=True)
                continue
            # round(mpf)는 float를 거치므로 mp.nint로 반올림한다.
            a, b = (Fraction(int(mp.nint(value * 10 ** digits)), 10 ** digits) for value in (a, b))
            if not (0 < a < b <= 1 and abs(a - Fraction(2 * i + 1, 2 * N)) <= Fraction(3, 2 * N)
                    and abs(b - Fraction(2 * j + 1, 2 * N)) <= Fraction(3, 2 * N)):
                print(f"cell a={i}/{N}, b={j}/{N}: findroot left the cell (a={float(a):.6f}, b={float(b):.6f})",
                      flush=True)
                continue  # 해가 그 칸(이웃 칸까지) 밖이거나 0 < a < b <= 1 밖이다.
            if any(abs(Fraction(p['a']) - a) + abs(Fraction(p['b']) - b) < Fraction(1, 10 ** (digits // 2))
                   for p in points):
                continue  # 이웃 칸에서 같은 해로 수렴했다.

            # digits자리 (a, b)에서 두 분산과 가중치를 Fraction으로 정확히 다시 계산한다.
            node = draw_variance.plane_nodes(n, a, b, None)
            point = {'a': str(a), 'b': str(b), 'a_decimal': f'{float(a):.17g}', 'b_decimal': f'{float(b):.17g}',
                     'cell': [int(i), int(j)]}
            for kind in grids:
                E, weights, labels = calculate_extremum.exact_extrema_and_weights(node, active=True, kind=kind)
                point |= {f'common_{kind}': float(E[1]),
                          f'variance_{kind}': float(calculate_variance.exact_variance(E[1:-1])),
                          f'labels_{kind}': ",".join(labels[1:-1])}
            point['min_weight'] = float(min(weights))
            points.append(point)

    result = {'settings': settings, 'cells': len(cells), 'points': points}
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(result, indent=1), encoding='utf-8')
    os.replace(temp, path)
    return result


def draw_zero_curves(curves, *, n=4, common=None, cache_dir=None, N=None, zoom=None):
    """최댓값 분산 0 곡선(M_1 = M_2)과 최솟값 분산 0 곡선(m_1 = m_2)을 한 (a, b) 평면에 겹쳐 그린다.

    curves = {kind: zero_variance_curve(..., kind=kind)의 점 목록}, common = common_zero_points의 points.
    cache_dir, N이 있으면 그 격자(최솟값 격자)에서 음수 가중치 영역을 회색으로 칠한다.
    zoom = (a_low, a_high, b_low, b_high)이면 그 범위를 확대한 작은 그림을 곡선이 없는 왼쪽 아래에 넣는다.
    """
    import matplotlib.pyplot as plt

    style = {"max": dict(color="tab:blue", label="max: M_1 = M_2"),
             "min": dict(color="tab:red", label="min: m_1 = m_2")}
    fig, ax = plt.subplots(figsize=(8, 7))
    if cache_dir is not None and N is not None:
        _, positive = load_extremum_plane(N, n=n, cache_dir=cache_dir, kind="min")
        defined = np.triu(np.ones((N + 1, N + 1), dtype=bool), 1)
        defined[0] = False
        negative = np.where(defined & ~positive, 1.0, np.nan)
        ax.pcolormesh(np.arange(N + 1) / N, np.arange(N + 1) / N, negative.T, cmap="Greys", vmin=0, vmax=3,
                      shading="nearest", zorder=0)
        ax.plot([], [], color="0.8", linewidth=6, label="negative weight")
    axes = [ax]
    if zoom is not None:
        inset = ax.inset_axes([0.07, 0.07, 0.42, 0.42], xlim=zoom[:2], ylim=zoom[2:])
        ax.indicate_inset_zoom(inset, edgecolor="black")
        axes.append(inset)
    for target in axes:
        for kind, points in curves.items():
            a = np.array([float(Fraction(p['a'])) for p in points])
            b = np.array([float(Fraction(p['b'])) for p in points])
            target.scatter(a, b, s=5, **style[kind] | {'label': f"{style[kind]['label']}  ({len(points)} points)"})
        for p in common or []:
            target.scatter([float(p['a_decimal'])], [float(p['b_decimal'])], s=90, marker="*", color="gold",
                           edgecolors="black", zorder=5,
                           label=f"both: a={float(p['a_decimal']):.10f}, b={float(p['b_decimal']):.10f}")
    ax.plot([0, 1], [0, 1], color="0.7", linewidth=0.8)
    ax.set(xlabel="a", ylabel="b", xlim=(0, 1), ylim=(0, 1),
           title=f"n={n}: variance = 0 curves of the maxima and the minima")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    return fig, ax
