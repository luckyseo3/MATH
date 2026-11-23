"""끝 노드를 고정하지 않은 0-symmetric 노드에서 분산이 0인 곡선을 구한다.

n=4 노드 (-b, -a, a, b), n=5 노드 (-b, -a, 0, a, b)에서 대칭이므로 서로 다른 내부 최댓값은 M_1, M_2 두 개이다.
따라서 분산 = 0  <=>  M_1 = M_2 이고, 식 하나에 변수 (a, b) 두 개이므로 해는 곡선을 이룬다.

방법: draw_variance.variance_plane(endpoint_variable=None)이 저장한 1/N 격자에서, 행 a = i/N마다
g(b) = M_1 - M_2의 부호가 바뀌는 이웃 격자점 쌍(양 끝 모두 양수 가중치)을 찾아 그 사이에서 b를 이분법으로 푼다.
g는 calculate_maximum.precise_maximums_and_weights(mpmath)로 계산하고,
마지막 b는 digits자리 십진수로 반올림한 Fraction에서 calculate_maximum / calculate_variance로 분산을 정확히 다시 구한다.
결과는 cache_dir/n{n}/N{N}/zero_curve.json에 저장하며, 설정이 같으면 다시 계산하지 않는다.
"""
import json
import multiprocessing
import os
from fractions import Fraction
from pathlib import Path

import numpy as np

import calculate_maximum
import calculate_variance
import draw_variance


def zero_curve_path(cache_dir, n, N):
    return Path(cache_dir) / f'n{n}' / f'N{N}' / 'zero_curve.json'


def sign_change_brackets(N, n, cache_dir):
    """행 i마다 g = M_1 - M_2의 부호가 바뀌는 (j, j+1) 쌍. 두 점 모두 양수 가중치여야 한다."""
    brackets = {}
    for i in draw_variance.plane_rows(N, n, cache_dir):
        with np.load(draw_variance.volume_slab_path(cache_dir, n, N, i)) as data:
            M, positive = data['M'], data['positive']
        g = M[:, 1] - M[:, 2]
        ok = positive & np.isfinite(g)
        change = ok[:-1] & ok[1:] & (np.sign(g[:-1]) != np.sign(g[1:]))
        pairs = [(i + 1 + int(k), i + 2 + int(k)) for k in np.nonzero(change)[0]]
        if pairs:
            brackets[i] = pairs
    return brackets


def solve_row(task):
    """a = i/N에서 각 구간 [j/N, (j+1)/N] 안의 M_1 = M_2 근 b를 이분법으로 digits자리까지 구한다."""
    from mpmath import mp
    n, N, i, pairs, digits = task
    dps = digits + 20
    a = Fraction(i, N)

    def g(b):
        M, _ = calculate_maximum.precise_maximums_and_weights(draw_variance.plane_nodes(n, a, b, None), dps)
        return M[1] - M[2]

    points = []
    with mp.workdps(dps):
        for j_low, j_high in pairs:
            low, high = Fraction(j_low, N), Fraction(j_high, N)
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
            b = Fraction(round((low + high) / 2 * 10 ** digits), 10 ** digits)

            # digits자리 b에서 Fraction으로 정확히 다시 계산한다.
            node = draw_variance.plane_nodes(n, a, b, None)
            M, weights, labels = calculate_maximum.exact_maximums_and_weights(node, active=True)
            variance = calculate_variance.exact_variance(M[1:-1])
            points.append({'a': str(a), 'b': str(b), 'b_decimal': f'{float(b):.17g}',
                           'common_max': float(M[1]), 'variance': float(variance),
                           'min_weight': float(min(weights)), 'labels': ",".join(labels[1:-1])})
    return points


def zero_variance_curve(N, *, n=4, cache_dir, digits=30, workers=None):
    """분산이 0인 곡선 위의 점 목록을 구한다 (a = i/N마다, a 오름차순). 저장된 결과가 있으면 읽는다.

    각 점은 a, b(Fraction 문자열), common_max(= M_1 = M_2), variance(정확한 분산의 float),
    min_weight, labels(F_1~F_{n-1}의 최댓값 후보)를 가진다.
    """
    if n not in (4, 5):
        raise ValueError("M_1 = M_2 curve is for n = 4, 5 with free endpoints")
    path = zero_curve_path(cache_dir, n, N)
    settings = {'n': n, 'N': N, 'digits': digits, 'method': 'bisection on M_1 - M_2 (mpmath), exact check'}
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['settings'] == settings:
            return saved['points']

    brackets = sign_change_brackets(N, n, cache_dir)
    tasks = [(n, N, i, pairs, digits) for i, pairs in sorted(brackets.items())]
    print(f"n={n}: solving M_1 = M_2 on {sum(len(t[3]) for t in tasks)} brackets in {len(tasks)} rows", flush=True)
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        points = [point for row in pool.map(solve_row, tasks) for point in row]

    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps({'settings': settings, 'points': points}, indent=1), encoding='utf-8')
    os.replace(temp, path)
    return points


def draw_zero_curve(points, *, n=4):
    """분산이 0인 곡선 b(a)와, 곡선을 따라간 공통 최댓값 M_1 = M_2를 그린다. 점 색은 최댓값 후보 조합이다."""
    import matplotlib.pyplot as plt

    a = np.array([float(Fraction(p['a'])) for p in points])
    b = np.array([float(Fraction(p['b'])) for p in points])
    common = np.array([p['common_max'] for p in points])
    labels = [p['labels'] for p in points]
    names = list(dict.fromkeys(labels))
    colors = {name: plt.get_cmap("tab10")(k % 10) for k, name in enumerate(names)}

    fig, (ax_curve, ax_max) = plt.subplots(1, 2, figsize=(13, 5.5))
    for name in names:
        mask = np.array([label == name for label in labels])
        text = draw_variance.label_text(name)
        ax_curve.scatter(a[mask], b[mask], s=6, color=colors[name], label=f"{text}  {int(mask.sum())}")
        ax_max.scatter(a[mask], common[mask], s=6, color=colors[name], label=text)
    best = int(np.argmin(common))
    for ax, y in ((ax_curve, b[best]), (ax_max, common[best])):
        ax.scatter([a[best]], [y], marker="o", s=20, color="magenta", edgecolors="black", zorder=5)
    ax_curve.set(xlabel="a", ylabel="b", xlim=(0, 1), ylim=(0, 1),
                 title=f"n={n}: variance = 0 curve (M_1 = M_2), {len(points)} points")
    ax_curve.plot([0, 1], [0, 1], color="0.7", linewidth=0.8)
    ax_curve.legend(fontsize=8, title="argmax of (F_1, ..., F_{n-1})", title_fontsize=8)
    ax_max.set_yscale("log")  # a → 0에서 공통 최댓값이 크게 커지므로 로그 축으로 그린다.
    ax_max.set(xlabel="a", ylabel="common maximum M_1 = M_2 (log scale)",
               title=f"min common maximum {common[best]:.10f}\nat a={a[best]:.3f}, b={b[best]:.10f}")
    ax_max.grid(True)
    fig.tight_layout()
    return fig, (ax_curve, ax_max)
