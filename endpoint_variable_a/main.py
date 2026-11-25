"""끝 노드를 고정하지 않은 노드의 기저를 a, b 슬라이더로 조절해 본다.

1️⃣ 4개 노드 (-a, -b, b, a)
2️⃣ 5개 노드 (-a, -b, 0, b, a)
3️⃣ n=4 분산 0 곡선에서 공통 최댓값 M이 가장 작은 a*를 중심으로 a* ± 0.25를 0.001 간격으로 고른 기저를
   a 슬라이더로 넘겨 본다 (곡선이 끝나는 쪽은 곡선 끝에서 자른다).
4️⃣ n=5에서 같은 방식의 분산 0 기저 슬라이더
분산 0 슬라이더는 키보드 ←/→로 한 칸, ↑/↓로 열 칸, Home/End로 처음/끝 점으로 옮길 수 있다.
네 그림은 서로 다른 창에 그린다.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
VARIANCE_CACHE = pathlib.Path(__file__).resolve().parent / 'variance_cache'  # 이 폴더의 결과

from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from matplotlib.widgets import Slider

import construct_basis_fraction
import draw_basis as draw
import zero_variance_curve

GRID_N = 1000        # 분산 0 곡선을 구한 격자 (variance_cache/n{n}/N1000/zero_curve.json)
WINDOW_RADIUS = Fraction(1, 4)    # 최솟값 점 a*를 중심으로 볼 반지름 (a* ± 0.25)
WINDOW_STEP = Fraction(1, 1000)   # 그 안에서 고르는 a 간격 (0.001 = 곡선의 모든 점)
ARROW_STEPS = {'left': -1, 'right': 1, 'down': -10, 'up': 10}  # 분산 0 슬라이더의 키보드 이동 칸 수


def curve_nodes(n, a, b):
    """분산 0 곡선의 노드: n=4는 (-b, -a, a, b), n=5는 (-b, -a, 0, a, b) (a 안쪽, b 바깥쪽)."""
    return [-b, -a, *([0] if n % 2 else []), a, b]


def zero_variance_bases(n, radius=WINDOW_RADIUS, step=WINDOW_STEP):
    """n=4, 5 분산 0 곡선에서 공통 최댓값 M이 가장 작은 점 a*를 중심으로 [a* - radius, a* + radius]의
    점을 step 간격으로 골라 각 점의 기저를 만든다. 곡선 밖으로 나가는 쪽은 곡선의 끝에서 자른다.

    곡선은 n4.py, n5.py와 같은 표기로 노드 curve_nodes(n, a, b), 0 < a < b <= 1이며 a = i/GRID_N마다 점이 하나다.
    반환값은 (samples, center)이다. samples는 점마다 {'a', 'b', 'nodes', 'basis', 'common_max'}이고
    basis는 F_0~F_n의 Fraction 계수 목록(낮은 차수부터, 오차 없음), center는 a*의 samples 안 위치이다.
    """
    points = zero_variance_curve.zero_variance_curve(GRID_N, n=n, cache_dir=VARIANCE_CACHE)
    a_grid = [Fraction(point['a']) for point in points]
    stride = Fraction(step) * GRID_N
    if stride.denominator != 1 or stride < 1:
        raise ValueError(f'step must be a positive multiple of 1/{GRID_N}.')
    best = min(range(len(points)), key=lambda k: points[k]['common_max'])
    a_star, radius = a_grid[best], Fraction(radius)
    picks = [k for k in range(len(points))
             if abs(a_grid[k] - a_star) <= radius and (k - best) % int(stride) == 0]
    samples = []
    for k in picks:
        a_value, b_value = a_grid[k], Fraction(points[k]['b'])
        node = curve_nodes(n, a_value, b_value)
        samples.append({'a': float(a_value), 'b': float(b_value), 'nodes': node,
                        'basis': construct_basis_fraction.get_new_basis(node),
                        'common_max': points[k]['common_max']})
    return samples, picks.index(best)


def draw_zero_variance_basis(samples, n, center=0):
    """분산 0 곡선 위 기저를 a 슬라이더(또는 화살표 키)로 넘겨 본다. 처음에는 samples[center]를 보여 준다.

    점선은 공통 최댓값 M이며 내부 F_1~F_{n-1}의 최댓값이 모두 여기에 닿는다.
    """
    t = np.linspace(-1, 1, 1000)
    curves = [[np.polynomial.polynomial.polyval(t, [float(c) for c in F_i]) for F_i in sample['basis']]
              for sample in samples]
    node_text = '(-b, -a, 0, a, b)' if n % 2 else '(-b, -a, a, b)'

    fig, ax = plt.subplots(figsize=(10, 6.5))
    fig.subplots_adjust(bottom=0.18)
    lines = [ax.plot(t, curve, label=rf'$F_{{{i}}}$')[0] for i, curve in enumerate(curves[0])]
    node_markers = ax.scatter([float(v) for v in samples[0]['nodes']], np.zeros(n),
                              color='black', zorder=5, label='Nodes')
    max_line = ax.axhline(samples[0]['common_max'], color='gray', linestyle='--', linewidth=1,
                          label=rf'$M = M_1 = \cdots = M_{{{n - 1}}}$')
    ax.axhline(0, linewidth=0.8)
    ax.set(xlabel='t', ylabel='Basis value')
    ax.grid(True)
    ax.legend(loc='upper right')

    a_values = np.array([s['a'] for s in samples])
    slider_ax = fig.add_axes([0.2, 0.06, 0.65, 0.03])
    # valinit의 빨간 표시가 M이 가장 작은 점 a*를 가리킨다.
    slider = Slider(slider_ax, 'a', a_values[0], a_values[-1], valinit=a_values[center], valstep=a_values)
    current = center

    def show(index):
        nonlocal current
        current = index
        sample = samples[index]
        for line, curve in zip(lines, curves[index]):
            line.set_ydata(curve)
        node_markers.set_offsets(np.column_stack(([float(v) for v in sample['nodes']], np.zeros(n))))
        M = sample['common_max']
        max_line.set_ydata([M, M])
        ax.set_title(rf'n={n}, Var = 0 : nodes ${node_text}$, '
                     rf'$a={sample["a"]:.3f}$, $b={sample["b"]:.6f}$, $M={M:.6f}$'
                     f'   [{index + 1}/{len(samples)}]' + ('  min M' if index == center else ''))
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in ARROW_STEPS:
            index = current + ARROW_STEPS[event.key]
        elif event.key == 'home':
            index = 0
        elif event.key == 'end':
            index = len(samples) - 1
        else:
            return
        slider.set_val(a_values[min(max(index, 0), len(samples) - 1)])  # 슬라이더가 show를 부른다.

    slider.on_changed(lambda value: show(int(np.argmin(np.abs(a_values - value)))))
    fig.canvas.mpl_connect('key_press_event', on_key)
    show(center)
    fig._zero_variance_controls = slider
    return fig, ax, slider


if __name__ == "__main__":
    # ←/→는 matplotlib 기본 단축키(이전/다음 보기)이므로 빼서 분산 0 슬라이더에 쓴다.
    for keymap in ('keymap.back', 'keymap.forward'):
        plt.rcParams[keymap] = [key for key in plt.rcParams[keymap] if key not in ARROW_STEPS]
    a, b = sp.symbols("a b")
    a_initial, b_initial = 0.90, 0.40
    # {심볼: (최솟값, 최댓값, 초기값, 간격)}
    parameters = {
        a: (0.02, 1.0, a_initial, 0.01),
        b: (0.01, 0.99, b_initial, 0.01),
    }

    ##################################################################################################################
    # 1️⃣ 4개 노드 (-a, -b, b, a)의 기저
    nodes_4 = [-a, -b, b, a]
    fig_4, ax_4, sliders_4 = draw.draw_adjustable_basis(nodes_4, parameters)
    fig_4.canvas.manager.set_window_title('n=4 : (-a, -b, b, a)')
    ax_4.set_title(r'n=4, Basis functions for nodes $(-a, -b, b, a)$')
    ax_4.legend(loc="upper right")

    ##################################################################################################################
    # 2️⃣ 5개 노드 (-a, -b, 0, b, a)의 기저
    nodes_5 = [-a, -b, 0, b, a]
    fig_5, ax_5, sliders_5 = draw.draw_adjustable_basis(nodes_5, parameters)
    fig_5.canvas.manager.set_window_title('n=5 : (-a, -b, 0, b, a)')
    ax_5.set_title(r'n=5, Basis functions for nodes $(-a, -b, 0, b, a)$')
    ax_5.legend(loc="upper right")

    ##################################################################################################################
    # 3️⃣ n=4 분산 0 곡선에서 M이 가장 작은 a*의 ± 0.25, 0.001 간격 기저 (노드 (-b, -a, a, b), a 안쪽 · b 바깥쪽)
    # 곡선은 variance_cache/n4/N1000/zero_curve.json (n4.py가 만든 결과)에서 읽는다.
    zero_samples_4, zero_center_4 = zero_variance_bases(4)
    fig_zero_4, ax_zero_4, zero_slider_4 = draw_zero_variance_basis(zero_samples_4, 4, zero_center_4)
    fig_zero_4.canvas.manager.set_window_title('n=4 : Var = 0 curve bases')

    ##################################################################################################################
    # 4️⃣ n=5 분산 0 곡선에서 M이 가장 작은 a*의 ± 0.25, 0.001 간격 기저 (노드 (-b, -a, 0, a, b), a 안쪽 · b 바깥쪽)
    # 곡선은 variance_cache/n5/N1000/zero_curve.json (n5.py가 만든 결과)에서 읽는다.
    zero_samples_5, zero_center_5 = zero_variance_bases(5)
    fig_zero_5, ax_zero_5, zero_slider_5 = draw_zero_variance_basis(zero_samples_5, 5, zero_center_5)
    fig_zero_5.canvas.manager.set_window_title('n=5 : Var = 0 curve bases')

    plt.show()
