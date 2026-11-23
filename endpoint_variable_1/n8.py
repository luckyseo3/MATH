import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
from fractions import Fraction
ENDPOINT_VARIABLE = Fraction(1)  # 끝 노드 a: 노드의 양 끝은 ±a
VARIANCE_CACHE = pathlib.Path(__file__).resolve().parent / 'variance_cache'  # 이 a로 구한 결과
import numpy as np
from fractions import Fraction
import matplotlib.pyplot as plt
import sympy as sp

import nodes
import construct_basis_sympy as basis
import draw_basis as draw
import calculate_maximum
import calculate_variance
import draw_variance as variance_plot

# 격자 간격은 1/GRID_N이다.
GRID_N = 500

if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL":nodes.precise_GL_nodes,
        "ESC": nodes.ESC_nodes(8, dps=40),
        "case1": np.array([Fraction("-1.00000000000000000000"), Fraction("-0.84645906502257362014"), Fraction("-0.56485449735052128961"), Fraction("-0.20261236737965795940"), Fraction("0.20261236737965809818"), Fraction("0.56485449735052162268"), Fraction("0.84645906502257350912"), Fraction("1.00000000000000000000")], dtype=object),
        "case2": np.array([Fraction("-1.00000000000000000000"), Fraction("-0.83869806502259058867"), Fraction("-0.54418149735056398075"), Fraction("-0.19062836737964597433"), Fraction("0.19062836737964611311"), Fraction("0.54418149735056431382"), Fraction("0.83869806502259047765"), Fraction("1.00000000000000000000")], dtype=object),
    }

    ##################################################################################################################
    # 1️⃣ basis를 차수와 종류별로 나누어 관찰
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=8,
        n_min=3,
        n_max=15,
        initial_basis="GL",
    )

    ##################################################################################################################
    # 2️⃣
    a, b, c = sp.symbols("a b c")
    adjustable_nodes = [-1, -c, -b, -a, a, b, c, 1]
    adjustable_basis = basis.get_New_Basis(
        adjustable_nodes, basis_name='adjustable'
    )
    a_initial = 0.40
    b_initial = 0.50
    c_initial = 0.60
    a_scope = [0.1, 0.9, a_initial, 0.01]
    b_scope = [0.1, 0.9, b_initial, 0.01]
    c_scope = [0.1, 0.9, c_initial, 0.01]

    fig_a, ax_a, a_sliders = draw.draw_adjustable_basis(
        adjustable_nodes,
        {a: tuple(a_scope), b: tuple(b_scope), c: tuple(c_scope)},
    )
    ax_a.set_title(rf'n=8, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")

    ##################################################################################################################
    # 3️⃣ 노드 (-1, -c, -b, -a, a, b, c, 1)에서 a, b, c에 따른 내부 F_1~F_7 최댓값 분산
    # 유리수 격자점에서 Fraction으로 오차 없이 계산하며, 분산이 작을수록 점이 크다.
    # 결과는 variance_cache/n8/N{GRID_N}/에 저장되어 중단해도 이어서 계산한다.
    variance_plot.variance_volume(GRID_N, n=8, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)

    ##################################################################################################################
    # 4️⃣ 지정한 (a, b, c)에서 시작해 분산을 줄이는 compass search
    # 결과는 variance_cache/n8/compass.json에 저장되어 같은 탐색은 다시 계산하지 않는다.
    search_starts = {
        "outer 1": (Fraction("0.1"), Fraction("0.2"), Fraction("0.9")),
        "outer 2": (Fraction("0.1"), Fraction("0.2"), Fraction("0.5")),
        "inner 1": (Fraction("0.1"), Fraction("0.2"), Fraction("0.7")),
        "outer 3": (Fraction("0.3"), Fraction("0.4"), Fraction("0.5")),
        "inner 2": (Fraction("0.3"), Fraction("0.5"), Fraction("0.7")),
        "outer 4": (Fraction("0.3"), Fraction("0.8"), Fraction("0.9")),
        "outer 5": (Fraction("0.5"), Fraction("0.6"), Fraction("0.7")),
        "outer 6": (Fraction("0.5"), Fraction("0.7"), Fraction("0.8")),
        "inner 3 (ESC)": tuple(nodes.ESC_nodes(8, dps=40)[4:7].tolist()),
        "inner 4": (Fraction("0.03"), Fraction("0.21"), Fraction("0.68")),
        "inner 5": (Fraction("0.26"), Fraction("0.32"), Fraction("0.76")),
        "inner 6": (Fraction("0.03"), Fraction("0.47"), Fraction("0.68")),
        "inner 7": (Fraction("0.31"), Fraction("0.77"), Fraction("0.88")),
    }
    histories = {}
    # print(f"{'start':>14} {'positive':>8} {'steps':>6} {'a':>12} {'b':>12} {'c':>12} {'variance':>22} {'positive':>8}")
    # for name, params in search_starts.items():
    #     initial_node = variance_plot.symmetric_nodes(8, params)
    #     history = variance_plot.cached_compass_search(8, Fraction("0.001"), Fraction("0.000001"), initial_node)
    #     histories[name] = history

    #     # 탐색은 numpy로 하고 시작 노드와 최종 노드만 Fraction으로 다시 확인한다.
    #     _, start_weights = calculate_maximum.exact_maximums_and_weights(initial_node)
    #     final_node = history[0][-1]
    #     final_maxima, final_weights = calculate_maximum.exact_maximums_and_weights(final_node)
    #     print(f"{name:>14} {str(min(start_weights) > 0):>8} {len(history[0]) - 1:>6} "
    #           f"{final_node[-4]:>12.8f} {final_node[-3]:>12.8f} {final_node[-2]:>12.8f} "
    #           f"{float(calculate_variance.exact_variance(final_maxima[1:-1])):>22.15e} {str(min(final_weights) > 0):>8}")

    # SHOW_HISTORIES로 탐색 경로의 처음 상태를 정하고, 그림 오른쪽의 "history" 체크박스로 켜고 끈다.
    SHOW_HISTORIES = True
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_volume(
        GRID_N, n=8, cache_dir=VARIANCE_CACHE, histories=histories, show_histories=SHOW_HISTORIES
    )

    plt.show()
