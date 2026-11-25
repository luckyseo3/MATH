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
import calculate_extremum
import calculate_variance
import s_compass_search as compass_search
import draw_variance as variance_plot


if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL": nodes.precise_GL_nodes,
        "GLo": nodes.precise_GLo_nodes,
        "Chev": nodes.precise_Chebyshev_nodes,
        "ESC": nodes.ESC_nodes(6, dps=40),
        "ESCLo": nodes.precise_ESCLo_nodes,
        "uniform": nodes.precise_uniform_nodes,
        "case1": np.array([Fraction("-1.0"), Fraction("-0.2"), Fraction("-0.1"), Fraction("0.1"), Fraction("0.2"), Fraction("1.0")], dtype=object),
        "case2": np.array([Fraction("-1.0"), Fraction("-0.8"), Fraction("-0.6"), Fraction("0.6"), Fraction("0.8"), Fraction("1.0")], dtype=object),
        "case3": np.array([Fraction("-1.0"), Fraction("-0.9"), Fraction("-0.3"), Fraction("0.3"), Fraction("0.9"), Fraction("1.0")], dtype=object),
        "case3 optimize": np.array([Fraction("-1.00000000000000000000"), Fraction("-0.87528199999999956038"), Fraction("-0.32611000000001294508"), Fraction("0.32611000000001294508"), Fraction("0.87528199999999956038"), Fraction("1.00000000000000000000")], dtype=object),
        "ESC optimize": np.array([Fraction("-1.0"), Fraction("-0.72591380756887968762"), Fraction("-0.25409619243112524511"), Fraction("0.25409619243112535614"), Fraction("0.72591380756887979864"), Fraction("1.00000000000000000000")], dtype=object)
    }

    ##################################################################################################################
    # 1️⃣ basis를 차수와 종류별로 나누어 관찰
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=3,
        n_min=3,
        n_max=15,
        initial_basis="GL",
    )

    ##################################################################################################################
    # 2️⃣ 4개 노드 (-1, -a, a, 1)의 기저를 a 슬라이더로 조절한다.
    a, b = sp.symbols("a b")
    adjustable_nodes = [-1, -b, -a, a, b, 1]
    adjustable_basis = basis.get_New_Basis(
        adjustable_nodes, basis_name='adjustable'
    )
    a_initial = 0.40
    b_initial = 0.50
    a_scope = [0.1, 0.9, a_initial, 0.01]
    b_scope = [0.1, 0.9, b_initial, 0.01]

    fig_a, ax_a, a_sliders = draw.draw_adjustable_basis(
        adjustable_nodes,
        {a: tuple(a_scope), b: tuple(b_scope)},
    )
    ax_a.set_title(rf'n=6, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")


    ##################################################################################################################
    # 5️⃣ 6개 노드 (-1, -b, -a, a, b, 1)에서 a, b에 따른 내부 F_1~F_5 최댓값 분산 (a < b인 점만 계산)
    # 마우스 드래그 또는 elev/azim 슬라이더로 그래프를 돌려볼 수 있다.
    # 간격 1/GRID_N 격자에서 Fraction으로 오차 없이 계산한다.
    # 결과는 variance_cache/n6/N{GRID_N}/에 저장되어 중단해도 이어서 계산한다.
    GRID_N = 1000
    variance_plot.variance_plane(GRID_N, n=6, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)
    # 각 격자점에서 F_1~F_5의 최댓값을 주는 후보(t*, ±1, tau_j)도 정확히 구한다.
    # 후보가 바뀌는 곳이 분산이 미분 불가능할 수 있는 경계이며, 곡면 위에 흰 선으로 표시된다.
    variance_plot.active_plane(GRID_N, n=6, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)

    ##################################################################################################################
    # 6️⃣ 끝점을 포함하는 노드들에서 시작해 분산을 줄이는 compass search
    histories = {}
    for name in ("case1", "case2", "case3", "ESC"):
        history = compass_search.compass_search_ver3(
            6, Fraction("0.001"), Fraction("0.000001"), basis_sources[name].tolist()
        )
        histories[name] = history
        node_history, min_history, h_history = history
        # print(f"\n[{name}]")
        # print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
        # for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
        #     node_text = ", ".join(f"{float(v):.20f}" for v in node)
        #     print(f"{step:>4} {float(step_h):>20.15f} {float(value):>20.15e}  [{node_text}]")

        # 탐색은 mpmath로 하고 최종 노드만 Fraction 정확값으로 다시 확인한다.
        exact_maxima, _ = calculate_extremum.exact_extrema_and_weights(node_history[-1], kind="max")
        exact_variance = calculate_variance.exact_variance(exact_maxima[1:-1])
        print(f"exact check: {float(exact_variance):.15e}")

    # 분산 곡면 위에 탐색 경로들을 색깔별로 표시한다. zoom 슬라이더나 마우스 휠로 최솟값 근처를 확대한다.
    # 최댓값 후보가 바뀌는 경계(미분 불가능할 수 있는 곳)는 흰 선이다.
    # SHOW_BOUNDARIES, SHOW_HISTORIES로 처음 상태를 정하고, 그림의 "boundaries" / "history" 체크박스로 켜고 끈다.
    SHOW_BOUNDARIES = True
    SHOW_HISTORIES = True
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_plane(
        GRID_N, n=6, cache_dir=VARIANCE_CACHE, histories=histories, show_boundaries=SHOW_BOUNDARIES, show_histories=SHOW_HISTORIES
    )
    # (a, b) 평면을 최댓값 후보 조합별로 색칠한 지도 (오른쪽 아래 "history" 체크박스)
    fig_regions, ax_regions = variance_plot.draw_active_regions(
        GRID_N, n=6, cache_dir=VARIANCE_CACHE, histories=histories, show_histories=SHOW_HISTORIES
    )

    plt.show()
