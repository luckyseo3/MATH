import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

import generate_nodes as nodes
import generate_basis as basis
import draw_basis as draw
import maximum
import draw_variance as variance_plot


if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL":nodes.GL_nodes,
        "ESC": nodes.ESC_nodes(7),
        "case1": np.array([-1, -0.95, -0.1, 0, 0.1, 0.95, 1]),
        "case2": np.array([-1, -0.95, -0.5, 0, 0.5, 0.95, 1]),
        "case3": np.array([-1, -0.80, -0.1, 0, 0.1, 0.80, 1]),
        "case4": np.array([-1, -0.80, -0.1, 0, 0.1, 0.80, 1]),
        "case5": np.array([-1, -0.90, -0.8, 0, 0.8, 0.90, 1]),
        "case6": np.array([-1, -0.50, -0.1, 0, 0.1, 0.50, 1]),
        "ESC optimal": np.array([-1.00000000000000000000, -0.79278273580483960004, -0.42562886791263354613, 0.00000000000000000000, 0.42562886791263365716, 0.79278273580484004412, 1.00000000000000000000]),
        "case1 optimal": np.array([-1.00000000000000000000, -0.87551500000000004320, -0.28348400000000023580, 0.00000000000000000000, 0.28348400000000023580, 0.87551500000000004320, 1.00000000000000000000]),
        "case2 optimal": np.array([-1.00000000000000000000, -0.91087000000000273303, -0.49232000000000986084, 0.00000000000000000000, 0.49232000000000986084, 0.91087000000000273303, 1.00000000000000000000]),
        "case3 optimal": np.array([-1.00000000000000000000, -0.75110099999999946352, -0.23670000000000004925, 0.00000000000000000000, 0.23670000000000004925, 0.75110099999999946352, 1.00000000000000000000]),
        "case4 optimal": np.array([-1.00000000000000000000, -0.75110099999999946352, -0.23670000000000004925, 0.00000000000000000000, 0.23670000000000004925, 0.75110099999999946352, 1.00000000000000000000]),
        # "case5 optimal": np.array([-1.00000000000000000000, -0.75110099999999946352, -0.23670000000000004925, 0.00000000000000000000, 0.23670000000000004925, 0.75110099999999946352, 1.00000000000000000000]),
        # "case6 optimal": np.array([-1.00000000000000000000, -0.75110099999999946352, -0.23670000000000004925, 0.00000000000000000000, 0.23670000000000004925, 0.75110099999999946352, 1.00000000000000000000]),
    }

    ##################################################################################################################
    # 1️⃣ basis를 차수와 종류별로 나누어 관찰
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=7,
        n_min=3,
        n_max=15,
        initial_basis="GL",
    )

    ##################################################################################################################
    # 2️⃣
    a, b = sp.symbols("a b")
    adjustable_nodes = [-1, -b, -a, 0, a, b, 1]
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
    ax_a.set_title(rf'n=7, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")


    ##################################################################################################################
    # 5️⃣ 7개 노드 (-1, -b, -a, 0, a, b, 1)에서 a, b에 따른 내부 F_1~F_6 최댓값 분산 (a < b인 점만 계산)
    # 간격 1/GRID_N 격자에서 Fraction으로 오차 없이 계산한다.
    # 결과는 variance_cache/n7/N{GRID_N}/에 저장되어 중단해도 이어서 계산한다.
    GRID_N = 1000
    variance_plot.variance_plane(GRID_N, n=7)
    # 각 격자점에서 F_1~F_6의 최댓값을 주는 후보(t*, ±1, tau_j)도 정확히 구한다.
    # 후보가 바뀌는 곳이 분산이 미분 불가능할 수 있는 경계이며, 곡면 위에 흰 선으로 표시된다.
    variance_plot.active_plane(GRID_N, n=7)

    # ##################################################################################################################
    # # 6️⃣ 끝점을 포함하는 노드들에서 시작해 분산을 줄이는 compass search
    # histories = {}
    # for name in ("ESC","case1","case2","case3","case4","case5","case6"):
    #     initial_node = basis_sources[name]
    #     initial_node[0], initial_node[3], initial_node[-1] = -1, 0, 1
    #     history = maximum.compass_search(
    #         7, 0.001, 0.000001, initial_node.tolist()
    #     )
    #     histories[name] = history
    #     node_history, min_history, h_history = history
    #     print(f"\n[{name}]")
    #     print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
    #     for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
    #         node_text = ", ".join(f"{float(v):.20f}" for v in node)
    #         print(f"{step:>4} {step_h:>20.15f} {float(value):>20.15e}  [{node_text}]")

    #     # 탐색은 numpy로 하고 최종 노드만 sympy로 다시 확인한다.
    #     sympy_variance, _ = maximum.variance_of_maximums_of_F(node_history[-1])
    #     print(f"sympy check: {float(sympy_variance):.15e}")
    # (탐색 경로를 그리려면 아래 두 함수에 histories=histories를 넘긴다.)

    # 분산 곡면 위에 최댓값 후보가 바뀌는 경계를 흰 선으로 표시한다. zoom 슬라이더나 마우스 휠로 확대한다.
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_plane(GRID_N, n=7)
    # (a, b) 평면을 최댓값 후보 조합별로 색칠한 지도
    fig_regions, ax_regions = variance_plot.draw_active_regions(GRID_N, n=7)

    plt.show()
