import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

import generate_nodes as nodes
import generate_basis as basis
import draw_basis as draw
import generalized_Lebegue_constant as generalized_lebesgue
import maximum
import draw_variance as variance_plot


if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL": nodes.GL_nodes,
        "GLo": nodes.GLo_nodes,
        "Chev": nodes.Chebyshev_nodes,
        "ESC": nodes.ESC_nodes(6),
        "ESCLo": nodes.ESCLo_nodes,
        "uniform": nodes.uniform_nodes,
        "case1": np.array([-1.0, -0.2, -0.1, 0.1, 0.2, 1.0]),
        "case2": np.array([-1.0, -0.8, -0.6, 0.6, 0.8, 1.0]),
        "case3": np.array([-1.0, -0.9, -0.3, 0.3, 0.9, 1.0]),
        "case3 optimize": np.array([-1.00000000000000000000, -0.87528199999999956038, -0.32611000000001294508, 0.32611000000001294508, 0.87528199999999956038, 1.00000000000000000000]),
        "ESC optimize": np.array([-1.0, -0.72591380756887968762, -0.25409619243112524511, 0.25409619243112535614, 0.72591380756887979864, 1.00000000000000000000])
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
    # 결과는 variance_cache/n6/variance.json에 저장되어 같은 점은 다시 계산하지 않는다.
    variance_a_values = [sp.Rational(k, 100) for k in range(1, 98)]   # 0.01 ~ 0.97
    variance_b_values = [sp.Rational(k, 100) for k in range(2, 99)]  # 0.02 ~ 0.98
    var_inner = variance_plot.variance_grid(variance_a_values, variance_b_values)

    ##################################################################################################################
    # 6️⃣ 끝점을 포함하는 노드들에서 시작해 분산을 줄이는 compass search
    histories = {}
    for name in ("case1", "case2", "case3", "ESC"):
        history = maximum.compass_search(
            6, 0.001, 0.000001, basis_sources[name].tolist()
        )
        histories[name] = history
        node_history, min_history, h_history = history
        # print(f"\n[{name}]")
        # print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
        # for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
        #     node_text = ", ".join(f"{float(v):.20f}" for v in node)
        #     print(f"{step:>4} {step_h:>20.15f} {float(value):>20.15e}  [{node_text}]")

        # 탐색은 numpy로 하고 최종 노드만 sympy로 다시 확인한다.
        sympy_variance, _ = maximum.variance_of_maximums_of_F(node_history[-1])
        print(f"sympy check: {float(sympy_variance):.15e}")

    # 분산 곡면 위에 탐색 경로들을 색깔별로 표시한다.
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_surface(
        variance_a_values, variance_b_values, var_inner, histories=histories
    )

    plt.show()
