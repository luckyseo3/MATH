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
        "ESC": nodes.ESC_nodes(5),
        "case1": np.array([-1, -0.1, 0, 0.1, 1]),
        "case2": np.array([-1, -0.3, 0, 0.3, 1]),
        "case3": np.array([-1, -0.5, 0, 0.5, 1]),
        "case4": np.array([-1, -0.7, 0, 0.7, 1]),
        "case5": np.array([-1, -0.9, 0, 0.9, 1]),
        "case1 optimize": np.array([-1.00000000000000000000, -0.36343097600000007219, 0.00000000000000000000, 0.36343097600000007219, 1.00000000000000000000]),
        "case2 optimize": np.array([-1.00000000000000000000, -0.36343097599999990566, 0.00000000000000000000, 0.36343097599999990566, 1.00000000000000000000]),
        "case3 optimize": np.array([-1.00000000000000000000, -0.61525633500000032061, 0.00000000000000000000, 0.61525633500000032061, 1.00000000000000000000]),
        "case4 optimize": np.array([-1.00000000000000000000, -0.61525633500000009857, 0.00000000000000000000, 0.61525633500000009857, 1.00000000000000000000]),
        "case5 optimize": np.array([-1.00000000000000000000, -0.82468554100000035501, 0.00000000000000000000, 0.82468554100000035501, 1.00000000000000000000]),

    }

    ##################################################################################################################
    # 1️⃣ basis를 차수와 종류별로 나누어 관찰
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=5,
        n_min=3,
        n_max=15,
        initial_basis="GL",
    )

    ##################################################################################################################
    # 2️⃣ 5개 노드 (-1, -a, 0, a, 1)의 기저를 a 슬라이더로 조절한다.
    a = sp.symbols("a")
    adjustable_nodes = [-1, -a, 0, a, 1]
    adjustable_basis = basis.get_New_Basis(
        adjustable_nodes, basis_name='adjustable'
    )
    a_initial = 0.40
    a_scope = [0.1, 0.9, a_initial, 0.01]

    fig_a, ax_a, a_sliders = draw.draw_adjustable_basis(
        adjustable_nodes,
        {a: tuple(a_scope)},
    )
    ax_a.set_title(rf'n=5, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")


    ##################################################################################################################
    # 5️⃣ 5개 노드 (-1, -a, 0, a, 1)에서 a에 따른 내부 F_1~F_4 최댓값 분산
    # 결과는 variance_cache/n5/variance.json에 저장되어 같은 점은 다시 계산하지 않는다.
    variance_a_values = [sp.Rational(k, 1000) for k in range(1, 1000)]   # 0.001 ~ 0.999
    var_inner = variance_plot.variance_curve(variance_a_values, n=5)

    ##################################################################################################################
    # 6️⃣ 끝점을 포함하는 노드들에서 시작해 분산을 줄이는 compass search
    histories = {}
    for name in ("ESC","case1","case2","case3","case4","case5"):
        initial_node = basis_sources[name]
        initial_node[0], initial_node[2], initial_node[-1] = -1, 0, 1
        history = maximum.compass_search(
            5, 0.001, 0.000000001, initial_node.tolist()
        )
        histories[name] = history
        node_history, min_history, h_history = history
        print(f"\n[{name}]")
        print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
        for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
            node_text = ", ".join(f"{float(v):.20f}" for v in node)
            print(f"{step:>4} {step_h:>20.15f} {float(value):>20.15e}  [{node_text}]")

        # 탐색은 numpy로 하고 최종 노드만 sympy로 다시 확인한다.
        sympy_variance, _ = maximum.variance_of_maximums_of_F(node_history[-1])
        print(f"[{name}] a={node_history[-1][-2]:.8f}, sympy check: {float(sympy_variance):.15e}")

    # 분산 곡선 위에 탐색 경로들을 색깔별로 표시한다.
    fig_variance, ax_variance = variance_plot.draw_variance_curve(
        variance_a_values, var_inner, n=5, histories=histories
    )

    plt.show()
