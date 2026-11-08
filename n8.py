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
        "ESC": nodes.ESC_nodes(8),
        "case1": np.array([-1.00000000000000000000, -0.84645906502257362014, -0.56485449735052128961, -0.20261236737965795940, 0.20261236737965809818, 0.56485449735052162268, 0.84645906502257350912, 1.00000000000000000000]),
        "case2": np.array([-1.00000000000000000000, -0.83869806502259058867, -0.54418149735056398075, -0.19062836737964597433, 0.19062836737964611311, 0.54418149735056431382, 0.83869806502259047765, 1.00000000000000000000]),
        "case3": np.array([-1, -0.80, -0.1, 0, 0.1, 0.80, 1]),
        "case4": np.array([-1, -0.80, -0.1, 0, 0.1, 0.80, 1]),
        "case5": np.array([-1, -0.90, -0.8, 0, 0.8, 0.90, 1]),
        "case6": np.array([-1, -0.50, -0.1, 0, 0.1, 0.50, 1]),
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
    adjustable_nodes = [-1, -c, -b, -a, 0, a, b, c, 1]
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

    # initial_node = nodes.ESC_nodes(8)
    # history = maximum.compass_search(8, 0.01, 1e-6, initial_node.tolist())

    # node_history, min_history, h_history = history
    # print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
    # for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
    #     node_text = ", ".join(f"{float(v):.20f}" for v in node)
    #     print(f"{step:>4} {step_h:>20.15f} {float(value):>20.15e}  [{node_text}]")

    plt.show()