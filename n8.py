import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

import generate_nodes as nodes
import generate_basis as basis
import draw_basis as draw
import draw_variance as variance_plot

# 격자 간격은 1/GRID_N이다.
GRID_N = 500

if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL":nodes.GL_nodes,
        "ESC": nodes.ESC_nodes(8),
        "case1": np.array([-1.00000000000000000000, -0.84645906502257362014, -0.56485449735052128961, -0.20261236737965795940, 0.20261236737965809818, 0.56485449735052162268, 0.84645906502257350912, 1.00000000000000000000]),
        "case2": np.array([-1.00000000000000000000, -0.83869806502259058867, -0.54418149735056398075, -0.19062836737964597433, 0.19062836737964611311, 0.54418149735056431382, 0.83869806502259047765, 1.00000000000000000000]),
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
    variance_plot.variance_volume(GRID_N, n=8)
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_volume(GRID_N, n=8)

    plt.show()
