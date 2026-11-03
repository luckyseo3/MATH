import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

import generate_nodes as nodes
import generate_basis as basis
import draw_basis as draw
import draw_curve as curve
import evaluate_basis as eval
import generalized_Lebegue_constant as generalized_lebesgue


if __name__ == "__main__":
    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL": nodes.GL_nodes,
        "GLo": nodes.GLo_nodes,
        "Chev": nodes.Chebyshev_nodes,
        "ESC": nodes.ESC_nodes,
        "ESCLo": nodes.ESCLo_nodes,
        "uniform": nodes.uniform_nodes
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
    a = sp.symbols("a")
    adjustable_nodes = [-1, -a, a, 1]
    adjustable_basis = basis.get_New_Basis(adjustable_nodes, basis_name='adjustable')
    a_initial = 0.5
    a_scope = [0.01, 0.99, a_initial, 0.01]

    fig_a, ax_a = plt.subplots(figsize=(10, 6))
    t_values_a, basis_funcs_a, lines_a = draw.draw_basis_1val(
        adjustable_basis, ax_a, a_initial
    )
    ax_a.set_title(rf'n=4, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")
    a_slider = draw.add_basis_slider(
        fig_a, ax_a, t_values_a, basis_funcs_a, lines_a, a_initial, a_scope
    )
    fig_a._adjustable_basis_slider = a_slider

    fig_adjustable_generalized_lebesgue, ax_adjustable_generalized_lebesgue, adjustable_generalized_lebesgue_slider = (
        draw.draw_adjustable_generalized_lebesgue_function(adjustable_basis, a_initial=a_initial)
    )

    ##################################################################################################################
    # 3️⃣ 일반화 르베그 함수를 차수와 basis 종류로 관찰
    fig_generalized_lebesgue, ax_generalized_lebesgue, generalized_lebesgue_checks, generalized_lebesgue_slider = draw.draw_generalized_lebesgue_function_selector(
        basis_sources, n_initial=3, n_min=3, n_max=15, initial_basis="GL"
    )

    ##################################################################################################################
    # 4️⃣ 5개 노드 (-1, -a, 0, a, 1)의 기저와 일반화 르베그 함수.
    five_node_basis = basis.get_New_Basis([-1, -a, 0, a, 1], basis_name='adjustable')
    fig_five, ax_five = plt.subplots(figsize=(10, 6))
    t_values_five, basis_funcs_five, lines_five = draw.draw_basis_1val(
        five_node_basis, ax_five, a_initial
    )
    ax_five.set_title(rf'n=5,Basis functions for $a={a_initial:.2f}$')
    ax_five.legend(loc="upper right")
    fig_five._adjustable_basis_slider = draw.add_basis_slider(
        fig_five, ax_five, t_values_five, basis_funcs_five, lines_five,
        a_initial, a_scope
    )
    fig_five_generalized_lebesgue, ax_five_generalized_lebesgue, five_generalized_lebesgue_slider = (
        draw.draw_adjustable_generalized_lebesgue_function(five_node_basis, a_initial=a_initial)
    )

    ##################################################################################################################
    # # 5️⃣ GL 노드 개수와 대칭 노드의 a에 따른 일반화 르베그 상수.

    fig, ax, checks, slider = draw.draw_generalized_lebesgue_constant(
        basis_sources,
        n_initial=20,
        n_min=3,
        n_max=20,
        initial_basis="GL",
    )

    fig_a4_constant, ax_a4_constant = generalized_lebesgue.generalized_Lebegue_constant_withA(4)
    fig_a5_constant, ax_a5_constant = generalized_lebesgue.generalized_Lebegue_constant_withA(5)
    
    plt.show()

    
