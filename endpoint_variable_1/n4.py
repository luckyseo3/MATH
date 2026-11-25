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
        "GL":nodes.precise_GL_nodes,
        "ESC": nodes.ESC_nodes(4, dps=40),
        "case1": np.array([Fraction("-1"), Fraction("-0.1"), Fraction("0.1"), Fraction("1")], dtype=object),
        "case2": np.array([Fraction("-1"), Fraction("-0.3"), Fraction("0.3"), Fraction("1")], dtype=object),
        "case3": np.array([Fraction("-1"), Fraction("-0.5"), Fraction("0.5"), Fraction("1")], dtype=object),
        "case4": np.array([Fraction("-1"), Fraction("-0.7"), Fraction("0.7"), Fraction("1")], dtype=object),
        "case5": np.array([Fraction("-1"), Fraction("-0.9"), Fraction("0.9"), Fraction("1")], dtype=object),
    }

    ##################################################################################################################
    # 1️⃣ basis를 차수와 종류별로 나누어 관찰
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=4,
        n_min=3,
        n_max=15,
        initial_basis="GL",
    )

    ##################################################################################################################
    # 2️⃣ 4개 노드 (-1, -a, a, 1)의 기저를 a 슬라이더로 조절한다.
    a = sp.symbols("a")
    adjustable_nodes = [-1, -a, a, 1]
    adjustable_basis = basis.get_New_Basis(
        adjustable_nodes, basis_name='adjustable'
    )
    a_initial = 0.40
    a_scope = [0.1, 0.9, a_initial, 0.01]

    fig_a, ax_a, a_sliders = draw.draw_adjustable_basis(
        adjustable_nodes,
        {a: tuple(a_scope)},
    )
    ax_a.set_title(rf'n=4, Basis functions for $a={a_initial:.2f}$')
    ax_a.legend(loc="upper right")


    ##################################################################################################################
    # 5️⃣ 4개 노드 (-1, -a, a, 1)에서 a에 따른 내부 F_1~F_3 최댓값 분산
    # 결과는 variance_cache/n4/variance.json에 저장되어 같은 점은 다시 계산하지 않는다.
    variance_a_values = [sp.Rational(k, 1000) for k in range(1, 1000)]   # 0.001 ~ 0.999
    var_inner = variance_plot.variance_curve(variance_a_values, n=4, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)

    ##################################################################################################################
    # 6️⃣ 끝점을 포함하는 노드들에서 시작해 분산을 줄이는 compass search
    histories = {}
    for name in ("ESC","case1","case2","case3","case4","case5"):
        initial_node = basis_sources[name]
        initial_node[0], initial_node[-1] = -1, 1
        history = compass_search.compass_search_ver3(
            4, Fraction("0.001"), Fraction("0.000000001"), initial_node.tolist()
        )
        histories[name] = history
        node_history, min_history, h_history = history
        print(f"\n[{name}]")
        print(f"{'step':>4} {'next h':>20} {'variance':>20}  nodes")
        for step, (node, value, step_h) in enumerate(zip(node_history, min_history, h_history)):
            node_text = ", ".join(f"{float(v):.20f}" for v in node)
            print(f"{step:>4} {float(step_h):>20.15f} {float(value):>20.15e}  [{node_text}]")

        # 탐색은 mpmath로 하고 최종 노드만 Fraction 정확값으로 다시 확인한다.
        exact_maxima, _ = calculate_extremum.exact_extrema_and_weights(node_history[-1], kind="max")
        exact_variance = calculate_variance.exact_variance(exact_maxima[1:-1])
        print(f"[{name}] a={float(node_history[-1][-2]):.8f}, exact check: {float(exact_variance):.15e}")

    # 분산 곡선 위에 탐색 경로들을 색깔별로 표시한다.
    fig_variance, ax_variance = variance_plot.draw_variance_curve(
        variance_a_values, var_inner, n=4, endpoint_variable=ENDPOINT_VARIABLE, histories=histories
    )

    plt.show()
