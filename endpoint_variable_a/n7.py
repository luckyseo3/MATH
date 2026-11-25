"""끝 노드를 고정하지 않은 0-symmetric 노드 (-c, -b, -a, 0, a, b, c)에서 내부 F_1~F_6 최댓값 분산을 관찰한다.

endpoint_variable_1/n9.py에서 양 끝 ±1을 뺀 경우이다. 끝 노드 c도 변수이므로 자유변수가 (a, b, c) 세 개이다.
0 < a < b < c <= 1인 간격 1/GRID_N 격자에서 분산을 Fraction으로 오차 없이 계산해 그린다.
c = 1인 면은 끝점이 ±1인 기존 경우(endpoint_variable_1/n7.py)와 같다.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
ENDPOINT_VARIABLE = None  # 끝 노드를 고정하지 않는다: 끝 노드 c가 격자 변수이다.
VARIANCE_CACHE = pathlib.Path(__file__).resolve().parent / 'variance_cache'  # 이 폴더의 결과

import matplotlib.pyplot as plt
import sympy as sp

import nodes
import draw_basis as draw
import draw_variance as variance_plot
import zero_variance_curve

GRID_N = 500         # 격자 간격 1/GRID_N = 0.002 (끝점 고정 n8.py, n9.py와 같다)


if __name__ == "__main__":

    # 생성 함수를 전달하면 선택한 방식과 n에 맞춰 기저를 새로 계산한다.
    # n개 노드로 최대 n차인 n+1개의 기저함수를 생성한다.
    basis_sources = {
        "GL": nodes.precise_GL_nodes,
        "ESC": nodes.ESC_nodes(7, dps=40),
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
    # 2️⃣ 7개 노드 (-c, -b, -a, 0, a, b, c)의 기저를 a, b, c 슬라이더로 조절한다.
    a, b, c = sp.symbols("a b c")
    adjustable_nodes = [-c, -b, -a, 0, a, b, c]
    a_initial = 0.20
    b_initial = 0.60
    c_initial = 0.90
    a_scope = [0.01, 0.99, a_initial, 0.01]
    b_scope = [0.01, 0.99, b_initial, 0.01]
    c_scope = [0.01, 1.00, c_initial, 0.01]

    fig_a, ax_a, a_sliders = draw.draw_adjustable_basis(
        adjustable_nodes,
        {a: tuple(a_scope), b: tuple(b_scope), c: tuple(c_scope)},
    )
    ax_a.set_title(r'n=7, Basis functions for nodes $(-c, -b, -a, 0, a, b, c)$')
    ax_a.legend(loc="upper right")

    ##################################################################################################################
    # 3️⃣ 노드 (-c, -b, -a, 0, a, b, c)에서 a, b, c에 따른 내부 F_1~F_6 최댓값 분산 (0 < a < b < c <= 1)
    # 유리수 격자점에서 Fraction으로 오차 없이 계산하며, 분산이 작을수록 점이 크다.
    # 결과는 variance_cache/n7/N{GRID_N}/에 a 값(슬랩)마다 저장되어 중단해도 이어서 계산한다.
    # 1/GRID_N = 0.002이면 약 2070만 점이다 (코어 7개로 약 1.7시간, 디스크 약 11GB).
    variance_plot.variance_volume(GRID_N, n=7, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)

    # 탐색 경로를 그리려면 histories={이름: compass_search의 반환값}을 넘긴다.
    # 그러면 "history" 체크박스가 생기고, show_histories=False로 처음에 숨길 수 있다.
    # 4️⃣ 최댓값 분산이 0인 곡선 (M_1 = M_2 = M_3): 식 두 개, 변수 세 개이므로 해는 곡선이다.
    # a = i/GRID_N마다 g_1 = M_1 - M_2, g_2 = M_2 - M_3의 부호가 모두 바뀌는 (b, c) 칸에서 Newton으로 풀고
    # Fraction으로 분산을 정확히 확인한다. 결과는 variance_cache/n7/N{GRID_N}/zero_curve.json에 저장된다.
    curve = zero_variance_curve.zero_variance_curve_volume(GRID_N, n=7, cache_dir=VARIANCE_CACHE)
    fig_curve, ax_curve = zero_variance_curve.draw_zero_curve_volume(curve, n=7)

    # 곡선(검은 점)을 3D 분산 그림에 겹쳐 그리며 "curve" 체크박스로 켜고 끈다.
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_volume(
        GRID_N, n=7, cache_dir=VARIANCE_CACHE, curve=curve
    )

    plt.show()
