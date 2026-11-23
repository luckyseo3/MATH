"""끝 노드를 고정하지 않은 0-symmetric 노드 (-b, -a, 0, a, b)에서 내부 F_1~F_4 최댓값 분산을 관찰한다.

n4.py의 노드 가운데에 0을 추가한 경우이다. 끝 노드 b도 변수이므로 자유변수가 (a, b) 두 개이고,
대칭이므로 서로 다른 내부 최댓값은 M_1, M_2 두 개, 방정식(내부 최댓값이 모두 같다)은 하나이다.
그래서 분산이 0인 점은 한 점이 아니라 (a, b) 평면 위의 곡선을 이룰 수 있다.
0 < a < b <= 1인 간격 1/GRID_N 격자에서 분산을 Fraction으로 오차 없이 계산해 그린다.
b = 1인 줄은 끝점이 ±1인 기존 경우(endpoint_variable_1/n5.py)와 같다.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
ENDPOINT_VARIABLE = None  # 끝 노드를 고정하지 않는다: 끝 노드 b가 격자 변수이다.
VARIANCE_CACHE = pathlib.Path(__file__).resolve().parent / 'variance_cache'  # 이 폴더의 결과

from fractions import Fraction

import matplotlib.pyplot as plt

import draw_variance as variance_plot
import zero_variance_curve

GRID_N = 1000        # 격자 간격 1/GRID_N = 0.001
NEAR_ZERO = 1e-6     # 분산이 이보다 작은 양수 가중치 격자점을 자홍색 점으로 표시한다.


if __name__ == "__main__":

    ##################################################################################################################
    # 1️⃣ 5개 노드 (-b, -a, 0, a, b)에서 a, b에 따른 내부 F_1~F_4 최댓값 분산 (0 < a < b <= 1)
    # 결과는 variance_cache/n5/N{GRID_N}/에 a 값(행)마다 저장되어 중단해도 이어서 계산한다.
    variance_plot.variance_plane(GRID_N, n=5, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)
    # 각 격자점에서 F_1~F_4의 최댓값을 주는 후보(t*, ±1, tau_j)도 정확히 구한다.
    variance_plot.active_plane(GRID_N, n=5, endpoint_variable=ENDPOINT_VARIABLE, cache_dir=VARIANCE_CACHE)

    ##################################################################################################################
    # 2️⃣ 분산이 0인 곡선: 대칭이므로 분산 = 0  <=>  M_1 = M_2.
    # 격자의 행 a = i/GRID_N마다 M_1 - M_2의 부호가 바뀌는 곳에서 b를 30자리까지 이분법으로 풀고,
    # 그 (a, b)에서 분산을 Fraction으로 정확히 다시 계산한다. 결과는 variance_cache/n5/N{GRID_N}/zero_curve.json.
    curve_points = zero_variance_curve.zero_variance_curve(GRID_N, n=5, cache_dir=VARIANCE_CACHE)
    curve = [(float(Fraction(p['a'])), float(Fraction(p['b']))) for p in curve_points]
    # 왼쪽: 곡선 b(a)를 최댓값 후보 조합별 색으로, 오른쪽: 곡선을 따라간 공통 최댓값 M_1 = M_2
    fig_curve, ax_curve = zero_variance_curve.draw_zero_curve(curve_points, n=5)

    ##################################################################################################################
    # 3️⃣ log10(분산) 3차원 곡면
    # - 자홍색 점: 분산 < NEAR_ZERO (양수 가중치)와 격자 국소 최솟값
    # - 흰 점: 음수 가중치 점
    # - 흰 선: 최댓값 후보가 바뀌는 경계 (분산이 미분 불가능할 수 있는 곳), "boundaries" 체크박스로 켜고 끈다.
    # zoom 슬라이더나 마우스 휠로 확대하고, 마우스 드래그 또는 elev/azim 슬라이더로 돌려 본다.
    fig_variance, ax_variance, variance_sliders = variance_plot.draw_variance_plane(
        GRID_N, n=5, cache_dir=VARIANCE_CACHE, show_boundaries=True,
        negative_color="white", near_zero=NEAR_ZERO,
    )

    ##################################################################################################################
    # 4️⃣ 최댓값 출현 지도: (a, b) 평면을 F_1~F_4의 최댓값이 나타나는 후보(t*, ±1, tau_j) 조합별로 색칠한다.
    # 흰 선은 후보가 바뀌는 경계, 흐린 영역은 음수 가중치, 검은 실선은 분산이 0인 곡선이다.
    # (검은 점선 계열은 가중치가 0인 곡선이다.)
    fig_regions, ax_regions = variance_plot.draw_active_regions(
        GRID_N, n=5, cache_dir=VARIANCE_CACHE, show_boundaries=True, near_zero=NEAR_ZERO, curve=curve,
    )

    plt.show()
