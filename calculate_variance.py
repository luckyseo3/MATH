"""GL type polynomials F_i 최댓값의 분산. 분산을 줄이는 탐색은 s_compass_search.py, s_newton_method.py, s_gradient_descent.py에 있다.

- exact_variance: 값 목록(예: calculate_extremum.exact_extrema_and_weights의 M[1:-1])의 모분산을 반올림 없이.
- precise_variance_of_maximums_of_F: 노드에서 calculate_extremum.precise_extrema_and_weights로 최댓값을 구해
  분산 (내부, 끝점 포함)을 mpmath dps자리로 계산한다.
"""
import math
from fractions import Fraction

from calculate_extremum import precise_extrema_and_weights

def exact_variance(values, gcd=True):
    """모분산을 반올림 없이 구한다.

    gcd=True이면 약분한 Fraction을, gcd=False이면 약분하지 않은 (분자, 분모) 정수 튜플을 반환한다.
    튜플의 분모는 항상 양수이다.
    큰 n에서는 Fraction의 약분(gcd)이 분산 계산 시간의 약 40%를 차지하지만 거의 줄지 않는다
    (n = 20에서 약 77만 비트 중 12비트). 비교만 할 때는 a/b < c/d  <=>  a*d < c*b 로 충분하다.
    """
    values = [Fraction(value) for value in values]
    n = len(values)

    # 공통분모 D = n * lcm(분모들). 분모의 곱보다 자릿수가 작고(실제 최댓값에서 약 절반),
    # n을 곱해 두었으므로 sum(numer)가 n으로 나누어떨어져 평균도 정수가 된다.
    common_denom = n * math.lcm(*(v.denominator for v in values))
    numer = [v.numerator * (common_denom // v.denominator) for v in values]  # D * v_i (정수)
    mean_numer = sum(numer) // n  # D * mean (정수, 나머지 없음)

    # Var = sum((D v_i - D mean)^2) / (n D^2). 끝까지 정수 연산이다.
    numerator = sum((x - mean_numer) ** 2 for x in numer)
    denominator = n * common_denom**2
    return Fraction(numerator, denominator) if gcd else (numerator, denominator)

def precise_variance_of_maximums_of_F(node, dps=40):
    """precise_extrema_and_weights로 구한 분산 (내부, 끝점 포함)을 float로 반환한다.

    분산 자체를 dps자리로 계산한 뒤 float로 바꾸므로 상대오차는 약 1e-16이다.
    """
    import mpmath as mp

    with mp.workdps(dps):
        M, _ = precise_extrema_and_weights(node, dps, kind="max")

        def variance(values):
            mean = sum(values) / len(values)
            return sum((value - mean) ** 2 for value in values) / len(values)

        return float(variance(M[1:-1])), float(variance(M))
