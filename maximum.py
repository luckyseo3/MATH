"""F_i 최댓값과 그 분산. 분산을 줄이는 탐색은 s_compass_search.py, s_gradient_descent.py에 있다.

최댓값은 두 가지 방법으로 구한다. 후보점(±1, 다른 노드, t*)과 알고리즘은 같고 수의 종류만 다르다.
- exact_maximums_and_weights (+ exact_variance): Fraction으로 오차 없이. 정확도용 (격자 정답지, 후보 지도, 최종 확인)
- precise_maximums_and_weights (+ precise_variance_of_maximums_of_F): mpmath dps자리로. 속도용 (탐색)
"""
import math
from fractions import Fraction

def exact_maximums_and_weights(node, active=False):
    """유리수 노드에서 F_i의 최댓값과 가중치를 Fraction으로 오차 없이 구한다.

    후보점은 ±1, 다른 노드, t*이다.
    다항식 계수는 낮은 차수부터 저장한다.
    active=True이면 (M, weights, labels)를 반환한다. labels[i]는 F_i의 최댓값을 주는
    후보점 전부의 이름("±1", "tau2"~"tau{n-1}", "t*")을 정렬해 "/"로 이은 문자열이다.
    같은 위치의 후보는 하나로 세며, 좌우 대칭인 가운데 F(n이 짝수일 때
    F_{n/2})는 t와 -t의 값이 항상 같으므로 양수 쪽 이름으로 합친다.
    """
    import basis_fraction as basis
    value_at = basis.value_at

    node = [Fraction(value) for value in node]  # float 입력은 원래 의도한 십진수가 아닌 이미 반올림된 값을 정확히 변환한다.
    n = len(node)

    L = basis.lagrange_basis(node)
    w = basis.weights(L)
    N = basis.normalized_lagrange_basis(L, w)
    G = basis.cumulative_lagrange_basis(N)
    F = basis.new_basis(G)

    # L_i의 최고차항 계수는 1 / D_i (D_i = prod_{j != i} (tau_i - tau_j))이다.
    D = [Fraction(1) / L_i[-1] for L_i in L]

    # 내부 F_i는 양 끝점에서 항상 0이므로 -1과 1은 같은 후보 "±1"로 본다.
    names = {value: f"tau{j + 1}" for j, value in enumerate(node)}
    names[node[0]] = names[node[-1]] = "±1"
    M = []
    labels = []

    for i, F_i in enumerate(F):
        # 구간의 양 끝점도 전역 최댓값의 후보이다.
        candidates = [Fraction(-1), Fraction(1)]

        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            candidates.extend(node[1:])

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            candidates.extend(node[:-1])

        else:
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            alpha = w[i] / D[i-1]
            beta = w[i-1] / D[i]
            candidates.extend(node[j] for j in range(n) if j not in (i-1, i))

            difference = alpha - beta
            if difference != 0:
                star = (alpha * node[i] - beta * node[i-1]) / difference
                if -1 <= star <= 1:
                    candidates.append(star)

        values = {point: value_at(F_i, point) for point in candidates}
        M.append(max(values.values())) # 이 부분을 min으로 바꾸면 최솟값으로 볼 수 있음
        if active:
            mirror = 2 * i == n
            best = {abs(point) if mirror else point for point, value in values.items() if value == M[-1]}
            labels.append("/".join(sorted({names.get(point, "t*") for point in best})))

    return (M, w, labels) if active else (M, w)

def exact_variance(values):
    """Return the population variance of Fractions without rounding."""
    values = [Fraction(value) for value in values]
    n = len(values)

    # 공통분모 D = n * lcm(분모들). 분모의 곱보다 자릿수가 작고(실제 최댓값에서 약 절반),
    # n을 곱해 두었으므로 sum(numer)가 n으로 나누어떨어져 평균도 정수가 된다.
    common_denom = n * math.lcm(*(v.denominator for v in values))
    numer = [v.numerator * (common_denom // v.denominator) for v in values]  # D * v_i (정수)
    mean_numer = sum(numer) // n  # D * mean (정수, 나머지 없음)

    # Var = sum((D v_i - D mean)^2) / (n D^2). 끝까지 정수 연산이고 마지막에 한 번만 약분한다.
    return Fraction(sum((x - mean_numer) ** 2 for x in numer), n * common_denom**2)

def precise_maximums_and_weights(node, dps=40):
    """exact_maximums_and_weights와 같은 후보에서 mpmath(dps자리)로 최댓값과 가중치를 구한다.

    dps=40이면 n=20에서도 Fraction 정확값과 15자리 이상 일치하고, Fraction보다 약 4 ~ 7배 빠르다 (n = 10 ~ 20).
    float 노드는 그 float 값 그대로(반올림 없이) mpf로 바꾼다.
    """
    import basis_mpmath as basis
    mp = basis.mp
    value_at = basis.value_at

    with mp.workdps(dps):
        node = [mp.mpf(value.numerator) / value.denominator
                if isinstance(value, Fraction) else mp.mpf(value) for value in node]
        n = len(node)

        L = basis.lagrange_basis(node)
        w = basis.weights(L)
        N = basis.normalized_lagrange_basis(L, w)
        G = basis.cumulative_lagrange_basis(N)
        F = basis.new_basis(G)

        # L_i의 최고차항 계수는 1 / D_i (D_i = prod_{j != i} (tau_i - tau_j))이다.
        D = [1 / L_i[-1] for L_i in L]

        M = []

        for i, F_i in enumerate(F):
            # 구간의 양 끝점도 전역 최댓값의 후보이다.
            candidates = [mp.mpf(-1), mp.mpf(1)]

            if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
                candidates.extend(node[1:])

            elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
                candidates.extend(node[:-1])

            else:
                # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
                alpha = w[i] / D[i-1]
                beta = w[i-1] / D[i]
                candidates.extend(node[j] for j in range(n) if j not in (i-1, i))

                difference = alpha - beta
                if difference != 0:
                    star = (alpha * node[i] - beta * node[i-1]) / difference
                    if -1 <= star <= 1:
                        candidates.append(star)

            M.append(max(value_at(F_i, point) for point in candidates))

        return M, w

def precise_variance_of_maximums_of_F(node, dps=40):
    """precise_maximums_and_weights로 구한 분산 (내부, 끝점 포함)을 float로 반환한다.

    분산 자체를 dps자리로 계산한 뒤 float로 바꾸므로 상대오차는 약 1e-16이다.
    """
    import mpmath as mp

    with mp.workdps(dps):
        M, _ = precise_maximums_and_weights(node, dps)

        def variance(values):
            mean = sum(values) / len(values)
            return sum((value - mean) ** 2 for value in values) / len(values)

        return float(variance(M[1:-1])), float(variance(M))
