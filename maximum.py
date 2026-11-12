from fractions import Fraction

import generate_basis as basis
import numpy as np
import sympy as sp


def maximums_of_F(node):
    """Return maxima on [-1, 1] for distinct, ordered numeric nodes."""
    node = [sp.sympify(value) for value in node]  # float 입력은 이미 반올림된 값을 SymPy Float로 유지하므로 이후 연산도 근사될 수 있다.
    n = len(node)

    if n < 3:
        raise ValueError("At least three node is required")
    for value in node:
        if not value.is_number or value.is_real is not True:
            raise ValueError("Nodes must be real numbers")
        if (value < -1) is sp.true or (value > 1) is sp.true:
            raise ValueError("Nodes must be real and lie in [-1, 1]")
    for left, right in zip(node, node[1:]):
        if sp.simplify(left >= right) is sp.true:
            raise ValueError("Nodes must be distinct and strictly increasing")

    lagrange = basis.lagrange_basis(node)
    weights, _, _ = basis.weights(lagrange)
    if any(w.equals(0) is True for w in weights):
        raise ValueError("Quadrature weights must be nonzero")
    
    F = basis.New_Basis(basis.Cumulative_Lagrange(
        basis.Normalized_Lagrange(lagrange, weights)))
    
    t = basis.x

    def evaluate(expr, point):
        # 무리수 노드에서는 정확히 0인 값도 근호 식으로 남아 수치 평가로
        # 부호를 판정할 수 없으므로(sp.Max가 거부), 기호적으로 먼저 정리한다.
        return sp.simplify(expr.subs(t, point))

    M = []  # F_i의 최댓값들을 저장할 리스트

    for i, F_i in enumerate(F):
        
        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            # 구간의 양 끝점도 전역 최댓값의 후보이다.
            maximum = sp.Max(evaluate(F_i, -1), evaluate(F_i, 1))
            for j in range(1, n):
                value = evaluate(F_i, node[j])
                maximum = sp.Max(maximum, value)

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            # 구간의 양 끝점도 전역 최댓값의 후보이다.
            maximum = sp.Max(evaluate(F_i, -1), evaluate(F_i, 1))
            for j in range(0, n-1):
                value = evaluate(F_i, node[j])
                maximum = sp.Max(maximum, value)

        else:
            maximum = sp.Max(evaluate(F_i, -1), evaluate(F_i, 1))
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            D_i = sp.prod(node[i-1] - node[j] for j in range(n) if j != i-1)
            D_next = sp.prod(node[i] - node[j] for j in range(n) if j != i)
            alpha = weights[i] / D_i
            beta = weights[i-1] / D_next

            for j in range(0, n):
                if j == i-1 or j == i:
                    continue
                value = evaluate(F_i, node[j])
                maximum = sp.Max(maximum, value)

            difference = alpha - beta
            if difference.equals(0) is not True:
                star = sp.simplify(
                    (alpha * node[i] - beta * node[i-1]) / difference)
                value = evaluate(F_i, star)
                if -1 <= star <= 1:
                    maximum = sp.Max(maximum, value)
        M.append(maximum)

    return M

def variance_of_maximums_of_F(node):
    """Return population variances (without, with) the endpoint functions."""
    maxima = maximums_of_F(node)

    def variance(values):
        mean = sum(values) / len(values)
        return sp.simplify(sum((value - mean) ** 2 for value in values) / len(values))

    return variance(maxima[1:-1]), variance(maxima)

def exact_maximums_and_weights(node, active=False):
    """유리수 노드에서 F_i의 최댓값과 가중치를 Fraction으로 오차 없이 구한다.

    후보점은 maximums_of_F와 같고, sympy를 거치지 않아 격자 계산에 쓸 수 있다.
    다항식 계수는 낮은 차수부터 저장한다.
    active=True이면 (M, weights, labels)를 반환한다. labels[i]는 F_i의 최댓값을 주는
    후보점 전부의 이름("±1", "tau2"~"tau{n-1}", "t*")을 정렬해 "/"로 이은 문자열이다.
    같은 위치의 후보는 하나로 세며, 좌우 대칭인 가운데 F(n이 짝수일 때
    F_{n/2})는 t와 -t의 값이 항상 같으므로 양수 쪽 이름으로 합친다.
    """
    node = [Fraction(value) for value in node]  # float 입력은 원래 의도한 십진수가 아닌 이미 반올림된 값을 정확히 변환한다.
    n = len(node)

    def value_at(polynomial, point):
        result = Fraction(0)
        for coefficient in reversed(polynomial):
            result = result * point + coefficient
        return result

    weights = []
    D = []
    G = []
    for i in range(n):
        polynomial = [Fraction(1)]
        D_i = Fraction(1)
        for j in range(n):
            if j == i:
                continue
            shifted = [Fraction(0), *polynomial]  # polynomial * t
            for k, coefficient in enumerate(polynomial):
                shifted[k] -= coefficient * node[j]
            polynomial = shifted
            D_i *= node[i] - node[j]

        integral = [Fraction(0)] + [c / (D_i * (k + 1)) for k, c in enumerate(polynomial)]
        weight = value_at(integral, Fraction(1)) - value_at(integral, Fraction(-1))
        if weight == 0:
            raise ValueError("Quadrature weights must be nonzero")

        G_i = [c / weight for c in integral]
        G_i[0] -= value_at(G_i, Fraction(-1)) + Fraction(1, 2)
        weights.append(weight)
        D.append(D_i)
        G.append(G_i)

    G_extended = [[Fraction(1, 2)], *G, [Fraction(-1, 2)]]

    # 내부 F_i는 양 끝점에서 항상 0이므로 -1과 1은 같은 후보 "±1"로 본다.
    names = {value: f"tau{j + 1}" for j, value in enumerate(node)}
    names[node[0]] = names[node[-1]] = "±1"
    M = []
    labels = []
    for i in range(n + 1):
        left, right = G_extended[i], G_extended[i + 1]
        F_i = [(left[k] if k < len(left) else 0) - (right[k] if k < len(right) else 0)
               for k in range(max(len(left), len(right)))]
        # 구간의 양 끝점도 전역 최댓값의 후보이다.
        candidates = [Fraction(-1), Fraction(1)]

        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            candidates.extend(node[1:])

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            candidates.extend(node[:-1])

        else:
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            alpha = weights[i] / D[i-1]
            beta = weights[i-1] / D[i]
            candidates.extend(node[j] for j in range(n) if j not in (i-1, i))

            difference = alpha - beta
            if difference != 0:
                star = (alpha * node[i] - beta * node[i-1]) / difference
                if -1 <= star <= 1:
                    candidates.append(star)

        values = {point: value_at(F_i, point) for point in candidates}
        M.append(max(values.values()))
        if active:
            mirror = 2 * i == n
            best = {abs(point) if mirror else point for point, value in values.items() if value == M[-1]}
            labels.append("/".join(sorted({names.get(point, "t*") for point in best})))

    return (M, weights, labels) if active else (M, weights)

def exact_variance(values):
    """Return the population variance of Fractions without rounding."""
    mean = sum(values) / len(values)  # values가 Fraction 등 정확한 수일 때만 아래 연산에 부동소수점 오차가 없다.
    return sum((value - mean) ** 2 for value in values) / len(values)

def numeric_weights(node):
    """Return float quadrature weights of the nodes."""
    node = np.array(node, dtype=float)  # float 변환 시 입력값이 반올림될 수 있다.
    weights = []
    for i in range(len(node)):
        others = np.delete(node, i)
        integral = np.polyint(np.poly(others) / np.prod(node[i] - others))  # 계수 전개, 노드 차이의 곱, 나눗셈과 항별 적분에서 반올림이 누적된다.
        weights.append(np.polyval(integral, 1) - np.polyval(integral, -1))  # 양 끝값 평가와 뺄셈에서 오차가 생기며 가까운 값의 차이는 소거에 취약하다.
    return weights

def numeric_maximums_of_F(node):
    """Return float maxima on [-1, 1] with the same candidates as maximums_of_F."""
    node = np.array(node, dtype=float)  # 입력을 float로 변환하면서 정밀도가 손실될 수 있다.
    n = len(node)

    weights = []
    D = []
    G = []  # G[i]의 계수 (높은 차수부터)
    for i in range(n):
        others = np.delete(node, i)
        D.append(np.prod(node[i] - others))  # 가까운 노드의 차이와 그 곱은 반올림 오차에 민감하다.
        L_i = np.poly(others) / D[i]  # 다항식 계수 전개와 나눗셈에서 오차가 누적될 수 있다.
        integral = np.polyint(L_i)  # 항별 적분의 나눗셈과 양 끝값 평가 및 뺄셈에도 반올림이 발생한다.
        weights.append(np.polyval(integral, 1) - np.polyval(integral, -1))
        if abs(weights[i]) < 1e-12:
            raise ValueError("Quadrature weights must be nonzero")

        G_i = integral / weights[i]  # 작은 가중치로 나누면 기존 오차가 증폭될 수 있다.
        G_i[-1] -= np.polyval(G_i, -1) + 0.5  # 다항식 평가와 상수항 보정에서도 반올림 및 소거가 발생할 수 있다.
        G.append(G_i)

    G_extended = [np.array([0.5]), *G, np.array([-0.5])]

    M = []  # F_i의 최댓값들을 저장할 리스트

    for i in range(n + 1):
        F_i = np.polysub(G_extended[i], G_extended[i + 1])  # 비슷한 계수의 뺄셈은 소거로 상대오차를 키울 수 있다.
        # 구간의 양 끝점도 전역 최댓값의 후보이다.
        candidates = [-1.0, 1.0]

        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            candidates.extend(node[1:])

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            candidates.extend(node[:-1])

        else:
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            alpha = weights[i] / D[i-1]  # 나눗셈과 차이 계산의 오차는 임계점과 후보 판정에 영향을 준다.
            beta = weights[i-1] / D[i]
            candidates.extend(np.delete(node, [i-1, i]))

            difference = alpha - beta  # 특히 차이가 0에 가까우면 star가 매우 민감하며 정확한 0 비교도 불안정하다.
            if difference != 0:
                star = (alpha * node[i] - beta * node[i-1]) / difference
                if -1 <= star <= 1:
                    candidates.append(star)

        M.append(np.max(np.polyval(F_i, candidates)))  # 후보별 다항식 평가의 반올림으로 최댓값이나 선택된 후보가 달라질 수 있다.

    return M

def numeric_variance_of_maximums_of_F(node):
    """Return float population variances (without, with) the endpoint functions."""
    maxima = numeric_maximums_of_F(node)
    return np.var(maxima[1:-1]), np.var(maxima)  # 평균, 편차, 제곱과 합산에서 반올림되며 입력 최댓값의 오차도 전파된다.

# 변화한 노드의 weight이 0이 되는 경우 예외처리 해야함
def compass_search(n, # 노드의 크기
                   h, # 간격
                   h_limit, # 종료 조건
                   initial_node # 초기 노드
                   ):
    
    min_node = initial_node
    min_value, _ = numeric_variance_of_maximums_of_F(initial_node)
    is_coordinatewise_stationary = True

    node_history = [min_node]
    min_history = [min_value]
    h_history = [h]
    history = [node_history, min_history, h_history]

    while(h > h_limit):

        for idx in range((n-2)//2):
            left = 1 + idx
            right = n-2 - idx
            for sign in [-1, 1]:
                node = initial_node[:]
                n_ = node[right] + sign*h  # float 입력이면 노드 이동과 경계 비교가 반올림에 영향을 받을 수 있다.

                if n_ <= 0 or n_ <= node[right - 1]: continue
                if n_ >= node[right + 1]: continue

                node[right] += sign*h
                node[left]  -= sign*h
                var, _ = numeric_variance_of_maximums_of_F(node)
                if(var < min_value):   # 분산 차이가 수치 오차 수준이면 개선 여부가 뒤바뀔 수 있다.
                    min_value = var
                    min_node = node[:]
                    is_coordinatewise_stationary = False

        if(is_coordinatewise_stationary):
            h /= 10  # float 간격을 반복해서 줄일 때도 반올림이 누적될 수 있다.
        else:
            initial_node = min_node[:]
            is_coordinatewise_stationary = True

        node_history.append(min_node)
        min_history.append(min_value)
        h_history.append(h)

    return history

