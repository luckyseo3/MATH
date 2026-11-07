import generate_basis as basis
import numpy as np
import sympy as sp


def maximums_of_F(node):
    """Return maxima on [-1, 1] for distinct, ordered numeric nodes."""
    node = [sp.sympify(value) for value in node]
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

def numeric_weights(node):
    """Return float quadrature weights of the nodes."""
    node = np.array(node, dtype=float)
    weights = []
    for i in range(len(node)):
        others = np.delete(node, i)
        integral = np.polyint(np.poly(others) / np.prod(node[i] - others))
        weights.append(np.polyval(integral, 1) - np.polyval(integral, -1))
    return weights

def numeric_maximums_of_F(node):
    """Return float maxima on [-1, 1] with the same candidates as maximums_of_F."""
    node = np.array(node, dtype=float)
    n = len(node)

    weights = []
    D = []
    G = []  # G[i]의 계수 (높은 차수부터)
    for i in range(n):
        others = np.delete(node, i)
        D.append(np.prod(node[i] - others))
        L_i = np.poly(others) / D[i]
        integral = np.polyint(L_i)
        weights.append(np.polyval(integral, 1) - np.polyval(integral, -1))
        if abs(weights[i]) < 1e-12:
            raise ValueError("Quadrature weights must be nonzero")

        G_i = integral / weights[i]
        G_i[-1] -= np.polyval(G_i, -1) + 0.5
        G.append(G_i)

    G_extended = [np.array([0.5]), *G, np.array([-0.5])]

    M = []  # F_i의 최댓값들을 저장할 리스트

    for i in range(n + 1):
        F_i = np.polysub(G_extended[i], G_extended[i + 1])
        # 구간의 양 끝점도 전역 최댓값의 후보이다.
        candidates = [-1.0, 1.0]

        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            candidates.extend(node[1:])

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            candidates.extend(node[:-1])

        else:
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            alpha = weights[i] / D[i-1]
            beta = weights[i-1] / D[i]
            candidates.extend(np.delete(node, [i-1, i]))

            difference = alpha - beta
            if difference != 0:
                star = (alpha * node[i] - beta * node[i-1]) / difference
                if -1 <= star <= 1:
                    candidates.append(star)

        M.append(np.max(np.polyval(F_i, candidates)))

    return M

def numeric_variance_of_maximums_of_F(node):
    """Return float population variances (without, with) the endpoint functions."""
    maxima = numeric_maximums_of_F(node)
    return np.var(maxima[1:-1]), np.var(maxima)

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
                n_ = node[right] + sign*h

                if n_ <= 0 or n_ <= node[right - 1]: continue
                if n_ >= node[right + 1]: continue

                node[right] += sign*h
                node[left]  -= sign*h
                var, _ = numeric_variance_of_maximums_of_F(node)
                if(var < min_value): 
                    min_value = var
                    min_node = node[:]
                    is_coordinatewise_stationary = False

        if(is_coordinatewise_stationary):
            h /= 10
        else:
            initial_node = min_node[:]
            is_coordinatewise_stationary = True

        node_history.append(min_node)
        min_history.append(min_value)
        h_history.append(h)

    return history

