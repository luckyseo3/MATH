import generate_basis as basis
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
        if value.is_real is False or (value < -1) is sp.true or (value > 1) is sp.true:
            raise ValueError("Nodes must be real and lie in [-1, 1]")
    for left, right in zip(node, node[1:]):
        if sp.simplify(left >= right) is sp.true:
            raise ValueError("Nodes must be distinct and strictly increasing")

    lagrange = basis.lagrange_basis(node)
    weights = basis.weights(lagrange)
    if any(w.equals(0) is True for w in weights):
        raise ValueError("Quadrature weights must be nonzero")
    
    F = basis.New_Basis(basis.Cumulative_Lagrange(
        basis.Normalized_Lagrange(lagrange, weights)))
    
    t = basis.x
    
    M = []  # F_i의 최댓값들을 저장할 리스트

    for i, F_i in enumerate(F):
        
        if i == 0:  # F_0의 임계점: node[1], ..., node[n-1]
            # 구간의 양 끝점도 전역 최댓값의 후보이다.
            maximum = sp.Max(F_i.subs(t, -1), F_i.subs(t, 1))
            for j in range(1, n):
                value = F_i.subs(t, node[j])
                maximum = sp.Max(maximum, value)

        elif i == n:  # F_n의 임계점: node[0], ..., node[n-2]
            # 구간의 양 끝점도 전역 최댓값의 후보이다.
            maximum = sp.Max(F_i.subs(t, -1), F_i.subs(t, 1))
            for j in range(0, n-1):
                value = F_i.subs(t, node[j])
                maximum = sp.Max(maximum, value)

        else:
            maximum = sp.Max(F_i.subs(t, -1), F_i.subs(t, 1))
            # 문서의 i, i+1은 Python의 i-1, i에 대응한다.
            D_i = sp.prod(node[i-1] - node[j] for j in range(n) if j != i-1)
            D_next = sp.prod(node[i] - node[j] for j in range(n) if j != i)
            alpha = weights[i] / D_i
            beta = weights[i-1] / D_next

            for j in range(0, n):
                if j == i-1 or j == i:
                    continue
                value = F_i.subs(t, node[j])
                maximum = sp.Max(maximum, value)

            difference = alpha - beta
            if difference.equals(0) is not True:
                star = (alpha * node[i] - beta * node[i-1]) / difference
                value = F_i.subs(t, star)
                if -1 <= star <= 1:
                    maximum = sp.Max(maximum, value)
        M.append(maximum)

    return M

def variance_of_maximums_of_F(node):
    """Return population variances with and without the endpoint functions."""
    maxima = maximums_of_F(node)

    def variance(values):
        mean = sum(values) / len(values)
        return sp.simplify(sum((value - mean) ** 2 for value in values) / len(values))

    return variance(maxima), variance(maxima[1:-1])


