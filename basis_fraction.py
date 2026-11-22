import math
from fractions import Fraction

##########################################
# 1. 노드 설정: n개
# 2. Lagrange basis: n개, 실제 차수 n-1
# 3. 가중치: 
# 4. Normalized Lagrange: n개
# 5. Cumulative Lagrange: n개
# 6. GL basis: n+1개, 실제 차수 n
# 7. 조절점 정의: n+1 개
##########################################

# node는 문자열 배열이다.
def lagrange_poly(i, node):

    node = [Fraction(value) for value in node] # 문자열 노드를 fracions 자료형으로 변환
    n = len(node)

    # 계수 리스트를 관리하는 방식
    coefficients  = [Fraction(1)]
    D_i = Fraction(1)

    for j in range(n):
        if j == i :
            continue

        result = [Fraction(0)] * (len(coefficients) + 1)

        # t - \tau_j 를 곱했을 때의 효과
        for k, coefficient in enumerate(coefficients):
            result[k] -= node[j] * coefficient # -\tau_j를 모두 곱하기
            result[k+1] += coefficient # 차수 하나씩 올리기

        coefficients = result

        D_i *= node[i] - node[j]

    L_i = [coeff / D_i for coeff in coefficients]

    return L_i

# node는 문자열 배열이다.
def lagrange_basis(node):
    return [lagrange_poly(i, node) for i in range(len(node))]

# 다항함수의 함숫값을 구하는 것은 a*t + b 에 t를 곱하고 다음 계수를 더하는 과정을 반복하는 것
# point도 문자열이라고 가정
#
# Fraction으로 Horner를 하면 매 단계 약분(gcd)을 하므로 느리다. 그래서 정수로 계산하고 마지막에 한 번만 약분한다.
#   계수 c_k = a_k / C  (C = 분모들의 lcm),  point = p / q  이면
#   P(p/q) = (sum_k a_k p^k q^(d-k)) / (C q^d)
#   분자는 정수 Horner로 구한다: r = a_d,  r <- r*p + a_k*q^(d-k)  (k = d-1, ..., 0)
def value_at(polynomial, point):

    point = Fraction(point)
    coefficients = [Fraction(coefficient) for coefficient in polynomial]
    if not coefficients:
        return Fraction(0)

    p, q = point.numerator, point.denominator
    C = math.lcm(*(c.denominator for c in coefficients))
    a = [c.numerator * (C // c.denominator) for c in coefficients]

    result = a[-1]
    q_power = 1
    for a_k in reversed(a[:-1]):
        q_power *= q
        result = result * p + a_k * q_power
    return Fraction(result, C * q_power)

def weights(L):

    w = []

    for L_i in L:
        # integral = [Fraction(0)] + [coeff / Fraction(k + 1) for k, coeff in enumerate(L_i)] # k = 계수의 차수
        # w_i = value_at(integral, "1") - value_at(integral, "-1")
        
        w_i = Fraction(0)
        for k, coeff in enumerate(L_i):
            if k % 2 == 0 :
                w_i += 2 * (coeff / (k+1)) 

        w.append(w_i)

    return w

# N_i = L_i / w_i 이므로 w_i = 0이면 정의되지 않는다.
def normalized_lagrange_basis(L, w):

    if any(w_i == 0 for w_i in w):
        raise ValueError("Quadrature weights must be nonzero")

    return [[coeff / w_i for coeff in L_i] for L_i, w_i in zip(L, w)]

def cumulative_lagrange_basis(N):

    G = [[Fraction(1, 2)]] # G_0 = 1/2 (상수 다항식)

    for N_i in N:
        integral = [Fraction(0)] + [coeff / Fraction(k + 1) for k, coeff in enumerate(N_i)]
        integral[0] -= value_at(integral, "-1") + Fraction(1, 2) # 상수항에만 -(적분값(-1)) - 1/2를 더한다
        G.append(integral)

    G.append([Fraction(-1, 2)]) # G_{n+1} = -1/2 (상수 다항식)

    return G

# F_i = G_i - G_{i+1}, i = 0, ..., n
# G_0 = 1/2, G_{n+1} = -1/2는 상수이고, 가운데 G_1, ..., G_n은 길이가 같은 계수 리스트이다.
def new_basis(G):
    n = len(G) - 2
    F = []

    for i in range(n + 1):
        if i == 0:  # F_0 = 1/2 - G_1
            F_i = [-coeff for coeff in G[1]]
            F_i[0] += G[0][0]

        elif i == n:  # F_n = G_n - (-1/2)
            F_i = G[n][:]
            F_i[0] -= G[n + 1][0]

        else:  # F_i = G_i - G_{i+1}
            F_i = [left - right for left, right in zip(G[i], G[i + 1])]

        F.append(F_i)

    return F

def get_new_basis(node):

    L = lagrange_basis(node)
    w = weights(L)
    N = normalized_lagrange_basis(L, w)
    G = cumulative_lagrange_basis(N)
    F = new_basis(G)

    return F