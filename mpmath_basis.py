"""mpmath 기반 기저 계산. 계수는 상수항부터 오름차순으로 저장한다.

문자열 노드를 권장한다. Fraction 버전과 달리 지정 정밀도의 근사값을 반환한다.
전체 계산: get_new_basis(node, dps=50)
개별 함수 계산: with mp.workdps(50): ...
"""

from mpmath import mp

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

    node = [mp.mpf(value) for value in node]
    n = len(node)

    # 계수 리스트를 관리하는 방식
    coefficients  = [mp.mpf(1)]
    D_i = mp.mpf(1)

    for j in range(n):
        if j == i :
            continue

        result = [mp.mpf(0)] * (len(coefficients) + 1)

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

    node = [mp.mpf(value) for value in node]
    if not node:
        raise ValueError("At least one node is required")
    if not all(mp.isfinite(value) for value in node):
        raise ValueError("Nodes must be finite")
    if len(set(node)) != len(node):
        raise ValueError("Nodes must be distinct at the working precision")
    
    return [lagrange_poly(i, node) for i in range(len(node))]

# 다항함수의 함숫값을 구하는 것은 a*t + b 에 t를 곱하고 다음 계수를 더하는 과정을 반복하는 것
# point도 문자열이라고 가정
def value_at(polynomial, point):

    point = mp.mpf(point)
    result = mp.mpf(0)
    for coefficient in reversed(polynomial):
        result = result * point + coefficient
    return result

def weights(L):

    w = []

    for L_i in L:
        # [-1, 1]에서 홀수 차수 항의 적분은 0이다.
        w_i = 2 * mp.fsum(coeff / (k + 1) for k, coeff in enumerate(L_i) if k % 2 == 0)
        w.append(w_i)

    return w

# N_i = L_i / w_i 이므로 w_i = 0이면 정의되지 않는다.
def normalized_lagrange_basis(L, w):

    if any(w_i == 0 for w_i in w):
        raise ValueError("Quadrature weights must be nonzero")

    return [[coeff / w_i for coeff in L_i] for L_i, w_i in zip(L, w)]

def cumulative_lagrange_basis(N):

    G = [[mp.mpf("0.5")]] # G_0 = 1/2 (상수 다항식)

    for N_i in N:
        integral = [mp.mpf(0)] + [coeff / mp.mpf(k + 1) for k, coeff in enumerate(N_i)]
        integral[0] -= value_at(integral, "-1") + mp.mpf("0.5") # 상수항에만 -(적분값(-1)) - 1/2를 더한다
        G.append(integral)

    G.append([mp.mpf("-0.5")]) # G_{n+1} = -1/2 (상수 다항식)

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

def get_new_basis(node, *, dps=None):
    """n개 노드로 n+1개 기저를 반환한다. dps=None이면 현재 mp.dps 사용.

    dps를 지정하면 입력 변환부터 전체 계산에 적용하고 종료 후 정밀도를 복원한다.
    """
    if dps is not None:
        if isinstance(dps, bool) or not isinstance(dps, int) or dps < 1:
            raise ValueError("dps must be a positive integer")
        with mp.workdps(dps):
            return get_new_basis(node)


    L = lagrange_basis(node)
    w = weights(L)
    N = normalized_lagrange_basis(L, w)
    G = cumulative_lagrange_basis(N)
    F = new_basis(G)

    return F