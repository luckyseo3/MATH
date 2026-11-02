import sympy as sp

##########################################
# 1. 노드 설정: n개
# 2. Lagrange basis: n개, 실제 차수 n-1
# 3. 가중치: 
# 4. Normalized Lagrange: n개
# 5. Cumulative Lagrange: n개
# 6. GL basis: n+1개, 실제 차수 n
# 7. 조절점 정의: n+1 개
##########################################

# 0. 전역 심볼 정의
x, a, b = sp.symbols('x a b')

# 1~2. 라그랑주 기저 L_i^n
def lagrange_poly(i, nodes):
    L = sp.Integer(1)
    for j in range(len(nodes)):
        if i != j:
            L *= (x - nodes[j]) / (nodes[i] - nodes[j])
    return sp.expand(L)

def lagrange_basis(nodes):
    return [lagrange_poly(i, nodes) for i in range(len(nodes))]

# 3. 가중치 w_i = integral L_i
def weights(basis):
    return [sp.integrate(L_i,(x,-1,1)) for L_i in basis]

# 4. Normalized Lagrange
# N_i = L_i / w_i
def Normalized_Lagrange(basis, w):
    return [sp.simplify(basis[i] / w[i]) for i in range(len(basis))]

# 5. Cumulative Lagrange
# G_i(x) = integral_{-1}^{x} N_i(s) ds - 1/2

def Cumulative_Lagrange(N):
    s = sp.symbols('s')
    G = []
    for N_i in N:
        integrand = N_i.subs(x, s) #x를 s로 대체
        G_i = sp.integrate(
            integrand, 
            (s, -1, x) # s에 대한 함수를 -1부터 x까지 적분
        ) - sp.Rational(1, 2)
        G.append(sp.expand(G_i))
    return G

# 6. GL basis F_i = G_i - G_{i+1}
#
# 논문의 인덱스:
# G_0 = 1/2
# G_1, ..., G_n = 위에서 계산한 G[0], ..., G[n-1]
# G_{n+1} = -1/2
def New_Basis(G):

    G_extended = [
        sp.Rational(1, 2),
        *G,
        -sp.Rational(1, 2)
    ]

    F = []

    for i in range(len(G)+1):
        F_i = sp.expand(G_extended[i] - G_extended[i + 1])
        F.append(F_i)
    return F

# 노드를 넣으면 GL과 비슷하게 basis를 구성하는 함수
def get_New_Basis(nodes, *, basis_name=None):
    import basis_cache

    nodes = list(nodes)
    try:
        identity = basis_cache.cache_identity(
            nodes, (lagrange_poly, lagrange_basis, weights, Normalized_Lagrange,
                    Cumulative_Lagrange, New_Basis, get_New_Basis)
        )
    except (OSError, TypeError, ValueError):
        identity = None
    if identity is not None:
        cached = basis_cache.load(identity, basis_name)
        if cached is not None:
            return cached
        
    lag = lagrange_basis(nodes)
    w = weights(lag)
    # print(w)
    N = Normalized_Lagrange(lag, w)
    G = Cumulative_Lagrange(N)
    F = New_Basis(G)
    
    if identity is not None:
        basis_cache.save(identity, F, basis_name)
    return F
