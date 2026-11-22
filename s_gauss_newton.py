import basis_fraction
import basis_mpmath
import maximum
from fractions import Fraction

# node 리스트를 입력하면 variable 리스트로 반환
def node_to_variable(node):
    n = len(node)
    return node[(n+1)//2 : -1]

# variable 리스트를 입력하면 node 리스트로 반환
def variable_to_node(variables, n):
    right = [Fraction(v) for v in variables]
    left = [-v for v in reversed(right)]
    middle = [Fraction(0)] if n%2 else []
    return [Fraction(-1)] + left + middle + right + [Fraction(1)]

# 주어진 node에서 i번째 free variable을 h만큼 움직임
def move_ith_variable(node, i, h):
    node = [Fraction(value) for value in node]
    n = len(node)

    if n%2 == 0 :
        # 짝수인 경우
        right = n // 2
        left = right - 1
    else:
        right = n // 2 + 1
        left = right - 2

    value_after_move = node[right + i] + h

    if value_after_move <= 0                   : raise ValueError("under 0")
    if node[right + i + 1] <= value_after_move : raise ValueError("break the order")
    if value_after_move <= node[right + i - 1] : raise ValueError("break the order")
    if 1 <= value_after_move                   : raise ValueError("over 1")

    node[right + i] = value_after_move
    node[left - i] = -value_after_move
    return node

# J P = f를 Fraction으로 오차 없이 푼다 (가우스-조르당 소거)
# J는 m×m (J[j][i] = ∂f_j/∂x_i), f는 길이 m. 특이행렬이면 ValueError
def solve_linear(J, f):
    m = len(f)
    A = [list(J[r]) + [f[r]] for r in range(m)]  # 확대행렬 [J | f]

    for col in range(m):
        # 정확한 수이므로 0이 아닌 피벗이면 어느 것이든 된다.
        pivot = next((r for r in range(col, m) if A[r][col] != 0), None)
        if pivot is None:
            raise ValueError(f"singular Jacobian (column {col})")
        A[col], A[pivot] = A[pivot], A[col] # pivot을 대각선에 위치하기

        p = A[col][col] # pivot값을 담을 변수 p
        A[col] = [a / p for a in A[col]] # pivot으로 열을 나누기
        for r in range(m):
            if r != col and A[r][col] != 0:
                factor = A[r][col]
                A[r] = [a - factor * b for a, b in zip(A[r], A[col])] # 열의 나머지 원소를 0으로 만들기

    return [A[r][m] for r in range(m)]

# node가 T+에 속하는지 검사: 순서, 가중치 > 0, 내부 F_i의 argmax가 모두 t*
# (M, reason)을 반환한다.
# - 속하면 M = 최댓값 리스트, reason = None
# - 속하지 않으면 M = None, reason = 실패 이유
# "t*/tau4"처럼 동점인 경우는 영역 경계(꺾이는 점)이므로 실패로 본다.
def in_T_plus(node):
    if not all(a < b for a, b in zip(node, node[1:])):
        return None, "order"
    try:
        M, w, labels = maximum.exact_maximums_and_weights(node, active=True)
    except (ValueError, ZeroDivisionError):
        return None, "weight = 0"
    if min(w) <= 0:
        return None, "weight <= 0"
    bad = [(i, label) for i, label in enumerate(labels[1:-1], start=1) if label != "t*"]
    if bad:
        return None, bad  # 예: [(2, 'tau4')] → F_2의 최댓값이 tau4에서 나옴
    return M, None

# ---------------------------------------------------------------
# 해석적 J (md "Gradient 구하기"). 인덱스는 Python 기준(0부터)이다.
# ---------------------------------------------------------------

# C_ik(t) = prod_{j != i, k} (t - τ_j) / (τ_i - τ_j)  의 계수 리스트 (상수항부터, 길이 n-1)
# C_ik는 τ_k에 무관하고  L_i(t) = C_ik(t) · (t - τ_k) / (τ_i - τ_k)  이다.
# 따라서 L_i를 (t - τ_k)로 합성 나눗셈한 몫에 (τ_i - τ_k)를 곱하면 된다. → O(n)
# (일차식 n-2개를 처음부터 곱하면 O(n^2)이다.)
# L_i(τ_k) = 0이므로 나머지는 0이다. Fraction이면 정확히 0, mpmath면 반올림 오차 수준이다.
# L_i, tau_i, tau_k는 같은 종류의 수(Fraction 또는 mpf)이다.
def C_polynomial(L_i, tau_i, tau_k):
    quotient = [0] * (len(L_i) - 1)

    # 합성 나눗셈: 최고차항부터 내려오며 r = c_d + τ_k · r
    r = L_i[-1]
    for d in range(len(L_i) - 2, -1, -1):
        quotient[d] = r
        r = L_i[d] + tau_k * r
    # 여기서 r은 나머지 L_i(τ_k)이다.

    return [q * (tau_i - tau_k) for q in quotient]

# 상수항이 0인 부정적분의 계수 리스트
def integral_polynomial(p):
    zero = p[0] * 0
    return [zero] + [c / (d + 1) for d, c in enumerate(p)]

# Φ_jk(t) = ∫_{-1}^t ∂L~_j/∂τ_k (s) ds  를 points의 각 t에서 계산한다. (k != j)
#   Φ_jk(t) = (A · B(t) - B · A(t)) / (B - τ_k A)^2
#   A(t) = ∫_{-1}^t C_jk,  B(t) = ∫_{-1}^t s C_jk,  A = A(1),  B = B(1)
# (k = j이면 L~_j가 τ_j에 무관하므로 Φ_jj = 0이다.)
def Phi(L_j, tau_j, tau_k, points):
    value_at = basis_mpmath.value_at
    C = C_polynomial(L_j, tau_j, tau_k)
    sC = [C[0] * 0] + C               # s · C_jk(s)
    IA = integral_polynomial(C)
    IB = integral_polynomial(sC)

    A_m1, B_m1 = value_at(IA, -1), value_at(IB, -1)
    A = value_at(IA, 1) - A_m1
    B = value_at(IB, 1) - B_m1
    denominator = (B - tau_k * A) ** 2   # = ((τ_j - τ_k) ω_j)^2 != 0

    values = []
    for t in points:
        A_t = value_at(IA, t) - A_m1
        B_t = value_at(IB, t) - B_m1
        values.append((A * B_t - B * A_t) / denominator)
    return values

# 해석적 J (mp.workdps(dps) 안에서 호출). J[r, x] = ∂f_r/∂x_x,  f_r = M_1 - M_{r+2}  (r, x = 0, ..., cnt_x-1)
# tau: mpf 노드, L, w: tau의 Lagrange 기저와 가중치 (basis_mpmath)
#
# 1. 포락선 정리: ∂M_i/∂τ_k = ∂F_i/∂τ_k (t*_i)  (t*_i는 미분하지 않는다)
# 2. Python 인덱스에서 F_i = G_i - G_{i+1}은 L~_{i-1} - L~_i를 적분한 것이므로
#      ∂M_i/∂τ_k = Φ_{i-1,k}(t*_i) - Φ_{i,k}(t*_i)
# 3. 자유변수 x는 τ_{right+x} = x, τ_{left-x} = -x 이므로 ∂/∂x = ∂/∂τ_{right+x} - ∂/∂τ_{left-x}
# 비용: (j, k) 쌍 O(n^2)개 × 다항식 연산 O(n) = O(n^3)
def analytic_jacobian(tau, L, w):
    from mpmath import mp
    n = len(tau)
    cnt_x = n//2 - 1
    if n%2 == 0:  # move_ith_variable과 같은 인덱스
        right = n // 2
        left = right - 1
    else:
        right = n // 2 + 1
        left = right - 2

    # 대칭이므로 서로 다른 내부 최댓값 M_1, ..., M_{cnt_x+1}의 t*만 필요하다. (maximum.py와 같은 식)
    star = {}
    for i in range(1, cnt_x + 2):
        alpha = w[i] * L[i-1][-1]
        beta = w[i-1] * L[i][-1]
        star[i] = (alpha * tau[i] - beta * tau[i-1]) / (alpha - beta)

    # 움직이는 노드 k
    K = [right + x for x in range(cnt_x)] + [left - x for x in range(cnt_x)]

    # dM[i][k] = ∂M_i/∂τ_k
    dM = {i: {k: mp.mpf(0) for k in K} for i in star}
    for j in range(cnt_x + 2):
        # L~_j는 F_{j+1}에 +로, F_j에 -로 들어간다.
        uses = [(i, sign) for i, sign in ((j + 1, 1), (j, -1)) if i in star]
        for k in K:
            if k == j:
                continue  # Φ_jj = 0
            values = Phi(L[j], tau[j], tau[k], [star[i] for i, _ in uses])
            for (i, sign), value in zip(uses, values):
                dM[i][k] += sign * value

    J = mp.matrix(cnt_x, cnt_x)
    for x in range(cnt_x):
        dM_dx = {i: dM[i][right + x] - dM[i][left - x] for i in dM}
        for r in range(cnt_x):
            J[r, x] = dM_dx[1] - dM_dx[r + 2]
    return J

# 중앙 차분 J (mp.workdps(dps) 안에서 호출). 해석적 J를 검증할 때 쓴다. 비용 O(n^4)
def finite_difference_jacobian(node, dps):
    from mpmath import mp
    cnt_x = len(node)//2 - 1
    # 중앙 차분 오차 ~ h^2 + 10^-dps / h 이므로 h = 10^(-dps*2/5) 정도로 둔다. (dps=100이면 10^-40)
    h = Fraction(1, 10**(dps * 2 // 5))
    J = mp.matrix(cnt_x, cnt_x)

    for i in range(cnt_x):
        node_plus = move_ith_variable(node, i, h)
        node_minus = move_ith_variable(node, i, -h)
        M_plus = maximum.precise_maximums_and_weights(node_plus, dps)[0]
        M_minus = maximum.precise_maximums_and_weights(node_minus, dps)[0]

        # f_j = M_1 - M_{j+2} (j = 0, ..., cnt_x-1)이므로 J[j, i] = ∂f_j/∂x_i
        for j in range(cnt_x):
            J[j, i] = ((M_plus[1] - M_plus[j + 2]) - (M_minus[1] - M_minus[j + 2])) / (2 * mp.mpf(h.numerator) / h.denominator)
    return J

# 분산 (분자, 분모) 튜플 (분모 > 0) 끼리의 비교와 float 변환
def ratio_less(a, b):
    return a[0] * b[1] < b[0] * a[1]  # a/b < c/d  <=>  a*d < c*b

def ratio_to_float(a):
    return a[0] / a[1]  # 큰 정수끼리의 / 도 정확히 반올림된 float를 준다.

def gauss_newton(
    initial_node,
    variance_tolerance,
    digits=40, # 매 걸음 뒤 자유변수를 10^-digits 격자에 반올림 (분산 하한 ~ 10^(-2*digits))
    dps=100,   # J, f, P를 계산할 mpmath 자릿수
    callback=None  # 걸음을 받아들일 때마다 callback(history)를 호출 (중간 저장용)
    ):
    # 걸음 P는 근삿값이어도 되므로 mpmath로 빠르게 구하고,
    # 걸음을 받아들일지 판정(T+ 검사, 분산 비교)과 결과 분산은 Fraction으로 정확히 계산한다.
    from mpmath import mp

    n = len(initial_node)
    cnt_x = n//2 -1  # 자유변수의 개수
    node = [Fraction(v) for v in initial_node]

    # 시작 노드가 T+ 밖이면 J가 의미 없으므로 시작하지 않는다.
    M, reason = in_T_plus(node)
    if M is None:
        raise ValueError(f"initial node is not in T+ ({reason})")
    # 분산은 약분하지 않은 (분자, 분모) 튜플로 다루고 교차곱으로 비교한다 (Fraction 약분 비용을 피한다).
    var = maximum.exact_variance(M[1:-1], gcd=False)
    tolerance = Fraction(variance_tolerance)
    tolerance = (tolerance.numerator, tolerance.denominator)

    # compass_search와 같은 형식: history[0][k], history[1][k]는 k번째 반복 후의 노드와 분산 (k = 0은 시작점)
    # 분산은 float로 기록한다. 정확한 값이 필요하면 노드에서 maximum.exact_variance로 다시 구한다.
    node_history = [node]
    variance_history = [ratio_to_float(var)]
    history = [node_history, variance_history]

    while ratio_less(tolerance, var):

        # 1. JP = f(X_{n}) : mpmath(dps자리)로 계산
        with mp.workdps(dps):
            # 1.1. J 구하기
            # 1.1.1. 현재 node의 L, w (mpmath, 반복마다 한 번: O(n^3))
            tau = [mp.mpf(v.numerator) / v.denominator for v in node]
            L = basis_mpmath.lagrange_basis(tau)
            w = basis_mpmath.weights(L)

            # 1.1.2. 해석적 J: O(n^3)  (검증할 때는 finite_difference_jacobian(node, dps)와 비교)
            J = analytic_jacobian(tau, L, w)

            # 1.2.  f(X_{n})} 구하기
            # 대칭이므로 서로 다른 내부 최댓값은 M_1, ..., M_{cnt_x+1}이다.
            # M은 직전 반복(또는 시작 전 검사)에서 Fraction으로 정확히 구한 현재 node의 최댓값이다.
            f = mp.matrix([mp.mpf((M[1] - M[j]).numerator) / (M[1] - M[j]).denominator
                           for j in range(2, cnt_x + 2)])

            # 1.3. JP = f 풀기
            P = mp.lu_solve(J, f)
            # 이후 걸음 줄이기는 Fraction으로 하므로 P를 dps자리 십진 문자열을 거쳐 Fraction으로 바꾼다.
            P = [Fraction(mp.nstr(p, dps)) for p in P]

        # 2 ~ 4. X_{n+1} = X_{n} - λP (걸음 줄이기)
        # λ = 1, 1/2, 1/4, ... 순서로 줄이며 T+에 속하고 분산이 줄어드는 첫 λ를 받아들인다.
        X = node_to_variable(node)
        lam = Fraction(1)
        while True:
            # 반올림하지 않으면 노드의 분모가 반복마다 커져 계산이 끝나지 않는다.
            X_new = [Fraction(round((x - lam * p) * 10**digits), 10**digits) for x, p in zip(X, P)]
            node_new = variable_to_node(X_new, n)

            # 3. X_{n+1} 검사 : 변화한 자유변수 노드가 T+에 속하는지 확인
            M_new, reason = in_T_plus(node_new)
            if M_new is not None:
                # 4. Var(X_{n+1})
                var_new = maximum.exact_variance(M_new[1:-1], gcd=False)
                if ratio_less(var_new, var):
                    break
                reason = "variance did not decrease"

            lam /= 2
            if lam < Fraction(1, 2**20):
                raise RuntimeError(f"no acceptable step (last reason: {reason})")

        node, M, var = node_new, M_new, var_new
        node_history.append(node)
        variance_history.append(ratio_to_float(var))
        if callback is not None:
            callback(history)

    return history


