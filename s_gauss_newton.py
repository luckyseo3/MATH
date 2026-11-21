import basis_fraction
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
        A[col], A[pivot] = A[pivot], A[col]

        p = A[col][col]
        A[col] = [a / p for a in A[col]]
        for r in range(m):
            if r != col and A[r][col] != 0:
                factor = A[r][col]
                A[r] = [a - factor * b for a, b in zip(A[r], A[col])]

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
    var = maximum.exact_variance(M[1:-1])

    # compass_search와 같은 형식: history[0][k], history[1][k]는 k번째 반복 후의 노드와 분산 (k = 0은 시작점)
    node_history = [node]
    variance_history = [var]
    history = [node_history, variance_history]

    while(var > variance_tolerance):

        # 1. JP = f(X_{n}) : mpmath(dps자리)로 계산
        with mp.workdps(dps):
            # 1.1. J 구하기
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
                var_new = maximum.exact_variance(M_new[1:-1])
                if var_new < var:
                    break
                reason = "variance did not decrease"

            lam /= 2
            if lam < Fraction(1, 2**20):
                raise RuntimeError(f"no acceptable step (last reason: {reason})")

        node, M, var = node_new, M_new, var_new
        node_history.append(node)
        variance_history.append(var)
        if callback is not None:
            callback(history)

    return history


