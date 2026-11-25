"""GL type polynomials F_i의 최댓값 또는 최솟값(kind="max"/"min")과 가중치. 최댓값의 분산은 calculate_variance.py에 있다.

최댓값은 두 가지 방법으로 구한다. 후보점(±1, 다른 노드, t*)과 알고리즘은 같고 수의 종류만 다르다.
- exact_extrema_and_weights: Fraction으로 오차 없이. 정확도용 (격자 정답지, 후보 지도, 최종 확인)
- precise_extrema_and_weights: mpmath dps자리로. 속도용 (탐색)
"""
from fractions import Fraction

def exact_extrema_and_weights(node, active=False, kind="max"):
    """유리수 노드에서 F_i의 최댓값과 가중치를 Fraction으로 오차 없이 구한다.

    후보점은 구간의 양 끝 ±1, 다른 노드, t*이다. 노드가 ±1을 포함하지 않아도 된다.
    다항식 계수는 낮은 차수부터 저장한다.
    active=True이면 (M, weights, labels)를 반환한다. labels[i]는 F_i의 최댓값을 주는
    후보점 전부의 이름("±1", "tau2"~"tau{n-1}", "t*")을 정렬해 "/"로 이은 문자열이다.
    같은 위치의 후보는 하나로 세며, 좌우 대칭인 가운데 F(n이 짝수일 때
    F_{n/2})는 t와 -t의 값이 항상 같으므로 양수 쪽 이름으로 합친다.
    """
    import construct_basis_fraction as basis
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

    # 내부 F_i는 구간의 양 끝점에서 항상 0이므로 -1과 1은 같은 후보 "±1"로 본다.
    # 끝 노드가 ±1이면(endpoint_variable a = 1) tau1, tau{n}도 "±1"이 되고, a < 1이면 tau 이름을 그대로 쓴다.
    names = {value: f"tau{j + 1}" for j, value in enumerate(node)}
    names[Fraction(-1)] = names[Fraction(1)] = "±1"
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

        if kind=="max" : 
            M.append(max(values.values()))
        elif kind=="min" :
            M.append(min(values.values()))
        else : 
            raise ValueError("Enter correct parameter")
        
        if active:
            mirror = 2 * i == n
            best = {abs(point) if mirror else point for point, value in values.items() if value == M[-1]}
            labels.append("/".join(sorted({names.get(point, "t*") for point in best})))

    return (M, w, labels) if active else (M, w)

def precise_extrema_and_weights(node, dps=40, kind="max"):
    """exact_extrema_and_weights와 같은 후보에서 mpmath(dps자리)로 최댓값(kind="max") 또는 최솟값(kind="min")과 가중치를 구한다.

    dps=40이면 n=20에서도 Fraction 정확값과 15자리 이상 일치하고, Fraction보다 약 4 ~ 7배 빠르다 (n = 10 ~ 20).
    float 노드는 그 float 값 그대로(반올림 없이) mpf로 바꾼다.
    """
    import construct_basis_mpmath as basis
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

            values = [value_at(F_i, point) for point in candidates]

            if kind=="max" :
                M.append(max(values))
            elif kind=="min" :
                M.append(min(values))
            else :
                raise ValueError("Enter correct parameter")

        return M, w
