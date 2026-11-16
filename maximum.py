"""F_i 최댓값과 그 분산, 그리고 분산을 줄이는 compass search.

최댓값은 두 가지 방법으로 구한다. 후보점(±1, 다른 노드, t*)과 알고리즘은 같고 수의 종류만 다르다.
- exact_maximums_and_weights (+ exact_variance): Fraction으로 오차 없이. 정확도용 (격자 정답지, 후보 지도, 최종 확인)
- precise_maximums_and_weights (+ precise_variance_of_maximums_of_F): mpmath dps자리로. 속도용 (탐색)
"""
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

def precise_maximums_and_weights(node, dps=40):
    """exact_maximums_and_weights와 같은 후보에서 mpmath(dps자리)로 최댓값과 가중치를 구한다.

    dps=40이면 n=20에서도 Fraction 정확값과 15자리 이상 일치하고, Fraction보다 약 4 ~ 7배 빠르다 (n = 10 ~ 20).
    float 노드는 그 float 값 그대로(반올림 없이) mpf로 바꾼다.
    """
    import mpmath as mp

    with mp.workdps(dps):
        node = [mp.mpf(value) for value in node]
        n = len(node)

        def value_at(polynomial, point):
            result = mp.mpf(0)
            for coefficient in reversed(polynomial):
                result = result * point + coefficient
            return result

        weights, D, G = [], [], []
        for i in range(n):
            polynomial = [mp.mpf(1)]
            D_i = mp.mpf(1)
            for j in range(n):
                if j == i:
                    continue
                shifted = [mp.mpf(0), *polynomial]  # polynomial * t
                for k, coefficient in enumerate(polynomial):
                    shifted[k] -= coefficient * node[j]
                polynomial = shifted
                D_i *= node[i] - node[j]
            integral = [mp.mpf(0)] + [c / (D_i * (k + 1)) for k, c in enumerate(polynomial)]
            weight = value_at(integral, 1) - value_at(integral, -1)
            if weight == 0:
                raise ValueError("Quadrature weights must be nonzero")
            G_i = [c / weight for c in integral]
            G_i[0] -= value_at(G_i, -1) + mp.mpf(1) / 2
            weights.append(weight)
            D.append(D_i)
            G.append(G_i)

        G_extended = [[mp.mpf(1) / 2], *G, [-mp.mpf(1) / 2]]
        M = []
        for i in range(n + 1):
            left, right = G_extended[i], G_extended[i + 1]
            F_i = [(left[k] if k < len(left) else 0) - (right[k] if k < len(right) else 0)
                   for k in range(max(len(left), len(right)))]
            candidates = [mp.mpf(-1), mp.mpf(1)]
            if i == 0:
                candidates.extend(node[1:])
            elif i == n:
                candidates.extend(node[:-1])
            else:
                alpha = weights[i] / D[i-1]
                beta = weights[i-1] / D[i]
                candidates.extend(node[j] for j in range(n) if j not in (i-1, i))
                difference = alpha - beta
                if difference != 0:
                    star = (alpha * node[i] - beta * node[i-1]) / difference
                    if -1 <= star <= 1:
                        candidates.append(star)
            M.append(max(value_at(F_i, point) for point in candidates))
        return M, weights

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

# 변화한 노드의 weight이 0이 되는 경우 예외처리 해야함
def compass_search(n, # 노드의 크기
                   h, # 간격
                   h_limit, # 종료 조건
                   initial_node, # 초기 노드
                   variance_function=precise_variance_of_maximums_of_F # 분산 계산 함수
                   ):

    min_node = initial_node
    min_value, _ = variance_function(initial_node)
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
                var, _ = variance_function(node)
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

def compass_search_ver2(n, # 노드의 크기
                   h, # 간격
                   h_limit, # 종료 조건
                   initial_node, # 초기 노드
                   variance_function=precise_variance_of_maximums_of_F # 분산 계산 함수
                   ):

    """compass_search에 다음을 더한 버전.

    - 이동은 (쌍 번호 idx, 부호 sign)의 튜플로 나타낸다. 쌍 idx는 오른쪽 노드 n-2-idx와
      대칭인 왼쪽 노드 1+idx이며, 오른쪽 노드를 sign*h, 왼쪽 노드를 -sign*h만큼 옮긴다.
    - 1) 직전에 성공한 이동을 먼저 시도하고, 좋아지면 다른 후보를 보지 않고 채택한다.
      2) 축 방향(쌍 하나) 후보 중 가장 좋은 것을 고른다.
      3) 축 방향이 모두 실패하면 대각선(쌍 두 개를 동시에, 모든 쌍 조합) 후보 중 가장 좋은 것을 고른다.
      4) 그것도 실패하면 h를 절반으로 줄인다.
    - 직전 이동의 반대 방향(왔던 점)은 분산이 더 크다는 것을 알고 있으므로 계산하지 않는다.
    - 대각선 후보는 쌍 개수 k에 대해 C(k,2)*4개이다 (n=10: 24개, n=20: 144개).
    """
    from itertools import combinations

    k = (n - 2) // 2
    axis_moves = [((idx, sign),) for idx in range(k) for sign in (-1, 1)]
    diagonal_moves = [((i, s), (j, t)) for i, j in combinations(range(k), 2)
                      for s in (-1, 1) for t in (-1, 1)]

    def moved(node, moves):
        """moves를 적용한 노드. 순서가 깨지면 None (짝수 n의 가운데 쌍이 0을 넘는 경우도 포함)."""
        node = node[:]
        for idx, sign in moves:
            node[n-2 - idx] += sign*h  # float 입력이면 노드 이동과 경계 비교가 반올림에 영향을 받을 수 있다.
            node[1 + idx]   -= sign*h
        if all(node[i] < node[i + 1] for i in range(n - 1)):
            return node
        return None

    def reverse(moves):
        return tuple(sorted((idx, -sign) for idx, sign in moves))

    def best_of(candidates, skip):
        """candidates 중 min_value보다 작은 것 가운데 가장 좋은 (분산, 노드, 이동). 없으면 None.
        skip에 든 이동(왔던 방향, 이미 실패한 직전 이동)은 계산하지 않는다."""
        best = None
        for moves in candidates:
            if tuple(sorted(moves)) in skip:
                continue
            node = moved(min_node, moves)
            if node is None:
                continue
            var, _ = variance_function(node)
            if var < min_value and (best is None or var < best[0]):  # 분산 차이가 수치 오차 수준이면 개선 여부가 뒤바뀔 수 있다.
                best = (var, node, moves)
        return best

    min_node = initial_node[:]
    min_value, _ = variance_function(initial_node)
    last_moves = None  # 직전에 성공한 이동

    node_history = [min_node]
    min_history = [min_value]
    h_history = [h]
    history = [node_history, min_history, h_history]

    while(h > h_limit):
        skip = set()
        found = None
        if last_moves is not None:
            skip.add(reverse(last_moves))
            found = best_of([last_moves], skip)
            skip.add(tuple(sorted(last_moves)))
        if found is None:
            found = best_of(axis_moves, skip)
        if found is None:
            found = best_of(diagonal_moves, skip)

        if found is None:
            h /= 2  # float 간격을 반복해서 줄일 때도 반올림이 누적될 수 있다.
            last_moves = None
        else:
            min_value, min_node, last_moves = found

        node_history.append(min_node)
        min_history.append(min_value)
        h_history.append(h)

    return history

def compass_search_ver3(n, # 노드의 크기
                   h, # 간격
                   h_limit, # 종료 조건
                   initial_node, # 초기 노드
                   variance_function=precise_variance_of_maximums_of_F, # 분산 계산 함수
                   state=None, # 이어서 계산할 때: 이전 실행의 on_step이 받은 상태 (initial_node는 무시된다)
                   on_step=None # 반복마다 on_step(상태)를 호출한다. 상태를 저장해 두면 state로 이어서 계산할 수 있다.
                   ):

    """compass_search_ver2를 first improvement로 바꾼 버전.

    - ver2는 축 방향(또는 대각선) 후보를 모두 계산하고 가장 좋은 것으로 움직이지만,
      ver3는 분산이 작아지는 후보를 찾는 즉시 움직인다.
    - 후보를 항상 같은 순서로 보면 앞쪽 쌍만 움직이므로, 성공한 이동은 그 목록의 맨 앞으로 옮긴다
      (move-to-front). 그래서 직전에 성공한 이동이 다음에 가장 먼저 시도된다.
    - 1) 축 방향 후보를 순서대로 보고, 좋아지는 것이 있으면 바로 움직인다.
      2) 축 방향이 모두 실패하면 대각선(모든 쌍 조합) 후보를 같은 방식으로 본다.
      3) 그것도 실패하면 h를 절반으로 줄인다.
    - 직전 이동의 반대 방향(왔던 점)은 분산이 더 크다는 것을 알고 있으므로 계산하지 않는다.
    - 실패한 단계(h를 줄이기 직전)는 ver2와 마찬가지로 후보를 모두 계산해야 한다.
    - 상태는 JSON으로 저장할 수 있는 dict이다 (이동은 [[idx, sign], ...] 리스트).
      {'h', 'last_moves', 'axis_moves', 'diagonal_moves', 'node_history', 'min_history', 'h_history'}
    """
    from itertools import combinations

    k = (n - 2) // 2
    axis_moves = [((idx, sign),) for idx in range(k) for sign in (-1, 1)]
    diagonal_moves = [((i, s), (j, t)) for i, j in combinations(range(k), 2)
                      for s in (-1, 1) for t in (-1, 1)]

    def moved(node, moves):
        """moves를 적용한 노드. 순서가 깨지면 None (짝수 n의 가운데 쌍이 0을 넘는 경우도 포함)."""
        node = node[:]
        for idx, sign in moves:
            node[n-2 - idx] += sign*h  # float 입력이면 노드 이동과 경계 비교가 반올림에 영향을 받을 수 있다.
            node[1 + idx]   -= sign*h
        if all(node[i] < node[i + 1] for i in range(n - 1)):
            return node
        return None

    def reverse(moves):
        return tuple(sorted((idx, -sign) for idx, sign in moves))

    def first_of(candidates, skip):
        """candidates를 순서대로 보며 min_value보다 작아지는 첫 (분산, 노드, 이동). 없으면 None.
        찾으면 그 이동을 candidates의 맨 앞으로 옮긴다."""
        for position, moves in enumerate(candidates):
            if tuple(sorted(moves)) == skip:
                continue
            node = moved(min_node, moves)
            if node is None:
                continue
            var, _ = variance_function(node)
            if var < min_value:  # 분산 차이가 수치 오차 수준이면 개선 여부가 뒤바뀔 수 있다.
                candidates.insert(0, candidates.pop(position))
                return var, node, moves
        return None

    def to_moves(saved):
        return tuple((idx, sign) for idx, sign in saved)

    if state is None:
        min_node = initial_node[:]
        min_value, _ = variance_function(initial_node)
        last_moves = None  # 직전에 성공한 이동
        node_history = [min_node]
        min_history = [min_value]
        h_history = [h]
    else:
        h = state['h']
        last_moves = to_moves(state['last_moves']) if state['last_moves'] is not None else None
        axis_moves = [to_moves(moves) for moves in state['axis_moves']]
        diagonal_moves = [to_moves(moves) for moves in state['diagonal_moves']]
        node_history = [list(node) for node in state['node_history']]
        min_history = list(state['min_history'])
        h_history = list(state['h_history'])
        min_node = node_history[-1]
        min_value = min_history[-1]
    history = [node_history, min_history, h_history]

    while(h > h_limit):
        skip = reverse(last_moves) if last_moves is not None else None
        found = first_of(axis_moves, skip)
        if found is None:
            found = first_of(diagonal_moves, skip)

        if found is None:
            h /= 2  # float 간격을 반복해서 줄일 때도 반올림이 누적될 수 있다.
            last_moves = None
        else:
            min_value, min_node, last_moves = found

        node_history.append(min_node)
        min_history.append(min_value)
        h_history.append(h)

        if on_step is not None:
            on_step({
                'h': h,
                'last_moves': [list(move) for move in last_moves] if last_moves is not None else None,
                'axis_moves': [[list(move) for move in moves] for moves in axis_moves],
                'diagonal_moves': [[list(move) for move in moves] for moves in diagonal_moves],
                'node_history': node_history,
                'min_history': min_history,
                'h_history': h_history,
            })

    return history

