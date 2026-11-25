"""분산을 줄이는 gradient descent (대칭 노드).

- 변수: 양수 쪽 내부 노드 x = (tau_{(n+1)//2+1}, ..., tau_{n-1}) (0과 끝 노드 a 제외, 오름차순).
  끝 노드 ±a(endpoint_variable)는 initial_node의 양 끝에서 읽고 움직이지 않는다.
  compass search처럼 대칭 쌍을 함께 움직이므로 노드는 항상 원점 대칭이다.
- 목적함수: 내부 F_1 ~ F_{n-1} 최댓값의 모분산 (calculate_extremum.precise_extrema_and_weights, mpmath dps자리).
  비교와 차분이 float 반올림(상대오차 1e-16)에 묻히지 않도록 분산을 mpf 그대로 쓴다.
- gradient: 기본은 중앙차분 (f(x+h e_j) - f(x-h e_j)) / 2h, 오차 O(h^2).
  해석적 gradient를 만들면 gradient_function=으로 바꿔 끼운다. (x, n, dps, a) -> mpf 리스트.
- step 크기: 줄어들면 받아들이기 (compass search의 "좋아지면 이동, 아니면 h를 줄이기"와 같은 생각).
  f(x - t g) < f(x) 이면 받아들이고, 아니면 t를 절반으로 줄여 다시 시도한다.
  다음 반복은 직전에 받아들인 t의 2배에서 시작한다.
- 반환값은 compass search와 같은 모양의 [node_history, variance_history, step_history]이다.
  draw_variance의 histories에 그대로 넘길 수 있다.

최적해 근처에서 분산의 Hessian 조건수는 n=10에서 약 6e3, n=15에서 약 6e4이므로
한 번의 반복마다 오차가 약 (1 - 1/조건수)배씩만 줄어든다. 분산을 아주 작게 만들려면 반복이 많이 필요하다.
"""
from fractions import Fraction

from mpmath import mp

import calculate_extremum
import calculate_variance


def full_node(x, n, a=1):
    """양수 쪽 노드 x와 끝 노드 a로 대칭 노드 전체를 만든다."""
    a = mp.mpf(a.numerator) / a.denominator if isinstance(a, Fraction) else mp.mpf(a)
    return [-a, *[-v for v in reversed(x)], *([mp.mpf(0)] if n % 2 else []), *x, a]


def half_node(node, n):
    """대칭 노드 전체에서 양수 쪽 노드 x를 꺼낸다."""
    return [mp.mpf(v) for v in node[(n + 1) // 2:-1]]


def variance(x, n, dps=40, a=1):
    """내부 최댓값의 모분산을 mpf로 반환한다."""
    with mp.workdps(dps):
        M, _ = calculate_extremum.precise_extrema_and_weights(full_node(x, n, a), dps, kind="max")
        M = M[1:-1]
        mean = mp.fsum(M) / len(M)
        return mp.fsum((value - mean) ** 2 for value in M) / len(M)


def central_difference_gradient(x, n, dps=40, a=1, h=mp.mpf('1e-12')):
    """중앙차분 gradient. dps=40에서 h=1e-12이면 절단오차 ~h^2, 반올림오차 ~1e-38/h 모두 충분히 작다."""
    with mp.workdps(dps):
        gradient = []
        for j in range(len(x)):
            plus, minus = x[:], x[:]
            plus[j] += h
            minus[j] -= h
            gradient.append((variance(plus, n, dps, a) - variance(minus, n, dps, a)) / (2 * h))
        return gradient


def is_valid(x, a=1):
    """0 < x_1 < ... < x_k < a 이어야 노드가 서로 다르고 순서가 유지된다."""
    return 0 < x[0] and all(p < q for p, q in zip(x, x[1:])) and x[-1] < a


def gradient_descent(n, # 노드의 크기
                     initial_node, # 초기 노드 (대칭, 양 끝은 ±a)
                     max_iterations=1000, # 최대 반복 횟수
                     variance_tolerance=1e-30, # 분산이 이보다 작아지면 멈춘다 (15자리 격자 노드의 바닥 근처)
                     step_limit=1e-30, # 줄어드는 step을 찾다가 step이 이보다 작아지면 멈춘다
                     initial_step=1e-3, # 첫 반복에서 시도하는 step
                     dps=40, # mpmath 자릿수
                     gradient_function=central_difference_gradient, # (x, n, dps, a) -> gradient
                     on_step=None # 반복마다 on_step(반복 번호, 노드, 분산, step)을 호출한다
                     ):

    with mp.workdps(dps):
        end = initial_node[-1]  # 끝 노드 a (endpoint_variable)
        a = mp.mpf(end.numerator) / end.denominator if isinstance(end, Fraction) else mp.mpf(end)
        x = half_node(initial_node, n)
        if not is_valid(x, a):
            raise ValueError("initial_node must be symmetric with 0 < positive nodes < a in increasing order")
        value = variance(x, n, dps, a)
        step = mp.mpf(initial_step)

        node_history = [full_node(x, n, a)]
        variance_history = [value]
        step_history = [mp.mpf(0)]
        history = [node_history, variance_history, step_history]

        for iteration in range(1, max_iterations + 1):
            if value < variance_tolerance:
                break
            gradient = gradient_function(x, n, dps, a)

            # 분산이 줄어드는 step을 찾을 때까지 절반으로 줄인다.
            while step >= step_limit:
                candidate = [v - step * g for v, g in zip(x, gradient)]
                if is_valid(candidate, a):
                    candidate_value = variance(candidate, n, dps, a)
                    if candidate_value < value:
                        break
                step /= 2
            else:
                break  # 더 줄일 수 있는 step이 없다.

            x, value = candidate, candidate_value
            node_history.append(full_node(x, n, a))
            variance_history.append(value)
            step_history.append(step)
            if on_step is not None:
                on_step(iteration, node_history[-1], value, step)

            step *= 2  # 다음 반복은 조금 더 큰 step부터 시도한다.

        return history


def to_fraction(value):
    """mpf를 같은 값의 Fraction으로 오차 없이 바꾼다."""
    from mpmath.libmp import to_rational  # man_exp는 부호를 버리므로 쓰지 않는다.

    numerator, denominator = to_rational(mp.mpf(value)._mpf_)
    return Fraction(numerator, denominator)


if __name__ == "__main__":
    import nodes

    ##################################################################################################################
    # ESC 노드에서 시작하는 gradient descent (중앙차분 gradient, 줄어들면 받아들이기)
    N = 10
    ITERATIONS = 10000000
    DPS = 40

    half = [float(v) for v in nodes.ESC_nodes(N)[(N + 1) // 2:-1]]
    initial_node = [-1.0, *[-v for v in reversed(half)], *([0.0] if N % 2 else []), *half, 1.0]

    def report(iteration, node, value, step):
        if iteration % 20 == 0:
            print(f"{iteration:>5} {float(value):>14.6e} {float(step):>10.2e}", flush=True)

    print(f"{'iter':>5} {'variance':>14} {'step':>10}")
    node_history, variance_history, step_history = gradient_descent(
        N, initial_node, max_iterations=ITERATIONS, dps=DPS, on_step=report)

    final = node_history[-1]
    print(f"\nn={N}: {len(node_history) - 1} iterations, "
          f"variance {float(variance_history[0]):.6e} -> {float(variance_history[-1]):.6e}")
    print("노드 (양수 쪽): " + ", ".join(mp.nstr(v, 20) for v in half_node(final, N)))

    # 탐색은 mpmath로 하고 최종 노드만 Fraction 정확값으로 다시 확인한다.
    exact_maxima, exact_weights = calculate_extremum.exact_extrema_and_weights([to_fraction(v) for v in final], kind="max")
    print(f"exact check: {float(calculate_variance.exact_variance(exact_maxima[1:-1])):.15e}, "
          f"min weight {float(min(exact_weights)):.4g}")
