import numpy as np
import sympy as sp

# 모든 함수는 "차수 n" 입력, "nodes 배열" 반환을 규칙으로 한다.

def GLo_nodes(n):
    # P'(n)의 근의 배열을 만둚
    interior_nodes = np.polynomial.legendre.Legendre.basis(n-1).deriv().roots()
    # 근의 배열에 양 끝점 -1, 1을 추가함
    nodes = np.concatenate(([-1.0], interior_nodes, [1.0]))
    return nodes
    
def uniform_nodes(n):
    """[-1, 1]에 등간격 노드 n개를 생성한다(간격은 2/(n-1)).

    n <= 0이면 빈 배열, n == 1이면 중앙 노드 [0]을 반환한다.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)):
        raise ValueError("n must be an integer.")
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([0.0])
    return np.linspace(-1.0, 1.0, n)


def GL_nodes(n):
    nodes, _ = np.polynomial.legendre.leggauss(n)
    return nodes

def Chebyshev_nodes(n): # 2종 체비셰프 노드를 만드는 함수
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([0.0])

    indices = np.arange(n)
    nodes = np.cos(np.pi * indices / (n - 1))
    return nodes[::-1]

def ESC_nodes(n, *, dps=None): # 1종 체비셰프 노드를 스케일링한 함수
    # 정확 계산용: float를 거치지 않고 고정밀 십진 문자열로 반환한다.
    if dps is not None:
        from mpmath import mp
        with mp.workdps(dps + 10):
            if n <= 0:
                return np.array([], dtype=object)
            if n == 1:
                return np.array(["0"], dtype=object)
            scale = mp.cos(mp.pi / (2 * n))
            values = [mp.nstr(mp.cos((2 * i + 1) * mp.pi / (2 * n)) / scale, dps)
                      for i in reversed(range(n))]
            values[0], values[-1] = "-1", "1"
            if n % 2:
                values[n // 2] = "0"
            from fractions import Fraction
            return np.array([Fraction(value) for value in values], dtype=object)
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([0.0])

    indices = np.arange(n)
    nodes = np.cos((2 * indices + 1) * np.pi / (2 * n))
    right_end = np.cos(np.pi / (2 * n))
    return (nodes / right_end)[::-1]

def Chebyshev_1st_nodes(n): # 1종 체비셰프 노드 함수
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([0.0])

    indices = np.arange(n)
    nodes = np.cos((2 * indices + 1) * np.pi / (2 * n))
    # right_end = np.cos(np.pi / (2 * n))
    return (nodes)[::-1]

def ESCLo_nodes(n):
    """2종 노드 n+2개에서 양 끝점을 제거하고 확대해 n개 노드를 반환한다.

    남은 오른쪽 끝점 cos(pi / (n+1))로 나누어 양 끝을 ±1로 맞춘다.
    n은 최종 노드 개수다. n <= 0이면 빈 배열, n == 1이면 [0]을 반환한다.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)):
        raise ValueError("n must be an integer.")
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([0.0])

    interior_nodes = Chebyshev_nodes(n + 2)[1:-1]
    return interior_nodes / interior_nodes[-1]


def crossing_legendre(n):
    if n%2 == 0:
        return False

    m = int(n/2)
    nodes = GL_nodes(m)
    cross_nodes = []

    cross_nodes.append(-1.0)
    cross_nodes.append(nodes[0])
    for i in range(m-1):

        # 새로운 노드 구하는 법
        # 1. 평균
        cross_nodes.append((nodes[i]+nodes[i+1])/2.0)
        cross_nodes.append(nodes[i+1])
    cross_nodes.append(1.0)

    return cross_nodes





def precise_GL_nodes(n):
    """GL 노드를 40자리 십진수에서 Fraction으로 변환한다."""
    from mpmath import mp
    from fractions import Fraction
    with mp.workdps(50):
        values, _ = mp.gauss_quadrature(n, "legendre")
        return np.array([Fraction(mp.nstr(v, 40)) if i != n // 2 or n % 2 == 0
                         else Fraction(0) for i, v in enumerate(values)], dtype=object)


def precise_GLo_nodes(n):
    from mpmath import mp
    from fractions import Fraction
    if n < 2:
        raise ValueError("GLo nodes require n >= 2.")
    with mp.workdps(50):
        def derivative(x):
            return mp.diff(lambda t: mp.legendre(n - 1, t), x)
        # 기존 float 근은 초기 추정에만 사용하고, 실제 노드는 고정밀로 다시 구한다.
        roots = [mp.findroot(derivative, str(v)) for v in GLo_nodes(n)[1:-1]]
        values = [Fraction(-1), *(Fraction(mp.nstr(v, 40)) for v in roots), Fraction(1)]
        if n % 2:
            values[n // 2] = Fraction(0)
        return np.array(values, dtype=object)


def precise_uniform_nodes(n):
    from fractions import Fraction
    if n <= 0:
        return np.array([], dtype=object)
    return np.array([Fraction(0)] if n == 1 else
                    [Fraction(2 * i, n - 1) - 1 for i in range(n)], dtype=object)


def precise_Chebyshev_nodes(n):
    from mpmath import mp
    from fractions import Fraction
    if n <= 1:
        return precise_uniform_nodes(n)
    with mp.workdps(50):
        values = [Fraction(mp.nstr(mp.cos(mp.pi * i / (n - 1)), 40))
                  for i in reversed(range(n))]
        values[0], values[-1] = Fraction(-1), Fraction(1)
        if n % 2:
            values[n // 2] = Fraction(0)
        return np.array(values, dtype=object)


def precise_ESCLo_nodes(n):
    from mpmath import mp
    from fractions import Fraction
    if n <= 1:
        return precise_uniform_nodes(n)
    with mp.workdps(50):
        scale = mp.cos(mp.pi / (n + 1))
        values = [Fraction(mp.nstr(mp.cos(mp.pi * i / (n + 1)) / scale, 40))
                  for i in reversed(range(1, n + 1))]
        values[0], values[-1] = Fraction(-1), Fraction(1)
        if n % 2:
            values[n // 2] = Fraction(0)
        return np.array(values, dtype=object)
