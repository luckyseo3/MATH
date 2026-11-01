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

def ESC_nodes(n): # 1종 체비셰프 노드를 스케일링한 함수
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



