import sympy as sp

# n+1 개의 조절점
x = sp.symbols('x')

# 정의된 nodes, Px, Py 배열을 입력하면 Px와 Py의 길이가 nodes로 맞는지 확인하고 Px와 Py를 묶어 리턴하는 함수
def get_control_points(nodes, Px, Py):
    n = len(nodes)
    if len(Px) != n + 1 or len(Py) != n + 1:
        raise ValueError("The number of control points must be n + 1 for both x and y.")
    else:
        return [Px, Py]
    
# curve
# 묶인 조절점 배열과 baiss의 배열을 입력 받아 곡선 배열?을 반환하는 함수
def get_curve(P, F):
    Px, Py = P[0], P[1]

    if len(Px) != len(F) or len(Py) != len(F):
        raise ValueError(
            "The number of control points must equal "
            "the number of basis functions."
        )
    
    curve_x = sp.expand(
        sum(Px[i] * F[i] for i in range(len(F)))
    )

    curve_y = sp.expand(
        sum(Py[i] * F[i] for i in range(len(F)))
    )

    curve_x_func = sp.lambdify(
        x, curve_x, modules='numpy'
    )

    curve_y_func = sp.lambdify(
        x, curve_y, modules='numpy'
    )

    return [curve_x_func, curve_y_func]
