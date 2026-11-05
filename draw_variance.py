"""노드 {-1, -b, -a, a, b, 1}에서 a, b에 따른 내부 F_i 최댓값 분산을 3차원으로 그린다."""
import numpy as np
import sympy as sp

import variance_cache
from variance_cache import inner_variance  # noqa: F401  (기존 호출 호환)


def variance_grid(a_values, b_values, *, n=6, workers=None):
    """(a, b) 격자에서 양 끝을 제외한 분산을 variance_cache를 거쳐 구한다.

    n=6이면 노드 [-1, -b, -a, a, b, 1], n=7이면 [-1, -b, -a, 0, a, b, 1].
    a < b가 아니거나 가중치가 0이 되어 정의되지 않는 점은 NaN으로 둔다.
    """
    if variance_cache.parameter_count(n) != 2:
        raise ValueError("3D surface needs exactly two parameters (n = 6 or 7)")
    pairs = [(a, b) for a in a_values for b in b_values]
    values = variance_cache.get_variances(n, pairs, workers=workers)
    return np.array(values, dtype=float).reshape(len(a_values), len(b_values))


def draw_variance_surface(a_values, b_values, Z, *, elev=25, azim=-130,
                          max_points=400, max_mesh=60):
    """마우스 드래그와 elev/azim 슬라이더로 돌려볼 수 있는 3차원 그래프.

    반환값은 (figure, axes, sliders)이다. 슬라이더가 가비지 컬렉션되면
    동작하지 않으므로 호출한 쪽에서 sliders를 변수로 들고 있어야 한다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    A, B = np.meshgrid([float(v) for v in a_values], [float(v) for v in b_values], indexing="ij")
    mask = ~np.isnan(Z)
    # 분산이 0.003 ~ 50 정도로 퍼져 있어 log10 축으로 그린다.
    L = np.log10(Z)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_axes([0.05, 0.18, 0.9, 0.77], projection="3d")
    # 드래그/슬라이더 회전 시 매 프레임 전체를 다시 그리므로, 촘촘한 격자에서는
    # (0.01 간격 약 4,700점에서 프레임당 0.75초) 화면이 멈춘 것처럼 보인다.
    # 점이 적으면 삼각분할 곡면과 격자점을, 많으면 최대 max_mesh x max_mesh로 줄인
    # 격자 곡면만 그려 회전을 가볍게 한다. 최솟값 표시는 항상 전체 데이터 기준이다.
    if mask.sum() <= max_points:
        ax.plot_trisurf(A[mask], B[mask], L[mask], cmap="viridis", alpha=0.85, edgecolor="none")
        ax.scatter(A[mask], B[mask], L[mask], color="black", s=10)
    else:
        # 정의되지 않은 점(a >= b)은 마스킹하고, 불투명 격자 곡면만 그린다.
        ax.plot_surface(A, B, np.ma.masked_invalid(L), cmap="viridis", edgecolor="none",
                        rcount=min(max_mesh, A.shape[0]), ccount=min(max_mesh, A.shape[1]),
                        antialiased=False)
    p, q = np.unravel_index(np.nanargmin(Z), Z.shape)
    ax.scatter([A[p, q]], [B[p, q]], [L[p, q]], color="red", s=60,
               label=f"min {Z[p, q]:.4g} at a={A[p, q]:.3g}, b={B[p, q]:.3g}")
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("log10(variance)")
    ax.set_title("Variance of max F_i (i = 1..5), nodes {-1, -b, -a, a, b, 1}, a < b")
    ax.legend(loc="upper left")
    ax.view_init(elev=elev, azim=azim)

    # 마우스로 드래그해도 회전하며, 슬라이더로 정확한 각도를 지정할 수도 있다.
    elev_slider = Slider(fig.add_axes([0.15, 0.08, 0.7, 0.03]), "elev", -90, 90, valinit=elev)
    azim_slider = Slider(fig.add_axes([0.15, 0.03, 0.7, 0.03]), "azim", -180, 180, valinit=azim)

    def update(_):
        ax.view_init(elev=elev_slider.val, azim=azim_slider.val)
        fig.canvas.draw_idle()

    elev_slider.on_changed(update)
    azim_slider.on_changed(update)

    def sync_sliders(_):
        # 마우스로 돌린 뒤 슬라이더 값도 현재 각도에 맞춘다.
        for slider, value in ((elev_slider, ax.elev), (azim_slider, ax.azim)):
            if abs(slider.val - value) > 1e-9:
                slider.eventson = False
                slider.set_val(((value + 180) % 360) - 180 if slider is azim_slider else value)
                slider.eventson = True

    fig.canvas.mpl_connect("button_release_event", sync_sliders)
    return fig, ax, (elev_slider, azim_slider)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    a_values = [sp.Rational(k, 10) for k in range(1, 9)]   # 0.1 ~ 0.8
    b_values = [sp.Rational(k, 10) for k in range(2, 10)]  # 0.2 ~ 0.9
    Z = variance_grid(a_values, b_values)

    print("a\\b " + " ".join(f"{float(b):>9.1f}" for b in b_values))
    for a, row in zip(a_values, Z):
        print(f"{float(a):.1f} " + " ".join("        -" if np.isnan(v) else f"{v:9.5f}" for v in row))

    fig, ax, sliders = draw_variance_surface(a_values, b_values, Z)
    fig.savefig("variance_surface.png", dpi=130)
    plt.show()
