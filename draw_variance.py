"""0 대칭·끝점 포함 노드에서 매개변수에 따른 내부 F_i 최댓값 분산을 3차원으로 그린다."""
import json
from pathlib import Path

import numpy as np
import sympy as sp

import maximum


def symmetric_nodes(n, params):
    return [-1, *(-p for p in reversed(params)), *([0] if n % 2 else []), *params, 1]


def has_negative_weight(n, params):
    return min(maximum.numeric_weights(symmetric_nodes(n, params))) < 0


def load_points(n):
    path = Path(__file__).resolve().parent / 'variance_cache' / f'n{n}' / 'variance.json'
    return json.loads(path.read_text(encoding='utf-8'))['points'] if path.exists() else {}


def save_points(n, points):
    path = Path(__file__).resolve().parent / 'variance_cache' / f'n{n}' / 'variance.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'points': points}, indent=1), encoding='utf-8')


def cached_variance(n, params, points):
    """params의 분산을 points에서 읽고, 없으면 계산해 추가한다. 정의되지 않으면 NaN."""
    key = ",".join(str(p) for p in params)
    if key not in points:
        try:
            value, _ = maximum.numeric_variance_of_maximums_of_F(symmetric_nodes(n, params))
            points[key] = {'value': float(value), 'exact': None}
        except ValueError:  # 가중치가 0인 점
            points[key] = {'value': None, 'exact': None}
    value = points[key]['value']
    return np.nan if value is None else value


def variance_curve(a_values, *, n=5):
    """노드 [-1, -a, 0, a, 1](n=5) 또는 [-1, -a, a, 1](n=4)에서 a에 따른 분산을 구한다.

    0 < a < 1이 아니거나 가중치가 0이 되어 정의되지 않는 점은 NaN으로 둔다.
    """
    if n not in (4, 5):
        raise ValueError("3D curve needs exactly one parameter (n = 4 or 5)")
    points = load_points(n)
    computed = len(points)

    V = np.full(len(a_values), np.nan)
    for p, a in enumerate(a_values):
        a = sp.Rational(str(a))
        if 0 < a < 1:
            V[p] = cached_variance(n, (a,), points)

    if len(points) != computed:
        save_points(n, points)
    return V


def variance_grid(a_values, b_values, *, n=6):
    """(a, b) 격자에서 양 끝을 제외한 분산을 variance_cache/n{n}/variance.json을 거쳐 구한다.

    n=6이면 노드 [-1, -b, -a, a, b, 1], n=7이면 [-1, -b, -a, 0, a, b, 1].
    a < b가 아니거나 가중치가 0이 되어 정의되지 않는 점은 NaN으로 둔다.
    """
    if n not in (6, 7):
        raise ValueError("3D surface needs exactly two parameters (n = 6 or 7)")
    points = load_points(n)
    computed = len(points)

    Z = np.full((len(a_values), len(b_values)), np.nan)
    for p, a in enumerate(a_values):
        for q, b in enumerate(b_values):
            a, b = sp.Rational(str(a)), sp.Rational(str(b))
            if 0 < a < b < 1:
                Z[p, q] = cached_variance(n, (a, b), points)

    if len(points) != computed:
        save_points(n, points)
    return Z


def add_view_sliders(fig, ax, elev, azim):
    from matplotlib.widgets import Slider

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
    return elev_slider, azim_slider


def draw_variance_curve(a_values, V, *, n=5, histories=None):
    """매개변수가 a 하나일 때 x=a, y=log10(분산)인 2차원 곡선을 그린다.

    음수 가중치가 있는 구간은 점선으로 표시한다.
    histories는 draw_variance_surface와 같고, 반환값은 (figure, axes)이다.
    """
    import matplotlib.pyplot as plt

    A = np.array([float(v) for v in a_values])
    V = np.asarray(V, dtype=float)
    mask = np.isfinite(V) & (V > 0)
    if not mask.any():
        raise ValueError("At least one positive finite variance is required")
    L = np.full_like(V, np.nan)
    L[mask] = np.log10(V[mask])

    negative = np.array([mask[p] and has_negative_weight(n, (A[p],)) for p in range(len(A))])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(A, np.where(negative, np.nan, L), color="black", linewidth=1, label="variance")
    if negative.any():
        ax.plot(A, np.where(negative, L, np.nan), color="black", linewidth=1, linestyle="--",
                label="negative weight")
    for index, (name, history) in enumerate((histories or {}).items()):
        node_history, min_history, _ = history
        history_a = [node[-2] for node in node_history]
        ax.scatter(history_a, np.log10(np.array(min_history, dtype=float)),
                   color=f"C{index % 10}", s=4, linewidths=0, label=name)
    ax.set_xlabel("a")
    ax.set_ylabel("log10(variance)")
    ax.set_title(f"Variance of max F_i, n={n}")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", markerscale=2)
    return fig, ax


def draw_variance_surface(a_values, b_values, Z, *, n=6, histories=None, elev=25, azim=-130,
                          max_points=400, max_mesh=60):
    """마우스 드래그와 elev/azim 슬라이더로 돌려볼 수 있는 3차원 그래프.

    음수 가중치가 있는 격자점은 작은 검은 점으로 표시한다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 작은 점으로 표시한다.
    반환값은 (figure, axes, sliders)이다. 슬라이더가 가비지 컬렉션되면
    동작하지 않으므로 호출한 쪽에서 sliders를 변수로 들고 있어야 한다.
    """
    import matplotlib.pyplot as plt

    A, B = np.meshgrid([float(v) for v in a_values], [float(v) for v in b_values], indexing="ij")
    Z = np.asarray(Z, dtype=float)
    mask = np.isfinite(Z) & (Z > 0)
    if not mask.any():
        raise ValueError("At least one positive finite variance is required")
    # 분산이 0.003 ~ 50 정도로 퍼져 있어 log10 축으로 그린다.
    L = np.full_like(Z, np.nan)
    L[mask] = np.log10(Z[mask])

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_axes([0.05, 0.18, 0.9, 0.77], projection="3d", computed_zorder=False)
    # 드래그/슬라이더 회전 시 매 프레임 전체를 다시 그리므로, 촘촘한 격자에서는
    # (0.01 간격 약 4,700점에서 프레임당 0.75초) 화면이 멈춘 것처럼 보인다.
    # 점이 적으면 삼각분할 곡면과 격자점을, 많으면 최대 max_mesh x max_mesh로 줄인
    # 격자 곡면만 그려 회전을 가볍게 한다.
    if mask.sum() <= max_points:
        if mask.sum() >= 3:
            ax.plot_trisurf(A[mask], B[mask], L[mask], cmap="viridis", alpha=0.85, edgecolor="none")
        ax.scatter(A[mask], B[mask], L[mask], color="gray", s=10)
    else:
        # 정의되지 않은 점(a >= b)은 마스킹하고, 불투명 격자 곡면만 그린다.
        ax.plot_surface(A, B, np.ma.masked_invalid(L), cmap="viridis", edgecolor="none",
                        rcount=min(max_mesh, A.shape[0]), ccount=min(max_mesh, A.shape[1]),
                        antialiased=False)
    # 점을 격자 곡면과 같은 간격으로 솎아 회전이 느려지지 않게 한다.
    stride = max(1, int(np.ceil(max(A.shape) / max_mesh)))
    negative = np.zeros_like(mask)
    for p in range(0, A.shape[0], stride):
        for q in range(0, A.shape[1], stride):
            negative[p, q] = mask[p, q] and has_negative_weight(n, (A[p, q], B[p, q]))
    if negative.any():
        ax.scatter(A[negative], B[negative], L[negative], color="black", s=1, linewidths=0,
                   depthshade=False, label="negative weight")
    for index, (name, history) in enumerate((histories or {}).items()):
        node_history, min_history, _ = history
        history_a = [node[-3] for node in node_history]
        history_b = [node[-2] for node in node_history]
        ax.scatter(history_a, history_b, np.log10(np.array(min_history, dtype=float)),
                   color=f"C{index % 10}", s=1, linewidths=0, depthshade=False, label=name)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("log10(variance)")
    ax.set_title(f"Variance of max F_i, n={n}, a < b")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper left", markerscale=4)
    ax.view_init(elev=elev, azim=azim)

    return fig, ax, add_view_sliders(fig, ax, elev, azim)


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
