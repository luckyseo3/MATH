"""0 대칭·끝점 포함 노드에서 매개변수에 따른 내부 F_i 최댓값 분산을 3차원으로 그린다."""
from fractions import Fraction
import gzip
import hashlib
import inspect
import json
import multiprocessing
import os
from pathlib import Path
import time

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


def volume_slab_path(n, N, i):
    return (Path(__file__).resolve().parent / 'variance_cache' / f'n{n}' / f'N{N}'
            / f'a{i:04d}.npz')


def compute_volume_slab(task):
    """a=i/N로 고정한 모든 (b, c)의 결과를 Fraction으로 정확히 계산해 저장한다.

    a{i}.npz: V(내부 분산), V_all(끝점 포함 분산), M(F_0~F_n의 최댓값), weights, positive.
    실수 값은 정확한 유리수를 가장 가까운 배정밀도로 반올림한 것이다.
    a{i}.exact.txt.gz: 내부 분산의 정확한 값 "분자/분모"를 한 줄에 하나씩 저장한다.
    (b, c)는 저장하지 않으며 np.triu_indices의 순서로 복원한다. 정의되지 않은 점은 NaN이다.
    """
    n, N, i = task
    path = volume_slab_path(n, N, i)
    count = (N - 1 - i) * (N - 2 - i) // 2
    V = np.full(count, np.nan)
    V_all = np.full(count, np.nan)
    M_values = np.full((count, n + 1), np.nan)
    weight_values = np.full((count, n), np.nan)
    positive = np.zeros(count, dtype=bool)

    path.parent.mkdir(parents=True, exist_ok=True)
    exact_path = path.with_name(path.stem + '.exact.txt.gz')
    exact_temp = path.with_name(path.stem + '.exact.tmp')
    a = Fraction(i, N)
    position = 0
    with gzip.open(exact_temp, 'wt', encoding='ascii', newline='\n') as exact_file:
        for j in range(i + 1, N - 1):
            for k in range(j + 1, N):
                try:
                    M, weights = maximum.exact_maximums_and_weights(
                        symmetric_nodes(n, (a, Fraction(j, N), Fraction(k, N))))
                except ValueError:
                    exact_file.write('nan\n')
                else:
                    inner = maximum.exact_variance(M[1:-1])
                    V[position] = float(inner)
                    V_all[position] = float(maximum.exact_variance(M))
                    M_values[position] = [float(value) for value in M]
                    weight_values[position] = [float(value) for value in weights]
                    positive[position] = min(weights) > 0
                    exact_file.write(f'{inner.numerator}/{inner.denominator}\n')
                position += 1

    # npz가 완료 표시이므로 마지막에 만든다. 중단되어도 깨진 파일이 완료로 남지 않는다.
    os.replace(exact_temp, exact_path)
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, V=V, V_all=V_all, M=M_values, weights=weight_values, positive=positive)
    os.replace(temp_path, path)
    return count


def variance_volume(N, *, n=8, workers=None):
    """0 < a < b < c < 1인 간격 1/N 격자의 분산을 병렬로 계산한다.

    n=8이면 노드 [-1, -c, -b, -a, a, b, c, 1], n=9이면 가운데에 0이 추가된다.
    a 값 하나가 끝날 때마다 variance_cache/n{n}/N{N}/에 저장하므로,
    중단 후 다시 실행하면 끝난 a는 건너뛰고 이어서 계산한다.
    """
    if n not in (8, 9):
        raise ValueError("Volume needs exactly three parameters (n = 8 or 9)")
    tasks = [(n, N, i) for i in range(1, N - 2) if not volume_slab_path(n, N, i).exists()]
    if not tasks:
        return

    total = sum((N - 1 - i) * (N - 2 - i) // 2 for _, _, i in tasks)
    print(f"n={n}, 1/N={1 / N:g}: {len(tasks)} slabs, {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_volume_slab, tasks):
            done += count
            elapsed = time.perf_counter() - start
            remaining = elapsed * (total - done) / done if done else 0
            print(f"{done:,}/{total:,} points, elapsed {elapsed / 3600:.2f} h, "
                  f"remaining about {remaining / 3600:.2f} h", flush=True)


def load_variance_volume(N, *, n=8):
    """양수 가중치 격자점들의 (a, b, c 인덱스, 분산)과 격자 위의 국소 최솟값을 반환한다.

    국소 최솟값은 주변 26개 격자점 중 양수 가중치인 점들보다 분산이 작은 점이다.
    interior가 False이면 주변에 양수 가중치가 아니거나 a<b<c를 벗어나는 점이 있다는 뜻이다.
    """
    def layer(i):
        grid = np.full((N + 1, N + 1), np.nan)
        if 1 <= i <= N - 3:
            with np.load(volume_slab_path(n, N, i)) as data:
                V, positive = data['V'], data['positive']
            j, k = np.triu_indices(N - 1 - i, 1)
            keep = positive & (V > 0)
            grid[j[keep] + i + 1, k[keep] + i + 1] = V[keep]
        return grid

    offsets = [(dj, dk) for dj in (-1, 0, 1) for dk in (-1, 0, 1)]
    points = {'a': [], 'b': [], 'c': [], 'V': []}
    minima = {'a': [], 'b': [], 'c': [], 'V': [], 'interior': []}
    previous, current, following = layer(0), layer(1), layer(2)
    for i in range(1, N - 2):
        finite = np.isfinite(current)
        j, k = np.nonzero(finite)
        points['a'].append(np.full(len(j), i, dtype=np.uint16))
        points['b'].append(j.astype(np.uint16))
        points['c'].append(k.astype(np.uint16))
        points['V'].append(current[finite])

        # 인덱스 0과 N은 항상 NaN이므로 roll로 생기는 경계 값은 비교에 영향이 없다.
        is_minimum = finite.copy()
        interior = finite.copy()
        for depth, other in enumerate((previous, current, following)):
            for dj, dk in offsets:
                if depth == 1 and dj == 0 and dk == 0:
                    continue
                neighbor = np.roll(other, (dj, dk), axis=(0, 1))
                is_minimum &= ~(neighbor <= current)
                interior &= np.isfinite(neighbor)
        j, k = np.nonzero(is_minimum)
        minima['a'].append(np.full(len(j), i, dtype=np.uint16))
        minima['b'].append(j.astype(np.uint16))
        minima['c'].append(k.astype(np.uint16))
        minima['V'].append(current[is_minimum])
        minima['interior'].append(interior[is_minimum])

        previous, current, following = current, following, layer(i + 2)

    return ({name: np.concatenate(parts) for name, parts in points.items()},
            {name: np.concatenate(parts) for name, parts in minima.items()})


def compute_plane_row(task):
    """n=6, 7에서 a=i/N로 고정한 모든 b의 결과를 Fraction으로 정확히 계산해 저장한다.

    저장 형식은 compute_volume_slab과 같고, j번째 값은 b=(i+1+j)/N이다.
    """
    n, N, i = task
    path = volume_slab_path(n, N, i)
    count = N - 1 - i
    V = np.full(count, np.nan)
    V_all = np.full(count, np.nan)
    M_values = np.full((count, n + 1), np.nan)
    weight_values = np.full((count, n), np.nan)
    positive = np.zeros(count, dtype=bool)

    path.parent.mkdir(parents=True, exist_ok=True)
    exact_path = path.with_name(path.stem + '.exact.txt.gz')
    exact_temp = path.with_name(path.stem + '.exact.tmp')
    a = Fraction(i, N)
    with gzip.open(exact_temp, 'wt', encoding='ascii', newline='\n') as exact_file:
        for position, j in enumerate(range(i + 1, N)):
            try:
                M, weights = maximum.exact_maximums_and_weights(
                    symmetric_nodes(n, (a, Fraction(j, N))))
            except ValueError:
                exact_file.write('nan\n')
                continue
            inner = maximum.exact_variance(M[1:-1])
            V[position] = float(inner)
            V_all[position] = float(maximum.exact_variance(M))
            M_values[position] = [float(value) for value in M]
            weight_values[position] = [float(value) for value in weights]
            positive[position] = min(weights) > 0
            exact_file.write(f'{inner.numerator}/{inner.denominator}\n')

    os.replace(exact_temp, exact_path)
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, V=V, V_all=V_all, M=M_values, weights=weight_values, positive=positive)
    os.replace(temp_path, path)
    return count


def variance_plane(N, *, n=6, workers=None):
    """0 < a < b < 1인 간격 1/N 격자의 분산을 Fraction으로 정확히 계산한다.

    n=6이면 노드 [-1, -b, -a, a, b, 1], n=7이면 가운데에 0이 추가된다.
    a 값 하나가 끝날 때마다 variance_cache/n{n}/N{N}/에 저장하므로 중단해도 이어서 계산한다.
    """
    if n not in (6, 7):
        raise ValueError("Plane needs exactly two parameters (n = 6 or 7)")
    tasks = [(n, N, i) for i in range(1, N - 1) if not volume_slab_path(n, N, i).exists()]
    if not tasks:
        return

    total = sum(N - 1 - i for _, _, i in tasks)
    print(f"n={n}, 1/N={1 / N:g}: {len(tasks)} rows, {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_plane_row, tasks, chunksize=4):
            done += count
            elapsed = time.perf_counter() - start
            remaining = elapsed * (total - done) / done if done else 0
            print(f"\r{done:,}/{total:,} points, elapsed {elapsed / 60:.1f} min, "
                  f"remaining about {remaining / 60:.1f} min", end="", flush=True)
    print()


def active_row_path(n, N, i):
    return volume_slab_path(n, N, i).with_suffix('.active.npz')


def compute_active_row(task):
    """a=i/N로 고정한 모든 b에서 내부 F_1~F_{n-1}의 최댓값을 주는 후보 이름을 정확히 구해 저장한다.

    label은 F_i별 후보 이름을 ","로 이은 문자열이며, 가중치가 0인 점은 빈 문자열이다.
    """
    n, N, i = task
    path = active_row_path(n, N, i)
    labels = []
    a = Fraction(i, N)
    for j in range(i + 1, N):
        try:
            _, _, active = maximum.exact_maximums_and_weights(
                symmetric_nodes(n, (a, Fraction(j, N))), active=True)
        except ValueError:
            labels.append('')
        else:
            labels.append(",".join(active[1:-1]))
    temp_path = path.with_name(path.stem + '.tmp')
    with open(temp_path, 'wb') as out:
        np.savez(out, labels=np.array(labels))
    os.replace(temp_path, path)
    return len(labels)


def active_plane(N, *, n=6, workers=None):
    """variance_plane과 같은 격자에서 최댓값 후보 라벨을 계산한다. 중단해도 이어서 계산한다."""
    tasks = [(n, N, i) for i in range(1, N - 1) if not active_row_path(n, N, i).exists()]
    if not tasks:
        return
    total = sum(N - 1 - i for _, _, i in tasks)
    print(f"active candidates n={n}, 1/N={1 / N:g}: {total:,} points to compute")
    start = time.perf_counter()
    done = 0
    with multiprocessing.Pool(workers or max(1, os.cpu_count() - 1)) as pool:
        for count in pool.imap_unordered(compute_active_row, tasks, chunksize=4):
            done += count
            elapsed = time.perf_counter() - start
            print(f"\r{done:,}/{total:,} points, elapsed {elapsed / 60:.1f} min, remaining about "
                  f"{elapsed * (total - done) / done / 60:.1f} min", end="", flush=True)
    print()


def load_active_plane(N, *, n=6):
    """(N+1)x(N+1) 정수 라벨 격자(-1은 정의되지 않음)와 라벨 이름 목록을 반환한다."""
    ids = np.full((N + 1, N + 1), -1, dtype=np.int32)
    names = {}
    for i in range(1, N - 1):
        with np.load(active_row_path(n, N, i)) as data:
            row = data['labels']
        ids[i, i + 1:N] = [names.setdefault(label, len(names)) if label else -1 for label in row]
    return ids, list(names)


def load_variance_plane(N, *, n=6):
    """(N+1)x(N+1) 격자 V[i, j](a=i/N, b=j/N)와 positive를 반환한다. 정의되지 않으면 NaN."""
    V = np.full((N + 1, N + 1), np.nan)
    positive = np.zeros((N + 1, N + 1), dtype=bool)
    for i in range(1, N - 1):
        with np.load(volume_slab_path(n, N, i)) as data:
            V[i, i + 1:N] = data['V']
            positive[i, i + 1:N] = data['positive']
    return V, positive


def variance_colormap(low, high, color_max=0.0, red_at=4.0):
    """log10(분산)용 고정 색: color_max까지 viridis, 그 위는 주황, red_at에서 빨강, 최댓값은 검붉은색."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    split = max(color_max, low + 1.0)
    top = max(high, split + 1.0)
    red = min((red_at - split) / (top - split), 0.9)
    viridis = plt.get_cmap("viridis")
    cmap = LinearSegmentedColormap.from_list("variance", [
        *((0.5 * t, viridis(t)) for t in np.linspace(0, 1, 32)),
        (0.5 + 0.25 * red, "darkorange"),
        (0.5 + 0.5 * red, "red"),
        (1.0, "#3d0000"),
    ], N=1024)
    return cmap, TwoSlopeNorm(vmin=low, vcenter=split, vmax=top)


def cached_compass_search(n, h, h_limit, initial_node):
    """maximum.compass_search의 결과를 variance_cache/n{n}/compass.json을 거쳐 구한다.

    h, h_limit, 시작 노드와 탐색 코드가 모두 같을 때만 저장된 결과를 재사용한다.
    """
    source = "".join(inspect.getsource(function) for function in (
        maximum.compass_search, maximum.numeric_variance_of_maximums_of_F,
        maximum.numeric_maximums_of_F))
    code = hashlib.sha256(source.encode('utf-8')).hexdigest()[:16]
    key = "|".join([code, repr(h), repr(h_limit), *(repr(float(v)) for v in initial_node)])

    path = Path(__file__).resolve().parent / 'variance_cache' / f'n{n}' / 'compass.json'
    runs = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    if key in runs:
        return runs[key]

    node_history, min_history, h_history = maximum.compass_search(n, h, h_limit, initial_node)
    runs[key] = [[[float(v) for v in node] for node in node_history],
                 [float(v) for v in min_history], [float(v) for v in h_history]]
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + '.tmp')
    temp_path.write_text(json.dumps(runs), encoding='utf-8')
    os.replace(temp_path, path)
    return runs[key]


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


def draw_variance_volume(N, *, n=8, max_points=3000, points_limit=30000, max_size=10,
                         color_max=0.0, red_at=4.0, zoom_limit=50, elev=25, azim=-130,
                         histories=None):
    """(a, b, c) 위치에 양수 가중치 격자점을 찍고, 분산은 점의 크기와 색으로 나타낸다.

    log10(분산)이 작을수록 점이 크다. 슬라이더 값 이하인 점만 표시하며,
    회전이 느려지지 않도록 표시할 점은 max_points개 이내로 고르게 솎는다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 가는 선으로 표시한다.
    반환값은 (figure, axes, sliders)이다. 슬라이더가 가비지 컬렉션되면
    동작하지 않으므로 호출한 쪽에서 sliders를 변수로 들고 있어야 한다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.widgets import Slider

    points, minima = load_variance_volume(N, n=n)
    if len(points['V']) == 0:
        raise ValueError("No grid point has all positive weights")
    L = np.log10(points['V'])
    low, high = float(L.min()), float(L.max())
    span = (high - low) or 1.0
    # 색은 슬라이더와 무관하게 log10(분산) 값으로만 정해진다.
    # 빨강 이상은 발산에 가까운 영역으로 본다.
    cmap, norm = variance_colormap(low, high, color_max, red_at)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.02, 0.39, 0.82, 0.57], projection="3d", computed_zorder=False)
    ax.set_proj_type("ortho")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("c")
    ax.view_init(elev=elev, azim=azim)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.88, 0.44, 0.02, 0.48]),
                 label="log10(variance)")

    for interior, marker, color, size, label in ((True, "*", "magenta", 25, "local min"),
                                                 (False, "x", "gray", 8, "local min on boundary")):
        chosen = minima['interior'] == interior
        if chosen.any():
            ax.scatter(minima['a'][chosen] / N, minima['b'][chosen] / N, minima['c'][chosen] / N,
                       marker=marker, color=color, s=size, depthshade=False, zorder=3, axlim_clip=True,
                       label=f"{label} ({int(chosen.sum())})")
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        ax.plot(path[:, -4], path[:, -3], path[:, -2], color=plt.get_cmap("tab20")(index % 20), linewidth=0.7,
                marker="o", markersize=1, zorder=4, axlim_clip=True, label=name)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper left", fontsize=7, ncols=2)

    A, B, C = points['a'] / N, points['b'] / N, points['c'] / N
    # 줌 중심은 처음에 분산이 가장 작은 격자점이며, 슬라이더로 옮길 수 있다.
    best = int(np.argmin(L))
    scatter = [None]

    def update(_=None):
        if scatter[0] is not None:
            scatter[0].remove()
        width = 1 / zoom_slider.val
        center = np.array([slider.val for slider in center_sliders])
        lower = np.clip(center - width / 2, 0, 1 - width)
        upper = lower + width
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
        ax.set_zlim(lower[2], upper[2])

        # 줌을 해도 점이 드물어지지 않도록 화면 안의 점들 중에서 고른다.
        threshold = threshold_slider.val
        selected = np.nonzero((L <= threshold)
                              & (lower[0] <= A) & (A <= upper[0])
                              & (lower[1] <= B) & (B <= upper[1])
                              & (lower[2] <= C) & (C <= upper[2]))[0]
        total = len(selected)
        shown = int(points_slider.val)
        if total > shown:
            selected = selected[np.linspace(0, total - 1, shown, dtype=int)]
        scatter[0] = ax.scatter(A[selected], B[selected], C[selected], c=L[selected], cmap=cmap,
                                norm=norm, s=0.5 + (max_size - 0.5) * (high - L[selected]) / span,
                                linewidths=0, depthshade=False, zorder=2, axlim_clip=True)
        ax.set_title(f"n={n}, grid 1/{N}: {len(selected):,} of {total:,} points, "
                     f"log10(var) <= {threshold:.2f}, zoom x{zoom_slider.val:.1f}")
        fig.canvas.draw_idle()

    center_sliders = [Slider(fig.add_axes([0.15, y, 0.7, 0.025]), f"{name} center", 0, 1,
                             valinit=float(value[best]))
                      for name, y, value in (("a", 0.35, A), ("b", 0.305, B), ("c", 0.26, C))]
    zoom_slider = Slider(fig.add_axes([0.15, 0.215, 0.7, 0.025]), "zoom", 1, zoom_limit, valinit=1)
    points_slider = Slider(fig.add_axes([0.15, 0.17, 0.7, 0.025]), "points", 100,
                           max(points_limit, max_points), valinit=max_points, valstep=100)
    threshold_slider = Slider(fig.add_axes([0.15, 0.125, 0.7, 0.025]), "log10(var) <=", low, high,
                              valinit=float(np.percentile(L, 5)))
    for slider in (*center_sliders, zoom_slider, points_slider, threshold_slider):
        slider.on_changed(update)

    def on_scroll(event):
        if event.inaxes is ax:
            zoom_slider.set_val(float(np.clip(zoom_slider.val * 1.25 ** event.step, 1, zoom_limit)))

    fig.canvas.mpl_connect("scroll_event", on_scroll)
    update()

    return fig, ax, (*center_sliders, zoom_slider, points_slider, threshold_slider,
                     *add_view_sliders(fig, ax, elev, azim))


def label_text(label):
    return "(" + ", ".join(part.replace("tau", "τ") for part in label.split(",")) + ")"


def draw_active_regions(N, *, n=6, histories=None):
    """(a, b) 평면을 내부 F_1~F_{n-1}의 최댓값 후보 조합별로 색칠한 지도.

    색이 바뀌는 곳이 분산이 미분 불가능할 수 있는 경계이다. 검은 실선/점선/일점쇄선은
    가중치 ω_1, ω_2, ω_3 = 0인 곡선(분산이 발산하는 곳)이고, 음수 가중치 영역은 흐리게 칠한다.
    후보가 정확히 같은 값으로 겹친 점(라벨에 "/")은 빨간 원, 분산 최소 격자점은 자홍색 별이다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    ids, names = load_active_plane(N, n=n)
    V, positive = load_variance_plane(N, n=n)
    weights = np.full((N + 1, N + 1, n), np.nan)
    for i in range(1, N - 1):
        with np.load(volume_slab_path(n, N, i)) as data:
            weights[i, i + 1:N] = data['weights']

    order = sorted(range(len(names)), key=lambda k: -(ids == k).sum())
    rank = np.full(len(names), -1)
    rank[order] = np.arange(len(order))
    # 영역이 10개를 넘으면 tab10 색이 겹치므로 tab20을 쓴다.
    palette = "tab10" if len(order) <= 10 else "tab20"
    base = plt.get_cmap(palette)(np.arange(len(order)) % plt.get_cmap(palette).N)
    cmap = ListedColormap(base)
    cmap.set_bad("white")
    shown = np.ma.masked_less(np.where(ids >= 0, rank[ids], -1), 0).astype(float)

    fig, ax = plt.subplots(figsize=(11, 8))
    fig.subplots_adjust(left=0.07, right=0.62, top=0.93, bottom=0.08)
    extent = (-0.5 / N, 1 + 0.5 / N, -0.5 / N, 1 + 0.5 / N)
    ax.imshow(shown.T, origin="lower", extent=extent, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5,
              interpolation="nearest")
    # 음수 가중치 영역은 흰색을 덮어 흐리게 보이게 한다.
    fade = np.zeros((N + 1, N + 1, 4))
    fade[(ids >= 0) & ~positive] = (1, 1, 1, 0.55)
    ax.imshow(fade.transpose(1, 0, 2), origin="lower", extent=extent, interpolation="nearest")

    grid = np.arange(N + 1) / N
    for k, style in zip(range((n + 1) // 2), ("-", "--", "-.", ":")):
        W = np.ma.masked_invalid(weights[:, :, k])
        if np.nanmin(W) < 0 < np.nanmax(W):
            ax.contour(grid, grid, W.T, levels=[0], colors="black", linewidths=1, linestyles=style)
    ties = np.array([("/" in name) for name in names])
    tie_i, tie_j = np.nonzero((ids >= 0) & ties[np.maximum(ids, 0)])
    if len(tie_i):
        ax.scatter(tie_i / N, tie_j / N, s=60, facecolors="none", edgecolors="red", linewidths=1.2, zorder=5)
    P = np.where(positive & np.isfinite(V) & (V > 0), V, np.inf)
    best = np.unravel_index(np.argmin(P), P.shape)
    ax.scatter([best[0] / N], [best[1] / N], marker="*", s=120, color="magenta", edgecolors="black",
               linewidths=0.5, zorder=6)
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        ax.plot(path[:, -3], path[:, -2], color="black", linewidth=0.8, zorder=4)

    handles = [Patch(color=base[r], label=f"{label_text(names[k])}  {int((ids == k).sum()):,}")
               for r, k in enumerate(order)]
    handles += [plt.Line2D([], [], color="black", linestyle=style, label=f"ω_{k + 1} = 0" if 2 * k + 1 == n else f"ω_{k + 1} = ω_{n - k} = 0")
                for k, style in zip(range((n + 1) // 2), ("-", "--", "-.", ":"))]
    handles += [plt.Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor="red", linestyle="",
                           label=f"exact tie ({len(tie_i)})"),
                plt.Line2D([], [], marker="*", color="magenta", linestyle="", markersize=10, label="min variance")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8,
              title=f"argmax of (F_1, ..., F_{n - 1})", title_fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_title(f"n={n}, grid 1/{N}: active candidates (faded = negative weight)")
    return fig, ax


def draw_variance_plane(N, *, n=6, histories=None, max_mesh=100, zoom_limit=50,
                        color_max=0.0, red_at=4.0, elev=25, azim=-130):
    """간격 1/N 격자의 log10(분산) 곡면을 그리고, zoom 슬라이더로 최솟값 근처를 확대한다.

    회전이 느려지지 않도록 화면 안의 격자를 최대 max_mesh x max_mesh로 솎아 그리며,
    확대할수록 솎는 간격이 줄어 원래 1/N 격자가 그대로 보인다.
    색은 log10(분산) 값으로만 정해진다. 음수 가중치 점은 작은 검은 점,
    양수 가중치 점들 사이의 격자 국소 최솟값은 자홍색 별로 표시한다.
    histories는 {이름: compass_search의 반환값}이며, 경로마다 다른 색의 가는 선으로 표시한다.
    반환값은 (figure, axes, sliders)이다.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.widgets import Slider

    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    V, positive = load_variance_plane(N, n=n)
    ids = load_active_plane(N, n=n)[0] if active_row_path(n, N, 1).exists() else None
    defined = np.isfinite(V) & (V > 0)
    L = np.full_like(V, np.nan)
    L[defined] = np.log10(V[defined])
    good = defined & positive
    if not good.any():
        raise ValueError("No grid point has all positive weights")
    low, high = float(np.nanmin(L)), float(np.nanmax(L))
    cmap, norm = variance_colormap(low, high, color_max, red_at)

    # 양수 가중치 점들 사이의 격자 국소 최솟값 (주변 8점과 비교)
    P = np.where(good, L, np.nan)
    is_minimum = good.copy()
    interior = good.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                neighbor = np.roll(P, (di, dj), axis=(0, 1))
                is_minimum &= ~(neighbor <= P)
                interior &= np.isfinite(neighbor)
    is_minimum &= interior
    best = np.unravel_index(np.nanargmin(P), P.shape)

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_axes([0.02, 0.28, 0.82, 0.63], projection="3d", computed_zorder=False)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("log10(variance)")
    ax.view_init(elev=elev, azim=azim)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.88, 0.33, 0.02, 0.52]),
                 label="log10(variance)")

    i, j = np.nonzero(is_minimum)
    ax.scatter(i / N, j / N, L[i, j], marker="*", color="magenta", s=40, depthshade=False,
               zorder=5, axlim_clip=True, label=f"grid local min ({len(i)})")
    for index, (name, history) in enumerate((histories or {}).items()):
        path = np.array(history[0], dtype=float)
        ax.plot(path[:, -3], path[:, -2], np.log10(np.array(history[1], dtype=float)),
                color=plt.get_cmap("tab10")(index % 10), linewidth=0.8, marker="o", markersize=1,
                zorder=4, axlim_clip=True, label=name)
    ax.legend(loc="upper left", fontsize=7)

    drawn = []

    def update(_=None):
        while drawn:
            drawn.pop().remove()
        width = 1 / zoom_slider.val
        center = np.array([slider.val for slider in center_sliders])
        lower = np.clip(center - width / 2, 0, 1 - width)
        first = np.floor(lower * N).astype(int)
        last = np.minimum(np.ceil((lower + width) * N).astype(int), N)
        stride = max(1, int(np.ceil((last - first).max() / max_mesh)))
        rows = np.arange(first[0], last[0] + 1, stride)
        cols = np.arange(first[1], last[1] + 1, stride)
        window = L[np.ix_(rows, cols)]
        A, B = np.meshgrid(rows / N, cols / N, indexing="ij")
        drawn.append(ax.plot_surface(A, B, np.ma.masked_invalid(window), cmap=cmap, norm=norm,
                                     rstride=1, cstride=1, edgecolor="none", antialiased=False,
                                     zorder=1))
        negative = np.isfinite(window) & ~positive[np.ix_(rows, cols)]
        dots = max(1, int(np.ceil(max(window.shape) / 60)))
        thin = np.zeros_like(negative)
        thin[::dots, ::dots] = True
        negative &= thin
        if negative.any():
            drawn.append(ax.scatter(A[negative], B[negative], window[negative], color="black", s=1,
                                    linewidths=0, depthshade=False, zorder=2, axlim_clip=True))
        if ids is not None:
            # 이웃한 두 격자점의 최댓값 후보가 다르면 그 사이에 짧은 선분을 곡면 높이로 그린다.
            label = ids[np.ix_(rows, cols)]
            half = stride / (2 * N)
            segments = []
            for axis in (0, 1):
                first_part = (slice(None, -1), slice(None)) if axis == 0 else (slice(None), slice(None, -1))
                second_part = (slice(1, None), slice(None)) if axis == 0 else (slice(None), slice(1, None))
                change = ((label[first_part] != label[second_part]) & (label[first_part] >= 0)
                          & (label[second_part] >= 0)
                          & np.isfinite(window[first_part]) & np.isfinite(window[second_part]))
                for p, q in zip(*np.nonzero(change)):
                    a0, b0 = A[p, q], B[p, q]
                    z = (window[p, q] + window[p + (axis == 0), q + (axis == 1)]) / 2
                    if axis == 0:
                        segments.append([(a0 + half, b0 - half, z), (a0 + half, b0 + half, z)])
                    else:
                        segments.append([(a0 - half, b0 + half, z), (a0 + half, b0 + half, z)])
            if segments:
                lines = Line3DCollection(segments, colors="white", linewidths=1.2, zorder=3, axlim_clip=True)
                ax.add_collection3d(lines)
                drawn.append(lines)
        ax.set_xlim(rows[0] / N, rows[-1] / N)
        ax.set_ylim(cols[0] / N, cols[-1] / N)
        if np.isfinite(window).any():
            ax.set_zlim(np.nanmin(window), np.nanmax(window))
        ax.set_title(f"n={n}, grid 1/{N}, zoom x{zoom_slider.val:.1f}: "
                     f"{len(rows)}x{len(cols)} mesh, step {stride}/{N}\n"
                     f"min log10(var) = {L[best]:.3f} at a={best[0] / N:g}, b={best[1] / N:g}")
        fig.canvas.draw_idle()

    # 줌 중심은 처음에 분산이 가장 작은 격자점이며, 슬라이더로 옮길 수 있다.
    center_sliders = [Slider(fig.add_axes([0.15, y, 0.7, 0.03]), f"{name} center", 0, 1,
                             valinit=best[axis] / N)
                      for axis, (name, y) in enumerate((("a", 0.23), ("b", 0.18)))]
    zoom_slider = Slider(fig.add_axes([0.15, 0.13, 0.7, 0.03]), "zoom", 1, zoom_limit, valinit=1)
    for slider in (*center_sliders, zoom_slider):
        slider.on_changed(update)

    def on_scroll(event):
        if event.inaxes is ax:
            zoom_slider.set_val(float(np.clip(zoom_slider.val * 1.25 ** event.step, 1, zoom_limit)))

    fig.canvas.mpl_connect("scroll_event", on_scroll)
    update()

    return fig, ax, (*center_sliders, zoom_slider, *add_view_sliders(fig, ax, elev, azim))


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
