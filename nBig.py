"""n = 10 ~ 20에서 ESC 노드로 시작하는 compass search를 수행하고 결과를 분석한다.

- 초깃값: ESC 노드를 정확히 대칭으로 맞춘 것 (끝점 ±1, 홀수 n이면 가운데 0)
- 분산 계산: maximum.precise_variance_of_maximums_of_F (mpmath DPS자리)
  float(numeric_variance_of_maximums_of_F)는 n=16에서 최댓값 오차가 약 1e-7, n=20에서 약 1e-4이고,
  ESC의 분산(약 5e-5)과 오차가 비슷해져 개선 여부를 판정할 수 없으므로 쓰지 않는다.
- 결과는 variance_cache/nBig/n{n}.json에 저장되어, 설정이 같으면 다시 계산하지 않는다.
  n마다 별도 프로세스에서 병렬로 계산한다.
"""
import json
import multiprocessing
import os
from pathlib import Path
import time

import numpy as np
import matplotlib.pyplot as plt

import generate_nodes as nodes
import draw_basis as draw
import maximum

N_VALUES = range(10, 21)
H_INITIAL = 0.001   # ESC와 최적 노드의 거리는 0.01 정도이다.
H_LIMIT = 1e-8      # 위치를 1e-8까지 다듬는다. 분산은 DPS자리로 계산하므로 비교가 오차에 묻히지 않는다.
DPS = 40            # n=20에서 Fraction 정확값과 15자리 이상 일치한다.

CACHE_DIR = Path(__file__).resolve().parent / 'variance_cache' / 'nBig'


def esc_initial_node(n):
    """ESC 노드를 정확히 대칭으로 맞춘다. compass search는 대칭 쌍을 함께 움직이므로 대칭이 유지된다."""
    half = [float(v) for v in nodes.ESC_nodes(n)[(n + 1) // 2:-1]]
    return [-1.0, *[-v for v in reversed(half)], *([0.0] if n % 2 else []), *half, 1.0]


def variance(node):
    return maximum.precise_variance_of_maximums_of_F(node, DPS)


def settings(n):
    return {'n': n, 'h': H_INITIAL, 'h_limit': H_LIMIT, 'dps': DPS, 'initial': esc_initial_node(n)}


def run(n):
    """n 하나의 compass search를 수행하고 저장한다. 같은 설정의 결과가 있으면 그것을 읽는다."""
    path = CACHE_DIR / f'n{n}.json'
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['settings'] == settings(n):
            return saved

    start = time.perf_counter()
    node_history, min_history, h_history = maximum.compass_search(
        n, H_INITIAL, H_LIMIT, esc_initial_node(n), variance_function=variance)
    elapsed = time.perf_counter() - start

    final = node_history[-1]
    M, weights = maximum.precise_maximums_and_weights(final, DPS)
    result = {
        'settings': settings(n),
        'esc_variance': float(min_history[0]),
        'variance': float(min_history[-1]),
        'node': [float(v) for v in final],
        'steps': len(node_history) - 1,
        'elapsed': elapsed,
        'min_weight': float(min(weights)),
        'maxima': [float(v) for v in M[1:-1]],
        'history': {'variance': [float(v) for v in min_history], 'h': [float(v) for v in h_history]},
    }
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(result, indent=1), encoding='utf-8')
    os.replace(temp, path)
    print(f"n={n} done: {result['steps']} steps, {elapsed / 60:.1f} min", flush=True)
    return result


def print_summary(results):
    print()
    print(f"{'n':>3} {'ESC variance':>14} {'compass variance':>17} {'ratio':>10} {'steps':>6} "
          f"{'min weight':>11} {'max M - min M':>14} {'time':>8}")
    for n, r in results.items():
        spread = max(r['maxima']) - min(r['maxima'])
        print(f"{n:>3} {r['esc_variance']:>14.6e} {r['variance']:>17.6e} "
              f"{r['variance'] / r['esc_variance']:>10.2e} {r['steps']:>6} "
              f"{r['min_weight']:>11.4g} {spread:>14.3e} {r['elapsed'] / 60:>6.1f} m")
    print()
    print("compass search 결과 노드 (양수 쪽, 0과 1 제외)")
    for n, r in results.items():
        half = r['node'][(n + 1) // 2:-1]
        print(f"n={n:>2}: " + ", ".join(f"{v:.10f}" for v in half))


if __name__ == "__main__":

    ##################################################################################################################
    # 1️⃣ n = 10 ~ 20에서 ESC 노드로 시작하는 compass search (n마다 병렬, 결과는 캐시)
    with multiprocessing.Pool(min(len(N_VALUES), max(1, os.cpu_count() - 1))) as pool:
        # 오래 걸리는 큰 n부터 시작해 전체 시간을 줄인다.
        computed = dict(zip(reversed(N_VALUES), pool.map(run, reversed(N_VALUES))))
    results = {n: computed[n] for n in N_VALUES}

    ##################################################################################################################
    # 2️⃣ 차수별 분산 출력
    # 분산은 내부 F_1 ~ F_{n-1} 최댓값의 모분산이다. 가중치가 음수가 되면 min weight가 음수로 표시된다.
    print_summary(results)

    ##################################################################################################################
    # 3️⃣ 차수별 기저 함수: 슬라이더(또는 좌우 방향키)로 n을, 라디오 버튼으로 ESC / compass를 고른다.
    basis_sources = {
        "ESC": nodes.ESC_nodes,
        "compass": {n: np.array(r['node']) for n, r in results.items()},
    }
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=N_VALUES[0],
        n_min=N_VALUES[0],
        n_max=N_VALUES[-1],
        initial_basis="compass",
    )

    plt.show()
