"""Gauss-Newton (s_newton_method.gauss_newton) 결과를 variance_cache/newton/n{n}.json에 저장한다.

python newton_runner.py convert            : 스크래치패드의 n = 9..21 결과를 저장소 형식으로 저장
python newton_runner.py run START SECONDS  : n = START, START+1, ...을 병렬로 계산하고 SECONDS초 뒤 중단
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
from fractions import Fraction
ENDPOINT_VARIABLE = Fraction(1)  # 끝 노드 a: 노드의 양 끝은 ±a
VARIANCE_CACHE = pathlib.Path(__file__).resolve().parent / 'variance_cache'  # 이 a로 구한 결과
import json, os, sys, time
import multiprocessing as mp_
from fractions import Fraction
from pathlib import Path

OUT = VARIANCE_CACHE / 'newton'
SCRATCH = Path(__file__).with_name('results')
TOL = Fraction(1, 10**60)
DIGITS, DPS = 40, 100
WORKERS = 6



def write_json(path, data):
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=1), encoding='utf-8')
    os.replace(tmp, path)


def describe(node):
    """최종 노드에서 Fraction으로 최댓값, 가중치, 분산을 정확히 다시 계산한다."""
    import calculate_extremum
    import calculate_variance
    M, w = calculate_extremum.exact_extrema_and_weights(node, kind="max")
    inner = M[1:-1]
    # 정확한 분산은 분모가 수천 자리라 저장하지 않는다. 필요하면 free에서 다시 계산할 수 있다.
    return dict(variance=float(calculate_variance.exact_variance(inner)),
                common_max=float(sum(inner) / len(inner)),
                max_range=float(max(inner) - min(inner)),
                min_weight=float(min(w)),
                maxima=[float(v) for v in inner])


def settings(n, start):
    return dict(n=n, method='s_newton_method.gauss_newton', start=start,
                variance_tolerance='1e-60', digits=DIGITS, dps=DPS)


def history_json(history):
    import s_newton_method as g
    nodes, variances = history
    return [dict(free=[str(v) for v in g.node_to_variable(node)], variance=float(var))
            for node, var in zip(nodes, variances)]


def convert():
    import s_newton_method as g
    for n in range(9, 22):
        r = json.loads((SCRATCH / f'n{n}.json').read_text())
        free = [Fraction(x) for x in r['free']]
        node = g.variable_to_node(free, n, ENDPOINT_VARIABLE)
        data = dict(settings=settings(n, 'ESC (float)'), status='done', iterations=None,
                    elapsed=r['time'], free=r['free'], free_decimal=[f'{float(v):.20f}' for v in free],
                    **describe(node), history=None,
                    note='2026-11-21 첫 실행 결과. history는 기록하지 않았다.')
        write_json(OUT / f'n{n}.json', data)
        print(f'saved n={n}', flush=True)


def run_one(n):
    import nodes, s_newton_method as g
    t0 = time.time()
    partial = OUT / f'n{n}.partial.json'

    # 이어서 계산: 이전 실행의 partial이 있으면 그 마지막 노드에서 시작한다.
    previous, previous_elapsed, resumes = [], 0.0, 0
    if partial.exists():
        old = json.loads(partial.read_text(encoding='utf-8'))
        previous = old['history'][:-1]          # 마지막 항목은 새 history의 0번째(시작점)와 같다.
        previous_elapsed = old['elapsed']
        resumes = old.get('resumes', 0) + 1
        init = g.variable_to_node([Fraction(x) for x in old['history'][-1]['free']], n, ENDPOINT_VARIABLE)
    else:
        init = [Fraction(v) for v in nodes.ESC_nodes(n, dps=60, endpoint_variable=ENDPOINT_VARIABLE)]

    def merged(history):
        return previous + history_json(history)

    def elapsed():
        return previous_elapsed + time.time() - t0

    def save_partial(history):
        h = merged(history)
        write_json(partial, dict(settings=settings(n, 'ESC (dps=60)'), status='running',
                                 iterations=len(h) - 1, elapsed=elapsed(), resumes=resumes,
                                 history=h))

    try:
        history = g.gauss_newton(init, TOL, digits=DIGITS, dps=DPS, callback=save_partial)
        node = history[0][-1]
        free = g.node_to_variable(node)
        h = merged(history)
        data = dict(settings=settings(n, 'ESC (dps=60)'), status='done',
                    iterations=len(h) - 1, elapsed=elapsed(), resumes=resumes,
                    free=[str(v) for v in free], free_decimal=[f'{float(v):.20f}' for v in free],
                    **describe(node), history=h)
        write_json(OUT / f'n{n}.json', data)
        if partial.exists():
            partial.unlink()
        msg = f"n={n} done  var {data['variance']:.2e}  iters {data['iterations']}  {data['elapsed']:.0f}s"
    except Exception as e:
        write_json(OUT / f'n{n}.failed.json', dict(settings=settings(n, 'ESC (dps=60)'), status='failed',
                                                   error=f'{type(e).__name__}: {e}', elapsed=time.time() - t0))
        msg = f'n={n} FAILED {type(e).__name__}: {e}'
    print(time.strftime('%H:%M:%S'), msg, flush=True)
    return n


def run(start, seconds):
    deadline = time.time() + seconds
    pool = mp_.Pool(WORKERS)
    pending = {}
    n_next = start
    while time.time() < deadline:
        for n in [n for n, r in pending.items() if r.ready()]:
            del pending[n]
        while len(pending) < WORKERS:
            if (OUT / f'n{n_next}.json').exists():   # 이미 끝난 n은 건너뛴다.
                n_next += 1
                continue
            pending[n_next] = pool.apply_async(run_one, (n_next,))
            print(time.strftime('%H:%M:%S'), f'start n={n_next}', flush=True)
            n_next += 1
        time.sleep(5)
    pool.terminate()
    print(time.strftime('%H:%M:%S'), f'deadline reached; stopped n={sorted(pending)} (partial files kept)', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'convert':
        convert()
    else:
        run(int(sys.argv[2]), float(sys.argv[3]))
