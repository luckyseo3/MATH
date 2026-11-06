import sympy as sp
import numpy as np
import inspect


def _constant_identity(nodes, mode):
    import basis_cache
    import generate_basis as generator

    funcs = (generator.lagrange_poly, generator.lagrange_basis, generator.weights,
             generator.Normalized_Lagrange, generator.Cumulative_Lagrange,
             generator.New_Basis, generator.get_New_Basis)
    return {'format_version': 1, 'basis': basis_cache.cache_identity(nodes, funcs),
            'numpy_version': np.__version__, 'mode': mode,
            'x_range': [-1, 1], 'x_samples': 4001,
            'algorithm': basis_cache._digest([inspect.getsource(generalized_Lebegue_function),
                                               inspect.getsource(_sampled_constant),
                                               inspect.getsource(generalized_Lebegue_constant_withA)])}


def generalized_Lebegue_function(basis):
    """Return the generalized Lebesgue function as a SymPy expression in x."""
    return sp.Add(*(sp.Abs(f) for f in basis))  


def generalized_constant_for_nodes(nodes, basis_name=None):
    """노드에 대응하는 일반화 르베그 상수의 표본 근삿값을 반환한다."""
    import generate_basis as generator
    import generalized_lebesgue_cache as constant_cache

    nodes = np.asarray(nodes)
    if np.iscomplexobj(nodes) and np.any(np.abs(nodes.imag) > 1e-10):
        raise ValueError('Nodes must be real.')
    nodes = np.asarray(nodes.real, dtype=float)
    if nodes.ndim != 1 or not len(nodes) or not np.all(np.isfinite(nodes)) or len(np.unique(nodes)) != len(nodes):
        raise ValueError('Expected distinct finite real nodes.')
    identity = _constant_identity(nodes, 'fixed')
    result = constant_cache.load(identity, basis_name)
    if isinstance(result, (float, int)) and not isinstance(result, bool) and np.isfinite(result):
        return float(result)
    weights = np.asarray(generator.weights(generator.lagrange_basis(nodes))[0], dtype=float)
    if not np.all(np.isfinite(weights)) or np.any(np.abs(weights) < 1e-12):
        raise ValueError('A normalization weight is zero or numerically near zero.')
    functions = generator.get_New_Basis(nodes, basis_name=basis_name)
    result = _sampled_constant(functions, np.linspace(-1, 1, 4001))
    constant_cache.save(identity, basis_name, result)
    return result

def generalized_Lebegue_constant_withBasisAndDegree(basis, start, end):
    """노드 생성 함수 basis(n)에 대해 각 n의 상수 근삿값을 점 그래프로 그린다.

    n은 노드 개수이며 start와 end를 모두 포함한다. 최댓값은 [-1, 1]의
    4001개 표본점에서 계산한다. (figure, axes)를 반환하며 show는 호출하지 않는다.
    """
    import matplotlib.pyplot as plt
    import basis_cache

    if not callable(basis):
        raise TypeError('basis must be a node generator: basis(n) -> n nodes.')
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer))
           for v in (start, end)) or not 1 <= start <= end:
        raise ValueError('Require positive integers with start <= end.')
    constants = []
    for n in range(start, end + 1):
        try:
            nodes = np.asarray(basis(n), dtype=float)
            if (nodes.ndim != 1 or len(nodes) != n or not np.all(np.isfinite(nodes))
                    or len(np.unique(nodes)) != n):
                raise ValueError('Expected n distinct finite real nodes.')
            name = basis_cache.node_source_name(basis)
            constants.append(generalized_constant_for_nodes(nodes, basis_name=name))
        except Exception as error:
            raise ValueError(f'Cannot calculate the constant for n={n}: {error}') from error

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(list(range(start, end + 1)), constants)
    ax.set(xlabel='n (number of nodes)', ylabel='Generalized Lebesgue constant (approx.)',
           title=f'Generalized Lebesgue constant: {getattr(basis, "__name__", "basis")}')
    ax.grid(True)
    return fig, ax

def generalized_Lebegue_constant_withA(n):
    """n=4 또는 5인 대칭 노드에서 a에 따른 상수 근삿값을 그린다.

    a=0.01,...,0.99, x 구간의 4001개 점을 사용한다. 가중치가 0인
    위치에서는 선을 끊는다. (figure, axes)를 반환하며 show는 호출하지 않는다.
    """
    import matplotlib.pyplot as plt
    import generate_basis as generator
    import generalized_lebesgue_cache as constant_cache

    if isinstance(n, (bool, np.bool_)) or n not in (4, 5) or not isinstance(n, (int, np.integer)):
        raise ValueError('n must be 4 or 5.')
    a = generator.a
    nodes = [-1, -a, a, 1] if n == 4 else [-1, -a, 0, a, 1]
    poles = [1 / np.sqrt(3)] if n == 4 else [1 / np.sqrt(5), np.sqrt(3 / 5)]
    a_values = np.sort(np.concatenate((np.arange(1, 100) / 100, poles)))
    x_values = np.linspace(-1, 1, 4001)
    identity = _constant_identity(nodes, 'adjustable_a')
    identity['a_values'] = a_values.tolist()
    result = constant_cache.load(identity, 'adjustable')
    if (not isinstance(result, list) or len(result) != len(a_values)
            or any(v is not None and (not isinstance(v, (int, float)) or not np.isfinite(v)) for v in result)):
        functions = generator.get_New_Basis(nodes, basis_name='adjustable')
        evaluators = [sp.lambdify((generator.x, a), f, 'numpy') for f in functions]
        result = []
        for value in a_values:
            if any(np.isclose(value, pole, rtol=0, atol=1e-12) for pole in poles):
                result.append(None)
                continue
            values = np.zeros_like(x_values)
            for evaluate in evaluators:
                values += np.abs(np.broadcast_to(evaluate(x_values, value), x_values.shape))
            result.append(float(np.max(values)))
        constant_cache.save(identity, 'adjustable', result)
    constants = [np.nan if value is None else value for value in result]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(a_values, constants, marker='.', markersize=3)
    for pole in poles:
        ax.axvline(pole, color='gray', linestyle='--', alpha=0.7)
    ax.set(xlabel='a', ylabel='Generalized Lebesgue constant (approx.)',
           title=f'Generalized Lebesgue constant: n={n}', xlim=(0, 1))
    ax.grid(True)
    return fig, ax


def _sampled_constant(basis, x_values):
    """기저 함수 절댓값 합의 표본 최댓값을 반환한다."""
    import generate_basis as generator

    values = np.zeros_like(x_values)
    for f in basis:
        evaluate = sp.lambdify(generator.x, f, 'numpy')
        values += np.abs(np.broadcast_to(evaluate(x_values), x_values.shape))
    if not np.all(np.isfinite(values)):
        raise ValueError('Basis values must be finite.')
    return float(np.max(values))
