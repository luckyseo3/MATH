"""생성된 SymPy 기저를 실행 가능한 코드 없이 JSON으로 보관한다."""
import hashlib
import inspect
import json
import os
from pathlib import Path
import re
import tempfile

import sympy as sp

CACHE_DIR = Path(__file__).resolve().parent / 'basis_cache'
_FORMAT_VERSION = 1


def _encode(expr):
    expr = sp.sympify(expr)
    if isinstance(expr, sp.Symbol):
        return {'type': 'Symbol', 'name': expr.name, 'assumptions': expr.assumptions0}
    if isinstance(expr, sp.Integer):
        return {'type': 'Integer', 'value': str(expr)}
    if isinstance(expr, sp.Rational):
        return {'type': 'Rational', 'p': str(expr.p), 'q': str(expr.q)}
    if isinstance(expr, sp.Float):
        return {'type': 'Float', 'mpf': list(expr._mpf_), 'precision': expr._prec}
    if expr.func in (sp.Add, sp.Mul, sp.Pow):
        return {'type': expr.func.__name__, 'args': [_encode(arg) for arg in expr.args]}
    raise ValueError(f'Unsupported cache expression: {expr!s}')


def _decode(data):
    kind = data['type']
    if kind == 'Symbol':
        return sp.Symbol(data['name'], **data['assumptions'])
    if kind == 'Integer':
        return sp.Integer(data['value'])
    if kind == 'Rational':
        return sp.Rational(data['p'], data['q'])
    if kind == 'Float':
        return sp.Float._new(tuple(data['mpf']), data['precision'])
    constructors = {'Add': sp.Add, 'Mul': sp.Mul, 'Pow': sp.Pow}
    if kind in constructors:
        return constructors[kind](*(_decode(arg) for arg in data['args']))
    raise ValueError('Unknown cached expression type.')


def _digest(value):
    payload = json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()


def cache_identity(nodes, functions):
    """입력 노드와 기저 계산 함수의 소스가 같을 때만 동일한 키를 만든다."""
    return {'format_version': _FORMAT_VERSION, 'sympy_version': sp.__version__,
            'algorithm': _digest([inspect.getsource(function) for function in functions]),
            'nodes': [_encode(node) for node in nodes]}


def _cache_path(identity, name):
    if name is None:
        name = 'custom-' + _digest(identity['nodes'])[:16]
    if not isinstance(name, str) or not name:
        raise ValueError('Basis name must be a non-empty string.')
    if (not re.fullmatch(r'[A-Za-z0-9_-]+', name)
            or name.upper() in {'CON', 'PRN', 'AUX', 'NUL', *(f'COM{i}' for i in range(1, 10)),
                                *(f'LPT{i}' for i in range(1, 10))}):
        name = 'name-' + _digest(name)[:16]
    return CACHE_DIR / name / f"n{len(identity['nodes'])}.json"


def node_source_name(source):
    name = getattr(source, '__name__', None)
    aliases = {'GL_nodes': 'GL', 'GLo_nodes': 'GLo', 'Chebyshev_nodes': 'Chev',
               'ESC_nodes': 'ESC', 'ESCLo_nodes': 'ESCLo', 'uniform_nodes': 'uniform'}
    return aliases.get(name, name)


def load(identity, name=None):
    path = _cache_path(identity, name)
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
        encoded = data['basis']
        if (data['identity'] != identity or data['checksum'] != _digest(encoded)
                or not isinstance(encoded, list) or len(encoded) != len(identity['nodes']) + 1):
            return None
        return [_decode(expr) for expr in encoded]
    except (OSError, ValueError, TypeError, KeyError, IndexError, AttributeError,
            RecursionError, ArithmeticError):
        return None


def save(identity, basis, name=None):
    """캐시 기록 실패가 원래의 기저 계산을 막지 않도록 한다."""
    temp_path = None
    try:
        encoded = [_encode(expr) for expr in basis]
        data = {'identity': identity, 'basis': encoded, 'checksum': _digest(encoded)}
        path = _cache_path(identity, name)
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix='.basis-', suffix='.tmp', delete=False) as out:
            temp_path = Path(out.name)
            json.dump(data, out, ensure_ascii=True, indent=2)
        os.replace(temp_path, path)
    except (OSError, ValueError, TypeError, RecursionError, OverflowError):
        pass
    finally:
        if temp_path is not None:
            try:
                temp_path.unlink(missing_ok=True)
            except OSError:
                pass
