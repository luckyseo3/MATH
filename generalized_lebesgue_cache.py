"""일반화 르베그 상수의 수치 근삿값과 계산 조건을 JSON으로 보관한다."""
import json
import os
from pathlib import Path
import tempfile

import basis_cache

CACHE_DIR = Path(__file__).resolve().parent / 'generalized_lebegue_cache'


def _path(identity, name):
    relative = basis_cache._cache_path(identity['basis'], name).relative_to(basis_cache.CACHE_DIR)
    return CACHE_DIR / relative


def load(identity, name):
    try:
        data = json.loads(_path(identity, name).read_text(encoding='utf-8'))
        if data['identity'] != identity or data['checksum'] != basis_cache._digest(data['result']):
            return None
        return data['result']
    except (OSError, ValueError, TypeError, KeyError, RecursionError):
        return None


def save(identity, name, result):
    temp_path = None
    try:
        path = _path(identity, name)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {'identity': identity, 'result': result, 'checksum': basis_cache._digest(result)}
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix='.lebesgue-', suffix='.tmp', delete=False) as out:
            temp_path = Path(out.name)
            json.dump(data, out, indent=2, allow_nan=False)
        os.replace(temp_path, path)
    except (OSError, ValueError, TypeError, RecursionError):
        pass
    finally:
        if temp_path is not None:
            try:
                temp_path.unlink(missing_ok=True)
            except OSError:
                pass
