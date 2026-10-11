"""Check the distributed file set and SHA-256 digests using the standard library."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRECTORIES = {'.git', '__pycache__', '.venv', 'venv', '.pytest_cache', 'outputs', 'build', 'dist'}


def distributed_files():
    return {
        p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*')
        if p.is_file()
        and not any(part in IGNORED_DIRECTORIES for part in p.relative_to(ROOT).parts)
        and p.suffix not in ('.pyc', '.pyo', '.exe')
        and p.name not in ('SHA256SUMS', 'graph_exhaustive', '.DS_Store', 'Thumbs.db')
    }


def main():
    expected = {}
    for line in (ROOT / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest, relative = line.split('  ', 1)
        path = (ROOT / relative).resolve()
        if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Invalid digest in SHA256SUMS')
        if not path.is_relative_to(ROOT) or relative in expected:
            raise ValueError('Invalid or duplicate path in SHA256SUMS')
        expected[relative] = digest
    actual = distributed_files()
    if set(expected) != set(actual):
        raise ValueError(f'File set differs: missing={sorted(set(expected)-set(actual))}, extra={sorted(set(actual)-set(expected))}')
    for relative, path in sorted(actual.items()):
        if path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest() != expected[relative]:
            raise ValueError(f'Integrity check failed: {relative}')
    print(json.dumps({'all_passed': True, 'verified_files': len(expected)}, indent=2))


if __name__ == '__main__':
    main()
