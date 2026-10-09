"""Verify all saved public-code certificates without importing either solver.

Run with normal Python (not -O). The input matrices, their identities, direct
factors, and equation sums are reconstructed independently from the source data.
The original classification outputs and original timing records are read-only.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import output_directory

from datetime import datetime, timezone
import argparse
import hashlib
import json
import platform
import time
import numpy as np
import verify_pbb_certificates as verifier

HERE = Path(__file__).resolve().parent
XYZ_SUBSETS = {
    'X': [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
    'Y': [[1, 5, 9], [2, 6, 7], [3, 4, 8], [3, 5, 7], [2, 4, 9], [1, 6, 8]],
    'Z': [[1, 4, 7], [2, 5, 8], [3, 6, 9]],
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_digest(rows):
    # This exact serialization also defines input_sha256 in the input convention.
    return hashlib.sha256(json.dumps(list(map(hex, rows))).encode('utf-8')).hexdigest()


def xyz_rows():
    blocks = {}
    for pauli, subsets in XYZ_SUBSETS.items():
        values = []
        for subset in subsets:
            matrix = np.ones((1, 1), dtype=np.uint8)
            for i in range(1, 10):
                matrix = np.kron(matrix, np.ones((1, 2), dtype=np.uint8)
                                 if i in subset else np.eye(2, dtype=np.uint8))
            values.append(matrix)
        blocks[pauli] = np.vstack(values)
    assert all(not np.any((blocks[a] @ blocks[b].T) % 2)
               for a, b in [('X', 'Y'), ('Y', 'Z'), ('Z', 'X')])
    rows = []
    for pauli in ('X', 'Y', 'Z'):
        for row in blocks[pauli]:
            pure = sum(1 << (2*int(i)) for i in np.flatnonzero(row))
            rows.append(pure if pauli == 'X' else pure << 1 if pauli == 'Z' else pure ^ (pure << 1))
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: certificate checks use assertions.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=HERE / 'results')
    parser.add_argument('--output-dir', type=Path, default=HERE.parent / 'outputs/verification')
    args = parser.parse_args()
    out = output_directory(args.output_dir)
    start = time.perf_counter()
    source_manifest = json.loads((HERE/'source_manifest.json').read_text())
    source = HERE/'sources/qcode_discovery'
    for entry in source_manifest['files']:
        assert digest(source/entry['path']) == entry['sha256']
    xyz_manifest = json.loads((HERE/'xyz_source_manifest.json').read_text())
    assert digest(HERE/xyz_manifest['file']) == xyz_manifest['sha256']
    xyz_input = json.loads((HERE/xyz_manifest['file']).read_text(encoding='utf-8'))
    assert xyz_input['subsets'] == XYZ_SUBSETS and xyz_input['n'] == 512
    catalog_path = source/'results/campaign7_publication_merged.jsonl'
    result_path = args.results_dir/'pbb_classification.jsonl'
    records = list(map(json.loads, catalog_path.read_text().splitlines()))
    results = list(map(json.loads, result_path.read_text().splitlines()))
    assert len(records) == len(results) == 368
    assert sorted(r['catalog_line'] for r in results) == list(range(1, 369))
    mapping, positive, negative, sizes = [], 0, 0, []
    for result in results:
        line = result['catalog_line']
        source_record = records[line-1]
        assert result['code_id'] == source_record.get('code_id', f'catalog_line_{line:03d}')
        for key in ('n', 'k', 'ell', 'm', 'd', 'd_is_exact'):
            assert result[key] == source_record[key]
        assert result['source_lc_any'] == source_record['lc_any']
        assert type(result['lc_css']) is bool
        rows = verifier.reconstruct(source_record)
        sha = rows_digest(rows)
        assert sha == result['input_sha256']
        assert len(verifier.canonical(rows)) == result['rank'] == result['n']-result['k']
        if result['lc_css']:
            verifier.verify_positive(rows, result)
            positive += 1
        else:
            verifier.verify_negative(rows, result)
            negative += 1
            sizes.extend(len(c['equation_indices']) for c in result['negative_certificates'])
        mapping.append({'dataset': 'PBB', 'catalog_line': line, 'code_id': result['code_id'],
                        'n': result['n'], 'k': result['k'], 'input_sha256': sha,
                        'lc_css': result['lc_css']})
    assert (positive, negative) == (11, 357)
    assert all(r['lc_css'] == r['source_lc_any'] for r in results)

    xyz_path = args.results_dir/'xyz_is_512.json'
    xyz = json.loads(xyz_path.read_text())
    assert xyz['source'] == 'https://arxiv.org/html/2607.14988v1'
    assert xyz['equation'] == 102 and xyz['subsets'] == XYZ_SUBSETS
    rows = xyz_rows()
    assert list(map(hex, rows)) == xyz['input_rows_hex']
    assert len(rows) == 768 and xyz['n'] == 512 and xyz['k'] == 9
    assert len(verifier.canonical(rows)) == xyz['rank'] == 503
    assert xyz['lc_css'] is False
    verifier.verify_negative(rows, xyz)
    mapping.append({'dataset': 'XYZ-IS', 'source': xyz['source'], 'equation': 102,
                    'n': 512, 'k': 9, 'input_sha256': rows_digest(rows), 'lc_css': False})
    assert all(name not in sys.modules for name in ('screen_pbb_catalog', 'run_demo', 'lc_css'))
    report = {'all_passed': True, 'verified_inputs': 369, 'positive_witnesses': positive,
              'negative_inputs': negative+1, 'negative_branch_certificates': 2*(negative+1),
              'pbb_certificate_equation_count_range': [min(sizes), max(sizes)],
              'verifies_saved_xyz_without_solver': True, 'imports_recognizer': False,
              'elapsed_seconds': time.perf_counter()-start,
              'timing_scope': 'Complete independent verification, including source hashes and matrix reconstruction; excluding report writes.',
              'verified_at_utc': datetime.now(timezone.utc).isoformat(),
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                          'platform': platform.platform()},
              'catalog_sha256': digest(catalog_path), 'pbb_results_sha256': digest(result_path),
              'xyz_results_sha256': digest(xyz_path),
              'entrypoint_sha256': digest(Path(__file__)),
              'verifier_module_sha256': digest(Path(verifier.__file__))}
    (out/'input_certificate_map.json').write_text(json.dumps(mapping, indent=2)+'\n')
    (out/'all_certificates_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
