"""Reconstruct and classify the published XYZ-IS example in Eq. (102).

Source: Baldelli et al., arXiv:2607.14988v1, Eqs. (94)--(102).
The original 512-qubit example is used without an added LC deformation.
An independent Kronecker-product reconstruction verifies the input, then the
ordinary-annihilator verifier checks the saved negative certificates.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.reproducibility import output_directory

import argparse
import hashlib
import json
import time
import numpy as np
import screen_pbb_catalog as scan
import verify_pbb_certificates as verify

HERE = Path(__file__).resolve().parent
SUBSETS = {
    'X': [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
    'Y': [[1, 5, 9], [2, 6, 7], [3, 4, 8], [3, 5, 7], [2, 4, 9], [1, 6, 8]],
    'Z': [[1, 4, 7], [2, 5, 8], [3, 6, 9]],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=HERE.parent / 'outputs/published_codes')
    args = parser.parse_args()
    out = output_directory(args.output_dir)
    if (out / 'xyz_is_512.json').exists():
        parser.error('XYZ output already exists; select a fresh --output-dir')
    rows, checks = [], {}
    for kind, subsets in SUBSETS.items():
        blocks = []
        for subset in subsets:
            mask = sum(1 << (9-i) for i in subset)
            free = [i for i in range(9) if mask >> i & 1]
            for fixed in range(512):
                if fixed & mask:
                    continue
                support = [fixed | sum(((v >> j) & 1) << i for j, i in enumerate(free))
                           for v in range(1 << len(free))]
                row = sum(1 << (2*i) for i in support)
                rows.append(row if kind == 'X' else row << 1 if kind == 'Z' else row ^ (row << 1))
            matrix = np.ones((1, 1), dtype=np.uint8)
            for i in range(1, 10):
                block = np.ones((1, 2), dtype=np.uint8) if i in subset else np.eye(2, dtype=np.uint8)
                matrix = np.kron(matrix, block)
            blocks.append(matrix)
        checks[kind] = np.vstack(blocks)
    assert all(not np.any((checks[a] @ checks[b].T) % 2)
               for a, b in [('X', 'Y'), ('Y', 'Z'), ('Z', 'X')])
    independent_rows = []
    for kind in SUBSETS:
        for check in checks[kind]:
            value = sum(1 << (2*int(i)) for i in np.flatnonzero(check))
            independent_rows.append(value if kind == 'X' else value << 1 if kind == 'Z' else value ^ (value << 1))
    assert independent_rows == rows
    start = time.perf_counter()
    answer = scan.classify(rows, 512)
    seconds = time.perf_counter()-start
    answer.update({'n': 512, 'k': 9, 'source': 'https://arxiv.org/html/2607.14988v1',
                   'equation': 102, 'subsets': SUBSETS, 'recognition_seconds': seconds,
                   'input_rows_hex': list(map(hex, rows))})
    assert len(verify.canonical(rows)) == answer['rank'] == 503
    assert not answer['lc_css']
    verify.verify_negative(rows, answer)
    answer['independent_reconstruction_verified'] = True
    answer['independent_negative_certificates_verified'] = True
    answer['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out/'xyz_is_512.json').write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps({k: answer[k] for k in ['n', 'k', 'rank', 'lc_css', 'factor_sizes',
                       'recognition_seconds', 'independent_reconstruction_verified',
                       'independent_negative_certificates_verified']}, indent=2))


if __name__ == '__main__':
    main()
