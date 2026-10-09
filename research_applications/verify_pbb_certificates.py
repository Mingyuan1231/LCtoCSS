"""Verify the saved PBB classification without importing the recognizer.

Uses NumPy to reconstruct the published block matrices. Negative certificates
are checked against independently constructed ordinary annihilator equations.
The mathematical implication uses the manuscript's two-anchor theorem.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import output_directory

import argparse
import hashlib
import json
import time
import numpy as np

HERE = Path(__file__).resolve().parent


def bits(value):
    return [i for i in range(value.bit_length()) if value >> i & 1]


def canonical(rows):
    basis = {}
    for row in rows:
        while row:
            p = (row & -row).bit_length()-1
            if p not in basis:
                basis[p] = row
                break
            row ^= basis[p]
    for p in sorted(basis, reverse=True):
        for q in basis:
            if q < p and basis[q] >> p & 1:
                basis[q] ^= basis[p]
    return [basis[p] for p in sorted(basis)]


def belongs(row, basis):
    for g in basis:
        p = (g & -g).bit_length()-1
        if row >> p & 1:
            row ^= g
    return row == 0


def reconstruct(record):
    ell, m = record['ell'], record['m']
    size = ell*m
    blocks = []
    for name in 'ABCD':
        block = np.zeros((size, size), dtype=np.uint8)
        for i in range(ell):
            for j in range(m):
                for a, b in record[name+'_terms']:
                    block[i*m+j, ((i+a) % ell)*m+(j+b) % m] ^= 1
        blocks.append(block)
    a, b, c, d = blocks
    zero = np.zeros_like(a)
    x = np.block([[a, b], [zero, zero]])
    z = np.block([[c, d], [b.T, a.T]])
    assert not np.any((x @ z.T + z @ x.T) % 2)
    interleaved = np.empty((2*size, 4*size), dtype=np.uint8)
    interleaved[:, 0::2], interleaved[:, 1::2] = x, z
    return [sum(1 << int(j) for j in np.flatnonzero(row)) for row in interleaved]


def verify_negative(rows, result):
    full_basis = canonical(rows)
    comp = result['rejected_factor']
    n = len(comp)
    assert comp and all(type(i) is int for i in comp)
    assert comp == sorted(set(comp)) and all(0 <= i < result['n'] for i in comp)
    local = [int(h, 16) for h in result['local_rref_hex']]
    projected = [sum(((r >> (2*i)) & 3) << (2*j) for j, i in enumerate(comp)) for r in rows]
    assert canonical(projected) == local
    # The projected subspace embeds in the input, so this is a direct factor.
    for row in local:
        embedded = sum(((row >> (2*j)) & 3) << (2*i) for j, i in enumerate(comp))
        assert belongs(embedded, full_basis)
    # Connected RREF support verifies the factor is indecomposable.
    adjacency = [set() for _ in range(n)]
    for row in local:
        support = set(j//2 for j in bits(row))
        for i in support:
            adjacency[i].update(support)
    seen, todo = set(), [0]
    while todo:
        i = todo.pop()
        if i not in seen:
            seen.add(i)
            todo.extend(adjacency[i]-seen)
    assert len(seen) == n
    pivots = [(g & -g).bit_length()-1 for g in local]
    nonpivots = [j for j in range(2*n) if j not in pivots]
    annihilators = [1 << j for j in nonpivots]
    for beta, g in enumerate(local):
        for t, j in enumerate(nonpivots):
            if g >> j & 1:
                annihilators[t] ^= 1 << pivots[beta]
    assert all((g & h).bit_count() % 2 == 0 for g in local for h in annihilators)
    assert sorted(c['anchor_variable'] for c in result['negative_certificates']) == [1, 2]
    for certificate in result['negative_certificates']:
        total = 0
        indices = certificate['equation_indices']
        assert indices and all(type(index) is int and 0 <= index <= len(local)*len(nonpivots) for index in indices)
        assert len(indices) == len(set(indices))
        for index in indices:
            if index == 0:
                total ^= 1 << certificate['anchor_variable']
                continue
            alpha, t = divmod(index-1, len(nonpivots))
            assert 0 <= alpha < len(local)
            s, h = local[alpha], annihilators[t]
            for i in range(n):
                x, z = (s >> (2*i)) & 1, (s >> (2*i+1)) & 1
                u, v = (h >> (2*i)) & 1, (h >> (2*i+1)) & 1
                total ^= (x*u ^ z*v) << (3*i)
                total ^= (x*v) << (3*i+1)
                total ^= (z*u) << (3*i+2)
                total ^= (z*v) << (3*n)
        assert total == 1 << (3*n), 'Certificate must reduce to 0 = 1'


def verify_positive(rows, result):
    qs = result['q']
    assert len(qs) == result['n']
    assert all(len(q) == 4 and all(type(v) is int and v in (0, 1) for v in q) and (q[0]*q[3] ^ q[1]*q[2]) == 1 for q in qs)
    transformed = []
    for row in rows:
        value = 0
        for i, (a, b, c, d) in enumerate(qs):
            x, z = (row >> (2*i)) & 1, (row >> (2*i+1)) & 1
            value ^= (x*a ^ z*c) << (2*i)
            value ^= (x*b ^ z*d) << (2*i+1)
        transformed.append(value)
    hx, hz = result['css_x_supports'], result['css_z_supports']
    assert all(len(v) == len(set(v)) and all(0 <= i < result['n'] for i in v) for v in hx+hz)
    pure = [sum(1 << (2*i) for i in v) for v in hx]
    pure += [sum(1 << (2*i+1) for i in v) for v in hz]
    assert canonical(pure) == canonical(transformed)
    assert all(len(set(x) & set(z)) % 2 == 0 for x in hx for z in hz)
    assert max(map(len, hx+hz), default=0) == result['max_check_weight']
    assert max((sum(i in support for support in hx+hz) for i in range(result['n'])), default=0) == result['max_qubit_participation']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=HERE / 'results')
    parser.add_argument('--output-dir', type=Path, default=HERE.parent / 'outputs/verification')
    args = parser.parse_args()
    out = output_directory(args.output_dir)
    start = time.perf_counter()
    source = HERE/'sources/qcode_discovery'
    manifest = json.loads((HERE/'source_manifest.json').read_text())
    for entry in manifest['files']:
        assert hashlib.sha256((source/entry['path']).read_bytes()).hexdigest() == entry['sha256']
    records = list(map(json.loads, (source/'results/campaign7_publication_merged.jsonl').read_text().splitlines()))
    results_path = args.results_dir/'pbb_classification.jsonl'
    results = list(map(json.loads, results_path.read_text().splitlines()))
    assert len(records) == len(results) == 368
    assert sorted(r['catalog_line'] for r in results) == list(range(1, 369))
    positive = negative = 0
    for result in results:
        record = records[result['catalog_line']-1]
        for key in ('n', 'k', 'ell', 'm'):
            assert result[key] == record[key]
        assert result['code_id'] == record.get('code_id', f"catalog_line_{result['catalog_line']:03d}")
        assert type(result['lc_css']) is bool
        assert result['source_lc_any'] == record['lc_any']
        rows = reconstruct(record)
        assert hashlib.sha256(json.dumps(list(map(hex, rows))).encode()).hexdigest() == result['input_sha256']
        assert len(canonical(rows)) == result['rank'] == record['n']-record['k']
        if result['lc_css']:
            verify_positive(rows, result)
            positive += 1
        else:
            verify_negative(rows, result)
            negative += 1
    report = {'verified_entries': len(results), 'positive_witnesses': positive,
              'negative_entries': negative, 'negative_branch_certificates': 2*negative,
              'all_passed': True, 'seconds': time.perf_counter()-start,
              'imports_recognizer': False, 'equation_method': 'ordinary annihilator bilinear formula',
              'results_sha256': hashlib.sha256(results_path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'pbb_independent_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
