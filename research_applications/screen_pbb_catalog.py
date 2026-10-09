"""Exact LC-to-CSS audit of a pinned external PBB catalog.

No artificial Clifford deformation is applied. The raw polynomial tuple is
expanded in lexicographic Z_ell x Z_m coordinates. Negative outputs contain
XOR certificates of inconsistency for both affine anchor branches on a factor.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")

import collections
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import lc_css as core
from reproducibility import output_directory
import argparse

CATALOG = HERE / 'sources/qcode_discovery/results/campaign7_publication_merged.jsonl'


def matrices(record):
    ell, m = record['ell'], record['m']
    size = ell*m
    def polynomial(terms):
        return [sum(1 << (((i+a) % ell)*m+(j+b) % m) for a, b in terms)
                for i in range(ell) for j in range(m)]
    a, b, c, d = [polynomial(record[name+'_terms']) for name in 'ABCD']
    def transpose(rows):
        return [sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(size)]
    xs = [u | (v << size) for u, v in zip(a, b)] + [0]*size
    zs = [u | (v << size) for u, v in zip(c, d)]
    zs += [u | (v << size) for u, v in zip(transpose(b), transpose(a))]
    # Check valid exponent tuples and full stabilizer commutativity before recognition.
    for name in 'ABCD':
        terms = record[name+'_terms']
        assert len(set(map(tuple, terms))) == len(terms)
        assert all(0 <= i < ell and 0 <= j < m for i, j in terms)
    assert all(((x & v).bit_count() ^ (z & u).bit_count()) % 2 == 0
               for x, z in zip(xs, zs) for u, v in zip(xs, zs))
    return [core.join(x, z, 2*size) for x, z in zip(xs, zs)]


def solve_with_certificate(equations, width, anchor):
    mask, basis = (1 << width)-1, {}
    for index, value in enumerate(itertools.chain([1 << anchor], equations)):
        provenance = 1 << index
        while value & mask:
            pivot = ((value & mask) & -(value & mask)).bit_length()-1
            if pivot not in basis:
                basis[pivot] = value, provenance
                break
            value ^= basis[pivot][0]
            provenance ^= basis[pivot][1]
        else:
            if value >> width & 1:
                indices = list(core.positions(provenance))
                check = 0
                for j in indices:
                    check ^= (1 << anchor) if j == 0 else equations[j-1]
                assert check == 1 << width
                return None, {'anchor_variable': anchor, 'equation_indices': indices,
                              'contradiction_verified': True, 'equations_read': index+1}
    solution = 0
    for pivot in sorted(basis, reverse=True):
        value = basis[pivot][0]
        bit = (value >> width & 1) ^ ((value & solution).bit_count() & 1)
        solution |= bit << pivot
    assert solution >> anchor & 1 == 0
    assert all(((e & mask & solution).bit_count() & 1) == (e >> width & 1) for e in equations)
    return solution, {'anchor_variable': anchor, 'rank': len(basis)}


def classify(rows, n):
    rows = core.validate_rows(rows, n)
    g, _ = core.rref(rows, 2*n)
    components = core.factors(g, n)
    qs = [None]*n
    statistics = []
    for component in components:
        size = len(component)
        local, pivots = core.rref([core.restrict(row, component) for row in g], 2*size)
        eqs = core.preserving_equations(local, pivots, size)
        certs = []
        for anchor in (1, 2):
            solution, details = solve_with_certificate(eqs, 3*size, anchor)
            certs.append(details)
            if solution is not None:
                break
        if solution is None:
            return {'lc_css': False, 'rank': len(g), 'factor_sizes': list(map(len, components)),
                    'rejected_factor': component, 'local_rref_hex': list(map(hex, local)),
                    'negative_certificates': certs}
        for j, i in enumerate(component):
            qs[i] = core.witness(solution >> (3*j) & 1, solution >> (3*j+1) & 1,
                                 solution >> (3*j+2) & 1)
        statistics.append({'n': size, 'rank': len(local), 'branches': len(certs)})
    transformed = core.act(rows, qs)
    xmask = sum(1 << (2*i) for i in range(n))
    xs = list(dict.fromkeys(r & xmask for r in transformed if r & xmask))
    zs = list(dict.fromkeys(r & (xmask << 1) for r in transformed if r & (xmask << 1)))
    assert core.rref(xs+zs, 2*n)[0] == core.rref(transformed, 2*n)[0]
    hx, hz = [core.split(r, n)[0] for r in xs], [core.split(r, n)[1] for r in zs]
    assert all((x & z).bit_count() % 2 == 0 for x in hx for z in hz)
    return {'lc_css': True, 'rank': len(g), 'factor_sizes': list(map(len, components)),
            'q': qs, 'css_x_supports': [list(core.positions(r)) for r in hx],
            'css_z_supports': [list(core.positions(r)) for r in hz],
            'max_check_weight': max((r.bit_count() for r in hx+hz), default=0),
            'max_qubit_participation': max((sum(r >> i & 1 for r in hx+hz) for i in range(n)), default=0),
            'factors': statistics}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=HERE.parent / 'outputs/published_codes')
    args = parser.parse_args()
    out = output_directory(args.output_dir)
    if (out / 'pbb_classification.jsonl').exists():
        parser.error('Classification output already exists; select a fresh --output-dir')
    records = [json.loads(line) for line in CATALOG.read_text().splitlines()]
    for index, record in enumerate(records, start=1):
        record['catalog_line'] = index
        record.setdefault('code_id', f'catalog_line_{index:03d}')
    records.sort(key=lambda r: (r['n'], r['code_id']))
    result_path = out / 'pbb_classification.jsonl'
    completed = {}
    start = time.perf_counter()
    with result_path.open('x', encoding='utf-8') as stream:
        for record in records:
            before = time.perf_counter()
            rows = matrices(record)
            answer = classify(rows, record['n'])
            assert record['n']-answer['rank'] == record['k'], record['code_id']
            output = {k:record[k] for k in ['code_id','catalog_line','n','k','ell','m','d','d_is_exact','lc_any']}
            output['source_lc_any'] = output.pop('lc_any')
            output.update(answer)
            output['input_sha256'] = hashlib.sha256(json.dumps(list(map(hex, rows))).encode()).hexdigest()
            output['seconds_including_construction_and_checks'] = time.perf_counter()-before
            stream.write(json.dumps(output)+'\n'); stream.flush()
            completed[record['code_id']] = output
            if len(completed) % 20 == 0 or output['lc_css'] or output['source_lc_any']:
                print(f"{len(completed)}/{len(records)} {record['code_id']} n={record['n']} "
                      f"LC-CSS={output['lc_css']} source={output['source_lc_any']} "
                      f"{output['seconds_including_construction_and_checks']:.3f}s", flush=True)
    counts = []
    for n in sorted({r['n'] for r in records}):
        rr = [r for r in completed.values() if r['n']==n]
        counts.append({'n':n,'total':len(rr),'lc_css':sum(r['lc_css'] for r in rr),
                       'previous_positive':sum(r['source_lc_any'] for r in rr),
                       'new_positive':sum(r['lc_css'] and not r['source_lc_any'] for r in rr),
                       'negative':sum(not r['lc_css'] for r in rr)})
    summary = {'catalog_size':len(records),'processed':len(completed),'by_size':counts,
               'catalog_sha256':hashlib.sha256(CATALOG.read_bytes()).hexdigest(),
               'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'recognition_core_sha256':hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
               'session_seconds':time.perf_counter()-start,
               'total_seconds':sum(r['seconds_including_construction_and_checks'] for r in completed.values()),
               'new_positive_ids':[r['code_id'] for r in completed.values() if r['lc_css'] and not r['source_lc_any']],
               'positive_mismatches':[r['code_id'] for r in completed.values() if not r['lc_css'] and r['source_lc_any']]}
    (out/'pbb_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
