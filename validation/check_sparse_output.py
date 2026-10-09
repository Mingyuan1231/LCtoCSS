"""Check sparse CSS extraction on previously recorded accepted inputs.

These checks verify sparse extraction, row spans, weights, and incidences.
The general bounds are proved in the manuscript.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import report_output, save_report, utc_now

import hashlib
import json

HERE = Path(__file__).resolve().parent


def extend(pivots, value):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = value
            return True
        value ^= pivots[pivot]
    return False


def basis(rows):
    result = {}
    for row in rows:
        extend(result, row)
    return result


def contains(pivots, value):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return False
        value ^= pivots[pivot]
    return True


def transform(row, maps):
    result = 0
    # Recorded matrices act on columns; the paper uses their transposes.
    for i, (a, b, c, d) in enumerate(maps):
        x, z = (row >> (2*i)) & 1, (row >> (2*i+1)) & 1
        result |= ((a*x) ^ (b*z)) << (2*i)
        result |= ((c*x) ^ (d*z)) << (2*i+1)
    return result


def support(row, n):
    return {i for i in range(n) if (row >> (2*i)) & 3}


def bounds(rows, n):
    sets = [support(row, n) for row in rows]
    weight = max(map(len, sets), default=0)
    degree = max((sum(i in s for s in sets) for i in range(n)), default=0)
    return len(rows), weight, degree


def main():
    report_path = report_output('sparse_output.json')
    source = HERE / 'cases.json'
    records = json.loads(source.read_text(encoding='utf-8'))
    checked = 0
    for case in records['cases']:
        answer = case['projector']
        if not answer['equivalent']:
            continue
        n, rows = case['n'], case['input_rows']
        transformed = [transform(row, answer['local_matrices']) for row in rows]
        assert all(support(a, n) == support(b, n) for a, b in zip(rows, transformed)), case['name']
        xmask = sum(1 << (2*i) for i in range(n))
        zmask = xmask << 1
        pure = [part for row in transformed for part in (row & xmask, row & zmask) if part]
        input_basis, pure_basis = basis(transformed), basis(pure)
        assert len(input_basis) == len(pure_basis) == case['rank'], case['name']
        assert all(contains(input_basis, row) for row in pure), case['name']
        assert all(contains(pure_basis, row) for row in transformed), case['name']
        m, w, delta = bounds(rows, n)
        new_m, new_w, new_delta = bounds(pure, n)
        assert new_m <= 2*m and new_w <= w and new_delta <= 2*delta, case['name']
        selected, pivots = [], {}
        for row in pure:
            if extend(pivots, row):
                selected.append(row)
        assert len(selected) == case['rank'], case['name']
        assert all(row in pure for row in selected), case['name']
        sm, sw, sd = bounds(selected, n)
        assert sm <= new_m and sw <= new_w and sd <= new_delta, case['name']
        checked += 1
    result = {
        'executed_at_utc': utc_now(),
        'description': 'Sparse output extraction on previously accepted cases; recognition not rerun',
        'source_record_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'accepted_inputs_checked': checked,
        'original_support_preserved': True,
        'projected_checks_span_transformed_input': True,
        'check_count_weight_and_incidence_bounds_passed': True,
        'independent_subset_retains_original_projected_rows': True,
        'all_passed': True,
    }
    save_report(result, report_path)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
