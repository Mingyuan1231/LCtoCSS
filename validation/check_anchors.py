"""Direct finite checks for factorwise anchors and the two-branch criterion.

This enumerates local maps and row-space membership directly.  It does not
implement or call the manuscript's affine-system solver.  Run with Python 3.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import report_output, save_report


from itertools import product
import json
import time


IDENTITY = (1, 0, 0, 1)
TRACE_ONE = tuple((a, b, c, 1 ^ a) for a, b, c in product(range(2), repeat=3))
CLIFFORD = tuple(q for q in product(range(2), repeat=4)
                 if q[0] * q[3] ^ q[1] * q[2] == 1)


def label(word):
    return sum({"I": 0, "X": 1, "Z": 2, "Y": 3}[p] << (2 * i)
               for i, p in enumerate(word))


def word(v, n):
    return "".join("IXZY"[(v >> (2 * i)) & 3] for i in range(n))


def span(rows):
    result = {0}
    for r in rows:
        result |= {v ^ r for v in tuple(result)}
    return result


def apply(v, maps):
    out = 0
    for i, (a, b, c, d) in enumerate(maps):
        x, z = (v >> (2 * i)) & 1, (v >> (2 * i + 1)) & 1
        out |= ((a * x ^ c * z) | ((b * x ^ d * z) << 1)) << (2 * i)
    return out


def isotropic(rows, n):
    return all(sum(((u >> (2 * i)) & 1) * ((v >> (2 * i + 1)) & 1)
                   + ((u >> (2 * i + 1)) & 1) * ((v >> (2 * i)) & 1)
                   for i in range(n)) % 2 == 0
               for u in rows for v in rows)


def scalar_masks(rows, n):
    space = span(rows)
    result = []
    for d in product(range(2), repeat=n):
        mask = sum((3 * d[i]) << (2 * i) for i in range(n))
        if all((r & mask) in space for r in rows):
            result.append(d)
    return result


def preserving_trace_one_maps(rows, n):
    space = span(rows)
    return [local for local in product(TRACE_ONE, repeat=n)
            if all(apply(r, local) in space for r in rows)]


def css_count_under_all_local_cliffords(rows, n):
    x_mask = sum(1 << (2 * i) for i in range(n))
    count = 0
    for local in product(CLIFFORD, repeat=n):
        transformed = [apply(r, local) for r in rows]
        transformed_space = span(transformed)
        if all((r & x_mask) in transformed_space for r in transformed):
            count += 1
    return count


def main():
    report_path = report_output('anchors.json')
    started = time.perf_counter()
    bad_words = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    bad_rows = [label(w) for w in bad_words]
    assert isotropic(bad_rows, 5) and len(span(bad_rows)) == 16
    bad_maps = preserving_trace_one_maps(bad_rows, 5)
    assert len(bad_maps) == 2
    assert all(all(p[1] == p[2] == 1 for p in local) for local in bad_maps)
    bad_masks = scalar_masks(bad_rows, 5)
    assert bad_masks == [(0,) * 5, (1,) * 5]
    lc_count = css_count_under_all_local_cliffords(bad_rows, 5)
    assert lc_count == 0
    result = {
        "method": "Direct finite row spaces and exhaustive local-map enumeration; no affine solver",
        "nonmaximal_anchor_required": {
            "words": bad_words, "n": 5, "rank": 4,
            "trace_one_assignments_checked": 8 ** 5,
            "preserving_maps": bad_maps,
            "preserving_scalar_masks": bad_masks,
            "local_clifford_assignments_checked": 6 ** 5,
            "css_witness_count": lc_count,
        },
        "nonmaximal_both_branches_required": [],
    }

    n = 6
    base = [(1 << (2 * i)) ^ sum(2 << (2 * j) for j in range(n)
                                      if abs(i - j) == 1)
            for i in range(n - 1)]
    anchor = 2  # zero-based qubit 3
    for name, twist, expected in [
        ("base", IDENTITY, (8, 8)),
        ("upper_shear", (1, 1, 0, 1), (0, 8)),
        ("lower_shear", (1, 0, 1, 1), (8, 0)),
    ]:
        maps = [IDENTITY] * n
        maps[anchor] = twist
        rows = [apply(r, maps) for r in base]
        assert isotropic(rows, n) and len(span(rows)) == 32
        preservers = preserving_trace_one_maps(rows, n)
        assert len(preservers) == 8
        masks = scalar_masks(rows, n)
        assert masks == [(0,) * n, (1,) * n]
        counts = (sum(p[anchor][1] == 0 for p in preservers),
                  sum(p[anchor][2] == 0 for p in preservers))
        assert counts == expected
        # Undo the shear, then Hadamards on odd-numbered path vertices.
        undo = [IDENTITY] * n
        undo[anchor] = twist  # each shear is an involution
        restored = [apply(r, undo) for r in rows]
        hadamards = [(0, 1, 1, 0) if i % 2 == 0 else IDENTITY
                     for i in range(n)]
        css = [apply(r, hadamards) for r in restored]
        x_mask = sum(1 << (2 * i) for i in range(n))
        assert all((r & x_mask) in span(css) for r in css)
        result["nonmaximal_both_branches_required"].append({
            "case": name, "words": [word(r, n) for r in rows],
            "n": n, "rank": 5, "anchor_qubit_one_based": anchor + 1,
            "trace_one_assignments_checked": 8 ** n,
            "trace_one_preservers": len(preservers),
            "local_anchor_possibilities": sorted({p[anchor] for p in preservers}),
            "b_branch_solutions": counts[0], "c_branch_solutions": counts[1],
            "preserving_scalar_masks": masks,
            "explicit_CSS_witness_verified": True,
        })

    # A single anchor on a decomposable code does not control the other factor.
    decomposable_rows = [1] + [r << 2 for r in bad_rows]
    decomposable_map = ((1, 0, 0, 0),) + bad_maps[0]
    decomposable_space = span(decomposable_rows)
    assert all(apply(r, decomposable_map) in decomposable_space
               for r in decomposable_rows)
    assert decomposable_map[0][1] == 0
    assert scalar_masks(decomposable_rows, 6) == [
        (0, 0, 0, 0, 0, 0), (0, 1, 1, 1, 1, 1),
        (1, 0, 0, 0, 0, 0), (1, 1, 1, 1, 1, 1),
    ]
    result["factorwise_anchors_required"] = {
        "words": [word(r, 6) for r in decomposable_rows],
        "n": 6, "rank": 5,
        "preserving_trace_one_map_with_b1_zero": decomposable_map,
        "reason_still_not_LC_CSS": "Five-qubit direct factor fails all 7776 local Clifford assignments",
    }
    result["elapsed_seconds"] = time.perf_counter() - started
    result["all_assertions_passed"] = True
    save_report(result, report_path)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
