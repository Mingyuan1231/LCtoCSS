"""Independent row-convention recognition against exact local-Clifford oracles.

The mathematical checks use only standard-library code and independent solvers.
The manuscript algorithm uses affine equations. The exact oracle enumerates
actual rowspace elements and intersects their pure-X / pure-Z eligibility sets
while searching all six local label maps per site. No projector is used by the
oracle. Pruning is justified because each additional assigned qubit can only
shrink either intersection, so cannot increase the sum of its dimensions.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import report_output, save_report


import itertools
import json
import random
import time

SEED = 0x20261008
RNG = random.Random(SEED)
START = time.monotonic()
OUT = None


def basis(rows, width):
    a = [int(x) for x in rows if x]
    rank = 0
    piv = []
    for col in range(width):
        k = next((k for k in range(rank, len(a)) if a[k] >> col & 1), None)
        if k is None:
            continue
        a[rank], a[k] = a[k], a[rank]
        for j in range(len(a)):
            if j != rank and a[j] >> col & 1:
                a[j] ^= a[rank]
        piv.append(col)
        rank += 1
    return a[:rank], piv


def kernel(rows, width):
    rr, piv = basis(rows, width)
    return [(1 << j) | sum(((row >> j) & 1) << p for row, p in zip(rr, piv))
            for j in range(width) if j not in piv]


def affine_solve(equations, width):
    # Incremental triangular elimination; free variables set to zero.
    triangular = {}
    for lhs, rhs in equations:
        while lhs:
            p = (lhs & -lhs).bit_length() - 1
            if p in triangular:
                old, bit = triangular[p]
                lhs ^= old
                rhs ^= bit
            else:
                triangular[p] = (lhs, rhs)
                break
        if not lhs and rhs:
            return None
    solution = 0
    for p in sorted(triangular, reverse=True):
        lhs, rhs = triangular[p]
        if ((lhs & solution).bit_count() & 1) ^ rhs:
            solution |= 1 << p
    assert all(((a & solution).bit_count() & 1) == b for a, b in equations)
    return solution


def components(rows, n):
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for row in rows:
        pivot_site = ((row & -row).bit_length() - 1) // 2
        for i in range(n):
            if row >> (2 * i) & 3:
                parent[root(i)] = root(pivot_site)
    return [[i for i in range(n) if root(i) == j]
            for j in range(n) if root(j) == j]


def manuscript_algorithm(input_rows, n, anchor_last=False):
    rows, _ = basis(input_rows, 2 * n)
    maps = [None] * n
    used_branches = []
    parts = components(rows, n)
    for sites in parts:
        local = [sum(((g >> (2 * site)) & 3) << (2 * j)
                     for j, site in enumerate(sites))
                 for g in rows if ((g & -g).bit_length() - 1) // 2 in sites]
        N = len(sites)
        annihilator = kernel(local, 2 * N)
        equations = []
        for g in local:
            for h in annihilator:
                lhs = rhs = 0
                for i in range(N):
                    x, z = g >> (2 * i) & 1, g >> (2 * i + 1) & 1
                    u, v = h >> (2 * i) & 1, h >> (2 * i + 1) & 1
                    lhs |= ((x * u) ^ (z * v)) << (3 * i)
                    lhs |= (x * v) << (3 * i + 1)
                    lhs |= (z * u) << (3 * i + 2)
                    rhs ^= z * v
                equations.append((lhs, rhs))
        anchor = N - 1 if anchor_last else 0
        for branch in (1, 2):
            solution = affine_solve(equations + [(1 << (3 * anchor + branch), 0)], 3 * N)
            if solution is not None:
                used_branches.append(branch)
                break
        if solution is None:
            return False, None, {"failed_component": sites, "branches": used_branches}
        for j, site in enumerate(sites):
            a, b, c = ((solution >> (3 * j + k)) & 1 for k in range(3))
            P = ((a, b), (c, 1 ^ a))
            def column_action(v):
                return ((P[0][0] * v[0]) ^ (P[0][1] * v[1]),
                        (P[1][0] * v[0]) ^ (P[1][1] * v[1]))
            image = next(v for v in ((1, 0), (0, 1), (1, 1)) if column_action(v) == v)
            null = next(v for v in ((1, 0), (0, 1), (1, 1)) if column_action(v) == (0, 0))
            maps[site] = (image[0], null[0], image[1], null[1])
    return True, maps, {"components": parts, "branches": used_branches}


# Oracle-specific routines use no affine solver, annihilator, decomposition, or
# projector. Matrices (a,b,c,d) act on ROW labels (x,z) by (ax+cz,bx+dz).
LOCAL_MAPS = tuple(q for q in itertools.product((0, 1), repeat=4)
                  if (q[0] * q[3]) ^ (q[1] * q[2]))


def rowspace_elements(independent_rows):
    elements = [0]
    for row in independent_rows:
        elements += [element ^ row for element in elements]
    assert len(set(elements)) == len(elements)
    return elements


def actual_row_action(rows, local_maps):
    transformed = []
    for row in rows:
        out = 0
        for i, (a, b, c, d) in enumerate(local_maps):
            x, z = row >> (2 * i) & 1, row >> (2 * i + 1) & 1
            out |= ((x * a) ^ (z * c)) << (2 * i)
            out |= ((x * b) ^ (z * d)) << (2 * i + 1)
        transformed.append(out)
    return transformed


def direct_css_dimensions(independent_rows, n):
    elements = rowspace_elements(independent_rows)
    xmask = sum(1 << (2 * i) for i in range(n))
    zmask = xmask << 1
    pure_x = [g for g in elements if not g & zmask]
    pure_z = [g for g in elements if not g & xmask]
    dx, dz = len(pure_x).bit_length() - 1, len(pure_z).bit_length() - 1
    assert len(pure_x) == 1 << dx and len(pure_z) == 1 << dz
    return dx, dz


def exhaustive_local_oracle(independent_rows, n):
    rank = len(independent_rows)
    elements = rowspace_elements(independent_rows)
    all_elements = (1 << len(elements)) - 1
    filters = []
    for i in range(n):
        site_filters = []
        for a, b, c, d in LOCAL_MAPS:
            pure_x = pure_z = 0
            for j, element in enumerate(elements):
                x, z = element >> (2 * i) & 1, element >> (2 * i + 1) & 1
                transformed_x = (x * a) ^ (z * c)
                transformed_z = (x * b) ^ (z * d)
                if transformed_z == 0:
                    pure_x |= 1 << j
                if transformed_x == 0:
                    pure_z |= 1 << j
            site_filters.append((pure_x, pure_z))
        filters.append(site_filters)
    visited = leaves = 0
    def search(i, candidates_x, candidates_z, choices):
        nonlocal visited, leaves
        visited += 1
        if i == n:
            leaves += 1
            dx = candidates_x.bit_count().bit_length() - 1
            dz = candidates_z.bit_count().bit_length() - 1
            if dx + dz == rank:
                return choices
            return None
        for q, (mask_x, mask_z) in zip(LOCAL_MAPS, filters[i]):
            next_x, next_z = candidates_x & mask_x, candidates_z & mask_z
            dx, dz = next_x.bit_count().bit_length() - 1, next_z.bit_count().bit_length() - 1
            if dx + dz < rank:
                continue
            ans = search(i + 1, next_x, next_z, choices + [q])
            if ans is not None:
                return ans
        return None
    result = search(0, all_elements, all_elements, [])
    if result is not None:
        assert sum(direct_css_dimensions(actual_row_action(independent_rows, result), n)) == rank
    return result is not None, {"nodes": visited, "full_assignments": leaves}


def literal_unpruned_oracle(independent_rows, n):
    """Second oracle: materialize full transformed spaces, with no pruning."""
    full_assignments = 0
    for maps in itertools.product(LOCAL_MAPS, repeat=n):
        full_assignments += 1
        transformed = actual_row_action(independent_rows, maps)
        if sum(direct_css_dimensions(transformed, n)) == len(independent_rows):
            return True, full_assignments
    return False, full_assignments


def symp(a, b, n):
    result = 0
    for i in range(n):
        result ^= ((a >> (2 * i) & 1) & (b >> (2 * i + 1) & 1))
        result ^= ((a >> (2 * i + 1) & 1) & (b >> (2 * i) & 1))
    return result


def all_isotropic(n):
    width = 2 * n
    for r in range(n + 1):
        for pivots in itertools.combinations(range(width), r):
            free_positions = [(i, j) for i, p in enumerate(pivots)
                              for j in range(p + 1, width) if j not in pivots]
            for assignment in range(1 << len(free_positions)):
                rows = [1 << p for p in pivots]
                for bit, (i, j) in enumerate(free_positions):
                    rows[i] |= ((assignment >> bit) & 1) << j
                if all(not symp(a, b, n) for i, a in enumerate(rows) for b in rows[i + 1:]):
                    yield rows


def random_isotropic(n, rank):
    rows, elements = [], {0}
    while len(rows) < rank:
        g = RNG.randrange(1, 1 << (2 * n))
        if g not in elements and all(not symp(g, h, n) for h in rows):
            rows.append(g)
            elements |= {g ^ h for h in list(elements)}
    return rows


def rebase_with_redundancy(rows):
    out = rows[:]
    if len(out) > 1:
        for _ in range(6 * len(out)):
            a, b = RNG.sample(range(len(out)), 2)
            out[a] ^= out[b]
        RNG.shuffle(out)
    out += [0, 0]
    for _ in range(len(rows) + 2):
        combination = 0
        for row in rows:
            if RNG.randrange(2):
                combination ^= row
        out.append(combination)
    RNG.shuffle(out)
    return out


def hidden_css(n, rank):
    rx = RNG.randrange(rank + 1)
    xs = []
    while len(xs) < rx:
        trial = RNG.randrange(1, 1 << n)
        new, _ = basis(xs + [trial], n)
        if len(new) > len(xs):
            xs = new
    kz = kernel(xs, n)
    zs = []
    while len(zs) < rank - rx:
        trial = 0
        for row in kz:
            if RNG.randrange(2):
                trial ^= row
        new, _ = basis(zs + [trial], n)
        if len(new) > len(zs):
            zs = new
    css = [sum(((g >> i) & 1) << (2 * i) for i in range(n)) for g in xs]
    css += [sum(((g >> i) & 1) << (2 * i + 1) for i in range(n)) for g in zs]
    maps = [RNG.choice(LOCAL_MAPS) for _ in range(n)]
    return actual_row_action(css, maps)


REPORT = {"seed": SEED, "implementation": "Independent row-action implementation",
          "oracle": "All six local GL(2,2) matrices; actual rowspace pure-X/Z intersections; monotone-intersection pruning",
          "coverage": {}, "mismatches": [], "witness_failures": [], "cases": 0,
          "positive": 0, "negative": 0, "redundant_rebased_checks": 0,
          "alternate_anchor_checks": 0, "oracle_nodes": 0, "oracle_full_assignments": 0,
          "literal_unpruned_cases": 0, "literal_unpruned_negative_cases": 0,
          "literal_unpruned_assignments": 0,
          "branch_c_used": 0, "examples": []}


def record_case(rows, n, category, known=None, rebase=True, alternate=False):
    independent, _ = basis(rows, 2 * n)
    assert len(independent) <= n
    assert all(not symp(g, h, n) for i, g in enumerate(independent) for h in independent[i + 1:])
    actual, witness, detail = manuscript_algorithm(rows, n)
    expected, stats = exhaustive_local_oracle(independent, n)
    REPORT["cases"] += 1
    REPORT["positive" if actual else "negative"] += 1
    REPORT["oracle_nodes"] += stats["nodes"]
    REPORT["oracle_full_assignments"] += stats["full_assignments"]
    REPORT["branch_c_used"] += detail["branches"].count(2)
    cov = REPORT["coverage"].setdefault(category, {})
    key = f"n{n}_r{len(independent)}"
    cov[key] = cov.get(key, 0) + 1
    if actual != expected or (known is not None and actual != known):
        REPORT["mismatches"].append({"category": category, "n": n, "rows": rows,
                                      "algorithm": actual, "oracle": expected, "known": known})
    if actual:
        assert all((a * d) ^ (b * c) for a, b, c, d in witness)
        transformed = actual_row_action(independent, witness)
        dx, dz = direct_css_dimensions(transformed, n)
        if dx + dz != len(independent):
            REPORT["witness_failures"].append({"n": n, "rows": rows, "maps": witness, "dimensions": [dx, dz]})
        # Independently confirm pure-coordinate projected checks remain in the
        # transformed space and that the two classical check spaces commute.
        elements = set(rowspace_elements(transformed))
        xm = sum(1 << (2 * i) for i in range(n))
        xs, zs = [g & xm for g in transformed], [g & (xm << 1) for g in transformed]
        assert all(g in elements for g in xs + zs)
        assert all(not symp(x, z, n) for x in xs for z in zs)
    if rebase:
        rebased = rebase_with_redundancy(independent)
        answer2, maps2, _ = manuscript_algorithm(rebased, n)
        assert answer2 == actual, ("generating-list dependence", n, rows, rebased)
        if answer2:
            assert sum(direct_css_dimensions(actual_row_action(independent, maps2), n)) == len(independent)
        REPORT["redundant_rebased_checks"] += 1
    if alternate:
        answer3, maps3, _ = manuscript_algorithm(rows, n, anchor_last=True)
        assert answer3 == actual, ("anchor dependence", n, rows)
        if answer3:
            assert sum(direct_css_dimensions(actual_row_action(independent, maps3), n)) == len(independent)
        REPORT["alternate_anchor_checks"] += 1
    if (not actual and len(REPORT["examples"]) < 10) or (2 in detail["branches"] and not any(e.get("branch_c") for e in REPORT["examples"])):
        REPORT["examples"].append({"category": category, "n": n, "rank": len(independent),
                                   "rows": rows, "answer": actual, "branch_c": 2 in detail["branches"],
                                   "maps": witness, "detail": detail})
    return actual


def checkpoint(stage):
    REPORT["elapsed_seconds"] = round(time.monotonic() - START, 3)
    REPORT["stage"] = stage
    REPORT["all_passed"] = not REPORT["mismatches"] and not REPORT["witness_failures"]
    if stage == 'complete':
        save_report(REPORT, OUT)
    print(json.dumps({k: REPORT[k] for k in ("stage", "cases", "positive", "negative", "elapsed_seconds", "all_passed")}), flush=True)


def main():
    global OUT
    OUT = report_output('exhaustive.json')
    for n, expected_count in enumerate((1, 4, 31, 514, 19381)):
        count = 0
        seen = set()
        for rows in all_isotropic(n):
            frozen = tuple(rows)
            assert frozen not in seen
            seen.add(frozen)
            answer = record_case(rows, n, "exhaustive", alternate=True)
            # Recheck every small negative with literal 6**n enumeration, all
            # spaces through n=3, and deterministic positive samples at n=4.
            if n <= 3 or not answer or count % 97 == 0:
                literal, assignments = literal_unpruned_oracle(rows, n)
                assert literal == answer, ("pruned/literal oracle disagreement", n, rows)
                REPORT["literal_unpruned_cases"] += 1
                REPORT["literal_unpruned_negative_cases"] += int(not literal)
                REPORT["literal_unpruned_assignments"] += assignments
            count += 1
        assert count == expected_count, (n, count, expected_count)
        checkpoint(f"exhaustive_n{n}_complete")
    for n in range(5, 9):
        for rank in range(n + 1):
            for _ in range(20):
                record_case(random_isotropic(n, rank), n, "seeded_isotropic", alternate=True)
            for _ in range(10):
                record_case(hidden_css(n, rank), n, "hidden_css", known=True, alternate=True)
        checkpoint(f"random_and_hidden_n{n}_complete")
    # Tensor factors include known obstructed five-cycle graph states, hidden
    # CSS, rank-zero sites, and partial stabilizers with independent row mixing.
    cycle5 = [(1 << (2 * i)) | (1 << (2 * ((i - 1) % 5) + 1)) |
              (1 << (2 * ((i + 1) % 5) + 1)) for i in range(5)]
    for n in range(5, 9):
        for _ in range(20):
            right_n = n - 5
            rr = RNG.randrange(right_n + 1)
            rows = cycle5 + [g << 10 for g in hidden_css(right_n, rr)]
            maps = [RNG.choice(LOCAL_MAPS) for _ in range(n)]
            record_case(actual_row_action(rows, maps), n, "obstructed_tensor", known=False, alternate=True)
        for _ in range(20):
            split = RNG.randrange(1, n)
            left = random_isotropic(split, RNG.randrange(split + 1))
            right_n = n - split
            right = random_isotropic(right_n, RNG.randrange(right_n + 1))
            rows = left + [g << (2 * split) for g in right]
            record_case(rows, n, "mixed_tensor", alternate=True)
        checkpoint(f"tensor_n{n}_complete")
    checkpoint("complete")
    assert REPORT["all_passed"]


if __name__ == "__main__":
    main()
