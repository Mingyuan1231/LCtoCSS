"""Independent LC-to-CSS validation using only the Python standard library.

Run from the repository root:
    python validation/validation_reference.py

Pauli vectors use interleaved bits (x_0,z_0,x_1,z_1,...).
Local matrices act on column labels; transpose them for right row action.  The oracle
enumerates the six invertible local binary matrices, pruning only when the
monotone sum of the X- and Z-column ranks already exceeds the stabilizer rank.
It does not use the projector criterion or the factor decomposition.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.reproducibility import output_file, utc_now


import argparse
import hashlib
import itertools
import json
import random
import time
from collections import Counter


def rref(rows: list[int] | tuple[int, ...], width: int) -> tuple[int, ...]:
    work = [row for row in rows if row]
    rank = 0
    for column in range(width):
        pivot = next((j for j in range(rank, len(work)) if work[j] >> column & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for j in range(len(work)):
            if j != rank and work[j] >> column & 1:
                work[j] ^= work[rank]
        rank += 1
    return tuple(work[:rank])


def linear_solve(equations: list[int], width: int) -> tuple[int, tuple[int, ...]] | None:
    reduced = rref(equations, width + 1)
    if 1 << width in reduced:
        return None
    pivots = [(row & -row).bit_length() - 1 for row in reduced]
    pivot_set = set(pivots)
    particular = sum(((row >> width) & 1) << pivot for row, pivot in zip(reduced, pivots))
    null_basis = []
    for free in range(width):
        if free not in pivot_set:
            null_basis.append((1 << free) | sum(((row >> free) & 1) << pivot for row, pivot in zip(reduced, pivots)))
    return particular, tuple(null_basis)


def annihilator(rows: tuple[int, ...], width: int) -> tuple[int, ...]:
    solution = linear_solve(list(rows), width)
    assert solution is not None
    return solution[1]


def symplectic(left: int, right: int, n: int) -> int:
    return sum(((left >> (2 * i) & 1) & (right >> (2 * i + 1) & 1)) ^
               ((left >> (2 * i + 1) & 1) & (right >> (2 * i) & 1)) for i in range(n)) & 1


def contains(rows: tuple[int, ...], vector: int, width: int) -> bool:
    return len(rref((*rows, vector), width)) == len(rows)


def restrict(row: int, qubits: list[int]) -> int:
    return sum(((row >> (2 * old)) & 3) << (2 * new) for new, old in enumerate(qubits))


def rref_support_factors(rows: tuple[int, ...], n: int) -> list[list[int]]:
    """Separate decomposition check using the qubit supports of RREF rows."""
    parents = list(range(n))

    def root(i: int) -> int:
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    for row in rows:
        pivot_qubit = ((row & -row).bit_length() - 1) // 2
        for qubit in range(n):
            if row >> (2 * qubit) & 3:
                parents[root(qubit)] = root(pivot_qubit)
    classes: dict[int, list[int]] = {}
    for qubit in range(n):
        classes.setdefault(root(qubit), []).append(qubit)
    return list(classes.values())


def scalar_factors(rows: tuple[int, ...], n: int) -> list[list[int]]:
    """Compute coordinate-mask algebra and its evaluation classes."""
    equations = []
    for row in rows:
        for check in annihilator(rows, 2 * n):
            equations.append(sum((((row >> (2 * i)) & (check >> (2 * i)) & 3).bit_count() & 1) << i for i in range(n)))
    solution = linear_solve(equations, n)
    assert solution is not None
    mask_basis = solution[1]
    classes: dict[tuple[int, ...], list[int]] = {}
    for i in range(n):
        classes.setdefault(tuple(mask >> i & 1 for mask in mask_basis), []).append(i)
    factors = list(classes.values())
    assert factors == rref_support_factors(rows, n)
    total_rank = 0
    for factor in factors:
        coordinate_mask = sum(3 << (2 * i) for i in factor)
        projected = rref([row & coordinate_mask for row in rows], 2 * n)
        assert all(contains(rows, row, 2 * n) for row in projected)
        total_rank += len(projected)
    assert total_rank == len(rows)
    return factors


def preserving_equations(rows: tuple[int, ...], n: int) -> list[int]:
    equations = []
    for row in rows:
        for check in annihilator(rows, 2 * n):
            coefficients = 0
            constant = 0
            for i in range(n):
                x, z = (row >> (2 * i)) & 1, (row >> (2 * i + 1)) & 1
                hx, hz = (check >> (2 * i)) & 1, (check >> (2 * i + 1)) & 1
                coefficients |= ((hx * x) ^ (hz * z)) << (3 * i)
                coefficients |= (hx * z) << (3 * i + 1)
                coefficients |= (hz * x) << (3 * i + 2)
                constant ^= hz * z
            equations.append(coefficients | (constant << (3 * n)))
    return equations


def matrix_apply(matrix: tuple[int, int, int, int], vector: int) -> int:
    a, b, c, d = matrix
    x, z = vector & 1, vector >> 1 & 1
    return ((a * x) ^ (b * z)) | (((c * x) ^ (d * z)) << 1)


def transform(row: int, matrices: list[tuple[int, int, int, int]]) -> int:
    return sum(matrix_apply(matrix, (row >> (2 * i)) & 3) << (2 * i) for i, matrix in enumerate(matrices))


def css_rank(rows: tuple[int, ...], n: int) -> bool:
    xmask = sum(1 << (2 * i) for i in range(n))
    zmask = xmask << 1
    return len(rref([row & xmask for row in rows], 2 * n)) + len(rref([row & zmask for row in rows], 2 * n)) == len(rows)


def projector_decide(raw_rows: list[int], n: int) -> dict:
    rows = rref(raw_rows, 2 * n)
    assert all(symplectic(a, b, n) == 0 for a in rows for b in rows)
    factors = scalar_factors(rows, n)
    projectors: list[tuple[int, int, int, int] | None] = [None] * n
    branches = []
    for factor in factors:
        local_rows = rref([restrict(row, factor) for row in rows], 2 * len(factor))
        equations = preserving_equations(local_rows, len(factor))
        witness = None
        for root_variable, name in ((1, "b_root=0"), (2, "c_root=0")):
            witness = linear_solve(equations + [1 << root_variable], 3 * len(factor))
            if witness is not None:
                branches.append(name)
                break
        if witness is None:
            return {"equivalent": False, "factors": factors, "failed_factor": factor}
        assignment = witness[0]
        local_projectors = []
        for local_i, global_i in enumerate(factor):
            a, b, c = ((assignment >> (3 * local_i + j)) & 1 for j in range(3))
            assert b * c == 0
            projector = (a, b, c, a ^ 1)
            local_projectors.append(projector)
            projectors[global_i] = projector
        assert all(contains(local_rows, transform(row, local_projectors), 2 * len(factor)) for row in local_rows)
    assert all(projector is not None for projector in projectors)
    matrices = []
    for projector in projectors:
        assert projector is not None
        image_vector = next(v for v in (1, 2, 3) if matrix_apply(projector, v) == v)
        kernel_vector = next(v for v in (1, 2, 3) if matrix_apply(projector, v) == 0)
        ux, uz = image_vector & 1, image_vector >> 1
        vx, vz = kernel_vector & 1, kernel_vector >> 1
        matrix = (vz, vx, uz, ux)
        assert (matrix[0] * matrix[3]) ^ (matrix[1] * matrix[2]) == 1
        assert matrix_apply(matrix, image_vector) == 1
        assert matrix_apply(matrix, kernel_vector) == 2
        matrices.append(matrix)
    changed = rref([transform(row, matrices) for row in rows], 2 * n)
    assert len(changed) == len(rows) and css_rank(changed, n)
    xmask = sum(1 << (2 * i) for i in range(n))
    pure_rows = rref([row & mask for row in changed for mask in (xmask, xmask << 1)], 2 * n)
    assert len(pure_rows) == len(rows)
    assert all(contains(changed, row, 2 * n) for row in pure_rows)
    return {"equivalent": True, "factors": factors, "branches": branches,
            "local_matrices": matrices, "css_generators": list(pure_rows)}


def extend_basis(basis: tuple[int, ...], vector: int) -> tuple[tuple[int, ...], int]:
    for pivot in basis:
        vector = min(vector, vector ^ pivot)
    if vector == 0:
        return basis, 0
    return tuple(sorted((*basis, vector), reverse=True)), 1


def exhaustive_local_oracle(raw_rows: list[int], n: int) -> dict:
    """Exact six-choice-per-qubit search; pruning uses the rank CSS test only."""
    rows = rref(raw_rows, 2 * n)
    rank = len(rows)
    columns = []
    for i in range(n):
        x = sum(((row >> (2 * i)) & 1) << j for j, row in enumerate(rows))
        z = sum(((row >> (2 * i + 1)) & 1) << j for j, row in enumerate(rows))
        columns.append((x, z, x ^ z))
    counts = {"branches": 0, "rank_pruned": 0, "accepted_leaves": 0}
    choices = tuple(itertools.permutations(range(3), 2))

    def visit(i: int, xb: tuple[int, ...], zb: tuple[int, ...], xr: int, zr: int) -> bool:
        if i == n:
            assert xr + zr == rank
            counts["accepted_leaves"] += 1
            return True
        for x_index, z_index in choices:
            counts["branches"] += 1
            xn, dx = extend_basis(xb, columns[i][x_index])
            zn, dz = extend_basis(zb, columns[i][z_index])
            if xr + dx + zr + dz > rank:
                counts["rank_pruned"] += 1
                continue
            if visit(i + 1, xn, zn, xr + dx, zr + dz):
                return True
        return False

    answer = visit(0, (), (), 0, 0)
    return {"equivalent": answer, **counts}


def all_isotropic_spaces(n: int):
    by_rank: set[tuple[int, ...]] = {()}
    yield ()
    for rank in range(1, n + 1):
        next_rank = set()
        for rows in by_rank:
            for vector in range(1, 1 << (2 * n)):
                if all(symplectic(vector, row, n) == 0 for row in rows):
                    larger = rref((*rows, vector), 2 * n)
                    if len(larger) == rank:
                        next_rank.add(larger)
        yield from sorted(next_rank)
        by_rank = next_rank


def random_isotropic(n: int, rank: int, rng: random.Random) -> list[int]:
    rows: tuple[int, ...] = ()
    while len(rows) < rank:
        vector = rng.randrange(1, 1 << (2 * n))
        if all(symplectic(vector, row, n) == 0 for row in rows):
            candidate = rref((*rows, vector), 2 * n)
            if len(candidate) > len(rows):
                rows = candidate
    return list(rows)


def paulis(word: str) -> int:
    mapping = {"I": 0, "X": 1, "Z": 2, "Y": 3}
    return sum(mapping[letter] << (2 * i) for i, letter in enumerate(word))


def graph_rows(n: int, edges: list[tuple[int, int]]) -> list[int]:
    rows = [1 << (2 * i) for i in range(n)]
    for left, right in edges:
        rows[left] ^= 1 << (2 * right + 1)
        rows[right] ^= 1 << (2 * left + 1)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-per-rank", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261003)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "outputs/validation/reference.json")
    args = parser.parse_args()
    if args.samples_per_rank < 0:
        parser.error("samples-per-rank must be nonnegative")
    rng = random.Random(args.seed)
    start = time.perf_counter()
    records = []
    group_counts: dict[str, Counter] = {}
    rank_counts: dict[str, Counter] = {}
    all_oracle_counts = Counter()

    def check_case(name: str, n: int, supplied: list[int], group: str, expected: bool | None = None) -> None:
        rank = len(rref(supplied, 2 * n))
        projected = projector_decide(supplied, n)
        oracle = exhaustive_local_oracle(supplied, n)
        assert projected["equivalent"] == oracle["equivalent"], (name, n, supplied, projected, oracle)
        if expected is not None:
            assert projected["equivalent"] == expected, (name, projected)
        record = {"name": name, "group": group, "n": n, "rank": rank,
                  "input_rows": supplied, "projector": projected, "oracle": oracle}
        records.append(record)
        count = group_counts.setdefault(group, Counter())
        count["cases"] += 1
        count["css" if projected["equivalent"] else "non_css"] += 1
        rank_counts.setdefault(f"{group}:n={n}", Counter())[str(rank)] += 1
        all_oracle_counts.update({key: value for key, value in oracle.items() if key != "equivalent"})

    expected_counts = {0: 1, 1: 4, 2: 31, 3: 514}
    for n in range(4):
        count = 0
        for index, rows in enumerate(all_isotropic_spaces(n)):
            check_case(f"all_n{n}_{index}", n, list(rows), "exhaustive_spaces")
            count += 1
        assert count == expected_counts[n]
        print(f"Exhaustive n={n}: {count} isotropic spaces passed", flush=True)

    for n in range(4, 7):
        for rank in range(n + 1):
            for sample in range(args.samples_per_rank):
                rows = random_isotropic(n, rank, rng)
                # Every sampled input includes redundant/zero rows.
                supplied = rows + [0]
                if rows:
                    supplied += [rows[0], rows[-1], rows[0] ^ rows[-1]]
                rng.shuffle(supplied)
                check_case(f"random_n{n}_r{rank}_{sample}", n, supplied, "random_with_redundancy")
        print(f"Sampled n={n}: all ranks and redundant rows passed", flush=True)

    c5 = graph_rows(5, [(i, (i + 1) % 5) for i in range(5)])
    wheel = graph_rows(6, [(i, (i + 1) % 5) for i in range(5)] + [(i, 5) for i in range(5)])
    noncss4 = [paulis(word) for word in ("ZIXZ", "YXYI", "IZZX")]
    bell = [paulis("XX"), paulis("ZZ")]
    examples = [
        ("C5", 5, c5, False),
        ("wheel_C5_plus_hub", 6, wheel, False),
        ("Cross_Vandeth_4_qubit_nonCSS", 4, noncss4, False),
        ("C5_plus_unused_qubit", 6, c5, False),
        ("nonCSS4_tensor_Bell", 6, noncss4 + [row << 8 for row in bell], False),
        ("three_Bell_pairs", 6, [row << (4 * offset) for offset in range(3) for row in bell], True),
        ("one_Bell_and_four_unused", 6, bell, True),
        ("six_Y_product_state", 6, [3 << (2 * i) for i in range(6)], True),
        ("one_Y_and_five_unused", 6, [3], True),
        ("fully_unconstrained_six_qubits", 6, [], True),
    ]
    for name, n, rows, expected in examples:
        check_case(name, n, rows + [0] + rows, "structured_cases", expected)

    report = {
        "executed_at_utc": utc_now(),
        "description": "Independent general stabilizer-code LC-to-CSS criterion versus exhaustive local-Clifford rank oracle",
        "seed": args.seed,
        "samples_per_rank": args.samples_per_rank,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.perf_counter() - start, 6),
        "total_cases": len(records),
        "mismatches": 0,
        "positive_witnesses_verified": sum(record["projector"]["equivalent"] for record in records),
        "decomposition_cross_checks": "Scalar-mask evaluation classes agree with connected components of RREF row supports for every case.",
        "group_counts": {group: dict(counter) for group, counter in group_counts.items()},
        "rank_counts": {group: dict(counter) for group, counter in rank_counts.items()},
        "oracle_counts": dict(all_oracle_counts),
        "oracle_pruning": "Complete six-local-choice search; a subtree is excluded only when rank(X)+rank(Z)>rank(S), which cannot decrease by adding columns.",
        "positive_witness_checks": ["Each local matrix is invertible", "Each projector preserves the factor", "Every local projector has trace one and determinant zero", "Transformed group satisfies the CSS rank equation", "Pure-X and pure-Z generators span exactly the transformed group"],
        "cases": records,
    }
    output_file(args.output).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "cases"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
