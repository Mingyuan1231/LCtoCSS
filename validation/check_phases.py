"""Check Clifford direction, signed Pauli groups, and local phase correction."""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.reproducibility import report_output, save_report


from collections import Counter, deque
from hashlib import sha256
import json
import random

import numpy as np

import validation_reference as v
SOURCE = Path(v.__file__)
PAULI = (np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], dtype=complex),
         np.diag([1, -1]).astype(complex), np.array([[0, -1j], [1j, 0]]))


def cliffords():
    # Pauli labels: I=0, X=1, Z=2, Y=3; sign is an exponent of -1.
    ident = ((0, 0), (1, 0), (2, 0), (3, 0))
    h = ((0, 0), (2, 0), (1, 0), (3, 1))
    s = ((0, 0), (3, 0), (2, 0), (1, 1))
    queue, seen = deque([ident]), {ident}
    unitaries = {ident: PAULI[0]}
    gates = (np.array([[1, 1], [1, -1]]) / np.sqrt(2), np.diag([1, 1j]))
    while queue:
        old = queue.popleft()
        for gate, unitary in zip((h, s), gates):
            new = tuple((gate[label][0], sign ^ gate[label][1]) for label, sign in old)
            if new not in seen:
                seen.add(new)
                queue.append(new)
                unitaries[new] = unitary @ unitaries[old]
    result = {}
    for table in sorted(seen):
        x, z = table[1][0], table[2][0]
        matrix = (x & 1, z & 1, x >> 1, z >> 1)
        result.setdefault(matrix, []).append(table)
    assert len(seen) == 24 and len(result) == 6
    assert all(len(lifts) == 4 for lifts in result.values())
    for table, unitary in unitaries.items():
        for label in range(4):
            target, sign = table[label]
            assert np.allclose(unitary @ PAULI[label] @ unitary.conj().T,
                               (-1)**sign * PAULI[target], atol=1e-12)
    return result, unitaries


def tensor(factors):
    out = np.ones((1, 1), dtype=complex)
    for factor in factors:
        out = np.kron(out, factor)
    return out


def pauli_matrix(row, n):
    return tensor(PAULI[(row >> (2*i)) & 3] for i in range(n))


def numerical_signed_audit(lifts, unitaries):
    rng = random.Random(202610081)
    counts = Counter()
    largest_error = 0.0
    for n in range(5):
        for rank in range(n+1):
            for sample in range(4):
                rows = v.random_isotropic(n, rank, rng)
                answer = v.projector_decide(rows + [0] + rows, n)
                counts["candidate_groups"] += 1
                if not answer["equivalent"]:
                    counts["non_css_groups_skipped"] += 1
                    continue
                signs = rng.randrange(1 << rank)
                group = signed_group(rows, signs, n)
                tables = [lifts[tuple(matrix)][0] for matrix in answer["local_matrices"]]
                unitary = tensor(unitaries[table] for table in tables)
                changed = conjugate(group, tables, n)
                # Every group element checked as an actual 2^n by 2^n operator.
                for source, phase in group.items():
                    target_group = conjugate({source: phase}, tables, n)
                    target, target_phase = next(iter(target_group.items()))
                    actual = unitary @ ((1j**phase) * pauli_matrix(source, n)) @ unitary.conj().T
                    expected = (1j**target_phase) * pauli_matrix(target, n)
                    error = float(np.max(np.abs(actual - expected)))
                    largest_error = max(largest_error, error)
                    assert error < 1e-12
                    counts["actual_operator_conjugations"] += 1
                pure = answer["css_generators"]
                # Construct precisely C Omega t = epsilon, then solve it.
                equations = []
                for row in pure:
                    coefficients = sum(v.symplectic(row, 1 << j, n) << j for j in range(2*n))
                    equations.append(coefficients | ((changed[row] // 2) << (2*n)))
                solution = v.linear_solve(equations, 2*n)
                assert solution is not None
                correction = solution[0]
                correction_matrix = pauli_matrix(correction, n)
                for row in pure:
                    signed = (1j**changed[row]) * pauli_matrix(row, n)
                    actual = correction_matrix @ signed @ correction_matrix.conj().T
                    expected = pauli_matrix(row, n)
                    error = float(np.max(np.abs(actual - expected)))
                    largest_error = max(largest_error, error)
                    assert error < 1e-12
                    counts["actual_all_plus_generator_corrections"] += 1
                counts["accepted_signed_groups"] += 1
    return {"counts": dict(counts), "max_absolute_matrix_error": largest_error,
            "seed": 202610081, "qubit_range": [0, 4], "samples_per_rank": 4}


def multiply(left, right, n):
    a, phase_a = left
    b, phase_b = right
    phase = phase_a + phase_b
    for i in range(n):
        x, z = (a >> (2*i)) & 1, (a >> (2*i+1)) & 1
        u, w = (b >> (2*i)) & 1, (b >> (2*i+1)) & 1
        phase += x*z + u*w - (x ^ u)*(z ^ w) + 2*z*u
    return a ^ b, phase % 4


def signed_group(rows, signs, n):
    group = {0: 0}
    for j, row in enumerate(rows):
        generator = (row, 2*((signs >> j) & 1))
        for old in list(group.items()):
            label, phase = multiply(old, generator, n)
            assert phase in (0, 2)
            assert label not in group
            group[label] = phase
    assert group[0] == 0
    return group


def conjugate(group, tables, n):
    changed = {}
    for row, phase in group.items():
        new = 0
        for i in range(n):
            label, sign = tables[i][(row >> (2*i)) & 3]
            new |= label << (2*i)
            phase += 2*sign
        assert new not in changed
        changed[new] = phase % 4
    return changed


def main():
    report_path = report_output('phases.json')
    lifts, unitaries = cliffords()
    counts = Counter()
    for n in range(4):
        for rows in v.all_isotropic_spaces(n):
            answer = v.projector_decide(list(rows), n)
            counts["binary_spaces"] += 1
            if not answer["equivalent"]:
                continue
            matrices = answer["local_matrices"]
            tables = [lifts[tuple(matrix)][0] for matrix in matrices]
            pure_rows = answer["css_generators"]
            xmask = sum(1 << (2*i) for i in range(n))
            assert all(not (row & xmask) or not (row & (xmask << 1)) for row in pure_rows)
            for signs in range(1 << len(rows)):
                group = signed_group(rows, signs, n)
                changed = conjugate(group, tables, n)
                assert all(row in changed for row in pure_rows)
                epsilon = [changed[row] // 2 for row in pure_rows]
                correction = next(t for t in range(1 << (2*n)) if all(
                    v.symplectic(row, t, n) == sign
                    for row, sign in zip(pure_rows, epsilon)))
                assert all((changed[row] + 2*v.symplectic(row, correction, n)) % 4 == 0
                           for row in pure_rows)
                counts["signed_groups_and_pauli_corrections"] += 1
    example = v.projector_decide([11], 2)
    transposed = [(a, c, b, d) for a, b, c, d in example["local_matrices"]]
    wrong = [v.transform(11, transposed)]
    assert wrong == [6] and not v.css_rank(tuple(wrong), 2)
    report = {
        "validator_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "counts": dict(counts),
        "physical_clifford_label_maps": len(lifts),
        "physical_signed_clifford_maps": sum(map(len, lifts.values())),
        "actual_numpy_operator_audit": numerical_signed_audit(lifts, unitaries),
        "phase_convention": "Hermitian P(x,z) = product_i i^(x_i*z_i) X_i^x_i Z_i^z_i; phases modulo 4",
        "orientation_counterexample": {
            "input": "YZ", "input_integer": 11,
            "validator_column_matrices": example["local_matrices"],
            "correct_column_action": "XX",
            "incorrect_untransposed_row_action": "ZX",
            "required_row_action_conversion": "transpose every returned local matrix",
        },
        "scope": "Finite sign and witness checks, all isotropic spaces and every valid independent sign assignment through n=3; no analytic proof or asymptotic timing inference",
    }
    save_report(report, report_path)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
