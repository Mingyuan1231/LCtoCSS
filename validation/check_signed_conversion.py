"""Check appendix sign recovery and syndrome conversion with Pauli matrices.

These finite checks cover signed output and syndrome conversion.
Run with Python 3.12 and NumPy.
"""

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.reproducibility import report_output, save_report, utc_now

from itertools import product
import hashlib
import json

import numpy as np


HERE = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
Y = 1j * X @ Z
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)
GATES = (I, H, S, H @ S, S @ H, H @ S @ H)


def tensor(matrices):
    result = np.ones((1, 1), dtype=complex)
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def pauli(label):
    return tensor({(0, 0): I, (1, 0): X, (0, 1): Z, (1, 1): Y}[
        tuple(label[i:i + 2])] for i in range(0, len(label), 2))


def label_of(word):
    return np.array([bit for char in word for bit in
                     {'I': (0, 0), 'X': (1, 0), 'Z': (0, 1),
                      'Y': (1, 1)}[char]], dtype=np.uint8)


def close(left, right):
    return np.allclose(left, right, atol=1e-10, rtol=0)


def identify(matrix, labels, matrices):
    for label, candidate in zip(labels, matrices):
        for sign in (0, 1):
            if close(matrix, (-1) ** sign * candidate):
                return label.copy(), sign
    raise AssertionError('Matrix is not a signed Hermitian Pauli')


def multiply_labels(a, sign_a, b, sign_b):
    c = a ^ b
    nu = lambda row: sum(int(row[i]) * int(row[i + 1])
                         for i in range(0, len(row), 2))
    kappa = (nu(a) + nu(b) - nu(c) + 2 * sum(
        int(a[i + 1]) * int(b[i]) for i in range(0, len(a), 2))) % 4
    assert kappa % 2 == 0
    return c, (sign_a + sign_b + kappa // 2) % 2


def row_basis(rows):
    pivots, selected = {}, []
    for index, row in enumerate(rows):
        value = row.copy()
        for pivot, basis_row in sorted(pivots.items()):
            if value[pivot]:
                value ^= basis_row
        if value.any():
            pivots[int(np.flatnonzero(value)[0])] = value
            selected.append(index)
    return selected


def projector(generators, n):
    result = np.eye(2 ** n, dtype=complex)
    for generator in generators:
        result = result @ (np.eye(2 ** n) + generator) / 2
    return result


def main():
    report_path = report_output('signed_conversion.json')
    # Includes rank zero, unused qubits, encoded qubits, and stabilizer states.
    families = [(1, []), (2, ['XI']), (2, ['XX', 'ZZ']),
                (3, ['ZZI', 'IZZ']), (3, ['XXX', 'ZZI', 'IZZ'])]
    cases = syndrome_cases = recovered_rows = 0
    regression_seen = False
    for n, words in families:
        r = len(words)
        labels = [np.array(bits, dtype=np.uint8)
                  for bits in product((0, 1), repeat=2 * n)]
        matrices = [pauli(label) for label in labels]
        omega = np.zeros((2 * n, 2 * n), dtype=np.uint8)
        for i in range(n):
            omega[2 * i, 2 * i + 1] = omega[2 * i + 1, 2 * i] = 1
        mixtures = np.eye(r, dtype=np.uint8)
        for i in range(r - 1):
            mixtures[i] ^= mixtures[i + 1]
        for shift in range(6):
            unitary = tensor(GATES[(shift + i) % 6] for i in range(n))
            # Matrix-derived binary map; independent of the appendix gate table.
            q_rows = []
            for j in range(2 * n):
                unit = np.eye(2 * n, dtype=np.uint8)[j]
                image = unitary @ pauli(unit) @ unitary.conj().T
                q_rows.append(identify(image, labels, matrices)[0])
            q = np.array(q_rows, dtype=np.uint8)
            assert np.array_equal(q @ omega @ q.T % 2, omega)
            for signs in product((0, 1), repeat=r):
                css = [(-1) ** sign * pauli(label_of(word))
                       for word, sign in zip(words, signs)]
                transformed = []
                for row in mixtures:
                    check = np.eye(2 ** n, dtype=complex)
                    for include, generator in zip(row, css):
                        if include:
                            check = check @ generator
                    transformed.append(check)
                inputs = [unitary.conj().T @ check @ unitary
                          for check in transformed]
                g0 = np.array([identify(check, labels, matrices)[0]
                               for check in inputs], dtype=np.uint8).reshape(r, 2 * n)
                t_rows = g0 @ q % 2
                signed_t = [identify(check, labels, matrices)
                            for check in transformed]
                combinations = {}
                for bits in product((0, 1), repeat=r):
                    coefficients = np.array(bits, dtype=np.uint8)
                    row = coefficients @ t_rows % 2
                    direct = np.eye(2 ** n, dtype=complex)
                    folded, sign = np.zeros(2 * n, dtype=np.uint8), 0
                    for include, check, (check_row, check_sign) in zip(
                            bits, transformed, signed_t):
                        if include:
                            direct = direct @ check
                            folded, sign = multiply_labels(
                                folded, sign, check_row, check_sign)
                    assert np.array_equal(folded, row)
                    assert close(direct, (-1) ** sign * pauli(row))
                    combinations[tuple(row)] = (coefficients, sign, direct)
                # All group rows model valid redundant input generators; retain
                # their original projections rather than reduced replacements.
                pure = []
                for row in combinations:
                    for parity in (0, 1):
                        projected = np.zeros(2 * n, dtype=np.uint8)
                        projected[parity::2] = np.array(row)[parity::2]
                        if projected.any():
                            pure.append(projected)
                c = np.array(pure, dtype=np.uint8).reshape(len(pure), 2 * n)
                entries = [combinations[tuple(row)] for row in c]
                k = np.array([entry[0] for entry in entries], dtype=np.uint8).reshape(len(pure), r)
                epsilon = np.array([entry[1] for entry in entries], dtype=np.uint8)
                assert np.array_equal(c, k @ g0 @ q % 2)
                chosen = row_basis(c)
                assert len(chosen) == r
                correction = next(label for label in labels if np.array_equal(
                    c[chosen] @ omega @ label % 2, epsilon[chosen]))
                assert np.array_equal(c @ omega @ correction % 2, epsilon)
                v = pauli(correction) @ unitary
                expected = projector([pauli(c[j]) for j in chosen], n)
                assert close(v @ projector(inputs, n) @ v.conj().T, expected)
                for error, error_matrix in zip(labels, matrices):
                    incoming = g0 @ omega @ error % 2
                    mapped = c @ omega @ (error @ q % 2) % 2
                    assert np.array_equal(mapped, k @ incoming % 2)
                    physical_error = v @ error_matrix @ v.conj().T
                    for row, bit in zip(c, mapped):
                        check = pauli(row)
                        assert close(check @ physical_error,
                                     (-1) ** int(bit) * physical_error @ check)
                    syndrome_cases += 1
                if words == ['XX', 'ZZ'] and signs == (1, 0):
                    assert close(transformed[0], tensor([Y, Y]))
                    assert combinations[tuple(label_of('XX'))][1] == 1
                    regression_seen = True
                cases += 1
                recovered_rows += len(pure)
    assert regression_seen
    result = {'executed_at_utc': utc_now(), 'scope': 'Finite Pauli-matrix checks of signed output and syndrome conversion',
              'cases': cases, 'output_rows_checked': recovered_rows,
              'error_syndrome_cases': syndrome_cases, 'projection_sign_regression_passed': regression_seen,
              'signed_products_match_matrices': True, 'coefficient_relations_passed': True,
              'redundant_sign_correction_passed': True, 'code_projectors_match': True,
              'syndrome_conversion_passed': True, 'all_passed': True,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    save_report(result, report_path)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
