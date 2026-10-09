"""Independently verify the saved signed BB conversions without running LDPC."""
if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")

import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'research_applications'))
from reproducibility import report_output, save_report
from verify_pbb_certificates import verify_positive


def tables():
    identity = np.eye(2, dtype=complex)
    h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    s = np.diag([1, 1j])
    p = [identity, np.array([[0, 1], [1, 0]]), np.diag([1, -1]),
         np.array([[0, -1j], [1j, 0]])]
    gates = {'I': identity, 'H': h, 'S': s, 'HS': h@s, 'SH': s@h, 'HSH': h@s@h}
    conjugation, multiplication = {}, {}
    for name, unitary in gates.items():
        conjugation[name] = []
        for pauli in p:
            image = unitary @ pauli @ unitary.conj().T
            matches = [(i, sign) for i, target in enumerate(p) for sign in (0, 1)
                       if np.allclose(image, (-1)**sign * target)]
            assert len(matches) == 1
            conjugation[name].append(matches[0])
    for a in range(4):
        for b in range(4):
            multiplication[a, b] = next(k for k in range(4)
                                        if np.allclose(p[a] @ p[b], 1j**k * p[a ^ b]))
    return conjugation, multiplication


def multiply(a, b, n, multiplication):
    phase = sum(multiplication[(a >> (2*i)) & 3, (b >> (2*i)) & 3] for i in range(n)) % 4
    return a ^ b, phase


def transform(row, sign, names, conjugation):
    out = 0
    for i, name in enumerate(names):
        label, change = conjugation[name][(row >> (2*i)) & 3]
        out |= label << (2*i)
        sign ^= change
    return out, sign


def symplectic(a, b, n):
    return sum(((a >> (2*i)) & 1) * ((b >> (2*i+1)) & 1)
               + ((a >> (2*i+1)) & 1) * ((b >> (2*i)) & 1) for i in range(n)) % 2


def native_rows(case):
    ell, m = case['l'], case['m']
    size = ell*m
    a, b = np.zeros((size, size), dtype=np.uint8), np.zeros((size, size), dtype=np.uint8)
    for i in range(ell):
        for j in range(m):
            for di, dj in [(3, 0), (0, case['A_y_exponents'][0]), (0, case['A_y_exponents'][1])]:
                a[i*m+j, ((i+di) % ell)*m+(j+dj) % m] ^= 1
            for di, dj in [(0, 3), (1, 0), (2, 0)]:
                b[i*m+j, ((i+di) % ell)*m+(j+dj) % m] ^= 1
    hx, hz = np.hstack([a, b]), np.hstack([b.T, a.T])
    assert not np.any((hx @ hz.T) % 2)
    return [sum(1 << (2*int(i)) for i in np.flatnonzero(row)) for row in hx] + [
        sum(1 << (2*int(i)+1) for i in np.flatnonzero(row)) for row in hz]


def main():
    output = report_output('bb_witnesses.json')
    directory = ROOT / 'application_demo/results'
    summary = json.loads((directory / 'summary.json').read_text(encoding='utf-8'))
    conjugation, multiplication = tables()
    checked_rows = trials = 0
    details = []
    for case in summary['cases']:
        path = directory / case['instance_file']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == case['instance_sha256']
        record = json.loads(path.read_text(encoding='utf-8'))
        n = record['n']
        assert n == case['n'] and case['k'] == 12 and len(case['runs']) == 5
        native = native_rows(case)
        mixed = []
        for j in range(n//2):
            row, phase = multiply(native[j], native[n//2+j], n, multiplication)
            assert phase % 2 == 0
            mixed.append((row, phase//2))
        mixed.extend((row, 0) for row in native[n//2:])
        generated = [transform(row, sign, record['generation_gates_withheld_from_recognizer'], conjugation)
                     for row, sign in mixed]
        rows = [int(row, 16) for row in record['input_rows_hex']]
        assert generated == list(zip(rows, record['input_signs']))
        first = case['runs'][0]
        witness = {'n': n, 'q': record['recovery_row_matrices'],
                   'css_x_supports': record['css_x_supports'], 'css_z_supports': record['css_z_supports'],
                   'max_check_weight': first['max_weight'], 'max_qubit_participation': first['max_degree']}
        verify_positive(rows, witness)
        pure = [sum(1 << (2*i) for i in support) for support in record['css_x_supports']]
        pure += [sum(1 << (2*i+1) for i in support) for support in record['css_z_supports']]
        for name, q in zip(record['recovery_gates'], record['recovery_row_matrices']):
            x, z = conjugation[name][1][0], conjugation[name][2][0]
            assert q == [x & 1, x >> 1, z & 1, z >> 1]
        transformed = [transform(row, sign, record['recovery_gates'], conjugation)
                       for row, sign in zip(rows, record['input_signs'])]
        correction = int(record['pauli_correction_hex'], 16)
        indices = record['syndrome_map_input_row_indices']
        signs = record['css_signs_before_pauli_correction']
        assert len(pure) == len(indices) == len(signs)
        for target, selected, expected_sign in zip(pure, indices, signs):
            assert len(selected) == len(set(selected)) and all(0 <= j < len(rows) for j in selected)
            value = exponent = 0
            for j in selected:
                value, phase = multiply(value, transformed[j][0], n, multiplication)
                exponent += phase + 2*transformed[j][1]
            assert value == target and exponent % 4 == 2*expected_sign
            assert symplectic(target, correction, n) == expected_sign
            checked_rows += 1
        assert case['decoder_first_seed'] == record['decoder']
        decoder = record['decoder']
        assert decoder['single_qubit_trials'] == decoder['single_qubit_corrected'] == 3*n
        assert decoder['random_trials'] == decoder['random_corrected'] == 1000
        assert decoder['logical_negative_controls'] == 2
        assert decoder['syndrome_and_logical_transport_passed'] == 3*n + 1002
        trials += decoder['syndrome_and_logical_transport_passed']
        details.append({'n': n, 'signed_output_rows': len(pure), 'max_weight': first['max_weight']})
    assert trials == 4518 and len(details) == 3
    report = {'all_passed': True, 'saved_instances': 3, 'signed_output_rows': checked_rows,
              'recorded_decoder_trials': trials, 'decoder_rerun': False,
              'syndrome_map_identity_checked': True, 'cases': details}
    save_report(report, output)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
