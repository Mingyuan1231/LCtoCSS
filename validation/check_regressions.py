"""Check production boundary cases, physical signs, and certificate rejection."""
if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")

from copy import deepcopy
from itertools import product
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'research_applications'))
import lc_css as core
import screen_pbb_catalog as classifier
import verify_pbb_certificates as verifier
import check_signed_conversion as matrix
from tools.reproducibility import report_output, save_report


def rejects(function, *args):
    try:
        function(*args)
    except (ValueError, AssertionError):
        return
    raise AssertionError('Invalid input or certificate was accepted')


def main():
    output = report_output('regressions.json')
    table = core.sign_tables()
    for n in (0, 1, 5):
        for rows in ([], [0, 0]):
            answer = core.recognize(rows, n)
            assert answer is not None and answer['rank'] == 0 and answer['pure_rows'] == []
            converted = core.verify_and_transport(rows, [0]*len(rows), answer, n, table)
            assert converted['max_weight'] == converted['max_degree'] == 0
            classified = classifier.classify(rows, n)
            assert classified['lc_css'] and classified['max_check_weight'] == 0
            classified['n'] = n
            verifier.verify_positive(rows, classified)
    malformed = [([-1], 1), ([4], 1), ([1, 2], 1), ([True], 1), ([1.0], 1), ([], -1), ([], True)]
    for rows, n in malformed:
        rejects(core.recognize, rows, n)
        rejects(classifier.classify, rows, n)
    answer = core.recognize([1], 1)
    rejects(core.verify_and_transport, [1], [], answer, 1, table)
    rejects(core.verify_and_transport, [1], [2], answer, 1, table)
    rejects(core.verify_and_transport, [1, 1], [0, 1], core.recognize([1, 1], 1), 1, table)

    signed_cases = 0
    gate_names = list(core.GATES)
    for n, words in [(1, []), (2, ['XI']), (2, ['XX', 'ZZ']),
                     (3, ['ZZI', 'IZZ']), (3, ['XXX', 'ZZI', 'IZZ'])]:
        native = [sum(int(bit) << j for j, bit in enumerate(matrix.label_of(word))) for word in words]
        for shift in range(6):
            names = [gate_names[(shift+i) % 6] for i in range(n)]
            unitary = matrix.tensor(matrix.GATES[(shift+i) % 6] for i in range(n))
            for signs in product((0, 1), repeat=len(native)):
                rows = core.act(native, [core.GATES[name] for name in names])
                input_signs = core.transform_signs(native, signs, names, table)
                incoming = [(-1)**sign * matrix.pauli([(row >> j) & 1 for j in range(2*n)])
                            for row, sign in zip(rows, input_signs)]
                native_operators = [(-1)**sign * matrix.pauli(matrix.label_of(word))
                                    for word, sign in zip(words, signs)]
                assert all(matrix.close(u, unitary @ v @ unitary.conj().T)
                           for u, v in zip(incoming, native_operators))
                answer = core.recognize(rows, n)
                assert answer is not None
                converted = core.verify_and_transport(rows, input_signs, answer, n, table)
                recovery = matrix.tensor(matrix.GATES[gate_names.index(core.GATE_NAMES[q])] for q in answer['q'])
                correction = matrix.pauli([(converted['pauli_correction'] >> j) & 1 for j in range(2*n)])
                physical = correction @ recovery
                expected = matrix.projector([matrix.pauli([(row >> j) & 1 for j in range(2*n)])
                                             for row in answer['pure_rows']], n)
                assert matrix.close(physical @ matrix.projector(incoming, n) @ physical.conj().T, expected)
                signed_cases += 1

    records = [json.loads(line) for line in (ROOT / 'research_applications/results/pbb_classification.jsonl').read_text().splitlines()]
    catalog = [json.loads(line) for line in (ROOT / 'research_applications/sources/qcode_discovery/results/campaign7_publication_merged.jsonl').read_text().splitlines()]
    negative = next(r for r in records if not r['lc_css'])
    positive = next(r for r in records if r['lc_css'])
    bad_cases = []
    for mutation in ('missing_branch', 'negative_index', 'out_of_range', 'duplicate_index', 'empty_factor', 'wrong_factor', 'changed_sum'):
        changed = deepcopy(negative)
        if mutation == 'missing_branch': changed['negative_certificates'].pop()
        elif mutation == 'negative_index': changed['negative_certificates'][0]['equation_indices'] = [-1]
        elif mutation == 'out_of_range': changed['negative_certificates'][0]['equation_indices'] = [10**9]
        elif mutation == 'duplicate_index': changed['negative_certificates'][0]['equation_indices'] *= 2
        elif mutation == 'empty_factor': changed['rejected_factor'] = []
        elif mutation == 'wrong_factor': changed['local_rref_hex'][0] = '0x0'
        else: changed['negative_certificates'][0]['equation_indices'] = [0]
        rows = verifier.reconstruct(catalog[negative['catalog_line']-1])
        rejects(verifier.verify_negative, rows, changed)
        bad_cases.append(mutation)
    for mutation in ('singular_map', 'extra_map_entry', 'missing_checks', 'wrong_degree'):
        changed = deepcopy(positive)
        if mutation == 'singular_map': changed['q'][0] = [0, 0, 0, 0]
        elif mutation == 'extra_map_entry': changed['q'][0].append(0)
        elif mutation == 'missing_checks': changed['css_x_supports'] = []; changed['css_z_supports'] = []
        else: changed['max_qubit_participation'] += 1
        rows = verifier.reconstruct(catalog[positive['catalog_line']-1])
        rejects(verifier.verify_positive, rows, changed)
        bad_cases.append(mutation)
    report = {'all_passed': True, 'zero_space_cases': 6, 'invalid_binary_inputs': len(malformed),
              'signed_production_cases': signed_cases, 'corrupted_certificates_rejected': bad_cases}
    save_report(report, output)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
