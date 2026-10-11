"""Reproduce the LC-to-CSS bivariate-bicycle application.

Python 3.12; NumPy is needed for signs; LDPC is needed for the decoder experiment.
Binary rows use interleaved (x_0,z_0,...,x_{n-1},z_{n-1}) coordinates.
The recognizer receives only the deformed checks, never their construction.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import argparse
import hashlib
import importlib.metadata
import itertools
import json
import platform
import random
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lc_css import (
    GATES, GATE_NAMES, act, coordinates, inverse, join, omega, positions,
    recognize, rref, sign_tables, signed_product, split, transform_signs,
    verify_and_transport,
)
from tools.reproducibility import output_directory

HERE = Path(__file__).resolve().parent
CASES = [(6, 6, (1, 2), 6), (12, 6, (1, 2), 12),
         (12, 12, (2, 7), 18)]


def bicycle(l, m, ay):
    h = l*m
    a, b = [], []
    for i in range(l):
        for j in range(m):
            a.append((1 << (((i+3) % l)*m+j)) ^
                     (1 << (i*m+(j+ay[0]) % m)) ^
                     (1 << (i*m+(j+ay[1]) % m)))
            b.append((1 << (i*m+(j+3) % m)) ^
                     (1 << (((i+1) % l)*m+j)) ^
                     (1 << (((i+2) % l)*m+j)))
    def transpose(rows):
        return [sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(h)]
    hx = [u | (v << h) for u, v in zip(a, b)]
    hz = [u | (v << h) for u, v in zip(transpose(b), transpose(a))]
    n = 2*h
    return [join(v, 0, n) for v in hx] + [join(0, v, n) for v in hz]


def exhaustive(rows, n, node_limit, seconds):
    """Exact six-frame DFS with partial-rank pruning; cutoff means UNRESOLVED."""
    started = time.perf_counter()
    g, _ = rref(rows, 2*n)
    rank = len(g)
    columns = []
    for i in range(n):
        x = sum(((row >> (2*i)) & 1) << j for j, row in enumerate(g))
        z = sum(((row >> (2*i+1)) & 1) << j for j, row in enumerate(g))
        columns.append((x, z, x ^ z))
    choices = list(itertools.permutations(range(3), 2))
    counters = {"nodes": 0, "pruned": 0, "complete_assignments": 0}
    def extend(basis, vector):
        for pivot in basis:
            vector = min(vector, vector ^ pivot)
        if not vector:
            return basis
        return tuple(sorted((*basis, vector), reverse=True))
    class Cutoff(Exception):
        pass
    def search(i, xb, zb):
        if i == n:
            counters["complete_assignments"] += 1
            return True
        for a, b in choices:
            if counters["nodes"] >= node_limit:
                raise Cutoff("node_limit")
            if counters["nodes"] % 1024 == 0 and time.perf_counter()-started >= seconds:
                raise Cutoff("time_limit")
            counters["nodes"] += 1
            xx, zz = extend(xb, columns[i][a]), extend(zb, columns[i][b])
            if len(xx) + len(zz) > rank:
                counters["pruned"] += 1
            elif search(i+1, xx, zz):
                return True
        return False
    try:
        status = "CSS_found" if search(0, (), ()) else "not_LC_CSS"
    except Cutoff as exc:
        status = "unresolved_" + str(exc)
    return {**counters, "status": status, "seconds": time.perf_counter()-started,
            "node_limit": node_limit, "time_limit_seconds": seconds,
            "full_frame_space": str(6**n)}


def decode_trials(rows, qs, transport, n, seed, sample_count):
    import numpy as np
    from ldpc import BpOsdDecoder
    hx, hz = transport["hx"], transport["hz"]
    def array(checks):
        return np.array([[(row >> i) & 1 for i in range(n)] for row in checks], dtype=np.uint8)
    # Deliberately anisotropic, independent per-qubit physical Pauli noise.
    original_probs = [0.990, 0.001, 0.007, 0.002]  # I, X, Z, Y
    css_probs = []
    for q in qs:
        mapped = [0.0]*4
        for label, probability in enumerate(original_probs):
            mapped[act([label], [q])[0]] = probability
        css_probs.append(mapped)
    px = [p[1]+p[3] for p in css_probs]
    pz = [p[2]+p[3] for p in css_probs]
    options = dict(max_iter=50, bp_method="minimum_sum", ms_scaling_factor=0.625,
                   schedule="parallel", omp_thread_count=1, osd_method="OSD_0", osd_order=0)
    dx = BpOsdDecoder(array(hz), error_channel=px, **options)
    dz = BpOsdDecoder(array(hx), error_channel=pz, **options)
    bomega = [omega(row, n) for row in rows]
    comega = [omega(join(row, 0, n), n) for row in hx]
    comega += [omega(join(0, row, n), n) for row in hz]
    bg, bp = rref(rows, 2*n)
    cg, cp = rref(act(rows, qs), 2*n)
    inv = inverse(qs)
    def trial(error):
        css_error = act([error], qs)[0]
        sin = sum(((row & error).bit_count() & 1) << i for i, row in enumerate(bomega))
        adapted = [(combination & sin).bit_count() & 1 for combination in transport["syndrome_map"]]
        direct = [(row & css_error).bit_count() & 1 for row in comega]
        assert adapted == direct
        x = dx.decode(np.array(adapted[len(hx):], dtype=np.uint8))
        z = dz.decode(np.array(adapted[:len(hx)], dtype=np.uint8))
        estimate = join(sum(int(v) << i for i, v in enumerate(x)),
                        sum(int(v) << i for i, v in enumerate(z)), n)
        pulled = act([estimate], inv)[0]
        residual, css_residual = error ^ pulled, css_error ^ estimate
        assert all((row & residual).bit_count() % 2 == 0 for row in bomega)
        css_success = coordinates(css_residual, cg, cp) is not None
        input_success = coordinates(residual, bg, bp) is not None
        assert css_success == input_success
        return input_success
    single = [label << (2*i) for i in range(n) for label in (1, 2, 3)]
    rng = random.Random(seed)
    sampled = [sum(label << (2*i) for i, label in enumerate(
               rng.choices(range(4), weights=original_probs, k=n))) for _ in range(sample_count)]
    start = time.perf_counter()
    single_corrected = sum(trial(error) for error in single)
    sample_corrected = sum(trial(error) for error in sampled)
    # Zero-syndrome logical errors must fail the stabilizer-membership test.
    def logical_vector(checks, stabilizers):
        basis, pivots = rref(checks, n)
        sb, sp = rref(stabilizers, n)
        for free in range(n):
            if free in pivots:
                continue
            value = (1 << free) | sum(((row >> free) & 1) << p for row, p in zip(basis, pivots))
            if coordinates(value, sb, sp) is None:
                return value
        raise AssertionError("Expected an encoded qubit")
    logical_x = join(logical_vector(hz, hx), 0, n)
    logical_z = join(0, logical_vector(hx, hz), n)
    for logical in (logical_x, logical_z):
        assert not trial(act([logical], inv)[0])
    return {"decoder": "ldpc.BpOsdDecoder", "ldpc_version": importlib.metadata.version("ldpc"),
            "options": options, "physical_probabilities_IXZY": original_probs,
            "css_x_marginals": px, "css_z_marginals": pz,
            "single_qubit_trials": len(single), "single_qubit_corrected": single_corrected,
            "random_trials": sample_count, "random_corrected": sample_corrected,
            "logical_negative_controls": 2,
            "syndrome_and_logical_transport_passed": len(single)+sample_count+2,
            "random_seed": seed, "seconds_including_checks": time.perf_counter()-start}


def small_crosscheck():
    """Compare recognition with independent row equations and exact local search."""
    import sys
    sys.path.insert(0, str(HERE.parent / "validation"))
    import check_algorithms as reference
    import validation_reference as ref
    tested = 0
    for n in range(1, 6):
        edges = list(itertools.combinations(range(n), 2))
        rng = random.Random(600+n)
        masks = range(1 << len(edges)) if n < 5 else [rng.randrange(1 << len(edges)) for _ in range(40)]
        for mask in masks:
            rows = ref.graph_rows(n, [edge for j, edge in enumerate(edges) if mask >> j & 1])
            result = recognize(rows, n)
            expected = reference.row_decide(rows, n)
            oracle = ref.exhaustive_local_oracle(rows, n)
            assert (result is not None) == (expected is not None) == oracle["equivalent"]
            bounded = exhaustive(rows, n, 100000, 10)
            assert bounded["status"] == ("CSS_found" if oracle["equivalent"] else "not_LC_CSS")
            if result is not None:
                assert rref(act(rows, result["q"]), 2*n)[0] == rref(result["pure_rows"], 2*n)[0]
            tested += 1
    return tested


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--nodes", type=int, default=250000)
    parser.add_argument("--search-seconds", type=float, default=10.0)
    parser.add_argument("--sizes", type=int, nargs="*", default=[72, 144, 288])
    parser.add_argument("--output-dir", type=Path, default=HERE.parent / "outputs/bb")
    args = parser.parse_args()
    if args.seeds < 1 or args.samples < 0 or args.nodes < 1 or args.search_seconds <= 0:
        parser.error("seeds and nodes must be positive, samples nonnegative, and search-seconds positive")
    if not args.sizes or any(n not in (72, 144, 288) for n in args.sizes):
        parser.error("sizes must be selected from 72, 144, and 288")
    out = output_directory(args.output_dir)
    if any((out / name).exists() for name in ["summary.json", *[f"bb_{n}_instance.json" for n in args.sizes]]):
        parser.error("Output files already exist; select a fresh --output-dir")
    table = sign_tables()
    small_count = small_crosscheck()
    print(f"Independent small-instance crosschecks: {small_count}", flush=True)
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "python": platform.python_version(), "platform": platform.platform(),
               "processor": platform.processor(), "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
               "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "core_sha256": hashlib.sha256((HERE.parent / "lc_css" / "__init__.py").read_bytes()).hexdigest(),
               "small_crosscheck_cases": small_count, "cases": []}
    for l, m, ay, distance in CASES:
        n = 2*l*m
        if n not in args.sizes:
            continue
        native = bicycle(l, m, ay)
        half = n // 2
        # Hide the original X/Z generator partition as well as the local frame.
        # (g_X,j, g_Z,j) -> (g_X,j g_Z,j, g_Z,j) is invertible even for redundant checks.
        mixed = [signed_product(native, [0]*n, (1 << j) | (1 << (half+j)), n)
                 for j in range(half)]
        mixed_rows = [value for value, sign in mixed] + native[half:]
        mixed_signs = [sign for value, sign in mixed] + [0]*half
        runs, first = [], None
        for instance in range(args.seeds):
            seed = 20261009 + 1000*n + instance
            rng = random.Random(seed)
            names = [rng.choice(list(GATES)) for _ in range(n)]
            hidden = [GATES[name] for name in names]
            rows = act(mixed_rows, hidden)
            signs = transform_signs(mixed_rows, mixed_signs, names, table)
            start = time.perf_counter()
            result = recognize(rows, n)
            elapsed = time.perf_counter()-start
            assert result is not None
            transport = verify_and_transport(rows, signs, result, n, table)
            full_elapsed = time.perf_counter()-start
            pure_initial = sum(not (split(row, n)[0] and split(row, n)[1]) for row in rows)
            run = {"seed": seed, "seconds_recognition_and_sparse_output": elapsed,
                   "seconds_including_signed_transport_and_validation": full_elapsed,
                   "input_pure_checks": pure_initial, "rank": result["rank"],
                   "input_max_weight": max(sum(bool(row >> (2*i) & 3) for i in range(n)) for row in rows),
                   "input_max_degree": max(sum(bool(row >> (2*i) & 3) for row in rows) for i in range(n)),
                   "x_checks": len(transport["hx"]), "z_checks": len(transport["hz"]),
                   "max_weight": transport["max_weight"], "max_degree": transport["max_degree"],
                   "factors": result["factors"]}
            runs.append(run)
            print(f"n={n}, seed={seed}: {elapsed:.4f}s, rank={result['rank']}, "
                  f"checks={run['x_checks']}+{run['z_checks']}, w={run['max_weight']}, "
                  f"degree={run['max_degree']}", flush=True)
            if first is None:
                first = (rows, signs, result, transport, names)
        rows, signs, result, transport, names = first
        search = exhaustive(rows, n, args.nodes, args.search_seconds)
        print(f"n={n} search: {search}", flush=True)
        decoding = decode_trials(rows, result["q"], transport, n, 913000+n, args.samples)
        print(f"n={n} decoder: {decoding['single_qubit_corrected']}/{3*n} single; "
              f"{decoding['random_corrected']}/{args.samples} random", flush=True)
        instance_data = {"n": n, "label_order": "interleaved x,z; zero-based qubits",
                         "input_generator_mixing": "(g_X,j,g_Z,j) -> (g_X,j*g_Z,j,g_Z,j)",
                         "input_rows_hex": [hex(row) for row in rows], "input_signs": signs,
                         "generation_gates_withheld_from_recognizer": names,
                         "recovery_gates": [GATE_NAMES[q] for q in result["q"]],
                         "recovery_row_matrices": result["q"],
                         "css_x_supports": [list(positions(row)) for row in transport["hx"]],
                         "css_z_supports": [list(positions(row)) for row in transport["hz"]],
                         "syndrome_map_input_row_indices": [list(positions(row)) for row in transport["syndrome_map"]],
                         "pauli_correction_hex": hex(transport["pauli_correction"]),
                         "css_signs_before_pauli_correction": transport["precorrection_signs"],
                         "decoder": decoding}
        artifact = out / f"bb_{n}_instance.json"
        artifact.write_text(json.dumps(instance_data, indent=2)+"\n")
        case = {"n": n, "k": n-result["rank"], "published_distance": distance,
                "l": l, "m": m, "A_y_exponents": ay,
                "runs": runs, "median_seconds": statistics.median(run["seconds_recognition_and_sparse_output"] for run in runs),
                "median_seconds_with_signed_transport_and_validation": statistics.median(run["seconds_including_signed_transport_and_validation"] for run in runs),
                "minimum_seconds": min(run["seconds_recognition_and_sparse_output"] for run in runs),
                "maximum_seconds": max(run["seconds_recognition_and_sparse_output"] for run in runs),
                "search_first_seed": search, "decoder_first_seed": decoding,
                "instance_file": artifact.name, "instance_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        summary["cases"].append(case)
        (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(f"Saved {out / 'summary.json'} and complete witnesses.", flush=True)


if __name__ == "__main__":
    main()
