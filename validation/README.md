# Implementation validation

Run from the repository root with Python 3.12. NumPy is needed for signed,
phase, BB, and certificate checks.

```sh
python validation/check_algorithms.py
python validation/check_sparse_output.py
python validation/check_signed_conversion.py
python validation/check_bb_witnesses.py
python validation/check_regressions.py
python validation/check_anchors.py
python validation/check_phases.py
```

Each command writes a JSON report under `outputs/validation/` and exits nonzero
on failure. `--output PATH` selects another report. Output paths cannot replace
source files or recorded results. Programs with mathematical assertions reject
optimized Python modes, which would disable those checks.

## Coverage

| Check | Coverage |
| --- | --- |
| `check_algorithms.py` | Production and independent row-action recognition on 1,136 fixtures; exact local-search comparison and bipartite construction for all 1,100 graphs through five vertices; Clifford table and signed triangle-to-star map |
| `check_sparse_output.py` | Sparse projections on 1,034 accepted fixtures; spans, supports, check counts, weights, and incidences |
| `check_signed_conversion.py` | 114 signed inputs, 660 output rows, and 5,208 syndrome cases using Pauli matrices |
| `check_bb_witnesses.py` | Three saved BB conversions: local maps, signs, corrections, generator relations, sparse checks, and recorded decoder counts |
| `check_regressions.py` | Zero-rank and zero-qubit inputs, malformed inputs, signed production conversion, and corrupted certificate rejection |
| `check_anchors.py` | Direct enumeration showing why both branches and factorwise anchors are necessary |
| `check_phases.py` | All 24 signed single-qubit Clifford maps; 550 binary spaces and 2,566 sign assignments through three qubits; sampled complex operators through four qubits |

The fixture `cases.json` contains exact-search answers for 550 isotropic spaces
on zero through three qubits, 576 seeded redundant inputs on four through six
qubits, and ten structured examples. Its seed is 20261003, with 32 samples per
rank. Regenerate cases without changing the fixture:

```sh
python validation/validation_reference.py
```

The output is `outputs/validation/reference.json`. Its `cases` field must equal
the fixture's `cases`; timings, timestamps, and code hashes are run-specific.

## Extended exact checks

```sh
python validation/check_exhaustive.py
```

This independent implementation enumerates all isotropic spaces through four
qubits, adds seeded examples through eight qubits, checks generator changes and
alternate anchors, and compares against two exact local-map oracles. The full
suite has 20,991 cases and can take several minutes.

The optional C++17 suite compares the adjacency criterion with complete labeled
local-complementation orbits:

```sh
g++ -std=c++17 -O2 validation/graph_exhaustive.cpp -o graph_exhaustive
./graph_exhaustive 7 outputs/validation/graph_exhaustive.json
```

On Windows, use `-o graph_exhaustive.exe` and `./graph_exhaustive.exe`. The
executable takes the maximum order and an output path. The complete suite
covers 2,131,020 graphs on zero through seven vertices. Finite enumeration is an
implementation check, not an asymptotic complexity measurement.

## Matrix conventions

Production and `check_algorithms.py` use right row action. The independent
`validation_reference.py` uses column action: `(a,b,c,d)` sends `(x,z)^T` to
`(a*x+b*z,c*x+d*z)^T`. A column map `T` corresponds to the production map
`Q = T^T`. The reference projector swaps `b` and `c`; branch order may select
another valid witness. `check_algorithms.py` verifies this transpose rule.

Exact-search oracles do not use the two-branch criterion. Certificate
verification reconstructs ordinary-annihilator equations without importing the
recognizer. These are checks of finite identities; the paper supplies their
general mathematical justification.
