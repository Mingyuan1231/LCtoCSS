# Bivariate bicycle CSS recovery

The controlled examples use the `[[72,12,6]]`, `[[144,12,12]]`, and
`[[288,12,18]]` BB codes. Five local-Clifford deformations per size test
recognition, sparse signed conversion, and a CSS decoder interface.

From the repository root with Python 3.12:

```sh
python -m pip install -r application_demo/requirements.txt
python application_demo/run_demo.py
```

Results go to `outputs/bb/`; the output files must not already exist. Select
another run with `--output-dir PATH`. `--sizes`, `--seeds`, `--samples`, `--nodes`,
and `--search-seconds` control the experiment. Use the defaults for the recorded
coverage.

The recorded experiment is in `results/summary.json`; complete signed inputs
and witnesses are in `results/bb_*_instance.json`. Check these witnesses without
running LDPC:

```sh
python validation/check_bb_witnesses.py
```

## Construction and conversion

Each native pair is mixed as `(g_X,g_Z) -> (g_X*g_Z,g_Z)`, then a seeded random
Clifford is applied to each qubit. Recognition receives only binary rows and
the qubit count. Signs are supplied separately to signed conversion. Generating
gates are saved as input provenance and are not given to the recognizer.

`lc_css.recognize(rows,n)` returns sparse pure checks and local binary maps.
`lc_css.verify_and_transport(rows,signs,result,n,lc_css.sign_tables())` returns
the pure checks, their signs before correction, a Pauli correction, and generator
coefficients for syndrome conversion. The syndrome map is expressed against
the supplied redundant input rows; saved `syndrome_map_input_row_indices`
lists the selected coefficients per output check. Signs must describe a valid
commuting stabilizer group without minus identity.

## Recorded observations

- All 15 inputs are accepted. Input maximum weight and participation are 12 and
  9; output checks have weight six and participation six.
- Decoding uses independent biased Pauli noise and exact syndromes. Logical
  success is tested by stabilizer membership after recovery.
- The 4,518 interface trials comprise 1,512 single-qubit errors, 3,000 random
  errors, and six deliberate logical errors. The logical controls remain
  failures; every trial preserves syndrome and success relations.
- Bounded local-frame search uses rank pruning and stops at a witness, complete
  exhaustion, 250,000 nodes, or approximately ten seconds. Recorded cutoffs are
  unresolved outcomes.
- Recognition timing excludes signed transport and decoder trials. Full setup
  adds sign recovery and algebraic validation. Medians cover different inputs.

Parameters and distances come from Bravyi et al., Nature 627, 778 (2024),
[arXiv:2308.07915v2](https://arxiv.org/html/2308.07915v2), Table 3. Distances and
circuit-level thresholds are not measured here. The
[reproducibility note](../REPRODUCIBILITY_NOTE.md) gives seeds, noise settings,
and timing scope.
