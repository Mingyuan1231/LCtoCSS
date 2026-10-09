# Local Clifford recognition of CSS codes

Code and certificates for **Local Clifford equivalence to CSS codes:
Polynomial-time recognition and conversion**.

Developed by Mingyuan Wang.

The algorithm determines whether a binary stabilizer space can be converted to
Calderbank-Shor-Steane (CSS) form using single-qubit Clifford transformations.
Accepted inputs receive local maps and pure X/Z checks. The published-code
classification also supplies independently verifiable certificates for rejected
inputs.

## Quick start

Use Python 3.12. From the repository root:

```sh
python -m pip install -r requirements.txt
python check_release.py
python research_applications/verify_all_certificates.py
```

Expected: `all_passed: true`, `verified_inputs: 369`, `positive_witnesses: 11`,
and `negative_branch_certificates: 716`. The verifier reconstructs every input
and checks certificates without importing or running the recognizer. Reports
go to `outputs/verification/`.

## Results and reproduction

| Result | Input | Generator | Independent check | Expected output |
| --- | --- | --- | --- | --- |
| BB recovery and decoder interface | 72-, 144-, and 288-qubit BB codes, five deformations each | `application_demo/run_demo.py` | `validation/check_bb_witnesses.py`; signed and syndrome checks in the demo | 15 accepted inputs; weight-six checks; 4,518 interface trials |
| PBB classification | Pinned catalog of 368 polynomial tuples | `research_applications/screen_pbb_catalog.py` | `research_applications/verify_all_certificates.py` | 11 CSS-equivalent entries; 357 exclusions |
| XYZ-IS classification | Published 512-qubit subset construction | `research_applications/check_xyz_is.py` | `research_applications/verify_all_certificates.py` | Rank 503; two branch contradictions |
| Recognition and graph conversion | 1,136 stabilizer fixtures and all graphs through five vertices | `validation/validation_reference.py`; `validation/check_algorithms.py` | Exact six-choice local-Clifford search | Zero mismatches; valid CSS and bipartite witnesses |

The PBB result classifies CSS equivalence per input. It does not classify
pairwise LC equivalence among excluded entries. Published distances are source
metadata and are not recomputed. Finite validation supports the implementation;
the general correctness and complexity results are proved in the paper.

## Use the recognizer

Run Python from this directory, or add it to your Python module search path:

```python
from lc_css import recognize

# Bell stabilizers XX and ZZ, with interleaved (x0,z0,x1,z1) bits.
result = recognize([0b0101, 0b1010], n=2)
if result is not None:
    print(result["q"])          # One right-acting 2x2 map per qubit.
    print(result["pure_rows"])  # Packed X checks, then Z checks.
```

Rows are nonnegative Python integers; bit `2*i` is `x_i` and bit `2*i+1` is
`z_i`. Indices are zero-based. `q[i] = (a,b,c,d)` sends the row label `(x,z)` to
`(x*a+z*c, x*b+z*d)` over GF(2). `None` means no CSS conversion exists. Invalid
widths or noncommuting labels raise `ValueError`.

Redundant and zero rows, unused qubits, and the zero-qubit space are supported.
The sparse checks can be redundant; `x_count` separates the X and Z parts and
`rank` is the stabilizer rank. Binary recognition does not specify code-space
signs. For signed conversion, supply a consistent signed stabilizer group to
`verify_and_transport`; see the [BB example](application_demo/README.md).

Recognition uses the standard library. NumPy is needed for signed matrix checks
and independent input reconstruction; the optional decoder also needs LDPC.
Programs using mathematical assertions reject `python -O` and `python -OO`.

## Repository guide

- [Reproducibility](REPRODUCIBILITY_NOTE.md): constructions, seeds, timings, and commands.
- [Published-code applications](research_applications/README.md): inputs and reruns.
- [Certificate format](research_applications/CERTIFICATES.md): input identities and equation indices.
- [Validation](validation/README.md): algorithm, graph, phase, sparsity, and rejection checks.
- [Third-party notices](THIRD_PARTY_NOTICES.md): catalog provenance and license.

Application `results/` directories contain the recorded evidence. Reruns write
to ignored `outputs/` directories or an explicit external destination.
Classification and BB reruns require fresh output files. `SHA256SUMS` identifies
the distributed files; `python check_release.py` checks their integrity. The
software works independently of a manuscript checkout.

Original code and documentation use the [MIT License](LICENSE). The bundled
upstream PBB catalog retains its Apache-2.0 license.
