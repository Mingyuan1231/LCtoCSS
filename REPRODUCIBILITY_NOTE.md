# Reproducibility of the quantum-code applications

This note documents the controlled BB recovery tests, published PBB and XYZ-IS
classifications, independent certificate verification, and recorded runtimes in
**Local Clifford equivalence to CSS codes: Polynomial-time recognition and
conversion**. Paths below are relative to this software directory.

## Controlled BB recovery and decoder checks

The known positive inputs use the `[[72,12,6]]`, `[[144,12,12]]`, and
`[[288,12,18]]` bivariate bicycle codes in Table 3 of
[Bravyi et al., arXiv:2308.07915v2](https://arxiv.org/html/2308.07915v2).
Their published distances are inherited by local Clifford conversion.
The recognizer and signed conversion are implemented in `lc_css.py`.
The construction and decoder tests are implemented in
[`application_demo/run_demo.py`](application_demo/run_demo.py), with the
original run recorded in
[`application_demo/results/summary.json`](application_demo/results/summary.json).

For each size, native generator pairs are mixed as
`(g_X,j, g_Z,j) -> (g_X,j*g_Z,j, g_Z,j)`, then each qubit receives an
independently sampled local Clifford from `I,H,S,HS,SH,HSH`. The five Python
random seeds are `20261009 + 1000*n + j`, for `j=0,...,4`. Recognition receives
only the binary checks and qubit count. Input signs are supplied to the signed
conversion step. The generating transformations are retained as verification
metadata.

All 15 inputs are accepted. Each has maximum input check weight 12 and maximum
qubit participation 9. Every recovered output has weight-six checks and
participation six. The numbers of X/Z checks are respectively 36/36, 72/72,
and 144/144; their ranks are 30/30, 66/66, and 138/138. Orthogonality, equality
of the transformed input and output row spaces, Pauli signs, and redundant
generator relations are checked separately. The three saved `bb_*_instance.json`
files contain complete signed inputs, local maps, check supports, Pauli
corrections, and syndrome maps.

The decoder is `ldpc.BpOsdDecoder`, using LDPC 2.4.1 and NumPy 2.3.5. Its
settings are 50 belief-propagation iterations, minimum-sum updates with scaling
0.625, parallel scheduling, one OpenMP thread, and order-zero `OSD_0`
postprocessing. Independent physical Pauli errors have
`(p_X,p_Y,p_Z)=(0.001,0.002,0.007)`. Exact syndromes and the transformed
site-dependent marginal probabilities are supplied to separate X/Z decoders.
This use of marginal probabilities leaves within-qubit X/Z correlation unused.

For the first deformation at each size, the tests cover every single-qubit
Pauli error, 1000 random errors with seed `913000+n`, and two deliberately
inserted zero-syndrome logical errors. Across the three sizes, 1512 single-qubit
errors and 3000 sampled errors are corrected; six logical-error controls remain
failures. All 4518 trials preserve the syndrome relations and the
stabilizer-membership test of logical success between the two descriptions.
These are finite decoder-interface checks with exact syndromes.

### Original BB timing and bounded-search comparison

The original run used Python 3.12.14 on Windows 10 with an Intel Core i7-10700KF.
Recognition timing includes factor decomposition, affine solves, local-map
construction, and sparse-check extraction. It excludes input generation,
signed transport, decoder trials, and file I/O. Full setup adds the syndrome
map, sign correction, and algebraic validation. Each median is over five
different deformations.

| Data qubits | Median recognition (s) | Median full setup (s) | Bounded search on the first deformation |
| --- | --- | --- | --- |
| 72 | 0.018 | 0.067 | Unresolved after 250000 nodes, 2.74 s |
| 144 | 0.113 | 0.324 | Unresolved after 250000 nodes, 5.65 s |
| 288 | 1.270 | 2.215 | Unresolved after 224256 nodes, 10.01 s |

The exhaustive baseline uses depth-first enumeration of six local binary
Clifford choices, with partial-rank pruning, and stops at a witness, complete
exhaustion, 250000 attempted assignments, or approximately 10 seconds. The
recorded cutoffs are unresolved outcomes. Recognition timings for the
same first-deformation inputs are retained in `summary.json`; a median over
different inputs is not a matched timing measurement for the search input.

See [`application_demo/README.md`](application_demo/README.md) for reproduction
commands and its output options. Reruns write to `outputs/bb/` and require fresh output
files; the distributed results are read-only.

## Source data and input conventions

The PBB inputs are all 368 entries of the public merged publication catalog in
the [qcode-discovery snapshot](https://github.com/qiskit-community/qcode-discovery/tree/db051c576f51c43d2bf9923b6a87033e7f88f0fe),
at commit `db051c576f51c43d2bf9923b6a87033e7f88f0fe`.
The source file is `campaign7_publication_merged.jsonl`, saved under
`research_applications/sources/qcode_discovery/results/`.

Writing $A_0,B_0,C_0,D_0$ for the four $\ell m\times\ell m$ binary circulant
blocks reconstructed from each entry's polynomial tuples, the input in grouped
$(x\mid z)$ coordinates is

$$
\widetilde B_{\mathrm{PBB}}=
\begin{pmatrix}
A_0&B_0&C_0&D_0\\
0&0&B_0^{\mathsf T}&A_0^{\mathsf T}
\end{pmatrix},\qquad n=2\ell m.
$$

Qubits are ordered by block and then lexicographically on
$\mathbb Z_\ell\times\mathbb Z_m$. A monomial $x^ay^b$ places a one in column
$(i+a,j+b)$ of row $(i,j)$, with addition modulo $\ell$ and $m$.
The recognizer uses interleaved $x,z$ columns. Commutativity and the catalog
value of $k=n-\operatorname{rank}B$ are checked for every input.

The XYZ-IS input is reconstructed from Eqs. (96)-(102) of
[Baldelli et al., arXiv:2607.14988v1](https://arxiv.org/html/2607.14988v1).
Qubits follow the nine tensor coordinates of the subset construction, with the
last coordinate varying fastest. The resulting 768 weight-eight checks have
rank 503 and form one indecomposable factor on 512 qubits. A separate
Kronecker-product construction reproduces the same binary matrix.

The source manifests record the versions and SHA-256 hashes:

- [`research_applications/source_manifest.json`](research_applications/source_manifest.json)
- [`research_applications/xyz_source_manifest.json`](research_applications/xyz_source_manifest.json)

Each PBB certificate is associated with its catalog line and code identifier;
the XYZ-IS record specifies its source equation and ordered subset tuples.
Hashes of the reconstructed binary matrices bind these descriptions to the
saved certificates. The complete correspondence is recorded in
[`research_applications/results/input_certificate_map.json`](research_applications/results/input_certificate_map.json).
[`CERTIFICATES.md`](research_applications/CERTIFICATES.md) specifies qubit and
row ordering, binary packing, and equation indices.

## Independent certificate verification

For a positive certificate, the verifier checks $\det Q_i=1$ for every local
map, orthogonality of the returned pure $X$ and pure $Z$ check spaces, and
equality of their combined span with the transformed input space $SQ$.

A negative certificate specifies a qubit set $V$, its factor basis $G_V$, and
one inconsistency witness for each anchor branch. The verifier checks that
$G_V$ spans the input's projection onto $V$, that its embedded rows belong to
$S$, and that its RREF support graph is connected. These conditions establish
an indecomposable direct factor. Writing either affine branch as $L\theta=f$,
its inconsistency witness is a binary row $\lambda$ satisfying

$$
\lambda L=0,\qquad \lambda f=1.
$$

The support of $\lambda$ is stored as the indices of equations whose sum is
$0=1$. The verifier reconstructs these equations from $G_V$ and its ordinary
annihilator, following the canonical ordering in
[`CERTIFICATES.md`](research_applications/CERTIFICATES.md), and checks both
witnesses. The verification entry point imports no recognition routine and
does not use the solver's elimination history. It independently checks all
11 positive witnesses, 714 PBB branch certificates, and two XYZ-IS branch
certificates against their inputs.

From this software directory, run:

```powershell
python -m pip install -r requirements.txt
python -X utf8 research_applications/verify_all_certificates.py
```

The verifier writes its report to
`outputs/verification/all_certificates_verification.json` and writes a new
input map in that output directory. The classification files and original recognition timings are
not overwritten. The expected report has `all_passed: true`,
`verified_inputs: 369`, `positive_witnesses: 11`, and
`negative_branch_certificates: 716`.

## Repeating the computations

```sh
python research_applications/screen_pbb_catalog.py
python research_applications/check_xyz_is.py
python research_applications/verify_all_certificates.py --results-dir outputs/published_codes --output-dir outputs/rerun_verification
python -m pip install -r application_demo/requirements.txt
python application_demo/run_demo.py
```

The classifier outputs go to `outputs/published_codes/`; the BB outputs go to
`outputs/bb/`. Select fresh run directories with `--output-dir PATH`. The
distributed `results/` files are never the default output targets.

The PBB source snapshot consists of its unmodified catalog and upstream license.
The XYZ construction is explicitly listed in `research_applications/sources/xyz_is.json`;
its manifest links to the versioned paper and identifies the equation range.
The independent verifier reconstructs the subset construction from its own
constants and compares them with this input file. Verification requires no
network connection or manuscript checkout.

## Runtime environment and timing scope

The original computations used an Intel Core i7-10700KF running Windows 10,
Python 3.12.14, and NumPy 2.3.5, in a single Python process.

| Computation | Recorded time | Included | Excluded |
| --- | --- | --- | --- |
| Serial PBB screening, all 368 entries | 49.7 s | Matrix construction, commutativity and rank checks, recognition, internal certificate checks, and input hashing | Catalog and result-file I/O; subsequent independent verification |
| XYZ-IS recognition | 0.663 s | Recognition and its internal checks | Input construction; independent verification; output serialization |

The PBB figure is the sum of the per-entry intervals. These are individual
measured runs. The exact values, environment, and timing boundaries are retained
in [`runtime_and_timing.json`](research_applications/runtime_and_timing.json).
The independent verification report records its own run separately.

Code hashes in recorded measurements identify the source used for those
measurements. `SHA256SUMS` identifies the distributed implementation. Repeated
runs record their own code hashes and elapsed times; they do not replace the
recorded measurements. Timing equality across machines or software versions is
not expected.
