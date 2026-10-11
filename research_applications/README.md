# Classification of published quantum codes

The PBB inputs are all 368 entries at qiskit-community/qcode-discovery commit
`db051c576f51c43d2bf9923b6a87033e7f88f0fe`. Eleven are LC-equivalent to CSS,
matching the previously known positives; the remaining 357 are excluded.
The XYZ-IS input is the 512-qubit, nine-logical-qubit code in Baldelli et al.,
[arXiv:2607.14988v1](https://arxiv.org/html/2607.14988v1), Eq. (102). Its 768
weight-eight checks have rank 503, and it is not LC-equivalent to CSS.

Inputs are reconstructed from published descriptions without additional
deformations. Classification is per input and does not determine pairwise
equivalence among rejected codes. Distances are source metadata, not recomputed
by these programs.

## Verify distributed certificates

From the repository root with Python 3.12 and NumPy:

```sh
python research_applications/verify_all_certificates.py
```

Expected: 369 inputs, 11 positive witnesses, 358 negative inputs, and 716 negative
branch certificates, with `all_passed: true`. Reports go to
`outputs/verification/`. The verifier imports no recognition code and preserves
the distributed classifications and timings. `--output-dir PATH` selects a
report directory; `--results-dir PATH` selects the certificates to verify.

## Repeat recognition

```sh
python research_applications/screen_pbb_catalog.py
python research_applications/check_xyz_is.py
python research_applications/verify_all_certificates.py --results-dir outputs/published_codes --output-dir outputs/rerun_verification
```

Generators write to `outputs/published_codes/`. Each requires its output file
to be absent. For another run, pass the same fresh `--output-dir PATH` to both
generators, then use it as `--results-dir` for verification. Recognition always
runs for every input; saved classifications are not silently reused.

## Files

| File | Purpose |
| --- | --- |
| `screen_pbb_catalog.py` | PBB construction and certificate-producing recognition using `lc_css` |
| `check_xyz_is.py` | Subset and Kronecker-product construction, recognition, and verification |
| `verify_pbb_certificates.py` | Independent matrix, factor, annihilator, and witness checks; standalone PBB verification |
| `verify_all_certificates.py` | Independent verification of PBB and XYZ-IS certificates |
| `results/pbb_classification.jsonl` | 368 decisions and certificates |
| `results/pbb_summary.json` | Recorded counts and screening times |
| `results/xyz_is_512.json` | XYZ input and two negative certificates |
| `results/input_certificate_map.json` | Source identities and hashes for 369 inputs |
| `runtime_and_timing.json` | Recorded environment and timing boundaries |
| `source_manifest.json` | Pinned catalog and license URLs, sizes, and hashes |
| `xyz_source_manifest.json` | Versioned source reference and construction hash |
| `sources/` | Licensed PBB catalog and explicit XYZ subset tuples |

[CERTIFICATES.md](CERTIFICATES.md) specifies coordinates and equation indices.
[Reproducibility](../docs/REPRODUCIBILITY.md) describes constructions and
timings. The certificates verify finite algebraic identities; the paper's
two-branch theorem supplies their classification meaning.
