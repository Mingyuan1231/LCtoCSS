# Certificate format and correspondence to the input

All indices in the JSON files and Python programs are zero-based unless stated
otherwise. The manuscript's mathematical indices are one-based.

## Source identities

The PBB catalog is the exact file
`results/campaign7_publication_merged.jsonl` at repository commit
`db051c576f51c43d2bf9923b6a87033e7f88f0fe` of
<https://github.com/qiskit-community/qcode-discovery>.
Its SHA-256 is
`c0197f0cb5e03bfc94215555ed0fbfefaa72d1c0e3fb89a31fc5045c5f7e25ba`.
`source_manifest.json` records the pinned URLs and hashes of the source snapshot.
In every classification record, `catalog_line` is the **one-based** line in that
file, and `code_id` is the source identifier, or `catalog_line_NNN` when absent.
No claim of pairwise LC inequivalence between these catalog records is made.

The XYZ-IS input is exactly the nine-coordinate subset construction specified in
Eqs. (96)--(102) of Baldelli et al., arXiv:2607.14988v1. It is the example to which
Definition 12 (genuine under local Pauli relabeling) is applied in the manuscript.
`xyz_source_manifest.json` identifies the versioned paper and the saved subset construction. The three
ordered lists of subsets are stored in `results/xyz_is_512.json`; the independent
verifier also specifies them directly, rather than trusting a solver-selected
list. Source distance values are not recertified by this calculation.

## Qubits, rows, and binary encoding

PBB qubits are ordered first by block `b=0,1`, then by `(i,j)` with
`0 <= i < ell`, `0 <= j < m`. Their index is `b*ell*m + i*m + j`.
A polynomial term `(a,b)` contributes a one in row `(i,j)`, column
`((i+a) mod ell,(j+b) mod m)` of the corresponding circulant block.
For the published grouped-coordinate matrix

```text
        [ A  B | C    D   ]
H_PBB = [ 0  0 | B^T  A^T ],
```

the upper row block precedes the lower row block; rows inside each block follow
the same lexicographic `(i,j)` order. All additions are over GF(2).

XYZ-IS qubits are the lexicographic nine-bit strings `(t1,...,t9)`, with `t9`
varying fastest. For each source subset `S`, the check block is the Kronecker
product, in the order `i=1,...,9`, of `[1,1]` if `i` belongs to `S` and `I_2`
otherwise. Blocks follow the source's ordered subset lists. The saved generator
rows place all X checks first, then all Y checks, then all Z checks. There are
768 rows on 512 qubits. A separate Kronecker-product implementation reconstructs
this entire input without running recognition.

Both inputs use interleaved scalar columns `(x0,z0,x1,z1,...)`. In a packed
integer row, bit `2*i` stores `x_i` and bit `2*i+1` stores `z_i`. The input hash is
SHA-256 of the UTF-8 encoding of the exact Python serialization
`json.dumps([hex(row) for row in rows])`, with default JSON separators and no
trailing newline. This convention is checked independently, not inferred from
the solver's claimed answer. `results/input_certificate_map.json` associates
all 369 inputs with source identifiers, dimensions, input hashes, and outcomes.

## Positive certificates

`q[i] = [a,b,c,d]` stores the matrix `[[a,b],[c,d]]` acting on row labels:
`(x,z) -> (x*a+z*c, x*b+z*d)`. The verifier requires determinant one over GF(2).
`css_x_supports` and `css_z_supports` list the qubit supports of the pure checks.
The verifier constructs the corresponding binary rows and checks both
orthogonality and equality of their combined row space with the transformed
original input space. These are binary LC-to-CSS certificates; the manuscript's
phase-correction procedure supplies signs when a signed realization is needed.

## Negative certificates

`rejected_factor` is the increasing list of original qubit indices in the
rejected direct factor. `local_rref_hex` stores its nonzero RREF rows, with the
factor's qubits relabeled in that order. Pivot columns are ordered increasingly.
The verifier reconstructs the full input row space and checks:

1. Its projection onto the listed qubits has exactly the stored RREF.
2. Every stored factor row, embedded in the full coordinates, belongs to the
   original input space. Together with step 1 this establishes a direct factor.
3. The RREF support graph on those qubits is connected, establishing
   indecomposability by the manuscript's scalar-mask lemma.

For local RREF rows `g_alpha` with pivots `p_beta`, list the nonpivot scalar
columns `j_t` increasingly. Define the ordinary annihilator row `h_t` to have a
one in column `j_t`, entry `g_beta[j_t]` in each pivot column `p_beta`, and zero
elsewhere. At each qubit, let `(x_i,z_i)` and `(u_i,v_i)` be the entries of
`g_alpha` and `h_t`. The reconstructed affine equation has coefficients

```text
a_i: x_i*u_i + z_i*v_i
b_i: x_i*v_i
c_i: z_i*u_i
rhs: sum_i z_i*v_i
```

This is the ordinary-annihilator formula in the manuscript, rather than the
solver's membership-expression construction.

For `N` factor qubits and rank `R`, `equation_indices` uses this exact order:

- Index 0: the homogeneous anchor on the first factor qubit. The field
  `anchor_variable=1` means `b_0=0`; `anchor_variable=2` means `c_0=0`.
- Index `1 + alpha*(2*N-R) + t`: the equation for `(g_alpha,h_t)` above.

The listed indices are the support of the binary inconsistency witness lambda.
The verifier requires their XOR to have zero coefficients and right-hand side
one: `lambda L=0`, `lambda f=1`. It independently checks both anchors. The
solver's `contradiction_verified`, `equations_read`, and timing fields are not
used to establish inconsistency. The saved `lc_css` field selects which kind of
certificate must be checked; it is not accepted without that check.

These identities certify rejection using the two-branch theorem in the
manuscript. They remove the need to trust the recognizer's elimination routine
or repeat its search, while retaining the mathematical theorem and the declared
input construction as the basis for interpreting the certificate.

## Independent verification and timing

Run from this directory using Python 3.12 and NumPy 2.3.5:

```powershell
python -X utf8 verify_all_certificates.py
```

The entry point imports `verify_pbb_certificates.py` but neither recognition
script nor the BB demonstration. It verifies the source catalog and construction hashes, reconstructs
all matrices, and checks 11 positive witnesses and 716 negative branch
certificates on 358 rejected inputs. It exits with an error on a failed check
and refuses Python's `-O` mode because the verifier uses assertions. Verification
does not modify the original classification files or their timing records.

`runtime_and_timing.json` records the runtime, processor, implementation, and
exact inclusion/exclusion boundaries for the original 49.7-second PBB run and
0.663-second XYZ recognition run. Exact start timestamps of those original runs
were not recorded. The all-certificate verification report has its own UTC time,
elapsed time, code hashes, and input/result hashes; it is a separate measurement.

Source hashes establish which snapshot is being checked; the algebraic
certificates establish the classification. The distributed inputs and certificates can be verified offline.
