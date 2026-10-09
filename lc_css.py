"""Local-Clifford recognition and signed conversion of binary stabilizer codes.

Rows are nonnegative Python integers with interleaved x,z bits, least-significant
qubit first. Local matrices act on row labels from the right. Recognition uses
only the standard library; signed Clifford tables additionally require NumPy.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import itertools

GATES = {
    "I": (1, 0, 0, 1), "H": (0, 1, 1, 0),
    "S": (1, 1, 0, 1), "HS": (1, 1, 1, 0),
    "SH": (0, 1, 1, 1), "HSH": (1, 0, 1, 1),
}
GATE_NAMES = {v: k for k, v in GATES.items()}


def positions(v):
    while v:
        b = v & -v
        yield b.bit_length() - 1
        v ^= b


def rref(rows, width, track=False):
    work = list(rows)
    coeff = [1 << i for i in range(len(work))]
    pivots = []
    for j in range(width):
        k = len(pivots)
        p = next((i for i in range(k, len(work)) if work[i] >> j & 1), None)
        if p is None:
            continue
        work[k], work[p] = work[p], work[k]
        coeff[k], coeff[p] = coeff[p], coeff[k]
        for i in range(len(work)):
            if i != k and work[i] >> j & 1:
                work[i] ^= work[k]
                if track:
                    coeff[i] ^= coeff[k]
        pivots.append(j)
    if track:
        return work[:len(pivots)], pivots, coeff[:len(pivots)]
    return work[:len(pivots)], pivots


def coordinates(row, basis, pivots, coeff=None):
    combination = 0
    for i, (g, p) in enumerate(zip(basis, pivots)):
        if row >> p & 1:
            row ^= g
            combination ^= coeff[i] if coeff is not None else 1 << i
    return None if row else combination


def solve(equations, width):
    """Packed Gaussian elimination, with free variables set to zero."""
    piv = {}
    mask = (1 << width) - 1
    for value in equations:
        while value & mask:
            p = ((value & mask) & -(value & mask)).bit_length() - 1
            if p not in piv:
                piv[p] = value
                break
            value ^= piv[p]
        else:
            if value >> width & 1:
                return None, len(piv)
    answer = 0
    for p in sorted(piv, reverse=True):
        rhs = (piv[p] >> width & 1) ^ ((piv[p] & answer).bit_count() & 1)
        answer |= rhs << p
    return answer, len(piv)


def act(rows, qs):
    out = []
    for row in rows:
        value = 0
        for i, (a, b, c, d) in enumerate(qs):
            x, z = row >> (2*i) & 1, row >> (2*i+1) & 1
            value |= (x*a ^ z*c) << (2*i)
            value |= (x*b ^ z*d) << (2*i+1)
        out.append(value)
    return out


def inverse(qs):
    return [(d, b, c, a) for a, b, c, d in qs]


def split(row, n):
    x = sum(((row >> (2*i)) & 1) << i for i in range(n))
    z = sum(((row >> (2*i+1)) & 1) << i for i in range(n))
    return x, z


def join(x, z, n):
    return sum(((x >> i & 1) << (2*i)) | ((z >> i & 1) << (2*i+1))
               for i in range(n))


def omega(row, n):
    x, z = split(row, n)
    return join(z, x, n)


def symplectic(a, b, n):
    ax, az = split(a, n)
    bx, bz = split(b, n)
    return ((ax & bz).bit_count() ^ (az & bx).bit_count()) & 1


def factors(g, n):
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for row in g:
        support = sorted({j // 2 for j in positions(row)})
        for j in support[1:]:
            parent[root(j)] = root(support[0])
    components = {}
    for i in range(n):
        components.setdefault(root(i), []).append(i)
    return sorted(components.values(), key=lambda c: c[0])


def restrict(row, component):
    return sum(((row >> (2*i)) & 3) << (2*j) for j, i in enumerate(component))


def preserving_equations(g, pivots, n):
    nonpivots = [j for j in range(2*n) if j not in set(pivots)]
    contributors = [[p for row, p in zip(g, pivots) if row >> j & 1]
                    for j in nonpivots]
    eqs = []
    for row in g:
        gp = []
        for i in range(n):
            x, z = row >> (2*i) & 1, row >> (2*i+1) & 1
            gp.extend([(x << (3*i)) ^ (z << (3*i+2)),
                       (x << (3*i+1)) ^ (z << (3*i)) ^ (z << (3*n))])
        for j, used in zip(nonpivots, contributors):
            equation = gp[j]
            for p in used:
                equation ^= gp[p]
            eqs.append(equation)
    return eqs


def witness(a, b, c):
    assert b*c == 0
    def left(v):
        x, z = v
        return a*x ^ b*z, c*x ^ (1 ^ a)*z
    vectors = [(1, 0), (0, 1), (1, 1)]
    u = next(v for v in vectors if left(v) == v)
    w = next(v for v in vectors if left(v) == (0, 0))
    return u[0], w[0], u[1], w[1]


def recognize(rows, n):
    """Return a CSS conversion, or None when no local Clifford map exists.

    The input may contain redundant and zero rows. Invalid binary widths or
    noncommuting labels raise ValueError. Output pure_rows retains nonzero
    projections of the supplied rows; these checks can be redundant.
    """
    rows = validate_rows(rows, n)
    g, _ = rref(rows, 2*n)
    qs = [None]*n
    stats = []
    for component in factors(g, n):
        size = len(component)
        local, pivots = rref([restrict(row, component) for row in g], 2*size)
        equations = preserving_equations(local, pivots, size)
        branches = 0
        for anchor in (1, 2):
            branches += 1
            answer, eq_rank = solve(itertools.chain(equations, [1 << anchor]), 3*size)
            if answer is not None:
                break
        if answer is None:
            return None
        for j, i in enumerate(component):
            qs[i] = witness(answer >> (3*j) & 1, answer >> (3*j+1) & 1,
                            answer >> (3*j+2) & 1)
        stats.append({"qubits": size, "stabilizer_rank": len(local),
                      "variables": 3*size, "equations_without_anchor": len(equations),
                      "branches_solved": branches, "affine_rank": eq_rank})
    transformed = act(rows, qs)
    # Retain nonzero projections of the supplied sparse checks.
    xmask = sum(1 << (2*i) for i in range(n))
    zmask = xmask << 1
    pure_x = list(dict.fromkeys(row & xmask for row in transformed if row & xmask))
    pure_z = list(dict.fromkeys(row & zmask for row in transformed if row & zmask))
    return {"q": qs, "pure_rows": pure_x + pure_z, "x_count": len(pure_x),
            "rank": len(g), "factors": stats}


def sign_tables():
    """Derive signed Clifford conjugation from actual 2x2 unitary matrices."""
    import numpy as np
    eye = np.eye(2, dtype=complex)
    h = np.array([[1, 1], [1, -1]]) / 2**0.5
    s = np.diag([1, 1j])
    gates = {"I": eye, "H": h, "S": s, "HS": h@s, "SH": s@h, "HSH": h@s@h}
    paulis = [eye, np.array([[0, 1], [1, 0]]), np.diag([1, -1]),
              np.array([[0, -1j], [1j, 0]])]
    table = {}
    for name, u in gates.items():
        labels = act(list(range(4)), [GATES[name]])
        table[name] = []
        for label, target in enumerate(labels):
            val = u @ paulis[label] @ u.conj().T
            sign = 0 if np.allclose(val, paulis[target]) else 1
            assert np.allclose(val, (-1)**sign * paulis[target])
            table[name].append(sign)
    return table


def transform_signs(rows, signs, names, table):
    return [sign ^ (sum(table[name][row >> (2*i) & 3] for i, name in enumerate(names)) & 1)
            for row, sign in zip(rows, signs)]


def signed_product(rows, signs, combination, n):
    label, exponent = 0, 0
    for i in positions(combination):
        x, z = split(label, n)
        u, v = split(rows[i], n)
        exponent += 2*signs[i] + (x & z).bit_count() + (u & v).bit_count()
        exponent += -((x ^ u) & (z ^ v)).bit_count() + 2*(z & u).bit_count()
        label ^= rows[i]
    assert exponent % 2 == 0
    return label, exponent % 4 // 2


def verify_and_transport(rows, input_signs, recovered, n, table):
    rows = validate_rows(rows, n)
    if len(input_signs) != len(rows) or any(type(s) is not int or s not in (0, 1) for s in input_signs):
        raise ValueError("Supply one binary sign for every input row")
    qs, pure = recovered["q"], recovered["pure_rows"]
    trans = act(rows, qs)
    bg, bp, bc = rref(trans, 2*n, track=True)
    assert rref(pure, 2*n)[0] == bg
    split_rows = [split(row, n) for row in pure]
    hx = [x for x, z in split_rows if x]
    hz = [z for x, z in split_rows if z]
    assert all((x & z).bit_count() % 2 == 0 for x in hx for z in hz)
    assert len(rref(hx, n)[0]) + len(rref(hz, n)[0]) == recovered["rank"]
    # K is expressed against the supplied redundant input rows directly.
    k = [coordinates(c, bg, bp, bc) for c in pure]
    assert all(v is not None for v in k)
    trans_signs = transform_signs(rows, input_signs, [GATE_NAMES[q] for q in qs], table)
    output_signs = []
    for c, combination in zip(pure, k):
        value, sign = signed_product(trans, trans_signs, combination, n)
        assert value == c
        output_signs.append(sign)
    correction, _ = solve([omega(c, n) | (sign << (2*n))
                           for c, sign in zip(pure, output_signs)], 2*n)
    assert correction is not None
    assert all(symplectic(c, correction, n) == sign for c, sign in zip(pure, output_signs))
    # Verify the resulting signed group against every redundant original row.
    cg, cp, cc = rref(pure, 2*n, track=True)
    for row, sign in zip(trans, trans_signs):
        combination = coordinates(row, cg, cp, cc)
        rebuilt, rebuilt_sign = signed_product(pure, [0]*len(pure), combination, n)
        assert rebuilt == row
        assert rebuilt_sign == sign ^ symplectic(row, correction, n)
    degrees = [sum(bool(row >> (2*i) & 3) for row in pure) for i in range(n)]
    weights = [sum(bool(row >> (2*i) & 3) for i in range(n)) for row in pure]
    return {"hx": hx, "hz": hz, "syndrome_map": k, "pauli_correction": correction,
            "max_weight": max(weights, default=0), "weight_set": sorted(set(weights)),
            "max_degree": max(degrees, default=0), "degree_set": sorted(set(degrees)),
            "x_rank": len(rref(hx, n)[0]), "z_rank": len(rref(hz, n)[0]),
            "precorrection_signs": output_signs}


def validate_rows(rows, n):
    """Validate packed binary labels and isotropy, preserving the supplied rows."""
    if type(n) is not int or n < 0:
        raise ValueError("n must be a nonnegative integer")
    rows = list(rows)
    if any(type(row) is not int or row < 0 or row.bit_length() > 2*n for row in rows):
        raise ValueError("Each row must be a nonnegative integer fitting in 2*n bits")
    basis, _ = rref(rows, 2*n)
    xmask = sum(1 << (2*i) for i in range(n))
    paired = [((row & xmask) << 1) | ((row >> 1) & xmask) for row in basis]
    if any((a & paired[j]).bit_count() & 1
           for i, a in enumerate(basis) for j in range(i+1, len(basis))):
        raise ValueError("Input stabilizer labels must commute")
    return rows
