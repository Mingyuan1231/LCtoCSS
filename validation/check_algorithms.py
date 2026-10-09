"""Check the manuscript's row convention, graph construction, and phases.

Run from the repository root: python validation/check_algorithms.py
Uses the standard library. The row-form implementation below is separate
from the reference implementation, whose exact-search answers are compared.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O or -OO: mathematical checks require assertions.")


import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproducibility import report_output, save_report

import hashlib
import itertools
import json
import platform
from datetime import datetime, timezone
import validation_reference as ref

HERE = Path(__file__).resolve().parent


def reduce_rows(rows, width):
    a = [list(row) for row in rows]
    pivots = []
    for j in range(width):
        p = next((i for i in range(len(pivots), len(a)) if a[i][j]), None)
        if p is None:
            continue
        k = len(pivots)
        a[k], a[p] = a[p], a[k]
        for i in range(len(a)):
            if i != k and a[i][j]:
                a[i] = [x ^ y for x, y in zip(a[i], a[k])]
        pivots.append(j)
    return a[:len(pivots)], pivots


def bits(rows, width):
    return [[(v >> j) & 1 for j in range(width)] for v in rows]


def packed(rows):
    return [sum(x << j for j, x in enumerate(row)) for row in rows]


def solve(equations, width):
    a, pivots = reduce_rows(bits(equations, width + 1), width + 1)
    if width in pivots:
        return None
    x = [0] * width
    for row, p in zip(a, pivots):
        x[p] = row[-1]
    return x


def factors(g, n):
    adjacency = [set() for _ in range(n)]
    for row in g:
        p = row.index(1) // 2
        for j, value in enumerate(row):
            if value:
                adjacency[p].add(j // 2)
                adjacency[j // 2].add(p)
    unseen = set(range(n))
    result = []
    while unseen:
        todo = [min(unseen)]
        component = set()
        while todo:
            i = todo.pop()
            if i in component:
                continue
            component.add(i)
            todo.extend(adjacency[i] - component)
        unseen -= component
        result.append(sorted(component))
    return result


def affine_system(g, n):
    """Each GP entry is an affine expression, packed with its constant last."""
    pivots = [row.index(1) for row in g]
    nonpivots = [j for j in range(2*n) if j not in pivots]
    eqs = []
    for row in g:
        gp = []
        for i in range(n):
            x, z = row[2*i:2*i+2]
            gp += [(x << (3*i)) ^ (z << (3*i+2)),
                   (x << (3*i+1)) ^ (z << (3*i)) ^ (z << (3*n))]
        for j in nonpivots:
            expression = gp[j]
            for beta, p in enumerate(pivots):
                if g[beta][j]:
                    expression ^= gp[p]
            eqs.append(expression)
    # Compare the exact row equations with the ordinary annihilator formula.
    other = []
    for row in g:
        for j in nonpivots:
            h = [0]*(2*n)
            h[j] = 1
            for beta, p in enumerate(pivots):
                h[p] = g[beta][j]
            expression = 0
            for i in range(n):
                x, z = row[2*i:2*i+2]
                u, v = h[2*i:2*i+2]
                expression ^= ((x*u) ^ (z*v)) << (3*i)
                expression ^= (x*v) << (3*i+1)
                expression ^= (z*u) << (3*i+2)
                expression ^= (z*v) << (3*n)
            other.append(expression)
    assert eqs == other
    return eqs


def multiply(q, r):
    a,b,c,d = q
    e,f,g,h = r
    return (a*e ^ b*g, a*f ^ b*h, c*e ^ d*g, c*f ^ d*h)


def act(rows, qs):
    out = []
    for row in rows:
        changed = []
        for i, (a,b,c,d) in enumerate(qs):
            x,z = row[2*i:2*i+2]
            changed += [x*a ^ z*c, x*b ^ z*d]
        out.append(changed)
    return out


def witness(a,b,c):
    p = (a,b,c,1^a)
    columns = [(1,0),(0,1),(1,1)]
    def left(v):
        x,z = v
        return (a*x ^ b*z, c*x ^ (1^a)*z)
    u = next(v for v in columns if left(v) == v)
    w = next(v for v in columns if left(v) == (0,0))
    q = (u[0],w[0],u[1],w[1])
    assert multiply(p,q) == multiply(q,(1,0,0,0))
    return q


def check_css(g, qs, n):
    t, _ = reduce_rows(act(g,qs),2*n)
    hx, px = reduce_rows([row[::2] for row in t],n)
    hz, pz = reduce_rows([row[1::2] for row in t],n)
    pure = [[v for x in row for v in (x,0)] for row in hx]
    pure += [[v for z in row for v in (0,z)] for row in hz]
    assert reduce_rows(pure,2*n)[0] == t
    assert len(px)+len(pz) == len(g)
    assert all(sum(x*z for x,z in zip(a,b)) % 2 == 0 for a in hx for b in hz)
    return hx,hz


def row_decide(raw, n):
    g, _ = reduce_rows(bits(raw,2*n),2*n)
    qs = [None]*n
    for component in factors(g,n):
        local, _ = reduce_rows([[row[j] for i in component for j in (2*i,2*i+1)] for row in g],2*len(component))
        eqs = affine_system(local,len(component))
        x = solve(eqs+[1 << 1],3*len(component))
        if x is None:
            x = solve(eqs+[1 << 2],3*len(component))
        if x is None:
            return None
        for j,i in enumerate(component):
            a,b,c = x[3*j:3*j+3]
            assert b*c == 0
            qs[i] = witness(a,b,c)
    check_css(g,qs,n)
    return qs


def graph_decide(a):
    n = len(a)
    qs = [None]*n
    rows = [sum((int(i==j) << (2*j)) | (a[i][j] << (2*j+1)) for j in range(n)) for i in range(n)]
    g, _ = reduce_rows(bits(rows,2*n),2*n)
    unseen = set(range(n))
    components = []
    while unseen:
        todo = [min(unseen)]
        component = set()
        while todo:
            i = todo.pop()
            if i in component:
                continue
            component.add(i)
            todo.extend(j for j in range(n) if a[i][j] and j not in component)
        unseen -= component
        components.append(sorted(component))
    for component in components:
        av = [[a[i][j] for j in component] for i in component]
        size = len(component)
        eqs = []
        for i in range(size):
            for j in range(i+1,size):
                eq = aij = av[i][j]
                eq = (aij << i) ^ (aij << j) ^ (aij << (2*size))
                for k in range(size):
                    eq ^= (av[i][k]*av[k][j]) << (size+k)
                eqs.append(eq)
        x = solve(eqs+[1 << size],2*size)
        if x is None:
            anchor = sum(av[0][j] << (size+j) for j in range(size))
            x = solve(eqs+[anchor],2*size)
        if x is None:
            return None
        aa, u = x[:size], x[size:]
        b = [sum(av[i][j]*u[j] for j in range(size)) % 2 for i in range(size)]
        for j,i in enumerate(component):
            qs[i] = witness(aa[j],b[j],u[j])
    hx,hz = check_css(g,qs,n)
    # RREF H_X is systematic on its pivot information set.
    hx, information = reduce_rows(hx,n)
    complement = [i for i in range(n) if i not in information]
    b = [[0]*n for _ in range(n)]
    for alpha,i in enumerate(information):
        for j in complement:
            b[i][j] = b[j][i] = hx[alpha][j]
    total = [multiply(q,(0,1,1,0)) if i in complement else q for i,q in enumerate(qs)]
    expected = [[v for j in range(n) for v in (int(i==j),b[i][j])] for i in range(n)]
    assert reduce_rows(act(g,total),2*n)[0] == reduce_rows(expected,2*n)[0]
    return {"adjacency": b, "local_row_matrices": total, "information_set": information}


def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def dagger(a):
    return [[a[j][i].conjugate() for j in range(2)] for i in range(2)]


def near(a,b):
    return all(abs(a[i][j]-b[i][j]) < 1e-12 for i in range(2) for j in range(2))


def phase_checks():
    identity = [[1+0j,0j],[0j,1+0j]]
    xx = [[0j,1+0j],[1+0j,0j]]
    zz = [[1+0j,0j],[0j,-1+0j]]
    yy = [[0j,-1j],[1j,0j]]
    pauli = {(0,0):identity,(1,0):xx,(0,1):zz,(1,1):yy}
    h = [[1/2**.5,1/2**.5],[1/2**.5,-1/2**.5]]
    s = [[1+0j,0j],[0j,1j]]
    table = [("I",identity,(1,0,0,1)),("H",h,(0,1,1,0)),("S",s,(1,1,0,1)),
             ("HS",mm(h,s),(1,1,1,0)),("SH",mm(s,h),(0,1,1,1)),("HSH",mm(mm(h,s),h),(1,0,1,1))]
    signed_images = {"I":((1,0,0),(0,1,0)),"H":((0,1,0),(1,0,0)),
                     "S":((1,1,0),(0,1,0)),"HS":((1,1,1),(1,0,0)),
                     "SH":((0,1,0),(1,1,0)),"HSH":((1,0,0),(1,1,1))}
    for name,u,q in table:
        for p,(x,z,sign) in zip((xx,zz),signed_images[name]):
            target = [[(-1)**sign*v for v in row] for row in pauli[(x,z)]]
            assert near(mm(mm(u,p),dagger(u)),target)
        for (x,z),p in pauli.items():
            ex,ez = act([[x,z]],[q])[0]
            transformed = mm(mm(u,p),dagger(u))
            target = pauli[(ex,ez)]
            assert near(transformed,target) or near(transformed,[[-v for v in row] for row in target])
        if name in ("H","S"):
            for (x,z),p in pauli.items():
                ex,ez = act([[x,z]],[q])[0]
                target = [[(-1)**(x*z)*v for v in row] for row in pauli[(ex,ez)]]
                assert near(mm(mm(u,p),dagger(u)),target)
    for (x,z),p in pauli.items():
        for (u,v),q in pauli.items():
            exponent = x*z+u*v-(x^u)*(z^v)+2*z*u
            target = [[1j**(exponent % 4)*val for val in row] for row in pauli[(x^u,z^v)]]
            assert near(mm(p,q),target)
    return {"clifford_label_cases":24,"table_signed_images":12,"elementary_signed_cases":8,"pauli_product_cases":16}


def main():
    report_path = report_output('algorithms.json')
    rerun = json.loads((HERE/'cases.json').read_text(encoding='utf-8'))
    import lc_css as core
    for case in rerun['cases']:
        n = case['n']
        result = row_decide(case['input_rows'],n)
        production = core.recognize(case['input_rows'], n)
        assert (production is not None) == case['oracle']['equivalent'], case['name']
        if production is not None:
            assert core.rref(core.act(case['input_rows'], production['q']), 2*n)[0] == core.rref(production['pure_rows'], 2*n)[0]
        assert (result is not None) == case['oracle']['equivalent'],case['name']
        if result is not None:
            g,_ = reduce_rows(bits(case['input_rows'],2*n),2*n)
            transposed = [(a,c,b,d) for a,b,c,d in case['projector']['local_matrices']]
            check_css(g,transposed,n)
    graphs = []
    for n in range(6):
        edges = list(itertools.combinations(range(n),2))
        count = positive = 0
        for mask in range(1 << len(edges)):
            selected = [edge for i,edge in enumerate(edges) if mask >> i & 1]
            a = [[0]*n for _ in range(n)]
            for i,j in selected:
                a[i][j] = a[j][i] = 1
            result = graph_decide(a)
            oracle = ref.exhaustive_local_oracle(ref.graph_rows(n,selected),n)
            assert (result is not None) == oracle['equivalent'],(n,mask)
            positive += result is not None
            count += 1
        graphs.append({"n":n,"cases":count,"positive":positive})
    triangle = [[0,1,1],[1,0,1],[1,1,0]]
    triangle_q = [(1,0,1,1),(1,1,1,0),(1,1,1,0)]
    triangle_g = bits(ref.graph_rows(3,[(0,1),(0,2),(1,2)]),6)
    check_css(triangle_g,triangle_q,3)
    assert graph_decide(triangle) is not None
    c5 = [[int((i-j)%5 in (1,4)) for j in range(5)] for i in range(5)]
    assert graph_decide(c5) is None
    # The explicit triangle graph witness includes the stated Pauli correction.
    # Full state-vector conjugation verifies phases, not just binary labels.
    h = [[1/2**.5,1/2**.5],[1/2**.5,-1/2**.5]]
    s = [[1+0j,0j],[0j,1j]]
    z = [[1+0j,0j],[0j,-1+0j]]
    gates = [mm(mm(h,s),h),mm(z,s),mm(z,s)]
    psi = [(-1)**sum(((v>>i)&1)*((v>>j)&1) for i,j in ((0,1),(0,2),(1,2)))/8**.5 for v in range(8)]
    for i,u in enumerate(gates):
        out = [0j]*8
        for v,amp in enumerate(psi):
            for bit in (0,1):
                w = (v & ~(1<<i)) | (bit<<i)
                out[w] += u[bit][(v>>i)&1]*amp
        psi = out
    target = [(-1)**sum(((v>>i)&1)*((v>>j)&1) for i,j in ((0,1),(0,2)))/8**.5 for v in range(8)]
    phase = psi[0]/target[0]
    assert all(abs(x-phase*y)<1e-12 for x,y in zip(psi,target))
    report = {"executed_at_utc":datetime.now(timezone.utc).isoformat(),"python":platform.python_version(),
              "reference_cases_checked":len(rerun['cases']),"row_convention_cases":len(rerun['cases']),
              "reference_transpose_witnesses_checked":rerun['positive_witnesses_verified'],
              "graph_cases":graphs,"phase_checks":phase_checks(),"triangle_signed_graph_witness_verified":True,
              "triangle_css_row_matrices":triangle_q,"triangle_graph_output":graph_decide(triangle),
              "mismatches":0,"source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'validation_reference.py',HERE/'cases.json',HERE.parent/'lc_css.py']}}
    save_report(report, report_path)
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
