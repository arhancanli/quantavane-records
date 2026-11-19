#!/usr/bin/env python3
"""Exact checker for Ramsey lower bounds: a colouring of the edges of K_n with no forbidden monochromatic clique.

    python3 checker.py FILE s t                 two colours: R(s, t) > n (no clique of size s in one colour and none
                                                of size t in the other, either way round)
    python3 checker.py FILE --forbid 1:3,2:3,0:16    any number of colours: colour symbol c has no clique of size s_c

FILE: n lines of n symbols (the n x n colour matrix; the diagonal is ignored but must be the same symbol 0); the
layout of github.com/arhancanli/quantavane-ramsey-records. The matrix must be symmetric. When it is circulant
(entry (i, j) depends only on j - i mod n), translations are automorphisms of every colour class, so it is enough to
search cliques through vertex 0; the program checks circulance before using it. Two separately written routines:
  A  branch and bound on bitsets with greedy-colouring bounds (in the style of Tomita's MCQ): is there a clique of
     size s in this colour?
  B  Bron-Kerbosch with pivoting on Python sets, stopped at the first clique of size s; for s = 3 it checks every
     edge's two neighbourhoods for a common vertex instead.
The repository above also has check.c, a third, compiled checker for two colours. Exit 0 = both PASS.
Pure Python: quick for triangles and small cliques; large cliques in dense colours (say K_15 on 280 vertices) can
take hours, where check.c takes minutes. MIT licence (LICENSE-CODE).
"""
import sys


def read(path):
    rows = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    n = len(rows)
    if any(len(r) != n for r in rows):
        raise ValueError("not an n x n matrix")
    for i in range(n):
        for j in range(i + 1, n):
            if rows[i][j] != rows[j][i]:
                raise ValueError("not symmetric at (%d, %d)" % (i, j))
    return rows


def circulant(rows):
    n = len(rows)
    return all(rows[i][j] == rows[0][(j - i) % n] for i in range(n) for j in range(n) if i != j)


def colour_bits(rows, c):
    n = len(rows)
    return [sum(1 << j for j in range(n) if j != i and rows[i][j] == c) for i in range(n)]


# ------------------------------------------------------------------------------------------- A: bitset B&B
def has_clique_a(adj, cand, need):
    """is there a clique of size need inside the vertex set cand (a bitmask)?"""
    if need <= 0:
        return True
    if bin(cand).count("1") < need:
        return False
    # greedy colouring of cand: colour classes are independent sets
    order, bound = [], []
    U, k = cand, 0
    while U:
        k += 1
        Q = U
        while Q:
            v = (Q & -Q).bit_length() - 1
            Q &= ~(1 << v)
            Q &= ~adj[v]
            U &= ~(1 << v)
            order.append(v)
            bound.append(k)
    if k < need:
        return False
    for idx in range(len(order) - 1, -1, -1):
        if bound[idx] < need:
            return False
        v = order[idx]
        if has_clique_a(adj, cand & adj[v], need - 1):
            return True
        cand &= ~(1 << v)
    return False


def check_a(rows, forbid, circ):
    n = len(rows)
    for c, s in forbid.items():
        adj = colour_bits(rows, c)
        starts = [0] if circ else range(n)
        allv = (1 << n) - 1
        for v in starts:
            # cliques through v, among vertices after v (or all neighbours of 0 when circulant)
            cand = adj[v] if circ else adj[v] & ~((1 << (v + 1)) - 1)
            if has_clique_a(adj, cand & allv, s - 1):
                return False, "colour %s has a clique of size %d (through vertex %d)" % (c, s, v)
    return True, "no forbidden clique (%s)" % ", ".join("colour %s: < %d" % kv for kv in sorted(forbid.items()))


# ------------------------------------------------------------------------------------------- B: Bron-Kerbosch
def check_b(rows, forbid, circ):
    n = len(rows)
    sys.setrecursionlimit(100000)
    for c, s in forbid.items():
        N = [set(j for j in range(n) if j != i and rows[i][j] == c) for i in range(n)]
        if s == 3:
            for u in range(n):
                for v in N[u]:
                    if v > u and N[u] & N[v]:
                        return False, "colour %s has a triangle %d, %d, %d" % (c, u, v, min(N[u] & N[v]))
            continue
        found = []

        def bk(R, P, X):
            if len(R) >= s:
                found.append(list(R))
                return True
            if len(R) + len(P) < s:
                return False
            pivot = max(P | X, key=lambda u: len(N[u] & P))
            for v in list(P - N[pivot]):
                if bk(R + [v], P & N[v], X & N[v]):
                    return True
                P = P - {v}
                X = X | {v}
            return False
        if circ:
            hit = bk([0], set(N[0]), set())
        else:
            hit = bk([], set(range(n)), set())
        if hit:
            return False, "colour %s has a clique of size %d: %s" % (c, s, found[0][:s])
    return True, "no forbidden clique (Bron-Kerbosch)"


def main(argv):
    try:
        path = argv[0]
        rows = read(path)
        symbols = sorted({ch for i, r in enumerate(rows) for j, ch in enumerate(r) if i != j})
        if "--forbid" in argv:
            spec = argv[argv.index("--forbid") + 1]
            options = [{a: int(b) for a, b in (x.split(":") for x in spec.split(","))}]
        else:
            s, t = int(argv[1]), int(argv[2])
            if len(symbols) != 2:
                raise ValueError("two colours expected, found %s" % symbols)
            options = [{symbols[1]: s, symbols[0]: t}, {symbols[1]: t, symbols[0]: s}]
    except (IndexError, ValueError) as e:
        print("bad input: %s (see the header of this file)" % (e,))
        return 2
    n = len(rows)
    circ = circulant(rows)
    print("colouring of K_%d, colours %s, %scirculant" % (n, ",".join(symbols), "" if circ else "not "))
    for forbid in options:
        res = []
        for name, fn in (("A", check_a), ("B", check_b)):
            ok, why = fn(rows, forbid, circ)
            res.append(ok)
            print("%s %s  %s" % (name, "PASS" if ok else "FAIL", why))
        if all(res):
            print("lower bound n + 1 = %d: PASS" % (n + 1))
            return 0
        if any(res):
            print("FAIL  the routines disagree")
            return 1
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
