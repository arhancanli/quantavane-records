#!/usr/bin/env python3
"""Exact checker for circulant weighing matrices CW(n, k).

A CW(n, k) is a circulant n x n matrix W with entries in {0, 1, -1} and W W^T = k I. It is proper when its first row
is not x^s w(x^d) for a CW(n/d, k) row w and some d > 1 dividing n, that is, when for every prime d dividing n the
support is not inside one residue class mod d. The table (Gordon's La Jolla repository) records which (n, k) exist.

    python3 checker.py FILE n k [--proper]

FILE: the first row (n integers in {0, 1, -1} on one line), or the whole matrix (n such lines).
  A  builds the n x n circulant (or reads the matrix and checks it is circulant) and computes W W^T entry by entry.
  B  written separately: the periodic autocorrelation of the first row, sum_i w_i w_(i+s), is 0 for s = 1..n-1 and
     k for s = 0, computed as the product w(x) w(x^-1) in Z[x]/(x^n - 1).
Both check the entries and that exactly k entries of the row are nonzero. --proper also checks properness.
Exit 0 = both PASS. MIT licence (LICENSE-CODE).
"""
import sys


def read(path, n):
    rows = [[int(x) for x in l.replace(",", " ").split()] for l in open(path) if l.strip() and not l.startswith("#")]
    if any(len(r) != n for r in rows) or not rows or any(x not in (-1, 0, 1) for r in rows for x in r):
        raise ValueError("expected rows of %d entries in {0, 1, -1}" % n)
    return rows


def check_a(rows, n, k):
    if len(rows) == 1:
        r = rows[0]
        M = [[r[(j - i) % n] for j in range(n)] for i in range(n)]
    elif len(rows) == n:
        M = rows
        if any(M[i][j] != M[0][(j - i) % n] for i in range(n) for j in range(n)):
            return False, "the matrix is not circulant"
    else:
        return False, "give one row or all n rows"
    for i in range(n):
        for j in range(i, n):
            s = sum(a * b for a, b in zip(M[i], M[j]))
            if s != (k if i == j else 0):
                return False, "(W W^T)[%d][%d] = %d" % (i, j, s)
    return True, "W W^T = %d I" % k


def check_b(rows, n, k):
    w = rows[0]
    # w(x) * w(x^-1) mod x^n - 1: coefficient of x^s is sum_i w_(i+s) w_i
    prod = [0] * n
    sup = [i for i, x in enumerate(w) if x]
    for i in sup:
        for j in sup:
            prod[(j - i) % n] += w[i] * w[j]
    ok = prod[0] == k and not any(prod[1:])
    return ok, "periodic autocorrelation %s" % ("0 at every shift, %d at 0" % k if ok else "nonzero at shift %d"
                                                % next(s for s in range(1, n) if prod[s]) if any(prod[1:]) else
                                                "%d at 0" % prod[0])


def proper(w, n):
    sup = [i for i, x in enumerate(w) if x]
    for d in range(2, n + 1):
        if n % d == 0 and all(d % p for p in range(2, int(d ** 0.5) + 1)):
            if len({i % d for i in sup}) == 1:
                return False, d
    return True, None


def main(argv):
    try:
        path, n, k = argv[0], int(argv[1]), int(argv[2])
        rows = read(path, n)
    except (IndexError, ValueError) as e:
        print("bad input: %s\n%s" % (e, __doc__))
        return 2
    ok = sum(1 for x in rows[0] if x) == k
    if not ok:
        print("FAIL  the row has %d nonzero entries, not %d" % (sum(1 for x in rows[0] if x), k))
    for name, fn in (("A", check_a), ("B", check_b)):
        res, why = fn(rows, n, k)
        print("%s %s  %s" % (name, "PASS" if res else "FAIL", why))
        ok &= res
    if "--proper" in argv:
        pr, d = proper(rows[0], n)
        print("%s  %s" % ("PASS" if pr else "FAIL", "proper" if pr else "support in one class mod %d: not proper" % d))
        ok &= pr
    print("CW(%d, %d): %s" % (n, k, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
