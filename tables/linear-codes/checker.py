#!/usr/bin/env python3
"""Exact checker for linear codes: an [n, k, d]_q code (minimum distance at least d).

    python3 checker.py FILE q n k d [--modulus c0,c1,...] [--exact]

FILE: a generator matrix, k lines of n field elements. For q prime the elements are 0..q-1; for q = p^e they are
integers c0 + p c1 + p^2 c2 + ... giving the coefficients in the polynomial basis of the Conway polynomial (built in
for q = 4, 8, 9, 16, 25, 27, 32, 49, 64; --modulus gives the low coefficients of another monic polynomial). The table
(codetables.de, MinT, Chen's database) keeps the largest d for each (q, n, k). Two separately written routines:
  A  every one of the q^k messages, in odometer order, the codeword updated by adding the change in one coordinate
     times its generator row (field tables by shift-and-add); a nonzero message giving the zero word means the rows
     are dependent and fails.
  B  one message per projective point (first nonzero coordinate 1), each codeword computed from scratch (field
     tables by schoolbook products); scalar multiples have equal weight, so this covers every nonzero codeword.
Both report the minimum weight; the check passes when both find it at least d (--exact: equal to d).
Exit 0 = both PASS. Pure Python: practical up to q^k of a few million. MIT licence (LICENSE-CODE).
"""
import sys

CONWAY = {4: (2, [1, 1]), 8: (2, [1, 1, 0]), 16: (2, [1, 1, 0, 0]), 32: (2, [1, 0, 1, 0, 0]),
          64: (2, [1, 1, 0, 1, 1, 0]), 9: (3, [2, 2]), 27: (3, [1, 2, 0]), 25: (5, [2, 4]), 49: (7, [3, 6])}


def _prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def field(q, modulus=None):
    """(add, mul, p) tables of GF(q); elements 0..q-1 in the coding c0 + p c1 + ..."""
    if _prime(q):
        return [[(a + b) % q for b in range(q)] for a in range(q)], [[a * b % q for b in range(q)] for a in range(q)], q
    if modulus is None:
        if q not in CONWAY:
            raise ValueError("give --modulus for q = %d" % q)
        p, f = CONWAY[q]
    else:
        f = modulus
        p = next(p for p in range(2, q + 1) if _prime(p) and q % p == 0)
    e = len(f)
    assert p ** e == q, "modulus degree does not match q"

    def dig(x):
        return [(x // p ** i) % p for i in range(e)]

    def num(ds):
        return sum(c * p ** i for i, c in enumerate(ds))
    add = [[num([(x + y) % p for x, y in zip(dig(a), dig(b))]) for b in range(q)] for a in range(q)]
    # multiplication: a * x^i by shifting and reducing x^e = -(f0 + f1 x + ...), then linearity in b's digits
    shifts = []
    for a in range(q):
        row, cur = [], dig(a)
        for _ in range(e):
            row.append(num(cur))
            top = cur[-1]
            cur = [0] + cur[:-1]
            cur = [(c - top * fi) % p for c, fi in zip(cur, f)]
        shifts.append(row)
    mul = [[0] * q for _ in range(q)]
    for a in range(q):
        for b in range(q):
            acc = 0
            for i, bi in enumerate(dig(b)):
                for _ in range(bi):
                    acc = add[acc][shifts[a][i]]
            mul[a][b] = acc
    return add, mul, p


def field_b(q, modulus=None):
    """multiplication table of GF(q) for routine B, written separately: schoolbook product of the digit polynomials,
    then reduction of the high terms from the top"""
    if _prime(q):
        return [[a * b % q for b in range(q)] for a in range(q)]
    p, f = (CONWAY[q] if modulus is None else (next(p for p in range(2, q + 1) if _prime(p) and q % p == 0), modulus))
    e = len(f)
    mul = [[0] * q for _ in range(q)]
    for a in range(q):
        da = [(a // p ** i) % p for i in range(e)]
        for b in range(q):
            db = [(b // p ** i) % p for i in range(e)]
            prod = [0] * (2 * e - 1)
            for i, x in enumerate(da):
                for j, y in enumerate(db):
                    prod[i + j] = (prod[i + j] + x * y) % p
            for d in range(2 * e - 2, e - 1, -1):         # x^d = x^(d-e) * x^e = -x^(d-e) (f0 + f1 x + ...)
                c, prod[d] = prod[d], 0
                for i in range(e):
                    prod[d - e + i] = (prod[d - e + i] - c * f[i]) % p
            mul[a][b] = sum(c * p ** i for i, c in enumerate(prod[:e]))
    return mul


def read(path, q, n, k):
    rows = [[int(x) for x in l.replace(",", " ").split()] for l in open(path)
            if l.split() and l.split()[0].lstrip("-").isdigit()]
    if len(rows) != k or any(len(r) != n for r in rows) or any(not 0 <= x < q for r in rows for x in r):
        raise ValueError("expected %d rows of %d elements in 0..%d" % (k, n, q - 1))
    return rows


def check_a(G, F, q, n, k):
    add, mul, p = F
    neg = [next(b for b in range(q) if add[a][b] == 0) for a in range(q)]
    # multiples[i][c] = c * row i
    mult = [[[mul[c][x] for x in row] for c in range(q)] for row in G]
    msg = [0] * k
    cw = [0] * n
    best = n + 1
    total = q ** k
    for step in range(1, total):
        i = 0                                    # odometer: the lowest coordinate that does not wrap
        while msg[i] == q - 1:
            old = mult[i][q - 1]
            cw = [add[a][neg[b]] for a, b in zip(cw, old)]
            msg[i] = 0
            i += 1
        old, new = mult[i][msg[i]], mult[i][msg[i] + 1]
        cw = [add[add[a][neg[b]]][c] for a, b, c in zip(cw, old, new)]
        msg[i] += 1
        wt = n - cw.count(0)
        if wt == 0:
            return False, 0, "a nonzero message gives the zero word: the rows are dependent"
        if wt < best:
            best = wt
    return True, best, "%d codewords, minimum weight %d" % (total - 1, best)


def check_b(G, mulb, q, n, k):
    p = next(p for p in range(2, q + 1) if q % p == 0)
    e = 1
    while p ** e < q:
        e += 1

    def addf(a, b):
        return sum((((a // p ** i) + (b // p ** i)) % p) * p ** i for i in range(e))
    add = [[addf(a, b) for b in range(q)] for a in range(q)]
    best, count = n + 1, 0
    for lead in range(k):
        for rest in range(q ** (k - 1 - lead)):
            coeffs = [0] * lead + [1]
            r = rest
            for _ in range(k - 1 - lead):
                coeffs.append(r % q)
                r //= q
            cw = [0] * n
            for c, row in zip(coeffs, G):
                if c:
                    mc = mulb[c]
                    cw = [add[a][mc[x]] for a, x in zip(cw, row)]
            wt = n - cw.count(0)
            count += 1
            if wt == 0:
                return False, 0, "a nonzero message gives the zero word"
            best = min(best, wt)
    return True, best, "%d projective codewords, minimum weight %d" % (count, best)


def main(argv):
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    try:
        path, q, n, k, d = args[0], int(args[1]), int(args[2]), int(args[3]), int(args[4])
        mod = [int(x) for x in argv[argv.index("--modulus") + 1].split(",")] if "--modulus" in argv else None
        F, Fb = field(q, mod), field_b(q, mod)
        G = read(path, q, n, k)
    except (IndexError, ValueError, AssertionError) as e:
        print("bad input: %s (see the header of this file)" % (e,))
        return 2
    ok = True
    exact = "--exact" in argv
    res = {}
    for name, fn, fl in (("A", check_a, F), ("B", check_b, Fb)):
        good, dmin, why = fn(G, fl, q, n, k)
        passed = good and (dmin == d if exact else dmin >= d)
        print("%s %s  %s" % (name, "PASS" if passed else "FAIL", why))
        res[name] = dmin
        ok &= passed
    if res["A"] != res["B"]:
        print("FAIL  the routines disagree")
        ok = False
    print("[%d,%d,%d]_%d: %s" % (n, k, d, q, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
