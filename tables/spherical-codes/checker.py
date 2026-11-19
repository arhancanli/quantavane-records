#!/usr/bin/env python3
"""Exact checker for spherical codes: N points on the unit sphere in R^d whose largest pairwise cosine is at most C.

    python3 checker.py FILE d N C

FILE: one point per line, d decimal coordinates separated by commas (other lines are ignored, as on
spherical-codes.org). Points need not be normalised: the cosine of x and y is <x, y> / (|x| |y|). C is the claimed
largest cosine as a decimal with E digits after the point; the claim holds when the true largest cosine is below
C + 10^-E, that is, when C is the true value rounded to its last stated digit or above it. The table keeps, for each
(d, N), the smallest largest cosine.

  A  exact: every coordinate is an exact decimal, so after scaling by a power of ten each point is an integer vector;
     for each pair it decides <x, y> / (|x| |y|) < C + 10^-E in integers, comparing squares after the signs (for
     C + 10^-E > 0: <x, y> <= 0 or <x, y>^2 10^(2E) < (C 10^E + 1)^2 |x|^2 |y|^2), and reports the largest cosine to
     30 digits.
  B  60-digit decimal arithmetic, written separately: normalises every point, takes every inner product and compares
     the largest with C + 10^-E.

Both also check that there are N points of dimension d, none zero and no two equal. Exit 0 = both PASS.
MIT licence (LICENSE-CODE).
"""
import sys
from decimal import Decimal, getcontext


def read_points(path, d):
    pts = []
    for line in open(path):
        parts = [p.strip() for p in line.strip().split(",")]
        if len(parts) != d:
            continue
        try:
            pts.append([Decimal(p) for p in parts])
        except Exception:
            continue
    return pts


def check_a(pts, C):
    # scale every coordinate to an integer
    D = max(max(-x.as_tuple().exponent, 0) for p in pts for x in p)
    X = [[int(x.scaleb(D)) for x in p] for p in pts]
    for p, q in zip(pts, X):                                  # the scaling is exact
        if any(Decimal(b).scaleb(-D) != a for a, b in zip(p, q)):
            return False, "a coordinate is not an exact decimal"
    E = max(-C.as_tuple().exponent, 0)
    Cn = int(C.scaleb(E))
    lim = (Cn + 1) ** 2
    s10 = 10 ** (2 * E)
    nn = [sum(a * a for a in x) for x in X]
    if min(nn) == 0:
        return False, "a zero vector"
    best = None                                               # (dot, |x|^2 |y|^2, i, j) of the largest cosine
    n = len(X)
    for i in range(n):
        xi, ni = X[i], nn[i]
        for j in range(i + 1, n):
            dot = sum(a * b for a, b in zip(xi, X[j]))
            den = ni * nn[j]
            # the cosine dot / sqrt(den) must be below B = (Cn + 1) / 10^E; compared through squares by sign
            if Cn + 1 > 0:
                bad = dot > 0 and dot * dot * s10 >= lim * den
            else:
                bad = dot >= 0 or dot * dot * s10 <= lim * den
            if bad:
                return False, "points %d and %d: cosine at least C + 10^-%d" % (i + 1, j + 1, E)
            if best is None or _bigger(dot, den, best[0], best[1]):
                best = (dot, den, i + 1, j + 1)
    getcontext().prec = 50
    val = Decimal(best[0]) / Decimal(best[1]).sqrt()
    return True, "largest cosine %s at points %d, %d (below C + 10^-%d)" % (format(val, ".30f"), best[2], best[3], E)


def _bigger(a, den, b, n0):
    """is a / sqrt(den) larger than b / sqrt(n0)? exact: signs first, then squares"""
    if (a >= 0) != (b >= 0):
        return a >= 0
    if a >= 0:
        return a * a * n0 > b * b * den
    return a * a * n0 < b * b * den


def check_b(pts, C):
    getcontext().prec = 60
    unit = []
    for p in pts:
        norm = sum(x * x for x in p).sqrt()
        if norm == 0:
            return False, "a zero vector"
        unit.append([x / norm for x in p])
    E = max(-C.as_tuple().exponent, 0)
    bound = C + Decimal(1).scaleb(-E)
    best = Decimal(-2)
    for i in range(len(unit)):
        u = unit[i]
        for j in range(i + 1, len(unit)):
            c = sum(a * b for a, b in zip(u, unit[j]))
            if c > best:
                best = c
    ok = best < bound
    return ok, "largest cosine %s (%s C + 10^-%d)" % (format(best, ".30f"), "below" if ok else "NOT below", E)


def main(argv):
    try:
        path, d, n, C = argv[0], int(argv[1]), int(argv[2]), Decimal(argv[3])
        pts = read_points(path, d)
    except (IndexError, ValueError, ArithmeticError) as e:
        print("bad input: %s\n%s" % (e, __doc__))
        return 2
    ok = True
    if len(pts) != n:
        print("FAIL  %d points of dimension %d, %d claimed" % (len(pts), d, n))
        ok = False
    if len({tuple(p) for p in pts}) != len(pts):
        print("FAIL  two points are equal")
        ok = False
    for name, fn in (("A", check_a), ("B", check_b)):
        res, why = fn(pts, C)
        print("%s %s  %s" % (name, "PASS" if res else "FAIL", why))
        ok &= res
    print("%d points in dimension %d with largest cosine <= %s: %s" % (n, d, C, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
