#!/usr/bin/env python3
"""Exact checker for covering arrays and Sherwood covering perfect hash families.

A covering array CA(N; t, k, v) is an N x k array over {0, ..., v-1} in which every N x t subarray contains every
t-tuple at least once. The table keeps, for each (t, k, v), the fewest rows CAN(t, k, v).

    python3 checker.py ca   FILE t v [--rows N] [--cols K]
    python3 checker.py cphf  FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]
    python3 checker.py scphf FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]

ca:   FILE has N lines of k symbols separated by spaces.
  A  for every t columns, the set of t-tuples read down the rows has v^t elements.
  B  for every symbol of every column, the set of rows holding it as a bitmask; for every t columns and every
     t-tuple, the AND of the t masks must be nonzero.
cphf: FILE has n lines of k integers (the layout of the ca-tools table; lines not starting with a digit are skipped): each integer is a vector of GF(q)^t, written
     in base q with the most significant digit first, each digit a field element in the Conway coding
     (c0 + p c1 + p^2 c2 + ... for the polynomial basis). Every t columns must have a row whose t vectors are linearly
     independent over GF(q); then Sherwood's construction gives CA(n(q^t - 1) + 1; t, k, q).
  A  rank of the t x t matrix by Gaussian elimination.
  B  written separately: no nonzero combination sum c_i x_i (first nonzero c_i = 1, so (q^t - 1)/(q - 1) of them)
     vanishes; its field tables come from schoolbook products of the digit polynomials (A's from shift-and-add).
scphf: the Sherwood families of github.com/arhancanli/quantavane-covering-arrays: each integer is a point of
     AG(t-1, q) (t-1 digits base q, most significant first); every t columns must have a row whose points span
     AG(t-1, q) (the vectors (1, point) independent); the CA has n(q^t - q) + q rows. Same routines A and B.
--rows N (all modes) states the claimed number of rows.
--cols K checks only the first K columns (deleting columns keeps either property). Fields: q prime, or q in the
built-in Conway table (4, 8, 9, 16, 25, 27, 32, 49, 64), or any q = p^e with --modulus giving the low coefficients
of the monic polynomial. Exit 0 = both PASS. Pure Python: practical up to a few million t-sets; the large families
in github.com/arhancanli/quantavane-covering-arrays have C checkers in its checker/ folder.
MIT licence (LICENSE-CODE).
"""
import itertools
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


# ------------------------------------------------------------------------------------------- covering arrays
def read_ca(path):
    rows = [list(map(int, l.split())) for l in open(path) if l.split() and l.split()[0].isdigit()]
    if not rows or len({len(r) for r in rows}) != 1:
        raise ValueError("rows of unequal length")
    return rows


def ca_a(rows, t, v, K):
    cols = list(zip(*rows))[:K]
    need = v ** t
    for S in itertools.combinations(range(K), t):
        if len(set(zip(*(cols[c] for c in S)))) != need:
            return False, "columns %s miss a t-tuple" % (S,)
    return True, "every one of the %d t-sets of columns covers all %d tuples" % (_comb(K, t), need)


def ca_b(rows, t, v, K):
    masks = [[0] * v for _ in range(K)]
    for r, row in enumerate(rows):
        for c in range(K):
            masks[c][row[c]] |= 1 << r
    for S in itertools.combinations(range(K), t):
        for tup in itertools.product(range(v), repeat=t):
            m = masks[S[0]][tup[0]]
            for c, x in zip(S[1:], tup[1:]):
                m &= masks[c][x]
                if not m:
                    break
            if not m:
                return False, "columns %s miss %s" % (S, tup)
    return True, "every t-tuple found in every t-set of columns by row masks"


def _comb(n, r):
    out = 1
    for i in range(r):
        out = out * (n - i) // (i + 1)
    return out


# ------------------------------------------------------------------------------------------- CPHF
def read_cphf(path, t, q, affine=False):
    rows = [list(map(int, l.split())) for l in open(path) if l.split() and l.split()[0].isdigit()]
    if not rows or len({len(r) for r in rows}) != 1:
        raise ValueError("rows of unequal length")

    def vec(x):
        if affine:                                        # a point of AG(t-1, q): the vector (1, digits)
            if not 0 <= x < q ** (t - 1):
                raise ValueError("entry %d is not a point of AG(%d, %d)" % (x, t - 1, q))
            return [1] + [(x // q ** (t - 2 - i)) % q for i in range(t - 1)]
        if not 0 <= x < q ** t:
            raise ValueError("entry %d is not a vector of GF(%d)^%d" % (x, q, t))
        return [(x // q ** (t - 1 - i)) % q for i in range(t)]
    return [[vec(x) for x in r] for r in rows]


def cphf_a(fam, t, F, K):
    add, mul, p = F
    q = len(add)
    inv = [0] * q
    for a in range(1, q):
        inv[a] = next(b for b in range(1, q) if mul[a][b] == 1)
    neg = [next(b for b in range(q) if add[a][b] == 0) for a in range(q)]

    def rank(M):
        M = [list(r) for r in M]
        rk = 0
        for c in range(t):
            piv = next((i for i in range(rk, t) if M[i][c]), None)
            if piv is None:
                continue
            M[rk], M[piv] = M[piv], M[rk]
            ip = inv[M[rk][c]]
            M[rk] = [mul[ip][x] for x in M[rk]]
            for i in range(t):
                if i != rk and M[i][c]:
                    co = neg[M[i][c]]
                    M[i] = [add[x][mul[co][y]] for x, y in zip(M[i], M[rk])]
            rk += 1
        return rk
    for S in itertools.combinations(range(K), t):
        if not any(rank([row[c] for c in S]) == t for row in fam):
            return False, "columns %s have no row with independent vectors" % (S,)
    return True, "all %d t-sets of columns have an independent row (Gaussian elimination)" % _comb(K, t)


def cphf_b(fam, t, mul, K):
    q = len(mul)
    p = next(p for p in range(2, q + 1) if q % p == 0)
    e = 1
    while p ** e < q:
        e += 1

    def addf(a, b):                                       # digit-wise addition mod p
        return sum((((a // p ** i) + (b // p ** i)) % p) * p ** i for i in range(e))
    add = [[addf(a, b) for b in range(q)] for a in range(q)]
    # combinations with first nonzero coefficient 1, one per projective point of GF(q)^t
    combos = []
    for c in itertools.product(range(q), repeat=t):
        nz = next((x for x in c if x), None)
        if nz == 1:
            combos.append(c)

    def independent(vs):
        for c in combos:
            s = [0] * t
            for ci, v in zip(c, vs):
                if ci:
                    s = [add[a][mul[ci][b]] for a, b in zip(s, v)]
            if not any(s):
                return False
        return True
    for S in itertools.combinations(range(K), t):
        if not any(independent([row[c] for c in S]) for row in fam):
            return False, "columns %s: every row has a vanishing combination" % (S,)
    return True, "no t-set of columns has a dependent combination in every row (%d combinations each)" % len(combos)


def main(argv):
    def opt(name, conv=int):
        return conv(argv[argv.index(name) + 1]) if name in argv else None
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    try:
        mode, path, t, v = args[0], args[1], int(args[2]), int(args[3])
        K = opt("--cols")
        if mode == "ca":
            rows = read_ca(path)
            k = len(rows[0])
            if any(x < 0 or x >= v for r in rows for x in r):
                raise ValueError("a symbol outside 0..%d" % (v - 1))
            K = min(K or k, k)
            want = opt("--rows")
            print("CA with %d rows, %d columns (checking %d), t = %d, v = %d" % (len(rows), k, K, t, v))
            if want is not None and want != len(rows):
                print("FAIL  %d rows, %d claimed" % (len(rows), want))
                return 1
            runs = (("A", lambda: ca_a(rows, t, v, K)), ("B", lambda: ca_b(rows, t, v, K)))
        elif mode in ("cphf", "scphf"):
            mod = opt("--modulus", lambda s: [int(x) for x in s.split(",")])
            F = field(v, mod)
            fam = read_cphf(path, t, v, affine=mode == "scphf")
            k = len(fam[0])
            K = min(K or k, k)
            N = len(fam) * (v ** t - v) + v if mode == "scphf" else len(fam) * (v ** t - 1) + 1
            print("%s(%d; %d, %d, %d) (checking %d columns) -> CA(%d; %d, %d, %d)"
                  % (mode.upper(), len(fam), k, v, t, K, N, t, K, v))
            want = opt("--rows")
            if want is not None and want != N:
                print("FAIL  gives %d rows, %d claimed" % (N, want))
                return 1
            Fb = field_b(v, mod)
            runs = (("A", lambda: cphf_a(fam, t, F, K)), ("B", lambda: cphf_b(fam, t, Fb, K)))
        else:
            raise ValueError("mode must be ca, cphf or scphf")
    except (IndexError, ValueError, AssertionError) as e:
        print("bad input: %s\n%s" % (e, __doc__))
        return 2
    ok = True
    for name, fn in runs:
        res, why = fn()
        print("%s %s  %s" % (name, "PASS" if res else "FAIL", why))
        ok &= res
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
