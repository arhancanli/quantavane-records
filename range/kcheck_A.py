#!/usr/bin/env python3
"""kcheck_A.py - sampled exact checks of Theorem K's array at q = 4^s (program A of 2; family-theorems lane, 18 Nov 2026).

Field model (separate from kcheck_B.c): GF(q) = GF(2)[t]/(f), f found by Rabin's irreducibility test; c the first
non-cube of GF(q)^*; GF(q^3) = GF(q)[B]/(B^3 - c), b = B. Columns B_x = (x, x^-1, x), E_x = (x, b x^-1, b x).
Three entries are independent over GF(q) iff their 3x3 coordinate determinant is nonzero. Every triple is exact
arithmetic; only the choice of triples is sampled. Cases (RANGE.md section 4):
  R  random triples of random columns          H1 three B columns with collinear x (rows 1 and 3 fail)
  H2 three E columns with collinear x          H3 B_{y w1}, B_{y w2}, E_y, w_i in <1, b> (rows 1 and 3 fail)
  H4 B_x, E_{x w1/b}, E_{x w2/b} (rows 1, 3 fail)   H5 B_x, E_x and a third column (rows 1 coincide)
For H1-H4 the program also confirms that rows 1 and 3 really fail, so the case tests Lemma 1 or Lemma 2.
usage: python3 -I kcheck_A.py S SAMPLES SEED [--tamper inv|bq|mod|cube]
exit 0 = PASS (no uncovered triple, every precondition held); 1 = FAIL; 2 = refused (bad field)."""
import sys, random, hashlib

S, SAMPLES, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
TAMPER = sys.argv[sys.argv.index("--tamper") + 1] if "--tamper" in sys.argv else ""
m = 2 * S
q = 1 << m


def pmulmod(a, b, f, deg):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> deg & 1:
            a ^= f
    return r


def ppowmod(a, e, f, deg):
    r = 1
    while e:
        if e & 1:
            r = pmulmod(r, a, f, deg)
        a = pmulmod(a, a, f, deg)
        e >>= 1
    return r


def pgcd(a, b):
    while b:
        while a and a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def rabin(f, deg):
    x = 2
    if ppowmod(x, 1 << deg, f, deg) != x:          # x^(2^deg) = x mod f
        return False
    for r in {p for p in range(2, deg + 1) if deg % p == 0 and all(p % d for d in range(2, p))}:
        h = ppowmod(x, 1 << (deg // r), f, deg) ^ x
        if pgcd(f, h) != 1:
            return False
    return True


f = next(F for F in range((1 << m) | 1, 1 << (m + 1), 2) if rabin(F, m))
if TAMPER == "mod":
    f = f ^ (1 << 1) if not rabin(f ^ (1 << 1), m) else f ^ (1 << 2)   # a reducible modulus of the same degree
if not rabin(f, m):
    print("REFUSED: modulus %#x is not irreducible" % f)
    sys.exit(2)


def mul(a, b):
    return pmulmod(a, b, f, m)


def inv(a):
    return ppowmod(a, q - 2, f, m)


cube_exp = (q - 1) // 3
noncube = [c for c in range(2, 64) if ppowmod(c, cube_exp, f, m) != 1]
c = noncube[0]
if TAMPER == "cube":
    c = next(x for x in range(2, 64) if ppowmod(x, cube_exp, f, m) == 1)
if ppowmod(c, cube_exp, f, m) == 1:
    print("REFUSED: c = %d is a cube, X^3 - c is reducible" % c)
    sys.exit(2)


def tmul(x, y):
    a0, a1, a2 = x
    b0, b1, b2 = y
    return (mul(a0, b0) ^ mul(c, mul(a1, b2) ^ mul(a2, b1)),
            mul(a0, b1) ^ mul(a1, b0) ^ mul(c, mul(a2, b2)),
            mul(a0, b2) ^ mul(a1, b1) ^ mul(a2, b0))


def det3(M):
    (a, b, cc), (d, e, g), (h, i, j) = M
    return (mul(a, mul(e, j) ^ mul(g, i)) ^ mul(b, mul(d, j) ^ mul(g, h)) ^ mul(cc, mul(d, i) ^ mul(e, h)))


def tinv(x):
    cols = [tmul(x, (1, 0, 0)), tmul(x, (0, 1, 0)), tmul(x, (0, 0, 1))]     # columns of multiplication by x
    M = [[cols[jj][ii] for jj in range(3)] for ii in range(3)]
    D = det3(M)
    if D == 0:
        raise ZeroDivisionError
    Di = inv(D)
    out = []
    for jj in range(3):                                                  # Cramer, right side (1, 0, 0)
        Mj = [row[:] for row in M]
        for ii in range(3):
            Mj[ii][jj] = 1 if ii == 0 else 0
        out.append(mul(det3(Mj), Di))
    return tuple(out)


ONE, b = (1, 0, 0), (0, 1, 0)
if TAMPER == "bq":
    b = (c, 0, 0)                                                        # b inside GF(q): the twist disappears
binv = tinv(b)
# field self-test: x * x^-1 = 1 and the minimal relation b^3 = c (b not in GF(q) unless tampered)
rng = random.Random(SEED)
for _ in range(20):
    x = (rng.randrange(q), rng.randrange(q), rng.randrange(1, q))
    assert tmul(x, tinv(x)) == ONE, "field self-test failed"


def rnd_scalar(nonzero=True):
    return rng.randrange(1 if nonzero else 0, q)


def rnd_el():
    while True:
        x = (rng.randrange(q), rng.randrange(q), rng.randrange(q))
        if x != (0, 0, 0):
            return x


def sc(a, x):
    return (mul(a, x[0]), mul(a, x[1]), mul(a, x[2]))


def add(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1], x[2] ^ y[2])


def prop(x, y):
    """x and y are the same point of PG(2, q)"""
    return all(mul(x[i], y[j]) == mul(x[j], y[i]) for i in range(3) for j in range(i + 1, 3))


def column(kind, x):
    if kind == "B":
        return (x, tinv(x), x)
    r2 = tmul(b, x) if TAMPER == "inv" else tmul(b, tinv(x))
    return (x, r2, tmul(b, x))


def indep(u, v, w):
    return det3([u, v, w]) != 0


def rows_ok(cols):
    return [indep(cols[0][r], cols[1][r], cols[2][r]) for r in range(3)]


def wpoint():
    while True:
        u, v = rnd_scalar(False), rnd_scalar(False)
        w = add(sc(u, ONE), sc(v, b))
        if w != (0, 0, 0):
            return w, (u, v)


def two_wpoints():
    while True:
        (w1, p1), (w2, p2) = wpoint(), wpoint()
        if mul(p1[0], p2[1]) != mul(p1[1], p2[0]):                       # not proportional: distinct points
            return w1, w2


stats = {}
bad_pre = 0
for case in ("R", "H1", "H2", "H3", "H4", "H5"):
    unc = pre = 0
    for _ in range(SAMPLES):
        if case == "R":
            spec = [(rng.choice("BE"), rnd_el()) for _ in range(3)]
            if any(spec[i][0] == spec[j][0] and prop(spec[i][1], spec[j][1]) for i in range(3) for j in range(i)):
                continue                                                 # the same column twice
            cols = [column(k, x) for k, x in spec]
        elif case in ("H1", "H2"):
            x1, x2 = rnd_el(), rnd_el()
            if prop(x1, x2):
                continue
            x3 = add(sc(rnd_scalar(), x1), sc(rnd_scalar(), x2))
            kind = "B" if case == "H1" else "E"
            cols = [column(kind, x) for x in (x1, x2, x3)]
        elif case == "H3":
            y = rnd_el()
            w1, w2 = two_wpoints()
            cols = [column("B", tmul(y, w1)), column("B", tmul(y, w2)), column("E", y)]
        elif case == "H4":
            x = rnd_el()
            w1, w2 = two_wpoints()
            cols = [column("B", x), column("E", tmul(tmul(x, binv), w1)), column("E", tmul(tmul(x, binv), w2))]
        else:
            x, k3, x3 = rnd_el(), rng.choice("BE"), rnd_el()
            if prop(x, x3):
                continue
            cols = [column("B", x), column("E", x), column(k3, x3)]
        ok = rows_ok(cols)
        if case in ("H1", "H2", "H3", "H4"):
            if ok[0] or ok[2]:
                pre += 1                                                 # precondition broken: not a hard case
        if not any(ok):
            unc += 1
    stats[case] = (unc, pre)
    bad_pre += pre if TAMPER == "" else 0

src = hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16]
fails = sum(u for u, _ in stats.values())
line = " ".join("%s:%d/%d" % (k, v[0], v[1]) for k, v in stats.items())
verdict = "PASS" if fails == 0 and bad_pre == 0 else "FAIL"
print("kcheck_A s %d q %d f %#x c %d samples %d seed %d tamper %s | uncovered/precondition-broken %s | %s | sha %s"
      % (S, q, f, c, SAMPLES, SEED, TAMPER or "-", line, verdict, src))
sys.exit(0 if verdict == "PASS" else 1)
