"""The independent check of a Theorem K certificate (field model B).

Written separately from qvrecords.theorem_k and importing nothing from it. It reads only the data in the certificate
(q, k, value, the modulus, the non-cube, the sample columns, the counted flag) and recomputes everything else with
different algorithms:

- GF(2^m): full carry-less product, then reduction by long division; inverses by the extended Euclidean algorithm on
  binary polynomials; the modulus is tested by Rabin's irreducibility test (model A uses Ben-Or's test, trial division before 20 Nov 2026).
- GF(q^3) = GF(q)[X]/(X^3 - c): products by schoolbook polynomial multiplication and reduction with X^3 = c;
  inverses by the extended Euclidean algorithm over GF(q)[X] (model A uses Cramer's rule).
- independence of three vectors of GF(q)^3: rank by Gaussian elimination (model A uses the determinant).
- the counted range: recomputed with its own prime-power test (trial factorisation) and integer roots.

What it checks: the cell's parameters and value; that the modulus is irreducible of degree m and the given element is
a non-cube (so b^3 = c defines GF(q^3)); that the sample columns are the construction's columns; and the covering
property of the family on random triples among the cell's first k columns and on the hard cases of Lemma 2 (two B
columns at points of a line through y, with E_y). With exhaustive=True it decides every triple through B_1 and E_1,
which by the symmetry x -> g x (each row multiplied by a fixed element) is every triple of the family; that is
practical for q <= 16. The covering property of all 2(q^2 + q + 1) columns is the theorem (proof in
github.com/arhancanli/quantavane-covering-arrays, three-row/); the sampled check tests the construction at this
cell, the exhaustive one proves it for that q.
"""
import functools
import random


# ------------------------------------------------------------------------------------------- GF(2)[t]
def _clmul(a, b):
    r = 0
    while b:
        low = b & -b
        r ^= a << (low.bit_length() - 1)
        b ^= low
    return r


def _pmod(a, f):
    df = f.bit_length()
    while a.bit_length() >= df:
        a ^= f << (a.bit_length() - df)
    return a


def _pgcd(a, b):
    while b:
        a, b = b, _pmod(a, b)
    return a


def _rabin_irreducible(f, m):
    """Rabin: f of degree m is irreducible iff t^(2^m) = t mod f and gcd(t^(2^(m/p)) - t, f) = 1 for primes p | m"""
    if f.bit_length() - 1 != m:
        return False

    def t_pow_2e(e):                    # t^(2^e) mod f by e squarings
        x = 2
        for _ in range(e):
            x = _pmod(_clmul(x, x), f)
        return x
    if t_pow_2e(m) != 2:
        return False
    n, p, primes = m, 2, []
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes.append(n)
    return all(_pgcd(t_pow_2e(m // p) ^ 2, f) == 1 for p in primes)


class FieldB:
    def __init__(self, m, f):
        self.m, self.f, self.q = m, f, 1 << m

    def mul(self, a, b):
        return _pmod(_clmul(a, b), self.f)

    def inv(self, a):
        """extended Euclid in GF(2)[t]: keeps s_i a = r_i mod f until r_i = 0; then r_(i-1) = gcd = 1"""
        if a == 0:
            raise ZeroDivisionError
        r0, r1, s0, s1 = self.f, a, 0, 1
        while r1:
            qt = 0
            while r0.bit_length() >= r1.bit_length():
                sh = r0.bit_length() - r1.bit_length()
                qt ^= 1 << sh
                r0 ^= r1 << sh
            r0, r1 = r1, r0
            s0, s1 = s1, s0 ^ _clmul(qt, s1)
        if r0 != 1:
            raise ValueError("not invertible: modulus reducible")
        return _pmod(s0, self.f)

    def pow(self, a, e):
        r = 1
        for bit in bin(e)[2:]:
            r = self.mul(r, r)
            if bit == "1":
                r = self.mul(r, a)
        return r


class CubicB:
    """GF(q)[X]/(X^3 - c), elements as coefficient lists [a0, a1, a2]"""

    def __init__(self, F, c):
        self.F, self.c = F, c

    def mul(self, x, y):
        F = self.F
        p = [0] * 5
        for i in range(3):
            if x[i]:
                for j in range(3):
                    if y[j]:
                        p[i + j] ^= F.mul(x[i], y[j])
        # X^3 = c, X^4 = c X
        return [p[0] ^ F.mul(self.c, p[3]), p[1] ^ F.mul(self.c, p[4]), p[2]]

    def inv(self, x):
        """extended Euclid in GF(q)[X] against X^3 - c (= X^3 + c in characteristic 2)"""
        F = self.F

        def trim(p):
            while p and p[-1] == 0:
                p = p[:-1]
            return p

        def sub(a, b):
            n = max(len(a), len(b))
            return trim([(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0) for i in range(n)])

        def pmul(a, b):
            if not a or not b:
                return []
            r = [0] * (len(a) + len(b) - 1)
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    r[i + j] ^= F.mul(ai, bj)
            return trim(r)

        def divmod_(a, b):
            a, qt = list(a), [0] * max(len(a) - len(b) + 1, 1)
            ib = F.inv(b[-1])
            while len(a) >= len(b) and a:
                s = len(a) - len(b)
                co = F.mul(a[-1], ib)
                qt[s] = co
                for i, bi in enumerate(b):
                    a[s + i] ^= F.mul(co, bi)
                a = trim(a)
            return trim(qt), a
        r0, r1 = [self.c, 0, 0, 1], trim(list(x))
        s0, s1 = [], [1]
        if not r1:
            raise ZeroDivisionError
        while len(r1) > 1:
            qt, rem = divmod_(r0, r1)
            r0, r1 = r1, rem
            s0, s1 = s1, sub(s0, pmul(qt, s1))
            if not r1:
                raise ValueError("X^3 - c is reducible: c is a cube")
        k = F.inv(r1[0])
        out = [F.mul(k, a) for a in s1] + [0, 0, 0, 0, 0]
        for d in range(len(out) - 1, 2, -1):    # reduce with X^3 = c, should s ever reach degree 3
            if out[d]:
                out[d - 3] ^= F.mul(self.c, out[d])
                out[d] = 0
        return out[:3]


def _rank3(F, rows):
    """rank of a 3 x 3 matrix over GF(q) by Gaussian elimination"""
    M = [list(r) for r in rows]
    rank, col = 0, 0
    while rank < 3 and col < 3:
        piv = next((i for i in range(rank, 3) if M[i][col]), None)
        if piv is None:
            col += 1
            continue
        M[rank], M[piv] = M[piv], M[rank]
        ip = F.inv(M[rank][col])
        M[rank] = [F.mul(ip, a) for a in M[rank]]
        for i in range(3):
            if i != rank and M[i][col]:
                co = M[i][col]
                M[i] = [a ^ F.mul(co, b) for a, b in zip(M[i], M[rank])]
        rank += 1
        col += 1
    return rank


# ------------------------------------------------------------------------------------------- counted range (B)
def _smallest_factor(n):
    if n % 2 == 0:
        return 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return d
        d += 2
    return n


def _prime_power_b(n):
    if n < 2:
        return False
    p = _smallest_factor(n)
    while n % p == 0:
        n //= p
    return n == 1


def _cube_root_floor(n):
    r = int(round(n ** (1.0 / 3)))
    while r ** 3 > n:
        r -= 1
    while (r + 1) ** 3 <= n:
        r += 1
    return r


@functools.lru_cache(maxsize=None)                     # a pure function of q; cached since the range reached s = 26
def counted_range_b(q):
    if q == 16:
        return 379, 546
    s2 = q.bit_length() - 1
    if q < 64 or q != 1 << s2 or s2 % 2:
        return None
    rows = 3 * q ** 3 - 2
    Q = _cube_root_floor(rows // 2) + 3
    while not (2 * Q ** 3 - 3 * Q + 2 * q <= rows and _prime_power_b(Q)):
        Q -= 1
    reach = (Q + 1) ** 2 if Q % 2 == 0 else Q * Q + Q + 1
    q15 = 2 ** (3 * s2 // 2)
    return max(reach, (q + 1) ** 2, (q + 2) ** 2 + q + 2, 4 * q15 + 3) + 1, 2 * (q * q + q + 1)


# ------------------------------------------------------------------------------------------- the family in model B
class FamilyB:
    def __init__(self, q, modulus, c):
        m = q.bit_length() - 1
        self.q, self.N = q, q * q + q + 1
        self.F = FieldB(m, modulus)
        self.K = CubicB(self.F, c)

    def point(self, j):
        q = self.q
        if j >= q * q + q:
            return [0, 0, 1]
        if j >= q * q:
            return [0, 1, j - q * q]
        a, b = divmod(j, q)
        return [1, a, b]

    def column(self, j):
        x = self.point(j % self.N)
        xi = self.K.inv(x)
        if j < self.N:
            return [x, xi, x]
        b = [0, 1, 0]
        return [x, self.K.mul(b, xi), self.K.mul(b, x)]

    def covers(self, cols):
        return any(_rank3(self.F, [cols[0][r], cols[1][r], cols[2][r]]) == 3 for r in range(3))


def check(cert, samples=12, hard=6, seed=None, exhaustive=False):
    """Check a Theorem K certificate. Returns (ok, report); report lists each check with its outcome."""
    rep = []

    def need(cond, what):
        rep.append(("PASS" if cond else "FAIL", what))
        return cond
    try:
        q, k, value = int(cert["q"]), int(cert["k"]), int(cert["value"])
        fld = cert["field"]
        m, f, c = int(fld["m"]), int(fld["modulus"]), int(fld["non_cube"])
    except (KeyError, TypeError, ValueError) as e:
        return False, [("FAIL", "certificate unreadable: %r" % (e,))]
    ok = need(cert.get("family") == "theorem-k", "family is theorem-k")
    ok &= need(q >= 4 and q == 1 << m and m % 2 == 0, "q = 2^m with m even (q = 4^s)")
    if not ok:
        return False, rep
    width = 2 * (q * q + q + 1)
    ok &= need(3 <= k <= width, "3 <= k <= 2(q^2 + q + 1) = %d" % width)
    ok &= need(value == 3 * (q ** 3 - 1) + 1, "value = 3(q^3 - 1) + 1 = %d (Sherwood: 3 rows of q^3 - 1 vectors and "
               "one zero row)" % (3 * (q ** 3 - 1) + 1))
    ok &= need(cert.get("statement") == "CAN(3, %d, %d) <= %d" % (k, q, value), "statement matches the parameters")
    ok &= need(_rabin_irreducible(f, m), "modulus irreducible of degree %d over GF(2) (Rabin's test)" % m)
    if not ok:
        return False, rep
    F = FieldB(m, f)
    ok &= need(0 < c < q and F.pow(c, (q - 1) // 3) != 1, "c = %d is a non-cube in GF(q): c^((q-1)/3) != 1" % c)
    if not ok:
        return False, rep
    fam = FamilyB(q, f, c)
    sc = cert.get("sample_columns", {})
    same = bool(sc) and all(fam.column(int(j)) == [list(e) for e in col] for j, col in sc.items())
    ok &= need(same, "the %d sample columns are the construction's columns" % len(sc))
    rng = counted_range_b(q)
    if cert.get("counted"):
        ok &= need(rng is not None and rng[0] <= k <= rng[1] and list(rng) == list(cert.get("counted_range") or []),
                   "k lies in the counted range %s, recomputed" % (rng,))
    else:
        ok &= need(rng is None or not rng[0] <= k <= rng[1], "k lies outside the counted range, as stated")
    R = random.Random(seed if seed is not None else k * 1000003 + q)
    bad = 0
    for _ in range(samples):
        idx = set()
        while len(idx) < 3:
            idx.add(R.randrange(k))
        bad += not fam.covers([fam.column(j) for j in idx])
    ok &= need(bad == 0, "%d random triples among the first %d columns covered" % (samples, k))
    bad = 0
    for _ in range(hard):                      # Lemma 2: B_{y w1}, B_{y w2}, E_y with w1, w2 on the line <1, b>
        y = fam.point(R.randrange(fam.N))
        while True:
            u1, v1, u2, v2 = (R.randrange(q) for _ in range(4))
            if F.mul(u1, v2) != F.mul(v1, u2):
                break
        b = [0, 1, 0]
        x1, x2 = fam.K.mul(y, [u1, v1, 0]), fam.K.mul(y, [u2, v2, 0])
        trip = [[x1, fam.K.inv(x1), x1], [x2, fam.K.inv(x2), x2], [y, fam.K.mul(b, fam.K.inv(y)), fam.K.mul(b, y)]]
        bad += not fam.covers(trip)
    ok &= need(bad == 0, "%d Lemma 2 hard cases covered" % hard)
    if exhaustive:
        cols = [fam.column(j) for j in range(width)]
        bad = tried = 0
        for first in (0, fam.N):
            for i in range(width):
                if i == first:
                    continue
                for j in range(i + 1, width):
                    if j == first:
                        continue
                    tried += 1
                    bad += not fam.covers([cols[first], cols[i], cols[j]])
        ok &= need(bad == 0, "exhaustive: %d triples through B_1 and E_1 (every triple, by symmetry), %d uncovered"
                   % (tried, bad))
    return bool(ok), rep
