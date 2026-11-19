"""Theorem K: the construction routine (field model A).

For q = 4^s, let c be a non-cube of GF(q) and b a root of X^3 - c in GF(q^3). Points of PG(2, q) are the elements x
of GF(q^3)^* / GF(q)^*, with GF(q^3) read as GF(q)^3 in the basis (1, b, b^2). The 3 x 2(q^2+q+1) array with columns

    B_x = (x, x^-1, x)      and      E_x = (x, b x^-1, b x)

is a covering perfect hash family CPHF(3; 2(q^2+q+1), q, 3): any three columns have a row whose three entries are
linearly independent over GF(q). By Sherwood, Martirosyan and Colbourn (J. Combin. Des. 14 (2006) 202-213) it gives a
covering array CA(3q^3 - 2; 3, k, q) for every k <= 2(q^2+q+1), so CAN(3, k, q) <= 3q^3 - 2.

The construction and a sketch of its proof are public in github.com/arhancanli/quantavane-covering-arrays, folder
three-row/. This module builds the certificate of any cell; qvrecords.verify checks it with separately written code.

Field model A: GF(2^m) as bit vectors modulo the first irreducible polynomial of degree m found by trial division
(for m = 4 that is x^4 + x + 1, the field of the ca-tools files); products by shift-and-add with reduction at every
step; inverses by Fermat (a^(q-2)); GF(q^3) inverses by Cramer's rule on the 3 x 3 multiplication matrix.
"""

THEOREM = "Theorem K"
SOURCE = "https://github.com/arhancanli/quantavane-covering-arrays/tree/main/three-row"


def _irreducible_trial(f, m):
    """True when the binary polynomial f of degree m has no factor of degree 1 .. m // 2 (trial division)."""
    for d in range(1, m // 2 + 1):
        for g in range(1 << d, 1 << (d + 1)):
            a = f
            while a.bit_length() >= g.bit_length():
                a ^= g << (a.bit_length() - g.bit_length())
            if a == 0:
                return False
    return True


class FieldA:
    """GF(2^m) with the first irreducible modulus found by trial division."""

    def __init__(self, m):
        self.m = m
        self.q = 1 << m
        self.modulus = next(f for f in range((1 << m) | 1, 1 << (m + 1), 2) if _irreducible_trial(f, m))

    def mul(self, a, b):
        m, f, r = self.m, self.modulus, 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> m & 1:
                a ^= f
        return r

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def inv(self, a):
        assert a, "zero has no inverse"
        return self.pow(a, self.q - 2)

    def first_non_cube(self):
        e = (self.q - 1) // 3
        return next(x for x in range(2, self.q) if self.pow(x, e) != 1)


class TheoremK:
    """The family at one q = 4^s: columns, points and the covering test, in field model A."""

    def __init__(self, q):
        if q < 4 or q & (q - 1) or (q.bit_length() - 1) % 2:
            raise ValueError("Theorem K needs q = 4^s with s >= 1; got %r" % (q,))
        self.q = q
        self.F = FieldA(q.bit_length() - 1)
        self.c = self.F.first_non_cube()
        self.N = q * q + q + 1                     # points of PG(2, q)
        self.width = 2 * self.N                    # columns of the CPHF

    # GF(q^3) = GF(q)[b]/(b^3 - c); elements are triples (a0, a1, a2) = a0 + a1 b + a2 b^2
    def tmul(self, x, y):
        m, c = self.F.mul, self.c
        return (m(x[0], y[0]) ^ m(c, m(x[1], y[2]) ^ m(x[2], y[1])),
                m(x[0], y[1]) ^ m(x[1], y[0]) ^ m(c, m(x[2], y[2])),
                m(x[0], y[2]) ^ m(x[1], y[1]) ^ m(x[2], y[0]))

    def det(self, u, v, w):
        m = self.F.mul
        return (m(u[0], m(v[1], w[2]) ^ m(v[2], w[1])) ^ m(u[1], m(v[0], w[2]) ^ m(v[2], w[0]))
                ^ m(u[2], m(v[0], w[1]) ^ m(v[1], w[0])))

    def tinv(self, x):
        """inverse in GF(q^3) by Cramer's rule: solve M y = e_0, M the matrix of multiplication by x"""
        cols = [self.tmul(x, e) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
        D = self.det(*cols)
        Di = self.F.inv(D)
        out = []
        for j in range(3):
            cj = list(cols)
            cj[j] = (1, 0, 0)
            out.append(self.F.mul(self.det(*cj), Di))
        return tuple(out)

    def point(self, j):
        """canonical order of PG(2, q): (1, a, b) for j < q^2, then (0, 1, b), then (0, 0, 1)"""
        q = self.q
        if j < q * q:
            return (1, j // q, j % q)
        if j < q * q + q:
            return (0, 1, j - q * q)
        return (0, 0, 1)

    def column(self, j):
        """column j: B_x for j < N, E_x for N <= j < 2N, x = point(j mod N)"""
        if not 0 <= j < self.width:
            raise IndexError("column %d outside 0..%d" % (j, self.width - 1))
        x = self.point(j % self.N)
        xi = self.tinv(x)
        if j < self.N:
            return (x, xi, x)
        b = (0, 1, 0)
        return (x, self.tmul(b, xi), self.tmul(b, x))

    def covered(self, cols):
        """the rows in which the three columns' entries are linearly independent"""
        return [r for r in range(3) if self.det(cols[0][r], cols[1][r], cols[2][r])]


# ------------------------------------------------------------------------------------------- the counted range
def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):   # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _is_prime_power(n):
    if n < 2:
        return False
    for e in range(n.bit_length(), 0, -1):
        r = _iroot(n, e)
        if r ** e == n and _is_prime(r):
            return True
    return False


def _iroot(n, e):
    lo, hi = 0, 1 << (n.bit_length() // e + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** e <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def counted_range(q):
    """The widths k at which the lab counts the cell CAN(3, k, q) <= 3q^3 - 2 as a record, as (k_from, k_to).

    q = 16: 379 <= k <= 546, the cells the keeper's table (ca-tools, PR #83, 17 November 2026) lists from this family.
    q = 4^s >= 64: k > R(q) up to 2(q^2+q+1), where R(q) is the largest width that any published construction class
    with at most 3q^3 - 2 rows reaches, by the lab's range analysis (refereed inside the lab, 18 November 2026; not yet
    peer reviewed). R(q) is the largest of: the class-I reach Q^2 + Q + 1 (Q odd) or (Q + 1)^2 (Q even), Q the largest
    prime power with 2Q^3 - 3Q + 2q <= 3q^3 - 2; (q + 1)^2; (q + 2)^2 + q + 2; and 4 q^(3/2) + 3.
    """
    if q == 16:
        return 379, 546
    if q < 64 or q & (q - 1) or (q.bit_length() - 1) % 2:
        raise ValueError("counted range defined for q = 16 and q = 4^s >= 64")
    NK = 3 * q ** 3 - 2
    Q = _iroot(NK // 2, 3) + 2
    while 2 * Q ** 3 - 3 * Q + 2 * q > NK or not _is_prime_power(Q):
        Q -= 1
    reach_I = Q * Q + Q + 1 if Q % 2 else (Q + 1) ** 2
    q32 = 1 << (3 * (q.bit_length() - 1) // 2)            # q^(3/2), exact for q = 4^s
    R = max(reach_I, (q + 1) ** 2, (q + 2) ** 2 + q + 2, 4 * q32 + 3)
    return R + 1, 2 * (q * q + q + 1)


# ------------------------------------------------------------------------------------------- certificates
_CACHE = {}


def family(q):
    if q not in _CACHE:
        if len(_CACHE) > 8:
            _CACHE.clear()
        _CACHE[q] = TheoremK(q)
    return _CACHE[q]


def certificate(q, k, show_columns=4):
    """value and certificate of the cell CAN(3, k, q) <= 3q^3 - 2"""
    K = family(q)
    if not 1 <= k <= K.width:
        raise ValueError("k must be in 1..%d for q = %d" % (K.width, q))
    value = 3 * q ** 3 - 2
    try:
        lo, hi = counted_range(q)
        counted = lo <= k <= hi
        rng = [lo, hi]
    except ValueError:
        counted, rng = False, None
    idx = sorted(set(list(range(min(k, show_columns))) + [k - 1]))
    cert = {
        "family": "theorem-k",
        "theorem": THEOREM,
        "statement": "CAN(3, %d, %d) <= %d" % (k, q, value),
        "t": 3, "k": k, "q": q, "value": value,
        "object": "the first %d columns of CPHF(3; %d, %d, 3), expanded by Sherwood's construction to CA(%d; 3, %d, %d)"
                  % (k, K.width, q, value, k, q),
        "field": {"m": K.F.m, "modulus": K.F.modulus, "modulus_bits": bin(K.F.modulus)[2:], "non_cube": K.c,
                  "extension": "GF(q^3) = GF(q)[b]/(b^3 - non_cube), coordinates in the basis (1, b, b^2)"},
        "columns_rule": "j < N: B_x = (x, x^-1, x); N <= j < 2N: E_x = (x, b x^-1, b x); x = point(j mod N); "
                        "N = q^2 + q + 1; points in the order (1, a0, a1) for j = a0 q + a1 < q^2, then (0, 1, j - q^2), "
                        "then (0, 0, 1); field elements as integers, bit i = coefficient of t^i",
        "sample_columns": {str(j): [list(e) for e in K.column(j)] for j in idx},
        "counted": counted,
        "counted_range": rng,
        "record_basis": ("q = 16: listed by the ca-tools table (keeper merged PR #83, 17 November 2026)" if q == 16 else
                         "past every published table (none lists strength 3 at q >= 64); k above the reach R(q) of "
                         "every published construction class with at most 3q^3 - 2 rows, by the lab's range analysis "
                         "(refereed inside the lab; not yet peer reviewed)") if counted else
                        "an upper bound from the construction; not counted as a record at this width",
        "source": SOURCE,
    }
    if q == 16 and k <= 468:
        cert["note"] = ("for k <= 468 the Sherwood restriction CA(12256; 3, 468, 16) of the same family is smaller "
                        "(three-row/families/CA_12256_3_468_16.cphf.txt); this cell is counted once")
    return value, cert


def ca_entry(q, k, row, col):
    """entry (row, col) of the covering array CA(3q^3 - 2; 3, k, q): Sherwood's expansion, row = r (q^3 - 1) + a - 1
    for CPHF row r and nonzero vector a in base q, and the last row all zero"""
    K = family(q)
    n = q ** 3 - 1
    if not 0 <= col < k <= K.width or not 0 <= row < 3 * n + 1:
        raise IndexError("entry outside the array")
    if row == 3 * n:
        return 0
    r, a = divmod(row, n)
    a += 1
    vec = (a // (q * q), (a // q) % q, a % q)
    h = K.column(col)[r]
    m = K.F.mul
    return m(vec[0], h[0]) ^ m(vec[1], h[1]) ^ m(vec[2], h[2])
