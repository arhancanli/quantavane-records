#!/usr/bin/env python3
"""range_cells_A.py - Theorem K cells past every table for q = 4^s, 3 <= s <= SMAX (program A of 2).

Family-theorems lane, 18 Nov 2026. Written separately from range_cells_B.c (that one uses int128 and trial
division; this one uses Python integers, Miller-Rabin with the first 13 prime bases, which is exact below 3.3e24,
and integer k-th roots for prime powers). Competitor classes and their bounds are proved in RANGE.md section 2.

Ours: N_K = 3q^3 - 2 rows for every k <= K = 2(q^2+q+1) (Theorem K).
 I    any covering array from a CPHF with at most 2 rows over a prime power q' >= q, symbols fused to q at
      2 rows per symbol: at least 2q'^3 - q'^2 - 2(q'-q) rows (Sherwood with a Colbourn-Lanus restricted pair); reach at most q'^2+q'+1 for odd q' and (q'+1)^2
      for even q' (Lemma R1, any 2-row CPHF, published or not).
      Competes when its rows <= N_K (a tie removes the cell). q1 = the largest such q'.
 IIa  PHF(3; k, v, 3) composed with an arc (v <= q+2): reach <= v^2 + v (Shangguan-Ge 2016, Theorem 1.5).
 IIb  random oversampling with 3 rows (any restriction): reach <= 4 q^(3/2) + 3 (Lemma R3).
 III  plications, PHF o RMS: at least 4q^3 - 5q^2 + 2q rows (checked > N_K below).
 IV   CPHF arrays with >= 4 rows, all pairs restricted: at least 4q^3 - 6q^2 + 3q rows (checked > N_K below).
 V    Roux-type recursions: reported as thresholds T2, T3 (assumption P, RANGE.md section 3).
Cells: k in (R, K], R = max reach over I, IIa, IIb.
usage: python3 -I range_cells_A.py SMAX [--nk-offset D]   (D != 0 is a tamper control)"""
import sys

BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


def is_prime(n):
    if n < 2:
        return False
    for p in BASES:
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def iroot(n, e):
    """largest x with x^e <= n"""
    lo, hi = 0, 1
    while hi ** e <= n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** e <= n:
            lo = mid
        else:
            hi = mid
    return lo


def is_prime_power(n):
    if n < 2:
        return False
    for e in range(n.bit_length(), 0, -1):
        r = iroot(n, e)
        if r ** e == n and is_prime(r):
            return True
    return False


def rows_I(qp, q):
    return 2 * qp ** 3 - qp * qp - 2 * (qp - q)


def main():
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    off = int(sys.argv[sys.argv.index("--nk-offset") + 1]) if "--nk-offset" in sys.argv else 0
    total = 0
    for s in range(3, smax + 1):
        q = 4 ** s
        NK = 3 * q ** 3 - 2 + off
        K = 2 * (q * q + q + 1)
        # class I: largest prime power q' with rows_I(q') <= NK; rows_I is increasing in q'
        top = iroot(NK // 2, 3) + 2
        while rows_I(top, q) <= NK:
            top += 1
        qp = top
        while not (qp >= q and rows_I(qp, q) <= NK and is_prime_power(qp)):
            qp -= 1
        q1 = qp
        q2 = q1 + 1                             # v2 (20 Nov 2026): scan up from q1 + 1; v1 began at `top`, which can lie
        while not is_prime_power(q2):           # one above the boundary and so skip a prime power at boundary + 1
            q2 += 1
        f = lambda x: x * x + x + 1 if x % 2 else (x + 1) ** 2
        e = q
        while 2 * e <= q1:
            e *= 2                              # the largest power of 2 in [q, q1], the only even prime powers there
        reachI = max(f(q1), f(e))
        marginI = rows_I(q2, q) - NK          # rows by which the next prime power's fused 2-row array exceeds ours
        v = q + 2
        reachIIa = v * v + v
        reachIIb = 4 * 8 ** s + 3              # q^(3/2) = 8^s exactly
        rowsIII = 4 * q ** 3 - 5 * q * q + 2 * q
        rowsIV = 4 * q ** 3 - 6 * q * q + 3 * q             # four rows, all six pairs restricted
        assert rowsIII > NK and rowsIV > NK, "a class III/IV construction would have fewer rows"
        R = max(reachI, reachIIa, reachIIb)
        assert R < K, "no cells"
        cells = K - R
        total += cells
        # class V thresholds (assumption P): doubling from a strength-3 ingredient of 2q^3 - q rows needs a
        # strength-2 array with <= T2 rows at width ceil((R+1)/2); with strength-2 ingredients of 2q^2 - q rows
        # it needs a strength-3 array with <= T3 rows at that width.
        T2 = (NK - (2 * q ** 3 - q)) // (q - 1)
        T3 = NK - (q - 1) * (2 * q * q - q)
        print("s %d q %d N_K %d K %d | I q1 %d reach %d next %d margin %d | IIa %d IIb %d | R %d k %d..%d cells %d total %d | V T2 %d T3 %d w %d"
              % (s, q, NK, K, q1, reachI, q2, marginI, reachIIa, reachIIb, R, R + 1, K, cells, total, T2, T3, (R + 2) // 2))


if __name__ == "__main__":
    main()
