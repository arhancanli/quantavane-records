# classcheck_C.py - a fourth, separately written re-derivation of the s-dependent quantities of RANGE.md sections 2 and 5
# (Quantavane, 20 Nov 2026): q1, the next prime power, the fusion margin per fused symbol, R(s), the cells and
# the totals, compared with out_A.txt; plus the class II-IV and P inequalities at each s. usage: python3 -I classcheck_C.py out_A.txt
# Independent re-derivation of RANGE.md's s-dependent quantities for s = 3..20 (19 Nov 2026), read against out_A.txt.
import sys
def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41):
        if n % p == 0: return n == p
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37,41):   # deterministic far beyond 1e13
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
def iroot(n, e):
    x = int(round(n ** (1.0 / e)))
    for y in (x - 2, x - 1, x, x + 1, x + 2):
        if y > 0 and y ** e == n: return y
    return None
def is_pp(n):
    if is_prime(n): return True
    for e in range(2, n.bit_length() + 1):
        r = iroot(n, e)
        if r and is_prime(r): return True
    return False
rows_I = lambda Q, q: 2*Q**3 - Q**2 - 2*(Q - q)          # RANGE.md class I floor with fusion
out = {}
for line in open(sys.argv[1]):
    w = line.split(); g = lambda k: int(w[w.index(k) + 1]); out[g("s")] = (g("q1"), g("next"), g("margin"), g("R"), g("cells"), g("total"), g("N_K"), g("K"))
tot = 0
for s in sorted(out):                                 # every s in the file (s = 3..20 on 20 Nov, 3..22 since)
    q = 4 ** s; NK = 3*q**3 - 2; K = 2*(q*q + q + 1)
    Q = int((NK / 2) ** (1/3)) + 3
    while rows_I(Q, q) > NK or not is_pp(Q): Q -= 1
    nx = Q + 1
    while not is_pp(nx): nx += 1
    assert rows_I(nx, q) > NK and rows_I(int((NK / 2) ** (1/3)) + 3, q) > NK
    R1 = Q*Q + Q + 1 if Q % 2 else (Q + 1)**2
    R = max(R1, (q + 1)**2)                        # class I incl. the power of 2 in [q, q1]
    IIa = (q + 2)**2 + q + 2; IIb = 4 * int(q ** 1.5) + 3 + 1; IId = q*q + 2*q + 2; IIe = q*q + 5*q
    margin = rows_I(nx, q) - NK
    per_sym = margin / (nx - q)
    cls3 = 4*q**3 - 5*q**2 + 2*q > NK; cls4 = 4*q**3 - 6*q**2 + 3*q > NK
    clsII_over = 3*nx**3 - 3*nx**2 - nx + 2*q > NK
    P_gap = (4*q**3 - 3*q**2) - NK
    cells = K - R; tot += cells
    o = out[s]
    ok = (o[0] == Q and o[1] == nx and o[2] == margin and o[3] == R and o[4] == cells and o[5] == tot and o[6] == NK and o[7] == K)
    print("s %2d q1 %d (odd %s) next %d | R %d = q1^2+q1+1: %s | IIa/R %.3f IIb/R %.2e IId/R %.3f IIe/R %.3f | margin/(next-q) = %.0f q | III %s IV %s II(Q>q) %s P gap %s | cells %d total %d | matches out_A: %s" % (
        s, Q, Q % 2 == 1, nx, R, R == Q*Q + Q + 1, IIa / R, IIb / R, IId / R, IIe / R, per_sym / q, cls3, cls4, clsII_over, P_gap > 0, cells, tot, ok))
