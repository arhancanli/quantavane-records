/* range_cells_B.c - Theorem K cells past every table, q = 4^s (program B of 2; family-theorems lane, 18 Nov 2026).
   Written separately from range_cells_A.py: fixed-width unsigned integers, prime powers by trial division (exact by
   definition), the class I bound found by scanning down from a floating cube root and the next prime power by
   scanning up from q1. Same classes and output format as program A (RANGE.md section 2), so the two outputs must be
   identical byte for byte.  build: cc -O2 -o range_cells_B range_cells_B.c   usage: range_cells_B SMAX [--nk-offset D]
   v2 (20 Nov 2026): 256-bit unsigned integers (four 64-bit words, every product checked for overflow) in place of
   unsigned __int128, so SMAX <= 30 (was 20: 4q^3 overflows 128 bits from s = 21). The logic and the output are those
   of v1 (kept as range_cells_B_v1_20261118.c); its output for SMAX = 20 must equal out_B.txt byte for byte.
   v2.1 (20 Nov 2026 18:3x, v2 kept as range_cells_B_v2_20261120.c): the cap is raised from SMAX <= 30 to SMAX <= 31,
   the real limit of this program: q1 and the next prime power must fit 64 bits for the trial division (small() exits
   otherwise), and q1 is about 1.14 * 4^s, below 2^64 only for s <= 31. Nothing else changes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;

typedef struct { u64 w[4]; } U;                        /* little-endian 256-bit unsigned integer */
static U uz(u64 x) { U r = {{x, 0, 0, 0}}; return r; }
static int ucmp(U a, U b) { for (int i = 3; i >= 0; i--) { if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1; } return 0; }
static U uadd(U a, U b) {
    U r; u128 c = 0;
    for (int i = 0; i < 4; i++) { c += (u128)a.w[i] + b.w[i]; r.w[i] = (u64)c; c >>= 64; }
    if (c) { fprintf(stderr, "overflow in add\n"); exit(3); }
    return r;
}
static U usub(U a, U b) {                              /* a - b, a >= b */
    if (ucmp(a, b) < 0) { fprintf(stderr, "negative in sub\n"); exit(3); }
    U r; u64 br = 0;
    for (int i = 0; i < 4; i++) { u128 t = (u128)a.w[i] - b.w[i] - br; r.w[i] = (u64)t; br = (t >> 64) ? 1 : 0; }
    return r;
}
static U umul(U a, U b) {                              /* exact or exit */
    u64 t[8] = {0};
    for (int i = 0; i < 4; i++) { u128 c = 0;
        for (int j = 0; j < 4; j++) { c += (u128)a.w[i] * b.w[j] + t[i + j]; t[i + j] = (u64)c; c >>= 64; }
        for (int k = i + 4; c && k < 8; k++) { c += t[k]; t[k] = (u64)c; c >>= 64; } }
    for (int k = 4; k < 8; k++) if (t[k]) { fprintf(stderr, "overflow in mul\n"); exit(3); }
    U r; memcpy(r.w, t, sizeof r.w); return r;
}
static U umuls(U a, u64 m) { return umul(a, uz(m)); }
static U udivs(U a, u64 d, u64 *rem) {                 /* a / d, d > 0 */
    U r; u128 cur = 0;
    for (int i = 3; i >= 0; i--) { cur = (cur << 64) | a.w[i]; r.w[i] = (u64)(cur / d); cur %= d; }
    if (rem) *rem = (u64)cur; return r;
}
static U ushl1(int k) { U r = uz(0); r.w[k / 64] = 1ULL << (k % 64); return r; }   /* 2^k, k < 256 */
static int uodd(U a) { return (int)(a.w[0] & 1); }
static long double uld(U a) { long double x = 0; for (int i = 3; i >= 0; i--) x = x * 18446744073709551616.0L + (long double)a.w[i]; return x; }
static U ufromld(long double x) {                      /* floor of a nonnegative value below 2^128 */
    U r = uz(0); long double hi = floorl(x / 18446744073709551616.0L);
    r.w[1] = (u64)hi; r.w[0] = (u64)(x - hi * 18446744073709551616.0L); return r;
}
static char *str(U x) {                                /* decimal, rotating static buffers */
    static char buf[32][96]; static int rr = 0; char *b = buf[rr++ & 31]; char t[96]; int n = 0; u64 d;
    if (ucmp(x, uz(0)) == 0) { strcpy(b, "0"); return b; }
    while (ucmp(x, uz(0)) > 0) { x = udivs(x, 10, &d); t[n++] = '0' + (int)d; }
    for (int i = 0; i < n; i++) b[i] = t[n - 1 - i];
    b[n] = 0; return b;
}

static int prime_power(u64 n) {          /* n = p^e, p prime, e >= 1 */
    if (n < 2) return 0;
    u64 p = 0;
    if (n % 2 == 0) p = 2;
    else for (u64 d = 3; d * d <= n; d += 2) if (n % d == 0) { p = d; break; }
    if (p == 0) return 1;                /* n prime */
    while (n % p == 0) n /= p;
    return n == 1;
}
static u64 small(U a) { if (a.w[1] || a.w[2] || a.w[3]) { fprintf(stderr, "value above 2^64\n"); exit(3); } return a.w[0]; }

static U rowsI(U qp, U q) {              /* restricted pair: 2Q^3 - Q^2 - 2(Q - q) */
    U q2 = umul(qp, qp), q3 = umul(q2, qp);
    return usub(usub(umuls(q3, 2), q2), umuls(usub(qp, q), 2));
}

int main(int argc, char **argv) {
    int smax = argc > 1 ? atoi(argv[1]) : 10; long long off = 0;
    for (int i = 2; i + 1 < argc; i++) if (!strcmp(argv[i], "--nk-offset")) off = atoll(argv[i + 1]);
    if (smax > 31) { fprintf(stderr, "SMAX <= 31 (q1 must fit 64 bits)\n"); return 2; }
    U total = uz(0), one = uz(1), two = uz(2);
    for (int s = 3; s <= smax; s++) {
        U q = ushl1(2 * s), qq = umul(q, q), q3 = umul(qq, q);
        U NK = usub(umuls(q3, 3), two);
        NK = off >= 0 ? uadd(NK, uz((u64)off)) : usub(NK, uz((u64)(-off)));
        U K = umuls(uadd(uadd(qq, q), one), 2);
        U c = uadd(ufromld(cbrtl(uld(udivs(NK, 2, 0)))), uz(4));   /* above the bound: checked by the loop */
        while (ucmp(rowsI(c, q), NK) <= 0) c = uadd(c, one);
        U q1 = usub(c, one);
        while (!(ucmp(rowsI(q1, q), NK) <= 0 && prime_power(small(q1)))) q1 = usub(q1, one);
        if (ucmp(q1, q) < 0) { fprintf(stderr, "no class I prime power\n"); return 1; }
        U q2 = uadd(q1, one);
        while (!prime_power(small(q2))) q2 = uadd(q2, one);
        /* Lemma R1: a 2-row CPHF over GF(x) has at most x^2+x+1 columns (x odd) or (x+1)^2 (x even) */
        U reachI = uodd(q1) ? uadd(uadd(umul(q1, q1), q1), one) : umul(uadd(q1, one), uadd(q1, one));
        U marginI = usub(rowsI(q2, q), NK);
        for (U x = q; ucmp(x, q1) <= 0; x = umuls(x, 2)) { U r = umul(uadd(x, one), uadd(x, one)); if (ucmp(r, reachI) > 0) reachI = r; }
        U v = uadd(q, two), reachIIa = uadd(umul(v, v), v);
        U reachIIb = uadd(umuls(ushl1(3 * s), 4), uz(3));
        U rowsIII = uadd(usub(umuls(q3, 4), umuls(qq, 5)), umuls(q, 2));
        U rowsIV = uadd(usub(umuls(q3, 4), umuls(qq, 6)), umuls(q, 3));        /* six restricted pairs */
        if (!(ucmp(rowsIII, NK) > 0 && ucmp(rowsIV, NK) > 0)) { fprintf(stderr, "class III/IV smaller\n"); return 1; }
        U R = reachI; if (ucmp(reachIIa, R) > 0) R = reachIIa; if (ucmp(reachIIb, R) > 0) R = reachIIb;
        if (ucmp(R, K) >= 0) { fprintf(stderr, "no cells\n"); return 1; }
        U cells = usub(K, R); total = uadd(total, cells);
        U T2 = udivs(usub(NK, usub(umuls(q3, 2), q)), small(usub(q, one)), 0);
        U T3 = usub(NK, umul(usub(q, one), usub(umuls(qq, 2), q)));
        printf("s %d q %s N_K %s K %s | I q1 %s reach %s next %s margin %s | IIa %s IIb %s | R %s k %s..%s cells %s total %s | V T2 %s T3 %s w %s\n",
               s, str(q), str(NK), str(K), str(q1), str(reachI), str(q2), str(marginI), str(reachIIa), str(reachIIb),
               str(R), str(uadd(R, one)), str(K), str(cells), str(total), str(T2), str(T3), str(udivs(uadd(R, two), 2, 0)));
    }
    return 0;
}
