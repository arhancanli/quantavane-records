/* range_cells_B.c - Theorem K cells past every table, q = 4^s (program B of 2; family-theorems lane, 18 Nov 2026).
   Written separately from range_cells_A.py: unsigned 128-bit integers, prime powers by trial division (exact by
   definition), the class I bound found by scanning down from a long-double cube root and the next prime power by
   scanning up from q1. Same classes and output format as program A (RANGE.md section 2), so the two outputs must be
   identical byte for byte.  build: cc -O2 -o range_cells_B range_cells_B.c   usage: range_cells_B SMAX [--nk-offset D] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;

static char *str(u128 x) {               /* decimal, rotating static buffers */
    static char buf[32][64]; static int r = 0; char *b = buf[r++ & 31]; char t[64]; int n = 0;
    if (x == 0) { strcpy(b, "0"); return b; }
    while (x) { t[n++] = '0' + (int)(x % 10); x /= 10; }
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

static u128 rowsI(u128 qp, u128 q) { return 2 * qp * qp * qp - qp * qp - 2 * (qp - q); }   /* restricted pair: 2Q^3 - Q^2 */

int main(int argc, char **argv) {
    int smax = argc > 1 ? atoi(argv[1]) : 10; long long off = 0;
    for (int i = 2; i + 1 < argc; i++) if (!strcmp(argv[i], "--nk-offset")) off = atoll(argv[i + 1]);
    if (smax > 20) { fprintf(stderr, "SMAX <= 20 (128-bit range)\n"); return 2; }
    u128 total = 0;
    for (int s = 3; s <= smax; s++) {
        u128 q = (u128)1 << (2 * s);
        u128 NK = 3 * q * q * q - 2 + (u128)(__int128)off;
        u128 K = 2 * (q * q + q + 1);
        u128 c = (u128)cbrtl((long double)(NK / 2)) + 4;   /* above the bound: checked by the loop */
        while (rowsI(c, q) <= NK) c++;
        u128 q1 = c - 1;
        while (!(rowsI(q1, q) <= NK && prime_power((u64)q1))) q1--;
        if (q1 < q) { fprintf(stderr, "no class I prime power\n"); return 1; }
        u128 q2 = q1 + 1;
        while (!prime_power((u64)q2)) q2++;
        /* Lemma R1: a 2-row CPHF over GF(x) has at most x^2+x+1 columns (x odd) or (x+1)^2 (x even) */
        u128 reachI = (q1 & 1) ? q1 * q1 + q1 + 1 : (q1 + 1) * (q1 + 1), marginI = rowsI(q2, q) - NK;
        for (u128 x = q; x <= q1; x <<= 1) { u128 r = (x + 1) * (x + 1); if (r > reachI) reachI = r; }
        u128 v = q + 2, reachIIa = v * v + v;
        u128 reachIIb = 4 * ((u128)1 << (3 * s)) + 3;
        u128 rowsIII = 4 * q * q * q - 5 * q * q + 2 * q, rowsIV = 4 * q * q * q - 6 * q * q + 3 * q;   /* six restricted pairs */
        if (!(rowsIII > NK && rowsIV > NK)) { fprintf(stderr, "class III/IV smaller\n"); return 1; }
        u128 R = reachI; if (reachIIa > R) R = reachIIa; if (reachIIb > R) R = reachIIb;
        if (R >= K) { fprintf(stderr, "no cells\n"); return 1; }
        u128 cells = K - R; total += cells;
        u128 T2 = (NK - (2 * q * q * q - q)) / (q - 1), T3 = NK - (q - 1) * (2 * q * q - q);
        printf("s %d q %s N_K %s K %s | I q1 %s reach %s next %s margin %s | IIa %s IIb %s | R %s k %s..%s cells %s total %s | V T2 %s T3 %s w %s\n",
               s, str(q), str(NK), str(K), str(q1), str(reachI), str(q2), str(marginI), str(reachIIa), str(reachIIb),
               str(R), str(R + 1), str(K), str(cells), str(total), str(T2), str(T3), str((R + 2) / 2));
    }
    return 0;
}
