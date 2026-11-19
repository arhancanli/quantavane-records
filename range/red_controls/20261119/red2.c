/* kcheck_B.c - sampled exact checks of Theorem K's array at q = 4^s (program B of 2; family-theorems lane, 18 Nov 2026).
   Written separately from kcheck_A.py. Field model: GF(Q) = GF(2)[x]/(F), Q = 2^(6s) = q^3, F the first polynomial of
   degree 6s passing Rabin's test (one 128-bit word, 6s <= 96). GF(q) is the subfield fixed by z -> z^q. b = h^((q^2+q+1)/3)
   for the first non-cube h of GF(Q)^*, so b^3 lies in GF(q) and b does not. Independence of u, v, w over GF(q) is the
   Moore determinant det[[u,v,w],[u^q,v^q,w^q],[u^q^2,v^q^2,w^q^2]] != 0 (no coordinates are used). Scalars of GF(q) are
   norms z^(q^2+q+1) of random z. Same cases as kcheck_A.py (R, H1-H5); hard cases confirm that rows 1 and 3 fail.
   build: cc -O2 -o kcheck_B kcheck_B.c    usage: kcheck_B S SAMPLES SEED [inv|bq|mod|cube|nuc]; SAMPLES = 0: exhaustive H3, H4, H5 by symmetry
   exit 0 PASS, 1 FAIL, 2 refused.
   v2 (19 Nov 2026): range 1 <= s <= 20 (was 16). The field arithmetic is unchanged. mulF reduces after every shift, so
   every intermediate has degree <= n and fits one 128-bit word whenever n = 6s <= 127 (s <= 21); the old cap of 16 was a
   guard, not a limit of the arithmetic. As a separate check of that, every run first compares mulF with a second
   multiplication written apart (mul2w: the full 256-bit carry-less product held in two 128-bit words, built from 64-bit
   blocks, then reduced by long division) on 4,000 random pairs and 16 edge pairs from its own generator, and tests
   z^(2^n) = z on 8 random z; it refuses on any mismatch and reports to stderr only, so stdout (the logs) is byte for
   byte that of v1 at s <= 16. The exhaustive mode (SAMPLES = 0) is refused above s = 8 (its counters and tables would
   overflow; it was run for s <= 6 only). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;
static int n;            /* 6s */
static u128 F;           /* modulus with bit n set */
static int S2;           /* 2s: z^q is 2s squarings */

static u128 mulF(u128 a, u128 b, u128 f, int deg) {     /* a*b mod f, deg(a), deg(b) < deg */
    u128 r = 0;
    for (int i = deg - 1; i >= 0; i--) {
        r <<= 1; if ((r >> deg) & 1) r ^= f;
        if ((b >> i) & 1) r ^= a;
    }
    return r;
}
static u128 mul(u128 a, u128 b) { return mulF(a, b, F, n); }
static u128 powF(u128 a, u128 e, u128 f, int deg) { u128 r = 1; while (e) { if (e & 1) r = mulF(r, a, f, deg); a = mulF(a, a, f, deg); e >>= 1; } return r; }
static u128 pw(u128 a, u128 e) { return powF(a, e, F, n); }
static int degree(u128 a) { int d = -1; while (a) { d++; a >>= 1; } return d; }
static u128 pgcd(u128 a, u128 b) {
    while (b) { int db = degree(b); while (a && degree(a) >= db) a ^= b << (degree(a) - db); u128 t = a; a = b; b = t; }
    return a;
}
static u128 x2k(int k, u128 f, int deg) { u128 r = 2; for (int i = 0; i < k; i++) r = mulF(r, r, f, deg); return r; }
static int rabin(u128 f, int deg) {
    if (x2k(deg, f, deg) != 2) return 0;
    for (int p = 2; p <= deg; p++) {
        int prime = 1; for (int d = 2; d * d <= p; d++) if (p % d == 0) prime = 0;
        if (!prime || deg % p) continue;
        if (pgcd(f, x2k(deg / p, f, deg) ^ 2) != 1) return 0;
    }
    return 1;
}
static u128 frob(u128 z) { for (int i = 0; i < S2; i++) z = mul(z, z); return z; }   /* z^q */
static u128 inv(u128 z) { return pw(z, (((u128)1) << n) - 2); }
static int indep(u128 u, u128 v, u128 w) {
    u128 u1 = frob(u), v1 = frob(v), w1 = frob(w), u2 = frob(u1), v2 = frob(v1), w2 = frob(w1);
    u128 d = mul(u, mul(v1, w2) ^ mul(w1, v2)) ^ mul(v, mul(u1, w2) ^ mul(w1, u2)) ^ mul(w, mul(u1, v2) ^ mul(v1, u2));
    return d != 0;
}
/* ---- v2: a second multiplication, written apart from mulF, used only to cross-check it ---- */
static u128 clmul64(u64 a, u64 b) {                     /* carry-less 64 x 64 -> 128 */
    u128 r = 0, A = a;
    for (int i = 0; i < 64; i++) if ((b >> i) & 1) r ^= A << i;
    return r;
}
typedef struct { u128 hi, lo; } w256;                    /* a polynomial of degree < 256 in two 128-bit words */
static w256 clmul128(u128 a, u128 b) {                  /* full carry-less product, no reduction */
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 p00 = clmul64(a0, b0), p11 = clmul64(a1, b1), mid = clmul64(a0, b1) ^ clmul64(a1, b0);
    w256 r; r.lo = p00 ^ (mid << 64); r.hi = p11; return r;
}
static int bit256(w256 x, int i) { return i >= 128 ? (int)((x.hi >> (i - 128)) & 1) : (int)((x.lo >> i) & 1); }
static w256 shl256(u128 f, int k) {                     /* f * x^k, 0 <= k < 256 - deg f */
    w256 r;
    if (k == 0) { r.hi = 0; r.lo = f; }
    else if (k < 128) { r.lo = f << k; r.hi = f >> (128 - k); }
    else { r.lo = 0; r.hi = f << (k - 128); }
    return r;
}
static u128 mul2w(u128 a, u128 b, u128 f, int deg) {    /* a*b mod f by full product then long division */
    w256 p = clmul128(a, b);
    for (int i = 255; i >= deg; i--) if (bit256(p, i)) { w256 t = shl256(f, i - deg); p.hi ^= t.hi; p.lo ^= t.lo; }
    if (p.hi) return ~(u128)0;                          /* cannot happen after a correct division */
    return p.lo;
}
static int mulcheck(int s) {                            /* 0 = agree; own generator, so `st` below is untouched */
    u64 x = 0x243F6A8885A308D3ULL ^ (u64)s;
    u128 mask = (((u128)1) << n) - 1;
    int bad = 0, cnt = 0;
    u128 edge[4] = {1, ((u128)1) << (n - 1), mask, mask ^ 1};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { cnt++; if (mulF(edge[i], edge[j], F, n) != mul2w(edge[i], edge[j], F, n)) bad++; }
    for (int i = 0; i < 4000; i++) {
        u128 a, b;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; a = (u128)x << 64;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; a = (a | x) & mask;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; b = (u128)x << 64;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; b = (b | x) & mask;
        cnt++; if (mulF(a, b, F, n) != mul2w(a, b, F, n)) bad++;
    }
    int fz = 0;
    for (int i = 0; i < 8; i++) {                       /* z^(2^n) = z in GF(2^n), by mul2w squarings */
        u128 z; x ^= x << 13; x ^= x >> 7; x ^= x << 17; z = ((u128)x << 64);
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; z = (z | x) & mask;
        u128 y = z; for (int k = 0; k < n; k++) y = mul2w(y, y, F, n);
        if (y != z) fz++;
    }
    fprintf(stderr, "kcheck_B mulcheck s %d n %d: mulF vs mul2w %d/%d disagree; z^(2^n) = z fails %d/8\n", s, n, bad, cnt, fz);
    return bad || fz;
}

static u64 st = 88172645463325252ULL;
static u64 rnd64(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static u128 rnd_el(void) { u128 z; do { z = (((u128)rnd64() << 64) | rnd64()) & ((((u128)1) << n) - 1); } while (!z); return z; }
static u128 NORM;         /* q^2 + q + 1 */
static u128 rnd_scalar(int nonzero) { if (!nonzero && rnd64() % 16 == 0) return 0; return pw(rnd_el(), NORM); }
static u128 b, binv; static const char *T = "";
static int same_point(u128 x, u128 y) { u128 r = mul(x, inv(y)); return frob(r) == r; }   /* x/y in GF(q) */
typedef struct { u128 r[3]; } col;
static col column(int E, u128 x) {
    col c;
    if (!E) { c.r[0] = x; c.r[1] = inv(x); c.r[2] = x; return c; }
    u128 d2 = strcmp(T, "nuc") ? b : (1 ^ b);   /* nuc: a row-2 twist off the conic that is not its nucleus */
    c.r[0] = x; c.r[1] = strcmp(T, "inv") ? mul(d2, inv(x)) : mul(b, x); c.r[2] = mul(b, x); return c;
}
static void two_w(u128 *w1, u128 *w2) {
    for (;;) {
        u128 a1 = rnd_scalar(0), c1 = rnd_scalar(0), a2 = rnd_scalar(0), c2 = rnd_scalar(0);
        *w1 = a1 ^ mul(c1, b); *w2 = a2 ^ mul(c2, b);
        if (*w1 && *w2 && mul(a1, c2) != mul(c1, a2)) return;
    }
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 2; }
    int s = atoi(argv[1]), samples = atoi(argv[2]); st ^= (u64)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    if (argc > 4) T = argv[4];
    if (s < 1 || s > 20) { fprintf(stderr, "1 <= s <= 20\n"); return 2; }      /* v2: was 16; n = 6s <= 120 < 128 */
    if (samples == 0 && s > 8) { printf("REFUSED: exhaustive mode needs s <= 8\n"); return 2; }
    n = 6 * s; S2 = 2 * s;
    for (u128 f = (((u128)1) << n) | 1; ; f += 2) if (rabin(f, n)) { F = f; break; }
    if (!strcmp(T, "mod")) { F ^= 2; if (rabin(F, n)) F ^= 6; }
    if (!rabin(F, n)) { printf("REFUSED: modulus not irreducible\n"); return 2; }
    if (mulcheck(s)) { printf("REFUSED: multiplication cross-check\n"); return 2; }
    u128 q = ((u128)1) << S2; NORM = q * q + q + 1;
    u128 cube = ((((u128)1) << n) - 1) / 3, h = 2;
    if (!strcmp(T, "cube")) { while (pw(h, cube) != 1) h++; } else { while (pw(h, cube) == 1) h++; }
    b = pw(h, NORM / 3);
    if (!strcmp(T, "bq")) b = pw(rnd_el(), NORM);                 /* an element of GF(q) */
    if (frob(b) == b && strcmp(T, "bq")) { printf("REFUSED: b lies in GF(q)\n"); return 2; }
    if (frob(mul(b, mul(b, b))) != mul(b, mul(b, b))) { printf("REFUSED: b^3 not in GF(q)\n"); return 2; }
    binv = inv(b);
    for (int i = 0; i < 20; i++) { u128 z = rnd_el(); if (mul(z, inv(z)) != 1) { printf("REFUSED: field self-test\n"); return 2; } }
    if (samples == 0) {                     /* exhaustive Lemma 1/2 cases via the symmetry x -> yx (y = 1) */
        u128 Q1 = (((u128)1) << n) - 1, z = 0, qm1 = q - 1;
        for (u128 t = 2; ; t++) {           /* z generates GF(q)^*: order exactly q - 1 */
            z = pw(t, Q1 / qm1); int ok = (z != 0);
            u128 m = qm1;
            for (u128 p = 2; p * p <= m; p++) if (m % p == 0) { if (pw(z, qm1 / p) == 1) ok = 0; while (m % p == 0) m /= p; }
            if (m > 1 && pw(z, qm1 / m) == 1) ok = 0;
            if (ok) break;
        }
        int Wn = (int)q + 1; u128 *W = malloc(sizeof(u128) * Wn), *Fq = malloc(sizeof(u128) * (size_t)q);
        Fq[0] = 0; u128 e = 1; for (u128 j = 1; j < q; j++) { Fq[j] = e; e = mul(e, z); }
        if (e != 1) { printf("REFUSED: GF(q) generator\n"); return 2; }
        for (int i = 0; i < (int)q; i++) W[i] = Fq[i] ^ b;              /* u + b */
        W[q] = 1;                                                        /* the point 1 */
        long long unc3 = 0, unc4 = 0, pre = 0, n3 = 0, unc5 = 0, n5 = 0;
        col B1 = column(0, 1), E1 = column(1, 1);
        for (int i = 0; i < Wn; i++) for (int j = i + 1; j < Wn; j++) {
            col c3[3] = {column(0, W[i]), column(0, W[j]), E1};
            col c4[3] = {B1, column(1, mul(binv, W[i])), column(1, mul(binv, W[j]))};
            int o3[3], o4[3];
            for (int r = 0; r < 3; r++) { o3[r] = indep(c3[0].r[r], c3[1].r[r], c3[2].r[r]); o4[r] = indep(c4[0].r[r], c4[1].r[r], c4[2].r[r]); }
            if (o3[0] || o3[2] || o4[0] || o4[2]) pre++;
            if (!o3[0] && !o3[1] && !o3[2]) unc3++;
            if (!o4[0] && !o4[1] && !o4[2]) unc4++;
            n3++;
        }
        /* H5: B_1, E_1 and every third column B_x / E_x (x != 1), x over all points: x = (1,a,b')-type via z-powers */
        long long Np = (long long)(q * q + q + 1);
        u128 g = 0;                         /* a generator of GF(Q)^* modulo GF(q)^*: points are g^i, i < q^2+q+1 */
        for (u128 t = 2; ; t++) { g = t; u128 m = NORM, ok = 1, mm = m; u128 gn = pw(g, qm1);
            for (u128 p = 2; p * p <= mm; p++) if (mm % p == 0) { if (pw(gn, m / p) == 1) ok = 0; while (mm % p == 0) mm /= p; }
            if (mm > 1 && pw(gn, m / mm) == 1) ok = 0; if (ok) break; }
        u128 x = g;                        /* the class of g has order q^2+q+1, so g^i, 0 < i < N, are the other points */
        for (long long i = 1; i < Np; i++, x = mul(x, g)) for (int k = 0; k < 2; k++) {
            col c5[3] = {B1, E1, column(k, x)}; int any = 0;
            for (int r = 0; r < 3; r++) any |= indep(c5[0].r[r], c5[1].r[r], c5[2].r[r]);
            if (!any) unc5++; n5++;
        }
        int okall = (unc3 == 0 && unc4 == 0 && unc5 == 0 && (pre == 0 || *T));
        printf("kcheck_B EXHAUSTIVE s %d n %d tamper %s | H3 %lld/%lld uncovered, H4 %lld/%lld, precondition broken %lld | H5 %lld/%lld | %s\n",
               s, n, *T ? T : "-", unc3, n3, unc4, n3, pre, unc5, n5, okall ? "PASS" : "FAIL");
        return okall ? 0 : 1;
    }
    const char *names[6] = {"R", "H1", "H2", "H3", "H4", "H5"}; int unc[6] = {0}, pre[6] = {0}, fails = 0, badpre = 0;
    for (int cs = 0; cs < 6; cs++) for (int it = 0; it < samples; it++) {
        col c[3];
        if (cs == 0) {
            int k[3]; u128 x[3];
            for (int j = 0; j < 3; j++) { k[j] = rnd64() & 1; x[j] = rnd_el(); }
            if ((k[0] == k[1] && same_point(x[0], x[1])) || (k[0] == k[2] && same_point(x[0], x[2])) || (k[1] == k[2] && same_point(x[1], x[2]))) continue;
            for (int j = 0; j < 3; j++) c[j] = column(k[j], x[j]);
        } else if (cs == 1 || cs == 2) {
            u128 x1 = rnd_el(), x2 = rnd_el(); if (same_point(x1, x2)) continue;
            u128 x3 = mul(rnd_scalar(1), x1) ^ mul(rnd_scalar(1), x2);
            c[0] = column(cs == 2, x1); c[1] = column(cs == 2, x2); c[2] = column(cs == 2, x3);
        } else if (cs == 3) {
            u128 y = rnd_el(), w1, w2; two_w(&w1, &w2);
            c[0] = column(0, mul(y, w1)); c[1] = column(0, mul(y, w2)); c[2] = column(1, y);
        } else if (cs == 4) {
            u128 x = rnd_el(), w1, w2; two_w(&w1, &w2);
            c[0] = column(0, x); c[1] = column(1, mul(mul(x, binv), w1)); c[2] = column(1, mul(mul(x, binv), w2));
        } else {
            u128 x = rnd_el(), x3 = rnd_el(); int k3 = rnd64() & 1; if (same_point(x, x3)) continue;
            c[0] = column(0, x); c[1] = column(1, x); c[2] = column(k3, x3);
        }
        int ok[3]; for (int r = 0; r < 3; r++) ok[r] = indep(c[0].r[r], c[1].r[r], c[2].r[r]);
        if (cs >= 1 && cs <= 4 && (ok[0] || ok[2])) pre[cs]++;
        if (!ok[0] && !ok[1] && !ok[2]) unc[cs]++;
    }
    for (int cs = 0; cs < 6; cs++) { fails += unc[cs]; if (!*T) badpre += pre[cs]; }
    unsigned long long Fhi = (unsigned long long)(F >> 64), Flo = (unsigned long long)F;
    printf("kcheck_B s %d n %d F %llx:%016llx samples %d seed %s tamper %s | uncovered/precondition-broken", s, n, Fhi, Flo, samples, argv[3], *T ? T : "-");
    for (int cs = 0; cs < 6; cs++) printf(" %s:%d/%d", names[cs], unc[cs], pre[cs]);
    printf(" | %s\n", (fails == 0 && badpre == 0) ? "PASS" : "FAIL");
    return (fails == 0 && badpre == 0) ? 0 : 1;
}
