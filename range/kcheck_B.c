/* kcheck_B.c - sampled exact checks of Theorem K's array at q = 4^s (program B of 2; family-theorems lane, 18 Nov 2026).
   Written separately from kcheck_A.py. Field model: GF(Q) = GF(2)[x]/(F), Q = 2^(6s) = q^3, F the first polynomial of
   degree 6s passing Rabin's test. GF(q) is the subfield fixed by z -> z^q. b = h^((q^2+q+1)/3)
   for the first non-cube h of GF(Q)^*, so b^3 lies in GF(q) and b does not. Independence of u, v, w over GF(q) is the
   Moore determinant det[[u,v,w],[u^q,v^q,w^q],[u^q^2,v^q^2,w^q^2]] != 0 (no coordinates are used). Scalars of GF(q) are
   norms z^(q^2+q+1) of random z. Same cases as kcheck_A.py (R, H1-H5); hard cases confirm that rows 1 and 3 fail.
   build: cc -O2 -o kcheck_B kcheck_B.c    usage: kcheck_B S SAMPLES SEED [inv|bq|mod|cube|nuc]; SAMPLES = 0: exhaustive H3, H4, H5 by symmetry
   exit 0 PASS, 1 FAIL, 2 refused.
   v2 (19 Nov 2026, kept as kcheck_B_v2_20261119.c): s <= 20 in one 128-bit word, with a cross-check of the multiply.
   v3 (20 Nov 2026): 1 <= s <= 30. Field elements, polynomials and exponents are multi-word (four 64-bit words; n = 6s
   <= 180 bits). The MAIN multiplication is now the full carry-less product (64 x 64-bit blocks into eight words) followed
   by long division by the modulus. Before anything else, every run compares it with a second multiplication written
   apart (shift-and-reduce after every bit, the method of v1/v2, on four words) on 4,024 pairs from its own generator,
   and tests z^(2^n) = z on 8 random z; it refuses on any mismatch and reports to stderr only. Random elements take
   ceil(n/64) words from the stream, highest word first (two words, as in v1/v2, for n <= 128), and every pair of
   random draws is now made in a fixed order, so stdout at s <= 20 must be that of v1/v2 byte for byte. The exhaustive
   mode (SAMPLES = 0) is refused above s = 8 (it was run for s <= 6 only).
   v3.1 (20 Nov 2026 18:3x, v3 kept as kcheck_B_v3_20261120.c): the cap is raised from s <= 30 to s <= 40; nothing else
   changes. The four-word layout is valid to n = 255 (s <= 42): elements, F and exponents fit four words, the shift-and-
   reduce accumulator needs n + 1 <= 256 bits, and the product buffer of eight words holds degree 2n - 2 <= 508. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;
#define NWMAX 4                                   /* words per element / polynomial / exponent */
typedef struct { u64 w[NWMAX]; } E;               /* little-endian: bit i of word i/64 = coefficient of x^i */
static int n;            /* 6s */
static E F;              /* modulus with bit n set */
static int S2;           /* 2s: z^q is 2s squarings */

static E ez(u64 x) { E r = {{x, 0, 0, 0}}; return r; }
static E e128(u128 x) { E r = {{(u64)x, (u64)(x >> 64), 0, 0}}; return r; }
static E exor(E a, E b) { for (int i = 0; i < NWMAX; i++) a.w[i] ^= b.w[i]; return a; }
static int eeq(E a, E b) { for (int i = 0; i < NWMAX; i++) if (a.w[i] != b.w[i]) return 0; return 1; }
static int eis0(E a) { for (int i = 0; i < NWMAX; i++) if (a.w[i]) return 0; return 1; }
static int eis1(E a) { return eeq(a, ez(1)); }
static int ebit(E a, int i) { return (int)((a.w[i >> 6] >> (i & 63)) & 1); }
static E esetbit(E a, int i) { a.w[i >> 6] |= 1ULL << (i & 63); return a; }
static int edeg(E a) { for (int i = NWMAX - 1; i >= 0; i--) if (a.w[i]) return 64 * i + 63 - __builtin_clzll(a.w[i]); return -1; }
static E eshl(E a, int k) {                       /* a * x^k, truncated to NWMAX words (callers keep degrees small) */
    E r = ez(0); int wq = k >> 6, bq = k & 63;
    for (int i = NWMAX - 1; i >= 0; i--) { if (i - wq < 0) break;
        u64 v = a.w[i - wq] << bq; if (bq && i - wq - 1 >= 0) v |= a.w[i - wq - 1] >> (64 - bq); r.w[i] = v; }
    return r;
}
static E emask(int k) { E r = ez(0); for (int i = 0; i < k; i++) r = esetbit(r, i); return r; }   /* 2^k - 1 */

/* ---- main multiplication: full carry-less product, then long division by f (deg f = deg) ---- */
static u128 clmul64(u64 a, u64 b) {               /* carry-less 64 x 64 -> 128 */
    u128 r = 0, A = a;
    for (int i = 0; i < 64; i++) if ((b >> i) & 1) r ^= A << i;
    return r;
}
static E mulM(E a, E b, E f, int deg) {            /* a*b mod f, deg(a), deg(b) < deg */
    u64 t[2 * NWMAX] = {0};
    int nw = (deg + 63) / 64;
    for (int i = 0; i < nw; i++) if (a.w[i]) for (int j = 0; j < nw; j++) if (b.w[j]) {
        u128 p = clmul64(a.w[i], b.w[j]); t[i + j] ^= (u64)p; t[i + j + 1] ^= (u64)(p >> 64); }
    for (int i = 2 * deg - 2; i >= deg; i--) if ((t[i >> 6] >> (i & 63)) & 1) {   /* subtract f * x^(i - deg) */
        int sh = i - deg, wq = sh >> 6, bq = sh & 63;
        for (int k = 0; k < NWMAX; k++) if (f.w[k]) {
            t[k + wq] ^= f.w[k] << bq;
            if (bq && k + wq + 1 < 2 * NWMAX) t[k + wq + 1] ^= f.w[k] >> (64 - bq);
        }
    }
    E r; for (int i = 0; i < NWMAX; i++) r.w[i] = t[i];
    for (int i = NWMAX; i < 2 * NWMAX; i++) if (t[i]) { E bad = ez(0); bad.w[NWMAX - 1] = ~0ULL; return bad; }  /* cannot happen */
    return r;
}
/* ---- second multiplication, for the cross-check only: shift and reduce after every bit (v1/v2's method) ---- */
static E mulS(E a, E b, E f, int deg) {
    E r = ez(0);
    for (int i = deg - 1; i >= 0; i--) {
        r = eshl(r, 1); if (ebit(r, deg)) r = exor(r, f);
        if (ebit(b, i)) r = exor(r, a);
    }
    return r;
}
static E mul(E a, E b) { return mulM(a, b, F, n); }
static E powM(E a, E e, E f, int deg) { E r = ez(1); int top = edeg(e); for (int i = 0; i <= top; i++) { if (ebit(e, i)) r = mulM(r, a, f, deg); a = mulM(a, a, f, deg); } return r; }
static E pw(E a, E e) { return powM(a, e, F, n); }
static E pwi(E a, u128 e) { return pw(a, e128(e)); }
static E pgcd(E a, E b) {
    while (!eis0(b)) { int db = edeg(b); while (!eis0(a) && edeg(a) >= db) a = exor(a, eshl(b, edeg(a) - db)); E t = a; a = b; b = t; }
    return a;
}
static E x2k(int k, E f, int deg) { E r = ez(2); for (int i = 0; i < k; i++) r = mulM(r, r, f, deg); return r; }
static int rabin(E f, int deg) {
    if (!eeq(x2k(deg, f, deg), ez(2))) return 0;
    for (int p = 2; p <= deg; p++) {
        int prime = 1; for (int d = 2; d * d <= p; d++) if (p % d == 0) prime = 0;
        if (!prime || deg % p) continue;
        if (!eis1(pgcd(f, exor(x2k(deg / p, f, deg), ez(2))))) return 0;
    }
    return 1;
}
static E frob(E z) { for (int i = 0; i < S2; i++) z = mul(z, z); return z; }   /* z^q */
static E inv(E z) { E e = emask(n); e.w[0] ^= 1; return pw(z, e); }       /* z^(2^n - 2) */
static int indep(E u, E v, E w) {
    E u1 = frob(u), v1 = frob(v), w1 = frob(w), u2 = frob(u1), v2 = frob(v1), w2 = frob(w1);
    E d = exor(exor(mul(u, exor(mul(v1, w2), mul(w1, v2))), mul(v, exor(mul(u1, w2), mul(w1, u2)))), mul(w, exor(mul(u1, v2), mul(v1, u2))));
    return !eis0(d);
}
static E ediv_small(E a, u64 d) { E r; u128 cur = 0; for (int i = NWMAX - 1; i >= 0; i--) { cur = (cur << 64) | a.w[i]; r.w[i] = (u64)(cur / d); cur %= d; } return r; }

static int mulcheck(int s) {                      /* 0 = agree; own generator, so `st` below is untouched */
    u64 x = 0x243F6A8885A308D3ULL ^ (u64)s;
    int nw = (n + 63) / 64, bad = 0, cnt = 0, fz = 0;
    E mask = emask(n), edge[4];
    edge[0] = ez(1); edge[1] = esetbit(ez(0), n - 1); edge[2] = mask; edge[3] = mask; edge[3].w[0] ^= 1;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { cnt++; if (!eeq(mulM(edge[i], edge[j], F, n), mulS(edge[i], edge[j], F, n))) bad++; }
    for (int i = 0; i < 4008; i++) {
        E a = ez(0), b = ez(0);
        for (int k = nw - 1; k >= 0; k--) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; a.w[k] = x; }
        for (int k = nw - 1; k >= 0; k--) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b.w[k] = x; }
        for (int k = 0; k < NWMAX; k++) { a.w[k] &= mask.w[k]; b.w[k] &= mask.w[k]; }
        if (i >= 4000) { if (i & 1) a = b; else b = a; }                 /* 8 squarings among the pairs */
        cnt++; if (!eeq(mulM(a, b, F, n), mulS(a, b, F, n))) bad++;
    }
    for (int i = 0; i < 8; i++) {                 /* z^(2^n) = z in GF(2^n), by squarings with the main multiplication */
        E z = ez(0);
        for (int k = nw - 1; k >= 0; k--) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; z.w[k] = x; }
        for (int k = 0; k < NWMAX; k++) z.w[k] &= mask.w[k];
        E y = z; for (int k = 0; k < n; k++) y = mul(y, y);
        if (!eeq(y, z)) fz++;
    }
    fprintf(stderr, "kcheck_B mulcheck s %d n %d: main (full product + division) vs shift-reduce %d/%d disagree; z^(2^n) = z fails %d/8\n", s, n, bad, cnt, fz);
    return bad || fz;
}

static u64 st = 88172645463325252ULL;
static u64 rnd64(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static E rnd_el(void) {                           /* highest word drawn first, as (rnd64() << 64) | rnd64() in v1/v2 */
    int nw = (n + 63) / 64; E mask = emask(n), z;
    if (nw < 2) nw = 2;                           /* v1/v2 drew two words also for n <= 64 */
    do { z = ez(0); for (int k = nw - 1; k >= 0; k--) z.w[k] = rnd64();
         for (int k = 0; k < NWMAX; k++) z.w[k] &= mask.w[k]; } while (eis0(z));
    return z;
}
static E NORM;            /* q^2 + q + 1 */
static E rnd_scalar(int nonzero) { if (!nonzero && rnd64() % 16 == 0) return ez(0); return pw(rnd_el(), NORM); }
static E b, binv; static const char *T = "";
static int same_point(E x, E y) { E r = mul(x, inv(y)); return eeq(frob(r), r); }   /* x/y in GF(q) */
typedef struct { E r[3]; } col;
static col column(int Ek, E x) {
    col c;
    if (!Ek) { c.r[0] = x; c.r[1] = inv(x); c.r[2] = x; return c; }
    E d2 = strcmp(T, "nuc") ? b : exor(ez(1), b);   /* nuc: a row-2 twist off the conic that is not its nucleus */
    c.r[0] = x; c.r[1] = strcmp(T, "inv") ? mul(d2, inv(x)) : mul(b, x); c.r[2] = mul(b, x); return c;
}
static void two_w(E *w1, E *w2) {
    for (;;) {
        E a1 = rnd_scalar(0), c1 = rnd_scalar(0), a2 = rnd_scalar(0), c2 = rnd_scalar(0);
        *w1 = exor(a1, mul(c1, b)); *w2 = exor(a2, mul(c2, b));
        if (!eis0(*w1) && !eis0(*w2) && !eeq(mul(a1, c2), mul(c1, a2))) return;
    }
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 2; }
    int s = atoi(argv[1]), samples = atoi(argv[2]); st ^= (u64)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    if (argc > 4) T = argv[4];
    if (s < 1 || s > 40) { fprintf(stderr, "1 <= s <= 40\n"); return 2; }      /* v3.1: n = 6s <= 240 < 64 NWMAX */
    if (samples == 0 && s > 8) { printf("REFUSED: exhaustive mode needs s <= 8\n"); return 2; }
    n = 6 * s; S2 = 2 * s;
    for (E f = esetbit(ez(1), n); ; f.w[0] += 2) if (rabin(f, n)) { F = f; break; }
    if (!strcmp(T, "mod")) { F.w[0] ^= 2; if (rabin(F, n)) F.w[0] ^= 6; }
    if (!rabin(F, n)) { printf("REFUSED: modulus not irreducible\n"); return 2; }
    if (mulcheck(s)) { printf("REFUSED: multiplication cross-check\n"); return 2; }
    E q = esetbit(ez(0), S2);
    NORM = esetbit(esetbit(ez(1), S2), 2 * S2);    /* q^2 + q + 1 */
    E cube = ediv_small(emask(n), 3), h = ez(2);
    if (!strcmp(T, "cube")) { while (!eis1(pw(h, cube))) h.w[0]++; } else { while (eis1(pw(h, cube))) h.w[0]++; }
    b = pw(h, ediv_small(NORM, 3));
    if (!strcmp(T, "bq")) b = pw(rnd_el(), NORM);                 /* an element of GF(q) */
    if (eeq(frob(b), b) && strcmp(T, "bq")) { printf("REFUSED: b lies in GF(q)\n"); return 2; }
    if (!eeq(frob(mul(b, mul(b, b))), mul(b, mul(b, b)))) { printf("REFUSED: b^3 not in GF(q)\n"); return 2; }
    binv = inv(b);
    for (int i = 0; i < 20; i++) { E z = rnd_el(); if (!eis1(mul(z, inv(z)))) { printf("REFUSED: field self-test\n"); return 2; } }
    if (samples == 0) {                     /* exhaustive Lemma 1/2 cases via the symmetry x -> yx (y = 1); s <= 8 */
        u128 qi = ((u128)1) << S2, Q1 = (((u128)1) << n) - 1, qm1 = qi - 1, normi = qi * qi + qi + 1; E z = ez(0);
        for (u128 t = 2; ; t++) {           /* z generates GF(q)^*: order exactly q - 1 */
            z = pwi(e128(t), Q1 / qm1); int ok = !eis0(z);
            u128 m = qm1;
            for (u128 p = 2; p * p <= m; p++) if (m % p == 0) { if (eis1(pwi(z, qm1 / p))) ok = 0; while (m % p == 0) m /= p; }
            if (m > 1 && eis1(pwi(z, qm1 / m))) ok = 0;
            if (ok) break;
        }
        int Wn = (int)qi + 1; E *W = malloc(sizeof(E) * Wn), *Fq = malloc(sizeof(E) * (size_t)qi);
        Fq[0] = ez(0); E e = ez(1); for (u128 j = 1; j < qi; j++) { Fq[j] = e; e = mul(e, z); }
        if (!eis1(e)) { printf("REFUSED: GF(q) generator\n"); return 2; }
        for (int i = 0; i < (int)qi; i++) W[i] = exor(Fq[i], b);       /* u + b */
        W[qi] = ez(1);                                                   /* the point 1 */
        long long unc3 = 0, unc4 = 0, pre = 0, n3 = 0, unc5 = 0, n5 = 0;
        col B1 = column(0, ez(1)), E1 = column(1, ez(1));
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
        /* H5: B_1, E_1 and every third column B_x / E_x (x != 1), x over all points */
        long long Np = (long long)normi;
        E g = ez(0);                        /* a generator of GF(Q)^* modulo GF(q)^*: points are g^i, i < q^2+q+1 */
        for (u128 t = 2; ; t++) { g = e128(t); u128 m = normi, ok = 1, mm = m; E gn = pwi(g, qm1);
            for (u128 p = 2; p * p <= mm; p++) if (mm % p == 0) { if (eis1(pwi(gn, m / p))) ok = 0; while (mm % p == 0) mm /= p; }
            if (mm > 1 && eis1(pwi(gn, m / mm))) ok = 0; if (ok) break; }
        E x = g;                           /* the class of g has order q^2+q+1, so g^i, 0 < i < N, are the other points */
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
    (void)q;
    const char *names[6] = {"R", "H1", "H2", "H3", "H4", "H5"}; int unc[6] = {0}, pre[6] = {0}, fails = 0, badpre = 0;
    for (int cs = 0; cs < 6; cs++) for (int it = 0; it < samples; it++) {
        col c[3];
        if (cs == 0) {
            int k[3]; E x[3];
            for (int j = 0; j < 3; j++) { k[j] = rnd64() & 1; x[j] = rnd_el(); }
            if ((k[0] == k[1] && same_point(x[0], x[1])) || (k[0] == k[2] && same_point(x[0], x[2])) || (k[1] == k[2] && same_point(x[1], x[2]))) continue;
            for (int j = 0; j < 3; j++) c[j] = column(k[j], x[j]);
        } else if (cs == 1 || cs == 2) {
            E x1 = rnd_el(), x2 = rnd_el(); if (same_point(x1, x2)) continue;
            E s1 = rnd_scalar(1); E s2 = rnd_scalar(1);                 /* drawn in this order (v1/v2: left operand first) */
            E x3 = exor(mul(s1, x1), mul(s2, x2));
            c[0] = column(cs == 2, x1); c[1] = column(cs == 2, x2); c[2] = column(cs == 2, x3);
        } else if (cs == 3) {
            E y = rnd_el(), w1, w2; two_w(&w1, &w2);
            c[0] = column(0, mul(y, w1)); c[1] = column(0, mul(y, w2)); c[2] = column(1, y);
        } else if (cs == 4) {
            E x = rnd_el(), w1, w2; two_w(&w1, &w2);
            c[0] = column(0, x); c[1] = column(1, mul(mul(x, binv), w1)); c[2] = column(1, mul(mul(x, binv), w2));
        } else {
            E x = rnd_el(), x3 = rnd_el(); int k3 = rnd64() & 1; if (same_point(x, x3)) continue;
            c[0] = column(0, x); c[1] = column(1, x); c[2] = column(k3, x3);
        }
        int ok[3]; for (int r = 0; r < 3; r++) ok[r] = indep(c[0].r[r], c[1].r[r], c[2].r[r]);
        if (cs >= 1 && cs <= 4 && (ok[0] || ok[2])) pre[cs]++;
        if (!ok[0] && !ok[1] && !ok[2]) unc[cs]++;
    }
    for (int cs = 0; cs < 6; cs++) { fails += unc[cs]; if (!*T) badpre += pre[cs]; }
    int nwF = (n + 64) / 64;
    printf("kcheck_B s %d n %d F ", s, n);
    if (nwF <= 2) printf("%llx:%016llx", F.w[1], F.w[0]);
    else { printf("%llx", F.w[nwF - 1]); for (int k = nwF - 2; k >= 0; k--) printf(":%016llx", F.w[k]); }
    printf(" samples %d seed %s tamper %s | uncovered/precondition-broken", samples, argv[3], *T ? T : "-");
    for (int cs = 0; cs < 6; cs++) printf(" %s:%d/%d", names[cs], unc[cs], pre[cs]);
    printf(" | %s\n", (fails == 0 && badpre == 0) ? "PASS" : "FAIL");
    return (fails == 0 && badpre == 0) ? 0 : 1;
}
