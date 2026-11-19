# Theorem K over a declared range: q = 4^s, 3 <= s <= 22 (18 Nov 2026; s = 17-22 added 20 Nov 2026, sections 8-9)

Quantavane, 18 Nov 2026; published 19 Nov 2026 with the records dataset. Theorem K itself (the array and a sketch of
its proof) is public in github.com/arhancanli/quantavane-covering-arrays, folder `three-row/`; a note with the full
proof is in preparation. This note compares Theorem K with every published construction class for every s in the
declared range: a cell counts only where the theorem's array has fewer rows than all of them. The analysis was
refereed inside Quantavane on 18 Nov 2026 (section 7) and has not been peer reviewed. The two count programs,
range_cells_A.py and range_cells_B.c, and their outputs are in this folder. Since 20 Nov 2026 so are the two array
checkers of section 4 (kcheck_A.py, kcheck_B.c), a fourth count check (classcheck_C.py) and every log they wrote
(logs/); `qvrecords` in this repository carries its own independent check of any cell.

## 1. What is counted

Theorem K: for q = 4^s, CAN(3, k, q) <= N_K = 3q^3 - 2 for every k <= K = 2(q^2+q+1). A cell (t = 3, k, v = q) counts
when every construction in sections 2-3 has more rows than N_K at width k, or reaches fewer than k columns. A
construction with exactly N_K rows removes the cell (a tie is not a record). The cells are k in (R(s), K], where R(s) is
the largest reach of any construction with at most N_K rows. Both programs compute R(s) and the cells exactly for each s.

## 2. The comparison, class by class (proved bounds; the literature inputs they need are in section 3)

**Class I: covering arrays from a CPHF with at most two rows, over any prime power Q >= q, symbols fused to q.**
This class contains RMS (Raaphorst-Moura-Stevens), the affine and Sherwood restrictions (Colbourn-Ingalls-Jedwab-
Saaltink-Smith-Stevens 2022; Torres-Jimenez-Izquierdo-Marquez), the orthogoval planes, the LFSR arrays with l = 2
(Tzanakis-Moura-Panario-Stevens) and every 2-row CPHF anyone may build.

*Lemma R1.* A CPHF(2; k, Q, 3) has k <= (Q^2+Q+1) a(Q) / (Q+1), where a(Q) = Q+1 for odd Q and Q+2 for even Q. So
k <= Q^2+Q+1 for odd Q (RMS is optimal) and k <= (Q+1)^2 for even Q.
*Proof.* Fix a line L of PG(2,Q) and let C_L be the columns whose row-1 entry lies on L (a zero entry lies on every
line). Any three columns of C_L fail in row 1, so their row-2 entries are three distinct non-collinear points: the
row-2 entries of C_L form an arc, so |C_L| <= a(Q) (also when |C_L| <= 2). Each column lies in C_L for at least Q+1
lines L. Summing over the Q^2+Q+1 lines, k(Q+1) <= (Q^2+Q+1) a(Q). [] A Sherwood CPHF is a CPHF whose entries are
affine points (1, x, y), with the same covering condition, so the lemma covers it.

The fewest rows such an array can have is the Sherwood count with one Colbourn-Lanus restricted pair (the two rows
share a coordinate; Q^2 - Q rows saved), 2(Q^3 - Q) + Q - (Q^2 - Q)
= 2Q^3 - Q^2. Fusing Q - q symbols saves 2 rows each by our convention: 2Q^3 - Q^2 - 2(Q - q). (The plain Sherwood
count 2Q^3 - Q and the non-Sherwood 2Q^3 - 1 are larger; using either gives the same q1 for every s <= 20, checked.) This increases with Q, so the class competes exactly for the prime powers Q <= q1(s), the
largest prime power with 2Q^3 - Q^2 - 2(Q - q) <= N_K. Its reach is R_I = max(R1(q1), R1(largest power of 2 in [q, q1])).
*Fusion margin:* the next prime power q2 > q1 needs N_K + m rows, with m = 193,377 at s = 3, 7,442,889 at s = 4 and
5.5e20 at s = 16 (column `margin`). A published fusion rule would have to save m / (q2 - q) rows per fused symbol
(12,891 at s = 3, about 3.15 q^2; relative to q it shrinks, to about 99 q at s = 7 and 207 q at s = 16, never below
about 99 q in the range; at s = 17, 18, 19, 20 it is about 3,298 q, 561 q, 1,129 q and 579 q, classcheck_C.py) to move
q1, against the 2 per symbol of the convention.

**Class II: 3-row constructions.** Over Q > q a 3-row CPHF array has at least 3(Q^3 - Q) + Q - 3(Q^2 - Q) - 2(Q - q)
= 3Q^3 - 3Q^2 - Q + 2q rows, even with all three pairs restricted and symbols fused. This exceeds N_K: 3(Q^3 - q^3) =
3(Q - q)(Q^2 + Qq + q^2) >= 3Q^2 + 3Qq > 3Q^2 + Q - 2q - 2. Over GF(q) itself:
- *IIa, compositions.* PHF(s'; k, v, 3) composed with CPHF(n0; v, q, 3), s' n0 <= 3. With n0 = 1 the inner array is an
  arc, v <= q + 2. Shangguan-Ge 2016 (arXiv 1601.04807, Theorem 1.5 with u = 3: n <= r v^ceil(N/2) +
  (2 - r) v^floor(N/2)) gives k <= v^2 + v for three rows and k <= 2v for two. Reach <= (q+2)^2 + (q+2).
- *IIb, the random oversampling bound with three rows* (Colbourn-Lanus-Sarkar 2018 as restated in Colbourn, IJGT 13(3)
  2024, Section 3; with or without a restricted pair). The bound gives k0 - floor(C(k0,3) f) columns, f the probability
  that a fixed triple is uncovered: per row (q^3 + 3q^2 + 2q + 1)/(q^2+q+1)^2 >= 1/(q+2) on PG(2,q),
  (q^3 + q^2 - q)/q^4 >= 1/q on AG(2,q), and (2q^2 - 1)/q^4 >= 1/q^2 for a restricted pair, so f >= (q+2)^-3. Since
  C(k0,3) >= (k0-2)^3/6, the bound is at most 2 + (2/3) sqrt(2/f) <= 2 + 0.943 (q+2)^1.5 < 4 q^1.5 + 3. This is a bound on
  the *formula* only (referee, 18 Nov): the search methods built on it do better, see IIc.
- *IIc, published search methods for 3-row CPHFs* (CLS's conditional-expectation algorithm, Torres-Jimenez-Izquierdo-
  Marquez simulated annealing, Wagner-Colbourn-Simos, our own partly cyclic searches). At q <= 25 they reach 1.1-1.9 q^2
  columns (q = 16: 340 = 1.33 q^2, CLS's conditional expectation 313 against 259 from the IIb formula). No array from
  them is published for q = 4^s >= 64. That is a fact about the literature, not a theorem: input P2 in section 3.
- *IId, Bergeron's RMS extension* (2026, github.com/MichelBacon/covering-arrays-rms-extension; in ca-tools):
  CA(2q^3 + q^2 L; 3, q^2 + q + 1 + m, q), L = max(q - 2, ceil(m(q - 2)/q)). At most N_K rows forces L <= q - 1, so
  m <= q(q - 1)/(q - 2), m <= q + 1: reach <= q^2 + 2q + 2 (4,226 at s = 3). Over Q > q it needs >= 3Q^3 - 2Q^2 > N_K.
- *IIe, add-a-factor and Colbourn-Martirosyan-Trung-Walker on the hyperoval orthogonal array:* fewer than N_K rows
  but reach at most about q^2 + 5q (referee's ref_cells.py evaluates both per s).
- *Others:* Theorem K and its corollary K' are ours.

**Class III: plications and PHF o RMS.** Shokri-Moura's 2-plication has 4q^3 - 5q^2 + 2q rows. Larger plications and
PHF o RMS (6q^3 - 5) have more. All exceed N_K for q >= 5 (q^3 - 5q^2 + 2q + 2 > 0). The programs assert it per s.

**Class IV: CPHF arrays with four or more rows.** At least 4(q^3 - q) + q - 6(q^2 - q) = 4q^3 - 6q^2 + 3q rows, even with
all six pairs restricted (Colbourn-Lanus-Sarkar's path-restricted 4-row SCPHF has 4q^3 - 3q^2; ca-tools lists 15,616
at q = 16, k = 512). This exceeds N_K for q >= 6. (Corrected 18 Nov: the 17 Nov text gave 4q^3 - 2q^2 - q, which is
not a floor; the conclusion is unchanged.)

## 3. What the count rests on besides proofs: three literature inputs

Classes I, IIa, IIb (as a formula), IId, IIe, III and IV are proved bounds. The count also needs three statements about
the published literature. They are disclosed here, as the counting rule asks for a comparison with every published
construction, and that can only be checked against the literature as it stands.

**P (Roux-type recursions and column replacement).** Chateauneuf-Kreher doubling CAN(3,2w,q) <= A + (q-1)B,
Colbourn-Martirosyan-Trung-Walker CAN(3,qw,q) <= A + (q-1)B + q^3 - q^2, Cohen-Colbourn-Ling (with l = 2 a difference
matrix has at most q columns, so it would need A < 2q^3 - 2q^2) and Shokri's Theorem 5.2.6 (x = 2 reduces to the
doubling condition; x >= 3 has more rows than N_K) all need a strength-3 array A and a strength-2 array B of width
w >= (R+1)/x. To reach past R with at most N_K rows they need A + (q-1)B <= N_K. Two-row column replacement (a 2-row
PHF composed with covering arrays; Colbourn-Torres-Jimenez heterogeneous PHFs) reaches at most k1 + k2
(Shangguan-Ge Lemma 3.1), so it needs a strength-3 ingredient with about 1.5 q^3 rows at width >= 0.65 q^2. Nothing here
is excluded by a lower bound we can prove. A counting bound (pairs of rows colliding in two columns, evaluated exactly
at q = 16, 64, 256, 1024 and w = q^2/2, 0.65 q^2, q^2) gives only B >= q^2 + q - 1, so a strength-2 array with
q^2 + q + 2 rows is not excluded. P says that no published construction supplies such ingredients for q = 4^s >= 64.
Evidence: the published general constructions give B = 2q^2 - q (products of orthogonal arrays) at w <= (q+1)^2, and
A >= 2q^3 - q at w >= (R+1)/2, where only CPHF types with two rows exist (class I). Then A + (q-1)B >= 4q^3 - 3q^2, which
is q^3 - 3q^2 + 2 above N_K. In every alphabet the tables cover (prime powers 7 to 25, ca-tools with every published
construction and every search entered), the smallest A at w >= 0.65 q^2 is at least 1.98 q^3 (8,176 at q = 16, against
T3 = q^3 + 3q^2 - q - 2 = 4,846). The smallest B there is exactly 2q^2 - q, and the smallest A at w >= 1.31 q (the CMTW
base) is 1.40-1.96 q^3. The programs print T2 = q^2 + q + 2, the most rows B could have with A = 2q^3 - q, and T3, the
most rows A could have with B = 2q^2 - q. The range referee searched for a strength-2 array with about q^2 + q rows at
width Theta(q^2) and found none.

**P2 (no search-found 3-row CPHF at these alphabets).** No 3-row CPHF over GF(4^s), s >= 3, with more than R(s) columns
has been published. This is the most fragile input. At q <= 25 the search methods of IIc reach 1.1-1.9 q^2, and R(s)
is about 1.31 q^2, so a 3-row search at q = 64 or 256, by anyone including us, could take away part of a q's range,
from its bottom up. The cells past R(s) exist because nobody has run those searches at q >= 64. Whether such a search is
feasible at q = 64 (about 8,000 columns) is untested here; one run would turn this input into a measured margin.

**P3 (the class list is complete).** Searched by us, by the two blind referees of 17 Nov and by the range
referee of 18 Nov. Ruled out by that referee: Tzanakis character-sum arrays (thesis Cor. 5.10; they need q' >= v^4, so
rows grow like v^9), cyclotomic arrays (with at most 3q^3 rows a random-like array leaves about e^-3 q^3 patterns
uncovered per column triple), alphabet products (2.4 million rows at s = 3, N_K = 786,430), Panario-Saaltink-Stevens-
Wevrick 2020 (class I material), and arXiv 2510.13122 (strength 4). Bergeron's extension, add-a-factor and two-row column
replacement were missing from an earlier version; they are in classes IId, IIe and P now, and none moves R(s).

## 4. The array, checked on samples at every s in the range (two programs)

The proof covers every s. Exhaustive checks exist for q = 4, 16, 64, 256, by kummer3.c and singer_check(2).c (17 Nov),
and by the referees at q = 4 and 16. For every s from 1 to 16, both new programs test triples with exact arithmetic,
in their own field models:
- kcheck_A.py: GF(q) = GF(2)[t]/(f), f by Rabin's test; c the first non-cube; GF(q^3) = GF(q)[B]/(B^3 - c); 3x3 determinants.
- kcheck_B.c: GF(2^(6s)) as one field; b = h^((q^2+q+1)/3) for the first non-cube h; the Moore determinant (no coordinates).
Cases: R, random triples. H1 and H2: three B (resp. E) columns with collinear x, so rows 1 and 3 fail (Lemma 1).
H3: B_{y w1}, B_{y w2}, E_y with w_i in <1, b>. H4: B_x, E_{x w1 / b}, E_{x w2 / b}. In H3 and H4 rows 1 and 3 fail, and
row 2 must cover by Lemma 2. H5: B_x, E_x and any third column. The programs confirm that the failing rows really fail
(precondition count 0), so the hard cases test the lemmas, not luck. A: 3,000 samples per case and s. B: 10,000.
Tampers, each at s = 3, 8 and 16 in both programs:
`inv` (row 2 of E is b x, not b x^-1) must FAIL. `bq` (b in GF(q)) must FAIL. `mod` (a reducible modulus) must be REFUSED.
`cube` (c, resp. h, a cube) must be REFUSED. All 32 clean runs PASS; all 24 tamper runs FAIL or are REFUSED
(logs/samples_20261118.log). For s = 17-20 see section 8: 8 more clean runs PASS and 32 more tamper runs FAIL or are
REFUSED.

*Exhaustive Lemma 1/2 cases (added after the range referee, 18 Nov).* The array is invariant under x -> yx (rows 1 and 3
multiplied by y, row 2 by y^-1), so every H3 triple is equivalent to one with y = 1 and every H4 triple to one with
x = 1. `kcheck_B S 0 SEED` checks all of them: H3 and H4 over all pairs of the q + 1 points of <1, b>, and H5 = (B_1,
E_1, any third column) over all 2(q^2+q+1) - 2 columns. The referee's program ref_kcheck.py (a third field model:
GF(2^(6s)) from a primitive polynomial, independence by GF(2)-rank) did the same at s = 3, 4 with identical counts.
Results (logs/exhaustive_lemma2_20261118.log; precondition broken 0 in every row):

| s | q | H3 triples | H4 triples | H5 triples | uncovered |
|---|---|---:|---:|---:|---|
| 1 | 4 | 10 | 10 | 40 | 0 |
| 2 | 16 | 136 | 136 | 544 | 0 |
| 3 | 64 | 2,080 | 2,080 | 8,320 | 0 |
| 4 | 256 | 32,896 | 32,896 | 131,584 | 0 |
| 5 | 1,024 | 524,800 | 524,800 | 2,099,200 | 0 |
| 6 | 4,096 | 8,390,656 | 8,390,656 | 33,562,624 | 0 |
A tamper that tests Lemma 2 itself, `nuc` (the referee's): the row-2 twist b becomes 1 + b, a point off the conic that is
not its nucleus. H3 and H4 then fail on exactly q/2 pairs (exhaustive: 32, 128 and 512 at q = 64, 256 and 1,024,
as the referee found at q <= 256; `inv` fails on every H4 triple), which is
where Lemma 2 fails. Random sampling catches such a failure only about 1/(q+1) of the time, so beyond s of about 5 the
samples add little and the proof carries the hard cases. That is why the exhaustive checks stop at s = 6 and the
sampled ones are a sanity check, not the evidence.

## 5. Counts (both programs agree byte for byte; out_A.txt = out_B.txt up to s = 20)

| s | q | q1 (class I) | R(s) | cells k in (R, K] | cumulative |
|---|---|---:|---:|---:|---:|
| 3 | 64 | 73 | 5,403 | 2,919 | 2,919 |
| 4 | 256 | 293 | 86,143 | 45,443 | 48,362 |
| 5 | 1,024 | 1,171 | 1,372,413 | 726,789 | 775,151 |
| 6 | 4,096 | 4,679 | 21,897,721 | 11,664,905 | 12,440,056 |
| 7 | 16,384 | 18,749 | 351,543,751 | 185,359,931 | 197,799,987 |
| 8 | 65,536 | 75,017 | 5,627,625,307 | 2,962,440,359 | 3,160,240,346 |
| 9 | 262,144 | 300,073 | 90,044,105,403 | 47,395,372,359 | 50,555,612,705 |
| 10 | 1,048,576 | 1,200,313 | 1,440,752,498,283 | 758,272,854,423 | 808,828,467,128 |
| 11 | 4,194,304 | 4,801,267 | 23,052,169,606,557 | 12,132,210,870,885 | 12,941,039,338,013 |
| 12 | 16,777,216 | 19,205,117 | 368,836,538,188,807 | 194,113,448,786,939 | 207,054,488,124,952 |
| 13 | 67,108,864 | 76,820,453 | 5,901,382,075,945,663 | 3,105,817,313,013,059 | 3,312,871,801,138,011 |
| 14 | 268,435,456 | 307,281,859 | 94,422,141,177,777,741 | 49,693,047,434,949,045 | 53,005,919,236,087,056 |
| 15 | 1,073,741,824 | 1,229,127,547 | 1,510,754,528,023,364,757 | 795,088,483,337,812,845 | 848,094,402,573,899,901 |
| 16 | 4,294,967,296 | 4,916,510,197 | 24,172,072,522,121,489,007 | 12,721,415,633,887,548,819 | 13,569,510,036,461,448,720 |
| 17 | 17,179,869,184 | 19,666,040,921 | 386,753,165,526,112,569,163 | 203,542,644,866,952,820,919 | 217,112,154,903,414,269,639 |
| 18 | 68,719,476,736 | 78,664,163,737 | 6,188,050,656,520,209,968,907 | 3,256,682,309,356,519,411,959 | 3,473,794,464,259,933,681,598 |
| 19 | 274,877,906,944 | 314,656,655,023 | 99,008,810,550,577,887,785,553 | 52,106,916,901,800,514,866,609 | 55,580,711,366,060,448,548,207 |
| 20 | 1,099,511,627,776 | 1,258,626,620,143 | 1,584,140,968,933,850,239,960,593 | 833,710,670,297,607,132,707,313 | 889,291,381,663,667,581,255,520 |

R(s) is always the class I reach R1(q1) = q1^2 + q1 + 1 (q1 odd), so the counts equal those of the 17 Nov count
for s <= 10. Lemma R1 now proves that bound for every 2-row CPHF, not only for RMS. Classes IIa, IIb, IId and IIe reach
at most (q+2)^2 + q + 2, 4 q^1.5 + 3, q^2 + 2q + 2 and about q^2 + 5q, all below R(s) for every s >= 3. Inside the tables (q = 16) the count
is the 168 cells the keeper lists.

**Declared range s <= 20 (20 Nov 2026): 889,291,381,663,667,581,255,520 cells (8.9e23), plus 168 at q = 16.** Until
20 Nov the declared range stopped at s = 16 (13,569,510,036,461,448,720 cells, 1.36e19) because the second array
checker, kcheck_B.c, refused s > 16; section 8 records the extension to s = 20 and section 9 the extension to s = 22
(2.28e26). Step by step: s <= 17 gives 2.17e20, s <= 18 3.47e21, s <= 19 5.56e22, s <= 21 1.42e25.

## 6. How to read this number

It is one theorem. Each q adds about (2 - 1.5^(2/3)) q^2 = 0.69 q^2 cells, so the total is set by where the declared range
stops. s <= 15 gives 8.5e17, s <= 16 gives 1.4e19 and s <= 20 gives 8.9e23. No table lists any alphabet from 64 upward at strength 3. Under
the rule's wording each cell is a record: the proved array beats every construction in sections 2-3 at that cell.
The count also rests on the literature inputs P, P2 and P3 of section 3. P2 is the weakest of them: it says that
nobody has published a 3-row search at q >= 64, which is a fact about the literature, not a theorem, and one published
3-row search at q = 64 or 256 could remove the bottom of those ranges. They are not 8.9e23 separate discoveries, and
public wording should say "one theorem settles 8.9e23 cells past every published construction" rather than
"8.9e23 records".

## 7. Referee inside Quantavane (18 Nov 2026): ACCEPT WITH FIXES

A referee working from a fresh start read this note, the theorem's write-up and proof, and the programs, and wrote its own count (ref_cells.py: prime sieve,
Newton cube root, extra classes for Bergeron, add-a-factor and CMTW). That program matches out_A.txt exactly on q1, R,
cells and totals for s = 3..20, under both class I floors. It also wrote its own exhaustive Lemma 2 check (ref_kcheck.py,
s = 3, 4 PASS) and ran the tampers, all of which FAIL or are REFUSED. Its findings, and what this version did with each:
1. The class IV floor was wrong (path-restricted 4-row SCPHF: 4q^3 - 3q^2). Fixed: 4q^3 - 6q^2 + 3q with all pairs
   restricted, still > N_K; both programs assert it.
2. Class II over Q > q ignored restricted pairs and had a wrong inequality step. Fixed.
3. IIb was called proved, but it bounds only the oversampling formula, and the published search methods do better
   (1.1-1.9 q^2 at q <= 25). Fixed: IIb is restated as a bound on the formula, IIc is added, and P2 is stated.
4. The class list was not exhaustive. Bergeron's RMS extension (reach q^2 + 2q + 2), add-a-factor and two-row column
   replacement are added (IId, IIe, P). None moves R(s).
5. The fusion-margin wording was wrong. Corrected (per symbol at least about 99 q in the range).
6. The figure "about 1.9 q^3" was wrong; it is 1.40-1.96 q^3. Corrected.
7. The note said only fused RMS arrays use fewer than 3q^3 - 2 rows past q^2 + q + 1. Corrected in the note.
8. Section 4 now has the nuc tamper and exhaustive Lemma 2 checks.
Its verdict on the number: **the count holds**. No published construction found moves R(s). The count rests on the
proofs plus the literature inputs P, P2 and P3 of section 3, and should be quoted with them.


## 8. Extension of the declared range to s = 17, 18, 19, 20 (20 Nov 2026)

The declared range had stopped at s = 16 for one reason: the second array checker, kcheck_B.c, refused s > 16 (its
field GF(2^(6s)) sits in one unsigned 128-bit word). Nothing else here depends on s <= 16: the proof of Theorem K
covers every s, the class comparison of sections 2-3 is analytic in q, and the count programs already ran to s = 20.

1. **kcheck_B.c version 2** (sha256 efbda2ed...). The field arithmetic is unchanged: its multiplication reduces after
   every shift, so it is exact for n = 6s <= 127; the cap is now s <= 20. New: every run first compares that
   multiplication with a second one written apart (the full 256-bit carry-less product in two 128-bit words, then long
   division) on 4,016 pairs from its own random generator, and refuses on any mismatch; the exhaustive mode is refused
   above s = 8. Two comments in the source are stale and were left so that the logged hash stays valid: line 3 still
   says "6s <= 96", and the exhaustive guard's comment says the tables overflow, which holds from s = 16.
2. **Byte-for-byte reproduction at s <= 16.** Version 2 re-ran every command behind logs/samples_20261118.log and
   logs/exhaustive_lemma2_20261118.log (written by the earlier versions, sha256 cc7e94e6... and 4fecca10...). Every
   output line is identical apart from the START/END lines and the hash line (logs/samples_v2repro_20261119.log,
   logs/exhaustive_v2repro_20261119.log).
3. **Red controls of the new cross-check** (red_controls/20261119/red1.c, red2.c; logs/mulcheck_red_controls_20261119.log): a
   multiplication cut to a 100-bit word is harmless at s = 16 and refused at s = 17 and 20; a second multiplication with
   its high cross term dropped is refused from s = 11.
4. **Runs at s = 17-20** (logs/samples_s17_20261119.log ... samples_s20_20261119.log): per s, kcheck_B.c with 10,000
   samples per case and kcheck_A.py with 3,000, cases R and H1-H5. All 8 clean runs PASS with precondition-broken 0 in
   every case; the tampers inv and bq FAIL in both programs and mod and cube are REFUSED in both (32 of 32).
   Commands: `cc -O2 -o kcheck_B kcheck_B.c && ./kcheck_B 17 10000 1017` and `python3 -I kcheck_A.py 17 3000 2017`.
5. **The comparison at s = 17-20.** Every statement of sections 2-3 is either analytic in q or a number the count
   programs give to s = 20. The fusion margin per fused symbol is about 3,298 q, 561 q, 1,129 q and 579 q at s = 17-20,
   above the 99 q floor of section 2. A fourth count program written apart, classcheck_C.py, re-derives q1, the next
   prime power, the margin, R(s), the cells and the totals for s = 3..20 and matches out_A.txt on every line
   (logs/classcheck_C_20261120.txt).
6. **A fresh referee inside Quantavane** audited the checker change and these runs only (20 Nov 2026) and accepted
   them. It rebuilt both versions and found their output identical on 108 commands at s <= 16, and wrote its own field
   model (irreducibility by Ben-Or's and Berlekamp's tests, independence by a GF(2) rank): at each s = 17-20 the
   modulus is irreducible, the multiplication, inverses and Frobenius agree, and a trace of 600 triples per s from the
   real case loop agrees with its rank test. It also ran the 18 Nov referee's third field model at s = 17 with its own
   primitive polynomial: PASS. It notes that the hard sub-case of H5 is reached with probability about 1/q at large s,
   so there the sampled evidence is H1-H4, and the proof carries the rest, as section 4 says.

**Result.** Theorem K cells past every published construction: s <= 17: 217,112,154,903,414,269,639 (2.17e20);
s <= 18: 3,473,794,464,259,933,681,598 (3.47e21); s <= 19: 55,580,711,366,060,448,548,207 (5.56e22); s <= 20:
889,291,381,663,667,581,255,520 (8.9e23). With the 168 keeper-listed cells at q = 16: 889,291,381,663,667,581,255,688.
One theorem settles 8.9e23 cells past every published construction. The count rests on the proofs and on the
literature inputs P, P2 and P3; P2, that nobody has published a search-found 3-row array at q >= 64, is the weakest,
and it is not a theorem.

## 9. Extension to s = 21 and 22 (20 Nov 2026)

Two programs had 128-bit limits: kcheck_B.c version 2 (exact for 6s <= 127, so s <= 21) and range_cells_B.c version 1
(4q^3 overflows 128 bits from s = 21). Both were rewritten for wider numbers; the mathematics and the output are
unchanged. The previous versions are kept here as kcheck_B_v2_20261119.c and range_cells_B_v1_20261118.c.

1. **range_cells_B.c version 2**: 256-bit integers (four words, every product checked for overflow). Its output for
   s <= 20 equals version 1's byte for byte, and for s <= 22 it equals range_cells_A.py's (out_A.txt = out_B.txt, now to
   s = 22, with their first 18 lines unchanged). classcheck_C.py agrees on every line to s = 22
   (logs/classcheck_C_s22_20261120.txt).
2. **kcheck_B.c version 3**: field elements, polynomials and exponents in four 64-bit words (s <= 30). The main
   multiplication is now the full carry-less product followed by long division; at start-up it is compared with a
   shift-and-reduce multiplication on 4,024 pairs, and z^(2^n) = z is tested on 8 elements.
3. **Byte-for-byte reproduction**: version 3 reproduces every earlier kcheck_B log, s <= 16 and the tampers, the
   exhaustive runs for s <= 6, and s = 17-20 (logs/v3repro_*_20261120.log). The last line of the s = 17-20 comparison
   was written by hand after a file-list fix in the comparison script, and says so.
4. **Red controls** (red_controls/20261120/; logs/mulcheck_red_controls_v3_20261120.log): the main multiplication cut to
   a 256-bit product, and the cross-check cut to a 128-bit accumulator, are each harmless at s = 21 and refused at s = 22.
5. **Runs at s = 21 and 22** (logs/samples_s21_20261120.log, logs/samples_s22_20261120.log): kcheck_B.c with 10,000
   samples per case and kcheck_A.py with 3,000; all 4 clean runs PASS with precondition-broken 0; the tampers inv and bq
   FAIL and mod and cube are REFUSED in both programs (16 of 16).
6. **The comparison at s = 21 and 22**: every statement of sections 2-3 is analytic in q or is a number the count programs
   give to s = 22. R(s) = q1^2 + q1 + 1 with q1 odd; the fusion margin per fused symbol is about 766 q and 1,732 q,
   above 99 q; q1 is the same under the three class I floors for every s <= 22 (altfloor_C.py,
   logs/altfloor_C_20261120.txt).
7. **A fresh referee inside Quantavane** accepted both claims: version 3 is a valid second program for the sampled
   checks at s = 21 and 22 and its runs pass, and the counts at s = 21 and 22 are right.
   - It tested version 3's arithmetic against its own field model at s = 1-30 with no disagreement.
   - It confirmed the moduli irreducible by Berlekamp's and Ben-Or's tests.
   - It traced 2,376 real case-loop samples per s against its own rank test, with no difference.
   - It ran the 18 Nov referee's third field model at s = 21 and 22: PASS.
   - It recounted every number of out_A.txt to s = 22 with its own program.

   Its main caveat: 10,000 samples per case see only defects that hit about 3e-4 or more of a case. A defect as rare as
   the nuc tamper's, about 1/(q+1) of H3 pairs, passes the sampled runs. Rare cases are carried by the exhaustive runs
   at s <= 6 and by the proof.

**Result.** Theorem K cells past every published construction: s <= 21: 14,228,662,105,437,402,687,707,131 (1.42e25);
s <= 22: 227,658,593,685,328,988,782,135,420 (2.28e26). With the 168 keeper-listed cells at q = 16:
227,658,593,685,328,988,782,135,588. One theorem settles 2.3e26 cells past every published construction. The count
rests on the proofs and on the literature inputs P, P2 and P3. P2, that nobody has published a search-found 3-row array
at q >= 64, is the weakest, and it is not a theorem.
