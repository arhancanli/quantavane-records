# Checker: linear codes

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE q n k d [--modulus c0,c1,...] [--exact]
```

## What it checks

```text
Exact checker for linear codes: an [n, k, d]_q code (minimum distance at least d).

    python3 checker.py FILE q n k d [--modulus c0,c1,...] [--exact]

FILE: a generator matrix, k lines of n field elements. For q prime the elements are 0..q-1; for q = p^e they are
integers c0 + p c1 + p^2 c2 + ... giving the coefficients in the polynomial basis of the Conway polynomial (built in
for q = 4, 8, 9, 16, 25, 27, 32, 49, 64; --modulus gives the low coefficients of another monic polynomial). The table
(codetables.de, MinT, Chen's database) keeps the largest d for each (q, n, k). Two separately written routines:
  A  every one of the q^k messages, in odometer order, the codeword updated by adding the change in one coordinate
     times its generator row (field tables by shift-and-add); a nonzero message giving the zero word means the rows
     are dependent and fails.
  B  one message per projective point (first nonzero coordinate 1), each codeword computed from scratch (field
     tables by schoolbook products); scalar multiples have equal weight, so this covers every nonzero codeword.
Both report the minimum weight; the check passes when both find it at least d (--exact: equal to d).
Exit 0 = both PASS. Pure Python: practical up to q^k of a few million. MIT licence (LICENSE-CODE).
```

## Limits

Both routines enumerate codewords: practical up to q^k of about ten million (a minute or two). Larger codes need a compiled distance program.

## Runs on our own certificates (19 November 2026)

- [100,4,89]_16 (a MinT cell): A PASS on all 65,535 codewords, B PASS on 4,369 projective ones, minimum weight 89 (under a second).
- [87,8,63]_7 (listed in Chen's database; sent to MinT): A PASS on 5,764,800 codewords, B on 960,800, minimum weight 63 (90 s).
- The Hamming code [7,4,3]_2 and the hexacode [6,3,4]_4 pass; claiming d = 4 for the Hamming code fails.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
