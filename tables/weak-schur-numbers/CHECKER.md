# Checker: weak Schur numbers

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE k n
```

## What it checks

```text
Exact checker for weak Schur numbers: WS(k) >= n.

A partition of {1, ..., n} into k sets is weakly sum-free when no set holds x, y, z with x + y = z and x != y.
WS(k) is the largest n for which such a partition exists; the table (Leaps in Bounds) keeps lower bounds.

    python3 checker.py FILE k n

FILE: k lines, each the numbers of one set separated by spaces or commas; or the colours of 1, 2, ..., n as one
string of base-36 digits 0..k-1 (line breaks ignored); or one line of n colours separated by spaces. The program checks the sets partition 1..n and runs two separately written routines:
  A  bitsets: for each set S as an integer M and each x in S, (M >> x) & M has a bit at y exactly when y and x + y
     are both in S; it must be empty apart from y = x.
  B  written separately: the sumset of each set from the square of the polynomial sum_{x in S} X^x, computed as one
     big-integer product with X = 2^24; the coefficient of X^z for z in S must be 0, after removing the pair x = y.
Exit 0 = both PASS. Both handle n in the millions (about a minute at n = 818,354). MIT licence (LICENSE-CODE).
```

## Limits

About two minutes at n = 818,354 (both routines).

## Runs on our own certificates (19 November 2026)

- WS(7) >= 2154, WS(8) >= 7002, WS(9) >= 22620 (Leaps in Bounds lists them) and WS(12) >= 818354: A and B PASS (WS(12): 123 s).
- WS(7) with 2 and 3 recoloured to the colour of 1: A and B FAIL (1 + 2 = 3).
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
