# Checker: circulant weighing matrices

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE n k [--proper]
```

## What it checks

```text
Exact checker for circulant weighing matrices CW(n, k).

A CW(n, k) is a circulant n x n matrix W with entries in {0, 1, -1} and W W^T = k I. It is proper when its first row
is not x^s w(x^d) for a CW(n/d, k) row w and some d > 1 dividing n, that is, when for every prime d dividing n the
support is not inside one residue class mod d. The table (Gordon's La Jolla repository) records which (n, k) exist.

    python3 checker.py FILE n k [--proper]

FILE: the first row (n integers in {0, 1, -1} on one line), or the whole matrix (n such lines).
  A  builds the n x n circulant (or reads the matrix and checks it is circulant) and computes W W^T entry by entry.
  B  written separately: the periodic autocorrelation of the first row, sum_i w_i w_(i+s), is 0 for s = 1..n-1 and
     k for s = 0, computed as the product w(x) w(x^-1) in Z[x]/(x^n - 1).
Both check the entries and that exactly k entries of the row are nonzero. --proper also checks properness.
Exit 0 = both PASS. MIT licence (LICENSE-CODE).
```

## Limits

Quadratic in n; instant for n in the hundreds.

## Runs on our own certificates (19 November 2026)

- CW(146, 64) and CW(242, 81), both proper: A, B and the properness test PASS.
- The lab's damaged CW(146, 64) row (one -1 moved): A and B FAIL.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
