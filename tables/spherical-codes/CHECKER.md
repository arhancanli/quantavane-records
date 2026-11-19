# Checker: spherical codes

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE d N C
```

## What it checks

```text
Exact checker for spherical codes: N points on the unit sphere in R^d whose largest pairwise cosine is at most C.

    python3 checker.py FILE d N C

FILE: one point per line, d decimal coordinates separated by commas (other lines are ignored, as on
spherical-codes.org). Points need not be normalised: the cosine of x and y is <x, y> / (|x| |y|). C is the claimed
largest cosine as a decimal with E digits after the point; the claim holds when the true largest cosine is below
C + 10^-E, that is, when C is the true value rounded to its last stated digit or above it. The table keeps, for each
(d, N), the smallest largest cosine.

  A  exact: every coordinate is an exact decimal, so after scaling by a power of ten each point is an integer vector;
     for each pair it decides <x, y> / (|x| |y|) < C + 10^-E in integers, comparing squares after the signs (for
     C + 10^-E > 0: <x, y> <= 0 or <x, y>^2 10^(2E) < (C 10^E + 1)^2 |x|^2 |y|^2), and reports the largest cosine to
     30 digits.
  B  60-digit decimal arithmetic, written separately: normalises every point, takes every inner product and compares
     the largest with C + 10^-E.

Both also check that there are N points of dimension d, none zero and no two equal. Exit 0 = both PASS.
MIT licence (LICENSE-CODE).
```

## Other checkers

`checker/sph_check.py` in github.com/arhancanli/quantavane-spherical-codes (50-digit decimal arithmetic) is a third checker.

## Limits

Quadratic in N; under a second for a few hundred points.

## Runs on our own certificates (19 November 2026)

- All 443 spherical codes here, from the files in github.com/arhancanli/quantavane-spherical-codes against the largest cosine in records.csv: A and B PASS (each under 1.1 s, 2 minutes in all).
- The regular tetrahedron passes with C = -0.33333333333333333333 and fails with C = -0.35; the square passes with C = 0 and fails with C = -0.1.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
