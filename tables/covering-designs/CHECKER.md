# Checker: covering designs

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE v k t [m] [--blocks b]
```

## What it checks

```text
Exact checker for covering designs C(v, k, t) and C(v, k, t, m).

A C(v, k, t, m) covering is a set of k-subsets (blocks) of {1, ..., v} such that every m-subset meets some block in at
least t points; C(v, k, t) is the case m = t (every t-subset lies in a block). The table keeps the fewest blocks.

    python3 checker.py FILE v k t [m] [--blocks B]

FILE: one block per line, points 1..v separated by spaces. The program runs two separately written routines and
passes only when both pass:

  A  m = t: generates the t-subsets inside every block and counts them against C(v, t).
     m > t: a point-by-point search for an m-subset meeting every block in at most t - 1 points (a counterexample);
     points with the same blocks form runs entered in order, and a run's contribution is bounded by its spare
     capacity. None may exist.
  B  written separately: a search over how many points of each type (the set of blocks holding a point) a
     counterexample would use, with a counting bound; for m = t as well. None may exist.
  github.com/arhancanli/quantavane-covering-records has verify.py, a third checker of the same kind as B.

Both also check the file: every block has k distinct points in 1..v, and (with --blocks) the number of blocks.
Exit status 0 = both PASS, 1 = a FAIL, 2 = bad input. MIT licence (LICENSE-CODE).
```

## Other checkers

`verify.py` in github.com/arhancanli/quantavane-covering-records is a third, separately written checker (an exact search for an m-subset that meets every block in fewer than t points).

## Limits

m = t is fast. For m > t both routines search: most coverings here take seconds, six take more than a minute and the hardest 17 minutes.

## Runs on our own certificates (19 November 2026)

- All 141 coverings in github.com/arhancanli/quantavane-covering-records (the 132 records here and 9 more files there): A and B PASS, 51 minutes in all; 135 took under a minute each, the longest C(36,9,5,11) <= 68 (17 minutes) and C(31,10,8,19) <= 37 (12 minutes).
- The Fano plane passes as C(7,3,2) and C(7,3,2,3); with one line removed it fails (tests/test_checkers.py).
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
