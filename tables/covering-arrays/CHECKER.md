# Checker: covering arrays

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py ca    FILE t v [--rows N] [--cols K]
python3 checker.py cphf  FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]
python3 checker.py scphf FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]
```

## What it checks

```text
Exact checker for covering arrays and Sherwood covering perfect hash families.

A covering array CA(N; t, k, v) is an N x k array over {0, ..., v-1} in which every N x t subarray contains every
t-tuple at least once. The table keeps, for each (t, k, v), the fewest rows CAN(t, k, v).

    python3 checker.py ca   FILE t v [--rows N] [--cols K]
    python3 checker.py cphf  FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]
    python3 checker.py scphf FILE t q [--rows N] [--cols K] [--modulus c0,c1,...]

ca:   FILE has N lines of k symbols separated by spaces.
  A  for every t columns, the set of t-tuples read down the rows has v^t elements.
  B  for every symbol of every column, the set of rows holding it as a bitmask; for every t columns and every
     t-tuple, the AND of the t masks must be nonzero.
cphf: FILE has n lines of k integers (the layout of the ca-tools table; lines not starting with a digit are skipped): each integer is a vector of GF(q)^t, written
     in base q with the most significant digit first, each digit a field element in the Conway coding
     (c0 + p c1 + p^2 c2 + ... for the polynomial basis). Every t columns must have a row whose t vectors are linearly
     independent over GF(q); then Sherwood's construction gives CA(n(q^t - 1) + 1; t, k, q).
  A  rank of the t x t matrix by Gaussian elimination.
  B  written separately: no nonzero combination sum c_i x_i (first nonzero c_i = 1, so (q^t - 1)/(q - 1) of them)
     vanishes; its field tables come from schoolbook products of the digit polynomials (A's from shift-and-add).
scphf: the Sherwood families of github.com/arhancanli/quantavane-covering-arrays: each integer is a point of
     AG(t-1, q) (t-1 digits base q, most significant first); every t columns must have a row whose points span
     AG(t-1, q) (the vectors (1, point) independent); the CA has n(q^t - q) + q rows. Same routines A and B.
--rows N (all modes) states the claimed number of rows.
--cols K checks only the first K columns (deleting columns keeps either property). Fields: q prime, or q in the
built-in Conway table (4, 8, 9, 16, 25, 27, 32, 49, 64), or any q = p^e with --modulus giving the low coefficients
of the monic polynomial. Exit 0 = both PASS. Pure Python: practical up to a few million t-sets; the large families
in github.com/arhancanli/quantavane-covering-arrays have C checkers in its checker/ folder.
MIT licence (LICENSE-CODE).
```

## Other checkers

For whole families of thousands of columns, use the compiled checkers of github.com/arhancanli/quantavane-covering-arrays (`checker/scphf_check.c`, `scphf4_check.c` to `scphf6_check.c`, `fullcheck3.c` and their fast variants), which decide every t-set of columns, and the ca-tools table's own `tools/cphf/cphfcheck.c`.

## Limits

Pure Python decides up to a few million t-sets of columns in minutes; whole families of thousands of columns need the compiled checkers named above.

## Runs on our own certificates (19 November 2026)

- github.com/arhancanli/quantavane-covering-arrays, three-row/: CPHF(3; 62, 5, 3), every one of its 37,820 column triples, and CPHF(3; 546, 16, 3) on its first 45 columns: A and B PASS (2 s and 6 s).
- Sherwood families of that repository on their first columns: CA_604_3_2876_4 (50 columns), CA_6611_3_2146_11 (40), CA_13117_3_10000_13 (40), CA_114205_4_145_13 (24, strength 4): A and B PASS, and each gives the row count its name states.
- Theorem K's family exported from qvrecords (q = 4, all 42 columns; q = 16, 40 columns): PASS; with one column replaced by a copy of another: FAIL (tests/test_checkers.py).
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
