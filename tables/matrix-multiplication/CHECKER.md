# Checker: fast matrix multiplication

A scheme of rank r for <n x m x p> is a list of r triads (alpha, beta, gamma) such that the sum of alpha (x) beta (x) gamma is the matrix multiplication tensor. The checkers are the three exact programs of github.com/arhancanli/quantavane-fmm-records, folder `check/`, each written separately:

- `fmm_check.py`: the tensor sum compared with the matrix multiplication tensor entry by entry, in numpy int64 slices for integer schemes and Python fractions for rational ones;
- `fmm_check2.c` (C): the full expansion after clearing denominators, modulo primes near 2^61, with enough primes that a height bound turns agreement modulo every prime into equality over the rationals, plus a random evaluation;
- `fmm_check3.py`: every coefficient sum_i P_i[a][b] Q_i[c][d] R_i[f][e] = [b = c][a = e][d = f], with entries scaled to integers by the least common denominator.

`sh check/run_all.sh` in that repository runs all three on every scheme there and on two damaged copies of each (one coefficient changed; the last term dropped), which must fail.

```sh
cc -O2 -o fmm_check2 check/fmm_check2.c -lm
./fmm_check2 SCHEME.mpl
python3 check/fmm_check.py SCHEME.mpl --expect-rank r --expect-format NxMxP
```

## Runs on our own certificates (19 November 2026)

- fmm_check2.c, compiled and run here on the two record schemes: <9 x 10 x 11> rank 640 and <9 x 11 x 11> rank 708, PASS (full expansion of 23,849,607 and 27,239,328 products agrees modulo three primes near 2^61).
- The repository's run_all.sh (14 November 2026, check/results.txt there) passed every scheme with all three programs and rejected every damaged copy.
