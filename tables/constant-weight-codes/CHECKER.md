# Checker: constant weight codes

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE n d w [--words M]
```

## What it checks

```text
Exact checker for constant weight codes: A(n, d, w) >= M.

A binary constant weight code is a set of M words of length n, each with exactly w ones, any two at Hamming distance
at least d. The table (Andries Brouwer's tables of constant weight codes) keeps the largest known M.

    python3 checker.py FILE n d w [--words M]

FILE: one word per line, as n characters 0/1, or as the positions of its ones (1..n, separated by spaces or commas),
or in the layout of Brouwer's files: a first line $BASE=16 (or 2, 8, 10, 32) and then each word as one number in that
base whose binary digits are the word. The program checks every word has length at most n and weight w, the words are
distinct, and runs two separately written routines:
  A  every pair of words: the population count of their XOR is at least d.
  B  the packing form: two words of weight w are at distance 2(w - |a & b|), so the code has distance >= d exactly
     when no set of s = w - ceil(d/2) + 1 positions lies in two words; every s-subset of every word goes into one
     dictionary, which must never see a subset twice.
Exit 0 = both PASS. A is quadratic in M (about a minute at M = 12,000). MIT licence (LICENSE-CODE).
```

## Limits

Routine A is quadratic in M: about 20 seconds at M = 12,285, 4 minutes at M = 40,740 and 15 minutes at M = 80,896.

## Runs on our own certificates (19 November 2026)

- 28 of the 39 codes here, from Andries Brouwer's own files of them: A and B PASS, from under a second to 15 minutes (A(58,4,5) >= 80,896, the largest). The other 11 were not run: no local copy of the file.
- A(27,4,6) >= 12,285 with one word replaced by a copy of another: FAIL in both routines.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
