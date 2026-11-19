# Checker: Ramsey numbers

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py FILE s t
python3 checker.py FILE --forbid 1:3,2:3,0:16
```

## What it checks

```text
Exact checker for Ramsey lower bounds: a colouring of the edges of K_n with no forbidden monochromatic clique.

    python3 checker.py FILE s t                 two colours: R(s, t) > n (no clique of size s in one colour and none
                                                of size t in the other, either way round)
    python3 checker.py FILE --forbid 1:3,2:3,0:16    any number of colours: colour symbol c has no clique of size s_c

FILE: n lines of n symbols (the n x n colour matrix; the diagonal is ignored but must be the same symbol 0); the
layout of github.com/arhancanli/quantavane-ramsey-records. The matrix must be symmetric. When it is circulant
(entry (i, j) depends only on j - i mod n), translations are automorphisms of every colour class, so it is enough to
search cliques through vertex 0; the program checks circulance before using it. Two separately written routines:
  A  branch and bound on bitsets with greedy-colouring bounds (in the style of Tomita's MCQ): is there a clique of
     size s in this colour?
  B  Bron-Kerbosch with pivoting on Python sets, stopped at the first clique of size s; for s = 3 it checks every
     edge's two neighbourhoods for a common vertex instead.
The repository above also has check.c, a third, compiled checker for two colours. Exit 0 = both PASS.
Pure Python: quick for triangles and small cliques; large cliques in dense colours (say K_15 on 280 vertices) can
take hours, where check.c takes minutes. MIT licence (LICENSE-CODE).
```

## Other checkers

`check.c` in github.com/arhancanli/quantavane-ramsey-records is a third, compiled checker for two colours (branch and bound with colouring bounds, from vertex 0 for circulant colourings).

## Limits

Quick for triangles and small cliques; large cliques in dense colours can take hours in Python, where check.c takes minutes.

## Runs on our own certificates (19 November 2026)

- R(3,18) >= 101 (100 points, not circulant): A and B PASS, 124 s.
- R(5,15) >= 283 (a circulant colouring of 282 points): A and B PASS, 95 s.
- R(4,18) >= 211 (210 points, not circulant) and R(8,12) >= 523 (522 points) were stopped after several minutes, to keep the machine's load down; for colourings this size check.c in quantavane-ramsey-records is the practical checker.
- Paley colourings of K_5 and K_17 and the three-colouring of K_16 from GF(16) pass; K_17 with one edge recoloured fails.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
