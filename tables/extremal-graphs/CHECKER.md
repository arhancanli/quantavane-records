# Checker: extremal graphs without short cycles

`checker.py` in this folder, pure Python 3 (standard library only, MIT licence). It prints one line per routine and a verdict, and exits 0 only when both routines pass.

```sh
python3 checker.py extremal FILE n t [--edges m]
```

## What it checks

```text
Exact checker for graphs of given girth: cages n(k, g) and extremal numbers ex(n; t).

    python3 checker.py cage     FILE k g [--vertices N]    a k-regular graph of girth exactly g (n(k, g) <= N)
    python3 checker.py extremal FILE n t [--edges m]       a graph on n vertices with m edges and no cycle of
                                                            length t or less (ex(n; t) >= m)

FILE: graph6 or sparse6 (one graph; .gz allowed) or an edge list ("u v" per line, vertices 0..n-1). The program
checks the graph is simple (no loops, no repeated edges) and runs two separately written routines:

  A  breadth-first search from every vertex, stopped at the depth that can still close a cycle shorter than the
     best found; gives the girth exactly (or proves it exceeds the bound that matters).
  B  the branch test: for every vertex v and every two neighbours x != y, the non-backtracking walks of length
     1..a that start v x and those of length 1..b that start v y end at disjoint sets of vertices, with
     a = floor((L-1)/2), b = ceil((L-1)/2) for L the shortest forbidden length; two such walks meeting give a closed
     walk, hence a cycle, of length <= L - 1, and every cycle that short makes two of them meet. For a cage it also
     finds a cycle of length g the same way (walks of length floor(g/2) and ceil(g/2)), so the girth is exactly g.
Exit 0 = both PASS. Pure Python: fine to tens of thousands of edges; larger graphs take hours.
MIT licence (LICENSE-CODE).
```

## Limits

Seconds for the sizes in this table.

## Runs on our own certificates (19 November 2026)

- All 59 extremal graphs the lab holds for this table (the 24 records here among them): A and B PASS; each of the 59 damaged copies the lab keeps (one edge added) is rejected by both routines.
- The tests in `tests/test_checkers.py`: a small known object passes and a damaged copy fails.
