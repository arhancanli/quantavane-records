#!/usr/bin/env python3
"""Exact checker for covering designs C(v, k, t) and C(v, k, t, m).

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
"""
import itertools
import math
import sys


def read_blocks(path, v, k):
    blocks = []
    for ln, line in enumerate(open(path), 1):
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        b = [int(x) for x in s]
        if len(b) != k or len(set(b)) != k or min(b) < 1 or max(b) > v:
            raise ValueError("line %d is not a %d-subset of 1..%d" % (ln, k, v))
        blocks.append(tuple(sorted(b)))
    return blocks


def check_a(blocks, v, k, t, m):
    if m == t:
        seen = set()
        for b in blocks:
            seen.update(itertools.combinations(b, t))
        return len(seen) == math.comb(v, t), "%d of %d t-subsets covered" % (len(seen), math.comb(v, t))
    # m > t: a point-by-point search for an m-subset S that meets every block in at most t - 1 points. Points are
    # taken in runs of equal type (the set of blocks holding them); points of one type are interchangeable, so a run
    # is entered in order (skipping one point skips the rest of its run). Bound: from each remaining run at most
    # min(points left in the run, the least spare capacity t - 1 - |S cap B| of its blocks).
    typ = {p: tuple(i for i, b in enumerate(blocks) if p in b) for p in range(1, v + 1)}
    order = sorted(range(1, v + 1), key=lambda p: (typ[p], -p))
    run_end = [0] * v
    for idx in range(v - 1, -1, -1):
        same = idx + 1 < v and typ[order[idx + 1]] == typ[order[idx]]
        run_end[idx] = run_end[idx + 1] if same else idx + 1
    spare = [t - 1] * len(blocks)
    S = []

    def bound(idx):
        tot, j = 0, idx
        while j < v:
            e = run_end[j]
            ty = typ[order[j]]
            tot += min([e - j] + [spare[i] for i in ty])
            j = e
        return tot

    def dfs(idx):
        if len(S) == m:
            return True
        if idx >= v or len(S) + bound(idx) < m:
            return False
        p = order[idx]
        ty = typ[p]
        if all(spare[i] > 0 for i in ty):
            for i in ty:
                spare[i] -= 1
            S.append(p)
            found = dfs(idx + 1)
            if found:
                return True
            S.pop()
            for i in ty:
                spare[i] += 1
        return dfs(run_end[idx])                  # leave p out, and with it the rest of its run
    sys.setrecursionlimit(10000)
    bad = dfs(0)
    return (not bad), ("counterexample m-subset %s" % sorted(S) if bad else "no m-subset misses every block")


def check_b(blocks, v, k, t, m):
    """written separately: points with the same blocks are interchangeable, so a bad m-subset is a choice of how many
    points x_T of each type T to take, with sum over the types inside each block at most t - 1 and sum x_T = m.
    Search over the types (largest first), x_T from 0 up, with the bound sum of min(count, capacity) over the rest."""
    masks = [frozenset(i for i, b in enumerate(blocks) if p in b) for p in range(1, v + 1)]
    count = {}
    for ty in masks:
        count[ty] = count.get(ty, 0) + 1
    free = count.pop(frozenset(), 0)          # points in no block: always allowed
    types = sorted(count, key=lambda ty: (-len(ty), sorted(ty)))
    caps = [t - 1] * len(blocks)
    nodes = 0
    sys.setrecursionlimit(10000)

    def most(ty):
        return min([count[ty]] + [caps[i] for i in ty])

    def go(j, size):
        nonlocal nodes
        nodes += 1
        if size >= m:
            return True
        if j == len(types) or size + sum(most(ty) for ty in types[j:]) < m:
            return False
        ty = types[j]
        for x in range(0, most(ty) + 1):
            for i in ty:
                caps[i] -= x
            ok = go(j + 1, size + x)
            for i in ty:
                caps[i] += x
            if ok:
                return True
        return False
    bad = go(0, free) if m > t or t > 1 else free >= m
    if m == t and not bad:
        return True, "no t-subset outside every block (%d nodes)" % nodes
    return (not bad), ("some m-subset meets every block in fewer than t points" if bad else
                       "no m-subset misses every block (%d nodes)" % nodes)


def main(argv):
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    try:
        path, v, k, t = args[0], int(args[1]), int(args[2]), int(args[3])
        m = int(args[4]) if len(args) > 4 else t
        assert 1 <= t <= m <= v and t <= k <= v
        blocks = read_blocks(path, v, k)
    except (IndexError, ValueError, AssertionError) as e:
        print("bad input: %s\n%s" % (e, __doc__))
        return 2
    want = None
    if "--blocks" in argv:
        want = int(argv[argv.index("--blocks") + 1])
    ok = True
    if want is not None and len(blocks) != want:
        print("FAIL  file has %d blocks, %d claimed" % (len(blocks), want))
        ok = False
    if len(set(blocks)) != len(blocks):
        print("note  the file repeats a block")
    for name, fn in (("A", check_a), ("B", check_b)):
        res, why = fn(blocks, v, k, t, m)
        print("%s %s  %s" % (name, "PASS" if res else "FAIL", why))
        ok &= res
    print("C(%d,%d,%d%s) <= %d: %s" % (v, k, t, "" if m == t else ",%d" % m, len(blocks), "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
