#!/usr/bin/env python3
"""Exact checker for graphs of given girth: cages n(k, g) and extremal numbers ex(n; t).

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
"""
import gzip
import sys
from collections import deque


# ------------------------------------------------------------------------------------------- reading graphs
def _n_and_rest(data):                              # data: bytes minus 63; a byte 126 is 63 here
    if data[0] < 63:
        return data[0], data[1:]
    if data[1] < 63:
        return (data[1] << 12) | (data[2] << 6) | data[3], data[4:]
    n = 0
    for c in data[2:8]:
        n = (n << 6) | c
    return n, data[8:]


def parse_graph6(line):
    data = [c - 63 for c in line]
    n, rest = _n_and_rest(data)
    bits = [(c >> (5 - i)) & 1 for c in rest for i in range(6)]
    edges, pos = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[pos]:
                edges.append((i, j))
            pos += 1
    return n, edges


def parse_sparse6(line):
    data = [c - 63 for c in line[1:]]
    n, rest = _n_and_rest(data)
    k = max(1, (n - 1).bit_length())
    bits = [(c >> (5 - i)) & 1 for c in rest for i in range(6)]
    edges, v, pos = [], 0, 0
    while pos + 1 + k <= len(bits):
        b = bits[pos]
        x = 0
        for i in range(k):
            x = (x << 1) | bits[pos + 1 + i]
        pos += 1 + k
        if b:
            v += 1
        if x >= n or v >= n:
            break
        if x > v:
            v = x
        else:
            edges.append((x, v))
    return n, edges


def read_graph(path):
    raw = gzip.open(path, "rb").read() if path.endswith(".gz") else open(path, "rb").read()
    lines = [l.strip() for l in raw.splitlines() if l.strip()]
    first = lines[0]
    if first.startswith(b">>"):
        first = first.split(b"<<", 1)[1] if b"<<" in first else first
    if first.startswith(b":"):
        return parse_sparse6(first)
    parts = first.split()
    if len(parts) == 1 and all(63 <= c <= 126 for c in first):
        return parse_graph6(first)
    edges = [tuple(map(int, l.split()[:2])) for l in lines if l.split()[0].isdigit()]
    n = max(max(e) for e in edges) + 1
    return n, edges


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        if u == v:
            raise ValueError("loop at %d" % u)
        if v in adj[u]:
            raise ValueError("repeated edge %d-%d" % (u, v))
        adj[u].add(v)
        adj[v].add(u)
    return [sorted(a) for a in adj]


# ------------------------------------------------------------------------------------------- A: BFS girth
def girth_a(adj, cap):
    """the girth if it is at most cap, else cap + 1 (no cycle of length <= cap)"""
    best = cap + 1
    for s in range(len(adj)):
        dist, par = {s: 0}, {s: -1}
        dq = deque([s])
        while dq:
            u = dq.popleft()
            du = dist[u]
            if 2 * du + 1 >= best:
                break
            for w in adj[u]:
                if w == par[u]:
                    continue
                if w in dist:
                    best = min(best, du + dist[w] + 1)
                else:
                    dist[w], par[w] = du + 1, u
                    dq.append(w)
    return best


# ------------------------------------------------------------------------------------------- B: branch test
def _branch(adj, v, x, depth):
    """end vertices of non-backtracking walks of length 1..depth starting v -> x"""
    seen = {x}
    frontier = [(x, v)]
    for _ in range(depth - 1):
        nxt = []
        for u, prev in frontier:
            for w in adj[u]:
                if w != prev:
                    seen.add(w)
                    nxt.append((w, u))
        frontier = nxt
    return seen


def meets(adj, L):
    """a vertex v and neighbours x != y whose branches of depths floor(L/2), ceil(L/2) meet, if any: then there is a
    cycle of length <= L; if none anywhere, there is no cycle of length <= L"""
    a, b = L // 2, L - L // 2
    for v in range(len(adj)):
        nb = adj[v]
        if len(nb) < 2:
            continue
        short = {x: _branch(adj, v, x, a) for x in nb} if a > 0 else {x: set() for x in nb}
        long_ = {x: _branch(adj, v, x, b) for x in nb} if b != a else short
        for x in nb:
            for y in nb:
                if x != y and short[x] & long_[y]:
                    return v, x, y
    return None


def main(argv):
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    try:
        mode, path, p1, p2 = args[0], args[1], int(args[2]), int(args[3])
        n, edges = read_graph(path)
        adj = adjacency(n, edges)
    except (IndexError, ValueError, OSError) as e:
        print("bad input: %s\n%s" % (e, __doc__))
        return 2
    ok = True
    m = sum(len(a) for a in adj) // 2
    print("graph: %d vertices, %d edges, degrees %d..%d" % (n, m, min(map(len, adj)), max(map(len, adj))))
    if mode == "cage":
        k, g = p1, p2
        want = int(argv[argv.index("--vertices") + 1]) if "--vertices" in argv else None
        if any(len(a) != k for a in adj):
            print("FAIL  not %d-regular" % k)
            ok = False
        if want is not None and n != want:
            print("FAIL  %d vertices, %d claimed" % (n, want))
            ok = False
        ga = girth_a(adj, g)
        resa = ga == g
        print("A %s  girth %s" % ("PASS" if resa else "FAIL", ga if ga <= g else "> %d" % g))
        short = meets(adj, g - 1)
        exact = meets(adj, g) if short is None else None
        resb = short is None and exact is not None
        print("B %s  %s" % ("PASS" if resb else "FAIL", "no cycle of length <= %d; a cycle of length %d through vertex %d"
                            % (g - 1, g, exact[0]) if resb else ("a cycle of length <= %d at vertex %d" % (g - 1, short[0])
                                                                 if short else "no cycle of length %d found" % g)))
        ok &= resa and resb
        print("n(%d,%d) <= %d: %s" % (k, g, n, "PASS" if ok else "FAIL"))
    elif mode == "extremal":
        nn, t = p1, p2
        want = int(argv[argv.index("--edges") + 1]) if "--edges" in argv else None
        if n > nn:                                   # fewer is fine: isolated vertices need not be listed
            print("FAIL  %d vertices, at most %d allowed" % (n, nn))
            ok = False
        if want is not None and m != want:
            print("FAIL  %d edges, %d claimed" % (m, want))
            ok = False
        ga = girth_a(adj, t)
        resa = ga > t
        print("A %s  %s" % ("PASS" if resa else "FAIL", "no cycle of length <= %d" % t if resa else "girth %d" % ga))
        short = meets(adj, t)
        resb = short is None
        print("B %s  %s" % ("PASS" if resb else "FAIL", "no cycle of length <= %d" % t if resb else
                            "a cycle of length <= %d at vertex %d" % (t, short[0])))
        ok &= resa and resb
        print("ex(%d;%d) >= %d: %s" % (nn, t, m, "PASS" if ok else "FAIL"))
    else:
        print("mode must be cage or extremal")
        return 2
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
