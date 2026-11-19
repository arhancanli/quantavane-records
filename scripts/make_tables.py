#!/usr/bin/env python3
"""Write tables/<field>/TABLE.md and SUBMIT.md from data/records.csv.

    python3 scripts/make_tables.py

CHECKER.md and checker.py in each folder are written by hand. MIT licence (LICENSE-CODE).
"""
import collections
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = "https://github.com/arhancanli/quantavane-records"

STATUS = {"keeper-accepted": "listed by the keeper", "submitted": "sent to the keeper",
          "verified-not-sent": "verified, not yet sent"}

FIELDS = {
    "covering-arrays": dict(
        title="Covering arrays",
        reads="CAN(t, k, v) <= N: a covering array CA(N; t, k, v), an N x k array over v symbols in which every t "
              "columns hold every t-tuple. Fewer rows is better.",
        keeper="Christopher Smolen's ca-tools table, the successor of Charles Colbourn's covering array tables "
               "(https://chrispbsmolen.github.io/ca-tools/). Rows marked \"past every table\" are widths or alphabets "
               "no published table reaches; there our array is smaller than every published construction we found.",
        cert="a covering array (N lines of k symbols 0..v-1), or a Sherwood covering perfect hash family in the "
             "ca-tools layout (n lines of k integers; give its kind, cphf or scphf, and the field if q is not prime)",
        summary=True),
    "covering-designs": dict(
        title="Covering designs",
        reads="C(v, k, t) <= b and C(v, k, t, m) <= b: b blocks of size k on v points such that every m-subset meets "
              "some block in at least t points (m = t for C(v, k, t)). Fewer blocks is better.",
        keeper="the Covering Repository (https://coveringrepository.com), which continues Daniel Gordon's La Jolla "
               "Covering Repository.",
        cert="the blocks, one per line, points 1..v"),
    "spherical-codes": dict(
        title="Spherical codes",
        reads="N points on the unit sphere in R^d with largest pairwise cosine C. A smaller largest cosine (a larger "
              "minimum angle) is better.",
        keeper="Henry Cohn's table of spherical codes (https://spherical-codes.org).",
        cert="the points, one per line, d comma-separated decimals, and the claimed largest cosine"),
    "constant-weight-codes": dict(
        title="Constant weight codes",
        reads="A(n, d, w) >= M: M binary words of length n and weight w, any two at distance at least d. More words "
              "is better.",
        keeper="Andries Brouwer's tables of constant weight codes (https://aeb.win.tue.nl/codes/Andw.html).",
        cert="the words, one per line (0/1 strings, positions of the ones, or Brouwer's $BASE= layout)"),
    "linear-codes": dict(
        title="Linear codes",
        reads="[n, k, d]_q: a linear code of length n and dimension k over GF(q) with minimum distance d. For given "
              "(q, n, k) a larger d is better.",
        keeper="Markus Grassl's codetables.de (https://www.codetables.de), MinT at the University of Salzburg "
               "(http://mint.sbg.ac.at/) and Eric Chen's database of quasi-twisted codes. The same code can hold the "
               "best value in two of these tables; each table's cell is a separate row here.",
        cert="a generator matrix, k lines of n field elements (Conway coding for prime powers)"),
    "cages": dict(
        title="Cages",
        reads="n(k, g) <= N: a k-regular graph of girth g on N vertices. Fewer vertices is better.",
        keeper="the House of Graphs cage table (https://houseofgraphs.org/meta-directory/cages), kept by Jan "
               "Goedgebeur, Jorik Jooken and Tibo Van den Eede. Since autumn 2026 it lists a first value for an empty "
               "cell only for graphs of at most 10,000 vertices.",
        cert="the graph in graph6 or sparse6 (gzip allowed) or as an edge list"),
    "extremal-graphs": dict(
        title="Extremal graphs without short cycles",
        reads="ex(n; t) >= m: a graph on n vertices with m edges and no cycle of length t or less. More edges is "
              "better.",
        keeper="the Combinatorics Wiki tables of extremal numbers ex(n; t) (http://combinatoricswiki.org/wiki/Ex(n;t)).",
        cert="the graph in graph6 or sparse6 or as an edge list"),
    "weighing-matrices": dict(
        title="Circulant weighing matrices",
        reads="CW(n, k): a circulant n x n matrix over {0, 1, -1} with W W^T = k I. The table records which (n, k) "
              "exist, and whether a proper one exists.",
        keeper="Daniel Gordon's La Jolla Circulant Weighing Matrix repository "
               "(https://github.com/dmgordo/circulant-weighing-matrices).",
        cert="the first row (n entries in {0, 1, -1})"),
    "weak-schur-numbers": dict(
        title="Weak Schur numbers",
        reads="WS(k) >= n: 1..n split into k sets, none holding x + y = z with x != y. A larger n is better.",
        keeper="Leaps in Bounds, kept by Daniel Munro (https://github.com/daniel-munro/leaps-in-bounds).",
        cert="the k sets, or the colours of 1..n"),
    "ramsey-numbers": dict(
        title="Ramsey numbers",
        reads="R(s, t) >= n (and more colours, and graph Ramsey numbers): a colouring of the edges of K_(n-1) with "
              "no forbidden monochromatic subgraph. A larger lower bound is better. Submissions are checked for "
              "cliques (R(s, t) and R(s_1, ..., s_r)); graph Ramsey numbers R(G, H) are listed but not yet open for "
              "submissions.",
        keeper="Stanislaw Radziszowski's dynamic survey Small Ramsey Numbers "
               "(https://www.combinatorics.org/ojs/index.php/eljc/article/view/DS1).",
        cert="the n x n colour matrix, one line per vertex"),
    "matrix-multiplication": dict(
        title="Fast matrix multiplication",
        reads="<n x m x p> <= r: a scheme computing the product of an n x m and an m x p matrix with r "
              "multiplications. Fewer is better.",
        keeper="the Collection of fast matrix multiplication algorithms, Universite de Lille "
               "(https://fmm.univ-lille.fr/; repository sedoglavic/fmm_digest).",
        cert="the scheme in the catalogue's Maple layout (Triad entries)", external=True),
    "degree-diameter": dict(
        title="Graphs and digraphs of given degree and diameter",
        reads="The largest known graphs or digraphs of maximum degree d and diameter k, in four tables: "
              "vertex-symmetric digraphs, Cayley graphs, bipartite graphs and arc-transitive graphs. A larger order "
              "is better. Submissions are open for vertex-symmetric digraphs given as Cayley or coset digraphs of "
              "PSL(2, q) or PGL(2, q) (times a cyclic group); the other three tables are listed but not yet open.",
        keeper="the Combinatorics Wiki degree-diameter tables "
               "(http://combinatoricswiki.org/wiki/The_Degree_Diameter_Problem_for_Vertex_Symmetric_Digraphs and its "
               "sister pages).",
        cert="a JSON spec with the group and the generators (the layout of github.com/arhancanli/quantavane-vsd)",
        external=True),
}


def esc(s):
    return s.replace("|", "\\|")


def table_rows(rs):
    out = ["| Cell | Ours | Previous best (source) | Status | Certificate |", "|---|---|---|---|---|"]
    for r in rs:
        prev = r["previous"]
        if r["previous_source"] and prev not in ("none listed",):
            prev = "%s (%s)" % (prev, r["previous_source"]) if prev else r["previous_source"]
        cert = r["certificate"]
        if cert.startswith("http"):
            cert = "[file](%s)" % cert
        st = STATUS[r["status"]] + (", " + r["status_date"] if r["status_date"] else "")
        if r["status_detail"].startswith("unlisted"):
            st += "; not listed under the table's size rule"
        if r["status_detail"].startswith("listed by Leaps in Bounds"):   # listed by another table than this one's keeper
            st = r["status_detail"]
        out.append("| %s | %s | %s | %s | %s |" % (esc(r["parameters"]), esc(r["ours"]), esc(prev), st, esc(cert)))
    return out


def ca_summary(rs):
    by = collections.defaultdict(lambda: [0, 0, 0, None, None])
    for r in rs:
        m = re.match(r"CAN\((\d+),k,(\d+)\)", r["parameters"])
        key = (int(m.group(1)), int(m.group(2)), r["table"].startswith("none"))
        b = by[key]
        b[0] += 1
        b[1] += int(r["cells"])
        b[2] += r["status"] == "keeper-accepted"
        lo, hi = int(r["k_from"]), int(r["k_to"])
        b[3] = lo if b[3] is None else min(b[3], lo)
        b[4] = hi if b[4] is None else max(b[4], hi)
    out = ["| t | v | where | record rows | widths covered | listed by the keeper | k from | k to |",
           "|---:|---:|---|---:|---:|---:|---:|---:|"]
    for (t, v, past), (n, cells, listed, lo, hi) in sorted(by.items()):
        out.append("| %d | %d | %s | %d | %s | %d | %s | %s |" % (
            t, v, "past every table" if past else "ca-tools", n, "{:,}".format(cells), listed,
            "{:,}".format(lo), "{:,}".format(hi)))
    return out


def main():
    with open(os.path.join(ROOT, "data", "records.csv"), newline="") as fh:
        recs = list(csv.DictReader(fh))
    by_field = collections.defaultdict(list)
    for r in recs:
        by_field[r["field"]].append(r)
    for f, cfg in FIELDS.items():
        rs = by_field.get(f, [])
        d = os.path.join(ROOT, "tables", f)
        os.makedirs(d, exist_ok=True)
        cnt = collections.Counter(r["status"] for r in rs)
        other = sum(1 for r in rs if r["status_detail"].startswith("listed by Leaps in Bounds"))
        cnt["keeper-accepted"] -= other                  # listed by another table, counted apart below
        cells = sum(int(r["cells"]) for r in rs)
        lines = ["# %s: living record table" % cfg["title"], "",
                 "**What a cell holds.** %s" % cfg["reads"], "",
                 "**The table this follows.** %s That table remains the reference for every cell; this one lists the "
                 "cells where Quantavane holds the best value we know, with the previous best and its source, and "
                 "takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md))." % (cfg["keeper"][0].upper() + cfg["keeper"][1:]), "",
                 "**Checker.** [CHECKER.md](CHECKER.md): " + (
                     "two separately written exact routines in [checker.py](checker.py)." if not cfg.get("external")
                     else "two separately written exact programs from a public Quantavane repository, named there."), "",
                 "**Now.** %s record rows covering %s cells: %s." % (
                     "{:,}".format(len(rs)), "{:,}".format(cells),
                     ", ".join(["{:,} {}".format(cnt[s], STATUS[s]) for s in STATUS if cnt[s]]
                               + (["{:,} listed by Leaps in Bounds".format(other)] if other else []))), "",
                 "Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best "
                 "is the value the table listed, or the best published value we found, when the record was set. "
                 "Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper "
                 "lists, the keeper's table holds the object.", ""]
        if cfg.get("summary"):
            lines += ["## By strength and alphabet", "", "Each record row is one array that holds a range of widths "
                      "k; data/records.csv has every row with its range, previous best and certificate.", ""]
            lines += ca_summary(rs)
            lines += ["", "## Rows past every table", ""]
            lines += table_rows([r for r in rs if r["table"].startswith("none")])
        else:
            lines += ["## Cells", ""]
            lines += table_rows(rs)
        open(os.path.join(d, "TABLE.md"), "w").write("\n".join(lines) + "\n")
        sub = ["# Submitting to the %s table" % cfg["title"].lower().replace("ramsey", "Ramsey").replace("schur", "Schur"), "",
               ("1. Check your construction yourself with [checker.py](checker.py) (see [CHECKER.md](CHECKER.md)). It "
                "is pure Python 3 with no other packages." if not cfg.get("external") else
                "1. Check your construction yourself with the programs named in [CHECKER.md](CHECKER.md)."),
               "2. Open an issue at %s/issues titled `%s: <cell> <value>`, and attach the certificate: %s. Say how "
               "the construction was found and how you would like to be credited." % (REPO, f, cfg["cert"]),
               "3. Within 72 hours we run both checkers on it and post their output in the issue.",
               "4. If both pass and your value beats the one in [TABLE.md](TABLE.md), we update the table, with your "
               "name and your source, in the next commit. If it beats the keeper's table too, send it to the keeper "
               "as well: the keeper's table is the reference, and we do not send other people's results on their "
               "behalf unless asked.",
               "5. If a check fails, the issue says which routine failed and why (for example the first uncovered "
               "tuple or the first short cycle), so the certificate can be fixed.", "",
               "Corrections to our own entries are welcome the same way. Certificates are published here under CC BY "
               "4.0, with credit to their authors.", ""]
        open(os.path.join(d, "SUBMIT.md"), "w").write("\n".join(sub))
    print("wrote TABLE.md and SUBMIT.md for %d fields" % len(FIELDS))


if __name__ == "__main__":
    main()
