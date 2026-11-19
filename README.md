# Quantavane records

Every record Quantavane has released, one row each, with its parameters, our value, the previous best and its source,
its status with the table's keeper, where its certificate is and how it was checked. A small Python package,
`qvrecords`, computes and checks any cell of a family theorem on demand, so the cells a theorem settles are queryable
without being stored. And living record tables, each with an exact checker, take better constructions from anyone.

Released means listed by a table's keeper, sent to a keeper, or published on quantavane.com or in a public
Quantavane repository. Results that are held, or that exist only in the lab's own files, are not here. Neither are
results whose second check is not on file.

## What is here

| Path | What it holds |
|---|---|
| `data/records.csv` | 3,161 records, one row each (columns below) |
| `data/families.csv` | the family theorem rows: 21 ranges of Theorem K, counted apart from the records |
| `qvrecords/` | the package: `cell`, `check`, `cell_at`, `family_cells`, `records`, `record`, `ca_record` |
| `tables/` | 12 living record tables: `TABLE.md`, `CHECKER.md`, `checker.py`, `SUBMIT.md` in each |
| `tests/` | the 1,000-cell family test and the checker tests |
| `scripts/make_tables.py` | writes each `TABLE.md` and `SUBMIT.md` from `data/records.csv` |
| `range/` | the range analysis behind the family cells (`RANGE.md`), its count programs with their outputs, and the two array checkers with their logs |

## The numbers (19 November 2026)

**Records** (`data/records.csv`): 3,161 rows covering 83,038,040 cells: 1,536,296 cells of keepers' tables and
81,501,744 covering-array cells past every published table.

| status | rows | cells |
|---|---:|---:|
| keeper-accepted: listed in the keeper's table | 2,655 | 1,533,823 |
| submitted: sent to the keeper, not yet listed | 330 | 2,323 |
| verified-not-sent: checked and published by Quantavane, not (yet) with a keeper | 176 | 81,501,894 |

By kind, in cells: 188,906 beat the best value listed before them; 1,248,484 are first values for cells the
table had left empty; 98,906 are project-scheduling instances (each counted once, for its best result); and
81,501,744 are covering-array cells past every published table, where our array is smaller than every published
construction we found (`kind` = `past-tables`; github.com/arhancanli/quantavane-covering-arrays, `wide/README.md`,
"Which cells count"). The records are in 21 fields; `python3 -m qvrecords summary` prints them by field and status.

**Family cells** (`data/families.csv`, served by `qvrecords`): 227,658,593,685,328,988,782,135,588 cells of the form
CAN(3, k, q) <= 3q^3 - 2, from Theorem K over q = 16 and q = 4^s for 3 <= s <= 22. The 168 cells at q = 16
(379 <= k <= 546) are listed in the ca-tools table (the keeper merged them on 17 November 2026). The others lie past
every published table: no table lists strength 3 at q >= 64. In words: one theorem settles 2.3e26 cells past every
published construction. That is one result, not 2.3e26 separate discoveries; the count is set by where the declared
range stops (s <= 16, the range until 20 November 2026, gives 1.4e19; s = 17 to 22 were added on 20 November 2026,
`range/RANGE.md` sections 8 and 9).

## Counting rules

1. A record is one cell of one table where our own object (a design, code, graph, colouring, array, scheme or
   schedule) gives a better value than the table listed when we found it. A first value for an empty cell is counted
   apart from a record that beats a listed value (`kind`).
2. One row can hold many cells. A covering array with k columns gives every width up to k (deleting columns keeps a
   covering array), so each covering-array row covers the range `k_from`..`k_to`, and `cells` says how many. A
   scheduling row covers one dataset of the Solutions Update, and `cells` is the number of its instances we hold.
3. Records and family cells are kept apart. A family cell is counted only inside the range where the lab's comparison
   with every published construction class shows the theorem's value is smaller (see "Family cells" below). Family
   cells are never added to the record counts, and record counts never include family cells.
4. Not counted, and not here: ties with a listed value; cells a keeper derives from our arrays by fusing symbols or
   other rules; cells where published work is smaller; superseded values; results held back; results without a
   second check on file.
5. The same object can hold the best value in two tables (for example a linear code in codetables.de and in MinT).
   Each table's cell is its own row.

## The columns of records.csv

`id` (the lab's ledger id), `part` (S: one object, one cell; R: one object, a range of cells), `field`, `table` (the
keeper's table, or "none: past every table"), `table_url`, `parameters`, `cells`, `k_from`, `k_to`, `ours`,
`previous` and `previous_source` (the value the table listed, or the best published value we found, and where it is
from), `kind`, `found` (when the object was found, where the ledger has it), `status`, `status_date`,
`status_detail` (the ledger's own words), `certificate`, `checkers`.

`certificate` is a URL when a public file holds the object (2,725 rows: our public repositories, or the keeper's
page); otherwise it says "with the keeper's table; no public Quantavane file yet" (71), "sent to the keeper; not yet
public" (222) or "not yet public" (143). `checkers` names the checks on file: two or more separately written
programs, or a keeper's own validation on listing, with red controls where the lab ran them.

## How to query

Python 3.8 or later, standard library only.

```sh
python3 -m qvrecords summary                       # counts by status and field, family cells by q
python3 -m qvrecords record S-ramsey-r-5-15         # one record
python3 -m qvrecords ca 3 11 1820                   # the record row holding CAN(3, 1820, 11)
python3 -m qvrecords cell theorem-k 64 6000 --check # one family cell, its certificate and the independent check
python3 -m qvrecords index 1000000000000 --check    # the 10^12-th counted family cell
```

```python
from qvrecords import cell, check, cell_at, family_cells
value, cert = cell("theorem-k", {"q": 4**10, "k": 2 * 10**12})
ok, report = check(cert)                 # (True, [("PASS", "modulus irreducible ..."), ...])
total, by_q = family_cells()             # 227658593685328988782135588
q, k = cell_at(123456789)                # the cell with that index
```

`data/records.csv` reads in any CSV tool.

## Family cells

Theorem K: for q = 4^s, let c be a non-cube of GF(q) and b a root of X^3 - c. With the points x of PG(2, q) as
elements of GF(q^3), the 3 x 2(q^2 + q + 1) array of columns B_x = (x, x^-1, x) and E_x = (x, b x^-1, b x) is a
covering perfect hash family, so by Sherwood, Martirosyan and Colbourn (J. Combin. Des. 14 (2006) 202-213)
CAN(3, k, q) <= 3q^3 - 2 for every k <= 2(q^2 + q + 1). The construction and a sketch of the proof are public in
github.com/arhancanli/quantavane-covering-arrays, folder `three-row/`; a note with the full proof is in preparation.

`cell("theorem-k", {"q": q, "k": k})` builds the certificate of one cell: the field (its modulus and non-cube), the
rule that gives every column, sample columns, the counted range and its basis. `theorem_k.ca_entry(q, k, row, col)`
gives any entry of the covering array itself. `check(cert)` is written separately and imports nothing from the
construction: it tests the modulus by Rabin's test, the non-cube, the value and the counted range (recomputed), the
sample columns (recomputed with other field arithmetic), and the covering property on random triples among the first
k columns and on the hard cases of the proof's Lemma 2; with `exhaustive=True` it decides every triple of the family,
which is practical for q <= 16. The sampled check tests the construction at the cell; the proof covers every triple.

The counted range at q = 4^s >= 64 is k > R(q), where R(q) is the largest width that any published construction class
with at most 3q^3 - 2 rows reaches, by the lab's range analysis (`qvrecords.theorem_k.counted_range` gives the
formula). That analysis, `range/RANGE.md`, was refereed inside the lab on 18 November 2026 and has not been peer
reviewed; published constructions we did not find could reach further. Besides its proofs it rests on three stated
facts about the literature, and the most fragile is P2: nobody has published a 3-row search at q >= 64, and such a
search could remove the bottom of a q's range. P2 is the weakest input of the count, and we say so plainly: it is a
statement about what has been published, not a theorem. The number is one theorem settling many cells past every
published construction, not that many separate discoveries. At q = 16 the counted range is the one the keeper lists.

Tests: `python3 -m unittest discover -s tests` draws 1,000 cells uniformly from all 2.3e26 counted cells and checks
each with both routines (the construction's and the independent one must agree on every column they both compute,
both find every sampled triple covered, and the certificate must pass); it also tests 20 cells at each q, the ends of
every range, eight damaged certificates (all fail) and the whole covering array at q = 4 by its definition.
`QV_SLOW=1` adds the exhaustive check at q = 16 (296,480 triples, about 6 seconds).

## Living record tables and submissions

| Table | Record rows here | Checker |
|---|---:|---|
| [covering arrays](tables/covering-arrays/TABLE.md) | 1,807 | covering arrays and Sherwood families |
| [covering designs](tables/covering-designs/TABLE.md) | 132 | C(v, k, t) and C(v, k, t, m) |
| [spherical codes](tables/spherical-codes/TABLE.md) | 443 | exact largest cosine |
| [constant weight codes](tables/constant-weight-codes/TABLE.md) | 39 | distance and packing |
| [linear codes](tables/linear-codes/TABLE.md) | 39 | minimum distance over GF(q) |
| [Ramsey numbers](tables/ramsey-numbers/TABLE.md) | 8 | cliques in each colour (graph Ramsey cells listed, not yet open) |
| [cages](tables/cages/TABLE.md) | 16 | regularity and girth |
| [extremal graphs](tables/extremal-graphs/TABLE.md) | 24 | edges and no short cycle |
| [degree and diameter](tables/degree-diameter/TABLE.md) | 57 | vertex-symmetric digraphs only, for now |
| [weak Schur numbers](tables/weak-schur-numbers/TABLE.md) | 7 | weakly sum-free partitions |
| [circulant weighing matrices](tables/weighing-matrices/TABLE.md) | 2 | W W^T = k I and properness |
| [fast matrix multiplication](tables/matrix-multiplication/TABLE.md) | 2 | the three exact checkers of quantavane-fmm-records |

Each table follows a keeper's table, which stays the reference; ours lists the cells where we hold the best value we
know, with the previous best and its source. To submit a better construction, open an issue here with the
certificate in the format the table's `CHECKER.md` names. We run both checkers within 72 hours and post their output
in the issue; a construction that passes and beats the listed value goes into the table with its author's name.

Fields with records here but no living table yet, and why: prime gaps (the endpoints are probable primes, so the
check is not exact in our sense; the prime gap list checks every gap), quantum codes (an exact distance check for
these lengths needs a compiled program we have not published), digital nets and ordered orthogonal arrays, permutation
codes and covering codes (same reason), project scheduling and vehicle routing (the checks need the benchmark
instance files; the Solutions Update and CVRPLIB check every submission themselves), atomic clusters (energies are
floating-point, not exact), and public benchmark leaderboards (the keepers score them).

## Licence

Data, tables and this README: CC BY 4.0 (`LICENSE`). Code (`qvrecords/`, `scripts/`, `tests/`, every `checker.py`):
MIT (`LICENSE-CODE`). The previous-best values are quoted from the tables and papers named in each row; please credit
those tables and their keepers too.

## How to cite

Cite this repository for the records, the family cells and the checkers, and for a single record also the table that
lists it.

```bibtex
@misc{quantavane-records,
  author       = {Canli, Arhan},
  title        = {Quantavane records: released records with certificates, family cells on demand, and living record tables},
  year         = {2026},
  howpublished = {\url{https://github.com/arhancanli/quantavane-records}},
  note         = {Quantavane; data CC BY 4.0, code MIT}
}
```

[CITATION.cff](CITATION.cff) gives the same reference for GitHub's "Cite this repository" button.

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane), who takes responsibility for the dataset.

The searches, constructions, checks, programs and text here were produced by AI agents (Claude, by Anthropic) working
for Quantavane under Arhan Canli's direction. The records marked keeper-accepted are listed in their keepers' tables;
the others are not listed by a keeper yet, and the range analysis behind the family cells has not been reviewed
outside Quantavane.
