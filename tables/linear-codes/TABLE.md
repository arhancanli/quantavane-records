# Linear codes: living record table

**What a cell holds.** [n, k, d]_q: a linear code of length n and dimension k over GF(q) with minimum distance d. For given (q, n, k) a larger d is better.

**The table this follows.** Markus Grassl's codetables.de (https://www.codetables.de), MinT at the University of Salzburg (http://mint.sbg.ac.at/) and Eric Chen's database of quasi-twisted codes. The same code can hold the best value in two of these tables; each table's cell is a separate row here. That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact routines in [checker.py](checker.py).

**Now.** 39 record rows covering 39 cells: 8 listed by the keeper, 31 sent to the keeper.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## Cells

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| [220,13,125] over GF(3) | 125 | 124 (the lower bound listed at codetables.de) | sent to the keeper, 2026-10-24 | sent to the keeper; not yet public |
| [153,10,99] over GF(4) | 99 | 98 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [182,13,114] over GF(4) | 114 | 112 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [195,13,122] over GF(4) | 122 | 121 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [39,12,19] over GF(5) | 19 | 18 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [63,14,33] over GF(5) | 33 | 32 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [70,12,40] over GF(5) | 40 | 39 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [93,13,55] over GF(5) | 55 | 54 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [57,10,36] over GF(7) | 36 | 35 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [66,10,43] over GF(7) | 43 | 42 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [66,12,40] over GF(7) | 40 | 39 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [87,8,63] over GF(7) | 63 | 62 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [90,11,59] over GF(7) | 59 | 58 (the lower bound listed at codetables.de) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| [102,10,73] over GF(9) | 73 | 72 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [130,10,96] over GF(9) | 96 | 95 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [130,13,90] over GF(9) | 90 | 89 (the lower bound listed at codetables.de) | listed by the keeper, 2026-11-10 | with the keeper's table; no public Quantavane file yet |
| [168,13,93]_3 | 93 | 89 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [176,14,96]_3 | 96 | 94 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [219,13,125]_3 | 125 | 122 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [220,13,125]_3 | 125 | 123 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [140,14,82]_4 | 82 | 80 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [153,10,99]_4 | 99 | 98 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [182,13,114]_4 | 114 | 110 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [195,13,122]_4 | 122 | 116 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [33,10,17]_5 | 17 | 15 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [39,12,19]_5 | 19 | 18 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [63,14,33]_5 | 33 | 32 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [70,12,40]_5 | 40 | 37 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [93,13,55]_5 | 55 | 52 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [57,10,36]_7 | 36 | 34 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [66,10,43]_7 | 43 | 39 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [66,12,40]_7 | 40 | 37 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [87,8,63]_7 | 63 | 62 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [90,11,59]_7 | 59 | 53 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [92,12,61]_8 | 61 | 56 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [102,10,73]_9 | 73 | 72 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [104,10,75]_9 | 75 | 74 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [130,10,96]_9 | 96 | 95 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
| [130,13,90]_9 | 90 | 81 (MinT's best known value) | sent to the keeper, 2026-11-10 | sent to the keeper; not yet public |
