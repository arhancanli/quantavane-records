# Weak Schur numbers: living record table

**What a cell holds.** WS(k) >= n: 1..n split into k sets, none holding x + y = z with x != y. A larger n is better.

**The table this follows.** Leaps in Bounds, kept by Daniel Munro (https://github.com/daniel-munro/leaps-in-bounds). That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact routines in [checker.py](checker.py).

**Now.** 7 record rows covering 7 cells: 3 listed by the keeper, 4 verified, not yet sent.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## Cells

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| WS(12) | 818354 | 815314 (Ageron and coauthors, 2021) | verified, not yet sent | not yet public |
| WS(13) | 2659652 | 2649772 (Ageron and coauthors, 2021) | verified, not yet sent | not yet public |
| WS(14) | 8592720 | 8560800 (Ageron and coauthors, 2021) | verified, not yet sent | not yet public |
| WS(15) | 27106320 | 27074400 (Ageron and coauthors, 2021) | verified, not yet sent | not yet public |
| WS(7) | 2154 | 2146 (Rowley, 2020) | listed by the keeper, 2026-11-09 | with the keeper's table; no public Quantavane file yet |
| WS(8) | 7002 | 6976 (Rowley, 2020) | listed by the keeper, 2026-11-09 | with the keeper's table; no public Quantavane file yet |
| WS(9) | 22620 | 22536 (Ageron and coauthors, 2021) | listed by the keeper, 2026-11-09 | with the keeper's table; no public Quantavane file yet |
