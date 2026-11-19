# Circulant weighing matrices: living record table

**What a cell holds.** CW(n, k): a circulant n x n matrix over {0, 1, -1} with W W^T = k I. The table records which (n, k) exist, and whether a proper one exists.

**The table this follows.** Daniel Gordon's La Jolla Circulant Weighing Matrix repository (https://github.com/dmgordo/circulant-weighing-matrices). That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact routines in [checker.py](checker.py).

**Now.** 2 record rows covering 2 cells: 2 sent to the keeper.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## Cells

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| CW(146, 8^2) | exists | open in the table (circulant weighing matrix tables (Gordon)) | sent to the keeper, 2026-11-08 | sent to the keeper; not yet public |
| CW(242, 9^2) | exists | open in the table (circulant weighing matrix tables (Gordon)) | sent to the keeper, 2026-11-08 | sent to the keeper; not yet public |
