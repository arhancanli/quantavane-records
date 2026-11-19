# Fast matrix multiplication: living record table

**What a cell holds.** <n x m x p> <= r: a scheme computing the product of an n x m and an m x p matrix with r multiplications. Fewer is better.

**The table this follows.** The Collection of fast matrix multiplication algorithms, Universite de Lille (https://fmm.univ-lille.fr/; repository sedoglavic/fmm_digest). That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact programs from a public Quantavane repository, named there.

**Now.** 2 record rows covering 2 cells: 2 listed by the keeper.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## Cells

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| ⟨9×10×11⟩ | 640 | 644 (Perminov, 2026) | listed by the keeper, 2026-11-13 | [file](https://github.com/arhancanli/quantavane-fmm-records/tree/main/schemes/9x10x11) |
| ⟨9×11×11⟩ | 708 | 710 (Stoch, 2026) | listed by the keeper, 2026-11-13 | [file](https://github.com/arhancanli/quantavane-fmm-records/tree/main/schemes/9x11x11) |
