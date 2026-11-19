# Ramsey numbers: living record table

**What a cell holds.** R(s, t) >= n (and more colours, and graph Ramsey numbers): a colouring of the edges of K_(n-1) with no forbidden monochromatic subgraph. A larger lower bound is better. Submissions are checked for cliques (R(s, t) and R(s_1, ..., s_r)); graph Ramsey numbers R(G, H) are listed but not yet open for submissions.

**The table this follows.** Stanislaw Radziszowski's dynamic survey Small Ramsey Numbers (https://www.combinatorics.org/ojs/index.php/eljc/article/view/DS1). That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact routines in [checker.py](checker.py).

**Now.** 8 record rows covering 8 cells: 4 sent to the keeper, 2 verified, not yet sent, 2 listed by Leaps in Bounds.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## Cells

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| R(3,3,16) | 293 | 291 (Exoo and Tatarevic, 2015) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| R_8(3) | 5384 | 5288 (revision 18 of the survey, 2026) | sent to the keeper, 2026-11-09 | sent to the keeper; not yet public |
| R(B5,B9) | 29 | 28 (Van Overberghe, 2024) | verified, not yet sent | not yet public |
| R(P,P) | 22 | 21 (Burr's bound, 1981, as given by Si, 2026) | verified, not yet sent | not yet public |
| R(3,18) ≥ 101 | 101 | 100 (a triangle-free graph on 99 points in Google's ramsey_number_bounds repository, 2026, the value revision 18 (24 April 2026) of Small Ramsey Numbers lists; the best circulant graph has 97 points (Coniglio et al., 2026)) | sent to the keeper, 2026-11-09 | [file](https://github.com/arhancanli/quantavane-ramsey-records/blob/main/colourings/R3_18_K100.txt) |
| R(4,18) ≥ 211 | 211 | 209 (209, from a colouring of 208 points in Google's ramsey_number_bounds repository, 2026 (revision 18 (24 April 2026) of Small Ramsey Numbers lists 205); our own first step past it was 210, earlier the same day) | sent to the keeper, 2026-11-09 | [file](https://github.com/arhancanli/quantavane-ramsey-records/blob/main/colourings/R4_18_K210.txt) |
| R(5,15) ≥ 283 | 283 | 275 (275, from a general inequality, as revision 18 (24 April 2026) of Stanisław Radziszowski's survey Small Ramsey Numbers prints it; our own circulant colouring of 278 points (R(5,15) ≥ 279, 24 October 2026) was the first step past it, and Radziszowski said it will go into the next revision of his survey) | listed by Leaps in Bounds 2026-11-09 (verified there); sent to Radziszowski's survey Small Ramsey Numbers, whose current version (24 April 2026) predates it | [file](https://github.com/arhancanli/quantavane-ramsey-records/blob/main/colourings/R5_15_K282.txt) |
| R(8,12) ≥ 523 | 523 | 511 (511, from a colouring of 510 points in Google's ramsey_number_bounds repository, 2026; our own steps past it were 512 (25 October 2026) and 519 (3 November 2026, sent to Stanisław Radziszowski for his survey that day)) | listed by Leaps in Bounds 2026-11-09 (verified there); sent to Radziszowski's survey Small Ramsey Numbers, whose current version (24 April 2026) predates it | [file](https://github.com/arhancanli/quantavane-ramsey-records/blob/main/colourings/R8_12_K522.txt) |
