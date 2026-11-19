# Covering arrays: living record table

**What a cell holds.** CAN(t, k, v) <= N: a covering array CA(N; t, k, v), an N x k array over v symbols in which every t columns hold every t-tuple. Fewer rows is better.

**The table this follows.** Christopher Smolen's ca-tools table, the successor of Charles Colbourn's covering array tables (https://chrispbsmolen.github.io/ca-tools/). Rows marked "past every table" are widths or alphabets no published table reaches; there our array is smaller than every published construction we found. That table remains the reference for every cell; this one lists the cells where Quantavane holds the best value we know, with the previous best and its source, and takes better constructions from anyone (see [SUBMIT.md](SUBMIT.md)).

**Checker.** [CHECKER.md](CHECKER.md): two separately written exact routines in [checker.py](checker.py).

**Now.** 1,807 record rows covering 82,937,799 cells: 1,774 listed by the keeper, 33 verified, not yet sent.

Built from [data/records.csv](../../data/records.csv) by scripts/make_tables.py; the previous best is the value the table listed, or the best published value we found, when the record was set. Where no public Quantavane file holds the certificate yet, the column says so; for cells a keeper lists, the keeper's table holds the object.

## By strength and alphabet

Each record row is one array that holds a range of widths k; data/records.csv has every row with its range, previous best and certificate.

| t | v | where | record rows | widths covered | listed by the keeper | k from | k to |
|---:|---:|---|---:|---:|---:|---:|---:|
| 3 | 4 | ca-tools | 10 | 90,398 | 10 | 2,872 | 100,000 |
| 3 | 4 | past every table | 1 | 65,004 | 0 | 100,001 | 165,004 |
| 3 | 5 | ca-tools | 24 | 90,920 | 23 | 56 | 100,000 |
| 3 | 5 | past every table | 1 | 54,375 | 0 | 100,001 | 154,375 |
| 3 | 7 | ca-tools | 22 | 91,339 | 22 | 1,474 | 100,000 |
| 3 | 7 | past every table | 1 | 37,618 | 0 | 100,001 | 137,618 |
| 3 | 8 | ca-tools | 22 | 91,003 | 22 | 291 | 100,000 |
| 3 | 8 | past every table | 1 | 16,455 | 0 | 100,001 | 116,455 |
| 3 | 9 | ca-tools | 23 | 92,195 | 23 | 365 | 100,000 |
| 3 | 9 | past every table | 1 | 121,150 | 0 | 100,001 | 221,150 |
| 3 | 11 | ca-tools | 17 | 91,768 | 17 | 551 | 100,000 |
| 3 | 11 | past every table | 1 | 107,658 | 0 | 100,001 | 207,658 |
| 3 | 13 | ca-tools | 16 | 91,379 | 16 | 774 | 100,000 |
| 3 | 13 | past every table | 1 | 34,368 | 0 | 100,001 | 134,368 |
| 3 | 16 | ca-tools | 18 | 91,376 | 18 | 314 | 100,000 |
| 3 | 16 | past every table | 1 | 222,146 | 0 | 100,001 | 322,146 |
| 3 | 17 | ca-tools | 14 | 91,575 | 14 | 1,329 | 100,000 |
| 3 | 17 | past every table | 1 | 6,533 | 0 | 100,001 | 106,533 |
| 3 | 19 | ca-tools | 15 | 92,181 | 15 | 402 | 100,000 |
| 3 | 19 | past every table | 1 | 60,839 | 0 | 100,001 | 160,839 |
| 3 | 23 | ca-tools | 7 | 86,823 | 7 | 2,493 | 100,000 |
| 3 | 23 | past every table | 1 | 202,074 | 0 | 100,001 | 302,074 |
| 3 | 25 | ca-tools | 7 | 82,893 | 7 | 2,977 | 100,000 |
| 3 | 25 | past every table | 1 | 308,135 | 0 | 100,001 | 408,135 |
| 3 | 27 | past every table | 2 | 436,669 | 0 | 21,197 | 534,635 |
| 3 | 29 | past every table | 2 | 609,399 | 0 | 26,131 | 684,678 |
| 3 | 31 | past every table | 2 | 811,081 | 0 | 31,777 | 860,049 |
| 3 | 32 | past every table | 2 | 910,718 | 0 | 33,825 | 944,542 |
| 3 | 37 | past every table | 2 | 1,486,309 | 0 | 265,264 | 1,751,572 |
| 3 | 41 | past every table | 1 | 260,627 | 0 | 383,199 | 643,825 |
| 3 | 43 | past every table | 1 | 2,453,971 | 0 | 454,386 | 2,908,356 |
| 3 | 61 | past every table | 1 | 8,177,605 | 0 | 1,580,390 | 9,757,994 |
| 3 | 67 | past every table | 1 | 7,404,363 | 0 | 2,205,466 | 9,609,828 |
| 3 | 73 | past every table | 2 | 16,252,780 | 0 | 2,989,885 | 19,242,664 |
| 3 | 97 | past every table | 2 | 19,816,766 | 0 | 8,178,895 | 27,995,660 |
| 3 | 109 | past every table | 2 | 21,645,101 | 0 | 12,349,850 | 33,994,950 |
| 4 | 5 | ca-tools | 10 | 6,450 | 10 | 247 | 7,853 |
| 4 | 7 | ca-tools | 27 | 5,051 | 27 | 641 | 9,904 |
| 4 | 8 | ca-tools | 32 | 9,245 | 32 | 68 | 10,000 |
| 4 | 9 | ca-tools | 46 | 9,052 | 46 | 79 | 10,000 |
| 4 | 11 | ca-tools | 38 | 8,683 | 38 | 100 | 10,000 |
| 4 | 13 | ca-tools | 46 | 9,022 | 46 | 125 | 10,000 |
| 4 | 16 | ca-tools | 78 | 8,529 | 78 | 159 | 10,000 |
| 4 | 17 | ca-tools | 61 | 8,932 | 61 | 176 | 10,000 |
| 4 | 19 | ca-tools | 46 | 8,796 | 46 | 200 | 10,000 |
| 4 | 23 | ca-tools | 17 | 98,646 | 17 | 266 | 100,000 |
| 4 | 25 | ca-tools | 16 | 98,642 | 16 | 314 | 100,000 |
| 5 | 4 | ca-tools | 1 | 10 | 1 | 601 | 610 |
| 5 | 5 | ca-tools | 1 | 115 | 1 | 782 | 896 |
| 5 | 7 | ca-tools | 37 | 789 | 37 | 695 | 1,483 |
| 5 | 8 | ca-tools | 38 | 1,071 | 38 | 623 | 1,693 |
| 5 | 9 | ca-tools | 34 | 829 | 34 | 819 | 1,647 |
| 5 | 11 | ca-tools | 103 | 9,579 | 103 | 341 | 10,000 |
| 5 | 13 | ca-tools | 53 | 6,148 | 53 | 273 | 10,000 |
| 5 | 16 | ca-tools | 27 | 589 | 27 | 208 | 1,028 |
| 5 | 17 | ca-tools | 38 | 5,725 | 38 | 408 | 10,000 |
| 5 | 19 | ca-tools | 18 | 5,251 | 18 | 471 | 10,000 |
| 5 | 23 | ca-tools | 13 | 5,477 | 13 | 345 | 10,000 |
| 5 | 25 | ca-tools | 17 | 5,334 | 17 | 407 | 10,000 |
| 6 | 7 | ca-tools | 40 | 2,000 | 40 | 801 | 2,800 |
| 6 | 8 | ca-tools | 71 | 2,933 | 71 | 239 | 3,171 |
| 6 | 9 | ca-tools | 66 | 3,739 | 66 | 212 | 3,950 |
| 6 | 11 | ca-tools | 96 | 3,991 | 96 | 210 | 4,200 |
| 6 | 13 | ca-tools | 90 | 4,149 | 90 | 252 | 4,400 |
| 6 | 16 | ca-tools | 110 | 4,532 | 110 | 269 | 4,800 |
| 6 | 17 | ca-tools | 76 | 4,631 | 76 | 270 | 4,900 |
| 6 | 19 | ca-tools | 76 | 4,797 | 76 | 204 | 5,000 |
| 6 | 23 | ca-tools | 88 | 4,735 | 88 | 266 | 5,000 |
| 6 | 25 | ca-tools | 70 | 4,733 | 70 | 268 | 5,000 |

## Rows past every table

| Cell | Ours | Previous best (source) | Status | Certificate |
|---|---|---|---|---|
| CAN(3,k,4), 100001 <= k <= 165004 | 1024 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_1024_3_165004_4/CA_1024_3_165004_4.cphf.txt) |
| CAN(3,k,25), 100001 <= k <= 408135 | 109225 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_109225_3_408135_25/CA_109225_3_408135_25.cphf.txt) |
| CAN(3,k,27), 21197 <= k <= 108143 | 117963 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_117963_3_108143_27/CA_117963_3_108143_27.cphf.txt) |
| CAN(3,k,11), 100001 <= k <= 207658 | 11891 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_11891_3_207658_11/CA_11891_3_207658_11.cphf.txt) |
| CAN(3,k,27), 184914 <= k <= 534635 | 137619 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_137619_3_534635_27/CA_137619_3_534635_27.cphf.txt) |
| CAN(3,k,29), 26131 <= k <= 135764 | 146189 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_146189_3_135764_29/CA_146189_3_135764_29.cphf.txt) |
| CAN(3,k,61), 1580390 <= k <= 9757994 | 1588501 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_1588501_3_9757994_61/CA_1588501_3_9757994_61.cphf.txt) |
| CAN(3,k,5), 100001 <= k <= 154375 | 1685 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_1685_3_154375_5/CA_1685_3_154375_5.cphf.txt) |
| CAN(3,k,29), 184914 <= k <= 684678 | 170549 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_170549_3_684678_29/CA_170549_3_684678_29.cphf.txt) |
| CAN(3,k,13), 100001 <= k <= 134368 | 17485 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_17485_3_134368_13/CA_17485_3_134368_13.cphf.txt) |
| CAN(3,k,31), 31777 <= k <= 167721 | 178591 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_178591_3_167721_31/CA_178591_3_167721_31.cphf.txt) |
| CAN(3,k,32), 33825 <= k <= 184913 | 196448 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_196448_3_184913_32/CA_196448_3_184913_32.cphf.txt) |
| CAN(3,k,31), 184914 <= k <= 860049 | 208351 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_208351_3_860049_31/CA_208351_3_860049_31.cphf.txt) |
| CAN(3,k,67), 2205466 <= k <= 9609828 | 2104939 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_2104939_3_9609828_67/CA_2104939_3_9609828_67.cphf.txt) |
| CAN(3,k,32), 184914 <= k <= 944542 | 229184 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_229184_3_944542_32/CA_229184_3_944542_32.cphf.txt) |
| CAN(3,k,73), 2989885 <= k <= 10998546 | 2722681 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_2722681_3_10998546_73/CA_2722681_3_10998546_73.cphf.txt) |
| CAN(3,k,73), 10998547 <= k <= 19242664 | 2722681 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_2722681_3_19242664_73/CA_2722681_3_19242664_73.cphf.txt) |
| CAN(3,k,16), 100001 <= k <= 322146 | 32656 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_32656_3_322146_16/CA_32656_3_322146_16.cphf.txt) |
| CAN(3,k,17), 100001 <= k <= 106533 | 34289 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_34289_3_106533_17/CA_34289_3_106533_17.cphf.txt) |
| CAN(3,k,37), 647633 <= k <= 1751572 | 354349 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_354349_3_1751572_37/CA_354349_3_1751572_37.cphf.txt) |
| CAN(3,k,37), 265264 <= k <= 647632 | 354349 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_354349_3_647632_37/CA_354349_3_647632_37.cphf.txt) |
| CAN(3,k,7), 100001 <= k <= 137618 | 3703 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_3703_3_137618_7/CA_3703_3_137618_7.cphf.txt) |
| CAN(3,k,19), 100001 <= k <= 160839 | 47899 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_47899_3_160839_19/CA_47899_3_160839_19.cphf.txt) |
| CAN(3,k,41), 383199 <= k <= 643825 | 482201 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_482201_3_643825_41/CA_482201_3_643825_41.cphf.txt) |
| CAN(3,k,8), 100001 <= k <= 116455 | 5048 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_5048_3_116455_8/CA_5048_3_116455_8.cphf.txt) |
| CAN(3,k,43), 454386 <= k <= 2908356 | 556291 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_556291_3_2908356_43/CA_556291_3_2908356_43.cphf.txt) |
| CAN(3,k,97), 8178895 <= k <= 13008618 | 6388129 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_6388129_3_13008618_97/CA_6388129_3_13008618_97.cphf.txt) |
| CAN(3,k,97), 13008619 <= k <= 27995660 | 6388129 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_6388129_3_27995660_97/CA_6388129_3_27995660_97.cphf.txt) |
| CAN(3,k,9), 100001 <= k <= 221150 | 7209 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_7209_3_221150_9/CA_7209_3_221150_9.cphf.txt) |
| CAN(3,k,23), 100001 <= k <= 302074 | 85031 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_85031_3_302074_23/CA_85031_3_302074_23.cphf.txt) |
| CAN(3,k,109), 12349850 <= k <= 15996200 | 9064549 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_9064549_3_15996200_109/CA_9064549_3_15996200_109.cphf.txt) |
| CAN(3,k,109), 15996201 <= k <= 33994950 | 9064549 | none listed | verified, not yet sent | [file](https://github.com/arhancanli/quantavane-covering-arrays/blob/main/wide/CA_9064549_3_33994950_109/CA_9064549_3_33994950_109.cphf.txt) |
