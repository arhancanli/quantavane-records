# Submitting to the Ramsey numbers table

1. Check your construction yourself with [checker.py](checker.py) (see [CHECKER.md](CHECKER.md)). It is pure Python 3 with no other packages.
2. Open an issue at https://github.com/arhancanli/quantavane-records/issues titled `ramsey-numbers: <cell> <value>`, and attach the certificate: the n x n colour matrix, one line per vertex. Say how the construction was found and how you would like to be credited.
3. Within 72 hours we run both checkers on it and post their output in the issue.
4. If both pass and your value beats the one in [TABLE.md](TABLE.md), we update the table, with your name and your source, in the next commit. If it beats the keeper's table too, send it to the keeper as well: the keeper's table is the reference, and we do not send other people's results on their behalf unless asked.
5. If a check fails, the issue says which routine failed and why (for example the first uncovered tuple or the first short cycle), so the certificate can be fixed.

Corrections to our own entries are welcome the same way. Certificates are published here under CC BY 4.0, with credit to their authors.
