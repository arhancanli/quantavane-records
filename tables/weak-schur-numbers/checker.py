#!/usr/bin/env python3
"""Exact checker for weak Schur numbers: WS(k) >= n.

A partition of {1, ..., n} into k sets is weakly sum-free when no set holds x, y, z with x + y = z and x != y.
WS(k) is the largest n for which such a partition exists; the table (Leaps in Bounds) keeps lower bounds.

    python3 checker.py FILE k n

FILE: k lines, each the numbers of one set separated by spaces or commas; or the colours of 1, 2, ..., n as one
string of base-36 digits 0..k-1 (line breaks ignored); or one line of n colours separated by spaces. The program checks the sets partition 1..n and runs two separately written routines:
  A  bitsets: for each set S as an integer M and each x in S, (M >> x) & M has a bit at y exactly when y and x + y
     are both in S; it must be empty apart from y = x.
  B  written separately: the sumset of each set from the square of the polynomial sum_{x in S} X^x, computed as one
     big-integer product with X = 2^24; the coefficient of X^z for z in S must be 0, after removing the pair x = y.
Exit 0 = both PASS. Both handle n in the millions (about a minute at n = 818,354). MIT licence (LICENSE-CODE).
"""
import sys

def read(path, k, n):
    raw = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    lines = [l.replace(",", " ").split() for l in raw]
    joined = "".join(raw)
    if all(len(l) == 1 for l in lines) and len(joined) == n:      # a colour string, one base-36 digit per number
        col = [int(c, 36) for c in joined]
        sets = [[i + 1 for i, c in enumerate(col) if c == j] for j in range(k)]
    elif len(lines) == 1 and len(lines[0]) == n:
        col = [int(x) for x in lines[0]]
        sets = [[i + 1 for i, c in enumerate(col) if c == j] for j in range(k)]
    else:
        sets = [sorted(int(x) for x in l) for l in lines]
    if len(sets) != k:
        raise ValueError("%d sets, %d claimed" % (len(sets), k))
    allv = sorted(x for s in sets for x in s)
    if allv != list(range(1, n + 1)):
        raise ValueError("the sets do not partition 1..%d" % n)
    return sets


def check_a(sets, n):
    for j, S in enumerate(sets):
        M = 0
        for x in S:
            M |= 1 << x
        for x in S:
            hit = (M >> x) & M & ~(1 << x)
            if hit:
                y = (hit & -hit).bit_length() - 1
                return False, "set %d holds %d + %d = %d" % (j + 1, x, y, x + y)
    return True, "every set weakly sum-free (bitsets)"


def check_b(sets, n):
    """sumsets by polynomial squaring: P(X) = sum_{x in S} X^x with X = 2^24 (Kronecker substitution); the
    coefficient of X^z in P^2 counts the ordered pairs (x, y) of S with x + y = z, the pair x = y once"""
    for j, S in enumerate(sets):
        if len(S) >= 1 << 23:
            return None, "set too large for 24-bit coefficients"
        buf = bytearray(3 * (n + 1))
        for x in S:
            buf[3 * x] = 1
        P = int.from_bytes(bytes(buf), "little")
        sq = (P * P).to_bytes(3 * (2 * n + 2), "little")
        inS = bytes(buf)
        for z in S:
            c = sq[3 * z] | sq[3 * z + 1] << 8 | sq[3 * z + 2] << 16
            if z % 2 == 0 and inS[3 * (z // 2)]:
                c -= 1
            if c:
                x = next(x for x in S if x < z - x and inS[3 * (z - x)])
                return False, "set %d holds %d + %d = %d" % (j + 1, x, z - x, z)
    return True, "every set weakly sum-free (sumsets by polynomial squaring)"


def main(argv):
    try:
        path, k, n = argv[0], int(argv[1]), int(argv[2])
        sets = read(path, k, n)
    except (IndexError, ValueError) as e:
        print("bad input: %s (see the header of this file)" % (e,))
        return 2
    ok = True
    for name, fn in (("A", check_a), ("B", check_b)):
        res, why = fn(sets, n)
        print("%s %s  %s" % (name, {True: "PASS", False: "FAIL", None: "SKIP"}[res], why))
        ok &= res is not False
    print("WS(%d) >= %d: %s" % (k, n, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
