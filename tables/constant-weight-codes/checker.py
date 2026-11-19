#!/usr/bin/env python3
"""Exact checker for constant weight codes: A(n, d, w) >= M.

A binary constant weight code is a set of M words of length n, each with exactly w ones, any two at Hamming distance
at least d. The table (Andries Brouwer's tables of constant weight codes) keeps the largest known M.

    python3 checker.py FILE n d w [--words M]

FILE: one word per line, as n characters 0/1, or as the positions of its ones (1..n, separated by spaces or commas),
or in the layout of Brouwer's files: a first line $BASE=16 (or 2, 8, 10, 32) and then each word as one number in that
base whose binary digits are the word. The program checks every word has length at most n and weight w, the words are
distinct, and runs two separately written routines:
  A  every pair of words: the population count of their XOR is at least d.
  B  the packing form: two words of weight w are at distance 2(w - |a & b|), so the code has distance >= d exactly
     when no set of s = w - ceil(d/2) + 1 positions lies in two words; every s-subset of every word goes into one
     dictionary, which must never see a subset twice.
Exit 0 = both PASS. A is quadratic in M (about a minute at M = 12,000). MIT licence (LICENSE-CODE).
"""
import itertools
import sys


def read(path, n):
    lines = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    base = None
    if lines and lines[0].upper().startswith("$BASE="):
        base = int(lines[0].split("=")[1])
        lines = lines[1:]
    words = []
    for l in lines:
        if base:
            words.append(int(l, base))
        elif set(l) <= {"0", "1"} and len(l) == n:
            words.append(int(l[::-1], 2))
        else:
            pos = [int(x) for x in l.replace(",", " ").split()]
            if min(pos) < 1 or max(pos) > n:
                raise ValueError("a position outside 1..%d" % n)
            words.append(sum(1 << (p - 1) for p in pos))
    if any(x >> n for x in words):
        raise ValueError("a word longer than %d" % n)
    return words


def check_a(words, n, d, w):
    for i, a in enumerate(words):
        for b in words[i + 1:]:
            if bin(a ^ b).count("1") < d:
                return False, "two words at distance %d" % bin(a ^ b).count("1")
    return True, "all %d pairs at distance >= %d" % (len(words) * (len(words) - 1) // 2, d)


def check_b(words, n, d, w):
    s = w - (d + 1) // 2 + 1
    if s <= 0:
        return True, "d > 2w is impossible to violate"
    seen = {}
    for idx, x in enumerate(words):
        pos = [i for i in range(n) if x >> i & 1]
        for sub in itertools.combinations(pos, s):
            if sub in seen:
                return False, "positions %s lie in words %d and %d" % ([p + 1 for p in sub], seen[sub] + 1, idx + 1)
            seen[sub] = idx
    return True, "no %d positions in two words (%d subsets)" % (s, len(seen))


def main(argv):
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    try:
        path, n, d, w = args[0], int(args[1]), int(args[2]), int(args[3])
        words = read(path, n)
    except (IndexError, ValueError) as e:
        print("bad input: %s (see the header of this file)" % (e,))
        return 2
    ok = True
    if "--words" in argv and len(words) != int(argv[argv.index("--words") + 1]):
        print("FAIL  %d words, %s claimed" % (len(words), argv[argv.index("--words") + 1]))
        ok = False
    if any(bin(x).count("1") != w for x in words):
        print("FAIL  a word of weight other than %d" % w)
        ok = False
    if len(set(words)) != len(words):
        print("FAIL  a repeated word")
        ok = False
    for name, fn in (("A", check_a), ("B", check_b)):
        res, why = fn(words, n, d, w)
        print("%s %s  %s" % (name, "PASS" if res else "FAIL", why))
        ok &= res
    print("A(%d,%d,%d) >= %d: %s" % (n, d, w, len(words), "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
