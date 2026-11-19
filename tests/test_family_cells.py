"""Tests of the family-cell API.

  python3 -m unittest discover -s tests -v              (about two minutes)
  QV_SLOW=1 python3 -m unittest discover -s tests -v    (adds the exhaustive check at q = 16)

The main test draws 1,000 cells uniformly from the 227,658,593,685,328,988,782,135,588 counted Theorem K cells and checks each
with both routines: the construction routine (field model A, qvrecords.theorem_k) and the independent check
(field model B, qvrecords.verify), which must agree on every column they both compute, both find every sampled triple
covered, and pass the certificate.
"""
import csv
import os
import random
import sys
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import qvrecords                                     # noqa: E402
from qvrecords import theorem_k, verify             # noqa: E402

SEED = 20261119


class FamilyCells(unittest.TestCase):

    def test_counts_match_families_csv(self):
        """the API's counted cells per q are the counts in data/families.csv"""
        with open(os.path.join(qvrecords.DATA, "families.csv"), newline="") as fh:
            rows = list(csv.DictReader(fh))
        declared = {}
        for r in rows:
            q = 16 if r["id"] == "F-K-q16" else 4 ** int(r["id"][len("F-K-s"):])
            declared[q] = (int(r["cells"]), int(r["k_from"]), int(r["k_to"]))
        total, by_q = qvrecords.family_cells()
        self.assertEqual(set(declared), set(by_q))
        for q, (n, lo, hi) in declared.items():
            self.assertEqual(by_q[q], n, q)
            self.assertEqual(theorem_k.counted_range(q), (lo, hi), q)
            if q != 16:
                self.assertEqual(verify.counted_range_b(q), (lo, hi), q)
        self.assertEqual(total, sum(n for n, _, _ in declared.values()))
        self.assertEqual(total, 227658593685328988782135588)

    def test_1000_random_cells_both_routines(self):
        total, _ = qvrecords.family_cells()
        R = random.Random(SEED)
        t0, triples, by_q = time.time(), 0, {}
        for n in range(1000):
            i = R.randrange(total)
            q, k = qvrecords.cell_at(i)
            by_q[q] = by_q.get(q, 0) + 1
            value, cert = qvrecords.cell("theorem-k", {"index": i})
            self.assertEqual((cert["q"], cert["k"], value), (q, k, 3 * q ** 3 - 2))
            self.assertTrue(cert["counted"])
            A = theorem_k.family(q)
            B = verify.FamilyB(q, cert["field"]["modulus"], cert["field"]["non_cube"])
            for _ in range(3):                   # the same triples through both routines
                idx = set()
                while len(idx) < 3:
                    idx.add(R.randrange(k))
                ca = [A.column(j) for j in idx]
                cb = [B.column(j) for j in idx]
                self.assertEqual([[list(e) for e in c] for c in ca], cb, (q, k, idx))
                self.assertTrue(A.covered(ca), (q, k, idx))
                self.assertTrue(B.covers(cb), (q, k, idx))
                triples += 1
            ok, rep = qvrecords.check(cert, samples=4, hard=2, seed=n)
            self.assertTrue(ok, (q, k, [x for x in rep if x[0] == "FAIL"]))
        print("\n  1,000 cells (by q: %s), %d triples through both routines plus check(): %.1f s"
              % (", ".join("%d: %d" % kv for kv in sorted(by_q.items())), triples, time.time() - t0))

    def test_every_q_stratified(self):
        """20 random counted cells at each of the 21 values of q, so the small fields are tested too"""
        R = random.Random(SEED + 1)
        for q, lo, hi in qvrecords._blocks():
            for n in range(20):
                k = R.randint(lo, hi)
                value, cert = qvrecords.cell("theorem-k", {"q": q, "k": k})
                self.assertTrue(cert["counted"])
                ok, rep = qvrecords.check(cert, samples=6, hard=3, seed=n)
                self.assertTrue(ok, (q, k, [x for x in rep if x[0] == "FAIL"]))

    def test_range_ends(self):
        """first and last counted cell of every q pass; the cell just below a range passes as not counted"""
        for q, lo, hi in qvrecords._blocks():
            for k in (lo, hi):
                ok, rep = qvrecords.check(qvrecords.cell("theorem-k", {"q": q, "k": k})[1], samples=4, hard=2)
                self.assertTrue(ok, (q, k, rep))
            v, cert = qvrecords.cell("theorem-k", {"q": q, "k": lo - 1})
            self.assertFalse(cert["counted"])
            self.assertTrue(qvrecords.check(cert, samples=4, hard=2)[0], (q, lo - 1))

    def test_index_roundtrip(self):
        total, by_q = qvrecords.family_cells()
        self.assertEqual(qvrecords.cell_at(0), (16, 379))
        self.assertEqual(qvrecords.cell_at(167), (16, 546))
        self.assertEqual(qvrecords.cell_at(168), (64, 5404))
        self.assertEqual(qvrecords.cell_at(total - 1), (4 ** 22, 2 * (4 ** 44 + 4 ** 22 + 1)))
        with self.assertRaises(IndexError):
            qvrecords.cell_at(total)

    def test_red_controls(self):
        """damaged certificates fail the check"""
        _, good = qvrecords.cell("theorem-k", {"q": 256, "k": 100000})
        self.assertTrue(qvrecords.check(good)[0])

        def damaged(**changes):
            c = {k: (dict(v) if isinstance(v, dict) else v) for k, v in good.items()}
            for key, val in changes.items():
                if key in ("modulus", "non_cube"):
                    c["field"][key] = val
                else:
                    c[key] = val
            return c
        m = good["field"]["m"]
        F = verify.FieldB(m, good["field"]["modulus"])
        cube = F.mul(F.mul(7, 7), 7)
        col = {j: [list(e) for e in v] for j, v in good["sample_columns"].items()}
        j0 = next(iter(col))
        col[j0][1][0] ^= 1
        for bad in (damaged(value=good["value"] - 1),
                    damaged(statement="CAN(3, 100000, 256) <= 50331645"),
                    damaged(modulus=(1 << m) | 1),                  # t^m + 1: divisible by t + 1
                    damaged(non_cube=cube),
                    damaged(sample_columns=col),
                    damaged(k=5000),                                # below the counted range, still marked counted
                    damaged(q=128),
                    damaged(k=2 * (256 ** 2 + 256 + 1) + 1)):
            self.assertFalse(qvrecords.check(bad)[0], bad)

    def test_covering_test_has_teeth(self):
        """model B's covering test rejects a family with one column repeated (exhaustive at q = 4)"""
        B = verify.FamilyB(4, 0b10011, 2)
        cols = [B.column(j) for j in range(42)]
        self.assertTrue(all(B.covers([cols[0], cols[i], cols[j]]) for i in range(1, 42) for j in range(i + 1, 42)))
        cols[7] = cols[8]
        self.assertFalse(B.covers([cols[7], cols[8], cols[0]]))

    def test_exhaustive_q4(self):
        _, cert = qvrecords.cell("theorem-k", {"q": 4, "k": 42})
        ok, rep = qvrecords.check(cert, exhaustive=True)
        self.assertTrue(ok, rep)

    def test_q4_covering_array_by_definition(self):
        """expand CA(190; 3, 42, 4) by Sherwood's construction and check every 3 columns see all 64 triples"""
        q, k = 4, 42
        rows = 3 * q ** 3 - 2
        A = [[theorem_k.ca_entry(q, k, r, c) for c in range(k)] for r in range(rows)]
        cols = list(zip(*A))
        for a in range(k):
            for b in range(a + 1, k):
                for c in range(b + 1, k):
                    seen = {(x * q + y) * q + z for x, y, z in zip(cols[a], cols[b], cols[c])}
                    self.assertEqual(len(seen), q ** 3, (a, b, c))

    @unittest.skipUnless(os.environ.get("QV_SLOW"), "set QV_SLOW=1")
    def test_exhaustive_q16(self):
        _, cert = qvrecords.cell("theorem-k", {"q": 16, "k": 546})
        t0 = time.time()
        ok, rep = qvrecords.check(cert, exhaustive=True)
        print("\n  q = 16 exhaustive: %s (%.0f s)" % (rep[-1][1], time.time() - t0))
        self.assertTrue(ok, rep)


if __name__ == "__main__":
    unittest.main()
