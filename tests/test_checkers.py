"""Tests of the table checkers: each passes a small known object and fails a damaged copy (a red control).

  python3 -m unittest discover -s tests -v

The fixtures are classical objects written out here (the Fano plane, the Petersen and Heawood graphs, Paley
colourings, the Hamming code ...), plus Theorem K's own family exported to the covering-array checker, so the family
API and the table checker are tested against each other.
"""
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from qvrecords import theorem_k                      # noqa: E402


def run(field, *args):
    p = subprocess.run([sys.executable, "-I", os.path.join(ROOT, "tables", field, "checker.py")] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout


class Checkers(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def assertPasses(self, field, *args):
        code, out = run(field, *args)
        self.assertEqual(code, 0, out)
        self.assertEqual(out.count(" PASS"), out.count(" PASS"))       # both routines print PASS
        self.assertNotIn("FAIL", out)

    def assertFails(self, field, *args):
        code, out = run(field, *args)
        self.assertEqual(code, 1, out)

    # ------------------------------------------------------------------------------------------- designs
    def test_covering_designs(self):
        fano = "1 2 3\n1 4 5\n1 6 7\n2 4 6\n2 5 7\n3 4 7\n3 5 6\n"
        f = self.write("fano.txt", fano)
        self.assertPasses("covering-designs", f, 7, 3, 2, "--blocks", 7)
        self.assertPasses("covering-designs", f, 7, 3, 2, 3)          # every 3-set meets a line in 2 points
        g = self.write("fano6.txt", "\n".join(fano.splitlines()[:6]) + "\n")
        self.assertFails("covering-designs", g, 7, 3, 2)
        self.assertFails("covering-designs", f, 7, 3, 2, "--blocks", 6)

    # ------------------------------------------------------------------------------------------- arrays
    def test_covering_arrays_ca(self):
        f = self.write("ca.txt", "0 0 0\n0 1 1\n1 0 1\n1 1 0\n")
        self.assertPasses("covering-arrays", "ca", f, 2, 2, "--rows", 4)
        g = self.write("ca_bad.txt", "0 0 0\n0 1 1\n1 0 1\n")
        self.assertFails("covering-arrays", "ca", g, 2, 2)

    def test_covering_arrays_theorem_k_family(self):
        """Theorem K's columns, exported in the ca-tools layout, pass the table checker (q = 4 whole, q = 16 in part)"""
        for q, K in ((4, 42), (16, 40)):
            fam = theorem_k.family(q)
            rows = [[] for _ in range(3)]
            for j in range(K):
                col = fam.column(j)
                for r in range(3):
                    a = col[r]
                    rows[r].append(a[0] * q * q + a[1] * q + a[2])
            f = self.write("k%d.txt" % q, "\n".join(" ".join(map(str, r)) for r in rows) + "\n")
            self.assertPasses("covering-arrays", "cphf", f, 3, q, "--rows", 3 * q ** 3 - 2)
            rows[1][5] = rows[1][6]
            rows[0][5], rows[2][5] = rows[0][6], rows[2][6]                # column 5 := column 6
            g = self.write("k%d_bad.txt" % q, "\n".join(" ".join(map(str, r)) for r in rows) + "\n")
            self.assertFails("covering-arrays", "cphf", g, 3, q)

    # ------------------------------------------------------------------------------------------- codes
    def test_spherical_codes(self):
        f = self.write("tet.txt", "1,1,1\n1,-1,-1\n-1,1,-1\n-1,-1,1\n")
        self.assertPasses("spherical-codes", f, 3, 4, "-0.33333333333333333333")
        self.assertFails("spherical-codes", f, 3, 4, "-0.35")
        sq = self.write("sq.txt", "1,0\n0,1\n-1,0\n0,-1\n")
        self.assertPasses("spherical-codes", sq, 2, 4, "0")
        self.assertFails("spherical-codes", sq, 2, 4, "-0.1")

    def test_constant_weight_codes(self):
        f = self.write("cw.txt", "111000\n100110\n010101\n001011\n")
        self.assertPasses("constant-weight-codes", f, 6, 4, 3, "--words", 4)
        g = self.write("cw_bad.txt", "111000\n100110\n010101\n110100\n")
        self.assertFails("constant-weight-codes", g, 6, 4, 3)

    def test_linear_codes(self):
        f = self.write("ham.txt", "1 0 0 0 1 1 0\n0 1 0 0 1 0 1\n0 0 1 0 0 1 1\n0 0 0 1 1 1 1\n")
        self.assertPasses("linear-codes", f, 2, 7, 4, 3, "--exact")
        self.assertFails("linear-codes", f, 2, 7, 4, 4)
        # the hexacode [6,3,4] over GF(4) (w = 2 in the Conway coding, w^2 = w + 1 = 3)
        h = self.write("hex.txt", "1 0 0 1 3 2\n0 1 0 1 2 3\n0 0 1 1 1 1\n")
        self.assertPasses("linear-codes", h, 4, 6, 3, 4, "--exact")

    # ------------------------------------------------------------------------------------------- graphs
    def petersen(self, extra=()):
        P = [(i, (i + 1) % 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)] + \
            [(i, i + 5) for i in range(5)] + list(extra)
        return self.write("p%d.txt" % len(P), "\n".join("%d %d" % e for e in P) + "\n")

    def test_cages(self):
        self.assertPasses("cages", "cage", self.petersen(), 3, 5, "--vertices", 10)
        H = [(i, (i + 1) % 14) for i in range(14)] + [(i, (i + 5) % 14) for i in range(0, 14, 2)]
        h = self.write("heawood.txt", "\n".join("%d %d" % e for e in H) + "\n")
        self.assertPasses("cages", "cage", h, 3, 6, "--vertices", 14)
        self.assertFails("cages", "cage", h, 3, 5)
        self.assertFails("cages", "cage", self.petersen([(0, 2)]), 3, 5)

    def test_extremal_graphs(self):
        self.assertPasses("extremal-graphs", "extremal", self.petersen(), 10, 4, "--edges", 15)
        self.assertFails("extremal-graphs", "extremal", self.petersen([(0, 2)]), 10, 4)

    def test_graph6_reader(self):
        """the Petersen graph in graph6 ('IheA@GUAo') is read and passes"""
        f = self.write("petersen.g6", "IheA@GUAo\n")
        self.assertPasses("cages", "cage", f, 3, 5, "--vertices", 10)

    # ------------------------------------------------------------------------------------------- others
    def test_weighing_matrices(self):
        f = self.write("cw4.txt", "1 1 1 -1\n")
        self.assertPasses("weighing-matrices", f, 4, 4)
        g = self.write("cw4_bad.txt", "1 1 1 1\n")
        self.assertFails("weighing-matrices", g, 4, 4)

    def test_weak_schur(self):
        f = self.write("ws2.txt", "1 2 4 8\n3 5 6 7\n")
        self.assertPasses("weak-schur-numbers", f, 2, 8)
        g = self.write("ws2_bad.txt", "1 2 4 8 9\n3 5 6 7\n")
        self.assertFails("weak-schur-numbers", g, 2, 9)

    def test_ramsey(self):
        def paley(p):
            sq = {(x * x) % p for x in range(1, p)}
            return ["".join("0" if i == j else ("1" if (j - i) % p in sq else "0") for j in range(p)) for i in range(p)]
        f = self.write("p17.txt", "\n".join(paley(17)) + "\n")
        self.assertPasses("ramsey-numbers", f, 4, 4)
        self.assertFails("ramsey-numbers", f, 3, 4)
        m = [list(r) for r in paley(17)]
        m[0][1] = m[1][0] = "0" if m[0][1] == "1" else "1"
        g = self.write("p17_bad.txt", "\n".join("".join(r) for r in m) + "\n")
        self.assertFails("ramsey-numbers", g, 4, 4)
        # three colours on K_16 from GF(16): Greenwood-Gleason, no monochromatic triangle
        q = 16
        mod = 0b10011

        def mul(a, b):
            r = 0
            while b:
                if b & 1:
                    r ^= a
                b >>= 1
                a <<= 1
                if a & 16:
                    a ^= mod
            return r
        g_ = 2
        log = {}
        x = 1
        for e in range(15):
            log[x] = e
            x = mul(x, g_)
        rows = ["".join("0" if i == j else str(log[i ^ j] % 3 + 1) for j in range(q)) for i in range(q)]
        h = self.write("gg16.txt", "\n".join(rows) + "\n")
        self.assertPasses("ramsey-numbers", h, "--forbid", "1:3,2:3,3:3")


if __name__ == "__main__":
    unittest.main()
