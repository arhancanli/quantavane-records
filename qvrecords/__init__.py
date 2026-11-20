"""qvrecords: Quantavane's records, queryable.

Two kinds of entry, kept apart:

- records: the single objects in data/records.csv (each a code, design, graph, colouring, array, schedule ... that
  set a table's best value), loaded by records(), record(id) and ca_record(t, v, k);
- family cells: cells settled by a proved construction over a whole range, computed on demand, never stored.
  cell(family, params) returns (value, certificate); check(certificate) verifies it with separately written code.

    >>> from qvrecords import cell, check
    >>> value, cert = cell("theorem-k", {"q": 64, "k": 6000})
    >>> value
    786430
    >>> check(cert)[0]
    True

Families served: "theorem-k" (CAN(3, k, q) <= 3q^3 - 2 for q = 4^s, any k <= 2(q^2+q+1)). family_cells() counts
the cells the lab counts as records and cell_at(i) names the i-th of them.
"""
import csv
import os

from . import theorem_k
from .verify import check as _check_k

__all__ = ["cell", "check", "families", "family_cells", "cell_at", "records", "record", "ca_record"]
__version__ = "0.1.0"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

S_MIN, S_MAX = 3, 26     # the declared range of Theorem K: q = 4^s for 3 <= s <= 26, plus q = 16 (s <= 16 until 20 Nov 2026)


def families():
    return {"theorem-k": {
        "statement": "CAN(3, k, q) <= 3q^3 - 2 for q = 4^s and 1 <= k <= 2(q^2 + q + 1)",
        "params": {"q": "a power of 4, at least 4", "k": "1 .. 2(q^2 + q + 1)"},
        "counted": "q = 16 with 379 <= k <= 546 (keeper-listed), and q = 4^s, 3 <= s <= 26, with k in counted_range(q)",
        "source": theorem_k.SOURCE}}


def cell(family, params):
    """(value, certificate) of one cell. params: {"q": q, "k": k}, or {"index": i} for the i-th counted cell."""
    if family != "theorem-k":
        raise ValueError("unknown family %r; served: %s" % (family, ", ".join(families())))
    if "index" in params:
        q, k = cell_at(int(params["index"]))
    else:
        q, k = int(params["q"]), int(params["k"])
    return theorem_k.certificate(q, k)


def check(cert, **kw):
    """(ok, report) for a certificate returned by cell(); kw: samples, hard, seed, exhaustive"""
    if not isinstance(cert, dict) or cert.get("family") != "theorem-k":
        return False, [("FAIL", "not a theorem-k certificate")]
    return _check_k(cert, **kw)


def _blocks():
    """the counted cells in their global order: q = 16 first, then s = 3 .. 26"""
    out = [(16,) + theorem_k.counted_range(16)]
    for s in range(S_MIN, S_MAX + 1):
        q = 4 ** s
        out.append((q,) + theorem_k.counted_range(q))
    return out


def family_cells():
    """number of counted Theorem K cells, by q and in all"""
    by_q = {q: hi - lo + 1 for q, lo, hi in _blocks()}
    return sum(by_q.values()), by_q


def cell_at(i):
    """(q, k) of the i-th counted cell, 0 <= i < family_cells()[0]"""
    if i < 0:
        raise IndexError("negative index")
    for q, lo, hi in _blocks():
        n = hi - lo + 1
        if i < n:
            return q, lo + i
        i -= n
    raise IndexError("past the last counted cell")


# ------------------------------------------------------------------------------------------- single records
_RECORDS = None


def records(path=None):
    """every row of data/records.csv, as dicts"""
    global _RECORDS
    if _RECORDS is None or path:
        with open(path or os.path.join(DATA, "records.csv"), newline="") as fh:
            _RECORDS = list(csv.DictReader(fh))
    return _RECORDS


def record(rid):
    return next((r for r in records() if r["id"] == rid), None)


def ca_record(t, v, k):
    """the covering-array record row whose range holds the cell CAN(t, k, v), or None"""
    import re
    for r in records():
        if r["field"] != "covering-arrays" or not r["k_from"]:
            continue
        m = re.match(r"CAN\((\d+),k,(\d+)\)", r["parameters"])
        if m and int(m.group(1)) == t and int(m.group(2)) == v and int(r["k_from"]) <= k <= int(r["k_to"]):
            return r
    return None
