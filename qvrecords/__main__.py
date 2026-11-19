"""python3 -m qvrecords: query records and family cells from the command line.

  python3 -m qvrecords summary                      counts: records by status and field, family cells
  python3 -m qvrecords record ID                    one row of data/records.csv
  python3 -m qvrecords ca T V K                     the covering-array record holding CAN(T, K, V)
  python3 -m qvrecords cell theorem-k Q K [--check] [--exhaustive]
  python3 -m qvrecords index I [--check]            the I-th counted family cell
"""
import collections
import json
import sys

from . import ca_record, cell, check, family_cells, record, records


def _show(value, cert, do_check, exhaustive=False):
    print(json.dumps(cert, indent=1))
    if do_check:
        ok, rep = check(cert, exhaustive=exhaustive)
        for st, what in rep:
            print("%s  %s" % (st, what))
        print("CHECK", "PASS" if ok else "FAIL")
        return 0 if ok else 1
    return 0


def main(a):
    if not a or a[0] == "summary":
        rs = records()
        by = collections.Counter((r["field"], r["status"]) for r in rs)
        cells = collections.Counter()
        for r in rs:
            cells[r["status"]] += int(r["cells"])
        print("records: %d rows covering %d cells" % (len(rs), sum(cells.values())))
        for st, n in sorted(collections.Counter(r["status"] for r in rs).items()):
            print("  %-18s %6d rows %12d cells" % (st, n, cells[st]))
        for (f, st), n in sorted(by.items()):
            print("  %-28s %-18s %5d" % (f, st, n))
        tot, by_q = family_cells()
        print("family cells (theorem-k, counted apart): %d" % tot)
        for q, n in by_q.items():
            print("  q = %-12d %d" % (q, n))
        return 0
    if a[0] == "record":
        r = record(a[1])
        print(json.dumps(r, indent=1) if r else "no record %s" % a[1])
        return 0 if r else 1
    if a[0] == "ca":
        r = ca_record(int(a[1]), int(a[2]), int(a[3]))
        print(json.dumps(r, indent=1) if r else "no record holds CAN(%s, %s, %s)" % (a[1], a[3], a[2]))
        return 0 if r else 1
    if a[0] == "cell":
        value, cert = cell(a[1], {"q": int(a[2]), "k": int(a[3])})
        return _show(value, cert, "--check" in a, "--exhaustive" in a)
    if a[0] == "index":
        value, cert = cell("theorem-k", {"index": int(a[1])})
        return _show(value, cert, "--check" in a)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
