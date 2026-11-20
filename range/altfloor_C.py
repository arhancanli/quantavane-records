# q1 under the three class I floors of RANGE.md section 2 (restricted pair, plain Sherwood, non-Sherwood; fused), s = 3..SMAX
# usage: python3 -I altfloor_C.py classcheck_C.py [SMAX=22]   (v2, 20 Nov: exact integer cube root, start asserted above the boundary)
import importlib.util, sys
spec = importlib.util.spec_from_file_location("cc", sys.argv[1]); cc = importlib.util.module_from_spec(spec)
src = open(sys.argv[1]).read().split("out = {}")[0]          # only the helper functions of classcheck_C.py
exec(compile(src, "cc", "exec"), cc.__dict__)
for s in range(3, int(sys.argv[2]) + 1 if len(sys.argv) > 2 else 23):
    q = 4 ** s; NK = 3 * q ** 3 - 2
    res = []
    for name, rows in (("restricted", lambda Q: 2*Q**3 - Q**2 - 2*(Q - q)), ("sherwood", lambda Q: 2*Q**3 - Q - 2*(Q - q)), ("plain", lambda Q: 2*Q**3 - 1 - 2*(Q - q))):
        Q = cc.icbrt(NK // 2) + 3
        assert rows(Q) > NK, 'start below the boundary'      # v2: the scan must start above the boundary
        while rows(Q) > NK or not cc.is_pp(Q): Q -= 1
        res.append(Q)
    print("s %d q1 restricted %d sherwood %d non-sherwood %d same %s" % (s, res[0], res[1], res[2], len(set(res)) == 1))
