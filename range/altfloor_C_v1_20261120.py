# q1 under the three class I floors of RANGE.md section 2 (restricted pair, plain Sherwood, non-Sherwood; fused), s = 3..22
import importlib.util, sys
spec = importlib.util.spec_from_file_location("cc", sys.argv[1]); cc = importlib.util.module_from_spec(spec)
src = open(sys.argv[1]).read().split("out = {}")[0]          # only the helper functions of classcheck_C.py
exec(compile(src, "cc", "exec"), cc.__dict__)
for s in range(3, 23):
    q = 4 ** s; NK = 3 * q ** 3 - 2
    res = []
    for name, rows in (("restricted", lambda Q: 2*Q**3 - Q**2 - 2*(Q - q)), ("sherwood", lambda Q: 2*Q**3 - Q - 2*(Q - q)), ("plain", lambda Q: 2*Q**3 - 1 - 2*(Q - q))):
        Q = int((NK / 2) ** (1 / 3)) + 3
        while rows(Q) > NK or not cc.is_pp(Q): Q -= 1
        res.append(Q)
    print("s %d q1 restricted %d sherwood %d non-sherwood %d same %s" % (s, res[0], res[1], res[2], len(set(res)) == 1))
