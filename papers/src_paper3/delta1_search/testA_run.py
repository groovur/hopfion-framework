import time
import math
import numpy as np
from fractions import Fraction
from testA_engine import enumerate_exprs, eval_expr_float
import mpmath as mp
from common import D1, phi, pi

MAXC = 30
TOL = 6.15e-5

t0 = time.time()
res = enumerate_exprs(MAXC)
t1 = time.time()
print(f"Enumerated {len(res)} expressions up to complexity {MAXC} in {t1-t0:.2f}s")

D1_f = float(mp.nstr(D1, 20))

vals = np.empty(len(res), dtype=np.float64)
for i, (a, b, c, tot) in enumerate(res):
    vals[i] = eval_expr_float(a, b, c)
t2 = time.time()
print(f"Evaluated in {t2-t1:.2f}s")

relerr = np.abs(vals - D1_f) / abs(D1_f)
hit_idx = np.where(relerr <= TOL)[0]
print(f"Number of hits within tol {TOL}: {len(hit_idx)}")

hits = []
for i in hit_idx:
    a, b, c, tot = res[i]
    hits.append((tot, a, b, c, relerr[i]))
hits.sort(key=lambda x: x[0])

candidate_found = False
with open("testA_hits.txt", "w") as f:
    f.write(f"MAXC={MAXC} TOL={TOL}\n")
    f.write(f"Total enumerated: {len(res)}\n")
    f.write(f"Total hits: {len(hits)}\n\n")
    for tot, a, b, c, err in hits:
        line = f"comp={tot:3d}  a={str(a):>8s} b={str(b):>8s} c={str(c):>8s}  relerr={err:.3e}"
        f.write(line + "\n")
        if a == Fraction(15,8) and b == Fraction(4,1) and c == Fraction(-1,1):
            candidate_found = True

print("Candidate (15/8, 4, -1) found among hits:", candidate_found)

# save sorted vals array for calibration
np.save("testA_vals_sorted.npy", np.sort(vals))
np.save("testA_vals_raw.npy", vals)

t3 = time.time()
print(f"Total time so far: {t3-t0:.2f}s")
