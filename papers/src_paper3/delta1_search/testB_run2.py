"""Test B (vectorized with numpy, double precision for the search; final frontier
entries re-verified at 60-digit mpmath precision).
"""
import time
import numpy as np
from fractions import Fraction
from math import gcd, isqrt
import mpmath as mp
from common import D1, phi, pi

mp.mp.dps = 60

BASIS_NAMES = ["pi","phi","sqrt5","pi^2","phi^2","pi*phi","pi/phi","1/pi","1/phi","log(phi)","log(2)"]
BASIS_MP = {
    "pi": pi, "phi": phi, "sqrt5": mp.sqrt(5), "pi^2": pi**2, "phi^2": phi**2,
    "pi*phi": pi*phi, "pi/phi": pi/phi, "1/pi": 1/pi, "1/phi": 1/phi,
    "log(phi)": mp.log(phi), "log(2)": mp.log(2),
}
BASIS_F = {k: float(v) for k, v in BASIS_MP.items()}
D1_F = float(D1)

MAXDEN = 12
MAXNUM = 48
MAXCOMP = 40

def comp_rat(f: Fraction):
    return 0 if f == 0 else abs(f.numerator) + f.denominator

def rationals_den_le(maxden, maxnum):
    out = [Fraction(0, 1)]
    seen = {Fraction(0, 1)}
    for q in range(1, maxden + 1):
        for p in range(-maxnum, maxnum + 1):
            if p == 0:
                continue
            if gcd(abs(p), q) != 1:
                continue
            fr = Fraction(p, q)
            if fr not in seen:
                seen.add(fr)
                out.append(fr)
    return out

ALL_RATS = rationals_den_le(MAXDEN, MAXNUM)
ALL_COMP = np.array([comp_rat(r) for r in ALL_RATS])
ALL_F = np.array([float(r) for r in ALL_RATS])
# keep only those usable (comp < MAXCOMP, since q0 alone could be 0)
mask = ALL_COMP < MAXCOMP
RATS = [r for r, m in zip(ALL_RATS, mask) if m]
RCOMP = ALL_COMP[mask]
RF = ALL_F[mask]
N = len(RATS)
print("Usable rationals (comp < 40):", N)

# denominators 1..MAXDEN for q0 best-rational-approx (vectorized)
Q0_DENS = np.arange(1, MAXDEN + 1)

def best_q0_vec(remaining):
    """remaining: ndarray of any shape. Returns (q0_num, q0_den, q0_comp, q0_val, err) arrays
    giving best rational approx per element, minimizing |q0_val - remaining| over den 1..MAXDEN,
    numerator unrestricted by MAXNUM (q0 has no explicit numerator cap stated beyond den<=12 rule
    - we still cap |numerator| at MAXNUM for consistency with the rational grammar)."""
    shape = remaining.shape
    flat = remaining.reshape(-1)
    best_err = np.full(flat.shape, np.inf)
    best_num = np.zeros(flat.shape, dtype=np.int64)
    best_den = np.ones(flat.shape, dtype=np.int64)
    for q in Q0_DENS:
        p = np.round(flat * q)
        p = np.clip(p, -MAXNUM, MAXNUM)
        val = p / q
        err = np.abs(val - flat)
        better = err < best_err
        best_err[better] = err[better]
        best_num[better] = p[better].astype(np.int64)
        best_den[better] = q
    return best_num.reshape(shape), best_den.reshape(shape), best_err.reshape(shape)

def reduce_frac(num, den):
    if num == 0:
        return Fraction(0, 1)
    g = gcd(abs(int(num)), int(den))
    return Fraction(int(num)//g, int(den)//g)

results = []  # (comp, relerr, q0str, biname, q1str, bjname, q2str)

t0 = time.time()

# --- single-term: q0 + q1*B_i ---
single_count = 0
for bi in BASIS_NAMES:
    Bi = BASIS_F[bi]
    remaining = D1_F - RF * Bi  # shape (N,)
    num, den, err = best_q0_vec(remaining)
    q0_comp = np.abs(num) + den
    # handle num==0 -> comp 0
    q0_comp = np.where(num == 0, 0, q0_comp)
    total_comp = q0_comp + RCOMP
    ok = total_comp <= MAXCOMP
    idxs = np.where(ok)[0]
    single_count += len(idxs)
    for i in idxs:
        q0 = reduce_frac(num[i], den[i])
        q1 = RATS[i]
        val = float(q0) + float(q1) * Bi
        relerr = abs(val - D1_F) / abs(D1_F)
        results.append((int(total_comp[i]), relerr, str(q0), bi, str(q1), None, None))

t1 = time.time()
print(f"Single-term: {single_count} kept, time {t1-t0:.2f}s")

# --- pair-term: q0 + q1*B_i + q2*B_j, i<j ---
pair_count = 0
for ii in range(len(BASIS_NAMES)):
    for jj in range(ii + 1, len(BASIS_NAMES)):
        bi, bj = BASIS_NAMES[ii], BASIS_NAMES[jj]
        Bi, Bj = BASIS_F[bi], BASIS_F[bj]
        # budget: RCOMP[i1] + RCOMP[i2] < MAXCOMP  (q0 could be 0)
        # build grid but prefilter by comp sum to save memory
        comp_sum = RCOMP[:, None] + RCOMP[None, :]  # (N,N)
        keep = comp_sum < MAXCOMP
        if not keep.any():
            continue
        i1, i2 = np.where(keep)
        q1f = RF[i1]
        q2f = RF[i2]
        remaining = D1_F - q1f * Bi - q2f * Bj
        num, den, err = best_q0_vec(remaining)
        q0_comp = np.abs(num) + den
        q0_comp = np.where(num == 0, 0, q0_comp)
        total_comp = q0_comp + RCOMP[i1] + RCOMP[i2]
        ok = total_comp <= MAXCOMP
        sel = np.where(ok)[0]
        pair_count += len(sel)
        for k in sel:
            q0 = reduce_frac(num[k], den[k])
            q1 = RATS[i1[k]]
            q2 = RATS[i2[k]]
            val = float(q0) + float(q1) * Bi + float(q2) * Bj
            relerr = abs(val - D1_F) / abs(D1_F)
            results.append((int(total_comp[k]), relerr, str(q0), bi, str(q1), bj, str(q2)))

t2 = time.time()
print(f"Pair-term: {pair_count} kept, time {t2-t1:.2f}s")
print(f"Total results: {len(results)}")

# Pareto frontier
results.sort(key=lambda r: (r[0], r[1]))
frontier = []
best_err_so_far = None
for r in results:
    comp_, relerr = r[0], r[1]
    if best_err_so_far is None or relerr < best_err_so_far - 1e-18:
        frontier.append(r)
        best_err_so_far = relerr

# Re-verify frontier entries at 60-digit precision
def mp_val(q0s, biname, q1s, bjname, q2s):
    q0 = mp.mpf(Fraction(q0s).numerator) / mp.mpf(Fraction(q0s).denominator)
    v = q0
    if biname is not None:
        q1f = Fraction(q1s)
        v += (mp.mpf(q1f.numerator)/mp.mpf(q1f.denominator)) * BASIS_MP[biname]
    if bjname is not None:
        q2f = Fraction(q2s)
        v += (mp.mpf(q2f.numerator)/mp.mpf(q2f.denominator)) * BASIS_MP[bjname]
    return v

frontier_verified = []
for r in frontier:
    comp_, relerr_f, q0, bi, q1, bj, q2 = r
    v = mp_val(q0, bi, q1, bj, q2)
    relerr_mp = abs(v - D1) / abs(D1)
    frontier_verified.append((comp_, float(relerr_mp), q0, bi, q1, bj, q2))

with open("testB_frontier.txt", "w") as f:
    f.write(f"Single-term kept: {single_count}\n")
    f.write(f"Pair-term kept: {pair_count}\n")
    f.write(f"Total candidate results (pre-Pareto): {len(results)}\n")
    f.write("Pareto frontier (mpmath 60-digit verified relerr):\n")
    for comp_, relerr_mp, q0, bi, q1, bj, q2 in frontier_verified:
        expr = f"{q0} + {q1}*{bi}"
        if bj is not None:
            expr += f" + {q2}*{bj}"
        f.write(f"comp={comp_:3d}  relerr={relerr_mp:.3e}  {expr}\n")
    # find where candidate C = 15/8 + 4*pi - phi sits
    f.write("\nKnown candidate C = 15/8 + 4*pi - 1*phi, complexity 30, relerr vs D1:\n")
    from common import C
    f.write(f"{float(abs(C-D1)/abs(D1)):.3e}\n")

print(f"Frontier size (verified): {len(frontier_verified)}")
for row in frontier_verified[:15]:
    print(row)

t3 = time.time()
print(f"Total Test B time: {t3-t0:.2f}s")
