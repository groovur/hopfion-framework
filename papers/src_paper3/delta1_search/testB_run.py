"""Test B: E = q0 + q1*B_i + q2*B_j, basis of 11 constants, q with den<=12, |num|<=48,
total complexity <= 40. For each (B_i,B_j,q1,q2), compute needed q0 = D1 - q1*B_i - q2*B_j,
find best rational approx with den<=12 via continued fractions, check total complexity and error.
"""
import time
from fractions import Fraction
import mpmath as mp
from common import D1, phi, pi

mp.mp.dps = 60

BASIS = {
    "pi": pi,
    "phi": phi,
    "sqrt5": mp.sqrt(5),
    "pi^2": pi**2,
    "phi^2": phi**2,
    "pi*phi": pi*phi,
    "pi/phi": pi/phi,
    "1/pi": 1/pi,
    "1/phi": 1/phi,
    "log(phi)": mp.log(phi),
    "log(2)": mp.log(2),
}
NAMES = list(BASIS.keys())

MAXDEN = 12
MAXNUM = 48
MAXCOMP = 40

def comp_rat(f: Fraction):
    if f == 0:
        return 0
    return abs(f.numerator) + f.denominator

def rationals_den_le(maxden, maxnum):
    """All reduced fractions p/q, 1<=q<=maxden, |p|<=maxnum, plus 0."""
    from math import gcd
    out = [Fraction(0,1)]
    seen = {Fraction(0,1)}
    for q in range(1, maxden+1):
        for p in range(-maxnum, maxnum+1):
            if p == 0:
                continue
            if gcd(abs(p), q) != 1:
                continue
            fr = Fraction(p, q)
            if fr not in seen:
                seen.add(fr)
                out.append(fr)
    return out

RATS = rationals_den_le(MAXDEN, MAXNUM)
print("Number of candidate rationals for q1,q2:", len(RATS))

def best_rational_approx(x_mpf, maxden):
    """Find best rational p/q with q<=maxden approximating x_mpf (mpmath float), via mpmath's
    continued fraction / Stern-Brocot search. Returns Fraction."""
    # Use mpmath's pslq-free approach: try all denominators 1..maxden, round numerator.
    best = None
    best_err = None
    for q in range(1, maxden+1):
        p = mp.nint(x_mpf * q)
        cand = Fraction(int(p), q)
        err = abs(mp.mpf(cand.numerator)/cand.denominator - x_mpf)
        if best_err is None or err < best_err:
            best_err = err
            best = cand
    return best

t0 = time.time()
results = []  # (total_comp, relerr, q0, name_i, q1, name_j, q2)

# term configurations: 0 terms (just q0 alone - trivial, skip since that's just rational != D1 unless huge denom)
# 1 basis term: q0 + q1*B_i
# 2 basis terms: q0 + q1*B_i + q2*B_j (i<j, i can equal... "at most 2 distinct basis elements")

pair_count = 0
single_count = 0

def try_expr(q0, q1_name, q1, q2_name, q2):
    val = q0
    if q1_name is not None:
        val = val + Fraction_to_mp(q1) * BASIS[q1_name]
    if q2_name is not None:
        val = val + Fraction_to_mp(q2) * BASIS[q2_name]
    relerr = abs(val - D1) / abs(D1)
    comp = comp_rat(q0) + (comp_rat(q1) if q1_name else 0) + (comp_rat(q2) if q2_name else 0)
    return val, relerr, comp

def Fraction_to_mp(f):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)

# Single basis term: q0 + q1*B_i
for bi in NAMES:
    for q1 in RATS:
        if q1 == 0:
            continue
        c_q1 = comp_rat(q1)
        if c_q1 >= MAXCOMP:
            continue
        remaining = D1 - Fraction_to_mp(q1) * BASIS[bi]
        q0 = best_rational_approx(remaining, MAXDEN)
        c_q0 = comp_rat(q0)
        total_comp = c_q0 + c_q1
        if total_comp > MAXCOMP:
            continue
        val, relerr, comp2 = try_expr(q0, bi, q1, None, None)
        results.append((total_comp, float(relerr), str(q0), bi, str(q1), None, None))
        single_count += 1

t1 = time.time()
print(f"Single-term search: {single_count} evaluated in {t1-t0:.2f}s")

# Two basis terms: q0 + q1*B_i + q2*B_j, i < j (distinct)
for ii in range(len(NAMES)):
    for jj in range(ii+1, len(NAMES)):
        bi, bj = NAMES[ii], NAMES[jj]
        for q1 in RATS:
            if q1 == 0:
                continue
            c_q1 = comp_rat(q1)
            if c_q1 >= MAXCOMP:
                continue
            for q2 in RATS:
                if q2 == 0:
                    continue
                c_q2 = comp_rat(q2)
                if c_q1 + c_q2 >= MAXCOMP:
                    continue
                remaining = D1 - Fraction_to_mp(q1)*BASIS[bi] - Fraction_to_mp(q2)*BASIS[bj]
                q0 = best_rational_approx(remaining, MAXDEN)
                c_q0 = comp_rat(q0)
                total_comp = c_q0 + c_q1 + c_q2
                if total_comp > MAXCOMP:
                    continue
                val, relerr, _ = try_expr(q0, bi, q1, bj, q2)
                results.append((total_comp, float(relerr), str(q0), bi, str(q1), bj, str(q2)))
                pair_count += 1

t2 = time.time()
print(f"Pair-term search: {pair_count} evaluated in {t2-t1:.2f}s")
print(f"Total results: {len(results)}")

# Build Pareto frontier: sort by complexity, keep track of min relerr seen so far
results.sort(key=lambda r: (r[0], r[1]))
frontier = []
best_err_so_far = None
for r in results:
    comp, relerr = r[0], r[1]
    if best_err_so_far is None or relerr < best_err_so_far:
        frontier.append(r)
        best_err_so_far = relerr
    if best_err_so_far is not None and best_err_so_far < 1e-12:
        pass

with open("testB_frontier.txt", "w") as f:
    f.write(f"Total single-term evaluated: {single_count}\n")
    f.write(f"Total pair-term evaluated: {pair_count}\n")
    f.write(f"Total results: {len(results)}\n")
    f.write("Pareto frontier (complexity, relerr, expression):\n")
    for r in frontier:
        comp, relerr, q0, bi, q1, bj, q2 = r
        expr = f"{q0} + {q1}*{bi}"
        if bj is not None:
            expr += f" + {q2}*{bj}"
        f.write(f"comp={comp:3d}  relerr={relerr:.3e}  {expr}\n")

print(f"Frontier size: {len(frontier)}")
for r in frontier[:20]:
    print(r)

t3 = time.time()
print(f"Total Test B time: {t3-t0:.2f}s")
