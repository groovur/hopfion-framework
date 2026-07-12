"""Test C: residual series hunt for delta = D1 - C.
delta ~ q * pi**a * phi**b, a in {-2..2}, b in {-30..10}, q rational den<=12 |num|<=48
found via continued-fraction best rational approximation.
"""
from fractions import Fraction
import mpmath as mp
from common import delta, pi, phi

mp.mp.dps = 60

MAXDEN = 12
MAXNUM = 48

def comp_rat(f: Fraction):
    return 0 if f == 0 else abs(f.numerator) + f.denominator

def best_rational_cf(x, maxden, maxnum):
    """Best rational approx to mpmath value x with denominator <= maxden via continued fractions,
    also capped at |numerator|<=maxnum. Returns Fraction."""
    # mpmath has mp.pslq or nstr; use continued fraction expansion manually.
    from mpmath import mpf
    a0 = mp.floor(x)
    # continued fraction convergents
    convergents = []
    x0 = x
    h_prev2, h_prev1 = 0, 1
    k_prev2, k_prev1 = 1, 0
    val = x
    for _ in range(30):
        a = mp.floor(val)
        h = a*h_prev1 + h_prev2
        k = a*k_prev1 + k_prev2
        convergents.append((int(h), int(k)))
        if k > maxden or abs(h) > maxnum:
            break
        frac_part = val - a
        if abs(frac_part) < mp.mpf(10)**-50:
            break
        val = 1/frac_part
        h_prev2, h_prev1 = h_prev1, h
        k_prev2, k_prev1 = k_prev1, k
    # pick best convergent (and semiconvergents) within bounds
    best = None
    best_err = None
    for h, k in convergents:
        if k == 0:
            continue
        if abs(h) > maxnum or k > maxden or k < 1:
            continue
        cand = Fraction(h, k)
        err = abs(mp.mpf(cand.numerator)/cand.denominator - x)
        if best_err is None or err < best_err:
            best_err = err
            best = cand
    # also brute-force small denominators to be safe (cheap, maxden<=12)
    for k in range(1, maxden+1):
        p = mp.nint(x * k)
        if abs(p) > maxnum:
            continue
        cand = Fraction(int(p), k)
        err = abs(mp.mpf(cand.numerator)/cand.denominator - x)
        if best_err is None or err < best_err:
            best_err = err
            best = cand
    return best

results = []  # (a, b, q, relerr)
for a in range(-2, 3):
    for b in range(-30, 11):
        base = pi**a * phi**b
        r = delta / base
        q = best_rational_cf(r, MAXDEN, MAXNUM)
        if q is None or q == 0:
            continue
        approx = mp.mpf(q.numerator)/mp.mpf(q.denominator) * base
        relerr = abs(approx - delta) / abs(delta)
        results.append((a, b, q, float(relerr), comp_rat(q)))

results.sort(key=lambda r: r[3])

with open("testC_hits.txt", "w") as f:
    f.write("Top 20 hits ranked by relative error:\n")
    for a, b, q, relerr, c in results[:20]:
        f.write(f"a={a:3d} b={b:4d} q={str(q):>8s} comp(q)={c:3d}  relerr={relerr:.3e}\n")

    f.write("\nTop 10 ranked by complexity(q) among rel err < 1e-2:\n")
    sub = [r for r in results if r[3] < 1e-2]
    sub.sort(key=lambda r: (r[4], r[3]))
    for a, b, q, relerr, c in sub[:10]:
        f.write(f"a={a:3d} b={b:4d} q={str(q):>8s} comp(q)={c:3d}  relerr={relerr:.3e}\n")

    f.write("\nFlagged (relerr < 1e-4 and comp(q) <= 10):\n")
    flagged = [r for r in results if r[3] < 1e-4 and r[4] <= 10]
    if not flagged:
        f.write("NONE\n")
    else:
        for a, b, q, relerr, c in flagged:
            f.write(f"a={a:3d} b={b:4d} q={str(q):>8s} comp(q)={c:3d}  relerr={relerr:.3e}\n")

print("Top 10:")
for r in results[:10]:
    print(r)
print("\nFlagged (relerr<1e-4, comp<=10):", [r for r in results if r[3] < 1e-4 and r[4] <= 10])

# specific framework constants
specific = {
    "1/(9*phi**6)": 1/(9*phi**6),
    "3/(2*pi)": mp.mpf(3)/(2*pi),
    "1/360": mp.mpf(1)/360,
    "1/phi**15": 1/phi**15,
    "1/(4*pi*phi**6)": 1/(4*pi*phi**6),
    "(15/8-phi)/phi**12": (mp.mpf(15)/8 - phi)/phi**12,
}
with open("testC_specific.txt", "w") as f:
    f.write(f"delta = {mp.nstr(delta, 50)}\n\n")
    for name, val in specific.items():
        ratio = delta / val
        relerr = abs(val - delta)/abs(delta)
        f.write(f"{name:25s} = {mp.nstr(val,30)}  delta/this = {mp.nstr(ratio,20)}  relerr(this vs delta) = {float(relerr):.3e}\n")

print("\nSpecific constants:")
for name, val in specific.items():
    ratio = delta / val
    print(name, mp.nstr(val, 20), "ratio delta/val =", mp.nstr(ratio, 15))
