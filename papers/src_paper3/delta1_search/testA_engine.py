"""Test A engine: E = a + b*pi + c*phi, a,b,c rational, comp = |p|+q each, comp(0)=0.
Enumerate all distinct (a,b,c) with total complexity <= MAXC.
"""
import mpmath as mp
from fractions import Fraction
from rationals import rationals_by_complexity

mp.mp.dps = 60
PI = mp.pi
PHI = (1 + mp.sqrt(5)) / 2

def frac_to_mp(f: Fraction):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)

def enumerate_exprs(max_comp):
    """Yield (a,b,c,total_comp) for all combos with total_comp <= max_comp."""
    buckets = rationals_by_complexity(max_comp)
    # flat list of (frac, comp) for comp 0..max_comp
    flat = []
    for c in range(max_comp + 1):
        for r in buckets[c]:
            flat.append((r, c))
    seen = set()
    results = []
    for a, ca in flat:
        if ca > max_comp:
            continue
        rem1 = max_comp - ca
        for b, cb in flat:
            if cb > rem1:
                continue
            rem2 = rem1 - cb
            for c, cc in flat:
                if cc > rem2:
                    continue
                total = ca + cb + cc
                key = (a, b, c)
                if key in seen:
                    continue
                seen.add(key)
                results.append((a, b, c, total))
    return results

def eval_expr(a, b, c):
    return frac_to_mp(a) + frac_to_mp(b) * PI + frac_to_mp(c) * PHI

import math as _math
PI_F = _math.pi
PHI_F = (1 + _math.sqrt(5)) / 2

def eval_expr_float(a, b, c):
    return float(a) + float(b) * PI_F + float(c) * PHI_F

if __name__ == "__main__":
    import time
    t0 = time.time()
    res = enumerate_exprs(30)
    t1 = time.time()
    print("count:", len(res), "time:", t1 - t0)
