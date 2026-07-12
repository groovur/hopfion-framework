"""Enumerate reduced rationals grouped by complexity comp(p/q) = |p| + q, comp(0) = 0."""
from fractions import Fraction
from math import gcd

def rationals_by_complexity(max_comp):
    """Return dict: complexity -> list of Fraction (reduced), including negatives, and 0 at comp 0."""
    buckets = {c: [] for c in range(max_comp + 1)}
    buckets[0].append(Fraction(0, 1))
    for q in range(1, max_comp + 1):
        for p in range(-(max_comp - q), (max_comp - q) + 1):
            if p == 0:
                continue
            if gcd(abs(p), q) != 1:
                continue
            c = abs(p) + q
            if c > max_comp:
                continue
            buckets[c].append(Fraction(p, q))
    return buckets

def all_rationals_upto(max_comp):
    buckets = rationals_by_complexity(max_comp)
    out = []
    for c in range(max_comp + 1):
        out.extend((r, c) for r in buckets[c])
    return out

if __name__ == "__main__":
    b = rationals_by_complexity(10)
    for c in range(11):
        print(c, len(b[c]), b[c][:5])
