#!/usr/bin/env python3
r"""
Bijection chase (2026-09-03): do the 6 trefoil midpoint-breaking patterns (Paper XVIII,
P18:conj:six_quarks) map to the 6 E_8 Coxeter exponents / 6 quark flavours?

Setup:
  - Trefoil T(2,3) has 3 crossing vertices = 3 COLOURS (Z_3, SU(3)_1; Papers XV/XVI), at z=+r_0.
  - 3 midpoint torsion-twist sites at z=0 are the breaking sites.
  - Six breaking patterns: C(3,1)=3 single + C(3,2)=3 double  (breaking all 3 -> vacuum; 0 -> baryon).
  - E_8 route (Paper V/XVII): 6 exponents {1,7,11,13,17,19}; up {1,11,19}=(t,c,u), down {7,13,17}=(b,s,d).
    Within an isospin type the 3 exponents are the 3 GENERATIONS (Z_3-BROKEN: masses span ~5 dex).

The test: can the Z_3-symmetric breaking triplets carry the Z_3-broken generation spread?
"""
from itertools import combinations

sites = (0, 1, 2)                                   # 3 midpoints = the trefoil's Z_3 (colour) triad
single = list(combinations(sites, 1))               # C(3,1) = 3
double = list(combinations(sites, 2))               # C(3,2) = 3

def z3_orbit(patterns):
    """cyclic Z_3 rotation r: i -> (i+1) mod 3; return orbit sizes."""
    def rot(p): return tuple(sorted((i+1) % 3 for i in p))
    seen, orbits = set(), []
    for p in patterns:
        if p in seen: continue
        orb = {p}; q = rot(p)
        while q not in orb: orb.add(q); q = rot(q)
        seen |= orb; orbits.append(orb)
    return orbits

print("="*74); print("SIX BREAKING PATTERNS vs SIX E_8 EXPONENTS: the bijection chase"); print("="*74)

print("\n[A] Breaking patterns and their Z_3 (colour) structure:")
print(f"    single-breakings C(3,1) = {single}   Z_3 orbits: {[len(o) for o in z3_orbit(single)]}")
print(f"    double-breakings C(3,2) = {double}   Z_3 orbits: {[len(o) for o in z3_orbit(double)]}")
print("    => each triplet is a SINGLE Z_3 orbit of size 3: the 3 members are Z_3-EQUIVALENT")
print("       (related by the colour rotation), hence degenerate under the unbroken Z_3.")
print("    => the 6 breakings factor as  2 (single/double = DEPTH) x 3 (Z_3 site = COLOUR).")

print("\n[B] The six quark flavours factor differently:")
up   = {'t':1,'c':11,'u':19}; down = {'b':7,'s':13,'d':17}
print(f"    up-type   E_8 exponents (=GENERATION): {up}   (t,c,u)")
print(f"    down-type E_8 exponents (=GENERATION): {down} (b,s,d)")
print("    => 6 flavours factor as  2 (isospin) x 3 (GENERATION).")
print("       The generation triplet is Z_3-BROKEN: exponents 1,11,19 (up) give masses spanning")
print("       ~5 orders of magnitude -- the opposite of the degenerate breaking triplet.")

print("\n[C] The obstruction (why the naive bijection fails):")
print("    A bijection breaking-pattern <-> exponent must map a Z_3-SYMMETRIC triplet (the 3")
print("    Z_3-equivalent single-breakings) onto a Z_3-BROKEN triplet (the 3 generations).")
print("    Nothing in the breaking geometry breaks Z_3 -- the sites are colour-related, and colour")
print("    is UNBROKEN (confinement). So the '3' in the breakings (COLOUR) is not the '3' in the")
print("    flavours (GENERATION). They coincide in COUNT (3=3) but are DIFFERENT quantum numbers.")

print("\n[D] Resolution / the idea this opens:")
print("    - 6 breaking patterns  = 2 isospin x 3 COLOUR   (geometry: writhe/depth x crossing site)")
print("    - 6 quark flavours     = 2 isospin x 3 GENERATION (algebra: writhe x E_8 exponent)")
print("    P18:conj:six_quarks identifies these because both are 6, but it conflates COLOUR with")
print("    GENERATION. Generation is NOT geometric (no Z_3-breaking in the trefoil); it is the")
print("    ALGEBRAIC E_8 label -- consistent with the no-saddle result: the generation/mass spread")
print("    is algebraic, not a configuration energy.  Three orthogonal quantum numbers, three")
print("    structures:  ISOSPIN (writhe, Paper XVIII) x COLOUR (crossings/SU(3)_1, Papers XV/XVI)")
print("                 x GENERATION (E_8 exponents, Papers V/XVII).")

print("\n[E] Secondary tension flagged for Paper XVIII:")
print("    P18 carries TWO isospin mechanisms -- single/double DEPTH (conj:six_quarks) and writhe")
print("    z-LEVEL (conj:isospin) -- but all breakings sit at z=0 midpoints, so 'depth' and 'z-level'")
print("    cannot both give up/down naively. Reconciling them is a further open point.")
print("="*74)
