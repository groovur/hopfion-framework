#!/usr/bin/env python3
r"""
Generation = ribbon TWIST: map each E_8 exponent to a framing/twist and test mirror = twist reversal.
(2026-09-03). See notes/quark_generation_e8_ribbon_twist.md.

Chain (standard CFT + Calugareanu): the WZW T-matrix is the framing anomaly (Dehn twist / topological
spin). The E_8 exponent m enters the quark mass via the framing phase T^(m)=(4m-7)/360. By
Calugareanu SL=Wr+Tw, with the writhe Wr fixed by the trefoil (colour), the framing phase IS the twist Tw.
Test: is the E_8 mirror duality m <-> h-m (h=30) exactly twist reversal Tw <-> -Tw?
"""
from fractions import Fraction as F

h = 30
exps = [1,7,11,13,17,19,23,29]
def T(m): return F(4*m-7, 360)                      # framing phase of exponent m

print("="*72); print("E_8 EXPONENT -> RIBBON TWIST, and the mirror-reversal test  (h=30)"); print("="*72)

# mirror-pair sum is constant -> the symmetry centre of the framing phase
pairs = [(1,29),(7,23),(11,19),(13,17)]
sums = {(a,b): T(a)+T(b) for a,b in pairs}
c0 = list(sums.values())[0]/2                        # centre = half the (constant) mirror-pair sum
print("\n[A] Mirror-pair framing-phase sums T^(m)+T^(30-m):")
for (a,b),s in sums.items():
    print(f"    ({a:>2},{b:>2}): {str(s):>7}")
allconst = len(set(sums.values()))==1
print(f"    all equal? {allconst}  ->  centre c0 = (sum)/2 = {str(c0)} = 53/360")

# twist = framing phase relative to the centre
def dT(m): return T(m) - c0
print("\n[B] Twist  Tw(m) = T^(m) - c0 :")
print("     m    T^(m)      Tw = T^(m)-c0    Tw*90    check (m-15)/90")
for m in exps:
    val = dT(m)
    print(f"    {m:>2}   {str(T(m)):>7}   {str(val):>9}     {str(val*90):>4}       {str(F(m-15,90)):>7}  "
          f"{'OK' if val==F(m-15,90) else 'X'}")
print("    => Tw(m) = (m-15)/90  EXACTLY: the twist is LINEAR in the exponent, centred at m=15=h/2.")

print("\n[C] MIRROR = TWIST REVERSAL test  (Tw(30-m) == -Tw(m)) :")
ok = all(dT(30-m) == -dT(m) for m in exps)
for a,b in pairs:
    print(f"    Tw({a:>2})={str(dT(a)):>7}   Tw({b:>2})={str(dT(b)):>7}   sum={str(dT(a)+dT(b)):>3}  "
          f"{'reversal OK' if dT(a)==-dT(b) else 'FAIL'}")
print(f"    => mirror m<->30-m is EXACT twist reversal Tw<->-Tw for all pairs: {ok}")

print("\n[D] Twist per quark (Tw in units of 1/90; equivalently m-15):")
assign = {'u':(1,19),'c':(2,11),'t':(3,1),'d':(1,17),'s':(2,13),'b':(3,7)}
pdg = {'u':2.16,'d':4.67,'s':93.4,'c':1275,'b':4180,'t':172570}
print("     q  gen  m   Tw=(m-15)   |Tw|    isospin   PDG mass(MeV)")
for q,(g,m) in sorted(assign.items(), key=lambda kv: -kv[1][1]):
    iso = 'up (+2/3)' if q in ('u','c','t') else 'down(-1/3)'
    print(f"    {q:>2}   {g}  {m:>2}    {m-15:>+3}       {abs(m-15):>2}    {iso}   {pdg[q]:>8.1f}")
print("    unused exponents 23,29 -> Tw=+8,+14 (the POSITIVE mirrors of b(-8), t(-14)).")

print("\n[E] READ:")
print("    - The E_8 exponent IS a twist: m = 15 + (twist), centred on the untwisted point m=15.")
print("    - The mirror duality m<->30-m is EXACTLY twist reversal (orientation/chirality flip).")
print("    - Strange (Tw=-2) is the LEAST-twisted quark -> nearest the untwisted centre -> the anchor,")
print("      consistent with strange = most lepton-like (a lepton is an unknotted, untwisted sector).")
print("    - |Tw| tracks GENERATION: gen-1,2 small (|Tw|=2,4), gen-3 large (|Tw|=8,14). Isospin (up/down)")
print("      is the orthogonal Z_2, NOT |Tw| -- exactly the three-independent-structures picture.")
print("    - Only the NEGATIVE-twist extremes (t=-14,b=-8) are used; the positive mirrors (29,23) are")
print("      unused -> a handedness/orientation asymmetry at the third generation (open: why one sign).")
print("="*72)
