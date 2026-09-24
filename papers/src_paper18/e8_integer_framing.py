#!/usr/bin/env python3
r"""
Integer-framing step: is Tw=m-15 a genuine Calugareanu framing of the trefoil T(2,3), and is the
E_8 mirror m<->30-m literal ribbon orientation/framing reversal? (2026-09-03).
See notes/quark_generation_e8_ribbon_twist.md.

Calugareanu-White-Fuller:  SL = Wr + Tw
  SL = self-linking (integer framing = linking of the two ribbon edges)
  Wr = writhe of the core knot  (trefoil standard diagram: +3 for right-handed)
  Tw = total twist of the ribbon about the core
Generation twist (e8_twist_framing.py):  Tw(m) = m-15 (integer, EVEN since exponents are odd).
Framing phase:  T^(m) = (4m-7)/360.
"""
from fractions import Fraction as F

Wr = 3                                              # right-handed trefoil writhe (blackboard framing SL=Wr)
exps = [1,7,11,13,17,19,23,29]
assign = {'t':1,'b':7,'c':11,'s':13,'d':17,'u':19}  # exponent per quark
def Tw(m): return m-15
def SL(m): return Wr + Tw(m)                        # = m-12
def Tphase(m): return F(4*m-7,360)

print("="*72); print("INTEGER FRAMING OF THE TREFOIL: SL = Wr + Tw,  Wr=3"); print("="*72)

print("\n[1] Framing per quark  (Tw=m-15 even, SL=Wr+Tw=m-12 integer):")
print("     q   m   Tw    SL=m-12   Tw even?   SL parity")
for q,m in sorted(assign.items(), key=lambda kv: Tw(kv[1])):
    print(f"    {q:>2}  {m:>2}  {Tw(m):>+3}    {SL(m):>+3}      {'yes' if Tw(m)%2==0 else 'NO':>3}       {'odd' if SL(m)%2 else 'even'}")
print("    unused: m=23 Tw=+8 SL=11 ; m=29 Tw=+14 SL=17")
print("    => every Tw is an EVEN integer, every SL an integer (odd, since Wr=3 odd + Tw even).")
print("       So Tw=m-15 is a genuine integer ribbon framing, not a half-integer/phase artefact.")

print("\n[2] Framing-phase / twist consistency:")
for m in [13,17,15,19,1]:
    lhs = Tphase(m); rhs = F(53,360) + F(Tw(m),90)
    tag = "(centre, m=15 not an exponent)" if m==15 else ""
    print(f"    m={m:>2}: T^(m)={str(lhs):>7}   53/360 + Tw/90 = {str(rhs):>7}   {'OK' if lhs==rhs else 'X'} {tag}")
print("    => T^(m) = 53/360 + Tw/90 EXACTLY: the framing phase advances by 1/90 = 1/(k*h(E_8)) per unit")
print("       twist, about the centre 53/360. The centre (Tw=0) is m=15 -- NOT an E_8 exponent -> the")
print("       untwisted ribbon is the generation-less (leptonic) limit, no quark sits there.")

print("\n[3] MIRROR m<->30-m as a ribbon operation:")
print("    Tw(30-m) = (30-m)-15 = -(m-15) = -Tw(m)          -> TWIST REVERSAL (exact).")
print("    SL(30-m) = Wr + Tw(30-m) = 3 - Tw(m) = 6 - SL(m) -> reflection of SL about Wr=3.")
for a,b in [(1,29),(7,23),(11,19),(13,17)]:
    print(f"      ({a:>2},{b:>2}): Tw {Tw(a):>+3}/{Tw(b):>+3} (sum {Tw(a)+Tw(b)});  SL {SL(a):>+3}/{SL(b):>+3} (sum {SL(a)+SL(b)}=2Wr)")
print("    => the E_8 mirror is FRAMING CONJUGATION: it reverses the twist about the blackboard framing")
print("       (Tw=0), holding the KNOT (writhe/crossings = COLOUR) FIXED. It is a pure GENERATION")
print("       operation, NOT a full knot-mirror (which would also flip Wr, i.e. flip the crossings).")

print("\n[4] VERDICT:")
print("    - generation = ribbon TWIST is now LITERAL: Tw=m-15 is a genuine even integer framing of the")
print("      trefoil, SL=m-12 the self-linking, centred on the untwisted (leptonic) framing m=15.")
print("    - the E_8 duality = framing conjugation (twist reversal, knot fixed) = a generation symmetry")
print("      that leaves colour untouched -- exactly the orthogonality the three-structure picture needs.")
print("    - the algebra->geometry map is closed at the level 2I<->two-tube: E_8 exponent <-> even ribbon")
print("      framing; E_8 Coxeter number h=30 <-> the twist range/period; mirror <-> framing conjugation.")
print("    CAVEAT [O]: the framing-anomaly phase is 1/90 per twist; confirming this is the trefoil's")
print("    ACTUAL Chern-Simons framing anomaly (topological spin) -- not just an E_8-derived integer that")
print("    is framing-shaped -- is the remaining step to make it fully first-principles.")
print("="*72)
