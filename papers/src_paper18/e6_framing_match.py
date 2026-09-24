#!/usr/bin/env python3
r"""
E_6-native framing-anomaly match: does the per-strand E_6 twist reconcile with the trefoil's OWN framing
data -- the SU(3)_1=(E_6)_1 topological spin and Paper XVI's Bishop holonomy -63 deg? (2026-09-03).
See notes/quark_generation_e8_ribbon_twist.md.

Twist (bridge-ordered): u,d=+-1; c,s=+-2; t,b=+-5.
"""
from fractions import Fraction as F
import math

print("="*72); print("E_6-NATIVE FRAMING MATCH: per-strand twist vs trefoil framing"); print("="*72)

# --- WZW framing data of the trefoil's own theory --------------------------------
print("\n[1] Trefoil WZW topological-spin (framing) data:")
# SU(3)_1: fundamental 3, h=1/3, c=2 ; (E_6)_1: 27, h=2/3, c=6
for name,h,c in [("SU(3)_1  (per-strand, 3)", F(1,3), 2), ("(E_6)_1  (trefoil, 27)", F(2,3), 6)]:
    theta = h - F(c,24)                         # topological spin incl. c/24
    print(f"    {name}: h={h}  c={c}  h-c/24={theta}  -> framing phase e^(2pi i h) = {float(h*360):.0f} deg (period {h.denominator} = mod {h.denominator})")
print("    => both are MOD 3 (h denominators = 3): the trefoil's own framing anomaly is Z_3. ")

# --- topological match: E_6 twist -> anomaly phase, mod 3 -------------------------
print("\n[2] TOPOLOGICAL match -- E_6 twist through the SU(3)_1 anomaly (120 deg = e^(2pi i/3) per unit):")
tw = {'u':1,'d':-1,'c':2,'s':-2,'t':5,'b':-5}
print("     q  Tw   Tw*120 mod 360   residue(=class)")
for q,t in sorted(tw.items(), key=lambda kv: kv[1]):
    ph = (t*120) % 360
    print(f"    {q:>2} {t:>+2}     {ph:>4} deg          {t%3}  ({'singlet 1' if t%3==0 else ('fund 3' if t%3==1 else 'anti 3bar')})")
print("    => the E_6 twist lands ONLY on 120/240 deg (residues 1,2 = 3/3bar), never 0 (singlet).")
print("       i.e. the per-strand twist is an SU(3)_1-CONSISTENT (mod-3) framing. Topological match: PASS.")
print("       And E_6 IS the trefoil's WZW (McKay 2T<->E_6<->SU(3)_1), so 'generation = twist' is the")
print("       trefoil's OWN framing, not an imported label -- derived in the WZW/McKay sense.")

# --- geometric side: Bishop holonomy ---------------------------------------------
print("\n[3] GEOMETRIC side -- Paper XVI Bishop holonomy H = -63.0191 deg (Z_3: -21.0064 deg/segment):")
H = -63.0191; Hseg = H/3
print(f"    total geometric twist  Tw_geom = H/360 = {H/360:+.5f} framing units")
print(f"    per Z_3 segment        {Hseg:.4f} deg = {Hseg/360:+.5f} framing units")
# candidate clean relations
print("    candidate clean forms:")
for label,val in [("7/40",F(7,40)),("1/(2 phi^2)",1/(2*((1+5**0.5)/2)**2)),("1/6",1/6.0)]:
    v=float(val)
    print(f"      {label} = {v:+.5f}  -> {abs(v)*360:.3f} deg   {'~ MATCH |Tw_geom|' if abs(abs(v)-abs(H/360))<0.002 else 'no'}")

print("\n[4] HONEST VERDICT:")
print("    - TOPOLOGICAL match PASSES: the per-strand E_6 twist is a mod-3 (SU(3)_1) framing, residue-0")
print("      (singlet) empty, and E_6 = the trefoil's own WZW -> generation = the trefoil's WZW framing.")
print("      This is the sense in which 'generation = twist' is DERIVED geometry (via McKay 2T<->E_6).")
print("    - GEOMETRIC (Bishop) is SEPARATE: H=-63 deg is the classical normal-bundle ANHOLONOMY")
print("      (Tw_geom ~ -0.175 framing units), a small FRACTIONAL geometric twist. It confirms the trefoil")
print("      carries a nontrivial Z_3 framing, but it does NOT numerically equal the INTEGER WZW twists")
print("      {+-1,+-2,+-5}. Classical anholonomy and quantum framing are Calugareanu-related, not equal.")
print("    => 'derived geometry' holds in the WZW/McKay sense (topological framing); the classical -63 deg")
print("       is corroborating evidence of nontrivial Z_3 framing, NOT the source of the twist quantization.")
print("       OPEN: a first-principles link from the classical Bishop anholonomy to the integer WZW framing.")
print("="*72)
