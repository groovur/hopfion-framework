#!/usr/bin/env python3
r"""
Hypothesis (user, 2026-09-03): E_8 is the GLOBAL trefoil approximation of the PER-STRAND E_6 twist.
Per-strand = E_6 (the quark's own group, 2T<->E_6<->SU(3)_1); E_8 (2I, lepton) = the coarse global/
lepton-facing description. Test whether E_6 gives the cleaner 6-quark twist structure.

E_6: Coxeter number h=12, exponents {1,4,5,7,8,11} (exactly SIX -> six quarks, none spare).
E_8: h=30, exponents {1,7,11,13,17,19,23,29} (eight; only 6 used, 23,29 unused).
"""
from fractions import Fraction as F

E6 = {'h':12, 'exps':[1,4,5,7,8,11]}
E8 = {'h':30, 'exps':[1,7,11,13,17,19,23,29]}

def pairs(h, exps):  return [(m, h-m) for m in exps if m < h-m]
def twists(h, exps): return {m: m - h//2 for m in exps}   # Tw = m - h/2

print("="*70); print("PER-STRAND E_6 vs GLOBAL E_8 for the quark generation twist"); print("="*70)

for name,G in (('E_6',E6),('E_8',E8)):
    h=G['h']; exps=G['exps']; Tw=twists(h,exps); pr=pairs(h,exps)
    used = 6
    print(f"\n[{name}]  h={h}  exps={exps}  ({len(exps)} exponents)")
    print(f"    mirror pairs (m+m'={h}): {pr}   ({len(pr)} pairs)")
    print(f"    twists Tw=m-{h//2}: {[Tw[m] for m in exps]}")
    print(f"    parity of Tw: {['even' if Tw[m]%2==0 else 'odd' for m in exps]}")
    spare = len(exps)-used
    print(f"    exponents for 6 quarks: {'ALL '+str(len(exps))+' used, NONE spare -> clean' if len(exps)==6 else str(used)+' used, '+str(spare)+' SPARE (23,29) -> the handedness puzzle'}")

print("\n[COMPARISON]")
print("    E_6: 6 exponents = 6 quarks, 3 COMPLETE mirror pairs, twists {+-1,+-2,+-5}. No spare.")
print("         -> DISSOLVES the E_8 unused-exponent (23,29) / gen-3 handedness puzzle: it was an")
print("            artefact of using the 8-exponent GLOBAL group for a 6-strand object.")
print("    E_8: 8 exponents, 2 spare (23,29). The 'handedness' was those 2 spares = E_8 having more")
print("         d.o.f. than the 6 strands need. Exactly what a global OVER-description would show.")

print("\n[WHY E_8 STILL APPEARS -- the architecture makes it the lepton bridge, not the fundamental group]")
print("    Paper XVI: m_q = m_ell[2I/E_8] * colour-correction[2T/E_6/SU(3)_1]. The mass REFERENCES the")
print("    lepton (2I/E_8), so E_8 enters as the lepton-facing BRIDGE; the strange<->muon anchor is an")
print("    E_8 bridge coincidence (a quark T-phase hitting a LEPTON T-phase). The INTRINSIC per-strand")
print("    twist is E_6 (the trefoil's own SU(3)_1 framing anomaly). => 'E_8 = global approximation of")
print("    per-strand E_6' is consistent with the architecture AND resolves the spare-exponent puzzle.")

print("\n[E_6 twist structure -> 6 quarks, cleanly]")
h=12; Tw=twists(h,E6['exps'])
print("    3 mirror pairs = 3 twist magnitudes {1,2,5} x 2 signs:")
for a,b in pairs(h,E6['exps']):
    print(f"      ({a:>2},{b:>2})  Tw = {Tw[a]:>+2} / {Tw[b]:>+2}   |Tw|={abs(Tw[a])}")
print("    => |Tw| in {1,2,5} = 3 generations; sign = isospin (up/down). All 6 used, mirror-complete.")
print("    (h/2=6 is the untwisted/leptonic centre; not an E_6 exponent -- same clean 'lepton=untwisted'.)")

print("\n[HONEST -- what must still be checked]")
print("    - Reproduce the STRANGE ANCHOR in E_6: is there an E_6 T-phase = a lepton T-phase? (the anchor")
print("      is cross-group E_6<->E_8, so it may live only in the E_8 bridge -- consistent, but verify.)")
print("    - The MASS values from E_6 exponents (the E_8 route gave down-sector ratios right; redo in E_6).")
print("    - The FRAMING-ANOMALY MATCH becomes E_6-native: trefoil SU(3)_1=(E_6)_1, c=6, Bishop holonomy")
print("      -63deg -- match the per-strand E_6 twist to THAT, not to E_8. This is the clean target.")
print("="*70)
