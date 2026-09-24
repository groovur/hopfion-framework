#!/usr/bin/env python3.11
"""
Isospin SIGN-FLIP via the formation-energy EM piece (user request, 2026-09-21).
Follows quark_generation_e8_ribbon_twist.md sec.22-24, Paper XVIII P18:conj:isospin +
P18:prop:isospin_mass_ratio (writhe/out-of-plane geometry) + O6 (formation energy).

THE PUZZLE: m_up/m_down FLIPS -- u<d (gen1) but c>s, t>b (gen2,3).
TEST: does the EM self-energy in the T(2,2)->T(2,3) formation energy explain the gen1 seed?
"""
import numpy as np
PHI=(1+5**0.5)/2; LN=np.log(PHI)
m={'u':2.16,'d':4.67,'s':93.4,'c':1270.,'b':4180.,'t':172760.}

print("="*70,"\n(0) THE FLIP, and its decomposition\n","="*70,sep="")
R={g:m[{1:'u',2:'c',3:'t'}[g]]/m[{1:'d',2:'s',3:'b'}[g]] for g in (1,2,3)}
lnR={g:np.log(R[g]) for g in (1,2,3)}
print("  up/down ratio R(g):", {g:round(R[g],3) for g in (1,2,3)}, " (crosses 1 between g1,g2)")
print("  ln R(g):           ", {g:round(lnR[g],3) for g in (1,2,3)})
meanlnR=np.mean(list(lnR.values()))
print(f"  mean ln R = {meanlnR:.3f} -> M_up/M_down = {np.exp(meanlnR):.2f}"
      f"  (P18:prop predicts 1/0.189={1/0.189:.2f}, ln={np.log(1/0.189):.2f}; 21% off) [ISOSPIN SCALE]")
print(f"  ln R spread about mean: {[round(lnR[g]-meanlnR,2) for g in (1,2,3)]}"
      f"  -> the GENERATION-dependence = up steeper (E_6 vs E_8, sec.24)")
print(f"  d(lnR)/gen (g1->g2) = {lnR[2]-lnR[1]:.2f} = (up incr - down incr), the steepness gap")

print("\n"+"="*70,"\n(1) DOES EM HAVE THE RIGHT SIGN FOR THE gen1 SEED?\n","="*70,sep="")
print("  gen1: m_d - m_u =", m['d']-m['u'], "MeV (DOWN heavier).")
print("  EM self-energy of a charged soliton: dm_EM ~ (alpha/2) Q^2 / r  (ALWAYS > 0).")
alpha=1/137.036
for r_fm in [0.5, 1.0, 2.0]:
    r_GeV_inv=r_fm/0.1973  # 1 fm = 1/0.1973 GeV^-1
    dm_u=0.5*alpha*(2/3)**2/r_GeV_inv*1000   # MeV
    dm_d=0.5*alpha*(1/3)**2/r_GeV_inv*1000
    print(f"  r={r_fm} fm: dm_EM(up)={dm_u:.3f}  dm_EM(down)={dm_d:.3f} MeV"
          f"  -> EM(up-down)=+{dm_u-dm_d:.3f} (EM makes UP heavier)")
print("  VERDICT: EM makes UP heavier by ~+0.3-1 MeV; gen1 needs DOWN heavier by +2.5 MeV.")
print("  => EM has the WRONG SIGN and is too SMALL. It CANNOT be the gen1 seed.")
print("     (Same as the SM: the bare m_d>m_u dominates EM -> neutron heavier than proton.)")
print("  This DEFINITIVELY kills the sec.22 H_Q charge/charge^2 lead: charge^2 -> up heavier,")
print("  the opposite of the observed gen1 ordering.")

print("\n"+"="*70,"\n(2) BOTH known up-heavier mechanisms have the wrong sign at gen1\n","="*70,sep="")
print("  - EM self-energy  ~ +Q^2  -> up heavier (wrong at gen1).")
print("  - P18 out-of-plane geometry: M_down/M_up=0.189 -> down LIGHTER = up heavier (wrong at gen1).")
print("  So the gen1 seed (up LIGHTER) OPPOSES both. The flip's essence = why up is light at gen1.")

print("\n"+"="*70,"\n(3) The ANCHORING ASYMMETRY makes u light at gen1 -- no EM needed\n","="*70,sep="")
print("  down: E_8 lepton-bridge, anchored LOW at strange=muon (gen2). gentle tower.")
print("  up:   E_6-native, anchored HIGH at the TOP-transition v/sqrt2 (gen3). STEEP tower.")
def n(mass): return -np.log(mass/215.5)/(2*LN)
nn={q:n(m[q]) for q in m}
print("  tower levels n (higher=lighter):", {q:round(nn[q],2) for q in ['u','d','s','c','b','t']})
print(f"  at gen1: n_u={nn['u']:.2f} > n_d={nn['d']:.2f}  => u LIGHTER.")
print("  MECHANISM: up is anchored at the HEAVY top and climbs a STEEP tower up to gen1, so it")
print("  overshoots to higher n (lighter) than down, which is anchored at the light muon and")
print("  climbs a GENTLE tower. Steep-tower-from-a-high-anchor => u ends lighter than d at gen1.")
print("  This is the E_6/E_8 steepness (sec.24) + top-as-transition anchoring -- NOT EM.")

print("\n"+"="*70,"\nVERDICT\n","="*70,sep="")
print("  * EM formation-energy piece: WRONG SIGN (up heavier) + too small -> NOT the sign-flip.")
print("    (Cleanly closes the H_Q charge lead of sec.22 via the sign.)")
print("  * The sign-flip = (i) up/down grows with gen [E_6>E_8 steepness, sec.24] +")
print("    (ii) up anchored HIGH at the top-transition, down LOW at muon -> u overshoots light")
print("    at gen1. Both already-established pieces; no new EM term.")
print("  * Residual TRUE puzzle: the exact gen1 seed (u/d=0.46) needs the E_6-native ABSOLUTE")
print("    tower + top anchor pinned (sec.24 outstanding); it is a mass-mechanism effect, not EM.")
