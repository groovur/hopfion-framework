#!/usr/bin/env python3.11
"""
Step 3, first calc (2026-09-22): does a sigma*R + barrier FORMATION-energy correction reproduce the
non-perturbative quark-mass residuals that survive RG running (Paper XVI colour-shift route,
P16:tab:rg_running)? Tests the sec.43 hypothesis. Honest: check the SCALE first (is the residual the
formation/constituent scale ~sigma*R, or something else?).
"""
import numpy as np
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; R=R0+r0     # confinement scale (sec.19)
# Paper XVI colour-shift route, run to comparison scale (P16:tab:rg_running)
q=['u','d','s','c','b','t']
m_static={'u':1.300,'d':0.497,'s':93.32,'c':265.7,'b':1930.,'t':3554.}   # MeV, post-RG
m_pdg   ={'u':2.16,'d':4.67,'s':93.40,'c':1275.,'b':4180.,'t':172570.}
# framework confinement/formation scales
sigma=646.5/9        # MeV/unit length (sec.19: sigma_tot = E_baryon/(3R0), fixed)
m_q0=215.5           # constituent scale
print("="*70,"\nPost-RG colour-shift residuals (Paper XVI) & the 'missing' energy\n","="*70,sep="")
print(f"  {'q':>2} {'m_static':>9} {'m_PDG':>9} {'ratio':>6} {'Dm=PDG-stat':>11} {'factor':>7}")
for k in q:
    dm=m_pdg[k]-m_static[k]; fac=m_pdg[k]/m_static[k]
    print(f"  {k:>2} {m_static[k]:9.2f} {m_pdg[k]:9.2f} {m_static[k]/m_pdg[k]:6.3f} {dm:11.2f} {fac:7.2f}")

print("\n"+"="*70,"\nis the missing energy the sigma*R FORMATION scale?\n","="*70,sep="")
print(f"  framework formation/confinement scales: sigma*R = {sigma*R:.0f} MeV, sigma*3R0 = {sigma*3*R0:.0f}"
      f" MeV, constituent m_q0 = {m_q0} MeV")
dms=[m_pdg[k]-m_static[k] for k in q]
print(f"  'missing' energy Dm spans {min(dms):.2f} to {max(dms):.0f} MeV -- SIX orders of magnitude.")
print(f"  sigma*R ~ {sigma*R:.0f} MeV is a ~CONSTANT (a fixed formation/confinement scale).")
print(f"  => a constant sigma*R canNOT match Dm spanning 0.08 - 169000 MeV. The current-mass residual")
print(f"     is NOT an additive formation energy.")

print("\n"+"="*70,"\nSCALE DIAGNOSIS: two DIFFERENT non-perturbative objects\n","="*70,sep="")
print("  (1) FORMATION/CONFINEMENT energy (sigma*R ~ 278 MeV): sets the CONSTITUENT / HADRON scale.")
print("      sec.19: proton = sigma*(3R0+R) = 924.8 MeV = 98.6% m_p. THIS is where sigma*R lives -- in")
print("      HADRON masses, largely DONE. It is NOT a correction to CURRENT quark masses.")
print("  (2) CURRENT-quark-mass residuals (this calc): the static/colour-shift CURRENT masses undershoot")
print("      by MULTIPLICATIVE factors (u1.7 d9.4 s1.0 c4.8 b2.2 t49) -- route-dependent, non-monotonic,")
print("      NOT a fixed sigma*R. A DISTINCT non-perturbative effect on the current-mass tower.")
print("  => sec.43 CONFLATED these. sigma*R (formation) is the HADRON scale (done, sec.19); the CURRENT-")
print("     mass residuals (sec.42's ~20-30%) are a SEPARATE non-perturbative correction, NOT the")
print("     formation energy. FIRST-CALC RESULT: formation-energy(sigma*R) hypothesis for the current-")
print("     mass residuals is FALSE on scale grounds.")
print("\n  NOTE: the two static routes DISAGREE on the residuals -- E_8 route (sec.24): down ~1.2 OVER,")
print("  up under; colour-shift (here): strange exact, ALL others under. So 'the residual' is not even")
print("  route-invariant -> the current-mass non-perturbative correction is not yet a well-posed single")
print("  quantity. Step 3 must FIRST fix the static route, THEN define the residual it corrects.")
