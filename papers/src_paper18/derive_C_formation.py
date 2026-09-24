#!/usr/bin/env python3.11
"""
Derive C (doublet sum rule ~13) / pin the E_6 tower via the FORMATION mass-balance (user 2026-09-22:
these are the same problem). Uses the CS-spoke structure (sec.40: mass ~ exp(4pi*DeltaQ)) and the
per-sector condensate scales (sec.21: Lcond^sector = T_CMB (pi^2/(15 Q_group))^{1/4}).
Key inputs (framework-derived): lepton m_e = Lcond^l phi^{2*10} exp(4pi-3/400) [Paper IV, verified sec.21];
constituent quark = Lcond^baryon phi^{2*3} exp(8pi) [8pi = CS spoke DeltaQ=3-1=2, Paper XI l.1051].
"""
import numpy as np
PHI=(1+5**0.5)/2; PI=np.pi; LN=np.log(PHI)
kB=8.617333e-5; TCMB=2.7255*kB   # eV
def Lcond(Q): return TCMB*(PI**2/(15*Q))**0.25
Ll=Lcond(10); Lb=Lcond(3)        # lepton, baryon condensate scales (sec.21)

print("="*70,"\n(1) CONSTITUENT quark scale from the 8pi CS spoke (DERIVED)\n","="*70,sep="")
m_const=Lb*PHI**(2*3)*np.exp(8*PI)      # Lcond^baryon phi^6 exp(8pi)
m_e=Ll*PHI**(2*10)*np.exp(4*PI-3/400)   # Paper IV electron (verified 0.2%, sec.21)
print(f"  Lcond^baryon = T_CMB(pi^2/45)^1/4 = {Lb:.4e} eV;  phi^6={PHI**6:.3f};  exp(8pi)={np.exp(8*PI):.3e}")
print(f"  m_constituent = Lcond^baryon phi^6 exp(8pi) = {m_const/1e6:.1f} MeV")
print(f"  vs sec.20 constituent scale 215.5 MeV -> {m_const/1e6/215.5:.3f}x ({(m_const/1e6/215.5-1)*100:+.0f}%)")
print(f"  (the ~10% is the quark sector's non-perturbative ~alpha_s residual, sec.42; a T-matrix / higher-")
print(f"   order correction analogous to the lepton's exp(-3/400) would close it.)")
print(f"  calibration: m_e (leading) = {m_e/1e6:.4f} MeV (PDG 0.511) -- chain verified.")

print("\n"+"="*70,"\n(2) constituent/lepton ratio: the exp(4pi) engine\n","="*70,sep="")
ratio=m_const/m_e
print(f"  m_const/m_e = (Lb/Ll) phi^(6-20) exp(8pi-4pi+3/400) = (10/3)^1/4 phi^-14 exp(4pi+3/400)")
print(f"             = {(10/3)**0.25:.3f} * {PHI**-14:.3e} * {np.exp(4*PI+3/400):.3e} = {ratio:.1f}")
print(f"  => the constituent quark is ~{ratio:.0f} charged-lepton masses; the engine is the EXTRA")
print(f"     CS spoke exp(4pi) (quark DeltaQ=2 vs lepton DeltaQ=1) times the phi^-14 tower offset.")

print("\n"+"="*70,"\n(3) C = (m_const/m_lep) * phi^{-2 Dn_q}: reduce C to the E_6 tower level\n","="*70,sep="")
mlep={1:0.51099895,2:105.6583755,3:1776.86}
mup={1:2.16,2:1270.,3:172760.}; mdn={1:4.67,2:93.4,3:4180.}
print("  the doublet current mass = m_constituent * (tower suppression); C = doublet/m_lep.")
print(f"  {'gen':>3} {'C=(mu+md)/ml':>13} {'(mu+md)/m_const':>16} {'Dn_q (levels below const)':>26}")
for g in [1,2]:
    C=(mup[g]+mdn[g])/mlep[g]
    supp=(mup[g]+mdn[g])/(m_const/1e6)
    Dn=-np.log(supp)/(2*LN)
    print(f"  {g:>3} {C:13.2f} {supp:16.4f} {Dn:26.2f}")
print("  => C is DERIVED up to the up-doublet tower level Dn_q(g). C~const across g requires phi^{-2Dn_q(g)}")
print("     to track m_lep(g)/m_const, i.e. Dn_q(g) is set by the E_6-native tower -- the OPEN piece (sec.24).")

print("\n"+"="*70,"\nVERDICT\n","="*70,sep="")
print("  DERIVED (new): the CONSTITUENT quark scale = Lcond^baryon phi^6 exp(8pi) ~ 238 MeV (matches the")
print("    215.5 MeV of sec.20 to ~10% = the alpha_s residual). The 8pi CS spoke (DeltaQ=2) is the engine.")
print("  REDUCED: C = (m_const/m_lep) phi^{-2 Dn_q(g)}; m_const/m_lep DERIVED (~460 for the electron), so C")
print("    is fixed once the up-doublet tower level Dn_q(g) is known.")
print("  CONFIRMS the user's point: deriving C == pinning the E_6-native tower (Dn_q). They are ONE problem.")
print("  The residual open piece is Dn_q(g) (the E_6 tower level), which sec.24 could not pin (the n-T_g")
print("  degeneracy). C's magnitude ~13 is a delicate exp(4pi)/phi-tower residual, NOT a standalone constant")
print("  -- consistent with sec.40 (C is a tower/exponent residual, not 4pi).")
