#!/usr/bin/env python3.11
"""
C in (m_up+m_down)=C*m_lep: 4pi-family vs group (13) family, weighted by REAL mass uncertainties.
User (2026-09-22): 4pi is natural/phi-related/framework-native (electron mass e^{4pi}); the error
pattern differs between candidates. Key point: g1 (u,d) is ~10% uncertain -> nearly unconstraining;
g2 (c,s) is ~1.5% -> the RELIABLE anchor. Let the reliable point decide. Flag post-hoc fits (CLAUDE.md 3).
"""
import numpy as np
PI=np.pi; phi=(1+5**0.5)/2
me,mmu=0.51099895,105.6583755
# PDG with 1-sigma (light quarks dominate g1 error)
mu1,dmu1=2.16,0.4; md1,dmd1=4.67,0.4          # u,d  (~10% each)
mc,dmc=1270.,20.; ms,dms=93.4,1.5             # c,s  (~1.5%)

C1=(mu1+md1)/me; dC1=np.sqrt(dmu1**2+dmd1**2)/me
C2=(mc+ms)/mmu;  dC2=np.sqrt(dmc**2+dms**2)/mmu
print("="*64,"\nC per doublet WITH uncertainties\n","="*64,sep="")
print(f"  g1 (u,d/e):  C = {C1:.2f} +- {dC1:.2f}   -> range [{C1-dC1:.2f}, {C1+dC1:.2f}]  (light quarks: LOOSE)")
print(f"  g2 (c,s/mu): C = {C2:.3f} +- {dC2:.3f}  -> range [{C2-dC2:.2f}, {C2+dC2:.2f}]  (RELIABLE anchor)")

print("\n"+"="*64,"\ncandidates vs the RELIABLE g2 (C=%.3f +- %.3f)\n"%(C2,dC2),"="*64,sep="")
cands=[("4*pi",4*PI),("4*pi + 1/3  (=4pi + h_SU(3)_1)",4*PI+1/3),
       ("13  (Q_b+Q_l = 3+10)",13.0),("4*pi + 1/phi",4*PI+1/phi),
       ("mean(C1,C2)",(C1+C2)/2)]
for lbl,v in cands:
    sig=(v-C2)/dC2
    print(f"  {lbl:34s} = {v:8.4f}   ({sig:+.2f} sigma from g2)")
print(f"  residual C2 - 4pi = {C2-4*PI:.4f}  (1/3 = {1/3:.4f}; 1/phi = {1/phi:.4f})")

print("\n"+"="*64,"\nREAD (honest)\n","="*64,sep="")
print(f"  - g1 is ~unconstraining (±{dC1:.1f}); g2 (±{dC2:.2f}) is the anchor and gives C = 12.90.")
print(f"  - pure 4pi=12.566 sits {(4*PI-C2)/dC2:+.1f} sigma from g2 -> MILDLY DISFAVOURED (a bit low).")
print(f"  - g2 is centered almost exactly on 4pi + 1/3 = {4*PI+1/3:.3f} (0.03%), and 1/3 = h_SU(3)_1,")
print(f"    the quark's conformal weight -> C = 4pi (solid-angle, as in m_e ~ e^{{4pi}}) + quark spin 1/3.")
print(f"    FRAMEWORK-NATIVE and fits the reliable point -- but POST-HOC (residual found by eye) and")
print(f"    a SINGLE reliable point; 13 also sits within g2's +1sigma edge. NOT decidable from 2 points.")
print(f"  - So the reliable data leans to C ~ 12.9 (4pi+1/3 or 13-low), consistent with the user's 4pi")
print(f"    pull once corrected by the quark spin. The arbiter is still the O6 formation mass-balance.")
