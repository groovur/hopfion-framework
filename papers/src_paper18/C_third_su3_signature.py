#!/usr/bin/env python3.11
"""
The '1/3' in C = 4pi + 1/3: is it the quark's SU(3) three-ness? (user 2026-09-22: fractional charge
-1/3/+2/3, baryon number 1/3.) All the 1/3's are the SAME SU(3) signature. Check the DOUBLET net
charge (the sum rule sums the doublet) and the O6 target. Not a derivation of 4pi -- a sharpening.
"""
import numpy as np
PI=np.pi
me,mmu=0.51099895,105.6583755
C1=(2.16+4.67)/me; C2=(1270.+93.4)/mmu

print("="*66,"\nThe three 1/3's are ONE SU(3) signature\n","="*66,sep="")
print("  quark electric charges: up +2/3, down -1/3")
print("  DOUBLET (up,down) NET charge = +2/3 + (-1/3) = +1/3   <- the sum rule sums the DOUBLET")
print("  baryon number per quark     = +1/3  (3 quarks -> baryon 1)")
print("  conformal weight h_SU(3)_1  = 1/3   (fundamental rep)")
print("  => all = the quark's SU(3) THREE-ness (charge quantized in thirds, 3 colours, 3 quarks/baryon).")
print("  The DOUBLET NET CHARGE (+1/3) is the cleanest match: it is exactly what the sum rule (m_u+m_d)")
print("  is a budget for, and it is GENERATION-INDEPENDENT (every (u,d) doublet has net charge +1/3).")

print("\n"+"="*66,"\nC = 4pi + 1/3 : consistency across doublets\n","="*66,sep="")
Cth=4*PI+1/3
print(f"  C_theory = 4pi + 1/3 = {Cth:.4f}")
print(f"  g1 (u,d/e):  measured {C1:.2f}   dev {abs(C1-Cth)/Cth*100:.1f}% (u,d ~10% uncertain -> loose)")
print(f"  g2 (c,s/mu): measured {C2:.3f}  dev {abs(C2-Cth)/Cth*100:.2f}% (reliable)")
print(f"  every (up,down) doublet has net charge +1/3 -> C = 4pi + 1/3 is GEN-INDEPENDENT, matching the")
print(f"  observed near-constant C. (baryon# would give doublet 2/3, worse; net charge 1/3 is the fit.)")

print("\n"+"="*66,"\nO6 TARGET (what the formation mass-balance must produce)\n","="*66,sep="")
print("  doublet budget:  (m_up + m_down) = [ 4pi + Q_doublet ] * m_lepton,  Q_doublet = +1/3")
print("  interpretation:  4pi   = the GEOMETRIC lepton-relaxation factor (solid angle of the director")
print("                          S^2; the quark doublet 'spreads' over the full sphere vs the lepton it")
print("                          conforms to -- the sparse-vacuum target).")
print("                   +1/3  = the quark SU(3) signature (net doublet charge / baryon-third / h_SU(3)).")
print("  So O6 (Q_H=2 lepton + Q_H=1 nu -> Q_H=3 doublet + byproduct) must give a mass-balance = a")
print("  4pi geometric (lepton) piece + a charge/SU(3) 1/3 piece. The charge piece is framework-native")
print("  (charge=writhe; the EM Q*Phi coupling, sec.26). The 4pi geometric piece is the OPEN part.")

print("\n"+"="*66,"\nHONEST STATUS\n","="*66,sep="")
print("  - The 1/3 = the quark's SU(3) three-ness (net doublet charge +1/3 the cleanest form) -- a")
print("    UNIFYING identification (charge, baryon#, conformal weight all the same 1/3). Group picture.")
print("  - 4pi = geometric (solid angle). C = 4pi + 1/3 = geometric + SU(3): the convergence the user saw.")
print("  - This SHARPENS the O6 target but does NOT derive 4pi (why exactly 4pi lepton masses). The O6")
print("    formation-geometry (director spreading over S^2) is the remaining piece; the 1/3 is now")
print("    understood (SU(3)). Caveat: still 1 reliable point + post-hoc; 13 within 1sigma (sec.36).")
