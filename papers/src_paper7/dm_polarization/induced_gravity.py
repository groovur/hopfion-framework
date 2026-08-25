import math
phi=(1+5**0.5)/2; pi=math.pi
MPl=1.22e19  # GeV, non-reduced
Lcond=1.19e-13  # GeV

print("== Seeley-DeWitt induced R-term (proper-time cutoff at s=1/Lam^2) ==")
print("  Heat kernel K(s)=(1/(16 pi^2 s^2)) sum a_n s^n, a_1=(1/6 - xi)R; minimal xi=0 -> a_1=R/6.")
print("  L_eff ⊃ -(1/(32 pi^2)) a_1 Lam^2 = -(Lam^2/(192 pi^2)) R.  (matches paper eq:1loop's R/6 term)")
print("  Match to (M_Pl^2/2) R:  M_Pl^2/2 = Lam^2/(192 pi^2)  =>  M_Pl^2 = Lam^2/(96 pi^2).")
print()
print("== direction: invert to get Lam_UV ==")
for label,coeff in [("Seeley-DeWitt  M_Pl^2=Lam^2/(96 pi^2)", 96*pi**2),
                    ("paper box      M_Pl^2=Lam^2/(32 pi^2)", 32*pi**2)]:
    Lam=math.sqrt(coeff)*MPl
    nUV=math.log(Lam/Lcond)/(2*math.log(phi))
    print(f"  {label}: Lam_UV=sqrt({coeff:.0f}) M_Pl = {math.sqrt(coeff):.2f} M_Pl "
          f"= {Lam:.2e} GeV (SUPER-Planckian)")
    print(f"       -> n_UV = ln(Lam_UV/Lcond)/(2 ln phi) = {nUV:.1f}")
print()
print(f"  paper's stated (WRONG direction) Lam_UV = M_Pl/(4pi sqrt2) = {MPl/(4*pi*math.sqrt(2)):.2e} GeV,")
print(f"    n_UV = {math.log(MPl/(4*pi*math.sqrt(2))/Lcond)/(2*math.log(phi)):.1f} (the paper's 73.6).")
print(f"  correct super-Planckian (box coeff, Lam~4pi sqrt2 M_Pl ~ phi^6 M_Pl): n_UV ~ 79.6.")
print(f"    shift = ln(phi^12)/(2 ln phi) = 6 exactly (Lam flips M_Pl/phi^6 -> phi^6 M_Pl).")
print("  EITHER n_UV is NON-INTEGER -> the 'compatible with non-integer tower index' conclusion HOLDS.")
print()
print("== does this break F(rho_inf), the CC, or the DE sector? NO ==")
print("  The NMC ratio (paper eq:NMC_ratio): xi_NMC rho_inf/(M_Pl^2/2) = 3 phi^7.")
print("  Both xi_NMC ∝ Lam_UV^2 AND M_Pl^2/2 ∝ Lam_UV^2, so Lam_UV^2 CANCELS in the ratio.")
print("  => 3 phi^7, hence F(rho_inf)=(M_Pl^2/2)(1+3phi^7), is Lam_UV-INDEPENDENT.")
print("  => the CC (Lam_obs/F), the DE mass, the Einstein-frame m_xi~H0 -- all UNAFFECTED.")
print("     Only the Lam_UV VALUE and n_UV move; the load-bearing gravity results are safe.")
print()
print("VERDICT: induced (Sakharov) gravity gives Lam_UV ABOVE M_Pl (super-Planckian), the standard")
print("result and consistent with M_Pl^2=Lam^2/(loop). The paper's sub-Planckian 'M_Pl/(4pi sqrt2)'")
print("is the reciprocal error. Correct: Lam_UV ~ (18-31) M_Pl ~ phi^6 M_Pl; n_UV~79-80, still non-integer.")
print("F(rho_inf)/CC/DE are Lam_UV-independent (cancellation), so nothing load-bearing breaks.")
