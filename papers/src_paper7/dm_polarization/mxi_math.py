"""
Run down the m_xi = phi^10 sqrt3 H0 vs "m_xi ~ H0 ~ 1e-33 eV" inconsistency FROM THE MATH.
Parent action: kinetic K(rho)(d rho)^2 with K(rho)=1/[phi^6(1+beta rho)]; V(rho)=Lam0 e^{-phi^6 rho}.
Physical scalar mass of fluctuations: m_xi^2 = V''(rho_inf)/K(rho_inf) (non-canonical norm).
"""
import numpy as np
phi=(1+5**0.5)/2
H0=2.27e-18; hbar_eV_s=6.582e-16; H0_eV=hbar_eV_s*H0   # H0 in eV

print("== the phi^20 is REAL: it comes from V'' and the non-canonical K ==")
print("V(rho)=Lam0 e^{-phi^6 rho} => V'=-phi^6 V, V''=phi^12 V.  At rho_inf: V''=phi^12 Lam_obs.")
print("K(rho_inf)=1/[phi^6(1+beta rho_inf)] , beta rho_inf=phi , 1+phi=phi^2 => K=1/(phi^6 phi^2)=1/phi^8.")
print(f"m_xi^2 = V''/K = (phi^12 Lam_obs)(phi^8) = phi^20 Lam_obs.   phi^12={phi**12:.1f}, phi^8={phi**8:.1f}, phi^20={phi**20:.1f}")
print("=> the phi^20 = phi^12 (potential curvature) x phi^8 (kinetic normalisation). GENUINE, not a typo.\n")

print("== numerical m_xi, using Lam_obs -> de Sitter H0^2 (Lam=3H^2) ==")
# de Sitter: Lam_obs (dimensionless, =Lam/M_Pl^... ) enters as m_xi = phi^10 sqrt(Lam_obs) M_Pl = phi^10 sqrt3 H0
mxi_eV = phi**10*np.sqrt(3)*H0_eV
print(f"H0 = {H0_eV:.3e} eV")
print(f"m_xi = phi^10 sqrt3 H0 = {phi**10*np.sqrt(3):.1f} H0 = {mxi_eV:.3e} eV")
print(f"  vs the paper's stated 'm_xi ~ H0 ~ 1e-33 eV': WRONG by factor phi^10 sqrt3 = {phi**10*np.sqrt(3):.0f}")
print(f"  correct order is ~1e-31 eV, NOT 1e-33 eV.\n")

print("== dimensional check of eq:mxi as written ==")
print("eq:mxi: m_xi^2 = phi^20 Lam_obs M_Pl^4, with Lam_obs/M_Pl^4 ~ 1e-122 (Lam_obs is energy density E^4).")
print("Then m_xi^2 = phi^20 (E^4)(E^4)=E^8 -> m_xi ~ E^4. DIMENSIONALLY WRONG.")
print("The number-correct form is m_xi^2 = phi^20 (Lam_obs/M_Pl^4) M_Pl^2  [dimensionless CC x M_Pl^2],")
print("i.e. m_xi = phi^10 sqrt(Lam_obs/M_Pl^4) M_Pl = phi^10 sqrt3 H0.  So 'M_Pl^4' in eq:mxi is a slip.\n")

print("== consequence 1: thawing DE ==")
# field thaws when H(z) ~ m_xi. H(z)=H0 sqrt(Om(1+z)^3+OL), Om=0.3, OL=0.7
r=phi**10*np.sqrt(3)   # m_xi/H0
Om,OL=0.3,0.7
zthaw=((r**2-OL)/Om)**(1/3)-1
print(f"m_xi/H0 = {r:.0f}. Field thaws at H(z)=m_xi => (1+z)^3=(r^2-OL)/Om => z_thaw = {zthaw:.0f}.")
print("A thawing field must be FROZEN until ~now (m_xi<~H0). z_thaw~50 means it thawed in the")
print("matter era and has ROLLED since -> tensions the proved thawing relation w_a=-3(1+w0).")
print("For genuine thawing NOW you need m_xi ~ H0 (NO phi^10). CONFLICT.\n")

print("== consequence 2: WP4 DM Frank constant (v=c/phi) ==")
c=2.99792458e8;G=6.674e-11;eV=1.602e-19;kB=1.381e-23;hbar=1.0546e-34;hbarc=hbar*c
Lam_cond=(2.7255*kB/eV)*(np.pi**2/150)**0.25; rho_cond=10*(Lam_cond*eV)**4/hbarc**3
K_RC=c**2*(200e3)**2/(4*np.pi*G); Kband=0.16
for lab,mfac in [("m_xi=H0 (loose)",1.0),("m_xi=phi^10 sqrt3 H0 (correct)",r)]:
    ell=(c/phi)/(mfac*H0); Kp=Kband*rho_cond*ell**2
    print(f"  {lab:32s}: K_phys/K_RC = {Kp/K_RC:.2e} ({np.log10(Kp/K_RC):+.1f} dex)")
print("=> with the CORRECT m_xi~213 H0, the DM Frank constant undershoots by ~3-4 dex.")
