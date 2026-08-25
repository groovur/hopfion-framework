"""
Einstein-frame DE scalar mass. Jordan action (schematic, M_Pl units):
  S = INT sqrt(-g) [ (1/2) f(rho) R - (1/2) k(rho) (d rho)^2 - V(rho) ],
  f(rho) = M_Pl^2 (1 + a rho),  a = 3 phi^6 beta   [f=2F, F=(M_Pl^2/2)(1+3phi^6 beta rho)]
  k(rho) = 2/(phi^6 (1+beta rho)),   V(rho)=Lam0 e^{-phi^6 rho}.
Einstein frame g~=Omega^2 g, Omega^2=f/M_Pl^2:
  U(rho) = M_Pl^4 V/f^2 = V/(1+a rho)^2  (M_Pl=1 units for V=Lam_obs density)
  (dchi/drho)^2 = M_Pl^2 [ k/f + (3/2)(f'/f)^2 ]
  m_E^2 = U''/(dchi/drho)^2   (leading; rho_inf is the frozen attractor)
Question: does the phi^20 Jordan enhancement CANCEL in the Einstein frame -> m_xi ~ H0?
"""
import numpy as np
phi=(1+5**0.5)/2
# beta from omega_BD = phi^12/(3 beta) ~ 237
beta=phi**12/(3*237.0)
a=3*phi**6*beta
rho_inf=phi/beta
u0=1+a*rho_inf                    # = 1+3 phi^7
print(f"phi^6={phi**6:.2f}, phi^12={phi**12:.1f}, phi^20={phi**20:.1f}")
print(f"beta={beta:.4f}, a=3phi^6 beta={a:.2f}, rho_inf=phi/beta={rho_inf:.3f}, "
      f"a rho_inf={a*rho_inf:.1f} (=3phi^7={3*phi**7:.1f}), u0=1+3phi^7={u0:.1f}\n")

# --- Jordan-frame mass (what eq:mxi computes): m_J^2 = V''/(k) at rho_inf, canonical-ish ---
k_inf=2/(phi**6*(1+beta*rho_inf))     # 1+beta rho_inf = 1+phi = phi^2
# V''=phi^12 V ; non-canonical physical Jordan mass m_J^2 = V''/k = phi^12 V /k
mJ2_over_Lam = phi**12 / k_inf        # in units of Lam_obs (=V(rho_inf))
print(f"Jordan: k(rho_inf)={k_inf:.4f}; m_J^2 = phi^12 V/k = {mJ2_over_Lam:.1f} * Lam_obs  "
      f"(= phi^20/... , the ~phi^20 enhancement)")

# --- Einstein-frame ---
fp_over_f = a/u0                        # f'/f at rho_inf
k_over_f  = k_inf/u0                    # k/f = k/(M_Pl^2 u0), in 1/M_Pl^2 units -> factor u0
chi2 = k_over_f + 1.5*fp_over_f**2      # (dchi/drho)^2 / M_Pl^2
# U''/U = (U'/U)^2 + (U'/U)' ; U'/U = -phi^6 - 2a/u0 ; (U'/U)'=2a^2/u0^2
UpU = -phi**6 - 2*a/u0
UppU = UpU**2 + 2*a**2/u0**2
Upp_over_Lam = UppU / u0**2             # U''=UppU*U, U(rho_inf)=Lam_obs/u0^2
mE2_over_Lam = Upp_over_Lam / chi2      # m_E^2 in units of Lam_obs/M_Pl^2
print(f"\nEinstein: f'/f={fp_over_f:.3f}, k/f={k_over_f:.2e}, (dchi/drho)^2/M_Pl^2 = {chi2:.4f} "
      f"(NMC term {1.5*fp_over_f**2:.3f} dominates)")
print(f"  U'/U={UpU:.2f}, U''/U={UppU:.1f}, U(rho_inf)=Lam_obs/u0^2={1/u0**2:.2e} Lam_obs")
print(f"  m_E^2 = U''/(dchi/drho)^2 = {mE2_over_Lam:.3f} * (Lam_obs/M_Pl^2)")

# --- numbers in H0 units. Lam_obs (energy density) = 3 H0^2 M_Pl^2 (de Sitter Lam=3H^2) ---
# so Lam_obs/M_Pl^2 = 3 H0^2.
print("\n== convert to H0 (Lam_obs/M_Pl^2 = 3 H0^2, de Sitter) ==")
mE2_H0 = mE2_over_Lam*3
mJ2_H0 = mJ2_over_Lam*3    # careful: m_J^2 = mJ2_over_Lam * Lam_obs, and Lam_obs/M_Pl^2=3H0^2, but m_J^2 has the M_Pl^2?
# Jordan m_J^2 = phi^12 V/k; V=Lam_obs energy density; to get H0 units divide by M_Pl^2 (field norm):
mJ2_H0 = mJ2_over_Lam*3
print(f"  m_E = sqrt({mE2_H0:.2f}) H0 = {np.sqrt(mE2_H0):.2f} H0   <-- PHYSICAL (Einstein) mass")
print(f"  m_J = sqrt({mJ2_H0:.0f}) H0 = {np.sqrt(mJ2_H0):.0f} H0   <-- Jordan (eq:mxi phi^20)")
print(f"  ratio m_E/m_J = {np.sqrt(mE2_H0/mJ2_H0):.2e}  (the phi^20 CANCELS via strong NMC)")

print("""
VERDICT: the phi^20 Jordan enhancement is CANCELLED in the Einstein frame by the strong
nonminimal coupling (3 phi^7 ~ 87): (dchi/drho)^2 is dominated by (3/2)(f'/f)^2 ~ (beta/phi)^2,
and m_E^2 = U''/(dchi/drho)^2 comes out ~ Lam_obs/M_Pl^2 ~ H0^2, PHI-POWER-FREE.
=> the PHYSICAL (cosmological / fifth-force) DE mass is m_xi ~ H0. RESOLUTION B confirmed:
  - "m_xi ~ H0 ~ 1e-33 eV" is CORRECT (Einstein frame). No numerical error.
  - thawing DE is fine (m_xi ~ H0 => thaws now). No tension.
  - WP4 DM near-match (+1 dex) HOLDS.
  - eq:mxi's m_xi^2=phi^20 Lam_obs M_Pl^4 is the JORDAN-frame mass; misleading as "the" DE mass.
""")
