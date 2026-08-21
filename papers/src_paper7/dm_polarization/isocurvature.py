"""Isocurvature number for the self-organised director texture.

The splay texture is DM only once the director UNFREEZES and self-organises.
An ultralight director (mass m_xi ~ H_0, the DE scalar mass at the attractor,
Eq. m_xi^2 = phi^20 Lam_obs M_Pl^4) is Hubble-frozen while H >> m_xi and only
starts to evolve/order when H(z) ~ m_xi. So the texture-DM density present at a
given epoch is gated by z_order.

beta_iso at the CMB is therefore (texture DM fraction present at recombination)
x (frozen director fluctuation)^2 / adiabatic power. The first factor is the
decisive one.
"""
import numpy as np

# cosmology
Om, OL = 0.315, 0.685
def Hratio(z):            # H(z)/H_0
    return np.sqrt(Om*(1+z)**3 + OL)

# --- when does the ultralight director unfreeze / order? H(z_order) = m_xi ---
# framework: m_xi ~ H_0 (paper: m_xi ≈ H_0 ≈ 1e-33 eV). scan c1 = m_xi/H_0 ~ O(1).
print("Director unfreezing (H(z)=m_xi):")
for c1 in (0.84, 1.0, 1.2, 2.0):
    # solve Hratio(z)=c1
    if c1 <= np.sqrt(OL):
        z = np.nan  # never (H asymptotes to sqrt(OL) H0)
    else:
        z = ((c1**2 - OL)/Om)**(1/3) - 1
    print(f"   m_xi/H_0 = {c1:4.2f} -> z_order = {z if np.isnan(z) else round(z,3)}"
          + ("  (never; H floor = %.2f H0)"%np.sqrt(OL) if np.isnan(z) else ""))

print("\n=> for m_xi ~ H_0, the director unfreezes at z_order ~ 0 (now/near future).")
print("   The self-organised texture is a PRESENT-EPOCH phenomenon, essentially")
print("   absent at recombination (z=1100).")

# --- frozen inflationary director fluctuation (the would-be isocurvature seed) ---
# H_inf from P_s = H_inf^2/(8 pi^2 eps M_Pl^2), Starobinsky eps = 2/N_e^2.
Ps, Ne = 2.10e-9, 57.73
eps = 2/Ne**2
Hinf_over_MPl = np.sqrt(Ps * 8*np.pi**2 * eps)
phi = (1+5**0.5)/2
f_dir_over_MPl = 1/phi**6         # director decay constant ~ condensate scale Lam_UV ~ M_Pl/phi^6
dtheta = Hinf_over_MPl/(2*np.pi*f_dir_over_MPl)
print(f"\nFrozen director fluctuation seed:")
print(f"   eps={eps:.2e}  H_inf/M_Pl={Hinf_over_MPl:.2e}  f_dir/M_Pl={f_dir_over_MPl:.3f}")
print(f"   delta_theta ~ H_inf/(2 pi f_dir) = {dtheta:.2e}")
print(f"   (adiabatic seed sqrt(P_s) = {np.sqrt(Ps):.2e} for comparison)")

# --- beta_iso at CMB ---
# texture DM density present at z=1100 relative to today ~ 0 (forms at z~0):
f_tex_at_cmb = 0.0   # to the precision of the ultralight unfreezing argument
beta_iso_cmb = f_tex_at_cmb * (dtheta/np.sqrt(Ps))**2
print(f"\nbeta_iso(CMB) = f_tex(z=1100) x (delta_theta/sqrt(P_s))^2 = {beta_iso_cmb:.1e}")
print(f"   Planck bound: beta_iso < 0.038.  SATISFIED trivially (texture is late).")

print("\n== HONEST VERDICT ==")
print("beta_iso(CMB) ~ 0: the ultralight director (m_xi~H_0) unfreezes at z~0, so")
print("the texture forms NOW and carries no CMB isocurvature -> safely below Planck.")
print("COST (the same lateness): this mechanism supplies NO DM at recombination, so")
print("it CANNOT be the CMB-era cold dark matter (peak heights, early growth). The")
print("director-strain DM is a LATE-TIME, galactic-scale component; early-universe")
print("CDM is a separate, still-open account. The isocurvature safety and the")
print("early-structure gap are the SAME fact (lateness), seen from two sides.")
print("NB: the frozen delta_theta ~ %.0e is itself comparable to sqrt(P_s); it is"%dtheta)
print("harmless ONLY because it does not gravitate as DM until z~0. If any texture")
print("DM were present at recombination, beta_iso would be O(1) -- excluded. The")
print("lateness is thus REQUIRED, not incidental.")
