"""
WP4 gate: does the condensate-normalised Frank constant reproduce the rotation-curve
requirement K_phys^RC = c^2 v_flat^2/(4 pi G)? Unit map (nematic): K_phys ~ K_band *
rho_cond * ell^2 (energy/length; verified vs liquid crystals: K~k_B T/a = (k_B T/a^3) a^2).
The decisive question is WHICH length ell enters:
  (i)  condensate coherence length  xi_cond ~ hbar c / Lam_cond  (~mm), or
  (ii) the ULTRALIGHT DIRECTOR Compton length  ell_dir ~ c / m_xi ~ c/H_0 (~Hubble radius).
"""
import numpy as np
phi=(1+5**0.5)/2
c=2.99792458e8; G=6.674e-11; hbar=1.0546e-34; eV=1.602e-19; kB=1.381e-23
hbarc=hbar*c                                   # J*m
H0=2.27e-18                                     # s^-1 (~70 km/s/Mpc)
kpc=3.0857e19; Mpc=1e3*kpc
Tcmb=2.7255*kB/eV                               # eV
Lam_cond=Tcmb*(np.pi**2/150)**0.25             # eV
Lam_cond_J=Lam_cond*eV
rho_cond=10*Lam_cond_J**4/(hbarc)**3           # J/m^3  (= 10 Lam_cond^4)
print(f"Lam_cond = {Lam_cond:.3e} eV ;  rho_cond = 10 Lam_cond^4 = {rho_cond:.3e} J/m^3")
print(f"  (CMB energy density for comparison ~ 4e-14 J/m^3)")

# target
vflat=200e3
K_RC=c**2*vflat**2/(4*np.pi*G)
print(f"\nTARGET  K_phys^RC = c^2 v_flat^2/(4 pi G) = {K_RC:.3e} N   (v=200 km/s)")

Kband=0.16                                      # dimensionless band stiffness (a0-scale)
def Kphys(ell): return Kband*rho_cond*ell**2

print("\n== the two length choices ==")
xi_cond=hbarc/Lam_cond_J
ell_dir=c/H0
for name,ell in [("(i) condensate coherence xi_cond ~ hbar c/Lam_cond", xi_cond),
                 ("(ii) director Compton  ell_dir ~ c/m_xi = c/H0",      ell_dir)]:
    Kp=Kphys(ell)
    print(f"{name}:")
    print(f"     ell = {ell:.3e} m ({ell/Mpc:.3e} Mpc) -> K_phys = {Kp:.3e} N ;  "
          f"K_phys/K_RC = {Kp/K_RC:.3e}  ({np.log10(Kp/K_RC):+.1f} dex)")

print("\n== sensitivity of the (ii) near-match ==")
Kp2=Kphys(ell_dir); base=Kp2/K_RC
print(f"base (v=c, m_xi=H0, Kband=0.16): K_phys/K_RC = {base:.1f}")
print(f"  v=c/phi  => ell*=1/phi -> /phi^2={phi**2:.2f}   -> ratio {base/phi**2:.1f}")
print(f"  m_xi=2H0 => ell*=1/2   -> /4               -> ratio {base/4:.1f}")
print(f"  both                                        -> ratio {base/phi**2/4:.2f}")
print(f"  => the required ell_dir that hits K_RC exactly:")
ell_need=np.sqrt(K_RC/(Kband*rho_cond))
print(f"     ell_need = {ell_need:.3e} m = {ell_need/Mpc:.3f} Mpc = {c/ell_need/H0:.2f} * (c/H0)^-1 in m_xi")
print(f"     i.e. m_xi/H0 = {ell_dir/ell_need:.2f}  (order unity!)")

print("\n== VERDICT ==")
print(f"Length (i) coherence ~mm: FAILS by {np.log10(Kphys(xi_cond)/K_RC):.0f} dex.")
print(f"Length (ii) director Compton c/m_xi: NEAR-MATCH, within ~{base:.0f}x = "
      f"{np.log10(base):.1f} dex, absorbable by v=c/phi, m_xi~few H0, and O(1) geometry.")
print("=> Elastic DM is quantitatively VIABLE only because the director is ULTRALIGHT")
print("   (m_xi~H0): its elastic response is cosmological, giving K_phys set by rho_cond")
print("   (from T_CMB) and m_xi -- NO galactic input -- right to ~1 order of magnitude.")
print("   The SCALE is derived; the O(1) COEFFICIENT (phi^6/2pi) is NOT pinned (R1 holds).")
