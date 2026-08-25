import numpy as np
phi=(1+5**0.5)/2; H0=1.494e-33  # eV
r3=np.sqrt(3)
print("(a) does V have a minimum / is V''=phi^12 V?")
print("  V(rho)=Lam0 e^{-phi^6 rho}: V'=-phi^6 V (never 0), V''=+phi^12 V. PURE RUNAWAY, NO minimum.")
print("  'potential minimum' = the ATTRACTOR rho_inf=phi/beta (frozen by Hubble friction, not V'=0).")
print(f"  V''(rho_inf)=phi^12 Lam_obs. phi^12={phi**12:.0f}. No flattening from V itself.")
print("  Only an Einstein-frame / nonminimal-coupling redefinition could change V''; the bare V does not.\n")
print("(b) canonical vs non-canonical DE mass (does the 1/phi^8 = K^{-1} belong?):")
print("  L=K(drho)^2 - V, K(rho_inf)=1/phi^8. Canonical chi: dchi=sqrt(2K) drho.")
print("  m_chi^2 = V''/(2K) = phi^12 Lam_obs * phi^8/2 = phi^20 Lam_obs/2  -> the phi^8 DOES enter.")
for lab,fac in [("non-canonical m_xi=phi^10 r3 H0", phi**10*r3),
                ("canonical-only (drop phi^8): phi^6 r3 H0", phi**6*r3)]:
    print(f"    {lab:42s} = {fac:6.1f} H0 = {fac*H0:.2e} eV")
print("  => even DROPPING the phi^8 (pure canonical) gives ~30 H0, still >> H0. Neither rescues m_xi~H0.\n")
print("VERDICT: (a) and (b) both CONFIRM the enhancement. m_xi is 30-210 H0 ~ 1e-31 eV, robustly >> H0.")
print("So consequences 1 (value) and 3 (thawing) are REAL and COUPLED, not clean typos.")
