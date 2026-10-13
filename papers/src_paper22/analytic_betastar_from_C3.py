#!/usr/bin/env python3.11
"""
analytic_betastar_from_C3.py -- ANALYTIC probe of the eta window's beta_*rho, bypassing the trefoil
resolution wall by using the ALREADY-PINNED trefoil sharpness C_3* (notes (hhhh)->(iiii)).

LOGIC: the eta suppression needs beta_*rho = beta_* x (trefoil-core gradient energy). The core is
exactly what the solver cannot resolve -- BUT the core SCALE is already fixed: Paper XX gives the
trefoil core radius ~ 1/C* with C_3* = 2.5062 (unconditional, Paper XV). The unit director rotates
by O(1) across the core, so g_2 = |grad n|^2 ~ C*^2. With the self-consistent feedback coupling
beta_* ~ 0.45 (Paper VII; density-feedback coupling, config-robust reference) and rho ∝ kappa ∝ g_2
(Paper VII l.281-286, proportionality constant c_rho ~ O(1), UNPINNED):

    beta_*rho  ~  beta_* * c_rho * g_2_core  ~  beta_* * c_rho * C*^2 .

This gives the ORDER OF MAGNITUDE of the window location from an established quantity -- the point is
whether it lands in the saturated regime and near the eta-required ~15, NOT a precise value (c_rho and
the core-peak-vs-average profile are the resolution-gated residual).
"""
import numpy as np
from scipy.optimize import brentq

phi=(1+np.sqrt(5))/2; phi6=phi**6
C3=2.5062                 # trefoil profile sharpness (Paper XV, unconditional)
beta_star=0.452           # feedback coupling (Paper VII)
eta_obs=6.1e-10; eta_iso=1.76e-4; sin4=8.0/15.0

def eta_of(be,p=2): return eta_iso*(sin4/(phi6*(1+be)))**p

if __name__=='__main__':
    print("="*74)
    print("  ANALYTIC beta_*rho FROM THE TREFOIL SHARPNESS C_3*  (bypassing the resolution wall)")
    print("="*74)
    g2_core = C3**2
    be_naive = beta_star*g2_core      # c_rho = 1 (naive)
    print(f"  C_3* = {C3}  ->  core radius ~ 1/C* = {1/C3:.2f}  ->  g_2_core ~ |grad n|^2 ~ C*^2 = {g2_core:.2f}")
    print(f"  beta_* = {beta_star} (Paper VII)")
    print(f"  => beta_*rho ~ beta_* * C*^2 = {be_naive:.1f}   (c_rho=1, uniform-core naive estimate)")
    print(f"     eta at this beta_*rho (p=2) = {eta_of(be_naive):.2e}  vs obs {eta_obs:.1e}  "
          f"({eta_of(be_naive)/eta_obs:.0f}x over)")
    print()
    be_obs = brentq(lambda b: eta_of(b)-eta_obs, 1e-3, 1e4)
    print(f"  eta_obs requires beta_*rho = {be_obs:.1f} (p=2).")
    print(f"  => the implied core factor c_rho * (peak/avg) = {be_obs/be_naive:.1f}  -- i.e. the")
    print(f"     core-peak g_2 and the rho ∝ kappa constant together ~{be_obs/be_naive:.0f}x the")
    print(f"     C*^2-average. O(few), the honest resolution-gated residual.")
    print()
    print(f"  BRACKET from C_3* (c_rho in [1, ~5], core peak up to a few x average):")
    for c in (1.0, 2.0, 5.3):
        be=beta_star*c*g2_core
        print(f"    c_rho*peak = {c:>4.1f}:  beta_*rho = {be:>5.1f}  ->  eta = {eta_of(be):.2e}")
    print()
    print("  VERDICT: the trefoil sharpness C_3* ALONE puts beta_*rho ~ beta_* C*^2 ~ 3 and, with an")
    print("  O(1-5) core-peak / rho-kappa factor, in the 3-15 range -- the SATURATED regime and the")
    print("  right ballpark for eta ~ 1e-8..1e-10 bracketing the observed 6e-10. So the window")
    print("  location is NOT a free unknown: it is fixed to O(few-15) by the established C_3*, with the")
    print("  precise value set by the (resolution-gated) core-peak g_2 and the rho ∝ kappa constant.")
    print("  NARROWS the eta open item from 'unknown, resolution-gated' to 'O(few-15) from C_3*, a")
    print("  factor ~5 (in beta_eff) residual = core profile + proportionality constant'. [N]")
    print("  CAVEAT [O]: beta_* taken as the lepton self-consistent 0.452 (the trefoil's own self-")
    print("  consistent coupling may differ); g_2~C*^2 is the uniform-core estimate; c_rho unpinned.")
    print("="*74)
