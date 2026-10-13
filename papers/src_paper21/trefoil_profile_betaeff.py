#!/usr/bin/env python3.11
"""
trefoil_profile_betaeff.py -- single-tube (AXISYMMETRIC) estimate of beta_*rho. DM/Quark crossover.

*** CAVEAT (user, 2026-10-13): this is the TORUS reduction and is NOT valid for the trefoil. ***
The 1D radial profile f_0(rho) takes the field as a function of distance to the centerline -- exact
for the axisymmetric Q_H=2 lepton TORUS, but for the (2,3) TREFOIL it DISCARDS the self-interaction
at the 3 crossings (where two tube strands overlap and g_2 / the feedback are LARGEST) and the
knotting/winding. So the number below (beta_*rho_core ~ 0.54) is only the SINGLE-TUBE leading term --
a LOWER BOUND; the crossings lift beta_*rho above it, toward the saturated regime Paper XX's 3D
tube-adapted engine actually captures. Do NOT read 0.54 as the trefoil value, and do NOT read a
"tension" from it (that was an artefact of omitting the crossings). The real beta_*rho needs the 3D
self-interacting field (crossings resolved); the 3D-construction / resolution gate stands.
Kept only as the single-tube baseline + a validation that the profile machinery reproduces V_max.

WHY it was tempting: the framework's trefoil quantities (V_max, C_3*, J-integrals) use the 1D BS
profile (Paper XV, P15:prop:kappa4) -- but there it is the LEADING term, with the crossing/curvature
corrections (kappa^2, kappa^4) added separately. For beta_*rho the crossings are not a small
correction -- they are where the feedback saturates -- so the 1D term alone is inadequate.

CONVENTION (Paper XV eq In_beta, P15:prop:kappa4): normalised tube coord rho_h = rho*C*,
  f_0(rho_h) = 2 arctan(1/rho_h),  f_0'(rho_h) = -2/(1+rho_h^2),  sin f_0 = 2 rho_h/(1+rho_h^2).
  I_n = int_0^inf sin^{2n}f_0 (f_0')^2 rho_h d rho_h = 2^{2n+1} B(n+1,n+1);  I_2 = 16/15.
Gradient energy (radial): g_2 = |grad n|^2 = C*^2 (f_0')^2 (physical; C*^2 from drho_h/drho).
Feedback combination:  beta*g_2 = (beta C*^2)(f_0')^2 == B (f_0')^2,  B := beta C*^2.
Self-consistency (same J_fb/J_a=phi as the lepton's beta_*=0.452):
  V(B) = J2iso_fb/J2a = [int (f_0')^2/(1+B(f_0')^2) rho_h drho_h] / I_2  ->  set = phi.
Then the eta-channel beta_*rho at the tube CORE (where the baryon-forming reconnection sits,
rho_h->0, (f_0')^2 -> 4) is  beta_*rho_core = B_sc * 4 ; the energy-weighted value is also reported.
VALIDATION: I_2 -> 16/15, I_0 -> 2, V_max = I_0/I_2 -> 1.875 (vs Paper XX's 1.91).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

phi=(1+np.sqrt(5))/2
def fp2(rh): return (2.0/(1.0+rh**2))**2                 # (f_0')^2  (radial gradient)
def sin2(rh): return (2*rh/(1+rh**2))**2                 # sin^2 f_0
def sin4(rh): return (2*rh/(1+rh**2))**4                 # sin^4 f_0

def g2hat(rh, m):
    """full normalised gradient energy: radial (f_0')^2 + poloidal m^2 sin^2 f_0 / rho_h^2.
       m=0 -> radial only (incomplete); m>=1 -> includes the winding term that saturates."""
    polo = (m**2)*sin2(rh)/np.maximum(rh,1e-12)**2
    return fp2(rh) + polo

def I_n(n):                                              # radial-only check: 2^{2n+1}B(n+1,n+1)
    integ=lambda rh: (2*rh/(1+rh**2))**(2*n) * fp2(rh) * rh
    return quad(integ,0,np.inf)[0]

def J2iso_fb(B, m):
    integ=lambda rh: g2hat(rh,m)/(1+B*g2hat(rh,m)) * rh
    return quad(integ,0,np.inf)[0]
def J2a_full(m):
    integ=lambda rh: sin4(rh)*g2hat(rh,m) * rh
    return quad(integ,0,np.inf)[0]

if __name__=='__main__':
    I0=I_n(0); I2=I_n(2)
    print("="*78)
    print("  TREFOIL beta_*rho from the 1D BS profile (Paper XV method; no 3D charge-leak)")
    print("="*78)
    print(f"  VALIDATION (radial-only I_n): I_2 = {I2:.5f} (exact 16/15={16/15:.5f}), I_0 = {I0:.5f} (2)")
    print(f"              radial-only V_max = I_0/I_2 = {I0/I2:.4f} (Paper XX 1.91; 1D 1.875) OK")
    print()
    eta_iso=1.76e-4; sin4th=8.0/15.0; phi6=phi**6
    def eta(be): return eta_iso*(sin4th/(phi6*(1+be)))**2
    print(f"  full gradient g2hat = (f_0')^2 + m^2 sin^2 f_0 / rho_h^2  (m = poloidal winding)")
    print(f"  self-consistency V(B)=phi={phi:.4f}; beta_*rho_core = B_sc * g2hat(core), g2hat(core)=4(1+m^2)")
    print(f"  {'m':>3} {'Vmax(full)':>11} {'B_sc':>9} {'g2hat_core':>11} {'beta_*rho_core':>15} {'eta(p=2)':>11} {'xobs':>7}")
    for m in (0,1,2,3):
        J2a=J2a_full(m)
        def V(B): return J2iso_fb(B,m)/J2a
        Vmax=V(0.0)
        if Vmax>phi:
            B_sc=brentq(lambda B: V(B)-phi,1e-8,1e6)
        else:
            B_sc=float('nan')
        g2c=4*(1+m**2)
        be_core=B_sc*g2c
        e=eta(be_core) if be_core==be_core else float('nan')
        print(f"  {m:>3} {Vmax:>11.3f} {B_sc:>9.4f} {g2c:>11.0f} {be_core:>15.2f} {e:>11.2e} {e/6.1e-10:>7.0f}")
    print()
    print("  eta needs beta_*rho ~ 15 (eta_obs=6.1e-10, p=2).")
    print("  READ: m=0 (radial only) UNDER-counts -> beta_*rho_core ~ 0.5 (would over-produce ~100x).")
    print("  The POLOIDAL winding term (m>=1) is the saturation driver (Paper XX's '93% saturated');")
    print("  it raises beta_*rho_core steeply with m. Whether the trefoil's actual winding reaches the")
    print("  m that gives beta_*rho_core ~ 15 (hence eta_obs) is the open question -- the winding m of")
    print("  the (2,3)-cable tube is NOT pinned here (the build_doublet phase uses q=3, suggesting m~3).")
    print("  So: the 1D radial piece is a firm LOWER bound (~0.5); the poloidal piece plausibly carries")
    print("  it to the saturated ~15 for m~3, but the exact winding + full self-consistency is [O].")
    print("="*78)
