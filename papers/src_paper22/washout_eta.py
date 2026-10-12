#!/usr/bin/env python3.11
"""
washout_eta.py -- Gamma-vs-H washout estimate for the reconnection baryogenesis eta.
DM/Quark crossover (src_paper21), step (d) redirect (notes (cc)): the 2+1->3 fusion is
EXOTHERMIC (downhill, no barrier), so eta's ~1e-10 is NOT a barrier -- it is set by the
competition between the reconnection rate Gamma and the Hubble rate H at the QCD epoch.

eta ~ eps * (n_carrier/s) * efficiency,  eps = CP asymmetry/event ~ O(delta_CP) ~ 0.1.
  Gamma = n sigma v  (reconnection = a STRONG-interaction rate: n~T^3, sigma~1/Lambda_QCD^2, v~1).
  H(T)  = 1.66 sqrt(g*) T^2 / M_Pl.
  Gamma>>H => deep equilibrium => strong washout, efficiency ~ H/Gamma (surviving fraction).
  Gamma<<H => frozen, efficiency ~ 1.
Brackets eta between the two limits and compares to observed eta_obs = 6.1e-10.
HONEST: order-of-magnitude; the true efficiency needs the full Boltzmann + the confinement-
transition (threading-irreversibility) dynamics, which a fixed-scale Gamma/H does not capture.
"""
import numpy as np

GeV = 1.0
MPl   = 1.22e19 * GeV        # Planck mass
eta_obs = 6.1e-10
eps   = 0.1                  # CP asymmetry per event ~ O(delta_CP) (update v)
n_over_s = 1e-2              # relativistic carrier number/entropy ~ 0.1/g*

def H(T, gstar=17.25):   return 1.66*np.sqrt(gstar)*T**2/MPl
def Gamma(T, Lam=0.2, gdof=40):
    n = 0.122*gdof*T**3          # relativistic number density n = zeta(3)/pi^2 g T^3
    sigma = 1.0/Lam**2           # geometric hadron cross-section ~ 1/Lambda_QCD^2
    v = 1.0
    return n*sigma*v

print("="*74)
print("  Gamma-vs-H WASHOUT ESTIMATE for reconnection baryogenesis eta")
print("="*74)
print(f"  observed eta = {eta_obs:.1e};  eps (CP/event) ~ {eps};  n/s ~ {n_over_s}")
print(f"\n  {'T (GeV)':>9} {'H (GeV)':>11} {'Gamma (GeV)':>12} {'Gamma/H':>10} {'regime':>14}")
for T in [1.0, 0.2, 0.15, 0.1, 0.05]:
    h=H(T); g=Gamma(T); r=g/h
    reg = 'equilibrium' if r>1 else 'frozen'
    print(f"  {T:>9.2f} {h:>11.2e} {g:>12.3e} {r:>10.2e} {reg:>14}")

# freeze-out T* where Gamma/H = 1
from scipy.optimize import brentq
Tstar = brentq(lambda T: Gamma(T)-H(T), 1e-25, 1e3)
print(f"\n  freeze-out scale (Gamma=H):  T* = {Tstar:.3e} GeV")
print(f"  (reconnection stays in EQUILIBRIUM until this T; T_QCD~0.15 GeV >> T*)")

# eta brackets at T_QCD
Tq=0.15; r=Gamma(Tq)/H(Tq)
eta_nowash = eps                      # transition-frozen, no washout: ~ eps
eta_fullwash = eps*n_over_s/r         # strong washout: efficiency ~ H/Gamma
print(f"\n  AT T_QCD=0.15 GeV:  Gamma/H = {r:.2e}  (DEEP equilibrium)")
print(f"  eta brackets:")
print(f"    no-washout (transition-frozen)  eta ~ eps          = {eta_nowash:.1e}   ({eta_nowash/eta_obs:.1e} x obs)")
print(f"    full strong-washout             eta ~ eps n/s H/G  = {eta_fullwash:.1e}   ({eta_fullwash/eta_obs:.1e} x obs)")
print(f"    OBSERVED                                            = {eta_obs:.1e}")
print("\n  READ: Gamma/H ~ 1e19 at QCD = DEEP equilibrium. eta is BRACKETED -- no-washout way")
print("  OVER-produces (eps~0.1 >> 6e-10), full-washout way UNDER-produces (<< 6e-10). Hitting")
print("  6e-10 needs a SPECIFIC intermediate efficiency the fixed-scale Gamma/H does NOT fix ->")
print("  eta NOT naturally pinned at the QCD scale; sensitive to the confinement-transition dynamics.")
