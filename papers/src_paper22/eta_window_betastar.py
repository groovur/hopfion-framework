#!/usr/bin/env python3.11
"""
eta_window_betastar.py -- the eta lynchpin (Paper XXII sec:open item i: the physical value of
beta_eff = beta_* rho at the QCD epoch). CORRECTED after user input (2026-10-13).

TWO REJECTED SHORTCUTS (recorded so they are not re-tried):
  * beta_* = phi^6 : WRONG. Paper VII fixes beta_* ~= 0.452 (self-consistency J_fb/J_a=phi). The
    e^{-phi^6 rho} is a SEPARATE (dark-energy) term. "beta_eff~phi^6" is a numerical coincidence.
  * "93% saturated -> beta_eff ~ 13" : WRONG. Paper XX's 93% is a SPATIAL/ENERGY fraction (the
    part of the config with beta g_2 >> 1), explicitly ENGINE-DEPENDENT (P20:rem:feedback_
    saturation), NOT the value of beta_eff.

CORRECTED PICTURE (user: "it's not phi^6, however it IS multiplied by phi^6 in the denominator,
and beta_* would be modified by rho in the early universe"):
  (1) beta_* is NOT a universal constant. The 0.452 is its self-consistent value at the VACUUM/
      lepton configuration (J_fb/J_a=phi solved there). At the high-density QCD/baryon epoch the
      self-consistency re-solves -> beta_* is density-modified. So "beta_eff = 0.452 * rho_huge
      over-suppresses" (the previous conclusion) is RETRACTED: it used a constant beta_*.
  (2) The per-power suppression is 1/[phi^6 (1+beta_* rho)]. The phi^6 (Bogomolny/WZW, physical,
      frame-robust) is FIXED and always present; it does the BULK of the work.
  (3) The baryon (Q_H=3 trefoil) feedback is DEEPLY SATURATED (beta g_2 >> 1; Paper XX
      P20:rem:feedback_saturation -- the saturation that denies the quark a static saddle). So the
      eta-channel beta_eff sits in the SATURATED regime, O(10), not at a tuned point.

DECOMPOSITION (p=2): eta = eta_iso * [sin^4 / (phi^6 (1+beta_eff))]^2.
"""
import numpy as np
from scipy.optimize import brentq

phi=(1+np.sqrt(5))/2; phi6=phi**6; phi12=phi**12
eta_obs=6.1e-10; eta_iso=1.76e-4; sin4=8.0/15.0

def eta_of(beta_eff,p): return eta_iso*(sin4/(phi6*(1+beta_eff)))**p

if __name__=='__main__':
    print("="*74)
    print("  ETA LYNCHPIN -- corrected: fixed phi^{6p} bulk + SATURATED feedback (Paper XX)")
    print("="*74)
    # (2) the fixed phi^12 does the bulk
    base=eta_iso*(sin4/phi6)**2    # beta_eff=0, p=2 : only the fixed phi^12
    print(f"  p=2: fixed phi^12 = {phi12:.0f}  -> eta_iso*(sin4/phi6)^2 = {base:.2e}")
    print(f"       (the physical Bogomolny normalisation ALONE supplies all but a factor")
    print(f"        {base/eta_obs:.0f} of the needed suppression -- the user's 'multiplied by phi^6').")
    # (3) the saturated feedback supplies the remainder
    be=brentq(lambda b: eta_of(b,2)-eta_obs,1e-3,1e6)
    print(f"  remaining factor {base/eta_obs:.0f} = (1+beta_eff)^2  ->  beta_eff = {be:.1f}")
    print(f"       beta_eff = {be:.1f} >> 1 : SATURATED. KEY (user): baryogenesis IS trefoil formation")
    print(f"       (2+2->3+1 assembles the baryon), so the feedback on the baryon channel is the")
    print(f"       trefoil's OWN -- the very config Paper XX (P20:rem:feedback_saturation) saturates.")
    print(f"       Paper XX's gradient reading (beta*g_2) and the eta density reading (beta_*rho) are")
    print(f"       two descriptions of the SAME trefoil texture, related by rho ∝ kappa (Paper VII) --")
    print(f"       ONE saturation, not a coincidence of two. What is NOT fixed is the MAGNITUDE: the")
    print(f"       self-consistent trefoil beta_eff is determined in principle but resolution-gated,")
    print(f"       and rho ∝ kappa carries an unpinned constant. 'Saturated'(>>1) is necessary not")
    print(f"       sufficient -- eta ~ 1/beta_eff^2 needs the specific ~15.")
    print()
    print("  VERDICT: eta suppression = (physical, FIXED phi^12 ~322x) x (feedback (1+beta_eff)^2 ~255x,")
    print("  beta_eff~15). The phi^6 does the bulk; the feedback is the forming-trefoil's OWN (same")
    print("  object as Paper XX, via rho ∝ kappa), deeply saturated. beta_*rho~15 is the value the")
    print("  saturated trefoil config MUST RETURN, not an independent input. STILL OPEN [O]: (i) that")
    print("  PRECISE saturated value -- gated by the trefoil resolution wall + the rho ∝ kappa")
    print("  constant; (ii) the heavy-band power p. Reframed 'free window location' -> 'same saturated")
    print("  trefoil object, magnitude resolution-gated' -- a real narrowing, NOT a closed prediction.")
    print("="*74)
