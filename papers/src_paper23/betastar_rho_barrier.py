#!/usr/bin/env python3.11
"""
betastar_rho_barrier.py -- Paper XXIII component F: calibrate the engine's feedback to the
S_eff convention and report beta_*rho for the trefoil -> the eta magnitude. DM/Quark crossover.

#############################################################################################
# RETRACTED NUMBER (notes (xxxx)->(zzzz), 2026-10-14). The CALIBRATION LOGIC below is correct #
# (beta_*rho = beta_* kappa = engine beta*g2, c cancels). BUT this script evaluates it at     #
# beta=0.452 on the RAW grid g2, which is NOT self-consistent: it gives V=0.30, not V=phi.    #
# The ~15-40 it prints is therefore a UNITS ARTIFACT, not the self-consistent value. The      #
# self-consistent ISOLATED-knot feedback is ~0.5 = 8 b* (see knot_in_medium.py), a robust     #
# medium-invariant FLOOR. The eta-relevant beta_*rho~15 is NOT this intrinsic value -- it is  #
# the ENVIRONMENT (chameleon) feedback of the dense QCD-epoch condensate, bracketed by the    #
# first-order window and freeze-out-selected at T_c (src_paper22/beta_rho_window.py; Paper    #
# XXII P22:sec:open, recast). Use knot_in_medium.py for the corrected picture. This script is #
# kept only to document the calibration logic and the artifact; do NOT cite its ~15 as a      #
# computed eta result.                                                                        #
#############################################################################################

THE CALIBRATION (resolves notes (kkkk)'s "units don't match S_eff"):
  Paper VII P7:rem:rho_extension -- the FUNDAMENTAL (Papers I-IV) feedback factor is (1+beta_* kappa),
  kappa = the Faddeev-Niemi gradient density = |grad n|^2. This is EXACTLY what the engine computes:
  (1+beta*g2), g2=|grad n|^2, beta=beta_*. The Paper VII k-essence rewrite S_eff uses (1+beta rho) with
  rho=c*kappa and beta=beta_*/c (P7:eq:Lambda_obs l.1390), so
        1 + beta*rho  =  1 + (beta_*/c)(c kappa)  =  1 + beta_* kappa     -- THE c CANCELS.
  => the S_eff feedback ARGUMENT  beta_*rho  =  beta_* kappa  =  the engine's  beta*g2  DIRECTLY.
  The normalisation c = phi^3 + p0^{Q_H} (P7:prop:cc_norm, p0=1/phi^2) enters ONLY the separate
  cosmological-constant exponent Lambda_obs = Lambda_0 e^{-phi^7 c/beta_*}, NOT the feedback factor.
  (This was the (kkkk) confusion: c was wrongly thought to rescale beta_*rho; it does not.)

THE ETA CONNECTION (Paper XXII P22:sec:open, P22:eq:seff):
  S_eff = sin^4theta/[phi^6 (1+beta_*rho)]; eta carries S_eff^2 (rate ~ |amplitude|^2, p=2, Paper XXI)
  -> suppression = phi^{12} (1+beta_*rho)^2. Reconnection baryogenesis OVERPRODUCES by ~1e4-1e5, so a
  suppression ~1e4-1e5 is needed; the paper's self-consistent beta_*rho~15 gives eta ~ 6e-10.

THIS SCRIPT combines the two computed loci (feedback_energy_3d = single-tube WALL; threading_barrier =
the formation BARRIER) under the calibration and reports beta_*rho and the implied eta suppression,
with the honest error budget (resolution-creep + unrelaxed ansatz).
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feedback_energy_3d import run as wall_run
from threading_barrier import run_d as barrier_run

phi = (1+np.sqrt(5))/2
ETA_OBS = 6.0e-10

def c_norm(Q_H):
    p0 = 1/phi**2
    return phi**3 + p0**Q_H

def eta_suppression(bstar_rho):
    """the framework's eta-suppression denominator phi^12 (1+beta_*rho)^2 (Paper XXII)."""
    return phi**12 * (1 + bstar_rho)**2

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=96)
    ap.add_argument('--L', type=float, default=6.0)
    ap.add_argument('--Lb', type=float, default=2.5, help='box for the threading model')
    ap.add_argument('--C_star', type=float, default=2.5062); ap.add_argument('--beta', type=float, default=0.452)
    ap.add_argument('--d_barrier', type=float, default=0.8)
    a = ap.parse_args()

    print("="*86)
    print("  (F) beta_*rho FOR THE TREFOIL -> eta  [calibration: beta_*rho = engine beta*g2; c CANCELS]")
    print("="*86)
    c = c_norm(3)
    print(f"  normalisation c = phi^3 + p0^Q_H (P7:prop:cc_norm), Q_H=3, p0=1/phi^2 = {1/phi**2:.4f}")
    print(f"     c(Q_H=3) = {c:.4f}  (enters Lambda_obs ONLY; CANCELS in the feedback factor -> NOT in beta_*rho)")
    print(f"  coupling beta_* = {a.beta} (FN self-consistent, V=phi); feedback arg beta_*rho = beta_* * |grad n|^2")
    print()

    # (A) single-tube wall (relaxed trefoil self-avoids)
    ow = wall_run(a.N, a.L, 3.0, 0.874, a.C_star, a.beta)
    bstar_rho_wall = ow['bg2_wall']
    # (B) threading barrier
    ob = barrier_run(a.d_barrier, a.N, a.Lb, a.C_star, a.beta)
    bstar_rho_barrier = ob['bg2_overlap']
    enh = ob['enh']

    print(f"  LOCUS                    beta_*rho    phi^12(1+b_*rho)^2    eta ~ eta_raw / suppression")
    for name, br in [("(A) single-tube wall", bstar_rho_wall),
                     ("(B) threading barrier", bstar_rho_barrier),
                     ("    paper XXII target", 15.0)]:
        sup = eta_suppression(br)
        # eta_raw (unsuppressed reconnection) overproduces by ~1e4-1e5 => eta_raw ~ ETA_OBS * 1e5 (mid)
        eta_mid = (ETA_OBS*8.0e4)/sup
        print(f"  {name:<24} {br:>8.2f}   {sup:>18.3e}   {eta_mid:>10.2e}")
    print()
    print(f"  barrier/wall enhancement = {enh:.2f}x")
    print("="*86)
    print("  *** RETRACTED NUMBER (see header + knot_in_medium.py). The ~15-40 above uses beta=0.452 on")
    print("  raw g2, which is NOT self-consistent (V=0.30, not phi); it is a UNITS ARTIFACT, not a")
    print("  computed eta result. The CALIBRATION LOGIC (beta_*rho = beta_* kappa = engine beta*g2, c")
    print("  cancels) is [E], grounded in P7:rem:rho_extension + prop:cc_norm -- but must be evaluated at")
    print("  the SELF-CONSISTENT coupling. Doing so (knot_in_medium.py): the isolated-knot beta_*rho is a")
    print("  robust, medium-invariant FLOOR ~0.5 = 8 b*, ~30x below 15. The eta-relevant beta_*rho~15 is")
    print("  the ENVIRONMENT (chameleon) feedback of the dense QCD-epoch condensate, bracketed by the")
    print("  first-order window and freeze-out-selected at T_c (src_paper22/beta_rho_window.py; Paper XXII")
    print("  P22:sec:open, recast) -- NOT the intrinsic trefoil value this script evaluates. ***")
    print("="*86)
