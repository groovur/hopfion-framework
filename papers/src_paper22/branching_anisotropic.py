#!/usr/bin/env python3.11
"""
branching_anisotropic.py -- REDO the baryon branching b_ch WITHOUT isotropising the sin^4 theta
(user catch; claude-hopfion s10 trap #1). notes (tt)->(uu). DM/Quark crossover src_paper21.

branching_eta.py took b_ch~0.2 as a direction-BLIND channel count. WRONG: the baryon (E_6, Q=3)
is the HEAVY band, and the heavy band CARRIES the directional suppression (claude-hopfion s3)
  S_eff = sin^4(theta) / [ phi^6 (1 + beta_* rho) ]        (Axiom A2)
while the competing light-band channels (elastic 2+2; the neutrino byproduct Q=1) do NOT. So the
baryon-forming reconnection (light leptons -> HEAVY baryon) projects onto the heavy band and its
rate is suppressed by S_eff. Isotropically averaging sin^4 (=8/15) washes out exactly the phi^6
(and density-feedback) suppression that the estimate was missing.

RESOLVED branching:  b_ch^res = b_ch^iso * S_eff^p ,  p = # heavy-band projections in the final
state. Final state 3+1 has ONE heavy object (the baryon) + one light (nu) => p=1 (natural); p=2
shown as a bracket. phi^6 is a PHYSICAL/topological power (WZW/Bogomolny lambda=phi^6, claude-
hopfion s6 discriminator -- NOT a frame artifact), so this suppression is real. (1+beta_* rho) is
the density feedback at the QCD epoch -- cutoff-set, spans orders (FLAG); it is the lever.
Report eta(iso) vs eta(resolved) across (1+beta_* rho) and p; locate where eta_obs sits.
"""
import numpy as np
phi=(1+np.sqrt(5))/2; phi6=phi**6                  # 17.944 -- WZW/Bogomolny power (physical)
eta_obs=6.1e-10
# from branching_eta.py (isotropic assembly):
eta_iso=1.76e-4           # A_abs=1 case; over-produces by ~2.9e5
b_ch_iso=0.20
sin4_avg=8.0/15.0         # <sin^4 theta> over solid angle (heavy band sits near the hard axis)

def S_eff(oneplus_brho, sin4=sin4_avg): return sin4/(phi6*oneplus_brho)

if __name__=='__main__':
    print("="*78)
    print("  ANISOTROPIC baryon branching: b_ch carries S_eff=sin^4/[phi^6(1+beta_* rho)] (heavy band)")
    print("="*78)
    print(f"  phi^6 = {phi6:.3f} (physical WZW/Bogomolny power, s6);  <sin^4> = {sin4_avg:.3f}")
    print(f"  ISOTROPIC (wrong): b_ch={b_ch_iso}, eta={eta_iso:.2e} => over-produces by {eta_iso/eta_obs:.1e}\n")
    print("  RESOLVED: eta^res = eta_iso * S_eff^p   (S_eff = <sin^4>/[phi^6 (1+beta_* rho)])")
    print(f"  {'1+beta*rho':>12} {'S_eff':>11} {'eta(p=1)':>11} {'ratio':>9} {'eta(p=2)':>11} {'ratio':>9}")
    for brho in [1.0, 5.0, 25.0, 1e2, 1e3, 1e4]:
        Se=S_eff(brho)
        e1=eta_iso*Se; e2=eta_iso*Se**2
        print(f"  {brho:>12.0f} {Se:>11.3e} {e1:>11.2e} {e1/eta_obs:>9.1e} {e2:>11.2e} {e2/eta_obs:>9.1e}")

    # where does eta_obs land?
    # p=1: eta_iso*<sin4>/(phi6*brho) = eta_obs  -> brho = eta_iso*<sin4>/(phi6*eta_obs)
    brho1=eta_iso*sin4_avg/(phi6*eta_obs)
    # p=2: eta_iso*(<sin4>/(phi6*brho))^2 = eta_obs -> brho = <sin4>/phi6 * sqrt(eta_iso/eta_obs)
    brho2=(sin4_avg/phi6)*np.sqrt(eta_iso/eta_obs)
    print(f"\n  eta_obs reached at:  p=1 -> (1+beta*rho) = {brho1:.1f}   p=2 -> (1+beta*rho) = {brho2:.1f}")

    print("\n  READ:")
    print("  * The ISOTROPIC estimate washed out the phi^6 directional suppression (user catch).")
    print(f"    Resolving it, the baryon channel (HEAVY band) carries S_eff -- one factor of 1/phi^6")
    print(f"    = {1/phi6:.3f} even at zero density feedback, pulling eta DOWN by ~{phi6:.0f}x per power p.")
    print("  * This CONVERTS the ~1e4-1e5 over-production into a DENSITY-FEEDBACK-controlled quantity:")
    print(f"    eta_obs sits at (1+beta*rho) ~ {brho1:.0f} (p=1) or ~ {brho2:.1f} (p=2) -- modest, NOT the")
    print("    extreme early-epoch beta*rho. If beta*rho is LARGER (claude-hopfion/memory flags it can")
    print("    be huge, ~rho~T^4), S_eff->0 and the channel UNDER-produces: the SAME density screening")
    print("    that turns the director/DM sector off in dense regions also throttles baryon formation.")
    print("  * So the honest picture: anisotropy supplies a REAL phi^6 suppression (physical power),")
    print("    recasting eta from 'over-produces' to 'set by beta*rho at the QCD epoch' -- a cutoff-set")
    print("    lever spanning orders. eta still unpinned, but now via ONE physical knob (the epoch")
    print("    density feedback), not a generic transport fudge. [O] beta*rho(T_QCD) is the quantity.")
