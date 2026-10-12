#!/usr/bin/env python3.11
"""
beta_rho_window.py -- the feedback strength beta_eff = beta_* rho is ONE variable that sets BOTH
(i) whether the confinement transition is first-order (hh: intermediate WINDOW; over-saturation
-> crossover) AND (ii) the baryon-branching suppression S_eff=sin^4/[phi^6(1+beta_eff)] -> eta
(uu). Overlay them on the SAME axis and find the beta_eff range where BOTH hold. notes (uu)->(vv).

User point: rho modifies beta HEAVILY in the QCD regime (beta_eff=beta_* rho, rho~T^4 ~1e47x
early; (gg)/(hh)). That heavy (1+beta_eff) IS the source of the baryon-branching suppression that
cures the (tt) over-production -- AND the first-order requirement (hh) caps beta_eff BELOW
over-saturation. So the two requirements CO-CONSTRAIN the same beta_eff(T_QCD).

(i) first-order barrier(beta_eff) from firstorder_transition.py (toy Landau, b=gamma=1 units).
(ii) eta(beta_eff) = eta_iso * S_eff^p,  S_eff = <sin^4>/[phi^6 (1+beta_eff)], p in {1,2}.
Report the overlap. HONEST: the toy's window LOCATION in physical units needs gamma/b from real
E_fb (still not extracted -- the residual open input); this shows the STRUCTURE + the consistency.
"""
import numpy as np
import firstorder_transition as fo

phi=(1+np.sqrt(5))/2; phi6=phi**6
eta_obs=6.1e-10; eta_iso=1.76e-4; sin4=8.0/15.0

def eta_of(beta_eff, p):
    S=sin4/(phi6*(1+beta_eff)); return eta_iso*S**p

if __name__=='__main__':
    b=g=1.0
    betas=np.array([0.0,0.5,1.0,2.0,5.0,10.0,1e2,1e3,1e4,1e5,1e6,1e7])
    print("="*82)
    print("  beta_eff = beta_* rho : ONE knob for BOTH first-order AND baryon-branching eta")
    print("="*82)
    print(f"  {'beta_eff':>10} {'1st-order?':>11} {'barrier':>9} | {'S_eff':>10} {'eta(p=1)':>10} {'eta(p=2)':>10}")
    fo_win=[]
    for be in betas:
        best=fo.find_Tc_and_barrier(b,g,be)
        is_fo = best is not None and best['barrier']>1e-4
        bar = best['barrier'] if best else 0.0
        if is_fo: fo_win.append(be)
        S=sin4/(phi6*(1+be))
        e1,e2=eta_of(be,1),eta_of(be,2)
        tag='FIRST-ORDER' if is_fo else 'crossover'
        print(f"  {be:>10.0f} {tag:>11} {bar:>9.4f} | {S:>10.2e} {e1:>10.2e} {e2:>10.2e}")

    # eta_obs crossing for each p
    from scipy.optimize import brentq
    be1=brentq(lambda x: eta_of(x,1)-eta_obs, 1e-3, 1e12)
    be2=brentq(lambda x: eta_of(x,2)-eta_obs, 1e-3, 1e12)
    print(f"\n  eta = eta_obs at:  beta_eff = {be1:.1f} (p=1)   beta_eff = {be2:.1f} (p=2)")
    if fo_win:
        print(f"  first-order WINDOW (toy b=g=1): beta_eff in [{min(fo_win):.1f}, {max(fo_win):.0f}]")
        in1 = min(fo_win)<=be1<=max(fo_win); in2 = min(fo_win)<=be2<=max(fo_win)
        print(f"  eta_obs point inside the first-order window?  p=1: {in1}   p=2: {in2}")

    print("\n  READ:")
    print("  * ONE variable beta_eff=beta_* rho does BOTH jobs. The heavy rho-enhancement at QCD")
    print("    (user point) is NOT a problem to dodge -- it is the SOURCE of the baryon suppression")
    print("    (large 1+beta_eff -> small S_eff) that cures the (tt) over-production.")
    print("  * It cannot be TOO heavy: beta_eff > ~1e3-1e4 OVER-saturates -> crossover (hh), so")
    print("    first-order itself CAPS beta_eff. Toy first-order window ~ [2, 1000].")
    print("  * CONSISTENCY is p-dependent (honest):")
    print(f"    - p=2: eta_obs at beta_eff~{be2:.0f} -- INSIDE the window; the window BRACKETS eta_obs")
    print(f"      (over at low beta_eff, under at high). A self-consistent beta_eff(T_QCD) EXISTS.")
    print(f"    - p=1 (the naive power, 1 heavy object): eta_obs at beta_eff~{be1:.0f} -- OUTSIDE/at the")
    print(f"      over-saturated EDGE; within the window p=1 OVER-produces by ~10-3000x. Mild TENSION.")
    print("  * So: a self-consistent solution EXISTS (p=2) but is NOT forced (p=1 mildly over-produces).")
    print("    NOT a derivation: window LOCATION in physical units needs gamma/b from real E_fb (still")
    print("    open, (hh)); p unpinned; beta_eff(T_QCD) co-constrained, not independently computed. [O]")
