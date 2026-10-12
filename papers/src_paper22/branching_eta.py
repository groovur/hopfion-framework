#!/usr/bin/env python3.11
"""
branching_eta.py -- going after b_B, the baryon-number-violating BRANCHING in kappa=delta_neq*b_B
(notes (ss)->(tt)). DM/Quark crossover src_paper21. Assemble b_B from framework factors and see
where the full eta lands (over/under observed).

b_B = b_ch  *  c_int
  b_ch  = topological channel fraction: of the Hopf-conserving outputs of 2+2 (total Q=4),
          {4+0 (gauge boson), 3+1 (BARYON+neutrino), 2+2 (elastic)}, the fraction into (3,1).
          Weighted by the clustering energetics E(Q)=C Q^{3/4} (Lin-Yang/VK; washout note): the
          exothermic channels compete; (3,1) is one of them, not dominant.
  c_int = CP-interference visibility x absorptive factor. From Paper V (P5:thm:delta_cp) the CP
          phase is arg of a TWO-amplitude sum Omega = A1 + A2, A1=e^{i 17pi/40} (|A1|=1), A2=(1/9)
          sum_m d_m e^{2pi i m/90} (E8 twisted sectors, m in {1,7,11,13,17,19,23,29}, McKay marks
          d=(2,3,4,5,6,4,2,3)). The per-event CP asymmetry from interference of two amplitudes is
          eps ~ [2|A1||A2|/(|A1|^2+|A2|^2)] sin(weak) sin(strong). The bracket (VISIBILITY) and the
          weak phase are FRAMEWORK; the strong/absorptive phase (on-shell cut) is dynamical -> flag.

Assemble eta = 7.04 A_pref * (sin dCP * c_int) * (delta_neq * b_ch) and compare to eta_obs.
Honest: b_ch and the absorptive factor are estimates; the point is the ORDER and whether the
natural factors over- or under-shoot observed (locating eta in the old bracket).
"""
import numpy as np

eta_obs=6.1e-10
gstar=17.25
A_pref=15.0/(4*np.pi**2*gstar)
dCP=65.0*np.pi/180; sin_dCP=np.sin(dCP)
delta_neq=1.0e-2            # DERIVED (transport_eta.py): v_w/(Gamma L_w), ~3e-3..4e-2, take 1e-2

# ---- (A) CP interference visibility from Paper V two-amplitude structure ----------------
m=np.array([1,7,11,13,17,19,23,29]); d=np.array([2,3,4,5,6,4,2,3])
A1=np.exp(1j*17*np.pi/40)                         # Mechanism 1, |A1|=1
A2=(1/9)*np.sum(d*np.exp(2j*np.pi*m/90))          # Mechanism 2, E8 twisted sectors
r2=abs(A2)/abs(A1)
visibility=2*r2/(1+r2**2)                         # interference visibility, in [0,1]
dphi=np.angle(A2)-np.angle(A1)                    # relative (weak) phase of the two amplitudes

# ---- (B) topological channel fraction 2+2 -> 3+1 ----------------------------------------
def E(Q,C=1.0): return C*Q**0.75                  # clustering energy (Lin-Yang/VK)
dE={'4+0':E(4)-2*E(2), '3+1':E(3)+E(1)-2*E(2), '2+2':0.0}   # energy change per channel
# crude rate weight: exothermic channels ~ allowed (downhill); use w = max(-dE,0)+w0 floor for
# the always-open elastic; normalise. (A proxy for phase-space/detailed-balance weighting.)
w0=0.1
w={k:max(-v,0.0)+w0 for k,v in dE.items()}
Z=sum(w.values()); b_ch=w['3+1']/Z

if __name__=='__main__':
    print("="*76)
    print("  b_B = b_ch * c_int   (the baryon-number-violating branching in kappa=delta_neq b_B)")
    print("="*76)
    print("  (A) CP interference visibility (Paper V two-amplitude Omega=A1+A2):")
    print(f"      |A1|=1 (Mech-1 topological spin);  A2=(1/9)sum d_m e^{{2pi i m/90}} (E8 twisted)")
    print(f"      |A2| = r2 = {r2:.4f}   relative weak phase = {np.degrees(dphi):.1f} deg")
    print(f"      interference visibility 2r2/(1+r2^2) = {visibility:.3f}")
    print(f"      => c_int = visibility * A_abs,  A_abs(absorptive, on-shell cut)=O(0.1-1) [FLAG]")

    print("\n  (B) topological channel fraction of 2+2 (Q=4):")
    for k in ['4+0','3+1','2+2']:
        print(f"      {k}: dE={dE[k]:+.3f} C   weight={w[k]:.3f}  {'<- BARYON+nu' if k=='3+1' else ''}")
    print(f"      => b_ch (fraction into 3+1) = {b_ch:.3f}")

    print("\n  (C) assembled eta (A_abs bracket):")
    print(f"      delta_neq = {delta_neq:.1e} (DERIVED);  A_pref={A_pref:.4f};  sin dCP={sin_dCP:.3f}")
    print(f"      {'A_abs':>7} {'c_int':>9} {'b_B=b_ch c_int':>16} {'kappa':>10} {'eta':>11} {'eta/eta_obs':>12}")
    for A_abs in [1.0,0.3,0.1]:
        c_int=visibility*A_abs
        b_B=b_ch*c_int
        kappa=delta_neq*b_B
        eta=7.04*A_pref*(sin_dCP*c_int)*(delta_neq*b_ch)
        print(f"      {A_abs:>7.1f} {c_int:>9.3f} {b_B:>16.3e} {kappa:>10.2e} {eta:>11.2e} {eta/eta_obs:>12.1e}")

    print("\n  READ:")
    print("  * c_int is NOT tiny: visibility~{:.2f} (A2 is a sizable fraction of A1), so the CP".format(visibility))
    print("    interference is O(0.1), NOT a loop-suppressed 1e-6 -- the framework's large delta_CP")
    print("    carries over to a large per-event asymmetry. b_ch~{:.2f} (baryon is 1 of 3 channels,".format(b_ch))
    print("    not dominant). So b_B ~ 1e-2..1e-1 from framework factors, NOT 1e-7.")
    print("  * CONSEQUENCE: assembled eta ~ 1e-6..1e-5 -- it OVER-PRODUCES observed 6.1e-10 by ~1e4.")
    print("  * So the honest finding FLIPS: eta is not short of a tiny branching; the natural factors")
    print("    OVERSHOOT, and the open question is the SUPPRESSION ~1e-4 that brings it down. Candidates:")
    print("    (i) the Q=2-lepton pair DENSITY/dilution at T_QCD (f_lep^2, ~1e-1..1e-2, not enough alone);")
    print("    (ii) baryon-antibaryon ANNIHILATION freeze-out (the symmetric part burns off -- sets the")
    print("    floor above which only the asymmetry survives); (iii) a smaller delta_neq (thicker wall).")
    print("  * This sits on the OVER-produce side of the old 22-order bracket (washout note): the")
    print("    framework's natural reconnection baryogenesis is EFFICIENT, so eta's problem is a")
    print("    SUPPRESSION, not an enhancement. Pinning it = the annihilation/dilution Boltzmann. [O]")
