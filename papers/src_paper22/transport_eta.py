#!/usr/bin/env python3.11
"""
transport_eta.py -- ANATOMY of the single open quantity: the transport efficiency kappa in
eta = 7.04 A_pref eps_CP kappa (wall_eta.py, notes (rr)->(ss)). DM/Quark crossover src_paper22.

Goal: NOT to re-derive the Gamma/H bracket (washout_eta.py: Gamma/H~3e19, deep equilibrium) but
to DECOMPOSE kappa into derived vs open factors, so the open part is a single physical number.

PHYSICS. Gamma/H ~ 3e19 = DEEP equilibrium. By CPT, n_B = 0 in EXACT equilibrium regardless of
C,CP. So the surviving asymmetry is set by the DEPARTURE from equilibrium at the bubble wall,
then by the fraction of reconnections that actually change baryon number:
  kappa = delta_neq * b_B
  delta_neq = (wall perturbation rate)/(re-equilibration rate) = (v_w/L_w)/Gamma   [DERIVED]
      -- the system equilibrates at Gamma; the wall sweeps a point in tau_wall=L_w/v_w; since
         Gamma*tau_wall>>1 the departure is small, delta_neq = 1/(Gamma tau_wall) = v_w/(Gamma L_w).
  b_B = branching: fraction of CP-biased reconnections that are baryon-number-VIOLATING
      (Q=2+2->3+1 threadings) among ALL condensate reconnections, x the net CP-interference
      efficiency c_int.                                                               [OPEN]
INPUTS (derived upstream): Gamma/H (washout_eta), v_w~0.15-0.2 (wall_eta, from alpha~0.05),
delta_CP=65 deg (Paper V), A_pref=15/(4pi^2 g*). L_w = c_w/T, c_w in [2,20] = DIMENSIONLESS wall
thickness, set by the cutoff-dependent stiffness of the L field (claude-hopfion s10: cutoff-set,
NOT derived -> flagged). Output: delta_neq (derived, bounds kappa from above), and the OPEN
branching b_B that observed eta requires. Verdict: eta's magnitude reduces to ONE branching ratio.
"""
import numpy as np

GeV=1.0; MPl=1.22e19; eta_obs=6.1e-10
Tq=0.15                      # T_QCD ~ Lcond(QCD)
gstar=17.25
A_pref=15.0/(4*np.pi**2*gstar)      # ~0.022
sin_dCP=np.sin(65.0*np.pi/180)      # 0.906

def H(T): return 1.66*np.sqrt(gstar)*T**2/MPl
def Gamma(T,Lam=0.2,gdof=40): return (0.122*gdof*T**3)*(1.0/Lam**2)   # n sigma v, v=1

if __name__=='__main__':
    G=Gamma(Tq); h=H(Tq)
    print("="*76)
    print("  ANATOMY of the transport efficiency kappa   (eta = 7.04 A_pref eps_CP kappa)")
    print("="*76)
    print(f"  T_QCD={Tq} GeV   Gamma={G:.3f} GeV   H={h:.2e} GeV   Gamma/H={G/h:.2e}  (deep eq.)")
    print(f"  A_pref=15/(4pi^2 g*)={A_pref:.4f}   sin(delta_CP)={sin_dCP:.3f}   eta_obs={eta_obs:.1e}")
    print(f"  CPT: n_B=0 in EXACT equilibrium -> only the DEPARTURE at the wall survives.\n")

    print("  (1) DERIVED departure-from-equilibrium delta_neq = v_w/(Gamma L_w), L_w=c_w/T:")
    print(f"      {'v_w':>5} {'c_w(=L_w*T)':>12} {'tau_wall*Gamma':>15} {'delta_neq':>12}")
    rows=[]
    for v_w in [0.15,0.20]:
        for c_w in [2,5,10,20]:
            Lw=c_w/Tq                      # wall thickness (GeV^-1)
            tau_wall=Lw/v_w                # time to sweep a point (GeV^-1)
            Gtau=G*tau_wall
            delta_neq=1.0/Gtau             # = v_w/(Gamma L_w)
            rows.append((v_w,c_w,Gtau,delta_neq))
            print(f"      {v_w:>5.2f} {c_w:>12} {Gtau:>15.1f} {delta_neq:>12.3e}")
    dn=[r[3] for r in rows]
    print(f"      => delta_neq in [{min(dn):.1e}, {max(dn):.1e}]  (DERIVED; this BOUNDS kappa above)")

    print("\n  (2) OPEN branching b_B that observed eta requires  (kappa = delta_neq * b_B):")
    print("      eta = 7.04 A_pref (sin dCP * c_int) (delta_neq * b_B)")
    for v_w,c_w,Gtau,delta_neq in [rows[1],rows[2],rows[5],rows[6]]:
        kappa_need=eta_obs/(7.04*A_pref*sin_dCP)
        bB_cint=kappa_need/delta_neq
        print(f"      v_w={v_w}, c_w={c_w}: delta_neq={delta_neq:.2e}  => needed b_B*c_int = {bB_cint:.2e}")

    print("\n  (3) VERDICT:")
    print("  * kappa is NOT a free fudge: it factorises into a DERIVED piece and ONE open piece.")
    print("    - delta_neq ~ 1e-2..1e-3 (DERIVED from v_w, Gamma, and the wall thickness) -- the")
    print("      departure from equilibrium; caps kappa <~ 1e-2, so eta cannot be O(1) (good:")
    print("      resolves the 'no-washout overproduces' edge of the old bracket WITHOUT tuning).")
    print("    - b_B*c_int ~ 1e-6..1e-7 = the fraction of reconnections that VIOLATE baryon number")
    print("      (Q=2+2->3+1 threadings among all condensate reconnections) x CP interference.")
    print("  * So eta's magnitude reduces to ONE physical ratio: the baryon-number-violating")
    print("    BRANCHING ~1e-6..1e-7. Plausible (B-violating threadings ARE rare vs generic")
    print("    reconnections) but NOT derived -- it needs the reconnection Boltzmann + the")
    print("    Q-selection rule on the 2I vertex. The wall thickness c_w is cutoff-set (flag).")
    print("  * The old 22-order bracket is now a ~1-order window on a BRANCHING RATIO: the")
    print("    exponential washout sensitivity is GONE because confinement caps the washout window")
    print("    (rate shuts off behind the wall), leaving a power-law (delta_neq), not exp(-Gamma tau).")
