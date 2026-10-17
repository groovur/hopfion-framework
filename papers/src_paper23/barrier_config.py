#!/usr/bin/env python3.11
"""
barrier_config.py -- Paper XXIII component B: the pion emission/absorption BARRIER and a reaction
coordinate from the Hopf link (Q_H=2) to the trefoil (Q_H=3). DM/Quark crossover (src_paper23).

GROUNDING: the Brieskorn family {a^2 + b^n = 0} IS the satellite twist family (Paper XVIII
P18:prop:satellite_twist_family): n=2 Hopf link T(2,2) (Q_H=2), n=3 trefoil T(2,3) (Q_H=3), n=4 the
unfused 3-component BARRIER. So a single-valued REACTION COORDINATE lambda in [0,1] is
    u(lambda) = (a^2 + (1-lambda) b^2 + lambda * c3 * b^3) / a,
with lambda=0 -> Hopf link (Q_H=2), lambda=1 -> trefoil (Q_H=3, c3=3 from knotted_map.py). The reaction
2+2 -> 3+1 changes Q_H by +1 ACROSS the barrier; that +1 is the PION (Q_H=1), emitted/absorbed exactly
at the reconnection lambda* where the zero set {numerator=0} changes topology (two components fuse).

This script: scans lambda, reports Q_H(lambda) and V(lambda) via the gate, and LOCATES the barrier
lambda* (where Q_H transitions 2->3 and the energy / self-interaction peaks). The barrier config
n(lambda*) is what component F (betastar_rho_barrier) will evaluate the feedback saturation on --
the user's insight: the saturation sits "just before full absorption / just after emission begins."

CAVEAT: single-valued polynomial interpolation of the Brieskorn exponent; the barrier lambda* is read
from Q_H(lambda) + E(lambda), not assumed. Proportions (c3) inherited from the tuned knotted map.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from knotted_map import s3_coords
from validation import gate_field

def reaction_doublet(X,Y,Z, lam, c3=3.0):
    a,b=s3_coords(X,Y,Z)
    Z1=a
    Z2=a**2 + (1-lam)*(b**2) + lam*c3*(b**3)
    nrm=np.sqrt(np.abs(Z1)**2+np.abs(Z2)**2)+1e-30
    return Z1/nrm, Z2/nrm

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--L',type=float,default=5.0)
    ap.add_argument('--N',type=int,default=128); ap.add_argument('--c3',type=float,default=3.0)
    a=ap.parse_args()
    ax=np.linspace(-a.L,a.L,a.N); h=ax[1]-ax[0]
    X,Y,Z=np.meshgrid(ax,ax,ax,indexing='ij')
    print("="*76)
    print(f"  REACTION COORDINATE Hopf link (Q_H=2) -> trefoil (Q_H=3)  [L={a.L} N={a.N} c3={a.c3}]")
    print(f"  u(lam)=(a^2 + (1-lam)b^2 + lam*{a.c3}*b^3)/a ;  barrier = where Q_H jumps + E peaks")
    print("="*76)
    print(f"  {'lambda':>7} {'Q_H':>7} {'V(b=0)':>8} {'E_fb(b=.452)':>13} {'J4/J2a':>8}")
    rows=[]
    for lam in np.linspace(0,1,11):
        z1,z2=reaction_doublet(X,Y,Z,lam,a.c3)
        g0=gate_field(z1,z2,h,beta=0.0); g=gate_field(z1,z2,h,beta=0.452)
        rows.append((lam,g['QH'],g0['V'],g['E_fb'],g['J4_over_J2a']))
        print(f"  {lam:>7.2f} {g['QH']:>7.3f} {g0['V']:>8.3f} {g['E_fb']:>13.3e} {g['J4_over_J2a']:>8.2f}")
    # locate barrier: steepest dQ_H/dlam (the reconnection) and/or E_fb maximum
    QH=np.array([r[1] for r in rows]); E=np.array([r[3] for r in rows]); lams=np.array([r[0] for r in rows])
    dQ=np.gradient(QH,lams); i_recon=int(np.argmax(np.abs(dQ))); i_emax=int(np.argmax(E))
    print()
    print(f"  Q_H: {QH[0]:.2f} (Hopf link) -> {QH[-1]:.2f} (trefoil);  net change ~ {QH[-1]-QH[0]:+.2f} (= pion unit)")
    print(f"  steepest Q_H change (reconnection) at lambda ~ {lams[i_recon]:.2f}  (dQ_H/dlam = {dQ[i_recon]:.2f})")
    print(f"  E_fb maximum (self-interaction peak) at lambda ~ {lams[i_emax]:.2f}")
    print(f"  => BARRIER lambda* ~ {0.5*(lams[i_recon]+lams[i_emax]):.2f}  (where the feedback saturation lives;")
    print(f"     'just before absorption / after emission'). Component F evaluates beta_*rho here.")
    print("="*76)
