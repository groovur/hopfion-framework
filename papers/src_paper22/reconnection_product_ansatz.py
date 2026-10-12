#!/usr/bin/env python3.11
"""
reconnection_product_ansatz.py -- step (d), first milestone: does the baryon-formation
channel Q_H=2 + Q_H=1 -> Q_H=3 exist, and what is its interaction energy E(R)?
DM/Quark crossover (src_paper21). Dual-purpose: the SAME reconnection rate feeds the
quark baryogenesis eta AND the dark co-genesis (notes (t),(u),(aa)).

METHOD (product ansatz, Battye-Sutcliffe; cf. arXiv:1503.08027):
  Two SMOOTH rational-map Hopfions (rational_map_ic.py, VALIDATED Q_H-convergent) at
  separation R on the x-axis. Vacuum is u=0 (north pole), soliton core u=infty, so the
  ADDITIVE stereographic ansatz u_tot = u_1 + u_2 reduces to soliton_i near c_i and
  vacuum far away -> combines Hopf charges for well-separated solitons.
  Diagnostics vs separation R: Q_H (does it -> Q1+Q2=3?) and the FADDEEV (sum) energy
  E_F = int|grad n|^2 + lam * J4 (well-defined scale, unlike the product E_fb which is
  scale-free) -> the interaction potential: attractive (favours fusion) or a barrier.

SCOPE (honest): semiclassical, product-ansatz (not relaxed) -> E(R) is the collective-
coordinate effective potential, order-of-magnitude. First milestone = channel existence
+ barrier SCALE (the input to eta's efficiency), NOT a precise rate. Relative orientation
fixed (aligned) here; the orientation/phase scan (attractive channel) is the next refinement.
"""
import numpy as np, sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feedback_saddle_core import whitehead_QH, nhat_from_doublet

def u_shifted(X, Y, Z, p, q, c):
    x,y,z = X-c[0], Y-c[1], Z-c[2]; r2 = x*x+y*y+z*z
    a = (2*(x+1j*y))/(1+r2); b = ((1-r2)+2j*z)/(1+r2)
    return (a**p)/(b**q + 1e-30)

def faddeev_energy(n, h, lam):
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    g2 = sum((dn[ax]**2).sum(0) for ax in range(3))
    def cross(u,v): return np.stack([u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]])
    F2 = np.zeros_like(g2)
    for i in range(3):
        for j in range(i+1,3):
            F2 = F2 + (n*cross(dn[i],dn[j])).sum(0)**2
    dV = h**3
    E2 = dV*g2.sum(); E4 = dV*F2.sum()
    return E2 + lam*E4, E2, E4

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=80); ap.add_argument('--L', type=float, default=5.0)
    ap.add_argument('--lam', type=float, default=1.0)
    ap.add_argument('--p1', type=int, default=1); ap.add_argument('--q1', type=int, default=2)  # Q_H=2
    ap.add_argument('--p2', type=int, default=1); ap.add_argument('--q2', type=int, default=1)  # Q_H=1
    ap.add_argument('--Rlist', type=str, default='6,5,4,3,2.5,2,1.5,1,0.5,0.0')
    a = ap.parse_args()
    ax = np.linspace(-a.L, a.L, a.N); h = ax[1]-ax[0]
    X,Y,Z = np.meshgrid(ax,ax,ax, indexing='ij')
    print("="*76)
    print(f"  PRODUCT-ANSATZ RECONNECTION  Q_H={a.p1*a.q1}+Q_H={a.p2*a.q2} -> ?   N={a.N} L={a.L} lam={a.lam}")
    print(f"  (expect Q_H -> {a.p1*a.q1 + a.p2*a.q2} at large R if the fusion channel exists)")
    print("="*76)
    def doublet_from_u(u):
        z1 = 1.0/np.sqrt(1+np.abs(u)**2); return z1.astype(complex), (u*z1).astype(complex)
    def su2_product(Z1a,Z2a, Z1b,Z2b):
        """quaternion/SU(2) product of two unit doublets (the correct Hopfion product ansatz):
        combines Hopf charges; far from a soliton its doublet -> (1,0)=identity."""
        Z1 = Z1a*Z1b - np.conj(Z2a)*Z2b
        Z2 = Z2a*Z1b + np.conj(Z1a)*Z2b
        return Z1, Z2
    print(f"  {'R':>5} {'Q_H':>7} {'E_F':>10} {'E2(grad)':>10} {'E4(Skyrme)':>11}")
    rows=[]
    for R in [float(s) for s in a.Rlist.split(',')]:
        u1 = u_shifted(X,Y,Z, a.p1,a.q1, (-R/2,0,0))
        u2 = u_shifted(X,Y,Z, a.p2,a.q2, (+R/2,0,0))
        Z1a,Z2a = doublet_from_u(u1); Z1b,Z2b = doublet_from_u(u2)
        Z1,Z2 = su2_product(Z1a,Z2a, Z1b,Z2b)
        QH = whitehead_QH(Z1, Z2, h)
        n = nhat_from_doublet(Z1, Z2)
        EF, E2, E4 = faddeev_energy(n, h, a.lam)
        rows.append((R,QH,EF)); print(f"  {R:>5.2f} {QH:>7.3f} {EF:>10.2f} {E2:>10.2f} {E4:>11.2f}")
    # interaction potential relative to large-R asymptote
    E_inf = rows[0][2]
    print(f"\n  interaction energy  Delta E_F(R) = E_F(R) - E_F(R={rows[0][0]:.1f}):")
    for R,QH,EF in rows:
        print(f"    R={R:>5.2f}  dE_F={EF-E_inf:>+9.2f}  Q_H={QH:.2f}")
    print("\n  READ: Q_H->3 at large R = channel exists. dE_F<0 monotonic as R->0 = ATTRACTIVE (barrierless")
    print("        fusion); dE_F>0 bump = a BARRIER (its height ~ the reconnection-rate scale for eta).")
