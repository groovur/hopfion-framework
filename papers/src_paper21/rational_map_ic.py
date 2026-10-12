#!/usr/bin/env python3.11
"""
rational_map_ic.py -- SMOOTH rational-map Hopfion initial condition (step (a)).
DM/Quark crossover (src_paper21). Fixes the tube-lift IC failure (feedback_saddle_core
smoke: Q_H=0 because the charge sat in 0.5% of grid points -> effectively single-strand).

CONSTRUCTION (Hietarinta-Salo / Sutcliffe rational-map Hopfion):
  S^3 coords from R^3:  a = 2(x+iy)/(1+r^2),  b = ((1-r^2)+2iz)/(1+r^2),  |a|^2+|b|^2=1.
  Rational map  u = a^p / b^q  ->  doublet Z=(Z1,Z2)=(b^q, a^p)/norm  ->  n_hat = Z^dag sigma Z.
  This is the AXIALLY-SYMMETRIC A_{p,q} Hopfion: core loop = unit circle (b=0), field smooth
  and charge DISTRIBUTED over S^3 (not a thin tube) -> Q_H converges with resolution.
  Hopf charge expected Q = p*q (VALIDATE numerically: converge to the integer as N grows).

GOAL of (a): a validated SMOOTH IC whose Q_H converges to the right integer -- the prereq
for the saddle search / collective-coordinate dynamics (steps d). The genuine KNOTTED
trefoil T(2,3) (Sutcliffe 2013, arXiv:1304.6021) is a MODULATED rational map, the next
refinement on top of this validated smooth base.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feedback_saddle_core import whitehead_QH, nhat_from_doublet, Efb_diagnostics

def rational_map_doublet(X, Y, Z, p, q):
    r2 = X**2 + Y**2 + Z**2
    a = (2*(X + 1j*Y)) / (1 + r2)            # S^3 azimuthal coord
    b = ((1 - r2) + 2j*Z) / (1 + r2)         # S^3 poloidal coord
    Z1 = b**q                                 # u = Z2/Z1 = a^p / b^q
    Z2 = a**p
    nrm = np.sqrt(np.abs(Z1)**2 + np.abs(Z2)**2) + 1e-30
    return Z1/nrm, Z2/nrm

def localization(n, h, X, Y, Z):
    """energy-weighted radius + fraction of grid carrying 50% of |grad n|^2 (tube-lift was ~0.5%)."""
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    g2 = sum((dn[ax]**2).sum(0) for ax in range(3))
    w = g2/g2.sum(); r = np.sqrt(X**2+Y**2+Z**2)
    rbar = (w*r).sum()
    flat = np.sort(g2.ravel())[::-1]; cum = np.cumsum(flat)/flat.sum()
    frac50 = (np.searchsorted(cum, 0.5)+1)/g2.size
    return rbar, frac50

if __name__ == '__main__':
    import argparse, time
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=int, default=1); ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--L', type=float, default=4.0)
    ap.add_argument('--Nlist', type=str, default='40,56,72,96')
    a = ap.parse_args()
    print("="*72)
    print(f"  RATIONAL-MAP Hopfion IC  u = a^{a.p}/b^{a.q}   expected Q_H = p*q = {a.p*a.q}")
    print(f"  (box L={a.L}; smooth A_{{p,q}} Hopfion -- Q_H should CONVERGE to the integer)")
    print("="*72)
    print(f"  {'N':>4} {'h':>7} {'Q_H':>8} {'r_bar':>7} {'frac(50% energy)':>17}  (tube-lift was 0.005)")
    for N in [int(s) for s in a.Nlist.split(',')]:
        ax = np.linspace(-a.L, a.L, N); h = ax[1]-ax[0]
        X,Y,Z = np.meshgrid(ax,ax,ax, indexing='ij')
        t0=time.time()
        Z1,Z2 = rational_map_doublet(X,Y,Z,a.p,a.q)
        QH = whitehead_QH(Z1,Z2,h)
        n = nhat_from_doublet(Z1,Z2)
        rbar, frac50 = localization(n,h,X,Y,Z)
        print(f"  {N:>4} {h:>7.3f} {QH:>8.3f} {rbar:>7.3f} {frac50:>17.3f}  [{time.time()-t0:.1f}s]")
    print("  READ: Q_H -> p*q as N grows (smooth IC, charge distributed) = rational-map IC VALIDATED.")
