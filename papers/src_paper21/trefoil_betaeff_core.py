#!/usr/bin/env python3.11
"""
trefoil_betaeff_core.py -- SMOKE TEST of the collective-coordinate route to the eta magnitude
(Paper XXII beta_*rho) and a prereq for the existence search. DM/Quark crossover.

IDEA (notes, 2026-10-13): the eta window's beta_*rho is the trefoil's self-consistent feedback at
the CORE. Rather than a full-field high-res relaxation (the resolution wall), EVALUATE the knotted
T(2,3) ansatz (already in feedback_saddle_core: torus_knot_centerline + build_doublet) and extract:
  * the self-consistent feedback coupling beta: solve V(beta) = J2iso_fb(beta)/J2a = phi
    (the same J_fb/J_a=phi self-consistency that gives beta_*=0.452 for the lepton);
  * the CORE-PEAK gradient energy g2_peak = max|grad n|^2 (vs the average that sets V);
  * the local core saturation  beta_eff_core = beta * g2_peak  == the eta channel's beta_*rho.
The analytic probe used g2 ~ C*^2 ~ 6.3 (uniform) -> beta_*rho ~ 2.8; the CORE PEAK is higher, and
beta is fixed by self-consistency, so this returns the ACTUAL beta_eff_core -- the factor-~5 residual.

SMOKE GOALS (compute discipline): (a) does Q_H converge to 3? (b) does beta_eff_core converge with N?
(c) time each N -> extrapolate the N~200 cost. NOT a production run; validates the method + cost.
"""
import numpy as np, sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feedback_saddle_core import (torus_knot_centerline, build_doublet,
                                   nhat_from_doublet, whitehead_QH, mu_fb)
from rational_map_ic import rational_map_doublet
from scipy.optimize import brentq

phi = (1+np.sqrt(5))/2

def g2_field(n, h):
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    return sum((dn[ax]**2).sum(0) for ax in range(3))        # |grad n|^2 at every point

def analyse(N, L, R0, r0, C_star, mode='smooth', p=1, q=3):
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X,Y,Z = np.meshgrid(ax,ax,ax, indexing='ij')
    if mode == 'smooth':                          # charge-resolving A_{p,q} base (Q_H=p*q)
        z1,z2 = rational_map_doublet(X,Y,Z,p,q)
    else:                                         # sharp knotted tube (loses Q_H at low N)
        pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
        t,cl = torus_knot_centerline(2,3,R0,r0)
        z1,z2 = build_doublet(pts, cl, t, 3, C_star)
        z1=z1.reshape(X.shape); z2=z2.reshape(X.shape)
    QH = whitehead_QH(z1,z2,h)
    n = nhat_from_doublet(z1,z2)
    g2 = g2_field(n,h); dV=h**3
    s4 = (1-n[2]**2)**2
    J2a = dV*np.sum(s4*g2)
    def V_of(beta): return (dV*np.sum(g2/(1+beta*g2)))/J2a
    # self-consistency V(beta)=phi ; V(0)=V_max decreasing in beta
    Vmax = V_of(0.0)
    if Vmax > phi:
        beta_sc = brentq(lambda b: V_of(b)-phi, 1e-6, 1e6)
    else:
        beta_sc = np.nan                                      # ansatz already below phi (coarse)
    g2_peak = g2.max()
    # energy-weighted mean g2 (the 'average' that sets the bulk of V)
    w = g2/g2.sum(); g2_ew = (w*g2).sum()
    be_core = beta_sc*g2_peak
    be_ew   = beta_sc*g2_ew
    return dict(N=N,h=h,QH=QH,Vmax=Vmax,beta_sc=beta_sc,g2_peak=g2_peak,
                g2_ew=g2_ew,be_core=be_core,be_ew=be_ew)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--Nlist',type=str,default='48,64,96')
    ap.add_argument('--L',type=float,default=5.0)
    ap.add_argument('--R0',type=float,default=3.0)
    ap.add_argument('--r0',type=float,default=0.874)
    ap.add_argument('--C_star',type=float,default=2.5062)
    ap.add_argument('--mode',type=str,default='smooth',choices=['smooth','sharp'])
    ap.add_argument('--p',type=int,default=1); ap.add_argument('--q',type=int,default=3)
    a=ap.parse_args()
    print("="*96)
    print(f"  TREFOIL core beta_eff SMOKE  T(2,3) R0={a.R0} r0={a.r0} C*={a.C_star}  box L={a.L}")
    print(f"  self-consistent beta (V=phi={phi:.4f}); beta_eff_core=beta*g2_peak == eta channel beta_*rho")
    print("="*96)
    print(f"  {'N':>4} {'h':>6} {'Q_H':>7} {'Vmax':>6} {'beta_sc':>9} {'g2_peak':>9} {'g2_ew':>8} "
          f"{'be_core':>8} {'be_ew':>7} {'t(s)':>7}")
    times=[]
    for N in [int(s) for s in a.Nlist.split(',')]:
        t0=time.time(); d=analyse(N,a.L,a.R0,a.r0,a.C_star,a.mode,a.p,a.q); dt=time.time()-t0; times.append((N,dt))
        print(f"  {d['N']:>4} {d['h']:>6.3f} {d['QH']:>7.3f} {d['Vmax']:>6.3f} {d['beta_sc']:>9.4f} "
              f"{d['g2_peak']:>9.2f} {d['g2_ew']:>8.3f} {d['be_core']:>8.2f} {d['be_ew']:>7.3f} {dt:>7.1f}")
    print(f"\n  eta needs beta_*rho ~ 15 (p=2). READ: does be_core approach ~15 as N grows? Q_H -> 3?")
    if len(times)>=2:
        (N1,t1),(N2,t2)=times[-2],times[-1]
        scal=np.log(t2/t1)/np.log(N2/N1) if t2>t1 else 3.0
        print(f"  timing scales ~ N^{scal:.1f}; extrapolated N=200 single-eval: ~{t2*(200/N2)**scal:.0f} s")
    print("="*96)
