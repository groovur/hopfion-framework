#!/usr/bin/env python3.11
"""
feedback_saddle_core.py -- reusable density-feedback E_fb functional + knot IC +
diagnostics, for the SADDLE (self-consistency) solver shared by the QUARK (trefoil
T(2,3)) and DARK (dihedral) sectors. DM/Quark crossover (src_paper21).

WHY THIS (corrected design, notes (x)+(y)):
  * E_fb = K_fb * J4 is a PRODUCT -> NO Derrick minimum -> the Hopfion is a SADDLE
    (Paper I; qh3_feedback_saddle.py). So MINIMIZATION (BFGS/Adam descent) collapses;
    the object is located by the self-consistency SCORE and a stationary-point search.
    (The bare-Faddeev BFGS/Lin-Yang minimum is for the SUM form, not this product.)
  * The "tube-resolution wall" was UNDER-RESOLUTION: the old run used N=64 h=0.175
    (~2 pts/tube-radius); the Faddeev-Hopfion literature standard is N~200-240 h~0.08
    (~5 pts/tube). High-res is HPC/GPU, not a fundamental barrier.
This module provides the REUSABLE pieces (IC, E_fb, diagnostics, Q_H) validated at low
res; the high-res saddle search is the compute-gated next layer (torch autograd descent
+ score-snapshot, or a matrix-free stationary-point solve).

Functional (ported from qh3_feedback_saddle.py; n = unit director field on a grid):
  g2 = |grad n|^2 ;  s4 = (1-nz^2)^2
  J2a = int s4*g2 dV ;  J2iso_fb = int g2/(1+beta*g2) dV ;  K_fb = J2a + mu*J2iso_fb,
  mu = 3 - phi ;  J4 = int (Fxy^2+Fxz^2+Fyz^2) dV, F_ij = n.(d_i n x d_j n)
  E_fb = K_fb * J4 ;  V = J2iso_fb/J2a (target phi) ;  s_opt = (phi^6 * J4/K_fb)^(1/2)
  score = |V - phi| + |s_opt - 1|   (self-consistency; small+stable = the saddle)
Q_H via the VALIDATED Whitehead integral (src_paper16/whitehead_hopf_charge.py).
"""
import numpy as np
from scipy.spatial import KDTree

phi = (1 + 5**0.5)/2
mu_fb = 3 - phi

# ---------- knot-centerline IC (reusable: (p,q) torus knot; trefoil = (2,3)) ----------
def torus_knot_centerline(p, q, R0, r0, NT=8000):
    t = np.linspace(0, 2*np.pi, NT, endpoint=False)
    Gx = (R0 + r0*np.cos(q*t))*np.cos(p*t)
    Gy = (R0 + r0*np.cos(q*t))*np.sin(p*t)
    Gz = r0*np.sin(q*t)
    return t, np.stack([Gx, Gy, Gz], axis=1)

def build_doublet(pts, centerline, t_arr, q, C_star):
    """S3-lift Hopf doublet Z=(z1,z2) around a knot centerline (whitehead_hopf_charge.py
    construction, validated to give Q_H -> p for the (p,q) knot at increasing N)."""
    tree = KDTree(centerline)
    def f0(rh): return 2*np.arctan(np.maximum(rh, 1e-9)**(-C_star))
    dists, idxs = tree.query(pts, k=80)
    t1 = t_arr[idxs[:,0]]; d1 = dists[:,0]
    dt = np.abs(((t_arr[idxs] - t1[:,None] + np.pi) % (2*np.pi)) - np.pi)
    mask = dt > np.pi/2
    fi = np.argmax(mask, axis=1); has2 = mask.any(axis=1)
    d2 = np.where(has2, dists[np.arange(len(pts)), fi], 1e6)
    t2 = t_arr[idxs[np.arange(len(pts)), fi]]
    f1 = f0(np.clip(d1,1e-6,None)*C_star); f2 = f0(np.clip(d2,1e-6,None)*C_star)
    w1 = 1.0/np.clip(d1,1e-6,None)**2; w2 = 1.0/np.clip(d2,1e-6,None)**2
    z1 = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2 = w1*np.sin(f1/2)*np.exp(1j*q*t1) + w2*np.sin(f2/2)*np.exp(1j*q*t2)
    m = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2) + 1e-12
    return z1/m, z2/m

def nhat_from_doublet(z1, z2):
    n1 = 2*np.real(np.conj(z1)*z2); n2 = 2*np.imag(np.conj(z1)*z2)
    n3 = np.abs(z1)**2 - np.abs(z2)**2
    return np.stack([n1, n2, n3], axis=0)

# ---------- diagnostics ----------
def whitehead_QH(z1, z2, h):
    dz1 = [np.gradient(z1, h, axis=i) for i in range(3)]
    dz2 = [np.gradient(z2, h, axis=i) for i in range(3)]
    A = [np.imag(np.conj(z1)*dz1[i] + np.conj(z2)*dz2[i]) for i in range(3)]
    curl = [np.gradient(A[2],h,axis=1)-np.gradient(A[1],h,axis=2),
            np.gradient(A[0],h,axis=2)-np.gradient(A[2],h,axis=0),
            np.gradient(A[1],h,axis=0)-np.gradient(A[0],h,axis=1)]
    return h**3*np.sum(sum(A[i]*curl[i] for i in range(3)))/(4*np.pi**2)

def Efb_diagnostics(n, h, beta):
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    g2 = sum((dn[ax]**2).sum(0) for ax in range(3))            # |grad n|^2
    s4 = (1 - n[2]**2)**2
    def cross(u, v): return np.stack([u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]])
    F2 = np.zeros_like(g2)
    for i in range(3):
        for j in range(i+1,3):
            F2 = F2 + (n*cross(dn[i], dn[j])).sum(0)**2
    dV = h**3
    J2a = dV*np.sum(s4*g2)
    J2iso_fb = dV*np.sum(g2/(1 + beta*g2))
    J4 = dV*np.sum(F2)
    K_fb = J2a + mu_fb*J2iso_fb
    E_fb = K_fb*J4
    V = J2iso_fb/J2a if J2a > 0 else np.nan
    s_opt = np.sqrt(phi**6 * J4/K_fb) if K_fb > 0 else np.nan
    score = abs(V - phi) + abs(s_opt - 1)
    return dict(J2a=J2a, J2iso_fb=J2iso_fb, J4=J4, K_fb=K_fb, E_fb=E_fb,
                V=V, s_opt=s_opt, s_opt6=s_opt**6, score=score, J4_over_J2a=J4/J2a)

# ---------- smoke / validation ----------
if __name__ == '__main__':
    import argparse, time
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=48)
    ap.add_argument('--L', type=float, default=6.0)
    ap.add_argument('--p', type=int, default=2); ap.add_argument('--q', type=int, default=3)
    ap.add_argument('--R0', type=float, default=3.0); ap.add_argument('--r0', type=float, default=0.874)
    ap.add_argument('--C_star', type=float, default=2.5062); ap.add_argument('--beta', type=float, default=0.452)
    a = ap.parse_args()
    ax = np.linspace(-a.L, a.L, a.N); h = ax[1]-ax[0]
    X,Y,Z = np.meshgrid(ax,ax,ax, indexing='ij')
    pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()], axis=1)
    print("="*68)
    print(f"  feedback_saddle_core SMOKE  knot T({a.p},{a.q})  N={a.N} L={a.L} h={h:.3f} beta={a.beta}")
    print(f"  tube radius ~1/C* = {1/a.C_star:.3f} ; pts/tube = {(1/a.C_star)/h:.1f}  (lit wants ~5)")
    print("="*68)
    t,cl = torus_knot_centerline(a.p, a.q, a.R0, a.r0)
    t0=time.time(); z1,z2 = build_doublet(pts, cl, t, a.q, a.C_star)
    z1=z1.reshape(a.N,a.N,a.N); z2=z2.reshape(a.N,a.N,a.N)
    QH = whitehead_QH(z1,z2,h)
    n = nhat_from_doublet(z1,z2)
    d = Efb_diagnostics(n,h,a.beta)
    print(f"  Q_H = {QH:.3f}   (expect ~{a.p} with grid underestimate; trefoil p=2 -> ~1.5-2 at N=48)")
    print(f"  E_fb={d['E_fb']:.3e}  K_fb={d['K_fb']:.3e}  J4={d['J4']:.3e}  J4/J2a={d['J4_over_J2a']:.3f}")
    print(f"  V={d['V']:.3f} (target phi={phi:.3f})  s_opt={d['s_opt']:.3f} (target 1)  s_opt^6={d['s_opt6']:.3f}")
    print(f"  self-consistency score = {d['score']:.3f}   [{time.time()-t0:.1f}s]")
    print("  (SMOKE validates IC+functional+Q_H machinery; HIGH-RES saddle search is the next, compute-gated layer)")
