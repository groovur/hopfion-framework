#!/usr/bin/env python3.11
"""
tube_q2_compare.py -- validate the lepton(Q=2)/quark(Q=3) split (a & b).

(a) FEEDBACK SATURATION: does the sharp quark tube drive the density-feedback
    term into saturation (beta*g2>>1, no restoring stiffness) more than the
    smoother lepton torus?
(b) SELF-CONSISTENCY REACHABILITY: V=J2iso/J2a is scale-INVARIANT and bounded
    above by its beta=0 value V_max=J2iso_bare/J2a (feedback only lowers V). So
    a self-consistent config (V=phi) can exist ONLY IF V_max >= phi. This is a
    clean, beta-independent, shape-only diagnostic. Compare V_max for Q=2 vs Q=3.

Q=3 = trefoil (2,3) torus knot, blend+partition (tube_trefoil).
Q=2 = unknot circle torus, Hopf winding q=2 (the smooth, crossing-free lepton-
      like config), reusing the Phase-2 circular machinery (tube_curvilinear).
Both scanned over tube sharpness C* (radius 1/C*) to separate TOPOLOGY (matched
C*) from FATNESS (varying C*).

Run:  python3.11 tube_q2_compare.py
"""
import numpy as np
from tube_curvilinear import (MU, PHI, circular_frame, circular_field_grid,
                              curvilinear_energy)
from tube_trefoil import BlendField, trefoil_frame, partition_weight

LAM = PHI**6


def shape_diagnostics(Gamma, T, N1, N2, k1, k2, s_grid, rho_grid, psi_grid,
                      build_n, node_weight=None):
    """J2a, J2iso_bare, J4 and the per-node (g2, dV) for saturation analysis."""
    ds = s_grid[1]-s_grid[0]; drho = rho_grid[1]-rho_grid[0]; dpsi = psi_grid[1]-psi_grid[0]
    S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
    cosp, sinp = np.cos(PSI), np.sin(PSI)
    er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
    X = Gamma[:, None, None, :] + RHO[..., None]*er
    n = build_n(X, S, RHO, PSI); n /= np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-12)
    hs = 1.0 - RHO*(k1[:, None, None]*cosp + k2[:, None, None]*sinp)
    Dn_s = (np.roll(n,-1,0)-np.roll(n,1,0))/(2*ds)/hs[...,None]
    Dn_p = (np.roll(n,-1,2)-np.roll(n,1,2))/(2*dpsi)/RHO[...,None]
    inter=(n[:,2:]-n[:,:-2])/(2*drho); lo=(n[:,1:2]-n[:,0:1])/drho; hi=(n[:,-1:]-n[:,-2:-1])/drho
    Dn_r = np.concatenate([lo, inter, hi], axis=1)
    g2 = (Dn_s**2).sum(-1)+(Dn_r**2).sum(-1)+(Dn_p**2).sum(-1)
    s4 = np.clip(1-n[...,2]**2, 0, 1)**2
    dV = hs*RHO*ds*drho*dpsi
    if node_weight is not None:
        dV = dV*node_weight
    Fsr=(n*np.cross(Dn_s,Dn_r)).sum(-1); Fsp=(n*np.cross(Dn_s,Dn_p)).sum(-1); Frp=(n*np.cross(Dn_r,Dn_p)).sum(-1)
    J4 = ((Fsr**2+Fsp**2+Frp**2)*dV).sum()
    J2a = (s4*g2*dV).sum(); J2iso = (g2*dV).sum()
    return dict(J2a=J2a, J2iso=J2iso, J4=J4, g2=g2, dV=dV)


def q3_diag(Cst, Ns=300, Nr=30, Np=64, rho_max=1.1):
    bf = BlendField(Cst=Cst)
    s_grid, Gamma, T, N1, N2, k1, k2, L = trefoil_frame(Ns)
    rho_grid = (np.arange(Nr)+0.5)*(rho_max/Nr); psi_grid = np.arange(Np)*(2*np.pi/Np)
    S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
    cosp, sinp = np.cos(PSI), np.sin(PSI)
    er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
    Xg = Gamma[:, None, None, :] + RHO[..., None]*er
    W = partition_weight(bf, Xg, S, RHO, L, rho_max)
    return shape_diagnostics(Gamma, T, N1, N2, k1, k2, s_grid, rho_grid, psi_grid,
                             lambda X, S, R, P: bf.n_at(X), node_weight=W)


def q2_diag(Cst, R=3.0, q=2, Ns=300, Nr=30, Np=64, rho_max=1.1):
    s_grid, Gamma, T, N1, N2, k1, k2 = circular_frame(R, Ns)
    rho_grid = (np.arange(Nr)+0.5)*(rho_max/Nr); psi_grid = np.arange(Np)*(2*np.pi/Np)
    build = circular_field_grid(Cst, q, R)
    return shape_diagnostics(Gamma, T, N1, N2, k1, k2, s_grid, rho_grid, psi_grid, build)


def self_consistent_beta(d):
    """Unique physical beta at which this shape satisfies BOTH V=phi and sopt=1.
    V=phi fixes beta_eff=x_phi via I(x_phi)=phi*J2a; sopt=1 then fixes
    beta = phi^6 J4 x_phi /(J2a(1+mu phi)).  (Same-grid ratio Q3/Q2 is the robust
    comparison; absolute value is normalization-dependent.)  Returns (beta_sc, xphi)
    or (None,..) if V=phi unreachable (V_max<phi)."""
    g2, dV = d['g2'], d['dV']; J2a, J4 = d['J2a'], d['J4']
    target = PHI*J2a
    if (g2*dV).sum() < target:            # I(0)=J2iso_bare < phi*J2a  => unreachable
        return None, None
    lo, hi = 0.0, 1e6
    for _ in range(80):                    # bisection on I(x)=phi*J2a (I decreasing in x)
        x = 0.5*(lo+hi); Ix = (g2/(1+x*g2)*dV).sum()
        if Ix > target: lo = x
        else: hi = x
    xphi = 0.5*(lo+hi)
    beta_sc = LAM*J4*xphi/(J2a*(1+MU*PHI))
    return beta_sc, xphi


def report(tag, d):
    Vmax = d['J2iso']/d['J2a']
    bsc, xphi = self_consistent_beta(d)
    bstr = f"beta_sc={bsc:.4f}" if bsc is not None else "beta_sc=NONE (V_max<phi)"
    print(f"  {tag:>18}  V_max={Vmax:6.3f}  {'>=phi' if Vmax>=PHI else '<phi!':>5}   {bstr}")
    for beta in [0.35, 0.452]:
        g2, dV = d['g2'], d['dV']
        Jbare = (g2*dV).sum(); Jfb = (g2/(1+beta*g2)*dV).sum()
        fe_sat = ((g2*dV)[beta*g2 > 1].sum())/Jbare
        print(f"       beta={beta:.3f}: suppress {Jbare/Jfb:4.2f}x, {100*fe_sat:4.1f}% saturated")
    return Vmax, bsc


if __name__ == "__main__":
    print("="*74)
    print(" Lepton(Q=2 unknot torus) vs Quark(Q=3 trefoil): V_max & feedback saturation")
    print("="*74)
    print("\n(b) V_max is scale-invariant & beta-independent; V=phi reachable iff V_max>=phi.")
    print("(a) saturation = fraction of isotropic energy where feedback is pinned (beta*g2>1).\n")
    for Cst in [1.5, 2.0, 2.5062]:
        print(f"--- tube sharpness C*={Cst} (radius {1/Cst:.3f}) ---")
        report(f"Q=3 trefoil", q3_diag(Cst))
        report(f"Q=2 torus(q=2)", q2_diag(Cst))
        print()
    print("READ: matched-C* rows isolate TOPOLOGY (knot vs unknot); across C* isolates")
    print("FATNESS. If Q=2 has higher V_max and/or lower saturation, that separates lepton")
    print("from quark. If they match at matched C*, the split is STABILITY not shape ->")
    print("the |grad E|^2 saddle-search is the needed (b) validator.")
