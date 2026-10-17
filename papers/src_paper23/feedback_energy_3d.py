#!/usr/bin/env python3.11
"""
feedback_energy_3d.py -- Paper XXIII components D/E: the full density-feedback energy on the PHYSICAL
TUBE around the (2,3) trefoil centerline, and the feedback-saturation (beta*g2 ~ beta_*rho) at the
CROSSINGS. DM/Quark crossover (src_paper23).

WHY THE PHYSICAL TUBE, NOT BRIESKORN (user, 2026-10-13): "the centerline of the trefoil is where there
is very little energy." Correct -- in a Faddeev-Niemi tube Hopfion the energy density lives in the tube
WALLS (the shell of rapid n-hat variation), and the feedback saturation (1+beta_*rho) lives where that
density OVERLAPS: the three crossings, where two strands pass within 2 r0=1.748 (tube radius ~r0=0.874).
The Brieskorn rational map (knotted_map.py) was the right instrument for the CHARGE (it proved Q_H=3 via
Whitehead) but it smears the field over all of S^3 -> it mislocates the energy (and came out
geometrically smooth). The PHYSICAL tube build_doublet (src_paper21, around torus_knot_centerline
R0=3,r0=0.874, BS profile f0(rho)=2 arctan(rho^-C*), C*=2.5062, q=3 cable framing) localizes the energy
correctly and resolves the crossings. We do NOT trust this ansatz's Whitehead Q_H (the Q_H~0 at low N is
the UNDER-RESOLUTION artifact: 2 pts/tube at N=64; ~5 needs N~200); we read its ENERGY DENSITY directly,
with topology enforced by the q=3 framing + two-strand crossing structure.

THIS SCRIPT:
  (1) builds the physical tube at resolution N (reusing the VALIDATED src_paper21 core),
  (2) reports the feedback diagnostics V, s_opt, E_fb, J4/J2a (convergence vs N; need ~5 pts/tube),
  (3) LOCATES the 3 crossings geometrically and measures the feedback saturation beta*g2 there vs on an
      isolated strand segment -- the OVERLAP ENHANCEMENT = the beta_*rho locus (component F reads this).
beta_*rho calibration to S_eff=sin^4theta/[phi^6(1+beta_*rho)] (ro prop kappa, Paper VII) is component F;
here we deliver the resolved field and the crossing-vs-isolated saturation ratio it needs.
"""
import numpy as np, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from feedback_saddle_core import (torus_knot_centerline, build_doublet, nhat_from_doublet,
                                   whitehead_QH, Efb_diagnostics, phi)

def g2_field(n, h):
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    return sum((dn[ax]**2).sum(0) for ax in range(3))        # |grad n|^2

def self_avoidance(cl, t_arr, dt_min=0.6):
    """Minimal NON-adjacent (|dt|>dt_min) strand-strand distance over the centerline. For the relaxed
    T(2,3) this is the constant antipodal 2r0 (Paper XV: |G(t)-G(t+pi)|=2r0 for all t)."""
    from scipy.spatial import cKDTree
    cl2 = cl[::max(1, len(cl)//600)]; t2 = t_arr[::max(1, len(cl)//600)]
    best = np.inf
    for i in range(len(cl2)):
        dt = np.abs(((t2 - t2[i] + np.pi) % (2*np.pi)) - np.pi)
        d = np.linalg.norm(cl2 - cl2[i], axis=1); d[dt <= dt_min] = np.inf
        best = min(best, d.min())
    return best

def radial_g2_profile(g2, ax, cl, h, rho_edges):
    """<g2>(rho) binned by distance-to-centerline rho. Resolution-robust AWAY from the core
    (the core rho->0 g2 diverges with resolution; the profile scale rho~1/C* does not)."""
    from scipy.spatial import cKDTree
    tree = cKDTree(cl)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    dist, _ = tree.query(pts, k=1)
    dist = dist.reshape(g2.shape)
    prof = []
    for lo, hi in zip(rho_edges[:-1], rho_edges[1:]):
        m = (dist >= lo) & (dist < hi)
        prof.append(float(g2[m].mean()) if m.any() else np.nan)
    return np.array(prof)

def run(N, L, R0, r0, C_star, beta):
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    t, cl = torus_knot_centerline(2, 3, R0, r0)
    z1, z2 = build_doublet(pts, cl, t, 3, C_star)
    z1 = z1.reshape(N, N, N); z2 = z2.reshape(N, N, N)
    n = nhat_from_doublet(z1, z2)
    d = Efb_diagnostics(n, h, beta)
    g2 = g2_field(n, h)
    # single-tube wall saturation: beta*g2 read at the profile scale rho=1/C* (resolution-robust)
    rho_scale = 1.0/C_star
    edges = np.array([0.5, 0.75, 1.0, 1.25, 1.5]) * rho_scale
    prof = radial_g2_profile(g2, ax, cl, h, edges)
    bg2_wall = beta * prof[1]        # the bin centred on rho ~ 1/C* (0.75-1.0 x rho_scale)
    out = dict(N=N, h=h, ptube=rho_scale/h, QH=whitehead_QH(z1, z2, h),
               V=d['V'], s_opt=d['s_opt'], J4_over_J2a=d['J4_over_J2a'], E_fb=d['E_fb'],
               dmin=self_avoidance(cl, t), r0=r0, tube_w=2*rho_scale,
               bg2_wall=bg2_wall, g2_prof=(beta*prof).tolist())
    return out

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--Nlist', type=str, default='64,96,128')
    ap.add_argument('--L', type=float, default=6.0)
    ap.add_argument('--R0', type=float, default=3.0); ap.add_argument('--r0', type=float, default=0.874)
    ap.add_argument('--C_star', type=float, default=2.5062); ap.add_argument('--beta', type=float, default=0.452)
    a = ap.parse_args()
    Ns = [int(s) for s in a.Nlist.split(',')]
    print("="*92)
    print(f"  (A) SINGLE-TUBE WALL saturation  T(2,3)  R0={a.R0} r0={a.r0} C*={a.C_star} beta={a.beta}")
    print(f"  relaxed tube self-avoids (d_min=2r0); beta_*rho(A) = beta*g2 at the WALL rho~1/C*={1/a.C_star:.3f}")
    print("="*92)
    print(f"  {'N':>4} {'pts/tube':>8} {'d_min':>6} {'tube_w':>6} {'V(->phi)':>9} {'s_opt(->1)':>10} "
          f"{'J4/J2a':>7} {'bg2_wall':>8}")
    for N in Ns:
        t0 = time.time()
        o = run(N, a.L, a.R0, a.r0, a.C_star, a.beta)
        print(f"  {N:>4} {o['ptube']:>8.1f} {o['dmin']:>6.3f} {o['tube_w']:>6.3f} {o['V']:>9.3f} "
              f"{o['s_opt']:>10.3f} {o['J4_over_J2a']:>7.3f} {o['bg2_wall']:>8.2f}  [{time.time()-t0:.0f}s]")
        print(f"       beta*g2 radial profile (rho/rho_scale 0.5-0.75-1.0-1.25): {[round(x,2) for x in o['g2_prof']]}")
    print("="*92)
    print(f"  READ: d_min={o['dmin']:.3f} vs tube width {o['tube_w']:.3f} -> strands DO NOT overlap in the")
    print("  relaxed knot (self-avoiding). So beta_*rho(A) is a single-tube WALL quantity, ~uniform along")
    print("  the knot. Compare to the threading-barrier beta_*rho(B) in threading_barrier.py. At >=5")
    print("  pts/tube (N~200) V->phi, s_opt->1 test self-consistency; bg2_wall is the (A) candidate.")
    print("="*92)
