#!/usr/bin/env python3.11
"""
threading_barrier.py -- Paper XXIII component F (primary route): the feedback saturation beta_*rho at the
THREADING BARRIER, where the Q_H=1 strand pierces the host tube during trefoil formation (Paper XVIII
satellite threading). DM/Quark crossover (src_paper23).

WHY (user, 2026-10-14): the RELAXED trefoil tube self-avoids (strands 2r0=1.748 apart, tube width ~0.80
-> no overlap; feedback_energy_3d.py). So the self-interaction that saturates the feedback (1+beta_*rho)
is NOT in the relaxed knot -- it is TRANSIENT, at the formation BARRIER, where one strand threads THROUGH
the other's tube (the half-threaded n=4 state, Paper XVIII). This is the user's primary hypothesis for
where the eta-setting beta_*rho lives; (A) single-tube wall (feedback_energy_3d.py) is the comparison.

LOCAL 3D MODEL (keeps the self-interaction -- NOT the vetoed 1D torus reduction): two straight tube
segments carrying the VALIDATED BS profile (build_doublet, src_paper21), brought to minimal separation d:
  host:     axis = x-axis through origin
  threader: axis = z-direction line through (0, d, 0)   -> axis-axis min distance = d at the origin
As d shrinks (threading deepens), the two cores approach; the feedback-relevant condensate DENSITY
(rho prop field energy density prop g2 = |grad n|^2, Paper VII rho prop kappa) ADDS where they overlap:
  rho_total(x) ~ g2_host(x) + g2_threader(x)      (leading-order overlap; the feedback couples to the
                                                   local scalar density, which stacks -- n-hat itself is
                                                   S^2-valued and does not add, but the density does).
beta_*rho(d) = beta * max_overlap(rho_total). The BARRIER = "just before full absorption / after emission
begins" (user) = d ~ tube width, where the cores are about to merge (reconnection). Report beta_*rho vs d,
locate the barrier, and compare beta_*rho(barrier) to the eta requirement (~15, Paper XXII) and to (A).

CAVEAT: the additive-density overlap is leading-order; the fully-relaxed threaded field would need the
relaxation that hit the wall. The d-TREND and the barrier ENHANCEMENT over a single tube are the robust
outputs; the absolute number is calibration-gated (component F / notes (kkkk)).
"""
import numpy as np, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from feedback_saddle_core import build_doublet, nhat_from_doublet, phi

def g2_of(n, h):
    dn = [np.stack([np.gradient(n[c], h, axis=ax) for c in range(3)]) for ax in range(3)]
    return sum((dn[ax]**2).sum(0) for ax in range(3))

def straight_tube_g2(ax, h, axis_point, axis_dir, C_star, q=1, span=8.0, NT=4000):
    """g2 field of a single straight BS-profile tube (reusing the validated build_doublet) whose core is
    the line axis_point + s*axis_dir. Open segment -> build_doublet finds no 'return' strand -> single tube."""
    axis_dir = np.asarray(axis_dir, float); axis_dir /= np.linalg.norm(axis_dir)
    s = np.linspace(-span/2, span/2, NT)
    cl = np.asarray(axis_point, float)[None, :] + s[:, None]*axis_dir[None, :]
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)   # framing parameter along the segment
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    z1, z2 = build_doublet(pts, cl, t_arr, q, C_star)
    N = len(ax)
    n = nhat_from_doublet(z1.reshape(N, N, N), z2.reshape(N, N, N))
    return g2_of(n, h)

def run_d(d, N, L, C_star, beta, q=1):
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    g2_host = straight_tube_g2(ax, h, (0, 0, 0), (1, 0, 0), C_star, q, span=4*L)
    g2_thr  = straight_tube_g2(ax, h, (0, d, 0), (0, 0, 1), C_star, q, span=4*L)
    rho_tot = g2_host + g2_thr
    # overlap region = near the origin segment of closest approach (|x|,|z| small, y in [0,d])
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    near = (np.abs(X) < 1.0) & (np.abs(Z) < 1.0) & (Y > -0.3) & (Y < d + 0.3)
    bg2_overlap = beta * rho_tot[near].max() if near.any() else np.nan
    bg2_single  = beta * g2_host[near].max() if near.any() else np.nan  # one tube, same region
    return dict(d=d, h=h, ptube=(1/C_star)/h, bg2_overlap=float(bg2_overlap),
                bg2_single=float(bg2_single), enh=float(bg2_overlap/bg2_single) if bg2_single else np.nan)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=96)
    ap.add_argument('--L', type=float, default=2.5)
    ap.add_argument('--C_star', type=float, default=2.5062)
    ap.add_argument('--beta', type=float, default=0.452)
    ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--dlist', type=str, default='1.748,1.2,0.8,0.6,0.4,0.25')
    a = ap.parse_args()
    ds = [float(s) for s in a.dlist.split(',')]
    tube_w = 2.0/a.C_star
    print("="*84)
    print(f"  (B) THREADING-BARRIER feedback saturation  C*={a.C_star} beta={a.beta} q={a.q}  N={a.N} L={a.L}")
    print(f"  two tubes at axis separation d; tube width ~{tube_w:.3f}; barrier ~ d<=tube width (cores merge)")
    print("="*84)
    print(f"  {'d':>7} {'pts/tube':>8} {'beta*rho_overlap':>16} {'beta*rho_1tube':>15} {'enhancement':>12}")
    for d in ds:
        t0 = time.time()
        o = run_d(d, a.N, a.L, a.C_star, a.beta, a.q)
        tag = "  <- barrier (d~tube_w)" if abs(d - tube_w) < 0.25 else ""
        print(f"  {o['d']:>7.3f} {o['ptube']:>8.1f} {o['bg2_overlap']:>16.2f} {o['bg2_single']:>15.2f} "
              f"{o['enh']:>11.2f}x  [{time.time()-t0:.0f}s]{tag}")
    print("="*84)
    print("  READ: beta*rho_overlap rises as d shrinks (threading deepens); the BARRIER is d ~ tube width")
    print("  (cores about to merge). beta_*rho(barrier) is the (B) candidate for the eta magnitude -- compare")
    print("  to (A) single-tube wall and to the eta requirement ~15 (Paper XXII). Absolute value is")
    print("  calibration-gated; the enhancement over one tube and the d-trend are the robust outputs.")
    print("="*84)
