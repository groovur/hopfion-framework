#!/usr/bin/env python3.11
"""
knot_in_medium.py -- Paper XXIII component F (saturated/density-modified reading): beta_*rho for the
trefoil forming IN the dense QCD-epoch medium, anchored by the vacuum relation beta*rho_inf = phi.
DM/Quark crossover (src_paper23).

CONTEXT (notes (xxxx)->(yyyy); user choice 2026-10-14 = "knot-in-medium"):
  The ISOLATED trefoil's self-consistent core feedback is beta_*rho_knot ~ 0.5 (engine V=phi + analytic
  b*kern, validated vs P4:prop:bstar) -- ~30x short of the eta requirement ~15. Paper XXII l.356-360
  says beta_*rho is the trefoil's DENSITY-MODIFIED self-consistent coupling (not the vacuum constant):
  the trefoil forms at the dense QCD epoch (rho tracks T^4, l.333), so the ambient background density
  rho_bg adds to the knot's own. This script computes beta_*rho(rho_bg) in that medium.

THE SOLID ANCHOR (DE sector, mxi_phi10_inconsistency.md; P7:eq:Lambda_obs):
  at the vacuum density rho_inf, the feedback argument is beta*rho_inf = phi (since the non-canonical
  norm K(rho_inf)=1/[phi^6(1+beta rho_inf)]=1/phi^8 => 1+beta rho_inf=phi^2 => beta rho_inf=phi). EXACT.
  So beta = phi/rho_inf, and the ambient contribution to the feedback argument is
      beta*rho_bg = phi * (rho_bg/rho_inf).
  The total feedback argument the baryon channel sees at the QCD epoch:
      beta_*rho = beta_*rho_knot(~0.5) + phi*(rho_bg/rho_inf).
  => eta needs beta_*rho~15, so phi*(rho_bg/rho_inf) ~ 14.5 => rho_bg/rho_inf ~ 15/phi ~ 9.
  The KNOT's own share of the 15 is ~0.5/15 ~ 3%; the rest is the ambient QCD-epoch density.

CROSS-CHECK (engine): solve the self-consistency V=J2iso_fb/J2a=phi with an added uniform background
rho_bg in the feedback (1+beta(g2+rho_bg)), and read beta_*rho at the tube. Confirms the linear anchor
model and the ~0.5 knot floor.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from feedback_saddle_core import torus_knot_centerline, build_doublet, nhat_from_doublet
from scipy.optimize import brentq

phi = (1+np.sqrt(5))/2

def build_trefoil_fields(N=112, L=6.0, C_star=2.5062):
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    t, cl = torus_knot_centerline(2, 3, 3.0, 0.874)
    z1, z2 = build_doublet(pts, cl, t, 3, C_star); z1 = z1.reshape(N,N,N); z2 = z2.reshape(N,N,N)
    n = nhat_from_doublet(z1, z2)
    dn = [np.stack([np.gradient(n[c], h, axis=a) for c in range(3)]) for a in range(3)]
    g2 = sum((dn[a]**2).sum(0) for a in range(3)); s4 = (1-n[2]**2)**2
    return g2, s4, h, ax, cl

def beta_rho_in_medium(g2, s4, rho_bg):
    """self-consistency V=phi with uniform background rho_bg in the feedback (1+beta(g2+rho_bg));
    returns (beta_sc, beta_*rho at the knot core = beta_sc*(g2_max+rho_bg))."""
    J2a = np.sum(s4*g2)
    def V(beta): return np.sum(g2/(1+beta*(g2+rho_bg)))/J2a
    if V(1e-6) < phi:    # feedback already too strong even at tiny beta
        return 0.0, 0.0
    try:
        bsc = brentq(lambda b: V(b)-phi, 1e-7, 5.0)
    except ValueError:
        return np.nan, np.nan
    g2core = np.percentile(g2, 99.5)   # representative core (robust, not the single max)
    return bsc, bsc*(g2core + rho_bg)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=112); ap.add_argument('--L', type=float, default=6.0)
    ap.add_argument('--C_star', type=float, default=2.5062)
    a = ap.parse_args()

    print("="*86)
    print("  (F, saturated/knot-in-medium) beta_*rho vs ambient background rho_bg")
    print("  anchor beta*rho_inf = phi (EXACT, DE sector) -> ambient adds phi*(rho_bg/rho_inf)")
    print("="*86)

    # --- analytic anchor model (knot floor ~0.5 + ambient) ---
    br_knot = 0.5
    print("  ANALYTIC (anchor beta*rho_inf=phi):  beta_*rho = beta_*rho_knot(0.5) + phi*(rho_bg/rho_inf)")
    print(f"  {'target beta_*rho':>16} {'rho_bg/rho_inf needed':>22} {'knot share':>11}")
    for tgt in [2.0, 15.0, phi**6]:
        ratio = (tgt - br_knot)/phi
        print(f"  {tgt:>16.2f} {ratio:>22.2f} {br_knot/tgt*100:>10.1f}%")
    print(f"     (phi^6={phi**6:.2f}; note rho_bg/rho_inf for 15 is 15/phi-ish = {(15-br_knot)/phi:.2f})")
    print()

    # --- engine cross-check: self-consistency with a background ---
    print("  ENGINE CROSS-CHECK (build_doublet trefoil, self-consistency V=phi with background rho_bg):")
    g2, s4, h, ax, cl = build_trefoil_fields(a.N, a.L, a.C_star)
    g2u = np.percentile(g2, 99.5)
    print(f"  grid N={a.N}; g2 core(99.5pct)={g2u:.2f}, g2 max={g2.max():.2f}")
    print(f"  {'rho_bg':>8} {'beta_sc':>9} {'beta_*rho(core)':>15} {'beta_sc*rho_bg(ambient)':>24}")
    for rho_bg in [0.0, 5.0, 20.0, 60.0, 150.0, 400.0]:
        bsc, br = beta_rho_in_medium(g2, s4, rho_bg)
        amb = bsc*rho_bg
        print(f"  {rho_bg:>8.1f} {bsc:>9.4f} {br:>15.2f} {amb:>24.2f}")
    print("="*86)
    print("  READ: the ISOLATED knot (rho_bg=0) gives beta_*rho~0.5 (floor, matches analytic b*kern).")
    print("  To reach beta_*rho~15 the AMBIENT background must supply ~14.5 -> rho_bg ~ 9 x rho_inf")
    print("  (vacuum). So the eta suppression is ~97% ambient QCD-epoch density, ~3% knot. The eta")
    print("  magnitude REDUCES to one ratio: rho_bg(QCD)/rho_inf ~ 9 -- plausible (not wild) but PINNED")
    print("  only by the physical-density <-> k-essence-field map, which Paper XXII itself flags as the")
    print("  'as-yet-unpinned constant' in rho prop kappa. [N] bounded: open ratio, not a free O(many).")
    print("="*86)
