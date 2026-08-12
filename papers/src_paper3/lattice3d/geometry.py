"""Lattice geometry: site layout, T(2,2) tube curves, and the pre-solve
geometry validation (nearest-curve distances, Gauss linking integrals).
"""
import numpy as np
from scipy.spatial import cKDTree

from common3d import R0, r0, PHI

N_SIDE = 3  # coloring period along each axis -> 3x3x3-site unit cell

# rotation matrices sending the color-0 (normal-z) torus to the color-1
# (normal-x) and color-2 (normal-y) orientations.
ROT = {
    0: np.eye(3),
    1: np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]]),
    2: np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]]),
}


def site_color(i, j, k):
    return (i + j + k) % 3


def tube_curve(center, color, delta, npts):
    """Sample one T(2,2) tube curve (npts points, t in [0,2pi))."""
    t = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    rr = R0 + r0 * np.cos(t + delta)
    x = rr * np.cos(t)
    y = rr * np.sin(t)
    z = r0 * np.sin(t + delta)
    pts = np.stack([x, y, z], axis=1)
    pts = pts @ ROT[color].T
    pts = pts + np.asarray(center)
    return pts


def core_ring(center, color, npts=2000):
    """The R0-radius core circle of a hopfion (minor radius -> 0), used for
    the coarse tangency / nearest-neighbor distance check."""
    t = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    x = R0 * np.cos(t)
    y = R0 * np.sin(t)
    z = np.zeros_like(t)
    pts = np.stack([x, y, z], axis=1)
    pts = pts @ ROT[color].T
    pts = pts + np.asarray(center)
    return pts


def build_sites(a_c, n_side=N_SIDE):
    """All sites of the n_side^3 unit cell (fractional index i,j,k in
    [0,n_side)); returns list of dict(idx=(i,j,k), center, color)."""
    sites = []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                c = site_color(i, j, k)
                center = a_c * np.array([i, j, k], dtype=float)
                sites.append({'idx': (i, j, k), 'center': center, 'color': c})
    return sites


def build_curves(a_c, n_side=N_SIDE, npts=300):
    """All 2*n_side^3 tube curves of the unit cell, unwrapped (real-space)
    coordinates -- NOT wrapped into [0,L). Each entry:
    dict(site=(i,j,k), tube=+-1, color=c, pts=(npts,3))."""
    curves = []
    for s in build_sites(a_c, n_side):
        for tube, delta in ((+1, 0.0), (-1, np.pi)):
            pts = tube_curve(s['center'], s['color'], delta, npts)
            curves.append({'site': s['idx'], 'tube': tube, 'color': s['color'],
                            'center': s['center'], 'pts': pts})
    return curves


def nearest_neighbor_site_pairs(a_c, n_side=N_SIDE):
    """Axis-adjacent site pairs on the periodic n_side^3 lattice, with the
    correct periodic-image displacement (site index wraps mod n_side; the
    real-space neighbor position is offset by +-L along that axis, not
    wrapped back into the cell). Each axis-adjacent pair always differs in
    color by +-1 mod 3 (perpendicular ring orientation) by construction.
    Returns list of dict(a=(i,j,k), b_center=array, axis=0/1/2)."""
    L = n_side * a_c
    pairs = []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                idx = (i, j, k)
                for axis in range(3):
                    nb = list(idx)
                    nb[axis] += 1
                    wrapped = nb[axis] >= n_side
                    if wrapped:
                        nb[axis] -= n_side
                    b_center = a_c * np.array(nb, dtype=float)
                    if wrapped:
                        b_center[axis] += L
                    pairs.append({'a': idx, 'b_idx': tuple(nb), 'b_center': b_center,
                                  'axis': axis})
    return pairs


def curve_min_distance(pts1, pts2):
    tree = cKDTree(pts2)
    d, _ = tree.query(pts1, k=1)
    return d.min()


def gauss_linking(pts1, pts2):
    """Discretized Gauss linking integral between two closed curves sampled
    at pts1 (n1,3), pts2 (n2,3). O(n1*n2); chunked over pts1 to bound
    memory."""
    n1 = pts1.shape[0]
    n2 = pts2.shape[0]
    d1 = np.roll(pts1, -1, axis=0) - pts1  # (n1,3) tangent-like segment
    d2 = np.roll(pts2, -1, axis=0) - pts2  # (n2,3)
    total = 0.0
    chunk = 200
    for s in range(0, n1, chunk):
        p1 = pts1[s:s + chunk]          # (c,3)
        t1 = d1[s:s + chunk]            # (c,3)
        r = p1[:, None, :] - pts2[None, :, :]      # (c,n2,3)
        rnorm = np.linalg.norm(r, axis=2)
        rnorm3 = np.clip(rnorm, 1e-12, None) ** 3
        cross = np.cross(t1[:, None, :], d2[None, :, :])  # (c,n2,3)
        integrand = np.einsum('ijk,ijk->ij', cross, r) / rnorm3
        total += integrand.sum()
    return total / (4.0 * np.pi)


def validate_geometry(a_c, report=print):
    """Runs the pre-solve geometry validation (a) and (b) described in the
    task. Returns dict with the raw numbers so callers can assert / log."""
    n_side = N_SIDE
    pairs = nearest_neighbor_site_pairs(a_c, n_side)
    curves = build_curves(a_c, n_side, npts=4000)
    curve_lookup = {}
    for c in curves:
        curve_lookup.setdefault(c['site'], []).append(c)

    # (a) minimum inter-curve distance for every nearest-neighbor hopfion pair
    min_dists = []
    for p in pairs:
        a_curves = curve_lookup[p['a']]
        b_curves_raw = curve_lookup[p['b_idx']]
        shift = p['b_center'] - a_c * np.array(p['b_idx'], dtype=float)
        for ca in a_curves:
            for cb in b_curves_raw:
                pts_b = cb['pts'] + shift
                dmin = curve_min_distance(ca['pts'], pts_b)
                min_dists.append(dmin)
    min_dists = np.array(min_dists)
    report(f"[geometry a_c={a_c}] inter-curve min-distance distribution over "
           f"{len(min_dists)} nearest-neighbor tube-pairs (perpendicular by "
           f"construction):")
    report(f"  min={min_dists.min():.5f}  mean={min_dists.mean():.5f}  "
           f"max={min_dists.max():.5f}  (2*r0={2*r0:.5f} = tangency threshold "
           f"for the core-circle picture)")

    # (b) Gauss linking integral, core-ring model, a_c = 5.8, a few pairs
    ac_link = 5.8
    pairs58 = nearest_neighbor_site_pairs(ac_link, n_side)
    checked = 0
    lk_values = []
    for p in pairs58[:6]:
        color_a = site_color(*p['a'])
        color_b = site_color(*p['b_idx'])
        ring_a = core_ring(ac_link * np.array(p['a'], dtype=float), color_a, npts=2000)
        ring_b = core_ring(p['b_center'], color_b, npts=2000)
        lk = gauss_linking(ring_a, ring_b)
        lk_values.append(lk)
        report(f"  [linking a_c=5.8] pair a={p['a']} b={p['b_idx']} axis={p['axis']}: "
               f"Lk={lk:.4f}")
        checked += 1
    lk_values = np.array(lk_values)
    ok = np.allclose(np.abs(lk_values), 1.0, atol=0.05)
    report(f"  |Lk|=1 for all {checked} sampled perpendicular neighbor pairs: {ok}")
    return {'min_dists': min_dists, 'lk_values': lk_values, 'lk_ok': bool(ok)}


# --- Fallback motif -------------------------------------------------------
# The 3-coloring motif above fails the linking validation (see results.md):
# for a fixed pair of colors (c, c+1 mod 3), ALL THREE axis directions from a
# site reach a neighbor of color c+1 (since i+j+k always increases by 1 along
# any single-axis step), but only the axis lying in the intersection line of
# the two rings' planes actually links (Lk=-1); the other two give Lk~0.
# Fallback (documented, used for all downstream physics): alternating xy
# (normal z) / xz (normal y) Hopf-chain rings along x, period 2*a_c, stacked
# with period a_c in y and z (translated copies of the same chain, not
# linked to their y/z neighbors -- only consecutive chain links along x).

def fallback_site_color(i):
    return 0 if i % 2 == 0 else 2  # reuse ROT[0] (normal z) / ROT[2] (normal y)


def build_sites_fallback(a_c, n_stack=1):
    """i in {0,1} (chain period 2*a_c along x); j,k in [0,n_stack) (period
    a_c along y,z)."""
    sites = []
    for i in range(2):
        for j in range(n_stack):
            for k in range(n_stack):
                c = fallback_site_color(i)
                center = np.array([i * a_c, j * a_c, k * a_c])
                sites.append({'idx': (i, j, k), 'center': center, 'color': c})
    return sites


def build_curves_fallback(a_c, n_stack=1, npts=300):
    curves = []
    for s in build_sites_fallback(a_c, n_stack):
        for tube, delta in ((+1, 0.0), (-1, np.pi)):
            pts = tube_curve(s['center'], s['color'], delta, npts)
            curves.append({'site': s['idx'], 'tube': tube, 'color': s['color'],
                            'center': s['center'], 'pts': pts})
    return curves


def validate_geometry_fallback(a_c, report=print):
    """Chain-along-x linking check (period 2*a_c) plus the y/z stacking
    distance check (period a_c, unlinked)."""
    Lx = 2.0 * a_c
    # (b) linking: consecutive chain neighbors i=0 <-> i=1 <-> i=0(next cell)
    pairs = [((0, 0, 0), (1, 0, 0), np.array([0., 0., 0.])),
             ((1, 0, 0), (0, 0, 0), np.array([Lx, 0., 0.]))]
    lk_values = []
    for a_idx, b_idx, b_shift in pairs:
        ca = fallback_site_color(a_idx[0])
        cb = fallback_site_color(b_idx[0])
        ring_a = core_ring(a_c * np.array(a_idx, dtype=float), ca, npts=2000)
        ring_b = core_ring(a_c * np.array(b_idx, dtype=float) + b_shift, cb, npts=2000)
        lk = gauss_linking(ring_a, ring_b)
        lk_values.append(lk)
        report(f"  [fallback linking a_c={a_c}] chain pair {a_idx}<->{b_idx}: Lk={lk:.4f}")
    lk_values = np.array(lk_values)
    ok = np.allclose(np.abs(lk_values), 1.0, atol=0.05)
    report(f"  [fallback] |Lk|=1 for all chain neighbor pairs: {ok}")

    # (a) min inter-curve distance: chain neighbors, both periodic images
    # (site 0 <-> site 1 within the cell, and site 1 <-> site 0 of the next
    # cell, i.e. site-0 image shifted by +Lx)
    curves = build_curves_fallback(a_c, n_stack=1, npts=4000)
    by_site = {}
    for c in curves:
        by_site.setdefault(c['site'], []).append(c)
    min_dists = []
    for ca in by_site[(0, 0, 0)]:
        for cb in by_site[(1, 0, 0)]:
            min_dists.append(curve_min_distance(ca['pts'], cb['pts']))
    for ca in by_site[(1, 0, 0)]:
        for cb in by_site[(0, 0, 0)]:
            min_dists.append(curve_min_distance(ca['pts'], cb['pts'] + np.array([Lx, 0., 0.])))
    min_dists = np.array(min_dists)
    report(f"  [fallback a_c={a_c}] chain-neighbor min-distance: "
           f"min={min_dists.min():.5f} mean={min_dists.mean():.5f} "
           f"(2*r0={2*r0:.5f})")
    return {'lk_values': lk_values, 'lk_ok': bool(ok), 'min_dists': min_dists}


if __name__ == "__main__":
    print("=== primary 3-coloring motif ===")
    for a_c in (6.0, 5.9, 5.8):
        validate_geometry(a_c)
        print()
    print("=== fallback alternating-chain motif ===")
    for a_c in (6.0, 5.9, 5.8):
        validate_geometry_fallback(a_c)
        print()
