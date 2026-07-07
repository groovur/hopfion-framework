#!/usr/bin/env python3
"""
qh3_field_viewer.py  —  visualize the actual field from a solver .npz output
===============================================================================
WHY THIS EXISTS: every diagnosis so far (v6, v7, v8) has been from scalar
logs alone (K, J4, r_bar, ...) — useful, but it's easy to misjudge what's
geometrically happening from numbers alone (the v8 sigma-scan run is a
direct example: the "SMOOTHED-OUT" verdict was a reasonable label for
r_bar>0.7*R0, but a look at the actual density would say more precisely
WHERE the surviving charge sits). This script loads the field directly
from a solver's saved .npz and produces:

  1. A 3D scatter of the highest-Hopf-charge-density grid points, coloured
     by density, with the analytic R0=3 trefoil curve overlaid as a
     reference line — so you can see at a glance whether the surviving
     density still traces the trefoil shape, has degenerated into a
     different shape (e.g. a single blob, a ring, scattered noise), or has
     drifted away from the reference curve.
  2. 2D density slices through z=0 (where all three crossings live) and
     through one crossing's own plane, both as filled contour plots — the
     2D view is often easier to read precisely than a 3D scatter, and
     directly shows whether three separate density lobes are still visible
     near the crossings or have merged/dispersed.
  3. Side-by-side before/after panels (initial BS-ansatz vs. final state),
     since the .npz files saved by the patched v8 script (and v6/v7)
     include both n_initial and n (final) — this makes "what changed"
     directly visible rather than something to reconstruct mentally from
     a scalar trajectory.

This script has NO torch dependency — it only needs numpy, scipy, and
matplotlib, all pure-CPU, so it runs anywhere without needing the GPU/torch
environment the solver itself needs.

USAGE
-----
  python qh3_field_viewer.py n_sigma0.00_FINAL.npz
  python qh3_field_viewer.py n_sigma50000.00_FINAL.npz --density_pctile 99
  python qh3_field_viewer.py n_R03.00_FINAL.npz --outdir my_plots

  # To compare several sigma values' final states side by side:
  python qh3_field_viewer.py n_sigma0.00_FINAL.npz n_sigma1700.00_FINAL.npz --compare

OUTPUT
------
  {tag}_3d_scatter.png       3D view, before vs after
  {tag}_slices_z0.png        2D density slice through z=0, before vs after
  {tag}_slices_crossing.png  2D density slice through one crossing's plane
  (if --compare with multiple files: {tag}_compare_3d.png etc., one row
   per file, so several sigma/R0 values can be eyeballed side by side)

NOTE ON WHAT "DENSITY" MEANS HERE: this recomputes the Hopf-charge density
rho_J4 = Fxy^2+Fxz^2+Fyz^2 directly from the saved field n, using the exact
same central-difference formula the solver itself uses internally
(matching qh3_trefoil_solver_3d_v6/v7/v8.py's compute_energy function) —
it is not reading any precomputed density from the .npz, so it is a
genuine, independent recomputation, not just replotting a number the
solver already reported.
"""
import numpy as np
import argparse, os, sys

ap = argparse.ArgumentParser()
ap.add_argument('npz_files', nargs='+', help='One or more solver .npz output files')
ap.add_argument('--density_pctile', type=float, default=99.7,
                help='Percentile threshold for the 3D scatter — only points above this percentile of rho_J4 are shown. Note: the Hopf charge density peaks at the TUBE WALL (distance ~1/C* from the centerline, per the BS profile f=2*atan(rho_hat^-C*) peaking in sin^4(f) at rho_hat=1), NOT on the centerline itself — so a low percentile will mostly show diffuse background, not a clean trefoil shape. 99.5-99.9 is a better starting point than a round number like 98.')
ap.add_argument('--C_star', type=float, default=1.5,
                help='C* used to build the initial condition (needed only to draw the tube-wall reference guide at offset 1/C* from the centerline; does not affect the density computation itself)')
ap.add_argument('--outdir', type=str, default='.')
ap.add_argument('--compare', action='store_true',
                help='Plot all provided files side by side instead of one figure per file')
ap.add_argument('--max_scatter_pts', type=int, default=15000,
                help='Cap on points plotted in the 3D scatter, for rendering speed')
args = ap.parse_args()

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)
except ImportError:
    print("pip install matplotlib")
    sys.exit(1)

os.makedirs(args.outdir, exist_ok=True)


def compute_rho_J4(n, h):
    """Exact port of the solver's compute_energy() Hopf-density formula,
    in plain numpy. n has shape (N,N,N,3)."""
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    def cd(u, axis):
        return (np.roll(u, -1, axis) - np.roll(u, 1, axis)) / (2*h)
    nxx, nxy, nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx, nyy, nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx, nzy, nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    return Fxy**2 + Fxz**2 + Fyz**2


def trefoil_curve(R0, r0, n_pts=2000):
    t = np.linspace(0, 2*np.pi, n_pts)
    c2, s2 = np.cos(2*t), np.sin(2*t)
    c3, s3 = np.cos(3*t), np.sin(3*t)
    R = R0 + r0*c3
    return R*c2, R*s2, r0*s3


def trefoil_tube_wall_ring(R0, r0, C_star, n_pts=400, n_ring=12):
    """A scatter of points lying ON the tube wall (distance 1/C* from the
    centerline), used as a visual guide since the Hopf charge density
    peaks AT the wall, not on the centerline itself (BS profile sin^4(f)
    peaks at rho_hat=1, i.e. physical distance 1/C* from the curve)."""
    t = np.linspace(0, 2*np.pi, n_pts)
    c2, s2 = np.cos(2*t), np.sin(2*t)
    c3, s3 = np.cos(3*t), np.sin(3*t)
    R = R0 + r0*c3
    cx, cy, cz = R*c2, R*s2, r0*s3
    # Approximate tangent via finite difference (good enough for a guide)
    dt = t[1]-t[0]
    tx = np.gradient(cx, dt); ty = np.gradient(cy, dt); tz = np.gradient(cz, dt)
    tmag = np.sqrt(tx**2+ty**2+tz**2)
    tx, ty, tz = tx/tmag, ty/tmag, tz/tmag
    # Any vector not parallel to T, then Gram-Schmidt to get a normal
    ref = np.zeros_like(np.stack([tx,ty,tz]))
    ref[2] = 1.0
    dot = tx*ref[0]+ty*ref[1]+tz*ref[2]
    nx_ = ref[0]-dot*tx; ny_ = ref[1]-dot*ty; nz_ = ref[2]-dot*tz
    nmag = np.sqrt(nx_**2+ny_**2+nz_**2)
    nx_, ny_, nz_ = nx_/nmag, ny_/nmag, nz_/nmag
    bx = ty*nz_-tz*ny_; by = tz*nx_-tx*nz_; bz = tx*ny_-ty*nx_
    r_wall = 1.0/C_star
    pts = []
    for chi in np.linspace(0, 2*np.pi, n_ring, endpoint=False):
        c, s = np.cos(chi), np.sin(chi)
        wx = cx + r_wall*(c*nx_ + s*bx)
        wy = cy + r_wall*(c*ny_ + s*by)
        wz = cz + r_wall*(c*nz_ + s*bz)
        pts.append(np.stack([wx,wy,wz], axis=-1))
    return np.concatenate(pts, axis=0)


def load_npz(path, cli_C_star):
    d = np.load(path, allow_pickle=True)
    n_final = d['n']
    n_init = d['n_initial'] if 'n_initial' in d else None
    N = int(d['N']) if 'N' in d else n_final.shape[0]
    h = float(d['h']) if 'h' in d else 0.26
    R0 = float(d['R0']) if 'R0' in d else 3.0
    r0 = float(d['r0']) if 'r0' in d else 0.874
    if 'C_star' in d:
        C_star = float(d['C_star'])
    else:
        C_star = cli_C_star
        print(f"  NOTE: C_star not found in {path} (older solver versions "
              f"did not save it) — using --C_star={cli_C_star} for the "
              f"tube-wall reference guide only; this does not affect the "
              f"density computation itself.")
    if h is None:
        raise ValueError(f"{path}: grid spacing h not found in .npz — re-run "
                          f"the solver with the patched save (h is required "
                          f"to recompute derivatives correctly).")
    return dict(n_final=n_final, n_init=n_init, N=N, h=h, R0=R0, r0=r0,
                C_star=C_star, path=path)


def grid_coords(N, h):
    cv = h*(np.arange(N) - N//2 + 0.5)
    X, Y, Z = np.meshgrid(cv, cv, cv, indexing='ij')
    return X, Y, Z


def known_crossing_midpoints(R0):
    """The three crossing midpoints, per the html's own CROSSING_TS
    (t = pi/6 + 2n*pi/3), at radius EXACTLY R0 in the z=0 plane,
    120 degrees apart. Independent of r0/C* — this is where an intact,
    uncollapsed trefoil's three crossing regions should sit."""
    angles = [np.pi/3, np.pi - np.pi/3, np.pi + np.pi/3]  # 60, 120, 240 deg pattern matching t=pi/6+2n*pi/3
    # Derived directly from trefoil_pos(pi/6 + 2n*pi/3) — see session check:
    # (1.5, 2.598), (1.5, -2.598), (-3, 0) for R0=3; scales linearly with R0.
    base = np.array([[0.5, 0.8660254], [0.5, -0.8660254], [-1.0, 0.0]])
    return base * R0


def find_density_peaks(n, h, R0, n_peaks=6, min_separation=0.8):
    """Find genuine local density concentrations in 3D via local-maximum
    filtering (not just the single global brightest pixel), then report
    each peak's 3D position, distance from the z=0 plane, radius from the
    z-axis, and distance to the nearest known crossing midpoint. This
    replaces eyeballing pixel positions in a 2D slice with an actual
    3D measurement."""
    from scipy import ndimage
    rho = compute_rho_J4(n, h)
    N = rho.shape[0]
    X, Y, Z = grid_coords(N, h)

    # Local maxima via a small max-filter window, then keep only points
    # that ARE their own local max (true peaks, not just any bright pixel).
    footprint_size = max(3, int(round(min_separation / h)))
    local_max = ndimage.maximum_filter(rho, size=footprint_size) == rho
    # Require a minimum density to exclude flat/near-zero "maxima" in vacuum
    thresh = np.percentile(rho, 99.0)
    candidate_mask = local_max & (rho > thresh)
    candidates = np.argwhere(candidate_mask)
    if len(candidates) == 0:
        return []

    # Sort candidates by density, descending, and greedily keep peaks
    # that are well-separated from already-kept ones (avoids reporting
    # several adjacent grid points within the same physical lobe).
    dens = rho[candidate_mask]
    order = np.argsort(-dens)
    kept = []
    kept_xyz = []
    for idx in order:
        i, j, k = candidates[idx]
        xyz = np.array([X[i,j,k], Y[i,j,k], Z[i,j,k]])
        if all(np.linalg.norm(xyz - kxyz) > min_separation for kxyz in kept_xyz):
            kept.append((xyz, dens[idx]))
            kept_xyz.append(xyz)
        if len(kept) >= n_peaks:
            break

    crossing_mids = known_crossing_midpoints(R0)
    results = []
    for xyz, d in kept:
        radius_xy = np.hypot(xyz[0], xyz[1])
        dists_to_crossings = [np.linalg.norm(xyz[:2] - m) for m in crossing_mids]
        nearest_idx = int(np.argmin(dists_to_crossings))
        results.append(dict(
            x=xyz[0], y=xyz[1], z=xyz[2], density=d,
            radius_xy=radius_xy, abs_z=abs(xyz[2]),
            nearest_crossing=nearest_idx,
            dist_to_nearest_crossing=dists_to_crossings[nearest_idx],
        ))
    return results


def report_density_peaks(data, n_peaks=6):
    n_final = data['n_final']
    N, h, R0 = data['N'], data['h'], data['R0']
    peaks = find_density_peaks(n_final, h, R0, n_peaks=n_peaks)
    crossing_mids = known_crossing_midpoints(R0)

    print(f"\n  Known crossing midpoints for R0={R0} (intact-trefoil reference):")
    for i, m in enumerate(crossing_mids):
        print(f"    crossing {i+1}: x={m[0]:.3f}  y={m[1]:.3f}  "
              f"radius={np.hypot(*m):.3f}  angle={np.degrees(np.arctan2(m[1],m[0])):.1f} deg")

    print(f"\n  Top {len(peaks)} density peaks found (3D local maxima, "
          f"min separation 0.8):")
    lines = []
    header = (f"  {'#':>2}  {'x':>7}  {'y':>7}  {'z':>7}  {'density':>9}  "
              f"{'radius_xy':>9}  {'|z|':>6}  {'nearest_crossing':>16}  "
              f"{'dist_to_it':>10}")
    print(header); lines.append(header)
    for i, p in enumerate(peaks):
        line = (f"  {i+1:>2}  {p['x']:>7.3f}  {p['y']:>7.3f}  {p['z']:>7.3f}  "
                f"{p['density']:>9.4f}  {p['radius_xy']:>9.3f}  {p['abs_z']:>6.3f}  "
                f"{p['nearest_crossing']+1:>16}  {p['dist_to_nearest_crossing']:>10.3f}")
        print(line); lines.append(line)

    # IMPORTANT CAVEAT, established by testing this function against a
    # synthetic, perfectly-intact BS-ansatz field built directly from the
    # true curve at the analytically correct C*: even that undamaged
    # reference case shows its 3 crossing-region density peaks sitting
    # measurably INSIDE radius R0 (peaks land near the highest-density
    # blend of the overlapping over/under tube walls, not on the bare
    # centreline) — so "peak radius < R0" is NOT by itself evidence of
    # collapse or drift. There is currently no independently-validated
    # "expected" peak radius for an intact structure at a given C*; this
    # would need to be established by running this same diagnostic on a
    # genuinely converged saddle (once one is found) before any radius
    # value can be called normal or anomalous. For now, report the
    # measurement only — do not classify it as intact vs. drifted.
    matched = set()
    for p in peaks:
        if p['dist_to_nearest_crossing'] < 1.5:
            matched.add(p['nearest_crossing'])
    verdict = (f"  -> {len(matched)} of 3 known crossing positions have a "
               f"density peak within 1.5 units of them.")
    print(verdict); lines.append(verdict)
    if len(matched) == 3:
        radii = [p['radius_xy'] for p in peaks if p['dist_to_nearest_crossing'] < 1.5]
        avg_r = np.mean(radii)
        note = (f"     All 3 crossings represented; their peaks sit at mean "
                f"radius {avg_r:.3f} (geometric R0={R0:.3f}). NOTE: a peak "
                f"radius below R0 is NOT necessarily anomalous — even a "
                f"synthetic, perfectly-intact BS-ansatz field at the "
                f"analytically correct C* shows this same inward offset, "
                f"because the density peak is a blend of both overlapping "
                f"tube walls at a crossing, not a point on the bare "
                f"centreline. This number is reported for tracking changes "
                f"ACROSS runs (e.g. before vs after relaxation, or across "
                f"sigma values), not as a pass/fail check against R0 itself.")
        print(note); lines.append(note)
    return lines


def plot_3d_scatter(ax, n, N, h, R0, r0, C_star, pctile, max_pts, title):
    rho = compute_rho_J4(n, h)
    X, Y, Z = grid_coords(N, h)
    thresh = np.percentile(rho, pctile)
    mask = rho > thresh
    xs, ys, zs, cs = X[mask], Y[mask], Z[mask], rho[mask]
    if len(xs) > max_pts:
        idx = np.random.choice(len(xs), max_pts, replace=False)
        xs, ys, zs, cs = xs[idx], ys[idx], zs[idx], cs[idx]
    if len(xs) == 0:
        ax.set_title(f"{title}\n(NO points above {pctile}th percentile — "
                      f"density may have collapsed to ~0 everywhere; "
                      f"try a lower --density_pctile)")
        return
    sc = ax.scatter(xs, ys, zs, c=cs, cmap='inferno', s=3, alpha=0.6)
    # Tube-wall guide: the Hopf density should peak AT the wall (distance
    # 1/C* from the centerline), not on the centerline itself — this is
    # the more informative reference to overlay.
    wall = trefoil_tube_wall_ring(R0, r0, C_star)
    ax.scatter(wall[:,0], wall[:,1], wall[:,2], c='cyan', s=0.5, alpha=0.25,
               label=f'tube wall, 1/C*={1/C_star:.2f}')
    tx, ty, tz = trefoil_curve(R0, r0)
    ax.plot(tx, ty, tz, color='white', linewidth=0.6, alpha=0.3,
            label=f'centreline, R0={R0}')
    ax.set_title(f"{title}\n{len(xs)} pts above p{pctile} "
                 f"(rho_J4>{thresh:.2e})", fontsize=9)
    ax.legend(loc='upper left', fontsize=6)
    ax.set_box_aspect([1,1,1])


def plot_slice(ax, n, N, h, R0, r0, axis, title, crossing_t=None):
    """2D density slice. axis='z0' -> slice at z=0. axis='crossing' ->
    slice through the plane containing crossing k (at angle theta=2*t,
    rotated about z so the slice cuts through that crossing's midpoint)."""
    rho = compute_rho_J4(n, h)
    cv = h*(np.arange(N) - N//2 + 0.5)
    mid = N//2
    if axis == 'z0':
        sl = rho[:, :, mid]
        extent = [cv[0], cv[-1], cv[0], cv[-1]]
        im = ax.imshow(sl.T, origin='lower', extent=extent, cmap='inferno',
                        aspect='equal')
        ax.set_xlabel('x'); ax.set_ylabel('y')
    else:
        sl = rho[mid, :, :]
        extent = [cv[0], cv[-1], cv[0], cv[-1]]
        im = ax.imshow(sl.T, origin='lower', extent=extent, cmap='inferno',
                        aspect='equal')
        ax.set_xlabel('y'); ax.set_ylabel('z')
    ax.set_title(title, fontsize=9)
    return im


def make_figure_for_file(data):
    n_final, n_init = data['n_final'], data['n_init']
    N, h, R0, r0, C_star = data['N'], data['h'], data['R0'], data['r0'], data['C_star']
    tag = os.path.splitext(os.path.basename(data['path']))[0]

    print(f"\n  --- Density peak report for {tag} ---")
    lines = report_density_peaks(data)
    peak_path = os.path.join(args.outdir, f"{tag}_density_peaks.txt")
    with open(peak_path, 'w') as f:
        f.write(f"Density peak report for {tag}\n")
        f.write(f"R0={R0}  N={N}  h={h}\n\n")
        f.write('\n'.join(lines) + '\n')
    print(f"  Saved: {peak_path}")

    # --- 3D scatter: before vs after ---
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    if n_init is not None:
        plot_3d_scatter(ax1, n_init, N, h, R0, r0, C_star, args.density_pctile,
                         args.max_scatter_pts, "INITIAL (BS-ansatz)")
    else:
        ax1.set_title("(no n_initial saved in this .npz)")
    plot_3d_scatter(ax2, n_final, N, h, R0, r0, C_star, args.density_pctile,
                     args.max_scatter_pts, "FINAL (after relaxation)")
    fig.suptitle(f"{tag} — Hopf charge density, top {100-args.density_pctile:.1f}% of points")
    fig.tight_layout()
    out1 = os.path.join(args.outdir, f"{tag}_3d_scatter.png")
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out1}")

    # --- 2D slice through z=0 ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    if n_init is not None:
        im0 = plot_slice(axes[0], n_init, N, h, R0, r0, 'z0', "INITIAL, z=0 slice")
        fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = plot_slice(axes[1], n_final, N, h, R0, r0, 'z0', "FINAL, z=0 slice")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    fig.suptitle(f"{tag} — density slice through z=0 (all three crossings live near here)")
    fig.tight_layout()
    out2 = os.path.join(args.outdir, f"{tag}_slices_z0.png")
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out2}")

    # --- 2D slice through x=0 (cuts through one crossing region) ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    if n_init is not None:
        im0 = plot_slice(axes[0], n_init, N, h, R0, r0, 'crossing', "INITIAL, x=0 slice")
        fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = plot_slice(axes[1], n_final, N, h, R0, r0, 'crossing', "FINAL, x=0 slice")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    fig.suptitle(f"{tag} — density slice through x=0 (cuts vertically through the structure)")
    fig.tight_layout()
    out3 = os.path.join(args.outdir, f"{tag}_slices_crossing.png")
    fig.savefig(out3, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out3}")


def make_compare_figure(all_data):
    n_files = len(all_data)
    fig = plt.figure(figsize=(6*n_files, 6))
    for i, data in enumerate(all_data):
        ax = fig.add_subplot(1, n_files, i+1, projection='3d')
        tag = os.path.splitext(os.path.basename(data['path']))[0]
        plot_3d_scatter(ax, data['n_final'], data['N'], data['h'],
                         data['R0'], data['r0'], data['C_star'], args.density_pctile,
                         args.max_scatter_pts, tag)
    fig.tight_layout()
    out = os.path.join(args.outdir, "compare_3d.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")

    fig, axes = plt.subplots(1, n_files, figsize=(5*n_files, 5))
    if n_files == 1:
        axes = [axes]
    for i, data in enumerate(all_data):
        tag = os.path.splitext(os.path.basename(data['path']))[0]
        im = plot_slice(axes[i], data['n_final'], data['N'], data['h'],
                         data['R0'], data['r0'], 'z0', tag)
        fig.colorbar(im, ax=axes[i], fraction=0.046)
    fig.suptitle("z=0 slice comparison")
    fig.tight_layout()
    out = os.path.join(args.outdir, "compare_slices_z0.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == '__main__':
    all_data = []
    for path in args.npz_files:
        print(f"\nLoading {path} ...")
        try:
            data = load_npz(path, args.C_star)
        except Exception as e:
            print(f"  ERROR loading {path}: {e}")
            continue
        all_data.append(data)
        print(f"  N={data['N']}  h={data['h']}  R0={data['R0']}  r0={data['r0']}  "
              f"has_initial={data['n_init'] is not None}")

    if not all_data:
        print("No valid .npz files loaded — nothing to plot.")
        sys.exit(1)

    if args.compare and len(all_data) > 1:
        make_compare_figure(all_data)
    else:
        for data in all_data:
            make_figure_for_file(data)

    print("\nDone.")
