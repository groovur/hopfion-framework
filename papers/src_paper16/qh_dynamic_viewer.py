#!/usr/bin/env python3
"""
qh_dynamic_viewer.py  —  visualize time-series output from the integrator
==========================================================================
Produces:
  1. TRAJECTORY PLOTS: J4(t), KE(t), snapback(t), theta_min(t), and
     the critical J4/J4_vacuum ratio.
  2. SPATIAL SLICES at any chosen snapshot: rho_J4 density through
     z=0, x=0, and y=0 planes (vacuum vs snapshot side-by-side).
  3. MULTI-RUN COMPARISON: overlay J4 trajectories across runs.
  4. NPY EXPORT: for each snapshot, writes n_<tag>.npy and (if present)
     v_<tag>.npy as standalone (N,N,N,3) float32 arrays readable
     directly by hopfion_viewer.html.

No torch required — numpy, scipy, matplotlib only.

USAGE
-----
  # Slices + npy export for several snapshots:
  python qh_dynamic_viewer.py --traj run/trajectory.csv \\
      --snaps run/snapshots/n_t00000800.npz \\
              run/snapshots/n_t00001400.npz

  # Export npy only, no slice plots:
  python qh_dynamic_viewer.py \\
      --snaps run/snapshots/n_t00001400.npz --no_slices

  # Shared colorbar scale between vacuum and snapshot panels
  # (useful when fields are comparable in magnitude):
  python qh_dynamic_viewer.py --snaps run/snapshots/n_t00000800.npz --shared_scale

  # Log-scale density (useful when dynamic range is enormous):
  python qh_dynamic_viewer.py --snaps run/snapshots/n_t00000800.npz --log_scale

  # Compare multiple trajectories:
  python qh_dynamic_viewer.py \\
      --traj run_ks5/trajectory.csv run_ks10/trajectory.csv \\
      --labels "kick=5" "kick=10"

COLORBAR MODES
--------------
  Default: each panel scaled independently to its own 99.5th percentile.
    Good for: seeing internal structure when vacuum >> snapshot magnitude.
  --shared_scale: vacuum and snapshot share the same colorbar axis.
    Good for: direct amplitude comparison (e.g. is the ring as bright
    as the main torus?). Reveals under/over-exposure honestly.
  --log_scale: log10(1 + rho) displayed instead of raw rho.
    Good for: blown-up / saturated fields where linear scale loses detail.
  --density_pctile N: percentile for vmax clipping (default 99.5).
    Lower values (e.g. 95) reveal faint structure at cost of clipping peaks.
"""
import numpy as np, argparse, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

try:
    from scipy.interpolate import RegularGridInterpolator
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--traj',    nargs='+', default=None,
                help='trajectory.csv file(s) to plot')
ap.add_argument('--snaps',   nargs='+', default=None,
                help='snapshot .npz files to make spatial slice plots for')
ap.add_argument('--labels',  nargs='+', default=None,
                help='legend labels for each trajectory (same order as --traj)')
ap.add_argument('--outdir',  default=None,
                help='output directory for plots and npy files '
                     '(default: same as first traj/snap file)')
ap.add_argument('--density_pctile', type=float, default=99.5,
                help='percentile threshold for colorbar vmax (default 99.5)')
ap.add_argument('--no_vacuum_panel', action='store_true',
                help='show only the snapshot column, even when n_vacuum is present')
ap.add_argument('--no_slices', action='store_true',
                help='skip slice plots entirely (useful when only npy export is wanted)')
ap.add_argument('--shared_scale', action='store_true',
                help='vacuum and snapshot panels share the same colorbar scale, '
                     'enabling direct amplitude comparison')
ap.add_argument('--log_scale', action='store_true',
                help='display log10(1+rho) instead of raw rho — useful when the '
                     'field has blown up or has extreme dynamic range')
ap.add_argument('--export_npy', action='store_true', default=True,
                help='export n and v arrays from each snapshot as standalone '
                     '.npy files readable by hopfion_viewer.html (default: on)')
ap.add_argument('--no_export_npy', action='store_true',
                help='disable npy export')
args = ap.parse_args()

if not args.traj and not args.snaps:
    print("Usage: qh_dynamic_viewer.py --traj trajectory.csv [--snaps snap.npz ...]")
    sys.exit(1)

export_npy = args.export_npy and not args.no_export_npy

outdir = args.outdir
if outdir is None and args.traj:
    outdir = os.path.dirname(args.traj[0]) or '.'
elif outdir is None and args.snaps:
    outdir = os.path.dirname(args.snaps[0]) or '.'
else:
    outdir = outdir or '.'
os.makedirs(outdir, exist_ok=True)

phi6 = ((1+5**0.5)/2)**6

# ── Rho_J4 computation (identical to the solvers) ────────────────────
def compute_rho_J4(n, h):
    nx,ny,nz = n[...,0],n[...,1],n[...,2]
    def cd(u,a): return (np.roll(u,-1,a)-np.roll(u,1,a))/(2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    return Fxy**2+Fxz**2+Fyz**2

def maybe_log(rho):
    """Apply log10(1+rho) transform if --log_scale is set."""
    if args.log_scale:
        return np.log10(1.0 + rho)
    return rho

def display_label():
    return "log₁₀(1+ρ_J4)" if args.log_scale else "ρ_J4"

# ── 1. TRAJECTORY PLOTS ───────────────────────────────────────────────
if args.traj:
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Dynamical Integrator — Trajectory', fontsize=14)

    colours = plt.cm.Set2(np.linspace(0, 1, max(len(args.traj), 1)))
    labels  = args.labels or [f'run {i+1}' for i in range(len(args.traj))]

    all_data = []
    for i, traj_path in enumerate(args.traj):
        try:
            data = np.genfromtxt(traj_path, delimiter=',', names=True)
            all_data.append((data, labels[i] if i < len(labels) else f'run {i+1}',
                             colours[i]))
        except Exception as e:
            print(f"Could not load {traj_path}: {e}")

    for data, label, col in all_data:
        t = data['t']
        axes[0,0].plot(t, data['J4'], color=col, label=label)
        axes[0,0].set(title='J4(t)', xlabel='t', ylabel='J4')

        if 'J4_frac_vacuum' in data.dtype.names:
            axes[0,1].plot(t, data['J4_frac_vacuum'], color=col, label=label)
            axes[0,1].axhline(1.0, color='grey', ls='--', lw=0.8, label='vacuum level')
            axes[0,1].set(title='J4/J4_vacuum (1.0 = fully returned)', xlabel='t')

        if 'KE' in data.dtype.names:
            axes[0,2].plot(t, data['KE'], color=col, label=label)
            axes[0,2].set(title='Kinetic Energy KE(t)', xlabel='t', ylabel='KE')

        if 'snapback' in data.dtype.names:
            axes[1,0].plot(t, data['snapback'], color=col, label=label)
            axes[1,0].set(title='Snapback  ||n-n_vac||² (0 = perfect return)',
                          xlabel='t', ylabel='snapback')

        if 'theta_min_kick' in data.dtype.names:
            axes[1,1].plot(t, data['theta_min_kick'], color=col, label=label)
            axes[1,1].set(title='θ_min in kick region (rad)', xlabel='t', ylabel='θ_min')

        if 'K' in data.dtype.names:
            axes[1,2].plot(t, data['K'], color=col, label=label)
            axes[1,2].set(title='K_fb(t)', xlabel='t', ylabel='K')

    for ax in axes.flat:
        if ax.lines: ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    traj_out = os.path.join(outdir, 'trajectory_plot.png')
    plt.savefig(traj_out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {traj_out}")

# ── 2. SPATIAL SLICES + NPY EXPORT for each snapshot ─────────────────
if args.snaps:
    for snap_path in args.snaps:
        try:
            snap = np.load(snap_path, allow_pickle=True)
        except Exception as e:
            print(f"Could not load {snap_path}: {e}")
            continue

        n      = snap['n']
        h      = float(snap.get('h', 0.26))
        step   = int(snap.get('step', 0))
        t_val  = float(snap.get('t', 0.0))
        J4_val = float(snap.get('J4', 0.0))
        N      = n.shape[0]
        tag    = os.path.basename(snap_path).replace('.npz', '')

        # ── NPY EXPORT ────────────────────────────────────────────────
        # Writes n and v as standalone (N,N,N,3) float32 npy files,
        # directly loadable by hopfion_viewer.html (drag-and-drop).
        if export_npy:
            n_npy_path = os.path.join(outdir, f'n_{tag}.npy')
            np.save(n_npy_path, n.astype(np.float32))
            print(f"Exported: {n_npy_path}  shape={n.shape} dtype=float32")

            if 'v' in snap:
                v = snap['v']
                v_npy_path = os.path.join(outdir, f'v_{tag}.npy')
                np.save(v_npy_path, v.astype(np.float32))
                print(f"Exported: {v_npy_path}  shape={v.shape} dtype=float32")

            if 'n_vacuum' in snap:
                vac_npy_path = os.path.join(outdir, f'n_vacuum_{tag}.npy')
                np.save(vac_npy_path, snap['n_vacuum'].astype(np.float32))
                print(f"Exported: {vac_npy_path}  (vacuum reference)")

        if args.no_slices:
            continue

        # ── SLICE PLOTS ───────────────────────────────────────────────
        rho = compute_rho_J4(n, h)

        has_vac = 'n_vacuum' in snap and not args.no_vacuum_panel
        if has_vac:
            n_vac   = snap['n_vacuum']
            rho_vac = compute_rho_J4(n_vac, h)

        cv = h*(np.arange(N) - N//2 + 0.5)

        fig, axes = plt.subplots(3, 2 if has_vac else 1,
                                  figsize=(12 if has_vac else 6, 14))
        if not has_vac:
            axes = np.array(axes).reshape(3, 1)

        scale_note = ""
        if args.shared_scale: scale_note += " [shared scale]"
        if args.log_scale:    scale_note += " [log]"
        fig.suptitle(f'Step {step}  t={t_val:.3f}  J4={J4_val:.3f}{scale_note}',
                     fontsize=13)

        iz = N//2
        ix = N//2
        iy = N//2

        # Apply log transform if requested
        rho_z0 = maybe_log(rho[:, :, iz])
        rho_x0 = maybe_log(rho[ix, :, :])
        rho_y0 = maybe_log(rho[:, iy, :])

        if has_vac:
            vac_z0 = maybe_log(rho_vac[:, :, iz])
            vac_x0 = maybe_log(rho_vac[ix, :, :])
            vac_y0 = maybe_log(rho_vac[:, iy, :])

        def vmax_pair(snap_arr, vac_arr=None):
            """Compute vmax for a snapshot/vacuum pair respecting scale mode."""
            snap_max = np.percentile(snap_arr, args.density_pctile)
            if args.shared_scale and vac_arr is not None:
                vac_max = np.percentile(vac_arr, args.density_pctile)
                combined = max(snap_max, vac_max)
                return combined, combined   # same scale for both
            vac_max = np.percentile(vac_arr, args.density_pctile) \
                      if vac_arr is not None else None
            return snap_max, vac_max

        slices = [
            (rho_z0, vac_z0 if has_vac else None, 'x', 'y',
             f'VACUUM z=0', f'SNAPSHOT z=0 t={t_val:.3f}'),
            (rho_x0, vac_x0 if has_vac else None, 'y', 'z',
             f'VACUUM x=0', f'SNAPSHOT x=0 t={t_val:.3f}'),
            (rho_y0, vac_y0 if has_vac else None, 'x', 'z',
             f'VACUUM y=0 (kick dir)', f'SNAPSHOT y=0 t={t_val:.3f}'),
        ]

        cmap = 'inferno'
        for row, (snap_arr, vac_arr, lx, ly, title_v, title_s) in enumerate(slices):
            snap_vmax, vac_vmax = vmax_pair(snap_arr, vac_arr)

            if has_vac and vac_arr is not None:
                im_v = axes[row,0].contourf(cv, cv, vac_arr.T, levels=40,
                                             cmap=cmap, vmin=0, vmax=vac_vmax)
                axes[row,0].set(title=title_v, xlabel=lx, ylabel=ly, aspect='equal')
                cb = plt.colorbar(im_v, ax=axes[row,0])
                cb.set_label(display_label(), fontsize=8)

                im_s = axes[row,1].contourf(cv, cv, snap_arr.T, levels=40,
                                             cmap=cmap, vmin=0, vmax=snap_vmax)
                axes[row,1].set(title=title_s, xlabel=lx, ylabel=ly, aspect='equal')
                cb = plt.colorbar(im_s, ax=axes[row,1])
                cb.set_label(display_label(), fontsize=8)

                # Warn in title if snapshot looks blown up
                if snap_vmax > 10 * vac_vmax:
                    axes[row,1].set_title(title_s + ' ⚠ SATURATED', color='red')
            else:
                im = axes[row,0].contourf(cv, cv, snap_arr.T, levels=40,
                                           cmap=cmap, vmin=0, vmax=snap_vmax)
                axes[row,0].set(title=title_s, xlabel=lx, ylabel=ly, aspect='equal')
                cb = plt.colorbar(im, ax=axes[row,0])
                cb.set_label(display_label(), fontsize=8)

        plt.tight_layout()
        out = os.path.join(outdir, f'{tag}_slices.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved:    {out}")

print("Done.")
