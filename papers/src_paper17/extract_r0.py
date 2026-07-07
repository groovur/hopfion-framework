"""
extract_r0.py  --  Measure tube minor radius r0 from Hopfion field snapshots.

File format (from hopfion_viewer.html):
  shape (N, N, N, 3), dtype float32, C-order
  last axis = [nx, ny, nz] interleaved

Usage:
  python extract_r0.py n_000900.npy [n_000500.npy n_001000.npy] --h 0.0375

Measures r0 at the three crossing vertices of T(2,3), where the tube axis
lies exactly in the z=const plane so z is perpendicular to the tube axis.
Uses cubic-spline zero-crossing of nz(z) for sub-grid accuracy (~h^2/12).

Then runs the decisive test: lambda3 * r0 / (2*R0) == phi^2 ?
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import argparse, sys

PHI   = (1 + 5**0.5) / 2
R0    = 3.0
LAM3  = PHI**6        # proved Q_H=3 Bogomolny parameter (Paper XV)
R0_PHI4 = 2*R0/PHI**4  # = 0.8754, the target value

# ── npy loader ────────────────────────────────────────────────────────────────
def load_field(path, h):
    """Load (N,N,N,3) float32 array; return nz as (N,N,N) float64 + metadata."""
    raw = np.load(path)
    if raw.ndim == 4 and raw.shape[3] == 3:
        N = raw.shape[0]
        nz = raw[:, :, :, 2].astype(np.float64)   # nz component
    elif raw.ndim == 3:
        # Scalar field -- assume it's nz directly (less likely but handle it)
        N = raw.shape[0]
        nz = raw.astype(np.float64)
        print(f"  Warning: scalar array detected, treating as nz directly")
    else:
        raise ValueError(f"Unexpected shape {raw.shape}; expected (N,N,N,3)")

    coords = (np.arange(N) - N//2) * h
    print(f"  Loaded {path}: shape={raw.shape}, N={N}, h={h}, "
          f"domain=[{coords[0]:.3f}, {coords[-1]:.3f}]")
    return nz, N, coords

# ── crossing locations ────────────────────────────────────────────────────────
def crossing_sites(r0_guess=0.874):
    """
    Three crossing vertices of T(2,3) at t_n = pi/6 + 2*pi*n/3.
    Gamma(t_n) = (R0*cos(2t_n), R0*sin(2t_n), r0*sin(3t_n))
    where cos(3t_n)=0, sin(3t_n)=+/-1 alternately.
    The tube axis has tang_z = 0 exactly at each crossing (Paper XVI).
    """
    sites = []
    for n in range(3):
        t  = np.pi/6 + 2*np.pi*n/3
        xc = R0 * np.cos(2*t)
        yc = R0 * np.sin(2*t)
        sgn = np.sin(3*t)          # +1 or -1
        zc = r0_guess * sgn        # tube centre height = ±r0
        sites.append(dict(n=n, xc=xc, yc=yc, zc=zc, sgn=float(sgn)))
    return sites

# ── r0 extraction at one crossing ────────────────────────────────────────────
def measure_r0_at_crossing(nz, N, coords, h, site, verbose=True):
    """
    At crossing C_n the tube axis is in the z=const plane, so z is perpendicular
    to the tube.  nz(z) along the vertical line through the tube centre:
      nz ~ -1 at z = zc          (f = pi, tube interior)
      nz =  0 at z = zc ± r0     (f = pi/2, tube boundary)
    We find the zero crossing on the far side from z=0 using cubic spline.
    Returns r0_measured, z_boundary.
    """
    xc, yc, zc, sgn = site['xc'], site['yc'], site['zc'], site['sgn']

    # Nearest grid point in xy, averaged over 3x3 neighbourhood
    ix0 = int(round(xc / h)) + N//2
    iy0 = int(round(yc / h)) + N//2
    ix0 = max(1, min(N-2, ix0))
    iy0 = max(1, min(N-2, iy0))

    nz_line = np.zeros(N)
    cnt = 0
    for di in range(-1, 2):
        for dj in range(-1, 2):
            ii = max(0, min(N-1, ix0+di))
            jj = max(0, min(N-1, iy0+dj))
            nz_line += nz[ii, jj, :]
            cnt += 1
    nz_line /= cnt

    if verbose:
        iz_c = int(round(zc / h)) + N//2
        iz_c = max(0, min(N-1, iz_c))
        print(f"    C{site['n']}: grid ({ix0},{iy0}), "
              f"nz at tube centre (z={zc:.3f}) = {nz_line[iz_c]:.4f}  "
              f"(expect ≈ -1.0)")

    # Cubic spline of nz(z)
    cs = CubicSpline(coords, nz_line)

    # Search for zero crossing on the far side of tube centre from z=0
    # sgn > 0: centre at +zc, search upward (z > zc)
    # sgn < 0: centre at -|zc|, search downward (z < -|zc|)
    iz_c  = int(round(zc / h)) + N//2
    iz_c  = max(1, min(N-2, iz_c))
    step  = 1 if sgn > 0 else -1
    iz_end = min(N-2, iz_c + step*N//4) if sgn > 0 else max(1, iz_c + step*N//4)

    bracket = None
    for iz in range(iz_c, iz_end, step):
        iz_next = iz + step
        if 0 <= iz_next < N and nz_line[iz] * nz_line[iz_next] < 0:
            bracket = (coords[min(iz,iz_next)], coords[max(iz,iz_next)])
            break

    if bracket is None:
        if verbose:
            print(f"    C{site['n']}: *** zero crossing not found "
                  f"(field may not be at Phase 1 exit) ***")
        return None, None

    z_bdy = brentq(cs, bracket[0], bracket[1], xtol=1e-10)
    r0_meas = abs(z_bdy - zc)

    if verbose:
        print(f"    C{site['n']}: z_boundary = {z_bdy:.6f}, "
              f"r0 = {r0_meas:.6f}")

    return r0_meas, z_bdy

# ── J4 and K from the full field ─────────────────────────────────────────────
def compute_K_J4(raw, N, h):
    """Compute K and J4 matching the Paper XVI / solver convention."""
    # raw shape (N,N,N,3), float32
    nx = raw[:,:,:,0].astype(np.float64)
    ny = raw[:,:,:,1].astype(np.float64)
    nz_arr = raw[:,:,:,2].astype(np.float64)

    # Central-difference gradients with periodic wrap
    def grad(f):
        return (
            (np.roll(f,-1,0)-np.roll(f,1,0))/(2*h),
            (np.roll(f,-1,1)-np.roll(f,1,1))/(2*h),
            (np.roll(f,-1,2)-np.roll(f,1,2))/(2*h),
        )

    dnx = grad(nx); dny = grad(ny); dnz = grad(nz_arr)

    # |grad n|^2 = sum over components and spatial directions
    grad_n_sq = sum(d**2 for d in dnx) + sum(d**2 for d in dny) + sum(d**2 for d in dnz)

    # sin^4(f): from n=(sin f cos g, sin f sin g, cos f) -> cos f = nz -> f=arccos(nz)
    # sin^4(f) = (1 - nz^2)^2
    sin4f = (1 - nz_arr**2)**2

    h3 = h**3
    K  = float(np.sum(grad_n_sq) * h3 / 2)
    J4 = float(np.sum(sin4f * grad_n_sq) * h3)
    return K, J4

# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description='Extract r0 from Hopfion snapshots')
    ap.add_argument('files', nargs='+',
                    help='Snapshot .npy files (e.g. n_000500.npy n_000900.npy n_001000.npy)')
    ap.add_argument('--h',  type=float, default=0.0375, help='Grid spacing (default 0.0375)')
    ap.add_argument('--r0-guess', type=float, default=0.874,
                    help='Initial r0 for locating crossings (default 0.874)')
    args = ap.parse_args()

    print(f"\nh = {args.h},  R0 = {R0},  r0_guess = {args.r0_guess}")
    print(f"phi^4 = {PHI**4:.6f},  target r0 = 2*R0/phi^4 = {R0_PHI4:.6f}")
    print(f"decisive-test target: lambda3*r0/(2*R0) = phi^2 = {PHI**2:.6f}\n")
    print("="*65)

    all_results = []
    sites = crossing_sites(r0_guess=args.r0_guess)

    for fpath in args.files:
        print(f"\nFile: {fpath}")
        try:
            raw = np.load(fpath)
        except Exception as e:
            print(f"  Error loading: {e}")
            continue

        if raw.ndim != 4 or raw.shape[3] != 3:
            print(f"  Unexpected shape {raw.shape}. Expected (N,N,N,3).")
            continue

        N = raw.shape[0]
        h = args.h
        coords = (np.arange(N) - N//2) * h
        nz = raw[:,:,:,2].astype(np.float64)
        print(f"  Shape: {raw.shape}, dtype: {raw.dtype}, N={N}, "
              f"domain [{coords[0]:.3f}, {coords[-1]:.3f}]")

        # K and J4
        print("  Computing K and J4...")
        K, J4 = compute_K_J4(raw, N, h)
        print(f"  K = {K:.2f},  J4 = {J4:.2f},  K/J4 = {K/J4:.4f}")

        # r0 extraction
        print("  Measuring r0 at crossing vertices:")
        r0_vals = []
        for site in sites:
            r0_m, z_b = measure_r0_at_crossing(nz, N, coords, h, site)
            if r0_m is not None:
                r0_vals.append(r0_m)

        if not r0_vals:
            print("  No measurements obtained.")
            continue

        r0_mean = np.mean(r0_vals)
        r0_std  = np.std(r0_vals)
        lhs     = LAM3 * r0_mean / (2*R0)
        resid   = lhs - PHI**2

        print(f"\n  r0 = {r0_mean:.6f} ± {r0_std:.6f}  (from {len(r0_vals)} crossings)")
        print(f"  K_min = {K:.1f}")
        print(f"  J4    = {J4:.1f}")
        print()
        print(f"  Decisive test: lambda3 * r0 / (2*R0)")
        print(f"    = phi^6 * {r0_mean:.6f} / 6")
        print(f"    = {lhs:.6f}  (target phi^2 = {PHI**2:.6f})")
        print(f"    residual = {resid:+.6f}  ({resid/PHI**2*100:+.4f}%)")
        print()
        frac_from_phi4 = (r0_mean - R0_PHI4)/R0_PHI4 * 100
        frac_from_874  = (r0_mean - 0.874)/0.874 * 100
        print(f"  |r0 - 2R0/phi^4| / (2R0/phi^4)  = {abs(frac_from_phi4):.4f}%  "
              f"({'✓ <0.1%' if abs(frac_from_phi4)<0.1 else '~ <0.5%' if abs(frac_from_phi4)<0.5 else '✗'})")
        print(f"  |r0 - 0.874|    / 0.874           = {abs(frac_from_874):.4f}%")
        print("="*65)

        all_results.append(dict(file=fpath, K=K, J4=J4, r0=r0_mean, r0_std=r0_std))

    if len(all_results) > 1:
        print("\nSummary across snapshots:")
        print(f"  {'Step':>12}  {'K':>8}  {'J4':>8}  {'r0':>10}  {'±':>8}")
        for r in all_results:
            step = r['file']
            print(f"  {step:>12}  {r['K']:>8.1f}  {r['J4']:>8.1f}  "
                  f"{r['r0']:>10.6f}  {r['r0_std']:>8.6f}")
        r0s = [r['r0'] for r in all_results]
        print(f"\n  Trend: r0 = {r0s[0]:.6f} → {r0s[-1]:.6f} "
              f"(delta = {r0s[-1]-r0s[0]:+.6f})")

if __name__ == '__main__':
    main()
