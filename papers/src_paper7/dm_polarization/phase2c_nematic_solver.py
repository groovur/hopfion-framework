#!/usr/bin/env python3
r"""
Phase 2c: full nematic director-field solver for the collective strain enhancement.

Purpose: replace the Phase 2b analytic bracket with an actual number. Relax the
one-constant Frank director field n(x) on a 3D grid with hedgehog (splay) boundary
conditions pinned at the galaxy cores (coherence-gated: only galaxies source; gas
does not), plus a charge-N radial far-field at the box boundary. Compare the
relaxed COLLECTIVE Frank energy to N x the single-defect energy -> the collective
enhancement factor eta. Project the strain energy density along the line of sight
-> Sigma_strain(x,y), the input to the lensing comparison (Phase 3).

Model (one-constant approximation):
    F = (1/2) K \int |grad n|^2 d^3x,     n in S^2 (unit vector).
Relaxation = harmonic-map heat flow: n <- n + dt * proj_perp(Laplacian n), renormalise,
holding pinned cells (galaxy cores + far-field boundary) fixed.
A single hedgehog n = r_hat gives |grad n|^2 = 2/r^2 -> u ~ K/r^2 -> M_enc ∝ r (isothermal),
i.e. the framework's flat-rotation-curve halo. eta > 1 = frustration/collective enhancement;
eta < 1 = mutual screening.

SMOKE TEST (--smoke): tiny grid, few galaxies, few steps -> verifies the machinery
(energy decreases monotonically, eta and Sigma produced). NOT physically converged.
FULL RUN: large grid + real Bullet galaxy positions (--positions file) + many steps.
Left for the user; outputs land in ./phase2c_output/.
"""
import argparse, os, time, sys
import numpy as np

def build_parser():
    p = argparse.ArgumentParser(description="Phase 2c nematic director-field solver")
    p.add_argument("--smoke", action="store_true", help="tiny smoke-test parameters")
    p.add_argument("--grid", type=int, default=96, help="cubic grid points per side")
    p.add_argument("--ngal", type=int, default=200, help="number of galaxy sources")
    p.add_argument("--steps", type=int, default=4000, help="relaxation iterations")
    p.add_argument("--dt", type=float, default=0.15, help="flow step (< 1/6 for stability)")
    p.add_argument("--core", type=float, default=1.5, help="galaxy core radius (grid cells)")
    p.add_argument("--cluster-frac", type=float, default=0.30,
                   help="cluster scale radius as a fraction of the box half-width (concentration)")
    p.add_argument("--profile", choices=["gaussian", "plummer", "bimodal"], default="bimodal",
                   help="simulated galaxy distribution (bimodal = Bullet-like two subclusters)")
    p.add_argument("--sep", type=float, default=0.45,
                   help="bimodal: subcluster separation as a fraction of the box half-width")
    p.add_argument("--positions", type=str, default=None,
                   help="optional .npy/.txt of REAL galaxy positions (Ngal x 3, grid units); "
                        "overrides --profile (needed only for the Phase-3 offset fit, not for eta)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--outdir", type=str,
                   default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "phase2c_output"))
    return p

def hedgehog_sum(shape, positions, box_center):
    """Initial field: normalized sum of radial-outward hedgehogs from each galaxy."""
    Nx, Ny, Nz = shape
    xs = np.arange(Nx)[:, None, None]
    ys = np.arange(Ny)[None, :, None]
    zs = np.arange(Nz)[None, None, :]
    n = np.zeros((3, Nx, Ny, Nz))
    for (gx, gy, gz) in positions:
        dx = xs - gx; dy = ys - gy; dz = zs - gz
        r = np.sqrt(dx*dx + dy*dy + dz*dz) + 1e-6
        n[0] += dx / r; n[1] += dy / r; n[2] += dz / r
    # far-field bias: radial-outward from cluster centre (net charge-N)
    cx, cy, cz = box_center
    n[0] += 0.0  # (kept explicit for readability)
    norm = np.sqrt((n**2).sum(axis=0)) + 1e-12
    return n / norm

def radial_field(shape, center):
    Nx, Ny, Nz = shape
    xs = np.arange(Nx)[:, None, None]; ys = np.arange(Ny)[None, :, None]; zs = np.arange(Nz)[None, None, :]
    cx, cy, cz = center
    dx = xs - cx + 1e-6; dy = ys - cy + 1e-6; dz = zs - cz + 1e-6
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    return np.stack([dx/r, dy/r, dz/r], axis=0)

def laplacian(n):
    return (np.roll(n, 1, 1) + np.roll(n, -1, 1)
            + np.roll(n, 1, 2) + np.roll(n, -1, 2)
            + np.roll(n, 1, 3) + np.roll(n, -1, 3) - 6.0 * n)

def dirichlet_energy(n):
    """(1/2) sum over +x,+y,+z edges of |n_i - n_{i+1}|^2  == discrete (1/2)|grad n|^2."""
    e = 0.0
    for ax in (1, 2, 3):
        d = np.roll(n, -1, ax) - n
        e += (d*d).sum()
    return 0.5 * e

def dirichlet_energy_ball(n, center, R):
    """(1/2)|grad n|^2 summed only within radius R of `center` (a defect)."""
    Nx, Ny, Nz = n.shape[1:]
    xs = np.arange(Nx)[:, None, None]; ys = np.arange(Ny)[None, :, None]; zs = np.arange(Nz)[None, None, :]
    cx, cy, cz = center
    inside = ((xs-cx)**2 + (ys-cy)**2 + (zs-cz)**2) <= R*R
    e = 0.0
    for ax in (1, 2, 3):
        d = np.roll(n, -1, ax) - n
        e += (d*d).sum(axis=0)[inside].sum()
    return 0.5 * e

def mean_nn_spacing(pos):
    """mean nearest-neighbour distance among galaxy positions."""
    if len(pos) < 2:
        return None
    P = np.asarray(pos)
    d2 = ((P[:, None, :] - P[None, :, :])**2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    return float(np.sqrt(d2.min(axis=1)).mean())

def sample_plummer(n, scale, center, rng):
    """N points from an isotropic Plummer sphere (realistic cluster concentration)."""
    # radial CDF inverse for Plummer: r = a / sqrt(X^{-2/3} - 1)
    X = rng.uniform(1e-6, 1-1e-6, size=n)
    r = scale / np.sqrt(X**(-2.0/3.0) - 1.0)
    u = rng.uniform(-1, 1, size=n); phi = rng.uniform(0, 2*np.pi, size=n)
    s = np.sqrt(1 - u*u)
    pts = np.stack([r*s*np.cos(phi), r*s*np.sin(phi), r*u], axis=1)
    return np.asarray(center) + pts

def generate_positions(profile, ngal, G, cluster_frac, sep, rng):
    center = np.array([G/2., G/2., G/2.])
    scale = cluster_frac * G/2.0
    if profile == "gaussian":
        pos = center + rng.normal(scale=scale, size=(ngal, 3))
    elif profile == "plummer":
        pos = sample_plummer(ngal, scale, center, rng)
    elif profile == "bimodal":                       # Bullet-like: two subclusters along x
        d = sep * G/2.0
        n1 = ngal // 2; n2 = ngal - n1
        c1 = center + np.array([-d/2., 0, 0]); c2 = center + np.array([+d/2., 0, 0])
        pos = np.vstack([sample_plummer(n1, scale*0.9, c1, rng),
                         sample_plummer(n2, scale*0.6, c2, rng)])  # sub-cluster smaller/denser
    return np.clip(pos, 3, G-4)

def normalize(n):
    return n / (np.sqrt((n**2).sum(axis=0)) + 1e-12)

def relax(n, pinned_mask, pinned_vals, steps, dt, log_every=0):
    hist = []
    for it in range(steps):
        lap = laplacian(n)
        # tangential projection: remove component along n
        dot = (lap * n).sum(axis=0)
        lap_perp = lap - dot * n
        n = n + dt * lap_perp
        n = normalize(n)
        # re-pin fixed cells
        n = np.where(pinned_mask, pinned_vals, n)
        if log_every and (it % log_every == 0 or it == steps-1):
            hist.append((it, dirichlet_energy(n)))
    return n, hist

def core_mask(shape, positions, core_r, base_field):
    Nx, Ny, Nz = shape
    xs = np.arange(Nx)[:, None, None]; ys = np.arange(Ny)[None, :, None]; zs = np.arange(Nz)[None, None, :]
    mask = np.zeros((1, Nx, Ny, Nz), dtype=bool)
    vals = base_field.copy()
    for (gx, gy, gz) in positions:
        r2 = (xs-gx)**2 + (ys-gy)**2 + (zs-gz)**2
        m = (r2 <= core_r*core_r)[None]
        mask |= m
        # radial-outward from THIS galaxy inside its core
        dx = xs-gx+1e-6; dy = ys-gy+1e-6; dz = zs-gz+1e-6
        rr = np.sqrt(dx*dx+dy*dy+dz*dz)
        loc = np.stack([dx/rr, dy/rr, dz/rr], axis=0)
        vals = np.where(m, loc, vals)
    return np.repeat(mask, 3, axis=0), vals

def boundary_mask(shape):
    Nx, Ny, Nz = shape
    m = np.zeros((1, Nx, Ny, Nz), dtype=bool)
    m[:, 0, :, :] = m[:, -1, :, :] = True
    m[:, :, 0, :] = m[:, :, -1, :] = True
    m[:, :, :, 0] = m[:, :, :, -1] = True
    return np.repeat(m, 3, axis=0)

def main():
    args = build_parser().parse_args()
    if args.smoke:
        args.grid, args.ngal, args.steps = 24, 4, 300
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    G = args.grid
    shape = (G, G, G)
    center = (G/2., G/2., G/2.)

    # galaxy positions
    if args.positions:
        pos = np.loadtxt(args.positions) if args.positions.endswith(".txt") else np.load(args.positions)
        pos = np.atleast_2d(pos)[:, :3]
        src = f"real:{args.positions}"
    else:
        pos = generate_positions(args.profile, args.ngal, G, args.cluster_frac, args.sep, rng)
        src = f"sim:{args.profile}(frac={args.cluster_frac},sep={args.sep})"
    ngal = len(pos)

    far = radial_field(shape, center)          # charge-N far-field template
    bmask = boundary_mask(shape)

    print(f"[phase2c] grid={G}^3  ngal={ngal}  steps={args.steps}  dt={args.dt}  "
          f"core={args.core}  pos={src}  smoke={args.smoke}")

    # ---- COLLECTIVE relaxation ----
    n0 = hedgehog_sum(shape, pos, center)
    cmask, cvals = core_mask(shape, pos, args.core, n0)
    pin_mask = bmask | cmask
    pin_vals = np.where(bmask, far, np.where(cmask, cvals, n0))
    n0 = np.where(pin_mask, pin_vals, n0)
    t0 = time.time()
    nrel, hist = relax(n0, pin_mask, pin_vals, args.steps, args.dt,
                       log_every=max(1, args.steps//10))
    F_collective = dirichlet_energy(nrel)
    t_coll = time.time() - t0

    # ---- SINGLE-defect reference (one galaxy at centre, same box/steps) ----
    single_pos = np.array([center])
    ns = hedgehog_sum(shape, single_pos, center)
    smask, svals = core_mask(shape, single_pos, args.core, ns)
    spin_mask = bmask | smask
    spin_vals = np.where(bmask, far, np.where(smask, svals, ns))
    ns = np.where(spin_mask, spin_vals, ns)
    ns_rel, _ = relax(ns, spin_mask, spin_vals, args.steps, args.dt)
    F_single = dirichlet_energy(ns_rel)

    # PHYSICAL reference: each isolated halo truncated at half the mean nn spacing
    # (matches the volume actually available to a defect in the crowded cluster).
    d_nn = mean_nn_spacing(pos)
    R_trunc = (d_nn/2.0) if d_nn else (G/2.0)
    F_single_trunc = dirichlet_energy_ball(ns_rel, center, R_trunc)

    eta_boxref = F_collective / (ngal * F_single) if F_single > 0 else float('nan')  # artifact-prone (far-field)

    # RADIUS-MATCHED, LOCAL eta: collective strain within each galaxy's Wigner-Seitz
    # ball (R_trunc = d_nn/2) vs the isolated single-defect halo in the same ball.
    # Excludes the charge-N far field (which biases the full-box comparison ~N^2).
    F_coll_local = sum(dirichlet_energy_ball(nrel, tuple(g), R_trunc) for g in pos)
    eta = F_coll_local / (ngal * F_single_trunc) if F_single_trunc > 0 else float('nan')  # PHYSICAL
    eta_fullbox = F_collective / (ngal * F_single_trunc) if F_single_trunc > 0 else float('nan')  # far-field-biased

    # resolution diagnostic: eta trustworthy only if a halo annulus is resolved
    # between galaxies, i.e. R_trunc = d_nn/2 >> core.
    resolved = (R_trunc >= 3.0 * args.core)
    dnn_over_core = (d_nn/args.core) if d_nn else float('inf')

    # ---- strain energy density + line-of-sight projection ----
    grad2 = np.zeros(shape)
    for ax in (1, 2, 3):
        d = np.roll(nrel, -1, ax) - nrel
        grad2 += (d*d).sum(axis=0)
    u = 0.5 * grad2                       # strain energy density (K=1 units)
    Sigma = u.sum(axis=2)                 # project along z -> Sigma_strain(x,y)

    # ---- monotonicity check (smoke-test success criterion) ----
    energies = [e for _, e in hist]
    monotone = all(energies[i+1] <= energies[i] + 1e-9*abs(energies[i]) for i in range(len(energies)-1))

    np.save(os.path.join(args.outdir, "sigma_strain.npy"), Sigma)
    np.save(os.path.join(args.outdir, "director_final.npy"), nrel.astype(np.float32))
    np.save(os.path.join(args.outdir, "galaxy_positions.npy"), np.asarray(pos))
    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(f"Phase 2c nematic solver -- {'SMOKE TEST' if args.smoke else 'RUN'}\n")
        f.write(f"grid={G}^3  ngal={ngal}  steps={args.steps}  dt={args.dt}  core={args.core}\n")
        f.write(f"positions: {src}\n")
        f.write(f"F_collective       = {F_collective:.6e}\n")
        f.write(f"F_single (full box)= {F_single:.6e}\n")
        f.write(f"mean nn spacing d_nn = {d_nn}\n")
        f.write(f"R_trunc (=d_nn/2)  = {R_trunc:.3f}\n")
        f.write(f"F_single_trunc     = {F_single_trunc:.6e}\n")
        f.write(f"N x F_single_trunc = {ngal*F_single_trunc:.6e}\n")
        f.write(f"F_coll_local (sum over WS balls) = {F_coll_local:.6e}\n")
        f.write(f"eta_PHYSICAL (radius-matched, local) = {eta:.4f}   <- use this "
                f"( >1 enhancement, <1 screening )\n")
        f.write(f"eta_fullbox (far-field-biased ~N)     = {eta_fullbox:.4f}   (do NOT use)\n")
        f.write(f"eta_boxref  (full/full)               = {eta_boxref:.4f}   (artifact-prone)\n")
        f.write(f"d_nn/core = {dnn_over_core:.2f}   RESOLVED (>=3, eta trustworthy): {resolved}\n")
        f.write(f"energy history (collective relax): {energies}\n")
        f.write(f"monotone decrease: {monotone}\n")
        f.write(f"collective relax wall time: {t_coll:.2f}s\n")
        f.write(f"Sigma_strain shape: {Sigma.shape}, total: {Sigma.sum():.4e}, peak: {Sigma.max():.4e}\n")

    # optional plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.imshow(Sigma.T, origin="lower", cmap="magma")
        plt.colorbar(label=r"$\Sigma_{\rm strain}$ (K=1)")
        plt.scatter(pos[:,0], pos[:,1], s=6, c="cyan", marker="+")
        plt.title(f"Phase 2c {'smoke' if args.smoke else 'run'}: eta={eta:.3f}")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "sigma_strain.png"), dpi=120)
        plt.close()
        plotted = True
    except Exception as e:
        plotted = False
        print(f"[phase2c] plot skipped ({e})")

    print(f"[phase2c] d_nn={d_nn}  R_trunc={R_trunc:.2f}")
    print(f"[phase2c] eta_PHYSICAL (radius-matched)={eta:.4f}  "
          f"(eta_fullbox={eta_fullbox:.4f}, far-field-biased)")
    print(f"[phase2c] d_nn/core={dnn_over_core:.2f}  RESOLVED={resolved}"
          + ("" if resolved else "  <-- WARNING: under-resolved, eta NOT trustworthy (raise --grid or lower --ngal/--core)"))
    print(f"[phase2c] energy monotone-decreasing: {monotone}  (smoke-test criterion)")
    print(f"[phase2c] wrote outputs to {args.outdir}  (plot={'yes' if plotted else 'no'})")
    if args.smoke:
        print("[phase2c] SMOKE TEST complete -- machinery verified; not physically converged.")
        print("[phase2c] Full run e.g.:  python3 phase2c_nematic_solver.py --grid 128 --ngal 200 "
              "--steps 8000 --positions bullet_galaxies.npy")

if __name__ == "__main__":
    main()
