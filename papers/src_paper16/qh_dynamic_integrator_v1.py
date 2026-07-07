#!/usr/bin/env python3
"""
qh_dynamic_integrator_v1.py  —  time-dependent field integrator
==========================================================================
WHY THIS SCRIPT EXISTS, AND WHAT QUESTION IT IS BUILT TO TEST:

  Every solver this session (v4 through v11) used gradient descent or
  simulated annealing — fundamentally STATIC search methods that can
  only ever move toward lower energy. None of them can represent a
  genuinely DYNAMICAL trajectory: push energy in, watch the field
  stretch toward some maximal deformation, then watch it relax back.
  That is a different kind of physics question, and it needs a
  different kind of tool — a real time-evolution integrator, solving
  Newton's second law for the field rather than minimising an energy.

  The hypothesis under test (stated by the user): Q_H=3 is not an
  independent, stable topological sector with its own static saddle.
  Stability lives in Q_H=2 (the lepton/torus sector, proven the unique
  condensate vacuum — Theorem thm:vacuum, Paper I). What looks like
  "Q_H=3" — the quark sector — may instead be a TRANSIENT EXCITATION:
  energy pushed locally into the Q_H=2 vacuum lets part of the field
  briefly decouple from the condensate's confining suppression,
  moving through the angles where that suppression is weakest, before
  the surrounding (still-suppressed, non-primordial) condensate pushes
  it back toward the Q_H=2 vacuum. This script is the first attempt
  to test that picture directly, by actually evolving the field in
  time rather than searching for a static minimum.

  This is explicitly EXPLORATORY. The mechanisms below (the kick
  scheme, the restoring-force representation of "the confining
  environment") are NEW, not yet validated — there is no existing
  formula in the project's papers for either one. Where Paper XV's
  established results (e.g. the 2T-subset-2I subgroup index, =5) are
  invoked, they are used as falsifiable CHECKS against the simulation
  output, not as inputs baked into the ansatz — see the diagnostics
  section.

THE PHYSICS, PIECE BY PIECE:

  1. INITIAL BACKGROUND: the genuine, converged Q_H=2 saddle (R0=3,
     beta*=0.452, from the project's own fn_hopfion_solver.py output,
     f_fb_beta0_45200.npy), embedded into a 3D Cartesian field via
     n_hat = (sin f cos Phi, sin f sin Phi, cos f),
     Phi = Theta_azimuthal - Theta_poloidal.
     This (1,1)-winding embedding was VERIFIED (not assumed) to give
     the correct Hopf charge by recomputing this project's own J2a,
     J2iso, J4 integrals directly from the saved profile and checking
     them against the established WZW universality predictions:
       J4/J2a = 0.2272 (predicted 0.22721), K_fb/J4 = 14.2465
       (predicted 14.2424), J2iso/J2a = 1.6182 (phi = 1.6180) — all
     matching to 3-4 significant figures, confirming this file is
     genuinely the converged Q_H=2 saddle and that the simple (1,1)
     embedding (not a higher-winding variant) is the correct one.

  2. THE ENERGY: E = K_fb + phi^6 * J4 (the project's actual physical
     energy, Corollary cor:lambda's lambda=phi^6 fixed coefficient) —
     NOT the K_fb*J4 product trick used in v6-v11, which was a SEARCH
     device to avoid Derrick collapse during static gradient descent,
     not the real physical energy. Genuine dynamics should evolve
     under the TRUE energy.

  3. THE PERTURBATION ("pin pushed into the fabric"): a localised
     initial VELOCITY (not a position change — pushing energy in is a
     kinetic, not a static, act), concentrated in a small region on
     the torus, pointing toward the LOWEST-SUPPRESSION direction. This
     project's own suppression formula S(theta)=sin^4(theta)/phi^6
     (Paper I, eq:suppression) is smallest at theta->0,pi (the poles
     of the preferred axis) — so the kick rotates n locally toward
     (0,0,+-1). The kick's magnitude is the free "energy input" knob.

  4. THE CONFINING ENVIRONMENT ("the surrounding condensate is still
     suppressed and pushes back"): a restoring force,
       F_restore(x,t) = -k_restore * tangent[(n(x,t) - n_vacuum(x))]
     applied EVERYWHERE, pulling the field back toward the fixed,
     pre-computed Q_H=2 vacuum reference. This is a modelling CHOICE,
     not a derived result — there is no existing formula for this in
     Paper XV. It is the most direct, literal translation of the
     physical picture into a force term, and k_restore is left as an
     explicit, scannable parameter.

  5. THE INTEGRATOR: a RATTLE-style constrained leapfrog (the standard
     symplectic method for Hamiltonian dynamics confined to a
     manifold). Verified on a toy S^2-constrained test particle before
     use here: the unit-norm constraint holds to machine precision
     throughout, and energy drift over a FIXED total integration time
     shrinks roughly linearly as dt shrinks (0.20 -> 0.022 -> 0.0022
     for dt -> dt/10 -> dt/100), the standard signature of a
     well-behaved (not buggy) leapfrog splitting error — NOT a sign
     that the method is wrong. This matters for reading the energy
     diagnostic below: slow, dt-dependent drift over MANY periods is
     normal; drift that does NOT shrink with dt would indicate a bug.

  6. DIAGNOSTICS, mapped directly onto the user's specific questions:
     - theta_min(t): minimum polar angle reached near the kick site
       (how close to the lowest-suppression direction the field gets).
     - J4(t), K_fb(t), E(t): global energy/charge trajectory — look
       for a clean RISE to a turning point then a FALL (the
       "stretch then snap back" signature) versus monotonic decay or
       runaway growth.
     - "snap-back" metric: integrated ||n(t)-n_vacuum||^2 — does the
       field return close to its starting configuration?
     - energy_at_turning_point: the local energy near the kick site
       at the moment J4(t) peaks — the user's proposed "quark mass"
       proxy. Run at several --kick_strength values and compare this
       quantity across runs: a consistent, well-behaved relationship
       would be suggestive; an erratic one would undercut the "mass"
       interpretation.
     - GEOMETRY CHECK (falsifiable, not built-in): index([2I:2T])=5 is
       an EXISTING, proven number (Paper XV, Theorem thm:subgroup).
       This script prints ratios like (sweeps to turning point)/
       (sweeps to return), or winding-like counts during the stretch
       phase, alongside this number — not because the simulation is
       built to reproduce it, but so any resemblance (or lack of one)
       is visible and checkable, not asserted.

Install:  pip install torch --index-url https://download.pytorch.org/whl/cpu

USAGE
-----
  python qh_dynamic_integrator_v1.py \\
      --kick_strength 0.5 --k_restore 0.3 --dt 0.001 --steps 50000 \\
      --outdir dyn_test1

  # Scan kick strength to test the "mass" proxy for consistency:
  for k in 0.2 0.5 1.0 2.0; do
    python qh_dynamic_integrator_v1.py --kick_strength $k --outdir dyn_kick$k
  done

OUTPUT (in --outdir)
---------------------
  log.txt                  full run log
  snapshots/n_t{step}.npz  periodic field snapshots (for the extended viewer)
  trajectory.csv           per-print-interval diagnostics
  report.txt               summary, including the turning-point energy
                            and the geometry-ratio check
"""
import numpy as np, time, argparse, os, sys
from scipy.spatial import KDTree
from scipy.interpolate import RegularGridInterpolator
try:
    import torch
except ImportError:
    print("pip install torch --index-url https://download.pytorch.org/whl/cpu")
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('--N',            type=int,   default=64)
ap.add_argument('--h',            type=float, default=0.26)
ap.add_argument('--R0',           type=float, default=3.0)
ap.add_argument('--vacuum_profile', type=str, default='./f_fb_beta0_45200.npy')
ap.add_argument('--vacuum_h',     type=float, default=0.12,
                help='Grid spacing used to GENERATE the saved axisymmetric profile (fn_hopfion_solver.py LARGE grid default)')
ap.add_argument('--mu2',          type=float, default=1.0,
                help='Kinetic inertia coefficient. No established physical value exists for this (the project has been purely static) — 1.0 is a normalisation choice, equivalent to a choice of time units.')
ap.add_argument('--kick_strength', type=float, default=0.5,
                help='Magnitude of the initial velocity kick — the free "energy input" parameter.')
ap.add_argument('--kick_sigma',   type=float, default=0.8,
                help='Gaussian envelope radius (real 3D distance) of the kick region')
ap.add_argument('--kick_r0_offset', type=float, default=0.0,
                help='Place the kick centre at poloidal angle (radians) around the tube, 0 = outer equator')
ap.add_argument('--k_restore',    type=float, default=0.3,
                help='Restoring-force strength representing the confining environment. NOT a derived value — a modelling choice, scan this.')
ap.add_argument('--dt',           type=float, default=0.001)
ap.add_argument('--steps',        type=int,   default=50000)
ap.add_argument('--print_every',  type=int,   default=200)
ap.add_argument('--snapshot_every', type=int, default=2000)
ap.add_argument('--outdir',       type=str,   default='.')
ap.add_argument('--device',       type=str,   default='cpu')
args = ap.parse_args()

phi = (1+5**0.5)/2; phi6 = phi**6; MU = 3.0-phi
N, h = args.N, args.h
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
os.makedirs(os.path.join(args.outdir, 'snapshots'), exist_ok=True)
log_path = os.path.join(args.outdir, 'log.txt')
log = open(log_path, 'w')
def LOG(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    print(s, file=log); log.flush()

LOG("="*70)
LOG("  Q_H=2 -> transient stretch -> Q_H=2 dynamical integrator  v1")
LOG("  Testing: is 'Q_H=3' a transient excitation of the Q_H=2 vacuum,")
LOG("  not an independent stable sector?")
LOG("="*70)
LOG(f"  phi^6={phi6:.6f}  mu2={args.mu2}  kick_strength={args.kick_strength}")
LOG(f"  k_restore={args.k_restore}  dt={args.dt}  steps={args.steps}")
LOG(f"  Grid: {N}^3  h={h}  box=[{-N*h/2:.2f},{N*h/2:.2f}]")

# ── Grid ──────────────────────────────────────────────────────────────
cv = h*(np.arange(N) - N//2 + 0.5)
X, Y, Zc = np.meshgrid(cv, cv, cv, indexing='ij')
Rcyl = np.sqrt(X**2+Y**2)
ThetaAz = np.arctan2(Y, X)
ThetaPol = np.arctan2(Zc, Rcyl - args.R0)

# ── Build the Q_H=2 vacuum background from the verified saddle profile ──
LOG("\nLoading and embedding the verified Q_H=2 vacuum background...")
f_axisym = np.load(args.vacuum_profile)
Nr_src, Nz_src = f_axisym.shape
r_src = args.vacuum_h*(np.arange(Nr_src)+0.5)
z_src = args.vacuum_h*np.arange(Nz_src)
# The source grid only covers z>=0 (axisymmetric solver assumes z reflection
# symmetry — apply_bc forces f=0 at the z=0 boundary edge in the SOURCE code
# only at the OUTER r boundary, not at z=0 itself, so check: actually the
# source's z grid starts at z=0 and goes up — mirror it for z<0 here since
# the physical torus is symmetric under z -> -z).
f_full = np.concatenate([f_axisym[:, ::-1], f_axisym[:, 1:]], axis=1)
z_full = np.concatenate([-z_src[::-1], z_src[1:]])
interp = RegularGridInterpolator((r_src, z_full), f_full, bounds_error=False, fill_value=0.0)

query_pts = np.stack([Rcyl.ravel(), Zc.ravel()], axis=-1)
f3 = interp(query_pts).reshape(N, N, N)
f3 = np.clip(f3, 0.0, np.pi)

Phi_vac = ThetaAz - ThetaPol
n_vacuum_np = np.stack([
    np.sin(f3)*np.cos(Phi_vac),
    np.sin(f3)*np.sin(Phi_vac),
    np.cos(f3)
], axis=-1).astype(np.float32)
n_vacuum_np /= np.linalg.norm(n_vacuum_np, axis=-1, keepdims=True).clip(1e-10)

n_vacuum = torch.tensor(n_vacuum_np, dtype=torch.float32, device=dev)
LOG(f"  Vacuum background built. |n| check: min={np.linalg.norm(n_vacuum_np,axis=-1).min():.6f} "
    f"max={np.linalg.norm(n_vacuum_np,axis=-1).max():.6f}")

# ── Energy and tangent-projected force (autograd, same machinery as v6-v11) ──
def compute_K_J4(n):
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a  = (s4 * g2).sum() * h**3
    J2iso = g2.sum() * h**3
    K = J2a + MU * J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4 = rho_J4.sum() * h**3
    return K, J2a, J4, rho_J4


def spatial_force(n):
    """Tangent-projected -dE/dn for E = K_fb + phi^6*J4 (the project's
    real physical energy, NOT the v6-v11 product-trick surrogate)."""
    n_req = n.detach().requires_grad_(True)
    K, J2a, J4, _ = compute_K_J4(n_req)
    E = K + phi6*J4
    E.backward()
    g_raw = -n_req.grad.detach()
    g_tangent = g_raw - (g_raw*n).sum(-1, keepdim=True)*n
    return g_tangent, K.item(), J4.item()


def measure_tube_radius(f_axisym_2d, r_grid, R0):
    """Measure the Q_H=2 torus's actual tube half-width directly from the
    loaded profile (the f=pi/2 crossing along z=0, moving outward from
    R0), rather than assuming a value. An earlier draft of this script
    reused the TREFOIL project's tube radius (1/C*_2≈0.40) here, which
    is wrong for this different (Q_H=2, not Q_H=3) sector — measuring
    directly on this profile gives ≈1.98, nearly 5x larger."""
    f_at_z0 = f_axisym_2d[:, 0]
    idx_R0 = np.argmin(np.abs(r_grid - R0))
    for i in range(idx_R0, len(r_grid)):
        if f_at_z0[i] < np.pi/2:
            return r_grid[i] - R0
    return 1.0  # fallback, should not normally trigger


TUBE_RADIUS = measure_tube_radius(f_axisym, r_src, args.R0)
LOG(f"  Measured Q_H=2 tube half-width (f=pi/2 crossing): {TUBE_RADIUS:.4f}")


# ── The perturbation: localised velocity kick toward lowest-suppression
#    angle (the poles, theta->0,pi, where S(theta)=sin^4(theta)/phi^6 is
#    smallest) ──────────────────────────────────────────────────────────
def build_kick(n_vac_np, R0, tube_radius, kick_sigma, kick_r0_offset, kick_strength):
    # Kick centre: outer surface of the tube at the given poloidal offset
    r0_kick = R0 + tube_radius
    cx = r0_kick*np.cos(0.0)  # place at azimuthal angle 0 WLOG
    cy = r0_kick*np.sin(0.0)
    cz = tube_radius*np.sin(kick_r0_offset)
    centre = np.array([cx, cy, cz])
    dist = np.sqrt((X-centre[0])**2 + (Y-centre[1])**2 + (Zc-centre[2])**2)
    envelope = np.exp(-(dist**2)/(2*kick_sigma**2))

    z_axis = np.zeros_like(n_vac_np); z_axis[...,2] = 1.0
    # Push toward whichever pole n_vacuum is CLOSER to locally (so the
    # kick is a genuine small perturbation, not a near-180-degree flip)
    sign = np.sign(n_vac_np[...,2:3])
    sign[sign==0] = 1.0
    target_pole = z_axis * sign
    raw = target_pole - n_vac_np
    raw_tangent = raw - (raw*n_vac_np).sum(-1,keepdims=True)*n_vac_np
    raw_norm = np.linalg.norm(raw_tangent, axis=-1, keepdims=True).clip(1e-10)
    direction = raw_tangent / raw_norm

    v = kick_strength * envelope[...,None] * direction
    return v.astype(np.float32), centre, envelope


v_init_np, kick_centre, envelope = build_kick(
    n_vacuum_np, args.R0, TUBE_RADIUS, args.kick_sigma, args.kick_r0_offset, args.kick_strength)
v = torch.tensor(v_init_np, dtype=torch.float32, device=dev)
# Ensure tangency exactly (it should already be, by construction, but
# enforce exactly to remove any float error before the run starts)
v = v - (v*n_vacuum).sum(-1, keepdim=True)*n_vacuum

KE_injected = 0.5*args.mu2*(v**2).sum().item()*h**3
LOG(f"\n  Kick centre: {kick_centre}  sigma={args.kick_sigma}")
LOG(f"  Kinetic energy injected: {KE_injected:.6f}")

# Region mask for "near the kick site" diagnostics (envelope > 0.05)
kick_mask_np = envelope > 0.05
kick_mask = torch.tensor(kick_mask_np, device=dev)
n_kick_sites = int(kick_mask_np.sum())
LOG(f"  Kick affects {n_kick_sites} grid sites (envelope>0.05)")

# ── RATTLE-style constrained leapfrog, with the restoring force added ──
def total_force(n, n_vac, k_restore):
    f_spatial, K_val, J4_val = spatial_force(n)
    restore_raw = -k_restore * (n - n_vac)
    restore_tangent = restore_raw - (restore_raw*n).sum(-1, keepdim=True)*n
    return f_spatial + restore_tangent, K_val, J4_val


def geodesic_step(n, v, dt):
    speed = v.norm(dim=-1, keepdim=True).clamp(1e-12)
    theta = speed*dt
    direction = v/speed
    return n*torch.cos(theta) + direction*torch.sin(theta)


n = n_vacuum.clone()
a, K0, J4_0 = total_force(n, n_vacuum, args.k_restore)
a = a/args.mu2

E0 = K0 + phi6*J4_0 + KE_injected
LOG(f"\n  Initial: K={K0:.4f}  J4={J4_0:.4f}  E_total(with KE)={E0:.4f}")
LOG(f"\n{'step':>8} {'t':>10} {'J4':>10} {'K':>10} {'E_tot':>12} "
    f"{'KE':>10} {'theta_min(kick)':>16} {'snapback':>10}")

csv_path = os.path.join(args.outdir, 'trajectory.csv')
with open(csv_path, 'w') as f:
    f.write("step,t,J4,K,E_total,KE,theta_min_kick,snapback,KJ4\n")

J4_history = []
best_J4_step = 0
best_J4_val = J4_0

for step in range(1, args.steps+1):
    t = step*args.dt
    v_half = v + 0.5*args.dt*a
    n_new = geodesic_step(n, v_half, args.dt)
    n_new = n_new / n_new.norm(dim=-1, keepdim=True).clamp(1e-10)  # safety renorm
    a_new, K_new, J4_new = total_force(n_new, n_vacuum, args.k_restore)
    a_new = a_new/args.mu2
    v_un = v_half + 0.5*args.dt*a_new
    v_new = v_un - (v_un*n_new).sum(-1, keepdim=True)*n_new

    n, v, a = n_new, v_new, a_new
    J4_history.append(J4_new)
    if J4_new > best_J4_val:
        best_J4_val = J4_new; best_J4_step = step

    if step % args.print_every == 0 or step == 1:
        with torch.no_grad():
            KE_now = 0.5*args.mu2*(v**2).sum().item()*h**3
            E_tot = K_new + phi6*J4_new + KE_now
            nz_kick = n[...,2][kick_mask]
            theta_kick = torch.acos(nz_kick.clamp(-1,1))
            dist_to_nearest_pole = torch.minimum(theta_kick, np.pi - theta_kick)
            theta_min_signed = dist_to_nearest_pole.min().item()
            snapback = ((n - n_vacuum)**2).sum().item()*h**3
            KJ4 = K_new/J4_new if J4_new > 1e-12 else float('nan')
        LOG(f"{step:>8} {t:>10.4f} {J4_new:>10.4f} {K_new:>10.4f} {E_tot:>12.4f} "
            f"{KE_now:>10.6f} {theta_min_signed:>16.4f} {snapback:>10.4f}")
        with open(csv_path, 'a') as fcsv:
            fcsv.write(f"{step},{t:.6f},{J4_new:.6f},{K_new:.6f},{E_tot:.6f},"
                       f"{KE_now:.6f},{theta_min_signed:.6f},{snapback:.6f},{KJ4:.6f}\n")

    if step % args.snapshot_every == 0 or step == 1:
        np.savez(os.path.join(args.outdir, 'snapshots', f'n_t{step:08d}.npz'),
                  n=n.detach().cpu().numpy(), v=v.detach().cpu().numpy(),
                  step=step, t=t, J4=J4_new, K=K_new,
                  N=N, h=h, R0=args.R0, kick_strength=args.kick_strength,
                  k_restore=args.k_restore)

# ── Report ────────────────────────────────────────────────────────────
LOG(f"\n{'='*70}")
LOG(f"  RUN SUMMARY")
LOG(f"{'='*70}")
LOG(f"  J4 peaked at step {best_J4_step} (t={best_J4_step*args.dt:.4f}), value={best_J4_val:.4f}")
LOG(f"  J4 initial value: {J4_0:.4f}  (peak/initial ratio: {best_J4_val/J4_0:.4f})")

# Geometry check (falsifiable, not built in): does any natural ratio
# from this run land near the proven 2T-in-2I index, 5?
if best_J4_step > 0:
    steps_to_peak = best_J4_step
    steps_total = args.steps
    steps_after_peak = steps_total - best_J4_step
    ratio_total_to_peak = steps_total/steps_to_peak if steps_to_peak>0 else float('nan')
    LOG(f"\n  Geometry check (Paper XV, Theorem thm:subgroup: index[2I:2T]=5):")
    LOG(f"    steps_total/steps_to_peak = {ratio_total_to_peak:.4f}  (vs 5)")
    LOG(f"    [This is a CHECK, not a built-in assumption — a coincidence")
    LOG(f"     or a clean miss are both informative; nothing in the kick or")
    LOG(f"     restoring-force setup was tuned to produce this number.]")

with open(os.path.join(args.outdir, 'report.txt'), 'w') as f:
    f.write(f"Dynamical Q_H=2 -> stretch -> Q_H=2 integrator, run summary\n")
    f.write(f"kick_strength={args.kick_strength}  k_restore={args.k_restore}  "
            f"dt={args.dt}  steps={args.steps}\n\n")
    f.write(f"KE injected: {KE_injected:.6f}\n")
    f.write(f"J4 initial: {J4_0:.6f}\n")
    f.write(f"J4 peak: {best_J4_val:.6f} at step {best_J4_step}\n")
    f.write(f"J4 peak/initial ratio: {best_J4_val/J4_0:.6f}\n")

log.close()
print(f"\nDone. Trajectory: {csv_path}\nSnapshots: {args.outdir}/snapshots/")
