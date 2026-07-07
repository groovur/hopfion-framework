#!/usr/bin/env python3
"""
qh_dynamic_integrator_v2.py  —  unified parameterised integrator
==========================================================================
DESIGN GOAL: a single script that can reproduce v1 exactly and also
explore every extension, controlled by explicit flags. Every combination
composes freely:

  --vacuum on/off      whether to load the Q_H=2 background at all
  --k_restore_mode     flat (v1) | suppression_weighted (v2) | off
  --gamma              0 = pure RATTLE (v1), >0 = Langevin damping
  --absorb_width       0 = off (v1), >0 = sponge boundary layer
  --resume             continue any existing npz snapshot

EXACT v1 REPRODUCTION (all three variants):
  # v1 kick=0.5 kr=0.3
  python qh_dynamic_integrator_v2.py \
      --gamma 0 --absorb_width 0 --k_restore 0.3 --k_restore_mode flat \
      --kick_strength 0.5 --steps 50000 --outdir repro_v1_k03

  # v1 kick=0.5 kr=1.0
  python qh_dynamic_integrator_v2.py \
      --gamma 0 --absorb_width 0 --k_restore 1.0 --k_restore_mode flat \
      --kick_strength 0.5 --steps 100000 --outdir repro_v1_k10

  # Pure RATTLE, no restoring force, just the energy functional:
  python qh_dynamic_integrator_v2.py \
      --gamma 0 --absorb_width 0 --k_restore_mode off \
      --kick_strength 0.5 --steps 50000 --outdir rattle_pure

EXPLORATION EXAMPLES:
  # Suppression-weighted restore, no damping (new physics, conservative):
  python qh_dynamic_integrator_v2.py \
      --gamma 0 --absorb_width 0 --k_restore 0.3 --k_restore_mode suppression_weighted \
      --kick_strength 0.5 --steps 80000 --outdir sw_nodamp

  # Suppression-weighted restore + Langevin damping outside kick only:
  python qh_dynamic_integrator_v2.py \
      --gamma 0.1 --gamma_profile outside_kick --absorb_width 0 \
      --k_restore 0.3 --k_restore_mode suppression_weighted \
      --kick_strength 0.5 --steps 80000 --outdir sw_damp01

  # Full v2: suppression-weighted + damping + sponge:
  python qh_dynamic_integrator_v2.py \
      --gamma 0.1 --gamma_profile outside_kick --absorb_width 0.15 \
      --k_restore 0.3 --k_restore_mode suppression_weighted \
      --kick_strength 0.5 --steps 80000 --outdir sw_full

  # Resume a v1 snapshot and continue with damping switched on:
  python qh_dynamic_integrator_v2.py \
      --resume dyn_v1/snapshots/n_t00050000.npz \
      --gamma 0.1 --gamma_profile outside_kick --k_restore_mode suppression_weighted \
      --steps 50000 --outdir resumed_with_damp

PARAMETERS
----------
  --vacuum_profile   path to f_fb_beta0_45200.npy (the Q_H=2 saddle).
                     If the file does not exist, the script exits with a
                     clear message rather than silently using a wrong field.

  --k_restore_mode   flat               : same as v1 — constant k_restore
                                          in ALL directions, pins field toward
                                          n_vacuum regardless of theta.
                     suppression_weighted: k_restore * sin^4(theta) pointwise —
                                          zero resistance at poles, full at
                                          equator. Follows S(theta)=sin^4/phi^6.
                     off                 : no restoring force. Pure energy
                                          dynamics from E=K_fb+phi^6*J4.

  --gamma            Langevin damping rate (1/time). 0 = pure conservative
                     RATTLE (v1-style). >0 = energy dissipates at rate gamma.
                     Half-life of KE = ln(2)/gamma time units.

  --gamma_profile    uniform             : gamma everywhere.
                     outside_kick        : 0 inside kick envelope, gamma outside.
                     boundary_only       : 0 everywhere except sponge layer.

  --absorb_width     Fraction of box half-width used as absorbing sponge
                     (0 = off = v1-style periodic with no absorption).
                     Typical: 0.10–0.20. Only meaningful if gamma_boundary > 0.

  --gamma_boundary   Damping strength specifically in the sponge layer
                     (added on top of gamma_profile, independently).

WHAT TO WATCH IN OUTPUT
-----------------------
  J4/J4_vac   (new column) — 1.0 at start, falls during stretch,
               returns toward 1.0 if dissipation is working.
               In v1 this stabilised at ~0.22–0.44.
  snapback     ||n(t) - n_vacuum||^2 — 0 means full return to vacuum.
               In v1 this peaked at ~140 and stayed there.
  KE           should decay toward 0 when gamma>0. Timescale ~ 1/gamma.

Install:  pip install torch --index-url https://download.pytorch.org/whl/cpu
"""
import numpy as np, argparse, os, sys
from scipy.interpolate import RegularGridInterpolator
try:
    import torch
except ImportError:
    print("pip install torch --index-url https://download.pytorch.org/whl/cpu")
    sys.exit(1)

# ── Arguments ─────────────────────────────────────────────────────────
ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
# Grid
ap.add_argument('--N',              type=int,   default=64)
ap.add_argument('--h',              type=float, default=0.26)
ap.add_argument('--R0',             type=float, default=3.0)
# Vacuum background
ap.add_argument('--vacuum_profile', type=str,
                default='/mnt/user-data/uploads/f_fb_beta0_45200.npy',
                help='Path to the saved Q_H=2 axisymmetric profile (.npy). '
                     'Required for vacuum reference and suppression-weighted restore.')
ap.add_argument('--vacuum_h',       type=float, default=0.12,
                help='Grid spacing used to generate the vacuum profile file.')
# Physics
ap.add_argument('--mu2',            type=float, default=1.0,
                help='Kinetic inertia. 1.0 is a normalisation choice (time units).')
ap.add_argument('--k_restore',      type=float, default=0.3,
                help='Restoring force strength (ignored when --k_restore_mode off).')
ap.add_argument('--k_restore_mode', type=str,   default='flat',
                choices=['flat', 'suppression_weighted', 'off'],
                help='flat=v1 behaviour; suppression_weighted=sin^4(theta)*k_restore; '
                     'off=no restoring force at all.')
# Kick
ap.add_argument('--kick_axis',      type=str,   default='pole',
                help='Direction the kick pushes n toward. Options: '
                     '"pole" (toward nearest ±z pole, original behaviour); '
                     '"x", "y", "z" (fixed Cartesian axis); '
                     '"a,b,c" (explicit unit vector, e.g. "0.707,0,0.707"); '
                     '"random" (random unit vector, seed controlled by --kick_seed).')
def _float_or_random(s):
    """Allow --kick_azimuth to be a float (radians) or the string 'random'."""
    if s is None:
        return 0.0
    if s.strip().lower() == 'random':
        return 'random'
    return float(s)

ap.add_argument('--kick_azimuth',   type=_float_or_random, default=0.0,
                help='Azimuthal angle (radians) of the kick centre around the '
                     'torus ring. 0=original position at (R0+r_tube, 0, 0). '
                     '"random"=draw uniformly from [0, 2π) using --kick_seed.')
def _int_or_none(s):
    """Allow --kick_seed to be omitted (uses auto-generated seed) or an int."""
    if s is None or s.lower() in ('none', 'auto', ''):
        return None
    return int(s)

ap.add_argument('--kick_seed',      type=_int_or_none, default=None,
                help='Integer seed for the kick RNG. If omitted (default), a seed is '
                     'generated automatically from the OS entropy and saved in the '
                     'snapshot, so every run is always reproducible from its output. '
                     'Pass an explicit integer to reproduce a specific kick.')
ap.add_argument('--kick_strength',  type=float, default=0.5)
ap.add_argument('--kick_sigma',     type=float, default=0.8)
ap.add_argument('--kick_r0_offset', type=float, default=0.0,
                help='Poloidal offset of kick centre on tube surface (radians).')
# Damping
ap.add_argument('--gamma',          type=float, default=0.0,
                help='Langevin damping rate. 0=pure conservative RATTLE (v1). '
                     '>0=energy dissipation at rate gamma.')
ap.add_argument('--gamma_profile',  type=str,   default='uniform',
                choices=['uniform', 'outside_kick', 'boundary_only'],
                help='Spatial profile of gamma damping field.')
ap.add_argument('--gamma_boundary', type=float, default=1.0,
                help='Extra damping strength in the absorbing sponge layer.')
ap.add_argument('--absorb_width',   type=float, default=0.0,
                help='Sponge layer thickness as fraction of box half-width. '
                     '0=off (v1-style). Typical: 0.10-0.20.')
# Icosahedral condensate correction
ap.add_argument('--ico_epsilon',    type=float, default=0.0,
                help='Strength of the icosahedral correction to the suppression '
                     'weight. 0=off (default, recovers exact suppression_weighted '
                     'behaviour). >0 adds a 5-fold azimuthal modulation to the '
                     'suppression landscape, representing the discrete icosahedral '
                     'geometry of the 2I condensate rather than its continuous '
                     'azimuthal average.\n'
                     'Only active when --k_restore_mode suppression_weighted.\n'
                     'The correction is delta_S_ico(n) = nz*(nx^5-10*nx^3*ny^2+'
                     '5*nx*ny^4) = Re[(nx+i*ny)^5*nz], the unique A_5-invariant '
                     'degree-6 polynomial on S^2 with 5-fold azimuthal symmetry. '
                     'It equals 0 at the poles, +0.256 at all 10 icosahedral '
                     'non-polar vertices (verified numerically), and averages to '
                     'zero over all azimuthal angles (so macroscopic observables '
                     'like sin_W are unchanged). Safe range: 0 <= epsilon < 2.')
# Run control
ap.add_argument('--dt',             type=float, default=0.001)
ap.add_argument('--steps',          type=int,   default=50000)
ap.add_argument('--print_every',    type=int,   default=200)
ap.add_argument('--snapshot_every', type=int,   default=2000)
ap.add_argument('--resume',         type=str,   default=None,
                help='Path to snapshot .npz to resume from. Kick is NOT re-applied. '
                     'New gamma/k_restore_mode/dt from THIS call are used.')
ap.add_argument('--outdir',         type=str,   default='.')
ap.add_argument('--device',         type=str,   default='cpu')
args = ap.parse_args()

phi = (1+5**0.5)/2; phi6 = phi**6; MU = 3.0-phi
N, h = args.N, args.h
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
os.makedirs(os.path.join(args.outdir, 'snapshots'), exist_ok=True)
log_path = os.path.join(args.outdir, 'log.txt')
log = open(log_path, 'a')

def LOG(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    print(s, file=log); log.flush()

# ── Determine active features for this run ───────────────────────────
USE_DAMPING   = args.gamma > 0 or args.absorb_width > 0
USE_RESTORE   = args.k_restore_mode != 'off'
USE_VACUUM    = USE_RESTORE or args.k_restore_mode != 'off'

# v1-exact flag: print clearly so output is self-documenting
is_v1_exact = (args.gamma == 0 and args.absorb_width == 0
               and args.k_restore_mode == 'flat' and args.resume is None)

LOG("="*70)
LOG("  Q_H=2 dynamical integrator  v2  (unified parameterised)")
LOG("="*70)
LOG(f"  N={N}  h={h}  R0={args.R0}  dt={args.dt}  steps={args.steps}")
LOG(f"  mu2={args.mu2}  kick_strength={args.kick_strength}  kick_sigma={args.kick_sigma}")
LOG(f"  kick_axis={args.kick_axis}  kick_azimuth={args.kick_azimuth}  "
    f"kick_seed={args.kick_seed}")
LOG(f"  k_restore_mode={args.k_restore_mode}  k_restore={args.k_restore}")
LOG(f"  gamma={args.gamma}  gamma_profile={args.gamma_profile}")
LOG(f"  absorb_width={args.absorb_width}  gamma_boundary={args.gamma_boundary}")
LOG(f"  ico_epsilon={args.ico_epsilon}"
    + (" [active: suppression_weighted + icosahedral correction]"
       if args.ico_epsilon != 0 and args.k_restore_mode == 'suppression_weighted'
       else " [inert: only active with --k_restore_mode suppression_weighted]"
       if args.ico_epsilon != 0 else ""))
if is_v1_exact:
    LOG("  ** v1-EXACT mode: gamma=0, no sponge, flat k_restore **")
if args.resume:
    LOG(f"  RESUMING from: {args.resume}")

# ── Grid ──────────────────────────────────────────────────────────────
cv = h*(np.arange(N) - N//2 + 0.5)
X, Y, Zc = np.meshgrid(cv, cv, cv, indexing='ij')
Rcyl = np.sqrt(X**2+Y**2)
ThetaAz  = np.arctan2(Y, X)
ThetaPol = np.arctan2(Zc, Rcyl - args.R0)

# ── Absorbing sponge ──────────────────────────────────────────────────
box_half = N*h/2
if args.absorb_width > 0:
    dist_from_wall = box_half - np.abs(
        np.stack([X, Y, Zc], axis=-1)).max(axis=-1)
    sponge_depth   = args.absorb_width * box_half
    sponge_weight  = np.clip(1.0 - dist_from_wall /
                              np.maximum(sponge_depth, 1e-6), 0.0, 1.0)
    gamma_boundary_field = args.gamma_boundary * sponge_weight**2
else:
    gamma_boundary_field = np.zeros((N,N,N), dtype=np.float32)

# ── Load Q_H=2 vacuum background ──────────────────────────────────────
n_vacuum_np = None
n_vacuum    = None
TUBE_RADIUS = 1.98   # measured default; overwritten below if file loads

if args.k_restore_mode != 'off' or args.resume:
    if not os.path.exists(args.vacuum_profile):
        LOG(f"\nERROR: vacuum_profile not found: {args.vacuum_profile}")
        LOG("  Set --k_restore_mode off to run without the vacuum reference,")
        LOG("  or provide the correct path to f_fb_beta0_45200.npy.")
        sys.exit(1)

    LOG("\nLoading Q_H=2 vacuum background...")
    f_axisym = np.load(args.vacuum_profile)
    Nr_src, Nz_src = f_axisym.shape
    r_src = args.vacuum_h*(np.arange(Nr_src)+0.5)
    z_src = args.vacuum_h*np.arange(Nz_src)
    f_full = np.concatenate([f_axisym[:, ::-1], f_axisym[:, 1:]], axis=1)
    z_full = np.concatenate([-z_src[::-1], z_src[1:]])
    interp = RegularGridInterpolator(
        (r_src, z_full), f_full, bounds_error=False, fill_value=0.0)
    f3 = interp(np.stack([Rcyl.ravel(), Zc.ravel()],
                          axis=-1)).reshape(N,N,N)
    f3 = np.clip(f3, 0.0, np.pi)
    Phi_vac = ThetaAz - ThetaPol
    n_vacuum_np = np.stack([np.sin(f3)*np.cos(Phi_vac),
                              np.sin(f3)*np.sin(Phi_vac),
                              np.cos(f3)], axis=-1).astype(np.float32)
    n_vacuum_np /= np.linalg.norm(
        n_vacuum_np, axis=-1, keepdims=True).clip(1e-10)
    n_vacuum = torch.tensor(n_vacuum_np, dtype=torch.float32, device=dev)

    # Measure tube radius from profile (not hardcoded)
    f_z0 = f_axisym[:,0]
    idx  = np.argmin(np.abs(r_src - args.R0))
    for i in range(idx, len(r_src)):
        if f_z0[i] < np.pi/2:
            TUBE_RADIUS = r_src[i] - args.R0; break
    LOG(f"  Tube half-width: {TUBE_RADIUS:.4f}  |n_vac| OK")

# ── Energy and force ──────────────────────────────────────────────────
def compute_K_J4(n):
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a))/(2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 +
          nyx**2+nyy**2+nyz**2 +
          nzx**2+nzy**2+nzz**2)
    K  = (s4*g2).sum()*h**3 + MU*(g2.sum()*h**3)
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4  = (Fxy**2+Fxz**2+Fyz**2).sum()*h**3
    return K, J4


def spatial_force(n):
    """
    Tangent-projected -dE/dn for E = K_fb + phi^6*J4, plus the
    chosen restoring force.  k_restore_mode controls the restore:

      off                  → no restoring term
      flat                 → -k_restore*(n - n_vac)        [v1 behaviour]
      suppression_weighted → -k_restore * S(n) * (n - n_vac)
                             where S(n) = sin^4(θ) + epsilon*delta_S_ico(n)
                             sin^4(θ) = (1-nz²)²: continuous azimuthal average
                             delta_S_ico(n) = nz*(nx^5-10*nx^3*ny^2+5*nx*ny^4)
                               = Re[(nx+i*ny)^5*nz]: A_5-invariant icosahedral
                               correction, 5-fold azimuthal symmetry, zero mean.
                             epsilon=0 (default): recovers pure sin^4 weighting.
                             epsilon>0: adds discrete 5-fold icosahedral structure.
                             Only ico_epsilon>0 with k_restore_mode='suppression_weighted'
                             activates the icosahedral correction. In flat or off
                             modes, ico_epsilon is inert (correct physics: flat mode
                             ignores the suppression landscape entirely).
    """
    n_r = n.detach().requires_grad_(True)
    K, J4 = compute_K_J4(n_r)
    (K + phi6*J4).backward()
    g_s = -(n_r.grad.detach())
    g_s = g_s - (g_s*n).sum(-1, keepdim=True)*n   # project to tangent

    if args.k_restore_mode == 'flat':
        restore_raw = -args.k_restore*(n - n_vacuum)
        g_s = g_s + restore_raw - (restore_raw*n).sum(-1,keepdim=True)*n

    elif args.k_restore_mode == 'suppression_weighted':
        # Base suppression: sin^4(theta), azimuthally symmetric
        S = (1.0 - n[...,2:3]**2).clamp(0,1)**2   # (N,N,N,1)

        if args.ico_epsilon != 0.0:
            # Icosahedral correction: delta_S_ico = Re[(nx+i*ny)^5 * nz]
            # = nz * (nx^5 - 10*nx^3*ny^2 + 5*nx*ny^4)
            # This is the unique A_5-invariant degree-6 polynomial on S^2
            # with 5-fold azimuthal symmetry (verified: equal at all 10
            # non-polar icosahedral vertices, zero mean over azimuth).
            nx, ny, nz = n[...,0:1], n[...,1:2], n[...,2:3]
            delta_S = nz * (nx**5 - 10.0*nx**3*ny**2 + 5.0*nx*ny**4)
            S = S + args.ico_epsilon * delta_S
            # Clamp to non-negative (guaranteed for epsilon<2, but safe regardless)
            S = S.clamp(min=0.0)

        restore_raw = -args.k_restore * S * (n - n_vacuum)
        g_s = g_s + restore_raw - (restore_raw*n).sum(-1,keepdim=True)*n

    # 'off': nothing to add — ico_epsilon also has no effect here
    return g_s, K.item(), J4.item()


# ── Gamma field ───────────────────────────────────────────────────────
def build_gamma_field(envelope_np):
    """Returns a (N,N,N,1) float32 tensor, or None if all-zero (RATTLE)."""
    g = np.zeros((N,N,N), dtype=np.float32)
    if args.gamma > 0:
        if args.gamma_profile == 'uniform':
            g[:] = args.gamma
        elif args.gamma_profile == 'outside_kick':
            g[envelope_np <= 0.05] = args.gamma
        # boundary_only: only the sponge contributes below
    g += gamma_boundary_field
    if g.max() == 0:
        return None   # no damping → pure RATTLE
    return torch.tensor(g, dtype=torch.float32, device=dev).unsqueeze(-1)


# ── Kick ──────────────────────────────────────────────────────────────
def resolve_kick_axis(n_init_np, rng):
    """
    Returns a unit vector (3,) giving the direction to push n toward.

    kick_axis modes:
      'pole'    — toward nearest ±z pole (original behaviour, v1-compatible)
      'x','y','z' — fixed Cartesian axis
      'a,b,c'   — explicit direction, e.g. '0,0.707,0.707'; normalised here
      'random'  — random unit vector drawn from rng (seeded by --kick_seed)

    In all modes except 'pole', the SAME axis is used at every site in the
    kick region. 'pole' uses a site-by-site sign flip to always push toward
    the closer of ±kick_axis, giving a smooth kick even where n_z changes sign.
    """
    ka = args.kick_axis.strip().lower()
    if ka == 'x':
        return np.array([1.0, 0.0, 0.0])
    elif ka == 'y':
        return np.array([0.0, 1.0, 0.0])
    elif ka == 'z' or ka == 'pole':
        return np.array([0.0, 0.0, 1.0])   # sign flip applied per-site below
    elif ka == 'random':
        v = rng.standard_normal(3)
        return v / np.linalg.norm(v)
    else:
        # parse 'a,b,c'
        try:
            parts = [float(x) for x in args.kick_axis.split(',')]
            assert len(parts) == 3
            v = np.array(parts)
            return v / np.linalg.norm(v)
        except Exception:
            raise ValueError(f"Cannot parse --kick_axis '{args.kick_axis}'. "
                             "Use 'pole', 'x', 'y', 'z', 'random', or 'a,b,c'.")


def build_kick(n_init_np):
    # Always use a concrete seed so the run is reproducible from the snapshot.
    # If the user didn't supply one, generate from OS entropy and record it.
    if args.kick_seed is None:
        args.kick_seed = int(np.random.SeedSequence().entropy & 0xFFFFFFFF)
    rng = np.random.default_rng(args.kick_seed)

    # Kick centre: azimuthal angle around the torus, plus poloidal offset
    az = rng.uniform(0.0, 2*np.pi) if args.kick_azimuth == 'random' else float(args.kick_azimuth)
    r0_kick = args.R0 + TUBE_RADIUS
    cx = r0_kick * np.cos(az)
    cy = r0_kick * np.sin(az)
    cz = TUBE_RADIUS * np.sin(args.kick_r0_offset)
    centre = np.array([cx, cy, cz])
    dist   = np.sqrt((X-cx)**2 + (Y-cy)**2 + (Zc-cz)**2)
    env    = np.exp(-(dist**2) / (2*args.kick_sigma**2))

    # Kick direction
    axis = resolve_kick_axis(n_init_np, rng)   # (3,) unit vector

    if args.kick_axis.strip().lower() == 'pole':
        # Site-by-site: push toward the ±axis that n is already closer to.
        # This is the original 'pole' behaviour: axis = z, sign flip per site.
        dot_with_n = (n_init_np * axis).sum(-1, keepdims=True)   # (N,N,N,1)
        sign = np.sign(dot_with_n); sign[sign == 0] = 1.0
        target = axis * sign        # (N,N,N,3), sign-flipped per site
    else:
        # All sites pushed toward the same axis direction
        target = np.broadcast_to(axis, n_init_np.shape).copy()

    raw       = target - n_init_np
    raw_t     = raw - (raw * n_init_np).sum(-1, keepdims=True) * n_init_np
    direction = raw_t / np.linalg.norm(raw_t, axis=-1, keepdims=True).clip(1e-10)
    v = (args.kick_strength * env[..., None] * direction).astype(np.float32)
    return v, env, axis, centre, az


# ── Integrator steps ──────────────────────────────────────────────────
def step_first_half(n, v, a, gamma_field):
    """Half-velocity update + geodesic position step.
    With gamma_field=None this is identical to v1's RATTLE half-step."""
    exp_h = 1.0 if gamma_field is None else torch.exp(-gamma_field*(args.dt/2))
    v_half = (v + 0.5*args.dt*a) * exp_h
    v_half = v_half - (v_half*n).sum(-1,keepdim=True)*n
    speed  = v_half.norm(dim=-1, keepdim=True).clamp(1e-12)
    th     = speed*args.dt
    d      = v_half/speed
    n_new  = n*torch.cos(th) + d*torch.sin(th)
    n_new  = n_new / n_new.norm(dim=-1, keepdim=True).clamp(1e-10)
    return n_new, v_half

def step_second_half(n_new, v_half, a_new, gamma_field):
    """Second half-velocity update + tangent re-projection."""
    exp_h  = 1.0 if gamma_field is None else torch.exp(-gamma_field*(args.dt/2))
    v_un   = (v_half + 0.5*args.dt*a_new) * exp_h
    return v_un - (v_un*n_new).sum(-1,keepdim=True)*n_new


# ── Initialise ────────────────────────────────────────────────────────
start_step = 0
if args.resume:
    LOG(f"\nResuming from: {args.resume}")
    snap       = np.load(args.resume, allow_pickle=True)
    n_np       = snap['n']
    v_np       = snap.get('v', np.zeros_like(snap['n']))
    start_step = int(snap.get('step', 0))
    if int(snap.get('N', N)) != N:
        LOG("  WARNING: snapshot N differs — field loaded as-is.")
    n = torch.tensor(n_np, dtype=torch.float32, device=dev)
    v = torch.tensor(v_np, dtype=torch.float32, device=dev)
    v = v - (v*n).sum(-1,keepdim=True)*n
    KE_injected = 0.5*args.mu2*(v**2).sum().item()*h**3
    LOG(f"  Resumed step={start_step}  KE={KE_injected:.6f}")
    envelope = np.zeros((N,N,N), dtype=np.float32)
    kick_az  = float(snap.get('kick_azimuth', 0.0))   # from original run
else:
    # Fresh start — use n_vacuum if loaded, else flat (north-pole) field
    if n_vacuum_np is not None:
        n_init_np = n_vacuum_np
    else:
        n_init_np = np.zeros((N,N,N,3), dtype=np.float32)
        n_init_np[...,2] = 1.0   # uniform north-pole field
    n = torch.tensor(n_init_np.copy(), dtype=torch.float32, device=dev)
    v_np, envelope, kick_axis_vec, kick_centre, kick_az = build_kick(n_init_np)
    v = torch.tensor(v_np, dtype=torch.float32, device=dev)
    if n_vacuum is not None:
        v = v - (v*n_vacuum).sum(-1,keepdim=True)*n_vacuum
    KE_injected = 0.5*args.mu2*(v**2).sum().item()*h**3
    LOG(f"\n  Kick centre: [{kick_centre[0]:.3f} {kick_centre[1]:.3f} {kick_centre[2]:.3f}]"
        f"  azimuth={kick_az:.4f} rad"
        + (" (was 'random')" if args.kick_azimuth == 'random' else ""))
    LOG(f"  Kick axis: {args.kick_axis}"
        + (f" → [{kick_axis_vec[0]:.4f} {kick_axis_vec[1]:.4f} {kick_axis_vec[2]:.4f}]"
           if args.kick_axis not in ('x','y','z') else ""))
    LOG(f"  KE injected: {KE_injected:.6f}  "
        f"sites: {int((envelope>0.05).sum())} (envelope>0.05)")

gamma_field = build_gamma_field(envelope)
kick_mask   = torch.tensor(envelope > 0.05, device=dev)

a, K0, J4_0 = spatial_force(n)
a = a/args.mu2
E0 = K0 + phi6*J4_0 + KE_injected
J4_vacuum = J4_0

LOG(f"\n  Start: K={K0:.4f}  J4={J4_0:.4f}  E_total={E0:.4f}")
LOG(f"  Integrator: {'RATTLE (conservative)' if gamma_field is None else 'Langevin (dissipative)'}")
LOG(f"  Sponge: {'off' if args.absorb_width==0 else f'width={args.absorb_width}'}")

# ── Trajectory CSV ────────────────────────────────────────────────────
csv_path    = os.path.join(args.outdir, 'trajectory.csv')
write_hdr   = not os.path.exists(csv_path)
csv_f       = open(csv_path, 'a')
if write_hdr:
    csv_f.write("step,t,J4,K,E_total,KE,theta_min_kick,"
                 "snapback,KJ4,J4_frac_vacuum\n")

LOG(f"\n{'step':>8} {'t':>9} {'J4':>9} {'J4/J4_0':>8} {'K':>10} "
    f"{'KE':>9} {'θ_min':>7} {'snapback':>10}")
LOG("-"*80)

best_J4 = J4_0; best_J4_step = start_step

# ── Main loop ─────────────────────────────────────────────────────────
for step in range(start_step+1, start_step+args.steps+1):
    t = step*args.dt

    n_new, v_half = step_first_half(n, v, a, gamma_field)
    a_new, K_new, J4_new = spatial_force(n_new)
    a_new = a_new/args.mu2
    v_new = step_second_half(n_new, v_half, a_new, gamma_field)
    n, v, a = n_new, v_new, a_new

    if J4_new < best_J4:
        best_J4 = J4_new; best_J4_step = step

    if step % args.print_every == 0 or step == start_step+1:
        with torch.no_grad():
            KE     = 0.5*args.mu2*(v**2).sum().item()*h**3
            E_tot  = K_new + phi6*J4_new + KE
            if kick_mask.any():
                th_k  = torch.acos(n[...,2][kick_mask].clamp(-1,1))
                th_min = torch.minimum(th_k, np.pi-th_k).min().item()
            else:
                th_min = float('nan')
            sb     = ((n - n_vacuum)**2).sum().item()*h**3 \
                     if n_vacuum is not None else float('nan')
            KJ4    = K_new/J4_new if J4_new > 1e-12 else float('nan')
            J4_frac = J4_new/J4_vacuum if J4_vacuum > 1e-12 else float('nan')
        LOG(f"{step:>8} {t:>9.3f} {J4_new:>9.3f} {J4_frac:>8.4f} "
            f"{K_new:>10.3f} {KE:>9.5f} {th_min:>7.4f} {sb:>10.3f}")
        csv_f.write(f"{step},{t:.6f},{J4_new:.6f},{K_new:.6f},{E_tot:.6f},"
                    f"{KE:.6f},{th_min:.6f},{sb:.6f},{KJ4:.6f},{J4_frac:.6f}\n")
        csv_f.flush()

    if step % args.snapshot_every == 0 or step == start_step+1:
        sv = dict(n=n.detach().cpu().numpy(),
                  v=v.detach().cpu().numpy(),
                  step=step, t=t, J4=J4_new, K=K_new,
                  N=N, h=h, R0=args.R0,
                  kick_strength=args.kick_strength,
                  kick_axis=args.kick_axis,
                  kick_azimuth=kick_az,       # resolved float, not 'random'
                  kick_seed=args.kick_seed,
                  k_restore=args.k_restore,
                  k_restore_mode=args.k_restore_mode,
                  gamma=args.gamma,
                  gamma_profile=args.gamma_profile,
                  absorb_width=args.absorb_width,
                  ico_epsilon=args.ico_epsilon)
        if n_vacuum_np is not None:
            sv['n_vacuum'] = n_vacuum_np
        np.savez(os.path.join(args.outdir, 'snapshots',
                               f'n_t{step:08d}.npz'), **sv)

csv_f.close()

# ── Summary ───────────────────────────────────────────────────────────
LOG(f"\n{'='*70}")
LOG(f"  SUMMARY  steps {start_step}→{start_step+args.steps}")
LOG(f"  J4_start = {J4_vacuum:.4f}")
LOG(f"  J4_min   = {best_J4:.4f}  at step {best_J4_step} "
    f"(t={best_J4_step*args.dt:.3f})")
LOG(f"  J4_min/J4_start = {best_J4/J4_vacuum:.4f}")
LOG(f"  J4_final = {J4_new:.4f}  ({J4_new/J4_vacuum*100:.1f}% of start)")
if n_vacuum is not None:
    LOG(f"  snapback_final = {sb:.3f}  (0=perfect return)")
steps_to = max(1, best_J4_step - start_step)
LOG(f"\n  Geometry check (index[2I:2T]=5):")
LOG(f"    steps_to_min={steps_to}  "
    f"total/steps_to_min={args.steps/steps_to:.3f}  (vs 5)")
log.close()
