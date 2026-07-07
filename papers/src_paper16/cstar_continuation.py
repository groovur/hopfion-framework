#!/usr/bin/env python3
"""
cstar_continuation.py
======================
Prototype C*-continuation solver for the Q_H=3 trefoil sector.

STRATEGY (per user direction): rather than annealing once at the
physical C*=2.5062 (where the existing grid badly under-resolves the
tube, Proposition prop:j4_resolution_artifact), perform a sequence of
SHORT anneals at a path of increasing C* values, starting from a
small, well-resolved C* and stepping up toward 2.5062, using EACH
step's relaxed field as the WARM START (not a fresh random/BS-ansatz
restart) for the NEXT C* value. This is a predictor-corrector-style
numerical continuation: track how the energy-minimising configuration
("saddle" in Paper XV's sense -- a stationary point of the full
E_geom=K*J4 functional) MOVES as C* is incremented, rather than
re-discovering it from scratch at each value.

This also tests a genuine PHYSICS question: does a smoothly-varying,
Q_H=3-preserving local minimum exist along the WHOLE path from small
C* to the physical value, or does something break down partway
(bifurcation, sudden energy jump, charge-changing event)? This is the
solver's REAL dynamics version of the earlier (purely kinematic,
diagnostic-only) C* sweep of Section sec:cstar_threshold.

GRID DESIGN NOTE (read this before changing N, h): per the diagnosed
resolution requirement, a FIXED grid must remain adequate for the
FINAL, smallest tube radius (1/C*_target). This prototype uses a
SMALL, CHEAP grid by default so the machinery can be validated
end-to-end quickly; production runs should significantly increase N
and reduce h (see --N, --h) once the pipeline is confirmed correct.

The per-strand Phi=chi+3t construction and its compensated Bishop
frame (Construction constr:compensated_frame, constr:perstrand) are
ported directly from perstrand_s3lift_v2.py / bishop_frame_v2.py.
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, default=48, help='grid points per axis')
ap.add_argument('--h', type=float, default=0.18, help='grid spacing')
ap.add_argument('--R0', type=float, default=3.0)
ap.add_argument('--r0', type=float, default=0.874)
ap.add_argument('--C_star_start', type=float, default=0.5)
ap.add_argument('--C_star_target', type=float, default=2.5062)
ap.add_argument('--n_steps', type=int, default=8, help='number of C* increments')
ap.add_argument('--sweeps_per_step', type=int, default=60,
                 help='Metropolis sweeps per C* increment (short, for prototype)')
ap.add_argument('--T0_frac', type=float, default=0.003,
                 help='initial temperature as fraction of starting E_geom, per step')
ap.add_argument('--cooling_rate', type=float, default=0.95)
ap.add_argument('--rotation_scale', type=float, default=0.2)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--outdir', type=str, default='.')
ap.add_argument('--device', type=str, default='cpu')
ap.add_argument('--blend_alpha', type=float, default=0.5,
                 help='predictor blend weight: 1.0 = use fresh reference field at each step '
                      '(fully replaces relaxed structure), 0.0 = ignore the new C* entirely '
                      '(keeps previous relaxed field unchanged). Slerp between previous '
                      'relaxed field and fresh per-strand reference at new C*.')
args = ap.parse_args()

torch.manual_seed(args.seed)
np.random.seed(args.seed)

PHI = (1+5**0.5)/2
MU = 3.0 - PHI
N, h = args.N, args.h
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
log_path = os.path.join(args.outdir, 'cstar_continuation_log.json')

print(f"{'='*70}")
print(f"  C*-CONTINUATION PROTOTYPE")
print(f"  Grid: N={N}, h={h}, box=[{-N*h/2:.2f},{N*h/2:.2f}]")
print(f"  C* path: {args.C_star_start} -> {args.C_star_target} in {args.n_steps} steps")
print(f"  Sweeps per step: {args.sweeps_per_step}")
print(f"{'='*70}")

# ── Grid ──────────────────────────────────────────────────────────────
cv = h*(np.arange(N) - N//2 + 0.5)
pts_np = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), axis=-1).reshape(-1, 3).astype(np.float32)
dist_from_origin = torch.tensor(
    np.linalg.norm(pts_np, axis=-1).reshape(N, N, N), dtype=torch.float32, device=dev)
pts_t = torch.tensor(pts_np, dtype=torch.float32, device=dev)

box = N*h
print(f"  Tube radius at target C*={args.C_star_target}: {1/args.C_star_target:.3f}  "
      f"(grid pts across diam: {2/args.C_star_target/h:.1f})")
print(f"  Tube radius at start  C*={args.C_star_start}: {1/args.C_star_start:.3f}  "
      f"(grid pts across diam: {2/args.C_star_start/h:.1f})")

needed_half_width = args.R0 + args.r0 + 1.0/min(args.C_star_start, args.C_star_target) + 0.5
actual_half_width = box/2
print(f"  Box half-width: {actual_half_width:.3f}  (needed for full C* path: {needed_half_width:.3f})")
if actual_half_width < needed_half_width:
    print(f"  FATAL: box too small for the fattest tube on this C* path "
          f"(needs half-width >= {needed_half_width:.3f}). Increase N or h. Aborting.")
    sys.exit(1)

# ── Per-strand construction machinery (numpy, built once per C*) ─────
R0, r0 = args.R0, args.r0
NT_frame = 20000
t_frame, T_frame, N1_frame, N2_frame, H_holonomy = build_compensated_frame_arclength(NT=NT_frame)
print(f"  Compensated Bishop frame holonomy: {np.degrees(H_holonomy):.4f} deg (closed)")

NT = 4000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
lobe_indices = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees = [KDTree(Gamma_pts[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr[li] for li in lobe_indices]

def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

def nearest_two_strands_robust(qpts):
    d_per_lobe, t_per_lobe = [], []
    for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
        d, idx = tree_l.query(qpts)
        d_per_lobe.append(d); t_per_lobe.append(t_l[idx])
    d_stack = np.stack(d_per_lobe, axis=1)
    t_stack = np.stack(t_per_lobe, axis=1)
    order = np.argsort(d_stack, axis=1)
    d_sorted = np.take_along_axis(d_stack, order, axis=1)
    t_sorted = np.take_along_axis(t_stack, order, axis=1)
    return t_sorted[:,0], d_sorted[:,0], t_sorted[:,1], d_sorted[:,1]

def chi_for_t(qpts, t_query):
    N1_pts, N2_pts = frame_at_t(t_query)
    curve_pts = curve_at_t(t_query)
    rel = qpts - curve_pts
    chi = np.arctan2(np.sum(rel*N2_pts,axis=1), np.sum(rel*N1_pts,axis=1))
    dist = np.linalg.norm(rel, axis=1)
    return chi, dist

# precompute (t1,d1,t2,d2) ONCE -- these don't depend on C*, only on
# the fixed grid and curve geometry, so building the field at a new
# C* is cheap (just re-evaluates f0 and the doublet, not the strand search)
print("  Precomputing per-point strand assignment (independent of C*)...")
t0 = time.time()
t1_g, d1_g, t2_g, d2_g = nearest_two_strands_robust(pts_np)
chi1_g, _ = chi_for_t(pts_np, t1_g)
chi2_g, _ = chi_for_t(pts_np, t2_g)
Phi1_g = chi1_g + 3*t1_g
Phi2_g = chi2_g + 3*t2_g
print(f"  done ({time.time()-t0:.1f}s)")

d1_t = torch.tensor(np.clip(d1_g,1e-6,None), dtype=torch.float32, device=dev)
d2_t = torch.tensor(np.clip(d2_g,1e-6,None), dtype=torch.float32, device=dev)
Phi1_t = torch.tensor(Phi1_g, dtype=torch.float32, device=dev)
Phi2_t = torch.tensor(Phi2_g, dtype=torch.float32, device=dev)

def build_n_C_star(C_star):
    """Build the per-strand n field at a given C*, using the PRECOMPUTED
    (C*-independent) strand assignment and poloidal phases."""
    f1 = 2*torch.atan(torch.pow(d1_t*C_star, -C_star))
    f2 = 2*torch.atan(torch.pow(d2_t*C_star, -C_star))
    w1 = 1.0/(d1_t**2)
    w2 = 1.0/(d2_t**2)
    z1a = torch.cos(f1/2).to(torch.complex64)
    z2a = (torch.sin(f1/2)*torch.exp(1j*Phi1_t))
    z1b = torch.cos(f2/2).to(torch.complex64)
    z2b = (torch.sin(f2/2)*torch.exp(1j*Phi2_t))
    z1 = w1*z1a + w2*z1b
    z2 = w1*z2a + w2*z2b
    norm = torch.sqrt(torch.abs(z1)**2 + torch.abs(z2)**2)
    z1, z2 = z1/norm, z2/norm
    nx = 2*torch.real(torch.conj(z1)*z2)
    ny = 2*torch.imag(torch.conj(z1)*z2)
    nz = torch.abs(z1)**2 - torch.abs(z2)**2
    n = torch.stack([nx,ny,nz], dim=-1).reshape(N,N,N,3)
    n = n / n.norm(dim=-1, keepdim=True).clamp(1e-10)
    return n

print("\nQUICK SELF-CHECK: build_n_C_star matches standalone perstrand_s3lift_v2 at C*=2.5062?")
nn = build_n_C_star(2.5062)
print(f"  shape={nn.shape}, |n| range=[{nn.norm(dim=-1).min().item():.5f},{nn.norm(dim=-1).max().item():.5f}]")

# ── Solver's actual energy functional, ported verbatim ───────────────
def compute_global(n):
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
    r_bar_num = (rho_J4 * dist_from_origin).sum()
    r_bar_den = rho_J4.sum()
    r_bar = (r_bar_num/r_bar_den).item() if r_bar_den.item() > 1e-12 else float('nan')
    return K.item(), J2a.item(), J4.item(), r_bar


def local_KJ4_contribution(n, mask):
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a_density  = s4 * g2
    J2iso_density = g2
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    return J2a_density, J2iso_density, rho_J4


def random_rotation_perturb(n, mask, scale):
    idx = mask.nonzero(as_tuple=True)
    if idx[0].shape[0] == 0:
        return n[mask], idx
    v = n[mask]
    rand_axis = torch.randn(v.shape[0], 3, device=dev)
    rand_axis = rand_axis - (rand_axis*v).sum(-1,keepdim=True)*v
    axis_norm = rand_axis.norm(dim=-1, keepdim=True).clamp(1e-10)
    rand_axis = rand_axis/axis_norm
    angle = (torch.rand(v.shape[0],1, device=dev)*2-1) * scale
    v_new = v*torch.cos(angle) + torch.cross(rand_axis, v, dim=-1)*torch.sin(angle)
    v_new = v_new / v_new.norm(dim=-1, keepdim=True).clamp(1e-10)
    return v_new, idx


# Precompute the 64 colour masks (mod-4 on each axis), exactly as in
# qh3_trefoil_solver_3d_v11.py (verified-exact stencil decomposition).
ar = torch.arange(N, device=dev)
I, J_, K_ = torch.meshgrid(ar, ar, ar, indexing='ij')
COLOUR_MASKS = []
for ci in range(4):
    for cj in range(4):
        for ck in range(4):
            COLOUR_MASKS.append((I%4==ci) & (J_%4==cj) & (K_%4==ck))


def anneal_sweep(n, K_cur, J4_cur, T, rotation_scale):
    n_accept = 0
    n_total = 0
    for mask in COLOUR_MASKS:
        v_new, idx = random_rotation_perturb(n, mask, rotation_scale)
        if idx[0].shape[0] == 0:
            continue
        J2a_d0, J2iso_d0, rho0 = local_KJ4_contribution(n, mask)
        n_trial = n.clone()
        n_trial[mask] = v_new
        J2a_d1, J2iso_d1, rho1 = local_KJ4_contribution(n_trial, mask)
        dJ2a = (J2a_d1 - J2a_d0).sum() * h**3
        dJ2iso = (J2iso_d1 - J2iso_d0).sum() * h**3
        dK = (dJ2a + MU*dJ2iso).item()
        dJ4 = ((rho1 - rho0).sum() * h**3).item()
        K_new = K_cur + dK
        J4_new = J4_cur + dJ4
        dE = (K_new*J4_new) - (K_cur*J4_cur)
        accept = dE < 0 or np.random.rand() < np.exp(-dE/max(T, 1e-12))
        n_total += idx[0].shape[0]
        if accept:
            n = n_trial
            K_cur, J4_cur = K_new, J4_new
            n_accept += idx[0].shape[0]
    accept_rate = n_accept / max(n_total, 1)
    return n, K_cur, J4_cur, accept_rate

print("\nEnergy functional and annealing machinery loaded.")
print(f"Colour classes: {len(COLOUR_MASKS)}")


def probe_dE(n, K_cur, J4_cur, mask, rotation_scale):
    """Dry-run a single colour-class proposal (no mutation) to measure
    the typical |dE_geom| scale for THIS field -- ported verbatim from
    the real solver's calibration mechanism (see calibrate_T0 docstring
    for why this is necessary: guessing T0 as a fraction of the TOTAL
    E_geom=K*J4 is wrong by orders of magnitude, since a single colour-
    class move only ever touches 1/64 of the grid)."""
    v_new, idx = random_rotation_perturb(n, mask, rotation_scale)
    if idx[0].shape[0] == 0:
        return 0.0
    J2a_d0, J2iso_d0, rho0 = local_KJ4_contribution(n, mask)
    n_trial = n.clone()
    n_trial[mask] = v_new
    J2a_d1, J2iso_d1, rho1 = local_KJ4_contribution(n_trial, mask)
    dJ2a = (J2a_d1 - J2a_d0).sum() * h**3
    dJ2iso = (J2iso_d1 - J2iso_d0).sum() * h**3
    dK = (dJ2a + MU*dJ2iso).item()
    dJ4 = ((rho1 - rho0).sum() * h**3).item()
    K_new = K_cur + dK
    J4_new = J4_cur + dJ4
    dE = (K_new*J4_new) - (K_cur*J4_cur)
    return abs(dE)


def calibrate_T0(n, K_cur, J4_cur, rotation_scale, n_probes=20):
    """Ported verbatim from qh3_trefoil_solver_3d_v11.py: measure the
    ACTUAL typical scale of a single colour-class's |dE_geom| on THIS
    field, rather than guessing T0 as a fraction of total E_geom.

    BUG FIX (v2): original used T0 = median_dE / ln(2), which sets 50%
    initial acceptance on a median-sized move — far too hot, causes the
    chain to accept almost everything early and scramble the field before
    cooling can help. Changed to T0 = median_dE / 4, giving ~19% initial
    acceptance on a median move, closer to the standard SA target of
    10-20% at the start of cooling."""
    probe_colours = [COLOUR_MASKS[i] for i in
                      np.random.choice(len(COLOUR_MASKS), size=min(n_probes, len(COLOUR_MASKS)), replace=False)]
    dEs = [probe_dE(n, K_cur, J4_cur, mask, rotation_scale) for mask in probe_colours]
    dEs = [d for d in dEs if d > 0]
    if not dEs:
        return 1.0
    median_dE = float(np.median(dEs))
    T0 = median_dE / 4.0   # ~19% acceptance on median move at start
    print(f"    [calibration] probed {len(dEs)} colour-class moves: "
          f"median|dE|={median_dE:.4e}  max|dE|={max(dEs):.4e}  -> T0={T0:.4e}")
    return T0

# ── C*-CONTINUATION MAIN LOOP ──────────────────────────────────────────
def slerp_fields(n_a, n_b, alpha):
    """Spherical linear interpolation between two unit-vector fields,
    pointwise. alpha=0 -> n_a, alpha=1 -> n_b. Falls back to linear
    blend (renormalised) where n_a, n_b are nearly parallel or
    antiparallel, to avoid the slerp formula's 0/0 singularity there."""
    dot = (n_a*n_b).sum(-1, keepdim=True).clamp(-1+1e-6, 1-1e-6)
    theta = torch.acos(dot)
    sin_theta = torch.sin(theta).clamp(1e-6)
    near_degenerate = (sin_theta < 1e-4)
    w_a = torch.sin((1-alpha)*theta)/sin_theta
    w_b = torch.sin(alpha*theta)/sin_theta
    out = w_a*n_a + w_b*n_b
    # fallback: plain linear blend + renormalise, where slerp is ill-conditioned
    lin = (1-alpha)*n_a + alpha*n_b
    lin = lin / lin.norm(dim=-1, keepdim=True).clamp(1e-10)
    out = torch.where(near_degenerate, lin, out)
    return out / out.norm(dim=-1, keepdim=True).clamp(1e-10)


# ── C*-CONTINUATION MAIN LOOP (predictor-corrector) ───────────────────
C_star_path = np.linspace(args.C_star_start, args.C_star_target, args.n_steps)
print(f"\nC* path: {[f'{c:.4f}' for c in C_star_path]}")
print(f"Blend alpha (predictor weight toward fresh reference field): {args.blend_alpha}")

history = []
n_field = None
t_run0 = time.time()

for step_i, C_star in enumerate(C_star_path):
    t_step0 = time.time()
    n_ref = build_n_C_star(C_star)
    if n_field is None:
        n_field = n_ref
        warm = False
        blend_used = None
    else:
        warm = True
        blend_used = args.blend_alpha
        n_field = slerp_fields(n_field, n_ref, args.blend_alpha)

    K_cur, J2a_cur, J4_cur, r_bar_cur = compute_global(n_field)
    E_cur = K_cur * J4_cur
    print(f"\n--- Step {step_i+1}/{args.n_steps}  C*={C_star:.4f}  "
          f"(warm_start={warm}, blend_alpha={blend_used}) ---")
    print(f"  post-predictor:  K={K_cur:.3f} J4={J4_cur:.3f} E={E_cur:.3f} r_bar={r_bar_cur:.3f}")

    # BUG FIX: calibrate T0 on n_ref (the pure analytic construction at this
    # C*), NOT on n_field (the slerp blend, which at steps 2+ contains structure
    # from the previous step's potentially-scrambled post-corrector field).
    # Calibrating on n_ref measures the actual energy landscape we're about to
    # explore; calibrating on n_field measured a corrupted proxy that gave T0
    # values ~50x too large, causing the corrector to accept almost all moves
    # and scramble the field rather than relax it.
    K_ref, _, J4_ref, _ = compute_global(n_ref)
    T0 = calibrate_T0(n_ref, K_ref, J4_ref, args.rotation_scale)
    T = T0
    for sweep in range(args.sweeps_per_step):
        n_field, K_cur, J4_cur, acc = anneal_sweep(n_field, K_cur, J4_cur, T, args.rotation_scale)
        T *= args.cooling_rate

    K_f, J2a_f, J4_f, r_bar_f = compute_global(n_field)
    E_f = K_f*J4_f
    dt_step = time.time()-t_step0
    print(f"  post-corrector:  K={K_f:.3f} J4={J4_f:.3f} E={E_f:.3f} r_bar={r_bar_f:.3f}  "
          f"(accept_rate(last)={acc:.3f}, {dt_step:.1f}s)")

    history.append(dict(step=step_i, C_star=float(C_star), warm_start=warm, blend_alpha=blend_used,
                         K_pred=K_cur, J4_pred=J4_cur, E_pred=E_cur, r_bar_pred=r_bar_cur,
                         K_post=K_f, J4_post=J4_f, E_post=E_f, r_bar_post=r_bar_f,
                         T0=T0, accept_rate_last=acc, wall_time=dt_step))

print(f"\nTotal continuation wall time: {time.time()-t_run0:.1f}s")

with open(log_path, 'w') as f:
    json.dump(history, f, indent=2)
print(f"History written to {log_path}")

print(f"\n{'C*':>8}  {'K_post':>10}  {'J4_post':>10}  {'E_post':>14}  {'r_bar_post':>11}")
for hrec in history:
    print(f"{hrec['C_star']:>8.4f}  {hrec['K_post']:>10.3f}  {hrec['J4_post']:>10.3f}  "
          f"{hrec['E_post']:>14.3f}  {hrec['r_bar_post']:>11.3f}")
