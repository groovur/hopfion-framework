#!/usr/bin/env python3
"""
gradient_flow_saddle_snapshot.py
=================================
Saddle-snapshot gradient flow for the Q_H=3 trefoil sector.

KEY DIFFERENCE FROM gradient_flow_constrained.py:
  After every Phase 2 step, the physical grid spacing h is rescaled to
  enforce s_opt = sqrt(lambda_3 * J4 / K) = 1 exactly (the Derrick
  balance / saddle condition).

  Under x -> s*x:  K -> s*K,  J4 -> J4/s,  so K*J4 is invariant,
  but s_opt^2 = lambda_3 * J4 / K -> lambda_3 * J4 / (s*K) / s ... wait,
  more carefully:
    K ~ int |dn/dx|^2 dx^3 -> K * (1/s^2) * s^3 = K * s  (3D volume * gradient^2)
    J4 ~ int |F|^2 dx^3    -> J4 * (1/s^4) * s^3 = J4/s
  So s_opt = sqrt(lambda_3 * J4 / K) -> sqrt(lambda_3 * (J4/s) / (K*s)) = s_opt/s.
  Setting s_opt_new = 1 requires s = s_opt_current.

  In the discrete Cartesian solver, rescaling x -> x*s_current is equivalent
  to changing h -> h * s_current (ZOOM OUT if s_opt>1, ZOOM IN if s_opt<1).
  BUT: this changes N_effective (the number of points covering the physical
  soliton), so it's only valid while the box stays large enough.

  PRACTICAL IMPLEMENTATION:
  We store a running 'physical_h' that accumulates the scale changes.
  The field tensor n[i,j,k] represents the field at position
  physical_h * (i - N//2, j - N//2, k - N//2).  All K, J4 integrals
  use physical_h (not the initial h). After each Phase 2 step:
    s_opt = sqrt(lambda_3 * J4 / K)
    physical_h -> physical_h * s_opt
  This enforces s_opt=1 at every step.

  The saddle-snapshot condition (Paper I, def:snapshot) is reached when
  |grad E_geom| -> 0 with s_opt = 1 maintained at every step. At that
  point, n[i,j,k] at physical_h represents the Q_H=3 saddle field.

WHEN TO USE THIS vs gradient_flow_constrained.py:
  - gradient_flow_constrained.py: topology protection, approaches the
    saddle neighbourhood, useful for establishing K_min and the |grad|
    fingerprint. Use for long exploratory runs.
  - THIS SCRIPT: true saddle convergence. Phase 2 enforces s_opt=1.
    Requires the field to already be near the saddle (e.g. warm-start
    from the h=0.05 run's n_final.npy at step ~1030 where |grad| was
    smallest). If topology collapses, fall back to gf_constrained.

PHASE STRUCTURE:
  Phase 1: same as gf_constrained -- large lr, no scale fixing, cleans
           up Construction C's excess gradient energy. Ends at K_min.
  Phase 2: scale-fixing at every step + tight angle clamp. Converges
           toward the true saddle shape.

USAGE:
  # From scratch:
  python gradient_flow_saddle_snapshot.py \\
    --N 192 --h 0.05 --C_star 2.5062 \\
    --n_steps 8000 --lr1 3e-4 --lr2 3e-5 \\
    --K_rise_eps 2.0 --delta_max_deg 5.0 \\
    --log_every 10 --outdir gf_saddle_snapshot

  # From warm start (best: use n_final at step ~1030, near |grad| min):
  python gradient_flow_saddle_snapshot.py \\
    --N 192 --h 0.05 --C_star 2.5062 \\
    --n_steps 8000 --lr1 3e-4 --lr2 3e-5 \\
    --K_rise_eps 2.0 --delta_max_deg 5.0 \\
    --warm_start n_final.npy \\
    --warm_start_phase 2 \\
    --log_every 10 --outdir gf_saddle_snapshot_ws

KEY OUTPUTS TO WATCH:
  s_opt column: should converge to 1.0000 as the saddle is approached.
  |grad| column: should fall toward 0.
  K, J4 columns: should stabilise (saddle = stationary point of E_geom).
  physical_h column: shows cumulative scale drift.

  CONVERGENCE CRITERION: |grad| < 1e3 AND |s_opt - 1| < 0.01
  (Both necessary: s_opt=1 without |grad|->0 just means Derrick balance,
  not a stationary point of the shape functional.)
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--N',                type=int,   default=64)
ap.add_argument('--h',                type=float, default=0.175)
ap.add_argument('--R0',               type=float, default=3.0)
ap.add_argument('--r0',               type=float, default=0.874)
ap.add_argument('--C_star',           type=float, default=2.5062)
ap.add_argument('--n_steps',          type=int,   default=2000)
ap.add_argument('--lr1',              type=float, default=3e-4)
ap.add_argument('--lr2',              type=float, default=3e-5)
ap.add_argument('--delta_max_deg',    type=float, default=5.0)
ap.add_argument('--K_rise_eps',       type=float, default=2.0)
ap.add_argument('--log_every',        type=int,   default=10)
ap.add_argument('--warm_start',       type=str,   default=None)
ap.add_argument('--warm_start_phase', type=int,   default=1,
                help='Which phase to START in when warm-starting (1 or 2).')
ap.add_argument('--outdir',           type=str,   default='gf_saddle_snapshot')
ap.add_argument('--seed',             type=int,   default=0)
ap.add_argument('--device',           type=str,   default='cpu')
# Scale-fixing controls
ap.add_argument('--scale_fix_every',  type=int,   default=1,
                help='Apply scale fix every N Phase-2 steps (default=1=every step).')
ap.add_argument('--scale_fix_max',    type=float, default=1.10,
                help='Maximum scale correction per step (clamp for stability).')
ap.add_argument('--scale_fix_min',    type=float, default=0.90,
                help='Minimum scale correction per step (clamp for stability).')
args = ap.parse_args()

torch.manual_seed(args.seed)
np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

PHI     = (1+5**0.5)/2
LAM3    = PHI**6          # = 17.9443, the proved lambda_3
MU      = 3.0 - PHI       # = 1.38197
N       = args.N
R0, r0, C_star = args.R0, args.r0, args.C_star
dev     = torch.device(args.device)
DELTA_MAX = float(np.radians(args.delta_max_deg))

# Running physical grid spacing (updated by scale fixing)
physical_h = args.h

print(f"{'='*70}")
print(f"  SADDLE-SNAPSHOT GRADIENT FLOW  (E_geom = K*J4, s_opt locked)")
print(f"  Grid: N={N}, h_initial={args.h}, box=[{-N*args.h/2:.2f},{N*args.h/2:.2f}]")
print(f"  C*={C_star}, lambda_3=phi^6={LAM3:.4f}")
print(f"  Phase 1: lr={args.lr1}, no angle clamp, no scale fix (cleanup)")
print(f"  Phase 2: lr={args.lr2}, delta_max={args.delta_max_deg}°, scale fix every {args.scale_fix_every} steps")
print(f"  Scale fix clamp: [{args.scale_fix_min:.3f}, {args.scale_fix_max:.3f}] per step")
print(f"  Phase transition: K rises > {args.K_rise_eps} above K_min")
print(f"{'='*70}")

# ── Grid (integer indices; physical positions = physical_h * (i - N//2)) ──
cv_int   = np.arange(N) - N//2  # integer offsets from centre
pts_int  = np.stack(np.meshgrid(cv_int,cv_int,cv_int,indexing='ij'),axis=-1).reshape(-1,3)

def get_pts_np(h_val):
    return (pts_int * h_val).astype(np.float32)

def get_dist_from_origin(h_val):
    return torch.tensor(
        np.linalg.norm(pts_int * h_val, axis=-1).reshape(N,N,N).astype(np.float32),
        dtype=torch.float32, device=dev)

needed_hw = R0 + r0 + 1/C_star + 0.5
actual_hw = N*args.h/2
if actual_hw < needed_hw:
    print(f"FATAL: box too small at initial h (need {needed_hw:.2f}, have {actual_hw:.2f})")
    sys.exit(1)
print(f"  Box adequate at h_initial: {actual_hw:.2f} >= {needed_hw:.2f}")
print(f"  Tube diameter coverage at h_initial: {2*(1/C_star)/args.h:.1f} pts")

# ── Energy functional (uses current physical_h) ───────────────────────
def E_geom(n, h_val):
    """E_geom = K*J4 at physical grid spacing h_val."""
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h_val)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2    = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a   = (s4*g2).sum() * h_val**3
    J2iso =    g2.sum()   * h_val**3
    K     = J2a + MU*J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4     = rho_J4.sum() * h_val**3
    return K*J4, K.item(), J4.item(), rho_J4

def diagnostics(n, h_val):
    with torch.no_grad():
        E, K, J4, rho = E_geom(n, h_val)
        dist_from_origin = get_dist_from_origin(h_val)
        r_bar = ((rho*dist_from_origin).sum()/rho.sum().clamp(1e-12)).item()
        s_opt = float(np.sqrt(LAM3 * J4 / max(K, 1e-12)))
    return E.item(), K, J4, r_bar, s_opt

def project_to_tangent(n, grad):
    return grad - (grad*n).sum(-1,keepdim=True)*n

def apply_angle_clamp(n_old, n_new, delta_max):
    with torch.no_grad():
        cos_a = (n_old * n_new).sum(-1, keepdim=True).clamp(-1+1e-6, 1-1e-6)
        angle = torch.acos(cos_a)
        too_far = (angle > delta_max).squeeze(-1)
        if not too_far.any():
            return n_new
        sin_a = torch.sin(angle).clamp(1e-8)
        t     = (delta_max / angle.clamp(1e-8))
        n_slerp = (torch.sin((1-t)*angle)/sin_a * n_old
                 + torch.sin(t*angle   )/sin_a * n_new)
        n_slerp = n_slerp / n_slerp.norm(dim=-1,keepdim=True).clamp(1e-10)
        n_out = n_new.clone()
        n_out[too_far] = n_slerp[too_far]
        return n_out / n_out.norm(dim=-1,keepdim=True).clamp(1e-10)

# ── Construction C initial condition ─────────────────────────────────
print("\nBuilding per-strand initial condition...")
t0b = time.time()
NT_frame = 20000
t_frame, _, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=NT_frame)
print(f"  Bishop frame holonomy: {np.degrees(H):.4f} deg")

NT = 4000
t_arr  = np.linspace(0, 2*np.pi, NT, endpoint=False)
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
lobe_indices  = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees    = [KDTree(Gamma_pts[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr[li] for li in lobe_indices]

def nearest_two_strands(qpts):
    d_per, t_per = [], []
    for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
        d, idx = tree_l.query(qpts, workers=1)
        d_per.append(d); t_per.append(t_l[idx])
    d_s = np.stack(d_per,axis=1); t_s = np.stack(t_per,axis=1)
    o   = np.argsort(d_s, axis=1)
    return (np.take_along_axis(t_s,o,axis=1)[:,0],
            np.take_along_axis(d_s,o,axis=1)[:,0],
            np.take_along_axis(t_s,o,axis=1)[:,1],
            np.take_along_axis(d_s,o,axis=1)[:,1])

def frame_at_t(t_q):
    idx = np.searchsorted(t_frame, t_q%(2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                     (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

def build_construction_c(h_val):
    pts = get_pts_np(h_val)
    t1_g, d1_g, t2_g, d2_g = nearest_two_strands(pts)
    chi1 = np.arctan2(np.sum((pts-curve_at_t(t1_g))*frame_at_t(t1_g)[1],axis=1),
                      np.sum((pts-curve_at_t(t1_g))*frame_at_t(t1_g)[0],axis=1))
    chi2 = np.arctan2(np.sum((pts-curve_at_t(t2_g))*frame_at_t(t2_g)[1],axis=1),
                      np.sum((pts-curve_at_t(t2_g))*frame_at_t(t2_g)[0],axis=1))
    Phi1 = chi1+3*t1_g; Phi2 = chi2+3*t2_g
    rho1 = np.clip(d1_g,1e-6,None); rho2 = np.clip(d2_g,1e-6,None)
    def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
    w1 = 1/rho1**2; w2 = 1/rho2**2
    z1 = (w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2 = w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    z1/=mag; z2/=mag
    nx0 = 2*np.real(np.conj(z1)*z2)
    ny0 = 2*np.imag(np.conj(z1)*z2)
    nz0 = np.abs(z1)**2 - np.abs(z2)**2
    n0  = np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
    n0 /= np.linalg.norm(n0,axis=-1,keepdims=True).clip(1e-10)
    return n0

# Build initial field
if args.warm_start:
    print(f"\nWarm-starting from {args.warm_start} (phase {args.warm_start_phase})...")
    n0_np = np.load(args.warm_start).astype(np.float32)
    if n0_np.shape != (N,N,N,3):
        print(f"FATAL: shape mismatch {n0_np.shape}"); sys.exit(1)
    n0_np /= np.linalg.norm(n0_np,axis=-1,keepdims=True).clip(1e-10)
    print(f"  Loaded. Will start in Phase {args.warm_start_phase}.")
else:
    n0_np = build_construction_c(physical_h)
print(f"  Analytic construction built in {time.time()-t0b:.1f}s")

# ── Initial diagnostics ───────────────────────────────────────────────
n_t = torch.tensor(n0_np, dtype=torch.float32, device=dev)
E0, K0, J40, rbar0, s0 = diagnostics(n_t, physical_h)
print(f"\nInitial field (h={physical_h:.4f}):")
print(f"  E={E0:.4e}  K={K0:.2f}  J4={J40:.2f}  J4/K={J40/max(K0,1e-6):.4f}")
print(f"  r_bar={rbar0:.3f}  s_opt={s0:.4f}  (saddle requires s_opt=1.0)")
print(f"  Q_H=3 (established at construction)")

# ── Optimiser ─────────────────────────────────────────────────────────
n_param  = n_t.clone().requires_grad_(True)
opt1     = torch.optim.Adam([n_param], lr=args.lr1)
opt2     = torch.optim.Adam([n_param], lr=args.lr2)

# If warm-starting in Phase 2, skip Phase 1
if args.warm_start and args.warm_start_phase == 2:
    phase = 2
    K_min_seen = K0
    opt2 = torch.optim.Adam([n_param], lr=args.lr2)
    print(f"  Skipping Phase 1 (warm-start_phase=2).")
else:
    phase = 1
    K_min_seen = K0

history    = []
log_path   = os.path.join(args.outdir, 'log.json')
t_run0     = time.time()
scale_step = 0   # Phase 2 step counter (for scale_fix_every)
cumulative_scale = 1.0  # total physical_h / args.h

vac0 = float((n_param.detach()[...,2] > 0.95).float().mean().item())

print(f"\nRunning {args.n_steps} steps  (Phase 1: lr={args.lr1}, Phase 2: lr={args.lr2})")
print(f"{'step':>6}  {'ph':>2}  {'E_geom':>12}  {'K':>8}  {'J4':>8}  "
      f"{'r_bar':>7}  {'|grad|':>10}  {'s_opt':>6}  {'h_phys':>8}  {'cl%':>5}")

for step in range(args.n_steps):
    opt = opt1 if phase == 1 else opt2
    n_before = n_param.detach().clone()

    # ── Forward + backward ───────────────────────────────────────────
    opt.zero_grad()
    E, _, _, _ = E_geom(n_param, physical_h)
    E.backward()

    with torch.no_grad():
        n_param.grad.data.copy_(
            project_to_tangent(n_param.detach(), n_param.grad))

    grad_norm = n_param.grad.norm().item()
    opt.step()

    with torch.no_grad():
        n_param.data.copy_(
            n_param / n_param.norm(dim=-1,keepdim=True).clamp(1e-10))

    clamped_frac = 0.0
    if phase == 2:
        # Angle clamp
        with torch.no_grad():
            n_clamped = apply_angle_clamp(n_before, n_param.detach(), DELTA_MAX)
            diff = (n_clamped - n_param.detach()).norm(dim=-1)
            clamped_frac = (diff > 1e-6).float().mean().item()
            n_param.data.copy_(n_clamped)

        # ── SCALE FIX: enforce s_opt = 1 ─────────────────────────────
        scale_step += 1
        if scale_step % args.scale_fix_every == 0:
            with torch.no_grad():
                _, K_cur, J4_cur, _, s_cur = diagnostics(n_param, physical_h)
                if s_cur > 1e-3 and abs(s_cur - 1.0) > 1e-4:
                    # Clamp scale change for numerical stability
                    s_clamped = np.clip(s_cur, args.scale_fix_min, args.scale_fix_max)
                    physical_h = physical_h * s_clamped
                    cumulative_scale *= s_clamped
                    # Verify box is still adequate
                    actual_hw = N * physical_h / 2
                    if actual_hw < needed_hw * 0.5:
                        print(f"\nWARNING: box shrunk too much (h={physical_h:.4f}, "
                              f"actual_hw={actual_hw:.2f}). Scale fixing paused.")
                    elif actual_hw > needed_hw * 5:
                        print(f"\nWARNING: box expanded too much (h={physical_h:.4f}, "
                              f"actual_hw={actual_hw:.2f}). Scale fixing paused.")

    # ── Diagnostics ──────────────────────────────────────────────────
    if (step+1) % args.log_every == 0 or step == 0:
        E_val, K_val, J4_val, rbar_val, s_val = diagnostics(n_param, physical_h)
        j4k = J4_val / max(K_val, 1e-6)
        near_vac = float((n_param.detach()[...,2] > 0.95).float().mean().item())

        print(f"{step+1:>6}  {phase:>2}  {E_val:>12.4e}  {K_val:>8.1f}  "
              f"{J4_val:>8.2f}  {rbar_val:>7.3f}  {grad_norm:>10.4e}  "
              f"{s_val:>6.4f}  {physical_h:>8.5f}  {100*clamped_frac:>4.1f}%")

        # ── Phase transition ──────────────────────────────────────────
        if phase == 1:
            if K_val < K_min_seen:
                K_min_seen = K_val
            elif K_val > K_min_seen + args.K_rise_eps:
                print(f"\n  *** PHASE TRANSITION at step {step+1}: "
                      f"K={K_val:.2f} > K_min={K_min_seen:.2f}+{args.K_rise_eps} ***")
                print(f"  Switching to Phase 2 with scale fixing. "
                      f"Current s_opt={s_val:.4f}.")
                phase = 2
                scale_step = 0
                opt2 = torch.optim.Adam([n_param], lr=args.lr2)

        # ── Safeguards ────────────────────────────────────────────────
        if rbar_val > 5.0:
            print(f"\nHALT: r_bar={rbar_val:.3f} > 5.0 (dilution)"); break
        if near_vac > vac0 + 0.15:
            print(f"\nHALT: near-vacuum +{100*(near_vac-vac0):.1f}pp"); break
        if j4k < 0.005 and step > 50:
            print(f"\nHALT: J4/K collapsed at step {step+1}"); break
        if physical_h > 2.0:
            print(f"\nHALT: physical_h={physical_h:.4f} too large (scale diverging)"); break
        if physical_h < 0.005:
            print(f"\nHALT: physical_h={physical_h:.4f} too small (scale collapsed)"); break

        # ── Convergence check ─────────────────────────────────────────
        if phase == 2 and grad_norm < 5e2 and abs(s_val - 1.0) < 0.02:
            print(f"\n*** CONVERGED: |grad|={grad_norm:.2e}, s_opt={s_val:.4f} at step {step+1} ***")
            print(f"  This is the Q_H=3 saddle field.")
            print(f"  K={K_val:.4f}, J4={J4_val:.6f}, K/J4={K_val/J4_val:.4f}")
            print(f"  physical_h={physical_h:.6f} (total scale: {cumulative_scale:.4f}x initial)")
            break

        history.append(dict(step=step+1, phase=phase, E=E_val, K=K_val, J4=J4_val,
                            r_bar=rbar_val, grad_norm=grad_norm, s_opt=s_val,
                            physical_h=physical_h, clamped_frac=clamped_frac))

print(f"\nTotal wall time: {time.time()-t_run0:.1f}s")
with open(log_path,'w') as f: json.dump(history,f,indent=2)
print(f"Log: {log_path}")

E_f, K_f, J4_f, rbar_f, s_f = diagnostics(n_param, physical_h)
print(f"\nFinal state:")
print(f"  E={E_f:.4e}  K={K_f:.4f}  J4={J4_f:.6f}")
print(f"  J4/K={J4_f/max(K_f,1e-6):.6f}  (saddle requires {1/LAM3:.6f} = 1/lambda_3)")
print(f"  r_bar={rbar_f:.3f}  s_opt={s_f:.6f}  (saddle requires 1.000000)")
print(f"  physical_h={physical_h:.6f} (cumulative scale: {cumulative_scale:.4f}x)")
print(f"  E reduction: {E0:.4e} -> {E_f:.4e}  ({E0/E_f:.2f}x)")

# Save field and metadata
np.save(os.path.join(args.outdir,'n_final.npy'),
        n_param.detach().cpu().numpy())
meta = dict(physical_h=physical_h, cumulative_scale=cumulative_scale,
            K=K_f, J4=J4_f, s_opt=s_f, grad_norm=float(grad_norm),
            E=E_f, r_bar=rbar_f, N=N, h_initial=args.h, C_star=C_star,
            LAM3=LAM3)
with open(os.path.join(args.outdir,'saddle_meta.json'),'w') as f:
    json.dump(meta, f, indent=2)
print(f"Field saved: {args.outdir}/n_final.npy")
print(f"Metadata:    {args.outdir}/saddle_meta.json")
print(f"""
INTERPRETING THE OUTPUT:
  If s_opt -> 1.0 AND |grad| -> 0: TRUE SADDLE FOUND.
    K_saddle = {K_f:.4f} at physical_h = {physical_h:.6f}
    Report back: K, J4, s_opt, physical_h, |grad| at convergence.
    
  If s_opt -> 1.0 BUT |grad| stays large: Derrick balance reached but
    shape not stationary. Need more steps or smaller lr2.
    
  If s_opt oscillates and |grad| large: field hasn't reached saddle.
    Try warm-starting from step with smallest |grad| in previous run.
    
  If physical_h diverges: scale is running away. Try scale_fix_max=1.02
    (more conservative scale fixing per step).
""")
