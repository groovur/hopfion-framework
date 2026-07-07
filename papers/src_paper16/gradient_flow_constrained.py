#!/usr/bin/env python3
"""
gradient_flow_constrained.py
============================
Constrained gradient-flow minimisation of E_geom = K*J4 in the Q_H=3
topological sector, using per-point step-size control (Option A) to
prevent discrete topology change events.

TOPOLOGY PROTECTION MECHANISM:
  The prior unconstrained run (gradient_flow_solver.py) escaped from
  Q_H=3 around step 711: the Adam optimiser was making per-grid-point
  rotations of O(20 radians) per step, large enough to jump over the
  topological wall in a single move. The fix is a POST-STEP PROJECTION:
  after each Adam update, compute the arc-distance each grid point
  actually moved on S^2 (angle = arccos(n_old . n_new)), and slerp any
  over-rotating points back to exactly delta_max. This is the standard
  projected-gradient method: the optimiser runs unconstrained, then the
  result is projected back onto the feasible set.

TWO-PHASE SCHEDULE:
  Phase 1 (cleanup): large lr, no angle clamp. Removes excess gradient
    energy from the initial Construction C ansatz efficiently. Ends
    automatically when K stops falling (the K-minimum transition).
  Phase 2 (saddle approach): smaller lr, tight angle clamp (delta_max).
    Careful approach to the Q_H=3 energy minimum with topology preserved.
  Transition detected automatically: K_current > K_min_seen + epsilon.

SNAPSHOT SYSTEM:
  --snapshots "step1,step2,..."  saves n_<step>.npy at the listed steps.
  Also saves automatically at Phase 1 exit (K_min moment) and at the
  final step regardless of --snapshots. Each snapshot file is named
  n_<step>.npy (e.g. n_0850.npy for step 850, the Phase 1 exit).
  The Phase 1 exit snapshot is the cleanest view of the Q_H=3 family
  boundary; Phase 2 snapshots show the vacuum-drift evolution through
  the 12 fraying sites.

CONVERGENCE CRITERION:
  Both |grad E| approaching zero AND J4/K stabilising. The saddle
  requires d/dn(K*J4) = 0, i.e. K*(dJ4/dn) + J4*(dK/dn) = 0.
  J4/K flattening is a necessary (not sufficient) condition.
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--N',            type=int,   default=64)
ap.add_argument('--h',            type=float, default=0.175)
ap.add_argument('--R0',           type=float, default=3.0)
ap.add_argument('--r0',           type=float, default=0.874)
ap.add_argument('--C_star',       type=float, default=2.5062)
ap.add_argument('--n_steps',      type=int,   default=2000)
ap.add_argument('--lr1',          type=float, default=3e-4,
                help='Phase 1 learning rate (cleanup, no angle clamp)')
ap.add_argument('--lr2',          type=float, default=1e-5,
                help='Phase 2 learning rate (saddle approach, with clamp)')
ap.add_argument('--delta_max_deg',type=float, default=9.0,
                help='Max rotation per grid point per step in Phase 2 (degrees)')
ap.add_argument('--K_rise_eps',   type=float, default=2.0,
                help='K rise above K_min to trigger Phase 2 transition')
ap.add_argument('--log_every',    type=int,   default=10)
ap.add_argument('--warm_start',   type=str,   default=None,
                help='path to n_final.npy to continue from')
ap.add_argument('--outdir',       type=str,   default='gf_constrained')
ap.add_argument('--seed',         type=int,   default=0)
ap.add_argument('--device',       type=str,   default='cpu')
ap.add_argument('--snapshots',    type=str,   default='',
                help='Comma-separated list of steps at which to save n_<step>.npy '
                     'snapshots. E.g. --snapshots "200,500,1000,2000". '
                     'The Phase 1 exit step and final step are always saved.')
args = ap.parse_args()

# Parse snapshot steps
snapshot_steps = set()
if args.snapshots.strip():
    for s in args.snapshots.split(','):
        s = s.strip()
        if s:
            snapshot_steps.add(int(s))

torch.manual_seed(args.seed)
np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

PHI    = (1+5**0.5)/2
MU     = 3.0 - PHI
N, h   = args.N, args.h
R0, r0, C_star = args.R0, args.r0, args.C_star
dev    = torch.device(args.device)
DELTA_MAX = float(np.radians(args.delta_max_deg))

print(f"{'='*70}")
print(f"  CONSTRAINED GRADIENT FLOW  (E_geom = K*J4, topology-safe)")
print(f"  Grid: N={N}, h={h}, box=[{-N*h/2:.2f},{N*h/2:.2f}]")
print(f"  C*={C_star}")
print(f"  Phase 1: lr={args.lr1}, no angle clamp (cleanup)")
print(f"  Phase 2: lr={args.lr2}, delta_max={args.delta_max_deg}° per point")
print(f"  Phase transition: K rises > {args.K_rise_eps} above its minimum")
if snapshot_steps:
    print(f"  Snapshots at steps: {sorted(snapshot_steps)} (+ Phase 1 exit + final)")
else:
    print(f"  Snapshots: Phase 1 exit + final only")
print(f"{'='*70}")

# ── Grid ──────────────────────────────────────────────────────────────
cv      = h*(np.arange(N) - N//2 + 0.5)
pts_np  = np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),axis=-1).reshape(-1,3).astype(np.float32)
dist_from_origin = torch.tensor(
    np.linalg.norm(pts_np,axis=-1).reshape(N,N,N), dtype=torch.float32, device=dev)

needed_hw = R0 + r0 + 1/C_star + 0.5
actual_hw = N*h/2
if actual_hw < needed_hw:
    print(f"FATAL: box too small (need {needed_hw:.2f}, have {actual_hw:.2f})")
    sys.exit(1)
print(f"  Box adequate: half-width {actual_hw:.2f} >= {needed_hw:.2f}")
print(f"  Tube diameter coverage at C*={C_star}: {2*(1/C_star)/h:.1f} pts")

# ── Energy functional ─────────────────────────────────────────────────
def E_geom(n):
    nx,ny,nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u,a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2    = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a   = (s4*g2).sum() * h**3
    J2iso =    g2.sum()   * h**3
    K     = J2a + MU*J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4    = rho_J4.sum() * h**3
    return K*J4, K.item(), J4.item(), rho_J4

def diagnostics(n):
    with torch.no_grad():
        E,K,J4,rho = E_geom(n)
        r_bar = ((rho*dist_from_origin).sum()/rho.sum().clamp(1e-12)).item()
    return E.item(), K, J4, r_bar

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

def save_snapshot(n_param, step, label=''):
    """Save n_<step>.npy and log a line to stdout."""
    fname = os.path.join(args.outdir, f'n_{step:06d}.npy')
    np.save(fname, n_param.detach().cpu().numpy())
    tag = f' [{label}]' if label else ''
    print(f"  >> SNAPSHOT saved: {fname}{tag}")

# ── Per-strand Construction C initial condition ───────────────────────
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

t1_g, d1_g, t2_g, d2_g = nearest_two_strands(pts_np)
chi1 = np.arctan2(np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[1],axis=1),
                  np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[0],axis=1))
chi2 = np.arctan2(np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[1],axis=1),
                  np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[0],axis=1))
Phi1 = chi1+3*t1_g; Phi2 = chi2+3*t2_g
rho1 = np.clip(d1_g,1e-6,None); rho2 = np.clip(d2_g,1e-6,None)

def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
w1 = 1/rho1**2; w2 = 1/rho2**2
z1_np = (w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
z2_np = w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
mag   = np.sqrt(np.abs(z1_np)**2+np.abs(z2_np)**2)
z1_np/=mag; z2_np/=mag
nx0 = 2*np.real(np.conj(z1_np)*z2_np)
ny0 = 2*np.imag(np.conj(z1_np)*z2_np)
nz0 = np.abs(z1_np)**2 - np.abs(z2_np)**2
n0_np = np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
n0_np /= np.linalg.norm(n0_np,axis=-1,keepdims=True).clip(1e-10)
print(f"  Analytic construction built in {time.time()-t0b:.1f}s")

if args.warm_start:
    print(f"\nWarm-starting from {args.warm_start}...")
    n0_np = np.load(args.warm_start).astype(np.float32)
    if n0_np.shape != (N,N,N,3):
        print(f"FATAL: shape mismatch {n0_np.shape} vs ({N},{N},{N},3)"); sys.exit(1)
    n0_np /= np.linalg.norm(n0_np,axis=-1,keepdims=True).clip(1e-10)
    print(f"  Loaded.")

# ── Initial diagnostics ───────────────────────────────────────────────
n_t   = torch.tensor(n0_np, dtype=torch.float32, device=dev)
E0, K0, J40, rbar0 = diagnostics(n_t)
_, _, vac0 = (None, None, ((n_t[...,2]>0.95).float().mean()).item())
print(f"\nInitial field:")
print(f"  E={E0:.4e}  K={K0:.2f}  J4={J40:.2f}  J4/K={J40/max(K0,1e-6):.4f}  r_bar={rbar0:.3f}")
print(f"  Q_H=3 (established at construction)")

# ── Optimiser and state ───────────────────────────────────────────────
n_param = n_t.clone().requires_grad_(True)
opt1 = torch.optim.Adam([n_param], lr=args.lr1)
opt2 = torch.optim.Adam([n_param], lr=args.lr2)

phase      = 1
K_min_seen = K0
history    = []
log_path   = os.path.join(args.outdir, 'log.json')
t_run0     = time.time()
phase1_exit_step = None  # will be set at phase transition

print(f"\nRunning {args.n_steps} steps  (Phase 1: lr={args.lr1}, Phase 2: lr={args.lr2})")
print(f"{'step':>6}  {'ph':>3}  {'E_geom':>12}  {'K':>8}  {'J4':>8}  "
      f"{'r_bar':>7}  {'|grad|':>10}  {'J4/K':>7}  {'clamped%':>9}")

for step in range(args.n_steps):
    opt = opt1 if phase == 1 else opt2
    n_before = n_param.detach().clone()

    # ── Forward + backward ───────────────────────────────────────────
    opt.zero_grad()
    E, _, _, _ = E_geom(n_param)
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
        with torch.no_grad():
            n_clamped = apply_angle_clamp(n_before, n_param.detach(), DELTA_MAX)
            diff = (n_clamped - n_param.detach()).norm(dim=-1)
            clamped_frac = (diff > 1e-6).float().mean().item()
            n_param.data.copy_(n_clamped)

    current_step = step + 1  # 1-indexed, matches log output

    # ── Snapshot at requested steps ───────────────────────────────────
    if current_step in snapshot_steps:
        save_snapshot(n_param, current_step)

    # ── Diagnostics ──────────────────────────────────────────────────
    if current_step % args.log_every == 0 or step == 0:
        E_val, K_val, J4_val, rbar_val = diagnostics(n_param)
        j4k = J4_val / max(K_val, 1e-6)
        near_vac = (n_param.detach()[...,2] > 0.95).float().mean().item()

        print(f"{current_step:>6}  {phase:>3}  {E_val:>12.4e}  {K_val:>8.1f}  "
              f"{J4_val:>8.2f}  {rbar_val:>7.3f}  {grad_norm:>10.4e}  "
              f"{j4k:>7.4f}  {100*clamped_frac:>8.2f}%")

        # ── Phase transition check ────────────────────────────────────
        if phase == 1:
            if K_val < K_min_seen:
                K_min_seen = K_val
            elif K_val > K_min_seen + args.K_rise_eps:
                print(f"\n  *** PHASE TRANSITION at step {current_step}: "
                      f"K={K_val:.2f} > K_min={K_min_seen:.2f}+{args.K_rise_eps} ***")
                print(f"  Switching to Phase 2: lr={args.lr2}, "
                      f"delta_max={args.delta_max_deg}°")
                phase = 2
                phase1_exit_step = current_step
                # Save Phase 1 exit snapshot automatically
                save_snapshot(n_param, current_step, label='Phase1_exit')
                opt2 = torch.optim.Adam([n_param], lr=args.lr2)

        # ── Dilution safeguards ───────────────────────────────────────
        if rbar_val > 4.5:
            print(f"\nHALT: r_bar={rbar_val:.3f} > 4.5 (dilution step {current_step})")
            break
        if near_vac > vac0 + 0.10:
            print(f"\nHALT: near-vacuum grew +{100*(near_vac-vac0):.1f}pp (dilution)")
            break
        if j4k < 0.01 and step > 50:
            print(f"\nHALT: J4/K={j4k:.5f} collapsed (topology lost step {current_step})")
            break

        history.append(dict(step=current_step, phase=phase, E=E_val, K=K_val,
                            J4=J4_val, r_bar=rbar_val, grad_norm=grad_norm,
                            J4_over_K=j4k, clamped_frac=clamped_frac))

print(f"\nTotal wall time: {time.time()-t_run0:.1f}s")
with open(log_path,'w') as f: json.dump(history,f,indent=2)
print(f"Log: {log_path}")

E_f, K_f, J4_f, rbar_f = diagnostics(n_param)
print(f"\nFinal: E={E_f:.4e}  K={K_f:.2f}  J4={J4_f:.2f}  "
      f"J4/K={J4_f/max(K_f,1e-6):.4f}  r_bar={rbar_f:.3f}")
print(f"E reduction: {E0:.4e} -> {E_f:.4e}  ({E0/E_f:.2f}x)")

# Always save final snapshot
final_step = min(args.n_steps, step+1)
save_snapshot(n_param, final_step, label='final')
np.save(os.path.join(args.outdir,'n_final.npy'),
        n_param.detach().cpu().numpy())
print(f"Field also saved as: {args.outdir}/n_final.npy")
