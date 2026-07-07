#!/usr/bin/env python3
"""
gradient_flow_constrained_q.py
===================================================================
Generic, q_pol-parameterised version of
papers/src_paper16/gradient_flow_constrained.py's topology-protected
constrained gradient flow.

WHY THIS EXISTS
---------------
papers/src_paper18/hopf_link_construction_v2.py's plain (unconstrained
Adam) gradient flow was found to let the field escape its starting
topological sector: per-grid-point rotations of O(10-20 radians) in a
single step are large enough to jump over the topological wall on the
discretised grid, producing the "radiating blob" energy collapse seen
in the q_pol=3-5 sweep (E dropping ~98% by step 2400, visible as a
detached, drifting cluster of vector energy in the rendered field, in
BOTH the target construction and the dedicated reference trefoil).
This is very likely the same numerical failure mode already diagnosed
and fixed in src_paper16/gradient_flow_constrained.py (built for the
q_pol=3 trefoil specifically): a two-phase schedule (large-lr cleanup,
then small-lr saddle approach) combined with a post-step angle-clamp
(slerp projection capping per-point rotation per step) that keeps the
optimiser from jumping the topological wall.

This script generalises that same protected optimiser to any real
q_pol (reusing hopf_link_construction_v2.py's general 2-arc
build_torus_knot_field() construction, NOT the odd-integer-only N-lobe
one -- most q_pol values of interest here, e.g. the golden-section
points ~3.764/~4.236 and the q=4 link sector, are non-integer or even
and cannot be built by the N-lobe method at all), so the q_pol=3-5
stability question can be tested with a tool that actually protects
the topological sector, rather than one that silently escapes it
partway through the run. Direct inspection of relaxed fields (see
README) already confirmed the general construction converges to the
correct topology despite an unrelaxed-frame symmetry artifact, so
there is no correctness reason to prefer the N-lobe method even where
it would apply.

DIFFERENCE FROM src_paper16/gradient_flow_constrained.py
----------------------------------------------------------
Only ONE field (the target q_pol construction) is built and relaxed
per run -- the original always also built+relaxed the fixed q=3
reference trefoil for comparison every single time, doubling cost.
That reference behavior is already on record from prior runs; dropping
it here roughly halves per-point wall time, which matters given these
runs can take hours to days depending on q_pol and n_steps.

CAVEAT: the Q_H=n identification (T(2,q) <-> Q_H=q) used in comments
below follows the pattern established for q=1,2,3 (Papers I/XVII/
XVIII); q=4,5 are acknowledged but topologically unestablished sectors
in the framework (Paper XVIII, rem:higher_charge_incoming) -- nothing
here proves stability/instability of any such sector, it only tests
whether THIS energy functional's gradient flow, run correctly (without
numerical topology escape), holds a given q_pol construction together
over long timescales, and if so, near what energy.
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--q_pol',        type=float, required=True,
                help='Curve-geometry parameter, any real value >= 2. '
                     'Uses the general 2-arc T(2,q_pol) construction '
                     '(valid for integer, non-integer, and even q).')
ap.add_argument('--alpha_wind',   type=float, default=1.0,
                help='Phase-winding interpolation (0=Q_H=2 field, '
                     '1=full q_pol-type winding). Default 1.0 (full '
                     'commitment), matching the q_pol sweep so far.')
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
ap.add_argument('--outdir',       type=str,   default='gf_constrained_q')
ap.add_argument('--seed',         type=int,   default=0)
ap.add_argument('--device',       type=str,   default='cpu')
ap.add_argument('--snapshots',    type=str,   default='',
                help='Comma-separated list of steps at which to save n_<step>.npy '
                     'snapshots. The Phase 1 exit step and final step are '
                     'always saved.')
args = ap.parse_args()

q_pol = args.q_pol
alpha = args.alpha_wind
wind_exp = 2.0 + alpha  # matches hopf_link_construction_v2.py's convention;
                        # alpha=1.0 default -> wind_exp=3.0 (full commitment,
                        # same as every q_pol point in the earlier sweeps)
if q_pol < 2.0:
    print(f"FATAL: q_pol must be >= 2 (got {q_pol})")
    sys.exit(1)

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
print(f"  CONSTRAINED GRADIENT FLOW, q_pol={q_pol}, wind_exp={wind_exp}  (E_geom = K*J4, topology-safe)")
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

# ── Energy functional (identical to src_paper16/gradient_flow_constrained.py) ──
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
    fname = os.path.join(args.outdir, f'n_{step:06d}.npy')
    np.save(fname, n_param.detach().cpu().numpy())
    tag = f' [{label}]' if label else ''
    print(f"  >> SNAPSHOT saved: {fname}{tag}")

# ── General 2-arc construction (build_torus_knot_field, any real q_pol) ─
print(f"\nBuilding q_pol={q_pol} initial condition (general 2-arc)...")
t0b = time.time()

NT = 8000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0 + r0*np.cos(q_pol*t_arr))*np.cos(2*t_arr)
Gy = (R0 + r0*np.cos(q_pol*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(q_pol*t_arr)
G_full = np.stack([Gx, Gy, Gz], axis=1)

mask1 = t_arr < np.pi
mask2 = t_arr >= np.pi
t1_arr = t_arr[mask1]; G1 = G_full[mask1]
t2_arr = t_arr[mask2]; G2 = G_full[mask2]
tau1_arr = 2*t1_arr
tau2_arr = 2*(t2_arr - np.pi)

tree1 = KDTree(G1)
tree2 = KDTree(G2)

dt = t_arr[1] - t_arr[0]
dGdt   = np.gradient(G_full, dt, axis=0)
T_arr  = dGdt / np.linalg.norm(dGdt, axis=1, keepdims=True).clip(1e-10)
d2Gdt2 = np.gradient(dGdt, dt, axis=0)
N_raw  = d2Gdt2 - (d2Gdt2*T_arr).sum(1,keepdims=True)*T_arr
N_norm = np.linalg.norm(N_raw, axis=1, keepdims=True).clip(1e-10)
N_arr  = N_raw / N_norm
B_arr  = np.cross(T_arr, N_arr)

N1_1 = N_arr[mask1];  N2_1 = B_arr[mask1]
N1_2 = N_arr[mask2];  N2_2 = B_arr[mask2]

d1, idx1 = tree1.query(pts_np, workers=1)
d2, idx2 = tree2.query(pts_np, workers=1)

def chi_angle(pts, curve_pts, nidx, N1_arr, N2_arr):
    rel  = pts - curve_pts[nidx]
    chi  = np.arctan2((rel*N2_arr[nidx]).sum(1), (rel*N1_arr[nidx]).sum(1))
    return chi

chi_1 = chi_angle(pts_np, G1, idx1, N1_1, N2_1)
chi_2 = chi_angle(pts_np, G2, idx2, N1_2, N2_2)
tau_1 = tau1_arr[idx1]
tau_2 = tau2_arr[idx2]

Phi1 = chi_1 + wind_exp * tau_1
Phi2 = chi_2 + wind_exp * tau_2

def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
rho1 = np.clip(d1, 1e-6, None); rho2 = np.clip(d2, 1e-6, None)
f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
w1 = 1/rho1**2;       w2 = 1/rho2**2

z1_np = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
z2_np = w1*np.sin(f1/2)*np.exp(1j*Phi1) + w2*np.sin(f2/2)*np.exp(1j*Phi2)
mag = np.sqrt(np.abs(z1_np)**2 + np.abs(z2_np)**2)
z1_np /= mag; z2_np /= mag

nx0 = 2*np.real(np.conj(z1_np)*z2_np)
ny0 = 2*np.imag(np.conj(z1_np)*z2_np)
nz0 = np.abs(z1_np)**2 - np.abs(z2_np)**2
n0_np = np.stack([nx0, ny0, nz0], axis=-1).reshape(N, N, N, 3).astype(np.float32)
n0_np /= np.linalg.norm(n0_np, axis=-1, keepdims=True).clip(1e-10)
print(f"  Construction (q_pol={q_pol}) built in {time.time()-t0b:.1f}s")

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
vac0 = ((n_t[...,2]>0.95).float().mean()).item()
print(f"\nInitial field (q_pol={q_pol}):")
print(f"  E={E0:.4e}  K={K0:.2f}  J4={J40:.2f}  J4/K={J40/max(K0,1e-6):.4f}  r_bar={rbar0:.3f}")

# ── Optimiser and state ───────────────────────────────────────────────
n_param = n_t.clone().requires_grad_(True)
opt1 = torch.optim.Adam([n_param], lr=args.lr1)
opt2 = torch.optim.Adam([n_param], lr=args.lr2)

phase      = 1
K_min_seen = K0
history    = []
log_path   = os.path.join(args.outdir, 'log.json')
t_run0     = time.time()
phase1_exit_step = None

print(f"\nRunning {args.n_steps} steps  (Phase 1: lr={args.lr1}, Phase 2: lr={args.lr2})")
print(f"{'step':>6}  {'ph':>3}  {'E_geom':>12}  {'K':>8}  {'J4':>8}  "
      f"{'r_bar':>7}  {'|grad|':>10}  {'J4/K':>7}  {'clamped%':>9}")

for step in range(args.n_steps):
    opt = opt1 if phase == 1 else opt2
    n_before = n_param.detach().clone()

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

    current_step = step + 1

    if current_step in snapshot_steps:
        save_snapshot(n_param, current_step)

    if current_step % args.log_every == 0 or step == 0:
        E_val, K_val, J4_val, rbar_val = diagnostics(n_param)
        j4k = J4_val / max(K_val, 1e-6)
        near_vac = (n_param.detach()[...,2] > 0.95).float().mean().item()

        print(f"{current_step:>6}  {phase:>3}  {E_val:>12.4e}  {K_val:>8.1f}  "
              f"{J4_val:>8.2f}  {rbar_val:>7.3f}  {grad_norm:>10.4e}  "
              f"{j4k:>7.4f}  {100*clamped_frac:>8.2f}%")

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
                save_snapshot(n_param, current_step, label='Phase1_exit')
                opt2 = torch.optim.Adam([n_param], lr=args.lr2)

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
print(f"E reduction: {E0:.4e} -> {E_f:.4e}  ({E0/max(E_f,1e-12):.2f}x)")

final_step = min(args.n_steps, step+1)
save_snapshot(n_param, final_step, label='final')
np.save(os.path.join(args.outdir,'n_final.npy'),
        n_param.detach().cpu().numpy())
print(f"Field also saved as: {args.outdir}/n_final.npy")
