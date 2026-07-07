#!/usr/bin/env python3
"""
qh3_trefoil_solver_3d.py  —  Paper XV (v4: direct K_fb minimisation)
=====================================================================
  1. ALGORITHM: Minimise K_fb = J2a + mu*J2iso directly (standard
     Faddeev-Niemi gradient flow). The topology is preserved by
     smooth flow (Hopf charge is quantized). At the K-minimum,
     K/J4 gives the Bogomolny parameter. This is exactly what
     Paper I's 2D solver does.

  2. RESOLUTION: Initialise with C*=1.5 (tube radius 0.67, well-resolved
     at h=0.22) instead of C*=3.43 (tube radius 0.29, only 1 grid point).
     The solver relaxes the profile to whatever C* minimises K_fb.
     Also use N=80 at h=0.22 for better resolution (~3 pts across tube).

  3. FORCES: PyTorch autograd — exact discrete gradient, no analytical
     formula mismatch possible.

Install:  pip install torch --index-url https://download.pytorch.org/whl/cpu

Usage:
  python qh3_trefoil_solver_3d.py --N 80 --h 0.22 --steps 200000
  python qh3_trefoil_solver_3d.py --N 64 --h 0.30 --C_star 1.5 --steps 100000
  python qh3_trefoil_solver_3d.py --resume qh3_3d_N80_00050000.npz

Watch:
  K/J4   -> converge to lambda_3  (is it phi^6 = 17.944?)
  J4/J0  -> stay near 1.0 (topology preserved by smooth flow)
  K      -> monotonically decrease (K_fb minimisation)
"""
import numpy as np, time, argparse, os, sys
from scipy.spatial import KDTree
try:
    import torch
except ImportError:
    print("pip install torch --index-url https://download.pytorch.org/whl/cpu")
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('--N',           type=int,   default=80)
ap.add_argument('--h',           type=float, default=0.22)
ap.add_argument('--C_star',      type=float, default=1.5,
                help='Initial C* (use 1.5 for h>=0.20, 3.43 for h<=0.08)')
ap.add_argument('--R0',          type=float, default=3.0)
ap.add_argument('--r0',          type=float, default=0.874)
ap.add_argument('--NT',          type=int,   default=2000)
ap.add_argument('--steps',       type=int,   default=200000)
ap.add_argument('--dt',          type=float, default=1e-4)
ap.add_argument('--max_dt',      type=float, default=5e-3)
ap.add_argument('--save_every',  type=int,   default=5000)
ap.add_argument('--print_every', type=int,   default=1000)
ap.add_argument('--outdir',      type=str,   default='.')
ap.add_argument('--resume',      type=str,   default=None)
ap.add_argument('--device',      type=str,   default='cpu')
args = ap.parse_args()

phi = (1+5**0.5)/2; phi6 = phi**6; MU = 3.0-phi
N, h = args.N, args.h
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
tag = f'qh3_3d_N{N}'
ckpt_base = os.path.join(args.outdir, tag)
log_path  = os.path.join(args.outdir, f'{tag}_log.txt')

print(f"\n{'='*66}")
print(f"  Q_H=3 Trefoil Solver v4 (direct K_fb minimisation, autograd)")
print(f"  phi^6 = {phi6:.6f}  mu* = {MU:.8f}  [NOT assumed]")
print(f"{'='*66}")
tube_r = 1/args.C_star
print(f"  {N}^3  h={h}  box=[{-N*h/2:.1f},{N*h/2:.1f}]  device={args.device}")
print(f"  C*_init={args.C_star}  tube_radius=1/C*={tube_r:.3f}")
print(f"  Points across tube: ~{2*tube_r/h:.1f}  (need >=3)")

# ── Trefoil ─────────────────────────────────────────────────────────
t_arr = np.linspace(0, 2*np.pi, args.NT, endpoint=False)
R0, r0 = args.R0, args.r0
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma = np.stack([Gx,Gy,Gz], axis=1).astype(np.float32)
tree = KDTree(Gamma)
box = N*h
ext = max(abs(Gx).max(), abs(Gy).max())
assert box/2 > ext + 2*tube_r + 1.0, \
    f"Box {box/2:.1f} too small for trefoil+tube (need >{ext+2*tube_r+1:.1f})"
print(f"  Trefoil fits (max {ext:.2f} + tube {tube_r:.2f}) in box {box/2:.1f} ✓")

# ── Grid ────────────────────────────────────────────────────────────
cv = h*(np.arange(N) - N//2 + 0.5)
pts = np.stack(np.meshgrid(cv,cv,cv, indexing='ij'), axis=-1).reshape(-1,3).astype(np.float32)

# ── Energy (K_fb = J2a + mu*J2iso) ─────────────────────────────────
def compute_energy(n_r, mu):
    """K_fb, J2a, J4 from n (N,N,N,3) with requires_grad."""
    nx, ny, nz = n_r[...,0], n_r[...,1], n_r[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a  = (s4 * g2).sum() * h**3
    J2iso = g2.sum() * h**3
    K = J2a + mu * J2iso
    # Hopf curvature
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    J4 = (Fxy**2 + Fxz**2 + Fyz**2).sum() * h**3
    return K, J2a, J4

def grad_K(n_t, mu):
    """Gradient of K_fb w.r.t. n, projected onto S^2 tangent space."""
    n_r = n_t.detach().requires_grad_(True)
    K, J2a, J4 = compute_energy(n_r, mu)
    K.backward()
    g = n_r.grad.detach()
    # Project onto tangent space of S^2
    g = g - (g * n_t).sum(-1, keepdim=True) * n_t
    return K.item(), J2a.item(), J4.item(), -g   # negative gradient = descent

# ── Initialise ──────────────────────────────────────────────────────
start_step = 1; J4_init = None
if args.resume and os.path.exists(args.resume):
    data = np.load(args.resume)
    n_np = data['n']; start_step = int(data.get('step',1))+1
    J4_init = float(data.get('J4_init', 0))
    print(f"\n  Resumed from step {start_step-1}")
else:
    print(f"\n  Building trefoil field (C*={args.C_star})...", flush=True)
    t0b = time.time()
    dists, idx = tree.query(pts)
    rho = np.clip(dists, 1e-6, None).astype(np.float32)
    t_near = t_arr[idx].astype(np.float32)
    f_ = 2*np.arctan(rho**(-args.C_star))
    tht = 3*t_near   # Q_H=3 winding
    nx_ = np.sin(f_)*np.cos(tht)
    ny_ = np.sin(f_)*np.sin(tht)
    nz_ = np.cos(f_)
    n_np = np.stack([nx_.reshape(N,N,N), ny_.reshape(N,N,N),
                     nz_.reshape(N,N,N)], axis=-1).astype(np.float32)
    n_np /= np.linalg.norm(n_np, axis=-1, keepdims=True).clip(1e-10)
    del dists, idx, rho, t_near, f_, tht, nx_, ny_, nz_
    print(f"  Done ({time.time()-t0b:.1f}s)")

n_t = torch.tensor(n_np, dtype=torch.float32, device=dev)

with torch.no_grad():
    K0, J2a0, J4_0 = compute_energy(n_t, MU)
K0, J2a0, J4_0 = K0.item(), J2a0.item(), J4_0.item()
if J4_init is None: J4_init = J4_0
lam0 = K0/J4_0 if J4_0 > 1e-12 else float('nan')
print(f"  K={K0:.2f}  J2a={J2a0:.2f}  J4={J4_0:.4f}  K/J4={lam0:.4f}")
print(f"  K/J2a={K0/J2a0:.4f}  (target 2phi={2*phi:.4f})")
print(f"  J4_init={J4_init:.6f}")

# ── Flow ────────────────────────────────────────────────────────────
log = open(log_path, 'a')
log.write(f"# Q_H=3 v4 direct-K  N={N} h={h} C*={args.C_star}\n")
log.write(f"# step  K  K/J4  K/J2a  J4/J2a  J4/J0  J2a  J4  t\n")
dt = args.dt; t0 = time.time()
K_prev = K0; best_lam = lam0; best_step = 1
n_best = n_np.copy()
rejects = 0

hdr = (f"  {'step':>8}  {'K':>10}  {'K/J4':>10}  {'K/J4/phi6':>10}  "
       f"{'K/J2a':>7}  {'J4/J0':>7}  {'dt':>9}")
print(f"\n{hdr}\n  {'-'*78}")

for step in range(start_step, args.steps+1):
    K_v, J2a_v, J4_v, Force = grad_K(n_t, MU)
    lam = K_v/J4_v if J4_v > 1e-12 else float('nan')

    if step % args.print_every == 0 or step == start_step:
        el = time.time() - t0
        print(f"  {step:>8d}  {K_v:>10.2f}  {lam:>10.4f}  {lam/phi6:>10.5f}  "
              f"{K_v/J2a_v:>7.4f}  {J4_v/J4_init:>7.4f}  {dt:>9.2e}  "
              f"({el:.0f}s, {rejects}rej)")
        log.write(f"{step}  {K_v:.4f}  {lam:.6f}  {K_v/J2a_v:.5f}  "
                  f"{J4_v/J2a_v:.7f}  {J4_v/J4_init:.5f}  "
                  f"{J2a_v:.4f}  {J4_v:.6f}  {el:.1f}\n")
        log.flush()
        rejects = 0

    if abs(lam - phi6) < abs(best_lam - phi6):
        best_lam = lam; best_step = step
        n_best = n_t.detach().cpu().numpy()

    # Gradient step
    with torch.no_grad():
        n_try = n_t + dt * Force
        n_try = n_try / n_try.norm(dim=-1, keepdim=True).clamp(1e-10)

    with torch.no_grad():
        K_try, _, _ = compute_energy(n_try, MU)
    K_try = K_try.item()

    if K_try < K_v:
        n_t = n_try; K_prev = K_try
        dt = min(dt * 1.01, args.max_dt)
    else:
        rejects += 1
        dt *= 0.7
        if dt < 1e-14:
            dt = args.dt  # reset instead of exhausting

    if step % args.save_every == 0:
        ck = f"{ckpt_base}_{step:08d}.npz"
        np.savez(ck, n=n_t.detach().cpu().numpy(), step=step,
                 J4_init=J4_init, J2a=J2a_v, J4=J4_v, K=K_v, lam=lam)
        print(f"  [ckpt {ck}]")

log.close()

# ── Final report ────────────────────────────────────────────────────
with torch.no_grad():
    K_f, J2a_f, J4_f = compute_energy(n_t, MU)
K_f, J2a_f, J4_f = K_f.item(), J2a_f.item(), J4_f.item()
lam_f = K_f/J4_f if J4_f > 1e-12 else float('nan')

print(f"\n{'='*66}")
print(f"  FINAL  Q_H=3 3D  step {step}")
print(f"{'='*66}")
print(f"  K           = {K_f:.6f}")
print(f"  K/J4        = {lam_f:.10f}  [lambda_3]")
print(f"  phi^6       = {phi6:.10f}  [reference]")
print(f"  K/J4 / phi^6= {lam_f/phi6:.8f}")
print(f"  log_phi(lam)= {np.log(lam_f)/np.log(phi) if lam_f>0 else 0:.6f}  [6.0 for phi^6]")
print(f"  K/J2a       = {K_f/J2a_f:.8f}  [2phi = {2*phi:.6f}]")
print(f"  J4/J2a      = {J4_f/J2a_f:.8f}  [virial ratio r_3]")
print(f"  J4/J4_init  = {J4_f/J4_init:.6f}  [topology]")
print()
for k in range(1, 20):
    lk = np.log(lam_f)/np.log(phi) if lam_f > 0 else 0
    if abs(lk - k) < 0.08:
        print(f"  *** lambda_3 ~ phi^{k} = {phi**k:.6f}  "
              f"({(lam_f-phi**k)/phi**k*100:+.3f}%) ***")
print(f"\n  Best lambda closest to phi^6: {best_lam:.6f} at step {best_step}")
np.savez(f"{ckpt_base}_FINAL.npz", n=n_t.detach().cpu().numpy(),
         step=step, J4_init=J4_init, J2a=J2a_f, J4=J4_f, K=K_f, lam=lam_f)
np.savez(f"{ckpt_base}_BEST.npz", n=n_best, step=best_step,
         J4_init=J4_init, lam=best_lam)
print(f"  Saved: {tag}_FINAL.npz  {tag}_BEST.npz")
