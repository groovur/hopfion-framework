#!/usr/bin/env python3
"""
qh3_trefoil_solver_jam.py  —  JAM variant of qh3_trefoil_solver_3d.py (2026-09-04)
==================================================================================
HANDLE (2) of the QGP-jam formation-energy test (notes/quark_generation_e8_ribbon_twist.md sec 19).

Identical physics to qh3_trefoil_solver_3d.py (K_fb = J2a + mu*J2iso Faddeev-Niemi gradient
flow, periodic box via torch.roll, autograd forces) with TWO changes:

  1. scipy.spatial.KDTree -> numpy brute-force nearest-point (base conda env has no scipy;
     see notes/pyenv.md). Initial-field construction only; flow physics unchanged.

  2. The "box too small" ASSERTION is relaxed to a notice. The whole point of the jam test
     is a TIGHT box: box = 2R with R = R0+r0 = 3.874 (item 1) puts the trefoil edge-to-edge
     with its periodic image => a LATTICE of jammed baryons at inter-baryon spacing d = R
     (the deconfinement threshold). Periodic BC = neighbours share field = the "overlap"
     mechanism for the completion energy.

WHY: the ISOLATED (large-box) flow ESCAPES the topological sector (J4/J4_init -> 0), which is
the long-standing saddle-finder blocker. The jam HYPOTHESIS: confinement to a cell of size R
walls off the escape route, so the tight-box flow should PRESERVE topology (J4/J4_init ~ 1) and
CONVERGE to a stable saddle. That alone un-blocks the saddle finder.

PASS/FAIL (from the note, user's reframe 2026-09-04):
  - PASS: J4/J4_init stays ~1 (no escape) AND K converges => jam stabilises the saddle.
  - The converged energy is the JAMMED energy; it should sit ABOVE the isolated core, and
    (per the sigma*R handle-1 estimate) map to ~924.8 MeV = ~98.6% of m_p, NOT the full 938.
    A full-938 (or runaway) result would be the WORRY (would imply free protons in the QGP).
  - K/J4 near phi^6 = 17.944 indicates a genuine Bogomolny-type saddle.

Marginal-jam box (d = R): box/2 = R = 3.874 -> box = 7.748. e.g. --N 36 --h 0.2152.
Denser jam (d < R, higher completion energy): shrink the box, e.g. --N 30 --h 0.20 (box=6.0).

Usage (base env has torch+numpy):
  /usr/local/Caskroom/miniconda/base/bin/python qh3_trefoil_solver_jam.py \
      --N 36 --h 0.2152 --C_star 1.5 --steps 60000 --outdir jam_marginal
"""
import numpy as np, time, argparse, os, sys
try:
    import torch
except ImportError:
    print("Use base conda env: /usr/local/Caskroom/miniconda/base/bin/python  (see notes/pyenv.md)")
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('--N',           type=int,   default=36)
ap.add_argument('--h',           type=float, default=0.2152)   # box = N*h = 7.747 = 2R (marginal jam)
ap.add_argument('--C_star',      type=float, default=1.5,
                help='Initial C* (1.5 for h>=0.20). Solver relaxes it.')
ap.add_argument('--R0',          type=float, default=3.0)
ap.add_argument('--r0',          type=float, default=0.874)
ap.add_argument('--NT',          type=int,   default=2000)
ap.add_argument('--steps',       type=int,   default=60000)
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
tag = f'qh3_jam_N{N}'
ckpt_base = os.path.join(args.outdir, tag)
log_path  = os.path.join(args.outdir, f'{tag}_log.txt')

print(f"\n{'='*66}")
print(f"  Q_H=3 Trefoil Solver -- JAM variant (tight periodic box)")
print(f"  phi^6 = {phi6:.6f}  mu* = {MU:.8f}  [NOT assumed]")
print(f"{'='*66}")
tube_r = 1/args.C_star
R = args.R0 + args.r0
print(f"  {N}^3  h={h}  box=[{-N*h/2:.2f},{N*h/2:.2f}]  device={args.device}")
print(f"  C*_init={args.C_star}  tube_radius=1/C*={tube_r:.3f}")
print(f"  Points across tube: ~{2*tube_r/h:.1f}  (need >=3)")

# --- numpy brute-force nearest point on a curve (replaces scipy KDTree) --------
def nearest(pts, Gamma, chunk=4096):
    n = pts.shape[0]
    dists = np.empty(n, np.float32); idx = np.empty(n, np.int64)
    for s in range(0, n, chunk):
        p = pts[s:s+chunk]                                  # (c,3)
        d2 = ((p[:, None, :] - Gamma[None, :, :])**2).sum(-1)   # (c,NT)
        j = d2.argmin(1)
        idx[s:s+chunk] = j
        dists[s:s+chunk] = np.sqrt(d2[np.arange(len(j)), j])
    return dists, idx

# --- Trefoil -------------------------------------------------------------------
t_arr = np.linspace(0, 2*np.pi, args.NT, endpoint=False)
R0, r0 = args.R0, args.r0
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma = np.stack([Gx,Gy,Gz], axis=1).astype(np.float32)
box = N*h
ext = max(abs(Gx).max(), abs(Gy).max())
if box/2 < ext + 2*tube_r + 1.0:
    print(f"  [JAM] Tight box: box/2={box/2:.2f} < isolated-fit {ext+2*tube_r+1:.2f}."
          f"  Trefoil (ext={ext:.2f}) is CONFINED, d~{box/2:.2f} vs R={R:.2f}. Intended.")
else:
    print(f"  Trefoil fits (max {ext:.2f} + tube {tube_r:.2f}) in box {box/2:.2f} (isolated regime).")

# --- Grid ----------------------------------------------------------------------
cv = h*(np.arange(N) - N//2 + 0.5)
pts = np.stack(np.meshgrid(cv,cv,cv, indexing='ij'), axis=-1).reshape(-1,3).astype(np.float32)

# --- Energy (K_fb = J2a + mu*J2iso) -------------------------------------------
def compute_energy(n_r, mu):
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
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    J4 = (Fxy**2 + Fxz**2 + Fyz**2).sum() * h**3
    return K, J2a, J4

def grad_K(n_t, mu):
    n_r = n_t.detach().requires_grad_(True)
    K, J2a, J4 = compute_energy(n_r, mu)
    K.backward()
    g = n_r.grad.detach()
    g = g - (g * n_t).sum(-1, keepdim=True) * n_t
    return K.item(), J2a.item(), J4.item(), -g

# --- Initialise ----------------------------------------------------------------
start_step = 1; J4_init = None
if args.resume and os.path.exists(args.resume):
    data = np.load(args.resume)
    n_np = data['n']; start_step = int(data.get('step',1))+1
    J4_init = float(data.get('J4_init', 0))
    print(f"\n  Resumed from step {start_step-1}")
else:
    print(f"\n  Building trefoil field (C*={args.C_star}) [numpy nearest]...", flush=True)
    t0b = time.time()
    dists, idx = nearest(pts, Gamma)
    rho = np.clip(dists, 1e-6, None).astype(np.float32)
    t_near = t_arr[idx].astype(np.float32)
    f_ = 2*np.arctan(rho**(-args.C_star))
    tht = 3*t_near
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

# --- Flow ----------------------------------------------------------------------
log = open(log_path, 'a')
log.write(f"# Q_H=3 JAM  N={N} h={h} box={box:.3f} C*={args.C_star}\n")
log.write(f"# step  K  K/J4  K/J2a  J4/J2a  J4/J4init  J2a  J4  t\n")
dt = args.dt; t0 = time.time()
K_prev = K0; best_lam = lam0; best_step = 1
n_best = n_np.copy()
rejects = 0

hdr = (f"  {'step':>8}  {'K':>10}  {'K/J4':>10}  {'K/J4/phi6':>10}  "
       f"{'K/J2a':>7}  {'J4/J4in':>8}  {'dt':>9}")
print(f"\n{hdr}\n  {'-'*80}")

for step in range(start_step, args.steps+1):
    K_v, J2a_v, J4_v, Force = grad_K(n_t, MU)
    lam = K_v/J4_v if J4_v > 1e-12 else float('nan')

    if step % args.print_every == 0 or step == start_step:
        el = time.time() - t0
        print(f"  {step:>8d}  {K_v:>10.2f}  {lam:>10.4f}  {lam/phi6:>10.5f}  "
              f"{K_v/J2a_v:>7.4f}  {J4_v/J4_init:>8.4f}  {dt:>9.2e}  "
              f"({el:.0f}s, {rejects}rej)")
        log.write(f"{step}  {K_v:.4f}  {lam:.6f}  {K_v/J2a_v:.5f}  "
                  f"{J4_v/J2a_v:.7f}  {J4_v/J4_init:.5f}  "
                  f"{J2a_v:.4f}  {J4_v:.6f}  {el:.1f}\n")
        log.flush()
        rejects = 0

    if abs(lam - phi6) < abs(best_lam - phi6):
        best_lam = lam; best_step = step
        n_best = n_t.detach().cpu().numpy()

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
            dt = args.dt

    if step % args.save_every == 0:
        ck = f"{ckpt_base}_{step:08d}.npz"
        np.savez(ck, n=n_t.detach().cpu().numpy(), step=step,
                 J4_init=J4_init, J2a=J2a_v, J4=J4_v, K=K_v, lam=lam)
        print(f"  [ckpt {ck}]")

log.close()

# --- Final report --------------------------------------------------------------
with torch.no_grad():
    K_f, J2a_f, J4_f = compute_energy(n_t, MU)
K_f, J2a_f, J4_f = K_f.item(), J2a_f.item(), J4_f.item()
lam_f = K_f/J4_f if J4_f > 1e-12 else float('nan')

print(f"\n{'='*66}")
print(f"  FINAL  Q_H=3 JAM  step {step}  box={box:.3f} (d~{box/2:.2f} vs R={R:.2f})")
print(f"{'='*66}")
print(f"  K            = {K_f:.6f}")
print(f"  K/J4         = {lam_f:.10f}  [lambda_3;  phi^6 = {phi6:.6f}]")
print(f"  K/J4 / phi^6 = {lam_f/phi6:.8f}")
print(f"  K/J2a        = {K_f/J2a_f:.8f}  [2phi = {2*phi:.6f}]")
print(f"  J4/J4_init   = {J4_f/J4_init:.6f}  [TOPOLOGY: ~1 = jam held it; ->0 = escaped]")
print()
if J4_f/J4_init > 0.5:
    print("  => TOPOLOGY PRESERVED: the tight box held the sector (jam PASS on stabilisation).")
else:
    print("  => TOPOLOGY LOST: the config escaped the sector even in the tight box (jam did NOT hold).")
print(f"\n  Best lambda closest to phi^6: {best_lam:.6f} at step {best_step}")
np.savez(f"{ckpt_base}_FINAL.npz", n=n_t.detach().cpu().numpy(),
         step=step, J4_init=J4_init, J2a=J2a_f, J4=J4_f, K=K_f, lam=lam_f, box=box)
np.savez(f"{ckpt_base}_BEST.npz", n=n_best, step=best_step,
         J4_init=J4_init, lam=best_lam, box=box)
print(f"  Saved: {tag}_FINAL.npz  {tag}_BEST.npz")
