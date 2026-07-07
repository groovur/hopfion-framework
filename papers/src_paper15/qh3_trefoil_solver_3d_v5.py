#!/usr/bin/env python3
"""
qh3_trefoil_solver_3d_v5.py  v5 — Definitive lambda_3 test
=========================================================
Minimises E_test = K_fb + lambda * J4 with lambda FIXED.
  The saddle of K×J4 satisfies: dK/dn + lambda_3 * dJ4/dn = 0
  This is exactly the stationarity condition for E_test = K + lambda*J4
  with lambda = lambda_3.

  If lambda_3 = phi^6: the minimum of E_test at lambda=phi^6 should
  have K/J4 = phi^6 and K/J2a = 2phi (virial). Starting from the BEST
  checkpoint (step 385 with K/J4=17.90≈phi^6), the force is nearly
  zero and the field should barely move — converging to the nearby saddle.

  If lambda_3 ≠ phi^6: the field will drift from step 385, and the
  converged K/J4 will reveal the true lambda_3.

WHAT TO WATCH:
  K/J4   → should converge to phi^6 if lambda_3=phi^6, else to lambda_3
  K/J2a  → should converge to 2phi=3.236 at the true saddle
  J4/J0  → should stabilise > 0 (topology preserved at saddle)
  E_test → should decrease monotonically

USAGE:
  # Definitive test from BEST checkpoint:
  python qh3_trefoil_solver_3d.py \\
      --resume ./qh3_best_resume/qh3_3d_N80_BEST.npz \\
      --lam 17.9443 --N 80 --h 0.22 --steps 20000

  # Fresh run (also valid — will evolve to saddle if lambda is right):
  python qh3_trefoil_solver_3d_v5.py \\
      --N 80 --h 0.22 --C_star 1.5 --lam 17.9443 --steps 10000

  # Test other lambda values:
  python qh3_trefoil_solver_3d_v5.py --resume BEST.npz --lam 11.090  # phi^5
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
ap.add_argument('--C_star',      type=float, default=1.5)
ap.add_argument('--R0',          type=float, default=3.0)
ap.add_argument('--r0',          type=float, default=0.874)
ap.add_argument('--NT',          type=int,   default=2000)
ap.add_argument('--lam',         type=float, default=None,
                help='Fixed Bogomolny parameter (default: phi^6 = 17.9443)')
ap.add_argument('--steps',       type=int,   default=20000)
ap.add_argument('--dt',          type=float, default=1e-4)
ap.add_argument('--max_dt',      type=float, default=5e-3)
ap.add_argument('--save_every',  type=int,   default=2000)
ap.add_argument('--print_every', type=int,   default=500)
ap.add_argument('--outdir',      type=str,   default='.')
ap.add_argument('--resume',      type=str,   default=None)
ap.add_argument('--device',      type=str,   default='cpu')
args = ap.parse_args()

phi = (1+5**0.5)/2; phi6=phi**6; MU=3-phi
N, h = args.N, args.h
LAM = args.lam if args.lam is not None else phi6
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
tag = f'qh3_v5_lam{LAM:.3f}_N{N}'
log_path = os.path.join(args.outdir, f'{tag}_log.txt')
ckpt_base= os.path.join(args.outdir, tag)

print(f"\n{'='*65}")
print(f"  Q_H=3 Solver v5 — Definitive lambda test")
print(f"  Minimising E = K_fb + {LAM:.6f} * J4  [lambda FIXED]")
print(f"  phi^6 = {phi6:.6f}  {'← testing this' if abs(LAM-phi6)<0.01 else f'← testing phi^k={np.log(LAM)/np.log(phi):.2f}'}")
print(f"  At saddle: K/J4 should converge to {LAM:.4f}")
print(f"  At virial: K/J2a should converge to 2phi = {2*phi:.4f}")
print(f"{'='*65}")

# ── Trefoil ────────────────────────────────────────────────────────
t_arr = np.linspace(0,2*np.pi,args.NT,endpoint=False)
Gx=(args.R0+args.r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy=(args.R0+args.r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz=args.r0*np.sin(3*t_arr)
Gamma=np.stack([Gx,Gy,Gz],axis=1).astype(np.float32)
tree=KDTree(Gamma)
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),axis=-1).reshape(-1,3).astype(np.float32)

# ── Energy ─────────────────────────────────────────────────────────
def compute(n_r, mu):
    nx,ny,nz=n_r[...,0],n_r[...,1],n_r[...,2]
    s4=(1-nz**2).clamp(0,1)**2
    def cd(u,a): return (torch.roll(u,-1,a)-torch.roll(u,1,a))/(2*h)
    nxx,nxy,nxz=cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz=cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz=cd(nz,0),cd(nz,1),cd(nz,2)
    g2=nxx**2+nxy**2+nxz**2+nyx**2+nyy**2+nyz**2+nzx**2+nzy**2+nzz**2
    J2a=(s4*g2).sum()*h**3; J2iso=g2.sum()*h**3; K=J2a+mu*J2iso
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*h**3
    return K, J2a, J4

def step_force(n_t, mu, lam):
    """Gradient of K + lam*J4 projected onto S^2 tangent space."""
    n_r=n_t.detach().requires_grad_(True)
    K,J2a,J4=compute(n_r,mu)
    E_test=K+lam*J4
    E_test.backward()
    g=n_r.grad.detach()
    g=g-(g*n_t).sum(-1,keepdim=True)*n_t  # project to tangent space
    return K.item(),J2a.item(),J4.item(),E_test.item(),-g

# ── Initialise ─────────────────────────────────────────────────────
start_step=1; J4_init=None
if args.resume and os.path.exists(args.resume):
    data=np.load(args.resume)
    n_np=data['n']; start_step=int(data.get('step',1))+1
    J4_init=float(data.get('J4_init',0))
    print(f"  Resumed from step {start_step-1}")
else:
    print(f"  Building trefoil field (C*={args.C_star})...",flush=True)
    dists,idx=tree.query(pts)
    rho=np.clip(dists,1e-6,None).astype(np.float32)
    t_near=t_arr[idx].astype(np.float32)
    f_=2*np.arctan(rho**(-args.C_star)); tht=3*t_near
    nx_=np.sin(f_)*np.cos(tht); ny_=np.sin(f_)*np.sin(tht); nz_=np.cos(f_)
    n_np=np.stack([nx_.reshape(N,N,N),ny_.reshape(N,N,N),
                   nz_.reshape(N,N,N)],axis=-1).astype(np.float32)
    n_np/=np.linalg.norm(n_np,axis=-1,keepdims=True).clip(1e-10)

n_t=torch.tensor(n_np,dtype=torch.float32,device=dev)
with torch.no_grad():
    K0,J2a0,J4_0=compute(n_t,MU)
K0,J2a0,J4_0=K0.item(),J2a0.item(),J4_0.item()
if J4_init is None: J4_init=J4_0
E0=K0+LAM*J4_0
lam0=K0/J4_0 if J4_0>1e-12 else float('nan')
print(f"  K={K0:.3f}  J2a={J2a0:.3f}  J4={J4_0:.4f}")
print(f"  K/J4={lam0:.4f} (vs lambda={LAM:.4f}; diff={abs(lam0-LAM)/LAM*100:.2f}%)")
print(f"  K/J2a={K0/J2a0:.4f} (vs 2phi={2*phi:.4f})")
print(f"  E_test = K + {LAM:.4f}*J4 = {E0:.3f}")

# ── Flow ──────────────────────────────────────────────────────────
log=open(log_path,'a')
log.write(f"# v5 lambda={LAM:.6f}  N={N} h={h}\n")
log.write(f"# step  K/J4  K/J2a  J4/J0  E_test  K  J4  t\n")
dt=args.dt; t0=time.time()
best_lam=lam0; best_step=start_step; n_best=n_np.copy()
E_prev=E0

print(f"\n  {'step':>7}  {'K/J4':>9}  {'K/J4/lam':>10}  "
      f"{'K/J2a':>7}  {'J4/J0':>7}  {'E_test':>10}  {'dt':>8}")
print(f"  {'-'*72}")

for step in range(start_step, args.steps+1):
    K_v,J2a_v,J4_v,E_v,Force=step_force(n_t,MU,LAM)
    lam=K_v/J4_v if J4_v>1e-12 else float('nan')

    if step%args.print_every==0 or step==start_step:
        el=time.time()-t0
        print(f"  {step:>7d}  {lam:>9.4f}  {lam/LAM:>10.5f}  "
              f"{K_v/J2a_v:>7.4f}  {J4_v/J4_init:>7.4f}  "
              f"{E_v:>10.3f}  {dt:>8.2e}  ({el:.0f}s)")
        log.write(f"{step}  {lam:.6f}  {K_v/J2a_v:.5f}  "
                  f"{J4_v/J4_init:.5f}  {E_v:.4f}  {K_v:.4f}  {J4_v:.6f}  {el:.1f}\n")
        log.flush()

    if abs(lam-LAM) < abs(best_lam-LAM):
        best_lam=lam; best_step=step; n_best=n_t.detach().cpu().numpy()

    with torch.no_grad():
        n_try=n_t+dt*Force
        n_try=n_try/n_try.norm(dim=-1,keepdim=True).clamp(1e-10)
    with torch.no_grad():
        K_try,_,J4_try=compute(n_try,MU)
    E_try=K_try.item()+LAM*J4_try.item()

    if E_try < E_v:
        n_t=n_try; dt=min(dt*1.008,args.max_dt)
    else:
        dt*=0.7
        if dt<1e-14: dt=args.dt

    if step%args.save_every==0:
        ck=f"{ckpt_base}_{step:08d}.npz"
        np.savez(ck,n=n_t.detach().cpu().numpy(),step=step,
                 J4_init=J4_init,K=K_v,J4=J4_v,lam=lam)
        print(f"  [ckpt {ck}]")

log.close()
with torch.no_grad():
    Kf,J2af,J4f=compute(n_t,MU)
Kf,J2af,J4f=Kf.item(),J2af.item(),J4f.item()
lamf=Kf/J4f if J4f>1e-12 else float('nan')

print(f"\n{'='*65}")
print(f"  RESULT  v5  lambda_test={LAM:.6f}")
print(f"{'='*65}")
print(f"  K/J4  = {lamf:.8f}  [converged; lambda_test = {LAM:.6f}]")
print(f"  Diff  = {abs(lamf-LAM)/LAM*100:.4f}%")
print(f"  K/J2a = {Kf/J2af:.8f}  [2phi = {2*phi:.6f}]")
print(f"  J4/J0 = {J4f/J4_init:.6f}  [>0 = topology preserved]")
print(f"  log_phi(K/J4) = {np.log(lamf)/np.log(phi):.6f}  [6.000 for phi^6]")
for k in range(1,20):
    if abs(np.log(lamf)/np.log(phi)-k)<0.07:
        print(f"  *** K/J4 ≈ phi^{k} = {phi**k:.6f} ({(lamf-phi**k)/phi**k*100:+.3f}%) ***")
print(f"\n  Best K/J4 = {best_lam:.6f} at step {best_step}")
np.savez(f"{ckpt_base}_FINAL.npz",n=n_t.detach().cpu().numpy(),
         step=step,J4_init=J4_init,K=Kf,J4=J4f,lam=lamf)
