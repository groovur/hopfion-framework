#!/usr/bin/env python3
"""
junction_flow.py
================
Gradient-flow minimisation of E_geom = K*J4 in the Y-junction domain:
a small box [-r_out, r_out]^3 centred on the Fermat-Steiner point (the
origin), with Dirichlet boundary conditions on the outer sphere set by
Construction C.

PHYSICAL SETUP:
  The trefoil T_{2,3} with R0=3, r0=0.874 has its Fermat-Steiner point
  exactly at the origin (verified: centroid of three crossing midpoints
  = (0,0,0) to machine precision).  The three Y-junction arms each have
  length exactly R0=3, meeting at 120° in the z=0 plane.

  The junction domain is the ball |x| < r_outer, embedded in a Cartesian
  box [-r_outer, r_outer]^3. Boundary points (|x| >= r_boundary) are
  set to Construction C and held FIXED throughout. Interior points
  (|x| < r_boundary) are free to evolve under gradient flow.

  r_boundary < r_outer < R0 ensures:
  - The arms are not inside the box (arm endpoints at distance R0=3)
  - The junction core IS inside the box
  - The tube profiles are represented on the boundary via Construction C

KEY DIFFERENCES FROM gradient_flow_constrained.py:
  1. DIRICHLET (not periodic) boundary conditions: finite differences at
     boundary points clamp to the fixed Construction C value.
  2. BOUNDARY MASK: points within r_boundary of the origin are interior
     (free); points outside are boundary (fixed). The gradient is zeroed
     at boundary points before each optimizer step.
  3. SMALL DOMAIN: r_outer ~ 1.5-2.0 << R0=3, so the box is much
     cheaper than the full trefoil solver.

WHAT E_junction MEASURES:
  E_junction = E_geom computed ONLY over interior points (|x| < r_junc),
  where r_junc < r_boundary is an inner measurement radius that excludes
  the boundary layer. This isolates the genuine junction energy from
  the arm-body contribution.

GENERATION DEPENDENCE:
  E_junction depends on the boundary condition, which comes from
  Construction C. At different WZW levels k=1,2,3 (quark generations),
  the effective tube profile changes through C* and the phi^6 coupling,
  making E_junction generation-dependent.
  Currently: C* = C*_3 = 2.5062 (physical saddle value, Paper XIV).
  Future: vary C* to probe generation dependence.
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--N',         type=int,   default=64,
                help='Grid points per side (box is [-r_out, r_out]^3)')
ap.add_argument('--r_outer',   type=float, default=2.0,
                help='Half-width of Cartesian box in natural units')
ap.add_argument('--r_boundary',type=float, default=1.6,
                help='Radius beyond which points are held fixed (BC layer)')
ap.add_argument('--r_junc',    type=float, default=0.8,
                help='Inner measurement radius for E_junction')
ap.add_argument('--R0',        type=float, default=3.0)
ap.add_argument('--r0',        type=float, default=0.874)
ap.add_argument('--C_star',    type=float, default=2.5062)
ap.add_argument('--n_steps',   type=int,   default=2000)
ap.add_argument('--lr',        type=float, default=3e-5)
ap.add_argument('--log_every', type=int,   default=50)
ap.add_argument('--outdir',    type=str,   default='junction_run')
ap.add_argument('--device',    type=str,   default='cpu')
args = ap.parse_args()

torch.manual_seed(0)
os.makedirs(args.outdir, exist_ok=True)

PHI  = (1+5**0.5)/2
MU   = 3.0 - PHI
N    = args.N
h    = 2*args.r_outer / N         # grid spacing
R0, r0, C_star = args.R0, args.r0, args.C_star
dev  = torch.device(args.device)

print(f"{'='*65}")
print(f"  JUNCTION GRADIENT FLOW  (Dirichlet BC from Construction C)")
print(f"  Box: [{-args.r_outer:.2f},{args.r_outer:.2f}]^3, N={N}, h={h:.4f}")
print(f"  r_boundary={args.r_boundary:.2f}, r_junc={args.r_junc:.2f}, r_outer={args.r_outer:.2f}")
print(f"  C*={C_star}, lr={args.lr}, n_steps={args.n_steps}")
print(f"{'='*65}")

if args.r_boundary >= args.r_outer:
    print("FATAL: r_boundary must be < r_outer"); sys.exit(1)
if args.r_boundary <= args.r_junc:
    print("FATAL: r_junc must be < r_boundary"); sys.exit(1)
if args.r_outer >= R0:
    print(f"WARNING: r_outer={args.r_outer} >= R0={R0}: arm endpoints inside box!")

# ── Grid ──────────────────────────────────────────────────────────────
cv   = -args.r_outer + h*(np.arange(N) + 0.5)
X,Y,Z = np.meshgrid(cv,cv,cv, indexing='ij')
pts_np = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1).astype(np.float32)
r_from_origin = np.linalg.norm(pts_np, axis=1).reshape(N,N,N)
r_t = torch.tensor(r_from_origin, dtype=torch.float32, device=dev)

# Masks (True = interior/free, False = boundary/fixed)
interior_mask = torch.tensor(r_from_origin < args.r_boundary,
                              dtype=torch.bool, device=dev)  # (N,N,N)
junc_mask     = torch.tensor(r_from_origin < args.r_junc,
                              dtype=torch.bool, device=dev)

n_interior = interior_mask.sum().item()
n_boundary = (~interior_mask).sum().item()
n_junc     = junc_mask.sum().item()
print(f"  Interior points (free):   {n_interior:>8d}  ({100*n_interior/N**3:.1f}%)")
print(f"  Boundary points (fixed):  {n_boundary:>8d}  ({100*n_boundary/N**3:.1f}%)")
print(f"  Junction measurement pts: {n_junc:>8d}  (|x|<{args.r_junc})")

# ── Energy functional (Dirichlet: uses clamp-to-boundary differences) ──
# For interior points: standard central differences, possibly hitting
# boundary points whose values are fixed.
# For boundary points: we never update them, so their gradient is zeroed.
def E_geom_dirichlet(n):
    """E_geom = K*J4 using standard central differences.
    Boundary points are fixed (their gradient will be zeroed later).
    n: (N,N,N,3) unit tensor."""
    nx,ny,nz = n[...,0],n[...,1],n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u,a):
        # Central difference: at boundaries this uses the fixed BC values,
        # which is correct for Dirichlet -- the gradient of interior points
        # near the boundary automatically incorporates the BC.
        return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2   = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a  = (s4*g2).sum()*h**3
    J2iso = g2.sum()*h**3
    K    = J2a + MU*J2iso
    Fxy  = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz  = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz  = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho  = Fxy**2 + Fxz**2 + Fyz**2
    J4   = rho.sum()*h**3
    return K*J4, K.item(), J4.item(), rho

def diagnostics(n):
    with torch.no_grad():
        E,K,J4,rho = E_geom_dirichlet(n)
        # E_junction: K * integral of rho over junction region
        E_junc_approx = float(rho[junc_mask].sum() * h**3 * K)
    return E.item(), K, J4, E_junc_approx

# ── Build Construction C boundary condition ───────────────────────────
print("\nBuilding Construction C initial/boundary field...")
t0b = time.time()
NT_frame = 20000
t_frame, _, N1f, N2f, H = build_compensated_frame_arclength(NT=NT_frame)
print(f"  Bishop frame holonomy: {np.degrees(H):.4f} deg")

NT = 4000
t_arr = np.linspace(0,2*np.pi,NT,endpoint=False)
Gx=(R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy=(R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz=r0*np.sin(3*t_arr)
Gamma_pts=np.stack([Gx,Gy,Gz],axis=1)
arc_starts=[0,2*np.pi/3,4*np.pi/3]
lobe_indices=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees=[KDTree(Gamma_pts[li]) for li in lobe_indices]
lobe_t_arrays=[t_arr[li] for li in lobe_indices]

def nearest_two(qpts):
    d_per,t_per=[],[]
    for tr,tl in zip(lobe_trees,lobe_t_arrays):
        d,idx=tr.query(qpts,workers=1); d_per.append(d); t_per.append(tl[idx])
    d_s=np.stack(d_per,axis=1); t_s=np.stack(t_per,axis=1)
    o=np.argsort(d_s,axis=1)
    return (np.take_along_axis(t_s,o,axis=1)[:,0],
            np.take_along_axis(d_s,o,axis=1)[:,0],
            np.take_along_axis(t_s,o,axis=1)[:,1],
            np.take_along_axis(d_s,o,axis=1)[:,1])

def frame_at(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NT_frame
    return N1f[idx],N2f[idx]
def curve_at(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                     (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)],axis=-1)

t1,d1,t2,d2 = nearest_two(pts_np)
chi1=np.arctan2(np.sum((pts_np-curve_at(t1))*frame_at(t1)[1],axis=1),
                np.sum((pts_np-curve_at(t1))*frame_at(t1)[0],axis=1))
chi2=np.arctan2(np.sum((pts_np-curve_at(t2))*frame_at(t2)[1],axis=1),
                np.sum((pts_np-curve_at(t2))*frame_at(t2)[0],axis=1))
Phi1=chi1+3*t1; Phi2=chi2+3*t2
rho1=np.clip(d1,1e-6,None); rho2=np.clip(d2,1e-6,None)
f1=2*np.arctan(np.maximum(rho1*C_star,1e-9)**(-C_star))
f2=2*np.arctan(np.maximum(rho2*C_star,1e-9)**(-C_star))
w1=1/rho1**2; w2=1/rho2**2
z1_np=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
z2_np=w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
mag=np.sqrt(np.abs(z1_np)**2+np.abs(z2_np)**2)
z1_np/=mag; z2_np/=mag
nx0=2*np.real(np.conj(z1_np)*z2_np)
ny0=2*np.imag(np.conj(z1_np)*z2_np)
nz0=np.abs(z1_np)**2-np.abs(z2_np)**2
n0_np=np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
n0_np/=np.linalg.norm(n0_np,axis=-1,keepdims=True).clip(1e-10)
print(f"  Built in {time.time()-t0b:.1f}s")

# ── Initial diagnostics ───────────────────────────────────────────────
n_t = torch.tensor(n0_np, dtype=torch.float32, device=dev)
E0, K0, J40, Ejunc0 = diagnostics(n_t)
print(f"\nInitial field (Construction C throughout):")
print(f"  E_total = {E0:.4e}  K = {K0:.2f}  J4 = {J40:.2f}")
print(f"  E_junction (|x|<{args.r_junc}) ~ {Ejunc0:.4e}")

# ── Setup optimiser ───────────────────────────────────────────────────
n_param = n_t.clone().requires_grad_(True)
opt = torch.optim.Adam([n_param], lr=args.lr)

# Store boundary values -- these never change
boundary_vals = n_t.detach().clone()  # (N,N,N,3) - full field
# boundary_vals[interior_mask] will be overwritten each step;
# boundary_vals[~interior_mask] are the fixed Dirichlet values

def project_tangent(n, grad):
    return grad - (grad*n).sum(-1,keepdim=True)*n

history = []
log_path = os.path.join(args.outdir,'log.json')
t0 = time.time()

print(f"\nRunning {args.n_steps} steps (lr={args.lr})...")
print(f"{'step':>6}  {'E_total':>12}  {'K':>8}  {'J4':>8}  "
      f"{'E_junc':>12}  {'|grad|':>10}  {'BC_err':>8}")

for step in range(args.n_steps):
    opt.zero_grad()
    E, _, _, _ = E_geom_dirichlet(n_param)
    E.backward()

    with torch.no_grad():
        # Project gradient to tangent space of S^2
        n_param.grad.data.copy_(project_tangent(n_param.detach(), n_param.grad))
        # ZERO GRADIENT AT BOUNDARY POINTS -- they must not move
        n_param.grad.data[~interior_mask] = 0.0

    grad_norm = n_param.grad[interior_mask].norm().item()
    opt.step()

    with torch.no_grad():
        # Renormalise
        n_param.data.copy_(
            n_param / n_param.norm(dim=-1,keepdim=True).clamp(1e-10))
        # ENFORCE DIRICHLET BC: restore boundary points exactly
        n_param.data[~interior_mask] = boundary_vals[~interior_mask]

    if (step+1) % args.log_every == 0 or step == 0:
        E_val, K_val, J4_val, Ejunc_val = diagnostics(n_param)
        # Boundary condition error: max deviation of boundary points from BC
        bc_err = (n_param.detach()[~interior_mask] -
                  boundary_vals[~interior_mask]).norm(dim=-1).max().item()

        print(f"{step+1:>6}  {E_val:>12.4e}  {K_val:>8.1f}  {J4_val:>8.2f}  "
              f"{Ejunc_val:>12.4e}  {grad_norm:>10.4e}  {bc_err:>8.2e}")

        history.append(dict(step=step+1, E=E_val, K=K_val, J4=J4_val,
                            E_junc=Ejunc_val, grad_norm=grad_norm, bc_err=bc_err))

        # Convergence: gradient norm near zero
        if grad_norm < 1e-2 * (history[0]['grad_norm'] if len(history)>0 else 1):
            print(f"\nCONVERGED: |grad| < 1% of initial at step {step+1}")
            break

        # Abort if BC is violated (should never happen)
        if bc_err > 1e-5:
            print(f"\nWARNING: BC violation {bc_err:.2e} at step {step+1}")

print(f"\nTotal wall time: {time.time()-t0:.1f}s")
with open(log_path,'w') as f: json.dump(history,f,indent=2)

E_f, K_f, J4_f, Ejunc_f = diagnostics(n_param)
print(f"\nFinal:")
print(f"  E_total  = {E_f:.4e}  (initial: {E0:.4e}, reduction: {E0/E_f:.6f}x)")
print(f"  E_junc   = {Ejunc_f:.4e}  (initial: {Ejunc0:.4e}, reduction: {Ejunc0/max(Ejunc_f,1e-30):.4f}x)")
print(f"  E_junc / E_total = {Ejunc_f/max(E_f,1e-30):.4e}")
print(f"  E_junc (solver units) vs E_profile* (thm:norm3):")
print(f"    E_junc = {Ejunc_f:.6f}")
print(f"    Note: unit conversion to L_cond units needed for physical comparison")

np.save(os.path.join(args.outdir,'n_junction.npy'),
        n_param.detach().cpu().numpy())
print(f"Field saved: {args.outdir}/n_junction.npy")
print(f"Log:  {log_path}")
