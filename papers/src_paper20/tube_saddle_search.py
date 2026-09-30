#!/usr/bin/env python3.11
"""
tube_saddle_search.py -- PHASE 5 / step (b): does the isolated Q_H=3 self-consistent
configuration actually EXIST as a critical point, and is it STABLE or a SADDLE?

Pure E-descent cannot answer this: E=K*J4 is scale-invariant so descent collapses,
and descent structurally slides off any saddle. Instead:
  1. SADDLE-SEARCH: minimise G(n) = 1/2 |P grad E|^2 (tangent-projected gradient
     norm squared). Critical points of E are the zeros of G; minimising G finds
     them regardless of stability. Monitor J4 to tell a NONTRIVIAL critical point
     (J4 stays large) from drift to the trivial vacuum (J4 -> 0).
  2. STABILITY: at the located critical point, the lowest tangent Hessian
     eigenvalue lambda_min (via Hessian-vector products + shifted power iteration)
     classifies it: lambda_min < 0 (clearly) => UNSTABLE SADDLE (quark-like, the
     object is not held open in isolation); lambda_min ~ 0 with the rest > 0 =>
     stable up to the marginal SCALE mode (which self-consistency fixes) = a
     genuine object (lepton-like).

Geometry (frame, h_s, partition W, pinned nodes, construction IC) reused from the
validated tube pipeline (tube_trefoil). The energy mirrors tube_relax / tube_curvilinear.

Smoke:  python3.11 tube_saddle_search.py --Ns 120 --Nr 16 --Np 32 --n_steps 60 --beta 0.24
"""
import numpy as np, torch, time, json, os, argparse
from tube_curvilinear import MU, PHI
from tube_trefoil import BlendField, trefoil_frame, partition_weight

LAM = PHI**6
ap = argparse.ArgumentParser()
ap.add_argument('--Ns', type=int, default=300)
ap.add_argument('--Nr', type=int, default=30)
ap.add_argument('--Np', type=int, default=64)
ap.add_argument('--rho_max', type=float, default=1.1)
ap.add_argument('--Cstar', type=float, default=2.5062)
ap.add_argument('--beta', type=float, default=0.24)      # trefoil self-consistent beta_sc
ap.add_argument('--n_steps', type=int, default=1500)
ap.add_argument('--lr', type=float, default=2e-3)
ap.add_argument('--log_every', type=int, default=25)
ap.add_argument('--hess_iters', type=int, default=60)
ap.add_argument('--outdir', type=str, default='tube_saddle_search')
ap.add_argument('--device', type=str, default='cpu')
args = ap.parse_args()
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)


def build_geometry():
    print(f"Precomputing geometry (Ns,Nr,Np={args.Ns},{args.Nr},{args.Np}, C*={args.Cstar})...")
    t0 = time.time()
    bf = BlendField(Cst=args.Cstar)
    s_grid, Gamma, T, N1, N2, k1, k2, L = trefoil_frame(args.Ns)
    rho_grid = (np.arange(args.Nr)+0.5)*(args.rho_max/args.Nr)
    psi_grid = np.arange(args.Np)*(2*np.pi/args.Np)
    ds = s_grid[1]-s_grid[0]; drho = rho_grid[1]-rho_grid[0]; dpsi = psi_grid[1]-psi_grid[0]
    S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
    cosp, sinp = np.cos(PSI), np.sin(PSI)
    er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
    Xg = Gamma[:, None, None, :] + RHO[..., None]*er
    n0 = bf.n_at(Xg); n0 /= np.linalg.norm(n0, axis=-1, keepdims=True).clip(1e-12)
    W = partition_weight(bf, Xg, S, RHO, L, args.rho_max)
    hs = 1.0 - RHO*(k1[:, None, None]*cosp + k2[:, None, None]*sinp)
    dV = hs*RHO*ds*drho*dpsi*W
    pinned = np.zeros((args.Ns, args.Nr, args.Np), bool)
    pinned[:, -1, :] = True; pinned |= (W < 0.999)
    print(f"  done ({time.time()-t0:.1f}s). pinned {100*pinned.mean():.1f}%")
    to = lambda a: torch.tensor(a, dtype=torch.float64, device=dev)   # float64 for Hessian accuracy
    return dict(n0=to(n0), dV=to(dV), hs=to(hs), RHO=to(RHO),
                pinned=torch.tensor(pinned, device=dev), ds=ds, drho=drho, dpsi=dpsi)


G = build_geometry()


def energy(n, beta):
    hs, RHO, dV = G['hs'], G['RHO'], G['dV']
    Dn_s = (torch.roll(n,-1,0)-torch.roll(n,1,0))/(2*G['ds'])/hs[...,None]
    Dn_p = (torch.roll(n,-1,2)-torch.roll(n,1,2))/(2*G['dpsi'])/RHO[...,None]
    inter=(n[:,2:]-n[:,:-2])/(2*G['drho']); lo=(n[:,1:2]-n[:,0:1])/G['drho']; hi=(n[:,-1:]-n[:,-2:-1])/G['drho']
    Dn_r = torch.cat([lo,inter,hi],1)
    g2 = (Dn_s**2).sum(-1)+(Dn_r**2).sum(-1)+(Dn_p**2).sum(-1)
    s4 = (1-n[...,2]**2).clamp(0,1)**2
    J2a=(s4*g2*dV).sum(); J2iso=(g2/(1+beta*g2)*dV).sum(); K=J2a+MU*J2iso
    Fsr=(n*torch.cross(Dn_s,Dn_r,dim=-1)).sum(-1)
    Fsp=(n*torch.cross(Dn_s,Dn_p,dim=-1)).sum(-1)
    Frp=(n*torch.cross(Dn_r,Dn_p,dim=-1)).sum(-1)
    J4=((Fsr**2+Fsp**2+Frp**2)*dV).sum()
    return K*J4, K, J4, J2a, J2iso


def tangent(n, v):
    return v - (v*n).sum(-1, keepdim=True)*n


def diag(n, beta):
    with torch.no_grad():
        E,K,J4,J2a,J2iso = energy(n,beta)
        V = (J2iso/J2a).item(); sopt=((LAM*J4/K)**0.5).item()
    return E.item(), K.item(), J4.item(), V, sopt


def grad_E_tan(n, beta, create_graph=False):
    E = energy(n, beta)[0]
    g = torch.autograd.grad(E, n, create_graph=create_graph)[0]
    return tangent(n, g)


def run():
    beta = args.beta; n0 = G['n0']; pinned = G['pinned']
    E0,K0,J40,V0,so0 = diag(n0, beta)
    gt0 = grad_E_tan(n0.clone().requires_grad_(True), beta).norm().item()
    print(f"\n{'='*70}\n saddle-search  beta={beta}  (target: |grad E|->0 with J4 large)\n{'='*70}")
    print(f"  init: E={E0:.3e} J4={J40:.2f} V={V0:.4f}(phi={PHI:.4f}) sopt={so0:.4f} |gradE|={gt0:.3e}")
    n = n0.clone().requires_grad_(True)
    opt = torch.optim.Adam([n], lr=args.lr)
    print(f"  {'step':>6} {'G=½|gradE|²':>13} {'|gradE|':>10} {'J4':>9} {'V':>7} {'sopt':>7}")
    hist=[]
    for step in range(args.n_steps):
        opt.zero_grad()
        gtan = grad_E_tan(n, beta, create_graph=True)
        Gval = 0.5*(gtan[~pinned]**2).sum()   # free nodes only; pinned are BCs, not DOF
        Gval.backward()
        with torch.no_grad():
            n.grad.copy_(tangent(n.detach(), n.grad)); n.grad[pinned] = 0.0
        opt.step()
        with torch.no_grad():
            n.data.copy_(n/n.norm(dim=-1,keepdim=True).clamp(1e-12))
            n.data[pinned] = n0[pinned]
        if (step+1) % args.log_every == 0 or step==0:
            gnorm = (2*Gval.item())**0.5
            E,K,J4,V,so = diag(n, beta)
            print(f"  {step+1:>6} {Gval.item():>13.4e} {gnorm:>10.3e} {J4:>9.2f} {V:>7.4f} {so:>7.4f}")
            hist.append(dict(step=step+1, G=Gval.item(), gradE=gnorm, J4=J4, V=V, sopt=so, E=E))
    # ── stability: lowest tangent-Hessian eigenvalue via shifted power iteration ──
    print("\n  Hessian stability at the located point (shifted power iteration):")
    n_star = n.detach().clone().requires_grad_(True)
    def Hv(v):
        g = grad_E_tan(n_star, beta, create_graph=True)
        hv = torch.autograd.grad((g*v).sum(), n_star, retain_graph=False)[0]
        return tangent(n_star.detach(), hv)
    torch.manual_seed(0)
    v = tangent(n_star.detach(), torch.randn_like(n_star)); v[pinned]=0; v/=v.norm()
    # 1) largest |eigenvalue| (for shift)
    lam_top=0
    for _ in range(30):
        w = Hv(v); w[pinned]=0; lam_top=(v*w).sum().item(); nw=w.norm()
        if nw<1e-30: break
        v = w/nw
    shift = abs(lam_top)*1.2 + 1e-6
    # 2) most-negative eigenvalue: dominant of (shift*I - H)
    v = tangent(n_star.detach(), torch.randn_like(n_star)); v[pinned]=0; v/=v.norm()
    lam_min=0
    for _ in range(args.hess_iters):
        w = shift*v - Hv(v); w[pinned]=0; nw=w.norm()
        if nw<1e-30: break
        v = w/nw
        lam_min = (v*Hv(v)).sum().item()
    print(f"    lambda_top≈{lam_top:.4e}  lambda_min≈{lam_min:.4e}")
    verdict = ("UNSTABLE SADDLE (lambda_min<0: an unstable shape mode -> quark not held open in isolation)"
               if lam_min < -1e-6*shift else
               "STABLE up to marginal mode (lambda_min>=0: object held open -> lepton-like)")
    print(f"    VERDICT: {verdict}")
    np.save(os.path.join(args.outdir, f'n_crit_beta{beta:.3f}.npy'), n.detach().cpu().numpy())
    json.dump(dict(history=hist, lam_top=lam_top, lam_min=lam_min, verdict=verdict, beta=beta),
              open(os.path.join(args.outdir,f'search_beta{beta:.3f}.json'),'w'), indent=2)
    print(f"  saved -> {args.outdir}/")


if __name__ == "__main__":
    run()
