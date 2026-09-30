#!/usr/bin/env python3.11
"""
tube_relax.py -- PHASE 4 of the tube-adapted Q_H=3 saddle build
(notes/qh3_saddle_infrastructure_scope.md). Relaxes the director field
n(s,rho,psi) on the trefoil tube grid via autograd on the curvilinear
Faddeev energy (validated in tube_curvilinear.py Phases 1-2 and tube_trefoil.py
Phase 3), with the density-feedback isotropic softening and the Paper-I
saddle-snapshot. Mirrors the proven qh3_feedback_saddle.py relaxation
(two-phase lr, tangent-projected Adam, angle-clamp topology protection,
best-self-consistency-score snapshot, plateau/ringing tail verdict) but on the
DOF-efficient tube grid (~1M DOF) instead of a uniform Cartesian grid (~30M).

GEOMETRY IS PRECOMPUTED ONCE (§8): frame, h_s, partition weight W, lab positions,
and the construction-C initial field. The relaxation variable is n on the grid;
each step only updates n and re-evaluates the curvilinear energy -> cheap. The
OUTER rho shell and the crossing (overlap, W<1) nodes are PINNED to the
construction (anchors the winding; crossings carry <1% of the energy, Phase 3).

Self-consistency targets (Paper I / Q_H=2): V = J2iso_fb/J2a -> phi ;
sopt = (phi^6 J4 / K_fb)^1/2 -> 1 ; score = |V-phi| + |sopt-1|.

Smoke:  python3.11 tube_relax.py --Ns 120 --Nr 16 --Np 32 --n_steps 40 --beta 0.35
Real:   python3.11 tube_relax.py --Ns 420 --Nr 40 --Np 80 --n_steps 3000 \
                --beta_list 0.0,0.2,0.35,0.452 --outdir tube_saddle_scan
"""
import numpy as np, torch, time, json, os, argparse
from tube_curvilinear import MU, PHI
from tube_trefoil import BlendField, trefoil_frame, partition_weight, R0, r0, C_STAR

LAM = PHI**6

ap = argparse.ArgumentParser()
ap.add_argument('--Ns', type=int, default=300)
ap.add_argument('--Nr', type=int, default=30)
ap.add_argument('--Np', type=int, default=64)
ap.add_argument('--rho_max', type=float, default=1.1)
ap.add_argument('--Cstar', type=float, default=C_STAR)
ap.add_argument('--beta', type=float, default=0.35)
ap.add_argument('--beta_list', type=str, default='')
ap.add_argument('--n_steps', type=int, default=2000)
ap.add_argument('--lr1', type=float, default=3e-4)
ap.add_argument('--lr2', type=float, default=1e-5)
ap.add_argument('--delta_max_deg', type=float, default=9.0)
ap.add_argument('--K_rise_eps', type=float, default=2.0)
ap.add_argument('--log_every', type=int, default=25)
ap.add_argument('--outdir', type=str, default='tube_saddle')
ap.add_argument('--device', type=str, default='cpu')
args = ap.parse_args()
DELTA_MAX = np.radians(args.delta_max_deg)
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)


# ── Precompute geometry + construction-C initial field ONCE ─────────────────
def build_geometry():
    print(f"Precomputing tube geometry (Ns,Nr,Np={args.Ns},{args.Nr},{args.Np}, "
          f"rho_max={args.rho_max}, C*={args.Cstar})...")
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
    dV = hs*RHO*ds*drho*dpsi*W                              # fixed volume measure
    rlab = np.linalg.norm(Xg, axis=-1)                      # |X| for r_bar
    # pinned = outer rho shell OR crossing-overlap node (anchors topology)
    pinned = np.zeros((args.Ns, args.Nr, args.Np), bool)
    pinned[:, -1, :] = True
    pinned |= (W < 0.999)
    print(f"  done ({time.time()-t0:.1f}s). pinned {100*pinned.mean():.1f}% of nodes; "
          f"h_s in [{hs.min():.3f},{hs.max():.3f}]")
    to = lambda a: torch.tensor(a, dtype=torch.float32, device=dev)
    return dict(n0=to(n0), dV=to(dV), hs=to(hs), RHO=to(RHO), rlab=to(rlab),
                pinned=torch.tensor(pinned, device=dev), ds=ds, drho=drho, dpsi=dpsi)


G = build_geometry()


def energy(n, beta):
    """Curvilinear Faddeev energy (torch, grad-enabled). Mirrors tube_curvilinear."""
    hs, RHO, dV = G['hs'], G['RHO'], G['dV']
    Dn_s   = (torch.roll(n, -1, 0) - torch.roll(n, 1, 0))/(2*G['ds'])   / hs[..., None]
    Dn_psi = (torch.roll(n, -1, 2) - torch.roll(n, 1, 2))/(2*G['dpsi']) / RHO[..., None]
    interior = (n[:, 2:] - n[:, :-2])/(2*G['drho'])
    low  = (n[:, 1:2] - n[:, 0:1])/G['drho']
    high = (n[:, -1:] - n[:, -2:-1])/G['drho']
    Dn_rho = torch.cat([low, interior, high], dim=1)
    g2 = (Dn_s**2).sum(-1) + (Dn_rho**2).sum(-1) + (Dn_psi**2).sum(-1)
    s4 = (1 - n[..., 2]**2).clamp(0, 1)**2
    J2a      = (s4*g2*dV).sum()
    J2iso_fb = (g2/(1.0 + beta*g2)*dV).sum()
    K_fb     = J2a + MU*J2iso_fb
    Fsr = (n*torch.cross(Dn_s,   Dn_rho, dim=-1)).sum(-1)
    Fsp = (n*torch.cross(Dn_s,   Dn_psi, dim=-1)).sum(-1)
    Frp = (n*torch.cross(Dn_rho, Dn_psi, dim=-1)).sum(-1)
    rho_J4 = Fsr**2 + Fsp**2 + Frp**2
    J4 = (rho_J4*dV).sum()
    return K_fb*J4, K_fb, J4, rho_J4, J2a, J2iso_fb


def diagnostics(n, beta):
    with torch.no_grad():
        E, K_fb, J4, rho, J2a, J2iso_fb = energy(n, beta)
        r_bar = ((rho*G['dV']*G['rlab']).sum()/(rho*G['dV']).sum().clamp(1e-12)).item()
        K_fb, J4, J2a, J2iso_fb = K_fb.item(), J4.item(), J2a.item(), J2iso_fb.item()
        V    = J2iso_fb/J2a if J2a > 1e-12 else 0.0
        sopt = (LAM*J4/K_fb)**0.5 if K_fb > 1e-12 else 0.0
    return dict(E=E.item(), K_fb=K_fb, J4=J4, J2a=J2a, J2iso_fb=J2iso_fb, r_bar=r_bar,
                V=V, sopt=sopt, sopt6=sopt**6, score=abs(V-PHI)+abs(sopt-1.0),
                J4_over_J2a=(J4/J2a if J2a>1e-12 else 0.0),
                Kfb_over_J4=(K_fb/J4 if J4>1e-12 else 0.0))


def project_to_tangent(n, grad):
    return grad - (grad*n).sum(-1, keepdim=True)*n


def apply_angle_clamp(n_old, n_new, dmax):
    cos_a = (n_old*n_new).sum(-1, keepdim=True).clamp(-1+1e-6, 1-1e-6)
    angle = torch.acos(cos_a)
    too_far = (angle > dmax).squeeze(-1)
    if not too_far.any():
        return n_new
    sin_a = torch.sin(angle).clamp(1e-8); t = dmax/angle.clamp(1e-8)
    n_slerp = (torch.sin((1-t)*angle)/sin_a*n_old + torch.sin(t*angle)/sin_a*n_new)
    n_slerp = n_slerp/n_slerp.norm(dim=-1, keepdim=True).clamp(1e-10)
    out = n_new.clone(); out[too_far] = n_slerp[too_far]
    return out/out.norm(dim=-1, keepdim=True).clamp(1e-10)


def run_beta(beta):
    n0 = G['n0']; pinned = G['pinned']
    d0 = diagnostics(n0, beta)
    print(f"\n{'-'*72}\n  beta={beta:.5f}\n{'-'*72}")
    print(f"  init: E={d0['E']:.3e} K_fb={d0['K_fb']:.2f} J4={d0['J4']:.2f} "
          f"V={d0['V']:.4f}(phi={PHI:.4f}) sopt={d0['sopt']:.4f} score={d0['score']:.4f} "
          f"r_bar={d0['r_bar']:.3f}")
    n_param = n0.clone().requires_grad_(True)
    opt = torch.optim.Adam([n_param], lr=args.lr1)
    phase = 1; K_min = d0['K_fb']
    best_score = d0['score']; n_best = n0.clone().cpu().numpy(); best = dict(d0, step=0)
    history = []; t0 = time.time()
    print(f"  {'step':>6} {'ph':>2} {'E':>11} {'K_fb':>8} {'J4':>8} {'V':>7} {'sopt':>7} "
          f"{'score':>7} {'r_bar':>6} {'|grad|':>9} {'clmp%':>6}")
    for step in range(args.n_steps):
        n_before = n_param.detach().clone()
        opt.zero_grad()
        E, _, _, _, _, _ = energy(n_param, beta)
        E.backward()
        with torch.no_grad():
            n_param.grad.copy_(project_to_tangent(n_param.detach(), n_param.grad))
            n_param.grad[pinned] = 0.0                       # freeze pinned nodes
        grad_norm = n_param.grad.norm().item()
        opt.step()
        with torch.no_grad():
            n_param.data.copy_(n_param/n_param.norm(dim=-1, keepdim=True).clamp(1e-10))
            clamped = 0.0
            if phase == 2:
                n_cl = apply_angle_clamp(n_before, n_param.detach(), DELTA_MAX)
                clamped = ((n_cl-n_param.detach()).norm(dim=-1) > 1e-6).float().mean().item()
                n_param.data.copy_(n_cl)
            n_param.data[pinned] = n0[pinned]                # hard-pin

        cs = step+1
        if cs % args.log_every == 0 or step == 0:
            d = diagnostics(n_param, beta); d['grad_norm'] = grad_norm; d['step'] = cs
            print(f"  {cs:>6} {phase:>2} {d['E']:>11.3e} {d['K_fb']:>8.1f} {d['J4']:>8.2f} "
                  f"{d['V']:>7.4f} {d['sopt']:>7.4f} {d['score']:>7.4f} {d['r_bar']:>6.3f} "
                  f"{grad_norm:>9.3e} {100*clamped:>5.1f}%")
            history.append(d)
            if d['score'] < best_score:
                best_score = d['score']; n_best = n_param.detach().cpu().numpy().copy(); best = dict(d)
            if phase == 1:
                if d['K_fb'] < K_min: K_min = d['K_fb']
                elif d['K_fb'] > K_min + args.K_rise_eps:
                    print(f"    *** PHASE 2 at step {cs}: K_fb={d['K_fb']:.2f}>{K_min:.2f}+eps "
                          f"(lr={args.lr2}, clamp={args.delta_max_deg}deg) ***")
                    phase = 2; opt = torch.optim.Adam([n_param], lr=args.lr2)
            if d['J4_over_J2a'] < 0.01 and step > 50:
                print(f"    HALT: J4/J2a collapsed at step {cs} (topology lost)"); break
    wall = time.time()-t0
    print(f"  best score={best_score:.5f} at step {best['step']}: V={best['V']:.5f} "
          f"sopt={best['sopt']:.5f} J4/J2a={best['J4_over_J2a']:.5f} "
          f"K_fb/J4={best['Kfb_over_J4']:.4f} r_bar={best['r_bar']:.3f}")
    print(f"  wall {wall:.1f}s ({wall/max(cs,1)*1000:.1f} ms/step)")

    tail = history[max(0, 3*len(history)//4):]
    if len(tail) >= 4:
        def st(k):
            a = np.array([h[k] for h in tail], float)
            dif = np.diff(a)
            return (a[-1]-a[0])/(abs(a[0])+1e-12), (a.max()-a.min())/(abs(a.mean())+1e-12), \
                   int(np.sum(dif[:-1]*dif[1:] < 0))
        j4 = st('J4'); sc = st('score')
        if j4[2] >= 3 and j4[1] > 0.05:      v = "RINGING (J4 oscillating)"
        elif abs(j4[0]) < 0.03 and sc[1] < 0.05: v = "PLATEAU -- metastable saddle candidate"
        elif j4[0] < -0.10:                  v = "SLOW COLLAPSE (J4 bleeding down)"
        else:                                v = "AMBIGUOUS -- extend steps"
        print(f"  tail (steps {tail[0]['step']}-{tail[-1]['step']}): J4 drift {j4[0]:+.2%} "
              f"pk-pk {j4[1]:.2%} flips {j4[2]}  VERDICT: {v}")
        best['tail_verdict'] = v

    tag = f"beta{beta:.3f}"
    np.save(os.path.join(args.outdir, f'n_best_{tag}.npy'), n_best)
    np.save(os.path.join(args.outdir, f'n_final_{tag}.npy'), n_param.detach().cpu().numpy())
    with open(os.path.join(args.outdir, f'history_{tag}.json'), 'w') as f:
        json.dump(history, f, indent=2)
    best['beta'] = beta; best['wall_s'] = wall
    return best


if __name__ == "__main__":
    betas = [float(x) for x in args.beta_list.split(',')] if args.beta_list else [args.beta]
    summary = [run_beta(b) for b in betas]
    with open(os.path.join(args.outdir, 'scan_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary -> {os.path.join(args.outdir, 'scan_summary.json')}")
    print(f"  {'beta':>7} {'score':>8} {'V':>8} {'sopt':>8} {'J4/J2a':>8} {'verdict'}")
    for s in summary:
        print(f"  {s['beta']:>7.3f} {s['score']:>8.4f} {s['V']:>8.4f} {s['sopt']:>8.4f} "
              f"{s['J4_over_J2a']:>8.4f} {s.get('tail_verdict','')}")
