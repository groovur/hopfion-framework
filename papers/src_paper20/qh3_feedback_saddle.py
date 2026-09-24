#!/usr/bin/env python3
"""
qh3_feedback_saddle.py  —  density-feedback saddle-snapshot solver for Q_H=3
===========================================================================
WHY THIS SCRIPT EXISTS
----------------------
The existing 3D Q_H=3 trefoil solvers (gradient_flow_solver.py,
gradient_flow_constrained.py, qh3_trefoil_solver_3d_v11.py) all minimise the
*bare* geometric energy  E_geom = K * J4  with  K = J2a + mu*J2iso  and NO
density feedback (they use the plain isotropic term J2iso, i.e. beta=0), by
pure descent / annealing. They are missing BOTH ingredients that Paper I's
axisymmetric Q_H=2 solver (src_paper1/fn_hopfion_solver.py) uses:

  (1) the DENSITY-FEEDBACK term:  J2iso_fb = integral kern/(1+beta*kern) dV
      (the isotropic stiffness is softened where the local energy density
      kern is high; beta* ~ 0.452 in Paper I units, Paper VIII/IX), and

  (2) the SADDLE-SNAPSHOT location method. Paper I is explicit
      (fn_hopfion_solver.py l.37-41): E_fb has no Derrick minimum in R^3 ---
      gradient descent ALWAYS collapses toward the trivial field. THE
      HOPFION IS A SADDLE-POINT of E_fb, not a minimum. It is located NOT by
      waiting for descent to converge, but by tracking the self-consistency
      score  score = |V - phi| + |sopt - 1|  throughout the descent and
      reverting to the best-scoring intermediate snapshot.

So "the isolated Q_H=3 saddle is transient" is not a tooling failure to fix
by a better minimiser --- it is the established nature of the object, TRUE
ALREADY FOR Q_H=2. The correct instrument is Paper I's: descend on E_fb WITH
the feedback term, score the self-consistency conditions, and snapshot the
best score. This script ports that instrument from the axisymmetric profile
formulation (Paper I) to the full 3D director field n(x) of the trefoil.

WHAT IS PORTED, AND THE 3D GENERALISATION
-----------------------------------------
  Paper I (axisymmetric)          ->   here (3D director n)
  kern = fr^2+fz^2 + sin^2(f)*A   ->   kern = g2 = |grad n|^2   (isotropic
                                        gradient energy density)
  J2a  = int sin^4(f)*kern dV     ->   J2a  = int s4 * g2 dV,  s4=(1-nz^2)^2
  J2iso_fb = int kern/(1+b*kern)  ->   J2iso_fb = int g2/(1+beta*g2) dV
  K_fb = J2a + mu*J2iso_fb        ->   same,  mu = 3 - phi
  J4   = int (Fxy^2+Fxz^2+Fyz^2)  ->   same (Hopf/Skyrme term)
  E_fb = K_fb * J4                ->   same
  V    = J2iso_fb/J2a  (-> phi)   ->   same
  sopt = (phi^6 * J4/K_fb)^(1/2)  ->   same   (Derrick target sopt=1;
                                        WZW universality reports sopt^6)
  score= |V-phi| + |sopt-1|       ->   same; revert to best-score snapshot

Autograd computes the descent force from E_fb directly (the feedback term
1/(1+beta*g2) is smooth in n), so no hand-coded Euler-Lagrange force is
needed --- unlike Paper I's numpy implementation.

TOPOLOGY PROTECTION: reused verbatim from gradient_flow_constrained.py ---
per-point post-step angle clamp (slerp any grid point rotating more than
delta_max back onto the feasible set) so the descent cannot jump the
topological wall in a single Adam step. Plus the dilution safeguards.

TWO OPEN, HONESTLY-FLAGGED CALIBRATION POINTS (the physics the local run
must settle, NOT assumed here):
  * beta transfer: Paper I's kern (with its geometric A(r,z) term and its
    discretisation) and this 3D g2 have DIFFERENT magnitudes, so beta*=0.452
    is NOT guaranteed to be the physical coupling here. The physical beta is
    the one at which the self-consistency score is minimised --- hence
    --beta_list runs a warm-started scan (as Paper I does) and reports
    score(beta); pick the minimiser, do not assume 0.452.
  * self-consistency TARGET: V=phi and sopt=1 are the Q_H=2 conditions.
    Whether the Q_H=3 trefoil saddle sits at the SAME targets, or at
    knot-specific ones, is open --- the script REPORTS V, sopt, sopt^6,
    J4/J2a, K_fb/J4 as diagnostics rather than asserting the Q=3 target.
  * WHAT WOULD MAKE THIS SUCCEED (the actual test): a beta at which the
    descent finds a genuine best-score plateau (score small and stable over
    many steps, not a fleeting dip on the way to collapse) at fixed J4/J0~1
    and r_bar~R0. That converged self-consistent saddle is the object whose
    geometry (winding, per-strand twist) is then physical and readable ---
    the prerequisite for pinning the E_6 tower level Dn_q (Paper XX, O2).

SMOKE TEST (fast, mechanical --- confirms the machinery runs; NOT physics):
  python3.11 qh3_feedback_saddle.py --N 24 --h 0.45 --NT_frame 2000 \
      --n_steps 60 --beta 0.452 --log_every 5 --outdir smoke_fb

REAL RUN (local; a beta scan is the informative one):
  python3.11 qh3_feedback_saddle.py --N 64 --h 0.175 \
      --beta_list "0.0,0.2,0.35,0.452,0.6,0.9" --n_steps 4000 \
      --outdir fb_saddle_scan
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
ap.add_argument('--n_steps',      type=int,   default=4000)
ap.add_argument('--beta',         type=float, default=0.452,
                help='Density-feedback coupling. J2iso_fb = int g2/(1+beta*g2) dV. '
                     'beta=0 recovers the feedback-free K=J2a+mu*J2iso. '
                     'Paper I value 0.452 is in Paper I kern units and may NOT '
                     'transfer to 3D g2 --- calibrate via --beta_list.')
ap.add_argument('--beta_list',    type=str,   default='',
                help='Comma-separated betas to scan, warm-started in sequence '
                     '(Paper I chain). Overrides --beta. The physical beta is '
                     'the one minimising the self-consistency score.')
ap.add_argument('--lr1',          type=float, default=3e-4,
                help='Phase 1 learning rate (cleanup, no angle clamp)')
ap.add_argument('--lr2',          type=float, default=1e-5,
                help='Phase 2 learning rate (saddle approach, with clamp)')
ap.add_argument('--delta_max_deg',type=float, default=9.0,
                help='Max rotation per grid point per step in Phase 2 (degrees)')
ap.add_argument('--K_rise_eps',   type=float, default=2.0,
                help='K_fb rise above its minimum to trigger Phase 2 transition')
ap.add_argument('--NT_frame',     type=int,   default=20000,
                help='Bishop-frame resolution for the construction-C IC. Lower '
                     '(e.g. 2000) for a fast smoke test.')
ap.add_argument('--log_every',    type=int,   default=10)
ap.add_argument('--warm_start',   type=str,   default=None,
                help='path to n_final.npy / n_best.npy to continue from')
ap.add_argument('--outdir',       type=str,   default='fb_saddle')
ap.add_argument('--seed',         type=int,   default=0)
ap.add_argument('--device',       type=str,   default='cpu')
args = ap.parse_args()

torch.manual_seed(args.seed)
np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

PHI    = (1+5**0.5)/2
MU     = 3.0 - PHI
LAM    = PHI**6
N, h   = args.N, args.h
R0, r0, C_star = args.R0, args.r0, args.C_star
dev    = torch.device(args.device)
DELTA_MAX = float(np.radians(args.delta_max_deg))

if args.beta_list.strip():
    BETAS = [float(b) for b in args.beta_list.split(',') if b.strip() != '']
else:
    BETAS = [args.beta]

print(f"{'='*72}")
print(f"  DENSITY-FEEDBACK SADDLE-SNAPSHOT SOLVER  (Q_H=3 trefoil)")
print(f"  E_fb = K_fb * J4,  K_fb = J2a + mu*J2iso_fb,")
print(f"  J2iso_fb = int g2/(1+beta*g2) dV   (Paper I feedback, ported to 3D)")
print(f"  saddle located by best score = |V-phi| + |sopt-1|  (Paper I method)")
print(f"  Grid: N={N}, h={h}, box=[{-N*h/2:.2f},{N*h/2:.2f}]   C*={C_star}")
print(f"  phi={PHI:.6f}  mu=3-phi={MU:.6f}  lam=phi^6={LAM:.5f}")
print(f"  beta(s) to run: {BETAS}")
print(f"{'='*72}")

# ── Grid ──────────────────────────────────────────────────────────────
cv      = h*(np.arange(N) - N//2 + 0.5)
pts_np  = np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),axis=-1).reshape(-1,3).astype(np.float32)
dist_from_origin = torch.tensor(
    np.linalg.norm(pts_np,axis=-1).reshape(N,N,N), dtype=torch.float32, device=dev)

needed_hw = R0 + r0 + 1/C_star + 0.5
actual_hw = N*h/2
if actual_hw < needed_hw:
    print(f"WARNING: box half-width {actual_hw:.2f} < recommended {needed_hw:.2f} "
          f"(fine for a smoke test, too small for physics).")
else:
    print(f"  Box adequate: half-width {actual_hw:.2f} >= {needed_hw:.2f}")

# ── Energy functional with density feedback ───────────────────────────
def E_fb(n, beta):
    """E_fb = K_fb*J4 with feedback-softened isotropic stiffness.
    Returns (E, K_fb, J4, rho_J4, J2a, J2iso_fb) --- all torch, E has grad."""
    nx,ny,nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u,a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2    = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a      = (s4*g2).sum() * h**3
    # density feedback: local energy density kern = g2 softens the isotropic term
    J2iso_fb = (g2/(1.0 + beta*g2)).sum() * h**3
    K_fb     = J2a + MU*J2iso_fb
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4    = rho_J4.sum() * h**3
    return K_fb*J4, K_fb, J4, rho_J4, J2a, J2iso_fb

def diagnostics(n, beta):
    """Self-consistency observables at the current field (no grad)."""
    with torch.no_grad():
        E,K_fb,J4,rho,J2a,J2iso_fb = E_fb(n, beta)
        r_bar = ((rho*dist_from_origin).sum()/rho.sum().clamp(1e-12)).item()
        K_fb, J4, J2a, J2iso_fb = K_fb.item(), J4.item(), J2a.item(), J2iso_fb.item()
        V     = J2iso_fb/J2a if J2a > 1e-12 else 0.0
        sopt  = (LAM*J4/K_fb)**0.5 if K_fb > 1e-12 else 0.0
        score = abs(V - PHI) + abs(sopt - 1.0)
    return dict(E=E.item(), K_fb=K_fb, J4=J4, J2a=J2a, J2iso_fb=J2iso_fb,
                r_bar=r_bar, V=V, sopt=sopt, sopt6=sopt**6, score=score,
                J4_over_J2a=(J4/J2a if J2a>1e-12 else 0.0),
                Kfb_over_J4=(K_fb/J4 if J4>1e-12 else 0.0))

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

# ── Per-strand Construction C initial condition (from gradient_flow_constrained) ──
def build_initial_field():
    print("\nBuilding per-strand construction-C initial condition "
          f"(NT_frame={args.NT_frame})...")
    t0b = time.time()
    t_frame, _, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=args.NT_frame)
    NTf = args.NT_frame
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
            d, idx = tree_l.query(qpts, workers=-1)
            d_per.append(d); t_per.append(t_l[idx])
        d_s = np.stack(d_per,axis=1); t_s = np.stack(t_per,axis=1)
        o   = np.argsort(d_s, axis=1)
        return (np.take_along_axis(t_s,o,axis=1)[:,0],
                np.take_along_axis(d_s,o,axis=1)[:,0],
                np.take_along_axis(t_s,o,axis=1)[:,1],
                np.take_along_axis(d_s,o,axis=1)[:,1])

    def frame_at_t(t_q):
        idx = np.searchsorted(t_frame, t_q%(2*np.pi)) % NTf
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
    return n0_np

# ── One beta: descend on E_fb, snapshot the best self-consistency score ──
def run_beta(beta, n0_np):
    n_t = torch.tensor(n0_np, dtype=torch.float32, device=dev)
    d0 = diagnostics(n_t, beta)
    vac0 = ((n_t[...,2]>0.95).float().mean()).item()
    print(f"\n{'─'*72}\n  beta = {beta:.5f}\n{'─'*72}")
    print(f"  init: E={d0['E']:.3e} K_fb={d0['K_fb']:.2f} J4={d0['J4']:.2f} "
          f"V={d0['V']:.4f}(phi={PHI:.4f}) sopt={d0['sopt']:.4f} "
          f"score={d0['score']:.4f} r_bar={d0['r_bar']:.3f}")

    n_param = n_t.clone().requires_grad_(True)
    opt1 = torch.optim.Adam([n_param], lr=args.lr1)
    opt2 = torch.optim.Adam([n_param], lr=args.lr2)
    phase = 1
    K_min_seen = d0['K_fb']
    best_score = d0['score']; n_best = n_t.clone().cpu().numpy(); best = dict(d0, step=0)
    history = []
    t_run0 = time.time()

    print(f"  {'step':>6} {'ph':>2} {'E_fb':>11} {'K_fb':>8} {'J4':>8} {'V':>7} "
          f"{'sopt':>7} {'score':>7} {'r_bar':>6} {'|grad|':>9} {'clmp%':>6}")
    for step in range(args.n_steps):
        opt = opt1 if phase == 1 else opt2
        n_before = n_param.detach().clone()
        opt.zero_grad()
        E, _, _, _, _, _ = E_fb(n_param, beta)
        E.backward()
        with torch.no_grad():
            n_param.grad.data.copy_(project_to_tangent(n_param.detach(), n_param.grad))
        # tangent |grad E_fb|: a saddle is a critical point (grad->0), so a true
        # saddle-approach shows |grad| DIP toward zero then RISE as descent flows
        # off the unstable direction. A monotone-small |grad| = genuine stationary
        # plateau; |grad| oscillating = ringing; |grad| large & falling = still moving.
        grad_norm = n_param.grad.norm().item()
        opt.step()
        with torch.no_grad():
            n_param.data.copy_(n_param / n_param.norm(dim=-1,keepdim=True).clamp(1e-10))

        clamped_frac = 0.0
        if phase == 2:
            with torch.no_grad():
                n_clamped = apply_angle_clamp(n_before, n_param.detach(), DELTA_MAX)
                diff = (n_clamped - n_param.detach()).norm(dim=-1)
                clamped_frac = (diff > 1e-6).float().mean().item()
                n_param.data.copy_(n_clamped)

        cs = step + 1
        if cs % args.log_every == 0 or step == 0:
            d = diagnostics(n_param, beta)
            d['grad_norm'] = grad_norm
            print(f"  {cs:>6} {phase:>2} {d['E']:>11.3e} {d['K_fb']:>8.1f} "
                  f"{d['J4']:>8.2f} {d['V']:>7.4f} {d['sopt']:>7.4f} "
                  f"{d['score']:>7.4f} {d['r_bar']:>6.3f} {grad_norm:>9.3e} "
                  f"{100*clamped_frac:>5.1f}%")
            d['step'] = cs; history.append(d)

            # saddle-snapshot: keep the best self-consistency score
            if d['score'] < best_score:
                best_score = d['score']
                n_best = n_param.detach().cpu().numpy().copy()
                best = dict(d)

            # phase transition (K_fb rises off its minimum) --- Paper I two-phase
            if phase == 1:
                if d['K_fb'] < K_min_seen:
                    K_min_seen = d['K_fb']
                elif d['K_fb'] > K_min_seen + args.K_rise_eps:
                    print(f"    *** PHASE 2 at step {cs}: K_fb={d['K_fb']:.2f} "
                          f"> K_min={K_min_seen:.2f}+{args.K_rise_eps} "
                          f"(lr={args.lr2}, clamp={args.delta_max_deg} deg) ***")
                    phase = 2
                    opt2 = torch.optim.Adam([n_param], lr=args.lr2)

            # dilution / collapse safeguards
            if d['r_bar'] > 4.5:
                print(f"    HALT: r_bar={d['r_bar']:.3f} > 4.5 (dilution) at step {cs}"); break
            if (n_param.detach()[...,2] > 0.95).float().mean().item() > vac0 + 0.10:
                print(f"    HALT: near-vacuum grew >+10pp (dilution) at step {cs}"); break
            if d['J4_over_J2a'] < 0.01 and step > 50:
                print(f"    HALT: J4/J2a collapsed at step {cs} (topology lost)"); break

    wall = time.time() - t_run0
    print(f"  best score={best_score:.5f} at step {best['step']}: "
          f"V={best['V']:.5f} sopt={best['sopt']:.5f} sopt6={best['sopt6']:.4f} "
          f"J4/J2a={best['J4_over_J2a']:.5f} K_fb/J4={best['Kfb_over_J4']:.4f} "
          f"r_bar={best['r_bar']:.3f}")
    print(f"  wall time {wall:.1f}s  ({wall/max(args.n_steps,1)*1000:.1f} ms/step)")

    # ── Plateau vs drift/ringing verdict over the last 25% of logged steps ──
    # Distinguishes (a) genuine metastable plateau, (b) slow monotone collapse,
    # (c) oscillation/ringing --- the thing 4000 steps could not resolve.
    tail = history[max(0, (3*len(history))//4):]
    if len(tail) >= 4:
        def stats(key):
            a = np.array([h[key] for h in tail], dtype=float)
            frac_drift = (a[-1]-a[0]) / (abs(a[0])+1e-12)   # net, signed
            frac_range = (a.max()-a.min()) / (abs(np.mean(a))+1e-12)  # peak-to-peak
            # zero-crossings of the step-to-step difference = oscillation count
            dif = np.diff(a); sign_flips = int(np.sum(dif[:-1]*dif[1:] < 0))
            return frac_drift, frac_range, sign_flips
        gN = np.array([h.get('grad_norm', np.nan) for h in tail], dtype=float)
        print(f"  ── tail diagnostics (last {len(tail)} logged pts, "
              f"steps {tail[0]['step']}–{tail[-1]['step']}) ──")
        print(f"     {'obs':>8} {'net drift':>11} {'pk-pk range':>12} {'sign flips':>11}")
        verdicts = {}
        for key in ('J4', 'V', 'score', 'r_bar'):
            fd, fr, sf = stats(key); verdicts[key] = (fd, fr, sf)
            print(f"     {key:>8} {fd:>+10.2%} {fr:>11.2%} {sf:>11}")
        print(f"     |grad| tail: {np.nanmin(gN):.3e} .. {np.nanmax(gN):.3e} "
              f"(trend {'FALLING' if gN[-1]<gN[0] else 'RISING'})")
        j4d = abs(verdicts['J4'][0]); j4flip = verdicts['J4'][2]
        if j4flip >= 3 and verdicts['J4'][1] > 0.05:
            verdict = "RINGING (J4 oscillating in the tail) — inertial-like, inspect period"
        elif j4d < 0.03 and verdicts['score'][1] < 0.05:
            verdict = "PLATEAU (J4 & score flat in tail) — metastable saddle candidate"
        elif verdicts['J4'][0] < -0.10:
            verdict = "SLOW COLLAPSE (J4 bleeding down in tail) — not yet stationary"
        else:
            verdict = "AMBIGUOUS — extend steps / inspect history"
        print(f"     VERDICT: {verdict}")
        best['tail_verdict'] = verdict
        best['tail_J4_drift'] = verdicts['J4'][0]

    tag = f"beta{beta:.3f}"
    np.save(os.path.join(args.outdir, f'n_best_{tag}.npy'), n_best)
    np.save(os.path.join(args.outdir, f'n_final_{tag}.npy'),
            n_param.detach().cpu().numpy())
    with open(os.path.join(args.outdir, f'history_{tag}.json'),'w') as f:
        json.dump(history, f, indent=2)
    best['beta'] = beta; best['wall_s'] = wall
    return n_best, best

# ── Main: (warm-started) beta scan ────────────────────────────────────
if args.warm_start:
    print(f"\nWarm-starting from {args.warm_start}...")
    n_seed = np.load(args.warm_start).astype(np.float32)
    if n_seed.shape != (N,N,N,3):
        print(f"FATAL: shape mismatch {n_seed.shape} vs ({N},{N},{N},3)"); sys.exit(1)
    n_seed /= np.linalg.norm(n_seed,axis=-1,keepdims=True).clip(1e-10)
else:
    n_seed = build_initial_field()

scan = []
for beta in BETAS:
    n_best, best = run_beta(beta, n_seed)
    scan.append(best)
    n_seed = n_best  # warm-start the next beta from this saddle (Paper I chain)

print(f"\n{'='*72}\n  SELF-CONSISTENCY SCAN vs beta "
      f"(physical beta = the one minimising score)\n{'='*72}")
print(f"  {'beta':>7} {'score':>8} {'V':>8} {'(phi)':>8} {'sopt':>7} {'sopt6':>7} "
      f"{'J4/J2a':>8} {'K_fb/J4':>8} {'r_bar':>6}")
for b in scan:
    print(f"  {b['beta']:>7.3f} {b['score']:>8.5f} {b['V']:>8.5f} {PHI:>8.5f} "
          f"{b['sopt']:>7.4f} {b['sopt6']:>7.4f} {b['J4_over_J2a']:>8.5f} "
          f"{b['Kfb_over_J4']:>8.4f} {b['r_bar']:>6.3f}")
best_overall = min(scan, key=lambda d: d['score'])
print(f"\n  -> best self-consistency at beta={best_overall['beta']:.3f} "
      f"(score={best_overall['score']:.5f}). Saddle field: "
      f"n_best_beta{best_overall['beta']:.3f}.npy")
print(f"  Reminder: V->phi and sopt->1 are the Q_H=2 targets; whether the Q_H=3\n"
      f"  saddle sits at the same targets is the open physics of this run.")
with open(os.path.join(args.outdir,'scan_summary.json'),'w') as f:
    json.dump(scan, f, indent=2)
print(f"\n  Scan summary: {os.path.join(args.outdir,'scan_summary.json')}")
