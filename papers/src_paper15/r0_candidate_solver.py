#!/usr/bin/env python3
"""
r0_candidate_solver.py  (v2 — corrected against the real solver files)
========================================================================
Empirical test of competing closed-form candidates for the trefoil minor
radius r0, by running Phase-1-only gradient flow (the ansatz-cleanup phase,
NOT the full saddle search) at each candidate r0 and comparing the resulting
K_min, J4, and E_geom = K*J4.

THIS VERSION IS CORRECTED against the actual project files:
  - bishop_frame_v2.py:  the Bishop frame construction (arc-length-uniform
    holonomy compensation) is reproduced VERBATIM from the real module,
    just parametrised over (R0, r0) instead of using module-level globals
    (the real bishop_frame_v2.py hardcodes R0=3.0, r0=0.874 at import time,
    which is exactly the thing we need to vary -- so the algorithm is
    copied faithfully but the hardcoded constants are turned into
    parameters; nothing about the transport/compensation method is changed).
  - gradient_flow_constrained.py: the energy functional E_geom = K*J4,
    K = J2a + mu*J2iso, the per-strand Construction C field, the Phase 1
    cleanup logic (lr1, no angle clamp), the K_rise_eps phase-transition
    criterion, and the three dilution safeguards (r_bar, near-vacuum
    fraction, J4/K collapse) are reproduced EXACTLY as in the real file,
    not reimplemented from memory.

RATIONALE FOR PHASE-1-ONLY (see conversation record):
  The Q_H=3 profile-direction landscape is FLAT (rational Jones polynomial
  at q_5: no irrational BPS anchor). Running into Phase 2, where Z3
  symmetry breaks spontaneously and stochastically, would inject seed-
  dependent noise into the r0 comparison for no benefit -- we already know
  Phase 2 just drifts toward (2+1) splitting regardless of r0. Phase 1
  ALONE removes the analytic-ansatz-specific noise (arctan profile,
  w=1/rho^2 weights, two-strand superposition artefacts) while leaving the
  GEOMETRY (R0, r0 -- gradient flow moves the FIELD, not the curve) and
  the TOPOLOGY (Q_H=3) untouched. This gives a genuine, ansatz-independent
  K_min(r0) at a small fraction of the cost of full two-phase runs.

CANDIDATES TESTED (edit CANDIDATES dict below to add more):
  control_0.874     : r0 = 0.874              (value used throughout the papers)
  sqrt2_over_phi    : r0 = sqrt(2)/phi        = 0.874032  (tightest numerical fit)
  R0_over_C2star    : r0 = R0/C2*             = 0.874172  (C2*=3.4318 derived
                       from Paper I's modified-BPS cubic thm:modBPS, set
                       against the independently-proved target 2^(4/3)/phi^5;
                       NOT dependent on the weaker Richardson g_infty estimate)
  phi2_over_3       : r0 = phi^2/3            = 0.872678
  seven_eighths     : r0 = 7/8                = 0.875000
  twoR0_over_phi4   : r0 = 2*R0/phi^4         = 0.875388  (Paper XVII solar-
                       angle match; equiv. to 2*R0/(3*phi+2))
  bracket_low/high  : 0.86, 0.90              (coarse bracket for overall shape)

USAGE:
  # Quick sanity check on a small grid first (a few minutes):
  python r0_candidate_solver.py --N 64 --h 0.175 --max_phase1_steps 600

  # Production run, matching the documented N=192 Phase-1-exit behaviour:
  python r0_candidate_solver.py --N 192 --h 0.05 --max_phase1_steps 1500

  # Subset of candidates only:
  python r0_candidate_solver.py --N 96 --h 0.075 \\
      --candidates control_0.874,R0_over_C2star,sqrt2_over_phi,twoR0_over_phi4

OUTPUT:
  Prints a table of K_min, J4_final, E_geom, exit_step for every candidate,
  sorted by r0. Saves results to r0_scan_results.json and (if matplotlib is
  available) a plot r0_scan.png.
"""

import numpy as np
import torch
import json
import time
import argparse
import sys
from scipy.spatial import KDTree

PHI = (1 + 5**0.5) / 2
MU  = 3.0 - PHI

# ── Candidate r0 values ──────────────────────────────────────────────────────
def build_candidates(R0):
    C2_star = 3.431820  # root of thm:modBPS cubic vs target 2^(4/3)/phi^5
    return {
        "control_0.874":     0.874000,
        "sqrt2_over_phi":    np.sqrt(2) / PHI,
        "R0_over_C2star":    R0 / C2_star,
        "phi2_over_3":       PHI**2 / 3,
        "seven_eighths":     7 / 8,
        "twoR0_over_phi4":   2 * R0 / PHI**4,
        "bracket_low_0.86":  0.860,
        "bracket_high_0.90": 0.900,
    }


# ══════════════════════════════════════════════════════════════════════════
# Bishop frame: faithful reproduction of bishop_frame_v2.py, parametrised
# over (R0, r0) instead of module-level globals.
# ══════════════════════════════════════════════════════════════════════════
def _Gamma(t, R0, r0):
    return np.array([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)])

def _Gamma_prime(t, R0, r0, h=1e-6):
    return (_Gamma(t+h, R0, r0) - _Gamma(t-h, R0, r0)) / (2*h)

def _build_raw_bishop_frame(R0, r0, NT=30000):
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
    tangents = np.array([_Gamma_prime(t, R0, r0) for t in t_arr])
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    T0 = tangents[0]
    arbitrary = np.array([0, 0, 1.0]) if abs(T0[2]) < 0.9 else np.array([1.0, 0, 0])
    N1_0 = arbitrary - np.dot(arbitrary, T0) * T0
    N1_0 /= np.linalg.norm(N1_0)
    N1 = np.zeros((NT, 3)); N1[0] = N1_0
    for i in range(1, NT):
        Tprev, Tcur = tangents[i-1], tangents[i]
        v = np.cross(Tprev, Tcur); s = np.linalg.norm(v); c = np.dot(Tprev, Tcur)
        if s < 1e-12:
            N1[i] = N1[i-1]; continue
        v_unit = v/s; angle = np.arctan2(s, c)
        vec = N1[i-1]
        rotated = (vec*np.cos(angle) + np.cross(v_unit, vec)*np.sin(angle)
                   + v_unit*np.dot(v_unit, vec)*(1-np.cos(angle)))
        N1[i] = rotated - np.dot(rotated, Tcur)*Tcur
        N1[i] /= np.linalg.norm(N1[i])
    N2 = np.cross(tangents, N1)
    return t_arr, tangents, N1, N2

def build_compensated_frame_arclength(R0, r0, NT=30000):
    """Verbatim port of bishop_frame_v2.build_compensated_frame_arclength,
    parametrised over (R0, r0) rather than module-level globals."""
    t_arr, T, N1_raw, N2_raw = _build_raw_bishop_frame(R0, r0, NT=NT)

    pts = np.array([_Gamma(t, R0, r0) for t in t_arr])
    seg = np.linalg.norm(np.diff(np.vstack([pts, pts[:1]]), axis=0), axis=1)
    cumArc = np.concatenate([[0], np.cumsum(seg)])
    L_total = cumArc[-1]

    cos_angle = np.clip(np.dot(N1_raw[0], N1_raw[-1]), -1, 1)
    sin_angle = np.dot(np.cross(N1_raw[0], N1_raw[-1]), T[0])
    H = np.arctan2(sin_angle, cos_angle)

    comp_angle = -H * cumArc[:NT] / L_total
    c, s = np.cos(comp_angle), np.sin(comp_angle)
    N1 = c[:, None]*N1_raw + s[:, None]*N2_raw
    N1 /= np.linalg.norm(N1, axis=1, keepdims=True)
    N2 = np.cross(T, N1)
    return t_arr, T, N1, N2, H


# ══════════════════════════════════════════════════════════════════════════
# Per-strand Construction C: verbatim port of the field-building section of
# gradient_flow_constrained.py, parametrised the same way the real CLI does.
# ══════════════════════════════════════════════════════════════════════════
def build_construction(N, h, R0, r0, C_star, NT_frame=20000, NT_curve=4000):
    t_frame, _, N1_frame, N2_frame, H = build_compensated_frame_arclength(R0, r0, NT=NT_frame)
    print(f"    Bishop frame holonomy: {np.degrees(H):.4f} deg")

    cv = h*(np.arange(N) - N//2 + 0.5)
    pts_np = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), axis=-1).reshape(-1, 3).astype(np.float32)

    t_arr = np.linspace(0, 2*np.pi, NT_curve, endpoint=False)
    arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
    Gx = (R0 + r0*np.cos(3*t_arr))*np.cos(2*t_arr)
    Gy = (R0 + r0*np.cos(3*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(3*t_arr)
    Gamma_pts = np.stack([Gx, Gy, Gz], axis=1)
    lobe_indices = [np.where((t_arr >= s) & (t_arr < s + 2*np.pi/3))[0] for s in arc_starts]
    lobe_trees = [KDTree(Gamma_pts[li]) for li in lobe_indices]
    lobe_t_arrays = [t_arr[li] for li in lobe_indices]

    def nearest_two_strands(qpts):
        d_per, t_per = [], []
        for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
            d, idx = tree_l.query(qpts, workers=1)
            d_per.append(d); t_per.append(t_l[idx])
        d_s = np.stack(d_per, axis=1); t_s = np.stack(t_per, axis=1)
        o = np.argsort(d_s, axis=1)
        return (np.take_along_axis(t_s, o, axis=1)[:, 0],
                np.take_along_axis(d_s, o, axis=1)[:, 0],
                np.take_along_axis(t_s, o, axis=1)[:, 1],
                np.take_along_axis(d_s, o, axis=1)[:, 1])

    def frame_at_t(t_q):
        idx = np.searchsorted(t_frame, t_q % (2*np.pi)) % NT_frame
        return N1_frame[idx], N2_frame[idx]

    def curve_at_t(t):
        return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                          (R0+r0*np.cos(3*t))*np.sin(2*t),
                          r0*np.sin(3*t)], axis=-1)

    t1_g, d1_g, t2_g, d2_g = nearest_two_strands(pts_np)
    chi1 = np.arctan2(np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[1], axis=1),
                       np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[0], axis=1))
    chi2 = np.arctan2(np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[1], axis=1),
                       np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[0], axis=1))
    Phi1 = chi1 + 3*t1_g; Phi2 = chi2 + 3*t2_g
    rho1 = np.clip(d1_g, 1e-6, None); rho2 = np.clip(d2_g, 1e-6, None)

    def f0(r): return 2*np.arctan(np.maximum(r, 1e-9)**(-C_star))
    f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
    w1 = 1/rho1**2; w2 = 1/rho2**2
    z1_np = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2_np = w1*np.sin(f1/2)*np.exp(1j*Phi1) + w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag = np.sqrt(np.abs(z1_np)**2 + np.abs(z2_np)**2)
    z1_np /= mag; z2_np /= mag
    nx0 = 2*np.real(np.conj(z1_np)*z2_np)
    ny0 = 2*np.imag(np.conj(z1_np)*z2_np)
    nz0 = np.abs(z1_np)**2 - np.abs(z2_np)**2
    n0_np = np.stack([nx0, ny0, nz0], axis=-1).reshape(N, N, N, 3).astype(np.float32)
    n0_np /= np.linalg.norm(n0_np, axis=-1, keepdims=True).clip(1e-10)
    return n0_np


# ══════════════════════════════════════════════════════════════════════════
# Energy functional: EXACT port of E_geom/diagnostics from
# gradient_flow_constrained.py (same variable names, same formulas).
# ══════════════════════════════════════════════════════════════════════════
def E_geom(n, h, dist_from_origin):
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    s4 = (1 - nz**2).clamp(0, 1)**2
    def cd(u, a): return (torch.roll(u, -1, a) - torch.roll(u, 1, a)) / (2*h)
    nxx, nxy, nxz = cd(nx, 0), cd(nx, 1), cd(nx, 2)
    nyx, nyy, nyz = cd(ny, 0), cd(ny, 1), cd(ny, 2)
    nzx, nzy, nzz = cd(nz, 0), cd(nz, 1), cd(nz, 2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a = (s4*g2).sum() * h**3
    J2iso = g2.sum() * h**3
    K = J2a + MU*J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4 = rho_J4.sum() * h**3
    return K*J4, K, J4, rho_J4

def diagnostics(n, h, dist_from_origin):
    with torch.no_grad():
        E, K, J4, rho = E_geom(n, h, dist_from_origin)
        r_bar = ((rho*dist_from_origin).sum()/rho.sum().clamp(1e-12)).item()
    return E.item(), K.item(), J4.item(), r_bar


# ══════════════════════════════════════════════════════════════════════════
# Phase 1 only: same loop structure, learning rate, K_rise_eps phase-
# transition criterion, and dilution safeguards as the real solver, but
# the loop simply STOPS at the moment the real solver would switch to
# Phase 2 (instead of continuing into Phase 2's angle-clamped approach).
# ══════════════════════════════════════════════════════════════════════════
def run_phase1_only(n0_np, h, R0, r0, lr1=3e-4, max_steps=1500,
                     K_rise_eps=2.0, log_every=10, device='cpu', seed=0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    dev = torch.device(device)
    N = n0_np.shape[0]

    cv = h*(np.arange(N) - N//2 + 0.5)
    pts_np = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), axis=-1).reshape(-1, 3).astype(np.float32)
    dist_from_origin = torch.tensor(
        np.linalg.norm(pts_np, axis=-1).reshape(N, N, N), dtype=torch.float32, device=dev)

    n_t = torch.tensor(n0_np, dtype=torch.float32, device=dev)
    E0, K0, J40, rbar0 = diagnostics(n_t, h, dist_from_origin)
    vac0 = ((n_t[..., 2] > 0.95).float().mean()).item()

    n_param = n_t.clone().requires_grad_(True)
    opt1 = torch.optim.Adam([n_param], lr=lr1)

    K_min_seen = K0
    K_min_step = 0
    history = []
    halt_reason = None

    for step in range(max_steps):
        opt1.zero_grad()
        E, _, _, _ = E_geom(n_param, h, dist_from_origin)
        E.backward()

        with torch.no_grad():
            grad = n_param.grad
            grad_proj = grad - (grad*n_param).sum(-1, keepdim=True)*n_param
            n_param.grad.data.copy_(grad_proj)

        opt1.step()

        with torch.no_grad():
            n_param.data.copy_(n_param / n_param.norm(dim=-1, keepdim=True).clamp(1e-10))

        current_step = step + 1

        if current_step % log_every == 0 or step == 0:
            E_val, K_val, J4_val, rbar_val = diagnostics(n_param, h, dist_from_origin)
            j4k = J4_val / max(K_val, 1e-6)
            near_vac = (n_param.detach()[..., 2] > 0.95).float().mean().item()
            history.append(dict(step=current_step, E=E_val, K=K_val, J4=J4_val,
                                 r_bar=rbar_val, J4_over_K=j4k))

            if K_val < K_min_seen:
                K_min_seen = K_val
                K_min_step = current_step
            elif K_val > K_min_seen + K_rise_eps:
                # This is exactly the real solver's Phase 1 -> Phase 2 trigger.
                # We stop here: Phase 1 exit reached.
                halt_reason = 'phase1_exit'
                break

            if rbar_val > 4.5:
                halt_reason = f'dilution_rbar_{rbar_val:.2f}'
                break
            if near_vac > vac0 + 0.10:
                halt_reason = f'dilution_vacuum_+{100*(near_vac-vac0):.1f}pp'
                break
            if j4k < 0.01 and step > 50:
                halt_reason = f'topology_lost_J4K_{j4k:.5f}'
                break

    with torch.no_grad():
        n_final = n_param / n_param.norm(dim=-1, keepdim=True).clamp(1e-10)
        E_f, K_f, J4_f, rbar_f = diagnostics(n_final, h, dist_from_origin)

    return dict(K_min=K_min_seen, K_min_step=K_min_step,
                K_final=K_f, J4_final=J4_f, E_geom_final=E_f,
                r_bar_final=rbar_f, exit_step=step+1,
                halt_reason=halt_reason or 'max_steps_reached',
                history=history)


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(description='Empirical r0 candidate scan (Phase 1 only, faithful to real solver)')
    ap.add_argument('--N', type=int, default=64, help='Grid size (matches real solver default; use 96-192 for production)')
    ap.add_argument('--h', type=float, default=0.175, help='Grid spacing (matches real solver default)')
    ap.add_argument('--R0', type=float, default=3.0)
    ap.add_argument('--C_star', type=float, default=2.5062)
    ap.add_argument('--lr1', type=float, default=3e-4)
    ap.add_argument('--max_phase1_steps', type=int, default=1500)
    ap.add_argument('--K_rise_eps', type=float, default=2.0,
                     help='Same default as the real solver: K rise above its minimum to trigger Phase 1 exit')
    ap.add_argument('--log_every', type=int, default=10)
    ap.add_argument('--device', type=str, default='cpu')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--candidates', type=str, default='',
                     help='Comma-separated subset of candidate names to run (default: all)')
    ap.add_argument('--outfile', type=str, default='r0_scan_results.json')
    args = ap.parse_args()

    candidates = build_candidates(args.R0)
    if args.candidates.strip():
        keep = set(s.strip() for s in args.candidates.split(','))
        candidates = {k: v for k, v in candidates.items() if k in keep}

    print(f"{'='*70}")
    print(f"  EMPIRICAL r0 CANDIDATE SCAN (Phase 1 only, faithful port)")
    print(f"  Grid: N={args.N}, h={args.h}, box=[{-args.N*args.h/2:.2f},{args.N*args.h/2:.2f}]")
    print(f"  R0={args.R0}, C*={args.C_star}, lr1={args.lr1}, "
          f"K_rise_eps={args.K_rise_eps}, max_steps={args.max_phase1_steps}")
    print(f"  Candidates: {list(candidates.keys())}")
    print(f"{'='*70}\n")

    results = {}
    for name, r0_val in candidates.items():
        print(f"\n--- Candidate: {name}  (r0 = {r0_val:.6f}) ---")
        needed_hw = args.R0 + r0_val + 1/args.C_star + 0.5
        actual_hw = args.N*args.h/2
        if actual_hw < needed_hw:
            print(f"  SKIPPED: box too small (need {needed_hw:.2f}, have {actual_hw:.2f})")
            continue

        t0 = time.time()
        print("  Building per-strand construction...")
        n0_np = build_construction(args.N, args.h, args.R0, r0_val, args.C_star)
        print(f"  Construction built in {time.time()-t0:.1f}s")

        t1 = time.time()
        res = run_phase1_only(n0_np, args.h, args.R0, r0_val,
                               lr1=args.lr1, max_steps=args.max_phase1_steps,
                               K_rise_eps=args.K_rise_eps, log_every=args.log_every,
                               device=args.device, seed=args.seed)
        print(f"  Phase 1 done in {time.time()-t1:.1f}s "
              f"({res['exit_step']} steps, halt: {res['halt_reason']}): "
              f"K_min={res['K_min']:.2f} at step {res['K_min_step']}, "
              f"J4_final={res['J4_final']:.2f}, "
              f"E_geom_final={res['E_geom_final']:.4e}")

        results[name] = dict(r0=r0_val, **{k: v for k, v in res.items() if k != 'history'})
        results[name]['history'] = res['history']

    print(f"\n{'='*70}")
    print("SUMMARY (sorted by r0)")
    print(f"{'='*70}")
    print(f"{'Candidate':<22}{'r0':>10}{'K_min':>10}{'J4_final':>10}{'E_geom':>12}{'exit_step':>10}  halt_reason")
    for name, r in sorted(results.items(), key=lambda kv: kv[1]['r0']):
        print(f"{name:<22}{r['r0']:>10.6f}{r['K_min']:>10.2f}{r['J4_final']:>10.2f}"
              f"{r['E_geom_final']:>12.4e}{r['exit_step']:>10d}  {r['halt_reason']}")

    with open(args.outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved full results to {args.outfile}")

    try:
        import matplotlib.pyplot as plt
        r0s = [r['r0'] for r in results.values()]
        Es = [r['E_geom_final'] for r in results.values()]
        Ks = [r['K_min'] for r in results.values()]
        order = np.argsort(r0s)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].plot(np.array(r0s)[order], np.array(Ks)[order], 'o-')
        axes[0].set_xlabel('r0'); axes[0].set_ylabel('K_min'); axes[0].set_title('K_min(r0)')
        axes[1].plot(np.array(r0s)[order], np.array(Es)[order], 'o-')
        axes[1].set_xlabel('r0'); axes[1].set_ylabel('E_geom_final'); axes[1].set_title('E_geom(r0)')
        for ax in axes:
            ax.axvline(0.874, color='gray', linestyle='--', alpha=0.5, label='control 0.874')
        plt.tight_layout()
        plt.savefig('r0_scan.png', dpi=120)
        print("Saved plot to r0_scan.png")
    except ImportError:
        print("(matplotlib not available; skipping plot)")


if __name__ == '__main__':
    main()
