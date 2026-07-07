#!/usr/bin/env python3
"""
crossing_transition_v2.py
=========================
Revised Q_H=2 → Q_H=3 crossing transition simulation.

WHAT WAS WRONG IN V1
---------------------
v1 used a local Rodrigues rotation to "un-wind" a crossing.
This is fundamentally broken: the Hopf charge Q_H is a
HOMOTOPY INVARIANT — ANY smooth field deformation preserves it.
Every λ in the v1 scan produced an identical Q_H=3 field.

CORRECT APPROACH
----------------
The near-trefoil Q_H=2 initial condition is built by changing
the phase winding in the per-strand doublet construction:

  Q_H=3 (full trefoil):    Phi = chi + 3*t   [standard]
  Q_H=2 (near-trefoil):    Phi = chi + 2*t   [one winding short]

Both use the same T(2,3) trefoil curve geometry but different
fiber phase winding. Changing q from 3→2 reduces Q_H by one
unit: the tube has the trefoil shape but the third crossing has
not been committed in the fiber. This is precisely the elastic
in the near-trefoil state just before the last positive crossing.

WHITEHEAD CHARGE
----------------
The Whitehead integral:
    Q_H = (1/4π²) ∫ A·(∇×A) d³x,   A_i = Im(Z† ∂_i Z)
is computed inline for any field via the inverse Hopf map
n(x) → Z(x) with gauge z1 ∈ ℝ≥0:
    z1(x) = √((1 + nz)/2)
    z2(x) = √((1 - nz)/2) · exp(i · arctan2(ny, nx))
This works on ARBITRARY n(x) fields, not just analytic constructions.
It is the only reliable way to track topology during gradient flow.

CROSSING EVENT SIGNATURE
------------------------
When Q_H=2 → Q_H=3 occurs in gradient flow (if it does), expect:
  - Whitehead Q_H jumps from ~2 to ~3 (possibly through non-integer)
  - J4 changes rapidly (crossing contributes new topological flux)
  - ρ_J4 at the target crossing vertex (C2) spikes transiently
  - r_bar changes as the flux redistributes

MODES
-----
  --mode compare   [DEFAULT] Build both Q_H=2 and Q_H=3 initial
                   conditions, compute Whitehead Q_H for each,
                   and run gradient flow. Compare trajectories.
                   Shows whether Q_H=2 spontaneously transitions.

  --mode winding2  Run only the Q_H=2 near-trefoil initial condition.
  --mode winding3  Run only the Q_H=3 full trefoil (reference).

USAGE
-----
  python crossing_transition_v2.py                     # compare both
  python crossing_transition_v2.py --mode winding2 --n_steps 2000
  python crossing_transition_v2.py --N 96 --n_steps 1000 --snap_every 50

  # Warm-start from Phase-1 exit of main run:
  python crossing_transition_v2.py --mode winding2 \\
         --warm_start gf_constrained/n_000850.npy
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

# ── CLI ───────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument('--mode',          type=str,   default='compare',
                choices=['compare','winding2','winding3'])
ap.add_argument('--N',             type=int,   default=64)
ap.add_argument('--h',             type=float, default=0.175)
ap.add_argument('--R0',            type=float, default=3.0)
ap.add_argument('--r0',            type=float, default=0.874)
ap.add_argument('--C_star',        type=float, default=2.5062)
ap.add_argument('--n_steps',       type=int,   default=800)
ap.add_argument('--lr',            type=float, default=3e-4)
ap.add_argument('--snap_every',    type=int,   default=100,
                help='Save n_XXXXXX.npy every N steps (0=disable)')
ap.add_argument('--log_every',     type=int,   default=10)
ap.add_argument('--whitehead_N',   type=int,   default=32,
                help='Grid size for inline Whitehead Q_H estimate '
                     '(32=fast ~2s, 48=medium ~10s, 64=slow ~30s). '
                     'Computed at step 0, n_steps/4, n_steps/2, n_steps.')
ap.add_argument('--outdir',        type=str,   default='crossing_v2')
ap.add_argument('--warm_start',    type=str,   default=None)
ap.add_argument('--device',        type=str,   default='cpu')
ap.add_argument('--seed',          type=int,   default=0)
args = ap.parse_args()

torch.manual_seed(args.seed)
np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

PHI = (1+5**0.5)/2
MU  = 3.0 - PHI
N, h = args.N, args.h
R0, r0, C_star = args.R0, args.r0, args.C_star
dev  = torch.device(args.device)

# Crossing vertex positions (where the ρ_J4 spike should appear)
T_CROSS = [np.pi/6 + 2*np.pi*k/3 for k in range(3)]
CROSS_PTS = np.array([
    [(R0+r0*np.cos(3*tc))*np.cos(2*tc),
     (R0+r0*np.cos(3*tc))*np.sin(2*tc),
      r0*np.sin(3*tc)]
    for tc in T_CROSS], dtype=np.float32)

# Whitehead Q_H snapshots at these fractional steps
WH_FRACS = [0.0, 0.25, 0.5, 1.0]

cv     = h*(np.arange(N) - N//2 + 0.5)
pts_np = np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),axis=-1).reshape(-1,3).astype(np.float32)
dist_from_origin = torch.tensor(
    np.linalg.norm(pts_np,axis=-1).reshape(N,N,N), dtype=torch.float32, device=dev)

print(f"{'='*68}")
print(f"  CROSSING TRANSITION v2 — Q_H=2 near-trefoil vs Q_H=3 trefoil")
print(f"  Mode: {args.mode.upper()}")
print(f"  Grid: N={N}, h={h}, box=[{-N*h/2:.2f},{N*h/2:.2f}]")
print(f"  C*={C_star}, R0={R0}, r0={r0}")
print(f"  Steps: {args.n_steps}, lr={args.lr}")
print(f"  Whitehead Q_H at steps: "
      f"{[int(f*args.n_steps) for f in WH_FRACS]}"
      f"  (grid N_wh={args.whitehead_N})")
print(f"{'='*68}")


# ══════════════════════════════════════════════════════════════════════
# ENERGY FUNCTIONAL
# ══════════════════════════════════════════════════════════════════════
def E_geom(n):
    nx,ny,nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u,a): return (torch.roll(u,-1,a)-torch.roll(u,1,a))/(2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2  = (nxx**2+nxy**2+nxz**2+nyx**2+nyy**2+nyz**2+nzx**2+nzy**2+nzz**2)
    K   = (s4*g2).sum()*h**3 + MU*g2.sum()*h**3
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho = Fxy**2+Fxz**2+Fyz**2
    J4  = rho.sum()*h**3
    return K*J4, K.item(), J4.item(), rho

def diag(n):
    with torch.no_grad():
        E,K,J4,rho = E_geom(n)
        rb = ((rho*dist_from_origin).sum()/rho.sum().clamp(1e-12)).item()
    return E.item(), K, J4, rb

def ptangent(n, g): return g - (g*n).sum(-1,keepdim=True)*n


# ══════════════════════════════════════════════════════════════════════
# WHITEHEAD Q_H — from arbitrary n(x) field
# ══════════════════════════════════════════════════════════════════════
def whitehead_from_n(n_np, N_wh=32, L_wh=6.5):
    """
    Compute Hopf charge Q_H from an arbitrary n(x) field via the
    Whitehead Berry–Chern–Simons integral.

    Inverse Hopf map with gauge z1 ∈ ℝ≥0:
        z1 = √((1 + nz)/2)
        z2 = √((1 - nz)/2) · exp(i · arctan2(ny, nx))

    This gauge has a line singularity at nz=−1 (south pole of S²),
    but the integrand A·(∇×A) is smooth away from that set
    and the integral converges.

    n_np: shape (N_field, N_field, N_field, 3) float32
    """
    h_wh  = 2*L_wh/(N_wh - 1)
    ax_wh = np.linspace(-L_wh, L_wh, N_wh)
    X_wh,Y_wh,Z_wh = np.meshgrid(ax_wh,ax_wh,ax_wh,indexing='ij')
    query_pts = np.stack([X_wh.ravel(),Y_wh.ravel(),Z_wh.ravel()],axis=1)

    # Interpolate n_np onto the whitehead grid using nearest-neighbour
    # (trilinear would be better but NN is fast and sufficient for Q_H estimate)
    N_f = n_np.shape[0]
    h_f = h  # field grid spacing

    def world_to_idx(qpts, N_field, h_field):
        """Map world coords to nearest grid indices (clipped)."""
        cv_f = h_field*(np.arange(N_field) - N_field//2 + 0.5)
        idx = np.round((qpts - cv_f[0]) / h_field).astype(int)
        return np.clip(idx, 0, N_field-1)

    ix = world_to_idx(query_pts[:,0], N_f, h_f)
    iy = world_to_idx(query_pts[:,1], N_f, h_f)
    iz = world_to_idx(query_pts[:,2], N_f, h_f)
    n_flat = n_np.reshape(-1, 3)
    n_q = n_flat[ix*N_f*N_f + iy*N_f + iz]   # (N_wh^3, 3)

    nx_q = n_q[:,0]; ny_q = n_q[:,1]; nz_q = n_q[:,2]

    # Inverse Hopf map
    nz_c = np.clip(nz_q, -1+1e-7, 1-1e-7)
    z1 = np.sqrt((1 + nz_c)/2)                        # real, ≥0
    z2 = np.sqrt((1 - nz_c)/2)*np.exp(1j*np.arctan2(ny_q, nx_q))

    z1 = z1.reshape(N_wh,N_wh,N_wh)
    z2 = z2.reshape(N_wh,N_wh,N_wh)

    # Berry connection A_i = Im(Z† ∂_i Z) = Im(z1*∂_i z1 + z2*conj∂_i z2)
    #  For real z1: Im(z1*∂_i z1) = z1 * Im(∂_i z1) = 0 (since z1 real ⟹ ∂_i z1 real)
    #  So A_i = Im(conj(z2)*∂_i z2)
    A = [np.imag(np.conj(z1)*np.gradient(z1.astype(complex),h_wh,axis=a)
               + np.conj(z2)*np.gradient(z2,h_wh,axis=a)) for a in range(3)]

    curlA = [
        np.gradient(A[2],h_wh,axis=1) - np.gradient(A[1],h_wh,axis=2),
        np.gradient(A[0],h_wh,axis=2) - np.gradient(A[2],h_wh,axis=0),
        np.gradient(A[1],h_wh,axis=0) - np.gradient(A[0],h_wh,axis=1),
    ]
    integrand = sum(A[i]*curlA[i] for i in range(3))
    return float(h_wh**3 * np.sum(integrand) / (4*np.pi**2))


# ══════════════════════════════════════════════════════════════════════
# PER-STRAND CONSTRUCTION — parameterised by winding q
# ══════════════════════════════════════════════════════════════════════
def build_construction(q_winding):
    """
    Build the per-strand initial field on the T(2,3) trefoil curve.

    q_winding=3: Phi = chi + 3*t → Q_H=3 (standard trefoil)
    q_winding=2: Phi = chi + 2*t → Q_H=2 (near-trefoil, one winding short)

    The Q_H value is verified by the Whitehead integral at startup.
    """
    t0 = time.time()
    label = f'Q_H={q_winding} (Phi=chi+{q_winding}*t)'
    print(f"\n  Building construction: {label}")

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
    lobe_idx  = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
    lobe_trees = [KDTree(Gamma_pts[li]) for li in lobe_idx]
    lobe_t    = [t_arr[li] for li in lobe_idx]

    def nearest_two(qpts):
        d_per, t_per = [], []
        for tree_l, tl in zip(lobe_trees, lobe_t):
            d, idx = tree_l.query(qpts, workers=1)
            d_per.append(d); t_per.append(tl[idx])
        ds = np.stack(d_per,axis=1); ts = np.stack(t_per,axis=1)
        o  = np.argsort(ds, axis=1)
        return (np.take_along_axis(ts,o,axis=1)[:,0],
                np.take_along_axis(ds,o,axis=1)[:,0],
                np.take_along_axis(ts,o,axis=1)[:,1],
                np.take_along_axis(ds,o,axis=1)[:,1])

    def frame_at_t(tq):
        idx = np.searchsorted(t_frame, tq%(2*np.pi)) % NT_frame
        return N1_frame[idx], N2_frame[idx]

    def curve_at_t(t):
        return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                         (R0+r0*np.cos(3*t))*np.sin(2*t),
                          r0*np.sin(3*t)], axis=-1)

    t1_g, d1_g, t2_g, d2_g = nearest_two(pts_np)
    chi1 = np.arctan2(np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[1],axis=1),
                      np.sum((pts_np-curve_at_t(t1_g))*frame_at_t(t1_g)[0],axis=1))
    chi2 = np.arctan2(np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[1],axis=1),
                      np.sum((pts_np-curve_at_t(t2_g))*frame_at_t(t2_g)[0],axis=1))

    # ── KEY CHANGE: q_winding controls Q_H ──────────────────────────
    Phi1 = chi1 + q_winding * t1_g
    Phi2 = chi2 + q_winding * t2_g

    def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    rho1 = np.clip(d1_g,1e-6,None); rho2 = np.clip(d2_g,1e-6,None)
    f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
    w1 = 1/rho1**2;       w2 = 1/rho2**2
    z1 = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2 = w1*np.sin(f1/2)*np.exp(1j*Phi1) + w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2)
    z1 /= mag; z2 /= mag
    nx0 = 2*np.real(np.conj(z1)*z2)
    ny0 = 2*np.imag(np.conj(z1)*z2)
    nz0 = np.abs(z1)**2 - np.abs(z2)**2
    n0  = np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
    n0 /= np.linalg.norm(n0,axis=-1,keepdims=True).clip(1e-10)
    print(f"  Built in {time.time()-t0:.1f}s")
    return n0, label


# ══════════════════════════════════════════════════════════════════════
# CROSSING SITE ρ_J4 MONITOR
# ══════════════════════════════════════════════════════════════════════
def crossing_site_rho(rho_t, radius=2):
    """J4 density integrated in a ball of `radius` grid cells around each crossing."""
    rho_np = rho_t.detach().cpu().numpy().reshape(-1)
    vals = []
    for cp in CROSS_PTS:
        ix = int(np.round((cp[0]-cv[0])/h)); ix = np.clip(ix,0,N-1)
        iy = int(np.round((cp[1]-cv[0])/h)); iy = np.clip(iy,0,N-1)
        iz = int(np.round((cp[2]-cv[0])/h)); iz = np.clip(iz,0,N-1)
        s = 0.0
        for dx in range(-radius,radius+1):
            for dy in range(-radius,radius+1):
                for dz in range(-radius,radius+1):
                    if dx**2+dy**2+dz**2 > radius**2: continue
                    s += float(rho_np[((ix+dx)%N)*N*N+((iy+dy)%N)*N+(iz+dz)%N])
        vals.append(s*h**3)
    return vals


# ══════════════════════════════════════════════════════════════════════
# SINGLE RUN
# ══════════════════════════════════════════════════════════════════════
def run(n0_np, label, q_winding, outdir):
    os.makedirs(outdir, exist_ok=True)
    n_np = n0_np.copy()
    if args.warm_start:
        print(f"\n  Warm-starting from {args.warm_start}")
        n_np = np.load(args.warm_start).astype(np.float32)
        n_np /= np.linalg.norm(n_np,axis=-1,keepdims=True).clip(1e-10)

    n_t = torch.tensor(n_np, dtype=torch.float32, device=dev)

    # ── Whitehead Q_H of initial construction ──────────────────────
    print(f"\n  Computing Whitehead Q_H of initial field "
          f"(N_wh={args.whitehead_N})...")
    t_wh = time.time()
    QH_init = whitehead_from_n(n_np, N_wh=args.whitehead_N)
    print(f"  Whitehead Q_H = {QH_init:.3f}  "
          f"(expected {q_winding}; t={time.time()-t_wh:.1f}s)")

    # ── Initial gradient-flow diagnostics ──────────────────────────
    E0, K0, J40, rb0 = diag(n_t)
    vac0 = (n_t[...,2] > 0.95).float().mean().item()
    with torch.no_grad():
        _, _, _, rho0 = E_geom(n_t)
        C0 = crossing_site_rho(rho0)

    print(f"\n  Initial: E={E0:.4e}  K={K0:.1f}  J4={J40:.2f}  "
          f"J4/K={J40/max(K0,1e-9):.4f}  r_bar={rb0:.3f}")
    print(f"  ρ_J4 at crossings: C0={C0[0]:.4f}  C1={C0[1]:.4f}  C2={C0[2]:.4f}")

    # ── Optimiser ──────────────────────────────────────────────────
    n_param = n_t.clone().requires_grad_(True)
    opt = torch.optim.Adam([n_param], lr=args.lr)

    history = []
    wh_steps = sorted(set(int(f*args.n_steps) for f in WH_FRACS))
    wh_results = {0: QH_init}   # step → Q_H value

    print(f"\n  {'step':>6}  {'E':>12}  {'K':>8}  {'J4':>8}  "
          f"{'J4/K':>7}  {'r_bar':>7}  {'ρ_J4@C2':>9}  {'|grad|':>10}"
          f"  {'Q_H_wh':>7}")
    print(f"  {'-'*88}")

    t_run0 = time.time()
    for step in range(args.n_steps):
        opt.zero_grad()
        E, _, _, _ = E_geom(n_param)
        E.backward()
        with torch.no_grad():
            n_param.grad.data.copy_(ptangent(n_param.detach(), n_param.grad))
        gn = n_param.grad.norm().item()
        opt.step()
        with torch.no_grad():
            n_param.data.copy_(n_param/n_param.norm(dim=-1,keepdim=True).clamp(1e-10))

        cur = step + 1   # actual step count (1-indexed)

        # Snapshot
        if args.snap_every > 0 and cur % args.snap_every == 0:
            fname = os.path.join(outdir, f'n_{cur:06d}.npy')
            np.save(fname, n_param.detach().cpu().numpy())
            print(f"    >> snapshot: {fname}")

        # Inline Whitehead Q_H at requested steps
        wh_str = ''
        if cur in wh_steps and cur not in wh_results:
            tw = time.time()
            qh = whitehead_from_n(n_param.detach().cpu().numpy(), N_wh=args.whitehead_N)
            wh_results[cur] = qh
            wh_str = f'{qh:>7.3f}'
            print(f"    [Whitehead Q_H at step {cur}: {qh:.3f}  ({time.time()-tw:.1f}s)]")
        else:
            wh_str = f"{'---':>7}"

        # Diagnostics
        if cur % args.log_every == 0 or step == 0:
            Ev,Kv,J4v,rbv = diag(n_param)
            with torch.no_grad():
                _,_,_,rhov = E_geom(n_param)
                Cv = crossing_site_rho(rhov)
                vac = (n_param[...,2]>0.95).float().mean().item()

            print(f"  {cur:>6}  {Ev:>12.4e}  {Kv:>8.1f}  {J4v:>8.2f}  "
                  f"{J4v/max(Kv,1e-9):>7.4f}  {rbv:>7.3f}  "
                  f"{Cv[2]:>9.4f}  {gn:>10.4e}  {wh_str}")

            history.append(dict(step=cur, E=Ev, K=Kv, J4=J4v,
                                J4_over_K=J4v/max(Kv,1e-9), r_bar=rbv,
                                rho_C0=Cv[0],rho_C1=Cv[1],rho_C2=Cv[2],
                                grad_norm=gn, near_vac=vac))

            # Halt checks
            if rbv > 4.5:
                print(f"\n  HALT: r_bar={rbv:.3f} > 4.5 (dilution)")
                break
            if vac > vac0 + 0.10:
                print(f"\n  HALT: near-vacuum grew +{100*(vac-vac0):.1f}pp")
                break
            if J4v/max(Kv,1e-9) < 0.005 and step > 30:
                print(f"\n  HALT: J4/K collapsed (topology lost at step {cur})")
                break

    Ef,Kf,J4f,rbf = diag(n_param)
    QH_final = whitehead_from_n(n_param.detach().cpu().numpy(), N_wh=args.whitehead_N)
    wh_results[cur] = QH_final

    print(f"\n  Final: E={Ef:.4e}  K={Kf:.1f}  J4={J4f:.2f}  "
          f"J4/K={J4f/max(Kf,1e-9):.4f}  r_bar={rbf:.3f}")
    print(f"  Whitehead Q_H: initial={QH_init:.3f} → final={QH_final:.3f}")
    print(f"  Wall time: {time.time()-t_run0:.1f}s")

    # Check for crossing event
    QH_vals = [wh_results[s] for s in sorted(wh_results.keys())]
    dQH = max(QH_vals) - min(QH_vals)
    if dQH > 0.3:
        print(f"\n  *** TOPOLOGY CHANGE DETECTED: ΔQ_H = {dQH:.3f} ***")
        print(f"  Q_H trajectory: "
              + " → ".join(f"{wh_results[s]:.2f}@{s}" for s in sorted(wh_results.keys())))
    else:
        print(f"\n  No topology change detected (ΔQ_H = {dQH:.3f} < 0.3)")

    # Save outputs
    np.save(os.path.join(outdir, 'n_final.npy'), n_param.detach().cpu().numpy())
    with open(os.path.join(outdir, 'log.json'), 'w') as f:
        json.dump(history, f, indent=2)
    with open(os.path.join(outdir, 'whitehead_trajectory.json'), 'w') as f:
        json.dump({str(k): v for k,v in wh_results.items()}, f, indent=2)
    print(f"  Saved: {outdir}/n_final.npy, log.json, whitehead_trajectory.json")

    return dict(label=label, QH_init=QH_init, QH_final=QH_final, dQH=dQH,
                J4_f=J4f, K_f=Kf, r_bar_f=rbf, steps=cur)


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
windings = []
if args.mode == 'winding2':  windings = [2]
elif args.mode == 'winding3': windings = [3]
else:                         windings = [2, 3]   # compare both

results = []
for q in windings:
    n0, lbl = build_construction(q)
    outd = os.path.join(args.outdir, f'q{q}')
    r = run(n0, lbl, q, outd)
    results.append(r)

if len(results) == 2:
    r2, r3 = results[0], results[1]
    print(f"\n{'='*68}")
    print(f"  COMPARISON SUMMARY")
    print(f"{'='*68}")
    print(f"  {'':20}  {'Q_H=2 near-trefoil':>22}  {'Q_H=3 trefoil':>16}")
    print(f"  {'Whitehead Q_H init':20}  {r2['QH_init']:>22.3f}  {r3['QH_init']:>16.3f}")
    print(f"  {'Whitehead Q_H final':20}  {r2['QH_final']:>22.3f}  {r3['QH_final']:>16.3f}")
    print(f"  {'ΔQ_H':20}  {r2['dQH']:>22.3f}  {r3['dQH']:>16.3f}")
    print(f"  {'J4 final':20}  {r2['J4_f']:>22.2f}  {r3['J4_f']:>16.2f}")
    print(f"  {'J4/K final':20}  {r2['J4_f']/max(r2['K_f'],1e-9):>22.4f}  "
          f"{r3['J4_f']/max(r3['K_f'],1e-9):>16.4f}")
    print(f"  {'r_bar final':20}  {r2['r_bar_f']:>22.3f}  {r3['r_bar_f']:>16.3f}")
    print()
    if abs(r2['QH_init'] - 2) > 0.5:
        print(f"  WARNING: Q_H=2 construction gave Q_H={r2['QH_init']:.3f}.")
        print(f"  Try --whitehead_N 48 for better resolution.")
    if r2['dQH'] > 0.3:
        print(f"  *** CROSSING EVENT DETECTED in Q_H=2 run ***")
        print(f"  The near-trefoil Q_H=2 field spontaneously transitioned during flow.")
        print(f"  Check whitehead_trajectory.json for the step-by-step Q_H history.")
    else:
        print(f"  No crossing event detected in Q_H=2 run ({r2['steps']} steps).")
        print(f"  The Q_H=2 near-trefoil is stable. More steps may be needed, or")
        print(f"  the Q_H=2→Q_H=3 transition requires external energy input.")
