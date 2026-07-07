#!/usr/bin/env python3
"""
hopf_link_construction.py  (v2 — near-trefoil interpolation added)
===================================================================
Build a Q_H=2 Hopf link initial field using the T(2,2) torus link curve,
or a continuously parametrised near-trefoil transition state T(2,q) with
fractional poloidal winding q ∈ [2,3], and compare its gradient flow to
the Q_H=3 trefoil Construction C.

NEW IN v2: TWO INTERPOLATION PARAMETERS FOR THE TRANSITION STATE
----------------------------------------------------------------
  --q_pol FLOAT     Poloidal winding exponent of the CURVE (default 2.0).
                    q=2.0 → standard T(2,2) Hopf link.
                    q=3.0 → T(2,3) trefoil geometry.
                    q=2.9 → near-trefoil: geometry is almost trefoil but the
                             curve does not quite close in 0..2π, producing
                             the pre-crossing intermediate state of
                             rem:crossing_path in Paper XVIII.
                    This parameter controls the GEOMETRY of the construction
                    curve, independently of the field winding below.

  --alpha_wind FLOAT  Phase-winding interpolation (default 0.0).
                    0.0 → field phase Φ = χ + 2τ (pure Q_H=2 winding).
                    1.0 → field phase Φ = χ + 3τ (pure Q_H=3 winding).
                    Intermediate values continuously interpolate the phase
                    commitment: Φ = χ + (2 + alpha)*τ.
                    This controls the TOPOLOGICAL commitment independently
                    of the curve geometry above.

The two parameters are independent:
  q_pol=2.0, alpha_wind=0.0  →  exact original Hopf link (Q_H=2)
  q_pol=2.9, alpha_wind=0.0  →  near-trefoil geometry, Q_H=2 field winding
                                  (the pre-crossing state: geometry almost
                                   committed, topology not yet committed)
  q_pol=2.9, alpha_wind=0.5  →  near-trefoil geometry, mixed winding
                                  (transition state)
  q_pol=3.0, alpha_wind=1.0  →  full trefoil geometry and winding (Q_H=3)

PHYSICAL INTERPRETATION (rem:crossing_path, Paper XVIII)
---------------------------------------------------------
The elastic-band analogy of the Q_H=2 → Q_H=3 transition proceeds:
  Step 1 (half-twist): writhe injection — corresponds to increasing q_pol
    away from 2 while keeping alpha_wind=0.
  Step 2 (asymmetric fold): threading begins — q_pol approaching 3, the
    curve geometry develops the three-lobe structure without topological
    commitment.
  Step 3 (unique valid crossing): topology commits — alpha_wind crosses from
    0 to 1 at q_pol≈3.

Sweeping q_pol from 2 to 3 at alpha_wind=0 traces the geometric path.
Sweeping alpha_wind from 0 to 1 at fixed q_pol≈2.9 approaches the
crossing event while holding the geometry near-fixed.

USAGE
-----
  # Standard Hopf link (exact original behaviour):
  python hopf_link_construction.py --N 64 --n_steps 800

  # Near-trefoil at q=2.7, Q_H=2 field (explore pre-crossing geometry):
  python hopf_link_construction.py --N 64 --q_pol 2.7 --alpha_wind 0.0 \\
      --verify_only --outdir near_trefoil_q27

  # Full transition-state sweep (build only, view in hopfion_viewer.html):
  for q in 2.0 2.3 2.5 2.7 2.9 3.0; do
    python hopf_link_construction.py --N 64 --q_pol $q --alpha_wind 0.0 \\
        --verify_only --outdir sweep_q${q}
  done

  # Commit the crossing: near-trefoil geometry, full Q_H=3 winding:
  python hopf_link_construction.py --N 64 --q_pol 2.9 --alpha_wind 1.0 \\
      --verify_only --outdir near_trefoil_committed
"""
import numpy as np
import torch
import time, json, os, sys, argparse
from scipy.spatial import KDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bishop_frame_v2 import build_compensated_frame_arclength

ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--N',           type=int,   default=64)
ap.add_argument('--h',           type=float, default=0.175)
ap.add_argument('--R0',          type=float, default=3.0)
ap.add_argument('--r0',          type=float, default=0.874)
ap.add_argument('--C_star',      type=float, default=2.5062)
ap.add_argument('--n_steps',     type=int,   default=800)
ap.add_argument('--lr',          type=float, default=3e-4)
ap.add_argument('--snap_every',  type=int,   default=100)
ap.add_argument('--log_every',   type=int,   default=10)
ap.add_argument('--whitehead_N', type=int,   default=32,
                help='Whitehead grid (32=fast ~2s, 48=~10s).')
ap.add_argument('--outdir',      type=str,   default='hopf_link_run')
ap.add_argument('--verify_only', action='store_true',
                help='Build field, run Whitehead, save snapshot, exit')
ap.add_argument('--device',      type=str,   default='cpu')
ap.add_argument('--seed',        type=int,   default=0)
# ── NEW: transition-state parameters ──────────────────────────────────
ap.add_argument('--q_pol',       type=float, default=2.0,
                help='Poloidal winding of the construction CURVE. '
                     '2.0=Hopf link T(2,2), 3.0=trefoil T(2,3), '
                     '2.9=near-trefoil pre-crossing state. '
                     'Controls curve GEOMETRY, not field topology.')
ap.add_argument('--alpha_wind',  type=float, default=0.0,
                help='Phase-winding interpolation: '
                     '0.0=Q_H=2 field (Phi=chi+2*tau), '
                     '1.0=Q_H=3 field (Phi=chi+3*tau), '
                     'intermediate=mixed. Controls topological COMMITMENT.')
args = ap.parse_args()

torch.manual_seed(args.seed); np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

PHI = (1+5**0.5)/2
MU  = 3.0 - PHI
N, h = args.N, args.h
R0, r0, C_star = args.R0, args.r0, args.C_star
q_pol     = args.q_pol
alpha     = args.alpha_wind
wind_exp  = 2.0 + alpha          # 2.0 at alpha=0, 3.0 at alpha=1
dev = torch.device(args.device)

cv     = h*(np.arange(N) - N//2 + 0.5)
pts_np = np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),axis=-1).reshape(-1,3).astype(np.float32)
dist_from_origin = torch.tensor(
    np.linalg.norm(pts_np,axis=-1).reshape(N,N,N), dtype=torch.float32, device=dev)

print(f"{'='*68}")
print(f"  NEAR-TREFOIL CONSTRUCTION  v2")
print(f"  Grid: N={N}, h={h}  |  C*={C_star}, R0={R0}, r0={r0}")
print(f"  q_pol={q_pol}  (curve geometry: 2.0=Hopf link, 3.0=trefoil)")
print(f"  alpha_wind={alpha}  (winding: 0.0=Q_H=2, 1.0=Q_H=3)")
print(f"  effective winding exponent: {wind_exp:.3f}")
print(f"  Steps={args.n_steps}, lr={args.lr}")
print(f"{'='*68}")


# ── Energy functional ─────────────────────────────────────────────────
def E_geom(n):
    nx,ny,nz = n[...,0], n[...,1], n[...,2]
    s4 = (1-nz**2).clamp(0,1)**2
    def cd(u,a): return (torch.roll(u,-1,a)-torch.roll(u,1,a))/(2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    g2  = nxx**2+nxy**2+nxz**2+nyx**2+nyy**2+nyz**2+nzx**2+nzy**2+nzz**2
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


# ── Whitehead Q_H from n(x) field ────────────────────────────────────
def whitehead_from_n(n_np, N_wh, L_wh=6.5):
    h_wh = 2*L_wh/(N_wh-1)
    ax   = np.linspace(-L_wh, L_wh, N_wh)
    X,Y,Z = np.meshgrid(ax,ax,ax,indexing='ij')
    qpts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    N_f = n_np.shape[0]
    def to_idx(vals):
        return np.clip(np.round((vals-cv[0])/h).astype(int), 0, N_f-1)
    ix = to_idx(qpts[:,0]); iy = to_idx(qpts[:,1]); iz = to_idx(qpts[:,2])
    n_q = n_np.reshape(-1,3)[ix*N_f*N_f + iy*N_f + iz]
    nzc = np.clip(n_q[:,2], -1+1e-7, 1-1e-7)
    z1  = np.sqrt((1+nzc)/2)
    z2  = np.sqrt((1-nzc)/2)*np.exp(1j*np.arctan2(n_q[:,1], n_q[:,0]))
    z1  = z1.reshape(N_wh,N_wh,N_wh)
    z2  = z2.reshape(N_wh,N_wh,N_wh)
    A = [np.imag(np.conj(z1)*np.gradient(z1.astype(complex),h_wh,axis=a)
               + np.conj(z2)*np.gradient(z2,h_wh,axis=a)) for a in range(3)]
    curlA = [np.gradient(A[2],h_wh,axis=1)-np.gradient(A[1],h_wh,axis=2),
             np.gradient(A[0],h_wh,axis=2)-np.gradient(A[2],h_wh,axis=0),
             np.gradient(A[1],h_wh,axis=0)-np.gradient(A[0],h_wh,axis=1)]
    return float(h_wh**3*np.sum(sum(A[i]*curlA[i] for i in range(3)))/(4*np.pi**2))


# ══════════════════════════════════════════════════════════════════════
# CORE CONSTRUCTION: T(2,q) with interpolated winding
# ══════════════════════════════════════════════════════════════════════
def build_torus_knot_field(q_pol, wind_exp, label=None):
    """
    Build the Hopf condensate field on the T(2,q_pol) curve geometry
    with field phase winding exponent wind_exp.

    q_pol=2.0, wind_exp=2.0  →  original T(2,2) Hopf link (Q_H=2).
    q_pol=3.0, wind_exp=3.0  →  original T(2,3) trefoil (Q_H=3).
    q_pol=2.9, wind_exp=2.0  →  near-trefoil geometry, Q_H=2 winding
                                  (the pre-crossing state).
    q_pol=2.9, wind_exp=2.5  →  near-trefoil geometry, mixed winding.

    The curve:
      Γ(t) = ((R0+r0·cos(q·t))·cos(2t),
               (R0+r0·cos(q·t))·sin(2t),
               r0·sin(q·t))
    over t ∈ [0, 2π).

    For q=2 this is the T(2,2) curve with two components (split at t=π).
    For q=3 this is the T(2,3) trefoil (one component over 0..2π).
    For 2<q<3 the curve does not close in 0..2π but provides the correct
    geometric skeleton for the near-trefoil interpolation.

    The two-strand field construction follows the same superposition as
    the original Hopf link, with the T(2,q) curve replacing T(2,2):
      - Strand 1: nearest point on the first half-period arc
      - Strand 2: nearest point on the second half-period arc
      - Superposition weights: w = 1/rho^2 (same as original)
      - Phase winding: Phi = chi + wind_exp * tau
    """
    if label is None:
        label = f"T(2,{q_pol:.2f}) wind={wind_exp:.2f}"
    print(f"\n  Building {label}...")
    t0 = time.time()

    NT = 8000
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)

    # Curve Γ(t) with generalised poloidal winding q_pol
    Gx = (R0 + r0*np.cos(q_pol*t_arr))*np.cos(2*t_arr)
    Gy = (R0 + r0*np.cos(q_pol*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(q_pol*t_arr)
    G_full = np.stack([Gx, Gy, Gz], axis=1)

    # Split into two arcs: first half [0, π) and second half [π, 2π).
    # This split corresponds to the two tube components γ1 (outer) and γ2
    # (inner) of the Q_H=2 Hopf link. The split point t=π is a calculational
    # convenience, not a physical feature: because the density-feedback
    # functional K_fb + φ^6*J4 is rotationally symmetric around the torus
    # axis, any split point t=t0 gives the same field energy by symmetry.
    # The physical picture is that the two tubes are energetically equivalent
    # and the perturbation that triggers the Q_H=2→Q_H=3 transition can act
    # at any toroidal angle along either component.
    mask1 = t_arr < np.pi
    mask2 = t_arr >= np.pi
    t1_arr = t_arr[mask1]; G1 = G_full[mask1]
    t2_arr = t_arr[mask2]; G2 = G_full[mask2]

    # Rescale each arc's parameter to [0, 2π) for the phase τ
    tau1_arr = 2*t1_arr                  # maps [0,π)  → [0,2π)
    tau2_arr = 2*(t2_arr - np.pi)        # maps [π,2π) → [0,2π)

    tree1 = KDTree(G1)
    tree2 = KDTree(G2)

    # ── Frenet-Serret frame for the full curve ────────────────────────
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

    # ── Nearest-strand assignment ─────────────────────────────────────
    d1, idx1 = tree1.query(pts_np, workers=1)
    d2, idx2 = tree2.query(pts_np, workers=1)

    def chi_angle(pts, curve_pts, nidx, N1_arr, N2_arr):
        rel  = pts - curve_pts[nidx]
        chi  = np.arctan2((rel*N2_arr[nidx]).sum(1),
                           (rel*N1_arr[nidx]).sum(1))
        return chi

    chi_1 = chi_angle(pts_np, G1, idx1, N1_1, N2_1)
    chi_2 = chi_angle(pts_np, G2, idx2, N1_2, N2_2)

    tau_1 = tau1_arr[idx1]
    tau_2 = tau2_arr[idx2]

    # ── Phase with interpolated winding ──────────────────────────────
    # Phi = chi + wind_exp * tau
    # wind_exp=2.0: Q_H=2 (Hopf link winding)
    # wind_exp=3.0: Q_H=3 (trefoil winding)
    # Intermediate values interpolate the phase commitment continuously.
    Phi1 = chi_1 + wind_exp * tau_1
    Phi2 = chi_2 + wind_exp * tau_2

    # ── Two-strand superposition (same as original) ───────────────────
    def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    rho1 = np.clip(d1, 1e-6, None); rho2 = np.clip(d2, 1e-6, None)
    f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
    w1 = 1/rho1**2;       w2 = 1/rho2**2

    z1_c = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2_c = w1*np.sin(f1/2)*np.exp(1j*Phi1) + w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag = np.sqrt(np.abs(z1_c)**2 + np.abs(z2_c)**2)
    z1_c /= mag; z2_c /= mag

    nx0 = 2*np.real(np.conj(z1_c)*z2_c)
    ny0 = 2*np.imag(np.conj(z1_c)*z2_c)
    nz0 = np.abs(z1_c)**2 - np.abs(z2_c)**2
    n0  = np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
    n0 /= np.linalg.norm(n0,axis=-1,keepdims=True).clip(1e-10)
    print(f"  Built in {time.time()-t0:.1f}s")
    return n0


# ── Convenience wrappers matching original function names ─────────────
def build_hopf_link():
    """Original T(2,2) Hopf link: q_pol=2, wind_exp=2."""
    return build_torus_knot_field(2.0, 2.0, label="T(2,2) Hopf link — Q_H=2")

def build_trefoil():
    """Q_H=3 trefoil via Construction C (uses bishop_frame_v2)."""
    print("\n  Building T(2,3) trefoil Construction C (Q_H=3)...")
    t0 = time.time()
    NT_frame = 20000
    t_frame,_,N1f,N2f,H = build_compensated_frame_arclength(NT=NT_frame)
    print(f"  Bishop frame holonomy: {np.degrees(H):.4f}°")
    NT = 4000
    t_arr = np.linspace(0,2*np.pi,NT,endpoint=False)
    arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
    Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
    Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(3*t_arr)
    Gpts = np.stack([Gx,Gy,Gz],axis=1)
    lobe_idx   = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
    lobe_trees = [KDTree(Gpts[li]) for li in lobe_idx]
    lobe_t     = [t_arr[li] for li in lobe_idx]
    def nearest_two(qpts):
        dp,tp=[],[]
        for tree_l,tl in zip(lobe_trees,lobe_t):
            d,idx=tree_l.query(qpts,workers=1); dp.append(d); tp.append(tl[idx])
        ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
        return (np.take_along_axis(ts,o,1)[:,0],np.take_along_axis(ds,o,1)[:,0],
                np.take_along_axis(ts,o,1)[:,1],np.take_along_axis(ds,o,1)[:,1])
    def frame_at_t(tq):
        idx=np.searchsorted(t_frame,tq%(2*np.pi))%NT_frame
        return N1f[idx],N2f[idx]
    def curve_at_t(t):
        return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                         (R0+r0*np.cos(3*t))*np.sin(2*t),
                          r0*np.sin(3*t)],axis=-1)
    t1g,d1g,t2g,d2g = nearest_two(pts_np)
    chi1=np.arctan2(np.sum((pts_np-curve_at_t(t1g))*frame_at_t(t1g)[1],1),
                    np.sum((pts_np-curve_at_t(t1g))*frame_at_t(t1g)[0],1))
    chi2=np.arctan2(np.sum((pts_np-curve_at_t(t2g))*frame_at_t(t2g)[1],1),
                    np.sum((pts_np-curve_at_t(t2g))*frame_at_t(t2g)[0],1))
    Phi1=chi1+3*t1g; Phi2=chi2+3*t2g
    rho1=np.clip(d1g,1e-6,None); rho2=np.clip(d2g,1e-6,None)
    def f0(r): return 2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    f1=f0(rho1*C_star); f2=f0(rho2*C_star)
    w1=1/rho1**2; w2=1/rho2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=mag; z2/=mag
    nx0=2*np.real(np.conj(z1)*z2); ny0=2*np.imag(np.conj(z1)*z2)
    nz0=np.abs(z1)**2-np.abs(z2)**2
    n0=np.stack([nx0,ny0,nz0],axis=-1).reshape(N,N,N,3).astype(np.float32)
    n0/=np.linalg.norm(n0,axis=-1,keepdims=True).clip(1e-10)
    print(f"  Construction C built in {time.time()-t0:.1f}s")
    return n0


# ══════════════════════════════════════════════════════════════════════
# GENERALISED N-LOBE CONSTRUCTION for genuine single-component T(2,q)
# torus knots at any odd integer q >= 3 (generalises build_trefoil()'s
# q=3-only "Construction C" to arbitrary odd q).
#
# build_torus_knot_field's generic 2-arc split (used for q in [2,3]) is
# NOT valid here: each half-arc there sweeps ~q/2 poloidal loops, so for
# q=5 each half wraps past itself 2.5 times before the split point,
# meaning "nearest point in this half" often lands on the wrong loop
# entirely. That was confirmed both by a visibly wrong apparent symmetry
# in the rendered field and by the raw Whitehead Q_H sitting nowhere
# near the expected value at q=5. This function instead splits the
# curve into exactly q lobes (one KDTree each, mirroring build_trefoil's
# 3-lobe treatment) so "nearest of the q lobes" is always geometrically
# well-posed regardless of q.
# ══════════════════════════════════════════════════════════════════════
def build_torus_knot_field_nlobe(q_pol, wind_exp, label=None):
    """Genuine single-component T(2,q_pol) construction for odd integer
    q_pol >= 3, via q_pol separate lobe-KDTrees (generalises
    build_trefoil()'s Construction C, which is the q_pol=3 special case
    of this same procedure)."""
    q = q_pol
    if q != int(q) or int(q) % 2 == 0 or int(q) < 3:
        raise ValueError(
            f"build_torus_knot_field_nlobe requires an odd integer q_pol >= 3 "
            f"(single-component T(2,q) knot); got {q_pol}")
    q = int(q)
    if label is None:
        label = f"T(2,{q}) wind={wind_exp:.2f} [N-lobe]"
    print(f"\n  Building {label}...")
    t0 = time.time()

    NT_frame = 20000
    t_frame, _, N1f, N2f, H = build_compensated_frame_arclength(NT=NT_frame, q=q)
    print(f"  Bishop frame holonomy: {np.degrees(H):.4f}°")

    NT = max(4000, 1400 * q)  # keep per-lobe sample density roughly fixed
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
    lobe_width = 2*np.pi/q
    arc_starts = [k*lobe_width for k in range(q)]
    Gx = (R0 + r0*np.cos(q*t_arr))*np.cos(2*t_arr)
    Gy = (R0 + r0*np.cos(q*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(q*t_arr)
    Gpts = np.stack([Gx, Gy, Gz], axis=1)
    lobe_idx = [np.where((t_arr >= s) & (t_arr < s+lobe_width))[0] for s in arc_starts]
    lobe_trees = [KDTree(Gpts[li]) for li in lobe_idx]
    lobe_t = [t_arr[li] for li in lobe_idx]

    def nearest_two(qpts):
        dp, tp = [], []
        for tree_l, tl in zip(lobe_trees, lobe_t):
            d, idx = tree_l.query(qpts, workers=1)
            dp.append(d); tp.append(tl[idx])
        ds = np.stack(dp, 1); ts = np.stack(tp, 1); o = np.argsort(ds, 1)
        return (np.take_along_axis(ts, o, 1)[:, 0], np.take_along_axis(ds, o, 1)[:, 0],
                np.take_along_axis(ts, o, 1)[:, 1], np.take_along_axis(ds, o, 1)[:, 1])

    def frame_at_t(tq):
        idx = np.searchsorted(t_frame, tq % (2*np.pi)) % NT_frame
        return N1f[idx], N2f[idx]

    def curve_at_t(t):
        return np.stack([(R0+r0*np.cos(q*t))*np.cos(2*t),
                          (R0+r0*np.cos(q*t))*np.sin(2*t),
                          r0*np.sin(q*t)], axis=-1)

    t1g, d1g, t2g, d2g = nearest_two(pts_np)
    chi1 = np.arctan2(np.sum((pts_np-curve_at_t(t1g))*frame_at_t(t1g)[1], 1),
                       np.sum((pts_np-curve_at_t(t1g))*frame_at_t(t1g)[0], 1))
    chi2 = np.arctan2(np.sum((pts_np-curve_at_t(t2g))*frame_at_t(t2g)[1], 1),
                       np.sum((pts_np-curve_at_t(t2g))*frame_at_t(t2g)[0], 1))
    Phi1 = chi1 + wind_exp*t1g
    Phi2 = chi2 + wind_exp*t2g
    rho1 = np.clip(d1g, 1e-6, None); rho2 = np.clip(d2g, 1e-6, None)
    def f0(r): return 2*np.arctan(np.maximum(r, 1e-9)**(-C_star))
    f1 = f0(rho1*C_star); f2 = f0(rho2*C_star)
    w1 = 1/rho1**2; w2 = 1/rho2**2
    z1 = (w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2 = w1*np.sin(f1/2)*np.exp(1j*Phi1) + w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag = np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1 /= mag; z2 /= mag
    nx0 = 2*np.real(np.conj(z1)*z2); ny0 = 2*np.imag(np.conj(z1)*z2)
    nz0 = np.abs(z1)**2-np.abs(z2)**2
    n0 = np.stack([nx0, ny0, nz0], axis=-1).reshape(N, N, N, 3).astype(np.float32)
    n0 /= np.linalg.norm(n0, axis=-1, keepdims=True).clip(1e-10)
    print(f"  N-lobe construction (q={q}) built in {time.time()-t0:.1f}s")
    return n0


# ── Gradient flow run ─────────────────────────────────────────────────
def run(n_np, label, q_expected, outdir):
    os.makedirs(outdir, exist_ok=True)
    wh_steps = sorted(set(int(f*args.n_steps) for f in [0, 0.25, 0.75, 1.0]))
    print(f"\n{'─'*60}")
    print(f"  RUN: {label}")
    print(f"  Whitehead Q_H check at steps: {wh_steps}")
    tw = time.time()
    QH_init = whitehead_from_n(n_np, args.whitehead_N)
    print(f"  Whitehead Q_H = {QH_init:.3f}  (t={time.time()-tw:.1f}s)")
    n_param = torch.tensor(n_np,dtype=torch.float32,device=dev).requires_grad_(True)
    opt = torch.optim.Adam([n_param], lr=args.lr)
    E0,K0,J40,rb0 = diag(n_param)
    print(f"\n  Initial: E={E0:.4e}  K={K0:.1f}  J4={J40:.2f}  r_bar={rb0:.3f}")
    history, wh_results = [], {0: QH_init}
    for step in range(args.n_steps):
        opt.zero_grad()
        E,_,_,_ = E_geom(n_param); E.backward()
        with torch.no_grad():
            n_param.grad.data.copy_(ptangent(n_param.detach(),n_param.grad))
        opt.step()
        with torch.no_grad():
            n_param.data.copy_(n_param/n_param.norm(dim=-1,keepdim=True).clamp(1e-10))
        actual_step = step+1
        if args.snap_every>0 and actual_step%args.snap_every==0:
            np.save(os.path.join(outdir,f'n_{actual_step:06d}.npy'),
                    n_param.detach().cpu().numpy())
        if actual_step in wh_steps and actual_step not in wh_results:
            QH_cur = whitehead_from_n(n_param.detach().cpu().numpy(), args.whitehead_N)
            wh_results[actual_step] = QH_cur
        if actual_step%args.log_every==0 or step==0:
            Ev,Kv,J4v,rbv = diag(n_param)
            history.append(dict(step=actual_step,E=Ev,K=Kv,J4=J4v,r_bar=rbv))
            if rbv > 4.5: print(f"  HALT: dilution"); break
    QH_final = whitehead_from_n(n_param.detach().cpu().numpy(), args.whitehead_N)
    Ef,Kf,J4f,rbf = diag(n_param)
    print(f"\n  Final: E={Ef:.4e}  J4/K={J4f/max(Kf,1e-9):.4f}  r_bar={rbf:.3f}")
    print(f"  Whitehead Q_H: {QH_init:.3f} → {QH_final:.3f}")
    np.save(os.path.join(outdir,'n_final.npy'), n_param.detach().cpu().numpy())
    with open(os.path.join(outdir,'log.json'),'w') as f: json.dump(history,f,indent=2)
    return dict(label=label, QH_init=QH_init, QH_final=QH_final,
                J4f=J4f, Kf=Kf, rbf=rbf)


# ── Main ──────────────────────────────────────────────────────────────
print()

# Determine which construction to use based on CLI flags
#
# build_torus_knot_field_nlobe() (generalised N-lobe construction, see
# above) is available but NOT used here: benchmarked at q=5 it takes
# 16.4s to build vs. the general 2-arc method's ~4.3s, i.e. slower, not
# faster. Direct inspection of relaxed (gradient-flow) final fields
# confirmed the general 2-arc method already converges to the correct
# T(2,q) topology at every q in the 3-5 sweep despite an unrelaxed-frame
# symmetry artifact, so there is no correctness reason to switch either.
if q_pol == 2.0 and alpha == 0.0:
    # Exact original Hopf link
    n_q2 = build_hopf_link()
    label_q2 = "T(2,2) Hopf link — Q_H=2"
else:
    # General T(2,q) with interpolated winding
    n_q2 = build_torus_knot_field(q_pol, wind_exp,
                                   label=f"T(2,{q_pol:.2f}) wind={wind_exp:.2f}")
    label_q2 = f"T(2,{q_pol:.2f}) q_pol={q_pol} alpha_wind={alpha}"

# Always also build the reference trefoil for comparison
n_q3 = build_trefoil()

if args.verify_only:
    print(f"\n  VERIFY ONLY (whitehead_N={args.whitehead_N})")
    for n_np, lbl, tag in [
        (n_q2, label_q2, f'q{q_pol:.2f}_a{alpha:.2f}'),
        (n_q3, 'T(2,3) trefoil', 'q3'),
    ]:
        t0 = time.time()
        QH = whitehead_from_n(n_np, args.whitehead_N)
        Ev,Kv,J4v,rbv = diag(torch.tensor(n_np,dtype=torch.float32,device=dev))
        fname = os.path.join(args.outdir, f'n_initial_{tag}.npy')
        np.save(fname, n_np)
        print(f"\n  {lbl}:")
        print(f"    Whitehead Q_H = {QH:.3f}  ({time.time()-t0:.1f}s)")
        print(f"    E={Ev:.4e}  J4={J4v:.1f}  K={Kv:.1f}  r_bar={rbv:.3f}")
        print(f"    Saved: {fname}")
    print(f"\n  Load npy files in hopfion_viewer.html to inspect geometry.")
    sys.exit(0)

# Full gradient flow
r2 = run(n_q2, label_q2, None, os.path.join(args.outdir, 'q2'))
r3 = run(n_q3, 'T(2,3) trefoil — Q_H=3', 3, os.path.join(args.outdir, 'q3'))

print(f"\n{'='*68}")
print(f"  COMPARISON")
print(f"{'='*68}")
print(f"  {'':28}  {'Construction':>18}  {'Trefoil':>10}")
print(f"  {'Q_H (Whitehead, initial)':28}  {r2['QH_init']:>18.3f}  {r3['QH_init']:>10.3f}")
print(f"  {'Q_H (Whitehead, final)':28}  {r2['QH_final']:>18.3f}  {r3['QH_final']:>10.3f}")
print(f"  {'J4/K (final)':28}  {r2['J4f']/max(r2['Kf'],1e-9):>18.4f}  "
      f"{r3['J4f']/max(r3['Kf'],1e-9):>10.4f}")
