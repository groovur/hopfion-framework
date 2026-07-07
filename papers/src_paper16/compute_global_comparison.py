"""
Evaluate the SOLVER's actual energy functional (compute_global, ported
verbatim from qh3_trefoil_solver_3d_v11.py) on three constructions:
  (A) the original, chi-free, single-strand ansatz (build_ic, the
      solver's own existing initial condition -- has the proven
      pi-discontinuity at crossings)
  (B) the chi-free GLOBAL S3-LIFT ansatz (smooth, bounded energy
      locally, but Q_H=0 -- Proposition prop:s3lift_charge_zero)
  (C) the corrected PER-STRAND Phi=chi+3t construction (smooth,
      Q_H=3 -- Theorem thm:perstrand_charge_three)

at MATCHED grid settings (the solver's own defaults: N=64, h=0.26,
C*=2.4987, R0=3.0, r0=0.874), to directly compare total K, J2a, J4,
r_bar, and the actual minimised objective E_geom = K*J4.
"""
import numpy as np
from scipy.spatial import KDTree
import sys, time
sys.path.insert(0, '.')
from bishop_frame_v2 import build_compensated_frame_arclength

# ── Solver's own defaults ────────────────────────────────────────────
N = 64
h = 0.26
C_star = 2.4987
R0, r0 = 3.0, 0.874
NT = 1500
MU = 3.0 - (1+5**0.5)/2   # solver's own mu_param = 3 - phi

cv = h*(np.arange(N) - N//2 + 0.5)
pts = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), axis=-1).reshape(-1,3).astype(np.float32)
dist_from_origin = np.linalg.norm(pts, axis=-1).reshape(N,N,N).astype(np.float32)

print(f"Grid: N={N}, h={h}, box=[{-N*h/2:.2f},{N*h/2:.2f}], C*={C_star}")
print(f"{pts.shape[0]} points total")

# ── Port compute_global to pure numpy (no torch needed for this test) ──
def compute_global_np(n):
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = np.clip(1 - nz**2, 0, 1)**2
    def cd(u, a): return (np.roll(u,-1,axis=a) - np.roll(u,1,axis=a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a  = float((s4 * g2).sum() * h**3)
    J2iso = float(g2.sum() * h**3)
    K = J2a + MU * J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4 = float(rho_J4.sum() * h**3)
    r_bar_num = float((rho_J4 * dist_from_origin).sum())
    r_bar_den = float(rho_J4.sum())
    r_bar = r_bar_num/r_bar_den if r_bar_den > 1e-12 else float('nan')
    return K, J2a, J2iso, J4, r_bar

# ── Construction A: original chi-free single-strand ansatz ──────────
def build_ansatz_A():
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
    Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
    Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(3*t_arr)
    Gamma = np.stack([Gx,Gy,Gz], axis=1).astype(np.float32)
    tree = KDTree(Gamma)
    dists, idx = tree.query(pts)
    rho = np.clip(dists, 1e-6, None).astype(np.float32)
    t_near = t_arr[idx].astype(np.float32)
    f_ = 2*np.arctan(rho**(-C_star))
    tht = 3*t_near
    nx_ = np.sin(f_)*np.cos(tht); ny_ = np.sin(f_)*np.sin(tht); nz_ = np.cos(f_)
    n_np = np.stack([nx_.reshape(N,N,N), ny_.reshape(N,N,N), nz_.reshape(N,N,N)], axis=-1).astype(np.float32)
    n_np /= np.linalg.norm(n_np, axis=-1, keepdims=True).clip(1e-10)
    return n_np

print("\n--- Building and evaluating Construction A (original, chi-free) ---")
t0=time.time()
nA = build_ansatz_A()
K,J2a,J2iso,J4,r_bar = compute_global_np(nA)
print(f"  K={K:.4f}  J2a={J2a:.4f}  J2iso={J2iso:.4f}  J4={J4:.4f}  r_bar={r_bar:.4f}  E_geom=K*J4={K*J4:.4f}  ({time.time()-t0:.1f}s)")

# ── Construction B: chi-free global S3-lift (Q_H=0, but smooth/bounded locally) ──
NT_b = 6000
t_arr_b = np.linspace(0, 2*np.pi, NT_b, endpoint=False)
Gx_b = (R0+r0*np.cos(3*t_arr_b))*np.cos(2*t_arr_b)
Gy_b = (R0+r0*np.cos(3*t_arr_b))*np.sin(2*t_arr_b)
Gz_b = r0*np.sin(3*t_arr_b)
Gamma_pts_b = np.stack([Gx_b,Gy_b,Gz_b], axis=1)
tree_b = KDTree(Gamma_pts_b)

def f0(rho_hat): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

def nearest_two_b(qpts, k=80, min_sep=np.pi/2):
    dists, idxs = tree_b.query(qpts, k=k)
    t1 = t_arr_b[idxs[:,0]]; d1 = dists[:,0]
    dt = np.abs(((t_arr_b[idxs]-t1[:,None]+np.pi)%(2*np.pi))-np.pi)
    mask = dt>min_sep; fi = np.argmax(mask,axis=1); has2=mask.any(axis=1)
    d2 = np.where(has2, dists[np.arange(len(qpts)),fi], 1e6)
    t2 = t_arr_b[idxs[np.arange(len(qpts)),fi]]
    return t1,d1,t2,d2

def build_ansatz_B():
    t1,d1,t2,d2 = nearest_two_b(pts)
    f1 = f0(np.clip(d1,1e-6,None)*C_star); f2 = f0(np.clip(d2,1e-6,None)*C_star)
    w1 = 1.0/np.clip(d1,1e-6,None)**2; w2 = 1.0/np.clip(d2,1e-6,None)**2
    z1u = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2u = w1*np.sin(f1/2)*np.exp(3j*t1) + w2*np.sin(f2/2)*np.exp(3j*t2)
    mag = np.sqrt(np.abs(z1u)**2+np.abs(z2u)**2)
    z1,z2 = z1u/mag, z2u/mag
    nx = 2*np.real(np.conj(z1)*z2); ny = 2*np.imag(np.conj(z1)*z2); nz = np.abs(z1)**2-np.abs(z2)**2
    n_np = np.stack([nx.reshape(N,N,N), ny.reshape(N,N,N), nz.reshape(N,N,N)], axis=-1).astype(np.float32)
    n_np /= np.linalg.norm(n_np, axis=-1, keepdims=True).clip(1e-10)
    return n_np

print("\n--- Building and evaluating Construction B (chi-free S3-lift, Q_H=0) ---")
t0=time.time()
nB = build_ansatz_B()
K,J2a,J2iso,J4,r_bar = compute_global_np(nB)
print(f"  K={K:.4f}  J2a={J2a:.4f}  J2iso={J2iso:.4f}  J4={J4:.4f}  r_bar={r_bar:.4f}  E_geom=K*J4={K*J4:.4f}  ({time.time()-t0:.1f}s)")

# ── Construction C: corrected per-strand Phi=chi+3t S3-lift (Q_H=3) ──
NT_frame = 30000
t_frame, T_frame, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=NT_frame)

NT_c = 6000
t_arr_c = np.linspace(0, 2*np.pi, NT_c, endpoint=False)
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
Gx_c = (R0+r0*np.cos(3*t_arr_c))*np.cos(2*t_arr_c)
Gy_c = (R0+r0*np.cos(3*t_arr_c))*np.sin(2*t_arr_c)
Gz_c = r0*np.sin(3*t_arr_c)
Gamma_pts_c = np.stack([Gx_c,Gy_c,Gz_c], axis=1)
lobe_indices = [np.where((t_arr_c>=s)&(t_arr_c<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees = [KDTree(Gamma_pts_c[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr_c[li] for li in lobe_indices]

def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

def nearest_two_strands_robust(qpts):
    d_per_lobe, t_per_lobe = [], []
    for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
        d, idx = tree_l.query(qpts)
        d_per_lobe.append(d); t_per_lobe.append(t_l[idx])
    d_stack = np.stack(d_per_lobe, axis=1)
    t_stack = np.stack(t_per_lobe, axis=1)
    order = np.argsort(d_stack, axis=1)
    d_sorted = np.take_along_axis(d_stack, order, axis=1)
    t_sorted = np.take_along_axis(t_stack, order, axis=1)
    return t_sorted[:,0], d_sorted[:,0], t_sorted[:,1], d_sorted[:,1]

def chi_for_t(qpts, t_query):
    N1_pts, N2_pts = frame_at_t(t_query)
    curve_pts = curve_at_t(t_query)
    rel = qpts - curve_pts
    chi = np.arctan2(np.sum(rel*N2_pts,axis=1), np.sum(rel*N1_pts,axis=1))
    dist = np.linalg.norm(rel, axis=1)
    return chi, dist

def build_ansatz_C():
    t1, d1, t2, d2 = nearest_two_strands_robust(pts)
    chi1, rho1 = chi_for_t(pts, t1)
    chi2, rho2 = chi_for_t(pts, t2)
    f1 = f0(np.clip(rho1,1e-6,None)*C_star); f2 = f0(np.clip(rho2,1e-6,None)*C_star)
    Phi1 = chi1 + 3*t1; Phi2 = chi2 + 3*t2
    z1a = np.cos(f1/2).astype(complex); z2a = np.sin(f1/2)*np.exp(1j*Phi1)
    z1b = np.cos(f2/2).astype(complex); z2b = np.sin(f2/2)*np.exp(1j*Phi2)
    w1 = 1.0/np.clip(rho1,1e-6,None)**2; w2 = 1.0/np.clip(rho2,1e-6,None)**2
    z1 = w1*z1a + w2*z1b; z2 = w1*z2a + w2*z2b
    norm = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    z1,z2 = z1/norm, z2/norm
    nx = 2*np.real(np.conj(z1)*z2); ny = 2*np.imag(np.conj(z1)*z2); nz = np.abs(z1)**2-np.abs(z2)**2
    n_np = np.stack([nx.reshape(N,N,N), ny.reshape(N,N,N), nz.reshape(N,N,N)], axis=-1).astype(np.float32)
    n_np /= np.linalg.norm(n_np, axis=-1, keepdims=True).clip(1e-10)
    return n_np

print("\n--- Building and evaluating Construction C (corrected per-strand, Q_H=3) ---")
t0=time.time()
nC = build_ansatz_C()
K,J2a,J2iso,J4,r_bar = compute_global_np(nC)
print(f"  K={K:.4f}  J2a={J2a:.4f}  J2iso={J2iso:.4f}  J4={J4:.4f}  r_bar={r_bar:.4f}  E_geom=K*J4={K*J4:.4f}  ({time.time()-t0:.1f}s)")
