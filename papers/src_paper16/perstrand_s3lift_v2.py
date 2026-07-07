"""
Corrected per-strand Phi_i = chi_i + 3*t_i ansatz, S3-lift blended,
using:
  (a) the arc-length-uniform compensated Bishop frame (bishop_frame_v2)
  (b) the ROBUST 3-lobe-partition strand search (fixes the k=80
      KDTree silent-failure bug found and confirmed this session)
"""
import numpy as np
from scipy.spatial import KDTree
from bishop_frame_v2 import build_compensated_frame_arclength, Gamma

R0, r0, C_star = 3.0, 0.874, 2.5062

NT_frame = 30000
t_frame, T_frame, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=NT_frame)

# Robust lobe partition (validated earlier this session for the C* sweep)
NT = 6000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
lobe_idx_ranges = [(s, s+2*np.pi/3) for s in arc_starts]

def lobe_of_t(t):
    """Which lobe (0,1,2) does parameter t belong to."""
    tm = t % (2*np.pi)
    return int(tm // (2*np.pi/3)) % 3

def f0(rho_hat): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

# build per-lobe KD-trees (robust strand search)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
lobe_indices = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees = [KDTree(Gamma_pts[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr[li] for li in lobe_indices]

def nearest_two_strands_robust(pts):
    """Find distance+parameter to nearest point in EACH of the 3 lobes,
    then return the 2 smallest (robust: always finds genuinely
    distinct strands, no silent k-too-small failure)."""
    d_per_lobe = []
    t_per_lobe = []
    for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
        d, idx = tree_l.query(pts)
        d_per_lobe.append(d)
        t_per_lobe.append(t_l[idx])
    d_stack = np.stack(d_per_lobe, axis=1)   # (npts,3)
    t_stack = np.stack(t_per_lobe, axis=1)   # (npts,3)
    order = np.argsort(d_stack, axis=1)
    d_sorted = np.take_along_axis(d_stack, order, axis=1)
    t_sorted = np.take_along_axis(t_stack, order, axis=1)
    return t_sorted[:,0], d_sorted[:,0], t_sorted[:,1], d_sorted[:,1]

def chi_for_t(pts, t_query):
    N1_pts, N2_pts = frame_at_t(t_query)
    curve_pts = curve_at_t(t_query)
    rel = pts - curve_pts
    chi = np.arctan2(np.sum(rel*N2_pts,axis=1), np.sum(rel*N1_pts,axis=1))
    dist = np.linalg.norm(rel, axis=1)
    return chi, dist

def Z_perstrand_s3lift_v2(pts):
    t1, d1, t2, d2 = nearest_two_strands_robust(pts)
    chi1, rho1 = chi_for_t(pts, t1)
    chi2, rho2 = chi_for_t(pts, t2)
    f1 = f0(np.clip(rho1,1e-6,None)*C_star)
    f2 = f0(np.clip(rho2,1e-6,None)*C_star)
    Phi1 = chi1 + 3*t1
    Phi2 = chi2 + 3*t2
    z1a = np.cos(f1/2).astype(complex); z2a = np.sin(f1/2)*np.exp(1j*Phi1)
    z1b = np.cos(f2/2).astype(complex); z2b = np.sin(f2/2)*np.exp(1j*Phi2)
    w1 = 1.0/np.clip(rho1,1e-6,None)**2
    w2 = 1.0/np.clip(rho2,1e-6,None)**2
    z1 = w1*z1a + w2*z1b
    z2 = w1*z2a + w2*z2b
    norm = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    return z1/norm, z2/norm

if __name__ == "__main__":
    # Sanity: robust search at the same s=-0.05 point that exposed the bug
    t1c = np.pi/6
    p_over = curve_at_t(np.array([t1c]))[0]
    p_under = curve_at_t(np.array([t1c+np.pi]))[0]
    midpoint = (p_over+p_under)/2
    direction = (p_under-p_over)/np.linalg.norm(p_under-p_over)
    pt = (midpoint - 0.05*direction).reshape(1,3)
    t1u,d1,t2u,d2 = nearest_two_strands_robust(pt)
    print(f"At s=-0.05 (the point that exposed the old bug):")
    print(f"  t1={t1u[0]:.4f}, d1={d1[0]:.6f}, t2={t2u[0]:.4f}, d2={d2[0]:.6f}")
    print(f"  (BUG CHECK: d2 should be a real, finite, sensible distance, NOT 1e6)")
