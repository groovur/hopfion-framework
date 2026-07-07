"""
Apply the complex-superposition (S3-lift) treatment to the corrected
per-strand doublets Z_i = (cos(f_i/2), sin(f_i/2)*exp(i*(chi_i+3t_i))),
and test smoothness across the crossing midline.
"""
import numpy as np
from scipy.spatial import KDTree
from bishop_frame_compensated import build_compensated_frame

R0, r0, C_star = 3.0, 0.874, 2.5062

NT_frame = 30000
t_frame, T_frame, N1_frame, N2_frame, H = build_compensated_frame(NT=NT_frame)

NT = 8000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
tree = KDTree(Gamma_pts)

def f0(rho_hat): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

def chi_for_t(pts, t_query):
    """chi of pts relative to the curve point and frame AT a SPECIFIED
    t_query (not necessarily the nearest point) -- needed so we can
    evaluate strand 'over' and strand 'under' doublets at the SAME
    spatial point x using EACH strand's own (possibly non-nearest) t."""
    N1_pts, N2_pts = frame_at_t(t_query)
    curve_pts = curve_at_t(t_query)
    rel = pts - curve_pts
    chi = np.arctan2(np.sum(rel*N2_pts,axis=1), np.sum(rel*N1_pts,axis=1))
    dist = np.linalg.norm(rel, axis=1)
    return chi, dist

def nearest_two_t(pts, k=80, min_sep=np.pi/2):
    dists, idxs = tree.query(pts, k=k)
    t1 = t_arr[idxs[:,0]]; d1 = dists[:,0]
    dt = np.abs(((t_arr[idxs]-t1[:,None]+np.pi)%(2*np.pi))-np.pi)
    mask = dt>min_sep; fi = np.argmax(mask,axis=1); has2=mask.any(axis=1)
    d2 = np.where(has2, dists[np.arange(len(pts)),fi], 1e6)
    t2 = t_arr[idxs[np.arange(len(pts)),fi]]
    return t1,d1,t2,d2

def Z_perstrand_s3lift(pts):
    t1, d1, t2, d2 = nearest_two_t(pts)
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
    t1c = np.pi/6
    p_over = curve_at_t(np.array([t1c]))[0]
    p_under = curve_at_t(np.array([t1c+np.pi]))[0]
    midpoint = (p_over+p_under)/2
    direction = (p_under-p_over)/np.linalg.norm(p_under-p_over)

    print("="*70)
    print("SMOOTHNESS CHECK: PER-STRAND PHI=CHI+3T DOUBLETS, S3-LIFT BLENDED")
    print("="*70)
    print(f"{'s':>8}  {'n_x':>10}  {'n_y':>10}  {'n_z':>10}")
    for s in np.linspace(-0.3,0.3,13):
        pt = (midpoint + s*direction).reshape(1,3)
        z1,z2 = Z_perstrand_s3lift(pt)
        z1,z2 = z1[0],z2[0]
        nx=2*np.real(np.conj(z1)*z2); ny=2*np.imag(np.conj(z1)*z2); nz=np.abs(z1)**2-np.abs(z2)**2
        print(f"{s:>8.4f}  {nx:>10.5f}  {ny:>10.5f}  {nz:>10.5f}")
