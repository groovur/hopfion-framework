"""
whitehead_perstrand_charge.py
==============================
Compute the Hopf charge Q_H of the CORRECTED per-strand S3-lift
ansatz (Phi_i = chi_i + 3*t_i, using the arc-length-uniform
compensated Bishop frame and the robust 3-lobe strand search) via
the Whitehead (Berry-Chern-Simons) integral:

    Q_H = (1 / 4*pi^2) * int A.(curl A) d^3x

where A_i = Im(Z-dagger d_i Z).

Validation: the same inverse-stereographic Q_H=1 reference doublet
used in whitehead_hopf_charge.py is included for comparison at the
same grid resolution.

Usage:
    python whitehead_perstrand_charge.py            # coarse: N=40
    python whitehead_perstrand_charge.py --N 80      # medium
    python whitehead_perstrand_charge.py --N 120     # fine

F. Manfredi / verification, June 2026
"""
import sys, time
import numpy as np
from scipy.spatial import KDTree
from bishop_frame_v2 import build_compensated_frame_arclength

N = 40
for i, arg in enumerate(sys.argv[1:]):
    if arg == '--N' and i+1 < len(sys.argv)-1:
        N = int(sys.argv[i+2])

print("="*68)
print(f"WHITEHEAD INTEGRAL FOR PER-STRAND (Phi=chi+3t) S3-LIFT  (N={N})")
print("="*68)

R0, r0, C_star = 3.0, 0.874, 2.5062

t0 = time.time()
NT_frame = 30000
t_frame, T_frame, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=NT_frame)
print(f"Compensated frame built, holonomy H={np.degrees(H):.4f} deg ({time.time()-t0:.1f}s)")

NT = 6000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
lobe_indices = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees = [KDTree(Gamma_pts[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr[li] for li in lobe_indices]

def f0(rho_hat): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]

def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

def nearest_two_strands_robust(pts):
    d_per_lobe, t_per_lobe = [], []
    for tree_l, t_l in zip(lobe_trees, lobe_t_arrays):
        d, idx = tree_l.query(pts)
        d_per_lobe.append(d); t_per_lobe.append(t_l[idx])
    d_stack = np.stack(d_per_lobe, axis=1)
    t_stack = np.stack(t_per_lobe, axis=1)
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

def Z_perstrand(pts, C_star_unused=None):
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

def Z_q1_reference(pts, C_star_unused=None):
    x,y,z = pts[:,0],pts[:,1],pts[:,2]
    rr2 = x**2+y**2+z**2
    denom = 1+rr2
    Z1 = ((1-rr2) + 2j*z)/denom
    Z2 = 2*(x+1j*y)/denom
    return Z1, Z2

def whitehead(Zfunc, N, L):
    h = 2*L/(N-1)
    ax = np.linspace(-L,L,N)
    X,Y,Z = np.meshgrid(ax,ax,ax,indexing='ij')
    pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    z1,z2 = Zfunc(pts)
    z1=z1.reshape(N,N,N); z2=z2.reshape(N,N,N)
    dz1=[np.gradient(z1,h,axis=a) for a in range(3)]
    dz2=[np.gradient(z2,h,axis=a) for a in range(3)]
    A=[np.imag(np.conj(z1)*dz1[a]+np.conj(z2)*dz2[a]) for a in range(3)]
    curlA=[np.gradient(A[2],h,axis=1)-np.gradient(A[1],h,axis=2),
           np.gradient(A[0],h,axis=2)-np.gradient(A[2],h,axis=0),
           np.gradient(A[1],h,axis=0)-np.gradient(A[0],h,axis=1)]
    integrand = sum(A[i]*curlA[i] for i in range(3))
    return h**3*np.sum(integrand)/(4*np.pi**2)

print(f"\n--- VALIDATION: Q_H=1 inverse-stereographic Hopf doublet ---")
t1 = time.time()
qh_ref = whitehead(Z_q1_reference, N, 6.0)
print(f"  Q_H (reference) = {qh_ref:.5f}  (expected 1.0; t={time.time()-t1:.1f}s)")

print(f"\n--- MAIN: per-strand Phi=chi+3t, S3-lift ansatz ---")
t2 = time.time()
qh_main = whitehead(Z_perstrand, N, 7.0)
print(f"  Q_H (per-strand S3-lift) = {qh_main:.5f}  (t={time.time()-t2:.1f}s)")
print(f"\nTotal time: {time.time()-t0:.1f}s")
