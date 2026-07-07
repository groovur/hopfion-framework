"""
Check whether the near-tube-core J4 excess of Construction C shrinks
(integrable, resolution artifact) or stays/grows (genuine new
divergence) under grid refinement, restricted to a small box around
ONE crossing region (matching the methodology already validated for
Theorem thm:s3_bounded_energy / Proposition prop:perstrand_smooth).
"""
import numpy as np
from scipy.spatial import KDTree
import sys, time
sys.path.insert(0, '.')
from bishop_frame_v2 import build_compensated_frame_arclength

R0, r0, C_star = 3.0, 0.874, 2.5062  # use the SAME C* as the rest of this
                                       # session's local checks (2.5062, not
                                       # the solver default 2.4987 -- close
                                       # enough that conclusions transfer)

NT_frame = 30000
t_frame, T_frame, N1_frame, N2_frame, H = build_compensated_frame_arclength(NT=NT_frame)
def frame_at_t(t_query):
    idx = np.searchsorted(t_frame, t_query % (2*np.pi)) % NT_frame
    return N1_frame[idx], N2_frame[idx]
def curve_at_t(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                      (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], axis=-1)

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

def build_n_C(qpts):
    t1, d1, t2, d2 = nearest_two_strands_robust(qpts)
    chi1, rho1 = chi_for_t(qpts, t1)
    chi2, rho2 = chi_for_t(qpts, t2)
    f1 = f0(np.clip(rho1,1e-6,None)*C_star); f2 = f0(np.clip(rho2,1e-6,None)*C_star)
    Phi1 = chi1 + 3*t1; Phi2 = chi2 + 3*t2
    z1a = np.cos(f1/2).astype(complex); z2a = np.sin(f1/2)*np.exp(1j*Phi1)
    z1b = np.cos(f2/2).astype(complex); z2b = np.sin(f2/2)*np.exp(1j*Phi2)
    w1 = 1.0/np.clip(rho1,1e-6,None)**2; w2 = 1.0/np.clip(rho2,1e-6,None)**2
    z1 = w1*z1a + w2*z1b; z2 = w1*z2a + w2*z2b
    norm = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    z1,z2 = z1/norm, z2/norm
    nx = 2*np.real(np.conj(z1)*z2); ny = 2*np.imag(np.conj(z1)*z2); nz = np.abs(z1)**2-np.abs(z2)**2
    n = np.stack([nx,ny,nz],axis=-1)
    return n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)

def rho_J4_density(n, h):
    nx,ny,nz = n[...,0], n[...,1], n[...,2]
    def cd(u,a): return (np.roll(u,-1,axis=a)-np.roll(u,1,axis=a))/(2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    return Fxy**2+Fxz**2+Fyz**2

t1c = np.pi/6
p_over = curve_at_t(np.array([t1c]))[0]
p_under = curve_at_t(np.array([t1c+np.pi]))[0]
midpoint = (p_over+p_under)/2

L = 1.5*r0
print("Refinement scan: J4 near a single crossing, Construction C")
print(f"{'Ngrid':>6}  {'h':>9}  {'J4_local':>12}  {'maxJ4dens':>12}")
for Ngrid in [41, 61, 81, 121, 161]:
    h_ = 2*L/(Ngrid-1)
    ax = np.linspace(-L,L,Ngrid)
    X,Y,Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    qpts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    n = build_n_C(qpts).reshape(Ngrid,Ngrid,Ngrid,3)
    rho_J4 = rho_J4_density(n, h_)
    margin = max(5,Ngrid//16)
    core = slice(margin,Ngrid-margin)
    rj_core = rho_J4[core,core,core]
    J4_local = rj_core.sum()*h_**3
    print(f"{Ngrid:>6}  {h_:>9.5f}  {J4_local:>12.4f}  {rj_core.max():>12.2f}")
