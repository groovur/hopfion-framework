"""
Test the prediction: a genuine THREE-strand complex superposition
(rather than two) should have a SMALLER Q_H deviation from zero near
C*~0.3-0.4, since the previously-omitted third lobe is now properly
included with its own correct weight, rather than acting as an
unaccounted-for source of asymmetry/circulation.
"""
import numpy as np
from scipy.spatial import KDTree
import time

R0, r0 = 3.0, 0.874
NT = 6000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
tree = KDTree(np.stack([Gx,Gy,Gz], axis=1))

def f0(rho_hat, C_star): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

def nearest_three(pts, k=400, min_sep=np.pi/2):
    dists, idxs = tree.query(pts, k=k)
    t_all = t_arr[idxs]
    results_d, results_t = [], []
    for strand_num in range(3):
        if strand_num == 0:
            d_this = dists[:,0]; t_this = t_all[:,0]
        else:
            excl = np.zeros((pts.shape[0], k), dtype=bool)
            for prev_t in results_t:
                dt = np.abs(((t_all - prev_t[:,None] + np.pi)%(2*np.pi))-np.pi)
                excl |= (dt < min_sep)
            valid = ~excl
            fi = np.argmax(valid, axis=1)
            has_valid = valid.any(axis=1)
            d_this = np.where(has_valid, dists[np.arange(pts.shape[0]), fi], 1e6)
            t_this = t_all[np.arange(pts.shape[0]), fi]
        results_d.append(d_this); results_t.append(t_this)
    return results_d, results_t

def Z_s3lift_3strand(pts, C_star):
    ds, ts = nearest_three(pts)
    z1 = np.zeros(pts.shape[0], dtype=complex)
    z2 = np.zeros(pts.shape[0], dtype=complex)
    for d,t in zip(ds, ts):
        f = f0(np.clip(d,1e-6,None)*C_star, C_star)
        w = 1.0/np.clip(d,1e-6,None)**2
        z1 += w*np.cos(f/2)
        z2 += w*np.sin(f/2)*np.exp(3j*t)
    mag = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    return z1/mag, z2/mag

def whitehead(Zfunc, C_star, N, L):
    h = 2*L/(N-1)
    ax = np.linspace(-L,L,N)
    X,Y,Z = np.meshgrid(ax,ax,ax,indexing='ij')
    pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    z1,z2 = Zfunc(pts,C_star)
    z1=z1.reshape(N,N,N); z2=z2.reshape(N,N,N)
    dz1=[np.gradient(z1,h,axis=a) for a in range(3)]
    dz2=[np.gradient(z2,h,axis=a) for a in range(3)]
    A=[np.imag(np.conj(z1)*dz1[a]+np.conj(z2)*dz2[a]) for a in range(3)]
    curlA=[np.gradient(A[2],h,axis=1)-np.gradient(A[1],h,axis=2),
           np.gradient(A[0],h,axis=2)-np.gradient(A[2],h,axis=0),
           np.gradient(A[1],h,axis=0)-np.gradient(A[0],h,axis=1)]
    integrand = sum(A[i]*curlA[i] for i in range(3))
    return h**3*np.sum(integrand)/(4*np.pi**2)

print(f"{'C*':>8}  {'tube_r':>8}  {'Q_H (3-strand)':>14}")
t0=time.time()
for C_star in [1.5, 1.1442, 0.8, 0.6, 0.4, 0.3, 0.2]:
    qh = whitehead(Z_s3lift_3strand, C_star, 50, 10.0)
    print(f"{C_star:>8.4f}  {1/C_star:>8.3f}  {qh:>14.5f}  (t={time.time()-t0:.0f}s)")
