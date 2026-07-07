"""
Sweep the GLOBAL S3-lift ansatz's Whitehead-integral Hopf charge
across C* values straddling C*_crit = 1/r0, to see whether Q_H
remains pinned at 0 throughout, or develops a nonzero value as the
two-strand blending region grows past its originally-intended local
scope.
"""
import sys, time
import numpy as np
from scipy.spatial import KDTree

R0, r0 = 3.0, 0.874
C_star_crit = 1.0/r0
NT = 6000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
tree = KDTree(np.stack([Gx,Gy,Gz], axis=1))

def f0(rho_hat, C_star):
    return 2*np.arctan(np.maximum(rho_hat, 1e-9)**(-C_star))

def nearest_two(pts, k=80, min_sep=np.pi/2):
    dists, idxs = tree.query(pts, k=k)
    t1 = t_arr[idxs[:,0]]; d1 = dists[:,0]
    dt = np.abs(((t_arr[idxs] - t1[:,None] + np.pi) % (2*np.pi)) - np.pi)
    mask = dt > min_sep
    fi = np.argmax(mask, axis=1)
    has2 = mask.any(axis=1)
    d2 = np.where(has2, dists[np.arange(len(pts)), fi], 1e6)
    t2 = t_arr[idxs[np.arange(len(pts)), fi]]
    return t1, d1, t2, d2

def Z_s3lift(pts, C_star):
    t1, d1, t2, d2 = nearest_two(pts)
    f1 = f0(np.clip(d1, 1e-6, None)*C_star, C_star)
    f2 = f0(np.clip(d2, 1e-6, None)*C_star, C_star)
    w1 = 1.0/np.clip(d1, 1e-6, None)**2
    w2 = 1.0/np.clip(d2, 1e-6, None)**2
    z1u = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2u = w1*np.sin(f1/2)*np.exp(3j*t1) + w2*np.sin(f2/2)*np.exp(3j*t2)
    mag = np.sqrt(np.abs(z1u)**2 + np.abs(z2u)**2)
    return z1u/mag, z2u/mag

def whitehead(Zfunc, C_star, N, L):
    h = 2*L/(N-1)
    ax = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    z1, z2 = Zfunc(pts, C_star)
    z1 = z1.reshape(N,N,N); z2 = z2.reshape(N,N,N)
    dz1 = [np.gradient(z1, h, axis=a) for a in range(3)]
    dz2 = [np.gradient(z2, h, axis=a) for a in range(3)]
    A = [np.imag(np.conj(z1)*dz1[a] + np.conj(z2)*dz2[a]) for a in range(3)]
    curlA = [
        np.gradient(A[2], h, axis=1) - np.gradient(A[1], h, axis=2),
        np.gradient(A[0], h, axis=2) - np.gradient(A[2], h, axis=0),
        np.gradient(A[1], h, axis=0) - np.gradient(A[0], h, axis=1),
    ]
    integrand = sum(A[i]*curlA[i] for i in range(3))
    QH = h**3 * np.sum(integrand) / (4*np.pi**2)
    return QH

N = 50
L = 7.0
print(f"C*_crit = {C_star_crit:.4f}")
print(f"Grid: N={N}, L={L}, h={2*L/(N-1):.4f}")
print(f"\n{'C*':>8}  {'tube_r':>8}  {'Q_H (S3-lift)':>14}")

C_star_values = [2.5062, 2.0, 1.5, 1.3, C_star_crit, 1.0, 0.8, 0.6, 0.4, 0.2]
t0 = time.time()
for C_star in C_star_values:
    qh = whitehead(Z_s3lift, C_star, N, L)
    print(f"{C_star:>8.4f}  {1/C_star:>8.3f}  {qh:>14.5f}   (t={time.time()-t0:.0f}s)")
