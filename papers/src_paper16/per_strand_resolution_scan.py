import numpy as np
from perstrand_s3lift_test import Z_perstrand_s3lift, curve_at_t

t1c = np.pi/6
p_over = curve_at_t(np.array([t1c]))[0]
p_under = curve_at_t(np.array([t1c+np.pi]))[0]
midpoint = (p_over+p_under)/2

def build_n(pts):
    z1,z2 = Z_perstrand_s3lift(pts)
    nx=2*np.real(np.conj(z1)*z2); ny=2*np.imag(np.conj(z1)*z2); nz=np.abs(z1)**2-np.abs(z2)**2
    n = np.stack([nx,ny,nz],axis=-1)
    return n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)

def gradient_energy_density(n, h):
    g2 = np.zeros(n.shape[:3])
    for axis in range(3):
        dn = (np.roll(n,-1,axis=axis)-np.roll(n,1,axis=axis))/(2*h)
        g2 += np.sum(dn**2, axis=-1)
    return g2

r0 = 0.874
L = 1.5*r0
print("Resolution scan: per-strand Phi=chi+3t, S3-lift blended")
print(f"{'Ngrid':>6}  {'h':>9}  {'maxg2':>10}  {'E_total':>10}")
for Ngrid in [41, 61, 81, 121]:
    h = 2*L/(Ngrid-1)
    ax = np.linspace(-L,L,Ngrid)
    X,Y,Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    n = build_n(pts).reshape(Ngrid,Ngrid,Ngrid,3)
    g2 = gradient_energy_density(n,h)
    margin = max(5,Ngrid//16)
    core = slice(margin,Ngrid-margin)
    g2c = g2[core,core,core]
    print(f"{Ngrid:>6}  {h:>9.5f}  {g2c.max():>10.2f}  {(g2c.sum()*h**3):>10.2f}")
