import numpy as np
import sys
sys.path.insert(0,'.')
from bishop_frame_compensated import build_compensated_frame

print("Building a much finer compensated frame (NT=150000) to test the")
print("N=81 spatial-grid energy bump for frame-resolution sensitivity...")

# Monkey-patch a finer frame into the perstrand_s3lift_test module
import perstrand_s3lift_test as mod
t_frame_fine, T_frame_fine, N1_frame_fine, N2_frame_fine, H_fine = build_compensated_frame(NT=150000)
mod.t_frame = t_frame_fine
mod.N1_frame = N1_frame_fine
mod.N2_frame = N2_frame_fine

def frame_at_t_fine(t_query):
    idx = np.searchsorted(t_frame_fine, t_query % (2*np.pi)) % 150000
    return N1_frame_fine[idx], N2_frame_fine[idx]
mod.frame_at_t = frame_at_t_fine

from perstrand_s3lift_test import Z_perstrand_s3lift, curve_at_t
r0 = 0.874
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

L = 1.5*r0
Ngrid = 81
h = 2*L/(Ngrid-1)
ax = np.linspace(-L,L,Ngrid)
X,Y,Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
n = build_n(pts).reshape(Ngrid,Ngrid,Ngrid,3)
g2 = gradient_energy_density(n,h)
margin = max(5,Ngrid//16)
core = slice(margin,Ngrid-margin)
g2c = g2[core,core,core]
print(f"With fine frame (NT=150000): N=81, maxg2={g2c.max():.2f}, E_total={(g2c.sum()*h**3):.2f}")
print("(compare to coarse-frame result: maxg2=72.99, E_total=159.88)")
