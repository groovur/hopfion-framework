"""
Resolution scan for the COMPLEX-DOUBLET-SUPERPOSITION (genuine S^3
lift) ansatz, compared against ansatz A (discontinuous, current
solver) and ansatz B (real-phase blend, also discontinuous) from
crossing_resolution_scan.py.

Ansatz C: build each strand's own Hopf doublet Z_over, Z_under
(complex, S^3-valued) using ITS OWN local profile f(d_strand), with
a FIXED internal phase theta per strand (theta_over=3*t1,
theta_under=3*t2). Combine via inverse-square-distance-weighted
COMPLEX superposition (not real-angle averaging), normalise, then
project to n via the Hopf map. This is a genuine two-channel S^3
construction, distinct from the single-nearest-point ansatz (A) and
the real-phase blend (B).
"""
import numpy as np
from scipy.spatial import KDTree

R0, r0, C_star = 3.0, 0.874, 2.5062

def f0(rho_hat):
    return 2*np.arctan(np.maximum(rho_hat, 1e-9)**(-C_star))

NT = 20000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
tree = KDTree(Gamma_pts)

t1 = np.pi/6
p_over = np.array([(R0+r0*np.cos(3*t1))*np.cos(2*t1),
                    (R0+r0*np.cos(3*t1))*np.sin(2*t1),
                    r0*np.sin(3*t1)])
t2 = t1+np.pi
p_under = np.array([(R0+r0*np.cos(3*t2))*np.cos(2*t2),
                     (R0+r0*np.cos(3*t2))*np.sin(2*t2),
                     r0*np.sin(3*t2)])
midpoint = (p_over+p_under)/2
theta_over_fixed, theta_under_fixed = 3*t1, 3*t2

def build_n_A(pts):
    dist, idx = tree.query(pts)
    rho = np.clip(dist, 1e-6, None)
    t_near = t_arr[idx]
    f = f0(rho*C_star)
    theta = 3*t_near
    n = np.stack([np.sin(f)*np.cos(theta), np.sin(f)*np.sin(theta), np.cos(f)], axis=-1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-10)

def build_n_C_S3lift(pts):
    """Genuine S^3-lift two-channel ansatz."""
    d_over = np.linalg.norm(pts - p_over, axis=1)
    d_under = np.linalg.norm(pts - p_under, axis=1)
    f_over = f0(np.clip(d_over,1e-6,None)*C_star)
    f_under = f0(np.clip(d_under,1e-6,None)*C_star)

    z1_over = np.cos(f_over/2).astype(complex)
    z2_over = np.sin(f_over/2)*np.exp(1j*theta_over_fixed)
    z1_under = np.cos(f_under/2).astype(complex)
    z2_under = np.sin(f_under/2)*np.exp(1j*theta_under_fixed)

    w_over = 1.0/np.clip(d_over,1e-6,None)**2
    w_under = 1.0/np.clip(d_under,1e-6,None)**2

    z1 = w_over*z1_over + w_under*z1_under
    z2 = w_over*z2_over + w_under*z2_under
    norm = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2)
    z1, z2 = z1/norm, z2/norm

    nx = 2*np.real(np.conj(z1)*z2)
    ny = 2*np.imag(np.conj(z1)*z2)
    nz = np.abs(z1)**2 - np.abs(z2)**2
    n = np.stack([nx,ny,nz], axis=-1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-10)

def gradient_energy_density(n, h):
    g2 = np.zeros(n.shape[:3])
    for axis in range(3):
        dn = (np.roll(n, -1, axis=axis) - np.roll(n, 1, axis=axis)) / (2*h)
        g2 += np.sum(dn**2, axis=-1)
    return g2

print("="*70)
print("RESOLUTION SCAN: ANSATZ A (discontinuous) vs ANSATZ C (S^3 lift)")
print("="*70)
print(f"{'Ngrid':>6}  {'h':>10}  {'E_A':>12}  {'E_C':>12}  {'ratio A/C':>10}  {'maxg2_A':>10}  {'maxg2_C':>10}  {'maxratio':>9}")

L = 1.5*r0
for Ngrid in [41, 61, 81, 121, 161, 241]:
    h = 2*L/(Ngrid-1)
    ax = np.linspace(-L, L, Ngrid)
    X, Y, Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    n_A = build_n_A(pts).reshape(Ngrid,Ngrid,Ngrid,3)
    n_C = build_n_C_S3lift(pts).reshape(Ngrid,Ngrid,Ngrid,3)
    g2_A = gradient_energy_density(n_A, h)
    g2_C = gradient_energy_density(n_C, h)

    margin = max(5, Ngrid//16)
    core = slice(margin, Ngrid-margin)
    g2_A_core = g2_A[core,core,core]
    g2_C_core = g2_C[core,core,core]

    EA = g2_A_core.sum()*h**3
    EC = g2_C_core.sum()*h**3
    maxA = g2_A_core.max()
    maxC = g2_C_core.max()
    print(f"{Ngrid:>6}  {h:>10.5f}  {EA:>12.4f}  {EC:>12.4f}  {EA/EC:>10.4f}  {maxA:>10.2f}  {maxC:>10.2f}  {maxA/maxC:>9.4f}")
