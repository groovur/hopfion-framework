"""
Resolution scan: does the energy gap between ansatz A (discontinuous)
and ansatz B (smoothed) GROW as grid resolution increases (the
signature of a genuine discontinuity, whose gradient energy formally
diverges in the continuum limit), or does it stay BOUNDED (the
signature of a real but finite physical feature, e.g. a steep but
smooth transition)?
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

def build_n_A(pts):
    dist, idx = tree.query(pts)
    rho = np.clip(dist, 1e-6, None)
    t_near = t_arr[idx]
    f = f0(rho*C_star)
    theta = 3*t_near
    n = np.stack([np.sin(f)*np.cos(theta), np.sin(f)*np.sin(theta), np.cos(f)], axis=-1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-10)

def build_n_B(pts):
    d_over = np.linalg.norm(pts - p_over, axis=1)
    d_under = np.linalg.norm(pts - p_under, axis=1)
    rho = np.minimum(d_over, d_under)
    f = f0(np.clip(rho,1e-6,None)*C_star)
    theta_over, theta_under = 3*t1, 3*t2
    w_over = 1.0/np.clip(d_over,1e-6,None)**2
    w_under = 1.0/np.clip(d_under,1e-6,None)**2
    ex = w_over*np.cos(theta_over) + w_under*np.cos(theta_under)
    ey = w_over*np.sin(theta_over) + w_under*np.sin(theta_under)
    theta_blend = np.arctan2(ey, ex)
    n = np.stack([np.sin(f)*np.cos(theta_blend), np.sin(f)*np.sin(theta_blend), np.cos(f)], axis=-1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-10)

def gradient_energy_density(n, h):
    g2 = np.zeros(n.shape[:3])
    for axis in range(3):
        dn = (np.roll(n, -1, axis=axis) - np.roll(n, 1, axis=axis)) / (2*h)
        g2 += np.sum(dn**2, axis=-1)
    return g2

print("="*70)
print("RESOLUTION SCAN: DOES THE A/B ENERGY GAP GROW WITH GRID DENSITY?")
print("="*70)
print(f"{'Ngrid':>6}  {'h':>10}  {'E_A':>12}  {'E_B':>12}  {'ratio':>8}  {'maxg2_A':>10}  {'maxg2_B':>10}  {'maxratio':>9}")

L = 1.5*r0
for Ngrid in [41, 61, 81, 121, 161, 241]:
    h = 2*L/(Ngrid-1)
    ax = np.linspace(-L, L, Ngrid)
    X, Y, Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    n_A = build_n_A(pts).reshape(Ngrid,Ngrid,Ngrid,3)
    n_B = build_n_B(pts).reshape(Ngrid,Ngrid,Ngrid,3)
    g2_A = gradient_energy_density(n_A, h)
    g2_B = gradient_energy_density(n_B, h)

    margin = max(5, Ngrid//16)
    core = slice(margin, Ngrid-margin)
    g2_A_core = g2_A[core,core,core]
    g2_B_core = g2_B[core,core,core]

    EA = g2_A_core.sum()*h**3
    EB = g2_B_core.sum()*h**3
    maxA = g2_A_core.max()
    maxB = g2_B_core.max()
    print(f"{Ngrid:>6}  {h:>10.5f}  {EA:>12.4f}  {EB:>12.4f}  {EA/EB:>8.4f}  {maxA:>10.2f}  {maxB:>10.2f}  {maxA/maxB:>9.4f}")
