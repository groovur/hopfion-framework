"""
Sweep C* from above to well below C*_crit = 1/r0, at FIXED grid
resolution and FIXED box, and track:
  (a) total gradient energy of the single-nearest-point ansatz
  (b) the SIZE of the "ambiguity set" (fraction of near-trefoil
      volume where the gap between 1st and 2nd nearest-distinct-
      strand distances is small) that is physically WEIGHTED by the
      field (i.e. weighted by sin^2(f0(rho)), so it only counts
      where the field actually has support)
  (c) whether the ambiguity region remains 3 disjoint blobs or
      merges into a connected whole (via connected-components count)
"""
import numpy as np
from scipy.spatial import KDTree
import networkx as nx

R0, r0 = 3.0, 0.874
C_star_crit = 1.0/r0
print(f"C*_crit = 1/r0 = {C_star_crit:.4f}")

NT = 4000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
tree = KDTree(Gamma_pts)

def f0(rho_hat, C_star):
    return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

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

def build_n(pts, C_star):
    t1, d1, t2, d2 = nearest_two(pts)
    rho = np.clip(d1, 1e-6, None)
    f = f0(rho*C_star, C_star)
    theta = 3*t1
    n = np.stack([np.sin(f)*np.cos(theta), np.sin(f)*np.sin(theta), np.cos(f)], axis=-1)
    n = n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
    return n, d1, d2

def gradient_energy_density(n, h):
    g2 = np.zeros(n.shape[:3])
    for axis in range(3):
        dn = (np.roll(n, -1, axis=axis) - np.roll(n, 1, axis=axis)) / (2*h)
        g2 += np.sum(dn**2, axis=-1)
    return g2

# FIXED box and grid, generous enough to contain the field's support
# even at the smallest C* tested
L_xy = R0 + 4.0
L_z = 4.0
N = 90
h = (2*L_xy)/(N-1)
xs = np.linspace(-L_xy, L_xy, N)
ys = np.linspace(-L_xy, L_xy, N)
zs = np.linspace(-L_z, L_z, N)
X,Y,Z = np.meshgrid(xs,ys,zs, indexing='ij')
pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
print(f"Fixed grid: N={N}, h={h:.4f}, box [-{L_xy:.1f},{L_xy:.1f}]^2 x [-{L_z:.1f},{L_z:.1f}]")

C_star_values = [2.5062, 2.0, 1.5, 1.3, C_star_crit, 1.0, 0.7, 0.5, 0.3, 0.2, 0.1]

print(f"\n{'C*':>8}  {'tube_r':>8}  {'E_total':>12}  {'maxg2':>10}  {'ambig_frac':>11}")
for C_star in C_star_values:
    n, d1, d2 = build_n(pts, C_star)
    n3 = n.reshape(N,N,N,3)
    g2 = gradient_energy_density(n3, h)
    E_total = g2.sum()*h**3
    maxg2 = g2.max()

    # ambiguity measure weighted by field support: sin^2(f0(d1)) * [d2/d1 close to 1]
    rho = np.clip(d1,1e-6,None)
    f = f0(rho*C_star, C_star)
    field_weight = np.sin(f)**2
    ambig = (d2 - d1)/np.clip(d1,1e-6,None)  # small near discontinuity
    is_ambig = ambig < 0.5  # within 50% -- a generous proxy threshold
    ambig_weighted_frac = (field_weight*is_ambig).sum() / np.clip(field_weight.sum(),1e-10,None)

    tube_r = 1.0/C_star
    print(f"{C_star:>8.4f}  {tube_r:>8.3f}  {E_total:>12.3f}  {maxg2:>10.2f}  {ambig_weighted_frac:>11.4f}")
