import numpy as np
from scipy.spatial import KDTree

print("="*70)
print("DENSE, FOCUSED SCAN OF |Z| RIGHT AT ALL THREE CROSSING REGIONS")
print("="*70)

R0, r0, C_star = 3.0, 0.874, 2.5062
def f0(rho_hat): return 2*np.arctan(np.maximum(rho_hat,1e-9)**(-C_star))

NT = 8000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma_pts = np.stack([Gx,Gy,Gz], axis=1)
tree = KDTree(Gamma_pts)

def nearest_two_distinct(pts, min_param_sep=np.pi/2, k=80):
    dists, idxs = tree.query(pts, k=k)
    t1 = t_arr[idxs[:,0]]
    d1 = dists[:,0]
    t_all = t_arr[idxs]
    dt = np.abs(((t_all - t1[:,None] + np.pi) % (2*np.pi)) - np.pi)
    mask_distinct = dt > min_param_sep
    first_distinct = np.argmax(mask_distinct, axis=1)
    has_distinct = mask_distinct.any(axis=1)
    idx2 = idxs[np.arange(len(pts)), first_distinct]
    d2 = dists[np.arange(len(pts)), first_distinct]
    t2 = t_arr[idx2]
    d2 = np.where(has_distinct, d2, 1e6)
    return t1, d1, t2, d2, has_distinct

def Z_unnormalized_magnitude(pts):
    t1, d1, t2, d2, has2 = nearest_two_distinct(pts)
    f1 = f0(np.clip(d1,1e-6,None)*C_star)
    f2 = f0(np.clip(d2,1e-6,None)*C_star)
    theta1, theta2 = 3*t1, 3*t2
    z1a = np.cos(f1/2).astype(complex); z2a = np.sin(f1/2)*np.exp(1j*theta1)
    z1b = np.cos(f2/2).astype(complex); z2b = np.sin(f2/2)*np.exp(1j*theta2)
    w1 = 1.0/np.clip(d1,1e-6,None)**2
    w2 = 1.0/np.clip(d2,1e-6,None)**2
    z1 = w1*z1a + w2*z1b
    z2 = w1*z2a + w2*z2b
    mag = np.sqrt(np.abs(z1)**2+np.abs(z2)**2)
    return mag

t1s = [np.pi/6, 5*np.pi/6, 3*np.pi/2]
def Gamma(t): return np.array([(R0+r0*np.cos(3*t))*np.cos(2*t), (R0+r0*np.cos(3*t))*np.sin(2*t), r0*np.sin(3*t)])

overall_min = np.inf
for i, t1c in enumerate(t1s):
    p_over = Gamma(t1c)
    p_under = Gamma(t1c+np.pi)
    midpoint = (p_over+p_under)/2
    # dense local box around this crossing's midpoint
    L = 1.2*r0
    N = 60
    ax = np.linspace(-L, L, N)
    X,Y,Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    pts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    mag = Z_unnormalized_magnitude(pts)
    print(f"Crossing {i+1}: min|Z| in local box = {mag.min():.6f}  (at {pts[mag.argmin()]})")
    overall_min = min(overall_min, mag.min())

print(f"\nOverall minimum |Z| across all three crossing regions: {overall_min:.6f}")
print(f"Zero found: {overall_min < 1e-6}")
