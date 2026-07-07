import numpy as np
from scipy.spatial import KDTree

R0, r0, C_star = 3.0, 0.874, 2.5062
NT = 8000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
Gamma = np.stack([Gx,Gy,Gz], axis=1).astype(np.float32)
tree = KDTree(Gamma)

def build_n_A(qpts):
    dists, idx = tree.query(qpts)
    rho = np.clip(dists, 1e-6, None)
    t_near = t_arr[idx]
    f_ = 2*np.arctan(rho**(-C_star))
    tht = 3*t_near
    nx = np.sin(f_)*np.cos(tht); ny = np.sin(f_)*np.sin(tht); nz = np.cos(f_)
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
def Gamma_at(t): return np.array([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)])
p_over = Gamma_at(t1c); p_under = Gamma_at(t1c+np.pi)
midpoint = (p_over+p_under)/2

L = 1.5*r0
print("Refinement scan: J4 near a single crossing, Construction A (original, discontinuous)")
print(f"{'Ngrid':>6}  {'h':>9}  {'J4_local':>12}  {'maxJ4dens':>12}")
for Ngrid in [41, 61, 81, 121, 161]:
    h_ = 2*L/(Ngrid-1)
    ax = np.linspace(-L,L,Ngrid)
    X,Y,Z = np.meshgrid(midpoint[0]+ax, midpoint[1]+ax, midpoint[2]+ax, indexing='ij')
    qpts = np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1)
    n = build_n_A(qpts).reshape(Ngrid,Ngrid,Ngrid,3)
    rho_J4 = rho_J4_density(n, h_)
    margin = max(5,Ngrid//16)
    core = slice(margin,Ngrid-margin)
    rj_core = rho_J4[core,core,core]
    J4_local = rj_core.sum()*h_**3
    print(f"{Ngrid:>6}  {h_:>9.5f}  {J4_local:>12.4f}  {rj_core.max():>12.2f}")
