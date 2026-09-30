#!/usr/bin/env python3.11
"""
local_bending_stiffness.py -- the STIFFNESS half of the breathing residual (2026-09-29).

omega=sqrt(k/mu); the breathing INERTIA differential is done (mu_up/mu_down=0.886,
breathing_inertia_crossing_vs_midpoint.py). This computes the Faddeev STIFFNESS k under a transverse
(radial = confinement-direction) displacement of the CROSSING(up) network vs the MIDPOINT(down) network:
displace the centerline Gamma(t) -> Gamma(t)+A*w(t)*rhat(t) with w localized on each network (Z_3-coherent),
rebuild the director, and read k = d^2 E_Faddeev/dA^2. If k_up/k_down ~1 the field bending stiffness is
uniform (differential lives in the string/junction coupling, cutoff-set); if k_up>k_down the field itself
carries the other half. Combined M_dn/M_up multiplier = sqrt(k_up mu_down/(k_down mu_up)); target 0.83.
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062; MU=3.0-PHI
N=int(sys.argv[1]) if len(sys.argv)>1 else 64
h=(2*(R0+r0+1/C_star+0.6))/N
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)
NTf=20000; t_frame,_,N1f,N2f,H=build_compensated_frame_arclength(NT=NTf)
def frame_at_t(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NTf; return N1f[idx],N2f[idx]
NT=4000; t_arr=np.linspace(0,2*np.pi,NT,endpoint=False)
def base_curve(t): return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)],-1)
def rhat(t):  return np.stack([np.cos(2*t),np.sin(2*t),np.zeros_like(t)],-1)   # radial (xy) = confinement dir

# Z_3-coherent windows on each network (sum of Gaussians in t, period 2pi)
def window(t, centers, sw=0.20):
    w=np.zeros_like(t)
    for c in centers:
        dt=np.abs(((t-c+np.pi)%(2*np.pi))-np.pi); w+=np.exp(-(dt/sw)**2)
    return w
cross_c=[np.pi/6+2*np.pi*k/3 for k in range(3)]           # 3t=pi/2 (z=+r0)
mid_c  =[k*np.pi/3 for k in range(6)]                     # 3t=0,pi,... (z=0)

def energy_for_displacement(A, centers):
    """Displace centerline radially by A*window on the given network; rebuild director; return E_Faddeev."""
    disp_t = base_curve(t_arr) + (A*window(t_arr,centers))[:,None]*rhat(t_arr)
    lobe_idx=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in [0,2*np.pi/3,4*np.pi/3]]
    trees=[KDTree(disp_t[li]) for li in lobe_idx]; tls=[t_arr[li] for li in lobe_idx]
    dp,tp=[],[]
    for tr,tl in zip(trees,tls):
        d,i=tr.query(pts,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    t1=np.take_along_axis(ts,o,1)[:,0]; d1=np.take_along_axis(ds,o,1)[:,0]
    t2=np.take_along_axis(ts,o,1)[:,1]; d2=np.take_along_axis(ds,o,1)[:,1]
    def disp_curve_at(tq):   # displaced curve at param tq (interp the radial bump)
        return base_curve(tq)+(A*window(tq,centers))[:,None]*rhat(tq)
    def chi(tq):
        n1,n2=frame_at_t(tq); rel=pts-disp_curve_at(tq); return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
    P1=chi(t1)+3*t1; P2=chi(t2)+3*t2; r1=np.clip(d1,1e-6,None); r2=np.clip(d2,1e-6,None)
    f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    f1=f0(r1*C_star); f2=f0(r2*C_star); w1=1/r1**2; w2=1/r2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
    n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),np.abs(z1)**2-np.abs(z2)**2],-1).reshape(N,N,N,3)
    n/=np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
    nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
    dd=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]
    g2=sum(dd[c][a]**2 for c in range(3) for a in range(3))
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=dd
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*h**3; J2a=(( (1-nz**2).clip(0,1)**2)*g2).sum()*h**3; J2iso=g2.sum()*h**3
    return (J2a+MU*J2iso)*J4, n.reshape(-1,3)

t0=time.time(); A=0.05
res={}
for name,centers,nsite in [('crossing(up)',cross_c,3),('midpoint(down)',mid_c,6)]:
    E0,n0=energy_for_displacement(0.0,centers); Ep,npf=energy_for_displacement(A,centers); Em,nmf=energy_for_displacement(-A,centers)
    k=(Ep+Em-2*E0)/A**2
    # SAME-mode inertia: dn/dA of the radial-bending deformation (tangent-projected)
    dn=(npf-nmf)/(2*A); dn=dn-(dn*n0).sum(-1,keepdims=True)*n0
    mu=(dn**2).sum()*h**3
    res[name]=(k/nsite, mu/nsite)
    print(f"  {name:16s}: k/site={k/nsite:.1f}  mu_bend/site={mu/nsite:.2f}  omega~sqrt(k/mu)={np.sqrt((k/nsite)/(mu/nsite)):.3f}")
print(f"  (built in {time.time()-t0:.1f}s, N={N}, A={A})")
ku,mu_u=res['crossing(up)']; kd,mu_d=res['midpoint(down)']
print(f"\n  CONSISTENT radial-bending mode (same k and mu):")
print(f"  k_up/k_down = {ku/kd:.4f}   mu_bend_up/mu_bend_down = {mu_u/mu_d:.4f}")
wu_wd=np.sqrt((ku/kd)/(mu_u/mu_d))
print(f"  omega_up/omega_down = {wu_wd:.4f}  ->  M_dn/M_up multiplier = {1/wu_wd:.4f}   (target 0.83)")
print(f"  [<1 = right direction (down lighter); >1 = wrong direction]")
