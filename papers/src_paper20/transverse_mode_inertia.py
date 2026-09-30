#!/usr/bin/env python3.11
"""
transverse_mode_inertia.py -- the breathing INERTIA for the PINNED transverse mode (2026-09-29).

pin_displacement_direction.py pinned the confined-oscillation displacement to TRANSVERSE-to-the-tube
(force ratio F_up/F_down=1.375, direction-independent in the {Nhat,Bhat} plane). The Airy residual
E0~(F^2/mu)^{1/3} needs mu for the SAME mode. Earlier mu_up/mu_down=0.886 was the DILATION value
(breathing_inertia) and 0.885 the radial-xy rigid value (local_bending) -- robust, but not the pure
transverse (Frenet-normal) rigid displacement. This computes mu = INT |dn/dA|^2 dV (tangent-projected on
S^2) for a rigid displacement A*w(t)*Nhat(t) of the crossing(up) vs midpoint(down) network, to lock the
coefficient. Nhat(t) built by finite difference on a fine grid, nearest-index lookup at query params.
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062
N=int(sys.argv[1]) if len(sys.argv)>1 else 64
h=(2*(R0+r0+1/C_star+0.6))/N
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)
NTf=20000; t_frame,_,N1f,N2f,H=build_compensated_frame_arclength(NT=NTf)
def frame_at_t(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NTf; return N1f[idx],N2f[idx]
def base_curve(t): return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)],-1)
# --- Frenet normal Nhat(t) on a fine grid, nearest-index lookup ---
NF=200000; tf=np.linspace(0,2*np.pi,NF,endpoint=False); dtf=tf[1]-tf[0]
Gf=base_curve(tf)
dG=(np.roll(Gf,-1,0)-np.roll(Gf,1,0))/(2*dtf); spf=np.linalg.norm(dG,axis=-1,keepdims=True); That=dG/spf
dThat=(np.roll(That,-1,0)-np.roll(That,1,0))/(2*dtf)/spf
kap=np.linalg.norm(dThat,axis=-1,keepdims=True); Nf=dThat/np.clip(kap,1e-12,None)
def nhat(t): return Nf[(np.searchsorted(tf,t%(2*np.pi)))%NF]        # Frenet normal at param t

NT=4000; t_arr=np.linspace(0,2*np.pi,NT,endpoint=False)
cross_c=[np.pi/6+2*np.pi*k/3 for k in range(3)]
mid_c  =[k*np.pi/3 for k in range(6)]
def window(t, centers, sw=0.20):
    w=np.zeros_like(t)
    for c in centers:
        dt=np.abs(((t-c+np.pi)%(2*np.pi))-np.pi); w+=np.exp(-(dt/sw)**2)
    return w

def director_for_displacement(A, centers):
    disp_t = base_curve(t_arr) + (A*window(t_arr,centers))[:,None]*nhat(t_arr)
    lobe_idx=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in [0,2*np.pi/3,4*np.pi/3]]
    trees=[KDTree(disp_t[li]) for li in lobe_idx]; tls=[t_arr[li] for li in lobe_idx]
    dp,tp=[],[]
    for tr,tl in zip(trees,tls):
        d,i=tr.query(pts,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    t1=np.take_along_axis(ts,o,1)[:,0]; d1=np.take_along_axis(ds,o,1)[:,0]
    t2=np.take_along_axis(ts,o,1)[:,1]; d2=np.take_along_axis(ds,o,1)[:,1]
    def disp_curve_at(tq): return base_curve(tq)+(A*window(tq,centers))[:,None]*nhat(tq)
    def chi(tq):
        n1,n2=frame_at_t(tq); rel=pts-disp_curve_at(tq); return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
    P1=chi(t1)+3*t1; P2=chi(t2)+3*t2; r1=np.clip(d1,1e-6,None); r2=np.clip(d2,1e-6,None)
    f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    f1=f0(r1*C_star); f2=f0(r2*C_star); w1=1/r1**2; w2=1/r2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
    n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),np.abs(z1)**2-np.abs(z2)**2],-1)
    return n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10), t1

t0=time.time(); A=0.05
res={}
for name,centers,ns in [('crossing(up)',cross_c,3),('midpoint(down)',mid_c,6)]:
    n0,t1=director_for_displacement(0.0,centers)
    npf,_=director_for_displacement(A,centers); nmf,_=director_for_displacement(-A,centers)
    dn=(npf-nmf)/(2*A); dn=dn-(dn*n0).sum(-1,keepdims=True)*n0        # tangent-projected mode shape
    mu=(dn**2).sum()*h**3/ns
    res[name]=mu
    print(f"  {name:16s}: mu_transverse/site = {mu:.3f}")
print(f"  (built in {time.time()-t0:.1f}s, N={N}, A={A}, displacement = Frenet normal)")
mu_u=res['crossing(up)']; mu_d=res['midpoint(down)']
print(f"\n  mu_transverse_up/mu_transverse_down = {mu_u/mu_d:.4f}")
print(f"  (compare: dilation 0.886, radial-xy rigid 0.885 -- expect ~0.886 if robust)")
F=1.375   # pinned tension force ratio F_up/F_down (pin_displacement_direction.py)
mult=1/((F**2/(mu_u/mu_d))**(1/3))
print(f"\n  PINNED Airy coefficient: M_dn/M_up mult = 1/(F^2/mu_ratio)^(1/3) = 1/({F}^2/{mu_u/mu_d:.4f})^(1/3) = {mult:.4f}")
print(f"  => M_dn/M_up = 0.189 x {mult:.4f} = {0.189*mult:.4f}   (measured 0.157; target mult 0.831)")
