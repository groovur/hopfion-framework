#!/usr/bin/env python3.11
"""
breathing_inertia_crossing_vs_midpoint.py -- the STRETCH/breathing mode of the confined trefoil (2026-09-29).

Context: E=K*J4 is SCALE-INVARIANT (this session) -> the trefoil's SIZE lambda is a flat zero mode
(the "breathing"/stretch coordinate, quark_sector_confinement_cherenkov.md sec.2-3: string tension
snaps it back; confined quark OSCILLATES, E_q/E_c=phi). The breathing has ZERO Faddeev stiffness
(scale-flat) -> its restoring force is the external STRING; its frequency is omega=sqrt(k_string/mu).
So the Faddeev-computable, chi-independent-in-the-ratio piece is the breathing INERTIA
    mu = INT |d n/d eps|^2 dV,   with the dilation deformation n_eps(x)=n(x/(1+eps)) -> dn/deps = -(x.grad)n.
This script computes mu, its CROSSING(up)/MIDPOINT(down) split (the ~1.2 residual test: does the breathing
inertia differ between the up and down networks?), and the dilation-derivative structure. chi (cutoff-set,
sec.3 'magnitudes open') cancels in the up/down ratio if uniform.
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062
N=int(sys.argv[1]) if len(sys.argv)>1 else 64
h=(2*(R0+r0+1/C_star+0.6))/N
cv=h*(np.arange(N)-N//2+0.5)
X=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1)          # (N,N,N,3) lab coords about the trefoil centre
pts=X.reshape(-1,3)
NTf=20000; t_frame,_,N1f,N2f,H=build_compensated_frame_arclength(NT=NTf)
def curve_at_t(t): return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)],-1)
def frame_at_t(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NTf; return N1f[idx],N2f[idx]
NT=4000; t_arr=np.linspace(0,2*np.pi,NT,endpoint=False); Gam=curve_at_t(t_arr)
lobe_idx=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in [0,2*np.pi/3,4*np.pi/3]]
lobe_trees=[KDTree(Gam[li]) for li in lobe_idx]; lobe_t=[t_arr[li] for li in lobe_idx]
def nearest_two(q):
    dp,tp=[],[]
    for tr,tl in zip(lobe_trees,lobe_t):
        d,i=tr.query(q,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    return (np.take_along_axis(ts,o,1)[:,0],np.take_along_axis(ds,o,1)[:,0],
            np.take_along_axis(ts,o,1)[:,1],np.take_along_axis(ds,o,1)[:,1])
def build_director(P):
    t1,d1,t2,d2=nearest_two(P)
    def chi(tq):
        n1,n2=frame_at_t(tq); rel=P-curve_at_t(tq); return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
    P1=chi(t1)+3*t1; P2=chi(t2)+3*t2; r1=np.clip(d1,1e-6,None); r2=np.clip(d2,1e-6,None)
    f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star))
    f1=f0(r1*C_star); f2=f0(r2*C_star); w1=1/r1**2; w2=1/r2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
    n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),np.abs(z1)**2-np.abs(z2)**2],-1)
    return n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10), t1

t0=time.time()
n0,t1=build_director(pts)                                   # eps=0
# dilation derivative dn/deps = -(x.grad)n : build n at scaled coords x/(1+eps) for +-eps, central diff
eps=0.01
np_p,_=build_director(pts/(1+eps)); np_m,_=build_director(pts/(1-eps))
dn=(np_p-np_m)/(2*eps)                                      # d n / d eps  (tangent-projected below)
dn=dn-(dn*n0).sum(-1,keepdims=True)*n0                     # project to S^2 tangent (physical mode)
inertia=(dn**2).sum(-1)                                     # |dn/deps|^2 density  (chi factored out)
print(f"  built in {time.time()-t0:.1f}s   N={N} h={h:.3f}")
mu=inertia.sum()*h**3
print(f"  breathing inertia mu = INT |dn/deps|^2 dV = {mu:.3f}  (x chi, cutoff-set)")

# crossing(up)/midpoint(down) partition by nearest-strand site signature sin(3 t1)
sin3=np.sin(3*t1)
up=(sin3>0.5); dn_b=(np.abs(sin3)<0.5)
fu=(np.sin(3*t_arr)>0.5).mean(); fd=(np.abs(np.sin(3*t_arr))<0.5).mean()
inr=inertia.reshape(-1)
mu_up=inr[up].sum()*h**3/fu; mu_dn=inr[dn_b].sum()*h**3/fd
print(f"\n  === breathing inertia per equal arc, crossing(up) vs midpoint(down) ===")
print(f"  mu_up(cross)/arc = {mu_up:.3f}   mu_down(mid)/arc = {mu_dn:.3f}")
print(f"  mu_up/mu_down = {mu_up/mu_dn:.4f}")
# omega ~ sqrt(k_string/mu); zero-point E0 ~ hbar omega ~ mu^-1/2 (equal string stiffness).
# => M_up/M_down from breathing ~ (mu_down/mu_up)^{1/2}; residual multiplier on M_dn/M_up ~ (mu_up/mu_down)^{1/2}
r=mu_up/mu_dn
print(f"\n  IF omega~mu^-1/2 with equal string stiffness:")
print(f"    M_dn/M_up breathing-multiplier = (mu_up/mu_down)^(1/2) = {r**0.5:.4f}")
print(f"    (Airy/linear omega~mu^-1/3 -> {r**(1/3):.4f})")
print(f"  target = 0.157/0.189 = {0.157/0.189:.4f}  (need this from the breathing differential)")
print(f"\n  READ: mu_up/mu_down != 1 => the breathing/stretch dynamics DOES differ between up and down")
print(f"  networks (unlike the static feedback). Compare the multiplier to 0.83.")
