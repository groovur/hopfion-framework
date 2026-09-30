#!/usr/bin/env python3.11
"""
isospin_feedback_residual.py -- does DENSITY FEEDBACK break the crossing/midpoint energy equality,
and can it be the ~1.2 isospin residual (M_down/M_up predicted 0.189 vs measured 0.157)?

faddeev_energy_crossing_vs_midpoint.py found the BARE (beta=0) Faddeev energy density EQUAL on the
crossing (up, z=+r0) and midpoint (down, z=0) networks (quartic 0.94, quadratic 1.00). But the crossing
is SHARP (two strands meet, high g2) and the midpoint is GENTLE (single strand, lower g2). The density
feedback g2/(1+beta*g2) is NONLINEAR: it SATURATES (caps at 1/beta) where g2 is high (crossing) and
stays RESPONSIVE where g2 is low (midpoint). So even with equal SUM of g2, the feedback-softened isotropic
energy can DIFFER between the networks. This tests whether that differential is ~1.2 and in the direction
that makes DOWN lighter (M_down/M_up 0.189 -> 0.157).
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062; MU=3.0-PHI
N=int(sys.argv[1]) if len(sys.argv)>1 else 56
h=(2*(R0+r0+1/C_star+0.6))/N
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)
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
    return (np.take_along_axis(ts,o,1)[:,0],np.take_along_axis(ds,o,1)[:,0])
t1,d1=nearest_two(pts)
def chi(tq):
    n1,n2=frame_at_t(tq); rel=pts-curve_at_t(tq); return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
_,d1b=nearest_two(pts)  # ensure order; reuse t1,d1
Phi1=chi(t1)+3*t1; rho1=np.clip(d1,1e-6,None)
# second strand for the blend
def nearest_two_full(q):
    dp,tp=[],[]
    for tr,tl in zip(lobe_trees,lobe_t):
        d,i=tr.query(q,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    return (np.take_along_axis(ts,o,1)[:,0],np.take_along_axis(ds,o,1)[:,0],
            np.take_along_axis(ts,o,1)[:,1],np.take_along_axis(ds,o,1)[:,1])
t1,d1,t2,d2=nearest_two_full(pts)
Phi1=chi(t1)+3*t1; Phi2=chi(t2)+3*t2; rho1=np.clip(d1,1e-6,None); rho2=np.clip(d2,1e-6,None)
f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star))
f1=f0(rho1*C_star); f2=f0(rho2*C_star); w1=1/rho1**2; w2=1/rho2**2
z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
z2=w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),np.abs(z1)**2-np.abs(z2)**2],-1).reshape(N,N,N,3)
n/=np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
nx,ny,nz=n[...,0],n[...,1],n[...,2]
dd=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]
g2=sum(dd[c][a]**2 for c in range(3) for a in range(3)).reshape(-1)
s4=((1-nz**2).clip(0,1)**2).reshape(-1)

sin3=np.sin(3*t1)
up=(sin3>0.5); dn=(np.abs(sin3)<0.5)
fu=(np.sin(3*t_arr)>0.5).mean(); fd=(np.abs(np.sin(3*t_arr))<0.5).mean()
def ratio(dens):  # per-equal-arc up/down
    return (dens[up].sum()/fu)/(dens[dn].sum()/fd)
print(f"N={N} h={h:.3f}  MU={MU:.4f}")
print(f"  BARE quadratic g2 up/down = {ratio(g2):.4f}   (script: ~1.00)")
print(f"  BARE J2a (s4*g2)  up/down = {ratio(s4*g2):.4f}")
print(f"\n  === feedback-softened ISOTROPIC density g2/(1+beta g2), up/down per arc ===")
print(f"  {'beta':>6} {'gfb up/down':>12} {'shift vs bare':>14} {'-> M_dn/M_up x':>16}")
for beta in [0.0,0.1,0.2,0.35,0.452,0.7,1.0]:
    gfb=g2/(1+beta*g2)
    r=ratio(gfb); shift=r/ratio(g2)
    # if mass ~ isotropic energy, up heavier by 'r' => M_dn/M_up multiplied by 1/r relative to bare
    print(f"  {beta:>6.3f} {r:>12.4f} {shift:>14.4f} {1/shift:>16.4f}")
print(f"\n  target: measured/predicted = 0.157/0.189 = {0.157/0.189:.4f} (need down lighter by this)")
print("  READ: if 'M_dn/M_up x' reaches ~0.83 for a physical beta and >1 gfb ratio (up heavier),")
print("  the feedback differential is a candidate for the ~1.2 residual. If gfb up/down stays ~1,")
print("  feedback does NOT break the crossing/midpoint equality -> not the residual.")
