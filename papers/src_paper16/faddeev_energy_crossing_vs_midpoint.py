#!/usr/bin/env python3.11
"""
Option 1 (focused field-energy calc): the Faddeev energy DENSITY localized on the CROSSING network
(up-type, z=+r0) vs the MIDPOINT/distal network (down-type, z=0) of the golden trefoil.
Session 2026-09-21. Tests sec.27/28 at the ENERGY level: does the crossing (up) network carry a
distinct Faddeev energy from the midpoint (down) network? (The twist-SCALING is deferred: sec.29's
mass~|Tw|^4.2 is NOT a geometric-twist energy, which scales ~Tw^2 -- the |Tw|->config mapping is open.)

Reuses the validated Construction-C trefoil ansatz + the E_geom functional of
gradient_flow_constrained.py. Single field evaluation (no gradient flow): light.
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength

PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062; MU=1.0
N=int(sys.argv[1]) if len(sys.argv)>1 else 40
h=(2*(R0+r0+1/C_star+0.6))/N
print(f"N={N} h={h:.4f} box_hw={N*h/2:.2f} (need {R0+r0+1/C_star+0.5:.2f})")

# ---- grid ----
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)

# ---- trefoil curve + compensated Bishop frame ----
t0=time.time()
NTf=20000
t_frame,_,N1f,N2f,H=build_compensated_frame_arclength(NT=NTf)
def curve_at_t(t): return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                                    (R0+r0*np.cos(3*t))*np.sin(2*t), r0*np.sin(3*t)],-1)
def frame_at_t(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NTf
    return N1f[idx],N2f[idx]
NT=4000; t_arr=np.linspace(0,2*np.pi,NT,endpoint=False)
Gam=curve_at_t(t_arr)
arc_starts=[0,2*np.pi/3,4*np.pi/3]
lobe_idx=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees=[KDTree(Gam[li]) for li in lobe_idx]; lobe_t=[t_arr[li] for li in lobe_idx]
def nearest_two(q):
    dp,tp=[],[]
    for tr,tl in zip(lobe_trees,lobe_t):
        d,i=tr.query(q,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    return (np.take_along_axis(ts,o,1)[:,0],np.take_along_axis(ds,o,1)[:,0],
            np.take_along_axis(ts,o,1)[:,1],np.take_along_axis(ds,o,1)[:,1])
t1,d1,t2,d2=nearest_two(pts)

# ---- Construction-C director (Phi = chi + 3 t) ----
def chi(tq):
    n1,n2=frame_at_t(tq); rel=pts-curve_at_t(tq)
    return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
Phi1=chi(t1)+3*t1; Phi2=chi(t2)+3*t2
rho1=np.clip(d1,1e-6,None); rho2=np.clip(d2,1e-6,None)
f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star))
f1=f0(rho1*C_star); f2=f0(rho2*C_star); w1=1/rho1**2; w2=1/rho2**2
z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
z2=w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),
            np.abs(z1)**2-np.abs(z2)**2],-1).reshape(N,N,N,3)
n/=np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
print(f"  ansatz built in {time.time()-t0:.1f}s")

# ---- Faddeev energy densities (central differences, matching E_geom) ----
def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
nx,ny,nz=n[...,0],n[...,1],n[...,2]
d=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]  # d[comp][axis]
g2=sum(d[c][a]**2 for c in range(3) for a in range(3))
(nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=d
Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
rhoJ4=(Fxy**2+Fxz**2+Fyz**2)
s4=(1-nz**2).clip(0,1)**2
rho_g2=g2; rho_J2a=s4*g2
J4=rhoJ4.sum()*h**3; J2iso=g2.sum()*h**3; J2a=(s4*g2).sum()*h**3
print(f"  J4={J4:.3f}  J2iso={J2iso:.3f}  J2a={J2a:.3f}  E=K*J4={(J2a+MU*J2iso)*J4:.2f}")

# ---- localize by sin(3t*): UP=crossing z=+r0 (sin3t>0.5); DOWN=midpoint/distal z=0 (|sin3t|<0.5);
#      the z=-r0 third (sin3t<-0.5) is the under-strand side (trefoil is NOT z->-z symmetric). ----
sin3=np.sin(3*t1)                          # nearest-strand site signature (flat, matches pts order)
up_bin  =(sin3> 0.5)                       # crossings, z=+r0
dn_bin  =(np.abs(sin3)<0.5)                # midpoints/distal, z=0
un_bin  =(sin3<-0.5)                       # under-strand side, z=-r0
rj=rhoJ4.reshape(-1); gg=rho_g2.reshape(-1)
frac_up=(np.sin(3*t_arr)>0.5).mean(); frac_dn=(np.abs(np.sin(3*t_arr))<0.5).mean()
E_up=rj[up_bin].sum()*h**3; E_dn=rj[dn_bin].sum()*h**3; E_un=rj[un_bin].sum()*h**3
g_up=gg[up_bin].sum()*h**3; g_dn=gg[dn_bin].sum()*h**3
print("\n  === Faddeev QUARTIC (J4) energy by z-network (equal ~1/3 arcs) ===")
print(f"  arc fractions: up(z=+r0)={frac_up:.3f} down(z=0)={frac_dn:.3f}")
print(f"  E_up(cross)={E_up:.2f}  E_down(mid)={E_dn:.2f}  E_under(z=-r0)={E_un:.2f}")
print(f"  UP/DOWN quartic-energy ratio (per equal arc) = {(E_up/frac_up)/(E_dn/frac_dn):.3f}")
print("  === quadratic (g2) ===")
print(f"  g_up={g_up:.2f}  g_down={g_dn:.2f}  UP/DOWN = {(g_up/frac_up)/(g_dn/frac_dn):.3f}")
print(f"\n  INTERPRETATION: ratio>1 => crossing (up) network carries MORE Faddeev energy per unit arc")
print(f"  than midpoint (down) -> up-type energetically distinct (grounds sec.27/28 at the ENERGY")
print(f"  level, beyond tangent orientation). ratio~1 => no energy asymmetry.")
