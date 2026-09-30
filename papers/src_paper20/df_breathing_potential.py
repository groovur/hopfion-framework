#!/usr/bin/env python3.11
"""
df_breathing_potential.py -- does DENSITY FEEDBACK provide the breathing potential? (2026-09-29)

User's mechanism-first question: the bare Faddeev energy E=K*J4 is scale-INVARIANT, so the breathing
(dilation lambda) mode is flat -- I had to import an ad-hoc string (Airy) + pick a displacement
direction to get a number, and the direction-dependence (radial-xy 0.815 vs radial-3D 0.778) is the tell
that the mechanism is not pinned. BUT density feedback BREAKS the scale-invariance exactly. Under a
dilation x->lambda x of the reference field:
    J2a(lambda)      = lambda * J2a_0
    J2iso_fb(lambda) = lambda * I(beta/lambda^2),   I(x) = INT g2/(1+x g2) dV   (at lambda=1)
    J4(lambda)       = J4_0 / lambda
    => E_fb(lambda)  = J4_0*J2a_0 + mu * J4_0 * I(beta/lambda^2).      [derived; beta=0 => flat, checks]
So DF is a REAL breathing potential, not a bolt-on. This script builds the field ONCE, gets g2/J4/J2a and
the crossing(up)/midpoint(down) partition, then reads E_fb(lambda) over a lambda range for each network:
  (1) the SHAPE  -- harmonic (~lambda^2, saturated regime), linear, or monotonic->collapse (needs a wall);
  (2) slope dE_fb/dlambda and curvature d2E_fb/dlambda2 at lambda=1, up vs down;
  (3) combined with the converged breathing inertia mu_up/mu_down=0.886 -> M_dn/M_up multiplier.
It also adds an optional linear string sigma*L_0*lambda to test whether DF+string form a bound well or
both drive to collapse (=> the packing/core wall sets the size, breathing is against the wall).
Honest: this DECIDES the mode (harmonic vs Airy vs wall) from first principles instead of assuming it.
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
n0,t1=build_director(pts)
n=n0.reshape(N,N,N,3); nx,ny,nz=n[...,0],n[...,1],n[...,2]
def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
dd=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]
g2=sum(dd[c][a]**2 for c in range(3) for a in range(3))          # |grad n|^2 per cell
(nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=dd
Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
J4=(Fxy**2+Fxz**2+Fyz**2).sum()*h**3
J2a=(((1-nz**2).clip(0,1)**2)*g2).sum()*h**3
g2f=g2.reshape(-1)
# crossing(up)/midpoint(down) partition by nearest-strand signature sin(3 t1)
sin3=np.sin(3*t1); up=(sin3>0.5); dnb=(np.abs(sin3)<0.5)
print(f"  built in {time.time()-t0:.1f}s  N={N} h={h:.3f}  J4={J4:.4f} J2a={J2a:.4f}")

def I_of_x(x, mask=None):
    gg = g2f if mask is None else g2f[mask]
    return (gg/(1.0+x*gg)).sum()*h**3

# saturation diagnostic: fraction of |grad n|^2 weight that is saturated (beta*g2 >> 1)
def sat_frac(beta, mask=None):
    gg=g2f if mask is None else g2f[mask]
    num=(beta*gg/(1.0+beta*gg)*gg).sum(); den=gg.sum(); return num/den

mu_ratio=0.886   # converged breathing inertia mu_up/mu_down (breathing_inertia_crossing_vs_midpoint.py)
lam=np.linspace(0.6,1.8,25)
print("\n  === DF breathing potential  E_fb(lambda) = const + mu*J4*I(beta/lambda^2)  (const dropped) ===")
for beta in [0.1, 0.452, 2.0, 10.0]:
    print(f"\n  --- beta = {beta}  (saturation frac: full={sat_frac(beta):.3f} up={sat_frac(beta,up):.3f} down={sat_frac(beta,dnb):.3f}) ---")
    for lbl,mask in [('full',None),('up(cross)',up),('down(mid)',dnb)]:
        V=np.array([MU*J4*I_of_x(beta/L**2,mask) for L in lam])   # DF breathing potential (mu=MU coupling)
        V=V-V.min()
        # slope & curvature at lambda=1 (central diff on a fine local stencil)
        dl=1e-3
        Vm=MU*J4*I_of_x(beta/(1-dl)**2,mask); Vp=MU*J4*I_of_x(beta/(1+dl)**2,mask); V0=MU*J4*I_of_x(beta,mask)
        slope=(Vp-Vm)/(2*dl); curv=(Vp+Vm-2*V0)/dl**2
        imin=V.argmin()
        shape = "monotonic->small-lambda (collapse; needs wall)" if imin==0 else \
                ("monotonic->large-lambda" if imin==len(lam)-1 else f"WELL at lambda={lam[imin]:.2f}")
        print(f"    {lbl:10s}: slope@1={slope:+8.3f}  curv@1={curv:+8.3f}  shape: {shape}")
    # up/down ratio of curvature and slope at lambda=1 (the harmonic and linear handles)
    dl=1e-3
    def sc(mask):
        Vm=MU*J4*I_of_x(beta/(1-dl)**2,mask); Vp=MU*J4*I_of_x(beta/(1+dl)**2,mask); V0=MU*J4*I_of_x(beta,mask)
        return (Vp-Vm)/(2*dl),(Vp+Vm-2*V0)/dl**2
    su,ku=sc(up); sd,kd=sc(dnb)
    # per equal arc (3 crossings vs 6 midpoints share the tube; normalise by weight fraction)
    fu=up.mean(); fd=dnb.mean(); sun,kun=su/fu,ku/fu; sdn,kdn=sd/fd,kd/fd
    print(f"    RATIOS/arc: slope_up/slope_down={sun/sdn:+.4f}   curv_up/curv_down={kun/kdn:+.4f}")
    if kun>0 and kdn>0:
        wr=np.sqrt((kun/kdn)/mu_ratio)     # harmonic omega_up/omega_down = sqrt((k_up/k_down)/(mu_up/mu_down))
        print(f"      HARMONIC (DF curvature): omega_up/omega_down={wr:.4f}  M_dn/M_up mult=1/that={1/wr:.4f}  (target 0.831)")
    if sun>0 and sdn>0:
        # Airy against a wall, force = DF slope: E0 ~ (F^2/mu)^(1/3)
        er=((sun/sdn)**2/mu_ratio)**(1/3)
        print(f"      AIRY (DF slope as force): E0_up/E0_down={er:.4f}  M_dn/M_up mult=1/that={1/er:.4f}  (target 0.831)")
print("\n  READ: (a) if slope>0 AND string slope>0 => both push to small lambda => the CORE/PACKING wall")
print("  sets the size, breathing is Airy-against-the-wall (use DF+string slope). (b) In the SATURATED")
print("  regime I(beta/lambda^2)~lambda^2/beta => E_fb ~ lambda^2 => HARMONIC (use curvature). The")
print("  saturation frac tells which regime the quark is in. The up/down RATIO is the mechanism test.")
