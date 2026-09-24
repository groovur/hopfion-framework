#!/usr/bin/env python3.11
"""
Per-strand twist Faddeev-energy SCALING (the honest version of sec.30; user 2026-09-21).
Question: does the Faddeev energy of a twisted director tube scale ~Tw^2 (quadratic only) or does the
QUARTIC term give a ~Tw^4 piece? If ~Tw^4, then sec.29's mass~|Tw|^4.2 is mass~energy (LINEAR), and the
retracted mass~energy^2 / phi^3=sqrt(phi^6) story is NOT needed.

Model: a straight periodic director tube n=(sin f cos Phi, sin f sin Phi, cos f), Phi = m*psi + q*(2pi z/L),
f(r) smooth pi->0. Two windings: m = MERIDIAN (director wraps the tube cross-section), q = LONGITUDINAL
(framing/twist along the tube). Compute K=(J2a+mu J2iso) and J4 (quartic) via the SAME E_geom finite
differences as gradient_flow_constrained.py. Report scaling of K, J4 vs each winding.
CAVEAT: on a single tube, varying q changes the GLOBAL framing (SL) -- the true PER-STRAND generation
twist (colour SL=3 fixed) is a 2-strand-cable relative twist; this straight-tube calc gives the LOCAL
energy SCALINGS, which is what fixes the mass-twist power. Construction of the SL-preserving cable twist
is the remaining step.
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=1.0
a=1.0                      # tube radius
N=64; Nz=64; L=6.0         # box: [-2a,2a]^2 x [0,L], periodic z
xs=np.linspace(-2*a,2*a,N,endpoint=False); zs=np.linspace(0,L,Nz,endpoint=False)
hx=xs[1]-xs[0]; hz=zs[1]-zs[0]
X,Y,Z=np.meshgrid(xs,xs,zs,indexing='ij')
r=np.sqrt(X**2+Y**2); psi=np.arctan2(Y,X)
def prof(r):                # f(0)=pi (n=-z core), f(a)=0 (n=+z outside): smooth, compact
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)
f=prof(r); sinf=np.sin(f); cosf=np.cos(f)

def energy(m,q):
    Phi=m*psi+q*(2*np.pi*Z/L)
    n=np.stack([sinf*np.cos(Phi), sinf*np.sin(Phi), cosf],-1)
    nrm=np.linalg.norm(n,axis=-1,keepdims=True); n=n/np.clip(nrm,1e-12,None)
    nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,ax,hh): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*hh)
    hs=[hx,hx,hz]
    d=[[cd(c,ax,hs[ax]) for ax in range(3)] for c in (nx,ny,nz)]
    g2=sum(d[c][ax]**2 for c in range(3) for ax in range(3))
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=d
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    dv=hx*hx*hz
    J2iso=g2.sum()*dv
    s4=(1-nz**2).clip(0,1)**2; J2a=(s4*g2).sum()*dv
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*dv
    return J2a+MU*J2iso, J4

def scaling(label, configs, xvals):
    Ks,J4s=[],[]
    for (m,q) in configs:
        K,J4=energy(m,q); Ks.append(K); J4s.append(J4)
    Ks=np.array(Ks); J4s=np.array(J4s); x=np.array(xvals,float)
    # log-log slope over the nonzero-x points
    msk=x>0
    def slope(y):
        yy=np.array(y)[msk]; xx=x[msk]
        return np.polyfit(np.log(xx),np.log(np.clip(yy,1e-30,None)),1)[0] if msk.sum()>1 else np.nan
    print(f"\n{label}")
    print(f"  {'param':>6} {'K':>12} {'J4':>12}")
    for xi,K,J in zip(xvals,Ks,J4s): print(f"  {xi:6.1f} {K:12.4f} {J:12.4f}")
    print(f"  log-log slope:  K ~ x^{slope(Ks):.2f}   J4 ~ x^{slope(J4s):.2f}   "
          f"E=K*J4 ~ x^{slope(Ks)+slope(J4s):.2f}")

print("="*66,"\nHow does the twisted-tube Faddeev energy scale? (E_geom = K*J4)\n","="*66,sep="")
scaling("(A) LONGITUDINAL winding q (framing/twist), m=1 fixed:",
        [(1,q) for q in [1,2,3,4,5,6]], [1,2,3,4,5,6])
scaling("(B) MERIDIAN winding m, q=1 fixed:",
        [(m,1) for m in [1,2,3,4,5,6]], [1,2,3,4,5,6])
scaling("(C) BOTH m=q=T (twist raises both):",
        [(T,T) for T in [1,2,3,4,5,6]], [1,2,3,4,5,6])
print("\n"+"="*66,"\nREAD\n","="*66,sep="")
print("  If any realization gives E=K*J4 ~ T^~4, then sec.29's mass~|Tw|^4.2 = mass~energy (LINEAR),")
print("  and the retracted mass~energy^2/phi^3 story is unnecessary. If only ~T^2, the mapping is")
print("  elsewhere. (Absolute magnitudes are profile/grid-dependent; only the SCALING powers matter.)")
