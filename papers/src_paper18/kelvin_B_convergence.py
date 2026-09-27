#!/usr/bin/env python3.11
"""
kelvin_B_convergence.py -- resolve the bending stiffness B via resolution convergence.
The coarse run gave dE/k^2 drifting DOWN ~4% (fit B<0), consistent with a finite-difference
artifact: central differences underestimate high-k gradients by ~(kh)^2, which fakes a
negative k^4 term. Test: fit y(k)=dE/(A^2 Lz/4)=T k^2 + B k^4 at increasing resolution.
  * B -> 0 as h->0  => tube is genuinely PURE TENSION (linear), crossover ell=sqrt(B/T) sub-tube.
  * B -> positive const => real bending stiffness; ell is the linear->quadratic crossover length.
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=3.0-PHI; a=1.0; m_mer=1; Lxy=2.5*a; Lz=8.0; A=0.10*a

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def E_bent(A,k,Nxy,Nz):
    xs=np.linspace(-Lxy,Lxy,Nxy,endpoint=False); hx=xs[1]-xs[0]
    zs=np.linspace(0,Lz,Nz,endpoint=False); hz=zs[1]-zs[0]
    X,Y,Zc=np.meshgrid(xs,xs,zs,indexing='ij')
    xc=A*np.sin(k*Zc); r=np.sqrt((X-xc)**2+Y**2); psi=np.arctan2(Y,X-xc)
    f=prof(r); Phi=m_mer*psi
    nx,ny,nz=np.sin(f)*np.cos(Phi),np.sin(f)*np.sin(Phi),np.cos(f)
    hs=[hx,hx,hz]
    def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*hs[ax])
    d_=[[cd(c,ax) for ax in range(3)] for c in (nx,ny,nz)]
    g2=sum(d_[c][ax]**2 for c in range(3) for ax in range(3))
    s4=(1-nz**2).clip(0,1)**2; dv=hx*hx*hz
    K=(s4*g2).sum()*dv+MU*g2.sum()*dv
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=d_
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*dv
    return K*J4

ks=[2*np.pi*n_/Lz for n_ in [1,2,3,4]]
print("="*72,"\n Bending stiffness B: resolution convergence (T k^2 + B k^4 fit)\n","="*72,sep="")
print(f"  {'Nxy':>4} {'Nz':>4} {'hx':>6} {'T':>10} {'B':>9} {'B/T':>9} {'ell=sqrt(B/T)':>13}")
for Nxy,Nz in [(80,64),(112,96),(144,128),(180,160)]:
    hx=2*Lxy/Nxy
    E0=E_bent(0,0,Nxy,Nz)
    y=np.array([(E_bent(A,k,Nxy,Nz)-E0)/(A**2*Lz/4) for k in ks])
    k2=np.array([k**2 for k in ks])
    T,B=np.linalg.lstsq(np.vstack([k2,k2**2]).T,y,rcond=None)[0]
    ell=np.sqrt(B/T) if (B>0 and T>0) else float('nan')
    print(f"  {Nxy:>4} {Nz:>4} {hx:6.3f} {T:10.1f} {B:9.1f} {B/T:9.4f} {ell:13.3f}")
print("\n  READ: if B climbs toward 0 (from negative) as h shrinks -> discretization artifact,")
print("  tube is PURE TENSION (linear Kelvin) at all physical wavelengths. If B settles positive")
print("  -> real bending; ell/a is the linear->quadratic (framework->classical-Kelvin) crossover.")