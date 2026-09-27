#!/usr/bin/env python3.11
"""
cable_twist_kinetic.py -- nail the sine-Gordon GRADIENT (kinetic) term of the cable relative
twist, fixing the earlier bug: cable_twist_solitonwidth.py used NON-integer q with periodic z,
so Dtheta didn't close over the period (boundary discontinuity -> the spurious q^0.9). Here q is
INTEGER (periodic) and Lz is large, so dDtheta/dz=2pi q/Lz spans a small-gradient range.

Analytic expectation: E=K*J4 with a z-gradient adding I2 (dDtheta/dz)^2 to both K and J4 ->
E-E0 = K_eff (2pi q/Lz)^2 * Lz = K_eff (2pi)^2 q^2 / Lz  -> QUADRATIC in q (standard SG gradient).
Extract K_eff, then the soliton width xi=sqrt(K_eff/V0) with a CLEAN K_eff, and re-test xi~R.
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=3.0-PHI; a=1.0; d=1.6*a; m_mer=1
Lxy=3.0*a; Nxy=88; Lz=16.0; Nz=96
xs=np.linspace(-Lxy,Lxy,Nxy,endpoint=False); hx=xs[1]-xs[0]
zs=np.linspace(0,Lz,Nz,endpoint=False); hz=zs[1]-zs[0]
X,Y,Z=np.meshgrid(xs,xs,zs,indexing='ij')

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def director(dth):                              # dth = relative-twist field Dtheta(z)
    def strand(cx,th):
        r=np.sqrt((X-cx)**2+Y**2); psi=np.arctan2(Y,X-cx)
        return 1.0/np.clip(r,1e-3,None)**2, prof(r), m_mer*psi+th
    w1,f1,P1=strand(-d/2,+dth/2); w2,f2,P2=strand(+d/2,-dth/2)
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    mag=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=mag; z2/=mag
    return np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),
                     np.abs(z1)**2-np.abs(z2)**2],-1)

def E_geom(n):
    nx,ny,nz=n[...,0],n[...,1],n[...,2]; hs=[hx,hx,hz]
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

print("="*70,"\n Cable twist PHONON dispersion: small oscillation about the vacuum\n","="*70,sep="")
E0=E_geom(director(np.pi+0*Z))                  # vacuum Dtheta=pi
eps=0.10                                         # small oscillation amplitude
print(f"  grid {Nxy}^2 x {Nz}, Lz={Lz}, d={d}a, eps={eps};  E0={E0:.1f}")
print(f"  Dtheta = pi + eps*sin(2pi n z/Lz):  E-E0 = (K_eff k^2 + V'') eps^2 Lz/4")
print(f"  {'n':>3} {'k':>7} {'E-E0':>11} {'y=(E-E0)/(eps^2 Lz/4)':>22}")
ns=[1,2,3,4]; ys=[]; k2=[]
for n_ in ns:
    k=2*np.pi*n_/Lz
    dE=E_geom(director(np.pi + eps*np.sin(k*Z)))-E0
    y=dE/(eps**2*Lz/4); ys.append(y); k2.append(k**2)
    print(f"  {n_:>3} {k:7.3f} {dE:11.2f} {y:22.2f}")
# y = V'' + K_eff k^2  -> linear fit in k^2
K_eff,Vpp=np.polyfit(k2,ys,1)                    # slope=K_eff, intercept=V''
print(f"\n  linear fit y = V'' + K_eff k^2:")
print(f"    K_eff (gradient/kinetic coeff) = {K_eff:.2f}")
print(f"    V'' (potential curvature at vacuum = mass^2 gap) = {Vpp:.2f}")
# cross-check V'' against the uniform-twist potential amplitude V0 (V=V0(1+cosDth) -> V''(pi)=V0)
Etop=E_geom(director(0*Z)); V0=(Etop-E0)/(2*Lz)
print(f"    cross-check: V0 from uniform-twist scan = {V0:.2f}  (should ~ V'')  ratio V''/V0={Vpp/max(V0,1e-9):.2f}")
xi=np.sqrt(K_eff/max(Vpp,1e-9)); R=3.0+np.sqrt(2)/PHI
print(f"\n  soliton width xi = sqrt(K_eff/V'') = {xi:.3f}   xi/a={xi/a:.2f}  xi/R={xi/R:.3f} (R={R:.2f})")
print("\n  READ: a CLEAN linear y-vs-k^2 => STANDARD sine-Gordon (quadratic gradient K_eff k^2 +")
print("  massive mode V''); V''~V0 confirms the potential; xi=sqrt(K_eff/V'') the twist-soliton width.")