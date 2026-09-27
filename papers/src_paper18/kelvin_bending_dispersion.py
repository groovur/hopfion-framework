#!/usr/bin/env python3.11
"""
kelvin_bending_dispersion.py -- the tube BENDING (Kelvin) dispersion SHAPE, magnitude sector.

The Kelvin wave = transverse displacement of the tube core, x_c(z)=A sin(kz) (gapless bending
mode; distinct from the P18 waveguide transverse-cavity cutoff E_c which is gapped). Its energy
per unit length is  dE = (T k^2 + B k^4) A^2/4 : T = LINE TENSION (long-lambda, linear dispersion
omega=c_s k), B = BENDING STIFFNESS (short-lambda, quadratic omega ~ k^2, the classical Kelvin
regime). Framework ripples propagate at c_s=c/phi (P18), so omega^2 = c_s^2 k^2 (1 + (B/T) k^2).
KEY TEST: does the framework tube have a nonzero TENSION T (=> a LINEAR long-lambda regime, from
the confinement string), unlike a classical thin-vortex Kelvin wave (pure quadratic, T=0)?
Compute T, B from the Faddeev energy of a bent tube; report the crossover length ell=sqrt(B/T).
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=3.0-PHI
a=1.0; m_mer=1
Lxy=2.5*a; Nxy=80; Lz=8.0; Nz=64
xs=np.linspace(-Lxy,Lxy,Nxy,endpoint=False); hx=xs[1]-xs[0]
zs=np.linspace(0,Lz,Nz,endpoint=False); hz=zs[1]-zs[0]
X,Y,Z=np.meshgrid(xs,xs,zs,indexing='ij')

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def E_bent(A,k):
    xc=A*np.sin(k*Z)                       # core displaced in x, bending wavenumber k
    r=np.sqrt((X-xc)**2+Y**2); psi=np.arctan2(Y,X-xc)
    f=prof(r); Phi=m_mer*psi
    n=np.stack([np.sin(f)*np.cos(Phi),np.sin(f)*np.sin(Phi),np.cos(f)],-1)
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

print("="*70,"\n Tube BENDING (Kelvin) dispersion shape: line tension T + bending B\n","="*70,sep="")
E0=E_bent(0.0,0.0)
A=0.12*a
ks=[2*np.pi*n_/Lz for n_ in [1,2,3,4]]
print(f"  grid {Nxy}^2 x {Nz}, a={a}, A={A}, E0={E0:.1f}")
print(f"  {'k':>7} {'dE=E-E0':>12} {'dE/(A^2 Lz/4)':>14}")
lhs=[]
for k in ks:
    dE=E_bent(A,k)-E0; y=dE/(A**2*Lz/4); lhs.append(y)
    print(f"  {k:7.3f} {dE:12.3f} {y:14.3f}")
# dE/(A^2 Lz/4) = T k^2 + B k^4  -> linear fit in k^2
k2=np.array([k**2 for k in ks]); y=np.array(lhs)
B,T=np.polyfit(k2,y,1)                      # y = B*(k^2)^... wait: y=T*k2 + B*k2^2
# fit y = T*k2 + B*k2^2 :
M=np.vstack([k2,k2**2]).T
T,B=np.linalg.lstsq(M,y,rcond=None)[0]
print(f"\n  fit y = T k^2 + B k^4:   T(line tension)={T:.2f}   B(bending)={B:.2f}")
if T>0 and B>0:
    ell=np.sqrt(B/T); R=3.0+np.sqrt(2)/PHI
    print(f"  crossover length ell=sqrt(B/T)={ell:.3f} (a={a}, R=R0+r0={R:.2f}); ell/a={ell/a:.2f}")
    print(f"  => omega^2 = c_s^2 k^2 (1 + (ell k)^2), c_s=c/phi:")
    print(f"     LINEAR (sound, tension) for k<<1/ell; QUADRATIC (bending, classical-Kelvin) for k>>1/ell.")
print(f"\n  T>0 (line tension present): {'YES' if T>0 else 'NO'} -> framework Kelvin has a LINEAR")
print(f"  long-wavelength regime from the confinement string, ABSENT in classical thin-vortex Kelvin")
print(f"  (pure quadratic). That linear-vs-quadratic crossover at k~1/ell is the distinguishing test.")
