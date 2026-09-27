#!/usr/bin/env python3.11
"""
cstar_selfconsistency.py -- the 1-day C* check: does the JAMMED (density-feedback) Q_H=3 saddle
prefer a FATTER tube than the isolated C*=2.5062, i.e. is it resolvable at moderate N?

Tube radius = 1/C* (profile f=2 arctan((rho C*)^{-C*}), transition at rho=1/C*). The density
feedback S_eff ~ 1/(1+beta*rho) SOFTENS the isotropic stiffness -> should SPREAD the director ->
FATTEN the tube. Test: scan C*, compute the Paper I self-consistency virial V=J2iso_fb/J2a (target
phi) and the Derrick sopt=sqrt(phi^6 J4/K_fb) (target 1), with feedback (beta=0.452) and without.
The self-consistent C* is where V=phi. If that C* is small (fat tube, radius~1), it's resolvable at
N~150 with existing solvers -> the expensive tube-adapted build is UNNECESSARY.

CAVEAT: on a moderate grid the THIN C* end is under-resolved (its V,sopt are unreliable); the FAT
end is accurate. So a V=phi crossing at FAT C* is trustworthy; run --N 96 and 128 to check the
crossing is resolution-stable.

Usage (smoke): python3.11 cstar_selfconsistency.py --N 48 --Cstar_list "1.0,1.5,2.0,2.5" --beta 0.452
Real:         python3.11 cstar_selfconsistency.py --N 128 --Cstar_list "0.8,1.0,1.2,1.5,1.8,2.2,2.6,3.0" --beta 0.452
"""
import numpy as np, argparse
from scipy.spatial import KDTree
PHI=(1+5**0.5)/2; MU=3.0-PHI; PHI6=PHI**6
ap=argparse.ArgumentParser()
ap.add_argument('--N',type=int,default=96)
ap.add_argument('--box',type=float,default=11.2)
ap.add_argument('--R0',type=float,default=3.0)
ap.add_argument('--r0',type=float,default=np.sqrt(2)/PHI)
ap.add_argument('--Cstar_list',type=str,default="1.0,1.5,2.0,2.5")
ap.add_argument('--beta',type=float,default=0.452)
ap.add_argument('--NT',type=int,default=6000)
args=ap.parse_args()
N=args.N; h=args.box/N
cv=h*(np.arange(N)-N//2+0.5)
pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)
t=np.linspace(0,2*np.pi,args.NT,endpoint=False)
R0,r0=args.R0,args.r0
G=np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)],1)
dist,idx=KDTree(G).query(pts); tn=t[idx]; rho=np.clip(dist,1e-6,None)

def observables(Cstar,beta):
    f=2*np.arctan((rho*Cstar)**(-Cstar)); th=3*tn
    n=np.stack([np.sin(f)*np.cos(th),np.sin(f)*np.sin(th),np.cos(f)],-1).reshape(N,N,N,3)
    n/=np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
    nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,a): return (np.roll(u,-1,a)-np.roll(u,1,a))/(2*h)
    d=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]
    g2=sum(d[c][a]**2 for c in range(3) for a in range(3))
    s4=(1-nz**2).clip(0,1)**2; dv=h**3
    J2a=(s4*g2).sum()*dv
    J2iso_fb=(g2/(1+beta*g2)).sum()*dv
    K_fb=J2a+MU*J2iso_fb
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=d
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*dv
    V=J2iso_fb/J2a if J2a>1e-12 else 0.0
    sopt=(PHI6*J4/K_fb)**0.5 if K_fb>1e-12 else 0.0
    return V,sopt,K_fb*J4,J4

Cs=[float(x) for x in args.Cstar_list.split(',')]
print("="*76,f"\n C* self-consistency scan  (N={N}, h={h:.3f}, beta={args.beta})\n","="*76,sep="")
print(f"  tube radius=1/C*; pts across tube radius at this N = (1/C*)/h")
print(f"  {'C*':>5} {'radius=1/C*':>11} {'pts/rad':>8} {'V=J2iso_fb/J2a':>14} {'(phi)':>7} {'sopt':>7} {'E=K_fb J4':>11}")
for C in Cs:
    V,sopt,E,J4=observables(C,args.beta)
    print(f"  {C:>5.2f} {1/C:>11.3f} {(1/C)/h:>8.1f} {V:>14.4f} {PHI:>7.4f} {sopt:>7.3f} {E:>11.1f}")
print(f"\n  Self-consistent C* = where V crosses phi={PHI:.4f}. Then tube radius=1/C*; check pts/radius")
print(f"  at N=150 (h={args.box/150:.3f}): resolvable if (1/C*)/{args.box/150:.3f} >~ 8-10.")
print(f"  If V=phi lands at FAT C* (~1, radius~1, well-resolved here) -> jammed tube resolvable at")
print(f"  N~150, expensive tube-adapted build UNNECESSARY. If only at thin C*~2.5 (under-resolved")
print(f"  here) -> re-run at N=128 to confirm; thin tube -> build needed.")
