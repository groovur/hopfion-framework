#!/usr/bin/env python3.11
"""
cable_twist_chiral.py -- does CHIRALITY split the two twist directions (isospin from twist)?

Test of the user's lead (2026-09-28): the framing-twist sine-Gordon has kink vs antikink
(twist two ways) = isospin. A PURE (achiral) cable gives V=V0(1+cos dtheta), EVEN in dtheta
=> kink/antikink degenerate = isospin-symmetric (cable_twist_potential.py confirmed this, 2D
z-uniform, achiral by construction). The (2,3)-cable is CHIRAL: its two strands wind
helically around each other. This script builds that 3D helical cable with winding rate tau
(sign = chirality) and asks: does V(dtheta) develop an ODD (sin) component?
  - odd component present, flipping sign with tau  => chirality SPLITS the twist directions:
    isospin from twist, mass-splitting sign set by the cable chirality (the writhe). DERIVED.
  - stays even for both tau  => the splitting is NOT elastic; it needs a parity-violating term
    (weak/chiral condensate) and the sign remains a formation selection. Clean negative.
PARITY CHECK: mirror sends (dtheta, tau) -> (-dtheta, -tau), so E(dtheta,tau)=E(-dtheta,-tau).
Hence at fixed tau an odd-in-dtheta part IS allowed and MUST reverse with tau -- that is the
discriminator computed here.

Run:  python3.11 cable_twist_chiral.py [--N 48 --Nz 32]
"""
import numpy as np, argparse
PHI=(1+5**0.5)/2; MU=3.0-PHI
ap=argparse.ArgumentParser()
ap.add_argument('--N',type=int,default=48)      # xy grid
ap.add_argument('--Nz',type=int,default=32)     # z grid (one helix pitch, periodic)
ap.add_argument('--a',type=float,default=1.0)   # tube radius
ap.add_argument('--dfac',type=float,default=1.6)# strand separation d = dfac*a
ap.add_argument('--ndt',type=int,default=25)
args=ap.parse_args()
a=args.a; d=args.dfac*a; Rc=d/2
L=3.0*a; N=args.N; Nz=args.Nz
xs=np.linspace(-L,L,N,endpoint=False); h=xs[1]-xs[0]
# one full helical pitch over z in [0,Pz); pitch chosen so the winding is O(1) over ~knot scale
Pz=2*np.pi          # tau=+-1 -> one turn per 2*pi in z (representative); scan-independent of Pz choice for the ODD-vs-EVEN question
zs=np.linspace(0,Pz,Nz,endpoint=False); hz=zs[1]-zs[0]
X,Y,Z=np.meshgrid(xs,xs,zs,indexing='ij')

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def strand(cx,cy,theta):
    r=np.sqrt((X-cx)**2+(Y-cy)**2); psi=np.arctan2(Y-cy,X-cx)
    f=prof(r); Phi=psi+theta            # m_mer=1 meridian winding + framing theta
    w=1.0/np.clip(r,1e-3,None)**2
    return w,f,Phi

def director(q,tau):
    # helical inter-strand winding (CHIRAL background, sign=tau) + a relative-twist KINK that
    # WINDS along z with charge q (=+1 kink / -1 antikink): phi(z)=q*2*pi*z/Pz. The chiral
    # coupling tau*d_z(phi) is nonzero only for a winding twist -> this is what a uniform
    # dtheta could not see.
    ang=tau*Z
    phi=q*2*np.pi*Z/Pz
    c1x,c1y= Rc*np.cos(ang), Rc*np.sin(ang)
    c2x,c2y=-Rc*np.cos(ang),-Rc*np.sin(ang)
    w1,f1,P1=strand(c1x,c1y,+phi/2)
    w2,f2,P2=strand(c2x,c2y,-phi/2)
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    mag=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=mag; z2/=mag
    nx=2*np.real(np.conj(z1)*z2); ny=2*np.imag(np.conj(z1)*z2); nz=np.abs(z1)**2-np.abs(z2)**2
    return np.stack([nx,ny,nz],-1)

def energy(q,tau):
    n=director(q,tau); nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,ax,hx): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*hx)
    nxx,nxy,nxz=cd(nx,0,h),cd(nx,1,h),cd(nx,2,hz)
    nyx,nyy,nyz=cd(ny,0,h),cd(ny,1,h),cd(ny,2,hz)
    nzx,nzy,nzz=cd(nz,0,h),cd(nz,1,h),cd(nz,2,hz)
    g2=nxx**2+nxy**2+nxz**2+nyx**2+nyy**2+nyz**2+nzx**2+nzy**2+nzz**2
    s4=(1-nz**2).clip(0,1)**2; dv=h*h*hz
    K=(s4*g2).sum()*dv+MU*g2.sum()*dv
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*dv
    return K*J4

if __name__=="__main__":
    print("="*72,f"\n chiral cable: kink(q=+1) vs antikink(q=-1) energy on helical (tau) background\n","="*72,sep="")
    print(f"  N={N} Nz={Nz} d={d}a   (split E(+q)-E(-q) at fixed tau = the isospin splitting)")
    print(f"  {'tau':>6} {'E(kink q=+1)':>14} {'E(antikink q=-1)':>16} {'split':>12} {'rel':>10}")
    E={}
    for tau in [0.0,+1.0,-1.0]:
        for q in [+1,-1]:
            E[(tau,q)]=energy(q,tau)
    for tau in [0.0,+1.0,-1.0]:
        ep,em=E[(tau,+1)],E[(tau,-1)]; split=ep-em
        print(f"  {tau:>6.1f} {ep:>14.2f} {em:>16.2f} {split:>12.3e} {split/((ep+em)/2):>10.2e}")
    print()
    # parity checks (must hold): E(+tau,+q)=E(-tau,-q)
    p1=abs(E[(+1.0,+1)]-E[(-1.0,-1)]); p2=abs(E[(+1.0,-1)]-E[(-1.0,+1)])
    print(f"  parity check E(+tau,+q)=E(-tau,-q): |diff|={p1:.2e}, {p2:.2e} (should be ~0)")
    split_tau1=E[(+1.0,+1)]-E[(+1.0,-1)]; split_tau0=E[(0.0,+1)]-E[(0.0,-1)]
    scale=(E[(+1.0,+1)]+E[(+1.0,-1)])/2
    print(f"  kink/antikink split at tau=0 (achiral): {split_tau0:.3e}")
    print(f"  kink/antikink split at tau=+1 (chiral): {split_tau1:.3e}  (rel {split_tau1/scale:.2e})")
    if abs(split_tau0)<1e-4*scale and abs(split_tau1)>10*abs(split_tau0)+1e-3*abs(scale):
        print("  => VERDICT: CHIRALITY SPLITS kink vs antikink (elastic isospin-from-twist).")
        print("     Sign of the split set by the cable chirality (writhe). Compare rel-split to 0.157/0.189.")
    elif abs(split_tau1)<1e-3*scale:
        print("  => VERDICT: NO elastic split even on the chiral cable -> the isospin splitting is NOT")
        print("     in the elastic/Faddeev sector; it needs a parity-violating term (weak/chiral")
        print("     condensate) and the sign remains a formation selection. Clean negative.")
    else:
        print("  => VERDICT: ambiguous; refine resolution.")
