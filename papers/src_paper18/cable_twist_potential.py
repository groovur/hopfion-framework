#!/usr/bin/env python3.11
"""
cable_twist_potential.py -- the MISSING sine-Gordon potential V(Delta-theta) for the
2-strand cable relative twist (the "remaining step" flagged in perstrand_twist_energy.py).

MACHINERY (not masses): the tube framing theta(z,t) obeys chi d2t theta - K d2z theta
+ V'(theta) = 0 (sine-Gordon), c_s=sqrt(K/chi). On a SINGLE tube V is flat (framing =
Goldstone). The periodic potential comes from the CABLE: two director strands whose
overlap energy depends on their RELATIVE framing Delta-theta = theta_1 - theta_2.
This script computes E(Delta-theta) via the framework's S^3-lift blend + Faddeev energy
(same as gradient_flow_constrained), on a 2D cross-section (z-uniform -> per unit length),
and asks: is it periodic ~ (1 - cos N*Delta-theta)?  If yes -> sine-Gordon, and we read
off the potential amplitude and period N.
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=3.0-PHI
a=1.0                          # tube radius
d=1.6*a                        # cable strand separation (overlapping tubes)
m_mer=1                        # meridian winding per strand
L=3.0*a; N=128                 # 2D cross-section grid
xs=np.linspace(-L,L,N,endpoint=False); h=xs[1]-xs[0]
X,Y=np.meshgrid(xs,xs,indexing='ij')

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def strand(cx, theta):
    r=np.sqrt((X-cx)**2+Y**2); psi=np.arctan2(Y,X-cx)
    f=prof(r); Phi=m_mer*psi+theta
    w=1.0/np.clip(r,1e-3,None)**2
    return w, f, Phi

def director(dtheta):
    w1,f1,P1=strand(-d/2, +dtheta/2)
    w2,f2,P2=strand(+d/2, -dtheta/2)
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    mag=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=mag; z2/=mag
    nx=2*np.real(np.conj(z1)*z2); ny=2*np.imag(np.conj(z1)*z2); nz=np.abs(z1)**2-np.abs(z2)**2
    return np.stack([nx,ny,nz],-1)

def energy(dtheta):
    n=director(dtheta); nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
    # 2D (x,y) derivatives; z-uniform so d/dz = 0
    nxx,nxy=cd(nx,0),cd(nx,1); nyx,nyy=cd(ny,0),cd(ny,1); nzx,nzy=cd(nz,0),cd(nz,1)
    g2=nxx**2+nxy**2+nyx**2+nyy**2+nzx**2+nzy**2
    s4=(1-nz**2).clip(0,1)**2
    J2a=(s4*g2).sum()*h*h; J2iso=g2.sum()*h*h
    K=J2a+MU*J2iso
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    J4=(Fxy**2).sum()*h*h
    return K, J4, K*J4

print("="*70,"\n Cable relative-twist potential V(Delta-theta)  [2-strand, per unit length]\n","="*70,sep="")
dths=np.linspace(0,2*np.pi,25)
Ks,J4s,Es=[],[],[]
for dt in dths:
    K,J4,E=energy(dt); Ks.append(K); J4s.append(J4); Es.append(E)
Ks,J4s,Es=map(np.array,(Ks,J4s,Es))
print(f"  d={d}a  m={m_mer}  grid {N}^2  h={h:.3f}")
print(f"  {'dtheta':>7} {'K':>10} {'J4':>10} {'E=K*J4':>12}")
for dt,K,J,E in zip(dths,Ks,J4s,Es):
    print(f"  {dt:7.3f} {K:10.3f} {J:10.3f} {E:12.2f}")

# Fourier decomposition V(dtheta) = c0 + sum_k c_k cos(k dtheta) (symmetric -> cosines only)
V=Es-Es.mean()
print("\n  Fourier content of V (harmonic amplitudes c_k, cos(k*dtheta)):")
for k in [1,2,3,4]:
    ck=2*np.mean(V*np.cos(k*dths))
    print(f"    k={k}: c_k={ck:9.2f}  ({'DOMINANT' if k==1 else f'{abs(ck/ (2*np.mean(V*np.cos(dths)))+1e-9)*100:.1f}% of k=1'})")
c1=2*np.mean(V*np.cos(dths)); recon=c1*np.cos(dths)
resid=np.sqrt(np.mean((V-recon)**2))/(np.abs(V).max()+1e-9)
print(f"  pure k=1 (single cosine) rel.resid = {resid:.4f}  -> vacuum at dtheta=pi (anti-aligned strands)")
print(f"  V(dtheta) = V0(1+cos dtheta),  V0 = {(Es.max()-Es.min())/2:.1f} per unit length")

print("\n"+"="*70,"\n  robustness: potential amplitude V0 vs strand separation d\n","="*70,sep="")
print(f"  {'d/a':>5} {'V0(peak-peak/2)':>16} {'k=1 resid':>10}")
for dfac in [1.2,1.4,1.6,1.8,2.0,2.4]:
    globals()['d']=dfac*a
    Es2=np.array([energy(dt)[2] for dt in dths]); V2=Es2-Es2.mean()
    r2=np.sqrt(np.mean((V2-(2*np.mean(V2*np.cos(dths)))*np.cos(dths))**2))/(np.abs(V2).max()+1e-9)
    print(f"  {dfac:5.1f} {(Es2.max()-Es2.min())/2:16.1f} {r2:10.4f}")
print("  -> V0 should FALL with separation (overlap-driven); a clean k=1 across d = robust sine-Gordon.")
