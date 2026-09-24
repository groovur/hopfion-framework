#!/usr/bin/env python3.11
"""
O6: director S^2 solid-angle / geometry across the Q_H=2 -> 3 step.
Session 2026-09-22. Does a geometric quantity of the director field jump by ~4pi (or the mass-ratio
C~12.9) from Q_H=2 (lepton) to Q_H=3 (quark)? Standard hopfion ansatz (stereographic S^3->S^2),
Hopf charge verified by FFT. Look for a LINEAR ~4pi (the factor sec.38 said the CS/vac-pol 4pi's canNOT
supply). Honest test: if nothing is ~4pi/~12.9, that confirms sec.38 (no geometric 4pi factor either).
"""
import numpy as np
PI=np.pi
N=int(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 64
L=4.0
x=np.linspace(-L,L,N,endpoint=False); h=x[1]-x[0]
X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r2=X**2+Y**2+Z**2
# stereographic R^3 -> S^3 subset C^2
Z1=(2*(X+1j*Y))/(1+r2); Z2=(2*Z+1j*(r2-1))/(1+r2)

def director(A,B):
    # w = Z1^A / Z2^B : rational map S^3 -> CP^1 = S^2 ; Hopf charge ~ A*B
    w=(Z1**A)/(Z2**B + 1e-12)
    aw2=np.abs(w)**2
    n=np.stack([2*np.real(w)/(1+aw2), 2*np.imag(w)/(1+aw2), (aw2-1)/(1+aw2)],-1)
    nrm=np.linalg.norm(n,axis=-1,keepdims=True); return n/np.clip(nrm,1e-12,None)

def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
def geom(n):
    nx,ny,nz=n[...,0],n[...,1],n[...,2]
    d=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]
    g2=sum(d[c][a]**2 for c in range(3) for a in range(3))
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=d
    # F_i = n . (d_j n x d_k n) : the pullback of the S^2 area 2-form (Fx=F_yz etc.)
    Fx=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)   # yz
    Fy=nx*(nyz*nzx-nzz*nyx)+ny*(nzz*nxx-nxz*nzx)+nz*(nxz*nyx-nyz*nxx)   # zx
    Fz=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)   # xy
    F=np.stack([Fx,Fy,Fz],-1)
    J2=g2.sum()*h**3
    J4=(Fx**2+Fy**2+Fz**2).sum()*h**3
    solidangle=np.linalg.norm(F,axis=-1).sum()*h**3   # total |F| = S^2 area-form magnitude integrated
    # Hopf charge via FFT: A = curl^{-1} F ; Q = (1/(4pi)^2) INT A.F  (calibrated to charge-1 below)
    k=2*PI*np.fft.fftfreq(N,d=h)
    KX,KY,KZ=np.meshgrid(k,k,k,indexing='ij'); K2=KX**2+KY**2+KZ**2; K2[0,0,0]=1
    Fh=[np.fft.fftn(F[...,i]) for i in range(3)]
    # A_hat = i k x F_hat / k^2
    Ah=[1j*(KY*Fh[2]-KZ*Fh[1])/K2, 1j*(KZ*Fh[0]-KX*Fh[2])/K2, 1j*(KX*Fh[1]-KY*Fh[0])/K2]
    A=[np.real(np.fft.ifftn(Ah[i])) for i in range(3)]
    AF=(A[0]*Fx+A[1]*Fy+A[2]*Fz).sum()*h**3
    return J2,J4,solidangle,AF

print(f"N={N} L={L} h={h:.3f}")
print(f"{'(A,B)':>7} {'Q(FFT/calib)':>12} {'J2':>10} {'J4':>10} {'solidAng':>10} {'AF_raw':>12}")
res={}
AF1=None
for (A,B) in [(1,1),(1,2),(1,3)]:
    n=director(A,B); J2,J4,sa,AF=geom(n)
    if AF1 is None: AF1=AF   # calibrate Q=1 to (1,1)
    Q=AF/AF1
    res[(A,B)]=(Q,J2,J4,sa,AF)
    print(f"  ({A},{B}) {Q:12.3f} {J2:10.2f} {J4:10.2f} {sa:10.2f} {AF:12.2f}")

print("\n=== across the Q_H=2 -> 3 step (using (1,2) and (1,3)) ===")
q2=res[(1,2)]; q3=res[(1,3)]
for lbl,i in [('J2',1),('J4',2),('solidAngle',3)]:
    print(f"  {lbl:11s}: Q2={q2[i]:.2f}  Q3={q3[i]:.2f}   ratio Q3/Q2={q3[i]/q2[i]:.3f}   diff={q3[i]-q2[i]:.2f}")
print(f"  targets: 4pi={4*PI:.3f}, C~12.9 ; look for any ratio/diff ~ these.")
print("  (magnitudes are grid/box-dependent; RATIOS across the step are the physical content.)")
