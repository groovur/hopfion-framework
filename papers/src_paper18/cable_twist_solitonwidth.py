#!/usr/bin/env python3.11
"""
cable_twist_solitonwidth.py -- pin the framing sine-Gordon's soliton width and check
the sector.  The framing is a DIRECTOR (orientation) mode -> c_dir=sqrt(K/chi)->0, NOT
the magnitude c/phi.  So its speed and gap are both chi-suppressed, but their RATIO is
chi-independent: gap/c_dir = 1/xi, xi = sqrt(K_eff/V0) = the twist-soliton width.
Compute K_eff (z-gradient stiffness of the relative twist) and V0 (potential amplitude)
ON THE SAME 3D cable, form xi, compare to the tube radius a and knot scale R=R0+r0.
This tells us the twist-quantum size; and since c_framing=c_dir != c/phi, whether it
touches the P18 waveguide cutoff E_c~pi c_s/R (which uses c_s=c/phi) at all.
"""
import numpy as np
PHI=(1+5**0.5)/2; MU=3.0-PHI
a=1.0; d=1.6*a; m_mer=1
Lxy=3.0*a; Nxy=96; Lz=6.0; Nz=24
xs=np.linspace(-Lxy,Lxy,Nxy,endpoint=False); hx=xs[1]-xs[0]
zs=np.linspace(0,Lz,Nz,endpoint=False); hz=zs[1]-zs[0]
X,Y,Z=np.meshgrid(xs,xs,zs,indexing='ij')

def prof(r):
    u=np.clip(r/a,0,1); return np.pi*(1+np.cos(np.pi*u))/2*(r<a)

def strand(cx, theta):
    r=np.sqrt((X-cx)**2+Y**2); psi=np.arctan2(Y,X-cx)
    f=prof(r); w=1.0/np.clip(r,1e-3,None)**2
    return w, f, m_mer*psi+theta

def director(dtheta_field):
    w1,f1,P1=strand(-d/2, +dtheta_field/2)
    w2,f2,P2=strand(+d/2, -dtheta_field/2)
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

# --- V0: uniform relative twist, potential well depth (per unit length) ---
Evac = E_geom(director(np.pi+0*Z))          # vacuum Dtheta=pi
Etop = E_geom(director(0*Z))                # barrier Dtheta=0
V0 = (Etop-Evac)/(2*Lz)                     # V=V0(1+cos); peak-to-peak=2V0, over length Lz
print("="*70,"\n Framing sine-Gordon: soliton width and sector\n","="*70,sep="")
print(f"  grid {Nxy}^2 x {Nz},  d={d}a,  Lz={Lz}")
print(f"  V0 (potential amplitude, per unit length) = {V0:.3f}   [Evac={Evac:.1f} Etop={Etop:.1f}]")

# --- K_eff: small z-gradient of the relative twist about the vacuum ---
print(f"\n  z-gradient stiffness (Dtheta=pi + q*2pi z/Lz), is E-Evac ~ q^2?")
print(f"  {'q':>5} {'E-Evac':>12} {'/q^2':>10}")
qs=[0.25,0.5,1.0,1.5,2.0]; coeffs=[]
for q in qs:
    Eq=E_geom(director(np.pi + q*2*np.pi*Z/Lz))
    dE=Eq-Evac; coeffs.append(dE/q**2)
    print(f"  {q:5.2f} {dE:12.3f} {dE/q**2:10.3f}")
# E-Evac = K_eff * integral (dDtheta/dz)^2 dz = K_eff*(2pi q/Lz)^2 * Lz = K_eff*(2pi)^2 q^2/Lz
slope=np.polyfit(np.log(qs),np.log(np.clip([E_geom(director(np.pi+q*2*np.pi*Z/Lz))-Evac for q in qs],1e-9,None)),1)[0]
K_eff=np.median(coeffs)*Lz/(2*np.pi)**2
print(f"  log-log slope E-Evac ~ q^{slope:.2f}  (2.0 = clean quadratic gradient)")
print(f"  K_eff (grad coeff of (dDtheta/dz)^2, per unit length) = {K_eff:.3f}")

xi=np.sqrt(K_eff/max(V0,1e-9))
print(f"\n  SOLITON WIDTH  xi = sqrt(K_eff/V0) = {xi:.3f}  (tube radius a={a}, so xi/a={xi/a:.2f})")
R=3.0+np.sqrt(2)/PHI
print(f"  vs tube radius a={a} and knot scale R=R0+r0={R:.2f}:  xi/a={xi/a:.2f}, xi/R={xi/R:.3f}")
print("\n  SECTOR: c_framing = c_dir = sqrt(K/chi) -> 0 (director/orientation), NOT c/phi.")
print("  P18 waveguide cutoff E_c ~ pi c_s/R uses c_s=c/phi (magnitude). Different sector:")
print("  framing SG = cold director-twist (gap & speed both chi-suppressed, ratio=1/xi);")
print("  Kelvin/waveguide ripples = fast magnitude sector. The SG gap is NOT E_c.")
