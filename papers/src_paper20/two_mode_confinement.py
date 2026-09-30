#!/usr/bin/env python3.11
"""
two_mode_confinement.py -- the confined quark's TWO transverse modes, N-hat + B-hat (2026-09-29).

Lead (user): the ~2.5% residual (pinned N-hat mode gives 0.153 vs 0.157) is a missing COUNTER-force. The
transverse oscillation is 2D: the Frenet-NORMAL N-hat is restored by string TENSION (linear/Airy, force
~sigma*kappa, ratio 1.375 -> pushes M_dn/M_up DOWN), the BINORMAL B-hat has NO tension restoring (dL/dA=0)
and is restored by the Faddeev BENDING stiffness (harmonic, ratio from local_bending -> pushes M_dn/M_up UP).
The down network contains a NEAR-INFLECTION (inner-midpoint kappa~0.026, tau~11.3) where the tension force
nearly vanishes and the bending mode dominates -- the counter-force, localised. Total confined zero-point per
quark = E_Nhat(Airy) + E_Bhat(harmonic). Single knob: sigma (tension), cutoff-set (Faddeev E is scale-invariant
so tension is NOT derivable from it). This computes all geometric + field ingredients, forms E_N(sigma),
E_B(sigma), finds the sigma that lands M_dn/M_up on 0.157, and reports whether the required tension/bending
balance is physical (order-1) rather than tuned.
Units: code units, hbar=1. Airy ground state (linear V=F|x| against wall): E0 = 2.338*(F^2/(2 mu))^(1/3).
Harmonic: E0 = (1/2) sqrt(k/mu).
"""
import numpy as np, sys, time
from scipy.spatial import KDTree
from scipy.optimize import brentq
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI; C_star=2.5062; MU=3.0-PHI
N=int(sys.argv[1]) if len(sys.argv)>1 else 64
h=(2*(R0+r0+1/C_star+0.6))/N
cv=h*(np.arange(N)-N//2+0.5); pts=np.stack(np.meshgrid(cv,cv,cv,indexing='ij'),-1).reshape(-1,3)
NTf=20000; t_frame,_,N1f,N2f,H=build_compensated_frame_arclength(NT=NTf)
def frame_at_t(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NTf; return N1f[idx],N2f[idx]
def base_curve(t): return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),(R0+r0*np.cos(3*t))*np.sin(2*t),r0*np.sin(3*t)],-1)
# --- Frenet frame T,N,B and curvature on a fine grid; nearest-index lookup ---
NF=200000; tf=np.linspace(0,2*np.pi,NF,endpoint=False); dtf=tf[1]-tf[0]; Gf=base_curve(tf)
dG=(np.roll(Gf,-1,0)-np.roll(Gf,1,0))/(2*dtf); spf=np.linalg.norm(dG,axis=-1,keepdims=True); Tf=dG/spf
dTf=(np.roll(Tf,-1,0)-np.roll(Tf,1,0))/(2*dtf)/spf; kapf=np.linalg.norm(dTf,axis=-1,keepdims=True)
Nf=dTf/np.clip(kapf,1e-12,None); Bf=np.cross(Tf,Nf)
def lut(arr,t): return arr[(np.searchsorted(tf,t%(2*np.pi)))%NF]
def Nhat(t): return lut(Nf,t)
def Bhat(t): return lut(Bf,t)
cross_c=[np.pi/6+2*np.pi*k/3 for k in range(3)]; mid_c=[k*np.pi/3 for k in range(6)]
def window(t,centers,sw=0.20):
    w=np.zeros_like(t)
    for c in centers:
        d=np.abs(((t-c+np.pi)%(2*np.pi))-np.pi); w+=np.exp(-(d/sw)**2)
    return w

# ---------- GEOMETRY (pure curve): g_N = INT w kappa ds ; h_B,h_N = d2L/dA2 ----------
NTg=40000; tg=np.linspace(0,2*np.pi,NTg,endpoint=False); dtg=tg[1]-tg[0]; Gg=base_curve(tg)
kap_g=lut(kapf[:,0],tg); sp_g=lut(spf[:,0],tg)      # [:,0] -> 1D (kapf,spf are (NF,1) from keepdims)
def gN(centers,ns): return (window(tg,centers)*kap_g*sp_g).sum()*dtg/ns    # INT w kappa ds  per site
def d2L(centers,dfun,ns,A=0.01):
    w=window(tg,centers)[:,None]
    def L(a):
        Gd=Gg+a*w*dfun(tg); dGa=(np.roll(Gd,-1,0)-np.roll(Gd,1,0))/(2*dtg); return (np.linalg.norm(dGa,axis=-1)*dtg).sum()
    return (L(A)+L(-A)-2*L(0))/A**2/ns
gN_u=gN(cross_c,3); gN_d=gN(mid_c,6)
hB_u=d2L(cross_c,Bhat,3); hB_d=d2L(mid_c,Bhat,6)
hN_u=d2L(cross_c,Nhat,3); hN_d=d2L(mid_c,Nhat,6)
print(f"  GEOMETRY/site: gN(INT w k ds) up={gN_u:.4f} dn={gN_d:.4f} (ratio {gN_u/gN_d:.3f})")
print(f"                 d2L/dA2 |Bhat  up={hB_u:.4f} dn={hB_d:.4f} ;  |Nhat up={hN_u:.4f} dn={hN_d:.4f}")

# ---------- FIELD: Faddeev bending stiffness k_Fadd and inertia mu along a direction ----------
NT=4000; t_arr=np.linspace(0,2*np.pi,NT,endpoint=False)
def energy_and_dir(A,centers,dfun):
    disp=base_curve(t_arr)+(A*window(t_arr,centers))[:,None]*dfun(t_arr)
    lobe=[np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in [0,2*np.pi/3,4*np.pi/3]]
    trees=[KDTree(disp[li]) for li in lobe]; tls=[t_arr[li] for li in lobe]
    dp,tp=[],[]
    for tr,tl in zip(trees,tls):
        d,i=tr.query(pts,workers=-1); dp.append(d); tp.append(tl[i])
    ds=np.stack(dp,1); ts=np.stack(tp,1); o=np.argsort(ds,1)
    t1=np.take_along_axis(ts,o,1)[:,0]; d1=np.take_along_axis(ds,o,1)[:,0]
    t2=np.take_along_axis(ts,o,1)[:,1]; d2=np.take_along_axis(ds,o,1)[:,1]
    def dcurve(tq): return base_curve(tq)+(A*window(tq,centers))[:,None]*dfun(tq)
    def chi(tq):
        n1,n2=frame_at_t(tq); rel=pts-dcurve(tq); return np.arctan2(np.sum(rel*n2,1),np.sum(rel*n1,1))
    P1=chi(t1)+3*t1; P2=chi(t2)+3*t2; r1=np.clip(d1,1e-6,None); r2=np.clip(d2,1e-6,None)
    f0=lambda r:2*np.arctan(np.maximum(r,1e-9)**(-C_star)); f1=f0(r1*C_star); f2=f0(r2*C_star)
    w1=1/r1**2; w2=1/r2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*P1)+w2*np.sin(f2/2)*np.exp(1j*P2)
    m=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=m; z2/=m
    n=np.stack([2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),np.abs(z1)**2-np.abs(z2)**2],-1)
    n/=np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)
    ng=n.reshape(N,N,N,3); nx,ny,nz=ng[...,0],ng[...,1],ng[...,2]
    def cd(u,ax): return (np.roll(u,-1,ax)-np.roll(u,1,ax))/(2*h)
    dd=[[cd(c,a) for a in range(3)] for c in (nx,ny,nz)]; g2=sum(dd[c][a]**2 for c in range(3) for a in range(3))
    (nxx,nxy,nxz),(nyx,nyy,nyz),(nzx,nzy,nzz)=dd
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4=(Fxy**2+Fxz**2+Fyz**2).sum()*h**3; J2a=(((1-nz**2).clip(0,1)**2)*g2).sum()*h**3; J2iso=g2.sum()*h**3
    E=(J2a+MU*J2iso)*J4
    return E, n

def k_and_mu(centers,dfun,ns,A=0.05):
    E0,n0=energy_and_dir(0.0,centers,dfun); Ep,npf=energy_and_dir(A,centers,dfun); Em,nmf=energy_and_dir(-A,centers,dfun)
    k=(Ep+Em-2*E0)/A**2/ns
    dn=(npf-nmf)/(2*A); dn=dn-(dn*n0).sum(-1,keepdims=True)*n0; mu=(dn**2).sum()*h**3/ns
    return k,mu

t0=time.time()
kB_u,muB_u=k_and_mu(cross_c,Bhat,3); kB_d,muB_d=k_and_mu(mid_c,Bhat,6)
kN_u,muN_u=k_and_mu(cross_c,Nhat,3); kN_d,muN_d=k_and_mu(mid_c,Nhat,6)
print(f"  FIELD (built {time.time()-t0:.1f}s): k_Fadd|Bhat/site up={kB_u:.1f} dn={kB_d:.1f} (ratio {kB_u/kB_d:.3f})  mu|Bhat up={muB_u:.2f} dn={muB_d:.2f}")
print(f"                                    k_Fadd|Nhat/site up={kN_u:.1f} dn={kN_d:.1f}  mu|Nhat up={muN_u:.2f} dn={muN_d:.2f}")

# ---------- ASSEMBLE the two-mode zero-point as a function of sigma (tension, code units) ----------
def E_N(sigma,gN_,mu_):    # Airy: linear force F=sigma*gN against the jam wall
    F=sigma*gN_; return 2.338*(F**2/(2*mu_))**(1/3)
def E_B(sigma,hB_,kF_,mu_): # harmonic: k = sigma*d2L/dA2|B + Faddeev bending
    k=sigma*hB_+kF_; return 0.5*np.sqrt(max(k,1e-12)/mu_)
def corr(sigma):
    Eu=E_N(sigma,gN_u,muN_u)+E_B(sigma,hB_u,kB_u,muB_u)
    Ed=E_N(sigma,gN_d,muN_d)+E_B(sigma,hB_d,kB_d,muB_d)
    return Ed/Eu
target=0.157/0.189
print(f"\n  target correction = {target:.4f}  (M_dn/M_up=0.157);  pure-Nhat corr={E_N(1,gN_d,muN_d)/E_N(1,gN_u,muN_u):.4f}")
print(f"  sigma ->  correction  ->  M_dn/M_up   [r=E_B/E_N up]   [k_B decomp: sigma*hB vs k_Fadd]")
for sg in [0.5,1,2,5,10,20,50,100,200,500]:
    r=E_B(sg,hB_u,kB_u,muB_u)/E_N(sg,gN_u,muN_u)
    print(f"   sigma={sg:6.1f}:  {corr(sg):.4f}   {0.189*corr(sg):.4f}      r={r:.4f}    sigma*hB={sg*hB_u:.1f} vs kFadd={kB_u:.1f}")
try:
    sg_star=brentq(lambda s: corr(s)-target, 0.01, 1e5)
    r_star=E_B(sg_star,hB_u,kB_u,muB_u)/E_N(sg_star,gN_u,muN_u)
    print(f"\n  sigma* (lands on 0.157) = {sg_star:.3f} code units;  bending admixture r={r_star:.4f}")
    print(f"  at sigma*: tension part sigma*hB_up={sg_star*hB_u:.2f} vs Faddeev bending kFadd_up={kB_u:.2f}"
          f"  (ratio {sg_star*hB_u/kB_u:.3f})")
    # tension vs bending magnitude sanity: Airy force sigma*gN vs Faddeev bending k as an energy over length R~ (R0)
    print(f"  N-hat linear force sigma*gN_up={sg_star*gN_u:.2f};  is sigma* ~ order of the field stiffness scale?")
except Exception as e:
    print("  no root:",e)
