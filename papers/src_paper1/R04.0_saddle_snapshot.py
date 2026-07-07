import numpy as np, time

phi=(1+5**0.5)/2; lam=phi**6; mu=3-phi; ep=1e-14
N=256; h=0.12; R0=4.0
r_=h*(np.arange(N)+0.5); z_=h*np.arange(N)
R,Z=np.meshgrid(r_,z_,indexing='ij')
vol=2*np.pi*2*R*h**2; vol[:,0]*=0.5
D0=(R-R0)**2+Z**2+ep

def bc(f):
    f=np.clip(f,0,np.pi); f[0,:]=0; f[-1,:]=0; f[:,-1]=0
    return f

def compute_geom(f):
    """E_geom = (J2a + mu*J2iso)*J4 -- same functional as fn_hopfion_plot.py"""
    f=np.clip(f,0,np.pi)
    fr=np.gradient(f,h,axis=0); fz=np.gradient(f,h,axis=1); fz[:,0]=0
    sf=np.sin(f); cf=np.cos(f); s2=sf**2; s3=s2*sf; s4=s2*s2
    A=1/D0+1/R**2; g2=fr**2+fz**2; fDG=fr*(R-R0)+fz*Z
    F13=-s2/D0*fDG; F12=s2/R*fr; F23=s2/R*fz
    kern=g2+s2*A
    J2a=float(np.sum(s4*kern*vol))
    J2iso=float(np.sum(kern*vol))
    J4=float(np.sum((F13**2+F12**2+F23**2)*vol))
    K=J2a+mu*J2iso; E=K*J4
    sopt=(lam*J4/K)**0.5 if K>0 else 0
    # Force
    loc_K=((4*s3*cf+mu*2*sf*cf)*kern+(s4+mu)*2*sf*cf*A)*vol
    loc_J4=2*(F13*(-2*sf*cf/D0*fDG)+F12*(2*sf*cf/R*fr)+F23*(2*sf*cf/R*fz))*vol
    fK_r=(2*s4+2*mu)*fr*vol; fK_z=(2*s4+2*mu)*fz*vol; fK_z[:,0]=0
    fJ4_r=2*(F13*(-s2/D0*(R-R0))+F12*(s2/R))*vol
    fJ4_z=2*(F13*(-s2/D0*Z)+F23*(s2/R))*vol; fJ4_z[:,0]=0
    def div(a,b): return np.gradient(a,h,axis=0)+np.gradient(b,h,axis=1)
    Force=-(( loc_K-div(fK_r,fK_z))*J4 + K*(loc_J4-div(fJ4_r,fJ4_z)))
    return J2a, J2iso, J4, E, sopt, Force

# IC: same formula as fn_hopfion_plot.py but for R0=4
f = bc(np.clip(np.pi/(1+D0/(0.6*R0)**2), 0, np.pi))

print(f"Running E_geom saddle-snapshot for R0={R0}, 10,000 steps")
print(f"(same algorithm as fn_hopfion_plot.py)")
print(f"{'step':>7} {'sopt':>8} {'J2iso/J2a':>11} {'J4/J2a':>10}")
print("-"*45)

best_dist=1e30; best_f=f.copy(); best_sopt=0; best_g=0; best_V=0
dt=1e-7; t0=time.time()

for step in range(1, 10001):
    J2a,J2iso,J4,E,sopt,Force = compute_geom(f)
    dist=abs(sopt-1.0)
    if dist<best_dist:
        best_dist=dist; best_f=f.copy()
        best_sopt=sopt; best_g=J4/J2a; best_V=J2iso/J2a
    if step%1000==0 or step==1:
        print(f"{step:7,} {sopt:8.5f} {J2iso/J2a:11.5f} {J4/J2a:10.6f}  [{time.time()-t0:.0f}s]")
    f_try=bc(f+dt*Force)
    if compute_geom(f_try)[3]<E:
        f=f_try; dt=min(dt*1.01,5e-4)
    else:
        dt*=0.7
        if dt<1e-18:
            print(f"  dt exhausted at step {step}"); break

print()
print(f"Best snapshot: sopt={best_sopt:.5f}, dist={best_dist:.5f}")
print(f"  J2iso/J2a (V at beta=0) = {best_V:.5f}")
print(f"  J4/J2a                  = {best_g:.6f}")
print(f"  Target J4/J2a (WZW)     = {2**(4/3)/phi**5:.6f}")
print(f"  Ratio                   = {best_g/(2**(4/3)/phi**5):.4f}")
print()

# Now check: what beta* does this profile need for V=phi?
from scipy.optimize import brentq
J2a_b=float(np.sum(np.sin(best_f)**4*(np.gradient(best_f,h,axis=0)**2+
    np.gradient(best_f,h,axis=1)**2+np.sin(best_f)**2*(1/D0+1/R**2))*vol))
fr2=np.gradient(best_f,h,axis=0); fz2=np.gradient(best_f,h,axis=1)
kern2=fr2**2+fz2**2+np.sin(best_f)**2*(1/D0+1/R**2)
def vres(b): return float(np.sum(kern2/(1+b*kern2)*vol))/J2a_b - phi
try:
    b_star=brentq(vres,0.001,5.0)
    print(f"  beta* for V=phi on this profile: {b_star:.5f}")
    print(f"  (cf. R0=3 beta*=0.452)")
except Exception as e:
    print(f"  beta* search: {e}")

np.save('./f_R04.0_egeom.npy', best_f)
print(f"\nProfile saved: /tmp/f_R04_egeom.npy")
