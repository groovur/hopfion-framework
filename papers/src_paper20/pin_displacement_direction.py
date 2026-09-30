#!/usr/bin/env python3.11
"""
pin_displacement_direction.py -- pin the confined-oscillation displacement direction (2026-09-29).

The Airy residual E0 ~ (F^2/mu)^{1/3} uses the string tension FORCE F = sigma*dL/dA for a displacement of
a network by A*dhat(t). The up/down ratio F_up/F_down was 1.28 (radial-xy) / 1.37 (radial-3D) -- but those
directions were CHOSEN. This pins the direction from the physics: for a tensioned tube the restoring is
sigma*kappa*Nhat (tension pulls along the curvature vector = Frenet NORMAL). Analytically, for a windowed
rigid displacement along a constant dhat,
    dL/dA = -INT w kappa (Nhat . dhat) ds     (tangent part = reparametrisation, drops to 1st order),
so the force PROJECTS the local curvature onto dhat, weighted by the window. Crossing kappa~0.40, midpoint
kappa~0.03 (15x): so F_up/F_down is HIGHLY direction-sensitive -- it runs from ~1 (dhat perp to Nhat) to
~kappa-ratio (dhat || Nhat). This script builds the Frenet frame and sweeps dhat across the transverse
{Nhat,Bhat} plane, locating radial-xy, radial-3D, Nhat, Bhat, and reports the resulting M_dn/M_up so we can
SEE whether the number is pinned by physics or is an artifact of the direction choice.
"""
import numpy as np
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI
NT=200000; t=np.linspace(0,2*np.pi,NT,endpoint=False); dt=t[1]-t[0]
def curve(tt): return np.stack([(R0+r0*np.cos(3*tt))*np.cos(2*tt),(R0+r0*np.cos(3*tt))*np.sin(2*tt),r0*np.sin(3*tt)],-1)
G=curve(t)
# Frenet frame by finite difference (periodic)
def d_ds_setup(G):
    dG=(np.roll(G,-1,0)-np.roll(G,1,0))/(2*dt)
    sp=np.linalg.norm(dG,axis=-1,keepdims=True)          # |dG/dt| = ds/dt
    That=dG/sp
    dThat=(np.roll(That,-1,0)-np.roll(That,1,0))/(2*dt)/sp   # dT/ds
    kappa=np.linalg.norm(dThat,axis=-1,keepdims=True)
    Nhat=dThat/np.clip(kappa,1e-12,None)
    Bhat=np.cross(That,Nhat)
    return That,Nhat,Bhat,sp[:,0],kappa[:,0]
That,Nhat,Bhat,sp,kappa=d_ds_setup(G)
def rhat_xy(tt): return np.stack([np.cos(2*tt),np.sin(2*tt),np.zeros_like(tt)],-1)
def rhat_3d(tt): GG=curve(tt); return GG/np.linalg.norm(GG,axis=-1,keepdims=True)
def window(tt,centers,sw=0.20):
    w=np.zeros_like(tt)
    for c in centers:
        d=np.abs(((tt-c+np.pi)%(2*np.pi))-np.pi); w+=np.exp(-(d/sw)**2)
    return w
cross_c=[np.pi/6+2*np.pi*k/3 for k in range(3)]        # z=+r0 crossings (up)
mid_c  =[k*np.pi/3 for k in range(6)]                  # z=0 midpoints (down)
wc=window(t,cross_c); wm=window(t,mid_c)

def force_num(centers,dvec):
    """numerical dL/dA for rigid displacement A*w(t)*dvec(t); central diff in A."""
    A=0.01; w=window(t,centers)[:,None]
    def L(a):
        Gd=G+a*w*dvec; dG=(np.roll(Gd,-1,0)-np.roll(Gd,1,0))/(2*dt); return (np.linalg.norm(dG,axis=-1)*dt).sum()
    return (L(A)-L(-A))/(2*A)

mu_ratio=0.886   # converged breathing inertia mu_up/mu_down (robust across dilation & radial modes)
nsite={'cross':3,'mid':6}
def report(lbl,dvec):
    fu=force_num(cross_c,dvec)/nsite['cross']; fd=force_num(mid_c,dvec)/nsite['mid']
    if abs(fd)<1e-9:
        print(f"  {lbl:26s}: F_down~0 (degenerate)"); return None
    r=fu/fd
    # tangent fraction of dvec on each network (unphysical reparam part)
    tf_u=(((dvec*That).sum(-1))**2 * wc).sum()/wc.sum(); tf_d=(((dvec*That).sum(-1))**2 * wm).sum()/wm.sum()
    er=((abs(r))**2/mu_ratio)**(1/3); mult=1/er*np.sign(r)
    print(f"  {lbl:26s}: F_up/F_down={r:+7.3f}  (tan.frac up={tf_u:.2f} dn={tf_d:.2f})  -> M_dn/M_up mult={mult:+.4f}")
    return r

print("  kappa: crossing avg={:.3f}  midpoint avg={:.3f}  (ratio {:.1f})".format(
    (kappa*wc).sum()/wc.sum(),(kappa*wm).sum()/wm.sum(),((kappa*wc).sum()/wc.sum())/((kappa*wm).sum()/wm.sum())))
print("  target M_dn/M_up mult = 0.157/0.189 = {:.4f}\n".format(0.157/0.189))
print("  === named candidate directions ===")
report('radial-xy',rhat_xy(t))
report('radial-3D (from origin)',rhat_3d(t))
report('Frenet normal Nhat',Nhat)
report('Frenet binormal Bhat',Bhat)
# transverse projections of the radial directions (drop tangent = reparam)
for lbl,dv in [('radial-xy _|_tube',rhat_xy(t)),('radial-3D _|_tube',rhat_3d(t))]:
    dp=dv-((dv*That).sum(-1,keepdims=True))*That; dp/=np.linalg.norm(dp,axis=-1,keepdims=True).clip(1e-12)
    report(lbl,dp)

print("\n  === sweep dhat = cos(psi) Nhat + sin(psi) Bhat across the transverse plane ===")
rr=[]
for psi in np.linspace(0,np.pi,13):
    dv=np.cos(psi)*Nhat+np.sin(psi)*Bhat
    fu=force_num(cross_c,dv)/nsite['cross']; fd=force_num(mid_c,dv)/nsite['mid']
    r=fu/fd if abs(fd)>1e-9 else np.nan; rr.append(r)
    mult=1/(((abs(r))**2/mu_ratio)**(1/3))*np.sign(r) if np.isfinite(r) else np.nan
    print(f"    psi={psi*180/np.pi:5.0f} deg: F_up/F_down={r:+8.3f}  M_dn/M_up mult={mult:+.4f}")
rr=np.array([x for x in rr if np.isfinite(x)])
print(f"\n  RANGE of F_up/F_down over transverse directions: [{np.nanmin(rr):+.2f}, {np.nanmax(rr):+.2f}]")
print("  READ: if the range is WIDE (F depends strongly on psi) the ~2% match is a direction ARTIFACT,")
print("  not a pinned number -> honest O3 = mechanism (Airy vs jam wall, tangent-geometry differential),")
print("  coefficient needs the jammed normal-mode problem. If NARROW & near 0.83, the direction is pinned.")
