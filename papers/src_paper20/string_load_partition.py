#!/usr/bin/env python3.11
"""
string_load_partition.py -- the CONFINEMENT (string) stiffness half of the O3 residual (2026-09-29).

The breathing/scale mode is Faddeev-scale-flat; its restoring force is the STRING V=sigma*L(config).
For a confined quark the string dominates the restoring, so the oscillation stiffness k ~ sigma*d2L/dA2
(L = tube/arc length; A = radial displacement amplitude of a network). This is PURE CURVE GEOMETRY:
displace the crossing(up) vs midpoint(down) network radially (confinement dir rhat) and read how much the
arc length curves -> the string stiffness per network. sigma cancels in the up/down RATIO.
Combined with the (already-computed, converged) breathing inertia ratio mu_up/mu_down=0.886:
    omega_up/omega_down = sqrt( (k_string_up/k_string_down) / (mu_up/mu_down) )
    M_dn/M_up multiplier = 1/that.  target 0.83  =>  needs k_string_up/k_string_down ~ 1.28.
"""
import numpy as np
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI
NT=200000; t=np.linspace(0,2*np.pi,NT,endpoint=False); dt=t[1]-t[0]
def curve(tt): return np.stack([(R0+r0*np.cos(3*tt))*np.cos(2*tt),(R0+r0*np.cos(3*tt))*np.sin(2*tt),r0*np.sin(3*tt)],-1)
def rhat_xy(tt): return np.stack([np.cos(2*tt),np.sin(2*tt),np.zeros_like(tt)],-1)      # radial in xy
def rhat_3d(tt):                                                                         # full radial from origin
    G=curve(tt); return G/np.linalg.norm(G,axis=-1,keepdims=True)
def window(tt,centers,sw=0.20):
    w=np.zeros_like(tt)
    for c in centers:
        d=np.abs(((tt-c+np.pi)%(2*np.pi))-np.pi); w+=np.exp(-(d/sw)**2)
    return w
cross_c=[np.pi/6+2*np.pi*k/3 for k in range(3)]        # z=+r0 crossings (up)
mid_c  =[k*np.pi/3 for k in range(6)]                  # z=0 midpoints (down)

def arclen(A,centers,rfun):
    G=curve(t)+(A*window(t,centers))[:,None]*rfun(t)
    dG=(np.roll(G,-1,0)-np.roll(G,1,0))/(2*dt)
    return (np.linalg.norm(dG,axis=-1)*dt).sum()

def stiffness(centers,rfun,A=0.02):
    L0=arclen(0.0,centers,rfun); Lp=arclen(A,centers,rfun); Lm=arclen(-A,centers,rfun)
    return (Lp+Lm-2*L0)/A**2, (Lp-Lm)/(2*A)            # d2L/dA2 (stiffness), dL/dA (tension force)

mu_ratio=0.886
for lbl,rfun in [('radial-xy',rhat_xy),('radial-3D (from origin)',rhat_3d)]:
    ku,fu=stiffness(cross_c,rfun); kd,fd=stiffness(mid_c,rfun)
    ku_s=ku/3; kd_s=kd/6                                # per site (3 crossings, 6 midpoints)
    r=ku_s/kd_s
    print(f"=== displacement dir: {lbl} ===")
    print(f"  crossing(up):  d2L/dA2/site = {ku_s:8.3f}   dL/dA = {fu/3:8.3f}")
    print(f"  midpoint(dn):  d2L/dA2/site = {kd_s:8.3f}   dL/dA = {fd/6:8.3f}")
    print(f"  k_string_up/k_string_down = {r:.4f}   (need ~1.28)")
    if r>0 and mu_ratio>0:
        wu_wd=np.sqrt(r/mu_ratio); print(f"  -> omega_up/omega_down=sqrt(r/mu_ratio)={wu_wd:.4f}  M_dn/M_up mult=1/that={1/wu_wd:.4f}  (target 0.83)")
    print()
print("READ: if k_string_up/k_string_down ~ 1.28 (crossing stiffer against the string), then combined")
print("with the breathing inertia it drives M_dn/M_up -> ~0.83 (0.189 -> 0.157). The junction/arc geometry")
print("attaches the string at the CROSSINGS, so up>down is expected; magnitude is the test.")
