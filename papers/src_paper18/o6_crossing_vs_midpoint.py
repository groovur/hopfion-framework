#!/usr/bin/env python3.11
"""
Option A (O6 analytic): crossing-vs-midpoint FORMATION energy on the golden trefoil.
Session 2026-09-21. Tests sec.27: up-type = CROSSING (shared) network, down-type = MIDPOINT/distal
(individual) network. Does a proper formation-energy proxy predict M_up/M_down better than Paper
XVIII's tangent-fraction (P18:prop:isospin_mass_ratio, 0.189 -> M_up/M_down=5.29 vs measured 6.38)?

Geometry (Paper XV/XVI, note sec.20): Gamma(t) = ((R0 + r0 cos3t) cos2t, (R0+r0 cos3t) sin2t,
r0 sin3t), R0=3, r0=sqrt2/phi. Verify against P18's exact tangent data (Gz=0 crossings; -0.525,
+0.321 midpt/distal; <Gz^2>=0.189).
"""
import numpy as np
PHI=(1+5**0.5)/2; R0=3.0; r0=np.sqrt(2)/PHI
alpha=1/137.036; Lam=215.5

t=np.linspace(0,2*np.pi,200001)[:-1]
c3,s3,c2,s2=np.cos(3*t),np.sin(3*t),np.cos(2*t),np.sin(2*t)
rho=R0+r0*c3
G=np.stack([rho*c2, rho*s2, r0*s3],1)                       # curve
# derivatives (analytic)
rhod=-3*r0*s3
Gd=np.stack([rhod*c2-2*rho*s2, rhod*s2+2*rho*c2, 3*r0*c3],1) # Gamma'
speed=np.linalg.norm(Gd,axis=-1)
That=Gd/speed[:,None]                                        # unit tangent
Tz=That[:,2]
# curvature kappa = |G' x G''| / |G'|^3
rhodd=-9*r0*c3
Gdd=np.stack([rhodd*c2-2*rhod*s2-2*(rhod*s2+2*rho*c2),
              rhodd*s2+2*rhod*c2+2*(rhod*c2-2*rho*s2),
              -9*r0*s3],1)
cross=np.cross(Gd,Gdd)
kappa=np.linalg.norm(cross,axis=-1)/speed**3

print("="*72,"\n(0) VERIFY geometry vs Paper XVIII P18:prop:isospin_mass_ratio\n","="*72,sep="")
# crossings: z=+r0 => sin3t=+1 (3t=pi/2). midpoints/distal: z=0 => sin3t=0 (3t=0,pi)
z=G[:,2]
def near(arr,val,tol): return np.where(np.abs(arr-val)<tol)[0]
# sample the three site categories at their exact t
def siteinfo(t3target):  # 3t = target
    tt=t3target/3.0
    i=np.argmin(np.abs(t-tt))
    return t[i], z[i], Tz[i], kappa[i]
print("  site (3t)          z        T_z       kappa")
for lbl,t3 in [('crossing (3t=pi/2)',np.pi/2),('z=0 A (3t=0)',1e-6),('z=0 B (3t=pi)',np.pi)]:
    ti,zi,tzi,ki=siteinfo(t3)
    print(f"  {lbl:20s} {zi:+.3f}   {tzi:+.4f}   {ki:.4f}")
# P18 midpoint/distal Tz = -0.525, +0.321; check which z=0 sites give these
tz_at_z0=[]
for t3 in [0,np.pi]:
    _,_,tzi,_=siteinfo(t3); tz_at_z0.append(tzi)
print(f"  z=0 T_z values: {tz_at_z0[0]:+.4f}, {tz_at_z0[1]:+.4f}  (P18: -0.525, +0.321)")
oop=0.5*(tz_at_z0[0]**2+tz_at_z0[1]**2)
print(f"  <T_z^2>_(z=0) = {oop:.4f}  (P18: 0.189)   crossing T_z^2 = {siteinfo(np.pi/2)[2]**2:.4f}")

print("\n"+"="*72,"\n(1) formation-energy PROXIES on crossing (UP) vs midpoint/distal (DOWN)\n","="*72,sep="")
# UP = crossing network (z=+r0). DOWN = the two z=0 networks (midpoint + distal).
# proxy 1: P18 in-plane tangent fraction (1 - T_z^2). up: 1-0 =1 ; down: 1-<T_z^2>=1-0.189
# proxy 2: bending energy ~ kappa^2 at the site.
# proxy 3: inter-strand proximity (sharing): min distance to a non-adjacent point on the curve.
def minstrand_dist(i):
    d=np.linalg.norm(G-G[i],axis=-1)
    # mask out points near i (within 5% of param) to exclude the local arc
    w=int(0.03*len(t))
    lo,hi=(i-w)%len(t),(i+w)%len(t)
    m=np.ones(len(t),bool)
    if lo<hi: m[lo:hi]=False
    else: m[lo:]=False; m[:hi]=False
    return d[m].min()
i_cross=np.argmin(np.abs(t-np.pi/6))       # 3t=pi/2
i_mid  =np.argmin(np.abs(t-0.0))           # 3t=0
i_dist =np.argmin(np.abs(t-np.pi/3))       # 3t=pi
k_cross=kappa[i_cross]; k_down=0.5*(kappa[i_mid]+kappa[i_dist])
d_cross=minstrand_dist(i_cross); d_down=0.5*(minstrand_dist(i_mid)+minstrand_dist(i_dist))
print(f"  curvature kappa:  crossing={k_cross:.4f}  down(mid/distal avg)={k_down:.4f}")
print(f"  inter-strand min dist: crossing={d_cross:.4f}  down={d_down:.4f}  (small=shared)")
print()
# integrated bending energy over each network's arc: partition the curve by nearest site-category.
# crossings at 3t=pi/2 (+2pi k/3); z=0 sites at 3t=0,pi (+2pi k/3). Assign each t to nearest category.
phase=(3*t)%(2*np.pi)
d_to_cross=np.minimum(np.abs(phase-np.pi/2),2*np.pi-np.abs(phase-np.pi/2))
d_to_z0=np.minimum(np.minimum(phase,2*np.pi-phase),np.abs(phase-np.pi))
is_cross=d_to_cross<d_to_z0
ds=speed*(t[1]-t[0])
Ebend_cross=np.sum((kappa**2*ds)[is_cross]); Ebend_down=np.sum((kappa**2*ds)[~is_cross])
print(f"  integrated bending E=INT kappa^2 ds:  crossing-arc={Ebend_cross:.3f}  z0-arc={Ebend_down:.3f}")

print("\n  M_up/M_down predicted by each proxy (up=crossing, down=midpt/distal):")
r1=(1-siteinfo(np.pi/2)[2]**2)/oop         # P18: M_up/M_down = in-plane|_C / <out-of-plane>_(z0) = 1/0.189
r2=k_cross/k_down                          # local curvature ratio
r2b=Ebend_cross/Ebend_down                 # integrated bending-arc ratio
r3=d_down/d_cross                          # sharing: shared(small d)=more energy => up/down ~ d_down/d_cross
for lbl,r in [('P18 out-of-plane tangent fraction (1/0.189)',r1),('local curvature kappa ratio',r2),
              ('integrated bending INT kappa^2 ds',r2b),('inter-strand proximity (sharing)',r3)]:
    print(f"    {lbl:42s}: M_up/M_down = {r:6.3f}")
# measured, two ways for the top
Mup=(2.16*1270*172760)**(1/3); Mdn=(4.67*93.4*4180)**(1/3)
Mup_f=(2.16*1270*334000)**(1/3)
print(f"  MEASURED M_up/M_down = {Mup/Mdn:.3f} (measured top) / {Mup_f/Mdn:.3f} (full-formation top 334 GeV)")
print(f"  P18 residual: measured/P18 = {(Mup/Mdn)/r1:.3f}  <- ~1.2 = the SAME uniform running factor")
print(f"  the down-sector E_8 fit carries (sec.24)? [flag, CLAUDE.md 3 -- one ratio, not claimed]")

print("\n"+"="*72,"\n(2) Q*Phi (gen1 EM) does NOT touch the triplet mean\n","="*72,sep="")
print("  Q*Phi ~ few MeV shifts only gen1 (sec.26): M_up,M_down are geometric means dominated by")
print("  c,t / s,b (GeV) -> the isospin SCALE is set by the crossing/midpoint geometry, gen1 EM is")
print("  a ~MeV correction invisible in the triplet mean. So (1) is the clean isospin-scale test.")

print("\n"+"="*72,"\n(3) VERDICT\n","="*72,sep="")
print("  Compare the proxies to measured M_up/M_down=6.38 and P18's 5.29 (tangent fraction, 21% off).")
print("  A proxy landing NEAR 6.38 with no free parameter = crossing/midpoint FORMATION energy is the")
print("  right isospin-scale object (sec.27). If a POWER is needed, report which and flag it.")
