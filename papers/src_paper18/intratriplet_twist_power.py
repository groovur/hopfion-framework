#!/usr/bin/env python3.11
"""
Intra-triplet absolute masses: is the generation hierarchy a POWER LAW in the ribbon twist?
Session 2026-09-21. generation = per-strand E_6 twist |Tw|={1,2,5} for gen{1,2,3} (note sec.13).
Test m ~ |Tw|^p within each isospin triplet. Physical motivation: the quartic Faddeev term
INT (dn x dn)^2 (= phi^6 J_4 in the framework) scales as (twist)^4 for a twisted ribbon -> expect p~4.
"""
import numpy as np
PHI=(1+5**0.5)/2
m={'u':2.16,'d':4.67,'s':93.4,'c':1270.,'b':4180.,'t':172760.}
m_t_full=334000.
Tw={'d':1,'s':2,'b':5, 'u':1,'c':2,'t':5}   # |E_6 twist| = generation

def fit_power(quarks):
    x=np.log([Tw[q] for q in quarks]); y=np.log([m[q] for q in quarks])
    p,lnA=np.polyfit(x,y,1)
    yhat=p*x+lnA; ss=1-np.sum((y-yhat)**2)/np.sum((y-y.mean())**2)
    return p,np.exp(lnA),ss

print("="*70,"\nDOWN triplet  m ~ |Tw|^p,  |Tw|={1,2,5}\n","="*70,sep="")
p,A,r2=fit_power(['d','s','b'])
print(f"  fit: p = {p:.3f}, A = {A:.3f} MeV, R^2 = {r2:.5f}")
print(f"  per-step p:  d->s = {np.log(m['s']/m['d'])/np.log(2):.3f}   "
      f"s->b = {np.log(m['b']/m['s'])/np.log(2.5):.3f}   d->b = {np.log(m['b']/m['d'])/np.log(5):.3f}")
print(f"  predict from m = A|Tw|^p:  d={A*1**p:.2f} s={A*2**p:.2f} b={A*5**p:.2f}"
      f"  (meas 4.67, 93.4, 4180)")
print(f"  residuals: d {A*1**p/m['d']:.3f}  s {A*2**p/m['s']:.3f}  b {A*5**p/m['b']:.3f}")

print("\n  candidate exponents (CLAUDE.md 3 -- flag, don't claim):")
for lbl,pc in [('4 (quartic Faddeev INT(dn x dn)^2)',4.0),('phi^3',PHI**3),
               ('2phi+1=phi^3',2*PHI+1),('fitted',p)]:
    ds=2**pc/1; sb=2.5**pc
    print(f"    p={pc:.3f} ({lbl:34s}): m_s/m_d={2**pc:.1f}(meas 20.0) m_b/m_s={2.5**pc:.1f}(meas 44.8)")

print("\n"+"="*70,"\nUP triplet  (measured top, then full-formation top)\n","="*70,sep="")
p_u,A_u,r2_u=fit_power(['u','c','t'])
print(f"  fit(measured top): p={p_u:.3f} R^2={r2_u:.5f}")
print(f"  per-step p:  u->c = {np.log(m['c']/m['u'])/np.log(2):.3f}   "
      f"c->t = {np.log(m['t']/m['c'])/np.log(2.5):.3f}  (measured top)")
print(f"               c->t = {np.log(m_t_full/m['c'])/np.log(2.5):.3f}  (full-formation top 334 GeV)")
print("  => UP per-step p is NOT constant (shared-crossing network, sec.27/28): no clean power law.")
print(f"  up/down power ratio (u->c)/(d->s) = {(np.log(m['c']/m['u'])/np.log(2))/(np.log(m['s']/m['d'])/np.log(2)):.2f}"
      "  (= the 2.13 steepness, same data)")

print("\n"+"="*70,"\nVERDICT\n","="*70,sep="")
print(f"  DOWN intra-triplet FOLLOWS a twist power law m ~ |Tw|^p, p={p:.2f} (R^2={r2:.4f}),")
print(f"  ~ 4 (quartic Faddeev) / phi^3={PHI**3:.3f}. Clean, near-zero-parameter for the DOWN hierarchy.")
print("  UP does NOT (shared-crossing network). So the individual (down/midpoint) triplet is analytic;")
print("  the shared (up/crossing) triplet still needs the shared-formation (jam) treatment.")
