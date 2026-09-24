#!/usr/bin/env python3.11
"""
Pin the E_6-native absolute tower with Q*Phi (linear EM) folded in, keeping the T(2,n) geometry
in view. Session 2026-09-21. Follows sec.24 (E_6 phase pi/4), sec.26 (EM linear -Q*Phi),
Paper XVIII T(2,n) tower + Jones moduli, confinement note sec.6.

STEP 1: EM-subtract -> topological masses m_top = m_meas + Q*Phi. Pin Phi by requiring the
        gen1 TOPOLOGICAL ordering to be up-heavier (consistent with the topological tower being
        up-heavier at all gens; sec.25/26).
STEP 2: test the E_6-native tower (phase pi/4, twist |Tw|={1,2,5}, top-transition anchor) on the
        topological up masses. Report honestly what closes and what does not.
STEP 3: tabulate the T(2,n) Jones tower (geometry the user asked to keep in mind).
"""
import numpy as np
PHI=(1+5**0.5)/2; LN=np.log(PHI); PI=np.pi; alpha=1/137.036
m={'u':2.16,'d':4.67,'s':93.4,'c':1270.,'b':4180.,'t':172760.}
Q={'u':2/3,'c':2/3,'t':2/3,'d':-1/3,'s':-1/3,'b':-1/3}
Tw={'u':1,'c':2,'t':5,'d':1,'s':2,'b':5}   # |E_6 per-strand twist| = generation (note sec.13)
Lam=215.5

print("="*72,"\nSTEP 1: fold Q*Phi in (EM-subtract), pin Phi by gen1 topological ordering\n","="*72,sep="")
print(" m_top = m_meas + Q*Phi (since m_meas = m_top - Q*Phi).  up gains +2/3 Phi, down loses 1/3 Phi.")
for Phi in [alpha*Lam, 2.5, 3.77]:
    mt={q:m[q]+Q[q]*Phi for q in m}
    flip = mt['u']>mt['d']
    print(f"  Phi={Phi:5.2f} MeV: m_top(u)={mt['u']:.2f} m_top(d)={mt['d']:.2f} -> "
          f"gen1 {'UP heavier (topological ordering consistent)' if flip else 'down still heavier'}")
# Phi that exactly equalises gen1 (u=d) then anything above flips it:
Phi_eq=(m['d']-m['u'])/(Q['u']-Q['d'])   # m_u+Q_u Phi = m_d+Q_d Phi
print(f"  => gen1 flips to up-heavier for Phi > {Phi_eq:.2f} MeV. alpha*Lam={alpha*Lam:.2f};"
      f" O(1) factor {Phi_eq/(alpha*Lam):.1f} -> Phi~a few MeV. Take Phi={Phi_eq:.2f} (gen1 marginal).")
Phi=Phi_eq*1.2   # slightly above threshold so gen1 is up-heavier
mt={q:m[q]+Q[q]*Phi for q in m}
print(f"  chosen Phi={Phi:.2f} MeV -> topological masses (MeV):",
      {q:round(mt[q],2) for q in ['u','d','s','c','b','t']})

print("\n"+"="*72,"\nSTEP 2: E_6-native tower on TOPOLOGICAL up masses (phase pi/4, top anchor)\n","="*72,sep="")
print(" base formula m ~ exp(-2pi Q_group T_g), Q_group=3; up anchored at top-transition (gen3).")
# (a) fixed tower level n, generation entirely in T_g: solve T_g each up quark needs
print(" (a) if all up at fixed n, generation in T_g: required (T_g - T_top):")
for q in ['c','u']:
    dT = -np.log(mt[q]/mt['t'])/(6*PI)
    print(f"     {q}: T-T_top = {dT:+.3f}  (twist gap from top |Tw_t|-|Tw|={5-Tw[q]})")
rc=-np.log(mt['c']/mt['t'])/(6*PI); ru=-np.log(mt['u']/mt['t'])/(6*PI)
print(f"     ratio (T_c-T_top)/(T_u-T_top) = {rc/ru:.3f};  twist-gap ratio (5-2)/(5-1)=3/4=0.75")
print(f"     -> {'MATCH' if abs(rc/ru-0.75)<0.05 else 'NO MATCH'}: T_g is NOT linear in the twist gap (n-T degeneracy).")
# (b) mass linear in |Tw| or exponent? (test on topological masses)
print(" (b) is ln(m_top) linear in |Tw|={1,2,5}? up:")
lnm=[np.log(mt[q]) for q in ['u','c','t']]; tw=[1,2,5]
s1=(lnm[1]-lnm[0])/(tw[1]-tw[0]); s2=(lnm[2]-lnm[1])/(tw[2]-tw[1])
print(f"     slope u->c (per |Tw|) = {s1:.2f};  c->t = {s2:.2f}  -> {'linear' if abs(s1-s2)<0.3 else 'NOT linear'}")
print(" => the absolute up tower does NOT reduce to a clean exp(k|Tw|) or single T_g law.")
print("    generation STEEPNESS works (E_6 phase 2.25, sec.24) but ABSOLUTE u,c do not close here.")

print("\n"+"="*72,"\nSTEP 3: the T(2,n) geometry (kept in view)\n","="*72,sep="")
print(" T(2,n) tower (P18): fused ODD n = golden/confined; unfused EVEN n = rational.")
print("   n=3 trefoil  Q_H=3 baryon   |J(q5)|=1/phi")
print("   n=4 unfused                 |J|=1")
print("   n=5 cinquefoil (Q_H=5)      |J|=phi   <- next fused sector, energy-input pushes here")
print("   n=6 unfused                 |J|=-1  (period 10, mirror |J^n|=|J^{10-n}|)")
print(" GENERATIONS = per-strand twist {1,2,5} of the FIXED T(2,3) (3 crossings=colour). The higher")
print(" T(2,n) (cinquefoil) are higher Q_H SECTORS, not generations. User's 'energy input -> higher")
print(" T' = transient access to Q_H=5 during (de)construction; 'shared lobes/oscillating' = the JAM")
print(" (sec.19) + confined-quark oscillation (confinement-note sec.3). Bearing on mass: the up-sector")
print(" formation is NOT a clean isolated trefoil -> no clean isolated tower should fit -> O6 (shared/")
print(" jammed formation Delta E), not an isolated-saddle mass, is the right object.")
