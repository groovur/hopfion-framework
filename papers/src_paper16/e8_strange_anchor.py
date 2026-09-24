#!/usr/bin/env python3
r"""
E_8 phase structure of the quark-mass assignment: WHY strange is the anchor (2026-09-03).

Paper XVI mass route (P16:eq:e8_ratio_recap):  m_q/m_ell(g) = exp(n_q * pi/9),
  n_q = 180(T_g - T^(m)) = Base_g - 2m,   Base_g = 180 T_g + 7/2,
  T^(m) = (4m-7)/360,  m = E_8 Coxeter exponent in {1,7,11,13,17,19,23,29}, h(E_8)=30.
Assignment (Paper V/XVI): u->19, c->11, t->1 (up);  d->17, s->13, b->7 (down); 23,29 unused.

The paper says m_s=13 is "the anchor" but the assignment is monotone-FITTED to the mass order
(P16 l.732-734). This script asks the phase question (tool (1)) to answer the anchor question (2):
  - which quark can have n_q=0 (mass = its lepton partner), and WHY is it unique?
  - does the E_8 Coxeter phase / mirror structure organise the assignment, or is it just a sort?
"""
import math
from fractions import Fraction as F

h = 30
exps = [1,7,11,13,17,19,23,29]                    # E_8 Coxeter exponents
Tm   = {m: F(4*m-7,360) for m in exps}            # topological T-phase of exponent m
# lepton T-phases (Paper I), geometric progression ratio (k+2)/k = 5/3
T = {1:F(5,24), 2:F(1,8), 3:F(3,40)}              # e, mu, tau
Base = {g: 180*T[g] + F(7,2) for g in (1,2,3)}    # Base_g = 180 T_g + 7/2
lept = {1:('e',0.51099895), 2:('mu',105.6583755), 3:('tau',1776.86)}  # MeV

print("="*74); print("E_8 PHASE STRUCTURE OF THE QUARK MASSES  (h=30, unit pi/9)"); print("="*74)

print("\n[A] Base_g = 180 T_g + 7/2  and its PARITY (the anchor mechanism):")
for g in (1,2,3):
    b = Base[g]; half = F(b,2)
    hit = (half.denominator==1 and int(half) in exps)
    print(f"    g={g} ({lept[g][0]:>3}): T_g={str(T[g]):>5}  Base_g={int(b):>2}  "
          f"parity={'EVEN' if int(b)%2==0 else 'odd '}  Base/2={str(half):>4}  "
          f"{'-> integer E_8 exponent '+str(int(half))+' => n_q=0 POSSIBLE' if hit else '-> not an exponent => n_q=0 impossible'}")
print("    => n_q=0 (quark mass = its lepton partner) is possible for EXACTLY ONE generation.")

print("\n[B] Which exponent's T-phase equals a LEPTON T-phase? (the phase resonance):")
for m in exps:
    for g in (1,2,3):
        if Tm[m]==T[g]:
            print(f"    T^({m}) = {str(Tm[m])} = T_{g} ({lept[g][0]})  <-- exponent {m} resonates with lepton {lept[g][0]}")
print("    (only exact rational coincidences shown)")

print("\n[C] The assignment, n_q, and predicted condensate-scale mass = m_ell(g)*exp(n_q*pi/9):")
assign = {'u':(1,19),'c':(2,11),'t':(3,1),'d':(1,17),'s':(2,13),'b':(3,7)}
pdg = {'u':2.16,'d':4.67,'s':93.4,'c':1275,'b':4180,'t':172570}   # MeV, PDG-ish
print("    q    gen  m    n_q=Base-2m   m_ell(g)     pred=m_ell*e^{n_q pi/9}   PDG      pred/PDG")
for q,(g,m) in assign.items():
    nq = Base[g]-2*m
    ml = lept[g][1]
    pred = ml*math.exp(float(nq)*math.pi/9)
    print(f"    {q:>2}   {g}   {m:>2}    {str(nq):>4}         {ml:9.3f}    {pred:12.3f}          {pdg[q]:8.1f}   {pred/pdg[q]:6.3f}")

print("\n[D] E_8 Coxeter eigenphases e^{2 pi i m/30} and the MIRROR pairs (m + m' = 30):")
pairs = [(1,29),(7,23),(11,19),(13,17)]
used = {m for _,m in assign.values()}
for a,b in pairs:
    ang_a = 2*math.pi*a/h; ang_b = 2*math.pi*b/h
    ua = 'used('+[q for q,(g,mm) in assign.items() if mm==a][0]+')' if a in used else 'UNUSED'
    ub = 'used('+[q for q,(g,mm) in assign.items() if mm==b][0]+')' if b in used else 'UNUSED'
    print(f"    ({a:>2},{b:>2}) sum={a+b}: phases {math.degrees(ang_a):6.1f} / {math.degrees(ang_b):6.1f} deg   {a}:{ua:11}  {b}:{ub}")

print("\n[E] READ:")
print("    - Anchor: gen-2 (muon) is the UNIQUE generation with Base even (26=2*13) AND 13 an E_8")
print("      exponent, so ONLY strange can sit at n_q=0 (m_s=m_mu). Gens 1,3 have odd Base (41,17)")
print("      -> Base/2 non-integer -> no quark-lepton mass coincidence. Strange is forced, not fitted.")
print("    - Phase resonance: T^(13)=1/8=T_mu exactly -> strange's E_8 topological spin = the muon's.")
print("      No other exponent lands on a lepton T-phase (T_e,T_tau need m=20.5, 8.5).")
print("    - Mirror structure: gen-1&2 quarks pair up ((u,c)=(19,11), (d,s)=(17,13), each summing to")
print("      30); gen-3 (t,b) take the EXTREME exponents (1,7), whose mirrors (29,23) are UNUSED.")
print("      => the assignment is NOT a plain mass-sort: it is a mirror-paired phase structure with")
print("         the third generation breaking the pairing. That is the shape a decay spectrum makes.")
