#!/usr/bin/env python3.11
"""
dark_parity_z2.py -- IS THE D-TYPE DARK PARITY DYNAMICAL?  (the open leg from notes (yyy):
Paper XXI's absolute stability rests on the dicyclic Z_2 being a CONSERVED DYNAMICAL symmetry,
not merely an outer automorphism -- note (f) flagged this as unproven.)

METHOD (reuses the validated CFT machinery already in the repo):
  * h-values / central charges: so2r_1(r), su2k(k) from dtype_cft_tower_engine.py (src_paper21).
  * the DYNAMICAL test is the SAME one Paper XV/XVII use for the neutrino's Majorana Z_2
    (paper15_simple_current.py): a Z_2 is dynamical iff it is a SIMPLE-CURRENT MONODROMY
    GRADING -- a conserved fusion selection rule -- not an abstract automorphism. The
    neutrino's parity is exactly this: the SU(2)_1 primary j=1/2 is an order-2 simple current
    (1/2 x 1/2 = 0), grading the theory into even {0} / odd {1/2}; the neutrino sits in the
    odd class => Majorana, conserved (P17:cor:majorana).

  For a simple current J the monodromy charge of a primary phi is
       Q_J(phi) = h_J + h_phi - h_{J x phi}   (mod 1),
  valued in (1/N)Z for an order-N current; for N=2 it is {0, 1/2} = a Z_2 grading that FUSION
  CONSERVES. If the D-type dark state sits in the odd class and the vacuum in the even class,
  the parity is dynamical in exactly the neutrino's sense.

RESULT (what the script shows): YES. In the SO(2r)_1 assignment (r odd = our case: SO(10),
SO(14)), the VECTOR v is an order-2 simple current whose monodromy grading is the SPINOR
PARITY: {1, v} even, {s, c} odd. The D-type dark state IS the spinor (the SAME odd-rank
complex spinor h=r/8 that gives sterility + the Dirac property) -> it is Z_2-ODD. So the dark
parity is NOT a new assumption: it is the spinor parity, a conserved monodromy grading, the
odd-rank spinor's FOURTH job (after no-closure/sterile, complex/Dirac).
REMAINING [O]: this conserves parity WITHIN the D-type CFT; absolute stability against decay
to the SM sectors needs those sectors to be EVEN under the SAME grading -- the cross-sector
modular embedding (all primaries as 2I-vacuum reps). Physically a spinor cannot be built from
non-spinor (A/E) reps, so the decay is spinor-parity-forbidden; the formal embedding is the gap.
"""
import sympy as sp
from dtype_cft_tower_engine import so2r_1, su2k


def monodromy_charge(hJ, hphi, hJphi):
    """Q_J(phi) = h_J + h_phi - h_{J x phi}  (mod 1), as an exact rational in [0,1)."""
    return sp.Rational(hJ + hphi - hJphi) % 1


# ---- SO(2r)_1 center = simple-current group. For r ODD it is Z_4 = {1,s,v,c} with
#      s^2=v, s^4=1 (spinors order 4, vector v order 2). Fuse by Z_4 addition. ----
# center-element index: 1<->0, s<->1, v<->2, c<->3 ; fusion = (i+j) mod 4.
CENTER_LABEL = {0: '1', 1: 's', 2: 'v', 3: 'c'}
LABEL_CENTER = {v: k for k, v in CENTER_LABEL.items()}

def so2r_fuse(a, b):
    return CENTER_LABEL[(LABEL_CENTER[a] + LABEL_CENTER[b]) % 4]

def analyse_so2r(r, name):
    c, h = so2r_1(r)   # h = {'1':0,'v':1/2,'s':r/8,'c':r/8}
    print(f"\n{'='*74}\n  {name}:  SO({2*r})_1,  c={c},  primaries 1,v,s,c  h = "
          f"0, 1/2, {sp.Rational(r,8)}, {sp.Rational(r,8)}\n{'='*74}")
    # 1. identify the order-2 simple current: v (v x v = 1)
    vv = so2r_fuse('v', 'v')
    print(f"  Vector v: v x v = {vv}  -> order { 2 if vv=='1' else '?'} simple current "
          f"(spinors s,c have order 4: s x s = {so2r_fuse('s','s')})  [center = Z_4, r odd]")
    # 2. monodromy grading under v
    print(f"  Monodromy grading under the vector current v  (Q_v = h_v+h_phi-h_(v x phi) mod 1):")
    parity = {}
    for phi in ('1', 'v', 's', 'c'):
        Q = monodromy_charge(h['v'], h[phi], h[so2r_fuse('v', phi)])
        parity[phi] = 'EVEN' if Q == 0 else ('ODD' if Q == sp.Rational(1,2) else f'Q={Q}')
        print(f"    phi={phi:>2}:  v x {phi} = {so2r_fuse('v',phi):>2},  Q_v = {Q}  -> {parity[phi]}")
    even = [p for p in parity if parity[p]=='EVEN']; odd = [p for p in parity if parity[p]=='ODD']
    print(f"  => Z_2 grading: EVEN {{{','.join(even)}}} , ODD {{{','.join(odd)}}}.")
    print(f"     The D-type DARK STATE is the SPINOR (s/c) -> Z_2-ODD; vacuum 1 is EVEN.")
    print(f"     Fusion conserves this charge -> the parity is DYNAMICAL (a selection rule),")
    print(f"     the SAME spinor (h=r/8) that gives sterility + Dirac. FOURTH job.")
    return parity


def analyse_neutrino():
    c, h = su2k(1)   # SU(2)_1: {'l=0':0, 'l=1':1/4}
    print(f"\n{'='*74}\n  NEUTRINO reference:  SU(2)_1,  c={c},  primaries j=0 (h=0), j=1/2 (h=1/4)"
          f"\n{'='*74}")
    hJ, hphi = h['l=1'], h['l=1']        # J = phi = j=1/2
    Q_self = monodromy_charge(hJ, hphi, h['l=0'])   # 1/2 x 1/2 = 0
    print(f"  j=1/2 is an order-2 simple current (1/2 x 1/2 = 0).")
    print(f"  Monodromy Q(1/2 on 1/2) = h_(1/2)+h_(1/2)-h_0 = 1/4+1/4-0 = {Q_self}"
          f"  -> j=1/2 is Z_2-ODD, j=0 EVEN.")
    print(f"  => the neutrino sits in the ODD class = MAJORANA parity (P17:cor:majorana),")
    print(f"     a conserved monodromy grading. THIS is 'dynamical Z_2'. The D-type test above")
    print(f"     is the SAME construct -> same sense of 'dynamical'.")


if __name__ == '__main__':
    print("="*74)
    print("  IS THE D-TYPE DARK PARITY DYNAMICAL?  (simple-current monodromy test)")
    print("="*74)
    analyse_neutrino()
    analyse_so2r(5, "2D_3 dark sector")
    analyse_so2r(7, "2D_5 dark sector")
    print(f"\n{'='*74}")
    print("  VERDICT: YES -- the dark parity is a DYNAMICAL simple-current monodromy grading")
    print("  (the SPINOR PARITY, generated by the vector v), in the SAME sense as the")
    print("  neutrino's Majorana Z_2. The Dirac dark spinor is Z_2-ODD, the vacuum EVEN; fusion")
    print("  conserves it. => upgrades Paper XXI from 'absolute stability IFF Z_2 dynamical' to")
    print("  'the Z_2 IS dynamical (spinor parity)'. The dark parity is the odd-rank spinor's")
    print("  FOURTH job (sterility, Dirac, now stability), not a separate assumption.")
    print("  REMAINING [O]: cross-sector embedding -- that the A/E (SM) sectors are EVEN under")
    print("  the SAME grading (a spinor cannot be built from non-spinor reps -> physically")
    print("  forbidden; formal 2I-vacuum modular embedding is the gap).")
    print("="*74)
