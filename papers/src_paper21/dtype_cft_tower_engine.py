#!/usr/bin/env python3.11
"""
dtype_cft_tower_engine.py -- D-type (binary-dihedral) sector via T-matrix closure.
Companion to e6_cft_tower_engine.py. Tests whether a dicyclic dark-Hopfion sector
(2D_3, 2D_5 -- the recombination-CDM candidates, notes/fabric_cmb_hypothesis_and_plan.md)
admits a mass tower.

Q_group = min integer Q>=1 s.t. Q*T_j is an ODD multiple of 1/4 for EVERY primary j
(P15:def:closes_tmatrix), T_j = h_j - c/24. The mass tower
m = Lcond * phi^{-2n} * exp(-2 pi Q_group T_g)  (Paper VI / XVII) NEEDS Q_group.

Two natural D-type WZW assignments tested:
  (a) (D_{n+2})_1 = SO(2(n+2))_1  [quark convention: (McKay Lie algebra)_1]:
        2D_3 -> SO(10)_1,  2D_5 -> SO(14)_1.
  (b) SU(2)_k with the D-type modular invariant [lepton convention: SU(2)_k]:
        2D_3 -> D_5 inv -> SU(2)_6,  2D_5 -> D_7 inv -> SU(2)_10   (k=2n).

RESULT (2026-10-05): BOTH fail closure. (a) spinor obstruction (h=r/8, r odd forces
Q*r odd vs vacuum Q*r even); (b) even-k obstruction (all even k fail; only ODD k --
the SM's SU(2)_1, SU(2)_3 -- close). => the D-type carries NO quantized charge (the
precise/3rd confirmation of darkness) and has NO Q_group => the mass tower exponent is
undefined => coldness is NOT determinable here; mass is set only by the undetermined
thermodynamic tower level n (like, but weaker than, the neutrino's unfixed n_nu).
Charge-less <=> mass-unpinned: the dark sector lacks both quantization conditions.
"""
import sympy as sp

def su2k(k):
    c = sp.Rational(3*k, k+2)
    h = {f"l={l}": sp.Rational(l*(l+2), 4*(k+2)) for l in range(k+1)}
    return c, h

def e6_1():
    return sp.Integer(6), {'1': sp.Integer(0), '27': sp.Rational(2,3), '27b': sp.Rational(2,3)}

def so2r_1(r):
    """SO(2r)_1: c=r; primaries 1,v,s,c with h = 0, 1/2, r/8, r/8."""
    return sp.Integer(r), {'1': sp.Integer(0), 'v': sp.Rational(1,2),
                           's': sp.Rational(r,8), 'c': sp.Rational(r,8)}

def closes(Q, Tlist):
    # Q*T is an odd multiple of 1/4  <=>  4*Q*T is an odd integer (exact rationals)
    for T in Tlist:
        x = sp.Rational(4*Q*T)
        if x.q != 1 or x.p % 2 == 0:
            return False
    return True

def Qgroup(c, h, Qmax=5000):
    T = [sp.Rational(hj - c/24) for hj in h.values()]
    for Q in range(1, Qmax+1):
        if closes(Q, T):
            return Q, T
    return None, T

def line(name, c, h, expect=None):
    Q, T = Qgroup(c, h)
    tag = "" if expect is None else f"  (expect {expect}: {'OK' if Q==expect else 'FAIL'})"
    print(f"  {name:28} c={str(c):4}  Q_group={Q}{tag}")
    return Q

if __name__ == '__main__':
    print("="*74)
    print("  VALIDATION: known sectors must reproduce Q_group {3,6,10}")
    print("="*74)
    line("(E_6)_1  [quark 2T]",       *e6_1(),   expect=3)
    line("SU(2)_1  [neutrino A/Z2]",  *su2k(1),  expect=6)
    line("SU(2)_3  [lepton 2I]",      *su2k(3),  expect=10)

    print("\n"+"="*74)
    print("  D-TYPE candidates, convention (a): (D_{n+2})_1 = SO(2(n+2))_1")
    print("="*74)
    line("SO(10)_1 [2D_3, D_5]", *so2r_1(5))
    line("SO(14)_1 [2D_5, D_7]", *so2r_1(7))

    print("\n"+"="*74)
    print("  D-TYPE candidates, convention (b): SU(2)_k with D-invariant (k=2n)")
    print("  + even/odd-k pattern (SM sectors are ODD k)")
    print("="*74)
    for k in (1, 2, 3, 4, 6, 10):
        line(f"SU(2)_{k}", *su2k(k))

    print("\n  CONCLUSION: D-type fails T-matrix closure by BOTH conventions ->")
    print("  no Q_group -> mass tower exponent undefined -> coldness undeterminable.")
    print("  Same spinorial/even structure that makes it DARK removes it from the tower.")
