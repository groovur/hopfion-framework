#!/usr/bin/env python3
"""
jones_T23_satellite_check.py
===================================================================
Checks Paper XVIII's Theorem 2.3 (thm:jones_eval / P18:thm:jones_eval):
    V_{T(2,3)}(q_5) = -1,   q_5 = e^{2 pi i / 5}
against direct computation, and computes the mathematically more
appropriate quantity -- the Jones polynomial of the trefoil as it is
actually constructed physically in Paper XVIII (Prop. 5.2 / satellite
construction): a (2,3)-cable of one component of the Q_H=2 Hopf link,
with the companion component still present and linked.

FINDINGS (see printed output below for the numbers)
----------------------------------------------------
1. Direct evaluation of the free-standing trefoil's standard Jones
   polynomial at q_5 does NOT give -1. It gives a genuinely irrational,
   non-real complex number of modulus sqrt(4 - phi) (phi = golden
   ratio), i.e. Theorem 2.3 does not check out as stated.

2. V_K(i) = (-1)^Arf(K) (Lickorish-Millett 1986) is where the "-1"
   claim actually comes from: the trefoil's Arf invariant is 1, giving
   V(i) = -1 exactly -- but at the 4th root of unity (SU(2)_2 /
   Ising-type point), not the 5th root (q_5) that Papers XVII-XIX use
   throughout for the k=3 WZW condensate. This is a coincidence of a
   different, unrelated piece of classical knot theory (the Arf
   invariant), not a WZW/quantum-dimension fact.

3. The trefoil is not actually a free-standing knot in this framework
   -- Prop. 5.2 builds it as the (2,3)-cable of one component of the
   Q_H=2 Hopf link, with the other Hopf-link component (the "companion")
   still present. The properly-computed satellite invariant (companion
   included, via the Kauffman-bracket/Temperley-Lieb cabling formula)
   is ALSO irrational at q_5, but a THIRD distinct value:
       V(q_5) = q_5^3 / phi,   |V(q_5)| = 1/phi
   This was cross-checked against an independent attempt at the
   cabling substitution which was caught to be WRONG (see the
   "known-bad" section below) -- kept in this file deliberately, to
   show the failure mode and why the accepted result should be
   trusted instead of just asserted.

CAVEAT (not resolved by this script -- a physics/reaction-bookkeeping
question, not a math one): whether the companion Hopf-link component
genuinely survives as a spectator after the Q_H=2+Q_H=1 -> Q_H=3
(Channel 2) reaction is not clearly pinned down in Papers
XVIII/XIX -- Channel 2's reaction equation shows no leftover particle,
unlike Channel 1. If the companion is actually consumed/absorbed
rather than surviving, the free-trefoil value (#1 above) may be the
physically relevant one instead of the satellite value (#3). Either
way, phi does not disappear at Q_H=3, contradicting P18:cor:no_phi's
claim that "the golden ratio phi is absent from every WZW amplitude of
the Q_H=3 sector at the k=3 level."

Requires: sympy (validated against python3.11 + sympy 1.14).
"""

import sympy as sp
from kauffman_bracket_engine import (
    A, identity_tangle, e_tangle, crossing_tangle,
    compose_tangles, close_up_bracket,
)

t = sp.symbols('t')
q5 = sp.exp(2 * sp.pi * sp.I / 5)
phi = (1 + sp.sqrt(5)) / 2


def jones_from_bracket(bracket, writhe):
    """V(t) = (-A^3)^{-w} * <bracket>, then substitute A = t^{-1/4}."""
    expr = sp.expand(sp.simplify((-A**3)**(-writhe) * bracket))
    return sp.simplify(sp.expand(expr.subs(A, t**sp.Rational(-1, 4))))


def permutation_after_word(n, word):
    perm = list(range(n))
    for i in word:
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
    return perm


def count_components(n, perm):
    seen = [False] * n
    comps = 0
    for i in range(n):
        if seen[i]:
            continue
        comps += 1
        j = i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
    return comps


print("=" * 70)
print("STEP 0 -- engine validation against known Jones polynomials")
print("=" * 70)

# Trefoil = closure of sigma_1^3 on 2 strands
n = 2
tangle = identity_tangle(n)
for _ in range(3):
    tangle = compose_tangles(crossing_tangle(n, 0, sign=+1), tangle, n)
bracket_tref = close_up_bracket(tangle, n)
V_tref_free = sp.simplify(jones_from_bracket(bracket_tref, writhe=3))
print("Trefoil (free, closure of sigma_1^3): V(t) =", V_tref_free)
print("  expected (standard tables):        V(t) = -t**4 + t**3 + t")
assert sp.simplify(V_tref_free - (-t**4 + t**3 + t)) == 0
print("  MATCH -- engine validated on trefoil.\n")

# Hopf link = closure of sigma_1^2 on 2 strands
tangle2 = identity_tangle(n)
for _ in range(2):
    tangle2 = compose_tangles(crossing_tangle(n, 0, sign=+1), tangle2, n)
bracket_hopf = close_up_bracket(tangle2, n)
V_hopf = sp.simplify(jones_from_bracket(bracket_hopf, writhe=2))
print("Hopf link (closure of sigma_1^2): V(t) =", V_hopf)
print("  expected (standard tables):     V(t) = -t**(1/2) - t**(5/2)")
assert sp.simplify(V_hopf - (-t**sp.Rational(1, 2) - t**sp.Rational(5, 2))) == 0
print("  MATCH -- engine validated on Hopf link.\n")

print("=" * 70)
print("STEP 1 -- Theorem 2.3's claim, checked directly")
print("=" * 70)
V_at_q5 = complex(V_tref_free.subs(t, q5).evalf(30))
print("Free trefoil V(q_5) =", V_at_q5, " modulus =", abs(V_at_q5))
print("Claimed by P18:thm:jones_eval: V(q_5) = -1")
print(f"  --> does NOT match (modulus {abs(V_at_q5):.6f} != 1)\n")

# closed form check: |V(q_5)|^2 = 2 + 1/phi^2 = 4 - phi
lhs = sp.nsimplify(abs(V_at_q5)**2, [sp.sqrt(5)])
print("Exact closed form of |V(q_5)|^2 for the free trefoil: 4 - phi =",
      sp.N(4 - phi), " (matches numeric:", abs(V_at_q5)**2, ")\n")

print("Where -1 DOES occur (Lickorish-Millett / Arf invariant check):")
q4 = sp.I  # 4th root of unity
V_at_q4 = sp.simplify(V_tref_free.subs(t, q4))
print("  Free trefoil V(i) =", V_at_q4, " (exact, no rounding)")
print("  This is the classical identity V_K(i) = (-1)^Arf(K); trefoil has Arf=1.")
print("  It is UNRELATED to the k=3 / q_5 WZW structure used throughout the papers.\n")

print("=" * 70)
print("STEP 2 -- the properly-computed satellite (companion included)")
print("=" * 70)
print("Physical picture (Prop. 5.2): T(2,3) = (2,3)-cable of ONE component of")
print("the Q_H=2 Hopf link, with the OTHER component (companion) still linked.\n")

n = 3
# cable twist: 3 crossings between the two cable strands (positions 0,1)
cable_word = [0, 0, 0]
# companion threading: verified by hand, strand-by-strand, to give exactly
# 2 same-sign crossings with EACH cable strand (linking number 1 with each,
# matching the validated Hopf-link building block above), and to return
# every strand to its starting position (a clean "through and back").
thread_word = [1, 0, 0, 1]
full_word = cable_word + thread_word

perm = permutation_after_word(n, full_word)
ncomp = count_components(n, perm)
print("Full word:", full_word)
print("Induced permutation:", perm, " -> number of link components:", ncomp)
assert ncomp == 2, "expected 2 components: fused cable pair + separate companion"
print("  OK: 2 components (cable pair fuses via the odd 3-twist; companion separate)\n")

tangle = identity_tangle(n)
for i in full_word:
    tangle = compose_tangles(crossing_tangle(n, i, sign=+1), tangle, n)
bracket = close_up_bracket(tangle, n)
writhe = len(full_word)
V_sat = sp.simplify(jones_from_bracket(bracket, writhe))
print("Satellite (cable + companion) Jones polynomial: V(t) =", V_sat)

val = complex(V_sat.subs(t, q5).evalf(30))
print("V(q_5) numeric:", val, " modulus:", abs(val))
closed_form = q5**3 / phi
print("Closed form check q_5^3/phi:", complex(sp.N(closed_form, 30)))
assert abs(val - complex(sp.N(closed_form, 30))) < 1e-15
print("  MATCH: V(q_5) = q_5^3 / phi  exactly, |V(q_5)| = 1/phi =", float(1 / phi), "\n")

print("=" * 70)
print("KNOWN-BAD ATTEMPT (kept deliberately -- caught, not hidden)")
print("=" * 70)
print("A naive 'systematic cabling substitution' (blindly mapping each of the")
print("Hopf link's 2 crossings to a fixed 2-crossing pattern without re-tracking")
print("where the doubled strand had moved) gives a DIFFERENT, wrong answer:")
bad_word = [0, 0, 0] + [1, 0, 1, 0]
bad_perm = permutation_after_word(n, bad_word)
bad_ncomp = count_components(n, bad_perm)
print("  bad word:", bad_word, " permutation:", bad_perm,
      " components:", bad_ncomp)
print("  The permutation shows the companion swapping places with a cable")
print("  strand instead of returning cleanly -- i.e. this word secretly encodes")
print("  a DIFFERENT link, not the intended one. Diagnosed via the same")
print("  permutation/component-count check used to validate the accepted word.")
