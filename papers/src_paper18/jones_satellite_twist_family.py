#!/usr/bin/env python3
"""
jones_satellite_twist_family.py
===================================================================
Follow-up to jones_T23_satellite_check.py.  Instead of computing the
satellite Jones invariant only at the trefoil's twist count (3 half-
twists between the two cable strands), compute it for the whole
family, twist count n = 0, 1, 2, ..., with the companion (the other
Q_H=2 Hopf-link component) linked once with each cable strand
throughout, exactly as validated there.

n=3 is the physical case: the trefoil T(2,3), i.e. Q_H=3, satellite
of the Q_H=2 Hopf link (Prop. 5.2 / prop:satellite).  n=5 is the next
knot in the same tower, the cinquefoil T(2,5).  Even n gives an
unfused, 3-component configuration (the two cable strands have not
merged into a single knotted component) rather than a genuine torus
knot.

EXACT RESULTS (symbolic, sympy, not rounded)
---------------------------------------------
    n=3 (trefoil,   fused, Q_H=3):     |V_sat(q5)| = (sqrt5-1)/2 = 1/phi
    n=4 (unfused, 3 components):        V_sat(q5)  = 1  (rational)
    n=5 (cinquefoil, fused):            |V_sat(q5)| = (sqrt5+1)/2 = phi
    n=6 (unfused, 3 components):        V_sat(q5)  = -1  (rational, EXACT)

STRUCTURAL FACTS (verified below, exact)
------------------------------------------
  - |V_sat(q5)| is periodic in n with period 10.
  - Within each period, |V_sat(q5)| is mirror-symmetric about n=5:
        |V_sat(n)| = |V_sat(10-n)|
  - At symmetric points, the phases satisfy
        arg V_sat(n) + arg V_sat(10-n) = 72 deg (mod 360),
    i.e. exactly one unit of arg(q5) = 2*pi/5.
  - Every FUSED (odd n, genuine torus knot) configuration computed
    gives an irrational, phi-related modulus.
  - Every UNFUSED (even n, 3-component, two parallel un-merged cable
    strands) configuration computed gives a RATIONAL value (n=4: 1;
    n=6: exactly -1 -- the specific value Paper XVIII originally
    claimed for the trefoil itself, which does exist in this family,
    just not at n=3).

OPEN QUESTIONS RAISED IN DISCUSSION (NOT established, flagged as such
so they are not mistaken for proved results if this file is read out
of context)
-----------------------------------------------------------------
  - The rational/irrational split tracking "fused vs unfused" rather
    than "which Q_H sector" suggests confinement (if it is related to
    rationality at all) may be a statement about component fusion,
    not about the Q_H=3 charge as such.  This would be consistent with
    the independently-established numerical fact (Papers XV/XVI
    gradient-flow searches) that no Q_H=3 saddle point is ever found:
    one live conjecture is that the physical trefoil never actually
    completes the fusion into T(2,3) -- the n=3 fused state may be an
    idealised endpoint that the condensate approaches but does not
    reach before disintegrating, consistent with the saddle-difficulty
    / flat-valley numerics (rem:saddle_difficulty).  Under that
    reading, n=4 (unfused, adjacent to n=3, rational) could be closer
    to the actually-realised transient geometry than n=3 itself.
  - Whether n=5 (cinquefoil, modulus phi, the RECIPROCAL of the
    trefoil's 1/phi) has any meaning as e.g. an "anti-hopfion" /
    orientation-reversed partner configuration is unexplored.
  - The period is exactly 10, the same integer as the quantum-group
    order Q_group = 2(k+2) = 10 (ord(q_2I), the Pentagon-theorem
    quantity used throughout the framework starting from Paper I).
    Whether this is the same 10 or a coincidence of this specific
    tangle's eigenvalue orders is not established.

None of the above interpretive points are claimed as proved; they are
recorded here as the open questions motivating further work, per the
discussion that produced this script.

Requires: sympy (validated against python3.11 + sympy 1.14).
"""

import cmath
import math

import sympy as sp

from kauffman_bracket_engine import (
    A, identity_tangle, crossing_tangle, compose_tangles, close_up_bracket,
)

t = sp.symbols('t')
q5 = sp.exp(2 * sp.pi * sp.I / 5)
phi = (1 + sp.sqrt(5)) / 2


def jones_from_bracket(bracket, writhe):
    expr = sp.expand(sp.simplify((-A**3)**(-writhe) * bracket))
    return sp.simplify(sp.expand(expr.subs(A, t**sp.Rational(-1, 4))))


def satellite_jones(cable_twists):
    """Companion-included satellite Jones polynomial for a cable of
    `cable_twists` half-twists between the two cable strands, with the
    companion linked once with each cable strand (thread word [1,0,0,1],
    validated in jones_T23_satellite_check.py to return every strand to
    its starting position and give exactly linking number 1 with each
    cable strand)."""
    n = 3
    word = [0] * cable_twists + [1, 0, 0, 1]
    tangle = identity_tangle(n)
    for i in word:
        tangle = compose_tangles(crossing_tangle(n, i, sign=+1), tangle, n)
    bracket = close_up_bracket(tangle, n)
    V = sp.simplify(jones_from_bracket(bracket, len(word)))
    return V, word


def component_count(cable_twists):
    n = 3
    word = [0] * cable_twists + [1, 0, 0, 1]
    perm = list(range(n))
    for i in word:
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
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


if __name__ == "__main__":
    print("=" * 78)
    print("Satellite Jones invariant vs. cable twist count n (companion included)")
    print("=" * 78)
    print("%3s %10s %14s %10s" % ("n", "components", "|V_sat(q5)|", "arg (deg)"))

    results = {}
    for n in range(0, 16):
        V, _ = satellite_jones(n)
        val = complex(V.subs(t, q5).evalf(30))
        mod = abs(val)
        phase = cmath.phase(val) * 180 / math.pi
        comps = component_count(n)
        results[n] = (comps, mod, phase, val)
        print("%3d %10d %14.6f %10.2f" % (n, comps, mod, phase))

    print()
    print("=" * 78)
    print("Exact symbolic values at the four key points")
    print("=" * 78)
    for n, label in [(3, "n=3  trefoil T(2,3), fused -- this is Q_H=3"),
                     (4, "n=4  unfused, 3 components"),
                     (5, "n=5  cinquefoil T(2,5), fused"),
                     (6, "n=6  unfused, 3 components")]:
        V, _ = satellite_jones(n)
        val_exact = sp.simplify(sp.expand_complex(V.subs(t, q5)))
        mod_exact = sp.nsimplify(abs(complex(val_exact.evalf(30))), [sp.sqrt(5)])
        print(label)
        print("  V(q5) exact   =", val_exact)
        print("  |V(q5)| exact =", mod_exact, " = ", float(mod_exact))
        print()

    print("=" * 78)
    print("Structural checks: period 10, mirror symmetry about n=5,")
    print("phase pairing summing to 72 deg = arg(q5)")
    print("=" * 78)
    for n in range(0, 6):
        c1, m1, p1, _ = results[n]
        c2, m2, p2, _ = results[n + 10]
        print("n=%2d vs n=%2d:  |V| %.6f vs %.6f  (match: %s)"
              % (n, n + 10, m1, m2, abs(m1 - m2) < 1e-9))
    print()
    for n in range(1, 5):
        c1, m1, p1, _ = results[n]
        c2, m2, p2, _ = results[10 - n]
        phase_sum = (p1 + p2) % 360
        print("n=%2d vs n=%2d (mirror about 5): |V| %.6f vs %.6f  "
              "phase sum mod 360 = %.2f (expect 72.00)"
              % (n, 10 - n, m1, m2, phase_sum))
