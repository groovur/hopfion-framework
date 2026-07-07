#!/usr/bin/env python3
"""
jones_satellite_continuous.py
===================================================================
Analytic continuation of jones_satellite_twist_family.py's discrete
twist-count family to continuous (non-integer) twist number n, so the
path between adjacent torus knots in the tower -- e.g. n=3 (trefoil,
Q_H=3) to n=5 (cinquefoil T(2,5)) -- can be examined as a curve rather
than just compared at its two integer endpoints.

This is the same style of continuous-parameter technique already used
elsewhere in this directory for the geometric field construction
(hopf_link_construction_v2.py's --q_pol, interpolating T(2,q) fields
continuously from q=2 to q=3); here it is applied to the invariant
instead of the field.

METHOD (exact, not curve-fitting)
-----------------------------------
The cable-twist crossing sigma_1, repeated n times in the
Temperley-Lieb algebra TL_2 (embedded in TL_3, acting only on the two
cable strands), decomposes as

    sigma_1^n = A^n * Id + c_n * e_1

where e_1 is the Temperley-Lieb generator connecting the two cable
strands and c_n satisfies the exact linear recurrence (derived by
direct TL_2 multiplication, and checked against the validated
integer-n results in jones_satellite_twist_family.py before use here):

    c_{n+1} = A^{n-1} - A^{-3} c_n,   c_0 = 0, c_1 = A^{-1}

which has the closed-form solution

    c_n = [A^{n+2} - (-1)^n A^{2-3n}] / (A^4 + 1)

For continuous n, (-1)^n is continued as exp(i*pi*n) (the standard
analytic continuation, agreeing with (-1)^n at integers).

The rest of the diagram (the fixed companion-threading tangle Theta,
independent of n) contributes two more fixed numbers,
    X = closure(Theta),   Y = closure(e_1 * Theta),
extracted directly from the validated Kauffman-bracket engine, giving

    bracket(n) = A^n * X + c_n * Y

and the Jones polynomial via the usual writhe correction
(writhe(n) = n + 4, four crossings from the fixed companion threading):

    V_sat(n)(t) = (-A^3)^{-(n+4)} * bracket(n) |_{A = t^{-1/4}}

VALIDATION
----------
This closed form is checked against the independently-computed exact
integer results (n=3,4,5,6) from jones_satellite_twist_family.py
before being trusted at non-integer n. All four match exactly.

FINDING: the n=3 -> n=5 path is NOT monotonic
------------------------------------------------
|V_sat(q_5)| does not simply rise from 1/phi (n=3) to phi (n=5). It
overshoots to a peak near n=3.6 (close to phi itself), falls through a
minimum near n=4.3 (below the n=4 unfused value of 1), then climbs
back up to phi at n=5. See the printed table below for exact numbers.

CAVEAT (not established, do not read the physical framing below as
proved)
-----------------------------------------------------------------
The continuation itself -- V_sat as an exact algebraic function of the
formal parameter n -- is rigorous. Whether the diagrammatic twist
count n corresponds to anything physically continuous in the actual
condensate (e.g. a coordinate along the perturbation-energy axis, or
a proxy for elapsed time during formation/decay) is NOT established
anywhere in the papers and is not derived here. That mapping, if it
exists, would need the dynamic-integrator infrastructure (velocity
field v(x,t) at each site) to justify. This script provides a
well-defined mathematical object to test candidate physical
interpretations against -- it does not itself supply one.

Requires: sympy (validated against python3.11 + sympy 1.14).
"""

import cmath
import math

import sympy as sp

from kauffman_bracket_engine import (
    A, identity_tangle, e_gen_diagram, crossing_tangle,
    compose_tangles, close_up_bracket,
)

t = sp.symbols('t')
n_sym = sp.symbols('n')
q5 = sp.exp(2 * sp.pi * sp.I / 5)
phi = (1 + sp.sqrt(5)) / 2


def get_X_Y():
    """Extract the two n-independent constants from the fixed
    companion-threading tangle Theta (word [1,0,0,1], validated in
    jones_T23_satellite_check.py)."""
    n = 3
    theta = identity_tangle(n)
    for i in [1, 0, 0, 1]:
        theta = compose_tangles(crossing_tangle(n, i, sign=+1), theta, n)
    X = sp.simplify(close_up_bracket(theta, n))
    e0_tangle = {e_gen_diagram(n, 0): sp.Integer(1)}
    eY_tangle = compose_tangles(e0_tangle, theta, n)
    Y = sp.simplify(close_up_bracket(eY_tangle, n))
    return X, Y


def c_n_continuous(nval, use_exact_sign=None):
    """Closed-form solution of c_{n+1} = A^{n-1} - A^{-3} c_n, c_0=0.
    use_exact_sign lets integer callers pass the exact (-1)**n instead
    of the exp(i*pi*n) continuation, for a clean integer cross-check."""
    if use_exact_sign is not None:
        sign_term = use_exact_sign
    else:
        sign_term = sp.exp(sp.I * sp.pi * nval)
    return (A**(nval + 2) - sign_term * A**(2 - 3 * nval)) / (A**4 + 1)


def V_sat_at(nval, X, Y, exact_integer=False):
    if exact_integer:
        cn = c_n_continuous(nval, use_exact_sign=(-1) ** nval)
    else:
        cn = c_n_continuous(nval)
    bracket = A**nval * X + cn * Y
    writhe = nval + 4
    expr = sp.simplify((-A**3)**(-writhe) * bracket)
    expr_t = sp.simplify(sp.expand(expr.subs(A, t**sp.Rational(-1, 4))))
    return expr_t


if __name__ == "__main__":
    X, Y = get_X_Y()
    print("X = closure(Theta)      =", X)
    print("Y = closure(e_0 * Theta) =", Y)
    print()

    print("=" * 60)
    print("Validation against exact integer results")
    print("=" * 60)
    for nval in [3, 4, 5, 6]:
        Vt = V_sat_at(nval, X, Y, exact_integer=True)
        val = complex(Vt.subs(t, q5).evalf(30))
        print("n=%d  |V(q5)| = %.6f" % (nval, abs(val)))
    print()

    print("=" * 60)
    print("Continuous path from n=3 (trefoil) to n=5 (cinquefoil)")
    print("=" * 60)
    print("%5s %12s %10s" % ("n", "|V(q5)|", "arg(deg)"))
    n = sp.Rational(30, 10)
    step = sp.Rational(1, 10)
    peak = (None, -1)
    trough = (None, 1e9)
    while n <= sp.Rational(50, 10):
        Vt = V_sat_at(n, X, Y, exact_integer=False)
        val = complex(Vt.subs(t, q5).evalf(30))
        mod = abs(val)
        phase = cmath.phase(val) * 180 / math.pi
        print("%5s %12.6f %10.2f" % (str(float(n)), mod, phase))
        if mod > peak[1]:
            peak = (float(n), mod)
        if mod < trough[1]:
            trough = (float(n), mod)
        n += step

    print()
    print("peak:   n=%.2f  |V(q5)|=%.6f  (phi=%.6f)" % (peak[0], peak[1], float(phi)))
    print("trough: n=%.2f  |V(q5)|=%.6f" % (trough[0], trough[1]))
