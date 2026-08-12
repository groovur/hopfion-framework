"""Exact checks for the anyonic identification of the golden-angle
attractor (main_paper3.tex, "The attractor as an anyonic observable"
and "The angular quantum from the orbifold spectrum").

V2 : Tr[rho U_mono] = p0 e^{-i 108 deg} + p1 e^{+i 36 deg} = 1/phi^2,
     exactly real, and equals the interferometry form
     S_{1/2,1/2} S_00 / S_{0,1/2}^2.
V3 : S_{1/2,1/2}/S_00 = 1.
Q  : phase group of the physical modular T-exponents (h - c/24):
     untwisted SU(2)_3 -> Z_40; adding the 2I-orbifold twisted sector
     (Coxeter exponents, h_m = m/90; Paper V) -> Z_360 = Z_{k|2I|},
     generator the m = 7 primary (exponent exactly 1/360); without the
     -c/24 shift the twisted exponents generate only Z_90.

Run: python3.11 trace_check.py   (all assertions must pass)
"""
from fractions import Fraction as F
from math import lcm

import sympy as sp

phi = (1 + sp.sqrt(5)) / 2

# ---- V2: monodromy expectation -------------------------------------
p0, p1 = 1 / phi**2, 1 / phi
M0 = sp.exp(-sp.I * sp.pi * sp.Rational(3, 5))   # -108 deg
M1 = sp.exp(sp.I * sp.pi / 5)                    # +36 deg
tr = p0 * M0 + p1 * M1
diff = sp.expand_complex(sp.expand(tr - 1 / phi**2, complex=True))
assert sp.simplify(sp.re(diff)) == 0 and sp.simplify(sp.im(diff)) == 0
print("V2: Tr[rho U_mono] = 1/phi^2 =", sp.nsimplify(sp.simplify(sp.re(tr)), [sp.sqrt(5)]), "(imaginary part 0 exactly)")

# interferometry form
S = lambda a, b: sp.sqrt(sp.Rational(2, 5)) * sp.sin((2*a+1)*(2*b+1)*sp.pi/5)
half = sp.Rational(1, 2)
assert sp.simplify(S(half, half) * S(0, 0) / S(0, half)**2 - 1 / phi**2) == 0
print("V2: equals S_{1/2 1/2} S_00 / S_{0 1/2}^2")

# ---- V3 -------------------------------------------------------------
assert sp.simplify(S(half, half) / S(0, 0) - 1) == 0
print("V3: S_{1/2 1/2}/S_00 = 1")

# ---- Phase-group orders --------------------------------------------
c24 = F(9, 5) / 24                                # c/24 = 27/360
untwisted = [F(0), F(3, 20), F(2, 5), F(3, 4)]
coxeter = [1, 7, 11, 13, 17, 19, 23, 29]
un = [h - c24 for h in untwisted]
tw = [F(m, 90) - c24 for m in coxeter]

order = lambda xs: lcm(*[x.denominator for x in xs])
assert order(un) == 40
assert order(un + tw) == 360
assert tw[1] == F(1, 360)                          # m = 7 generator
assert order([F(m, 90) for m in coxeter]) == 90    # without -c/24
assert 360 == 3 * 120                              # k |2I|
print("Q: untwisted Z_40; with 2I twisted sector Z_360 = Z_{k|2I|};")
print("   m=7 exponent = 1/360; without -c/24 shift only Z_90")
print("all checks passed")
