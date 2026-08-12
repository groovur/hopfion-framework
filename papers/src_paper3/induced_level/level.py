"""
Exact symbolic exploration of twisted-sector quantum dimensions and the
induced abelian Chern-Simons level for the condensate CFT SU(2)_3 with a
2I-orbifold twisted sector.

All results are exact: rationals, elements of Q(sqrt5), or (where the
twisted-sector candidate forces it) exact algebraic numbers of higher
degree, certified either by an explicit denested radical or by a minimal
polynomial. Floating point (mpmath, high precision) is used only to
recognize / cross-check closed forms -- never as the reported value.

Run: python3.11 level.py
"""

import sympy as sp
import mpmath as mp

sp.init_printing(use_unicode=False)

# ----------------------------------------------------------------------
# Framework data (exact)
# ----------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
c = sp.Rational(9, 5)
c24 = c / 24  # = 3/40

# Untwisted SU(2)_3 primaries
h_un = [sp.Integer(0), sp.Rational(3, 20), sp.Rational(2, 5), sp.Rational(3, 4)]
d_un = [sp.Integer(1), phi, phi, sp.Integer(1)]
D_untw2 = sp.radsimp(sp.simplify(sum(d ** 2 for d in d_un)))  # = 5 + sqrt(5)

# Twisted sector: E8 Coxeter exponents, h_m = m/90
m_list = [1, 7, 11, 13, 17, 19, 23, 29]
h_tw = [sp.Rational(m, 90) for m in m_list]

print("=" * 78)
print("Untwisted-sector data")
print("=" * 78)
print("h_un  =", h_un)
print("d_un  =", d_un, " (phi =", phi, ")")
print("D_untw^2 = sum d_a^2 =", D_untw2, " = 2 + 2*phi^2, numeric:", sp.N(D_untw2, 30))
print()

# ----------------------------------------------------------------------
# Candidate A1: d_m = 1  (abelian twisted sector)
# ----------------------------------------------------------------------
d_A1 = [sp.Integer(1)] * 8

# ----------------------------------------------------------------------
# Candidate A2: finite-E8 Dynkin marks {2,3,4,5,6,4,2,3}, matched
# positionally to the exponent list m_list = [1,7,11,13,17,19,23,29] in
# the order both are given in the task statement (increasing-exponent
# order <-> stated mark order). Explicit pairing:
#   m=1 -> mark 2,  m=7  -> mark 3,  m=11 -> mark 4,  m=13 -> mark 5,
#   m=17 -> mark 6, m=19 -> mark 4,  m=23 -> mark 2,  m=29 -> mark 3.
# (Sanity check: these marks are exactly the dimensions of the 8
# nontrivial irreps of the binary icosahedral group 2I -- the finite
# nodes of the affine E8 Dynkin diagram under McKay correspondence --
# and sum to h-1 = 29 with h = 30 the E8 Coxeter number.)
# ----------------------------------------------------------------------
marks = [2, 3, 4, 5, 6, 4, 2, 3]
A2_pairing = list(zip(m_list, marks))
d_A2 = [sp.Integer(x) for x in marks]

# ----------------------------------------------------------------------
# Candidate A3: d_m = sin(pi m/30) / sin(pi/30)  (quantum dims at the
# E8 Coxeter root of unity q = exp(i pi/30), h = 30).
# These are honest algebraic numbers but NOT elements of Q(sqrt5): each
# nonzero, non-unit value has degree 4 over Q. Certified below by exact
# minimal polynomials (computed symbolically, not fit numerically).
# ----------------------------------------------------------------------
_sin1 = sp.sin(sp.pi / 30).rewrite(sp.sqrt)


def d3_exact(m):
    s = sp.sin(sp.pi * m / 30).rewrite(sp.sqrt)
    return sp.radsimp(sp.simplify(s / _sin1))


d_A3 = [sp.nsimplify(sp.N(d3_exact(m), 40)) if False else d3_exact(m) for m in m_list]

print("=" * 78)
print("Candidate A3: exact radical values and minimal polynomials")
print("=" * 78)
_seen = {}
for m in m_list:
    val = d_A3[m_list.index(m)]
    key = sp.N(val, 15)
    if key not in _seen:
        mp_ = sp.minimal_polynomial(val, sp.Symbol("x"))
        _seen[key] = mp_
        print(f"  m={m:2d}: d_m minimal polynomial = {mp_}   (degree {sp.degree(mp_)})")
    print(f"        d_{m} = {sp.N(val, 25)}")
print()

CANDIDATES = {
    "A1": d_A1,
    "A2": d_A2,
    "A3": d_A3,
}

# ----------------------------------------------------------------------
# Question A: total D^2 and modularity gate (exact where tractable,
# else certified via 60-digit numerics -- the gaps below are orders of
# magnitude too large to be finite-precision artifacts).
# ----------------------------------------------------------------------
mp.mp.dps = 60


def mp_d(sym_val):
    return mp.mpf(str(sp.N(sym_val, 60)))


mp_phi = (1 + mp.sqrt(5)) / 2
mp_h_un = [mp.mpf(0), mp.mpf(3) / 20, mp.mpf(2) / 5, mp.mpf(3) / 4]
mp_d_un = [mp.mpf(1), mp_phi, mp_phi, mp.mpf(1)]
mp_h_tw = [mp.mpf(m) / 90 for m in m_list]

results_A = {}

print("=" * 78)
print("Question A: total D^2 and modularity gate, per candidate")
print("=" * 78)

for label, d_tw in CANDIDATES.items():
    # exact D^2 (works cleanly for A1, A2; A3 stays exact via sympy too,
    # just algebraically larger)
    D2_tw = sp.radsimp(sp.simplify(sum(d ** 2 for d in d_tw)))
    D2_total = sp.radsimp(sp.simplify(D_untw2 + D2_tw))

    # high-precision numeric cross-check of D^2, and full Gauss sums
    d_tw_mp = [mp_d(d) for d in d_tw]
    h_all_mp = mp_h_un + mp_h_tw
    d_all_mp = mp_d_un + d_tw_mp

    D2_mp = sum(d ** 2 for d in d_all_mp)
    D_mp = mp.sqrt(D2_mp)

    Thp = mp.mpc(0)
    for h, d in zip(h_all_mp, d_all_mp):
        Thp += d ** 2 * mp.exp(2j * mp.pi * h)
    Thm = Thp.conjugate()
    absThp = abs(Thp)
    ThpThm = (Thp * Thm).real
    gate_pass = abs(absThp - D_mp) < mp.mpf(10) ** (-40)  # generous; true fails are O(1)-O(100)
    arg = mp.arg(Thp)
    c_eff = (4 * arg / mp.pi) % 8

    results_A[label] = dict(
        D2_total_exact=D2_total,
        D2_total_mp=D2_mp,
        D_mp=D_mp,
        Thp=Thp,
        absThp=absThp,
        ThpThm=ThpThm,
        gate_pass=gate_pass,
        c_eff=c_eff,
    )

    print(f"--- Candidate {label} ---")
    print(f"  D^2_total (exact)      = {D2_total}")
    print(f"  D^2_total (60dp check) = {mp.nstr(D2_mp, 40)}")
    print(f"  D  (=sqrt D^2)         = {mp.nstr(D_mp, 40)}")
    print(f"  Theta_+                = {mp.nstr(Thp, 30)}")
    print(f"  |Theta_+|              = {mp.nstr(absThp, 40)}")
    print(f"  Theta_+ * Theta_-      = {mp.nstr(ThpThm, 40)}")
    print(f"  D^2 (for comparison)   = {mp.nstr(D2_mp, 40)}")
    print(f"  gate |Theta_+| == D ?  -> diff = {mp.nstr(absThp - D_mp, 30)}  PASS={gate_pass}")
    print(f"  arg(Theta_+) (deg)     = {mp.nstr(arg * 180 / mp.pi, 20)}")
    print(f"  c_eff mod 8 (raw phase readout, only meaningful if gate passes) = {mp.nstr(c_eff, 20)}")
    print()

# ----------------------------------------------------------------------
# Question B: induced U(1) level
# ----------------------------------------------------------------------
print("=" * 78)
print("Question B: induced abelian CS level k_ind")
print("=" * 78)

targets = {
    "360": sp.Integer(360),
    "360/phi^2": sp.radsimp(sp.simplify(360 / phi ** 2)),
    "90": sp.Integer(90),
    "40": sp.Integer(40),
}
print("Targets (exact):")
for k, v in targets.items():
    print(f"  {k:>10} = {v}  = {sp.N(v, 20)}")
print()

results_B = {}

for label, d_tw in CANDIDATES.items():
    h_all = h_un + h_tw
    d_all = d_un + d_tw
    D2 = sp.radsimp(sp.simplify(sum(d ** 2 for d in d_all)))

    row = {}
    for qname, qfun in [("q=h", lambda h: h), ("q=h-c/24", lambda h: h - c24)]:
        S2 = sp.radsimp(sp.simplify(sum(d ** 2 * qfun(h) ** 2 for d, h in zip(d_all, h_all))))
        S2n = sp.radsimp(sp.simplify(S2 / D2))
        matches_raw = [name for name, t in targets.items() if sp.simplify(S2 - t) == 0]
        matches_norm = [name for name, t in targets.items() if sp.simplify(S2n - t) == 0]
        row[qname] = dict(S2=S2, S2n=S2n, matches_raw=matches_raw, matches_norm=matches_norm)
        print(f"[{label}, {qname}]")
        print(f"   S2 (raw)     = {S2}   = {sp.N(S2, 20)}   exact target match: {matches_raw or 'none'}")
        print(f"   S2/D2 (norm) = {S2n}  = {sp.N(S2n, 20)}   exact target match: {matches_norm or 'none'}")
    results_B[label] = row
    print()

# A3 needs high-precision route for S2 as well (algebraic but not Q(sqrt5))
print("[A3 high-precision cross-check of S2, matching the exact sympy route above]")
d_A3_mp = [mp_d(d) for d in d_A3]
h_all_mp = mp_h_un + mp_h_tw
d_all_mp = mp_d_un + d_A3_mp
D2_A3_mp = sum(d ** 2 for d in d_all_mp)
for qname, qfun in [("q=h", lambda h: h), ("q=h-c/24", lambda h: h - mp.mpf(c24))]:
    S2 = sum(d ** 2 * qfun(h) ** 2 for d, h in zip(d_all_mp, h_all_mp))
    S2n = S2 / D2_A3_mp
    print(f"   [{qname}]  S2 = {mp.nstr(S2, 30)}   S2/D2 = {mp.nstr(S2n, 30)}")
print()

# ----------------------------------------------------------------------
# Full Z_360 abelian lattice, maximal-refinement reference
# ----------------------------------------------------------------------
print("=" * 78)
print("Full Z_360 abelian-lattice reference (all d=1, q_n = n/360)")
print("=" * 78)
N = 360
S2_360 = sp.Rational(sum(n ** 2 for n in range(N)), N ** 2)
D2_360 = sp.Integer(N)  # sum of d_n^2 = 360*1
S2_360_over_360 = sp.simplify(S2_360 / 360)
S2_360_over_D2 = sp.simplify(S2_360 / D2_360)  # identical to S2/360 since D2=360
print(f"  S2 = sum_{{n=0}}^{{359}} (n/360)^2 = {S2_360} = {sp.N(S2_360, 20)}")
print(f"  D^2 (this lattice) = {D2_360}")
print(f"  S2/360 = {S2_360_over_360} = {sp.N(S2_360_over_360, 20)}")
print(f"  S2/D^2 = {S2_360_over_D2} = {sp.N(S2_360_over_D2, 20)}  (identical to S2/360 here, D^2=360)")
matches = [name for name, t in targets.items() if sp.simplify(S2_360 - t) == 0]
matches_n = [name for name, t in targets.items() if sp.simplify(S2_360_over_360 - t) == 0]
print(f"  exact target match (raw): {matches or 'none'}")
print(f"  exact target match (normalized): {matches_n or 'none'}")
print()

# ----------------------------------------------------------------------
# Question C
# ----------------------------------------------------------------------
print("=" * 78)
print("Question C: does alpha_inv = 360/phi^2 read as k_ind * <U> with k_ind=360")
print("             for any (candidate, normalization) pair from B?")
print("=" * 78)
alpha_inv = sp.radsimp(sp.simplify(360 / phi ** 2))
print(f"  alpha_inv = 360/phi^2 = {alpha_inv} = {sp.N(alpha_inv, 20)}")
found_360 = []
for label, row in results_B.items():
    for qname, d in row.items():
        if sp.simplify(d["S2"] - 360) == 0:
            found_360.append((label, qname, "raw"))
        if sp.simplify(d["S2n"] - 360) == 0:
            found_360.append((label, qname, "normalized"))
if sp.simplify(S2_360 - 360) == 0:
    found_360.append(("full-Z360-lattice", "q=n/360", "raw"))
if sp.simplify(S2_360_over_360 - 360) == 0:
    found_360.append(("full-Z360-lattice", "q=n/360", "normalized"))
print(f"  (candidate, charge convention, normalization) giving k_ind = 360 EXACTLY: {found_360 or 'NONE'}")
print()
print("Done.")
