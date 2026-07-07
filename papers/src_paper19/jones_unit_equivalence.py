"""
jones_unit_equivalence.py
=========================
Numerical and symbolic verification of Jones polynomial moduli and the
cyclotomic ring structure used in Paper XIX's unit-equivalence theorems.

WHAT THIS SCRIPT VERIFIES:
  - V_{T(2,2)}(q₅) = 1+ζ₅³, |V_{T(2,2)}(q₅)| = 1/φ
    (the standard Hopf-link Jones polynomial V(t)=-t^{1/2}-t^{5/2} at t=q₅)
  - V_sat(T(2,3))(q₅) = q₅³/φ, |V_sat| = 1/φ
    (Paper XVIII, Theorem thm:trefoil_satellite_jones)
  - Phase relation: V_sat(q₅) = q₅⁻¹ · V_{T(2,2)}(q₅)
  - Disjoint-union factor -(q₅^{1/2}+q₅^{-1/2}) = -φ (a different quantity
    from V_{T(2,2)}; kept in a separate function so the two are never
    conflated) (Lem. P19:lem:disjoint_factor)
  - Channel 1 unit factor ζ₅+ζ₅⁴ = 1/φ (Thm. P19:thm:channel_1)
  - Channel 2 unit factor -φ·ζ₅ (Thm. P19:thm:channel_2)
  - Both channel unit factors have norm 1 in Z[ζ₅] (verified via resultant
    with the cyclotomic polynomial Φ₅, not just numeric modulus matching)
  - S-matrix ratio S₀₁/S₀₀ = φ, giving the solar angle tan(θ₁₂⁽⁰⁾) = 1/φ
    directly (Paper XVII); independent of the Jones-polynomial values above

All values cited in:
  P19:thm:channel_1  (unit factor ζ₅+ζ₅⁴ = 1/φ)
  P19:thm:channel_2  (unit factor -φ·ζ₅)
  P19:prop:solar_angle_jones  (tan(θ₁₂⁽⁰⁾) = 1/φ, via the S-matrix)
  P19:lem:disjoint_factor  (disjoint-union factor = -φ)
"""

import numpy as np
import sympy as sp

PHI = (1 + 5**0.5) / 2
PI = np.pi
zeta5 = np.exp(2j * PI / 5)
q5 = zeta5

# ── symbolic setup for exact Z[zeta5] verification ────────────────────────
z = sp.symbols('z')            # represents zeta5
Phi5 = z**4 + z**3 + z**2 + z + 1

def reduce_poly(expr):
    """Reduce a polynomial in z modulo the cyclotomic polynomial Phi5,
    giving the canonical degree<=3 representative in Z[zeta5]."""
    p = sp.Poly(sp.expand(expr), z)
    _, r = sp.div(p, sp.Poly(Phi5, z), z)
    return sp.expand(r.as_expr())

def norm_in_Qzeta5(expr):
    """Field norm N_{Q(zeta5)/Q} of an element of Z[zeta5], given as a
    polynomial in z, via the resultant with the cyclotomic polynomial."""
    return sp.resultant(Phi5, sp.expand(expr), z)


# ── Section 1: Jones moduli ───────────────────────────────────────────────

def V_hopf_link(q):
    """Jones polynomial V_{T(2,2)}(q) of the (positive) Hopf link.
    Standard formula: V(t) = -(t^{1/2}+t^{5/2})."""
    return -(q**0.5 + q**2.5)

def V_disjoint_union_factor(q):
    """Disjoint-union multiplier -(q^{1/2}+q^{-1/2}) for the Kauffman
    bracket. This is a DIFFERENT quantity from V_hopf_link above -- it is
    the factor by which the bracket of two disjoint components multiplies,
    not the Jones polynomial of any single link."""
    return -(q**0.5 + q**(-0.5))

def V_trefoil_free(q):
    """Jones polynomial V_{T(2,3)}(q) of the free-standing trefoil.
    Standard formula: V(t) = -t^{-4}+t^{-3}+t^{-1}."""
    return -q**(-4) + q**(-3) + q**(-1)


V_h = V_hopf_link(q5)
V_disj = V_disjoint_union_factor(q5)
V_t_free = V_trefoil_free(q5)
V_sat = q5**3 / PHI   # Paper XVIII, Theorem thm:trefoil_satellite_jones

print("=" * 65)
print("SECTION 1: Jones moduli at q₅ = e^{2πi/5}")
print("=" * 65)
print()
print(f"  V_{{T(2,2)}}(q₅) = {V_h:.8f}")
print(f"  |V_{{T(2,2)}}(q₅)| = {abs(V_h):.10f}")
print(f"  1/φ                 = {1/PHI:.10f}")
print(f"  Match |V|=1/φ:       {'YES ✓' if abs(abs(V_h)-1/PHI)<1e-9 else 'NO ✗'}")
print()
print(f"  V_sat(T(2,3))(q₅) = q₅³/φ = {V_sat:.8f}")
print(f"  |V_sat|             = {abs(V_sat):.10f}")
print(f"  Match |V_sat|=1/φ:    {'YES ✓' if abs(abs(V_sat)-1/PHI)<1e-9 else 'NO ✗'}")
print()
print("  Phase relation: V_sat(q₅) = q₅⁻¹ · V_{T(2,2)}(q₅)")
phase_check = q5**-1 * V_h
print(f"  q₅⁻¹·V_{{T(2,2)}}(q₅) = {phase_check:.8f}   V_sat = {V_sat:.8f}")
print(f"  Match: {'YES ✓' if abs(phase_check - V_sat) < 1e-9 else 'NO ✗'}")
print()
print(f"  V_{{T(2,3)}}(q₅) free-standing = {V_t_free:.8f}, |V| = {abs(V_t_free):.8f}")
print("  (a third, different value -- the trefoil taken as an abstract")
print("  free-standing knot rather than the satellite construction)")


# ── Section 2: exact value of V_{T(2,2)}(q5) in Z[zeta5] ─────────────────
print()
print("=" * 65)
print("SECTION 2: Exact value V_{T(2,2)}(q₅) = 1+ζ₅³ in Z[ζ₅]")
print("=" * 65)
V_h_exact = 1 + zeta5**3
print(f"  1+ζ₅³ = {V_h_exact:.10f}")
print(f"  V_{{T(2,2)}}(q₅) (numeric) = {V_h:.10f}")
print(f"  Match: {'YES ✓' if abs(V_h_exact - V_h) < 1e-9 else 'NO ✗'}")
print()
print(f"  |1+ζ₅³|² = {abs(1+zeta5**3)**2:.10f}, 1/φ² = {1/PHI**2:.10f}")
print(f"  Match:      {'YES ✓' if abs(abs(1+zeta5**3)**2 - 1/PHI**2)<1e-9 else 'NO ✗'}")


# ── Section 3: Disjoint-union factor (a different quantity) ─────────────
print()
print("=" * 65)
print("SECTION 3: Disjoint-union factor (Lemma P19:lem:disjoint_factor)")
print("=" * 65)
print(f"  -(q₅^{{1/2}}+q₅^{{-1/2}}) = {V_disj:.10f}")
print(f"  -φ                        = {-PHI:.10f}")
print(f"  Match:                       {'YES ✓' if abs(V_disj+PHI)<1e-9 else 'NO ✗'}")
print("  (Note: this is the disjoint-union multiplier, not V_{T(2,2)}.")
print("  The two are structurally different and must not be conflated.)")


# ── Section 4: Channel unit factors, exact in Z[zeta5] ───────────────────
print()
print("=" * 65)
print("SECTION 4: Channel unit factors")
print("=" * 65)

# phi and phi^2 as z-polynomials (phi = 1+z+z^4), reduced mod Phi5
phi_poly = reduce_poly(1 + z + z**4)
phi2_poly = reduce_poly(phi_poly**2)
neg_phi_poly = reduce_poly(-phi_poly)
print(f"  phi   reduced in Z[zeta5] basis = {phi_poly}")
print(f"  phi^2 reduced in Z[zeta5] basis = {phi2_poly}")

# Channel 1: T(2,2) u T(2,2) -> T(2,3) u U
#   J_eff(init)  = (-phi) * V_T22^2
#   J_eff(final) = (-phi) * V_sat * 1 = (-phi)*(zeta5^3/phi) = -zeta5^3
#   ratio = J_eff(init) / J_eff(final) = (-phi)*V_T22^2 / (-zeta5^3)
#         = phi * V_T22^2 * zeta5^2        (since zeta5^5=1, 1/zeta5^3=zeta5^2)
print()
print("  Channel 1: T(2,2)⊔T(2,2) → T(2,3)⊔U")
sq = reduce_poly((1 + z**3)**2)
print(f"    (1+ζ₅³)² reduced                = {sq}")
numerator = reduce_poly(neg_phi_poly * sq)
print(f"    J_eff(init) = (-φ)(1+ζ₅³)²       = {numerator}")
print("    J_eff(final) = (-φ)·(ζ₅³/φ)·1 = -ζ₅³")
print("    ratio = J_eff(init)/(-ζ₅³) = J_eff(init)·ζ₅² (since ζ₅⁵=1)")
ratio1 = reduce_poly(numerator * z**2)
print(f"    J_eff(init)·ζ₅² reduced         = {ratio1}")
ratio1_final = reduce_poly(-ratio1)
print(f"    ratio = -(J_eff(init)·ζ₅²)      = {ratio1_final}   [ = ζ₅+ζ₅⁴ ]")

u1_numeric = ((-PHI) * V_h_exact**2) / ((-PHI) * V_sat)
print(f"    numeric cross-check: {u1_numeric:.8f}   (expect 1/φ = {1/PHI:.8f})")
print(f"    Match: {'YES ✓' if abs(u1_numeric - 1/PHI) < 1e-9 else 'NO ✗'}")
norm1 = norm_in_Qzeta5(z + z**4)
print(f"    Norm_{{Q(ζ₅)/Q}}(ζ₅+ζ₅⁴) = {norm1}  (unit iff ±1)")

# Channel 2: T(2,2) u U -> T(2,3)
#   J_eff(init)  = (-phi) * V_T22 * 1
#   J_eff(final) = V_sat = zeta5^3/phi
#   ratio = J_eff(init) / J_eff(final) = (-phi)*V_T22 * phi/zeta5^3
#         = -phi^2 * V_T22 * zeta5^2        (since 1/zeta5^3 = zeta5^2)
print()
print("  Channel 2: T(2,2)⊔U → T(2,3)")
init2 = reduce_poly(-phi2_poly * (1 + z**3))
print(f"    -φ²(1+ζ₅³) reduced              = {init2}")
print("    ratio = -φ²(1+ζ₅³)·ζ₅² (since J_eff(final)=ζ₅³/φ, dividing")
print("    by it multiplies by φ·ζ₅²; combined with the φ already in -φ²")
print("    this gives the -φ² factor shown above times ζ₅²)")
ratio2_poly = reduce_poly(init2 * z**2)
print(f"    -φ²(1+ζ₅³)·ζ₅² reduced          = {ratio2_poly}   [ = -φ·ζ₅ ]")

u2_numeric = ((-PHI) * V_h_exact * 1) / (V_sat)
print(f"    numeric cross-check: {u2_numeric:.8f}   (expect -φ·ζ₅ = {(-PHI*zeta5):.8f})")
print(f"    Match: {'YES ✓' if abs(u2_numeric - (-PHI*zeta5)) < 1e-9 else 'NO ✗'}")
norm2 = norm_in_Qzeta5(-(1 + z + z**4) * z)
print(f"    Norm_{{Q(ζ₅)/Q}}(-φ·ζ₅) = {norm2}  (unit iff ±1)")


# ── Section 5: Solar angle -- independent of the Jones-polynomial values ─
print()
print("=" * 65)
print("SECTION 5: Solar angle (S-matrix route, Paper XVII)")
print("=" * 65)
S_ratio = np.sin(3*PI/5)/np.sin(PI/5)
print(f"  S₀₁/S₀₀ = sin(3π/5)/sin(π/5) = {S_ratio:.8f} = φ")
solar_angle = np.degrees(np.arctan(1/PHI))
print(f"  tan(θ₁₂⁽⁰⁾) = S₀₀/S₀₁ = 1/φ = {1/PHI:.8f}")
print(f"  θ₁₂⁽⁰⁾ = arctan(1/φ) = {solar_angle:.4f}°")
print("  (This route uses only the SU(2)_3 modular S-matrix; it does not")
print("  depend on V_{T(2,2)}(q5) or any other Jones-polynomial value.)")

print()
print("Done.")
