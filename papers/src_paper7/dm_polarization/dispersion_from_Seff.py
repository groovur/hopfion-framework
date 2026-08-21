"""Derive the condensate director dispersion from the semi-Dirac two-band model
and test whether its angular structure IS the Hopfion suppression S_eff = sin^4.

Semi-Dirac point: one axis quadratic (heavy), the perpendicular axis linear
(light). Write it with the heavy axis PERPENDICULAR to the radial preferred
direction e_r (the physical 'hard perpendicular, easy radial'):
    H = A k_perp^2 sigma_x + v k_par sigma_y,
with k_par = k cos(theta) (along e_r), k_perp = k sin(theta).  d = (d_x,d_y),
E^2 = d_x^2 + d_y^2.  We check the coefficient of the heavy (k^4) term against
S_eff = sin^4(theta), and compare with the paper's written form
E^2 = A^2 k^4 cos^4 + v^2 k^2 sin^2 (heavy ALONG e_r).
"""
import sympy as sp

A, v, k, th = sp.symbols('A v k theta', positive=True)

# heavy axis PERPENDICULAR to e_r:
kpar, kperp = k*sp.cos(th), k*sp.sin(th)          # par = radial, perp = transverse
dx = A*kperp**2                                    # quadratic (heavy) band, perp
dy = v*kpar                                        # linear (light) band, radial
E2 = sp.expand_trig(sp.simplify(dx**2 + dy**2))
print("heavy-perp model:  E^2 =", E2)

heavy_coeff = sp.simplify(dx**2 / (A**2 * k**4))   # angular factor of the k^4 term
light_coeff = sp.simplify(dy**2 / (v**2 * k**2))   # angular factor of the k^2 term
print("  heavy (k^4) angular factor =", heavy_coeff, "   (S_eff = sin^4?)",
      sp.simplify(heavy_coeff - sp.sin(th)**4) == 0)
print("  light (k^2) angular factor =", light_coeff)

# the '(sin^2)^2' structure: the heavy band d_x is itself quadratic in k_perp,
# d_x = A k^2 sin^2(theta); the energy squares it -> sin^4.
print("  d_x =", sp.simplify(dx), "  ~ sin^2(theta) ;  d_x^2 ~ sin^4(theta) = (sin^2)^2")

# physical limits
print("  theta=0 (radial):   E^2 =", sp.simplify(E2.subs(th,0)),   "-> linear/light/EASY")
print("  theta=pi/2 (perp):  E^2 =", sp.simplify(E2.subs(th,sp.pi/2)), "-> quadratic/heavy/HARD")

print("\n--- paper's written form (heavy ALONG e_r): E^2 = A^2 k^4 cos^4 + v^2 k^2 sin^2 ---")
E2_paper = A**2*k**4*sp.cos(th)**4 + v**2*k**2*sp.sin(th)**2
paper_heavy = sp.simplify(E2_paper.coeff(A**2)/k**4)
print("  heavy angular factor =", paper_heavy, " -> equals S_eff=sin^4 ?",
      sp.simplify(paper_heavy - sp.sin(th)**4) == 0)
print("  theta=0 (radial):", sp.simplify(E2_paper.subs(th,0)), "-> HEAVY at radial (contradicts 'easy radial')")

print("\nVERDICT: S_eff = sin^4(theta) IS the heavy-band angular factor of the")
print("semi-Dirac dispersion -- but ONLY with the heavy axis PERPENDICULAR to e_r")
print("(E^2 = A^2 k^4 sin^4 + v^2 k^2 cos^2). It is intrinsic to the band structure")
print("(the transverse quadratic band k_perp^2, squared in E^2 -> (sin^2)^2 = sin^4),")
print("no external field. The paper's 'heavy ALONG e_r' / cos^4 form is the SAME")
print("dispersion with the angle measured from the other axis, and its wording is")
print("backwards: it must read heavy=perpendicular, light=radial to match S_eff.")
