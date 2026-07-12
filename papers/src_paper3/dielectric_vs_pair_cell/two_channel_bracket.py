"""Two-channel (polarization x vertex-suppression) bracket for Delta_1.

Computes aInv_obs/aInv_geo = eps_eff * <1/(1+b X)>^{-1} for a contact-packed
condensate of Q_H=2 two-tube (T(2,2)) hopfions and compares with the
golden/geo target 112.5/phi^10 = 0.914695.

All inputs are fixed by the framework; there is no free density:
  b*   : V(b*) = phi, V(b) = (15/8) arctan(sqrt(8b))/sqrt(8b)  [Paper III]
  C*^2 = (3/4) phi^5 / 2^(4/3)  (thin-torus tube-width unit)   [Paper II]
  contact spacing 2(R0+1) = 8 condensate units, R0 = 3
  <alpha_pair> : 2D anisotropic-dielectric solver, pair separation
  s = 2 R0/(C2* C*) = 0.96231 solver units (pair_cell/, h=0.1, eta=1e-2;
  real parts grid-converged at the ~1% level, imaginary parts not
  converged there and excluded).

Vertex channel: <1/(1+bX)> = 1 - n_tube * D with per-tube deficit
D = 2 pi b * Integral_0^inf X/(1+bX) t dt = 2 pi b * 4 arctan(x*)/x*,
x* = sqrt(8 b*) (same integral as J_fb, Paper III Prop. virial_exact).
The vertex normalization (exponent 1 on the suppression average) is an
assumption of this bracket, not derived from the Paper IV vertex.
"""
import numpy as np
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R0 = 3.0

b_star = brentq(lambda b: (15/8)*np.arctan(np.sqrt(8*b))/np.sqrt(8*b) - PHI,
                1e-6, 10.0, xtol=1e-14, rtol=1e-15)
x_star = np.sqrt(8*b_star)
Cstar2 = (3/4)*PHI**5 / 2**(4/3)
C2star = brentq(lambda C: C*(2*C+1)/((C*C+1)*(3*C-1)) - 2**(4/3)/PHI**5,
                1.0, 10.0, xtol=1e-14)
s_pair = 2*R0/(C2star*np.sqrt(Cstar2))
target = 112.5/PHI**10

# <alpha_pair> (real part) from pair_cell/results.md, h=0.1, eta=1e-2
alpha_pair = {"A": -1.850339, "B": -2.952718}

# contact packing of pairs, converted to solver units (1 solver unit = C*)
n_pair = {"square": (1/(2*(R0+1))**2) * Cstar2,
          "triangular": (2/(np.sqrt(3)*(2*(R0+1))**2)) * Cstar2}

D = 2*np.pi*b_star * 4*np.arctan(x_star)/x_star   # vertex deficit per tube

print(f"b* = {b_star:.8f}   x* = {x_star:.8f}   C*^2 = {Cstar2:.6f}")
print(f"C2* = {C2star:.8f}  s_pair = {s_pair:.8f}  D = {D:.5f}")
print(f"target golden/geo = 112.5/phi^10 = {target:.6f}\n")
for lat, npair in n_pair.items():
    vertex = 1/(1 - 2*npair*D)
    for m, a in alpha_pair.items():
        t = npair*a
        eps_eff = 1 + 2*t/(1 - t)          # 2D Maxwell-Garnett
        prod = eps_eff*vertex
        print(f"{lat:10s} model {m}: eps_eff={eps_eff:.5f}  "
              f"vertex={vertex:.5f}  product={prod:.5f}  "
              f"dev={(prod/target-1)*100:+.2f}%")
