"""Shared constants and single-tube dielectric-profile kernel for the 3D
lattice homogenization build (Paper III, src_paper3/lattice3d).

Everything here is pinned to the pre-registered task specification: constants
are derived from their defining equations (never hand-copied), to 1e-14
where the defining equation is transcendental.
"""
import numpy as np
from scipy.optimize import brentq

PHI = (1.0 + np.sqrt(5.0)) / 2.0

# Target for the whole build.
T_TARGET = 112.5 / PHI**10

# --- Torus geometry -----------------------------------------------------
R0 = 3.0


def _find_C2star():
    def g(C):
        return C * (2.0 * C + 1.0) / ((C**2 + 1.0) * (3.0 * C - 1.0)) - 2.0**(4.0 / 3.0) / PHI**5
    return brentq(g, 1.01, 10.0, xtol=1e-14, rtol=1e-14)


C2_STAR = _find_C2star()
r0 = R0 / C2_STAR  # minor radius of the T(2,2) tube


# --- Dielectric kernel constants ----------------------------------------
def _find_bstar():
    def f(x):
        # x = sqrt(8b); solves (15/8) arctan(x)/x = phi
        return (15.0 / 8.0) * np.arctan(x) / x - PHI
    x_star = brentq(f, 1e-6, 10.0, xtol=1e-14, rtol=1e-14)
    return x_star**2 / 8.0


B_STAR = _find_bstar()
XSTAR2 = 8.0 * B_STAR  # dimensionless kernel amplitude, ~0.53721045

CSTAR2 = (3.0 / 4.0) * PHI**5 / 2.0**(4.0 / 3.0)
CSTAR = np.sqrt(CSTAR2)  # ~1.81682488, tube-width unit (Paper II)


def u_profile(rho):
    """u_i(rho) = xstar2 / (1+(rho/Cstar)^2)^2, the per-tube dimensionless
    dielectric kernel (beta*X_i in the validated 2D single-tube chain)."""
    return XSTAR2 / (1.0 + (rho / CSTAR) ** 2) ** 2


def eps_tan_of_U(U, eta):
    return 1.0 / (1.0 + U) ** 2 + 1j * eta


def eps_rad_of_U(U, eta):
    return (1.0 - 3.0 * U) / (1.0 + U) ** 3 + 1j * eta


# validated single-tube transverse polarizability (2D chain, eta-> extrapolated;
# see src_paper3/dielectric_cell/results.md, R1 beta=beta_star row) -- used only
# for the V3 chain-of-custody cross-check, not recomputed here.
ALPHA_V3 = -0.97107 + 0.0135j


if __name__ == "__main__":
    print(f"PHI = {PHI:.14f}")
    print(f"T_TARGET = {T_TARGET:.10f}")
    print(f"C2_STAR = {C2_STAR:.10f}  r0 = {r0:.10f}")
    print(f"B_STAR = {B_STAR:.14f}  XSTAR2 = {XSTAR2:.10f}")
    print(f"CSTAR2 = {CSTAR2:.10f}  CSTAR = {CSTAR:.10f}")
    print("cross-check against task-spec rounded literals:")
    print(f"  C2* spec=3.43182008 computed={C2_STAR:.8f} diff={abs(C2_STAR-3.43182008):.2e}")
    print(f"  r0  spec=0.874172   computed={r0:.6f}       diff={abs(r0-0.874172):.2e}")
    print(f"  xstar2 spec=0.53720965 computed={XSTAR2:.8f} diff={abs(XSTAR2-0.53720965):.2e}")
    print(f"  Cstar spec=1.81682488  computed={CSTAR:.8f}  diff={abs(CSTAR-1.81682488):.2e}")
