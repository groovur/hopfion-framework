#!/usr/bin/env python3
"""
jacobi_bvp_solver.py
====================
Solves the Sturm-Liouville BVP for the thick-torus profile correction
f2(rho) of the density-feedback Faddeev-Niemi Hopfion.

BVP (Paper XIV, Theorem 4.1, eq. 4.5):
    L0[f2](rho) = S2(rho),    f2(0) = f2(inf) = 0

where:
    L0[f2] = f2'' - V''(f0) f2
    V''(f0) = C*^2 cos(2f0) / rho^2           (Jacobi potential)
    S2(rho) = (C*/(2R0^2)) sin(f0)(1-C*cos(f0))   (toroidal curvature source)
    f0(rho) = 2 arctan(rho^{-C*})              (modified BPS profile)

Discretisation (Paper XIV, Remark 4.5, eq. 4.13):
    Finite differences in log-space t = log(rho).
    Setting A_i = C*^2 cos(2f_{0,i}) and B_i = rho_i^2 S2(rho_i),
    the scheme at interior point i is:

        (1/dt^2 - 1/(2dt)) f_{i+1}
      + (-2/dt^2 - A_i)    f_i
      + (1/dt^2 + 1/(2dt)) f_{i-1}  = B_i

    forming a tridiagonal system with N = 8000 interior points
    on t in [-12, 12] (rho in [e^{-12}, e^{12}], dt ~ 0.003).

Parameters:
    C*  = 3.4318   (modified BPS exponent, Paper I Theorem 7.3)
    R0  = 3        (torus major radius in condensate units)
    phi = (1+sqrt(5))/2

Output:
    Table of f0(rho) and f2(rho) at representative rho values,
    reproducing Paper XIV Remark 4.5.

Reference:
    F. Manfredi, "The Thick-Torus Profile Correction to the
    Density-Feedback Hopfion" (Paper XIV), Section 4.

Usage:
    python jacobi_bvp_solver.py
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


# ═══════════════════════════════════════════════════════════════════
# Physical constants
# ═══════════════════════════════════════════════════════════════════
PHI = (1 + np.sqrt(5)) / 2          # golden ratio
CSTAR = 3.4318                       # modified BPS exponent
R0 = 3.0                             # torus major radius


# ═══════════════════════════════════════════════════════════════════
# Profile functions
# ═══════════════════════════════════════════════════════════════════
def f0(rho):
    """Modified BPS profile: f0(rho) = 2 arctan(rho^{-C*})."""
    return 2.0 * np.arctan(rho**(-CSTAR))


def sin_f0(rho):
    """sin(f0) = 2 rho^{C*} / (rho^{2C*} + 1)."""
    t = rho**(2 * CSTAR)
    return 2.0 * rho**CSTAR / (t + 1)


def cos_f0(rho):
    """cos(f0) = (rho^{2C*} - 1) / (rho^{2C*} + 1)."""
    t = rho**(2 * CSTAR)
    return (t - 1) / (t + 1)


def cos_2f0(rho):
    """cos(2f0) = 1 - 2 sin^2(f0)."""
    s = sin_f0(rho)
    return 1.0 - 2.0 * s * s


# ═══════════════════════════════════════════════════════════════════
# BVP components
# ═══════════════════════════════════════════════════════════════════
def jacobi_potential_times_rho2(rho):
    """A(rho) = C*^2 cos(2f0(rho)).  [V''(f0) = A/rho^2]"""
    return CSTAR**2 * cos_2f0(rho)


def source_times_rho2(rho):
    """B(rho) = rho^2 S2(rho) = rho^2 (C*/(2R0^2)) sin(f0)(1 - C* cos(f0))."""
    return rho**2 * (CSTAR / (2.0 * R0**2)) * sin_f0(rho) * (1.0 - CSTAR * cos_f0(rho))


# ═══════════════════════════════════════════════════════════════════
# Solver
# ═══════════════════════════════════════════════════════════════════
def solve_bvp(N=8000, t_min=-12.0, t_max=12.0):
    """
    Solve L0[f2] = S2 by finite differences in log-space t = log(rho).

    Parameters
    ----------
    N : int
        Number of interior grid points.
    t_min, t_max : float
        Endpoints of the log-rho grid (rho_min = e^{t_min}, rho_max = e^{t_max}).

    Returns
    -------
    rho : ndarray, shape (N+2,)
        Grid points in rho-space (including boundaries).
    f2 : ndarray, shape (N+2,)
        Profile correction f2(rho) at each grid point.
    """
    # Uniform grid in t = log(rho)
    t_grid = np.linspace(t_min, t_max, N + 2)
    dt = t_grid[1] - t_grid[0]
    rho_grid = np.exp(t_grid)

    # Evaluate coefficients at interior points (indices 1..N)
    rho_int = rho_grid[1:-1]
    A = jacobi_potential_times_rho2(rho_int)
    B = source_times_rho2(rho_int)

    # Tridiagonal matrix coefficients (Paper XIV eq. 4.13)
    coeff_upper = 1.0 / dt**2 - 1.0 / (2.0 * dt)     # f_{i+1} coefficient
    coeff_diag = -2.0 / dt**2 - A                      # f_i coefficient
    coeff_lower = 1.0 / dt**2 + 1.0 / (2.0 * dt)     # f_{i-1} coefficient

    # Build sparse tridiagonal system
    diag_main = coeff_diag
    diag_upper = np.full(N - 1, coeff_upper)
    diag_lower = np.full(N - 1, coeff_lower)

    matrix = sparse.diags(
        [diag_lower, diag_main, diag_upper],
        offsets=[-1, 0, 1],
        shape=(N, N),
        format='csc'
    )

    # Solve with boundary conditions f2(t_min) = f2(t_max) = 0
    f2_interior = spsolve(matrix, B)

    # Assemble full solution
    f2_full = np.zeros(N + 2)
    f2_full[1:-1] = f2_interior

    return rho_grid, f2_full


# ═══════════════════════════════════════════════════════════════════
# Verification
# ═══════════════════════════════════════════════════════════════════
def verify_solution(rho, f2):
    """Check BVP residual and asymptotic exponents."""
    t = np.log(rho)
    dt = t[1] - t[0]
    N = len(rho)

    # Reconstruct f2'' from the finite-difference stencil
    # In log-space: f2''(rho) = (g''(t) - g'(t)) / rho^2
    # where g(t) = f2(e^t)
    g = f2.copy()
    residuals = []
    for i in range(1, N - 1):
        gpp = (g[i + 1] - 2 * g[i] + g[i - 1]) / dt**2
        gp = (g[i + 1] - g[i - 1]) / (2 * dt)
        f2pp = (gpp - gp) / rho[i]**2
        V_pp = CSTAR**2 * cos_2f0(rho[i]) / rho[i]**2
        S2 = (CSTAR / (2 * R0**2)) * sin_f0(rho[i]) * (1 - CSTAR * cos_f0(rho[i]))
        lhs = f2pp - V_pp * f2[i]
        residuals.append(abs(lhs - S2))

    residuals = np.array(residuals)

    # Asymptotic exponents
    p_plus = (1 + np.sqrt(1 + 4 * CSTAR**2)) / 2
    p_minus = (1 - np.sqrt(1 + 4 * CSTAR**2)) / 2

    # Near rho=0: f2 ~ rho^{p+}
    mask_small = (rho > 0.01) & (rho < 0.03)
    if np.any(mask_small) and np.all(np.abs(f2[mask_small]) > 0):
        idx = np.where(mask_small)[0]
        log_rho_s = np.log(rho[idx])
        log_f2_s = np.log(np.abs(f2[idx]))
        slope_small = np.polyfit(log_rho_s, log_f2_s, 1)[0]
    else:
        slope_small = np.nan

    # Near rho=inf: f2 ~ rho^{2-C*}
    mask_large = (rho > 50) & (rho < 200)
    if np.any(mask_large) and np.all(np.abs(f2[mask_large]) > 0):
        idx = np.where(mask_large)[0]
        log_rho_l = np.log(rho[idx])
        log_f2_l = np.log(np.abs(f2[idx]))
        slope_large = np.polyfit(log_rho_l, log_f2_l, 1)[0]
    else:
        slope_large = np.nan

    return {
        'max_residual': np.max(residuals[10:-10]),
        'rms_residual': np.sqrt(np.mean(residuals[10:-10]**2)),
        'slope_small': slope_small,
        'slope_small_theory': p_plus,
        'slope_large': slope_large,
        'slope_large_theory': 2 - CSTAR,
        'p_plus': p_plus,
        'p_minus': p_minus,
    }


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════
def main():
    print("Jacobi BVP Solver for the Thick-Torus Profile Correction")
    print("=" * 60)
    print(f"  C*  = {CSTAR}")
    print(f"  R0  = {R0}")
    print(f"  phi = {PHI:.10f}")
    print()

    # Solve
    rho, f2 = solve_bvp(N=8000, t_min=-12, t_max=12)

    # Verify
    v = verify_solution(rho, f2)
    print("Verification:")
    print(f"  Max BVP residual (interior):  {v['max_residual']:.2e}")
    print(f"  RMS BVP residual (interior):  {v['rms_residual']:.2e}")
    print(f"  Small-rho exponent: {v['slope_small']:.3f}  (theory: {v['slope_small_theory']:.3f})")
    print(f"  Large-rho exponent: {v['slope_large']:.3f}  (theory: {v['slope_large_theory']:.3f})")
    print()

    # Key properties
    peak_idx = np.argmax(np.abs(f2))
    sign_changes = np.where(np.diff(np.sign(f2)))[0]
    sign_change_rho = []
    for sc in sign_changes:
        if abs(f2[sc]) > 1e-12:
            r_cross = rho[sc] - f2[sc] * (rho[sc + 1] - rho[sc]) / (f2[sc + 1] - f2[sc])
            sign_change_rho.append(r_cross)

    print("Key properties:")
    print(f"  Peak |f2|:      {np.abs(f2[peak_idx]):.4f} at rho = {rho[peak_idx]:.4f}")
    print(f"  f2 < 0 for rho < {sign_change_rho[0]:.2f} (core flattening)")
    print(f"  f2 > 0 for rho > {sign_change_rho[0]:.2f} (wing extension)")
    print()

    # Table (Paper XIV Remark 4.5)
    rho_table = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    labels = [
        "deep core, near f=pi",
        "inner core",
        "inner core",
        "approaching transition",
        "transition point, f=pi/2",
        "outer transition",
        "outer region",
        "sign change at rho~2.5",
        "wing extension",
    ]

    print(f"{'rho':>8s}  {'f0(rho)':>10s}  {'f2(rho)':>14s}  {'meaning'}")
    print(f"{'-'*8}  {'-'*10}  {'-'*14}  {'-'*30}")
    for r, label in zip(rho_table, labels):
        idx = np.argmin(np.abs(rho - r))
        f0_val = f0(rho[idx])
        f2_val = f2[idx]
        print(f"{rho[idx]:8.1f}  {f0_val:10.4f}  {f2_val:+14.2e}  {label}")

    print()
    print(f"Indicial exponents: p_+ = {v['p_plus']:.4f}, p_- = {v['p_minus']:.4f}")
    print(f"Wronskian: W = p_+ - p_- = sqrt(1+4C*^2) = {v['p_plus'] - v['p_minus']:.4f}")


if __name__ == '__main__':
    main()
