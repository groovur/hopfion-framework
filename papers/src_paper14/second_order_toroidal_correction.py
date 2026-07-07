#!/usr/bin/env python3
"""
second_order_toroidal_correction.py
====================================
Computes the O(1/R0^4) profile correction to the thick-torus
ratio r_3D = J4^{3D}/J2a^{3D} from the Battye-Sutcliffe profile
perturbation f2(rho)*cos(2chi)/R0^2.

Method:
  1. Solve the Sturm-Liouville BVP  L0[f2] = S2(rho)
     with Jacobi potential V''(f0) and toroidal curvature source S2.
  2. Compute the four Paper XIV integrals (I4, I2a, I4_tor, I2a_tor)
     by full 2D (rho, chi) Gauss quadrature with f = f0 + f2*cos2chi/R0^2.
  3. Assemble J4^{3D}, J2a^{3D} via the cross-coupled toroidal formula
     and extract the corrected ratio r_3D.

Reference: Paper XIV, Section 6 (Second-order toroidal correction).

Usage:
  python second_order_toroidal_correction.py [--R0 3.0] [--Cstar 3.4318] [--N 8000] [--Nchi 64]

Author: F. Manfredi
"""
import argparse
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq


# ═══════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════
PHI = (1 + np.sqrt(5)) / 2
R_WZW = 2**(4/3) / PHI**5           # 0.22721...
QT2 = PHI**8 + 1                     # 47.979...


# ═══════════════════════════════════════════════════════════════════
# Profile functions for f0(rho) = 2*arctan(rho^{-C*})
# satisfying the modified BPS equation f0' = -C* sin(f0)/rho
# ═══════════════════════════════════════════════════════════════════
def f0_val(rho, Cs):
    return 2 * np.arctan(rho**(-Cs))

def f0_prime(rho, Cs):
    return -2 * Cs * rho**(Cs - 1) / (rho**(2*Cs) + 1)

def sin_f0(rho, Cs):
    return 2 * rho**Cs / (rho**(2*Cs) + 1)

def cos_f0(rho, Cs):
    return (rho**(2*Cs) - 1) / (rho**(2*Cs) + 1)

def cos_2f0(rho, Cs):
    s = sin_f0(rho, Cs)
    return 1 - 2 * s**2


# ═══════════════════════════════════════════════════════════════════
# BVP: L0[f2] = S2(rho)
# L0[f2] = f2'' - V''(f0)*f2,  V''(f0) = C*^2 cos(2f0)/rho^2
# S2(rho) = (C*/(2R0^2)) sin(f0) (1 - C* cos(f0))
# ═══════════════════════════════════════════════════════════════════
def solve_bvp(rho_grid, Cs, R0):
    """Solve L0[f2] = S2 on a log-spaced grid. Returns (f2, f2')."""
    N = len(rho_grid)
    log_rho = np.log(rho_grid)
    h = log_rho[1] - log_rho[0]  # uniform in log(rho)

    # Coefficients in log-space: g'' - g' - rho^2 V'' g = rho^2 S2
    V_pp = Cs**2 * cos_2f0(rho_grid, Cs) / rho_grid**2
    S2 = (Cs / (2 * R0**2)) * sin_f0(rho_grid, Cs) * (1 - Cs * cos_f0(rho_grid, Cs))

    n_int = N - 2
    diag_a = np.full(n_int, 1/h**2 + 1/(2*h))
    diag_b = -2/h**2 - rho_grid[1:-1]**2 * V_pp[1:-1]
    diag_c = np.full(n_int, 1/h**2 - 1/(2*h))
    rhs = rho_grid[1:-1]**2 * S2[1:-1]

    A = sparse.diags([diag_a[1:], diag_b, diag_c[:-1]], [-1, 0, 1],
                      shape=(n_int, n_int), format='csc')
    g = np.zeros(N)
    g[1:-1] = spsolve(A, rhs)

    # f2 = g (in log-rho space, f2(rho) = g(log rho))
    f2 = g.copy()

    # f2'(rho) = g'(u)/rho
    gp = np.zeros(N)
    for i in range(1, N-1):
        gp[i] = (g[i+1] - g[i-1]) / (2*h)
    f2p = gp / rho_grid

    return f2, f2p


# ═══════════════════════════════════════════════════════════════════
# Azimuthally-averaged integrals via Gauss quadrature in chi
# ═══════════════════════════════════════════════════════════════════
def compute_integrals(rho_grid, Cs, R0, f2=None, f2p=None, N_chi=64):
    """
    Compute I4, I2a, I4_tor, I2a_tor.
    If f2 is provided, uses f = f0 + f2*cos(2chi)/R0^2 and averages over chi.
    If f2 is None, uses f = f0 (no chi-dependence).

    Returns dict with keys: I4, I2a, I4_tor, I2a_tor
    """
    N = len(rho_grid)

    # Gauss-Legendre quadrature on [0, 2*pi]
    pts, wts = np.polynomial.legendre.leggauss(N_chi)
    chi_pts = np.pi * (pts + 1)
    chi_wts = wts * np.pi

    I4_int = np.zeros(N)
    I2a_int = np.zeros(N)
    I4t_int = np.zeros(N)
    I2at_int = np.zeros(N)

    for i in range(N):
        r = rho_grid[i]
        fp0 = f0_prime(r, Cs)

        if f2 is not None:
            # Average over chi
            g4 = g2a = g4t = g2at = 0.0
            for j in range(N_chi):
                eps = f2[i] * np.cos(2 * chi_pts[j]) / R0**2
                epsp = f2p[i] * np.cos(2 * chi_pts[j]) / R0**2

                f_val = f0_val(r, Cs) + eps
                fp_val = fp0 + epsp
                sf = np.sin(f_val)

                s4fp2 = sf**4 * fp_val**2
                s2r2 = sf**2 / r**2

                g4 += chi_wts[j] * s4fp2 / r
                g2a += chi_wts[j] * sf**4 * (fp_val**2 + s2r2) * r

                if r < R0:
                    w = r / np.sqrt(R0**2 - r**2)
                    g4t += chi_wts[j] * s4fp2 * w
                    g2at += chi_wts[j] * sf**6 * w

            I4_int[i] = g4 / (2*np.pi)
            I2a_int[i] = g2a / (2*np.pi)
            I4t_int[i] = g4t / (2*np.pi)
            I2at_int[i] = g2at / (2*np.pi)
        else:
            sf = sin_f0(r, Cs)
            s4 = sf**4
            fp2 = fp0**2

            I4_int[i] = s4 * fp2 / r
            I2a_int[i] = s4 * (fp2 + sf**2 / r**2) * r

            if r < R0:
                w = r / np.sqrt(R0**2 - r**2)
                I4t_int[i] = s4 * fp2 * w
                I2at_int[i] = sf**6 * w

    from numpy import trapezoid as trapz
    mask = rho_grid < R0 * (1 - 1e-6)

    return {
        'I4':     trapz(I4_int, rho_grid),
        'I2a':    trapz(I2a_int, rho_grid),
        'I4_tor': trapz(I4t_int[mask], rho_grid[mask]),
        'I2a_tor':trapz(I2at_int[mask], rho_grid[mask]),
    }


def assemble_3D(integrals, R0):
    """Compute J4^{3D}, J2a^{3D}, and ratio from Paper XIV cross-coupled formula."""
    J4 = R0 * integrals['I4'] + QT2 * integrals['I2a_tor']
    J2a = R0 * integrals['I2a'] + QT2 * integrals['I4_tor']
    return J4, J2a, J4 / J2a


def F_BPS(C):
    """Modified BPS equation: F(C*) gives the WZW ratio J4/J2a."""
    return C * (2*C + 1) / ((C**2 + 1) * (3*C - 1))


# ═══════════════════════════════════════════════════════════════════
# Main computation
# ═══════════════════════════════════════════════════════════════════
def main(R0=3.0, Cs=3.4318, N=8000, N_chi=64):
    print(f"Parameters: R0={R0}, C*={Cs}, N={N}, N_chi={N_chi}")
    print(f"Constants:  phi={PHI:.10f}, Qt^2={QT2:.3f}, r_WZW={R_WZW:.10f}")

    # Grid
    rho = np.exp(np.linspace(np.log(1e-4), np.log(50), N))

    # Step 1: Solve BVP
    f2, f2p = solve_bvp(rho, Cs, R0)
    peak_idx = np.argmax(np.abs(f2))
    print(f"\nBVP: peak |f2| = {np.abs(f2[peak_idx]):.6f} at rho = {rho[peak_idx]:.4f}")

    # Step 2: Base integrals (f0 only)
    base = compute_integrals(rho, Cs, R0, f2=None, N_chi=N_chi)
    J4_0, J2a_0, r3D_0 = assemble_3D(base, R0)

    # Step 3: Corrected integrals (f0 + f2)
    corr = compute_integrals(rho, Cs, R0, f2=f2, f2p=f2p, N_chi=N_chi)
    J4_c, J2a_c, r3D_c = assemble_3D(corr, R0)

    # Step 4: Yukawa chain
    C_thin = 3.392  # Paper IV thin-torus saddle
    k = 3
    C_3D_0 = brentq(lambda C: F_BPS(C) - r3D_0, 2, 5)
    C_3D_c = brentq(lambda C: F_BPS(C) - r3D_c, 2, 5)
    dy_0 = np.exp((C_3D_0 - C_thin) / (2*k)) - 1
    dy_c = np.exp((C_3D_c - C_thin) / (2*k)) - 1
    vEW_raw = 0.69  # Paper IV thin-torus residual
    vEW_0 = vEW_raw - dy_0 * 100
    vEW_c = vEW_raw - dy_c * 100

    # Report
    gap_0 = (r3D_0 - R_WZW) / R_WZW * 100
    gap_c = (r3D_c - R_WZW) / R_WZW * 100

    print(f"\n{'='*65}")
    print(f"RESULTS")
    print(f"{'='*65}")

    print(f"\n  Integral shifts (f0+f2 vs f0):")
    for key in ['I4', 'I2a', 'I4_tor', 'I2a_tor']:
        shift = (corr[key] - base[key]) / base[key] * 100
        print(f"    {key:>8s}: {base[key]:.10f} -> {corr[key]:.10f}  ({shift:+.6f}%)")

    print(f"\n  3D ratio:")
    print(f"    r_WZW          = {R_WZW:.10f}")
    print(f"    r3D (f0)       = {r3D_0:.10f}  gap: {gap_0:+.6f}%")
    print(f"    r3D (f0+f2)    = {r3D_c:.10f}  gap: {gap_c:+.6f}%")

    print(f"\n  Yukawa correction:")
    print(f"    C*_3D (f0)   = {C_3D_0:.6f}  dy/y = {dy_0*100:+.4f}%  v_EW = {vEW_0:+.4f}%")
    print(f"    C*_3D (f0+f2)= {C_3D_c:.6f}  dy/y = {dy_c*100:+.4f}%  v_EW = {vEW_c:+.4f}%")

    print(f"\n  Scaling: gap * R0^4 = {abs(gap_c)*R0**4:.4f}% (O(1) check)")

    return {
        'r3D_base': r3D_0, 'r3D_corr': r3D_c,
        'gap_base': gap_0, 'gap_corr': gap_c,
        'vEW_base': vEW_0, 'vEW_corr': vEW_c,
        'f2_peak': np.abs(f2[peak_idx]),
        'f2_rho_peak': rho[peak_idx],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Second-order toroidal correction')
    parser.add_argument('--R0', type=float, default=3.0)
    parser.add_argument('--Cstar', type=float, default=3.4318)
    parser.add_argument('--N', type=int, default=8000)
    parser.add_argument('--Nchi', type=int, default=64)
    args = parser.parse_args()
    main(R0=args.R0, Cs=args.Cstar, N=args.N, N_chi=args.Nchi)
