#!/usr/bin/env python3.11
"""
su3_colorsum_coeffs.py -- SU(3) colour-sum of the semi-Dirac holonomy potential ->
Landau coefficients b2,b3,b4 of the Z3 confinement potential. DM/Quark crossover
(src_paper21), notes (nn)->(oo). Toward the exact first-order strength + eta.

Single-phase holonomy potential (semidirac_thermal.py): f(phi) = -sum c_n cos(n phi),
  c_n = 1/n^s  with  s=3 (SEMI-DIRAC, B3, epoch geometry)  vs  s=4 (relativistic, B4, std QCD).
SU(3) along the real-L direction: eigenvalue phases (phi,0,-phi) (traceless), Polyakov loop
  L = (1 + 2 cos phi)/3  (real; L=1 at phi=0 deconfined, L=0 at phi=2pi/3 confined).
  V(L) = f(phi)+f(0)+f(-phi) = const + 2*sum_n c_n [cos(n phi)-1],  cos phi = (3L-1)/2.
Expand V(L) = (b2/2)L^2 - (2 b3/3)L^3 + (b4/4)L^4 + ... (real-axis form) -> read b2,b3,b4.
COMPARE semi-Dirac (s=3) vs relativistic (s=4): is the cubic b3 ENHANCED relative to b4
(=> MORE strongly first-order) for the semi-Dirac epoch geometry?
First-order strength (cubic Landau, real axis V=(b2/2)l^2-(2b3/3)l^3+(b4/4)l^4): barrier/
two-minima exist when the cubic dominates; a scale-free measure is R = b3^2/(b2 b4).
"""
import numpy as np

def g_s(phi, s, nmax=2000):
    n = np.arange(1, nmax+1)
    return np.sum(np.cos(np.outer(phi,n))/n**s, axis=1)

def V_of_L(L, s):
    # cos phi = (3L-1)/2 ; phi in [0, pi]; V(L) = 2*(g_s(phi)-g_s(0))  (drop overall const & sign)
    x = (3*L - 1)/2.0
    x = np.clip(x, -1.0, 1.0)
    phi = np.arccos(x)
    g0 = g_s(np.array([0.0]), s)[0]
    return 2*(g_s(phi, s) - g0)

def coeffs(s):
    L = np.linspace(0.0, 1.0, 400)
    V = V_of_L(L, s)
    # fit V(L) to a polynomial; read L^2,L^3,L^4 coefficients (V(0)=0, V'(0)=0 expected near L=0)
    p = np.polyfit(L, V, 6)          # highest power first
    c = p[::-1]                       # c[k] = coeff of L^k
    c2, c3, c4 = c[2], c[3], c[4]
    # map to b2,b3,b4 via real-axis form (b2/2)L^2-(2b3/3)L^3+(b4/4)L^4
    b2 = 2*c2; b3 = -1.5*c3; b4 = 4*c4
    return b2, b3, b4

if __name__ == '__main__':
    print("="*72)
    print("  SU(3) COLOUR-SUM Landau coefficients from the holonomy potential")
    print("  semi-Dirac (epoch geometry, s=3/B3)  vs  relativistic (std QCD, s=4/B4)")
    print("="*72)
    print(f"  {'dispersion':>18} {'b2':>9} {'b3':>9} {'b4':>9} {'b3/b4':>8} {'R=b3^2/(b2 b4)':>16}")
    for name, s in [("relativistic (B4)",4), ("SEMI-DIRAC (B3)",3)]:
        b2,b3,b4 = coeffs(s)
        R = b3**2/(b2*b4) if (b2*b4)!=0 else float('nan')
        print(f"  {name:>18} {b2:>9.4f} {b3:>9.4f} {b4:>9.4f} {b3/b4:>8.3f} {R:>16.3f}")
    print("\n  READ: larger b3/b4 (and R) = more cubic-dominated = MORE strongly first-order.")
    print("  If semi-Dirac(B3) > relativistic(B4): the epoch geometry ENHANCES first-order (derived).")
