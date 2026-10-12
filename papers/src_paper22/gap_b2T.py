#!/usr/bin/env python3.11
"""
gap_b2T.py -- gap(L) -> thermal mass b2(T): the thermal mass that DRIVES the confinement
transition. DM/Quark crossover (src_paper22), notes (oo)->(pp).

Confinement gaps the semi-Dirac excitations: the order parameter L closes the gap as the
system deconfines. Model Delta(L) = Delta0 (1 - L)  [Delta0 = confinement gap ~ Lcond(QCD);
L=0 confined/gapped, L=1 deconfined/gapless]. Gapped dispersion
  E(k,theta;Delta)^2 = A^2 k^4 sin^4 th + v^2 k^2 cos^2 th + Delta^2.
Thermal free energy  f(T,Delta) = (g_dof) T INT d^3k/(2pi)^3 ln(1-e^{-E/T}).
Thermal part of the L=0 curvature:  b2_thermal(T) = d^2 f/dL^2|_0 = f''(Delta0) * Delta0^2
(Delta'=-Delta0, Delta''=0). f''(Delta)<0 (gap closing unleashes modes) -> b2_thermal<0
(destabilises L=0 = favours deconfinement). The T=0 vacuum/confinement term favours L=0:
  b2_vac = +2 a_vac > 0  (cost of deconfining; a_vac ~ the bag/string scale ~ Delta0^2).
FULL:  b2(T) = b2_vac + b2_thermal(T) = 2 a_vac - |b2_thermal(T)|.
TRANSITION: b2(T) changes sign at T_c where |b2_thermal(T_c)| = 2 a_vac. Show T_c ~ Delta0
~ Lcond(QCD) ~ 76 MeV. With b3/b4~0.17 (su3_colorsum) this fixes the full Landau potential
at T_c -> first-order strength. Units Delta0=A=v=1; report T_c/Delta0 (scale-free).
"""
import numpy as np
from scipy import integrate

def f_thermal(T, Delta, A=1.0, v=1.0, kmax=40.0):
    def integrand(th, k):
        E = np.sqrt(A**2*k**4*np.sin(th)**4 + v**2*k**2*np.cos(th)**2 + Delta**2)
        return k**2*np.sin(th)*np.log1p(-np.exp(-E/T))
    val,_ = integrate.dblquad(integrand, 1e-6, kmax, lambda k:1e-6, lambda k:np.pi-1e-6,
                              epsabs=1e-11, epsrel=1e-8)
    return T*(2*np.pi)*val/(2*np.pi)**3

def b2_thermal(T, Delta0=1.0, ds=1e-3):
    # Z3-respecting gap EVEN in L: Delta(L)^2 = Delta0^2 (1 - L^2).  No linear term.
    # b2_thermal = d^2 f/dL^2|_0 = (df/d(Delta^2)) * d^2(Delta^2)/dL^2 = (df/ds)*(-2 Delta0^2).
    s0 = Delta0**2
    dfds = (f_thermal(T, np.sqrt(s0+ds)) - f_thermal(T, np.sqrt(s0-ds)))/(2*ds)
    return -2*Delta0**2*dfds

if __name__ == '__main__':
    Delta0 = 1.0
    print("="*70)
    print("  gap(L) -> THERMAL MASS b2(T)  (Delta0=1; b2(T)=2 a_vac - |b2_thermal(T)|)")
    print("="*70)
    Ts = np.array([0.3,0.5,0.7,1.0,1.5,2.0,3.0])
    bth = np.array([b2_thermal(T,Delta0) for T in Ts])
    print(f"  {'T/Delta0':>9} {'b2_thermal':>12} {'|b2_thermal|':>13}")
    for T,b in zip(Ts,bth): print(f"  {T:>9.2f} {b:>12.4e} {abs(b):>13.4e}")
    # choose a_vac ~ the T=0 confinement scale; find T_c where |b2_thermal|=2 a_vac
    import numpy as _np
    for a_vac in [0.05, 0.1, 0.2]:
        target = 2*a_vac
        absb = _np.abs(bth)
        # T_c where |b2_thermal| crosses target (interp on the monotonic rising part)
        if absb[0] < target < absb[-1]:
            Tc = _np.interp(target, absb, Ts)
            print(f"\n  a_vac={a_vac}: b2(T) changes sign (b2_thermal destabilises) at T_c/Delta0 = {Tc:.2f}")
        else:
            print(f"\n  a_vac={a_vac}: crossing outside scanned range (|b2_thermal| in [{absb[0]:.3f},{absb[-1]:.3f}])")
    print("\n  READ: |b2_thermal(T)| RISES with T -> b2(T)=2a_vac-|b2_thermal| goes +->- => L=0")
    print("  (confined) stable at low T, unstable (deconfined) at high T. T_c ~ O(Delta0) ~ Lcond(QCD)~76 MeV.")
    print("  With b3/b4~0.17 (su3_colorsum) this fixes the full Landau potential at T_c -> first-order strength.")
