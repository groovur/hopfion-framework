#!/usr/bin/env python3.11
"""
firstorder_z3.py -- Z3 Polyakov-loop confinement free energy with the FRAMEWORK's inputs,
first-order check. DM/Quark crossover (src_paper21), notes (kk).

V(L) = (b2/2)|L|^2 - (b3/3)(L^3+L*^3) + (b4/4)|L|^4 - h(L+L*)      (Landau-Ginzburg, Z3)
On the Z3-aligned real axis L=l (real):
  V(l) = (b2/2) l^2 - (2 b3/3) l^3 + (b4/4) l^4 - 2 h l
FRAMEWORK INPUTS:
  * b3 != 0 GUARANTEED: trefoil 3 colours admit the Z3 cubic 3x3x3->1 (colour-singlet baryon)
    = the QCD-center cubic. (DERIVED from 2T/trefoil structure.) Normalise b3=b4=1.
  * explicit Z3 breaking h ~ e^{-F_q/T}, F_q = single-STRAND free energy. The framework has NO
    free light quarks -- a strand exists only confined, at the CONSTITUENT scale ~300 MeV. At
    T_QCD~150 MeV: h ~ e^{-300/150} = e^{-2} ~ 0.14 (WEAK breaking). Real QCD: free light quarks
    F_q~m_q (few MeV) -> h~1 (strong -> crossover, the Columbia plot).
  * b2(T): the thermal mass term, swept through the transition (high T: deconfined l!=0 favoured;
    low T: confined l~0). T_c ~ Lcond(QCD) ~ 76 MeV.
CHECK: for the framework h~0.14, is the transition FIRST-ORDER (two coexisting minima + barrier
at some b2) or CROSSOVER (single min throughout)? Compare h=0 (pure gauge) and h~1 (real QCD).
"""
import numpy as np

def V(l, b2, b3, b4, h):
    return 0.5*b2*l**2 - (2*b3/3.0)*l**3 + 0.25*b4*l**4 - 2*h*l

def minima_maxima(b2, b3, b4, h, lmax=4.0, N=40000):
    l = np.linspace(-0.5, lmax, N)
    f = V(l, b2, b3, b4, h); df = np.gradient(f, l); s = np.sign(df)
    mins, maxs = [], []
    for i in np.where(np.diff(s) != 0)[0]:
        if s[i] < 0 and s[i+1] > 0: mins.append((l[i], f[i]))
        elif s[i] > 0 and s[i+1] < 0: maxs.append((l[i], f[i]))
    return mins, maxs

def transition_type(b3, b4, h, b2lo=-2.0, b2hi=2.0, Nb=800):
    """sweep b2 (temperature); FIRST-ORDER if some b2 has 2 minima with ~degenerate depth + a
    barrier between; else CROSSOVER. Return the coexistence b2 and barrier if first-order."""
    best=None
    for b2 in np.linspace(b2hi, b2lo, Nb):
        mins, maxs = minima_maxima(b2, b3, b4, h)
        if len(mins) >= 2 and len(maxs) >= 1:
            m = sorted(mins, key=lambda x: x[0])
            f_lo, f_hi = m[0][1], m[-1][1]
            gap = abs(f_lo - f_hi)
            barrier = max(mx[1] for mx in maxs) - min(f_lo, f_hi)
            if best is None or gap < best['gap']:
                best = dict(b2=b2, gap=gap, barrier=barrier,
                            l_conf=m[0][0], l_deconf=m[-1][0],
                            dl=abs(m[-1][0]-m[0][0]))
    return best

if __name__ == '__main__':
    b3=b4=1.0
    print("="*76)
    print("  Z3 CONFINEMENT FREE ENERGY -- first-order check with framework inputs")
    print("  V(l)=(b2/2)l^2 - (2/3)l^3 + (1/4)l^4 - 2h*l   (b3=b4=1; h=explicit Z3 breaking)")
    print("="*76)
    print(f"  {'case':>26} {'h':>6} {'transition':>13} {'barrier':>10} {'order-param jump':>17}")
    cases = [("pure gauge (h=0)",0.0),
             ("FRAMEWORK (constituent)",np.exp(-2.0)),
             ("h=0.3",0.3),
             ("real QCD (light q)",1.0)]
    for name,h in cases:
        best = transition_type(b3,b4,h)
        if best and best['barrier']>1e-4 and best['dl']>0.05:
            print(f"  {name:>26} {h:>6.3f} {'FIRST-ORDER':>13} {best['barrier']:>10.4f} {best['dl']:>17.3f}")
        else:
            print(f"  {name:>26} {h:>6.3f} {'crossover':>13} {'--':>10} {'--':>17}")
    # locate the critical endpoint h_c (first-order -> crossover)
    hs = np.linspace(0, 1.0, 101); h_c=None
    for h in hs:
        best=transition_type(b3,b4,h)
        fo = best and best['barrier']>1e-4 and best['dl']>0.05
        if not fo: h_c=h; break
    print(f"\n  critical endpoint: first-order for h < h_c ~ {h_c:.2f}  (crossover above).")
    print(f"  framework h~0.14 {'<' if (h_c and 0.14<h_c) else '>'} h_c  ->  "
          f"{'FIRST-ORDER' if (h_c and 0.14<h_c) else 'crossover'}.")
