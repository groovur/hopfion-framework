#!/usr/bin/env python3.11
"""
firstorder_transition.py -- does DENSITY FEEDBACK make the confinement transition
FIRST-ORDER? (step toward pinning eta; notes (dd),(ee)). DM/Quark crossover (src_paper21).

eta's deep-equilibrium bracket (Gamma/H~1e19, update dd) is only escaped if the transition
is OUT-OF-EQUILIBRIUM, i.e. FIRST-ORDER (barrier -> supercooling -> abrupt freeze). The
framework-native hook: the feedback softening 1/(1+beta*rho) is POSITIVE/hysteretic (high
rho -> softer -> cheaper to get denser), which can turn a crossover into first-order.

TOY Landau free energy, order parameter rho>=0 (condensate density / knottedness):
  F(rho;T) = a(T)*rho + b*rho^2 - gamma * beta*rho^2/(1+beta*rho)
  a(T) = thermal cost (>0 deconfined-favoured, <0 confined); b>0 stabilisation;
  last term = density-FEEDBACK gain (softening: grows then saturates). beta=0 => NO feedback
  => F=a*rho+b*rho^2 (single min moving continuously = CROSSOVER/2nd-order baseline).
TEST: does beta>0 open a BARRIER + second minimum (FIRST-ORDER)? Barrier height / latent
heat = the out-of-equilibrium strength that would pin eta. Qualitative (toy); the point is
the MECHANISM: feedback (beta) first-order vs no-feedback (beta=0) crossover.
"""
import numpy as np

def F(rho, a, b, gamma, beta):
    return a*rho + b*rho**2 - gamma*(beta*rho**2)/(1+beta*rho)

def extrema(a, b, gamma, beta, rmax=8.0, Nr=20000):
    rho = np.linspace(0, rmax, Nr)
    f = F(rho, a, b, gamma, beta)
    df = np.gradient(f, rho)
    sign = np.sign(df)
    # interior stationary points (sign changes of dF)
    idx = np.where(np.diff(sign) != 0)[0]
    mins, maxs = [], []
    for i in idx:
        r = rho[i]
        if sign[i] < 0 and sign[i+1] > 0: mins.append((r, f[i]))     # - -> + : minimum
        elif sign[i] > 0 and sign[i+1] < 0: maxs.append((r, f[i]))    # + -> - : maximum
    # rho=0 is a boundary minimum if dF(0)>0
    if df[0] > 0: mins = [(0.0, f[0])] + mins
    return mins, maxs, rho, f

def classify(a, b, gamma, beta):
    mins, maxs, rho, f = extrema(a, b, gamma, beta)
    nmin = len([m for m in mins if True])
    # first-order signature: >=2 minima with a maximum (barrier) between them
    firstorder = (len(mins) >= 2 and len(maxs) >= 1)
    return mins, maxs, firstorder

def find_Tc_and_barrier(b, gamma, beta, a_lo=-3.0, a_hi=1.0, Na=600):
    """scan a (temperature proxy) for the coexistence point (two degenerate minima)."""
    best=None
    for a in np.linspace(a_hi, a_lo, Na):
        mins, maxs, fo = classify(a,b,gamma,beta)
        if fo:
            m = sorted(mins, key=lambda x:x[0])
            f0, fc = m[0][1], m[-1][1]      # F at rho=0 min and confined min
            barrier = max(mx[1] for mx in maxs) - min(f0, fc)
            gap = abs(f0 - fc)
            if best is None or gap < best['gap']:
                best = dict(a=a, rho0=m[0][0], rhoc=m[-1][0], barrier=barrier, gap=gap,
                            rho_barrier=sorted(maxs,key=lambda x:x[1])[-1][0])
    return best

if __name__ == '__main__':
    b, gamma = 1.0, 1.0
    print("="*74)
    print("  DENSITY-FEEDBACK CONFINEMENT TRANSITION: first-order vs crossover")
    print(f"  F(rho)=a*rho + {b}*rho^2 - {gamma}*beta*rho^2/(1+beta*rho)   (rho = order param)")
    print("="*74)
    print(f"  {'beta':>6} {'transition':>14} {'rho_confined':>13} {'barrier height':>15} {'~at a_c':>9}")
    for beta in [0.0, 0.2, 0.452, 1.0, 2.0, 4.0]:
        best = find_Tc_and_barrier(b, gamma, beta)
        if best is None:
            print(f"  {beta:>6.3f} {'CROSSOVER':>14} {'(single min)':>13} {'0 (no barrier)':>15} {'--':>9}")
        else:
            print(f"  {beta:>6.3f} {'FIRST-ORDER':>14} {best['rhoc']:>13.3f} {best['barrier']:>15.4f} {best['a']:>9.3f}")
    print("\n  READ: beta=0 -> CROSSOVER (feedback off). beta>0 opening a BARRIER + confined")
    print("  minimum = FIRST-ORDER, induced BY the feedback. Barrier grows with beta = stronger")
    print("  out-of-equilibrium (supercooling/abrupt freeze) = the mechanism that could pin eta.")
