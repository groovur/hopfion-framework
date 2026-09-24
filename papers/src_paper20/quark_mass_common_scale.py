#!/usr/bin/env python3.11
"""
quark_mass_common_scale.py -- RG-run the quark masses to a COMMON scale and re-test
the E_6/CFT tower prediction m_u/m_c = exp(-2pi) (e6_cft_tower_engine.py, probe C).

WHY: PDG quotes m_u at 2 GeV but m_c at m_c=1.27 GeV (different scales), so the raw
ratio 2.16/1270 mixes scales. The physical, RG-invariant comparison is m_u(mu)/m_c(mu)
at a COMMON mu. Because the QCD mass anomalous dimension gamma_m(alpha_s) is FLAVOUR-
INDEPENDENT, the ratio m_i(mu)/m_j(mu) is the SAME at every common mu (both masses
scale by the identical factor) -- so "common scale" yields ONE well-defined number.

METHOD: numerically integrate the MS-bar 3-loop beta-function (alpha_s) and 2-3-loop
mass anomalous dimension, a = alpha_s/(4pi):
  d a /d ln mu^2 = -a^2 (b0 + b1 a + b2 a^2)
  d ln m/d ln mu^2 = -(g0 a + g1 a^2 + g2 a^3)
anchored at alpha_s(M_Z)=0.1179, with n_f thresholds at m_b, m_c.

VALIDATION (must reproduce published values before the numbers are trusted -- sec.1):
  * alpha_s(2 GeV) ~ 0.30   (PDG)
  * m_c(2 GeV) ~ 1.09 GeV   (running PDG m_c(m_c)=1.27 up to 2 GeV; widely quoted)
  * m_b(m_b)=4.18 -> m_b(2 GeV) ~ 4.87 GeV
Only if these match do we trust the m_u/m_c common-scale number.
"""
import numpy as np
from scipy.integrate import solve_ivp

MZ = 91.1876
ALPHAS_MZ = 0.1179
ZETA3 = 1.2020569

def bcoef(nf):
    b0 = 11 - 2*nf/3
    b1 = 102 - 38*nf/3
    b2 = 2857/2 - 5033*nf/18 + 325*nf**2/54
    return b0, b1, b2

def gcoef(nf):
    g0 = 4.0
    g1 = 202/3 - 20*nf/9
    g2 = 1249 - (2216/27 + 160*ZETA3/3)*nf - 140*nf**2/81
    return g0, g1, g2

def run_alpha(mu_hi, a_hi, mu_lo, nf):
    """Integrate a=alpha_s/(4pi) from mu_hi (a_hi) down to mu_lo at fixed nf. t=ln mu."""
    b0, b1, b2 = bcoef(nf)
    def rhs(t, y):
        a = y[0]
        # d a/d ln mu = 2 * d a/d ln mu^2 = -2 a^2 (b0 + b1 a + b2 a^2)
        return [-2*a**2*(b0 + b1*a + b2*a**2)]
    sol = solve_ivp(rhs, [np.log(mu_hi), np.log(mu_lo)], [a_hi],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    return sol.y[0, -1]

def alpha_s(mu):
    """alpha_s(mu) in MS-bar with thresholds at m_b=4.18, m_c=1.27 (mass=threshold)."""
    mb, mc = 4.18, 1.27
    a = ALPHAS_MZ/(4*np.pi)
    if mu >= mb:
        a = run_alpha(MZ, a, mu, 5); return 4*np.pi*a
    a = run_alpha(MZ, a, mb, 5)          # to m_b, nf=5
    if mu >= mc:
        a = run_alpha(mb, a, mu, 4); return 4*np.pi*a
    a = run_alpha(mb, a, mc, 4)          # to m_c, nf=4
    a = run_alpha(mc, a, mu, 3); return 4*np.pi*a

def run_mass(m_hi, mu_hi, mu_lo, nf):
    """Run a mass from mu_hi to mu_lo at fixed nf (no threshold crossing here)."""
    b0, b1, b2 = bcoef(nf); g0, g1, g2 = gcoef(nf)
    def rhs(t, y):
        a, lnm = y
        da   = -2*a**2*(b0 + b1*a + b2*a**2)
        dlnm = -2*(g0*a + g1*a**2 + g2*a**3)
        return [da, dlnm]
    a_hi = alpha_s(mu_hi)/(4*np.pi)
    sol = solve_ivp(rhs, [np.log(mu_hi), np.log(mu_lo)], [a_hi, np.log(m_hi)],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    return np.exp(sol.y[1, -1])

# ── VALIDATION ────────────────────────────────────────────────────────────────
print("="*74)
print("  VALIDATION (reproduce published values before trusting the ratio)")
print("="*74)
a2 = alpha_s(2.0); a3 = alpha_s(3.0)
print(f"  alpha_s(2 GeV)  = {a2:.4f}   (PDG ~0.30)   {'OK' if 0.28<a2<0.32 else 'CHECK'}")
print(f"  alpha_s(3 GeV)  = {a3:.4f}   (PDG ~0.25)   {'OK' if 0.24<a3<0.26 else 'CHECK'}")
mc_2  = run_mass(1.27, 1.27, 2.0, 4)     # m_c(m_c)=1.27 -> m_c(2 GeV)
mc_3  = run_mass(1.27, 1.27, 3.0, 4)
mb_2  = run_mass(4.18, 4.18, 2.0, 4)     # m_b(m_b)=4.18 -> m_b(2 GeV)
print(f"  m_c(2 GeV)      = {mc_2:.4f} GeV   (~1.09)   {'OK' if 1.06<mc_2<1.12 else 'CHECK'}")
print(f"  m_c(3 GeV)      = {mc_3:.4f} GeV   (~0.99)   {'OK' if 0.96<mc_3<1.02 else 'CHECK'}")
print(f"  m_b(2 GeV)      = {mb_2:.4f} GeV   (~4.87)   {'OK' if 4.75<mb_2<5.0 else 'CHECK'}")

# ── COMMON-SCALE RATIOS ─────────────────────────────────────────────────────────
print("\n"+"="*74)
print("  m_u/m_c AT A COMMON SCALE vs the CFT prediction exp(-2pi)")
print("="*74)
PRED = np.exp(-2*np.pi)
mu_2GeV = 2.16e-3     # m_u(2 GeV) GeV
# m_c already at 2 GeV = mc_2 ; ratio at 2 GeV:
r_2 = mu_2GeV/mc_2
# also at M_Z (ratio is scale-invariant; cross-check by running both to M_Z)
mu_MZ = run_mass(mu_2GeV, 2.0, MZ, 4)  # crude: light u run in nf=4 (ignores b threshold, tiny for a ratio)
mc_MZ = run_mass(mc_2,    2.0, MZ, 4)
r_MZ = mu_MZ/mc_MZ
print(f"  CFT prediction         m_u/m_c = exp(-2pi)      = {PRED:.4e}")
print(f"  raw PDG (mixed scales) 2.16/1270                = {2.16/1270:.4e}  ({(2.16/1270)/PRED:.3f}x pred)")
print(f"  common scale @2 GeV    2.16/{mc_2*1e3:.0f} MeV          = {r_2:.4e}  ({r_2/PRED:.3f}x pred)")
print(f"  common scale @M_Z      (RG-invariant check)     = {r_MZ:.4e}  ({r_MZ/PRED:.3f}x pred)")
print(f"\n  m_u PDG uncertainty ~ +0.49/-0.26 MeV on 2.16 (~+23%/-12%) DOMINATES:")
for mu_val,lab in [(2.16-0.26,'m_u low 1.90'),(2.16,'m_u cent 2.16'),(2.16+0.49,'m_u high 2.65')]:
    r = (mu_val*1e-3)/mc_2
    print(f"    {lab} MeV -> m_u/m_c@2GeV = {r:.4e}  ({r/PRED:.3f}x pred, "
          f"{'pred INSIDE' if 0.85<PRED/r<1.18 else 'pred outside'} band)")
print(f"\n  => is exp(-2pi) within the experimental band for m_u/m_c at a common scale?")
