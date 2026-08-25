"""
#2: does the Vainshtein r_V (and the 'r_V ~ phi^-8/3 topological invariant' corollary)
survive the Jordan->Einstein frame finding (physical m_xi ~ H0, phi-power-free)?

Paper: r_V^3 propto omega_BD / m_xi^2, with omega_BD = phi^12/(3 beta), m_xi^2 = phi^20 Lam_obs.
   => r_V propto (phi^12/phi^20)^{1/3} = phi^{-8/3}.  ("topological", 8=20-12, 20=2 Q_cond).
The phi^20 is the JORDAN-frame (non-canonical) scalar mass. The PHYSICAL (Einstein/canonical)
mass is phi-power-free ~ H0. Test the two robust questions.
"""
import numpy as np
phi=(1+5**0.5)/2
H0_eV=1.494e-33; AU=1.496e11; Rsun=6.96e8; c=2.998e8
r3=np.sqrt(3)

print("== Q1: is Cassini resolved robustly (independent of which mass)? ==")
# r_V propto m^{-2/3}. Paper: r_V(Jordan, m~213 H0) = 9.3e6 AU, gamma_PPN=(Rsun/r_V)^{1/2}.
rV_J=9.3e6*AU
g_J=np.sqrt(Rsun/rV_J)
# physical mass m_E ~ H0 -> m smaller by 213 -> r_V larger by 213^{2/3}
fac=213**(2/3)
rV_E=rV_J*fac
g_E=np.sqrt(Rsun/rV_E)
print(f"  Jordan mass (~213 H0): r_V={rV_J/AU:.1e} AU, gamma_PPN={g_J:.1e}")
print(f"  physical mass (~H0)  : r_V={rV_E/AU:.1e} AU (x{fac:.0f}), gamma_PPN={g_E:.1e}")
print(f"  Cassini bound 2.3e-5: BOTH pass. r_V >> solar system either way. RESOLUTION ROBUST.\n")

print("== Q2: does the phi^{-8/3} 'topological invariant' scaling survive? ==")
print("  r_V^3 propto omega_BD/m^2 = phi^12/m^2.")
print(f"  Jordan m^2 = phi^20 Lam_obs  -> r_V propto phi^(12-20)/3 = phi^-8/3 = {phi**(-8/3):.3f}")
print(f"  Physical m^2 ~ Lam_obs (phi-free) -> r_V propto phi^12/3 = phi^4 = {phi**4:.3f}")
print("  => the ENTIRE phi-dependence of r_V flips (phi^-8/3 -> phi^+4) when the physical")
print("     (canonical) mass replaces the Jordan mass. The '8 = 2 Q_cond - 12' bookkeeping,")
print("     hence 'r_V fixed by the topological charge Q_cond=10', rests on the phi^20 that")
print("     the strong NMC CANCELS in the canonical (physical) mass. So the corollary's")
print("     golden-ratio scaling is a JORDAN-FRAME ARTIFACT, not a physical invariant.\n")

print("== the caveat that could rescue the corollary ==")
print("  IF the screening is CHAMELEON (mass density-dependent) and the LOCAL (solar,")
print("  high-density) scalar mass -- not the cosmological ~H0 -- enters r_V, the local")
print("  mass could carry the phi^20 and the corollary could stand. The paper cites BOTH")
print("  Vainshtein AND chameleon screening. Which mass enters r_V (cosmological Einstein ~H0,")
print("  or a local high-density mass) is the open screening-mechanism question.")
print("\nHONEST NET: Cassini resolution ROBUST; the phi^-8/3 topological-invariant COROLLARY")
print("is frame-fragile -- it needs the physical (or local-chameleon) mass to be re-derived.")
