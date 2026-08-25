import numpy as np
phi=(1+5**0.5)/2
beta=0.452
print("== the two channels are DIFFERENT physics ==")
print("CHAMELEON (rem:chameleon): density-dependent COUPLING, NOT a mass.")
print("  omega_BD^eff = omega_BD (1+beta rho); at high density (beta rho>>1) -> phi^12 rho/3,")
print("  BETA-INDEPENDENT and m_xi-INDEPENDENT.")
oBD=phi**12/(3*beta)
brho=1e6   # beta rho_sun ~ 1e6 (solar-system density, paper l.1058)
oBD_eff=oBD*(1+brho)
print(f"  omega_BD={oBD:.0f}; at beta rho_sun~{brho:.0e}: omega_BD^eff={oBD_eff:.1e} >> Cassini 43000")
print(f"  margin ~ {oBD_eff/43000:.0e}. Resolves Cassini at the framework coupling. ROBUST, no m_xi.\n")

print("VAINSHTEIN (thm:cassini): r_V=(G_N M omega_BD/m_xi^2)^(1/3), uses the COSMOLOGICAL m_xi.")
print("  r_V ~ omega_BD/m_xi^2 = phi^12/m_xi^2. The phi^-8/3 corollary needs m_xi^2 = phi^20 (Jordan).")
print("  Physical cosmological m_xi ~ H0 (phi-free, Einstein) -> r_V ~ phi^4, NOT phi^-8/3.\n")

print("== does the chameleon RESCUE the phi^-8/3 corollary? NO ==")
print("  The chameleon is a COUPLING suppression (beta_eff=beta/(1+beta rho)); it carries NO")
print("  phi^20 mass. The only place phi^20 lives is the Jordan cosmological m_xi in r_V, which")
print("  the Einstein frame cancels. So no local chameleon mass supplies phi^20 to r_V.")
print("  => phi^-8/3 'topological invariant' is a Jordan-frame artifact; chameleon does not save it.\n")

print("== NET ==")
print("  Cassini resolution: ROBUST via the CHAMELEON (coupling, m_xi-independent, beta=0.452).")
print("  Vainshtein is 'complementary' (paper l.1083); its phi^-8/3 corollary is frame-fragile.")
print("  Demoting the corollary does NOT threaten Cassini -- the chameleon carries it.")
