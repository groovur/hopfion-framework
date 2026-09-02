# GATE 1: equation of state (a-scaling) of each fabric component. Does ANY sub-component
# scale as 1/a^3 (matter/CDM) at recombination? Scaling arguments (Phase 0: Lam_cond ~ 1/a).
print("Phase-0 inputs: Lam_cond ~ T_CMB ~ 1/a ; xi_cond = 1/Lam_cond ~ a ; rho_cond ~ Lam_cond^4 ~ 1/a^4")
print()
print("Frank constant K ~ rho_cond * xi_cond^2  (units energy/length):")
print("  K ~ (1/a^4)(a^2) = 1/a^2")
print()
print("== component-by-component a-scaling ==")
rows=[
 ("thermal condensate rho_cond ~ Lam_cond^4", "1/a^4", "RADIATION (= CMB reservoir; it IS the temperature)"),
 ("frozen-attractor vacuum Lam_obs (CC)", "const", "DARK ENERGY (exp-suppressed, frozen)"),
 ("baryonic Hopfions (knots, rest-mass)", "1/a^3", "MATTER -- but only Omega_b~0.05"),
 ("DISORDERED director stress (recomb candidate):", "", ""),
 ("   u ~ K |grad n|^2, |grad n|~1/xi_cond~1/a, K~1/a^2", "1/a^4", "RADIATION -> FAILS as CDM"),
 ("   (= rho_cond: disordered flucts on coherence scale carry ~rho_cond)", "", ""),
 ("GALACTIC disclination strain (attached to baryons)", "1/a^3", "MATTER -- but LATE (orders z~few), ABSENT at recomb"),
]
for name,scale,note in rows:
    if scale: print(f"  {name:52s} {scale:>7}  {note}")
    else:     print(f"  {name}")
print()
print("== the disordered-stress check (the crux) ==")
print("  Disordered director fluctuations live on the coherence scale xi_cond and are RELATIVISTIC")
print("  (director/magnitude modes ~ c or c/phi, ultralight m_xi~H0 => not non-relativistic).")
print("  u_disordered ~ K/xi_cond^2 ~ rho_cond ~ 1/a^4.  RADIATION-like, NOT matter.")
print("  For matter (1/a^3) you need NON-relativistic massive excitations = the Hopfion knots (baryons),")
print("  which are only Omega_b~0.05. No EXTRA cold-matter fabric component exists at recombination.")
print()
print("== GATE 1 VERDICT: RED for route 1 ==")
print("  NO fabric sub-component is matter-like (1/a^3) at recombination:")
print("   - disordered stress = radiation (1/a^4)")
print("   - frozen vacuum = DE (const)")
print("   - the only matter-like strain (galactic disclination) is LATE (z~few), absent at z~1100.")
print("  => option 4 in the 'fabric bulk/stress = CDM' sense (route 1) FAILS.")
print("  => 'no separate CDM' can ONLY survive via ROUTE 2: modified pre-recombination gravity")
print("     (density-feedback + induced Einstein eqs) reproducing the CMB peaks with radiation+DE+baryons")
print("     and NO cold matter -- a radical, whole-program claim. Otherwise OPTION 3 (Omega_DM is an input/gap).")
