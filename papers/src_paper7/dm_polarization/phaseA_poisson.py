# PHASE A (first pass): the framework's modified Poisson equation, and the ONE number that decides route 2.
# Framework gravity = induced GR (G_N) + Brans-Dicke fifth force (density-feedback scalar), chameleon-screened.
#   grad^2 Phi = 4 pi G_eff * delta_rho_source,   G_eff = G_N [1 + 2 beta_eff^2],  beta_eff = beta/(1+beta rho).
# The fifth-force ENHANCEMENT over standard gravity is the 2 beta_eff^2 term. To REPLACE CDM, gravity on the
# baryons must be boosted by ~ (Omega_m/Omega_b) so baryons gravitate like the full matter budget.
import math
phi=(1+5**0.5)/2
beta=0.452
Om, Ob = 0.31, 0.049   # total matter, baryons

print("== the enhancement the framework CAN provide ==")
# background attractor: beta rho_inf = phi (order 1). beta_eff = beta/(1+phi) = beta/phi^2.
beta_eff = beta/phi**2
enh = 2*beta_eff**2
print(f"  beta={beta}, beta rho_inf=phi={phi:.2f} (attractor) => beta_eff=beta/phi^2={beta_eff:.4f}")
print(f"  G_eff/G_N = 1 + 2 beta_eff^2 = 1 + {enh:.4f} = {1+enh:.4f}  (a ~{enh*100:.1f}% enhancement, UNscreened)")
print(f"  If the plasma density screens it (beta rho>>1): beta_eff->0, G_eff->G_N (0% enhancement, standard).")
print(f"  => the framework's fifth-force boost is AT MOST ~{enh*100:.0f}%, and >=0% (screened).")
print()
print("== the enhancement route 2 NEEDS ==")
need = Om/Ob
print(f"  To replace CDM, baryon gravity must mimic the full matter budget: boost ~ Omega_m/Omega_b = {need:.1f}x")
print(f"  i.e. G_eff/G_N ~ {need:.1f} (a ~{(need-1)*100:.0f}% enhancement).")
print()
print("== VERDICT ==")
print(f"  framework can supply: ~{enh*100:.1f}% (best case) to 0% (screened).")
print(f"  route 2 needs:        ~{(need-1)*100:.0f}%.")
print(f"  SHORTFALL: off by ~{(need-1)/enh:.0f}x. The density-feedback fifth force (beta~0.45, /phi^4 suppressed)")
print(f"  is ~2 orders of magnitude too WEAK to replace CDM at the linear level.")
print()
print("== plus the timing check (independent) ==")
print("  The full coherence transition is at z~few (T_ord~Lam_cond<<T_rec), NOT at recombination -> the peaks")
print("  cannot be the ordering signature; at z~1100 the fabric is disordered/radiation-like and only mildly")
print("  (<=6%) enhances gravity, or is screened to standard.")
print()
print("PHASE-A FIRST-PASS CLOSURE: route 2 FAILS at the linear level. The modified-gravity enhancement is")
print("~6% at most (600% needed). => OPTION 3 stands: recombination Omega_DM is a genuine input the framework")
print("does NOT derive. Its modified gravity is far too weak to eliminate CDM. No large Boltzmann run needed.")
