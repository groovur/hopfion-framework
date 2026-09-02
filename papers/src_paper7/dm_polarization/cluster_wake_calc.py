#!/usr/bin/env python3
r"""
Director WAKE of a galaxy moving at v through the condensate -- the cluster-deficit lever.
Careful setup (moving-source director response), which OVERTURNS the optimistic f~0.5 estimate.

Director dynamics: chi d_t^2 n = K grad^2 n  (c_dir^2 = K/chi). Moving source ~ delta(x - v t).
Fourier: response ~ source / [chi(omega^2 - c_dir^2 k^2)], on-shell omega = k.v.
Denominator vanishes when v cos(theta_k) = c_dir  -> a MACH CONE at cos(theta_k)=c_dir/v.
For v>c_dir (supersonic) the source radiates director waves into the cone (Cherenkov) = the wake.

TWO candidate deposition channels -- BOTH fail to ENHANCE beyond the RAR value:

(A) ENERGY / radiated wake.
  Cherenkov drag force F_drag ~ (radiated power)/v. The radiated power scales with the medium's
  ability to CARRY AWAY energy, i.e. with c_dir. As c_dir -> 0 the medium is FROZEN: no propagating
  modes, no radiation, NO dissipation -> F_drag -> 0. A frozen medium exerts (elastic) restoring, not
  drag. => radiated-wake deposition -> 0 as c_dir->0. Also the mass is energy/c^2 ~ (v/c)^2-suppressed
  (v/c~1e-3) -> ~1e-6 M_halo even if it didn't vanish. Negligible either way.

(B) CONFIGURATION / halo left behind.
  c_dir->0 -> the heavy director cannot follow the galaxy -> the galaxy leaves its (frozen) strain
  halo behind = it becomes intracluster strain. GOOD -- this deposits ~M_halo. BUT re-forming a new
  halo takes t_form ~ R_halo/c_dir -> infinity as c_dir->0. So each galaxy deposits its ONE halo ONCE
  and then goes halo-less; it CANNOT re-source to deposit again. => total deposited ~ N_gal * M_halo
  = the RAR value (~5x stellar). NO enhancement over what the halos already provide. Redistribution,
  not amplification.

=> NEITHER channel enhances the cluster mass beyond the RAR/MOND-short value. The colder the director
   (c_dir->0, the very feature that makes DM 'cold'), the MORE it forbids the wake enhancement:
   frozen medium = no dissipation (A) and infinite re-formation time (B).
"""
import math

# scaling check for (A): Cherenkov/relaxational drag vanishes with c_dir
# F_drag ~ m_eff * v * (relaxation rate) ; relaxation rate ~ c_dir / L  -> F_drag ~ m_eff v c_dir / L
def drag_scaling(m_eff, v, c_dir, L):
    return m_eff * v * c_dir / L

# (B): total deposited = sum of the galaxies' one-time halo shedding = RAR value
def deposited_configuration(N_gal, M_halo):
    return N_gal * M_halo   # each galaxy deposits its single frozen halo ONCE (t_form->inf blocks re-forming)

N_gal, M_halo = 200, 5.0        # M_halo in stellar-mass units (~5x stars)
Mstar_total = 200.0
DM_needed = 28.0 * Mstar_total

print("=== Director wake of a moving galaxy -- careful calc ===\n")
print("(A) radiated-wake drag ~ m_eff v c_dir / L  ->  scales WITH c_dir:")
for cdir in (1.0, 0.1, 0.01, 0.0):
    print(f"    c_dir={cdir:>4}: F_drag/(m_eff v/L) = {cdir:.3f}   -> c_dir->0 gives F_drag->0 (frozen, no radiation)")
print("    + radiated mass = energy/c^2 ~ (v/c)^2 ~ 1e-6 -> negligible regardless.\n")

dep = deposited_configuration(N_gal, M_halo)
print("(B) configuration shedding (galaxy leaves its frozen halo behind):")
print(f"    total deposited = N_gal * M_halo = {dep:.0f}  (stellar units) = {dep/Mstar_total:.0f}x stellar")
print(f"    DM needed = {DM_needed:.0f} = {DM_needed/Mstar_total:.0f}x stellar  -> STILL {DM_needed/dep:.1f}x SHORT")
print("    each galaxy deposits its ONE halo ONCE; t_form~R_halo/c_dir->inf blocks re-forming -> no accumulation.\n")

print("=== VERDICT ===")
print("The wake lever does NOT close the cluster deficit. c_dir->0 (the coldness that makes DM 'cold') is exactly")
print("what forbids the enhancement: frozen medium -> NO radiated drag (A), and INFINITE halo re-formation time")
print("-> single deposition = the RAR value, no amplification (B). The optimistic f~0.5 estimate assumed repeated")
print("shedding-and-resourcing, which c_dir->0 forbids. So the cluster deficit stands: the framework shares the")
print("MOND cluster problem, and the wake/drag lever -- carefully analysed -- does not resolve it.")
print("HONEST: uncertainties remain (exact director dynamics, halo formation history), but the direction is clear")
print("-- colder director = LESS wake enhancement, opposite to the optimistic hope.")
