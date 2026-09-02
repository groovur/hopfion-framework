#!/usr/bin/env python3
"""
Phase 2b: collective strain enhancement for the Bullet Cluster -- analytic bracket.

Physics: the splay strain gives an isothermal halo M_enc(r) = k V^2 r / G (M ∝ r,
the same mechanism as galactic flat rotation curves), with V the velocity scale of
the source. The question is which V the strain 'sees':
  - INDIVIDUAL (linear): each galaxy sources its own halo at v_gal ~ 200 km/s.
    Total = sum of N such halos.
  - COLLECTIVE (isothermal): the N-galaxy system sources ONE halo at the cluster
    dispersion sigma_cl ~ 1000 km/s. Since M_iso = sigma^2 R / G IS the cluster
    virial (=observed) mass, full collectivity matches the cluster DM by construction.
The enhancement HEADROOM between the two is (sigma_cl/v_gal)^2. Compare to the
enhancement NEEDED to turn the linear strain into the observed DM.
Everything here is order-of-magnitude; the true collective efficiency needs the
full nematic solver (Phase 2c). This brackets it.
"""
import math

# --- cluster budget (standard fractions; Bullet is a massive merger) ---
f_DM, f_gas, f_star = 0.85, 0.12, 0.03     # of total lensing mass
strain_per_star = 5.0                       # rotation-curve calibration: strain ~5x stellar mass

# --- velocity scales ---
v_gal   = 200.0    # km/s, L* member galaxy circular velocity
sigma_cl = 1000.0  # km/s, Bullet main-cluster velocity dispersion
N_gal   = 200      # luminous member galaxies (order of magnitude)

print("=== Phase 2b: collective strain enhancement (Bullet Cluster) ===\n")

# LINEAR (individual) strain as a fraction of the cluster DM
strain_linear = strain_per_star * f_star          # fraction of total lensing mass
frac_of_DM    = strain_linear / f_DM
short_linear  = 1.0/frac_of_DM
print("LINEAR (each galaxy sources its own ~5x-stellar halo):")
print(f"  strain = {strain_per_star:.0f} x stars = {strain_linear*100:.0f}% of total; DM = {f_DM*100:.0f}%")
print(f"  strain / DM = {frac_of_DM:.2f}  ->  SHORT by {short_linear:.1f}x\n")

# ENHANCEMENT NEEDED to close the gap
enh_needed = short_linear
print(f"ENHANCEMENT NEEDED (linear -> observed DM): {enh_needed:.1f}x\n")

# ENHANCEMENT HEADROOM available if strain goes fully collective (isothermal at sigma_cl)
headroom = (sigma_cl/v_gal)**2
print("COLLECTIVE HEADROOM (strain sees cluster sigma, not galaxy v):")
print(f"  (sigma_cl/v_gal)^2 = ({sigma_cl:.0f}/{v_gal:.0f})^2 = {headroom:.0f}x  available")
print(f"  full collectivity gives an isothermal M = sigma^2 R/G = the cluster VIRIAL mass (matches by construction)\n")

# required fractional collectivity
f_needed = enh_needed/headroom
print(f"REQUIRED COLLECTIVE EFFICIENCY = needed/headroom = {enh_needed:.1f}/{headroom:.0f} = {f_needed*100:.0f}%\n")

print("=== VERDICT ===")
print(f"- Linear (independent) strain gives ~{frac_of_DM*100:.0f}% of the cluster DM -> ~{short_linear:.1f}x SHORT.")
print(f"  The amplitude shortfall is REAL and confirmed.")
print(f"- BUT the collective lever has AMPLE headroom: ~{headroom:.0f}x available vs ~{enh_needed:.1f}x needed.")
print(f"  Closing the gap needs only ~{f_needed*100:.0f}% collective efficiency -- the strain field")
print(f"  responding partly to the cluster potential (sigma) rather than each galaxy (v).")
print("- So the Bullet amplitude is NOT a hard wall: the mechanism has the CAPACITY to close it.")
print("  Whether it achieves ~20% collectivity is the UNRESOLVED number -> needs the nematic solver (2c).")
print("  (Bracket: 0% collectivity -> 5.7x short [fails]; 100% -> isothermal = matches; truth in between.)")
