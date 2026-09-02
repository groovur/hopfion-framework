#!/usr/bin/env python3
r"""
Cluster DM deficit vs polaron WAKE DEPOSITION (framework-specific, MOND-distinct).

Density-feedback (screening) gives the RAR + planets/dwarfs/cored profiles (direction correct).
The cluster TOTAL-MASS deficit (MOND ~2-5x short) is separate. Framework lever: c_dir->0 makes
every galaxy SUPERSONIC (v >> c_dir), so a moving galaxy sheds a strain WAKE it can't drag; since
c_dir->0 the shed strain RELAXES ~never -> it PERSISTS as intracluster strain and ACCUMULATES over
the cluster's history. MOND has no analog. Question: can the accumulated wake reach the deficit?

Rough scaling (order-of-magnitude; the deposition fraction f is the real unknown -> needs the wake calc):
  intracluster strain M_icl ~ N_gal * f * M_halo_per_gal * N_crossings
"""
import math

# --- cluster / galaxy parameters (typical rich cluster) ---
N_gal        = 200            # luminous member galaxies
Mstar_gal    = 1.0            # per-galaxy stellar mass (units: set Mstar_gal=1 -> everything in these units)
Mhalo_per_gal= 5.0 * Mstar_gal  # galactic director-strain halo ~5x stars (rotation-curve calibration)
sigma_cl     = 1000.0         # km/s cluster velocity dispersion
R_cl_Mpc     = 1.0            # cluster radius
t_cluster_Gyr= 10.0           # cluster age

# --- DM deficit to close ---
# cluster budget DM/gas/stars ~ 85/12/3 -> DM ~ 28 x stellar; MOND-like strain gives ~ (v/sigma)^2 of it.
DM_needed_over_stellar = 28.0
Mstar_total = N_gal * Mstar_gal
DM_needed   = DM_needed_over_stellar * Mstar_total

# single-object strain (RAR, universal a0) supplies ~5x stellar total (the MOND-short piece)
M_strain_rar = Mhalo_per_gal * N_gal      # = 5 * stellar_total
deficit      = DM_needed - M_strain_rar   # the ~2-5x shortfall to make up
deficit_factor = DM_needed / M_strain_rar

# --- wake accumulation ---
# crossing time ~ R/sigma ; number of crossings over cluster age
kms_per_Mpc_Gyr = 1.0/ (3.086e19/ (1e3*3.156e16))   # (km/s) per (Mpc/Gyr)  -- convert
t_cross_Gyr = (R_cl_Mpc * 3.086e19) / (sigma_cl) / (3.156e16)   # Mpc->km /(km/s) -> s -> Gyr
N_crossings = t_cluster_Gyr / t_cross_Gyr

print("=== Cluster DM deficit vs wake deposition (order-of-magnitude) ===\n")
print(f"stellar total           = {Mstar_total:.0f}  (units of per-galaxy stellar mass)")
print(f"DM needed (~28x stellar)= {DM_needed:.0f}")
print(f"RAR single-object strain= {M_strain_rar:.0f}  (~5x stellar)")
print(f"=> DEFICIT              = {deficit:.0f}   (deficit factor {deficit_factor:.1f}x)\n")
print(f"crossing time  ~ R/sigma = {t_cross_Gyr:.2f} Gyr")
print(f"N crossings over {t_cluster_Gyr:.0f} Gyr = {N_crossings:.1f}\n")

print("Accumulated intracluster wake strain  M_icl ~ N_gal * f * M_halo * N_crossings :")
for f in (0.01, 0.05, 0.1, 0.3):
    M_icl = N_gal * f * Mhalo_per_gal * N_crossings
    closes = M_icl / deficit
    print(f"  deposition fraction f={f:>4}: M_icl={M_icl:8.0f}  = {closes:5.2f} x deficit "
          f"({'CLOSES' if closes>=1 else 'short'})")

print("\n=== READ ===")
print("Deficit factor ~5.6x (need ~28x stellar; RAR single-object strain gives ~5x).")
print("Wake accumulation over ~10 crossings closes it for a deposition fraction f ~ 0.4-0.5 PER CROSSING")
print("(f=0.3 -> 0.67x deficit). NOT 'few percent' -- SUBSTANTIAL: each galaxy sheds ~half its")
print("(continuously re-sourced) halo per supersonic crossing. Plausible for v>>c_dir (c_dir->0: the halo")
print("cannot follow at all, so the loosely-bound outer halo detaches), but not small.")
print("VERDICT: RIGHT ORDER OF MAGNITUDE (not orders off), unlike screening (wrong sign) -- a live lever,")
print("but it needs a LARGE f. THE UNKNOWN is f from the actual director wake calc (v>>c_dir response).")
print("Promising + framework-specific (MOND has no wake); magnitude is borderline-reachable, not guaranteed.")
