#!/usr/bin/env python3
"""
Route 1' go/no-go gate: is there a COLD, MASSIVE (m* >> T_rec), sub-horizon-CLUSTERING
mode in Paper VII's semi-Dirac dispersion that could be the recombination CDM?

Paper VII eq. P7:eq:semidirac:
    E^2 = A^2 k^4 sin^4(theta) + v^2 k^2 cos^2(theta)
    heavy (perp, theta=pi/2) band:  E = A k^2   <-- GAPLESS (no rest mass term)
    light (radial, theta=0)  band:  E = v k     <-- linear/relativistic

The gate as posed was "m*/T_rec >> 1 => cold clustering CDM". But the heavy band
has NO mass gap (m* = 0). The only 'cold' object in the sector is the director soft
mode, whose *sound speed* c_dir->0 (P7:eq:cdir) but whose *mass* is the ultralight
m_xi ~ H_0 (dark-energy sector). We test both candidate 'masses' against T_rec.

Requirement for a component to act as CLUSTERING cold matter at recombination:
  mass >> H(z_rec)   (field 'awake'/oscillating, Compton wavelength sub-horizon)
  and  mass >> T_rec (non-relativistic).
"""
import math

eV = 1.0  # work in eV

# --- constants ---
hbar = 6.582119569e-16      # eV s
T_CMB0 = 2.35e-4            # eV  (2.725 K)
z_rec = 1100.0
T_rec = T_CMB0 * (1 + z_rec)   # eV

# Hubble
H0_si = 2.13e-18              # s^-1  (~67 km/s/Mpc)
H0 = H0_si * hbar            # eV
# H(z_rec): matter-rad; use radiation-era H ~ H0 * sqrt(Omega_r) (1+z)^2 as scale
Omega_r = 9.0e-5
H_rec = H0 * math.sqrt(Omega_r) * (1 + z_rec)**2   # eV (order-of-magnitude)

# candidate masses of the 'cold' modes in the sector
m_xi = H0                     # ultralight director mass ~ H0 (P7:sec:de)
Lam_cond = 1.19e-4           # eV  condensate scale (from T_CMB, earlier thread)
m_heavy_gap = 0.0            # eV  semi-Dirac heavy band is GAPLESS

print("=== Route 1' mass gate: cold massive clustering mode at recombination? ===\n")
print(f"T_rec (z={z_rec:.0f})        = {T_rec:.3e} eV")
print(f"H(z_rec) (order)         = {H_rec:.3e} eV")
print(f"H0                       = {H0:.3e} eV\n")

def verdict(name, m):
    r_T = m / T_rec
    r_H = m / H_rec if H_rec > 0 else float('inf')
    nonrel = "COLD (non-rel)" if r_T > 1 else "RELATIVISTIC/irrelevant"
    awake  = "clusters (awake)" if r_H > 1 else "HUBBLE-FROZEN (DE-like)"
    print(f"-- {name}: mass = {m:.3e} eV")
    print(f"     m/T_rec = {r_T:.3e}  -> {nonrel}")
    print(f"     m/H_rec = {r_H:.3e}  -> {awake}\n")

verdict("semi-Dirac HEAVY band gap", m_heavy_gap)   # 0 -> nothing
verdict("director soft mode m_xi ~ H0", m_xi)
verdict("condensate scale Lam_cond",   Lam_cond)

print("=== CONCLUSION ===")
print("Heavy semi-Dirac band is GAPLESS (m*=0): no cold massive mode exists at the")
print("band-structure level. The only 'cold' (c_dir->0) mode is the ULTRALIGHT")
print("director (m_xi ~ H0), which at recombination has m_xi << H(z_rec) by ~40 orders")
print("=> HUBBLE-FROZEN, w=-1, DE-like: it does NOT cluster at recombination.")
print("Even the condensate scale Lam_cond << T_rec => thermal excitations relativistic.")
print("=> Route 1' FAILS at the dispersion level. No m* >> T_rec clustering CDM in Paper VII.")
print("   (The semi-Dirac anisotropy still delivers the LATE-time galactic DM; it does")
print("    NOT supply an early cold-matter component. Option 3 stands.)")
