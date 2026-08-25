"""
Test the stiff-core / floppy-edge picture: does the director-wave speed c_dir
vary with galactic radius r? (User's idea: outer edge at a 'higher cutoff'.)

Frank elasticity is HARMONIC: energy = (1/2) K |grad n|^2. The wave speed
c_dir = sqrt(K/chi) is built from BULK band integrals (K, chi) that carry NO
r-dependence. The disclination background n=e_r enters the ENERGY DENSITY
(u ~ K/r^2) but not the PROPAGATION SPEED of perturbations. So at harmonic
order c_dir(r) = const. Radial structure lives in the background, and in
nonlinear corrections that scale as (core/r)^2 -> matter at the CORE, not the edge.
This script makes that explicit for a flat-curve galaxy (v = const).
"""
import numpy as np
c = 2.99792458e8; kpc = 3.0857e19
v = 200e3                         # flat-curve velocity (const in r)
r = np.array([0.5,1,2,5,10,20,50])*kpc
ell_core = 1e-3*kpc              # disclination core ~ microscopic (illustrative)

# c_dir is a BULK constant (band K, chi). Take the a0-scale value ~0.03c for scale.
c_dir = 0.03*c

print(f"{'r[kpc]':>7} {'|grad n|~1/r':>12} {'u~K/r^2 (rel)':>13} {'shear v/r':>11} "
      f"{'c_dir(r)/c':>11} {'t_el/t_dyn':>11} {'NL (core/r)^2':>13}")
for rr in r:
    strain = 1/rr
    u_rel  = (kpc/rr)**2                    # u ~ 1/r^2 (relative units)
    shear  = v/rr
    cdir_r = c_dir                          # <-- CONSTANT: harmonic Frank
    coh    = v/c_dir                        # t_el/t_dyn = (r/c_dir)/(r/v) = v/c_dir, r-INDEP
    nl     = (ell_core/rr)**2               # nonlinear splay correction scale
    print(f"{rr/kpc:7.1f} {strain*kpc:12.3f} {u_rel:13.3e} {shear:11.3e} "
          f"{cdir_r/c:11.3e} {coh:11.3e} {nl:13.1e}")

print()
print("READ-OFF:")
print("- c_dir(r)/c is CONSTANT -> the director wave speed does NOT vary with radius.")
print("- t_el/t_dyn = v/c_dir is CONSTANT -> coherence is UNIFORM across the disk")
print("  (the whole medium is coherent-or-not together, set by the CUTOFF, not by r).")
print("- Radial variation is entirely in the BACKGROUND (strain ~1/r, u ~1/r^2), not")
print("  in the perturbation speed. Nonlinear corrections ~ (core/r)^2 -> O(1) only at")
print("  r ~ core (microscopic): a CORE effect, not an EDGE effect.")
print()
print("=> The 'outer edge at a higher cutoff' picture is NOT supported: c_dir is uniform.")
print("   Edge ripples are a FREE-BOUNDARY effect (medium terminates at R_halo; splay")
print("   instability lives on that surface), not a radial wave-speed/cutoff gradient.")
