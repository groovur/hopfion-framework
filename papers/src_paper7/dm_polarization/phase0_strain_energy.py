#!/usr/bin/env python3
"""
Phase 0.1: baryon-sourced Frank strain energy at recombination vs what CDM needs.

Robust bound: the elastic strain energy density u_strain = (1/2) K |grad n|^2 CANNOT
exceed the condensate energy density rho_cond (you cannot store more elastic energy in
a medium than the medium's own energy). With K ~ Lam_cond^2 and |grad n|_max ~ 1/xi_cond
= Lam_cond, u_strain_max ~ Lam_cond^4 ~ rho_cond. So rho_cond is the CEILING.

Key structural fact (Gate 1): rho_cond ~ 10 Lam_cond^4, Lam_cond = 0.506 T_CMB, so
rho_cond ~ rho_gamma (radiation-like). At recombination rho_b ~ rho_gamma too, so the
strain ceiling ~ rho_b -- but CDM needs Omega_cdm/Omega_b ~ 5.4 x rho_b.
"""
import math

# constants / cosmology
T0 = 2.35e-4          # eV, T_CMB today
z_rec = 1100.0
kfac = (math.pi**2/150)**0.25   # Lam_cond = kfac * T_CMB
Omega_b = 0.0493
Omega_cdm = 0.265
Omega_gamma = 5.38e-5

T_rec = T0*(1+z_rec)                     # eV
Lam_cond_rec = kfac * T_rec              # eV
rho_cond = 10.0 * Lam_cond_rec**4        # eV^4  (P7/Paper XII)
rho_gamma = (math.pi**2/15) * T_rec**4   # eV^4

# baryon energy density at recombination
rho_b_over_gamma_0 = Omega_b/Omega_gamma
rho_b_over_gamma_rec = rho_b_over_gamma_0/(1+z_rec)   # rho_b ~ (1+z)^3, rho_gamma ~ (1+z)^4
rho_b = rho_b_over_gamma_rec * rho_gamma

# required strain to mimic CDM
rho_strain_req = (Omega_cdm/Omega_b) * rho_b
# ceiling on strain
rho_strain_max = rho_cond   # elastic energy <= medium energy

print("=== Phase 0.1: strain energy at recombination (z=%.0f) ===\n" % z_rec)
print(f"T_rec              = {T_rec:.4f} eV")
print(f"Lam_cond(z_rec)    = {Lam_cond_rec:.4f} eV   (= 0.506 T_CMB)")
print(f"T_rec/Lam_cond     = {T_rec/Lam_cond_rec:.2f}   <-- NOTE: ~2, constant in time (Lam_cond tracks T)\n")
print(f"rho_gamma(z_rec)   = {rho_gamma:.3e} eV^4")
print(f"rho_cond(z_rec)    = {rho_cond:.3e} eV^4   (= strain CEILING; ~rho_gamma, radiation-like)")
print(f"rho_b(z_rec)       = {rho_b:.3e} eV^4\n")
print(f"rho_cond/rho_b     = {rho_cond/rho_b:.2f}     (strain ceiling in units of baryon density)")
print(f"REQUIRED (CDM)     = {Omega_cdm/Omega_b:.2f} x rho_b = {rho_strain_req:.3e} eV^4\n")

shortfall = rho_strain_req/rho_strain_max
print(f"ceiling / required = {rho_strain_max/rho_strain_req:.3f}")
print(f"SHORTFALL          = {shortfall:.1f}x  (even the ABSOLUTE ceiling falls this far short)\n")

print("=== VERDICT ===")
print("Even the ABSOLUTE ceiling (ALL condensate energy stored as coherent strain, S=1,")
print("distortion at the correlation scale -- unphysical at recombination where the medium")
print(f"is only marginally ordered) gives rho_strain <= rho_cond ~ {rho_cond/rho_b:.1f} rho_b,")
print(f"but CDM needs {Omega_cdm/Omega_b:.1f} rho_b => short by ~{shortfall:.0f}x at the ceiling, MORE realistically.")
print("STRUCTURAL REASON: strain energy <= rho_cond, and rho_cond is RADIATION-like (Gate 1),")
print("so at recombination rho_cond ~ rho_gamma ~ rho_b. The strain is capped near rho_b; CDM")
print("needs ~5x more. The channel does NOT close the gap. (MOND-frontier tension confirmed.)")
