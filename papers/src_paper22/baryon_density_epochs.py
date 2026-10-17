#!/usr/bin/env python3.11
"""
baryon_density_epochs.py -- baryon NUMBER density and chemical potential mu_B at the QCD epoch vs
recombination, to assess whether the cosmological transition can reach the high-mu_B QCD critical point
(RHIC-BES/FAIR/NICA) and whether E_sb/h needs a mu_B (net baryon density) dependence. src_paper22.

QUESTION (user, 2026-10-14): to extend the framework's first-order transition to a (T, mu_B) prediction,
figure out the baryon density at the QCD soup vs recombination.

KEY PHYSICS the numbers expose: in the EARLY UNIVERSE the net baryon asymmetry is fixed, eta = n_B/n_gamma
~ 6e-10, at ALL epochs. So the baryon NUMBER density n_B ~ eta n_gamma ~ a^-3 (dilutes), but the
CHEMICAL POTENTIAL stays mu_B/T ~ eta ~ 1e-10 << 1 at BOTH QCD and recombination. The hot QCD "soup" is
DENSE (thermal q+qbar ~ g_* T^3, huge) but nearly SYMMETRIC (mu_B ~ 0). RHIC/FAIR is the OPPOSITE corner:
cold-ish, HIGH net baryon density, mu_B ~ 200-800 MeV (mu_B/T ~ 1-5). The early universe never goes there.
"""
import numpy as np
from scipy.special import zeta

# constants / conversions
eV_per_MeV = 1e6
hbarc_cm = 1.97327e-5          # eV*cm  (hbar c); 1 eV = 1/hbarc_cm cm^-1
def eV3_to_cm3(n_eV3): return n_eV3 / hbarc_cm**3   # eV^3 -> cm^-3

ZETA3 = zeta(3)                # 1.2020569
eta_B = 6.1e-10               # baryon-to-photon ratio (Planck)

def n_gamma(T_eV):             # photon number density, eV^3
    return (2*ZETA3/np.pi**2) * T_eV**3

def report(name, T_eV, g_star):
    ng = n_gamma(T_eV)
    nB = eta_B * ng
    # thermal (symmetric) particle density ~ (g_star/2) * n_gamma-like (q+qbar+g soup)
    n_thermal = (g_star/2.0) * (2*ZETA3/np.pi**2) * T_eV**3
    # mu_B/T for a near-symmetric relativistic gas: n_net ~ (g_B/6) mu T^2 -> mu/T ~ 6 n_net/(g_B T^3)
    # order-of-magnitude: mu_B/T ~ n_B / n_thermal (net/total)
    muB_over_T = nB / n_thermal
    print(f"  {name:<16} T = {T_eV:.3e} eV")
    print(f"     n_gamma   = {ng:.3e} eV^3 = {eV3_to_cm3(ng):.3e} cm^-3")
    print(f"     n_B (net) = {nB:.3e} eV^3 = {eV3_to_cm3(nB):.3e} cm^-3   (= eta * n_gamma)")
    print(f"     n_thermal = {n_thermal:.3e} eV^3  (q+qbar+g soup, g_*={g_star})")
    print(f"     mu_B/T    ~ {muB_over_T:.2e}      -> mu_B ~ {muB_over_T*T_eV:.2e} eV")
    return nB

if __name__ == '__main__':
    print("="*78)
    print("  BARYON DENSITY & CHEMICAL POTENTIAL: QCD epoch vs recombination")
    print("="*78)
    T_qcd = 150*eV_per_MeV          # 150 MeV
    T_rec = 0.259                   # eV  (z~1100, T~3000 K)
    nB_qcd = report("QCD soup", T_qcd, g_star=51.0)
    print()
    nB_rec = report("recombination", T_rec, g_star=3.36)   # photons + nu (matter era)
    print()
    print(f"  n_B(QCD)/n_B(recomb) = {nB_qcd/nB_rec:.3e}   (= a^-3 dilution = (T_QCD/T_rec)^3 = {(T_qcd/T_rec)**3:.3e})")
    print("="*78)
    print("  RHIC-BES/FAIR/NICA critical-point search:  mu_B ~ 200-800 MeV, T ~ 150-160 MeV")
    print(f"     -> mu_B/T ~ 1.3-5.3   (vs cosmological mu_B/T ~ {eta_B/25:.1e} ~ eta)")
    print("="*78)
    print("  CONCLUSION:")
    print("  * The QCD 'soup' density that breaks strings / sets E_sb is the THERMAL (symmetric)")
    print("    q+qbar+gluon density ~ g_* T^3 -- captured ENTIRELY by T in h = e^{-E_sb/T}.")
    print("  * The NET baryon density (mu_B) is ~1e10x smaller: mu_B/T ~ eta ~ 1e-10 at BOTH QCD and")
    print("    recombination (and all cosmological epochs). The early universe sits on the mu_B ~ 0 AXIS.")
    print("  * RHIC/FAIR scans mu_B/T ~ 1-5 -- a DIFFERENT corner, ~1e10x higher mu_B. Cosmology never")
    print("    reaches it. So E_sb/h needs NO mu_B dependence for the framework's domain (baryogenesis is")
    print("    at mu_B~0); and the high-mu_B critical point is NOT accessible to the cosmological calc.")
    print("="*78)
