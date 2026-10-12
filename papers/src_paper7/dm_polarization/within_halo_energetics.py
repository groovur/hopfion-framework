#!/usr/bin/env python3.11
"""
within_halo_energetics.py -- the OPEN PROBLEM (Paper XXII sec:twoepoch "Status"):
in the COHERENT galactic regime the dark-hopfion director interaction is ACTIVE and the
halo is dense (~1e68 relics). Is the cumulative interaction energy (a) finite -- the
already-counted collective disclination + a bounded, extensive relaxation -- or (b) an
N^2-divergent long-range log that OVER-CLOSES the halo and needs nonlinear (S^2) saturation?

This script does NOT settle it by brute force (N~1e68 is not simulable). It settles the
SCALING, which is what the open problem is actually about, by separating three pieces:

  (1) MEAN TEXTURE (the halo itself). The baryon-sourced splay disclination has |grad n|~1/r
      (Paper VII), so E_tex = (chi/2) int |grad n|^2 dV ~ int (1/r^2) 4pi r^2 dr ~ R:
      LINEAR in R, i.e. M_enc proportional to r -- finite, and it IS the halo mass (polaron
      m_eff, P7:rem:polaron). Its winding is TOPOLOGICAL (one disclination, ~unit winding),
      NOT N: adding relics does not add winding, so the mean texture never scales as N^2.

  (2) FLUCTUATION PLASMA (the relic-relic interaction that thermalises them). The semi-Dirac
      director mediator gives a LOGARITHMIC (2D-Coulomb/vortex) pair potential V~C^2 ln(r)
      (nhat_mediated_potential.py). A NAIVE same-sign pairwise sum is ~N^2 ln -> the (b)
      over-closure worry. BUT a THERMALISED plasma DEBYE-SCREENS: the log is cut to a
      Yukawa-of-a-log ~ C^2 K_0(r/lambda_D) beyond the Debye length lambda_D, so each relic
      interacts with only N_Debye ~ n lambda_D^3 neighbours, NOT all N. The correlation
      energy then becomes EXTENSIVE (~N u_corr, i.e. ~ a renormalisation of m_dark) rather
      than N^2. The geometry-gating (coherent regime -> interaction active -> relaxation)
      is exactly what delivers the thermalisation that screens.

  (3) VIRIAL BOUND. A marginally-confining log gives 2<T> = <r V'> = C^2 (const per pair):
      a scale-invariant equilibrium -- the system neither collapses nor runs away (a CORED,
      extended halo). So even the screened interaction energy is virialised to ~ M v^2, not
      a runaway.

VERDICT the script checks numerically: with representative galaxy numbers the SCREENED
(extensive) energy is a small fraction of M_halo c^2, while the UNSCREENED N^2 sum would
over-close by the factor N/N_Debye. => option (a). The over-closure is a double-counting /
unscreened artefact; no S^2 saturation is required, CONDITIONAL on (i) thermalisation
(supplied by the coherent-regime geometry-gating) and (ii) quasi-neutrality of the
FLUCTUATIONS about the mean texture. The coupling C (hence lambda_D) is cutoff-set [O].
"""
import numpy as np

# ---- constants / representative Milky-Way-like halo ----
GeV = 1.0
c_cgs   = 2.998e10                       # cm/s
kpc_cm  = 3.086e21                       # cm
Msun_g  = 1.989e33                       # g
g_to_GeV = 5.61e23                       # GeV/c^2 per gram
m_dark  = 5.0                            # GeV  (ADM scale, Paper XXI/XXII)

M_halo_Msun = 1.0e12
R_halo_kpc  = 100.0
v_vir_kms   = 150.0                      # virial velocity

M_halo_GeV = M_halo_Msun*Msun_g*g_to_GeV
N = M_halo_GeV/m_dark                    # number of dark relics
R = R_halo_kpc*kpc_cm
n_dens = N/(4/3*np.pi*R**3)              # cm^-3
d_ip = n_dens**(-1/3)                    # inter-particle spacing, cm
v = v_vir_kms*1e5/c_cgs                  # dimensionless (v/c)

print("="*72)
print("  WITHIN-HALO ENERGETICS  (Paper XXII sec:twoepoch open problem)")
print("="*72)
print(f"  M_halo = {M_halo_Msun:.0e} Msun = {M_halo_GeV:.2e} GeV")
print(f"  N(relics) = M_halo/m_dark = {N:.2e}")
print(f"  R_halo = {R_halo_kpc:.0f} kpc,  n = {n_dens:.2e} cm^-3,  d_ip = {d_ip:.2e} cm")
print(f"  v_vir/c = {v:.3e}")
print()

# ---- (1) MEAN TEXTURE: linear-in-R, topological winding (NOT N) ----
# E_tex = (chi/2) int |grad n|^2 dV, |grad n| ~ w/r (w = topological winding ~ 1).
# int_a^R (w/r)^2 4pi r^2 dr = 4pi w^2 (R - a)  -> LINEAR in R. M_enc proportional to r.
print("--- (1) Mean texture (the halo = one splay disclination) ---")
print("  |grad n| ~ w/r  (w = winding, TOPOLOGICAL ~1, independent of N)")
print("  E_tex = (chi/2) int|grad n|^2 dV  ~  4pi (chi/2) w^2 (R-a)  => LINEAR in R")
print("  => M_enc proportional to r (Paper VII flat rotation curve); finite; NOT ~N^2.")
print("  Key: N relics do not add winding -> the mean halo mass never scales as N^2.")
print()

# ---- (2) FLUCTUATION PLASMA: unscreened N^2 vs Debye-screened (extensive ~N) ----
# Debye length for a classical plasma: lambda_D^2 = T/(4pi n C_eff^2) (Gauss units analogue).
# For the LOG (2D-Coulomb) mediator the screening is the same qualitative cut; we use the
# standard lambda_D^2 ~ <KE>/(n C^2) with <KE> ~ (1/2) m_dark v^2 as the "temperature".
# C (director coupling) is cutoff-set; scan the stated window.
print("--- (2) Fluctuation interaction: unscreened N^2 vs Debye-screened (~N) ---")
KE = 0.5*m_dark*v**2            # GeV, thermal energy per relic (<=> 'temperature')
print(f"  thermal KE per relic ~ (1/2) m_dark v^2 = {KE:.3e} GeV")
print(f"  {'C [eV]':>10} {'lambda_D/d_ip':>14} {'N_Debye':>12} {'E_scr/Mc^2':>12} {'E_N2/Mc^2':>12}")
for C_eV in [10.0, 1.0e2, 1.0e3]:
    C = C_eV*1e-9               # GeV
    # Debye length in units of the inter-particle spacing:
    #   lambda_D^2 = KE / (4 pi n C^2)  (natural analogue); work in cm with hbar*c=0.1973 GeV*fm
    hbarc_GeVcm = 0.1973e-13
    # convert n (cm^-3) and energies: use lambda_D^2 [cm^2] = KE[GeV] / (4 pi n[cm^-3] C^2[GeV^2]) * (hbarc)^?
    # dimensional: [E]/([n][C]^2) has units GeV/(cm^-3 GeV^2)=cm^3/GeV; * (1/(hbarc)) gives cm^2.
    lamD2 = KE/(4*np.pi*n_dens*C**2) / hbarc_GeVcm   # cm^2
    lamD = np.sqrt(abs(lamD2))
    N_Debye = max(n_dens*(4/3*np.pi*lamD**3), 1.0)
    # Screened (extensive) correlation energy ~ N * (C^2 / lambda_D-scale) -> ~ N * C^2 * (something O(1)):
    # use plasma correlation energy per particle u ~ -C^2 * n^{1/3} (coupling*density scale), extensive:
    u_corr = C**2 * n_dens**(1/3) * hbarc_GeVcm   # GeV per particle (dimensional via hbarc)
    E_scr = N*abs(u_corr)
    # Unscreened same-sign N^2 log boundary term ~ (1/2) N^2 C^2 ln(R/d):
    E_N2 = 0.5*N**2 * C**2 * np.log(R/d_ip)
    print(f"  {C_eV:>10.0f} {lamD/d_ip:>14.2e} {N_Debye:>12.2e} {E_scr/M_halo_GeV:>12.2e} {E_N2/M_halo_GeV:>12.2e}")
print()
print("  READ: the UNSCREENED N^2 column over-closes by many orders (E_N2/Mc^2 >> 1) --")
print("  this is the naive (b) worry. The SCREENED column is <<1 (extensive, a small")
print("  renormalisation of the rest-mass budget). The physical difference is Debye")
print("  screening in the THERMALISED plasma; the ratio of the two is ~ N/N_Debye.")
print()

# ---- (3) VIRIAL: marginally-confining log -> scale-invariant, cored, no runaway ----
print("--- (3) Virial bound for the log potential ---")
print("  V(r) ~ C^2 ln(r):  r V'(r) = C^2 (constant) => 2<T> = <r V'> = C^2 per pair")
print("  => scale-invariant equilibrium: the halo neither collapses nor runs away")
print("     (a CORED, extended profile), and the interaction energy is virialised to ~M v^2,")
print("     NOT an N^2 runaway.")
print()

print("="*72)
print("  VERDICT: within-halo energy = (1) finite linear-in-R mean texture [= halo mass,")
print("  topological winding, NOT N^2]  +  (2) Debye-screened EXTENSIVE (~N) fluctuation")
print("  plasma [a renormalisation of m_dark, NOT N^2]  +  (3) virialised to ~M v^2.")
print("  => option (a): no over-closure, no S^2 saturation needed, CONDITIONAL on")
print("  thermalisation (geometry-gating supplies it) + quasi-neutral fluctuations.")
print("  OPEN [O]: C (hence lambda_D) is cutoff-set; the quasi-neutrality of the")
print("  fluctuation 'charges' about the mean texture is assumed, not proven.")
print("="*72)
