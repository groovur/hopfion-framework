"""
Cheap scaling pass: does the Frank-elastic director strain give the DM medium
enough SHEAR rigidity to rotate coherently at galactic scales, and where do the
"edge ripples" (splay-unstable boundary modes) sit?

All physics is order-of-magnitude / scaling. Inputs taken from the paper + the
existing dm_polarization scripts (director_speed.py, cutoff_embedding.py):
  - director-wave speed  c_dir^2 = K/chi,  scaling c_dir^2 ~ Lam^-2.09 (bend)
  - anchor: at Lam=120 (band units, band velocity = c), c_dir ~ 3e-3 c
  - a0-preferred physical cutoff Lam* ~ 11-12 band units (cutoff_embedding.py)
  - splay branch unstable: c_dir^2_splay < 0, |c_dir^2|_splay ~ Lam^-0.82
  - disclination embedding n_hat = e_r : |grad n|^2 = 2/r^2, u = K/r^2 = rho_DM c^2
"""
import numpy as np

# ---- constants (SI) ----
c   = 2.99792458e8      # m/s
G   = 6.674e-11         # m^3 kg^-1 s^-2
kpc = 3.0857e19         # m
km  = 1.0e3

print("="*72)
print("1)  Frank constant K is FIXED by the rotation curve (not free)")
print("="*72)
# disclination: u = K/r^2 = rho_DM c^2  => M_enc(r)=int rho 4pi r'^2 dr' = 4pi K r/c^2
#   => v_flat^2 = G M_enc/r = 4 pi G K / c^2   =>   K = c^2 v_flat^2 /(4 pi G)
for vkms in (150, 200, 250):
    v = vkms*km
    K = c**2 * v**2 / (4*np.pi*G)          # units: N (= J/m, a Frank constant)
    print(f"  v_flat={vkms:3d} km/s -> K = c^2 v^2/(4piG) = {K:.2e} N   "
          f"(v/c = {v/c:.2e})")
print("  => K is pinned by v_flat via the disclination model; M_enc=4piKr/c^2 (flat).")

print()
print("="*72)
print("2)  Coherence criterion:  elastic wave-crossing vs. rotation period")
print("="*72)
# tau_elastic / tau_rot = (R/c_dir)/(2piR/v) = v/(2 pi c_dir).
# Coherent (rigid-like) medium  <=>  c_dir >> v/(2pi)  <=>  elastic waves win.
# c_dir(Lam) from anchor c_dir^2(120)= (3e-3)^2, scaling ~ Lam^-2.09  (band vel = c)
def c_dir_over_c(Lam):
    cdir2_120 = (3.0e-3)**2
    return np.sqrt(cdir2_120 * (Lam/120.0)**(-2.09))

vgal = 200*km
vgal_c = vgal/c
print(f"  galactic v_flat = 200 km/s  = {vgal_c:.2e} c ;  threshold v/(2pi) = {vgal_c/(2*np.pi):.2e} c")
print(f"  {'Lam(band)':>10} {'c_dir/c':>10} {'c_dir/v_gal':>12} {'tau_el/tau_rot':>14}  regime")
for Lam in (11, 12, 30, 60, 120):
    cd = c_dir_over_c(Lam)
    ratio = cd/vgal_c
    tau = vgal_c/(2*np.pi*cd)
    regime = "COHERENT (elastic)" if tau < 1 else "winds up (fluid)"
    print(f"  {Lam:>10} {cd:>10.2e} {ratio:>12.1f} {tau:>14.2e}  {regime}")
print("  Criterion c_dir > v/(2pi) ~ 1e-4 c holds by 30-500x across the whole")
print("  plausible cutoff range -> the coherence verdict is ROBUST to the cutoff,")
print("  unlike a0 (which needs the exact O(1) coefficient).")

print()
print("="*72)
print("3)  Shear rigidity vs. collisionless CDM: static shear modulus")
print("="*72)
# static Frank energy density u = K/r^2 = rho_DM c^2  is the gravitating mass.
# effective STATIC shear modulus for a texture shear over scale ell~r: mu ~ K/r^2 = rho c^2.
# but shear-WAVE speed = c_dir (dynamic), tiny, because rotational inertia chi is huge:
#   c_dir^2 = mu_dyn/chi  with chi ~ Lam^3.5 (heavy director).
# => the medium is a "heavy elastic solid": large static rigidity, slow shear waves.
print("  static shear modulus  mu_static ~ K/r^2 = rho_DM c^2  (= the gravitating")
print("  mass-energy density itself); shear-wave speed c_dir = sqrt(mu/chi) is small")
print("  only because the rotational inertia chi ~ Lam^3.5 is large ('heavy director').")
print("  Net: a heavy elastic SOLID, not a collisionless gas -> it resists shear and")
print("  stays coherent, suppressing the substructure a particle halo would develop.")

print()
print("="*72)
print("4)  Edge ripples: splay-unstable boundary modes")
print("="*72)
# splay branch: c_dir^2_splay < 0 (imaginary freq = instability). growth rate ~ |c_dir|_splay * k.
# fastest growth at the largest k the coherent core allows -> wavelength ~ core/gradient scale.
def cdir2_splay(Lam):  # band units, anchor |cdir2|_splay(10)=8.3e-3, ~Lam^-0.82
    return 8.3e-3 * (Lam/10.0)**(-0.82)
for Lam in (11, 30, 120):
    print(f"  Lam={Lam:>3}: |c_dir^2|_splay = {cdir2_splay(Lam):.2e} (band)  "
          f"-> growth-rate scale ~ |c_dir|_splay/ell")
print("  => a splay (Frederiks-like) instability lives at the disclination BOUNDARY,")
print("     wavelength set by the core/gradient scale (cutoff-set magnitude).")
print("     Observable counterpart: disk warps / corrugations / rings as ELASTIC")
print("     boundary modes rather than tidal features.")

print()
print("="*72)
print("5)  Payoff for the cutoff hunt")
print("="*72)
print("  The a0 magnitude AND the ripple wavelength scale with the same Lam*.")
print("  But the COHERENCE verdict (robust, 30-500x margin) is a cutoff-INDEPENDENT")
print("  qualitative prediction: elastic DM = smooth, coherent, substructure-poor halos.")
print("  A measured coherence/warp scale would give a SECOND handle on Lam*, alongside a0.")
