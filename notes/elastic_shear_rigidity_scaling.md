# Elastic shear rigidity of the director-strain DM — cheap scaling pass

STATUS: interpretive / scaling. Tags each result as robust vs cutoff-set.
Scripts: dm_polarization/{elastic_shear_scaling,crossover_match,crossover_phi6_check,
energy_cutoff_pin}.py; inputs from Paper 7 §dark_matter + {director_speed,cutoff_embedding}.py.
PAPER STATUS: the QUALITATIVE heavy-elastic-solid prediction (coherence + substructure
suppression + intrinsic warps, signs robust) is now Remark P7:rem:elastic_solid in
Paper 7. Everything QUANTITATIVE below (K value, coherence margins, a0/phi^6/cutoff)
stays scratch-only — cutoff-set and, for the phi^6 pinning, a negative result.

## The question
Standard collisionless CDM is a gas of particles on independent orbits — no shear
rigidity. This framework's DM is a FRANK-ELASTIC director strain (stiffness K), so
it can support SHEAR like a membrane/solid. Does that rigidity matter at galactic
scales, and where do the "edge ripples" (splay-unstable boundary modes) sit?

## STANDALONE SUMMARY — the director strain as a heavy elastic solid
One structural fact drives everything: the director strain is an ELASTIC medium with
genuine shear rigidity (Frank K), not a collisionless gas. It is "heavy" — large static
rigidity (mu_static ~ K/r^2 = rho_DM c^2, the gravitating mass density itself), but slow
shear waves (c_dir -> 0) because the rotational inertia chi ~ Lam^3.5 outgrows K. Three
consequences, robust part first:
  (a) COHERENCE (CONDITIONAL on a small cutoff): elastic waves cross the galaxy faster
      than it rotates (c_dir > v/(2pi)) ONLY if the relevant cutoff is small (Lam~11-120,
      the a0 scale). At the large gravitational UV cutoff Lam_UV=M_Pl/phi^6 (~1e30 band
      units) c_dir->0 far below v/(2pi), the elastic dynamics vanish, and the strain
      reduces to STATIC cold DM (dynamically like collisionless CDM). Which cutoff applies
      is unresolved (Result 2; cutoff_derivation_scope.md).
  (b) SUBSTRUCTURE SUPPRESSION (qualitative): a shear-rigid medium resists fragmentation,
      tending to fewer subhalos than collisionless CDM — the direction the missing-
      satellite / too-big-to-fail observations point. (No subhalo mass function computed;
      tendency only. NOTE: the central profile is rho ~ 1/r^2 (SIS-like, from M_enc ∝ r),
      which is CUSPY — so this does NOT claim to solve cusp-core; only substructure.)
  (c) INTRINSIC DISK MODES (interpretive, magnitude cutoff-set): the splay instability
      (K_splay<0) at the disclination edge makes warps / corrugations / rings elastic
      boundary modes of the medium rather than tidal accidents. Wavelength = cutoff-set.
Robust = the SIGNS (shear rigidity; c_dir > v_gal; splay<0). Cutoff-set = every magnitude
(K value, coherence scale, warp wavelength). This is exactly what P7:rem:elastic_solid
states, at the qualitative level; the numbers below are the scratch backing.

## Result 1 — K is fixed by the rotation curve (not free)
Disclination embedding n_hat=e_r: |grad n|^2=2/r^2, Frank density u=K/r^2, and the
DM mass-energy IS that elastic energy: rho_DM c^2 = u = K/r^2. Then
  M_enc(r) = int rho_DM 4pi r'^2 dr' = 4pi K r / c^2   (∝ r, flat curves)
  v_flat^2 = G M_enc/r = 4piGK/c^2   =>   K = c^2 v_flat^2 /(4piG).
So the Frank constant is PINNED by v_flat. For v=200 km/s, K≈4.3e36 N. (Exact
within the disclination model; it is the same M_enc∝r already in P7:sec:flat_curves,
read backwards to isolate K.)

## Result 2 — coherence criterion (holds ONLY at small cutoff; FAILS at large UV cutoff)
Elastic-wave crossing time vs rotation period:
  tau_el/tau_rot = (R/c_dir)/(2piR/v) = v/(2pi c_dir).
Coherent (rigid-like) medium <=> c_dir > v/(2pi) ~ 1e-4 c (for v=200 km/s).
Using c_dir^2 ~ Lam^-2.09 (director_speed.py) anchored at c_dir(120)≈3e-3 c:
  Lam=11 (a0-preferred): c_dir≈3.6e-2 c, c_dir/v_gal≈55, tau_el/tau_rot≈3e-3  COHERENT
  Lam=120              : c_dir≈3.0e-3 c, c_dir/v_gal≈4.5, tau_el/tau_rot≈3.5e-2 COHERENT
CRITICAL CORRECTION (very large UV cutoff): c_dir=v/(2pi) at Lam~3e3; per-rotation
coherence FAILS above that, and even Hubble-time substructure suppression at ~kpc
(needs c_dir > kpc*H_0/c ~ 2.5e-7 c) fails above Lam~1e6. The cutoff that governs c_dir
is the DIRECTOR BAND cutoff (condensate coherence scale / amplitude-mode gap), which is
<= the gravitational Lam_UV=M_Pl/phi^6 (~1e30 band units) but NOT equal to it -- M_Pl/phi^6
is only the EFT ceiling. IF the director cutoff is large (near that ceiling), c_dir->0
(results.md: "at the larger physical cutoff comfortably cold") and ALL elastic dynamics
vanish: the strain is STATIC cold pressureless DM, indistinguishable from collisionless
CDM. Coherence holds ONLY if the director cutoff is small (near the a0 scale Lam~11-120).
WHERE the director band cutoff actually sits is the unresolved question (WP3 of
cutoff_derivation_scope.md); the natural |k|=1/xi_cond reading could put it low enough to
keep the elastic dynamics. The earlier "robust to the cutoff" claim assumed the small-cutoff
regime and is WITHDRAWN until the director cutoff is derived.

## Result 3 — heavy elastic SOLID (the distinguishing prediction)
Static shear modulus mu_static ~ K/r^2 = rho_DM c^2 (large — it is the gravitating
mass density itself). Shear waves are slow (c_dir small) ONLY because the rotational
inertia chi ~ Lam^3.5 is huge ("heavy director"). Net: a heavy elastic solid, not a
collisionless gas. Physical consequence — it resists shear and stays coherent, so it
SUPPRESSES the small-scale substructure a particle halo develops. This runs opposite
to collisionless CDM's substructure/cusp tensions (missing satellites, too-big-to-fail,
cusp-core): elastic DM => smoother, more coherent, substructure-poor halos.

## Result 4 — edge ripples = splay-unstable boundary modes
Splay branch c_dir^2_splay<0 (imaginary freq = instability), |c_dir^2|_splay ~ Lam^-0.82.
Lives at the disclination boundary (coherent medium meets vacuum); a Frederiks-like
instability with wavelength set by the core/gradient scale (magnitude cutoff-set).
Observable counterpart: disk warps / corrugations / rings as ELASTIC boundary modes
rather than tidal accidents. (Real disks show exactly these; here they'd be intrinsic.)

## Result 4b — is c_dir radius-dependent? NO (cdir_radial.py) — tests "edge at higher cutoff"
User hypothesis: the rippling outer edge is at a higher UV cutoff. TEST: Frank elasticity
is HARMONIC, u=(1/2)K|grad n|^2, so c_dir=sqrt(K/chi) is built from BULK band integrals
(K, chi) carrying NO r-dependence. The disclination background n=e_r enters the ENERGY
DENSITY (u~K/r^2) but not the perturbation SPEED. Result (flat curve, v=const):
  c_dir(r)      = CONST   (wave speed uniform)
  t_el/t_dyn    = v/c_dir = CONST  (coherence UNIFORM across the disk -- coherent-or-not
                  together, set by the global CUTOFF, not by r)
  background:   strain ~1/r, u ~1/r^2 (all radial variation is here, not in the speed)
  nonlinear:    (core/r)^2 ~ 1e-6 at kpc -> O(1) only at r~core (microscopic) = a CORE
                effect, not an EDGE effect.
VERDICT: the "outer edge at a higher cutoff" picture is NOT supported at harmonic order --
c_dir is uniform. The edge ripples are a FREE-BOUNDARY effect (the medium terminates at
R_halo; the splay instability lives on that surface), NOT a radial wave-speed/cutoff
gradient. Caveat: the observed warps are at the STELLAR-disk edge while the DM disclination
extends to R_halo, so mapping the splay-boundary mode to observed warps needs an extra
assumption about where the DM boundary sits.

## Result 5 — payoff for the cutoff hunt
- a0 magnitude and ripple wavelength both scale with the same Lam*.
- The coherence verdict is NOT cutoff-independent (correction): it holds only in the
  small-cutoff regime (Lam~11-120) and FAILS at the large gravitational UV cutoff, where
  the medium is static cold DM. So "elastic DM = smooth/coherent/substructure-poor" is a
  prediction of the SMALL-cutoff regime, contingent on the unresolved IR-vs-UV question.
- A measured halo-coherence scale (substructure suppression) or intrinsic warp
  wavelength would then be a discriminator BETWEEN the two regimes (small cutoff =>
  elastic signatures present; large cutoff => none), i.e. an observational handle on
  where the relevant cutoff sits.

## Cutoff status (what the papers already fix vs what's open)
Fixed/derived: Lam_UV = M_Pl/(4pi sqrt2) ≈ M_Pl/phi^6 ≈ 0.0557 M_Pl (Seeley-DeWitt,
P7:cor:luv_tower); Lam_cond = T_CMB(pi^2/150)^{1/4}; m_xi ~ H_0 (DE sector).
Open = ONE O(1) ratio: m_xi / (semi-Dirac crossover E*), claimed ≈ phi^6. The a0
work (cutoff_embedding.py) reduced it to "does m_xi sit at Lam*≈11-12 band units?"
(vs phi^6≈17.9 — a factor ~1.5 gap, "order unity" per the paper).

## Crossover-matching RESULT (crossover_match.py + crossover_phi6_check.py, RUN)
Semi-Dirac dispersion E=sqrt((kx^2/2)^2 + ky^2) (band units): quadratic heavy axis
(kx^2/2), linear light axis (ky). Crossover k*=2, E*=sqrt8=2.828.
Regularised band_K(Lam)=A Lam^a, A=5.226e-3, a=1.357 (smooth cutoff). a0=band_K*cH0.
a0=cH_0/2pi requires band cutoff Lam*=12.4 => m_xi/E* in [12.4 (light axis), 38.4
(heavy axis)] -- the physical ratio is bracketed by the dispersion anisotropy.

TEST of the S_eff claim m_xi/E* = phi^6: propagate phi^n through both axes ->
  phi^5 (11.1): a0 in [0.47,0.93]e-10  -> obs 1.2e-10 ABOVE bracket   DISFAVOURED
  phi^6 (17.9): a0 in [0.65,1.79]e-10  -> obs INSIDE, geomean 1.08e-10 (0.90x obs)  BEST
  phi^7 (29.0): a0 in [0.90,3.44]e-10  -> obs inside but near low edge (1.46x)  allowed
So the a0 data FAVOUR phi^6: obs sits at the geometric centre of the phi^6 light/heavy
bracket (within 10%), phi^5 is excluded (a0 above), phi^7 off-centre. phi^6 log-position
in the a0-required [12.4,38.4] bracket = 0.33.

STATUS: this CONFIRMS the independently-motivated S_eff ratio phi^6 reproduces a0
(within the dispersion anisotropy) and singles it out over phi^5; it does NOT uniquely
DERIVE phi^6 (bracket ~2.7x wide still admits phi^7). Model-dependent inputs: the
band_K power law (regularised integral), H_0, and the O(1) 3D-embedding geom factor
folded into band_K. Upgrade path to a unique derivation: fix the O(1) embedding factor
and narrow the axis-map ambiguity (proper angular average of the anisotropic cutoff
instead of the light/heavy extremes) -> would collapse the 2.7x bracket.

## PINNING ATTEMPT — negative result (energy_cutoff_pin.py, conv_check.py)
Goal: collapse the 2.7x light/heavy bracket by replacing the CIRCULAR momentum
cutoff exp(-2|k|^2/Lam^2) (ambiguous as an energy for the anisotropic band) with
the ISO-ENERGY cutoff exp(-2(E(k)/m_xi)^2), single-valued in m_xi. Cartesian
integrator validated vs the polar reference (K_bend within 3-20%, best at large scale).

RESULT: the bracket collapses, but the required ratio MOVES AWAY from phi^6, not onto it:
  circular cutoff (geom_coeff=1): m_xi/E* in [12.4, 38.4]  (bracket contains phi^6=17.9)
  ENERGY cutoff  (geom_coeff=1): m_xi/E* ~ 8-10 (grid-sensitive, not converged)
  => under the correct single-valued cutoff, phi^6 OVERSHOOTS a0 by ~3x.
Convergence is poor (m_req/E* = 8.6, 9.3, 10.4 at dk=0.06,0.15,0.10; exponent swings
1.5-1.9) -> the elongated energy-cutoff domain is not numerically settled here.

HONEST VERDICT (supersedes the "phi^6 favoured" line above):
- phi^6 is NOT pinned. Its viability hinges on an UNDER-SPECIFIED modelling choice:
  is the physical condensate cutoff on ENERGY (bandwidth -> energy cutoff -> ratio ~9,
  phi^6 fails) or on MOMENTUM (coherence length -> circular cutoff -> bracket [12,38],
  phi^6 fits)? The framework does not fix this, and geom_coeff (the 3D-embedding O(1)
  factor) is undetermined on top of it.
- The earlier "a0 data favour phi^6" was an ARTIFACT of the circular-cutoff bracket
  happening to straddle phi^6; it does not survive the energy-cutoff prescription.
- ROBUST content is only order-of-magnitude: a0 ~ cH_0 (band_K ~ O(0.1-1), a0 ~ 1e-10
  with no tuning) and required ratio ~ O(10). The exact coefficient (the 2pi, the phi^6)
  is not derivable at this modelling level.
- To actually pin it one must FIRST derive, from the condensate, (i) whether the cutoff
  is on E or |k|, and (ii) the absolute geom_coeff (the dimensionless-band -> physical-K
  normalisation via the condensate energy density). Both are beyond a scaling pass; this
  is the "UV cutoff regularisation" the paper already flags as open. No paper edit.

## Coherence channel: a coarse DISCRIMINATOR, not a fine constraint (honest correction)
At the a0-matched Lam*=12.4, c_dir~0.03c, coherence margin ~48x. Coherence holds for
any Lam up to ~1e3 (per-rotation) or ~1e6 (kpc substructure), then FAILS -- and the
gravitational UV cutoff M_Pl/phi^6 is ~1e30, deep in the failure regime. So the yes/no
coherence verdict does not finely constrain Lam*, but it DOES cleanly separate two
regimes: small cutoff (Lam < ~1e3-1e6) => elastic signatures present; large cutoff
(the full UV scale) => none, static cold DM. An observed presence/absence of the
elastic signatures (substructure suppression, intrinsic warps) therefore tells us which
regime Nature sits in -- a coarse but real handle on where the relevant cutoff lies.
Earlier optimism that coherence is cutoff-robust is withdrawn.
