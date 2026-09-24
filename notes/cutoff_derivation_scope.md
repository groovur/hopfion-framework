# Scoping the cutoff derivation — closing a0 and the elastic-solid magnitudes

GOAL: derive a_0 (and the elastic-solid coherence/warp magnitudes) from the
condensate with NO free coefficient, by fixing the two blockers surfaced by the
pinning attempt: (U1) the UV cutoff's nature (E vs |k|) and value, and (U2) the
geom_coeff that converts the dimensionless band stiffness K_band -> physical K_phys.
STATUS: scope only. Includes the WP1 gate analysis (done now, cheap). No paper edit.

## RAISED STAKES (from the large-UV-cutoff + c_dir(r) findings, this session)
Two later results make U1 (nature AND VALUE of the cutoff) more load-bearing than "the
O(1) factor in a0":
- The cutoff VALUE decides whether P7:rem:elastic_solid's dynamics exist at all. c_dir^2
  ~ Lam^-2.09 -> 0, so at the small (a0-scale) cutoff the medium is a dynamically coherent
  elastic solid (substructure suppression, intrinsic warps), while at the large gravitational
  UV cutoff Lam_UV=M_Pl/phi^6 (~1e30 band units) c_dir->0 and it is STATIC cold DM (no elastic
  signatures). So deriving the cutoff value is not just a0's coefficient -- it decides the
  remark's entire physics. Coherence FAILS above Lam~1e3 (per-rotation) / ~1e6 (kpc substructure).
- The cutoff is a SINGLE GLOBAL scale, not a field: c_dir(r)=const (cdir_radial.py; harmonic
  Frank -> bulk K, chi, no r-dependence). So there is exactly ONE number to pin, and it
  sets both a0's O(1) factor and whether the elastic dynamics survive -- one target, two payoffs.
IMPORTANT distinction the WP must respect: the cutoff governing c_dir is the DIRECTOR BAND
cutoff (where the semi-Dirac director EFT breaks -- set by the condensate's own coherence
length / amplitude-mode gap), which is <= the GRAVITATIONAL cutoff Lam_UV=M_Pl/phi^6 but is
NOT equal to it. M_Pl/phi^6 is only the absolute EFT ceiling; the physical director cutoff
is the condensate scale and is likely far lower. So "c_dir cut off at M_Pl/phi^6 -> static
cold DM" is the UPPER-BOUND worst case, not the expected value. The real question WP3 must
answer: where does the director band cutoff actually sit -- near the a0 scale (Lam~11-12,
elastic dynamics survive) or much higher (static cold DM)? The E-vs-|k| choice feeds this:
|k|=1/xi_cond (coherence length) is the natural moderate candidate. WP3 now decides the
remark's physics, not just a0's coefficient.

## What is already derived (band units, dimensionless)
- Semi-Dirac dispersion from S_eff (dispersion_from_Seff.py):
    E^2 = A^2 k^4 sin^4(theta) + v^2 k^2 cos^2(theta),  band units A=1/2, v=1.
    S_eff = sin^4(theta) IS the heavy-band angular factor (intrinsic; no ext field).
- Crossover k*=2, E*=sqrt8=2.83.  Frank K(bend/splay), chi, c_dir^2=K/chi: validated
  band integrals with a RADIAL DISK cutoff |k|<Lam (frank_run.py), Goldstone-gated.
- Physical anchors: v(light) -> c ; Lam_cond = T_CMB(pi^2/150)^{1/4} = 1.19e-4 eV ;
  rho_cond ~ 10 Lam_cond^4 ; director mass m_xi ~ H_0 ; Lam_UV = M_Pl/phi^6.

## WP1 (DONE, gate analysis) — the problem REFRAMES
Two independent facts change the target:

(1) a0 ~ cH_0 is an INFRARED coincidence, not a band result.
    a_0 = v_flat^2/r_*, r_* = c/m_xi = c/H_0 (Hubble radius). The scale a_0 ~ cH_0
    comes entirely from the director mass m_xi ~ H_0 (IR), the standard MOND-cosmology
    coincidence. The band (UV) does NOT set the ~1e-10 magnitude; it only supplies the
    O(1) factor (the 1/2pi). So "phi^6" was never needed for the a_0 SCALE — only,
    possibly, for the O(1) factor.

(2) geom_coeff is NOT O(1) — the "band_K = 0.159" match is a convention.
    cutoff_embedding sets a_0/cH_0 = band_K with geom_coeff:=1. But an INDEPENDENT
    determination of the Frank constant from the rotation curve (Result 1 of the
    elastic-solid note) gives
        K_phys^RC = c^2 v_flat^2 /(4 pi G) ~ 4.3e36 N   (v=200 km/s).
    The band gives a DIMENSIONLESS number K_band ~ 0.1-0.8. Converting requires
    K_phys = K_band * [condensate energy density] * [length]^p. Setting geom_coeff=1
    hides that conversion factor. The two determinations of K_phys agree only if that
    factor works out; it is a specific (and likely large, structured) combination of
    condensate scales, NOT a free O(1). The pinning attempt failed precisely because
    it varied the cutoff shape at fixed geom_coeff=1 -- an unjustified convention.

REFRAMED UNKNOWNS:
  U1 (cutoff nature): frank_run already ARGUES for |k| (radial disk) on Goldstone
     grounds; an anisotropic-in-k energy cutoff is suspect. Needs: verify the gate
     under an energy cutoff, and derive the physical |k|-cutoff = 1/xi_cond (healing
     length). Likely |k|, value = inverse condensate coherence length.
  U2 (geom_coeff): the real open number. = the condensate-energy-density normalisation
     that turns K_band into K_phys. Test: does it reproduce K_phys^RC ~ 4.3e36 N?

## Work packages (after the WP1 reframe)
WP2 -- Physical (A, v) from S_eff.  v=c (given). Extract A in condensate units from
   the phi^6-normalised S_eff -> crossover length ell*=A/v and physical E*, and the
   healing length xi_cond from the condensate profile (Paper I/II). Cost ~1-2h
   (symbolic + one normalisation). GATE: does 1/xi_cond land at a sensible |k|-cutoff?

WP3 -- Cutoff nature (U1).  Check dE(q->0) gate under an energy cutoff vs the disk
   cutoff numerically (reuse frank_run integrand); decide E vs |k| physically from
   what terminates the semi-Dirac band (healing length -> |k|). Cost ~1h (band reuse).

WP4 -- geom_coeff (U2), the crux.  Build the full unit map: K_phys = K_band * rho_cond
   * ell*^p (correct dimensions), using rho_cond ~ 10 Lam_cond^4 and (A,v) from WP2.
   Compare K_phys to K_phys^RC = c^2 v_flat^2/(4 pi G). Cost ~1-2h. GATE: agreement
   within O(1) => geom_coeff derived; a large unexplained gap => scheme-dependent,
   a_0 stays order-of-magnitude (honest stop).

WP5 -- Assemble.  a_0 prediction with propagated error; state whether the O(1) factor
   (1/2pi) and/or phi^6 is derived or remains modelling. Cost ~30min.

## Cost & compute
Total ~1 focused day. Compute is LIGHT: the heavy band integrals are done and cached
(frank_run checkpoint.jsonl); WP3/WP4 reuse them. No long runs. (Per the compute-cost
rule: no multi-hour sweeps needed; the expense is analysis, not CPU.)

## Risks (honest)
R1 (moderate-HIGH): U2 may carry an irreducible scheme-dependent O(1) (or worse) that
   no derivation removes -> a_0 stays order-of-magnitude (a0 ~ cH_0), the 1/2pi and
   phi^6 NOT derived. The WP1 reframe already leans this way. This is the likeliest
   outcome and it is a legitimate result: "a0 ~ cH_0 robust; O(1) factor open."
R2 (moderate): A from the 2D band may not normalise cleanly without the full 3D
   condensate profile (the 2D semi-Dirac is a reduction). WP2 may need Paper I's profile.
R3 (low-moderate): the phi^6-in-S_eff may govern the ANISOTROPY (v vs A ratio), not the
   cutoff ratio -- in which case phi^6 never enters a_0 and the earlier target was wrong.

## Go / no-go
WP1 (done) already reframes: the phi^6-pins-a_0 program is probably ill-posed (a_0 is
IR; the O(1) factor is the only place a UV number could enter, via geom_coeff). RECO:
run WP4 FIRST as the true gate -- compute K_phys from the condensate and compare to
K_phys^RC ~ 4.3e36 N. If they agree within O(1), geom_coeff is derived and the rest
follows; if the gap is large/structured, stop and report "a0 ~ cH_0 robust, O(1) factor
scheme-dependent" as the honest closure. WP2/WP3 are only worth doing if WP4 gates green.

## WP4 -- DONE (wp4_geomcoeff_gate.py). Result: ORDER-OF-MAGNITUDE SUCCESS.
Unit map (nematic, verified vs liquid crystals K~k_B T/a): K_phys ~ K_band * rho_cond * ell^2,
rho_cond = 10 Lam_cond^4 = 4.18e-14 J/m^3 (= CMB energy density, from T_CMB). Target
K_phys^RC = c^2 v_flat^2/(4 pi G) = 4.29e36 N. The decisive input is WHICH length ell:
  (i)  condensate coherence xi_cond = hbar c/Lam_cond ~ 1.7 mm -> K_phys = 1.8e-20 N,
       -56 dex. CATASTROPHIC FAIL.
  (ii) ULTRALIGHT DIRECTOR Compton length ell_dir = c/m_xi ~ c/H_0 ~ 4.3 Gpc
       -> K_phys = 1.2e38 N, K_phys/K_RC = 27 (+1.4 dex). NEAR-MATCH.
The 27x residual is fully absorbable: v=c/phi (/phi^2=2.6), m_xi=2H_0 (/4), both -> 2.6x;
exact match needs m_xi/H_0 = 5.2 (order unity). So:
- The DM Frank constant is DERIVED (scale, not coefficient) from rho_cond (T_CMB) and
  m_xi~H_0 (DE sector), with NO galactic input, right to ~1 order of magnitude.
- This works ONLY because the director is ULTRALIGHT: its elastic response is COSMOLOGICAL
  (ell~c/H_0), which is what supplies the enormous galactic Frank constant. The condensate
  coherence length (mm) fails by 56 orders. So the elastic-DM identification is
  quantitatively viable, and viable *because* m_xi~H_0 -- a real cross-sector consistency
  (it is the a0~cH_0 / rho_cond~rho_CMB coincidence in disguise).
- The O(1) COEFFICIENT (phi^6, 2pi) is NOT pinned: the 27x depends on v (c vs c/phi),
  m_xi/H_0 (~5), and O(1) geometry. R1 CONFIRMED. HONEST CLOSURE: "a0 ~ cH_0 and the DM
  Frank constant are derived to ~1 dex from T_CMB + m_xi~H_0; the exact coefficient is a
  structured unit-map value (v, m_xi/H_0, geometry), not a pinned pure phi-power."
NET across WP1-WP4: the pinning of phi^6 is NOT achieved (and WP2 showed why -- phi^6 is in
the unit map, not the band, and the map carries a factor-phi v-ambiguity + m_xi/H_0). But
WP4 upgrades the a0/elastic-magnitude status from "order-of-magnitude, coefficient free" to
"SCALE DERIVED to ~1 dex from T_CMB and the ultralight director, coefficient open." That is
a stronger and honest closure. No paper edit beyond what P7:rem:elastic_solid already says
(magnitudes cutoff/condensate-set); could add one sentence that the scale is condensate-derived.

## WP3 -- DONE (wp3_v_and_reconcile.py). v resolved; m_xi is the real lever.
WP3 resolves the light-axis velocity: v = c/phi (the attractor k-essence sound speed
c_s=1/phi, P7:rem:cs_de -- the DM medium sits at the attractor rho_inf; v=c was the
low-density O(beta^0) assumption). Effect on the WP4 gate: only factor phi^2=2.6.
RECONCILED WP4 gate (K_phys/K_RC), which the m_xi question forced open:
             m_xi=H0    m_xi=phi^10 H0   m_xi=phi^10 sqrt3 H0
  v=c        +1.4 dex   -2.7 dex         -3.2 dex
  v=c/phi    +1.0 dex   -3.2 dex         -3.6 dex
Exact match wants m_xi = 3.2 H0 (v=c/phi) or 5.2 H0 (v=c) -- i.e. a FEW x H0.
=> The gate essentially DEMANDS m_xi ~ few H0. That matches the paper's OPERATIONAL
"m_xi ~ H0" (used in the Vainshtein sector) but TENSIONS the explicit FORMULA
m_xi = phi^10 sqrt3 H0 (~213 H0, P7 l.1011), which would undershoot by ~3 dex.
FLAG: the paper states m_xi inconsistently -- formula phi^10 sqrt3 H0 (~213 H0) vs
repeated "m_xi ~ H0 ~ 1e-33 eV". WP4/WP3 independently favour m_xi ~ few H0. Worth a
dedicated check of the phi^10 prefactor in P7:eq:mxi (is it phi^10 or is "~H0" right?).

## OVERALL CLOSURE (WP1-WP4 + m_xi)
1. phi^6 is NOT pinned: it lives in the unit map, not the band shape (WP2); the map carries
   a factor-phi (v, WP3-resolved to c/phi) and, dominantly, the m_xi prefactor.
2. The elastic-DM identification is QUALITATIVELY VINDICATED: the DM Frank constant needs a
   COSMOLOGICAL elastic length, which the ULTRALIGHT director (ell~c/m_xi) supplies -- 50+ dex
   better than the condensate coherence length (mm). Powered by rho_cond(T_CMB) + m_xi(~H0).
3. QUANTITATIVELY it lands within ~1 dex IF m_xi ~ few H0 (operational), or ~3 dex short if
   m_xi = phi^10 sqrt3 H0 (formula). The O(1..100) coefficient is NOT pinned. R1 stands.
4. Honest paper-ready statement: "the dark-matter scale is set by the same condensate energy
   density (rho_cond ~ 10 Lam_cond^4, from T_CMB) and the ultralight director mass m_xi ~ H0
   that fix dark energy; it reproduces the galactic Frank constant to ~1 order of magnitude
   with no galactic input, the residual coefficient being tied to the precise m_xi/H0 and O(1)
   geometry." NO stronger claim; and resolve the m_xi phi^10 tension before quoting a number.
