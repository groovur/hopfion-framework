# RESEARCH PROGRAM: does the framework's elastic mechanism, extended to recombination, match the CMB — and where does it differ from ΛCDM?

Capstone program (2026-08-31). Consolidates the session's DM/CMB thread into a
falsifiable research plan. Everything here is a PLAN; nothing enters papers unless
derived + verified. Companion: `fabric_cmb_hypothesis_and_plan.md` (why option 3
stands), `nondeweld_director_scale_program.md` (Part A blocked), base file §9.

## CENTRAL QUESTION
The framework's density-feedback + Frank-elastic gravity reproduces the GALACTIC data
particle-free (flat curves, RAR/a0, no BAO shift) but is a LOW-density/LOW-accel
mechanism that SCREENS to ≈GR at recombination -> currently reduces to GR+baryons ->
underproduces the 3rd CMB peak. Q: extended INTO the high-density recombination regime,
does its own mechanism reproduce the observed C_l (positions + heights + polarization),
and WHERE does it predict something DIFFERENT from ΛCDM (falsifiable)?
Framing: judge the framework on the OBSERVATIONS (gravitational: rotation curves, lensing/
Bullet, CMB peaks, LSS, BAO), NOT on reproducing ΛCDM's cold-particle construct. Divergences
from ΛCDM are TESTS, not failures. Outcome may be SUCCESS, or an honest TENSION/falsification
(the relativistic-MOND frontier); both are valid science.

## THE CRUX (what this session's dead-ends sharpened — the unexamined channel)
Two DIFFERENT things have been conflated under "screening":
  - the FIFTH FORCE (modification to how mass gravitates): screens at high rho (phaseA:
    ~6% -> ~0%). Established.
  - the STRAIN ENERGY (the Frank elastic energy density sourced by baryon perturbations):
    it is REAL energy, gravitates via E=mc^2 as an ordinary stress-energy source, and does
    NOT obviously screen. Its magnitude is the "heavy elastic solid" magnitude (P7:rem:
    elastic_solid: condensate-derived, order-of-magnitude, CUTOFF-DEPENDENT, "not fixed here").
NEVER EVALUATED at recombination. If the baryon-sourced strain energy density at z~1100 is a
NON-negligible fraction of / comparable to the baryon energy (or larger, if the medium is
"heavy"), it acts as EXTRA gravitating, baryon-tracking, adiabatic clustering matter -> could
drive the 3rd peak WITHOUT a cold particle and WITHOUT a new director scale (Part A moot).
This is the program's make-or-break, and it is CHEAP to bound.

## PHASE 0 — The decisive cheap estimate: strain-energy magnitude + clustering at recombination
GOAL: bound the baryon-sourced Frank strain energy density rho_strain at z~1100 relative to
rho_baryon, and its effective (w, c_s^2, clustering).
TASKS:
  0.1 From the Frank energy F=1/2∫[K(∇n̂)^2] with K set by rho_cond(T_CMB) and the boundary
      conditions imposed by baryon density perturbations delta_b: estimate rho_strain(delta_b)
      at recombination densities. Is rho_strain/rho_b ~ O(0.01)? O(1)? O(5)? (need ~5 to mimic
      Omega_m/Omega_b).
  0.2 Does the strain energy SCREEN? Separate cleanly from the fifth force. If the strain is
      the gravitating source (not the coupling modification), argue whether S_eff/(1+beta rho)
      suppresses it or only the force.
  0.3 EOS/clustering: is the baryon-sourced strain adiabatic (tracks delta_b -> adiabatic, GOOD)
      and pressureless (c_dir->0 -> clusters, GOOD)? Confirm it is NOT the isocurvature splay
      texture (that's the independent, late one).
GATE 0:
  - rho_strain/rho_b too small (<<1) at recombination -> framework underproduces 3rd peak ->
    HONEST TENSION with CMB (a real, reportable result; the MOND frontier). STOP or go to Phase D
    to characterize the tension precisely.
  - rho_strain/rho_b ~ O(1) or tunably ~5 (cutoff-dependent!) -> LIVE: escalate to Phase A/B.
COST: analytic + a short script. Days. DO THIS FIRST.

## PHASE A — Framework linear perturbation theory (the hard prerequisite; only if Gate 0 live)
Derive from the density-feedback + Frank action:
  A.1 modified potential eqn: ∇^2 Phi = 4πG_eff(rho) [delta_rho_b + delta_rho_strain], with
      G_eff + screening from S_eff/(1+beta rho).
  A.2 strain perturbation eqn: delta_rho_strain[delta_b] response (elastic, c_dir->0).
  A.3 each component's (rho, p, c_s^2) through recombination; couple to baryon+photon+neutrino
      Boltzmann hierarchy.
Load-bearing; no shortcut. This is the real intellectual work.

## PHASE B — Level 0 analytic peak estimate (cheap; may decide it before any big run)
Hu-Sugiyama semi-analytic acoustic-peak formulas with the framework's effective driving
(strain energy + any residual force from Phase A). Estimate 1st/2nd/3rd peak heights + the
odd/even alternation. Compare to Planck. GATE B: in range -> Phase C; far off -> tension (Phase D).

## PHASE C — Level 1 parametric Boltzmann (existing tools, no code surgery)
CLASS/CAMB with Omega_cdm=0 + an EFFECTIVE baryon-tracking strain component (parametrised
rho_strain(a, k), EOS, clustering from Phase A). Scan; can it reproduce the observed C_l
(TT/TE/EE, positions + 3rd peak + damping)? Tests viability with existing tools.

## PHASE D — WHERE IT DIFFERS FROM ΛCDM (the falsifiable predictions — the scientific payoff)
Systematically, wherever the mechanism (baryon-tracking elastic strain, no cold particle,
density-feedback, scale-dependent G_eff) departs from a cold collisionless particle:
  D.1 SMALL-SCALE P(k): strain TRACKS baryons -> possibly SUPPRESSED small-scale power / different
      substructure vs ΛCDM. Tests: Lyman-alpha, satellite counts, small-scale lensing.
  D.2 GROWTH (sigma8/S8, fsigma8): elastic/modified growth history -> a specific deviation.
      Relevant to the S8 tension (does the framework predict lower S8 naturally?).
  D.3 BULLET CLUSTER (the sharp one): ΛCDM lensing mass follows collisionless galaxies, offset
      from gas. Does baryon-sourced strain follow GAS (most baryon mass) or galaxies? If forced
      to track total baryons (gas-dominated) -> lensing at the gas -> TENSION; OR a distinctive
      offset signature to derive. CRITICAL TEST — could falsify or distinguish.
  D.4 ISOCURVATURE: any splay-texture admixture -> CMB isocurvature signature (ΛCDM adiabatic;
      Planck constrains tightly).
  D.5 PEAK–Omega_b RELATION: modified driving may give a different 3rd-peak/baryon-density
      dependence than cold wells -> a specific pattern in peak ratios.
  D.6 SCALE-DEPENDENT G_eff(k, rho) from density feedback -> k-dependent growth / CMB-lensing feature.
These make the program FALSIFIABLE and are worth stating even if Phases A-C are hard.

## PHASE E — Full framework Boltzmann (LARGE; only if B/C warrant)
Implement Phase-A eqns in a Boltzmann code (CLASS source modification); compute C_l; fit Planck.
Months. Defer until the cheap phases justify it.

## SUCCESS / FAILURE CRITERIA (state up front — honesty)
  SUCCESS: extended mechanism reproduces C_l (positions + 3rd peak + polarization phase) within
    Planck errors, with distinct falsifiable D-predictions.
  HONEST FAILURE: strain energy too weak / screens completely -> underpredicts 3rd peak ->
    the framework's DM works for galaxies but not the CMB (like relativistic MOND). A VALID,
    reportable outcome that would disfavor the framework's cosmological completeness. Going in,
    accept this is a real possible result.

## COST / RISK / ORDERING (compute is not free — cheap-first)
  Phase 0 = analytic + short script (DAYS) — DECISIVE, DO FIRST.
  Phase B = analytic (days) — second cheap gate.
  Phase A = derivation (cheap compute, hard physics) — only if Gate 0 live.
  Phase C = existing tools (weeks). Phase E = large (months) — LAST.
  RISK: (i) Gate 0 may kill it cheaply (strain too weak). (ii) even if magnitude works, coherence/
  adiabaticity/positions must ALL come out (high bar). (iii) Bullet Cluster (D.3) is an independent
  potential falsifier regardless of the CMB.

## FIRST CONCRETE STEP
Phase 0.1: estimate rho_strain/rho_b at recombination from the Frank energy with K(rho_cond, T_CMB)
and baryon-perturbation boundary conditions. One script + the P7 elastic magnitude. That single
number gates the entire program.

## PHASE 0.1 — DONE (2026-08-31, phase0_strain_energy.py). GATE 0: FAILED (honestly). + a correction.
CEILING: strain energy u = (1/2)K|grad n|^2 <= rho_cond (elastic energy <= medium energy; K~Lam_cond^2,
|grad n|_max ~ 1/xi_cond = Lam_cond -> u_max ~ Lam_cond^4 ~ rho_cond). Numbers at z_rec=1100:
  rho_cond ~ rho_gamma ~ 2.95e-3 eV^4 (radiation-like, Gate 1); rho_b ~ 2.45e-3 eV^4.
  strain ceiling / rho_b = 1.2;  CDM needs Omega_cdm/Omega_b = 5.4 x rho_b.
  => even the ABSOLUTE (unphysical, S=1) ceiling is ~4.5x SHORT; realistically much more.
STRUCTURAL: strain <= rho_cond, rho_cond radiation-like -> at recombination rho_cond ~ rho_b, so strain
capped near rho_b; CDM wants 5x. Channel does NOT close the gap. MOND-frontier tension, now QUANTIFIED.

CORRECTION (important, caught here): the earlier session claim "T_rec/Lam_cond ~ 2000, medium 2000x too
hot at recombination" was WRONG -- it used Lam_cond TODAY. Lam_cond = 0.506 T_CMB TRACKS T, so
T/Lam_cond ~ 2 at ALL epochs. The medium is perpetually MARGINALLY disordered (near-critical), including
at recombination; pretransitional/transient order is NOT exponentially suppressed there. Long-range
ordering is HUBBLE-gated (m_xi~H0 -> z~0), not thermal. => the user's transient/oscillatory angular-order
picture IS physically live at recombination. But it doesn't rescue the magnitude: the strain ENERGY is
capped at rho_cond ~ rho_b. The gap survives on ENERGY, not temperature.

ON OSCILLATION (user): strain tracks + comoves with baryons -> oscillates IN PHASE with them -> adds to
baryon LOADING, not the NON-oscillating CDM driving the 3rd peak needs. Near-critical fluctuations are
also incoherent -> featureless hump unless globally phase-locked (Bell/§7). So oscillatory strain is the
wrong behavior for the gap regardless of magnitude.
ON GALACTIC FLUTTER (user): edge warps/corrugations/rings are elastic MODES of the ORDERED medium (absent
at recombination -> no long-range order). The elastic CONSTANT K~Lam_cond^2 carries over (what we used) ->
gives the 4.5x ceiling. Good cross-scale calibration idea; lands on the same ceiling.

VERDICT: Phase 0 gate FAILED. The strain-energy channel underproduces the recombination matter budget by
>=4.5x, structurally (strain <= rho_cond, radiation-like). This is the honest MOND-frontier result. Next
options: (D) characterize the tension precisely as a falsifiable prediction (how far off is the 3rd peak?),
or accept option 3. Phases A/C/E (full perturbation theory / Boltzmann) are NOT warranted for a channel
that fails the magnitude gate by ~5x. Recommend Phase D framing over a large run.

## PHASE D.5 — the falsifiable prediction stated (2026-08-31): the THIRD PEAK.
Planck data (user-confirmed, correct): peak1/peak2 ~ 2.3 (D_l: ~5750/2500 uK^2); peak3/peak2 ~ 0.98
(~2450/2500) -- third peak nearly EQUAL to second; ~7 peaks fit by 6-param LambdaCDM.
KEY discriminator logic:
  - first-to-second ratio is a BARYON diagnostic (odd/even alternation from baryon loading). BOTH
    LambdaCDM and no-CDM fit ~2.3 by dialing baryons. NOT a CDM discriminator. The baryon density
    used (Omega_b h^2~0.022) is INDEPENDENTLY fixed by BBN + primordial deuterium -> not a free fudge.
  - THE THIRD PEAK is the discriminator. Once baryons are raised to fit peak 2, a NO-CDM model predicts
    peak 3 COLLAPSES (potentials decay w/o cold matter; radiation driving) -> peak3/peak2 << 1 (decaying
    sequence). CDM sustains the potentials -> peak 3 stays high. Observed peak3/peak2~0.98 = CDM signature.
FRAMEWORK PREDICTION: strain capped at ~rho_b + oscillates IN PHASE with baryons (Phase 0.1) => behaves
like a HIGH-BARYON NO-CDM model at recombination => predicts a DECAYING sequence, peak3/peak2 << 1,
vs observed ~0.98. THIS IS THE FALSIFIABLE DIVERGENCE FROM LambdaCDM, and it currently goes AGAINST the
framework. Exact predicted ratio needs a Boltzmann run (Phase C/E, not warranted for a 5x-short channel),
but the DIRECTION (suppressed 3rd peak) is robust textbook physics.
CORRECT a user misconception (gently, recorded): the peaks are SPATIAL scales (k~l/D) at ONE epoch (last
scattering), NOT a time sequence. Peak 3 is a smaller patch than peak 2, not a later moment. "Structure
forms at peak 2, available at peak 3" does not map onto the physics (no time elapses between peaks).
HUMILITY: a modified-gravity theory COULD lift peak 3 with an UNSCREENED potential-sustaining mechanism;
the framework screens to ~GR at recombination + strain too weak (5x) -> cannot, as it stands.
=> Phase D.5 = the honest, quantified MOND-frontier prediction: framework predicts a SUPPRESSED 3rd peak;
   Planck shows a PROMINENT one. A real falsifiable statement to hold against data. Program stands PARKED
   at option-3/tension; Phase D framing (this) is the scientific deliverable, not a large run.

## PHASE D.3 — SCOPED (2026-08-31): the Bullet Cluster (independent of the CMB). Verdict: likely FAILS (amplitude).
OBSERVATION (1E 0657-558, Clowe+ 2006): galaxies (collisionless) sailed through -> on the two sides;
hot gas (collisional, DOMINANT baryon mass ~4x stars) shocked -> left in the MIDDLE; lensing mass -> on the
GALAXIES, OFFSET from the gas at ~8 sigma. Gravitating mass follows the SUBDOMINANT collisionless galaxies,
NOT the dominant gas. Classic killer for baryon-sourced gravity.
THE DOUBLE BIND for the framework's baryon-sourced, comoving strain (both options fail):
  - strain tracks TOTAL baryon mass (gas-dominated): amplitude OK (~5x baryons ~ lensing mass) but LOCATION
    on the GAS -> fails the offset.
  - strain tracks GALAXIES only (compact structures): LOCATION right (galaxies) but amplitude ~5x STARS,
    while clusters need ~6x TOTAL baryons ~ ~30x stars -> ~6x SHORT -> fails amplitude (= MOND cluster problem).
  => a baryon-sourced strain cannot put the RIGHT AMOUNT in the RIGHT PLACE (right amount lives with the gas,
     right place is the galaxies). LambdaCDM satisfies both (collisionless DM: abundant AND follows galaxies).
THREE SUB-TESTS:
  1. TRACKING/location: framework-specific ESCAPE -- each galaxy carries its pre-formed director-strain halo
     (galactic DM), which travels with the galaxy -> lensing follows galaxies -> could PASS the offset
     (better than naive MOND). But forces you into the amplitude-fail row.
  2. AMPLITUDE: strain ~5x stars, need ~30x stars -> ~6x short (MOND cluster problem). FAILS.
  3. COLLISIONALITY: elastic medium w/ shear rigidity ("heavy elastic solid") -> may SHOCK/lag like gas,
     not sail through like galaxies -> could WORSEN the offset. Uncertain, leans against.
NARROW ESCAPE (blocked): a scale/environment-dependent strain, ~6x larger per galaxy in clusters than in
isolation. Blocked by (a) rotation curves CALIBRATE strain to ~5x baryons (cranking 6x risks galaxy fits),
(b) density-feedback SCREENS in dense environments (S_eff->0) -> pushes strain DOWN in clusters, wrong way.
VERDICT: independent of the CMB, likely FAILS on amplitude (~6x short), even if the galaxy-halo feature saves
the offset. CONSISTENT THEME: galaxies work; clusters AND the CMB both need ~5-6x more mass than a
baryon-sourced strain gives. Same ~5x shortfall, three arenas (rotation curves OK; CMB 3rd peak; clusters).
WHAT WOULD SETTLE IT: the actual Bullet lensing + X-ray + optical maps + a concrete strain-sourcing model
(mass-weighted vs gradient/structure-weighted) -> compute strain-induced convergence kappa, compare offset +
amplitude to data. A real (data-driven) calc if pursued; the double-bind argument already gives the direction.

## PHASE D.3 — REVISED (2026-08-31, user correction): local per-source strain dissolves the "double bind."
CORRECTION: the double-bind framing was too binary. Density feedback is LOCAL + PER-SOURCE: each galaxy
sources its OWN feedback and carries its OWN strain halo BOUND to it. Not one global field choosing gas-lump
vs galaxy-lump. This changes 2 of 3 sub-tests in the framework's favor:
  1. OFFSET: per-galaxy halos follow the galaxies IF the diffuse/shocked/incoherent GAS sources little
     COHERENT strain (plausible: splay disclination needs a coherent ordered seed; hot ICM has none).
     -> lensing at galaxies -> PLAUSIBLY PASSES the famous offset. BETTER than naive MOND (which sources
     gravity from gas mass -> lensing at gas -> fails).
  3. COLLISIONALITY: strain halo is BOUND to each galaxy -> travels with it -> collisionless-like -> PASSES.
  => the Bullet Cluster is NOT a qualitative "wrong place" failure. Reduces to a single QUANTITATIVE question.
REMAINING = AMPLITUDE only: per-galaxy strain ~5x STELLAR mass (rotation-curve calibrated, optical radius);
cluster needs ~30x stellar mass (DM 85/gas 12/stars 3). Galaxy halos ~5x stars -> ~15% of cluster -> total
~30% -> ~3x low. = MOND cluster-mass deficit, now ISOLATED from the offset. Same ~5x amplitude theme.
CHECKABLE CAVEATS (determine 3x-short vs closer):
  (a) does hot/shocked/unordered ICM seed a splay disclination (source strain) or not? Offset pass depends on
      "gas sources ~0 coherent strain." If gas sources strain ∝ mass -> lensing back at gas -> offset degrades.
  (b) FULL halos + STRIPPED strain: cluster galaxy halos are tidally stripped -> stripped strain -> intracluster
      strain STILL following the galaxy distribution -> raises amplitude above naive 5x-stellar. Needs modeling;
      could close much of the gap.
REVISED VERDICT: local sourcing RESOLVES the offset + collisionality (a real framework strength, understated
before). Remaining is a QUANTITATIVE amplitude question, leaning short (MOND clusters) but genuinely UNCERTAIN
pending (a)+(b). Much better standing than "fails the Bullet Cluster." Supersedes the double-bind verdict above.

## PHASE D.3 — REFINED again (2026-08-31, user, on the two open flags). Both lean the framework's way.
FLAG (a) gas strain-sourcing: it's the CONDENSATE's local order ("gluon soup"=disordered director), NOT the
gas's EM temperature directly (director is dark, decoupled from the 10^8 K plasma). Shocked/diffuse/turbulent
ICM = no coherent disclination seed + collision stirs the condensate -> condensate stays DISORDERED -> little
strain. Bound galaxy = coherent seed -> strain. => gas sources little strain -> OFFSET RESOLVES. Mechanism =
coherence, not gas heat, but the user's picture is right.
FLAG (b) amplitude, "closer galaxies -> more strain": this is ELASTIC INTERACTION between defects, genuinely
separation-dependent (NOT naive). Overlapping/frustrated director distortions between nearby galaxies -> strain
energy ABOVE the sum of isolated halos -> a dense galaxy field (cluster) sources COLLECTIVE strain > naive
5x/galaxy -> attacks the amplitude gap. Screening does NOT kill it: (1) the halo/ICM environment is DIFFUSE
(beta rho small) -> weakly screened (screening is a galaxy-CORE effect); (2) it's the strain ENERGY (gravitates
via E=mc^2), not the fifth FORCE that screening suppresses.
=> UPDATED: offset likely RESOLVED; amplitude has a real LEVER (collective, weakly-screened strain) to close
the ~3x gap -- how much is the open number. Materially better than "fails on amplitude." The remaining question
is a CALCULATION: strain energy of the N-galaxy configuration at the Bullet's actual galaxy positions vs the
lensing map. Well-posed; the user's two flags are its two physical inputs (coherence-gated sourcing + collective
separation-dependent strain).

## PHASE D.3 — CALCULATION SCOPED (2026-08-31): collective strain vs Bullet lensing map.
GOAL: predict kappa_model(theta) = [Sigma_baryons + Sigma_strain]/Sigma_crit for the Bullet Cluster; compare
OFFSET (kappa peak at galaxies?) + AMPLITUDE (peak kappa value?) to the published lensing map. Decide if
COLLECTIVE strain closes the ~3x amplitude gap.
MODEL: Sigma_baryons from data (optical galaxies + X-ray gas). Sigma_strain from director field n(x) with
COHERENCE-GATED BCs (bound galaxies source; incoherent shocked gas does not -- Flag a). Single-galaxy halo =
SIS (M_enc ∝ r -> rho ∝ 1/r^2), normalized by rotation-curve/stellar mass. COLLECTIVE term = extra Frank energy
from interacting galaxy-induced defects (Flag b).
PHASES (cheap-first):
  1. single-galaxy strain profile rho_strain(r; M_star): SIS, core + tidal cutoff. Cheap (days).
  2a. LINEAR baseline: sum individual SIS at Bullet galaxy positions -> the floor, ~3x short (number to beat). Cheap.
  2b. COLLECTIVE analytic estimate -- THE DECISIVE CHEAP STEP: galaxy splay-defects interact like 2D-Coulomb,
      E_int ~ pi K s^2 ln(R/d). SELF/linear energy ~ N; INTERACTION ~ N^2 -> for N~100s cluster galaxies the
      collective term can dominate. Estimate the N^2 enhancement factor for the Bullet's actual galaxy count/
      spacing. GATE: reaches ~3x -> amplitude plausibly closes; ~1.5x or WRONG SIGN (mutual screening not
      frustration) -> honest shortfall. This single estimate gates the expensive solver.
  2c. FULL nematic solver (weeks, ONLY if 2b promising/borderline): solve n(x) with all galaxy BCs, u=1/2 K|grad n|^2,
      project to Sigma_strain. Multi-defect nematic elasticity PDE. Precision only.
  3. assemble kappa_model, compare OFFSET + AMPLITUDE to the published map. Cheap once 1,2 done.
DATA (all published): Bullet galaxy catalog + stellar masses (optical); X-ray gas map (Chandra); weak-lensing
kappa (Clowe+ 2006 / Bradac+). No new observations.
KEY UNCERTAINTIES: (1) normalization under tidal STRIPPING (does stripped strain still count, following galaxies?);
(2) SIGN of the collective term (frustration enhances vs mutual screening reduces -- 2b must get the sign);
(3) coherence-gating assumption (gas sources ~0 strain; if ∝ mass -> offset degrades).
BOTTOM LINE: Phase 2b (analytic N^2-scaling collective-strain estimate) is the whole ballgame and is CHEAP;
the solver (2c) is gated behind it. Falsifiable either way: collective reaches 3x (Bullet passes) or too small/
wrong-signed (honest amplitude failure). NEXT: run Phase 2b -- the N^2 collective enhancement factor.

## PHASE 2b — DONE (2026-08-31, phase2b_collective_strain.py). Bullet amplitude: SOFTENED, genuinely OPEN.
Framing: splay strain gives isothermal M=k V^2 r/G (same mechanism as galactic flat curves); enhancement =
which V the strain sees. INDIVIDUAL: each galaxy at v_gal~200 km/s. COLLECTIVE: N-galaxy system at cluster
sigma~1000 km/s; full collectivity -> M=sigma^2 R/G = cluster VIRIAL (=observed) mass -> matches BY CONSTRUCTION.
NUMBERS (N~200, v_gal=200, sigma_cl=1000, fractions DM/gas/star=85/12/3, strain~5x stars):
  - LINEAR strain = 5x stars = 15% of total -> 18% of DM -> 5.7x SHORT (shortfall REAL, confirmed).
  - HEADROOM (collective) = (sigma/v)^2 = 25x available; NEEDED = 5.7x -> required efficiency ~23%.
MOND ANCHOR: MOND is partly collective (boost ~ local accel = collective potential) and lands clusters ~2x short
(not 5.7x). So MOND-like partial collectivity softens 5.7x -> ~2x RESIDUAL (the known MOND cluster problem),
better but not fully closed unless framework collectivity exceeds MOND's.
VERDICT: amplitude NOT a hard wall -- ~25x headroom vs 5.7x gap, ~23% efficiency closes it. Likely lands at a
RESIDUAL ~2x (MOND-like), genuinely OPEN. Exact collective efficiency needs the nematic solver (Phase 2c).
BULLET CLUSTER STANDING (net, this session): OFFSET resolved (local per-source, coherence-gated, Flags a);
COLLISIONALITY resolved (halo bound to galaxy); AMPLITUDE softened from 5.7x-fail to ~2x-residual-open (Flag b,
collective isothermal lever w/ ample headroom). Much better than "fails the Bullet Cluster." NOT closed; the
~2x residual is the honest open number, same MOND cluster-scale tension, pinnable only by Phase 2c.

## PHASE 2c — SOLVER BUILT + SMOKE-TESTED (2026-08-31). Awaiting user's full run.
Script: papers/src_paper7/dm_polarization/phase2c_nematic_solver.py. One-constant nematic director-field
relaxation (harmonic-map heat flow): pin hedgehog cores at galaxy positions (coherence-gated; gas not sourced)
+ charge-N radial far-field boundary; relax n(x); F=(1/2)K∫|grad n|^2; eta = F_collective/(N*F_single_trunc),
F_single_trunc = single-defect energy within R=d_nn/2 (spacing-truncated reference -- the PHYSICAL comparison;
a naive full-box reference is an artifact, both directions). Projects u(x) along LOS -> Sigma_strain(x,y) for
Phase 3 lensing. eta>1 = collective enhancement (frustration); eta<1 = mutual screening.
SMOKE TEST (24^3, 4 gal, 300 steps): PASSES machinery (energy monotone-decreasing, outputs written). eta values
are NOT physical (tiny/packed config, R_trunc ~ core). Needs python3.11 (default python3=3.14 has no numpy).
FULL RUN (user to do): python3.11 phase2c_nematic_solver.py --grid 128 --ngal 200 --steps 8000 --positions
bullet_galaxies.npy. REQUIRED real input: Bullet galaxy positions (Ngal x 3, grid units) via --positions (mock
otherwise). matplotlib absent -> PNG skipped, .npy saved. Outputs -> phase2c_output/summary.txt (eta_PHYSICAL etc).
INTERPRETATION when done: compare eta_PHYSICAL to the ~23% collective efficiency Phase 2b said closes the ~5.7x
amplitude gap (i.e. eta needs to reach ~5.7 to fully close; ~1-2 = MOND-like residual; <1 = screening, fails).

## PHASE 2c — SMOKE TEST FINDING (2026-08-31): dimensionless eta is ILL-POSED; needs calibrated mass profile.
Machinery VERIFIED (relaxation converges, monotone energy, outputs, resolution diagnostic d_nn/core works;
simulated bimodal Bullet-like positions fine -- REAL positions only needed for the Phase-3 offset fit, not eta).
BUT: a resolved run (G=64, N=16, d_nn/core=14.5) gives eta spanning 0.4->12 depending on definition:
  - eta_local (radius-matched, WS ball) = 0.44  -> near-galaxy SCREENING (neighbours align, reduce local splay).
  - eta_fullbox (vs truncated ref) = 11.9        -> charge-N FAR-FIELD ~N^2 energy (outside cluster, spurious).
  - eta_boxref = 0.03                            -> artifact.
=> the dimensionless "collective enhancement" is DEFINITION-SENSITIVE (which radius / which linear reference) and
   does NOT have a clean answer in the one-constant model. Reason: energy isn't additive for overlapping fields;
   near-field screens while far-field (charge-N) grows ~N^2; the physical total-within-cluster is aperture-dependent.
ROBUST REFORMULATION (the real Phase 2c): compute the ABSOLUTE strain MASS profile M_strain(r) for the cluster,
CALIBRATE K from a single galaxy's rotation curve (M_strain,gal(r_opt)=~5 M_star), convert the cluster's
dimensionless energy profile to mass, compare M_strain,cl(R_cluster) to the observed cluster DM. Report the RATIO
(not a dimensionless eta). ALSO: the one-constant HARMONIC model omits the framework's real physics -- splay
instability K_splay<0, nonlinear stabilization, density-feedback screening -- all of which shift the number.
STATUS: Phase 2c is a genuine solver-development + calibration PROJECT, not a parameter plug-in. The Bullet
AMPLITUDE therefore stays where Phase 2b left it: OPEN, ~2x MOND-like residual, ample analytic headroom, but the
actual efficiency needs the calibrated NONLINEAR treatment to pin. Do NOT quote the one-constant eta as physical.
DECISION: invest in the calibrated-mass + nonlinear reformulation, or park the Bullet amplitude at the honest
~2x-residual-open verdict (Phase 2b) with Phase 2c machinery in place for later.

## PHASE 2c / A1 — KEY FINDING (2026-08-31): M~r halo is NOT a minimum of positive Frank energy.
Built A1 disclination/hedgehog solver (phase2c_A1_disclination.py), reusing gradient_flow_constrained.py's
semi-Dirac Frank energy (K=J2a+mu* J2iso, sin^4 theta), S^2-projection, slerp angle-clamp, + multi-galaxy
oriented (oblate-hedgehog) sources + isotropic-vs-oriented N=2 comparison + an M(r)-slope validation gate.
Smoke on base-env torch: machinery RUNS (autograd, monotone energy, comparison). BUT the validation GATE FAILS:
single-galaxy M(r) slope = 0.2 (collapsed) to 2.0 (boundary-forced), never ~1.
DIAGNOSIS (physics, not a bug): a hedgehog (radial splay, degree-1) SHRINKS under the positive Frank/Dirichlet
energy in 3D (Derrick's theorem) -> relaxation drives it toward uniform/collapse, AWAY from the extended M~r
profile. Adding J4 (Option 2) stabilises a FIXED-SIZE soliton, still not M~r to galactic radii.
=> M~r is DRIVEN by the SPLAY INSTABILITY K_splay<0 (negative splay energy extends the radial texture until a
   cutoff stops it). That cutoff = the acceleration scale, "NOT fixed here" (Paper VII). So the halo energy is
   UNBOUNDED/CUTOFF-DEPENDENT by construction, and the M~r EXTENT -- hence the CLUSTER AMPLITUDE -- is intrinsically
   tied to the undetermined cutoff. It is NOT a single derivable number; it is cutoff-parametrised.
CONSEQUENCE: Option 1 (positive Frank) is physically INCAPABLE of the halo (wrong energy). A correct solver needs
the K_splay<0 energy (frank_run.py) + an EXPLICIT cutoff, and its output is a FUNCTION of the cutoff, not a number.
=> The deep dive CONFIRMS (does not resolve) the earlier verdict: the cluster amplitude cannot be pinned without
   specifying the acceleration-scale cutoff Paper VII leaves open -- exactly the "~2x MOND-like residual, not fixed
   here" status. phase2c_A1_disclination.py is kept (machinery + the Derrick finding); the physical A1 would be a
   K_splay<0 + cutoff solver whose result is cutoff-parametrised.
DECISION for the user: (a) build the K_splay<0 + explicit-cutoff solver (harder; answer is a cutoff-family, not a
   single number), or (b) accept the amplitude is cutoff-set (Paper VII's own position) and hold the honest verdict.
