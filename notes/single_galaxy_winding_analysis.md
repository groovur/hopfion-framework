# Single-galaxy director-strain winding — analytic derivation + the sign tension

Goal: pin the M~r halo physics before more solver code (A1). Result: the winding
profile + M~r check out, but the MEASURED stiffnesses expose a genuine tension in the
framework's own energetics that must be resolved to define the gravitating mass.

## Setup: +1 disclination with character angle psi
A radial-BC galaxy sources a +1 disclination. Character psi = angle(director, radial):
psi=0 radial (splay), psi=pi/2 azimuthal (bend). For alpha = phi + psi(rho):
  splay = cos(psi)/rho,   bend = sin(psi)/rho.
Frank energy density:
  u = (1/2)[ K_splay cos^2(psi) + K_bend sin^2(psi) ] / rho^2.
Energy per unit length:
  E(psi) = pi [ K_splay cos^2(psi) + K_bend sin^2(psi) ] ln(R/a).
=> u ~ 1/rho^2 -> M_enc ~ r (the paper's theta~ln r, |grad n|~1/r, flat rotation curves). CHECKS OUT.

## The tension (measured stiffnesses)
frank_run.py test: K_splay = -8.45 (UNSTABLE), K_bend = +0.534 (STABLE). Ratio ~ -16.
Minimising E(psi): K_splay < K_bend and NEGATIVE -> minimum at psi=0 (radial splay) with
  E = pi K_splay ln(R/a) < 0  -> NEGATIVE Frank energy -> negative mass excess -> ANTI-gravitates.
The positive-energy option (psi=pi/2, bend, E=pi K_bend ln>0) is DISFAVOURED (splay is what's unstable).
=> Cannot get positive gravitating mass AND K_splay<0 from the same (quadratic) Frank energy.
Precise conflict:
  - frank_run: K_splay<0 (good: explains WHY the texture forms spontaneously);
  - paper flat_curves: u ~ +K/r^2 (positive mass);
  - energy-min texture with K_splay<0 is NEGATIVE-energy radial splay.

## Three resolutions (the real open question)
1. GRAVITATING MASS != Frank-energy-relative-to-uniform. It is the coupling energy of holding the
   baryon BC, or the absolute reconfiguration cost (positive), not the (possibly negative) Frank term.
   -> the paper's "the Frank energy gravitates" needs sharpening.
2. K_splay<0 is only a FORMATION driver; the halo is a METASTABLE positive-energy bend winding (psi=pi/2),
   held by topology/dynamics, not the ground state. Mass = +K_bend (positive) but not the energy minimum.
3. A higher-order/stabilising term (J4, density feedback) makes the STABILISED texture's TOTAL energy
   positive though the quadratic splay stiffness is negative -- like a negative mass^2 field whose
   soliton is positive-energy after the quartic stabilises it.

## Assessment / recommendation
Most likely (1) or (3): the quadratic stiffness being negative is an INSTABILITY, consistent with the
STABILISED object having POSITIVE energy. Then the DM mass = the full stabilised energy (Frank + J4/feedback),
POSITIVE, and K_splay<0 is the seeding instability, NOT the sign of the halo mass.
=> If the framework confirms this, the A1 solver energy must be the STABILISED functional and the mass is its
   positive total -- which also fixes exactly what to compute. This is the precise decision that unblocks A1.
OPEN QUESTION FOR THE FRAMEWORK: is the DM mass the (positive) energy of the fully-stabilised texture (with
K_splay<0 as the seeding instability), rather than the (negative) quadratic Frank energy? Needs a framework-
level answer (Paper VII / DM papers), not a search.

## RESOLUTION (2026-08-31, user): the POLARON picture. Sign paradox dissolved; mass = drag.
User: resolution (1) makes sense, and the negative energy may relate to non-visible DRAG forces. -> This IS
the polaron: a baryon dressed by the strain distortion it induces.
  - K_splay<0 -> baryon induces a FAVOURABLE splay -> NEGATIVE binding energy (like a polaron).
  - BUT the polaron's EFFECTIVE INERTIAL MASS is LARGER than bare (it must DRAG its distortion cloud when it
    moves). Negative binding + increased inertia coexist -- the defining polaron signature.
=> the DM "mass" is the POSITIVE extra INERTIAL mass of baryon+strain-cloud, NOT the (negative) static Frank
   energy. Sign paradox dissolved (inertial mass positive by construction). Resolution (1), made precise.
DRAG = the same dressing, dynamical:
  - accelerating dressed baryon drags the cloud -> extra inertia -> the rotation-curve "dark mass";
  - baryon moving through the medium -> strain cloud LAGS (c_dir->0 slow mode) -> wake -> DRAG =
    DYNAMICAL FRICTION (satellite orbital decay, sinking globulars, LMC orbit -- standardly attributed to DM).
  - matches P7:rem:elastic_solid ("heavy elastic solid, shear rigidity, c_dir->0"): large effective inertia +
    strong drag. "Dark mass" and "invisible drag" = two faces of one strain dressing.
IMPLICATION for A1: compute the EFFECTIVE INERTIAL MASS (dynamical response of the strain to baryon motion),
NOT the static Frank energy. Positive by construction -> sidesteps the sign problem, but it is a DYNAMICAL calc
(resistance to acceleration), harder than a static energy. This is why the static-energy solver fought the sign.
PHENOMENOLOGY (Phase-D, falsifiable): polaron/drag DM is NOT pure collisionless CDM -> velocity/acceleration-
DEPENDENT effects, esp. modified DYNAMICAL FRICTION + orbital-decay signatures, testable vs the systems where
dynamical friction is measured. A genuine CDM-distinguishing prediction. [SPECULATIVE lead; not yet derived.]

## POLARON EFFECTIVE MASS + DRAG — analytic derivation (2026-08-31). Sign paradox RESOLVED rigorously.
Director dynamics (framework's own c_dir^2=K/chi): chi d_t^2 n = K grad^2 n + source.
Comoving distortion n(x-Vt) => d_t n = -(V.grad)n. Field momentum (pi=chi d_t n):
  P_i = -INT pi d_i n d^3x = chi V_j INT (d_i n . d_j n) d^3x  ==  M_ij V_j.
EFFECTIVE INERTIAL MASS TENSOR:  M_ij = chi INT (d_i n . d_j n) d^3x ;  m_eff = (chi/3) INT |grad n|^2 d^3x.
KEY RESULTS:
  - POSITIVE BY CONSTRUCTION: INT|grad n|^2 is a sum of squares, chi>0 (susceptibility). SIGN OF K_splay NEVER
    ENTERS. The negative Frank energy is irrelevant to the mass. Sign paradox dead.
  - Instability (splay: c_dir^2=K_splay/chi<0) is K<0, NOT chi<0 -> config GROWS (halo forms), inertia chi>0 ->
    POSITIVE mass. Formation and mass cleanly separated.
  - Flat curves: |grad n|^2~1/r^2 -> INT|grad n|^2 d^3x ~ r -> m_eff(r) ~ r -> M_enc ~ r. Positive mass.
  - Heavy: chi=K/c_dir^2 large as c_dir->0 -> huge dressing -> "heavy elastic solid" (P7:rem:elastic_solid).
=> DM mass = chi * (positive Dirichlet integral of the comoving distortion). This is EXACTLY what the static
   solver already computes (INT|grad n|^2 = 2 J2iso). A1 REVIVED + tractable: relax the extended winding, compute
   INT|grad n|^2 (positive), x chi/3. Cutoff still sets MAGNITUDE (halo extent); SIGN resolved.
DRAG (same cloud, dynamical):
  - RIGOROUS: the effective-mass halo gravitates -> STANDARD gravitational dynamical friction (=CDM orbital
    decay). Framework reproduces observed drag via the halo; no tension.
  - EXTRA (needs gamma_rot, unpinned): finite c_dir -> cloud lags/dissipates -> elastic drag F~gamma_rot V.
    c_dir->0 exactly (lossless) -> rigid ride -> pure mass, NO extra drag. So elastic drag in [0, small], open.
  - DISTINGUISHING: gravitational friction F~1/v^2 (Chandrasekhar); elastic/viscous F~V. Different v-scaling.
HONEST vs GR+CDM: standard dynamical friction is well-explained; elastic drag is a PREDICTION, not a required
anomaly. Puzzles where LESS friction than CDM is seen (Fornax GC timing; "fast bars") are the place to look, but
debated (cored halos also explain them). [SPECULATIVE until gamma_rot pinned.]

## A1 chi-FREE RESULTS (2026-09-01, phase2c_A1_disclination.py, base-torch). Config NAILED; isotropic SCREENS.
Polaron mass proxy = INT|grad n|^2 (positive, chi/3 cancels in ratios).
1. M~r CONFIG NAILED: single-galaxy slope = 1.04-1.06 (grid 48-72) with MINIMAL relaxation. KEY: FULL relaxation
   of the positive Frank energy DEGRADES it (slope->1.7-2.1) = the Derrick collapse, seen directly. Analytic
   hedgehog is correct; over-relaxing the (wrong-sign) positive energy was the bug. Use minimal relaxation.
2. COLLECTIVE ENHANCEMENT (chi-free, ISOTROPIC hedgehogs), eta = m(N2)/(2 m(N1)), stable across steps:
     sep=8 ->0.559 ; sep=14 ->0.632 ; sep=22 ->0.728 ; sep=32 ->0.819.
   => eta<1 ALWAYS (SCREENING); MORE screening when CLOSER; ->1 (additive) as sep->large. Two radial hedgehogs
   point oppositely between them -> gradients CANCEL -> closer = more cancellation.
IMPLICATION: isotropic clustering REDUCES the mass (wrong way for the amplitude). Enhancement CANNOT come from
isotropic sources -> must come from ORIENTATION. Confirms "isotropic averaging misses it" (base file section 10)
NUMERICALLY. => chi_splay/chi_bend is NOT optional: it is REQUIRED, the only thing that can flip screening into
enhancement (oriented disks add where spherical hedgehogs cancel). User's "chi nice-to-have" was too modest.
NEXT: chi_splay/chi_bend from the semi-Dirac band (extend frank_run.py: static rotational susceptibility), then
the oriented (oblate, aligned-vs-perp) comparison -> does orientation reverse the screening?

## CHI + ORIENTED COMPARISON (2026-09-01, frank_chi.py + phase2c_A1_disclination.py). Screening NOT reversed.
CHI (frank_chi.py, diamagnetic q=0 term of frank_run): chi(G1, validated generator) = 601 (cutoff-dep, Lam=20).
Anisotropy check chi(G2)/chi(G1) ~ 2.69 -- SUGGESTS chi is anisotropic, BUT G2 is a MODELLING GUESS for the 2nd
rotation DOF (framework's true 2nd generator not confirmed). So: isotropic chi solid; chi ANISOTROPY tentative.
ORIENTED COMPARISON (config-based, chi-free at leading order; N=2, sep=14, oblate q=0.3, minimal relaxation):
  eta = m_polaron(N2)/(2 m_polaron(N1)), single-galaxy slope clean (~1.05 sphere; 1.38 oblate, config less clean):
    isotropic (sphere)    : eta_iso     = 0.632
    oblate disks, aligned : eta_aligned = 0.567
    oblate disks, perp    : eta_perp    = 0.497
  => orientation does NOT reverse the screening -- it slightly WORSENS it (perp worst). ALL eta<1.
IMPLICATION: at the level we can compute (config + isotropic chi), collective baryon-sourced strain SCREENS
(clustering REDUCES per-galaxy mass) regardless of orientation. This REFUTES the Phase-2b analytic optimism
("ample headroom, ~23% efficiency"): the actual computation shows the collective effect goes the WRONG way.
=> the cluster amplitude gap does NOT close via collective strain -- it worsens. Option 3 / amplitude-short
   verdict stands, now with the collective effect SHOWN to screen. ONLY remaining lever = the chi ANISOTROPY
   (tentative ~2.7), which needs the framework's true 2nd rotation generator to settle -- the "only get that far"
   boundary the user anticipated.
CAVEATS: N=2, one separation, minimal relaxation, oblate model (q=0.3) is a modelling choice with slope~1.38
   (less clean), chi anisotropy unconfirmed. Suggestive, not definitive; but the screening is robust across
   sep/steps for the isotropic case and the oriented cases don't reverse it.
NOTE: oblate single-disk absolute mass (1842) > sphere (700) -- oblate concentrates gradients -> MORE mass per
   galaxy (model-dependent), even though the enhancement RATIO screens. Absolute amplitude depends on the
   (uncertain) single-galaxy config; the ROBUST result is the screening RATIO.

## SHADOW / SCREENING TEST (2026-09-01, user's point). Shadow WORSENS the screening.
User: are we accounting for galaxies BLOCKING field energy (density-feedback screening S_eff->0 inside high-rho
galaxies -- "shadow of the anisotropy")? Real omission -- point-source model missed it. Tested via mass-weighting
(exclude screened interior: S_eff=(1-exp(-(r/rgal)^2)) per galaxy):
  rgal=0: eta_iso=0.632 eta_al=0.567 eta_pe=0.497
  rgal=4: eta_iso=0.589 eta_al=0.545 eta_pe=0.474
  rgal=6: eta_iso=0.546 eta_al=0.521 eta_pe=0.453
=> shadow LOWERS eta (worse screening). Reason: the near-galaxy interior holds each halo's STRONGEST,
   NON-cancelling |grad n|^2 (the good individual mass); the cancellation is in the BETWEEN region, OUTSIDE the
   screened zones. Removing interiors strips the good mass, leaving the cancellation relatively stronger.
SCREENING IS ROBUST: isotropic + oriented + screened-interior ALL give eta<1. Deep reason: two radial
baryon-sourced halos point OPPOSITELY between them -> gradients cancel there -- TOPOLOGICAL, not an artifact.
UNTESTED (harder, uncertain): (a) CONFIG-level shadow (source at galaxy EDGE, interior free, field routes AROUND)
-- changes the config not just the weight; skeptical it reverses (between-galaxy fields still oppose). (b)
DIRECTIONAL shadow (semi-Dirac anisotropy has a direction -> directional blocking) -- couples to the chi
anisotropy, needs the 2nd generator (unavailable).
BOTTOM LINE: collective baryon-sourced strain SCREENS robustly. Cluster amplitude gap does NOT close via
collective strain (worsens). Remaining levers both need the unavailable 2nd generator / chi anisotropy.

## RELAXATION EXPLORATION (2026-09-01): minimal relaxation IS correct; screening is robust.
Q (user): are we relaxing the right way? Over-relaxing positive Frank -> Derrick collapse. Tested the PHYSICAL
K_splay<0 (+ J4) energy (aniso_energy in phase2c_A1_disclination.py, --aniso; K_splay=-8.45, K_bend=+0.534):
  positive-Frank, FULL relax (steps400): slope 2.19 (COLLAPSE).
  K_splay<0 + J4, FULL relax: slope 2.99, m blows up to 1.3e5 (RUNAWAY -- unbounded splay concentrates).
  J4 scan lam=0.1..8: slope ~3.05-3.08 ALWAYS, m~7e4 -- J4 does NOT bound the runaway at any strength.
  minimal relaxation (analytic hedgehog): slope ~1.05 (clean M~r).
FINDING: the M~r halo is a SADDLE of BOTH energies (positive Frank collapses; K_splay<0 runs away). It is NOT a
minimum of any simple local energy. It is stabilised by the ACCELERATION-SCALE CUTOFF (halo extends until g~a0,
where S_eff screens -- position-dependent, = the user's SHADOW/screening physics at the halo EDGE), which our
energies lack. Without it K_splay<0 just runs away.
=> MINIMAL RELAXATION IS CORRECT, not a hack: we IMPOSE the physical M~r config (analytic hedgehog) because the
stabilising physics lives OUTSIDE the simple energy; relaxing a simple energy MOVES AWAY (collapse/runaway).
=> the SCREENING result (eta<1) came from the CORRECT (imposed, minimal-relax) configs -> robust, NOT an artifact
of wrong relaxation. Screening = geometric: two M~r halos' radial gradients oppose between galaxies -> cancel.
IMPLICATION: to relax "properly" one would need the position-dependent acceleration-cutoff (S_eff at the halo
edge) in the energy -- a bigger model. But it would not change the between-galaxy cancellation. Screening stands.

## ISOTROPIC-INJUSTICE CHECK (2026-09-01, user's point 1): NOT a significant injustice; screening robust.
Decomposed the config into splay^2/twist^2/bend^2 (frank_decompose) -> anisotropic mass = chi_s(splay+twist) +
chi_b*bend; scanned rho=chi_b/chi_s and solved eta=1. (grid64, sep14, minimal relax, isotropic hedgehogs):
  single: splay=1.13e4  twist=1.9e1  bend=8.2e2
  N2    : splay=1.15e4  twist=2.1e1  bend=8.6e2
  CHANNEL enhancements: eta_splay = 11460/(2*11320) = 0.506 ; eta_bend = 855/(2*817) = 0.523. BOTH <1 (screen).
=> the anisotropic mass is a weighted average of 0.506 and 0.523 -> ALWAYS ~0.51 <1 for ANY positive chi_b/chi_s.
   NO chi anisotropy reverses the screening. Reversal would need chi_b/chi_s ~300 (clean bend) -- unreachable.
REASON: config is SPLAY-DOMINATED (splay ~11000 >> bend ~800). Two radial halos cancel in SPLAY, generating
   little physical BEND -> the anisotropy has nothing to reweight. => isotropic mass was NOT a significant
   injustice; the screening is robust to the splay/bend anisotropy.
CAVEAT: perfect hedgehog has bend=0 (pure splay); the ~817 single bend is partly NUMERICAL (discretization +
   pinned core) -> bend channel precision limited. But the SPLAY channel (0.51, clean) dominates the screening,
   so the conclusion is robust. A cleaner test needs better numerics (physical bend cleanly separated).

## SECOND GENERATOR = EDGE SCREENING (2026-09-01, user's point 2). UNIFICATION.
S_eff = sin^4(theta)/[phi^6 (1+beta rho)] (Axiom A2); theta = director tilt from e_r=grad rho/|grad rho|.
Two director-rotation DOF split by action on theta:
  G_phi (rotation ABOUT e_r, azimuthal): does NOT change theta -> S_eff invariant -> SOFT Goldstone = frank_run's
    generator (chi=601). Easy (no anisotropy cost).
  G_theta (POLAR TILT, the 2nd generator): CHANGES theta -> changes S_eff=sin^4(theta) directly. chi_theta set by
    the sin^4 curvature AND modulated by 1/(1+beta rho). => chi_theta is DENSITY-DEPENDENT: large where UNSCREENED
    (low rho, halo EDGE), suppressed where screened (high rho, galaxy interior). Unlike chi_phi (density-blind),
    chi_theta LIVES AT THE HALO EDGE.
=> the 2nd generator is not merely related to the edge screening -- it IS the same structure. UNIFICATION of
   previously-separate walls:
   - 2nd generator chi_theta  == the tilt tied to sin^4 theta
   - edge screening 1/(1+beta rho) == the density modulation of that same chi_theta
   - M~r cutoff (acceleration scale) == where the edge screening kicks in
   - missing relaxation stabiliser == the same density-dependent edge term
   ALL ONE THING: the halo-edge density-dependent anisotropy. It (a) stabilises M~r (stops the K_splay<0 runaway
   at the edge), (b) IS chi_theta, (c) IS the shadow/screening. Explains why simple energies couldn't stabilise
   the halo (lacked the edge term) and why isotropic mass missed the anisotropy (chi_theta is EDGE-localised;
   the splay-dominated bulk config had little edge structure).
HONEST LIMIT: cannot compute chi_theta from frank_run directly -- needs the framework's H dependence on the POLAR
   TILT (how theta enters d; a sigma_z-like gap?), NOT confirmed in the repo. Guessed G2 (~2.7) was a stand-in.
WELL-POSED NEXT TASK (unified): model the halo-edge density-dependent anisotropy -- sin^4(theta)/(1+beta rho) as a
   POSITION-DEPENDENT weight in the energy. One ingredient simultaneously (a) stabilises M~r in full relaxation,
   (b) supplies the real chi_theta, (c) implements the shadow/screening. The three walls were one wall.

## EDGE-ENERGY BUILD (2026-09-01): edge term does NOT stabilise M~r locally; stabilisation is GLOBAL.
Built the ê_r-keyed aniso energy (K_splay<0) + EDGE cutoff (splay_wt: K_splay<0 active only inside R_c). Modeling
anchor (user): halos are SPHERICAL as observed w/ GR -> local-radial ê_r (disc-compression would conflict). Test:
does the edge term stabilise M~r under FULL relaxation?
  aniso NO edge, full relax: m=1.26e5 slope=2.99 (runaway).
  aniso EDGE R_c=15, full relax: m=1.11e5 slope=2.88 (STILL runaway).
FINDING: the edge term does NOT fix it. WHY: K_splay<0 drives |grad n| up EVERYWHERE (a high-gradient "splay foam",
slope3 = ~uniform energy density) because the splay energy is UNBOUNDED below; cutting it at large r just moves the
foam INSIDE R_c. The runaway is DISTRIBUTED, not an edge over-extension -> a LOCAL edge term addresses the wrong end.
DEEPER: the extended M~r halo is held by a GLOBAL self-consistent condition -- the halo extends until its OWN
gravity's acceleration g=GM(r)/r^2 ~1/r drops to a0 (cutoff radius depends on enclosed mass -> non-local). A LOCAL
gradient-descent solver CANNOT produce it (positive Frank collapses; K_splay<0 foams; J4/edge-cutoff don't tame it).
=> the user's edge-screening insight is PHYSICALLY correct (the acceleration edge IS the stabiliser) but it is a
   GLOBAL self-consistency condition, NOT a local energy term. Hence IMPOSING the analytic hedgehog (minimal
   relaxation) is the CORRECT approach (confirmed again) -- the local solver can only distort the global-held config.
NET (whole halo thread): M~r = saddle of all local energies, held by the global acceleration scale (impose it);
halo spherical (local-radial e_r, matches round observed halos); collective strain SCREENS robustly (isotropic +
oriented + shadow + anisotropy-weighted); 2nd generator = edge screening = M~r cutoff = one halo-edge physics, but
GLOBAL not a local generator. To do it "properly" needs a SELF-CONSISTENT solver (halo gravity sets the cutoff),
a much bigger model. phase2c_A1_disclination.py carries the machinery + all these findings.

## FRAMING CORRECTION (2026-09-01, user): screening was OVER-EMPHASISED; dominant DM = RAR, not the inverse.
User: clusters have MORE DM; the screening seems to predict the INVERSE of observed. Correct catch -> I
over-emphasised a SUBDOMINANT effect. Two contributions were conflated:
  (1) SINGLE-OBJECT strain (each galaxy's own M~r halo) = the DOMINANT DM, universal acceleration scale a0 =
      the RADIAL ACCELERATION RELATION (RAR). Matches galaxies AND dwarfs.
  (2) GALAXY-GALAXY collective strain (the eta<1 screening I computed) = a SUBDOMINANT correction on top of (1).
I let (2) drive the narrative; it is NOT the prediction.
DOMINANT prediction (1) gives the OBSERVED direction, not the inverse:
  - LOW-acceleration systems have the MOST DM. Dwarf spheroidals (low accel, low density) have the HIGHEST DM
    fractions (DM/stars ~10-1000). Universal-a0 strain: lower accel -> more boost -> more effective DM. => the
    framework gets dwarfs RIGHT, OPPOSITE of "denser->less". RAR matches galaxies + dwarfs.
CLUSTERS: UNDER-prediction, NOT inversion. Universal-a0 strain calibrated at galactic v~200 applied to a cluster
  sigma~1000 falls short by ~(200/1000)^2 = the MOND CLUSTER PROBLEM (~2-5x short). A DEFICIT, not a reversal --
  clusters get SOME strain, just not enough.
NET across mass: dwarfs OK (high DM matched), galaxies OK (RAR), clusters SHORT (MOND problem). The genuine open
issue is the CLUSTER DEFICIT (shortfall in the dominant single-object strain), NOT the galaxy-galaxy screening.
LESSON: the eta<1 screening is real but a SMALL correction to a SUBDOMINANT term; do not present it as the
framework predicting less DM in dense environments. The framework predicts the RAR (more DM at low accel) + the
cluster deficit.

## DENSITY-FEEDBACK CONSISTENCY + CLUSTER WAKE (2026-09-01, user: "cluster cores like planets").
CONSISTENCY (user insight): density-feedback screening (S_eff->0 at high rho -> DECOUPLE -> GR) is ONE monotonic
law spanning ~30 dex in density:
  planets/solar system (rho~1 g/cc, beta rho~1e6): STRONGLY screened -> NO DM (solar-system requirement MET by
    the same mechanism). Cluster cores/disks: mildly screened. Dwarfs/outskirts/voids (lowest rho): unscreened ->
    MAX DM (dwarfs = highest DM fraction). => lands both anchors (solar-system no-DM, dwarf max-DM) with one law.
  EXTENDS to a prediction: screening dense CORES -> CORED DM profiles (not NFW cusps) -> matches the CORE-CUSP
    problem (dwarfs+clusters show flat cores). Distinctive, right direction (debated: baryonic feedback also cores).
  => density-feedback gets the DM DIRECTION + PROFILE right (planets, dwarfs, cores, RAR). This is the framework's
     GALAXY-SCALE SUCCESS. Separate from the cluster TOTAL-MASS deficit (magnitude).
CLUSTER WAKE (cluster_wake_estimate.py): the deficit lever (point 1, drag/wake), MOND-distinct.
  Mechanism: c_dir->0 -> every galaxy SUPERSONIC (v>>c_dir) -> sheds a strain WAKE it can't drag; c_dir->0 -> shed
  strain ~never relaxes -> PERSISTS as intracluster strain -> ACCUMULATES over cluster history.
  Numbers: deficit ~5.6x (need ~28x stellar; RAR gives ~5x). M_icl ~ N_gal*f*M_halo*N_crossings, N_crossings~10.
  Closes the deficit for deposition fraction f~0.4-0.5 PER CROSSING (NOT few percent -- SUBSTANTIAL: each galaxy
  sheds ~half its continuously-re-sourced halo per crossing). Plausible for v>>c_dir (halo can't follow at all)
  but LARGE. => RIGHT ORDER OF MAGNITUDE (not orders off), unlike screening (wrong sign). A LIVE lever, borderline-
  reachable, framework-specific. THE UNKNOWN = f, needs the actual director wake calc (v>>c_dir response).
NET: density-feedback = galaxy success (RAR, planets, dwarfs, cored profiles); cluster deficit = open, with the
wake/drag as a promising MOND-distinct lever (needs f from the wake calc). NOT the inverse of observed -- the
framework has the right DM direction across scales; only the cluster magnitude is short, with a live mechanism.

## WAKE CALC -- CAREFUL (2026-09-01, cluster_wake_calc.py): OVERTURNS the optimistic estimate. Deficit STANDS.
Set up the moving-source director response (chi d_t^2 n = K grad^2 n; Mach cone at cos theta_k=c_dir/v). Two
deposition channels, BOTH fail to ENHANCE beyond the RAR value:
  (A) RADIATED wake (Cherenkov): F_drag ~ m_eff v c_dir/L -> scales WITH c_dir -> as c_dir->0 the medium is
      FROZEN, no propagating modes, NO radiation, NO dissipation -> F_drag->0. + radiated mass=E/c^2~(v/c)^2~1e-6.
      Negligible either way.
  (B) CONFIGURATION (galaxy leaves its frozen halo behind): deposits ~M_halo, GOOD -- but re-forming a new halo
      takes t_form~R_halo/c_dir->inf. So each galaxy deposits its ONE halo ONCE, then goes halo-less; no re-source,
      no accumulation. Total ~ N_gal*M_halo = RAR value (~5x stellar) -> STILL 5.6x SHORT. Redistribution, not
      amplification.
KEY (elegant): c_dir->0 -- the very COLDNESS that makes the DM 'cold' -- is EXACTLY what forbids the wake
enhancement (frozen -> no dissipation (A); infinite re-forming time -> single deposition (B)). The optimistic
f~0.5 (prev turn) assumed repeated shed-and-resource, which c_dir->0 FORBIDS. Self-corrected.
=> WAKE LEVER FAILS. Cluster deficit STANDS -- the framework shares the MOND cluster problem, and neither
   framework-specific lever resolves it: density-feedback SCREENS clusters (galaxy success, not cluster fix);
   wake/drag does NOT enhance (c_dir->0 forbids it).

## OVERALL DM SCORECARD (honest, end of this investigation)
GALAXIES/DWARFS/PLANETS: SUCCESS. Density-feedback = one law, right direction across ~30 dex: solar-system no-DM,
  dwarf max-DM, RAR, cored profiles. NOT the inverse of observed. The polaron effective mass gives positive M~r
  (P7:rem:polaron). Galaxy-scale DM is a genuine particle-free account.
CLUSTERS: SHORT (open). ~5.6x total-mass deficit = the universal-a0/MOND cluster problem. Framework-specific
  levers explored and FAILED: (a) collective strain SCREENS (not enhances); (b) density-feedback SCREENS clusters;
  (c) wake/drag does NOT enhance (c_dir->0 forbids). No framework fix found. Honest: shared MOND cluster deficit.
RECOMBINATION CDM: input (option 3, separate thread; fabric_cmb_hypothesis_and_plan.md).
NET: the framework's DM is a real success at galaxy/dwarf/planet scales (RAR + density-feedback direction), with a
genuine, unresolved cluster-magnitude deficit shared with MOND. Not the inverse; a shortfall at one scale.

## FRAMEWORK vs MOND + CLUSTER-DEFICIT RECOVERY (2026-09-02, user).
FRAMEWORK != MOND: MOND is a FORCE LAW (boost slaved to local baryonic g_N, universal a0). The framework is a
physical FIELD (director strain) with DOF MOND lacks: (i) condensate-state-dependent response (density-feedback +
anisotropy, not universal a0); (ii) the strain can DECOUPLE from baryons (relic/fossil deformation = "stress
without the population"); (iii) INERTIA (polaron effective mass); (iv) TOPOLOGY (disclinations from K_splay<0);
(v) self-interaction. CMB: framework has a real mass component (still 3rd-peak-open, but structurally unlike MOND).
WHY the explored levers failed the cluster: collective strain SCREENS, density-feedback SCREENS, wake DOESN'T
accumulate -- ALL are the baryon-TRACKING channel, where framework ~= MOND, so they inherit MOND's deficit.
THE MOND-IMPOSSIBLE CHANNEL (uncomputed): TOPOLOGICALLY-TRAPPED DISCLINATIONS from turbulent assembly.
  K_splay<0 SPONTANEOUSLY forms disclinations (instability, not response). Quiescent galaxy -> 1 disclination
  (the halo, RAR). CLUSTER built by violent hierarchical merging -> turbulent reordering TRAPS MANY disclinations
  (Kibble-Zurek: rapid symmetry-breaking during mergers freezes defects). They are TOPOLOGICAL (persistent,
  c_dir->0), GRAVITATE (Frank energy=mass), and are NOT baryon-tracking (medium defects seeded by DYNAMICS).
  => extra cluster mass = trapped-defect population, scaling with MERGER/ASSEMBLY HISTORY, not instantaneous
  baryon g_N. MOND (local force law) cannot encode assembly history; LambdaCDM halo is set by collapse not recent
  merger state. SAME object as the recombination-CDM primordial-defect lead (baryon-independent coherent
  deformation; seeded early by ordering, late by mergers) = the CDM ontology, two epochs.
TESTABLE DISCRIMINATOR: DISTURBED/recently-merged clusters should show LARGER mass excess than RELAXED clusters
  of equal baryon content. MOND: no difference (boost tracks baryons). LambdaCDM: excess tracks virialized halo,
  not recent merger state. Framework: excess tracks disturbance. -> relaxed-vs-disturbed cluster comparison
  (excess mass vs dynamical-disturbance indicator) SEPARATES the three theories. Framework-distinct, falsifiable.
STATUS: (A) distinctions = ESTABLISHED framework features. (B) defect recovery = SPECULATIVE (mechanism sound +
  framework-specific; magnitude [defect density x Frank energy =? ~2-3x deficit] UNCOMPUTED; the ONE lever left
  after the MOND-like ones screened). (C) prediction = genuine discriminator IF (B) holds.
REFRAME: cluster problem = "baryon-tracking channel shares MOND's deficit; the TOPOLOGICAL channel is where the
  framework departs, and it makes a falsifiable disturbed-cluster prediction." NEXT (if pursued): Kibble-Zurek
  defect-density estimate for cluster assembly x Frank energy per defect -> does it reach the ~2-3x deficit?

## KIBBLE-ZUREK CLUSTER-DEFECT ESTIMATE (2026-09-02, cluster_kibble_zurek.py). LIVE + physically-scaled.
Trapped-disclination mass from cluster assembly. K_splay<0 spontaneously makes disclinations; c_dir->0 -> merger-
trapped defects never heal -> persistent network. KZ: ~1 line per xi_turb^2 (xi_turb = merger-reordering coherence).
SCALING (analytic, K cancels): rho_def ~ K/xi_turb^2 ; M_def(R_cl) ~ K R_cl^3/xi_turb^2 ; M_RAR ~ K R_cl (isothermal)
  => M_def/M_RAR ~ (R_cl/xi_turb)^2.
NUMBERS: deficit needs M_def/M_RAR ~ 4.6 (cluster ~28x stellar; RAR ~5x). => required xi_turb ~ 0.47 R_cl.
  xi/R=1.0->1x ; 0.7->2x ; 0.5->4x (CLOSES) ; 0.3->11x ; 0.1->100x (over).
KEY: cluster mergers ARE subcluster-scale (subclusters ~0.3-0.7 R_cl), so the REQUIRED xi_turb ~0.5 R_cl MATCHES
the actual merger scale -- NOT tuned. => RIGHT ORDER OF MAGNITUDE + PHYSICALLY-SCALED. FIRST lever that is LIVE
(all MOND-like channels -- collective, density-feedback, wake -- SCREENED/failed; this topological one survives).
CAVEATS (uncomputed): (i) strongly SCALE-SENSITIVE -- fine turbulence (xi<<R) OVER-produces ((R/xi)^2 blows up),
so needs xi_turb NOT much below ~0.5 R_cl; (ii) defect ANNIHILATION (opposite-charge disclinations heal; c_dir->0
slows it but like/unlike mix sets the NET density) -- the real unknown; (iii) K->a0 normalization. Promising, open.
DISCRIMINATOR (framework-distinct, falsifiable): DISTURBED/recently-merged clusters -> MORE mass excess than
RELAXED clusters of equal baryon content (excess tracks merger state, not baryon g_N [MOND: none] or virial halo
[LambdaCDM: tracks halo]). A relaxed-vs-disturbed cluster comparison separates the three theories.
NET (cluster deficit): the baryon-tracking channel shares MOND's deficit (screens); the TOPOLOGICAL defect channel
is LIVE, framework-specific, right-order, physically-scaled -- the honest candidate recovery, with a falsifiable
disturbed-cluster test. Next precision step (if pursued): net defect density after annihilation (a defect-network
scaling-solution calc) x Frank energy, with the K=a0/G normalization -> pin M_def/M_RAR for realistic mergers.

## DEFECT-NETWORK SCALING SOLUTION (2026-09-02, cluster_defect_scaling.py). PINNED: recovers deficit, no free norm.
Refines the KZ estimate by adding ANNIHILATION. Driven-dissipative scaling: a defect network stirred at injection
scale xi_turb with pair-annihilation self-regulates to xi_net ~ xi_turb (creation at xi_turb balances annihilation;
density neither runs away nor drains). => ANNIHILATION RESOLVES THE RUNAWAY: (R/xi)^2 is bounded by the merger
scale, NOT a fine-turbulence blowup. c_dir->0 modification: annihilation is ADVECTION-only (merger flow), so
regulation happens DURING mergers then FREEZES (no coarsening after) -> current density set by the last major merger.
RESULT: M_def/M_RAR ~ (R_cl/xi_turb)^2, xi_turb = subcluster merger scale.
  xi/R = 0.7->2.0 ; 0.6->2.8 ; 0.5->4.0 (CLOSES) ; 0.4->6.2 ; 0.3->11 (over).
  Plausible merger bracket (xi/R=0.3-0.7): M_def/M_RAR in [2,11] -> BRACKETS the needed 4.6x; recovers at
  xi_turb ~ 0.47 R_cl (= subcluster scale, PHYSICAL not tuned).
ABSOLUTE NORM (K=a0/G): K cancels in the ratio; M_RAR = the a0-calibrated deep-MOND value = known ~5x stellar;
  M_def = ratio x M_RAR -> at 4.6x: M_def ~23x stellar, + RAR 5x = 28x stellar = observed cluster DM. MATCHES,
  NO FREE NORMALIZATION (anchored to a0/G).
OPEN: the c_dir->0 advection-only annihilation is NOT textbook (self-moving) scaling -> a real defect-network
  simulation WITH advection is needed to confirm xi_net~xi_turb and fix the O(1) coefficient. [scaling estimate,
  not a rigorous network solution.]
DISCRIMINATOR (sharpened): frozen network -> current density = last-major-merger xi_turb -> VIOLENT/recent mergers
  (finer xi_turb) => MORE excess; relaxed clusters keep the fossil network (c_dir->0, no coarsening). Falsifiable:
  excess vs merger-disturbance indicator. MOND: none; LambdaCDM: excess tracks virial halo not merger state.
=> RESTING POINT for the cluster thread: the topological-defect channel RECOVERS the cluster deficit at the
   physical merger scale, annihilation-regulated, with NO free normalization (a0/G anchor) -- a framework-specific,
   falsifiable answer. Remaining precision: a defect-network sim (advection + c_dir->0) to fix the O(1).

## ============ CONTINUATION SUMMARY (post-compaction state, 2026-09-02) ============
WHERE THE DM INVESTIGATION STANDS (one ontology: CDM = coherent director DEFORMATION without the baryon population):
  - GALAXIES/DWARFS/PLANETS: SUCCESS. density-feedback (S_eff=sin^4θ/[φ^6(1+βρ)]) screens high-ρ -> RAR, right
    direction across ~30 dex (solar-system no-DM, dwarf max-DM), cored profiles. Polaron effective mass
    m_eff=(χ/3)∫|∇n|^2 (POSITIVE, sign paradox resolved) -> M∝r. In paper: P7:rem:polaron, P7:rem:cdm_ontology.
  - CLUSTERS: baryon-tracking channel shares MOND's deficit (~4.6x short; collective strain SCREENS, wake doesn't
    accumulate, density-feedback screens -- all the MOND-like channels FAILED). The MOND-IMPOSSIBLE channel that
    WORKS: merger-trapped DISCLINATION network (K_splay<0 spontaneously makes disclinations; cluster mergers trap
    them; c_dir->0 -> frozen, persist). Defect-network scaling (annihilation self-regulates ξ_net~ξ_turb) gives
    M_def/M_RAR ~ (R_cl/ξ_turb)^2 ~ 4.6x at ξ_turb~0.47 R_cl = the PHYSICAL subcluster-merger scale, NO free
    normalization (K=a0/G anchors M_RAR). RECOVERS the deficit. Falsifiable discriminator: disturbed/recently-
    merged clusters -> MORE excess (MOND: none; LambdaCDM: tracks halo not merger state).
  - RECOMBINATION CDM: input (option 3). Open lead = primordial COHERENT deformation (same defect ontology, seeded
    by early ordering not mergers). See fabric_cmb_hypothesis_and_plan.md.
KEY SCRIPTS (papers/src_paper7/dm_polarization/): phase2c_A1_disclination.py (M∝r config + collective screening),
  frank_chi.py (χ=601, anisotropy tentative), frank_run.py (K_splay=-8.45, K_bend=+0.534), phase0_strain_energy.py,
  cluster_wake_calc.py, cluster_kibble_zurek.py, cluster_defect_scaling.py.
IMMEDIATE NEXT TASK (this turn): defect_network_sim.py -- 2D complex-field (CGL/XY) disclination dynamics with a
  merger-like STIRRING FLOW (advection) in the c_dir->0 (slow-relaxation / high-Peclet) limit. GOAL: measure the
  steady/frozen net defect spacing ξ_net vs the stirring scale ξ_turb -> CONFIRM the scaling ξ_net~ξ_turb and FIX
  the O(1) coefficient -> pin M_def/M_RAR=(R_cl/ξ_net)^2. Setup + smoke here; user runs the full parameter scan
  (ξ_turb, Peclet) on their machine. Uses base-env python (/usr/local/Caskroom/miniconda/base/bin/python, has numpy).
OPEN Q_H LABELLING (deferred): notes/paper_todo_qh_labelling.md -- Paper VII "Q=2"=vacuum vs ladder Q_H=2=lepton;
  don't write specific Q_H into Papers I/VII until the vacuum-winding question is settled with the user.

## DEFECT-NETWORK SIM BUILT + SMOKE-TESTED (2026-09-02, defect_network_sim.py).
2D complex-field (CGL/XY) disclination dynamics: d_t psi = D lap psi + (1-|psi|^2)psi/tau - (v.grad)psi, with
incompressible merger-like STIRRING v at scale xi_turb (refreshed each eddy-turnover), c_dir->0 = small D
(annihilation slow, advection-dominated; Peclet=amp*xi_turb/D). Counts vortices -> xi_net=L/sqrt(N_v) ->
SCALING COEFFICIENT c=xi_net/xi_turb. Cluster: M_def/M_RAR=(R_cl/(c*xi_turb))^2 at xi_turb=0.47 R_cl.
SMOKE (N=48, xi=8, D=0.15, Peclet=32, 400 steps -- POOR STATS, N_v still drifting 58->32, NOT converged):
  N_v~40, xi_net=7.63, c = xi_net/xi_turb = 0.95 (~1 -> scaling assumption xi_net~xi_turb CONFIRMED in smoke),
  M_def/M_RAR = 4.98 -> RECOVERS the ~4.6x deficit. [PRELIMINARY SMOKE, not the pinned value.]
=> encouraging: the scaling xi_net~xi_turb (c~1) holds even at smoke resolution, and c~1 recovers the deficit.
FULL RUN (user, on base-env python w/ numpy): python defect_network_sim.py --N 256 --xi 16 --steps 20000
  (well-resolved ~256 defects); then SCAN Peclet (--D 0.05/0.15/0.3 = c_dir->0 is small D/high Peclet) and
  xi (8/16/32) to check c is SCALE-INDEPENDENT and get c(Peclet). c(Peclet) is the pinned O(1) -> M_def/M_RAR.
  If c stays ~1 in the c_dir->0 (high-Peclet) limit -> the topological channel RECOVERS the cluster deficit,
  pinned, no free normalization. If c drifts far from 1 -> revise.

## DEFECT-SIM FULL RUN + DYNAMICS CATCH (2026-09-02). Overdamped=WRONG dynamics; c=2.77 is an ARTIFACT.
User ran defect_network_sim.py --N 256 --xi 16 --steps 20000 (D=0.15, Peclet=64): N_v~33, xi_net=44.4,
c=xi_net/xi_turb=2.77, M_def/M_RAR=0.59 -> does NOT recover (and N_v still drifting down 46->30 = still coarsening).
CATCH: the sim uses OVERDAMPED dynamics (d_t psi = D lap psi + ...). But the framework's director is INERTIAL:
chi d_t^2 n = K grad^2 n, c_dir^2=K/chi, c_dir->0 = HEAVY (large chi). Overdamped: defects have mobility ~D, move
together, ANNIHILATE, network COARSENS (xi_net ~ sqrt(D t) grows -> sparse -> c large -> fails). Inertial-heavy:
defects have MASS ~chi, move at ~c_dir->0, OSCILLATE not overdamped-relax -> DON'T self-annihilate -> FROZEN ->
xi_net~xi_turb -> c~1 -> recovers. So c=2.77 is an OVERDAMPED ARTIFACT; wrong dynamics for a c_dir->0 director.
(The earlier smoke c~0.95 was under-relaxed = accidentally less-coarsened = closer to frozen, right-ish for the
wrong reason.)
FIX (in progress): rewrite the sim as a DAMPED WAVE: chi d_t^2 psi = K lap psi + (1-|psi|^2)psi/tau - gamma d_t psi
+ stirring, with c_dir=sqrt(K/chi) SMALL (large chi) and SMALL gamma. Then defects are heavy/frozen; measure c.
HONEST open: even inertial-heavy has SOME damping gamma -> over Gyr (many eddy-turnovers) slow annihilation may
still coarsen. Recovery requires annihilation slow ENOUGH over the cluster lifetime. The inertial run decides.
STATUS: the topological-recovery verdict is UNSETTLED pending the correct (inertial) sim -- do NOT quote c=2.77
as the framework's answer (wrong dynamics); do NOT quote c~1 as settled either (needs the inertial run).

## INERTIAL DYNAMICS IMPLEMENTED + RUN (2026-09-02). Freeze-vs-coarsen CONFIRMED; but frozen c~2, not ~1.
defect_network_sim.py now has --inertial (chi d_t^2 psi = c_dir^2 lap psi + (1-|psi|^2)psi/tau - gamma d_t psi
- adv), velocity-Verlet stepping, --cdir (default 0.1), --gamma (default 0.02), --no-stir (pure coarsening test).
IC: dense random phases smoothed ~xi. Overdamped path retained (default) for comparison; it is the WRONG dynamics.

NO-STIR (pure coarsening, N=96 xi=8, 40000 steps) -- the clean freeze-vs-coarsen discriminator:
  OVERDAMPED (D=0.4): N_v 766 -> 8 -> ... -> 0 by step 32k. Defects have mobility, coarsen, ANNIHILATE AWAY. c->inf.
  INERTIAL (c_dir=0.1,g=0.02): N_v 870 -> 32 (fast shedding of sub-xi IC noise) -> FROZEN at ~24-26 rest of run.
  => CONFIRMED: heavy (c_dir->0) defects freeze into a persistent fossil network; overdamped ones disappear.
     The dynamics catch was right -- overdamped coarsening was an artifact.
STIRRED (merger injection, N=128 xi=16, 30000 steps):
  OVERDAMPED (D=0.15): c=3.47, M_def/M_RAR=0.38 -- reproduces the ~2.77 artifact family (fails, keeps coarsening).
  INERTIAL (c_dir=0.1,g=0.02): N_v=5323, xi_net=1.75, c=0.11 -- GRID-SCALE NOISE ARTIFACT (opposite failure):
    c_dir->0 makes the elastic term (c_dir^2 lap = 0.01*lap) too weak to heal the windings advection injects, so
    the field shreds to ~1 defect per 3 cells. "M_def/M_RAR=377" is nonsense. NOT a physical network.
  => BOTH stirred numbers are artifacts in OPPOSITE directions (overdamped over-coarsens c~3.5; inertial
     under-heals c~0.1). The 2D-proxy sim CANNOT pin the stirred steady-state coefficient. The ad hoc
     advection-in-acceleration coupling is the culprit (see caveat below).

## VERDICT (2026-09-02): sim does NOT settle clusters; earlier "c~1 recovers 4.6x" is RETRACTED.
TRUSTWORTHY result (no-stir freeze test only): inertial c_dir->0 defects FREEZE into a persistent fossil network
(overdamped ones coarsen to zero). The KZ/freeze picture is CONFIRMED -- a real, framework-specific,
MOND-impossible feature (recent mergers -> more trapped defects -> more excess: the falsifiable discriminator).
BUT the frozen network sits at c~2 (~2x COARSER than the injection scale) -> M_def/M_RAR ~ 1 -> an O(1) boost on
top of the RAR, NOT the ~4.6x needed to close the cluster deficit.
The stirred steady-state magnitude is NOT reliably computable with this 2D proxy (c swings 0.1<->3.5 with
dynamics/coupling). So the coefficient is genuinely UNPINNED by this tool.
=> HONEST STATUS: the topological fossil-network channel is REAL and persistent, but the simulation does NOT
   support the deficit-closing claim. Best trustworthy estimate = an O(1) partial boost (helps, likely does NOT
   fully close the cluster deficit). The earlier cluster_defect_scaling.py "brackets 4.6x, recovers at
   xi_turb~0.47 R_cl" ASSUMED c~1; the sim says c is not 1 in any trustworthy regime -> that optimism is RETRACTED.
   Pinning it properly needs a genuine 3D inertial director sim with RIGOROUS advection (material derivative
   chi D^2psi/Dt^2, not -adv in the acceleration) -- beyond this 2D proxy. Parked there, honestly.

HONEST NUMBER: even in the CORRECT inertial dynamics, the frozen network sits at c~2 (no-stir plateau: N_v~24 in
96^2 seeded at xi=8 -> xi_net~19.6 -> c~2.4), NOT c~1. There is a fast initial annihilation BURST (the heavy
defects still let very-close opposite pairs annihilate via the core term) before the freeze. c~2 gives
M_def/M_RAR = 1/(c*0.47)^2 ~ 1.1 -- a MINOR boost (~1x RAR), NOT the ~4.6x needed to close the cluster deficit.
=> TRENDING NEGATIVE for full recovery: the topological fossil network is real and persistent (freeze confirmed),
   but freezes ~2x too COARSE to supply the cluster deficit on its own. The stirred-inertial steady state (pending)
   is the last lever -- continuous merger injection may hold it denser than free decay. If stirred c also ~2 ->
   the topological channel gives only an O(1) boost, not the deficit; honest verdict = does NOT close clusters.
CAVEAT on the sim (state honestly): advection is added to the ACCELERATION (-adv in accel), which is an ad hoc
coupling -- physically advection is first-order transport (material derivative), not a force. So the stirred
inertial number is a PROXY, not rigorous. The no-stir freeze result is clean; the stirred magnitude is soft.

## 3D INERTIAL SIM (2026-09-02, defect_network_sim_3d.py). c_dir->0 => c->1 (freezes at injection scale).
Built the 3D version (disclination LINES, not points) with the advection bug FIXED: stirring is modelled as a
Kibble-Zurek RE-QUENCH (rapidly reorient at xi_turb, then freeze), not the ad hoc advection-in-acceleration.
IC/quench = FFT-Gaussian-correlated unit field at xi_turb (clean, no sub-xi grid noise -> no spurious initial
annihilation burst, unlike the 2D roll-smoothed IC). Measure = total defect-LINE length via plaquette windings
in all 3 planes -> xi_net = sqrt(V/P). REPORT c = xi_net/xi_net0 (frozen/imprinted) -- CALIBRATION-FREE (the
geometric O(1) cancels): how much coarser the frozen network is than the freshly-quenched xi_turb tangle.
Benchmarks (numpy, no scipy): 64^3 8.5min/20k steps, 96^3 14min, 128^3 47min.

c_dir SCAN (48^3, xi=8, gamma=0.02, 1500-step smoke), c=xi_net/xi_net0, P0=394:
  c_dir=0.3 : c=2.81 (P 394->50, coarsens -- waves propagate fast enough to find+annihilate)
  c_dir=0.1 : c=1.10 (P 394->324, nearly frozen)
  c_dir=0.03: c=1.03 (P 394->372, frozen)
  c_dir=0.01: c=1.005 (P 394->390, frozen)
=> CLEAR monotonic trend: as c_dir->0 (THE FRAMEWORK'S REGIME, chi diverges -> heavy director), the KZ tangle
   FREEZES at the injection scale, c->1. c=1 gives M_def/M_RAR = 1/(1*0.47)^2 = 4.5 -> RECOVERS the ~4.6x deficit.
   This REVERSES the 2D pessimism (c~2): the 2D c~2 was partly a roll-smoothed-IC grid-noise annihilation burst;
   the clean FFT IC + true c_dir->0 limit freezes at c~1.

HONEST CAVEAT (do not overstate): with ANY gamma>0 and c_dir>0 the true ground state is defect-FREE (P->0
eventually) -- "frozen" is a TIMESCALE statement, not a hard plateau. The physical content: coarsening RATE ->0
as c_dir->0, and the framework sits at c_dir->0 (chi~diverges), so the coarsening time >> cluster age -> c~1 TODAY.
DECISIVE TESTS RUNNING (background): (A) 48^3 long-time (20000-step) coarsening-RATE scan c_dir=0.15/0.08/0.04/0.02
-- does rate ->0 as c_dir->0? (B) 96^3 higher-res confirmation (c_dir=0.05, 15000 steps) -- is c~1 resolution-robust?
STATUS: TRENDING POSITIVE for topological cluster recovery in the CORRECT (c_dir->0 inertial) dynamics -- pending
the rate-scan + higher-res confirmation. The overdamped c=2.77 and the 2D c~2 are both superseded.

## RATE SCAN VERDICT (2026-09-02, 48^3, 20000 steps). NO plateau -- freezing is a TIMESCALE effect. CONDITIONAL.
Long-time (20k-step) coarsening scan, P0=394, P(t) history:
  c_dir=0.15: P->0 by ~5000 steps    (fully coarsens)
  c_dir=0.08: P->0 by ~11000 steps   (fully coarsens, slower)
  c_dir=0.04: P=116@20k, still falling -> 0 at ~29000 (extrapolated)
  c_dir=0.02: P=314@20k, still falling -> 0 at ~102000 (extrapolated)
=> DECISIVE: there is NO true plateau at any nonzero c_dir. P drains monotonically to 0 in every case -- the
   defect-FREE ground state always wins EVENTUALLY. c_dir controls ONLY THE RATE. The short-smoke "c->1" was
   "hasn't coarsened YET," not a frozen state. So "frozen" is STRICTLY a timescale statement (my caveat was right,
   the interim optimism was premature).
COARSENING TIME vs c_dir (t_zero in steps): 5k, 11k, 29k, 102k for c_dir=0.15,0.08,0.04,0.02.
  log-log slope: t_coarsen ~ 1/c_dir^1.5-1.7 (steepening toward small c_dir) -> DIVERGES as c_dir->0.

PHYSICAL MAPPING (order-of-magnitude, tag [N]/[O] -- NOT rigorously derived from the sim):
  Freezing over a cluster age requires the network to coarsen by << xi_turb in t_cluster. Ballistic defect speed
  ~ c_dir -> condition c_dir * t_cluster << xi_turb, i.e. c_dir << xi_turb/t_cluster.
  With xi_turb ~ 0.5 Mpc (subcluster merger scale), t_cluster ~ 10 Gyr: c_dir << ~50 km/s.
  KEY CONSISTENCY: this is the SAME "cold" condition the framework ALREADY asserts (§4: c_dir->0 = cold,
  pressureless, no Jeans obstruction -- what lets the director-DM cluster at all). So the topological-network
  freezing is NOT a new assumption: it is TIED to the coldness that gives the DM sector its clustering. The
  coldness that makes it dark-matter-like is the same c_dir->0 that freezes the fossil network.

HONEST VERDICT (tagged):
  [E-sim] In the correct inertial dynamics, at FREEZE-IN the KZ tangle sits at c=xi_net/xi_turb~1 (recovers the
    ~4.6x deficit) -- BUT only transiently; it coarsens away at rate ~c_dir^1.5-2.
  [N] Coarsening time diverges as c_dir->0; the framework HAS c_dir->0, so the fossil network is parametrically
    arrestable, and the freezing condition (c_dir << xi_turb/t_cluster ~ 50 km/s) COINCIDES with the coldness the
    framework already requires for director-DM clustering.
  [O] The ACTUAL numerical value of c_dir is a CUTOFF-SET magnitude the framework does NOT derive (CLAUDE.md
    §10). So recovery is CONDITIONAL: it holds IFF c_dir is below the (age/dynamical-time)-set threshold. Cannot
    be asserted as settled without pinning c_dir (or an independent bound that c_dir < ~50 km/s in clusters).
  [E] FALSIFIABLE DISCRIMINATOR stands and is now sharper: since the network coarsens (slowly), RELAXED old
    clusters have had longer to coarsen -> LESS excess; recently-DISTURBED/merged clusters -> freshly re-quenched,
    denser network -> MORE excess. MOND: none; LambdaCDM: excess tracks the halo, not the merger/relaxation age.
=> NET: the topological channel is a REAL, framework-specific, parametrically-viable cluster-recovery mechanism
   whose sufficiency hinges on one un-derived magnitude (c_dir), and which makes a distinct falsifiable prediction.
   Honest status = CONDITIONAL RECOVERY, not settled recovery. Supersedes both the c=2.77 "fails" and the naive
   "c~1 closes it." Awaiting 96^3 higher-res confirmation of the rate/robustness.

## 96^3 HIGHER-RES CONFIRMATION (2026-09-02). Resolution-robust; verdict stands.
96^3, xi_turb=12, c_dir=0.05, gamma=0.02, 15000 steps (wall 725s): P0=1072 (scales ~volume/xi^3 vs 394@48^3 ->
measurement machinery sound), P=738 after 15k steps, c=1.20 and STILL DECLINING (1072->738 monotonic, no plateau).
=> CONFIRMS the 48^3 rate scan at higher resolution: small c_dir -> slow-but-steady coarsening, c->1 only
   transiently, no true freeze. The CONDITIONAL-RECOVERY verdict is resolution-robust.

## HEALING-RATE REFINEMENT (2026-09-02, user). Survival is set by ANNIHILATION rate, not signal speed.
Key distinction (folded into P7:rem:defect_network, 3rd paragraph): the quantity that decides whether the
fossil network survives is the defect HEALING (annihilation) rate, which is DISTINCT from the director signal
speed c_dir. THREE rates now:
  (1) signal speed c_dir=sqrt(K/chi) -- how fast orientation disturbances propagate (->0).
  (2) healing rate -- how fast a trapped defect DISAPPEARS: a MOBILITY (line tension K driving annihilation vs
      rotational viscosity gamma_rot resisting + heavy inertia chi), NOT the wave speed. gamma_rot is the SAME
      viscosity the paper already flags as unpinned in P7:rem:polaron (elastic drag).
  (3) magnitude/Bell mode c_s=c/phi -- fast smooth/global reconfiguration + core-melting energy scale.
TOPOLOGICAL PROTECTION (the payoff): a topological defect CANNOT be removed by smooth reconfiguration however
fast -- only by annihilation (partner) or core-melting. So even if the shared-orientation (Bell) geometry makes
smooth reconfiguration GLOBAL and FAST, that mode carries no topological charge and heals smooth strain, NOT the
trapped network. => the "Bell makes reconfiguration fast" worry does NOT threaten the fossil network. Survival
needs the ANNIHILATION rate (tension/viscosity, topologically bottlenecked) slow over a cluster age -- met in the
same cold-director limit, contingent on the same unpinned quantities (gamma_rot, c_dir).
HONEST: the 3D sim's coarsening IS a healing rate (annihilation via the core term) at ONE fixed damping
(gamma=0.02); it found healing slow at small c_dir but never zero -- topology bottlenecks healing, does NOT forbid
it. Real gamma_rot unpinned -> healing rate unpinned, but it is the physically-correct survival criterion
(sharper than "c_dir vs assembly velocity").
