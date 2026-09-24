# STANDING HYPOTHESIS: the CMB peaks as the fabric coherence transition (+ modeling plan)
Consolidated so it survives context compaction. Speculative (route 2); NOT in the paper.

## THE HYPOTHESIS (user's, refined)
The CMB acoustic peaks are the SIGNATURE of the fabric's cyclical coherence transition through
EXTENDED recombination -- NOT the imprint of a pre-existing external CDM. Recombination is not instant
(last-scattering has real thickness, dz~100+); the fabric clusters where it can, is torn back by
radiation pressure, clusters again -- and that cyclical clustering-vs-tearing IS the acoustic oscillation.
Same KIND of process as an electron reconfiguring under a photon (Paper XIII): a geometric reconfiguration
in progress, cyclical. => "no separate CDM"; Omega_DM~0.26 that LambdaCDM infers is an artifact of assuming
standard GR + a permanent cold component.

## WHY IT'S NOT CRAZY
The cyclical clustering-vs-radiation-pressure IS the standard acoustic-oscillation mechanism; the hypothesis
only re-attributes the GRAVITY SOURCE (fabric coherence transition, density-feedback) instead of external CDM.
Peak POSITIONS (harmonic ladder l~220,440,660) come from the sound horizon -> free either way.

## THE DECISIVE TEST: the THIRD PEAK
Peak HEIGHTS, esp. the third peak, are the direct fingerprint of the CDM density. Observed: HIGH third peak
=> the driving fluid fell into DEEP, PERSISTENT wells during the oscillations. The hypothesis must produce a
high third peak from TRANSIENT, partial clustering (no permanently-clustered component).
MECHANISM LEAN (against): (i) the fabric is radiation-like before ordering (Gate 1, u_disordered~rho_cond~1/a^4);
(ii) the density-feedback SCREENS in overdense regions (omega_BD^eff=omega_BD(1+beta rho) grows) -> weaker
driving where you need it. Both push the third peak LOWER than even no-CDM standard gravity. Double strike.
HONEST LIMIT: the third-peak height from transient-clustering + screened gravity is a BOLTZMANN-level calc;
cannot be settled by scaling/mechanism argument alone. Genuine open hypothesis, not a proven failure.

## KEY THREAD FINDINGS (context, so they survive)
- Ordering epoch CLOSED (ordering_epoch.py): T_ord ~ Lam_cond ~ T_CMB(today) ~ 1e-4-1e-3 eV, ~200-2600x below
  T_rec~0.26 eV. Medium DISORDERED at recombination; orders at z~O(1)-few. Director-strain DM is LATE-time
  galactic only. Cutoff-independent (set by Lam_cond from T_CMB).
- Phase 0 (fabric under expansion): F2 (stretch); framework fixes it via rho_CMB=10 Lam_cond^4 => Lam_cond~1/a,
  xi_cond~a, rho_cond~1/a^4. 2I PRESERVED under isotropic stretch (scale grows, symmetry intact).
- Gate 1 EOS (gate1_eos.py): NO fabric sub-component is matter-like (1/a^3) at recombination. Thermal condensate
  ~1/a^4 (radiation); frozen vacuum ~const (DE); disordered director stress ~1/a^4 (radiation); galactic
  disclination strain ~1/a^3 (matter) but LATE. => option-4 route-1 ("fabric bulk=CDM") FAILS.
- Route 2 (modified gravity, no CDM): faces the peak-alternation/third-peak test + the screening obstacle.
  This hypothesis is its strongest, testable form.
- Fallback = OPTION 3: recombination Omega_DM is a genuine input the framework does not yet derive.

## MODELING PLAN (start CHEAP; the big run is last)
PREREQUISITE -- Phase A: framework LINEAR PERTURBATION THEORY (does not exist yet; load-bearing).
  Derive from the density-feedback action: (A1) modified Poisson/Einstein for delta_rho -> potential Phi with
  the density-feedback + chameleon screening; (A2) the fabric coherence-transition dynamics (how ordering turns
  on through extended recombination, the transient clustering); (A3) the perturbation EOS of each component.
  This is the real intellectual work; everything downstream needs it.
LEVEL 0 -- analytic peak estimate (CHEAP, DO FIRST, may already decide it):
  Hu-Sugiyama semi-analytic acoustic-peak formulas: peak heights vs driving (well depth) + baryon loading +
  radiation driving. Plug in the framework's EFFECTIVE driving (transient fabric clustering, screened) from
  Phase A. Estimate the 3rd-peak height; compare to observed. GATE: if driving is far too weak (screening),
  3rd peak too low -> route 2 likely fails (confirms lean). If in range -> escalate. NO large run.
LEVEL 1 -- parametric modified Boltzmann (MEDIUM): CLASS/CAMB with Omega_cdm=0 + an EFFECTIVE fabric component
  (parametrised time-dependent clustering/EOS that turns on at recombination). Scan; can ANY config reproduce
  the observed peaks (esp. 3rd)? Tests viability with existing tools, no code-surgery.
LEVEL 2 -- full framework Boltzmann (LARGE; only if L0/L1 warrant): implement Phase-A perturbation eqs in
  CLASS source; compute C_l; fit Planck. Months of work; needs Phase A complete.
COST/RISK: Phase A = derivation (cheap compute, hard physics). L0 = cheap + fairly decisive. L1 = existing
  tools. L2 = large run the user worried about -- DEFER until L0/L1 justify. Do Phase A + L0 FIRST.

## PHASE A FIRST PASS -- DONE (phaseA_poisson.py). ROUTE 2 CLOSED (negative). No large run needed.
Modified Poisson: grad^2 Phi = 4 pi G_eff delta_rho, G_eff = G_N[1 + 2 beta_eff^2], beta_eff = beta/(1+beta rho).
TWO independent, decisive failures of route 2:
(1) MAGNITUDE too small by ~89x. Fifth-force boost G_eff/G_N = 1 + 2(beta/phi^2)^2 ~ 1.06 (6%, unscreened) to
    1.0 (0%, screened). To replace CDM, baryon gravity must mimic the full matter budget: G_eff/G_N ~ Omega_m/
    Omega_b ~ 6.3 (533%). The density-feedback fifth force (beta~0.45, /phi^4 suppressed) is ~2 orders too weak.
(2) STRUCTURAL: even a large boost on baryons cannot provide CDM's NON-OSCILLATING wells -- boosted baryons
    still oscillate (radiation-locked). CDM's role is the persistent well, which no baryon-boost supplies.
Plus TIMING (independent): the coherence transition is at z~few (T_ord~Lam_cond<<T_rec), NOT recombination.
=> ROUTE 2 FAILS. OPTION 3 is the closure: recombination Omega_DM~0.26 is a GENUINE INPUT the framework does
   not derive; its modified gravity is far too weak (and structurally wrong) to eliminate CDM.

## FINAL CLOSURE OF THE WHOLE DM/CDM THREAD
- Late-time galactic DM: DERIVED (director strain, orders z~few). Firm.
- Dark energy: DERIVED (frozen vacuum, Lambda_obs). Firm.
- Recombination CDM (Omega_DM~0.26): NOT derived -- a genuine input/gap. Route 1 (fabric-stress=CDM) failed at
  Gate 1 (fabric is radiation-like, no matter-like component at recomb). Route 2 (modified gravity, no CDM)
  failed at Phase A (fifth force ~89x too weak + can't make non-oscillating wells).
- The paper's dm_status/dm_ic honestly flags this (Omega_DM an input). The ordering-epoch sharpening (paper_todo)
  can be added: the director strain is definitively LATE (z~few), so it is NOT the recombination CDM.
STATUS: dark-sector cosmology thread CLOSED at option 3. The Boltzmann run (Level 1/2) is NOT needed -- the
linear enhancement estimate already excludes route 2 by ~2 orders of magnitude.

## UPDATE 2026-08-30 -- BELL-MECHANISM REFRAME: the "CLOSED" above is for the LOCAL route only; a distinct NONLOCAL route is OPEN
Verified in repo (Paper VIII/IX), not guessed. This REOPENS a route the phaseA closure did not test.

GROUNDING (established in the framework already):
- Paper VIII (main_paper8.tex:334, :478, :633) ESTABLISHES the no-signaling GEOMETRIC mechanism: the condensate
  orientation n_hat is a SINGLE SHARED field quantity, a boundary condition "fixed at pair creation"; the
  correlation arises because both sites probe the SAME geometry, NOT from FTL communication; measurement
  back-action propagates at c_s = c/phi; the currency is literally "the geometric cost of the H->C" reduction
  (:478). Realistic + deterministic at field level, no signaling.
- Paper IX "post-quantum window" (P9:rem:nosignaling) = the untested LEPTON Bell CHSH ceiling (no experiment to
  date), NOT a causality hedge. The no-signaling mechanism itself is established for the quantum regime.
- So "nonlocal correlation, c-limited visibility, no signaling" is the framework's OWN mechanism, resting on its
  FOUNDING premise: ONE CONNECTED 2I CONDENSATE FROM THE START (the very reason for the hopfion). The cosmology
  INHERITS this premise; it is not a new bolt-on.

THE REFRAME (user's): CMB structure is globally correlated because it is ONE connected geometric configuration
sharing a boundary condition (Bell mechanism at cosmic scale); the photons are the c-limited "smoke-ring"
readout -- the geometric cost of the knots reconfiguring to match the new geometry. CDM is then the BOOKKEEPING
ARTIFACT of forcing that globally-correlated readout through FIXED, LOCAL gravity. Bonus: same answer-shape as
Bell also addresses the HORIZON PROBLEM (super-horizon correlations without inflation-as-causal-contact).

WHY phaseA_poisson.py DOES NOT CLOSE THIS:
- phaseA tested the LOCAL FIFTH-FORCE route: G_eff = G_N[1+2 beta_eff^2] boosting BARYON gravity. Failed twice
  (magnitude ~89x too weak; STRUCTURALLY boosted baryons still OSCILLATE, radiation-locked -> no non-oscillating
  wells).
- The Bell/global-geometry route is DISTINCT: the non-oscillating "wells" are the GLOBAL GEOMETRIC TEMPLATE (a
  boundary condition, like n_hat fixed at the source), NOT oscillating baryons => ESCAPES the structural prong.
  It does not rely on a local fifth force => the local-Poisson magnitude bound may not apply.

THE PIVOT (the real open question -- checkable BEFORE any Boltzmann run):
  Is the framework's gravity FUNDAMENTALLY LOCAL (modified Poisson, grad^2 Phi = 4 pi G_eff delta_rho_local),
  or does the density-feedback / geometric-cost mechanism carry a genuinely NONLOCAL source term (potential
  sourced by the GLOBAL geometric configuration, not local stress-energy)?
   - LOCAL   => phaseA's closure STANDS; the Bell reframe does not rescue it; option 3 (Omega_DM input) holds.
   - NONLOCAL => phaseA does not capture it; this route (call it 2') is genuinely OPEN, needs new perturbation th.
  phaseA ASSUMED local modified-Poisson. The Bell mechanism (n_hat as shared boundary condition IS a nonlocal
  correlation) SUGGESTS the nonlocal term exists. This assumption was NEVER checked against Paper VII's actual
  gravitational sector. NEXT STEP: determine whether Paper VII induced/density-feedback gravity has a nonlocal
  (global-geometric) potential source, or is strictly local. That decides local-Poisson closure vs open.

THE CALCULATION STILL NEEDED (the user's named target; unchanged in form):
  Even if nonlocal, the global geometric template must reproduce the C_l from the geometric-cost dynamics:
  (i)  PEAK POSITIONS -- the sound-horizon / harmonic-ladder scale (l~220,440,660): a boundary-condition template
       is NOT automatically this DYNAMICAL scale; the geometric-cost dynamics must GENERATE it.
  (ii) THIRD-PEAK HEIGHT -- the non-oscillating template must have the right depth/shape.
  Requires the (still nonexistent) Phase A framework linear perturbation theory, but now with the
  GLOBAL-GEOMETRIC-CORRELATION source, not merely a local fifth force.

REVISED STATUS (supersedes "CLOSED at option 3" for the nonlocal route):
  - Route 1 (fabric bulk/stress = CDM): CLOSED negative (Gate 1 -- radiation-like, no matter comp. at recomb).
  - Route 2 LOCAL (fifth-force-boosted baryons): CLOSED negative (phaseA -- 89x weak + oscillating).
  - Route 2' NONLOCAL (Bell/global-geometric template, no local CDM): OPEN. Hinges on the LOCAL-vs-NONLOCAL
    gravity pivot; NOT excluded by existing calcs.
  - Option 3 (Omega_DM a genuine input) remains the honest FALLBACK if the pivot resolves LOCAL.
  PAPER STANCE UNCHANGED (Omega_DM input honestly flagged). Do NOT put route 2' in the paper -- speculative.

## UPDATE 2026-08-30 (2) -- PIVOT RESOLVES *LOCAL*; but that opens ROUTE 1' (GEOMETRIC CDM), which Gate 1 did NOT test
User's intuition: gravity IS fundamentally local (local geometric configuration holds the density). This does NOT
force option 3 -- it opens a THIRD route distinct from route 1 (killed) and route 2' (nonlocal, now deprioritised).

ROUTE 1' -- CDM = a COLD, PRESSURELESS, ANISOTROPIC geometric configuration of the vacuum substrate.
Not the isotropic thermal bulk (Gate 1 killed that: relativistic -> 1/a^4). Instead a population of local geometric
textures (director-polarisation configurations; "LCD-polarisation" analogy) that are:
  - PRESSURE-DECOUPLED from radiation by the SEMI-DIRAC anisotropy (easy-radial / hard-perpendicular). Radiation
    pressure acts along the easy-radial axis; the configuration sits in the hard-perpendicular axis -> photons do
    not push it -> it does NOT oscillate. This is CDM's defining "no radiation pressure" property, geometrically.
  - NON-RELATIVISTIC along the hard-perpendicular axis: the semi-Dirac band is QUADRATIC (massive, E~p_perp^2/2m*)
    there, LINEAR (relativistic) along easy-radial. A mode carried on the quadratic axis is COLD/massive -> w~0,
    matter-like (1/a^3). Gate 1 averaged ISOTROPICALLY and saw only the relativistic (linear) behaviour; it never
    tested the anisotropic quadratic mode. => GENUINE ESCAPE HATCH from Gate 1.
  - LOCALLY GRAVITATING + PERSISTENT: unaffected by radiation pressure, these wells do not oscillate or wash out;
    they deepen by density feedback (the "knots" / debris-caught-in-eddies picture). Standard CDM cosmology is then
    REPRODUCED, not overturned: persistent wells, matter-radiation equality ~z3400, baryons free-fall in at
    recombination (knots stabilise, photons decouple along easy-radial -> CMB). Framework inherits the CMB fit.

THE CRUX NUMBER (cheap, decisive, DO FIRST): the EFFECTIVE MASS m* of the hard-perpendicular (quadratic) band.
  COLD (viable CDM) iff m* >> T_rec ~ 0.26 eV (so the mode is non-relativistic at recombination and earlier).
  If m* ~ Lam_cond ~ 1e-4 eV -> still relativistic -> fails exactly like route 1.
  m* = curvature of the quadratic band; estimate from the framework's semi-Dirac dispersion (find where defined --
  DM papers / Paper VII excitation spectrum). ONE number largely gates route 1' before ANY full model.

DOWNSTREAM (only if m* passes the gate):
  (a) ABUNDANCE: number density x m* must give Omega ~ 0.26 at recombination and equality at ~z3400.
  (b) SEEDING: primordial fluctuations from the framework's OWN inflation (Paper VII flux/tower).
  (c) CONTINUITY vs late-time strain: is this cold config the PROTO-form of the galactic disclination strain
      (one component, cold from early), or a DISTINCT anisotropic mode? The ordering epoch (T_ord~Lam_cond, z~few)
      was the NEMATIC coherence transition -- possibly a DIFFERENT DOF from this cold massive mode. RECONCILE:
      does route 1' imply TWO director-sector components (early cold anisotropic CDM + late nematic strain), and is
      that consistent? Open.

COST HONESTY: user is right that a FULL CMB reconstruction needs the whole model (geometric-stress temperature +
local density + free-flowing baryon/photon + vacuum sim). But the GO/NO-GO does NOT: the m* estimate is cheap and
decisive. Do the cheap gate FIRST; only invest in the full model if m* >> T_rec.

REVISED ROUTE MAP:
  - Route 1  (isotropic fabric bulk = CDM): CLOSED (Gate 1, relativistic).
  - Route 1' (ANISOTROPIC cold geometric config = CDM, LOCAL gravity): OPEN, most promising + most conservative.
    Gate = effective mass m* vs T_rec. NOT tested by Gate 1 (isotropic average missed the quadratic mode).
  - Route 2  (local fifth-force baryons): CLOSED (phaseA).
  - Route 2' (nonlocal Bell template): deprioritised (pivot resolves LOCAL per user intuition), not formally excluded.
  - Option 3 (Omega_DM input): fallback iff route 1' fails the m* gate.

## UPDATE 2026-08-30 (3) -- ROUTE 1' GATE RUN (route1prime_mass_gate.py). VERDICT: FAILS at dispersion level.
Read Paper VII eq P7:eq:semidirac directly: E^2 = A^2 k^4 sin^4(theta) + v^2 k^2 cos^2(theta).
  - HEAVY (perp, theta=pi/2) band: E = A k^2 -> GAPLESS. No rest-mass term. m* = 0. So the "m* >> T_rec cold
    massive mode" the gate looked for DOES NOT EXIST at the band-structure level. "Heavy/cold" here = vanishing
    GROUP VELOCITY (and, for the director, vanishing SOUND SPEED c_dir=sqrt(K/chi)->0, P7:eq:cdir), NOT a gap.
  - The only genuinely cold object is the DIRECTOR soft mode, mass = ultralight m_xi ~ H0 (DE sector).
CLUSTERING-CDM requires mass in the window: m >> T_rec (non-relativistic) AND m >> H(z_rec) (awake, not frozen).
Numbers (T_rec=0.259 eV, H_rec~1.6e-29 eV):
  - heavy-band gap  m=0            -> nothing.
  - director m_xi ~ H0 =1.4e-33 eV -> m/T_rec~5e-33 (relativistic-irrelevant) AND m/H_rec~9e-5 (HUBBLE-FROZEN,
    w=-1, DE-like). Cold in c_dir but FROZEN at recombination: does NOT cluster (fuzzy-DM onset is z~0, not z~1100).
  - condensate Lam_cond=1.2e-4 eV  -> m/H_rec~7e24 (awake) BUT m/T_rec~5e-4 (RELATIVISTIC) -> radiation-like (Gate1).
  NOTHING in the sector has a mass anywhere near the required >~0.26 eV. The two available scales straddle it
  uselessly: the director is cold-but-frozen; the condensate is awake-but-hot. No cold+awake mode exists.
=> ROUTE 1' FAILS. The semi-Dirac anisotropy genuinely delivers the LATE-time galactic DM (splay instability,
   flat curves, radiation-pressure decoupling) but is GAPLESS, so it supplies NO early cold-matter component at
   recombination. This is the sharpest statement of the gap: read off Paper VII's own dispersion, not a scaling.

## FINAL VERDICT OF THE WHOLE THREAD (all routes exhausted)
  Route 1  (isotropic fabric bulk = CDM):        CLOSED (Gate 1, radiation-like).
  Route 1' (anisotropic cold geometric = CDM):   CLOSED (this gate -- heavy band gapless; no m* >> T_rec mode).
  Route 2  (local fifth-force baryons):          CLOSED (phaseA, ~89x weak + oscillating).
  Route 2' (nonlocal Bell template):             OPEN in principle, but pivot resolves LOCAL (user) -> not pursued.
  => OPTION 3 IS THE HONEST STATUS: recombination Omega_DM~0.26 is a GENUINE INPUT the framework does not derive.
     The DM sector explains LATE-time galactic DM (director strain, z~few) and DE (frozen vacuum), NOT the
     recombination cold-matter budget. Paper stance (Omega_DM flagged as input) is correct and unchanged.
     Nothing speculative goes in the paper.

## UPDATE 2026-08-30 (4) -- METHOD CORRECTION (user): mass is EMERGENT from geometry, not a precondition.
The route-1' "mass gate" asked the WRONG question. Reading a gap off the LINEAR dispersion tests small
FLUCTUATIONS (phonons/Goldstone), not the gravitating energy of a large-amplitude COHERENT configuration
(winding/defect/texture), whose energy is accumulated Frank/gradient energy (E=mc^2 = coherent geometric
stress), NOT a band gap. A gapless dispersion does NOT preclude massive coherent textures. Conceded.

CORRECT geometry-first gate: does the geometry support a COHERENT, LOCALIZED, NUMBER-CONSERVED configuration
at z~1100 beyond the baryons? Key equivalence: "redshifts as a^-3" <=> "has a fixed energy scale that does
NOT track the temperature." Thermal (scale ~ T ~ 1/a) -> a^-4 radiation; frozen -> const DE; only a
temperature-INDEPENDENT conserved localized lump (or field with m>>H) -> a^-3 matter. That fixed scale IS the
emergent mass -- so mass is the downstream SIGNATURE of a coherent lump, not an input.

Applying it (geometry-first, no dispersion gap invoked):
  - Bulk structured vacuum: energy set by Lam_cond ~ T ~ 1/a (Phase 0) -> tracks the bath -> a^-4 (thermal) or
    const (frozen CC/DE). Substrate coherence does NOT buy a^-3 (its scale is not T-independent).
  - Localized conserved lumps = HOPFIONS = the baryons (Omega_b~0.05, already counted, already a^-3).
  - Extra non-baryonic coherent textures (director strain/disclinations) need the director ORDERED to exist as
    defects; ordering onsets at T_ord~Lam_cond -> z~few (LATE). At recombination the director is DISORDERED ->
    only incoherent fluctuations (ride the bath -> radiation-like).
  => Still OPTION 3, but the reason is now geometric+sharp: a^-3 requires a T-independent number-conserved
     coherent geometric config; the only ones at recombination are baryonic Hopfions; everything else tracks the
     bath (radiation) or is frozen (DE); non-baryonic textures form too late (z~few).

THE ONE GEOMETRY-FIRST ESCAPE this reframing opens (better question than the mass gate):
  A PRIMORDIAL COHERENT-DEFECT population -- number-conserved localized geometric configs frozen in by an
  EARLIER, higher-T transition (the 2I structure crystallizing), DISTINCT from the late nematic ordering. Such
  lumps would be present AT recombination with fixed rest energy -> genuine a^-3 cold matter, geometry-sourced,
  no separate particle. CONCRETE OPEN QUESTION: does the framework's high-T history contain a symmetry-breaking
  transition leaving a stable, non-baryonic, number-conserved defect/texture network of the right abundance?
  Not currently derived -> today still option 3, but THIS is the direction to probe if reopening.
  (route1prime_mass_gate.py numbers stand as a fluctuation-level check but do NOT settle the texture question.)

## UPDATE 2026-08-30 (5) -- CHECKED the primordial-defect lead against the repo. Reframed to a SPECTRUM question.
Repo grounding (grep + reads):
  - Target space = Faddeev-Niemi director S²: pi_2(S²)=Z (hedgehog monopoles), pi_3(S²)=Z (HOPFIONS = the
    particles). The SO(3)=RP^3 in Paper XIII:1345 is the knot's INTERNAL rotational config space (spin-statistics,
    pi_1=Z_2 -> half-integer spin), NOT a cosmological order-parameter manifold. Paper VII's DM "disclination" is a
    splay TEXTURE (energetic, splay-unstable), forming at the LATE nematic ordering (z~few).
  - grep for relic abundance / stable-neutral dark state / cosmological monopole|wall|defect abundance / freeze-out
    / dark hopfion: NOTHING. Framework has no such object or calculation. Undeveloped, not already solved.
THREE obstacles to the DEFECT-NETWORK version:
  (1) formation epoch: director defects form at the LATE nematic transition (z~few); no high-T Kibble transition
      leaving a network is specified anywhere -> nothing at recombination.
  (2) active-source CMB exclusion: a defect NETWORK as an active/incoherent source gives a broad hump, NOT the
      sharp harmonic acoustic peaks (standard, model-independent for active sources) -> can't be the coherent seed.
  (3) generic Kibble overabundance (~1 defect/correlation volume) = monopole/overclosure problem, not tuned 0.26.
THE REFRAME (correct translation into framework language): "stable, non-baryonic, number-conserved texture" in a
soliton theory = a STABLE, NEUTRAL, MASSIVE, NON-SM HOPFION relic ("dark Hopfion"). This IS standard CDM realised
as the framework's own soliton: mass = geometry (coherent stress), number-conserved, a^-3, clusters. It sidesteps
all three obstacles (passive relic, not active source -> peaks survive; a particle cold from freeze-out, not the
late nematic ordering; abundance = freeze-out number, not Kibble density). Coheres with the whole geometry-first
thread (coherent lumps = mass = matter).
=> The lead survives as a SPECTRUM question, not a defect question:
   (a) does the Q_H tower contain a STABLE, ELECTRICALLY NEUTRAL, NON-BARYONIC, MASSIVE state? -- REPO question,
       Papers XIII/XV/XVII (the particle spectrum). CHEAP next step, answerable from the papers.
   (b) does its cosmological production give Omega~0.26? -- NOT in the repo; a freeze-out/production calc. Do NOT
       fabricate a Kibble-Zurek number: the framework does not yet specify the forming transition to ground one.
STATUS UNCHANGED: still option 3 today. But the open lead is now a well-posed spectrum question (dark-Hopfion
relic), not a vague defect-network hope. Next: read the Q_H spectrum for a stable neutral non-SM state.

## UPDATE 2026-08-30 (6) -- "dark Hopfion" -> "vacuum structure, no feedback" (user). Repo CONFIRMS option 3, structurally.
User correctly rejected "dark Hopfion": a Hopfion is a knot in the MAGNITUDE rho -> it DOES density feedback,
not truly dark. DM in the framework lives in the ORIENTATION n_hat, not rho (Paper VII). Right instinct.
Read Paper VII IC structure (P7:rem:reheating_IC / P7:sec:dm_ic / BAO passage l.1975-2009) -- it ALREADY reasons
exactly this and states the gap:
  - Two cold director-sector components, each disqualified for recombination CDM:
    (i) BARYON-SOURCED strain: adiabatic (inherits P_s via baryons) but TRACKS baryons -> only Omega_b~0.05
        (magnitude fail).
    (ii) SPLAY self-organised texture: independent, but director is ultralight (m_xi~H0) -> Hubble-frozen until
        H~m_xi -> "orders only at z~0, absent at recombination" AND is ISOCURVATURE. Acoustic peaks require
        ADIABATIC cold matter (isocurvature -> out-of-phase peaks, excluded). Two independent disqualifiers.
  - "No density feedback" DOES exist (l.549): high-rho decoupling S_eff=S/(1+beta rho)->0, "preferred direction
    decouples from local physics." But that is SCREENING -- the director effect turns OFF where rho is high ->
    SUPPRESSES clustering. Wrong sign for CDM wells (anti-CDM), not the user's clustering mechanism.
PRECISE MISSING INGREDIENT: a non-baryonic, UNSCREENED, NON-ultralight (m >> H(z_rec)), ADIABATIC orientation/
vacuum configuration -> would cluster as a^-3 at recombination. The framework's orientation sector has ONE mass
scale, m_xi ~ H0, and that ultralightness is WELDED to the DE sector (it is what makes the strain cold and ties
a0 to cH0). The same feature that gives the elegant late-DM + DE unification FORBIDS an early cold component.
Structural, not a loose end. => OPTION 3 CONFIRMED from the paper's own IC structure; it is already the paper's
stated position (P7:sec:dm_ic). The user's "vacuum structure, no feedback" is the correct SHAPE; the framework's
realization of it is ultralight + late + isocurvature by construction, so it cannot be the recombination CDM.
(Supersedes the (5) "dark-Hopfion relic" spectrum lead: the correct home is orientation/vacuum, not a rho-knot
particle; and that home is structurally barred by ultralightness. No spectrum read needed.)

## UPDATE 2026-08-31 (7) -- REFRAME (user, methodological): file the gap as a REGIME to extend, not a particle to find.
User's point: CDM is a THEORY; "cold particle" is its construct, not an observation (50 yr null direct detection).
The OBSERVATIONS are gravitational (rotation curves, cluster dispersions, lensing/Bullet, CMB peaks, LSS, BAO).
"Third peak -> Omega_cdm~0.26" is an interpretation INSIDE GR, not a datum. Don't chase the construct; reproduce
the observations via the framework's OWN mechanism, and where it differs from LambdaCDM that's a TEST, not a failure.
CORRECT (I had been importing the LambdaCDM framing). Honest map:
  - Framework's Frank-elastic mechanism ALREADY reproduces (particle-free): flat rotation curves (M_enc∝r splay),
    RAR / MOND a0 (strain tracks baryons), no BAO shift. Real wins on the GALACTIC data.
  - It is a LOW-acceleration / LOW-density mechanism. At recombination densities it SCREENS to ≈GR (phaseA: 5th
    force ~6% -> ~0% screened). So at z~1100 the framework has NO active distinctive mechanism -> reduces to
    GR+baryons -> underproduces the 3rd peak. Not "missing particle" -- the mechanism is TURNED OFF at that epoch.
  - This is the relativistic-MOND frontier EXACTLY: nails galaxies, struggles at the CMB 3rd peak (recombination
    is not the modified-gravity regime). A shared hard problem of ALL medium/elastic DM theories, not a unique flaw.
    Paper VII itself hedges here ("not fixed here", cutoff-dependent) => status OPEN, not FALSIFIED.
REFRAMED OPEN PROBLEM (replaces "find the recombination CDM"): does the framework's elastic/density-feedback
mechanism, EXTENDED into the high-density recombination regime (where it currently screens off), reproduce the
observed C_l (esp. 3rd peak)? The concrete task = find an UNSCREENED, recombination-ACTIVE piece of the mechanism.
Every prior lead (Bell coherence, transient order, heavy-elastic dressing) was reaching for exactly this; each is
blocked by a specific recorded obstacle (missing scale [Part A]; screening; isocurvature), NOT ruled out.
=> File the gap as "a regime the elastic mechanism does not yet reach," not "a cold particle the framework lacks."
   Option 3 language in the paper (Omega_DM input) stays as the honest placeholder; the RESEARCH target is regime-
   extension + a falsifiable CMB prediction that may DIFFER from LambdaCDM.

## UPDATE 2026-09-01 -- CDM ONTOLOGY (user): CDM = coherent deformation WITHOUT a knot ("wake without the heat").
Reframing that clarifies what CDM IS across the whole sector, and sharpens the recombination question.
FOUR states of the ONE condensate:
  baryons   = coherent KNOTS (Hopfions; localised, EM-coupled).
  DM/CDM    = coherent DEFORMATION (director strain/wake) -- persists INDEPENDENTLY of the knot that sourced it
              (c_dir->0: a deformation relaxes only over cosmological times, survives after the baryon moves on).
              "stress without the heat" = COHERENT geometric strain, cold.
  radiation = INCOHERENT (thermal) condensate excitations (the "heat"), redshift as radiation.
  dark energy = the frozen vacuum (CC).
KEY: CDM is the COHERENT part of the deformation, NOT the incoherent thermal part. This is the mass-vs-temperature
(coherent-vs-incoherent geometric stress) distinction, applied to NAME CDM: cold, coherent, EM-neutral, baryon-
DECOUPLED fabric deformation.
RE-OPENS RECOMBINATION CDM, sharpened: the recombination cold-matter component = a PRIMORDIAL COHERENT DEFORMATION
(a frozen fabric texture left from the early ordering "gluon soup"->structured), persisting (c_dir->0), gravitating,
WITHOUT baryons. Gate 1 earlier found the primordial INCOHERENT (thermal) strain is RADIATION-like -- but that is
the HEAT, which this framing EXCLUDES. The decisive, still-OPEN question: is there a primordial COHERENT (frozen,
ordered) deformation at z~1100, distinct from the thermal condensate, that gravitates as COLD MATTER? = the
strongest form of the primordial-coherent-texture lead. NOT closed by Gate 1 (which addressed the incoherent part).
CLUSTER DEFICIT: the fossil-deformation picture is appealing (accumulated frozen wakes = intracluster CDM) but the
wake calc (each baryon deposits ONE halo, t_form->inf; collective strain SCREENS) gives no accumulation beyond the
RAR value -> clusters stay short. The reframing clarifies WHAT cluster CDM is, not the AMOUNT.
NEXT (recombination): can the early-universe ordering leave a COHERENT frozen strain surviving to recombination as
cold matter? Distinct question from Gate 1's thermal (radiation) result. The primordial-coherent-deformation
derivation is the open lead.
