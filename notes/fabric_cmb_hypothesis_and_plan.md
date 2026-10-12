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

## UPDATE 2026-10-05 -- FOLLOWED the "baryon-seeded topological-amplification" lead (user). CLOSES (turbulence engine is late).
Lead (from the "matter = topologically-LOCKED stress; lock != mass" reframe): a BARYON-SEEDED (->adiabatic)
disclination network, amplified by the SAME cluster-channel multiplicity M_def/M_RAR~(R/xi)^2 (~5x, 28/5), to lift
Omega_b -> Omega_DM AT RECOMBINATION. FOLLOWED; it CLOSES, three ways (first is decisive):
1. THE AMPLIFICATION ENGINE IS TURBULENCE, WHICH IS LATE. The cluster (R/xi)^2 uses xi_net~xi_turb~1/2 R_cl --
   the network is injected at the MERGER-TURBULENCE scale. Turbulence is nonlinear + late (z~few): needs collapsed/
   virializing structure to stir the director and trap disclination lines at small xi. At recombination baryons are
   LINEAR (delta_b~1e-5), no collapse, no turbulence -> NO injection -> the dense network (small xi) simply does not
   form. The amplification cannot be ported to z~1100 because its driver is absent. (Decisive on its own.)
2. If a dense network is postulated anyway: amplification needs structure at xi < baryon-clustering scale = MORE
   structure than baryons carry, but baryons have NO small-scale power at recombination (Silk damping) -> the
   small-scale network CANNOT be baryon-correlated -> ISOCURVATURE -> out-of-phase peaks, excluded.
3. A coarsening (topologically-persistent but evolving) network is an ACTIVE SOURCE -> broad hump, not sharp
   harmonic peaks (standard model-independent defect-CMB result). Baryon-seeding fixes the initial phase, not the
   active evolution.
THE FORK (airtight): adiabatic => baryon-traced => no amplification => capped at Omega_b (~5x short); OR amplified
=> extra small-scale structure => isocurvature/active-source => excluded by peaks. The lead needs BOTH arms; they
are mutually exclusive. The topological LOCK buys "cold + persistent"; the AMPLIFICATION needs frustration/
turbulence, which breaks "adiabatic." Neither supplies the note's named missing ingredient (non-ultralight +
ADIABATIC). => OPTION 3 CONFIRMED from the amplification side. This is the same wall as obstacles (2)-(3) of the
2026-08-30(5) defect-network analysis, now reached via the cluster-channel-at-recombination route. DO NOT re-walk:
porting the cluster (R/xi)^2 amplification to recombination fails because its turbulence engine is a late-time,
nonlinear phenomenon. The lock/ripple/knot intuition is a good REFRAME (lock != mass) but does not cross this gap.

## UPDATE 2026-10-05 (b) -- MODIFIED-GRAVITY route gated (Gates 0+1). CLOSES; option 3 confirmed from the GRAVITY side.
Distinct from the component routes (1,1',2,2' = DM as a component, all closed by the m* gate). This is the
relativistic-MOND frontier: can density-feedback MODIFIED GRAVITY + baryons (no CDM component) reproduce the peaks?
GATE 0 (z_eq without CDM): z_eq(LCDM)=3430 >> z_rec~1090 (matter-dom at recomb, frozen potentials, correct boost);
  z_eq(baryons only)=538 < z_rec => RADIATION-dominated at recombination => potentials DECAY => acoustic driving
  wrong. Confirmed: with only Omega_b, matter-domination is not reached by recombination.
GATE 1 (screening at recomb): condensate rho(z_rec)/rho_today ~ (1+z)^4 = 1.4e12; S_eff ~ 1/(1+beta* ratio) ~
  1.6e-12 => director OFF by ~6e11 AND the chameleon screens the fifth force to GR to ~1 part in 1e12. NO
  modified-gravity rescue: at recombination the framework = GR + baryons + radiation, no CDM.
VERDICT: modified-gravity route CLOSES. Face-value framework predicts a definite, LARGE, LCDM-DISCREPANT CMB
(radiation-dom at recomb: suppressed/shifted 3rd peak, decaying-potential driving, enhanced early-ISW) -- almost
certainly in strong tension with Planck. => OPTION 3 confirmed from the gravity side too; gap robust from BOTH
component and modified-gravity angles. Only remaining non-input avenue = route 2/B (peaks = fabric coherence
transition), a different + speculative mechanism. Script: inline python (Gates 0/1), 2026-10-05.
KEY for the quark<->DM "shared sigma/K" intuition (user): the CDM gap is set by SCREENING+GAPLESS DISPERSION, NOT
by the stiffness MAGNITUDE (Gate 1's 1e-12 is the density ratio (1+z)^4, independent of sigma/K_Faddeev). So the
quark residual's sigma/K CANNOT move the CDM verdict. Also the two phenomena are at OPPOSITE cosmic epochs: quark
deconfinement (QGP) at T>T_c~150 MeV (z>>1e9, long before recomb -- quarks confined by recomb) vs DM activation at
z~few. Same feedback clock, opposite density extremes; a thematic/two-cooling link, NOT a shared recomb-era handle.
sigma/K stays valid only as an INTERNAL-CONSISTENCY cross-check (does knot-scale K match galactic-scale K under
running?) -- worth doing on its own merits, does not touch the CDM gap. [[o3_isospin_residual]], [[reconnection-barrier-redirect]].

## UPDATE 2026-10-05 (c) -- GATE 2: semi-analytic CMB acoustic-peak estimate. no-CDM CMB is GROSSLY EXCLUDED (not merely "tension").
Script: scratchpad/gate2_cmb_peaks.py (tight-coupling photon-baryon + nu + CDM + Phi in conf-Newtonian gauge,
Phi=Psi i.e. neglect nu anisotropic stress). ENGINE VALIDATED on LCDM geometry to 0.1%: r_s(z*)=144.4 Mpc,
D_A=13891 Mpc, theta*=0.01040, ell_A=302 -- all match Planck (theta*~0.01041, ell_A~301.7). So peak POSITIONS
are trustworthy; absolute peak HEIGHTS are NOT (sharp projection + fluid-nu + common damping -> LCDM 3rd peak
comes out too tall), heights usable only QUALITATIVELY.
HEADLINE (robust, pure r_s/D_A integrals): no-CDM (keep omega_b + flat, Lambda->0.95): r_s=195.1 Mpc,
D_A=27728 Mpc, theta*=0.00703, ell_A=447 -- theta* DOWN 48% vs observed. Whole peak series shifts outward;
first peak lands ell~325-340 (after the known ~0.27 ell_A phase-shift that maps the engine's 295 -> true 220)
instead of observed 220. Planck pins ell_A to ~0.03% => conflict at OVERWHELMING significance.
DYNAMICAL confirmation (robust): potential at first-peak scale phi(z*)/phi_prim = 0.67 (LCDM) vs 0.24 (no-CDM)
-> well decayed to 1/4 of primordial by recombination (and still decaying) = the decaying-potential / wrong
acoustic driving Gate 0 predicted, now a number.
HEIGHTS (qualitative): LCDM strong odd/even + tall 3rd peak (H3/H2=1.7, the CDM signature); no-CDM peaks
nearly equal-height, slowly declining (H2/H1=1.09, H3/H1=1.01), washed-out alternation, NO 3rd-peak boost.
NO RESCUE by retuning: tuning H0 to drag theta* back to 0.0104 then breaks z_eq (535 vs 3400), heights, AND
BAO/distance-ladder simultaneously -- the degeneracy CDM breaks; no single-parameter move restores all.
VERDICT: Gates 0/1/2 close the modified-gravity route from THREE independent angles (equality, screening,
peaks). Option 3 robust. The ONLY avenue that could give a VIABLE (because DIFFERENT) prediction is route 2/B
(peaks = fabric coherence transition), not GR-with-missing-matter. [[o3_isospin_residual]].

## UPDATE 2026-10-05 (d) -- CORRECTION to the 2026-08-30 (6) "dark Hopfion" closure. The NEUTRAL-HEAVY-KNOT door is OPEN, not shut.
Trigger: user re-raised "maybe there is a dark Hopfion, one not visible." On close re-examination the (6) closure
("option 3 structurally FORCED") was TOO STRONG: it was forced only for the rho-knot and the smooth splay texture,
NOT for a neutral heavy orientation KNOT. Both pillars of (6) conflate distinct things:
  PILLAR 1 -- "a Hopfion does density feedback -> not truly dark." CONFLATES gravitating with EM-visible. Feedback =
    the knot has a density and GRAVITATES (what we want); "dark" = no EM charge. COUNTEREXAMPLE in the ladder: the
    NEUTRINO (Q_H=1) is EM-neutral, gravitates, IS a knot. neutral+gravitating+knot is already consistent.
  PILLAR 2 -- "orientation sector has ONE mass scale m_xi~H0 -> ultralight -> late." TRUE of the SMOOTH director
    mode (splay texture / DE-welded Goldstone / galactic strain). FALSE for KNOTS: the charged leptons and baryons
    ARE orientation-sector (n_hat) Hopfions with MeV-GeV masses (Faddeev/phi-tower energy, NOT m_xi). So
    "orientation => ultralight" fails for a localized knot.
A NEUTRAL, HEAVY (m>>T_rec=0.26 eV), STABLE orientation Hopfion ("a heavy version of the neutrino's neutrality")
therefore clears BOTH (6) disqualifiers AND the isocurvature one (a freeze-out relic in prior plasma contact is
ADIABATIC, like any thermal relic; also the one-condensate ontology gives adiabaticity structurally). It would be
a^-3 (localized, number-conserved, T-independent coherent stress = mass-is-geometry, claude-hopfion s2), neutral,
cold, early -- EXACTLY the missing ingredient. This REVIVES the (5) spectrum lead that (6) short-circuited.
STILL GENUINELY OPEN (unanswered, not closed):
  (a) EXISTENCE: does the Q_H spectrum PERMIT or FORCE a stable EM-neutral heavy non-baryonic knot? Known ladder has
      ONE neutral knot = the LIGHT neutrino -> such a state is BEYOND the current SM-knot dictionary (allowed by
      ontology, not yet shown forced). Natural candidate: a HEAVY/STERILE partner in the Q_H=1 (neutrino) sector --
      not a new sector. [SPECTRUM READ of Papers XIII/XV/XVII -- in progress this session.]
  (b) ABUNDANCE: freeze-out -> Omega~0.26? Uncalculated; do NOT fabricate (forming transition not specified).
  (c) STABILITY: what conserved charge forbids decay to light-nu + photon?
FRAMING: this IDENTIFIES CDM with a framework soliton (does NOT eliminate it) -> it IS positing a particle, in
tension with the (7) "regime not particle" methodological stance. Both bets live. KEY: "option 3 structurally
forced" downgraded to "option 3 holds for rho-knot + smooth texture; neutral-heavy-knot route OPEN, well-posed,
unexplored." [[o3_isospin_residual]].

## UPDATE 2026-10-05 (e) -- SPECTRUM READ + ontology refinement (user) + the D-TYPE dark-Hopfion lead.
ONTOLOGY REFINEMENT (user): "it's still not really a particle, it's still topology if analogous to a particle."
Correct and it DISSOLVES the (d) tension with note (7): a dark Hopfion is NOT an external point-quantum (WIMP)
added to the ontology -- it is the ONE condensate in a knotted configuration the fabric already holds; mass =
geometric stress as always. "Behaves like a CDM particle" (localized, number-conserved, a^-3, clusters) is the
PHENOMENOLOGY; the ONTOLOGY stays pure geometry. So CDM-as-topology is the SHARPEST form of geometry-first, not a
retreat to a particle. The (7) "regime not particle" stance was against importing an external quantum; internal
topology is compatible.
SPECTRUM READ (Papers XIII/XV/XVII), ESTABLISHED [E]:
  - Matter sectors are McKay/ADE groups; SM fills only A (cyclic SU(2)_1 -> QH=1 nu), E_6 (2T -> QH=3 quark),
    E_8 (2I -> QH=2 lepton). XVII l.411-424 remark: "No D-type group (binary dihedral) appears in the three matter
    sectors; D-type carry their own Z2 outer automorphism, NOT the Z3 colour, and are expected to appear only in
    hypothetical sectors BEYOND THE STANDARD MODEL." => D-series (+ E_7=2O) are RESERVED, unused, colourless slots.
  - Composite tower (XIII) exhausted by SM bosons: Q=0 Higgs, Q=4 gamma/W/Z/gluons, Q=8 graviton -- no dark state.
  - Neutrino neutrality+stability via a Z2 (Cor P17:majorana): j=1/2 primary of SU(2)_1 self-conjugate (j*=k-j=1/2)
    + generates Z2 (1/2x1/2=0) => Majorana. The framework ALREADY uses Z2-self-conjugacy for stability.
D-TYPE DARK-HOPFION LEAD [N/I -- my synthesis, NOT in papers]: a D-type (binary-dihedral) knot hits every CMB
requirement: non-baryonic (no Z3 colour [E]); EM-neutral (plausible: own Z2 decoupled from SM Z6=Z2xZ3 centre --
[O] needs charge projection); STABLE via its OWN conserved Z2 = a DARK PARITY (lightest D-state cannot decay to
A/E sectors -- the textbook DM-stability architecture, R-parity/KK-parity-like; and the framework already uses a
Z2-self-conjugacy for the neutrino, one McKay-branch over -- answers the old open-(c) stability question
structurally); cold/heavy via its OWN knot tower Lcond*phi^{2Q}*e^CS (NOT m_xi~H0 -- it's a knot, not the smooth
mode) [O]; adiabatic+early (freeze-out relic) [generic]; clusters a^-3 (localized coherent stress) [E s2].
FORCED-vs-ALLOWED FORK (decides prediction vs accommodation) -- HYPOTHESIS TO TEST: if matter sectors = finite
SUBGROUPS of the 2I vacuum (|2I|=120=2^3*3*5), then by Lagrange: 2T (order 24 | 120, E_6) available; binary
dihedral 2D_n (order 4n, n|30: 2D_3=12, 2D_5=20 | 120, D-type) AVAILABLE; cyclic (A) available; but 2O (order
48 does NOT divide 120, E_7) UNAVAILABLE. This EXACTLY matches SM usage {A,E_6,E_8} + D reserved + E_7 absent ==>
would FORCE the D-type dark sector as a subgroup the vacuum already contains, and FORCE E_7 absent. Strong
prediction-grade argument IF the papers ground sectors as 2I-subgroups. [VERIFY against papers -- in progress.]
STILL [O]: (a) is the 2I-subgroup grounding actually the framework's principle (forced) or just an ADE labelling
(allowed)? (b) D-type EM charge projection; (c) D-type mass tower -> cold?; (d) abundance Omega~0.26 (freeze-out;
do NOT fabricate). [[o3_isospin_residual]].

## UPDATE 2026-10-05 (f) -- FORCED-vs-ALLOWED resolved + STABILITY (user flagged both). Verdict: PARTIALLY forced, NOT clean.
GROUNDING IS REAL [E]: matter sectors ARE 2I-subgroups. P15 Thm (2T subset 2I, index [2I:2T]=5=pentagon=k+2@k=3;
5 cosets = 5 lepton primaries; the restriction 2I->2T IS the lepton->quark descent). So subgroup arithmetic is
the framework's actual selection mechanism, not an imposition.
FORCED consequences [E/Lagrange]: (i) binary-dihedral (dicyclic) 2D_n EXISTS in 2I (2D_3 o12, 2D_5 o20, Q_8 o8 all
| 120); (ii) E_7/2O EXCLUDED (o48 does NOT divide 120 -> 2O not subset 2I) = a CLEAN standalone prediction: no
binary-octahedral matter sector, ever. Independent of the DM question.
INDEX PROBE (is realization forced, pentagon-style?) -> NEGATIVE. [2I:2D_5]=6, [2I:2D_3]=10. BUT 6,10 are ALREADY
Q_group=2(k+2) of the NEUTRINO (k=1->6, P17:prop:qh1_Qgroup) and LEPTON (k=3->10). So the index COINCIDES with
assigned values -> NO independent coset<->primary mechanism (the pentagon's actual content). s3 trap: number
matching != forcing. => realization ALLOWED, NOT forced.
STABILITY [O] -- CORRECTION to (e): the Z2-parity was OVER-CLAIMED. An outer automorphism is NOT automatically a
conserved dynamical symmetry (contrast the neutrino's genuine FUSION-category Z2, P17:cor:majorana). The real
discriminator is SUBGROUP NESTING:
  - Q_8=2D_2 (o8) IS subset 2T (it is "Set A" of the quark 24-cell, P15:eq:Q8_vertices) -> nested in a realized SM
    sector -> relaxes continuously into the quark, NO barrier -> UNSTABLE. Not a DM candidate.
  - 2D_3 (o12) and 2D_5 (o20) are subset 2I but NOT subset 2T (12,20 do not divide 24; A_4 has no order-6 subgroup
    so neither dicyclic lifts into 2T) -> "SIDEWAYS BRANCHES" -> no symmetric descent into quark/lepton -> decay to
    SM needs RECONNECTION (barriered, Paper XVIII T(2,2)->T(2,3) wall) -> at least METASTABLE; ABSOLUTELY stable
    only if the Z2 is additionally a genuine conserved symmetry (unproven).
  => stability is PLAUSIBLE, UNPROVEN, and ONLY for the sideways dicyclics 2D_3/2D_5 (not a generic D-type property).
NET VERDICT (user's 2-part instinct vindicated): the dark-Hopfion route is PARTIALLY forced (sub-symmetry EXISTS;
E_7 EXCLUDED) but NOT a clean prediction -- realization is UNFORCED (index probe negative) and stability is OPEN
(plausible only for 2D_3/2D_5, metastable-pending-barrier at best). Concrete candidates now NAMED: 2D_3 (o12) and
2D_5 (o20). Next tests if pursued: their KNOT TYPE (distinct from unknot/torus/trefoil?), EM-charge projection,
mass tower -> cold?, and whether the dicyclic Z2 is dynamically conserved. [[o3_isospin_residual]].

## UPDATE 2026-10-05 (g) -- PAPER XVII edit (no-E_7 prediction) + DARKNESS TEST PASSES (sterile).
PAPER EDIT: the no-E_7 forced prediction added to P17:rem:ade_structure (main_paper17.tex). Cites Paper XV
Thm 9.1 (2T subset 2I) -- the "9.1" is INFERRED by section-count, VERIFY against compiled PDF (CLAUDE s11c/d).
DARKNESS TEST -> AFFIRMATIVE/STRONG. Mechanism (P15:thm:charge_generator): electric charge Q=T_3L+T_8L/sqrt3 is a
generator of SU(3)_L (electroweak factor of trinification E_6 ⊃ SU(3)_C x SU(3)_L x SU(3)_R), and charge is
QUANTIZED in units of 1/3 BY THE Z_3 CENTRE of SU(3)_L (all SM charges in (1/3)Z require that Z_3). A D-type group
carries its OWN Z_2, NOT the Z_3 (XVII remark) -> no Z_3 centre -> cannot host SU(3)_L reps -> Q has no SU(3)_L to
act on -> NO electric charge. Same mechanism that makes the neutrino (A-type/SU(2)_1) EM-neutral. GOES FURTHER
than nu: the neutrino is EM-neutral but WEAKLY charged because PAIRED with the charged lepton in the (nu_L,ell_L)
SU(2)_L doublet (XVII l.219); the D-type has NO SM partner (separate McKay branch, own Z_2, orthogonal to the SM
trinification gauge structure) -> no SU(2)_L weak, no colour, no EM -> FULL SM SINGLET = STERILE, couples only via
gravity + its own Z_2 (= the natural dark parity). IDEAL CDM profile.
CAVEATS [O]: (i) the "Z_2 not Z_3 => no SU(3)_L embedding" step rests on the charge-quantization/centre mechanism
(solid) but a from-scratch rep-theory proof that the binary-dihedral WZW admits NO SU(3)_L embedding is NOT done;
(ii) a D-type could carry its OWN hidden-U(1) Z_2-quantized "dark charge" -> hidden-sector/self-interacting DM,
still SM-dark.
CANDIDATE SCORECARD now: non-baryonic [E], DARK/sterile [this test, strong], topologically protected [knot off-
ladder], stability plausible [sideways dicyclic]. REMAINING [O]: stable SOLITON exists in dihedral class?
(construction/relaxation, likely hits the O2 tube-resolution wall), MASS tower -> cold?, ABUNDANCE ~0.26.
[[o3_isospin_residual]], [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (h) -- MASS TOWER for the D-type: NO Q_group (T-matrix fails to close) -> coldness UNDETERMINABLE.
Engine: papers/src_paper20/dtype_cft_tower_engine.py (companion to e6_cft_tower_engine.py; sympy-exact; VALIDATED:
reproduces (E_6)_1->3, SU(2)_1->6, SU(2)_3->10). Q_group = min Q s.t. Q*T_j = odd/4 for every primary (T_j=h_j-c/24,
P15:def:closes_tmatrix). The mass tower m=Lcond*phi^{-2n}*exp(-2pi Q_group T_g) (Paper VI/XVII) NEEDS Q_group.
RESULT: the D-type FAILS T-matrix closure by BOTH natural WZW assignments:
  (a) (D_{n+2})_1 = SO(2r)_1 [quark convention]: SO(10)_1 (2D_3), SO(14)_1 (2D_5) -> Q_group=None. Obstruction =
      the SPINOR primary h=r/8 (r odd forces Q*r ODD, vacuum T=-r/24 forces Q*r EVEN -> irreconcilable).
  (b) SU(2)_k with D-invariant [lepton convention]: 2D_3->SU(2)_6, 2D_5->SU(2)_10 -> Q_group=None. GENERAL pattern
      (verified): ALL EVEN-k SU(2) fail closure; only ODD k close -- and the SM's SU(2)_1, SU(2)_3 are odd. The
      D-type invariants live at even k (6,10).
INTERPRETATION: T-matrix closure IS the charge-quantization condition. Its failure = the D-type carries NO quantized
SM charge = a THIRD independent confirmation of darkness (after Z_2-not-Z_3 and no-SU(3)_L-embedding). CONSEQUENCE:
with no Q_group the mass-tower EXPONENT is undefined; the mass is set only by the undetermined thermodynamic tower
level n -> COLDNESS IS NOT DETERMINABLE from the current framework (weaker even than the neutrino, which at least
has Q_group+T_g; its n_nu is only bounded). ELEGANT: charge-less <=> mass-unpinned -- the dark sector lacks BOTH
quantization conditions together. This does NOT kill the candidate (it did not come out definitively light): the
framework PERMITS a dark D-type of ANY mass; CDM needs it cold (m>>0.26 eV), which is ALLOWED but UNPREDICTED. The
mass is a free parameter at this stage; no established structure pins it (and the CS-spoke route is blocked too --
the D-type knot is off-ladder, its type/existence open). [[dark-hopfion-cdm-lead]], [[o3_isospin_residual]].

## UPDATE 2026-10-05 (i) -- FORMING TRANSITION: QCD-epoch co-formation ("gluon soup") + the MASS reframe (user).
TWO connected user points that reframe the open items:
(1) FORMING TRANSITION = the QCD/confinement epoch ("holding together the gluon soup", T~150 MeV, z~1e12 -- EARLY,
    good for a cold relic). PRECEDENT: baryon formation (two Hopf fibers -> trefoil T(2,3)) topologically FORCES a
    byproduct = the PION (Paper XVIII Thm 5.3). So "assembly sheds a forced topological byproduct" is ESTABLISHED,
    not new. HOOK: Q_8=2D_2 subset 2T (24-cell "Set A") -> dihedral structure is physically present in the soup.
    CAVEAT that sharpens it: the Q_8 INSIDE 2T is the UNSTABLE one (relaxes to quark); the STABLE candidates
    2D_3/2D_5 are NOT subset 2T -- sideways branches of the 2I VACUUM. So the clean picture is "the QCD-epoch
    reorganization of the global 2I vacuum sheds a stable dihedral vacuum-defect" (same epoch, sourced from vacuum
    not baryon). PAYOFF [I, speculative]: if the dark relic CO-FORMS with baryons, its abundance is TIED to the
    baryon abundance -> attacks the Omega_DM/Omega_b~5 "coincidence" (co-genesis at the QCD scale naturally gives
    O(few)). First handle in this thread that could pin the ABUNDANCE rather than leave it free. HOPE not result.
(2) MASS REFRAME: 0.26 eV is the recombination TEMPERATURE (non-rel-at-recomb floor), NOT a mass/average. Cosmology
    pins the ABUNDANCE (Omega_DM h^2~0.12) + "cold", NOT the mass (CDM mass unmeasured: keV->PBH all work). So the
    framework leaving mass free is NOT a defect -- it matches reality. KEY: "cold" differs thermal vs non-thermal.
    THERMAL relic needs m >~ few keV (free-streaming/Lyman-alpha; 0.26 eV thermal = HOT, ruled out). NON-THERMAL /
    TOPOLOGICAL relic (formed cold from a phase transition, axion/defect-like) is COLD AT ANY MASS. So point (1)'s
    QCD-epoch topological formation => non-thermal => COLD BY CONSTRUCTION, regardless of mass => the unpinned mass
    (update (h)) is then CONSISTENT with viable CDM. The two questions resolve each other.
GRADIENT-DESCENT scope decision (user): a relaxation tests SOLITON EXISTENCE/STABILITY (the open structural item),
NOT mass/abundance. HIGH-res -> hits the O2 tube-resolution wall (don't). COARSE existence test IS tractable:
dihedral (2D_3) symmetric IC -> watch Q_H + energy + localization -> qualitative outcome {stable lump / saddle /
unwinds-or-relaxes-to-SM}. Scripts -> NEW DIR src_paper21 (DM/Quark crossover; too big for Paper XX if it pans out).
[[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (j) -- COARSE dihedral relaxation RAN -> hits the resolution wall; INCONCLUSIVE (not a no).
Scripts: src_paper21/dihedral_soliton_relaxation.py (Faddeev E=kappa|grad n|^2+lambda F^2, torch Adam, S^2; Q_H via
the VALIDATED Whitehead integral A_i=Im(Z^dag d_i Z), Q=1 ref ->0.92 at N=48 OK). IC = 3 Q=1 cores on a ring
(genuinely D_3: A3=cos3phi energy-density amplitude = 0.32 nonzero; a naive phase-winding e^{3i phi} IC was
AXIALLY symmetric A3=0 and was discarded).
MAIN RUN (N=48, 3000 steps, 268s): config SHRANK (r_bar 3.2->1.65) then the CHARGE FELL THROUGH THE LATTICE ->
E->0.01 (=vacuum, NOT a finite VK plateau), Q_H->non-integer 0.33, field DELOCALIZED to box edge (r_bar->5.8),
A3->0. = classic lattice-unwinding = the SAME tube-resolution wall that blocked qh3_feedback_saddle (memory).
LAMBDA SCAN (N=48, 300 steps, lam=1/6/25) -- DIAGNOSTIC separating "wall artifact" vs "IC not near a soliton":
larger lambda (bigger soliton) SLOWS charge loss + keeps r_bar larger (lam=1: Q_H->0.38, r->1.65; lam=25:
Q_H->1.59, r~2.2 at step 300). That lambda-dependence = the sub-grid-shrink FINGERPRINT => the unwinding is
(largely) the RESOLUTION WALL, NOT a bad IC. (Mild, do-not-over-read: A3 stays finite while charge is held --
the 3-fold does NOT spontaneously symmetrize to axial before unwinding.)
VERDICT: INCONCLUSIVE on dihedral-soliton EXISTENCE. The coarse free-flow route CANNOT answer it (it hits the wall,
as predicted). It is NOT a "no". Next if pursued: a Q_H-CONSTRAINED relaxation (project/penalise to FIX the charge,
so the question becomes "is there a finite-E localized D_3 config at fixed Q_H" -- the repo gradient_flow_solver
approach), OR accept the wall (as the quark program did) and rest on the ALGEBRAIC results (dark/sterile confirmed
3 ways; no-E_7 forced [in P17]; mass unpinned). Compute discipline: do NOT brute-force higher N -- that IS the wall.
[[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (k) -- KEY (user): the dark soliton SIZE is UNPINNED (like the mass) -> "larger" DISSOLVES the wall.
Faddeev soliton size mu* ~ sqrt(lambda/kappa). For lepton/quark lambda=phi^6 is FIXED by their WZW simple-current
(s_opt, SU(2)_3 Fibonacci). The dark sector has NO closing WZW / NO Q_group -> lambda is NOT forced to phi^6 ->
the dark SIZE is free, SAME under-determination that leaves the dark MASS free (update (h)). So the dark Hopfion
could be LARGER than the (thin) lepton/quark ones -- nothing forbids it. Hints it MIGHT be: off-ladder D-branch
(dihedral) knots are generically fatter than tight T(2,q) tubes; unpinned mass could be LIGHT -> large Compton/
de Broglie (fuzzy-DM-like at the extreme). CRUCIAL: the resolution WALL exists BECAUSE the known solitons are THIN
-- "larger dark soliton" = "larger lambda/kappa" = exactly the regime the wall does NOT bite (the lam-scan already
showed bigger lambda -> charge held far longer). So the wall was an ARTIFACT of assuming lepton-thinness, which the
dark sector is under no obligation to share. => reopens the existence probe as a CONDITIONAL test: "IF the dark
soliton is large (large lambda/kappa), does the dihedral class support a stable soliton?" (sidesteps the wall).
CAVEAT: large lambda is a CHOICE here, consistent with the dark under-determination but NOT derived -- a LIVES
result shows the dihedral class CAN hold a stable soliton, not that the real dark lambda is large. RUN: N=48,
lam=30, 4000 steps -- in progress (bxy3yjj10). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (l) -- FAT-SOLITON run done -> free flow is the WRONG INSTRUMENT (can't even hold the trefoil). Existence NOT numerically accessible.
Fat run (N=48, lam=30, 4000 steps, 352s): vs the thin run, the user's size-caveat WORKED -- config stayed LOCALIZED
(r_bar 1.4-2.3, NEVER ran to box edge) and 3-FOLD PERSISTED (A3~0.37-0.47 throughout, NEVER symmetrised to axial).
BUT the CHARGE still drained (Q_H 2.6->1.0->0->-0.17, E->0.08=vacuum); larger lambda only SLOWED it (zero-cross
~step 2500 vs ~400 thin), never a finite-E plateau.
REAL DIAGNOSIS (clincher): the drain is NOT size/dihedral -- FREE gradient flow cannot CONSERVE Q_H on a lattice
(S^2-projection preserves |n|=1 but NOT the Hopf charge; the repo gradient_flow_solver SAYS this and MONITORS+HALTS
on drift). DECISIVE: that repo solver applied to the KNOWN trefoil ALSO hit the tube-resolution wall (memory, O2).
If the method CAN'T HOLD THE TREFOIL WE KNOW EXISTS, its failure on the dihedral config says NOTHING about dihedral
existence. => the free-flow probe is CATEGORICALLY the wrong instrument at available resolution; soliton existence
is NOT numerically accessible this way. NOT a "no".
POSITIVE MORSEL (real): at large lambda the D_3 structure is ROBUST under flow (stays localized + 3-fold, does not
symmetrise to axial) -- the dihedral config does not WANT to become axial; only the topology-leak (method) ends it.
NET: the gradient-descent route is EXHAUSTED honestly. Existence needs a Q_H-CONSERVING/CONSTRAINED method (real
work; same class that the trefoil itself needs). Recommendation: REST ON THE ALGEBRA (dark/sterile 3 ways; no-E_7
forced [P17]; mass+size unpinned but non-thermal-formation => cold OK). [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (m) -- BIG REFRAME (user): DM = the GRAVITATIONAL RESIDUE of a TRANSIENT dark Hopfion ("wake, knot gone").
User: "if the dark Hopfion has a way to be created it would have existed, and as its charge bled to zero it would
have left gravitational residue, but no longer be with us." This is the framework's OWN deepest DM ontology
(claude-hopfion s4; note 2026-09-01: "CDM = coherent DEFORMATION without a knot, wake without the heat, persisting
via c_dir->0 after the sourcing knot is gone") -- now applied to the DARK sector.
WHY IT'S STRONG:
  * DISSOLVES the stability problem entirely: the dark Hopfion NEED NOT survive -- its TRANSIENCE is the mechanism
    (transient knot -> permanent frozen gravitational residue). All the "is the dihedral soliton stable" numerics
    become MOOT (the charge-drain runs are the FIRST HALF of the mechanism, not a failure).
  * REPAIRS the s4 Omega_b-cap failure: there the strain is BARYON-sourced -> tracks baryons -> only Omega_b. Here
    the sourcing population is DARK HOPFIONS (non-baryonic) -> residue NOT capped at Omega_b -> can reach Omega_DM.
THREE CRUXES [O]:
  (1) SCALING a^-3 vs a^-2 (THE crux): smooth frozen texture = a^-2 (= the route-1' failure). Survives ONLY if the
      residue is LOCALIZED, NUMBER-CONSERVED frozen lumps (each a^-3, = the s2 matter criterion). Localized-lumps
      vs smooth-texture on unwinding is the open question.
  (2) TIMING: unwinding must COMPLETE before recombination. A still-active texture at z~1100 = ACTIVE SOURCE ->
      broad CMB features NOT sharp peaks (standard texture-DM exclusion = note obstacle (2)). QCD-epoch formation
      (z~1e12) gives a huge runway -> plausible, and is WHY the forming epoch matters.
  (3) FREEZING needs c_dir->0, which the gradient-flow SIM LACKS: in the sim the charge-drain sent E->0 (dissipated
      to vacuum, NO residue) because pure gradient flow has no trapping. Physical c_dir->0 would FREEZE the strain
      before full dissipation. So the sim neither supports nor tests the residue -- it is a DIFFERENT mechanism
      (dissipation vs freezing); testing it needs the freezing dynamics, not relaxation.
ANALOG: unwinding cosmic TEXTURES (Turok) leaving gravitational imprints -- real, CMB-constrained as active source,
constraint EVADED by pre-recombination completion. STATUS: most promising dark-sector version yet; NOT a result
(scaling + timing unestablished, abundance uncomputed). This REFRAMES the thread: stop seeking a STABLE dark relic;
the target is a TRANSIENT dark Hopfion + its FROZEN a^-3 residue. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (n) -- EoS of the residue (point 1) + INTERACTIONS decide it (user) -> FIRST dark route NOT structurally barred.
EoS (scaling): residue = frozen director strain. TWO regimes:
  * LOCALIZED number-conserved RIGID lumps (fixed proper size+energy; frozen by c_dir->0, held rigid by the strain's
    genuine STATIC SHEAR RIGIDITY, P7:rem:elastic_solid "static rigidity is large" even as c_dir->0): n ~ a^-3,
    E fixed -> rho ~ a^-3, w=0 => CDM. 
  * SMOOTH space-filling TEXTURE: grad n stretches -> rho ~ a^-2, w=-1/3 => the route-1' FAILURE.
THE DECIDER = INTERACTIONS (user's N-body point): localized-vs-coarsened is set by how the lumps interact. s4:
"as c_dir->0 the tangle freezes but NO true plateau; defect-free ground state wins EVENTUALLY; c_dir sets the
COARSENING RATE t_coarsen ~ 1/c_dir^{1.5-2}". So: a^-3 EARLY (pre-coarsening) -> drifts to a^-2. Is a^-3 AT
RECOMBINATION iff (a) formed as localized lumps EARLY (QCD epoch -> huge runway) AND (b) c_dir below the coldness
threshold (~50 km/s, s4) so t_coarsen >> t_universe. (b) is the SAME coldness the sector already requires -- NO new
ingredient. ESCAPES the s9 a^-2 wall because that wall = the LATE coarsened smooth end-state; the residue is a^-3
in its EARLY localized phase = where recombination sits if formation is early. => FIRST dark-sector route that is
NOT structurally barred. Open piece: localized-lumps-vs-texture on unwinding (= the freezing/coarsening dynamics).
INTERACTIONS (user, physics): one SHARED n_hat (s7 monism) => inherently MANY-BODY, not independent relics.
  * dark-dark: overlapping strain dressings -> attract/repel/capture ("whirlpools, stick or fly off"); bound dark
    composites (Paper XIII Q=4 mechanism, dark sector); => SELF-INTERACTING DM (SIDM) = falsifiable (core-cusp,
    cluster offsets), a FEATURE.
  * dark-baryon: via shared n_hat -> enables QCD co-formation (abundance hook) BUT risks partial baryon-tracking
    (toward Omega_b cap); not fully capped (dark Hopfions are their own source). Real tension, flag it.
  * N-BODY = the s4 KIBBLE-ZUREK frozen disclination NETWORK (the cluster-mass channel), now dark-Hopfion-sourced
    instead of merger-turbulence-sourced. SAME object -> the dark CDM residue and the cluster-mass network may be
    ONE frozen-defect-network physics. Full sim = N interacting field-solitons each at the resolution wall ->
    intractable; but the NETWORK physics is already in the framework.
STATUS: strongest result in the thread -- clears the a^-3 hurdle that killed every earlier director route, on the
SAME coldness the sector already has. NOT a derivation (scaling-level; a^-3 is long-lived METASTABLE "no true
plateau"; localized-lump formation + abundance unestablished). [[dark-hopfion-cdm-lead]], [[o3_isospin_residual]].

## UPDATE 2026-10-05 (o) -- POINT 2 (freezing dynamics) via the CORRECT INERTIAL tool -> CONFIRMS the a^-3 EoS.
KEY method point: my src_paper21 relaxation was OVERDAMPED (Adam=steepest descent) -> dissipated E->0 (no freeze).
WRONG dynamics. The framework is INERTIAL (chi d_t^2 n = K lap n); the right tool ALREADY EXISTS:
papers/src_paper7/dm_polarization/defect_network_sim_3d.py (inertial velocity-Verlet, c_dir=sqrt(K/chi); measures
c=xi_net/xi_net0: c~1 freezes localized=a^-3, c>>1 coarsens=a^-2; --requench = merger re-seeding).
RAN (smoke N=48, single-quench freeze test), c_dir scan: c_dir=0.08 -> c=1.06 (STAYS LOCALIZED = a^-3); c_dir=0.2
-> c=1.16; c_dir=0.5 -> c=21, P->0 (coarsens away). Sharp threshold 0.2->0.5. => CONFIRMS point-1 EoS from the
correct dynamics: residue is a^-3 (localized frozen) IFF c_dir below the coldness threshold (s4 ~50 km/s, the SAME
coldness the sector already needs). My earlier "dissipated to vacuum" was the overdamped artifact; inertial FREEZES.
REQUENCH (user's merger/whirlpool interactions): 3 mergers @ c_dir=0.2 -> P re-spikes at each merger, ends 378 vs
330 single-quench, c=1.08 vs 1.16 => INTERACTIONS ARE STABILIZING: ongoing dark-Hopfion activity re-seeds the
network, fights coarsening, holds it localized (a^-3). A feature, not just a complication.
UNIFICATION: EARLY passive frozen residue (recombination CDM, single-quench-then-freeze, must be PASSIVE by z~1100
or it's an active source = wrong peaks) and LATE merger-sustained network (= the s4 cluster-mass channel, ongoing
mergers help) are the SAME frozen-defect-network physics at two epochs.
CAVEATS: generic KZ-quench PROXY (not dark-Hopfion-specific seeding -- that needs the dark Hopfion = resolution
wall); smoke/coarse (full s4 runs N=96/8000 gave the "no true plateau, t_coarsen~1/c_dir^{1.5-2}" result); c_dir in
code units (physical threshold mapping from s4). NET: points 1+2 agree -- dark-Hopfion-wake CDM is a^-3-VIABLE
(conditionally: cold enough c_dir + passive-by-recombination + localized seeding). Strongest dark route; still not
a derivation of Omega_DM. [[dark-hopfion-cdm-lead]], [[o3_isospin_residual]].

## UPDATE 2026-10-05 (p) -- ABUNDANCE via QCD co-formation (ADM) -> PREDICTS m_dark ~ 5 GeV, EXPLAINS Omega_DM/Omega_b.
Engine: src_paper21/adm_abundance.py. Residue gravitates with ~ dark-Hopfion rest energy (c_dir->0 traps the
unwound energy, little radiated -> E_res ~ f*m_dark, f~1). If the dark Hopfion is a FORCED topological byproduct of
baryon formation (like the PION one-per-baryon, P18 Thm 5.3), its number tracks the baryon ASYMMETRY:
n_dark/n_b = r (topological, not a free knob). Then Omega_DM/Omega_b = r*f*(m_dark/m_p) -- pure ADM arithmetic.
OBS Omega_DM/Omega_b = 5.36 -> require r*f*m_dark = 5.03 GeV. FRAMEWORK-NATIVE (r=1 one-per-baryon, f~1 full trap):
=> m_dark ~ 5.0 GeV = the CANONICAL ADM mass (~5 m_p). This simultaneously (1) EXPLAINS the Omega_DM/Omega_b~5
coincidence with NO fine-tuning (an O(few)*QCD-scale object, which QCD-epoch co-formation gives for free), (2) PINS
the previously-unpinned dark mass (update (h)), (3) ~5 GeV + STERILE -> evades 50yr direct-detection nulls (feature).
SELF-CONSISTENT timing: 5 GeV can't be THERMALLY made at T_QCD~150 MeV (Boltzmann-killed) but it's NON-THERMAL /
TOPOLOGICAL (like the proton 0.94 GeV forming at T~150 MeV via confinement, not thermal) -- the non-thermal
co-formation the whole picture needed.
CAVEATS: (a) r = topological byproduct ratio ASSUMED O(1), NOT derived -- now THE key open piece (is the dark
Hopfion a FORCED byproduct of baryon formation, how many per baryon? -- a topological-bookkeeping question the
framework could answer but hasn't); (b) f~1 assumes near-total trapping; (c) eta (baryon asymmetry) NOT explained
(baryogenesis, SM-open too -- ADM TIES dark to eta, doesn't derive eta); (d) must be ASYMMETRIC (symmetric
annihilating dark Hopfions -> thermal relic DECOUPLED from Omega_b -> coincidence NOT explained).
NET: the dark-Hopfion-wake CDM now has a coherent END-TO-END (conditional) story: dark/sterile (3 ways) + a^-3
viable (pts 1+2) + abundance via ADM pins m_dark~5 GeV explaining the coincidence. FALSIFIABLE: ~5 GeV sterile ADM;
SIDM (dark-dark); recently-merged clusters denser. The remaining UNCOMPUTED pieces are r (the forced-byproduct
ratio) and eta. [[dark-hopfion-cdm-lead]], [[o3_isospin_residual]].

## UPDATE 2026-10-05 (q) -- BOOKKEEPING CLOSES THE LYNCHPIN: dihedral byproduct NOT forced -> ADM coincidence FAILS.
Template (P18 thm:pion_forced): baryon formation 2+2 -> 3+1 forces the PION via TWO saturated conservation laws:
(i) Hopf charge/linking 4=3+1 (minimal decomp with a confined Q_H=3); (ii) component count (4 circles, 3 fuse to
trefoil, 4th is the forced spectator). The pion SATURATES both books -- nothing left over.
APPLY to a dihedral byproduct: for a D-type object to be FORCED, the reaction must carry a conserved quantity only
a dihedral object can absorb. Enumerate: Hopf charge (balanced=pion), components (balanced=spectator), EM/colour/
E/p (balanced by quark+pion), and the D-type's distinguishing charge = its OWN Z_2 (dihedral outer automorphism)
is carried by NONE of the ladder objects (2I lepton, 2T baryon, SU(2)_1 pion) -> conserved TRIVIALLY 0=0 -> forces
NOTHING. Also un-assemblable: a dihedral knot is a different SYMMETRY CLASS, not a different Q_H -> cannot be built
by charge-conserving fusion of torus knots/unknots (you get torus knots + unknot spectators, never a dihedral
linkage; many-body soup reactions don't change this). => r_forced = 0. The dark Hopfion is NOT a forced byproduct
of baryon formation. (Subtler 2D_n-not-subset-2T "expulsion" idea: not a clean conservation law, needs the vacuum
locally dihedral, not grounded -- not leaned on.)
CONSEQUENCE -- ABUNDANCE CLOSES IN THE NEGATIVE: the ADM explanation of Omega_DM/Omega_b~5 REQUIRED r~O(1) tied to
baryon number (forced co-formation). The books say that tie does NOT exist. => the ~5 GeV prediction (update (p))
EVAPORATES (it rested entirely on forced co-formation). If the dark Hopfion forms, it forms as a KIBBLE-ZUREK
vacuum defect (density set by correlation length, NOT tied to Omega_b) = the framework's OWN obstacle (3): generic
KZ -> OVERCLOSURE/UNTUNED, not 0.26. So the Omega_DM/Omega_b coincidence is NOT explained; abundance NOT derivable.
NET (the "close it out one way or another" the user asked): the co-formation/ADM coincidence story FAILS the
framework's own bookkeeping. SURVIVING: dark/sterile (3 ways), a^-3-viable IF formed+cold (pts 1+2), no-E_7 forced
[P17]. The CDM ABUNDANCE remains UNEXPLAINED -- back to option 3 on the recombination budget, now with a fully
characterised dark sector and one firm prediction. Abundance tag: UNSUPPORTED-HOPE (was the exciting ADM lead).
[[dark-hopfion-cdm-lead]], [[o3_isospin_residual]].

## UPDATE 2026-10-05 (r) -- CORRECTION (user): (q) used the DILUTE regime; check the SATURATED SOUP too. Both fail, differently.
Two user points: (1) "partially forced" (e/f) = the dihedral SECTOR is AVAILABLE (2D_n subset 2I forced; E_7 absent)
= it CAN form; "r_forced=0" (q) = not forced as a per-baryon BYPRODUCT. AVAILABLE != FORCED-TO-BE-PRODUCED -> defects
CAN exist (abundance != 0), just not shackled 1:1 to baryons. (q) was sloppy to read as a flat death.
(2) ENVIRONMENT: (q)'s reaction 2+2->3+1 is the DILUTE/LATE isolated-object regime -- NOT the saturated gluon soup.
In the dense plasma the bookkeeping differs: high feedback, overlapping preimages, COLLECTIVE Kibble-Zurek freeze-out
(not clean 2-body); a dense tangle HAS all symmetry classes locally and DOES shed dihedral defects -> the dilute
"nothing left over" does NOT apply. Real gap; (q) only covered the dilute channel.
BUT evaluating the SOUP channel (honest): KZ defect density ~ 1/corr-vol ~ (few T_c)^3 ~ O(n_gamma); baryon density
= the ASYMMETRY n_b = eta n_gamma, eta~6e-10. => n_dark/n_b ~ 1e7-1e10 -> Omega_dark/Omega_b ~ 1e7-1e10 at ~GeV =
CATASTROPHIC OVERCLOSURE (not ~5). To avoid: defects must annihilate/coarsen -> frozen (still too much) OR scaling
(ACTIVE CMB source -> wrong peaks) = the framework's own obstacle (3), now via the dense-regime physics.
NET: the ADM coincidence needed n_dark ~ the baryon ASYMMETRY. DILUTE channel: no forced tie to n_b (bookkeeping).
DENSE channel: sets n_dark by KZ DENSITY ~n_gamma = ~1e10x the asymmetry -> overproduce. NEITHER gives Omega_DM/
Omega_b~5. Conclusion of (q) HOLDS but now BOTH regimes checked (the user's correction). The dihedral SECTOR is
available (can form dark matter) but its abundance is NOT naturally 0.26: too-untied (dilute) or too-much (dense).
Still OPTION 3 on abundance. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (s) -- DARK BARYOGENESIS route: the D-type is a STERILE but DIRAC sector (complex spinor). Open direction + SO(10) hook.
Framing fix: a Z_2 = STABILITY (dark parity), NOT an asymmetry. An asymmetry needs a conserved U(1) (dark number).
Q: does the D-type have one? YES. Candidates 2D_3<->D_5<->SO(10) (r=5), 2D_5<->D_7<->SO(14) (r=7). Reality of
SO(2r) spinors (textbook): r ODD -> spinor COMPLEX (s*=c distinct) -> Dirac-like -> conserved dark-number U(1) ->
ASYMMETRY POSSIBLE. Both candidates r odd -> complex (SAME fact that makes SO(10)'s 16 complex/chiral). Contrast
the neutrino (SU(2)_1 Majorana, self-conjugate -> NO asymmetry). ELEGANT COHERENCE: the odd-r spinor (h=r/8) that
FAILED T-matrix closure (=> sterile, update (h)) is the SAME one that is COMPLEX (=> asymmetry-capable). Sterility
and Dirac-ness share one origin. FIRM (rep theory), new structural result.
SO(10) HOOK: 2D_3's algebra = D_5 = SO(10) = THE GUT group, home of B-violation + the chiral 16. Natural site for
CO-GENESIS (correlated baryon+dark asymmetries). LEAP (flag): framework uses SO(10) as McKay-correspondent ALGEBRA,
NOT a gauged GUT with dynamics -- treating it as a real B-violating GUT is unbuilt.
OBSTACLES (why direction not result): (1) NO baryogenesis foundation -- framework doesn't explain eta; can't
correlate two asymmetries when one isn't derived. (2) STERILITY vs SHARING (standard ADM tension): what makes it
dark makes tying its asymmetry to the baryon one hard; need an early-active/late-decoupled CONNECTOR -- the shared
2I VACUUM (both 2T,2D_3 subset 2I; a vacuum-level CP-violating reorg could seed both) is the natural framework
candidate, but model-building. (3) MAGNITUDE n_dark~n_b model-dependent, not derived.
VERDICT: OPEN direction with a REAL structural opening (sterile-yet-Dirac sector on the SO(10) algebra, 2I vacuum
as connector) -- NOT a result; pins abundance only CONDITIONAL on physics not yet in the framework (dark co-genesis
+ a baryogenesis mechanism). FOR PAPER 21: confirms the call -- firm paper-ready content is the STRUCTURE, now
PLUS this (D-type = sterile DIRAC sector; sterility+Dirac share one origin; SO(10) is its algebra). DM abundance
stays an OPEN APPLICATION with a named next step (dark co-genesis on 2I), not a dead end. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (t) -- BARYOGENESIS mechanism (user, foundation step): chiral reconnection, all 3 Sakharov conds framework-native.
Repo check: NO baryogenesis mechanism in the papers (only derived CKM phase delta_CP, Paper 5 thm:delta_cp, two
competing WZW; neutrino CP Paper 10) -> eta genuinely unexplained. INGREDIENTS all present though.
MATTER/ANTIMATTER = trefoil CHIRALITY: T(2,3) is a chiral knot (writhe +3 != mirror -3); baryon = trefoil wound
one way, antibaryon = mirror winding. P18:rem:wound_state: trefoil is a WOUND STATE "encoding the history of its
formation" (assembled by chiral threading, not a minimum). Reactions REVERSIBLE (P18 inverse 3+1->2+2) => in
equilibrium balance, no net B => asymmetry needs OUT-OF-EQUILIBRIUM.
SAKHAROV MAP (all framework-native, assembled not postulated):
  (1) B-VIOLATION: reconnection 2+2->3+1 (& reverse) -- baryons FORM FROM LEPTONS by threading; B = knot-sector
      count changed by reconnection, not fundamental. [P18]
  (2) C&CP: derived WZW phase delta_CP [Paper 5] + chiral writhe/Chern-Simons -> biases trefoil-vs-antitrefoil
      threading direction.
  (3) OUT-OF-EQ: confinement transition FREEZES the reversible reconnections ("confinement as threading
      irreversibility", P18) -> CP-biased ratio frozen in.
NET: chiral CP-violating reconnection preferentially winds one trefoil handedness; confinement freezes the excess
before washout -> that excess = eta.
FAVORABLE: EVADES the "SM CKM CP ~10 orders too small for baryogenesis" objection -- that is specific to EW
SPHALERONS; here B-violation is TOPOLOGICAL RECONNECTION at the confinement scale (different, plausibly far more
efficient channel, where the CP phase is O(1)).
HONEST LIMIT: the STRUCTURE (CP-biased chiral freeze-out) maps; the MAGNITUDE eta~6e-10 is NOT derived (needs
reconnection rates x CP interference x freeze-out dilution, unbuilt -- the hard step in every baryogenesis model).
Framework handle to try: the writhe-3 framing phase theta_{3/2}^3 = i (P15, a genuine chiral phase in the
reconnection amplitude). NOT faked.
CONNECTOR (next step, previewed): ONE mechanism, not two -- the chiral CP reconnection happens IN THE 2I VACUUM
(contains BOTH 2T baryon and 2D_3 dark); the SAME biased freeze-out biases the dark handedness too, and the
D-type's COMPLEX/Dirac spinor (update (s)) gives it the dark-number U(1) to carry the bias -> CO-GENESIS ->
n_dark ~ n_b naturally. BUILD NEXT. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (u) -- DARK-BARYON CONNECTOR: ADM co-genesis resolves the (r) overproduction. m_dark~5 GeV re-emerges (conditional).
CONNECTOR = the SHARED WZW phase: delta_CP is a property of the 2I WZW T-matrix (Paper 5), common to the whole
vacuum -> the SAME delta_CP biasing the chiral trefoil(baryon) freeze-out biases any chiral 2D_3(dark) freeze-out
in the same vacuum -> correlated asymmetries by construction (CO-GENESIS).
RESOLVES the (r) OVERPRODUCTION via standard ADM, all pieces framework-native:
  1. KZ abundant SYMMETRIC formation (n_dark ~ n_gamma).
  2. SELF-ANNIHILATION removes the symmetric part -- dark Hopfions are SIDM (dark-dark via shared n_hat); dark+
     anti-dark annihilate.
  3. CP-biased ASYMMETRIC residue survives = the DM; size ~ delta_CP x (branching into dark sector).
  4. if dark asym ~ baryon asym (same delta_CP, comparable subgroup branching 2D_3:2T = order 12:24) ->
     n_dark^res ~ n_b -> ADM -> m_dark ~ 5 GeV RE-EMERGES.
=> the (r) overproduction is what ADM WANTS (abundant symmetric production, burned to the asymmetric residue by
annihilation). The dark SIDM (flagged as falsifiable feature) IS the required annihilation channel -- pieces
interlock.
HONEST CATCH: this is the STABLE-RELIC version -> RE-RAISES the existence/stability question (ADM needs the
asymmetric residue to SURVIVE as stable dark Hopfions, Z_2-protected lightest). DIFFERENT BET than the transient-
WAKE reframe (m): wake AVOIDED needing stability; ADM REQUIRES it. Can't have both; ADM picks stability (Z_2-
plausible, numerically wall-blocked). Compounding unresolved: branching ~O(1) ASSUMED (comparable subgroups); SIDM
annihilation efficiency ASSUMED; rests on the UNBUILT eta magnitude (update t).
NET: abundance now has a COHERENT framework-native MECHANISM (ADM co-genesis: one chiral CP event in the shared 2I
vacuum -> correlated asymmetries; SIDM burns the symmetric overproduction; asymmetric stable residue = DM at ~5
GeV). UPGRADE from (r) "both regimes fail". But a CHAIN of plausible O(1) assumptions on an unbuilt eta + unproven
stability = a DEVELOPED CONDITIONAL SCENARIO, not a derivation. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (v) -- PAPER 21 DRAFT + eta ATTEMPT (lynchpin stays open, but the obstacle is now precisely located).
PAPER 21 drafted: papers/main_paper21.tex (checkpoint). Thesis = STRUCTURE (firm): ADE completion + no-E_7 theorem;
D-type sterile (T-matrix closure fails both WZW ways, engine-validated) + Dirac (complex spinor r odd); sterility
+Dirac share the odd-r-spinor origin; SO(10) algebra. DM = CONDITIONAL application (a^-3 EoS; ADM co-genesis ->
m_dark~5 GeV) with 4 open conditions (sec 8). Epistemic scope stated up front. Cites I,V,VII,XV-XX,CIZ; convention
fixes done (no backslashed underscores; \path for filenames).
eta ATTEMPT (the lynchpin): handle = writhe-3 framing phase theta_{3/2}^3 = i (P15). KEY SUBTLETY: a SINGLE CP
phase gives NO asymmetry (trefoil writhe +3 -> phase i, mirror -3 -> i*=-i; |i|=|-i| -> equal rates at tree level).
Asymmetry needs TWO interfering channels with different weak AND strong phases (eps ~ sin d_weak sin d_strong).
FRAMEWORK HAS THEM: Paper 5 derives delta_CP from TWO COMPETING WZW MECHANISMS -> the interference exists, biased
by the framing phase i -> eps (CP asym per event) plausibly O(delta_CP) ~ O(0.1), NOT tiny/forbidden. BUT eta !=
eps: eta ~ eps x (carrier/entropy) x (washout efficiency). Landing eta~6e-10 from eps~0.1 needs a ~1e-9 dilution/
efficiency factor = a full non-eq reconnection-rate BOLTZMANN calc the framework has NOT built.
VERDICT: eta NOT derivable now. BUT the obstacle is precisely located: eps (the CP asymmetry) is NOT the problem
(framework handle exists, O(0.1)); the missing piece is the EFFICIENCY/dilution (~1e-9) = reconnection-rate
Boltzmann. The lynchpin stays open, now with the exact missing calculation named. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (w) -- CORRECTION (user caught a conflation): STABILITY is topological (solid), only EXISTENCE is wall-gated. Charge does NOT physically decay.
User: "I thought the large dark Hopfion was stable but the charge decayed -- is it full conservation we want, or
stability that decays due to charge decay?" This exposed a conflation of TWO pictures I slid between:
  A = STABLE KNOT (ADM relic): Q_H conserved, knot survives, net asymmetry = DM. (what update (u)/5 GeV needs)
  B = TRANSIENT WAKE (update m): knot unwinds (charge decays), frozen wake is the DM.
RESOLUTION: the charge is the HOPF charge Q_H (topological). Q_H conservation = stability against UNWINDING
(SAME protection as the trefoil); Z_2 dark parity = stability against decay to SM sectors; Dirac/complex spinor
(r odd, update s) => dark Hopfion != conjugate => NET Q_H can be nonzero = the conserved dark-number the ADM
asymmetry lives in. So stability + asymmetry BOTH from Q_H-conservation + Dirac, and the charge does NOT physically
decay. The two pictures COLLAPSE to A; Q_H conservation FORBIDS the B unwinding. (B is also ADM-incompatible: if
the knot unwinds, its net dark-number is lost -> 5 GeV evaporates.)
THE SIM'S CHARGE DECAY WAS AN ARTIFACT: free gradient flow on a lattice doesn't respect Q_H conservation (leaks it
through the grid) -- KNOWN method failure because the SAME method couldn't hold the trefoil either (which exists &
is stable). NOT a real unwinding.
RECORD CORRECTION: earlier updates listed "stability" as unproven/wall-blocked -- that CONFLATED stability with
EXISTENCE. STABILITY (if a dihedral Q_H soliton exists, does it last?) = topological + Z_2 = SOLID, as firm as the
trefoil's, NOT open. EXISTENCE (is there a finite-E minimizer in the dihedral Q_H sector?) = the genuinely open,
wall-gated question. OPEN LEDGER SHRINKS to: (1) existence (charge-conserving relaxation, wall) + (2) eta
efficiency (reconnection-rate Boltzmann). Stability is off the open list. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (x) -- PAST THE WALL: literature says it's mostly a WRONG-TOOL gap + one real physics question. Methods identified.
Web lit review (Faddeev-Niemi/Skyrme Hopfions -- a MATURE numerical field):
STANDARD MINIMIZATION = rational-map ansatz IC + BFGS/quasi-Newton on 180^3-240^3 lattices (spacing ~0.08)
[Battye-Sutcliffe; Sutcliffe]. The repo used ADAM/OVERDAMPED on ~48^3 -- the one method guaranteed to leak charge
+ dissipate. So our "wall" is largely WRONG TOOL (steepest descent, coarse) vs the field's good-ansatz+quasi-Newton
+fine-grid. Existence reframed: blocked -> an established computation not yet run with the right tooling.
EXISTENCE is partly a THEOREM: Lin-Yang (CMP 2004) -- energy minimizers exist for infinitely many Hopf charges in
the Faddeev model; E <= C|Q|^{3/4} (Vakulenko-Kapitanskii). Firm analytic footing for the base model.
CONSTRUCTION MATCHES FRAMEWORK: Sutcliffe "Torus knots as Hopfions" (2013, arXiv:1304.6021) builds (P,Q)-torus-knot
Hopfions as |Q| SINE-GORDON kink strings twisted P/Q -- exactly the framework's trefoil-as-torus-knot (+ resonates
w/ the sine-Gordon tube-mode notes). Clean symmetric ICs.
eta / RECONNECTION RATES: COLLECTIVE-COORDINATE / MODULI-SPACE methods + product ansatz for interactions +
instanton-holonomy moduli space [isospinning hopfions 1301.6803; scattering IOP 2012; instanton-holonomy
hep-th/0108082]. Semi-analytic reconnection/scattering rates WITHOUT full-grid per-knot relaxation = the route to
eta's efficiency factor, and addresses the user's N-body/whirlpool dynamics.
THE REAL HARD PART (honest, not just tooling): standard Faddeev minimizers have AT MOST discrete symmetry, in
PRACTICE only C_2; first knot (trefoil) at Q=7 NOT Q=3. => (a) a DIHEDRAL (D_3) minimizer is NOT guaranteed (base
model favors LOW symmetry); (b) the density-feedback functional DIFFERS from standard Faddeev (assigns trefoil to
Q_H=3) so the standard catalog does NOT transfer -- framework needs its OWN minimization w/ its OWN functional.
BUT: the framework ALREADY relies on density-feedback to get the trefoil's 3-fold at low Q -> the dihedral question
is the SAME symmetry-imposition mechanism, not a new leap.
PATH FORWARD: (1) replace solver: rational-map IC + BFGS/arrested-Newton (charge-conserving) + finer grid + the
density-feedback functional -> direct route to existence (240^3 big but routine, NOT a fundamental wall). (2) impose
D_3 in the ansatz (symmetric minimization, torus-knot/sine-Gordon). (3) eta rate via collective coordinates/moduli
space (not full-grid). (4) analytic existence: adapt Lin-Yang to density-feedback + D_3-restricted. NET: wall =
tooling gap + one real physics Q (does density-feedback admit a D_3 minimizer?), both DOABLE research programs.
[[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (y) -- SOLVER BUILD STARTED (shared quark+dark): E_fb functional DONE; IC needs rational-map (tube-lift fails). + saddle correction.
Key design CORRECTION (from reading qh3_feedback_saddle.py): E_fb = K_fb*J4 is a PRODUCT -> NO Derrick minimum ->
the density-feedback Hopfion is a SADDLE (Paper I; true for Q_H=2 & 3). So BFGS/Adam MINIMIZATION COLLAPSES it; the
right tool is a STATIONARY-POINT / self-consistency SADDLE search (score = |V-phi|+|s_opt-1|, snapshot best), NOT
minimization. My earlier "rational-map + BFGS" was right on ansatz+resolution, WRONG on optimizer (BFGS is for the
bare SUM-form Faddeev w/ a minimum; the framework's PRODUCT form has none).
RESOLUTION WALL QUANTIFIED: old qh3 run N=64 h=0.175 = ~2 pts/tube-radius; lit standard N~200-240 h~0.08 = ~5
pts/tube. Wall = under-resolution (~1/4 of needed), fixable at N~200 but HPC/GPU, NOT fundamental.
BUILD (src_paper21/feedback_saddle_core.py, reusable quark+dark): E_fb functional + diagnostics (J2a,J2iso_fb,J4,
K_fb,V,s_opt,score) + Whitehead Q_H = BUILT & CORRECT (Q_H routine validated -> 0.92 on the Q=1 inverse-stereo
reference at N=48). IC (knotted Hopfion) = NOT YET: the ported nearest-strand tube-lift gives Q_H=0 -- only 0.5% of
grid points find the 2nd strand, |z2|~0.03 => effectively SINGLE-strand => Q_H=0 identically. The knot's charge
lives in the THIN 2-strand crossing regions, which a tube-lift under-samples to nothing. THIS IS WHY THE FIELD USES
THE RATIONAL-MAP ANSATZ (smooth, charge properly distributed) -- NOT a tube-lift. So the IC is the real next build
step: implement Sutcliffe "torus knots as Hopfions" rational map (trefoil = T(2,3)), validate Q_H at N~150-200,
THEN the high-res saddle search (HPC). Once validated on the trefoil it RE-OPENS the O2 quark saddle (pins E_6 tower
Dn_q) AND sets up the dark (dihedral-IC still TBD -- which knot, off the torus ladder). [[o2-feedback-saddle-solver]], [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (z) -- COURSE-CORRECTION: prior O2 record shows the resolution wall is PHYSICAL; numerical route settled.
Re-read [[o2-feedback-saddle-solver]] (detailed v6-v11 history). It TEMPERS this session's "rational-map+high-res
re-opens O2" enthusiasm. ESTABLISHED there (not tentative): (1) the jammed density-feedback tube is THIN (radius
~0.36-0.40 at EVERY beta; the "does feedback fatten it" check = NO) -> no IC changes the relaxed object's thinness.
(2) BOTH FD and Berg-Lüscher solid-angle J4 estimators FAIL to converge; solid-angle WORSE; need ~10-15 pts/tube
-> h<=0.03 -> N>=300 (27M+ pts, hours/1000 steps) IMPRACTICAL. (3) tube-adapted curvilinear solver ALREADY built
(tube_{curvilinear,trefoil,relax}.py); |grad E|^2 saddle search attempted, ill-conditioned. (4) CONCLUSION: quark
O2 pinned ALGEBRAICALLY (e6_cft_tower_engine, validated), NOT numerically. So "rational-map+BFGS re-opens O2" =
re-litigating settled ground; my fresh lit read (rational-map+BFGS routine) is for the BARE SUM Faddeev (has a
minimum), NOT the framework's PRODUCT E_fb (thin jammed saddle, physical wall).
SALVAGE (honest): feedback_saddle_core.py functional+diagnostics+Q_H are correct & reusable. The smoke Q_H=0 was a
real IC bug (tube-lift -> single-strand); the RATIONAL-MAP ansatz WOULD fix the IC charge -- but as an IC for a
future rented-HW run, NOT a wall-breaker (relaxed saddle still thin). DARK sector is WORSE: the algebraic escape
that saved the quark (CFT tower) is BLOCKED there (D-type has no Q_group, T-matrix doesn't close, update h). So
dark EXISTENCE is stuck BOTH ways: numerical (thin-tube wall) + algebraic (no tower).
NET: I was over-optimistic. Numerical field-relaxation route is wall-blocked by the PHYSICAL thin tube (my own
record). Options: (a) finish rational-map IC as a clean reusable IC for future HPC only (modest); (b) tube-adapted
curvilinear route (the record's actual requirement, real HPC); (c) accept existence wall-gated both sectors, rest
on algebraic/structural results. Defer to prior direct work over fresh inference (CLAUDE s4,s5,s9).
[[o2-feedback-saddle-solver]], [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (aa) -- COLLECTIVE-COORDINATE DYNAMICS is DUAL-PURPOSE (quark + dark) + fills a stated Paper XVIII gap.
Checked reuse for the quark sector (user). FINDINGS: (1) NEW direction -- no collective-coordinate/moduli-space/
product-ansatz work in the repo. (2) Paper XVIII EXPLICITLY flags it as missing: develops the reaction TOPOLOGY
(2+2->3+1 bookkeeping, pion emission) but l.134 "what remains missing is the DYNAMICAL content of confinement",
l.273 decay is "a condition, not a full dynamical proof of decay, which requires the [dynamics]". Collective-
coordinate dynamics = exactly that missing content -> for the quark it COMPLETES a stated open end, not a bolt-on.
(3) SHARED LYNCHPIN: the reconnection RATE feeds BOTH -- quark/baryogenesis (rate x CP bias x washout = eta's
efficiency factor) AND dark co-genesis (same reconnection in the shared 2I vacuum). ONE build, both sectors. eta
IS a baryon-formation-rate calc; the dark abundance rides on it.
QUARK-SECTOR PAYOFFS: baryon-formation dynamics (2+2->3+1) as a process; pion-emission RATE; the dynamical proof of
Q=1->Q=2 decay XVIII only has as a condition; confinement string-breaking (d_q=1) made dynamical.
CAVEATS: semiclassical (structure of rates, not precise; order-of-mag for violent reconnection); still needs the
rational-map approximate soliton (step a); DON'T conflate with the settled isospin-residual reconnection-barrier
([[reconnection-barrier-redirect]] -- that was a DIRECT-FIELD barrier, wall-blocked, for a DIFFERENT quantity, the
2.5% tower-exponent residual; the collective-coordinate RATE is a different method for eta). SIDESTEPS the static-
resolution wall (collective coords on approximate solitons, not full-field).
VERDICT: a->d is well-motivated on BOTH sectors. The reconnection rate = one calc serving baryogenesis-eta + dark
co-genesis, filling XVIII's acknowledged missing dynamical content. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (bb) -- STEP (a) DONE: rational-map IC VALIDATED (smooth, Q_H converges) -- fixes the tube-lift Q_H=0.
src_paper21/rational_map_ic.py: smooth rational-map Hopfion u=a^p/b^q -> doublet (b^q,a^p) -> n_hat; Q_H via the
validated Whitehead. RESULT: Q_H converges MONOTONICALLY to the integer with N: Q=1 (1,1) -> 0.936/0.960/0.970/
0.977 at N=40/56/72/96; Q=3 (1,3) -> 2.509/2.711/2.800/2.863. Fast (0.7s @ N=96). WHY it works where the tube-lift
gave 0: charge DISTRIBUTED -- 50% of |grad n|^2 spread over ~4-6% of the grid vs the tube-lift's 0.5% (single-
strand -> Q=0). This is the smooth, well-resolved IC the saddle/dynamics need.
SCOPING (honest): validated object = the SMOOTH AXIALLY-SYMMETRIC A_{p,q} Hopfion (core = wound loop), Q_H=3 = a
legit IC in the Q_H=3 SECTOR. The framework's quark is the genuine TREFOIL KNOT T(2,3) = a MODULATED rational map
(Sutcliffe 2013) -- the remaining refinement. Matters because the framework's Q_H=3 minimizer is the KNOT (non-
standard; standard Faddeev knots start at Q=7), so a trefoil-shaped IC biases the saddle toward it, not the axial
loop. STATUS: rational-map METHOD + smooth convergent-Q_H IC = DONE (the tube-lift fix); knotted-trefoil modulation
= next sub-step; then (d) collective-coordinate dynamics. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (cc) -- STEP (d) first milestone: product-ansatz numerics HARD; but robust result = EXOTHERMIC channel + eta redirect.
Built src_paper22/reconnection_product_ansatz.py (two rational-map Hopfions, product ansatz, E(R), Q_H(R)). It hit
GENUINE documented-hard subtleties, not a quick bug: (1) the grid Hopf charge is DOUBLET-GAUGE-dependent -- the
(b^q,a^p) convention gives correct single-soliton Q_H (validated +1 in (a)) but far-field (-1)^q != identity, so
can't feed the SU(2) product (needs vacuum=identity); the identity-vacuum (1,u) convention combines but THROWS
AWAY b^q's phase -> corrupts Q_H (gives -0.68 for Q=1). The two requirements CONFLICT on a finite grid. (2) the
Hopfion product ansatz is itself imperfect (literature spends papers on it). (3) box/resolution tradeoff for
separated solitons. => a VALIDATED interaction-energy calc needs a GAUGE-INVARIANT n_hat-based Hopf charge (solve
curl A = f from the solid-angle 2-form, not the doublet) + a careful ansatz = real research-numerics, NOT a quick
milestone. Do NOT claim it works -- it doesn't yet.
ROBUST RESULT (no hard numerics): by VK/Lin-Yang E ~ C Q^{3/4} (clustering), fusion 2+1->3 is EXOTHERMIC:
E(2)+E(1)=C(2^0.75+1)=2.68C > E(3)=C 3^0.75=2.28C (dE~-0.4C). The baryon-formation channel EXISTS and is DOWNHILL
(knots beat separated multisolitons). First milestone answered via energetics.
ETA REDIRECT (the useful part): fusion is DOWNHILL/not barriered -> the eta efficiency (~1e-9) CANNOT come from a
reconnection barrier; it comes from COSMOLOGICAL WASHOUT (reconnection rate Gamma vs Hubble H at the QCD epoch --
how much CP-biased asymmetry survives expansion). So eta-efficiency = a cosmological Boltzmann question (Gamma vs
H), NOT a soliton-barrier question -> more tractable (no thin-tube wall, no product ansatz) and the physically
correct home for the ~1e-9. NEXT (if pursued): the Gamma-vs-H washout estimate, OR a gauge-invariant n_hat-charge
product-ansatz if the interaction geometry is wanted. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-05 (dd) -- Gamma-vs-H WASHOUT ESTIMATE: eta is DEEP-EQUILIBRIUM-bracketed over ~22 orders, NOT pinned.
src_paper22/washout_eta.py. Reconnection = strong-interaction rate Gamma = n sigma v (n~T^3, sigma~1/Lambda_QCD^2,
v~1). RESULT: at T_QCD=0.15 GeV, Gamma/H ~ 3e19 = DEEP EQUILIBRIUM (freeze-out T* ~ 1e-20 GeV, i.e. NEVER in any
relevant epoch). => eta is BRACKETED over ~22 orders: no-washout (transition-frozen) eta~eps~0.1 (1.6e8x too big);
full strong-washout eta~eps(n/s)(H/Gamma)~3e-23 (5e-14x too small). Observed 6.1e-10 sits BETWEEN but hitting it
needs a SPECIFIC intermediate efficiency (~1e-8) that a fixed-scale Gamma/H does NOT fix -- set by the detailed
confinement-transition dynamics the framework does not specify.
HONEST VERDICT: the eta lynchpin is NOT resolved, and the reason is STRUCTURAL (not effort): this is the SAME
problem that keeps standard baryogenesis OFF the QCD scale -- low T -> deep equilibrium -> no clean out-of-eq
freeze -> eta hypersensitive to the transition. "Threading irreversibility" would have to beat Gamma/H~1e19, which
a rate estimate can't deliver (crossover QCD transition doesn't abruptly shut off washout). So the baryogenesis ->
co-genesis -> 5 GeV chain is STRUCTURALLY COHERENT but QUANTITATIVELY UNPINNED AT ITS ROOT: eta is bracketed, not
predicted, and the dark co-genesis riding on it inherits that. The washout estimate LOCATED eta precisely as a
transition-dynamics-sensitive quantity spanning 22 orders -- not a number that falls out. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-05 (ee) -- CONFIRMED (user): n_hat/solid-angle charge RULED OUT; product-ansatz wouldn't pin eta anyway -> TRANSITION DYNAMICS is the eta direction.
(1) n_hat RULED OUT: the geometric/solid-angle (Berg-Luscher) Hopf-charge estimator was already found WORSE than
FD and non-convergent (notes (z) l.919; O2 record) -- the thin tube makes the director wrap ~O(1)/plaquette ->
solid-angle Omega SATURATES. So update (cc)'s "gauge-invariant n_hat-charge via the solid-angle 2-form" invoked a
ruled-out method. (Nuance: that ruling was for the THIN static tube; on SMOOTH rational-map fields it might not
saturate -- but MOOT, see (2).)
(2) DECISIVE: eta is NOT a reconnection-GEOMETRY quantity -- it's a transition-DYNAMICS quantity (washout (dd):
Gamma/H~1e19 deep equilibrium, eta bracketed over 22 orders, set by WHEN/HOW ABRUPTLY reconnection freezes). So
even a perfect product-ansatz interaction (any charge estimator) would NOT pin eta. product-ansatz/n_hat = WRONG
instrument for the lynchpin.
(3) => the SPECIFIED (first-order-like) CONFINEMENT-TRANSITION DYNAMICS is the right and only direction for eta:
eta lives in the OUT-OF-EQUILIBRIUM FREEZE of the transition. If confinement is effectively FIRST-ORDER (bubble
nucleation, latent heat, abrupt freeze), the CP-biased asymmetry freezes BEFORE washout completes -> efficiency
(hence eta) set by nucleation/freeze dynamics, not left bracketed. FRAMEWORK HOOK: DENSITY FEEDBACK could MAKE the
confinement transition first-order -- the (1+beta* rho) softening is a positive-feedback/hysteretic loop on
clustering, exactly what turns a crossover into first-order. SM QCD crossover is NOT out-of-eq enough (why SM
baryogenesis isn't at QCD); the density-feedback version MIGHT be. That is the eta-pinning mechanism to pursue.
[[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-06 (ff) -- FEEDBACK -> FIRST-ORDER confinement transition CONFIRMED (above a threshold). eta-pinning reduces to ONE ratio gamma/b.
src_paper22/firstorder_transition.py: toy Landau F(rho)=a rho + b rho^2 - gamma beta rho^2/(1+beta rho). RESULT:
beta=0 -> CROSSOVER (baseline, feedback off). beta>0 opens a BARRIER (FIRST-ORDER) ABOVE a threshold in (beta,
gamma/b): at framework beta*~0.45, first-order needs gamma/b >~ 4 (gamma=1,2 crossover; gamma=4 barrier 0.027;
gamma=8 -> 2.0; gamma=16 -> 29). beta-threshold drops with gamma (g=1: beta>1.15; g=2: beta>0.55; g=4: beta>0.3).
So density feedback INDUCES a first-order transition once the feedback GAIN gamma exceeds ~4x the knot-core
STABILIZATION b (at beta~0.45); below -> crossover; above -> robustly out-of-equilibrium (barrier grows fast).
WHY THIS MATTERS: FIRST direction in the eta chase that is (a) VIABLE, (b) where eta ACTUALLY LIVES (first-order ->
supercool -> abrupt freeze = the out-of-eq the washout (dd) showed is the decider), (c) NOT resolution-wall-blocked
-- it reduces to ONE ratio gamma/b (feedback-gain : stabilization), read ALGEBRAICALLY from the framework's own
E_fb = K_fb J4 (K_fb = J2a + mu J2iso_fb). => converts "eta bracketed over 22 orders" into a SHARP question: is the
framework's gamma/b >~ 4? YES -> first-order -> eta pinnable by nucleation/freeze dynamics; NO -> crossover -> eta
stays bracketed, QCD-scale baryogenesis genuinely fails (clean negative). CAVEATS: toy Landau (mechanism, not the
exact coefficients); gamma/b and the beta-mapping must come from the real E_fb; even first-order here is initially
WEAK near threshold (modest supercooling) -> strength matters for how well eta is pinned. NEXT: extract gamma/b
from E_fb. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-06 (gg) -- RESOLUTION (user): feedback strength = beta*·rho, and rho was ~10^47x larger early -> transition DEEPLY first-order at QCD.
User: "hopfions are stable at beta~0.452 NOW, but beta was different early." CORRECT, refined: the feedback factor
is (1+beta*·rho) -- the STRENGTH is beta*·rho, NOT beta* alone. rho (condensate ~ T^4; Lcond~T, claude-hopfion s5)
was rho_QCD/rho_today ~ T_QCD^4/(10 Lcond^4) ~ (150 MeV)^4/(10*(1.19e-4 eV)^4) ~ 2.5e47 (~10^47) larger at the QCD
epoch than today. So beta_eff = beta*·rho was ~10^47x STRONGER early -> DEEPLY SATURATED feedback -> STRONGLY
FIRST-ORDER confinement transition -> robustly out-of-equilibrium -> the eta-pinning mechanism was ON exactly when
baryons formed. TODAY beta*·rho small -> unsaturated -> CROSSOVER -> Hopfions stable (the beta*~0.452 value).
RESOLVES the (ff) worry: at beta*~0.45 the TOY gave crossover, but that used TODAY's strength; the transition
happened at QCD where beta*·rho was ~10^47 above any plausible first-order threshold -> deeply first-order WHEN IT
MATTERED, out-of-eq GENERIC (not tuned). Self-consistent: baryogenesis/confinement early (strong feedback, first-
order) + stable baryons today (weak feedback, crossover). CAVEATS: 10^47 ratio robust; the exact threshold mapping
is toy (needs real E_fb) but QCD swamps it by dozens of orders. DEEPLY first-order != automatically eta=6e-10:
strong first-order guarantees STRONG out-of-eq (escapes the 22-order washout bracket of (dd)) but the EXACT eta
needs the nucleation/supercooling/freeze calc. VIABILITY of the transition-dynamics eta-route = RESOLVED (mechanism
was on); the NUMBER is the quantitative follow-up. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-06 (hh) -- CORRECTION to (gg): first-order is an intermediate WINDOW; huge beta_eff OVER-saturates to crossover. gamma/b extraction STILL needed (now the crux).
User asked "do we still need gamma/b given (gg)?" -> checking the toy shows (gg) OVER-READ. firstorder_transition
beta-scan (b=gamma=1): beta<~1 crossover; beta~2-1e4 FIRST-ORDER (barrier peaks ~beta=10 at 0.0095); beta>~1e6
CROSSOVER AGAIN (barrier->0). First-order is an INTERMEDIATE-strength WINDOW, NOT monotonic. PHYSICAL: deep
saturation 1/(1+beta* rho)->0 = stiffness softened to ZERO -> no barrier -> continuous/crossover. So (gg)'s
"beta* rho~1e47 -> deeply first-order" is WRONG: 1e47 is OVER-saturated -> CROSSOVER. Feedback was strong early
(true) but "stronger" OVERSHOOTS the window, doesn't land in it.
CONSEQUENCE: first-order requires beta_eff IN THE WINDOW at the transition. Universe cools through T_QCD with
beta_eff=beta*·rho(T) DECREASING, so the question = where does beta_eff(T_QCD) sit vs the window? Needs: (1) window
LOCATION in physical units = the gamma/b ratio from real E_fb [STILL NEEDED, now the crux]; (2) beta_eff(T_QCD)=
beta*·rho(T_QCD) in those units [where the transition sits]. Together -> first-order vs crossover at QCD. (gg)
established only that feedback was STRONG early, NOT first-order (likely past the window, over-saturated).
NET: ANSWER TO USER = YES, gamma/b extraction still needed (the crux). The high early density pushes toward the
UPPER (saturated/crossover) edge, so first-order is NOT guaranteed -- it needs beta_eff(T_QCD) to land in the
intermediate window. NEXT: extract gamma/b (window location) + estimate beta_eff(T_QCD). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (ii) -- gamma/b extraction does NOT cleanly close: the framework lacks the BULK confinement free energy. eta chase terminus.
Attempted the gamma/b extraction (hh's crux). It does NOT cleanly map from E_fb, for a structural reason: TOY =
bulk Landau F(rho)=a rho+b rho^2-(feedback), explicit stabilization b, first-order in an intermediate WINDOW.
FRAMEWORK E_fb=K_fb*J4 = PRODUCT, scale-free, NO bulk stabilization (scale fixed by self-consistency sopt=1), and
a SINGLE-KNOT functional -- NOT the bulk confinement free energy a phase transition needs. The algebraic ratios
that ARE specified (K_fb=2phi J2a, J2iso_fb=phi J2a, J4=0.227 J2a, Paper I) are single-soliton, not bulk-thermo;
the unsoftened J2iso (needed for the feedback GAIN) isn't cleanly available; the condensate potential V(rho)=
Lam0 e^{-phi^6 rho} the framework gives is the DARK-ENERGY one, not the confinement free energy.
CONCLUSION (honest terminus of the eta chase): the first-order-transition route is VIABLE IN PRINCIPLE (mechanism
real: feedback -> first-order in an intermediate window) but pinning whether the framework's confinement transition
SITS in that window needs the BULK CONFINEMENT THERMODYNAMICS (free energy + nucleation), which the framework does
NOT specify (it gives a single-knot energy + a cosmological condensate potential, neither = the confinement
transition free energy). So eta stays UNPINNED, obstruction precisely located: NOT the resolution wall (static
saddle), NOT the washout bracket per se (first-order could escape it), but a MISSING FRAMEWORK INPUT -- the bulk
confinement free energy/nucleation dynamics. This is a "needs new framework development" boundary, not a grind-able
computation. ETA CHASE TERMINATES HERE honestly. [[dark-hopfion-cdm-lead]], [[o2-feedback-saddle-solver]].

## UPDATE 2026-10-06 (jj) -- AUDIT (user): did early-universe calcs wrongly use TODAY's T_CMB? NO error; surfaces a useful consistency.
User flagged: the framework input is T_CMB(today); did we misuse it for the early QCD universe? AUDIT: NO concrete
result corrupted -- where the epoch mattered we used EARLY scales/scalings: Gates used rho(z)/rho_today~(1+z)^4;
washout (dd) used T_QCD=0.15 GeV for H,Gamma and Lambda_QCD~0.2 GeV for sigma (NOT today's Lcond); rho-ratio (gg)
used T_QCD^4 on the early side, Lcond(today) only on the TODAY side; m_dark~5 GeV = ADM ratio x m_p (today budget).
T_CMB(today) only anchors (a) TODAY's Lcond and (b) the DIMENSIONLESS constants (phi-powers, alpha=360/phi^2,
beta*, mu=3-phi, McKay/WZW) = GEOMETRIC/topological -> T-INDEPENDENT, and the 2I structure is PRESERVED at all
epochs (uniform scaling; T/Lcond~2 ALWAYS, near-critical) -> legit at QCD.
GENUINE CONSISTENCY SURFACED: Lcond TRACKS T (Lcond=0.506 T), so Lcond(QCD)=0.506*150MeV ~ 76 MeV ~ T_QCD ~
Lambda_QCD -- the condensate scale, temperature, and CONFINEMENT scale ALL COINCIDE at the QCD epoch (vs today's
Lcond~1e-4 eV, 12 orders smaller). Self-consistent: confinement happens WHEN Lcond reaches the QCD scale. Using
today's Lcond for an early confinement calc would be a ~12-order error -- AVOIDED by scaling.
SHARPENS THE ETA TERMINUS: the MISSING bulk confinement free energy ((ii)) has its SCALE at Lcond(QCD)~76 MeV ~
T_QCD ~ Lambda_QCD (coinciding), NOT today's scale; and near-critical (T/Lcond~2) is exactly where first-order-vs-
crossover is delicate (consistent with the "window" finding (hh)). So the user's point LOCATES eta's missing input
more precisely (a near-critical confinement free energy at Lcond~76 MeV), doesn't overturn results. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (kk) -- RESEARCH+BRAINSTORM: the bulk confinement free energy = a Z3 Polyakov-loop Landau-Ginzburg; trefoil Z3 FAVORS first-order.
Lit (Polyakov-loop effective potential): the deconfinement free energy is a Z3-symmetric Landau-Ginzburg for the
order parameter L: V(L)=(b2(T)/2)|L|^2 - (b3/3)(L^3+c.c.) + (b4/4)|L|^4. The CUBIC b3 (allowed ONLY by Z3) makes
pure-gauge SU(3) FIRST-ORDER; light dynamical quarks break Z3 -> crossover (Columbia plot). So first-order-vs-
crossover = how intact the Z3 is.
FRAMEWORK MAP (apt): the TREFOIL has an exact Z3 (three colours/crossings = Z3 outer automorphism of 2T, P15/XVI)
= the analog of the QCD CENTER Z3. => the framework's bulk confinement free energy = a Z3 Landau-Ginzburg in a
trefoil-condensation order parameter, trefoil-Z3 -> cubic -> generically FIRST-ORDER. AND plausibly stays first-
order where real QCD is crossover: real QCD crossover BECAUSE light dynamical quarks break Z3; in the framework
the "quarks" are the STRANDS OF THE TREFOIL (not separate light fields) -> Z3 carried by the knot, far more INTACT
-> first-order favoured. Density feedback (1+beta* rho) adds to b2/b3/b4 (toy) -> 2nd independent first-order driver.
SCALE: coefficients at Lcond(QCD)~76 MeV, near-critical (T/Lcond~2) -- the T_CMB-audit (jj) regime.
NUCLEATION/eta: standard first-order -- bubble nucleation (bounce action, supercooling T_nuc<T_c), CP-biased
RECONNECTION at the BUBBLE WALLS (EW-baryogenesis-style at the confinement scale). eta ~ eps x wall-efficiency,
computable once V(L) fixed. CAVEAT: strong first-order + FAST walls SUPPRESSES baryogenesis -> OPTIMAL wall speed
(first-order necessary NOT sufficient for eta=6e-10).
NET: converts the missing input from "unspecified" to a CONCRETE TRACTABLE TEMPLATE -- derive b2,b3,b4 (+feedback)
from the trefoil/2T structure, Z3 Landau-Ginzburg at ~76 MeV, bubble-nucleation baryogenesis at the walls. CAVEATS:
Polyakov loop is an ANALOGY (framework not a gauge theory; L=trefoil-condensation, Z3=trefoil colour, not a Wilson
holonomy); coefficients still need framework derivation; efficiency has an optimal-wall-speed sweet spot.
[[dark-hopfion-cdm-lead]]. Sources: Polyakov-loop effective potentials (hep-ph/0011193; arXiv:2405.17335; BNL
Polyakov modeling); first-order QCD baryogenesis/nucleation (arXiv:2407.17549; hep-ph/0304275; arXiv:2003.08892).

## UPDATE 2026-10-06 (ll) -- FIRST-ORDER CHECK RUN (firstorder_z3.py): framework confinement = FIRST-ORDER (favored), weak/marginal. Proposition -> checked.
V(l)=(b2/2)l^2-(2/3)l^3+(1/4)l^4-2h*l (Z3 Polyakov-loop, b3=b4=1 normalized, h=explicit Z3 breaking). RESULTS:
pure-gauge h=0 FIRST-ORDER (barrier 0.050, jump 1.33); framework h~0.02-0.14 FIRST-ORDER (barrier 0.0005-0.03);
h=0.3 crossover; real QCD h~1 crossover. Critical endpoint h_c~0.21 (first-order below).
DERIVED (firm): cubic b3!=0 -- trefoil 3 colours -> Z3 cubic 3x3x3->1 (colour-singlet baryon) = the QCD-center
cubic. First-order DRIVER present by construction.
KEY framework-vs-QCD difference (robust): explicit Z3 breaking h is SMALL because confinement forbids FREE LIGHT
quarks -- only STRING-BREAKING breaks Z3, costing ~2 m_const ~600 MeV >> T_QCD~150 -> h~e^{-600/150}~0.02 (string-
breaking) to ~0.14 (single constituent), WELL BELOW h_c~0.21. Real QCD crossover ONLY because free light quarks
give h~1; the framework has NONE (quarks = confined trefoil strands) -> lands FIRST-ORDER. => the proposition is
now a CHECKED result: framework confinement is first-order, differing from QCD by the absence of free light quarks.
CAVEATS: coefficients NORMALIZED not fully derived (Z3 structure/cubic derived; exact h_c needs real b3/b4 from the
framework thermal free energy) -> "first-order FAVORED" not "certain"; WEAK first-order (barrier << pure-gauge, near
endpoint); b2(T) is a temperature proxy.
ETA UPGRADE: a first-order transition (even weak) gives OUT-OF-EQUILIBRIUM at bubble walls -> eta ESCAPES the
22-order deep-equilibrium washout bracket (dd) -- a genuine freeze-out exists. BUT weak first-order = modest out-of-
eq, often INSUFFICIENT to fully suppress washout -> eta is "favored nonzero + partially pinned", NOT robustly 6e-10.
Weakness is the limiting factor. NEXT if pursued: derive b3/b4 (real h_c) from the 2T thermal free energy; the
bubble-wall efficiency for eta. [[dark-hopfion-cdm-lead]]. Scripts: src_paper22/firstorder_z3.py.

## UPDATE 2026-10-06 (mm) -- BANKED: first-order confinement FAVORED. Next = derive b3/b4 from the THERMAL FREE ENERGY OF THE EPOCH'S GEOMETRY (user).
Banking (ll): framework confinement transition = FIRST-ORDER (favored, weak/marginal) via Z3 cubic (trefoil 3
colours) + small explicit breaking (no free light quarks, only expensive string-breaking) -> h<h_c. eta escapes
the 22-order washout bracket; weak first-order -> partially pinned, not robustly 6e-10.
NEXT (user's directive + key insight): derive the EXACT b3/b4 (hence h_c, first-order strength, eta efficiency)
from the condensate THERMAL FREE ENERGY -- but "of the GEOMETRY OF THE EPOCH": the thermal FE must be computed for
the QCD-epoch condensate state (Lcond(QCD)~76 MeV, 2I structure preserved, near-critical T/Lcond~2, trefoils
forming), NOT today's Lcond~1e-4 eV. CHAIN: eta <- b3/b4 <- thermal FE F(T,L) <- the epoch's condensate geometry.
APPROACH: F(T,L) = T INT d^3k/(2pi)^3 ln(1-e^{-E(k;L)/T}) over the SEMI-DIRAC excitations (E^2=A^2k^4 sin^4 th +
v^2 k^2 cos^2 th, claude-hopfion s3), with an L-dependent gap (confined gapped, deconfined gapless) and the Z3
(three-colour) structure giving the cubic b3. Expand in L -> b2(T),b3,b4. Scale set by Lcond(QCD)~76 MeV.
[[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (nn) -- THERMAL FREE ENERGY OF THE EPOCH'S GEOMETRY DERIVED: semi-Dirac -> T^3, CUBIC-Bernoulli (B3) Weiss potential.
src_paper22/semidirac_thermal.py. The QCD-epoch condensate excitations have the SEMI-DIRAC dispersion
E^2=A^2k^4 sin^4 th + v^2 k^2 cos^2 th (claude-hopfion s3), NOT relativistic. Computed the thermal FE f(T)=
T INT d^3k/(2pi)^3 ln(1-e^{-E/T}) and its holonomy (Polyakov-phase) dependence. RESULTS (clean numerics):
(1) f(T) ~ T^p, FITTED p=3.000 -> f ~ T^3/(vA) [vs relativistic T^4] -- the semi-Dirac anisotropy (k^4 hard + k
    easy) reduces effective phase space.
(2) holonomy Fourier coeffs c_n = 1/(pi n^3) EXACTLY (c_n*n^3 = 0.31831 = 1/pi for every n) -> f(phi) =
    -(T^3/pi) sum cos(n phi)/n^3 = a CUBIC Bernoulli B3(phi) [vs relativistic 1/n^4 -> quartic B4].
=> the epoch-geometry Weiss potential is LOWER-DEGREE, CUBIC-DOMINATED (B3) vs standard-QCD (B4). The cubic IS the
first-order driver -> the semi-Dirac epoch geometry PLAUSIBLY ENHANCES first-order vs a relativistic condensate.
STATUS: thermal FE OF THE EPOCH'S GEOMETRY = DERIVED in structure (T^3/B3, framework-specific, firm numerics) --
the user's flagged prerequisite MET. INTERPRETATION "B3 -> stronger first-order" plausible, needs the b3/b4 mapping.
REMAINING for exact b3/b4 -> eta: SU(3) 3-COLOR SUM (3 eigenvalue phases -> Z3 order param L) + GAP(L) relation
(confinement gaps the semi-Dirac modes) + explicit breaking h -> h_c, first-order strength, bubble-wall efficiency.
Chain to eta now ONE well-defined calc away (color-sum+gap), not open-ended. (A,v at epoch set the SCALE; the B3
SHAPE is dispersion-fixed.) [[dark-hopfion-cdm-lead]]. Papers: 21 (structure), 22 (genesis dynamics) drafted.

## UPDATE 2026-10-06 (oo) -- SU(3) COLOUR-SUM done: b3/b4 DERIVED (~0.17), semi-Dirac MODESTLY enhances first-order. Paper scoping decided.
PAPER SCOPING (user): KEEP SEPARATE -- 21 = firm STRUCTURAL sector (no-E7, sterile Dirac); 22 = conditional
GENESIS BRIDGE (quark sector -> DM via co-genesis at the first-order confinement transition). Reason: epistemic
cleanliness (firm vs open) + the bridge is a distinct contribution. DONE: 21 now DEFERS genesis to 22 (cross-ref +
bibitem added, main_paper21.tex sec:abundance); Paper 22 drafted (main_paper22.tex).
SU(3) COLOUR-SUM (su3_colorsum_coeffs.py): single-phase holonomy f(phi)=-sum c_n cos(n phi), c_n=1/n^s (s=3 semi-
Dirac B3 / s=4 relativistic B4); SU(3) real-L slice phases (phi,0,-phi), L=(1+2cos phi)/3; V(L)=const+2 sum c_n
[cos(n phi)-1], expand -> b2,b3,b4. RESULT: relativistic b3/b4=0.159, R=b3^2/(b2 b4)=0.29; SEMI-DIRAC b3/b4=0.169,
R=0.50. => b3!=0 (cubic/first-order driver DERIVED, from trefoil Z3); semi-Dirac epoch geometry ENHANCES first-
order MODESTLY (b3/b4 ~6% larger; R ~1.7x). First-order realized when b2(T)<b3^2/b4 (as pure-gauge SU(3), similar
b3/b4, IS first-order) -> weakly-to-moderately first-order, consistent with (ll), now with DERIVED coefficients.
STILL NEEDED for eta (correcting my "one calc" over-compression): (a) gap(L) -> b2(T) [thermal mass: how
confinement gaps the semi-Dirac modes, drives b2(T) through 0 at T_c -> transition T + first-order STRENGTH];
(b) bubble-wall efficiency -> eta. CHAIN: eta <- wall-eff <- [b2(T) from gap(L)] + [b3/b4 DERIVED] <- thermal FE
(T^3/B3 DERIVED) <- epoch geometry (DERIVED). Two links derived this turn; gap(L)->b2(T) and wall-eff remain.
[[dark-hopfion-cdm-lead]]. Scripts: src_paper21/{semidirac_thermal,su3_colorsum_coeffs}.py.

## UPDATE 2026-10-06 (pp) -- gap(L) -> b2(T) DERIVED: thermal mass drives the transition; T_c ~ 2 Lcond ~ T_QCD (near-critical match).
src_paper22/gap_b2T.py. Z3-respecting gap EVEN in L: Delta(L)^2=Delta0^2(1-L^2) (Delta0=confinement gap~Lcond(QCD);
no linear term). Gapped semi-Dirac modes E^2=A^2k^4sin^4+v^2k^2cos^2+Delta^2. b2_thermal(T)=d^2f/dL^2|_0 =
-2 Delta0^2 (df/d(Delta^2)). RESULT: b2_thermal<0 at all T, |b2_thermal| RISES MONOTONICALLY with T (4e-4 at
T/Delta0=0.3 -> 0.15 at T/Delta0=3). So b2(T)=2a_vac-|b2_thermal(T)| goes +->- : L=0 CONFINED stable at low T,
DECONFINED (unstable) at high T -- THERMAL MASS DRIVES THE TRANSITION (derived). Sign change at T_c~2.3 Delta0
(a_vac~0.05). With Delta0~Lcond(QCD)~76 MeV: T_c ~ 2 Lcond ~ 150 MeV ~ T_QCD -- AND T_c~2Lcond EXACTLY matches the
framework near-critical condition T/Lcond~2 (claude-hopfion s9). Not an input, a CONSISTENCY (lands on QCD scale +
near-criticality).
(Fixed: first tried Delta(L)=Delta0(1-L) LINEAR -> gave a spurious linear/explicit-breaking term + non-monotonic
f''; the physical Z3-respecting gap is EVEN, Delta^2~(1-L^2).)
CHAIN: epoch geometry [v] -> thermal FE T^3/B3 [v] -> b3/b4~0.17 [v] -> b2(T) drives transition, T_c~2Lcond [v,
this turn]. Confinement transition FIRST-ORDER (b3!=0) at the QCD scale -- derived. REMAINING: (1) NORMALIZATION
reconciliation -> first-order STRENGTH: b2(T)[gap, Delta0 scale] and b3/b4[color-sum, Weiss T^3 scale] are in
DIFFERENT normalizations; both from the SAME semi-Dirac FE -> compute the FULL F(T,L) (holonomy+gap together) in
ONE calc -> consistent b2,b3,b4 -> barrier/latent heat at T_c. (2) bubble-wall efficiency -> eta. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (qq) -- UNIFIED F(T,L) done: the transition is STRONGLY first-order (one normalization). Normalization step CLOSED.
src_paper22/unified_FTL.py. The full free energy in ONE normalization (Matsubara/holonomy sum):
  F(T,L) = -T sum_n (1/n) I_n(T,Delta(L)) [2 cos(n phi)+1] + a_vac L^2,  I_n=INT d^3k/(2pi)^3 e^{-nE/T},
  E=sqrt(semi-Dirac^2 + Delta^2), SU(3) real-L phases (phi,0,-phi) with cos(n phi)=T_n((3L-1)/2), gap EVEN
  Delta(L)^2=Delta0^2(1-L^2), a_vac L^2 = T=0 confinement cost. Holonomy+gap TOGETHER = exactly the b2/b3/b4
  normalization reconciliation that (pp) flagged. Fast fixed-grid (Gauss-Legendre) vectorised integral.
RESULT -- the GLOBAL minimum of F(L) JUMPS DISCONTINUOUSLY from the confined branch (small L) to the deconfined
well (L=1) at T_c (true equal-depth coexistence), a LARGE order-parameter jump => STRONGLY first-order:
    a_vac   T_c/Delta0   L_conf   jump dL=1-L_conf   barrier/max|F|
    0.1       0.92        0.348        0.652              0.95
    0.2       1.17        0.406        0.594              0.54
    0.4       1.48        0.454        0.546              0.35
So dL~0.55-0.65 (LARGE) and barrier a substantial fraction of the FE scale => genuinely (NOT merely weakly)
first-order. Driven by the Z3 trefoil cubic (n=3k holonomy harmonics = 3-fold colour structure), as in standard
pure-gauge SU(3). T_c/Delta0 ~ O(1-1.5) ~ O(Lcond(QCD)) (coexistence T_c; the (pp) spinodal T_c~2Delta0 is the
higher curvature-sign-change estimate -- ordering consistent, both O(1-2)Delta0 ~ QCD scale).
KEY RECONCILIATION (two earlier facts now consistent): the quartic-Landau b4 fit went NEGATIVE/ill-defined
(mid-session) PRECISELY BECAUSE the jump is large -- a strongly first-order transition lives BEYOND the quartic
truncation's validity, so the honest instrument is the FULL F(T,L), not a local b2/b3/b4 fit. The (ll) "weak/
marginal" verdict was the Landau-truncation under-reading; the gap coupling + full holonomy harmonics make it
strong. [[o2-e6-cft-tower]].
CAVEATS (claude-hopfion s6/s10): (i) TOY semi-analytic matrix+Landau model; A=v=Delta0=1 units; a_vac (T=0
confinement/bag scale) is PARAMETRIZED not derived (cutoff-set magnitude -- flag). (ii) confined min sits at
L_conf~0.35-0.45 NOT 0 (incomplete center symmetry in the toy) -> the PRECISE dL/strength are indicative, the
first-order CHARACTER + O(1) jump are robust. (iii) strength DECREASES with a_vac (0.95->0.35) -> "strongly"
clearest at small a_vac; even at a_vac=0.4 strength~0.35 is substantial.
CHAIN now: epoch geometry [v] -> thermal FE T^3/B3 [v] -> b3/b4~0.17 [v] -> b2(T) drives it, T_c~O(Lcond) [v] ->
UNIFIED F(T,L): STRONGLY first-order, dL~0.6, barrier/FE~0.5 [v, this turn]. Out-of-equilibrium Sakharov condition
STRONGLY MET -> large departure from equilibrium at the bubble wall (good for baryo/co-genesis; escapes the
washout equilibrium of (ii)). REMAINING for eta: ONE step -- bubble-wall velocity/efficiency x the derived WZW
delta_CP -> the actual baryon + dark asymmetry. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (rr) -- FINAL eta step done (wall_eta.py): all 3 Sakharov conditions DERIVED/MET; eta magnitude = a TRANSPORT number; co-genesis ratio is the falsifiable output.
src_paper22/wall_eta.py. Folds the DERIVED strong first-order transition (qq) + the DERIVED WZW CP phase into eta.
CP INPUT (actual numbers, Paper V P5:thm:delta_cp, NOT from memory): delta_CP = 65.0 deg (LO 17pi/40 = 76.5 deg;
subleading -> 65.0 deg, within 0.12sigma of PDG 65.4+-3.3). sin(delta_CP)=0.906 = O(1) => CP is NOT the bottleneck
(confirms washout (v)).
LATENT HEAT / alpha extracted FROM F(T,L) (entropy s=-dF/dT across T_c, per branch):
    a_vac   T_c/Delta0   L_heat     rho_rad    alpha    v_w(deflag..deton)
    0.1       0.92       2.52e-1    4.01e0     0.063    [0.20, 0.99]
    0.2       1.17       4.79e-1    1.05e1     0.046    [0.18, 0.99]
    0.4       1.48       9.16e-1    2.74e1     0.033    [0.15, 0.99]
NUANCE (reconciles (qq)): the transition is STRONG in ORDER-PARAMETER space (dL~0.6, big L-jump) but ENERGETICALLY
MILD (alpha = L_heat/rho_rad ~ 0.03-0.06, i.e. latent heat only ~3-6% of the radiation bath) -- because the FE
barrier scale (~0.01-0.07 Delta0^4) is small next to g* T^4. Both true, different measures; alpha~0.05 is typical
for a QCD-scale transition. v_w ~ 0.15-0.2 (deflagration) up to ~1 (detonation branch). The departure from
equilibrium is genuine; what matters for baryogenesis is washout-avoidance (below), not alpha size.
WASHOUT-AVOIDANCE IS AUTOMATIC + STRONG here: confinement itself removes the free strands that B-violating
reconnection needs, so the "sphaleron shutoff in the broken phase" needs NO v/T>1 tuning (unlike EWBG) -- the
broken(confined) phase simply cannot reconnect. => eta is NOT Sakharov-blocked.
eta MAGNITUDE: eta_B = 7.04 n_B/s, n_B/s ~ A_pref eps_CP kappa_transport, A_pref=15/(4pi^2 g*)~2.2e-2,
eps_CP=sin(delta_CP) c_int. eta_max(kappa=1) ~ 1.4e-2 .. 1.4e-1 (OVER-produces vs eta_obs=6.1e-10). Matching eta_obs
needs net transport efficiency kappa_transport ~ 4e-9 .. 4e-8. ** This is the ONE remaining non-derived factor, now
cleanly a TRANSPORT/branching number (full reconnection Boltzmann + diffusion), NOT a missing Sakharov condition.**
Exactly the obstacle (ii)/(v) located -- now ISOLATED and quantified, with both the CP input (65 deg) and the
out-of-equilibrium condition (strong 1st-order) DERIVED, so the gap is strictly the efficiency.
ROBUST FALSIFIABLE OUTPUT (kappa-INDEPENDENT -- this is THE prediction, not eta): baryon + dark asymmetries share
the SAME delta_CP and freeze at the SAME transition (co-genesis) -> Omega_dark/Omega_b set by the 2D_3:2T subgroup
branching (|2D_3|:|2T| = 12:24) ~ O(5), near the observed ~5.4. eta's absolute value stays a transport calc; the
RATIO is the framework output. [[dark-hopfion-cdm-lead]] [[o2-e6-cft-tower]].
STATUS: the Paper 22 (co-genesis) program is now internally complete at the STRUCTURAL level -- all Sakharov
conditions met, CP derived, transition derived, ratio predicted; the single open quantity (eta magnitude) is a
transport efficiency, flagged [O]. This is the honest terminus of the eta chase (do not re-walk as a Sakharov
problem).

## UPDATE 2026-10-06 (ss) -- ANATOMY of the transport efficiency kappa: it FACTORISES into a DERIVED departure-from-equilibrium x ONE open branching ratio. Exponential-washout sensitivity GONE.
src_paper22/transport_eta.py. Opened up kappa (the sole open factor in eta=7.04 A_pref eps_CP kappa). PHYSICS: Gamma/H
~3e19 at T_QCD = DEEP equilibrium; by CPT n_B=0 in EXACT equilibrium regardless of C,CP -> ONLY the departure at the
wall survives. kappa factorises:
  kappa = delta_neq * b_B ,
  delta_neq = v_w/(Gamma L_w) = 1/(Gamma tau_wall)   [DERIVED]  -- leading departure from equilibrium; the system
    re-equilibrates at Gamma, the wall sweeps a point in tau_wall=L_w/v_w, and Gamma*tau_wall~30-370>>1 so the
    departure is SMALL and POWER-LAW. With v_w~0.15-0.2 (from alpha~0.05, wall_eta) and L_w=c_w/T, c_w in [2,20]
    (dimensionless wall thickness, cutoff-set via the L-field stiffness -- FLAG, claude-hopfion s10):
    delta_neq in [2.7e-3, 3.6e-2]. This CAPS kappa <~ 1e-2 -> eta CANNOT be O(1): the "no-washout overproduces"
    EDGE of the old 22-order bracket (washout note (cc)/(dd)) is resolved WITHOUT tuning.
  b_B = baryon-number-VIOLATING branching: fraction of CP-biased reconnections that are Q=2+2->3+1 threadings among
    ALL condensate reconnections, x CP-interference c_int.   [OPEN, interpretive]
RESULT: with the DERIVED prefactors (A_pref=0.022, sin delta_CP=0.906, delta_neq~1e-2), observed eta=6.1e-10 needs
  b_B * c_int ~ 3e-7 .. 8e-7.
So eta's magnitude REDUCES TO ONE PHYSICAL RATIO: a B-violating branching ~1e-6..1e-7. [I] plausible (B-violating
threadings ARE rare vs generic condensate reconnections) but NOT derived -- needs the reconnection Boltzmann + the
Q-selection rule on the 2I vertex. NB this ~1e-7 is DEFINITIONAL (b_B := whatever is left); the physics content is
the IDENTIFICATION of the leftover with a branching ratio + the argument that ~1e-7 is a plausible magnitude, NOT a
second derivation of eta (do not cite it as independent corroboration -- CLAUDE.md s3).
KEY STRUCTURAL GAIN: the old bracket spanned ~22 orders because it kept an EXPONENTIAL washout exp(-Gamma tau).
That is REMOVED: confinement shuts the B-violating rate off behind the wall (automatic, (rr)), so the washout
window is capped and the suppression is the POWER-LAW delta_neq, not exp(-Gamma tau). The remaining uncertainty is
a ~1-order window on the DERIVED delta_neq times ONE open branching ratio -- a well-posed (hard) quantity, not a
22-order mystery. Caveat: this is a SCHEMATIC leading-order anatomy (kappa=delta_neq*b_B), not the full transport
calc, which could reshuffle O(1)-O(10) factors. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (tt) -- going after b_B FLIPS the finding: framework factors OVER-produce eta by ~1e4-1e5. The open quantity sharpens to the CROSS-CFT fusion amplitude.
src_paper22/branching_eta.py. Tried to assemble b_B = b_ch * c_int from framework structure; result CORRECTS the
(ss) framing ("need a tiny ~1e-7 branching, plausible" -- WRONG on inspection; CLAUDE.md s9).
(A) CP interference c_int: from Paper V (P5:thm:delta_cp) the CP phase = arg(Omega), Omega=A1+A2, A1=e^{i 17pi/40}
(|A1|=1, Mech-1 topological spin), A2=(1/9)sum_m d_m e^{2pi i m/90} (E8 twisted, m in {1,7,11,13,17,19,23,29},
marks (2,3,4,5,6,4,2,3)). COMPUTED: |A2|=2.85 (so the two contributions are COMPARABLE -- Mech-2 is "subleading"
only as a PHASE correction, not in modulus), relative weak phase -15.5 deg, interference VISIBILITY 2r2/(1+r2^2)=
0.63. => c_int ~ 0.6 x absorptive(O(0.1-1)) = NOT loop-suppressed. The large delta_CP (65 deg) carries straight
through to a large per-event CP asymmetry. [Caveat: mapping Paper V's A1,A2 to the two baryogenesis amplitudes is
an ASSUMPTION; robust claim = the two contributions are comparable so no small CP factor exists.]
(B) topological channel fraction b_ch: Hopf-conserving outputs of 2+2 (Q=4) = {4+0 gauge boson, 3+1 BARYON+nu, 2+2
elastic}; energetics E=C Q^{3/4} (dE: 4+0=-0.54C most exothermic, 3+1=-0.08C, 2+2=0) -> crude b_ch~0.20 (baryon is
1 of 3, not dominant).
(C) ASSEMBLED: b_B = b_ch*c_int ~ 1e-2..1e-1 (NOT 1e-7). eta = 7.04 A_pref (sin dCP c_int)(delta_neq b_ch) ~
1e-5..1e-4 -> OVER-PRODUCES observed 6.1e-10 by ~1e4-1e5.
** THE FINDING FLIPS: reconnection baryogenesis here is EFFICIENT; eta's problem is a SUPPRESSION ~1e-6, NOT a
smallness. ** Localisation: the suppression is NOT in CP (c_int~0.6 computed) and NOT in delta_neq (~1e-2 derived)
-- the ONLY place left is b_ch, specifically the CROSS-CFT FUSION amplitude. Lepton=SU(2)_3, baryon=(E_6)_1 are
DIFFERENT CFTs, so lepton+lepton->baryon is a cross-sector fusion; the crude energetic b_ch~0.2 is BLIND to the
WZW/T-matrix selection rule governing it. That selection rule is the SAME machinery as Q_group + the D-type
darkness closure (dtype/e6_cft_tower_engine) => FRAMEWORK-COMPUTABLE. If the cross-sector amplitude is strongly
suppressed (~1e-7) eta lands; if not, the framework over-produces. So the open quantity sharpens from "a transport
efficiency" to "the SU(2)_3 -> (E_6)_1 reconnection fusion amplitude" -- a WZW/T-matrix number. NEXT if pursued:
compute that cross-CFT fusion amplitude (branching of two SU(2)_3 j primaries into the E_6 baryon channel) with the
T-matrix engine. [[o2-e6-cft-tower]] [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (uu) -- USER CATCH: branching_eta ISOTROPISED the sin^4 (s10 trap #1). Resolved -> baryon carries S_eff; over-production recast as a DENSITY-FEEDBACK knob; UNIFIES with DM screening.
src_paper22/branching_anisotropic.py. (tt)'s b_ch~0.2 was a direction-BLIND channel count -- WRONG: the baryon
(E_6, Q=3) is the HEAVY band, which CARRIES the directional suppression S_eff=sin^4(theta)/[phi^6(1+beta_* rho)]
(claude-hopfion s3, Axiom A2), while the competing light-band channels (elastic 2+2; the Q=1 neutrino byproduct)
do NOT. Isotropically averaging sin^4(=8/15) washed out exactly the phi^6 (+density) suppression. phi^6=17.94 is a
PHYSICAL/topological power (WZW/Bogomolny lambda=phi^6, s6 discriminator -- NOT a frame artifact) so the
suppression is real.
RESOLVED b_ch^res = b_ch^iso * S_eff^p, p = # heavy-band projections in final state 3+1 (ONE heavy baryon => p=1
natural; p=2 bracketed). RESULT (eta_iso=1.76e-4, over by 2.9e5):
  phi^6 alone (beta*rho=1): over-production 2.9e5 -> 1.6e4 (p=1) / 250 (p=2).
  eta_obs reached at (1+beta*rho) ~ 8600 (p=1) or ~16 (p=2).
So resolving the anisotropy CONVERTS the ~1e4-1e5 over-production into a quantity set by ONE physical knob:
beta_*rho at the QCD epoch. UNIFICATION (the payoff): the SAME S_eff denominator (1+beta_*rho) that turns the
director/DM sector OFF in dense regions (s4 screening) ALSO throttles heavy-band baryon formation at the dense QCD
epoch -- one Axiom-A2 mechanism, right direction (dense epoch -> suppressed baryogenesis). NOT a coincidence: S_eff
is one object in two sectors.
CAVEATS [O]: (1) beta_*rho(T_QCD) is cutoff-set + PIVOTAL: modest (~16, p=2) lands eta_obs, but the naive rho~T^4
value (notes/memory flag ~1e47) drives S_eff->0 => UNDER-produces AND would push the transition to CROSSOVER
(killing (qq) first-order). So beta_*rho(T_QCD) now controls BOTH the first-order character AND baryogenesis
efficiency -- the open quantity CONSOLIDATES to one number, with a real large-beta_*rho tension. (2) power p (1 vs 2)
not pinned. NET: eta still unpinned, but "generic transport fudge ~1e-7" -> "ONE physical density-feedback knob
beta_*rho(T_QCD), shared with the first-order diagnosis and DM screening." Do NOT re-isotropise the sin^4 (s10).
[[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-06 (vv) -- USER: rho modifies beta heavily at QCD. The SAME beta_eff=beta_*rho sets first-order AND baryon-branching; overlay -> a self-consistent window EXISTS (p=2), mild tension (p=1).
src_paper22/beta_rho_window.py. User pointed back to (gg)/(hh): the feedback strength is beta_eff=beta_*rho (NOT
beta_* alone), and rho (condensate ~T^4, Lcond~T) was ~1e47x larger early -> beta heavily rho-modified at QCD. KEY:
that SAME beta_eff is the (1+beta_eff) in S_eff=sin^4/[phi^6(1+beta_eff)] of the baryon branching (uu). So one
variable does BOTH jobs, and the heavy rho-enhancement is NOT a problem to dodge -- it is the SOURCE of the
baryon-branching suppression that cures the (tt) over-production. Overlaid barrier(beta_eff) [firstorder_transition.py]
and eta(beta_eff)=eta_iso*S_eff^p on the SAME axis:
  toy first-order WINDOW (b=gamma=1 units): beta_eff in ~[2, 1000] (barrier peaks ~beta_eff=10; >~1e4 over-saturates
    -> crossover, (hh)).
  eta=eta_obs at: beta_eff~15 (p=2), beta_eff~8574 (p=1).
CONSISTENCY is p-DEPENDENT (honest, and I had to FIX the script's text which wrongly claimed p=1 was inside --
its own output said False; CLAUDE.md s1/s9):
  - p=2: eta_obs point beta_eff~15 is INSIDE the first-order window; the window BRACKETS eta_obs (over at low
    beta_eff, under at high) -> a self-consistent beta_eff(T_QCD) EXISTS giving BOTH first-order AND ~observed eta.
  - p=1 (naive: 1 heavy baryon in final state): eta_obs point beta_eff~8574 is OUTSIDE/at the over-saturated edge;
    within the window p=1 OVER-produces by ~10-3000x. Mild TENSION.
So a self-consistent solution EXISTS but is NOT forced; p=1 vs p=2 (how many heavy-band projection factors the
light+light->heavy+light vertex carries) is the pivot and is unpinned. NET for eta: from "22-order transport
mystery" (dd) -> "ONE co-constrained variable beta_eff(T_QCD)", with a self-consistent window (p=2) and a mild p=1
over-production. RESIDUAL OPEN [O]: (1) window LOCATION in physical units = gamma/b from the real feedback energy
E_fb (STILL not extracted -- the needs-development boundary of (ii)/(hh), unchanged); (2) the power p; (3)
beta_eff(T_QCD) independently (co-constrained here, not computed). Do NOT re-walk eta as Sakharov/transport -- it is
now this one density-feedback window question. [[dark-hopfion-cdm-lead]] [[o2-e6-cft-tower]].

## UPDATE 2026-10-06 (ww) -- co-genesis ratio r = n_dark/n_b COMPUTED from group-theoretic sector weights (was assumed O(1)). r=0.5-0.9 -> m_dark~4-11 GeV. Still a consistency.
src_paper22/cogenesis_branching.py. Pinned the status of r (user pushed; checked ALL code+notes): the abundance's
"r~O(1)" had NEVER been positively computed -- only (q) r_forced=0 (the dihedral is NOT a forced per-baryon
byproduct; P18 pion saturates the 2+2->3+1 slot -- the one prior DERIVATION, negative) and (u)'s "comparable
subgroups 12:24" HEURISTIC (adm_abundance.py: "ASSUMED O(1), key uncomputed piece"). The quantum-dim machinery
(cc_p0_linking: <W_j>=S_0j/S_00=d_j) existed but was never applied. The computed BARYON-channel branching b_ch
(branching_eta/anisotropic) is a DIFFERENT quantity (baryon-making fraction in the eta chain), NOT r.
COMPUTED r here, from the group-theoretic SECTOR WEIGHTS of the two McKay subgroups (baryon 2T order 24 McKay E_6;
dark 2D_3 order 12 McKay D_5; verified |G|=sum d^2): PRINCIPLE (stated, ASSUMED) = in the shared-2I co-genesis the
asymmetry branches in proportion to each sector's vacuum weight.
  (a) group ORDER   r=|2D_3|/|2T|=12/24 = 0.50  (regular-rep dim = total DOF)
  (b) #IRREPS       r=6/7              = 0.86  (species count = McKay affine nodes)
  (c) level-1 prim. r=4/3              = 1.33  [CFT-assignment dependent; cross-check]
=> r in [0.5,1.3], ALL O(1); pure group-theory (a,b) give r=0.5-0.9. Omega_dm/Omega_b=5.36 (f=1) -> m_dark =
5.36 m_p/r = 3.8-10.1 GeV. CORRECTION to the "canonical r=1 -> 5 GeV": the group-ORDER weight gives r=0.5 ->
m_dark~10 GeV; 5 GeV was the r=1 idealization. ADM few-GeV scale ROBUST (factor ~2-3), not uniquely 5.
STATUS unchanged in kind: r is now a COMPUTED sector-weight O(1) (not assumed), but (i) branching~static weight is
an ASSUMPTION (true r needs freeze-out reconnection DYNAMICS, uncomputed); (ii) m_dark UNPINNED (D-type no Q_group
-> no mass tower) so the observed ratio FIXES m_dark, not predicts the abundance. Net: Omega_dm/Omega_b stays a
CONSISTENCY; the upgrade is r assumed->computed, and the clean record of r_forced=0 (not a forced byproduct) vs the
O(1) (a shared-vacuum branching weight). Paper 22 sec:abund + eq:adm should carry the computed r=0.5-0.9 + m_dark
~4-11 GeV + the f factor (to match Paper 21's eq:adm). [[dark-hopfion-cdm-lead]].
IMPLICATIONS OF HIGHER m_dark (~10 GeV, r=0.5 vs ~5 GeV, r~1) -- order-of-magnitude check (inline python, knot size
~1/Lcond(QCD), Lcond=76-150 MeV, geometric sigma=pi R^2):
  RETRACTED (my first answer, wrong; CLAUDE.md s9): "higher m_dark stresses the ADM symmetric burn-off". The burn-off
    needs <sigma v> >> 3e-26 cm^3/s; even UNITARITY-capped, <sigma v> ~1.5e-18 (10 GeV) / 5.9e-18 (5 GeV) cm^3/s, i.e.
    ~1e8x the requirement at BOTH masses (geometric ~5e-16..2e-15, above unitarity -- composite, like N-Nbar). So
    annihilation is NOT a discriminator and the ADM mechanism has NO internal preference within 4-11 GeV.
  ACTUAL (mild): (1) same mass density, ~half the number (energy budget into the dark sector fixed by Omega_dm);
    (2) COLDER: non-thermal formation at T_QCD -> v_form ~ T_QCD/m ~0.015c (10 GeV) vs 0.03c (5 GeV); cold either way
    -> partly resolves the earlier "coldness undeterminable" flag: once the ADM match forces m_dark ~ few GeV,
    coldness follows for any m_dark in 4-11 GeV; (3) elastic sigma/m ~3e-3..1.2e-2 cm^2/g (10 GeV) vs ~6e-3..2.4e-2
    (5 GeV) -- BOTH below the 0.1-1 cm^2/g SIDM band, so observable SIDM needs the director-mediated (long-range
    n_hat) enhancement regardless of mass; (4) each knot carries ~10 m_p rest energy -> plausibly more wound ->
    existence (Paper 21 open item) somewhat harder [qualitative; off-ladder dicyclic, E~Q^{3/4} not assumed];
    (5) unchanged: sterility / direct-detection evasion, ADM "coincidence" (r<1 is WHY heavier), no mass-tower
    selection (no Q_group).
  NET: higher m_dark is MILD and creates NO framework tension. Paper statement: m_dark ~4-11 GeV, spread from the
  sector-weight measure, no internal preference. Estimates order-of-magnitude (sigma, unitarity cap, v_form).
WHY f IS ABSENT FROM PAPER 22 (user asked; traced). f = E_res/m_dark was defined in (p)/adm_abundance.py under the
  TRANSIENT-WAKE picture B (m): fraction of an UNWINDING knot's energy trapped as frozen director strain. (w) collapsed
  A/B to the STABLE KNOT A (Q_H conservation forbids unwinding; B also loses the ADM dark number) -> gravitating mass
  = m_dark, f = 1 IDENTICALLY. Paper 22 was written post-(w) with r*m_dark/m_p (since first commit 47dd0fa) --
  CORRECT. Paper 21's eq:adm had imported the wake-era r*f form from adm_abundance.py -- the stale one. FIXED: Paper 21
  eq:adm -> r*m_dark/m_p with a stable-knot sentence; adm_abundance.py docstring marked wake-era / superseded.
  RETRACT my earlier "add f back to Paper 22" (wrong direction).
  Independent supporting argument (user, [I] -- plausible, NOT established): in the wake picture f depends on director
  trapping, but the director coupling S_eff ~ 1/(1+beta_* rho) is SCREENED off at QCD density (s4) -- the same large
  beta_eff used in (uu)/(vv) -- so even an unwinding knot's trapping (f) would be throttled at the Paper 22 epoch.
  Not decisive alone (a wake could freeze later at lower rho); topological reason (w) is the established one.
  Structural Paper 21 vs dynamical Paper 22 split correlates with this: the structural paper inherited the old
  formula; the dynamical one was built on (w). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-07 (xx) -- WAKE ONTOLOGY STANDS; only picture B retired. DM is TWO-COMPONENT. New open question: late-time DOUBLE-COUNTING.
User: "is the wake picture out? I thought geometry-first was fundamental." CLARIFIED (checked P7 text):
  STANDS (foundational, P7:rem:cdm_ontology, P7:rem:polaron): DM = coherent DEFORMATION of the director that a knot
    induces, persisting independently of its source (c_dir->0), "stress without heat"; each knot is polaron-DRESSED.
    Geometry-first: one connected condensate reconfigures; light on the easy axis; knots move slower, filling and
    reshaping the geometry (claude-hopfion s1,s4,s7). Baryon-sourced galactic strain / RAR / cluster fossil network:
    untouched.
  RETIRED by (w) = only picture B (m): the DARK knot unwinds and its frozen residue is the recombination CDM (Q_H
    conservation forbids that unwinding; B also loses the ADM dark number). f was "fraction of an unwinding knot's
    energy left behind" -> meaningless when nothing unwinds.
  Picture A is STILL geometry-first: the stable dark knot is itself a knotted configuration of the same condensate
    (topology, not an added particle -- same status as the proton, D-type slot of the 2I ADE structure). It also
    DRESSES the director like every knot: gravitating mass = m_dark + dressing (factor >=1, NOT the old f). Dressing
    screened at formation density (S_eff->0), may switch on late as rho falls -> dark knots could also source late
    strain. [O] uncomputed.
  "Back to a particle?" PHENOMENOLOGICALLY yes (localized, stable, number-conserved, a^-3, ADM), ONTOLOGICALLY no.
  STERILITY (Paper 21 cor:sterile): closure fails -> no quantised SM charge; no SU(3)_L rep for Q (no EM); no weak
    doublet partner; not in 2T (no colour) -> full SM singlet. Couples to: gravity, its own Z_2, and the shared
    director n_hat (SIDM + dressing). WORDING FIX OWED: cor:sterile says "interacting ONLY gravitationally and
    through its own Z_2" but Paper 21 sec:abundance says "self-interacting through the shared director field" --
    compatible (director coupling is geometric-sector, not SM) but the corollary should name it.
  => DM is TWO-COMPONENT: (1) EARLY stable dark knots fill the recombination Omega_DM gap [conditional: existence
    open, abundance = consistency, m_dark~4-11 GeV]; (2) LATE baryon-sourced director strain = galactic dynamics/RAR
    [established, P7]; (+ possible dark-knot dressing).
  NEW OPEN QUESTION (not previously in notes, checked): late-time DOUBLE-COUNTING. Cold knots with Omega~0.26 fall
    into CDM-like halos around galaxies; P7 already reproduces the RAR from baryon-sourced strain ALONE -> knot halo
    + strain could exceed the RAR. Candidate resolutions: (1) density screening (S_eff ~ 1/(1+beta_* rho)) suppresses
    strain where knot halos are dense; (2) knot halos stay diffuse/cored at galactic scales; (3) the empirical a_0 /
    strain normalisation already absorbs the knot halo. Assessed next. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-07 (yy) -- DOUBLE-COUNTING ASSESSED against P7's actual text: (3) is the only structural escape; (1) only in dense cores; (2) not available. Real tension.
P7 FACTS USED (read, not recalled): (a) source: "visible matter ... imposes orientation boundary conditions on the
director" (P7:sec:director_strain) -- the mechanism is KNOTS imposing BCs, nothing baryon-specific; (b) magnitude:
condensate-derived "to within about an order of magnitude", a_0 coefficient needs a cutoff/crossover ratio ~phi^6
and "is not derived here" (P7:rem:elastic_solid, P7:sec:accel_scale); (c) screening: S_eff = sin^4/[phi^6(1+beta_*
rho)]; (beta_* rho_sun)^-1 ~ 1e-6 at Solar-System densities, unsuppressed in voids (P7:rem:lorentz); planets->dwarfs
spans ~30 dex (s4); (d) at large cutoff the strain is "phenomenologically indistinguishable from collisionless CDM"
(P7:rem:elastic_solid).
(1) SCREENING -- UNLIKELY as the resolution. Knot halos ARE galactic-density matter (cosmic mean ~5x baryons, the
  dwarf regime where P7 has MAX strain). A factor ~5 in rho on a ~30-dex screening lever arm barely moves beta_* rho.
  Can only matter in genuinely dense environments (cluster cores, already screened in P7).
(2) DIFFUSE/CORED KNOT HALOS -- NOT AVAILABLE from established structure. Knots are cold (v_form~0.015-0.03c, ~0 today),
  number-conserved, a^-3, with bare sigma/m ~1e-3..1e-2 cm^2/g (ww) -- below the 0.1-1 coring band -> they cluster as
  collisionless CDM (cuspy, galaxy-scale halos). Only a director-mediated interaction could change this, and P7 makes
  dynamical elastic response conditional on the (unfixed) cutoff. Would need new physics.
(3) NORMALISATION ABSORBS -- the ONLY escape the structure leaves room for: the strain/a_0 coefficient is
  explicitly underived (O(1)-O(10) slack), so the observed excess can be SHARED (knot halo + reduced strain). COST:
  the RAR's TIGHTNESS is what the strain explains and CDM halos don't; a cosmic-ratio knot halo does not track
  baryons point-by-point -> adds scatter. Viable only if the knot halo is sub-dominant inside rotation-curve radii in
  the deep-a_0 regime (dwarfs/LSBs) -- i.e. CDM's cusp/diversity problems partly reappear and must be checked.
STRUCTURE MAKES IT WORSE, NOT BETTER (fact a): dark knots are knots, so they should ALSO source strain (dressing,
  (xx)) once unscreened -> total = baryons + knot halo + baryon-strain + knot-strain. "Strain tracks baryons" becomes
  "strain tracks all knots"; the RAR (vs g_bar only) is then protected only if knot-sourced strain is small or
  follows baryons. Open.
UPSIDE (same coin): a cosmic-ratio knot component naturally supplies CLUSTER mass (ΛCDM works for clusters), with
  strain screened in cores -> addresses the s4 cluster deficit (strain ~5x vs lensing ~28x) WITHOUT needing the
  conditional fossil disclination network. So the knots help clusters exactly where they hurt galaxies.
VERDICT: most likely = (3), with (1) acting only in dense cluster cores; (2) unsupported. It is a REAL tension, not
  resolved: the framework now has the classic MOND(galaxies)/CDM(clusters) complementarity built in, and the RAR
  tightness test is decisive. CONCRETE NEXT TEST: toy rotation-curve fit -- baryons + cosmic-ratio NFW knot halo +
  strain with FREE coefficient (and knot-sourced strain on/off) -> can a reduced strain coefficient reproduce the RAR's
  small scatter across dwarf -> massive galaxies? [O] [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-07 (zz) -- USER: "there is anisotropy missing; MOND has none." (yy) ISOTROPISED (s10 trap #1, 2nd time). Prior axis-resolved work bears directly.
(yy) treated knot halo, strain AND screening as isotropic = MOND+CDM in disguise. notes/single_galaxy_winding_analysis.md
(2026-09-01) already has the axis-resolved machinery; re-reading it changes the verdict:
  PRIOR RESULT 1 -- collective hedgehog CANCELLATION (topological, robust: isotropic, oriented, shadowed, splay/bend
    reweighted all give eta<1): two radial strain sources point OPPOSITELY between them -> gradients cancel; closer =
    more cancellation. => knot-sourced strain is NOT additive. A swarm of dark knots sourcing radial hedgehogs cancels
    between knots, leaving at most the aggregate's outer field. This DEFUSES (yy)'s "structure makes it worse"
    (baryons + knot halo + baryon-strain + knot-strain): the knot-strain term largely cancels.
  PRIOR RESULT 2 -- the M~r halo is held by a GLOBAL self-consistent condition: it extends until g = G M_enc/r^2 ~ a_0,
    M_enc = TOTAL enclosed mass. A knot halo raises M_enc -> pushes the strain's M~r onset OUTWARD -> strain and knot
    halo are NON-ADDITIVE (strain only operates where TOTAL g < a_0). A 4th option missed in (yy).
  WHERE THE ANISOTROPY ENTERS (MOND has none):
    (A) the preferred direction e_r = grad rho/|grad rho| (P7:rem:lorentz) is set by the TOTAL density -> the knot halo
        does not just add mass, it RE-ORIENTS the axis the strain couples along; S_eff = sin^4(theta)/[phi^6(1+beta_*
        rho)] with theta = director tilt from e_r. (yy)'s option (1) kept only the 1/(1+beta_* rho) magnitude and dropped
        the sin^4(theta) direction.
    (B) the inertial (polaron) mass is a TENSOR M_ij = chi INT d_i n . d_j n (P7:eq:meff): for a radial texture d_r n = 0,
        only TANGENTIAL components survive -> the dressing acts on tangential motion (rotation curves) and not on radial/
        vertical motion; lensing sees the gravitating Frank energy. The knot halo is direction-blind.
    (C) two rotation DOF: G_phi (about e_r, soft, density-blind) vs G_theta (polar tilt, changes sin^4 theta, chi_theta
        density-dependent, edge-localised) -- the halo-edge anisotropy (prior note).
  CONSEQUENCE: "double counting" is DIRECTION-DEPENDENT, not one normalisation. The kinematic RAR is a TANGENTIAL relation
    (rotation ⊥ radial e_r, sin^4 theta = 1, M_ij tangential -> strain maximal there); vertical-force (K_z), radial/
    pressure-supported dynamics and lensing weight the isotropic knot halo relatively more. => the framework predicts
    a PROBE-DEPENDENT dark excess (kinematic-RAR vs lensing-RAR vs vertical dynamics disagree in a definite pattern) --
    absent in BOTH MOND (isotropic g) and CDM (isotropic halo). Tension -> falsifiable discriminator.
  STATUS: (yy)'s verdict "(3) only escape" is SUPERSEDED: non-additive global cutoff (result 2) + collective knot-strain
    cancellation (result 1) both reduce the overlap, and the residual overlap is axis-dependent (A,B). Not yet
    quantified. UNRESOLVED from the prior note: theta's exact role for an aligned radial hedgehog (n = e_r -> sin^4=0)
    vs tilts; chi_theta not computable from frank_run (needs H's polar-tilt dependence). NEXT: axis-resolved toy --
    disk + knot halo + strain with the global a_0 cutoff on TOTAL M_enc; compare tangential (rotation) vs vertical vs
    lensing excess. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-07 (aaa) -- SEMI-ANALYTIC DOUBLE-COUNT TEST (anisotropy resolved): the global a_0 onset SELF-CORRECTS most of it; anisotropy -> a lensing/rotation PREDICTION that MEASURES the knot share.
src_paper7/dm_polarization/anisotropic_double_count.py (new; extends phase2c_A1 machinery semi-analytically, since the
M~r config must be IMPOSED). Toy: exponential-disk baryons (M_b 1e7..3e11), NFW knot halo M_200=f_knot*M_b (c=10),
strain g_s = zeta*g_src*[1/(1-exp(-sqrt(g_tot/a_0)))-1] (amplitude from SOURCE mass, onset from TOTAL Newtonian g =
prior result 2), knot-sourced strain fraction kappa_k (prior result 1). Anisotropy (radial hedgehog): probe at angle psi
to r_hat sees strain weight B: 1.5 sin^2 psi (M_ij), A: (15/8) sin^4 psi; rotation psi=90, radial psi=0, lensing = scalar 1.
Rotation fit fixes zeta_eff = w_t*zeta -> anisotropy PREDICTS lensing and radial-orbit excess.
GATE: first stand-in (sharp onset at outermost a_0 crossing, linear M growth) FAILED (zeta_eff 0.62, rms 0.18 dex,
  massive -0.14 dex) -- stand-in artifact, discarded. Smooth onset passes (trivially: f_knot=0 = the empirical RAR).
RESULTS (kinematic RAR fit; observed intrinsic scatter ~0.06 dex):
   f_knot kappa | zeta_eff rms   dwarf  massive | lens/rot B   A    ISO | radial/rot (A,B)
   0      0     |  1.000  0.000  0.00   0.00    |   0.67     0.53  1.00 |  0.00
   1      0/0.1 |  1.05-1.09 0.03 -0.01 +0.03   |   0.69     0.57  1.00 |  0.08
   5.36   0/0.1 |  1.01-1.13 0.054-0.065 -0.02 +0.075 | 0.76-0.77 0.67-0.68 1.00 | 0.29-0.31
READ:
 (1) DOUBLE-COUNTING LARGELY SELF-CORRECTS through the non-additive onset: even a COSMIC-ratio knot halo leaves RAR
     scatter ~0.05-0.065 dex (~observed) with zeta_eff ~1 -- the knot halo raises g_tot and pushes the strain out
     where knots matter. Residual: massive galaxies over-predicted ~+0.075 dex (~19%), dwarfs ~-0.02 dex. Mild
     tension, toy-level. (yy)'s "(3) only escape, real tension" is SUPERSEDED: it was the additive/isotropic framing.
 (2) ANISOTROPY PREDICTION (absent in MOND and CDM): at g_bar~1e-12 m/s^2 lensing sees 0.67-0.77 (B) / 0.53-0.68 (A)
     of the rotation-curve dark excess (~ -0.1 to -0.15 dex in g_obs); purely radial orbits see 0-0.3 of it.
 (3) The ratio DEPENDS on f_knot (knot halo is isotropic): lens/rot -> 2/3 (B) or 8/15 (A) if strain dominates,
     -> 1 if knots dominate. => lensing-RAR vs kinematic-RAR (and radial vs tangential orbit dynamics, e.g. dSph/
     elliptical Jeans anisotropy) MEASURES the knot share. A framework-specific, falsifiable discriminator.
CAVEATS [O]: toy (spherical baryons, fixed c=10, cosmic f relative to galaxy baryons, no feedback/baryon loss, no
  input scatter). M_ij->per-probe mapping is the B HYPOTHESIS; lensing != dynamics is effectively a GRAVITATIONAL SLIP
  -> must be checked against P7's WEP/PPN sector. CHECK AGAINST DATA (not done): the weak-lensing RAR of isolated
  galaxies (KiDS; reported broadly consistent with the kinematic RAR) -- a -0.1..-0.15 dex lensing deficit at
  g_bar~1e-12 may already be constrained; needs the actual error bars. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (bbb) -- GRAVITATIONAL-SLIP + KiDS check: the director texture is a GLOBAL MONOPOLE (lenses, no Newtonian force). Retracts (aaa)'s slip; surfaces a real P7 galaxy-sector tension; KiDS constrains it.
src_paper7/dm_polarization/monopole_stress_check.py (sympy-verified). Followed research-log skill (verify vs code;
retract cleanly).
RETRACTION: (aaa)'s lensing/rotation slip (0.53-0.77 via the polaron INERTIAL tensor M_ij) is WITHDRAWN -- wrong
object. The polaron m_eff (P7:rem:polaron) is the SOURCE baryon's own inertia (P_i=M_ij V_j), not the halo that test
stars/photons respond to. Both rotation stars and lensing photons respond to the GRAVITATING Frank energy; the slip,
if any, is set by THAT energy's anisotropic STRESS, not the source inertia.
[E] (verified, sympy, both one-constant Frank AND pure splay, any sign of K): the baryon-sourced director hedgehog
n=r_hat has field-theory stress  rho=+c/r^2, p_r=-c/r^2 (radial TENSION), p_t=0  ->  rho + sum p_i = 0. This is the
Barriola-Vilenkin GLOBAL MONOPOLE: ZERO Komar/active gravitational mass (weak-field R_00=4piG(rho+sum p_i)=0 -> NO
Newtonian force on test particles) but a solid-angle deficit that DEFLECTS light (lenses). So a FREE director texture
LENSES WITHOUT PULLING -- the opposite extreme of normal matter.
[O] TENSION with P7's flagship galaxy result: P7's flat rotation curves use M_enc=INT u dV ~ r with u=K/r^2 as a
NEWTONIAN gravitating density -- but the proper GR active mass (rho+sum p) VANISHES for this texture, so the strain
exerts no Newtonian pull. P7's "pressureless" (c_dir->0) is absence of JEANS/kinetic pressure SUPPORT; the STATIC
elastic stress p_r=-u is a different object and is present (P7:rem:elastic_solid even calls it a shear-rigid elastic
SOLID). P7 never computes the strain's galactic-scale lensing or active mass (only Solar-System gamma_PPN via
Vainshtein, and the bullet remark which ASSUMES the strain is the lensing mass). So this is a genuine gap, not
previously in the notes.
LIKELY RESOLUTION [I/O]: the halo is NOT a free texture -- it is pinned by baryon BCs and held by the GLOBAL a0
self-consistency condition (single_galaxy_winding_analysis). A BOUND system's active mass includes the binding/pinning
stress (the baryon force that holds the winding), which can restore rho+sum p != 0. Unmodeled; needs the pinned-config
stress, not the free-texture one. The slip is then set by the bound-restored fraction (lenses+pulls equally) vs
residual free-texture (lenses only).
[E, data] KiDS-1000 weak-lensing RAR (Brouwer et al. 2021, A&A 650 A113): g_obs from LENSING follows the SAME RAR as
rotation curves, extended ~2 decades lower in acceleration; agrees with MOND/EG predictions. => real galaxy halos LENS
and PULL equally (small slip). This EMPIRICALLY excludes the pure-texture (global-monopole) limit as the whole DM:
a pure texture would show lensing mass with ~no rotation-curve mass. The data REQUIRE the active mass to be mostly
restored -- i.e. bound/pinned strain, and/or the KNOTS (normal a^-3 matter: rho>0, ordinary pressure, positive active
mass -> pull AND lens normally) carrying the dynamical+lensing DM. Favors the two-component picture (xx) with the
dynamical mass NOT a free texture.
NET vs the double-count thread (yy)/(zz)/(aaa): the slip is NOT the (aaa) polaron number. The physical question is
whether the strain's active mass is free-texture (lenses only) or bound-restored (normal). KiDS says ~normal (no large
slip), so the framework must deliver a bound-restored strain OR lean on the knots for the dynamical RAR. CONCRETE OPEN
TASK: compute the PINNED director-halo stress-energy (baryon BC + a0 binding) and its rho+sum p -- does binding restore
the Newtonian pull P7 assumes? That is the real unresolved calc, and it underlies P7's rotation-curve mechanism itself,
not just the slip. CAVEATS: P7 DOES derive Einstein-frame GR (so weak-field R_00 applies to the emergent metric); the
escape is the binding stress or a (screened, per P7) fifth force, not a non-GR force law. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (ccc) -- PINNED-HALO ACTIVE-MASS calc: in standard GR the strain LENSES flat (M~r) but PULLS Keplerian (concentrated). Flags an underived director->metric coupling in P7's flagship galaxy result.
src_paper7/dm_polarization/active_mass_profiles.py (sympy-verified). SOURCE-CHECK first (research-log skill): the
prior "monopole" discards in paper7_audit_changes.md were (a) the SCALAR-magnitude rho hedgehog -> Keplerian (replaced
by the director) and (b) direction-only sin^4θ inert for spherical sources -- NEITHER is the director strain's GR
ACTIVE MASS. So (bbb)/(ccc) are new, not re-walking. (Audit DID flag the flat-curve config as an ASSUMED disclination,
imposed-not-relaxed -- adjacent.)
CALC: slow-matter dynamics couples to the Komar active-mass density rho+sum p_i (weak-field R_00=4piG(rho+p_r+2p_t));
lensing couples to energy density + solid-angle deficit. General (Euler on a (grad n)^{2k} energy): rho+sum p_i =
(2k-2) f. VERIFIED (sympy, hedgehog n=r_hat):
  Frank (k=1): rho=1/r^2, rho+sum p = 0            -> global monopole: M_enc[energy]=INT rho 4pi r^2 dr ~ R
     (flat-curve-amplitude LENSING) but ZERO Newtonian active mass -> NO pull.
  Faddeev J4 (k=2, the Skyrme STABILISER, not optional): f4=1/r^4, rho+sum p = 2 f4 > 0 -> the ONLY dynamical pull,
     but M_active=INT 2 f4 4pi r^2 dr -> CONST (concentrated) -> Keplerian dynamics, NOT flat.
RESULT: in standard GR the director strain gives flat-curve LENSING (M~r) but CONCENTRATED/Keplerian DYNAMICS ->
lensing >> dynamics. This CONTRADICTS (a) P7's flat-ROTATION-curve (dynamical) claim and (b) KiDS, where the
weak-lensing RAR = kinematic RAR (Brouwer 2021). The Frank result rho+sum p=0 is CONFIG-INDEPENDENT (holds for any
f=1/2 W(n)|grad n|^2, incl. the sin^4θ-weighted P7 energy -- W is an overall factor, cancels in the trace), so
relaxation/pinning of the Frank part cannot restore the pull; only a NON-quadratic term (J4, or a direct
director-matter force) carries active mass.
RESOLUTION [O] (the real gap): this applies the standard-GR MATTER-SOURCE formula to the director. P7 asserts the
Frank energy "gravitates (E=mc^2)" as energy density rho (pressureless), but (i) the strain is NOT pressureless as a
STRESS source (p_r=-rho; P7's "pressureless" is Jeans/kinetic, (bbb)), and (ii) the director is part of P7's
EMERGENT-gravity sector and P7 NEVER DERIVES how the director strain enters the emergent Einstein equations (as rho,
or rho+sum p). If it sources Phi via rho (not rho+sum p) -- e.g. a direct director-matter coupling / fifth force --
dynamics is flat; but P7 SCREENS the fifth force at galactic radii (G_eff/G_N=1+6e-5, P7:sec:halo_screening) and says
flat curves are NOT a G_eff variation. So P7 is caught between: (A) strain gravitates as rho+sum p (standard GR) ->
no flat dynamics; (B) strain pulls matter directly (fifth force) -> flat dynamics but contradicts P7's own screening.
This is an INTERNAL inconsistency / underived coupling in P7's galaxy sector, now precisely located.
KiDS CONSTRAINT: lensing RAR = kinematic RAR means real halos lens AND pull equally -> whatever carries the galactic
dynamical mass must ALSO lens equally (normal-matter-like). The pure strain (lens>>pull) fails; NORMAL matter does it.
=> favors the KNOTS (a^-3, rho>0, ordinary pressure -> pull+lens equally) as the galactic dynamical DM, with the strain
at most a subdominant lensing-only piece -- a reframe of P7's director-strain-as-DM thesis. [[dark-hopfion-cdm-lead]].
NOT fatal-claimed: hinges on the underived director->metric coupling; P7 derives GR for the SCALAR sector but not the
director sector's gravitational coupling, which is the thing to derive next. Config caveat: bare hedgehog; true pinned
relaxation shifts the J4 profile (not the Frank rho+sum p=0, which is config-independent).

## UPDATE 2026-10-11 (ddd) -- USER REFRAME ("halo follows the lens, but drags?") = P7's OWN polaron picture. CORRECTS (ccc) over-claim + (bbb) over-retraction: P7's dynamical mass is the COMOVING momentum mass, not the static Komar mass.
User asked: "so the halo follows the lens, but drags?" -- this is P7:rem:polaron almost verbatim ("a moving baryon
DRAGS its distortion cloud, so the halo produces gravitational dynamical friction"). Forced a re-examination (skill
s5: a pushback is a cue to redo, not defend). I had computed the WRONG gravitating quantity.
P7 HAS TWO inequivalent "gravitating" quantities it never cleanly separates:
  (i) STATIC Frank energy density u (sec:director_strain, flat_curves): M_enc=INT u dV ~ r. My (bbb)/(ccc) computed
      its Komar mass rho+sum p = 0 (global monopole) -> lensing only.
  (ii) POLARON momentum mass (rem:polaron): the dressing COMOVING with the baryon, n_hat(x-Vt), momentum P_i=M_ij V_j,
      M_ij = chi INT d_i n_hat d_j n_hat d^3x,  m_eff=(chi/3)INT|grad n|^2 ~ r.  "The quantity that gravitates and acts
      dynamically." This is what P7's flat ROTATION curves actually rest on -- a COMOVING/dynamical quantity, NOT the
      static texture.
CORRECTION to (ccc): "in GR the strain lenses flat but pulls Keplerian -> P7's flat curves fail" was TOO STRONG. It
refuted only reading (i) (static energy as a GR matter source). P7's mechanism is (ii) the polaron momentum mass. The
two DIVERGE (static Komar = 0 vs polaron ~ r); (ccc) showed (i) fails, which just means P7 MUST lean on (ii).
CORRECTION to (bbb): the M_ij tangential tensor I used in (aaa) and then dismissed as "source-baryon inertia, not the
halo" IS P7's gravitating halo (rem:polaron identifies M_ij as what gravitates). Partial un-retract. BUT the per-probe
slip numbers (aaa, 0.53-0.77) still don't stand as given: shell-averaged M_ij ~ INT dOmega (dij - r_i r_j) = (8pi/3)dij
= ISOTROPIC, so the SPHERICAL halo's mass is isotropic (no net lensing-vs-rotation slip from the tensor); the
anisotropy survives only as a radial-vs-tangential ORBIT distinction at a point (dSph/elliptical Jeans anisotropy),
not an angular halo asymmetry.
KEY (favorable) CONSEQUENCE: m_eff ~ INT|grad n|^2 ~ INT u ~ the LENSING mass -> the polaron predicts dynamics ~
lensing ~ r (both from the same gradient integral) -> FLAT curves AND lensing=kinematic RAR -> CONSISTENT with KiDS
(Brouwer 2021). So P7's own mechanism naturally gives the KiDS consistency my (ccc) said it failed.
THE REAL OPEN QUESTION (sharpened, [O]): is the polaron COMOVING MOMENTUM mass a legitimate CONSERVATIVE gravitating
halo (sources a central potential -> stable flat curves), or does its "drag" character make it DISSIPATIVE (decays
orbits, cannot sustain a curve)? Equivalently: does the inertial M_ij (from field momentum) equal the active
gravitational mass? The static texture says active mass = 0; the comoving dressing carries momentum/kinetic
stress-energy (chi (d_t n)^2) ABSENT from the static texture, which is where a nonzero active mass could live. This is
the underived director->metric coupling, now precisely: compute the gravitating mass of the COMOVING (not static)
distortion -> conservative central (flat curves, KiDS-ok) or dissipative drag? CONFIG caveat unchanged. Credit: user's
reframe located this. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (eee) -- USER picture resolves (ccc)/(ddd): the DYNAMICAL (gravitating) halo mass is the ROTATING director's KINETIC term, conservative not dissipative. Static lenses, moving pulls -- one mechanism, = P7's polaron.
User: "the geometry rotates locally; the lens is the lightest thing that follows the geometry; the halo drags after
-- maybe not dissipation; and not sure if baryon-tracking or dark-hopfions or a combination." Pursued it; it holds and
resolves the (ccc) worry.
[E] (sympy, rotating_active_mass.py; Hilbert T_munu): the Komar active-mass density of the director is
   rho + sum p_i = 2 chi (d_t n)^2.
STATIC texture (d_t n=0): 0  -> global monopole, LENSES (deficit) but NO Newtonian pull (= (bbb)/(ccc); the gradient
part CANCELS). ROTATING/moving director (d_t n != 0): +2 chi (d_t n)^2 > 0 -> PULLS. So ALL the dynamical pull comes
from the TIME-DEPENDENCE -- the orbiting/dispersing baryons (and/or dark hopfions) DRAGGING the director. Exactly the
user's "lens follows the static geometry, halo drags after."
NOT DISSIPATIVE (user's "maybe not dissipation", CONFIRMED established): c_dir->0 => medium FROZEN, no propagating
modes, no radiation, F_drag->0; "a frozen medium exerts (elastic) RESTORING, not drag" (cluster_wake_calc.py). So the
drag is REACTIVE/elastic (conservative) -> it SUSTAINS orbits (flat curves), does not decay them. This answers (ddd)'s
open conservative-vs-dissipative question: CONSERVATIVE.
= P7's POLARON, unified: for a comoving distortion n(x-Vt), d_t n = -V.grad n, so INT 2chi(d_t n)^2 = 2 M_ij V_iV_j,
M_ij=chi INT d_i n d_j n (P7:rem:polaron). The polaron momentum mass IS the rotating-director active mass. So P7's two
"gravitating" quantities reconcile: static energy = LENSING (global monopole), polaron/kinetic = DYNAMICAL pull.
RESOLVES (ccc) over-claim: the strain does NOT fail flat curves; the flat-curve pull is the kinetic/polaron term
(nonzero), the static Komar mass (0) is the LENSING-deficit piece. And both scale as the same INT|grad n|^2 ->
lensing ~ dynamics ~ r -> CONSISTENT with KiDS (lensing RAR = kinematic RAR, Brouwer 2021).
PROFILE (self-consistency sketch [I], NOT derived): active density ~ 2 chi <v^2> |grad n|^2, <v^2>=baryon kinetic
(rotation v_c^2 OR dispersion sigma^2); |grad n|^2~1/r^2 -> M_active~r (flat) = the lensing profile. chi=K/c_dir^2->inf
(cold) so slow motion suffices.
OPEN [O]: (1) does the director co-move at the baryon kinetic scale, and how much (sets the a0 normalization)? (2) the
r-profile and a0 are self-consistent, not derived. (3) ROTATION- vs DISPERSION-supported: <v^2> must be the LOCAL
kinetic scale (so dSph, DM-rich + low rotation, work via sigma^2), NOT bulk rotation -- else dSphs mispredicted;
needs checking. (4) SOURCE stays OPEN (user): baryon-tracking and/or dark-hopfion co-moving dressing -- both drag the
director; which dominates undetermined. NEXT: the rotating/dispersing-halo PROFILE calc (M_active(r) and lensing(r)
from 2chi<v^2>|grad n|^2 with a realistic director + v(r)) -> does it give flat curves + lensing=dynamics + the a0
scale, for disks AND dispersion systems? [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (fff) -- SEGUE: (eee) derivation written into Paper 7 (rem:active_mass); Koide-exact m_tau prediction located (in papers as K=2/3; sharp falsifiability NOT). Next: profile calc.
PAPER 7: added Remark P7:rem:active_mass after rem:polaron (+ BarriolaVilenkin1989 bibitem): the static director
texture has rho+sum p=0 (global monopole: deficit LENSING, no pull); the dynamical pull is 2chi(d_t n)^2 = the comoving
polaron mass; reactive/conservative (c_dir->0); both lensing & dynamics ~ INT(grad n)^2 -> common RAR; profile/a0 =
cutoff-set, not fixed (matches P7's existing hedge). Verified: cites/labels/envs balanced. States current content
plainly (s7), no change-narration.
KOIDE (user side-check, verified /tmp): K=2/3 EXACT is in Paper I (P1:thm:koide, from Z3; m_tau/m_e 0.03%). The SHARP
falsifiability framing (exact Koide -> m_tau=1776.969 MeV; PDG 1776.86+/-0.12 is 0.91 sigma low; measure to +/-0.05
MeV to test Q=2/3 exact) is NOT in the papers -- a clean lepton-sector falsifiable prediction worth adding to Paper
I/II later (NOT DM; parked). Numbers reproduced exactly (m_tau=1776.969, 0.91 sigma, Q=0.6666667).
NEXT (back to DM): the rotating/dispersing-HALO PROFILE calc -- M_active(r)=INT 2chi(d_t n)^2 and the lensing mass
M_lens(r) from a director |grad n|^2~2/r^2 with a matter kinetic profile <v^2>(r) (rotation v_c^2 OR dispersion
sigma^2) -> do both go ~r (flat curves + extended lensing halo), is the dynamics/lensing ratio O(1) (KiDS), and does
a0 emerge? Key subtlety to handle: the LENSING source for the kinetic/stiff anisotropic stress (both lensing &
dynamics are kinetic-dominated for cold director, so ratio ~O(1), not the naive c_dir->0 blowup). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (ggg) -- ROTATING-HALO PROFILE: flat curves + dSph halos EMERGE self-consistently; universal a0 does NOT (scale-free toy -> a0 stays cutoff-set, = P7's standing). KiDS-shape OK.
src_paper7/dm_polarization/profile_rotating_halo.py. Self-consistent closure for the comoving-director active mass
(eee): rho_act=4 chi eta <v^2>/r^2, dM_act/dr=(A/r)(M_b+M_act), A=16 pi chi eta G (dimensionless self-consistency
amplitude). Results:
(1) FLAT CURVES [E,toy]: A~O(1) gives an extended ~flat v_c(r) (log-sd 0.003 dex over r>3Rd at A=1, v_flat~240 km/s).
    The rotating-director bootstrap DOES produce flat rotation curves. A<~0.5 under-produces, A~2 runs away.
(2) DISPERSION systems [E,toy]: sigma=const -> M_act=16 pi chi eta sigma^2 r ~ r directly. dSphs (DM-rich, ~no
    rotation) get an M~r halo from sigma^2 alone -> the mechanism needs LOCAL <v^2>, not bulk rotation. dSph concern
    (eee) RESOLVED: dispersion-supported systems work.
(3) UNIVERSAL a0 -- NO [O]: crossover (M_act=M_b) is at FIXED r_cross/Rd~0.43 for all M_b -> g_cross ∝ GMb/Rd^2 ∝
    Mb^0.4 (Rd∝Mb^0.3), NOT a universal acceleration (g_cross/a0 = 0.15,0.37,0.94,2.36 for Mb=1e8..1e11). EXPECTED:
    the toy's A is DIMENSIONLESS (scale-free: only G,Mb,Rd), so it structurally cannot pick an absolute acceleration.
    A universal a0 requires the DIMENSIONFUL cutoff coupling (chi=K/c_dir^2 carrying cH0); that is P7's a0=cH0/(2pi)
    "cutoff-set, not derived" -- OUTSIDE the toy's scope, neither confirmed nor refuted here.
LENSING: both dynamics (INT 2chi(d_t n)^2) and the kinetic-dominated lensing energy ~ r -> common RAR SHAPE
(KiDS-consistent in shape); the O(1) normalization ratio needs the anisotropic-stress lensing kernel (flagged, not
pinned).
NET (DM thread (bbb)->(ggg)): the rotating/comoving-director mechanism (USER's "geometry rotates, lens follows, halo
drags, not dissipation") CONSISTENTLY gives flat curves + extended lensing + dSph halos, resolving the (ccc)
global-monopole tension and matching KiDS in shape. REMAINING (both pre-existing P7 gaps, not new): (a) the a0
NORMALIZATION + universality = the cutoff coupling (P7's "cutoff-set"); (b) the exact lensing/dynamics O(1) factor;
(c) SOURCE open (user): baryon-tracking and/or dark-hopfion dressing. Written into Paper 7 as rem:active_mass (the
static-lenses/comoving-pulls derivation; profile/a0 left cutoff-set, consistent). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (hhh) -- CUTOFF COUPLING (source-check: already done -> a0~cH0) + LENSING/DYNAMICS slip COMPUTED = 1+(c_dir/v)^2/4 ~ 1 (KiDS), tied to the sector's coldness.
TWO asks (user), both resolved:
(1) CUTOFF COUPLING -- ALREADY DONE (source-check, results.md + crossover_match.py/cutoff_embedding.py): the absolute
scale my (ggg) toy lacked IS supplied by the director mass m_xi~H0. a0 = (band_K(Lam*)) cH0, band_K=5.23e-3 Lam^1.36;
a0 ~ cH0 RIGHT ORDER, NO tuning (robust). Exact a0=cH0/(2pi) needs the cutoff/crossover hierarchy m_xi/E* ~ 12-38.
SUGGESTIVE [I] (framework-native): m_xi/E* = phi^6=17.94 (the S_eff normalization itself) brackets a0 within ~1.5x
(light-axis 1.49x obs, heavy 0.54x); does NOT nail the 2pi. So a0~cH0 robust, exact coeff cutoff-set (phi^6 lead,
unproven). The NEW comoving-kinetic active mass (eee) inherits the SAME |grad n|^2 geometry + m_xi cutoff -> same
a0~cH0 (the extra kinetic v^2/c_dir^2 factor folds into the same cutoff-set coefficient). (ggg)'s missing universal
a0 was the dimensionless-toy artifact; the dimensionful m_xi supplies it. a0 NOT independently derived (=P7 standing).
(2) LENSING/DYNAMICS -- COMPUTED [E] (lensing_dynamics_slip.py, sympy-checked stress decomposition):
   STATIC texture: anisotropic (p_r=-rho, p_t=0) -> dynamics rho+p_r+2p_t=0 (global monopole), LENSING-only (deficit).
   COMOVING kinetic: isotropic stiff (p_r=p_t=rho_kin) -> Phi=Psi -> lensing mass = dynamics mass (NO slip).
   Sum: dynamics=4 rho_kin (kinetic only pulls); lensing=4 rho_kin + rho_static (anisotropic texture adds to lensing
   only). With d_t n~v grad n, chi=K/c_dir^2: rho_static/rho_kin=c_dir^2/v^2, so
       lensing/dynamics = 1 + (1/4)(c_dir/v)^2.
   => controlled by the SAME coldness c_dir/v the director-DM sector already needs (c_dir->0, pressureless
   clustering). Cold director (c_dir<<v_gal) -> ratio->1 = weak-lensing RAR = kinematic RAR (KiDS/Brouwer 2021);
   residual is a DEFINITE-SIGN prediction (lensing slightly ABOVE dynamics) ~ (c_dir/v)^2, a future precision test.
   KiDS consistency thus REQUIRES the sector's coldness -- not a new assumption. CAVEAT: c_dir/v cutoff-set, not
   pinned; the isotropic-vs-anisotropic split is the clean physics.
NET (galaxy thread (bbb)->(hhh) COMPLETE at mechanism level): rotating/comoving director -> flat curves (ggg) +
dSph halos (ggg) + a0~cH0 (crossover_match) + lensing=dynamics to O((c_dir/v)^2) (hhh). Static lenses, moving pulls,
reactive not dissipative. All consistent with the galaxy data (RAR, KiDS) at the mechanism level. REMAINING opens
(all pre-existing, cutoff-set): exact a0 coefficient (phi^6 lead), the numeric c_dir/v, and the SOURCE (baryon vs
dark-hopfion, user-open). Paper 7 carries rem:active_mass. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (iii) -- Koide tau prediction added to Paper I; phi^6 lead checked: it's a VACUUM-director quantity (source-blind), but the TIGHT baryon-RAR decides the galactic source = baryon-tracking (dark hopfions = separate early/cluster role).
PAPER I: added Remark P1:rem:koide_tau after the Koide table: exact K=2/3 (P1:thm:koide) + known m_e,m_mu => m_tau=1776.969
MeV; PDG 1776.86+/-0.12 is 0.91 sigma low; residual ~6e-6 in K too small for any known correction; falsifiable by
m_tau to +/-0.05 MeV. Verified (refs P1:thm:koide, PDG2024 resolve; label; env balanced). Lepton-sector; done.
PHI^6 LEAD (user: "why would it fit in; might help decide the source?"): CHECKED. phi = SU(2)_3 WZW = Fibonacci anyon
quantum dimension (established, external_literature_connections); phi^6 = lambda = Bogomolny parameter ~17.944
(P7: S_eff=sin^4/[phi^6(1+beta rho)], induced G_N=G_src^2/(4pi phi^6), DE potential e^{-phi^6 rho}). So phi^6 is a
VACUUM/CONDENSATE-director coupling set by the 2I/SU(2)_3 structure -- common to EVERY knot (baryon E_6, dark SO(10),
lepton) and to induced gravity. => a0 = band_K cH0 with m_xi/E* ~ phi^6 is a DIRECTOR-INTRINSIC scale; it does NOT
discriminate the strain SOURCE (baryon vs dark-hopfion). [I]
BUT the source IS decided -- by the DATA, not phi^6: the RAR is TIGHT (dark accel tracks the BARYON distribution
point-by-point). A baryon-sourced director strain (baryons impose the BCs, P7:sec:director_strain) gives exactly that.
Dark-hopfion-sourced strain would track the DARK distribution, NOT baryons -> would BREAK the tight baryon-RAR. So:
  GALACTIC strain (flat curves, RAR, a0) = BARYON-TRACKING (forced by the tight RAR; P7's established mechanism).
  DARK HOPFIONS = a SEPARATE component, for the EARLY (recombination CDM) + CLUSTER mass roles -- different epoch/scale,
  not the galactic RAR source.
So the two "sources" are not competing for one role: baryon-tracking strain does galaxies; dark hopfions do early/cluster.
This resolves the (eee)/(ggg)/(hhh) open "source" question for the GALACTIC case. phi^6's contribution: it makes a0 a
baryon-independent vacuum scale, consistent with (not source-deciding for) the baryon-tracking galactic strain.
Caveat: "dark hopfions could add some galactic strain" is not excluded, only bounded to be subdominant by the RAR
tightness (same logic as the (yy)/(zz) double-count: a non-baryon-tracking component adds RAR scatter). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (jjj) -- SYNTHESIS (user): TWO dark phenomena at TWO epochs, one condensate. Dark hopfion = early CDM + clusters (closes the s9 option-3 gap); strain = late galactic RAR. Reproduces the observed LCDM-CMB/clusters + MOND-galaxies split. Galactic-overlap = the one open tension.
User: "could be BOTH -- dark hopfions (stable, sterile, baryon-independent) dominant EARLY, no [dominant] effect now;
DM research seeks ONE unifying force when it could be TWO at two times." Pursued; it is the framework's actual
structure. Refinement: not "no effect now" (cold a^-3 relics still gravitate/cluster today) but SUBDOMINANT in
galaxies (RAR-tightness-bounded), dominant early + in clusters. And not "two forces" -- TWO MODES of the ONE
condensate: a stable KNOT (dark hopfion) vs a baryon-dragged director TEXTURE (strain). One fabric, two
configurations, two epochs.
EARLY role -- CLOSES s9 OPTION 3 [the payoff]: claude-hopfion s9 names the missing recombination-CDM ingredient as
"non-baryonic, UNSCREENED, NON-ultralight (m>>H_rec), ADIABATIC, a^-3." The D-type dark hopfion is ALL of these:
non-baryonic (D-branch, not a baryon knot); unscreened (a sterile KNOT carrying its own Q_H, NOT the director texture
that S_eff screens); m_dark~5-10 GeV >> H(z_rec); adiabatic (ADM CO-GENESIS -> density correlated with baryons, NOT
Kibble-Zurek isocurvature); cold a^-3. s9 had DISMISSED "a dark Hopfion (rho-knot): does feedback, isn't dark" -- but
that is a SCALAR-rho knot; the D-type is a sterile ORIENTATION knot (no T-matrix closure, Paper 21), so the dismissal
does NOT apply. => the dark hopfion is EXACTLY the component s9 said was missing -> closes option 3 (recombination
Omega_DM), which the strain (LATE z~few, isocurvature) structurally CANNOT. [I/O] conditional on existence
(wall-gated) + Omega~0.26 (ADM) + adiabaticity (co-genesis).
CLUSTER role: dark hopfion at cosmic Omega_dm/Omega_b=5.4 -> M_tot/M_bar~6.4 = the observed cluster ratio (LCDM-like,
fair samples); strain (baryon-tracking, ~5x STELLAR << gas) subdominant in gas-rich clusters. So clusters =
dark-hopfion-dominated, closing the s4 lensing deficit (strain alone ~5x, needs ~28x vs stellar) WITHOUT the
conditional fossil-network channel.
GALACTIC role: baryon-tracking STRAIN gives the tight RAR (hhh: lensing=dynamics); dark hopfion present but must be
SUBDOMINANT in the inner RAR region or it detracks the RAR (the (yy)/(zz) double-count, now sharpened: cosmic-ratio
CDM in galaxies would be ~5x baryons -> CDM-like, not RAR). THE OPEN TENSION. Candidate resolution [I, unverified]:
VELOCITY-DEPENDENT SIDM -- the dark hopfion self-interacts via the shared n_hat (claude-hopfion s4); if sigma/m rises
at LOW v (galaxies ~100 km/s) vs HIGH v (clusters ~1000 km/s), it cores/diffuses in galaxies (subdominant inner,
strain dominates the RAR) but stays concentrated in clusters (dominant) -- the standard vel-dependent-SIDM
phenomenology, here n_hat-mediated. Would deliver galaxy-subdominant/cluster-dominant automatically. Needs sigma(v)
+ halo modeling (not done).
NET -- the "relativistic-MOND frontier" (claude-hopfion s9) resolved as a HYBRID, both components from one condensate:
 * CMB (z~1100): dark-hopfion CDM (adiabatic, cold) -> LCDM-like peaks. [closes option 3]
 * clusters: dark-hopfion CDM (cosmic ratio) -> LCDM-like. [closes s4 deficit]
 * galaxies (z~0): baryon-tracking strain -> tight RAR/MOND. [hhh]
This reproduces the OBSERVED split (LCDM wins CMB+clusters, MOND wins galaxies) that neither pure theory explains --
via two condensate modes at two epochs. Distinct from LCDM (galaxies = strain, not CDM halos -> tight RAR, no
cusp/diversity problem) and from MOND (CMB+clusters = real dark-hopfion CDM). FALSIFIABLE discriminators: (1) the
dark hopfion is sterile ~5-10 GeV ADM (direct-detection-evading, SIDM); (2) galaxy dynamics = RAR with LENSING
slightly > kinematics ~(c_dir/v)^2 (hhh); (3) vel-dependent SIDM core-sizes galaxies vs clusters; (4) recently-merged
clusters denser (s4). OPEN (all pre-existing): dark-hopfion EXISTENCE (wall-gated), Omega~0.26 (ADM conditional),
the galactic-overlap sigma(v). NEXT: the galactic-overlap / vel-dependent-SIDM check (does n_hat-mediated sigma(v)
make the dark hopfion galaxy-subdominant, cluster-dominant?). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (kkk) -- GALACTIC-OVERLAP explored: n_hat-mediated vel-dependent SIDM gives the galaxy(thermalised)/cluster(collisionless) split ROBUSTLY in direction, but TIGHT within sigma/m bounds. Resolvable via massive-mediator n_hat-SIDM; direction solid, exact RAR open.
src_paper7/dm_polarization/galactic_overlap_sidm.py. GAP: a cosmic-ratio dark hopfion, if COLLISIONLESS CDM in
galaxies, forms halos that do NOT track baryons tightly -> would add scatter and DETRACK the tight galactic RAR (the
strain's). RESOLUTION mechanism: the dark hopfion self-interacts via the shared n_hat (claude-hopfion s4); scatterings/
particle N=(rho/m)sigma v t. N>~few -> THERMALISED -> DM relaxes to the (baryon-dominated inner) potential -> TRACKS
baryons (SIDM-baryon connection) -> RAR-preserving; N<~1 -> collisionless NFW -> dominant (clusters, CMB).
RESULT (N across systems, sigma/m=0.5(v_cl/v)^n cm^2/g, cluster-normalised to the Bullet bound):
  split is ROBUST: N ∝ rho/v^(n-1) -> galaxies (dense+SLOW) >> clusters (diffuse+FAST) for any n>=0. Driven by rho/v.
  n=1: dwarf N=17 (sigma/m=20, in SIDM band), spiral N=1.7 (marginal), cluster N=0.1 (collisionless). OK-ish.
  n=2: dwarf N=675 but sigma/m=800 (EXCLUDED, dwarf bound ~tens); spiral N=17 (good); cluster 0.1.
  n=4 (massless n_hat/Rutherford): wildly overshoots dwarf sigma/m (1e6). EXCLUDED.
VERDICT: the galaxy-thermalised / cluster-collisionless split is in the RIGHT DIRECTION and robust (rho/v). BUT a
single power law is TIGHT: n=1 keeps dwarfs in-bound but leaves spirals only marginal; n=2 thermalises spirals but
overshoots dwarf sigma/m. => needs a MASSIVE-mediator (Yukawa) n_hat-SIDM with a turnover (the standard vel-dependent
SIDM that fits dwarf->cluster), NOT the c_dir->0 MASSLESS limit (n=4, excluded). So the n_hat interaction must be a
massive-mediator SIDM, sigma/m(cluster)<~0.5, mild/turnover vel-dependence -- a SPECIFIC, FALSIFIABLE condition.
Builds on (ww): the bare GEOMETRIC knot sigma/m~1e-3 is below the SIDM band, so this REQUIRES the n_hat-mediated
(enhanced, vel-dependent) interaction -- which is exactly the sector's flagged SIDM feature, now with a required
sigma(v) form. OPEN [O]: (1) compute the actual n_hat-mediated sigma(v) (needs the mediator mass; is it Yukawa with
n~1-2? the c_dir/m_xi structure); (2) does thermalised dark hopfion + baryon-tracking strain give the EXACT observed
RAR normalization + 0.05-dex scatter (halo modeling, not done)? NET: the overlap is PLAUSIBLY resolvable, direction
robust, pending the n_hat sigma(v) and the combined-RAR fit. The two-epoch hybrid (jjj) holds together IF the n_hat
SIDM is massive-mediator with cluster sigma/m<~0.5. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (lll) -- n_hat-mediated sigma(v): RESOLVING the anisotropy (user) gives a LOGARITHMIC (vortex) interaction, NOT the isotropized Yukawa. (kkk)'s n=4-excluded was an isotropization artifact. Log -> plasma-like relaxation -> natural galaxy/cluster split. Bell-consistent.
src_paper7/dm_polarization/nhat_mediated_potential.py. User flagged the field anisotropy + Bell. BOTH decisive.
MEDIATOR = the SEMI-DIRAC director mode E^2(k)=v^2 k_par^2 + A^2 k_perp^4 (k_par along e_r = easy/LINEAR, k_perp =
hard/QUADRATIC; claude-hopfion s3), NOT an isotropic massless Yukawa. Static potential V(r)=INT d^3k e^{ik.r}/E^2(k):
inner k-integral ~ C/|k| on both axes -> V(r) ∝ INT cos(kr)/k dk = Ci(kmax r)-Ci(kmin r) = LOGARITHMIC (2D-Coulomb/
vortex), IR-regulated by kmin~m_xi (saturates at r~1/m_xi ~ HUBBLE scale). VERIFIED: dV/dln r -> 1.0 across
1/kmax<<r<<1/kmin. Force F=-dV/dr ~ 1/r (Coulomb-like force).
CORRECTS (kkk): the n=4 (massless-Yukawa Rutherford, EXCLUDED) was an ISOTROPIZATION ARTIFACT (claude-hopfion s10
trap #1). The real n_hat interaction is LOGARITHMIC -- qualitatively different, NOT excluded.
IMPLICATION for the galaxy/cluster split: a log potential (1/r force) gives PLASMA-LIKE relaxation, t_relax ~
v^3/(n C^2 ln L) -- the SAME v^3/n scaling as gravitational two-body relaxation. Dense+SLOW galaxies relax
(thermalise -> baryon-tracking -> RAR-preserving); diffuse+FAST clusters do NOT (collisionless -> dominant CDM). So
the overlap split is NATURAL for the log interaction (even more robust than (kkk)'s tuned power-law) -- IF the
director coupling C makes t_relax(galaxy) < t_age < t_relax(cluster). C (the dark-hopfion<->n_hat coupling strength)
NOT computed [O]; gravitational-strength C would NOT relax galaxies, so it needs C >> gravity (plausible: the dark
hopfion IS a knot in n_hat, strongly coupled, unlike a test mass).
BELL (Papers VIII/IX, user's lead -- apt): the interaction is through the SINGLE SHARED n_hat ("the geometric cost
of the shared geometry"), NOT particle exchange; its long-range, m_xi/Hubble-regulated (globally coherent) log form
is exactly the shared-single-field mechanism (VIII:334,478,633; s7 "structure is global, readout causal"). So this
SIDM is a shared-geometry back-action -- the Bell n_hat at galactic/cluster scale. Conceptually unifies the DM
self-interaction with the Bell nonlocality mechanism (one shared director).
NET: the galactic-overlap (kkk) is PLAUSIBLY resolved, and better than via tuned SIDM: the anisotropy-forced
LOGARITHMIC n_hat interaction gives plasma-like relaxation that naturally thermalises galaxies (RAR-preserving) and
leaves clusters collisionless (dominant). RESCUES from the n=4 exclusion. OPEN [O]: (1) compute C (dark-hopfion<->
n_hat coupling) -> is t_relax(galaxy)<t_age<t_relax(cluster)? (2) the long-range (Hubble-regulated) log is unusual
(cosmological coupling / global coherence, Bell-like) -- needs care it doesn't violate bounds. (3) the exact RAR
from thermalised-dark-hopfion+strain still needs halo modeling. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (mmm) -- PLAN + deferred paper edits. Next: compute the n_hat coupling C -> does the relaxation split galaxies from clusters? Paper additions HELD.
Direction (user): compute (1) the dark-hopfion<->n_hat coupling C and the relaxation-time split t_relax(galaxy) vs
t_relax(cluster); then (2) the long-range/bounds check and (3) the combined-RAR halo modeling ONLY if needed for (1).
PAPER EDITS DEFERRED until the result is solid: the dark-hopfion material (D-type as early/cluster CDM; the
two-epoch hybrid; the n_hat-SIDM / log-interaction resolution) belongs in PAPER 22 (the D-type/co-genesis paper),
NOT Paper 7. Paper 7 would later get FORWARD REFERENCES to the dark-hopfion paper (P22) where the two-epoch division
and the galactic strain's subdominant-dark-hopfion companion are invoked -- added only once (1) [+(2)/(3) if needed]
are solid. Nothing added to any paper yet this sub-thread. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (nnn) -- (1) coupling C: galaxy/cluster thermalisation window C~60-4000 eV (weak, sub-keV; cutoff-set but plausible; split robust). (2) long-range log is dark-only (no baryon fifth force) BUT raises a recombination dark-pressure concern. (3) not needed.
(1) COMPUTED (per-encounter thermalisation criterion for the log interaction (lll): C lnLambda >~ (1/2) m_dark v^2;
cumulative/plasma soft scatters only lower the needed C). With m_dark=5 GeV, lnLambda~10:
  C_thresh = (1/2)m_dark v^2/lnLambda: dwarf 2.5 eV, spiral 62 eV, group 694 eV, cluster 4000 eV.
  WINDOW (galaxy-thermalised, cluster-collisionless): C in ~[60, 4000] eV; width = (v_cl/v_spiral)^2 ~ 64; dwarfs
  thermalise for C>~2.5 eV. Required C ~ 10-1000 eV = WEAK (C/m_dark~1e-8) soft-director-Goldstone coupling.
  C ~ Frank K (the log coefficient) is CUTOFF-SET (not independently pinned, same boundary as a0), BUT the window is
  a modest sub-keV scale and the v^2 split is ROBUST -> galaxy/cluster thermalisation split is ACHIEVABLE for a
  plausible weak C. (1) ANSWERED: resolvable, gated on the cutoff-set C landing in [~60,4000] eV.
  DENSITY note: full t_relax ~ m^2 v^3/(n C^2 lnL) -> also collisionless at RECOMBINATION (low cosmic-mean n, smooth)
  -> dark hopfion clusters as CDM then; thermalisation only switches on in collapsed HALOS. So the three regimes
  (recomb collisionless, galaxy thermalised, cluster collisionless) are consistent via (n,v).
(2) LONG-RANGE CHECK: the log is IR-regulated only at m_xi~H0 -> HUBBLE-RANGE force. GOOD: it is DARK-SECTOR ONLY
  (baryons carry no dark Q_H -> NO baryonic fifth force -> evades WEP/fifth-force bounds). CONCERN [O]: a long-range
  dark-dark force can give the dark hopfion a dark-sector PRESSURE / sound speed; if non-negligible at recombination
  it would modify the CDM clustering -> the CMB peaks (jjj's clean-CDM role). The weak C (~1e-8 m_dark) and the smooth
  low-density recomb medium plausibly keep it pressureless-enough, but this is the genuine remaining tension: the SAME
  long-range log that resolves the galactic overlap must NOT spoil the recombination-CDM role. NEEDS the dark-fluid
  sound-speed estimate at z~1100 (not done). This is the one check that could break the two-epoch hybrid.
(3) halo-RAR modeling: NOT needed for (1)/(2); deferred.
NET: the galactic-overlap is resolvable (log interaction, weak cutoff-set C in a robust window); the two-epoch hybrid
holds IFF the long-range log stays pressureless-enough at recombination (2, open). Paper edits still DEFERRED (mmm):
dark-hopfion/two-epoch -> Paper 22; Paper 7 forward-refs later -- pending (2). [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (ooo) -- (2) recombination sound speed: c_s^2=C/m~1e-8 -> lambda_J<<CMB peaks (SAFE) -- BUT only if the long-range log is SCREENED (else divergent cosmological self-energy). Screening (Debye/Bell) is the deep open. Anisotropy averages out at recombination.
src_paper7/dm_polarization (sound-speed estimate). Dark pressure from the log n_hat force: P_int~C n -> c_s^2=dP/drho
=C/m_dark. For the window C~60-4000 eV, m=5 GeV: c_s^2~1e-8..1e-6, c_s~7-270 km/s, comoving Jeans lambda_J~0.005-0.2
Mpc << CMB-peak scales (>10 Mpc) => the dark hopfion is CDM-like for the PEAKS. (2) PASSES at the sound-speed level;
the SAME weak C does double duty (thermalise galaxies C>~60 eV AND stay pressureless for peaks c_s^2=C/m~1e-8).
DEEPER ISSUE SURFACED (honest -- not a clean pass): the log is IR-cut only at m_xi~H0 -> HUBBLE-range. A long-range
UNSCREENED repulsive force on a cosmological density of same-(dark)-charge hopfions has a DIVERGENT self-energy
(u_int~n^2 C r_max^3, per-particle energy >> m_dark c^2) -- UNPHYSICAL. So the picture is consistent ONLY if the
interaction is SCREENED to a local range. Two candidate screenings, both plausible: (a) DEBYE -- the mobile cold
dark-hopfion medium rearranges to screen the log at lambda_D << Hubble (standard plasma screening) -> local -> c_s^2
=C/m valid, self-energy finite; (b) BELL (Papers VIII/IX, user's lead): the shared-n_hat is a geometric CORRELATION
(boundary condition), NOT a propagating force carrying energy -- "no signal travels" -> it may not produce a
conventional long-range self-energy/pressure AT ALL. (a) and (b) both save it; which (and whether a Bell correlation
can still drive the galactic RELAXATION, which I derived as a force) is THE deep open question now -- the NATURE of
the n_hat interaction (screened force vs shared-geometry correlation) underlies BOTH the galactic thermalisation and
the recombination safety.
ANISOTROPY (user flag for (2)): at recombination the 2I condensate is ISOTROPIC on cosmological scales (uniform
expansion preserves 2I, random local director; claude-hopfion s1) -> the anisotropic log AVERAGES to an isotropic
<C> -> anisotropy is NOT decisive at recombination (unlike galaxies, coherent radial director). c_s set by <C>.
FALSIFIABLE: the UPPER window (C~keV, c_s~270 km/s, lambda_J~0.2 Mpc) suppresses SUB-Mpc power (Lyman-alpha/
satellites) -- a small-scale-suppression prediction distinct from pure CDM (like warm/SIDM), while the CMB peaks stay
CDM-like.
NET: (2) passes at the sound-speed level (hybrid holds for the peaks) IFF the long-range log is screened (Debye or
Bell-correlation; plausible). The hybrid's remaining crux is the NATURE of the n_hat interaction. Paper edits:
dark-hopfion/two-epoch -> Paper 22 is now coherent enough to DRAFT, with the n_hat-nature/screening flagged [O];
Paper 7 forward-refs after. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (ppp) -- USER: density feedback too, not just anisotropy. CHECKED: S_eff=sin^4/[phi^6(1+beta_*rho)] -- the (1+beta_*rho) is NEGLIGIBLE at all diffuse DM densities (unscreened), screens only stars/planets/knot-cores. Doesn't change (ooo); SUPPORTS the hybrid.
User flagged the SECOND factor in S_eff (I'd emphasized sin^4 anisotropy). Sourced P7:rem:lorentz (l.551):
(beta_*rho_sun)^-1 ~ 1e-6 -> threshold 1/beta_* ~ 1e-6 g/cm^3 (below=UNscreened beta_*rho<1, above=SCREENED).
beta_*rho across densities: recomb cosmic mean (3e-21) 2e-15; galaxy halo (1e-24) 7e-19; ISM/disk (2e-24) 1e-18;
cluster gas (3e-26) 2e-20 -- ALL << 1 => UNSCREENED (density feedback NEGLIGIBLE in all diffuse DM). SCREENED only:
star/Solar System (rho~1.4) beta_*rho~1e6; dark-hopfion KNOT CORE (~nuclear 2e14) ~1e20.
IMPLICATIONS: (1) does NOT change (ooo) self-energy/c_s (the long-range log lives in the diffuse unscreened medium)
nor the two-epoch diffuse-density results. (2) SUPPORTS the hybrid: (a) no DM in the Solar System (dense, screened =
Cassini); (b) in dense cluster CORES the director STRAIN screens off (claude-hopfion s4 "wrong sign for the deficit")
but the DARK-HOPFION KNOT (sterile, gravitates as a knot, not director-strain) is NOT subject to S_eff screening ->
provides the core mass -> the two-epoch division (strain=diffuse galactic, dark-hopfion knot=dense cluster core) is
CONSISTENT with the screening; (c) the knot cores being screened means the dark hopfion's dense core decouples from
the director, only its diffuse tail couples long-range. (3) The RAR within a (diffuse, unscreened) halo is NOT shaped
by density feedback (uniform max-DM); the "30 dex planets->dwarfs" (s4) is the screened EXTREME (dense bodies off) +
the a0/acceleration law (diffuse). Both sin^4 anisotropy AND (1+beta_*rho) to be carried in the Paper 22 draft, with
the density factor noted sub-dominant for diffuse DM. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (qqq) -- DRAFTED into Paper 22: new section "The D-type sector as two-epoch dark matter" (sec:twoepoch), both anisotropy + density feedback carried, opens flagged. Paper 7 forward-refs still pending.
papers/main_paper22.tex: added \section sec:twoepoch after sec:abund (+ Paper8/Bell bibitem). Content: (1) two modes
of one condensate at two epochs; (2) relic = recombination CDM (non-baryonic, sterile, m>>H_rec, adiabatic via
co-genesis) + cluster mass at cosmic ratio -- closes the strain's late/isocurvature gap; (3) galaxies: strain RAR +
relic thermalised via the SHARED-n_hat interaction; the semi-Dirac propagator gives a LOGARITHMIC (vortex) potential,
S_eff carries BOTH sin^4 AND (1+beta_*rho) [density factor negligible in diffuse DM, screens stars/knot-cores];
plasma relaxation t~v^3/(nC^2) -> galaxy-thermalised/cluster-collisionless; weak C~10-1e3 eV, c_s^2=C/m~1e-8 safe for
peaks; (4) reproduces LCDM(CMB+clusters)/MOND(galaxies); (5) CONDITIONAL, central open = the NATURE of the n_hat
interaction (screened force vs Bell no-signaling correlation, Paper VIII) + the self-energy; falsifiable: sterile
few-GeV SIDM, sub-Mpc power suppression, lensing>kinematics by O((c_dir/v)^2). Verified: cites/refs/envs balanced,
Paper8 bibitem added. STATE: Paper 22 states it as a conditional program (matches its abstract); the "what changed"
narration lives here, not the paper (s7). Paper 7 FORWARD-REFS to P22 (two-epoch; the galactic strain's relic
companion) DEFERRED until the n_hat-nature crux is settled. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (rrr) -- n_hat-INTERACTION NATURE (Paper VIII grounded): it is a DILEMMA. Bell = Aharonov-Bohm HOLONOMY (no energy -> can't relax); the DM force = CONFINING log (relaxes but OVERCLOSES unless nonlinearly screened). The (lll)-(nnn) n_hat-SIDM resolution is CONDITIONAL, not clean. Paper 22 sec:twoepoch SIDM claim needs a caveat.
Read Paper VIII (source-check). MECHANISM: the shared n_hat does Bell correlations as an AHARONOV-BOHM HOLONOMY -- a
single global, source-fixed BOUNDARY CONDITION (fixed at pair creation); "the correlation is a holonomy of the shared
condensate geometry, NOT a transmitted influence"; no signal, no energy transfer (only the density back-action
delta_rho propagates, at c_s=c/phi). So the Bell regime carries NO energy.
APPLYING to the dark-hopfion overlap resolution (lll-nnn), TWO natures, both problematic:
 (A) HOLONOMY (Bell): no energy transfer -> CANNOT thermalise/relax a halo (relaxation needs energy redistribution).
     Also it is the entangled-PAIR regime (n_hat a boundary condition at pair creation), NOT a dense dynamical medium.
     => if the dark-hopfion interaction is holonomy-like, it does NOT relax galaxies -> the overlap is UNRESOLVED
     (dark hopfion stays collisionless -> detracks the RAR).
 (B) FORCE (director distortion, the (lll) LOG): relaxes, BUT the log is CONFINING (V~C ln r GROWS with r). For a
     dense cosmological medium (N_in_range ~ n R^3 ~ 1e67 in a halo) the cumulative interaction energy ~ n^2 C R^3 ln
     per volume >> rho_dm -> OVERCLOSES. A confining potential cannot be an unscreened cosmological self-interaction.
     => the (lll) log-force OVERCLOSES unless SCREENED to a finite range.
ESCAPES (plausible, unestablished): the director S^2 NONLINEARITY (|n_hat|=1 saturates the distortion -> caps the log
at a saturation length) and/or the director COHERENCE LENGTH (set by the baryon-sourced strain -- coherent in a halo,
random/cancelling cosmologically) could SCREEN the confining-log to a finite range -> finite energy + local relaxation.
If the mode is effectively MASSIVE (Yukawa, falling) instead of gapless-log, it reverts to the (kkk) tight power-law
SIDM window. So the n_hat-SIDM overlap resolution is CONDITIONAL on a nonlinear/coherence screening of the confining
interaction -- the refined crux, NOT the clean resolution (lll)-(nnn) implied.
CONCEPTUAL PAYOFF (grounded, user's Bell lead): the ONE shared n_hat does BOTH -- Bell nonlocality as an A-B HOLONOMY
(source-fixed pairs, no energy) AND the dark-matter self-interaction as a director FORCE (dynamical medium). Same
field, two regimes. The DM is the force regime; its viability hinges on the confining-log being nonlinearly screened.
CORRECTION owed: Paper 22 sec:twoepoch presents the plasma relaxation as resolving the overlap; it should be SOFTENED
to "conditional on a nonlinear/coherence screening of the confining n_hat interaction" (the log overcloses otherwise).
Flag for the next Paper 22 revision (committed; not edited now). RETRACTS the (lll)/(nnn) optimism that the log
"naturally" resolved the overlap -- it does so only if screened, and the holonomy alternative cannot relax at all.
NET: the galactic OVERLAP is the two-epoch hybrid's genuine unresolved weak point; CMB+cluster roles (dark-hopfion
CDM, no competing strain there) stand. [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-11 (sss) -- USER: the Bell paper predicts GEOMETRY-DEPENDENCE (planar/coherent=quantum, isotropic/spherical=classical). This GATES the DM n_hat-interaction and RESOLVES the (rrr) cosmological-overclosure dilemma: active in coherent galaxies, WASHED OUT in isotropic cosmology. Grounded, not ad hoc. Shared falsifiable signature.
Paper VIII's three geometry predictions (verified, P8 Summary): PLANAR (all existing Bell tests) CHSH~2.6 (quantum,
~Tsirelson); BIASED-3D degrades 2.6->2.0 (sigma~30 deg); ISOTROPIC-3D CHSH~1.7 (meas) / ~1.0 (weak) = CLASSICAL. QM
predicts 2.828 for all. MECHANISM: the SAME S_eff=sin^4(theta)/phi^6 anisotropy is RESOLVED in a coherent (planar)
geometry -> strong correlations, but WASHES OUT under isotropic averaging -> classical. (= the "never isotropize"
principle, now a falsifiable Bell prediction.)
CONNECTION (user, resolves (rrr)): the DM n_hat-interaction is governed by the SAME S_eff anisotropy, so it inherits
the SAME geometry-dependence ->
  GALAXIES: baryon-sourced strain makes n_hat COHERENT (radial) -> the anisotropic interaction is RESOLVED/ACTIVE ->
    relaxes the dark hopfions (thermalise -> track baryons -> RAR-safe). (kkk-lll relaxation, now gated-ON here.)
  COSMOLOGY/recombination: director ISOTROPIC/random (uniform 2I expansion, pre-structure, no coherent strain) ->
    the anisotropic interaction WASHES OUT -> effectively INACTIVE -> NO long-range confining-log -> NO cosmological
    overclosure (the (rrr) fatal concern) -> dark hopfion is collisionless CDM (peaks safe, (ooo)).
So the (rrr) COSMOLOGICAL overclosure is RESOLVED: the confining-log is geometry-GATED OFF where the director is
isotropic (cosmology), ON where coherent (galaxies). The "screening" (rrr) needed is the isotropic WASHING-OUT --
GROUNDED in Paper VIII's established geometry-dependence, NOT ad hoc. The gate = the baryon-sourced director coherence.
UNIFICATION: ONE S_eff anisotropy -> Bell violations (planar quantum / spherical classical) AND DM n_hat-interaction
(coherent-galaxy active / isotropic-cosmology washed). One principle, two arenas. SHARED FALSIFIABLE SIGNATURE: the
geometry-dependence is testable in Bell experiments (planar vs spherical CHSH; satellite/3D networks, Paper VIII) AND
underlies the DM behaviour -- a cross-domain test of the same sin^4 structure.
UPGRADES (rrr): the galactic-overlap resolution is now GROUNDED (geometry/coherence-gating = Bell's own prediction),
not the ad hoc screening (rrr) flagged. REMAINING [O] (refined, honest): the WITHIN-halo energetics -- in the COHERENT
galactic regime the log is active + dense (1e67 hopfions); is the cumulative confining energy just the accounted
coherent halo-DM (strain) background + a FINITE perturbative mutual relaxation, or does it overclose within the halo?
(Plausibly the former: the coherent distortion IS the halo DM already counted; relaxation is the perturbation on top.)
Not fully settled. Also the relaxation RATE and coherence length. NET: the two-epoch hybrid's galactic overlap is now
PLAUSIBLY resolved by the grounded geometry-gating; cosmological overclosure resolved; within-halo energetics the
refined open. Paper 22 sec:twoepoch: revise the overlap resolution to the GEOMETRY-GATING (Bell-grounded) + the shared
falsifiable geometry-dependence, rather than (rrr)'s "conditional on ad hoc screening." [[dark-hopfion-cdm-lead]].

## UPDATE 2026-10-12 (ttt) -- Paper IX read (full); sec:twoepoch REVISED to the Bell-grounded geometry-gating. [E] the geometry-dependence is now backed by BOTH Bell papers (VIII kinematic, IX division-algebra); paper edit done.
Read main_paper9.tex in full (2410 lines) to ground the (sss) geometry-gating in the division-algebra Bell paper, not
just Paper VIII's kinematics. CONFIRMS/DEEPENS (sss): IX derives p_cond = cos^4(t/2)/[cos^4+sin^4] = Born rule at the
SQUARED argument (Hopf S^3->S^2, t->t^2); S(theta)=sin^4(theta)/phi^6 is the SAME density-feedback J4 kernel (so the
Bell correlation and the DM n_hat-interaction share one kernel -- not an analogy, the same S_eff). Two ceilings: lepton
p_cond -> 8/3~=2.667 (P9:prop:projective, post-quantum, structural gap 2sqrt2-8/3~=0.161); photon Stokes -> 2sqrt2
(Tsirelson, matches Hensen/Poh). 3 lepton generations <-> C,H,O; C is the unique division algebra supporting composite
entanglement (Tsirelson gap = H->C projection + Fourier mismatch). NET: the geometry-dependence (coherent=strong,
isotropic=washed) is a property of the sin^4 kernel itself, established in both Bell papers -- the (sss) gating rests on
established physics, not one paper's kinematic model.
EDIT (main_paper22.tex): sec:twoepoch "Galaxies" + "Status" paragraphs REWRITTEN. Galaxies: overlap resolution is now
the GEOMETRY-GATING -- the relic interaction inherits the S_eff sin^4 geometry-dependence of the Bell correlations
(Papers VIII-IX), RESOLVED/active where the baryon-sourced strain makes n_hat coherent (galaxy -> relaxes, tracks
baryons, RAR-safe), WASHED OUT where n_hat is isotropic (recombination medium + inter-galactic -> inert, collisionless
CDM, no overclosure). Replaces (rrr)'s "conditional on ad hoc screening/Debye." Status: central open MOVED from "the
nature of the n_hat interaction" to the WITHIN-HALO energetics only (cosmological overclosure resolved by isotropic
washing); added the shared cross-domain falsifiable signature (planar-vs-spherical Bell geometry-dependence probes the
same kernel that gates the dark sector). Added \bibitem{Paper9} (Division Algebras..., doi 10.5281/zenodo.20075227),
reused verbatim from Paper VIII's bibliography per CLAUDE.md s12.
CHAIN: grounds (sss) [[dark-hopfion-cdm-lead]]; closes the (rrr) overclosure dilemma via the (sss) gating now backed by
IX. CAVEATS [O]: within-halo energetics (coherent regime: log active + dense ~1e67 hopfions) -- cumulative confining
energy plausibly = the already-counted coherent halo strain + finite perturbative relaxation, but not proven vs an
S^2-saturation requirement; relaxation RATE and coherence length still uncomputed; C (director coupling) uncomputed
(10-10^3 eV window); c_dir cutoff-set (not pinned). STATUS: two-epoch hybrid's galactic overlap PLAUSIBLY resolved
(grounded geometry-gating), cosmological overclosure resolved, within-halo energetics the single refined open; the
paper now states this as current content. No commit yet (awaiting user).

## UPDATE 2026-10-12 (vvv) -- USER intuition captured: the "clapping hand" = the three-speed separation. Faithful; maps onto the computed lensing/dynamics slip. Two refinements recorded.
USER framing: both MOND and CDM are "right" in different regimes; geometry changes, then light, then baryons move to
fill, and the interconnectedness changes the geometry as it reconfigures -- like a clap, where a distant observer SEES
the clap (light) before HEARING it (sound).
VERDICT [E on the mechanism, I on the framing]: a faithful intuition for the THREE SPEEDS (claude-hopfion s4): light
c > magnitude/Bell mode c_s=c/phi (~0.62c, fast smooth GLOBAL) > director strain c_dir->0 (slow, the halo). The clap's
"see-before-hear" IS the computed lensing-vs-dynamics slip (lensing/dynamics = 1 + (1/4)(c_dir/v)^2,
src_paper7/dm_polarization/lensing_dynamics_slip.py): the LENS responds to the fast-settled geometry, the dynamical
halo DRAGS behind at c_dir. = the earlier "the lens is the lightest thing that follows the geometry, then the halo
drags after" (notes session). So the analogy is not loose -- it names a real, falsifiable framework number.
REFINEMENT 1 (not "both theories right"): the framework is NEITHER MOND nor CDM; it REPRODUCES each phenomenology in
its regime (collisionless clustering early/clusters; RAR/a0 in galaxies) and DIFFERS where they diverge (merger-
dependent cluster mass; planar-vs-spherical Bell). Stronger than "both right in their lane" -- a third thing that
PREDICTS the boundary.
REFINEMENT 2 (the "sound" freezes and is interconnected): unlike a real clap, the slow director reconfiguration does
not pass -- it SETS (c_dir->0, "stress without heat") into a persistent halo (the fossil texture). And "inter-
connectedness" = the Bell shared-n_hat (claude-hopfion s7): the reconfiguration is globally COORDINATED (one shared
field), though back-action still propagates at c_s=c/phi -> NO signaling. The passive distant observer slightly under-
sells the shared-geometry coordination.
PAPER: added a one-sentence physics-register version to Paper XXII sec:twoepoch (fast-settled lens vs c_dir-lagging
halo = origin of the slip); the clapping metaphor stays in the notes only. [[dark-hopfion-cdm-lead]]

## UPDATE 2026-10-12 (www) -- WITHIN-HALO energetics (the refined open from (sss)/(ttt)) PLAUSIBLY RESOLVED: Debye screening makes it extensive (~N), not N^2. Mechanism identified, scaling checked.
Script: src_paper7/dm_polarization/within_halo_energetics.py (scaling, not N-body; N~1e68 unsimulable). The worry: in the
COHERENT galactic regime the dark-hopfion director log-interaction is ACTIVE and dense -> does the cumulative energy go
as N^2 (over-close) and need nonlinear S^2 saturation? ANSWER: NO, it splits into three finite pieces --
 (1) MEAN TEXTURE = the halo = ONE splay disclination, |grad n|~1/r -> E=(chi/2)int|grad n|^2 ~ LINEAR in R (M_enc prop
     r, Paper VII). Winding is TOPOLOGICAL (~1), NOT N -> adding relics never makes the halo mass ~N^2.
 (2) FLUCTUATION PLASMA (the relic-relic log that thermalises them): a THERMALISED plasma DEBYE-SCREENS the log ->
     Yukawa-of-log ~K_0(r/lambda_D), each relic sees only N_Debye~n*lambda_D^3 neighbours -> correlation energy
     EXTENSIVE (~N, a renormalisation of m_dark), NOT N^2. CHECKED: for the stated C~10-10^3 eV window, lambda_D ~ 1e-9
     kpc << R_halo=100 kpc (lambda_D/d_ip~1e9-1e11, but << R/d_ip~4e22) -> screening is EFFECTIVE inside the halo;
     E_screened/Mc^2 ~ 1e-32..1e-28 (negligible) vs E_unscreened(N^2)/Mc^2 ~ 1e53..1e57 (absurd). The ratio = N/N_Debye.
 (3) VIRIAL: log is marginally confining, r V'=C^2 const -> 2<T>=C^2/pair -> scale-invariant equilibrium (CORED halo),
     interaction virialised to ~M v^2, no runaway.
VERDICT: option (a) -- NO over-closure, NO S^2 saturation needed. UPGRADES (sss)/(ttt): the within-halo open is now
"plausibly resolved via Debye screening in the thermalised coherent-regime plasma." REMAINING [O]: C (hence lambda_D)
cutoff-set; and the resolution is CONDITIONAL on quasi-neutrality of the FLUCTUATION director-charges about the mean
texture (assumed, not proven -- if the fluctuations were net same-sign, the N^2 boundary term returns). Thermalisation
itself is supplied by the (sss) geometry-gating (coherent -> interaction active -> relax). [[dark-hopfion-cdm-lead]]

## UPDATE 2026-10-12 (xxx) -- DARK-HOPFION EXISTENCE (Paper XXI open-(ii)): pinned the TARGET KNOT TYPE (prism/D_n-graph-link) and laid out the existence logic. Status: ALLOWED, not FORCED. One script bug caught+fixed.
Script: src_paper21/dark_hopfion_existence.py (group/topology facts checked; existence logic stated, not a sim -- the
free-flow route is exhausted/wrong-instrument, notes (j)-(l)). TWO tractable advances:
 A. KNOT TYPE [N proposal]: the A-type ladder is cyclic Z_q -> lens space -> T(2,q) torus knot (nu/lepton/baryon =
    T(2,1/2/3)). The D-type (McKay D-hat) Kleinian singularity C^2/2D_n has LINK S^3/2D_n = a PRISM MANIFOLD, Seifert
    over the orbifold S^2(2,2,n), whose plumbing is the D_n Dynkin TREE (chain + two-prong fork), NOT a torus knot ->
    confirms "off the T(2,Q_H) ladder." (2,2,n) is NOT pairwise coprime (gcd(2,2)=2) -> not a homology sphere -> H_1
    nontrivial. CHECKED: dicyclic abelianization H_1(2D_3)=H_1(2D_5)=Z_4 (n odd) CONTAINS a Z_2 = the sector's own DARK
    PARITY. So the dark parity is literally the torsion homology of the prism manifold. This is the TARGET a charge-
    conserving existence search should aim at.
 B. EXISTENCE LOGIC: PURE Faddeev (sum form) -> Vakulenko-Kapitanskii E>=c|Q|^{3/4} + Lin-Yang (2004) => minimisers
    EXIST in nontrivial Hopf classes (existence in principle). DENSITY-FEEDBACK (product E_fb=K_fb*J4) -> SADDLE not min
    (feedback_saddle_core) -- BUT the BARYON/trefoil has the same saddle structure and EXISTS (proton, confinement-
    stabilised), so "saddle != doesn't exist." The sterile dark hopfion has NO confinement -> its STABILISER is the open
    physics: (i) lambda not pinned to phi^6 (no Q_group) -> possibly LARGE -> a FAT genuine minimiser (notes (k)); (ii)
    its Z_2 parity + prism/graph-link rigidity may protect a saddle. Both plausible, unforced.
VERDICT: existence ALLOWED (nontrivial Faddeev class -> Lin-Yang minimiser for the pure functional), NOT FORCED (index
probe negative; feedback-saddle needs a sterile stabiliser; free flow is the wrong instrument). Numerics need a
Q_H-CONSERVING/constrained high-res search (rational_map_ic.py smooth base + saddle-score, GPU-gated). Paper XXI's
"existence open, stability settled" stands; now with a NAMED target topology + the stabiliser question as the handle.
BUG CAUGHT (CLAUDE.md s1): first draft used divisibility (12|24) to test 2D_3 in 2T -> printed "True", contradicting
the text; fixed to the real criterion (A_4 has no order-6 subgroup -> D_3,D_5 both out). [[dark-hopfion-cdm-lead]]

## UPDATE 2026-10-12 (yyy) -- Cross-referenced the PRIOR existence thread (e)-(m) into (xxx): the knot-type proposal ANSWERS a long-open item; stability REFINED; fat-soliton route reinforced. dark_hopfion_existence.py extended.
Read back the full prior existence thread (notes (e)-(m), 2026-10-05). Four developments folded into (xxx) /
src_paper21/dark_hopfion_existence.py:
 1. KNOT TYPE -- the (xxx) prism/D_n-graph-link proposal ANSWERS the "knot type? (distinct from unknot/torus/trefoil)"
    question that (f),(h) repeatedly posed and left OPEN. Not a restatement -- a resolution of a standing item.
 2. SYNTHESIS (new): the D_n-graph-link/prism topology is intrinsically MORE DISTRIBUTED ("fatter") than a thin T(2,q)
    tube -> it IS the large-lambda/large-size regime of the fat-soliton route (k). So knot-type (A) and fat-soliton
    existence (B) are ONE proposal, mutually reinforcing: the proposed topology is exactly the regime that sidesteps
    the thin-tube resolution wall.
 3. EMPIRICAL MORSEL folded in (notes (l)): the fat (lam=30) D_3 free-flow run stayed LOCALIZED + 3-fold-ROBUST
    (A3~0.37-0.47, never symmetrized to axial); only the lattice topology-leak (which drains the KNOWN trefoil too)
    ended it. => the dihedral config is dynamically robust, does NOT want to go axial; the drain is method, not
    instability. Positive support for existence (not a proof).
 4. STABILITY REFINED (reconciles (xxx)'s "settled" with the prior (e)->(f) correction): unwinding = topological,
    settled; DECAY splits into (1) the SIDEWAYS-BRANCH reconnection barrier [SOLID: 2D_3,2D_5 NOT subset 2T -> decay
    needs a barriered T(2,2)->T(2,3)-type reconnection, not continuous relaxation -> metastable; contrast Q_8=2D_2
    subset 2T = unstable] and (2) the Z_2 dark parity [UNPROVEN as a dynamical symmetry -- (e) over-claimed it, (f)
    corrected: an outer automorphism is not automatically conserved]. So honest stability = "unwinding-settled +
    decay-metastable (solid barrier) -> absolute IFF Z_2 dynamical (open)."
 5. CONTEXT: the TRANSIENT-wake route (m) was explored and RETIRED by the stable-knot picture (w) (the ADM dark NUMBER
    needs a conserved relic); so stable existence is the correct target -- don't resurrect transient.
PAPER XXI EDIT [DONE, user approved]: open-problem (ii) and the abstract RE-ORDERED. Stability now splits into
unwinding (topological, settled) vs decay: the SIDEWAYS-BRANCH reconnection barrier (solid -> metastable; 2D_3,2D_5
NOT subset 2T, ref P21:rem:dtype_slot + Paper XVIII T(2,2)->T(2,3) wall) is foregrounded, and the Z_2 dark parity is
presented as the additional route to ABSOLUTE stability, explicitly conditional on the Z_2 being a conserved dynamical
symmetry "rather than merely an outer automorphism, which we do not establish here." Abstract clause softened to match
("topological against unwinding and barrier-protected (metastable) against decay, absolute only if its Z_2 parity is
dynamical"). Fixes the over-claim that the prior (f) correction had already flagged internally.
STATUS: existence ALLOWED-not-forced with a NAMED target topology (prism/D_n-graph-link, Z_2 = its torsion homology)
that doubles as the fat-soliton regime; stability = unwinding-settled + decay-metastable(solid) + absolute-iff-Z_2-
dynamical. Numerics still need a Q_H-conserving high-res search (GPU-gated). [[dark-hopfion-cdm-lead]]
