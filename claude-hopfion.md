# Hopfion Framework — Geometric Base (read before re-deriving foundations)

Stable, load-bearing facts about the framework's geometry, so they are not
re-worked each session. Epistemic tags per CLAUDE.md §6:
**[E]** established/verified · **[I]** interpretive stance of the framework ·
**[N]** derived this program, in `notes/` (not yet in a paper) · **[O]** open.
Verify against the cited paper before building new results on any single line
(CLAUDE.md §1). Sources are Paper VII (`gravity/`) unless noted.

## 1. Ontology — one structured condensate
- **[I]** Monistic geometry: the vacuum **is** a structured condensate ("the
  fabric"). There is no separate metric/substrate — the fabric is the thing the
  knots are made from; its geometry changes through its own self-interaction
  (density feedback). Gravity is a function of clustered density, not a fixed
  always-on background.
- **[E]** The vacuum carries **binary-icosahedral (2I)** structure. Under
  isotropic cosmic expansion the 2I symmetry is **preserved** (uniform scaling;
  only the scale grows) — no elongation/breaking, which would need anisotropic
  stretch (Phase 0, `notes/fabric_eos_plan.md`). **[N]**
- **[E]** One connected condensate **from the start** is the founding premise and
  the reason for the hopfion. The cosmology **inherits** this; it is not a new
  assumption (Paper VIII/IX; see §7).

## 2. Mass and temperature are geometry (emergent, not preconditions)
- **[I]** **Mass is output, not input.** `E=mc²` = **coherent geometric stress**:
  a stable, localized, conserved geometric configuration (a knot/texture). Do
  **not** treat "rest mass" as a precondition to look up — it is the downstream
  signature of a coherent lump. When testing whether something is matter-like,
  ask about the *configuration*, not a band-structure gap.
- **[I]** **Temperature = incoherent geometric stress** — disordered fluctuation
  energy that forms no stable knot (distributed, releasable). Cosmic evolution =
  geometric stress cohering from hot/soup to cold/structured as the fabric
  dilutes. (Working interpretation; no geometric *redefinition* of T is in a
  paper yet — `notes/recombination_cdm_brainstorm.md`.)
- **[E]** Cosmological-scaling equivalence used repeatedly: **a⁻³ (matter) ⇔ a
  temperature-INDEPENDENT, number-conserved coherent configuration** (fixed rest
  energy, or a field with m≫H). Anything whose energy scale rides the bath
  (∝T∝1/a) is **a⁻⁴ radiation**; anything frozen is **const (DE)**. So "does it
  cluster as CDM" is a question about temperature-independent coherent lumps, not
  about a mass parameter. **[N]** derivation in `notes/fabric_eos_plan.md`.

## 3. The semi-Dirac dispersion (the anisotropy engine)
- **[E]** Paper VII `P7:eq:semidirac`:
  `E² = A²k⁴ sin⁴θ + v²k² cos²θ`, from `H = A k_⊥² σ_x + v k_∥ σ_y`,
  `k_∥ = k cosθ` (along n̂=ê_r), `k_⊥ = k sinθ`.
  - **easy / radial / light axis:** linear, `E = v k` → relativistic.
  - **hard / perpendicular axis:** quadratic, `E = A k²` → "heavy" = vanishing
    **group velocity**, **GAPLESS** (no rest-mass term; E→0 as k→0). "Heavy" is
    NOT a mass gap.
  - The heavy band carries the `sin⁴θ` of the directional suppression
    `S_eff = sin⁴θ / [φ⁶(1+β_*ρ)]` (Axiom A2). The same `sin⁴=(sin²)²` structure,
    via the pentagon `k+2=5`, fixes the **three fermion generations**.
- **[I]** This anisotropy does double duty cosmologically: radiation pressure /
  photon free-streaming run along the **easy-radial** axis; a configuration that
  avoids that pressure sits **hard-perpendicular**.

## 4. Director sector = dark matter (soft mode)
- **[E]** DM lives in the **orientation (director) field n̂** (Axiom A1), the soft
  mode — not in the scalar magnitude ρ. DM = **director strain**: the Frank
  elastic energy `F = ½∫[K_splay(∇·n̂)² + K_bend(n̂×(∇×n̂))²]` of the distortion.
  It gravitates (`E=mc²`) but has no EM coupling → invisible. Being **baryon-
  sourced**, the strain **tracks baryons** automatically (recovers the radial-
  acceleration relation with no tuned halo).
- **[E]** Stiffnesses have **opposite signs**: `K_bend>0` (stable),
  **`K_splay<0`** (unstable) → uniform director spontaneously develops a radial
  **splay/disclination** texture → `|∇n̂|~1/r`, `u~K/r²`, **`M_enc ∝ r`** → flat
  rotation curves (the signature of the splay-driven disclination; a defect-free
  director gives the Keplerian monopole).
- **[E]** **Coldness** is dynamical, not a gap: director speed
  `c_dir² = K/χ`; in the semi-Dirac medium χ (rotational susceptibility) diverges
  faster than K → **`c_dir → 0`** → cold, pressureless, no Jeans obstruction →
  clusters on all scales. (Splay channel: `c_dir²<0` = the instability.)
- **[E]** **`m_ξ ~ H₀`** — the director mass is ultralight (dark-energy sector).
  Consequence: at recombination `m_ξ ≪ H(z_rec)` → **Hubble-frozen, w=−1, DE-
  like** → it does NOT cluster then; the director-strain DM is **LATE (z~few)**.
- **[E]** **Acceleration scale:** strain overtakes Newton at `~cH₀`, matching MOND
  `a₀≈cH₀/(2π)`; the ratio (cutoff / semi-Dirac crossover) is of order **φ⁶** (the
  same normalization as `S_eff`). Magnitude is condensate-derived (from ρ_cond and
  m_ξ) to ~an order of magnitude, no galactic input; exact coefficient not derived.
- **[E]** **Polaron mass — resolves the `K_splay<0` sign paradox** (`P7:rem:polaron`).
  The gravitating/dynamical DM mass is the **effective inertial mass** of a baryon
  dressed by its comoving director distortion:
  `m_eff = (χ/3)∫|∇n̂|² d³x` — **POSITIVE by construction** (χ>0, integrand a sum of
  squares), so **the sign of K_splay never enters**. `K_splay<0` drives the texture's
  FORMATION (`c_dir²<0` in the splay channel), not the sign of the halo mass;
  formation and gravitating mass are cleanly separated. `|∇n̂|²~1/r²` → `m_eff∝r` →
  `M_enc∝r`, positive. **Drag** = the same dressing dynamically: the halo gives the
  standard gravitational **dynamical friction**; a further **elastic** drag (finite
  c_dir, set by the rotational viscosity γ_rot) would scale `F~V`, distinct from
  Chandrasekhar `1/v²` — a CDM-distinguishing signature, unpinned (γ_rot not fixed).
  [Was the A1 static-energy sign tension; `notes/single_galaxy_winding_analysis.md`.]
- **[I]** **CDM ontology — four states of the one condensate** (`P7:rem:cdm_ontology`):
  **baryon** = coherent **knot** (Hopfion, carries Hopf charge, EM-coupled); **dark
  matter** = coherent **deformation** of the director carrying **no Hopf charge** — a
  texture of the medium, not the knot — persisting via `c_dir→0` after the sourcing
  baryon moves on ("**stress without heat**"); **radiation** = *incoherent* thermal
  condensate excitations; **dark energy** = frozen vacuum. DM is not a particle but the
  cold, coherent, EM-neutral **deformation field** of the fabric. (Hopf charge here is
  relative to the ground state — §8.)
- **[E]** **Density-feedback screening → the RAR across ~30 dex** (`S_eff` denominator
  `(1+β_*ρ)`, §3): high-ρ → `S_eff→0` → director OFF → **no DM** (planets, Solar
  System, and cluster **cores**); low-ρ → max director → **max DM** (dwarfs); one law
  spans planets→dwarfs and gives **cored** (not cuspy) profiles. The same feedback that
  turns DM off in dense regions **screens clusters** — the wrong sign for the cluster
  deficit — so the galaxy/RAR mechanism is screening, but the **cluster** fix is the
  topological channel below, NOT screening.
- **[I/E]** **vs MOND** (`P7:rem:bullet`, `P7:rem:defect_network`): MOND = a force law
  slaved to the *local* baryonic `g_N` with universal `a₀`; the framework = a physical
  field with **inertia, topology, and decoupling** that merely *coincides* with MOND in
  the deep-`a₀` galactic regime. Same `a₀`, but the framework can carry mass the
  local-`g_N` law structurally cannot (next bullet). Bullet-cluster lensing/gas offset
  is the *expected* outcome (each galaxy's splay halo rides through collisionlessly;
  shocked incoherent gas seeds no splay) — `P7:rem:bullet`.
- **[N]** **Cluster-mass deficit — the topological disclination-network channel**
  (`P7:rem:defect_network`; `notes/single_galaxy_winding_analysis.md`). Baryon-tracking
  strain (RAR, ~5× stellar) underestimates cluster lensing (~28×). Framework-specific,
  baryon-INDEPENDENT recovery: hierarchical/turbulent merger assembly traps disclination
  lines (**Kibble–Zurek**), frozen by `c_dir→0` → a persistent **fossil network**.
  `M_def/M_RAR ~ (R_cl/ξ_net)²`; recovers when the network freezes at the injection
  scale `ξ_net~ξ_turb~½R_cl` (subcluster merger scale), no free normalization beyond
  the `a₀` anchor.
  - **CONDITIONAL [O]**: 3D inertial sim (`defect_network_sim_3d.py`) — as `c_dir→0` the
    KZ tangle freezes at `c=ξ_net/ξ_turb→1` (recovers), but there is **no true plateau**:
    the defect-free ground state always wins eventually; `c_dir` sets only the coarsening
    **RATE** (`t_coarsen ~ 1/c_dir^{1.5–2}`, diverges as `c_dir→0`). "Frozen" is a
    **timescale** statement. Survival over a cluster age needs `c_dir ≪ ξ_turb/t_cl ~
    50 km/s` — the **SAME coldness** the sector already requires (§4). Recovery holds IFF
    `c_dir` is below that threshold; `c_dir` is cutoff-set, NOT derived. (Overdamped
    `c=2.77` and 2D `c≈2` were wrong-dynamics artifacts — the framework is INERTIAL,
    `χ∂²_t n̂=K∇²n̂`, not overdamped/diffusive.)
  - **[E] Three speeds — don't conflate**: (i) **signal** `c_dir=√(K/χ)→0`; (ii)
    **healing/annihilation** rate = a **mobility** (line tension K vs rotational
    viscosity γ_rot + inertia χ) — *this*, not the signal speed, sets network survival;
    (iii) **magnitude/Bell** mode `c_s=c/φ` (fast, smooth, global). **Topological
    protection**: a defect heals only by annihilation or core-melting, so the fast smooth
    (Bell "one geometric object", §7) reconfiguration heals smooth strain but **cannot
    annihilate a disclination** — the fossil network is safe from it. γ_rot unpinned.
  - **[E] Falsifiable discriminator**: because the network slowly coarsens,
    recently-merged/disturbed clusters carry a **denser** network → **more** mass excess
    than relaxed clusters of equal baryon content. MOND: none; ΛCDM: excess tracks the
    halo, not the merger/relaxation state.

## 5. Magnitude sector = dark energy; scales
- **[E]** Scalar-magnitude perturbation sound speed **`c_s = 1/φ`** is the **dark-
  energy** (magnitude-mode) speed (`P7:rem:cs_de`) — relativistic; **NOT** dark
  matter. Do not attribute `1/φ` to DM.
- **[E]** **`Λ_cond = T_CMB (π²/150)^{1/4} ≈ 1.19×10⁻⁴ eV`**; condensate energy
  `ρ_cond ~ 10 Λ_cond⁴` (Paper XII, fixed by T_CMB). `Λ_cond ~ T ~ 1/a`, so the
  condensate bulk tracks **radiation (a⁻⁴)**; the frozen attractor value is the
  **CC / dark energy (const)**.

## 6. Frames — Jordan vs Einstein (avoid the recurring trap)
- **[E]** Physical (Einstein-frame) DE scalar mass is **`m_ξ ~ H₀`**; the Jordan
  `m_ξ² = φ²⁰ Λ_obs M_Pl⁴` (~213 H₀) is cancelled by strong NMC (`1+3φ⁷≈88`).
  Use the physical mass for cosmological dynamics.
- **[E]** **φ-power discriminator:** *physical/topological* powers are real
  (WZW/Bogomolny `λ=φ⁶`, `α⁻¹=360/φ²`, `ω_BD∝φ¹²`); *scalar-tensor-norm* powers
  are **frame artifacts** (`m_ξ φ²⁰`, `r_V∝φ⁻⁸ᐟ³`, `m_g∝φ⁻⁸` — Jordan-only, no
  observable depends on them); `Λ_UV~φ⁶=4π√2` is a ~1% coincidence.
  Cassini is secured frame-robustly by the **chameleon**, not by r_V.
- **[E]** Induced gravity is **flux-primary**: `G_N = G_src²/(4πφ⁶)`;
  `Λ_UV = M_Pl/φ⁶ ≈ 0.0557 M_Pl` is a **tower** relation (gap = 3 steps = φ⁶), and
  `Λ_UV` sub-Planckian is load-bearing (inflation `P_s ∝ (Λ_UV/M_Pl)⁸`).

## 7. Nonlocality without signaling (Bell mechanism) — Papers VIII, IX
- **[E]** Paper VIII (`main_paper8.tex:334, :478, :633`): the condensate
  orientation **n̂ is a single SHARED field quantity**, a boundary condition fixed
  at pair creation. Correlations arise because both sites probe the **same
  geometry** — **no signal travels**; back-action propagates at `c_s=c/φ`; the
  currency is "the **geometric cost** of the ℍ→ℂ reduction." Realistic +
  deterministic at field level. Apparent nonlocality = shared geometry, not FTL.
- **[E]** Paper IX "**post-quantum window**" (`P9:rem:nosignaling`) = the untested
  **lepton** Bell CHSH ceiling (no experiment to date) — **not** a causality
  hedge. The no-signaling mechanism itself is established for the quantum regime.
- **[I]** Cosmological use: a globally-correlated fabric configuration + c-limited
  visibility (the "smoke-ring") is the SAME mechanism (would also address the
  horizon problem). Structure is global; the readout is causal.

## 8. Sector ladder (charges) — see Papers XV, XVI, XVII
- **[E]** Ladder, **verified against Paper XVII's abstract** (l.46–60) and the
  Paper XV/XVI titles:
  - `Q_H=0` — the **bare condensate**: no winding, no topological charge; the
    structured space, the condensate at rest.
  - `Q_H=1` — **neutrino** sector (Paper XVII, l.91, 216): the lightest
    topological excitation once embedded in the density-feedback medium. This is
    the object the standard Faddeev–Niemi literature calls the *minimal*
    Hopfion, realised here as the lightest **fermion**, not as the ground state.
  - `Q_H=2` — **charged-lepton** sector (e/μ/τ tower): torus, **2I**, WZW
    `SU(2)_3`, `Q_group=10`. Stabilised **because of** the density-feedback
    term, not despite it. **[E]/[I] Geometry (user, 2026-09-02): TWO TUBES —
    an INNER and an OUTER tube** — the fusion of two single-tube `Q=1` quanta,
    `j_½×j_½ = j₀⊕j₁`, into the `j₀` singlet (Paper I `P1:thm:vacuum` (iii),
    l.1318). The **single tube is unstable and collapses** (`Q=1` "not a stable
    ground state ... decays to `Q=2`"; also l.224, l.1599). So the physical
    `Q=2` object is intrinsically two-tube, not a doubled single tube.
  - `Q_H=3` — **baryon** sector (three-quark-like, distinct from the lepton
    tower): trefoil `T(2,3)`, **2T**, `(E_6)_1 ⊃ SU(3)_1^{×3}`, `Q_group=3`
    (Papers XV/XVI). **[I] (user, 2026-09-02):** the two objects that combine are
    the **two pre-images (Hopf fibers) reconfiguring into the `T(2,3)` trefoil**;
    and **Hopf-charge conservation topologically forces a neutrino-like
    byproduct** — identified in the hadronic context with the **pion** (Paper
    XVIII Thm 5.3, used in Paper XIX l.1489, 1583).
  - Gauge bosons/photon appear as `Q=4, J=1` composites of two `Q=2` (Paper
    XIII) — not re-verified here.
  - Physically `Q=2` is genuinely **two tubes**
    (inner + outer) — `P1:thm:vacuum` (iii) fuses two `Q=1` single tubes into the
    `Q=2` singlet and the single tube is unstable/collapses — while the footnote
    at l.323–326 is only a **solver-measure convention** (the `z→−z` bilateral
    factor in `vol=2π·2r·h²` *labels* a single tube `Q=2` in the axisymmetric
    code). Convention vs physics, not a contradiction. **What REMAINS open is the
    cross-paper label:** `Q_H=2` names the **two-tube condensate vacuum** in I/VII
    but the **charged lepton** in XV–XVII (both carry the same 2I / `SU(2)_3` /
  `j_0` data), and `Q_H=0` (no winding) vs Paper I's winding-carrying vacuum.
  — `Q_H` is the Hopf charge
    RELATIVE to the condensate ground state.** `Q_H=0` = the wound `Q=2` vacuum
    (Paper I Thm 6.1); `Q_H=1,2,3` = ν/ℓ/baryon excitations. Edits: XVII
    definition rewritten (l.278–307, `n→n_vac` at infinity, `Q_H=Q_tot−Q_vac`);
    Paper I remark `P1:rem:qh_relative`; Paper VII note in `P7:rem:cdm_ontology`;
    XVI l.1047/2052 clarifiers; XV footnote.
  - **[N]** **Quark-sector check (2026-09-02): the reconciliation is a NAMING choice
    with NO physics at stake.** Whether the ground state carries winding does **not**
    change anything in Papers XV–XIX: XVII l.138–141 states each sector `Q_H=N` gets
    its identity from a **triple** (McKay group, WZW CFT, knot type), and every quark
    result flows from that triple — colour exclusion (group coprimality `gcd(2,3)=1`),
    generations/masses (coset `M(5,6)` CFT), confinement/fractional charge (trefoil rep
    theory), `Q_group{3,6,10}`, `Q_coset=h(E_8)=30` (knot `(p,q)=(2,3)` + Brieskorn
    `Σ(2,3,5)`) — **never** from the absolute integer or the base being winding-free.
    `Q_H` and physics-winding are already decoupled (`Q_H=3` but trefoil `(p,q)=(2,3)`;
    `Q_group≠Q_H`). The CFT `j_0` the physics uses is the **identity primary**, not the
    geometric ground state, so a winding-carrying vacuum never enters it. Read `Q_H` as a
    **sector index / excitation count above the condensate** and Paper I (2I/`Q=2`
    winding-carrying ground state) and XV–XIX are consistent. 
  - **Q_H=1 paper, XVII, close read** it is not *purely* a label —
    XVII's topology section carries an **explicit trivial-vacuum assumption**:
    `Q_H` is defined as the **absolute** Hopf charge (π₃ generator, l.278–284) with
    the finite-energy BC "**n** approaches a **constant** at spatial infinity"
    (l.305–307). That constant-at-infinity IS the winding-free background, and it is
    the **opposite framing** from XVI l.1047 ("vacuum = `Q_H=2`") — a real framing
    inconsistency between the two ladder papers. Still **physics-inert**: redefining
    `Q_H` relative to a wound background gives the same ℤ with the same knot
    generators, changing only the absolute integers. Likely **no contradiction** at
    all — XVII's "constant at infinity" is the *local* soliton reference; Paper I's
    `Q=2` is the *global* condensate texture; restating the BC as "approaches the
    condensate ground state at infinity" (a relative Hopf charge) reconciles them.

## 9. Standing cosmological status (the recurring "CDM gap")
- **[I]** **Frame the gap as a REGIME, not a missing particle.** CDM is a theory;
  "cold particle" is its construct (50 yr null direct detection). The data are
  gravitational (rotation curves, lensing/Bullet, CMB peaks, LSS, BAO); "3rd peak →
  Ω_cdm" is a GR-internal interpretation, not a datum. The framework's Frank-elastic
  mechanism *already* reproduces the galactic data particle-free (flat curves, RAR/
  a₀, no BAO shift) but is a **low-density/low-accel** mechanism that **screens to
  ≈GR at recombination** (phaseA) — so it goes quiet exactly where the 3rd peak
  lives. This is the **relativistic-MOND frontier** (nails galaxies, hard at the CMB
  peak), a shared problem of all medium/elastic DM theories — OPEN, not falsified
  (Paper VII hedges "not fixed here"). Real target: extend the elastic mechanism
  into the high-density regime (an unscreened, recombination-active piece), yielding
  a falsifiable CMB prediction that may DIFFER from ΛCDM — don't chase Ω=0.26.
- **[E]** Framework **derives**: late-time galactic DM (director strain, z~few),
  dark energy (frozen vacuum, Λ_obs), baryonic matter (Hopfions, Ω_b~0.05),
  radiation (thermal condensate).
- **[O]** Framework does **NOT** derive the **recombination cold-matter budget**
  (`Ω_DM≈0.26` at z~1100). All internal routes are exhausted at the calculation
  level (`notes/fabric_cmb_hypothesis_and_plan.md`): route 1 (isotropic bulk,
  radiation-like), route 1′ (no temperature-independent coherent non-baryonic lump
  at recombination — textures order too late), route 2 (local fifth-force ~89×
  weak). **Current honest status = OPTION 3: `Ω_DM` is an input the DM sector does
  not explain at recombination.** 
- **[E]** **Why the gap is structural (do not re-litigate).** Paper VII's own IC
  section (`P7:rem:reheating_IC`, `P7:sec:dm_ic`) already states it: the two cold
  director-sector components each fail for recombination — the **baryon-sourced
  strain** is adiabatic but tracks baryons (only Ω_b≈0.05), and the independent
  **splay texture** is ultralight → Hubble-frozen → "orders only at z~0, absent at
  recombination," *and* is **isocurvature** (peaks need adiabatic). The "no
  feedback" that exists (`P7:rem:lorentz`, l.549) is high-ρ **screening**
  (S_eff→0) — it turns the director *off* in dense regions, the wrong sign for
  clustering. The missing ingredient is a non-baryonic, unscreened, **non-
  ultralight** (m≫H_rec), adiabatic orientation/vacuum configuration; the sector's
  only mass scale is m_ξ~H₀, welded to the DE sector. The feature that unifies
  late-DM + DE (ultralightness) is exactly what forbids an early cold component.
  A "dark Hopfion" (a ρ-knot) does feedback and isn't dark; a smooth vacuum
  configuration is DE-like (const) or radiation, not a⁻³. **So option 3 is the
  paper's own position, for a structural reason — not an unexplored lead.** Full
  chain in `notes/fabric_cmb_hypothesis_and_plan.md`.
- **[I]** **Two-cooling picture (one-line frame for the gap):** as the one
  condensate cools it holds structure in stages — the **charge sector** locks up at
  recombination (z≈1100, α-governed: atoms bind, photons decouple → Ω_b≈0.05), the
  **director sector** long-range-orders only much later (z≈0–few). CORRECTION: Λ_cond
  TRACKS T (Λ_cond=0.506 T_CMB), so T/Λ_cond≈2 at ALL epochs — the medium is
  perpetually MARGINALLY disordered (near-critical), NOT "2000× too hot" at
  recombination (an earlier error that used Λ_cond today). Long-range ordering is
  Hubble-gated (m_ξ~H₀→z≈0), not thermal. So transient/pretransitional order IS live
  at recombination — but the baryon-sourced strain ENERGY is capped at ρ_cond (~ρ_b,
  radiation-like) while CDM needs ~5 ρ_b (Phase 0.1, `phase0_strain_energy.py`: ≥4.5×
  short). The gap survives on ENERGY, not on temperature. The sectors are related at
  the **root** (one 2I/WZW knot — Paper III, "two constants from one knot"), but
  root-related ≠ one supplies the other's scale: α is a dimensionless 2I angle,
  the director mass a dimensionful scale welded to DE. **[E]** the welding is
  structural and sharp: the DM/orientation sector's **entire UV cutoff** is
  m_ξ~H₀, and m_ξ²=φ²^Q·Λ_obs·M_Pl⁴ (Q=10) anchors it to Λ_obs at the field-
  equation level — the whole sector lives at ≲H₀, so there is **no** higher
  director scale to find without severing the DE anchor. Pursuing a second, non-DE
  orientation scale needs new physics (a second orientation sector); program
  **parked** pending that — `notes/nondeweld_director_scale_program.md`.

## 10. How to test in this framework — recurring traps (pick the right instrument)
General verification discipline is in CLAUDE.md §§1,3,4,6. These are the
framework-SPECIFIC ways a test gives a wrong answer through the wrong instrument;
each cost a real detour this session.
- **Never isotropically average an anisotropic (semi-Dirac / director) medium.**
  The physics lives on the resolved axes (easy-radial vs hard-perpendicular; the
  `sin⁴θ` structure; directional pressure-decoupling). Averaging over θ washes out
  exactly what carries the result. *This session:* Gate 1's angle-average gave
  "radiation-like"; the anisotropic reconsideration (route 1′) was the correct
  frame. Resolve the axes first.
- **A linear-dispersion gap tests FLUCTUATIONS, not COHERENT configurations.**
  `E(k)` describes small modes (phonons/Goldstone); it does NOT bound the
  gravitating energy of a large-amplitude coherent texture/defect (its energy is
  Frank/gradient energy, non-perturbative). *This session:* the route-1′ "mass
  gate" read the gapless band and under-answered. For "does it gravitate as
  matter," ask about **coherent configurations** (solitons/defects/windings),
  never a band gap. (See §2: mass is emergent — test the configuration, not a rest
  mass; for matter-scaling test temperature-independence + number-conservation.)
- **Don't assume gravity is LOCAL modified-Poisson without checking Paper VII's
  actual source structure.** Density feedback / geometric-cost / shared-n̂ (§7) may
  carry a nonlocal source term; a local `∇²Φ=4πG_eff δρ_local` calc silently
  excludes it. *This session:* `phaseA` assumed local; the Bell reframe implied
  nonlocal. Establish local-vs-nonlocal before trusting a local perturbation calc.
- **Flag cutoff-SET magnitudes as not-derived.** Several Paper VII results are
  explicitly cutoff-dependent ("not fixed here"): the acceleration scale, whether
  the elastic strain dynamics act on galactic timescales, the strain magnitude.
  Signs/scalings are robust; specific numbers are not. Don't quote a cutoff-set
  value as a derived prediction.
- **Run the Jordan/Einstein discriminator on any φ-power before trusting it**
  (§6). A scalar-tensor-norm power (`φ²⁰`, `φ⁻⁸ᐟ³`, `φ⁻⁸`) is likely a frame
  artifact; a topological/WZW power (`φ⁶`, `φ¹²`, `360/φ²`) is physical.
