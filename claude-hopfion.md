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

## 8. Sector ladder (charges) — see Papers XIII, XVII
- **[E]** Topological-charge ladder: `Q_H=0` vacuum/**photon** (ripple, Paper
  XVII), `1` neutrino, `2` lepton, `3` quark; gauge bosons/photon appear as
  `Q=4, J=1` composites of two `Q=2` (Paper XIII). Verify the exact assignment in
  the cited paper before use.

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
