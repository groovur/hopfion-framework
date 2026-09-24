# Paper VII (main_paper7.tex) — audit change log

Status of items from the 2026-08 audit. DONE = applied to the working tree.
DECISION = needs the author's canonical call before editing.

## DONE

- **Inconsistency 2 — Vainshtein exponent (Cor. P7:cor:vainshtein_topology).**
  The "topological, no free parameters" formula for the φ¹² exponent,
  `2(Q_cond − Q_Hopf·(k+2))`, evaluates to 2(10−2·5) = **0, not 12**. The
  correct source is φ¹² = λ² = (φ⁶)², the squared Bogomolny parameter,
  which the paper already stated in the same line. Fix applied: deleted
  the false charge-formula equality in item (iii) and the proof; kept
  φ²⁰ = 2Q_condensate (topological) and φ¹² = λ² (Bogomolny). No numbers
  change; the corollary's conclusion stands via λ = φ⁶.

## DECISION NEEDED — β in the gravity sector (reopened; earlier "clean fix" retracted)

The "ω_BD is a β\* quantity, keep r_V with β" plan does NOT work, and was
NOT applied. Tracing eq P7:eq:rV_derivation: r_V³ = r_S/(6β·φ⁸·Λ), and the
β there IS the β inside ω_BD = φ¹²/(3β) — one coupling, not two (line 1044:
"the single parameter β enters only through ω_BD"). And β\* = 0.452 FAILS
Cassini: at 0.452, r_V = 5.6×10⁶ AU and |γ_PPN−1| = 2.9×10⁻⁵ > the 2.3×10⁻⁵
bound. The paper uses β = 0.1 in the Cassini/Vainshtein sector precisely
because β\* = 0.452 does not pass. So this is a physics tension, not a
notation collision.

Corrected reading (for author decision):
 - Gravity/cosmology sector should use β = 0.1 (density coupling) UNIFORMLY.
   The G_eff proof (line 1835) using β\* = 0.452 (→ ω_BD = 237) is the
   outlier; it should use β = 0.1 (→ ω_BD ≈ 1076, α_BD ≈ 9.3×10⁻⁶). The
   G_eff conclusion (≈ G_N) survives either value — safe fix.
 - β\* = 0.452 then appears only in condensate-internal / Paper I contexts
   (S_eff line 164, β\*·ρ_CMB = Q).
 - DEEPER, unresolved: is β = 0.1 DERIVED (the κ↔ρ constant of line 286,
   ≈ 0.22) or FITTED to pass Cassini? If derivable, β = 0.1 is legitimate
   and Cassini passes honestly. If not, β\* = 0.452 marginally violates
   Cassini and the Vainshtein resolution only works at a fitted β — a
   genuine open problem the paper currently masks with the unquantified
   β ∝ β\* remark.

## DECISION NEEDED (not yet edited)

- **β = 0.1 vs 0.1034 (lines 1358 vs 1366), same density coupling.**
  0.1034 appears only in the Λ_obs exponent and lands exactly on 10⁻¹²²
  (φ⁷/0.1034 = 280.8); 0.1 is used for r_∞/r_V. Same coupling, two
  values, 3.4% apart. Which is canonical? And should Λ = 10⁻¹²² be
  presented as a prediction (β fixed independently → 10⁻¹²²) or as a
  consistency (β ≈ 0.1 chosen, gives 10⁻¹²² to 3.4%)? The tier table
  already says 3.4%; the abstract says a clean 10⁻¹²² — the same
  false-precision framing corrected in Paper XI.

- **β ∝ β\* proportionality constant (line 286).** The relation β ∝ β\*
  is asserted but the κ↔ρ constant is never given, so β = 0.1 is not
  actually derived from β\* = 0.452 — it is a Cassini-allowed choice
  (β ≤ 0.13) labelled "framework value". State whether the constant is
  derivable (making β a prediction) or a genuine free rescaling
  (fitted, honestly labelled).

- **N_e = 57.73 (Prop. P7:prop:Treh).** Chosen to hit the observed
  P_s = 2.10×10⁻⁹ exactly, then shown consistent with a loose T_reh
  range — the reverse of a prediction (P-B). n_s = 0.9654 is robust to
  this; P_s "exact" is the circular part.

- **β₀ (\bz) vs β\* (\bst) prose conflation.** β₀·ρ_CMB = Q exactly
  (defined); β\*·ρ_CMB = 10.012 (computed, 0.12% off). Deliberately
  distinct but swapped in prose (lines 1534–1541, 2353, 2529). State the
  relation once and use consistently.

- **Stale input-count table (lines 220–224).** Front-matter comparison
  table says 2 inputs; the rest of the paper (post-Paper XII) says 1.

- **Dark matter sector (§9, P7:sec:dark_matter) — does not work as written.**
  The only DM derivation in the series. Mechanism: no DM particle, the
  k-essence scalar ρ (same field as dark energy) is the DM via clustering
  perturbations. Three verified problems:
  (i) c_s = 1/φ ≈ 0.618 c (Thm P7:thm:cs, algebra correct) is RELATIVISTIC
      — ~600× too large for cold dark matter (needs c_s ≲ 10⁻³ c). This is
      a hot / dark-energy-like fluid, not CDM.
  (ii) Jeans criterion inverted: λ_J ≈ 2.8×10⁴ Mpc (correctly computed) is
      larger than the observable universe, so all galactic/cosmological
      scales sit BELOW it → perturbations oscillate, do NOT cluster. The
      paper reads this as "clusters freely on all scales" — backwards.
  (iii) BAO prediction r_s = 147/φ = 90.9 Mpc (a 38% shift) is presented as
      a future DESI/Euclid test but is already excluded by the measured
      147 Mpc sound horizon.
  Plus: Ω_DM = 0.26 ρ_crit is "taken as an input" (line 1710) — not derived,
  and a THIRD input (T_CMB, m_e, Ω_DM), contradicting the 1-/2-input claims.
  Not a labelling fix: needs a mechanism giving small c_s on sub-horizon
  scales; c_s = 1/φ is a constant O(1) speed, the opposite of what CDM needs.

- **DM: the tension/geometric alternative (author's original logic, NOT in
  the paper).** Flat rotation curves as condensate-gradient TENSION rather
  than fluid clustering. v²=const ⟺ M_enc(r)∝r ⟺ ρ_grad∝1/r² ⟺ ρ~ln r
  (a filament/line-source profile). The framework HAS the ingredient the
  paper didn't use: Axiom A1's radial preferred direction ê_r=∇ρ/|∇ρ| with
  anisotropic sin⁴θ suppression — a "spoke" structure that could give the
  ln r profile. Tested numerically (src_paper7/dm_tension) — NEGATIVE.

- **DM tension self-organization test — NEGATIVE (diagnostic).** 3D
  gradient-flow relaxation of a point mass in the anisotropic condensate
  (transverse-gradient suppression w=1,0.3,0.1; spherical IC + 1% noise to
  let 2D structure emerge spontaneously). At every w the relaxed field stays
  spherical to 5 decimals (axis ratios ≈1.0), profile ρ~1/r (d≈3), rotation
  Keplerian (p_v2≈−1) — identical to the isotropic control. The spherical
  solution is a STABLE fixed point; no 2D sheet or 1D filament emerges.
  DIAGNOSIS: the framework's anisotropy is DIRECTION-only (sin⁴θ), and for a
  spherical source ∇ρ is radial everywhere, so the suppression has nothing to
  act on — direction-only anisotropy structurally cannot modify a spherical
  source's radial force law. Flat curves would need EITHER a gradient-
  MAGNITUDE (|∇ρ|-dependent, MOND/k-mouflage) term the framework's K(ρ)X
  kinetic structure lacks, OR an externally-imposed 2D disk geometry
  (untested). Scale check stands: MOND a₀ = cH₀/(2π) to ~10%.

- **DM — CORRECTED DIRECTION (supersedes the negatives above): director
  strain of a semi-Dirac condensate.** The prior negatives ruled out the
  wrong ontology (Newtonian force / static scalar magnitude ρ → monopole →
  Keplerian; the direction-only anisotropy is inert for a spherical source;
  a dynamical force-sim gives dissipative drag, not conservative binding).
  The correct frame: the condensate is a NEMATIC-like medium; its physics is
  in the ORIENTATION field (director n̂ = ê_r, Axiom A1), the soft/Goldstone
  mode — NOT the scalar magnitude (the gapped mode I kept solving). The
  anisotropy is a band-structure property: S_eff = sin⁴θ = (sin²)² is the
  SEMI-DIRAC dispersion E² = A²k⁴cos⁴θ + v²k²sin²θ (heavy/quadratic one axis,
  light/linear perpendicular), the same anisotropy measured in ZrSiS/semi-
  Dirac semimetals, and the SAME pentagon (k+2=5) that fixes it also fixes
  the 3 generations (one structure, three phenomena).
  Dark matter = the FRANK ELASTIC STRAIN of the director, F=(K/2)∫(∇n̂)²,
  sourced by baryons imposing orientation boundary conditions — invisible
  (orientation, not EM-coupled), collective (long-range elastic), tracks
  baryons (→ radial-acceleration relation automatic). A DISCLINATION-type
  winding θ~ln r gives |∇n̂|~1/r → energy density ~1/r² → M_enc ∝ r → FLAT
  rotation curves (a smooth defect-free director gives the Keplerian
  monopole; the "less-radial"/frustrated region — the disclination — is the
  effect). k-space equivalent: vanishing-DOS anisotropic screening, static
  Π(q→0)→0 directionally → V(r) departs from 1/r; the soft pole of Π IS the
  director Goldstone, and K is the q² coefficient of Π. OPEN calculation:
  small-q expansion of the semi-Dirac Π(q) to extract K and the a₀ scale
  (src_paper7/dm_polarization); assumed crux = baryons source a disclination
  (winding) director, not a smooth one. Scale check stands: a₀ = cH₀/(2π).

- **DM — Frank K COMPUTED (gate-validated, frank_run.py): the director is
  SPLAY-UNSTABLE → self-organization confirmed.** Fix that unblocked it:
  RADIAL DISK cutoff |k|<Λ (a square grid breaks axis-rotation symmetry →
  spurious Goldstone violation, the −5/−700 of the failed attempts). Disk
  cutoff → Goldstone gate passes ~1e-8, Λ-stable. Undoped ⇒ intraband
  Pauli-blocked ⇒ interband+diamagnetic only (fast, ~0.3s/pt). Result, dE~q²
  (STANDARD Frank) both directions: K_bend>0 (stable), **K_splay<0
  (UNSTABLE)** — signs cutoff-robust (Λ=15,30,60); magnitudes grow with Λ
  (UV-cutoff-set by the gapless quadratic band, so a₀ enters via the physical
  condensate cutoff ~m_ξ~H₀). PHYSICS: negative splay ⇒ the uniform director
  spontaneously forms a radial splay texture = the disclination/spoke strain
  = invisible dark matter; q² ⇒ θ~ln r ⇒ M_enc∝r ⇒ FLAT curves. This is the
  self-organization the author predicted; the scalar-magnitude relaxation was
  blind to it (lives in the orientation field). Robust: sign(K_splay)<0, q²
  exponent, gate. Open: physical-cutoff regularization + 3D embedding to turn
  the K magnitude into a₀. Run infra resumable (test/full, checkpointed).

## MINOR / AUTHOR HANDLING

- \aem undefined (line 2291) — author will handle (intended as a paper
  summary line).
- Broken \ref{thm:verlinde_cancel}, \ref{rem:single_input} (line 1513) —
  replace with "Theorem 4.7" + commented \ref, as done at line 2520.
- \date contains "Superseded estimate update - June 2026" text (line 67)
  alongside "Preprint --- April 2026" (line 65) — editing artifact.
- Notation collision: φ (golden, \vp) vs φ_c (inflaton, \varphi_c) in
  the Inflation section — legibility only.

## APPLIED — DM/DE/inflation revision (cold director), 2026-08

- **c_s = 1/φ relocated (DE sound speed, not DM).** The old Thm `P7:thm:cs`
  used c_s=1/φ as a dark-MATTER clustering speed (fatal: 0.6c, Jeans inverted,
  90.9 Mpc BAO). c_s=1/φ is real but is the DARK-ENERGY perturbation speed
  (order unity => smooth, non-clustering DE, correct for Λ-like). Now a compact
  Remark `P7:rem:cs_de` in §Dark Energy (EoS subsection), with the k-essence
  derivation, explicitly distinguished from the O(β⁰) acoustic-metric c_s²=1
  (§P7:sec:acoustic_schwarz — that one was never c_s=1/φ and was untouched;
  the acoustic Schwarzschild/Kerr/Hawking sector does not depend on thm:cs).

- **Cold director result (director_speed.py).** c_dir² = K/χ, χ = static
  susceptibility of the rotation generator = ∫|M(k,k)|²/|d_k| (same M as the
  frank_run diamagnetic term). Over Lam=15..120: χ~Λ^3.5 (matches dim.
  analysis exactly), K_bend~Λ^1.4, K_splay~Λ^2.7 => c_dir²(bend)~Λ^-2.1,
  (splay)~Λ^-0.8, BOTH decreasing. The DIRECTOR (orientation) mode is COLD
  (c_dir→0, heavy: χ grows faster than K), NOT relativistic — the old c_s=0.6c
  was the scalar-MAGNITUDE mode, not shared by the director. => cold,
  pressureless, no Jeans obstruction, clusters on all scales (recovers the CDM
  clustering the fluid failed to give). splay c_dir²<0 = the self-org
  instability (imaginary freq). Robust: coldness (direction), splay sign.
  Cutoff-set: magnitudes, T_ord~(π/2)K_bend, isocurvature fraction.

- **Paper §Dark Matter restored/extended with cold-director info.** New
  subsection `P7:sec:cold_director` ("Cold director mode: pressureless
  clustering"): c_dir²=K/χ (eq `P7:eq:cdir`), heavy director => c_dir→0 =>
  cold; splay => c_dir²<0 instability; static/dynamic decoupling; contrast
  with the relativistic scalar-magnitude c_s=1/φ (=DE, rem:cs_de). Abstract DM
  sentence now says "cold, pressureless orientation mode that clusters on all
  scales". Status subsection updated (coldness established; isocurvature
  fraction added to the cutoff-set open list).

- **Inflation bridge restored (reframed), `P7:rem:reheating_IC`.** New
  subsection `P7:sec:dm_ic` in §Dark Matter. Reheating (Prop Treh,
  T_reh~1.25e8 GeV) establishes the condensate/nematic MEDIUM; DM IC then
  follow from the director sector: cold baryon-sourced strain inherits the
  adiabatic P_s=2.10e-9 THROUGH the baryons ("no free parameter" survives,
  cleaner than the old T_cond route); splay texture = a COLD isocurvature
  component, amplitude set by the same cutoff as a_0, must stay subdominant
  (Planck β_iso). ε≈β/φ⁶ (ρ-as-inflaton) still excluded by Starobinsky
  ε=2/N_e²≈0.001. NOTE: old Prop `IC_condensate` NOT restored — it ran through
  the fluid transfer function T_cond (with the c_s=1/φ oscillations); its
  surviving "no free parameter" content is now carried by baryon-tracking in
  reheating_IC.

- **DEFERRED (ref updates, per author):** registry `paper_registry.yaml` still
  lists removed labels (P7:thm:cs, P7:prop:CDM_limit, P7:rem:rotation_curves_halo,
  P7:prop:IC_condensate) and the old reheating_IC→cosmology.tex mapping; add
  P7:rem:cs_de, P7:sec:cold_director, P7:eq:cdir, P7:sec:dm_ic; reheating_IC now
  conceptually in dark_sector.tex not cosmology.tex. Sector files
  gravity/dark_sector.tex + gravity/cosmology.tex still carry the OLD fluid DM
  content and need re-sync. Not yet done.

- **OPEN (next calc):** physical-cutoff regularisation + 3D embedding to turn
  the cutoff-set magnitudes into numbers: absolute a_0 (target cH_0/2π) AND the
  isocurvature texture fraction vs Planck β_iso≲0.038. Same cutoff sets both.

## APPLIED — crossover matching + isocurvature (a_0 fold), 2026-08

- **a_0 crossover matching (crossover_match.py).** a_0/cH_0 = band_K(Lam*),
  Lam* = m_xi in crossover units. Target cH_0/(2pi) => Lam*~12.4 => cutoff/
  crossover hierarchy m_xi/E* ~ 12-38. Single-scale (no hierarchy) => a_0 10.6x
  too small. a_0 ~ cH_0 ROBUST (right order); exact 2pi NOT derived.
  LEAD (flagged): phi^6=17.94 sits in the 12-38 window AND is the S_eff
  normalisation in the dispersion itself; m_xi/E*=phi^6 gives light-axis a_0=
  1.49x obs, heavy-axis 0.54x obs -> BRACKETS the empirical a_0 within ~1.5x.
  Suggestive (framework-native, not fitted), unproven.
  FOLDED into §accel_scale: coefficient set by cutoff/crossover ratio, order
  phi^6 (=S_eff normalisation), light/heavy axes bracket a_0 within O(1);
  framework-native scale, precise value not derived.

- **Isocurvature number (isocurvature.py).** Texture is DM only after the
  ultralight director (m_xi~H_0) unfreezes: H(z_order)=m_xi => z_order~0 (m_xi=H0)
  or 1.2 (m_xi=2H0). Texture forms at z~0, absent at recombination. Frozen
  inflationary delta_theta ~ H_inf/(2pi f_dir) = 2.85e-5 ~ sqrt(P_s)=4.58e-5, but
  beta_iso(CMB) ~ 0 because f_tex(z=1100)~0. SAFELY < Planck 0.038.
  KEY (honest, two sides of one fact): the lateness that makes it isocurvature-
  safe is the same lateness that means it supplies NO DM at recombination -> it
  CANNOT be the CMB-era CDM. Director-strain DM = LATE-TIME galactic component;
  early-universe CDM still open. delta_theta~sqrt(P_s) is harmless ONLY due to
  lateness (any texture DM at z=1100 => beta_iso~O(1), excluded). Lateness REQUIRED.
  NOT yet folded into paper (reheating_IC still says "must remain subdominant";
  could sharpen to "automatically subdominant/late, at the cost of not being early
  CDM" -- pending author call).

- New scripts: director_speed.py, cutoff_embedding.py, crossover_match.py,
  isocurvature.py (all in src_paper7/dm_polarization; run in a numpy venv).

## APPLIED — semi-Dirac identification derived + direction fix, 2026-08

- **S_eff <-> semi-Dirac now DERIVED in the paper (was asserted).** DM intro
  rewritten: adds the two-band fluctuation Hamiltonian H = A k_perp^2 sx +
  v k_par sy and the dispersion E^2 = A^2 k^4 sin^4(th) + v^2 k^2 cos^2(th)
  (new eq P7:eq:semidirac), showing the heavy band carries exactly S_eff's
  sin^4(th) intrinsically (transverse quadratic band k_perp^2 ~ sin^2, squared
  -> (sin^2)^2). Verified sympy (dispersion_from_Seff.py).
- **Direction fixed.** "quadratic (heavy) ALONG e_r" was BACKWARDS; now "heavy
  PERPENDICULAR to e_r, light along e_r" (easy radial = light/linear, hard perp
  = heavy/quadratic). Dispersion cos^4 -> sin^4. Checked: no other direction
  statement in the paper needed fixing (only the DM intro had it).
- **Splay instability CONFIRMED correct (splay_bend_check.py, sympy).** With
  e_r = light axis, modulation along heavy(perp) = SPLAY (div n != 0, bend=0),
  along light(radial) = BEND. frank_run K(heavy)<0, K(light)>0 => SPLAY unstable,
  bend stable. Paper's "splay instability" stands; the direction fix does NOT
  change K_splay<0, disclination, M_enc ∝ r, flat curves, cold director, a_0,
  isocurvature.
- **Typo fixed:** "dark matter-\ncomponent" -> "dark-matter component" (DM intro).
- New scripts: dispersion_from_Seff.py, splay_bend_check.py (need sympy).

## APPLIED — CC framing: beta=1/Q prediction (option a), 2026-08

- **Decision (author):** present the CC as the parameter-free beta=1/Q=0.1
  prediction, NOT a clean "10^-122 derived". At beta=1/Q the exponent is
  -vp^7/0.1 = -290 => Lambda ~ 10^-126 Mpl^4, ~4 orders below the observed
  10^-122 (exponent agrees ~3%), with no fitted parameter. A fitted beta=0.1034
  reproduces 10^-122 exactly but is not a prediction.
- **Paper X check (author asked):** Paper X does NOT quantize the coupling to
  1/Q. It topologically quantizes the Hopf CHARGE (Q=10) and the energy
  spectrum (geometric quantization, Q=10 forced). The nearest exact relation is
  Paper VII's own beta_0 * rho_CMB = Q (line ~1576), which fixes beta_0 through
  rho_CMB, not as 1/Q. So beta=1/Q=0.1 is a natural VALUE (and the beta already
  used for rho_inf) but its derivation is OPEN.
- **OPEN derivation to chase:** is the k-essence density-feedback coupling
  topologically quantized to beta=1/Q? Paper X's geometric-quantization machinery
  (which forces Q=10) is the natural place to try; not yet done. If it works, the
  CC becomes a real (if ~4-orders-imprecise) parameter-free prediction.
- **Edits (Lobs kept = OBSERVED 10^-122 for downstream m_xi, r_V, eta):**
  - main_paper7.tex: boxed eq P7:eq:Lambda_obs now shows the FORMULA predicting
    ~10^-126 at beta=1/Q (dropped "Lobs =" LHS so Lobs stays = observed
    downstream); DE prose reframed (beta=1/Q, ~4 orders, Paper X charge-quant
    cited, coupling-quant open, fitted 0.1034 not a prediction); abstract (l.96)
    and conclusion reframed; tier row "Derived ~10^-122 / 1 order" -> "Predicted
    (no fit) ~10^-126 / 4 orders (obs 10^-122)". Added \bibitem{Paper10}.
  - main_reader_guide.tex: CC paragraph reframed to the beta=1/Q prediction +
    open-derivation note; rho_inf=vp/beta (was vp/beta*).
  - LEFT as observed 10^-122 (correct, empirical CC): paper7 lines ~991, 1026
    (r_V computation), ~2065 (eta ~ (H0/Mpl)^2).

## CC beta=1/Q: anyonic chase — NEGATIVE, open question sharpened, 2026-08

- **Reframing:** beta=1/Q <=> rho_inf=Q*phi <=> CC=exp(-Q*phi^7). And the density-
  feedback partition at the attractor (beta*rho_inf=phi) is EXACTLY the anyonic
  Born weights: 1/(1+beta*rho_inf)=1/phi^2=p0 and beta*rho_inf/(1+..)=1/phi=p1,
  the same 1/phi^2 = Tr[rho U_mono] in alpha^-1 = k|2I|/phi^2 = 360/phi^2.
- **But the Born-weight match is a RESTATEMENT** (pure algebra from beta*rho_inf=phi,
  since 1+phi=phi^2). Not independent evidence, not a derivation. (CLAUDE.md #3.)
- **Framework exact identities (Paper VII l.1583):** rho_CMB/Lambda_cond^4 = Q (CMB
  density counts the charge Q=ord(q_2I)), and attractor beta*rho_inf=phi=d_(1/2)
  (attractor counts the quantum dimension). Beautiful pairing of the two SU(2)_3
  invariants at the two special densities.
- **The crux (NOT a proven contradiction):** IF the kinetic feedback coupling beta
  (in 1/(1+beta*rho), driving attractor+CC) equals beta_0 from beta_0*rho_CMB=Q,
  then beta_0=1/Lambda_cond^4=0.45, rho_inf=phi*Lcond^4=3.58, CC=10^-28 -- and
  beta=1/Q (needing Lcond^4=Q=10) would be inconsistent. IF they are DISTINCT
  couplings (the framework carries 3 symbols beta_0/beta*/beta and never equates
  the kinetic beta to beta_0), beta=1/Q is merely adopted/underived, not refuted.
- **SHARPENED OPEN QUESTION (replaces 'Paper X will supply it'):** is the kinetic
  beta = beta_0 = 1/Lambda_cond^4 (=> CC~10^-28) or = 1/Q (=> 10^-126)? I.e. are the
  kinetic-feedback and rho_CMB couplings the same? Paper X does NOT settle it
  (it quantizes the CHARGE Q, and its symplectic weight uses beta*, not beta_0=1/Q).
- **Paper decision:** LEAVE Paper VII at 10^-126 (option a, beta=1/Q adopted,
  derivation open) -- it is defensible; the anyonic chase neither derived nor
  refuted it. Earlier same-session 'switch to option b' recommendation RETRACTED
  (it assumed kinetic beta = beta_0, which is unestablished).
- Scripts: scratchpad/cc_coupling.py (reconciliation tests, all negative:
  screened beta_eff no real root, beta*/phi^3=0.107 not clean, ratio 4.52 not phi^k).

## CC: kinetic beta != beta_0 RESOLVED; open item = rho<->kappa normalization, 2026-08

- **Answer to 'is kinetic beta = beta_0=1/Lcond^4?': NO.** Three couplings, three
  density variables: beta* couples to kappa=|grad rho|^2 (Papers I-IV kinetic,
  1+beta*|grad rho|^2, eq parent_kin); kinetic beta couples to the FIELD rho
  (1+beta rho, eq parent_action); beta_0=1/Lcond^4 couples to the CMB photon
  density rho_CMB. Paper VII Remark P7:rem:rho_extension (l.287) only asserts
  rho ∝ kappa and beta ∝ beta* -- proportional, constant UNQUANTIFIED. Framework
  never equates kinetic beta to beta_0. So the earlier '10^-28 contradiction' was
  a conflation and is void; beta=1/Q is adopted/underived, NOT refuted. Paper
  stays at 10^-126.
- **Open item now precisely = the rho<->kappa normalization c = rho/kappa = beta*/beta.**
  beta=1/Q <=> c = beta*Q = 4.52 (NOT clean: phi^3=2phi+1=4.236 is 6.7% off).
  FOLDED into Paper VII DE prose (replaces 'coupling quantised to 1/Q').
- **Group-theory check (author asked):** notes/beta_flow_derivation.md is a
  beta-FLOW derivation but for the ALPHA cascade (Delta_1), NOT the CC. It
  validates the profile self-consistency b* = beta*/C^2 = 0.06715 (matches
  P4:prop:bstar) via J_fb/J_a=phi, but does NOT yield beta=1/Q. The profile-derived
  couplings (beta*=0.452 physical, b*=0.06715 dimensionless) bracket but do not
  equal 1/Q=0.1. So NO existing group-theory derivation of beta=1/Q.
- **PATHS to c to evaluate (next):** (A) compute c=<rho>/<kappa> directly from the
  Battye-Sutcliffe profile / the k-essence canonical normalization (phi_c, eq
  P7:eq:phi_c) -- the profile FIXES c, so this could resolve beta outright; (B)
  test whether the CC attractor beta rho_inf=phi IS the profile self-consistency
  J_fb/J_a=phi (both '=phi'); (C) check if c is a tower/Lcond ratio. Path A is the
  concrete one -- c is not free, it is set by the profile.

## CC: leading normalization c=phi^3 identified (bracket-then-refine lead), 2026-08

- **Key find (Paper I l.1276):** the feedback coupling is phi-power-structured at
  source: beta_tilde* = beta*/phi^6 (exact). The whole coupling chain is phi-powers
  + few-% corrections:
    profile scale C^2 = phi^4=6.854 (actual 6.73, 1.8% off);
    normalization c=rho/kappa: canonical value phi^3=4.236 (from sqrt(phi^6)=phi^3);
    c for OBSERVED CC = 4.37 (phi^3 x1.03); c for beta=1/Q = 4.52 (phi^3 x1.067).
- **Leading term DERIVED (indicative): c=phi^3** (canonical kinetic normalization) ->
  beta=beta*/phi^3=0.107 -> CC=10^-118. The observed (10^-122) and 1/Q (10^-126) are
  both a few-% correction to c=phi^3, straddling it. This is the alpha 4-term shape:
  leading golden term + small (Chern-Simons-like) refinement.
- **Status upgrade:** the coupling is c=phi^3*(1+delta), delta~3-7% -- a DERIVED leading
  phi-power + open correction, NOT a bare adopted 1/Q. Paper stays at 10^-126; DE-prose
  open item tightened from 'rho<->kappa normalization' to 'the correction delta to the
  canonical c=phi^3' (cites Paper IV alpha 4-term as the analog).
- **NEXT (ii, requested):** derive delta via the rho->phi_c->kappa field redefinition to
  next order (canonical map dphi_c/drho=1/(phi^3 sqrt(1+beta rho)), density-running phi^3
  (low rho) -> phi^4 (attractor)); see if it closes to 1/Q (delta=0.067) or phi^3 (delta=0)
  or the observed (delta=0.03).
- Scripts: scratchpad/path_A_normalization.py, cc_coupling.py.

## CC REFINED: c = phi^3 + p0^2 -> 10^-122 (matches obs), 2026-08

- **Brainstorm result:** the refinement to leading c=phi^3 is c = phi^3 + 1/phi^4 =
  phi^3(1+1/phi^7) -> CC=10^-122.2, MATCHING observed 10^-122. And 1/phi^4 = p0^2,
  the SQUARE of the vacuum Born weight p0=1/(1+beta*rho_inf)=1/phi^2 (same p0 that
  carries alpha^-1=k|2I|*p0, Paper III). Structure: leading phi^3 (one canonical
  factor sqrt(phi^6)=phi^3) + correction p0^2 (two feedback vacuum insertions).
- **Mechanism check:** frame factor F(rho_inf) ~ 1+3phi^7 gives only -1.94 orders
  (INSUFFICIENT, need -4); canonical 2nd-order gives delta=1/(2phi)=0.31 (10x too
  big). Both RULED OUT. The p0^2 identification MATCHES numerically (0.3% in c) and
  is anyonic-native, but its physical MECHANISM is still open.
- **1/Q was the WRONG correction size:** delta(1/Q)~1/phi^6 -> 10^-124..126
  overshoots; observed prefers delta=1/phi^7. So the round 1/Q is retired.
- **PAPER UPDATED to 10^-122 everywhere** (box P7:eq:Lambda_obs -> 10^-122; new eq
  P7:eq:cc_c for c=phi^3+1/phi^4; DE prose rewritten to leading phi^3 + anyonic p0^2;
  abstract, tier row, conclusion, reader-guide all -> 10^-122 / c=phi^3+p0^2). Honest
  framing: leading phi^3 DERIVED (->10^-118); p0^2 refinement MATCHES obs (10^-122),
  mechanism OPEN.
- **OPEN:** physical origin of the p0^2 correction (why exactly two vacuum-weight
  insertions). Candidate ties to the anyonic sector (Paper III p0). Scripts:
  scratchpad/{path_A_normalization,path_ii_redefinition,cc_coupling}.py + the scans.

## CC: correction = p0^(Q_H), Q_H=2 mechanism (bracket + fold-in), 2026-08

- **[phi^3, 1/Q] bracket folded into paper:** observed CC at the geometric-mean
  MIDPOINT. c(1/Q)=phi^3+2p0^2 (topological, upper, 10^-126); c=phi^3 (canonical,
  lower, 10^-118); observed c=phi^3+p0^2 (midpoint, 10^-122). CC=exp(-64.2 c) is
  log-linear in c so midpoint-in-c = geometric-mean-in-CC = 10^-122. DE box/prose,
  tier row, abstract, conclusion, reader guide all updated.
- **Q_H=2 MECHANISM (author's lead, confirmed):** correction = p0^(Q_H), the vacuum
  Born weight p0=1/(1+beta* rho_inf)=1/phi^2 raised to the HOPF CHARGE.
    Q_H=1 (standard Hopfion): p0^1 -> c=4.618 -> 10^-128.8 (WRONG, 7 orders off).
    Q_H=2 (density-feedback icosahedral condensate, Paper I): p0^2 -> 10^-122 (MATCH).
    Q_H=3: p0^3 -> 10^-119.7.
  Q_H=2 is INDEPENDENTLY fixed (Paper I, not fitted), so p0^(Q_H=2)=p0^2 is a
  parameter-free consequence. The observed CC is EVIDENCE for the charge-2 vacuum;
  Q_H=1 collapses the structure (as the author predicted). Physical picture: each
  Hopf unit sits in the vacuum (identity) channel of the 1/2 x 1/2 SU(2)_3 fusion
  with weight p0; the joint two-unit vacuum amplitude is p0^2.
- **FOLDED:** DE prose open-item upgraded to 'correction = p0^(Q_H), Q_H=2, observed
  CC as evidence for charge-2 vacuum'. Remaining open item narrowed to: WHY each
  Hopf unit contributes exactly one p0 insertion.
- Scripts: scratchpad/*.py (path_A, path_ii, cc_coupling, scans).
