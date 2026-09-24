# Non-DE-welded director scale program (SPECULATIVE — not in any paper)

Goal: find whether the director/orientation sector admits a SECOND mass scale,
NOT the ultralight m_xi~H0 welded to dark energy — one high enough that its
ordering (or near-threshold transient order) operates near RECOMBINATION rather
than z~few. If it exists, the recombination-CDM gap (option 3, `fabric_cmb_
hypothesis_and_plan.md`, base file §9) could close as a THIRD constant from the
one knot (Paper III: "two constants from one knot" -> alpha + linking scale).

All below is [SPECULATIVE]. Tag epistemic status; nothing enters papers unless derived.

## The two COUPLED questions
A. THE SCALE: is there a non-DE-welded director scale, and can it be >~ eV (so its
   ordering threshold sits near T_rec~0.26 eV, not Lam_cond~1.2e-4 eV = z~few)?
B. THE MECHANISM: near that threshold, does transient/oscillating angular order
   produce gravitating wells that SEED structure, without a permanent cold component?
They are coupled: the mechanism (B) needs the scale (A). At recombination we are
~2000x above the CURRENT threshold (T_rec/Lam_cond ~ 0.26/1.2e-4), where pretransitional
order is exponentially negligible. So B can only operate at recombination if A supplies
a higher director scale.

## The mechanism (user's hypothesis, captured faithfully)
alpha's POTENTIAL (the 2I angular structure / value 360/phi^2) exists always in the
vacuum geometry, but "heat" prevents the angular structure from HOLDING. Near threshold,
transient 2I/director order flickers into existence (fast, UV-timescale ~ "attoseconds,"
like an electron reconfiguring), momentarily creating a "vacuum structure" (a well) to be
filled, but the surrounding hot soup disperses it before it coheres globally (c_dir->0 =
slow global locking). Further cooling -> more/longer-lived order -> at T_ord, permanent.
=> the transient wells, if biased by density/temperature, could seed structure.
PHYSICS STATUS: this is standard PRETRANSITIONAL / short-range order above a nematic
transition (real: cybotactic clusters; critical slowing-down as T->T_ord). Sound as far
as it goes. The flicker being fast-local / slow-global reconciles "rapid" with c_dir->0.

## THE CRUX for B: COHERENT vs INCOHERENT (this decides peaks vs hump)
Transient order-parameter fluctuations are an ACTIVE/INCOHERENT source unless something
phase-locks them. Active/incoherent sources (like generic defect networks) give a broad
featureless CMB hump, NOT the sharp harmonic acoustic peaks -> EXCLUDED as the structure
seed. The ONLY way the mechanism survives the peaks is if the transient ordering is
GLOBALLY COHERENT (phase-locked across space). Candidate coherence source: the Bell /
shared-geometry mechanism (base file §7, Paper VIII) -- one connected 2I configuration
whose ordering is a global boundary condition, read out causally (c-limited). If the
flicker is synchronized by that global geometry (a shared "clock"/cooling), it could
imprint COHERENT (peaked) structure. This ties together the whole session's threads:
transient order (this note) + oscillation + Bell global correlation (§7).
=> DECISION GATE for B: is the transient director ordering coherent (global, phase-locked
   -> peaks) or incoherent (independent local fluctuations -> hump, excluded)?

## Part A: candidate non-DE-welded director scales (to investigate)
- **Paper III LINKING SCALE (the 2nd "constant from one knot").** The Hopf LINKING is
  itself an angular/orientation structure (linking of field lines). If the director's
  angular order is tied to the linking, the linking scale -- set by the WZW/2I structure,
  same as alpha, NOT by DE -- is the natural non-DE-welded director scale. FIRST TASK:
  read its VALUE in main_paper3.tex; is it >~ eV (threshold near T_rec)? This is the most
  promising candidate and directly realizes "third constant from the knot."
- m_WZW (P4:thm:compton, beta*=1/m_WZW^2): the WZW quantum mass scale of the Q_H=2 knot.
- R_0 / the Hopf-polar / linking-number scales (P3, P4).
- CHECK each: (i) non-DE-welded? (ii) magnitude vs T_rec? (iii) does it govern the
  ORIENTATION (director) sector specifically, or only the charge sector?

## Obstacles (honest, must each be cleared)
1. SCALE-MECHANISM COUPLING: recombination is ~2000x above the current (Lam_cond) threshold;
   need A to deliver a >~eV director scale or B is exponentially suppressed at z~1100.
2. COHERENCE (the crux): incoherent transient order is peak-excluded; needs global phase-lock
   (Bell/§7) to give sharp peaks. Unproven.
3. ADIABATIC vs ISOCURVATURE: the existing splay texture is isocurvature (base §9); a new
   seed must be adiabatic (or the peaks fail) -- check what the transient-order seed produces.
4. NO DOUBLE-COUNTING: a new director scale must not re-count energy already assigned to DE
   (m_xi~H0 / Lam_obs) or to radiation (Lam_cond bulk). Separate bookkeeping cleanly.
5. ABUNDANCE: must yield Omega~0.26 at recombination.
6. WHY TWO director scales? Physical justification for a 2nd scale decoupled from m_xi~H0
   (which is welded to DE by the a0<->cH0 unification, §4). What breaks the welding?

## Connections
- Base file §9 (two-cooling picture; the "third constant from the knot" open direction).
- §7 Bell mechanism (the coherence source for the crux).
- `fabric_cmb_hypothesis_and_plan.md` (why option 3 stands; what a rescue must beat).
- `delta1-dielectric-program` / `OP4_delta1_investigation.md` (the same one-knot scale-origin
  problem in the charge sector; shared actor = feedback coupling beta).

## FIRST CONCRETE TASK
Read the LINKING SCALE value + definition in main_paper3.tex (the 2nd constant). Determine:
is it non-DE-welded (set by WZW/2I, not H0/Lam_obs)? What is its energy value? Could its
ordering/near-threshold behavior sit near T_rec~0.26 eV? That single read gates whether Part A
has a live candidate before any mechanism (B) work.

### FIRST TASK — DONE (2026-08-31). Result: linking scale is DIMENSIONLESS -> eliminated as the mass.
main_paper3.tex: the "Linking Scale" is the Morse/virial parameter mu* = sqrt(k+2)/phi = 3-phi
~ 1.382 (P3:eq lines 255-264), and the Linking Scale Conjecture is the energy-integral RATIO
J_4/J_a = 2^{4/3}/phi^5 ~ 0.227 (P3, l.336). BOTH "constants from one knot" (alpha, linking scale)
are PURE NUMBERS. => the linking scale is the WRONG TYPE to be a dimensionful director mass (same
number-vs-mass trap as mu_UV/mu_IR). It could at most be a dimensionless FACTOR inside such a mass;
the dimensionful anchor must come from elsewhere. Candidate 1 ELIMINATED as the scale itself.

### REDIRECT — Part A now points to the DIMENSIONFUL scales, and one specific question
Dimensionful scales in the framework: m_WZW (P4:thm:compton, beta*=1/m_WZW^2 -- the Q_H=2 charge-
sector knot mass, ABOVE Lam_cond, plausibly near/above T_rec), R_0 (knot size), Lam_cond (sets the
CURRENT director stiffness K~Lam_cond^2 -> late ordering), Lam_UV, m_xi~H0 (DE-welded).
THE REAL CRUX (obstacle #6, sharpened): the director stiffness is currently Lam_cond-set (-> z~few).
Is there a physical route by which the DIRECTOR sector couples instead to a HIGHER, charge-sector
scale (m_WZW)? If the orientation stiffness were tied to m_WZW rather than Lam_cond, its threshold
could sit near/above T_rec -> ordered/cold at recombination. NO current reason it does -- but this is
now the precise, well-posed question, and it is the honest content of "why two director scales /
what breaks the DE-welding." NEXT: read what physically sets K (the Frank stiffness) in Paper VII/
the DM papers -- is it forced to be Lam_cond^2, or could an internal (m_WZW/R_0) scale enter?

### SECOND TASK — DONE (2026-08-31). Result: PART A BLOCKED at the framework level (negative).
Read what fixes the orientation scale in main_paper7.tex:
  1. m_xi^2 = phi^{2Q}·Lambda_obs·MPl^4 with Q=10 (condensate Hopf charge), l.1107/1496.
     => m_xi = phi^10·sqrt(Lambda_obs)·MPl ~ sqrt(3)·H0. The scalar mass is ANCHORED TO
     Lambda_obs (dark energy) AT THE FIELD-EQUATION LEVEL: m_xi ∝ sqrt(cosmological constant).
     This IS the welding, and it is structural.
  2. Frank stiffness K depends on m_xi (P7:rem:elastic_solid) -> K inherits the DE anchor.
  3. DECISIVE: the semi-Dirac dispersion's OWN UV cutoff -- the DM-sector condensate bandwidth --
     IS m_xi ~ H0 (l.1898-1900); the crossover is BELOW it (cutoff/crossover ~ phi^6 ->
     crossover ~ m_xi/phi^6). So the ENTIRE orientation/DM sector lives at energies <~ H0.
CONSEQUENCE: Part A has no candidate and the sector's construction FORBIDS one. The orientation
sector's ceiling is its UV cutoff m_xi~H0 ~ 1e-33 eV. A recombination-era director scale
(T_rec~0.26 eV) is ~30 orders of magnitude ABOVE the sector's own bandwidth. The higher framework
scales (tower M_n=Lam_cond·phi^{2n}; m_WZW) are all MAGNITUDE/CHARGE-sector, NOT orientation. To
raise the orientation cutoff above H0 you must sever m_xi from Lambda_obs -- the SAME anchor the
dark-energy derivation rests on (load-bearing, not a knob).
=> PART A BLOCKED. Part B (transient order) is moot without a scale to place it near recombination.
   The recombination-CDM gap (option 3) stands FIRMER: not "the director mass is ultralight" but
   "the entire orientation sector's UV cutoff is m_xi~H0 ∝ sqrt(Lambda_obs); no higher director
   scale exists without remaking the DE anchor." Only way forward = NEW PHYSICS: a second
   orientation sector with its own higher, non-Lambda_obs cutoff -- a major structural addition
   that ripples into the DE derivation, not a small fix. PROGRAM PARKED pending that.
