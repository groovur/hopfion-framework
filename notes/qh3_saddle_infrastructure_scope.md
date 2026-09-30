# Scope: Hopf-density infrastructure for the jammed Q_H=3 saddle (2026-09-27)

The blocked Paper XX problems (O1 absolute anchor, O2 gen-1 Φ, O3 tangent functional) all need the
converged jammed/self-consistent Q_H=3 saddle geometry (n, tangent, energy). Blocked by the RESOLUTION
WALL. This scopes what a better Hopf-charge-density discretization would take, and whether it can reach
the saddle at tractable N.

## THE CORE NUMERICAL PROBLEM
Faddeev J4 = ∫F² d³x, F_ij = n·(∂_i n × ∂_j n) — QUARTIC in first derivatives. The trefoil tube core is
thin (~1/C*~0.40, C*=2.5062). Uniform grid: ~10–15 pts across the core radius needs h~0.033 → over the
±5.6 box, **N~320 (≈33M pts)**. Measured cost scaling (~200 ms/step at N=96=0.88M) → ~7 s/step at N=320
→ ~8 h per relaxation, and a self-consistency (β/C*) scan → days. AND J4 convergence is still marginal
there (the N-convergence test: only ~5.7 pts/tube at N=160, J4 still drifting).

## WHAT DOESN'T WORK (established/reasoned this session)
- **Solid-angle (Berg–Lüscher) J4**: TRIED — WORSE than FD (grows faster). The field is genuinely
  UNDER-RESOLVED (few pts across tube); no LOCAL estimator recovers continuum J4 from insufficient data.
- **High-order / spectral FD**: limited by the SHARP arctan tube profile (Gibbs at the core); helps the
  (kh)^p error at the margin but not the fundamental under-resolution. Not a silver bullet.
- **Brute-force N~320**: possible but expensive (~days for a scan) AND marginal convergence.

## THE REAL FIX = RESOLVE THE TUBE EFFICIENTLY (concentrate DOF where the field varies)
- **(3a) TUBE-ADAPTED curvilinear coordinates [recommended for the true saddle].** Represent the field
  in (arc-length s, cross-section ρ, ψ) following the trefoil (Bishop frame — repo HAS bishop_frame_v2.py).
  Cross-section resolved cheaply (ρ,ψ ~15×15) → arc s~570 → **~1.3×10⁵ DOF vs ~3×10⁷ (≈260× fewer)**.
  Transform K, J4 into the curvilinear metric (careful but standard). CHALLENGE = the **3 CROSSINGS**:
  the knot self-approaches, so the field there depends on TWO arc positions — single-tube local coords
  break down. Fix: HYBRID (tube coords away from crossings + small local Cartesian patches at the 3
  crossings), or a superposition/blend. This is the main build effort.
- **(3b) REDUCED collective-coordinate ansatz + semi-analytic energy [faster, riskier first attempt].**
  Parametrize the config (C*, R0, r0, framing, crossing structure), minimise energy + self-consistency
  (V=φ) over the FEW parameters; compute the cross-section profile energy finely in 1D + local crossing
  patches (sidesteps the 3D grid). Cheap. RISK: the true jammed saddle may not lie in the ansatz family
  (the construction-C relaxation escaped — the saddle deforms).

## C* CHECK — RESOLVED (2026-09-27): the jammed tube does NOT fatten [the escape is closed]
Before running a new construction-level C* scan (`src_paper20/cstar_selfconsistency.py`), checked the
EXISTING compute. The relaxed FEEDBACK saddle already exists: `src_paper20/data/fb_saddle_scan/`
(`qh3_feedback_saddle.py`, N=64, h=0.175, C*=2.5062) annealed the full E=K_fb·J4 over β∈{0,0.2,0.35,
0.452,0.6,0.9} and saved the converged fields (`n_best_beta*.npy`). Measured the effective tube radius
(where f=π/2, i.e. mean n_z crosses 0, binned by distance to the trefoil) directly from the relaxed
fields. CONTROL: the analytic construction on the SAME grid gives exactly 0.399=1/C* (binning calibrated).
RESULT: relaxed radius ≈ **0.36–0.37 at EVERY β** (0.90–0.93× the construction 0.399) — the jam **thins
the tube slightly, does NOT fatten it**. The hypothesis (feedback softens stiffness → spreads director →
fatter resolvable tube) is FALSE. Feedback drives V→isotropy (V: 2.70→~1.1 as β:0→0.452, see scan_summary)
but does not spread the poloidal profile; spreading over MORE points would be easier to resolve, so nothing
prevented it — the energy just doesn't want it. **So the jammed tube is thin (~C*≈2.8), the "resolvable at
N~150, build unnecessary" branch is CLOSED, and the tube-adapted build (3a) is genuinely required.**
CAVEAT (honest): relaxed fields are at ~2.1 pts/radius (under-resolved), but the control is clean and the
TREND is thinning not fattening — no fat-tube preference to find at higher N. `cstar_selfconsistency.py`
(construction-level C* scan) is therefore redundant and its thin end unreliable — NOT RUN.

## a&b (lepton/quark split) CHECK — shape hypothesis REFUTED (2026-09-27, `tube_q2_compare.py`)
Tested the user's hypothesis "density feedback stable at lepton (Q=2) not quark (Q=3); the sharp quark
pushes the feedback out." Compared Q=3 trefoil (blend+partition) vs Q=2 unknot torus (winding q=2, a
SMOOTH crossing-free proxy — NOT Paper I's true two-tube inner+outer lepton) at matched tube sharpness C*.
Key β-independent diagnostic: V=J2iso/J2a is scale-invariant and ≤ its β=0 value V_max=J2iso_bare/J2a
(feedback only lowers V), so V=φ is reachable IFF V_max≥φ. RESULT (C*=2.5062): Q3 V_max=1.907, Q2 1.877;
saturation ~95% BOTH; self-consistent β (unique β making V=φ AND sopt=1) β_sc≈0.24(Q3) vs 0.22(Q2). So:
- **Saturation is TOPOLOGY-BLIND** — knot and unknot saturate feedback identically; saturation is a
  tube-SHARPNESS (C*) effect, not a crossings/knotting effect. "Quark pushes out feedback more" = FALSE.
- **BOTH have a self-consistent config** (V_max>φ, similar β_sc). The isolated quark is NOT "no config".
- ⇒ against the UNKNOT PROXY, neither first-order diagnostic separates lepton from quark.
- **BUT the proxy was misleading — the TRUE two-tube Q=2 REVIVES the split (2026-09-27).** Compared
  against Paper I's actual converged axisymmetric lepton (`src_paper1/f_fb_beta0.45200.npy`, reproduced
  its self-consistency: V=1.6179=φ, sopt=1.12, J4/J2a=0.226 vs target 0.227). Each sector at ITS OWN
  self-consistent β (the physical coupling where V=φ, sopt=1):
  - **True two-tube lepton: V_max=2.215, only 12.9% saturated → feedback RESPONSIVE/active.**
  - **Trefoil quark: V_max=1.907, ~92.7% saturated (β_sc≈0.24) → feedback SATURATED/pinned.**
  - Unknot proxy (V_max=1.877, ~95% sat) sat with the QUARK because it too is a single thin tube — that
    is why the proxy comparison wrongly read "topology-blind." The two-tube (inner+outer) structure
    spreads the lepton field → low g2 → feedback stays responsive; the sharp single quark tube (C*=2.5,
    radius 0.4) → high g2 → feedback pinned at 1/β. **(a) VALIDATED**: feedback is an active stabilizer
    for the spread lepton, self-saturated (ineffective) for the sharp quark = the user's "exploding quark
    pushes out the feedback." CAVEAT (§4): lepton on Paper I's 2D axisym engine, quark on the 3D tube
    engine → absolute saturation % not perfectly cross-comparable, but the 13%-vs-93% chasm and the
    engine-robust V_max (2.215 vs 1.907) both confirm it; saturation is sharpness-driven and the quark's
    sharpness is its own saddle value. [§9 note: hypothesis went right→"refuted"(proxy)→revived(true
    geometry); the refutation was proxy-specific. Building the real geometry, per user, was decisive.]
- **LABELING FIX + (b) recalibration (2026-09-27, user):** Paper I EXPLICITLY establishes single-tube
  instability — `P1:thm:vacuum`(iii), l.1318–1325: "the Q=1 sector is … not a stable ground state … but
  decays to Q=2" — via a CONFORMAL argument (the j₀ fusion channel has lower conformal weight). So the
  quark trefoil (a single KNOTTED tube) inheriting instability → confinement is largely PRE-ESTABLISHED;
  (b)'s saddle-search CONFIRMS + QUANTIFIES for the trefoil specifically, it does not discover. (a)'s
  feedback-saturation is a COMPLEMENTARY (energetic) mechanism, distinct from Paper I's conformal one —
  not redundant, not conflicting. CAUTION: Paper I's footnote l.323–326 says the axisymmetric solver's
  "Q=2" profile is itself a single tube (Battye–Sutcliffe Q=1) RELABELED Q=2 by the vol=2π·2r h² measure;
  the two-tube fusion is the THEOREM's picture. Tube-count label is OPEN (claude-hopfion §8). So the
  Paper XX remark's "two-tube torus" was FIXED → "extended, low-curvature torus" (label-neutral). The (a)
  comparison is label-INDEPENDENT (compared Paper I's actual converged STABLE lepton saddle, 13% sat, vs
  the trefoil, 93% sat) and stands. The lepton/quark distinction is really SPREAD (responsive feedback)
  vs SHARP-KNOT (saturated), not tube-count.
- **(b) STILL OPEN as a numerical confirmation** (is the quark's self-consistent saddle actually UNSTABLE?): needs the |∇E|² saddle
  search + Hessian. But (a) now supplies the mechanism for why it would not be stabilized (saturated
  feedback exerts no gradient-responsive restoring force). The "no self-consistent isolated quark" framing
  (from over-reading the β=0 descent) stays retracted — a self-consistent saddle DOES exist (V_max>φ); the
  question is only its stability.
- **(b) ATTEMPTED, INCONCLUSIVE — the |∇E|² saddle-search is ill-conditioned (2026-09-28).** Built
  `src_paper20/tube_saddle_search.py` (minimise G=½|P∇E|² via double-backward autograd; then lowest tangent
  Hessian eigenvalue by shifted power iteration). Runs at β=0.24 (Ns,Nr,Np=300,30,64, 2500 steps): |∇E|
  fell 41000→~8500 then PLATEAUED; a bug fix (sum G over FREE nodes only — pinned outer-shell/crossing
  nodes are BCs whose frozen nonzero ∇E floored G) lowered the plateau to ~5800 but it STILL stalls
  (decrements collapsing to ~0, J4 stays ~6100, V/sopt never leave ~0.30/~9.85 — never nears the
  self-consistent (φ,1)). So NO critical point is located, and the reported lam_min (+179 then +232, i.e.
  it even FLIPPED SIGN from the 60-step smoke's −1930) is evaluated far from a critical point and is NOT a
  valid stability statement. ROOT CAUSE (structural, not tuning): minimising |∇E|² squares the conditioning
  of E, and this landscape is brutally stiff (λ_top≈1.4e5) → gradient descent on G is hopelessly
  ill-conditioned. Cannot distinguish "no accessible critical point" from "method can't reach it" →
  INCONCLUSIVE for (b), NOT evidence of stability. VERDICT: rest (b) on Paper I's established single-tube
  instability (`P1:thm:vacuum` iii) + (a)'s saturation mechanism (both solid, (a) written into Paper XX
  `P20:rem:feedback_saturation`); the numerical confirmation was only ever quantification of an
  already-established fact and is not worth chasing with a stiff method. If ever revisited, use a
  FIXED-SCALE relaxation (descend E holding J4 or r_bar fixed → removes the marginal collapse mode,
  well-conditioned, read shape-stability directly), NOT |∇E|² minimisation.
CAVEAT: E=K·J4 is exactly SCALE-INVARIANT (K~λ, J4~1/λ), so descent collapses the size for ANY Q — the
scale is fixed by self-consistency, not minimization; "it collapses under descent" is therefore not the
quark-specific fact. NEXT to settle (b) + get O2 geometry: |∇E|² saddle-search + Hessian/stability (the
only tool that finds an UNSTABLE critical point and characterizes its stability); or build the true
two-tube Q=2 for an exact lepton contrast.

## BUILD PROGRESS — tube-adapted curvilinear solver (`src_paper20/tube_curvilinear.py`, 2026-09-27)
Phased, validation-gated (per user: build + validate locally before renting HW). KEY enabler: with a
BISHOP frame the tube metric is DIAGONAL — g=diag(h_s²,1,ρ²), √g=h_s·ρ, h_s=1−ρ(κ₁cosψ+κ₂sinψ), no
frame-twist cross term. So the Cartesian Faddeev energy carries over verbatim with scaled directional
derivatives D_s=(1/h_s)∂_s, D_ρ=∂_ρ, D_ψ=(1/ρ)∂_ψ and the h_s·ρ volume weight. Director n kept in LAB
frame so s4=(1−n_z²)² and the Hopf 2-form are unchanged.
- **PHASE 1 — straight tube (h_s=1): PASSED.** Curvilinear FD energy converges 2nd-order to an
  independent analytic-quadrature continuum truth, for K and J4, including s-winding (D_s operator).
  (J4 err quarters each doubling: 33→10→2.7→1.1%.) Gate caught a hand-truth error (spurious F_sψ; d_s n
  ∥ d_ψ n here → F_sψ=0) that the FD engine got right — §1 discipline confirmed.
- **PHASE 2 — circular tube (unknot, h_s=1+(ρ/R)sinψ): PASSED.** Validated against a SECOND independent
  engine (Cartesian FD port of qh3_feedback_saddle.E_fb) on the same analytic field. Both converge to
  K≈1314, J4≈2777; Cartesian approaches from below (dJ4% −12.9→−5.6→−2.8 as N:120→240). Validates the
  h_s CURVATURE factor. Demonstrated the payoff: curvilinear (320,80,160)≈4M DOF is BETTER converged
  than Cartesian N=240=13.8M pts. A trefoil eval ≈1M DOF (Ns~570×~30×64) → seconds; relaxation may be
  LOCAL-tractable, rented HW for final β/C* scans + resolution push.
- **PHASE 3 — trefoil + 3 crossings: PASSED, and the crossings turned out to be a NON-ISSUE**
  (`src_paper20/tube_trefoil.py`). Reused `bishop_frame_v2.py` (arc-length-resampled + Bishop curvatures
  κ1,κ2 from dT/ds). Crossing self-approach handled by "blend first" (user choice): the existing
  inverse-square two-strand spinor superposition, extended to a PARTITION OF UNITY
  W_s=ρ_s⁻²/Σρ_i⁻² so the self-overlapping chart is counted once (density is chart-invariant; only dV
  reweighted). VALIDATED against the independent Cartesian engine on the SAME blended field: converge
  toward one energy (finest dK=−1.5%, dJ4=−4.3%, Cartesian rising from below, same thin-tube signature
  as Phase 2). DECISION METRIC (`--diag`): overlap = 3.96% of nodes, but only **0.94% of K and 0.04%
  of J4** sit in crossings — the Faddeev term lives in the smooth single-strand tube wall, not the
  crossings. So **CARTESIAN PATCHES ARE NOT WORTH BUILDING** — the blend is sufficient. The feared hard
  part is closed.
- **PHASE 4 — relaxation/saddle-snapshot: BUILT + machinery VALIDATED, physics run OPEN**
  (`src_paper20/tube_relax.py`). torch/autograd relaxation of n(s,ρ,ψ) on the tube grid: geometry
  (frame, h_s, partition W, lab positions, construction IC) PRECOMPUTED once; outer-ρ shell + crossing
  (W<1) nodes PINNED to anchor topology; two-phase tangent-projected Adam + angle-clamp; V=φ/sopt=1
  best-score snapshot; plateau/ringing/collapse tail verdict. All mirrors qh3_feedback_saddle.
  - **VALIDATION (key):** at β=0 the tube grid gives V=1.91, sopt=6.14, MATCHING the WELL-RESOLVED
    Cartesian (Phase 3 N=200: sopt=(φ⁶·5851/2854)^½=6.06). The scan_summary N=64 values (V=2.70,
    sopt=4.03) were the UNDER-RESOLVED artifacts the memory flagged — the tube grid corrects them.
  - **DYNAMICS healthy:** E decreases monotonically, J4/J2a~5 (topology held), r_bar pinned at R0, no
    collapse/dilution; phase-2 transition fires.
  - **OPEN (the physics, needs real compute — do NOT pre-guess §10):** (1) on the RESOLVED grid,
    feedback pushes sopt AWAY from 1 (K_fb shrinks with β) and the construction starts far from target
    (sopt~9.9) — so self-consistency must be reached by RELAXING n, not by β alone; whether descent
    passes through a metastable min-score plateau BEFORE the saddle collapses (J4→0) is the empirical
    question. (2) Schedule tuning needed: the K_rise phase-2 trigger fired early (step 200) and froze
    progress at lr2=1e-5 — raise --K_rise_eps (stay in phase 1) and --n_steps. (3) β and the V/sopt
    TARGETS themselves may be knot-specific (memory's standing open calibration).
  - **COST (§8):** ~30 ms/step at 128k nodes → ~330 ms/step at the real 1.34M grid (420,40,80) → ~16
    min/β for 3000 steps + ~2 min precompute. A 4-β scan ~1.2 hr — LOCAL-tractable; rented HW only for
    resolution push / wider scans.
Modules uncommitted, self-tests inline: `tube_curvilinear.py [--phase2]`, `tube_trefoil.py [--diag]`,
`tube_relax.py --Ns .. --Nr .. --Np .. --n_steps .. --beta_list ..`.

## HIGHEST-ROI FIRST STEP (do BEFORE the big build) [SUPERSEDED by the C* check above — kept for context]
**Determine the self-consistent tube width (the jammed C*).** C* is set by the Derrick/Bogomolny (J2 vs
J4) balance, so the jammed saddle HAS a self-consistent value — which MAY differ from the isolated
crossing-analysis C*=2.5. If the jammed tube is FATTER (C*~1, radius~1), it is resolvable at moderate
N~150 with EXISTING solvers, and the whole expensive build is unnecessary. Cheap 1-parameter energy
minimisation over C* (in the jam, with feedback). Tempered: it may instead confirm C*~2.5 (thin, build
needed) — but either way it's the essential cheap first step that decides the whole effort.

## EFFORT / RISK / VERDICT
- C* check: ~1 day, HIGH ROI (might avoid everything or confirm the need).
- Tube-adapted build (3a): weeks (curvilinear Faddeev energy + crossing patches + relaxation +
  validation against a known Hopfion). Moderate risk (crossings). Best shot at the TRUE saddle.
- Reduced ansatz (3b): days–week. Higher risk (ansatz richness). Fast if it works.
**Could it reach the saddle at tractable N?** If jammed C* fatter → YES trivially. Else tube-adapted →
plausibly YES (~10⁵ DOF), but crossings are a real risk. Honest: a real multi-week effort, moderate
risk, NOT a guaranteed win. START with the C* check.
