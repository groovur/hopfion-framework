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
