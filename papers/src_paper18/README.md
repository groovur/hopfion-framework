# src_paper18 — scripts

Two independent lines of investigation live here: (1) exact Jones-polynomial
computation via a validated Kauffman-bracket/Temperley-Lieb engine, and
(2) numerical field construction / gradient-flow relaxation of the
condensate itself. Requires `python3.11`. The field-construction scripts
need `torch` and `scipy`; the Jones scripts need `sympy`.

## Jones polynomial / knot invariant scripts

```bash
# Engine + validation: checks Theorem 2.3's Jones_{T(2,3)}(q5)=-1 claim
# against direct computation, and computes the correct satellite invariant
# (companion Hopf-link component included, via cabling).
python3.11 kauffman_bracket_engine.py      # reusable TL engine, no output
python3.11 jones_T23_satellite_check.py    # the core check + satellite value

# Extends to the whole twist-count family n=0..15 (n=3 trefoil, n=5
# cinquefoil, even n = unfused 3-component links). Exact symbolic values,
# period-10 / mirror-symmetry structural checks.
python3.11 jones_satellite_twist_family.py

# Analytic continuation of the twist count to continuous n, so the path
# between n=3 (trefoil) and n=5 (cinquefoil) can be examined as a curve.
# Finds a non-monotonic path: peak near n~3.6, trough near n~4.3.
python3.11 jones_satellite_continuous.py

# Channel unit-factor verification (Theorems channel_1/channel_2 in
# main_paper19.tex): exact cyclotomic reduction, norm-in-Z[zeta5]
# verification via resultant with Phi_5.
python3.11 ../src_paper19/jones_unit_equivalence.py
```

All of these are read-only checks (print exact/numeric results, no files
written) except where noted.

## Field construction / gradient-flow scripts

```bash
# Bishop frame (parallel-transported, holonomy-compensated) for the
# T(2,q) curve. q=3 default (trefoil); generalised to any odd q so it
# also supports q=5 (cinquefoil) etc.
python3.11 bishop_frame_v2.py

# Builds Q_H=2/Q_H=3 initial conditions via phase-winding interpolation
# (Phi = chi + wind_exp*tau) and relaxes with gradient flow. This is the
# ORIGINAL winding-interpolation approach (fixed T(2,3) curve geometry,
# only the fiber phase winds continuously 2->3).
python3.11 crossing_transition_v2.py --mode compare --n_steps 2000

# Builds T(2,q) fields with the curve GEOMETRY itself continuously
# interpolated (q_pol), independently of the fiber phase winding
# (alpha_wind). Two constructions:
#   - build_torus_knot_field(): generic 2-arc nearest-strand split.
#     Valid for q_pol in [2,3] (where each half-arc covers ~1.5 poloidal
#     loops, a reasonable local nearest-strand proxy). NOT reliable much
#     beyond q_pol~3.5-4: each half-arc then covers >2 loops, so
#     "nearest point in this half" can land on the wrong loop, aliasing
#     the true symmetry (visually confirmed: q_pol=5 unrelaxed field
#     showed a spurious ~3-fold pattern instead of 5-fold).
#   - build_torus_knot_field_nlobe(): generalises build_trefoil()'s
#     3-lobe Construction C to any odd integer q_pol>=3, using q_pol
#     separate lobe-KDTrees so "nearest lobe" is always well-posed.
#     Written to fix the above, but NOT currently wired into the CLI
#     routing: benchmarked slower than the generic method at q=5 (16.4s
#     vs 4.3s build time), AND direct inspection of relaxed
#     (gradient-flow final) fields showed the generic method already
#     converges to the correct T(2,q) topology at every q in 3-5 despite
#     the unrelaxed-frame artifact -- so there was no correctness gain to
#     justify the slower build. Kept available (call directly) in case a
#     future check needs a construction whose INITIAL, unrelaxed frame is
#     also topologically clean.
python3.11 hopf_link_construction_v2.py --N 64 --q_pol 3.5 --alpha_wind 1.0 \
    --verify_only --outdir some_dir          # build + diagnostics only
python3.11 hopf_link_construction_v2.py --N 64 --q_pol 5.0 --alpha_wind 1.0 \
    --n_steps 2400 --whitehead_N 64 --outdir some_dir   # full gradient flow
```

`--verify_only` builds the field, runs the Whitehead $Q_H$ check once, and
exits (fast, ~10-30s). Without it, full gradient flow runs (minutes, since
each invocation also builds+relaxes the reference $T(2,3)$ trefoil for
comparison — this doubles the cost of every run but keeps a fixed,
known-good baseline alongside the target construction).

Caveat re: the raw Whitehead $Q_H$ diagnostic — it is unreliable on
*unrelaxed* initial fields (observed values scattered between 0.4-1.4
regardless of the target $Q_H$) and needs adequate `--whitehead_N`
resolution to be trusted even on relaxed fields (32 is fast but coarse;
64 costs only ~0.15s more per call and is cheap enough to always use).

### IMPORTANT: plain gradient flow loses topology over long runs

`hopf_link_construction_v2.py`'s `run()` uses plain, unconstrained Adam.
At long step counts (~2400+) this was found to let the field escape its
starting topological sector: per-grid-point rotations of O(10-20 rad) in
a single step are large enough to jump the topological wall on the
discretised grid. Symptom: energy collapses ~98% and a visible chunk of
the field's vector energy detaches and drifts away from the main tube
(confirmed in rendered output, in BOTH the target construction and the
dedicated reference trefoil — not q-dependent). **Do not trust long
(>~1000 step) runs of `hopf_link_construction_v2.py`'s plain flow for
stability questions** — use `gradient_flow_constrained_q.py` (below)
instead, which was built specifically to prevent this.

```bash
# Topology-protected (angle-clamped, two-phase) constrained gradient
# flow, generalising papers/src_paper16/gradient_flow_constrained.py
# (originally q=3-only) to any real q_pol via the general 2-arc
# construction above. Only builds/relaxes ONE field per run (no redundant
# reference-trefoil rebuild), roughly halving cost vs. the original.
#
# IMPORTANT: pass --delta_max_deg 2.0, NOT the script's own default of
# 9.0. Paper XVI (main_paper16.tex, lines 3240-3297) found the 9 deg
# clamp still lets floating-point moment-estimate drift in Adam
# accumulate into a spurious (2+1) Z3-symmetry-breaking pattern across
# the trefoil's 3 midpoints (visually: 2 of 3 midpoints densify, 1
# depletes) -- confirmed to be a NUMERICAL ARTEFACT, not real physics,
# because tightening to 2 deg suppresses it (fraying stays Z3-symmetric
# throughout Phase 2 at 2 deg). The CLI default is left at 9.0 only for
# consistency with gradient_flow_constrained.py's own default; always
# override it for any run whose result you intend to trust.
python3.11 gradient_flow_constrained_q.py --q_pol 4.0 --n_steps 20000 \
    --delta_max_deg 2.0 --outdir data/gf_constrained_q4.0
```

Two-phase schedule: Phase 1 (large `--lr1`, no angle clamp) does cheap
cleanup; Phase 2 (small `--lr2`, `--delta_max_deg` clamp per step)
carefully approaches the energy minimum once $K$ starts rising again,
without letting any grid point jump the topological wall in one step.
`--n_steps 20000` is a large budget, not an expectation — halt
conditions (dilution, near-vacuum growth, $J_4/K$ collapse) will stop a
run early if it genuinely destabilizes even under the constraint.

## Sweep orchestration scripts

```bash
# --verify_only energy sweep across q_pol=2..5 (11 points), no relaxation.
# Checks for a barrier/trough shape in the raw constructed-field energy
# between the trefoil (q=3) and cinquefoil (q=5) regions.
python3.11 sweep_beyond_trefoil.py
# -> data/sweep_beyond_trefoil/summary.json + n_initial_*.npy snapshots

# Full gradient-flow relaxation at 5 key q_pol points (3.0, 3.5, 4.0, 4.5,
# 5.0), testing whether any point shows a distinguishing stability
# signature (small energy/Q_H drift = near equilibrium) vs. the others.
# Motivating question: does the trefoil (Q_H=3) sit at a genuine
# potential-energy stop, or can the same formation energy push the
# condensate on toward Q_H=4 (acknowledged but unestablished in Paper
# XVIII, rem:higher_charge_incoming) and Q_H=5 (T(2,5), unexplored
# anywhere in the framework)?
python3.11 sweep_gradient_flow.py
# -> data/gradient_flow_sweep/<q>/... + summary.json
```

Edit `Q_VALUES`, `N_STEPS`, `WHITEHEAD_N` at the top of
`sweep_gradient_flow.py` to change the sweep. At `N=64` grid resolution,
budget roughly 0.5s/gradient-step *per construction*, and each point runs
two constructions (target + reference trefoil) — e.g. 2400 steps is
roughly 35-40 minutes per point, several hours for a 5-point sweep. Run
with `run_in_background` (or equivalent) rather than waiting inline.

### Results so far (see `data/` subfolders for full logs)

- `sweep_beyond_trefoil/`: raw constructed-field energy vs. $q_{\rm pol}$
  is non-monotonic between the trefoil and cinquefoil — a local peak near
  $q\approx3.5$, a trough near $q\approx4.0$ (close to the $Q_H=4$ sector),
  then rising again toward $q\approx4.75$-$5.0$. Unrelaxed-ansatz energy
  only; not yet confirmed to survive relaxation.
- `gradient_flow_sweep_800steps/`: an earlier, shorter (800-step) full
  relaxation sweep showed no distinguishing stability signature — every
  point dropped by a uniform ~68-71% in energy with no plateau anywhere,
  meaning 800 steps was not enough for any point to approach equilibrium
  (not a null result on the underlying physics question).
- `gradient_flow_sweep/`: current, longer (2400-step, `whitehead_N=64`)
  rerun of the same 5-point sweep, in progress / most recent.

## Reference images (project root, not this folder)

- `knots.png` — $T(2,1)$, $T(2,3)$, $T(2,5)$ standard torus-knot renders.
- `cinqfoil.png` — constructed-field renders at $q_{\rm pol}=4, 4.25, 4.5,
  5$ (the unfused/transitional region between trefoil and cinquefoil).
