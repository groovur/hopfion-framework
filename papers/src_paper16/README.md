# Paper XVI — Numerical Dataset and Script Reference

**Paper:** *Quark Generation Masses in the Density-Feedback Hopfion:
A Survey of Ruled-Out Mechanisms and the $Q_H=3$ Gradient-Flow Landscape*
(Paper XVI in the density-feedback Faddeev–Niemi Hopfion series)

This repository contains the Python scripts, field snapshots (`.npy`), and
run logs (`.json`) underlying the numerical results reported in Paper XVI.
All scripts require Python ≥ 3.11, PyTorch ≥ 2.0, NumPy, and SciPy.

---

## Dependency graph

```
bishop_frame_v2.py          ← standalone, no project deps
    └── perstrand_s3lift_v2.py
    └── whitehead_perstrand_charge.py
    └── gradient_flow_constrained.py
    └── arc_segment_analytic.py
    └── junction_flow.py
```

`bishop_frame_v2.py` must be in the same directory as every other script
that imports it.

---

## Scripts

### `bishop_frame_v2.py` — Holonomy-compensated Bishop frame

Builds the arc-length-uniform, holonomy-compensated parallel-transport
frame along the trefoil T(2,3) at R₀=3, r₀=0.874. The frame has total
holonomy −63.019° (confirmed to < 3×10⁻⁶ degrees). Used as a dependency
by all other scripts.

**No CLI.** Import and call directly:

```python
from bishop_frame_v2 import build_compensated_frame_arclength
t, T, N1, N2, H = build_compensated_frame_arclength(NT=20000)
# H is the holonomy angle in radians
```

**Paper reference:** Sec. `sec:bishop_holonomy`, Proposition `prop:holonomy_value`.

---

### `perstrand_s3lift_v2.py` — Per-strand Construction C (diagnostic)

Standalone diagnostic that builds the per-strand S³-lift ansatz
(Construction C, `constr:perstrand`) using the corrected 3-lobe
strand search. Provides the `Z_perstrand_s3lift_v2(pts)` function
returning the Hopf doublet (z₁, z₂) at arbitrary query points. Also
confirms the crossing-midline field is smooth.

**To run the crossing diagnostic:**

```bash
python perstrand_s3lift_v2.py
```

Output confirms that the robust lobe-partition strand search eliminates
the KDTree silent-failure bug at crossing midlines, and that the field
passes through S³ continuously.

**Paper reference:** Construction `constr:perstrand`,
Propositions `prop:phase_jump`, `prop:s3lift_charge_zero`,
`prop:perstrand_smooth`.

---

### `whitehead_perstrand_charge.py` — Hopf charge verification

Computes the Hopf charge Q_H of the per-strand Construction C field via
the Whitehead (Berry–Chern–Simons) integral:

```
Q_H = (1 / 4π²) ∫ A·(∇×A) d³x
```

Includes a Q_H=1 reference field for validation. Convergence toward
Q_H=3 as N increases is the key output (Theorem `thm:perstrand_charge_three`).

**Usage:**

```bash
# Coarse (fast, ~2 min)
python whitehead_perstrand_charge.py

# Medium resolution
python whitehead_perstrand_charge.py --N 80

# Fine (paper result, ~20 min)
python whitehead_perstrand_charge.py --N 120
```

**Expected output at N=120:** Q_H ≈ 2.97 (converging toward 3.0).

**Paper reference:** Theorem `thm:perstrand_charge_three`,
Table in Sec. `sec:charge_check`.

---

### `arc_segment_analytic.py` — Arc-segment J₄ partition

Computes ρ_{J₄} arc-segment energy integrals on the analytic
Construction C field at two grid resolutions (N=64 and N=96),
assigning each grid point to its nearest arc-segment via KDTree.
Produces Table `tab:arc_segments` and the density ratios ρ_B/ρ_A.

**Usage:**

```bash
python arc_segment_analytic.py
```

Runs both N=64 (h=0.175) and N=96 (h=0.113) automatically.
Wall time: ~5 min total.

**Expected output:**
- B-segments (C→M, 30°): ~7.0% each of total J₄
- A-segments (M→C, 90°): ~26.1% each of total J₄
- Density ratio B/A ≈ 0.81

**Paper reference:** Sec. `sec:arc_segment`, Table `tab:arc_segments`,
Proposition `prop:z3_symmetry`.

---

### `gradient_flow_constrained.py` — Production gradient flow solver

The main Q_H=3 gradient-flow solver. Minimises E_geom = K·J₄ on a
periodic Cartesian grid starting from Construction C, using two-phase
Adam optimisation with per-point angle clamping (Option A topology
protection).

**Dependencies:** `bishop_frame_v2.py` must be in the same directory.

**To reproduce the Paper XVI gradient flow results (Fig. gf_fine):**

```bash
python gradient_flow_constrained.py \
  --N 192 --h 0.05 --C_star 2.5062 \
  --n_steps 4000 --lr1 3e-4 --lr2 1e-4 \
  --K_rise_eps 2.0 --delta_max_deg 9.0 \
  --log_every 10 --outdir gf_fine
```

Wall time: ~3.8 hours on CPU. Outputs: `gf_fine/log.json`,
`gf_fine/n_final.npy`.

**Key result:** Phase transition at step ~840–850, K_min ≈ 2307.5.
Phase 2 shows J₄ decaying monotonically (topology drift); the
K-minimum at Phase 1 exit is the closest numerical approach to the
Q_H=3 sector boundary.

**To reproduce the coarser N=144 run (referenced in Sec. `sec:gradient_flow`):**

```bash
python gradient_flow_constrained.py \
  --N 144 --h 0.075 --C_star 2.5062 \
  --n_steps 4000 --lr1 3e-4 --lr2 1e-4 \
  --K_rise_eps 2.0 --delta_max_deg 9.0 \
  --log_every 10 --outdir gf_coarse
```

Wall time: ~1.5 hours. K_min ≈ 2313.9 at step ~711.

**With cascade snapshots (for fraying site analysis):**

```bash
python gradient_flow_constrained.py \
  --N 192 --h 0.05 --C_star 2.5062 \
  --n_steps 4000 --lr1 3e-4 --lr2 5e-5 \
  --K_rise_eps 1.0 --delta_max_deg 5.0 \
  --log_every 10 \
  --snapshots "50,100,200,300,500,750,1000,1500,2000" \
  --outdir gf_cascade_study
```

Snapshots are saved as `n_<step>.npy` plus an automatic Phase 1 exit
snapshot tagged `[Phase1_exit]`. The Phase 1 exit snapshot is the most
theoretically significant: it represents the Q_H=3 field at the
K-minimum before vacuum drift begins.

**To warm-start from a saved checkpoint:**

```bash
python gradient_flow_constrained.py \
  --N 192 --h 0.05 --C_star 2.5062 \
  --n_steps 6000 --lr1 3e-4 --lr2 5e-5 \
  --K_rise_eps 2.0 --delta_max_deg 5.0 \
  --warm_start gf_fine/n_final.npy \
  --log_every 10 --outdir gf_continuation
```

**Full argument reference:**

| Argument | Default | Description |
|---|---|---|
| `--N` | 64 | Grid points per side |
| `--h` | 0.175 | Grid spacing (solver units; R₀=3) |
| `--C_star` | 2.5062 | Profile sharpness parameter |
| `--n_steps` | 2000 | Total gradient steps |
| `--lr1` | 3e-4 | Phase 1 (cleanup) learning rate |
| `--lr2` | 1e-5 | Phase 2 (constrained) learning rate |
| `--delta_max_deg` | 9.0 | Per-point max rotation per step in Phase 2 (degrees) |
| `--K_rise_eps` | 2.0 | K rise above K_min to trigger Phase 2 |
| `--log_every` | 10 | Steps between log lines |
| `--warm_start` | — | Path to `.npy` checkpoint to continue from |
| `--outdir` | `gf_constrained` | Output directory |
| `--snapshots` | — | Comma-separated steps for `.npy` snapshots |

**Paper reference:** Sec. `sec:saddle_prerequisite` (initial condition),
Sec. `sec:gradient_flow` (run results and interpretation),
Proposition `prop:gf_phase1`, Proposition `prop:gf_topology_loss`,
Remark `rem:gf_twophase`.

---

### `fraying_site_geometry.py` — Geometric data for the 12 fraying sites

Computes and verifies the exact geometry of all four site types that
make up the 12 fraying sites of the Q_H=3 condensate: crossings,
midpoints, distal lobes, and crossing approach/departure zones.

No external dependencies beyond NumPy. Standalone — does not require
`bishop_frame_v2.py`.

**Usage:**

```bash
python fraying_site_geometry.py
```

Wall time: ~2 s.

**Output includes:**

- Tabulated R_xy, z, φ_az, κ, tang_z for all three Z3-related points in
  each category
- 3D equidistance verification: each distal lobe lies at d=3.6263 from its
  two adjacent crossings, d=6.9293 from the far crossing
- κ ordering: crossing (0.399) > lobe (0.349) > midpoint (0.026), confirming
  S_eff hierarchy
- Tangent z-components: tang_z = 0 at crossings (exactly horizontal),
  tang_z ≈ −0.525 at midpoints (downward), tang_z ≈ +0.321 at lobes (upward)
- Full 12-site cascade ordering and Z3 coplanarity structure

**Fixes relative to inline snippet:** the original in-paper snippet
incorrectly identified the far crossing (C2, φ=180°) as adjacent to lobe L0
(φ=0°). The correct adjacent crossings are C0 and C1 at φ=±60°. The
equidistance result is unchanged but the labelling in the comment was wrong.

**Paper reference:** Proposition `prop:geom_data`, Remark `rem:torsion_contrast`,
Remark `rem:arcseg_physics`, Remark `rem:rhoJ4_anticorrelation`.

---

### `junction_flow.py` — Y-junction energy solver

Minimises E_geom = K·J₄ in the Y-junction domain (a ball centred on the
Fermat–Steiner origin) with Dirichlet boundary conditions from
Construction C. Resolves the junction energy E_∪ to seven significant
figures (Proposition `prop:ejunc_zero`).

**Dependencies:** `bishop_frame_v2.py` must be in the same directory.

**To reproduce the Paper XVI junction results:**

```bash
# Physical saddle value C*=2.5062 (converged run, ~4 min)
python junction_flow.py \
  --N 60 --C_star 2.5062 \
  --n_steps 10000 --lr 1e-4 \
  --log_every 50 --outdir junction_production_2.5062

# C*=2.0 (2000 steps, ~3 min)
python junction_flow.py \
  --N 60 --C_star 2.0 \
  --n_steps 2000 --lr 1e-4 \
  --log_every 50 --outdir junction_production_2.0

# C*=3.0 (2000 steps, ~3 min)
python junction_flow.py \
  --N 60 --C_star 3.0 \
  --n_steps 2000 --lr 1e-4 \
  --log_every 50 --outdir junction_production_3.0
```

**Expected results:**
- C*=2.5062: E_∪/E_total = 6.9×10⁻⁸ (converged at step ~2900)
- C*=2.0:    E_∪/E_total = 1.8×10⁻⁶
- C*=3.0:    E_∪/E_total = 1.8×10⁻⁹
- Scaling: E_∪ ∝ C*⁻¹⁶

**Full argument reference:**

| Argument | Default | Description |
|---|---|---|
| `--N` | 64 | Grid points per side (box is [-r_out, r_out]³) |
| `--r_outer` | 2.0 | Half-width of Cartesian box |
| `--r_boundary` | 1.6 | BC layer radius (points beyond are held fixed) |
| `--r_junc` | 0.8 | Inner measurement radius for E_junction |
| `--C_star` | 2.5062 | Profile sharpness parameter |
| `--n_steps` | 2000 | Gradient steps |
| `--lr` | 3e-5 | Learning rate |
| `--log_every` | 50 | Steps between log lines |
| `--outdir` | `junction_run` | Output directory |

**Paper reference:** Sec. `sec:what_is_needed`, Proposition `prop:ejunc_zero`,
Remark `rem:ejunc_consequence`, Eq. `eq:baryon_final`.

---

## Output files

Each solver run produces:

| File | Description |
|---|---|
| `<outdir>/log.json` | Step-by-step diagnostics (E, K, J₄, r̄, \|∇\|, clamped%) |
| `<outdir>/n_final.npy` | Final field array, shape (N,N,N,3), dtype float32 |
| `<outdir>/n_<step>.npy` | Snapshot at requested step (if `--snapshots` used) |
| `<outdir>/n_<step>.npy [Phase1_exit]` | Automatic snapshot at K-minimum |

Field arrays contain the unit vector field **n** : [0,N)³ → S², with the
south pole (0,0,−1) representing the vacuum. To load and inspect:

```python
import numpy as np
n = np.load("gf_fine/n_final.npy")   # shape (192,192,192,3)
print(n.shape, n.dtype)
# Vacuum fraction (should be low for a healthy Q_H=3 field):
print("vacuum fraction:", (n[...,2] > 0.95).mean())
```

To visualise ρ_{J₄} (the Hopf flux density), the companion HTML viewer
`hopfion_viewer.html` renders the `.npy` field using 1.5×10⁵ particles
coloured by ρ_{J₄} intensity.

---

## Physical constants used throughout

| Symbol | Value | Description |
|---|---|---|
| φ | (1+√5)/2 = 1.6180… | Golden ratio |
| λ₃ | φ⁶ = 17.9443… | Proved saddle ratio K/J₄ at Q_H=2 BPS point |
| μ | 3−φ = 1.3820… | Density-feedback coupling |
| R₀ | 3.0 | Trefoil tube major radius (solver units) |
| r₀ | 0.874 | Trefoil tube minor radius (solver units) |
| C* | 2.5062 | Physical saddle profile sharpness |
| E*_profile | 0.005041 | Analytic profile energy (L_cond=1 units, Thm. norm3) |
| E_baryon | 1.04405 · E*_profile | Baryon energy after string tension correction |

---

## Reproducing Paper XVI figures in order

1. **Sec. `sec:bishop_holonomy`** — no run needed; `bishop_frame_v2.py`
   prints holonomy = −63.0190° when imported.

2. **Sec. `sec:perstrand_winding` / `sec:charge_check`** (Q_H=3 verification):
   ```bash
   python whitehead_perstrand_charge.py --N 120
   ```

3. **Sec. `sec:gradient_flow`** (K-minimum landscape):
   ```bash
   python gradient_flow_constrained.py \
     --N 192 --h 0.05 --C_star 2.5062 \
     --n_steps 4000 --lr1 3e-4 --lr2 1e-4 \
     --K_rise_eps 2.0 --delta_max_deg 9.0 \
     --log_every 10 --outdir gf_fine
   ```

4. **Prop. `prop:ejunc_zero`** (junction energy):
   ```bash
   python junction_flow.py --N 60 --C_star 2.5062 \
     --n_steps 10000 --lr 1e-4 --outdir junction_production_2.5062
   ```

5. **Sec. `sec:arc_segment`** (arc-segment table):
   ```bash
   python arc_segment_analytic.py
   ```

6. **Prop. `prop:geom_data`** (fraying site geometry table):
   ```bash
   python fraying_site_geometry.py
   ```

7. **Sec. `sec:arc_segment` / 12-site fraying** (cascade experiment):
   ```bash
   python gradient_flow_constrained.py \
     --N 192 --h 0.05 --C_star 2.5062 \
     --n_steps 4000 --lr1 3e-4 --lr2 5e-5 \
     --K_rise_eps 1.0 --delta_max_deg 5.0 \
     --snapshots "50,100,200,300,500,750,1000" \
     --outdir gf_cascade_study
   ```
   Load the Phase 1 exit snapshot `n_<step>.npy [Phase1_exit]` into
   `hopfion_viewer.html` to visualise the 12-site fraying pattern.
