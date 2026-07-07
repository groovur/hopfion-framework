# Density-Feedback FN Hopfion — Reproduction Package

Code and data to reproduce all numerical results and figures in:

> *The Density-Feedback Faddeev–Niemi Hopfion: Exact Algebraic Predictions
> for Standard-Model Parameters* (doi:[10.5281/zenodo.19342027], 2026)
> and its companion paper *BPS Structure and Fixed-Point Theorem…* (doi:[10.5281/zenodo.19363491], 2026)

---

## Contents

| File | Purpose |
|------|---------|
| `fn_hopfion_solver.py` | Primary solver — finds the self-consistent saddle at a given R0 |
| `R0_saddle_snapshot.py` | E_geom warm-start generator — fast (~1 min) initial profile for any R0 |
| `fn_hopfion_plot.py` | Produces **Figure 1** (four-panel profile plot, any R0) |
| `fn_hopfion_plot_symmetric.py` | Produces **Figure 2** (symmetric publication-quality plot, R0=3 only) |
| `weinberg_rge_twoloop.py` | Reproduces the two-loop SM RGE check (Table in §Weinberg angle) |
| `f_fb_beta0_45200.npy` | **Converged R0=3 profile** (256×256, β\*=0.452, 24 sessions) — reference dataset |

---

## Requirements

```
python >= 3.9
numpy
scipy
matplotlib
```

Install with:
```bash
pip install numpy scipy matplotlib
```

---

## Quick start — verify the key numbers in 5 minutes

Load the provided converged profile and print all paper quantities:

```python
import numpy as np

phi = (1 + 5**0.5) / 2
lam = phi**6
mu  = 3 - phi

f    = np.load('f_fb_beta0_45200.npy')
N=256; h=0.12; R0=3.0; beta=0.452

r_ = h*(np.arange(N)+0.5); z_ = h*np.arange(N)
R, Z = np.meshgrid(r_, z_, indexing='ij')
vol = 2*np.pi*2*R*h**2; vol[:,0] *= 0.5
D0  = (R-R0)**2 + Z**2 + 1e-14
A   = 1/D0 + 1/R**2

fr = np.gradient(f,h,axis=0); fz = np.gradient(f,h,axis=1); fz[:,0]=0
sf = np.sin(f); s2=sf**2; s4=s2**2
kern = fr**2 + fz**2 + s2*A
fb   = 1/(1 + beta*kern)
fDG  = fr*(R-R0) + fz*Z
F13  = -s2/D0*fDG; F12 = s2/R*fr; F23 = s2/R*fz

J2a  = np.sum(s4*kern*vol)
J2fb = np.sum(kern*fb*vol)
J4   = np.sum((F13**2+F12**2+F23**2)*vol)
K_fb = J2a + mu*J2fb

print(f"V  = J2iso_fb/J2a = {J2fb/J2a:.8f}  (phi = {phi:.8f})")
print(f"J4/J2a  = {J4/J2a:.8f}  (WZW = {2**(4/3)/phi**5:.8f})")
print(f"K_fb/J4 = {K_fb/J4:.8f}  (WZW = {lam/2**(1/3):.8f})")
print(f"sopt^6  = {(lam*J4/K_fb)**3:.8f}  (WZW = 2.0)")
```

Expected output:
```
V  = J2iso_fb/J2a = 1.61803399  (phi = 1.61803399)
J4/J2a  = 0.22716825  (WZW = 0.22721402)
K_fb/J4 = 14.24648   (WZW = 14.24238)
sopt^6  = 1.99828     (WZW = 2.0)
```

Gaps from WZW predictions: J4/J2a −0.020%, K_fb/J4 +0.03%, sopt^6 −0.09% —
all consistent with O(h²) discretisation error at h=0.12.

---

## Reproducing Figure 1

`fn_hopfion_plot.py` runs a 10,000-step E_geom saddle-snapshot flow from a BS-ansatz
initial condition, saves the profile, and produces the four-panel figure.

```bash
# R0=3 (paper figure, ~2 minutes)
python3 fn_hopfion_plot.py --R0 3.0

# R0=4 or R0=5 (same script, different R0)
python3 fn_hopfion_plot.py --R0 4.0
python3 fn_hopfion_plot.py --R0 5.0
```

Outputs (in working directory):
- `hopfion_saddle_profile_R0{R0}.png` — the figure
- `hopfion_saddle_profile.npy` — the E_geom profile (usable as warm-start)
- `hopfion_saddle_stats.txt` — key scalar quantities

> **Note on the display:** The plot uses z ∈ [0, 10] on the horizontal axis and
> r ∈ [0, 12] on the vertical axis. The torus centre is at (z=0, r=R0). If R0 > 3
> the Hopf tube appears shifted right in the vertical panel — this is correct, not a bug.

---

## Reproducing Figure 2

`fn_hopfion_plot_symmetric.py` produces the symmetric publication-quality figure
(paper Figure 2). It reads `hopfion_saddle_profile.npy` from the working directory.

```bash
# First generate the profile (if not already done):
python3 fn_hopfion_plot.py --R0 3.0

# Then produce Figure 2:
python3 fn_hopfion_plot_symmetric.py
```

Output: `hopfion_profile_clean.png`

This script is hard-coded to R0=3 and β\*=0.452. It does not run the solver
itself — it reads the `.npy` file and computes derived fields.

---

## Reproducing the R0=3 converged result (Table 1, §Numerics)

The provided `f_fb_beta0_45200.npy` is the fully-converged reference profile
(24 solver sessions, ~200 hours total). To reproduce it from scratch:

### Step 1 — Generate E_geom warm-start (~1 minute)

```bash
python3 R0_saddle_snapshot.py --R0 3.0
```

Output: `f_R03.0_egeom.npy`

This runs 10,000 steps of gradient descent on the geometric functional
E_geom = (J2a + μ\*·J2iso)·J4 (no feedback). The E_geom saddle provides
a well-shaped initial profile that prevents collapse in the feedback solver.

### Step 2 — Run the feedback solver

**Verification run** (~45 minutes, single beta):
```bash
python3 fn_hopfion_solver.py \
    --R0 3.0 \
    --warmstart f_R03.0_egeom.npy \
    --beta_list "0.452" \
    --outdir outputs_R03_verify
```

Expected result in `outputs_R03_verify/report.txt`:
```
V    = 1.61803399  (phi = 1.61803399)
sopt = 1.000...
J4/J2a  ≈ 0.22717  (WZW 0.22721)
K_fb/J4 ≈ 14.246   (WZW 14.242)
score   < 0.001
```

**Full production run** (~8 hours, 11-beta warm-start chain):
```bash
python3 fn_hopfion_solver.py \
    --R0 3.0 \
    --warmstart f_R03.0_egeom.npy \
    --beta_list "0.3,0.35,0.4,0.42,0.44,0.44721,0.452,0.47,0.47213,0.5,0.55" \
    --outdir outputs_R03_full
```

The warm-start chain is important: each β step reshapes the profile gradually
toward the equilibrium shape at the next β. Jumping directly to β=0.452 from
a cold BS-ansatz initial condition causes the profile to collapse toward f=0.
The chain discovered β\*≈0.452 on the first run; subsequent runs used
`--beta_list "0.452"` for verification.

---

## Reproducing the R0-universality dataset (Table 2, §Universality)

### R0=4

```bash
# Step 1: E_geom warm-start (~1 minute)
python3 R0_saddle_snapshot.py --R0 4.0
# Output: f_R04.0_egeom.npy

# Step 2: Feedback solver, first session (~3 hours)
python3 fn_hopfion_solver.py \
    --R0 4.0 \
    --warmstart f_R04.0_egeom.npy \
    --beta_list "0.46,0.47,0.475,0.48,0.483,0.485,0.488,0.49,0.495,0.5" \
    --outdir outputs_R04_s1

# Step 3: Second session from best profile (~3 hours, optional)
python3 fn_hopfion_solver.py \
    --R0 4.0 \
    --warmstart outputs_R04_s1/f_beta_best.npy \
    --beta_list "0.46,0.47,0.475,0.48,0.483,0.485,0.488,0.49,0.495,0.5" \
    --outdir outputs_R04_s2
```

Expected best result: J4/J2a ≈ 0.2174 (ratio to WZW ≈ 0.957), β\* ≈ 0.483

### R0=5

```bash
# Step 1: E_geom warm-start (~1 minute)
python3 R0_saddle_snapshot.py --R0 5.0
# Output: f_R05.0_egeom.npy

# Step 2: Feedback solver (~3 hours)
python3 fn_hopfion_solver.py \
    --R0 5.0 \
    --warmstart f_R05.0_egeom.npy \
    --beta_list "1.05,1.08,1.1,1.12,1.14,1.15,1.153,1.154,1.157,1.16" \
    --outdir outputs_R05_s1

# Step 3: Second session (optional)
python3 fn_hopfion_solver.py \
    --R0 5.0 \
    --warmstart outputs_R05_s1/f_beta_best.npy \
    --beta_list "1.05,1.08,1.1,1.12,1.14,1.15,1.153,1.154,1.157,1.16" \
    --outdir outputs_R05_s2
```

Expected best result: J4/J2a ≈ 0.2316 (ratio to WZW ≈ 1.019), β\* ≈ 1.154

> **Why R0=4 and R0=5 froze after 2–3 sessions:** At larger R0 the E_geom saddle
> *is* the E_fb saddle — the feedback makes no further change to the profile shape.
> This is consistent with the self-similarity argument: as R0 increases the Hopf
> tube moves away from the r=0 axis and the 1/r² correction becomes negligible,
> so E_geom captures the full physics. The β\* values are NOT universal
> (β\* grows with R0), but J4/J2a is predicted to be universal by the WZW algebra.

---

## Reproducing the Weinberg-angle two-loop RGE check

```bash
python3 weinberg_rge_twoloop.py
```

Runs in seconds. Outputs sin²θW(M_Z) from two-loop SM running (Machacek–Vaughn 1983)
and compares to the WZW prediction 3/(8φ). Expected output:

```
sin²θW(M_Z) from two-loop SM RGE: 0.23122  (observed)
WZW prediction 3/(8phi):          0.23176
Gap:                               +0.234%  (three-loop level)
```

This establishes that Route B (RGE derivation of sin²θW from a GUT boundary
condition) does not independently confirm 3/(8φ) at two-loop order.

---

## Using R0_saddle_snapshot.py profiles with the plot scripts

Both `fn_hopfion_plot.py` and `fn_hopfion_plot_symmetric.py` can use profiles
from `R0_saddle_snapshot.py` as input for different-R0 plots. The naming
convention is:

- `R0_saddle_snapshot.py --R0 X` saves `f_R0X.0_egeom.npy`
- `fn_hopfion_plot.py` saves `hopfion_saddle_profile.npy` (always this name)
- `fn_hopfion_plot_symmetric.py` reads `hopfion_saddle_profile.npy` (hard-coded)

To plot a different R0 with `fn_hopfion_plot_symmetric.py`:
```bash
# 1. Run R0_saddle_snapshot for the desired R0
python3 R0_saddle_snapshot.py --R0 4.0
# Output: f_R04.0_egeom.npy

# 2. Rename/copy to the expected filename
cp f_R04.0_egeom.npy hopfion_saddle_profile.npy

# 3. Run the symmetric plot (note: axis labels in suptitle are hard-coded to R0=3)
python3 fn_hopfion_plot_symmetric.py
```

> **Known display issue:** `fn_hopfion_plot_symmetric.py` has R0=3 hard-coded
> in its `suptitle` and the `--R0` argument from `fn_hopfion_plot.py` is not
> threaded into it. The plots will be geometrically correct for any R0 but the
> title will still say β\*=0.452. Use `fn_hopfion_plot.py --R0 X` directly
> for R0 ≠ 3 figures with correct metadata.

---

## Solver options reference

```
fn_hopfion_solver.py:
  --R0          Torus major radius (default: 3.0)
  --outdir      Output directory (default: outputs_R0{R0})
  --warmstart   Path to .npy warm-start profile
                (auto-detects f_fb_beta0_45200.npy or previous best)
  --beta_list   Comma-separated β values, e.g. "0.40,0.45,0.452"
  --max_steps   Max gradient-descent steps per β (default: 800,000)
  --grid        SMALL=128×128 h=0.20 (testing, ~10 min)
                LARGE=256×256 h=0.12 (production, default)

R0_saddle_snapshot.py:
  --R0          Torus major radius (default: 3.0)
  (no other options; always runs 10,000 steps on 256×256)

fn_hopfion_plot.py:
  --R0          Torus major radius (default: 3.0)
  --outdir      Output directory (default: .)
```

---

## Output files

`fn_hopfion_solver.py` writes to `--outdir`:

| File | Contents |
|------|---------|
| `log.txt` | Full run log with per-step diagnostics |
| `results.csv` | One row per β: V, sopt, sopt6, J4/J2a, K_fb/J4, score, … |
| `report.txt` | Human-readable summary table |
| `f_beta{value}.npy` | Converged profile at each β |
| `f_beta_best.npy` | Profile at lowest score (use as warm-start for next session) |

---

## Expected runtimes

| Task | Grid | Time |
|------|------|------|
| E_geom warm-start (any R0) | 256×256 | ~1 min |
| Single-β verification run | 256×256 | ~45 min |
| 11-β production run (R0=3) | 256×256 | ~8 hours |
| R0=4 or R0=5 universality (2 sessions) | 256×256 | ~6 hours |
| Quick test | 128×128 (--grid SMALL) | ~10 min |
| Two-loop RGE check | — | seconds |
| Figure 1 or 2 | 256×256 | ~2 min |

All times are on a single CPU core. The solver is pure NumPy and is
not parallelised; it scales linearly with the number of β values and
grid size.

---

## Numerical values in the paper

All key numbers can be verified from `f_fb_beta0_45200.npy` using the
quick-start snippet above. The full table from Paper 1 §Numerics:

| Quantity | Numerical | WZW prediction | Gap |
|----------|-----------|----------------|-----|
| V = J2iso_fb/J2a | 1.61803399 | φ = 1.61803399 | < 10⁻⁸ |
| J4/J2a | 0.22716825 | 2^(4/3)/φ⁵ = 0.22721402 | −0.020% |
| K_fb/J4 | 14.2465 | λ/2^(1/3) = 14.2424 | +0.030% |
| sopt⁶ | 1.99828 | 2 | −0.086% |
| β\* | 0.45200 | — | — |

All gaps are consistent with O(h²) discretisation error at h=0.12
(expected ~0.04–0.1% for the 256×256 grid).
