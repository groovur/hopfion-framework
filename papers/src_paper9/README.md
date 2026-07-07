# `chsh_transition_curve.py`

Validation and exploration tool for the CHSH$(f,\beta)$ transition curve of the
density-feedback Faddeev–Niemi Hopfion condensate framework.

Reproduces and validates all numerical claims of:

> F. Manfredi, *"FN Hopfion Condensate — Paper IX: Division Algebras, the Born
> Rule, and the Tsirelson Bound from Anisotropic Flux Suppression"* (2026).

---

## Requirements

```
python >= 3.9
numpy
scipy
matplotlib  (optional, for --plot)
```

Install dependencies:

```bash
pip install numpy scipy matplotlib
```

---

## Quick start

```bash
# Reproduce Paper IX Table 1 (f transition, no feedback)
python chsh_transition_curve.py

# Same table with density feedback β* = 0.452 (Paper VIII value)
python chsh_transition_curve.py --beta 0.452

# Validate ALL numerical claims in Paper IX
python chsh_transition_curve.py --validate

# Show photon/Stokes results and the post-quantum conjecture value
python chsh_transition_curve.py --photon
```

---

## What the script computes

### The model

The $f$-parameterised condensate probability:

```
p_f(θ) = 1 / (1 + tan^(2f)(θ/2))
```

| Value | Physical meaning |
|-------|-----------------|
| `f = 1` | Standard QM Born rule (`ℂ`-valued, `p = cos²(θ/2)`) |
| `f = 2` | Condensate pcond (`ℍ`-valued, quaternionic) |
| `f → ∞` | Projective measurement; pcond lepton ceiling `8/3 ≈ 2.667` |

The signed correlation: `x_f(θ) = 2·p_f(θ) − 1`

With density feedback `β > 0`:

```
f_eff(θ; f₀, β) = f₀ / (1 + β·cos θ)
```

- `θ ≈ 0` (aligned): `f_eff = f₀/(1+β) < f₀` — feedback softens measurement
- `θ ≈ π` (anti-aligned): `f_eff = f₀/(1−β) > f₀` — no feedback, sharpens
- `θ = π/2`: `f_eff = f₀` — no correction at perpendicular

**Requires `β < 1`** for positivity at all angles. Paper VIII value: `β* = 0.452`.

### Rescaling conventions

The script returns **two** rescaled CHSH values:

| Name | Formula | Used in |
|------|---------|---------|
| Fixed rescaling | `CHSH_raw / (3/8)` | Paper IX Table 1, all `f` |
| `f`-dependent | `CHSH_raw / ⟨sin^{2f} θ⟩_2D` | Theorem `7.2` statement |

`⟨sin^{2f} θ⟩_2D = Γ(f+½) / (√π · Γ(f+1))` is computed exactly via Gamma functions.

At `f = 2`: both agree (both equal `8/3`). At other `f`, they differ.

### Photon/Stokes results

The Stokes/Malus observable `cos(2(φ_n − α))` gives (Theorem `thm:stokes`):

```
E^Stokes(ψ) = −½·cos(2ψ)     [Parseval factor ½ for the m=2 mode of L²(S¹)]
```

CHSH raw = `√2`; after standard QM normalisation (`×2`): CHSH = `2√2 = Tsirelson`.

The post-quantum conjecture (`9.18`) predicts:

```
CHSH = 4·x_cond(π/4) = 8√2/3 ≈ 3.771   [> Tsirelson, no-signaling satisfied]
```

---

## Command-line reference

```
python chsh_transition_curve.py [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--validate` | off | Reproduce and check all Paper IX numerical claims |
| `--scan` | off | Scan CHSH over a `(f, β)` grid |
| `--photon` | off | Print Stokes/Malus and post-quantum conjecture values |
| `--f F[,F,...]` | `1.0,...,10.0` | Comma-separated `f` values |
| `--beta B` | `0.0` | Single density-feedback `β` value |
| `--beta-values B[,B,...]` | `0.0,0.452` | Multiple `β` values for `--scan` |
| `--N N` | `1_000_000` | Monte Carlo sample count |
| `--seed S` | `42` | Random seed |
| `--tol T` | `0.005` | Tolerance for `--validate` pass/fail |
| `--plot` | off | Save `chsh_curve.png` |

---

## Usage examples

### Reproduce Paper IX Table 1

```bash
python chsh_transition_curve.py
```

Output matches the table in §7.3 (pcond rescaling `×8/3`, `β = 0`).

### Paper VIII comparison point

```bash
python chsh_transition_curve.py --f 4.0 --beta 0.452
```

Reproduces `CHSH_rescaled ≈ 2.60` (Paper VIII Monte Carlo: `≈ 2.611`).

### Full validation suite

```bash
python chsh_transition_curve.py --validate
```

Checks 28 claims from Paper IX (12 analytical, 16 Monte Carlo) at `N = 1 000 000`.
All 28 pass at tolerance `0.005`. For tighter checks:

```bash
python chsh_transition_curve.py --validate --N 5000000 --tol 0.002
```

### Scan the (f, β) surface

```bash
python chsh_transition_curve.py --scan --f 2.0,4.0,6.0,10.0 --beta-values 0.0,0.1,0.2,0.3,0.452
```

### Photon sector results

```bash
python chsh_transition_curve.py --photon
```

Prints the two CHSH ceilings (lepton `8/3` and photon `2√2`), the Stokes
exact derivation, and the post-quantum conjecture value `8√2/3 ≈ 3.771`.

### Save a plot

```bash
python chsh_transition_curve.py --plot
python chsh_transition_curve.py --beta 0.452 --plot
```

Saves `chsh_curve.png` with the `f`-transition curve, Tsirelson bound,
lepton ceiling `8/3`, and classical bound labelled.

---

## Paper IX numerical targets

All values reproduced to within Monte Carlo noise at `N = 1 000 000`:

| Paper claim | Value | Source |
|------------|-------|--------|
| `A₁ = 4 − 2√2` | `1.1716` | Theorem `3.3` |
| Parseval sum `= 1/√2` | `0.70711` | Theorem `3.3` |
| `C(π/8)` at `f=2` | `0.6415` | Theorem `4.1` |
| `C(3π/8)` at `f=2` | `0.2445` | Theorem `4.1` |
| CHSH_raw at `f=2` | `0.886` | Theorem `4.1` |
| CHSH_rescaled at `f=2` | `2.363` | Theorem `4.1` |
| CHSH_raw at `f=1` | `0.6533` | Theorem `7.2` |
| CHSH_rescaled (`f`-dep) at `f=1` | `1.307` | Theorem `7.2`(iii) |
| CHSH_rescaled (fixed) at `f=4` | `2.613` | Table 1 |
| `C(π/8)` at `f=2, β*` | `0.6353` | Theorem `8.3`(ii) |
| `C(3π/8)` at `f=2, β*` | `0.2432` | Theorem `8.3`(ii) |
| CHSH_rescaled at `f=2, β*` | `2.341` | Theorem `8.3`(ii) |
| CHSH_rescaled at `f=4, β*` | `2.60` | Theorem `8.3`(iii) |
| Stokes CHSH raw | `√2` | Theorem `thm:stokes` (exact) |
| Stokes CHSH (`×2`) | `2√2` | Theorem `thm:stokes` (exact) |
| Post-quantum CHSH | `8√2/3 ≈ 3.771` | Conjecture `9.18` (exact) |

---

## Notes on conventions

**Fourier coefficient ratio sign.** Theorem `3.3` in Paper IX states
the geometric ratio as `A_{2k+3}/A_{2k+1} = -(3−2√2) = 2√2−3 ≈ −0.172`.
This is confirmed numerically (`A₃/A₁ = −0.17157` from direct integration).
The expression `(√2−3)^k ≈ (−1.586)^k` that appears in some drafts is incorrect;
the correct ratio is `(2√2−3)^k ≈ (−0.172)^k`.

**Rescaling.** The paper's Table 1 consistently uses the fixed pcond rescaling
`÷(3/8) = ×(8/3)` for all `f`. The `f`-dependent rescaling
`÷⟨sin^{2f}⟩_2D` is used only in the mathematical statement of Theorem `7.2`.
These agree at `f = 2` (both give `8/3`) and diverge at other `f`.

**Density feedback sign.** The correct form is `f_eff = f₀/(1 + β·cos θ)`
with **signed** cosine, matching Paper VIII and reproduced by the proof values
`C(π/8) = 0.6353` in Theorem `8.3`. The `|cos θ|` form gives
`C(π/8) ≈ 0.588`, which does not match the paper.

---

## Archived version

The code for this script is permanently archived at Zenodo:
DOI: [10.5281/zenodo.20075227](https://doi.org/10.5281/zenodo.20075227)
