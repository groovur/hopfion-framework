# Measurement drift in the SU(2)_3 probe-pair monodromy coupling

Script: `drift.py`. Random seed fixed (`20260812`). Raw numeric output: `results_raw.json`.

## Model

Two `j=1/2` anyons in the SU(2)_3 (Fibonacci-like) sector fuse to channels
`c in {0,1}`, quantum dimensions `d_0=1`, `d_1=phi`, Born weights
`p0 = 1/phi^2 = 0.381966`, `p1 = 1/phi = 0.618034`. Monodromy eigenvalues
`M_0 = exp(-i*108deg)`, `M_1 = exp(+i*36deg)`. The static identification is
`alpha_inv = 360 * Re Tr[rho U_mono]`.

For the Born-diagonal state, `Tr[rho_born U] = p0*M0 + p1*M1` was checked at
50-digit precision: its imaginary part vanishes exactly and its real part
equals `1/phi^2` exactly, giving `alpha_inv_born = 137.50776405...`.

A single interferometer pass with fringe parameter `t in (0,1]` measures the
monodromy via Kraus operators diagonal in the fusion basis,

```
K_+ = diag( sqrt((1+t*Re M0)/2), sqrt((1+t*Re M1)/2) )
K_- = diag( sqrt((1-t*Re M0)/2), sqrt((1-t*Re M1)/2) )
```

`K_+^2 + K_-^2 = I` was verified at `t = 1, 0.5, 0.1, 0.01`. Because `rho`
starts diagonal and every Kraus operator here is diagonal, `rho_n` stays
diagonal for all `n`: the full dynamics collapses onto a scalar classical
process on `p0(n) = Pr(channel 0)`, with `Re Tr[rho_n U] = b + p0(n)*(a-b)`
(`a = Re M0 = -0.309017`, `b = Re M1 = 0.809017`), an affine, monotone
function of `p0(n)`. All reported statistics of `360*Re x_n` are this affine
map applied to the statistics of `p0(n)`.

Relaxation ("Tacoma" competition): after each pass, `p0 -> (1-g)*p0 + g*p0_born`.

## R1 — collapse statistics (no relaxation)

`N = 10^4` trajectories, `n = 1000` passes, `t in {1, 0.5, 0.1, 0.01}`.

Martingale check: `E[p0(n)]` (equivalently `E[Re x_n]`) is a mathematical
identity here, not an approximation — proved exactly in R4 below, for any
`t`, at every single pass. The Monte Carlo values below fluctuate around
`alpha_inv_born = 137.508` by Monte Carlo noise only.

Collapse-endpoint check: a trajectory collapses to a corner (`p0 -> 0` or
`p0 -> 1`) at a rate controlled by `t`. `mean_dist_to_boundary =
mean(min(p0, 1-p0))` measures how far trajectories are from having collapsed
after 1000 passes:

| t | mean_dist_to_boundary @n=1000 | frac -> channel 0 | frac -> channel 1 | Born p0 / p1 |
|---|---|---|---|---|
| 1.0 | 3e-16 (fully collapsed) | 0.3776 | 0.6224 | 0.38197 / 0.61803 |
| 0.5 | 5e-48 (fully collapsed) | 0.3855 | 0.6145 | 0.38197 / 0.61803 |
| 0.1 | 0.038 (partially collapsed) | 0.3715 | 0.6285 | — (not converged yet) |
| 0.01 | 0.375 (essentially uncollapsed) | 0.0828 | 0.9172 | — (not converged yet) |

For `t=1.0` and `t=0.5`, 1000 passes are enough to fully collapse essentially
every trajectory, and the observed channel-0 frequencies (0.3776, 0.3855)
agree with the Born value 0.38197 to within Monte Carlo error (`sqrt(p0
p1/N) = 0.0049`, i.e. 0.9 sigma and 0.5 sigma respectively) — this is the
expected verification. For `t=0.1` and especially `t=0.01`, 1000 passes are
*not* enough to collapse the ensemble (number of passes to collapse scales
roughly as `1/t^2` for weak measurement), so the raw `p0>0.5` frequencies at
`n=1000` are not yet meaningful Born-rule checks for those two columns; they
are listed for completeness, not as a discrepancy.

## R2 — measured-value drift, no relaxation

Ensemble mean and variance of `360*Re x_n` vs `n` (same runs as R1). Collapse
endpoint values: `360*Re M0 = 360*cos(108deg) = -111.246`, `360*Re M1 =
360*cos(36deg) = 291.246`.

**t = 1.0** (collapses almost entirely within ~20 passes):

| n | 0 | 1 | 5 | 20 | 100 | 1000 |
|---|---|---|---|---|---|---|
| mean | 137.51 | 137.39 | 137.67 | 139.30 | 139.27 | 139.27 |
| var | 0 | 13224 | 30324 | 37867 | 38073 | 38073 |

**t = 0.1** (still collapsing at n=1000):

| n | 0 | 1 | 5 | 20 | 100 | 1000 |
|---|---|---|---|---|---|---|
| mean | 137.51 | 137.51 | 137.66 | 138.47 | 139.14 | 141.55 |
| var | 0 | 113 | 561 | 2159 | 9007 | 33228 |

**t = 0.01** (barely moved from the Born state by n=1000):

| n | 0 | 1 | 5 | 20 | 100 | 1000 |
|---|---|---|---|---|---|---|
| mean | 137.51 | 137.50 | 137.49 | 137.42 | 137.47 | 137.71 |
| var | 0 | 1.13 | 5.55 | 22.4 | 111.7 | 1106.2 |

In every case the mean wanders within Monte Carlo noise around 137.508 (no
drift), while the variance grows monotonically from 0 and saturates once the
ensemble has collapsed (saturation value `38073`-`38376`, consistent across
`t`, set by the two-point Born-weighted distribution over the endpoints
`{-111.2, 291.2}`: `Var = p0*p1*(360*(a-b))^2 = 0.382*0.618*402.49^2 =
38260`, matching the fully-collapsed columns above to Monte Carlo precision).
**Variance growth-then-saturation, not a mean shift, is the entire
measurement effect with no relaxation** — exactly as the martingale property
requires.

## R3 — rethermalization competition (measurement + relaxation)

`N = 2000` trajectories per `(t,g)` cell, run for `min(20/g, 60000)` passes,
steady-state mean/std estimated by averaging over the last 30% of the run.

| t | g | steps | steady mean | bias (mean-137.508) | steady std | bias/SE |
|---|---|---|---|---|---|---|
| 1.0 | 0.1 | 200 | 141.740 | 4.232 | 138.207 | 1.37 |
| 1.0 | 0.01 | 2000 | 139.025 | 1.517 | 183.799 | 0.37 |
| 1.0 | 0.001 | 20000 | 138.341 | 0.834 | 193.623 | 0.19 |
| 1.0 | 0.0001 | 60000 | 142.122 | 4.614 | 194.080 | 1.06 |
| 0.5 | 0.1 | 200 | 137.143 | -0.365 | 91.195 | -0.18 |
| 0.5 | 0.01 | 2000 | 140.470 | 2.962 | 161.348 | 0.82 |
| 0.5 | 0.001 | 20000 | 141.349 | 3.841 | 188.315 | 0.91 |
| 0.5 | 0.0001 | 60000 | 138.842 | 1.334 | 194.301 | 0.31 |
| 0.1 | 0.1 | 200 | 137.359 | -0.149 | 21.756 | -0.31 |
| 0.1 | 0.01 | 2000 | 137.491 | -0.017 | 66.729 | -0.01 |
| 0.1 | 0.001 | 20000 | 138.799 | 1.291 | 138.953 | 0.42 |
| 0.1 | 0.0001 | 60000 | 137.404 | -0.104 | 182.787 | -0.03 |
| 0.01 | 0.1 | 200 | 137.568 | 0.060 | 2.202 | 1.22 |
| 0.01 | 0.01 | 2000 | 137.696 | 0.188 | 7.505 | 1.12 |
| 0.01 | 0.001 | 20000 | 137.384 | -0.124 | 23.459 | -0.24 |
| 0.01 | 0.0001 | 60000 | 135.813 | -1.695 | 67.944 | -1.12 |

"bias/SE" uses `SE = std/sqrt(N)` with `N=2000`, a lower bound on the true
standard error (steady-state snapshots are autocorrelated, so real SE is
larger and true significance is even smaller than shown). **Every bias
value across the full 4x4 grid is within about 1.4 estimated standard
errors of zero — no statistically significant bias anywhere on the grid.**

Scaling of the fluctuation amplitude (std), fit on the small-`t` cells where
the perturbative regime applies cleanly (see R4):

- `std` vs `t` at fixed `g=0.1`, using `t in {0.01, 0.1}`: log-log slope
  `0.995` (i.e. `std ~ t^1`).
- `std` vs `g` at fixed `t=0.01`: log-log slope `-0.496` (i.e.
  `std ~ g^-1/2`).

At larger `t` (0.5, 1.0) the fluctuation amplitude *saturates* rather than
continuing to grow — it approaches the bounded two-point maximum
(`~360*|a-b|/2 ~ 201`, the largest possible std of any distribution
supported on `[-111.2, 291.2]`), because the per-pass measurement pull then
dominates the per-pass relaxation pull strongly enough to nearly fully
collapse the state between relaxation kicks. A naive slope fit across the
*whole* `t` grid at `g=0.1` gives `0.91` and across the whole `g` grid at
`t=1.0` gives `-0.05` — both artifacts of this saturation, not the
underlying `t^1 g^-1/2` law, which only holds in the small-`t`, small-signal
regime.

## R4 — exact / small-t expansion (sympy)

**Exact martingale property.** `E[p0'] - p0`, computed symbolically to all
orders in `t` for the measurement step alone (not a truncated series),
simplifies to exactly `0`. This is an algebraic identity, true for every
`t in (0,1]`, not just to `O(t^2)`.

**Exact fixed point under measurement + relaxation.** Since one measurement
pass preserves `E[p0]` exactly (any `t`) and relaxation maps
`E[p0] -> (1-g)*E[p0] + g*p0_born`, the fixed point of the mean map is the
solution of `(1-g)*p0 + g*p0_born = p0`, which sympy solves to `p0 = p0_born`
for any `g != 0` and any `t`. **The bias of the mean is exactly zero at
every order in `t`, not just to leading order** — measurement (martingale)
plus linear relaxation toward the Born state cannot produce a biased steady
state, by construction of the model.

**Leading-order variance injected per pass.** `Var(p0' | p0)`, series-expanded
in `t`, has vanishing `O(1)` and `O(t)` terms; the leading term is

```
Var(p0' | p0) = t^2 * p0^2 * (1-p0)^2 * (a-b)^2 + O(t^4)
```

(sympy-verified exact factorization: `t2_coeff_factored = p0**2*(a-b)**2*(p0-1)**2`).
At `p0 = p0_born` this evaluates to `0.069660...` (the coefficient of `t^2`).

**Steady-state variance balance.** Treating variance injection (per pass,
`O(t^2)`) against relaxation damping (`(1-g)^2` per pass) as a linear
recursion gives the closed form

```
Var_ss(p0) = (1-g)^2 * f(p0_born,t) / (1 - (1-g)^2),   f = 0.069660... * t^2
```

whose small-`g` leading behavior is `Var_ss(p0) ~ f(p0_born,t) / (2g) =
0.03483 * t^2 / g`. Converting to `360*Re x` via `Var(360 Re x) = 360^2
(a-b)^2 Var(p0)` gives

```
sigma_ss(360 Re x)  ~  360*|a-b| * sqrt(0.03483 * t^2 / g)
                     =  402.49 * 0.1866 * t / sqrt(g)
                     ~  75.1 * t / sqrt(g)      (leading order, small t, small g)
```

This closed-form prediction reproduces both empirical slopes above (`t^1`,
`g^-1/2`) and matches the small-`t` numeric table to ~10-20%: e.g. at
`t=0.01, g=1e-4`, predicted `sigma = 75.1`, observed `67.9`; at `t=0.01,
g=1e-1`, predicted `2.19`, observed `2.20`. **No bias term of any order in
`t` appears anywhere in this derivation** — consistent with the exact
martingale result above, not an independent check of it (the two are the
same fact: `E[p0'']=p0_born` for all `n` is a linear-algebra identity that
this variance calculation does not touch).

## Answer

No. Repeated measurement of the fusion-channel monodromy does not bias the
measured coupling `alpha_inv = 360*Re Tr[rho U]` toward any secular drift or
resonant pumping, under any of the conditions probed (`t` from 0.01 to 1,
with or without Born-state relaxation at rates `g` from `1e-4` to `1e-1`):
the ensemble mean of `360*Re x_n` is an exact martingale, proved
symbolically to all orders in `t`, and stays pinned at `137.508` before,
during, and after relaxation reaches steady state — every bias measured
across the full 16-cell `(t,g)` grid is statistically consistent with zero
(under ~1.4 sigma). What *does* happen is a pure fluctuation effect: with no
relaxation, the variance of the measured value grows from a delta function
at the Born value toward a saturation set by the two Born-weighted collapse
endpoints (`-111.2`, `291.2`); with relaxation competing against
measurement, the ensemble reaches a genuine non-Born-collapsed steady state
whose fluctuation amplitude obeys `sigma ~ t/sqrt(g)` in the small-signal
regime (closed form and exponents both derived and numerically confirmed),
saturating to a bounded maximum when the per-pass measurement pull
dominates the relaxation pull (`t` of order 1, small `g`). The "Tacoma
Narrows" mechanism as posed — a resonant *shift* of the measured coupling —
does not occur in this linear model; producing an actual bias would require
a nonlinear vacuum response (relaxation that is not a simple linear pull
toward the Born state), which is outside the scope of the model tested
here. For context only, no fitting performed: the fluctuation amplitudes
found here (order 1-200 in `alpha_inv` units, i.e. absolute scale `1e0` to
`1e2`) are many orders of magnitude larger than the quoted experimental
`alpha_inv` residual scale (`~1e-6` absolute) or the Cs/Rb tension
(`~1e-9` relative) — this model, as specified, says nothing quantitative
about those numbers; it answers only the qualitative bias-vs-fluctuation
question posed.
