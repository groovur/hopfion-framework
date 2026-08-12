# Nonlinear rectification test ("Tacoma Narrows" hypothesis)

Script: `rectify.py`. Raw numeric output: `results_raw.json`. Shooting solver
and background functions copied verbatim from `../../dielectric_cell/common.py`
(not modified); `beta_star = 0.06715131` and `rho_c(beta_star) = 0.51913431`
reproduced and asserted against that file's validated values before use.

## Model

Background single tube `X0(rho) = 8/(1+rho^2)^2`. A measurement event is a
weak oscillatory probe driving the k-essence fluctuation,
`delta X(rho,tau) = A*f(rho)*cos(omega*tau)`, `A << 1`. Time-averaging the
nonlinear permittivities `eps_tan(X) = 1/(1+bX)^2`,
`eps_rad(X) = (1-3bX)/(1+bX)^3` over one drive period gives, to `O(A^2)`,

```
<eps(X0 + delta X)>_tau = eps(X0) + (A^2/4) f(rho)^2 eps''(X0)
```

(the `O(A)` term vanishes under `<cos>_tau = 0`; `<cos^2>_tau = 1/2` produces
the quarter-amplitude coefficient). This is the rectification term tested.

## R1 — closed-form second derivatives and sign

With `u = b*X` (sympy-verified against `diff(eps, X, 2)` directly, not just
asserted):

```
eps_tan''(X)  =  6 b^2 / (1+bX)^4
eps_rad''(X)  =  6 b^2 (5 - 3bX) / (1+bX)^5
```

**Signs.**

- `eps_tan''(X) > 0` for every `X` on the physical domain (`1+bX>0` always
  makes the denominator positive; numerator `6b^2>0`). No sign change.
- `eps_rad''(X)` changes sign at `u = bX = 5/3`: positive for `u<5/3`,
  negative for `u>5/3`. On this background `X0` is maximal at `rho=0`,
  giving `u_max = b*X0(0) = 0.53721`, well below `5/3 = 1.6667`. **`eps_rad''`
  therefore stays positive over the entire physical domain here** — same
  sign as `eps_tan''`, everywhere.
- At the `eps_rad=0` shell `rho_c = 0.51913` (`u = bX = 1/3` by definition of
  that shell), `eps_rad''(rho_c) = 0.025682 > 0`: the zero of `eps_rad`
  itself is *not* a sign-change locus of `eps_rad''`. The two loci (`eps_rad
  = 0` at `u=1/3`, `eps_rad''=0` at `u=5/3`) are distinct and, for this
  background/`b*` combination, the second is never reached.

So both rectification shifts `<delta eps_tan>` and `<delta eps_rad>` are
**strictly positive** (for any `A^2 f(rho)^2 > 0`) everywhere on the domain
probed here — nonlinear rectification raises both permittivities on average,
it does not flip sign near the resonant shell.

## R2 — response shift, profile-shaped drive

Drive `f(rho) = X0(rho)/8` (natural profile shape), `eta = 1e-4`,
`beta = beta_star`. Baseline (`A=0`): `alpha_pol(0) = -0.97104087 +
0.01366438j`, matching the independently validated `dielectric_cell` run at
the same `(beta,eta)` to `<1e-6`.

| A^2 | Re alpha_pol(A) | Im alpha_pol(A) | Re r(A) | Im r(A) |
|---|---|---|---|---|
| 1e-2 | -0.97101975 | 0.01365801 | -2.1842e-5 | 6.2455e-6 |
| 1e-3 | -0.97103876 | 0.01366374 | -2.1842e-6 | 6.2461e-7 |
| 1e-4 | -0.97104066 | 0.01366431 | -2.1842e-7 | 6.2460e-8 |

`r(A) = [alpha_pol(A)-alpha_pol(0)]/alpha_pol(0)`. Log-log slope of `|r|` vs
`A^2` over the three points: **1.0000** — confirms `r ~ A^2` cleanly (pure
quadratic rectification, as the perturbative construction requires).

```
C_rect (profile-shaped)  =  r/A^2  =  -0.0021842 + 0.0006246 j     (complex, ~constant across all three A^2)
```

`Re(C_rect) < 0`: the profile-shaped rectification *reduces* `Re(alpha_pol)`
in magnitude (pushes it toward zero, since `alpha_pol(0)<0`).

## R3 — resonant-shell drive

Drive centered on the `eps_rad=0` shell, `f(rho) =
exp(-(rho-rho_c)^2/(2*0.1^2))`, `rho_c = 0.51913`, same `A^2` grid, same
`eta`. Log-log slope of `|r|` vs `A^2`: **1.0000** (again clean `A^2`
scaling).

```
C_rect (shell-centered)  =  -0.0034378 + 0.0011812 j
```

```
enhancement = |C_rect_shell| / |C_rect_profile| = 1.600
```

**Enhancement is 1.6x, not >10x.** The resonant shell amplifies rectification
somewhat (as expected — the shell is where `eps_rad` is most sensitive) but
by an O(1) factor, not by orders of magnitude. **This does not support a
literal "Tacoma Narrows" resonant blow-up**: `eps_rad''` itself is smooth,
finite, and non-singular at `rho_c` (see R1 — the zero of `eps_rad` is not a
pole or sign-change of `eps_rad''`), so there is no resonant divergence for
a drive centered there to amplify; the shell is simply a mildly more
"productive" region than the profile average, order-1 same as anywhere else
on the tube.

## R4 — sign chain and implied amplitude

Sign chain: `eps_eff(n,alpha) = 1 + 2 n alpha/(1-n alpha)`, so
`d(eps_eff)/d(alpha) = 2n/(1-n alpha)^2 > 0` for any real `alpha` (denominator
is a square) — `eps_eff`, and hence `alpha_inv` (proportional to `eps_eff` by
construction), is monotonically *increasing* in `alpha_pol`.

Residual: `alpha_inv(formula) - alpha_inv(measured) = 137.035998 -
137.035999177 = -1.177e-6` (formula low). To move the rectification-corrected
value up to match measured, need `Delta(alpha_pol) > 0` (since
`d(eps_eff)/d(alpha)>0`). But `alpha_pol(0) = Re(alpha0) = -0.97104 < 0`
(anti-screening) and `Delta(alpha_pol) = alpha_pol(0)*C_rect*A^2`, so
`Delta(alpha_pol)>0` requires `alpha_pol(0)*C_rect > 0`, i.e. (since
`alpha_pol(0)<0`) **`C_rect` must be negative**.

The R2 profile-shaped `C_rect` computed above is `Re(C_rect) = -0.0021842`
— **negative, the correct sign** for rectification to move the coupling in
the direction needed to explain the observed residual. This is a sign
consistency check only, not a fit.

Using the fractional-proportionality relation (`alpha_inv ∝ eps_eff`, so
relative shifts are equal regardless of the unknown proportionality
constant) and `n = 0.0459`:

```
eps_eff0             = 0.914662
d(eps_eff)/d(alpha)  = 0.084133
frac_needed           = [alpha_inv(measured)-alpha_inv(formula)] / alpha_inv(formula) = 8.589e-9
coefficient           = d(eps_eff)/d(alpha) * alpha_pol(0) * Re(C_rect) / eps_eff0 = 1.9509e-4   (per unit A^2)

A^2_implied  =  frac_needed / coefficient  =  4.403e-5
A_implied    =  sqrt(A^2_implied)          =  0.00664
```

**This inferred `A` is physically tiny (`A ~ 6.6e-3 << 1`)** — well inside the
weak-drive perturbative regime this whole calculation assumes is valid (the
`O(A^4)` corrections neglected above would be `O(A^2) ~ 4e-5` relative to the
`O(A^2)` term itself, i.e. utterly negligible at this amplitude). It is
labeled here explicitly as an **inference, not a fit**: `A` was solved for
given the computed `C_rect`, not tuned to reproduce the residual by
adjusting any model parameter.

## Answer

**Yes, nonlinear rectification biases the coupling — in the correct
direction — but the effect requires only a tiny drive amplitude, and the
"Tacoma" resonance is not singular.** Concretely:

1. The nonlinear vacuum response *does* rectify a periodic measurement drive
   into a genuine `O(A^2)` DC shift of both `eps_tan` and `eps_rad`, confirmed
   both in closed form (R1) and numerically in the full shooting solver (R2,
   R3: clean `A^2` scaling, log-log slope `1.0000` in both cases) — this is
   the qualitative mechanism the linear (`../drift.py`) analysis could not
   produce and explicitly flagged as the remaining open possibility.
2. The rectification coefficient `C_rect` for the physically motivated
   profile-shaped drive is `-0.002184 + 0.000625j`; its real part is negative,
   which is *exactly* the sign needed to shift `alpha_inv` upward (toward the
   measured value) given the anti-screening sign of `alpha_pol` and the
   monotonic `eps_eff(alpha_pol)` relation — the mechanism points the right
   way, not by construction but as a nontrivial output of the sign chain.
3. Resonant enhancement is real but modest: centering the drive on the
   `eps_rad=0` shell increases `|C_rect|` by only `1.6x`, not an order of
   magnitude or more, because `eps_rad''` is smooth and non-singular at that
   shell (its own zero is at `u=1/3`, but `eps_rad''` only changes sign at
   the different, unreached point `u=5/3`) — there is no genuine resonant
   divergence for a "Tacoma Narrows" blow-up to ride on, only a mild
   geometric enhancement.
4. Given the computed `C_rect` and `n=0.0459`, explaining the entire observed
   four-term-formula residual (`1.177e-6` in `alpha_inv`, i.e.
   `8.6e-9` relative) via this mechanism alone would require drive amplitude
   `A ~ 0.0066` — physically tiny, safely inside the perturbative regime used
   to derive the `O(A^2)` formula. This is reported as a scale check only
   (not a fit, not a claim that this mechanism *is* the source of the
   residual): it shows the required amplitude is not absurd, so the
   mechanism is quantitatively plausible as a contributor, but nothing here
   establishes that it is the actual or complete explanation.
