# coherent_pair — coherent two-tube (T(2,2) cross-section) dielectric cell

Standalone solver (single file, `coherent_pair.py`, no imports from the
sibling suites) for the 2D dipole polarizability of the coherent
two-tube k-essence cell, its vertex-suppression deficit, and the
two-channel contact-packing product against the golden/geo target
`112.5/phi^10 = 0.914695`.

## Physics

The coherent background is the exact degree-2 CP^1 map

    w(z) = lam^2 / ((z - z1)(z - z2)),   z1,2 = -+ s/2,   lam^2 = s,

so each lump reduces to the unit tube `f = 2 arctan(1/rho)` of the 1D
cell (`../dielectric_cell/`) in the far-separated limit. Kinetic
density `X = 8|w'|^2/(1+|w|^2)^2` (equals `8/(1+rho^2)^2` exactly for a
single lump). Dielectric tensor of the fluctuation medium of
`L(X) = X/(1+bX)`:

    eps_tan = 1/(1+bX)^2,   eps_rad = (1-3bX)/(1+bX)^3,
    eps_ij = eps_tan delta_ij + (eps_rad - eps_tan) ghat_i ghat_j.

Anisotropy prescriptions (`--prescription`):

- `gradf` (default): `ghat` along `grad f`, `f = 2 arctan|w|`. This is
  the unique coherent generalization that reduces exactly to the
  validated 1D single-tube cell.
- `pullback`: isotropic `eps = (1-bX)/(1+bX)^3` (the pullback metric of
  a holomorphic map is conformal). Does NOT reduce to the 1D cell;
  comparison only.
- `A` / `B`: the toy superpositions of `../dielectric_cell/pair_cell/`
  (incoherent `X1+X2` / coherent gradient sum), for cross-checks.

Framework constants are computed internally to machine precision:
`b* = 0.06715131` (from `V(b*) = phi`), `C*^2 = (3/4) phi^5 / 2^(4/3)`,
`C2* = 3.43182008`, default separation `s = 2 R0/(C2* C*) = 0.96230706`.
Limit absorption `eps -> eps + i*eta`.

## Numerics and the instability

Method identical to the validated `pair_cell` suite: 9-point
cell-centered FD for `div(eps grad u) = 0`, harmonic face averaging,
Dirichlet `u = -E.r` at `|x|,|y| = L = 40`, dipole from an m=1 ring fit
on `r in [15,25]` with a growing-mode term absorbing finite-L leakage.
One LU factorization serves both field orientations.

`eps_rad` crosses zero on a contour near each core (`rho_c = 0.519`
per isolated tube at `b*`). The limit-absorption physics is resolved
only when the `|eps_rad| < eta` shell is covered by several cells;
its width is about `2 eta / |d eps_rad / d rho|`. Criterion table
(h = width/4):

| eta   | required h |
|-------|-----------|
| 1e-2  | 0.007     |
| 3e-3  | 0.002     |
| 1e-3  | 0.0007    |

Take h -> 0 before eta -> 0. On under-resolved grids `Im(alpha)`
collapses toward 0 linearly in eta or destabilizes entirely; `Re(alpha)`
is far less sensitive and grid-converges at the ~1% level already at
h = 0.1. Treat Re as the deliverable and Im as criterion-limited unless
the table above is satisfied.

Every completed solve appends one JSON line to `checkpoint.jsonl`
immediately and is skipped on rerun (safe to interrupt and resume).
Run levels strictly sequentially; do not launch grids in parallel
(BLAS contention makes everything slower). Memory: the full-domain
complex LU needs roughly 1 GB at h=0.1, 6-8 GB at h=0.05, 30+ GB at
h=0.025 (do not attempt h=0.025 full-domain on a 32 GB machine).

## Usage

Validation gates (run first, ~1 min):

    python3.11 coherent_pair.py --validate --h 0.2

Expected: V1 uniform |alpha| < 1e-9; V3 single-tube limit within ~1% of
the 1D reference at eta=1e-2 (and identical to the pair_cell suite's
single-tube value at the same h, eta); V3b coherent X matches the
single-tube X to ~3e-2 at s=30 within rho < 3 (the finite-s residue
correction of the product map is O(rho/s), so this deviation falls off
as 1/s and does not affect the s -> inf reduction).

Coarse physics pass (~2 min):

    nice -n 10 python3.11 coherent_pair.py --h 0.2 --eta 1e-2 3e-3 1e-3

Production ladder (hours; sequential; run overnight):

    nice -n 10 python3.11 coherent_pair.py --h 0.1 0.05 --eta 1e-2 3e-3 \
        --two-channel

`--two-channel` prints the vertex deficit of the same background and
the contact-packing products (square and triangular lattices) against
the golden/geo target, using Re<alpha> at the finest completed
(h, eta). Variants: `--prescription pullback|A|B`, `--form sum
--alpha-phase <rad>` for the sum-form background with a relative lump
phase, `--s <val>` / `--beta <val>` for off-framework parameters.

## Interpreting the output

`<alpha>/(2 alpha_single_1D)` near 1 means the tubes polarize
independently; a strong suppression (|ratio| well below 1/2) would be
required for contact packing to close the dilute Maxwell-Garnett route
on its own. The two-channel product tests instead whether polarization
(eps_eff) and vertex suppression (`<1/(1+bX)>`) jointly reproduce
`112.5/phi^10` at contact packing with no free density. The toy
prescriptions A and B bracket the target at about +6% / -5%
(see `../../paper3_support/two_channel_bracket.py`); the `gradf`
coherent background is the physically preferred single value.
