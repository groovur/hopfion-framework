# pair_cell: 2D anisotropic-dielectric polarizability of a two-tube pair

2D electrostatics solver for div(eps(x,y) . grad u) = 0 with a uniform
applied field, computing the induced dipole polarizability alpha (per unit
length) of k-essence tube backgrounds. Extends the 1D radial shooting code
in `../` (see `../common.py`, `../results.md`) to non-circular geometries,
specifically the two-tube pair at separation s = 2*R0/(C2*.C*) = 0.96230706.

Run everything with `python3.11` (numpy/scipy).

## Physics conventions

- Background: X(rho) = 8/(1+rho^2)^2 per tube; eps_tan = 1/(1+beta*X)^2,
  eps_rad = (1-3*beta*X)/(1+beta*X)^3, beta = beta* = 0.06715131 (root of
  (15/8)arctan(sqrt(8b))/sqrt(8b) = phi). Overall 1/phi^6 normalization is
  dropped (eps -> 1 at infinity), matching the 1D code.
- Tensor: eps_ij = eps_tan delta_ij + (eps_rad - eps_tan) ghat_i ghat_j,
  ghat = g/|g| with g = sqrt(X1) rhat1 + sqrt(X2) rhat2 for the pair
  (radial for a single tube). Where |g| = 0 the tensor is isotropic
  eps_tan.
- Pair models: Model A (incoherent) uses X_A = X1 + X2 in eps_tan/eps_rad;
  Model B (coherent gradient) uses X_B = |g|^2, which includes the cross
  term 2*sqrt(X1 X2) rhat1.rhat2.
- Limit absorption: eps_tan and eps_rad each get +i*eta; alpha(eta) is
  computed at finite eta.
- Normalization gate: isolated step cylinder (eps1, radius R) gives
  alpha = R^2 (eps1-1)/(eps1+1); eps1=3, R=1 gives +0.5.

## Numerical method

Cell-centered finite differences on [-L,L]^2 (L=40), 9-point stencil with
harmonic face averaging of the diagonal tensor components and arithmetic
averaging of eps_xy; complex sparse direct solve (SuperLU). Dirichlet BC
u = -E.r on the outer boundary. The dipole is read from a ring
r in [15,25] by least-squares m=1 angular fit that includes 1/r, 1/r^3,
and a growing r term; the growing term absorbs the finite-L Dirichlet
leakage and the 1/r coefficient is the dipole moment.

The pair problem is mirror-symmetric about both axes, so production runs
on fine grids use a quadrant-reduced solver (4x fewer unknowns): E along
the pair axis has u odd in x/even in y, E transverse has u even in x/odd
in y, and the mirror ghost cells fold back into the matrix with the
corresponding sign.

## Files

- `common2d.py` - solver library: tensor constructors (single tube, pair
  models A/B), full-domain assembly (`assemble_system`, LU reuse across
  field directions via `solve_system_multi`), quadrant assembly
  (`assemble_system_quadrant`, `solve_quadrant`), dipole extraction
  (`extract_dipole`, `extract_dipole_quadrant`), constants (BETA_STAR).
- `checkpoint.py` - append-only JSONL result store (`checkpoint.jsonl`).
  Every completed solve is written and fsynced immediately; runners skip
  already-checkpointed points, so all runs are interruption-safe and
  resumable.
- `validate2d.py` - gates V1 (uniform medium, alpha = 0) and V2 (step
  cylinder, alpha = +0.5, grid ladder h = 0.4/0.2/0.1).
- `phase1.py` - coarse full-domain runs: pair eta sweep at h = 0.2
  (models A/B, eta = 1e-2/3e-3/1e-3, both orientations), pair h = 0.1 at
  eta = 1e-2, far-separation check s = 10 (gate V4), single-tube ladder
  h = 0.2/0.1/0.05 at eta = 1e-2 (gate V3 against the 1D shooting
  result).
- `phase2.py` - quadrant-solver ladder: cross-validation points at
  h = 0.2/0.1, then pair + matching single-tube runs at h = 0.05 and
  0.025, eta = 1e-2 and 3e-3. These etas are the ones whose resonant-
  annulus width is resolvable at affordable h; smaller eta needs
  h <~ 0.007 (see results.md). The h = 0.0125 level is skipped (sparse LU
  projected at 25-40 GB against a 20 GB ceiling).
- `make_results.py` - regenerates `results.md` from `checkpoint.jsonl`:
  validation-gate tables, pair tables (full-domain and quadrant),
  same-grid ratios <alpha>/(2 x single), Richardson h->0 at fixed eta
  with a cautious Re-only eta->0 statement, negative-eps geometry of the
  two models, annulus-resolution criterion, caveats, timings. Also
  recomputes 1D shooting references at arbitrary eta by importing
  `../common.py`.
- `run_pair.py` - entry point chaining validate2d -> phase1 -> phase2 ->
  make_results. Safe to re-run at any time; completed solves are skipped.
- `checkpoint.jsonl` - one JSON record per completed solve: task, grid h,
  L, N, model, eta, s, orientation, Re/Im alpha, wall seconds, timestamp.
- `results.md` - generated report (do not edit by hand; edit
  make_results.py and regenerate).
- `phase1.log`, `phase2.log` - run logs.

## Usage

```
cd pair_cell
nice -n 10 python3.11 run_pair.py        # full pipeline, resumable
python3.11 make_results.py               # regenerate results.md only
```

Runs are strictly sequential by design (one sparse factorization at a
time); do not launch phases concurrently, as BLAS thread oversubscription
degrades all of them.
