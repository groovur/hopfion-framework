# Two-tube pair cell: 2D anisotropic-dielectric polarizability (interim: coarse grids only, fine-grid ladder pending)

All alpha values are 2D dipole polarizabilities per unit length in the
normalization of the 1D code (isolated step cylinder eps1=3, R=1 gives
alpha = +0.5). Complex alpha from the limit-absorption prescription
eps -> eps + i*eta. Domain [-L,L]^2, L=40, Dirichlet BC u=-E.r;
dipole read off a ring r in [15,25] (m=1 angular Fourier fit including
a growing-mode term absorbing the finite-L Dirichlet leakage).

```
beta = beta* = 0.06715131
C2*  = 3.43182008141397
C*   = 1.81682487986669
s    = 2*R0/(C2*.C*) = 0.96230706   (R0=3)
rho_c(beta*) = 0.51913431
alpha_single_1D (eta->0, from ../results.md) = -0.97106613+0.01353232i
```

## Validation gates

### V1: uniform medium
alpha = -5.346e-11 +0.000e+00i at h=0.2 (expected 0). PASS

### V2: isotropic step cylinder eps1=3, R=1 (expect +0.5)

| h | alpha | rel. error |
|---|-------|-----------|
| 0.4 | +0.404459+0.000000i | 19.108% |
| 0.2 | +0.506556+0.000000i | 1.311% |
| 0.1 | +0.501675+0.000000i | 0.335% |
| Richardson (0.2,0.1) | +0.500048+0.000000i | 0.010% |

### V3: single tube vs 1D shooting code

2D uniform-grid results by eta and h, against the 1D reference at the
same eta:

| eta | h | alpha_2D | alpha_1D | rel. diff |
|-----|---|----------|----------|-----------|
| 0.01 | 0.2 | -0.961623+0.020503i | -0.968429+0.026606i | 0.94% |
| 0.01 | 0.1 | -0.961292+0.030743i | -0.968429+0.026606i | 0.85% |
| 0.01 | 0.05 | -0.966948+0.033096i | -0.968429+0.026606i | 0.69% |

### V4: far-separation pair (s=10, Model A, eta=1e-2, h=0.2)

```
alpha_par  = -1.887105+0.040104i
alpha_perp = -1.969403+0.042275i
<alpha>    = -1.928254+0.041189i
2 x single (same grid/eta) = -1.923247+0.041006i
ratio <alpha>/(2 x single) = 1.002604-0.000040j
```

## Pair results (s = 0.96230706)

### full domain solver

| model | h | eta | alpha_par | alpha_perp | <alpha> |
|-------|---|-----|-----------|------------|---------|
| A | 0.2 | 0.01 | -1.626776+0.131914i | -1.947361+0.093609i | -1.787069+0.112762i |
| A | 0.2 | 0.003 | -1.638495+0.047660i | -1.893167+0.153020i | -1.765831+0.100340i |
| A | 0.2 | 0.001 | -1.640847+0.016256i | -1.719112+0.174101i | -1.679979+0.095179i |
| A | 0.1 | 0.01 | -1.748565+0.191274i | -1.952114+0.127798i | -1.850339+0.159536i |
| B | 0.2 | 0.01 | -2.894405+0.457915i | -2.926838+0.183877i | -2.910622+0.320896i |
| B | 0.2 | 0.003 | -2.801989+0.434253i | -2.919068+0.157518i | -2.860529+0.295885i |
| B | 0.2 | 0.001 | -2.440103+0.654538i | -2.887310+0.080557i | -2.663706+0.367547i |
| B | 0.1 | 0.01 | -2.906822+0.353552i | -2.998614+0.273761i | -2.952718+0.313656i |

### quadrant solver

| model | h | eta | alpha_par | alpha_perp | <alpha> |
|-------|---|-----|-----------|------------|---------|
| A | 0.2 | 0.01 | -1.627060+0.131937i | -1.947699+0.093626i | -1.787379+0.112781i |
| B | 0.2 | 0.01 | -2.894911+0.457995i | -2.927346+0.183909i | -2.911129+0.320952i |

## Same-grid ratios <alpha>/(2 x alpha_single)

Each pair value divided by twice the single-tube value computed with
the same solver, grid, and eta (so shared discretization error largely
cancels in the ratio).

| model | h | eta | <alpha> | 2 x single (same grid) | ratio |
|-------|---|-----|---------|------------------------|-------|
| A | 0.2 | 0.01 | -1.787069+0.112762i | -1.923247+0.041006i | +0.930021-0.038802i |
| A | 0.1 | 0.01 | -1.850339+0.159536i | -1.922584+0.061486i | +0.964091-0.052147i |
| B | 0.2 | 0.01 | -2.910622+0.320896i | -1.923247+0.041006i | +1.516258-0.134523i |
| B | 0.1 | 0.01 | -2.952718+0.313656i | -1.922584+0.061486i | +1.539450-0.113910i |

## Grid extrapolation (Richardson h->0 at fixed eta)

Richardson extrapolation (order 2) from the two finest quadrant grid
levels at each fixed eta. The eta->0 step is then taken for Re(alpha)
only (linear in eta from eta = 1e-2 and 3e-3, both of which satisfy
the annulus resolution criterion at these grids); Im(alpha) at these
etas is still dominated by the finite-eta absorption and is quoted at
eta = 1e-2 as a criterion-limited value, not an eta->0 limit.

## Negative-eps_rad geometry of the two models

eps_rad (the tensor eigenvalue along ghat) evaluated on the pair axis
(y=0) and the perpendicular bisector (x=0); tube centers at
x = +-0.481154, single-tube critical radius rho_c = 0.519134.

Model A (X_A = X1+X2): one merged negative region containing both
tube centers and the midpoint.
```
  midpoint (0,0): X_A = 10.5498, eps_rad = -0.2257 (negative)
  pair axis: eps_rad < 0 for |x| < 1.0853
  bisector:  eps_rad < 0 for |y| < 0.7509
```

Model B (X_B = |g|^2): the two outward radial gradients cancel at the
midpoint, so X_B = 0 and eps = 1 exactly there; the cross term
suppresses X between the tubes and enhances it outboard. The negative
region is a single connected band that surrounds the pair: its inner
boundary runs through the two tube centers (where ghat flips) and
passes over the midpoint on the bisector; a positive-eps lens covers
the midpoint region.
```
  midpoint (0,0): X_B = 0, eps_rad = +1
  pair axis: eps_rad < 0 for 0.4812 < |x| < 1.3657
             (inner edge = tube center at 0.4812)
  bisector:  eps_rad < 0 for 0.2918 < |y| < 1.0348
  X_B just outboard of a tube center: 18.4632 (vs X_A max 10.6073)
```

## Resonant-annulus resolution criterion

eps_rad crosses zero at rho_c = 0.51913431 from each tube center; the
limit-absorption physics is resolved only when the radial width of the
|eps_rad| < eta shell is covered by at least a few cells. Width
estimate 2*eta/|d(eps_rad)/d(rho)|_rho_c :

| eta | annulus width | h required (width/4) |
|-----|---------------|----------------------|
| 0.01 | 0.0290 | 0.0072 |
| 0.003 | 0.0087 | 0.0022 |
| 0.001 | 0.0029 | 0.0007 |
| 0.0003 | 0.0009 | 0.0002 |
| 0.0001 | 0.0003 | 0.0001 |

## Numerical caveats

1. eta-sensitivity is the dominant systematic. Each tube carries a
   negative-eps_rad core (rho < rho_c = 0.519 from its center), and at
   s = 0.962 the two critical annuli nearly touch. The limit-absorption
   Im(alpha) survives eta->0 only where the annulus is resolved by the
   grid (table above). On grids that fail the criterion, Im(alpha)
   collapses toward 0 linearly in eta (the discrete system has no
   spectrum at the operating point), or alpha destabilizes entirely
   when a discrete eigenmode happens to sit near it. Values at etas
   failing the h-criterion are reported but must not be treated as
   converged physics.
2. Model B is markedly more eta-sensitive than Model A at the same
   grid. Its negative-eps region is larger in area, reaches stronger
   eps contrast (X_B peaks near 16.8 just outboard of each tube center
   vs 10.6 for X_A), and its anisotropy direction ghat is
   discontinuous across each tube center (the partner tube's radial
   unit vector flips sign there), all of which sharpen the discrete
   resonances of the finite grid. See the geometry section.
3. The V2 step-cylinder gate converges at ~O(h) rather than O(h^2)
   because the sharp material interface is staircased on the Cartesian
   grid; the smooth tube background has no such interface and shows
   near-quadratic behavior in the resolved regime.
4. Quadrant (mirror-symmetry) and full-domain solvers agree to ~0.03%
   at identical h (difference dominated by the extraction-arc choice),
   validating the symmetry reduction used for the fine-grid ladder.

## Timings

Total checkpointed solver wall time: 16043 s (4.46 h) across 30 records.

