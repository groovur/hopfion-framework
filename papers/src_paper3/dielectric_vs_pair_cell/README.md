# paper3_support

Support material for the two-channel Delta_1 candidate discussed in
`papers/main_paper3.tex` (Remark P3:rem:delta1_candidate).

## two_channel_bracket.py

Computes `aInv_obs/aInv_geo = eps_eff * <1/(1+beta X)>^{-1}` for a
contact-packed condensate of Q_H=2 two-tube (T(2,2)) hopfions and compares
the result with the golden/geo target `112.5/phi^10 = 0.914695`. All inputs
(b*, C*^2, contact spacing, pair separation, `<alpha_pair>`) are fixed by the
framework; there is no free density parameter. Running it with Python 3.11
reproduces the table in the tex comment block:

| lattice    | model | eps_eff | vertex  | product | dev     |
|------------|-------|---------|---------|---------|---------|
| square     | A     | 0.82576 | 1.17679 | 0.97178 | +6.24%  |
| square     | B     | 0.73568 | 1.17679 | 0.86571 | -5.36%  |
| triangular | A     | 0.80148 | 1.20988 | 0.96974 | +6.02%  |
| triangular | B     | 0.70090 | 1.20988 | 0.84797 | -7.29%  |

## dielectric_cell/ and pair_cell/

Snapshot copies of the solver suites whose originals live in
`papers/src_paper3/dielectric_cell/` (and its `pair_cell/` subdirectory).
The `pair_cell/checkpoint.jsonl` and `pair_cell/results.md` snapshots here
were taken while the fine-grid ladder was still appending to the original
checkpoint file; they reflect the solver state at copy time, not the final
converged run.

## Dependency chain

1. 1D shooting solver (`dielectric_cell/`) -> single-tube alpha
2. 2D pair solver (`pair_cell/`) -> `<alpha_pair>`
3. `<alpha_pair>` -> `two_channel_bracket.py`
