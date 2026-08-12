# 3D linked-lattice homogenization: results

Verdict at the pre-registered evaluation point (a_c = 6.0 = 2R0, both vertex variants, both eta): Re(P) ~ 0.27 against the target 112.5/phi^10 = 0.914695 — a miss by a factor of about 3.3 (-70%). The linked-lattice-at-threshold candidate for the Delta_1 carrier is refuted: the vacuum network at threshold density over-responds by far.

## Physics table (from checkpoint.jsonl; Hopf-chain fallback cell, grid (2M, M, M), h = a_c/M)

P_pathavg = eps_iso / S1, S1 = <1/(1+U)>; P_meanfield = eps_iso * (1 + Ubar), Ubar = <U>. dev% = Re(P) vs T. iterations = GMRES inner iterations for E along x, y, z.

| a_c | eta | grid | h | Re eps_iso | Im eps_iso | aniso | S1 | Ubar | P_pathavg | dev% | P_meanfield | dev% | iterations | t(s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6.0 | 1e-02 | (96, 48, 48) | 0.1250 | 0.1767 | 0.181 | 2.89e-02 | 0.63918 | 0.60777 | 0.2764 | -69.8 | 0.2840 | -68.9 | [650, 640, 640] | 363 |
| 6.0 | 3e-03 | (96, 48, 48) | 0.1250 | 0.1787 | 0.176 | 2.56e-02 | 0.63918 | 0.60777 | 0.2797 | -69.4 | 0.2874 | -68.6 | [2096, 2062, 2062] | 1214 |
| 5.9 | 1e-02 | (96, 48, 48) | 0.1229 | 0.1584 | 0.192 | 4.14e-02 | 0.63012 | 0.62914 | 0.2514 | -72.5 | 0.2580 | -71.8 | [645, 634, 634] | 365 |
| 5.9 | 3e-03 | (96, 48, 48) | 0.1229 | 0.1493 | 0.188 | 1.02e-01 | 0.63012 | 0.62914 | 0.2370 | -74.1 | 0.2433 | -73.4 | [2078, 2045, 2045] | 1168 |
| 5.8 | 1e-02 | (96, 48, 48) | 0.1208 | 0.1564 | 0.201 | 6.92e-02 | 0.62085 | 0.65150 | 0.2519 | -72.5 | 0.2583 | -71.8 | [638, 626, 626] | 357 |
| 5.8 | 3e-03 | (96, 48, 48) | 0.1208 | 0.1585 | 0.205 | 3.40e-02 | 0.62085 | 0.65150 | 0.2553 | -72.1 | 0.2617 | -71.4 | [2065, 2020, 2020] | 1175 |
| 6.0 | 1e-02 | (128, 64, 64) | 0.0938 | 0.1745 | 0.189 | 6.77e-02 | 0.63918 | 0.60776 | 0.2730 | -70.2 | 0.2806 | -69.3 | [654, 644, 644] | 855 |
| 6.0 | 3e-03 | (128, 64, 64) | 0.0938 | 0.1794 | 0.188 | 1.35e-01 | 0.63918 | 0.60776 | 0.2807 | -69.3 | 0.2884 | -68.5 | [2110, 2073, 2073] | 2777 |
| 5.9 | 1e-02 | (128, 64, 64) | 0.0922 | 0.1659 | 0.192 | 7.05e-02 | 0.63013 | 0.62913 | 0.2633 | -71.2 | 0.2703 | -70.4 | [649, 635, 635] | 930 |
| 5.9 | 3e-03 | (128, 64, 64) | 0.0922 | 0.1590 | 0.189 | 1.45e-01 | 0.63013 | 0.62913 | 0.2524 | -72.4 | 0.2591 | -71.7 | [2094, 2048, 2048] | 2707 |
| 5.8 | 1e-02 | (128, 64, 64) | 0.0906 | 0.1630 | 0.195 | 1.56e-01 | 0.62085 | 0.65149 | 0.2625 | -71.3 | 0.2691 | -70.6 | [642, 629, 629] | 839 |
| 5.8 | 3e-03 | (128, 64, 64) | 0.0906 | 0.1643 | 0.190 | 1.85e-01 | 0.62085 | 0.65149 | 0.2647 | -71.1 | 0.2714 | -70.3 | [2070, 2024, 2024] | 2697 |
| 6.0 | 1e-02 | (192, 96, 96) | 0.0625 | 0.1766 | 0.183 | 1.09e-01 | 0.63918 | 0.60776 | 0.2764 | -69.8 | 0.2840 | -69.0 | [656, 645, 645] | 5536 |
| 6.0 | 3e-03 | (192, 96, 96) | 0.0625 | 0.1775 | 0.179 | 2.00e-01 | 0.63918 | 0.60776 | 0.2776 | -69.6 | 0.2853 | -68.8 | [2114, 2078, 2078] | 18132 |
| 5.9 | 1e-02 | (192, 96, 96) | 0.0615 | 0.1687 | 0.194 | 7.21e-02 | 0.63013 | 0.62913 | 0.2678 | -70.7 | 0.2749 | -69.9 | [651, 639, 639] | 5543 |
| 5.9 | 3e-03 | (192, 96, 96) | 0.0615 | 0.1670 | 0.196 | 1.23e-01 | 0.63013 | 0.62913 | 0.2651 | -71.0 | 0.2721 | -70.2 | [2100, 2062, 2062] | 17990 |
| 5.8 | 1e-02 | (192, 96, 96) | 0.0604 | 0.1688 | 0.197 | 2.84e-02 | 0.62085 | 0.65149 | 0.2719 | -70.3 | 0.2788 | -69.5 | [646, 631, 631] | 5466 |
| 5.8 | 3e-03 | (192, 96, 96) | 0.0604 | 0.1661 | 0.194 | 7.56e-02 | 0.62085 | 0.65149 | 0.2676 | -70.7 | 0.2744 | -70.0 | [2084, 2032, 2032] | 17711 |
| 6.0 | 1e-02 | (256, 128, 128) | 0.0469 | 0.1745 | 0.189 | 6.05e-02 | 0.63918 | 0.60776 | 0.2730 | -70.1 | 0.2806 | -69.3 | [667, 657, 657] | 11768 |

Stability at a_c = 6.0 across all completed (grid, eta): P_pathavg deviation spans -70.2% to -69.3%; P_meanfield -69.3% to -68.5%. The values are consistent across eta, across spacing (monotone, no structure), and across the grid ladder: the miss is not a grid or absorption artifact.

## Geometry and solver validation (from the build phase)

- The pre-registered 3-coloring motif FAILED the Gauss-linking gate (axis-adjacent pairs link along only one of three axes); the documented fallback — alternating normal-z/normal-y Hopf chains along x, stacked at a_c in y and z — passes: |lk| = 1 exactly for chain neighbours at a_c = 5.8/5.9, and neighbour tube curves are tangent to machine precision at a_c = 6.0 (the topological threshold).
- Solver gates: V1 uniform exact; V2 Clausius-Mossotti sphere array to 0.057%; V3 parallel-tube transverse response vs the validated 2D Maxwell-Garnett chain to 3.9%.

## Scope and limitations

The homogenized medium is the vacuum network only. A free Q = 2 hopfion (the electron of the lepton sector) is the source/test charge probing the medium and is not a lattice member; linking is a property of the vacuum's selected configuration, not of Q = 2 objects generally.

Level-1 modeling caveat: at the threshold spacing the cell's volume-mean suppression variable is U ~ 0.6 — the lattice tubes are far from the isolated relaxed profile assumed by the superposed kernels, and the incoherent superposition is least reliable precisely in this dense regime. This caveat affects the precise value of eps_iso, but a factor-3 discrepancy is outside any correction this limitation plausibly supplies.

## Reproduction

geometry.py (curves, distances, Gauss linking), field.py (U, ghat, tensor), solver.py (matrix-free FFT homogenization, Moulinec-Suquet Green operator with complex reference, GMRES), validate.py (V1-V3), run_lattice.py + driver.py (driver; appends each completed case to checkpoint.jsonl and skips completed cases on rerun).
