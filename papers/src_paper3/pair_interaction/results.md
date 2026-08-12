# Metric interaction energy between two Q_H=2 hopfion tubes vs separation d

Computation: `interaction.py` in this directory. Superposed-kernel interaction
only (kappa = kappa_1 + kappa_2), no field relaxation. R0 = 3, mustar =
3 - phi = 1.3819660113. Two configurations (perp, coax), two beta values
(0.452, 0.06715131). d grid: 3.5, 4, 4.5, 5, 5.5, 5.8, 6.0, 6.2, 6.5, 7, 8, 10.

## Numerics and cost

Midpoint-rule 3D quadrature, box [-12, d+12] along the separation axis and
[-12, 12] in the two transverse axes, looped in slabs along the elongated
axis to bound memory. Grid ladder h = 0.25, 0.125, 0.0625 (up to ~8e7 points
at h = 0.0625). Coarse pass (h = 0.25, all 24 d/config combinations) took
0.78 s; projected total for the full three-level ladder was ~57 s; actual
total wall time was 61 s. All three levels were run for the full grid; no
subsetting was needed.

## Convergence

`Delta_J2a` (linear functional, exact cancellation expected):

| config | d | Delta_J2a (h=0.0625) |
|---|---|---|
| perp | 3.5 | -3.4e-13 |
| perp | 6.0 | 1.1e-13 |
| perp | 10.0 | 5.7e-13 |
| coax | 3.5 | -9.1e-13 |
| coax | 6.0 | 1.7e-13 |
| coax | 10.0 | 8.5e-13 |

Machine-precision zero, as required by linearity of J2a — confirms the
sum/difference bookkeeping in the code is correct, independent of quadrature.

`Delta_J4` across the grid ladder (perp, values at h = 0.25 / 0.125 / 0.0625):

| d | perp | coax |
|---|---|---|
| 3.5 | 43.2883 / 43.2872 / 43.2870 | 23.6945 / 23.6944 / 23.6944 |
| 6.0 | 94.4985 / 94.4985 / 94.4985 | 3.17815 / 3.17815 / 3.17815 |
| 10.0 | 3.62519 / 3.62519 / 3.62519 | 0.410262 / 0.410261 / 0.410261 |

`E_pair(d; beta=0.452)` across the grid ladder:

| d | perp | coax |
|---|---|---|
| 3.5 | 1130076.46 / 1130072.06 / 1130071.26 | 1105187.30 / 1105184.28 / 1105183.65 |
| 6.0 | 1215719.97 / 1215716.84 / 1215716.19 | 1083602.18 / 1083599.30 / 1083598.69 |
| 10.0 | 1085391.08 / 1085388.37 / 1085387.81 | 1081345.30 / 1081342.45 / 1081341.85 |

Relative change from h=0.125 to h=0.0625 is ~1e-6, so h=0.0625 results below
are converged well past the precision needed for sign/shape conclusions.

Single-circle thin-tube check: analytic thin-tube estimate
J2a_thin = 16*pi^2*R0 = 473.741. Numeric single-circle J2a at h=0.0625:
477.174 (perp box) / 476.919 (coax box) — agreement at the 0.7-0.9% level,
consistent with the expected curvature correction to the thin-tube
approximation (the box shape differs slightly between configs, which is why
the two numbers differ from each other at the ~0.05% level; both are well
converged in h).

Far-separation check at d = 10 (largest sampled separation): Delta_J4 =
3.625 (perp), 0.410 (coax); Delta_Jfb(beta=0.452) = -3.073 (perp), -0.382
(coax). These are small compared to their peak values near d~5-6 (94.5 and
23.7 for Delta_J4) but not yet negligible — the kernel tail (~1/rho^4) makes
the cross-term integral decay slowly, so full vanishing would require d
well beyond 10. The qualitative point holds: all Delta F magnitudes decrease
monotonically toward 0 as d increases across the sampled range.

## Delta_Jfb(d) and Delta_J4(d) tables (h = 0.0625)

### perp

| d | Delta_Jfb(beta=0.452) | Delta_Jfb(beta=0.06715131) | Delta_J4 |
|---|---|---|---|
| 3.5 | -26.7399 | -8.8715 | 43.2870 |
| 4.0 | -26.9131 | -9.7393 | 49.2994 |
| 4.5 | -27.8784 | -11.4106 | 61.3217 |
| 5.0 | -28.8282 | -13.4330 | 78.0245 |
| 5.5 | -28.8959 | -14.8754 | 92.9554 |
| 5.8 | -28.2194 | -14.9898 | 96.1321 |
| 6.0 | -27.4236 | -14.6567 | 94.4985 |
| 6.2 | -26.3586 | -14.0001 | 89.7982 |
| 6.5 | -24.3203 | -12.5145 | 78.2250 |
| 7.0 | -20.1485 | -9.3653 | 54.4000 |
| 8.0 | -11.7334 | -4.1505 | 20.9265 |
| 10.0 | -3.0726 | -0.7924 | 3.6252 |

### coax

| d | Delta_Jfb(beta=0.452) | Delta_Jfb(beta=0.06715131) | Delta_J4 |
|---|---|---|---|
| 3.5 | -17.3273 | -5.0378 | 23.6944 |
| 4.0 | -11.7408 | -3.1974 | 14.8026 |
| 4.5 | -8.0466 | -2.0953 | 9.6044 |
| 5.0 | -5.6057 | -1.4158 | 6.4497 |
| 5.5 | -3.9784 | -0.9842 | 4.4663 |
| 5.8 | -3.2684 | -0.8013 | 3.6304 |
| 6.0 | -2.8778 | -0.7020 | 3.1782 |
| 6.2 | -2.5414 | -0.6174 | 2.7929 |
| 6.5 | -2.1206 | -0.5124 | 2.3166 |
| 7.0 | -1.5902 | -0.3818 | 1.7249 |
| 8.0 | -0.9379 | -0.2237 | 1.0100 |
| 10.0 | -0.3816 | -0.0908 | 0.4103 |

Delta_Jfb < 0 everywhere for both configurations, both betas: Jfb is a
concave, saturating function of kappa, so Jfb[k1+k2] < Jfb[k1] + Jfb[k2]
whenever the two kernels have overlapping support, regardless of separation.
This sign is a generic consequence of the functional's concavity, not by
itself the sign of the physical interaction. Delta_J4 > 0 everywhere, also
generic: Delta_J4 = (1/2) Integral kappa_1 kappa_2 dV exactly (the cross
term of the square), which is positive-definite for any nonzero overlap.

## E_pair(d) tables (h = 0.0625)

E_pair(d) = Kfb[kappa_1+kappa_2; beta] * J4[kappa_1+kappa_2], the physically
clean quantity (E evaluated directly on the superposed field; not a
difference against the isolated single-hopfion energies, since the product
form Kfb*J4 is not additive even when the two kernels have disjoint
support).

### perp

| d | E_pair(beta=0.452) | E_pair(beta=0.06715131) |
|---|---|---|
| 3.5 | 1130071.26 | 1404393.41 |
| 4.0 | 1140111.59 | 1416223.71 |
| 4.5 | 1159439.23 | 1439749.83 |
| 5.0 | 1186566.20 | 1472580.91 |
| 5.5 | 1211565.39 | 1502220.83 |
| 5.8 | 1217620.54 | 1508754.61 |
| 6.0 | 1215716.19 | 1505732.39 |
| 6.2 | 1208937.54 | 1496669.41 |
| 6.5 | 1191598.70 | 1474172.84 |
| 7.0 | 1155634.46 | 1427755.52 |
| 8.0 | 1106958.08 | 1362916.21 |
| 10.0 | 1085387.81 | 1329904.53 |

### coax

| d | E_pair(beta=0.452) | E_pair(beta=0.06715131) |
|---|---|---|
| 3.5 | 1105183.65 | 1366490.04 |
| 4.0 | 1095263.60 | 1349677.81 |
| 4.5 | 1089805.01 | 1339876.34 |
| 5.0 | 1086662.96 | 1333950.76 |
| 5.5 | 1084775.79 | 1330243.81 |
| 5.8 | 1084010.03 | 1328692.76 |
| 6.0 | 1083598.69 | 1327851.61 |
| 6.2 | 1083253.03 | 1327137.04 |
| 6.5 | 1082840.30 | 1326264.42 |
| 7.0 | 1082338.18 | 1325185.20 |
| 8.0 | 1081763.67 | 1323904.53 |
| 10.0 | 1081341.85 | 1322883.35 |

Central-difference estimate of dE_pair/dd at d = 6.0, using the samples at
d = 5.8 and 6.2 (spacing 0.4):

| config | beta | dE_pair/dd at d=6 |
|---|---|---|
| perp | 0.452 | -21707 |
| perp | 0.06715131 | -30215 |
| coax | 0.452 | -1892 |
| coax | 0.06715131 | -3890 |

All four are negative: at d = 6, E_pair is decreasing in d in every
configuration/beta combination — the interaction is repulsive there.

## A geometric feature specific to the perp configuration (important caveat)

For the perp geometry, the minimum Euclidean distance between the two
defining circles (not their kernel-broadened tubes, the bare curves) is
|6 - d| for d in this range (verified numerically by direct grid search over
both circle parameters). At d = 6.0 = 2*R0 the two circles are exactly
tangent (zero separation) — a geometric coincidence of this parametrization
given R0 = 3, not a topological effect. This produces a real, non-monotonic
bump in Delta_J4(d) and Delta_Jfb(d) (and hence in E_pair(d)) centered near
d approx 5.7-5.9, i.e. very close to, but on the near side of, d = 6:
E_pair(perp, beta=0.452) rises from d=3.5 to a local maximum at d=5.8
(1217620.5) then falls through d=6.0 (1215716.2), d=6.2 (1208937.5), and
onward monotonically down to the d=10 floor. The same shape appears for
beta=0.06715131. In other words, for perp: the interaction is attractive
(dE_pair/dd > 0, energy lowered by decreasing d) for d roughly between 3.5
and 5.7, and repulsive (dE_pair/dd < 0) from about d=5.8 out through d=10.
By d = 6.0 itself the derivative is already negative in every case checked
above, so "repulsive at d near 6" holds, but the extremum sits close enough
to d=6 that it should be read as a genuine, geometry-driven feature of the
"circle-touches-circle" parametrization at 2*R0 = 6, not as evidence of any
distinguished physical role for d=6 as such. This is flagged explicitly
because it is the opposite of a null result: there IS a local feature near
d=6 in the perp configuration, and it traces to bare-curve tangency, not to
the field-overlap physics being blind to topology.

For the coax configuration, the minimum circle-circle distance is exactly d
(no tangency anywhere in this range), and correspondingly E_pair(d) is
monotonically decreasing over the entire sampled range 3.5 <= d <= 10, for
both betas, with no local extremum, no kink, and no special behavior at
d = 6 — this configuration is a clean confirmation that in the absence of
a coincidental curve-tangency, the superposed-kernel energy has no feature
at d = 6.

## Plain answer

In both configurations and both beta values, at d near 6 the metric
(field-overlap) interaction is repulsive: E_pair(d) is decreasing in d
there, so increasing separation lowers the superposed-kernel energy. For
coax, this repulsion holds smoothly across the entire range 3.5 <= d <= 10
with no extremum and no feature at d=6, exactly as expected for a kernel
overlap that is blind to topology. For perp, the same repulsive conclusion
holds at d=6 and beyond, but there is a local maximum of E_pair(d) near
d approx 5.8 (just below d=6) below which the interaction becomes
attractive as d decreases further toward 3.5; that maximum is traced
directly to the two defining circles of the perp construction becoming
exactly tangent at d = 2*R0 = 6, a coordinate/geometric coincidence of
this specific setup (R0=3) rather than a topological signature. No interior
minimum of E_pair(d) exists over the sampled range 3.5 <= d <= 10 for
either configuration; E_pair(d) is bounded below by an asymptotic floor
approached monotonically as d increases past the (perp) local maximum or
(coax) from d=3.5 onward, and that floor is not yet fully reached at d=10
(still declining slowly there, consistent with the slow 1/rho^4 kernel
tail). So: repulsive near and above d=6 in every case tested; the "no
feature at d=6" expectation is confirmed for coax and contradicted for
perp, where the near-d=6 feature is real but attributable to bare-curve
tangency, not topology.
