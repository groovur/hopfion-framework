# The 2R0 linking threshold and the linked-lattice vacuum candidate

Status: step 1 of the pre-registered linked-lattice program (threshold
proof + equilibrium argument). Steps 2 (interaction sign) and 3
(periodic-cell response at spacing exactly 2R0) are separate.

## Theorem (linking bound)

If two circles C1, C2 of common radius R in R^3 have linking number
lk(C1, C2) != 0, then their centers satisfy |c1 - c2| < 2R, for every
relative orientation. The bound is sharp: the perpendicular Hopf
configuration is linked for every center distance d < 2R, and at
d = 2R the circles touch.

Proof. Let D1 be the flat disk bounded by C1. Since lk(C1, C2) != 0,
the curve C2 meets every Seifert surface of C1, in particular D1: there
is a point p in C2 with p in D1, so |p - c1| < R (if |p - c1| = R the
circles intersect, a degenerate configuration). Since p lies on the
circle C2, |p - c2| = R. By the triangle inequality
|c1 - c2| <= |c1 - p| + |p - c2| < R + R = 2R.

Sharpness. Take C1 in the xy-plane centered at the origin, C2 in the
xz-plane centered at (d, 0, 0). C2 crosses the xy-plane at x = d - R
and x = d + R; it threads C1 iff |d - R| < R, i.e. 0 < d < 2R. So
separations arbitrarily close to 2R are attained, only by
configurations approaching the perpendicular one. QED.

Two corollaries relevant to the condensate:

1. A network of pairwise-linked R0-circles cannot dilate beyond
   nearest-neighbor spacing 2R0 without unlinking events (tube
   crossings), which cost the crossing energy of the field.
2. Near the threshold the relative orientation of linked neighbors is
   forced toward perpendicular: tilted configurations unlink at
   strictly smaller separations.

## Equilibrium argument (conditional)

Assume (A1) the vacuum condensate is a network of mutually linked
Q_H = 2 hopfions, and (A2) the inter-hopfion interaction is repulsive
at separations of order 2R0 (step 2 computes this). Then the network is
pushed outward by (A2) and blocked at spacing 2R0 by the theorem: the
equilibrium spacing is the linking threshold itself,

    s_vacuum = 2R0 = 6  (condensate units, R0 = 3),

with locally perpendicular neighbor orientations (corollary 2). No
fitted parameter enters: the spacing is topology + repulsion.

Under this spacing, the combined polarization x Axiom-A2-vertex
requirement for the golden-angle step solves to spacing 6.027 (+0.45%
of 2R0) in the crude model (Maxwell-Garnett at n alpha ~ 0.18,
mean-field vertex, parallel-tube cross-section statistics); the
periodic-cell computation of step 3 replaces all three approximations
and is the decisive test.

## What is assumed, what is derived

Derived here: the threshold value 2R0 and the perpendicularity of
near-threshold neighbors (unconditional geometry); the equilibrium at
the threshold GIVEN (A1) + (A2).

Assumed: (A1) the linked topology of the vacuum network. Individual
hopfions are already stabilized by their own Hopf charge, so mutual
linking is not needed for single-soliton stability; (A1) is a statement
about which multi-hopfion configuration the condensate selects. A
formation argument (dense-phase entanglement locked in by the crossing
barrier as the condensate relaxes: the network can be diluted only to
the threshold, not through it) is plausible but not proved. No paper in
the series currently asserts or excludes (A1).

(A2) is computable and is step 2 (src_paper3/pair_interaction/). Note
the superposed-kernel computation of step 2 measures the metric
(field-overlap) interaction only; the crossing barrier that enforces
corollary 1 is a separate, topological effect and is not needed for the
equilibrium argument, only repulsion is.

## Pre-registered step 3

Periodic-cell (Bloch q -> 0) electromagnetic response of the linked
lattice at spacing exactly 2R0 (no fitted density), times the
Axiom-A2 vertex factor with exponent exactly 1, compared with
112.5/phi^10 = 0.914695. Sub-percent agreement identifies the Delta_1
carrier; a clear miss removes the linked-lattice candidate. The
orientation statistics of the cross-section must reflect perpendicular
neighbors (corollary 2), not parallel tubes.
