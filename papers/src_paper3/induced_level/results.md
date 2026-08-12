# Twisted-sector quantum dimensions and the induced U(1) level

Condensate CFT: SU(2)_3 (untwisted primaries h in {0, 3/20, 2/5, 3/4}, quantum
dims {1, phi, phi, 1}) with a 2I-orbifold twisted sector labeled by the E8
Coxeter exponents m in {1,7,11,13,17,19,23,29}, h_m = m/90, c = 9/5.
D_untw^2 = 5 + sqrt(5).

All values below are exact: rationals or elements of Q(sqrt5), except where
noted (candidate A3's twisted-sector quantities are honest algebraic numbers
of degree 4 over Q, certified by minimal polynomial, not elements of
Q(sqrt5)). Numerics (mpmath, 50-60 digit precision) are used only to
cross-check the exact symbolic results; every "no exact match" verdict below
is a symbolic zero-test (`sympy.simplify(difference) == 0`), not a numeric
tolerance call, and every gate/target failure is confirmed by a gap of order
1-400 at 60-digit precision -- far outside any possibility of being a
precision artifact.

## Question A: twisted-sector quantum dimensions d_m

Three pre-registered candidates were tested; no others.

- **A1**: d_m = 1 for all m (abelian twisted sector).
- **A2**: d_m = finite-E8 Dynkin marks {2,3,4,5,6,4,2,3}, matched
  positionally to the exponent list (1,7,11,13,17,19,23,29) in the order
  both are given in the task statement:
  m=1->2, m=7->3, m=11->4, m=13->5, m=17->6, m=19->4, m=23->2, m=29->3.
  (These marks are the dimensions of the 8 nontrivial irreps of the binary
  icosahedral group 2I, i.e. the finite nodes of the affine E8 Dynkin
  diagram under McKay correspondence; they sum to h-1=29 with E8 Coxeter
  number h=30.)
- **A3**: d_m = sin(pi m/30)/sin(pi/30) (E8-Coxeter-root-of-unity quantum
  dims, q = exp(i*pi/30)).

### A3 exact values (not in Q(sqrt5))

| m | d_m minimal polynomial (degree) | d_m (25 dp) |
|---|---|---|
| 1, 29 | x - 1 (deg 1) | 1.0000000000000000000000000 |
| 7, 23 | x^4 - 6x^3 - 4x^2 + 9x + 1 (deg 4) | 6.401420105502707915065975 |
| 11, 19 | x^4 - 11x^3 + 21x^2 - 11x + 1 (deg 4) | 8.739681318220424342718522 |
| 13, 17 | x^4 - 9x^3 - 4x^2 + 6x + 1 (deg 4) | 9.357715306970319190923109 |

Each nontrivial A3 value is a genuine quartic algebraic number, not
expressible in Q(sqrt5) alone (30 is a constructible-angle denominator so
d_m is a real square-root radical, but its minimal polynomial over Q has
degree 4, not 2).

### D^2 = D_untw^2 + sum d_m^2

| Candidate | D^2 (exact) | D^2 (numeric, 40 dp) |
|---|---|---|
| A1 | 13 + sqrt(5) | 15.2360679774997896964091736687312762354 |
| A2 | 124 + sqrt(5) | 126.2360679774997896964091736687312762354 |
| A3 | degree-4 algebraic (twisted part has minimal polynomial x^4-432x^3+8424x^2-52488x+104976) | 419.0901573324838085054635290782549334148 |

### Gauss sums Theta_pm = sum_a d_a^2 exp(+-2*pi*i*h_a) over ALL 12 primaries, and modularity gate

Gate: |Theta_+| = D (equivalently Theta_+ * Theta_- = D^2), and
Theta_+/|Theta_+| = exp(2*pi*i*c_eff/8).

| Candidate | \|Theta_+\| | D | \|Theta_+\| - D | Gate \|Theta_+\|=D | c_eff mod 8 (raw phase, only meaningful if gate holds) |
|---|---|---|---|---|---|
| A1 | 9.279078765960998310 | 3.903340617663258669 | 5.375738148297739641 | **FAIL** | 1.465849873412633592 |
| A2 | 112.0432701505338971 | 11.23548254315317735 | 100.8077876073807198 | **FAIL** | 1.378690286762413075 |
| A3 | 393.4093478883388624 | 20.47169160896294394 | 372.9376562793759184 | **FAIL** | 1.336453274173248493 |

All three candidates fail the modularity gate, by large margins (5.4, 101,
373 respectively, at 60-digit precision -- not near-misses). None of the
raw c_eff readouts is a recognizably simple rational (e.g. none equals
c mod 8 = 9/5 = 1.8), consistent with the gate already failing: c_eff is
only physically meaningful for a genuine modular tensor category, which
none of these three candidates form. In particular the physical target
c mod 8 = 9/5 is not reproduced by any candidate.

**Per the task's own note**: A1 (the abelian/Z_360 reading) is retained as
a live candidate for Question B despite the gate failure, since an abelian
Z_360 phase lattice need not be the modular closure of an honest MTC.

## Question B: induced U(1) level

Framework charge rule: q_a = h_a (native), and variant q_a = h_a - c/24 = h_a
- 3/40. S2 = sum_a d_a^2 q_a^2 (raw), and S2/D^2 (normalized). Targets
(pre-registered, exact-or-nothing): 360, 360/phi^2 = 540 - 180*sqrt(5)
(numeric 137.50776405...), 90, 40.

| Candidate | charge convention | S2 (raw, exact) | S2 numeric | S2/D^2 (exact) | S2/D^2 numeric | exact target match |
|---|---|---|---|---|---|---|
| A1 | q=h | 73*sqrt(5)/800 + 73069/64800 | 1.331649227638 | 19*sqrt(5)/53136 + 230083/2656800 | 0.087401108318 | none |
| A1 | q=h-c/24 | 89*sqrt(5)/1600 + 99077/129600 | 0.888864305940 | 312989/5313600 - 67*sqrt(5)/265680 | 0.058339481502 | none |
| A2 | q=h | 73*sqrt(5)/800 + 313957/64800 | 5.049056635046 | 83851*sqrt(5)/199208160 + 38901103/996040800 | 0.039996941571 | none |
| A2 | q=h-c/24 | 89*sqrt(5)/1600 + 2909/1296 | 2.368980046681 | 75377*sqrt(5)/249010200 + 7207111/398416320 | 0.018766269297 | none |
| A3 | q=h | (degree-4 algebraic) | 13.564875980157 | (degree-4 algebraic) | 0.032367441093 | none |
| A3 | q=h-c/24 | (degree-4 algebraic) | 5.297418077206 | (degree-4 algebraic) | 0.012640282728 | none |

No exact match to any pre-registered target for any (candidate, charge
convention, normalization) triple. The raw and normalized values are not
even close to 360, 137.5..., 90, or 40 in most cases (closest raw value is
A3's 13.56, still off by an order of magnitude from the smallest target 40).

### Full Z_360 abelian lattice (maximal-refinement reference)

All 360 phases exp(2*pi*i*n/360), d_n = 1, q_n = n/360, n = 0..359.

- S2 = sum_{n=0}^{359} (n/360)^2 = 258121/2160 = 119.50046296296296296 (exact rational)
- D^2 (this lattice) = 360 (exact, since all d_n = 1)
- S2/360 = S2/D^2 = 258121/777600 = 0.33194573045267489712 (exact rational; these two normalizations coincide here because D^2 = 360 exactly)

No exact match to any target (raw or normalized).

## Question C

alpha_inv = 360/phi^2 = 540 - 180*sqrt(5) (exact, numeric 137.50776405...).
Checking every (candidate, charge convention, normalization) triple computed
in Question B (six spectrum-derived triples plus the full-Z_360-lattice
reference, in both raw and normalized form) for an exact k_ind = 360:
**none** of them equals 360 exactly. The closest raw value obtained anywhere
is A3's 13.56 (q=h) -- off from 360 by more than an order of magnitude, and
every normalized value is well below 1.

## Plain-language summary

None of the three pre-registered twisted-sector quantum-dimension
assignments (uniform d_m=1, finite-E8 Dynkin marks, or E8-Coxeter-root-of-
unity dims sin(pi m/30)/sin(pi/30)) yields a modular tensor category: the
Gauss-sum gate |Theta_+| = D fails for all three, by large, precisely
computed margins, not near-misses. Correspondingly, no induced U(1) level
computed from any of these candidates -- under either the native charge
rule q_a = h_a or the phase-exponent variant q_a = h_a - c/24, raw or
D^2-normalized -- reproduces any of the four pre-registered targets
(360, 360/phi^2, 90, 40) exactly, nor does the maximal-refinement Z_360
abelian-lattice reference. In particular, alpha_inv = 360/phi^2 does not
admit a reading as k_ind * <U> with k_ind = 360 emerging from any spectrum
constructed here: that decomposition, if wanted, has to be imposed by hand
(k_ind := 360, <U> := 1/phi^2) rather than derived from summing over the
stated primary content. The one genuinely new structural fact surfaced here
is that candidate A3's twisted-sector quantum dimensions are not elements
of Q(sqrt5) but honest degree-4 algebraic numbers (minimal polynomials
given above), so any downstream use of A3 needs to track a strictly larger
number field than the rest of the framework.
