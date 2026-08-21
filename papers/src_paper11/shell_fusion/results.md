# Open-shell multiplet structure vs. SU(2)_3 WZW fusion data

`shell.py` builds the p^2 and d^2 two-electron electrostatic problem from
scratch (single-particle m_l, m_s states, antisymmetrized Slater
determinants, exact Condon-Shortley coefficients from Wigner 3j symbols)
and diagonalizes the Coulomb interaction exactly in sympy. No floats, no
tolerances anywhere; every claim below is either an exact symbolic
equality (verified in code) or a proof (stated as such).

Independent validation note: the first implementation of the exchange
integral used the naive index order c^k(m_a,m_b)*c^k(m_b,m_a) for the
off-diagonal ("exchange") term. This silently produced a NEGATIVE
exchange contribution for some (m,m') pairs (since c^k(m,m') =
(-1)^(m-m') c^k(m',m) for l=1, an antisymmetry, not a symmetry), and gave
term degeneracies {5,9,1} instead of the required {9,5,1} -- 3P and 1D
swapped. This was caught only by the hard requirement that the exact
diagonalization reproduce the known degeneracy multiset, not by any
individual number looking wrong. The fix: the correct formula pairs
bra/ket indices as c^k(m_a,m_c)*c^k(m_d,m_b) (second factor index order
reversed relative to the first), which reduces the exchange integral to
sum_k [c^k(m_a,m_b)]^2 F^k -- manifestly non-negative, as a real Coulomb
kernel requires. This is documented in `shell.py` at `raw_integral()`.

## Part A -- validated multiplet tables

### p^2 (equivalently p^4 by particle-hole symmetry)

Basis: 6 single-particle states (m_l in {-1,0,1}) x (m_s = +-1/2),
15 antisymmetrized two-electron states = C(6,2). Diagonalized exactly in
each (M_L,M_S) block; term labels assigned by the exact Racah-theorem
peeling algorithm (a term (L,S) must have identical energy at every one
of its (2L+1)(2S+1) microstates, verified against the exact eigenvalues,
not assumed).

| Term | Degeneracy (2L+1)(2S+1) | Energy (raw F0,F2) | Energy (scaled, F2_scaled=F2_raw/25) |
|------|------|------|------|
| 3P   | 9  | F0 - F2/5   | F0 - 5 F2_scaled  |
| 1D   | 5  | F0 + F2/25  | F0 + 1 F2_scaled  |
| 1S   | 1  | F0 + 2F2/5  | F0 + 10 F2_scaled |

Total degeneracy: 9+5+1 = 15 = C(6,2). EXACT match to the textbook
result E = F0 + {-5,+1,+10} F2. Ground term (lowest energy, F2>0):
3P -- matches Hund's rule. Both facts asserted and checked in `shell.py`
(`run_p2`).

### d^2

Basis: 10 single-particle states (m_l in {-2,...,2}) x (m_s=+-1/2),
45 antisymmetrized states = C(10,2). Same exact diagonalization +
peeling procedure, extended to k=0,2,4 (F0,F2,F4).

| Term | Degeneracy | Energy (raw F0,F2,F4) | Energy (Racah A,B,C) |
|------|------|------|------|
| 3F | 21 | F0 - 8F2/49 - F4/49  | A - 8B |
| 3P | 9  | F0 + F2/7 - 4F4/21   | A + 7B |
| 1G | 9  | F0 + 4F2/49 + F4/441 | A + 4B + 2C |
| 1D | 5  | F0 - 3F2/49 + 4F4/49 | A - 3B + 2C |
| 1S | 1  | F0 + 2F2/7 + 2F4/7   | A + 14B + 7C |

(Racah scaling: F2_raw = 49 F2_scaled, F4_raw = 441 F4_scaled,
A = F0 - 49 F4_scaled, B = F2_scaled - 5 F4_scaled, C = 35 F4_scaled.)

Total degeneracy: 21+9+9+5+1 = 45 = C(10,2). All five term energies match
the textbook Racah A,B,C formulas EXACTLY (verified symbolically in
`shell.py`, `run_d2`; this is a 5-way check against 3 free parameters --
overdetermined and non-trivial). Ground term: 3F -- matches Hund's rule.

## Part B -- dimensionless structural invariants (F0,A independent)

p-shell (F2 units): term count = 3; degeneracy multiset = {9,5,1}; terms
with S>0: 1 (only 3P).

r1 = (E(1D)-E(3P)) / (E(1S)-E(3P)) = 6/15 = **2/5** (exact, rational)
r2 = (E(1S)-E(1D)) / (E(1D)-E(3P)) = 9/6 = **3/2** (exact, rational)

d-shell: term count = 5; degeneracy multiset = {21,9,9,5,1}; terms with
S>0: 2 (3F, 3P). All pairwise energy-difference ratios among the 5 terms
reduce to ratios of small integers once expressed in A,B,C (e.g.
(E(3P)-E(3F))/B = 15, a pure SO(3)/Racah integer, no C dependence since
both are S=1 terms).

## Part C -- pre-registered fingerprint tests

**Test 1 -- phi-in-energies: FAIL to find any (PASS as a proven null).**
r1=2/5, r2=3/2, and every other ratio produced by Parts A/B are ratios of
integers (rational numbers), because every diagonalized term energy is a
rational linear combination of F0, F2 (and F4). phi=(1+sqrt5)/2 and sqrt5
both have degree-2 minimal polynomials over Q (x^2-x-1 and x^2-5
respectively) and are therefore NOT rational. A rational number cannot
equal an irrational number. Hence no term-energy ratio in p^2 or d^2 can
equal any element of {phi, 1/phi, phi^2, 1/phi^2, 2/phi, phi/2, sqrt5,
1/sqrt5, phi/sqrt5} -- this is a theorem given the exact rationality of
the diagonalized spectrum, not a scan that happened to come up empty.

**Test 2 -- fusion-multiplicity map: FAIL (no exact match; one shared
small integer).** p-shell term count (3) equals SU(2)_3's non-vacuum
primary count (3). That is the only exact integer coincidence found.
Nothing else lines up: the degeneracy multiset {9,5,1} does not match
the primary count (4), the non-vacuum count (3), or any quantum-dimension
combination (dims are {1,phi,phi,1}, D^2=2 phi sqrt5, both irrational --
and by Test 1's theorem no p/d-shell rational number can ever equal
them). No fusion coefficient N_ab^c (all 0 or 1 in SU(2)_3) or S-matrix
entry is reproduced by any p- or d-shell quantity. A-priori target count:
roughly 10-15 small integer/rational SU(2)_3 quantities exist to test
against (primary count 4, non-vacuum count 3, a handful of small
dimension-squared sums, ~6 independent fusion coefficients that are
themselves just 0/1). Matching exactly one of them (3=3) among that many
candidates, using small integers that are common on both sides regardless
of physical content, is unremarkable and expected by chance.

**Test 3 -- truncation test: FAIL (no genuine k=3 cutoff).** Max S=1 and
max L=2 for p^2 are the ordinary Pauli-exclusion / SO(3) triangle-rule
bounds for 2 electrons in an l=1 shell, not a Hill-bound-style saturation
at 3 independent cooperative modes. The term count (3) is a Racah/SO(3)
counting fact (number of L,S combinations compatible with the Pauli
principle for l=1), unrelated in origin and mechanism to the WZW level
k=3 truncation (which forbids j>3/2 primaries in the fusion category).
There is no "number of independent modes" object in the single-atom
multiplet structure that saturates at exactly 3 the way the fusion
category's primary spin does.

**Test 4 -- the pentagon/5: appears, wrong origin.** The 1D term's
degeneracy is 5, and this IS meaningful -- but as 2L+1 with L=2, an
ordinary SO(3) orbital-degeneracy fact, not as k+2=5 (the WZW truncation
index). r1=2/5 is likewise a ratio of Slater-Condon integers (6/15
reduced), rational, hence by Test 1 provably distinct from any
k+2-related irrational invariant of the fusion category. The numeral 5
appears in both contexts for unrelated reasons: SO(3) representation
dimension vs. anyonic level truncation.

## Part D -- chance accounting and verdict

Enumerating the small structural invariants produced by the p+d shells:
term counts (3, 5), degeneracy multisets ({9,5,1}, {21,9,9,5,1}), S>0
counts (1, 2), and O(10) pairwise ratios, essentially all reducible to
single-digit-to-two-digit integers or simple rational ratios thereof.
Against these, SU(2)_3 offers O(10-15) low-complexity integer/rational
targets (primary count, non-vacuum count, small dimension combinations,
0/1 fusion coefficients) plus a family of provably-irrational targets
(quantum dims, D^2, S-matrix entries) that can never match by Test 1's
theorem. With that many candidates on each side, 1-2 coincidental integer
matches (3 terms = 3 non-vacuum primaries; 5 = 2L+1 happens to equal
k+2) are the statistically expected outcome of an essentially unconstrained
search over small integers, not evidence of shared structure.

**Verdict: (i).** The open-shell multiplet structure of p^2/p^4 and d^2
is fully accounted for by ordinary SO(3)/Racah angular-momentum algebra
and the Pauli principle. Every energy, degeneracy, and ratio computed
here is exactly rational, and Test 1 proves as a theorem that no such
quantity can carry a phi- or sqrt5-type signature; Tests 2-3 find no
exact correspondence to SU(2)_3 fusion data beyond shared small integers
that are the generic residue of any small combinatorial enumeration; Test
4's "5" is SO(3)'s 2L+1, not the fusion category's k+2. There is no
condensate-specific (phi, k=3 truncation, pentagon-identity) signature in
single-atom light-element open-shell multiplets beyond coincidental
shared integers.

This localizes any condensate-specific structure entirely to the
many-body COOPERATIVE regime -- the Hill-bound n_H <= k=3 claim,
open problem [3] of Paper~XI~\cite{Paper11} -- and NOT to single-atom
multiplet term structure. The 2p^4 correlation/exchange energy that sets
IE(O) is, on this evidence, generic Racah physics: whatever condensate
signature exists in the oxygen/light-atom ionization-energy anomalies
must arise from how the multiplet energies are USED or COMBINED across
many electrons/atoms (a cooperative-mode counting problem), not from the
single-configuration term algebra itself, which this computation shows
carries no exact SU(2)_3 fingerprint.
