# L3 draft: the coupling as vacuum-channel monodromy fraction

Claim to derive: alpha^-1 (golden-angle attractor) = (full monodromy
cycle) x (vacuum-channel fusion probability) = 360 x p_0 = 360/phi^2.

Status: derivation sketch with verification points; not established.
Pre-registered failure modes at the end.

## Setup

Photon exchange between two charges is, in the condensate, an exchange
process between two Q_H = 2 flux knots threading the SU(2)_3 vacuum:
topologically a Hopf-link event T(2,2) between two j = 1/2 anyonic
worldlines (each Q_H = 2 hopfion carries the fundamental anyon charge
of the fusion sector; Papers XVIII/XIX satellite and reaction
machinery). The coupling is the statistical phase accumulated per
exchange, weighted by the channel in which the pair propagates.

## Step 1 — cycle closure (L1, exact)

One braiding of the pair in channel c = 1 advances the pair phase by
the monodromy angle 2*pi*(h_1 - 2h_{1/2}) = 2*pi/10 = 36 deg. The
condensate charge is Q = 10 = ord(q_{2I}) (Paper I). A charge-Q cycle
therefore accumulates exactly 360 deg: the coupling's angular
normalization is the closure of Q monodromies, not a unit convention.

## Step 2 — channel weighting (L2, exact)

The exchange pair fuses as 1/2 x 1/2 = 0 + 1 with Born weights
p_c = d_c/(d_{1/2})^2: p_0 = 1/phi^2, p_1 = 1/phi (sum = 1 exactly).
The vacuum (c = 0) channel is the one in which the exchanged pair
annihilates back to the condensate ground state — the channel that
defines a completed photon exchange between external charges, as
opposed to the c = 1 channel in which the pair remains a propagating
excitation.

## Step 3 — the bridge (the open part)

Proposed statement: the measured coupling counts, per full cycle, the
fraction of exchange events that complete (return to vacuum):

    alpha^-1 = (cycle) x P(vacuum channel) = 360 x 1/phi^2.

Physical reading: alpha = phi^2/360 is the probability-weighted phase
cost per completed exchange; strong coupling would mean every braiding
completes, alpha^-1 = 360; the golden ratio suppression is the
Born weight of annihilation.

What a real derivation must supply:
(a) Why the coupling is LINEAR in p_0 (an amplitude-squared/Born
    argument at the level of the exchange density matrix), rather than
    p_1, p_0 - p_1, or |amplitude| = 1/phi.
(b) Why the phase unit is the c = 1 monodromy (the propagating
    channel's angle) while the weight is the c = 0 probability —
    i.e. the exchange advances phase while propagating and completes
    by annihilating. This mixed structure is the crux; it must drop
    out of a single well-defined TQFT quantity, e.g.
    Tr[rho_exchange * U_monodromy] with rho the fusion density matrix:
    NOTE Tr[rho U] = p_0*e^{-i*108 deg} + p_1*e^{+i*36 deg} — compute
    its phase/modulus exactly (verification point V2 below) and check
    whether any natural reading gives 360/phi^2 before inventing a
    bespoke weighting.
(c) The degree normalization: Step 1 fixes 360 as Q monodromies; the
    derivation must show the coupling is measured in units of
    (cycle/360), i.e. per-monodromy, so that the angle appears in
    degrees. This is where the framework's alpha-as-angle
    identification (P3:rem:packing_layers) becomes a theorem or fails.

## RESOLUTION OF THE CRUX (V2 computed, exact)

The mixed structure of (b) dissolves. The single natural TQFT
quantity — the monodromy expectation in the exchange state —
evaluates to the vacuum-channel weight exactly:

    Tr[rho U_mono] = p_0 e^{-i 108 deg} + p_1 e^{+i 36 deg}
                   = (3 - sqrt 5)/2 = 1/phi^2   (exactly real; sympy)

and equals the standard interferometry identity
S_{1/2,1/2} S_00 / S_{0,1/2}^2 = 1/(2 cos(pi/5))^2 = 1/phi^2
(verified exact). This is the monodromy scalar of anyonic
interferometry (the fringe-visibility factor of interfering anyonic
charges): the imaginary parts of the two channels cancel identically
and the expectation is the real Born weight of the vacuum channel.
No channel projection is imposed by hand.

The bridge statement therefore sharpens to:

    alpha^-1 = (closure cycle) x Tr[rho U_mono] = 360 x 1/phi^2,

i.e. the coupling is the fundamental braiding cycle weighted by the
monodromy scalar of the exchange pair. Item (a) is answered (the
linearity in p_0 is the reality of the monodromy expectation, not an
assumption); item (b) is answered (no mixed weighting exists — one
operator, one state). Note also that Q = 10 closes BOTH channels
(c = 1 winds +1 turn: 10 x 36 = 360; c = 0 winds -3 turns:
10 x (-108) = -1080), so the 360 deg cycle is the winding-one
fundamental closure and is channel-independent as a closure
condition. Item (c) — why the coupling is measured in cycle units
(degrees), the alpha-as-angle identification — remains the open
axiom-level step.

## Verification points (bracket engine, exact)

V1: |Jones_{T(2,2)}(q_5)| = 1/phi (proved, Paper XIX): the Hopf-link
    amplitude modulus is 1/phi = p_1 = sqrt(p_0). Determine which of
    p_0 = |J|^2 (amplitude-squared -> Born) or p_1 = |J| the bridge
    actually needs; amplitude-squared favors (a).
V2: Tr[rho U_mono] = (1/phi^2) e^{-i 108 deg} + (1/phi) e^{i 36 deg}
    — exact value, phase, modulus (sympy).
V3: The S-matrix normalized Hopf amplitude S_{1/2,1/2}/S_00 = 1
    exactly — check and interpret (it says the linked exchange is
    unsuppressed relative to vacuum; the suppression lives entirely in
    the fusion weights).

## Pre-registered failure modes

If the natural TQFT quantity yields 360 x p_1 = 222.5, 360 x |J| =
222.5, 360 x (p_1 - p_0) = 85.0, or the phase of Tr[rho U] rather
than 360 x p_0 = 137.5 — the identification is wrong and this note
records a negative. No reweighting after the fact: the bridge must
select p_0 by argument (a), not by matching the target.
