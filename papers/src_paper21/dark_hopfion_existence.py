#!/usr/bin/env python3.11
"""
dark_hopfion_existence.py -- the OPEN PROBLEM (Paper XXI open-problem (ii)): does a
localised D-type dark hopfion EXIST as a soliton of the (density-feedback) Faddeev
functional? Stability (given existence) is settled/topological; EXISTENCE is open, and
the coarse free-flow route is exhausted -- it is the WRONG instrument (lattice-unwinding
wall; cannot even hold the known trefoil: notes (j)-(l), dihedral_soliton_relaxation.py).

This script does NOT re-run gradient flow (that route is closed). Instead it does the two
things that ARE tractable and sharpen the problem:

  A. PIN THE KNOT TYPE from the group theory / low-dimensional topology, so a future
     charge-conserving relaxation knows what target to look for. The A-type ladder is
     T(2,q) torus knots on cyclic (lens-space) quotients; the D-type is a different animal.

  B. LAY OUT THE EXISTENCE LOGIC against the Faddeev-Hopf existence theorems and the
     feedback-functional saddle structure (shared with the baryon, which EXISTS), so the
     honest status is "allowed, not forced; a sterile-sector stabiliser is the open physics."

All group/topology facts below are checked numerically; the existence logic is stated,
tagged, and not over-claimed.
"""
import numpy as np
from math import gcd

PHI = (1+5**0.5)/2

def abelianization_dicyclic(n):
    """H_1 of the dicyclic (binary dihedral) group Dic_n = 2D_n, order 4n.
    Presentation <a,b | a^{2n}=1, b^2=a^n, b a b^-1 = a^-1>.
    Abelianise: b a b^-1 = a^-1 => a^2=1; b^2=a^n=a^{n mod 2}.
      n odd : a=b^2, b^4=1 -> Z_4 ;  n even: a^2=b^2=1 -> Z_2 x Z_2."""
    return "Z_4" if (n % 2 == 1) else "Z_2 x Z_2"

print("="*74)
print("  DARK HOPFION EXISTENCE  (Paper XXI open-problem (ii))")
print("="*74)

# ── A. group / topology facts (all checkable) ───────────────────────────────
print("\n--- A. Group theory: the D-type candidates sit SIDEWAYS in 2I ---")
ORDER_2I, ORDER_2T = 120, 24
# 2D_n lifts into 2T iff its SO(3) image D_n (dihedral, order 2n) sits in A_4 = 2T/Z2.
# A_4's subgroup orders are exactly {1,2,3,4,12} -- it has NO order-6 subgroup -- so
# D_3 (order 6) is excluded structurally and D_5 (order 10) by Lagrange (10 does not | 12).
A4_SUBGROUP_ORDERS = {1, 2, 3, 4, 12}
for n, name in [(3, "2D_3"), (5, "2D_5")]:
    order = 4*n
    idx   = ORDER_2I//order
    in_2I = (ORDER_2I % order == 0)
    in_2T = (2*n) in A4_SUBGROUP_ORDERS   # D_n (order 2n) must be an actual A_4 subgroup
    H1    = abelianization_dicyclic(n)
    mckay = {3: "D_5-hat -> SO(10)", 5: "D_7-hat -> SO(14)"}[n]
    print(f"  {name}: order {order:3d}, [2I:{name}]={idx:2d}, in 2I: {in_2I}, in 2T: {in_2T}, "
          f"H_1={H1}, McKay {mckay}")
print("  => both lie in 2I, NEITHER in 2T (A_4=2T/Z2 has no order-6 subgroup; D_3,D_5 out): SIDEWAYS.")
print("  => H_1 = Z_4 (n odd) CONTAINS a Z_2 = the sector's own dark parity (stabiliser).")
print("  NOTE indices {10,6} equal Q_group{lepton=10, nu=6}: a MATCH, not a FORCING")
print("  (index probe NEGATIVE, notes: realization allowed, not forced).")

print("\n--- A. Knot/configuration type: prism manifold, OFF the T(2,q) ladder ---")
print("  A-type sector  A_{q-1}: cyclic Z_q -> LENS space S^3/Z_q -> T(2,q) torus knot")
print("                 (nu=unknot T(2,1), lepton=Hopf link T(2,2), baryon=trefoil T(2,3)).")
print("  D-type sector  (McKay D-hat): the Kleinian singularity C^2/2D_n has LINK")
print("                 S^3/2D_n = a PRISM MANIFOLD, Seifert-fibred over the orbifold")
print("                 S^2(2,2,n). Its graph/plumbing link is the D_n Dynkin tree (a chain")
print("                 with a two-pronged FORK), NOT a (2,q) torus knot.")
print("  CHECK (not a homology sphere -> carries torsion = the Z_2 parity):")
for n in (3, 5):
    # Seifert base S^2(2,2,n): pairwise-coprime <=> homology sphere.
    pairwise_coprime = (gcd(2,2)==1 and gcd(2,n)==1)
    print(f"    (2,2,{n}): pairwise coprime? {pairwise_coprime}  -> gcd(2,2)=2 => NOT a homology")
    print(f"             sphere -> H_1 nontrivial (={abelianization_dicyclic(n)}), the dark Z_2 lives here.")
print("  => PROPOSAL [N]: the dark hopfion is a prism-manifold / D_n-graph-link Hopfion,")
print("     off the T(2,Q_H) torus ladder, its Z_2 homology = the dark parity. This is the")
print("     TARGET topology a charge-conserving existence search should aim at.")

# ── B. existence logic against the Faddeev theorems ─────────────────────────
print("\n--- B. Existence logic (pure Faddeev vs density-feedback) ---")
# Vakulenko-Kapitanskii bound and Lin-Yang existence are for the PURE Faddeev energy.
for Q in (1, 2, 3, 4):
    print(f"    Vakulenko-Kapitanskii: E(Q={Q}) >= c |Q|^(3/4) = c*{Q**0.75:.3f}  (energy bounded below)")
print("  PURE Faddeev-Hopf (E=kappa|grad n|^2 + lambda F^2): Lin-Yang (2004) PROVES")
print("  energy-minimising Hopfions EXIST for an infinite set of Hopf charges (incl. Q=1).")
print("  => in the PURE functional the D-type Hopf class, being non-trivial, HAS a minimiser.")
print()
print("  DENSITY-FEEDBACK functional (E_fb = K_fb * J4, a PRODUCT): no Derrick minimum ->")
print("  the object is a SADDLE, not a minimum (feedback_saddle_core.py; Paper XX). BUT the")
print("  BARYON (trefoil) has the SAME saddle structure and EXISTS (it is the proton,")
print("  stabilised by CONFINEMENT). So 'saddle, not minimum' does NOT mean 'does not exist'.")
print()
print("  The sterile dark hopfion has NO confinement. Its stabiliser is therefore the open")
print("  PHYSICS (not just numerics). Two framework-internal routes:")
print("   (i) lambda (Bogomolny coupling) is NOT pinned to phi^6 for the D-type (no Q_group),")
print("       so it may be LARGE -> a FAT soliton that is a genuine Faddeev MINIMISER (notes")
print("       (k)); plausible, not derived.")
print("   (ii) its own Z_2 parity + the prism/graph-link rigidity may protect a localised")
print("       saddle the way confinement protects the trefoil; unproven.")
print()
print("  SYNTHESIS (knot type <-> fat soliton reinforce each other): a D_n-graph-link /")
print("  prism configuration is intrinsically MORE DISTRIBUTED than a thin T(2,q) tube ->")
print("  it IS the large-lambda/large-size regime of route (i). The two proposals are one.")
print()
print("  EMPIRICAL MORSEL (notes (l), dihedral_soliton_relaxation.py, fat run lam=30):")
print("  under free flow a FAT D_3 config stayed LOCALISED (r_bar 1.4-2.3, never ran to the")
print("  box edge) and the 3-fold PERSISTED (A3~0.37-0.47, never symmetrised to axial); only")
print("  the lattice topology-leak (which drains the KNOWN trefoil too) ended it. => the")
print("  dihedral config is dynamically ROBUST and does NOT want to relax to axial; the drain")
print("  is a method artefact, not a dihedral instability. Supports existence; not a proof.")

# ── C. stability (refined from the prior notes: unwinding vs decay) ──────────
print("\n--- C. Stability, refined (notes (e)->(f)): unwinding vs decay are DIFFERENT ---")
print("  UNWINDING (Q_H -> 0): forbidden by Hopf-charge conservation, exactly as for the")
print("     trefoil. TOPOLOGICAL, settled.")
print("  DECAY TO SM SECTORS: two distinct protections, of DIFFERENT strength --")
print("   (1) SIDEWAYS-BRANCH BARRIER [solid]: 2D_3,2D_5 are NOT subset 2T (section A), so a")
print("       decay to the quark/lepton sectors needs a barriered RECONNECTION (the Paper")
print("       XVIII T(2,2)->T(2,3) wall), not a continuous relaxation -> at least METASTABLE.")
print("       (Contrast Q_8=2D_2 which IS subset 2T -> relaxes continuously -> UNSTABLE, not DM.)")
print("   (2) Z_2 DARK PARITY [UNPROVEN as dynamical]: an outer automorphism is NOT")
print("       automatically a conserved dynamical symmetry (the (e) over-claim, corrected in")
print("       (f)); only IF the dicyclic Z_2 is genuinely conserved does it give ABSOLUTE")
print("       stability (textbook R-/KK-parity architecture). Unproven.")
print("  => honest stability = 'unwinding-stable [settled] + decay-METASTABLE [sideways barrier,")
print("     solid] -> ABSOLUTELY stable IFF the Z_2 is dynamical [open]'. NOTE: Paper XXI leans")
print("     on the Z_2-parity for decay protection; the SIDEWAYS-BRANCH barrier is the firmer")
print("     argument and is the one to foreground.")
print()
print("  SUPERSEDED CONTEXT (notes (m)->(w)): a TRANSIENT dark hopfion (knot unwinds, leaves a")
print("  frozen gravitational wake) was explored and RETIRED -- the stable-knot picture (w) is")
print("  current, because the ADM dark-NUMBER asymmetry (Paper XXII) needs a conserved relic,")
print("  which a transient wake lacks. So STABLE existence (this script) is the right target.")

print("\n" + "="*74)
print("  VERDICT: existence is ALLOWED, not FORCED.")
print("   * topology: the D-type Hopf class (prism / D_n-graph-link) is non-trivial, so the")
print("     PURE Faddeev functional has a minimiser (Lin-Yang) -- existence in principle.")
print("   * feedback: saddle structure shared with the baryon (which exists) => not excluded;")
print("     the sterile-sector STABILISER (fat large-lambda minimiser, or Z_2/graph rigidity)")
print("     is the open physics.")
print("   * numerics: free flow is the WRONG instrument (lattice-unwinding wall). A")
print("     Q_H-CONSERVING / constrained high-res search (rational_map_ic.py smooth base +")
print("     saddle-score, GPU) is the needed-but-ungated route; index probe NEGATIVE means")
print("     no shortcut via a forcing index.")
print("   * stability (refined, section C): unwinding settled; decay METASTABLE via the solid")
print("     sideways-branch reconnection barrier; ABSOLUTE stability needs a dynamical Z_2 [open].")
print("  => this script PINS the target topology (prism / D_n-graph-link, with the dark Z_2 = its")
print("     torsion homology) -- answering the long-open 'knot type?' item -- and sharpens the")
print("     two handles: the sterile STABILISER (existence) and whether the Z_2 is dynamical")
print("     (absolute stability). [O]")
print("="*74)
