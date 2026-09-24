# PAPER TODO (later): Q_H labelling conflict in Papers I & VII

Deferred fix, flagged by the user 2026-09-02. Full diagnosis in claude-hopfion.md §8 [O].

THE CONFLICT: Papers I & VII predate the Q_H={0..3} ladder (Papers XV-XVII) and use "Q=2" for the
condensate VACUUM (P1:thm:vacuum; P7 l.88,145,1099,1102,1117,1410,1435). The current ladder uses:
  Q_H=0 bare condensate (no winding) / 1 neutrino / 2 charged lepton (e,mu,tau) / 3 baryon (+pion byproduct).
So "Q_H=2" names the VACUUM in I/VII and the CHARGED LEPTON in XV-XVII; and the ladder's Q_H=0 (no winding)
contradicts Paper I's vacuum carrying winding.

GEOMETRY CLARIFIED (user, 2026-09-02): Q_H=2 is genuinely TWO TUBES (an inner + an outer tube) -- P1:thm:vacuum
item (iii) fuses two Q=1 single tubes into the Q=2 singlet (j_1/2 x j_1/2 = j_0 + j_1), and the single tube is
UNSTABLE/collapses (l.224, l.1318, l.1599). Q_H=3 = the two pre-images reconfiguring into the T(2,3) trefoil.
The old "Paper I internal inconsistency on the single tube" is RESOLVED: the footnote l.323-326 is only a
SOLVER-MEASURE CONVENTION (z->-z bilateral factor labels a single tube "Q=2" in the axisymmetric code), not the
physics. Convention vs physics, not a contradiction. The CROSS-PAPER label (vacuum vs lepton; Q_H=0 vs winding)
remains the open item.

UNDERLYING PHYSICS QUESTION (author decision, not derivable): does the vacuum CARRY WINDING (Q=2, per I/VII)
or NOT (Q_H=0 bare, per the ladder)? This must be settled before renumbering.

STATUS OF P7:rem:cdm_ontology (added this session): uses the LADDER convention (Hopf charge distinguishes the
knot/matter from the medium/DM), number-free -- so it sits in latent tension with P7's body ("Q=2 = vacuum")
but with NO explicit contradiction (the remark cites no Q value). Acceptable interim.

QUARK-SECTOR CHECK (2026-09-02): is "bare condensate = Q_H=0 = no winding" load-bearing in XV-XIX, or would
anything change if the ground state carries winding (Paper I: 2I two-tube Q=2)? ANSWER: NOT load-bearing --
it is a LABELLING ORIGIN. [E, verified in XV/XVI/XVII]:
  - XVII l.138-141 states the framework assigns each sector Q_H=N its identity via a TRIPLE (McKay group, WZW
    CFT, knot type). Every quark result flows from that triple, never from the absolute integer or base winding:
    colour exclusion = group coprimality gcd(2,3)=1; generations/masses = WZW/coset M(5,6) CFT; confinement/
    fractional charge = trefoil rep theory; Q_group{3,6,10}, Q_coset=h(E8)=30 = knot (p,q)=(2,3) + Brieskorn
    Sigma(2,3,5).
  - Q_H and physics-winding are ALREADY decoupled: at Q_H=3 the trefoil is (p,q)=(2,3) [label != winding];
    Q_group(3,6,10) != Q_H(3,1,2).
  - "Vacuum" is overloaded 4 ways (CFT j_0 identity primary / geometric 2I Q=2 ground state / internal n=0
    fully-confined trefoil / XVI l.1047 relative "Q_H=2 plays the vacuum role"). The physics uses the CFT j_0,
    which is NOT the geometric ground state -> a winding-carrying geometric vacuum never enters it.
  - XVI already thinks RELATIVELY (2I->2T breaking; confinement hierarchy Q_H=1 subset 2 subset 3 = successive
    excitations).
CONCLUSION [N]: nothing in XV-XIX changes if the ground state carries winding. The fix is purely the LABEL
convention -- read Q_H as a SECTOR INDEX / excitation-count above the condensate (relative), not an absolute
Hopf charge from a winding-free zero. Then Paper I (2I/Q=2 winding-carrying ground state) and XV-XIX (sectors by
their triples) are consistent. Residual author-decision [O]: define Q_H absolute (additive homotopy invariant ->
integers offset by the vacuum winding = a renaming) vs relative (drop the "no winding" gloss on Q_H=0) -- NO
downstream physics either way.

Q_H=1 PAPER (XVII) CLOSE READ (2026-09-02) -- CORRECTS the "purely a label, no assumption" claim above:
XVII DOES contain an EXPLICIT trivial-vacuum assumption, in the topology section (the sharpest place, since Q_H=1
is the lowest excitation):
  - l.278-284: Q_H DEFINED as the ABSOLUTE Hopf charge (linking number of preimages; generator of pi_3(S^2)=Z).
    Not a relative index.
  - l.305-307: "Finite energy forces n to approach a CONSTANT at spatial infinity" (one-point compactification).
    This constant-at-infinity BC IS the winding-free-vacuum assumption: the soliton is counted against a trivial
    background. So XVII frames Q_H as absolute-on-trivial-vacuum -- OPPOSITE to XVI l.1047 ("vacuum = Q_H=2").
    => a real FRAMING inconsistency between the two ladder papers; the assumption is NOT merely a label.
BUT physics-inert [E for derivations, N for invariance]: all neutrino results flow from McKay group + WZW/coset
CFT + knot (Z2<->A1<->SU(2)_1; colour excl gcd(2,3)=1; Q_group=6; k=1 forced by McKay Q_group=2(k+2), l.1038;
dim_q=phi; coset M(5,6); masses; PMNS) -- none use the absolute Hopf number or the vacuum winding. Redefining Q_H
RELATIVE to a wound background (relative homotopy: linking of excitation preimages vs the background, instead of
compactify-to-constant) still gives Z with the SAME knot generators (unknot 1 / torus 2 / trefoil 3); only the
absolute integers + framing shift.
RECONCILIATION [I/N]: the two "vacua" may be DIFFERENT objects -> possibly NO contradiction. XVII's "constant at
infinity" = the LOCAL asymptotic reference for a localized particle soliton; Paper I's "Q=2 vacuum" = the GLOBAL
space-filling condensate ground-state texture (2I two-tube). A localized neutrino locally approaches the constant
while the global condensate is wound. Fix = restate XVII's "constant at infinity" as "approaches the condensate
ground state at infinity" = a RELATIVE Hopf charge. Precise lines to touch if reconciling: XVII l.278-307.

========================================================================================================
RESOLVED (2026-09-02): convention adopted -- Q_H is the Hopf charge RELATIVE to the condensate ground state.
User approved (plan: /Users/frederick/.claude/plans/twinkling-noodling-seahorse.md). Q_H=0 = the wound Q=2
ground state (Paper I Thm 6.1); Q_H=1,2,3 = neutrino/lepton/baryon excitations. Forward refs to XV-XVII from
I/VII kept in PROSE (no bibitems added, per user); back-refs use \cite{Paper1}.
EDITS MADE (all verified: env balance, refs/cites resolve, no backslashed underscores):
  - XVII l.278-307: definition rewritten relative-to-ground-state (n -> n_vac at infinity; Q_H=Q_tot-Q_vac);
    proof BC "constant at infinity" -> "condensate ground state at infinity". P17:prop:qh1_unknot unchanged.
  - Paper I: new remark P1:rem:qh_relative after P1:rem:symmetry_breaking_noncircular (states the convention;
    disambiguates from the z->-z solver-measure footnote).
  - Paper VII: cross-ref sentence at end of P7:rem:cdm_ontology (Hopf charge counted relative to the Q=2
    ground state; ladder on same footing).
  - Paper XVI: l.1047 "effective vacuum of 2I->2T transition, not the condensate ground state"; l.2052
    "bare (unwound) trivial vacuum".
  - Paper XV: footnote at first Q_H mention (convention + back-cite Paper I Thm 6.1).
COMPILE: pdflatex not available in-session; structural checks pass; USER TO COMPILE each edited paper.
Reference number resolved: P1:thm:vacuum = Theorem 6.1 of Paper I (from main_paper3.tex:1013 commented \ref).

VERIFICATION -- does the vacuum winding affect XV/XVI/XVII results? (2026-09-02, full trace). ANSWER: NO,
and the check is SHARPER than the earlier group/CFT argument. Traced every Q_H / pi_3 / winding use:
  - Every QUANTITATIVE use of Q_H is the LOCALIZED (relative) charge (nu=1, ell=2, baryon=3):
    * XVII l.1161-1163: J_{2a}=32 Q_H/15, J_4=8 Q_H/5 ("Q_H times the single-circle result, one preimage per
      unit of Hopf charge"). XVII l.1232-1233: "N=2 for Q_H=2 (two preimage circles in the Hopf link), N=1 for
      Q_H=1 (single preimage circle)" -- explicitly the localized soliton's winding.
    * XVI l.2588: "charge-N Hopfion", N = localized toroidal winding. XVI l.3223: Faddeev-Niemi floor
      E >= C|Q_H|^{3/4}, bounding the localized trefoil vs decay to the LOCAL vacuum.
  - Every pi_3(S^2)=Z use is DISTINCTNESS-only (XV Step 3 l.2297-2304: each sector = one non-degenerate class ->
    forces multiplicity-1 j=3/2). The VALUE never enters.
  - NO additive vacuum-charge usage anywhere (only Q_tot-Q_vac is our own new XVII def line l.287). No VK bound
    uses an absolute total-field charge.
CONSEQUENCES:
  (1) No result changes under the relative convention -- confirmed by direct inspection, not just the triple
      argument. phi^9 sqrt(J_a J_4)=6,10; mass ratios; energy floor; sector assignments all use localized charges.
  (2) The convention fix is CORRECTIVE, not cosmetic: the COMPUTATIONS were already relative (localized preimage
      counts 1,2,3); only XVII's topology-section DEFINITION said "absolute, n->constant at infinity" -- which was
      inconsistent with the paper's own J-integrals. Our edit aligns the definition with the computations.
  (3) The wound vacuum MAKES the relative reading necessary: a space-filling wound ground state has infinite
      background INT|grad n|^2; the finite particle energies the papers compute ARE the excess-over-background
      (relative) quantities. Relative-Q_H is what makes the energies finite.
