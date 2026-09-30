# Paper XX O2: the gen-1 field Phi (isospin sign-flip) -- state after the (c) run (2026-09-28)

O2 = derive Phi, the charge-linear field that flips ONLY generation 1 (u<d, while c>s, t>b).
Mechanism (settled, P20:prop:gen1_flip): dm = -Q*Phi; up(+2/3) down-shifts, down(-1/3) up-shifts,
dm_down-dm_up = +Phi > 0 (down heavier); Phi ~ MeV, generation-independent, so ~50% of m_u at gen1,
negligible at GeV -> flips gen1 only. Isospin SCALE (M_up/M_down=6.38) separately CLOSED (P18 tangent
5.29 x ~1.2 running, zero params). See quark_generation secs 25-29 for the full development.

## (c) RESULT -- MAGNITUDE is a product of two DERIVED scales [N, this session]
Phi = alpha_em * Lambda_jam, and BOTH factors are framework-derived (no new parameter):
  - alpha_em^{-1} = 360/phi^2 = 137.5  (icosahedral 2I/WZW value; measured 137.036, 0.3% off)
  - Lambda_jam = sigma(3R0+R)/3 = 308.3 MeV  (jam per-quark = m_p^jam/3, confinement energetics --
    the SAME jammed scale that anchors O1)
  => Phi = alpha_em * sigma(3R0+R)/3 = 2.24 MeV = 89-90% of observed m_d-m_u = 2.51 MeV.
This upgrades Phi from "alpha_em x a motivated formation scale" to "a product of two independently-
derived framework scales." Residual ~10% = formation-scale/running ambiguity (jam 308 vs the ~344 MeV
that would saturate the splitting; no principled scale lands on 344 -> not number-matched).

## SIGN -- confirmed a BOUNDARY CONDITION (not derivable from geometry alone)
The trefoil geometry supplies the z-level asymmetry a signed field couples to: up-type sites at z=+r0
(<z>=+0.87), down-type at z=0 -- genuinely asymmetric (NOT up=+r0/down=-r0), so a z-linear field
distinguishes them by construction. BUT no framework quantity has a clean definite SIGN there that fixes
which member is raised (the T_z=-0.52 seen at first was a biased distal-only sample; over all z=0 sites
<T_z>~0). So the sign stays the sector's ONE isospin input -- strictly less than the SM, which inputs the
FULL doublet splitting (sign AND magnitude of m_d-m_u) as a Yukawa. Reframed in P20:rem:gen1_sign.

## WHAT A FULL ACTION-LEVEL DERIVATION STILL NEEDS (the genuine open piece)
Show the density-feedback action, WITH the Paper XIII EM sector (photon = Q=4,J=1 composite), PRODUCES
the charge-linear z-level coupling at strength alpha_em*Lambda -- i.e. derive the coupling, not motivate
it. That would fix both the residual 10% and the sign. Multi-paper (needs the EM coupling to the director
at the z-levels); NOT attempted here beyond confirming the magnitude is compositional and the sign is a BC.

## ACTION-LEVEL + SIGN-DERIVATION HUNT (2026-09-28) -- structure clarified, NOT cracked
- **Action level NOT derivable with tools at hand.** Both the charge-linear -Q*Phi and the
  mass-proportional-to-tangent-orientation are MOTIVATED not derived (P18:rem:isospin_mass_ratio_scope
  says so). TRAP flagged (CLAUDE.md sec.10): the semi-Dirac dispersion (hard axis=heavy) is the WRONG
  instrument -- it bounds FLUCTUATIONS, not the coherent-configuration (mass) energy; do NOT chase it.
  The config-energy candidate (curvature/bending) goes the WRONG way (sec.28: integrated bending 0.44,
  inverted). No clean action-level derivation available; needs Paper XIII's EM coupling to the z-levels
  explicitly (not worked out). Genuine multi-step open derivation.
- **The sign is a Z_2 SSB DIRECTION = the TREFOIL CHIRALITY [key reframe].** P18:conj:isospin is explicit:
  WITHOUT the external field both I_3 values are equally probable (isospin-symmetric); the field SELECTS.
  What is selected is geometrically definite: the Z_2 writhe asymmetry puts the crossing vertices at
  z=+r0 (not -r0), and that orientation IS the trefoil handedness (mirror image -> z=-r0). So "the isospin
  input" = the trefoil CHIRALITY, ONE discrete bit, not a continuous parameter -- strictly sharper than
  the SM (which inputs the number m_d-m_u). Magnitude+structure then follow (this note, top).
- **LEAD [O, speculative, unverified]:** that one bit may not be independent. The framework has other
  discrete-symmetry directions -- weak parity (left-handed doublets) and matter/antimatter (Q_H=+3 not -3).
  If the trefoil handedness is fixed by the SAME choice as weak parity, m_d>m_u is a CONSEQUENCE of parity's
  direction, not a separate input -> collapses several SM sign-inputs into one Z_2. The framework does NOT
  currently connect handedness<->parity, so this is a direction for future work, not a result. Do not
  overclaim in the papers; the committed statement stays "sign = single isospin input" (P20:rem:gen1_sign).

## PAPER XX EDITS MADE (2026-09-28)
- P20:prop:gen1_flip: Phi range 1.6-2.3 MeV (formation scale 237-308), jam value 2.3 = best (was a bare
  1.6 that mismatched the remark).
- P20:rem:gen1_sign: reframed -- the sign is the sector's SINGLE isospin input, < SM Yukawa (which inputs
  sign+magnitude); geometry gives the z-asymmetry, not the sign; full action derivation would remove it.
- P20:rem:gen1_np: magnitude = product of two derived scales (alpha_em=360/phi^2, Lambda_jam=sigma(3R0+R)/3),
  no new parameter, 2.24 MeV.
- O2 open-problem: sharpened to "the action-level step" + Paper13 cite (bibitem added verbatim from
  main_paper5/14, inserted between Paper12 and Paper15). No dangling refs; spans clean.

## ISOSPIN-FROM-TWIST-SIGN (user, 2026-09-28) -- best lead on the sign, mechanism identified
User: "trefoils twist two ways... on formation there may be an isospin from the twist." Developed:
- MAPS to the framing-twist SINE-GORDON already in the paper (P18:prop:twist_sinegordon,
  V=V0(1+cos phi), kink=2pi twist). A sine-Gordon has KINK vs ANTIKINK = twist two ways.
  => **generation = |Tw| in {1,2,5} (magnitude); isospin = kink/antikink (SIGN of the twist).**
  Generation + isospin = magnitude + sign of ONE object, the framing twist. Framework-native.
- CONSISTENT: a PURE sine-Gordon has kink/antikink DEGENERATE (mirror images) = "isospin-symmetric
  without a field" (P18:conj:isospin). Verified: the isolated-cable Faddeev computation
  (cable_twist_potential.py) gave ONLY V0(1+cos phi), NO linear term -> statically symmetric. Correct.
- SPLITTING MECHANISM: lifting the kink/antikink degeneracy needs a term LINEAR in the twist gradient
  (a chiral/cholesteric bias lambda d_z phi), which shifts kink/antikink by +-2 pi lambda -> splitting
  LINEAR in twist, SIGNED by lambda = the -Q*Phi structure. Source of lambda = "on formation": the
  writhe injection with definite (2,3)-cable-framing handedness (P18 l.536,853, "minimal-action framing,
  not free") IN the chiral 2I condensate (vacuum 2I is chiral). So P18:conj's abstract "external field
  that selects" = the FORMATION writhe injection; its sign = the condensate/framing chirality.
- WHAT THIS DERIVES vs INPUTS: derives isospin=twist-sign, gen=|Tw|, both signs exist (doublet),
  splitting linear-in-twist, ordering sign = condensate chirality. STILL input: the 2I condensate
  chirality itself (ONE discrete choice, SAME as the writhe/framing, candidate = weak parity). So it
  reduces "the isospin input" to the vacuum chirality, not a new parameter.
- COMPUTED [E, 2026-09-28] `cable_twist_chiral.py`: the elastic split is a CLEAN NEGATIVE.
  Built a 3D helical two-strand cable (winding rate tau = chirality) and compared kink (q=+1) vs
  antikink (q=-1) twist energy. RESULT: EXACTLY degenerate (split ~1e-16 relative) at tau=0,+1,-1;
  parity checks E(+tau,+q)=E(-tau,-q) hold to machine zero. The Faddeev energy does NOT split the twist
  directions, even on a chiral cable. REASON (a symmetry, robust beyond the model): E=K*J4 is
  PARITY-INVARIANT, relating (kink,right-knot)<->(antikink,left-knot); a same-knot kink/antikink split
  would need parity broken WITHOUT flipping the knot, which the parity-invariant elastic energy cannot do
  (in the symmetric cable it is strand-exchange = 180deg z-rotation; in general it is parity).
- CONCLUSION: isospin = twist direction (kink/antikink) is GEOMETRICALLY correct, but the SPLITTING is
  NOT elastic. It must come from the CHARGE (EM -Q*Phi, which breaks the strand symmetry via up/down
  charge -- exactly P20:prop:gen1_flip) or an explicit parity-violating term. The elastic writhe-twist
  route is RULED OUT. And since the elastic sector is parity-symmetric, the sign is inseparable from
  PARITY VIOLATION -> ties the isospin ordering to the same parity input as the weak sector (the earlier
  unification lead, now sharpened: not an elastic chirality, a genuine parity-violation link).
- Magnitude relation to Phi=alpha_em*Lambda: OPEN whether the chiral-twist bias and the EM -Q*Phi are
  the SAME coupling (Q depends on I_3=twist-sign) or complementary. Do not conflate yet.

## PARITY-VIOLATING STRUCTURE = CHERN-SIMONS (user+synthesis, 2026-09-28) [lead, framework-native]
User: static energy is a mirror pair (parity) -> only charge breaks it; the parity violation must be
dynamical/topological; "perhaps during the Q_H=2->3 (T(2,2->3)) transition?"; + the 8-critical-point
chart (real eigenvalues=achiral source/sink/saddle vs COMPLEX=chiral spirals). Synthesis:
- The parity-violating structure is the **Chern-Simons term** -- ALREADY in the framework
  (P20:rem:charge_cs: SU(3)_1 CS, chiral edge modes carry the Z_3 triality charge = the fractionally-
  charged quark excitations; EM Hall response fractionalises charge into thirds; the 8pi CS spoke sets
  the constituent scale). CS properties that fit EVERY constraint:
  (1) PARITY-VIOLATING -- the archetypal 3D parity-odd topological term (chiral edge modes / QHE).
  (2) NOT in the static Faddeev energy K*J4 -- which is EXACTLY why cable_twist_chiral found exact
      kink/antikink degeneracy. CS is a separate topological term, invisible to the elastic test.
  (3) FIRST-ORDER in time / GYROSCOPIC -- the Hall response is a chiral rotational flow -> turns the
      achiral real-eigenvalue saddle into a SPIRAL (complex eigenvalues = the user's critical-point chart).
  (4) WHERE THE CHARGE ALREADY IS -- the edge modes ARE the fractional charge. So "only charge breaks it"
      and "the parity violation" are the SAME object: CS supplies charge AND chirality. -Q*Phi = the CS
      charge-field coupling; sign(Phi) = sign(CS level).
- LOCATION: the CS/Hall chiral flow acts during FORMATION (Q_H=2->3 writhe injection), spiraling with a
  handedness set by the CS level -> selects the isospin. Matches "on formation" and the transition guess.
- UNIFICATION (one signed integer, the CS level = chiral-edge handedness = Hopf orientation Q_H>0):
  isospin sign <-> Chern-Simons chirality <-> fractional charge <-> matter(Q_H>0) <-> weak parity.
- STATUS: LEAD/synthesis, framework-native, resolves the degeneracy puzzle and locates the parity
  violation, but NOT a derivation. VERIFY: does the CS/Hall term, during the Q_H=2->3 transition, give
  (a) a spiral (complex-eigenvalue) flow, and (b) an isospin selection of the right SIGN and ~MeV
  magnitude? Needs the CS coupling in the transition dynamics (Paper XIII EM + the SU(3)_1 CS edge
  dynamics) -- the same action-level step O2 already flags, now with a concrete parity-odd term to use.
  Do not put in the papers as derived; committed statement stays "sign = single isospin/parity input."

## ALGEBRAIC LOCK of the isospin sign to the CS level (2026-09-28) [best route, per "O2 is algebraic"]
CS parity flips k->-k (and conjugates topological spins). The framework sits at a DEFINITE level:
SU(3)_1, k=+1, c=+2>0, topological spin e^{+2pi i/3} (h=1/3, not the conjugate), framing period mod 3
(P20 l.217). The isospin ordering = sign(-Q*Phi) with Phi the CS/Hall response, so:
  **isospin sign = sign(k) = sign(c) = phase of topological spin = chiral-edge direction =
   fractional-charge sign = Q_H>0 (matter).** ONE signed CS datum.
- ACHIEVES: the isospin sign is NOT independent -- it is the CS level sign, LOCKED to the fractional
  charge, topological spin, framing-mod-3, and matter. Collapses "a free bit" -> "the same vacuum-chirality
  bit already spent on charge/matter." This is the algebraic form of the sign result.
- CANNOT: derive k=+1 vs -1 (that IS the vacuum-chirality parity input). But it is demonstrably ONE
  shared bit, not a separate input. Honest ceiling.
- e8_handedness.py (prior): the GENERATION twist sign = mass hierarchy, NOT chirality; gen-1<->gen-2 are
  twist-mirrors at SAME isospin -> confirms isospin is the separate CS/writhe Z_2, not the elastic twist.

## BUILD SCOPE -- CS/Hall spiral test (gated on Phase 0)
Goal: show the CS/Hall term turns the achiral Q_H=3 saddle into a CHIRAL SPIRAL (complex eigenvalues),
handedness selecting isospin. Foundation: `crossing_transition_v2.py` (Q_H=2->3 via fiber winding
chi+2t->chi+3t, gradient flow + Whitehead Q_H tracking) = the achiral baseline (gradient flow -> REAL
eigenvalues, no spiral).
- **Phase 0 [GATING CRUX, algebraic]:** derive the explicit CS/Hall coupling term in the director
  dynamics from SU(3)_1 CS + EM sector. NOT worked out in the framework; everything downstream needs it.
  Main risk. THIS is also the algebraic sign derivation (show CS/Hall -> -Q*Phi with sign(Phi)=sign(k)).
- **Phase 1:** reuse crossing_transition_v2 + tube pipeline; confirm pure-gradient transition has REAL
  eigenvalues (baseline, cheap).
- **Phase 2:** add CS term (d_t n = -grad E + lambda_CS[Hall]); linearize at saddle; eigenvalues REAL
  (CS inert) vs COMPLEX (spiral, chiral). sign(Im lambda) = handedness.
- **Phase 3:** does spiral handedness select isospin with right sign + ~MeV magnitude?
- EFFORT: multi-week, front-loaded on Phase 0 (itself algebra). RECOMMENDATION: do Phase 0 as PURE ALGEBRA
  first (CS-level derivation of -Q*Phi's sign); if it closes, the sign is algebraically reduced to the CS
  level and the spiral sim is confirmation, not load-bearing. Matches the framework's O1/up-ratio pattern
  (algebra closed them, not field simulation).

## PHASE 0 DONE (2026-09-28) -- ALGEBRAIC CLOSURE + PARITY NO-GO -> SKIP THE BUILD
Worked the CS algebra. SU(3)_1: c=8/(1+3)=2, primaries 1/3/3bar with h=0,1/3,1/3, topological spin
theta=e^{2pi i/3}, framing anomaly e^{2pi i c/24}=e^{i pi/6}, framing period mod 3; parity: k->-k, theta->conj.
- STRUCTURE DERIVED: formation = writhe injection = +1 unit of self-linking = +1 CS FLUX quantum with the
  (2,3)-cable handedness; a charge Q in that flux gets a CS/Aharonov-Bohm energy ~ Q*Phi_flux = the
  charge-linear -Q*Phi, sign(Phi)=sign(injected flux)=sign(k). The -Q*Phi structure falls out of CS. [E]
- PARITY NO-GO on the absolute sign [the key result]: SU(3)_{+1} and SU(3)_{-1} are EXACT mirror theories
  (both unitary, related by parity); the two cable handednesses are mirror images with EQUAL action, so
  minimal-action CANNOT choose. Nothing parity-even can: the isospin ordering is parity-ODD, and the
  framework's parity-even sector (elastic energy, computed DEGENERATE) structurally cannot fix a parity-odd
  sign. The only parity-odd object is the CS term, whose sign = the vacuum chirality BY DEFINITION.
- CLOSURE: isospin sign = sign(k), (1) NOT independent -- locked to fractional charge Y=t/3, topological
  spin e^{+2pi i/3}, framing anomaly c=+2, matter Q_H>0 (one shared vacuum-chirality bit); (2) NOT
  derivable -- that bit IS the parity input; a parity no-go forbids deriving it from the parity-even sector,
  and the CS mirror theories are equally valid. (3) => the SPIRAL SIMULATION (Phases 1-3) CANNOT do better;
  it would only reconfirm the lock. **SKIP THE BUILD.** O2 sign is closed: one shared vacuum-chirality bit,
  provably not an independent parameter and provably not derivable within the framework's chirality budget.
- Honest ceiling reached; matches the SM (which also takes a chirality/sign as input) but with LESS input
  (one shared bit vs per-doublet Yukawa sign+magnitude; magnitude is compositional, this note top).

## STATUS TAGS (CLAUDE.md sec.6)
- mechanism (linear Q*Phi, gen1-only): [E] established. isospin scale 6.38: [E] closed.
- magnitude Phi=alpha_em*Lambda_jam=2.24 (90%): [N] compositional derivation, two derived scales.
- sign: [O]/boundary-condition (one isospin input). full action-level Phi: [O] (needs Paper XIII EM).
