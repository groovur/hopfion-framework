# Quark generation, the E_8 anchor, and generation = ribbon twist

Session 2026-09-03. Thread: the confined-quark mass is algebraic (no saddle) -> E_8 route ->
strange anchor -> isospin closure -> bijection chase -> three orthogonal quantum numbers ->
generation = ribbon twist (framing). Tags per CLAUDE.md sec.6.
Scripts: papers/src_paper16/e8_strange_anchor.py, papers/src_paper18/breaking_exponent_bijection.py,
         papers/src_paper18/e8_twist_framing.py (this step).
Papers: added P17:rem:quark_mass_algebraic (XVII); pointer at P16:rem:e8_negative_result (XVI).

## 1. The confined-quark mass is ALGEBRAIC, not a saddle [in papers now]
The Q_H=3 quark is metastable/transient with NO stable field saddle (Paper XV/XVI/XVIII). So its
mass cannot be a configuration (saddle) energy, as lepton masses are -- it must be ALGEBRAIC.
Paper V's E_8 route is that algebraic mechanism: m_q/m_ell(g) = exp(n_q*pi/9), n_q = Base_g - 2m,
Base_g = 180 T_g + 7/2, m = E_8 Coxeter exponent, pi/9 = pi Q/(k h(E_8)) (Q=10,k=3,h=30).
This reframes Paper XVI's "E_8 = failed mechanism": what failed was a SADDLE-based DERIVATION of the
assignment, not the algebraic architecture, which the instability makes MANDATORY.

## 2. Strange anchor: PARITY-forced, not fitted [derived this session]
n_q=0 (m_q=m_ell) requires m = Base_g/2. Base_1=41(odd), Base_2=26(EVEN), Base_3=17(odd). Only gen-2
gives an integer, and 26/2=13 is an E_8 exponent -> ONLY strange can anchor (m_s=m_mu). Equivalently
T^(13)=(4*13-7)/360=1/8=T_mu: strange's E_8 spin = the muon's -> most lepton-like quark. (Same 26=2*13
= two-loop QCD b_1 in Paper V.)

## 3. Mirror structure + isospin closure [derived this session]
E_8 exponents pair m+m'=h=30: (13,17),(11,19) complete (used), (1,29),(7,23) half (gen-3 uses 1,7;
23,29 unused). Mirror-pure by isospin: (13,17)=(s,d) down, (11,19)=(c,u) up, (t,b)=(1,7) mixed extremes.
ISOSPIN SPLIT reduces to one bit -> fixed by anchor: strange is DOWN because m_s=m_mu and the muon is
the I_3=-1/2 member of (nu_mu,mu-)_L; weak isospin inherited. No UP-type anchor possible (would need a
neutrino-mass ~meV resonance). This RESOLVES the up/down direction P18:conj:isospin leaves free.
Open: lift I_3-inheritance from a symmetry argument to a derivation from the action.

## 4. Bijection chase: color != generation [derived this session]
6 breaking patterns (P18:conj:six_quarks) factor as 2(single/double) x 3(Z_3 site = COLOUR, degenerate).
6 flavours factor as 2(isospin) x 3(GENERATION = E_8 exponents, Z_3-BROKEN, ~5-dex mass spread).
The two "3"s coincide in count but are DIFFERENT quantum numbers: the breaking Z_3 is COLOUR (crossings,
unbroken by confinement), so it cannot carry the Z_3-broken generation spread. => P18:conj:six_quarks
CONFLATES colour with generation. Generation is NOT breaking-geometry; it is the E_8/algebraic label.
Secondary tension: P18 has TWO isospin mechanisms (breaking DEPTH vs writhe z-LEVEL) -- both breakings
sit at z=0 midpoints, so they can't both give up/down naively. Needs reconciling.

## 5. THREE ORTHOGONAL QUANTUM NUMBERS, three geometric homes [the frontier]
User's principle: algebra needs a geometric structure to map onto (as 2I <-> two-tube hopfion). So the
E_8 generation label MUST map to a physical trefoil feature. The clean home (Calugareanu SL = Wr + Tw):
  - COLOUR      = WRITHE (Wr) = the 3 crossings (Z_3, SU(3)_1). Fixed by knot type T(2,3).
  - ISOSPIN     = the Z_2 up-down z-asymmetry (crossings z=+r0 vs midpoints z=0). Parity, orthogonal.
  - GENERATION  = TWIST (Tw) = the ribbon FRAMING. A fixed trefoil carries a family of framings (all
                 T(2,3), differing by twist) -> generation is the framing, varying while colour stays.
The algebra->geometry bridge already exists in the mass formula: the WZW T-MATRIX IS the framing anomaly
(Dehn twist / topological spin), and the E_8 exponent enters mass via T^(m)=(4m-7)/360, a framing phase.
So m -> T^(m) -> ribbon twist is a standard CFT chain: the T-phase the mass formula uses IS the twist.
Generation != breaking-geometry (that's colour x isospin); generation = framing-geometry (the twist).

## 6. RESULT: E_8 exponent = ribbon twist; mirror = EXACT twist reversal (e8_twist_framing.py)
Framing phase T^(m)=(4m-7)/360; mirror-pair sum T^(m)+T^(30-m)=53/180 CONSTANT -> centre c0=53/360.
Twist relative to centre:  Tw(m) = T^(m) - c0 = (m-15)/90  EXACTLY  -> the exponent is a LINEAR twist
coordinate, m = 15 + 90*Tw, centred on the untwisted point m=15=h/2.
MIRROR TEST PASSES EXACTLY: Tw(30-m) = -Tw(m) for all four pairs -> the E_8 duality m<->30-m IS ribbon
orientation/twist reversal (algebra<->geometry map, like 2I<->two-tube).
Twist per quark (Tw=m-15): u(+4) d(+2) s(-2) c(-4) b(-8) t(-14).
  - centre (Tw~0) = leptonic; strange (Tw=-2) near centre = most lepton-like (anchor status itself is the
    Base-parity COMBINATION Base_2-30=2*Tw, not twist magnitude alone -- d ties |Tw|=2).
  - gen-1 <-> gen-2 are TWIST-MIRROR partners: (u,c)=(+-4), (d,s)=(-+2), same |Tw| opposite sign; gen-3
    (t,b) the large-twist outlier. So |Tw| does NOT linearly order generations by mass -- it is structured
    (gen1/gen2 mirror, gen3 extreme). Isospin (up/down) is the orthogonal Z_2, not |Tw|.
  - HANDEDNESS: only the NEGATIVE-twist extremes are used (t=-14,b=-8); positive mirrors (23=+8,29=+14)
    UNUSED. Gen-3 picks one twist chirality -> a particle/antiparticle or L/R-looking asymmetry. NEW OPEN Q.
CAVEAT: this "twist" is the fractional FRAMING PHASE (T-matrix), not yet an integer Calugareanu framing
with the trefoil's actual Wr; making it a genuine integer ribbon framing is the next geometric step.

## 7. Handedness check (e8_handedness.py): deflates to mass hierarchy, but yields two real things
The gen-3 "one twist sign" is NOT a fundamental chirality. Sorting by twist: t(-14) b(-8) c(-4) s(-2)
d(+2) u(+4) is MONOTONE in mass (more-negative twist = heavier). Twist sign = light(+)/heavy(-) = mass
order. Unused +8,+14 (exps 23,29) = "ultra-light" slots below u,d that don't exist = the mass FLOOR.
So the handedness is the flavour hierarchy restated, not L/R. [initial chirality/CP excitement RETRACTED]
BUT the check paid off:
  (a) INTEGER-FRAMING DEPTH: E_8 exponents are all ODD -> Tw=m-15 is always EVEN. Even framing = ribbon
      closes with consistent orientation/spin structure. => the integer-framing step should use EVEN
      Calugareanu framings. SL = Wr + Tw = 3 + (m-15) = m-12 (odd, since Wr=3 odd + Tw even).
  (b) LEAD [speculative]: gen-1<->gen-2 are exact twist-mirrors ((u,c)=(+-4),(d,s)=(-+2)); gen-3 (t,b)
      UNPAIRED (mirrors 23,29 empty). Parallels CKM: (1,2) mix strongly (Cabibbo ~13deg), 3rd weakly.
      Twist-paired<->strongly-mixed, unpaired<->weakly-mixed. Caveat: twist-mirror is WITHIN isospin,
      CKM is up-down misalignment -- not airtight. Checkable against actual CKM angles later.
NEXT: integer-framing step with the EVEN constraint -- SL=Wr+Tw, Wr=trefoil writhe(=3), Tw even; test
whether the SL values {t:-11,b:-5,c:-1,s:1,d:5,u:7} (=m-12) are physically meaningful framings.

## 8. E_8 ratio vs twist correlation (e8_ratio_twist_correlation.py): isospin, not chirality
Q: does the E_8 mass ratio (pred/PDG at condensate scale) correlate with twist chirality? ANSWER: NO
robust twist correlation. Raw Pearson |Tw|=+0.67 is ENTIRELY top-driven (leave-out-top -> +0.02);
non-monotone (|Tw|=4 c,u dips below |Tw|=2 s,d); charm is an outlier. NOT a twist law.
ROBUST signal = ISOSPIN VARIANCE: down-type ratios {d:1.26,s:1.13,b:1.21} cluster TIGHT (x1.11 spread)
= correct up to a single uniform ~1.2 factor (running); up-type {u:0.67,c:0.34,t:1.93} scatter (x5.78)
= genuinely wrong ratios. Pearson-of-mean misses it (means 1.20 vs 0.98 close); the signal is variance.
WHY: anchor is DOWN (strange=muon, n_s=0) -> down sector calibrated -> ratios right; UP unanchored
(no up-quark = a lepton; needs ~meV neutrino resonance) -> scattered; up hierarchy steeper (harder).
So twist's role in ACCURACY is INDIRECT (twist->isospin->which sector anchored). Twist gives the
geometric LABELS; the anchor/isospin governs the ACCURACY -- cleanly separate. => integer-framing is a
statement about the labels and will NOT by itself fix the up-sector masses. PARTIAL WIN: the E_8 route
gets the whole DOWN sector right up to one uniform (running) factor; UP sector is the open problem.

## 9. INTEGER-FRAMING STEP (e8_integer_framing.py): PASSES -- generation = literal ribbon framing
Calugareanu SL=Wr+Tw, Wr=3 (right-handed trefoil, blackboard framing SL=Wr). Tw(m)=m-15, SL=m-12.
  - Every Tw is an EVEN integer, every SL an integer (odd) -> Tw=m-15 is a genuine integer ribbon
    framing, not a half-integer/phase artefact. Framing per quark: SL={t:-11,b:-5,c:-1,s:+1,d:+5,u:+7}.
  - Framing phase is exactly framing-linear: T^(m) = 53/360 + Tw/90 (verified). Phase advances by
    1/90 = 1/(k h(E_8)) per unit twist, about centre 53/360. Centre Tw=0 <-> m=15 (NOT an exponent)
    = the untwisted, generation-less LEPTONIC limit -- no quark sits there (a lepton is unknotted AND
    untwisted). Nice.
  - MIRROR = FRAMING CONJUGATION: Tw->-Tw (exact), SL->6-SL (reflection about Wr=3), every pair
    SL+SL'=2Wr=6. Flips the twist, holds the KNOT (crossings=COLOUR) FIXED -> a PURE GENERATION
    symmetry, not a full knot-mirror. Exactly the colour/generation orthogonality the picture needs.
  - MAP CLOSED at the 2I<->two-tube level: E_8 exponent <-> even ribbon framing; h=30 <-> twist range;
    mirror <-> framing conjugation. Generation-as-twist is now LITERAL geometry, not CFT-algebra.
CAVEAT [O]: not yet proven that 1/90-per-twist is the trefoil's ACTUAL Chern-Simons framing anomaly
(topological spin) vs an E_8 integer that is merely framing-shaped. Remaining first-principles step:
compute the trefoil's framing anomaly independently and match to 1/90.

## 10. RECON of existing framing content (Papers 16,18) -- a correction + two gifts (2026-09-03)
User asked to check for existing colour/ribbon twist + the 2T<->E_6 structure BEFORE the framing match.
(1) EXISTING FRAMING:
   - Paper XVIII l.466,536: trefoil GLOBAL self-linking = 3, the (2,3)-cable framing, "NOT a free choice:
     the MINIMAL-ACTION framing" (BPS-dictated). This is the COLOUR/knot framing, FIXED.
   - Paper XVI P16:prop:bishop_holonomy (l.2811): trefoil Bishop-frame HOLONOMY H = -63.0191 deg
     (-1.09989 rad), exactly Z_3-symmetric (-21.0064 deg/segment) -- the intrinsic parallel-transport
     ribbon twist, DISTINCT from the theta=3t coordinate winding.
   CORRECTION [important]: my SL=Wr+Tw=m-12 RE-FRAMED the whole trefoil, but the global SL is FIXED at 3
   (colour). So generation CANNOT be a re-framing of the whole knot -> it must be a PER-STRAND / internal
   twist, distinct from the fixed global colour framing. RETRACT the SL=m-12 (global) reading; KEEP
   Tw=m-15 + mirror=twist-reversal as a PER-STRAND twist. Two twists were blurred: global colour framing
   (SL=3, fixed) vs per-strand generation twist.
(2) 2T<->E_6 -- NOT a tension, a clean split:
   - Paper XVI l.581: m_q = m_ell[2I/E_8] * colour-correction[2T/E_6/SU(3)_1]. E_8 exponents enter as the
     LEPTON-inherited part (generation); E_6 is COLOUR. Exactly the three-structure split, by construction.
   - Paper XVII l.1511 (Paper V): CKM mixing from Coxeter exponents of E_6, E_7 AND E_8 (chain 2T<2O<2I).
     So mixing uses all three E-groups, richer than just E_8.
TWO GIFTS for the framing-anomaly match:
   (a) TARGET EXISTS: don't recompute the trefoil framing -- Paper XVI has it: Bishop holonomy -63.02 deg,
       Z_3-symmetric. Match = does our per-strand E_8-twist (phase 1/90 = 4 deg/twist) reconcile with -63
       deg + the SU(3)_1 topological spin (down-quark h=1/3 -> e^{2 pi i/3})?
   (b) MIXING MACHINERY EXISTS: the gen1<->gen2 twist-mirror / CKM lead -> check against Paper V's existing
       E_6/E_7/E_8 CKM construction, not a new one.

## 11. STRUCTURAL SHIFT (user, 2026-09-03): E_8 = global approx of per-strand E_6 twist
Hypothesis: the intrinsic per-strand generation twist is E_6 (2T/SU(3)_1, the quark's OWN group); E_8
(2I, lepton) is the coarse GLOBAL/lepton-facing approximation. Strongly supported (e6_perstrand_vs_e8_global.py):
  - E_6 has EXACTLY 6 Coxeter exponents {1,4,5,7,8,11} = 6 quarks, ALL used, NONE spare, 3 COMPLETE
    mirror pairs (1,11),(4,8),(5,7) -> twists {+-1,+-2,+-5}. |Tw| in {1,2,5} = generation; sign = isospin.
  - This DISSOLVES the gen-3 handedness/unused-exponent puzzle: the 2 spare E_8 exponents (23,29) were an
    artefact of using the 8-exponent GLOBAL group (E_8) for a 6-STRAND object. E_6 over-describes nothing.
  - E_8 still appears legitimately as the LEPTON BRIDGE: m_q = m_ell[2I/E_8]*colour[2T/E_6] references the
    lepton (E_8), and the strange<->muon anchor is an E_8 bridge coincidence (quark T-phase = LEPTON T-phase).
    Intrinsic per-strand twist = E_6. So the E_8 picture (anchor, mass) and the E_6 picture (per-strand twist,
    framing) are the SAME structure at two levels -- consistent with the architecture.
TRADE-OFF [honest]: the clean EVEN-framing result was E_8-SPECIFIC (E_8 exps all odd -> even twists). E_6
twists {+-1,+-2,+-5} are MIXED parity. E_6 buys 6=6 completeness, loses all-even framing. Needs interpretation.
LEAD [speculative]: E_6 |Tw|max = 5 = the pentagon k+2=5 that fixes 3 generations (claude-hopfion sec.3);
3rd generation = maximal twist = pentagon. Flag, not a claim.
REFRAMES the framing-anomaly match -> E_6-NATIVE: match the per-strand twist to the trefoil's OWN WZW
SU(3)_1=(E_6)_1 (c=6) + Paper XVI Bishop holonomy -63 deg -- NOT to E_8. Cleaner, more physical target.
TO CHECK: (i) does the strange anchor survive as a cross-group E_6<->E_8 coincidence? (ii) do the E_6
exponents reproduce the (down-sector) masses the E_8 route got right? (iii) the E_6-native framing match.

## 12. FRAMING PARITY CORRECTED (2026-09-03): mod-2 is a red herring; quark=mod-3, lepton=mod-5
User: shouldn't E_6 framing be even (per-strand) while E_8 is the lepton crossover? ANSWER: NO -- "even/odd"
(mod 2) is NOT the physical framing invariant. RETRACT my earlier "even framing = spin structure" claim:
the E_8 even-ness was pure exponent ARITHMETIC (E_8 exps all odd, h/2=15 odd -> Tw=m-15 even), not physics.
The framing anomaly is e^{2pi i h} per unit framing, so the periodicity = denominator of the topological spin h:
  - QUARK (per-strand, SU(3)_1=(E_6)_1, fundamental 3): h=1/3 -> period 3 (MOD 3).
  - LEPTON (SU(2)_3, 2I j=1): h=2/5 -> period 5 (MOD 5).
Neither is mod 2. The RIGHT invariant (mod 3) is CLEAN and confirms itself: E_6 twists mod 3 give residue
counts 0->0, 1->3, 2->3. residue-0 (colour-SINGLET framing) is EMPTY = no quark is a colour singlet (quarks
carry colour) -- exactly right. The 6 quarks split 3-3 into the 3/3bar (fundamental/anti-fundamental) classes.
=> the mod-2 "even" story tested the quark with the LEPTON's ruler. Sharpened two-level picture: framing
periods are 3 (quark/colour) and 5 (lepton/pentagon k+2=5) -- the framework's two structural numbers; the
E_6->E_8 crossover is mod-3 -> mod-5. [verified inline; residue-0-empty is a genuine consistency check.]

## 13. E_6 BRIDGE ORDERING (e6_bridge_ordering.py): |Tw|=generation, mirror=isospin doublet
Assigned the 6 E_6 exponents to quarks via the E_8 bridge. Result (up=+Tw convention):
  u:m7(Tw+1,gen1) d:m5(-1,gen1) | c:m8(+2,gen2) s:m4(-2,gen2) | t:m11(+5,gen3) b:m1(-5,gen3).
FORCED, two ways:
  - |Tw|=generation: {1,2,5}=gen{1,2,3}, distinct magnitudes (unlike E_8 where gen1,2 SHARE |Tw|={2,4}).
    Fixed by the E_8 bridge (gen-3 = largest twist in both).
  - MIRROR PAIR = ISOSPIN DOUBLET: E_6 pairs (5,7)=(u,d),(4,8)=(c,s),(1,11)=(t,b) = the PHYSICAL generation
    doublets. Contrast E_8 pairs = cross-generation (s-d,c-u), gen-3 unpaired. => STRONGEST confirmation E_6
    is the per-strand group: its mirror symmetry IS the up-down doublet structure; its |Tw| ARE the 3 gens.
NOT forced: which doublet member is up vs down (needs electric-charge/writhe); the mod-3 residue 3/3bar
meaning (residue-0 stays empty = no colour singlets, correct; 1-vs-2 split doesn't align with isospin/gen).
Masses (item 3) come via the E_8 bridge, not the E_6 twist. READY for the E_6-native framing match.

## 14. E_6-NATIVE FRAMING MATCH (e6_framing_match.py): topological PASS; classical Bishop separate
TOPOLOGICAL match PASSES: every E_6 twist through the SU(3)_1 anomaly (120 deg/unit) lands only on
120/240 deg (residues 1,2 = 3/3bar), NEVER 0 (singlet). So the per-strand twist is a mod-3,
SU(3)_1-consistent framing with the singlet framing empty = the trefoil's own Z_3 framing anomaly.
And E_6 IS the trefoil's WZW (McKay 2T<->E_6<->SU(3)_1), so generation=twist is the trefoil's OWN framing
-- DERIVED at the framework's standard (same as 2I<->two-tube, a McKay/WZW correspondence, not classical geom).
WZW spins: SU(3)_1 (3) h=1/3->120deg; (E_6)_1 (27) h=2/3->240deg; both mod 3.
GEOMETRIC (Bishop) SEPARATE + honest: H=-63.0191 deg -> Tw_geom = -0.17505 framing units. A near-coincidence
7/40=0.175->63.000 deg flagged by the script, but REJECTED: off by 0.019 deg, Bishop is converged-numerical,
and it is the total torsion of THEIR specific trefoil embedding (NOT a knot invariant) so a clean fraction
would be an embedding artefact. [CLAUDE.md sec.3 -- coincidence, not claimed.] The classical anholonomy
(-0.175, fractional) is Calugareanu-related to but NOT equal to the integer WZW twists {+-1,+-2,+-5}. It
corroborates a nontrivial Z_3 framing, not the quantization.
VERDICT: generation=twist DERIVED in the WZW/McKay sense (framework standard). OPEN (now precise): a
first-principles link from the classical Bishop anholonomy (-63 deg) to the INTEGER WZW framing -- i.e. show
the classical geometry QUANTIZES to the E_6 twists. That is the one rung left.

## 15. BISHOP RUNG -- CORRECTED (2026-09-03, user pushback): NOT accidental; root-siblings, not identity
FIRST took the Bishop holonomy as an embedding-dependent ACCIDENT -- WRONG. The trefoil geometry is
GOLDEN-DETERMINED: R0=3, r0=sqrt(2)/phi=0.874032 (Paper XV l.640, r0=R0/C2* from the virial fixed point;
r0_candidate_solver.py: sqrt(2)/phi is the tightest fit). So the Bishop holonomy is DETERMINED, not arbitrary.
HIGH-PRECISION (total torsion at golden r0, converged NT 50k-800k): H = -1.099975 rad = -63.0239 deg.
Near-forms ALL close but RESOLVED non-exact: 1.1 rad (=3*11/30, off 2.5e-5), 7/40 turn (off 6.6e-5),
11/30 rad/segment (off 8.7e-6). => H is a genuine TRANSCENDENTAL (torsion integral of the golden torus knot),
structural in ORIGIN (golden) but not in VALUE (no clean closed form). [CLAUDE.md sec.3 -- flagged, not claimed.]
CORRECTED VERDICT:
  - NOT accidental (user right): the Bishop holonomy and the E_6/WZW framings share the SAME ROOT -- both flow
    from the golden/BPS structure that fixes r0=sqrt2/phi and the trefoil. Consistent, not independent. This
    STRENGTHENS generation=twist: its geometry is golden-determined, not arbitrary.
  - NOT a clean quantisation either: the value is a transcendental ~-1.1 rad; Calugareanu ties this geometric
    torsion to the SINGLE SL=3 framing, not to the six integer generation framings {+-1,+-2,+-5} (separate
    WZW/E_6 shifts). A transcendental != an integer; classical torsion does not PRODUCE the topological framings.
  => root-siblings (both golden/BPS-determined), NOT accident and NOT identity. generation=twist grounded in
     the golden geometry at ROOT level (stronger than 'just McKay'); the classical value is transcendental torsion,
     not the integer framing quantised. OPEN: whether the near-1.1-rad/11-30 structure is exact under a better
     (analytic) torsion computation, or genuinely transcendental.

## 16. PAPER REVISIT of committed E_8 remarks (2026-09-03): SOUND, no changes
Re-checked P17:rem:quark_mass_algebraic vs all later work. Every claim holds (no-saddle->algebraic->E_8;
parity anchor strange=13; isospin via muon I_3; O(1) masses; I_3-derivation open). The E_6 work CONFIRMS
the anchor is genuinely E_8-level (exp 13 beyond E_6 range), so 'mass/anchor = E_8' is right, not superseded.
E_6 reframes the INTERPRETATION (E_8 = lepton bridge/global approx; unused 23,29 = over-description) but
contradicts NO claim. RECOMMENDATION: leave committed remarks as-is; E_6/generation-as-twist is a deeper,
still-note-only layer -- when it firms up, add a forward POINTER ('generation has a deeper per-strand E_6
realisation; E_8 here is the lepton bridge'), not a rewrite. Do NOT commit the E_6 picture to papers yet.

## 17. UP-SECTOR MASSES (up_sector_problem.py): REAL failure, STRUCTURAL cause (2026-09-04)
Scale-invariant mass-RATIO test (ratios are ~RG-invariant, so this bypasses the scale ambiguity):
  DOWN ratios good: m_s/m_d pred/PDG 0.90, m_b/m_s 1.07 (~10%).
  UP ratios FAIL: m_c/m_u 0.50, m_t/m_c 5.78 -- off ~2x in OPPOSITE directions (u->c gap too small,
  c->t too big). Scale-invariant => the up-sector failure is REAL, not a running/scale artefact.
CAN'T be patched in E_8: required up-type exponents ~{18.4, 9.4, 1.95}; 9.4 is NOT an E_8 exponent
(odd: 7,11,13,...), so charm has no clean exponent giving its ratio. Not a fitting error -- an inadequacy.
STRUCTURAL CAUSE: the down sector is anchored to the CHARGED LEPTON (strange=muon, both I_3=-1/2), which
calibrates the down tower. Up-type is I_3=+1/2; the only same-isospin lepton is the NEUTRINO (~meV,
~massless). So the up sector has NO massive lepton anchor -- a charged-lepton-referenced formula
STRUCTURALLY anchors down and ORPHANS up. The up-sector failure is the shadow of neutrinos being light.
=> fixes collapse: (a) re-fit E_8 exponents = DEAD (no exponent for charm); (c) top-special = PARTIAL
(u,c still off by ~3 in n without top); (b) reference the neutrino = RIGHT isospin but neutrinos too light
to anchor a MeV-GeV tower -> needs a genuinely NEW mechanism. This IS Paper XVI's 'generation-dependent
correction, not yet identified', now given a REASON: it is the isospin sector with no massive lepton partner.
STATUS: up-sector masses = real, STRUCTURAL open problem. Framework predicts down (anchored), not up.

## 18. UP-SECTOR REFRAME (user, 2026-09-04): formation energy + the bare top -- TOP SOLVED, u,c open
User insight: (i) a quark can't exist free in the 2I vacuum without energy input; it forms in collision
T(2,2)->T(2,3) (Paper XVIII production thresholds); mass = FORMATION energy (kinetic + barrier), not bare mass.
(ii) the TOP decays weakly FASTER than QCD hadronisation, so it never forms a bound T(2,3) -- doesn't exist
as other particles do.
QUANTIFIED (up_sector code): tau_top/tau_QCD = 0.15 (top decays 6.8x faster than it can hadronise) -- the ONLY
quark that never hadronises (all others, incl. b~1e-12 s, do). 
CLEAN RESULT (top): E_8 predicts the FULL T(2,3) formation m_t=m_tau*exp(15pi/9)=334 GeV; measured 173 GeV =
v/sqrt2 (EW scale, y_t~1). The top decays before COMPLETING the trefoil -> stops at the EW scale, never reaches
full formation -> pred/meas=1.93. The top's 2x 'over-prediction' is the INCOMPLETE-FORMATION deficit, not a
formula failure. Turns the worst up-sector miss into a feature.
STRUCTURAL: down {d,s,b} ALL confined -> homogeneous -> fits; up {u,c,t} has u,c confined but top BARE ->
inhomogeneous -> no single formula should fit. The up-sector 'failure' is largely the top being a different object.
HONEST RESIDUAL: even excluding the bare top, confined u,c still miss ~2x (m_c/m_u pred 293 vs 590; n_c-n_u req
3 vs formula 1). So TOP = explained (bare/incomplete formation); u,c = still open. The formation-energy FRAME is
right; the u->c gap for the two confined up quarks is the remaining piece.
NEXT: compute the T(2,2)->T(2,3) FORMATION energy (E(Q_H=3)-E(Q_H=2) + 8pi Chern-Simons barrier + kinetic) =
Paper XVIII open O6, and test whether formation energy (not bare mass) closes the confined-u,c gap.
This supersedes sec.17's 'structurally unanchored' framing: the mechanism is bare-top + formation-energy.

## 19. >>> RESUME HERE (2026-09-04) -- QGP jamming model, pinning R (item 1) <<<
ARC SO FAR (this note, secs 1-18): generation-as-twist DERIVED at WZW/McKay level -- generation = per-strand
E_6 twist (2T<->E_6<->SU(3)_1), E_8 = lepton bridge/global approx, mirror=isospin doublets, framing mod-3.
Committed to papers: P17:rem:quark_mass_algebraic (E_8 mass/anchor), P16 pointer -- SOUND, leave as-is (E_6
is deeper note-only layer). Bishop holonomy -63deg = golden-determined (r0=sqrt2/phi) transcendental, root-
sibling of the framings, NOT their quantisation. UP-SECTOR is the live open problem.

UP-SECTOR STATE: down {d,s,b} predicted (anchored strange=muon, ratios right up to uniform ~1.2 factor).
UP {u,c,t} fails. Cause (user): (i) TOP is bare -- decays weakly 6.8x faster than QCD hadronisation, never
forms T(2,3); E_8 predicts full-formation 334 GeV, measured 173=v/sqrt2 (incomplete formation). TOP SOLVED.
(ii) confined u,c still miss ~2x (m_c/m_u pred 293 vs 590) -- OPEN. Static formation-energy route BLOCKED
(Q_H=3 saddle finder escapes sector = the instability). Algebraic route = the failing E_8 route.

CURRENT TASK -- the QGP JAMMING analytical model (user's reframe; define analytics BEFORE numbers):
Quarks are unstable isolated (no saddle) but 'stable' in the gluon soup because (a) DECAY-BLOCK: un-threading
T(2,3)->T(2,2) (baryon->lepton) needs a strand to sweep ~R to escape; in the dense soup neighbours occupy that
room, so decay is blocked when inter-quark spacing d < R (topological JAM, not an energy barrier); (b) continuous
RECOMBINATION (strand exchange) pins the population in Q_H=3. Deconfinement = jamming release at d ~ R.
Mass = the JAMMED-STRAND formation energy (well-defined BECAUSE the jam stabilises it), not the isolated saddle.
THREE THINGS TO DEFINE (user's plan): 1. pin R (un-threading scale) from C*=2.5062, r0=sqrt2/phi, R0=3 -- DOING
NOW (item 2 needs it). 2. the jammed formation energy = E(T(2,3) in jam) - E(T(2,2)) [ΔQ=1 step WITH the jam,
not the blocked isolated saddle], vs master formula m ∝ phi^{2N} exp(4pi ΔQ). 3. recombination-vs-unthreading
rate ratio vs d/R -- may fall out of 1+2. Then test: does the JAMMED formation energy close the confined u,c gap.
ITEM 1 DONE -- R PINNED (e6... inline calc, src_paper16): R = R0 + r0 = 3 + sqrt2/phi = 3.874 (trefoil OUTER
extent = the un-threading/jamming scale; a neighbour within R blocks a strand escaping T(2,3)->T(2,2)).
Framework-determined: R0=3 (from 3R0=N3^2=9, N3=3=baryon Q_group), r0=sqrt2/phi=R0/C2*. Dimensionless
R/r0 = 1 + 3phi/sqrt2 = 4.43. Jam condition: d<R jammed/deconfined, d>R confined; deconfinement at d~R,
density n~1/R^3 (~1/fm^3~Lambda_QCD^3 physically, R~1fm). PRIMARY R=R0+r0=3.87 (global knot size);
ALTERNATIVE if un-threading is a LOCAL reconnection: strand gap=1.52 or tube 1/C*=0.399 (item 2 decides which).
Other scales: Rg=3.12, R_xy=R0-r0=2.126.
ITEM 2 -- HANDLE (1) DONE (cheap sigma*R test, 2026-09-04): KEY correction -- E_baryon = sigma_tot*3R0
ENTIRELY (paper16 l.1393-1429: sigma_int*3R0=profile, sigma_ext*3R0=the 4.4%/1.044). So sigma_tot is
FIXED = 646.5/9 = 71.83 MeV per unit length, NOT free. Missing E_cup=291.8 MeV (the 31% paper XVI
couldn't source geometrically, l.1586-94 'QCD binding via string network') <-> definite extra length
L_cup = 291.8/71.83 = 4.062 units. JAM CELL R=R0+r0=3.874 gives sigma*R = 278.3 MeV = 95.4% of E_cup.
Parameter-free proton: m_p ~ sigma*(3R0 + R) = sigma*(4R0+r0) = 924.8 MeV vs 938.3 (-1.44%). Y-junction
(3 arms) + ONE jam cell. NOTHING fitted. Discarded circular candidate 3R0*(mp-Ebar)/Ebar (tautology).
Did NOT chase the 0.19-unit residual (1/2 tube-radius or L_arc/pi would fit = fishing, not claimed).
USER REFRAME (2026-09-04, important): the 1.4% deficit is EXPECTED and REASSURING, not a failure -- if
the jam recovered the FULL proton (938) it would imply free protons EXIST in the QGP, which we DON'T see.
Constrained overlap + winding-switching (recombination) during a jam CANNOT complete a full isolated
proton; the last ~1.4% is the de-jamming/completion (localisation) energy, paid only when the jam RELEASES
(d>R, hadronisation). So sigma*(3R0+R)=924.8 = the JAMMED energy; m_p=938 = the RELEASED energy. This
SHARPENS handle (2): the tight-box flow should land ~924.8 (jammed, ~98.6%), and should NOT reach 938 --
full recovery would be the worry. Deficit magnitude 13.5 MeV = light-quark/current-mass ORDER (proton
current content 2m_u+m_d~9 MeV; not an exact hit, noted not claimed) -- consistent with completion =
the current/localisation piece. A perfect fit would be a red flag (CLAUDE.md 3/6); 98.6% with a physical
reason for the residual is stronger. HANDLE (2) NEXT: constrained flow in tight box 2R~7.7, background it.
HANDLE (2) LAUNCHED (2026-09-04): src_paper16/qh3_trefoil_solver_jam.py (JAM variant of
qh3_trefoil_solver_3d.py: numpy nearest-point instead of scipy KDTree; box assertion relaxed to
allow tight box; periodic BC = lattice of jammed baryons). Run: --N 36 --h 0.2152 (box=7.747=2R,
marginal jam d=R) --C_star 1.5 --steps 200000, outdir src_paper16/jam_marginal, bg ID bdktedf9b,
~30 min. EARLY (200 steps): K falling, K/J2a=3.31 vs 2phi=3.236, K/J4 descending toward phi^6=17.94,
J4/J4_init=1.015 (TOPOLOGY HELD -- opposite of the established large-box escape, CLAUDE.md 4). PASS =
J4/J4_init~1 AND K converges (jam stabilises the saddle = un-blocks the 3D saddle finder). Converged
energy should map to ~924.8 MeV (jammed, 98.6% of m_p), NOT 938 (full recovery would be the WORRY).
K/J4->phi^6 = genuine Bogomolny saddle. CONTROL (same-solver large box escaping) deferred (est. hrs);
large-box escape is established in papers. MeV conversion (K natural units -> vs 646.5/924.8) TBD from
converged K vs isolated-profile K. ENV: base conda has torch+numpy NO scipy; torch_intel is EMPTY
(user's 'conda activate torch_intel' note was wrong -- verified). See notes/pyenv.md.
HANDLE (2) RESULT = BLOCKED on the open saddle-finder (2026-09-04): the run ESCAPED -- J4/J4init
1.46->0.42 monotonically (in-place Hopf DILUTION), K/J2a drifting off 2phi. This is the KNOWN v6-v8
failure of the UNCONSTRAINED K_fb flow (I mistakenly used the v4-lineage solver). Repo history
(qh3_trefoil_solver_3d_v11.py docstring, gradient_flow_constrained.py): unconstrained descent ALWAYS
dilutes (v6-v8) or spikes (v9); topology needs post-step slerp projection (gradient_flow_constrained.py,
two-phase) or Metropolis (v11); the saddle-finder is OPEN, paper16 'a literal prerequisite'. Tight box
does NOT help -- dilution is in-place on T^3, not spreading. So there is NO stable Q_H=3 saddle to find
(dilution IS the instability = the thread's finding). Killed PID 29938. MISTAKES: (a) double-backgrounded
(nohup& inside run_in_background) -> lost harness tracking + false 'completed' notification; (b) read
'topology holding' at 200 steps too early (turned over by 2000) -- CLAUDE.md 9. HANDLE (1) analytic result
(sigma*R=924.8=98.6% m_p, deficit-expected reframe) STANDS independent of the numerical saddle.
DECISION POINT: (a) accept handle-1 + move to generation-dependence / item 3 (the real u,c prize) --
RECOMMENDED; or (b) invest in gradient_flow_constrained.py (topology-protected) tight-vs-large ENERGY
comparison (the one clean numerical jam test, but nontrivial -- constrained solver has its own
convergence issues, saddle-finder open). Do NOT re-run the unconstrained jam solver.

### WHEN REVISITING HANDLE (2) -- checklist (2026-09-04, user + model)
Goal: measure the JAM COMPRESSION energy = E(confined, topology-held) - E(isolated, topology-held),
and its dependence on inter-baryon spacing d, WITHOUT the escape/dilution confounding it.
 1. SOLVER: use the TOPOLOGY-PROTECTED flow, never the unconstrained one. Options in src_paper16:
    gradient_flow_constrained.py (post-step slerp projection to delta_max + two-phase schedule) or
    qh3_trefoil_solver_3d_v11.py (Metropolis annealing). Unconstrained K_fb ALWAYS dilutes(v6-8)/spikes(v9).
 2. MULTIPLE BARYONS (user): put 2+ distinct trefoils in the box at controlled separation d -- NOT one
    trefoil with periodic self-images (at box=2R a single trefoil touches its own image and RECONNECTS,
    an artefact). Real sharing/overlap is between distinct baryons. Models the QGP lattice properly.
 3. TIGHTER THAN THE MULTIPLES (user): sweep d downward from d=R (marginal) to d<R (denser jam) -- i.e.
    box radii TIGHTER than integer multiples of the baryon size. Predict E_jam(d) RISES as d shrinks
    (compression -> completion energy grows -> in-medium mass DROP, a real QGP signature). d=R should
    give ~924.8 MeV (98.6% m_p); denser -> lower jammed mass, larger completion.
 4. FIXED-TOPOLOGY ENERGY: the deliverable is the energy DIFFERENCE at held J4 (large vs tight box), the
    compression. Both runs must hold topology (that's what the constrained solver buys).
 5. CONTROL: isolated large box, SAME constrained solver -- confirms the constraint (not the box) holds
    topology, isolates tooling from physics (CLAUDE.md 4). Establishes the E_isolated baseline.
 6. MeV CONVERSION: pin natural-units K -> MeV once via the master-formula profile (E_baryon=646.5 <->
    E_profile*, C*=2.5062). Then E_jam(d) is directly comparable to E_cup=291.8.
 7. NUMERICS: >=3 pts across tube (tube_r=1/C*); conservative dt/max_dt (my jam run's dt-adaptation
    thrashed -- big steps drove the dilution faster). Box slightly > 2R for single-baryon runs to avoid
    self-image reconnection; explicit multi-baryon preferred.
 8. THEN generation-dependence: add the ribbon twist to E_jam(d) and test the u,c ratio (option-a work).
 9. PROCESS: launch with run_in_background=true ALONE (no nohup&) so the harness tracks it and gives a
    real completion signal; nohup& inside it detaches and loses tracking (2026-09-04 mistake).

## 20. OPTION (a) -- generation-dependence: sigma*R is UNIVERSAL, gen lives in the TOWER (2026-09-04)
KEY FINDING (paper16 l.1596-1605, [E]): sigma*R is GENERATION-UNIVERSAL. All quarks are the SAME trefoil
CURVE T(2,3); generation = ribbon TWIST (framing), a ribbon not a curve property -> R=R0+r0 and sigma_tot
are curve properties, identical across generations. So the twist does NOT change the confinement length/
tension; sigma*R cannot be the u,c lever. Paper's OWN words: 'Current quark masses are lighter by factors
phi^{-2n_q^(g)} from the generation-dependent tower level n_q^(g)... exactly as leptons (Paper XII). The
derivation of n_q^base from a Q_H=3 thermodynamic matching (analogue of n_e=20) is the remaining OPEN
calculation.' So the up-sector RATIO problem = a TOWER problem, and the open piece = the tower BASE n_q^base.
TOWER TEMPLATE (paper17 l.812-882, P17:eq:mass_formula, Paper VI): m^(g) = Lcond * phi^{-2n} *
exp(-2pi * Q_group * T_g). Lepton base n_e = 20 = 2*Q_group^(l) = 2*10, from CMB energy identity Efb/V=rCMB
(Paper XII Thm E; the type-(i) 'CMB normalisation' integer, NOT the irrational absolute level ~50.25).
=> NATURAL BARYON ANALOGUE (candidate, [N] to verify): n_q^base = 2*Q_group^(baryon) = 2*3 = 6 (= the phi^6
already in the master formula P16:eq:mq_leading). Q_group triple {baryon3, nu6, lepton10} (claude-hopfion 8).
LIGHT/HEAVY regime [I]: constituent scale m_q^(0)=215.5 MeV. LIGHT quarks u,d,s (<215) sit BELOW it
(tower-suppressed phi^{-2n}, n>0); HEAVY c,b,t (>215) sit ABOVE (n<0, a different regime -- not 'light
constituent quarks'). Both up and down sectors cross this scale, so light/heavy alone != why down works/up
fails. NEXT: derive n_q^base via the CMB-energy-identity analogue for Q_H=3 (template gives 2*Q_group=6),
then build the quark spectrum m(g)=215.5*phi^{-2(n_q^(g))}*exp(...) and test the up-sector u,c ratio.
Down works via E8 route (up_sector_problem.py); recheck with the tower parametrisation + the E6 twist->T_g.
DIAGNOSTIC (2026-09-04, inline calc): required tower level n from m=215.5*phi^-2n:
  u+4.78 c-1.85 t-6.95 (up) | d+3.98 s+0.87 b-3.08 (down). +=below constituent(light), -=above(heavy).
  Mean gen increment: UP -5.87/gen, DOWN -3.53/gen -> UP tower ~1.6x STEEPER (ratio 1.66; individual
  increments scatter 27-30% so NOT claiming phi=1.618 or 5/3=1.667 -- numerology trap, CLAUDE.md 3/6).
  Isospin split n_up-n_down GROWS with gen: +0.80, -2.72, -3.87 (g=1,2,3). Light/heavy: u,d,s below
  215.5 (n>0), c,b,t above (n<0). CONCLUSION: up-sector problem is NOT confinement (sigma*R universal),
  NOT the base (n_q^base=6 same for both isospins) -- it is an ISOSPIN-DEPENDENT TOWER STEEPNESS
  (I_3=+1/2 up climbs faster than I_3=-1/2 down). => the lever = how isospin couples to the tower
  increment (T_g or effective Q_group split by I_3). Ties to the existing isospin-split thread (muon
  resonance, Z_2 z-asymmetry). NEXT FORK: derive the I_3-dependence of T_g / Q_group in the tower.
GEOMETRY CONSTANTS: R0=3, r0=sqrt(2)/phi=0.874032, C*=2.5062 (tube radius 1/C*=0.399), R_xy=R0-r0=2.126
(midpoint/crossing radius). Trefoil param Gamma(t)=((R0+r0 cos3t)cos2t, (R0+r0 cos3t)sin2t, r0 sin3t).
Scripts: src_paper16/ (e8_strange_anchor, up_sector_problem, e6_*, bishop_frame_v2), src_paper18/ (e6_*,
cherenkov_*, breaking_exponent_bijection, sweep_beyond_trefoil). Dynamic integrator: src_paper16/
qh_dynamic_integrator_v2.py (kick+Langevin, T(2,2) vacuum -> kick -> evolve; matches kinetic-formation but
gives a TRANSIENT shower/fragmentation, not a clean mass). Viz: hopfion_trefoil.html, sweep_beyond_trefoil.py.

## 21. n_q^base DERIVED via the CMB identity; u,c gap = isospin tower steepness (2026-09-21)
New machine, `python3.11` env (numpy/scipy/torch/sympy; see [[pyenv]]).
Script: `src_paper16/nq_base_tower.py`. Two results.

(A) **n_q^base = 6 is DERIVED at the CMB-arithmetic level [N->E for the arithmetic half].**
The lepton n_e=20 comes from the chain (Paper XII Thm A,E): rho_CMB = Q_group * Lcond^4 with
Lcond^(sector) = T_CMB (pi^2/(15 Q_group))^{1/4}, so rho_CMB/Lcond^4 = 15 Q_group/15 = Q_group
(an integer), and n_base = 2 Q_group. Paper XVII ALREADY generalises this from the lepton
(Q_group=10, Lcond=T_CMB(pi^2/150)^{1/4}) to the neutrino (P17:prop:qh1_cmb: Q_group=6,
Lcond=T_CMB(pi^2/90)^{1/4}, rho_CMB=6 Lcond^4). The BARYON is the SAME construction:
  **Lcond^(baryon) = T_CMB (pi^2/45)^{1/4} = 1.607e-4 eV, rho_CMB/Lcond^4 = 45/15 = 3 =
  Q_group^(baryon), => n_q^base = 2*3 = 6.** VERIFIED numerically (ratio 3.000000). Calibration:
  the lepton chain m_e = Lcond^(l) phi^20 e^{4pi-3/400} recovers 0.51212 MeV (0.2%).
So the "n_q^base = 2 Q_group = 6" of sec.20 is NOT just an analogy -- it is the arithmetic
(CMB thermodynamic) half of the Paper XII/XVII derivation, run for Q_group=3.
  **CAVEAT [O] -- the saddle half is BLOCKED.** The FULL two-route equivalence (Thm E) also
  needs the profile-normalisation route phi^9 sqrt(Ja[f*] J4[f*]) = Q_group, which integrates
  over the SADDLE PROFILE f*. Q_H=3 has NO stable saddle (the mass is algebraic -- the whole
  reason for the E_8 route, P17:rem:quark_mass_algebraic). So the profile half cannot be
  evaluated for the baryon; n_q^base=6 rests on the CMB-arithmetic route alone. This is the
  SAME no-saddle obstruction, reappearing -- it does not break the CMB derivation (that half is
  saddle-free) but it means the baryon lacks the second, independent confirmation the lepton and
  neutrino have. Sharpens sec.20's "remaining OPEN calculation" to: DONE (arithmetic), the
  saddle-route cross-check is what's structurally unavailable.

(B) **The u,c gap is an ISOSPIN-DEPENDENT TOWER STEEPNESS -- E_6 twist does NOT explain it.**
Empirical tower levels from m = 215.5 phi^{-2n} (reproduces sec.20 exactly): u+4.78 d+3.98
s+0.87 c-1.84 b-3.08 t-6.95. Findings:
  - DOWN {d,s,b} is a REGULAR tower: increments -3.11, -3.95 (~-3.5/gen). Works (anchored).
  - UP {u,c,t} climbs STEEPER: increments -6.63, -5.10. Confined slope ratio (u->c)/(d->s)=2.13
    (top excluded, bare/incomplete-formation sec.18). Isospin split n_up-n_down GROWS with gen:
    +0.80, -2.71, -3.87.
  - **The E_6 twist does NOT carry it**: the twist isospin-split {2,4,10} (gen 1,2,3) gives
    split/twist ratios {0.40, -0.68, -0.39} -- NOT constant. So the generation-as-twist geometry
    (secs 11-14) labels generations but does NOT generate the I_3-dependent steepness. Clean
    negative: twist gives labels, isospin steepness is a separate coupling (consistent w/ sec.8).
  - The u,c number: measured m_c/m_u = 588; the isospin-blind null (up stepped by the DOWN
    per-gen increment) = 20; up needs an EXTRA Delta n = -3.51 tower levels gen1->2. That extra
    -3.51 ~ the down mean increment (-3.53): the confined up tower climbs at ~2x the down rate
    gen1->2 (2.13, near 2 -- flagged NOT claimed as exactly 2, CLAUDE.md 3).
  CONFIRMS sec.17/sec.20: not the base (6 is isospin-blind), not confinement (sigma*R universal,
  sec.20), not the twist -- a genuinely NEW I_3->tower-increment coupling is required. The up
  sector has no massive lepton anchor (I_3=+1/2 partner = the ~meV neutrino), so nothing
  calibrates its steepness. LEVER still open: derive the I_3-dependence of the effective tower
  increment (T_g or Q_group split by I_3). NOT the phi^{-2n} base and NOT sigma*R.

## 22. THE I_3 -> TOWER LEVER (2026-09-21): worked, not closed -- charge is the lead, degeneracy the wall
Script: `src_paper16/i3_tower_lever.py`. Two independent routes agree on the shape; the lever is
identified but underdetermined. Tags per CLAUDE.md 6.

STRUCTURE OF THE GAP (both routes):
  - Lumped tower (m=215.5 phi^{-2n}): DOWN per-gen increments -3.11,-3.95; UP -6.63,-5.10.
    Confined slope ratio (u->c)/(d->s) = **2.13** (the ONLY top-free comparison).
  - E_8 route (unambiguous -- no n-T degeneracy): DOWN deficits (req-assigned) -0.66,-0.35,-0.55
    ~ uniform -0.5 = ln(1.2)*9/pi (the known uniform running factor, pred/meas~1.2). UP deficits
    +1.13,+3.12,+0.00. Confined up EXTRA beyond the uniform: u +1.65, c +3.64 -> ratio **2.21**.
  => CONVERGENT: the confined up sector needs an EXTRA generation-growing suppression of ~2.1-2.2
    per generation relative to down. (BONUS: E_8 top full-formation deficit = 0.00, pred/meas=1.00
    EXACTLY -- clean re-confirmation of sec.18: top stops at v/sqrt2, full T(2,3) would be 334 GeV.)

CANDIDATES TESTED:
  - H_Q (steepness ~ |electric charge|, ratio |2/3|/|1/3|=2): matches the steepness (2.13, 2.21 vs
    2.0) -- the ONLY candidate with the right shape. Charge is a real, mass-INDEPENDENT framework
    quantum number (writhe/self-linking, alpha=360/phi^2), and steepness~charge makes the up-down
    split GROW with generation automatically (steeper up tower, common-ish gen-1 start -> diverges),
    reproducing u<d at gen1 but c>s, t>b. BEST LEAD.
  - H_T (split ~ T-matrix phase T_g={5/24,1/8,0.075}): FAILS -- T_g DECREASES with gen while |split|
    GROWS (wrong trend); split/T_g not constant (+3.85,-21.7,-51.6).
  - H_tw (split ~ E_6 twist-split {2,4,10}): FAILS -- twist-split all +, but n-split flips sign
    (+0.80,-2.71,-3.87); split/tw not constant.

WHY H_Q IS NOT CLAIMED (honest caveats, CLAUDE.md 3/6):
  (1) ONE clean data point. gen2->3 is top-contaminated (bare top) in the lumped route and top-as-
    input in the E_8 route; so the "law" rests on the single u,c-vs-d,s comparison. 2.13~2.21 is not
    razor-clean and one point can't distinguish charge-2 from anything near 2.
  (2) Exp sensitivity: the 6% slope miss (2.13 vs 2) becomes a ~50% miss in m_c/m_u (pred 400 vs 588).
  (3) DECISIVE against the naive mechanism: the natural physical origin (EM self-energy in the
    T(2,2)->T(2,3) formation energy) scales as CHARGE-SQUARED -> ratio (2/3)^2/(1/3)^2 = **4**, NOT the
    observed ~2. So charge-LINEAR steepness has no ready mechanism; charge-SQUARED has a mechanism but
    the wrong number. Either way not settled. The charge match may be coincidence.
  (4) T-vs-n DEGENERACY (the core wall): m=Lcond phi^{-2n} exp(-2pi*3*T) fixes only 2n lnphi + 6pi T
    per mass. Fixing T from the E_6 twist (T=Tw/3) forces e.g. n_u=-1.75 vs lumped +4.78, n_d=+10.5 vs
    +3.98 -- the same masses with wildly different n. Isospin CANNOT be uniquely assigned to n vs T
    without an independent per-quark handle. This is why the lumped-vs-twist pictures disagree.

STATUS: lever IDENTIFIED (extra gen-growing up-suppression ~2/gen; charge is the only right-shaped
lead), NOT closed. NEXT (the one calc that could settle it): compute the actual T(2,2)->T(2,3)
FORMATION-energy EM contribution (Paper XVIII open O6, sec.18) and read off whether it scales as
charge^1 (ratio 2, would confirm H_Q) or charge^2 (ratio 4, would kill it and demand a new term).
That formation-energy computation is now the single highest-value step -- it simultaneously tests
H_Q, breaks the n-T degeneracy (gives an independent mass piece), and closes sec.18's O6.

## 23. LEVER REFRAME (user, 2026-09-21): isospin asymmetry = GROUP asymmetry (E_6 vs E_8), top = sector transition
User pushback on sec.22's E_8-heavy framing, both points CORRECT per the notes:
  (i) The quark sector is E_6-NATIVE (sec.11: E_6 has exactly 6 Coxeter exponents = 6 quarks;
      Paper XVI l.581 m_q = m_ell[2I/E_8]*colour[2T/E_6/SU(3)_1] -- E_8 is only the LEPTON-inherited
      piece). My sec.22 i3_tower_lever leaned on the E_8 route, i.e. the lepton BRIDGE, not the
      quark's own group.
  (ii) TOP = the SECTOR TRANSITION (sec.18: top stops at v/sqrt2, y_t~1, never completes T(2,3)).
      Sharpened: top sits at the fermion<->EW/gauge boundary (Q_H=3 -> Q_H=4 gauge composites,
      claude-hopfion 8), NOT a normal up-tower member. So the up "tower" is u,c ONLY; top is the
      upper boundary.

REFRAME -- this IMPROVES the lever (supersedes sec.22's H_Q as the lead):
The isospin asymmetry is NOT a new I_3 coupling. It is a GROUP asymmetry set by whether a massive
lepton anchor exists:
  - DOWN {d,s,b}: anchored by strange=muon -> uses the **E_8 lepton bridge** (h(E_8)=30), gentler tower.
  - UP {u,c}: no massive lepton partner (neutrino too light) -> **E_6-NATIVE** (h(E_6)=12), the quark's
    OWN group, STEEPER tower. Top = transition (v/sqrt2), the boundary, not a member.
PREDICTION (structural, no new coupling): phase unit ~ 1/h, so up/down tower-steepness ratio ~
h(E_8)/h(E_6) = 30/12 = **2.5**. OBSERVED (clean top-free comparison) = **2.13**. Same ballpark, ~15%
off -- comparable to H_Q's 6% (2.13 vs 2.0) BUT with real mechanism (only down has a lepton to
anchor to -> E_8; up is native -> E_6), needing NO extra term. This is now the LEAD over H_Q.
  - Why better than H_Q (charge): H_Q had no mechanism (charge-linear) or the wrong number
    (charge^2 -> 4). H_h (group asymmetry) is exactly the E_6/E_8 split the notes already established
    (secs 11-16) and explains the ORIGIN of the isospin dependence.
  - Same honest caveats survive: ONE clean data point (top-free), exp sensitivity, and the
    n-vs-T degeneracy. 2.13 vs 2.5 not razor-clean.

OUTSTANDING CALC (unchanged priority, now sharper): build the **E_6-NATIVE mass formula** -- the
analogue of the E_8 route (m_q/m_ref = exp(n_q * pi Q/(k h)), n_q = Base'_g - 2 m_E6) with h(E_6)=12,
the E_6 exponents {u:7,c:8,t:11 / d:5,s:4,b:1}, and the up sector anchored at the TOP-TRANSITION
(v/sqrt2) instead of a lepton -- and test whether it reproduces m_u, m_c (note sec.11 l.160 item ii,
NEVER DONE). NOT built here: guessing the phase unit / Base offset without grounding = numerology.
The formation-energy EM piece (sec.22, Paper XVIII O6) remains the other route; the two are
complementary (group-native spectrum vs dynamical formation energy).

## 24. E_6-NATIVE MASS FORMULA test (2026-09-21): phase-ratio EXPLAINS the steepness; spectrum does NOT fall out
Script: `src_paper16/e6_native_mass.py`. Grounding from papers: T_g=(6-g)/(8(g+2))={5/24,1/8,3/40}
(P17 l.601,627, same across sectors); Base_g=180 T_g+7/2={41,26,17}, 180=6 h(E_8); E_8 phase
pi/9 = pi Q/(k h) with Q=10,k=3,h(E_8)=30. E_6-native swaps in the QUARK's own numbers: Q=Q_group^
baryon=3, k=1 (SU(3)_1=(E_6)_1, c-match at P16 l.468), h(E_6)=12.

RESULT (1) -- WHY up is steeper, GROUNDED: E_6-native phase = pi*3/(1*12) = **pi/4**. Ratio to the
E_8-bridge phase pi/9 is **(pi/4)/(pi/9) = 9/4 = 2.25**. Observed up/down tower steepness = 2.13
(sec.22) -> match to **6%**. So the up sector being E_6-NATIVE (larger phase unit, quark's own group)
vs the down sector being E_8-lepton-BRIDGED gives the steepness asymmetry with the right magnitude
and NO new coupling. This is the best-grounded form of the sec.23 group-asymmetry lead (better than
the raw h-ratio 2.5 or charge 2.0): the actual phase unit pi Q/(k h) uses ALL three quark numbers.

RESULT (2) -- E_8 down sanity (confirms established): pred/meas d,s,b = {1.26,1.13,1.21} (uniform
~1.2 running); u,c = {0.67,0.34} (scattered, unanchored); top full-formation = 1.000 EXACTLY.

RESULT (3) -- but the NAIVE E_6-native formula does NOT reproduce the SPECTRUM. E_6 mirror pairs =
isospin doublets, so the naive form gives m_up/m_dn = exp(-2(m^up_E6 - m^dn_E6)*phase). The phase each
doublet REQUIRES is NOT constant and NOT pi/4: u/d -0.19, c/s +0.33, t/b +0.19-0.22. The SIGN even
flips (u lighter than d at gen1; c,t heavier at gen2,3) -- the naive single-phase mirror cannot
capture the gen1<->gen2 up-down inversion.

SYNTHESIS: E_6-native cleanly explains the **generation-steepness TREND** (2.25~2.13) -- that axis is
the phase unit pi Q/(k h). It does NOT by itself give the **isospin doublet splitting** (within-gen
up-vs-down, which flips sign gen1->gen2). These are ORTHOGONAL: generation steepness = E_6 phase unit
(SOLVED in shape); isospin splitting = the sign-flipping within-doublet piece (STILL the hard open
part, ties to the formation-energy EM piece sec.22 and the sec.4 "two isospin mechanisms").
OUTSTANDING to make E_6-native a real spectrum: pin (a) the reference mass (top-transition v/sqrt2 per
user, or the sector scale), (b) the Base' offset (E_8's 7/2 analogue), (c) the T_g source (lepton
phases vs the SU(3)_1 twist framing) -- NOT guessed here (would be numerology). The isospin sign-flip
is the deeper physics and is not an E_6-vs-E_8 question.

## 25. ISOSPIN SIGN-FLIP via the formation-energy EM piece (2026-09-21): EM is the WRONG SIGN
Script: `src_paper16/isospin_signflip_em.py`. User asked to work the sign-flip (u<d at gen1, but
c>s, t>b) via the EM formation energy. RESULT: **EM is not the mechanism -- it has the wrong sign.**

DECOMPOSITION of the flip. R(g)=m_up/m_down = {0.46, 13.60, 41.3} (crosses 1 between g1,g2).
ln R = {-0.77, +2.61, +3.72}. Splits cleanly into:
  - MEAN ln R = 1.85 -> M_up/M_down = 6.38. This is the ISOSPIN SCALE, and it is Paper XVIII's
    P18:prop:isospin_mass_ratio (writhe/out-of-plane tangent geometry, M_down/M_up=0.189=1/5.29,
    predicts ln=1.67 vs 1.85, 21% off). [existing framework result]
  - SPREAD about the mean {-2.62,+0.76,+1.87} = the GENERATION-dependence = up climbs steeper than
    down = the E_6-native(up)/E_8-bridge(down) phase-ratio 2.25 (sec.24). d(lnR)/gen = up_incr -
    down_incr = 3.38.

EM SIGN TEST (the requested piece). gen1 needs m_d - m_u = +2.51 MeV (DOWN heavier). EM self-energy
of a charged soliton dm_EM ~ (alpha/2) Q^2/r is ALWAYS POSITIVE and larger for the higher charge:
up(2/3) vs down(1/3) -> EM makes UP heavier by +0.1 to +0.5 MeV (r=0.5-2 fm). So EM:
  - has the WRONG SIGN (makes up heavier; gen1 needs down heavier), and
  - is TOO SMALL (~0.5 MeV vs the needed 2.5 MeV).
=> **EM CANNOT be the gen1 seed.** This is exactly the SM situation (the bare m_d>m_u dominates the
EM splitting -> the neutron is heavier than the proton). It also DEFINITIVELY kills the sec.22 H_Q
charge/charge^2 lead: charge^2 -> up heavier, the opposite of the observed gen1 ordering.

BOTH known "up-heavier" mechanisms (EM ~+Q^2; P18 out-of-plane M_down/M_up=0.189) have the WRONG SIGN
at gen1. So the flip's essence = why up is LIGHT at gen1, opposing both.

RESOLUTION (no EM): the ANCHORING ASYMMETRY. Down = E_8 lepton-bridge, anchored LOW at strange=muon
(gen2), gentle tower. Up = E_6-native, anchored HIGH at the TOP-transition v/sqrt2 (gen3), STEEP
tower. Because up is anchored at the heavy top and climbs a steep tower up to gen1, it OVERSHOOTS to
higher n (lighter) than down (anchored at the light muon, gentle tower): n_u=4.78 > n_d=3.98 -> u
lighter. So the sign-flip = E_6/E_8 steepness (sec.24) + top-as-transition anchoring (sec.23), TWO
already-established pieces, NO new EM term. The flip is a STRUCTURAL consequence of the two towers
being anchored at opposite ends (up high/top, down low/muon) with different (E_6 vs E_8) steepness.

STATUS [CORRECTED in sec.26 -- read that]: sec.25 concluded "EM route CLOSED." That over-generalised
from the SELF-energy (Q^2, which IS wrong-sign/too-small) to ALL EM. The user pushed back (CLAUDE.md
9); redone in sec.26: the LINEAR charge*field coupling (-Q*Phi) is signed, right-sign, and reopens EM
as the gen1-seed flipper. Keep: the Q^2 self-energy is not the flip; the decomposition (mean=P18,
spread=steepness) stands. Retract: "EM is a red herring" -- it is not.

## 26. EM REOPENED (user pushback, 2026-09-21): the LINEAR charge coupling, not the self-energy
Script: `src_paper16/isospin_em_linear.py`. User: don't close EM so fast; it may factor into pinning
the E_6-native absolute tower. CORRECT -- sec.25 tested only the positive-definite Q^2 SELF-energy and
wrongly generalised. EM enters TWO ways:
  (A) SELF-energy ~ +c Q^2: positive-definite, up heavier, ~0.5 MeV. sec.25's sign kill applies HERE
      only -- not the flip.
  (B) LINEAR coupling ~ -Q*Phi (charge x a background field): SIGNED. **P18:conj:isospin explicitly
      invokes such an external field** ("under any external field distinguishing the two z-levels...
      selects which quark charge emerges"). So the linear piece is the FRAMEWORK-NATIVE EM term, and
      it was never tested in sec.25.

THE LINEAR TERM WORKS -- right sign, right structure:
  - SIGN: dm=-Q*Phi -> up(+2/3) shifts DOWN, down(-1/3) shifts UP; dm_down - dm_up = +Phi > 0 =>
    DOWN heavier. RIGHT SIGN for the gen1 seed (opposite to the Q^2 self-energy).
  - MAGNITUDE: Phi ~ alpha*Lambda_const = (1/137)*215.5 = 1.57 MeV; dm_down-dm_up=+1.57 vs observed
    m_d-m_u=2.51 (order-of-magnitude match, an O(1) geometric factor closes it).
  - GENERATION STRUCTURE (the decisive part): Phi is a fixed formation-scale field -> the shift is
    a CONSTANT ~1.5 MeV. Fractional effect: gen1 u 48.5%, d 11%; gen2 s 0.6%, c 0.08%; gen3 ~0%.
    So it FLIPS ONLY gen1 and leaves c>s, t>b untouched. EXACTLY the sign-flip structure.

REVISED PICTURE (supersedes sec.25's closure):
  m_measured = m_topological(E_6/E_8 tower + anchoring + P18 out-of-plane) + EM.
  - Topological tower: up-heavier at ALL gens (the "bare" ordering).
  - EM LINEAR: constant down+/up- by ~1.5 MeV -> flips ONLY gen1 -> observed sign-flip.
  So EM is NOT a red herring; it is the ~MeV charge-linear offset that MUST enter the E_6-native
  absolute-tower pinning (fit the topological tower to m_measured + Q*Phi, i.e. EM-subtract first).
  This is why EM factors into pinning the tower (user's point): at gen1 EM is ~50% of the up mass,
  so ignoring it would badly bias the tower fit at the light end.

UNDETERMINED (honest): Phi's SIGN (fixed by observation, as strange=muon fixed the down direction --
P18:conj:isospin says the field selects the charge, it doesn't derive which) and its exact MAGNITUDE
(the O(1) geometric factor on alpha*Lambda_const). But the STRUCTURE -- signed, gen-independent,
O(alpha*Lam)~MeV, flips only gen1 -- is correct and framework-native.
NEXT: fold Q*Phi into the E_6-native absolute-tower pinning (sec.24 outstanding) as a charge-linear
EM offset; then the O6 topological Delta E carries the generation/steepness and EM carries the gen1
isospin seed. The two are complementary, not competing.

## 27. E_6 TOWER + Q*Phi PIN, and the GEOMETRIC root of the up-sector: up = the SHARED crossings (2026-09-21)
Script: `src_paper16/e6_tower_em_pin.py`. User: pin E_6 tower with Q*Phi, keep T(2,n) geometry in mind
(trefoil vs cinquefoil; energy input pushes higher T; trefoils recombine in the soup sharing lobes,
cooling -> individuality, oscillating). This session's most useful step -- it explains WHY up resists a
clean tower.

CONCRETE RESULTS:
  (1) Q*Phi EM-subtraction WORKS. m_top = m_meas + Q*Phi. gen1 flips to up-heavier (topological
      ordering consistent) for Phi > 2.51 MeV; alpha*Lam_const=1.57 MeV, so an O(1) geometric factor
      ~1.6 gives Phi~2.5-3 MeV. Chosen Phi=3 MeV: m_top(u)=4.17, m_top(d)=3.67 (up heavier), gen2,3
      unchanged. So the linear EM piece (sec.26) does its job at gen1 and is negligible above.
  (2) The ABSOLUTE up tower does NOT close (verified): required (T_c-T_top)/(T_u-T_top)=0.46 vs the
      twist-gap 3/4=0.75; ln(m_top) not linear in |Tw|={1,2,5} (slopes 5.72 vs 1.64). Generation
      STEEPNESS works (E_6 phase 2.25, sec.24); absolute u,c do not fall out of an isolated-trefoil
      E_6 tower. Consistent with the n-T degeneracy (sec.22) and every prior attempt.

THE GEOMETRIC ROOT (user's shared-lobe picture, grounded) [I on the mapping, E on the geometry]:
P18:conj:isospin fixes the isospin geometry: **up-type = the CROSSING-VERTEX network (z=+r_0); down-type
= the MIDPOINT/distal-lobe network (z=0).** The crossings are exactly WHERE THE TWO STRANDS MEET -- the
SHARED/interaction points; the midpoints/distal lobes are the OUTER, individual parts of the tube. So:
  - DOWN-type quarks live on the INDIVIDUAL (midpoint) network -> admit a clean tower (E_8 lepton-bridge,
    anchored strange=muon; works to ~1.2, sec.24). 
  - UP-type quarks live on the SHARED (crossing) network -> their mass is inherently a SHARING/formation
    quantity, NOT an isolated-particle tower value -> no clean isolated tower should fit. THIS is the
    geometric reason the up sector resists (deeper than sec.17's 'no lepton anchor'): up sits on the
    shared crossings.
This is the user's "trefoils share a lobe/parts of the structure" made precise: the up-type charge IS
the shared-crossing network. Down-type is the individual-lobe network.

UNIFICATION (three threads meet):
  - The EM Q*Phi (sec.26) is PHYSICALLY the external field coupling to the z-LEVEL (up z=+r_0 vs down
    z=0) -- P18:conj:isospin's "field distinguishing the two z-levels." So Q*Phi and the crossing/midpoint
    geometry are the SAME thing: the field acts on z, up and down sit at different z -> signed shift.
  - Higher T(2,n) (cinquefoil T(2,5), |J(q5)|=phi; the Q_H=5 sector) = MORE crossings = MORE sharing.
    User's "energy input -> higher T" = (de)construction transiently accessing more-shared, higher-Q_H
    configurations. The "shared mess -> individuality, oscillating" on cooling = the JAM (sec.19) +
    the confined-quark oscillation (confinement-note sec.3, E_q/E_c=phi marginal).
  - So the up-sector mass is an O6 SHARED-FORMATION Delta E on the crossing network, not an isolated
    saddle/tower. Down-sector (individual midpoints) is why down is cleanly towered.

STATUS: Q*Phi pinned (Phi~2.5-3 MeV, O(1)*alpha*Lam). Absolute up tower does NOT close as an isolated
trefoil -- and that is now UNDERSTOOD, not just observed: up-type IS the shared-crossing network. NEXT
(O6, sharpened): compute the T(2,2)->T(2,3) formation Delta E on the CROSSING (shared) network vs the
MIDPOINT (individual) network -- the up/down split should emerge as shared-vs-individual formation
energy, with Q*Phi the z-level EM piece and the generation steepness the E_6 phase. Test whether the
crossing-network (shared) formation reproduces the up masses that the isolated tower cannot.

## 28. OPTION A (O6 analytic, crossing vs midpoint): the ISOSPIN SCALE closes (2026-09-21)
Script: `src_paper16/o6_crossing_vs_midpoint.py`. Built the golden trefoil (R0=3, r0=sqrt2/phi),
verified geometry, tested formation-energy proxies for M_up/M_down (up=crossing network, down=
midpoint/distal, sec.27).

GEOMETRY VERIFIED vs P18:prop:isospin_mass_ratio: unit-tangent T_z at the z=0 sites = +0.3206,
-0.5249 (P18: +0.321, -0.525); <T_z^2>_(z=0) = 0.1891 (P18: 0.189); crossing T_z=0 exact. My
trefoil = the paper's.

PROXIES for M_up/M_down (measured 6.38 w/ measured top, 7.95 w/ full-formation top):
  - P18 out-of-plane tangent fraction  (1-T_z^2|_C)/<T_z^2>_(z0) = 1/0.189 = **5.29**  <- BEST
  - local curvature kappa ratio (crossing/down)                 = 2.13
  - integrated bending INT kappa^2 ds (crossing-arc/z0-arc)     = 0.44  (WRONG direction)
  - inter-strand min-distance (sharing proxy)                   = 1.08  (crossings NOT much closer)
KEY: P18's tangent-fraction (5.29) is the best proxy, and the residual **measured/P18 = 6.38/5.29 =
1.207 ~ the SAME uniform ~1.2 running factor** the down-sector E_8 fit carries (sec.24; P18 itself
flagged the 21% as "a generation-universal correction, not yet identified"). So:
  **the ISOSPIN SCALE (M_up/M_down) is essentially CLOSED: crossing/midpoint tangent geometry (5.29)
  x uniform running (~1.2) = ~6.35 ~ measured 6.38, ZERO free parameters.** [1.207~1.2 flagged, not
  claimed as exact -- CLAUDE.md 3, one ratio.]

WHAT THE PROXY TEST TELLS US (refines sec.27):
  - The working mechanism is the TANGENT ORIENTATION (in-plane vs out-of-plane EXIT at the site),
    NOT literal bending energy (curvature ratio 2.13, integrated bending 0.44 wrong-way) and NOT
    literal spatial sharing (inter-strand min-distance ratio 1.08 ~ crossings are NOT much closer to
    other strands in this embedding). So "up = shared crossings" is TOPOLOGICAL/interactional (the
    crossing is where strands CROSS in the knot diagram), not a spatial-proximity statement. Honest
    correction to sec.27's "shared" wording: the crossing network's distinction is its IN-PLANE
    tangent exit, not spatial closeness.
  - Q*Phi (gen1 EM) and the E_6 generation steepness (sec.24) do NOT touch the triplet MEAN (dominated
    by c,t/s,b at GeV) -> the isospin scale is a clean, separate, now-closed piece.

COMMITTED TO PAPER (2026-09-21): added `P18:rem:isospin_scale_residual` after
`P18:rem:isospin_mass_ratio_scope` (main_paper18.tex) -- records (a) the independent tangent-data
verification + the residual ~1.2 matching the Paper XVII anchored-sector overshoot, (b) the operative
quantity = tangent ORIENTATION not bending/proximity. Refs resolve, remark envs balance, no
backslashed underscores. USER TO COMPILE (no pdflatex in-session). Scripts o6_crossing_vs_midpoint.py
etc. are Paper-18 physics currently in src_paper16 (move to src_paper18 pending user OK; convention:
quark-mass scripts live in src_paper16 on the trefoil geometry, XVII already cites e8_strange_anchor.py there).

VERDICT: Option A CLOSES the isospin SCALE (up/down triplet mean) analytically: P18 tangent geometry
+ uniform running, no free parameters, ~0.5% after running. What REMAINS open = the INTRA-triplet
(per-generation) absolute masses (the E_6 tower, sec.24, still not absolutely pinned) + the gen1 seed
(Q*Phi, sec.26). Option B (jam MD) would target the intra-triplet shared-formation, but the isospin
SCALE no longer needs it.

## 29. INTRA-TRIPLET: the DOWN hierarchy is a TWIST POWER LAW m ~ |Tw|^~4.2 (2026-09-21) [big]
Script: `src_paper18/intratriplet_twist_power.py`. Attacked the intra-triplet (the last open piece,
sec.28) analytically FIRST (CLAUDE.md 8) before the heavy jam MD. It largely cracked -- for DOWN.

DOWN {d,s,b}, twist |Tw|={1,2,5} (E_6, gen 1/2/3): **m ~ |Tw|^p FITS with p=4.219, R^2=0.99987.**
  - per-step p: d->s=4.322, s->b=4.148, d->b=4.223 -- CONSISTENT (the middle point s is not forced
    by a 2-pt fit, so the collinearity in log-log is a real 3-point power law).
  - one-amplitude prediction m=A|Tw|^p, A=4.79 MeV: d=4.79, s=89.3, b=4262 (meas 4.67, 93.4, 4180);
    residuals within 5%. Near-ZERO-parameter for the entire down generation hierarchy.
PHYSICAL ORIGIN [I, motivated]: the QUARTIC Faddeev term INT(dn x dn)^2 (= phi^6 J_4, the framework's
own energy, E_geom=K*J4 in the solver) scales as (twist)^4 for a twisted ribbon -> m ~ Tw^4. The
same quartic "4" as sin^4 theta / J_4 / the "three generations from sin^4" (claude-hopfion 3).
EXPONENT [O, flag not claim, CLAUDE.md 3]: fit p=4.22. Nearest: **phi^3=4.236 (0.4% from fit, ratios
6-8%)** best; 4 (pure quartic, physical) underfits (p 5% low, ratios ~20% off: m_s/m_d=16 vs 20). So
the ROBUST result is the POWER LAW (R^2=0.9999); the exact exponent (~4.2, phi^3-ish, quartic-motivated)
is secondary and not pinned. Tower form: n(g)=n_1 - (p/2 ln phi) ln|Tw| -- tower level LOG in twist.

UP {u,c,t}: NO clean power law (u->c p=9.2, c->t p=5.4 meas / 6.1 full-top) -- per-step p DECREASES,
the shared-crossing network (sec.27/28). up/down power ratio (u->c)/(d->s)=2.13 = the steepness (same
data). So the INDIVIDUAL (down/midpoint) triplet is analytic; the SHARED (up/crossing) triplet is not.

STATUS: intra-triplet DOWN hierarchy SOLVED analytically (twist power law, R^2=0.9999, quartic/phi^3,
one amplitude). This was the "heavy" open piece -- half of it fell to a light calc. What REMAINS for
the heavy jam MD (Option B) = the UP (shared-crossing) triplet ONLY, whose steepness decreases with
generation (sharing strongest at low gen, de-shares toward the top-transition). NEXT: either (i) probe
UP analytically as down-power-law x a shared-crossing correction (light), or (ii) the jam MD for UP.

## 30. OPTION 1 (Faddeev energy, crossing vs midpoint): a clean NEGATIVE that confirms sec.28 (2026-09-21)
Script: `src_paper16/faddeev_energy_crossing_vs_midpoint.py` (reuses the validated Construction-C
trefoil ansatz + E_geom; single field eval, ~3s at N=64, no gradient flow -- light, not "heavy").

FIRST, a correction to sec.29's physical origin [CLAUDE.md 1/9, caught while designing]: I claimed
"quartic Faddeev -> twist^4", but a geometric ribbon twist enters the Faddeev energy as (d_s n)^2 and
(F_{s psi})^2, BOTH ~ Tw^2, NOT Tw^4. So the data's mass~|Tw|^4.2 (sec.29) is NOT a direct
geometric-twist energy. Either mass ~ energy^2 (touches the phi^3=sqrt(phi^6) hint -- an energy->mass
square root) or the E_6 framing integer |Tw| is not a literal geometric twist rate. Mapping OPEN.
[** RETRACTED in sec.31 (user flag): the "~Tw^2" here is a GLOBAL frame twist; the GENERATION twist is
PER-STRAND/colour-fixed, a different object whose energy I did NOT compute. So the mass~energy^2 / phi^3
inference is premature -- it compared the generation-twist mass to the wrong twist's energy. See sec.31. **]

RESULT (N=64): the CROSSING (up, z=+r0) and MIDPOINT (down, z=0) networks carry ~EQUAL Faddeev energy
density per unit arc: quartic ratio E_up/E_down = 0.94 (N=32: 0.91, converging ~0.95), quadratic = 1.00.
And E(z=+r0) = E(z=-r0) EXACTLY -> the energy density IS z->-z symmetric even though the tangent
orientation (Gamma_z) is not (P18:conj:isospin's asymmetry is in ORIENTATION, not energy).

INTERPRETATION (confirms sec.28 AND P18:conj:isospin at the field level):
  - The up/down mass-scale asymmetry (M_up/M_down~6.4) is NOT a static energy difference between the
    crossing and midpoint networks -- they are energetically equivalent. It is the tangent-ORIENTATION
    coupling to an external field (P18:conj:isospin), exactly sec.28's conclusion.
  - Equal energy on both networks IS P18:conj:isospin's premise made explicit: "absent a field, both
    isospin values equally probable" -- because the two z-networks cost the same energy. The field
    breaks the degeneracy; the geometry (tangent orientation) sets which way and by how much (sec.28).
  - CLOSES the "up costs more static energy" route (the crossing network does NOT carry more energy).

IMPLICATION for the intra-triplet UP masses: a STATIC field-energy calc cannot produce the up/down
difference (networks are energetically equivalent) -- so neither this calc NOR the heavier jam MD
(also static-ish field energies) will crack the UP sector by energy. The up masses are an
ORIENTATIONAL / external-field-selection + formation-dynamics quantity, not a static energy. This
REDIRECTS Option B: don't measure static shared-formation energy for up; model the field-selection
(P18:conj:isospin field) + the Q*Phi z-level coupling (sec.26) dynamically.

STATUS of the up-sector program (consolidated):
  - isospin SCALE: closed (P18 tangent orientation x running, sec.28).
  - DOWN intra-triplet: closed (twist power law m~|Tw|^4.2, sec.29).
  - gen1 seed: Q*Phi linear EM (sec.26).
  - UP intra-triplet: NOT a static-energy problem (this sec) -> the genuinely open piece is the
    external-field selection + formation dynamics on the crossing network; no static shortcut.

## 31. TWIST DISAMBIGUATION (user, 2026-09-21): the mass~energy^2/phi^3 lead is PREMATURE
User flagged: make sure the generation twist isn't confused with the internal COLOUR ribbon twist
(the §10 error recurring). Checked -- one real confusion found, in the sec.30 inference (not the scripts).

THE DISTINCT FRAMINGS (§5, §10), so they are never blurred again:
  - COLOUR = writhe Wr=3 = the 3 crossings = global self-linking SL=3 (BPS-fixed (2,3)-cable framing,
    P18 l.466,536). Z_3/SU(3)_1. FIXED for ALL quarks -- not a variable.
  - ISOSPIN = the z-level (crossing z=+r0 vs midpoint z=0). The sec.30 networks.
  - GENERATION = per-strand INTERNAL ribbon twist |Tw|={1,2,5} (E_6), at FIXED knot type T(2,3)
    (colour SL=3 stays put; §5 "a family of framings, all T(2,3), differing by twist"). The sec.29
    power-law variable.
  - meridian winding theta=3t (the "3" in the ansatz Phi=chi+3t) = the cable/colour coordinate,
    DISTINCT from the Bishop-frame ribbon twist (§10); NOT a free variable.
  - Bishop holonomy -63 deg = geometric torsion, root-sibling (§15).

CHECK of the recent scripts:
  - sec.29 (intratriplet_twist_power.py): |Tw|={1,2,5} = per-strand GENERATION twist, colour fixed.
    CLEAN, no confusion.
  - sec.30 (faddeev_energy): localized the ISOSPIN z-networks on the fixed-colour trefoil. Correct
    for isospin. Does NOT use the generation twist.

THE CONFUSION (in sec.30's inference, now CORRECTED): "geometric twist -> Faddeev energy ~Tw^2" was
computed for a GLOBAL FRAME rotation. The GENERATION twist is PER-STRAND / internal (colour-fixed) --
a DIFFERENT object whose energy scaling I have NOT computed. So the chain "mass~|Tw|^4.2 vs energy~Tw^2
=> mass~energy^2 => phi^3=sqrt(phi^6)" compared the generation-twist MASS to the WRONG twist's ENERGY.
=> **RETRACT the mass~energy^2 / phi^3=sqrt(phi^6) inference as premature.** If the PER-STRAND
generation-twist energy scales ~Tw^4.2 directly, then mass~energy (linear), no square root, no phi^6
connection. The sec.29 power law m~|Tw|^4.2 STANDS (it is the mass vs per-strand framing at fixed
colour); only the energy-scaling INTERPRETATION (sec.29/30) is retracted pending a correct calc.

THE ACTUAL OPEN QUESTION (sharpened): compute the Faddeev energy of a PER-STRAND (internal, colour-SL=3-
preserving) ribbon twist on T(2,3) as a function of |Tw|, and see whether it scales ~Tw^2 (global-frame-
like) or ~Tw^4.2 (matching the mass directly) or otherwise. Only THEN is the mass-vs-energy relation
(and any phi-power) meaningful. This is the correct next calc -- and it needs the per-strand-twist
construction (§10's open item), NOT the global frame rotation.

## 32. PER-STRAND TWIST ENERGY (2026-09-21): geometric energy ~Tw^{1.5-2.8}, NOT Tw^4.2 -> mass is ALGEBRAIC
Script: `src_paper16/perstrand_twist_energy.py`. Built the honest version of sec.30: a twisted director
tube n=(sin f cos Phi, sin f sin Phi, cos f), Phi=m*psi+q*(2pi z/L), and measured how the framework
energy E_geom=K*J4 scales with the winding (via the SAME finite-difference E_geom).

RESULT (log-log slopes of E=K*J4):
  (A) longitudinal winding q (framing/twist), m=1: E ~ q^1.49
  (B) meridian winding m, q=1:                     E ~ m^2.60
  (C) both m=q=T:                                  E ~ T^2.84
None is ~Tw^4. The physical saddle energy E_phys ~ sqrt(K*J4) (Derrick/Bogomolny E2=E4 at optimal
scale) is even lower, ~Tw^{0.75-1.4}. And no clean mass~energy^k: 4.2/2.84=1.48, 4.2/1.49=2.8 (neither
a clean power).

CONCLUSION (closes the thread, ties to sec.1):
  - The geometric Faddeev TWIST energy scales ~Tw^{1.5-2.8}, NOT the data's mass~|Tw|^4.2 (sec.29).
  - DEEPER: the quark has NO stable saddle and its mass is ALGEBRAIC (sec.1, P17:rem:quark_mass_algebraic).
    So relating the generation twist to the quark MASS via ANY Faddeev configuration energy was
    misframed -- there is no saddle energy to equate the mass to. (The LEPTON, Q_H=2, has a saddle and
    its mass IS a configuration energy; the QUARK does not.)
  - => the twist->quark-mass relation m~|Tw|^4.2 is ALGEBRAIC (the E_6 WZW T-matrix / tower level
    n ~ -ln|Tw|), NOT a geometric twist energy. This CONFIRMS sec.1 and DEFINITIVELY closes the
    mass~energy^2 / phi^3=sqrt(phi^6) thread (retracted sec.31): the mass is not the Faddeev energy at
    all, so there is no energy->mass power to carry a phi^6.
  - The geometric twist energy (~Tw^{1.5-2.8}) is a REAL quantity (relevant to the confinement/formation
    dynamics, and to the LEPTON sector where masses ARE saddle energies), just not the quark generation
    mass. The exponent 4.2 lives in the algebraic E_6 tower, to be understood there (why n ~ -ln|Tw|),
    not in geometry.
STATUS: the per-strand-twist ENERGY question is answered (scales ~Tw^{1.5-2.8}, real but not the quark
mass). The open piece is now purely ALGEBRAIC: derive n ~ -ln|Tw| (m~|Tw|^4.2) from the E_6 WZW tower
-- a WZW/CFT calc, not a field-energy one. No phi^6 / square-root motivation survives.

## 33. REP-THEORY VERDICT (2026-09-21): m~|Tw|^4.2 is EMERGENT, not a WZW identity
Script: `src_paper18/e6_tower_reptheory.py` (sympy exact rationals + transcendentals). Tested whether
the down power law (sec.29) is an algebraic identity from the E_6/E_8 WZW data or a numerical coincidence.

TELL 1 -- per-step exponent NOT exactly constant: p(d->s)=4.167, p(s->b)=4.223 -> 1.3% spread. A true
identity would give EXACTLY constant p (spread ~0). 1.3% = a good-but-inexact fit.
TELL 2 (decisive) -- the clean power law REQUIRES the measured leptons. The E_8 route is m_down(g) =
m_lepton(g) x exp(n_q pi/9); the lepton ratios (206.8, 16.8) x E_8 factors (0.087, 2.85) = (17.96, 47.9)
"conspire" to near-equal twist-powers. But dropping the measured leptons and using the rep-theory
formula pieces alone (exp(-2pi*10*T_g) x exp(n_q pi/9)) gives p(d->s)=4.03 vs p(s->b)=4.57 -- 13% spread,
NOT a clean power. So the |Tw|^4.2 regularity lives in [MEASURED leptons] x [E_8 exponents], not in the
pure WZW data.
TELL 3 -- p~4.2 vs phi^3=4.236 (0.9% off) is WITHIN the fit's own 1.3% scatter -> "p=phi^3" is not
established; it is inside the coincidence's noise. (2phi+1=phi^3 same; h(E_6)/e=4.41 worse.)

VERDICT: **m~|Tw|^4.2 is an EMERGENT numerical near-coincidence (good to ~1.3% over 3 points), NOT a
fundamental rep-theory identity.** The FUNDAMENTAL structure is the E_8 route (lepton mass x E_8-exponent
correction, established, works to the ~1.2 uniform running); the twist power law is a re-parametrization
that happens to fit. So: stop chasing a rep-theory / phi^3 / phi^6 "meaning" for the exponent 4.2 -- there
isn't one at the identity level. (R^2=0.9999 over 3 points with 2 fit params + near-collinearity is easy;
it does not certify a fundamental law -- the per-step spread and the leptons-required test do.)

CONSEQUENCE for the program: the DOWN hierarchy is DERIVED by the E_8 route (m_lep x exp(n_q pi/9),
sec.24), full stop; the |Tw|^4.2 (sec.29) is a nice empirical restatement, not a new law. This CLOSES
the twist->mass thread cleanly: the quark mass is algebraic (E_8 exponents + lepton reference, sec.1/32),
the geometric-twist and phi^6 interpretations are both dead, and the down sector needs no further
"exponent derivation" -- it is the E_8 route.

## 34. DOUBLET SUM RULE (user intuition, 2026-09-21): (m_up+m_down)=C*m_lepton, C~13 -- ties UP to leptons
User: down-from-lepton makes sense; up? "the full mass of a lepton must go somewhere." Tested
(src_paper18/doublet_sum_rule.py, up_sector_e6_vs_e8_lepton.py) -- the intuition WORKS, as a per-doublet
mass budget (NOT a per-quark relation, and NOT product/seesaw -- those fail).

RESULT: (m_up + m_down)/m_charged_lepton per generation:
  g1 (u,d / e):   13.366     g2 (c,s / mu):  12.904     g3 (t,b / tau):  99.6 [BARE TOP -- excluded]
  CONFINED doublets g1,g2: C = 13.37, 12.90 -> **3.5% spread**. C_bar ~ 13.1.
CROSS-PREDICTION (non-circular -- C from one doublet predicts the OTHER's up quark):
  C(g2)=12.90 -> m_u = C*m_e - m_d = 1.92 (meas 2.16, 11%);  C(g1)=13.37 -> m_c = C*m_mu - m_s = 1319
  (meas 1270, 3.9%). So the doublet sum rule predicts the confined up quarks to ~4-11% from leptons+downs.

WHY THIS MATTERS (resolves the sec.17 up-anchor puzzle):
  sec.17 said the UP sector has NO massive lepton anchor (its isospin partner is the ~meV neutrino). The
  sum rule shows UP IS anchored to the CHARGED lepton -- not directly, but through the DOUBLET: the
  doublet's total quark mass is ~13 m_lepton, and since DOWN is lepton-anchored (strange=muon, E_8), UP
  takes "the REST of the budget": **m_up = C*m_lepton - m_down**. This is exactly the user's "full lepton
  mass goes somewhere" -- the charged-lepton mass is the doublet's budget, split between up and down.
  It also explains why the top breaks it: the top is the bare transition (v/sqrt2, sec.23), not a confined
  doublet member, so it is not bound by the confined-doublet budget.

CAVEATS (CLAUDE.md 3/6): only 2 confined doublets -> C~13 rests on 2 points; not derived (sits between
4pi=12.566 and 13 = the E_8 strange-anchor exponent = 2 Q_group^nu+1; flagged, not claimed). The
PRODUCT (seesaw) m_u m_d / m_lep^2 does NOT hold (38.6, 10.6 -- scatters); only the SUM does. So this is a
strong LEAD, not a closed result.

STATUS: the open UP intra-triplet (sec.30/33) now has its first real handle -- a doublet mass-budget sum
rule tying up to the SAME charged leptons that anchor down, predicting m_c to <4% and m_u to ~11%. NEXT:
(i) derive C (is it 13 = anchor exponent, or 4pi?); (ii) a 3rd confined check is impossible (only 3
generations, top bare) -- so C stays 2-point unless it is derived; (iii) reconcile with the E_6 steepness
(sec.24) and the field-selection (sec.30) -- the sum rule may BE the "formation dynamics" answer for up.

## 35+36. DERIVING C: reliable data leans to 4pi + quark-spin, not pure 4pi or pure 13 (2026-09-22)
Scripts: `src_paper18/derive_C_budget.py`, `C_4pi_vs_group.py`. User: Q_b+Q_l=13 is the nicest GROUP
form, but 4pi is natural/phi-related/framework-native (electron mass e^{4pi}); the g1/g2 error pattern
differs between candidates. KEY (uncertainty-weighted): the two doublets have VERY different reliability:
  g1 (u,d/e):  C = 13.37 +- 1.11  (light quarks ~10% -> NEARLY UNCONSTRAINING)
  g2 (c,s/mu): C = 12.90 +- 0.19  (c,s ~1.5% -> the RELIABLE ANCHOR)
So the error asymmetry the user saw is mostly which point to TRUST. Against the reliable g2:
  pure 4pi=12.566        : -1.78 sigma  (mildly DISFAVOURED, a bit low)
  4pi + 1/3 = 12.900     : -0.02 sigma  (BULLSEYE; 1/3 = h_SU(3)_1 = the quark conformal weight)
  13 = Q_b+Q_l           : +0.51 sigma  (fine, within 1 sigma)
  4pi + 1/phi = 13.184   : +1.48 sigma
RESULT: the reliable data leans to **C ~ 12.9 = 4pi + (small quark correction)**, mildly disfavouring
PURE 4pi and consistent with 13. Nicest framework-native form: **C = 4pi + h_SU(3)_1 = 4pi + 1/3** =
(the solid-angle factor that sizes the lepton mass, e^{4pi} in m_e) + (the quark's own spin 1/3) --
matches the user's 'quark conforms to the lepton (4pi) + its own quark nature' picture, better than the
pure-group 3+10 or pure 4pi.
CAVEATS (CLAUDE.md 3, HARD): (i) 4pi+1/3 found POST-HOC (residual C2-4pi=0.3375 ~ 1/3=0.333); (ii) the
residual is ALSO ~ pi/9=0.349 (the quark E_8 phase unit) -- 1.3% vs 3.6%, so even the correction is not
unique; (iii) ONE reliable point (g2); (iv) 13 within 1 sigma. So this is a LEAD refining the user's 4pi,
NOT a derivation. Data cannot select 4pi+1/3 vs 13 vs 4pi+pi/9.
STATUS: C ~ 12.9, base = 4pi (vindicates the user's geometric instinct; echoes m_e ~ e^{4pi}), + a small
quark correction (~1/3 = h_SU(3), or pi/9). The ARBITER remains the O6 Q_H=2+1->3 formation mass-balance,
which would derive C (and pick the correction) rather than fit it.

## 38. O6 / the 4pi mechanism CHECKED (2026-09-22): framework's 4pi's don't give C's factor-of-4pi
User: maybe C's 4pi is like the fine-structure "1/(4pi) per unit WZW level" vacuum-polarisation. Checked
against the papers -- good instinct (the structure is REAL), but the magnitude/role does not match.
GROUNDED (the framework's 4pi's):
  - CS SPOKE amplitude C_CS = 4pi*Delta Q (P17:eq:CS_el, from HI[Q]=Q -> INT A^dA=16pi^2 Q, WZW norm
    1/(4pi)). Lepton (Q_H=2, Delta Q=1): 4pi. Quark constituent (Delta Q=3-1=2): 8pi (Paper XI l.1051).
    But these enter masses EXPONENTIATED: exp(4pi)=2.9e5, exp(8pi) astronomical. NOT C~12.9.
  - VACUUM-POLARISATION 1/(4pi) PER WZW LEVEL (Paper IV l.867-917, the user's reference): suppresses the
    QED vac.pol. by (1 - 1/(4pi)) ~ 0.92/level (~8% correction, the alpha^-1 4th term -1/(36 pi phi^6)).
    A SMALL correction, not a factor of 4pi.
  - MASS TOWER = phi^2 ~ 2.618 PER LEVEL (Paper VI/XI l.179), not 4pi.
CHECK: C ~ 4pi (~12.9) would need 4pi as a LINEAR factor-of-12.9 mass ratio. NONE of the three framework
4pi/per-level structures gives that: exp(4pi) too big, 1/(4pi) too small (~8%), phi^2 wrong value
(phi^{2n}=12.9 -> n=2.66, non-integer). So the vacuum-polarisation route does NOT produce C's 4pi.
CONSEQUENCE (honest): the appealing C=4pi+1/3 has NO derived mechanism among the framework's 4pi's. A
genuine 4pi-as-mass-factor would need a SOLID-ANGLE / phase-space 4pi (INT_{S^2} dOmega = 4pi, the director
spreading over the full sphere) -- geometrically plausible but NOT an established framework 4pi. Combined
with sec.36 (data can't separate 4pi+1/3 from 13), this tips slightly AGAINST the clean 4pi reading:
numerically attractive, mechanistically unsupported so far.
STATUS: the 1/3 = SU(3) charge (sec.37) is the SOLID part. The ~13 magnitude of C is genuinely UN-derived
-- neither 4pi (no factor mechanism) nor 13=Q_b+Q_l (a sum of group orders, unmotivated as a mass budget)
has a real derivation. O6 remains the arbiter, and it must produce the ~12.9 factor from formation
GEOMETRY (a solid-angle 4pi is the best remaining candidate, distinct from the CS/vac-pol 4pi's). Do NOT
re-chase the vacuum-polarisation 1/(4pi) route for C (ruled out on magnitude).

## 39. O6 SOLID-ANGLE computed (2026-09-22): NO 4pi jump -- the solid-angle route is FALSIFIED
Script: `src_paper18/o6_solid_angle_Q2_Q3.py`. Built standard hopfions (stereographic S^3->S^2, w=Z1^A/
Z2^B), Hopf charge verified by FFT, computed the director geometry across Q_H=2->3.
CHARGE CONVERGENCE (N=96): (1,1)->1.00, (1,2)->1.95, (1,3)->2.82 (converging to 1,2,3 -- genuinely the
Q_H=1,2,3 configs).
GEOMETRY ACROSS Q_H=2->3 (ratios, the physical content):
  J2 (Dirichlet)  Q3/Q2 = 1.21
  J4 (Faddeev)    Q3/Q2 = 1.80
  solid angle INT|F|  Q3/Q2 = 1.33
ALL modest (~order the charge ratio 3/2=1.5). NONE is 4pi (12.57) or C~12.9. So the director's S^2
solid angle (and J2, J4) does NOT jump by 4pi across the lepton->quark step -- it grows by ~1.3.
=> the "solid-angle 4pi" hypothesis (sec.37/38, the LAST candidate for C's 4pi) is FALSIFIED.

CONSOLIDATED VERDICT on C (sec.35-39): the doublet budget C~12.9 has NO derived mechanism:
  - CS spoke 4pi*DeltaQ: exponentiated (exp(4pi)~3e5) -- too big (sec.38).
  - vacuum-pol 1/(4pi) per WZW level: ~8% correction -- too small (sec.38).
  - director solid angle across Q_H=2->3: ~1.33 -- too small, not 4pi (this sec).
  - group 13=Q_b+Q_l: a clean integer sum but no MASS-budget mechanism (sec.35).
So C~12.9 (~4pi+1/3 numerically) is an EMPIRICAL constant; the ONLY principled part is the +1/3 =
SU(3) charge (sec.37). The 4pi is numerology-adjacent (fits, no mechanism). Down-sector derived (E_8),
isospin scale (sec.28), gen1 seed (Q*Phi, sec.26) all stand; the doublet SUM RULE (m_u+m_d)=C*m_lep is
a strong EMPIRICAL regularity (cross-predicts u,c to ~4-11%) with C un-derived.
STATUS: the O6/4pi chase is DONE (all mechanistic routes for C's magnitude exhausted -- negative). Do
NOT re-chase C's 4pi via CS, vacuum-pol, or solid-angle (all tried, all fail). If C is ever derived it
needs a genuinely new mass-budget principle for the Q_H=2+1->3 formation, not the geometric/CS 4pi's.

## 40. CONCEPTUAL (user, 2026-09-22): the "C ~ 4pi" reading CONFLATES Hopf charge with a mass ratio
User asked (a) are we conflating Hopf charge with the mass ratio, (b) why is m_e ~ exp(4pi-3/400)?
Read Paper IV (l.1166-1255). Both answered; the answer reframes sec.36-39.

(b) WHY exp(4pi-3/400): the electron is the j=1/2 BOUNDARY STATE of the SU(2)_3 WZW CFT. Its mass is an
open-string PARTITION-FUNCTION AMPLITUDE, Z ~ exp(-action) -> masses are EXPONENTIALS of an action. The
Chern-Simons action of the topological transition is Delta S_CS = 4pi (= 4pi*Delta Q, Delta Q=1; from
INT A^dA = 16pi^2 Q divided by the WZW norm 1/(4pi)). So the leading factor is e^{4pi} (universal CS
spoke). -3/400 = -T_{1/2}/Q = the WZW modular T-matrix (topological-spin) phase; phi^20=phi^{2Q}
(Verlinde/tower, Q=10); -alpha/pi = QED scheme. So 4pi is a CS ACTION IN THE EXPONENT; exponentiated
because mass = exp(-S) amplitude. NOT a linear 4pi.

(a) CONFLATION -- YES, and this is the key: HOPF CHARGE ENTERS MASS EXPONENTIALLY (e^{4pi*DeltaQ} + tower
phi^{2Q}). The natural Q_H-sector mass RATIO is ~ e^{4pi*DeltaQ} phi^{2 DeltaQ} = HUGE (e^{4pi}~2.9e5),
NOT ~13. C=(m_u+m_d)/m_lep~12.9 is small because it compares CURRENT quark masses (tower-suppressed far
below constituent) to the lepton, and the big e^{4pi} is COMMON to both absolute scales -> CANCELS in the
ratio. C is a RESIDUAL set by the phi-tower powers + E_8/E_6 exponents, NOT by e^{4pi}. So "C ~ 4pi"
compares a MASS RATIO (12.9) to the BARE ACTION VALUE (4pi) -- but mass ~ e^{4pi}, not 4pi. CATEGORY
ERROR. This is EXACTLY why sec.38/39 found no mechanism: there cannot be one (a mass ratio can't equal a
bare action when masses ~ exp(action)).
CONSEQUENCE (reframes sec.36-39): the "C = 4pi + 1/3" reading is a CONFLATION -- the 4pi is a numerical
COINCIDENCE (C ~ bare CS action value; no mechanism, and none possible). The +1/3=SU(3) charge CAN enter
linearly (charge enters linearly via the EM Q*Phi term, sec.26), so that small piece is defensible, but
it is a tiny correction to a coincidental 4pi. => DROP the 4pi interpretation of C. The doublet SUM RULE
(m_u+m_d)=C*m_lep STANDS as an empirical regularity (cross-predicts u,c ~4-11%), but C~12.9 is a
tower/exponent RESIDUAL, its magnitude un-derived, and its "~4pi" has no physical content.
NET for the up-sector: down derived (E_8); isospin scale (sec.28); gen1 seed (Q*Phi, sec.26); the doublet
sum rule is EMPIRICAL with C a residual (not 4pi, not Hopf). The one principled fragment is +1/3=SU(3)
charge (linear, EM). No 4pi/phi^6/solid-angle meaning survives (secs 31,32,33,38,39,40).

## 41. OPEN AVENUE (user, 2026-09-22): the CS DYNAMICS / edge modes / Hall response -- UN-EXPLORED
User: we've used the CS ACTION VALUE (4pi=CS action of the topological transition) but not its DYNAMICS
-- varying S_CS w.r.t. the external EM potential gives the induced/Hall CURRENT; the boundary supports
CHIRAL EDGE MODES that carry (fractional) electric charge (quantum-Hall edge transport). "Doing all this
without the equation of motion?" CORRECT -- and it may be where fractional charge really lives.

CONFIRMED (grep all papers): NO edge-mode / Hall-response / induced-current / chiral-edge / filling-
fraction / CS-equation-of-motion content anywhere. The framework does the CHARGE sector STATICALLY:
  - fractional charge <- Z_2 WRITHE asymmetry (P18:conj:isospin), kinematic assignment.
  - alpha^-1 = k|2I|/V*^2 (Paper XIX l.1135), WZW S-MATRIX modular data (V*=phi), static.
  - 4pi <- topological action value INT A^dA=16pi^2 Q, in the mass exponent.
We (and the framework) NEVER varied the action: no delta S/delta A, no induced current j=sigma_xy*(dual F),
no Hall response, no edge modes. CS theory without its dynamics.

WHY IT MATTERS: fractional charge is intrinsically an EDGE/HALL phenomenon in CS theory (bulk gapped/
topological; charge-carrying excitations on the boundary; FQH filling nu -> quasiparticle charge e*nu).
FRAMEWORK-NATIVE HOOK: the charge sector's CFT is SU(2)_3 -> level **k=3**; an abelian CS/FQH state at
level 3 gives nu=1/3 -> quasiparticle charge **e/3 = down-quark charge**, 2e/3 = up-quark charge. The
level-3 -> THIRD-integer charges is exactly FQH edge physics -- DYNAMICAL, not the static writhe
assignment. So the "+1/3" we kept finding (secs 34-40) may be the fractional HALL/EDGE charge nu=1/k
(k=3), NOT "the SU(3) conformal weight". The whole charge sector may be edge dynamics we skipped.

WHAT THIS CHANGES vs closing:
  - STILL SOUND to close: the "C~4pi" MASS-BUDGET reading (category error, sec.40). Independent of this.
  - MUST NOT close: the CS-DYNAMICAL sector (EOM, induced current, Hall, edge modes) is UN-EXPLORED and
    is the natural home for fractional charge (+ via edge-current energy + the linear Q*Phi coupling
    sec.26, a potential feed into the MASS). Reframes the open question from "derive C's magnitude" to
    "does the CS dynamics give the fractional charges (nu=1/k, k=3 -> 1/3,2/3) and their currents/
    energies, done only statically so far?"
CAVEAT (CLAUDE.md 3): SU(2)_3 is NON-abelian; the clean nu=1/k -> e/k is ABELIAN (U(1)_3). The EM-U(1)
Hall response of the condensate must be worked out = the un-built CS-EOM computation. So k=3->1/3 charge
is a SUGGESTIVE, framework-native LEAD, not a result.
NEXT (if pursued): the condensate's CS equation of motion / induced EM current / edge-mode structure --
does fractional charge emerge dynamically (nu=1/k=1/3), and does the chiral edge current carry energy
that contributes to the quark mass? This is a genuinely NEW sector, orthogonal to the (static) mass-tower
and charge-assignment work. Do NOT close the up-sector program without flagging this as the live avenue.

## 42. RESIDUAL DIAGNOSTIC (step 1, 2026-09-22): static=exact leading, residuals=the dynamical layer
User strategy CONFIRMED. Script: `src_paper18/residual_diagnostic.py`. Source: Paper VII summary table
(main_paper7.tex l.2590-2652) + this session's quark residuals. Classified every quantitative residual
by whether its sector is a STATIC SADDLE (weakly-coupled, stable) or DYNAMICAL/LOOP (strong/transient):
  EM / charged-lepton (STATIC SADDLE): alpha^-1 5e-9, lambda 2e-4, m_e 1.3e-4, m_H 1e-3, sin^2thW 2.3e-3,
    m_mu/m_tau 5e-3  -> ALL <= 0.5%, geomean 0.011%.
  cosmology (MIXED): rho_CMB/Lcond 1.2e-3, T_CMB chain 2.2e-3, rho_infty 3.4e-2, P_s 8e-2, g_W 9e-2
    -> 0.1-9%, geomean 1.5%.
  quark / neutrino (DYNAMICAL/LOOP): neutrino 28%, down 'running' 20%, isospin ratio 21%, up u,c ~factor
    -> ALL >= 20%, geomean 33%.
RESULT: a **~40x CLEAN GAP** between static-saddle (<=0.5%) and dynamical (>=20%) sectors; cosmology
(mixed static+dynamical inputs) sits between. The dynamical residual scale ~20-30% ~ **alpha_s** (strong
coupling: alpha_s(2GeV)~0.30, alpha_s(m_c)~0.35), vs the EM/lepton scale ~alpha_em/pi~2e-3. So RESIDUAL
SIZE TRACKS THE SECTOR'S DYNAMICAL COUPLING (alpha_em tiny -> sub-%; alpha_s~0.3 -> ~20-30%).

INTERPRETATION (confirms user, sec.40-41): the framework's static/topological results are the EXACT
LEADING ORDER; the residuals ARE the dynamical layer, LARGEST where the object has no static saddle
(the quark: metastable/transient, sec.1). Charged leptons (stable saddles) -> static is exact to sub-%;
quarks (QCD/confinement, no saddle) + neutrinos (seesaw loop) -> ~alpha_s/loop residuals. This is not a
collection of unrelated misses -- it is ONE dynamical correction whose size = the sector's coupling.
CONSEQUENCE for the program: don't chase individual quark residuals (C, running, isospin, up) as static
constants -- they are the STRONG/FORMATION dynamics (alpha_s scale), forced in the quark sector. The
right target is the DYNAMICAL LAYER:
  - quark MASS residuals (~20-30% ~ alpha_s): step (3) formation/confinement dynamics (sec.18/19 QGP-
    jamming, the transient formation energy) -- the alpha_s-scale correction to the static tower.
  - fractional CHARGE: step (2) CS-edge/Hall dynamics (sec.41) -- a distinct dynamical handle.
CAVEATS: sector classification has some judgment; session quark residuals have uncertainties. But the
40x split is robust to reclassification, and the alpha_s identification is order-of-magnitude (suggestive).
STATUS: step (1) DONE and POSITIVE. The framework is a static/topological theory exact at leading order,
with a systematic dynamical layer (~coupling per sector) that is UN-DEVELOPED and forced in the quark
sector. NEXT (highest value): develop the strong/formation dynamics (step 3) for the quark mass residuals
-- the alpha_s-scale correction -- and/or the CS-edge dynamics (step 2) for fractional charge.

## 43. STEP 3 SCOPE (2026-09-22): the residuals are NON-PERTURBATIVE (formation/confinement), -> Paper 20
Opening step 3 (develop the dynamical layer for the quark mass residuals). Grounding check changed the
scope IMPORTANTLY:
  - The PERTURBATIVE alpha_s layer is ALREADY BUILT: Paper V/XVI derive a framework-internal
    alpha_s(LUV) = C_F sin(pi/5)/(2 phi^{43/7}) ~ 0.0204 (C_F=4/3, phi, pi/5, d_eff=43/14; verified 0.3%
    vs PDG), and Paper XVI runs the quark masses via QCD RG (P16:tab:rg_running).
  - CRUCIAL (P16 l.1928): RG running changes entries by "at most a few percent, leaving the pattern
    (one near-exact match, FIVE large undershoots, generation-growing deviation) UNCHANGED." So the
    perturbative alpha_s running does NOT close the ~20-30% residuals -- they PERSIST.
=> The residuals are NON-PERTURBATIVE: the metastable-quark FORMATION/CONFINEMENT dynamics (no static
saddle, sec.1; the transient formation, sec.18/19). Consistent with sec.42's ~alpha_s~0.3 scale being
exactly where perturbation theory becomes O(1) and non-perturbative effects dominate. The undershoots
(static tower UNDER-predicts) = the static/perturbative skeleton MISSING the positive non-perturbative
formation/confinement energy.

STEP 3 = Paper 20: the NON-PERTURBATIVE formation/confinement dynamics of Q_H=3, as the source of the
persistent ~20-30% (generation-growing) quark-mass residuals. Pieces (all framework-native, partly begun):
  (A) FORMATION ENERGY: mass = static tower + transient formation energy (barrier + string sigma*R +
      kinetic) in the dynamical medium. sec.18 (top=incomplete formation, v/sqrt2), sec.19 (sigma*R
      proton=924.8=98.6% m_p). Needs the TOPOLOGY-PROTECTED flow (sec.19 checklist; the unconstrained
      solver dilutes -- the blocked piece).
  (B) CONFINEMENT KINEMATICS: c_s=c/phi Cherenkov, tube waveguide (confinement note secs 1-4) -- the
      framework's own strong-dynamics kinematics; the alpha_s-scale corrections.
  (C) ISOSPIN/CROSSING dynamics (sec.30): up/down 21% from crossing(shared)/midpoint(individual)
      field-selection -- the orientational, not-static piece.
FIRST CALC (tractable, analytic -- NOT the blocked solver): take Paper XVI's SURVIVING residuals (post-RG
undershoots, generation-growing) and test the FORMATION-ENERGY hypothesis analytically -- does a
sigma*R-type + barrier formation correction (sec.19 handle-1 style, using the derived sigma and R=R0+r0)
reproduce the undershoot pattern (positive, generation-growing)? This connects sec.42's diagnostic to
sec.18/19's formation energy WITHOUT the blocked gradient-flow solver. If it tracks, the residuals are
the formation energy; then (heavier) the topology-protected flow quantifies it.
CAVEAT: the full formation energy needs the topology-protected solver (open/blocked, sec.19). The first
analytic calc is a hypothesis test, not the full computation.
STATUS: step 3 scoped. It is a coherent NEW SECTOR (non-perturbative formation/confinement dynamics) =
Paper 20. The perturbative alpha_s (Paper V/XVI) and the static tower (E_8) are the DONE leading layers;
step 3 is the non-perturbative correction they leave.

## 44. STEP 3 FIRST CALC (2026-09-22): formation-energy hypothesis FALSE on scale; residual is ILL-POSED
Script: `src_paper18/step3_formation_energy_test.py`. Tested sec.43's hypothesis (residuals = sigma*R +
barrier formation correction) against Paper XVI's post-RG colour-shift residuals (P16:tab:rg_running).

RESULT 1 -- FORMATION-ENERGY (sigma*R) HYPOTHESIS IS FALSE (scale): the 'missing' energy Dm=m_PDG-m_static
spans 0.08 - 169000 MeV (SIX orders), but sigma*R ~ 278 MeV is a ~CONSTANT formation/confinement scale.
A constant sigma*R cannot match a 6-order spread. => the current-quark-mass residual is NOT the sigma*R
formation energy.
RESULT 2 -- sec.43 CONFLATED two non-perturbative objects (my error, caught by the calc):
  (a) FORMATION/CONFINEMENT energy sigma*R~278 MeV = the CONSTITUENT/HADRON scale. sec.19: proton =
      sigma*(3R0+R)=924.8=98.6% m_p. This is HADRON physics, largely DONE. NOT a current-mass correction.
  (b) CURRENT-quark-mass residuals (sec.42) = a DISTINCT non-perturbative correction to the current-mass
      tower, multiplicative & route-dependent (colour-shift: u1.7 d9.4 s1.0 c4.8 b2.2 t49), NOT sigma*R.
RESULT 3 (the deeper find) -- the CURRENT-MASS RESIDUAL IS NOT EVEN WELL-POSED: the framework has TWO
static quark-mass routes that DISAGREE:
  - E_8 Coxeter route (Paper V, sec.24): m_q/m_ell=exp(n_q pi/9); m_s(cond)=m_mu=105.66; down ~1.2 OVER,
    up under.
  - COLOUR-SHIFT route (Paper XVI l.581): m_q=m_ell*colour[2T/E_6/SU(3)]; m_s(LUV)=22.56->93.32 via RG;
    strange EXACT, ALL others UNDER.
Same strange-singled-out, but different numbers/patterns (P16:rem:two_strange_results notes they are two
independent calcs). So 'the residual' depends on which route -> ill-posed until the route is fixed.

CONSEQUENCE (reorders step 3): BEFORE computing any dynamical/non-perturbative correction, the framework
must RESOLVE its static quark-mass route (E_8 Coxeter vs colour-shift -- two mechanisms, disagreeing
condensate-scale masses). Only then is the residual a single well-defined quantity to correct. This is a
more foundational open problem than sec.43 assumed. sec.19's sigma*R (hadron formation) STANDS and is
separate/done; sec.42's diagnostic (residuals ~alpha_s, non-perturbative) STANDS; but step 3's target
(the current-mass residual) needs the static route settled first.
STEP 3 REORDERED: (i) reconcile/select E_8-Coxeter vs colour-shift static route [foundational, likely a
close read of Papers V/XVI + which is more derived]; (ii) define the residual of the chosen route;
(iii) THEN the non-perturbative correction (which is NOT sigma*R formation -- that is hadrons). Paper 20
scope now includes the route reconciliation as step 0.

## 45. CORRECTION to sec.44 (user, 2026-09-22): the route IS resolved; PAPERS ARE BEHIND THE NOTES
sec.44 claimed the static quark-mass route is "ill-posed -- two routes disagree, must reconcile." WRONG
-- I grounded it in the BEHIND-PAPERS (Paper V pure-E_8; Paper XVI colour-shift) instead of THIS
conversation's resolved position. User corrected:
  - Paper V's E_8-for-the-ENTIRE-scale is the EARLY, SUPERSEDED version.
  - THIS conversation resolved (sec.16, sec.24): **DOWN = E_8 lepton-bridge** (strange=muon anchor;
    pred/meas {1.26,1.13,1.21}, the uniform ~1.2); **UP = E_6-native** (steeper, phase pi/4, steepness
    2.25). E_6 = quark's native group; E_8 = lepton BRIDGE (P17:rem:quark_mass_algebraic is sound as the
    bridge). NOT ambiguous.
  - The Paper XVI colour-shift table (sec.44's data) is ANOTHER behind-version, NOT our route.
RETRACT sec.44's "ill-posed / reconcile the route." The route is settled in the notes: E_8-bridge(down)/
E_6-native(up). sec.44's scale result (sigma*R = HADRON scale, sec.19, NOT the current-mass residual)
STANDS; only the "two-route ambiguity" framing is retracted.

THE REAL FINDING (user, key): **THE PAPERS ARE WELL BEHIND THE NOTES.** Everything this conversation
established is NOTE-ONLY: E_6-native/E_8-bridge split (sec.11-16,24), Q_H relative convention (committed
earlier), isospin scale (sec.28), Q*Phi gen1 (sec.26), doublet sum rule (sec.34), residual diagnostic
(sec.42), the static-leading/dynamical-residual framing (sec.42-43), the retractions of 4pi/phi^6/solid-
angle (sec.31-40). The papers still carry Paper V pure-E_8, Paper XVI colour-shift, P17 E_8-only. So
re-deriving "from the papers" (as sec.44 did) trips on the STALE state.

CONSEQUENCE for step 3 / the residual (well-posed now, via the NOTES' route):
  - DOWN (E_8-bridge): ~1.2 OVER (uniform) -- consistent with condensate->observation running; small,
    arguably closed as running.
  - UP (E_6-native): UNDER (the steeper tower); absolute u,c don't close (sec.24). THIS is the genuine
    open non-perturbative piece -- NOT a universal sigma*R (sec.44), NOT the colour-shift table.
So step 3's target = the UP-sector E_6-native residual, in the notes' route.

RECOMMENDATION (reprioritise): the highest-value move is now CONSOLIDATION -- the resolved physics is
note-only and the papers are behind, so exploratory calcs keep tripping on stale papers. Write up the
resolved synthesis (E_8-down/E_6-up, isospin scale, Q*Phi, doublet sum rule, residual diagnostic, the
static-leading + dynamical-residual framing) as **Paper 20** (+ forward pointers into XV-XVII per sec.16),
rather than more calcs. The physics is largely settled in the notes; the gap is that it is not in the
papers. Do step 3's remaining calc (up-sector non-perturbative) AFTER/WITHIN that consolidation, from the
notes' route, not the stale papers.

## 46. CONSOLIDATION (2026-09-22): CLAUDE.md sec.12 added; Paper 20 DRAFT created
Per user: (1) added CLAUDE.md sec.12 "the notes lead the papers -- read notes/ before grounding in a
paper" (the notes supersede the papers on open threads; re-deriving from stale papers manufactures false
open problems, e.g. the sec.44 two-route error). (2) Drafted `papers/main_paper20.tex` (DRAFT, matching
XVIII/XIX format: theorem/prop/conj/remark/open-problems, bibliography). Paper XX = "The Dynamical Layer
of the Quark Sector: E_6-Native Generations, the Isospin Mass-Scale, and Static Residuals as Strong
Dynamics." Consolidates this conversation:
  - sec.groups: E_6-native(up)/E_8-bridge(down), phase-unit ratio 9/4=2.25 (Prop); generation=ribbon
    twist |Tw|={1,2,5} (Rem); the |Tw|^4.2 power law is emergent not fundamental (Rem, sec.33).
  - sec.isospin: tangent-orientation mass-scale M_dn/M_up=0.189 (Prop, sec.28); orientation-not-energy
    (Rem, sec.30).
  - sec.gen1: charge-linear Q*Phi flip (Prop, sec.26).
  - sec.sumrule: doublet sum rule C~13 (Conj, sec.34); retracted 4pi/phi^6/solid-angle/group readings of
    C (Rem, secs 38-40), only 1/3=SU(3) charge principled.
  - sec.residuals: the diagnostic -- static=leading (<=0.5%), dynamical=~alpha_s (~20-30%), ~40x split
    (Prop, sec.42); perturbative alpha_s already built + doesn't close (Rem, sec.43); sigma*R=hadron not
    current-mass (Rem, sec.44).
  - Open problems: up-type non-perturbative residual; C's magnitude (formation mass-balance); CS-edge/
    Hall dynamics (sec.41); Phi derivation; M_dn/M_up residual.
VERIFIED: env balance (15 thm-like begin/end, document/abstract/bib balanced), 0 backslashed underscores,
8 sections, all 10 cites have bibitems (Paper1 now cited). FLAGGED: Paper XVIII DOI is a PLACEHOLDER
(zenodo.XXXXXXXX) -- USER to fill; cross-paper \ref numbers left as %-comments per CLAUDE.md 11 (USER to
verify/number on compile). USER TO COMPILE (no pdflatex in-session). E_6 picture now IN a paper (draft) --
supersedes sec.16's "note-only" status once XX is finalised + forward pointers added to XV-XVII.

## 47. OPTION 1 (CS-edge fractional charge, 2026-09-22): DERIVED from the SU(3)_1 Z_3 center (triality)
Script: `src_paper18/cs_edge_fractional_charge.py`. Opened the CS-edge/Hall avenue (sec.41: fractional
charge as an edge/Hall, not static-writhe, quantity). REFINED the target: the fractional '3' is the COLOUR
SU(3)_1=(E_6)_1 Z_3 center (triality = the 3 crossings), NOT the SU(2)_3 lepton level k=3 (sec.41 was
mis-aimed at the lepton level).

DERIVATION (framework-native chain): T(2,3) -> 3 crossings -> SU(3)_1 colour (Papers XV-XVII) -> Z_3
center -> charge in thirds.
  - SU(3)_1 primaries: 1(triality 0,h=0), 3(t=1,h=1/3), 3bar(t=2,h=1/3), c=2.
  - topological spin exp(2pi i h) of the fundamental = e^{2pi i/3} = the Z_3 center element -> h=1/3 IS
    the triality charge (the conformal weight encodes the center).
  - Z_3-center consistency of SU(3)_colour (x) U(1)_em forces Y = triality/3 (mod 1). quark (t=1)->Y=1/3
    -> Q=I_3+Y/2 = {+2/3 up, -1/3 down} with I_3=+-1/2 the z-network isospin (sec.28). lepton (t=0) ->
    Y integer -> Q={0 nu, -1 charged}.
  => quarks fractional (carry colour triality), leptons integer (colour singlet) -- DERIVED from the
     trefoil's 3-fold (Z_3) topology. UPGRADES the static writhe conjecture P18:conj:isospin to a CFT
     derivation. Full charge = I_3 (z-network, sec.28) + Y/2 (SU(3)_1 triality, this sec); both trefoil-
     geometric.
CS-EDGE realisation [I, framing]: the SU(3)_1 CS edge modes carry the triality charge (bulk-edge
correspondence); charge fractionalises into e/3 via the Z_3 center. This is the DYNAMICAL picture of the
same center charge -- the "equation of motion" avenue of sec.41 -- but the substance is the group-theory
Z_3 derivation; the explicit edge-mode/Hall spectrum (guaranteed to agree by bulk-edge) was not computed
(would confirm, not add).
CAVEAT (CLAUDE.md 3): Y=t/3 center-consistency is the STANDARD SU(3)xU(1) quantisation (the SM Z_6
quotient). The FRAMEWORK-native input is that SU(3)_1 is DERIVED from T(2,3)/2T (Papers XV-XVII), so the
chain trefoil->3 crossings->Z_3->thirds is the framework's own; the charge quantisation given SU(3)_1 is
standard. So this DERIVES fractional charge WITHIN the framework's colour CFT, cleanly.
STATUS: Option 1 POSITIVE -- fractional charge derived from the colour Z_3 center, upgrading P18's writhe
conjecture; refines sec.41 (colour-3, not level-3). Feeds Paper XX open problem O3 (CS-edge charge): O3
is now largely ANSWERED (Z_3-center derivation), with the explicit edge spectrum the remaining technical
step.

COMMITTED TO PAPER XX (2026-09-22): added new section "Fractional Charge from the Colour Z_3 Center"
(P20:sec:charge) with P20:prop:frac_charge (Y=triality/3 -> Q={2/3,-1/3}; leptons integer) and
P20:rem:charge_cs (colour-3 not level-3; CS-edge realisation; standard SU(3)xU(1) quantisation caveat).
O3 rewritten from open avenue to "explicit CS edge spectrum" (charge now derived, edge spectrum is the
confirmatory remainder). Abstract result added as (3), renumbered; intro + summary updated. VERIFIED:
17/17 thm-likes balanced, 9 sections, abstract (1)-(6), 0 backslashed underscores, new labels defined.

PAPER XVIII UPDATE ASSESSMENT: NOT needed now. P18:conj:isospin assigns BOTH isospin (z-network->I_3) AND
fractional charge (2/3,-1/3) from writhe; Paper XX KEEPS the isospin part (uses I_3) and DERIVES the
fractional denominator (thirds) from the colour Z_3 center. CONSISTENT (I_3 + Y/2 = same charges), not
contradicted -- Paper XX supplies the deeper origin. A one-line POINTER on P18:conj:isospin (fractional
denominator derived from colour Z_3 in Paper XX; writhe gives isospin) is a FORWARD REF -> deferred to
finalisation per user (18/19/20 are drafts). Exact spot flagged: P18:conj:isospin / its scope remark.

## 48. DERIVE C via FORMATION (2026-09-22): constituent scale DERIVED; C reduced to the E_6 tower
Script: `src_paper18/derive_C_formation.py`. Attempted to derive C / pin the E_6 tower via the formation
mass-balance (user: these are ONE problem). Real progress + honest open remainder.

NEW RESULT [E, derived]: the CONSTITUENT quark scale = Lcond^baryon phi^6 exp(8pi) = 237 MeV, from
first principles: Lcond^baryon=T_CMB(pi^2/45)^{1/4} (sec.21), phi^6=phi^{2 Q_group^baryon}, and the 8pi
CS spoke (DeltaQ=3-1=2, quark from the Q_H=1 neutrino; Paper XI l.1051). Matches sec.20's EMPIRICAL
215.5 MeV to ~10% (= the alpha_s non-perturbative residual, sec.42; a T-matrix correction analogous to
the lepton's exp(-3/400) would close it). So the constituent scale, previously an input, is DERIVED.
The engine is the extra CS spoke: quark has exp(8pi) (DeltaQ=2) vs lepton exp(4pi) (DeltaQ=1).

REDUCTION of C: constituent/m_e = (Lb/Ll) phi^{-14} exp(4pi+3/400) = (10/3)^{1/4} phi^-14 exp(4pi+3/400)
= 463 (derived). Then C = (m_const/m_lep) * phi^{-2 Dn_q(g)}, so C is DERIVED up to the up-doublet tower
level Dn_q(g): back-out gives Dn_q(g1)=+3.69, Dn_q(g2)=-1.82 (gen2 negative -- charm 1270 > constituent
237, doublet above the constituent scale). Dn_q(g) is the E_6-native tower level.
=> deriving C == pinning Dn_q(g) == pinning the E_6 tower. CONFIRMS the user's insight: ONE problem.
The remaining open piece is Dn_q(g) (the E_6 tower level), blocked by the n-T_g degeneracy of sec.24.
C's magnitude ~13 is a DELICATE exp(4pi)/phi-tower residual (m_const/m_lep~463 times the tower
suppression), NOT a standalone constant -- consistent with sec.40 (C is a tower/exponent residual, not 4pi).

STATUS: Paper XX open problem O2 (C via formation mass-balance) PARTIALLY ANSWERED: (i) constituent scale
DERIVED (new, 8pi spoke); (ii) constituent/lepton ratio DERIVED (~463); (iii) C reduced to the E_6 tower
level Dn_q(g). What REMAINS = Dn_q(g), the E_6 absolute tower (sec.24 degeneracy). The constituent-scale
derivation is a clean new result worth adding to Paper XX (a proposition); the C-reduction updates O2.
NEXT (the true remaining knot): pin Dn_q(g) -- break the n-T_g degeneracy with an independent handle on
n or T_g. This is the one irreducible open computation for the up-sector masses.
COMMITTED TO PAPER XX (2026-09-22): P20:prop:constituent (m_const=Lcond^baryon phi^6 e^{8pi}~237 MeV;
m_const/m_e~463; C=(m_const/m_lep) phi^{-2 Dn_q}), P20:rem:C_residual; O2 rewritten to "pin the E_6 tower
level Dn_q"; abstract result (5) + summary updated. Verified: 19/19 thm-likes, 0 backslashed underscores.

## 49. EXTERNAL-LIT CHECK (user segue, 2026-09-22): 3 core structures ESTABLISHED + 1 mainstream program
Web search -- is the framework's CFT/topological content already established? Recorded fully in
[[external-literature-connections]]. Summary:
  (1) phi = SU(2)_3 WZW = FIBONACCI ANYON quantum dimension (established TQC). V*=phi is the Fibonacci
      qdim; the j0,j1 vs j1/2,3/2 split is the standard Fibonacci substructure. VALIDATES phi.
  (2) Solar angle theta_12=arctan(1/phi)=31.72deg (Paper XIX) is the ESTABLISHED A_5 (icosahedral)
      GOLDEN-RATIO MIXING (Everett-Stuart arXiv:0812.1057; modular A_5 King et al 2206.14869, JHEP 2019).
      NOT novel -> framework MUST cite it; its contribution is the TOPOLOGICAL origin, not the mixing.
  (3) MODULAR FLAVOR SYMMETRY (Feruglio 2017+): Yukawas=modular forms, hierarchy from tau near fixed
      points, finite groups incl. Gamma_5=A_5. The mainstream analog of the framework's CFT-for-masses +
      RG-fixed-point tower. Position the framework as its TOPOLOGICAL/condensate realisation.
  (4) Fractional charge (sec.47) = the standard SM Z_6=Z_3(colour)xZ_2xU(1)_Y quotient (Y=triality/3).
      Confirms sec.47's caveat: standard mechanism, framework-native part = SU(3)_1 from the trefoil.
NET: VALIDATION -- the framework's phi, icosahedral-golden-ratio mixing, and fractional charge are
MAINSTREAM (not ad hoc); it sits at TQC (x) modular-flavour (x) SM-Z_6, unified by a hopfion/McKay origin
+ gravity. TO-DO (finalisation): cite A_5 golden-ratio mixing + modular flavour symmetry; note phi=Fibonacci
qdim; verify which papers already cite these. Honest positioning: novelty = the topological UNIFICATION,
not the individual mixing/charge results (which are known).

## Epistemic ledger
[E] Calugareanu SL=Wr+Tw; WZW T-matrix = framing (standard CFT); colour=crossings, isospin=Z_2 (framework).
[N] derived this session: parity-gate anchor; mirror-pure isospin; colour!=generation; the 3-orthogonal
    separation; generation = ribbon twist (proposal, motivated by T-matrix = framing).
[O] open: lift I_3-inheritance to a derivation; reconcile the two isospin mechanisms; verify the
    exponent<->twist map (this step); does the mirror = orientation/twist reversal.
[E] (2026-09-21) n_q^base=6 via the CMB arithmetic (Lcond^(baryon)=T_CMB(pi^2/45)^{1/4},
    rho_CMB=3 Lcond^4), parallel to P17:prop:qh1_cmb; verified. m_e recovered 0.2%.
[O] (2026-09-21) the profile-normalisation (saddle) half of n_q^base is BLOCKED -- no stable
    Q_H=3 saddle; and the u,c isospin tower-steepness (~2x, NOT the E_6 twist) needs a new
    I_3->tower coupling. See [[coset5-2t-in-2i-geometry]] for the parked coset-5/pentagon-5 lead.
[N] (2026-09-21, sec.22-23) I_3 lever WORKED: confined up needs extra gen-growing suppression ratio
    ~2.1-2.2/gen (two routes agree); H_T and H_tw ruled out. LEAD (sec.23, user) = GROUP asymmetry:
    down E_8-lepton-bridged (h=30) / up E_6-native (h=12), ratio h_E8/h_E6=2.5 vs observed 2.13,
    NO new coupling. Supersedes the electric-charge lead. Top = sector transition (up tower = u,c only).
[O] (2026-09-21) NEXT = compute the T(2,2)->T(2,3) formation-energy EM piece (Paper XVIII O6):
    tests charge^1 vs charge^2, breaks the n-T degeneracy, closes sec.18 O6. Highest-value step.
[N] (2026-09-21, sec.24) E_6-native phase pi Q/(k h)=pi/4 (Q=3,k=1,h=12); ratio to E_8-bridge pi/9
    is 9/4=2.25 ~ observed steepness 2.13 (6%) -> GENERATION-steepness explained, grounded, no new
    coupling. But naive E_6 mirror formula does NOT give the ISOSPIN doublet splitting (phases not
    constant, sign flips gen1<->gen2). Steepness=E_6 phase (solved-in-shape).
[E] (2026-09-21, sec.25) the EM Q^2 SELF-energy is NOT the sign-flip: +Q^2 makes up heavier (wrong
    sign) and ~0.5 MeV (too small). [This kills the SELF-energy only -- see sec.26 correction.]
[N] (2026-09-21, sec.26, CORRECTS sec.25) EM REOPENED via the LINEAR charge*field coupling -Q*Phi
    (P18:conj:isospin's external field): SIGNED (down heavier, right sign), Phi~alpha*Lam_const~1.6
    MeV, generation-INDEPENDENT -> flips ONLY gen1 (48% of m_u, negligible at gen2,3). This IS the
    gen1-seed flipper and MUST enter the E_6 tower pinning (EM-subtract before fitting). H_Q (Q^2)
    stays retired; H_linear (Q^1) is the live EM hypothesis.
[E] decomposition stands: mean(up/down)=P18 out-of-plane (0.189, 21%); generation-spread=E_6/E_8
    steepness (sec.24). Flip = topological up-heavier tower + EM-linear gen1 flip.
[O] undetermined: Phi sign (obs-fixed) + magnitude (O(1) on alpha*Lam); fold Q*Phi into sec.24 tower.
[N] (2026-09-21, sec.27) Q*Phi pinned: Phi~2.5-3 MeV (O(1)*alpha*Lam_const) EM-subtracts to flip gen1
    up-heavier; gen2,3 untouched. Absolute up E_6 tower still does NOT close (T_g/twist non-linear).
[I/N] (sec.27, GEOMETRIC ROOT, user) up-type = CROSSING network (z=+r_0); down-type = MIDPOINT/distal
    network (z=0) [P18:conj:isospin]. Q*Phi = the z-level field coupling. Higher T(2,n)=more crossings
    (energy-input -> higher Q_H).
[N] (2026-09-21, sec.29) INTRA-TRIPLET DOWN = twist power law m~|Tw|^p (PER-STRAND generation twist,
    colour SL=3 fixed), |Tw|={1,2,5}, p=4.22, R^2=0.9999, one amplitude (~m_d). Power law ROBUST.
[E] (2026-09-21, sec.32) per-strand twist Faddeev ENERGY scales ~Tw^{1.5-2.8} (E=K*J4), NOT Tw^4.2.
    Quark mass is ALGEBRAIC/saddle-less (sec.1) so it is NOT any Faddeev energy. Geometric/phi^6 route DEAD.
[E] (2026-09-22, sec.44) STEP 3 first calc: formation-energy(sigma*R) hypothesis FALSE for current-mass
    residuals (sigma*R~278 MeV constant vs Dm spanning 6 orders). sigma*R = HADRON scale (sec.19, done),
    NOT a current-mass correction. [sec.44's "two-route ambiguity" RETRACTED in sec.45.]
[E] (2026-09-22, sec.45, user) ROUTE IS RESOLVED: E_8-bridge(down)/E_6-native(up), sec.16/24; Paper V
    pure-E_8 SUPERSEDED. sec.44 misread the BEHIND-papers. KEY: **papers are well behind the notes** --
    all of this session (E_6/E_8 split, isospin scale, Q*Phi, sum rule, residual diagnostic, static-
    leading framing) is NOTE-ONLY. RECOMMENDATION: CONSOLIDATE into Paper 20 (physics settled in notes;
    gap = not in papers) before more calcs. Step 3 target = UP-sector E_6-native residual (notes' route).
[E/O] (2026-09-22, sec.43) STEP 3 SCOPE -> Paper 20: perturbative alpha_s is DONE (Paper V/XVI derive
    alpha_s(LUV)=C_F sin(pi/5)/(2 phi^{43/7})~0.0204, RG-run) and does NOT close residuals (P16: pattern
    of undershoots PERSISTS). => residuals are NON-PERTURBATIVE formation/confinement (metastable quark,
    no saddle). Step 3 = the formation energy (sec.18/19, needs topology-protected solver) + confinement
    kinematics (c_s=c/phi) + isospin/crossing (sec.30). FIRST CALC: test analytically whether a sigma*R
    +barrier formation correction reproduces P16's surviving undershoot pattern (no blocked solver).
[E] (2026-09-22, sec.42) RESIDUAL DIAGNOSTIC (step 1, CONFIRMS user): residuals split ~40x by sector --
    static-saddle (EM/charged-lepton) ALL <=0.5% (geomean 0.011%); dynamical/loop (quark/neutrino) ALL
    >=20% (geomean 33% ~ alpha_s); cosmology (mixed) between. Static=EXACT leading order; residuals ARE
    the dynamical layer, size ~ sector coupling (alpha_em tiny vs alpha_s~0.3), largest where no static
    saddle (quark). => target the DYNAMICAL LAYER (strong/formation dynamics, alpha_s scale), not
    individual constants. Source: Paper VII table + session.
[E] (2026-09-22, sec.48) CONSTITUENT quark scale DERIVED = Lcond^baryon phi^6 exp(8pi) = 237 MeV
    (215.5 to ~10%=alpha_s); 8pi CS spoke (DeltaQ=2) is the engine. C = (m_const/m_lep) phi^{-2 Dn_q(g)},
    m_const/m_lep~463 derived; C reduced to the E_6 tower level Dn_q(g) (open, sec.24). deriving C ==
    pinning the E_6 tower (user's insight confirmed). Paper XX O2 partially answered. Remaining knot =
    Dn_q(g) via breaking the n-T_g degeneracy.
[N] (2026-09-22, sec.47) FRACTIONAL CHARGE DERIVED from the COLOUR SU(3)_1 Z_3 center (triality): quark
    Y=triality/3=1/3 -> Q=I_3+Y/2={2/3,-1/3}; leptons triality-0 -> integer. Chain: T(2,3)->3 crossings
    ->SU(3)_1->Z_3->thirds (framework-native; SU(3)_1 derived Papers XV-XVII). Upgrades P18 writhe
    conjecture. REFINES sec.41: the '3' is COLOUR SU(3)_1 (triality), NOT SU(2)_3 level. h=1/3=triality/3
    (topological spin = center phase). CS-edge = the dynamical realisation (bulk-edge), not explicitly
    computed. Answers Paper XX open O3 (charge). CAVEAT: Y=t/3 is standard SU(3)xU(1) quantisation.
[O] (2026-09-22, sec.41, user) OPEN AVENUE: CS DYNAMICS (EOM/induced current/Hall/edge modes) is
    UN-EXPLORED (grep: no edge/Hall/EOM content). Framework does charge STATICALLY (writhe conj, S-matrix
    alpha). Fractional charge is an EDGE/HALL phenomenon; SU(2)_3 k=3 -> nu=1/3 FQH edge charge e/3=down,
    2e/3=up (level-3 -> third-charges, framework-native, DYNAMICAL not static). The "+1/3" may be nu=1/k
    not conformal weight. CAVEAT: SU(2)_3 non-abelian vs abelian nu=1/k -> the EM Hall response (=the
    un-built CS-EOM) must be worked out. NEXT: condensate CS EOM/induced current/edge modes; does frac
    charge emerge dynamically + does edge-current energy feed the mass? Live avenue -- do NOT close.
[E] (2026-09-22, sec.40) CONFLATION resolved: Hopf charge enters mass EXPONENTIALLY (e^{4pi*DeltaQ} +
    tower); m_e~exp(4pi) because mass = exp(-S) WZW boundary-state amplitude, 4pi = the CS ACTION. So a
    mass RATIO between Q_H sectors ~ e^{4pi} (huge), not ~13. C~12.9 is a tower/exponent RESIDUAL (the
    e^{4pi} cancels, being common to both scales); "C~4pi" compares a mass ratio to a BARE ACTION =
    category error/coincidence. DROP the 4pi reading. Sum rule empirical; only +1/3=SU(3) charge (linear
    EM) principled. Explains why secs 38/39 found no mechanism (none possible).
[E] (2026-09-22, sec.39) O6 solid-angle COMPUTED (hopfions Q_H=1,2,3, charge FFT-verified): director
    geometry across Q_H=2->3 grows by J2 1.21, J4 1.80, solid-angle 1.33 -- NONE ~4pi/12.9. Solid-angle
    4pi hypothesis FALSIFIED. C~12.9 has NO mechanism (CS exp too big, vac-pol too small, solid-angle
    too small, group-13 unmotivated). C EMPIRICAL; only +1/3=SU(3) charge principled. O6/4pi chase DONE.
[E] (2026-09-22, sec.38) the framework's 4pi's do NOT give C's factor-of-4pi: CS spoke 4pi*DeltaQ is
    EXPONENTIATED (exp(4pi)=2.9e5); vac-pol 1/(4pi) per WZW level (Paper IV) is a ~8% correction; tower
    is phi^2/level. So C=4pi+1/3 has NO derived mechanism among them; a factor-4pi would need a solid-
    angle/phase-space 4pi (not established). 1/3=SU(3) charge SOLID (sec.37); ~13 magnitude UN-derived.
    Do NOT re-chase the vac-pol 1/(4pi) route. O6 must give ~12.9 from formation GEOMETRY.
[N] (2026-09-22, sec.35+36) DERIVING C: uncertainty-weighted, g2(c,s) is the reliable anchor (C=12.90
    +-0.19), g1(u,d) unconstraining (+-1.1). Reliable data DISFAVOURS pure 4pi (-1.8sig), lands on
    4pi+1/3=12.90 (h_SU(3)_1 quark spin; -0.02sig), 13 within +0.5sig. C = 4pi + quark-correction
    (~1/3 or pi/9): user's 4pi base vindicated + a quark quantum number. POST-HOC + 1 reliable point +
    correction not unique -> LEAD, not derived. O6 formation mass-balance is the arbiter.
[N] (2026-09-21, sec.34) DOUBLET SUM RULE (user): (m_up+m_down)=C*m_charged_lepton, C~13 for the 2
    CONFINED doublets (13.37, 12.90; 3.5%). Cross-predicts m_c to 3.9%, m_u to 11% from leptons+downs.
    Resolves sec.17: UP IS lepton-anchored via the DOUBLET budget (m_up=C*m_lep-m_down), not its neutrino
    partner. Top breaks it (bare transition). CAVEAT: 2 points, C not derived (~4pi or 13=anchor exp);
    product/seesaw fails, only SUM works. Strong LEAD for the open UP intra-triplet.
[E] (2026-09-21, sec.33) REP-THEORY VERDICT: m~|Tw|^4.2 is EMERGENT, NOT a WZW identity. Per-step p
    scatters 1.3% (4.167 vs 4.223); the clean power REQUIRES measured leptons (rep-theory-internal gives
    4.03 vs 4.57, 13%); phi^3=4.236 is within the fit noise. FUNDAMENTAL = the E_8 route (m_lep x
    exp(n_q pi/9), sec.24); |Tw|^4.2 is a re-parametrization. Twist->mass thread CLOSED: down = E_8 route,
    no exponent "meaning" to derive.
[E] (2026-09-21, sec.30) Faddeev energy density is EQUAL on crossing(up) and midpoint(down) networks
    (quartic 0.94, quadratic 1.00; E(+r0)=E(-r0) exact). => up/down is NOT a static energy difference;
    it is tangent-ORIENTATION x external field (confirms sec.28 + P18:conj:isospin premise). CLOSES the
    static-energy route for UP: neither this nor the jam MD cracks UP by energy -> UP = field-selection
    + formation dynamics, no static shortcut.
[N] (2026-09-21, sec.28) ISOSPIN SCALE CLOSED analytically: P18 out-of-plane tangent fraction gives
    M_up/M_down=5.29; residual to measured 6.38 = 1.207 ~ the uniform ~1.2 running (flagged). Geometry
    +running, ZERO free params. Working mechanism = TANGENT ORIENTATION (in/out-of-plane exit), NOT
    curvature or spatial proximity (min-dist ratio 1.08) -> "shared crossings" is TOPOLOGICAL, not
    spatial (corrects sec.27 wording). Open = intra-triplet absolute masses (E_6 tower sec.24) + gen1
    seed (Q*Phi sec.26). Option B jam-MD only needed for the intra-triplet, NOT the isospin scale.
