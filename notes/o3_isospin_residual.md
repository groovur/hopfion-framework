# Paper XX O3: the M_dn/M_up ~1.2 residual (2026-09-29)

Prediction (tangent geometry, P20:prop:isospin_scale): M_dn/M_up = 0.189. Measured: 0.157. Factor ~1.2.

## RG-RUNNING vs DYNAMICAL decomposition of the down-type ~1.2 (2026-09-29) -- caption corrected
Question (user): is the c/b or generation structure affected by the u/d breathing work? Answer: NO to both;
plus a decisive decomposition of the down-type ~1.2 over-prediction. KEY THEORY FACT: the mass anomalous
dimension is FLAVOUR-UNIVERSAL, so RG running preserves ALL mass ratios -- intra-triplet (m_s/m_d, sqrt-phi
ladder, e^{-2pi}) AND M_dn/M_up (ratio of geomeans). => running CANNOT cause any ratio residual.
Numbers (all brought to common 2 GeV where needed):
- M_dn/M_up over-pred = down_geomean(1.145)/up_geomean(0.939) = 1.220 = breathing 1/0.811=1.233 (1%). RG-INV
  => purely dynamical (confirms O3 "running cancels in the ratio" rigorously).
- ANTI-CORRELATION is the tell: down OVER-pred (1.145), up UNDER-pred (0.939). Flavour-universal running moves
  both the SAME way => cannot produce opposite deviations => the asymmetry is NOT running, it is the network
  (breathing) layer; sign right (crossing/up stiffer => up heavier => up under-pred, down over-pred).
- Genuine running is small+specific: b excess 1.21@m_b -> 1.05@2GeV (~15% was b-specific SCALE-QUOTING, not
  generation-universal); common all-6 offset = sqrt(1.145*0.939)=1.037 (~4%, anchor calibration).
- intra-triplet SHAPE (RG-inv, dynamical): m_s/m_d pred/meas=0.90, m_b/m_s=0.93 (the sqrt-phi/n_e=20 physics).
CONCLUSION: down-type ~1.2 is NOT "generation-universal RG running" -- it is the confined DYNAMICAL layer
(breathing), RG-invariant, = the isospin residual. They UNIFY (paper's thesis). Genuine running = only b's
scale-quoting (~15%, b-specific) + ~4% common offset. GENERATION STRUCTURE + c/b individual routes UNAFFECTED
by the breathing work (intra-triplet ratios RG-invariant; breathing is generation-universal, cancels within a
triplet; c via E_6-native, b via E_8-bridge, neither routes through M_dn/M_up). EDIT: Paper XX table caption
(l.733) corrected -- down/up asymmetry = dynamical layer (RG-inv, cites P20:rem:isospin_breathing), only b
scale-quoting + ~4% common as genuine running; c,b deviation line aligned (c ~13%, b ~5% at common 2 GeV).

## FINAL PINNED RESULT (2026-09-29) -- consistent transverse mode, FOLD INTO PAPER XX
The confined-mass mode is the TRANSVERSE oscillation restored by the tube tension. Only the Frenet-NORMAL
(Nhat) component has string restoring (dL/dA = -INT w kappa (Nhat.dhat) ds; the binormal Bhat gives dL/dA=0,
a flat/unconfined direction). So the mass mode = Nhat oscillation, and F AND mu must both be Nhat (consistency:
the earlier 0.147 wrongly mixed the transverse FORCE 1.375 with the DILATION inertia 0.886 -- different modes).
Both computed for the Nhat mode:
  - FORCE  F_up/F_down = 1.375  (`pin_displacement_direction.py`; = windowed-curvature ratio, direction-INDEP
    in the transverse plane; crossing tube more curved => stiffer string confinement => heavier up).
  - INERTIA mu_up/mu_down = 1.009  (`transverse_mode_inertia.py`; CONVERGED N=48/72/80 = 1.012/1.009/1.010;
    a rigid transverse translation changes the field nearly network-independently, unlike the dilation 0.886).
  - Airy (linear confinement, V=sigma*L): M_dn/M_up mult = 1/(F^2/mu_ratio)^(1/3) = 1/(1.375^2/1.009)^(1/3)
    = 0.811.  => M_dn/M_up = 0.189 x 0.811 = 0.1533  vs measured 0.157  => 2.4% (residual-on-residual 1.023).
CLEAN READ: the residual is ESSENTIALLY the tangent/CURVATURE geometry alone (mu~1 => M_dn/M_up mult ~ F^{-2/3}
= 1.375^{-2/3} = 0.810); the inertia is a ~1% correction. ONE geometric quantity (windowed curvature ratio),
the SAME tangent orientation that sets the leading 0.189, Airy-quantized against the jam wall. No free direction
(force direction-indep; mode = the string-confined Nhat). This is the pinned, mechanism-first coefficient: 0.153.
DIRECTION PINNED + INERTIA MATCHED => COHERENT => folding into Paper XX O3 with the data (user, 2026-09-29).

## TWO-MODE COUNTER-FORCE -- COMPUTED, OVERSHOOTS (`two_mode_confinement.py`, 2026-09-29) [NEGATIVE/partial]
Operationalised the counter-force lead: total confined transverse zero-point = E_Nhat(Airy, tension) +
E_Bhat(harmonic, Faddeev bending), single knob sigma. Analytic target: closing 0.153->0.157 needs only a 6.6%
bending admixture (r=E_B/E_N=0.066). DIRECT computation: gN ratio 1.375 (confirms N-hat); but k_Fadd|Bhat
up/down = 0.084 (DOWN ~12x stiffer in bending) and d2L/dA2 down=8.8 vs up=1.0 -- BOTH driven by the down
near-inflection (rigidly displacing the near-straight inner-midpoint rebuilds the field enormously). So the
bending mode is NOT a gentle 6.6% counter-term: it is STIFF and down-heavy, giving correction~3.4 => M_dn/M_up
~0.65 (down far too heavy), no physical sigma reaching 0.157 (E_B~52 >> E_N~0.55 at sigma=1; k_Fadd~3e5-3.8e6
field units >> any reasonable tension force sigma*gN~1). Adding the mode RUINS the pure-Nhat 0.153, not improves.
HONEST READ: (1) the counter-force is REAL and RIGHT-SIGNED -- the down near-inflection carries a transverse
(harmonic) stiffness the linear-curvature force misses, pushing M_dn/M_up UP toward 0.157 (user's intuition
correct). (2) Its MAGNITUDE is not a clean 6.6%: computed directly it OVERSHOOTS (~10x too big) because the
near-inflection makes the down quadratic stiffness ~10x the up. (3) The true size depends on sigma (cutoff-set;
Faddeev E scale-invariant => tension NOT derivable) AND on which quantisation wins at the inflection (Airy weak-
linear vs harmonic strong-quadratic, itself sigma-dependent). => the lead IDENTIFIES the right mechanism but
CANNOT pin the coefficient without the cutoff. Naive addition overshoots; honest estimate stays SINGLE-MODE
N-hat Airy = 0.153 (2.5%), the 2.5% genuinely cutoff-limited, NOT a solved 0%. CAVEATS: tension-only-quadratic
B-hat variant (drop k_Fadd, use sigma*h_B) untested but h_B is ALSO down>>up (8.8 vs 1.0), same overshoot
direction; the huge down k_Fadd may be partly a rigid-displacement/near-inflection artifact. DID NOT keep fitting
B-hat variants to hit 0.157 (would be fit-for-wrong-reason, CLAUDE.md sec.3). Paper XX O3 already carries the
honest 0.153/2% framing + "higher-order jam dynamics" for the remainder -- consistent, NO paper change needed.

## ALGEBRAIC FRENET DATA + THE MISSING COUNTER-FORCE (2026-09-29, user lead: "oscillating countering force")
Moved the force computation to algebra (sympy, exact). |G'|^2 = 9 r0^2 cos^2(3t) + 4(R0+r0 cos3t)^2 (closed).
Exact curvature/torsion at the site types (R0=3, r0=sqrt2/phi):
  crossing (up, 3t=pi/2):   kappa=0.399  tau=-0.168  |G'|=6.55   [kappa=sqrt(64R0^4+900R0^2 r0^2+2025 r0^4)/(4R0^2+9r0^2)^{3/2}]
  distal lobe (down,3t=0):  kappa=0.349  tau=+0.014  |G'|=8.18
  inner-mid (down,3t=pi):   kappa=0.026  tau=+11.28  |G'|=5.00   <-- NEAR-INFLECTION (kappa~0), torsion spike.
KEY: the DOWN network is HETEROGENEOUS -- one of its sites (inner-midpoint) is a near-inflection where the
N-hat TENSION force (F=sigma INT w kappa (Nhat.dhat) ds, PROPORTIONAL TO CURVATURE) nearly VANISHES. So the
pure-Nhat model treats those down quarks as too soft => UNDERESTIMATES M_down => M_dn/M_up too LOW (0.153<0.157)
-- exactly the sign of the 2.5% miss. Those quarks are still confined, by the force the curvature-term omits:
the Faddeev BENDING stiffness (Bhat mode, k_up/k_down=0.796, local_bending) + the TORSION coupling that mixes
Nhat<->Bhat (tau=11.3 at the inflection). That restoring RAISES M_down => raises M_dn/M_up toward 0.157 (right
dir). => THE MISSING FORCE (user's intuition, confirmed): the SECOND transverse mode (Bhat, bending-restored,
NOT tension) + torsion N-B mixing, localized at the down near-inflection where the leading curvature term is
blind. MAGNITUDE is CUTOFF-SET (sigma/K_Faddeev ratio, not derived -- same status as accel scale / strain mag),
so directionally derived, coefficient open; the 2.5% remainder is consistent with it. Algebra HELPED: exact
resolution-free force functional + revealed the near-inflection a single numeric 1.375 hid + gave tau(t) = the
N-B coupling carrier. NOTE the code's Bishop frame handles the inflection smoothly (numeric 1.375 trustworthy);
Frenet tau=11.3 is the frame's near-inflection blow-up but flags the real geometry. NEXT (optional): two-mode
zero-point E_Nhat(Airy)+E_Bhat(harmonic) with sigma/K_Faddeev as the one cutoff parameter -> does it hit 0.157.
CORRECTION to earlier status text below: the pinned consistent result is 0.153 (~2.5%), NOT 6% (6% = superseded
radial-3D-with-dilation-inertia mismatch).

## STATUS (2026-09-29) -- SETTLED SO FAR, NOT YET FOLDED INTO PAPER XX (user: "no fold in yet")
Settled position: (1) DF is OUT as the residual mechanism, both statically (isospin_feedback_residual, ~1%)
and dynamically (df_breathing_potential: raw ~1.01, per-arc wrong-direction). (2) The confined-mass MODE is
AIRY-against-the-jam-wall -- pinned from first principles (df_breathing_potential: E_fb(lambda) monotonic =>
DF+string both drive to collapse => core/packing wall sets the size). (3) The up/down differential is carried
by the STRING TANGENT GEOMETRY (F_up/F_down~1.28), the SAME tangent orientation that sets the leading 0.189.
DIRECTION NOW PINNED (df below, `pin_displacement_direction.py`): the confined oscillation is TRANSVERSE to
the tube (longitudinal = reparametrisation, zero to 1st order). The tension force dL/dA = -INT w kappa (Nhat.dhat)
ds projects the local CURVATURE onto the displacement, so F_up/F_down = the windowed-curvature ratio = 1.375
INDEPENDENT of the transverse direction (Nhat vs Bhat: sweep is FLAT at 1.375 for all psi) -- there is NO
direction freedom once the unphysical tangent part is dropped. => M_dn/M_up mult = 0.777 (pinned) => 0.189 x
0.777 = 0.147 vs measured 0.157 = 6% OVERSHOOT (over-corrects). The earlier radial-xy "2% match" (0.815) was
an ARTIFACT: radial-xy carries a 14% TANGENT (reparametrisation) fraction on the crossing network that
asymmetrically dilutes its curvature projection; the physical (transverse) number is 0.147, NOT 0.154.
So (user's principle, CLAUDE.md sec.3/9): the mechanistically-correct direction gives the WORSE number
(0.147, 6%) and that is the honest one; the closer 0.154 was fit for the wrong reason (tangent contamination).
Paper XX O3 stays on the honest "dynamical residual" framing until we decide to fold; the pinned coefficient is
0.147 (6% overshoot), mechanism = Airy-against-jam-wall + tangent/curvature geometry + breathing inertia 0.886.

## DIRECTION PIN (`pin_displacement_direction.py`, 2026-09-29) -- the residual coefficient IS pinnable
Feared the force ratio would be wildly direction-sensitive (curvature projection). IT IS NOT: because dL/dA
couples only to the Nhat-component (cos psi) of the displacement, the RATIO F_up/F_down cancels the cos psi
and equals the windowed-curvature ratio 1.375 for EVERY transverse direction (sweep flat, range [1.37,1.37]).
Named directions: radial-xy 1.278 (but 14% tangent-contaminated -> 0.816); radial-3D 1.373; Frenet Nhat 1.375;
radial-xy_|_ 1.359; radial-3D_|_ 1.455; all transverse-physical ones cluster 1.36-1.46 -> M_dn/M_up mult
0.748-0.783, tight around 0.777. kappa windowed: crossing 0.384 / midpoint 0.256 (ratio 1.5). CAVEAT: inertia
0.886 was the dilation/radial-xy value (local_bending 0.885, breathing 0.886 -- robust); the transverse-mode
inertia not separately recomputed (expected ~0.886). RESULT: direction pinned, coefficient 0.777 -> 0.147,
a genuine 6% (=0.157/0.147=1.068) residual-on-residual = the higher-order confined/jam dynamics + Airy approx.

## DENSITY FEEDBACK — TESTED, RULED OUT (`src_paper16/isospin_feedback_residual.py`)
Hypothesis (user): the sharp CROSSING (up) network saturates the feedback while the gentle MIDPOINT
(down) stays responsive, so the feedback-softened energy differs between them = the residual. TEST
(N=56, Construction-C trefoil): bare quadratic g2 up/down = 1.006 (matches faddeev_energy_crossing_vs_
midpoint.py's ~1.00); feedback-softened g2/(1+beta*g2) up/down peaks at 1.014 (beta=0.1), 1.008 (0.35),
1.006 (0.452) -> at most a ~1% differential -> M_dn/M_up multiplier 0.99-1.01. TARGET = 0.157/0.189 =
0.83 (a 17% shift). MISSED BY ~17x. RULED OUT.
REASON: g2 (field gradient) is dominated by the TUBE WALL (poloidal profile, same C*=2.5062 on both
networks), NOT the centerline curvature (kappa 0.399 crossing vs 0.026 midpoint). So the two networks
have near-identical g2 distributions and saturate the feedback almost equally -- which is also why the
bare energies came out equal. CAVEAT: N=56, crossing core somewhat under-resolved (could raise its
saturation a little), but the bare-sum equality is the tell; it is a ~1% effect, not 20%.

## THE CORRECT FRAMING (user, 2026-09-29) -- structural, not a missing mechanism
The measured quark masses are NOT the isolated-static object we compute. We compute the ISOLATED, STATIC
trefoil (WZW/tangent ideal) -- which for the quark does not even exist as a stable object (it is the
saddle that collapses in free space, the feedback-saturation result). The MEASURED masses are confined-
QCD constructs: current-quark parameters at a renormalisation scale, in a scheme, extracted only from
quarks bound in hadrons/jets. So the 0.189->0.157 gap (and every ~10-20% quark residual) is the
STRUCTURAL distance between "isolated topological ideal" and "confined, running, current-mass value" --
expected BY CONSTRUCTION, size ~alpha_s (the dynamical layer). This is why EVERY static-isolated mechanism
(DF ~1%, RG running cancels in the ratio, curvature, bending) fails to produce it: the residual is not in
the static sector; it is the confinement/formation dynamics that turns the ideal into the shared-jam value
(P20 O1, Section residuals). The negatives are POSITIVE evidence that the residual IS the dynamical layer.
CLEAN MIRROR of the lepton: lepton static=physical (free, responsive feedback) -> matches to 0.013%;
quark static!=physical (confined, saturated) -> ~20% residual. Same feedback-regime split, residual side.
HONEST LIMIT: explains the residual in KIND + sizes it (~alpha_s), does NOT compute 0.157 (needs the
confined/formation dynamics). Paper XX O3 + P20:prop:isospin_scale sharpened to this framing (2026-09-29).

## STRETCHED-DYNAMICS LEAD -- the BREATHING mode differentiates up/down (2026-09-29) [best residual lead]
Prior work (`notes/quark_sector_confinement_cherenkov.md` sec.2-3): confinement is ELASTIC (string tension
"snaps back" the stretched tube); the confined quark OSCILLATES, E_q/E_c=c/c_s=phi (marginal, R-indep).
Magnitudes (string tension, coupling) flagged OPEN there. NEW piece built this session: the STRETCH/
breathing mode = the flat scale zero-mode of E=K*J4 (scale-invariant). Zero Faddeev stiffness -> restoring
force is the external STRING; omega=sqrt(k_string/mu). The Faddeev-computable, chi-independent-in-ratio
piece is the breathing INERTIA mu=INT|dn/deps|^2 dV (dilation deformation n(x/(1+eps))).
- `src_paper16/breathing_inertia_crossing_vs_midpoint.py`: mu_up(crossing)/mu_down(midpoint)=0.886
  +-0.002, CONVERGED across N=48,64,72,88 (0.885/0.884/0.888/0.886 -- not a resolution artifact)
  -- a ~11% up/down differential (vs density feedback's ~1%, RULED OUT above). RIGHT DIRECTION: crossing
  (up) has LOWER breathing inertia -> higher omega -> heavier up -> M_dn/M_up shrinks 0.189 -> toward 0.157.
- Magnitude: M_dn/M_up multiplier = (mu_up/mu_down)^{1/2}=0.94 (harmonic) or ^{1/3}=0.96 (Airy/linear).
  Right sign, ~HALF the needed 0.83 (6% of the needed 17%). The OTHER HALF is in the string STIFFNESS
  k (assumed equal here): full M_dn/M_up=sqrt(k_up mu_down/(k_down mu_up)); reaching 0.83 needs
  k_up/k_down~1.28, plausible since the crossing exits IN-plane (T_z=0) and the midpoint OUT-of-plane
  (T_z=+-0.5) so the string couples differently. SAME tangent orientation that sets the leading 0.189 ->
  unification hint (leading ratio + dynamical residual both from tangent geometry).
- STATUS [lead, CLAUDE.md sec.3/6]: right direction + half the magnitude from first principles (inertia);
  string-stiffness half UNCOMPUTED (cutoff-set). Caveats: N=64; omega~mu^-1/2 assumes harmonic (real
  potential is linear string); uniform-dilation breathing (local oscillation may differ more).
- RADIAL-BENDING mode checked + RULED OUT as the residual mode (`local_bending_stiffness.py`, N=64):
  displaced the crossing/midpoint networks radially (confinement dir), consistent same-mode k AND mu:
  k_up/k_down=0.796 (crossing SOFTER), mu_bend_up/mu_bend_down=0.885 -> omega_up/omega_down=0.948 ->
  M_dn/M_up multiplier=1.055 = WRONG DIRECTION. RESOLUTION: the radial-bending mode is a DIFFERENT mode
  from the breathing/mass mode. The mass residual is the SCALE mode (the flat zero-direction confinement
  stabilises; its ground state IS the confined mass), which is Faddeev-SCALE-FLAT (no wrong-way Faddeev
  stiffness) -> its physics is INERTIA (0.886, right dir) + STRING coupling. The radial bending is a
  Faddeev EXCITATION independent of confinement (not the mass); its wrong-way result says: do NOT identify
  the confined mass with radial rattling -- use the breathing/scale mode. [Note the prior E_q~pi hbar c/R
  cavity-mode picture is thus an EXCITATION-spectrum statement, not the ground-state mass.]
- STRING LOAD-PARTITION DONE (`string_load_partition.py`, pure curve geometry): displaced crossing/midpoint
  networks radially, measured how the arc length L responds. KEY: the confinement potential is LINEAR
  (V=sigma*L) -> the confined quark is in a LINEAR well -> AIRY, not harmonic. The Airy ground state
  E0 ~ (F^2/mu)^{1/3} depends on the LINEAR tension FORCE F=sigma*dL/dA, NOT the curvature d2L/dA2. The
  earlier harmonic treatment (curvature stiffness, wrong direction 1.055) was the WRONG physics for a string.
  Linear force ratio: F_up/F_down = dL/dA crossing/midpoint = 1.278 (radial-xy) / 1.372 (radial-3D) -- the
  ~1.28 needed. Crossing (inner radius, in-plane tangent) stretches the string MORE per radial displacement.
- **RESIDUAL ~CLOSED (to ~2-6%)**: breathing mechanism = inertia mu_up/mu_down=0.886 + string force
  F_up/F_down=1.28, Airy: E0_up/E0_down=(F_up^2 mu_down/(F_down^2 mu_up))^{1/3}=(1.638*1.129)^{1/3}=1.227
  -> M_dn/M_up correction = 1/1.227 = 0.815 (radial-xy) / 0.778 (radial-3D). OBSERVED residual = 0.157/0.189
  = 0.831. So 0.189 x 0.815 = 0.154 ~ measured 0.157 (2%); radial-3D 0.147 (6%). The ~1.2 residual IS the
  confined breathing-oscillation (Airy: string tension force + stretch inertia), NOT a static effect.
- STATUS [strong lead, NOT established -- CLAUDE.md sec.3/6/9]: gives the FULL residual magnitude (not just
  half), right direction, from two computed quantities (inertia converged; force geometric). MODELING
  CHOICES (the uncertainty): (i) displacement direction (radial-xy 0.815 vs radial-3D 0.778, 2-6% spread);
  (ii) AIRY-dominance = the string linear term dominates the Faddeev harmonic term (k_Faddeev=0.80 computed;
  a mix would sit between Airy and harmonic; holds only if sigma large enough -- cutoff-set); (iii) the
  confined mode = local radial oscillation (a choice); (iv) N=64 inertia/force resolution. So: "captures the
  residual within the modeling uncertainty," not "derives 0.157 exactly." Strongest handle yet on O3.
- NEXT to solidify: (1) pin the displacement direction / mode physically (radial-xy vs 3D vs the actual
  confinement pull); (2) the Airy-vs-harmonic mix (does sigma dominate k_Faddeev at the confined scale?);
  (3) resolution. If (1)-(2) hold, O3 is mechanistically closed.
## DF-AS-BREATHING-POTENTIAL -- TESTED DIRECTLY, RULED OUT + PINS THE MODE (`df_breathing_potential.py`, 2026-09-29)
User's mechanism-first question: does density feedback + the breathing COMBINE to give the residual (even if
the number is worse than the ad-hoc Airy)? The insight is CORRECT in principle: the bare E=K*J4 is scale-flat
but DF BREAKS the scale-invariance EXACTLY. Under a dilation the DF energy reduces (derived, checked beta=0=>flat):
    E_fb(lambda) = J4_0*J2a_0 + mu*J4_0*I(beta/lambda^2),   I(x)=INT g2/(1+x g2) dV  (field built ONCE at lambda=1).
So DF is a REAL breathing potential, not a bolt-on. Computed it explicitly, crossing(up)/midpoint(down) split,
beta=0.1..10 (responsive->saturated, sat frac 0.66->0.99), N=48 and N=64 (stable):
- RAW up/down DF-potential slope = 1.01 (NEARLY IDENTICAL) -- same as the static g2 (~1.006, isospin_feedback_
  residual). Reason unchanged: g2 is TUBE-WALL-dominated (same C* both networks), not centerline-curvature. So
  DF does NOT differentiate the networks DYNAMICALLY either -- confirms the static ruling now for the breathing.
- Per-arc-normalised (3 crossings vs 6 midpoints) the DF differential is ~0.60 and goes the WRONG DIRECTION:
  M_dn/M_up multiplier = 1.17-1.37 (>1; harmonic 1.17-1.22, Airy 1.34-1.37), target 0.831. Either normalisation:
  DF is NOT the residual mechanism. RULED OUT dynamically (not just statically). [normalisation subtlety = why
  the number is unreliable; the RAW ~1.01 is the tell, as in the static case.]
- WHAT IT DID PIN (the real payoff -- decides the MODE from first principles instead of assuming it): E_fb(lambda)
  is MONOTONIC increasing at ALL beta => DF pushes to small lambda, AND the string sigma*L(lambda) pushes to small
  lambda => BOTH drive to collapse => the confined size is set by the CORE/PACKING (JAM) WALL, and the breathing is
  AIRY-AGAINST-THE-WALL. This VINDICATES the Airy container (not a smooth harmonic well) from first principles, and
  identifies the wall as the shared-boundary/jam (the confined problem we keep circling). [In the SATURATED quark
  regime I(beta/lambda^2)~lambda^2/beta => E_fb~lambda^2, harmonic-like, but still monotonic->wall.]
- CONCLUSION: the residual mechanism is the STRING TANGENT GEOMETRY (F_up/F_down=1.28, string_load_partition) --
  a pure tangent-ORIENTATION effect, the SAME tangent geometry that sets the leading 0.189 -- NOT a DF effect. DF
  is now cleanly OUT both statically AND dynamically; the mode is Airy-against-the-jam-wall (pinned, not assumed);
  the tangent-geometry force carries the differential. Mechanistically honest: one geometry (tangent orientation)
  for both the leading ratio and the residual; DF governs the LEPTON (responsive) not the residual sign of the quark.

## Prior brainstorm routes: Route 1 (string) + Route 3 (quantize breathing) here; Route 4 (decay) = the
## string-break threshold sigma*L(lambda_max)=m(Q_H=2)+m(Q_H=1) = the well's finite depth, for free.
