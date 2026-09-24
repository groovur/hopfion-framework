# PLAN: fabric equation of state & recombination CDM (research program)
STATUS: speculative research plan. Goal: decide option 4 (fabric supplies its own recombination CDM,
no separate species) by deriving the equation of state (redshift scaling) of the incoherent condensate
director stress -- done properly in FABRIC terms, not by importing metric-expansion assumptions.

## THE PIVOTAL FORK (decide this first -- everything downstream depends on it)
"Redshift = space stretching" + "space = the fabric" => expansion is something happening TO the fabric.
Which of these is the framework's fabric doing under expansion?
  (F1) CREATED: new condensate fills the growing volume => intrinsic density ~ const (DE-like). 
  (F2) STRETCHED: fixed amount of fabric thinned => intrinsic density ~ 1/a^3 (matter-like) or per-volume.
  (F3) METRIC-ON-FIXED-FABRIC: fabric fixed, the induced metric expands (standard FRW on the condensate).
This choice sets the EOS of the fabric's OWN stress, hence whether it can be CDM. It is the crux and must
be pinned from the condensate action, NOT assumed. (The paper derives Einstein eqs from the action --
Phase 0 checks which of F1/F2/F3 that derivation actually implies.)

## Phase 0 -- Foundations (get the definitions right; PREREQUISITE)
0.1 Read the existing cosmology derivation (P7:thm:einstein, frozen-attractor de Sitter, inflation a(t))
    and determine which of F1/F2/F3 the condensate action implies. Does the framework CREATE or STRETCH fabric?
0.2 Define redshift geometrically: light = fabric ripple; fabric stretch -> wavelength stretch. Check it
    reproduces z=a_0/a-1. Confirm the scale factor a's physical meaning in fabric terms (spacing/density/vol).
0.3 Define temperature geometrically (from prior brainstorm): T = incoherent stored fabric stress. Relate
    T ~ 1/a (cooling) to the fabric releasing stored stress as it stretches (CMB = released stress?).
    GATE 0: is the fabric picture self-consistent with the framework's own FRW/Einstein derivation? If a
    contradiction appears, that is itself a key result (the fabric picture needs the action reworked).

## Phase 1 -- EOS of the three components, DERIVED in fabric terms (the gate)
1.1 Radiation (fabric ripples/photons): derive rho ~ 1/a^4 from ripple stretch + number dilution. (recover standard)
1.2 Matter (fabric knots/Hopfions, baryons): derive rho ~ 1/a^3 from fixed rest-mass + 1/a^3 number. (recover standard)
1.3 THE OPEN ONE -- incoherent director stress (the CDM candidate): its EOS depends on F1/F2/F3 (Phase 0) AND
    on localisation. Sub-questions:
      - is the disordered pre-ordering director stress attached to conserved structures (baryons) => 1/a^3, matter?
      - delocalised/vacuum => const, DE (feeds Lambda, fails as CDM)?
      - relativistic fluctuation => 1/a^4, radiation (fails)?
    RULED OUT earlier: it is NOT the baryon-Hopfions' internal deformation (chameleon-screened, view (a)/Paper XIII);
    it must be the fabric DIRECTOR field's own incoherent stress. Derive its dilution law from Phase 0's F-choice.
    GATE 1 (go/no-go for option 4): matter-like (1/a^3) -> option 4 viable; radiation/DE -> option 4 fails,
    revert to option 3 (recombination CDM is an unexplained Omega_DM input).

## Phase 2 -- Recombination CDM test (only if Gate 1 = matter-like)
2.1 Magnitude: does the fabric stress give effective Omega ~ 0.26 at recombination? (from rho_cond, deformation scale)
2.2 Clustering: do the fabric-stress PERTURBATIONS track density (non-oscillating gravitating wells)? (mechanism
    exists via density-dependent deformation, but must be shown quantitatively)
2.3 CMB fit: cold + non-oscillating + right Omega + right scaling -> must reproduce the acoustic peaks. HIGH BAR.

## Phase 3 -- Consistency & closure
3.1 Continuity: the recombination incoherent stress must smoothly COHERE into the late-time director strain
    (z~few ordering). One fabric, two guises -- check the transition conserves the right budget.
3.2 No double-counting: the fabric stress used as CDM must NOT be the same energy already counted as the CC/DE
    (Lambda_obs). Separate the coherent-vacuum (DE) from incoherent-stress (CDM) bookkeeping cleanly.
3.3 VERDICT: option 4 works (fabric = its own CDM, unifying recomb-CDM + galactic-DM as one fabric) OR option 3
    (gap; Omega_DM stays an input).

## Critical path & risks
CRITICAL PATH: Phase 0 (F1/F2/F3) -> Phase 1.3 (director-stress EOS) -> Gate 1. Phases 2-3 only if green.
RISK A (high): the answer may be DE-like (const) -- the condensate is a vacuum, and its bulk energy is the CC;
  the incoherent stress may just feed Lambda, not cluster as CDM. Then option 4 fails.
RISK B: clustering (2.2) may not reproduce the CMB peaks even if the scaling is matter-like (the HIGH BAR).
RISK C: the fabric picture may conflict with the framework's own FRW derivation (Gate 0) -> reframe needed.
COST: Phase 0 ~ read + reason (cheap, but conceptually load-bearing); Phase 1 ~ the real derivation; Phases 2-3
  ~ heavier (perturbation theory / CMB). Do Phase 0 + Gate 1 FIRST; they decide whether the rest is worth it.

## PHASE 0 -- DONE. Verdict: F2 (fabric stretches), 2I preserved.
- F1 (creation) OUT: no reason for the constrained hot "soup" if fabric is created on demand.
- F3 (fixed fabric) OUT: the fabric IS the substrate; monistic geometry has no separate metric to expand.
- F2 (stretch) IN, and the framework FIXES it: rho_CMB = 10 Lam_cond^4 (Paper XII) with rho_CMB ~ 1/a^4
  => Lam_cond ~ T_CMB ~ 1/a, xi_cond = 1/Lam_cond ~ a, rho_cond ~ Lam_cond^4 ~ 1/a^4 (RADIATION-like).
- 2I under isotropic stretch: PRESERVED (uniform scaling keeps icosahedral symmetry); only the SCALE grows
  (R_0 ~ a). No "elongation"/breaking (that needs anisotropic stretch, not observed).
- KEY: the condensate scale IS the temperature scale, so the fabric bulk tracks RADIATION, not matter.

## GATE 1 -- DONE (gate1_eos.py). Verdict: RED for route 1.
a-scaling of every fabric component:
  thermal condensate rho_cond ~ Lam_cond^4       -> 1/a^4  RADIATION (it IS the CMB reservoir)
  frozen-attractor vacuum Lam_obs (CC)           -> const  DARK ENERGY
  baryonic Hopfions (rest-mass)                  -> 1/a^3  MATTER, but only Omega_b~0.05
  DISORDERED director stress (recomb candidate)  -> 1/a^4  RADIATION (u~K/xi^2~rho_cond; flucts relativistic)
  GALACTIC disclination strain (baryon-attached) -> 1/a^3  MATTER, but LATE (z~few), ABSENT at recombination
NO fabric sub-component is matter-like (1/a^3) at recombination. Disordered stress = radiation; the only
matter-like strain is late. => OPTION 4 route-1 ("fabric bulk/stress = CDM") FAILS.
Remaining paths for "no separate CDM": ROUTE 2 (modified pre-recombination gravity reproduces the CMB
peaks with radiation+DE+baryons, NO cold matter) -- radical, whole-program. Else OPTION 3 (Omega_DM input/gap).

## NET (this line of inquiry)
The fabric supplies radiation (thermal condensate), DE (frozen vacuum), baryonic matter (Hopfions), and a
LATE-time galactic DM (disclination strain, z~few). It does NOT supply a cold-matter component at recombination.
So the recombination Omega_DM~0.26 is either (route 2) unnecessary if the modified gravity fits the CMB, or
(option 3) a genuine input/gap. Route 1 is closed. Decision now: scope route 2, or accept the gap.
