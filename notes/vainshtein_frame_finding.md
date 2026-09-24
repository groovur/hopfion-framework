# Vainshtein r_V and the "phi^-8/3 topological invariant" corollary -- frame check

TRIGGER: the Einstein-frame m_xi finding (physical mass ~H0, phi-power-free; Jordan mass
=phi^20 ~213 H0). Does the Cassini result and the P7:cor:vainshtein_topology corollary survive?
Script: dm_polarization/vainshtein_frame.py.

## Two questions
Q1 (Cassini resolution): ROBUST. r_V propto m^{-2/3}; swapping Jordan m~213 H0 -> physical m~H0
   makes r_V ~36x larger (9.3e6 -> 3.3e8 AU), gamma_PPN 2.2e-5 -> 3.7e-6. BOTH << Cassini bound
   2.3e-5. r_V >> solar system either way. The Cassini/Vainshtein RESOLUTION does not depend on
   which mass is used. (P7:thm:cassini's conclusion stands.)

Q2 (the phi^-8/3 corollary): FRAME-FRAGILE.
   r_V^3 propto omega_BD/m^2 = phi^12/m^2.
   - Jordan m^2 = phi^20 Lam_obs -> r_V propto phi^{(12-20)/3} = phi^-8/3.  (paper: "8 = 2Q_cond-12")
   - Physical m^2 ~ Lam_obs (phi-free) -> r_V propto phi^{12/3} = phi^+4.
   The ENTIRE golden-ratio scaling of r_V flips (phi^-8/3 -> phi^+4) when the physical
   (canonical) mass replaces the Jordan one. The '8 = 2 Q_condensate - 12' bookkeeping -- hence
   "r_V is a topological invariant fixed by Q_condensate=10" -- rests on the phi^20 that the
   strong NMC CANCELS in the physical mass. So the corollary's phi-scaling is a JORDAN-FRAME
   ARTIFACT, not a physical invariant.

## The caveat that could rescue P7:cor:vainshtein_topology
IF the screening is CHAMELEON (density-dependent mass) and the LOCAL, solar-environment
(high-density) scalar mass -- not the cosmological ~H0 -- enters r_V, that local mass could
carry the phi^20 and the corollary could stand. The paper cites BOTH Vainshtein and chameleon
screening (l.122). WHICH mass legitimately enters r_V (cosmological Einstein ~H0, Jordan
non-canonical ~213 H0, or a local high-density chameleon mass) is the OPEN screening-mechanism
question that decides the corollary.

## Paper status
- P7:eq:mxi clarified (Jordan vs Einstein; my edit) and the Vainshtein sentence made frame-honest
  (uses the Jordan expression; Cassini robust either way). DONE.
- P7:cor:vainshtein_topology (r_V ~ phi^-8/3 "topological invariant") is NOT edited -- it needs
  the author to decide which mass enters r_V. If the cosmological/physical mass, the corollary's
  golden-ratio scaling does not survive; if a local chameleon mass carrying phi^20, it may.
  Flagged, no edit.

## CHAMELEON vs VAINSHTEIN resolved (chameleon_vs_vainshtein.py)
The two channels are DIFFERENT physics (paper rem:chameleon + thm:cassini):
- CHAMELEON (rem:chameleon): a density-dependent COUPLING, NOT a mass.
  omega_BD^eff = omega_BD(1+beta rho); at solar density beta rho_sun~1e6, omega_BD^eff~2.4e8
  >> Cassini 43000 (margin ~1e3-1e4). At high density this -> phi^12 rho/3, BETA- and
  m_xi-INDEPENDENT. This is the framework's PRIMARY, robust Cassini resolution (works at
  beta=0.452 directly, l.1050).
- VAINSHTEIN (thm:cassini): r_V=(G_N M omega_BD/m_xi^2)^{1/3}, uses the COSMOLOGICAL m_xi.
  Its phi^-8/3 corollary needs the Jordan m_xi^2=phi^20; the physical (Einstein) cosmological
  mass ~H0 gives r_V ~ phi^4. "Complementary" to the chameleon (l.1083).

DOES CHAMELEON RESCUE THE phi^-8/3 COROLLARY? NO. The chameleon carries NO phi^20 mass
(it is a coupling suppression). The only phi^20 is the Jordan cosmological m_xi in r_V, which
the strong NMC cancels in the physical frame. So there is no local chameleon mass to supply
phi^20 to r_V -- the earlier "chameleon could rescue it" caveat is CLOSED (it cannot).

## FINAL on #2
- Cassini resolution: ROBUST -- carried by the CHAMELEON (density-dependent coupling,
  m_xi- and frame-independent). P7:thm:cassini's CONCLUSION stands solidly.
- P7:cor:vainshtein_topology (r_V ~ phi^-8/3 "topological invariant fixed by Q_cond=10"):
  a JORDAN-FRAME ARTIFACT. Physical mass ~H0 -> r_V ~ phi^4. NOT rescued by the chameleon.
  Demoting/qualifying it does NOT threaten Cassini (the chameleon carries that independently).
  Recommend: add a frame caveat to the corollary (its golden-ratio scaling is a Jordan-frame
  quantity), or demote it from "topological invariant". Author call; NOT yet edited.

## PAPER EDITS DONE
- Added Remark P7:rem:rV_frame after the corollary proof: the phi^-8/3 exponent carries the
  Jordan-frame m_xi^2=phi^20; physical Einstein mass is phi-power-free (m_xi~H0), so r_V ~ phi^4
  (mass corrected) and further in a full Einstein treatment -> the golden-ratio scaling is a
  Jordan-frame quantity, "topological invariant" reading holds only in that frame; Cassini secured
  frame-robustly by the chameleon (rem:chameleon).
- Tier table (l.2436): "topologically determined / Proved / exact" -> "phi^-8/3 golden-ratio scaling
  / Proved (Cor), Jordan-frame (Rem) / Cassini secured by chameleon".
- Conclusion (l.2637): dropped "topological invariant... closing the open problem"; now "golden-ratio
  scaling phi^-8/3 ... Jordan-frame quantity, Cassini secured frame-robustly by the chameleon".
Env balance intact (remark 22/22, equation 89/89). #2 fully closed.
