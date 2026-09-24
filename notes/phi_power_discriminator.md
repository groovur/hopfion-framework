# Jordan->Einstein discriminator over the phi-power "invariants" (phi_power_discriminator.py)

METHOD: a phi-power that arises from the SCALAR-TENSOR normalisation (non-canonical kinetic
K=1/phi^8, or the conformal factor F) is a FRAME ARTIFACT -- it cancels when you compute the
physical (canonical / Einstein-frame) quantity. A phi-power of TOPOLOGICAL/WZW origin (Bogomolny
lambda, quantum dimension, WZW T-matrix, McKay order) is PHYSICAL and frame-free. (Established by
the m_xi case: phi^20 = phi^12(V'') x phi^8(1/K) cancelled to phi^0, m_xi~H0.)

## Results
PHYSICAL (topological/WZW, SURVIVE):
  - lambda = phi^6  (Bogomolny/BPS, Paper II) -- soliton property, no gravity frame.
  - alpha^-1 = 360/phi^2 (WZW/anyonic S-matrix) -- non-gravitational, matches 137.036.
  - omega_BD ~ phi^12 = lambda^2, and alpha_BD ~ phi^-12 = 1/lambda^2 -- the phi-content is the
    topological Bogomolny lambda^2. The physical fifth-force coupling ~ phi^-6 (=1/lambda) is real.

FRAME ARTIFACT (scalar-tensor normalisation, CANCEL/CHANGE in Einstein frame):
  - m_xi^2 ~ phi^20 = phi^12(V'') x phi^8(1/K)  -> Einstein frame phi^0 (m_xi~H0). CANCELS.
  - r_V ~ phi^-8/3 = (phi^12/phi^20)^{1/3}: MIXED. The phi^20(artifact) cancels; the phi^12(=lambda^2,
    physical) survives -> physical r_V ~ phi^4. So r_V HAS a topological scaling (phi^4 from lambda^2),
    just NOT the phi^-8/3 the corollary claimed (that was contaminated by the m_xi artifact).
  - m_g^2 ~ phi^-8 (graviton TT): tied to the NMC (S_eff/(1+beta rho)); Einstein ~ /Omega^2 -> ~phi^-15.
    Frame-dependent AND untested (only v_GW~c matters, holds for any tiny m_g).
  - F(rho_inf) ~ phi^7 (3phi^7 = NMC ratio): this IS the Jordan->Einstein conformal map, frame-defining.

COINCIDENCE (not a phi-power at all):
  - Lambda_UV = M_Pl/(4pi sqrt2) ~ M_Pl/phi^6: 4pi sqrt2=17.772 vs phi^6=17.944, 1.0% off.
    A numerical coincidence (within the framework's own 0.1-2% chance-rate caveat), not derived.

## RULE OF THUMB (the payoff)
A phi-power is PHYSICAL iff its ORIGIN is topological/WZW (Bogomolny lambda, quantum dimension,
T-matrix, McKay order), NOT the non-canonical kinetic norm 1/phi^8 or the conformal factor F.
=> The GRAVITY/scalar-tensor sector is where the frame artifacts live (m_xi, r_V exponent, m_g, F).
   The PARTICLE/topological sector (alpha^-1, lambda, quark masses, WZW) is frame-free and robust.
This is a concrete criterion for separating real phi-structure from normalisation artifacts, and it
is falsifiable per-quantity by the Einstein-frame computation.

## Paper status
DONE: m_xi (P7:eq:mxi clarified), r_V corollary (P7:rem:rV_frame + tier/conclusion softened),
and m_g^2 caveat added after the v_GW bound (phi^-8 is Jordan-frame/NMC, no observable depends on it).
Lambda_UV~phi^6 needs NO edit -- the paper ALREADY calls 4pi sqrt2 ~ phi^6 'a structural coincidence
of the framework' (l.582) and shows the 1% ratio (l.629). Discriminator CONFIRMED the paper's own framing.
The topological phi-powers (lambda, alpha^-1, omega_BD content) need NO change -- physical.
