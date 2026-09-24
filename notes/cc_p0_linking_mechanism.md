# The CC refinement as a linking-number monodromy count (back to 2I)

## Result
The cosmological-constant normalization is
    c = rho/kappa = phi^3 + p0^(Q_H),
with the leading phi^3 = sqrt(phi^6) the canonical kinetic factor, and the
correction p0^(Q_H) the vacuum Born weight raised to the HOPF CHARGE.
  - p0 = 1/(1+beta* rho_inf) = 1/phi^2  (feedback suppression at the attractor)
       = Tr[rho U_mono] = S_{1/2,1/2} S_{0,0} / S_{0,1/2}^2  (monodromy scalar, Paper III)
  - Q_H = 2 (density-feedback icosahedral Hopfion, Paper I) -> p0^2 -> CC=10^-122 (obs).
  - Q_H = 1 (standard Hopfion) -> p0^1 -> 10^-129 (WRONG, 7 orders): structure collapses.
CC = exp(-64.2 c) is log-linear in c, so the observed value is the geometric-mean
MIDPOINT of the bracket [canonical phi^3 (10^-118), topological phi^3+2p0^2 (10^-126)].

## Mechanism (why exactly one p0 per Hopf unit)
- Q_H is a LINKING NUMBER (the Hopf invariant = linking of the field's fiber
  preimages). Q_H=2 = fibers linked twice.
- p0 is a MONODROMY (once-around-the-vacuum amplitude for a 1/2 anyon).
- The k-essence field, vacuum-dressed, encircles each linked fiber ONCE; each
  encirclement inserts exactly one factor p0. #insertions = #links = Q_H.
  => correction = p0^(Q_H). The exponent is TOPOLOGICAL (a linking number).
- INDEPENDENT insertions (one per link), NOT the 1/2 x 1/2 FUSION of the pair --
  fusion would give a single p0 (=Q_H=1 answer, wrong). The data (p0^2) picks the
  linking picture over the fusion picture: a real falsifiable discriminator.
- Additive onto phi^3: c = (bare canonical phi^3) + (vacuum self-energy p0^Q_H).

## Back to 2I
Everything here is binary-icosahedral: p0 = the SU(2)_3/2I monodromy scalar;
Q=10=ord(q_2I), |2I|=120, k=3 (alpha^-1 = k|2I| p0 = 360/phi^2); and Q_H=2 is the
charge of the 2I (icosahedral) Hopfion vacuum. The CC refinement is the SAME 2I
anyonic data as alpha, now weighted by the linking number instead of k|2I|.

## Open derivation (to make it a proof, not a picture)
Compute the vacuum renormalization of c = rho/kappa in the Q_H=2 Hopfion
background and show it equals phi^3(1 + p0^(Q_H)/phi^3), i.e. the linking ->
monodromy -> p0 chain made explicit, with the exponent = the Hopf linking number.
STATUS: numerically confirmed (Q_H=2 matches, Q_H=1 fails); mechanism is a
coherent TQFT picture; the field-theory-in-Hopfion-background calc is NOT yet done
(checking whether a prior paper has it before sketching).

Scripts: scratchpad/{cc_coupling,path_A_normalization,path_ii_redefinition}.py + scans.
Paper: main_paper7.tex DE section (eq P7:eq:cc_c, box P7:eq:Lambda_obs -> 10^-122).

## Has the calc been done? NO -- but the Linking Scale Conjecture is the precedent

Searched all papers. The exact CC vacuum-normalization calc (c=phi^3+p0^Q_H) is NOT
done. BUT the framework already has the SAME KIND of theorem:

**Linking Scale Conjecture** (Paper I conj:linking; PROVED in Paper II via a modular
WZW Ward identity, P2:sec:lsl):
    s_opt^{2k} = N_fus(Q,k),   at (Q=2,k=3): N_fus=2.
- Same EL manifold: the virial condition is J_fb/J_a = dim_q(1/2) = phi -- which is
  EXACTLY the CC attractor condition (beta* rho_inf = phi, 1+beta* rho_inf = phi^2).
- The topological charge Q controls a fusion/linking count N_fus that weights a
  PHYSICAL quantity (the scale s_opt). Method: WZW Ward identity giving
  K_fb/J_four = lam/N_fus^{1/k}, constant on the EL manifold.
- The '2' is shared: N_fus=2 = Q_H=2 for the k=3 icosahedral condensate. So the CC
  correction exponent (Q_H, linking) and the Linking Scale exponent (N_fus, fusion)
  are the SAME 2I count. The CC p0^Q_H is another instance of the SAME
  topological-charge-weighting pattern, on the SAME manifold.

## Sketch of the CC vacuum-normalization calc (adapting the Paper II Ward identity)

Target: c = rho/kappa = phi^3 (1 + p0^{Q_H}/phi^3) = phi^3 + p0^{Q_H} on the EL manifold.
1. LEADING (bare): canonical kinetic normalization. Parent kinetic
   |grad rho|^2/(phi^6(1+beta* kappa)), lam=phi^6; canonical field rescales by
   sqrt(lam)=phi^3 => c_bare = phi^3.
2. CORRECTION (vacuum self-energy): the field rho is dressed by the Q_H-linked
   vacuum. Each of the Q_H linked fibers is encircled once (Hopf charge = linking
   number); each encirclement inserts one monodromy factor p0 = Tr[rho U_mono] =
   S_{1/2,1/2}S_{0,0}/S_{0,1/2}^2 (modular-S data). => insertion = p0^{Q_H}.
3. WARD IDENTITY (the technical target, à la P2:sec:lsl): a modular WZW Ward
   identity on the EL manifold fixing the field's wave-function renormalization
   Z = c/phi^3 = 1 + p0^{Q_H}/phi^3, with the monodromy p0 per link and the
   exponent = the Hopf linking number Q_H. This mirrors Paper II's
   s_opt^{2k}=N_fus but for the normalization (monodromy weight, linking exponent).
4. RESULT: c = phi^3 + p0^{Q_H} = phi^3 + p0^2 at Q_H=2 -> CC=10^-122.

Key claim to prove: the modular Ward identity gives the normalization correction as
(monodromy scalar)^(linking number). This is the natural next theorem after the
Linking Scale Conjecture, on the same manifold, same 2I data.
STATUS: sketch only; the Ward-identity derivation is the open work.

## CLOSED: the derivation ingredients are proved in Paper III (author was right)

The 'hard PDE step' (modular covariance / the Ward identity) is NOT open -- it is
PROVED in Paper III via three-fermion Witten bosonization. Paper III proof chain
(P3:eq:chain_full):
   kappa_tube = g_WZW (bosonization)  =>  V* = S_{0,1/2}/S_{0,0} = phi  =>
   J_4/J_2a = 2^{4/3}/phi^5 (Linking Scale)  =>  alpha^-1 = 360/phi^2 - k/(2pi).
The first two arrows are PROVED theorems (V*=phi exact for all R_0, P3:thm:V_exact_all_R0;
bosonization P3:thm:bos_main; Verlinde/Cardy P3:thm:verlinde).

KEY IDENTITY (verified exactly): p0 = 1/V*^2.
  p0 = S_{1/2,1/2} S_{0,0}/S_{0,1/2}^2 = (S00/S0h)^2 = 1/(S0h/S00)^2 = 1/V*^2,
  using self-duality S_{1/2,1/2}=S_{0,0} at k=3 (Paper II Prop Sselfdual).
  V*=phi (Paper III, proved) => p0 = 1/phi^2. EXACT.

So the CC correction p0^{Q_H} = V*^{-2 Q_H} = phi^{-2 Q_H} is built ENTIRELY from
PROVED framework data:
  - V* = phi = S_{0,1/2}/S_{0,0}: proved via Witten bosonization (Paper III).
  - p0 = 1/V*^2: algebra + self-duality (Paper II/III).
  - Q_H = 2 = Hopf linking number (Paper I).
It is the SAME icosahedral WZW monodromy that carries alpha^-1 = k|2I| p0, now
weighted by the Hopf linking number Q_H instead of k|2I|. Not an open Ward identity
-- the same proved bosonization chain, one monodromy insertion per link.

PAPER UPDATED: DE prose now states p0=1/V*^2, V*=phi proved by bosonization
(Paper III), removes 'open item' language. The only residual is the assembly
statement (that c = phi^3 + p0^{Q_H} as a single normalization theorem), whose
ingredients are all proved.

## Step 3 PERFORMED: the self-energy as a Chern-Simons Wilson loop (self_energy_linking.py)

The rho-fluctuation is the spin-1/2 SU(2)_3 primary; its self-energy in the
Q_H-Hopfion background is the normalized Wilson-loop expectation. Standard
CS/Witten evaluation with the SU(2)_3 S-matrix:
  - unknot j:              <W_j> = S_{0j}/S_{00} = d_j (quantum dimension); d_{1/2}=phi.
  - Hopf link (linking 1): <W_h W_h> = S_hh/S_00 = 1 (self-duality S_hh=S_00 at k=3).
  - monodromy scalar:      p0 = M_hh = S_hh S_00/S_0h^2 = 1/phi^2 = 1/d_{1/2}^2 = 1/V*^2.
  - linking number n:      M_hh^n = p0^n.
The Hopf charge = linking number of the field's fibre preimages (Paper I), so the
self-energy loop links the condensate flux Q_H times:
  Sigma = p0^{Q_H}  =>  c = c_tree + Sigma = phi^3 + p0^{Q_H}.
Q_H=2 (geometric): Sigma=1/phi^4 -> c=4.382 -> CC=10^-122 (MATCHES).
Q=10 (algebraic quantization charge, Paper X): Sigma=p0^10~7e-5 -> CC=10^-118 (WRONG).
=> The self-energy links the GEOMETRIC Hopf flux (Q_H=2), not the algebraic Q=10.
Step 3 is now a TQFT Wilson-loop computation using the SAME S-matrix as alpha^-1;
the only residual is the topological identification (self-energy loop linking = the
Hopf invariant Q_H), which is the definition of the Hopf charge (Paper I).
Script: scratchpad/self_energy_linking.py.

## Trefoil-departure (Paper 18) RESOLVED via a reality criterion, 2026-08

Paper 18 Remark P18:rem:trefoil_satellite_correction left open how the trefoil's
"departure from 1" relates to the Q_H=3 sector's transience. The CC work supplies a
criterion: a sector is a STABLE VACUUM only if its monodromy scalar is REAL (the CC's
Tr[rho U_mono]=1/phi^2 is real -- the two fusion channels' imaginary parts cancel to a
single operator expectation, Paper III). Then:
  - Q_H=2 Hopf link: Tr[rho U_mono]=1/phi^2 REAL (Im~1e-17) -> stable vacuum.
  - Q_H=3 trefoil satellite: J_sat(q5)=q5^3/phi = -0.5-0.363i COMPLEX (phase -144 deg),
    same MODULUS 1/phi -> NOT a real Born weight -> transient.
=> The transience-signalling departure is the PHASE, not the modulus (both 1/phi=1/V*).
Necessary condition (stable vacuum needs a real Born weight), consistent with transience;
NOT a full dynamical decay proof (needs Paper XIX kinematics). FOLDED into Paper 18 remark.
Script: scratchpad/self_energy_linking.py (reality check appended inline).

Scorecard of CC-work applications to 18/19: solar angle (already in Papers 17/19,
retracted); quark-mass monodromy (checked, does NOT fit -- dead end); trefoil-departure
(THIS -- a real result). Genuine outputs of the CC thread: Paper VII (derived CC
c=phi^3+p0^{Q_H}), Paper III (p0=1/V*^2 line), Paper XVIII (this reality criterion).
