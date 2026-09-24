# The m_xi phi^10 inconsistency (math run-down)

TRIGGER: WP4 needs m_xi ~ few H0, but P7:eq:mxi's formula gives m_xi = phi^10 sqrt3 H0 ~ 213 H0,
while the paper ALSO says "m_xi ~ H0 ~ 1e-33 eV" (l.677, 968, 1493). Which is right? -> the math.
Script: dm_polarization/mxi_math.py.

## The phi^10 is REAL (not a typo)
m_xi^2 = V''(rho_inf)/K(rho_inf), the standard non-canonical-scalar mass.
  V(rho)=Lam0 e^{-phi^6 rho} => V'' = phi^12 V ; at rho_inf, V''(rho_inf)=phi^12 Lam_obs.
  K(rho_inf)=1/[phi^6(1+beta rho_inf)] = 1/(phi^6 phi^2) = 1/phi^8   (beta rho_inf=phi, 1+phi=phi^2).
  m_xi^2 = V''/K = phi^12 * phi^8 * Lam_obs = phi^20 Lam_obs.
So phi^20 = phi^12 (potential curvature) x phi^8 (non-canonical kinetic norm). GENUINE.
=> m_xi = phi^10 sqrt(Lam_obs) M_Pl = phi^10 sqrt3 H0 = 213 H0 = 3.2e-31 eV.

## Three consequences
1. NUMERICAL ERROR in the paper: "m_xi ~ H0 ~ 1e-33 eV" (l.677,968,1493) is WRONG by phi^10 sqrt3
   = 213x. Correct order is ~1e-31 eV, not 1e-33 eV. ("~H0" is OK only as 'Hubble-SCALE/ultralight'.)
2. DIMENSIONAL slip in P7:eq:mxi as written: "m_xi^2 = phi^20 Lam_obs M_Pl^4" with Lam_obs an
   energy density (E^4) gives m_xi ~ E^4. The number-correct form is
   m_xi^2 = phi^20 (Lam_obs/M_Pl^4) M_Pl^2  [Lam_obs/M_Pl^4 = dimensionless 1e-122]. So "M_Pl^4" -> "M_Pl^2".
3. THAWING-DE TENSION (physics, real): m_xi=213 H0 => field thaws at H(z)=m_xi => z_thaw~52 (matter
   era) and has ROLLED since. But the proved thawing relation w_a=-3(1+w0) and w~-1 today need the
   field FROZEN until ~now, i.e. m_xi <~ H0. So m_xi=213 H0 tensions the DE sector's own thawing story.
   (Cassini/Vainshtein SURVIVES: r_V ~ m_xi^{-2/3} only shrinks by 213^{2/3}~35, still >> solar system.)

## Net
- The phi^10 is derived and correct; the paper's "m_xi ~ 1e-33 eV" throwaway is the error.
- With the correct m_xi~213 H0, WP4's DM Frank constant undershoots by ~3.6 dex (v=c/phi), NOT the
  +1 dex the loose "m_xi~H0" gave. So the elastic-DM quantitative match is ~3-4 dex short.
- AND thawing DE wants m_xi <~ H0, which the phi^10 formula violates. So BOTH the DM near-match and
  thawing DE would prefer m_xi ~ H0 (no phi^10) -- but the V''/K math gives phi^20. This is a genuine
  internal tension in Paper 7's DE/DM sector, for the author to adjudicate. NO paper edit made.
  Candidates to check: (a) is V really Lam0 e^{-phi^6 rho} at the minimum (or does something flatten
  V'' by phi^12)? (b) is the kinetic norm 1/phi^8 or should the DE mass use canonical M_Pl norm?

## Candidates (a) and (b) CHECKED (mxi_ab_check.py) -- both CONFIRM the enhancement
(a) V has NO minimum: V=Lam0 e^{-phi^6 rho} is a pure runaway (V'=-phi^6 V never 0). The
    "potential minimum" is the ATTRACTOR rho_inf (frozen by Hubble friction, not V'=0).
    V''(rho_inf)=phi^12 Lam_obs is correct; the bare V does not flatten it. (Only an
    Einstein-frame/nonminimal redefinition could, and that is not what the paper does.)
(b) The non-canonical norm K=1/phi^8 DOES belong: canonical chi has m^2=V''/(2K)=phi^20 Lam_obs/2.
    Even DROPPING phi^8 (pure canonical) gives m_xi=phi^6 sqrt3 H0 ~ 31 H0 -- still >> H0.
=> Neither candidate rescues m_xi ~ H0. m_xi is robustly 30-210 H0 ~ 1e-31 eV.

## Consequence for the paper edit (WHY I HELD OFF)
Because (a)/(b) confirm the enhancement, consequences 1 (the eV value) and 3 (thawing) are
COUPLED, not independent typos: if the paper is corrected to state m_xi ~ 200 H0, it would
STILL prove the thawing relation w_a=-3(1+w0) -- which needs m_xi <~ H0 -- so it becomes
INTERNALLY INCONSISTENT. Silently fixing the eV number makes the paper worse, not better.
This is a genuine DE-sector tension for the author to adjudicate, TWO coherent resolutions:
  A) phi^20 stands => m_xi ~ 200 H0 ~ 1e-31 eV. Then: fix all "m_xi~H0~1e-33 eV" to ~1e-31 eV
     AND qualify/revisit thawing (a field that thawed at z~50 is not w~-1-today thawing);
     WP4 DM undershoots ~3.6 dex; Vainshtein/Cassini SURVIVES (r_V only /35).
  B) m_xi ~ H0 is physically mandated (thawing + DM + operational usage) => the DE mass is
     NOT V''/K in the Jordan frame; it must be defined in the EINSTEIN frame, where the
     nonminimal coupling F(rho_inf) flattens V'' by ~phi^12 and restores m_xi~H0. Then fix
     eq:mxi's derivation, keep "m_xi~H0", thawing fine, WP4 near-match (+1 dex).
Resolution B is physically favoured (thawing is a proved theorem; DM near-match; operational
usage) and is a KNOWN subtlety (scalar-tensor masses are frame-dependent; the fifth-force /
cosmological mass is the Einstein-frame one). CHECK NEEDED: redo V''/K in the Einstein frame
with F(rho_inf)=M_Pl^2(39phi+25)/2 -- does the conformal factor remove the phi^12? If yes, B.
NO paper edit made pending this.

## RESOLVED: Einstein-frame mass (einstein_frame_mxi.py) -- m_xi ~ H0 is PHYSICALLY CORRECT
F(rho)=(M_Pl^2/2)(1+3phi^6 beta rho) (P7:eq:Feff, linear NMC; 3phi^7~87 at rho_inf -- STRONG).
Einstein frame g~=Omega^2 g, Omega^2=f/M_Pl^2, f=2F:
  U(rho)=V/(1+a rho)^2,  (dchi/drho)^2 = M_Pl^2[k/f + (3/2)(f'/f)^2],  a=3phi^6 beta.
At rho_inf the canonical norm is DOMINATED by the NMC mixing (3/2)(f'/f)^2 ~ (beta/phi)^2
(not k/f). The phi powers CANCEL symbolically:
  U'' ~ phi^12 U ~ phi^12 (Lam_obs/u0^2) ~ phi^12 Lam_obs/phi^14 = Lam_obs/phi^2
  (dchi/drho)^2 ~ M_Pl^2 (beta/phi)^2 = M_Pl^2 beta^2/phi^2
  m_E^2 = U''/(dchi/drho)^2 ~ Lam_obs/(M_Pl^2 beta^2) ~ H0^2/beta^2 ~ H0^2.   PHI-POWER-FREE.
Numerically: m_E = 1.07 H0 (physical, Einstein) vs m_J = 151 H0 (Jordan, eq:mxi phi^20).
The strong nonminimal coupling cancels the entire phi^20 Jordan enhancement.

## FINAL VERDICT -- the paper's physics is CORRECT; do NOT apply the "three fixes"
- Consequence 1 (m_xi ~ H0 ~ 1e-33 eV): CORRECT (Einstein/physical frame). NO error, NO fix.
- Consequence 3 (thawing DE): NO tension -- m_xi ~ H0 thaws now, consistent with w_a=-3(1+w0).
- WP4 DM near-match (+1 dex) HOLDS (m_xi ~ few H0).
- The ONLY residual: P7:eq:mxi "m_xi^2 = phi^20 Lam_obs M_Pl^4" is the JORDAN-frame mass and is
  quoted as "the DE scalar mass" -> misleading (and gives 151-213 H0, not H0). Optional CLARITY
  edit: note it is the Jordan-frame expression, whose strong-NMC conformal reduction gives the
  physical m_xi ~ H0. (Also the "M_Pl^4" is dimensionally the Jordan bookkeeping.) This is an
  author call, not an error in the physics. All "three fixes" WITHDRAWN.
Robustness: the phi-cancellation is symbolic (bookkeeping of phi powers), so m_xi ~ H0 is robust
to the O(1) corrections (the U' term at the non-extremal attractor shifts the O(1) coeff, not the
phi scaling). m_xi ~ (0.5-2) H0.
