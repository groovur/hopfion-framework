# Semi-Dirac condensate polarization: small-q setup and first result

Model: 2D semi-Dirac two-band H(k)=(k_x^2/2m)sigma_x+(v k_y)sigma_y (m=v=1),
E=|d|, d=(k_x^2/2, k_y). Static undoped interband Lindhard polarization
  Pi(q)=int d^2k/(2pi)^2 (1 - dhat_k . dhat_{k+q})/(E_k+E_{k+q}),  Pi(0)=0.

## Result: anisotropic small-q exponents (confirmed)

| q-direction | fitted Pi(q) ~ q^s | prediction (anisotropic scaling) |
|---|---|---|
| quadratic/heavy axis (q_x) | s = 0.96 | 1 (linear, Dirac-like) |
| linear/light axis (q_y)    | s = 0.50 | 1/2 (square-root)      |

(anisotropic scaling k_x~L^{1/2}, k_y~L, E~L, d^2k~L^{3/2} => Pi~L^{1/2},
giving Pi~q_x^1 and Pi~q_y^{1/2}. Confirmed numerically.)

## Interpretation

- Pi(0)=0: vanishing DOS -> NO Thomas-Fermi screening -> long-range
  interaction is modified, not exponentially cut off.
- Anisotropic: the sqrt(q) light-axis channel is more IR-singular than the
  linear heavy-axis channel, so the anomalous long-range response is
  concentrated PERPENDICULAR to the radial (heavy) axis -- i.e. in the
  "less-radial" regions. This is the director / disclination channel.
- The sqrt(q) soft mode is the director Goldstone; its being unscreened is
  why the Frank elastic strain (dark matter) is long-range.

## Next pieces (open)

1. Frank constant K: from the DIRECTOR (transverse/orientational) correlator,
   not the density Pi computed here -- K is the q^2 coefficient of the
   director susceptibility. Density Pi confirms the soft channel exists;
   K quantifies its stiffness.
2. 3D embedding + actual V(r): map the 2D semi-Dirac response into the 3D
   galactic geometry (anisotropy relative to the local radial), get V(r),
   test M_enc(r) ∝ r (flat curves). Real-space director/disclination
   argument already gives this cleanly (theta~ln r -> M∝r); the polarization
   is the microscopic justification.
3. a_0 scale: from the semi-Dirac crossover momentum (where quadratic and
   linear terms balance), mapped to the condensate scale m_xi ~ H_0.
   Target: a_0 = cH_0/(2pi) (verified as the MOND scale).

Run: python polarization.py

# Frank K — correct setup for the (compute-heavy) run

Three quick attempts (frank_k.py, k2, k3) FAILED the Goldstone gate
(dE(q->0)/theta^2 must -> 0; got ~ -5 to -14). Diagnosis: the naive
d(R_theta k) expansion at fixed k omits the INTRABAND paramagnetic term,
which is what enforces the Goldstone protection. Correct formula (per k,
lower band filled; dtheta H = -G.sigma with G=(kx*ky,-kx); d2theta H =
G'.sigma with G'=(ky^2-kx^2,-ky); |-(k)>=(1,-e^{i phi})/sqrt2):

  dE/theta0^2 = INT d2k/(2pi)^2 [
      (1/4)(-dhat.G')                                              # diamagnetic
    + (1/4) sum_{s=+-} |<+(k+sq)|G.sigma|-(k)>|^2 / (-|d_k|-|d_{k+sq}|)   # interband
    + (1/4) sum_{s=+-} |<-(k+sq)|G.sigma|-(k)>|^2 / (-|d_k|+|d_{k+sq}|) ] # intraband
  K = 2 (dE/theta0^2) / q^2 ,  along q||y (bend) and q||x (splay) separately.

WHY IT IS EXPENSIVE (the setup requirements):
1. The intraband denominator (-|d_k|+|d_{k+sq}|) ~ q.grad|d| -> 0 near the
   semi-Dirac node and along |d_{k+q}|=|d_k| contours. K is the SMALL
   RESIDUAL of large, near-cancelling +q / -q intraband contributions.
   A uniform grid + eta-regularisation does NOT resolve it (frank_k5 gives
   a spurious ~ -700). REQUIRES: adaptive/node-resolving integration
   (dense sampling where |grad|d|| is small), proper principal value
   (symmetric +-q pairing to cancel the 1/den pole analytically, keeping
   the O(q^2) residual), and an eta -> 0 extrapolation.
2. Convergence in the k-cutoff Kmax AND grid N must both be checked; the
   quadratic band (kx^2) makes the large-k tail slow.
3. NON-NEGOTIABLE VALIDATION GATE: dE(q->0)/theta^2 must extrapolate to 0
   (Goldstone). Do NOT trust any K unless the gate passes to a few %.
   Report the gate residual alongside K.
4. Deliverables: K_splay, K_bend (sign MUST be >0 for a stable director;
   if <0, the uniform director is unstable -> modulated ground state, a
   physical result to report honestly), their ratio (the anisotropy),
   and the q-exponent (2 = standard Frank; anomalous exponent would mean
   the disclination flat-curve argument needs redoing).
5. Scale: map the crossover momentum (where kx^2/2 ~ ky) to the condensate
   scale m_xi ~ H_0 to get a_0; target a_0 = c H_0/(2pi).

STATUS: formula derived and validated in structure; numerics NOT converged
(the intraband PV needs the careful setup above). No reliable K yet.

# Frank K — RESULT (gate-validated, frank_run.py)

FIX that unblocked everything: use a RADIAL DISK cutoff |k|<Lam (a square/
product grid breaks the axis-rotation symmetry -> spurious Goldstone
violation growing with cutoff, which was the -5/-63/-700 in the failed
attempts). With the disk cutoff the q=0 Goldstone gate passes to ~1e-8 and
is Lam-stable. Undoped semimetal => intraband Pauli-blocked => interband +
diamagnetic only (no near-singular PV; the run is FAST, ~0.3 s/point).

Result (K = 2 dE/theta^2 / q^2), dE ~ q^2 (STANDARD Frank) for both:
  bend  (q||y, gradient along light axis):  K_bend  > 0  (STABLE)
  splay (q||x, gradient along heavy/radial): K_splay < 0  (UNSTABLE)
Signs are cutoff-robust (Lam=15,30,60 all agree). Magnitudes GROW with Lam
(K_bend 0.36->0.94->2.5; K_splay -3.6->-26->-166): the stiffness is
UV-cutoff-dependent (gapless quadratic band contributes at all k), so its
MAGNITUDE is set by the physical condensate cutoff, not universal.

PHYSICS: negative splay stiffness => the uniform director is unstable and
SPONTANEOUSLY SELF-ORGANISES into a radial splay texture. Splay = axis
rotation with gradient along the radial (heavy) direction => radial texture
=> disclination/spoke strain => the invisible dark-matter strain. Standard
Frank (q^2) => disclination theta~ln r => M_enc ∝ r => flat rotation curves.
This is the self-organisation mechanism; the scalar-magnitude relaxation was
blind to it (it lives in the orientation field).

Robust: sign of K (splay-unstable), q^2 exponent, Goldstone gate.
Cutoff-set: magnitude of K (=> a_0 scale enters via the physical condensate
UV cutoff ~ m_xi ~ H_0; target a_0 = cH_0/(2pi)).

Run infra (frank_run.py): resumable (checkpoint.jsonl, partial-sum records,
skip-on-rerun), test/full modes. Full mode (fine grid, small-q, silent) is
available for precision but NOT needed for the sign/mechanism, which the
fast test already establishes. Next physics step: physical-cutoff
regularisation + 3D embedding to turn the magnitude into a_0.

# Director-wave speed, susceptibility, ordering scale (director_speed.py)

Goldstone relation for the broken-rotation director: c_dir^2 = K/chi, with
  K   = Frank stiffness (frank_run: K_bend>0 stable, K_splay<0 unstable),
  chi = static susceptibility of the rotation generator G,
      = INT_disk d2k/(2pi)^2 |M(k,k)|^2/|d_k|  (same M as the diamagnetic term).

Result over Lam = 15,30,60,120 (disk cutoff, band units v=1):
  chi     ~ Lam^3.51   (matches anisotropic-scaling dim. analysis L^3.5 exactly)
  K_bend  ~ Lam^1.42     K_splay ~ Lam^2.70
  => c_dir^2 (bend)  ~ Lam^-2.09   (splay) ~ Lam^-0.82   -- BOTH decreasing.

  cdir2_bend: 8.2e-4, 1.9e-4, 4.4e-5, 1.1e-5   (-> 0 in the UV)
  |cdir2|_splay: 8.3e-3, 5.3e-3, 2.9e-3, 1.5e-3

KEY PHYSICS (robust to normalization; only signs + Lam-scaling are used):
- The DIRECTOR (orientation) mode is COLD: c_dir^2 is small and DECREASES with
  cutoff (chi grows faster than K -- the director is "heavy", large rotational
  inertia). It is NOT relativistic. The old c_s = 1/phi ~ 0.6c that killed the
  fluid picture was the scalar-MAGNITUDE mode; the director does not share it.
  At Lam=120, c_dir ~ 3e-3 c and falling -> near/below the CDM bound c_s<~1e-3 c;
  at the (larger) physical condensate cutoff it is comfortably cold.
- => the director-strain DM is genuinely cold/pressureless: no Jeans obstruction,
  clusters on all scales -- recovering exactly the CDM clustering the fluid failed
  to give (there Jeans length > horizon). Static strain gravitates; the slow mode
  provides no pressure support (static/dynamic decoupling confirmed).
- splay c_dir^2 < 0 = imaginary frequency = the self-organization instability;
  |c_dir^2|_splay ~ Lam^-0.82 sets the growth-rate scale.

CUTOFF-SET (same open problem as a_0): absolute c_dir, the ordering temperature
T_ord ~ (pi/2)K_bend (band units, ~Lam^1.4), and the texture/isocurvature
fraction all scale with the physical condensate cutoff ~ m_xi. Robust = coldness
(direction), splay instability (sign); cutoff-set = magnitudes.

Run: python director_speed.py

# Physical-cutoff regularisation + 3D embedding (cutoff_embedding.py) -- SETUP

Smooth rotationally-symmetric form factor f(|k|)=exp(-|k|^2/Lam^2) replaces the
sharp disk cutoff (keeps the Goldstone gate: |k|-only). 3D embedding: n_hat=e_r,
|grad e_r|^2=2/r^2, u=(1/2)|K|(2/r^2), M_enc=INT u 4pi r^2 dr ∝ r (flat curves),
a_0 = v_flat^2/r_*, r_* = c/m_xi ~ Hubble radius.

RIGOROUS (confirmed under smooth cutoff):
- Goldstone gate = 0.0 exactly (form factor preserves it).
- Coldness robust to regularisation scheme: cdir2_bend 2.4e-3->5.3e-4->1.2e-4
  (Lam=10,20,40), splay cdir2<0. NOT a sharp-cutoff artifact.

a_0 EMBEDDING (target a_0/cH_0 = 1/2pi = 0.159):
  Lam=10: band_K=0.119 -> a0=8.1e-11 m/s^2   (obs 1.2e-10)
  Lam=20: band_K=0.302 -> a0=2.0e-10
  Lam=40: band_K=0.790 -> a0=5.4e-10
  a_0 comes out the RIGHT ORDER (~10^-10 m/s^2) with NO tuning, and the target
  cH_0/2pi corresponds to band_K=0.159, i.e. a physical cutoff Lam*~11-12 band
  units. So the open problem is REDUCED to one number: does m_xi sit at Lam*~11
  in band units? -> the crossover-matching step (match the band heavy/light
  crossover to the real condensate dispersion). This is the remaining physical
  input; the a_0 coefficient is otherwise cutoff-set (band_K grows with Lam).

ISOCURVATURE scaffold: T_ord ~ (pi/2)K_bend (band units, cutoff-set); rho_tex ~
|K|/xi^2, xi~1/Lam; f_iso = rho_tex/rho_DM needs the rho_DM(m_xi) normalisation
(same unit map as a_0). Planck bound beta_iso<~0.038 to check once normalised.

NEXT: the crossover-matching (fix Lam* = m_xi in band units) closes BOTH a_0 and
the isocurvature fraction -- they share the one cutoff.

Run: python cutoff_embedding.py

# Crossover matching -- RESULT (crossover_match.py)

Fix Lam* = m_xi in the band's crossover units (k*=2, E*=2.83); then
a_0/cH_0 = band_K(Lam*), band_K(Lam) = 5.23e-3 * Lam^1.36 (smooth cutoff).

- TARGET a_0 = cH_0/(2pi) (band_K=0.159) => Lam* = 12.4 crossover units,
  i.e. cutoff/crossover hierarchy m_xi/E* ~ 12 (light-axis map) to ~38 (heavy).
- SINGLE-SCALE (m_xi ~ E*, no hierarchy, Lam_band~k*=2): band_K=0.015 =>
  a_0 = 1.0e-11 m/s^2 = 10.6x BELOW observed (1.2e-10). a_0 ~ cH_0 order holds.

VERDICT (honest): a_0 ~ cH_0 is robust (right order, no tuning). The exact
cH_0/(2pi) is NOT derived -- it requires the condensate UV cutoff to sit
~10-40x above the semi-Dirac crossover. This is a specific, falsifiable
condition, not a free fit.

SPECULATIVE lead (flagged; framework-native-constants rule applies): the
required hierarchy window (12-38) contains phi^6 = 17.94 -- and phi^6 is NOT a
random number here, it is the S_eff normalisation in the dispersion itself
(S_eff = sin^4 theta/[phi^6(1+beta rho)]). Taking m_xi/E* = phi^6:
  light-axis map:  Lam*=17.9 -> a_0 = 1.79e-10 = 1.49x observed
  heavy-axis map:  Lam*= 8.5 -> a_0 = 6.5e-11  = 0.54x observed
So the framework's OWN phi^6 brackets the empirical a_0 within a factor ~1.5
either way -- suggestive (the hierarchy is physically motivated, not fitted),
but it does NOT nail the 2pi. Status: promising, unproven.

ISOCURVATURE: at Lam* the ordering scale T_ord ~ (pi/2)K_bend ~ 0.25 m_xi, i.e.
the director orders LATE (near the condensate scale, not at reheating). A late,
cold texture is naturally subdominant -> expected to satisfy Planck beta_iso <
0.038, but a number still needs the rho_DM(m_xi) normalisation.

Run: python crossover_match.py

# Isocurvature number (isocurvature.py) -- RESULT

The splay texture is DM only once the ultralight director (m_xi ~ H_0, Eq mxi)
UNFREEZES: H(z_order) = m_xi.  For m_xi/H_0 = 1, z_order = 0; for 2, z_order=1.2.
=> the texture forms at z ~ 0 (present epoch), essentially absent at z=1100.

Frozen inflationary director fluctuation (would-be isocurvature seed):
  H_inf/M_Pl = 1.0e-5 (P_s=2.1e-9, Starobinsky eps=2/N_e^2=6e-4),
  f_dir ~ Lam_UV ~ M_Pl/phi^6, delta_theta ~ H_inf/(2pi f_dir) = 2.85e-5,
  comparable to the adiabatic sqrt(P_s)=4.58e-5.

beta_iso(CMB) = f_tex(z=1100) x (delta_theta/sqrt(P_s))^2 ~ 0, since f_tex(1100)~0.
=> SAFELY below Planck beta_iso < 0.038.

HONEST VERDICT (the two sides of one fact):
- SAFE: beta_iso(CMB) ~ 0 because the texture is late (m_xi~H_0 => z_order~0).
- COST: the same lateness means the mechanism supplies NO DM at recombination
  -> it CANNOT be the CMB-era cold dark matter (peak heights, early growth).
  The director-strain DM is a LATE-TIME, galactic-scale component; early-universe
  CDM is a separate, still-open account.
- The frozen delta_theta ~ 3e-5 is itself ~ sqrt(P_s): harmless ONLY because the
  texture does not gravitate until z~0. Any texture DM at recombination would give
  beta_iso ~ O(1) (excluded). So the lateness is REQUIRED, not incidental.

Run: python isocurvature.py

# S_eff <-> semi-Dirac: DERIVED (dispersion_from_Seff.py, sympy-verified)

NOT previously done (the band model H=(kx^2/2)sx+ky sy was ASSUMED in all DM
scripts). Now derived symbolically:

Semi-Dirac with heavy axis PERPENDICULAR to e_r, k_par=k cos th (radial),
k_perp=k sin th:  H = A k_perp^2 sx + v k_par sy,
  E^2 = A^2 k^4 sin^4(th) + v^2 k^2 cos^2(th).
- heavy (k^4) angular factor = sin^4(th) == S_eff  [sympy: TRUE, exact].
- intrinsic: d_x = A k^2 sin^2(th) (transverse quadratic band); E^2 squares it
  -> sin^4 = (sin^2)^2. No external field. => user is RIGHT: sin^4 is the band
  structure, not a coincidence or analogy.
- limits: th=0 (radial) E^2=v^2 k^2 (linear/light/EASY); th=pi/2 (perp)
  E^2=A^2 k^4 (quadratic/heavy/HARD). Matches "easy radial, hard perpendicular".

DIRECTION FIX (the catch): the match REQUIRES heavy = PERPENDICULAR to e_r,
light = radial. The paper's DM intro says "quadratic (heavy) ALONG e_r" and the
note's dispersion E^2=A^2 k^4 cos^4+v^2 k^2 sin^2 puts heavy at radial -- BOTH
BACKWARDS (cos^4 != sin^4; sympy confirms paper form's heavy factor = cos^4, not
S_eff). Must read heavy=perpendicular, light=radial.

CONSEQUENCE TO CHECK: the scripts use the same backwards convention (results.md:
"splay = q||x = heavy/radial"). Frank convention: splay ~ q PERP to n_hat, bend
~ q ALONG n_hat. With n_hat=e_r and heavy=perp now fixed, the heavy/radial
assignment AND the splay/bend<->q-direction labels need a careful geometric
re-pass -- the UNSTABLE mode may be bend, not splay. The qualitative result
(one mode unstable -> self-organised texture -> q^2 -> flat curves) is
label-independent and survives; only WHICH mode is unstable is in question.

Run: python dispersion_from_Seff.py

# Splay/bend labelling — RESOLVED (splay_bend_check.py, sympy)

Director n_hat = e_r = LIGHT axis (from the S_eff derivation). Rotate by
theta(x,y); compute div n_hat (splay) and |n x curl n| (bend) explicitly:
  - modulation along HEAVY axis (x, PERP to e_r): splay != 0, bend = 0  => SPLAY
  - modulation along LIGHT axis (y, ALONG e_r):   splay = 0, bend != 0  => BEND
(consistent with the Frank convention: splay = grad PERP to n_hat, bend = grad
ALONG n_hat.)

frank_run: K(heavy-axis modulation) = -166 < 0 ; K(light-axis modulation) = +2.5 > 0.
=> SPLAY is the UNSTABLE mode; BEND is stable. The paper's "splay instability" is
CORRECT. My earlier worry that it might be bend was WRONG (verified, not assumed).

NET of the semi-Dirac check: the identification is DERIVED and the splay
instability stands. The ONLY correction to the paper is cosmetic-but-real DIRECTION
wording: "quadratic (heavy) along e_r" -> "heavy PERPENDICULAR to e_r, linear
(light) along e_r"; dispersion form A^2 k^4 cos^4 + v^2 k^2 sin^2 -> A^2 k^4 sin^4
+ v^2 k^2 cos^2. K_splay<0, disclination, M_enc ∝ r, flat curves all UNCHANGED.

Run: python splay_bend_check.py
