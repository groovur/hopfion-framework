# Quark sector: condensate-Cherenkov, confinement kinematics, the T(2,n) tower

Session 2026-09-02/03. Thread from the "snake.png" wake artifact → magnitude-mode Cherenkov →
confinement dynamics → the torus-knot tower. Epistemic tags per CLAUDE.md §6.
Scripts: `papers/src_paper18/cherenkov_confinement.py`, `.../bounded_motion_cherenkov.py`.
Paper: added `P18:rem:cherenkov_confinement` (after `P18:rem:production_conditions`) + Open Problem O7.

## 1. Condensate-Cherenkov (magnitude/DE mode) — the derivable kinematics [E from c_s=c/φ]
The condensate's scalar (magnitude/density) mode has sound speed c_s = c/φ (P7:rem:cs_de, §5).
A source (parton) at v>c_s radiates into it (Frank–Tamm / Cherenkov). All from c_s=1/φ:
  - threshold  β = 1/φ ≈ 0.618  (γ ≈ 1.272; KE threshold 0.27 mc²).
  - UR radiated fraction → 1 − 1/φ² = 1/φ  (golden identity 1/φ² = 1 − 1/φ). EXACT.
  - cone  cosθ_c = 1/(φβ) → arccos(1/φ) = 51.8° at β→1 (Mach angle 38.2°).
  - a light-speed parton has magnitude-Mach M = 1/c_s = φ.
Absolute dE/dx magnitude = coupling·cutoff, NOT derived (cutoff-set, §10). Threshold/fraction/cone ARE.

## 2. Confinement crossover: v = c/φ separates elastic ↔ radiative [E boundary; magnitudes O]
  - v < c/φ : deformation outruns the parton → tube STRETCHES quasi-statically → energy stored as
    STRING TENSION → snaps back = CONFINED (elastic regime). [user's "more input → more stretch → snaps back"]
  - v > c/φ : Mach cone → energy RADIATED (Cherenkov) → not recovered = de-confining.
Cherenkov does NOT set the string-tension magnitude (that's the tube's elastic modulus × cross-section,
cutoff-set). It sets the KINEMATIC boundary of clean confinement + the UR fraction 1/φ.

## 3. Bounded-motion resolution — why a CONFINED quark doesn't bleed [E scaling; I photon; O magnitudes]
Tension: hadronic partons are ultra-rel (v~c > c/φ) → naively always radiate → but hadrons don't bleed.
Resolution (this is the merit-worthy core, now in P18:rem:cherenkov_confinement):
  (a) Continuous Frank–Tamm cone needs UNBOUNDED coherent translation. A confined quark oscillates
      (no net translation) → the continuous cone never forms. Primary point.
  (b) The wound tube WAVEGUIDES the ripples the quark emits. Cutoff E_c ~ π ħ c_s/R; confined-quark
      scale E_q ~ π ħ c/R → RATIO IS R-INDEPENDENT:  E_q/E_c = 1/c_s = φ.  [the clean result]
      → confined quark sits a fixed factor φ ABOVE the sound cutoff: MARGINAL.
      Soft internal modes (below cutoff) = evanescent → TRAPPED as standing waves = internal hadron
      spectrum, no loss ("directed by the surrounding wound medium" — user). Hardest modes (~φ up) leak
      → FINITE WIDTHS (hadrons long-lived but not perfectly stable — physically right).
  (c) ESCAPE attempt (net translation, v>c/φ) → straight-line Cherenkov → ESCAPE DAMPING, opposing
      de-confinement (user's "also damped" — correct, but only on escape).
So continuous bleed is REPLACED by trapping (a,b) + escape damping (c).

## 4. Photon channel separation (the "photon" thought) [I reasoned]
Real photons escape the hadron (radiative decays) while sound-channel ripples stay trapped.
Reason = DECOUPLING, not cutoff: the photon is colour-NEUTRAL → doesn't couple to the colour/director
tube → escapes. (A cutoff argument gives the OPPOSITE: E_c^EM = φ·E_c^sound, faster EM mode trapped MORE
easily — so escape must be by decoupling.) "The hadron glows (photons) but doesn't hiss (sound trapped)."

## 5. Temperature axis — the QGP picture [O, user's hypothesis, NOT derived]
Two-axis de-confinement: SPEED (v>c/φ radiates) AND TEMPERATURE (order/disorder).
High-T early universe → condensate disordered → no stable trefoil → gluon soup; cools → order + "room"
for stability → trefoil stabilises → confinement onset ("never recover until cool enough" — user).
Already partly in the paper: P18:rem:production_conditions (two-threshold: energy + density; QGP through T_c)
and Open Problem O6 (reproduce T_c~150–160 MeV). The relation of the v>c/φ radiative boundary to T_c is
UNESTABLISHED (new Open Problem O7). The trefoil-stability temperature is the QCD/charge scale, NOT the DE
Λ_cond scale — mapping needs new work, can't inherit.

## 6. The torus-knot tower T(2,n) + the Q_H=2→3 satellite [E, Paper XVIII already has this]
User asked: higher trefoil configs? YES — Paper XVIII P18:def:torus_tower, P18:prop:preimage_tower:
  Q_H=1 → T(2,1) unknot (ν) · Q_H=2 → T(2,2) HOPF LINK (charged lepton = the TWO TUBES) ·
  Q_H=3 → T(2,3) trefoil (baryon; 3 crossings at 120°/Z_3 = the three quark colours).
  n=5 → cinquefoil T(2,5) (next fused knot), then n=7…
Jones tower {1/φ, 1, φ, −1, …}, period 10, mirror |J^(n)|=|J^(10−n)|; ODD n = fused genuine knot
(irrational, golden) = "confined"; EVEN n = unfused link (rational). 
USER'S KEY INTUITION CONFIRMED: "T(2,3) as a perturbance of the Q_H=2 field, 120° = 3 quarks" IS the
paper's SATELLITE construction (P18:prop:satellite): thread the Q_H=1 unknot through the T(2,2) Hopf link
with (2,3)-cable framing → the two-component link FUSES into the single-component trefoil. Panel (B):
T(2,2) and the T(2,3) satellite carry the SAME modulus 1/φ, differing by ONE unit of phase q_5 → Q_H=3 is
Q_H=2 wound up by one half-twist. This is the §8 "two pre-images reconfiguring into T(2,3)" made precise.
OPEN [O]: physical assignment of the higher tower — cinquefoil T(2,5) (5 crossings) pentaquark-shaped?
even-n unfused links = meson/multi-hadron? fused/unfused = single vs multi-hadron? period-10 = Q_group=10?
Paper XVIII computes the Jones invariants but does NOT assign the higher members — genuinely open.

## 7. Q_H=3 stability language (user's "unstable/decays" — confirmed) [E]
Q_H=3 is METASTABLE/confined, not stable: XV l.991 "not a stable isolated soliton"; XVIII l.948 metastable;
XIX l.72 "transient and unstable"; decay channel = "loses topological integrity at one/more crossing"
(XVIII l.2361, XIX l.444). Nuance: part of "no stable Q_H=3 saddle" was a gradient-flow TOOLING finding
(CLAUDE.md working-principles §4), not pure physics. Physical: topologically confined → long-lived, but
NOT the ground state (Q_H=2 is).

## THROUGH-LINE: one sound speed c_s=c/φ gives the whole confinement kinematics —
the deconfinement boundary (v>c/φ), the confined quark's Mach-φ, the marginal-widths factor φ, the UR
radiated fraction 1/φ, and (via colour-decoupling) the photon escape. Derived: kinematics. Open: coupling
magnitudes, photon-decoupling rigour, T_c link, tower particle-assignment.
