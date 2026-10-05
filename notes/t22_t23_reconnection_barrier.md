# T(2,2)→T(2,3) reconnection: viz, barrier scan, and where the residual actually lives (2026-10-02)

Investigation of the Q_H=2 -> Q_H=3 transition (forced pion emission, Paper XVIII), triggered by building a
three.js visualization and then asking whether the "break" is real and whether the small prediction residuals
are "energy drawn from the surrounding field on interaction/measurement" (user hypothesis).

## VISUALIZATION (`papers/src_paper18/t22_t23_transition.html`) -- DONE
Three.js morph of the exact centerline Gamma(t)=((R0+r0 cos(q t))cos2t,(R0+r0 cos(q t))sin2t, r0 sin(q t)),
smoothly interpolating the CLOSED T(2,2) (ring, q=2) and T(2,3) (trefoil, q=3) endpoints (gapless; exact
Gamma_q is open for non-integer q so we blend the two closed endpoints). Matches the img_paper18/
n_initial_q2.{00..90} progression (ring folds into trefoil). Features: fiber-winding rainbow overlay
(Phi=chi+q t) vs R/G/B quark-colour toggle; reconnection PINCH+flash at tau=0.5 (tube necks, white flash);
Q_H=1 pion (unknot ring) detaching AT the pinch, absorb(2+1->3)/emit(2+2->3+1) toggle; persistent field
TETHER (thins with separation, opacity floor -> one continuous condensate never fully severs); play/pause/
reverse/scrub. Reuses hopfion_trefoil.html conventions (R0=3, r0=0.874, Frenet/parallel-transport frame).

## PHYSICS OF THE "BREAK" (discussion, established)
- The break is REAL and LOCAL. Q_H is a homotopy invariant => NO smooth deformation of a closed field takes
  Q_H 2->3 (crossing_transition_v2.py's v1 bug). The transition MUST pass through a reconnection event (field
  degenerate/singular on a local sheet) with an energy barrier. Local, energy-costed, where the Q_H=1 is forced.
- The connection is ALSO REAL and CONTINUOUS. One shared condensate (Bell premise, n-hat single shared field):
  trefoil + pion are textures in the SAME fabric, coupled through overlapping gradient tails; connection decays
  with distance but never severs at the medium level. The animation's tether captures this.
- Net: a continuous energetic process punctuated by one localised topological reconnection.

## SMOKE TESTS (crossing_transition_v2.py) -- TWO WALLS FOUND
1. WHITEHEAD-CHARGE RESOLUTION WALL. The Whitehead integral Q_H=(1/4pi^2)INT A.(curl A) of a KNOWN Q_H=3
   trefoil reads: 0.090 (N=32,wh=24), 0.984 (N=48,wh=48), 1.361 (N=64,wh=64). Converging toward 3 but SLOWLY;
   extrapolates to needing N>=200 (impractical). => CANNOT reliably read the absolute topological charge, so
   cannot watch Q_H go 2->3 or timestamp the reconnection by charge. (Same tube-resolution wall as prior O2.)
2. NO SPONTANEOUS TRANSITION. Plain gradient flow from the Q_H=2 near-trefoil just relaxes (stable); script's
   own verdict "the transition requires external energy input." Consistent with Q_H=3 being a SADDLE (barrier
   between sectors) -- downhill flow rolls back to the Q_H=2 min. (Confirms the no-stable-isolated-quark story.)

## BARRIER SCAN (gradient_flow_constrained_q.py, N=64, 2000 steps, topology-safe) -- DONE
Scanned q_pol=2.0..3.0 (7 pts), relaxed each (two-phase + 9deg/step angle clamp), recorded E_geom=K*J4.
All held topology (r_bar ~3.1-3.2, no HALT, 0% clamped). Cheap: ~2 min/run, ~15 min total (NOT rented hw;
the "hours-days" docstring was for larger N / old doubled version).
  q_pol : E_geom(x1e6) : J4/K
  2.0   : 0.314 : 0.249   (ring / Hopf link)
  2.2   : 0.643 : 0.291
  2.4   : 1.535 : 0.324
  2.5   : 1.525 : 0.302
  2.6   : 1.511 : 0.293
  2.8   : 1.549 : 0.295
  3.0   : 1.565 : 0.291   (trefoil T(2,3))
FINDING: monotonic CLIMB from ring (0.31e6) to a trefoil-level PLATEAU (~1.5e6); NO barrier hump above the
trefoil along the geometric path. Trefoil ~5x the ring in bare shape energy; climb steepest at q=2.2-2.4.
Caveats: NOT fully converged (|grad|~1.5e5 at step 2000, scale-inv functional relaxes slowly) -> fine
structure (tiny q=2.4 bump) within noise, no reliable small barrier; geometric path != true MEP (upper bound);
E_geom is BARE (beta=0, no feedback). ROBUST result: ring->trefoil is uphill ~5x, no spontaneous route --
independent, topology-safe corroboration of "quark has no stable isolated saddle; formation needs energy input."

## THE CALIBRATION FINDING -- E_geom IS NOT THE MASS (the key redirect)
Tried to convert the formation energy to MeV via the baryon calibration. It does NOT convert, and finding that
is the result. The framework's mass formula (derive_C_formation.py) is
    m = Lcond^sector * phi^{2 Q_H} * exp(CS spoke),  Lcond^sector = T_CMB (pi^2/(15 Q_group))^{1/4}.
Baryon: Lcond^baryon ~ 1.61e-4 eV; constituent = Lcond^baryon phi^6 e^{8pi} ~ 236 MeV; completion = /phi^6 ~
13.2 MeV. The FADDEEV gradient-flow energy E_geom=K*J4 DOES NOT ENTER this formula at all. Proof they are not
proportional: ring/trefoil SHAPE-energy ratio = 0.20, but lepton/baryon MASS ratio = ~5e-4 -- off by ~400x.
So E_geom (scale-invariant shape number) and the physical mass are DIFFERENT CURRENCIES; no legitimate factor
converts the barrier to MeV. => THE MASS, AND THE RESIDUAL, LIVE IN THE Lcond*phi-tower*exp(CS) MACHINERY,
NOT IN THE GRADIENT-FLOW SHAPE ENERGY.

## HYPOTHESIS VERDICT: right in spirit, wrong location
User hypothesis (residual = energy drawn from the field/medium on interaction): CORRECT IN SPIRIT -- the
residual IS a dressing/interaction correction from the condensate. But it lives in the EXPONENT of the tower
mass formula, not in a reconnection barrier. derive_C_formation.py states it: "the ~10% is the quark sector's
non-perturbative ~alpha_s residual, a T-matrix/higher-order correction analogous to the lepton's exp(-3/400)."
So: lepton dressing = e^{-3/400 - alpha/pi} (radiative); quark dressing = the analogous alpha_s/WZW-T-matrix
correction to the exp(CS-spoke) exponent. The reconnection barrier and the mass residual are in DIFFERENT
sectors of the theory. For the isospin M_dn/M_up residual specifically, the breathing calc (0.153, see
o3_isospin_residual.md) is the right handle -- it lives in the tower/tangent geometry, not gradient flow.

## DO NOT RE-WALK (for the next session)
- Do NOT try to read Q_H via the Whitehead integral to track the transition: resolution wall (needs N>=200).
- Do NOT try to extract the mass residual from an E_geom gradient-flow barrier: wrong currency (E_geom != mass;
  masses are Lcond*phi^{2Q}*exp(CS)). The residual is a tower-EXPONENT (alpha_s/T-matrix) correction.
- Banked positives: (1) viz done; (2) topology-safe ring->trefoil formation energy ~5x uphill, no spontaneous
  transition (corroborates no-stable-isolated-quark); (3) the clean redirect above.
- Live lead for the residual (if pursued): compute the quark's alpha_s/WZW-T-matrix higher-order correction to
  exp(CS spoke), analogous to the lepton's e^{-3/400} -- this is the OPEN piece derive_C_formation.py names
  (the E_6 tower level Dn_q(g) / the T-matrix residual).
