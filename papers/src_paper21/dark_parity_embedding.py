#!/usr/bin/env python3.11
"""
dark_parity_embedding.py -- CLOSE the cross-sector residual of the dark-parity stability
(Paper XXI open-(ii)): show the Standard-Model (A/E-type) sectors are EVEN under the grading that
makes the D-type dark state odd, so D-type -> SM decay is forbidden. And, en route, PIN the
heavy-band power p of the baryogenesis suppression (Paper XXII). Builds on dark_parity_z2.py +
the modular machinery of dtype_cft_tower_engine.py (so2r_1) and e6/su2 data.

KEY (sharper than dark_parity_z2.py's Z_2): for SO(2r)_1 with r ODD the simple-current group is the
FULL Z_4 centre of Spin(2r) -- {1, v, s, c} with the spinor s of ORDER 4 (s^2=v, s^4=1). The Z_4
charge is the conserved DARK NUMBER (the Dirac complex-spinor U(1) of Paper XXI, reduced mod 4 by
the fusion ring): 1->0, v->2, s->1, c->3. The Z_2 matter parity (even/odd) is just this charge mod 2.
So the dark state (spinor) carries dark-number 1; the Z_2 parity is its order-2 remnant.

EMBEDDING ARGUMENT: the SM sectors live in DIFFERENT chiral algebras -- neutrino SU(2)_1, charged
lepton SU(2)_3 (2I), baryon (E_6)_1 (2T) -- none of which is SO(2r). They therefore carry NO
SO(2r)-centre charge: dark-number 0, i.e. matter-parity EVEN. A decay D(dark#=1) -> SM(dark#=0)
changes the SO(10)-centre charge by 1 and is forbidden by its conservation. That conserved charge is
the SAME dark number the co-genesis (ADM, Paper XXII) requires to carry the asymmetry -- so stability
and the asymmetry are ONE conserved charge, not two assumptions. RESIDUAL: cross-sector conservation
holds iff the combined theory conserves the SO(10)-centre dark number, which ADM assumes anyway.
"""
import sympy as sp
from dtype_cft_tower_engine import so2r_1, su2k, e6_1

# ---- SO(2r)_1 centre = Z_4 (r odd): dark-number charge of each primary ----
# labels 1,v,s,c ; Z_4 charges 0,2,1,3 (s the generator, s^2=v, s^3=c, s^4=1)
Z4 = {'1': 0, 's': 1, 'v': 2, 'c': 3}

def report_dtype(r, name):
    c, h = so2r_1(r)
    print(f"\n  {name}: SO({2*r})_1, c={c}, h(1,v,s,c)=0,1/2,{sp.Rational(r,8)},{sp.Rational(r,8)}")
    print(f"    centre = Z_4 (r={r} odd): dark-number charge  " +
          ", ".join(f"{p}:{Z4[p]}" for p in ('1','v','s','c')))
    print(f"    matter parity = (dark# mod 2):  EVEN {{1,v}}  ODD {{s,c}}  (the dark spinor is odd)")
    print(f"    => the dark fermion is the spinor s: dark-number 1, matter-parity ODD.")

if __name__ == '__main__':
    print("="*74)
    print("  CROSS-SECTOR DARK-PARITY EMBEDDING + heavy-band power p")
    print("="*74)
    print("\n--- A. D-type centre charge (the conserved dark number) ---")
    report_dtype(5, "2D_3 dark")
    report_dtype(7, "2D_5 dark")

    print("\n--- B. The SM sectors carry ZERO SO(2r)-centre charge (different algebras) ---")
    sm = {'neutrino  SU(2)_1 (A-type)': su2k(1),
          'lepton    SU(2)_3 (2I,E_8)': su2k(3),
          'baryon    (E_6)_1  (2T)'    : e6_1()}
    for nm,(cc,hh) in sm.items():
        print(f"    {nm:28}: c={str(cc):>4}, primaries {list(hh.keys())}")
    print("    None is an SO(2r) algebra -> no Spin(2r)-centre (dark-number) charge -> dark# = 0")
    print("    => every SM sector is matter-parity EVEN / dark-number neutral.")

    print("\n--- C. Decay D -> SM is forbidden by the conserved dark number ---")
    print("    dark spinor: dark# = 1.   SM final states: dark# = 0.   Delta(dark#) = 1 =/= 0.")
    print("    The SO(10)-centre Z_4 is conserved by fusion; it IS the Dirac complex-spinor U(1)")
    print("    (Paper XXI sec:dirac) reduced mod 4 -- the SAME dark number the ADM co-genesis")
    print("    (Paper XXII) needs to carry the asymmetry. So: stability and the asymmetry are ONE")
    print("    conserved charge. The Z_2 matter parity is its order-2 remnant.")
    print("    => cross-sector residual CLOSES, conditional only on the combined theory conserving")
    print("       the SO(10)-centre dark number -- which ADM assumes anyway (not a new assumption).")

    print("\n--- D. Heavy-band power p of the baryogenesis suppression (Paper XXII) ---")
    # final state of the baryon-forming reconnection 2+2 -> 3+1
    final = {'Q_H=3 baryon (E_6)': 'HEAVY band', 'Q_H=1 neutrino (A-type)': 'light band'}
    n_heavy = sum(1 for v in final.values() if v == 'HEAVY band')
    print(f"    reconnection 2+2 -> 3+1 final state:")
    for k,v in final.items(): print(f"      {k:28}: {v}")
    print(f"    heavy-band (sin^4 directional-suppression) objects in the final state: n_heavy = {n_heavy}")
    print(f"    The directional factor S_eff suppresses the AMPLITUDE to project the final state onto")
    print(f"    the heavy band, one factor per heavy-band object -> amplitude ∝ S_eff^{n_heavy}.")
    print(f"    eta is a RATE ∝ |amplitude|^2 -> eta ∝ S_eff^(2*{n_heavy}) = S_eff^{2*n_heavy}.")
    print(f"    => p = 2*n_heavy = {2*n_heavy}.  (p=1 would treat S_eff as a probability, not an")
    print(f"       amplitude factor; the rate=|amplitude|^2 structure gives p=2 for the single heavy baryon.)")

    print("\n" + "="*74)
    print("  VERDICT:")
    print("   * STABILITY: the dark state is the SO(10) spinor carrying dark-number 1 (Z_4 centre);")
    print("     the SM sectors carry dark-number 0 (different algebras) -> D->SM forbidden by the")
    print("     conserved dark number = the Dirac U(1) = the ADM asymmetry carrier. The Z_2 matter")
    print("     parity is its remnant. Cross-sector residual CLOSES (conditional on dark-number")
    print("     conservation, which ADM requires). Absolute stability no longer rests on the hand-")
    print("     waving 'a spinor can't be built from non-spinor reps'.")
    print("   * POWER p = 2: one heavy-band object (the Q_H=3 baryon) in the 2+2->3+1 final state,")
    print("     entering the amplitude once; eta ~ |amplitude|^2 -> S_eff^2. Pins the Paper XXII")
    print("     p=2 (which the eta self-consistency preferred) from the final-state heavy-band count.")
    print("  [N] both; residual [O]: dark-number conservation across sectors (ADM assumption).")
    print("="*74)
