#!/usr/bin/env python3.11
"""
o2_phi_chirality_probe.py -- SKETCH/probe for Paper XX O2 (the generation-1 flip Phi=m_d-m_u SIGN),
using the pion-emission chirality. DM/Quark crossover (src_paper23).

#############################################################################################
# RETRACTED (2026-10-15, advisor catch). This probe compared the WRONG PAIR. The z->-z mirror #
# flips Q_H (+2.770 -> -2.770 in the run below) -> it is the ANTI-trefoil (opposite-charge     #
# sector), NOT the down quark. Up and down are BOTH Q_H=+3, distinguished INTERNALLY: Paper    #
# XVIII/XX assign isospin to the crossing-vertex network (z=+r0, up-type) vs the midpoint/     #
# distal network (z=0, down-type) of the SAME T(2,3) (P20:prop:isospin_scale, M_dn/M_up=0.189).#
# So "writhe +/- = up/down" is UNSUPPORTED and the WZW-phase -> Phi-sign route is WITHDRAWN.    #
# Two reasons it fails: (i) the pair is particle/anti (parity), not isospin; the E_fb          #
# degeneracy only shows the Faddeev energy is PARITY-invariant, nothing about isospin. (ii) a   #
# CP-ODD phase (theta^3=i, or Im(A1 A2*)) gives a RATE asymmetry, not a MASS splitting -- CPT   #
# forces conjugate states to equal mass; mass is CP-EVEN. Also theta_j=e^{2pi i h_j} is ALREADY #
# in the e6 tower (T_g=h_j-c/24) -- not a new handle. O2's Phi stays as Paper XX states it      #
# (charge-linear alpha_em*Lambda, sign anchor-fixed). Kept only to document the dead lead; the  #
# only valid output below is the exact PARITY-invariance of E_fb (minor).                      #
#############################################################################################

THE LEAD (2026-10-15). O2 needs the up/down mass-splitting SIGN. The memory/Paper XX record that the
STATIC-energy route is REFUTED: the Faddeev energy is z->-z symmetric (P20:rem:orientation), so up and
down (the two writhe chiralities of the T(2,3) trefoil -- the Z_2 writhe-asymmetry isospin label of
Paper XVIII) are energetically DEGENERATE absent a chiral term. The Phi sign therefore cannot come from
|E_fb|; it must come from a CHIRAL quantity of the FORMATION -- i.e. the pion-emission reconnection, which
IS chiral (handed threading; the pion carries a twist quantum, Paper XXIII). The named candidate is the
writhe-3 WZW framing phase theta^3_{3/2}=i (Paper XV), already flagged in Paper XXII as "a genuine chiral
phase in the reconnection amplitude."

WHAT THIS PROBE SHOWS (cheap, reuses the engine; it does NOT compute Phi):
 (1) E_fb NULL: the trefoil field and its z->-z mirror have IDENTICAL E_fb, V, |Q_H| -> the sign is NOT
     in the (parity-even) energy. Re-confirms the refutation on the full reaction-coordinate field.
 (2) WRITHE is the chiral LABEL: the writhe sign flips (+ <-> -) between the two chiralities -> up/down
     = writhe +/- . This is the quantity that distinguishes them (the Z_2 writhe asymmetry).
 (3) REFRAME: so Phi's sign = (writhe-sign label) x (sign the WZW framing phase theta^3_{3/2}=i imparts to
     the up-vs-down formation AMPLITUDE). That is an ALGEBRAIC WZW question (Paper XV/III modular data) --
     it SIDESTEPS the resolution wall, exactly as the e6_cft_tower did for the mass-RATIO spacing. This
     is the tractable O2 target: evaluate the chiral framing phase, not relax a field.
CAVEAT [O]: this probe establishes (1)+(2) and PINPOINTS the route (3); it does NOT establish that
theta^3_{3/2}=i gives Phi>0 (down heavier). That is the next, algebraic step (a phase -> a real splitting
needs the amplitude-interference bookkeeping, as for delta_CP=65 deg in Paper XXII). Flagged, not claimed.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from knotted_map import brieskorn_doublet
from torsion_framing import trace_components, gauss_wr_lk
from feedback_saddle_core import nhat_from_doublet, Efb_diagnostics, whitehead_QH

def efb_of(X, Y, Z, h, c=3.0, beta=0.452):
    z1, z2 = brieskorn_doublet(X, Y, Z, p=2, q=3, i=1, j=0, c=c)
    n = nhat_from_doublet(z1, z2)
    d = Efb_diagnostics(n, h, beta)
    return d, whitehead_QH(z1, z2, h)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=96); ap.add_argument('--L', type=float, default=5.0)
    ap.add_argument('--Nb', type=int, default=400)
    a = ap.parse_args()
    ax = np.linspace(-a.L, a.L, a.N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')

    print("="*82)
    print("  O2 Phi-SIGN CHIRALITY PROBE (sketch): is the m_d-m_u sign in the formation chirality?")
    print("="*82)
    # (1) E_fb null: trefoil vs its z->-z mirror
    dp, qhp = efb_of(X, Y,  Z, h)     # chirality +  (up-candidate)
    dm, qhm = efb_of(X, Y, -Z, h)     # chirality -  (down-candidate = z->-z mirror)
    print("  (1) ENERGY NULL (parity-even Faddeev energy -> up/down degenerate):")
    print(f"      chirality +:  E_fb={dp['E_fb']:.4e}  V={dp['V']:.4f}  J4/J2a={dp['J4_over_J2a']:.4f}  Q_H={qhp:+.3f}")
    print(f"      chirality -:  E_fb={dm['E_fb']:.4e}  V={dm['V']:.4f}  J4/J2a={dm['J4_over_J2a']:.4f}  Q_H={qhm:+.3f}")
    rel = abs(dp['E_fb']-dm['E_fb'])/abs(dp['E_fb'])
    print(f"      |dE_fb|/E_fb = {rel:.2e}  -> {'DEGENERATE (sign NOT in the energy)' if rel<1e-6 else 'NOT degenerate (!)'}")
    # (2) writhe is the chiral label: trace the trefoil centerline, writhe +; mirror (z->-z) flips sign
    comps = trace_components(1.0, c3=3.0, Nb=a.Nb)      # trefoil (lambda=1)
    du = 2*np.pi/a.Nb
    wr_p = gauss_wr_lk(comps)[1]
    comps_m = [np.column_stack([c[:,0], c[:,1], -c[:,2]]) for c in comps]   # z->-z mirror
    wr_m = gauss_wr_lk(comps_m)[1]
    print("  (2) WRITHE = the chiral LABEL (flips sign between chiralities = up vs down):")
    print(f"      writhe(+) = {wr_p:+.3f}   writhe(-) = {wr_m:+.3f}   (Z_2 writhe asymmetry, Paper XVIII)")
    print("="*82)
    print("  *** RETRACTED (see header). The z->-z pair above is Q_H=+3 vs Q_H=-3 = trefoil vs")
    print("  ANTI-trefoil (a PARITY / particle-antiparticle pair), NOT up vs down (both Q_H=+3,")
    print("  distinguished internally by the crossing z=+r0 / midpoint z=0 site networks, Paper XVIII/")
    print("  XX, P20:prop:isospin_scale M_dn/M_up=0.189). So the only valid result here is that E_fb is")
    print("  exactly PARITY-invariant; 'writhe=up/down' and the WZW-phase->Phi-sign route are WITHDRAWN")
    print("  (a CP-odd phase gives a rate asymmetry, not a mass splitting -- CPT; and theta_j is already")
    print("  in the e6 tower). O2's Phi stays as Paper XX gives it (alpha_em charge-linear). ***")
    print("="*82)
