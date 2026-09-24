#!/usr/bin/env python3.11
"""
E_6-NATIVE quark mass formula test (user: quark sector is E_6, not E_8-bridge; top=transition).
Session 2026-09-21. Follows quark_generation_e8_ribbon_twist.md sec.23.

E_8 lepton-bridge route (established, works for DOWN; P17:eq:e8_quark_ratio):
   m_q/m_ell(g) = exp((Base_g - 2 m_E8) * pi/9),  Base_g = 180 T_g + 7/2,
   phase pi/9 = pi Q/(k h),  Q=10 (lepton Q_group), k=3, h(E_8)=30,  T_g=(6-g)/(8(g+2)).
E_6-NATIVE analogue (quark's own group): Q=Q_group^baryon=3, k=1 (SU(3)_1=(E_6)_1), h(E_6)=12
   -> phase_E6 = pi*3/(1*12) = pi/4.   Base'_g = 6 h(E_6) T_g + offset = 72 T_g + offset.

Tests: (1) phase-unit ratio vs observed steepness; (2) E_8 route reproduces down (sanity);
(3) does the naive E_6-native mirror formula reproduce the isospin doublet splittings?
"""
import numpy as np
PI=np.pi; PHI=(1+5**0.5)/2

# masses MeV; leptons; E_6 & E_8 exponent assignments (note sec.13,17)
m={'u':2.16,'d':4.67,'s':93.4,'c':1270.,'b':4180.,'t':172760.}
m_t_full=334000.
mlep={1:0.511,2:105.658,3:1776.86}
gen={'u':1,'d':1,'s':2,'c':2,'b':3,'t':3}
E8={'d':17,'s':13,'b':7,'u':19,'c':11,'t':1}
E6={'d':5,'s':4,'b':1,'u':7,'c':8,'t':11}     # h(E_6)=12; mirror pairs (5,7)(4,8)(1,11)=doublets
Tg={g:(6-g)/(8*(g+2)) for g in (1,2,3)}         # = 5/24, 1/8, 3/40

print("T_g =",{g:f'{Tg[g]:.4f}' for g in (1,2,3)}, " Base_g(E8)=",
      {g:180*Tg[g]+3.5 for g in (1,2,3)})

print("\n"+"="*70,"\n(1) PHASE-UNIT RATIO  (why up is steeper)\n","="*70,sep="")
ph_E8=PI*10/(3*30); ph_E6=PI*3/(1*12)
print(f"  E_8 bridge phase = pi Q/(k h) = pi*10/(3*30) = pi/9 = {ph_E8:.5f}")
print(f"  E_6 native phase = pi*3/(1*12) = pi/4          = {ph_E6:.5f}")
print(f"  ratio (E_6 native)/(E_8 bridge) = {ph_E6/ph_E8:.4f}  (= 9/4 = 2.25)")
print(f"  OBSERVED up/down tower steepness ratio (sec.22) = 2.13   -> match to {100*abs(2.25/2.13-1):.0f}%")
print("  => E_6-native up tower is intrinsically steeper by the phase-unit ratio. RIGHT SHAPE.")

print("\n"+"="*70,"\n(2) E_8 lepton-bridge reproduces DOWN (sanity)\n","="*70,sep="")
for q in ['d','s','b','u','c','t']:
    g=gen[q]; nq=180*Tg[g]+3.5-2*E8[q]
    pred=mlep[g]*np.exp(nq*PI/9)
    meas=m_t_full if q=='t' else m[q]
    print(f"  {q}: pred={pred:9.2f}  meas={meas:9.2f}  pred/meas={pred/meas:.3f}")

print("\n"+"="*70,"\n(3) NAIVE E_6-native: does mirror-exponent give the doublet splitting?\n","="*70,sep="")
print("  E_6 mirror pairs = isospin doublets; m_up/m_down = exp(-2(m_E6^up - m_E6^dn)*phase)")
print("  Solve for the phase each doublet REQUIRES (should be constant = pi/4 if the formula holds):")
for (up,dn) in [('u','d'),('c','s'),('t','b')]:
    dexp=E6[up]-E6[dn]
    for label,mup in [('meas',m[up])]+([('full',m_t_full)] if up=='t' else []):
        ratio=mup/m[dn]
        req_phase=np.log(ratio)/(-2*dexp) if ratio<1 else np.log(ratio)/(2*dexp)
        # signed: ln(m_up/m_dn) = +-2 dexp * phase
        req_phase=np.log(ratio)/(2*dexp)   # keep sign
        print(f"  {up}/{dn}: dexp={dexp:2d}  m_up/m_dn={ratio:7.3f} ({label})"
              f"  required phase = {req_phase:+.4f}   (pi/4={PI/4:.4f})")
print("  => if these are NOT ~pi/4 and NOT constant, the naive E_6 mirror formula FAILS")
print("     for absolute/isospin splittings (steepness trend != full spectrum).")

print("\n"+"="*70,"\nVERDICT\n","="*70,sep="")
print("  - phase-unit ratio 9/4=2.25 ~ observed steepness 2.13 (5.6%): E_6-native explains WHY")
print("    up is steeper (quark's own group, larger phase unit). BEST-grounded steepness lead.")
print("  - naive E_6 mirror formula does NOT cleanly give absolute masses / doublet splittings")
print("    (required phase not constant, not pi/4) -> E_6-native is not turnkey; the reference")
print("    mass, Base' offset, and T_g source (lepton phases vs SU(3)_1 twist) must be pinned.")
print("  - top=transition (v/sqrt2) is the up-sector's upper anchor, not a tower member.")
