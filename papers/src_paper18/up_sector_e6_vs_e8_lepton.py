#!/usr/bin/env python3.11
"""
UP sector: is E_6 still in the running / E_8 as a base, and does a lepton-mass bookkeeping supply the
'boost' the E_8 route misses? (user 2026-09-21: down-from-lepton makes sense; up? full lepton mass
must go somewhere.) Exploratory, honest -- flag coincidences (CLAUDE.md 3).
"""
import numpy as np
PI=np.pi; PHI=(1+5**0.5)/2
mlep={1:0.51099895,2:105.6583755,3:1776.86}          # e, mu, tau (MeV)
mup ={1:2.16,2:1270.,3:172760.}; mup_full3=334000.    # u,c,t ; top full-formation
mdn ={1:4.67,2:93.4,3:4180.}
# E_8 exponents (lepton-bridge): up {19,11,1}, down {17,13,7}; Base_g={41,26,17}
Base={1:41,2:26,3:17}; E8up={1:19,2:11,3:1}; E8dn={1:17,2:13,3:7}
E6up={1:7,2:8,3:11}                                   # E_6-native up exponents (note sec.13)

print("="*70,"\n(1) E_8 route for UP (references leptons e,mu,tau): too light?\n","="*70,sep="")
print("  g  m_up     E8_pred   pred/meas   boost=meas/pred")
boost={}
for g in [1,2,3]:
    nq=Base[g]-2*E8up[g]; pred=mlep[g]*np.exp(nq*PI/9)
    meas=mup[g]; boost[g]=meas/pred
    tag=" (top: full-form 334GeV, pred/full=%.2f)"%(pred/mup_full3) if g==3 else ""
    print(f"  {g}  {meas:8.1f}  {pred:8.2f}   {pred/meas:.3f}      {meas/pred:.3f}{tag}")
print(f"  => u,c UNDER-predicted (boost {boost[1]:.2f}, {boost[2]:.2f}); top OK at full-formation.")
print(f"  boost ratio c/u = {boost[2]/boost[1]:.2f}  (E_6/E_8 phase ratio 9/4=2.25? steepness 2.13?)")

print("\n"+"="*70,"\n(2) E_6-native UP (phase pi/4, E_6 exps), anchored at top-transition v/sqrt2\n","="*70,sep="")
# ratio form within up, anchored at top (g3): m(g)/m(t) = exp(-2*(m_E6(g)-m_E6(t))*phase)? test scaling
# Simpler: does E_6 exponent spacing linearize ln(m_up)? (sec.24 said no; re-confirm w/ top-transition)
for topmass,lbl in [(mup[3],'measured top'),(mup_full3,'full-form top')]:
    lm={1:np.log(mup[1]),2:np.log(mup[2]),3:np.log(topmass)}
    s12=(lm[2]-lm[1])/(E6up[2]-E6up[1]); s23=(lm[3]-lm[2])/(E6up[3]-E6up[2])
    print(f"  {lbl}: ln(m_up) slope per E_6 exp:  u->c={s12:.3f}  c->t={s23:.3f}  const? {abs(s12-s23)<0.3}")
print("  => E_6 exponent spacing does NOT linearize up masses either (confirms sec.24).")

print("\n"+"="*70,"\n(3) LEPTON BOOKKEEPING tests ('full lepton mass goes somewhere')\n","="*70,sep="")
print("  test conservation/seesaw relations between up, down, lepton:")
for g in [1,2,3]:
    e=mlep[g]; u=mup[g]; d=mdn[g]
    print(f"  g{g}: m_u*m_d/m_lep^2 = {u*d/e**2:9.2f}   (m_u+m_d)/m_lep = {(u+d)/e:8.2f}   "
          f"m_u*m_d/(m_lep*m_d0) [g] ...")
print("  m_u*m_d/m_lep^2 across g:", [round(mup[g]*mdn[g]/mlep[g]**2,2) for g in [1,2,3]],
      "-> NOT constant => no clean seesaw")
print("  (m_u+m_d)/m_lep across g:", [round((mup[g]+mdn[g])/mlep[g],2) for g in [1,2,3]],
      "-> NOT constant => no clean sum rule")
# test: is the UP BOOST itself a lepton quantity?
print(f"\n  up boost {{{boost[1]:.2f},{boost[2]:.2f},~1}} vs candidates:")
print(f"    m_lep(g)-based? sqrt(m_mu/m_e)={np.sqrt(mlep[2]/mlep[1]):.2f} (vs boost c/u={boost[2]/boost[1]:.2f})")
print(f"    electric-charge^2 up/down (2/3)^2/(1/3)^2=4? ; charge 2?; steepness 2.13?")

print("\n"+"="*70,"\nHONEST READ\n","="*70,sep="")
print("  - E_8 is the BASE for up too (references leptons), but UNDER-predicts u,c: the up sector is")
print("    steeper than the lepton-bridge gives -> E_6 (native, h=12) IS still in the running as the")
print("    STEEPNESS source (sec.24), but E_6 exponent spacing alone does NOT linearize the absolute")
print("    up masses. Top saturates the EW scale (v/sqrt2) = 'full mass' for the transition object.")
print("  - Simple lepton conservation/seesaw/sum bookkeeping does NOT hold (no constant ratio).")
print("  - So 'full lepton mass goes somewhere' is not a clean per-generation mass identity; the up")
print("    excess is the E_6 steepness + the (dynamical, sec.30) field-selection, not a lepton budget.")
