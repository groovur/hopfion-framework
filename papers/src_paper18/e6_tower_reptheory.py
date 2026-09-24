#!/usr/bin/env python3.11
"""
Representation-theory origin of the down-sector power law m ~ |Tw|^4.2 (sec.29/32).
Session 2026-09-21. Question: is it an algebraic IDENTITY from the E_6/E_8 WZW data, or an emergent
NUMERICAL coincidence of the E_8-corrected lepton masses? Tell: is the per-step exponent EXACTLY
constant? Uses sympy for exact rationals (T_g, exponents, pi/9) + mpmath for the transcendentals.
"""
import sympy as sp

phi=(1+sp.sqrt(5))/2; pi=sp.pi
# --- exact framework data (Paper XVII) ---
gens=[1,2,3]
Tg={1:sp.Rational(5,24),2:sp.Rational(1,8),3:sp.Rational(3,40)}     # T_g=(6-g)/(8(g+2))
Base={g:180*Tg[g]+sp.Rational(7,2) for g in gens}                    # 41,26,17
E8_down={1:17,2:13,3:7}                                              # d,s,b E_8 exponents (mirror pairs)
E6_down={1:5,2:4,3:1}                                                # d,s,b E_6 exponents
Tw={g:abs(E6_down[g]-6) for g in gens}                               # |E_6 twist| = 1,2,5
nq={g:Base[g]-2*E8_down[g] for g in gens}                            # E_8 correction exponent = 7,0,3
print("Tg   =",{g:Tg[g] for g in gens})
print("Base =",{g:Base[g] for g in gens},"  E8_down=",E8_down,"  n_q=",nq,"  |Tw|=",Tw)

# --- measured leptons (inputs; the E_8 route references them) & down masses ---
mlep={1:sp.Float('0.51099895'),2:sp.Float('105.6583755'),3:sp.Float('1776.86')}  # MeV
mdn ={g: mlep[g]*sp.exp(nq[g]*pi/9) for g in gens}
print("\nE_8-route down masses (MeV):",{g:float(mdn[g]) for g in gens},
      " (meas 4.67,93.4,4180)")

# --- power-law fit and, crucially, PER-STEP exponent constancy ---
import math
lnTw={g:math.log(float(Tw[g])) for g in gens}; lnm={g:math.log(float(mdn[g])) for g in gens}
p_ds=(lnm[2]-lnm[1])/(lnTw[2]-lnTw[1])     # d->s  (|Tw| 1->2)
p_sb=(lnm[3]-lnm[2])/(lnTw[3]-lnTw[2])     # s->b  (|Tw| 2->5)
p_db=(lnm[3]-lnm[1])/(lnTw[3]-lnTw[1])
print(f"\nper-step exponent p:  d->s={p_ds:.4f}  s->b={p_sb:.4f}  d->b={p_db:.4f}")
print(f"  constancy: |p_ds - p_sb| = {abs(p_ds-p_sb):.4f}  ({abs(p_ds-p_sb)/p_db*100:.1f}% spread)")
print("  => EXACTLY constant (~0) means algebraic IDENTITY; a % spread means NUMERICAL coincidence.")

# --- the decomposition: down ratio = lepton ratio x E_8 factor ---
print("\ndecomposition  m_down(g2)/m_down(g1) = [m_lep ratio] x exp(dn_q pi/9):")
for (a,b,tw) in [(1,2,'d->s'),(2,3,'s->b')]:
    lr=float(mlep[b]/mlep[a]); e8=float(sp.exp((nq[b]-nq[a])*pi/9)); prod=lr*e8
    twr=float(Tw[b])/float(Tw[a])
    print(f"  {tw}: lepton {lr:8.2f} x E8 {e8:7.4f} = {prod:7.3f} = (|Tw| ratio {twr:.2f})^{math.log(prod)/math.log(twr):.3f}")
print("  the wildly different lepton ratios x E_8 factors give near-equal twist-powers -- the")
print("  'conspiracy'. If exact, the lepton FORMULA must force it; if ~1%, it is a 3-point coincidence.")

# --- is p a rep-theory number? ---
print("\np vs rep-theory constants:")
for lbl,v in [('phi^3',float(phi**3)),('h(E_6)/e',12/math.e),('c(E_6)=6 -> 6*ln2/ln(2.5)?',None),
              ('2phi+1',float(2*phi+1))]:
    if v: print(f"  {lbl:22s} = {v:.4f}  (fit ~{p_db:.3f})")

# --- CROSS-CHECK: do the FORMULA lepton masses (framework tower) also give a clean power? ---
# lepton tower ratio m(g2)/m(g1) = phi^{-2(n_g2-n_g1)} exp(-2 pi*10*(T_g2-T_g1)); absolute n_g irrational.
# Use the framework's exp part alone (the T_g / topological-spin piece) as the rep-theory-internal test:
print("\nrep-theory-INTERNAL (drop measured leptons): down mass from T_g + exponents only")
print("  m_down(g) ~ exp(-2 pi Q T_g) * exp(n_q pi/9), Q=10 (E_8 bridge). check |Tw|^p:")
for Q in [10]:
    mint={g: sp.exp(-2*pi*Q*Tg[g])*sp.exp(nq[g]*pi/9) for g in gens}
    lmi={g:math.log(float(mint[g])) for g in gens}
    pa=(lmi[2]-lmi[1])/(lnTw[2]-lnTw[1]); pb=(lmi[3]-lmi[2])/(lnTw[3]-lnTw[2])
    print(f"    Q={Q}: p(d->s)={pa:.3f}  p(s->b)={pb:.3f}  (constant? {abs(pa-pb)<0.05})")
