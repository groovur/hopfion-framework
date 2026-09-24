#!/usr/bin/env python3.11
"""
Derive C in the doublet budget (m_up+m_down)=C*m_lepton (sec.34). User picture: the quark relaxes to
the stable environment -- sparse vacuum -> the lepton (Q_H=2). Test structural candidates for C against
PDG data WITH light-quark uncertainties, and the budget/split decomposition. sympy for exact forms.
"""
import numpy as np, sympy as sp
PI=sp.pi; phi=(1+sp.sqrt(5))/2
# PDG masses (MeV) with rough 1-sigma (light quarks are uncertain -> so is C)
mlep={1:0.51099895,2:105.6583755,3:1776.86}
# central, low, high for u,d,s (MSbar 2GeV), c,b
mu={1:(2.16,1.99,2.26),2:(1270.,1250,1290),3:(172760.,172000,173500)}
md={1:(4.67,4.51,4.79),2:(93.4,92.2,94.6),3:(4180.,4160,4200)}

def C_of(g, ui, di): return (mu[g][ui]+md[g][di])/mlep[g]
print("="*66,"\nC = (m_u+m_d)/m_lep, with light-quark uncertainty ranges\n","="*66,sep="")
for g in [1,2]:
    c_c=C_of(g,0,0); c_lo=C_of(g,1,1); c_hi=C_of(g,2,2)
    print(f"  g{g}: C = {c_c:.2f}  (range {c_lo:.2f} - {c_hi:.2f})")
c1=C_of(1,0,0); c2=C_of(2,0,0)
print(f"  confined mean = {(c1+c2)/2:.2f};  spread {abs(c1-c2)/((c1+c2)/2)*100:.1f}%")

print("\n"+"="*66,"\nSTRUCTURAL candidates for C (framework-native)\n","="*66,sep="")
Qb,Qnu,Ql=3,6,10        # Q_group triple {baryon, nu, lepton} (claude-hopfion 8)
cands=[("Q_b + Q_l = 3+10", Qb+Ql), ("2 Q_nu + 1 = 2*6+1", 2*Qnu+1),
       ("E_8 strange-anchor exponent", 13), ("4*pi", float(4*PI)),
       ("Q_l + Q_b (=13) vs 4pi", None), ("Q_nu + 7", Qnu+7), ("phi^2 * 5", float(phi**2*5))]
for lbl,v in cands:
    if v is None: continue
    dev1=abs(v-c1)/c1*100; dev2=abs(v-c2)/c2*100
    print(f"  {lbl:32s} = {v:7.3f}   dev g1 {dev1:4.1f}%  g2 {dev2:4.1f}%")
print("  NOTE: 3+10 = 2*6+1 = 13 = E_8 anchor exp all give 13 (multiple routes -> numerology hazard,")
print("  CLAUDE.md 3). 4pi=12.57 also within the ~few-% data spread. Data CANNOT distinguish (2 points).")

print("\n"+"="*66,"\nBUDGET / SPLIT decomposition\n","="*66,sep="")
print("  m_up = C*m_lep * f_up(g),  f_up = m_u/(m_u+m_d) = the isospin/generation SPLIT")
for g in [1,2,3]:
    fu=mu[g][0]/(mu[g][0]+md[g][0])
    tag=" (bare top)" if g==3 else ""
    print(f"  g{g}: f_up = {fu:.3f}  f_down = {1-fu:.3f}   (C*m_lep budget split){tag}")
print("  => budget C ~ const (lepton-tied); SPLIT f_up grows with gen = the E_6/E_8 steepness (sec.24/28).")
print("  So m_up(g) = C * m_lep(g) * f_up(g): C=lepton budget (new), f_up=established isospin structure.")

print("\n"+"="*66,"\nformation check: does Q_H=2+Q_H=1 -> Q_H=3 give C?\n","="*66,sep="")
print(f"  Q_group: lepton(Q_H=2)={Ql}, nu(Q_H=1)={Qnu}, baryon(Q_H=3)={Qb}.")
print(f"  C~13 = Q_b+Q_l={Qb+Ql} (quark's own order + the lepton it relaxes toward) is the cleanest")
print(f"  framework-native form and matches the user's 'conform to the lepton' picture -- BUT with the")
print(f"  multiple-13 degeneracy + 2-point data + light-quark errors, it is a MOTIVATED candidate, not")
print(f"  a derivation. A real derivation needs the Q_H=2+1->3 formation mass-balance (Paper XVIII O6).")

print("\n"+"="*66,"\nVERDICT\n","="*66,sep="")
print("  C = 13.1 +- ~5% (2 confined doublets, light-quark errors). Best framework-native form:")
print("  C = Q_group^baryon + Q_group^lepton = 3+10 = 13 (quark order + lepton target), matching the")
print("  'quark relaxes to the lepton' picture. NOT pinned vs 4pi / 2Q_nu+1 (all ~13, within errors).")
print("  Derivation route = the formation mass-balance (O6), not available from the 2-point fit.")
