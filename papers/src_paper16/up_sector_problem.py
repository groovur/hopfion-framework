#!/usr/bin/env python3
r"""
Characterise the up-sector mass problem (2026-09-03). Is the up-type scatter a scale/running artefact,
or a genuine failure of the E_8 route m_q/m_ell(g)=exp(n_q pi/9)?  Test with SCALE-INVARIANT mass RATIOS.
See notes/quark_generation_e8_ribbon_twist.md.
"""
import math
phi=(1+5**0.5)/2
# PDG-ish MS-bar masses (MeV); scales noted -- ratios within a sector are ~RG-invariant (same gamma).
m = {'u':2.16,'d':4.67,'s':93.4,'c':1275.0,'b':4180.0,'t':172570.0}
ml = {1:0.51099895, 2:105.6583755, 3:1776.86}      # e, mu, tau
# E_8 assignment: (generation, exponent) ; n_q = Base_g - 2m ; Base = {41,26,17}
E8 = {'u':(1,19),'c':(2,11),'t':(3,1),'d':(1,17),'s':(2,13),'b':(3,7)}
Base = {1:41,2:26,3:17}
K = 9/math.pi                                        # n = (9/pi) ln(ratio)

print("="*72); print("UP-SECTOR PROBLEM: scale-invariant mass-ratio test"); print("="*72)

def nreq(q):
    g=E8[q][0]; return K*math.log(m[q]/ml[g])
def npred(q):
    g,mm=E8[q]; return Base[g]-2*mm

print("\n[1] Required vs predicted n_q  (n_q = (9/pi) ln(m_q/m_ell(g))):")
print("     q  gen  m_exp  n_pred  n_req   miss")
for q in ['u','c','t','d','s','b']:
    print(f"    {q:>2}   {E8[q][0]}   {E8[q][1]:>3}    {npred(q):>+4}   {nreq(q):>+5.2f}   {nreq(q)-npred(q):>+5.2f}")

print("\n[2] SCALE-INVARIANT within-sector mass RATIOS (predicted vs PDG):")
def ratio_pred(qh,ql):   # m_qh/m_ql predicted = (m_l(gh)/m_l(gl)) exp((n_h-n_l) pi/9)
    gh=E8[qh][0]; gl=E8[ql][0]
    return (ml[gh]/ml[gl])*math.exp((npred(qh)-npred(ql))*math.pi/9)
for (qh,ql) in [('c','u'),('t','c'),('s','d'),('b','s')]:
    pr=ratio_pred(qh,ql); ob=m[qh]/m[ql]
    print(f"    m_{qh}/m_{ql}:  predicted {pr:9.2f}   PDG {ob:9.2f}   pred/PDG {pr/ob:.3f}")
print("    => DOWN ratios (s/d, b/s) within ~10%; UP ratios (c/u, t/c) off by ~2x. Ratios are")
print("       ~scale-invariant, so the up-sector failure is REAL, not a running/scale artefact.")

print("\n[3] What the up-sector REQUIRES (from the ratios):")
nu,nc,nt = nreq('u'),nreq('c'),nreq('t')
print(f"    required up-type n = {{{nu:.2f}, {nc:.2f}, {nt:.2f}}}   increments {{{nc-nu:+.2f}, {nt-nc:+.2f}}}")
print(f"    formula up-type  n = {{3, 4, 15}}                 increments {{+1, +11}}")
print(f"    required up-type exponents m=(Base-n)/2 = {{{(41-nu)/2:.2f}, {(26-nc)/2:.2f}, {(17-nt)/2:.2f}}}")
print("    => the fitted exponents (19,11,1) are wrong for the RATIOS; required ~{18.4, 9.4, 1.9}.")
print("       charm is the worst miss (needs ~9.4, got 11). The clean increments {~3,~6} are suggestive.")

print("\n[4] Down-type, for contrast:")
nd,ns,nb = nreq('d'),nreq('s'),nreq('b')
print(f"    required down n = {{{nd:.2f}, {ns:.2f}, {nb:.2f}}}   vs formula {{7,0,3}}  miss {{{nd-7:+.2f},{ns-0:+.2f},{nb-3:+.2f}}}")
print("    => down miss is ~UNIFORM (-0.5): a single overall shift (running/normalisation). Down WORKS.")

print("\n[5] READ / hypotheses to test:")
print("    The down sector is anchored (strange=muon) and its ratios/hierarchy come out right up to a")
print("    uniform shift. The up sector has NO anchor, its exponent assignment was fitted to mass ORDER")
print("    (not ratios), and the ratios fail by ~2x -- a genuine problem. Candidate resolutions:")
print("     (a) wrong up-type exponent assignment: fit exponents to RATIOS, see if a clean set exists;")
print("     (b) up-type references the WRONG lepton: I_3=+1/2 up-quarks vs the I_3=+1/2 doublet member")
print("         (neutrino), not the charged lepton -- a different base;")
print("     (c) the TOP is special (electroweak-scale, Yukawa ~1) and breaks the tower for gen-3 up.")
print("="*72)
