#!/usr/bin/env python3.11
"""
top_anchor_offset.py -- fix the anchor SCHEME/SCALE, then measure the top-tower
offset consistently and test whether it equals C (doublet sum rule) or is a
mixed-scale artifact.

CONTEXT (Paper XX, O2): the E_6-native up-tower SPACING is CFT-fixed
(P20:prop:up_ratio, m_u/m_c=e^{-2pi}). The ABSOLUTE anchor is open. Anchoring at
the top, the two light generations sat a factor ~13 below the tower -- BUT that
~13 used m_u at 2 GeV with m_t at the top scale (MIXED scales), exactly the
inconsistency that made the raw m_u/m_c look 9% off. The framework's masses are
MS-bar (electron m_e = Lcond phi^20 e^{4pi-3/400-alpha/pi}, where -3/400=T_{1/2}/Q
is the WZW pole->MSbar conversion, Paper IV Thm 8.11), so the anchor must be MS-bar
and at a COMMON scale. This script recomputes the offset consistently.

Because gamma_m is flavour-independent, m_i/m_j at a common scale is RG-invariant;
the physical top-offset is offset(g) = (m_g/m_t)_obs / (m_g/m_t)_tower at a common
scale. Test: is offset ~ C(g1)=(m_u+m_d)/m_e ~ 13 (a structural link), or does
scheme-consistency move it away (the 13 was a mixed-scale artifact)?

Validated 3-loop MS-bar RG runner (thresholds at m_c,m_b,m_t).
"""
import numpy as np
from scipy.integrate import solve_ivp

MZ = 91.1876; ALPHAS_MZ = 0.1179; ZETA3 = 1.2020569
MC, MB, MT = 1.27, 4.18, 162.5          # MS-bar m_q(m_q) thresholds

def bcoef(nf): return (11-2*nf/3, 102-38*nf/3, 2857/2-5033*nf/18+325*nf**2/54)
def gcoef(nf): return (4.0, 202/3-20*nf/9, 1249-(2216/27+160*ZETA3/3)*nf-140*nf**2/81)

def _run_a(mu_hi, a_hi, mu_lo, nf):
    b0,b1,b2 = bcoef(nf)
    s = solve_ivp(lambda t,y:[-2*y[0]**2*(b0+b1*y[0]+b2*y[0]**2)],
                  [np.log(mu_hi),np.log(mu_lo)],[a_hi],rtol=1e-11,atol=1e-15)
    return s.y[0,-1]

def alpha_s(mu):
    a = ALPHAS_MZ/(4*np.pi)
    # M_Z is between m_b and m_t -> nf=5 there
    if mu >= MT:
        a = _run_a(MZ,a,MT,5); a=_run_a(MT,a,mu,6); return 4*np.pi*a
    if mu >= MB:
        a=_run_a(MZ,a,mu,5); return 4*np.pi*a
    a=_run_a(MZ,a,MB,5)
    if mu >= MC:
        a=_run_a(MB,a,mu,4); return 4*np.pi*a
    a=_run_a(MB,a,MC,4); a=_run_a(MC,a,mu,3); return 4*np.pi*a

def _run_m_fixed(m,mu0,mu1,nf):
    b0,b1,b2=bcoef(nf); g0,g1,g2=gcoef(nf)
    a0=alpha_s(mu0)/(4*np.pi)
    def rhs(t,y):
        a,lnm=y
        return [-2*a**2*(b0+b1*a+b2*a**2), -2*(g0*a+g1*a**2+g2*a**3)]
    s=solve_ivp(rhs,[np.log(mu0),np.log(mu1)],[a0,np.log(m)],rtol=1e-11,atol=1e-15)
    return np.exp(s.y[1,-1])

def run_mass(m,mu0,mu1):
    """Run a quark mass between any two scales, crossing m_c,m_b,m_t with LO
    (mass-continuous) matching and the correct n_f in each interval."""
    thr=[MC,MB,MT]
    pts=sorted(set([mu0,mu1]+[t for t in thr if min(mu0,mu1)<t<max(mu0,mu1)]))
    if mu0>mu1: pts=pts[::-1]
    def nf_between(a,b):
        mid=np.sqrt(a*b); return 3+sum(mid>t for t in thr)
    for lo,hi in zip(pts[:-1],pts[1:]):
        m=_run_m_fixed(m,lo,hi,nf_between(lo,hi))
    return m

# ── VALIDATION ────────────────────────────────────────────────────────────────
print("="*76,"\n  VALIDATION\n","="*76,sep="")
print(f"  alpha_s(2GeV)={alpha_s(2.0):.4f} (~0.30)   alpha_s(MZ)={alpha_s(MZ):.4f} (0.1179)")
mc2=run_mass(MC,MC,2.0); print(f"  m_c(2GeV)={mc2:.4f} (~1.09)")
muMZ_chk=run_mass(2.16e-3,2.0,MZ)
print(f"  m_u(MZ)={muMZ_chk*1e3:.3f} MeV (~1.2-1.3)   [runs 2GeV->MZ across m_b]")

# ── COMMON-SCALE MASSES (MS-bar at M_Z, consistent scheme) ──────────────────────
print("\n"+"="*76,"\n  TOP-TOWER OFFSET at a COMMON scale (M_Z), scheme-consistent\n","="*76,sep="")
# up-type running masses at M_Z
m_u_MZ = run_mass(2.16e-3, 2.0, MZ)      # m_u(2GeV) -> M_Z
m_c_MZ = run_mass(MC,      MC,  MZ)      # m_c(m_c)  -> M_Z
m_t_MZ = run_mass(MT,      MT,  MZ)      # m_t(m_t)  -> M_Z  (MS-bar anchor)
print(f"  m_u(MZ)={m_u_MZ*1e3:.3f} MeV   m_c(MZ)={m_c_MZ:.4f} GeV   m_t(MZ)={m_t_MZ:.2f} GeV")

# E_6-native tower levels (from e6_cft_tower_engine.py)
PHI=(1+5**0.5)/2
Tg={1:5/24,2:1/8,3:3/40}; mE6={1:7,2:8,3:11}; A=72; ph=np.pi/4
n={g:A*Tg[g]-2*mE6[g] for g in (1,2,3)}
def tower_ratio(g):  # m_g/m_t predicted by the tower
    return np.exp(-(n[g]-n[3])*ph)

print(f"\n  tower levels n = {[round(n[g],3) for g in (1,2,3)]}  (n_u-n_t={n[1]-n[3]:.1f}, n_c-n_t={n[2]-n[3]:.1f})")
print(f"  {'q':>3} {'(m/m_t)_obs':>12} {'(m/m_t)_tower':>14} {'OFFSET':>8}")
off={}
for g,name in [(1,'u'),(2,'c')]:
    r_obs = (m_u_MZ if g==1 else m_c_MZ)/m_t_MZ
    r_tow = tower_ratio(g)
    off[g]=r_obs/r_tow
    print(f"  {name:>3} {r_obs:>12.3e} {r_tow:>14.3e} {off[g]:>8.3f}")
print(f"  offset(u)/offset(c) = {off[1]/off[2]:.3f}  (=1 would be a single uniform normalisation)")
off_mean=np.sqrt(off[1]*off[2])
print(f"  => uniform top-offset (geo-mean) = {off_mean:.2f}")

# ── COMPARE: is the MIXED-scale ~13 an artifact? and offset vs C ────────────────
print("\n"+"="*76,"\n  MIXED-SCALE vs COMMON-SCALE, and the C test\n","="*76,sep="")
off_mixed_u = (2.16e-3/172.76)/tower_ratio(1)     # m_u@2GeV, m_t=pole (the old ~13)
print(f"  MIXED-scale offset(u) [m_u@2GeV / m_t=172.76 pole] = {off_mixed_u:.2f}  (the old '~13')")
print(f"  COMMON-scale offset (M_Z, MS-bar)                  = {off_mean:.2f}")
# doublet sum-rule C
mlep={1:0.51099895,2:105.6583755}; md={1:4.67,2:93.4}; mu={1:2.16,2:1097.}  # c@2GeV
C1=(mu[1]+md[1])/mlep[1]; C2=(mu[2]+md[2])/mlep[2]
print(f"\n  C(g1)=(m_u+m_d)/m_e = {C1:.2f}   C(g2)=(m_c+m_s)/m_mu = {C2:.2f}   (2 GeV masses)")
print(f"  offset({off_mean:.1f}) vs C(~{C1:.0f}): {'CLOSE' if abs(off_mean-C1)/C1<0.1 else 'NOT close'} "
      f"-> the offset={off_mean:.1f} is {'~ C' if abs(off_mean-C1)/C1<0.1 else 'NOT C: the mixed-scale 13~C coincidence DISSOLVES under scheme discipline'}")
print(f"\n  disciplined candidates for {off_mean:.2f} (report only, do NOT claim):")
for lbl,v in [('2pi',2*np.pi),('e^2',np.e**2),('phi^4',PHI**4),('4pi',4*np.pi),
              ('C(g1)',C1),('2phi^3+...',None)]:
    if v: print(f"    {lbl:8s}={v:6.2f}  ({'within 5%' if abs(v-off_mean)/off_mean<0.05 else f'{(v/off_mean-1)*100:+.0f}%'})")
print("\n  anchor-scheme sensitivity of the offset (pole vs MS-bar top, ~7% span):")
for lbl,mt in [('MS-bar m_t(MZ)',m_t_MZ),('pole 172.76',172.76),('v/sqrt2 173.9',173.9)]:
    o=(m_u_MZ/mt)/tower_ratio(1); print(f"    anchor={lbl:20s}-> offset(u)={o:.2f}")
