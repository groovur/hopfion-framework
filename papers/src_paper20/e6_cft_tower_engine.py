#!/usr/bin/env python3.11
"""
e6_cft_tower_engine.py  --  exact WZW/CFT modular-data engine for the quark tower
=================================================================================
GOAL (O2, Paper XX): pin the E_6-native up-quark tower level ALGEBRAICALLY from
conformal data, instead of from a field saddle (which the resolution wall blocks).
This is the non-circular route flagged after the feedback-solver work: the
framework already DERIVES the quark CHARGES from (E_6)_1 conformal weights
(Paper XV Thm P15:thm:charge, h_(1,0)=1/3=|Q_d|, h_27=2/3=|Q_u|), so extending
"conformal weight -> observable" from charge to MASS is the natural program.

GROUNDED INPUTS (verified against the papers, NOT guessed):
  * T_g is a T-MATRIX PHASE (Paper XVII P17:prop, l.601-627):
      T_g = h_{1/2}(k_g) - c(k_g)/24  for SU(2)_{k_g},  g<->k_g=1,2,3
          = 3/(4(k+2)) - k/(8(k+2)) = (6-k)/(8(k+2)) = {5/24, 1/8, 3/40}.
  * Mass tower (Paper XVII l.785,817):  m ∝ exp(-2 pi Q_group T_g).
  * E_8-bridge (down) tower (Paper XVII l.1460):
      Base_g = 180 T_g + 7/2   (180 = 6 h(E_8) = 6*30),   n_q = Base_g - 2 m,
      m = E_8 Coxeter exponent of the quark;  down phase pi/9 = pi Q/(k h),
      Q=10, k=3, h(E_8)=30.  Anchored to the MEASURED leptons (E_8 = lepton bridge).
  * E_6 Coxeter exponents {1,4,5,7,8,11} (h(E_6)=12), split into isospin-doublet
      mirror pairs summing to h=12:  (d,u)=(5,7), (s,c)=(4,8), (b,t)=(1,11).
      down = {d:5,s:4,b:1}, up = {u:7,c:8,t:11}.
  * E_6-native (up) phase pi Q/(k h) = pi*3/(1*12) = pi/4  (Q=Q_group^baryon=3,
      k=1 [SU(3)_1=(E_6)_1], h(E_6)=12).  Steepness ratio (pi/4)/(pi/9)=9/4=2.25
      ~ observed up/down 2.13 (Paper XX, sec.24).

OPEN KNOBS (what O2 must pin): for the E_6-NATIVE up tower, the analogues of the
E_8 Base coefficient (180=6h(E_8) -> 6h(E_6)=72 ?), the offset (7/2 -> ?), and the
ANCHOR (up is anchored HIGH at the top-transition v/sqrt2 ~ 174 GeV ~ m_t, sec.23/24;
down is anchored LOW at the leptons).  KEY: in mass RATIOS to the top anchor, the
Base offset CANCELS, so m_u/m_t, m_c/m_t, m_u/m_c are PARAMETER-FREE predictions of
the grounded data -- the real test.

The engine is exact (sympy rationals; floats only at the final MeV comparison), and
VALIDATES itself (S unitarity, Verlinde fusion -> nonneg integers, known h/c, and it
must reproduce T_g={5/24,1/8,3/40} and the down d,s,b tower) BEFORE the open probe.
"""
import sympy as sp

phi = (1 + sp.sqrt(5)) / 2
pi  = sp.pi

# ══════════════════════════════════════════════════════════════════════════════
# 1. WZW modular data
# ══════════════════════════════════════════════════════════════════════════════
def su2k(k):
    """SU(2)_k modular data (exact). Primaries labelled by l=2j in 0..k."""
    c = sp.Rational(3*k, k+2)
    h = {l: sp.Rational(l*(l+2), 4*(k+2)) for l in range(k+1)}     # h_j=j(j+1)/(k+2), l=2j
    Tphase = {l: sp.nsimplify(h[l] - c/24) for l in range(k+1)}    # T-matrix phase (mod 1)
    # S-matrix  S_{l,l'} = sqrt(2/(k+2)) sin(pi (l+1)(l'+1)/(k+2))
    S = sp.Matrix(k+1, k+1, lambda a, b:
                  sp.sqrt(sp.Rational(2, k+2)) * sp.sin(pi*(a+1)*(b+1)/(k+2)))
    return dict(k=k, c=c, h=h, Tphase=Tphase, S=sp.simplify(S))

def su3_level1():
    """SU(3)_1 primaries: 1=(0,0), 3=(1,0), 3bar=(0,1). h=(p^2+pq+q^2+3p+3q)/12."""
    def hpq(p, q): return sp.Rational(p*p + p*q + q*q + 3*p + 3*q, 12)
    c = sp.Integer(2)
    prim = {'1': (0, 0), '3': (1, 0), '3bar': (0, 1)}
    h = {name: hpq(*pq) for name, pq in prim.items()}
    T = {name: sp.nsimplify(h[name] - c/24) for name in prim}
    return dict(c=c, h=h, Tphase=T)

def e6_level1():
    """(E_6)_1 primaries: 1, 27, 27bar. c=6, h_27=2/3 (Z_3 simple currents)."""
    c = sp.Integer(6)
    h = {'1': sp.Integer(0), '27': sp.Rational(2, 3), '27bar': sp.Rational(2, 3)}
    T = {name: sp.nsimplify(h[name] - c/24) for name in h}
    return dict(c=c, h=h, Tphase=T)

# Lie-theoretic Coxeter data
COX = dict(
    E6=dict(h=12, exponents=[1, 4, 5, 7, 8, 11]),
    E8=dict(h=30, exponents=[1, 7, 11, 13, 17, 19, 23, 29]),
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. VALIDATION  (must pass before any probe is trusted -- CLAUDE.md sec.1)
# ══════════════════════════════════════════════════════════════════════════════
def verlinde_fusion(S):
    """N_{ab}^c = sum_i S_ai S_bi conj(S_ci)/S_0i. Return as nested dict of ints."""
    n = S.shape[0]
    N = {}
    for a in range(n):
        for b in range(n):
            for cc in range(n):
                val = sum(S[a, i]*S[b, i]*sp.conjugate(S[cc, i])/S[0, i] for i in range(n))
                N[(a, b, cc)] = sp.simplify(sp.nsimplify(sp.re(sp.expand(val))))
    return N

def validate():
    print("="*78)
    print("  VALIDATION  (engine correctness before probes)")
    print("="*78)
    ok = True

    # (a) SU(2)_k: known c, h_{1/2}, and S unitarity + Verlinde integrality
    for k in (1, 2, 3):
        d = su2k(k)
        h_half = d['h'][1]                       # l=1 <-> j=1/2
        exp_h = sp.Rational(3, 4*(k+2))
        SSt = sp.simplify(d['S']*d['S'].T)
        unit = SSt == sp.eye(k+1)
        Nf = verlinde_fusion(d['S'])
        allint = all(v.is_integer and v >= 0 for v in Nf.values())
        print(f"  SU(2)_{k}: c={d['c']}, h_(1/2)={h_half} (exp {exp_h}: "
              f"{'OK' if h_half==exp_h else 'FAIL'}), S unitary: {unit}, "
              f"Verlinde N>=0 int: {allint}")
        ok &= (h_half == exp_h) and unit and allint

    # (b) reproduce Paper XVII T_g = {5/24, 1/8, 3/40} at j=1/2, k_g=1,2,3
    Tg = {g: su2k(g)['Tphase'][1] for g in (1, 2, 3)}
    exp_Tg = {1: sp.Rational(5, 24), 2: sp.Rational(1, 8), 3: sp.Rational(3, 40)}
    print(f"  T_g (SU(2)_k, j=1/2): {[str(Tg[g]) for g in (1,2,3)]}  "
          f"(Paper XVII {[str(exp_Tg[g]) for g in (1,2,3)]}: "
          f"{'OK' if Tg==exp_Tg else 'FAIL'})")
    ok &= (Tg == exp_Tg)

    # (c) SU(3)_1 / (E_6)_1 charges = conformal weights (Paper XV)
    s3, e6 = su3_level1(), e6_level1()
    print(f"  SU(3)_1: h_3={s3['h']['3']} (|Q_d|=1/3: {'OK' if s3['h']['3']==sp.Rational(1,3) else 'FAIL'}), "
          f"(E_6)_1: h_27={e6['h']['27']} (|Q_u|=2/3: {'OK' if e6['h']['27']==sp.Rational(2,3) else 'FAIL'})")
    ok &= (s3['h']['3'] == sp.Rational(1, 3)) and (e6['h']['27'] == sp.Rational(2, 3))

    # (d) E_6 exponents mirror-pair to h=12 (isospin doublets)
    pairs = [(5, 7), (4, 8), (1, 11)]
    mir = all(a+b == COX['E6']['h'] for a, b in pairs) and \
          sorted([x for p in pairs for x in p]) == COX['E6']['exponents']
    print(f"  E_6 exponents {COX['E6']['exponents']} mirror to h=12 as doublets "
          f"(d,u)(s,c)(b,t): {'OK' if mir else 'FAIL'}")
    ok &= mir

    print(f"\n  VALIDATION {'PASSED' if ok else 'FAILED'} -- "
          f"{'probes trustworthy' if ok else 'DO NOT trust probes'}")
    return ok

# ══════════════════════════════════════════════════════════════════════════════
# 3. PROBE B (validation of the grounded machinery): E_8-bridge DOWN tower
# ══════════════════════════════════════════════════════════════════════════════
def probe_B_down():
    print("\n"+"="*78)
    print("  PROBE B: E_8-bridge DOWN tower (reproduce d,s,b; grounded machinery)")
    print("="*78)
    Tg   = {g: su2k(g)['Tphase'][1] for g in (1, 2, 3)}
    Base = {g: 180*Tg[g] + sp.Rational(7, 2) for g in (1, 2, 3)}     # {41,26,17}
    mE8  = {1: 17, 2: 13, 3: 7}                                       # d,s,b E_8 exponents
    nq   = {g: Base[g] - 2*mE8[g] for g in (1, 2, 3)}                 # {7,0,3}
    mlep = {1: sp.Float('0.51099895'), 2: sp.Float('105.6583755'), 3: sp.Float('1776.86')}
    meas = {1: 4.67, 2: 93.4, 3: 4180.}
    print(f"  Base_g=180 T_g+7/2 = {[Base[g] for g in (1,2,3)]}   n_q=Base-2m = {[nq[g] for g in (1,2,3)]}")
    print(f"  {'q':>3} {'n_q':>4} {'m_pred (MeV)':>13} {'m_meas':>9} {'ratio':>6}")
    for g, name in [(1,'d'),(2,'s'),(3,'b')]:
        mp = float(mlep[g]*sp.exp(nq[g]*pi/9))
        print(f"  {name:>3} {int(nq[g]):>4} {mp:>13.3f} {meas[g]:>9.2f} {mp/meas[g]:>6.3f}")
    print("  (expect uniform ~1.2 = the down-sector RG running; confirms sec.24 RESULT 2)")

# ══════════════════════════════════════════════════════════════════════════════
# 4. PROBE C (THE OPEN ONE): E_6-native UP tower, PARAMETER-FREE ratios to top
# ══════════════════════════════════════════════════════════════════════════════
def probe_C_up():
    print("\n"+"="*78)
    print("  PROBE C: E_6-native UP tower -- parameter-free mass RATIOS (O2 test)")
    print("="*78)
    Tg    = {g: su2k(g)['Tphase'][1] for g in (1, 2, 3)}
    mE6up = {1: 7, 2: 8, 3: 11}                       # u,c,t E_6 Coxeter exponents
    A     = 6*COX['E6']['h']                          # 72 = 6 h(E_6), analogue of 180=6 h(E8)
    phase = pi*3/(1*COX['E6']['h'])                   # pi/4, E_6-native phase pi Q/(k h)
    # n_q(g) = A T_g + B - 2 m ; the offset B CANCELS in n_q(g)-n_q(3):
    nq_noB = {g: A*Tg[g] - 2*mE6up[g] for g in (1, 2, 3)}
    dnq    = {g: sp.nsimplify(nq_noB[g] - nq_noB[3]) for g in (1, 2, 3)}
    meas   = {1: 2.16, 2: 1270., 3: 172760.}          # u,c,t current masses (MeV)
    print(f"  A=6h(E_6)={A}, phase=pi Q/(k h)={phase}, up E_6 exponents {mE6up}")
    print(f"  n_q(g)-n_q(3) (B cancels) = {[str(dnq[g]) for g in (1,2,3)]}")
    print(f"\n  PARAMETER-FREE ratios (anchor = top, m_t={meas[3]/1e3:.1f} GeV = v/sqrt2):")
    print(f"  {'q':>3} {'Dn_q':>7} {'m/m_t pred':>12} {'m/m_t meas':>12} {'pred/meas':>10}")
    for g, name in [(1,'u'),(2,'c'),(3,'t')]:
        for sign, tag in [(-1, '')]:                  # m ∝ exp(sign * n_q * phase)
            r_pred = float(sp.exp(sign*dnq[g]*phase))
            r_meas = meas[g]/meas[3]
            print(f"  {name:>3} {str(dnq[g]):>7} {r_pred:>12.3e} {r_meas:>12.3e} "
                  f"{(r_pred/r_meas if r_meas else float('nan')):>10.3f}")
    # the internal u/c ratio (independent of the anchor entirely)
    r_uc_pred = float(sp.exp(-(dnq[1]-dnq[2])*phase))
    r_uc_meas = meas[1]/meas[2]
    print(f"\n  INTERNAL u/c ratio (anchor-free): pred {r_uc_pred:.3e} vs meas {r_uc_meas:.3e} "
          f"-> {r_uc_pred/r_uc_meas:.3f}")
    print("  Read: if u/c matches but m/m_t is off by a UNIFORM factor, the tower SHAPE is")
    print("  grounded and only the top-anchor normalisation (B / reference) needs one number.")

if __name__ == '__main__':
    if validate():
        probe_B_down()
        probe_C_up()
