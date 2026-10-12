#!/usr/bin/env python3.11
"""
cogenesis_branching.py -- COMPUTE the co-genesis ratio r = n_dark/n_b from the group-theoretic
SECTOR WEIGHTS of the two McKay subgroups, turning the abundance's "r ~ O(1) (comparable
subgroups 12:24)" heuristic (notes (u), adm_abundance.py "ASSUMED, key uncomputed piece") into an
actual number. DM/Quark crossover (src_paper22), notes (vv)->(ww).

STATUS of r before this (checked, code+notes):
  * (q): r_forced = 0  -- the dihedral is NOT a FORCED per-baryon byproduct (P18 pion saturates
    the 2+2->3+1 slot). The only prior DERIVATION, and negative.
  * (u): r ~ O(1) re-introduced as a HEURISTIC from the subgroup orders |2D_3|:|2T| = 12:24.
  * never positively computed.

PRINCIPLE (assumed, stated): in the shared-2I-vacuum co-genesis a single CP-biased reconnection
seeds BOTH sectors; the asymmetry branches in proportion to each sector's WEIGHT in the vacuum.
The natural, measure-independent-up-to-O(1) weights are group-theoretic:
  (a) group ORDER |G| (regular-rep dimension = total DOF of the sector);
  (b) number of IRREPS of G (= McKay/affine-Dynkin node count = number of species).
Baryon sector = 2T (binary tetrahedral, order 24, McKay E_6); dark = 2D_3 (dicyclic order 12,
McKay D_5). r = (dark weight)/(baryon weight). Then Omega_dm/Omega_b = r f (m_dark/m_p) fixes
m_dark (f~1 trap fraction). This is a STATIC sector-weight estimate: the true r needs the
freeze-out DYNAMICS (reconnection amplitudes x phase space), NOT computed here -- so it upgrades r
from "assumed O(1)" to "computed sector-weight O(1), bracketed", not to a derived dynamical rate.
"""
# ---- irrep dimensions (verify |G| = sum d^2) ----
G = {
  '2T  (baryon, McKay E6)': [1,1,1,2,2,2,3],         # binary tetrahedral, order 24
  '2D_3 (dark, McKay D5)' : [1,1,1,1,2,2],           # dicyclic order 12
}
Odm_over_Ob = 0.1200/0.02237    # Planck 2018 = 5.363
m_p = 0.9383                    # GeV

def order(dims): return sum(d*d for d in dims)
def nirrep(dims): return len(dims)

if __name__=='__main__':
    print("="*74)
    print("  CO-GENESIS RATIO r = n_dark/n_b from group-theoretic SECTOR WEIGHTS")
    print("="*74)
    for name,dims in G.items():
        print(f"  {name:28} irrep dims {dims}  |G|=sum d^2={order(dims)}  #irreps={nirrep(dims)}")
    oT = order(G['2T  (baryon, McKay E6)']); oD = order(G['2D_3 (dark, McKay D5)'])
    nT = nirrep(G['2T  (baryon, McKay E6)']); nD = nirrep(G['2D_3 (dark, McKay D5)'])
    assert oT==24 and oD==12, "group orders off"

    r_order  = oD/oT               # 12/24 = 0.50  (total DOF / regular-rep dim)
    r_nirrep = nD/nT               # 6/7  = 0.86   (number of species)
    r_cft    = 4/3                 # level-1 WZW primary count SO(10)_1(4) vs (E6)_1(3) [CFT-assignment dependent]
    print(f"\n  SECTOR-WEIGHT ESTIMATES of r = dark/baryon:")
    print(f"    (a) group order   r = |2D_3|/|2T|   = {oD}/{oT} = {r_order:.3f}")
    print(f"    (b) # irreps      r = {nD}/{nT}            = {r_nirrep:.3f}")
    print(f"    (c) level-1 prim. r = 4/3            = {r_cft:.3f}   [CFT-assignment dependent; cross-check]")
    rs = [r_order, r_nirrep, r_cft]
    print(f"    => r in [{min(rs):.2f}, {max(rs):.2f}]  -- ALL O(1); group-theory (a,b) give r = 0.5-0.9")

    print(f"\n  ABUNDANCE  (Omega_dm/Omega_b = {Odm_over_Ob:.3f}, f=1):  m_dark = (ratio) m_p / r")
    print(f"    {'r':>8} {'m_dark (GeV)':>14}")
    for r in sorted(rs):
        print(f"    {r:>8.2f} {Odm_over_Ob*m_p/r:>14.2f}")
    print(f"    => m_dark ~ {Odm_over_Ob*m_p/max(rs):.1f}-{Odm_over_Ob*m_p/min(rs):.1f} GeV (the ADM few-GeV scale)")

    print("\n  READ:")
    print("  * r is NOT arbitrary: group-theoretic sector weights give r = 0.5 (orders) to 0.9")
    print("    (#irreps), i.e. O(1) but < 1 -- the '12:24' heuristic made precise. The round 'r=1'")
    print("    of adm_abundance.py is the idealization; the group orders actually give r=0.5.")
    print("  * CONSEQUENCE: m_dark ~ 4-11 GeV (r<1 pushes it ABOVE the r=1 value of ~5 GeV). The")
    print("    canonical ADM few-GeV scale is robust to the weight measure (~factor-2-3 spread).")
    print("  * STILL a CONSISTENCY, not a prediction: (i) branching ~ static sector weight is an")
    print("    ASSUMPTION (the true r needs the freeze-out reconnection dynamics, uncomputed); (ii)")
    print("    m_dark is unpinned by the framework (D-type has NO Q_group -> no mass tower), so the")
    print("    observed ratio FIXES m_dark rather than the framework predicting Omega_dm/Omega_b.")
    print("  * What IS firm: r = O(1) is now computed (not assumed) from sector weights, and (q)'s")
    print("    r_forced=0 means it is NOT a forced per-baryon byproduct -- the O(1) is a shared-vacuum")
    print("    branching weight, not a one-per-baryon tie.")
