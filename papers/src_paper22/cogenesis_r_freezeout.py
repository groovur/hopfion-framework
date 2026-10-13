#!/usr/bin/env python3.11
"""
cogenesis_r_freezeout.py -- the FREEZE-OUT (thermal) refinement of the co-genesis ratio
r = n_dark/n_b, the piece cogenesis_branching.py explicitly left open ("the true r needs the
freeze-out DYNAMICS ... NOT computed here"). DM/Quark crossover (src_paper22), notes (ww)->(zzz).

cogenesis_branching.py weighted r by STATIC group data: group ORDER |G| (0.5), #irreps (0.86),
level-1 primary count (1.33). Those are ZERO-temperature DOF counts. At a FREEZE-OUT the asymmetry
in each co-forming sector is seeded by the SAME CP bias (shared-2I delta_CP) acting on the THERMAL
degrees of freedom actually populated at T_QCD. For a WZW/CFT sector the thermal DOF count is the
CENTRAL CHARGE c (Stefan-Boltzmann for a CFT: rho ~ c T^4; g_* proportional to c). So the leading
DYNAMICAL (thermal-DOF-weighted) estimate is

    r = n_dark/n_b  ~  c_dark / c_baryon      (common per-DOF CP bias cancels),

with c the WZW central charge c = k dim(g)/(k + h^v). This is SHARPER than the raw group order
(which counts the regular-rep / T=0 DOF) and is the physically correct freeze-out weight.

RESULT: c_baryon = c((E_6)_1) = 6, c_dark = c(SO(10)_1) = 5 (2D_3) or c(SO(14)_1) = 7 (2D_5):
    r = 5/6 = 0.83  (2D_3)   or   7/6 = 1.17  (2D_5),
converging with the #irreps estimate (0.86) and tightening the O(1) band. Then
Omega_dm/Omega_b = r (m_dark/m_p) fixes m_dark ~ 6.0 GeV (2D_3) / 4.3 GeV (2D_5) -- the ADM few-GeV
scale, now from a thermal-DOF weight rather than an assumed O(1).

HONEST [O]: this still assumes EQUAL PER-DOF CP asymmetry (eps_dark = eps_baryon) -- the common
delta_CP acting identically per DOF. The genuinely-uncomputed residual is the channel-dependent
reconnection AMPLITUDE ratio eps_dark/eps_baryon (the trefoil-forming vs dihedral-defect-forming
rates). So r is now a THERMAL-DOF estimate (c-ratio), not yet a full reconnection-amplitude rate;
and m_dark remains framework-unpinned (no D-type Q_group) -> the observed ratio FIXES m_dark.
"""

def c_wzw(dim_g, h_dual, k=1):
    """WZW central charge c = k dim(g)/(k + h^v)."""
    from fractions import Fraction as Fr
    return Fr(k*dim_g, k + h_dual)

# (algebra, dim, dual Coxeter number h^v)
SECTORS = {
    '(E_6)_1  baryon 2T' : ('E6',   78, 12),
    'SO(10)_1 dark 2D_3' : ('so10', 45,  8),
    'SO(14)_1 dark 2D_5' : ('so14', 91, 12),
}
Odm_over_Ob = 0.1200/0.02237   # Planck 2018 = 5.363
m_p = 0.9383                   # GeV

if __name__ == '__main__':
    print("="*74)
    print("  CO-GENESIS r FROM FREEZE-OUT THERMAL DOF (central charge)")
    print("="*74)
    c = {}
    for name,(alg,dim,hv) in SECTORS.items():
        cc = c_wzw(dim,hv)
        c[name] = cc
        print(f"  {name:22}  c = {dim}/{1+hv} = {cc} = {float(cc):.4f}")
    c_bar = c['(E_6)_1  baryon 2T']
    print(f"\n  Freeze-out weight r = c_dark/c_baryon  (common per-DOF CP bias cancels):")
    print(f"    {'sector':>22} {'r = c_dark/c_bar':>18} {'m_dark = (ratio)m_p/r [GeV]':>30}")
    for name in ('SO(10)_1 dark 2D_3','SO(14)_1 dark 2D_5'):
        r = c[name]/c_bar
        m_dark = Odm_over_Ob*m_p/float(r)
        print(f"    {name:>22} {str(r)+' = '+format(float(r),'.3f'):>18} {m_dark:>30.2f}")

    print(f"\n  CONVERGENCE of the weight estimates for 2D_3 (all O(1)):")
    print(f"    group order |2D_3|/|2T| = 12/24   = 0.500   (T=0 regular-rep DOF)")
    print(f"    # group irreps          = 6/7     = 0.857")
    print(f"    central charge c        = 5/6     = 0.833   <-- freeze-out thermal DOF (this script)")
    print(f"    level-1 primaries       = 4/3     = 1.333   (topological/quantum-dim count)")
    print(f"    => the THERMAL-DOF weight (0.83) is the freeze-out-physical one; it agrees with the")
    print(f"       #irreps count and gives m_dark ~ 6 GeV (2D_3), the canonical ADM scale.")
    print(f"\n  STATUS: upgrades r from 'assumed O(1)' / 'static group order 0.5' to a THERMAL-DOF")
    print(f"  (central-charge) weight r = 5/6 ~ 0.83. REMAINING [O]: the per-channel CP-asymmetry")
    print(f"  (reconnection-amplitude) ratio eps_dark/eps_baryon, assumed 1 here; and m_dark stays")
    print(f"  framework-unpinned (no D-type Q_group) so the observed ratio FIXES it.")
    print("="*74)
