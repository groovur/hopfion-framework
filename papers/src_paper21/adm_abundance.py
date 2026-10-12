#!/usr/bin/env python3.11
"""
adm_abundance.py -- abundance of the dark-Hopfion-wake CDM via QCD-epoch CO-FORMATION.
DM/Quark crossover (src_paper21). See notes/fabric_cmb_hypothesis_and_plan.md (i)-(o).

LOGIC: the residue gravitates with ~ the dark-Hopfion rest energy (c_dir->0 traps the
unwound energy as frozen strain, little radiated -> E_residue ~ f * m_dark, f~1). If
the dark Hopfion is a FORCED topological byproduct of baryon formation (like the PION,
one per baryon -- Paper XVIII Thm 5.3), its number tracks the baryon ASYMMETRY:
    n_dark / n_b = r   (r a TOPOLOGICAL ratio, not a free knob; r=1 = one-per-baryon).
Then the cosmic-coincidence ratio is pure arithmetic:
    Omega_DM / Omega_b = r * (E_residue / m_p) = r * f * (m_dark / m_p).
This is the ADM mechanism: the famous Omega_DM/Omega_b ~ 5 is EXPLAINED if r*f*m_dark ~
5 m_p -- i.e. an O(few)*QCD-scale object, which QCD-epoch co-formation gives for free.
NB: this EXPLAINS the coincidence + PINS the (previously free) dark mass; it does NOT
derive the baryon asymmetry eta itself (that is baryogenesis, unsolved in SM too).

SUPERSEDED FRAMING (notes (w),(ww)): this script was written under the TRANSIENT-WAKE picture
(notes (m)), where f = E_res/m_dark is the fraction of an UNWINDING knot's energy trapped as
frozen director strain. Update (w) adopted the STABLE-KNOT picture instead (Hopf-charge
conservation forbids unwinding; the wake also loses the ADM dark number). For the stable relic
the gravitating mass IS m_dark, so f = 1 identically and Omega_DM/Omega_b = r*(m_dark/m_p)
(Papers XXI, XXII). The f column below is kept only for the wake-era comparison. r is computed
from sector weights in src_paper22/cogenesis_branching.py (r = 0.5-0.9), no longer just assumed.
"""
m_p      = 0.9383                 # GeV
Ob_h2    = 0.02237                # Planck 2018
Odm_h2   = 0.1200
ratio_obs = Odm_h2 / Ob_h2
eta      = 6.1e-10                # baryon asymmetry n_b/n_gamma

print("="*70)
print("  DARK-HOPFION-WAKE ABUNDANCE via QCD co-formation (ADM framing)")
print("="*70)
print(f"  Observed:  Omega_DM h^2 = {Odm_h2},  Omega_b h^2 = {Ob_h2}")
print(f"             Omega_DM / Omega_b = {ratio_obs:.3f}   (the coincidence to explain)")
print(f"  Baryon asymmetry eta = n_b/n_gamma = {eta:.1e}")

print("\n  ADM prediction:  Omega_DM/Omega_b = r * f * (m_dark/m_p)")
print("  => required  r*f*m_dark = {:.2f} * m_p = {:.2f} GeV\n".format(ratio_obs, ratio_obs*m_p))

print(f"  {'r (dark/baryon #)':>18} {'f (trap frac)':>14} {'=> m_dark (GeV)':>16}")
for r in (1, 2, 3, 5):
    for f in (1.0, 0.5):
        m_dark = ratio_obs * m_p / (r * f)
        print(f"  {r:>18} {f:>14} {m_dark:>16.2f}")

print("\n  FRAMEWORK-NATIVE case: one forced byproduct per baryon (r=1, like the pion),")
print("  c_dir->0 traps ~all the unwound energy (f~1):")
m_canon = ratio_obs * m_p
print(f"     => m_dark ~ {m_canon:.2f} GeV  == the CANONICAL ADM mass (~5 m_p).")
print( "     This PINS the previously-unpinned dark mass (update (h)) AND explains the")
print( "     Omega_DM/Omega_b coincidence with NO fine tuning -- an O(few)*QCD-scale object.")
print( "     ~5 GeV + sterile (no SM coupling) => evades 50yr direct-detection nulls (a feature).")

print("\n  HONEST CAVEATS:")
print("   * r = topological byproduct ratio: ASSUMED O(1), not yet derived from the")
print("     confinement byproduct structure (the key uncomputed piece).")
print("   * f~1 assumes c_dir->0 traps nearly all unwound energy as frozen strain.")
print("   * eta (baryon asymmetry) itself is NOT explained (baryogenesis; SM-open too).")
print("   * symmetric (annihilating) dark Hopfions would instead give a thermal-relic")
print("     abundance DECOUPLED from Omega_b -> coincidence NOT explained. ADM (asymmetry")
print("     tracking the baryon asymmetry via co-formation) is what buys the coincidence.")
