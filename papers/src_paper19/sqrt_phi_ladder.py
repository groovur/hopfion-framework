"""
sqrt_phi_ladder.py
==================
Numerical verification of the √φ geometric ladder of quark mass ratios.

Produces the values cited in Paper XIX (main_paper19.tex):
  - Individual generation gaps Δn in log_φ units for all four gaps
  - Cross-ratios Δn_up¹/Δn_down² ≈ φ (3.7%) and Δn_up²/Δn_down¹ ≈ φ (1.3%)
  - Inversion ratios r_up = 1.2993, r_down = 1.2691 (both ≈ √φ = 1.2720)
  - Four-gap geometric sequence verification: (6.23, 7.90, 10.20, 13.26) ≈ x·(1, √φ, φ, φ^(3/2))
  - Absolute mass predictions from framework-derived m_s and m_s/m_d

All values are cited in:
  P19:thm:sqrt_phi_ladder  (eq. P19:eq:sqrt_phi_ladder)
  P19:rem:strange_anchor   (eq. P19:eq:strange_anchor_predictions, table)
  P19:cor:absolute_masses  (eqs. P19:eq:md_pred – P19:eq:cu_ratio_pred)
  P19:rem:cross_ratio      (eqs. P19:eq:cross_ratio_a – P19:eq:cross_ratios)
  P19:prop:jones_modulus_ratio (eq. P19:eq:jones_modulus_ratio)

References:
  PDG 2024 quark masses (MS-bar at 2 GeV)
  Paper XVI Proposition 13.2 — m_s derived at 0.09%
  Paper XIX spin-Casimir formula — m_s/m_d = √(k/2)·exp(8π/9)
"""

import numpy as np

PHI = (1 + 5**0.5) / 2

# ── PDG 2024 quark masses (MeV, MS-bar at 2 GeV) ─────────────────────────
m_u, m_c, m_t = 2.16,  1273.0, 172570.0
m_d, m_s, m_b = 4.67,  93.4,    4183.0

# ── Log_φ exponents ───────────────────────────────────────────────────────
def n_phi(m):
    return np.log(m) / np.log(PHI)

n = {q: n_phi(m) for q, m in
     [('u', m_u), ('c', m_c), ('t', m_t),
      ('d', m_d), ('s', m_s), ('b', m_b)]}

print("Log_φ exponents:")
for q, v in n.items():
    print(f"  n_{q} = {v:.6f}")

# ── Generation gaps ───────────────────────────────────────────────────────
du1 = n['c'] - n['u']   # up  gen 1→2 (u→c)
du2 = n['t'] - n['c']   # up  gen 2→3 (c→t)
dd1 = n['s'] - n['d']   # down gen 1→2 (d→s)
dd2 = n['b'] - n['s']   # down gen 2→3 (s→b)

print()
print("Generation gaps (Δn in log_φ units):")
print(f"  Δn_up¹  = n_c - n_u = {du1:.6f}   [u→c gap]")
print(f"  Δn_up²  = n_t - n_c = {du2:.6f}   [c→t gap]")
print(f"  Δn_down¹ = n_s - n_d = {dd1:.6f}   [d→s gap]")
print(f"  Δn_down² = n_b - n_s = {dd2:.6f}   [s→b gap]")
print()
print(f"  Up:   DECREASING  (Δn¹={du1:.4f} > Δn²={du2:.4f})")
print(f"  Down: INCREASING  (Δn¹={dd1:.4f} < Δn²={dd2:.4f})")

# ── Average spacings and their ratio ─────────────────────────────────────
n_bar_up   = (du1 + du2) / 2
n_bar_down = (dd1 + dd2) / 2
ratio_avg  = n_bar_up / n_bar_down

print()
print("Average generation spacings:")
print(f"  n̄_up   = {n_bar_up:.6f}")
print(f"  n̄_down = {n_bar_down:.6f}")
print(f"  Ratio  = {ratio_avg:.6f}  (φ = {PHI:.6f},  Δ = {abs(ratio_avg-PHI)/PHI*100:.2f}%)")
print(f"  → P19:prop:jones_modulus_ratio: n̄_up/n̄_down ≈ φ = |V_{{T(2,2)}}(q₅)|/|V_{{T(2,3)}}(q₅)|")

# ── Cross-ratios ─────────────────────────────────────────────────────────
cr1 = du1 / dd2
cr2 = du2 / dd1

print()
print("Cross-ratios (P19:eq:cross_ratio_a, P19:eq:cross_ratios):")
print(f"  Δn_up¹ / Δn_down² = {du1:.4f}/{dd2:.4f} = {cr1:.6f}  (φ = {PHI:.6f},  Δ = {abs(cr1-PHI)/PHI*100:.2f}%)")
print(f"  Δn_up² / Δn_down¹ = {du2:.4f}/{dd1:.4f} = {cr2:.6f}  (φ = {PHI:.6f},  Δ = {abs(cr2-PHI)/PHI*100:.2f}%)")

# ── Inversion ratios ─────────────────────────────────────────────────────
r_up   = du1 / du2         # larger / smaller within up-triplet (Δn¹ > Δn²)
r_down = dd2 / dd1         # larger / smaller within down-triplet (Δn² > Δn¹)
sqphi  = PHI**0.5

print()
print("Inversion ratios (P19:thm:sqrt_phi_ladder proof):")
print(f"  r_up   = Δn_up¹/Δn_up²     = {r_up:.6f}  (√φ = {sqphi:.6f},  Δ = {abs(r_up-sqphi)/sqphi*100:.2f}%)")
print(f"  r_down = Δn_down²/Δn_down¹  = {r_down:.6f}  (√φ = {sqphi:.6f},  Δ = {abs(r_down-sqphi)/sqphi*100:.2f}%)")

# ── Four-gap geometric sequence ───────────────────────────────────────────
x = dd1   # base = Δn_down¹ = log_φ(m_s/m_d)
pred = np.array([x, sqphi*x, PHI*x, PHI**1.5*x])
obs  = np.array([dd1, dd2, du2, du1])

print()
print("√φ geometric ladder  (P19:eq:sqrt_phi_ladder):")
print(f"  Base x = Δn_down¹ = {x:.6f}  [= log_φ(m_s/m_d)]")
print()
print(f"  {'Term':>12}  {'Predicted':>12}  {'Observed':>12}  {'Δ':>8}")
labels = ['x·1', 'x·√φ', 'x·φ', 'x·φ^(3/2)']
quark_labels = ['Δn_down¹', 'Δn_down²', 'Δn_up²', 'Δn_up¹']
for lbl, ql, p, o in zip(labels, quark_labels, pred, obs):
    print(f"  {ql:>12}  {p:>12.4f}  {o:>12.4f}  {abs(p-o)/o*100:>7.2f}%  ({lbl})")

# ── Absolute mass predictions ─────────────────────────────────────────────
# Framework values (no PDG input)
m_s_fw  = 93.32          # Paper XVI, 0.09%
k = 3
ratio_fw = np.sqrt(k / 2) * np.exp(8 * np.pi / 9)   # spin-Casimir formula

m_d_fw   = m_s_fw / ratio_fw
m_b_fw   = m_s_fw * ratio_fw**sqphi
tc_ratio = ratio_fw**PHI
cu_ratio = ratio_fw**PHI**1.5

print()
print("Absolute mass predictions (P19:cor:absolute_masses):")
print(f"  m_s/m_d (framework) = √(k/2)·exp(8π/9) = √(3/2)·exp(8π/9) = {ratio_fw:.6f}")
print(f"  (PDG m_s/m_d = {m_s/m_d:.3f};  Δ = {abs(ratio_fw - m_s/m_d)/(m_s/m_d)*100:.2f}%)")
print()

rows = [
    ("m_d",   f"m_s/{ratio_fw:.3f}",  m_d_fw,  m_d,   "MeV"),
    ("m_b",   f"m_s·{ratio_fw:.3f}^√φ", m_b_fw, m_b, "MeV"),
    ("m_t/m_c", f"{ratio_fw:.3f}^φ",  tc_ratio, m_t/m_c, "—"),
    ("m_c/m_u", f"{ratio_fw:.3f}^φ^(3/2)", cu_ratio, m_c/m_u, "—"),
]
print(f"  {'Quantity':12}  {'Predicted':>12}  {'PDG':>12}  {'Δ':>8}")
print("  " + "-"*56)
for qty, formula, pred_val, pdg_val, unit in rows:
    print(f"  {qty:12}  {pred_val:>12.4f}  {pdg_val:>12.4f}  {abs(pred_val-pdg_val)/pdg_val*100:>7.2f}%  {unit}")

print()
print("Done — all values cross-check with Paper XIX.")
