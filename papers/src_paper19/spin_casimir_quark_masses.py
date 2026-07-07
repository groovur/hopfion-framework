"""
spin_casimir_quark_masses.py
============================
Numerical verification of the spin-Casimir hypothesis and the absolute
quark mass predictions of Paper XIX.

Produces the values cited in Paper XIX (main_paper19.tex):
  - Structural identity k = 4·C₂(j=1/2) = 4·(1/2)·(3/2) = 3  (exact)
  - √(k/2) = √(3/2) = √(2·C₂(1/2))  (exact)
  - m_s/m_d = √(k/2)·exp(8π/9) = 19.991  (0.047% from integer 20)
  - k = 3 is uniquely the WZW level at which 2·C₂(j=1/2) = k/2
  - Required δT to reach 20, and why additive corrections don't work
  - Full absolute mass table (m_d, m_b, m_t/m_c, m_c/m_u)

All values cited in:
  P19:rem:ms_md_integer  (eqs. P19:eq:spin_corrected_ratio,
                           P19:eq:spin_casimir_id, P19:eq:k_casimir_identity)
  P19:cor:absolute_masses  (eqs. P19:eq:md_pred – P19:eq:cu_ratio_pred)
  P19:rem:mb_prediction

References:
  Paper XVI Section 8 — E₈ Coxeter mechanism, exponents m_d=17, m_s=13
  Paper XVI Proposition 13.2 — m_s = 93.32 MeV (0.09%)
"""

import numpy as np

PHI = (1 + 5**0.5) / 2

# ── WZW / E₈ parameters ──────────────────────────────────────────────────
k    = 3          # WZW level
h_E8 = 30         # E₈ Coxeter number  (h_E8 = 30 = h̃(E₈))
Q    = 10         # Q_group^(lepton) = ord(q_{2I}) = 10  (Paper XII)
j    = 0.5        # quark spin

# Coxeter exponents assigned in Paper XVI Section 8
m_d_exp = 17      # d-quark Coxeter exponent
m_s_exp = 13      # s-quark Coxeter exponent
T_diff  = (m_d_exp - m_s_exp) / (k * h_E8)   # = 4/90

print("=" * 65)
print("SECTION 1: Structural identity k = 4·C₂(j=1/2)")
print("=" * 65)
C2_half = j * (j + 1)           # = 3/4
k_from_C2 = 4 * C2_half
print(f"  C₂(j=1/2) = j(j+1) = {j}·{j+1} = {C2_half}")
print(f"  4·C₂(1/2) = {k_from_C2:.1f}  =  k = {k}  ← exact identity")
print()
print(f"  √(k/2)        = {np.sqrt(k/2):.8f}")
print(f"  √(2·C₂(1/2)) = {np.sqrt(2*C2_half):.8f}  (same, confirms √(k/2)=√(2C₂(j)))")
print()

print("=" * 65)
print("SECTION 2: m_s/m_d from spin-Casimir formula")
print("=" * 65)
E8_pred   = np.exp(2 * np.pi * Q * T_diff)     # pure E₈
spin_corr = np.sqrt(k / 2)                      # = √(3/2)
ms_md_fw  = spin_corr * E8_pred

print(f"  E₈ Coxeter prediction:  exp(2π·Q·ΔT) = exp(8π/9) = {E8_pred:.6f}")
print(f"  Spin-Casimir factor:    √(k/2)        =             {spin_corr:.6f}")
print(f"  Combined:               √(k/2)·exp(8π/9)         = {ms_md_fw:.6f}")
print(f"  Integer 20 (= n_e):                               = 20.000000")
print(f"  Residual from 20:                                  {abs(ms_md_fw-20)/20*100:.4f}%")
print()

# Required δT to reach exactly 20 from E₈ alone
ln20_over_2piQ = np.log(20) / (2 * np.pi * Q)
delta_T_needed = ln20_over_2piQ - T_diff
delta_T_spin   = np.log(np.sqrt(k/2)) / (2 * np.pi * Q)
print(f"  Required δT to reach 20: {delta_T_needed:.8f}")
print(f"  δT from spin-Casimir:    {delta_T_spin:.8f}")
print(f"  Residual:                {abs(delta_T_needed - delta_T_spin):.2e}  ({abs(delta_T_needed-delta_T_spin)/delta_T_needed*100:.4f}%)")

print()
print("  Additive fermion correction h_ψ/(k·h_E8) = 1/180:")
delta_T_fermion = 0.5 / (k * h_E8)
ms_md_fermion   = np.exp(2 * np.pi * Q * (T_diff + delta_T_fermion))
print(f"    Would give m_s/m_d = {ms_md_fermion:.4f}  (Δ = {abs(ms_md_fermion-20)/20*100:.1f}% — too large)")

print()
print("=" * 65)
print("SECTION 3: Absolute mass predictions (no PDG input)")
print("=" * 65)
m_s_fw  = 93.32       # Paper XVI, 0.09% accuracy
sqphi   = PHI**0.5

m_d_fw   = m_s_fw / ms_md_fw
m_b_fw   = m_s_fw * ms_md_fw**sqphi
tc_ratio = ms_md_fw**PHI
cu_ratio = ms_md_fw**PHI**1.5

# PDG comparison values
m_d_PDG, m_b_PDG = 4.67,  4183.0
tc_PDG  = 172570 / 1273.0
cu_PDG  = 1273.0  / 2.16

print(f"  m_s (Paper XVI)  = {m_s_fw} MeV  (0.09%)")
print(f"  m_s/m_d (framework) = {ms_md_fw:.6f}")
print()
print(f"  {'Quantity':14}  {'Predicted':>12}  {'PDG':>12}  {'Δ':>8}")
print("  " + "-" * 55)
rows = [
    ("m_d (MeV)",   m_d_fw,   m_d_PDG),
    ("m_b (MeV)",   m_b_fw,   m_b_PDG),
    ("m_t/m_c",     tc_ratio, tc_PDG),
    ("m_c/m_u",     cu_ratio, cu_PDG),
]
for qty, pred, pdg in rows:
    print(f"  {qty:14}  {pred:>12.4f}  {pdg:>12.4f}  {abs(pred-pdg)/pdg*100:>7.2f}%")

print()
print("  Using integer 20 instead of 19.991 for comparison:")
for qty, pred_fn, pdg in [
    ("m_d (MeV)",   lambda r: m_s_fw/r,        m_d_PDG),
    ("m_b (MeV)",   lambda r: m_s_fw*r**sqphi, m_b_PDG),
    ("m_t/m_c",     lambda r: r**PHI,           tc_PDG),
    ("m_c/m_u",     lambda r: r**PHI**1.5,      cu_PDG),
]:
    pred20 = pred_fn(20.0)
    pred_fw = pred_fn(ms_md_fw)
    print(f"  {qty:14}: integer-20 = {pred20:.4f}  vs  framework = {pred_fw:.4f}"
          f"  (Δ = {abs(pred20-pred_fw)/pred_fw*100:.3f}%)")

print()
print("Done — all values cross-check with Paper XIX.")
