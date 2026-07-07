"""
Strange-quark mass: framework prediction run down to PDG comparison scale.

Formula (Paper V, colour-level-shift route):
    m_q^(g) = v_EW * exp(-2*pi*Q*T_g) * phi^{-2*delta_n^(C)}

For the strange quark: g=2 (down-type), T_2=1/8, delta_n^(C)=1/2 (down-type
shift is phi^{-3} overall, per Paper V eq:down_quark_mass).

This predicts m_s at the condensate scale Lambda_UV ~ 6.8e17 GeV. To compare
with the PDG value (quoted at mu=2 GeV in the MS-bar scheme), the mass must
be run down through proper multi-loop, multi-threshold QCD RGEs. This script
performs that running at one-loop order, with flavor-threshold matching at
each quark mass (nf=6 -> 5 -> 4 -> 3), using alpha_s(M_Z)=0.1180 (PDG) as
the anchor for alpha_s's own running.

Result: m_s(2 GeV) = 93.32 MeV, vs PDG 93.40 MeV -- a 0.09% match.
"""

import math

# ── Framework constants ──────────────────────────────────────────────────
phi = (1 + math.sqrt(5)) / 2
Q = 10                      # lepton-sector closing integer (Paper I/IV)
v_EW = 246.22                # GeV, Paper XIV (0.015% residual vs measured 246.24 GeV... )
T2 = 1.0 / 8                  # T_g for generation g=2 (k_g=2), Paper I eq:T

# ── Step 1: framework prediction at the condensate scale Lambda_UV ────────
LUV = 6.8e17                 # GeV, Paper VI

m_s_LUV = v_EW * math.exp(-2 * math.pi * Q * T2) * phi**(-3)   # GeV, down-type shift
print(f"m_s(Lambda_UV) = {m_s_LUV*1000:.4f} MeV")

# ── Step 2: one-loop alpha_s running with flavor thresholds ───────────────
alpha_s_MZ = 0.1180           # PDG 2024
M_Z = 91.1876                 # GeV
m_t_thr = 172.57               # GeV, pole mass (top threshold)
m_b_thr = 4.18                 # GeV, MS-bar m_b(m_b) (bottom threshold)
m_c_thr = 1.27                 # GeV, MS-bar m_c(m_c) (charm threshold)
mu_target = 2.0                # GeV, standard light-quark comparison scale

def run_alpha_s_1loop(alpha0, mu0, mu, nf):
    """One-loop running of alpha_s from mu0 to mu with nf active flavors."""
    b0 = (33 - 2 * nf) / 3
    denom = 1 + alpha0 / (2 * math.pi) * b0 * math.log(mu / mu0)
    return alpha0 / denom

# Anchor alpha_s at LUV by running UP from M_Z through the top threshold
alpha_s_mt = run_alpha_s_1loop(alpha_s_MZ, M_Z, m_t_thr, nf=5)
alpha_s_LUV = run_alpha_s_1loop(alpha_s_mt, m_t_thr, LUV, nf=6)
print(f"alpha_s(Lambda_UV) = {alpha_s_LUV:.6f}  [anchored via alpha_s(M_Z)=0.1180]")

# Run alpha_s back DOWN from LUV through each flavor threshold
alpha_s_at_mt  = run_alpha_s_1loop(alpha_s_LUV,  LUV,     m_t_thr, nf=6)
alpha_s_at_mb  = run_alpha_s_1loop(alpha_s_at_mt, m_t_thr, m_b_thr, nf=5)
alpha_s_at_mc  = run_alpha_s_1loop(alpha_s_at_mb, m_b_thr, m_c_thr, nf=4)
alpha_s_at_2GeV = run_alpha_s_1loop(alpha_s_at_mc, m_c_thr, mu_target, nf=3)

print(f"alpha_s(m_t) = {alpha_s_at_mt:.6f}")
print(f"alpha_s(m_b) = {alpha_s_at_mb:.6f}")
print(f"alpha_s(m_c) = {alpha_s_at_mc:.6f}")
print(f"alpha_s(2 GeV) = {alpha_s_at_2GeV:.6f}  [real-world value ~0.30]")

# ── Step 3: run the quark mass itself through the same thresholds ─────────
gamma0 = 8.0   # universal one-loop QCD mass anomalous dimension

def run_mass_1loop(m0, alpha0, alpha_f, nf):
    """One-loop running of a quark mass between two scales with nf flavors."""
    b0 = (33 - 2 * nf) / 3
    exponent = gamma0 / (2 * b0)
    return m0 * (alpha_f / alpha0) ** exponent

m = m_s_LUV
m = run_mass_1loop(m, alpha_s_LUV,   alpha_s_at_mt, nf=6)   # LUV   -> m_t
m = run_mass_1loop(m, alpha_s_at_mt, alpha_s_at_mb, nf=5)   # m_t   -> m_b
m = run_mass_1loop(m, alpha_s_at_mb, alpha_s_at_mc, nf=4)   # m_b   -> m_c
m = run_mass_1loop(m, alpha_s_at_mc, alpha_s_at_2GeV, nf=3) # m_c   -> 2 GeV
m_s_2GeV = m

m_s_PDG = 93.4e-3   # GeV, PDG 2024 MS-bar at 2 GeV

print()
print(f"m_s(2 GeV) predicted = {m_s_2GeV*1000:.4f} MeV")
print(f"m_s(2 GeV) measured  = {m_s_PDG*1000:.2f} MeV  (PDG 2024)")
print(f"Ratio predicted/measured = {m_s_2GeV/m_s_PDG:.4f}  "
      f"({abs(1-m_s_2GeV/m_s_PDG)*100:.2f}% deviation)")
