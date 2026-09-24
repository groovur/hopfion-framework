#!/usr/bin/env python3.11
"""
n_q^base tower calculation for the Q_H=3 (baryon/quark) sector.
Session 2026-09-21. See notes/quark_generation_e8_ribbon_twist.md sec.20.

Two parts:
  PART A -- DERIVE n_q^base via the CMB-energy-identity analogue for Q_H=3,
            exactly parallel to Paper XII (lepton, Q_group=10, n=20) and
            P17:prop:qh1_cmb (neutrino, Q_group=6). Verify arithmetic +
            recover m_e as calibration. Flag the saddle-blocked half.
  PART B -- EMPIRICAL tower diagnostic: extract the tower level each quark
            mass requires, quantify the isospin-dependent steepness, test
            candidate structures for the isospin split, give the u,c result.
"""
import numpy as np

PHI   = (1 + 5**0.5) / 2
LNPHI = np.log(PHI)
PI    = np.pi

# --- data (PDG central; MSbar 2 GeV for u,d,s; m_c(m_c), m_b(m_b), m_t pole) ---
m = {'u':2.16, 'd':4.67, 's':93.4, 'c':1270.0, 'b':4180.0, 't':172760.0}  # MeV
lep = {'e':0.51100, 'mu':105.658, 'tau':1776.86}                           # MeV
I3  = {'u':+0.5,'c':+0.5,'t':+0.5, 'd':-0.5,'s':-0.5,'b':-0.5}
gen = {'u':1,'d':1, 'c':2,'s':2, 't':3,'b':3}
m_q0 = 215.5   # constituent baryon scale (paper16), the tower anchor

T_CMB_K  = 2.7255
kB       = 8.617333e-5          # eV/K
T_CMB_eV = kB * T_CMB_K

def Lcond(Qgroup):
    """Sector condensate scale: rho_CMB / Lcond^4 = Qgroup  =>  Lcond=T_CMB(pi^2/(15 Q))^{1/4}."""
    return T_CMB_eV * (PI**2 / (15.0 * Qgroup))**0.25

print("="*74)
print("PART A -- DERIVATION of n_q^base via the CMB energy identity")
print("="*74)
rho_CMB_over_TCMB4 = PI**2/15   # Stefan-Boltzmann photon energy density coefficient
for name, Q in [('lepton QH=2', 10), ('neutrino QH=1', 6), ('BARYON QH=3', 3)]:
    L  = Lcond(Q)
    ratio = rho_CMB_over_TCMB4 * T_CMB_eV**4 / L**4   # = rho_CMB / Lcond^4
    print(f"  {name:14s}  Q_group={Q:2d} | Lcond=T_CMB(pi^2/{15*Q:3d})^(1/4)={L:.4e} eV"
          f" | rho_CMB/Lcond^4={ratio:.6f} (=Q_group)  n_base=2Q={2*Q}")
print("  => Baryon: Lcond^(baryon)=T_CMB(pi^2/45)^(1/4);  n_q^base = 2*Q_group^(baryon)=6")

# calibration: recover m_e from lepton chain  m_e = Lcond^(l) * phi^20 * e^{4pi-3/400}
L_lep = Lcond(10)
m_e_pred = L_lep * PHI**20 * np.exp(4*PI - 3/400)   # eV
print(f"\n  CALIBRATION  m_e = Lcond^(l) phi^20 e^(4pi-3/400) = {m_e_pred/1e6:.5f} MeV"
      f"  (PDG 0.51100)  ratio={m_e_pred/1e6/lep['e']:.4f}")
print("  => chain verified for the STABLE (saddle-bearing) lepton sector.")
print("  CAVEAT [O]: the profile-normalisation half  phi^9 sqrt(Ja[f*]J4[f*])=Q_group")
print("  needs a SADDLE profile f*. Q_H=3 has NO stable saddle (mass is algebraic),")
print("  so that half cannot be evaluated. n_q^base=6 is the ARITHMETIC (CMB) half of")
print("  the derivation, parallel to P17:prop:qh1_cmb; the saddle half is BLOCKED.")

print()
print("="*74)
print("PART B -- EMPIRICAL tower diagnostic (m = m_q0 * phi^{-2 n})")
print("="*74)
def nreq(mass):  # required tower level from m = m_q0 phi^{-2n}  (n>0 => below constituent)
    return -np.log(mass/m_q0)/(2*LNPHI)
n = {q: nreq(m[q]) for q in m}
print(f"  constituent anchor m_q0 = {m_q0} MeV,  n = -ln(m/m_q0)/(2 ln phi)")
print("   q   m(MeV)     n_req   I3   gen")
for q in ['u','d','s','c','b','t']:
    print(f"   {q}  {m[q]:9.2f}  {n[q]:+6.2f}  {I3[q]:+.1f}   {gen[q]}")

print("\n  per-generation tower increment  Delta n (g -> g+1):")
for iso,(g1,g2,g3) in [('UP',('u','c','t')), ('DOWN',('d','s','b'))]:
    d12, d23 = n[g2]-n[g1], n[g3]-n[g2]
    print(f"    {iso:4s}: {g1}->{g2} = {d12:+.2f} | {g2}->{g3} = {d23:+.2f} | mean {(n[g3]-n[g1])/2:+.2f}/gen")
print(f"  confined-only slope ratio (u->c)/(d->s) = {(n['c']-n['u'])/(n['s']-n['d']):.2f}"
      "   (top excluded: bare, incomplete-formation per sec.18)")

print("\n  ISOSPIN SPLIT  n_up - n_down  by generation (grows => up steeper):")
split = {g: n[{1:'u',2:'c',3:'t'}[g]] - n[{1:'d',2:'s',3:'b'}[g]] for g in (1,2,3)}
for g in (1,2,3):
    print(f"    gen{g}: {split[g]:+.2f}")

# --- candidate structures for the isospin split (test, don't assume) ---
print("\n  TEST candidate laws for the isospin split (CLAUDE.md 3: coincidence != law):")
# E_6 twist isospin split = (Tw_up - Tw_down) = {2,4,10} for gen {1,2,3}
tw_split = {1:2, 2:4, 3:10}
print("    E_6 twist split (Tw_up-Tw_down) = {2,4,10}; split/twist_split =",
      {g: round(split[g]/tw_split[g],3) for g in (1,2,3)}, "-> NOT constant")
# ratio of increments
print(f"    up-tower steepness / down = {abs(n['t']-n['u'])/abs(n['b']-n['d']):.2f} (all 3 gen)")

print()
print("="*74)
print("PART B.2 -- the u,c ratio, down-calibrated tower")
print("="*74)
# down sector anchored (strange=muon). Predict up ratios from the SAME tower spacing
# (isospin-blind null): if up followed the down per-gen increment, what would m_c/m_u be?
dn_down_12 = n['s'] - n['d']         # down gen1->gen2 increment
m_c_null   = m['u'] * PHI**(-2*dn_down_12)   # up-quark stepped by the DOWN increment
print(f"  measured m_c/m_u          = {m['c']/m['u']:.2f}")
print(f"  isospin-blind null (up steps by down's Dn={dn_down_12:+.2f}): m_c/m_u = "
      f"{m_c_null/m['u']:.2f}")
print(f"  => up sector climbs a factor {(m['c']/m['u'])/(m_c_null/m['u']):.2f} steeper than the")
print(f"     down-calibrated tower predicts -- the residual I_3-dependent steepness.")
print(f"  extra tower levels up needs vs down, gen1->2: Dn_up-Dn_down = "
      f"{(n['c']-n['u'])-(n['s']-n['d']):+.2f}")
