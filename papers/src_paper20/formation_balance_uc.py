#!/usr/bin/env python3.11
"""
formation_balance_uc.py -- pin the ABSOLUTE scale of the confined u,c doublet from
the T(2,2)->T(2,3) formation balance (Paper XX O2, Paper XVIII O6).

The E_6-native RATIO is fixed: m_u/m_c=e^{-2pi} (P20:prop:up_ratio). What remains is
the single absolute anchor M_0 of the confined tower m_g = M_0 exp(-n_g pi/4),
n_g = 72 T_g - 2 m_g^{E6}. This script (a) fits M_0 from u,c; (b) compares it to the
framework's DERIVED formation scales (the QGP-jam constituent and the completion/
localisation energy, sec.19); (c) sanity-checks the ratio formula against the
committed master exp(-2pi Q_gr T_g).

Established formation scales (sec.19, all framework-derived):
  * sigma_tot = 646.5/9 = 71.83 MeV/unit (FIXED, not free)
  * R = R0+r0 = 3 + sqrt2/phi = 3.874 (jam/un-threading scale)
  * jammed proton  m_p^jam = sigma*(3R0+R) = 924.8 MeV (98.6% of m_p, parameter-free)
  * completion/localisation deficit = m_p - m_p^jam = 13.5 MeV (paid at hadronisation;
    sec.19: 'light-quark/current-mass order', 2m_u+m_d~9 MeV)
  * constituent scale = Lcond^baryon phi^6 e^{8pi} ~ 234 MeV (8pi CS spoke, sec.48)
"""
import numpy as np
PHI=(1+5**0.5)/2; PI=np.pi
kB=8.617333e-5; TCMB=2.7255*kB   # eV

# ── derived formation scales ────────────────────────────────────────────────
sigma=646.5/9                    # MeV/unit
R0=3.0; r0=np.sqrt(2)/PHI; R=R0+r0
m_p_jam=sigma*(3*R0+R)           # jammed proton
m_p=938.272
deficit=m_p-m_p_jam              # completion/localisation, total (3 quarks)
Lcond_b=TCMB*(PI**2/(15*3))**0.25    # baryon condensate scale (eV)
m_const=Lcond_b*PHI**6*np.exp(8*PI)/1e6   # MeV
print("="*74,"\n  DERIVED formation scales (sec.19, sec.48)\n","="*74,sep="")
print(f"  sigma={sigma:.2f} MeV/u   R=R0+r0={R:.3f}   jammed proton sigma(3R0+R)={m_p_jam:.1f} MeV")
print(f"  completion deficit m_p - jam = {deficit:.1f} MeV (3 quarks) -> per quark {deficit/3:.2f} MeV")
print(f"  constituent Lcond^baryon phi^6 e^8pi = {m_const:.1f} MeV   (jammed/hadron scale)")

# ── fit the confined-tower anchor M_0 from u,c ──────────────────────────────
Tg={1:5/24,2:1/8,3:3/40}; mE6={1:7,2:8,3:11}
n={g:72*Tg[g]-2*mE6[g] for g in (1,2,3)}         # bare levels: n_u=1, n_c=-7, n_t=-16.6
m_obs={1:2.16,2:1270.}                            # u,c current (MeV)
M0={g: m_obs[g]/np.exp(-n[g]*PI/4) for g in (1,2)}
M0_gm=np.sqrt(M0[1]*M0[2])
print("\n"+"="*74,"\n  CONFINED-TOWER ANCHOR M_0  (m_g = M_0 exp(-n_g pi/4))\n","="*74,sep="")
print(f"  bare levels n: u={n[1]:.1f} c={n[2]:.1f}   (n_u-n_c={n[1]-n[2]:.1f} -> m_u/m_c=e^-2pi)")
print(f"  M_0 from u = {M0[1]:.2f} MeV;  from c = {M0[2]:.2f} MeV;  geo-mean = {M0_gm:.2f} MeV")
print(f"  (u,c agree to {abs(M0[1]/M0[2]-1)*100:.0f}% = the e^-2pi residual)")

print("\n"+"="*74,"\n  WHAT IS M_0?  compare to the derived formation scales\n","="*74,sep="")
for lbl,val in [('completion/quark (m_p-jam)/3',deficit/3),
                ('constituent phi^6 e^8pi',m_const),
                ('jammed proton /3',m_p_jam/3),
                ('electron m_e',0.511),('down m_d(2GeV)',4.67)]:
    print(f"  M_0={M0_gm:.2f} vs {lbl:30s}={val:8.2f} MeV  -> ratio {M0_gm/val:.3f}")
print(f"\n  => M_0 ~ {M0_gm:.1f} MeV sits at the COMPLETION/localisation scale (proton deficit/3),")
print(f"     NOT the constituent ({m_const:.0f} MeV, the jammed/hadron scale). The confined current-")
print(f"     mass tower anchors on the de-jamming completion energy paid at hadronisation.")

# ── sanity: does my ratio formula, not the committed master, give m_u/m_c? ───
print("\n"+"="*74,"\n  RATIO-FORMULA consistency check (CLAUDE.md sec.1)\n","="*74,sep="")
r_mine=np.exp(-(n[1]-n[2])*PI/4)
r_master=np.exp(-2*PI*3*(Tg[1]-Tg[2]))           # exp(-2pi Q_gr T_g), Q_gr^baryon=3
r_data_raw=2.16/1270.; r_data_cs=2.16/1097.
print(f"  m_u/m_c:  mine (72 T_g-2m_E6, pi/4) = {r_mine:.3e} = e^-2pi")
print(f"            committed master exp(-2pi Q_gr T_g), Q=3 = {r_master:.3e}")
print(f"            data raw {r_data_raw:.3e}   common-scale {r_data_cs:.3e}")
print(f"  => my E_6-native formula matches data ({r_mine/r_data_cs:.2f}x); the bare master")
print(f"     exp(-2pi Q_gr T_g) does NOT ({r_master/r_data_cs:.0f}x) -- the E_6-native form carries")
print(f"     the extra 6h(E_6) coefficient + Coxeter-exponent term. FLAG: reconcile the E_6-native")
print(f"     generation formula with the master before over-building on the absolute anchor.")
