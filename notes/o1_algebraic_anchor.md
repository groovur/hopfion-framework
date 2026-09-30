# Paper XX O1: algebraic absolute anchor of the E_6 up-tower (2026-09-28)

O1 = the single ABSOLUTE normalisation M_0 of the confined u,c tower
m_g = M_0 exp(-n_g pi/4). The tower SPACING is parameter-free (m_u/m_c=e^{-2pi},
P20:prop:up_ratio). Only M_0 was open, previously gotten from the FRAGILE subtraction
    completion = m_p - sigma(3R0+R) = 938.27 - 924.8 = 13.5 MeV  (-> /3 = 4.5 MeV/quark),
a difference of two ~930 MeV numbers (a 1% error in sigma(3R0+R) -> ~70% in the 13.5).

## ALGEBRAIC ANCHOR [N, this session] -- no subtraction
**completion energy = Lambda_cond^baryon * e^{8pi} = 13.2 MeV** (total), M_0 = /3 = 4.41 MeV/quark.
This is the constituent scale WITHOUT the phi^6 tower factor:
  constituent = Lambda_cond^baryon * phi^{2Q_gr} * e^{8pi} = 237 MeV  (P20:prop:constituent, Q_gr=3 -> phi^6)
  completion  = Lambda_cond^baryon *            e^{8pi} =  13.2 MeV  = constituent / phi^6.
Lambda_cond^baryon = T_CMB (pi^2/45)^{1/4} = 1.607e-4 keV (baryon-sector condensate scale).
`formation_balance_uc.py` + the anchor test reproduce: direct 13.22 vs subtraction 13.49 (2.0% apart);
data anchor M_0(u)=2.16/e^{-pi/4}=4.74, M_0(c)=1270/e^{7pi/4}=5.20 MeV/quark (geo-mean 4.96); direct
4.41 matches to ~10% (within the e^{-2pi} u/c shape residual).

## PHYSICAL STORY (current vs constituent)
The e^{8pi} Chern-Simons spoke sets the bare formation scale; the phi^{2Q_gr}=phi^6 WZW/Bogomolny tower
factor is the TOPOLOGICAL DRESSING that lifts it to the jammed constituent mass. The current-quark
anchor (the de-jamming/localisation energy released at hadronisation) is the UN-dressed spoke scale;
constituent = anchor x phi^6. This is the framework's current-vs-constituent distinction: constituent/
completion = phi^6. [Interpretive; the numeric agreement is the evidence.]

## CROSS-CHECKS / STATUS (CLAUDE.md sec.3, sec.6)
- Direct (Lambda e^{8pi}) vs subtraction (m_p - sigma(3R0+R)): agree 2%. INDEPENDENCE: direct inputs
  {T_CMB, 8pi, phi^6} vs subtraction {m_p, sigma=646.5/9, R0, r0} are disjoint IFF sigma's natural-units
  ->MeV conversion is independent of Lambda_cond. The notes flag that conversion as "TBD via the
  master-formula profile" (quark_generation sec.19, l.329/365), so partial circularity is NOT fully
  excluded -- state as corroboration-modulo-that, not a clean independent check.
- Direct vs data anchor (4.7-5.2/quark): ~10%, = the expected shared-jam / e^{-2pi} residual.
- [N] lead: numerically solid (2%/10%), physically coherent; the phi^6 current/constituent split is
  interpretive. Tag accordingly in the paper.

## CONSEQUENCE FOR O1 / THE JAMMED SOLVER
The ABSOLUTE anchor is now ALGEBRAIC (Lambda e^{8pi} = constituent/phi^6) -- it does NOT need the jammed
Q_H=3 field configuration, correcting the old O1 text ("computing the completion energy ... requires the
jammed configuration"). Only the ~5-10% per-generation SHAPE residual (the shared-jam deviation, the u/c
and d-s splits at common scale) plausibly needs the jammed config; the anchor itself does not. Consistent
with the standing conclusion (do not build the jammed field solver; masses/anchor are algebraic).

## PRE-EXISTING FLAG -- RESOLVED (stale in the script)
`formation_balance_uc.py`'s printed warning ("E_6-native ratio e^{-2pi} != committed master
exp(-2pi Q_gr T_g), reconcile") is STALE. The notes already settled it (quark_generation secs 23-24):
the up sector is **E_6-NATIVE** -- phase pi/4 = pi Q_gr/(k h(E_6)) = 3pi/(1*12), grounded in the quark's
own numbers (Q_gr=3, k=1=(E_6)_1, h(E_6)=12) -- and the bare exp(-2pi Q_gr T_g) is Paper XVI's LEADING
form (P16:eq:mq_leading) which the E_6-native refinement SUPERSEDES for the up sector (with n_q=72 T_g-
2 m_E6 carrying the Coxeter-exponent structure). They are NOT meant to be equal; there is no contradiction
to reconcile -- the reconciliation IS the E_8-bridge(down)/E_6-native(up) split (l.1188). Up-ratio checks
out: n_u-n_c=72(5/24-1/8)-2(7-8)=8 -> m_u/m_c=e^{-8 pi/4}=e^{-2pi} (0.95x data). The script comment
predates the settlement; ignore it. (Still genuinely open, but SEPARATE and about O2 not O1: the isospin
doublet SPLITTING / gen-1 sign-flip, which the naive E_6-native form does not give, sec.24 RESULT(3).)
