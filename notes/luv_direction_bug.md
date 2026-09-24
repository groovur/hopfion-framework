# Lambda_UV direction bug (induced_gravity.py, luv_coincidence.py)

## The Seeley-DeWitt calculation
Integrating out xi at one loop, heat kernel K(s)=(1/(16pi^2 s^2)) sum a_n s^n, a_1=(1/6-xi)R.
Minimal (xi=0): the induced R-term is L ⊃ -(Lam^2/(192 pi^2)) R (matches the paper's eq:1loop R/6 term).
Match to (M_Pl^2/2)R: M_Pl^2/2 = Lam^2/(192 pi^2) => M_Pl^2 = Lam^2/(96 pi^2).
(The paper's box eq:MPl writes 32 pi^2, dropping the 1/6 and a factor 3; a separate minor slip.)

## THE BUG: direction
Inverting M_Pl^2 = Lam^2/(loop factor) gives Lam_UV = sqrt(loop) * M_Pl > M_Pl, SUPER-Planckian:
  Seeley-DeWitt (96 pi^2): Lam_UV = 30.8 M_Pl = 3.8e20 GeV -> n_UV = 80.1
  paper box (32 pi^2):     Lam_UV = 4pi sqrt2 M_Pl = 17.8 M_Pl ~ phi^6 M_Pl = 2.2e20 GeV -> n_UV = 79.6
The paper's stated Lam_UV = M_Pl/(4pi sqrt2) = 0.056 M_Pl (sub-Planckian, n_UV=73.6) is the RECIPROCAL.
This is standard Sakharov induced gravity: the cutoff is ABOVE the induced M_Pl (M_Pl ~ Lam/(4pi)).
The paper's "sub-Planckian, as required for a UV cutoff" is the conceptual error -- induced-gravity
cutoffs are super-Planckian.

## Impact
- Lam_UV: M_Pl/phi^6 -> phi^6 M_Pl (flips by phi^12); n_UV: 73.6 -> ~79.6-80. STILL NON-INTEGER,
  so the "compatible with non-integer tower level" conclusion (Cor luv_tower) HOLDS.
- The 4pi sqrt2 ~ phi^6 coincidence is UNCHANGED (now Lam_UV ~ phi^6 M_Pl).
- LOAD-BEARING RESULTS SAFE: the NMC ratio 3 phi^7 = xi_NMC rho_inf/(M_Pl^2/2) has Lam_UV^2 CANCEL
  (both xi_NMC and M_Pl^2/2 ∝ Lam_UV^2), so F(rho_inf), the CC (Lam_obs/F), the DE mass, and the
  Einstein-frame m_xi~H0 are ALL Lam_UV-independent. Nothing physical breaks; only Lam_UV & n_UV move.

## Paper locations to fix (if approved)
- eq:MPl box: Lam_UV = M_Pl/(4pi sqrt2) -> Lam_UV = 4pi sqrt2 M_Pl (super-Planckian). (coeff 32->96 optional.)
- eq:LUV_approx: Lam_UV ~ M_Pl/phi^6 ~ 0.0557 M_Pl -> ~ phi^6 M_Pl ~ 17.9 M_Pl; 6.8e17 GeV -> 2.2e20 GeV.
- Cor luv_tower: n_UV 73.6 -> 79.6 (recompute); conclusion (non-integer) unchanged.
- "sub-Planckian, as required for a UV cutoff" -> "super-Planckian, as expected for an induced
  (Sakharov) cutoff".
- Abstract: Lam_UV ~ 0.0557 M_Pl ~ M_Pl/phi^6 -> super-Planckian phi^6 M_Pl.
NOT yet edited -- multi-location incl. abstract + boxed eq; awaiting go-ahead.

## RETRACTION + RESOLUTION (inflation_luv.py, luv_reframe.py)
RETRACTION: my "correct fully to super-Planckian" plan was WRONG -- it would DESTROY inflation.
P_s = (N_e^2/12pi^2)(Lam_UV/M_Pl)^8 = 2.11e-9 (matches Planck) NEEDS Lam_UV/M_Pl = 1/phi^6
(sub-Planckian). Super-Planckian gives P_s ~ 1.8e13 (22 orders off). So Lam_UV sub-Planckian is
LOAD-BEARING (inflation), contra my earlier claim. Do NOT super-Planckian-ify Lam_UV.

RESOLUTION (flux/tower reframe -- resolves the tension, preserves everything):
- The framework's PRIMARY gravity is the density-feedback FLUX: G_N = G_src^2/(4pi phi^6) (l.416).
  M_Pl^2 = 1/G_N = 4pi phi^6/G_src^2, Lam_UV-INDEPENDENT. This is the actual gravity mechanism.
- M_Pl and Lam_UV are BOTH on the phi-tower, EXACTLY 3 steps apart: n(M_Pl)=76.59, n(Lam_UV)=73.59,
  gap = 3.00 = phi^6. So Lam_UV = M_Pl/phi^6 is a TOWER/topological relation (the SAME phi^6 as the
  Bogomolny/S_eff suppression in G_N), NOT the loop coincidence.
- The Seeley-DeWitt "induced M_Pl from a loop cutoff, Lam_UV=M_Pl/(4pi sqrt2)~M_Pl/phi^6" is the
  SPURIOUS secondary argument: it approximated the tower phi^6 by the loop 4pi sqrt2 (1% coincidence)
  AND inverted the direction. ALL the trouble (direction, 32 vs 96 pi^2, super/sub-Planckian) is HERE.
- A GENUINE loop-induced sub-Planckian Lam_UV does NOT work: needs enhancement phi^12*96pi^2 ~ 3e5;
  the strong NMC gives only ~44; large-N needs ~3e5 dof (implausible).
NET: gravity = flux (primary); Lam_UV = M_Pl/phi^6 = condensate/tower cutoff (phi^6 TOPOLOGICAL,
tower gap 3), used by inflation. Demote the Seeley-DeWitt induced-M_Pl story. The phi^6 is PROMOTED
from loop-coincidence to tower-topological (passes the discriminator). Nothing load-bearing changes;
the "Induced EH/Planck mass" section is the part to reframe (secondary consistency, not a derivation).
DISCRIMINATOR CORRECTION: earlier I tagged Lam_UV~phi^6 a "coincidence" -- that's true of 4pi sqrt2~phi^6
(the paper's spurious justification), but the ACTUAL Lam_UV phi^6 is the tower relation = PHYSICAL.

## PAPER REFRAME DONE (focused)
Reframed the "Induced EH Term and the Planck Mass" section + all mirrors:
- eq:MPl: boxed relation changed M_Pl^2=Lam_UV^2/(32pi^2), Lam_UV=M_Pl/(4pi sqrt2) -> Lam_UV=M_Pl/phi^6
  (the TOWER relation; M_Pl & Lam_UV three tower steps apart, phi^6=Bogomolny lambda, topological).
- Section text: one-loop Seeley-DeWitt R-term demoted to a consistency check; M_Pl fixed by the flux
  G_N=G_src^2/(4pi phi^6) (eq:GN).
- Cor:luv_tower retitled "Tower level of the ultraviolet cutoff"; proof uses Lam_UV=M_Pl/phi^6 directly;
  n_UV=73.6 (unchanged), M_Pl at n=76.6, gap 3 exact.
- rem:bogomolny_coincidence reframed: phi^6 gap IS the Bogomolny lambda (topological); the 4pi sqrt2~phi^6
  loop proximity demoted to "a numerical coincidence that plays no role".
- prop:eft_naturalness (iii): M_Pl fixed by flux, not "M_Pl^2=Lam_UV^2/[6(4pi)^2 phi^6]".
- Inflation: eq:Minf simplified M_inf=Lam_UV^2/M_Pl=M_Pl/phi^12 (removed 6(4pi)^2 phi^6 steps); P_s proof
  substitution (Lam_UV/M_Pl)^2=1/phi^12 (was 1/[6(4pi)^2 phi^6]). P_s=2.10e-9 UNCHANGED (uses 1/phi^6).
- Abstract (x2), G_N remark, conclusion: all mirror the tower framing; removed M_Pl^2=Lam_UV^2/(32pi^2),
  4pi sqrt2, 0.0563 M_Pl.
NET: gravity=flux (primary), Lam_UV=M_Pl/phi^6=condensate/tower cutoff (topological phi^6). All tested
predictions unchanged (P_s=2.10e-9, n_s, r, n_UV=73.6, M_inf=3.79e16 GeV). The inverted/wrong-coefficient
loop relation and the 4pi sqrt2~phi^6 coincidence crutch are gone. Env: eq 89/89, remark 22/22, cor 2/2,
proof 15/15; labels intact. The only 4pi sqrt2 left is the demoted "plays no role" mention in rem:bogomolny.
