# PAPER TODO (later) — sharpen dm_ic/bao/dm_status with the ordering-epoch closure

STATUS: the ordering epoch is now estimated (ordering_epoch.py, notes/ordering_epoch_closure.md):
T_ord ~ Lam_cond ~ T_CMB(today) ~ 1e-4-1e-3 eV, ~200-2600x below T_rec~0.26 eV -> medium DISORDERED at
recombination, orders at z~O(1)-few. Robust (cutoff-independent; set by Lam_cond from T_CMB).

## Edits to make (hedged as an order-of-magnitude estimate)
1. P7:sec:dm_ic / P7:rem:reheating_IC: change "whether the adiabatic baryon-sourced strain instead
   operates at recombination depends on the ordering epoch ... not fixed here" TO something like:
   "An order-of-magnitude estimate places the nematic ordering temperature at the condensate scale,
   T_ord ~ Lam_cond ~ T_CMB, so the medium orders only at z~O(1)-few, well after recombination: the
   director strain is a late-time (z<~few) galactic component and does not supply recombination-era
   cold matter. The recombination cold-dark-matter budget (Omega_DM, taken as an input) is a separate,
   as-yet-unidentified component."
2. P7:sec:bao: the 147-Mpc sentence — keep "no shift is firm", and replace "recovery ... is contingent
   on the ordering epoch, not settled here" with: "the standard r_s is set by the input recombination
   CDM; the late-time director strain (pressureless, comoving) shifts it not at all, but does not itself
   supply it."
3. P7:sec:dm_status: the "epoch at which the nematic medium orders -- which decides whether the sector
   contributes recombination-era CDM or only the late-time galactic component" -> resolve toward the
   late-time answer (z~few from Lam_cond), noting recombination CDM stays an input.
HEDGE: this is a heuristic T_ord ~ rho_cond xi_cond^3 estimate (O(1-10) factor), robust in DIRECTION
(T_ord << T_rec for any factor) but not a rigorous derivation. Present as "estimate places ... late",
not "proved". Pending author decision to fold in.
