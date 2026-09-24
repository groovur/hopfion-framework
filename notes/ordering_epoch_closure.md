# Ordering epoch CLOSED (ordering_epoch.py) -- director-strain DM is late-time

TRIGGER: the 147-Mpc question hinges on whether the nematic medium orders before recombination
(so the adiabatic baryon-sourced strain is present as recombination-era cold matter). Earlier we
thought this was cutoff-set/unknown. It is NOT -- the ordering temperature is set by the CONDENSATE
scale Lam_cond (fixed by T_CMB), not the unresolved UV cutoff.

## The estimate
Nematic ordering temperature ~ energy per coherence cell:
  T_ord ~ rho_cond * xi_cond^3.  rho_cond = 10 Lam_cond^4, xi_cond = 1/Lam_cond (natural units)
  => T_ord ~ 10 Lam_cond.  With the O(0.1-1) band stiffness and (pi/2) factor: T_ord ~ (0.1-10) Lam_cond.
Lam_cond = T_CMB(pi^2/150)^{1/4} = 1.19e-4 eV ~ T_CMB(today). So T_ord ~ 1e-4 to 1e-3 eV.
Recombination: T_rec ~ 0.26 eV (z~1100). => T_ord is 200-2600x BELOW T_rec, for ANY band factor.

## Result (robust to the O(1-10) uncertainty)
The nematic medium is thoroughly DISORDERED at recombination and orders only when the universe cools
to ~Lam_cond, i.e. z_order ~ O(1)-few:
  T_ord=Lam_cond -> z_order~-0.5 ; T_ord=10 Lam_cond -> z_order~4.
Empirical anchor: galactic DM (RAR/flat curves) is observed out to z~few, requiring T_ord > T(z~2)~7e-4 eV
=> the favoured value is T_ord ~ 10 Lam_cond, z_order ~ 4 (galaxy-formation epoch). Ordered for z<~4,
disordered before. Either way z_order << z_rec=1100.

## PHYSICAL POINT (why it's robust and cutoff-independent)
Lam_cond ~ T_CMB(today) by construction (Lam_cond = T_CMB(pi^2/150)^{1/4} ~ T_CMB/2). A condensate orders
when the universe cools to ~its own scale = ~now. At recombination the universe was ~1100x hotter than the
condensate scale, so thoroughly disordered. The director-strain DM is therefore INTRINSICALLY late-time
(tied to the condensate scale being ~today's temperature) -- not a recombination-era object. No dependence
on the UV cutoff; only on Lam_cond (from T_CMB).

## CLOSURE of the 147-Mpc / WP-line
- Ordering epoch: z ~ O(1)-few (LATE), robustly. NOT recombination.
- The adiabatic baryon-sourced strain is ABSENT at recombination -> the director strain does NOT supply
  recombination-era CDM. It is a LATE-TIME (z<~4) galactic component.
- The recombination CDM that sets the standard 147 Mpc and the CMB peaks (Omega_DM~0.26, already an INPUT)
  is a SEPARATE, unidentified component -- consistent with the notes' "early-universe CDM still open".
- 147 Mpc: set by the input recombination CDM (standard), NOT by the director strain; the strain (late,
  pressureless, comoving) shifts nothing. So "no shift" holds and the standard scale is set by the input,
  not recovered by this sector.
NET: the contingency is resolved NEGATIVELY -- the director strain is late-time; recombination CDM is a
separate input/gap. Paper's dm_ic/bao "not settled/contingent" can be sharpened to "late-time; ordering
epoch z~few from the condensate scale; recombination CDM is the separate Omega_DM input". Pending author call.
