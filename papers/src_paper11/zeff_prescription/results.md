# Independent Z_eff prescriptions vs tower-level placement

Question: can Z_eff be supplied independently of the measured
ionisation energy, so that the master formula of Paper XI predicts a
tower level rather than restating it?

Tower level: n(E) = n(E_Ry) - ln(E/E_Ry)/(2 ln phi), n(E_Ry) = -12.102.
Experimental first ionisation energies (NIST), Z = 1..20.

## Result

| prescription | mean abs dn (Z<=20) |
|---|---|
| one-electron orbital, IE = Z_eff^2 E_Ry/n_qn^2, Slater Z_eff | 0.946 |
| Slater total-energy difference, IE = E(ion) - E(atom) | 0.191 |

The one-electron orbital formula is not Slater's prescription: Slater's
rules were constructed for total energies
E = -E_Ry sum_i (Z_eff,i / n*_i)^2, and the ionisation energy is the
difference between the ion and the neutral atom. Used that way the
same rules are five times more accurate.

Oxygen, the case discussed in Paper XI Section 6:
Slater screening for the 2p electron is
sigma = 0.85*2 (1s) + 0.35*5 (other n=2) = 3.45, so Z_eff = 4.55
(2s and 2p are one Slater group, screening one another at 0.35).
One-electron orbital: IE = 4.55^2 E_Ry/4 = 70.4 eV, a factor 5.2 too
large. Total-energy difference: IE = 14.17 eV, 4% high, |dn| = 0.041.

## Bearing on the multi-layer spiral claims

Even the correct Slater prescription leaves a scatter of about 0.19
tower levels across Z <= 20 (worst cases: Al 0.61, Si 0.43, S 0.39).
That is roughly ten times the near-integer window used in the
shell-completion placements (0.017) and about ninety times the
near-degeneracy separations discussed (0.0021). No prescription
available here can certify structure at those scales; the placements
remain descriptive, and an independent Z_eff accurate to better than
0.01 tower levels would be required to make them predictive.

Run: python3.11 slater.py

## Authoritative Z_eff values (mendeleev package, Clementi et al. data)

Oxygen, from `mendeleev.element('O').zeff(n, o, method=...)`:

| orbital | Clementi | Slater | sigma (Clementi) | sigma (Slater) |
|---|---|---|---|---|
| 1s | 7.6579 | 7.70 | 0.342 | 0.30 |
| 2s | 4.4916 | 4.55 | 3.508 | 3.45 |
| 2p | 4.4532 | 4.55 | 3.547 | 3.45 |

Two points follow.

1. The 2s and 2p effective charges differ by only 0.038, so the two
   orbitals see essentially the same effective nuclear charge. This
   validates Slater's grouping of 2s and 2p into a single screening
   group, and rules out treating 2s as an inner shell screening 2p at
   0.85. Taking the 1s contribution at the standard 0.85 each, the
   Clementi value implies an intra-shell screening of 0.369 per n=2
   electron, against Slater's uniform 0.35.

2. Slater and Clementi agree on Z_eff(O, 2p) to 2% (4.55 vs 4.4532).
   Neither is near 2. The value Z_eff = 2 is not an orbital screening
   constant; it is the measured ionisation energy re-expressed through
   the one-electron formula IE = Z_eff^2 E_Ry / n_qn^2.

Screening ceiling: a single electron cannot screen more than one unit
of nuclear charge. Reaching sigma_tot = 6.00 with the 1s pair at 0.85
and 2p-2p at 0.35 would require each 2s electron to screen 1.625 units;
even granting the 1s pair the maximum 1.0 each, the requirement is
1.475. No admissible assignment of screening constants produces
sigma_tot = 6.00.
