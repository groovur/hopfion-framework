"""Independent Z_eff prescriptions tested against tower-level placement.

Compares three ways of supplying Z_eff without inverting the measured
ionisation energy:
  (a) one-electron orbital formula IE = Z_eff^2 E_Ry / n_qn^2 with
      Slater's rules (the naive route),
  (b) Slater total-energy difference IE = E(ion) - E(atom) with
      E = -E_Ry sum_i (Z_eff,i/n*_i)^2 (Slater's actual prescription),
  (c) quantum-defect / effective principal number n*, IE = E_Ry (q+1)^2/n*^2.
Tower level: n(E) = n(E_Ry) - ln(E/E_Ry)/(2 ln phi), n(E_Ry) = -12.102.
"""
import math
PHI = (1 + 5 ** 0.5) / 2
ERY = 13.6057
NREF = -12.102
L = 2 * math.log(PHI)
nlev = lambda E: NREF - math.log(E / ERY) / L

GROUPS = ["1s", "2s2p", "3s3p", "3d", "4s4p", "4d", "4f", "5s5p", "5d", "6s6p"]
NSTAR = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7, 5: 4.0, 6: 4.2}
PRINC = {"1s": 1, "2s2p": 2, "3s3p": 3, "3d": 3, "4s4p": 4, "4d": 4,
         "4f": 4, "5s5p": 5, "5d": 5, "6s6p": 6}
CAP = {"1s": 2, "2s2p": 8, "3s3p": 8, "3d": 10, "4s4p": 8, "4d": 10,
       "4f": 14, "5s5p": 8, "5d": 10, "6s6p": 8}
FILL = ["1s", "2s2p", "3s3p", "4s4p", "3d", "4s4p", "5s5p", "4d",
        "6s6p", "4f", "5d"]

def config(nelec):
    """Occupancy per group by Aufbau order (madelung, simplified)."""
    order = ["1s", "2s2p", "3s3p", "4s4p", "3d", "5s5p", "4d", "6s6p",
             "4f", "5d"]
    occ = {g: 0 for g in GROUPS}
    left = nelec
    for g in order:
        room = CAP[g] - occ[g]
        take = min(room, left)
        occ[g] += take
        left -= take
        if left == 0:
            break
    return occ

def zeff(occ, g):
    n = PRINC[g]
    s = 0.0
    same = occ[g] - 1
    s += same * (0.30 if g == "1s" else 0.35)
    for h in GROUPS:
        if h == g or occ[h] == 0:
            continue
        m = PRINC[h]
        if g in ("3d", "4d", "4f", "5d"):
            s += occ[h] * 1.00 if (m < n or (m == n and h < g)) else 0.0
        else:
            if m == n - 1:
                s += occ[h] * 0.85
            elif m < n - 1:
                s += occ[h] * 1.00
            elif m == n and h != g:
                s += occ[h] * 0.35
    return sum(occ.values()) and (ZNUC[0] - s)

def total_energy(Z, nelec):
    if nelec <= 0:
        return 0.0
    occ = config(nelec)
    global ZNUC
    ZNUC = [Z]
    E = 0.0
    for g in GROUPS:
        if occ[g] == 0:
            continue
        ze = zeff(occ, g)
        E += occ[g] * (ze / NSTAR[PRINC[g]]) ** 2
    return -ERY * E

def ie_slater_total(Z, k):
    """k-th ionisation energy by Slater total-energy difference."""
    return total_energy(Z, Z - k) - total_energy(Z, Z - k + 1)

# experimental first IEs (NIST), eV
EXP = {1: ("H", 13.598), 2: ("He", 24.587), 3: ("Li", 5.392),
       4: ("Be", 9.323), 5: ("B", 8.298), 6: ("C", 11.260),
       7: ("N", 14.534), 8: ("O", 13.618), 9: ("F", 17.423),
       10: ("Ne", 21.565), 11: ("Na", 5.139), 12: ("Mg", 7.646),
       13: ("Al", 5.986), 14: ("Si", 8.152), 15: ("P", 10.487),
       16: ("S", 10.360), 17: ("Cl", 12.968), 18: ("Ar", 15.760),
       19: ("K", 4.341), 20: ("Ca", 6.113)}

print(f"{'el':>3} {'IE_exp':>8} {'IE_1e':>8} {'dn_1e':>7} {'IE_tot':>8} {'dn_tot':>7}")
d1, dt = [], []
for Z, (el, ie) in EXP.items():
    occ = config(Z)
    global ZNUC
    ZNUC = [Z]
    g = [h for h in GROUPS if occ[h] > 0][-1]
    ze = zeff(occ, g)
    ie_1e = ze ** 2 * ERY / PRINC[g] ** 2
    ie_tot = ie_slater_total(Z, 1)
    a = nlev(ie) - nlev(ie_1e)
    b = nlev(ie) - nlev(ie_tot) if ie_tot > 0 else float("nan")
    d1.append(abs(a)); dt.append(abs(b))
    print(f"{el:>3} {ie:8.3f} {ie_1e:8.2f} {a:7.3f} {ie_tot:8.2f} {b:7.3f}")
print(f"\nmean |dn| one-electron orbital : {sum(d1)/len(d1):.3f}")
print(f"mean |dn| Slater total-energy : {sum(dt)/len(dt):.3f}")
print(f"\nO 2p check: screening 3.45 -> Z_eff {8-3.45:.2f}, "
      f"IE_1e = {(8-3.45)**2*ERY/4:.1f} eV vs 13.618 exp "
      f"(factor {(8-3.45)**2*ERY/4/13.618:.1f})")
