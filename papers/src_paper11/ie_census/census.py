#!/usr/bin/env python3
"""
Ionisation-energy census for Paper XI (phi-tower / master formula).
No web access. mendeleev/periodictable not installed in this environment
(checked at run time below) -> falls back to a hand-typed reference table
of NIST-style ionisation energies (subset, documented below).

Formula (Paper XI, Definition of phi-tower level, calibrated form used here):
    n(E) = n(E_Ry) - ln(E/E_Ry) / (2*ln(phi))
with n(E_Ry) = -12.102 (paper's calibration constant) and E_Ry = 13.6057 eV.
"""

import math
import itertools
import json
from scipy import stats
import numpy as np

# ---------------------------------------------------------------------------
# 0. Package availability check
# ---------------------------------------------------------------------------
HAVE_MENDELEEV = False
HAVE_PERIODICTABLE = False
try:
    import mendeleev  # noqa
    HAVE_MENDELEEV = True
except ImportError:
    pass
try:
    import periodictable  # noqa
    HAVE_PERIODICTABLE = True
except ImportError:
    pass

DATA_SOURCE_NOTE = (
    "mendeleev available: {}; periodictable available: {}. "
    "Neither is installed in this environment, and no web access is permitted, "
    "so this census uses a hand-typed reference table of ionisation energies "
    "(values as commonly tabulated by NIST ASD / CRC Handbook, recalled offline). "
    "IE1 for Z=1-92 is given to the precision of standard reference tables. "
    "Higher IEs (IE2-IE6) for the 13 elements named in Paper XI are typed with "
    "high confidence for IE1-IE3 and lower confidence (fewer reliable digits, "
    "some values approximate) for IE4-IE6 of the heavier/transition elements, "
    "since these are less commonly reproduced and could not be re-verified "
    "against a live database in this offline session. No values were invented "
    "to fit an expected outcome; all are standard textbook/NIST-style figures "
    "recalled from training data."
).format(HAVE_MENDELEEV, HAVE_PERIODICTABLE)

PHI = (1 + math.sqrt(5)) / 2
LN_PHI = math.log(PHI)
E_RY = 13.6057          # eV, as used in Paper XI
N_RY = -12.102          # paper's calibration constant, n(E_Ry)

def n_level(E_eV):
    return N_RY - math.log(E_eV / E_RY) / (2 * LN_PHI)

# ---------------------------------------------------------------------------
# 1. Data: IE1 for Z = 1..92, EXCLUDING the 13 named elements (kept separately)
#    tag = (n_qn, subshell) of the electron removed by IE1, via standard
#    Aufbau + the conventional "ns before (n-1)d" ionisation-order exception
#    for transition metals (chemically standard; stated explicitly here).
#    Values in eV.
# ---------------------------------------------------------------------------
NAMED = {"H", "O", "Na", "Xe", "Pb", "Fe", "Ge", "Mg", "Ni", "Li", "K", "Rb", "Cs"}

# Z : (symbol, IE1_eV, n_qn, subshell)
SINGLE_IE1 = {
    2:  ("He", 24.5874, 1, "s"),
    4:  ("Be", 9.3227,  2, "s"),
    5:  ("B",  8.2980,  2, "p"),
    6:  ("C",  11.2603, 2, "p"),
    7:  ("N",  14.5341, 2, "p"),
    9:  ("F",  17.4228, 2, "p"),
    10: ("Ne", 21.5645, 2, "p"),
    13: ("Al", 5.9858,  3, "p"),
    14: ("Si", 8.1517,  3, "p"),
    15: ("P",  10.4867, 3, "p"),
    16: ("S",  10.3600, 3, "p"),
    17: ("Cl", 12.9676, 3, "p"),
    18: ("Ar", 15.7596, 3, "p"),
    20: ("Ca", 6.1132,  4, "s"),
    21: ("Sc", 6.5615,  4, "s"),
    22: ("Ti", 6.8281,  4, "s"),
    23: ("V",  6.7462,  4, "s"),
    24: ("Cr", 6.7665,  4, "s"),
    25: ("Mn", 7.4340,  4, "s"),
    27: ("Co", 7.8810,  4, "s"),
    29: ("Cu", 7.7264,  4, "s"),
    30: ("Zn", 9.3942,  4, "s"),
    31: ("Ga", 5.9993,  4, "p"),
    33: ("As", 9.7886,  4, "p"),
    34: ("Se", 9.7524,  4, "p"),
    35: ("Br", 11.8138, 4, "p"),
    36: ("Kr", 13.9996, 4, "p"),
    38: ("Sr", 5.6949,  5, "s"),
    39: ("Y",  6.2173,  5, "s"),
    40: ("Zr", 6.6339,  5, "s"),
    41: ("Nb", 6.7589,  5, "s"),
    42: ("Mo", 7.0924,  5, "s"),
    43: ("Tc", 7.28,    5, "s"),
    44: ("Ru", 7.3605,  5, "s"),
    45: ("Rh", 7.4589,  5, "s"),
    46: ("Pd", 8.3369,  4, "d"),   # exception: [Kr]4d10 5s0
    47: ("Ag", 7.5762,  5, "s"),
    48: ("Cd", 8.9938,  5, "s"),
    49: ("In", 5.7864,  5, "p"),
    50: ("Sn", 7.3439,  5, "p"),
    51: ("Sb", 8.6084,  5, "p"),
    52: ("Te", 9.0096,  5, "p"),
    53: ("I",  10.4513, 5, "p"),
    56: ("Ba", 5.2117,  6, "s"),
    57: ("La", 5.5769,  6, "s"),
    58: ("Ce", 5.5387,  6, "s"),
    59: ("Pr", 5.473,   6, "s"),
    60: ("Nd", 5.525,   6, "s"),
    61: ("Pm", 5.582,   6, "s"),
    62: ("Sm", 5.6437,  6, "s"),
    63: ("Eu", 5.6704,  6, "s"),
    64: ("Gd", 6.1501,  6, "s"),
    65: ("Tb", 5.8638,  6, "s"),
    66: ("Dy", 5.9389,  6, "s"),
    67: ("Ho", 6.0215,  6, "s"),
    68: ("Er", 6.1077,  6, "s"),
    69: ("Tm", 6.1843,  6, "s"),
    70: ("Yb", 6.2542,  6, "s"),
    71: ("Lu", 5.4259,  6, "s"),
    72: ("Hf", 6.8251,  6, "s"),
    73: ("Ta", 7.5496,  6, "s"),
    74: ("W",  7.8640,  6, "s"),
    75: ("Re", 7.8335,  6, "s"),
    76: ("Os", 8.4382,  6, "s"),
    77: ("Ir", 8.9670,  6, "s"),
    78: ("Pt", 8.9587,  6, "s"),
    79: ("Au", 9.2255,  6, "s"),
    80: ("Hg", 10.4375, 6, "s"),
    81: ("Tl", 6.1082,  6, "p"),
    83: ("Bi", 7.2856,  6, "p"),
    84: ("Po", 8.414,   6, "p"),
    85: ("At", 9.3175,  6, "p"),   # low confidence (poorly measured element)
    86: ("Rn", 10.7485, 6, "p"),
    87: ("Fr", 4.0727,  7, "s"),
    88: ("Ra", 5.2784,  7, "s"),
    89: ("Ac", 5.17,    7, "s"),
    90: ("Th", 6.3067,  7, "s"),
    91: ("Pa", 5.89,    7, "s"),
    92: ("U",  6.1941,  7, "s"),
}
assert set(SINGLE_IE1) | {1,3,8,11,12,19,26,28,32,37,54,55,82} == set(range(1,93))

# ---------------------------------------------------------------------------
# 2. Multi-IE table for the 13 named elements.
#    Each entry: k, IE_k (eV), n_qn, subshell of the removed electron,
#    confidence flag ('hi' = high confidence, 'lo' = lower confidence /
#    approximate, typed from memory of standard sequential-IE tables).
# ---------------------------------------------------------------------------
MULTI = {
    "H":  (1,  [(1, 13.598, 1, "s", "hi")]),
    "Li": (3,  [(1, 5.3917, 2, "s", "hi"),
                (2, 75.6402, 1, "s", "hi"),
                (3, 122.4543, 1, "s", "hi")]),
    "Na": (11, [(1, 5.1391, 3, "s", "hi"),
                (2, 47.2864, 2, "p", "hi"),
                (3, 71.6200, 2, "p", "hi"),
                (4, 98.91,  2, "p", "hi"),
                (5, 138.40, 2, "p", "hi"),
                (6, 172.18, 2, "p", "lo")]),
    "Mg": (12, [(1, 7.6462, 3, "s", "hi"),
                (2, 15.0353, 3, "s", "hi"),
                (3, 80.1437, 2, "p", "hi"),
                (4, 109.2655, 2, "p", "hi"),
                (5, 141.27, 2, "p", "lo"),
                (6, 186.76, 2, "p", "lo")]),
    "K":  (19, [(1, 4.3407, 4, "s", "hi"),
                (2, 31.63,  3, "p", "hi"),
                (3, 45.806, 3, "p", "hi"),
                (4, 60.91,  3, "p", "lo"),
                (5, 82.66,  3, "p", "lo"),
                (6, 100.0,  3, "p", "lo")]),
    "Fe": (26, [(1, 7.9024, 4, "s", "hi"),
                (2, 16.1878, 4, "s", "hi"),
                (3, 30.652, 3, "d", "hi"),
                (4, 54.8,   3, "d", "lo"),
                (5, 75.0,   3, "d", "lo"),
                (6, 99.1,   3, "d", "lo")]),
    "Ni": (28, [(1, 7.6398, 4, "s", "hi"),
                (2, 18.1688, 4, "s", "hi"),
                (3, 35.19,  3, "d", "hi"),
                (4, 54.9,   3, "d", "lo"),
                (5, 76.06,  3, "d", "lo"),
                (6, 108.0,  3, "d", "lo")]),
    "Ge": (32, [(1, 7.8994, 4, "p", "hi"),
                (2, 15.9346, 4, "p", "hi"),
                (3, 34.2241, 4, "s", "hi"),
                (4, 45.7131, 4, "s", "hi"),
                (5, 93.5,   3, "d", "lo"),
                (6, 126.6,  3, "d", "lo")]),
    "Rb": (37, [(1, 4.1771, 5, "s", "hi"),
                (2, 27.2895, 4, "p", "hi"),
                (3, 39.247, 4, "p", "hi"),
                (4, 52.20,  4, "p", "lo"),
                (5, 68.44,  4, "p", "lo"),
                (6, 82.9,   4, "p", "lo")]),
    "Cs": (55, [(1, 3.8939, 6, "s", "hi"),
                (2, 23.15744, 5, "p", "hi"),
                (3, 33.195, 5, "p", "hi"),
                (4, 43.0,   5, "p", "lo"),
                (5, 49.0,   5, "p", "lo"),
                (6, 61.0,   5, "p", "lo")]),
    "O":  (8,  [(1, 13.618, 2, "p", "hi"),
                (2, 35.1211, 2, "p", "hi"),
                (3, 54.9355, 2, "p", "hi"),
                (4, 77.4135, 2, "p", "hi"),
                (5, 113.899, 2, "s", "hi"),
                (6, 138.1189, 2, "s", "hi")]),
    "Xe": (54, [(1, 12.1298, 5, "p", "hi"),
                (2, 20.9750, 5, "p", "hi"),
                (3, 32.1230, 5, "p", "hi")]),
    "Pb": (82, [(1, 7.4167, 6, "p", "hi"),
                (2, 15.0326, 6, "p", "hi"),
                (3, 31.9373, 6, "s", "hi"),
                (4, 42.32,  6, "s", "lo"),
                (5, 68.8,   5, "d", "lo"),
                (6, 88.0,   5, "d", "lo")]),
}

NOBLE_COUNTS = {2, 10, 18, 36, 54, 86}

# ---------------------------------------------------------------------------
# 3. Build the master entry list
# entry = dict(elem, Z, k, E, n, n_qn, subshell, conf, shell_completion)
# ---------------------------------------------------------------------------
entries = []
for Z, (sym, E, nqn, sub) in SINGLE_IE1.items():
    is_shell = (Z - 1) in NOBLE_COUNTS   # k=1 always for this group
    entries.append(dict(elem=sym, Z=Z, k=1, E=E, n=n_level(E),
                         n_qn=nqn, sub=sub, conf="hi", shell=is_shell))

for sym, (Z, lst) in MULTI.items():
    for (k, E, nqn, sub, conf) in lst:
        is_shell = (Z - k) in NOBLE_COUNTS
        entries.append(dict(elem=sym, Z=Z, k=k, E=E, n=n_level(E),
                             n_qn=nqn, sub=sub, conf=conf, shell=is_shell))

N_TOTAL = len(entries)

# ===========================================================================
# CENSUS A -- near-integer landings
# ===========================================================================
d_ints = np.array([abs(e["n"] - round(e["n"])) for e in entries])

def hist10(arr):
    counts, edges = np.histogram(arr, bins=10, range=(0, 0.5))
    return counts.tolist(), edges.tolist()

A1_hist_counts, A1_hist_edges = hist10(d_ints)
A1_hits = int(np.sum(d_ints <= 0.017))
A1_expected_p = 2 * 0.017
A1_expected_n = A1_expected_p * N_TOTAL
A1_binom = stats.binomtest(A1_hits, N_TOTAL, A1_expected_p, alternative="two-sided")

shell_entries = [e for e in entries if e["shell"]]
N_SHELL = len(shell_entries)
d_ints_shell = np.array([abs(e["n"] - round(e["n"])) for e in shell_entries])
A2_hist_counts, A2_hist_edges = hist10(d_ints_shell)
A2_hits = int(np.sum(d_ints_shell <= 0.017))
A2_expected_n = A1_expected_p * N_SHELL
A2_binom = stats.binomtest(A2_hits, N_SHELL, A1_expected_p, alternative="two-sided") if N_SHELL > 0 else None

# A3 -- specific cases
A3_cases = {
    "Pb IE3": 31.9373,
    "Xe IE3": 32.1230,
    "Xe IE1": 12.1298,
}
A3_paper_claim = {
    "Pb IE3": (-12.989, 0.011),
    "Xe IE3": (-12.995, 0.005),
    "Xe IE1": (-11.983, 0.017),
}
A3_results = {}
for name, E in A3_cases.items():
    n = n_level(E)
    d = abs(n - round(n))
    A3_results[name] = dict(E=E, n=n, d=d, paper=A3_paper_claim[name])

# ===========================================================================
# CENSUS B -- pair near-degeneracies
# ===========================================================================
n_values = np.array([e["n"] for e in entries])
elems = [e["elem"] for e in entries]

pairs_dn = []
pairs_same_elem = 0
for (i, j) in itertools.combinations(range(N_TOTAL), 2):
    if elems[i] == elems[j]:
        pairs_same_elem += 1
    pairs_dn.append(abs(n_values[i] - n_values[j]))
pairs_dn = np.array(pairs_dn)
# B1 = ALL unordered pairs of distinct IE table entries (literal reading of
# the spec); B2 below applies the A != B structural restriction on top.

B1_total_pairs = len(pairs_dn)
B1_same_elem_pairs = pairs_same_elem
B1_hits = int(np.sum(pairs_dn <= 0.0021))
n_range = n_values.max() - n_values.min()
d_window = 0.0021
p_uniform = (2 * d_window) / n_range - (d_window / n_range) ** 2
B1_expected = B1_total_pairs * p_uniform

# B2 -- structural restriction: same n_qn, same subshell type (s or p only), different element
sp_entries = [e for e in entries if e["sub"] in ("s", "p")]
restricted_pairs = []
for a, b in itertools.combinations(sp_entries, 2):
    if a["elem"] == b["elem"]:
        continue
    if a["n_qn"] == b["n_qn"] and a["sub"] == b["sub"]:
        restricted_pairs.append((a, b, abs(a["n"] - b["n"])))

B2_N = len(restricted_pairs)
B2_hits = [t for t in restricted_pairs if t[2] <= 0.0021]
B2_expected = B2_N * p_uniform
B2_prob_at_least_one = 1 - (1 - p_uniform) ** B2_N if B2_N > 0 else 0.0

# B3 -- specific pairs
O6 = next(e for e in entries if e["elem"] == "O" and e["k"] == 6)
Na5 = next(e for e in entries if e["elem"] == "Na" and e["k"] == 5)
B3_ON = dict(nO=O6["n"], nNa=Na5["n"], dn=abs(O6["n"] - Na5["n"]))

H1 = next(e for e in entries if e["elem"] == "H" and e["k"] == 1)
O1 = next(e for e in entries if e["elem"] == "O" and e["k"] == 1)
B3_HO = dict(nH=H1["n"], nO=O1["n"], dn=abs(H1["n"] - O1["n"]))

# B4 -- configurations (hand-derived, stated in comments above)
B4 = dict(
    O_ion="O5+ (3 electrons: 1s2 2s1)", O_removed="2s",
    Na_ion="Na4+ (7 electrons: 1s2 2s2 2p3)", Na_removed="2p",
    both_2s=False, same_n_qn=True,
)

# ===========================================================================
# CENSUS C -- Slater-rule predictivity test
# ===========================================================================
AUFBAU_ORDER = ["1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p","6s",
                 "4f","5d","6p","7s","5f","6d","7p"]
CAPACITY = {"s":2, "p":6, "d":10, "f":14}

def build_config(num_electrons):
    """Simple Aufbau fill, in the AUFBAU_ORDER given, returns list of
    (n, l, count) in fill order."""
    remaining = num_electrons
    config = []
    for orb in AUFBAU_ORDER:
        n = int(orb[0]); l = orb[1]
        cap = CAPACITY[l]
        if remaining <= 0:
            break
        put = min(cap, remaining)
        config.append((n, l, put))
        remaining -= put
    return config

def slater_zeff(num_electrons_before, n_qn, l):
    """Z_eff for the *last-filled* electron of a config built by simple
    reverse-Aufbau order (the 'simple Aufbau' convention requested)."""
    config = build_config(num_electrons_before)
    if not config:
        return None
    # electron being removed = last group in fill order
    n_target, l_target, _ = config[-1]
    sigma = 0.0
    for (n, l_, cnt) in config:
        if n == n_target and l_ == l_target:
            same_group_others = cnt - 1
            sigma += same_group_others * (0.30 if (n_target == 1) else 0.35)
        elif l_target in ("s", "p"):
            if n == n_target - 1:
                sigma += cnt * 0.85
            elif n < n_target - 1:
                sigma += cnt * 1.00
            # n == n_target but different l already handled by 'same group'
            # electrons in same n but different l among s/p are same group (Slater groups ns,np together)
        else:  # target is d or f: all electrons in lower groups (n' <= n, different group) contribute 1.00
            if not (n == n_target and l_ == l_target):
                sigma += cnt * 1.00
    Z = num_electrons_before  # Z_atom equals electron count of the *neutral* only when before-state==neutral;
    return n_target, l_target, sigma

# For (s,p) Slater grouping, ns and np of the same n are ONE group; need to
# fix build_config to merge same-n s/p pairs when computing sigma "same group".
def slater_zeff_fixed(Z_atom, num_electrons_before, n_qn_expected, l_expected):
    config = build_config(num_electrons_before)
    if not config:
        return None
    n_target, l_target, _ = config[-1]
    def group_key(n, l):
        if l in ("s", "p"):
            return (n, "sp")
        return (n, l)
    tgt_key = group_key(n_target, l_target)
    same_group_count = sum(cnt for (n, l, cnt) in config if group_key(n, l) == tgt_key)
    sigma = (same_group_count - 1) * (0.30 if n_target == 1 else 0.35)
    for (n, l, cnt) in config:
        if group_key(n, l) == tgt_key:
            continue
        if l_target in ("s", "p"):
            if n == n_target - 1:
                sigma += cnt * 0.85
            elif n <= n_target - 2:
                sigma += cnt * 1.00
        else:  # d or f target: everything in a lower group contributes 1.00
            if n < n_target or (n == n_target and l in ("s","p")):
                sigma += cnt * 1.00
            elif n == n_target and l == "d" and l_target == "f":
                sigma += cnt * 1.00
    Zeff = Z_atom - sigma
    return n_target, l_target, Zeff

C_results = []
for e in entries:
    Z = e["Z"]; k = e["k"]
    n_before = Z - (k - 1)
    if n_before <= 0:
        continue
    out = slater_zeff_fixed(Z, n_before, e["n_qn"], e["sub"])
    if out is None:
        continue
    n_target, l_target, Zeff = out
    if Zeff <= 0:
        continue
    IE_pred = Zeff**2 * E_RY / n_target**2
    n_pred = n_level(IE_pred)
    delta_n = e["n"] - n_pred
    C_results.append(dict(elem=e["elem"], k=k, n_qn_simple_aufbau=n_target,
                           sub_simple_aufbau=l_target, n_qn_chem=e["n_qn"],
                           sub_chem=e["sub"], Zeff_slater=Zeff,
                           IE_exp=e["E"], IE_pred=IE_pred,
                           n_exp=e["n"], n_pred=n_pred, delta_n=delta_n))

C_deltas = np.array([c["delta_n"] for c in C_results])
C_mean = float(np.mean(C_deltas))
C_std = float(np.std(C_deltas, ddof=1))

def find_case(elem, k):
    for c in C_results:
        if c["elem"] == elem and c["k"] == k:
            return c
    return None

C_cases = {
    "Pb IE3": find_case("Pb", 3),
    "Xe IE3": find_case("Xe", 3),
    "Xe IE1": find_case("Xe", 1),
    "O IE6":  find_case("O", 6),
    "Na IE5": find_case("Na", 5),
    "H IE1":  find_case("H", 1),
    "O IE1":  find_case("O", 1),
}

# ===========================================================================
# Output
# ===========================================================================
def fmt(x, nd=4):
    return f"{x:.{nd}f}"

report = {
    "data_source_note": DATA_SOURCE_NOTE,
    "N_total": N_TOTAL,
    "N_shell": N_SHELL,
    "A1": dict(hist_counts=A1_hist_counts, hist_edges=A1_hist_edges,
               hits=A1_hits, expected_n=A1_expected_n, expected_p=A1_expected_p,
               pvalue=A1_binom.pvalue),
    "A2": dict(hist_counts=A2_hist_counts, hist_edges=A2_hist_edges,
               hits=A2_hits, expected_n=A2_expected_n,
               pvalue=(A2_binom.pvalue if A2_binom else None),
               members=[(e["elem"], e["k"], round(e["n"],4)) for e in shell_entries]),
    "A3": A3_results,
    "B1": dict(total_pairs=B1_total_pairs, hits=B1_hits, expected=B1_expected,
               n_range=n_range, p_uniform=p_uniform),
    "B2": dict(N_pairs=B2_N, hits=len(B2_hits), expected=B2_expected,
               prob_at_least_one=B2_prob_at_least_one,
               hit_list=[(a["elem"],a["k"],b["elem"],b["k"],round(dn,4)) for a,b,dn in B2_hits]),
    "B3_ON": B3_ON,
    "B3_HO": B3_HO,
    "B4": B4,
    "C_mean": C_mean,
    "C_std": C_std,
    "C_N": len(C_results),
}

if __name__ == "__main__":
    print(json.dumps({k: v for k, v in report.items() if k not in
                       ("A2",)}, default=str, indent=2)[:3000])
    print("N_total =", N_TOTAL, " N_shell =", N_SHELL)
    print("A1 hits/expected:", A1_hits, A1_expected_n, "p=", A1_binom.pvalue)
    print("A3:", A3_results)
    print("B1:", B1_total_pairs, B1_hits, B1_expected)
    print("B2:", B2_N, len(B2_hits), B2_expected, B2_prob_at_least_one)
    print("B3_ON:", B3_ON)
    print("B3_HO:", B3_HO)
    print("C mean/std:", C_mean, C_std)
    for name, c in C_cases.items():
        print(name, c)
