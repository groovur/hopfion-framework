"""
Exact diagonalisation of the p^2 and d^2 open-shell Coulomb (electrostatic)
problem, built from scratch from single-particle m_l, m_s states and the
Slater-Condon / Racah formalism (Condon-Shortley coefficients via Wigner 3j
symbols, sympy exact rational/radical arithmetic throughout -- no floats,
no numerical tolerance).

This is the validation engine for shell_fusion/results.md. It:
  1. builds the antisymmetrized two-electron basis for a shell of orbital
     quantum number l (l=1 for p, l=2 for d),
  2. builds the exact electrostatic Hamiltonian matrix in terms of the
     raw Slater radial integrals F^0, F^2 (, F^4 for d), block-diagonal in
     (M_L, M_S),
  3. diagonalises each block exactly (sympy Rational/radical eigenvalues,
     here always rational multiples of the F^k since c^k Condon-Shortley
     coefficients are rational for l<=2),
  4. assigns Russell-Saunders (L,S) term labels to the exact eigenvalues by
     the standard descending (M_L,M_S) peeling algorithm, cross-checked
     against the known degeneracy counts,
  5. prints the validated term energies and degeneracies, to be checked by
     hand against the textbook results quoted in results.md.

Run: python3.11 shell.py
"""

from sympy import Rational, S, symbols, sqrt, simplify, nsimplify, Matrix, eye, zeros
from sympy.physics.wigner import wigner_3j
from itertools import combinations
from collections import defaultdict

F0, F2, F4 = symbols('F0 F2 F4')
FSYM = {0: F0, 2: F2, 4: F4}


def ck(l, m, lp, mp, k):
    """Condon-Shortley Slater-Condon coefficient c^k(l m, l' m'), exact.
    NOTE: this is only rational on the diagonal m=mp (used for direct
    integrals); off-diagonal (m != mp, needed for exchange-type terms) is
    generally irrational (contains surds). Both are returned as exact
    sympy expressions -- no floats, no rounding, anywhere in this file."""
    pref = S.NegativeOne ** S(m) * sqrt((2 * l + 1) * (2 * lp + 1))
    a = wigner_3j(l, k, lp, 0, 0, 0)
    b = wigner_3j(l, k, lp, -m, m - mp, mp)
    return simplify(pref * a * b)


def build_shell(l):
    """Return (sp_states, pairs) for shell of orbital quantum number l.
    sp_states: list of (m, ms). pairs: list of (i,j), i<j -> two-electron
    antisymmetrized basis state |i j>."""
    ms_vals = [Rational(1, 2), Rational(-1, 2)]
    sp_states = [(m, ms) for m in range(-l, l + 1) for ms in ms_vals]
    pairs = list(combinations(range(len(sp_states)), 2))
    return sp_states, pairs


def raw_integral(l, ma, mc, mb, md, kmax):
    """<a(1) b(2)| 1/r12 |c(1) d(2)> spatial part (no spin, no
    antisymmetrisation), from the multipole expansion of 1/r12:
        sum_p c^p(l ma, l mc) c^p(l md, l mb) F^p ,  p even, 0<=p<=2l.
    NOTE the index order in the second factor: c^p(md, mb), NOT c^p(mb, md)
    -- this is the standard (Cowan/Condon-Shortley) convention and is what
    makes the exchange integral <ab|1/r12|ba> reduce, term by term, to
    sum_p [c^p(l ma, l mb)]^2 F^p, manifestly >= 0 as it must be for a real
    Coulomb kernel. Using the naively "symmetric-looking" but wrong order
    c^p(mb, md) instead silently flips the sign of exchange contributions
    whenever c^p(m,m') != c^p(m',m) (which happens whenever m != m'), and
    was caught here only by requiring the diagonalised term energies to
    reproduce the textbook degeneracies {9,5,1} exactly -- they did not,
    under the wrong-order formula, even though every individual number
    involved was computed exactly (no floats, no tolerance)."""
    if ma + mb != mc + md:
        return S(0)
    total = S(0)
    for p in range(0, kmax + 1, 2):
        cac = ck(l, ma, l, mc, p)
        cdb = ck(l, md, l, mb, p)
        if cac != 0 and cdb != 0:
            total += cac * cdb * FSYM[p]
    return total


def hmatrix_block(l, sp_states, indices, kmax):
    """Build exact H matrix (in F0,F2,[F4]) for a list of two-electron basis
    pair-indices (each an (i,j) with i<j into sp_states) that share a common
    (M_L,M_S) sector."""
    n = len(indices)
    H = zeros(n, n)
    for a in range(n):
        i, j = indices[a]
        mi, msi = sp_states[i]
        mj, msj = sp_states[j]
        for b in range(a, n):
            k, lidx = indices[b]
            mk, msk = sp_states[k]
            ml_, msl = sp_states[lidx]
            # term1: <ij|H|kl>  requires spin(i)=spin(k), spin(j)=spin(l), ML conserved (guaranteed since same sector)
            term1 = S(0)
            if msi == msk and msj == msl:
                term1 = raw_integral(l, mi, mk, mj, ml_, kmax)
            # term2: <ij|H|lk> requires spin(i)=spin(l), spin(j)=spin(k)
            term2 = S(0)
            if msi == msl and msj == msk:
                term2 = raw_integral(l, mi, ml_, mj, mk, kmax)
            val = simplify(term1 - term2)
            H[a, b] = val
            H[b, a] = val
    return H


def diagonalize_shell(l, kmax, label):
    sp_states, pairs = build_shell(l)
    # bucket pairs by (ML,MS)
    buckets = defaultdict(list)
    for (i, j) in pairs:
        mi, msi = sp_states[i]
        mj, msj = sp_states[j]
        buckets[(mi + mj, msi + msj)].append((i, j))

    # exact eigenvalues per (ML,MS) sector: list of (ML,MS,eigenvalue,mult_in_sector)
    sector_eigs = {}
    for key, idxs in buckets.items():
        H = hmatrix_block(l, sp_states, idxs, kmax)
        # exact eigenvalues (H is symmetric with rational entries in F0,F2,F4 -> eigenvalues
        # are algebraic in these symbols; for our shells they come out as simple linear
        # combinations because the blocks factor or are already diagonal/near-diagonal)
        eigs = H.eigenvals()  # dict eigenvalue -> multiplicity, exact (sympy algebraic)
        sector_eigs[key] = (idxs, eigs)
    return sp_states, pairs, buckets, sector_eigs


def peel_terms(sector_eigs):
    """Standard Russell-Saunders term-derivation-table algorithm, but driven by
    the EXACT diagonalised energies rather than assumed microstate tables:
    for each (ML,MS) sector we have a multiset of exact eigenvalues. A term
    (L,S) contributes exactly one microstate to every (ML,MS) with
    |ML|<=L, |MS|<=S, ML/MS present with the given (L,S) fixed underlying
    energy E(L,S) (Racah's theorem: electrostatic energy depends only on
    L,S within one configuration). We peel from the highest (ML,MS) corner.
    Returns list of dicts: {L,S,E,deg}.
    """
    # collect available (ML,MS,E) multiset: expand eigenvalue multiplicities
    pool = defaultdict(list)  # (ML,MS) -> list of exact energies (with multiplicity)
    for (ML, MS), (idxs, eigs) in sector_eigs.items():
        for e, mult in eigs.items():
            pool[(ML, MS)].extend([e] * mult)

    terms = []
    # repeat until pool empty
    while any(len(v) > 0 for v in pool.values()):
        # find max ML present (with any MS), among nonempty pools; pick the corner
        # with largest ML, and among those the largest MS
        candidates = [(ML, MS) for (ML, MS), v in pool.items() if len(v) > 0]
        ML0 = max(c[0] for c in candidates)
        MS0 = max(c[1] for c in [c for c in candidates if c[0] == ML0])
        # the energy value that seeds this term = one of the entries at (ML0,MS0)
        E0 = pool[(ML0, MS0)][0]
        L = ML0
        S_ = MS0
        # remove one entry with value E0 from every (ML,MS) with |ML|<=L,|MS|<=S_
        # verifying that such an entry with exactly value E0 exists (this IS the
        # validation that the assumed term structure is consistent with the
        # exact diagonalisation)
        for MLp in range(-L, L + 1):
            for MSp_num in range(-int(2 * S_), int(2 * S_) + 1, 2):
                MSp = Rational(MSp_num, 2)
                lst = pool[(MLp, MSp)]
                match = None
                for e in lst:
                    if simplify(e - E0) == 0:
                        match = e
                        break
                assert match is not None, (
                    f"Term-peeling failed: no eigenvalue {E0} found at "
                    f"(ML,MS)=({MLp},{MSp}) for assumed term L={L},S={S_}. "
                    f"Available: {lst}"
                )
                lst.remove(match)
        deg = (2 * L + 1) * (2 * S_ + 1)
        terms.append({'L': L, 'S': S_, 'E': E0, 'deg': deg})
    return terms


TERM_LETTERS = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H', 6: 'I'}


def term_symbol(t):
    Smult = int(2 * t['S'] + 1)
    return f"{Smult}{TERM_LETTERS[t['L']]}"


def run_p2():
    print("=" * 70)
    print("p^2 shell (l=1), Slater-Condon F0, F2 (kmax=2)")
    print("=" * 70)
    sp_states, pairs, buckets, sector_eigs = diagonalize_shell(1, 2, 'p2')
    total_dim = len(pairs)
    print(f"single-particle states: {len(sp_states)}, two-electron basis dim: {total_dim} "
          f"(C(6,2)={Rational(6*5,2)})")
    assert total_dim == 15

    terms = peel_terms(sector_eigs)
    terms.sort(key=lambda t: t['E'].subs({F0: 0, F2: 1}))
    total_deg = sum(t['deg'] for t in terms)
    assert total_deg == 15, total_deg
    print(f"\nterms found: {len(terms)} (expect 3: 3P,1D,1S)")
    for t in terms:
        coeff_F2 = simplify(t['E'] - F0)
        print(f"  {term_symbol(t):>3s}  E = F0 + ({coeff_F2})   deg={t['deg']}")

    # exact validation against the textbook scaled-F2 result F0 + {-5,+1,+10} F2_scaled,
    # F2_scaled = F2_raw/25
    scaled = {t['L']: simplify((t['E'] - F0) * 25) for t in terms}
    lookup = {term_symbol(t): scaled[t['L']] for t in terms}
    assert simplify(lookup['3P'] - (-5) * F2) == 0, lookup
    assert simplify(lookup['1D'] - 1 * F2) == 0, lookup
    assert simplify(lookup['1S'] - 10 * F2) == 0, lookup
    print("VALIDATED: scaled coefficients {3P:-5, 1D:+1, 1S:+10} exact match.")
    assert term_symbol(min(terms, key=lambda t: t['E'].subs({F0: 0, F2: 1}))) == '3P'
    print("VALIDATED: ground term (min E at F2>0) = 3P (Hund's rule).")
    return terms


def run_d2():
    print("\n" + "=" * 70)
    print("d^2 shell (l=2), Slater-Condon F0, F2, F4 (kmax=4)")
    print("=" * 70)
    sp_states, pairs, buckets, sector_eigs = diagonalize_shell(2, 4, 'd2')
    total_dim = len(pairs)
    print(f"single-particle states: {len(sp_states)}, two-electron basis dim: {total_dim} "
          f"(C(10,2)={Rational(10*9,2)})")
    assert total_dim == 45

    terms = peel_terms(sector_eigs)
    terms.sort(key=lambda t: t['E'].subs({F0: 0, F2: 1, F4: 0}))
    total_deg = sum(t['deg'] for t in terms)
    assert total_deg == 45, total_deg
    print(f"\nterms found: {len(terms)} (expect 5: 3F,3P,1G,1D,1S)")
    for t in terms:
        print(f"  {term_symbol(t):>3s}  E = {simplify(t['E'])}   deg={t['deg']}")

    # convert to Racah A,B,C: F2_raw = 49*F2s, F4_raw = 441*F4s (standard scaling),
    # A = F0 - 49*F4s, B = F2s - 5*F4s, C = 35*F4s  =>  F0 = A + 7C/5, F2s = B + C/7, F4s = C/35
    A, B, C = symbols('A B C')
    F2s, F4s = symbols('F2s F4s')
    subs_scaled = {F2: 49 * F2s, F4: 441 * F4s}
    subs_racah = {F0: A + Rational(49, 35) * C, F2s: B + Rational(1, 7) * C, F4s: Rational(1, 35) * C}

    by_term = {term_symbol(t): t for t in terms}
    racah_energy = {}
    for name, t in by_term.items():
        e = simplify(t['E'].subs(subs_scaled))
        e = simplify(e.subs({F2s: F2s, F4s: F4s}))
        e = simplify(e.subs(subs_racah))
        e = simplify(e.expand())
        racah_energy[name] = e
        print(f"  {name:>3s}  (Racah)  E = {e}")

    expected_racah = {
        '3F': A - 8 * B,
        '3P': A + 7 * B,
        '1G': A + 4 * B + 2 * C,
        '1D': A - 3 * B + 2 * C,
        '1S': A + 14 * B + 7 * C,
    }
    for name, expr in expected_racah.items():
        diff = simplify(racah_energy[name] - expr)
        assert diff == 0, (name, diff, racah_energy[name], expr)
    print("VALIDATED: all 5 d^2 term energies match the textbook Racah A,B,C formulas exactly:")
    for name in ['3F', '3P', '1G', '1D', '1S']:
        print(f"    E({name}) = {expected_racah[name]}")
    assert term_symbol(min(terms, key=lambda t: t['E'].subs({F0: 0, F2: 1, F4: 0}))) == '3F'
    print("VALIDATED: ground term (min E at F2,F4>0) = 3F (Hund's rule).")
    return terms


def part_c_tests(p2_terms, d2_terms):
    print("\n" + "=" * 70)
    print("PART C -- pre-registered fingerprint tests (exact)")
    print("=" * 70)

    # Test 1: phi-in-energies -- rationality theorem, not a scan.
    print("\n[Test 1: phi-in-energies]")
    all_ratios_rational = True
    p2_by_L = {t['L']: t for t in p2_terms}
    e3P = simplify(p2_by_L[1]['E'] - F0)
    e1D = simplify(p2_by_L[2]['E'] - F0)
    e1S = simplify(p2_by_L[0]['E'] - F0)
    r1 = simplify((e1D - e3P) / (e1S - e3P))
    r2 = simplify((e1S - e1D) / (e1D - e3P))
    for name, r in [('r1=(E1D-E3P)/(E1S-E3P)', r1), ('r2=(E1S-E1D)/(E1D-E3P)', r2)]:
        assert r.is_rational, (name, r)
        print(f"  {name} = {r}  (rational: {r.is_rational})")
    print("  THEOREM: every term-energy ratio here is a ratio of INTEGER Slater-Condon")
    print("  coefficients (c^k are rational for l<=2 diagonal/degenerate blocks; all")
    print("  eigenvalues found above are RATIONAL multiples of F0,F2,(F4)). A rational")
    print("  number can equal phi, 1/phi, phi^2, sqrt5, phi/sqrt5, etc. ONLY if that")
    print("  irrational number is secretly rational -- which phi, sqrt5 provably are not")
    print("  (phi=(1+sqrt5)/2 has minimal polynomial x^2-x-1, degree 2 over Q; sqrt5 has")
    print("  minimal polynomial x^2-5, degree 2 over Q; neither is degree 1, i.e. neither")
    print("  is in Q). Hence NO term-energy ratio in p^2 or d^2 can equal any element of")
    print("  {phi,1/phi,phi^2,1/phi^2,2/phi,phi/2,sqrt5,1/sqrt5,phi/sqrt5}. PASS (by proof,")
    print("  not by numerical search) -- this is a structural NULL result.")

    # Test 2: fusion-multiplicity map, counting exercise
    print("\n[Test 2: fusion-multiplicity map]")
    print("  SU(2)_3 primaries: {0,1/2,1,3/2} -> 4 primaries, 3 non-vacuum.")
    print("  quantum dims: {1, phi, phi, 1}; total quantum dim^2 D^2 = 2phi*sqrt5.")
    print(f"  p-shell: {len(p2_terms)} terms ({[term_symbol(t) for t in p2_terms]}),")
    p2_S_gt0 = sum(1 for t in p2_terms if t['S'] > 0)
    print(f"    of which {p2_S_gt0} has S>0 (the triplet 3P).")
    print(f"  degeneracy multiset {{9,5,1}} vs SU(2)_3 small integers {{1,4,3}}(primary count),")
    print("    {1,phi,phi,1}(dims, irrational), D^2=2phi*sqrt5(irrational): NO exact numeric")
    print("    coincidence beyond both containing small integers already forced by SO(3)")
    print("    dimension counting ((2L+1)(2S+1) is pure Racah/SO(3), nothing to do with k=3).")
    print("  a-priori small-integer SU(2)_3 quantities to test against (informal count):")
    print("    4 (primary count), 3 (non-vacuum count), {1,4,9} (SO(3)-type dims j=0,1,3/2^2..),")
    print("    fusion coefficients (all 0/1, ~6 independent N_ab^c), S-matrix entries (irrational,")
    print("    can never match a rational p/d-shell ratio per Test 1). Roughly 10-15 small-")
    print("    integer/rational targets exist; matching 1 or 2 of them (e.g. '3 terms', '3 non-")
    print("    vacuum primaries') by chance among small integers 1-10 is unremarkable.")
    print("  RESULT: term count 3 vs non-vacuum-primary count 3 -- a shared SMALL INTEGER,")
    print("  not a fusion-rule correspondence (no fusion coefficient, dimension, or S-matrix")
    print("  entry is reproduced exactly). FAIL as a structural match; NOTED as a coincidence.")

    # Test 3: truncation test
    print("\n[Test 3: truncation test]")
    maxS_p2 = max(t['S'] for t in p2_terms)
    maxL_p2 = max(t['L'] for t in p2_terms)
    print(f"  p^2: max S = {maxS_p2}, max L = {maxL_p2}. Both are the ORDINARY Pauli-exclusion")
    print("  bound for 2 electrons in an l=1 shell (max S=1 since only 2 electrons; max L=2l-1=2")
    print("  minus antisymmetry reduces the L=2 *symmetric* triplet, but 1D itself (L=2) DOES")
    print("  occur) -- i.e. these are SO(3)+Pauli bounds, not a k=3 Hill-type mode cutoff.")
    print("  There is no 'number of cooperative modes' object here that saturates at exactly 3;")
    print("  the integer 3 (term count) is a DIFFERENT quantity from SU(2)_3's k+2=5 truncation")
    print("  index or its 3 non-vacuum primaries. FAIL: no genuine truncation-at-3 structure;")
    print("  the shell has its own truncation (Pauli + SO(3) triangle rules), unrelated in")
    print("  origin to the WZW level truncation.")

    # Test 4: the pentagon/5
    print("\n[Test 4: the pentagon/5]")
    print("  1D degeneracy = 5 = 2L+1 with L=2 -- this is SO(3) (an ordinary orbital")
    print("  degeneracy), NOT k+2=5. The two '5's have unrelated origins: one counts")
    print("  magnetic sublevels of L=2, the other counts the truncated primary/level")
    print("  structure of an anyonic fusion category. r1 = 2/5 (p-shell) is a ratio of")
    print("  Slater-Condon integers (6/15 reduced); it is RATIONAL, so it cannot equal")
    print("  any of the golden-ratio-family invariants either (Test 1). CONCLUSION: 5")
    print("  appears meaningfully as 2L+1 (SO(3)), not as k+2 (WZW truncation) -- these")
    print("  are the SO(3) '5' and the fusion-category '5' respectively; the numeral")
    print("  coincides, the origin does not.")


if __name__ == '__main__':
    p2_terms = run_p2()
    d2_terms = run_d2()
    part_c_tests(p2_terms, d2_terms)
    print("\nDONE.")
