"""
chsh_transition_curve.py
========================
CHSH(f,beta) transition curve for the f-parameterised condensate model.

Reproduces and validates all numerical claims of:
  F. Manfredi, "FN Hopfion Condensate — Paper IX:
  Division Algebras, the Born Rule, and the Tsirelson Bound
  from Anisotropic Flux Suppression" (2026).

The f-parameterised probability:
  p_f(theta) = cos^(2f)(theta/2) / (cos^(2f)(theta/2) + sin^(2f)(theta/2))
             = 1 / (1 + tan^(2f)(theta/2))

Special cases:
  f=1: p_QM(theta) = cos^2(theta/2)   [standard complex Born rule, C-valued]
  f=2: p_cond(theta) = quaternionic (H-valued) condensate probability
  f->inf: step function (projective measurement, ceiling 8/3)

With density feedback beta > 0:
  f_eff(theta; f0, beta) = f0 / (1 + beta * cos(theta))
  p_{f0,beta}(theta) = 1 / (1 + tan^(2*f_eff)(theta/2))

  At theta=0  (aligned):    f_eff = f0/(1+beta) < f0  [feedback softens]
  At theta=pi (anti-aligned): f_eff = f0/(1-beta) > f0  [no feedback, sharpens]

Paper IX numerical targets (see --validate):
  Thm 4.1:   C(pi/8)=0.6415, C(3pi/8)=0.2445, raw=0.886, rescaled=2.363
  Thm 7.2 f=1:      raw=0.6533, rescaled=1.307 (at standard angles)
  Thm 7.2 f=4:      raw=0.980,  rescaled=2.613
  Thm 8.3 f=2,b*:  C(pi/8)=0.6353, C(3pi/8)=0.2432, raw=0.878, rescaled=2.341
  Thm 8.3 f=4,b*:  C(pi/8)=0.728,  C(3pi/8)=0.250,  raw=0.975, rescaled=2.60
  Conj 9.18: xcond(pi/4) = sqrt(2)/(3/2) = 2*sqrt(2)/3, CHSH=8*sqrt(2)/3=3.771
  Thm 9.14:          CHSH_raw = sqrt(2), CHSH_rescaled (x2) = 2*sqrt(2) = Tsirelson
"""

import argparse
import math
import numpy as np
from scipy.special import gamma as gamma_fn


# ── Constants ─────────────────────────────────────────────────────────────────
TSIRELSON    = 2.0 * math.sqrt(2.0)   # 2√2 ≈ 2.828
BETA_STAR    = 0.452                   # density-feedback coupling (Paper IX/VIII)
PHI          = (1.0 + math.sqrt(5.0)) / 2.0   # golden ratio


# ── Core functions ─────────────────────────────────────────────────────────────

def x_f(theta, f):
    """Signed correlation x_f(theta) = 2*p_f(theta) - 1.

    Parameters
    ----------
    theta : array_like   mismatch angle(s)
    f     : float or array_like   exponent (scalar or same shape as theta)

    Returns
    -------
    array_like in [-1, 1]
    """
    theta = np.asarray(theta, dtype=float)
    c = np.cos(theta / 2.0) ** f
    s = np.sin(theta / 2.0) ** f
    denom = c ** 2 + s ** 2
    return (c ** 2 - s ** 2) / np.where(denom == 0, 1.0, denom)


def f_eff(theta, f0, beta):
    """Direction-dependent effective exponent with density feedback.

    f_eff(theta) = f0 / (1 + beta * cos(theta))

    Requires beta < 1 for positivity at all theta.
    At theta=0:   f_eff = f0/(1+beta)  < f0   (feedback softens near-aligned)
    At theta=pi:  f_eff = f0/(1-beta)  > f0   (no feedback, sharpens anti-aligned)
    At theta=pi/2: f_eff = f0           = f0   (no correction at perpendicular)
    """
    return f0 / (1.0 + beta * np.cos(theta))


def x_fb(theta, f0, beta):
    """Signed correlation with density feedback.

    Uses f_eff(theta; f0, beta) as the effective exponent.
    """
    fe = f_eff(theta, f0, beta)
    return x_f(theta, fe)


def sin2f_avg_exact(f):
    """Exact 2D planar average <sin^{2f}(theta)>_2D = Gamma(f+1/2)/(sqrt(pi)*Gamma(f+1)).

    Special values:
      f=1: 1/2
      f=2: 3/8
      f=3: 5/16
      f->inf: -> 0  (rescaling diverges, CHSH_raw -> 1 from below)
    """
    return gamma_fn(f + 0.5) / (math.sqrt(math.pi) * gamma_fn(f + 1.0))


def xcond(theta):
    """The condensate correlation xcond(theta) = 2*cos(theta)/(1 + cos^2(theta)).

    This is x_f at f=2 in closed form: the quaternionic (H-valued) observable.
    """
    c = np.cos(theta)
    return 2.0 * c / (1.0 + c ** 2)


def stokes_corr(psi):
    """Stokes/Malus correlation E^Stokes(psi) = -(1/2)*cos(2*psi).

    Ensemble average of cos(2*(phi-alpha)) * cos(2*(phi-beta)) over phi_n.
    The factor 1/2 is the Parseval factor for the m=2 mode of L^2(S^1).
    After the QM normalisation convention (x2): CHSH = 2*sqrt(2) = Tsirelson.
    """
    return -0.5 * np.cos(2.0 * psi)


# ── CHSH computation ───────────────────────────────────────────────────────────

def chsh_mc(f0, beta=0.0, N=1_000_000, seed=42,
            angles=None):
    """Compute CHSH_raw and CHSH_rescaled by Monte Carlo.

    Standard Bell-test angles (Paper IX):
      a=0, a'=pi/4, b=pi/8, b'=3pi/8

    The f-dependent rescaling denominator <sin^{2f}>_2D is exact (Gamma functions).
    For beta > 0 the rescaling still uses f0 (the geometric average is beta-independent,
    as established in Theorem 8.3 of Paper IX).

    Parameters
    ----------
    f0    : float   base measurement-strength exponent
    beta  : float   density-feedback coupling (0 = feedback-free, 0.452 = Paper VIII)
    N     : int     number of Monte Carlo samples
    seed  : int     random seed for reproducibility
    angles: tuple   (a1,a2,b1,b2) override; default (0, pi/4, pi/8, 3pi/8)

    Returns
    -------
    chsh_raw       : float
    chsh_rescaled  : float   = chsh_raw / <sin^{2f0}>_2D
    avg_sin2f      : float   exact rescaling denominator
    correlations   : dict    E11, E12, E21, E22 and autocorrelations C(pi/8), C(3pi/8)
    """
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2.0 * math.pi, N)

    if angles is None:
        a1, a2 = 0.0, math.pi / 4.0
        b1, b2 = math.pi / 8.0, 3.0 * math.pi / 8.0
    else:
        a1, a2, b1, b2 = angles

    def mismatch(angle):
        d = np.abs(angle - phi) % (2.0 * math.pi)
        return np.where(d > math.pi, 2.0 * math.pi - d, d)

    def E(alpha, beta_angle):
        tA = mismatch(alpha)
        tB = mismatch(beta_angle)
        if beta == 0.0:
            xA = x_f(tA, f0)
            xB = x_f(tB, f0)
        else:
            xA = x_fb(tA, f0, beta)
            xB = x_fb(tB, f0, beta)
        return -np.mean(xA * xB)

    E11 = E(a1, b1)
    E12 = E(a1, b2)
    E21 = E(a2, b1)
    E22 = E(a2, b2)

    chsh_raw = abs(E11 + E12) + abs(E21 - E22)
    avg_fdep  = sin2f_avg_exact(f0)     # f-dependent: Gamma(f+1/2)/...
    avg_fixed = 3.0 / 8.0               # fixed pcond rescaling (paper table)
    chsh_res_fixed = chsh_raw / avg_fixed   # used in paper Table 1
    chsh_res_fdep  = chsh_raw / avg_fdep    # used in Theorem 7.2

    # Autocorrelations C(pi/8) and C(3pi/8) for direct comparison with paper
    # C(psi) ≡ (1/2pi) int x(phi)*x(phi+psi) dphi ≈ -E(0, psi)
    C_pi8   = -E(0.0, math.pi / 8.0)
    C_3pi8  = -E(0.0, 3.0 * math.pi / 8.0)

    correlations = dict(E11=E11, E12=E12, E21=E21, E22=E22,
                        C_pi8=C_pi8, C_3pi8=C_3pi8)
    # Return fixed-rescaling by default (paper convention)
    return chsh_raw, chsh_res_fixed, avg_fixed, correlations, \
           chsh_res_fdep, avg_fdep


def chsh_stokes_exact():
    """Exact CHSH for the photon Stokes observable (Theorem 9.14, Paper IX).

    E^Stokes(psi) = -(1/2)*cos(2*psi)  [Parseval identity, m=2 mode]

    At optimal angles (0, pi/4, pi/8, -pi/8):
      All |a-b| = pi/8 except |a'-b'| = 3pi/8
      CHSH_raw = sqrt(2)
      After Stokes normalisation (x2): CHSH = 2*sqrt(2) = Tsirelson

    Returns
    -------
    chsh_raw : float   = sqrt(2)
    chsh_stokes : float   = 2*sqrt(2)  [after QM normalisation x2]
    """
    a1, a2 = 0.0, math.pi / 4.0
    b1, b2 = math.pi / 8.0, -math.pi / 8.0   # optimal photon angles

    E11 = stokes_corr(a1 - b1)   # psi = pi/8
    E12 = stokes_corr(a1 - b2)   # psi = pi/8
    E21 = stokes_corr(a2 - b1)   # psi = pi/8
    E22 = stokes_corr(a2 - b2)   # psi = 3pi/8

    raw = abs(E11 + E12) + abs(E21 - E22)
    return raw, 2.0 * raw


def chsh_postquantum_exact():
    """CHSH under the angle-doubling conjecture (Conjecture conj:9.18).

    E_cond,phot(psi) = -xcond(2*psi)
    At optimal angles (0, pi/4, pi/8, -pi/8):
      All |psi| = pi/8 or 3pi/8
      CHSH = 4 * xcond(pi/4) = 4 * 2*(1/sqrt(2)) / (1 + 1/2) = 8*sqrt(2)/3

    Returns
    -------
    xcond_pi4  : float   xcond(pi/4) = 2*(1/sqrt(2)) / (1+1/2) = 4*sqrt(2)/3 / 2
    chsh       : float   8*sqrt(2)/3 ≈ 3.771
    """
    xc = float(xcond(math.pi / 4.0))
    return xc, 4.0 * xc


# ── Exact analytical checks ────────────────────────────────────────────────────

def exact_checks():
    """Analytical cross-checks not requiring Monte Carlo."""
    s2   = math.sqrt(2.0)
    phi  = PHI

    checks = {}

    # A1 = 4 - 2*sqrt(2)
    checks['A1'] = 4.0 - 2.0 * s2

    # Parseval sum = 1/sqrt(2)
    r    = 3.0 - 2.0 * s2    # |ratio| = |-(3-2*sqrt(2))| = 3-2*sqrt(2) ≈ 0.172
    A1   = 4.0 - 2.0 * s2
    parseval = (A1 ** 2 / 2.0) / (1.0 - r ** 2)
    checks['parseval_sum'] = parseval   # should be 1/sqrt(2)

    # C(pi/8) and C(3pi/8) from exact Fourier series
    def C_exact(psi):
        return sum((A1 * (-r) ** k) ** 2 / 2.0 * math.cos((2 * k + 1) * psi)
                   for k in range(200))
    checks['C_pi8_exact']  = C_exact(math.pi / 8.0)
    checks['C_3pi8_exact'] = C_exact(3.0 * math.pi / 8.0)
    checks['CHSH_raw_f2']  = checks['C_pi8_exact'] + checks['C_3pi8_exact']
    checks['CHSH_res_f2']  = checks['CHSH_raw_f2'] / (3.0 / 8.0)

    # f=1 exact
    checks['CHSH_raw_f1'] = 0.5 * math.sqrt(1.0 + 1.0 / s2)
    checks['CHSH_res_f1'] = checks['CHSH_raw_f1'] / 0.5

    # Stokes
    stokes_raw, stokes_norm = chsh_stokes_exact()
    checks['stokes_raw']  = stokes_raw
    checks['stokes_norm'] = stokes_norm

    # Post-quantum
    xc, pq_chsh = chsh_postquantum_exact()
    checks['xcond_pi4']  = xc
    checks['pq_chsh']    = pq_chsh
    checks['8sqrt2_3']   = 8.0 * s2 / 3.0

    return checks


# ── CLI modes ──────────────────────────────────────────────────────────────────

def run_table(args):
    """Print the CHSH(f) transition table (original script behaviour)."""
    Ts = TSIRELSON
    print("CHSH(f) transition — Paper IX Theorem 7.1")
    print("=" * 72)
    print(f"{'f':>6}  {'beta':>6}  {'raw':>9}  {'<sin^2f>':>10}  "
          f"{'rescaled':>10}  {'%Tsir':>7}")
    print("-" * 72)

    beta = args.beta
    for f in args.f_values:
        raw, res, avg, _, res_fdep, avg_fdep = chsh_mc(f, beta=beta, N=args.N, seed=args.seed)
        print(f"{f:>6.2f}  {beta:>6.3f}  {raw:>9.5f}  {avg:>10.6f}  "
              f"{res:>10.5f}  {100.0 * res / Ts:>6.2f}%")

    print()
    print(f"Tsirelson 2*sqrt(2) = {Ts:.6f}")
    if beta > 0:
        print(f"Density feedback beta = {beta} (Paper VIII convention)")


def run_validate(args):
    """Reproduce and check all Paper IX numerical claims."""
    tol = args.tol
    N   = args.N
    s2  = math.sqrt(2.0)

    # ── Analytical checks (no MC needed) ──────────────────────────────────
    exact = exact_checks()
    rows_exact = [
        ("A1 = 4-2√2",           exact['A1'],           1.172,  0.001),
        ("Parseval sum = 1/√2",   exact['parseval_sum'], 1/s2,   1e-8),
        ("C(pi/8) f=2 exact",     exact['C_pi8_exact'],  0.6415, 0.001),
        ("C(3pi/8) f=2 exact",    exact['C_3pi8_exact'], 0.2445, 0.001),
        ("CHSH_raw f=2 exact",    exact['CHSH_raw_f2'],  0.886,  0.001),
        ("CHSH_res f=2 exact",    exact['CHSH_res_f2'],  2.363,  0.002),
        ("CHSH_raw f=1 exact",    exact['CHSH_raw_f1'],  0.6533, 0.001),
        ("CHSH_res f=1 exact",    exact['CHSH_res_f1'],  1.307,  0.001),
        ("Stokes CHSH_raw",       exact['stokes_raw'],   s2,     1e-10),
        ("Stokes CHSH (x2)",      exact['stokes_norm'],  2*s2,   1e-10),
        ("xcond(pi/4)",           exact['xcond_pi4'],    4*s2/3/2, 1e-8),
        ("PQ CHSH = 8√2/3",       exact['pq_chsh'],      8*s2/3, 1e-8),
    ]

    # ── Monte Carlo checks ────────────────────────────────────────────────
    print("Validating Paper IX numerical claims...")
    print(f"Monte Carlo N={N:,}  seed={args.seed}")
    print()

    # Thm 4.1: f=2, beta=0
    raw_f2, res_f2, _, corr_f2, *_ = chsh_mc(2.0, beta=0.0, N=N, seed=args.seed)
    # Thm 7.2: f=1
    raw_f1, _, _, corr_f1, res_f1_fdep, _ = chsh_mc(1.0, beta=0.0, N=N, seed=args.seed)
    # Thm 7.2: f=4
    raw_f4, res_f4, _, corr_f4, *_ = chsh_mc(4.0, beta=0.0, N=N, seed=args.seed)
    # Thm 8.3: f=2, beta*
    raw_fb2, res_fb2, _, corr_fb2, *_ = chsh_mc(2.0, beta=BETA_STAR, N=N, seed=args.seed)
    # Thm 8.3: f=4, beta*
    raw_fb4, res_fb4, _, corr_fb4, *_ = chsh_mc(4.0, beta=BETA_STAR, N=N, seed=args.seed)

    rows_mc = [
        # (label, computed, paper_target, tolerance)
        # Thm 4.1
        ("Thm 4.1  C(pi/8)",    corr_f2['C_pi8'],   0.6415, 0.002),
        ("Thm 4.1  C(3pi/8)",   corr_f2['C_3pi8'],  0.2445, 0.002),
        ("Thm 4.1  CHSH_raw",   raw_f2,             0.886,  0.003),
        ("Thm 4.1  CHSH_res",   res_f2,             2.363,  0.005),
        # Thm 7.2 f=1
        ("Thm 7.2  f=1  CHSH_raw",       raw_f1,              0.6533, 0.002),
        # f=1 rescaled uses f-dependent (1/0.5=2) per Thm 7.2(iii)
        ("Thm 7.2  f=1  CHSH_res (f-dep)", res_f1_fdep,      1.307,  0.004),
        # Thm 7.2 f=4
        ("Thm 7.2  f=4  CHSH_raw",     raw_f4,             0.980,  0.003),
        ("Thm 7.2  f=4  CHSH_res",     res_f4,             2.613,  0.005),
        # Thm 8.3 f=2, beta*
        ("Thm 8.3  f=2  C(pi/8)",     corr_fb2['C_pi8'],  0.6353, 0.002),
        ("Thm 8.3  f=2  C(3pi/8)",    corr_fb2['C_3pi8'], 0.2432, 0.002),
        ("Thm 8.3  f=2  CHSH_raw",    raw_fb2,            0.878,  0.003),
        ("Thm 8.3  f=2  CHSH_res",    res_fb2,            2.341,  0.005),
        # Thm 8.3 f=4, beta*
        ("Thm 8.3  f=4  C(pi/8)",     corr_fb4['C_pi8'],  0.728,  0.003),
        ("Thm 8.3  f=4  C(3pi/8)",    corr_fb4['C_3pi8'], 0.250,  0.003),
        ("Thm 8.3  f=4  CHSH_raw",    raw_fb4,            0.975,  0.003),
        ("Thm 8.3  f=4  CHSH_res",    res_fb4,            2.60,   0.01),
    ]

    all_rows = rows_exact + rows_mc
    passed = 0
    failed = 0

    print(f"  {'Claim':<42}  {'Computed':>10}  {'Target':>10}  {'Δ':>8}  {'':>4}")
    print("  " + "-" * 80)
    for label, computed, target, tol_check in all_rows:
        delta = abs(computed - target)
        ok = delta <= tol_check
        mark = "✓" if ok else "✗ FAIL"
        if ok: passed += 1
        else:  failed += 1
        print(f"  {label:<42}  {computed:>10.5f}  {target:>10.5f}  "
              f"{delta:>8.5f}  {mark}")

    print()
    print(f"Results: {passed} passed, {failed} failed "
          f"(tolerance: analytical=exact, MC≈{tol})")

    if failed > 0:
        print("WARNING: Some checks failed. Increase N (--N 7000000) for better MC accuracy.")

    # Stokes exact check
    stokes_raw, stokes_norm = chsh_stokes_exact()
    print()
    print(f"Thm stokes (exact, no MC):")
    print(f"  CHSH_raw  = {stokes_raw:.8f}  (target = sqrt(2) = {math.sqrt(2):.8f})")
    print(f"  CHSH (x2) = {stokes_norm:.8f}  (target = 2*sqrt(2) = {2*math.sqrt(2):.8f})")

    # Post-quantum conjecture
    xc, pq = chsh_postquantum_exact()
    print()
    print(f"Conj 9.18 (exact):")
    print(f"  xcond(pi/4) = {xc:.8f}  (= 2*(1/√2)/(1+1/2) = {8*math.sqrt(2)/3/4:.8f})")
    print(f"  CHSH        = {pq:.8f}  (target 8√2/3 = {8*math.sqrt(2)/3:.8f})")


def run_scan(args):
    """Scan CHSH over a grid of (f, beta) values."""
    print(f"CHSH(f, beta) surface scan  [N={args.N:,}]")
    print()
    print(f"{'f':>6}  {'beta':>6}  {'raw':>9}  {'rescaled':>10}  {'%Tsir':>7}")
    print("-" * 50)
    for f in args.f_values:
        for b in args.beta_values:
            raw, res, *_ = chsh_mc(f, beta=b, N=args.N, seed=args.seed)
            print(f"{f:>6.2f}  {b:>6.3f}  {raw:>9.5f}  {res:>10.5f}  "
                  f"{100.0*res/TSIRELSON:>6.2f}%")
        if len(args.beta_values) > 1:
            print()


def run_photon(args):
    """Print the Stokes and post-quantum photon results."""
    stokes_raw, stokes_norm = chsh_stokes_exact()
    xc, pq = chsh_postquantum_exact()
    s2 = math.sqrt(2.0)

    print("Photon polarization results (Paper IX Section 9)")
    print("=" * 60)
    print()
    print("Theorem thm:stokes (Stokes/Malus observable, m=2 Fourier mode):")
    print(f"  E^Stokes(psi) = -(1/2)*cos(2*psi)")
    print(f"  Parseval factor = 1/2  (m=2 mode of L^2(S^1))")
    print(f"  CHSH_raw        = {stokes_raw:.6f}  (= sqrt(2) = {s2:.6f})")
    print(f"  After x2 norm.  = {stokes_norm:.6f}  (= 2*sqrt(2) = Tsirelson)")
    print()
    print("Two CHSH ceilings:")
    print(f"  Lepton (pcond, Parseval 3/8, x8/3): ceiling = 8/3 = {8/3:.6f}")
    print(f"  Photon (Stokes, Parseval 1/2,  x2): ceiling = {stokes_norm:.6f} = Tsirelson")
    print(f"  Structural gap 2*sqrt(2) - 8/3 = {stokes_norm - 8/3:.6f}")
    print()
    print("Conjecture conj:9.18 (pcond coupling for spin-1):")
    print(f"  xcond(pi/4) = 2*(1/sqrt(2))/(1+1/2) = {xc:.6f}")
    print(f"  CHSH = 4*xcond(pi/4) = {pq:.6f}  (= 8*sqrt(2)/3 = {8*s2/3:.6f})")
    print(f"  Exceeds Tsirelson by {pq - TSIRELSON:.4f} ({100*(pq/TSIRELSON - 1):.2f}%)")
    print(f"  No-signaling: satisfied (integral of xcond(2*(phi-beta)) over phi = 0)")


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_float_list(s):
    return [float(x) for x in s.split(',')]


def main():
    parser = argparse.ArgumentParser(
        description="CHSH(f,beta) transition curve — Paper IX validation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Reproduce Paper IX Table (f transition, beta=0):
  python chsh_transition_curve.py

  # Run with density feedback beta* = 0.452:
  python chsh_transition_curve.py --beta 0.452

  # Validate ALL numerical claims in Paper IX:
  python chsh_transition_curve.py --validate

  # Higher accuracy validation (slower):
  python chsh_transition_curve.py --validate --N 5000000

  # Scan (f, beta) surface:
  python chsh_transition_curve.py --scan --f 1.0,2.0,4.0 --beta-values 0.0,0.2,0.452

  # Photon/Stokes results and post-quantum conjecture:
  python chsh_transition_curve.py --photon

  # Single (f, beta) pair:
  python chsh_transition_curve.py --f 4.0 --beta 0.452
"""
    )

    parser.add_argument('--validate', action='store_true',
                        help='Reproduce and check all Paper IX numerical claims')
    parser.add_argument('--scan', action='store_true',
                        help='Scan CHSH over a grid of (f, beta) values')
    parser.add_argument('--photon', action='store_true',
                        help='Print Stokes/Malus and post-quantum results')

    parser.add_argument('--f', type=parse_float_list,
                        default='1.0,1.2,1.5,1.8,2.0,2.5,3.0,4.0,6.0,10.0',
                        metavar='F[,F,...]',
                        help='Comma-separated f values (default: 1.0,1.2,...,10.0)')
    parser.add_argument('--beta', type=float, default=0.0,
                        help='Density-feedback coupling beta (default: 0.0; '
                             'Paper VIII value: 0.452)')
    parser.add_argument('--beta-values', type=parse_float_list,
                        default='0.0,0.452', metavar='B[,B,...]',
                        help='Beta values for --scan (default: 0.0,0.452)')
    parser.add_argument('--N', type=int, default=1_000_000,
                        help='Monte Carlo sample count (default: 1000000)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('--tol', type=float, default=0.005,
                        help='Tolerance for --validate checks (default: 0.005)')
    parser.add_argument('--plot', action='store_true',
                        help='Save a CHSH(f) curve plot as chsh_curve.png')

    args = parser.parse_args()
    args.f_values    = args.f
    args.beta_values = args.beta_values

    if args.validate:
        run_validate(args)
    elif args.scan:
        run_scan(args)
    elif args.photon:
        run_photon(args)
    else:
        run_table(args)

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            f_dense = np.linspace(1.0, 8.0, 60)
            rescaled_vals = []
            for f in f_dense:
                raw, res, *_ = chsh_mc(f, beta=args.beta, N=300_000, seed=args.seed)
                rescaled_vals.append(res)

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(f_dense, rescaled_vals, color='#1f77b4', lw=2,
                    label=r'CHSH$_{\rm rescaled}(f)$' +
                          (f', $\\beta={args.beta}$' if args.beta > 0 else ''))
            ax.axhline(TSIRELSON, color='red', ls='--', lw=1.5,
                       label=r'Tsirelson $2\sqrt{2}$')
            ax.axhline(8/3, color='purple', ls='-.', lw=1.2,
                       label=r'Lepton ceiling $8/3$')
            ax.axhline(2.0, color='gray', ls=':', lw=1.0,
                       label='Classical bound 2')
            ax.axvline(2.0, color='orange', ls=':', lw=1.2,
                       label=r'$f=2$ (pcond, $\mathbb{H}$)')
            ax.axvline(1.0, color='green', ls=':', lw=1.2,
                       label=r'$f=1$ (Born rule, $\mathbb{C}$)')
            ax.set_xlabel(r'Measurement exponent $f$', fontsize=12)
            ax.set_ylabel(r'CHSH$_{\rm rescaled}(f)$', fontsize=12)
            ax.set_title(
                'CHSH$(f)$ Transition Curve — Paper IX Theorem 7.1'
                + (f'\n$\\beta={args.beta}$' if args.beta > 0 else ''),
                fontsize=11)
            ax.legend(fontsize=9)
            ax.set_xlim(1, 8)
            ax.set_ylim(1.0, 3.2)
            plt.tight_layout()
            plt.savefig('chsh_curve.png', dpi=150, bbox_inches='tight')
            print("\nPlot saved: chsh_curve.png")
            plt.close()
        except ImportError:
            print("\n(matplotlib not available — install with: pip install matplotlib)")


if __name__ == "__main__":
    main()
