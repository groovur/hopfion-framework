#!/usr/bin/env python3
"""
fn_hopfion_solver.py
====================
Self-consistent solver for the density-feedback FN Hopfion.
Parametrised over torus major radius R0.

PHYSICS
-------
Minimises  E_fb = K_fb * J4  where

    K_fb     = J2a + mu* * J2iso_fb
    J2a      = integral  sin^4(f) * kern  dV
    J2iso_fb = integral  kern / (1+beta*kern)  dV
    J4       = integral  (F13^2 + F12^2 + F23^2)  dV
    kern     = |grad f|^2 + sin^2(f) * A(r,z)
    A(r,z)   = 1/D0 + 1/r^2,   D0 = (r-R0)^2 + z^2

    mu* = 3 - phi,   phi = (1+sqrt(5))/2,   lam = phi^6

Self-consistency conditions at the minimum:
    (A)  V    = J2iso_fb / J2a  =  phi           [virial balance]
    (B)  sopt = sqrt(lam*J4/K_fb)  =  1          [Derrick condition]

WZW UNIVERSALITY PREDICTION
----------------------------
Linking Scale Conjecture: at the self-consistent minimum,

    J4/J2a   = 2^(4/3) / phi^5   ≈  0.22721
    K_fb/J4  = lam / 2^(1/3)     ≈  14.2424
    sopt^6   = 2

These values should be INDEPENDENT of R0 (the universality test).

ALGORITHM (saddle-snapshot)
---------------------------
E_fb has no Derrick minimum in R^3 — gradient descent always collapses
toward f=0. The Hopfion is a saddle-point of E_fb. The solver finds it
by tracking  score = |V - phi| + |sopt - 1|  throughout each descent
and reverting to the best intermediate snapshot. A warm-start chain
(each beta seeds the next) keeps the profile near the saddle.

USAGE
-----
  # Original R0=3 run (reproduces the published result):
  python3 fn_hopfion_solver.py \
    --R0 3.0 \
    --beta_list "0.3,0.35,0.4,0.42,0.44,0.44721,0.452,0.47,0.47213,0.5,0.55" \
    --outdir outputs_R03

  # Universality test at R0=4:
  python3 fn_hopfion_solver.py --R0 4.0 --outdir outputs_R04

  # Universality test at R0=5:
  python3 fn_hopfion_solver.py --R0 5.0 --outdir outputs_R05

  # Quick test on small grid:
  python3 fn_hopfion_solver.py --R0 3.0 --grid SMALL

  # Custom beta scan:
  python3 fn_hopfion_solver.py --R0 3.0 --beta_list "0.40,0.45,0.452,0.47"

  # Override warm-start:
  python3 fn_hopfion_solver.py --R0 4.0 --warmstart my_profile.npy

OUTPUT (all in --outdir)
------------------------
  log.txt            full run log
  results.csv        one row per beta
  report.txt         human-readable summary
  f_beta{v}.npy      converged profile at each beta
"""

import numpy as np, time, os, csv, argparse, sys
from scipy.ndimage import zoom

# ── CLI ───────────────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--R0',         type=float, default=3.0,
               help='Torus major radius (default 3.0)')
p.add_argument('--outdir',     default=None,
               help='Output directory (default: outputs_R0{R0})')
p.add_argument('--warmstart',  default=None,
               help='Path to .npy warm-start profile (auto-detected if omitted)')
p.add_argument('--beta_list',  default=None,
               help='Comma-separated beta values, e.g. "0.40,0.45,0.452"')
p.add_argument('--max_steps',  type=int, default=None,
               help='Max gradient-descent steps per beta')
p.add_argument('--grid',       choices=['SMALL', 'LARGE'], default='LARGE',
               help='SMALL=128x128 h=0.20 (testing)  LARGE=256x256 h=0.12 (production)')
args = p.parse_args()

R0     = args.R0
OUTDIR = args.outdir or f'outputs_R0{R0:.1f}'
os.makedirs(OUTDIR, exist_ok=True)

# ── Grid ──────────────────────────────────────────────────────────────────────
if args.grid == 'SMALL':
    Nr, Nz, h, MAX_DEF, PEV = 128, 128, 0.20, 300_000, 10_000
else:
    Nr, Nz, h, MAX_DEF, PEV = 256, 256, 0.12, 800_000, 25_000

MAX_STEPS = args.max_steps or MAX_DEF

# ── Constants ─────────────────────────────────────────────────────────────────
phi     = (1.0 + 5.0**0.5) / 2.0   # 1.6180339887…
lam     = phi**6                    # 17.9442719…
mu      = 3.0 - phi                 # 1.3819660…
ep      = 1e-14

# WZW universality predictions
WZW_g    = 2.0**(4.0/3.0) / phi**5   # J4/J2a  ≈ 0.22721
WZW_leff = lam / 2.0**(1.0/3.0)      # K_fb/J4 ≈ 14.2424
WZW_s6   = 2.0                        # sopt^6

# ── Geometry ──────────────────────────────────────────────────────────────────
r_         = h * (np.arange(Nr) + 0.5)
z_         = h * np.arange(Nz)
R, Z       = np.meshgrid(r_, z_, indexing='ij')
vol        = 2.0 * np.pi * 2.0 * R * h**2
vol[:, 0] *= 0.5
D0         = (R - R0)**2 + Z**2 + ep   # rebuilt per R0 below

# ── Beta list: auto-scale around known beta*(R0=3)=0.452 ─────────────────────
def default_beta_list(R0):
    """
    Estimate beta* for a given R0 and return a scan list.

    Calibrated from the shifted-IC analysis:
      R0=3 → beta* ≈ 0.452  (confirmed)
      R0=4 → beta* ≈ 0.51   (shifted-IC estimate)
      R0=5 → beta* ≈ 0.55   (shifted-IC estimate)
    Linear fit: beta*(R0) ≈ 0.452 + 0.049*(R0 - 3)
    """
    beta_est = 0.452 + 0.049 * (R0 - 3.0)
    lo  = max(0.05, beta_est * 0.60)
    hi  = min(4.0,  beta_est * 1.80)
    # 8 coarse points across the range + 5 dense around estimate
    coarse = np.linspace(lo, hi, 8)
    dense  = np.linspace(beta_est * 0.85, beta_est * 1.15, 5)
    return sorted(set(np.round(np.concatenate([coarse, dense]), 5)))

if args.beta_list:
    BETAS = [float(x) for x in args.beta_list.split(',')]
else:
    BETAS = default_beta_list(R0)

# ── Logging ───────────────────────────────────────────────────────────────────
LOG = open(os.path.join(OUTDIR, 'log.txt'), 'w', buffering=1)

def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    print(s, file=LOG); LOG.flush()

log("=" * 70)
log("FN HOPFION — DENSITY-FEEDBACK SOLVER")
log("=" * 70)
log(f"R0       = {R0}")
log(f"Grid     = {Nr}x{Nz},  h={h}")
log(f"phi      = {phi:.10f}")
log(f"mu*      = {mu:.10f}")
log(f"lam      = {lam:.8f}")
log(f"MAX_STEPS= {MAX_STEPS:,}  per beta")
log(f"Beta list: {BETAS}")
log()
log("WZW universality predictions:")
log(f"  J4/J2a   = {WZW_g:.8f}  (= 2^(4/3)/phi^5)")
log(f"  K_fb/J4  = {WZW_leff:.8f}  (= lam/2^(1/3))")
log(f"  sopt^6   = {WZW_s6:.8f}")
log()

# ── BCs ───────────────────────────────────────────────────────────────────────
def apply_bc(f):
    f = np.clip(f, 0.0, np.pi)
    f[0, :]  = 0.0
    f[-1, :] = 0.0
    f[:, -1] = 0.0
    return f

# ── Build IC when no warm-start is available ──────────────────────────────────
def make_ic_bs(C=None):
    """
    BS-ansatz initial condition: f = 2*arctan(C / d).
    C controls tube radius; default C = 0.6*R0 (matches R0=3 original IC).
    """
    if C is None:
        C = 0.6 * R0
    d  = np.sqrt(D0)
    return apply_bc(2.0 * np.arctan(C / d))

# ── Warm-start loading ────────────────────────────────────────────────────────
WARMSTART_CANDIDATES = [
    args.warmstart,
    '/mnt/user-data/uploads/f_fb_beta0_45200.npy',   # R0=3 reference
    f'{OUTDIR}/f_beta_best.npy',                       # previous run in same outdir
]

def shift_to_R0(f_src, R0_src, R0_dst, h_grid):
    """
    Translate a profile in the r-direction so its torus centre moves
    from R0_src to R0_dst.  This is a rigid integer-pixel shift —
    the only correct warm-start strategy when changing R0.

    Under the shift: D0 = (r - R0_dst)^2 + z^2 is unchanged in
    shape; the profile torus simply sits at the new major radius.
    """
    shift = round((R0_dst - R0_src) / h_grid)
    N     = f_src.shape[0]
    f_out = np.zeros_like(f_src)
    if shift >= 0:
        f_out[shift:, :] = f_src[:N - shift, :]
    else:
        f_out[:N + shift, :] = f_src[-shift:, :]
    # re-apply BCs
    f_out[0, :] = 0.0; f_out[-1, :] = 0.0; f_out[:, -1] = 0.0
    return np.clip(f_out, 0.0, np.pi)


f_warm   = None
R0_ws    = None   # R0 the warm-start was generated at

for path in WARMSTART_CANDIDATES:
    if path and os.path.exists(path):
        raw = np.load(path)
        if raw.shape == (Nr, Nz):
            f_warm = np.clip(raw, 0.0, np.pi)
            log(f"Warm-start: {path}  (exact shape {raw.shape})")
        else:
            sc     = Nr / raw.shape[0]
            f_warm = np.clip(zoom(raw, sc, order=3)[:Nr, :Nz], 0.0, np.pi)
            log(f"Warm-start: {path}  (resized {raw.shape} → {f_warm.shape})")
        # Infer the R0 the warm-start was built for from the profile peak.
        # Use peak_r directly (not rounded to a grid multiple) so that a
        # profile already centred at R0 doesn't get a spurious ±1-pt shift.
        peak_r = r_[int(np.argmax(f_warm[:, 0]))]
        R0_ws  = peak_r                          # use actual peak location
        log(f"  Profile peak at r={peak_r:.3f}  →  R0_warmstart = {R0_ws:.3f}")
        break

if f_warm is None:
    log(f"No warm-start found — using BS ansatz IC  (C={0.6*R0:.2f}, R0={R0})")
    f_warm = make_ic_bs()
elif abs(R0_ws - R0) > 0.5 * h:
    # Warm-start is from a different R0: shift the profile to the new centre.
    log(f"Shifting warm-start from R0={R0_ws:.2f} → R0={R0:.2f} "
        f"(shift = {round((R0-R0_ws)/h)} grid pts)")
    f_warm = shift_to_R0(f_warm, R0_ws, R0, h)
    log(f"Shift complete.  New profile peak at r={r_[int(np.argmax(f_warm[:,0]))]:.3f}")
else:
    log(f"Warm-start R0={R0_ws:.2f} matches target R0={R0:.2f} — no shift needed.")

log()

# ── Core compute ──────────────────────────────────────────────────────────────
def compute(f, beta):
    """
    Returns J2a, J2iso, J2iso_fb, J4, E_fb, Force  at the given (f, beta).
    D0 is module-level, computed for the current R0.
    """
    f  = np.clip(f, 0.0, np.pi)
    fr = np.gradient(f, h, axis=0)
    fz = np.gradient(f, h, axis=1); fz[:, 0] = 0.0

    sf = np.sin(f); cf = np.cos(f)
    s2 = sf**2; s3 = s2*sf; s4 = s2*s2

    A   = 1.0/D0 + 1.0/R**2
    fDG = fr*(R - R0) + fz*Z

    F13 = -s2/D0 * fDG
    F12 =  s2/R  * fr
    F23 =  s2/R  * fz

    kern = fr**2 + fz**2 + s2*A

    J2a      = float(np.sum(s4  * kern * vol))
    J2iso    = float(np.sum(       kern * vol))
    J4       = float(np.sum((F13**2 + F12**2 + F23**2) * vol))

    fb       = 1.0 / (1.0 + beta * kern)
    fb2      = fb * fb
    J2iso_fb = float(np.sum(kern * fb * vol))

    K_fb = J2a + mu * J2iso_fb
    E_fb = K_fb * J4

    # ── EL force: -dE_fb/df ───────────────────────────────────────────────────
    loc_J2a      = (4*s3*cf * kern + s4 * 2*sf*cf * A) * vol
    fJ2a_r       =  2*s4 * fr * vol
    fJ2a_z       =  2*s4 * fz * vol;  fJ2a_z[:, 0] = 0.0

    loc_J2iso_fb = 2*sf*cf * A * fb2 * vol
    fJ2iso_fb_r  = 2 * fr      * fb2 * vol
    fJ2iso_fb_z  = 2 * fz      * fb2 * vol;  fJ2iso_fb_z[:, 0] = 0.0

    loc_J4 = 2*(F13*(-2*sf*cf/D0 * fDG)
              + F12*( 2*sf*cf/R  * fr)
              + F23*( 2*sf*cf/R  * fz)) * vol
    fJ4_r  = 2*(F13*(-s2/D0*(R-R0)) + F12*(s2/R)) * vol
    fJ4_z  = 2*(F13*(-s2/D0*Z)      + F23*(s2/R)) * vol;  fJ4_z[:, 0] = 0.0

    def div(a, b):
        return np.gradient(a, h, axis=0) + np.gradient(b, h, axis=1)

    dK_fb = (loc_J2a - div(fJ2a_r, fJ2a_z)) + \
            mu * (loc_J2iso_fb - div(fJ2iso_fb_r, fJ2iso_fb_z))
    dJ4   = loc_J4 - div(fJ4_r, fJ4_z)

    Force = -(dK_fb * J4 + K_fb * dJ4)

    return J2a, J2iso, J2iso_fb, J4, E_fb, Force


# ── Gradient-descent with saddle-snapshot ─────────────────────────────────────
def run_beta(beta, f_init):
    """
    Minimise E_fb at the given beta. Because the Hopfion is a saddle-point
    of E_fb (not a minimum), we track  score = |V-phi| + |sopt-1|  throughout
    and revert to the best intermediate snapshot at the end.
    """
    f      = apply_bc(f_init.copy())
    dt     = 1e-7
    Ep     = 1e30
    consec = 0
    t0     = time.time()

    best_score = 1e30
    f_best     = f.copy()
    best       = {}

    log(f"\n{'─'*68}")
    log(f"  beta = {beta:.6f}")
    log(f"{'─'*68}")

    for step in range(1, MAX_STEPS + 1):

        J2a, J2iso, J2iso_fb, J4, E_fb, Force = compute(f, beta)

        if step % PEV == 0 or step == 1:
            V     = J2iso_fb / J2a   if J2a  > 1e-12 else 0.0
            sopt  = (lam*J4/(J2a+mu*J2iso_fb))**0.5 if (J2a+mu*J2iso_fb) > 0 else 0.0
            sopt6 = sopt**6
            score = abs(V - phi) + abs(sopt - 1.0)

            log(f"  step={step:>8,}  E={E_fb:>12.3f}"
                f"  V={V:.5f}(φ={phi:.5f})"
                f"  sopt={sopt:.5f}"
                f"  sopt6={sopt6:.4f}"
                f"  score={score:.5f}"
                f"  dt={dt:.1e}"
                f"  [{time.time()-t0:.0f}s]")

            if score < best_score:
                best_score = score
                f_best     = f.copy()
                best = dict(J2a=J2a, J2iso=J2iso, J2iso_fb=J2iso_fb,
                            J4=J4, E_fb=E_fb, V=V, sopt=sopt, sopt6=sopt6,
                            score=score, step=step)

        f_try = apply_bc(f + dt * Force)
        _, _, _, _, E_try, _ = compute(f_try, beta)

        if E_try < E_fb:
            f  = f_try
            dt = min(dt * 1.005, 1e-3)
        else:
            dt *= 0.7
            if dt < 1e-18:
                log(f"  dt exhausted at step {step}.")
                break

        rel = abs(E_fb - Ep) / max(abs(E_fb), 1e-30)
        if rel < 1e-10 and step > 10_000:
            consec += 1
            if consec >= 5:
                log(f"  Converged at step {step:,}  [{time.time()-t0:.0f}s]")
                break
        else:
            consec = 0
        Ep = E_fb

    # ── Final eval & revert if needed ─────────────────────────────────────────
    J2a, J2iso, J2iso_fb, J4, E_fb, _ = compute(f, beta)
    V_f     = J2iso_fb/J2a   if J2a > 1e-12 else 0.0
    sopt_f  = (lam*J4/(J2a+mu*J2iso_fb))**0.5 if (J2a+mu*J2iso_fb) > 0 else 0.0
    score_f = abs(V_f - phi) + abs(sopt_f - 1.0)

    if best and score_f > best_score:
        log(f"  Reverting to snapshot at step {best['step']:,} "
            f"(score {best_score:.5f} < {score_f:.5f})")
        f         = f_best
        J2a       = best['J2a'];       J2iso    = best['J2iso']
        J2iso_fb  = best['J2iso_fb'];  J4       = best['J4']
        E_fb      = best['E_fb'];      V_f      = best['V']
        sopt_f    = best['sopt'];      score_f  = best['score']

    sopt6_f = sopt_f**6
    KJ4_f   = (J2a + mu*J2iso_fb) / J4 if J4 > 1e-12 else 0.0
    g_f     = J4 / J2a if J2a > 1e-12 else 0.0

    return f, dict(
        J2a=J2a, J2iso=J2iso, J2iso_fb=J2iso_fb, J4=J4, E_fb=E_fb,
        V=V_f, sopt=sopt_f, sopt6=sopt6_f, KJ4=KJ4_f, g=g_f,
        score=score_f,
        V_gap=V_f - phi,
        sopt6_gap_pct=100.0*(sopt6_f/WZW_s6 - 1.0),
        g_ratio=g_f/WZW_g if WZW_g > 0 else 0.0,
        KJ4_ratio=KJ4_f/WZW_leff if WZW_leff > 0 else 0.0,
    )


# ── Main loop ─────────────────────────────────────────────────────────────────
log("=" * 70)
log("STARTING BETA SCAN")
log("=" * 70)

FIELDS = ['beta', 'V', 'V_gap', 'sopt', 'sopt6', 'sopt6_gap_pct',
          'J4_J2a', 'g_ratio', 'Kfb_J4', 'KJ4_ratio',
          'J2a', 'J2iso', 'J2iso_fb', 'J4', 'E_fb', 'score']

CSV_PATH = os.path.join(OUTDIR, 'results.csv')
with open(CSV_PATH, 'w', newline='') as cf:
    csv.DictWriter(cf, FIELDS).writeheader()

results  = []
f_chain  = f_warm.copy()
best_run = None

for beta in BETAS:

    f_conv, res = run_beta(beta, f_chain)

    row = dict(beta=beta,
               V=res['V'], V_gap=res['V_gap'],
               sopt=res['sopt'], sopt6=res['sopt6'],
               sopt6_gap_pct=res['sopt6_gap_pct'],
               J4_J2a=res['g'], g_ratio=res['g_ratio'],
               Kfb_J4=res['KJ4'], KJ4_ratio=res['KJ4_ratio'],
               J2a=res['J2a'], J2iso=res['J2iso'],
               J2iso_fb=res['J2iso_fb'], J4=res['J4'],
               E_fb=res['E_fb'], score=res['score'])
    results.append(row)

    with open(CSV_PATH, 'a', newline='') as cf:
        csv.DictWriter(cf, FIELDS).writerow(row)

    np.save(os.path.join(OUTDIR, f'f_beta{beta:.5f}.npy'), f_conv)

    log(f"\n  ┌── RESULT  beta={beta:.6f} {'─'*38}┐")
    log(f"  │  V    = J2iso_fb/J2a = {res['V']:.8f}   (phi={phi:.8f}  gap={res['V_gap']:+.6f})")
    log(f"  │  sopt = sqrt(lam*J4/K_fb) = {res['sopt']:.8f}   (target 1.0)")
    log(f"  │  sopt^6              = {res['sopt6']:.8f}   (WZW target {WZW_s6:.1f}  gap={res['sopt6_gap_pct']:+.3f}%)")
    log(f"  │  J4/J2a              = {res['g']:.8f}   (WZW {WZW_g:.8f}  ratio={res['g_ratio']:.5f})")
    log(f"  │  K_fb/J4             = {res['KJ4']:.8f}   (WZW {WZW_leff:.8f}  ratio={res['KJ4_ratio']:.5f})")
    log(f"  │  score               = {res['score']:.8f}   (0 = perfect)")
    log(f"  └{'─'*66}┘")

    # Note crossing of V=phi
    if len(results) >= 2:
        d0 = results[-2]['V_gap']
        d1 = results[-1]['V_gap']
        if d0 * d1 < 0:
            b0, b1 = results[-2]['beta'], results[-1]['beta']
            b_x = b0 - d0*(b1-b0)/(d1-d0)
            log(f"\n  *** V crosses phi between beta={b0:.5f} and beta={b1:.5f}")
            log(f"  *** Interpolated beta_cross ≈ {b_x:.6f}")

    # Track global best
    if best_run is None or res['score'] < best_run['score']:
        best_run = dict(beta=beta, **res)
        np.save(os.path.join(OUTDIR, 'f_beta_best.npy'), f_conv)

    # Warm-start chain
    f_chain = f_conv.copy()


# ── Summary ───────────────────────────────────────────────────────────────────
log("\n\n" + "=" * 70)
log("FINAL SUMMARY")
log("=" * 70)
log(f"{'beta':>8}  {'V':>10}  {'V_gap':>9}  {'sopt':>8}  "
    f"{'sopt6':>8}  {'J4/J2a':>10}  {'K/J4':>10}  {'score':>10}")
log("─" * 88)

for r in results:
    flag = " ← BEST" if r['beta'] == best_run['beta'] else ""
    log(f"  {r['beta']:>6.5f}  {r['V']:>10.6f}  {r['V_gap']:>+9.5f}  {r['sopt']:>8.5f}  "
        f"{r['sopt6']:>8.5f}  {r['J4_J2a']:>10.6f}  {r['Kfb_J4']:>10.5f}  "
        f"{r['score']:>10.6f}{flag}")

log("=" * 70)
log()
log(f"Best self-consistent solution:  beta = {best_run['beta']:.6f}")
log(f"  V    = J2iso_fb/J2a = {best_run['V']:.8f}  (phi={phi:.8f}  gap={best_run['V_gap']:+.6f})")
log(f"  sopt = {best_run['sopt']:.8f}  (target 1.0)")
log(f"  sopt6= {best_run['sopt6']:.8f}  (WZW target {WZW_s6})")
log(f"  J4/J2a  = {best_run['g']:.8f}  (WZW {WZW_g:.8f}  ratio={best_run['g_ratio']:.5f})")
log(f"  K_fb/J4 = {best_run['KJ4']:.8f}  (WZW {WZW_leff:.8f}  ratio={best_run['KJ4_ratio']:.5f})")
log(f"  score   = {best_run['score']:.8f}")
log()
log("WZW UNIVERSALITY CHECK:")
log(f"  J4/J2a  ratio (obs/pred) = {best_run['g_ratio']:.5f}  {'✓ UNIVERSAL' if abs(best_run['g_ratio']-1)<0.05 else '✗ R0-DEPENDENT'}")
log(f"  K_fb/J4 ratio (obs/pred) = {best_run['KJ4_ratio']:.5f}  {'✓ UNIVERSAL' if abs(best_run['KJ4_ratio']-1)<0.05 else '✗ R0-DEPENDENT'}")
log(f"  sopt^6  = {best_run['sopt6']:.5f}  (target 2.0)  {'✓' if abs(best_run['sopt6']-2)<0.1 else '✗'}")

# ── Written report ────────────────────────────────────────────────────────────
rpt = os.path.join(OUTDIR, 'report.txt')
with open(rpt, 'w') as f:
    f.write(f"FN HOPFION DENSITY-FEEDBACK SOLVER\n")
    f.write(f"R0={R0}  Grid={Nr}x{Nz}  h={h}\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"WZW predictions:\n")
    f.write(f"  J4/J2a  = {WZW_g:.8f}\n")
    f.write(f"  K_fb/J4 = {WZW_leff:.8f}\n")
    f.write(f"  sopt^6  = {WZW_s6}\n\n")
    f.write(f"{'beta':>8}  {'V':>12}  {'V_gap':>9}  {'sopt':>8}  "
            f"{'sopt6':>8}  {'J4/J2a':>10}  {'K/J4':>10}  {'score':>10}\n")
    f.write("─" * 88 + "\n")
    for r in results:
        flag = " BEST" if r['beta'] == best_run['beta'] else ""
        f.write(f"  {r['beta']:>6.5f}  {r['V']:>12.8f}  {r['V_gap']:>+9.5f}  "
                f"{r['sopt']:>8.5f}  {r['sopt6']:>8.5f}  "
                f"{r['J4_J2a']:>10.6f}  {r['Kfb_J4']:>10.5f}  "
                f"{r['score']:>10.6f}{flag}\n")
    f.write("\n")
    f.write(f"Best: beta={best_run['beta']:.6f}\n")
    f.write(f"  V=      {best_run['V']:.8f}  (phi={phi:.8f})\n")
    f.write(f"  sopt=   {best_run['sopt']:.8f}\n")
    f.write(f"  sopt6=  {best_run['sopt6']:.8f}  (WZW 2.0)\n")
    f.write(f"  J4/J2a= {best_run['g']:.8f}  (WZW {WZW_g:.8f}  ratio={best_run['g_ratio']:.4f})\n")
    f.write(f"  K/J4=   {best_run['KJ4']:.8f}  (WZW {WZW_leff:.8f}  ratio={best_run['KJ4_ratio']:.4f})\n")
    f.write(f"  score=  {best_run['score']:.8f}\n")

log(f"\nAll outputs in: {OUTDIR}/")
log("  log.txt          — this log")
log("  results.csv      — data table")
log("  report.txt       — summary report")
log("  f_beta*.npy      — profiles per beta")
log("  f_beta_best.npy  — best profile (use as warm-start for next run)")
LOG.close()
