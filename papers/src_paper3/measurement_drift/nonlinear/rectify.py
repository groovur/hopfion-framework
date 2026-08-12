"""
Nonlinear-rectification ("Tacoma Narrows") test for the k-essence dielectric
cell, Paper III context.

Context (see ../drift.py, ../results.md): the LINEAR measurement model
(diagonal Kraus maps on the fusion-channel density matrix) was proved an
exact martingale -- repeated weak measurement cannot bias the mean of the
measured coupling, symbolically, to all orders in the fringe parameter t.
That result leaves open the stated remaining hypothesis: a NONLINEAR vacuum
response could rectify a periodic measurement drive into a DC (zero-
frequency) shift of the effective dielectric response, which would bias the
measured coupling even though no linear mechanism can.

This script tests that hypothesis directly in the k-essence dielectric-cell
model used throughout src_paper3 (see ../../dielectric_cell/common.py, whose
validated functions are copied, not modified, below).

Model
-----
k-essence Lagrangian L(X) = X/(1+bX), b = b* = 0.06715131 (root of
(15/8)*arctan(sqrt(8b))/sqrt(8b) = phi, matching common.BETA_STAR exactly).
Background single-tube profile X0(rho) = 8/(1+rho^2)^2.

A measurement event is modeled as a weak oscillatory probe driving the
scalar fluctuation about the background,
    delta X(rho, tau) = A * f(rho) * cos(omega * tau),   A << 1.
The tangential/radial effective permittivities are nonlinear functions of X,
    eps_tan(X) = 1/(1+bX)^2,          eps_rad(X) = (1-3bX)/(1+bX)^3,
so a period average of eps(X0 + delta X) over one drive cycle produces an
O(A^2) DC rectification term (Duffing-style):
    <eps(X0 + delta X)>_tau = eps(X0) + (1/4) A^2 f(rho)^2 eps''(X0) + O(A^4),
using <cos^2> = 1/2 and the O(A) term averaging to zero.

R1: closed-form eps_tan''(X), eps_rad''(X) (sympy) and sign of the
    rectified shift over rho.
R2: profile-shaped drive f(rho) = X0(rho)/8; rerun the validated shooting
    solver with eps -> eps + (A^2/4) f^2 eps'', extract alpha_pol(A) and the
    fractional shift r(A) = [alpha_pol(A) - alpha_pol(0)] / alpha_pol(0) for
    A^2 in {1e-2, 1e-3, 1e-4} at eta = 1e-4. Fit r ~ C_rect * A^2.
R3: same with a drive localized on the eps_rad = 0 resonant shell,
    f(rho) = exp(-(rho - rho_c)^2 / (2*0.1^2)), compare C_rect.
R4: sign chain from alpha_pol -> eps_eff -> alpha_inv, and the amplitude A^2
    the observed four-term-formula residual would imply under this
    mechanism (inference only, not a fit).
"""
import json
import os
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# =======================================================================
# Functions copied verbatim (unmodified) from ../../dielectric_cell/common.py
# =======================================================================
PHI = (1 + np.sqrt(5)) / 2


def X(rho):
    return 8.0 / (1.0 + rho**2)**2


def find_x_star():
    def f(x):
        return (15.0 / 8.0) * np.arctan(x) / x - PHI
    x_star = brentq(f, 1e-6, 10.0, xtol=1e-14, rtol=1e-14)
    return x_star


X_STAR = find_x_star()
BETA_STAR = X_STAR**2 / 8.0


def rho_c_of_beta(beta):
    val = 24.0 * beta
    if val < 1.0:
        return None
    return np.sqrt(np.sqrt(val) - 1.0)


def eps_tan_bg(rho, beta):
    Xr = X(rho)
    return 1.0 / (1.0 + beta * Xr)**2


def eps_rad_bg(rho, beta):
    Xr = X(rho)
    return (1.0 - 3.0 * beta * Xr) / (1.0 + beta * Xr)**3


def rhs_factory(eps_rad_fn, eps_tan_fn, eta):
    def f(rho, y):
        w, v = y
        er = eps_rad_fn(rho) + 1j * eta
        et = eps_tan_fn(rho) + 1j * eta
        dw = v / (rho * er)
        dv = et * w / rho
        return np.array([dw, dv], dtype=complex)
    return f


def shoot_alpha(eps_rad_fn, eps_tan_fn, eta, rho_min=1e-4, rho_max=200.0,
                 rtol=1e-10, atol=1e-12, dense=False):
    f = rhs_factory(eps_rad_fn, eps_tan_fn, eta)
    er0 = eps_rad_fn(rho_min) + 1j * eta
    et0 = eps_tan_fn(rho_min) + 1j * eta
    s = np.sqrt(et0 / er0)
    if s.real < 0:
        s = -s
    w0 = rho_min**s
    wp0 = s * rho_min**(s - 1.0)
    v0 = rho_min * er0 * wp0
    y0 = np.array([w0, v0], dtype=complex)

    sol = solve_ivp(f, (rho_min, rho_max), y0, method='RK45',
                     rtol=rtol, atol=atol, dense_output=dense)
    if not sol.success:
        return None, {'success': False, 'message': sol.message}

    wR, vR = sol.y[:, -1]
    erR = eps_rad_fn(rho_max) + 1j * eta
    wpR = vR / (rho_max * erR)
    R = rho_max
    c1 = (wR + R * wpR) / (2.0 * R)
    c2 = R * (wR - R * wpR) / 2.0
    alpha = -c2 / c1
    info = {'success': True, 'nfev': sol.nfev, 'c1': c1, 'c2': c2,
            'w_end': wR, 'wp_end': wpR}
    return alpha, info


def bg_eps_fns(beta):
    return (lambda rho: eps_rad_bg(rho, beta), lambda rho: eps_tan_bg(rho, beta))


# sanity: reproduce the validated BETA_STAR/rho_c values from common.py
assert abs(BETA_STAR - 0.06715131) < 1e-8
RHO_C = rho_c_of_beta(BETA_STAR)
assert abs(RHO_C - 0.51913431) < 1e-7

results = {}

# =======================================================================
# R1: closed-form second derivatives eps_tan''(X), eps_rad''(X) (sympy)
# =======================================================================
Xs, bs = sp.symbols('X b', positive=True)

eps_tan_expr = 1 / (1 + bs * Xs)**2
eps_rad_expr = (1 - 3 * bs * Xs) / (1 + bs * Xs)**3

eps_tan_pp = sp.simplify(sp.diff(eps_tan_expr, Xs, 2))
eps_rad_pp = sp.simplify(sp.diff(eps_rad_expr, Xs, 2))

# hand-derived closed forms (via u = b*X substitution), cross-checked below
us = sp.symbols('u', positive=True)  # u = b*X
eps_tan_pp_hand = 6 * bs**2 / (1 + bs * Xs)**4
eps_rad_pp_hand = 6 * bs**2 * (5 - 3 * bs * Xs) / (1 + bs * Xs)**5

assert sp.simplify(eps_tan_pp - eps_tan_pp_hand) == 0, "eps_tan'' hand form mismatch"
assert sp.simplify(eps_rad_pp - eps_rad_pp_hand) == 0, "eps_rad'' hand form mismatch"

eps_tan_pp_str = str(sp.nsimplify(eps_tan_pp_hand))
eps_rad_pp_str = str(sp.nsimplify(eps_rad_pp_hand))

eps_tan_pp_num = sp.lambdify((Xs, bs), eps_tan_pp_hand, 'numpy')
eps_rad_pp_num = sp.lambdify((Xs, bs), eps_rad_pp_hand, 'numpy')

b_val = BETA_STAR

# sign analysis over the physical range of u = b*X along the background
# u(rho=0) = b*X(0) = b*8 is the maximum value of u anywhere on the profile
# (X0(rho) is maximal and monotonically decreasing from rho=0)
u_max = b_val * X(0.0)
u_at_rho_c = b_val * X(RHO_C)  # should equal 1/3 exactly (eps_rad(rho_c)=0 defn)

# eps_tan'' = 6 b^2/(1+u)^4 > 0 for all u > -1: always positive (no sign change)
# eps_rad'' = 6 b^2 (5-3u)/(1+u)^5: positive for u < 5/3, negative for u > 5/3,
# zero at u = 5/3 (distinct from the eps_rad=0 shell at u=1/3)
u_sign_change_rad = sp.Rational(5, 3)

rho_grid = np.linspace(1e-3, 5.0, 2000)
X0_grid = X(rho_grid)
u_grid = b_val * X0_grid
eps_tan_pp_grid = eps_tan_pp_num(X0_grid, b_val)
eps_rad_pp_grid = eps_rad_pp_num(X0_grid, b_val)

results['R1'] = dict(
    eps_tan_pp_closed_form="6*b**2/(1+b*X)**4",
    eps_rad_pp_closed_form="6*b**2*(5-3*b*X)/(1+b*X)**5",
    u_max_at_rho0=float(u_max),
    u_at_rho_c=float(u_at_rho_c),
    u_sign_change_eps_rad_pp=float(u_sign_change_rad),
    eps_tan_pp_sign="positive everywhere on the physical domain (u > -1 always)",
    eps_rad_pp_sign_at_rho0=("positive" if eps_rad_pp_grid[0] > 0 else "negative"),
    eps_rad_pp_min_over_grid=float(np.min(eps_rad_pp_grid)),
    eps_rad_pp_max_over_grid=float(np.max(eps_rad_pp_grid)),
    eps_rad_pp_at_rho_c=float(eps_rad_pp_num(X(RHO_C), b_val)),
    note="u never exceeds u_max<5/3 on this background, so eps_rad'' stays "
         "positive (same sign as eps_tan'') over the ENTIRE physical domain; "
         "the eps_rad=0 shell at u=1/3 is not a sign-change locus for eps_rad''.",
)

# =======================================================================
# R2: profile-shaped drive, response shift alpha_pol(A)
# =======================================================================
def f_profile(rho):
    return X(rho) / 8.0


def make_perturbed_eps(beta, f_fn, A2):
    def er(rho):
        Xr = X(rho)
        base = eps_rad_bg(rho, beta)
        return base + (A2 / 4.0) * f_fn(rho)**2 * eps_rad_pp_num(Xr, beta)

    def et(rho):
        Xr = X(rho)
        base = eps_tan_bg(rho, beta)
        return base + (A2 / 4.0) * f_fn(rho)**2 * eps_tan_pp_num(Xr, beta)
    return er, et


ETA_R2 = 1e-4
A2_LIST = [1e-2, 1e-3, 1e-4]

er0, et0 = bg_eps_fns(BETA_STAR)
alpha0, info0 = shoot_alpha(er0, et0, eta=ETA_R2, rtol=1e-10, atol=1e-12)
# cross-check against dielectric_cell/results.md value at beta_star, eta=1e-4:
# Re(alpha)=-0.97104087, Im(alpha)=1.36643754e-02
assert abs(alpha0.real - (-0.97104087)) < 1e-6
assert abs(alpha0.imag - 1.36643754e-02) < 1e-6


def run_A2_scan(f_fn, label):
    rows = []
    for A2 in A2_LIST:
        er, et = make_perturbed_eps(BETA_STAR, f_fn, A2)
        alpha, info = shoot_alpha(er, et, eta=ETA_R2, rtol=1e-10, atol=1e-12)
        r = (alpha - alpha0) / alpha0
        rows.append(dict(A2=A2, alpha_re=alpha.real, alpha_im=alpha.imag,
                          r_re=r.real, r_im=r.imag, C_rect=(r / A2)))
    return rows


rows_profile = run_A2_scan(f_profile, "profile")

# log-log slope fit of |r| vs A^2 (expect slope ~1: r ~ A^2)
A2_arr = np.array([row['A2'] for row in rows_profile])
r_abs_arr = np.array([abs(complex(row['r_re'], row['r_im'])) for row in rows_profile])
slope_profile, intercept_profile = np.polyfit(np.log(A2_arr), np.log(r_abs_arr), 1)

C_rect_profile_vals = [row['C_rect'] for row in rows_profile]
C_rect_profile_mean = np.mean(C_rect_profile_vals)

results['R2'] = dict(
    eta=ETA_R2,
    alpha0_re=alpha0.real,
    alpha0_im=alpha0.imag,
    rows=[dict(A2=row['A2'], alpha_re=row['alpha_re'], alpha_im=row['alpha_im'],
               r_re=row['r_re'], r_im=row['r_im'],
               C_rect_re=row['C_rect'].real, C_rect_im=row['C_rect'].imag)
          for row in rows_profile],
    loglog_slope_r_vs_A2=float(slope_profile),
    C_rect_profile_mean_re=float(C_rect_profile_mean.real),
    C_rect_profile_mean_im=float(C_rect_profile_mean.imag),
)

# =======================================================================
# R3: resonant-shell drive, compare rectification coefficient
# =======================================================================
SIGMA_SHELL = 0.1


def f_shell(rho):
    return np.exp(-(rho - RHO_C)**2 / (2.0 * SIGMA_SHELL**2))


rows_shell = run_A2_scan(f_shell, "shell")
r_abs_arr_shell = np.array([abs(complex(row['r_re'], row['r_im'])) for row in rows_shell])
slope_shell, intercept_shell = np.polyfit(np.log(A2_arr), np.log(r_abs_arr_shell), 1)

C_rect_shell_vals = [row['C_rect'] for row in rows_shell]
C_rect_shell_mean = np.mean(C_rect_shell_vals)

enhancement = abs(C_rect_shell_mean) / abs(C_rect_profile_mean)

results['R3'] = dict(
    rho_c=RHO_C,
    sigma_shell=SIGMA_SHELL,
    rows=[dict(A2=row['A2'], alpha_re=row['alpha_re'], alpha_im=row['alpha_im'],
               r_re=row['r_re'], r_im=row['r_im'],
               C_rect_re=row['C_rect'].real, C_rect_im=row['C_rect'].imag)
          for row in rows_shell],
    loglog_slope_r_vs_A2=float(slope_shell),
    C_rect_shell_mean_re=float(C_rect_shell_mean.real),
    C_rect_shell_mean_im=float(C_rect_shell_mean.imag),
    enhancement_factor=float(enhancement),
    tacoma_supported_qualitatively=bool(enhancement > 10.0),
)

# =======================================================================
# R4: sign chain and implied amplitude
# =======================================================================
n_val = 0.0459
alpha_inv_formula = 137.035998
alpha_inv_measured = 137.035999177
residual = alpha_inv_formula - alpha_inv_measured  # < 0, formula low

# eps_eff(n, alpha) = 1 + 2 n alpha / (1 - n alpha);  d(eps_eff)/d(alpha) = 2n/(1-n alpha)^2 > 0
alpha_re0 = alpha0.real  # < 0 (anti-screening), the "bare"/formula value
eps_eff0 = 1.0 + 2.0 * n_val * alpha_re0 / (1.0 - n_val * alpha_re0)
d_eps_eff_d_alpha = 2.0 * n_val / (1.0 - n_val * alpha_re0)**2

# fractional-shift equality (alpha_inv proportional to eps_eff, so relative
# shifts are equal regardless of the unknown proportionality constant):
#   Delta(alpha_inv)/alpha_inv_formula = Delta(eps_eff)/eps_eff0
# We need the sign/magnitude of Delta(alpha_inv) = alpha_inv_measured - alpha_inv_formula = -residual (>0)
delta_alpha_inv_needed = alpha_inv_measured - alpha_inv_formula  # > 0
frac_needed = delta_alpha_inv_needed / alpha_inv_formula

# Delta(eps_eff) = d_eps_eff_d_alpha * Delta(alpha_pol),  Delta(alpha_pol) = alpha_re0 * C_rect * A^2
# => frac_needed = [d_eps_eff_d_alpha * alpha_re0 * C_rect * A^2] / eps_eff0
# Solve for A^2 using the measured (profile-shaped) C_rect (real part; this is the
# physically realized drive shape, R3's shell value is a hypothetical upper bound):
C_rect_used = C_rect_profile_mean.real

# sign of C_rect required to move alpha_inv the correct direction (up, since
# delta_alpha_inv_needed > 0 and d_eps_eff_d_alpha > 0 and alpha_re0 < 0):
#   need alpha_re0 * C_rect > 0  =>  since alpha_re0 < 0, need C_rect < 0
sign_C_rect_required = "negative" if (alpha_re0 * 1.0) < 0 else "positive"
# (the required sign of C_rect itself, given alpha_re0 < 0 and we need alpha_re0*C_rect>0)
required_C_rect_sign = "negative"

coefficient = d_eps_eff_d_alpha * alpha_re0 * C_rect_used / eps_eff0
A2_implied = frac_needed / coefficient if coefficient != 0 else float('nan')

results['R4'] = dict(
    n=n_val,
    alpha_inv_formula=alpha_inv_formula,
    alpha_inv_measured=alpha_inv_measured,
    residual_formula_minus_measured=residual,
    delta_alpha_inv_needed=delta_alpha_inv_needed,
    frac_needed=frac_needed,
    alpha_re0=alpha_re0,
    eps_eff0=eps_eff0,
    d_eps_eff_d_alpha=d_eps_eff_d_alpha,
    C_rect_profile_re_used=C_rect_used,
    required_C_rect_sign_for_correct_direction=required_C_rect_sign,
    actual_C_rect_profile_sign=("negative" if C_rect_used < 0 else "positive"),
    signs_consistent=bool((C_rect_used < 0) == (required_C_rect_sign == "negative")),
    coefficient_frac_alpha_inv_per_A2=coefficient,
    A2_implied=A2_implied,
    A_implied=(float(np.sqrt(abs(A2_implied))) if not np.isnan(A2_implied) else float('nan')),
    A_physically_tiny=bool(not np.isnan(A2_implied) and abs(A2_implied) < 1.0),
)

# =======================================================================
# dump
# =======================================================================
with open(os.path.join(OUTDIR, 'results_raw.json'), 'w') as fh:
    json.dump(results, fh, indent=2, default=str)

print("DONE")
print("BETA_STAR =", BETA_STAR, "RHO_C =", RHO_C)
print("R1: eps_tan_pp sign: always positive; eps_rad_pp sign at rho=0:",
      results['R1']['eps_rad_pp_sign_at_rho0'], " u_max=", u_max, " u_at_rho_c=", u_at_rho_c)
print("R2: alpha0 =", alpha0, " C_rect_profile_mean =", C_rect_profile_mean,
      " loglog_slope=", slope_profile)
print("R3: C_rect_shell_mean =", C_rect_shell_mean, " enhancement =", enhancement,
      " loglog_slope=", slope_shell)
print("R4: required C_rect sign =", required_C_rect_sign,
      " actual sign =", results['R4']['actual_C_rect_profile_sign'],
      " A2_implied =", A2_implied)
