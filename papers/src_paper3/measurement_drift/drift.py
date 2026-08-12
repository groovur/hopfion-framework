"""
Measurement-drift simulation for the SU(2)_3 probe-pair monodromy coupling
(Paper III context). "Tacoma Narrows" hypothesis test: does repeated weak
measurement of the fusion-channel monodromy drive the Born-diagonal state
away from its initial value (secular bias / resonant pumping), or does the
martingale structure of the measurement forbid this in the mean, leaving
only a variance effect (and, once damping/relaxation competes against it,
a bounded steady-state fluctuation)?

Physics setup
--------------
Fusion channels c in {0,1}, quantum dims d_0=1, d_1=phi, Born weights
p0 = 1/phi^2, p1 = 1/phi, monodromy eigenvalues
    M_0 = exp(-i*108deg),  M_1 = exp(+i*36deg).
alpha_inv = 360 * Re Tr[rho U_mono].

Because rho stays diagonal in the fusion basis under both the measurement
Kraus maps (diagonal by construction) and the Born-state relaxation
(convex combination of diagonal states), the entire dynamics reduces
exactly to a scalar classical stochastic process on p0(n) = Pr(channel 0)
for each trajectory. Tr[rho_n U] = p0(n)*M0 + (1-p0(n))*M1, so
Re Tr[rho_n U] = b + p0(n)*(a-b), with a = Re M0, b = Re M1 -- an affine,
strictly monotonic function of p0(n). All statistics of interest are
therefore statistics of the scalar p0(n) pushed through this affine map.

One interferometer pass (fringe parameter t in (0,1]):
    K_+ = diag( sqrt((1+t*a)/2), sqrt((1+t*b)/2) )
    K_- = diag( sqrt((1-t*a)/2), sqrt((1-t*b)/2) )
    K_+^2 + K_-^2 = I   (checked below).
Outcome probabilities: p(s) = p0*(1+s*t*a)/2 + (1-p0)*(1+s*t*b)/2.
Update: p0 -> p0 * (1+s*t*a)/2 / p(s).

Relaxation (Tacoma competition): after each measurement pass,
    p0 -> (1-g)*p0 + g*p0_born.
"""
import numpy as np
import sympy as sp
import mpmath as mp
import json
import os

OUTDIR = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260812)

# ---------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------
phi = (1 + 5**0.5) / 2
P0_BORN = 1 / phi**2
P1_BORN = 1 / phi
A = float(mp.cos(mp.radians(108)))   # Re M0
B = float(mp.cos(mp.radians(36)))    # Re M1

# exact high-precision check that Tr[rho_born U] is real and equals 1/phi^2
mp.mp.dps = 50
phi_mp = (1 + mp.sqrt(5)) / 2
p0_mp = 1 / phi_mp**2
p1_mp = 1 / phi_mp
M0_mp = mp.e**(-1j * mp.radians(108))
M1_mp = mp.e**(1j * mp.radians(36))
trU_mp = p0_mp * M0_mp + p1_mp * M1_mp
assert abs(mp.im(trU_mp)) < mp.mpf('1e-45'), "Born Tr[rho U] not real!"
assert abs(mp.re(trU_mp) - p0_mp) < mp.mpf('1e-45'), "Born Tr[rho U] != 1/phi^2!"
ALPHA_INV_BORN = 360.0 * float(mp.re(trU_mp))

# ---------------------------------------------------------------------
# core vectorized single-pass measurement update (over N trajectories)
# ---------------------------------------------------------------------
def measure_step(p0, t, rng):
    X = p0 * A + (1 - p0) * B
    p_plus = 0.5 + 0.5 * t * X
    u = rng.random(p0.shape)
    plus_mask = u < p_plus
    p0_plus = p0 * (1 + t * A) / 2 / p_plus
    p0_minus = p0 * (1 - t * A) / 2 / (1 - p_plus)
    return np.where(plus_mask, p0_plus, p0_minus)

def relax_step(p0, g):
    return (1 - g) * p0 + g * P0_BORN

def re_x_of_p0(p0):
    return B + p0 * (A - B)

results = {}

# ---------------------------------------------------------------------
# sanity: Kraus completeness at a few (p0, t) points
# ---------------------------------------------------------------------
for t_test in [1.0, 0.5, 0.1, 0.01]:
    Kp2_0 = (1 + t_test * A) / 2
    Km2_0 = (1 - t_test * A) / 2
    Kp2_1 = (1 + t_test * B) / 2
    Km2_1 = (1 - t_test * B) / 2
    assert abs((Kp2_0 + Km2_0) - 1) < 1e-14
    assert abs((Kp2_1 + Km2_1) - 1) < 1e-14

# =======================================================================
# R1 + R2 : collapse statistics and measured-value drift, no relaxation
# =======================================================================
N_R1 = 10_000
N_STEPS_R1 = 1000
T_LIST = [1.0, 0.5, 0.1, 0.01]
CHECKPOINTS = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]

r1r2 = {}
for t in T_LIST:
    p0 = np.full(N_R1, P0_BORN)
    mean_p0_series = []
    mean_360rex_series = []
    var_360rex_series = []
    step_recorded = []
    for n in range(N_STEPS_R1 + 1):
        if n in CHECKPOINTS:
            rex = re_x_of_p0(p0)
            mean_p0_series.append(float(np.mean(p0)))
            mean_360rex_series.append(float(np.mean(360 * rex)))
            var_360rex_series.append(float(np.var(360 * rex)))
            step_recorded.append(n)
        if n < N_STEPS_R1:
            p0 = measure_step(p0, t, rng)

    # collapse-endpoint frequencies at n = N_STEPS_R1
    frac_channel0 = float(np.mean(p0 > 0.5))
    frac_channel1 = float(np.mean(p0 < 0.5))
    frac_undecided = 1.0 - frac_channel0 - frac_channel1
    # how "collapsed" trajectories are (should sit at 0 or 1)
    mean_dist_to_boundary = float(np.mean(np.minimum(p0, 1 - p0)))

    r1r2[t] = dict(
        steps=step_recorded,
        mean_p0=mean_p0_series,
        mean_360rex=mean_360rex_series,
        var_360rex=var_360rex_series,
        frac_channel0=frac_channel0,
        frac_channel1=frac_channel1,
        frac_undecided=frac_undecided,
        mean_dist_to_boundary=mean_dist_to_boundary,
    )

results['R1R2'] = r1r2
results['endpoints'] = dict(
    val_360ReM0=360.0 * A,
    val_360ReM1=360.0 * B,
    alpha_inv_born=ALPHA_INV_BORN,
    p0_born=P0_BORN,
    p1_born=P1_BORN,
)

# =======================================================================
# R3 : rethermalization competition (measurement + relaxation)
# =======================================================================
N_R3 = 2000
G_LIST = [1e-1, 1e-2, 1e-3, 1e-4]
r3 = {}
for t in T_LIST:
    for g in G_LIST:
        steps = min(int(20 / g), 60000)
        burn_frac = 0.7  # use last 30% of run as "steady state" window
        p0 = np.full(N_R3, P0_BORN)
        snap_means = []
        snap_vars = []
        record_from = int(burn_frac * steps)
        for n in range(steps):
            p0 = measure_step(p0, t, rng)
            p0 = relax_step(p0, g)
            if n >= record_from:
                rex = 360 * re_x_of_p0(p0)
                snap_means.append(float(np.mean(rex)))
                snap_vars.append(float(np.var(rex)))
        ss_mean = float(np.mean(snap_means))
        ss_std = float(np.mean(np.sqrt(snap_vars)))  # avg of per-snapshot std
        ss_bias = ss_mean - ALPHA_INV_BORN
        r3[(t, g)] = dict(steps=steps, ss_mean=ss_mean, ss_bias=ss_bias, ss_std=ss_std)

results['R3'] = {f"{t}|{g}": v for (t, g), v in r3.items()}

# log-log slope fits: sigma vs t (at fixed g, largest g for best equilibration),
# and sigma vs g (at fixed t).
def loglog_slope(xs, ys):
    lx = np.log(np.array(xs))
    ly = np.log(np.array(ys))
    slope, intercept = np.polyfit(lx, ly, 1)
    return float(slope), float(intercept)

# sigma vs t at fixed g = 1e-1 (best equilibrated)
g_fixed = 1e-1
sig_vs_t = [r3[(t, g_fixed)]['ss_std'] for t in T_LIST]
slope_sigma_t, _ = loglog_slope(T_LIST, sig_vs_t)

# sigma vs g at fixed t = 1.0
t_fixed = 1.0
sig_vs_g = [r3[(t_fixed, g)]['ss_std'] for g in G_LIST]
slope_sigma_g, _ = loglog_slope(G_LIST, sig_vs_g)

# bias magnitude vs t, g (should be ~noise floor, not a real power law)
bias_vals = [abs(r3[(t, g)]['ss_bias']) for t in T_LIST for g in G_LIST]

results['R3_fits'] = dict(
    slope_sigma_vs_t_at_g0p1_fullgrid=slope_sigma_t,
    slope_sigma_vs_g_at_t1_fullgrid=slope_sigma_g,
    max_abs_bias=float(np.max(bias_vals)),
    median_abs_bias=float(np.median(bias_vals)),
)

# clean small-t (perturbative, unsaturated) slope fits, as used in results.md:
# sigma vs g at fixed small t=0.01 (full 4-point g grid, unsaturated regime)
sig_vs_g_smallt = [r3[(0.01, g)]['ss_std'] for g in G_LIST]
slope_g_smallt, _ = loglog_slope(G_LIST, sig_vs_g_smallt)
# sigma vs t at fixed g=0.1, restricted to unsaturated t in {0.01,0.1}
sig_vs_t_small = [r3[(t, 0.1)]['ss_std'] for t in [0.01, 0.1]]
slope_t_smallt, _ = loglog_slope([0.01, 0.1], sig_vs_t_small)
results['R3_fits']['slope_sigma_vs_g_at_t0p01_unsaturated'] = slope_g_smallt
results['R3_fits']['slope_sigma_vs_t_unsaturated_at_g0p1'] = slope_t_smallt

# =======================================================================
# R4 : exact / small-t symbolic analysis (sympy)
# =======================================================================
p0s, ts, gs, as_, bs = sp.symbols('p0 t g a b', real=True)

X = p0s * as_ + (1 - p0s) * bs
p_plus = sp.Rational(1, 2) + sp.Rational(1, 2) * ts * X
p_minus = 1 - p_plus
p0_plus = p0s * (1 + ts * as_) / 2 / p_plus
p0_minus = p0s * (1 - ts * as_) / 2 / p_minus

# --- exact martingale property of the measurement step (all orders in t) ---
E_p0_prime = sp.simplify(p_plus * p0_plus + p_minus * p0_minus)
martingale_exact = sp.simplify(E_p0_prime - p0s)
assert martingale_exact == 0, f"measurement step is NOT an exact martingale: {martingale_exact}"

# --- exact fixed point of the mean under measurement + relaxation ---
# E[p0''] = (1-g)*E[p0'] + g*p0_born = (1-g)*p0 + g*p0_born  (exact, any t)
p0_born_s = sp.symbols('p0_born', real=True)
mean_map = (1 - gs) * p0s + gs * p0_born_s
fixed_point_mean = sp.solve(sp.Eq(mean_map, p0s), p0s)
# -> should return [p0_born] for g != 0: bias in the MEAN is exactly zero.

# --- variance injected by one measurement pass, series in t ---
Var_p0_prime = sp.simplify(p_plus * p0_plus**2 + p_minus * p0_minus**2 - p0s**2)
Var_series = sp.series(Var_p0_prime, ts, 0, 4).removeO()
Var_series = sp.expand(Var_series)
# leading (t^2) coefficient, as a function of p0, a, b
t2_coeff = sp.simplify(Var_series.coeff(ts, 2))
t2_coeff_factored = sp.factor(t2_coeff)

# evaluate at p0 = p0_born symbolically-numerically
phi_s = sp.GoldenRatio
p0_born_expr = 1 / phi_s**2
a_val = sp.cos(sp.rad(108))
b_val = sp.cos(sp.rad(36))
t2_coeff_at_born = sp.nsimplify(t2_coeff_factored.subs({p0s: p0_born_expr, as_: a_val, bs: b_val}))
t2_coeff_at_born_float = float(t2_coeff_at_born.evalf())

# steady-state variance of p0 from the balance recursion:
# Var(n+1) = (1-g)^2 * ( Var(n) + f(p0,t) ),  f = Var_p0_prime evaluated near p0_born
# => Var_ss = (1-g)^2 * f / (1 - (1-g)^2)
f_at_born_t2 = t2_coeff_at_born_float  # coefficient of t^2 in f(p0_born, t)
g_sym = sp.symbols('g', positive=True)
t_sym = sp.symbols('t', positive=True)
f_leading = f_at_born_t2 * t_sym**2
Var_ss_expr = (1 - g_sym)**2 * f_leading / (1 - (1 - g_sym)**2)
Var_ss_small_g = sp.series(Var_ss_expr, g_sym, 0, 1).removeO()  # leading small-g behavior
Var_ss_small_g_simplified = sp.simplify(Var_ss_small_g)

# convert to sigma of 360*Re x: Var(360 Re x) = 360^2 (a-b)^2 Var(p0)
ab_diff2 = float((a_val - b_val)**2)
sigma_theory_coeff = 360.0 * abs(float(a_val - b_val))  # multiplies sqrt(Var_ss(p0))

results['R4'] = dict(
    martingale_exact_zero=(martingale_exact == 0),
    fixed_point_mean=[str(x) for x in fixed_point_mean],
    t2_coeff_factored=str(t2_coeff_factored),
    t2_coeff_at_p0_born=t2_coeff_at_born_float,
    Var_ss_small_g_law=str(Var_ss_small_g_simplified),
    sigma_360Rex_coeff=sigma_theory_coeff,
)

# theoretical sigma prediction vs t,g on the R3 grid, for comparison
theory_sigma = {}
for t in T_LIST:
    for g in G_LIST:
        var_ss_p0 = f_at_born_t2 * t**2 * (1 - g)**2 / (1 - (1 - g)**2)
        sigma_pred = sigma_theory_coeff * np.sqrt(var_ss_p0)
        theory_sigma[f"{t}|{g}"] = float(sigma_pred)
results['R4_theory_sigma_grid'] = theory_sigma

# ---------------------------------------------------------------------
# dump raw results for results.md to consume
# ---------------------------------------------------------------------
with open(os.path.join(OUTDIR, 'results_raw.json'), 'w') as fh:
    json.dump(results, fh, indent=2, default=str)

print("DONE")
print("alpha_inv_born =", ALPHA_INV_BORN)
print("martingale exact zero:", martingale_exact == 0)
print("fixed point of mean:", fixed_point_mean)
print("t^2 coeff of Var(p0') at p0_born:", t2_coeff_at_born_float)
print("R3 fits:", results['R3_fits'])
