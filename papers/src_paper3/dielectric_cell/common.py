import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2


def X(rho):
    return 8.0 / (1.0 + rho**2)**2


def find_x_star():
    def f(x):
        return (15.0/8.0) * np.arctan(x)/x - PHI
    # f(0+) = 15/8 - phi > 0 ; find where it crosses zero
    x_star = brentq(f, 1e-6, 10.0, xtol=1e-14, rtol=1e-14)
    return x_star


X_STAR = find_x_star()
BETA_STAR = X_STAR**2 / 8.0


def rho_c_of_beta(beta):
    # (1+rho_c^2)^2 = 24*beta  => rho_c = sqrt(sqrt(24*beta) - 1)
    val = 24.0 * beta
    if val < 1.0:
        return None  # no negative-eps region
    return np.sqrt(np.sqrt(val) - 1.0)


def eps_tan_bg(rho, beta):
    Xr = X(rho)
    return 1.0 / (1.0 + beta*Xr)**2


def eps_rad_bg(rho, beta):
    Xr = X(rho)
    return (1.0 - 3.0*beta*Xr) / (1.0 + beta*Xr)**3


def eps_tan_step(rho, eps1, R):
    return np.where(rho < R, eps1, 1.0)


def eps_rad_step(rho, eps1, R):
    return np.where(rho < R, eps1, 1.0)


def rhs_factory(eps_rad_fn, eps_tan_fn, eta):
    """Return f(rho, y) for y=[w, v], v = rho*eps_rad*w'.
    eps_rad_fn, eps_tan_fn: functions of rho (real-valued, background physics)
    eta: regularization added as +i*eta to both eps_rad and eps_tan.
    """
    def f(rho, y):
        w, v = y
        er = eps_rad_fn(rho) + 1j*eta
        et = eps_tan_fn(rho) + 1j*eta
        dw = v / (rho * er)
        dv = et * w / rho
        return np.array([dw, dv], dtype=complex)
    return f


def shoot_alpha(eps_rad_fn, eps_tan_fn, eta, rho_min=1e-4, rho_max=200.0,
                 rtol=1e-10, atol=1e-12, dense=False):
    """Shooting solve. Returns alpha_pol (complex), and solver info dict.

    Regular IC at rho_min: treating eps_rad, eps_tan as locally constant near
    the origin, the indicial (Frobenius) exponent for w ~ rho^s solutions of
    the ODE is s = +/- sqrt(eps_tan/eps_rad). When eps_tan == eps_rad (the
    isotropic case used in V1-V3) this reduces to s=1, i.e. w ~ A*rho exactly
    as stated. In the anisotropic background (R1/R2), and especially where
    eps_rad(0) < 0 (origin inside the negative-eps shell), s != 1 in general;
    using the literal w=rho, w'=1 IC there excites a spurious admixture of
    the irregular (growing-at-origin) branch that does not vanish as
    rho_min -> 0 (confirmed numerically: naive IC gives non-convergent
    rho_min sensitivity). We therefore use the correct local regular
    exponent s (principal branch, Re(s) >= 0) built from the *regularized*
    (eta-shifted) local eps values, which is the literal same physical
    regularity condition (finite, single-valued potential at the origin)
    applied consistently with the eta -> 0 limit-absorption prescription.
    """
    f = rhs_factory(eps_rad_fn, eps_tan_fn, eta)
    er0 = eps_rad_fn(rho_min) + 1j*eta
    et0 = eps_tan_fn(rho_min) + 1j*eta
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
    erR = eps_rad_fn(rho_max) + 1j*eta
    wpR = vR / (rho_max * erR)
    R = rho_max
    c1 = (wR + R*wpR) / (2.0*R)
    c2 = R*(wR - R*wpR) / 2.0
    alpha = -c2 / c1
    info = {'success': True, 'nfev': sol.nfev, 'c1': c1, 'c2': c2,
            'w_end': wR, 'wp_end': wpR}
    return alpha, info


def bg_eps_fns(beta):
    return (lambda rho: eps_rad_bg(rho, beta), lambda rho: eps_tan_bg(rho, beta))


def step_eps_fns(eps1, R):
    return (lambda rho: eps_rad_step(rho, eps1, R), lambda rho: eps_tan_step(rho, eps1, R))
