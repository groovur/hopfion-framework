"""V1-V3 solver validations (must pass before any physics run)."""
import time
import numpy as np

from common3d import CSTAR2, ALPHA_V3, u_profile, eps_tan_of_U, eps_rad_of_U
from solver import eps_eff_tensor


def v1_uniform(report=print):
    shape = (16, 8, 8)
    val = 2.5 + 0.3j
    eps6 = tuple(np.full(shape, val, dtype=complex) if i < 3 else
                 np.zeros(shape, dtype=complex) for i in range(6))
    eps_eff, info = eps_eff_tensor(eps6, h=0.5, eta=0.0, tol=1e-8, maxiter=500, report=lambda *a: None)
    err = np.abs(eps_eff - val * np.eye(3)).max()
    report(f"V1 uniform medium: eps={val}, max|eps_eff-eps*I|={err:.3e}  "
           f"PASS={err < 1e-6}")
    return err < 1e-6


def v2_dilute_spheres(report=print):
    L = 6.0
    M = 48
    h = L / M
    xc = (np.arange(M) + 0.5) * h
    X, Y, Z = np.meshgrid(xc, xc, xc, indexing='ij')
    c = L / 2.0
    r = np.sqrt((X - c)**2 + (Y - c)**2 + (Z - c)**2)
    eps1 = 3.0
    eps_scalar = np.where(r < 1.0, eps1, 1.0).astype(complex)
    zero = np.zeros_like(eps_scalar)
    eps6 = (eps_scalar, eps_scalar, eps_scalar, zero, zero, zero)
    eps_eff, info = eps_eff_tensor(eps6, h=h, eta=0.0, tol=1e-8, maxiter=1000, report=lambda *a: None)
    f = 4.0 * np.pi / 3.0 / L**3
    cm = 1.0 + 3.0 * f * (eps1 - 1.0) / (eps1 + 2.0 - f * (eps1 - 1.0))
    iso = np.trace(eps_eff).real / 3.0
    dev = 100.0 * (iso - cm) / cm
    report(f"V2 dilute spheres: f={f:.5f}  eps_eff_iso={iso:.6f}  CM={cm:.6f}  "
           f"dev={dev:.4f}%  PASS={abs(dev) < 1.0}")
    return abs(dev) < 1.0


def v3_parallel_tube(report=print):
    L = 18.0
    Mxy = 64
    Nz = 4
    h = L / Mxy
    xc = (np.arange(Mxy) + 0.5) * h
    zc = (np.arange(Nz) + 0.5) * h
    X, Y, Z = np.meshgrid(xc, xc, zc, indexing='ij')
    cx = cy = L / 2.0

    def wrap(d, Lp):
        return d - Lp * np.round(d / Lp)

    dx = wrap(X - cx, L)
    dy = wrap(Y - cy, L)
    rho = np.sqrt(dx**2 + dy**2)
    eta = 1e-2
    u = u_profile(rho)
    U = u
    safe = rho > 1e-12
    gx = np.sqrt(u) * np.where(safe, dx / np.where(safe, rho, 1.0), 0.0)
    gy = np.sqrt(u) * np.where(safe, dy / np.where(safe, rho, 1.0), 0.0)
    gmag = np.sqrt(gx**2 + gy**2)
    gsafe = gmag > 1e-12
    ghx = np.where(gsafe, gx / np.where(gsafe, gmag, 1.0), 0.0)
    ghy = np.where(gsafe, gy / np.where(gsafe, gmag, 1.0), 0.0)
    et = eps_tan_of_U(U, eta)
    er = eps_rad_of_U(U, eta)
    diff = er - et
    exx = et + diff * ghx**2
    eyy = et + diff * ghy**2
    ezz = et * np.ones_like(et)
    exy = diff * ghx * ghy
    zero = np.zeros_like(et)
    eps6 = (exx, eyy, ezz, exy, zero, zero)
    eps_eff, info = eps_eff_tensor(eps6, h=h, eta=eta, tol=1e-7, maxiter=1000, report=lambda *a: None)
    n = CSTAR2 / (L**2)
    alpha = ALPHA_V3
    mg = 1.0 + 2.0 * n * alpha / (1.0 - n * alpha)
    dev = 100.0 * (eps_eff[0, 0].real - mg.real) / mg.real
    report(f"V3 parallel tube: n={n:.6f}  alpha={alpha}  MG={mg:.6f}  "
           f"solver eps_eff_xx={eps_eff[0,0]:.6f}  dev={dev:.3f}%  "
           f"PASS(<=10%)={abs(dev) < 10.0}")
    return abs(dev) < 10.0


if __name__ == "__main__":
    t0 = time.time()
    ok1 = v1_uniform()
    ok2 = v2_dilute_spheres()
    ok3 = v3_parallel_tube()
    print(f"\nAll pass: {ok1 and ok2 and ok3}   (total time {time.time()-t0:.1f}s)")
