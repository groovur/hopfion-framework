"""Matrix-free FFT homogenization for the periodic cell problem
div(eps grad u) = 0 with imposed mean field E (Moulinec-Suquet Lippmann-
Schwinger formulation, complex reference medium eps_ref).

We solve the fixed-point equation
    e + Gamma0[(eps - eps_ref) e] = E        (*)
(e = E + grad u, the total local field; Gamma0 the periodic Green operator
of the homogeneous eps_ref medium) as a linear system via matrix-free GMRES
rather than plain Picard iteration on (*). This is a deliberate departure
from literal Eyre-Milton fixed-point coefficients (not reproduced here from
memory with confidence) that keeps the same Green operator and the same
Eyre-Milton-optimal complex reference medium eps_ref = sqrt(eps_min*eps_max),
but replaces the fixed-point sweep with a Krylov solve for robustness across
the sign-changing-real-part / eta-regularized contrast regions where plain
Picard iteration is not guaranteed to contract. Documented as such; not
literally "the Eyre-Milton scheme."
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


def _fftn(a):
    return np.fft.fftn(a)


def _ifftn(a):
    return np.fft.ifftn(a)


class Gamma0:
    def __init__(self, shape, h):
        Nx, Ny, Nz = shape
        kx = 2.0 * np.pi * np.fft.fftfreq(Nx, d=h)
        ky = 2.0 * np.pi * np.fft.fftfreq(Ny, d=h)
        kz = 2.0 * np.pi * np.fft.fftfreq(Nz, d=h)
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        self.KX, self.KY, self.KZ = KX, KY, KZ
        K2 = KX**2 + KY**2 + KZ**2
        K2[0, 0, 0] = 1.0
        self.K2 = K2
        self.shape = shape

    def apply(self, tau, eps_ref):
        """tau: (3,Nx,Ny,Nz) complex. Returns Gamma0(tau), same shape."""
        t0h = _fftn(tau[0])
        t1h = _fftn(tau[1])
        t2h = _fftn(tau[2])
        kdott = self.KX * t0h + self.KY * t1h + self.KZ * t2h
        factor = kdott / (eps_ref * self.K2)
        factor[0, 0, 0] = 0.0
        out = np.empty_like(tau)
        out[0] = _ifftn(self.KX * factor)
        out[1] = _ifftn(self.KY * factor)
        out[2] = _ifftn(self.KZ * factor)
        return out


def tensor_apply(eps6, e):
    """eps6 = (exx,eyy,ezz,exy,exz,eyz), each (Nx,Ny,Nz) complex.
    e: (3,Nx,Ny,Nz). Returns eps . e, same shape."""
    exx, eyy, ezz, exy, exz, eyz = eps6
    ex, ey, ez = e[0], e[1], e[2]
    out = np.empty_like(e)
    out[0] = exx * ex + exy * ey + exz * ez
    out[1] = exy * ex + eyy * ey + eyz * ez
    out[2] = exz * ex + eyz * ey + ezz * ez
    return out


def eyre_milton_eps_ref(eps6):
    """sqrt(eps_min*eps_max) using the diagonal-component real-part extremes
    (Eyre-Milton optimal reference medium)."""
    exx, eyy, ezz = eps6[0], eps6[1], eps6[2]
    diag_re = np.concatenate([exx.real.ravel(), eyy.real.ravel(), ezz.real.ravel()])
    eps_min = diag_re.min()
    eps_max = diag_re.max()
    diag_im = np.concatenate([exx.imag.ravel(), eyy.imag.ravel(), ezz.imag.ravel()])
    im = diag_im.mean()
    return np.sqrt(complex(eps_min, im) * complex(eps_max, im))


def solve_cell(eps6, h, Efield, eta, tol=1e-7, maxiter=2000, restart=40,
               eps_ref=None, callback_report=None):
    """Solve for the total local field e (3,Nx,Ny,Nz) given the tensor
    field eps6 and imposed mean field Efield (3,)."""
    shape = eps6[0].shape
    n = int(np.prod(shape))
    G = Gamma0(shape, h)
    if eps_ref is None:
        eps_ref = eyre_milton_eps_ref(eps6)

    E_uniform = np.zeros((3,) + shape, dtype=complex)
    for i in range(3):
        E_uniform[i] = Efield[i]

    I3 = [np.ones(shape, dtype=complex) * 0.0 for _ in range(1)]

    def matvec(x):
        e = x.reshape((3,) + shape)
        tau = tensor_apply(eps6, e)
        # subtract eps_ref * e
        tau = tau - eps_ref * e
        Ge = G.apply(tau, eps_ref)
        out = e + Ge
        return out.ravel()

    A = LinearOperator((3 * n, 3 * n), matvec=matvec, dtype=complex)
    resids = []

    def cb(rk):
        resids.append(rk)

    x0 = E_uniform.ravel()
    b = E_uniform.ravel()
    sol, info = gmres(A, b, x0=x0, rtol=tol, atol=0.0, restart=restart,
                       maxiter=maxiter, callback=cb, callback_type='pr_norm')
    e = sol.reshape((3,) + shape)
    n_iter = len(resids)
    final_resid = resids[-1] if resids else None
    if callback_report is not None:
        callback_report(info, n_iter, final_resid)
    return e, {'info': info, 'n_iter': n_iter, 'final_resid': final_resid,
               'eps_ref': eps_ref}


def eps_eff_tensor(eps6, h, eta, tol=1e-7, maxiter=2000, restart=40,
                    eps_ref=None, report=print):
    """Solve 3 field directions (x,y,z) and assemble the full eps_eff
    tensor from eps_eff[:,l] = <eps . e_l> (mean flux for unit mean field
    E=e_l)."""
    shape = eps6[0].shape
    n = float(np.prod(shape))
    eps_eff = np.zeros((3, 3), dtype=complex)
    solve_info = []
    if eps_ref is None:
        eps_ref = eyre_milton_eps_ref(eps6)
    for l, Efield in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
        e, info = solve_cell(eps6, h, Efield, eta, tol=tol, maxiter=maxiter,
                              restart=restart, eps_ref=eps_ref)
        D = tensor_apply(eps6, e)
        Dmean = np.array([D[0].sum(), D[1].sum(), D[2].sum()]) / n
        eps_eff[:, l] = Dmean
        info['Efield'] = Efield
        solve_info.append(info)
        report(f"    solve E={Efield}: gmres info={info['info']} "
               f"n_iter={info['n_iter']} final_resid={info['final_resid']} "
               f"eps_ref={eps_ref:.6f}")
    return eps_eff, solve_info
