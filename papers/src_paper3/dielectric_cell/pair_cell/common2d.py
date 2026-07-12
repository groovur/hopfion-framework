import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

PHI = (1.0 + np.sqrt(5.0)) / 2.0


def find_x_star():
    def f(x):
        return (15.0 / 8.0) * np.arctan(x) / x - PHI
    return brentq(f, 1e-6, 10.0, xtol=1e-14, rtol=1e-14)


X_STAR = find_x_star()
BETA_STAR = X_STAR ** 2 / 8.0


def tube_X(rho):
    return 8.0 / (1.0 + rho ** 2) ** 2


def eps_tan_of_X(Xval, beta, eta):
    return 1.0 / (1.0 + beta * Xval) ** 2 + 1j * eta


def eps_rad_of_X(Xval, beta, eta):
    return (1.0 - 3.0 * beta * Xval) / (1.0 + beta * Xval) ** 3 + 1j * eta


def harmonic_mean(a, b):
    return 2.0 * a * b / (a + b)


class Grid:
    """Cell-centered square grid [-L, L]^2 with a 1-cell ghost layer."""

    def __init__(self, L, h):
        self.L = L
        self.h = h
        N = int(round(2 * L / h))
        if N % 2 == 1:
            N += 1  # keep even so no cell center sits exactly at x=0
        self.N = N
        xc = -L + h / 2.0 + np.arange(N) * h
        self.xc = xc
        xc_pad = np.concatenate(([xc[0] - h], xc, [xc[-1] + h]))
        self.xc_pad = xc_pad
        self.Np = N + 2
        Xp, Yp = np.meshgrid(xc_pad, xc_pad, indexing='ij')
        self.Xp = Xp
        self.Yp = Yp


class QuadrantGrid:
    """Cell-centered grid on the quadrant [0, L]^2 with a 1-cell ghost
    layer on every side. The ghost layers at x<0 / y<0 sit at the mirror
    images (-h/2) of the first interior cells and are eliminated through
    mirror symmetry of the solution in assemble_system_quadrant."""

    def __init__(self, L, h):
        self.L = L
        self.h = h
        M = int(round(L / h))
        self.N = M
        xc = h / 2.0 + np.arange(M) * h
        self.xc = xc
        xc_pad = np.concatenate(([-h / 2.0], xc, [xc[-1] + h]))
        self.xc_pad = xc_pad
        self.Np = M + 2
        Xp, Yp = np.meshgrid(xc_pad, xc_pad, indexing='ij')
        self.Xp = Xp
        self.Yp = Yp


def stencil_offsets(grid, exx_p, eyy_p, exy_p):
    """9-point mixed-derivative stencil coefficients for
    div(eps grad u) = 0 on a cell-centered grid with harmonic averaging of
    the diagonal tensor components at faces and arithmetic averaging of the
    off-diagonal component. Returns dict (di,dj) -> (N,N) coefficient
    array (units: flux per unit potential; the common 1/h factor of the
    divergence is dropped since the equation RHS is zero)."""
    N = grid.N
    h = grid.h

    def sl(arr, di, dj):
        # slice of padded array shifted by (di,dj) relative to interior cells
        return arr[1 + di:1 + di + N, 1 + dj:1 + dj + N]

    exx_c = sl(exx_p, 0, 0)
    exx_e = sl(exx_p, 1, 0)
    exx_w = sl(exx_p, -1, 0)
    eyy_c = sl(eyy_p, 0, 0)
    eyy_n = sl(eyy_p, 0, 1)
    eyy_s = sl(eyy_p, 0, -1)
    exy_c = sl(exy_p, 0, 0)
    exy_e = sl(exy_p, 1, 0)
    exy_w = sl(exy_p, -1, 0)
    exy_n = sl(exy_p, 0, 1)
    exy_s = sl(exy_p, 0, -1)

    exx_R = harmonic_mean(exx_c, exx_e)
    exx_L = harmonic_mean(exx_w, exx_c)
    eyy_T = harmonic_mean(eyy_c, eyy_n)
    eyy_B = harmonic_mean(eyy_s, eyy_c)
    exy_R = 0.5 * (exy_c + exy_e)
    exy_L = 0.5 * (exy_w + exy_c)
    exy_T = 0.5 * (exy_c + exy_n)
    exy_B = 0.5 * (exy_s + exy_c)

    inv_h = 1.0 / h
    inv_4h = 1.0 / (4.0 * h)

    C = -exx_R * inv_h - exx_L * inv_h - eyy_T * inv_h - eyy_B * inv_h
    E_ = exx_R * inv_h + exy_T * inv_4h - exy_B * inv_4h
    W_ = exx_L * inv_h - exy_T * inv_4h + exy_B * inv_4h
    Nn = exy_R * inv_4h - exy_L * inv_4h + eyy_T * inv_h
    S_ = -exy_R * inv_4h + exy_L * inv_4h + eyy_B * inv_h
    NE = (exy_R + exy_T) * inv_4h
    SE = -(exy_R + exy_B) * inv_4h
    NW = -(exy_L + exy_T) * inv_4h
    SW = (exy_L + exy_B) * inv_4h

    offsets = {
        (0, 0): C, (1, 0): E_, (-1, 0): W_, (0, 1): Nn, (0, -1): S_,
        (1, 1): NE, (1, -1): SE, (-1, 1): NW, (-1, -1): SW,
    }
    return offsets


def assemble_system(grid, exx_p, eyy_p, exy_p):
    """
    Build the sparse matrix A once from the dielectric tensor (independent
    of the applied field). Returns (A, ghost_info); ghost_info carries the
    boundary-coupling bookkeeping used by build_rhs() to form the RHS for
    any applied field, so one LU factorization of A can be reused for
    multiple field directions.
    """
    N = grid.N
    offsets = stencil_offsets(grid, exx_p, eyy_p, exy_p)

    ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    idx = ii * N + jj

    rows_list = []
    cols_list = []
    data_list = []
    # per-offset ghost bookkeeping, reused later to build the RHS for any
    # applied field without recomputing the tensor/matrix
    ghost_info = []

    for (di, dj), coeff in offsets.items():
        ni = ii + di
        nj = jj + dj
        interior_mask = (ni >= 0) & (ni < N) & (nj >= 0) & (nj < N)
        r_i = idx[interior_mask]
        c_i = (ni[interior_mask]) * N + nj[interior_mask]
        d_i = coeff[interior_mask]
        rows_list.append(r_i)
        cols_list.append(c_i)
        data_list.append(d_i)
        ghost_mask = ~interior_mask
        if np.any(ghost_mask):
            gi = ii[ghost_mask] + di + 1  # padded index
            gj = jj[ghost_mask] + dj + 1
            r_g = idx[ghost_mask]
            d_g = coeff[ghost_mask]
            ghost_info.append((r_g, gi, gj, d_g))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    data = np.concatenate(data_list)

    A = sp.csr_matrix((data, (rows, cols)), shape=(N * N, N * N))
    return A, ghost_info


def assemble_system_quadrant(grid, exx_p, eyy_p, exy_p, parity, Efield):
    """
    Assemble and return (A, b) on a QuadrantGrid (domain [0,L]^2), folding
    the two mirror boundaries (x=0, y=0) into the matrix via the known
    parity of the solution instead of solving them as separate unknowns.

    parity = (sign_x, sign_y): sign_x = u(-x,y)/u(x,y), sign_y =
    u(x,-y)/u(x,y). For the pair problem (mirror-symmetric tensor field
    about both axes): E along the pair axis (x) gives parity (-1,+1)
    ("par": odd in x, even in y); E transverse (y) gives parity (+1,-1)
    ("perp": even in x, odd in y).

    Only the innermost ghost ring (mirror side) maps back into the domain;
    the outer ghost ring (x=L, y=L) is still ordinary Dirichlet -> RHS.
    """
    N = grid.N
    sign_x, sign_y = parity
    offsets = stencil_offsets(grid, exx_p, eyy_p, exy_p)

    ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    idx = ii * N + jj
    Ex, Ey = Efield
    Ubc = -(Ex * grid.Xp + Ey * grid.Yp)
    b = np.zeros(N * N, dtype=complex)

    rows_list = []
    cols_list = []
    data_list = []

    for (di, dj), coeff in offsets.items():
        ni = ii + di
        nj = jj + dj
        mirror_x = (ni == -1)
        mirror_y = (nj == -1)
        real_bd = (ni == N) | (nj == N)

        interior = (~mirror_x) & (~mirror_y) & (~real_bd)
        r_i = idx[interior]
        c_i = ni[interior] * N + nj[interior]
        d_i = coeff[interior]
        rows_list.append(r_i)
        cols_list.append(c_i)
        data_list.append(d_i)

        only_x = mirror_x & (~mirror_y) & (~real_bd)
        r_mx = idx[only_x]
        c_mx = 0 * N + nj[only_x]
        d_mx = coeff[only_x] * sign_x
        rows_list.append(r_mx)
        cols_list.append(c_mx)
        data_list.append(d_mx)

        only_y = mirror_y & (~mirror_x) & (~real_bd)
        r_my = idx[only_y]
        c_my = ni[only_y] * N + 0
        d_my = coeff[only_y] * sign_y
        rows_list.append(r_my)
        cols_list.append(c_my)
        data_list.append(d_my)

        both = mirror_x & mirror_y
        r_b = idx[both]
        c_b = np.zeros(r_b.shape, dtype=int)
        d_b = coeff[both] * sign_x * sign_y
        rows_list.append(r_b)
        cols_list.append(c_b)
        data_list.append(d_b)

        if np.any(real_bd):
            gi = ii[real_bd] + di + 1
            gj = jj[real_bd] + dj + 1
            ubc_vals = Ubc[gi, gj]
            r_rb = idx[real_bd]
            d_rb = coeff[real_bd]
            np.add.at(b, r_rb, -d_rb * ubc_vals)

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    data = np.concatenate(data_list)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N * N, N * N))
    return A, b


def solve_quadrant(grid, exx_p, eyy_p, exy_p, parity, Efield):
    A, b = assemble_system_quadrant(grid, exx_p, eyy_p, exy_p, parity, Efield)
    u_flat = spla.spsolve(A, b)
    return u_flat.reshape(grid.N, grid.N)


def build_rhs(grid, ghost_info):
    """Return a function Efield -> b (RHS vector) reusing the ghost
    bookkeeping from assemble_system, without rebuilding the matrix."""
    N = grid.N
    Xp, Yp = grid.Xp, grid.Yp

    def make_b(Efield):
        Ex, Ey = Efield
        Ubc = -(Ex * Xp + Ey * Yp)
        b = np.zeros(N * N, dtype=complex)
        for (r_g, gi, gj, d_g) in ghost_info:
            ubc_vals = Ubc[gi, gj]
            np.add.at(b, r_g, -d_g * ubc_vals)
        return b

    return make_b


def solve_system(grid, A, ghost_info, Efield):
    """Single-RHS solve (spsolve). For multiple Efield directions on the
    same tensor, prefer solve_system_multi to reuse one LU factorization."""
    make_b = build_rhs(grid, ghost_info)
    b = make_b(Efield)
    u_flat = spla.spsolve(A, b)
    return u_flat.reshape(grid.N, grid.N)


def solve_system_multi(grid, A, ghost_info, Efields):
    """Factorize A once (splu) and solve for each Efield in Efields.
    Returns a list of (N,N) complex solutions, one per Efield."""
    make_b = build_rhs(grid, ghost_info)
    lu = spla.splu(A.tocsc())
    us = []
    for Efield in Efields:
        b = make_b(Efield)
        u_flat = lu.solve(b)
        us.append(u_flat.reshape(grid.N, grid.N))
    return us


def assemble_and_solve(grid, exx_p, eyy_p, exy_p, Efield, verbose=False):
    """Convenience single-field-direction solve (assembles + factorizes +
    solves in one call)."""
    A, ghost_info = assemble_system(grid, exx_p, eyy_p, exy_p)
    return solve_system(grid, A, ghost_info, Efield)


def bilinear_interp(grid, field, xq, yq):
    """Bilinear interpolation of a real (N,N) array `field` at points (xq,yq).
    xq, yq: 1D arrays. Uses interior grid grid.xc (uniform spacing h)."""
    h = grid.h
    xc = grid.xc
    N = grid.N
    fi = (xq - xc[0]) / h
    fj = (yq - xc[0]) / h
    i0 = np.floor(fi).astype(int)
    j0 = np.floor(fj).astype(int)
    i0 = np.clip(i0, 0, N - 2)
    j0 = np.clip(j0, 0, N - 2)
    ti = fi - i0
    tj = fj - j0
    f00 = field[i0, j0]
    f10 = field[i0 + 1, j0]
    f01 = field[i0, j0 + 1]
    f11 = field[i0 + 1, j0 + 1]
    return (f00 * (1 - ti) * (1 - tj) + f10 * ti * (1 - tj) +
            f01 * (1 - ti) * tj + f11 * ti * tj)


def tensor_from_scalar(eps_val):
    """Isotropic tensor from a scalar (possibly complex) field eps_val
    (array shape (Np,Np)): exx=eyy=eps_val, exy=0."""
    exx = eps_val.astype(complex)
    eyy = eps_val.astype(complex)
    exy = np.zeros_like(exx)
    return exx, eyy, exy


def tensor_from_g(Xscalar, gx, gy, beta, eta):
    """Build the anisotropic tensor eps_ij = et*delta_ij + (er-et)*ghat_i*ghat_j
    from the scalar background field Xscalar (used in et,er) and the
    (possibly unnormalized) gradient-direction vector field (gx,gy).
    Where |g| ~ 0, ghat is set to (0,0), which makes the tensor isotropic
    (eps_ij = et*delta_ij) there -- exactly the "use isotropic eps_tan"
    fallback."""
    et = eps_tan_of_X(Xscalar, beta, eta)
    er = eps_rad_of_X(Xscalar, beta, eta)
    gmag = np.sqrt(gx ** 2 + gy ** 2)
    safe = gmag > 1e-12
    ghx = np.where(safe, gx / np.where(safe, gmag, 1.0), 0.0)
    ghy = np.where(safe, gy / np.where(safe, gmag, 1.0), 0.0)
    diff = (er - et)
    exx = et + diff * ghx ** 2
    eyy = et + diff * ghy ** 2
    exy = diff * ghx * ghy
    return exx, eyy, exy


def single_tube_tensor(Xp, Yp, beta, eta, center=(0.0, 0.0)):
    """Single k-essence tube centered at `center`; g is radial (unit vector
    pointing away from the tube center)."""
    dx = Xp - center[0]
    dy = Yp - center[1]
    rho = np.sqrt(dx ** 2 + dy ** 2)
    Xval = tube_X(rho)
    return tensor_from_g(Xval, dx, dy, beta, eta)


def pair_tensor(Xp, Yp, beta, eta, s, model):
    """Two identical tubes at (+-s/2, 0). model = 'A' (incoherent, X_A=X1+X2)
    or 'B' (coherent gradient, X_B=|g|^2, cross term included)."""
    c1 = (-s / 2.0, 0.0)
    c2 = (s / 2.0, 0.0)
    dx1 = Xp - c1[0]
    dy1 = Yp - c1[1]
    rho1 = np.sqrt(dx1 ** 2 + dy1 ** 2)
    dx2 = Xp - c2[0]
    dy2 = Yp - c2[1]
    rho2 = np.sqrt(dx2 ** 2 + dy2 ** 2)
    X1 = tube_X(rho1)
    X2 = tube_X(rho2)
    # unit radial vectors from each tube (guard rho=0, not hit on our grid
    # since N is kept even so no cell center sits at x=0, but tube centers
    # are offset from the origin as well for s>0, so this is just a safety
    # floor)
    rho1s = np.where(rho1 > 1e-12, rho1, 1e-12)
    rho2s = np.where(rho2 > 1e-12, rho2, 1e-12)
    r1x, r1y = dx1 / rho1s, dy1 / rho1s
    r2x, r2y = dx2 / rho2s, dy2 / rho2s
    sq1 = np.sqrt(X1)
    sq2 = np.sqrt(X2)
    gx = sq1 * r1x + sq2 * r2x
    gy = sq1 * r1y + sq2 * r2y
    if model == 'A':
        Xscalar = X1 + X2
    elif model == 'B':
        Xscalar = gx ** 2 + gy ** 2
    else:
        raise ValueError("model must be 'A' or 'B'")
    return tensor_from_g(Xscalar, gx, gy, beta, eta)


def extract_dipole_quadrant(grid, u, Efield, component, r_inner=15.0,
                             r_outer=25.0, n_radii=10, n_theta=90):
    """Dipole extraction for a QuadrantGrid solution. `component` is 'cos'
    (E along the pair/x axis: mirror symmetry forces the sin-theta
    coefficient to vanish exactly, so we fit only the cos-theta channel)
    or 'sin' (E transverse/y axis: cos-theta coefficient vanishes).
    Sampling is restricted to the open first-quadrant angular range
    (0, pi/2) since the domain only covers x,y > 0. Returns the complex
    dipole coefficient (px for 'cos', py for 'sin')."""
    Ex, Ey = Efield
    radii = np.linspace(r_inner, r_outer, n_radii)
    thetas = np.linspace(1e-3, np.pi / 2.0 - 1e-3, n_theta)
    R, TH = np.meshgrid(radii, thetas, indexing='ij')
    R = R.ravel()
    TH = TH.ravel()
    Xq = R * np.cos(TH)
    Yq = R * np.sin(TH)

    ur = bilinear_interp(grid, u.real, Xq, Yq)
    ui = bilinear_interp(grid, u.imag, Xq, Yq)
    uq = ur + 1j * ui
    resid = uq + (Ex * Xq + Ey * Yq)

    ang = np.cos(TH) if component == 'cos' else np.sin(TH)
    invr = 1.0 / R
    invr3 = 1.0 / R ** 3
    basis = np.stack([ang * invr, ang * invr3, ang * R], axis=1)
    coef, *_ = np.linalg.lstsq(basis, resid, rcond=None)
    return coef[0]


def extract_dipole(grid, u, Efield, r_inner=15.0, r_outer=25.0,
                    n_radii=6, n_theta=72):
    """Fit u + E.r on a ring r in [r_inner,r_outer] to (px cos + py sin)/r
    via least squares, returning complex (px, py).

    A finite-L Dirichlet boundary (u = -E.r imposed exactly at r=L, rather
    than the true u = -E.r + p.r/r^2) slightly shifts the effective field
    amplitude inside the domain away from the nominal E. This produces a
    small *r-growing* admixture (E_eff - E).r in the residual u + E.r,
    which competes with the decaying dipole term p.rhat/r at large r and
    corrupts a naive fit that assumes only the decaying piece is present.
    We therefore fit both pieces simultaneously (decaying 1/r, 1/r^3 dipole
    terms plus a growing r term absorbing the finite-L leakage) and read
    off the coefficient of the 1/r term as the true dipole moment.
    """
    Ex, Ey = Efield
    radii = np.linspace(r_inner, r_outer, n_radii)
    thetas = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    R, TH = np.meshgrid(radii, thetas, indexing='ij')
    R = R.ravel()
    TH = TH.ravel()
    Xq = R * np.cos(TH)
    Yq = R * np.sin(TH)

    ur = bilinear_interp(grid, u.real, Xq, Yq)
    ui = bilinear_interp(grid, u.imag, Xq, Yq)
    uq = ur + 1j * ui

    resid = uq + (Ex * Xq + Ey * Yq)

    cB = np.cos(TH)
    sB = np.sin(TH)
    invr = 1.0 / R
    invr3 = 1.0 / R ** 3
    basis = np.stack([cB * invr, sB * invr, cB * invr3, sB * invr3,
                       cB * R, sB * R], axis=1)
    coef, *_ = np.linalg.lstsq(basis, resid, rcond=None)
    px, py = coef[0], coef[1]
    return px, py
