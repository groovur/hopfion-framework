"""Grid construction and the per-grid-point dielectric tensor field.

Periodic images of every tube curve are built explicitly (not via cKDTree's
internal minimum-image wrap) so that the literal spec cutoff (10 units) is
honored even though the fallback unit cell is smaller than 2*cutoff in y,z:
all periodic shells needed to cover every image within `cutoff` of any
point in the base cell are included in each tube's search tree.
"""
import numpy as np
from scipy.spatial import cKDTree

from common3d import u_profile, eps_tan_of_U, eps_rad_of_U


class Grid:
    """Cell-centered grid on [0,Lx) x [0,Ly) x [0,Lz), uniform spacing h in
    all three axes (Nx = Lx/h etc., generally Nx=2*Ny=2*Nz for the fallback
    cell L=(2*a_c, a_c, a_c))."""

    def __init__(self, L, h):
        self.L = np.asarray(L, dtype=float)
        self.h = h
        self.N = tuple(int(round(Li / h)) for Li in L)
        xc = [ (np.arange(Ni) + 0.5) * h for Ni in self.N ]
        self.xc = xc
        Xc, Yc, Zc = np.meshgrid(xc[0], xc[1], xc[2], indexing='ij')
        self.pts = np.stack([Xc.ravel(), Yc.ravel(), Zc.ravel()], axis=1)
        self.shape = self.N


def periodic_images(pts, L, cutoff):
    """Return all periodic replicas of `pts` (n,3) needed so that every
    image within `cutoff` of any point of the base cell [0,L) is present."""
    Lx, Ly, Lz = L
    shells = [int(np.ceil(cutoff / Li)) for Li in (Lx, Ly, Lz)]
    out = []
    for si in range(-shells[0], shells[0] + 1):
        for sj in range(-shells[1], shells[1] + 1):
            for sk in range(-shells[2], shells[2] + 1):
                shift = np.array([si * Lx, sj * Ly, sk * Lz])
                out.append(pts + shift)
    return np.concatenate(out, axis=0)


def tube_field_contribution(grid_pts, curve_pts, L, cutoff, eta):
    """Nearest-point distance & unit direction from every grid point to one
    tube curve (via its periodic images within cutoff). Returns (rho, rhat)
    with rho = np.inf and rhat = 0 where nothing is within cutoff."""
    images = periodic_images(curve_pts, L, cutoff)
    tree = cKDTree(images)
    rho, idx = tree.query(grid_pts, k=1, distance_upper_bound=cutoff, workers=-1)
    finite = np.isfinite(rho)
    rhat = np.zeros_like(grid_pts)
    nearest = np.zeros_like(grid_pts)
    nearest[finite] = images[idx[finite]]
    disp = grid_pts - nearest
    rho_safe = np.where(finite & (rho > 1e-12), rho, 1.0)
    rhat[finite] = disp[finite] / rho_safe[finite, None]
    return rho, rhat, finite


def build_U_g(grid, curves, cutoff=10.0):
    """U = sum_i u_i(rho_i); g = sum_i sqrt(u_i) rhat_i, over all tube
    curves in `curves` (each a dict with key 'pts')."""
    n = grid.pts.shape[0]
    U = np.zeros(n)
    g = np.zeros((n, 3))
    for c in curves:
        rho, rhat, finite = tube_field_contribution(grid.pts, c['pts'], grid.L, cutoff, 0.0)
        u_i = np.zeros(n)
        u_i[finite] = u_profile(rho[finite])
        U += u_i
        g += np.sqrt(u_i)[:, None] * rhat
    return U, g


def build_eps_tensor(U, g, eta):
    """eps_jk = eps_tan*delta_jk + (eps_rad-eps_tan)*ghat_j*ghat_k, plus the
    volume averages S1=<1/(1+U)> and Ubar=<U> needed for the vertex
    factor."""
    n = U.shape[0]
    gmag = np.linalg.norm(g, axis=1)
    safe = gmag > 1e-12
    ghat = np.zeros_like(g)
    ghat[safe] = g[safe] / gmag[safe, None]

    eps_tan = eps_tan_of_U(U, eta)
    eps_rad = eps_rad_of_U(U, eta)
    diff = eps_rad - eps_tan

    # symmetric tensor stored as 6 components: xx,yy,zz,xy,xz,yz
    exx = eps_tan + diff * ghat[:, 0] ** 2
    eyy = eps_tan + diff * ghat[:, 1] ** 2
    ezz = eps_tan + diff * ghat[:, 2] ** 2
    exy = diff * ghat[:, 0] * ghat[:, 1]
    exz = diff * ghat[:, 0] * ghat[:, 2]
    eyz = diff * ghat[:, 1] * ghat[:, 2]

    S1 = np.mean(1.0 / (1.0 + U))
    Ubar = np.mean(U)
    return (exx, eyy, ezz, exy, exz, eyz), S1, Ubar
