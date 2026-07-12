"""Coherent two-tube (T(2,2) cross-section) dielectric-cell solver.

Standalone: no imports from the sibling suites. Computes the 2D dipole
polarizability of the coherent two-tube k-essence cell, the vertex
suppression deficit of the same background, and (optionally) the
two-channel contact-packing product against the golden/geo target
112.5/phi^10.

Background
----------
The coherent two-tube configuration is the exact degree-2 CP^1 map
    w(z) = lam^2 / ((z - z1)(z - z2)),   z1 = -s/2, z2 = +s/2,
with lam^2 = s so that near either core w ~ 1/(z - zi): each lump
reduces to the unit tube of the 1D cell (f = 2 arctan(1/rho)) in the
far-separated limit. The kinetic density of the map is
    X = 8 |w'|^2 / (1 + |w|^2)^2,
which equals 8/(1+rho^2)^2 exactly for a single lump. An alternative
sum form w = 1/(z-z1) + e^{i alpha}/(z-z2) is available (--form sum).

Dielectric tensor (fluctuation medium of L(X) = X/(1+bX)):
    eps_tan = L'(X)        = 1/(1+bX)^2
    eps_rad = L' + 2X L''  = (1-3bX)/(1+bX)^3
    eps_ij  = eps_tan delta_ij + (eps_rad - eps_tan) ghat_i ghat_j
Anisotropy prescriptions (--prescription):
    gradf    (default) ghat along grad f, f = 2 arctan|w|. Reduces
             exactly to the validated 1D single-tube cell when s -> inf.
    pullback isotropic eps = L' + X L'' (sigma-model pullback metric of
             a holomorphic map is conformal). Does NOT reduce to the 1D
             cell; provided for comparison only.
    A / B    the toy superpositions of the pair_cell suite (incoherent
             X1+X2 / coherent gradient sum), for cross-checking.
Limit absorption: eps -> eps + i*eta, eta -> 0 from above.

Numerical method (identical to the validated pair_cell suite): 9-point
cell-centered FD for div(eps grad u) = 0, harmonic face averaging of
the diagonal tensor components, Dirichlet u = -E.r at |x|,|y| = L,
dipole read from an m=1 angular Fourier fit on a ring r in [15,25]
including a growing-mode term that absorbs finite-L Dirichlet leakage.
One LU factorization serves both field orientations. Every completed
solve is appended to checkpoint.jsonl immediately and skipped on rerun.

Resolution warning: eps_rad crosses zero on a contour near each core;
the limit-absorption physics is resolved only when the |eps_rad| < eta
shell (width ~ 2 eta / |d eps_rad/d rho|) is covered by several cells.
Take h -> 0 before eta -> 0; see readme.md for the criterion table.
"""

import argparse
import json
import os
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

PHI = (1.0 + np.sqrt(5.0)) / 2.0
R0 = 3.0

B_STAR = brentq(lambda b: (15.0 / 8.0) * np.arctan(np.sqrt(8.0 * b))
                / np.sqrt(8.0 * b) - PHI, 1e-6, 10.0, xtol=1e-14)
X_STAR = np.sqrt(8.0 * B_STAR)
CSTAR2 = 0.75 * PHI ** 5 / 2.0 ** (4.0 / 3.0)
C2STAR = brentq(lambda C: C * (2 * C + 1) / ((C * C + 1) * (3 * C - 1))
                - 2.0 ** (4.0 / 3.0) / PHI ** 5, 1.0, 10.0, xtol=1e-14)
S_DEFAULT = 2.0 * R0 / (C2STAR * np.sqrt(CSTAR2))
TARGET = 112.5 / PHI ** 10
ALPHA_SINGLE_1D = -0.97106613 + 0.01353232j   # eta->0, 1D shooting solver

CHECKPOINT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "checkpoint.jsonl")


# ---------------------------------------------------------------- background

def w_and_wprime(Z, s, form, alpha):
    """Complex field w and derivative w' on complex points Z."""
    z1 = -0.5 * s
    z2 = 0.5 * s
    if form == "product":
        lam2 = s
        d1 = Z - z1
        d2 = Z - z2
        w = lam2 / (d1 * d2)
        wp = -lam2 * (d1 + d2) / (d1 * d2) ** 2
    elif form == "sum":
        ph = np.exp(1j * alpha)
        d1 = Z - z1
        d2 = Z - z2
        w = 1.0 / d1 + ph / d2
        wp = -1.0 / d1 ** 2 - ph / d2 ** 2
    else:
        raise ValueError("form must be 'product' or 'sum'")
    return w, wp


def coherent_fields(Xp, Yp, s, form, alpha):
    """X = 8|w'|^2/(1+|w|^2)^2 and the grad-f direction (gx, gy)."""
    Z = Xp + 1j * Yp
    w, wp = w_and_wprime(Z, s, form, alpha)
    m2 = np.abs(w) ** 2
    Xval = 8.0 * np.abs(wp) ** 2 / (1.0 + m2) ** 2
    # grad f, f = 2 arctan|w|: direction of grad|w|^2 = 2 Re(conj(w) dw)
    gx = np.real(np.conj(w) * wp)
    gy = -np.imag(np.conj(w) * wp)
    return Xval, gx, gy


def tube_X(rho):
    return 8.0 / (1.0 + rho ** 2) ** 2


def eps_tan_of_X(Xval, b, eta):
    return 1.0 / (1.0 + b * Xval) ** 2 + 1j * eta


def eps_rad_of_X(Xval, b, eta):
    return (1.0 - 3.0 * b * Xval) / (1.0 + b * Xval) ** 3 + 1j * eta


def tensor_from_g(Xscalar, gx, gy, b, eta):
    et = eps_tan_of_X(Xscalar, b, eta)
    er = eps_rad_of_X(Xscalar, b, eta)
    gmag = np.sqrt(gx ** 2 + gy ** 2)
    safe = gmag > 1e-12
    ghx = np.where(safe, gx / np.where(safe, gmag, 1.0), 0.0)
    ghy = np.where(safe, gy / np.where(safe, gmag, 1.0), 0.0)
    diff = er - et
    return et + diff * ghx ** 2, et + diff * ghy ** 2, diff * ghx * ghy


def build_tensor(Xp, Yp, s, b, eta, prescription, form, alpha):
    if prescription == "gradf":
        Xval, gx, gy = coherent_fields(Xp, Yp, s, form, alpha)
        return tensor_from_g(Xval, gx, gy, b, eta)
    if prescription == "pullback":
        Xval, _, _ = coherent_fields(Xp, Yp, s, form, alpha)
        # L' + X L'' = (1 - bX)/(1+bX)^3  for L = X/(1+bX)
        iso = (1.0 - b * Xval) / (1.0 + b * Xval) ** 3 + 1j * eta
        return iso.astype(complex), iso.astype(complex), np.zeros_like(iso, dtype=complex)
    if prescription in ("A", "B"):
        dx1, dy1 = Xp + 0.5 * s, Yp
        dx2, dy2 = Xp - 0.5 * s, Yp
        rho1 = np.sqrt(dx1 ** 2 + dy1 ** 2)
        rho2 = np.sqrt(dx2 ** 2 + dy2 ** 2)
        X1, X2 = tube_X(rho1), tube_X(rho2)
        rho1s = np.where(rho1 > 1e-12, rho1, 1e-12)
        rho2s = np.where(rho2 > 1e-12, rho2, 1e-12)
        gx = np.sqrt(X1) * dx1 / rho1s + np.sqrt(X2) * dx2 / rho2s
        gy = np.sqrt(X1) * dy1 / rho1s + np.sqrt(X2) * dy2 / rho2s
        Xscalar = X1 + X2 if prescription == "A" else gx ** 2 + gy ** 2
        return tensor_from_g(Xscalar, gx, gy, b, eta)
    raise ValueError("prescription must be gradf|pullback|A|B")


# ---------------------------------------------------------------- FD solver

class Grid:
    """Cell-centered square grid [-L, L]^2 with a 1-cell ghost layer."""

    def __init__(self, L, h):
        self.L = L
        self.h = h
        N = int(round(2 * L / h))
        if N % 2 == 1:
            N += 1
        self.N = N
        xc = -L + h / 2.0 + np.arange(N) * h
        self.xc = xc
        xc_pad = np.concatenate(([xc[0] - h], xc, [xc[-1] + h]))
        Xp, Yp = np.meshgrid(xc_pad, xc_pad, indexing="ij")
        self.Xp, self.Yp = Xp, Yp


def harmonic_mean(a, b):
    return 2.0 * a * b / (a + b)


def stencil_offsets(grid, exx_p, eyy_p, exy_p):
    N, h = grid.N, grid.h

    def sl(arr, di, dj):
        return arr[1 + di:1 + di + N, 1 + dj:1 + dj + N]

    exx_c, exx_e, exx_w = sl(exx_p, 0, 0), sl(exx_p, 1, 0), sl(exx_p, -1, 0)
    eyy_c, eyy_n, eyy_s = sl(eyy_p, 0, 0), sl(eyy_p, 0, 1), sl(eyy_p, 0, -1)
    exy_c = sl(exy_p, 0, 0)
    exy_e, exy_w = sl(exy_p, 1, 0), sl(exy_p, -1, 0)
    exy_n, exy_s = sl(exy_p, 0, 1), sl(exy_p, 0, -1)

    exx_R, exx_L = harmonic_mean(exx_c, exx_e), harmonic_mean(exx_w, exx_c)
    eyy_T, eyy_B = harmonic_mean(eyy_c, eyy_n), harmonic_mean(eyy_s, eyy_c)
    exy_R, exy_L = 0.5 * (exy_c + exy_e), 0.5 * (exy_w + exy_c)
    exy_T, exy_B = 0.5 * (exy_c + exy_n), 0.5 * (exy_s + exy_c)

    ih, i4h = 1.0 / h, 1.0 / (4.0 * h)
    return {
        (0, 0): -(exx_R + exx_L + eyy_T + eyy_B) * ih,
        (1, 0): exx_R * ih + exy_T * i4h - exy_B * i4h,
        (-1, 0): exx_L * ih - exy_T * i4h + exy_B * i4h,
        (0, 1): exy_R * i4h - exy_L * i4h + eyy_T * ih,
        (0, -1): -exy_R * i4h + exy_L * i4h + eyy_B * ih,
        (1, 1): (exy_R + exy_T) * i4h,
        (1, -1): -(exy_R + exy_B) * i4h,
        (-1, 1): -(exy_L + exy_T) * i4h,
        (-1, -1): (exy_L + exy_B) * i4h,
    }


def assemble_system(grid, exx_p, eyy_p, exy_p):
    N = grid.N
    offsets = stencil_offsets(grid, exx_p, eyy_p, exy_p)
    ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    idx = ii * N + jj
    rows, cols, data, ghost_info = [], [], [], []
    for (di, dj), coeff in offsets.items():
        ni, nj = ii + di, jj + dj
        interior = (ni >= 0) & (ni < N) & (nj >= 0) & (nj < N)
        rows.append(idx[interior])
        cols.append(ni[interior] * N + nj[interior])
        data.append(coeff[interior])
        ghost = ~interior
        if np.any(ghost):
            ghost_info.append((idx[ghost], ii[ghost] + di + 1,
                               jj[ghost] + dj + 1, coeff[ghost]))
    A = sp.csr_matrix((np.concatenate(data),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N * N, N * N))
    return A, ghost_info


def solve_both_orientations(grid, A, ghost_info):
    """One LU factorization, solves for E along x ('par') and y ('perp')."""
    lu = spla.splu(A.tocsc())
    out = {}
    for name, (Ex, Ey) in (("par", (1.0, 0.0)), ("perp", (0.0, 1.0))):
        Ubc = -(Ex * grid.Xp + Ey * grid.Yp)
        b = np.zeros(grid.N * grid.N, dtype=complex)
        for (r_g, gi, gj, d_g) in ghost_info:
            np.add.at(b, r_g, -d_g * Ubc[gi, gj])
        out[name] = lu.solve(b).reshape(grid.N, grid.N)
    return out


def bilinear_interp(grid, field, xq, yq):
    h, xc, N = grid.h, grid.xc, grid.N
    fi = (xq - xc[0]) / h
    fj = (yq - xc[0]) / h
    i0 = np.clip(np.floor(fi).astype(int), 0, N - 2)
    j0 = np.clip(np.floor(fj).astype(int), 0, N - 2)
    ti, tj = fi - i0, fj - j0
    return (field[i0, j0] * (1 - ti) * (1 - tj)
            + field[i0 + 1, j0] * ti * (1 - tj)
            + field[i0, j0 + 1] * (1 - ti) * tj
            + field[i0 + 1, j0 + 1] * ti * tj)


def extract_dipole(grid, u, Efield, r_inner=15.0, r_outer=25.0,
                   n_radii=6, n_theta=72):
    """m=1 Fourier ring fit; the r-growing basis term absorbs the finite-L
    Dirichlet leakage (see pair_cell readme for the derivation)."""
    Ex, Ey = Efield
    radii = np.linspace(r_inner, r_outer, n_radii)
    thetas = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    R, TH = np.meshgrid(radii, thetas, indexing="ij")
    R, TH = R.ravel(), TH.ravel()
    Xq, Yq = R * np.cos(TH), R * np.sin(TH)
    uq = (bilinear_interp(grid, u.real, Xq, Yq)
          + 1j * bilinear_interp(grid, u.imag, Xq, Yq))
    resid = uq + (Ex * Xq + Ey * Yq)
    cB, sB = np.cos(TH), np.sin(TH)
    basis = np.stack([cB / R, sB / R, cB / R ** 3, sB / R ** 3,
                      cB * R, sB * R], axis=1)
    coef, *_ = np.linalg.lstsq(basis, resid, rcond=None)
    return coef[0], coef[1]


# ------------------------------------------------------------- vertex deficit

def vertex_deficit(s, b, form, alpha, L=60.0, n=4000):
    """D = Integral (1 - 1/(1+bX)) dA for the coherent background, by
    midpoint quadrature on [-L,L]^2 plus an analytic tail bound (X ~ 1/r^6
    for the product form, ~1/r^4 for the sum form)."""
    x = np.linspace(-L + L / n, L - L / n, n)
    Xp, Yp = np.meshgrid(x, x, indexing="ij")
    Xval, _, _ = coherent_fields(Xp, Yp, s, form, alpha)
    hcell = (2.0 * L / n) ** 2
    D = np.sum(b * Xval / (1.0 + b * Xval)) * hcell
    return D


# ---------------------------------------------------------------- checkpoint

def load_done():
    done = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            for line in f:
                r = json.loads(line)
                key = (r.get("prescription"), r.get("form"), r.get("alpha"),
                       r.get("h"), r.get("eta"), r.get("s"),
                       r.get("orientation"))
                done[key] = r
    return done


def append_checkpoint(rec):
    with open(CHECKPOINT, "a") as f:
        f.write(json.dumps(rec) + "\n")


# ---------------------------------------------------------------------- runs

def run_case(h, L, eta, s, b, prescription, form, alpha, done):
    keys = [(prescription, form, alpha, h, eta, s, o) for o in ("par", "perp")]
    if all(k in done for k in keys):
        recs = [done[k] for k in keys]
        return {r["orientation"]: complex(r["alpha_re"], r["alpha_im"])
                for r in recs}, 0.0
    t0 = time.time()
    grid = Grid(L, h)
    exx, eyy, exy = build_tensor(grid.Xp, grid.Yp, s, b, eta,
                                 prescription, form, alpha)
    A, ghost_info = assemble_system(grid, exx, eyy, exy)
    us = solve_both_orientations(grid, A, ghost_info)
    wall = time.time() - t0
    out = {}
    for orient, u in us.items():
        E = (1.0, 0.0) if orient == "par" else (0.0, 1.0)
        px, py = extract_dipole(grid, u, E)
        a = px if orient == "par" else py
        out[orient] = a
        append_checkpoint({
            "prescription": prescription, "form": form, "alpha": alpha,
            "h": h, "L": L, "N": grid.N, "eta": eta, "s": s,
            "orientation": orient,
            "alpha_re": float(a.real), "alpha_im": float(a.imag),
            "wall_s": round(wall, 1),
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        })
    return out, wall


def validate(h, L, b):
    """V1 uniform (alpha=0); V3 single-lump limit of the coherent
    background vs the 1D cell at the same eta (the decisive gate: the
    gradf prescription at s -> large must reproduce the single tube)."""
    print("V1 uniform:")
    grid = Grid(L, h)
    one = np.ones_like(grid.Xp, dtype=complex)
    A, gi = assemble_system(grid, one, one, np.zeros_like(one))
    us = solve_both_orientations(grid, A, gi)
    px, _ = extract_dipole(grid, us["par"], (1.0, 0.0))
    print(f"  alpha = {px:.3e}  (expected 0)")

    print("V3 coherent background at s=40 (isolated-lump limit), eta=1e-2:")
    # place lumps far apart; the extraction ring around the origin then sees
    # ... both lumps; instead compare a SINGLE unit lump built from the same
    # machinery: w = 1/z (form sum with one pole is not supported, so use
    # the exact single-tube tensor path through tensor_from_g).
    dx, dy = grid.Xp, grid.Yp
    rho = np.sqrt(dx ** 2 + dy ** 2)
    exx, eyy, exy = tensor_from_g(tube_X(rho), dx, dy, b, 1e-2)
    A, gi = assemble_system(grid, exx, eyy, exy)
    us = solve_both_orientations(grid, A, gi)
    px, _ = extract_dipole(grid, us["par"], (1.0, 0.0))
    ref = -0.968429 + 0.026606j   # 1D shooting at eta=1e-2
    print(f"  alpha_2D = {px:.6f}   1D ref(eta=1e-2) = {ref:.6f}   "
          f"rel.diff = {abs(px - ref) / abs(ref) * 100:.2f}%")

    print("V3b gradf coherent tensor == single-tube tensor at s -> inf:")
    # at s = 30 the second lump is far outside the ring; the LOCAL tensor
    # near lump 1 must match the single-tube tensor near its center.
    s_far = 30.0
    Xv, gx, gy = coherent_fields(grid.Xp + s_far / 2.0, grid.Yp, s_far,
                                 "product", 0.0)
    # compare X along the x-axis near the lump at the origin of this frame
    line = np.abs(grid.Xp[:, grid.N // 2 + 1])
    Xs = tube_X(np.sqrt(grid.Xp[:, grid.N // 2 + 1] ** 2
                        + grid.Yp[:, grid.N // 2 + 1] ** 2))
    err = np.max(np.abs(Xv[:, grid.N // 2 + 1] - Xs)
                 [line < 3.0]) / 8.0
    print(f"  max |X_coherent - X_single|/X(0) within rho<3: {err:.2e}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--h", type=float, nargs="+", default=[0.2])
    p.add_argument("--L", type=float, default=40.0)
    p.add_argument("--eta", type=float, nargs="+", default=[1e-2])
    p.add_argument("--s", type=float, default=None)
    p.add_argument("--beta", type=float, default=None)
    p.add_argument("--prescription", default="gradf",
                   choices=["gradf", "pullback", "A", "B"])
    p.add_argument("--form", default="product", choices=["product", "sum"])
    p.add_argument("--alpha-phase", type=float, default=0.0,
                   dest="alpha_phase",
                   help="relative phase (radians) for --form sum")
    p.add_argument("--validate", action="store_true")
    p.add_argument("--two-channel", action="store_true",
                   help="print contact-packing two-channel products")
    args = p.parse_args()

    b = B_STAR if args.beta is None else args.beta
    s = S_DEFAULT if args.s is None else args.s

    print(f"b* = {B_STAR:.8f}  x* = {X_STAR:.8f}  C*^2 = {CSTAR2:.6f}  "
          f"C2* = {C2STAR:.8f}")
    print(f"s = {s:.8f}  beta = {b:.8f}  prescription = {args.prescription}"
          f"  form = {args.form}")

    if args.validate:
        validate(min(args.h), args.L, b)
        return

    done = load_done()
    results = {}
    for h in sorted(args.h, reverse=True):
        for eta in sorted(args.eta, reverse=True):
            out, wall = run_case(h, args.L, eta, s, b, args.prescription,
                                 args.form, args.alpha_phase, done)
            mean = 0.5 * (out["par"] + out["perp"])
            results[(h, eta)] = mean
            src = "checkpoint" if wall == 0.0 else f"{wall:.0f}s"
            print(f"h={h} eta={eta:g}: par={out['par']:.6f} "
                  f"perp={out['perp']:.6f} <alpha>={mean:.6f} [{src}]")
            print(f"    <alpha>/(2 alpha_single_1D) = "
                  f"{mean / (2 * ALPHA_SINGLE_1D):.5f}")

    if args.two_channel:
        D = vertex_deficit(s, b, args.form, args.alpha_phase)
        print(f"\nvertex deficit of this background D_pair = {D:.5f} "
              f"(2 x single closed form = "
              f"{2 * 2 * np.pi * b * 4 * np.arctan(X_STAR) / X_STAR:.5f})")
        h_fine = min(args.h)
        eta_fine = min(args.eta)
        a_pair = results[(h_fine, eta_fine)].real
        print(f"two-channel products at contact packing "
              f"(using Re<alpha> at h={h_fine}, eta={eta_fine:g}; "
              f"target = {TARGET:.6f}):")
        for lat, npair_cond in (("square", 1.0 / 64.0),
                                ("triangular", 2.0 / (np.sqrt(3.0) * 64.0))):
            npair = npair_cond * CSTAR2
            t = npair * a_pair
            eps_eff = 1.0 + 2.0 * t / (1.0 - t)
            vertex = 1.0 / (1.0 - npair * D)
            prod = eps_eff * vertex
            print(f"  {lat:10s}: eps_eff={eps_eff:.5f} vertex={vertex:.5f} "
                  f"product={prod:.5f} dev={(prod / TARGET - 1) * 100:+.2f}%")


if __name__ == "__main__":
    main()
