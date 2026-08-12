"""
Metric (field-overlap) interaction energy between two Q_H=2 hopfion tubes,
as a function of center separation d, in the density-feedback framework.

Each hopfion is idealized as a circular tube of radius R0 with kinetic-density
kernel kappa_i(x) = 8/(1+rho_i(x)^2)^2, rho_i(x) = Euclidean distance from x
to circle i (the curve, not its center). Two geometric configurations are
computed:

  perp: circle 1 in the xy-plane centered at the origin; circle 2 in the
        xz-plane centered at (d,0,0). This is the linking-relevant geometry.
  coax: circle 1 in the xy-plane centered at the origin; circle 2 parallel
        to circle 1, coaxial, centered at (0,0,d). Reference geometry.

Only the superposed-kernel interaction is computed here (kappa = kappa_1 +
kappa_2, no field relaxation).

Functionals of a kernel field kappa on R^3:
  J2a[kappa]        = Integral kappa dV
  Jfb[kappa;beta]   = Integral kappa/(1+beta*kappa) dV
  J4[kappa]         = Integral kappa^2/4 dV      (quartic proxy)
  Kfb[kappa;beta]   = J2a[kappa] + mustar*Jfb[kappa;beta]
  E[kappa;beta]     = Kfb[kappa;beta] * J4[kappa]

Interaction quantities reported, for F in {J2a, Jfb, J4}:
  Delta F(d) = F[kappa_1+kappa_2] - F[kappa_1] - F[kappa_2]

The product form E is not additive even for non-overlapping fields, so the
physically meaningful quantity there is
  E_pair(d) = Kfb[kappa_1+kappa_2] * J4[kappa_1+kappa_2]
i.e. E evaluated directly on the superposed field, as a function of d.

Numerics: midpoint-rule quadrature on a 3D box, looping over slabs along the
axis of separation so memory stays bounded (each slab is a 2D array in the
other two coordinates). kappa_1, kappa_2, and their sum are evaluated at
IDENTICAL grid points within a given (d, config, h) run, so quadrature bias
cancels in the differences Delta F.
"""

import json
import time
import numpy as np

R0 = 3.0
PAD = 12.0
MUSTAR = 3.0 - (1.0 + np.sqrt(5.0)) / 2.0  # = 3 - phi
BETAS = (0.452, 0.06715131)
D_GRID = [3.5, 4.0, 4.5, 5.0, 5.5, 5.8, 6.0, 6.2, 6.5, 7.0, 8.0, 10.0]
H_LADDER = [0.25, 0.125, 0.0625]
CONFIGS = ("perp", "coax")

assert abs(MUSTAR - 1.3819660113) < 1e-9, MUSTAR


def box_bounds(d, config, pad=PAD):
    """Return (xlo,xhi, ylo,yhi, zlo,zhi) for the box containing both circles."""
    if config == "perp":
        return (-pad, d + pad, -pad, pad, -pad, pad)
    elif config == "coax":
        return (-pad, pad, -pad, pad, -pad, d + pad)
    else:
        raise ValueError(config)


def rho1(x, y, z, r0=R0):
    """Distance from (x,y,z) to circle 1: radius r0, in xy-plane, centered at origin."""
    s = np.sqrt(x ** 2 + y ** 2)
    return np.sqrt((s - r0) ** 2 + z ** 2)


def rho2(x, y, z, d, config, r0=R0):
    """Distance from (x,y,z) to circle 2, depending on configuration."""
    if config == "perp":
        # circle 2: radius r0, in xz-plane, centered at (d,0,0)
        s2 = np.sqrt((x - d) ** 2 + z ** 2)
        return np.sqrt((s2 - r0) ** 2 + y ** 2)
    elif config == "coax":
        # circle 2: radius r0, in xy-plane, centered at (0,0,d)
        s2 = np.sqrt(x ** 2 + y ** 2)
        return np.sqrt((s2 - r0) ** 2 + (z - d) ** 2)
    else:
        raise ValueError(config)


def compute(d, h, config, betas=BETAS, mustar=MUSTAR, r0=R0, pad=PAD):
    """Midpoint-rule quadrature of all needed functionals at grid spacing h.

    Returns a dict of raw integrals: J2a_1,2,s ; J4_1,2,s ; Jfb_1,2,s[beta].
    Loops over the x-axis (slab = 2D array in y,z) for perp geometry — since
    the separation is along x there — and over z for coax (separation along
    z), so each slab is bounded in size regardless of d.
    """
    xlo, xhi, ylo, yhi, zlo, zhi = box_bounds(d, config, pad)
    Nx = int(round((xhi - xlo) / h))
    Ny = int(round((yhi - ylo) / h))
    Nz = int(round((zhi - zlo) / h))
    xs = xlo + (np.arange(Nx) + 0.5) * h
    ys = ylo + (np.arange(Ny) + 0.5) * h
    zs = zlo + (np.arange(Nz) + 0.5) * h
    dV = h ** 3

    acc = {
        "J2a_1": 0.0, "J2a_2": 0.0, "J2a_s": 0.0,
        "J4_1": 0.0, "J4_2": 0.0, "J4_s": 0.0,
    }
    for b in betas:
        acc[f"Jfb_1_{b}"] = 0.0
        acc[f"Jfb_2_{b}"] = 0.0
        acc[f"Jfb_s_{b}"] = 0.0

    # loop over the axis along which the box is elongated by d (x for perp,
    # z for coax); the other two axes form the 2D slab via broadcasting.
    if config == "perp":
        Y = ys.reshape(-1, 1)
        Z = zs.reshape(1, -1)
        for x in xs:
            r1 = rho1(x, Y, Z, r0)
            r2 = rho2(x, Y, Z, d, config, r0)
            k1 = 8.0 / (1.0 + r1 ** 2) ** 2
            k2 = 8.0 / (1.0 + r2 ** 2) ** 2
            ks = k1 + k2
            acc["J2a_1"] += k1.sum()
            acc["J2a_2"] += k2.sum()
            acc["J2a_s"] += ks.sum()
            acc["J4_1"] += (k1 ** 2).sum()
            acc["J4_2"] += (k2 ** 2).sum()
            acc["J4_s"] += (ks ** 2).sum()
            for b in betas:
                acc[f"Jfb_1_{b}"] += (k1 / (1.0 + b * k1)).sum()
                acc[f"Jfb_2_{b}"] += (k2 / (1.0 + b * k2)).sum()
                acc[f"Jfb_s_{b}"] += (ks / (1.0 + b * ks)).sum()
    else:  # coax
        X = xs.reshape(-1, 1)
        Y = ys.reshape(1, -1)
        for z in zs:
            r1 = rho1(X, Y, z, r0)
            r2 = rho2(X, Y, z, d, config, r0)
            k1 = 8.0 / (1.0 + r1 ** 2) ** 2
            k2 = 8.0 / (1.0 + r2 ** 2) ** 2
            ks = k1 + k2
            acc["J2a_1"] += k1.sum()
            acc["J2a_2"] += k2.sum()
            acc["J2a_s"] += ks.sum()
            acc["J4_1"] += (k1 ** 2).sum()
            acc["J4_2"] += (k2 ** 2).sum()
            acc["J4_s"] += (ks ** 2).sum()
            for b in betas:
                acc[f"Jfb_1_{b}"] += (k1 / (1.0 + b * k1)).sum()
                acc[f"Jfb_2_{b}"] += (k2 / (1.0 + b * k2)).sum()
                acc[f"Jfb_s_{b}"] += (ks / (1.0 + b * ks)).sum()

    out = {}
    out["J2a_1"] = acc["J2a_1"] * dV
    out["J2a_2"] = acc["J2a_2"] * dV
    out["J2a_s"] = acc["J2a_s"] * dV
    out["J4_1"] = 0.25 * acc["J4_1"] * dV
    out["J4_2"] = 0.25 * acc["J4_2"] * dV
    out["J4_s"] = 0.25 * acc["J4_s"] * dV
    for b in betas:
        out[f"Jfb_1_{b}"] = acc[f"Jfb_1_{b}"] * dV
        out[f"Jfb_2_{b}"] = acc[f"Jfb_2_{b}"] * dV
        out[f"Jfb_s_{b}"] = acc[f"Jfb_s_{b}"] * dV

    out["Nx"], out["Ny"], out["Nz"] = Nx, Ny, Nz
    out["npoints"] = Nx * Ny * Nz
    return out


def derived(out, betas=BETAS, mustar=MUSTAR):
    """Compute Delta J2a, Delta Jfb(beta), Delta J4, Kfb, E_pair, etc from raw integrals."""
    d = {}
    d["Delta_J2a"] = out["J2a_s"] - out["J2a_1"] - out["J2a_2"]
    d["Delta_J4"] = out["J4_s"] - out["J4_1"] - out["J4_2"]
    for b in betas:
        d[f"Delta_Jfb_{b}"] = out[f"Jfb_s_{b}"] - out[f"Jfb_1_{b}"] - out[f"Jfb_2_{b}"]
        Kfb_1 = out["J2a_1"] + mustar * out[f"Jfb_1_{b}"]
        Kfb_2 = out["J2a_2"] + mustar * out[f"Jfb_2_{b}"]
        Kfb_s = out["J2a_s"] + mustar * out[f"Jfb_s_{b}"]
        E_1 = Kfb_1 * out["J4_1"]
        E_2 = Kfb_2 * out["J4_2"]
        E_s = Kfb_s * out["J4_s"]
        d[f"Kfb_1_{b}"] = Kfb_1
        d[f"Kfb_2_{b}"] = Kfb_2
        d[f"Kfb_s_{b}"] = Kfb_s
        d[f"E_1_{b}"] = E_1
        d[f"E_2_{b}"] = E_2
        d[f"E_pair_{b}"] = E_s
        d[f"E_int_naive_{b}"] = E_s - E_1 - E_2
    return d


def thin_tube_J2a(r0=R0):
    """Analytic thin-tube estimate: 2*pi*R0 * Integral_0^inf 8/(1+r^2)^2 * 2*pi*r dr.
    Integral_0^inf r/(1+r^2)^2 dr = 1/2, so inner integral = 8*2*pi*(1/2) = 8*pi.
    Thin-tube J2a = 2*pi*R0*8*pi = 16*pi^2*R0.
    """
    return 16.0 * np.pi ** 2 * r0


def main():
    t_start = time.time()
    results = {}

    # --- coarse timing pass at h=0.25 over full (d, config) grid ---
    t0 = time.time()
    for config in CONFIGS:
        for d in D_GRID:
            key = (config, d, 0.25)
            results[key] = compute(d, 0.25, config)
    t_coarse = time.time() - t0
    n_runs_coarse = len(CONFIGS) * len(D_GRID)
    print(f"[coarse h=0.25] {n_runs_coarse} runs in {t_coarse:.2f} s "
          f"({t_coarse / n_runs_coarse:.3f} s/run)")

    # project cost of remaining ladder levels (h=0.125 ~8x points, h=0.0625 ~64x points
    # relative to h=0.25, since Nx*Ny*Nz scales as (1/h)^3)
    proj_0125 = t_coarse * 8
    proj_00625 = t_coarse * 64
    total_proj = t_coarse + proj_0125 + proj_00625
    print(f"[projection] h=0.125 full ladder ~ {proj_0125:.1f} s, "
          f"h=0.0625 full ladder ~ {proj_00625:.1f} s, "
          f"total projected ~ {total_proj:.1f} s ({total_proj/60:.1f} min)")

    run_h0125 = True
    run_h00625 = total_proj < 30 * 60

    if run_h0125:
        t0 = time.time()
        for config in CONFIGS:
            for d in D_GRID:
                key = (config, d, 0.125)
                results[key] = compute(d, 0.125, config)
        t_mid = time.time() - t0
        print(f"[h=0.125] {n_runs_coarse} runs in {t_mid:.2f} s")

    if run_h00625:
        t0 = time.time()
        for config in CONFIGS:
            for d in D_GRID:
                key = (config, d, 0.0625)
                results[key] = compute(d, 0.0625, config)
        t_fine = time.time() - t0
        print(f"[h=0.0625] {n_runs_coarse} runs in {t_fine:.2f} s")
    else:
        print("[h=0.0625] skipped: projected total exceeds 30 minutes")

    t_total = time.time() - t_start
    print(f"[total wall time] {t_total:.1f} s ({t_total/60:.2f} min)")

    # serialize raw + derived results
    serializable = {}
    for key, out in results.items():
        config, d, h = key
        skey = f"{config}|{d}|{h}"
        rec = dict(out)
        rec.update(derived(out))
        serializable[skey] = rec

    with open("raw_results.json", "w") as f:
        json.dump(serializable, f, indent=2)

    print("wrote raw_results.json")
    print(f"thin-tube analytic J2a (single circle, R0={R0}): {thin_tube_J2a():.6f}")


if __name__ == "__main__":
    main()
