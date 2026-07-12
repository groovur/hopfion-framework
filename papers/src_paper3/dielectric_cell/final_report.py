import time
import numpy as np
import common as c

ETAS = np.array([1e-2, 1e-3, 1e-4, 1e-5])

def run_beta_table(beta, rho_max=200.0, rho_min=1e-4, rtol=1e-10, atol=1e-12):
    er, et = c.bg_eps_fns(beta)
    rows = []
    for eta in ETAS:
        alpha, info = c.shoot_alpha(er, et, eta=eta, rho_min=rho_min, rho_max=rho_max,
                                     rtol=rtol, atol=atol)
        rows.append((eta, alpha.real, alpha.imag, info['nfev']))
    return rows

def extrapolate_eta0(rows, npts=2):
    """Linear least-squares fit alpha(eta) = alpha0 + k*eta using the
    `npts` smallest-eta points; returns (Re0, Im0, slope_re, slope_im)."""
    sub = sorted(rows, key=lambda r: r[0])[:npts]
    etas = np.array([r[0] for r in sub])
    res = np.array([r[1] for r in sub])
    ims = np.array([r[2] for r in sub])
    A = np.vstack([np.ones_like(etas), etas]).T
    (re0, kre), *_ = np.linalg.lstsq(A, res, rcond=None)
    (im0, kim), *_ = np.linalg.lstsq(A, ims, rcond=None)
    return re0, im0, kre, kim

def grid_doubling_check(beta, rho_max=200.0, eta=1e-4):
    er, et = c.bg_eps_fns(beta)
    out = []
    for rtol in [1e-9, 1e-10, 1e-11, 1e-12]:
        alpha, info = c.shoot_alpha(er, et, eta=eta, rho_min=1e-4, rho_max=rho_max,
                                     rtol=rtol, atol=rtol*1e-2)
        out.append((rtol, alpha, info['nfev']))
    return out

def rho_max_check(beta, eta=1e-4):
    er, et = c.bg_eps_fns(beta)
    out = []
    for rmax in [100.0, 200.0, 400.0]:
        alpha, info = c.shoot_alpha(er, et, eta=eta, rho_min=1e-4, rho_max=rmax)
        out.append((rmax, alpha))
    return out

def rho_min_check(beta, eta=1e-4, rho_max=200.0):
    er, et = c.bg_eps_fns(beta)
    out = []
    for rmin in [1e-3, 1e-4, 1e-5]:
        alpha, info = c.shoot_alpha(er, et, eta=eta, rho_min=rmin, rho_max=rho_max)
        out.append((rmin, alpha))
    return out


def main():
    lines = []
    def p(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        lines.append(s)

    t_start = time.time()

    p("x_star  =", format(c.X_STAR, ".8f"))
    p("beta_star =", format(c.BETA_STAR, ".8f"))
    p("rho_c(beta_star) =", format(c.rho_c_of_beta(c.BETA_STAR), ".8f"))
    p("rho_c(beta_star/2) =", c.rho_c_of_beta(c.BETA_STAR/2.0))
    p("rho_c(2*beta_star) =", format(c.rho_c_of_beta(2.0*c.BETA_STAR), ".8f"))
    p()

    # ---- validations ----
    p("=== V1: uniform medium (beta=0) ===")
    er, et = c.bg_eps_fns(0.0)
    a, info = c.shoot_alpha(er, et, eta=0.0)
    p(f"alpha = {a}  (expected 0)")
    p()

    p("=== V2: step inclusion eps1=3, R=1 ===")
    er, et = c.step_eps_fns(3.0, 1.0)
    a, info = c.shoot_alpha(er, et, eta=0.0)
    p(f"alpha = {a}  expected = 0.5  diff = {abs(a-0.5):.3e}")
    p()

    p("=== V3: step inclusion eps1=-3+i*1e-3, R=1 ===")
    eps1 = -3.0 + 1j*1e-3
    er, et = c.step_eps_fns(eps1, 1.0)
    a, info = c.shoot_alpha(er, et, eta=0.0)
    expected = (eps1-1)/(eps1+1)
    p(f"alpha = {a}  expected = {expected}  diff = {abs(a-expected):.3e}")
    p()

    # ---- physics runs ----
    betas = {
        "R1 beta=beta_star": c.BETA_STAR,
        "R2a beta=beta_star/2": c.BETA_STAR/2.0,
        "R2b beta=2*beta_star": 2.0*c.BETA_STAR,
    }
    extraps = {}
    for label, beta in betas.items():
        p(f"=== {label} (beta={beta:.8f}, rho_c={c.rho_c_of_beta(beta)}) ===")
        rows = run_beta_table(beta)
        for eta, re, im, nfev in rows:
            p(f"  eta={eta:.0e}: Re(alpha)={re:.8f}  Im(alpha)={im:.8e}  nfev={nfev}")
        re0, im0, kre, kim = extrapolate_eta0(rows, npts=2)
        p(f"  eta->0 (linear extrap. from eta=1e-4,1e-5): Re0={re0:.8f}  Im0={im0:.8e}")
        extraps[label] = (re0, im0)
        # grid/resolution doubling check
        gd = grid_doubling_check(beta)
        p("  rtol refinement (grid-doubling analog):")
        for rtol, alpha, nfev in gd:
            p(f"    rtol={rtol:.0e}: alpha={alpha}  nfev={nfev}")
        rmaxc = rho_max_check(beta)
        p("  rho_max sensitivity:")
        for rmax, alpha in rmaxc:
            p(f"    rho_max={rmax}: alpha={alpha}")
        rminc = rho_min_check(beta)
        p("  rho_min sensitivity:")
        for rmin, alpha in rminc:
            p(f"    rho_min={rmin}: alpha={alpha}")
        p()

    p(f"Total compute time: {time.time()-t_start:.3f}s")
    p()

    # ---- R3 dilute effective medium ----
    p("=== R3: 2D Maxwell-Garnett / Clausius-Mossotti ===")
    p("Formula used: eps_eff = 1 + 2*n*alpha_pol/(1 - n*alpha_pol)")
    p("Inverted:      n = (eps_eff-1) / ( (eps_eff+1) * alpha_pol )")
    re0, im0 = extraps["R1 beta=beta_star"]
    alpha_c = complex(re0, im0)
    alpha_mod_signed = complex(np.sign(re0)*abs(alpha_c), 0.0)  # modulus with sign of Re
    p(f"Using R1 (beta=beta_star) eta->0 extrapolated alpha = {alpha_c}")
    p(f"  |alpha| = {abs(alpha_c):.8f}   Re(alpha) = {re0:.8f}")
    targets = {"eps_eff=1.09326": 1.09326, "eps_eff=1/1.09326=0.914695": 1.0/1.09326}
    for tlabel, eps_eff in targets.items():
        n_re = (eps_eff - 1.0) / ((eps_eff + 1.0) * re0)
        n_mod = (eps_eff - 1.0) / ((eps_eff + 1.0) * alpha_mod_signed.real)
        n_full_complex = (eps_eff - 1.0) / ((eps_eff + 1.0) * alpha_c)
        p(f"  {tlabel}: n(using Re alpha)={n_re:.8f}   n(using |alpha|, signed)={n_mod:.8f}   n(full complex)={n_full_complex}")
    p()

    with open("results.md", "w") as f:
        f.write("# Dielectric-cell dipole polarizability: results\n\n```\n")
        f.write("\n".join(lines))
        f.write("\n```\n")

if __name__ == "__main__":
    main()
