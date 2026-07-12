import time
import numpy as np
import common as c

ETAS = [1e-2, 1e-3, 1e-4, 1e-5]

def run_beta(label, beta, rho_max=200.0, rho_min=1e-4, rtol=1e-10, atol=1e-12):
    er, et = c.bg_eps_fns(beta)
    print(f"\n=== {label}  (beta={beta:.8f}, rho_c={c.rho_c_of_beta(beta)}) ===")
    results = []
    for eta in ETAS:
        t0 = time.time()
        alpha, info = c.shoot_alpha(er, et, eta=eta, rho_min=rho_min, rho_max=rho_max,
                                     rtol=rtol, atol=atol)
        dt = time.time()-t0
        if alpha is None:
            print(f"  eta={eta:.0e}: FAILED -- {info.get('message')}")
            results.append((eta, None, None))
            continue
        print(f"  eta={eta:.0e}: Re(alpha)={alpha.real:.8f}  Im(alpha)={alpha.imag:.8e}  "
              f"nfev={info['nfev']}  t={dt:.3f}s")
        results.append((eta, alpha.real, alpha.imag))
    return results

if __name__ == "__main__":
    t_start = time.time()
    res_star = run_beta("R1: beta = beta_star", c.BETA_STAR)
    res_half = run_beta("R2a: beta = beta_star/2", c.BETA_STAR/2.0)
    res_double = run_beta("R2b: beta = 2*beta_star", 2.0*c.BETA_STAR)
    print(f"\nTotal physics-run time: {time.time()-t_start:.3f}s")

    # resolution check: tighten rtol/atol further at eta=1e-4 for beta_star to confirm convergence
    print("\n--- Grid/tolerance refinement check at eta=1e-4, beta=beta_star ---")
    er, et = c.bg_eps_fns(c.BETA_STAR)
    for rt in [1e-10, 1e-11, 1e-12]:
        alpha, info = c.shoot_alpha(er, et, eta=1e-4, rho_min=1e-4, rho_max=200.0,
                                     rtol=rt, atol=rt*1e-2)
        print(f"  rtol={rt}: alpha={alpha}  nfev={info['nfev']}")

    # rho_max sensitivity at beta_star, eta=1e-4
    print("\n--- rho_max sensitivity at eta=1e-4, beta=beta_star ---")
    for rmax in [100.0, 200.0, 400.0]:
        alpha, info = c.shoot_alpha(er, et, eta=1e-4, rho_min=1e-4, rho_max=rmax)
        print(f"  rho_max={rmax}: alpha={alpha}")

    # rho_min sensitivity
    print("\n--- rho_min sensitivity at eta=1e-4, beta=beta_star ---")
    for rmin in [1e-3, 1e-4, 1e-5]:
        alpha, info = c.shoot_alpha(er, et, eta=1e-4, rho_min=rmin, rho_max=200.0)
        print(f"  rho_min={rmin}: alpha={alpha}")
