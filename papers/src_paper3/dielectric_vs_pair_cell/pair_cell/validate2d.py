"""Validation gates V1 (uniform medium) and V2 (isotropic step cylinder).
Results are checkpointed to checkpoint.jsonl; V3 (single tube vs the 1D
shooting code) and V4 (far-separation pair) are produced by phase1.py /
phase2.py and assessed in make_results.py from the same checkpoint file."""
import time
import numpy as np
import common2d as c2
import checkpoint as ck

R_INNER, R_OUTER = 15.0, 25.0
L_DOMAIN = 40.0


def uniform_tensor(g):
    exx = np.ones((g.Np, g.Np), dtype=complex)
    eyy = np.ones((g.Np, g.Np), dtype=complex)
    exy = np.zeros((g.Np, g.Np), dtype=complex)
    return exx, eyy, exy


def step_tensor(g, eps1, R):
    rho = np.sqrt(g.Xp ** 2 + g.Yp ** 2)
    val = np.where(rho < R, eps1, 1.0 + 0j)
    return c2.tensor_from_scalar(val)


def main():
    done = ck.load()

    # V1: uniform medium -> alpha = 0
    if not ck.have(done, task="v1", h=0.2):
        t0 = time.time()
        g = c2.Grid(L=L_DOMAIN, h=0.2)
        exx, eyy, exy = uniform_tensor(g)
        u = c2.assemble_and_solve(g, exx, eyy, exy, (1.0, 0.0))
        px, py = c2.extract_dipole(g, u, (1.0, 0.0), R_INNER, R_OUTER)
        ck.append(dict(task="v1", h=0.2, L=L_DOMAIN, N=g.N,
                        alpha_re=float(px.real), alpha_im=float(px.imag),
                        wall_s=round(time.time() - t0, 1)))
        print(f"V1 h=0.2: alpha={px:.3e}", flush=True)
        del g

    # V2: step cylinder eps1=3, R=1 -> alpha = +0.5, grid ladder
    for h in [0.4, 0.2, 0.1]:
        if ck.have(done, task="v2", h=h):
            continue
        t0 = time.time()
        g = c2.Grid(L=L_DOMAIN, h=h)
        exx, eyy, exy = step_tensor(g, 3.0 + 0j, 1.0)
        u = c2.assemble_and_solve(g, exx, eyy, exy, (1.0, 0.0))
        px, py = c2.extract_dipole(g, u, (1.0, 0.0), R_INNER, R_OUTER)
        ck.append(dict(task="v2", h=h, L=L_DOMAIN, N=g.N,
                        alpha_re=float(px.real), alpha_im=float(px.imag),
                        wall_s=round(time.time() - t0, 1)))
        print(f"V2 h={h}: alpha={px:.6f} rel.err={abs(px-0.5)/0.5:.3%}",
              flush=True)
        del g


if __name__ == "__main__":
    main()
