"""Phase-1 runs: pair eta sweep at h=0.2, pair h=0.1 at eta=1e-2,
V4 far-separation check, single-tube h=0.05 V3 gate at eta=1e-2.
Strictly sequential; every solve is checkpointed to checkpoint.jsonl
the moment it completes."""
import time
import numpy as np
import common2d as c2
import checkpoint as ck
from scipy.optimize import brentq

L_DOMAIN = 40.0
R_INNER, R_OUTER = 15.0, 25.0
SINGLE_1D = complex(-0.97106613, 0.01353232)


def compute_s():
    PHI = c2.PHI

    def f(C):
        return C * (2 * C + 1) / ((C ** 2 + 1) * (3 * C - 1)) - 2 ** (4.0 / 3.0) / PHI ** 5

    C2star = brentq(f, 0.5, 10.0, xtol=1e-14, rtol=1e-14)
    Cstar = np.sqrt(0.75 * PHI ** 5 / 2 ** (4.0 / 3.0))
    return 2.0 * 3.0 / (C2star * Cstar), C2star, Cstar


S_PAIR, C2_STAR, C_STAR = compute_s()
BETA = c2.BETA_STAR


def run_pair_point(g, h, model, eta, s, done, task="pair"):
    """One factorization, both orientations; checkpoint each orientation."""
    if (ck.have(done, task=task, h=h, model=model, eta=eta, orientation="par", s=s)
            and ck.have(done, task=task, h=h, model=model, eta=eta,
                        orientation="perp", s=s)):
        print(f"skip (checkpointed): {task} h={h} model={model} eta={eta:g}",
              flush=True)
        return
    t0 = time.time()
    exx, eyy, exy = c2.pair_tensor(g.Xp, g.Yp, BETA, eta, s, model)
    A, gi = c2.assemble_system(g, exx, eyy, exy)
    u_par, u_perp = c2.solve_system_multi(g, A, gi, [(1.0, 0.0), (0.0, 1.0)])
    wall = time.time() - t0
    for orientation, u, E in [("par", u_par, (1.0, 0.0)),
                               ("perp", u_perp, (0.0, 1.0))]:
        px, py = c2.extract_dipole(g, u, E, r_inner=R_INNER, r_outer=R_OUTER)
        a = px if orientation == "par" else py
        ck.append(dict(task=task, h=h, L=L_DOMAIN, N=g.N, model=model, eta=eta,
                        s=s, orientation=orientation,
                        alpha_re=float(a.real), alpha_im=float(a.imag),
                        cross_re=float((py if orientation == "par" else px).real),
                        cross_im=float((py if orientation == "par" else px).imag),
                        wall_s=round(wall, 1)))
        print(f"{task} h={h} model={model} eta={eta:g} {orientation}: "
              f"alpha={a:.6f}  (wall {wall:.0f}s for pair of orientations)",
              flush=True)


def run_single_point(g, h, eta, done):
    if ck.have(done, task="single", h=h, eta=eta):
        print(f"skip (checkpointed): single h={h} eta={eta:g}", flush=True)
        return
    t0 = time.time()
    exx, eyy, exy = c2.single_tube_tensor(g.Xp, g.Yp, BETA, eta)
    A, gi = c2.assemble_system(g, exx, eyy, exy)
    u = c2.solve_system(g, A, gi, (1.0, 0.0))
    px, py = c2.extract_dipole(g, u, (1.0, 0.0), r_inner=R_INNER, r_outer=R_OUTER)
    wall = time.time() - t0
    ck.append(dict(task="single", h=h, L=L_DOMAIN, N=g.N, eta=eta,
                    alpha_re=float(px.real), alpha_im=float(px.imag),
                    wall_s=round(wall, 1)))
    print(f"single h={h} eta={eta:g}: alpha={px:.6f}  (wall {wall:.0f}s)",
          flush=True)


def main():
    t_start = time.time()
    done = ck.load()
    print(f"s = {S_PAIR:.8f}   beta = {BETA:.8f}", flush=True)

    # (a) pair h=0.2 full eta sweep, both models
    g02 = c2.Grid(L=L_DOMAIN, h=0.2)
    for model in ["A", "B"]:
        for eta in [1e-2, 3e-3, 1e-3]:
            run_pair_point(g02, 0.2, model, eta, S_PAIR, done)

    # (c) V4: s=10 model A h=0.2 eta=1e-2, plus single-tube reference at
    # the same grid/eta
    run_pair_point(g02, 0.2, "A", 1e-2, 10.0, done, task="pair_s10")
    run_single_point(g02, 0.2, 1e-2, done)
    del g02

    # (b) pair h=0.1 at eta=1e-2, both models
    g01 = c2.Grid(L=L_DOMAIN, h=0.1)
    for model in ["A", "B"]:
        run_pair_point(g01, 0.1, model, 1e-2, S_PAIR, done)
    run_single_point(g01, 0.1, 1e-2, done)
    del g01

    # (d) single tube h=0.05 at eta=1e-2 (V3 grid ladder point)
    g005 = c2.Grid(L=L_DOMAIN, h=0.05)
    run_single_point(g005, 0.05, 1e-2, done)
    del g005

    print(f"Phase 1 total wall time: {time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
