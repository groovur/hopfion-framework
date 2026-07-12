"""Phase-2 ladder: symmetry-reduced (quadrant) solves of the pair problem
on finer grids, plus matching single-tube references at every grid level.
Scope: primary deliverable is grid-converged Re(alpha) at eta = 1e-2 and
3e-3 (the etas whose resonant-annulus width is resolvable on affordable
uniform grids); Im(alpha) at smaller eta is criterion-limited and not
pursued here. Grid ladder h = 0.05, 0.025. The h = 0.0125 level (10.2M
unknowns per quadrant) is projected at 25-40 GB for the sparse LU and is
skipped under the 20 GB ceiling.

Strictly sequential; every solve is checkpointed to checkpoint.jsonl the
moment it completes, and the run is resumable (checkpointed points are
skipped)."""
import time
import numpy as np
import common2d as c2
import checkpoint as ck
from phase1 import compute_s

L_DOMAIN = 40.0
R_INNER, R_OUTER = 15.0, 25.0
BETA = c2.BETA_STAR
S_PAIR = compute_s()[0]

ETAS = [1e-2, 3e-3]
LADDER = [0.05, 0.025]

PARITY = {"par": ((-1, 1), (1.0, 0.0), "cos"),
          "perp": ((1, -1), (0.0, 1.0), "sin")}


def run_pair_quadrant(gq, h, model, eta, orientation, done, s=S_PAIR,
                       task="pair_quad"):
    if ck.have(done, task=task, h=h, model=model, eta=eta,
               orientation=orientation, s=s):
        print(f"skip: {task} h={h} {model} eta={eta:g} {orientation}",
              flush=True)
        return
    parity, E, comp = PARITY[orientation]
    t0 = time.time()
    exx, eyy, exy = c2.pair_tensor(gq.Xp, gq.Yp, BETA, eta, s, model)
    u = c2.solve_quadrant(gq, exx, eyy, exy, parity, E)
    a = c2.extract_dipole_quadrant(gq, u, E, comp,
                                    r_inner=R_INNER, r_outer=R_OUTER)
    wall = time.time() - t0
    ck.append(dict(task=task, h=h, L=L_DOMAIN, N=gq.N, model=model, eta=eta,
                    s=s, orientation=orientation,
                    alpha_re=float(a.real), alpha_im=float(a.imag),
                    wall_s=round(wall, 1)))
    print(f"{task} h={h} {model} eta={eta:g} {orientation}: alpha={a:.6f} "
          f"({wall:.0f}s)", flush=True)


def run_single_quadrant(gq, h, eta, done, task="single_quad"):
    if ck.have(done, task=task, h=h, eta=eta):
        print(f"skip: {task} h={h} eta={eta:g}", flush=True)
        return
    parity, E, comp = PARITY["par"]
    t0 = time.time()
    exx, eyy, exy = c2.single_tube_tensor(gq.Xp, gq.Yp, BETA, eta)
    u = c2.solve_quadrant(gq, exx, eyy, exy, parity, E)
    a = c2.extract_dipole_quadrant(gq, u, E, comp,
                                    r_inner=R_INNER, r_outer=R_OUTER)
    wall = time.time() - t0
    ck.append(dict(task=task, h=h, L=L_DOMAIN, N=gq.N, eta=eta,
                    alpha_re=float(a.real), alpha_im=float(a.imag),
                    wall_s=round(wall, 1)))
    print(f"{task} h={h} eta={eta:g}: alpha={a:.6f} ({wall:.0f}s)", flush=True)


def main():
    t_start = time.time()
    done = ck.load()
    print(f"phase2 start: s={S_PAIR:.8f} beta={BETA:.8f}", flush=True)

    # Gate: quadrant-vs-full cross-check points at h=0.2 and h=0.1
    # (full-domain values already checkpointed by phase 1).
    for h in [0.2, 0.1]:
        gq = c2.QuadrantGrid(L=L_DOMAIN, h=h)
        for model in ["A", "B"]:
            for orientation in ["par", "perp"]:
                run_pair_quadrant(gq, h, model, 1e-2, orientation, done)
        del gq

    # Ladder levels: pairs, then the matching single-tube reference, plus
    # the single-tube reference at the gate levels for same-grid ratios.
    for h in [0.2, 0.1]:
        gq = c2.QuadrantGrid(L=L_DOMAIN, h=h)
        for eta in ETAS:
            run_single_quadrant(gq, h, eta, done)
        del gq
    for h in LADDER:
        gq = c2.QuadrantGrid(L=L_DOMAIN, h=h)
        for model in ["A", "B"]:
            for eta in ETAS:
                for orientation in ["par", "perp"]:
                    run_pair_quadrant(gq, h, model, eta, orientation, done)
        for eta in ETAS:
            run_single_quadrant(gq, h, eta, done)
        del gq

    print(f"phase2 total wall: {time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
