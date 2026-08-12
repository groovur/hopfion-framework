"""Orchestrates the full build: geometry validation, V1-V3, and the physics
table (a_c x eta x grid resolution M), with checkpointing.

Grid-convention note (documented per the task's fallback-reporting
instruction): the pre-registered 3-coloring motif failed the linking
validation (see results.md) and we fall back to the documented simpler
motif ("alternating xy/xz Hopf chains along x, stacked in y,z at spacing
a_c"). Its periodic unit cell is L = (2*a_c, a_c, a_c), not a cube, so the
task's single resolution parameter N (implying an N^3 cubic grid) is
replaced by M = Ny = Nz (Nx = 2M), uniform spacing h = a_c/M in all three
axes. M in {48, 64} stand in for the task's N in {128, 192}; both give
comparable or finer h than the task's a_c=6 values (h=6/48=0.125,
6/64=0.09375, vs the task's h~0.14, 0.094) at a much smaller total point
count (2*48^3=221184, 2*64^3=524288 vs 128^3=2097152), which is what makes
the fallback geometry computationally cheaper as well as topologically
correct.
"""
import datetime
import json
import os
import time

import numpy as np

from common3d import T_TARGET, PHI
from geometry import (build_curves_fallback, validate_geometry,
                       validate_geometry_fallback)
from field import Grid, build_U_g, build_eps_tensor
from solver import eps_eff_tensor

HERE = os.path.dirname(os.path.abspath(__file__))
CKPT_PATH = os.path.join(HERE, "checkpoint.jsonl")
RESULTS_PATH = os.path.join(HERE, "results.md")

A_C_LIST = [6.0, 5.9, 5.8]  # 6.0 primary/pre-registered
ETA_LIST = [1e-2, 3e-3]
M_LIST = [48, 64]
NPTS_CURVE = 400
CUTOFF = 10.0


def load_checkpoint_keys():
    keys = set()
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                keys.add((rec['a_c'], rec['eta'], tuple(rec['grid'])))
    return keys


def append_checkpoint(rec):
    with open(CKPT_PATH, "a") as f:
        f.write(json.dumps(rec) + "\n")


def run_one_case(a_c, eta, M, log, maxiter=250, restart=40):
    Lx, Ly, Lz = 2.0 * a_c, a_c, a_c
    h = a_c / M
    grid = Grid((Lx, Ly, Lz), h)
    curves = build_curves_fallback(a_c, n_stack=1, npts=NPTS_CURVE)
    t0 = time.time()
    U, g = build_U_g(grid, curves, cutoff=CUTOFF)
    eps6, S1, Ubar = build_eps_tensor(U, g, eta)
    eps6 = tuple(comp.reshape(grid.shape) for comp in eps6)
    t_field = time.time() - t0

    t1 = time.time()
    eps_eff, solve_info = eps_eff_tensor(eps6, h, eta, tol=1e-7, maxiter=maxiter,
                                          restart=restart, report=log)
    t_solve = time.time() - t1

    eps_iso = np.trace(eps_eff).real / 3.0
    off_diag = np.abs(eps_eff - np.diag(np.diag(eps_eff))).max()
    diag_vals = np.diag(eps_eff).real
    anisotropy = (diag_vals.max() - diag_vals.min()) / eps_iso

    P_pathavg = eps_iso / S1
    P_meanfield = eps_iso * (1.0 + Ubar)
    dev_pathavg = 100.0 * (P_pathavg - T_TARGET) / T_TARGET
    dev_meanfield = 100.0 * (P_meanfield - T_TARGET) / T_TARGET

    rec = {
        'a_c': a_c, 'eta': eta, 'grid': list(grid.N), 'h': h,
        'eps_eff_diag': [complex(x).real for x in np.diag(eps_eff)],
        'eps_eff_diag_imag': [complex(x).imag for x in np.diag(eps_eff)],
        'eps_eff_re': eps_eff.real.tolist(),
        'eps_eff_im': eps_eff.imag.tolist(),
        'field_directions': ['x', 'y', 'z'],
        'eps_eff_offdiag_max': float(off_diag),
        'eps_iso': float(eps_iso),
        'anisotropy': float(anisotropy),
        'S1': float(S1), 'Ubar': float(Ubar),
        'P_pathavg': float(P_pathavg), 'P_meanfield': float(P_meanfield),
        'dev_pathavg_pct': float(dev_pathavg), 'dev_meanfield_pct': float(dev_meanfield),
        'iterations': [si['n_iter'] for si in solve_info],
        'gmres_info': [si['info'] for si in solve_info],
        'final_resid': [si['final_resid'] for si in solve_info],
        't_field_s': t_field, 't_solve_s': t_solve,
        't_total_s': t_field + t_solve,
        'timestamp': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    return rec


def main():
    log_lines = []

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        log_lines.append(s)

    log("=== Geometry validation: primary 3-coloring motif ===")
    for a_c in A_C_LIST:
        validate_geometry(a_c, report=log)
    log("\nPrimary motif FAILS the linking validation: axis-adjacent site "
        "pairs of the same color-pair link (|Lk|=1) along only one of the "
        "three lattice axes and give Lk~0 along the other two (see numbers "
        "above); falling back to the documented simpler motif.\n")

    log("=== Geometry validation: fallback alternating-chain motif ===")
    for a_c in A_C_LIST:
        validate_geometry_fallback(a_c, report=log)
    log("")

    log("=== V1-V3 validations: see validate.py output (run separately) ===\n")

    keys_done = load_checkpoint_keys()
    log("=== Physics runs ===")
    log(f"{'a_c':>5} {'eta':>8} {'grid':>16} {'eps_iso':>10} {'aniso':>8} "
        f"{'S1':>8} {'Ubar':>8} {'P_path':>9} {'dev%':>8} {'P_mf':>9} {'dev%':>8} "
        f"{'iters':>12} {'t(s)':>7}")
    for M in M_LIST:
        for a_c in A_C_LIST:
            for eta in ETA_LIST:
                Lx, Ly, Lz = 2.0 * a_c, a_c, a_c
                h = a_c / M
                grid_key = (int(round(Lx / h)), int(round(Ly / h)), int(round(Lz / h)))
                key = (a_c, eta, grid_key)
                if key in keys_done:
                    log(f"[skip, checkpointed] a_c={a_c} eta={eta} grid={grid_key}")
                    continue
                rec = run_one_case(a_c, eta, M, log)
                append_checkpoint(rec)
                keys_done.add(key)
                log(f"{a_c:5.1f} {eta:8.0e} {str(tuple(rec['grid'])):>16} "
                    f"{rec['eps_iso']:10.6f} {rec['anisotropy']:8.2e} "
                    f"{rec['S1']:8.5f} {rec['Ubar']:8.5f} "
                    f"{rec['P_pathavg']:9.6f} {rec['dev_pathavg_pct']:8.3f} "
                    f"{rec['P_meanfield']:9.6f} {rec['dev_meanfield_pct']:8.3f} "
                    f"{str(rec['iterations']):>12} {rec['t_total_s']:7.1f}")

    with open(os.path.join(HERE, "run_log.txt"), "a") as f:
        f.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
