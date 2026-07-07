#!/usr/bin/env python3
"""
sweep_gradient_flow.py
===================================================================
Runs actual gradient flow (not --verify_only) on hopf_link_construction_v2.py
at a set of q_pol starting points spanning the trefoil (q=3) to the
cinquefoil (q=5), to check whether relaxation settles to a stable
configuration (energy plateaus, Q_H stabilizes) or keeps drifting
(dilution / instability), and whether any such stability correlates
with phi-related values -- checked AFTER the run, not presupposed.

This is real computation (~12-13 min per point on CPU, since each
invocation also builds+relaxes the reference T(2,3) trefoil for
comparison, per the existing script's structure). Points chosen from
the prior --verify_only energy sweep: q=3.0 (trefoil baseline),
q=3.5 (energy peak), q=4.0 (energy trough / Q_H=4 sector), q=4.5,
q=5.0 (cinquefoil).

Output: data/gradient_flow_sweep/<q>/... (per-point snapshots + log.json)
        data/gradient_flow_sweep/summary.json
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'data', 'gradient_flow_sweep')
os.makedirs(OUTDIR, exist_ok=True)

Q_VALUES = [3.0, 3.5, 4.0, 4.5, 5.0]
N_STEPS = 2400
WHITEHEAD_N = 64

results = []
for q in Q_VALUES:
    point_dir = os.path.join(OUTDIR, f'q{q:.2f}')
    print(f"\n{'='*70}\nq_pol = {q}  (gradient flow, {N_STEPS} steps)\n{'='*70}", flush=True)
    cmd = [
        sys.executable, os.path.join(HERE, 'hopf_link_construction_v2.py'),
        '--N', '64', '--q_pol', str(q), '--alpha_wind', '1.0',
        '--n_steps', str(N_STEPS), '--snap_every', '200', '--log_every', '20',
        '--whitehead_N', str(WHITEHEAD_N), '--outdir', point_dir,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, timeout=1800)
    out = proc.stdout
    if proc.returncode != 0:
        print(out)
        print(proc.stderr)
        raise RuntimeError(f"run failed at q_pol={q}")

    with open(os.path.join(point_dir, 'run_stdout.txt'), 'w') as f:
        f.write(out)

    # Parse both blocks (target construction, then reference trefoil):
    # "Initial: E=... K=... J4=... r_bar=..." and
    # "Final: E=... J4/K=... r_bar=..." and "Whitehead Q_H: X -> Y"
    blocks = re.findall(
        r"RUN: (.+?)\n.*?Whitehead Q_H = (-?[\d.]+).*?\n\s*Initial: E=([\d.eE+-]+)\s+K=([\d.eE+-]+)\s+J4=([\d.eE+-]+)\s+r_bar=([\d.eE+-]+).*?"
        r"Final: E=([\d.eE+-]+)\s+J4/K=([\d.eE+-]+)\s+r_bar=([\d.eE+-]+)\n\s*Whitehead Q_H: (-?[\d.]+) . (-?[\d.]+)",
        out, re.DOTALL)
    parsed = []
    for b in blocks:
        label, qh0, E0, K0, J40, rb0, Ef, J4K_f, rbf, qh_i, qh_f = b
        parsed.append(dict(label=label, QH_init=float(qh0), E0=float(E0),
                            K0=float(K0), J40=float(J40), rbar0=float(rb0),
                            Ef=float(Ef), J4K_final=float(J4K_f), rbar_f=float(rbf),
                            QH_before=float(qh_i), QH_after=float(qh_f)))
    if len(parsed) != 2:
        print(out)
        raise RuntimeError(f"expected 2 blocks (construction+trefoil), got {len(parsed)} at q_pol={q}")

    entry = dict(q_pol=q, construction=parsed[0], trefoil_ref=parsed[1])
    results.append(entry)
    c = parsed[0]
    print(f"  construction: E {c['E0']:.4e} -> {c['Ef']:.4e}  "
          f"(delta {(c['Ef']-c['E0'])/c['E0']*100:+.2f}%)  "
          f"QH {c['QH_before']:.3f} -> {c['QH_after']:.3f}  r_bar_f={c['rbar_f']:.3f}")

with open(os.path.join(OUTDIR, 'summary.json'), 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "=" * 70)
print("SUMMARY: relaxation behavior, q_pol = 3 -> 5")
print("=" * 70)
print(f"{'q_pol':>6} {'E0':>12} {'Ef':>12} {'%change':>9} {'QH0':>7} {'QHf':>7} {'r_bar_f':>9}")
for r in results:
    c = r['construction']
    pct = (c['Ef'] - c['E0']) / c['E0'] * 100
    print(f"{r['q_pol']:>6.2f} {c['E0']:>12.4e} {c['Ef']:>12.4e} {pct:>8.2f}% "
          f"{c['QH_before']:>7.3f} {c['QH_after']:>7.3f} {c['rbar_f']:>9.3f}")
print("\n(Small |%change| in E and QH => relaxation is near a stable point.")
print(" Large drift or rbar_f blowing up (>4.5 triggers HALT) => unstable/diluting.)")
