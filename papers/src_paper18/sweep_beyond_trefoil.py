#!/usr/bin/env python3
"""
sweep_beyond_trefoil.py
===================================================================
Orchestrates hopf_link_construction_v2.py --verify_only over a sweep
of q_pol from 2 (Hopf link) through 3 (trefoil) and on to 5
(cinquefoil), at fixed alpha_wind=1.0 (full Q_H=3-type winding
commitment), to check whether the constructed field's energy
E_geom = K*J4 shows a barrier between the trefoil and the cinquefoil
-- i.e. whether there is independent, non-Jones-polynomial evidence
for or against the hypothesis that the trefoil is a waypoint on a
continuous deformation that can proceed further toward T(2,5) given
enough energy, rather than a fully independent, separate object.

This does NOT run gradient flow (no relaxation) -- it only evaluates
the energy of the constructed initial ansatz at each q_pol, which is
fast and a reasonable first check. Whether the (unrelaxed) ansatz
energy's q_pol-dependence tracks the true (relaxed) energy landscape
is itself an assumption, not verified here.

Output: data/sweep_beyond_trefoil/summary.json and a printed table.
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'data', 'sweep_beyond_trefoil')
os.makedirs(OUTDIR, exist_ok=True)

Q_VALUES = [2.0, 2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]

results = []
for q in Q_VALUES:
    print(f"\n{'='*70}\nq_pol = {q}\n{'='*70}")
    cmd = [
        sys.executable, os.path.join(HERE, 'hopf_link_construction_v2.py'),
        '--N', '64', '--q_pol', str(q), '--alpha_wind', '1.0',
        '--verify_only', '--outdir', OUTDIR, '--whitehead_N', '32',
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, timeout=600)
    out = proc.stdout
    if proc.returncode != 0:
        print(out)
        print(proc.stderr)
        raise RuntimeError(f"run failed at q_pol={q}")

    # Parse the FIRST "Whitehead Q_H = ... / E=... J4=... K=... r_bar=..." block
    # (the first block is always our target q_pol construction; the second
    # block, later in the output, is always the reference T(2,3) trefoil).
    m = re.search(r"Whitehead Q_H = ([\d.eE+-]+)[^\n]*\n\s+E=([\d.eE+-]+)\s+J4=([\d.eE+-]+)\s+K=([\d.eE+-]+)\s+r_bar=([\d.eE+-]+)", out)
    if not m:
        print(out)
        raise RuntimeError(f"could not parse output at q_pol={q}")
    QH, E, J4, K, rbar = (float(x) for x in m.groups())
    results.append(dict(q_pol=q, QH=QH, E=E, J4=J4, K=K, r_bar=rbar))
    print(f"  parsed: QH={QH:.3f} E={E:.4e} J4={J4:.1f} K={K:.1f} r_bar={rbar:.3f}")

with open(os.path.join(OUTDIR, 'summary.json'), 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "=" * 70)
print("SUMMARY: constructed-field energy vs q_pol (alpha_wind=1.0)")
print("=" * 70)
print(f"{'q_pol':>6} {'E':>14} {'J4':>10} {'K':>10} {'r_bar':>8}")
for r in results:
    print(f"{r['q_pol']:>6.2f} {r['E']:>14.4e} {r['J4']:>10.1f} {r['K']:>10.1f} {r['r_bar']:>8.3f}")

Es = [r['E'] for r in results]
qs = [r['q_pol'] for r in results]
imin, imax = Es.index(min(Es)), Es.index(max(Es))
print(f"\nmin E at q_pol={qs[imin]} (E={Es[imin]:.4e})")
print(f"max E at q_pol={qs[imax]} (E={Es[imax]:.4e})")
print("\n(Look for a local maximum strictly between q_pol=3 and q_pol=5")
print(" -- that would be the energy-barrier signature being tested.)")
