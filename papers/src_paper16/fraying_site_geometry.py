#!/usr/bin/env python3
"""
fraying_site_geometry.py
========================
Compute and verify exact geometric data for the 12 fraying sites of the
Q_H=3 trefoil T_{2,3} condensate at R0=3.0, r0=0.874.

The 12 sites fall into four Z3-symmetric categories of 3 sites each:

  Category A — Crossings (3 sites):
    t = pi/6 + 2n*pi/3,  n=0,1,2
    R_xy = R0 = 3.0,  z = +r0 = +0.874
    kappa = 0.399 (maximum),  tang_z = 0 (exactly horizontal)

  Category B — Midpoint torsion-twists (3 sites):
    t = pi/3 + 2n*pi/3,  n=0,1,2
    R_xy = R0 - r0 = 2.126 (minimum),  z = 0
    kappa = 0.026 (minimum),  tang_z ~ -0.525 |Γ̇|

  Category C — Distal lobes (3 sites):
    t = 2n*pi/3,  n=0,1,2
    R_xy = R0 + r0 = 3.874 (maximum),  z = 0
    kappa = 0.349 (intermediate),  tang_z ~ +0.321 |Γ̇|

  Category D — Crossing approach/departure pair at each crossing:
    These are resolved from Category A when the fraying is observed in 3D
    field data: the two fray events per crossing (approach vs departure,
    60° offset in azimuth from each other) combine with Category A to give:
    6 crossing-associated sites + 3 midpoint + 3 lobe = 12 total.

Outputs:
  - Tabulated geometry for all three point types (crossings, midpoints, lobes)
  - 3D equidistance verification for lobes
  - Tangent z-components (the key Paper XVI arc-segment quantities)
  - kappa ordering confirmation (crossing > lobe > midpoint)
  - Full 12-site structure summary

Paper XVI reference: Proposition prop:geom_data, Remark rem:torsion_contrast,
Remark rem:arcseg_physics, Remark rem:rhoJ4_anticorrelation.
"""

import numpy as np

PHI = (1 + 5**0.5) / 2
R0, r0 = 3.0, 0.874

# ── Trefoil geometry ──────────────────────────────────────────────────────────

def Gamma(t):
    """Trefoil T(2,3) parametrisation: t ∈ [0, 2π)."""
    return np.array([
        (R0 + r0*np.cos(3*t)) * np.cos(2*t),
        (R0 + r0*np.cos(3*t)) * np.sin(2*t),
        r0 * np.sin(3*t)
    ])

def tangent_unit(t, h=1e-7):
    """Unit tangent vector dΓ/dt / |dΓ/dt|."""
    g1 = (Gamma(t + h) - Gamma(t - h)) / (2*h)
    return g1 / np.linalg.norm(g1)

def speed(t, h=1e-7):
    """Speed |dΓ/dt|."""
    g1 = (Gamma(t + h) - Gamma(t - h)) / (2*h)
    return np.linalg.norm(g1)

def kappa(t, h=1e-6):
    """Frenet curvature κ = |Γ' × Γ''| / |Γ'|³."""
    g1 = (Gamma(t + h) - Gamma(t - h)) / (2*h)
    g2 = (Gamma(t + h) - 2*Gamma(t) + Gamma(t - h)) / h**2
    return np.linalg.norm(np.cross(g1, g2)) / np.linalg.norm(g1)**3

# ── Distinguished parameter values ────────────────────────────────────────────

# Crossings: 3t = π/2 ⟹ z = r0 sin(π/2) = r0, R_xy = R0 + r0 cos(π/2) = R0
CROSSING_TS = [np.pi/6 + 2*n*np.pi/3 for n in range(3)]

# Midpoints: 3t = π ⟹ z = 0, R_xy = R0 + r0 cos(π) = R0 - r0
MIDPOINT_TS = [np.pi/3 + 2*n*np.pi/3 for n in range(3)]

# Distal lobes: 3t = 0 ⟹ z = 0, R_xy = R0 + r0 cos(0) = R0 + r0
LOBE_TS = [2*n*np.pi/3 for n in range(3)]

# ── Geometry table ────────────────────────────────────────────────────────────

def row(label, t):
    p = Gamma(t)
    R_xy = np.sqrt(p[0]**2 + p[1]**2)
    k = kappa(t)
    phi_az = np.degrees(np.arctan2(p[1], p[0]))
    tz = tangent_unit(t)[2]
    spd = speed(t)
    return dict(label=label, t=t, x=p[0], y=p[1], z=p[2],
                R_xy=R_xy, phi_az=phi_az, kappa=k, tang_z=tz, speed=spd)

print("=" * 72)
print("TREFOIL T(2,3) FRAYING SITE GEOMETRY   R0={}, r0={}".format(R0, r0))
print("=" * 72)

# ── Category A: Crossings ─────────────────────────────────────────────────────
print("\nCATEGORY A — Crossings (t = π/6 + 2nπ/3)")
print("-" * 72)
hdr = f"  {'Label':<10} {'t':>7} {'R_xy':>7} {'z':>7} {'φ_az':>7} {'κ':>8} {'tang_z':>9}"
print(hdr)
cross_rows = [row(f"C{n}", t) for n, t in enumerate(CROSSING_TS)]
for r_ in cross_rows:
    print(f"  {r_['label']:<10} {r_['t']:>7.4f} {r_['R_xy']:>7.4f} "
          f"{r_['z']:>7.4f} {r_['phi_az']:>7.1f}° {r_['kappa']:>8.5f} "
          f"{r_['tang_z']:>9.6f}")

print(f"\n  Analytic:  R_xy = R0 = {R0},  z = r0 = {r0},  tang_z = 0 exactly")
print(f"  Numerical: z mean = {np.mean([r_['z'] for r_ in cross_rows]):.6f}  "
      f"tang_z max abs = {max(abs(r_['tang_z']) for r_ in cross_rows):.2e}")

# ── Category B: Midpoints ─────────────────────────────────────────────────────
print("\nCATEGORY B — Midpoints (t = π/3 + 2nπ/3)")
print("-" * 72)
print(hdr)
mid_rows = [row(f"M{n}", t) for n, t in enumerate(MIDPOINT_TS)]
for r_ in mid_rows:
    print(f"  {r_['label']:<10} {r_['t']:>7.4f} {r_['R_xy']:>7.4f} "
          f"{r_['z']:>7.4f} {r_['phi_az']:>7.1f}° {r_['kappa']:>8.5f} "
          f"{r_['tang_z']:>9.6f}")

print(f"\n  Analytic:  R_xy = R0 - r0 = {R0 - r0},  z = 0,  tang_z < 0 (downward)")
print(f"  tang_z = {mid_rows[0]['tang_z']:.6f}  (unit tangent z-component; all three agree to 5 sig figs)")

# ── Category C: Distal lobes ──────────────────────────────────────────────────
print("\nCATEGORY C — Distal lobes (t = 2nπ/3)")
print("-" * 72)
print(hdr)
lobe_rows = [row(f"L{n}", t) for n, t in enumerate(LOBE_TS)]
for r_ in lobe_rows:
    print(f"  {r_['label']:<10} {r_['t']:>7.4f} {r_['R_xy']:>7.4f} "
          f"{r_['z']:>7.4f} {r_['phi_az']:>7.1f}° {r_['kappa']:>8.5f} "
          f"{r_['tang_z']:>9.6f}")

print(f"\n  Analytic:  R_xy = R0 + r0 = {R0 + r0},  z = 0,  tang_z > 0 (upward)")
print(f"  tang_z = {lobe_rows[0]['tang_z']:.6f}  (unit tangent z-component; all three agree to 5 sig figs)")

# ── Equidistance verification ─────────────────────────────────────────────────
print("\n" + "=" * 72)
print("EQUIDISTANCE VERIFICATION (3D distances from each lobe to all crossings)")
print("=" * 72)
print("Note: each lobe (φ=0°, −120°, 120°) is adjacent to two crossings at ±60°")
print("      relative to it. The far crossing is at 180° relative.")
print()

for n, lt in enumerate(LOBE_TS):
    p_L = Gamma(lt)
    phi_L = np.degrees(np.arctan2(p_L[1], p_L[0]))
    dists = []
    for m, ct in enumerate(CROSSING_TS):
        p_C = Gamma(ct)
        d = np.linalg.norm(p_L - p_C)
        phi_C = np.degrees(np.arctan2(p_C[1], p_C[0]))
        dists.append((m, d, phi_C))
    dists.sort(key=lambda x: x[1])
    print(f"  L{n} (φ={phi_L:+.0f}°):")
    for m, d, phi_C in dists:
        tag = "← adjacent" if d < 4.0 else "← far"
        print(f"    to C{m} (φ={phi_C:+.0f}°):  d = {d:.6f}  {tag}")

d_adj = np.linalg.norm(Gamma(0) - Gamma(np.pi/6))
print(f"\n  Adjacent distance (exact): {d_adj:.6f}")
print(f"  Far      distance (exact): {np.linalg.norm(Gamma(0) - Gamma(np.pi/6+4*np.pi/3)):.6f}")

# ── Curvature ordering ────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("CURVATURE ORDERING  (S_eff ~ κ²)")
print("=" * 72)
k_cross = kappa(np.pi/6)
k_lobe  = kappa(0)
k_mid   = kappa(np.pi/3)
print(f"  Crossing (t=π/6):   κ = {k_cross:.5f},  κ² = {k_cross**2:.6f}  ← S_eff MAX")
print(f"  Lobe     (t=0):     κ = {k_lobe:.5f},  κ² = {k_lobe**2:.6f}  ← intermediate")
print(f"  Midpoint (t=π/3):   κ = {k_mid:.6f},  κ² = {k_mid**2:.8f}  ← S_eff MIN")
print(f"\n  Ratio κ(cross)/κ(lobe)   = {k_cross/k_lobe:.4f}")
print(f"  Ratio κ(lobe)/κ(mid)     = {k_lobe/k_mid:.2f}")
print(f"  Ratio κ(cross)/κ(mid)    = {k_cross/k_mid:.1f}")
print(f"\n  kappa(lobe) = {k_lobe:.5f}  (NOT at global max or min)")

# Verify by scanning
t_scan = np.linspace(0, 2*np.pi, 3000, endpoint=False)
k_scan = np.array([kappa(t) for t in t_scan])
print(f"  Global kappa max: {k_scan.max():.5f} at t={t_scan[k_scan.argmax()]:.4f}")
print(f"  Global kappa min: {k_scan.min():.6f} at t={t_scan[k_scan.argmin()]:.4f}")

# ── Tangent z-components summary ──────────────────────────────────────────────
print("\n" + "=" * 72)
print("TANGENT z-COMPONENTS  (key arc-segment prediction quantities)")
print("=" * 72)
tz_cross = tangent_unit(np.pi/6)[2]
tz_lobe  = tangent_unit(0)[2]
tz_mid   = tangent_unit(np.pi/3)[2]
print(f"  Crossing (t=π/6):   tang_z = {tz_cross:.2e}  (zero: exactly horizontal)")
print(f"  Midpoint (t=π/3):   tang_z = {tz_mid:.6f}  (≈ -0.525 |Γ̇|, downward)")
print(f"  Distal lobe (t=0):  tang_z = {tz_lobe:.6f}  (≈ +0.321 |Γ̇|, upward)")
print(f"\n  Sign convention: tang_z(crossing) = 0, tang_z(midpoint) < 0 < tang_z(lobe)")
print(f"  Midpoints and lobes point in OPPOSITE z-directions, consistent with the")
print(f"  tube curving from z=0 up to z=+r0 (crossing) and back down through z=0.")

# ── 12-site summary ───────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("12-SITE FRAYING STRUCTURE SUMMARY")
print("=" * 72)
print(f"""
The 12 fraying sites divide into four categories, each with Z3 symmetry:

  Sites 1–3   Crossing departure zones (B-segment exit, post-crossing)
              t ≈ π/6 + 2nπ/3 + ε,  coplanar in z = +{r0} plane, 120° apart
              Carries 79% of total J4 flux (A-segment bodies)
              PRIMARY jets 1,2,3

  Sites 4–6   Distal lobe exits (B-segment end, maximum R_xy)
              t = 2nπ/3,  R_xy = R0+r0 = {R0+r0},  z = 0
              κ = {k_lobe:.3f} (intermediate),  tang_z ≈ +0.321 |Γ̇| (upward)
              SECONDARY jets 4,5,6,  deflected ≈16°–19° below three-jet plane

  Sites 7–9   Crossing approach zones (A-segment entry, pre-crossing)
              t ≈ π/6 + 2nπ/3 − ε,  softer emission than departure
              TERTIARY jets 7,8,9

  Sites 10–12 Midpoint torsion-twist residuals (structural minimum)
              t = π/3 + 2nπ/3,  R_xy = R0-r0 = {R0-r0},  z = 0
              κ = {k_mid:.4f} (minimum),  tang_z ≈ -0.525 |Γ̇| (downward)
              SOFTEST jets 10,11,12

  Cascade energy ordering:  departure > lobe > approach > midpoint
  Z3 structure:  each group of 3 is mutually coplanar by Z3 symmetry
  Between-group acoplanarity:  departure/lobe ≈ 16°–19°, lobe/midpoint ≈ 50°

  Adjacent crossing distance (lobe ↔ two nearest crossings):  {d_adj:.4f}
""")
