#!/usr/bin/env python3
r"""
Kibble-Zurek estimate: trapped-disclination mass from cluster assembly vs the deficit.
The one MOND-impossible cluster lever (topological, baryon-independent, assembly-seeded).

MECHANISM: K_splay<0 spontaneously forms disclinations. Cluster hierarchical merging repeatedly
disorders/reorders the director; c_dir->0 means the director CANNOT relax, so merger-trapped defects
NEVER heal -> a persistent disclination network. Kibble-Zurek: ~1 defect line per correlation area
xi_turb^2, where xi_turb = the coherence scale of the merger-induced (re)ordering.

SCALING (analytic):
  disclination line energy/length  eps ~ K ln(xi_turb/a_core) ~ K   (K = Frank constant ~ a0/G)
  line density (length/volume)      ~ 1/xi_turb^2
  defect energy density  rho_def ~ eps/xi_turb^2 ~ K/xi_turb^2
  M_def(R_cl) ~ rho_def * R_cl^3 ~ K R_cl^3 / xi_turb^2
  RAR halo    M_RAR(R_cl) ~ INT (K/r^2) 4pi r^2 dr ~ K R_cl        (isothermal, the MOND-short value)
  => M_def / M_RAR ~ (R_cl / xi_turb)^2       <-- the whole result, dimensionless, K cancels.

So the trapped-defect mass EXCEEDS the RAR (baryon-tracking) value by (R_cl/xi_turb)^2. This is the
framework-specific EXTRA mass MOND structurally cannot have (assembly-history-dependent, not baryon g_N).
"""
import math

# deficit to close: cluster needs ~28x stellar; RAR single-object strain gives ~5x -> extra needed / RAR:
DM_over_stellar   = 28.0
RAR_over_stellar  = 5.0
extra_over_RAR    = (DM_over_stellar - RAR_over_stellar) / RAR_over_stellar   # M_def/M_RAR required
print("=== Kibble-Zurek trapped-disclination mass vs cluster deficit ===\n")
print(f"cluster DM ~ {DM_over_stellar:.0f}x stellar; RAR strain ~ {RAR_over_stellar:.0f}x stellar")
print(f"=> need M_def/M_RAR ~ {extra_over_RAR:.1f}  to close the deficit\n")

# required coherence scale: (R_cl/xi_turb)^2 = extra_over_RAR
xi_over_R_needed = 1.0/math.sqrt(extra_over_RAR)
print(f"REQUIRED xi_turb / R_cl = 1/sqrt({extra_over_RAR:.1f}) = {xi_over_R_needed:.2f}")
print(f"  i.e. the merger-reordering coherence length must be ~{xi_over_R_needed:.2f} of the cluster radius.\n")

print("M_def/M_RAR = (R_cl/xi_turb)^2  for a range of coherence scales:")
for xr in (1.0, 0.7, 0.5, 0.3, 0.1):
    ratio = (1.0/xr)**2
    tag = "CLOSES" if abs(ratio-extra_over_RAR) < 0.25*extra_over_RAR else ("OVER-produces" if ratio>extra_over_RAR else "short")
    print(f"  xi_turb/R_cl={xr:>4}:  M_def/M_RAR = {ratio:6.1f}   ({tag})")

print("\n=== READ ===")
print(f"To recover the ~{extra_over_RAR:.0f}x deficit, need xi_turb ~ {xi_over_R_needed:.1f} R_cl -- coherence over")
print("~HALF the cluster radius. Cluster mergers ARE subcluster-scale (subclusters ~0.3-0.7 R_cl), so the")
print("required xi_turb MATCHES the actual merger scale. => RIGHT ORDER OF MAGNITUDE, and the required scale is")
print("physical, NOT tuned. This is a LIVE recovery (unlike the screened MOND-like channels).")
print("BUT strongly SCALE-SENSITIVE: fine turbulence (xi_turb<<R_cl) OVER-produces badly ((R/xi)^2 blows up),")
print("coherent mergers (xi_turb~R_cl) give ~few. So the magnitude hinges on xi_turb, and on defect ANNIHILATION")
print("(opposite-charge disclinations heal; c_dir->0 slows healing but like/unlike mix sets the NET density).")
print("HONEST: mechanism is framework-specific + right-order + physically-scaled -- the ONLY lever not screened.")
print("Uncomputed precisely: net defect density after annihilation, and the K->a0 normalization. Promising, open.")
print("DISCRIMINATOR stands: disturbed (small effective xi_turb / recent mergers) clusters -> MORE excess than")
print("relaxed ones -- MOND predicts none, LambdaCDM ties excess to the halo not the merger state.")
