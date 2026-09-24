#!/usr/bin/env python3
r"""
Condensate-Cherenkov energy loss + the confinement crossover (exploratory, 2026-09-03).

IDEA (this session): a parton (a coloured knot-fragment) moving through the density-feedback
condensate radiates into the MAGNITUDE (density / DE) mode, whose sound speed is c_s = c/phi
(Paper VII, P7:rem:cs_de). By Frank--Tamm this is Cherenkov radiation above threshold v > c_s.

The point is NOT the (uncomputable, coupling-set) absolute energy-loss rate. It is that the SAME
sound speed c_s=c/phi that sets the Cherenkov threshold ALSO sets the elastic<->radiative crossover
of the confining flux tube:
  * v < c_s : the tube deformation propagates AHEAD of the parton -> the tube stretches
              quasi-statically, energy stored elastically as STRING TENSION -> SNAPS BACK -> CONFINED.
  * v > c_s : the deformation cannot keep up -> a Mach cone forms, energy is RADIATED (Cherenkov),
              lost to the medium, the tube cannot follow -> escape / de-confining.
So the confinement/deconfinement KINEMATIC boundary is v = c_s = c/phi (derived). Temperature is the
second, independent axis (order/disorder of the condensate): high-T early universe = disordered
condensate, no room for a stable trefoil -> gluon soup; cools -> order sets in -> trefoil stabilises
-> confinement onset. (Temperature axis is the user's hypothesis, NOT derived here.)

Frank--Tamm in the condensate:  dE/dx = g_eff^2 * INT_{v>c_s} omega [1 - c_s^2/(v^2)] domega
  - threshold beta_th = c_s = 1/phi
  - Cherenkov factor  F(beta) = 1 - c_s^2/beta^2 = 1 - 1/(phi^2 beta^2)
  - cone angle        cos theta_c = c_s/v = 1/(phi beta)
The coupling g_eff and the UV cutoff (parton size / tube cross-section) set the MAGNITUDE and are
NOT derived (cutoff-set, CLAUDE.md sec.10). The kinematics below ARE derived from c_s=1/phi.
"""
import math

phi = (1 + 5**0.5) / 2
cs  = 1/phi                       # magnitude sound speed in units of c

print("="*74)
print("CONDENSATE-CHERENKOV + CONFINEMENT CROSSOVER   (c_s = c/phi, phi = %.6f)" % phi)
print("="*74)

# --- 1. Threshold (DERIVED) --------------------------------------------------
beta_th  = cs
gamma_th = 1/math.sqrt(1 - beta_th**2)
KE_over_mc2 = gamma_th - 1
print("\n[1] Cherenkov threshold  (a parton radiates into the magnitude mode above this):")
print(f"    beta_th = c_s = 1/phi          = {beta_th:.6f} c")
print(f"    gamma_th = 1/sqrt(1-1/phi^2)   = {gamma_th:.6f}")
print(f"    kinetic-energy threshold (gamma-1) m c^2 = {KE_over_mc2:.4f} m c^2")
print("    -> partons in hadrons are ultra-relativistic (beta~1) => DEEP in the radiative regime.")

# --- 2. The golden identity at ultra-relativistic speed (DERIVED, exact) -----
F_ur = 1 - 1/phi**2
print("\n[2] Ultra-relativistic Cherenkov factor  F = 1 - 1/(phi^2 beta^2):")
print(f"    beta -> 1 :  F -> 1 - 1/phi^2 = {F_ur:.6f}")
print(f"    EXACT identity: 1 - 1/phi^2 = 1/phi = {1/phi:.6f}   (since 1/phi^2 = 1 - 1/phi)")
print("    -> the maximal radiated fraction saturates at exactly 1/phi.")

# --- 3. Cherenkov cone angle (DERIVED) ---------------------------------------
theta_ur = math.degrees(math.acos(1/phi))     # cos th = 1/(phi*beta), beta->1
print("\n[3] Cherenkov cone (magnitude mode), cos theta_c = 1/(phi beta):")
print(f"    beta -> 1 :  theta_c = arccos(1/phi) = {theta_ur:.2f} deg  (from the velocity axis)")
print(f"                 Mach angle (from the wavefront) = {90-theta_ur:.2f} deg")

# --- 4. F(beta) and cone across the relativistic range (DERIVED shape) -------
print("\n[4] Energy-loss SHAPE dE/dx ∝ F(beta) (magnitude = coupling*cutoff, NOT derived):")
print("     beta     gamma    F=1-1/(phi^2 b^2)   theta_c(deg)   regime")
for beta in (0.55, 0.618, 0.65, 0.75, 0.85, 0.95, 0.99, 0.999):
    if beta <= cs:
        print(f"    {beta:5.3f}   {1/math.sqrt(1-beta**2):5.3f}    {'--- below threshold ---':>18}      --        ELASTIC (confining)")
    else:
        F = 1 - 1/(phi**2 * beta**2)
        th = math.degrees(math.acos(1/(phi*beta)))
        print(f"    {beta:5.3f}   {1/math.sqrt(1-beta**2):5.3f}    {F:16.4f}   {th:8.2f}       RADIATIVE (Cherenkov)")

# --- 5. Confinement crossover (interpretation) -------------------------------
print("\n[5] CONFINEMENT <-> RADIATION crossover is EXACTLY v = c_s = c/phi:")
print("    v < c/phi : deformation outruns nothing -> stored as string tension -> snaps back = CONFINED")
print("    v > c/phi : Mach cone -> energy radiated (Cherenkov) -> not recovered = DE-CONFINING")
print("    Cherenkov does NOT set the string-tension MAGNITUDE (that is the tube's elastic modulus x")
print("    cross-section, cutoff-set). It sets the KINEMATIC BOUNDARY of clean confinement, and the")
print("    ultra-relativistic radiated fraction 1/phi.")

print("\n[6] TWO-AXIS de-confinement map (speed AND temperature):")
print("    SPEED axis      : v>c/phi -> radiative escape (this calc, derived boundary).")
print("    TEMPERATURE axis: high-T -> condensate disordered -> no stable trefoil -> gluon soup;")
print("                      cools -> order + 'room' -> trefoil stabilises -> confinement onset.")
print("                      [user's hypothesis; NOT derived -- the trefoil-stability temperature is")
print("                       the QCD/charge-sector scale, not the DE Lambda_cond scale.]")

print("\n" + "="*74)
print("DERIVED (from c_s=1/phi): threshold beta=1/phi, gamma=%.3f; UR factor=1/phi; cone=%.1f deg." % (gamma_th, theta_ur))
print("NOT DERIVED: dE/dx magnitude (coupling g_eff, UV cutoff); string-tension value; T_deconf.")
print("="*74)
