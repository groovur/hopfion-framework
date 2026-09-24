#!/usr/bin/env python3
r"""
Bounded-motion Cherenkov: why a CONFINED quark does not continuously bleed energy (2026-09-03).

The straight-line result (cherenkov_confinement.py): a parton at v>c_s=c/phi radiates into the
magnitude mode. Hadronic partons are ultra-relativistic (v~c > c/phi), so naively they should
radiate forever -- but hadrons don't bleed energy. Resolution (this session):

A confined quark is NOT in free translation; it oscillates INSIDE the wound trefoil tube. The tube
is a WAVEGUIDE for the magnitude ripples the quark emits. A waveguide of transverse size R passes
only modes shorter than a cutoff wavelength lambda_c ~ 2R; longer modes are EVANESCENT (trapped).
So:
  * emission BELOW cutoff  -> evanescent -> trapped as standing modes = the hadron's INTERNAL
    excitation spectrum, NOT radiated to infinity. No net loss. (= "directed by the wound medium".)
  * emission ABOVE cutoff  -> propagates out -> leakage = hadronic widths / transitions.
  * NET TRANSLATION (escape attempt, v>c/phi) -> straight-line Cherenkov -> radiates ->
    ESCAPE DAMPING (the bounded motion is 'also damped', as suspected, but only on escape).

The point below: the confinement scale and the waveguide cutoff are the SAME scale (both ~ c_s/R),
so the confined quark sits WITHIN A FACTOR ~phi of cutoff -> MARGINAL trapping. That is not a
coincidence to tune away; it is why hadrons are long-lived (trapped) yet have finite widths and
radiative transitions (marginal leakage). Magnitudes are cutoff-set (CLAUDE.md sec.10); the SCALING
and the marginality are the robust, derivable content.

Channel separation (connects to the 'photon' thought): the magnitude/sound ripples are density
perturbations of the SAME condensate that forms the tube, so they couple to it and are waveguided.
The photon (a colour-NEUTRAL easy-radial excitation, v=c) does NOT couple to the colour/director
tube, so it escapes freely -- which is why the hadron radiates real PHOTONS (radiative decays) while
its sound-channel ripples stay trapped. (This is the correct reason; a naive waveguide-cutoff
argument fails, since the faster EM mode has the HIGHER cutoff E_c=h c/2R = phi * E_c^sound, i.e. it
would be trapped MORE easily, not less -- the escape is by decoupling, not by cutoff.)
"""
import math

hbarc = 197.327            # MeV*fm
phi   = (1+5**0.5)/2
cs    = 1/phi              # magnitude sound speed / c

print("="*76)
print("BOUNDED-MOTION CHERENKOV: the tube as a waveguide for magnitude ripples")
print("c_s = c/phi = %.4f c ;  hbar c = %.2f MeV*fm" % (cs, hbarc))
print("="*76)

def cutoff_E(R):           # waveguide cutoff energy E_c = 2*pi*hbar*c_s/lambda_c, lambda_c~2R
    lam_c = 2*R
    return 2*math.pi*hbarc*cs/lam_c    # MeV

def quark_E(R):            # particle-in-a-tube ground scale ~ pi*hbar/R (ultra-rel quark), MeV
    return math.pi*hbarc/R

print("\n[1] Waveguide cutoff of the tube (lambda_c ~ 2R) vs the confined-quark scale (~pi hbar/R):")
print("      R(fm)   lambda_emit@300MeV(fm)   E_cutoff(MeV)   E_quark(MeV)   ratio  regime")
for R in (0.5, 0.75, 1.0, 1.25, 1.5):
    lam_emit = 2*math.pi*hbarc*cs/300.0     # wavelength of a 300 MeV magnitude ripple
    Ec = cutoff_E(R); Eq = quark_E(R)
    ratio = Eq/Ec
    reg = "trapped (below cutoff)" if ratio<1 else "marginal/leaky"
    print(f"    {R:5.2f}   {lam_emit:8.2f}                {Ec:8.0f}       {Eq:8.0f}     {ratio:5.2f}  {reg}")

print("\n[2] Why marginal is STRUCTURAL, not tuned:")
print("    E_cutoff  = 2*pi*hbar*c_s/(2R) = pi*hbar*c_s/R   (waveguide, ~ c_s/R)")
print("    E_quark   ~ pi*hbar/R                            (particle-in-tube, ~ c/R)")
print("    ratio E_quark/E_cutoff = 1/c_s = phi = %.4f   -> INDEPENDENT of R." % phi)
print("    => the confined quark always sits a factor phi ABOVE the sound-cutoff: marginal.")
print("       Soft internal modes (below cutoff) are trapped; the hardest sit just above => leak")
print("       as finite widths. Both are the SAME scale c/R -- no coincidence to remove.")

print("\n[3] The three regimes of the confined quark's magnitude emission:")
print("    (a) bounded, below cutoff : EVANESCENT -> standing modes = internal hadron spectrum. No loss.")
print("    (b) bounded, above cutoff : leaks out  -> hadronic widths / sound-channel transitions.")
print("    (c) net escape (v>c/phi)  : straight-line Cherenkov -> ESCAPE DAMPING (opposes de-confining).")
print("    -> continuous Cherenkov bleed is REPLACED by (a) trapping + (c) escape damping.")
print("       The Frank-Tamm cone (needs unbounded coherent translation) does NOT apply to (a).")

print("\n[4] Channel separation (the 'photon' connection) -- escape is by DECOUPLING, not cutoff:")
print("    magnitude/sound ripples: density modes of the SAME condensate as the tube -> couple ->")
print("                             waveguided/trapped (stay inside the hadron).")
print("    photon (colour-NEUTRAL) : does NOT couple to the colour/director tube -> escapes freely.")
print("    (A cutoff argument would say the opposite: E_c^EM = phi * E_c^sound, so the faster EM mode")
print("     is trapped MORE easily -- so the photon's escape must be decoupling, not kinematics.)")
print("    => a bounded, accelerating charged quark radiates real photons (radiative transitions)")
print("       while its sound-channel ripples remain trapped: the hadron glows but does not hiss.")

print("\n" + "="*76)
print("DERIVED (scaling): cutoff ~ c_s/R; E_quark/E_cutoff = phi (R-independent) => marginal trapping.")
print("                   bounded-below-cutoff = trapped (no loss); escape = Cherenkov-damped.")
print("REASONED [I]:      photon escapes by COLOUR-DECOUPLING from the tube; sound ripples couple/trapped.")
print("NOT DERIVED:       absolute widths / coupling; the exact waveguide cutoff constant (~2R).")
print("="*76)
