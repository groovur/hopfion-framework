#!/usr/bin/env python3
r"""
Static rotational susceptibility chi of the semi-Dirac condensate -- the diamagnetic
(q=0) term of frank_run.py, isolated. Enters c_dir^2 = K/chi and the polaron effective
mass m_eff = (chi/3) INT|grad n|^2 (Paper VII P7:rem:polaron, P7:eq:meff).

Same band + generator as frank_run.py:
  d = (kx^2/2, ky),  E_pm = +-|d|,  lower band |-(k)>=(1,-e^{i phi})/sqrt2, phi=atan2(d2,d1).
  generator G=(kx ky, -kx),  M(k,k) = <+(k)|G.sigma|-(k)>.
Static (van Vleck / diamagnetic) susceptibility:
  chi = INT_disk d^2k/(2pi)^2  |M(k,k)/2|^2 / (2|d_k|)      (positive).
This is frank_run's SUBTRACTION term with sign flipped -- q-INDEPENDENT, hence a single
value for generator G (no splay/bend split at q=0 from this generator alone).

To probe a possible chi anisotropy we also evaluate a SECOND generator (the orthogonal
director-rotation DOF), G2=(ky, kx*ky) [rotation about the other in-plane axis], as a
first check of whether chi_splay =? chi_bend. This second generator is a MODELLING choice
(the framework's exact second generator not confirmed here) -- flagged accordingly.
"""
import numpy as np, sys

def chi_integral(Lam, Nr, Nth, gen="G1"):
    rs = np.linspace(Lam/Nr, Lam, Nr); dr = rs[1]-rs[0]
    th = np.linspace(0, 2*np.pi, Nth, endpoint=False); dth = th[1]-th[0]
    C, S = np.cos(th), np.sin(th)
    acc = 0.0
    for r in rs:
        kx, ky = r*C, r*S
        d1, d2 = kx**2/2, ky
        nd = np.hypot(d1, d2)
        phi = np.arctan2(d2, d1)
        if gen == "G1":
            G1, G2 = kx*ky, -kx           # frank_run's generator
        else:                              # G2: orthogonal rotation DOF (modelling choice)
            G1, G2 = ky, kx*ky
        # M(k,k) = -(G1 - i G2) e^{i phi} + (G1 + i G2) e^{-i phi}
        z0 = -(G1 - 1j*G2)*np.exp(1j*phi) + (G1 + 1j*G2)*np.exp(-1j*phi)
        integrand = (np.abs(z0)/2)**2 / (2*nd)
        acc += np.sum(integrand) * r*dr*dth / (2*np.pi)**2
    return float(acc)

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "test"
    if mode == "test":
        Lam, Nr, Nth = 20.0, 400, 720
    else:
        Lam, Nr, Nth = 40.0, 2000, 2880
    chi1 = chi_integral(Lam, Nr, Nth, "G1")
    chi2 = chi_integral(Lam, Nr, Nth, "G2")
    print(f"[frank_chi] Lam={Lam} Nr={Nr} Nth={Nth}")
    print(f"  chi(G1, frank_run generator)      = {chi1:.6e}")
    print(f"  chi(G2, orthogonal DOF [model])   = {chi2:.6e}")
    print(f"  chi(G2)/chi(G1)                   = {chi2/chi1:.4f}  "
          f"(=1 -> chi isotropic -> orientation effect is CONFIG-only, chi-free)")
    print("  NOTE: G1 is the validated generator; G2 is a modelling guess for the 2nd rotation DOF.")
    print("  If chi is isotropic, the oblate aligned-vs-perp comparison is chi-free (config ratio).")

if __name__ == "__main__":
    main()
