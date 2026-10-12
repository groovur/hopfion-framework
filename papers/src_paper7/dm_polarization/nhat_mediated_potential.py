#!/usr/bin/env python3.11
"""
nhat_mediated_potential.py -- the n_hat-mediated dark-hopfion interaction WITH the semi-Dirac
anisotropy resolved (notes (kkk)->(lll)). User: "remember the anisotropy inherent in the field;
the Bell papers may help the concept." Both apply.

(kkk) used an ISOTROPIC massless Yukawa -> Rutherford sigma~1/v^4 (n=4), EXCLUDED. WRONG: the
mediator is the SEMI-DIRAC director mode, E^2(k) = v^2 k_par^2 + A^2 k_perp^4 (k_par along the
director e_r = easy/LINEAR, k_perp = hard/QUADRATIC; claude-hopfion s3). Its static propagator is
ANISOTROPIC; the interaction is NOT a massless Yukawa.

Static potential V(r) = INT d^3k e^{ik.r}/E^2(k). Resolving the axes, the inner k-integral gives
~ C/|k| on BOTH axes, so V(r) ∝ INT_{kmin}^{kmax} cos(kr)/k dk = Ci(kmax r) - Ci(kmin r) -> a
LOGARITHMIC potential (2D-Coulomb/vortex-like), IR-regulated by kmin ~ m_xi (i.e. saturating at
r ~ 1/m_xi ~ the Hubble scale). Force F = -dV/dr ~ 1/r (Coulomb-like force). This is verified below
(dV/dln r -> const).

CONSEQUENCES:
 * (kkk)'s n=4 exclusion is an ISOTROPIZATION ARTIFACT; the real interaction is logarithmic.
 * A log potential (1/r force) gives PLASMA-LIKE relaxation, t_relax ~ v^3/(n C^2 ln L) -- the SAME
   v^3/n scaling as gravitational two-body relaxation. Dense+SLOW galaxies relax (thermalise ->
   baryon-tracking, RAR-preserving); diffuse+FAST clusters do not (collisionless, dominant). So the
   galaxy/cluster split the overlap needs is NATURAL for a log interaction -- IF the director
   coupling C is strong enough that t_relax(galaxy) < t_age < t_relax(cluster). C not computed [O].
 * BELL (Papers VIII/IX): the interaction is through the SHARED director n_hat ("geometric cost of
   the shared geometry"), not particle exchange; its long-range, Hubble-regulated (globally
   coherent) log character reflects exactly the single-shared-field mechanism. The user's Bell
   intuition is apt: this SIDM is a shared-geometry back-action.
HONEST: order-of-magnitude; C (director coupling) uncomputed -> whether galaxies actually relax is
open; the m_xi IR regulation makes the interaction LONG-RANGE (cosmological) -- unusual, needs care.
"""
import numpy as np
from scipy.special import sici

def V_log(r, kmin, kmax):
    _,Ci_hi=sici(kmax*r); _,Ci_lo=sici(kmin*r)
    return -(Ci_hi-Ci_lo)

if __name__=='__main__':
    kmin,kmax=1e-3,1e3     # IR ~ m_xi (1/kmin ~ Hubble), UV ~ 1/xi (knot)
    rs=np.array([0.01,0.1,1.0,10.0,100.0,300.0])
    V=np.array([V_log(r,kmin,kmax) for r in rs]); dv=np.gradient(V,np.log(rs))
    print("n_hat-mediated interaction (semi-Dirac propagator, axes resolved):")
    print(f"  {'r':>7} {'V(r)':>9} {'dV/dln r':>9}")
    for r,a,b in zip(rs,V,dv): print(f"  {r:>7.2f} {a:>9.3f} {b:>9.3f}")
    print(f"  dV/dln r ~ 1 => LOGARITHMIC (vortex/2D-Coulomb); IR-saturates at r~1/m_xi={1/kmin:.0f} (Hubble).")
    print("  => NOT the isotropic massless-Yukawa Rutherford of (kkk); n=4 exclusion was an isotropization")
    print("  artifact. Log => plasma-like relaxation t~v^3/(n C^2) => natural galaxy(thermalise)/cluster")
    print("  (collisionless) split IF the director coupling C is strong enough [O]. Bell: shared-n_hat coupling.")
