"""
Pin the cutoff: replace the CIRCULAR momentum cutoff exp(-2|k|^2/Lam^2) (which is
ambiguous as an energy for the anisotropic semi-Dirac band -> the 2.7x light/heavy
bracket) with the physically correct ISO-ENERGY cutoff exp(-2(E(k)/m_xi)^2). That
is single-valued: m_xi is a true energy, so m_xi/E* is unambiguous and the bracket
collapses to one number. Then invert band_K(m_xi)=target for the REQUIRED m_xi/E*
and compare to phi^6.

Validation: the same Cartesian integrator with the circular cutoff must reproduce
cutoff_embedding.py's polar result K_bend(Lam=10,20,40) = 0.1189, 0.3017, 0.7900.

Band: H=(kx^2/2)sx + ky sy ; d1=kx^2/2, d2=ky, E=sqrt(d1^2+d2^2), phi=atan2(d2,d1).
G=(kx ky, -kx). z0 = -(G1-iG2)e^{i phi}+(G1+iG2)e^{-i phi}.
K_bend (q||y): 2/q^2 * sum f^2 * sum_s 0.25[|z|^2/4/(-E-Ep) - |z0|^2/4/(-2E)].
"""
import numpy as np
phi = (1+5**0.5)/2
Estar = np.sqrt(8.0)                 # crossover energy (band units)
q = 0.05

def band_K_bend(cutoff, scale, dk=0.15, ext_ky=None, ext_kx=None):
    """cutoff in {'circ','ener'}; scale = Lam (circ) or m_xi (ener)."""
    if cutoff == 'circ':
        kxm = kym = 2.6*scale
    else:  # energy cutoff: ky_max ~ 2.6 m_xi (light), kx_max ~ sqrt(2*2.6 m_xi) (heavy)
        kym = ext_ky if ext_ky else 2.6*scale
        kxm = ext_kx if ext_kx else np.sqrt(2*2.6*scale)
    kx1 = np.arange(-kxm+dk/2, kxm, dk)
    ky1 = np.arange(-kym+dk/2, kym, dk)
    kx, ky = np.meshgrid(kx1, ky1)
    d1, d2 = kx**2/2, ky
    E = np.hypot(d1, d2); ph = np.arctan2(d2, d1)
    G1, G2 = kx*ky, -kx
    z0 = -(G1-1j*G2)*np.exp(1j*ph) + (G1+1j*G2)*np.exp(-1j*ph)
    f2 = np.exp(-2*(kx**2+ky**2)/scale**2) if cutoff=='circ' else np.exp(-2*(E/scale)**2)
    Esafe = np.where(E>0, E, np.inf)
    tot = np.zeros_like(kx)
    for s in (+1,-1):
        kyp = ky + s*q                      # bend shift q||y (kx, d1 unchanged)
        d2p = kyp; Ep = np.hypot(d1, d2p); php = np.arctan2(d2p, d1)
        z = -(G1-1j*G2)*np.exp(1j*ph) + (G1+1j*G2)*np.exp(-1j*php)
        tot += 0.25*((np.abs(z)/2)**2/(-E-Ep) - (np.abs(z0)/2)**2/(-2*Esafe))
    integ = np.sum(f2*tot)*dk*dk/(2*np.pi)**2
    return 2*integ/q**2

print("== VALIDATION: Cartesian circular cutoff vs polar reference ==")
for L, ref in [(10,0.1189),(20,0.3017),(40,0.7900)]:
    val = band_K_bend('circ', L)
    print(f"  Lam={L:3d}: Cartesian K_bend={val:.4f}   polar ref={ref:.4f}   ratio={val/ref:.3f}")

print("\n== ISO-ENERGY cutoff: band_K(m_xi), single-valued (no light/heavy bracket) ==")
mxis = np.array([20.,30.,40.,50.,60.,70.])
Ks = np.array([band_K_bend('ener', m) for m in mxis])
aE, lnAE = np.polyfit(np.log(mxis), np.log(Ks), 1); AE = np.exp(lnAE)
print(f"  band_K(m_xi) = {AE:.4e} * m_xi^{aE:.3f}   (energy cutoff)")
for m,K in zip(mxis,Ks): print(f"    m_xi={m:5.1f}  (m_xi/E*={m/Estar:5.2f})  band_K={K:.4f}")

print("\n== INVERT for the REQUIRED m_xi/E* (geom_coeff=1 convention) ==")
for target,label in [(1/(2*np.pi),'cH_0/2pi'), (0.176,'a0_obs/cH_0')]:
    m_req = (target/AE)**(1/aE)
    print(f"  band_K={target:.4f} ({label:12s}) -> m_xi={m_req:5.2f}  m_xi/E*={m_req/Estar:5.2f}"
          f"   vs phi^6={phi**6:.2f}  (ratio {(m_req/Estar)/phi**6:.2f})")

print("\n== forward: what a0 does m_xi/E*=phi^n give (energy cutoff, geom_coeff=1) ==")
cH0 = 2.99792458e8*2.27e-18
for n in (5,6,7):
    m = phi**n*Estar
    K = band_K_bend('ener', m)
    print(f"  phi^{n} (m_xi/E*={phi**n:5.2f}): band_K={K:.4f}  a0={K*cH0:.2e}  (obs 1.2e-10, ratio {K*cH0/1.2e-10:.2f})")
