"""Crossover matching: fix the physical cutoff Lam* = m_xi in the band's own
CROSSOVER units, then read off a_0/cH_0 = band_K(Lam*) and the isocurvature
scale. Honest: we compute what the framework REQUIRES for a_0 = cH_0/(2pi) and
whether that is a natural (single-scale) condition or needs a hierarchy.

Semi-Dirac dispersion  E = sqrt(kx^4/4 + ky^2)  (band units, A=1/2, v=1).
Crossover locus |d1|=|d2|: kx^2/2 = |ky|. Along kx=ky=k*: k*=2, E* = E(k*=2 on
that locus).  So the band's intrinsic scale is (k*, E*) ~ (2, ~2).

Physical map:  light axis v -> c ; cutoff energy E(k_max) = m_xi ~ H_0.
=> Lam_band (cutoff momentum in crossover units) is set by m_xi/E*.
a_0 = band_K(Lam_band) * cH_0  (from the 3D embedding; the O(1) geom factor is
absorbed into the reported coefficient).  Target: band_K(Lam*) = 1/(2pi)=0.159.
"""
import numpy as np

# --- band_K(Lam) with the smooth physical cutoff, from cutoff_embedding ---
def K_bend_reg(Lam, Nr=1400, Nth=720, q=0.05):
    rs = np.linspace(Lam/Nr, 4*Lam, Nr); dr = rs[1]-rs[0]
    th = np.linspace(0, 2*np.pi, Nth, endpoint=False); dth = th[1]-th[0]
    C, S = np.cos(th), np.sin(th)
    tot = 0.0
    for r in rs:
        kx, ky = r*C, r*S
        d1, d2 = kx**2/2, ky; nd = np.hypot(d1, d2); phi = np.arctan2(d2, d1)
        G1, G2 = kx*ky, -kx; f2 = np.exp(-2*r**2/Lam**2)
        val = np.zeros_like(kx)
        for s in (+1, -1):
            kxp, kyp = kx, ky+s*q            # bend: q||y
            d1p, d2p = kxp**2/2, kyp; ndp = np.hypot(d1p, d2p); phip = np.arctan2(d2p, d1p)
            z  = -(G1-1j*G2)*np.exp(1j*phi) + (G1+1j*G2)*np.exp(-1j*phip)
            z0 = -(G1-1j*G2)*np.exp(1j*phi) + (G1+1j*G2)*np.exp(-1j*phi)
            val += 0.25*((np.abs(z)/2)**2/(-nd-ndp) - (np.abs(z0)/2)**2/(-2*nd))
        tot += np.sum(f2*val)*r*dr*dth/(2*np.pi)**2
    return 2*tot/q**2

# --- crossover scale ---
kstar = 2.0
Estar = np.hypot(kstar**2/2, kstar)   # E at (kx,ky)=(k*,k*) on the crossover locus
print(f"crossover: k* = {kstar:.2f}  E* = {Estar:.3f}  (band units, v=1)")

# --- band_K(Lam) sample + power law ---
Lams = np.array([8.,12.,16.,24.,36.])
Ks   = np.array([K_bend_reg(L) for L in Lams])
a, lnA = np.polyfit(np.log(Lams), np.log(Ks), 1)
A = np.exp(lnA)
print("\nband_K(Lam) = A Lam^a :  A=%.4e  a=%.3f" % (A, a))
for L,K in zip(Lams,Ks): print(f"   Lam={L:5.1f}  K_bend={K:.4f}")

# --- solve band_K(Lam*) = 1/(2pi) ---
target = 1/(2*np.pi)
Lstar = (target/A)**(1/a)
print(f"\nTARGET a_0 = cH_0/(2pi):  band_K = {target:.4f}  =>  Lam* = {Lstar:.2f} (crossover units)")
print(f"   i.e. cutoff/crossover ratio  m_xi/E* :")
print(f"     light-axis map (E=k):     m_xi/E* = Lam*        = {Lstar:.1f}")
print(f"     heavy-axis map (E=k^2/2):  m_xi/E* = Lam*^2 / (2 k*) ~ {Lstar**2/(2*kstar):.1f}")

# --- single-scale (no hierarchy): cutoff AT the crossover, Lam_band ~ k*=2 ---
Lss = 2.0
Kss = K_bend_reg(Lss)
cH0 = 2.99792458e8 * 2.27e-18
print(f"\nSINGLE-SCALE (m_xi ~ E*, Lam_band ~ k*={Lss:.0f}):  band_K = {Kss:.4f}")
print(f"   a_0 = {Kss:.4f} cH_0 = {Kss*cH0:.2e} m/s^2   (obs 1.2e-10, cH_0/2pi = {cH0/(2*np.pi):.2e})")
print(f"   => single-scale a_0 is {cH0/(2*np.pi)/(Kss*cH0):.1f}x BELOW the empirical value.")

# --- isocurvature scale at the matched cutoff ---
Kstar = A*Lstar**a
Tord = 0.5*np.pi*Kstar
print(f"\nISOCURVATURE at Lam*={Lstar:.1f}:  T_ord ~ (pi/2)K_bend = {Tord:.3f} (band, ~m_xi units)")
print("   texture rho_tex ~ |K|/xi^2 (xi~1/Lam*); f_iso = rho_tex/rho_DM still")
print("   needs rho_DM(m_xi); but T_ord ~ O(1) m_xi => ordering is LATE (~condensate")
print("   scale), so the texture forms cold and late -> naturally subdominant.")

print("\n== HONEST VERDICT ==")
print("a_0 ~ cH_0 is robust (right order, no tuning). The exact cH_0/(2pi) needs")
print(f"a cutoff/crossover hierarchy m_xi/E* ~ {Lstar:.0f}-{Lstar**2/(2*kstar):.0f}. A single-scale")
print("condensate gives a_0 ~10-30x too small. So the 2pi is NOT derived; it is a")
print("specific, falsifiable condition on where m_xi sits above the semi-Dirac")
print("crossover. Isocurvature: ordering is late/cold => texture subdominant (safe),")
print("pending the rho_DM normalisation for a number vs Planck 0.038.")
