"""
WP3: resolve the light-axis velocity v (c vs c/phi), then RECONCILE the WP4 gate
with the framework's actual m_xi (which is ~phi^10 H0, not H0).

WP3 argument for v: the DM medium is the condensate AT THE ATTRACTOR rho_inf=phi/beta.
Radial propagation there has the k-essence sound speed c_s = 1/phi (P7:rem:cs_de). The
v->c assumption is the O(beta^0) low-density limit, not the attractor. So the light-axis
band velocity for the DM medium is v = c/phi. (Subdominant: factor phi vs the m_xi factor.)

Unit map: K_phys = K_band * rho_cond * ell_dir^2, ell_dir = v/m_xi (director Compton length).
"""
import numpy as np
phi=(1+5**0.5)/2
c=2.99792458e8; G=6.674e-11; hbar=1.0546e-34; eV=1.602e-19; kB=1.381e-23
hbarc=hbar*c; H0=2.27e-18; kpc=3.0857e19; Mpc=1e3*kpc
Tcmb=2.7255*kB/eV; Lam_cond=Tcmb*(np.pi**2/150)**0.25
rho_cond=10*(Lam_cond*eV)**4/hbarc**3
Kband=0.16
vflat=200e3; K_RC=c**2*vflat**2/(4*np.pi*G)

def Kphys(v,mxi):
    ell=v/mxi; return Kband*rho_cond*ell**2
def report(vlabel,v,mxilabel,mxi):
    r=Kphys(v,mxi)/K_RC
    print(f"  v={vlabel:6s}  m_xi={mxilabel:12s}: K_phys/K_RC = {r:.2e}  ({np.log10(r):+.1f} dex)")

# H0 as an energy (for m_xi in the same units as v/ell): use m_xi as inverse-length via c/m_xi_len
# work directly with ell_dir = v / m_xi where m_xi has freq units (rad/s): ell=v/(m_xi_freq).
print("K_RC = %.2e N ;  rho_cond=%.2e J/m^3\n"%(K_RC,rho_cond))
print("WP4 gate, K_phys/K_RC, over v in {c, c/phi} and m_xi in {H0, phi^10 H0, phi^10 sqrt3 H0}:")
for vlab,v in [("c",c),("c/phi",c/phi)]:
    for mlab,m in [("H0",H0),("phi^10 H0",phi**10*H0),("phi^10 sqrt3 H0",phi**10*np.sqrt(3)*H0)]:
        report(vlab,v,mlab,m)

print("\nExact-match m_xi (K_phys=K_RC) for each v:")
for vlab,v in [("c",c),("c/phi",c/phi)]:
    ell_need=np.sqrt(K_RC/(Kband*rho_cond)); m_need=v/ell_need
    print(f"  v={vlab:6s}: m_xi_need = {m_need/H0:.2f} H0")

print("""
VERDICT (WP3 + reconciliation):
- WP3: v = c/phi (attractor k-essence sound speed) is the physical light-axis velocity for
  the DM medium; v=c was the low-density assumption. Effect on WP4: factor phi^2=2.6 only.
- The DOMINANT lever is m_xi. The framework FORMULA m_xi = phi^10 sqrt3 H0 (~213 H0) makes
  the gate UNDERSHOOT by ~3 dex; the operational "m_xi ~ H0" (used in the Vainshtein sector)
  gives ~1-2 dex. The exact match wants m_xi ~ few H0.
- So closure hinges on the phi^10 prefactor in m_xi, which the paper states INCONSISTENTLY
  (formula phi^10 sqrt3 H0 ~ 213 H0  vs  repeated "m_xi ~ H0 ~ 1e-33 eV"). This is a real
  internal tension to flag, NOT something WP3 resolves.
- Robust regardless: the elastic-DM scale is 50+ dex closer than the condensate-coherence
  (mm) length -- the ULTRALIGHT-director mechanism is qualitatively right (cosmological
  elastic length needed & supplied). Quantitatively it is ~1-3 dex off depending on m_xi,
  and the O(1..100) coefficient is NOT pinned. R1 stands.
""")
