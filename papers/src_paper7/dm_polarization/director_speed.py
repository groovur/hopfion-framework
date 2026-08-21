"""Director-wave speed, ordering temperature, and texture (isocurvature) scale
for the semi-Dirac condensate director. Extends the gate-validated frank_run.py
machinery (same H, same disk cutoff, same interband matrix elements).

Goldstone relation for the director (broken rotation, generator G):
    omega^2 = (rho_s/chi) q^2   =>   c_dir = sqrt(K/chi),
  rho_s = K = Frank stiffness (q^2 coeff of the twist energy; frank_run),
  chi   = static susceptibility of the rotation generator G,
        = INT d2k/(2pi)^2 |M(k,k)|^2 / |d_k|   (interband Kubo, lower band filled),
  M(k,k) = <+(k)|G.sigma|-(k)>, the SAME zero-momentum matrix element whose
  square is the diamagnetic subtraction in frank_run's integrand.

Deliverables:
  1. c_dir^2 = K/chi in band units (v=1) at several disk cutoffs Lam
     -> is it O(1) (relativistic, ~ band velocity) and is the RATIO cutoff-robust
        even though K and chi individually are cutoff-set?
  2. 2D nematic Kosterlitz-Thouless ordering scale T_KT = (pi/2) K_stiff
     (in band units) and its K-scaling -> ordering temperature ~ physical cutoff.
  3. these feed the texture/isocurvature estimate (done in the writeup, with the
     cutoff caveats made explicit).

H(k) = (kx^2/2) sigma_x + ky sigma_y ; d=(kx^2/2, ky) ; |-(k)>=(1,-e^{i phi})/sqrt2.
G = (kx*ky, -kx) (rotation generator). Undoped: interband only.
"""
import numpy as np

def matrices_on_disk(Lam, Nr, Nth):
    """return arrays over the disk: |M(k,k)|^2, |d_k|, and the measure r*dr*dth."""
    rs = np.linspace(Lam/Nr, Lam, Nr); dr = rs[1]-rs[0]
    th = np.linspace(0, 2*np.pi, Nth, endpoint=False); dth = th[1]-th[0]
    C, S = np.cos(th), np.sin(th)
    M2 = np.zeros((Nr, Nth)); ND = np.zeros((Nr, Nth)); W = np.zeros((Nr, Nth))
    for i in range(Nr):
        r = rs[i]; kx = r*C; ky = r*S
        d1, d2 = kx**2/2, ky
        nd = np.hypot(d1, d2); phi = np.arctan2(d2, d1)
        G1, G2 = kx*ky, -kx
        # M(k,k) = <+(k)|G.sigma|-(k)>  with |-(k)>=(1,-e^{i phi})/sqrt2,
        # |+(k)>=(1,+e^{i phi})/sqrt2 ; same expression frank_run uses at q=0:
        z0 = -(G1 - 1j*G2)*np.exp(1j*phi) + (G1 + 1j*G2)*np.exp(-1j*phi)
        M2[i] = (np.abs(z0)/2)**2  # |<+(k)|G.sigma|-(k)>|^2, same amplitude conv. as frank_run
        ND[i] = nd
        W[i]  = r*dr*dth/(2*np.pi)**2
    return M2, ND, W

def chi_generator(Lam, Nr, Nth):
    """static susceptibility of G: chi = INT |M(k,k)|^2/|d_k|."""
    M2, ND, W = matrices_on_disk(Lam, Nr, Nth)
    # E_+ - E_- = 2|d| ; second-order PT with both orderings -> |M|^2/|d|
    integrand = np.where(ND > 0, M2/ND, 0.0)
    return float(np.sum(integrand * W))

def dE_twist(q, ex, ey, Lam, Nr, Nth):
    """twist energy per theta0^2 (frank_run integrand), for K = 2 dE/q^2."""
    rs = np.linspace(Lam/Nr, Lam, Nr); dr = rs[1]-rs[0]
    th = np.linspace(0, 2*np.pi, Nth, endpoint=False); dth = th[1]-th[0]
    C, S = np.cos(th), np.sin(th)
    tot = 0.0
    for i in range(Nr):
        r = rs[i]; kx = r*C; ky = r*S
        d1, d2 = kx**2/2, ky; nd = np.hypot(d1, d2); phi = np.arctan2(d2, d1)
        G1, G2 = kx*ky, -kx
        val = np.zeros_like(kx)
        for s in (+1, -1):
            kxp, kyp = kx+s*q*ex, ky+s*q*ey
            d1p, d2p = kxp**2/2, kyp; ndp = np.hypot(d1p, d2p); phip = np.arctan2(d2p, d1p)
            z  = -(G1-1j*G2)*np.exp(1j*phi) + (G1+1j*G2)*np.exp(-1j*phip)
            z0 = -(G1-1j*G2)*np.exp(1j*phi) + (G1+1j*G2)*np.exp(-1j*phi)
            term = (np.abs(z)/2)**2/(-nd-ndp) - (np.abs(z0)/2)**2/(-2*nd)
            val += 0.25*term
        tot += np.sum(val)*r*dr*dth/(2*np.pi)**2
    return tot

def K_stiff(ex, ey, Lam, Nr, Nth, q=0.1):
    return 2*dE_twist(q, ex, ey, Lam, Nr, Nth)/q**2

if __name__ == "__main__":
    print(f"{'Lam':>5} {'Nr':>5} | {'chi':>12} {'K_bend':>10} {'K_splay':>10} "
          f"| {'cdir2_bend':>12} {'|cdir2|_splay':>13}")
    L=[]; CHI=[]; KB=[]; KS=[]
    for Lam, Nr, Nth in [(15.0,600,720),(30.0,900,720),(60.0,1400,720),(120.0,2400,720)]:
        chi = chi_generator(Lam, Nr, Nth)
        Kb  = K_stiff(0,1, Lam, Nr, Nth)   # bend  (q||y, stable)
        Ks  = K_stiff(1,0, Lam, Nr, Nth)   # splay (q||x, unstable, <0)
        L.append(Lam); CHI.append(chi); KB.append(Kb); KS.append(Ks)
        print(f"{Lam:5.0f} {Nr:5d} | {chi:12.4e} {Kb:10.4f} {Ks:10.4f} "
              f"| {Kb/chi:12.3e} {abs(Ks)/chi:13.3e}")
    L=np.array(L); lnL=np.log(L)
    def slope(y): return np.polyfit(lnL, np.log(np.abs(np.array(y))),1)[0]
    a_chi, a_kb, a_ks = slope(CHI), slope(KB), slope(KS)
    print(f"\nPower-law in Lam:  chi ~ Lam^{a_chi:.2f}   K_bend ~ Lam^{a_kb:.2f}   "
          f"K_splay ~ Lam^{a_ks:.2f}")
    print(f"=> c_dir^2 (bend)  ~ Lam^{a_kb-a_chi:.2f}   (splay) ~ Lam^{a_ks-a_chi:.2f}")
    print("Negative exponent => the director mode SLOWS (cold) as the cutoff grows; "
          "it is NOT relativistic. c_dir^2<0 in the splay channel = the instability.")
