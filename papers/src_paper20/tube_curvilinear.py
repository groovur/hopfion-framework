#!/usr/bin/env python3.11
"""
tube_curvilinear.py -- Faddeev energy of a density-feedback director field in
TUBE-ADAPTED (Bishop-frame) curvilinear coordinates (s, rho, psi).

Motivation (notes/qh3_saddle_infrastructure_scope.md, option 3a): the blocked
Paper XX O1-O3 need the converged jammed Q_H=3 saddle, but a uniform Cartesian
grid cannot resolve the thin trefoil tube (radius ~1/C*~0.4) without N~320 and
even then J4 is only marginally converged. Concentrating DOF inside the tube via
(s,rho,psi) coords cuts ~3e7 -> ~1e5 DOF. This file is the FOUNDATION: the
curvilinear energy functional and its validation gates. Crossings + relaxation
come in later phases.

WHY BISHOP FRAME => the metric is DIAGONAL. Position
    x(s,rho,psi) = Gamma(s) + rho (cos psi N1(s) + sin psi N2(s)).
Parallel transport gives N1' = -kappa1 T, N2' = -kappa2 T (NO N1<->N2 twist term),
so ds^2 metric is diag(h_s^2, 1, rho^2) with
    h_s = 1 - rho (kappa1 cos psi + kappa2 sin psi),   sqrt(g) = h_s * rho.
Orthonormal directional derivatives of the LAB-frame director n=(nx,ny,nz):
    D_s n = (1/h_s) d_s n,   D_rho n = d_rho n,   D_psi n = (1/rho) d_psi n.
Then, EXACTLY as in the Cartesian solver (qh3_feedback_saddle.py:E_fb) but with
these scaled derivatives and the h_s*rho volume weight:
    g2   = |D_s n|^2 + |D_rho n|^2 + |D_psi n|^2
    s4   = (1 - nz^2)^2
    J2a  = INT s4 g2 dV ;  J2iso_fb = INT g2/(1+beta g2) dV ;  K = J2a + MU J2iso_fb
    F_sr = n.(D_s n x D_rho n), F_sp = n.(D_s n x D_psi n), F_rp = n.(D_rho n x D_psi n)
    J4   = INT (F_sr^2 + F_sp^2 + F_rp^2) dV ;   E = K J4
n is kept in the FIXED lab frame so s4 and the Hopf structure are identical to
the Cartesian code (no rotation of the target S^2).

Run:  python3.11 tube_curvilinear.py    # runs the Phase-1 straight-tube gates
"""
import numpy as np

PHI = (1 + 5**0.5) / 2
MU  = 3.0 - PHI          # isotropic-stiffness coefficient, matches qh3_feedback_saddle


# ----------------------------------------------------------------------------
# Core: Faddeev energy on a (s,rho,psi) grid given a Bishop frame along Gamma(s).
# Convention: s and psi are PERIODIC; rho is a bounded radial coord, STAGGERED
# (rho_j = (j+0.5) drho) so no node sits on the singular axis rho=0.
# ----------------------------------------------------------------------------
def curvilinear_energy(Gamma, T, N1, N2, kappa1, kappa2,
                       s_grid, rho_grid, psi_grid, build_n_lab, beta=0.0,
                       node_weight=None):
    """node_weight: optional (Ns,Nr,Np) partition-of-unity weight multiplying the
    volume measure, for self-overlapping charts (crossings). None => all-ones (an
    injective chart, e.g. straight/circular). The energy DENSITY is chart-invariant
    (orthonormal frame); only dV is reweighted so overlap volume is counted once."""
    Ns, Nr, Np = len(s_grid), len(rho_grid), len(psi_grid)
    ds   = s_grid[1]   - s_grid[0]
    drho = rho_grid[1] - rho_grid[0]
    dpsi = psi_grid[1] - psi_grid[0]     # = 2 pi / Np

    S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
    cosp, sinp = np.cos(PSI), np.sin(PSI)

    # lab positions of every grid node
    er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
    X  = Gamma[:, None, None, :] + RHO[..., None]*er
    n  = build_n_lab(X, S, RHO, PSI)                    # (Ns,Nr,Np,3), unit lab vectors
    n  = n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-12)

    # metric scale factor h_s (>0 required: rho < 1/kappa everywhere)
    hs = 1.0 - RHO*(kappa1[:, None, None]*cosp + kappa2[:, None, None]*sinp)

    # orthonormal directional derivatives (central; periodic in s, psi)
    Dn_s   = (np.roll(n, -1, 0) - np.roll(n, 1, 0)) / (2*ds)   / hs[..., None]
    Dn_psi = (np.roll(n, -1, 2) - np.roll(n, 1, 2)) / (2*dpsi) / RHO[..., None]
    Dn_rho = np.empty_like(n)                                  # rho not periodic
    Dn_rho[:, 1:-1] = (n[:, 2:] - n[:, :-2]) / (2*drho)
    Dn_rho[:, 0]    = (n[:, 1]  - n[:, 0])   / drho            # one-sided at inner edge
    Dn_rho[:, -1]   = (n[:, -1] - n[:, -2])  / drho            # one-sided at outer edge

    g2 = (Dn_s**2).sum(-1) + (Dn_rho**2).sum(-1) + (Dn_psi**2).sum(-1)
    nz = n[..., 2]
    s4 = np.clip(1 - nz**2, 0, 1)**2

    dV       = hs * RHO * ds * drho * dpsi
    if node_weight is not None:
        dV = dV * node_weight
    J2a      = (s4 * g2 * dV).sum()
    J2iso_fb = (g2 / (1.0 + beta*g2) * dV).sum()
    K        = J2a + MU*J2iso_fb

    Fsr = (n * np.cross(Dn_s,   Dn_rho)).sum(-1)
    Fsp = (n * np.cross(Dn_s,   Dn_psi)).sum(-1)
    Frp = (n * np.cross(Dn_rho, Dn_psi)).sum(-1)
    J4  = ((Fsr**2 + Fsp**2 + Frp**2) * dV).sum()

    return dict(K=float(K), J4=float(J4), J2a=float(J2a),
                J2iso_fb=float(J2iso_fb), E=float(K*J4))


# ----------------------------------------------------------------------------
# Straight-tube frame: Gamma(s)=(0,0,s), constant orthonormal frame, kappa=0.
# ----------------------------------------------------------------------------
def straight_frame(L, Ns):
    s = np.arange(Ns) * (L / Ns)                     # periodic sampling of [0,L)
    Gamma = np.stack([np.zeros(Ns), np.zeros(Ns), s], 1)
    T  = np.tile(np.array([0, 0, 1.0]), (Ns, 1))
    N1 = np.tile(np.array([1.0, 0, 0]), (Ns, 1))
    N2 = np.tile(np.array([0, 1.0, 0]), (Ns, 1))
    kappa1 = np.zeros(Ns); kappa2 = np.zeros(Ns)
    return s, Gamma, T, N1, N2, kappa1, kappa2


# ----------------------------------------------------------------------------
# Analytic straight-tube test field and its INDEPENDENT continuum energy.
# n = (sin f cosPhi, sin f sinPhi, cos f),  f=f(rho),  Phi = psi + q k s  (k=2pi/L).
# For a straight tube (h_s=1) the integrand is s-independent in magnitude, so the
# continuum energy = L * (2D rho,psi integral), computed with analytic derivatives:
#   D_rho n : |.|^2 = f'^2
#   D_psi n : |.|^2 = sin^2 f / rho^2
#   D_s   n : |.|^2 = (q k)^2 sin^2 f
#   Hopf comps:  F_rp = f' sin f / rho ; F_sr = q k f' sin f ; F_sp = 0
# NOTE F_sp = n.(D_s n x D_psi n) = 0 here: Phi = psi + q k s enters n only through
# Phi, so d_s n = q k d_Phi n and d_psi n = d_Phi n are PARALLEL -> cross product 0.
# (An earlier hand-truth wrongly set F_sp = q k sin^2 f/rho; the FD engine returned
# ~0 and caught it -- exactly the §1 discipline. Only F_sr, F_rp survive.)
# ----------------------------------------------------------------------------
def make_field(Cst, q, L):
    k = 2*np.pi / L
    def f_of(rho):      return 2*np.arctan((rho*Cst)**(-Cst))
    def fp_of(rho):     # d f / d rho
        u = (rho*Cst)**(-Cst)
        return 2 * (1/(1+u**2)) * (-Cst) * (rho*Cst)**(-Cst-1) * Cst
    def build_n_lab(X, S, RHO, PSI):
        f   = f_of(RHO); Phi = PSI + q*k*S
        return np.stack([np.sin(f)*np.cos(Phi),
                         np.sin(f)*np.sin(Phi),
                         np.cos(f)], -1)
    return f_of, fp_of, k, build_n_lab


def continuum_truth(Cst, q, L, rho_max, Nr_fine=200000):
    """High-res 1D (rho) midpoint quadrature of the analytic straight-tube energy
    densities, x (2 pi) in psi x L in s. Independent of the FD engine."""
    f_of, fp_of, k, _ = make_field(Cst, q, L)
    rho = (np.arange(Nr_fine)+0.5) * (rho_max/Nr_fine)
    drho = rho_max/Nr_fine
    f  = f_of(rho); fp = fp_of(rho); sf = np.sin(f); nz = np.cos(f)
    g2 = fp**2 + sf**2/rho**2 + (q*k)**2 * sf**2
    s4 = (1-nz**2)**2
    # dV per unit (psi,s): rho drho ; then x 2pi (psi) x L (s)
    w = rho*drho * (2*np.pi) * L
    J2a      = (s4*g2*w).sum()
    J2iso    = (g2*w).sum()                       # beta=0 truth
    K        = J2a + MU*J2iso
    Frp = fp*sf/rho
    Fsr = q*k*fp*sf
    # F_sp = 0 : d_s n and d_psi n are parallel (both ~ d_Phi n) -> cross product 0
    J4  = ((Frp**2+Fsr**2)*w).sum()
    return dict(K=K, J4=J4, J2a=J2a, J2iso_fb=J2iso, E=K*J4)


def _gate(name, Cst, q, L, rho_max, resolutions):
    print(f"\n=== {name}: C*={Cst}, q={q}, L={L}, rho_max={rho_max} ===")
    truth = continuum_truth(Cst, q, L, rho_max)
    print(f"  continuum truth:  K={truth['K']:.5f}  J4={truth['J4']:.5f}  E={truth['E']:.4f}")
    _, _, _, build_n_lab = make_field(Cst, q, L)
    print(f"  {'(Ns,Nr,Np)':>16} {'K':>10} {'J4':>10} {'E':>12}  {'errK%':>7} {'errJ4%':>7}")
    prev = None
    for (Ns, Nr, Np) in resolutions:
        s_grid, Gamma, T, N1, N2, k1, k2 = straight_frame(L, Ns)
        rho_grid = (np.arange(Nr)+0.5) * (rho_max/Nr)
        psi_grid = np.arange(Np) * (2*np.pi/Np)
        r = curvilinear_energy(Gamma, T, N1, N2, k1, k2,
                               s_grid, rho_grid, psi_grid, build_n_lab, beta=0.0)
        eK  = 100*(r['K']-truth['K'])/truth['K']
        eJ4 = 100*(r['J4']-truth['J4'])/truth['J4']
        print(f"  {str((Ns,Nr,Np)):>16} {r['K']:>10.4f} {r['J4']:>10.4f} {r['E']:>12.3f}"
              f"  {eK:>7.3f} {eJ4:>7.3f}")
        prev = r
    return truth, prev


# ============================================================================
# PHASE 2 -- circular tube (unknot). Tests the CURVATURE metric factor
# h_s = 1 - rho(kappa1 cos psi + kappa2 sin psi) = 1 + (rho/R) sin psi.
# Circle Gamma(s)=R(cos(s/R),sin(s/R),0). Valid Bishop frame (zero holonomy):
#   N1=(0,0,1) const [out of plane], N2=(cos,sin,0) [outward radial]
#   => kappa1=0, kappa2=-1/R  (verified: T'=-N2/R = kappa2 N2).
# Independent cross-check: the SAME analytic field on a Cartesian grid, via a
# numpy port of qh3_feedback_saddle.E_fb. Two unrelated discretizations must
# converge to one number -> validates h_s (§1).
# ============================================================================
def circular_frame(R, Ns):
    s = np.arange(Ns) * (2*np.pi*R / Ns)
    a = s / R
    Gamma = np.stack([R*np.cos(a), R*np.sin(a), np.zeros(Ns)], 1)
    T  = np.stack([-np.sin(a), np.cos(a), np.zeros(Ns)], 1)
    N1 = np.tile(np.array([0, 0, 1.0]), (Ns, 1))
    N2 = np.stack([np.cos(a), np.sin(a), np.zeros(Ns)], 1)
    kappa1 = np.zeros(Ns); kappa2 = -np.ones(Ns)/R
    return s, Gamma, T, N1, N2, kappa1, kappa2


def circular_field_lab(Cst, q, R):
    """Analytic unknot field at arbitrary LAB points X (for the Cartesian engine),
    via the analytic nearest-point-on-circle map. Same n(f,Phi) as the curvilinear
    grid so the two engines integrate one continuous field."""
    def f_of(rho): return 2*np.arctan((rho*Cst)**(-Cst))
    def n_at(X):
        x, y, z = X[..., 0], X[..., 1], X[..., 2]
        a = np.arctan2(y, x)                       # nearest-point param on circle
        Gx, Gy = R*np.cos(a), R*np.sin(a)
        er_x, er_y = np.cos(a), np.sin(a)          # N2 (radial), N1=(0,0,1)
        dx, dy, dz = x-Gx, y-Gy, z-0.0
        comp2 = dx*er_x + dy*er_y                  # along N2
        comp1 = dz                                 # along N1
        rho = np.sqrt(comp1**2 + comp2**2)
        # engine node sits at rho(cos psi N1 + sin psi N2): comp1=rho cos psi,
        # comp2=rho sin psi  =>  psi = arctan2(comp2, comp1). Must match exactly.
        psi = np.arctan2(comp2, comp1)
        f = f_of(np.clip(rho, 1e-9, None)); Phi = psi + q*a
        return np.stack([np.sin(f)*np.cos(Phi), np.sin(f)*np.sin(Phi), np.cos(f)], -1)
    return n_at


def circular_field_grid(Cst, q, R):
    """Same field addressed by (s,rho,psi) grid coords directly (curvilinear engine).
    psi convention must match circular_field_lab: comp along N1 = rho sin psi."""
    def f_of(rho): return 2*np.arctan((rho*Cst)**(-Cst))
    def build(X, S, RHO, PSI):
        f = f_of(np.clip(RHO, 1e-9, None)); Phi = PSI + q*(S/R)
        return np.stack([np.sin(f)*np.cos(Phi), np.sin(f)*np.sin(Phi), np.cos(f)], -1)
    return build


def cartesian_energy(n_at, box, N, beta=0.0):
    """Numpy port of qh3_feedback_saddle.E_fb (periodic central-diff FD, torus BC).
    Field must decay to vacuum before the box edge for the periodic BC to be clean."""
    h = box / N
    cv = h*(np.arange(N) - N//2 + 0.5)
    X = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), -1)
    n = n_at(X); n = n/np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-12)
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    def cd(u, a): return (np.roll(u, -1, a) - np.roll(u, 1, a)) / (2*h)
    nxx, nxy, nxz = cd(nx, 0), cd(nx, 1), cd(nx, 2)
    nyx, nyy, nyz = cd(ny, 0), cd(ny, 1), cd(ny, 2)
    nzx, nzy, nzz = cd(nz, 0), cd(nz, 1), cd(nz, 2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    s4 = np.clip(1-nz**2, 0, 1)**2
    dv = h**3
    J2a = (s4*g2).sum()*dv
    J2iso_fb = (g2/(1+beta*g2)).sum()*dv
    K = J2a + MU*J2iso_fb
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4 = (Fxy**2+Fxz**2+Fyz**2).sum()*dv
    return dict(K=float(K), J4=float(J4), J2a=float(J2a), J2iso_fb=float(J2iso_fb), E=float(K*J4))


def gate2(Cst=2.5062, q=1, R=3.0, rho_max=1.6):
    print(f"\n=== PHASE 2 circular tube: C*={Cst}, q={q}, R={R}, rho_max={rho_max} ===")
    print("  Two independent engines must converge to ONE energy (validates h_s curvature).")
    build = circular_field_grid(Cst, q, R)
    print(f"\n  [curvilinear, h_s=1+(rho/R)sin psi]")
    print(f"  {'(Ns,Nr,Np)':>16} {'K':>10} {'J4':>10} {'E':>12}")
    cur = None
    for (Ns, Nr, Np) in [(96,24,48),(160,40,80),(240,60,120),(320,80,160)]:
        s_grid, Gamma, T, N1, N2, k1, k2 = circular_frame(R, Ns)
        rho_grid = (np.arange(Nr)+0.5)*(rho_max/Nr)
        psi_grid = np.arange(Np)*(2*np.pi/Np)
        cur = curvilinear_energy(Gamma, T, N1, N2, k1, k2,
                                 s_grid, rho_grid, psi_grid, build, beta=0.0)
        print(f"  {str((Ns,Nr,Np)):>16} {cur['K']:>10.4f} {cur['J4']:>10.4f} {cur['E']:>12.2f}")
    n_at = circular_field_lab(Cst, q, R)
    box = 2*(R + rho_max + 0.4)
    print(f"\n  [Cartesian FD cross-check, box={box:.1f}]")
    print(f"  {'N':>16} {'K':>10} {'J4':>10} {'E':>12}  {'dK%':>7} {'dJ4%':>7}")
    for N in [120, 180, 240]:
        car = cartesian_energy(n_at, box, N, beta=0.0)
        dK  = 100*(car['K']-cur['K'])/cur['K']
        dJ4 = 100*(car['J4']-cur['J4'])/cur['J4']
        print(f"  {N:>16} {car['K']:>10.4f} {car['J4']:>10.4f} {car['E']:>12.2f}"
              f"  {dK:>7.3f} {dJ4:>7.3f}")
    print("\n  Gate passes if the two engines' finest K,J4 agree (dK%,dJ4% small &"
          "\n  shrinking). Cartesian is the harder-resolved one (thin tube on a coarse"
          "\n  Cartesian grid), so expect it to approach the curvilinear value from below.")


if __name__ == "__main__":
    import sys
    if "--phase2" in sys.argv:
        gate2()
        sys.exit(0)
    print("PHASE 1 -- straight-tube validation gates (h_s=1). Curvilinear FD energy"
          "\nmust converge to the independent analytic-quadrature continuum truth.")
    # 1a: purely poloidal (no s-winding) -> tests rho,psi operators + Jacobian
    _gate("1a poloidal-only", Cst=2.5062, q=0, L=6.0, rho_max=2.0,
          resolutions=[(4,16,24),(4,32,48),(4,64,96),(4,128,192)])
    # 1b: add meridional winding q=3 -> also tests the d_s operator (F_sr,F_sp,J4)
    _gate("1b winding q=3", Cst=2.5062, q=3, L=6.0, rho_max=2.0,
          resolutions=[(24,16,24),(48,32,48),(96,64,96),(160,96,144)])
    print("\nGate passes if errK%, errJ4% -> 0 as resolution increases.")
