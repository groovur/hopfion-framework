#!/usr/bin/env python3.11
"""
tube_trefoil.py -- PHASE 3 of the tube-adapted Q_H=3 saddle build
(notes/qh3_saddle_infrastructure_scope.md). Puts the validated curvilinear
Faddeev engine (tube_curvilinear.py, Phases 1-2 passed) on the actual TREFOIL,
with the 3 self-crossings handled by "blend first" (user, 2026-09-27): the
existing inverse-square two-strand spinor superposition, extended to a
PARTITION OF UNITY so the self-overlapping tube chart is not double-counted:
    W_s(X) = rho_s^-2 / sum_i rho_i^-2    (sum over strands covering X)
Away from crossings only one strand covers -> W=1; in the overlap the two
charts' W sum to 1. The energy DENSITY is chart-invariant (orthonormal frame),
so only the volume measure dV is reweighted.

VALIDATION: total energy cross-checked against the Cartesian FD engine on the
SAME blended field (Cartesian integrates lab space once -> no double-count, it
is ground truth). Also reports the fraction of K and J4 inside the crossing
regions -> decides whether Cartesian patches are worth building (user's plan).

Run:  python3.11 tube_trefoil.py
"""
import numpy as np, time
from scipy.spatial import KDTree
from bishop_frame_v2 import build_compensated_frame_arclength
from tube_curvilinear import curvilinear_energy, cartesian_energy, MU, PHI

R0, r0 = 3.0, 0.874
C_STAR = 2.5062


def _trefoil_pts(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                     (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)], -1)


def trefoil_frame(Ns, NT=30000):
    """Arc-length-uniform Bishop frame on the trefoil + Bishop curvatures k1,k2."""
    t_arr, T_t, N1_t, N2_t, H = build_compensated_frame_arclength(NT=NT)
    G_t = _trefoil_pts(t_arr)
    seg = np.linalg.norm(np.diff(np.vstack([G_t, G_t[:1]]), axis=0), axis=1)
    cumArc = np.concatenate([[0.0], np.cumsum(seg)])         # len NT+1
    L = cumArc[-1]
    s_grid = np.arange(Ns) * (L/Ns)
    t_of_s = np.interp(s_grid, cumArc[:NT], t_arr)           # arc-length -> t
    def vinterp(V):                                          # per-component interp in t
        return np.stack([np.interp(t_of_s, t_arr, V[:, k]) for k in range(3)], 1)
    Gamma = _trefoil_pts(t_of_s)
    T  = vinterp(T_t);  N1 = vinterp(N1_t);  N2 = vinterp(N2_t)
    # re-orthonormalise after interpolation
    T  /= np.linalg.norm(T, axis=1, keepdims=True)
    N1 -= (N1*T).sum(1, keepdims=True)*T
    N1 /= np.linalg.norm(N1, axis=1, keepdims=True)
    N2  = np.cross(T, N1)
    ds = L/Ns
    dTds = (np.roll(T, -1, 0) - np.roll(T, 1, 0)) / (2*ds)   # periodic central diff
    kappa1 = (dTds*N1).sum(1)
    kappa2 = (dTds*N2).sum(1)
    return s_grid, Gamma, T, N1, N2, kappa1, kappa2, L


# ----------------------------------------------------------------------------
# Blended construction-C field at arbitrary lab points (used by BOTH engines).
# Nearest two DISTINCT strands (arc-separated), inverse-square spinor blend.
# ----------------------------------------------------------------------------
class BlendField:
    def __init__(self, Cst=C_STAR, NT_curve=6000, NT_frame=30000, min_sep_frac=1/6):
        self.Cst = Cst
        self.t = np.linspace(0, 2*np.pi, NT_curve, endpoint=False)
        self.G = _trefoil_pts(self.t)
        self.tree = KDTree(self.G)
        tf, _, N1f, N2f, _ = build_compensated_frame_arclength(NT=NT_frame)
        self.tf, self.N1f, self.N2f, self.NTf = tf, N1f, N2f, NT_frame
        self.min_sep = 2*np.pi*min_sep_frac          # arc-separation (in t) for "distinct"

    def _frame(self, tq):
        idx = np.searchsorted(self.tf, tq % (2*np.pi)) % self.NTf
        return self.N1f[idx], self.N2f[idx]

    def nearest_two(self, X, k=60):
        d, idx = self.tree.query(X, k=k, workers=-1)
        t1 = self.t[idx[:, 0]]; d1 = d[:, 0]
        dt = np.abs(((self.t[idx] - t1[:, None] + np.pi) % (2*np.pi)) - np.pi)
        far = dt > self.min_sep
        fi = np.argmax(far, axis=1)
        has2 = far.any(axis=1)
        rows = np.arange(len(X))
        d2 = np.where(has2, d[rows, fi], np.inf)
        t2 = self.t[idx[rows, fi]]
        return t1, d1, t2, d2

    def _phase(self, X, tq):
        N1q, N2q = self._frame(tq)
        rel = X - _trefoil_pts(tq)
        chi = np.arctan2((rel*N2q).sum(1), (rel*N1q).sum(1))
        return chi + 3*tq

    def n_at(self, X):
        Xf = X.reshape(-1, 3)
        t1, d1, t2, d2 = self.nearest_two(Xf)
        r1 = np.clip(d1, 1e-6, None); r2 = np.clip(d2, 1e-6, None)
        P1 = self._phase(Xf, t1); P2 = self._phase(Xf, t2)
        f0 = lambda r: 2*np.arctan((r*self.Cst)**(-self.Cst))
        f1 = f0(r1); f2 = f0(r2)
        w1 = 1/r1**2; w2 = np.where(np.isfinite(d2), 1/r2**2, 0.0)
        z1 = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
        z2 = w1*np.sin(f1/2)*np.exp(1j*P1) + w2*np.sin(f2/2)*np.exp(1j*P2)
        mag = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2).clip(1e-12)
        z1 /= mag; z2 /= mag
        n = np.stack([2*np.real(np.conj(z1)*z2),
                      2*np.imag(np.conj(z1)*z2),
                      np.abs(z1)**2 - np.abs(z2)**2], -1)
        return n.reshape(X.shape)


def partition_weight(bf, X_grid, s_of_node, rho_of_node, L, rho_max):
    """W = rho_home^-2 / (rho_home^-2 + rho_other^-2); rho_home = node's chart radius;
    rho_other = distance to nearest strand whose arc-param is FAR from the node's home s.
    W=1 where no other strand is within rho_max (single-chart coverage)."""
    Xf = X_grid.reshape(-1, 3)
    d, idx = bf.tree.query(Xf, k=60, workers=-1)
    tt = bf.t[idx]                                       # (M,k) params of neighbours
    # home arc param for each node: t at that node's s (map s->t via curve length)
    # cheaper: home strand = the neighbour nearest in arc to the node's own (s->t)
    t_home = np.interp(s_of_node.ravel(), np.linspace(0, L, len(bf.t), endpoint=False), bf.t)
    dt = np.abs(((tt - t_home[:, None] + np.pi) % (2*np.pi)) - np.pi)
    other = dt > bf.min_sep
    big = np.where(other, d, np.inf)
    rho_other = big.min(axis=1)                          # nearest OTHER-strand distance
    rho_home = np.clip(rho_of_node.ravel(), 1e-6, None)
    w_home = 1.0/rho_home**2
    w_other = np.where(rho_other < rho_max, 1.0/np.clip(rho_other, 1e-6, None)**2, 0.0)
    W = w_home/(w_home + w_other)
    return W.reshape(X_grid.shape[:-1])


def run_trefoil(Cst=C_STAR, rho_max=1.1):
    print(f"=== PHASE 3 trefoil (blend + partition of unity): C*={Cst}, rho_max={rho_max} ===")
    bf = BlendField(Cst=Cst)
    print("  [curvilinear tube chart + partition weight]")
    print(f"  {'(Ns,Nr,Np)':>16} {'K':>10} {'J4':>10} {'E':>12} {'Wmin':>6} {'olap%':>6} {'s':>5}")
    cur = None
    for (Ns, Nr, Np) in [(200,20,48),(300,30,64),(420,40,80)]:
        t0 = time.time()
        s_grid, Gamma, T, N1, N2, k1, k2, L = trefoil_frame(Ns)
        rho_grid = (np.arange(Nr)+0.5)*(rho_max/Nr)
        psi_grid = np.arange(Np)*(2*np.pi/Np)
        S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
        cosp, sinp = np.cos(PSI), np.sin(PSI)
        er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
        Xg = Gamma[:, None, None, :] + RHO[..., None]*er
        W = partition_weight(bf, Xg, S, RHO, L, rho_max)
        olap = 100*(W < 0.999).mean()
        cur = curvilinear_energy(Gamma, T, N1, N2, k1, k2, s_grid, rho_grid, psi_grid,
                                 lambda X, S, RHO, PSI: bf.n_at(X), beta=0.0, node_weight=W)
        print(f"  {str((Ns,Nr,Np)):>16} {cur['K']:>10.3f} {cur['J4']:>10.3f} {cur['E']:>12.1f}"
              f" {W.min():>6.3f} {olap:>6.2f} {time.time()-t0:>5.1f}")
    box = 2*(R0 + r0 + rho_max + 0.4)
    print(f"\n  [Cartesian FD cross-check on the SAME blended field, box={box:.1f}]")
    print(f"  {'N':>16} {'K':>10} {'J4':>10} {'E':>12}  {'dK%':>7} {'dJ4%':>7}")
    for N in [140, 200]:
        t0 = time.time()
        car = cartesian_energy(bf.n_at, box, N, beta=0.0)
        dK = 100*(car['K']-cur['K'])/cur['K']; dJ4 = 100*(car['J4']-cur['J4'])/cur['J4']
        print(f"  {N:>16} {car['K']:>10.3f} {car['J4']:>10.3f} {car['E']:>12.1f}"
              f"  {dK:>7.2f} {dJ4:>7.2f}  ({time.time()-t0:.1f}s)")
    print("\n  Gate: curvilinear+partition K,J4 should match the Cartesian engine (which"
          "\n  counts lab space once). Agreement validates the blend-crossing treatment.")


def crossing_energy_fraction(Cst=C_STAR, rho_max=1.1, Ns=420, Nr=40, Np=80):
    """How much of K and J4 sits in the overlap (crossing) regions -> patch decision.
    Recomputes the curvilinear densities inline and splits by W<0.999."""
    bf = BlendField(Cst=Cst)
    s_grid, Gamma, T, N1, N2, k1, k2, L = trefoil_frame(Ns)
    rho_grid = (np.arange(Nr)+0.5)*(rho_max/Nr); psi_grid = np.arange(Np)*(2*np.pi/Np)
    ds = s_grid[1]-s_grid[0]; drho = rho_grid[1]-rho_grid[0]; dpsi = psi_grid[1]-psi_grid[0]
    S, RHO, PSI = np.meshgrid(s_grid, rho_grid, psi_grid, indexing='ij')
    cosp, sinp = np.cos(PSI), np.sin(PSI)
    er = cosp[..., None]*N1[:, None, None, :] + sinp[..., None]*N2[:, None, None, :]
    Xg = Gamma[:, None, None, :] + RHO[..., None]*er
    n = bf.n_at(Xg); n /= np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-12)
    W = partition_weight(bf, Xg, S, RHO, L, rho_max)
    hs = 1.0 - RHO*(k1[:, None, None]*cosp + k2[:, None, None]*sinp)
    Dn_s = (np.roll(n,-1,0)-np.roll(n,1,0))/(2*ds)/hs[...,None]
    Dn_p = (np.roll(n,-1,2)-np.roll(n,1,2))/(2*dpsi)/RHO[...,None]
    Dn_r = np.empty_like(n)
    Dn_r[:,1:-1]=(n[:,2:]-n[:,:-2])/(2*drho); Dn_r[:,0]=(n[:,1]-n[:,0])/drho; Dn_r[:,-1]=(n[:,-1]-n[:,-2])/drho
    g2=(Dn_s**2).sum(-1)+(Dn_r**2).sum(-1)+(Dn_p**2).sum(-1)
    s4=np.clip(1-n[...,2]**2,0,1)**2
    Fsr=(n*np.cross(Dn_s,Dn_r)).sum(-1); Fsp=(n*np.cross(Dn_s,Dn_p)).sum(-1); Frp=(n*np.cross(Dn_r,Dn_p)).sum(-1)
    j4d=Fsr**2+Fsp**2+Frp**2
    dV=hs*RHO*ds*drho*dpsi*W
    Kd=s4*g2+MU*g2
    cross=W<0.999
    Ktot=(Kd*dV).sum(); J4tot=(j4d*dV).sum()
    Kx=(Kd*dV)[cross].sum(); J4x=(j4d*dV)[cross].sum()
    print(f"=== crossing-energy fraction (Ns,Nr,Np={Ns},{Nr},{Np}, rho_max={rho_max}) ===")
    print(f"  overlap nodes: {100*cross.mean():.2f}%   Wmin={W.min():.3f}")
    print(f"  K   in crossings: {100*Kx/Ktot:6.2f}%   (K_tot={Ktot:.1f})")
    print(f"  J4  in crossings: {100*J4x/J4tot:6.2f}%   (J4_tot={J4tot:.1f})")
    print("  -> small fraction => blend is sufficient, Cartesian patches NOT worth it;")
    print("     large fraction => J4 in crossings is under-resolved by the blend, patch.")


if __name__ == "__main__":
    import sys
    if "--diag" in sys.argv:
        crossing_energy_fraction()
    else:
        run_trefoil()
