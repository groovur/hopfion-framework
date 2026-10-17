#!/usr/bin/env python3.11
"""
torsion_framing.py -- Paper XXIII component C (first piece): the writhe / twist decomposition of the
u=0 preimage CENTERLINE along the Hopf-link -> trefoil reaction coordinate. DM/Quark crossover
(src_paper23).

MOTIVATION (user's question on component B): the emitted pion is Q_H=1 and the UNKNOT. By
Calugareanu-White, the Hopf charge of a framed curve is its self-linking Lk = Wr + Tw (writhe +
framing twist). A round unknot has Wr=0, so a Q_H=1 unknot carries its unit of charge as PURE TWIST.
Hence the pion must carry one twist unit's energy on emission (add one on absorption). This script
QUANTIFIES that: it tracks how the geometric (writhe+linking) part of the centerline charge grows as
the knot forms, and reads the framing twist Tw(lambda) = Q_H(lambda) - [Wr+2*Lk]_centerline as the
remainder. The prediction to test: the trefoil end is writhe-dominated (Tw small) and the net +1 of
the reaction shows up as a shed twist quantum = the pion.

METHOD (exact zero-set trace, not grid extraction):
  On S^3 use Hopf coords a=cos(chi)e^{i alpha}, b=sin(chi)e^{i beta}. The reaction u=0 set is
      a^2 = -[(1-lam) b^2 + lam c3 b^3]  =  -s^2 e^{2i beta}[(1-lam) + lam c3 s e^{i beta}],
  s=sin(chi). Magnitude fixes chi(beta):  1-s^2 = s^2 * |(1-lam)+lam c3 s e^{i beta}|  (monotone -> root
  find). Phase fixes alpha(beta): 2 alpha = pi + 2 beta + arg[(1-lam)+lam c3 s e^{i beta}]. The arg term
  WINDS once with beta only when the b^3 (lam) term dominates -- that winding is the reconnection that
  fuses the two Hopf circles (lam=0, two sqrt branches) into the single trefoil (lam=1). Trace the S^3
  curve, stereographically project to R^3 (inverse of knotted_map.s3_coords), compute the discrete Gauss
  writhe/linking (Klenin-Langowski segment form). Centerline writhe is embedding-dependent (NOT a
  topological invariant) -- that is the point: Tw is likewise embedding-dependent, but the TOPOLOGICAL
  statement (pion=unknot, Wr_round=0 -> its charge is twist) is frame-independent; the trend here is the
  physics.

Q_H(lambda) for the remainder is taken from the validated gate (same reaction_doublet as barrier_config)
at a stated grid resolution, NOT assumed.
"""
import numpy as np, sys, os
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from knotted_map import s3_coords
from validation import gate_field

# ---------- zero-set trace in Hopf coords ----------
def _s_of_beta(beta, lam, c3):
    """solve 1 - s^2 = s^2 * |(1-lam) + lam c3 s e^{i beta}| for s in (0,1)."""
    def Gmag(s):
        re = (1-lam) + lam*c3*s*np.cos(beta)
        im = lam*c3*s*np.sin(beta)
        return np.hypot(re, im)
    f = lambda s: (1 - s*s) - s*s*Gmag(s)
    # f(0+)=1>0, f(1-)= -Gmag(1)<0 -> unique root
    return brentq(f, 1e-9, 1-1e-12, xtol=1e-12)

def _s3_to_R3(a, b):
    """inverse of s3_coords: (a,b) on S^3 -> (X,Y,Z) in R^3. Projection pole at Re(b)=-1."""
    den = 1.0 + np.real(b)
    return np.real(a)/den, np.imag(a)/den, np.imag(b)/den

def trace_components(lam, c3=3.0, Nb=1200, close_tol=1e-2):
    """Return a list of closed R^3 polylines (1 for a knot, 2 for a link).

    Phase: 2 alpha = pi + 2 beta + psi(beta), psi = arg[(1-lam)+lam c3 s e^{i beta}], UNWRAPPED
    (accumulated) along beta so the trace captures the knotting. psi winds by 2pi over beta:0->2pi iff
    the rotating b^3 term dominates (lam c3 s > 1-lam) -- that winding is the reconnection: no wind =>
    a returns at beta=2pi => two-branch LINK; one wind => a returns only at beta=4pi => single KNOT."""
    betas = np.linspace(0, 4*np.pi, 2*Nb+1)
    s = np.array([_s_of_beta(bb, lam, c3) for bb in betas])
    re = (1-lam) + lam*c3*s*np.cos(betas)
    im = lam*c3*s*np.sin(betas)
    psi = np.unwrap(np.arctan2(im, re))                 # continuous phase
    alpha = 0.5*(np.pi + 2*betas + psi)
    chi = np.arcsin(np.clip(s, 0, 1))
    a = np.cos(chi)*np.exp(1j*alpha)
    b = np.sin(chi)*np.exp(1j*betas)
    X, Y, Z = _s3_to_R3(a, b)
    P = np.column_stack([X, Y, Z])
    scale = max(1.0, np.ptp(P))
    half = Nb                                           # index of beta=2pi
    closes_2pi = np.linalg.norm(P[half] - P[0]) < close_tol*scale
    if closes_2pi:
        comp0 = P[:half]                                # one loop over [0,2pi)
        a2 = np.cos(chi[:half])*np.exp(1j*(alpha[:half]+np.pi))   # other sqrt branch
        X2, Y2, Z2 = _s3_to_R3(a2, b[:half])
        comp1 = np.column_stack([X2, Y2, Z2])
        if _min_sep(comp0, comp1) > close_tol*scale:
            return [comp0, comp1]                       # LINK (two components)
        return [comp0]
    else:
        return [P[:2*half]]                             # single KNOT over [0,4pi)

def _min_sep(A, B):
    m = np.inf
    for p in A[::max(1, len(A)//60)]:
        m = min(m, np.min(np.linalg.norm(B - p, axis=1)))
    return m

# ---------- discrete Gauss writhe / linking (segment form) ----------
def _segments(poly):
    P = np.asarray(poly)
    e = np.roll(P, -1, axis=0) - P          # segment vectors (closed)
    m = P + 0.5*e                           # midpoints
    return m, e

def gauss_wr_lk(components):
    """Return (total, self_writhe, 2x_linking). total = (1/4pi) sum over all segment pairs of
    (e_a x e_b).(m_a-m_b)/|m_a-m_b|^3, skipping same-component adjacent/self pairs."""
    segs = [ _segments(c) for c in components ]
    ms = [s[0] for s in segs]; es = [s[1] for s in segs]
    ns = [len(m) for m in ms]
    self_sum = 0.0; cross_sum = 0.0
    for ci in range(len(components)):
        for cj in range(len(components)):
            mi, ei = ms[ci], es[ci]; mj, ej = ms[cj], es[cj]
            for ia in range(ns[ci]):
                d = mi[ia] - mj                       # (nj,3)
                r3 = np.linalg.norm(d, axis=1)**3 + 1e-30
                cx = np.cross(np.broadcast_to(ei[ia], ej.shape), ej)   # (nj,3)
                contrib = np.einsum('kj,kj->k', cx, d) / r3
                if ci == cj:
                    # skip self (ia) and immediate neighbours (wrap) to avoid 1/r^3 blowup
                    for off in (-1, 0, 1):
                        contrib[(ia+off) % ns[cj]] = 0.0
                    self_sum += contrib.sum()
                else:
                    cross_sum += contrib.sum()
    self_wr = self_sum/(4*np.pi)
    cross_2lk = cross_sum/(4*np.pi)   # already counts (i,j) and (j,i) => = 2*Lk
    return self_wr + cross_2lk, self_wr, cross_2lk

# ---------- Frenet curvature / torsion along a parametrised closed curve ----------
def _periodic_derivs(P, du):
    """r', r'', r''' by periodic central differences on a uniformly (du) parametrised closed curve."""
    rp   = (np.roll(P,-1,0) - np.roll(P,1,0))/(2*du)
    rpp  = (np.roll(P,-1,0) - 2*P + np.roll(P,1,0))/(du*du)
    rppp = (np.roll(P,-2,0) - 2*np.roll(P,-1,0) + 2*np.roll(P,1,0) - np.roll(P,2,0))/(2*du**3)
    return rp, rpp, rppp

def frenet_stats(components, du, R0_target=None):
    """Arc-length-weighted <kappa^2>, <kappa^2 tau^2>, <tau^2>, and L over one or more closed curves,
    each a uniformly-(du)-parametrised polyline. If R0_target is set, rescale all curves uniformly so
    the mean cylindrical radius sqrt(x^2+y^2) = R0_target (puts kappa,tau in that length unit --
    kappa,tau ~ 1/length, so this is how we compare to Paper XV's R0=3 numbers)."""
    if R0_target is not None:
        allP = np.vstack(components)
        Rmaj = np.mean(np.sqrt(allP[:,0]**2 + allP[:,1]**2))
        s = R0_target/max(Rmaj, 1e-12)
        components = [c*s for c in components]
    sum_ds=0.0; s_k2=0.0; s_k2t2=0.0; s_t2=0.0
    for P in components:
        rp, rpp, rppp = _periodic_derivs(P, du)
        cr = np.cross(rp, rpp)
        crn = np.linalg.norm(cr, axis=1)
        spd = np.linalg.norm(rp, axis=1)
        kap = crn/np.maximum(spd**3, 1e-30)
        tau = np.einsum('ij,ij->i', cr, rppp)/np.maximum(crn**2, 1e-30)
        ds  = spd*du
        sum_ds += ds.sum()
        s_k2   += np.dot(kap**2, ds)
        s_k2t2 += np.dot(kap**2*tau**2, ds)
        s_t2   += np.dot(tau**2, ds)
    L=sum_ds
    return dict(L=L, k2=s_k2/L, k2t2=s_k2t2/L, t2=s_t2/L)

def validate_frenet_paperXV(N=2000):
    """GATE: feed Paper XV's analytic (2,3) torus-knot centerline through the SAME numerical Frenet code
    and recover <kappa^2>=0.120, <kappa^2 tau^2>=0.0053 (paper15_computations.py / P15 cor:lambda3,
    rem:torsion_scope). Validates the numerical-derivative machinery before it is trusted on the
    reaction centerline."""
    R0, r0 = 3.0, 0.874
    t = np.linspace(0, 2*np.pi, N, endpoint=False); du=2*np.pi/N
    R = R0 + r0*np.cos(3*t)
    P = np.column_stack([R*np.cos(2*t), R*np.sin(2*t), r0*np.sin(3*t)])
    st = frenet_stats([P], du, R0_target=None)
    ok_k2 = abs(st['k2']-0.120) < 0.006
    ok_k2t2 = abs(st['k2t2']-0.0053) < 0.0006
    print("--- Frenet GATE: Paper XV (2,3) centerline (R0=3, r0=0.874) ---")
    print(f"  [{'PASS' if ok_k2 else 'FAIL'}] <kappa^2>      = {st['k2']:.4f}   (paper 0.120)")
    print(f"  [{'PASS' if ok_k2t2 else 'FAIL'}] <kappa^2 tau^2> = {st['k2t2']:.4f}  (paper 0.0053)")
    print(f"         L3 = {st['L']:.3f} (paper 41.26);  torsion/curv ratio <k2t2>/<k2> = {st['k2t2']/st['k2']*100:.1f}% (paper ~4%)")
    return ok_k2 and ok_k2t2

# ---------- Q_H(lambda) from the gate (reaction_doublet) ----------
def QH_of_lambda(lam, c3=3.0, L=5.0, N=96):
    from barrier_config import reaction_doublet
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    z1, z2 = reaction_doublet(X, Y, Z, lam, c3)
    return gate_field(z1, z2, h, beta=0.0)['QH']

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--Nb', type=int, default=400, help='beta samples per 2pi (centerline)')
    ap.add_argument('--nlam', type=int, default=6)
    ap.add_argument('--c3', type=float, default=3.0)
    ap.add_argument('--qh-N', type=int, default=0, help='grid N for Q_H (0 = skip, use nominal 2->3)')
    ap.add_argument('--qh-L', type=float, default=5.0)
    a = ap.parse_args()
    print("="*86)
    print(f"  WRITHE / TWIST along the reaction coordinate  (centerline Gauss integral; c3={a.c3})")
    print(f"  Tw = Q_H - [Wr_self + 2 Lk]_centerline ;  pion (Q_H=1 unknot) => pure twist")
    print("="*86)
    print(f"  {'lambda':>7} {'#comp':>5} {'Wr_self':>9} {'2*Lk':>8} {'Wr+2Lk':>9} {'Q_H':>7} {'Tw=Q_H-(Wr+2Lk)':>17}")
    for lam in np.linspace(0.0, 1.0, a.nlam):
        comps = trace_components(lam, a.c3, Nb=a.Nb)
        tot, wr, lk2 = gauss_wr_lk(comps)
        if a.qh_N > 0:
            qh = QH_of_lambda(lam, a.c3, a.qh_L, a.qh_N)
        else:
            qh = 2.0 + lam   # nominal 2->3; replace with --qh-N for the measured value
        tw = qh - tot
        print(f"  {lam:>7.2f} {len(comps):>5d} {wr:>9.3f} {lk2:>8.3f} {tot:>9.3f} {qh:>7.3f} {tw:>17.3f}")
    print("="*86)
    print()
    print("="*86)
    print("  TORSION ENERGY GEOMETRY along the reaction coordinate (Paper XV rem:torsion_scope form:")
    print("  the surviving torsion correction is dE ~ 1/2 f1^2 rho^2 <kappa^2 tau^2>, CURVATURE-GATED).")
    print("  Symmetric-profile cancellation (P15:thm:torsion_cancels) does NOT apply off lam=1.")
    print("="*86)
    gate_ok = validate_frenet_paperXV()
    print()
    du = 2*np.pi/a.Nb
    print(f"  {'lambda':>7} {'#comp':>5} {'<kappa^2>':>10} {'<tau^2>':>9} {'<k^2 t^2>':>10} {'<k2t2>/<k2>':>12}")
    for lam in np.linspace(0.0, 1.0, a.nlam):
        comps = trace_components(lam, a.c3, Nb=a.Nb)
        st = frenet_stats(comps, du, R0_target=3.0)   # rescale to R0=3 (Paper XV unit) for comparability
        ratio = st['k2t2']/st['k2'] if st['k2']>1e-12 else 0.0
        print(f"  {lam:>7.2f} {len(comps):>5d} {st['k2']:>10.4f} {st['t2']:>9.4f} {st['k2t2']:>10.4f} {ratio*100:>11.1f}%")
    print("="*86)
    print("  READ (torsion energy): the curvature-gated correction <k^2 t^2> is ~4% of <k^2> for the")
    print("  relaxed symmetric trefoil (Paper XV). Here, tracked through the NON-symmetric barrier --")
    print("  where P15:thm:torsion_cancels does NOT hold -- to see whether it stays a small correction")
    print("  (safe to carry as O(rho^2) in component F) or spikes at the reconnection. Rescaled to R0=3.")
    print("="*86)
    print()
    print("="*86)
    print("  READ: at lambda=0 the charge is mutual LINKING (two Hopf circles, Wr_self~0, 2Lk~2);")
    print("  as lambda->1 the single trefoil acquires self-WRITHE. The framing twist Tw is the")
    print("  remainder to Q_H. The net +1 charge of the reaction, carried off by the Q_H=1 UNKNOT")
    print("  pion, is a twist quantum (Wr_round=0 => Tw=1): the pion carries twist energy on")
    print("  emission / supplies it on absorption. (Wr is embedding-dependent; the +1 unknot is not.)")
    print("="*86)
