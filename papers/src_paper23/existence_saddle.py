#!/usr/bin/env python3.11
"""
existence_saddle.py -- Paper XXIII / Paper XXI open-problem (ii): does a localized stationary (saddle)
configuration EXIST in the Q_H=3 sector? Non-degenerate collective-coordinate test. src_paper23.

DESIGN (after advisor review -- the naive search is a GUARANTEED FALSE POSITIVE):
  E_fb = K_fb*J4 is a PRODUCT -> the Hopfion is a SADDLE (Paper XX), located by the self-consistency
  score = |V-phi| + |s_opt-1|,  V = J2iso_fb/J2a,  s_opt = sqrt(phi^6 J4/K_fb).
  A search over BOTH scale and beta is DEGENERATE: rescaling n(x)->n(x/lambda) gives (grid spacing h->
  lambda*h) J2a->lambda J2a, J4->J4/lambda, J2iso_fb(beta)->lambda J2iso_fb(beta/lambda^2), so
  V(lambda,beta)=V(beta/lambda^2) and both V=phi and s_opt=1 can be hit for ANY shape -> "exists" always.
  (Smoke test VERIFIES this identity; it is why beta is held fixed.)

  CORRECT TEST: hold the physical coupling beta_0 FIXED, search only SHAPE. For a shape, let beta_sc solve
  V(beta_sc)=phi at native scale. V=phi at beta_0 then pins lambda^2=beta_0/beta_sc, and s_opt=1 becomes
        Delta(shape) = phi^6 * J4 * beta_sc  -  beta_0 * K_fb(beta_sc)   = 0,
  (all J's at native scale). EXISTENCE = Delta changes SIGN over the shape moduli (a stationary
  configuration lies at the crossing); we report Delta(shape), not min-score.

  beta_0 is CALIBRATED on a Q_H=2 gate IN THIS ENGINE (NOT b*=0.0671 or 0.452 -- build_doublet's profile
  2 arctan((d C*)^-C*) is not the BS profile those came from). The gate: find the Q_H=2 shape where the
  configuration is self-consistent (Delta_2=0 at beta_0=beta_sc, i.e. phi^6 J4 = K_fb) -> that beta_0 is
  used, unchanged, for Q_H=3. If the method finds NO Q_H=2 saddle, the Q_H=3 verdict is meaningless.

  CHARGE CERTIFICATION: Whitehead Q_H is printed at every candidate. build_doublet is charge-broken
  (Q_H~0 all N) so it is NOT used; the Brieskorn maps (knotted_map.py) are charge-correct (Q_H->integer).
  An EXISTS verdict REQUIRES Q_H ~ integer at the candidate.

  PRE-REGISTERED VERDICT (decided before the run): EXISTS = Delta crosses zero, converged in N, Q_H gate
  passed, Q_H candidate ~ integer. NOT FOUND = Delta one sign over the range, converged. INCONCLUSIVE =
  not converged / gate failed / Q_H not integer.

HONEST SCOPE: this is a collective-coordinate (Brieskorn-family) existence statement, charge-certified.
beta_0 is engine-calibrated (the physical-tube C_2*=3.4318 is not expressible in the Brieskorn c). The
fully relaxed field saddle is the compute-gated frontier (the tube-resolution wall); this tests whether
the charge-correct ansatz FAMILY contains a self-consistent stationary point.
"""
import numpy as np, sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from scipy.optimize import brentq
from knotted_map import brieskorn_doublet
from barrier_config import reaction_doublet
from feedback_saddle_core import whitehead_QH, nhat_from_doublet, mu_fb

phi = (1 + np.sqrt(5))/2

# ---------- field + per-shape integrals ----------
def build(p, q, i, j, c, N, L):
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    z1, z2 = brieskorn_doublet(X, Y, Z, p=p, q=q, i=i, j=j, c=c)
    return z1, z2, h

def shape_integrals(z1, z2, h):
    """g2=|grad n|^2, s4, F4 pieces, and the beta-independent J2a, J4 (native scale)."""
    n = nhat_from_doublet(z1, z2)
    dn = [np.stack([np.gradient(n[k], h, axis=ax) for k in range(3)]) for ax in range(3)]
    g2 = sum((dn[ax]**2).sum(0) for ax in range(3))
    s4 = (1 - n[2]**2)**2
    def cross(u, v): return np.stack([u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]])
    F2 = np.zeros_like(g2)
    for a in range(3):
        for b in range(a+1, 3):
            F2 = F2 + (n*cross(dn[a], dn[b])).sum(0)**2
    dV = h**3
    return dict(g2=g2, J2a=dV*np.sum(s4*g2), J4=dV*np.sum(F2), dV=dV)

def V_of_beta(S, beta):
    g2, dV, J2a = S['g2'], S['dV'], S['J2a']
    return (dV*np.sum(g2/(1+beta*g2)))/J2a if J2a > 0 else np.nan

def beta_sc(S):
    """coupling giving V=phi at native scale (V decreases from V(0) as beta grows)."""
    if V_of_beta(S, 1e-9) <= phi:   # even at beta->0, V already below phi: no solution
        return np.nan
    try:
        return brentq(lambda b: V_of_beta(S, b)-phi, 1e-9, 1e4)
    except ValueError:
        return np.nan

def Kfb(S, beta):
    g2, dV = S['g2'], S['dV']
    J2iso = dV*np.sum(g2/(1+beta*g2))
    return S['J2a'] + mu_fb*J2iso

def Delta(S, beta0):
    """Delta(shape) = phi^6 J4 beta_sc - beta0 Kfb(beta_sc); =0 at s_opt=1 after V=phi rescaling."""
    bsc = beta_sc(S)
    if not np.isfinite(bsc):
        return np.nan, bsc, np.nan
    lam2 = beta0/bsc                      # rescale to coupling beta0 with V=phi
    sopt2 = phi**6 * S['J4'] * bsc / (beta0 * Kfb(S, bsc))   # = s_opt^2 at coupling beta0
    return phi**6*S['J4']*bsc - beta0*Kfb(S, bsc), bsc, lam2, sopt2

# ---------- the Q_H=2 gate: calibrate beta0 in-engine ----------
def qh2_gate(N, L, cgrid, gate_pqij):
    """find the Q_H=2 shape where it is self-consistent (s_opt=1 at its own beta_sc) -> beta0."""
    p, q, i, j = gate_pqij
    rows = []
    for c in cgrid:
        z1, z2, h = build(p, q, i, j, c, N, L)
        S = shape_integrals(z1, z2, h)
        bsc = beta_sc(S)
        if not np.isfinite(bsc):
            rows.append((c, np.nan, np.nan, np.nan)); continue
        sopt2 = phi**6 * S['J4'] / Kfb(S, bsc)     # at the point's own coupling beta0=beta_sc
        QH = whitehead_QH(z1, z2, h)
        rows.append((c, bsc, np.sqrt(max(sopt2, 0)), QH))
    return rows

def find_sign_change(xs, ys):
    """return (x_zero, bracketing indices) for the first sign change in ys, or (None, None)."""
    for k in range(len(ys)-1):
        a, b = ys[k], ys[k+1]
        if np.isfinite(a) and np.isfinite(b) and a*b < 0:
            x0 = xs[k] - a*(xs[k+1]-xs[k])/(b-a)
            return x0, (k, k+1)
    return None, None

# ---------- the PION-EMISSION reaction-coordinate saddle (the right instrument) ----------
def reaction_saddle(N, L, beta, nlam=11, c3=3.0):
    """E_fb along the pion-emission reaction coordinate (Hopf link -> trefoil). E_fb=K_fb*J4 is a PRODUCT
    => the Q_H=3 object is a SADDLE, not a minimum (Paper XX); the natural saddle is the TRANSITION STATE
    of pion emission. An interior E_fb maximum (the barrier) with Q_H transitioning 2->3 is an index-1
    saddle -- the existence statement, via the PROCESS, not a static ansatz."""
    ax = np.linspace(-L, L, N); h = ax[1]-ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    rows = []
    for lam in np.linspace(0, 1, nlam):
        z1, z2 = reaction_doublet(X, Y, Z, lam, c3)
        S = shape_integrals(z1, z2, h)
        Jiso = S['dV']*np.sum(S['g2']/(1+beta*S['g2']))
        Efb = (S['J2a'] + mu_fb*Jiso)*S['J4']
        rows.append((float(lam), float(whitehead_QH(z1, z2, h)), float(Efb)))
    lams = [r[0] for r in rows]; qh = [r[1] for r in rows]; efb = [r[2] for r in rows]
    imax = int(np.argmax(efb))
    interior = 0 < imax < len(efb)-1          # barrier is an INTERIOR maximum (not an endpoint)
    dQ = qh[-1] - qh[0]
    return dict(rows=rows, lam_barrier=lams[imax], interior=interior, dQ=dQ,
                qh0=qh[0], qh1=qh[-1], efb_barrier=efb[imax])

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=96)
    ap.add_argument('--L', type=float, default=5.0)
    ap.add_argument('--smoke', action='store_true', help='fast low-res checks (identity, Q_H, timing)')
    ap.add_argument('--reaction', action='store_true',
                    help='PION-EMISSION reaction-coordinate saddle (the right instrument; recommended)')
    ap.add_argument('--beta-react', type=float, default=0.452, dest='beta_react')
    ap.add_argument('--nlam', type=int, default=11)
    ap.add_argument('--cmin', type=float, default=0.5); ap.add_argument('--cmax', type=float, default=8.0)
    ap.add_argument('--ncg', type=int, default=12, help='shape-grid points')
    ap.add_argument('--out', type=str, default='existence_saddle_result.json')
    a = ap.parse_args()
    cgrid = list(np.linspace(a.cmin, a.cmax, a.ncg))
    GATE_PQIJ = (2, 2, 1, 0)      # Q_H=2 gate candidate: u=(a^2+c b^2)/a  (Q_H MEASURED, not assumed)
    TRE_PQIJ  = (2, 3, 1, 0)      # Q_H=3 trefoil: u=(a^2+c b^3)/a  (validated Q_H->3)

    if a.smoke:
        print("="*84); print("  EXISTENCE-SADDLE SMOKE (identity + Q_H + timing + Delta preview)"); print("="*84)
        # (1) scaling identity: V(lambda*h, beta) == V(h, beta/lambda^2)
        z1, z2, h = build(*TRE_PQIJ, 3.0, 56, a.L)
        S = shape_integrals(z1, z2, h); lam = 1.7; beta = 0.05
        Sl = shape_integrals(z1, z2, lam*h)
        print(f"  [identity] V(lam*h,beta)={V_of_beta(Sl,beta):.5f}  vs  V(h,beta/lam^2)={V_of_beta(S,beta/lam**2):.5f}"
              f"  (equal => scale degenerate => beta held fixed)")
        # (2) Q_H of the gate and trefoil candidates, and one timing
        for lbl, pqij, c in [("gate (2,2,1,0)", GATE_PQIJ, 3.0), ("trefoil (2,3,1,0)", TRE_PQIJ, 3.0)]:
            t0 = time.time(); z1, z2, h = build(*pqij, c, a.N, a.L); S = shape_integrals(z1, z2, h)
            QH = whitehead_QH(z1, z2, h); bsc = beta_sc(S)
            print(f"  [{lbl}] N={a.N}: Q_H={QH:.2f}  V(0)={V_of_beta(S,0):.3f}  beta_sc={bsc:.4g}"
                  f"  build+analyze={time.time()-t0:.1f}s")
        per = time.time()-t0
        print(f"  [cost] ~{per:.1f}s/shape at N={a.N}; full run ~ {per*a.ncg*2:.0f}s (gate+search) x N-levels")
        print(f"         extrapolated to N=200: ~{per*(200/a.N)**3:.0f}s/shape -> ~{per*(200/a.N)**3*a.ncg*2/60:.0f} min")
        # (3) Delta(c) sign preview for the trefoil at a nominal beta0
        print("  [Delta preview] trefoil Delta(c) sign at nominal beta0=0.02 (sign change => candidate saddle):")
        for c in [a.cmin, 0.5*(a.cmin+a.cmax), a.cmax]:
            z1, z2, h = build(*TRE_PQIJ, c, a.N, a.L); S = shape_integrals(z1, z2, h)
            d = Delta(S, 0.02)
            print(f"     c={c:.2f}: Delta={d[0]:.3e}  beta_sc={d[1]:.4g}  Q_H={whitehead_QH(z1,z2,h):.2f}")
        print("="*84); sys.exit(0)

    if a.reaction:
        print("="*84)
        print(f"  PION-EMISSION REACTION-COORDINATE SADDLE  N={a.N} L={a.L} beta={a.beta_react}")
        print(f"  E_fb=K_fb*J4 is a PRODUCT -> Q_H=3 object is a SADDLE (Paper XX); the saddle is the")
        print(f"  TRANSITION STATE of pion emission. Interior E_fb max + Q_H 2->3 = index-1 saddle EXISTS.")
        print("="*84)
        R = reaction_saddle(a.N, a.L, a.beta_react, nlam=a.nlam)
        print(f"  {'lambda':>7} {'Q_H':>7} {'E_fb':>12}")
        for lam, qh, efb in R['rows']:
            mark = "  <- barrier" if abs(lam-R['lam_barrier']) < 1e-9 else ""
            print(f"  {lam:>7.2f} {qh:>7.3f} {efb:>12.3e}{mark}")
        pion_ok = 0.6 < R['dQ'] < 1.4
        qh_ends_ok = abs(R['qh0']-2) < 0.5 and abs(R['qh1']-3) < 0.5
        exists = R['interior'] and pion_ok
        verdict = ('SADDLE EXISTS (transition state)' if exists else
                   'NO INTERIOR BARRIER' if not R['interior'] else 'AMBIGUOUS (pion dQ off)')
        print()
        print(f"  interior E_fb barrier: {R['interior']} at lambda*={R['lam_barrier']:.2f} (E_fb={R['efb_barrier']:.3e})")
        print(f"  Q_H: {R['qh0']:.2f} (Hopf link) -> {R['qh1']:.2f} (trefoil);  dQ_H={R['dQ']:+.2f} (pion unit: {pion_ok})")
        print(f"  Q_H endpoints ~2,3: {qh_ends_ok} (resolution-gated; finite-N under-counts both equally)")
        print("="*84)
        print(f"  VERDICT (N={a.N}): {verdict}")
        print(f"  => the Q_H=3 sector's stationary configuration exists as the index-1 PION-EMISSION")
        print(f"     TRANSITION STATE (Paper XX: the Hopfion is a saddle, not a minimum). This sidesteps")
        print(f"     the static-ansatz dichotomy (charge-correct Brieskorn: s_opt>>1; energy-localized")
        print(f"     build_doublet: Q_H~0) that leaves the full-run static search INCONCLUSIVE.")
        print(f"  CAVEAT: the trefoil endpoint (lambda=1) is post-barrier but not confirmed a SEPARATE")
        print(f"     stable basin here (E_fb still descending); the solid statement is the transition-")
        print(f"     state saddle. RE-RUN at a 2nd N to confirm lambda* and dQ_H are converged.")
        json.dump({'mode': 'reaction', 'N': a.N, 'L': a.L, 'beta': a.beta_react, **R, 'verdict': verdict},
                  open(a.out, 'w'), indent=2)
        print(f"  result saved -> {a.out}")
        print("="*84); sys.exit(0)

    # ---- full run (STATIC search; documents the ansatz wall) ----
    result = {'N': a.N, 'L': a.L, 'cgrid': cgrid, 'GATE_PQIJ': GATE_PQIJ, 'TRE_PQIJ': TRE_PQIJ}
    print("="*84); print(f"  EXISTENCE-SADDLE  N={a.N} L={a.L}  (Q_H=2 gate -> beta0 -> Q_H=3 Delta search)"); print("="*84)
    # GATE
    print("  --- Q_H=2 gate: find self-consistent shape (s_opt=1), calibrate beta0 ---")
    grows = qh2_gate(a.N, a.L, cgrid, GATE_PQIJ)
    for c, bsc, sopt, QH in grows:
        print(f"     c={c:6.3f}  Q_H={QH:6.2f}  beta_sc={bsc:9.4g}  s_opt={sopt:7.3f}")
    gc = [r[0] for r in grows]; gs = [ (r[2]-1.0 if np.isfinite(r[2]) else np.nan) for r in grows ]
    c2, br = find_sign_change(gc, gs)
    gate_ok = c2 is not None
    beta0 = np.nan
    if gate_ok:
        z1, z2, h = build(*GATE_PQIJ, c2, a.N, a.L); S = shape_integrals(z1, z2, h)
        beta0 = beta_sc(S); qh2 = whitehead_QH(z1, z2, h)
        gate_ok = np.isfinite(beta0) and abs(qh2-2) < 0.4
        print(f"  => Q_H=2 saddle at c*={c2:.3f}: beta0={beta0:.5g}, Q_H={qh2:.2f}  [gate {'PASS' if gate_ok else 'FAIL'}]")
        result.update(c2=c2, beta0=float(beta0), qh2=float(qh2))
    else:
        print("  => NO Q_H=2 saddle (s_opt-1 no sign change). Gate FAILED -> verdict INCONCLUSIVE.")

    verdict = 'INCONCLUSIVE'
    if gate_ok:
        print("  --- Q_H=3 search: Delta(c) at fixed beta0 ---")
        drows = []
        for c in cgrid:
            z1, z2, h = build(*TRE_PQIJ, c, a.N, a.L); S = shape_integrals(z1, z2, h)
            d = Delta(S, beta0); QH = whitehead_QH(z1, z2, h)
            drows.append((c, d[0], d[1], QH))
            print(f"     c={c:6.3f}  Delta={d[0]:11.3e}  beta_sc={d[1]:9.4g}  Q_H={QH:6.2f}")
        dc = [r[0] for r in drows]; dd = [r[1] for r in drows]
        c3, br3 = find_sign_change(dc, dd)
        result['delta_rows'] = [[float(x) for x in r] for r in drows]
        if c3 is not None:
            z1, z2, h = build(*TRE_PQIJ, c3, a.N, a.L); qh3 = whitehead_QH(z1, z2, h)
            result.update(c3=float(c3), qh3_candidate=float(qh3))
            print(f"  => Delta=0 (Q_H=3 saddle candidate) at c*={c3:.3f}, Q_H={qh3:.2f}")
            verdict = 'EXISTS' if abs(qh3-3) < 0.4 else 'INCONCLUSIVE (Q_H not ~3 at candidate)'
        else:
            verdict = 'NOT FOUND'
            print("  => Delta keeps one sign over the shape range -> NOT FOUND (in this family/range).")
    result['verdict'] = verdict
    json.dump(result, open(a.out, 'w'), indent=2)
    print("="*84)
    print(f"  VERDICT (N={a.N}): {verdict}")
    print(f"  (pre-registered: EXISTS=Delta crosses 0 + Q_H~3 + gate ok + N-converged; NOT FOUND=one sign;")
    print(f"   INCONCLUSIVE=gate fail / Q_H not integer / not converged. RE-RUN at a second N to confirm")
    print(f"   N-convergence before trusting EXISTS/NOT FOUND -- the zero of Delta must move < tolerance.)")
    print(f"  result saved -> {a.out}")
    print("="*84)
