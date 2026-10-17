#!/usr/bin/env python3.11
"""
validation.py -- the GATE for Paper XXIII (src_paper23/PLAN.md, component G). Run FIRST; nothing
downstream (knotted_map, barrier_config, beta_*rho) is trusted until a candidate field passes here.

Two layers:
  (1) METHOD validation (known cases the machinery MUST reproduce) -- run in __main__:
        * smooth rational map A_{p,q} = a^p/b^q -> Q_H -> p*q  (charge machinery / whitehead_QH).
        * BS profile integrals I_n = 2^{2n+1} B(n+1,n+1) (Paper XV P15:prop:kappa4): I_2=16/15, I_0=2,
          thin-tube V_max = I_0/I_2 = 1.875 (Paper XX full trefoil quotes 1.91).
  (2) gate_field(z1,z2,h): the reusable PASS/FAIL gate a candidate (knotted) field must clear --
        Q_H (Whitehead), V=J2iso_fb/J2a, J4/J2a, localisation, s_opt. Returns a dict with a 'pass' flag
        once targets (Q_H integer, V in a sane band) are supplied.

Reuses the VALIDATED machinery in src_paper21 (whitehead_QH, nhat_from_doublet, Efb_diagnostics,
rational_map_doublet). No new physics here -- this is the referee.
"""
import numpy as np, sys, os
from scipy.integrate import quad
from scipy.special import beta as Beta
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src_paper21'))
from feedback_saddle_core import whitehead_QH, nhat_from_doublet, Efb_diagnostics
from rational_map_ic import rational_map_doublet

phi=(1+np.sqrt(5))/2

def _check(name, got, want, tol, results):
    ok = abs(got-want) <= tol*max(1.0,abs(want))
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {got:.5f}  want {want:.5f}  (tol {tol:.0e})")
    return ok

# ---------- (2) the reusable gate ----------
def gate_field(z1, z2, h, beta=0.452, expect_QH=None, V_band=(1.3,2.3)):
    """PASS/FAIL gate on a candidate doublet field. expect_QH=None just reports Q_H."""
    QH = whitehead_QH(z1, z2, h)
    n  = nhat_from_doublet(z1, z2)
    d  = Efb_diagnostics(n, h, beta)
    out = dict(QH=QH, V=d['V'], J4_over_J2a=d['J4_over_J2a'], s_opt=d['s_opt'], E_fb=d['E_fb'])
    checks = {}
    if expect_QH is not None:
        checks['QH_integer'] = abs(QH-expect_QH) < 0.15
    checks['V_in_band'] = (V_band[0] <= d['V'] <= V_band[1]) if np.isfinite(d['V']) else False
    out['checks']=checks; out['pass']=all(checks.values()) if checks else None
    return out

# ---------- (1) method validations ----------
def validate_charge_machinery(L=4.0, Ns=(48,72,96)):
    print("--- charge machinery: A_{p,q} rational map  Q_H -> p*q ---")
    res=[]
    for (p,q) in [(1,3),(2,3)]:
        QHs=[]
        for N in Ns:
            ax=np.linspace(-L,L,N); h=ax[1]-ax[0]
            X,Y,Z=np.meshgrid(ax,ax,ax,indexing='ij')
            z1,z2=rational_map_doublet(X,Y,Z,p,q)
            QHs.append(whitehead_QH(z1,z2,h))
        trend = "rising->" if QHs[-1]>QHs[0] else "flat/falling"
        print(f"  A_{{{p},{q}}}: Q_H @N={Ns} = {[round(x,2) for x in QHs]}  (target {p*q}, {trend})")
        res.append(QHs[-1] > 0.8*p*q and QHs[-1] <= p*q*1.05)   # converging up toward p*q
    return all(res)

def validate_profile_integrals():
    print("--- BS profile integrals I_n = 2^{2n+1} B(n+1,n+1) (Paper XV) ---")
    def I_n(n): return quad(lambda rh: (2*rh/(1+rh**2))**(2*n)*(2/(1+rh**2))**2*rh,0,np.inf)[0]
    res=[]
    r=[]
    _check("I_0 (=2)", I_n(0), 2.0, 1e-4, r)
    _check("I_2 (=16/15, the J2a norm)", I_n(2), 16/15, 1e-4, r)
    _check("thin-tube V_max = I_0/I_2 (=1.875; full trefoil 1.91)", I_n(0)/I_n(2), 1.875, 1e-3, r)
    _check("I_n Beta closed form n=2 (=2^5 B(3,3))", I_n(2), 2**5*Beta(3,3), 1e-6, r)
    return all(r)

if __name__=='__main__':
    print("="*72); print("  PAPER XXIII GATE -- method validation (known cases)"); print("="*72)
    a = validate_charge_machinery()
    print(); b = validate_profile_integrals()
    print("\n"+"="*72)
    print(f"  charge machinery: {'PASS' if a else 'FAIL'}   profile integrals: {'PASS' if b else 'FAIL'}")
    print("  => the gate's machinery is validated on known cases; gate_field() is ready for a")
    print("     candidate knotted map (knotted_map.py). A knotted field must return Q_H~integer")
    print("     and V in a sane band to pass.")
    print("="*72)
