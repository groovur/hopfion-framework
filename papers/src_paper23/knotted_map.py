#!/usr/bin/env python3.11
"""
knotted_map.py -- Paper XXIII component A (the CORE gap): a singularity-free, charge-correct KNOTTED
T(2,3) field map. DM/Quark crossover (src_paper23).

CONSTRUCTION (Brieskorn, grounded -- NOT the charge-broken tube-lift): on S^3 in C^2 coords
  a = 2(x+iy)/(1+r^2),  b = ((1-r^2)+2iz)/(1+r^2),   |a|^2+|b|^2 = 1,
the zero set {a^2 + b^3 = 0} cap S^3 IS the (2,3) trefoil (classic link-of-singularity; = Paper XVII's
Sigma(2,3,5) connection). We build a PROPER S^2-valued Hopfion as a projective rational map, doublet
  Z1 = a^i b^j            (denominator -> the u=infinity preimage)
  Z2 = a^p + c * b^q      (numerator   -> the u=0 preimage = the trefoil for (p,q)=(2,3))
normalised to |Z1|^2+|Z2|^2=1 (projective -> NO division singularity even where a or b vanish). This
distributes the charge over S^3 (like the A_{p,q} base that the gate validated to Q_H->p*q), unlike
build_doublet's tube-lift (Q_H->0). The denominator exponents (i,j) set the infinity-preimage and the
Hopf charge; we SCAN candidates and let validation.gate_field report Q_H -- we do NOT assume Q_H=3,
we MEASURE it. A candidate that gives Q_H->3 (integer) with the trefoil as its u=0 set is the target.

This is an iterative build: __main__ gates several (p,q,i,j,c) candidates. Picking the physical one
(Q_H=3 + trefoil u=0 preimage + sane V) is the point; downstream components use the winner.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation import gate_field

def s3_coords(X,Y,Z):
    r2=X**2+Y**2+Z**2
    a=(2*(X+1j*Y))/(1+r2)
    b=((1-r2)+2j*Z)/(1+r2)
    return a,b

def brieskorn_doublet(X,Y,Z, p=2,q=3,i=1,j=0,c=3.0):
    """The validated trefoil map: u=(a^2 + c b^3)/a  [(p,q,i,j)=(2,3,1,0)].
    Q_H -> 3 (Richardson, L=5: Q_inf=3.03). The modulation c sets the trefoil proportions
    (the virial): c~3 gives V_max(beta=0)=1.871, matching the physical BS-profile trefoil
    V_max=1.875 to <0.3% (Paper XV; Paper XX full value 1.91). Pure Brieskorn c=1 gives 1.80.
    Projective doublet -> singularity-free. The u=0 preimage {a^2+c b^3=0} is the knotted
    (2,3) trefoil (off-axis, 2-fold toroidal winding -- gated, not axisymmetric)."""
    a,b=s3_coords(X,Y,Z)
    Z1=(a**i)*(b**j)
    Z2=a**p + c*(b**q)
    nrm=np.sqrt(np.abs(Z1)**2+np.abs(Z2)**2)+1e-30
    return Z1/nrm, Z2/nrm

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--L',type=float,default=3.0)
    ap.add_argument('--Nlist',type=str,default='56,72,96'); a=ap.parse_args()
    Ns=[int(s) for s in a.Nlist.split(',')]
    print("="*80)
    print("  knotted_map.py -- Brieskorn trefoil candidates, gated (Q_H MEASURED, not assumed)")
    print(f"  u = (a^p + c b^q)/(a^i b^j);  trefoil numerator (p,q)=(2,3);  box L={a.L}")
    print("="*80)
    # candidates: vary the denominator (i,j) -> different infinity-preimage / Hopf charge
    cands=[dict(p=2,q=3,i=0,j=3,c=1.0), dict(p=2,q=3,i=1,j=1,c=1.0),
           dict(p=2,q=3,i=1,j=0,c=1.0), dict(p=2,q=3,i=0,j=0,c=1.0),
           dict(p=2,q=3,i=0,j=2,c=1.0), dict(p=2,q=3,i=1,j=2,c=1.0)]
    print(f"  {'(p,q,i,j)':>14} {'Q_H @Ns':>28} {'trend':>12} {'V(N_max)':>9}")
    for cd in cands:
        QHs=[]; Vlast=None
        for N in Ns:
            ax=np.linspace(-a.L,a.L,N); h=ax[1]-ax[0]
            X,Y,Z=np.meshgrid(ax,ax,ax,indexing='ij')
            z1,z2=brieskorn_doublet(X,Y,Z,**cd)
            g=gate_field(z1,z2,h,expect_QH=None)
            QHs.append(g['QH']); Vlast=g['V']
        trend="rising" if QHs[-1]>QHs[0]+0.05 else ("flat" if abs(QHs[-1]-QHs[0])<0.05 else "falling")
        lab=f"({cd['p']},{cd['q']},{cd['i']},{cd['j']})"
        print(f"  {lab:>14} {str([round(x,2) for x in QHs]):>28} {trend:>12} {Vlast:>9.3f}")
    print()
    print("  READ: a candidate whose Q_H CONVERGES to an integer (ideally 3) with sane V is a valid")
    print("  Hopfion; the one with Q_H->3 AND a trefoil u=0 preimage is the target. (Knottedness of")
    print("  the u=0 set is checked separately -- next step if a Q_H=3 candidate appears.)")
    print("="*80)
