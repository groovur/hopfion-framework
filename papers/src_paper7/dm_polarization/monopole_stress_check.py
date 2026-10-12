#!/usr/bin/env python3.11
"""
monopole_stress_check.py -- gravitational-slip check for the director-strain halo (notes (bbb)).

Question: does the baryon-sourced director hedgehog (P7 flat-curve halo, n = r_hat, |grad n|~1/r)
lens the same as it pulls, or is there a lensing-vs-dynamics SLIP? The slip is set by the
ANISOTROPIC STRESS of the gravitating Frank energy -- NOT the polaron inertial tensor (that is the
SOURCE baryon's inertia, not the halo test stars orbit; (aaa)'s polaron slip is RETRACTED).

Field-theory stress-energy of a static director texture with energy density f(grad n):
   T00 = f,   Tij = (df/d(d_i n_a)) d_j n_a - dij f.
Compute rho = T00 and the sum of principal pressures sum p_i = T^i_i (the Komar/active-mass
density source: weak-field R_00 = 4 pi G (rho + sum p_i), so slow-matter acceleration a = -grad Phi
has source rho + sum p_i). Also the radial/tangential split p_r, p_t.

RESULT (verified, sympy): BOTH one-constant Frank and pure splay give
   rho = +c/r^2,  p_r = -c/r^2,  p_t = 0,  rho + sum p_i = 0.
This is the Barriola-Vilenkin GLOBAL MONOPOLE: zero active gravitational mass (NO Newtonian force
on test particles) but a solid-angle deficit that DEFLECTS light (lenses). So a FREE director
texture lenses WITHOUT exerting Newtonian force -- the opposite extreme of normal matter, and in
tension with P7's flat-rotation-curve mechanism (which needs the strain to PULL, treating rho=u as
a Newtonian source and ignoring the p_r=-u cancellation). P7's "pressureless" (c_dir->0) is the
absence of JEANS/kinetic pressure support; the STATIC elastic stress p_r=-u is a different thing
and is present. Likely resolution [O]: the halo is NOT a free texture -- it is pinned by baryon
boundary conditions and held by the global a0 condition (single_galaxy_winding_analysis); the
binding/pinning stress contributes to T_munu and can restore a nonzero active mass. Unmodeled.
"""
import sympy as sp

x,y,z = sp.symbols('x y z', real=True)
r = sp.sqrt(x*x+y*y+z*z); X=[x,y,z]
n = [x/r, y/r, z/r]                                   # hedgehog n = r_hat
dn = [[sp.diff(n[a],X[i]) for i in range(3)] for a in range(3)]   # dn[a][i] = d_i n_a
rh = [x/r, y/r, z/r]

def report(name, f, T):
    rho = sp.simplify(f)
    trT = sp.simplify(sum(T[i][i] for i in range(3)))
    pr  = sp.simplify(sum(rh[i]*rh[j]*T[i][j] for i in range(3) for j in range(3)))
    pt  = sp.simplify((trT-pr)/2)
    print(f"  {name:12}: rho={rho}  p_r={pr}  p_t={pt}  rho+sum_p={sp.simplify(rho+trT)}")

if __name__=='__main__':
    print("Director hedgehog stress-energy (field-theory T_munu):")
    K = sp.symbols('K', positive=True)
    f1 = sp.Rational(1,2)*K*sum(dn[a][i]**2 for a in range(3) for i in range(3))
    T1 = [[K*sum(dn[a][i]*dn[a][j] for a in range(3)) - (f1 if i==j else 0) for j in range(3)] for i in range(3)]
    report("one-constant", f1, T1)
    Ks = sp.symbols('K_s', real=True)
    divn = sum(dn[a][a] for a in range(3))
    fs = sp.Rational(1,2)*Ks*divn**2
    Ts = [[Ks*divn*dn[i][j] - (fs if i==j else 0) for j in range(3)] for i in range(3)]
    report("pure-splay", fs, Ts)
    print("\n  => rho+sum_p = 0 (global monopole): lenses (deficit), NO Newtonian force. See docstring.")
