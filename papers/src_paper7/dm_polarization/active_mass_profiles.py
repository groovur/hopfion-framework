#!/usr/bin/env python3.11
"""
active_mass_profiles.py -- the PINNED-HALO ACTIVE-MASS calculation (notes (ccc), follows (bbb)).

Question: in proper GR, does the baryon-sourced director-strain halo exert the NEWTONIAN PULL that
flat rotation curves need, or only LENS? Slow-matter dynamics couples to the Komar/active mass
density  rho + sum p_i  (weak-field R_00 = 4 pi G (rho + p_r + 2 p_t)); light (lensing) couples to
the energy density + the solid-angle deficit. So compare, for the director hedgehog n=r_hat:
   M_enc(r)    = INT rho 4 pi r^2 dr          (energy density = P7's flat-curve LENSING mass)
   M_active(r) = INT (rho + sum p_i) 4 pi r^2 dr   (the NEWTONIAN DYNAMICAL mass)

General result (verified, sympy; Euler on a (grad n)^{2k} term): a homogeneous degree-2k gradient
energy has  rho + sum p_i = (2k-2) f.
   Frank (k=1, quadratic):   rho + sum p = 0         -> global monopole: lenses, NO pull.
   Faddeev J4 (k=2, quartic): rho + sum p = 2 f4     -> positive active mass (the ONLY dynamical pull).
With the hedgehog scalings f2 ~ 1/r^2, f4 ~ 1/r^4:
   M_enc[Frank]  ~ INT (1/r^2) r^2 dr  ~ R        (flat-curve amplitude -> LENSING is flat-curve-like)
   M_active[J4]  ~ INT (1/r^4) r^2 dr  -> const   (CONCENTRATED -> DYNAMICS is Keplerian, NOT flat)

=> in standard GR the strain gives flat-curve LENSING but Keplerian DYNAMICS: lensing >> dynamics.
This CONTRADICTS (a) P7's flat-rotation-curve claim (dynamical) and (b) KiDS, where lensing RAR =
kinematic RAR (Brouwer 2021). CAVEAT [O]: this applies the standard-GR matter-source formula to the
director. P7 asserts the Frank energy "gravitates" as energy density rho (pressureless), but the
director is part of P7's EMERGENT-gravity sector, and P7 never derives how the director strain
enters the emergent Einstein equations (as rho, or as rho+sum p). The resolution lives in that
underived coupling; the result here flags the gap and shows the naive-GR answer fails the data.
"""
import sympy as sp

x,y,z = sp.symbols('x y z', real=True); r=sp.sqrt(x*x+y*y+z*z); X=[x,y,z]
n=[x/r,y/r,z/r]
d=[[sp.diff(n[a],X[i]) for i in range(3)] for a in range(3)]
def cross(u,v): return [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]

# Frank
f2=sp.Rational(1,2)*sum(d[a][i]**2 for a in range(3) for i in range(3))
T2=[[sum(d[a][i]*d[a][j] for a in range(3))-(f2 if i==j else 0) for j in range(3)] for i in range(3)]
s2=sp.simplify(f2+sum(T2[i][i] for i in range(3)))
# Faddeev J4
din=[[d[a][i] for a in range(3)] for i in range(3)]
f4=sp.simplify(sp.Rational(1,2)*sum((sum(n[a]*cross(din[i],din[j])[a] for a in range(3)))**2
                                    for i in range(3) for j in range(3)))
R=sp.symbols('R',positive=True); rr=sp.symbols('r',positive=True)
f2r=sp.simplify(f2.subs({y:0,z:0})); f4r=sp.simplify(f4.subs({y:0,z:0}))
Menc =sp.integrate(f2r.subs(x,rr)*4*sp.pi*rr**2,(rr,1,R))
Mact =sp.integrate(2*f4r.subs(x,rr)*4*sp.pi*rr**2,(rr,1,R))

if __name__=='__main__':
    print("Director hedgehog, standard-GR active mass vs energy (lensing) mass:")
    print(f"  Frank:   rho = {sp.simplify(f2)},  rho+sum_p = {s2}  (k=1 -> 0: lenses, no pull)")
    print(f"  Faddeev: f4  = {f4},  rho+sum_p = 2 f4 = {sp.simplify(2*f4)}  (k=2 -> active mass)")
    print(f"  M_enc[Frank energy]   = INT rho 4pi r^2 dr          = {sp.simplify(Menc)}   (~R: LENSING flat)")
    print(f"  M_active[J4]          = INT (rho+sum_p) 4pi r^2 dr  = {sp.simplify(Mact)}   (->const: DYNAMICS Keplerian)")
    print("\n  => lensing(flat, ~R) >> dynamics(concentrated). Tension w/ flat rotation curves + KiDS")
    print("     (lensing RAR = kinematic RAR). Resolution = P7's underived director->metric coupling. [O]")
