#!/usr/bin/env python3.11
"""
rotating_active_mass.py -- the active (gravitating) mass of a MOVING/rotating director (notes (eee)).
Follows the user's picture: the geometry rotates locally; the (light) lens follows the static
geometry; the halo "drags after" -- and that drag is REACTIVE/elastic (c_dir->0: frozen, no
radiation, NO dissipation; cluster_wake_calc.py), NOT viscous. So it can SUSTAIN orbits.

Komar active-mass density (slow-matter source): rho + sum p_i. For the director field with
kinetic chi(d_t n)^2 and Frank K(grad n)^2 (Hilbert T_munu, mostly-minus):
   rho + sum p_i = 2 chi (d_t n)^2.
   -> STATIC texture (d_t n = 0): 0   (global monopole: lenses via deficit, NO Newtonian pull).
   -> ROTATING/moving director (d_t n != 0): +2 chi (d_t n)^2 > 0   (PULLS; conservative).
The gradient (texture) part CANCELS in rho+sum p (the (bbb) global-monopole result); ALL the
dynamical pull comes from the TIME-DEPENDENCE -- the orbiting/dispersing baryons dragging the
director. This IS P7's polaron momentum mass: for a comoving distortion n(x-Vt), d_t n = -V.grad n,
so 2 chi (d_t n)^2 integrated = 2 chi INT (V.grad n)^2 = 2 M_ij V_i V_j with M_ij=chi INT d_i n d_j n
(P7:rem:polaron). Static lenses, moving pulls -- one mechanism.

PROFILE (self-consistency sketch, NOT a derivation): active density ~ 2 chi <v^2> |grad n|^2 where
<v^2> is the baryon kinetic scale (rotation v_c^2 OR dispersion sigma^2) and |grad n|^2 ~ 1/r^2
(disclination). Then M_active(r) = INT 2 chi <v^2>/r^2 4pi r^2 dr ~ r  -> M_enc proportional r =
FLAT curves, and proportional to the SAME |grad n|^2 integral as the lensing energy -> lensing =
kinematic RAR (KiDS, Brouwer 2021). chi=K/c_dir^2 -> infinity as c_dir->0, so the cold director
gives a large mass from even slow motion. CAVEATS [O]: assumes the director co-moves at the baryon
kinetic scale (why/how much?); the r-profile and the a0 normalization are not derived; rotation-
supported vs dispersion-supported systems (disks vs dSph/ellipticals) need the <v^2> to be the
local kinetic scale, not bulk rotation, or dSphs (DM-rich, low rotation) would be mispredicted.
"""
import sympy as sp

if __name__=='__main__':
    pt,px,py,pz,chi,K=sp.symbols('n_t n_x n_y n_z chi K',positive=True)
    grad2=px**2+py**2+pz**2
    L  = sp.Rational(1,2)*chi*pt**2 - sp.Rational(1,2)*K*grad2      # Lagrangian density
    rho= sp.Rational(1,2)*chi*pt**2 + sp.Rational(1,2)*K*grad2      # energy density T_00
    trT= K*grad2 + 3*L                                             # sum p_i = T^i_i (mostly-minus)
    aps= sp.simplify(rho+trT)
    print("Director active (Komar) mass density  rho + sum p_i =", aps)
    print("  static texture (n_t=0):   ", aps.subs(pt,0), "  -> global monopole: lenses, no pull")
    print("  rotating/moving (grad=0): ", sp.simplify(aps.subs({px:0,py:0,pz:0})), "  -> PULLS (conservative)")
    print("\n=> the pull is the TIME-DEPENDENT (drag) part = P7's polaron momentum mass; the static")
    print("   texture only lenses. Reactive (c_dir->0: no dissipation), so it sustains flat curves.")
