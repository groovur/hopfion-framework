#!/usr/bin/env python3.11
"""
lensing_dynamics_slip.py -- the lensing-vs-dynamics ratio for the director halo (notes (hhh)).
Follows (eee): static texture lenses (anisotropic, global monopole), comoving/rotating director pulls.

Decompose the halo stress into its two pieces and read the slip:
  STATIC texture (hedgehog): rho=1/2 K|grad n|^2, p_r=-rho, p_t=0  -> ANISOTROPIC.
     dynamics source rho+p_r+2p_t = 0 (global monopole: no pull);
     lensing: anisotropic stress -> Phi != Psi -> a solid-angle-deficit LENSING mass ~ rho_static.
  COMOVING/kinetic (d_t n != 0): rho_kin=1/2 chi (d_t n)^2, p_r=p_t=rho_kin  -> ISOTROPIC stiff (p=rho).
     isotropic pressure => Phi = Psi => LENSING mass = DYNAMICS mass (no slip; relativistic factor acts on both).

So, summing: dynamics = 4 rho_kin (only the isotropic kinetic part pulls); lensing = 4 rho_kin + rho_static
(the anisotropic static texture adds to LENSING only). With d_t n ~ v grad n and chi=K/c_dir^2,
rho_static/rho_kin = c_dir^2/v^2, hence
     lensing/dynamics = 1 + (1/4)(c_dir/v)^2.
The slip is controlled by the SAME coldness parameter c_dir/v that the director-DM sector already
requires (c_dir -> 0 for pressureless clustering): if c_dir << v_gal the slip -> 0 (weak lensing RAR =
kinematic RAR, as KiDS/Brouwer 2021 find); the residual is a DEFINITE-SIGN prediction (lensing slightly
ABOVE dynamics) of order (c_dir/v)^2 -- a future precision test. KiDS consistency thus REQUIRES the same
coldness the sector hinges on; it is not a new assumption. CAVEAT: c_dir is cutoff-set (not pinned); the
isotropic-vs-anisotropic split is the clean physics, the exact (c_dir/v) at galactic scales is not derived.
"""
import numpy as np

def slip(cdir_over_v): return 1.0 + 0.25*cdir_over_v**2

if __name__=='__main__':
    print("lensing/dynamics = 1 + (1/4)(c_dir/v)^2   [static=anisotropic lensing-only; kinetic=isotropic no-slip]")
    print(f"  {'c_dir/v':>9} {'lensing/dynamics':>18}")
    for x in [0.0,0.05,0.1,0.2,0.5,1.0,2.0]:
        print(f"  {x:>9.2f} {slip(x):>18.3f}")
    print("\n  cold director (c_dir << v_gal): ratio -> 1  => KiDS-consistent (lensing RAR = kinematic RAR).")
    print("  SAME coldness the DM sector needs for clustering; residual slip (lensing slightly >1) ~ (c_dir/v)^2,")
    print("  a definite-sign future precision test. c_dir/v cutoff-set (not pinned).")
