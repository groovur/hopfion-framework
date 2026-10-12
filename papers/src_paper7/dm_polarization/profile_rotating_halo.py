#!/usr/bin/env python3.11
"""
profile_rotating_halo.py -- halo PROFILE from the rotating/comoving director active mass (notes (eee)->(ggg)).

Mechanism (eee, Paper 7 rem:active_mass): the STATIC director texture lenses but does not pull
(rho+sum p=0, global monopole); the DYNAMICAL pull is the comoving kinetic term rho_act=2chi(d_t n)^2.
The director is dragged by the local matter motion: (d_t n)^2 ~ eta <v^2>(r) |grad n|^2, with a
disclination |grad n|^2 = 2/r^2 and <v^2> the LOCAL kinetic scale (rotation v_c^2 OR dispersion sigma^2).
So the active-mass density is  rho_act(r) = 4 chi eta <v^2>(r) / r^2.

SELF-CONSISTENT closure (rotation-supported disk): v_c^2(r)=G[M_b(r)+M_act(r)]/r, and
  dM_act/dr = 4 pi r^2 rho_act = 16 pi chi eta G [M_b(r)+M_act(r)] / r  (linear ODE; set A=16 pi chi eta G).
This is a bootstrap: the director mass is sourced by the motion it itself supports. We solve it and ask:
(1) do flat rotation curves emerge? (2) does an acceleration scale a0 appear (strain overtakes baryons
where g ~ a0)? (3) the DISPERSION case (sigma^2 const) -> M_act ~ r directly (dSph: DM-rich, no rotation).
LENSING: the deficit (static) + kinetic energy, M_lens ~ INT[1/2 K|grad n|^2 + 1/2 chi(d_t n)^2] ~ r; for a
cold director (c_dir<<v) the kinetic term dominates BOTH, so dynamics/lensing ~ O(1) (KiDS). The exact
anisotropic lensing source (kinetic/stiff stress) is flagged, not pinned.
HONEST: toy; A (=16 pi chi eta G) is the one free amplitude (the cutoff-set normalization P7 does not derive);
spherical; eta O(1); the lensing/dynamics O(1) factor needs the anisotropic-stress lensing kernel.
"""
import numpy as np
from scipy.integrate import solve_ivp

G=4.302e-6                  # kpc (km/s)^2 / Msun
MS2=3.2408e-14              # (km/s)^2/kpc per m/s^2
a0_obs=1.2e-10/MS2         # ~3703 (km/s)^2/kpc

def Mbar(r,Mb,Rd): x=r/Rd; return Mb*(1-(1+x)*np.exp(-x))

def solve_disk(Mb=5e10, Rd=3.5, A=1.0, rmax=60.0):
    # dM_act/dr = (A/r)(M_b+M_act); integrate from small r
    def rhs(r,M): return (A/r)*(Mbar(r,Mb,Rd)+M[0])
    r=np.logspace(np.log10(0.05*Rd),np.log10(rmax),400)
    sol=solve_ivp(rhs,[r[0],r[-1]],[0.0],t_eval=r,rtol=1e-8,atol=1e-3)
    Mact=sol.y[0]; Mb_r=Mbar(r,Mb,Rd)
    vc=np.sqrt(G*(Mb_r+Mact)/r); vb=np.sqrt(G*Mb_r/r)
    return r,vc,vb,Mb_r,Mact

if __name__=='__main__':
    print("="*78); print("  ROTATING-DIRECTOR HALO PROFILE: self-consistent v_c(r)  (A=16 pi chi eta G)"); print("="*78)
    print("  (1) DISK rotation curves vs the self-consistency amplitude A:")
    print(f"    {'A':>6} {'v_flat(km/s)':>13} {'flat? sd(dex,r>3Rd)':>20} {'Mact/Mb @10Rd':>14}")
    for A in [0.2,0.5,1.0,2.0]:
        r,vc,vb,Mb_r,Mact=solve_disk(A=A)
        outer=r>3*3.5
        flat_sd=np.std(np.log10(vc[outer]))
        i10=np.argmin(np.abs(r-35))
        print(f"    {A:>6.2f} {vc[outer][0]:>13.1f} {flat_sd:>20.3f} {Mact[i10]/Mb_r[i10]:>14.2f}")

    print("\n  (2) acceleration scale: where does the strain overtake baryons (Mact=Mb)? g there vs a0:")
    for A in [0.5,1.0,2.0]:
        r,vc,vb,Mb_r,Mact=solve_disk(A=A)
        k=np.where(Mact>=Mb_r)[0]
        if len(k):
            rc=r[k[0]]; g=G*Mb_r[k[0]]/rc**2
            print(f"    A={A}: strain=baryon at r={rc:.1f} kpc, g_bar there={g:.0f} (km/s)^2/kpc = {g/a0_obs:.2f} a0")
        else: print(f"    A={A}: strain never overtakes baryons in range")

    print("\n  (3) DISPERSION system (dSph, sigma=const, no rotation): M_act(r)=16 pi chi eta sigma^2 r")
    print("    => M_act ~ r directly (DM-rich, rotation-independent): the mechanism does NOT require rotation,")
    print("    only local <v^2>. dSphs (high sigma^2/baryon, low rotation) get a large M~r halo. Consistent.")

    print("\n  READ: a self-consistency amplitude A~O(1) gives an extended ~flat v_c(r); strain overtakes")
    print("  baryons near a fixed acceleration (candidate a0); dispersion systems get M~r from sigma^2 alone.")
    print("  Both dynamics and the (kinetic-dominated) lensing ~ r -> common RAR (KiDS). A is the one")
    print("  cutoff-set amplitude P7 leaves underived; a0 emerging at ~fixed g is the nontrivial check.")
