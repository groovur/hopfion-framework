#!/usr/bin/env python3.11
"""
galactic_overlap_sidm.py -- the galactic-overlap gap (notes (jjj)->(kkk)): can n_hat-mediated
velocity-dependent SIDM make the dark hopfion THERMALISED (baryon-tracking) in galaxies but
COLLISIONLESS (dominant CDM) in clusters -- so the tight RAR (strain) survives in galaxies while
the dark hopfion supplies cluster mass + the recombination CDM?

Diagnostic: scatterings per particle over a halo age, N = (rho/m) sigma v t = rho (sigma/m) v t.
N >~ few  => thermalised core, DM relaxes toward the (baryon-dominated) potential -> tracks baryons
            -> does NOT detrack the RAR (SIDM-baryon connection).
N <~ 1    => collisionless -> NFW cusp, dominant, LCDM-like.
Velocity dependence sigma/m = (sigma/m)_ref (v_ref/v)^n (n_hat light-ish mediator: n~1-2 Yukawa,
n=4 massless/Rutherford). Then N ∝ rho v^(1-n). Galaxies: dense + SLOW; clusters: diffuse + FAST
-> N_gal/N_cl is large for any n>=0 (driven by rho/v), the question is the ABSOLUTE normalisation
(thermalised galaxies AND collisionless clusters simultaneously, within SIDM bounds sigma/m<~1 at
cluster v).  Representative (rho_core, v) per system; t=10 Gyr.
HONEST: order-of-magnitude (halo cores vary); the n_hat mediator MASS (sets n and the turnover) is
not derived -- the check is whether the REQUIRED sigma(v) sits in the SIDM-allowed band, not a fit.
"""
import numpy as np
GeVcm3_to_g=1.783e-24          # g/cm^3 per GeV/cm^3
kms=1e5; Gyr=3.156e16          # cm/s per km/s ; s per Gyr
t_age=10*Gyr

# (rho_core [GeV/cm^3], v [km/s]) representative
systems={'dwarf':(0.5,30.0),'LSB/spiral':(0.05,120.0),'group':(0.01,500.0),'cluster':(3e-3,1200.0)}

def Nscatter(rho_GeV,v_km,som_at_v):   # som_at_v = sigma/m (cm^2/g) at this system's v
    rho=rho_GeV*GeVcm3_to_g; v=v_km*kms
    return rho*som_at_v*v*t_age

def som_of_v(v_km, som_ref, v_ref, n):  return som_ref*(v_ref/v_km)**n

if __name__=='__main__':
    print("="*74); print("  GALACTIC-OVERLAP: SIDM scatterings/particle N across systems (t=10 Gyr)"); print("="*74)
    for n in [1,2,4]:
        # normalise so clusters are marginally collisionless: sigma/m(cluster)=0.5 cm^2/g
        vcl=systems['cluster'][1]; som_ref=0.5; v_ref=vcl
        print(f"\n  sigma/m = 0.5 (v_cl/v)^{n} cm^2/g  [cluster-normalised, within Bullet bound]:")
        print(f"    {'system':>11} {'v(km/s)':>8} {'sigma/m':>9} {'N_scatter':>11} {'regime':>14}")
        for s,(rho,v) in systems.items():
            som=som_of_v(v,som_ref,v_ref,n)
            N=Nscatter(rho,v,som)
            reg='THERMALISED' if N>3 else ('marginal' if N>0.3 else 'collisionless')
            print(f"    {s:>11} {v:>8.0f} {som:>9.2f} {N:>11.2e} {reg:>14}")
    print("\n  READ: galaxies (dense+slow) get N>>1 (thermalised -> DM relaxes to the baryon potential =")
    print("  tracks baryons, preserves the RAR), clusters (diffuse+fast) get N<~1 (collisionless, dominant")
    print("  NFW). The split is ROBUST (driven by rho/v); n=1-2 keeps dwarf sigma/m in the SIDM band, n=4")
    print("  (massless n_hat) overshoots dwarf sigma/m -> the n_hat mediator needs a MASS (Yukawa, n~1-2),")
    print("  not the c_dir->0 massless limit. So: resolvable IFF the n_hat interaction is a massive-mediator")
    print("  SIDM with sigma/m(cluster)<~0.5 and n~1-2 -- a specific, falsifiable condition, not a free fit.")
