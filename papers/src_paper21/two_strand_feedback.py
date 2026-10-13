#!/usr/bin/env python3.11
"""
*** STATUS: FAILED numerically (2026-10-13). KEPT as a record, DO NOT trust its numbers. ***
The naive DOUBLET SUPERPOSITION (z1_A+z1_B, z2_A+z2_B) creates artificial SINGULARITIES where the
summed doublet vanishes -> n_hat ill-defined -> spurious g2 ~ 2700 (vs a real single-tube ~O(C*^2)),
and a flat/absurd ~1e5x "enhancement" that is an artefact, not self-interaction. The baseline tube
(d=10) also sits OUTSIDE the box. So this does NOT test the barrier idea. A proper two-strand field
needs a singularity-free multi-soliton construction (multiply rational maps, not add doublets) --
which is the SAME hard construction as the knotted field. Lesson re-logged (notes (nnnn)).
-------------------------------------------------------------------------------------------------
two_strand_feedback.py -- PROOF-OF-CONCEPT of the user's idea (2026-10-13): the feedback SATURATION
beta_*rho sits at the pion emission/absorption THRESHOLD = the n=4 reconnection BARRIER (Paper XVIII
P18:rem:wound_state), where two tube-strands are mid-threading and MAXIMALLY overlap. This captures the
SELF-INTERACTION the 1D single-tube profile discarded (notes (mmmm) retraction), and is LOCAL +
topologically pinned, not a global knot relaxation.

MODEL (local, honest): two straight parallel Hopf tubes along x_hat, at (y,z)=(0,+-d/2), each the
standard Hopf-tube doublet  z1=cos(f0/2), z2=sin(f0/2) e^{i m psi}  (rho,psi polar about its centre;
f0(rho)=2 arctan((rho C*)^{-1}), the BS profile; m = poloidal winding). Fields are SUPERPOSED (summed
doublets, renormalised -- the same combination build_doublet uses). The configuration is x-translation
invariant, so g2 = |d_y n|^2 + |d_z n|^2 on the (y,z) plane. We scan the strand separation d and read
the PEAK g2 and its ENHANCEMENT over a single isolated tube. If the overlap at the reconnection
separation (d ~ the tube core, where strands must meet to re-link) boosts g2 by a large factor, the
self-interaction CAN carry beta_*rho from the single-tube ~O(1) up to the saturated regime -- the
quantitative test of whether the barrier picture closes the eta magnitude.

CAVEATS (do not over-read): (i) two STRAIGHT PARALLEL tubes is a local model of the crossing, not the
exact n=4 half-threaded cable; (ii) the physical barrier separation d_barrier is set by the (2,3)-cable
framing (not pinned here -- we scan d); (iii) beta normalisation (self-consistent vs reference) still
enters the absolute beta_*rho. So this tests the MECHANISM + the enhancement factor, not a final number.
"""
import numpy as np

C_star=2.5062
def f0(rho): return 2*np.arctan((np.maximum(rho,1e-9)*C_star)**(-1.0))

def tube_doublet(Y,Z,y0,z0,m):
    rho=np.sqrt((Y-y0)**2+(Z-z0)**2); psi=np.arctan2(Z-z0,Y-y0)
    f=f0(rho)
    return np.cos(f/2).astype(complex), np.sin(f/2)*np.exp(1j*m*psi)

def nhat(z1,z2):
    n1=2*np.real(np.conj(z1)*z2); n2=2*np.imag(np.conj(z1)*z2); n3=np.abs(z1)**2-np.abs(z2)**2
    return np.stack([n1,n2,n3])

def g2_peak(d, m, L=3.0, N=240):
    ax=np.linspace(-L,L,N); h=ax[1]-ax[0]
    Y,Z=np.meshgrid(ax,ax,indexing='ij')
    za=tube_doublet(Y,Z,0,+d/2,m); zb=tube_doublet(Y,Z,0,-d/2,m)
    z1=za[0]+zb[0]; z2=za[1]+zb[1]
    nrm=np.sqrt(np.abs(z1)**2+np.abs(z2)**2)+1e-30; z1/=nrm; z2/=nrm
    n=nhat(z1,z2)
    dny=np.stack([np.gradient(n[c],h,axis=0) for c in range(3)])
    dnz=np.stack([np.gradient(n[c],h,axis=1) for c in range(3)])
    g2=(dny**2).sum(0)+(dnz**2).sum(0)
    return g2.max()

if __name__=='__main__':
    print("="*74)
    print("  TWO-STRAND CROSSING FEEDBACK (barrier picture; captures self-interaction)")
    print(f"  Hopf tubes, C*={C_star}, core radius ~1/C*={1/C_star:.2f}; scan separation d")
    print("="*74)
    m=1
    single=g2_peak(10.0,m)        # isolated tube (d huge)
    print(f"  single isolated tube: g2_peak = {single:.2f}  (baseline; = C*^2-ish)")
    print(f"  {'d':>6} {'d/core':>7} {'g2_peak':>9} {'enhancement':>12}")
    for d in [3.0,2.0,1.5,1.0,0.8,0.6,0.5,0.4,0.3,0.2]:
        g=g2_peak(d,m); print(f"  {d:>6.2f} {d*C_star:>7.2f} {g:>9.2f} {g/single:>11.2f}x")
    print()
    print("  READ: enhancement = g2_peak(two strands)/g2_peak(single). If it grows large as d -> core")
    print("  (strands overlapping at the reconnection barrier), the SELF-INTERACTION boosts the feedback")
    print("  g2 -> beta_*rho rises from the single-tube ~O(1) toward the saturated regime. That is the")
    print("  user's barrier picture: saturation at the emission/absorption threshold. The factor needed")
    print("  to reach beta_*rho~15 from the single-tube ~0.5 is ~30x. NOTE: straight-parallel-tube model")
    print("  + d_barrier not pinned + beta normalisation open -> this tests the MECHANISM, not a number.")
    print("="*74)
