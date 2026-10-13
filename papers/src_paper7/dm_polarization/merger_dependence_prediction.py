#!/usr/bin/env python3.11
"""
merger_dependence_prediction.py -- SHARPEN Paper VII's merger-dependent cluster-mass signature
(P7:rem:defect_network) from a qualitative "disturbed clusters carry a denser network -> larger
excess" into a QUANTITATIVE, falsifiable M_excess(relaxation-state) curve, and state the
MOND / LambdaCDM discriminators.

FRAMEWORK (P7:rem:defect_network): the fossil disclination network adds mass beyond the
baryon-tracking strain,
    M_excess  ==  M_def/M_strain  ~  (R_cl / xi_net)^2 ,
with xi_net the frozen-network spacing. A merger INJECTS the network at the subcluster scale
xi_inj ~ xi_turb ~ (1/2) R_cl -> at injection M_excess ~ (R_cl/(R_cl/2))^2 = 4 (order of the
observed cluster-to-RAR deficit ~28x/5x ~ 5-6: right ballpark, no tuned normalisation).

COARSENING: Paper VII stresses the network heals by ANNIHILATION (a defect mobility), NOT by the
signal speed c_dir. A line network coarsening by tension-driven annihilation against rotational
viscosity obeys the standard curvature law
    d(xi^2)/dt = 2 M_mob      ->   xi_net(t)^2 = xi_inj^2 + 2 M_mob t ,
so, writing tau_coarsen = xi_inj^2/(2 M_mob) (the e-fold-ish coarsening time),
    M_excess(t) = (R_cl/xi_net)^2 = 4 / (1 + t/tau_coarsen) .
t = time since the last major merger (the cluster's dynamical/relaxation age). SURVIVAL (the
network persists over a cluster age) needs tau_coarsen >~ t_cluster -- the SAME cold/heavy-director
limit the sector already requires; M_mob (hence tau_coarsen) is cutoff-set, NOT derived, so the
SLOPE's magnitude is uncertain while its SIGN and FORM (monotonic decline) are predicted.

PREDICTION (the test): measure the excess beyond the baryon-tracking RAR mass,
    M_excess = (M_lens - M_RAR)/M_RAR ,
versus a relaxation proxy (X-ray morphology, centroid shift w, BCG--X-ray-peak offset, or the
substructure-inferred time since last major merger). Framework: M_excess DECLINES monotonically
with relaxation time as 4/(1+t/tau_coarsen). Below, the curve for several tau_coarsen.
"""
import numpy as np

R_cl_Mpc   = 2.0                  # representative cluster radius (R_200-ish)
xi_inj     = 0.5*R_cl_Mpc         # injection scale ~ 1/2 R_cl (subcluster merger)
M_excess_0 = (R_cl_Mpc/xi_inj)**2 # = 4 at injection

def M_excess(t_Gyr, tau_Gyr):
    return M_excess_0/(1.0 + t_Gyr/tau_Gyr)

if __name__ == '__main__':
    print("="*74)
    print("  MERGER-DEPENDENT CLUSTER MASS EXCESS  (sharpening P7:rem:defect_network)")
    print("="*74)
    print(f"  xi_inj = (1/2)R_cl -> M_excess(t=0, just merged) = (R_cl/xi_inj)^2 = {M_excess_0:.0f}")
    print(f"  (matches the cluster deficit scale ~28x-lensing / ~5x-RAR ~ 5-6, no tuned norm.)")
    print(f"\n  M_excess(t) = {M_excess_0:.0f} / (1 + t/tau_coarsen)   [t = time since last major merger]\n")
    taus = [5.0, 10.0, 20.0]      # tau_coarsen (Gyr); survival needs tau >~ t_cluster
    ts   = [0.0, 1.0, 2.0, 4.0, 8.0]
    hdr = "  t_relax [Gyr] " + "".join(f"tau={tau:>4.0f}Gyr " for tau in taus)
    print(hdr); print("  " + "-"*(len(hdr)-2))
    for t in ts:
        row = f"  {t:>12.0f}  "
        for tau in taus:
            row += f"{M_excess(t,tau):>11.2f} "
        print(row)
    print("\n  (M_excess = mass beyond the ~5x baryon-tracking RAR; 4 -> ~0 as the cluster relaxes.)")

    print("\n  SLOPE near a merger (the headline number):")
    for tau in taus:
        slope = -M_excess_0/tau   # dM_excess/dt at t=0
        print(f"    tau_coarsen={tau:>4.0f} Gyr:  dM_excess/dt|_0 = {slope:+.2f} per Gyr "
              f"(~{100/tau:.0f}% of the peak excess lost per Gyr of relaxation)")

    print("\n" + "="*74)
    print("  DISCRIMINATORS (same baryonic content, varying dynamical state):")
    print("   * THIS FRAMEWORK: M_excess DECLINES monotonically with relaxation time,")
    print("     M_excess = 4/(1+t/tau); disturbed/merging clusters show the LARGEST excess.")
    print("   * MOND (force law slaved to local baryonic g_N): NO merger dependence at all --")
    print("     M_excess fixed by the baryon distribution, independent of dynamical state.")
    print("   * LambdaCDM (collisionless halo): lensing mass = total halo mass, set by the")
    print("     ASSEMBLY/mass, not the relaxation state; disturbed clusters show transient")
    print("     projection/substructure scatter, NOT a systematic excess-vs-relaxation trend")
    print("     (and major mergers can BIAS projected lensing mass UP transiently -- different")
    print("     sign/structure from the framework's smooth monotonic decline).")
    print("  => a clean three-way discriminator in existing lensing + X-ray-morphology samples")
    print("     (relaxed vs disturbed cluster lensing). SLOPE magnitude is cutoff-set (M_mob);")
    print("     SIGN + monotonic form 4/(1+t/tau) are the framework prediction. [O on slope]")
    print("="*74)
