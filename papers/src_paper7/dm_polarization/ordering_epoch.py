import math
kB=8.617e-5  # eV/K
Tcmb0=2.7255  # K
Tcmb0_eV=Tcmb0*kB           # 2.35e-4 eV
Lam_cond=Tcmb0_eV*(math.pi**2/150)**0.25   # condensate scale (= T_CMB(pi^2/150)^{1/4})
Trec=0.26   # eV, recombination (z~1100, T~3000 K)
zrec=1100
print(f"Tcmb(0) = {Tcmb0_eV:.3e} eV ;  Lam_cond = {Lam_cond:.3e} eV  (= condensate scale, from T_CMB)")
print(f"T_rec ~ {Trec} eV  (z~{zrec})")
print()
print("== nematic ordering temperature T_ord ~ (band factor) x Lam_cond ==")
print("  T_ord ~ energy per coherence cell = rho_cond * xi_cond^3.")
print("  rho_cond = 10 Lam_cond^4 (natural), xi_cond = 1/Lam_cond => rho_cond xi_cond^3 = 10 Lam_cond.")
print("  With the O(0.1-1) band stiffness K_band and (pi/2) factor: T_ord ~ (0.1 - 10) Lam_cond.")
print()
print(f"  {'band factor':>12} {'T_ord [eV]':>12} {'z_order (T(z)=T_ord)':>22} {'vs recomb':>12}")
for f in (0.1,1.0,10.0):
    Tord=f*Lam_cond
    zord=Tord/Tcmb0_eV-1
    verdict="DISORDERED at recomb" if Tord<Trec else "ordered at recomb"
    print(f"  {f:>12.1f} {Tord:>12.3e} {zord:>22.2f} {verdict:>22}")
print()
print(f"  For ANY band factor, T_ord < {Trec} eV = T_rec by {Trec/(10*Lam_cond):.0f}-{Trec/(0.1*Lam_cond):.0f}x.")
print("  => the medium is thoroughly DISORDERED at recombination (z~1100) in every case.")
print()
print("== empirical anchor: galactic DM works out to z~few (observed RAR/flat curves) ==")
for z in (0,2,4):
    print(f"    T(z={z}) = {Tcmb0_eV*(1+z):.3e} eV -> needs T_ord > this for medium ordered at z={z}")
print("  => T_ord ~ 10 Lam_cond ~ 1.2e-3 eV (z_order~4) is the empirically-favoured value:")
print("     medium orders at z~4 (galaxy-formation epoch), giving late-time galactic DM, disordered earlier.")
print()
print("== CLOSURE ==")
print("  T_ord ~ Lam_cond scale ~ 1e-4 to 1e-3 eV, set by the CONDENSATE scale (from T_CMB), NOT the UV cutoff.")
print("  This is ~200-2600x BELOW T_rec ~ 0.26 eV. So the nematic medium is DISORDERED at recombination and")
print("  orders only at z ~ O(1)-few (when T cools to ~Lam_cond). ROBUST to the O(1-10) band factor.")
print("  => the adiabatic baryon-sourced strain is ABSENT at recombination; the director-strain DM is a")
print("     LATE-TIME (z<~4) galactic component. It does NOT supply recombination-era CDM, so the standard")
print("     147 Mpc is NOT recovered by this sector -- the recombination CDM (Omega_DM~0.26, already an input)")
print("     is a separate/gap, not the same object.")
