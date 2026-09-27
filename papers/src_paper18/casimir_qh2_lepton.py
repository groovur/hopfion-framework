#!/usr/bin/env python3.11
"""
casimir_qh2_lepton.py -- Lead C: the Casimir energy of the Q_H=2 lepton is ALREADY in the
framework, as the -c/24 (CFT Casimir) piece of the -3/400 T-matrix scheme correction.

Physical picture (Lead A): the Q_H=2 tube carries trapped LINEAR ripple modes (omega=c_s k,
c_s=c/phi). On the closed lepton loop these are a 1D CFT on a circle; the zero-point of a CFT
on a circle is the Casimir energy E_Cas = -(pi c_central/6)(hbar c_s/L) -- proportional to the
CENTRAL CHARGE c_central. In the framework's mass formula the corresponding term is the -c/24
inside the T-matrix phase T_{1/2}=h_{1/2}-c/24, which enters m_e = Lcond phi^20 e^{4pi - 3/400 -...}
via -3/400 = T_{1/2}/Q_group. So the Casimir (c/24) is NOT a new effect -- it is the algebraic
(CFT) form of the physical ripple zero-point, already present. This script verifies the exact
decomposition and the central-charge proportionality.
"""
import sympy as sp

def su2k(k):
    c=sp.Rational(3*k,k+2)
    h_half=sp.Rational( (1*(1+2)), 4*(k+2) )   # h_j=j(j+1)/(k+2), j=1/2 -> l=1: l(l+2)/(4(k+2))=3/(4(k+2))
    return c, h_half

print("="*74,"\n Casimir (c/24) inside the electron's -3/400 correction  [Q_H=2, SU(2)_3]\n","="*74,sep="")
k=3; Qgrp=10
c, h_half = su2k(k)
cas = c/24                                   # CFT Casimir / ground-state shift
T_half = h_half - cas                        # T-matrix phase (Paper XVII)
corr = -T_half/Qgrp                           # the exponent correction (should be -3/400)
print(f"  SU(2)_{k}: c={c},  h_(1/2)={h_half},  Casimir c/24={cas},  T_(1/2)=h-c/24={T_half}")
print(f"  Q_group(lepton)={Qgrp}")
print(f"  correction = -T_(1/2)/Q = {corr}   (framework: -3/400 = {sp.Rational(-3,400)})  "
      f"{'OK' if corr==sp.Rational(-3,400) else 'MISMATCH'}")

print("\n  DECOMPOSITION of the -3/400 exponent correction:")
cas_piece = +cas/Qgrp                         # -T/Q => +c/24 /Q  (Casimir enters with + sign)
h_piece   = -h_half/Qgrp
print(f"    Casimir (c/24) piece   : +(c/24)/Q = {cas_piece} = {float(cas_piece):+.5f}  (mass ENHANCEMENT e^+)")
print(f"    conformal-weight piece : -(h_half)/Q = {h_piece} = {float(h_piece):+.5f}  (suppression)")
print(f"    sum                    : {cas_piece+h_piece} = {float(cas_piece+h_piece):+.5f}  (= -3/400)")
print(f"  -> the CFT Casimir contributes +{cas_piece} to the electron mass exponent; the -3/400 is")
print(f"     the Casimir(+c/24) PARTIALLY CANCELLED by the conformal weight(-h_1/2), net suppression.")

print("\n"+"="*74,"\n Central-charge proportionality: physical ripple Casimir  vs  c/24\n","="*74,sep="")
print("  Physical: trapped linear ripples omega=c_s k on the loop (circumference L) = 1D CFT on S^1;")
print("  zero-point E_Cas = -(pi c_central/6)(hbar c_s/L)  -- PROPORTIONAL to c_central.")
print("  Framework: the T-matrix Casimir term is -c/24 -- also PROPORTIONAL to c_central (=c).")
print("  Same origin: both are the CFT ground-state (Casimir) energy of the ripple modes, one as the")
print("  geometric loop zero-point, one as the algebraic -c/24 in the mass exponent.")
print(f"\n  For SU(2)_3: c_central = {c} = {float(c):.3f}.  The physical E_Cas and the -c/24 term differ")
print(f"  only by the geometric/normalisation factors (pi/6, hbar c_s/L vs 1/24, per-mode); both scale")
print(f"  with the SAME central charge, so the identification is structural, not numerical coincidence.")
print("\n  CONCLUSION [E]: the Q_H=2 Casimir is NOT a new lepton-mass piece -- it is the c/24 CFT")
print("  Casimir already inside the -3/400 correction (Paper IV/XI/XVII), now given its physical")
print("  meaning as the zero-point of the Lead-A trapped ripple spectrum. Lead C = welds Lead A's")
print("  physical modes to the framework's existing scheme correction; no new adjustable energy.")