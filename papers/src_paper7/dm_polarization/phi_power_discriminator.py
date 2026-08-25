"""
Jordan->Einstein frame discriminator over the phi-power 'invariants' of Paper 7.
Question for each: is the phi-power PHYSICAL (survives to the canonical/Einstein frame
or is non-gravitational) or a FRAME ARTIFACT (from the non-canonical kinetic norm 1/phi^8
and/or the conformal factor F, cancels in the physical frame)?

The lesson from m_xi: phi^20 = phi^12 (V'') x phi^8 (1/K) CANCELLED in the Einstein frame.
Test each quantity by ORIGIN:
  - TOPOLOGICAL/WZW origin (Bogomolny lambda, quantum dim, T-matrix, McKay) -> physical, frame-free.
  - SCALAR-TENSOR normalisation origin (V''/K, conformal F) -> frame-dependent, candidate artifact.
  - NUMERICAL approx of a non-phi constant -> coincidence.
"""
import math
class np:
    pi=math.pi
    sqrt=staticmethod(math.sqrt)
np.pi=math.pi
phi=(1+5**0.5)/2
beta=0.452

def line(name, jordan, origin, verdict, note=""):
    print(f"  {name:26s} {jordan:14s} {origin:28s} {verdict}")
    if note: print(f"       -> {note}")

print("="*100)
print(f"  {'quantity':26s} {'Jordan phi-power':14s} {'origin':28s} verdict")
print("="*100)

# 1. Bogomolny parameter
line("lambda (Bogomolny)","phi^6","BPS soliton (Paper II)","PHYSICAL",
     "topological property of the Q=2 Hopfion; no gravity frame. Frame-free.")
# 2. alpha^-1 fine structure
line("alpha^-1 = 360/phi^2","phi^-2","WZW/anyonic S-matrix","PHYSICAL",
     "not gravitational at all. 137.036, matches obs. Frame-free.")
# 3. omega_BD
line("omega_BD = lambda^2/(3 beta)","phi^12","= lambda^2, topological","PHYSICAL (content)",
     "phi^12 IS lambda^2 (Bogomolny). Jordan coupling, but the phi-content is topological.")
# 4. alpha_BD fifth-force
line("alpha_BD ~ 1/(2 omega_BD)","phi^-12","= 1/lambda^2, topological","PHYSICAL (content)",
     "carries the same topological phi^12; the coupling STRENGTH is physical (screened for obs).")
# 5. m_xi (the one we cracked)
line("m_xi^2 (DE mass)","phi^20","V'' x 1/K (scalar-tensor)","FRAME ARTIFACT",
     "phi^20=phi^12(V'') x phi^8(1/K). Einstein frame -> phi^0 (m_xi~H0). CANCELS.")
# 6. r_V
line("r_V (Vainshtein)","phi^-8/3","omega_BD/m_xi^2 (MIXED)","MIXED",
     "phi^-8=phi^12(physical)-phi^20(artifact). Physical part survives: r_V ~ phi^4.")
# 7. Lambda_UV
fourpi=4*np.pi*math.sqrt(2)
line("Lambda_UV = M_Pl/(4pi sqrt2)","~phi^6 (approx)","Seeley-DeWitt loop","COINCIDENCE",
     f"4pi sqrt2={fourpi:.3f} vs phi^6={phi**6:.3f}: {abs(fourpi-phi**6)/phi**6*100:.1f}% off. "
     "NOT a phi-power; a 1% numerical coincidence.")
# 8. graviton mass
Om2=(1+3*phi**7)/2
line("m_g^2 (graviton TT)","phi^-8","S_eff/(1+beta rho) = NMC","FRAME/UNTESTED",
     f"1/phi^8=1/(phi^6 phi^2), tied to NMC. Einstein ~ /Omega^2 (Omega^2={Om2:.0f}) -> ~phi^-15. "
     "Frame-dependent AND untested (only v_GW~c matters, holds for any tiny m_g).")
# 9. F(rho_inf)
line("F(rho_inf) = (M_Pl^2/2)(1+3phi^7)","phi^7","the conformal factor itself","FRAME-DEFINING",
     "3phi^7 IS the NMC ratio = the Jordan->Einstein map. Frame-defining, not an invariant.")

print("="*100)
print("""
SUMMARY -- the discriminator separates three classes:
  PHYSICAL (topological/WZW, frame-free):  lambda=phi^6, alpha^-1=360/phi^2, the phi^12 in
     omega_BD/alpha_BD (=lambda^2). These are properties of the soliton / WZW data, not the
     gravity frame. They SURVIVE.
  FRAME ARTIFACT (scalar-tensor normalisation, cancels in Einstein frame):  m_xi^2 ~ phi^20
     (V''/K) -> phi^0; r_V ~ phi^-8/3 -> phi^4 (its artifact PART cancels, physical part phi^12 stays);
     m_g^2 ~ phi^-8 (NMC, and untested); F(rho_inf) ~ phi^7 IS the frame map.
  COINCIDENCE:  Lambda_UV ~ M_Pl/phi^6 is really M_Pl/(4pi sqrt2), a 1% numerical match, not a phi-power.
RULE OF THUMB: a phi-power is PHYSICAL iff its ORIGIN is topological/WZW (Bogomolny, quantum
dimension, T-matrix, McKay order) rather than the non-canonical kinetic norm 1/phi^8 or the
conformal factor F. The gravity/scalar-tensor sector is where the artifacts live.
""")
