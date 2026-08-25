import math
phi=(1+5**0.5)/2; pi=math.pi
Ne=57.73
# P_s = (Ne^2/12pi^2)(Lam_UV/M_Pl)^8, needs P_s=2.10e-9 (matches Planck)
def Ps(x): return (Ne**2/(12*pi**2))*x**8*(18/19)**4
for lab,x in [("Lam_UV/M_Pl = 1/phi^6 (sub-Planckian, paper)", 1/phi**6),
              ("Lam_UV/M_Pl = 1/(4pi sqrt2) (exact induced, sub)", 1/(4*pi*math.sqrt(2))),
              ("Lam_UV/M_Pl = 4pi sqrt6 (SUPER-Planckian, 'corrected')", 4*pi*math.sqrt(6))]:
    print(f"  {lab:52s}: x={x:.4g}, P_s={Ps(x):.2e}")
print()
print(f"  Observed P_s ~ 2.1e-9.  Paper's 1/phi^6 gives {Ps(1/phi**6):.2e} (MATCHES).")
print(f"  Super-Planckian gives {Ps(4*pi*math.sqrt(6)):.1e} -- off by ~{math.log10(Ps(4*pi*math.sqrt(6))/2.1e-9):.0f} ORDERS.")
print()
print("=> P_s (hence n_s, r, the whole CMB normalisation, load-bearing, matches Planck)")
print("   DEPENDS on Lam_UV/M_Pl = 1/phi^6 (SUB-Planckian). My earlier 'nothing load-bearing")
print("   depends on Lam_UV' was WRONG: the CC/DE sector cancels Lam_UV, but INFLATION does NOT.")
print("   'Correcting' Lam_UV to super-Planckian would DESTROY the inflation predictions.")
