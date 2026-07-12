import numpy as np
from fractions import Fraction as Fr

k = 3
kp2 = k + 2                      # 5
phi = (1 + 5**0.5)/2
golden_angle = 360.0/phi**2      # = 360*(2-phi)

spins = [Fr(0), Fr(1,2), Fr(1), Fr(3,2)]

def h(j):                        # conformal weight / topological spin exponent
    return Fr(j.numerator*(j.numerator+j.denominator), (j.denominator**2)*kp2) if False else j*(j+1)/kp2

def qdim(j):                     # quantum dimension
    return np.sin((2*float(j)+1)*np.pi/kp2)/np.sin(np.pi/kp2)

print("== VALIDATION against known SU(2)_3 / Fibonacci-anyon data ==")
for j in spins:
    hj = j*(j+1)/kp2
    ang = float(hj)*360.0
    print(f"  j={str(j):>3}  h_j={str(hj):>6} (rational)  twist angle 360 h_j = {ang:8.3f} deg   d_j={qdim(j):.6f}")
print(f"  -> spin-1 twist should be Fibonacci tau = e^(4 pi i/5) = 144 deg : {float(Fr(2,5))*360:.1f} deg")
print(f"  -> d_(1/2)=d_1 should equal phi={phi:.6f}")

print("\n== MONODROMY (R^2) phases for fusing two spin-1/2 quanta ==")
a = b = Fr(1,2)
ha = a*(a+1)/kp2
for c in [Fr(0), Fr(1)]:
    hc = c*(c+1)/kp2
    exponent = hc - ha - ha            # rational
    ang = float(exponent % 1)*360.0
    print(f"  channel c={str(c):>3}: monodromy exp = {str(exponent):>7} * 2pi  ->  {ang:8.3f} deg   (rational multiple of 360)")

print("\n== TARGET ==")
print(f"  golden angle 360/phi^2 = 360*(2-phi) = {golden_angle:.6f} deg")
print(f"  is golden angle a rational multiple of 360?  360*(2-phi) with (2-phi)={2-phi:.6f} irrational -> NO")

print("\n== DIRECT COMPARISON: does any modular-data phase equal the golden angle? ==")
phases = {"twist j=0":0.0,"twist j=1/2":54.0,"twist j=1":144.0,"twist j=3/2":270.0,
          "monodromy c=0":float(Fr(-6,20)%1)*360,"monodromy c=1":float(Fr(2,20)%1)*360}
for name,ang in phases.items():
    print(f"  {name:>16}: {ang:8.3f} deg   | diff from golden {abs(ang-golden_angle):7.3f}")
