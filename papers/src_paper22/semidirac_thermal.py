#!/usr/bin/env python3.11
"""
semidirac_thermal.py -- thermal free energy of the EPOCH's GEOMETRY (user's key point):
the QCD-epoch condensate excitations have the SEMI-DIRAC dispersion (claude-hopfion s3),
NOT a relativistic one. This changes the Polyakov-loop (Weiss/GPY) potential -> the b2,b3,b4
of the confinement Landau potential. Toward pinning eta (notes (mm)).

Semi-Dirac:  E(k,theta)^2 = A^2 k^4 sin^4(theta) + v^2 k^2 cos^2(theta)
  easy axis (theta=0): E=v k (linear);  hard axis (theta=pi/2): E=A k^2 (quadratic).
Thermal free energy density:  f(T) = T INT d^3k/(2pi)^3 ln(1 - e^{-E/T}).
HOLONOMY (Polyakov phase phi, the color/confinement order parameter L~e^{i phi}):
  f(T,phi) = T INT d^3k/(2pi)^3 Re ln(1 - e^{-E/T} e^{i phi})
           = -T INT d^3k/(2pi)^3 sum_{n>=1} (1/n) e^{-nE/T} cos(n phi).
EXPECTED (scaling argument): INT d^3k e^{-nE/T} ~ (T/n)^2/(vA) for semi-Dirac (vs (T/n)^3
relativistic), so f(phi) ~ (T^3/vA) sum cos(n phi)/n^3 = a CUBIC Bernoulli B_3(phi) (vs the
relativistic QUARTIC B_4). The B3 structure directly feeds the cubic b3 of the Z3 potential.
This script CONFIRMS the T-power and the B3 structure numerically. Units: A=v=1 (dimensionless;
the shape, not the scale, is what sets b3/b4). The SU(3) 3-color sum + gap(L) -> exact b3/b4
is the remaining step (flagged).
"""
import numpy as np
from scipy import integrate

def E_semidirac(k, th, A=1.0, v=1.0):
    return np.sqrt(A**2*k**4*np.sin(th)**4 + v**2*k**2*np.cos(th)**2)

def f_thermal(T, A=1.0, v=1.0, kmax=60.0):
    """f(T) = T/(2pi)^3 * 2pi * INT k^2 dk sin th dth ln(1-e^{-E/T})."""
    def integrand(th, k):
        E = E_semidirac(k, th, A, v)
        return k**2*np.sin(th)*np.log1p(-np.exp(-E/T))
    val,_ = integrate.dblquad(integrand, 1e-6, kmax, lambda k:1e-6, lambda k:np.pi-1e-6,
                              epsabs=1e-10, epsrel=1e-8)
    return T*(2*np.pi)*val/(2*np.pi)**3

def f_holonomy_fourier(T, nmax=40, A=1.0, v=1.0, kmax=60.0):
    """coefficients c_n = T/(2pi)^3 2pi INT k^2 dk sin th dth (1/n) e^{-nE/T}, so
    f(phi) = - sum_n c_n cos(n phi). Return c_n (n=1..nmax) and the implied B-structure."""
    cs=[]
    for n in range(1,nmax+1):
        def integ(th,k):
            E=E_semidirac(k,th,A,v)
            return k**2*np.sin(th)*np.exp(-n*E/T)/n
        val,_=integrate.dblquad(integ,1e-6,kmax,lambda k:1e-6,lambda k:np.pi-1e-6,
                                epsabs=1e-12,epsrel=1e-8)
        cs.append(T*(2*np.pi)*val/(2*np.pi)**3)
    return np.array(cs)

if __name__=='__main__':
    print("="*74)
    print("  SEMI-DIRAC THERMAL FREE ENERGY (epoch geometry)  E^2=A^2k^4 sin^4 th + v^2k^2 cos^2 th")
    print("="*74)
    # (1) T-power: f(T) ~ T^p ; semi-Dirac expect p=3 (vs relativistic p=4)
    Ts=[0.5,1.0,2.0,4.0]
    fs=[f_thermal(T) for T in Ts]
    print("  T-scaling  f(T) ~ T^p :")
    for T,f in zip(Ts,fs): print(f"    T={T:<4} f={f:>12.5e}")
    import numpy as _np
    p=_np.polyfit(_np.log(Ts),_np.log(_np.abs(fs)),1)[0]
    print(f"    => fitted power p = {p:.3f}  (semi-Dirac expect ~3; relativistic would be 4)")
    # (2) holonomy Fourier coefficients c_n -> is it 1/n^3 (B3, cubic) not 1/n^4 (B4, quartic)?
    print("\n  holonomy Fourier coeffs c_n (f(phi)=-sum c_n cos n phi) at T=2:")
    cs=f_holonomy_fourier(2.0,nmax=8)
    for n,c in enumerate(cs,1):
        print(f"    n={n}  c_n={c:>12.5e}   c_n*n^3={c*n**3:>10.4e}   c_n*n^4={c*n**4:>10.4e}")
    # which power makes c_n*n^p ~ const?
    ns=_np.arange(1,len(cs)+1)
    q=-_np.polyfit(_np.log(ns),_np.log(_np.abs(cs)),1)[0]
    print(f"    => c_n ~ 1/n^q with q = {q:.3f}  (B3/cubic if q~3; B4/quartic if q~4)")
    print("\n  READ: p~3 and q~3 => semi-Dirac gives a T^3, CUBIC-Bernoulli Weiss potential")
    print("  (vs relativistic T^4/B4) -> a framework-specific b2,b3,b4. Next: SU(3) 3-color sum")
    print("  + gap(L) -> exact b3/b4 -> h_c, first-order strength, eta.")
