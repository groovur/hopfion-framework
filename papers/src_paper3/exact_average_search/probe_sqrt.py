import sympy as sp
s, c = sp.symbols('s c', positive=True)
eps_tan = s**4/(s**2+c**2)**2
eps_rad = (s**2-3*c**2)*s**4/(s**2+c**2)**3
integrand = sp.sqrt(eps_tan*eps_rad)
integrand = sp.simplify(integrand)
print("integrand:", integrand)
# domain where real: s^2 - 3c^2 >=0 i.e. s >= sqrt(3)*c
try:
    I = sp.integrate(integrand, (s, sp.sqrt(3)*c, sp.oo))
    print("integral:", I)
except Exception as e:
    print("failed:", e)
