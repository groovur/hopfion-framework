import sympy as sp

s, c = sp.symbols('s c', positive=True)

# test basic integral
expr = 1/(s**2+c**2)**2
I = sp.integrate(expr, (s, 1, sp.oo))
print("test1:", sp.simplify(I))

expr2 = s**4/(s**2+c**2)**3
I2 = sp.integrate(expr2, (s, 1, sp.oo))
print("test2:", sp.simplify(I2))
