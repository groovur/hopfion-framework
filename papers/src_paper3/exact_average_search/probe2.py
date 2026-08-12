import sympy as sp

s, c = sp.symbols('s c', positive=True)

def normalize_atan(expr):
    expr = expr.subs(sp.atan(1/c), sp.pi/2 - sp.atan(c))
    expr = sp.simplify(expr)
    return expr

expr = 1/(s**2+c**2)**2
I = sp.integrate(expr, (s, 1, sp.oo))
I = normalize_atan(I)
print("test1:", I)

expr2 = s**4/(s**2+c**2)**3
I2 = sp.integrate(expr2, (s, 1, sp.oo))
I2 = normalize_atan(I2)
print("test2:", I2)

# Now try to write as A0(c) + A1(c)*atan(c) form using collect
A = sp.Symbol('A')
def to_R0_R1(expr):
    e = expr.subs(sp.atan(c), A)
    e = sp.together(e)
    poly = sp.Poly(sp.numer(e), A) if e.is_rational_function(A) else None
    return e

print(to_R0_R1(I2))
