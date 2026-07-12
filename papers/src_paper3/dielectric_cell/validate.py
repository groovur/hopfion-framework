import time
import numpy as np
import common as c

def report(label, alpha, expected=None, info=None):
    print(f"{label}: alpha = {alpha}")
    if expected is not None:
        print(f"   expected = {expected}, |diff| = {abs(alpha-expected):.3e}")
    if info is not None:
        print(f"   nfev={info.get('nfev')}")

# V1: uniform medium beta=0 -> eps_rad=eps_tan=1
t0 = time.time()
er, et = c.bg_eps_fns(0.0)
alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=1e-4, rho_max=200.0)
report("V1 (uniform, beta=0)", alpha, expected=0.0, info=info)
print(f"  time={time.time()-t0:.3f}s")

# V2: isotropic step inclusion eps1=3, R=1
t0 = time.time()
er, et = c.step_eps_fns(3.0, 1.0)
alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=1e-4, rho_max=200.0)
expected = 1.0**2 * (3.0-1.0)/(3.0+1.0)
report("V2 (step eps1=3, R=1)", alpha, expected=expected, info=info)
print(f"  time={time.time()-t0:.3f}s")

# V3: step inclusion eps1 = -3 + i*eta, R=1
t0 = time.time()
eta3 = 1e-3
eps1 = -3.0 + 1j*eta3
er, et = c.step_eps_fns(eps1, 1.0)
alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=1e-4, rho_max=200.0)
expected = 1.0**2 * (eps1-1.0)/(eps1+1.0)
report("V3 (step eps1=-3+i*1e-3, R=1)", alpha, expected=expected, info=info)
print(f"  time={time.time()-t0:.3f}s")

# sensitivity checks: rho_min and rho_max variation for V2
print("\n--- V2 sensitivity ---")
for rmin in [1e-3, 1e-4, 1e-5]:
    er, et = c.step_eps_fns(3.0, 1.0)
    alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=rmin, rho_max=200.0)
    print(f"  rho_min={rmin}: alpha={alpha}")
for rmax in [50.0, 100.0, 200.0, 400.0]:
    er, et = c.step_eps_fns(3.0, 1.0)
    alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=1e-4, rho_max=rmax)
    print(f"  rho_max={rmax}: alpha={alpha}")

# tighter tolerance check (resolution / accuracy doubling analog for adaptive RK45)
print("\n--- V2 tolerance convergence ---")
for rt in [1e-8, 1e-10, 1e-12]:
    er, et = c.step_eps_fns(3.0, 1.0)
    alpha, info = c.shoot_alpha(er, et, eta=0.0, rho_min=1e-4, rho_max=200.0, rtol=rt, atol=rt*1e-2)
    print(f"  rtol={rt}: alpha={alpha}  nfev={info['nfev']}")
