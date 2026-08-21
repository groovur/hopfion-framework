"""Which director mode is the unstable one? Isolate the GEOMETRIC labelling
(splay vs bend for a given modulation direction) from the band numbers.

Established:
  - e_r (director n_hat) = the LIGHT axis of the semi-Dirac dispersion
    (S_eff derivation: radial = light = 'easy'; heavy = perpendicular).
  - frank_run (gate-validated): rotation texture modulated along the HEAVY axis
    has K < 0 (unstable); along the LIGHT axis K > 0 (stable).

Unknown to pin: is 'modulation along the heavy axis' a SPLAY or a BEND of e_r?
Frank definitions (unambiguous, real space):
    splay = div n_hat ;  bend = | n_hat x (curl n_hat) |.
Put e_r along y (light). Heavy axis = x. Rotate the director by a small
theta(x,y): n_hat = (sin theta, cos theta) ~ (theta, 1).
Test modulation along x (heavy) vs along y (light).
"""
import sympy as sp

x, y, qh, ql, th0 = sp.symbols('x y q_h q_l theta0', positive=True)

def modes(theta):
    # n_hat in the xy-plane, small-theta about e_r = y_hat
    nx, ny = sp.sin(theta), sp.cos(theta)
    div = sp.diff(nx, x) + sp.diff(ny, y)                 # splay = div n
    curl_z = sp.diff(ny, x) - sp.diff(nx, y)              # (curl n)_z
    # n x (curl n): n=(nx,ny,0), curl=(0,0,curl_z) -> (ny*curl_z, -nx*curl_z, 0)
    bend = sp.sqrt(sp.simplify((ny*curl_z)**2 + (nx*curl_z)**2))
    return sp.simplify(div), sp.simplify(bend)

print("=== modulation along HEAVY axis (x):  theta = theta0 cos(q_h x) ===")
div_h, bend_h = modes(th0*sp.cos(qh*x))
print("  splay = div n_hat =", sp.simplify(sp.series(div_h, th0, 0, 2).removeO()))
print("  bend  =", sp.simplify(sp.series(bend_h, th0, 0, 2).removeO()))

print("\n=== modulation along LIGHT axis (y):  theta = theta0 cos(q_l y) ===")
div_l, bend_l = modes(th0*sp.cos(ql*y))
print("  splay = div n_hat =", sp.simplify(sp.series(div_l, th0, 0, 2).removeO()))
print("  bend  =", sp.simplify(sp.series(bend_l, th0, 0, 2).removeO()))

print("""
READING:
- modulation along HEAVY (x, perpendicular to e_r): splay != 0, bend = 0  => SPLAY
- modulation along LIGHT (y, along e_r):            splay = 0, bend != 0  => BEND
(Standard Frank: splay = gradient PERPENDICULAR to n_hat; bend = gradient ALONG n_hat.)

Combine with frank_run:
  K(modulation along HEAVY) = -166 < 0  ==> SPLAY is UNSTABLE
  K(modulation along LIGHT) = +2.5 > 0  ==> BEND  is STABLE

=> The 'splay instability' label is CORRECT. The only error is the paper's
   direction wording ('heavy along e_r' -> should be heavy PERPENDICULAR to e_r);
   it does not change which mode is unstable.
""")
