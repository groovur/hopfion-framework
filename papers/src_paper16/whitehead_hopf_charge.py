"""
whitehead_hopf_charge.py
========================
Compute the Hopf charge Q_H of the global S3-lift ansatz via the
Whitehead / Berry-Chern-Simons integral:

    Q_H = (1 / 4*pi^2) * int A.(curl A) d^3x

where A_i = Im(Z† d_i Z) is the Berry connection of the normalised
Hopf doublet Z = (z1, z2), built by the global two-strand S3-lift
(Construction constr:global_s3lift, Paper XVI sec:global_nondegen).

Key facts established before running this script:
  * z2=0 on the ENTIRE midline curve M(t)=(G(t)+G(t+pi))/2, a closed
    vortex loop in the z=0 plane (not just at the three crossings).
  * The single-strand ansatz has A.(curl A) = 0 identically (A is a
    pure gradient => curl A = grad(sin^2) x grad(3t), then A.curl A = 0).
  * The S3-lift is smooth everywhere (|Z|_min >= 0.386, proved), so the
    Berry connection is bounded and the integral converges.

Validation: the inverse-stereographic Q_H=1 Hopf doublet
    Z1 = (1-r^2 + 2iz) / (1+r^2),  Z2 = 2(x+iy) / (1+r^2)
is smooth, |Z|=1 everywhere, and gives Q_H converging to 1.0 at
increasing N (0.72 at N=20, 0.86 at N=30, 0.92 at N=40, 0.95 at N=56).

Usage:
    python whitehead_hopf_charge.py               # coarse: N=40, ~30s
    python whitehead_hopf_charge.py --N 80        # medium: N=80, ~5min
    python whitehead_hopf_charge.py --N 120       # fine:  N=120, ~20min

F. Manfredi / verification, June 2026
"""

import sys, time
import numpy as np
from scipy.spatial import KDTree

# ── parse N from command line ──────────────────────────────────────────────
N = 40
for i, arg in enumerate(sys.argv[1:]):
    if arg == '--N' and i+1 < len(sys.argv)-1:
        N = int(sys.argv[i+2])

print("="*68)
print(f"WHITEHEAD INTEGRAL FOR S3-LIFT HOPF CHARGE  (N={N})")
print("="*68)

# ── trefoil constants (Paper XV) ──────────────────────────────────────────
R0, r0, C_star = 3.0, 0.874, 2.5062
NT = 8000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz = r0*np.sin(3*t_arr)
tree = KDTree(np.stack([Gx,Gy,Gz], axis=1))

def f0(rho_hat):
    return 2*np.arctan(np.maximum(rho_hat, 1e-9)**(-C_star))

def nearest_two(pts, k=80, min_sep=np.pi/2):
    dists, idxs = tree.query(pts, k=k)
    t1 = t_arr[idxs[:,0]]; d1 = dists[:,0]
    dt = np.abs(((t_arr[idxs] - t1[:,None] + np.pi) % (2*np.pi)) - np.pi)
    mask = dt > min_sep
    fi = np.argmax(mask, axis=1)
    has2 = mask.any(axis=1)
    d2 = np.where(has2, dists[np.arange(len(pts)), fi], 1e6)
    t2 = t_arr[idxs[np.arange(len(pts)), fi]]
    return t1, d1, t2, d2

# ── doublet builders ───────────────────────────────────────────────────────
def Z_s3lift(pts):
    """Global S3-lift (Construction constr:global_s3lift, Paper XVI)."""
    t1, d1, t2, d2 = nearest_two(pts)
    f1 = f0(np.clip(d1, 1e-6, None)*C_star)
    f2 = f0(np.clip(d2, 1e-6, None)*C_star)
    w1 = 1.0/np.clip(d1, 1e-6, None)**2
    w2 = 1.0/np.clip(d2, 1e-6, None)**2
    z1u = (w1*np.cos(f1/2) + w2*np.cos(f2/2)).astype(complex)
    z2u = w1*np.sin(f1/2)*np.exp(3j*t1) + w2*np.sin(f2/2)*np.exp(3j*t2)
    mag = np.sqrt(np.abs(z1u)**2 + np.abs(z2u)**2)
    return z1u/mag, z2u/mag

def Z_hopf1(pts):
    """Inverse-stereographic Q_H=1 reference (smooth, no vortex, |Z|=1)."""
    x, y, z = pts[:,0], pts[:,1], pts[:,2]
    r2 = x**2 + y**2 + z**2
    z1 = ((1 - r2) + 2j*z) / (1 + r2)
    z2 = 2*(x + 1j*y)   / (1 + r2)
    return z1, z2   # already normalised: |z1|^2+|z2|^2 = 1

# ── Whitehead integral ─────────────────────────────────────────────────────
def whitehead(Zfunc, N, L, label):
    """Q_H = (1/4pi^2) int A.(curl A) d^3x,  A_i = Im(Z† d_i Z)."""
    h = 2*L/(N-1)
    ax = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    t0 = time.time()
    print(f"\n  [{label}] building Z on {N}^3 = {N**3:,} grid points ...",
          flush=True)
    z1, z2 = Zfunc(pts)
    z1 = z1.reshape(N,N,N); z2 = z2.reshape(N,N,N)
    print(f"    Z built in {time.time()-t0:.1f}s", flush=True)

    # Berry connection  A_i = Im(Z† d_i Z)
    dz1 = [np.gradient(z1, h, axis=a) for a in range(3)]
    dz2 = [np.gradient(z2, h, axis=a) for a in range(3)]
    A = [np.imag(np.conj(z1)*dz1[a] + np.conj(z2)*dz2[a]) for a in range(3)]

    # curl A
    curlA = [
        np.gradient(A[2], h, axis=1) - np.gradient(A[1], h, axis=2),
        np.gradient(A[0], h, axis=2) - np.gradient(A[2], h, axis=0),
        np.gradient(A[1], h, axis=0) - np.gradient(A[0], h, axis=1),
    ]

    integrand = sum(A[i]*curlA[i] for i in range(3))
    QH = h**3 * np.sum(integrand) / (4*np.pi**2)
    print(f"    Q_H = {QH:.6f}  (raw integral = {h**3*np.sum(integrand):.4f})")
    print(f"    total time {time.time()-t0:.1f}s")
    return QH

# ── run ────────────────────────────────────────────────────────────────────
L = 6.0
print(f"\nBox: [-{L},{L}]^3,  h = {2*L/(N-1):.4f},  tube_radius ~ {1/C_star:.3f}")
print(f"Points per tube radius ~ {(1/C_star)/( 2*L/(N-1) ):.1f}")
print(f"(Expect ~5-10% underestimate at N=80; multiply by ~1.05 to correct)")

print("\n--- VALIDATION: Q_H=1 inverse-stereographic Hopf doublet ---")
q1 = whitehead(Z_hopf1, N, L, "Q_H=1 Hopf")
print(f"  Expected: 1.0000  Got: {q1:.4f}  "
      f"({'OK' if abs(q1-1)<0.15 else 'CHECK'})")

print("\n--- MAIN: Global S3-lift ansatz ---")
qs = whitehead(Z_s3lift, N, L, "S3-lift")
print(f"\n  S3-lift Q_H = {qs:.4f}")
print(f"  (Correction for grid underestimate: ~{qs * (1/q1 if q1>0.1 else 1):.4f})")
print()
# classify
q_corr = qs * (1/q1) if q1 > 0.1 else qs
nearest_int = round(q_corr)
print(f"  Nearest integer: {nearest_int}")
if abs(q_corr - nearest_int) < 0.25:
    print(f"  RESULT: Q_H = {nearest_int}  (within 25% of integer after correction)")
else:
    print(f"  RESULT: Q_H ambiguous at this resolution — increase N")
