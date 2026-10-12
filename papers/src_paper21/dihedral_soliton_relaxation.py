#!/usr/bin/env python3.11
"""
dihedral_soliton_relaxation.py -- COARSE existence/stability probe for a D-type
(dihedral) dark-Hopfion in the Faddeev functional. DM/Quark crossover (src_paper21).

SCOPE (deliberately limited; see notes/fabric_cmb_hypothesis_and_plan.md (g)-(i)):
  Tests ONLY soliton EXISTENCE/STABILITY -- does a 3-fold (dihedral-flavoured)
  charged field configuration, under Faddeev gradient flow, (a) settle to a stable
  localized lump [candidate LIVES], (b) relax to an axially-symmetric / torus-type
  config [relaxes toward the SM torus-knot ladder], or (c) unwind (Q_H->0, E->0)?
  It says NOTHING about mass or abundance (those need the forming transition +
  cosmology; the mass is unpinned, update (h)). Coarse by design -- a HIGH-res run
  would hit the O2 tube-resolution wall and is NOT the goal.

METHOD (reuses validated pieces):
  * Field = Hopf doublet Z=(Z1,Z2), |Z|=1; n_hat = Z^dag sigma Z.
  * Faddeev energy  E = kappa * int |grad n|^2 + lambda * int F^2,
      F_ij = n . (d_i n x d_j n)   (both terms -> Vakulenko-Kapitanskii stability).
  * Q_H via the VALIDATED Whitehead integral (src_paper16/whitehead_hopf_charge.py):
      Q_H = (1/4pi^2) int A.(curl A),  A_i = Im(Z^dag d_i Z).
  * Gradient flow: torch Adam on the doublet, renormalise |Z|=1 each step.
  * IC: (n,m) toroidal Hopfion, doublet (a^m, b^n) from S^3 inverse-stereographic
      coords a=2(x+iy)/(1+r^2), b=((1-r^2)+2iz)/(1+r^2).  (n,m)=(1,3) = 3-fold
      azimuthal (dihedral-flavoured, C_3); (3,1) = axial Q=3; (1,1) = Q=1 VALIDATION.

MONITORS: E, Q_H, r_bar (energy-weighted radius, localization), A3 (cos3phi amplitude
of the energy density = the 3-fold order parameter: persists=dihedral, decays=axial).
"""
import numpy as np, torch, time, argparse

ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, default=40)
ap.add_argument('--L', type=float, default=6.0)
ap.add_argument('--steps', type=int, default=400)
ap.add_argument('--lr', type=float, default=5e-3)
ap.add_argument('--kappa', type=float, default=1.0)
ap.add_argument('--lam', type=float, default=1.0)
ap.add_argument('--ic', type=str, default='dihedral', choices=['dihedral', 'toroidal'])
ap.add_argument('--ncore', type=int, default=3, help='dihedral IC: cores in a ring (D_ncore)')
ap.add_argument('--Rc', type=float, default=1.6, help='dihedral IC: core ring radius')
ap.add_argument('--n', type=int, default=1, help='toroidal IC: poloidal winding')
ap.add_argument('--m', type=int, default=3, help='toroidal IC: azimuthal winding')
ap.add_argument('--log_every', type=int, default=25)
A = ap.parse_args()

dev = 'cpu'
ax = np.linspace(-A.L, A.L, A.N)
h = ax[1]-ax[0]
X, Y, Zc = np.meshgrid(ax, ax, ax, indexing='ij')
s = X**2 + Y**2 + Zc**2
phi = np.arctan2(Y, X)

# ---- IC ----
def Z_core(cx, cy, cz):
    """Validated Q=1 inverse-stereographic Hopf doublet centred at (cx,cy,cz)."""
    dx, dy, dz = X-cx, Y-cy, Zc-cz
    rr = dx**2 + dy**2 + dz**2
    Z1 = ((1-rr) + 2j*dz)/(1+rr)
    Z2 = (2*(dx + 1j*dy))/(1+rr)
    return Z1, Z2, rr

if A.ic == 'toroidal':
    # phase-winding ansatz (AXIALLY symmetric density -- kept for comparison)
    a = (2*(X+1j*Y))/(1+s); b = ((1-s)+2j*Zc)/(1+s)
    Araw, Braw = a**A.m, b**A.n
    nrm = np.sqrt(np.abs(Araw)**2+np.abs(Braw)**2)+1e-12
    Z1_0, Z2_0 = (Braw/nrm).astype(np.complex128), (Araw/nrm).astype(np.complex128)
else:
    # GENUINELY DIHEDRAL: ncore Q=1 cores on a ring (D_ncore symmetry), weighted-
    # superposed like the validated S3-lift (w_i = 1/d_i^2) -> 3-lumpy density.
    Z1_0 = np.zeros_like(s, dtype=np.complex128)
    Z2_0 = np.zeros_like(s, dtype=np.complex128)
    for k in range(A.ncore):
        th = 2*np.pi*k/A.ncore
        Z1k, Z2k, rr = Z_core(A.Rc*np.cos(th), A.Rc*np.sin(th), 0.0)
        w = 1.0/(rr + 0.25)
        Z1_0 += w*Z1k; Z2_0 += w*Z2k
    nrm = np.sqrt(np.abs(Z1_0)**2+np.abs(Z2_0)**2)+1e-12
    Z1_0 /= nrm; Z2_0 /= nrm

def nhat_from_Z(Z1, Z2):
    n1 = 2*(Z1.real*Z2.real + Z1.imag*Z2.imag)   # 2 Re(conj Z1 Z2)
    n2 = 2*(Z1.real*Z2.imag - Z1.imag*Z2.real)   # 2 Im(conj Z1 Z2)
    n3 = (Z1.real**2 + Z1.imag**2) - (Z2.real**2 + Z2.imag**2)
    return n1, n2, n3

def whitehead_QH(Z1, Z2):
    """Q_H = (1/4pi^2) int A.(curl A), A_i = Im(Z^dag d_i Z). (validated formula)"""
    dZ1 = [np.gradient(Z1, h, axis=i) for i in range(3)]
    dZ2 = [np.gradient(Z2, h, axis=i) for i in range(3)]
    Ac = [np.imag(np.conj(Z1)*dZ1[i] + np.conj(Z2)*dZ2[i]) for i in range(3)]
    curl = [np.gradient(Ac[2],h,axis=1)-np.gradient(Ac[1],h,axis=2),
            np.gradient(Ac[0],h,axis=2)-np.gradient(Ac[2],h,axis=0),
            np.gradient(Ac[1],h,axis=0)-np.gradient(Ac[0],h,axis=1)]
    integ = sum(Ac[i]*curl[i] for i in range(3))
    return h**3*np.sum(integ)/(4*np.pi**2)

# ---- torch energy ----
def d(f, axis):   # central difference via roll (periodic; field ~ vacuum at boundary)
    return (torch.roll(f, -1, axis) - torch.roll(f, 1, axis)) / (2*h)

def energy_and_density(Zt):
    Z1 = Zt[0] + 1j*Zt[1]; Z2 = Zt[2] + 1j*Zt[3]
    nn = (Z1.conj()*Z1 + Z2.conj()*Z2).real.sqrt()
    Z1 = Z1/nn; Z2 = Z2/nn
    n1 = 2*(Z1.real*Z2.real + Z1.imag*Z2.imag)
    n2 = 2*(Z1.real*Z2.imag - Z1.imag*Z2.real)
    n3 = (Z1.real**2+Z1.imag**2) - (Z2.real**2+Z2.imag**2)
    n = torch.stack([n1, n2, n3])
    dn = [torch.stack([d(n[c], ax_) for c in range(3)]) for ax_ in range(3)]  # dn[ax][comp]
    grad2 = sum((dn[ax_]**2).sum(0) for ax_ in range(3))                       # |grad n|^2
    def cross(u, v): return torch.stack([u[1]*v[2]-u[2]*v[1],
                                         u[2]*v[0]-u[0]*v[2],
                                         u[0]*v[1]-u[1]*v[0]])
    F2 = torch.zeros_like(grad2)
    for i in range(3):
        for j in range(i+1, 3):
            Fij = (n*cross(dn[i], dn[j])).sum(0)
            F2 = F2 + Fij**2
    dens = A.kappa*grad2 + A.lam*F2
    return (h**3)*dens.sum(), dens

# order parameters
r = torch.tensor(np.sqrt(s), dtype=torch.float64)
cos3 = torch.tensor(np.cos(3*phi), dtype=torch.float64)
def monitors(dens):
    w = dens.detach(); tot = w.sum()
    rbar = (w*r).sum()/tot
    A3 = (w*cos3).sum()/tot            # cos(3phi) projection of energy density
    return float(rbar), float(A3), float(tot*h**3)

# ---- validate IC charge, then flow ----
print("="*74)
desc = f"dihedral {A.ncore}-core (D_{A.ncore})" if A.ic=='dihedral' else f"toroidal (n,m)=({A.n},{A.m})"
print(f"  DIHEDRAL SOLITON RELAXATION  IC={desc}  N={A.N} L={A.L} steps={A.steps}")
print("="*74)
QH0 = whitehead_QH(Z1_0, Z2_0)
print(f"  IC Q_H (Whitehead) = {QH0:.3f}   (grid underestimate ~ -10-30% at this N)")
# Q_H=1 validation
a1 = (2*(X+1j*Y))/(1+s); b1 = ((1-s)+2j*Zc)/(1+s)
print(f"  VALIDATION Q_H=1 ref = {whitehead_QH(b1, a1):.3f}  (expect ~0.9 at this N)")

Zt = torch.tensor(np.stack([Z1_0.real, Z1_0.imag, Z2_0.real, Z2_0.imag]),
                  dtype=torch.float64, requires_grad=True)
opt = torch.optim.Adam([Zt], lr=A.lr)
print(f"\n  {'step':>5} {'E':>12} {'Q_H':>7} {'r_bar':>7} {'A3(cos3phi)':>12}")
t0 = time.time()
for step in range(A.steps+1):
    opt.zero_grad()
    E, dens = energy_and_density(Zt)
    if step < A.steps:
        E.backward(); opt.step()
        with torch.no_grad():   # renormalise |Z|=1 (kill norm/gauge drift)
            nn = (Zt[0]**2+Zt[1]**2+Zt[2]**2+Zt[3]**2).sqrt()
            Zt.data = Zt.data/nn
    if step % A.log_every == 0 or step == A.steps:
        with torch.no_grad():
            rbar, A3, Eval = monitors(dens)
            Zc1 = (Zt[0]+1j*Zt[1]).numpy(); Zc2 = (Zt[2]+1j*Zt[3]).numpy()
            nn = np.sqrt(np.abs(Zc1)**2+np.abs(Zc2)**2)+1e-12
            qh = whitehead_QH(Zc1/nn, Zc2/nn)
        print(f"  {step:>5} {Eval:>12.4f} {qh:>7.3f} {rbar:>7.3f} {A3:>12.5f}")
print(f"\n  elapsed {time.time()-t0:.1f}s")
print("  READ: Q_H preserved + E plateau + A3 finite -> stable DIHEDRAL lump (candidate LIVES)")
print("        A3 -> 0 (axial) -> relaxes toward torus ladder ; Q_H->0 / E->0 -> unwinds")
