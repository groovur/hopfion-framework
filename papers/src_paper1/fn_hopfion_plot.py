#!/usr/bin/env python3
"""
fn_hopfion_plot.py  —  Recompute and plot the FN Hopfion saddle profile
========================================================================

Runs the saddle-snapshot flow at mu* = 3 - phi on the LARGE grid,
saves the profile, and produces a figure with:

  Panel A  — f(r,z) filled contours with f=pi/2 isoline (Hopf tube)
  Panel B  — |grad f|^2 energy density (where the soliton lives)
  Panel C  — sin^4(f) * |grad f|^2 anisotropic weight (J2a integrand)
  Panel D  — 1D radial slice f(rho) at fixed z=0, vs BS rational ansatz

Usage:
  python3 fn_hopfion_plot.py [--outdir PATH] --R0 3.0 ]

Output:
  hopfion_saddle_profile.npy   (the raw f array, LARGE grid)
  hopfion_saddle_profile.png   (the four-panel figure)
  hopfion_saddle_stats.txt     (key numbers: ratio, br, J2iso_hat, etc.)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import time, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', default='.')
parser.add_argument('--R0',         type=float, default=3.0,
               help='Torus major radius (default 3.0)')
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
phi   = (1 + 5**0.5) / 2
lam   = phi**6          # 17.9443
phi9  = lam**1.5        # 76.0132
BS_J2 = 0.6143
MU    = 3.0 - phi       # 1.38197  (exact: 1 + 1/phi^2)
ep    = 1e-14

# ── LARGE grid ────────────────────────────────────────────────────────────────
Nr, Nz = 256, 256
h      = 0.12
R0     = args.R0

r_   = h * (np.arange(Nr) + 0.5)
z_   = h * np.arange(Nz)
R, Z = np.meshgrid(r_, z_, indexing='ij')
vol        = 2 * np.pi * 2 * R * h**2
vol[:, 0] *= 0.5
D0   = (R - R0)**2 + Z**2 + ep
iR0  = max(1, int(round(R0 / h)) - 1)

print(f"Grid: {Nr}×{Nz},  h={h},  R0={R0},  R_max={Nr*h:.1f}")
print(f"mu*  = 3 - phi = {MU:.10f}")
print(f"lambda = phi^6 = {lam:.6f}")
print()

# ── Core ──────────────────────────────────────────────────────────────────────
def compute(f, mu):
    f  = np.clip(f, 0.0, np.pi)
    fr = np.gradient(f, h, axis=0)
    fz = np.gradient(f, h, axis=1); fz[:, 0] = 0.0
    sf = np.sin(f); cf = np.cos(f)
    s2 = sf**2; s3 = s2*sf; s4 = s2*s2; s5 = s4*sf
    D  = D0; A = 1.0/D + 1.0/R**2
    g2 = fr**2 + fz**2; fDG = fr*(R - R0) + fz*Z
    F13 = -s2/D*fDG; F12 = s2/R*fr; F23 = s2/R*fz
    kern  = g2 + s2*A
    J2a   = float(np.sum(s4*kern*vol))
    J2iso = float(np.sum(   kern*vol))
    J4    = float(np.sum((F13**2 + F12**2 + F23**2)*vol))
    K = J2a + mu*J2iso
    loc_K  = ((4*s3*cf + mu*2*sf*cf)*kern + (s4+mu)*2*sf*cf*A)*vol
    loc_J4 = 2*(F13*(-2*sf*cf/D*fDG)+F12*(2*sf*cf/R*fr)+F23*(2*sf*cf/R*fz))*vol
    fK_r   = (2*s4+2*mu)*fr*vol; fK_z = (2*s4+2*mu)*fz*vol; fK_z[:,0]=0
    fJ4_r  = 2*(F13*(-s2/D*(R-R0))+F12*(s2/R))*vol
    fJ4_z  = 2*(F13*(-s2/D*Z)+F23*(s2/R))*vol; fJ4_z[:,0]=0
    def div(a, b): return np.gradient(a, h, axis=0) + np.gradient(b, h, axis=1)
    dK  = loc_K  - div(fK_r,  fK_z)
    dJ4 = loc_J4 - div(fJ4_r, fJ4_z)
    Force = -(dK*J4 + K*dJ4)
    return J2a, J2iso, J4, K*J4, Force

def bc(f):
    f = np.clip(f, 0.0, np.pi)
    f[0,:]=0; f[-1,:]=0; f[:,-1]=0; f[iR0,0]=np.pi
    return f

def s_opt(J2a, J2iso, J4, mu):
    K = J2a + mu*J2iso
    return float(np.sqrt(lam*J4/K)) if K > 0 and J4 > 0 else 0.0

# ── Reference gm_iso ─────────────────────────────────────────────────────────
f_ref = bc(np.clip(np.pi/(1+D0/(0.82*R0)**2), 0, np.pi))
_,J2iso_r,J4_r,_,_ = compute(f_ref, mu=0)
GM_ISO = float(np.sqrt(J2iso_r*J4_r))
print(f"gm_iso = {GM_ISO:.6f}")

def br(J2a, J2iso, J4, mu):
    so    = s_opt(J2a, J2iso, J4, mu)
    ratio = J2a / J2iso if J2iso > 0 else 0
    # br_grid: formula used in all sessions (gm_iso normalised)
    br_g  = phi9 * np.sqrt(max(J2a*J4, 0)) * BS_J2 / GM_ISO
    # br_true: physically normalised — uses BS_J2 directly
    br_t  = phi9 * BS_J2 * np.sqrt(ratio*(ratio+mu)/lam) if J2iso > 0 else 0
    # J2iso_hat: the saddle profile's J2iso in BS natural units
    # Derived from br_grid = phi9 * J2iso_hat * sqrt(C/lam), so:
    C = ratio*(ratio+mu) if J2iso > 0 else 1e-30
    J2iso_hat_val = br_g / (phi9 * np.sqrt(C/lam)) if C > 0 else 0
    return br_g, br_t, so, ratio, J2iso_hat_val

# ── Flow ──────────────────────────────────────────────────────────────────────
print(f"\nRunning saddle-snapshot flow at mu* = 3-phi = {MU:.8f} ...")
f = bc(np.clip(np.pi/(1+D0/(0.6*R0)**2), 0, np.pi))
dt = 1e-7
MAX_STEPS = 10_000
PEV       = 1_000

best_dist  = 1e30
best_f     = f.copy()
best_stats = {}
t0 = time.time()

for step in range(1, MAX_STEPS+1):
    J2a, J2iso, J4, Eg, Force = compute(f, MU)
    br_g, br_t, so, ratio, J2isohat = br(J2a, J2iso, J4, MU)

    dist = abs(so - 1.0)
    if dist < best_dist:
        best_dist  = dist
        best_f     = f.copy()
        best_stats = dict(step=step, br_grid=br_g, br_true=br_t,
                          s_opt=so, ratio=ratio, J2iso_hat=J2isohat,
                          J2a=J2a, J2iso=J2iso, J4=J4)

    if step % PEV == 0 or step == 1:
        print(f"  step={step:>6d}  br={br_g:>7.4f}  s_opt={so:.4f}  "
              f"ratio={ratio:.4f}  J2iso_hat={J2isohat:.5f}  [{time.time()-t0:.0f}s]",
              flush=True)

    f_try = bc(f + dt*Force)
    Eg_try = compute(f_try, MU)[3]
    if Eg_try < Eg:
        f = f_try; dt = min(dt*1.01, 5e-4)
    else:
        dt *= 0.7
        if dt < 1e-18:
            print(f"  dt exhausted at step {step}"); break

s = best_stats
print(f"\nSaddle snapshot at step {s['step']}:")
print(f"  s_opt     = {s['s_opt']:.5f}")
print(f"  J2iso_hat = {s['J2iso_hat']:.5f}  (BS_J2 = {BS_J2})")
print(f"  ratio     = {s['ratio']:.5f}")
print(f"  br_grid   = {s['br_grid']:.5f}")
print(f"  br_true   = {s['br_true']:.5f}")

# ── Save profile ──────────────────────────────────────────────────────────────
npy_path = os.path.join(args.outdir, "hopfion_saddle_profile.npy")
np.save(npy_path, best_f)
print(f"\nProfile saved: {npy_path}")

# ── Derived fields ────────────────────────────────────────────────────────────
f = best_f
fr = np.gradient(f, h, axis=0)
fz = np.gradient(f, h, axis=1); fz[:,0] = 0.0
sf = np.sin(f); s2 = sf**2; s4 = s2*s2
D  = D0; A = 1.0/D + 1.0/R**2
g2 = fr**2 + fz**2
kern = g2 + s2*A           # full isotropic integrand (per unit vol / 2piR drdz)
anis = s4 * kern            # sin^4 weighted (J2a integrand)

# ── FIGURE ────────────────────────────────────────────────────────────────────
# Plot domain: r in [0, 12], z in [0, 10]  (soliton lives near r=R0=3)
ir_max = int(12.0/h); iz_max = int(10.0/h)
ir_max = min(ir_max, Nr); iz_max = min(iz_max, Nz)

r_plot = r_[:ir_max]
z_plot = z_[:iz_max]
f_plot = f[:ir_max, :iz_max]
kr_plot = kern[:ir_max, :iz_max]
an_plot = anis[:ir_max, :iz_max]

# Symmetric extension to negative z for visual clarity
z_full = np.concatenate([-z_plot[::-1][:-1], z_plot])
f_full = np.concatenate([f_plot[:, ::-1][:,:-1], f_plot], axis=1)
kr_full = np.concatenate([kr_plot[:, ::-1][:,:-1], kr_plot], axis=1)
an_full = np.concatenate([an_plot[:, ::-1][:,:-1], an_plot], axis=1)

fig, axes = plt.subplots(2, 2, figsize=(13, 10),
                          gridspec_kw={'hspace': 0.35, 'wspace': 0.35})

CMAP_F    = 'RdYlBu_r'
CMAP_DENS = 'inferno'
ISOLINE_COLOR = '#00e5ff'   # cyan — stands out on both colormaps
ISOLINE_LW    = 2.0

titles = [
    r'$f(r,z)$  [field angle]',
    r'$|\nabla f|^2 + \sin^2\!f \cdot A$  [isotropic energy density]',
    r'$\sin^4\!f \cdot (|\nabla f|^2 + \sin^2\!f \cdot A)$  [$J_{2a}$ integrand]',
    r'Radial slice  $f(\varrho)$  at $z=0$',
]

for ax, data, title, cmap in zip(
        axes.flat[:3],
        [f_full, kr_full, an_full],
        titles[:3],
        [CMAP_F, CMAP_DENS, CMAP_DENS]):

    if cmap == CMAP_F:
        vmin, vmax = 0, np.pi
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        label = r'$f$ (rad)'
        ticks = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
        tick_labels = [r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
    else:
        vmax = np.percentile(data, 99)
        norm = mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=vmax)
        label = r'arb. units'
        ticks = None; tick_labels = None

    im = ax.pcolormesh(z_full, r_plot, data,
                       cmap=cmap, norm=norm, shading='auto', rasterized=True)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(label, fontsize=9)
    if ticks is not None:
        cb.set_ticks(ticks); cb.set_ticklabels(tick_labels)

    # f = pi/2 isoline (Hopf tube cross-section)
    ax.contour(z_full, r_plot, data if cmap == CMAP_F else f_full,
               levels=[np.pi/2],
               colors=[ISOLINE_COLOR], linewidths=ISOLINE_LW,
               linestyles='solid')
    # f = pi isoline (soliton core)
    ax.contour(z_full, r_plot, f_full,
               levels=[np.pi * 0.95],
               colors=['white'], linewidths=1.0, linestyles='dashed', alpha=0.7)

    ax.axvline(0, color='gray', lw=0.5, ls=':')
    ax.set_xlabel(r'$z$', fontsize=11)
    ax.set_ylabel(r'$r$', fontsize=11)
    ax.set_title(title, fontsize=10, pad=6)
    ax.set_xlim(z_full[0], z_full[-1])
    ax.set_ylim(0, r_plot[-1])

    # Mark torus centre
    ax.plot(0, R0, '+', color='white', ms=8, mew=1.5)

# Panel D: 1D radial slice
ax = axes[1, 1]
# rho = distance from torus ring centre (R0, 0)
rho_slice = r_[:ir_max] - R0   # along z=0
f_slice = f[:ir_max, 0]        # z=0 column

# BS rational ansatz for comparison
a_BS = 0.82 * R0
rho_BS = np.linspace(-R0*0.95, r_plot[-1]-R0, 400)
r_BS = rho_BS + R0
f_BS = np.where(r_BS > 0, np.pi / (1 + (r_BS - R0)**2/a_BS**2), 0)

ax.plot(rho_BS, f_BS/np.pi, '--', color='#888888', lw=1.8,
        label=r'BS ansatz ($q=0.82$)')
ax.plot(rho_slice, f_slice/np.pi, color='#e63946', lw=2.2,
        label=r'Physical Hopfion ($\mu^*=3-\varphi$)')
ax.axvline(0, color='gray', lw=0.5, ls=':')
ax.axhline(0.5, color=ISOLINE_COLOR, lw=1.2, ls='--', alpha=0.8,
           label=r'$f = \pi/2$ (Hopf tube)')
ax.set_xlabel(r'$\varrho = r - R_0$', fontsize=11)
ax.set_ylabel(r'$f / \pi$', fontsize=11)
ax.set_title(titles[3], fontsize=10, pad=6)
ax.set_xlim(-R0*0.9, r_plot[-1]-R0)
ax.set_ylim(-0.05, 1.10)
ax.legend(fontsize=8, loc='upper right')
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))

# ── Global title ──────────────────────────────────────────────────────────────
fig.suptitle(
    rf'FN Hopfion saddle-point profile at $\mu^* = 3-\varphi = {MU:.5f}$'
    rf'  |  $\beta_0\rho_{{\mathrm{{CMB}}}}={s["br_true"]:.3f}$'
    rf'  |  ratio$^*={s["ratio"]:.4f}$'
    rf'  |  $\widehat{{J}}_{{2,\mathrm{{iso}}}}={s["J2iso_hat"]:.4f}$ (BS$={BS_J2}$)',
    fontsize=11, y=0.98
)

# ── Annotation: cyan = Hopf tube ──────────────────────────────────────────────
for ax in axes.flat[:3]:
    ax.annotate(r'$f=\pi/2$ (Hopf tube)',
                xy=(0, R0), xytext=(2.5, 7.5),
                fontsize=7.5, color=ISOLINE_COLOR,
                arrowprops=dict(arrowstyle='->', color=ISOLINE_COLOR,
                                lw=0.8, shrinkA=2, shrinkB=2))

filename = f"hopfion_saddle_profile_R0{R0}.png"
png_path = os.path.join(args.outdir, filename)
fig.savefig(png_path, dpi=160, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Figure saved: {png_path}")

# ── Stats file ────────────────────────────────────────────────────────────────
txt_path = os.path.join(args.outdir, "hopfion_saddle_stats.txt")
with open(txt_path, 'w') as fh:
    fh.write("FN Hopfion Saddle-Point Profile — Session 24\n")
    fh.write("="*50 + "\n\n")
    fh.write(f"Grid:           {Nr}×{Nz},  h={h},  R0={R0}\n")
    fh.write(f"mu*:            {MU:.10f}  (= 3 - phi = 1 + 1/phi^2)\n")
    fh.write(f"Snapshot step:  {s['step']}\n\n")
    fh.write(f"s_opt:          {s['s_opt']:.6f}  (Derrick balance)\n")
    fh.write(f"J2iso_hat:      {s['J2iso_hat']:.6f}  (BS_J2 = {BS_J2})\n")
    fh.write(f"ratio*:         {s['ratio']:.6f}  (J2a/J2iso)\n")
    fh.write(f"br_grid:        {s['br_grid']:.6f}  (gm_iso normalised)\n")
    fh.write(f"br_true:        {s['br_true']:.6f}  (BS_J2 normalised)\n\n")
    fh.write(f"lambda = phi^6: {lam:.6f}\n")
    fh.write(f"phi^9:          {phi9:.6f}\n")
    fh.write(f"gm_iso:         {GM_ISO:.6f}\n")
print(f"Stats saved:  {txt_path}")
print("\nDone.")
