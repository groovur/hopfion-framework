#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

phi = (1 + 5**0.5) / 2

f = np.load('hopfion_saddle_profile.npy')
N = 256; h = 0.12; R0 = 3.0

# Correct coordinates from the solver:
# r_ = h*(np.arange(N) + 0.5)  -> starts at h/2
# z_ = h*np.arange(N)          -> starts at 0 (torus equatorial plane)
r_ = h * (np.arange(N) + 0.5)
z_ = h * np.arange(N)
R, Z = np.meshgrid(r_, z_, indexing='ij')

ep = 1e-14
D0  = (R - R0)**2 + Z**2 + ep
A   = 1.0/D0 + 1.0/np.where(R < ep, ep, R)**2
sf  = np.sin(f); sf2 = sf**2; sf4 = sf**4
fr  = np.gradient(f, h, axis=0)
fz  = np.gradient(f, h, axis=1)
kern = fr**2 + fz**2 + sf2 * A

iso_dens = kern
J2a_dens = sf4 * kern

# Clip for display
def trim(arr, pct=99.5):
    return np.clip(arr, 0, np.percentile(arr[arr>0], pct)) if arr.max()>0 else arr

# Radial profile at z_=0 (j=0, torus equatorial plane)
j0 = 0
rho_axis = r_ - R0  # centred on torus
f_slice  = f[:, j0] / np.pi

# BPS comparison
rho_pos  = rho_axis[rho_axis > 0]
f_bs_pos = 2 * np.arctan(0.82 / rho_pos) / np.pi

# ── Figure ────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 8.5), facecolor='white')
gs  = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.38)

cmap_f = 'RdYlBu_r'
cmap_e = 'inferno'

# Panels A–C: show region 0..12 in r, 0..12 in z (where soliton lives)
rmax_plot = 12.0; zmax_plot = 12.0
ri = r_ <= rmax_plot; zi = z_ <= zmax_plot

def pcolour(ax, field, cmap, title, xlabel=True, vmax=None):
    Rp = R[np.ix_(ri, zi)]; Zp = Z[np.ix_(ri, zi)]
    Fp = field[np.ix_(ri, zi)]
    vmax = vmax or np.percentile(Fp[Fp>0], 99.5)
    im = ax.pcolormesh(Zp.T, Rp.T, Fp.T, cmap=cmap,
                       vmin=0, vmax=vmax, rasterized=True, shading='auto')
    ax.contour(Zp.T, Rp.T, f[np.ix_(ri,zi)].T,
               levels=[np.pi/2], colors='cyan', linewidths=1.6)
    ax.set_ylabel(r'$r$', fontsize=12)
    if xlabel:
        ax.set_xlabel(r'$z$', fontsize=12)
    ax.set_title(title, fontsize=10.5)
    ax.set_xlim(0, zmax_plot); ax.set_ylim(0, rmax_plot)
    return im

# Panel A: field angle
ax0 = fig.add_subplot(gs[0,0])
im0 = pcolour(ax0, f, cmap_f, r'(a) $f(r,z)$ — field angle', vmax=np.pi)
cb0 = fig.colorbar(im0, ax=ax0, pad=0.02)
cb0.set_ticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
cb0.set_ticklabels(['$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'], fontsize=8)
# label the Hopf tube
ax0.annotate(r'$f{=}\pi/2$', xy=(0.8, 3.2), xytext=(3.5, 6.5),
             color='cyan', fontsize=9,
             arrowprops=dict(arrowstyle='->', color='cyan', lw=1.2))

# Panel B: isotropic energy density
ax1 = fig.add_subplot(gs[0,1])
im1 = pcolour(ax1, trim(iso_dens), cmap_e,
              r'(b) $|\nabla f|^2+\sin^2\!f\cdot A$ — isotropic density')
fig.colorbar(im1, ax=ax1, pad=0.02).set_label('arb.\ units', fontsize=8)

# Panel C: J2a integrand
ax2 = fig.add_subplot(gs[1,0])
im2 = pcolour(ax2, trim(J2a_dens), cmap_e,
              r'(c) $\sin^4\!f\cdot(|\nabla f|^2{+}\sin^2\!f\cdot A)$ — $J_{2a}$ integrand')
fig.colorbar(im2, ax=ax2, pad=0.02).set_label('arb.\ units', fontsize=8)

# Panel D: radial slice
ax3 = fig.add_subplot(gs[1,1])
ax3.plot(rho_axis, f_slice, color='C3', lw=2.2,
         label=r'Physical Hopfion ($\mu^*=3{-}\varphi$)')
ax3.plot(rho_pos,  f_bs_pos, 'k--', lw=1.6,
         label=r'BS ansatz ($q{=}0.82$)')
ax3.axhline(0.5, color='cyan', lw=1.3, ls='--', label=r'$f{=}\pi/2$')
ax3.axvline(0.0, color='gray', lw=0.8, ls=':')
ax3.set_xlabel(r'$\varrho = r - R_0$', fontsize=12)
ax3.set_ylabel(r'$f/\pi$', fontsize=12)
ax3.set_title(r'(d) Radial slice at $z=0$', fontsize=10.5)
ax3.set_xlim(-2.5, 9); ax3.set_ylim(-0.02, 1.05)
ax3.legend(fontsize=8.5, loc='upper right')
ax3.grid(True, alpha=0.3)

fig.suptitle(
    r'FN Hopfion saddle point ($256\times256$, $h=0.12$, $\beta^*=0.452$): '
    r'$\mu^*=3{-}\varphi$, $\;\beta_0\rho_\mathrm{CMB}=10.012$, '
    r'$\;J_{2\mathrm{iso}}^\mathrm{fb}/J_{2a}=\varphi$, '
    r'$\;s_\mathrm{opt}=2^{1/6}$ to $0.010\%$',
    fontsize=10, y=0.997)

plt.savefig('hopfion_profile_clean.png',
            bbox_inches='tight', facecolor='white', dpi=200)
print("Done")
