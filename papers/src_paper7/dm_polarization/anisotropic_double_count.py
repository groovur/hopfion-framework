#!/usr/bin/env python3.11
"""
anisotropic_double_count.py -- semi-analytic test of the late-time DOUBLE-COUNTING question
(notes (xx),(yy),(zz)): stable dark knots (Omega~0.26, CDM-like halo) PLUS the director strain
(P7 flat-curve mechanism), with the director ANISOTROPY resolved per probe.

COMPONENTS (per galaxy, spherical toy):
  baryons : exponential disk, spherical-equivalent M_b(<r) = M_b[1-(1+x)e^-x], x=r/R_d,
            R_d = 3.5 kpc (M_b/5e10)^0.3.
  knots   : NFW halo, M_200 = f_knot * M_b, concentration c=10 (f_knot = 0 / 1 / cosmic 5.36).
  strain  : the imposed M~r disclination halo (local solvers cannot hold it -- prior notes). Two
            physical ingredients: AMPLITUDE set by the SOURCE mass, ONSET governed by the TOTAL Newtonian
            g (prior result 2: halo held until total g ~ a_0, so a knot halo pushes it outward):
              g_s = zeta * g_src * [ 1/(1 - exp(-sqrt(g_tot/a_0))) - 1 ],
              g_src = G (M_b + kappa_k M_knot)(<r)/r^2,  g_tot = G (M_b + M_knot)(<r)/r^2.
            Deep regime (g_tot = g_src << a_0): g_s -> sqrt(g_src a_0) = the M~r flat-curve halo.
            With no knots and zeta=1 this IS the empirical RAR (calibration = how P7 anchors a_0).
            kappa_k = surviving fraction of knot-SOURCED strain after the collective hedgehog
            cancellation (prior result 1; 0 = fully cancelled). zeta = the UNDERIVED strain/a_0
            coefficient (P7:sec:accel_scale) -> fitted.
            (A first version with a SHARP onset at the outermost a_0 crossing + linear M growth FAILED
            the gate: zeta_eff=0.62, rms 0.18 dex, massive -0.14 dex -- an artifact of the stand-in.)

ANISOTROPY (radial hedgehog n = r_hat, prior notes; P7:eq:meff): a probe moving at angle psi to r_hat
sees the strain with weight
    B (inertial tensor M_ij = (3/2) m_eff (I - r r)):  w(psi) = 1.5 sin^2 psi
    A (sin^4 coupling, normalised to isotropic mean 1): w(psi) = (15/8) sin^4 psi
    ISO (control = MOND+CDM, no anisotropy):           w = 1
  rotation curves (tangential, psi=90): w_t = 1.5 (B) / 1.875 (A) / 1
  radial-orbit dynamics (psi=0):          w_r = 0 (A,B) / 1
  lensing (gravitating Frank energy = scalar trace):  w_L = 1
The kinematic (rotation) RAR is fitted -> it fixes zeta_eff = w_t*zeta; the anisotropy then
PREDICTS lensing (strain x 1/w_t) and radial dynamics (no strain).

TARGET: empirical kinematic RAR g_obs = g_bar / (1 - exp(-sqrt(g_bar/a_0))), a_0=1.2e-10 m/s^2.
VALIDATION GATE: f_knot=0, kappa_k=0 must reproduce the RAR with zeta_eff~1 and small scatter.

HONEST SCOPE: toy (spherical baryons, fixed c, sharp onset, no feedback, no intrinsic scatter in
inputs). The M_ij -> per-probe mapping is the B HYPOTHESIS made quantitative: P7 uses only the scalar
m_eff, and direction-dependent dynamics vs isotropic gravitation must be checked against P7's WEP
sector. theta for an exactly aligned hedgehog (sin^4 -> 0) is unresolved (notes (zz)).
"""
import numpy as np
from scipy.optimize import minimize_scalar

G      = 4.302e-6                    # kpc (km/s)^2 / Msun
MS2    = 3.2408e-14                  # 1 (km/s)^2/kpc in m/s^2
A0     = 1.2e-10/MS2                 # a_0 in (km/s)^2/kpc  (~3703)
RHOC   = 136.0                       # Msun/kpc^3 (h=0.7)
F_COSMIC = 0.1200/0.02237            # 5.36

def mcgaugh(gb):  return gb/(1.0-np.exp(-np.sqrt(gb/A0)))
def M_bar(r,Mb,Rd): x=r/Rd; return Mb*(1-(1+x)*np.exp(-x))
def M_nfw(r,M200,c=10.0):
    if M200<=0: return np.zeros_like(r)
    R200=(3*M200/(4*np.pi*200*RHOC))**(1/3); rs=R200/c
    m=lambda y: np.log(1+y)-y/(1+y)
    return M200*m(r/rs)/m(c)

def galaxy(Mb, f_knot, kappa_k, r):
    Rd = 3.5*(Mb/5e10)**0.3
    Mb_r = M_bar(r,Mb,Rd); M200=f_knot*Mb; Mk_r=M_nfw(r,M200)
    gN_tot = G*(Mb_r+Mk_r)/r**2
    g_src  = G*(Mb_r+kappa_k*Mk_r)/r**2
    gs1 = g_src*(1.0/(1.0-np.exp(-np.sqrt(gN_tot/A0)))-1.0)   # strain at zeta=1
    return dict(gbar=G*Mb_r/r**2, gN=gN_tot, gs1=gs1, Rd=Rd)

MASSES = np.logspace(7,11.5,16)
def population(f_knot,kappa_k,kin=True):
    out=[]
    for Mb in MASSES:
        Rd=3.5*(Mb/5e10)**0.3
        r = np.logspace(np.log10(0.3*Rd), np.log10(15*Rd if kin else 300.0), 60)
        g=galaxy(Mb,f_knot,kappa_k,r); g['Mb']=Mb; g['r']=r; out.append(g)
    return out

def resid(pop,zeta_eff):
    return np.concatenate([np.log10((g['gN']+zeta_eff*g['gs1'])/mcgaugh(g['gbar'])) for g in pop])

def fit(pop):
    f=lambda z: np.mean(resid(pop,z)**2)
    z=minimize_scalar(f,bounds=(0.0,5.0),method='bounded').x
    return z, np.sqrt(f(z))

def mass_trend(pop,z):
    lo=[np.mean(np.log10((g['gN']+z*g['gs1'])/mcgaugh(g['gbar']))) for g in pop if g['Mb']<1e9]
    hi=[np.mean(np.log10((g['gN']+z*g['gs1'])/mcgaugh(g['gbar']))) for g in pop if g['Mb']>3e10]
    return np.mean(lo), np.mean(hi)

def probe_ratio(f_knot,kappa_k,zeta_eff,w_t,gb_target=1e-12/MS2):
    """At large radius (g_bar ~ 1e-12 m/s^2): lensing/rotation and radial/rotation dark-excess ratios."""
    pop=population(f_knot,kappa_k,kin=False); L=[]; R=[]
    for g in pop:
        i=np.argmin(np.abs(g['gbar']-gb_target))
        if abs(np.log10(g['gbar'][i]/gb_target))>0.3: continue
        rot = g['gN'][i]+zeta_eff*g['gs1'][i]
        lens= g['gN'][i]+(zeta_eff/w_t)*g['gs1'][i]
        rad = g['gN'][i]+(zeta_eff*(0.0 if w_t!=1.0 else 1.0))*g['gs1'][i]
        gb=g['gbar'][i]
        L.append((lens-gb)/(rot-gb)); R.append((rad-gb)/(rot-gb))
    return (np.median(L),np.median(R),len(L)) if L else (np.nan,np.nan,0)

if __name__=='__main__':
    print("="*92)
    print("  DOUBLE-COUNTING with director ANISOTROPY: knots (NFW) + strain (M~r, global a0 onset on TOTAL g)")
    print("="*92)
    print(f"  a0={A0:.0f} (km/s)^2/kpc;  galaxies M_b=1e7..3e11 Msun;  target = empirical kinematic RAR\n")
    print(f"  {'f_knot':>6} {'kappa_k':>7} | {'zeta_eff':>8} {'rms(dex)':>8} {'dwarf<1e9':>9} {'massive>3e10':>12} |"
          f" {'lens/rot B':>10} {'lens/rot A':>10} {'lens/rot ISO':>12} | {'radial/rot (A,B)':>16}")
    for f_knot in [0.0, 1.0, F_COSMIC]:
        for kappa_k in ([0.0] if f_knot==0 else [0.0,0.1]):
            pop=population(f_knot,kappa_k)
            z,rms=fit(pop); lo,hi=mass_trend(pop,z)
            LB,RB,_=probe_ratio(f_knot,kappa_k,z,1.5)
            LA,_,_ =probe_ratio(f_knot,kappa_k,z,1.875)
            LI,_,n =probe_ratio(f_knot,kappa_k,z,1.0)
            print(f"  {f_knot:>6.2f} {kappa_k:>7.2f} | {z:>8.3f} {rms:>8.3f} {lo:>+9.3f} {hi:>+12.3f} |"
                  f" {LB:>10.2f} {LA:>10.2f} {LI:>12.2f} | {RB:>16.2f}")
    print("\n  columns: zeta_eff = fitted tangential strain coefficient (w_t*zeta); rms = RAR scatter about")
    print("  the empirical relation (observed intrinsic ~0.06 dex); dwarf/massive = mean residual (dex) by mass;")
    print("  lens/rot = DARK-excess ratio lensing:rotation at g_bar~1e-12 m/s^2 (1 = isotropic, as MOND/CDM);")
    print("  radial/rot = dark excess seen by purely radial orbits relative to rotation.")
