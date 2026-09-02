#!/usr/bin/env python3
r"""
3D defect-network simulation, INERTIAL (c_dir->0) director dynamics -- pin the frozen line-network spacing.
Successor to defect_network_sim.py (2D, point defects). See notes/single_galaxy_winding_analysis.md.

WHY 3D + why this model:
  - Disclinations are LINES in 3D, not points. The 2D point proxy is where the c~2-vs-c~1 ambiguity lived.
  - The 2D stirred runs were ARTIFACTS: advection was added to the ACCELERATION (-adv as a force), which is
    physically wrong (advection is first-order transport). Here stirring is modelled the way the Kibble-Zurek
    physics actually works: a RE-QUENCH -- rapidly reorient the director at scale xi_turb (a merger), then let
    the inertial c_dir->0 dynamics FREEZE the trapped defect lines. No ad hoc advection force.

MODEL: 3D complex order parameter psi=|psi|e^{i theta} (director orientation as a phase; disclinations=vortex
lines, +-1 windings). INERTIAL (damped-wave) dynamics, velocity-Verlet:
    chi d_t^2 psi = K grad^2 psi + (1-|psi|^2) psi/tau - gamma d_t psi
  written as  d_t^2 psi = c_dir^2 grad^2 psi + (1-|psi|^2)psi/tau - gamma d_t psi,  c_dir^2 = K/chi.
  c_dir->0 (SMALL c_dir = heavy/large chi) => defects are massive, don't self-propagate -> FROZEN network.
  Overdamped comparison (--overdamped: d_t psi = D grad^2 psi + core) retained -- the WRONG dynamics (coarsens).

QUENCH: IC = random director correlated at xi_turb (FFT-filtered) -> a Kibble-Zurek defect tangle at ~xi_turb.
  --requench K : K merger kicks over the run (blend toward a fresh xi_turb field, strength --kick) = repeated
  mergers re-injecting defects. Default (K=0) = single quench = the clean freeze test.

MEASURE: total defect-LINE length via plaquette windings in all 3 coordinate planes -> line density ->
  xi_net = k*sqrt(V/L_line). The geometric k cancels in the calibration-free ratio c = xi_net/xi_imprint, which
  is what we report: how much COARSER the frozen network is than the freshly-imprinted xi_turb tangle.
  c~1 => network freezes at the injection scale (recovers cluster deficit); c>>1 => coarsens (fails).

SMOKE: --smoke (small box). FULL: e.g.  python defect_network_sim_3d.py --N 96 --xi 12 --steps 8000
Uses numpy only (no scipy): /usr/local/Caskroom/miniconda/base/bin/python or any python3 with numpy.
"""
import argparse, os, time
import numpy as np

def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=96, help="grid size (cubic box L=N cells = R_cl)")
    p.add_argument("--xi", type=float, default=12.0, help="quench correlation scale xi_turb (cells)")
    p.add_argument("--steps", type=int, default=8000)
    p.add_argument("--dt", type=float, default=0.1)
    p.add_argument("--tau", type=float, default=1.0, help="core relaxation time")
    p.add_argument("--cdir", type=float, default=0.3, help="inertial wave speed c_dir=sqrt(K/chi); SMALL=heavy")
    p.add_argument("--gamma", type=float, default=0.02, help="inertial damping (small=nearly lossless)")
    p.add_argument("--overdamped", action="store_true", help="use OVERDAMPED d_t psi=D lap psi+core (WRONG dyn.)")
    p.add_argument("--D", type=float, default=0.15, help="overdamped diffusion")
    p.add_argument("--requench", type=int, default=0, help="number of merger re-quench kicks over the run")
    p.add_argument("--kick", type=float, default=0.6, help="re-quench blend strength (0..1)")
    p.add_argument("--measure-every", type=int, default=500)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    return p

def correlated_unit_field(N, xi, rng):
    """Random complex field with correlation length ~xi (FFT Gaussian filter), normalized to unit modulus."""
    noise = rng.standard_normal((N, N, N)) + 1j*rng.standard_normal((N, N, N))
    k = 2*np.pi*np.fft.fftfreq(N)                       # cells^-1
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    kmag = np.sqrt(KX**2 + KY**2 + KZ**2)
    filt = np.exp(-0.5*(kmag*xi)**2)                    # correlation length ~ xi cells
    field = np.fft.ifftn(np.fft.fftn(noise)*filt)
    return field/(np.abs(field) + 1e-12)

def lap3(a):
    return (np.roll(a,1,0)+np.roll(a,-1,0)+np.roll(a,1,1)+np.roll(a,-1,1)
            +np.roll(a,1,2)+np.roll(a,-1,2) - 6*a)

def _wrap(x):
    return (x + np.pi) % (2*np.pi) - np.pi

def defect_line_length(theta):
    """Total defect-line length ~ sum of |plaquette winding| over all 3 coordinate-plane orientations."""
    tot = 0
    for (a0, a1) in ((0,1), (1,2), (0,2)):             # in-plane axis pairs (normal = the third axis)
        d0 = _wrap(np.roll(theta,-1,a0) - theta)       # edge diff along a0
        d1 = _wrap(np.roll(theta,-1,a1) - theta)       # edge diff along a1
        winding = d0 + np.roll(d1,-1,a0) - np.roll(d0,-1,a1) - d1
        tot += int(np.abs(np.round(winding/(2*np.pi)).astype(int)).sum())
    return tot                                          # pierced plaquettes summed over 3 orientations

def main():
    a = build_parser().parse_args()
    if a.smoke:
        a.N, a.xi, a.steps, a.measure_every = 48, 8.0, 1500, 300
    rng = np.random.default_rng(a.seed)
    N = a.N
    psi = correlated_unit_field(N, a.xi, rng)
    vel = np.zeros_like(psi)
    P0 = defect_line_length(np.angle(psi))             # freshly-imprinted tangle at xi_turb (calibration ref)
    xinet0 = np.sqrt(N**3 / max(P0, 1))
    kick_at = set(int((j+1)*a.steps/(a.requench+1)) for j in range(a.requench)) if a.requench > 0 else set()
    hist = []
    t0 = time.time()
    for it in range(a.steps):
        if it in kick_at:                              # merger re-quench: blend toward a fresh xi_turb field
            fresh = correlated_unit_field(N, a.xi, rng)
            psi = (1-a.kick)*psi + a.kick*fresh
            psi = psi/(np.abs(psi)+1e-12); vel[:] = 0
        if a.overdamped:
            psi = psi + a.dt*(a.D*lap3(psi) + (1.0-np.abs(psi)**2)*psi/a.tau)
        else:
            accel = a.cdir**2*lap3(psi) + (1.0-np.abs(psi)**2)*psi/a.tau - a.gamma*vel
            vel = vel + a.dt*accel
            psi = psi + a.dt*vel
        if it % a.measure_every == 0 or it == a.steps-1:
            P = defect_line_length(np.angle(psi))
            hist.append((it, P))
    P_frozen = hist[-1][1]
    xinet = np.sqrt(N**3 / max(P_frozen, 1))
    c = xinet / xinet0                                 # calibration-free: frozen spacing / imprinted spacing
    xi_turb_over_Rcl = 0.47
    Mdef_over_MRAR = 1.0/(c*xi_turb_over_Rcl)**2
    mode = (f"OVERDAMPED(D={a.D})" if a.overdamped else f"INERTIAL(c_dir={a.cdir},gamma={a.gamma})")
    mode += f"  requench={a.requench}(kick={a.kick})" if a.requench else "  single-quench(freeze test)"
    lines = [
        f"defect_network_sim_3d  N={N}^3  xi_turb={a.xi}  {mode}  {'SMOKE' if a.smoke else 'RUN'}",
        f"imprinted tangle:  P0={P0}  xi_net0=sqrt(V/P0)={xinet0:.2f} cells   (freshly quenched at xi_turb)",
        f"frozen network:    P ={P_frozen}  xi_net =sqrt(V/P) ={xinet:.2f} cells",
        f"COEFFICIENT c = xi_net/xi_net0 = {c:.3f}   (c~1: freezes at injection scale; c>>1: coarsened away)",
        f"=> cluster M_def/M_RAR = (R_cl/(c*xi_turb))^2 at xi_turb=0.47 R_cl = {Mdef_over_MRAR:.2f}  (need ~4.6)",
        f"P history (it,P_line): {hist}",
        f"wall {time.time()-t0:.1f}s",
    ]
    print("\n".join(lines))
    with open(os.path.join(a.outdir, "defect_network_3d_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
