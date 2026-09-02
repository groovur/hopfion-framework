#!/usr/bin/env python3
r"""
Defect-network simulation (advection + c_dir->0): pin the scaling coefficient c = xi_net/xi_turb.
Refines cluster_defect_scaling.py (which ASSUMED xi_net~xi_turb). See notes/single_galaxy_winding_analysis.md.

MODEL: 2D complex order parameter psi = |psi| e^{i theta} (the director orientation as a phase; disclinations
= vortices, +-1 windings). Complex Ginzburg-Landau relaxation + a merger-like STIRRING FLOW (advection):
    d_t psi = D grad^2 psi + (1-|psi|^2) psi / tau  -  (v . grad) psi
  - D grad^2  : orientation relaxation -> drives opposite vortices together = ANNIHILATION. c_dir->0 = SMALL D
                (heavy/slow director: annihilation is slow, advection-dominated). Peclet = amp*xi_turb/D.
  - (1-|psi|^2)/tau : core term, keeps |psi|~1 away from defect cores (sets core size).
  - v : incompressible stirring at coherence scale xi_turb (the merger/turbulent injection), refreshed each
        eddy-turnover (mergers are transient). This CREATES vortex pairs (shear) and ADVECTS them.

GOAL: run to statistical steady state, count vortices N_v, get xi_net = L/sqrt(N_v), and the SCALING COEFFICIENT
    c = xi_net / xi_turb.  Then for a cluster: M_def/M_RAR = (R_cl/xi_net)^2 = (R_cl/(c*xi_turb))^2, with
    xi_turb ~ 0.47 R_cl (subcluster merger scale) -> M_def/M_RAR ~ 4.6/c^2. c~1 => recovers the deficit.

SETUP + SMOKE here (small grid, few steps). FULL parameter scan (xi_turb, Peclet, box size) = user runs it.
Uses numpy: /usr/local/Caskroom/miniconda/base/bin/python  (or any python3 with numpy).
"""
import argparse, os, time
import numpy as np

def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=128, help="grid size (box L = N cells = R_cl)")
    p.add_argument("--xi", type=float, default=16.0, help="stirring coherence scale xi_turb (cells)")
    p.add_argument("--D", type=float, default=0.15, help="orientation relaxation (small = c_dir->0)")
    p.add_argument("--amp", type=float, default=0.6, help="stirring velocity amplitude")
    p.add_argument("--tau", type=float, default=1.0, help="core relaxation time")
    p.add_argument("--dt", type=float, default=0.1, help="time step")
    p.add_argument("--steps", type=int, default=4000, help="total steps")
    p.add_argument("--measure-from", type=float, default=0.5, help="start averaging after this fraction of steps")
    p.add_argument("--inertial", action="store_true",
                   help="INERTIAL (damped-wave) director dynamics chi d_t^2 psi = K lap psi - gamma d_t psi + ... "
                        "-- the framework's actual c_dir->0 regime (heavy, frozen defects). Default is OVERDAMPED "
                        "(diffusive, WRONG for the framework -- coarsens).")
    p.add_argument("--cdir", type=float, default=0.1, help="inertial: wave speed c_dir=sqrt(K/chi); SMALL=heavy")
    p.add_argument("--gamma", type=float, default=0.02, help="inertial: damping (small=nearly lossless)")
    p.add_argument("--no-stir", action="store_true", help="disable stirring (pure coarsening test from a dense IC)")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    return p

def incompressible_flow(N, xi, amp, rng):
    """Random incompressible flow (v = curl of streamfunction) with spectral power peaked at k ~ 1/xi."""
    k = np.fft.fftfreq(N)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    Kmag = np.sqrt(KX**2 + KY**2)
    k0 = 1.0/xi
    spec = np.exp(-((Kmag - k0)**2)/(2*(0.5*k0)**2)); spec[0, 0] = 0.0
    stream = np.real(np.fft.ifft2(spec * np.exp(2j*np.pi*rng.random((N, N)))))
    vx = 0.5*(np.roll(stream, -1, 1) - np.roll(stream, 1, 1))    # d stream / dy
    vy = -0.5*(np.roll(stream, -1, 0) - np.roll(stream, 1, 0))   # -d stream / dx
    vmax = np.sqrt(vx*vx + vy*vy).max() + 1e-12
    return vx/vmax*amp, vy/vmax*amp

def lap(a):
    return (np.roll(a,1,0)+np.roll(a,-1,0)+np.roll(a,1,1)+np.roll(a,-1,1) - 4*a)

def grad(a):
    ax = 0.5*(np.roll(a,-1,0) - np.roll(a,1,0))
    ay = 0.5*(np.roll(a,-1,1) - np.roll(a,1,1))
    return ax, ay

def count_vortices(theta):
    def wrap(x): return (x + np.pi) % (2*np.pi) - np.pi
    d_right = wrap(np.roll(theta,-1,0) - theta)     # +x edge diffs
    d_up    = wrap(np.roll(theta,-1,1) - theta)     # +y edge diffs
    winding = d_right + np.roll(d_up,-1,0) - np.roll(d_right,-1,1) - d_up
    charge  = np.round(winding/(2*np.pi)).astype(int)
    return int(np.abs(charge).sum())

def main():
    a = build_parser().parse_args()
    if a.smoke:
        a.N, a.xi, a.steps = 48, 8.0, 400
    rng = np.random.default_rng(a.seed)
    N = a.N
    # dense random initial defect network (both modes) at scale ~xi (smoothed random phases)
    psi = np.exp(2j*np.pi*rng.random((N, N)))
    for _ in range(int(a.xi)):                          # smooth to imprint coherence ~xi
        psi = psi + 0.25*lap(psi); psi /= np.abs(psi)+1e-12
    vel = np.zeros_like(psi)                            # d_t psi (inertial only)
    vx, vy = incompressible_flow(N, a.xi, a.amp, rng)
    t_eddy = max(1, int((a.xi/max(a.amp,1e-6))/a.dt))
    peclet = a.amp*a.xi/a.D
    nmeas, Nv_acc = 0, 0
    hist = []
    t0 = time.time()
    for it in range(a.steps):
        if (not a.no_stir) and it % t_eddy == 0 and it > 0:
            vx, vy = incompressible_flow(N, a.xi, a.amp, rng)
        if a.no_stir:
            adv = 0.0
        else:
            gx, gy = grad(psi); adv = vx*gx + vy*gy
        if a.inertial:                                  # chi d_t^2 psi = K lap psi + core - gamma d_t psi - adv
            accel = a.cdir**2*lap(psi) + (1.0-np.abs(psi)**2)*psi/a.tau - a.gamma*vel - adv
            vel = vel + a.dt*accel
            psi = psi + a.dt*vel
        else:                                           # OVERDAMPED (diffusive) -- wrong for c_dir->0
            psi = psi + a.dt*(a.D*lap(psi) + (1.0-np.abs(psi)**2)*psi/a.tau - adv)
        if it >= a.measure_from*a.steps:
            Nv = count_vortices(np.angle(psi))
            Nv_acc += Nv; nmeas += 1
            if it % max(1, a.steps//10) == 0:
                hist.append((it, Nv))
    Nv_mean = Nv_acc/max(nmeas,1)
    xi_net = N/np.sqrt(max(Nv_mean, 1e-9))
    c = xi_net/a.xi
    # cluster application: M_def/M_RAR = (R_cl/(c xi_turb))^2 with xi_turb ~ 0.47 R_cl
    xi_turb_over_Rcl = 0.47
    Mdef_over_MRAR = 1.0/(c*xi_turb_over_Rcl)**2
    mode = (f"INERTIAL(c_dir={a.cdir},gamma={a.gamma})" if a.inertial else f"OVERDAMPED(D={a.D})") \
           + ("  NO-STIR(coarsening test)" if a.no_stir else f"  stir(xi={a.xi},amp={a.amp})")
    lines = [
        f"defect_network_sim  N={N}^2  {mode}  {'SMOKE' if a.smoke else 'RUN'}",
        f"steady-state vortices N_v = {Nv_mean:.1f}   xi_net = L/sqrt(N_v) = {xi_net:.2f} cells",
        f"SCALING COEFFICIENT c = xi_net/xi_turb = {c:.3f}   (c~1 confirms xi_net~xi_turb assumption)",
        f"=> cluster M_def/M_RAR = (R_cl/(c*xi_turb))^2 at xi_turb=0.47 R_cl = {Mdef_over_MRAR:.2f}  "
        f"(need ~4.6 to close deficit)",
        f"N_v history (it,Nv): {hist}",
        f"wall {time.time()-t0:.1f}s",
    ]
    print("\n".join(lines))
    with open(os.path.join(a.outdir, "defect_network_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    if a.smoke:
        print("\n[smoke] machinery check only (small box, few steps, poor statistics). For the real c, run e.g.:")
        print("  python defect_network_sim.py --N 256 --xi 16 --D 0.15 --steps 20000   (well-resolved, ~256 defects)")
        print("  then scan Peclet: --D 0.05 / 0.3 (c_dir->0 is small D / high Peclet) and xi=8,16,32 to check c is")
        print("  scale-independent. c(Peclet) is the pinned O(1) -> M_def/M_RAR.")

if __name__ == "__main__":
    main()
