#!/usr/bin/env python3
r"""
Phase 2c / A1: baryon-sourced director-strain halo(s) via topologically-protected
splay disclinations, with an ISOTROPIC-vs-ORIENTED comparison at minimal N.

Reuses the validated semi-Dirac Frank energy from src_paper16/gradient_flow_constrained.py:
    s4 = (1 - n_z^2)^2 = sin^4(theta)            (semi-Dirac directional suppression)
    J2a  = INT s4 * |grad n|^2 dV                (anisotropic splay term)
    J2iso= INT      |grad n|^2 dV
    K    = J2a + mu* * J2iso,   mu* = 3 - phi    (the Frank energy; Option-1 objective)
Option 1 (Frank-only): minimise K with the galaxy cores PINNED to a disclination
pattern, so the winding is held (a bare hedgehog would relax to uniform and give a
localised texture, NOT M_enc ~ r). Validation gate: single-galaxy M(r) must be ~ r
(log-log slope ~1). S^2-projected gradient flow + slerp angle-clamp (topology-safe),
as in the paper-16 solver.

SOURCES:
  - hedgehog (ISOTROPIC): n = (x-g)/|x-g|            -- 3D radial, orientation-free.
  - disclination (ORIENTED): escaped splay about spin axis a:
        n = cos(alpha) a + sin(alpha) e_perp,  alpha = arctan(rho_perp/rc),
        e_perp = d_perp/|d_perp|  (radial in the plane perpendicular to a).
    On-axis n||a; far-field radial-in-plane (splay). Axis a = galaxy spin axis.

MINIMAL COMPARISON (the point): N=2 in three configs -- isotropic, aligned (a1||a2),
perpendicular (a1 _|_ a2) -- to see if relative orientation changes the collective
strain. Averaging isotropically would wash this out (base file section 10).

Run examples (torch_intel):
  python phase2c_A1_disclination.py --mode calibrate --grid 64 --steps 2000
  python phase2c_A1_disclination.py --mode compare   --grid 64 --steps 2000
"""
import argparse, os, time, math
import numpy as np
import torch

PHI = (1 + 5**0.5) / 2
MU  = 3 - PHI

def parse():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["calibrate", "compare"], default="calibrate")
    p.add_argument("--grid", type=int, default=64)
    p.add_argument("--h", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=2000)
    p.add_argument("--lr", type=float, default=0.05)
    p.add_argument("--core", type=float, default=2.0, help="pinned core radius (cells)")
    p.add_argument("--sep", type=float, default=16.0, help="galaxy separation (cells) for N=2")
    p.add_argument("--delta-max-deg", type=float, default=8.0, help="angle clamp per step")
    p.add_argument("--oblate", type=float, default=0.3, help="disk axis-ratio q<1 for oriented sources")
    p.add_argument("--rgal", type=float, default=0.0,
                   help="galaxy screening radius (cells); >0 weights mass by S_eff~0 inside each galaxy (the "
                        "density-feedback 'shadow'). 0 = point sources (old behaviour).")
    p.add_argument("--aniso", action="store_true",
                   help="relax the PHYSICAL K_splay<0 + J4 energy (favours the extended hedgehog, opposes "
                        "Derrick collapse) instead of the positive Frank energy. The correct relaxation.")
    p.add_argument("--lam-j4", type=float, default=0.02, help="J4 Skyrme stabiliser weight (sets halo scale)")
    p.add_argument("--rc", type=float, default=0.0,
                   help="edge radius (cells): K_splay<0 driving cut off beyond R_c (the acceleration/screening "
                        "edge) -> stabilises M~r under FULL relaxation. 0 = no edge (runs away).")
    p.add_argument("--rc-width", type=float, default=3.0, help="edge transition width (cells)")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                    "phase2c_A1_output"))
    return p.parse_args()

def grid_coords(N, h, dev):
    cv = h * (torch.arange(N, device=dev) - N//2 + 0.5)
    X, Y, Z = torch.meshgrid(cv, cv, cv, indexing="ij")
    return torch.stack([X, Y, Z], dim=-1)   # (N,N,N,3)

def hedgehog_ic(coords, g):
    d = coords - g
    return d / d.norm(dim=-1, keepdim=True).clamp(1e-6)

def oriented_hedgehog_ic(coords, g, axis, q):
    """Anisotropic (oblate-disk) hedgehog: n = normalize(dir of grad rho) for an
       oblate mass with symmetry axis `axis` and axis-ratio q = sigma_par/sigma_perp.
       q=1 -> spherical hedgehog (isotropic); q<1 -> disk: axial gradient enhanced.
       Recovers M_enc ~ r (radial splay) while carrying the galaxy orientation."""
    a = axis / axis.norm()
    d = coords - g
    d_par = (d * a).sum(-1, keepdim=True) * a
    d_perp = d - d_par
    e_dir = d_perp + d_par / (q*q)                    # grad rho direction (oblate)
    return e_dir / e_dir.norm(dim=-1, keepdim=True).clamp(1e-6)

def frank_energy(n, h):
    """K = J2a + mu* J2iso  (reused from gradient_flow_constrained E_geom); also energy density."""
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    s4 = (1 - nz**2).clamp(0, 1)**2
    def cd(u, a): return (torch.roll(u, -1, a) - torch.roll(u, 1, a)) / (2*h)
    comps = []
    for c in (nx, ny, nz):
        comps += [cd(c, 0), cd(c, 1), cd(c, 2)]
    g2 = sum(x**2 for x in comps)      # |grad n|^2  -> POLARON MASS density (chi-free proxy)
    u_a  = s4 * g2                     # anisotropic (sin^4) energy density
    dens = u_a + MU * g2               # total Frank energy density (signed sense via anisotropy)
    J2a  = u_a.sum() * h**3
    J2iso = g2.sum() * h**3
    K = J2a + MU * J2iso
    # polaron mass proxy = INT |grad n|^2 (chi/3 factor omitted; cancels in ratios)
    return K, dens, g2

# measured semi-Dirac Frank stiffnesses (frank_run.py test): splay UNSTABLE, bend stable
K_SPLAY = -8.45
K_BEND  =  0.534

def aniso_energy(n, h, lam_j4, splay_wt=None):
    """PHYSICAL Frank energy with K_splay<0 (favours the extended radial splay -> opposes Derrick
       collapse). splay_wt (field, 1 inside R_c -> 0 outside) = the EDGE term: the K_splay<0 driving
       acts only within the acceleration/screening edge R_c, so the halo extends to R_c and STOPS
       (stabilises M~r) instead of running away. + lam_j4 * J4 (mild Skyrme regulariser)."""
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    def cd(u, a): return (torch.roll(u, -1, a) - torch.roll(u, 1, a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    div = nxx + nyy + nzz                                   # splay = div n
    cx, cy, cz = nzy-nyz, nxz-nzx, nyx-nxy                  # curl n
    bx, by, bz = ny*cz-nz*cy, nz*cx-nx*cz, nx*cy-ny*cx      # bend = n x curl n
    bend2 = bx*bx + by*by + bz*bz
    ksp = K_SPLAY if splay_wt is None else K_SPLAY*splay_wt   # edge-cut splay driving
    E_frank = 0.5*(ksp*div*div + K_BEND*bend2).sum()*h**3
    Fxy = nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    J4 = (Fxy*Fxy+Fxz*Fxz+Fyz*Fyz).sum()*h**3
    return E_frank + lam_j4*J4

def frank_decompose(n, h):
    """Return (INT splay^2, INT twist^2, INT bend^2) densities as spatial fields, so the
       ANISOTROPIC polaron mass = chi_s*INT splay^2 + chi_t*INT twist^2 + chi_b*INT bend^2
       (isotropic chi -> sum = INT|grad n|^2). Tests the 'isotropic injustice'."""
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    def cd(u, a): return (torch.roll(u, -1, a) - torch.roll(u, 1, a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz = cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz = cd(nz,0),cd(nz,1),cd(nz,2)
    div = nxx + nyy + nzz
    cx, cy, cz = nzy-nyz, nxz-nzx, nyx-nxy
    twist = nx*cx + ny*cy + nz*cz
    bx, by, bz = ny*cz-nz*cy, nz*cx-nx*cz, nx*cy-ny*cx
    splay2 = div*div
    twist2 = twist*twist
    bend2  = bx*bx + by*by + bz*bz
    return splay2, twist2, bend2

def build_field(coords, gals, rc, dev):
    """gals: list of (pos(3), axis(3), q).  q=1 -> isotropic hedgehog; q<1 -> oriented disk.
       Returns initial n and a pinned-core mask + values."""
    N = coords.shape[0]
    def ic(gt, at, q):
        return oriented_hedgehog_ic(coords, gt, at, q)
    n = torch.zeros_like(coords)
    for (g, a, q) in gals:
        gt = torch.tensor(g, dtype=torch.float32, device=dev)
        at = torch.tensor(a, dtype=torch.float32, device=dev)
        w = 1.0 / ((coords - gt).norm(dim=-1, keepdim=True).clamp(1e-3))**2   # nearest dominates
        n = n + w * ic(gt, at, q)
    n = n / n.norm(dim=-1, keepdim=True).clamp(1e-6)
    mask = torch.zeros(N, N, N, dtype=torch.bool, device=dev)
    vals = n.clone()
    for (g, a, q) in gals:
        gt = torch.tensor(g, dtype=torch.float32, device=dev)
        at = torch.tensor(a, dtype=torch.float32, device=dev)
        r = (coords - gt).norm(dim=-1)
        m = r <= rc
        mask |= m
        vals = torch.where(m.unsqueeze(-1), ic(gt, at, q), vals)
    # OUTER BOUNDARY: pin the edge shell to radial-outward from the galaxy centroid
    # (the charge-N far field), so the extended M~r halo is not contaminated by the
    # periodic wrap. Without this the profile steepens (boundary artefact).
    gc = torch.tensor(np.mean([g for g, _, _ in gals], axis=0), dtype=torch.float32, device=dev)
    far = (coords - gc); far = far / far.norm(dim=-1, keepdim=True).clamp(1e-6)
    bnd = torch.zeros(N, N, N, dtype=torch.bool, device=dev)
    bnd[0,:,:] = bnd[-1,:,:] = bnd[:,0,:] = bnd[:,-1,:] = bnd[:,:,0] = bnd[:,:,-1] = True
    vals = torch.where(bnd.unsqueeze(-1), far, vals)
    mask |= bnd
    n = torch.where(mask.unsqueeze(-1), vals, n)
    return n, mask, vals

def slerp_clamp(n_old, n_new, dmax):
    cos_a = (n_old * n_new).sum(-1, keepdim=True).clamp(-1+1e-6, 1-1e-6)
    ang = torch.acos(cos_a)
    too = (ang > dmax).squeeze(-1)
    if not too.any():
        return n_new
    sa = torch.sin(ang).clamp(1e-8)
    t = (dmax / ang.clamp(1e-8))
    sl = (torch.sin((1-t)*ang)/sa) * n_old + (torch.sin(t*ang)/sa) * n_new
    sl = sl / sl.norm(dim=-1, keepdim=True).clamp(1e-10)
    out = n_new.clone(); out[too] = sl[too]
    return out / out.norm(dim=-1, keepdim=True).clamp(1e-10)

def relax(n0, mask, vals, h, steps, lr, dmax, efn):
    n = n0.clone().requires_grad_(True)
    hist = []
    for it in range(steps):
        E = efn(n)
        (grad,) = torch.autograd.grad(E, n)
        with torch.no_grad():
            grad = grad - (grad*n).sum(-1, keepdim=True)*n     # tangent projection
            n_new = n - lr*grad
            n_new = n_new / n_new.norm(dim=-1, keepdim=True).clamp(1e-10)
            n_new = slerp_clamp(n, n_new, dmax)
            n_new[mask] = vals[mask]                            # re-pin cores
        n = n_new.detach().requires_grad_(True)
        if it % max(1, steps//10) == 0 or it == steps-1:
            with torch.no_grad():
                hist.append(float(efn(n)))
    return n.detach(), hist

def mass_profile(dens, coords, center, rmax, nbin=24):
    r = (coords - center).norm(dim=-1)
    edges = torch.linspace(0, rmax, nbin+1, device=dens.device)
    M = []
    for e in edges[1:]:
        M.append(float(dens[r <= e].sum()))
    return edges[1:].cpu().numpy(), np.array(M)

def loglog_slope(r, M):
    ok = (M > 0) & (r > 0)
    if ok.sum() < 3: return float("nan")
    lr, lM = np.log(r[ok]), np.log(M[ok])
    hi = lr > np.median(lr)     # outer half (avoid core)
    if hi.sum() < 2: hi = ok[ok]
    return float(np.polyfit(lr[hi], lM[hi], 1)[0])

def main():
    a = parse()
    if a.smoke: a.grid, a.steps = 24, 300
    os.makedirs(a.outdir, exist_ok=True)
    dev = torch.device("cpu")
    N, h = a.grid, a.h
    coords = grid_coords(N, h, dev)
    c = N/2 * h
    center = torch.tensor([0., 0., 0.], device=dev)  # box centred at origin
    dmax = math.radians(a.delta_max_deg)
    rmax = 0.45 * N * h
    if a.aniso:
        if a.rc > 0:
            rr = coords.norm(dim=-1)                       # distance from box centre (galaxy at origin)
            splay_wt = 0.5*(1.0 - torch.tanh((rr - a.rc*h)/(a.rc_width*h)))   # 1 inside R_c -> 0 outside
            elabel = f"PHYSICAL aniso K_splay<0 (edge R_c={a.rc}) + {a.lam_j4}*J4"
        else:
            splay_wt = None
            elabel = f"PHYSICAL aniso K_splay<0 (no edge -> runaway) + {a.lam_j4}*J4"
        efn = lambda nn: aniso_energy(nn, h, a.lam_j4, splay_wt)
    else:
        efn = lambda nn: frank_energy(nn, h)[0]
        elabel = "positive-Frank (Derrick-collapses)"
    lines = [f"Phase 2c/A1 disclination halo  mode={a.mode}  grid={N}^3 h={h} steps={a.steps} "
             f"core={a.core}  energy={elabel}"]
    print(f"[A1] energy = {elabel}")

    def run(gals, tag):
        n0, mask, vals = build_field(coords, gals, a.core, dev)
        t0 = time.time()
        n, hist = relax(n0, mask, vals, h, a.steps, a.lr, dmax, efn)
        K, dens, g2 = frank_energy(n, h)
        g2 = g2.detach()
        # density-feedback SCREENING ('shadow'): S_eff ~ 0 inside each galaxy (r<rgal), ->1 outside
        if a.rgal > 0:
            w = torch.ones_like(g2)
            for (gpos, _, _) in gals:
                gt = torch.tensor(gpos, dtype=torch.float32, device=dev)
                ri = (coords - gt).norm(dim=-1)
                w = w * (1.0 - torch.exp(-(ri/a.rgal)**2))     # 0 at galaxy centre -> 1 far
            g2_eff = g2 * w
        else:
            g2_eff = g2
        gc = torch.tensor(np.mean([g for g, _, _ in gals], axis=0), dtype=torch.float32, device=dev)
        r, M = mass_profile(g2_eff, coords, gc, rmax)         # screened polaron mass profile
        slope = loglog_slope(r, M)
        m_tot = float(M[-1])                                  # polaron mass within Rcl (screened if rgal>0)
        # splay/twist/bend decomposition (screened), for the anisotropic (chi_b/chi_s) mass
        sp2, tw2, bd2 = frank_decompose(n, h)
        if a.rgal > 0:
            sp2, tw2, bd2 = sp2*w, tw2*w, bd2*w
        S = float(sp2.sum()*h**3); T = float(tw2.sum()*h**3); B = float(bd2.sum()*h**3)
        Sigma = g2_eff.sum(dim=2)
        mono = all(hist[i+1] <= hist[i] + 1e-6*abs(hist[i]) for i in range(len(hist)-1))
        np.save(os.path.join(a.outdir, f"sigma_{tag}.npy"), Sigma.cpu().numpy())
        line = (f"[{tag:16s}] m={m_tot:.4e} slope={slope:.2f}  splay={S:.3e} twist={T:.3e} bend={B:.3e} "
                f" t={time.time()-t0:.1f}s")
        print(line); lines.append(line)
        return dict(m=m_tot, M=M, slope=slope, S=S, T=T, B=B)

    z = [0., 0., 1.]; x = [1., 0., 0.]; q = a.oblate
    if a.mode == "calibrate":
        run([([0.,0.,0.], z, 1.0)], "single_isotropic")       # spherical hedgehog -> want slope ~1
        lines.append("GATE: single-galaxy M(r)=INT|grad n|^2 within r should scale ~r (slope ~1). "
                     "This is the chi-free polaron mass profile; if slope~1 the config is good.")
    else:  # 'compare' -> chi-FREE collective enhancement (isotropic hedgehogs, same orientation)
        s = a.sep * h
        g1 = [-s/2, 0., 0.]; g2p = [+s/2, 0., 0.]
        r1 = run([([0.,0.,0.], z, 1.0)],            "N1_single_iso")
        r2 = run([(g1, z, 1.0), (g2p, z, 1.0)],     "N2_isotropic")
        eta_iso = r2['m']/(2.0*r1['m']) if r1['m']>0 else float('nan')
        def eta_rho(rho):   # ANISOTROPIC mass = (splay+twist) + rho*bend ; rho = chi_bend/chi_splay
            m1 = (r1['S']+r1['T']) + rho*r1['B']; m2 = (r2['S']+r2['T']) + rho*r2['B']
            return m2/(2.0*m1) if m1>0 else float('nan')
        num = 2*(r1['S']+r1['T']) - (r2['S']+r2['T']); den = r2['B'] - 2*r1['B']
        rho_star = num/den if abs(den) > 1e-30 else float('inf')
        lines.append(f"ISOTROPIC-chi collective enhancement eta_iso = {eta_iso:.3f}  (<1 = screening)")
        lines.append(f"  single: splay={r1['S']:.3e} twist={r1['T']:.3e} bend={r1['B']:.3e}")
        lines.append(f"  N2    : splay={r2['S']:.3e} twist={r2['T']:.3e} bend={r2['B']:.3e}")
        lines.append("ANISOTROPIC eta(chi_bend/chi_splay) -- does bend-weighting reverse the screening?")
        for rho in (1, 2, 3, 5, 10, 20, 50):
            lines.append(f"  chi_b/chi_s={rho:>3} : eta={eta_rho(rho):.3f}")
        lines.append(f"REVERSAL THRESHOLD chi_b/chi_s (eta=1) = {rho_star:.2f}  "
                     f"(frank_chi tentative chi ratio ~2.7). Small+reachable -> isotropic WAS an injustice; "
                     f"huge/negative -> anisotropy cannot rescue the screening.")

    with open(os.path.join(a.outdir, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[A1] wrote {a.outdir}/summary.txt (+ sigma_*.npy)")

if __name__ == "__main__":
    main()
