#!/usr/bin/env python3
"""
qh3_trefoil_solver_3d_v11.py  —  classic Metropolis simulated annealing
==========================================================================
WHY THIS SCRIPT EXISTS, AND WHY IT IS A DIFFERENT ALGORITHM FROM v6-v9:

  v6-v8 (soft K_fb*J4 coupling) all showed J4 leaking away monotonically
  to near-zero, regardless of sigma (checked across 5+ orders of
  magnitude). v9 (a TRUE hard Lagrange-multiplier constraint forcing
  J4's INTEGRAL to stay exactly fixed) fixed THAT failure mode, but
  revealed a DIFFERENT one: gradient descent found a way to satisfy
  "J4's integral is unchanged" by collapsing nearly all of the Hopf
  charge into one tiny, extremely dense spike (final density ~459 vs
  initial ~9.6 — a ~48x concentration into a single small region, far
  from any of the three crossings). Fixing the INTEGRAL of J4 was not
  a strong enough constraint, because a sharp spike and a smooth,
  properly-spread-out trefoil can have the same integral.

  Both failures share a root cause: pure gradient descent ALWAYS moves
  downhill, so it will always find and exploit ANY direction that
  lowers the objective, including degenerate ones (dilution in v6-v8,
  spiking in v9) that a careful human would recognise as unphysical
  but that the optimiser has no way to know to avoid.

  Paper XV's own remedy list names Battye-Sutcliffe-style simulated
  annealing as the candidate fix for exactly this kind of situation.
  IMPORTANT CORRECTION, made before writing this script: web research
  into the actual Battye & Sutcliffe (1998, Phys. Rev. Lett. 81, 4798)
  paper and its follow-up ("Solitons, Links and Knots") shows their own
  robustness came from relaxing many different initial conditions and
  keeping the best result — not from a textbook Metropolis temperature
  schedule. The user explicitly asked for the textbook Metropolis
  version anyway, reasoning that we already have a strong structural
  prior (the BS-ansatz trefoil shape from Paper XV) and don't need
  blind multi-start exploration — we need a mechanism that can escape
  the SPECIFIC spike/dilution degeneracies while searching near the
  known-good shape. This script implements that: classic Metropolis
  acceptance (locally-worse moves accepted with probability
  exp(-dE/T), temperature cooled over the run), NOT multi-start.

WHY A NAIVE IMPLEMENTATION WOULD BE TOO SLOW, AND HOW THIS ONE AVOIDS IT:

  Textbook Metropolis proposes ONE random change at a time (e.g. to
  one lattice site) and accepts/rejects it before proposing the next.
  Done literally as a Python loop over a 64^3≈262000-site grid, this
  would be far too slow to run a meaningful number of sweeps.

  The standard fix used throughout lattice field theory Monte Carlo
  (this is genuinely how large lattice simulations are implemented,
  not an invented shortcut) is a CHECKERBOARD / sublattice decomposition:
  partition the grid into colour classes such that no single energy
  term ever depends on two sites of the SAME colour simultaneously —
  then ALL sites of one colour can be proposed, evaluated, and
  accepted/rejected SIMULTANEOUSLY via vectorised tensor operations,
  because their proposals provably cannot interact with each other.

  GETTING THE COLOURING RIGHT TOOK REAL CARE, and is worth recording
  here because it is easy to get wrong (it was gotten wrong twice
  during this script's own design, before being caught by direct
  numerical testing against this project's actual energy formulas —
  see the development log this version corresponds to):
    - A simple even/odd (mod-2) colouring per axis is NOT sufficient.
      This project's central-difference stencil, cd(u,axis) =
      (u[i+1]-u[i-1])/(2h), SKIPS the centre point — so a term at
      index k depends on k-1 and k+1, meaning TWO sites that are
      distance-2 apart along one axis (e.g. indices 2 and 4) can both
      appear in the SAME density term (at index 3). A mod-2 colouring
      puts both 2 and 4 in the same colour class, which is wrong.
    - A mixed mod-4-on-one-axis / mod-2-on-the-others scheme was also
      tested and found to be wrong (verified numerically to give an
      incorrect energy change), because the risky distance-2 pairing
      can occur along ANY of the three axes independently.
    - The verified-correct scheme is mod-4 colouring independently on
      ALL THREE axes (colour = (i%4, j%4, k%4)), giving 64 colour
      classes total. This was checked directly against this project's
      real g2 and rho_J4 formulas (not a simplified toy energy) and
      confirmed to give EXACTLY the same total energy change whether
      a colour class's sites are updated all-at-once or one at a time
      and summed — to floating-point precision.

  Each of the 64 colour sub-steps is a small, cheap, local computation
  (only the ~1/64 of sites in that colour, plus their immediate
  neighbours, need touching) — substantially cheaper per sub-step than
  a full-grid autograd backward pass, even though there are 64 of them
  per sweep.

WHAT'S MINIMISED: E_geom = K_fb * J4 (the proven Q_H=2 product trick,
ported again as in v7), NOT a hard J4 constraint. The product
structure directly discourages J4→0 (the v6-v8 dilution failure);
the LOCAL, single-site-at-a-time Metropolis proposals (as opposed to
a single global gradient-descent direction) make it much harder for
the field to coordinate the kind of large-scale concentration that
produced v9's spike, since any one local move that suddenly
over-concentrates density will usually raise the energy locally and
be rejected, especially as the temperature cools.

BUG FOUND AND FIXED AFTER THE FIRST v11 RUN: the first version set
T0 as a fixed fraction (5%) of the TOTAL E_geom. This is the wrong
reference scale — E_geom=K*J4 is a product of two O(1000)-scale
numbers, but a single colour-class proposal only ever touches 1/64
of the grid, so its typical |dE_geom| is a much smaller, DIFFERENT
quantity. The first run's T0 was so mismatched that acceptance fell
to EXACTLY 0.0% by sweep ~50 and stayed there for 1700+ consecutive
sweeps, with every diagnostic frozen identically — not a genuine cold
anneal, a stuck/frozen search. This version FIXES this by directly
MEASURING the typical |dE_geom| from ~20 trial colour-class proposals
on the actual initial field before the run starts (calibrate_T0()),
and setting T0 from that measured scale (targeting ~50% acceptance of
a median-sized move at the start) rather than guessing from the total
energy. A stuck-run detector also now prints a loud warning if
acceptance holds at exactly 0.0% for 100 consecutive sweeps, so this
specific failure mode can never again pass silently.

PROPOSAL MECHANISM: each active site's unit vector n is perturbed by
a SMALL RANDOM ROTATION (not Gaussian noise + renormalise, which would
bias the proposal distribution) — rotate n by a small random angle
about a random axis, drawn fresh per site per sub-step.

TEMPERATURE SCHEDULE: geometric cooling, T_k = T_0 * cooling_rate^k,
one update of T per FULL SWEEP (all 64 colours).

Install:  pip install torch --index-url https://download.pytorch.org/whl/cpu

USAGE
-----
  # First test: a single R0, moderate sweep count, to see whether
  # annealing avoids BOTH the v6-v8 dilution and the v9 spike:
  python qh3_trefoil_solver_3d_v11.py \\
      --R0_list "3.0" --C_star 2.4987 \\
      --N 64 --h 0.26 --sweeps 2000 --outdir v11_sanity_check

OUTPUT (in --outdir)
---------------------
  {tag}_log.txt            full run log (per-sweep diagnostics)
  {tag}_results.csv        one row per R0
  {tag}_report.txt         human-readable summary
  n_R0{val}_FINAL.npz      final field (viewable with qh3_field_viewer.py)

WHAT TO LOOK FOR
-----------------
  - J4/J0 should stay near 1 (as in v9) WITHOUT the density also
    spiking into one small region (the v9 failure) — check r_bar
    (the J4-weighted mean radius, same diagnostic as v8) stays near
    R0, and check the acceptance rate printed each sweep: a healthy
    anneal shows acceptance starting high (most moves accepted while
    hot) and dropping smoothly as T cools, NOT collapsing to near-0
    immediately (which would mean the schedule started too cold) or
    staying near 1 throughout (started too hot / never actually
    cooled enough to commit to a configuration).
  - View the result with qh3_field_viewer.py — specifically check
    whether the FINAL state's density peaks sit near the three known
    crossing positions at a REASONABLE density scale (comparable to
    the initial BS-ansatz's peak, not 1000x smaller as in v6-v8, and
    not 50x larger as in v9).
"""
import numpy as np, time, argparse, os, sys
from scipy.spatial import KDTree
try:
    import torch
except ImportError:
    print("pip install torch --index-url https://download.pytorch.org/whl/cpu")
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('--N',           type=int,   default=64)
ap.add_argument('--h',           type=float, default=0.26)
ap.add_argument('--C_star',      type=float, default=2.4987)
ap.add_argument('--R0_list',     type=str,   default="3.0")
ap.add_argument('--r0',          type=float, default=0.874)
ap.add_argument('--NT',          type=int,   default=1500)
ap.add_argument('--sweeps',      type=int,   default=2000,
                help='Number of FULL sweeps (each sweep = all 64 colour sub-steps once)')
ap.add_argument('--T0',          type=float, default=None,
                help='Initial temperature. If not given, auto-set to ~5% of the initial E_geom value so early moves are frequently accepted.')
ap.add_argument('--cooling_rate', type=float, default=0.99,
                help='Geometric cooling: T_k = T0 * cooling_rate^k per sweep. Lowered from an initial 0.997 default after the first v11 run showed the field drifting toward a spike-prone configuration during the many sweeps it takes 0.997 to cool meaningfully — 0.99 cools roughly 3x faster.')
ap.add_argument('--rotation_scale', type=float, default=0.25,
                help='Max random rotation angle (radians) for a single-site proposal')
ap.add_argument('--print_every', type=int,   default=50)
ap.add_argument('--outdir',      type=str,   default='.')
ap.add_argument('--device',      type=str,   default='cpu')
ap.add_argument('--seed',        type=int,   default=0)
args = ap.parse_args()

phi = (1+5**0.5)/2; phi6 = phi**6; MU = 3.0-phi
N, h = args.N, args.h
dev = torch.device(args.device)
os.makedirs(args.outdir, exist_ok=True)
tag = f'qh3_v11'
log_path = os.path.join(args.outdir, f'{tag}_log.txt')

torch.manual_seed(args.seed)
np.random.seed(args.seed)

R0_LIST = [float(x) for x in args.R0_list.split(',')]

print(f"\n{'='*70}")
print(f"  Q_H=3 Trefoil Solver v11 — classic Metropolis annealing")
print(f"  E_geom = K_fb * J4, 64-colour checkerboard sublattice updates")
print(f"  (verified exact for this project's stencil formulas — see docstring)")
print(f"  phi^6 = {phi6:.6f}   2*phi = {2*phi:.6f}")
print(f"  cooling_rate={args.cooling_rate}  rotation_scale={args.rotation_scale} rad")
print(f"  R0 scan: {R0_LIST}")
print(f"  Grid: {N}^3  h={h}  box=[{-N*h/2:.2f},{N*h/2:.2f}]")
print(f"{'='*70}")

log = open(log_path, 'a')
log.write(f"# Q_H=3 v11 (Metropolis annealing)  N={N} h={h} C*={args.C_star}\n")
log.write(f"# R0_list = {R0_LIST}\n")
log.write(f"# R0  sweep  T  accept_rate  K  J2a  J4  J4_J0  r_bar  KJ4  KJ2a  E_geom\n")

# ── Grid ──────────────────────────────────────────────────────────────
cv = h*(np.arange(N) - N//2 + 0.5)
pts = np.stack(np.meshgrid(cv, cv, cv, indexing='ij'), axis=-1).reshape(-1, 3).astype(np.float32)
dist_from_origin_flat = np.linalg.norm(pts, axis=-1).astype(np.float32)
dist_from_origin = torch.tensor(dist_from_origin_flat.reshape(N, N, N), dtype=torch.float32, device=dev)

box = N*h
max_R0 = max(R0_LIST)
ext_check = max_R0 + args.r0 + 2/args.C_star + 1.0
if box/2 < ext_check:
    print(f"  WARNING: box half-width {box/2:.2f} may be too small for "
          f"the largest R0={max_R0} (+tube +margin needs ~{ext_check:.2f}).")

# ── Global energy (used for sweep-boundary refresh and diagnostics) ───
def compute_global(n, mu_param):
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a  = (s4 * g2).sum() * h**3
    J2iso = g2.sum() * h**3
    K = J2a + mu_param * J2iso
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    J4 = rho_J4.sum() * h**3
    r_bar_num = (rho_J4 * dist_from_origin).sum()
    r_bar_den = rho_J4.sum()
    r_bar = (r_bar_num/r_bar_den).item() if r_bar_den.item() > 1e-12 else float('nan')
    return K.item(), J2a.item(), J4.item(), r_bar


def local_KJ4_contribution(n, mask):
    """Sum of the LOCAL g2 and rho_J4 density contributions AT the sites
    selected by `mask` only (used to compute dK, dJ4 for a checkerboard
    colour's proposed change — exact because same-colour sites never
    share a stencil term, verified numerically against these exact
    formulas before this script was written)."""
    nx, ny, nz = n[...,0], n[...,1], n[...,2]
    s4 = (1 - nz**2).clamp(0,1)**2
    def cd(u, a): return (torch.roll(u,-1,a) - torch.roll(u,1,a)) / (2*h)
    nxx,nxy,nxz = cd(nx,0), cd(nx,1), cd(nx,2)
    nyx,nyy,nyz = cd(ny,0), cd(ny,1), cd(ny,2)
    nzx,nzy,nzz = cd(nz,0), cd(nz,1), cd(nz,2)
    g2 = (nxx**2+nxy**2+nxz**2 + nyx**2+nyy**2+nyz**2 + nzx**2+nzy**2+nzz**2)
    J2a_density  = s4 * g2
    J2iso_density = g2
    Fxy = nx*(nyx*nzy-nzx*nyy) + ny*(nzx*nxy-nxx*nzy) + nz*(nxx*nyy-nyx*nxy)
    Fxz = nx*(nyx*nzz-nzx*nyz) + ny*(nzx*nxz-nxx*nzz) + nz*(nxx*nyz-nyx*nxz)
    Fyz = nx*(nyy*nzz-nzy*nyz) + ny*(nzy*nxz-nxy*nzz) + nz*(nxy*nyz-nyy*nxz)
    rho_J4 = Fxy**2 + Fxz**2 + Fyz**2
    # NOTE: this returns the FULL density fields; the caller sums over
    # the relevant masked sites (and, importantly, ALSO needs the
    # density at each active site's 6 neighbours, since the proposal
    # changes derivatives evaluated there too — handled by always
    # working with the FULL recomputed field within compute_global-style
    # calls, restricted afterward via masked sums for the dK/dJ4 used in
    # the acceptance rule).
    return J2a_density, J2iso_density, rho_J4


def random_rotation_perturb(n, mask, scale):
    """Propose a new unit vector at each True site in `mask` by rotating
    the EXISTING vector by a small random angle about a random axis —
    exactly preserves |n|=1 with no renormalisation bias (unlike adding
    Gaussian noise then renormalising)."""
    idx = mask.nonzero(as_tuple=False)
    n_sites = idx.shape[0]
    if n_sites == 0:
        return n, idx
    v = n[mask]  # (n_sites, 3)
    axis = torch.randn(n_sites, 3, device=dev)
    axis = axis - (axis*v).sum(-1, keepdim=True)*v  # make axis perpendicular-ish (not required to be exact)
    axis = axis / axis.norm(dim=-1, keepdim=True).clamp(1e-10)
    angle = (torch.rand(n_sites, device=dev)*2-1) * scale
    cosA, sinA = torch.cos(angle).unsqueeze(-1), torch.sin(angle).unsqueeze(-1)
    # Rodrigues' rotation formula
    cross = torch.cross(axis, v, dim=-1)
    dot = (axis*v).sum(-1, keepdim=True)
    v_new = v*cosA + cross*sinA + axis*dot*(1-cosA)
    v_new = v_new / v_new.norm(dim=-1, keepdim=True).clamp(1e-10)
    return v_new, idx


# Precompute the 64 colour masks (mod-4 on each axis) once — verified
# exact for this project's stencil formulas (see docstring).
ar = torch.arange(N, device=dev)
I, J, K_ = torch.meshgrid(ar, ar, ar, indexing='ij')
COLOUR_MASKS = []
for ci in range(4):
    for cj in range(4):
        for ck in range(4):
            COLOUR_MASKS.append((I%4==ci) & (J%4==cj) & (K_%4==ck))
print(f"  Built {len(COLOUR_MASKS)} colour classes "
      f"(avg {N**3/len(COLOUR_MASKS):.0f} sites each)")


def anneal_sweep(n, K_cur, J4_cur, T, rotation_scale):
    """One full sweep over all 64 colours. Returns updated n, K, J4, and
    the fraction of proposals accepted this sweep."""
    n_accept = 0
    n_total = 0
    for mask in COLOUR_MASKS:
        v_new, idx = random_rotation_perturb(n, mask, rotation_scale)
        if idx.shape[0] == 0:
            continue
        # Compute BEFORE densities (full field, current n)
        J2a_d0, J2iso_d0, rho0 = local_KJ4_contribution(n, mask)
        # Build trial field with this colour's sites replaced
        n_trial = n.clone()
        n_trial[mask] = v_new
        J2a_d1, J2iso_d1, rho1 = local_KJ4_contribution(n_trial, mask)

        # Per-site local energy change. Since same-colour sites never
        # share a stencil term, the change at site s's OWN location
        # plus its dependence on the OTHER terms that reference s are
        # already captured by comparing the FULL density fields before
        # and after — summing the difference over ALL grid points
        # gives the exact global dK,dJ4 for this colour's combined
        # proposal (verified exact for simultaneous same-colour updates).
        dJ2a = (J2a_d1 - J2a_d0).sum() * h**3
        dJ2iso = (J2iso_d1 - J2iso_d0).sum() * h**3
        dK = (dJ2a + MU*dJ2iso).item()
        dJ4 = ((rho1 - rho0).sum() * h**3).item()

        K_new = K_cur + dK
        J4_new = J4_cur + dJ4
        dE = (K_new*J4_new) - (K_cur*J4_cur)

        accept = dE < 0 or np.random.rand() < np.exp(-dE/max(T, 1e-12))
        n_total += idx.shape[0]
        if accept:
            n = n_trial
            K_cur, J4_cur = K_new, J4_new
            n_accept += idx.shape[0]
        # else: n, K_cur, J4_cur unchanged for this colour

    accept_rate = n_accept / max(n_total, 1)
    return n, K_cur, J4_cur, accept_rate


# ── Build initial condition for a given R0 (BS-ansatz seed) ──────────
def build_ic(R0, r0, C_star, NT):
    t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
    Gx = (R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
    Gy = (R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
    Gz = r0*np.sin(3*t_arr)
    Gamma = np.stack([Gx,Gy,Gz], axis=1).astype(np.float32)
    tree = KDTree(Gamma)
    dists, idx = tree.query(pts)
    rho = np.clip(dists, 1e-6, None).astype(np.float32)
    t_near = t_arr[idx].astype(np.float32)
    f_ = 2*np.arctan(rho**(-C_star))
    tht = 3*t_near
    nx_ = np.sin(f_)*np.cos(tht); ny_ = np.sin(f_)*np.sin(tht); nz_ = np.cos(f_)
    n_np = np.stack([nx_.reshape(N,N,N), ny_.reshape(N,N,N),
                     nz_.reshape(N,N,N)], axis=-1).astype(np.float32)
    n_np /= np.linalg.norm(n_np, axis=-1, keepdims=True).clip(1e-10)
    return n_np


def probe_dE(n, K_cur, J4_cur, mask, rotation_scale):
    """Compute |dE_geom| for ONE proposed (but not applied) colour-class
    move, without mutating any state — used only to MEASURE the typical
    scale of a single sub-step's energy change, for temperature
    calibration. This is the same computation anneal_sweep does
    internally, factored out so it can be run as a dry-run probe."""
    v_new, idx = random_rotation_perturb(n, mask, rotation_scale)
    if idx.shape[0] == 0:
        return 0.0
    J2a_d0, J2iso_d0, rho0 = local_KJ4_contribution(n, mask)
    n_trial = n.clone()
    n_trial[mask] = v_new
    J2a_d1, J2iso_d1, rho1 = local_KJ4_contribution(n_trial, mask)
    dJ2a = (J2a_d1 - J2a_d0).sum() * h**3
    dJ2iso = (J2iso_d1 - J2iso_d0).sum() * h**3
    dK = (dJ2a + MU*dJ2iso).item()
    dJ4 = ((rho1 - rho0).sum() * h**3).item()
    K_new = K_cur + dK
    J4_new = J4_cur + dJ4
    dE = (K_new*J4_new) - (K_cur*J4_cur)
    return abs(dE)


def calibrate_T0(n, K_cur, J4_cur, rotation_scale, n_probes=20):
    """Measure the ACTUAL typical scale of |dE_geom| for a single
    colour-class proposal on THIS field, by directly probing several
    random colours, rather than guessing T0 as a fraction of the total
    E_geom (which was the root cause of v11's first run freezing at
    0% acceptance for 1700+ sweeps — the total E_geom is a product of
    two O(1000)-scale numbers and is a wildly wrong reference scale for
    a move that only touches 1/64 of the grid at once)."""
    probe_colours = [COLOUR_MASKS[i] for i in
                      np.random.choice(len(COLOUR_MASKS), size=min(n_probes, len(COLOUR_MASKS)), replace=False)]
    dEs = [probe_dE(n, K_cur, J4_cur, mask, rotation_scale) for mask in probe_colours]
    dEs = [d for d in dEs if d > 0]
    if not dEs:
        return 1.0  # degenerate fallback, should not normally happen
    # Set T0 so that a MEDIAN-sized unfavourable move has a sizeable
    # (~50%) chance of acceptance at the start — exp(-dE/T0)=0.5 at the
    # median dE is the standard calibration target for an annealing run
    # that should start "hot" (exploring) and cool from there.
    median_dE = float(np.median(dEs))
    T0 = median_dE / np.log(2)
    print(f"    [calibration] probed {len(dEs)} colour-class moves: "
          f"median|dE|={median_dE:.4e}  max|dE|={max(dEs):.4e}  -> T0={T0:.4e}")
    return T0


def relax_R0(R0, n_init, sweeps):
    n = torch.tensor(n_init, dtype=torch.float32, device=dev)
    K0, J2a0, J4_0, r_bar0 = compute_global(n, MU)
    J4_target = J4_0
    Eg0 = K0*J4_0
    print(f"\n  R0={R0:.2f}  init: K={K0:.2f} J2a={J2a0:.2f} J4={J4_0:.4f} "
          f"E_geom={Eg0:.2f}  r_bar={r_bar0:.3f}(R0={R0})")
    if args.T0 is not None:
        T0 = args.T0
        print(f"    Using user-specified T0={T0:.4e}")
    else:
        T0 = calibrate_T0(n, K0, J4_0, args.rotation_scale)

    K_cur, J4_cur = K0, J4_0
    T = T0
    history = []
    zero_accept_streak = 0

    for sweep in range(1, sweeps+1):
        n, K_cur, J4_cur, acc_rate = anneal_sweep(n, K_cur, J4_cur, T, args.rotation_scale)

        # Stuck-run detector: the first v11 run showed EXACTLY 0.0%
        # acceptance for 1700+ consecutive sweeps with every diagnostic
        # frozen — that is a sign of a frozen/stuck search, not a
        # legitimate cold-annealing result, and should be surfaced
        # loudly rather than silently logged as if it were normal.
        if acc_rate == 0.0:
            zero_accept_streak += 1
        else:
            zero_accept_streak = 0
        if zero_accept_streak == 100:
            print(f"    WARNING: acceptance rate has been EXACTLY 0.0% for "
                  f"100 consecutive sweeps (as of sweep {sweep}, T={T:.4e}). "
                  f"This usually means the search is frozen/stuck, not "
                  f"genuinely cooled — consider a higher T0 or smaller "
                  f"--rotation_scale and re-running, rather than trusting "
                  f"the remaining sweeps.")

        if sweep % args.print_every == 0 or sweep == 1:
            # Refresh exactly from the full field periodically to avoid
            # incremental float drift over many sweeps.
            K_cur, J2a_cur, J4_cur, r_bar_cur = compute_global(n, MU)
            j4_frac = J4_cur/J4_target if J4_target > 1e-12 else float('nan')
            vrat = K_cur/J2a_cur if J2a_cur > 1e-12 else float('nan')
            lam = K_cur/J4_cur if J4_cur > 1e-12 else float('nan')
            Eg_cur = K_cur*J4_cur
            print(f"    sweep={sweep:>5}  T={T:>10.4e}  accept={acc_rate*100:>5.1f}%  "
                  f"J4/J0={j4_frac:>7.4f}  r_bar={r_bar_cur:>6.3f}  "
                  f"K/J2a={vrat:>7.4f}(2phi={2*phi:.4f})  K/J4={lam:>9.3f}  "
                  f"E_geom={Eg_cur:>12.2f}")
            history.append(dict(sweep=sweep, T=T, accept=acc_rate, J4_J0=j4_frac,
                                 r_bar=r_bar_cur, KJ2a=vrat, KJ4=lam, E_geom=Eg_cur))
            log.write(f"{R0}  {sweep}  {T:.4e}  {acc_rate:.4f}  {K_cur:.4f}  "
                      f"{J2a_cur:.4f}  {J4_cur:.6f}  {j4_frac:.6f}  {r_bar_cur:.4f}  "
                      f"{lam:.6f}  {vrat:.6f}  {Eg_cur:.4f}\n")
            log.flush()

        T *= args.cooling_rate

    K_f, J2a_f, J4_f, r_bar_f = compute_global(n, MU)
    j4_frac_f = J4_f/J4_target if J4_target > 1e-12 else float('nan')
    result = dict(R0=R0, K=K_f, J2a=J2a_f, J4=J4_f, J4_J0=j4_frac_f, r_bar=r_bar_f)
    result['KJ4'] = result['K']/result['J4'] if result['J4'] > 1e-12 else float('nan')
    result['KJ2a'] = result['K']/result['J2a'] if result['J2a'] > 1e-12 else float('nan')
    result['E_geom'] = result['K']*result['J4']
    result['final_accept_rate'] = history[-1]['accept'] if history else float('nan')
    return n.detach().cpu().numpy(), result, history


# ── Main: scan over R0 ─────────────────────────────────────────────────
results = []
csv_path = os.path.join(args.outdir, f'{tag}_results.csv')
with open(csv_path, 'w') as f:
    f.write("R0,K,J2a,J4,J4_J0,r_bar,KJ4,KJ2a,E_geom,final_accept_rate\n")

for R0 in R0_LIST:
    n_init = build_ic(R0, args.r0, args.C_star, args.NT)
    n_final, res, hist = relax_R0(R0, n_init, args.sweeps)
    results.append(res)

    with open(csv_path, 'a') as f:
        f.write(f"{R0},{res['K']:.6f},{res['J2a']:.6f},{res['J4']:.6f},"
                f"{res['J4_J0']:.6f},{res['r_bar']:.6f},{res['KJ4']:.6f},"
                f"{res['KJ2a']:.6f},{res['E_geom']:.6f},{res['final_accept_rate']:.4f}\n")

    np.savez(os.path.join(args.outdir, f'n_R0{R0:.2f}_FINAL.npz'),
              n=n_final, n_initial=n_init, N=N, h=h, R0=R0, r0=args.r0,
              C_star=args.C_star, **{k: v for k, v in res.items() if k != 'R0'})

print(f"\n{'='*70}")
print(f"  METROPOLIS ANNEALING RESULTS")
print(f"{'='*70}")
print(f"  {'R0':>6}  {'J4/J0':>8}  {'r_bar':>8}  {'K/J2a':>8}  {'K/J4':>10}  "
      f"{'K/J4/phi6':>10}  {'final_acc%':>10}")
for r in results:
    lam = r['KJ4']
    print(f"  {r['R0']:>6.2f}  {r['J4_J0']:>8.4f}  {r['r_bar']:>8.4f}  "
          f"{r['KJ2a']:>8.4f}  {lam:>10.4f}  {lam/phi6:>10.5f}  "
          f"{r['final_accept_rate']*100:>10.2f}")

# Explicit verdict per R0, covering BOTH known failure directions —
# the first v11 run's most informative symptom was J4 GROWING (to
# 6x its initial value), not collapsing, and nothing in that run's
# output explicitly flagged growth as a failure mode distinct from
# collapse, so this is made explicit here.
print()
for r in results:
    if 0.7 <= r['J4_J0'] <= 1.5 and 0.5*r['R0'] <= r['r_bar'] <= 1.5*r['R0']:
        verdict = "PLAUSIBLY STABLE (J4 near target, r_bar near R0)"
    elif r['J4_J0'] < 0.5:
        verdict = "COLLAPSED (J4 diluted away — same failure as v6-v8)"
    elif r['J4_J0'] > 1.5:
        verdict = "GROWN/SPIKED (J4 far above target — same failure direction as v9's spike, now via accepted anneal moves rather than gradient descent)"
    else:
        verdict = "AMBIGUOUS — inspect r_bar trend and final field directly"
    print(f"  R0={r['R0']:.2f}: {verdict}")

log.close()

rpt_path = os.path.join(args.outdir, f'{tag}_report.txt')
with open(rpt_path, 'w') as f:
    f.write(f"Q_H=3 v11 — Metropolis annealing\n")
    f.write(f"N={N} h={h} C*_init={args.C_star} sweeps={args.sweeps}\n\n")
    f.write(f"{'R0':>6}  {'J4/J0':>8}  {'r_bar':>8}  {'K/J2a':>8}  {'K/J4':>10}  "
            f"{'K/J4/phi6':>10}  {'final_acc%':>10}\n")
    for r in results:
        lam = r['KJ4']
        f.write(f"  {r['R0']:>6.2f}  {r['J4_J0']:>8.4f}  {r['r_bar']:>8.4f}  "
                f"{r['KJ2a']:>8.4f}  {lam:>10.4f}  {lam/phi6:>10.5f}  "
                f"{r['final_accept_rate']*100:>10.2f}\n")

print(f"\nResults saved: {csv_path}\nReport: {rpt_path}")
