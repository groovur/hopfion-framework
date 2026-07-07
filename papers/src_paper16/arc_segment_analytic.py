#!/usr/bin/env python3
"""
arc_segment_analytic.py
=======================
Compute rho_J4 arc-segment energy integrals on the ANALYTIC Construction C
field (exactly Z3-symmetric by construction).

6 segments in correct t order (monotone, no wraparound bug):
  Seg 0  B0: C0 -> M0   t in [pi/6,    pi/3]     30 deg  (crossing->midpoint)
  Seg 1  A0: M0 -> C1   t in [pi/3,    pi/6+2p/3] 90 deg  (midpoint->crossing)
  Seg 2  B1: C1 -> M1   t in [pi/6+2p/3, pi/3+2p/3] 30 deg
  Seg 3  A1: M1 -> C2   t in [pi/3+2p/3, pi/6+4p/3] 90 deg
  Seg 4  B2: C2 -> M2   t in [pi/6+4p/3, pi/3+4p/3] 30 deg
  Seg 5  A2: M2 -> C0   t in [pi/3+4p/3, pi/6+2pi]  90 deg (wraps)

B-segments (C->M, 30 deg): should be Z3-equal, SHORT, high energy density
A-segments (M->C, 90 deg): should be Z3-equal, LONG,  lower energy density
"""
import numpy as np, sys, time
sys.path.insert(0,'.')
from bishop_frame_v2 import build_compensated_frame_arclength
from scipy.spatial import KDTree

R0, r0, C_star = 3.0, 0.874, 2.5062
PHI = (1+5**0.5)/2
MU  = 3.0 - PHI

# Correct monotone boundary sequence
C0 = np.pi/6;         M0 = np.pi/3
C1 = np.pi/6+2*np.pi/3; M1 = np.pi/3+2*np.pi/3
C2 = np.pi/6+4*np.pi/3; M2 = np.pi/3+4*np.pi/3

BOUNDS  = np.array([C0, M0, C1, M1, C2, M2])          # monotone
SEG_HI  = np.array([M0, C1, M1, C2, M2, C0+2*np.pi])  # upper bound (last wraps)
SEG_NAMES = ['B0 (C0->M0, 30°)', 'A0 (M0->C1, 90°)',
             'B1 (C1->M1, 30°)', 'A1 (M1->C2, 90°)',
             'B2 (C2->M2, 30°)', 'A2 (M2->C0, 90°)']
SEG_TYPE  = ['B','A','B','A','B','A']

def assign_seg(t_arr):
    """Vectorised: assign each t in [0,2pi) to segment 0-5."""
    t = t_arr % (2*np.pi)
    seg = np.zeros(len(t), dtype=int)
    for s in range(5):
        seg[(t >= BOUNDS[s]) & (t < BOUNDS[s+1])] = s
    # segment 5: M2 -> C0 (wraps: t>=M2 OR t<C0)
    seg[(t >= M2) | (t < C0)] = 5
    return seg

# Build curve samples and KDTree
NT = 6000
t_arr = np.linspace(0, 2*np.pi, NT, endpoint=False)
Gx=(R0+r0*np.cos(3*t_arr))*np.cos(2*t_arr)
Gy=(R0+r0*np.cos(3*t_arr))*np.sin(2*t_arr)
Gz=r0*np.sin(3*t_arr)
curve_pts = np.stack([Gx,Gy,Gz],axis=1)
main_tree = KDTree(curve_pts)
seg_of_sample = assign_seg(t_arr)

# Verify spans
print("Segment spans (t-width):")
for s,(nm,lo,hi) in enumerate(zip(SEG_NAMES,BOUNDS,SEG_HI)):
    span = hi-lo
    print(f"  Seg {s} {nm}: {np.degrees(span):.1f} deg")
print()

# Lobe trees for Construction C
arc_starts = [0, 2*np.pi/3, 4*np.pi/3]
lobe_indices  = [np.where((t_arr>=s)&(t_arr<s+2*np.pi/3))[0] for s in arc_starts]
lobe_trees    = [KDTree(curve_pts[li]) for li in lobe_indices]
lobe_t_arrays = [t_arr[li] for li in lobe_indices]

NT_frame=20000
t_frame,_,N1f,N2f,_=build_compensated_frame_arclength(NT=NT_frame)

def nearest_two(pts):
    d_per,t_per=[],[]
    for tr,tl in zip(lobe_trees,lobe_t_arrays):
        d,idx=tr.query(pts,workers=1); d_per.append(d); t_per.append(tl[idx])
    d_s=np.stack(d_per,axis=1); t_s=np.stack(t_per,axis=1)
    o=np.argsort(d_s,axis=1)
    return (np.take_along_axis(t_s,o,axis=1)[:,0],
            np.take_along_axis(d_s,o,axis=1)[:,0],
            np.take_along_axis(t_s,o,axis=1)[:,1],
            np.take_along_axis(d_s,o,axis=1)[:,1])

def frame_at(tq):
    idx=np.searchsorted(t_frame,tq%(2*np.pi))%NT_frame
    return N1f[idx],N2f[idx]
def curve_at(t):
    return np.stack([(R0+r0*np.cos(3*t))*np.cos(2*t),
                     (R0+r0*np.cos(3*t))*np.sin(2*t),
                      r0*np.sin(3*t)],axis=-1)

def build_n(pts):
    t1,d1,t2,d2=nearest_two(pts)
    chi1=np.arctan2(np.sum((pts-curve_at(t1))*frame_at(t1)[1],axis=1),
                    np.sum((pts-curve_at(t1))*frame_at(t1)[0],axis=1))
    chi2=np.arctan2(np.sum((pts-curve_at(t2))*frame_at(t2)[1],axis=1),
                    np.sum((pts-curve_at(t2))*frame_at(t2)[0],axis=1))
    Phi1=chi1+3*t1; Phi2=chi2+3*t2
    rho1=np.clip(d1,1e-6,None); rho2=np.clip(d2,1e-6,None)
    f1=2*np.arctan(np.maximum(rho1*C_star,1e-9)**(-C_star))
    f2=2*np.arctan(np.maximum(rho2*C_star,1e-9)**(-C_star))
    w1=1/rho1**2; w2=1/rho2**2
    z1=(w1*np.cos(f1/2)+w2*np.cos(f2/2)).astype(complex)
    z2=w1*np.sin(f1/2)*np.exp(1j*Phi1)+w2*np.sin(f2/2)*np.exp(1j*Phi2)
    mag=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=mag; z2/=mag
    nx=2*np.real(np.conj(z1)*z2); ny=2*np.imag(np.conj(z1)*z2)
    nz=np.abs(z1)**2-np.abs(z2)**2
    n=np.stack([nx,ny,nz],axis=-1)
    return (n/np.linalg.norm(n,axis=-1,keepdims=True).clip(1e-10)).astype(np.float32)

def run(N, h):
    t0=time.time()
    cv=h*(np.arange(N)-N//2+0.5)
    X,Y,Z=np.meshgrid(cv,cv,cv,indexing='ij')
    pts=np.stack([X.ravel(),Y.ravel(),Z.ravel()],axis=1).astype(np.float32)
    n=build_n(pts).reshape(N,N,N,3)

    nx,ny,nz=n[...,0],n[...,1],n[...,2]
    def cd(u,a): return (np.roll(u,-1,a)-np.roll(u,1,a))/(2*h)
    nxx,nxy,nxz=cd(nx,0),cd(nx,1),cd(nx,2)
    nyx,nyy,nyz=cd(ny,0),cd(ny,1),cd(ny,2)
    nzx,nzy,nzz=cd(nz,0),cd(nz,1),cd(nz,2)
    Fxy=nx*(nyx*nzy-nzx*nyy)+ny*(nzx*nxy-nxx*nzy)+nz*(nxx*nyy-nyx*nxy)
    Fxz=nx*(nyx*nzz-nzx*nyz)+ny*(nzx*nxz-nxx*nzz)+nz*(nxx*nyz-nyx*nxz)
    Fyz=nx*(nyy*nzz-nzy*nyz)+ny*(nzy*nxz-nxy*nzz)+nz*(nxy*nyz-nyy*nxz)
    rho=(Fxy**2+Fxz**2+Fyz**2)
    J4_total=float(rho.sum()*h**3)

    _,nn_idx=main_tree.query(pts.reshape(-1,3))
    seg_of_pt=seg_of_sample[nn_idx].reshape(N,N,N)

    seg_J4=np.zeros(6)
    pt_count=np.zeros(6,dtype=int)
    for s in range(6):
        mask=(seg_of_pt==s)
        seg_J4[s]=float(rho[mask].sum()*h**3)
        pt_count[s]=int(mask.sum())

    elapsed=time.time()-t0
    return J4_total, seg_J4, pt_count, elapsed

print("="*65)
print("ARC-SEGMENT ENERGY — ANALYTIC CONSTRUCTION C")
print("="*65)

results = {}
for N,h in [(64,0.175),(96,0.113)]:
    print(f"\n--- N={N}, h={h:.3f} ({time.strftime('%H:%M:%S')}) ---")
    J4, seg_J4, pt_count, elapsed = run(N, h)
    results[(N,h)] = (J4, seg_J4)

    B_segs = seg_J4[[0,2,4]]   # C->M, 30 deg each
    A_segs = seg_J4[[1,3,5]]   # M->C, 90 deg each

    # Energy DENSITY per degree (normalise by arc-length in t)
    B_density = B_segs / (np.pi/6)   # per radian of t
    A_density = A_segs / (np.pi/2)

    print(f"  J4_total = {J4:.3f}  ({elapsed:.1f}s)")
    print(f"\n  {'Segment':<24} {'J4':>10} {'%tot':>7} {'pts':>8} {'J4/rad':>10}")
    for i,(nm,v,pc) in enumerate(zip(SEG_NAMES,seg_J4,pt_count)):
        rad = np.pi/6 if SEG_TYPE[i]=='B' else np.pi/2
        print(f"  {nm:<24} {v:>10.3f} {100*v/J4:>7.2f} {pc:>8d} {v/rad:>10.3f}")

    print(f"\n  B-segments (C->M, 30°, n=3): {B_segs}")
    print(f"    mean={B_segs.mean():.3f}  std/mean={B_segs.std()/B_segs.mean():.4f}")
    print(f"  A-segments (M->C, 90°, n=3): {A_segs}")
    print(f"    mean={A_segs.mean():.3f}  std/mean={A_segs.std()/A_segs.mean():.4f}")

    print(f"\n  Energy density ratio B/A (per radian): {B_density.mean()/A_density.mean():.4f}")
    print(f"  Raw integral  ratio B/A (per segment): {B_segs.mean()/A_segs.mean():.4f}")
    print(f"  3*(B-A) = {3*(B_segs.mean()-A_segs.mean()):.3f}")
    print(f"  J4/phi^6 = {J4/((1+5**0.5)/2)**6:.4f}")

import time
