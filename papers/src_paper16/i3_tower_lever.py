#!/usr/bin/env python3.11
"""
I3 -> tower-coupling lever for the up/down quark mass split.
Session 2026-09-21. Follows nq_base_tower.py (n_q^base=6 derived) and
notes/quark_generation_e8_ribbon_twist.md sec.20-21.

QUESTION: the up-sector tower climbs ~2x steeper than down (sec.21). Where does
isospin enter the mass formula  m = Lcond phi^{-2n} exp(-2pi Q_group T_g)?
Test candidate structures HONESTLY (report failures; flag single-data-point leads).

Candidates:
  H_Q  : steepness (per-gen increment) scales with |electric charge| (up 2/3, down 1/3 -> ratio 2)
  H_T  : split proportional to the T-matrix phase T_g (lepton phases 5/24,1/8,~0.075)
  H_tw : split proportional to the E_6 twist isospin-split {2,4,10}
Also: demonstrate the T-vs-n degeneracy that blocks a unique assignment.
"""
import numpy as np
PHI=(1+5**0.5)/2; LN=np.log(PHI); PI=np.pi
m_q0=215.5  # MeV constituent baryon scale (tower anchor)

# masses (MeV): current/MSbar. top: measured pole AND E_8 full-formation (334 GeV, sec.18)
m   = {'u':2.16,'d':4.67,'s':93.4,'c':1270.,'b':4180.,'t':172760.}
m_t_full = 334000.  # E_8 predicts full T(2,3) formation; measured 173 GeV = incomplete (bare top)
gen = {'u':1,'d':1,'s':2,'c':2,'b':3,'t':3}
Qem = {'u':2/3,'d':-1/3,'s':-1/3,'c':2/3,'b':-1/3,'t':2/3}   # electric charge
Tw  = {'u':+1,'d':-1,'c':+2,'s':-2,'t':+5,'b':-5}            # E_6 per-strand twist (note sec.13)
Tlep= {1:5/24,2:1/8,3:0.075}                                 # lepton/quark T-matrix phase per gen

def nreq(mass): return -np.log(mass/m_q0)/(2*LN)
n={q:nreq(m[q]) for q in m}
n_t_full=nreq(m_t_full)

print("="*72,"\nLUMPED-n tower levels  (m = m_q0 phi^{-2n})\n","="*72,sep="")
for q in ['u','d','s','c','b','t']:
    print(f"  {q}: n={n[q]:+6.3f}  gen{gen[q]}  Qem={Qem[q]:+.3f}")
print(f"  (top full-formation 334 GeV -> n={n_t_full:+.3f})")

print("\n"+"="*72,"\nPER-GENERATION INCREMENTS and slope ratio up/down\n","="*72,sep="")
up=['u','c','t']; dn=['d','s','b']
dU=[n[up[1]]-n[up[0]], n[up[2]]-n[up[1]]]           # measured top
dU_full=[n[up[1]]-n[up[0]], n_t_full-n[up[1]]]      # full-formation top
dD=[n[dn[1]]-n[dn[0]], n[dn[2]]-n[dn[1]]]
print(f"  DOWN incr: d->s={dD[0]:+.3f}  s->b={dD[1]:+.3f}")
print(f"  UP   incr: u->c={dU[0]:+.3f}  c->t={dU[1]:+.3f}   (measured top)")
print(f"  UP   incr: u->c={dU_full[0]:+.3f}  c->t={dU_full[1]:+.3f}   (E_8 full-formation top)")
print(f"  slope ratio up/down  gen1->2 = {dU[0]/dD[0]:.3f}   [ONLY clean comparison; top-free]")
print(f"  slope ratio up/down  gen2->3 = {dU[1]/dD[1]:.3f} (meas top, contaminated) / "
      f"{dU_full[1]/dD[1]:.3f} (full-form top)")

print("\n"+"="*72,"\nH_Q: steepness ~ |electric charge|  (ratio |Q_up|/|Q_down| = 2)\n","="*72,sep="")
print(f"  predicted slope ratio = |2/3|/|1/3| = 2.000")
print(f"  measured  slope ratio gen1->2 = {dU[0]/dD[0]:.3f}   (off {100*(dU[0]/dD[0]/2-1):+.0f}%)")
# predict m_c/m_u by stepping up-quark with 2x the down increment
mc_pred = m['u']*PHI**(-2*(2*dD[0]))
print(f"  predict m_c: step u by 2x(down incr)=2*({dD[0]:+.2f}) -> m_c={mc_pred:.0f} MeV "
      f"(meas 1270);  m_c/m_u pred {mc_pred/m['u']:.0f} vs meas {m['c']/m['u']:.0f}")
print("  => slope-level agreement decent (2.13 vs 2) but exp amplifies 6% -> mass off ~50%.")
print("     ONE clean data point (gen1->2); gen2->3 top-contaminated. NOT a confirmed law.")

print("\n"+"="*72,"\nH_T: split ~ T_g   /   H_tw: split ~ E_6 twist-split {2,4,10}\n","="*72,sep="")
split={g:n[{1:'u',2:'c',3:'t'}[g]]-n[{1:'d',2:'s',3:'b'}[g]] for g in (1,2,3)}
tw_split={1:2,2:4,3:10}
print("  gen   n-split   T_g     twist-split   split/T_g   split/tw")
for g in (1,2,3):
    print(f"   {g}   {split[g]:+6.3f}  {Tlep[g]:.4f}     {tw_split[g]:2d}        "
          f"{split[g]/Tlep[g]:+7.2f}    {split[g]/tw_split[g]:+.3f}")
print("  H_T: split/T_g NOT constant AND T_g DECREASES while |split| GROWS -> WRONG TREND.")
print("  H_tw: split/tw NOT constant AND twist-split all + while n-split flips sign -> FAILS.")

print("\n"+"="*72,"\nT-vs-n DEGENERACY (why the assignment is not unique)\n","="*72,sep="")
print("  m = Lcond phi^{-2n} exp(-2pi*3*T): a single mass fixes only the combination")
print("  2n*lnphi + 6pi*T. Fixing T from the E_6 twist (T=Tw/3, SU(3)_1 h=1/3) forces n:")
for q in ['u','d']:
    T=Tw[q]/3.0
    # ln m = ln Lcond -2n lnphi -6pi T  (absorb Lcond via m_q0 anchor at n=0,T=0)
    n_forced = -(np.log(m[q]/m_q0)+6*PI*T)/(2*LN)
    print(f"    {q}: T=Tw/3={T:+.3f} -> n_forced={n_forced:+.2f}  (vs lumped n={n[q]:+.2f})")
print("  => the same masses are reproduced with wildly different n once T carries the twist.")
print("     Without an INDEPENDENT fix of n OR T per quark, isospin cannot be uniquely")
print("     assigned to n vs T. This degeneracy is the core obstruction (sec.21 [O]).")

print("\n"+"="*72,"\nE_8-ROUTE CROSS-CHECK (unambiguous: no n-T degeneracy; sec.17)\n","="*72,sep="")
# m_q/m_ell(g) = exp(n_q pi/9), n_q = Base_g - 2m. Down works, up fails (sec.17).
mlep={1:0.511,2:105.658,3:1776.86}  # e,mu,tau
Base={1:41,2:26,3:17}
Ecox={'d':17,'s':13,'b':7,'u':19,'c':11,'t':1}   # assigned E_8 exponents (sec.17 remark)
mtop={'t':m_t_full}                              # top full-formation for the check
print("   q  g  E8_m  n_q(E8)  n_q(req)  deficit=req-E8   pred/meas")
for q in ['d','s','b','u','c','t']:
    g=gen[q]; nq_e8=Base[g]-2*Ecox[q]
    mq = m_t_full if q=='t' else m[q]
    nq_req = 9*np.log(mq/mlep[g])/PI
    pred = mlep[g]*np.exp(nq_e8*PI/9)
    print(f"   {q}  {g}   {Ecox[q]:2d}   {nq_e8:+5.1f}   {nq_req:+6.2f}    {nq_req-nq_e8:+6.2f}"
          f"          {pred/mq:.2f}")
# extra up correction beyond the ~ -0.5 uniform down (running) deficit
dd=[9*np.log(m[q]/mlep[gen[q]])/PI-(Base[gen[q]]-2*Ecox[q]) for q in ['d','s','b']]
du=[9*np.log((m_t_full if q=='t' else m[q])/mlep[gen[q]])/PI-(Base[gen[q]]-2*Ecox[q]) for q in ['u','c','t']]
print(f"  down deficits {[f'{x:+.2f}' for x in dd]} ~ uniform -0.5 (=ln(1.2)*9/pi, the running factor)")
print(f"  up   deficits {[f'{x:+.2f}' for x in du]}  (top~0: full-formation OK, sec.18)")
print(f"  confined up extra beyond uniform: u {du[0]-np.mean(dd):+.2f}, c {du[1]-np.mean(dd):+.2f}"
      f"  -> ratio {(du[1]-np.mean(dd))/(du[0]-np.mean(dd)):.2f}")
print("  => SAME picture: confined up needs an extra gen-growing suppression, ratio ~2-2.8/gen.")

print("\n"+"="*72,"\nVERDICT\n","="*72,sep="")
print("  - steepness ratio up/down = 2.13 (one clean comparison), NEAR electric-charge ratio 2")
print("    -> H_Q is the most STRUCTURED lead (charge is a real, mass-independent framework")
print("       quantum number; steepness~charge makes up-down splitting GROW with gen naturally),")
print("       but underdetermined: 1 data point + exp sensitivity. NOT claimed.")
print("  - H_T (T_g-scaling) and H_tw (twist-scaling) both FAIL (wrong trend / wrong sign).")
print("  - core obstruction = T-vs-n degeneracy: need an independent per-quark handle on n or T.")
print("  INDEPENDENT CHECK to settle H_Q: a 2nd top-free up/down slope comparison -- i.e. the")
print("  confined c->? increment without the bare top, OR derive steepness~|Q| from the writhe/")
print("  self-linking origin of electric charge (alpha=360/phi^2) in the formation energy.")
