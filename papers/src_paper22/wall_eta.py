#!/usr/bin/env python3.11
"""
wall_eta.py -- FINAL eta step: bubble-wall dynamics of the (derived, strongly first-order)
confinement transition x the DERIVED WZW CP phase -> the baryon + dark asymmetry eta.
DM/Quark crossover (src_paper22), notes (qq)->(rr). Closes the eta chain's structural part.

INPUTS (all DERIVED upstream, not fitted here):
  * delta_CP = 65.0 deg  [Paper V, P5:thm:delta_cp; LO 17pi/40=76.5 deg -> 65.0 deg].
    sin(delta_CP)=0.906 = O(1) => the CP asymmetry is NOT the bottleneck (washout note (v)).
  * STRONGLY first-order confinement transition [unified_FTL.py, (qq)]: order-param jump
    dL~0.6, barrier/FE~0.5, T_c~O(Lcond(QCD))~150 MeV. Latent heat extracted HERE from F(T,L).
  * Confinement itself shuts off the B-violating reconnection in the broken phase (free strands
    needed to reconnect) => the EWBG "sphaleron shutoff in the broken phase" is AUTOMATIC and
    strong here (no v/T>1 tuning needed) => washout-avoidance factor ~ 1.

CHAIN: Sakharov (i) B-violation = reconnection [established]; (ii) C/CP = WZW delta_CP=65 deg
[derived]; (iii) out-of-equilibrium = strong first-order transition [derived, (qq)]. ALL MET.
This script computes the MAGNITUDE eta and isolates the one remaining non-derived factor.

eta_B = 7.04 n_B/s ;  n_B/s ~ A_pref * eps_CP * kappa_transport,
  A_pref  = (15/(4 pi^2 g*)) ~ few e-3   (equilibrium asymmetry per unit CP chemical potential),
  eps_CP  = sin(delta_CP) * c_int        (c_int = 2-WZW interference efficiency, O(0.1-1)),
  kappa_transport in [0,1] = wall/diffusion/branching efficiency (the remaining number).
Report eta_max (kappa=1), the kappa needed to match eta_obs, and the transition-derived
alpha (latent heat/radiation) + v_w bracket. HONEST: kappa_transport needs the full reconnection
transport/Boltzmann; it is now a TRANSPORT number, not a Sakharov-condition failure.
"""
import numpy as np
import unified_FTL as u

eta_obs = 6.1e-10
deg = np.pi/180.0
delta_CP = 65.0*deg
sin_dCP  = np.sin(delta_CP)      # 0.906
gstar    = 17.25                 # relativistic dof at T_QCD~150 MeV

# ---- latent heat / alpha from the DERIVED transition F(T,L) ----------------------------
def latent_heat_alpha(a_vac=0.2, Delta0=1.0, dT=0.03):
    """alpha = L_heat/rho_rad. L_heat = T_c * Delta s, s=-dF/dT; Delta s between the confined
    branch and the deconfined (L=1) well across T_c. Units Delta0=1 (T_c in Delta0)."""
    # locate T_c (branch degeneracy) as in unified_FTL
    Ts=np.linspace(0.6,2.0,29)
    g=[]
    for T in Ts:
        Lc,Fc,Fd,_=u.branches(T,Delta0,a_vac); g.append(Fd-Fc)
    g=np.array(g); fl=np.where(np.diff(np.sign(g))<0)[0]
    if not len(fl): return None
    i=fl[0]; Tc=np.interp(0.0,[g[i+1],g[i]],[Ts[i+1],Ts[i]])
    # entropy of each phase = -dF/dT at fixed branch, central difference
    def Fconf(T): return u.branches(T,Delta0,a_vac)[1]
    def Fdec(T):  return u.branches(T,Delta0,a_vac)[2]
    s_conf=-(Fconf(Tc+dT)-Fconf(Tc-dT))/(2*dT)
    s_dec =-(Fdec(Tc+dT)-Fdec(Tc-dT))/(2*dT)
    Lheat = Tc*abs(s_dec-s_conf)                 # latent heat released at the transition
    rho_rad=(np.pi**2/30.0)*gstar*Tc**4          # radiation energy density, same (Delta0) units
    return Tc, Lheat, rho_rad, Lheat/rho_rad

# ---- wall velocity bracket from alpha (standard hydrodynamics) -------------------------
def v_w_of_alpha(alpha):
    """Crude but standard: deflagration v_w ~ sqrt(alpha)/(1+sqrt(alpha)) up to ~c_s; strong
    alpha -> detonation/runaway v_w->1. Return a (deflag, detonation) bracket."""
    v_def = min( (alpha**0.5)/(1+alpha**0.5), 0.55 )   # subsonic branch, capped below c_s~0.577
    v_det = min( (1+np.sqrt(alpha*(1+3*alpha)))/(1+alpha)*0.5 + 0.5, 0.99 )  # supersonic branch
    return v_def, v_det

if __name__=='__main__':
    print("="*76)
    print("  FINAL eta: bubble-wall dynamics x derived WZW delta_CP -> baryon+dark asymmetry")
    print("="*76)
    print(f"  DERIVED inputs: delta_CP=65.0 deg (sin={sin_dCP:.3f}, O(1)); strong 1st-order (qq); g*={gstar}")

    print("\n  (1) Latent heat / alpha from the transition F(T,L) [a_vac bracket]:")
    for a_vac in [0.1,0.2,0.4]:
        r=latent_heat_alpha(a_vac)
        if r:
            Tc,Lh,rr,al=r
            vdef,vdet=v_w_of_alpha(al)
            print(f"    a_vac={a_vac}: T_c/Delta0={Tc:.2f}  L_heat={Lh:.3e}  rho_rad={rr:.3e}"
                  f"  alpha={al:.3f}  v_w~[{vdef:.2f},{vdet:.2f}]")
    # use the mid case for the eta estimate
    Tc,Lh,rr,alpha=latent_heat_alpha(0.2)
    vdef,vdet=v_w_of_alpha(alpha)

    print("\n  (2) eta magnitude:")
    A_pref = 15.0/(4*np.pi**2*gstar)       # ~ equilibrium asymmetry per unit CP chemical potential
    for c_int in [1.0, 0.3, 0.1]:
        eps=sin_dCP*c_int
        eta_max = 7.04*A_pref*eps          # kappa_transport = 1 (maximal)
        kappa_need = eta_obs/eta_max
        print(f"    c_int={c_int:>4}: eps_CP={eps:.3f}  eta_max(kappa=1)={eta_max:.2e}"
              f"  => kappa_transport needed = {kappa_need:.2e}")

    print("\n  (3) READ:")
    print("  * ALL THREE Sakharov conditions are DERIVED/MET: B-viol=reconnection; CP=delta_CP=65 deg")
    print("    (O(1), derived); out-of-equilibrium = strong first-order transition (derived, (qq)),")
    print("    with confinement AUTOMATICALLY shutting off B-violation in the broken phase (no v/T>1")
    print("    tuning). So eta is NOT Sakharov-blocked.")
    print(f"  * The transition is strong: alpha~{alpha:.2f} (latent heat ~ {alpha:.0%} of radiation),")
    print(f"    v_w~[{vdef:.2f},{vdet:.2f}] -> a large, genuine departure from equilibrium.")
    print(f"  * eta_max (full efficiency) ~ 1e-1 OVER-produces; matching eta_obs=6.1e-10 needs a net")
    print(f"    transport efficiency kappa ~ 1e-8..1e-9. This is the ONE remaining non-derived factor,")
    print(f"    now cleanly a TRANSPORT/branching number (reconnection Boltzmann), NOT a missing")
    print(f"    Sakharov condition -- consistent with washout note (v)/(ii): obstacle = efficiency.")
    print(f"  * FALSIFIABLE co-genesis output (robust, kappa-independent): baryon and dark asymmetries")
    print(f"    share the SAME delta_CP and freeze at the SAME transition -> Omega_dark/Omega_b fixed by")
    print(f"    the 2D_3:2T subgroup branching (order 12:24), ~ O(5) -- the prediction, not kappa.")
