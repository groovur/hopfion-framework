# WP2 — physical dispersion scales (v, A, anisotropy) from S_eff

GOAL: pin the physical semi-Dirac parameters (v, A) and the condensate length so WP3
(cutoff nature/value) and WP4 (geom_coeff) can run. STATUS: partial success + one catch.

## Framework ingredients (Paper 7, exact)
Parent action (P7:eq:parent_action):
  S = INT d4x sqrt(-g) [ g^{mu nu} d_mu rho d_nu rho / (phi^6 (1+beta rho)) + Lam0 e^{-phi^6 rho} ]
      + INT S_eff dOmega,   S_eff = sin^4(theta) / [phi^6 (1+beta rho)].
Attractor rho_inf = phi/beta, so 1+beta rho_inf = 1+phi = phi^2.
Semi-Dirac (P7:eq:semidirac): H = A k_perp^2 sigma_x + v k_par sigma_y,
  E^2 = A^2 k^4 sin^4(theta) + v^2 k^2 cos^2(theta); heavy(perp)=sin^4=S_eff angular factor.

## Result A -- where phi^6 actually sits (SOLID, and it corrects the a0 framing)
The band-INTERNAL anisotropy A/v is O(1): band units A=1/2, v=1 => A/v=1/2, crossover
k*=v/A=2, E*=sqrt8. phi^6 is NOT in A/v. It is the S_eff SUPPRESSION MAGNITUDE
(sin^4/phi^6), i.e. it lives in the PHYSICAL NORMALISATION / unit map, not the band shape.
=> "m_xi/E* ~ phi^6" is a claim about the UNIT MAP (band E* -> physical energy), not about
the dispersion's internal structure. This is why band-shape studies (circular vs energy
cutoff) could never pin phi^6: phi^6 was never in the band shape. It is in the map. WP4.

## Result B -- the light-axis velocity v: c vs c/phi (CATCH, factor-phi, unresolved)
cutoff_embedding ASSUMED v -> c ("light-axis velocity v -> c", flagged MODELLED). But the
radial (light-axis) propagation of condensate fluctuations has a DERIVED speed: the
k-essence sound speed of the scalar/magnitude mode, c_s = 1/phi (P7:rem:cs_de, from
c_s^2 = 1 + rho K'/K = 1/phi^2 at the attractor). If the light-axis mode is that magnitude
mode, v = c/phi, NOT c -- a factor phi below the assumed value, which shifts A=c ell*,
E*=c/ell*, and the entire unit map by phi. Honest status: v is the DIRECTOR band's light-axis
velocity; c (relativistic kinetic scale) and c_s=c/phi (derived magnitude sound speed) are
both candidates and the framework does not cleanly identify which. FLAG for WP4: a factor phi
here directly changes the a0 coefficient and the m_xi/E* ratio.

## Result C -- absolute scale is set by R_0 = R_0(T_CMB) (=> WP4 needs it)
A = c ell*, E* = c/ell*: the crossover LENGTH ell* fixes A and E*. ell* is a condensate
length. The Linking-Scale results hold "for all R_0" (P7 l.320; Papers II/III) -- the
dimensionless ratios are R_0-independent, but the ABSOLUTE scale ell* is tied to the
condensate size R_0, which is fixed by the single input T_CMB (via Lam_cond =
T_CMB(pi^2/150)^{1/4} = 1.19e-4 eV). So the absolute crossover/healing scale is computable
in principle but needs R_0(T_CMB) and the profile -- that is WP4's unit map, not a pure number.

## Result D -- the one derived transverse scale
S_eff's second-order expansion gives the graviton (transverse TT) mass
  m_g^2 = 16 pi / (315 phi^6 (1+beta rho_inf)) = 16 pi / (315 phi^8)   (P7:eq:graviton_mass),
carrying the 1/phi^6 suppression x 1/phi^2 from 1+beta rho_inf. Physically m_g < 1e-29 eV
(P7 l.519). This is the TT sector; the director's transverse stiffness shares the 1/phi^6
normalisation structure but its absolute scale still rides on R_0 (Result C).

## WP2 verdict + handoff
FIXED: (i) phi^6 is in the unit map, not the band shape (Result A) -- so pinning must be a
unit-map calc, confirming WP4 is the crux. (ii) The transverse suppression normalisation is
1/phi^6, verified via m_g (Result D).
OPEN, handed to WP3/WP4: (i) v = c or c/phi (Result B) -- a factor-phi ambiguity that moves
the a0 coefficient; must be resolved by identifying the light-axis mode. (ii) the absolute
ell*/E*/geom_coeff needs R_0(T_CMB) + profile (Result C).
RISK UPDATE (R1 raised): with phi^6 in the map and the absolute scale tied to R_0 plus an
unresolved factor-phi in v, geom_coeff is a structured combination (R_0, phi-powers), not a
clean pure number -- so a unique phi^6 pinning is now LESS likely, and "a0 ~ cH_0 robust,
coefficient a specific-but-not-pure-number map" is the probable honest closure. WP4 decides.
