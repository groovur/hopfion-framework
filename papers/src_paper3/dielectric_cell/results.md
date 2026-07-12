# Dielectric-cell dipole polarizability: results
# Can you delegate building the actual coherent two-tube backgound solver to se where we land? 
# Can you update paper3.tex with these findings? (including the 

```
x_star  = 0.73294642
beta_star = 0.06715131
rho_c(beta_star) = 0.51913431
rho_c(beta_star/2) = None
rho_c(2*beta_star) = 0.89182102

=== V1: uniform medium (beta=0) ===
alpha = (-0+0j)  (expected 0)

=== V2: step inclusion eps1=3, R=1 ===
alpha = (0.49999999705789944-0j)  expected = 0.5  diff = 2.942e-09

=== V3: step inclusion eps1=-3+i*1e-3, R=1 ===
alpha = (1.9999997488784345+0.0004999998747190484j)  expected = (1.9999997500000621+0.0004999998750000312j)  diff = 1.122e-09

=== R1 beta=beta_star (beta=0.06715131, rho_c=0.5191343115235857) ===
  eta=1e-02: Re(alpha)=-0.96842863  Im(alpha)=2.66062182e-02  nfev=3266
  eta=1e-03: Re(alpha)=-0.97081256  Im(alpha)=1.48516682e-02  nfev=3662
  eta=1e-04: Re(alpha)=-0.97104087  Im(alpha)=1.36643754e-02  nfev=4034
  eta=1e-05: Re(alpha)=-0.97106360  Im(alpha)=1.35455265e-02  nfev=4382
  eta->0 (linear extrap. from eta=1e-4,1e-5): Re0=-0.97106613  Im0=1.35323211e-02
  rtol refinement (grid-doubling analog):
    rtol=1e-09: alpha=(-0.9710408792852563+0.01366437538117152j)  nfev=2762
    rtol=1e-10: alpha=(-0.9710408714230293+0.013664375392786346j)  nfev=4034
    rtol=1e-11: alpha=(-0.9710408704708181+0.01366437539686161j)  nfev=6338
    rtol=1e-12: alpha=(-0.9710408703623462+0.013664375396340975j)  nfev=9974
  rho_max sensitivity:
    rho_max=100.0: alpha=(-0.9708394625064171+0.013664355322675358j)
    rho_max=200.0: alpha=(-0.9710408714230293+0.013664375392786346j)
    rho_max=400.0: alpha=(-0.9710912322081781+0.013664380428818764j)
  rho_min sensitivity:
    rho_min=0.001: alpha=(-0.9710408718896343+0.013664378517982178j)
    rho_min=0.0001: alpha=(-0.9710408714230293+0.013664375392786346j)
    rho_min=1e-05: alpha=(-0.9710408713958077+0.013664375381739004j)

=== R2a beta=beta_star/2 (beta=0.03357565, rho_c=None) ===
  eta=1e-02: Re(alpha)=-0.50486515  Im(alpha)=5.81972798e-03  nfev=1412
  eta=1e-03: Re(alpha)=-0.50493584  Im(alpha)=5.82070530e-04  nfev=1406
  eta=1e-04: Re(alpha)=-0.50493655  Im(alpha)=5.82071508e-05  nfev=1406
  eta=1e-05: Re(alpha)=-0.50493655  Im(alpha)=5.82071518e-06  nfev=1406
  eta->0 (linear extrap. from eta=1e-4,1e-5): Re0=-0.50493655  Im0=1.10110530e-13
  rtol refinement (grid-doubling analog):
    rtol=1e-09: alpha=(-0.5049365536281237+5.820715150161027e-05j)  nfev=914
    rtol=1e-10: alpha=(-0.504936546425391+5.8207150793827586e-05j)  nfev=1406
    rtol=1e-11: alpha=(-0.5049365455368399+5.8207150713457936e-05j)  nfev=2198
    rtol=1e-12: alpha=(-0.5049365454354344+5.820715070071789e-05j)  nfev=3446
  rho_max sensitivity:
    rho_max=100.0: alpha=(-0.5048358407612968+5.8197080380925025e-05j)
    rho_max=200.0: alpha=(-0.504936546425391+5.8207150793827586e-05j)
    rho_max=400.0: alpha=(-0.5049617269176991+5.8209668803302586e-05j)
  rho_min sensitivity:
    rho_min=0.001: alpha=(-0.5049365464180773+5.8207150801165636e-05j)
    rho_min=0.0001: alpha=(-0.504936546425391+5.8207150793827586e-05j)
    rho_min=1e-05: alpha=(-0.5049365464242077+5.820715079330688e-05j)

=== R2b beta=2*beta_star (beta=0.13430261, rho_c=0.8918210192377795) ===
  eta=1e-02: Re(alpha)=-1.70386874  Im(alpha)=1.63556670e-01  nfev=2486
  eta=1e-03: Re(alpha)=-1.71348040  Im(alpha)=1.47889904e-01  nfev=2888
  eta=1e-04: Re(alpha)=-1.71444802  Im(alpha)=1.46310812e-01  nfev=3254
  eta=1e-05: Re(alpha)=-1.71454485  Im(alpha)=1.46152775e-01  nfev=3602
  eta->0 (linear extrap. from eta=1e-4,1e-5): Re0=-1.71455561  Im0=1.46135216e-01
  rtol refinement (grid-doubling analog):
    rtol=1e-09: alpha=(-1.7144480322561901+0.1463108112750957j)  nfev=2336
    rtol=1e-10: alpha=(-1.7144480232695616+0.14631081158600423j)  nfev=3254
    rtol=1e-11: alpha=(-1.7144480222227267+0.14631081161682036j)  nfev=5096
    rtol=1e-12: alpha=(-1.7144480221035652+0.1463108116210418j)  nfev=8006
  rho_max sensitivity:
    rho_max=100.0: alpha=(-1.714045213095403+0.1463107727799274j)
    rho_max=200.0: alpha=(-1.7144480232695616+0.14631081158600423j)
    rho_max=400.0: alpha=(-1.7145487442486915+0.14631082156938333j)
  rho_min sensitivity:
    rho_min=0.001: alpha=(-1.714448075891328+0.1463107847547844j)
    rho_min=0.0001: alpha=(-1.7144480232695616+0.14631081158600423j)
    rho_min=1e-05: alpha=(-1.7144480231622137+0.1463108110093789j)

Total compute time: 2.369s

=== R3: 2D Maxwell-Garnett / Clausius-Mossotti ===
Formula used: eps_eff = 1 + 2*n*alpha_pol/(1 - n*alpha_pol)
Inverted:      n = (eps_eff-1) / ( (eps_eff+1) * alpha_pol )
Using R1 (beta=beta_star) eta->0 extrapolated alpha = (-0.9710661264175603+0.013532321104294458j)
  |alpha| = 0.97116041   Re(alpha) = -0.97106613
  eps_eff=1.09326: n(using Re alpha)=-0.04588000   n(using |alpha|, signed)=-0.04587555   n(full complex)=(-0.04587109423874948-0.0006392380083672997j)
  eps_eff=1/1.09326=0.914695: n(using Re alpha)=0.04588000   n(using |alpha|, signed)=0.04587555   n(full complex)=(0.04587109423874947+0.0006392380083672995j)

```
