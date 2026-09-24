# The β-flow within Q_H=2: analytic derivation

Working notes toward Open Problem 4 (Paper III): testing whether the
still-undetermined scale ratio μ_UV/μ_IR in the α cascade can be built
from an action/energy difference *within* the Q_H=2 sector — the
condensate's feedback coupling β flowing from its bare value (β=0) to
its self-consistent fixed point (β=β*) — rather than a transition
between different Q_H sectors (which is the mechanism already used,
successfully, for the e^{4π} factor in m_e/Λ_cond).

Status: **derivation validated against established repo values; final
step (normalizing an energy difference into a scale-ratio logarithm)
still open.**

## 1. Setup

Battye–Sutcliffe profile (thin-torus, tube-only piece):
$$f(\rho) = 2\arctan(C/\rho), \qquad \mathrm{kern}(\rho) = (f')^2+\sin^2f/\rho^2 = \frac{2\sin^2f}{\rho^2}$$
using the locked-angle identity $(f')^2=\sin^2f/\rho^2$ (`P4:lem:locked_angle`).

Nondimensionalize with $x=\rho/C$ (so $C=1$ below; all quantities are
then pure numbers). Then $\sin f = 2x/(1+x^2)$ and
$$\mathrm{kern}(x) = \frac{8}{(1+x^2)^2}.$$

Paper I's own thin-torus integral definitions (`main_paper1.tex`,
"Thin-torus reduction"):
$$J_4 \approx (2\pi)^2R_0\int_0^\infty \sin^4\!f\,\frac{(f')^2}{\rho}\,d\rho \equiv (2\pi)^2R_0\,I_4,$$
$$J_a \approx (2\pi)^2R_0\int_0^\infty \sin^4\!f\left((f')^2+\frac{\sin^2f}{\rho^2}\right)\rho\,d\rho \equiv (2\pi)^2R_0\,I_{2a}.$$

The feedback functional (unweighted, parallel to $J_{2\mathrm{iso}}$):
$$J_{fb}(\beta) = \int_0^\infty \frac{\mathrm{kern}}{1+\beta\,\mathrm{kern}}\,\rho\,d\rho \equiv (2\pi)^2R_0\,I_{fb}(b), \qquad b\equiv\beta/C^2.$$

## 2. Closed forms (all derived by hand, then checked numerically —
## see `beta_flow_verify.py`)

**$I_{fb}(b)$.** Substituting $u=1+x^2$:
$$I_{fb}(b) = 4\int_1^\infty \frac{du}{u^2+8b} = \frac{4}{a}\arctan(a), \qquad a\equiv\sqrt{8b},$$
with $I_{fb}(0)=4$ (the $a\to0$ limit). **Matches `I_{2iso,tube}=4`
exactly** (established in `P4:lem:virial_factor`'s proof).

**$I_{2a}$** (β-independent, uses $\mathrm{kern}=2\sin^2f/\rho^2$ so the
integrand is $2\sin^6f/\rho$):
$$I_{2a} = 128\int_0^\infty \frac{x^5}{(1+x^2)^6}\,dx = 128\cdot\frac{1}{60} = \frac{32}{15}.$$
**Matches the established value exactly.**

**$I_4$** (β-independent, integrand $\sin^6f/\rho^3$):
$$I_4 = 64\int_0^\infty \frac{x^3}{(1+x^2)^6}\,dx = 64\cdot\frac{1}{40} = \frac{8}{5}.$$
**Ratio check:** $I_4/I_{2a} = (8/5)/(32/15) = 3/4$, matching
`P1:prop:bps`'s "BPS ratio" exactly. (Note: this 3/4 is *not* the
physical thick-torus target ratio $2^{4/3}/\varphi^5\approx0.227$ —
that discrepancy is already known and accounted for elsewhere in the
framework via the full 2D EL profile correction; it does not affect
this calculation, which only uses the simple radial ansatz for the
β-flow itself.)

## 3. Solving for b* — an independent check

The self-consistency condition is $J_{fb}(\beta^*)/J_a=\varphi$, i.e.
$I_{fb}(b^*) = \varphi\cdot(32/15) = 3.451806\ldots$. Solving
numerically:
$$b^* \equiv \beta^*/C^2 = 0.0671513\ldots$$

**This matches the framework's own quoted $b^*\approx0.06715$
(`P4:prop:bstar`) exactly.** This is an independent reconstruction —
built from the raw profile and Paper I's own stated integral
definitions, not looked up — landing precisely on an already-published
number. This validates the whole setup (profile, kern, integral
definitions, self-consistency condition) before trusting anything
built on top of it.

## 4. The energy functional across the flow

$$K_{fb}(\beta) = J_a + \mu^{*2}J_{fb}(\beta), \qquad \mu^{*2}=3-\varphi \text{ (established)},$$
$$E_{fb}(\beta) = K_{fb}(\beta)\cdot J_4.$$

Since $J_4$ and $J_a$ don't depend on $\beta$, define the reduced
(dimensionless, $(2\pi)^4R_0^2$-stripped) energy
$$\mathcal{E}(b) \equiv I_{2a} + (3-\varphi)\,I_{fb}(b).$$

At $b=0$: $\mathcal{E}(0) = \frac{32}{15}+4(3-\varphi) = 7.6612\ldots$

At $b=b^*$: $\mathcal{E}(b^*) = \frac{32}{15}+(3-\varphi)\varphi\cdot\frac{32}{15} = \frac{32}{15}(1+3\varphi-\varphi^2)$.
Using $\varphi^2=\varphi+1$: $1+3\varphi-\varphi^2 = 1+3\varphi-\varphi-1 = 2\varphi$, so
$$\mathcal{E}(b^*) = \frac{64\varphi}{15} = 6.9036\ldots$$

**This is exactly $2\varphi\cdot I_{2a}$ — i.e. it independently
reproduces the already-proved identity $K_{fb}=2\varphi J_a$** (Paper
I, Proposition 2.5). Another consistency check passed, not a new
result — but confirms the reduced-energy bookkeeping is right.

## 5. Where this stands: the magnitude problem

The energy *decreases* from bare to self-consistent
($\mathcal{E}(0)=7.66 \to \mathcal{E}(b^*)=6.90$), consistent with
$\beta^*$ being an energy-minimizing fixed point — a sensible sign.

But the **size** of this change is far too small to plausibly reduce
to the target log-ratio:
$$\Delta\mathcal{E} = \mathcal{E}(0)-\mathcal{E}(b^*) = 0.7576,\qquad
\ln\!\big(\mathcal{E}(0)/\mathcal{E}(b^*)\big) = 0.1041,$$
versus the needed $\log(\mu_{UV}/\mu_{IR})\approx26.86$ (Δ₁'s
requirement). Neither the raw difference nor the log-ratio of these
dimensionless reduced energies is anywhere near $26.86$ — off by two
orders of magnitude either way.

**This doesn't kill the "β-flow within Q_H=2" idea, but it does rule
out the naive version of it**: a bare energy difference of the reduced
functional $\mathcal{E}(b)$, taken at face value, is not the quantity
that produces Δ₁. If this mechanism is right at all, either (a) there's
a large physical prefactor still missing (analogous to how $\Delta
S_{CS}=4\pi$ needed an actual Chern–Simons normalization, not just a
raw energy functional value), or (b) the correct quantity to
exponentiate isn't $\mathcal{E}(0)-\mathcal{E}(b^*)$ but something else
built from the same flow (e.g. an integrated "RG time" $\int
d\beta\,(\partial S/\partial\beta)$ rather than the endpoint
difference, or a quantity involving $Q=10$ or $R_0=3$ multiplicatively).

**Next step, if pursued:** identify the actual Euclidean-action
normalization convention used to get $\Delta S_{CS}=4\pi$ in Paper IV
(`P4:thm:second_spoke`) and check whether the same normalization,
applied to $E_{fb}(\beta)$ across this flow, changes the order of
magnitude enough to approach $26.86$.
