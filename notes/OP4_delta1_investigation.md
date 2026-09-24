# Open Problem 4 (Paper III): deriving Δ₁ / μ_UV/μ_IR — working notes

Status: **open**. This file tracks a brainstorming session investigating whether
the golden-angle numerator (`360 = k·|2I| = 3×120`) has a group-theoretic origin,
and whether that connects to deriving the still-unresolved scale ratio in the
α cascade. Nothing below is established; this is a record of hypotheses tried
and their outcomes, so future sessions don't re-walk the same dead ends.

## The target

`P3:prop:cascade` (foundations/alpha.tex, main_paper3.tex):

```
φ⁶/sin⁴(π/5) = 150.332  --Δ₁=-12.824-->  360/φ² = 137.508  --Δ₂=-0.472-->  α⁻¹_phys = 137.036
```

Δ₂ = k/(2π) is derived (WZW one-loop, exact via Knizhnik–Zamolodchikov).
Δ₁ is **not** derived: it requires `log(μ_UV/μ_IR) = Δ₁·2π/k`, and solving
numerically gives:

```
log(μ_UV/μ_IR) = 26.86   (corrected from an internally-inconsistent "≈e^4.27"
                           previously in main_paper3.tex — now fixed)
μ_UV/μ_IR ≈ 4.62×10^11
```

This is Open Problem 4. **Finding: it is not actually resolved in Paper IV**,
despite main_paper3.tex's open-problems appendix claiming so (see "Dead ends"
below) — this has been corrected in main_paper3.tex.

## Established anchors available to build from

- `P1:thm:k3`: k=3 forced by ℤ₁₀ ⊂ 2I (published).
- `P1:thm:vacuum`: |2I|=120, established as a physical fact about the
  condensate ground state, not a bookkeeping device (published).
- `P3:thm:verlinde` / `P3:thm:bos_main`: V* = φ, the condensate's
  self-consistent virial fixed point, proved via the Verlinde S-matrix
  (published, exact). This likely already answers "why is the vacuum
  φ-ordered" — no separate derivation needed.
- `P4:thm:spoke_ground`: Q_H=2 is the singlet bound state of two Q_H=1
  quanta via j_{1/2}×j_{1/2}=j_0⊕j_1 (published).
- `P4:thm:spoke_mass_spectrum`: three generations (j=1/2,1,3/2), strictly
  mass-ordered, each with its own T-matrix phase (published) — a genuine,
  already-proved 3-element φ-ordered "tower."
- `P4:thm:compton` (Compton-wavelength identification): β* = 1/m_WZW².
- `P4:thm:mass_renorm`: m_WZW²/m₀² = (8φ/15)·(2R₀²/(2R₀²+1)), exact at
  physical R₀=3 using the locked Hopf angle (θ=π/4) and winding-virial
  factorization lemmas.
- `P3:thm:bos_3d` / `P3:cor:bos_valid`: Witten bosonization proves the
  classical geometric (soliton) picture and the WZW quantum picture are
  the *same physics*, up to O(R₀⁻²) — this is the standing bridge any
  "geometric ↔ group-theoretic" unification argument should use.
- `P3:rem:packing_layers` (newly labeled/extracted this session): the
  golden angle is described in the source prose as "the angle at which
  successive icosahedral packing layers are most efficiently arranged;
  the angle of minimum geometric stress in a φ-ordered condensate." This
  was previously un-extractable (no label/status/source tags) — now fixed
  in both main_paper3.tex and foundations/alpha.tex.

## Hypotheses tried, and their status

### H1 — static count: 360 = k×|2I| as "number of gauge orientations"
An earlier (external/pasted) attempt tried deriving α⁻¹ = k|2I|/dim_q² from
WZW current correlators, orbifold charge quantization, Verlinde sums, and
Chern–Simons Wilson loops. **All failed** — each gives the wrong power of
|2I| (bare correlator: no |2I| dependence at all; charge quantization:
|2I|² not |2I|; Verlinde/CS: no clean match).

**Ruled out on physical grounds, independent of the failed calculations**:
α runs with scale (real, established fact — the whole cascade exists
because α is scale-dependent). A static combinatorial count (a fixed
number of orientations/copies/microstates) is scale-independent by
construction and cannot be the mechanism, regardless of whether some
calculation can be made to numerically match it. This was the key insight
that redirected the whole investigation toward RG-flow interpretations.

### H2 — Q_H=1 charge-carrier framing (embedding degeneracy)
Considered and explicitly set aside by the user: "geometrical and
group-theoretical, without dropping to charge carriers." Not pursued
further.

### H3 — k copies of |2I|, one per generation
If "layer" = "generation," `k×|2I|` reads as "one full copy of the
icosahedral structure per generation" (k=3 already independently forced
as generation count, `P4:thm:three_gen`). Plausible fit to
`P3:rem:packing_layers`'s "successive layers" language, since
`P4:thm:spoke_mass_spectrum` gives an actual established 3-element
φ-ordered tower (generations), unlike:
  - the φ-spiral (P11) — explicitly a *display choice*, not physical
    (`P11:rem:Narms_choice`), ruled out as a candidate for "layers."
  - the Q_H=0,2,4,6,8 spin tower (P13) — a different index (angular
    momentum of higher composite charge sectors), not obviously the
    same "layers" as the golden-angle remark.
Status: plausible, **not tested against the numbers**, superseded in
priority by H5 below once the "α runs" constraint was raised.

### H4 — intrinsic 3-fold axes of 2I itself
Alternative to H3: the icosahedron has 10 native 3-fold rotation axes;
"k=3" might reflect structure *inside* 2I rather than 3 external copies
of it. Not distinguished from H3 — would need to identify which
observable actually counts one way vs the other. Not tested further.

### H5 — RG-flow / discrete-step reframing (current leading direction)
Following the "α runs, so it can't be a static count" objection: `k×|2I|`
would have to enter as a **rate, step-count, or trajectory length** in
the RG flow, not a snapshot number. Since Δ₂'s slope is bare k (no |2I|),
the natural place for |2I| to enter is the **length of the flow** — i.e.
literally the still-open μ_UV/μ_IR ratio, or the number of discrete
steps/images between the UV condensate scale and the IR icosahedral
fixed point. This reframes "derive k×|2I|" as *equivalent to* solving
Open Problem 4 correctly, rather than a separate formula to derive
independently. **This is the live thread.**

## Numerology checks run (both ruled out — dead ends, don't re-chase)

1. **Old paper text claimed μ_UV/μ_IR ≈ e^4.27.** Recomputing directly
   from Δ₁=(k/2π)ln(ratio) with k=3 gives ln(ratio)=26.86, not 4.27 — the
   old figure didn't follow from the framework's own stated equation.
   **Fixed in main_paper3.tex** (now reads `≈26.86`, `4.6×10^11`).
   Before the fix was known, `ln(72)=4.2767` (0.16% off the old, wrong
   4.27) was found to equal `k×|2T|` (|2T|=24, the *Q_H=3 baryon sector's*
   McKay group, not 2I) — an intriguing near-match, but it was a match to
   an erroneous number and is now moot. **Dead end.**

2. **Thin/thick-torus boundary-flux ratio vs α_geo.** From the Pohozaev
   variation-(d) calculation (main_paper3.tex §"Radial dilation"): thin-
   torus estimate `4π²·2/R₀ ≈ 26.32` (Q_math=2 ansatz) vs the value
   required to satisfy the normalization conjecture, `B_required(R₀=3) ≈
   0.175` (using Q=10). Their ratio, 150.23, is close to `α⁻¹_geo =
   φ⁶/sin⁴(π/5) = 150.332` (0.067% off). **Ruled out**: a 0.067% "match"
   is coincidence-sized, not identity-sized (real identities in this
   framework match to machine precision or are exactly provable); no
   algebraic path connects π²·2^(2/3)·φ^(23/2) (boundary-flux side) to
   φ⁶/sin⁴(π/5) (WZW-angle side); and both source numbers are already
   individually flagged as unreliable in their own context (thin-torus
   formula explicitly "fails by ≈150×"; B_required is the still-open
   piece of a different, separate conjecture). **Dead end, don't re-chase.**

## Dead end: the "Resolved in Paper IV" cross-reference

`main_paper3.tex`'s Open Problems appendix claimed (before this session's
edit) that deriving Δ₁/μ_UV/μ_IR was "Resolved in Paper IV... via the
Compton-wavelength identification and one-loop mass renormalization
(Theorems 3.1 and 4.3 therein)." **Checked directly against Paper IV and
found false**: `P4:thm:compton` proves β*=1/m_WZW² (nothing about μ_UV/
μ_IR); `P4:thm:mass_renorm` proves m_WZW²/m₀² (also nothing about μ_UV/
μ_IR). A full grep of main_paper4.tex for `12.824`, `150.332`, `Δ_1`
turns up the golden angle only as a *given input*, never derived. Paper
IV's actual route to higher precision is a different, independent
mechanism (T-matrix correction giving the three-term formula, then a
Chern-Simons correction giving the four-term formula) that bypasses
Δ₁'s derivation rather than supplying it. **Corrected in main_paper3.tex
this session** — the bullet has been reopened/reworded to state this
honestly.

## Resolved sub-question: the origin of the numerator 360

**$360 = Q\times(\pi/5\text{ in degrees}) = 10\times36°$, exactly.** $Q=10$
is the independently-established topological Hopf charge
(`P1:thm:charge_q10`); $36°=\pi/5$ is the *same* pentagon angle already
used in $\alpha^{-1}_{\rm geo}=\varphi^6/\sin^4(\pi/5)$ — not a separate
assumption. This is a genuinely different, better-grounded finding than
the $k\times|2I|=360$ conjecture (which never found a mechanism): it
uses only quantities already established elsewhere in the framework,
and connects $360$ directly to the angle already sitting in the
"solid ground" term, rather than to an unrelated group-order product.

This resolves the "why $360$, and is it just the Babylonian degree
convention" concern raised earlier (checked at the time: the
$\alpha^{-1}\approx360/\varphi^2$ match fails completely in radians or
turns, which was concerning if $360$ had no independent meaning). It
does not resolve $\Delta_1$ itself — knowing why $360$ is the right
numerator doesn't explain why $\varphi^6/\sin^4(\pi/5)$ comes out
numerically close to $360/\varphi^2\cdot(\ldots)$. But it does mean the
golden angle and $\alpha^{-1}_{\rm geo}$ are now both expressed using
*only* $\varphi$, $Q$, and the single angle $\pi/5$ — no free-floating
"$360$" left unexplained.

**Corrected in both `papers/main_paper3.tex` and `foundations/alpha.tex`**
(`P3:rem:no_fitted_params`, `P3:rem:why_golden_angle`) to state this
rather than "$360$ is the number of degrees in a full rotation."

## Live lead, tightened: $\Delta_1 \approx 4\pi + (15/8-\varphi)$

Starting point: $\Delta_1=12.824$ vs.\ $4\pi=12.566$ — a $2.05\%$ miss
on its own, but the first candidate all session in the *right order of
magnitude* (everything else was off by 10–100$\times$) and built from
an already-established quantity ($4\pi$ is the same Chern–Simons unit
that produces $m_e/\Lambda_{\rm cond}$ via $\Delta Q=1$), not one
invented to fit.

**Tightened significantly**: the residual $\Delta_1-4\pi=0.2577$ matches
$15/8-\varphi=0.2570$ almost exactly. Combined:
$$\Delta_1 \approx 4\pi + \left(\frac{15}{8}-\varphi\right) = 12.823337,
\quad\text{vs. actual } 12.824120 \;\; (0.0061\%\text{ match}).$$
Both terms are already-established, independently-derived quantities:
$4\pi$ (CS action, as above) and $15/8-\varphi$ (the exact distance the
condensate's virial function travels from its bare thin-torus value
$V(0)=15/8$ to the self-consistent fixed point $V(\beta^*)=\varphi$,
`P4:lem:virial_factor`, `P3:thm:V_exact_all_R0`). $0.0061\%$ is tighter
than several of the framework's own accepted physical residuals.

**Exact closed form derived**: $\Delta_1 = (2136\varphi-3392)/5$ (from
$\varphi^6=8\varphi+5$, $(3-\varphi)^2=5(2-\varphi)=5/\varphi^2$). This
proves $\Delta_1$ is algebraic of degree 2 over $\mathbb{Q}$ (an element
of $\mathbb{Q}(\sqrt5)$). Since $\pi$ is transcendental (Lindemann) and
algebraic numbers are closed under addition, $4\pi+(15/8-\varphi)$ is
provably transcendental. **An algebraic number and a transcendental
number can never be exactly equal, at any precision** — so the
$0.0061\%$ match above is, by proof, not an exact identity and never
can be. Checked whether the *log-ratio* $\Delta_1\cdot2\pi/k=26.86$
(rather than $\Delta_1$ itself) might be the right transcendental
target instead — it isn't: that quantity is $\pi$ times $\Delta_1$'s
own algebraic coefficients by construction ($\frac{2\pi}{15}(2136\varphi-3392)$),
so it's not an independent quantity, just a rescaling, and doesn't
resolve the mismatch. If a real derivation exists, it must either treat
the $4\pi+(15/8-\varphi)$ match as a genuine (unexplained) numerical
coincidence, or identify a *different*, adjacent quantity — not
$\Delta_1$ in its current exact form — as the true transcendental
target. Not yet identified.

**Consistency check performed**: using the $R_0=3$-corrected virial
$V(0;R_0)=\frac{15}{8}(1+\frac{1}{2R_0^2})$ instead of the thin-torus
value *worsens* the match to $0.81\%$. This is not evidence against the
candidate: $\Delta_1$ itself has zero $R_0$-dependence (both
$\varphi^6/\sin^4(\pi/5)$ and $360/\varphi^2$ are pure, $R_0$-independent
$\varphi$-expressions), so a match that holds at the thin-torus
(also $R_0$-independent) level and degrades once $R_0$-dependence is
introduced is exactly the pattern expected if the candidate is real.

**Proposed physical reading** (not yet a derivation): both terms are
properties of the vacuum/condensate itself, not of any specific soliton
state within it — consistent with the "Hopfion is the smoke, the
condensate is the air" framing from early in this investigation. $4\pi$
is the topological cost of a charge-sector transition in the vacuum's
own configuration space; $15/8-\varphi$ is the distance the vacuum's own
virial ratio travels during relaxation to self-consistency. Under this
reading, $\Delta_1$ measures two aspects of the vacuum's total behavior,
with the observed electromagnetic coupling as a downstream consequence.

**Status: unverified, logged in `main_paper3.tex`** (new remark
`P3:rem:delta1_candidate`, inserted directly after `P3:prop:cascade`,
explicitly flagged as a numerically-confirmed candidate sum, not a
derivation). **Next step, not yet attempted**: look for a single,
unified vacuum-level calculation (a combined free energy or effective
action spanning both the topological Chern–Simons sector structure and
the classical relaxation dynamics) that would produce
$4\pi+(15/8-\varphi)$ as one expression rather than two separately
derived pieces added together.

**Secondary thread (not yet resolved)**: whether the "$4\pi$" piece
could be a *second*, distinct $\Delta Q=1$ topological transition from
the one already used for $m_e$ (candidate: Paper XVII's conjectured
$Q_H=1$ [neutrino] $\to Q_H=2$ [charged lepton] transition, specifically
about acquiring electric charge). Caveat: $\Delta S_{CS}=4\pi\Delta Q$ is
purely topological (depends only on charge numbers, not species), so if
this is "the same transition" physically as the $m_e$ one, it is not
actually independent — unresolved. Also flagged: Paper IV calls its
$Q_H=1$ state "the electron sector" while Paper XVII (a conjecture)
calls $Q_H=1$ "the neutrino" — an unresolved inconsistency in how
$Q_H=1$ is physically identified across papers.

## Unified physical picture attempted: WZW action (kinetic + WZ split)

Motivation: $4\pi$ (SU(2) double cover, 2 full rotations) and $15/8$
(structurally meaningful — the bare thin-torus virial ratio, a real
geometric fact about the BS profile) and $\varphi$ (fixed point of the
self-referential map $x\mapsto1+1/x$, and literally what the
condensate's own feedback loop converges to) suggest a single coherent
story rather than an arbitrary sum:

- **Kinetic piece, no new derivation needed**: $\int_0^{\beta^*}dV =
  V(\beta^*)-V(0) = \varphi-15/8$, trivially the fundamental theorem of
  calculus applied to the virial function's known endpoints.
- **Topological piece, a specific hypothesis**: the "bare" ($\beta=0$,
  $V=15/8$) and "self-consistent" ($\beta=\beta^*$, $V=\varphi$) states
  are hypothesized to be the *unbound* and *bound* configurations of the
  **same** $Q_H{=}1\to Q_H{=}2$ fusion already used for $m_e$
  (`P4:thm:second_spoke`, `P4:thm:spoke_ground`) — not a second,
  independent transition. Under this reading $4\pi$ is the same
  physical event observed through a different quantity (the vacuum's
  virial ratio, rather than the electron's mass) — one microphysical
  event, two observable footprints. This resolves the "is this
  double-counting" caveat from the secondary thread above, if correct.
- **Why additive, structurally**: the standard $\mathrm{SU}(2)_k$ WZW
  action itself splits as $S_{WZW}=S_{\rm kinetic}+\frac{k}{12\pi}\int_B
  \mathrm{Tr}(g^{-1}dg)^3$ (smooth sigma-model term plus topological WZ
  term, added not multiplied) — exactly the shape being proposed here,
  and this framework already claims to *be* the $\mathrm{SU}(2)_3$ WZW
  model (`P3:thm:bosonization`), so this isn't an imported analogy.

Status: a coherent, motivated physical *picture*, not a derivation. Not
shown: that the WZW action, evaluated rigorously on this specific
trajectory, actually produces $4\pi+(15/8-\varphi)$ (or something that
reduces to it) rather than merely resembling that shape.

## Refinement: the $4\pi$ piece should probably reduce via a trig identity, not stand as a bare transcendental

Key observation: $\pi$ is transcendental, but $\sin$/$\cos$ of a
*rational multiple* of $\pi$ is always algebraic. This is not
hypothetical — it is the exact mechanism **already at work** in this
derivation: $\sin^4(\pi/5)$ reduces exactly to $(3-\varphi)^2/16$ (no
$\pi$ survives), and `P4:rem:sign_pattern` already establishes
$\cos(2\pi QT_1)=\cos(2\pi QT_{3/2})=0$ exactly for the established
$T$-matrix phases. So the "$4\pi$" in the candidate almost certainly
should not appear as a bare additive/exponential transcendental term —
if the WZW-action picture above is right, the topological contribution
should appear *inside* a trig/phase expression that reduces
algebraically, the same way the kinetic term already does.

**Checked**: $\cos(2\pi QT_j)$ for $j=1/2,1,3/2$ (using $T_{1/2}=3/40$,
$T_1=13/40$, $T_{3/2}=27/40$, established) all vanish identically —
this is a structural/arithmetic fact ($QT_j$ is always a quarter-integer
for these $T_j$, giving phase $\equiv\pi/2\pmod\pi$), not new
information, and doesn't provide a nonzero target to match.

**Computed the exact algebraic value the "topological piece" must take**,
assuming the kinetic piece is exactly $15/8-\varphi$ (forced, since
$\Delta_1$ is exactly known):
$$\text{topological}_{\rm exact} = \Delta_1 - (15/8-\varphi)
= \frac{17128\varphi-27211}{40} = 12.567154\ldots$$
This is provably algebraic (difference of two elements of
$\mathbb{Q}(\varphi)$), so it cannot equal $4\pi$ exactly — checked
against candidates:
- $4\pi = 12.566371$ (diff $0.0008$, $0.006\%$)
- $24\varphi^2/5 = 12.566563$ (diff $0.0006$, $0.005\%$) — found via
  noting $\text{topological}_{\rm exact}\cdot\sin^4(\pi/5)\approx3/2$,
  i.e. $\text{topological}_{\rm exact}\approx\frac{3/2}{\sin^4(\pi/5)}=
  \frac{24}{(3-\varphi)^2}=\frac{24\varphi^2}{5}$ (using
  $(3-\varphi)^2=5/\varphi^2$).

Neither is exact. Two competing "shapes" were considered: a
transcendental one ($4\pi$-flavored, physically motivated by the CS
action, provably cannot be exact) or an algebraic one
($24\varphi^2/5$-flavored, which *could* in principle be exact).

**$24\varphi^2/5$ investigated and downgraded.** Found the exact
rewriting $24\varphi^2/5 = k\cdot I_4\cdot\varphi^2$ ($k=3$,
$I_4=8/5$, both established profile/level quantities — verified with
exact fractions). Checked whether $k\cdot I_4$ alone (without $\varphi^2$)
has independent meaning: found $k\cdot I_4 = c+k$ ($c=3k/(k+2)=9/5$,
the central charge) — but this is **circular, not a finding**: $I_4=8/5$
is a fixed classical-profile constant with no $k$-dependence in its own
derivation, while $c/k+1=3/(k+2)+1$ only equals $8/5$ because $k=3$ is
plugged in; the "identity" $I_4=1+c/k$ reduces algebraically to "$k=3$",
which was already known. Same species of trap as the $360=k|2I|$
substitution into $\Delta_1$ earlier in this document. **Then checked
the $\varphi^2$ factor itself**: attempted to justify it via "$\varphi^2$
as the quantum-mechanical weight of the charged sector" (Born-rule-style,
paralleling $360/\varphi^2$'s own denominator). Traced this framing back
to its actual source: the original external/pasted H1 attempt at the
very start of this investigation, which relabeled $\varphi^2$ (already
present in the golden angle purely as the standard phyllotaxis
convention, $360°\times(2-\varphi)$, nothing to do with quantum
mechanics) as "$\dim_q^2$" and added the "quantum weight" gloss without
derivation. That framing was never established anywhere in the actual
papers (checked: `dim_q`/$\varphi$ never appears as a Born-rule weight
in any established Verlinde sum or partition function in Papers I, III,
IV, or XVII) and was already flagged unreliable in H1. Reusing it here to
justify $\varphi^2$ (rather than $\varphi^1$) in $k\cdot I_4\cdot\varphi^2$
imports that same ungrounded premise.

**Net effect**: $k\cdot I_4\cdot\varphi^2=24\varphi^2/5$ remains an exact
arithmetic identity (real, verified) but has **no surviving physical
justification** — downgraded to the same tier as the numerology already
killed elsewhere in this document ($\ln(72)=k|2T|$, $\beta^*\times400$):
a tight numerical near-miss with no mechanism, not a physically
motivated candidate. The $4\pi+(15/8-\varphi)$ candidate is unaffected
by this downgrade — its grounding (established CS action; established
virial relaxation) never relied on any Born-rule/quantum-weight framing.

## Resolution: the narrative was wrong, not the search

**$\pi/5$ is the pentagon angle (WZW-modular, fixed by $k$), not a
Hopf-polar location.** Checked directly: the only established,
physically-derived use of the suppression function $S(\theta)=
\sin^4\theta/\varphi^6$ (`P1:eq:suppression`) is in the Weinberg-angle
derivation (`P1:thm:amplitude_ratio`), where $\theta$ is a genuine
Hopf-fibration polar angle and the physically-justified value is
$\theta=\pi/2$ (the equatorial vacuum locus, `P1:thm:vacuum`) — not
$\pi/5$. $\pi/5=\pi/(k+2)$ is instead the WZW modular parameter
($q=e^{i\pi/(k+2)}$, `P1:thm:k3`), a fixed group-theoretic quantity, not
a dynamical location subject to any geometric ($R_0$-dependent)
correction. This retracts H8's $\theta'$ finding (was fitting a
continuous parameter, same species as the numerology already killed).

**Exact identity found**: $\sin^2(\pi/5)=(3-\varphi)/4$ (from
$\cos(\pi/5)=\varphi/2$ and $\varphi^2=\varphi+1$). And $3-\varphi=\mu^*
=\sqrt{k+2}/\varphi$ is *exact for the full physical profile at every
$R_0$* — not a BPS-ansatz leading-order approximation awaiting
correction (unlike $J_4/J_{2a}$, which genuinely does have that
BPS-vs-physical gap). This is established via `P3:thm:V_exact_all_R0`,
a fixed-point existence theorem: $V^*=J_{fb}/J_a=\varphi$ holds exactly
at every finite $R_0$, and the LSC-equivalence theorem
(`P3:thm:lsc_equivalence`) shows $\mu^*=3-\varphi$ follows from this
with no approximation. So:

$$\alpha^{-1}_{\rm geo} = \frac{\varphi^6}{\sin^4(\pi/5)} = \frac{16\varphi^6}{(3-\varphi)^2}$$

is **fully exact**, built from two independently-exact ingredients
($\varphi^6=\lambda$, $\mu^*=3-\varphi$). The golden angle $360/\varphi^2$
is *also* fully exact (a geometric/packing constant). Checked whether
there's a hidden algebraic identity making these equal (cross-multiplying
using $\varphi^8=21\varphi+13$): there isn't one — $336\varphi+208\neq
3600-1800\varphi$. So **$\Delta_1$ is the genuine, exact, computable
difference between two independently-derived exact quantities that
share the same $\varphi$/$k{=}3$ origin — not an RG-running effect, not
a scale ratio, not a residual awaiting a derivation.**

**Conclusion**: Paper III's "WZW anti-screening running from the
condensate UV scale to the golden-angle attractor" framing for $\Delta_1$
was a narrative laid over an unexplained numerical gap, not a derived
mechanism — consistent with every mechanism this session tried (H1
through H8) failing for a structural reason. **Narrative corrected in
`papers/main_paper3.tex` and `foundations/alpha.tex`** (this session):
`P3:prop:cascade`, the "WZW anti-screening mechanism" subsection (now
"WZW one-loop correction"), `P3:rem:residual_status`,
`P3:rem:two_gaps_independent`, the abstract/intro cascade diagrams, the
`P3:sec:sign` anti-screening-sign section, the open-problems list, and
the Summary section were all corrected to state plainly: only $\Delta_2=
k/(2\pi)$ (golden angle $\to$ physical value) is a genuine WZW
running/anti-screening effect; $\Delta_1$ ($\alpha^{-1}_{\rm geo}\to$
golden angle) is an exact difference between two exact quantities, not
a running effect, with Open Problem 4 reframed accordingly (from
"derive the scale ratio" to "is there a deeper reason these two exact
quantities are close, or is it coincidence").

**Not yet done**: `main_paper4.tex` and `main_reader_guide.tex` also
contain instances of the same conflation (e.g.
`main_paper4.tex:688`, "the WZW anti-screening that generated $\Delta_1$
and $\Delta_2$," and a table row spanning $\Lambda_{\rm cond}\to m_e$
as one continuous anti-screening effect) — these were not corrected
this session and remain inconsistent with the corrected Paper III
narrative. Flagged for a future pass.

## Next direction (H5, in progress)

Test whether the Compton-wavelength identification (β*=1/m_WZW²) and the
mass-renormalization theorem give a natural candidate for μ_UV and μ_IR
as actual physical mass/energy scales of the condensate (candidates:
Λ_cond, m₀, m_WZW, m_e), and check numerically whether any combination's
log-ratio reproduces 26.86, or relates to it via the already-known R₀=3
thick-torus correction factors (18/19, 48φ/95, etc.) or φ-powers.

**Tested and ruled out**: μ_UV=Λ_cond, μ_IR=m_WZW (via the spoke formula
chain) gives ln-ratio=22.24, missing the target by 4.62 in the log (a
factor of ~101) — a real, quantified miss, not close enough to call a
match. A tempting but unexplained near-match to ln(Q²)=ln(100)=4.605 on
that gap was noted but flagged as unverified pattern-matching with no
mechanism (same category as the numerology already killed below) —
raised in confidence by finding the *old, erroneous* "e^4.27" figure
that used to be in main_paper3.tex was itself just a units slip
(Δ₁/k, missing the 2π factor) rather than anything physical, which is
a caution against taking unexplained near-matches seriously without a
mechanism.

**Reading 2 test (geometric length-ratio, R₀ vs C\*)**: ruled out on
order-of-magnitude grounds. $\ln(R_0/C^*)=0.50$ (using the established
$C^{*2}=3\varphi^5/(8\cdot2^{1/3})$) — off from the target by a factor
of ~54 in the log. Every genuinely huge hierarchy in this framework
comes from an *exponentiated action* (e.g. $e^{4\pi}$ in $m_e/\Lambda_{\rm cond}$),
never a bare ratio of comparable-sized geometric lengths — so this
rules out $\mathcal{B}(R_0)$-type quantities even more firmly than the
earlier numerical near-miss did (wrong category of object, not just a
coincidental miss).

## H6 — β-flow within Q_H=2 (same topology, geometry changes)

Distinct from the $Q_H=1\to Q_H=2$ topology-changing CS action that
gives $e^{4\pi}$ (`P4:thm:second_spoke`): this hypothesis says the
$\mu_{UV}/\mu_{IR}$ ratio might come from an action difference *within*
the fixed $Q_H=2$ sector, as the feedback coupling β flows from its bare
value (β=0) to its self-consistent fixed point (β=β*=1/m_WZW²,
`P4:thm:compton`).

**Derivation and verification**: see `notes/beta_flow_derivation.md` and
`notes/beta_flow_verify.py` (run: exact match on all checks).
Built $I_{fb}(b)=\frac{4}{\sqrt{8b}}\arctan(\sqrt{8b})$ from scratch
(BS profile, Paper I's own thin-torus integral definitions) and
validated it three independent ways against numbers already in the
repo: $I_{fb}(0)=4$ (matches `I_{2iso,tube}`), $I_{2a}=32/15$ and
$I_4=8/5$ (both exact, giving the established BPS ratio $I_4/I_{2a}=3/4$),
and — the strongest check — solving $I_{fb}(b^*)/I_{2a}=\varphi$ gives
$b^*=0.0671513$, an **exact independent reproduction of the framework's
own quoted $b^*\approx0.06715$** (`P4:prop:bstar`). The reduced energy
$\mathcal{E}(b^*)=2\varphi\, I_{2a}$ also independently reproduces the
established identity $K_{fb}=2\varphi J_a$.

**Result: ruled out in its naive/bare form.** $\mathcal{E}(0)=7.66\to
\mathcal{E}(b^*)=6.90$ (energy decreases, sensible sign for an
EL-minimizing fixed point), but the magnitude is ~2 orders of magnitude
too small: $\Delta\mathcal{E}=0.76$, $\ln(\mathcal{E}(0)/\mathcal{E}(b^*))=0.10$,
versus the needed $26.86$. Doesn't kill the idea outright — it kills
"bare energy difference of the reduced functional, no further
normalization." If pursued further, the next step is finding the actual
Euclidean-action normalization convention behind $\Delta S_{CS}=4\pi$
(`P4:thm:second_spoke`) and checking whether applying the same
normalization to $E_{fb}(\beta)$ across this flow closes the ~2-orders
gap.

**Structural problem found for H6/Route B, likely fatal in this form.**
Traced exactly how $\Delta S_{CS}=4\pi$ is derived (`P4:eq:delta_cs`):
it comes from the **Chern–Simons topological term** $\int A\wedge dA=
16\pi^2Q$, which depends *only* on the integer Hopf charge $Q$, not on
any continuous profile parameter. $\Delta S_{CS}=4\pi\cdot\Delta Q$.
Every large exponential hierarchy in this framework ($e^{4\pi}$,
$e^{49\pi/6}$) comes from exactly this mechanism, scaled by $\Delta Q$.
The $\beta{=}0\to\beta^*$ flow holds $Q_H=2$ fixed ($\Delta Q=0$), so
this mechanism gives *zero* by construction — consistent with every
continuous-parameter effect in the framework being small ($O(R_0^{-2})$,
a few percent at most, never a large hierarchy). **Route B, in the
"action difference at fixed topology" form, is structurally disfavored**,
not just numerically small by chance — there is no precedent anywhere
in the repo for a fixed-topology continuous deformation producing a
large effect. Would need an entirely different, non-CS, non-perturbative
mechanism to revive it.

## H7 (Route A revisited) — chameleon screening machinery exists, wrong pair of scales

Route A's original idea (μ as a literal radial position/local-density
scale within the tube cross-section, holographic-RG style) turns out to
have real, already-published machinery behind it:
`quantum/chameleon.tex` has an established three-scale hierarchy using
the screening function $\varepsilon(\rho)=1/[\varphi^6(1+\beta\rho)]$ —
the *same* mathematical form as $\mathcal{E}(\beta)$ from H6 — evaluated
at three concrete environments:

| Environment | $\beta\rho$ | $\varepsilon$ |
|---|---|---|
| Cosmic attractor | $\varphi$ (exact) | $1/\varphi^8\approx0.021$ |
| Bohr radius (atomic) | $\approx2.6\times10^{25}$ | $\approx2.2\times10^{-27}$ |
| QCD/nucleon scale | $\approx9.3\times10^{48}$ | $\approx6\times10^{-51}$ |

Note the cosmic-attractor row uses $\beta\rho_\infty=\varphi$ — the
*same* self-consistency condition independently reconstructed in H6
($I_{fb}(\beta^*)/I_{2a}=\varphi$). Same fixed point, same formula shape.

**Checked and ruled out**: the established cosmic-vs-atomic ratio
($\beta\rho_{\rm atom}/\beta\rho_\infty\approx1.6\times10^{25}$, or its
4th root $\sim2\times10^6$, since these $\rho$'s scale as mass$^4/
\Lambda_{\rm cond}^4$) does **not** match the target $4.6\times10^{11}$
in either form (off by 14 orders of magnitude raw, or 5 orders as a 4th
root). This existing hierarchy screens across *different macroscopic
environments* (cosmic vacuum / atomic / nuclear density) — it is not
about position *within* a single soliton's cross-section, so this is
the wrong pair of scales, not a failure of the mechanism itself.

**What's actually required, now well-posed**: build the analogous map
using $\mathrm{kern}(\rho)$ (already derived exactly in H6:
$\mathrm{kern}(\rho)=8C^2/(\rho^2+C^2)^2$) as the *local* "ambient
density" in the same screening formula, instead of an external
environmental density — then ask what radial locations $\rho_{UV}$,
$\rho_{IR}$ within the *single* $Q_H=2$ profile correspond to the two
$\alpha$-cascade endpoints ($\alpha^{-1}_{\rm geo}=150.332$ and the
golden angle $137.508$), and check whether the resulting radial/scale
ratio reproduces $\Delta_1$. Attempted below (H8) — same magnitude wall.

## H8 — $\theta'$ angle-correction to $\sin^4(\pi/5)$: found, then retracted

Noted that $\alpha^{-1}_{\rm geo}=\varphi^6/\sin^4(\pi/5)$ has the same
algebraic form as $1/S(\pi/5)$, where $S(\theta)=\sin^4\theta/\varphi^6$
is the established condensate suppression function (`P1:eq:suppression`,
used with a genuine physical meaning in `P1:thm:amplitude_ratio`, the
Weinberg-angle derivation, where $\theta$ is the **polar angle on the
Hopf-fibration base $S^2$** and the physically-justified value is
$\theta=\pi/2$ — the equatorial vacuum locus where the $Q=2$ soliton
actually sits).

Solved for the angle $\theta'$ at which $\varphi^6/\sin^4(\theta')$
exactly hits the golden angle: $\theta'=36.944°$, vs. $\pi/5=36°$ — only
a $0.944°$ ($2.62\%$) shift. Matched this shift to
$1/38=1/(2(2R_0^2+1))$ at $R_0=3$ (an already-established combination
from the virial-factorization chain) to $0.35\%$, and checked the
resulting angle directly: $\theta'=\frac{\pi}{5}\cdot\frac{39}{38}$
gives $\alpha^{-1}=137.465$ vs. golden angle $137.508$ — a $0.031\%$
match, using only already-established $R_0=3$ quantities.

**Retracted.** This assumed $\theta=\pi/5$ is a genuine Hopf-polar
*location* (like the established $\theta=\pi/2$ vacuum angle) subject
to a thick-torus/$R_0$-dependent geometric correction. But $\pi/5$ is
much more plausibly **the pentagon angle** — a fixed, purely
group-theoretic quantity tied to $k+2=5$ (the same $5$ underlying
$\dim_q(1/2)=\sin(2\pi/5)/\sin(\pi/5)=2\cos(\pi/5)=\varphi$,
`P1:thm:dimq`, already proved), not a dynamical location that could
shift with $R_0$. A fixed group-theoretic angle has no principled
reason to receive an $R_0$-dependent correction — $R_0$ is a classical
profile parameter, unconnected to $k$ or the pentagon geometry. Given
that, the $0.031\%$ match is most likely another coincidental numerical
fit (same species as the numerology killed earlier in this
investigation — $\ln(72)=k|2T|$, $\beta^*\times400$ — just tighter
because it came from fitting a continuous parameter rather than
choosing from a short list of integers), not a genuine derived
correction. **Downgraded from "leading candidate" back to "unproven,
likely coincidental."** The real question — is $\sin^4(\pi/5)$ in
$\alpha^{-1}_{\rm geo}$ a WZW-modular quantity (pentagon angle, fixed,
$k$-determined) or a Hopf-polar quantity (location, $R_0$-dependent) —
now looks like it should resolve toward "WZW-modular," which is a
different derivation path than anything tried above (needs CFT/Verlinde
refinement, not profile-geometry refinement). Next: examine the
WZW-modular reading directly.

## 2026-08-31 — Cross-check: is OP4's scale ratio the recombination-CDM "non-ultralight orientation scale"? NO.
Asked whether the director-sector mass that recombination CDM would need (see
`fabric_cmb_hypothesis_and_plan.md`: a non-baryonic, unscreened, m>>H_rec, adiabatic
orientation config) relates to Paper III's μ_UV/μ_IR. Answer: no derived connection, and
two reasons not to force one:
  1. PREMISE STALE: the μ_UV/μ_IR *running* framing of Δ₁ is RETRACTED in current
     `main_paper3.tex` (`P3:prop:cascade` item 1: "the retracted running framing"). Δ₁ is
     the exact algebraic difference (2136φ−3392)/5 between two independently-exact icosahedral
     φ-quantities (φ⁶/sin⁴(π/5)=150.33, 360/φ²=137.508), NOT an energy-scale ratio. OP4 is now
     "does Δ₁'s value / the endpoint near-coincidence have a deeper origin." (`beta_flow_derivation.md`
     l.4 still uses the old "scale ratio" language — stale framing.)
  2. TYPE/SECTOR MISMATCH: Δ₁ is dimensionless, in the EM/charge sector (α, Q_H=2 WZW). The CDM
     scale is a dimensionful MASS in the director/DM sector, non-ultralight. A number ≠ a mass
     without a dimensionful anchor + a reason the DM sector gains a 2nd scale. Equating them =
     the coincidence-restated trap (CLAUDE.md §3).
Only honest thread: Δ₁'s dielectric candidate route uses the feedback coupling β
(`delta1-dielectric-program`), and the DM sector's screening also uses β — a shared ACTOR, not a
bridge. DM ultralightness (m_ξ~H₀) is welded to the DE scale, not α's UV endpoint. => the two open
problems are independent; do not conflate.
