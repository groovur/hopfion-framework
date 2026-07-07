# Density-Feedback Faddeev–Niemi Hopfion Framework
## Mathematical Knowledge Repository

Single input: **T_CMB = 2.7255 K**. All results below derive from this alone,
plus geometry and group theory — no fitted parameters.

This repository extracts, organises, and cross-checks every formal result
(theorems, propositions, corollaries, conjectures, constructions, remarks,
definitions, and labeled equations) from the Density-Feedback Hopfion paper
series (Papers I–XIX) into a single, internally consistent, queryable index.
It is the source of truth for a Q&A system, and for anyone who
wants to check a claim, a residual, or an open problem without reading
nineteen papers.

---

### Repository Structure

```
hopfion-framework/
├── paper_registry.yaml       # Every paper: status, DOI, and its result labels
├── master_table.tex          # All numbered outputs: value, residual, status, source
├── open_problems.tex         # Current open problems, per-paper (P#:OP#)
├── preamble.tex              # Shared macros, theorem environments, status commands
│
├── foundations/
│   ├── pentagon.tex          # Pentagon Theorem (k+2=5), 4 parts, 3 routes
│   ├── suppression.tex       # S(θ,ρ)=sin⁴θ/[φ⁶(1+β*ρ)], saddle point, WZW
│   ├── alpha.tex             # α⁻¹ three- and four-term formulas and all proofs
│   └── wzw_group.tex         # McKay correspondence, WZW group-theory backbone
│
├── electroweak/
│   ├── lepton_masses.tex     # T-phases, φ-tower, m_e chain
│   ├── weinberg.tex          # sin²θ_W, M_W/M_Z, v_EW
│   └── yukawa.tex            # y_e, MSbar→pole, Koide, thick-torus corrections
│
├── strong/
│   ├── colour.tex            # 2T⊂2I, McKay quiver, N_c=3, confinement
│   ├── qcd_coupling.tex      # d_eff=43/14, α_s, Λ_QCD, chameleon decoupling
│   ├── quark_masses.tex      # Quark Koide, CKM angles, CP phase, E8 assignment
│   ├── strong_cp.tex         # θ_QCD=0 from Hopf fibre = Peccei–Quinn
│   ├── topology_qh3.tex      # Q_H=3 sector topology, numerical construction
│   └── confinement_topology.tex # Jones-polynomial confinement, jets, quark-mass dynamics
│
├── gravity/
│   ├── newton.tex            # G_N, Einstein equations, Cassini bound, Kerr background
│   ├── cosmology.tex         # Inflation, n_s, r, P_s, N_e, reheating
│   ├── dark_sector.tex       # Λ_obs, w_a, c_s=1/φ, BAO shift, Vainshtein
│   └── wep_instanton.tex     # Weak equivalence principle, instanton suppression
│
├── quantum/
│   ├── hilbert_space.tex     # H=L²(C_Q^red,dμ), Born rule, φ²ⁿ spectrum
│   ├── schrodinger.tex       # Schrödinger equation derived from the condensate
│   ├── bell.tex              # Born rule, Malus, three CHSH regimes, Hurwitz algebras
│   ├── chameleon.tex         # Screening hierarchy: cosmic/atomic/QCD
│   ├── spin_statistics.tex   # Spin-1/2, spin tower, Pauli exclusion from Q-topology
│   └── neutrino.tex          # Q_H=1 sector, PMNS, Majorana mass, solar angle
│
├── chemistry/
│   ├── phi_spiral.tex        # φ-spiral, Z_eff(O), Rydberg tower level
│   ├── icosahedrite.tex      # Al₆₃Cu₂₄Fe₁₃: Al/Cu=φ², Fe=F₇
│   └── cooperativity.tex     # n_H ≤ k=3 from WZW fusion 3/2⊗3/2=0
│
├── profile/
│   ├── normalisation.tex     # φ⁹√(J_aJ₄)=10, x*=R₀/C*, one-input chain
│   ├── virial.tex            # Virial condition, saddle point, β*=0.452
│   └── bps.tex               # BPS structure, fixed-point theorem
│
└── scalar/
    └── higgs.tex             # Higgs mass and quartic coupling from the WZW condensate
```

New sectors are added to this table (and to `CLAUDE-MAINTAINER.md`'s sector map) the
moment a result doesn't fit an existing subfolder — see `CLAUDE-MAINTAINER.md` for the
filing rule.

---

### Source Papers

| Key | Title | Status | DOI |
|-----|-------|--------|-----|
| Foundational | Pentagon Theorem | published | [10.5281/zenodo.20173651](https://doi.org/10.5281/zenodo.20173651) |
| Reader's guide | — | published | [10.5281/zenodo.20173714](https://doi.org/10.5281/zenodo.20173714) |
| Paper I | The Density-Feedback Faddeev–Niemi Hopfion | published | [10.5281/zenodo.19342027](https://doi.org/10.5281/zenodo.19342027) |
| Paper II | BPS Structure and Fixed-Point Theorem | published | [10.5281/zenodo.19363491](https://doi.org/10.5281/zenodo.19363491) |
| Paper III | Two Constants from One Knot | published | [10.5281/zenodo.19478629](https://doi.org/10.5281/zenodo.19478629) |
| Paper IV | The Hopf Spoke, WZW Fermion Mass Renormalization | published | [10.5281/zenodo.19504729](https://doi.org/10.5281/zenodo.19504729) |
| Paper V | Colour Confinement, the QCD Scale | published | [10.5281/zenodo.19638857](https://doi.org/10.5281/zenodo.19638857) |
| Paper VI | The Golden Tower: Fermion Scale Hierarchy | published | [10.5281/zenodo.19646945](https://doi.org/10.5281/zenodo.19646945) |
| Paper VII | Emergent Gravity, Inflation, Dark Energy, WEP | published | [10.5281/zenodo.19646953](https://doi.org/10.5281/zenodo.19646953) |
| Paper VIII | Geometry-Dependent Bell Violations | published | [10.5281/zenodo.18737070](https://doi.org/10.5281/zenodo.18737070) |
| Paper IX | The Born Rule, the Tsirelson Bound, Hurwitz Algebras | published | [10.5281/zenodo.20075227](https://doi.org/10.5281/zenodo.20075227) |
| Paper X | Quantisation of the Hopfion Condensate | published | [10.5281/zenodo.20075279](https://doi.org/10.5281/zenodo.20075279) |
| Paper XI | The Golden-Spiral as a Universal Mass Hierarchy Map | published | [10.5281/zenodo.20173763](https://doi.org/10.5281/zenodo.20173763) |
| Paper XII | The Profile Normalisation Conjecture | published | [10.5281/zenodo.20210460](https://doi.org/10.5281/zenodo.20210460) |
| Paper XIII | Atomic Quantum Mechanics from the Hopfion Condensate | published | [10.5281/zenodo.20471482](https://doi.org/10.5281/zenodo.20471482) |
| Paper XIV | The Thick-Torus Profile Correction | published | [10.5281/zenodo.20479790](https://doi.org/10.5281/zenodo.20479790) |
| Paper XV | The Q_H=3 Sector of the Density-Feedback Hopfion | published | [10.5281/zenodo.20691001](https://doi.org/10.5281/zenodo.20691001) |
| Paper XVI | Quark Generation Masses: Survey of Ruled-Out Mechanisms | published | [10.5281/zenodo.21012650](https://doi.org/10.5281/zenodo.21012650) |
| Paper XVII | The Q_H=1 Sector: Topology, Group Structure, Colour Exclusion | published | [10.5281/zenodo.21013412](https://doi.org/10.5281/zenodo.21013412) |
| Paper XVIII | Preimage Topology and Confinement | **draft** | [10.5281/zenodo.21047750](https://doi.org/10.5281/zenodo.21047750) |
| Paper XIX | Dynamical Origin of Quark Masses | **draft** | [10.5281/zenodo.21225452](https://doi.org/10.5281/zenodo.21225452) |

`paper_registry.yaml` is the authoritative source for this table — it also
lists every extracted result label per paper. If this README and the
registry ever disagree, trust the registry.

---

### Viewing numerical field snapshots (Paper XV onward)

From Paper XV on, the numerical construction/gradient-flow scripts in each
`papers/src_paperN/` folder save condensate field snapshots as `.npy`
files. These can be inspected directly in a browser with
`papers/src_paper16/hopfion_viewer.html`. Drag a `.npy` file onto the page
(no server or install needed; it renders the field as 8000 - 3×10⁵ particles
coloured by $\rho_{J_4}$ intensity). There is a directional vector option to
visualize the field flow. The viewer lives in `src_paper16` but
is the shared tool referenced by later papers' scripts (XVII, XVIII, XIX)
too, there's only the one copy.

Some scripts (e.g. `qh_dynamic_integrator_v2.py`) save `.npz` archives
instead of bare `.npy` — the viewer only accepts `.npy` directly, so pull
the field array `.npy`s out first with an archiver utility or:

```python
import numpy as np
np.save('snapshot.npy', np.load('snapshot.npz')['n'])
```

**Static, self-contained visualisations** (no data file needed — just open
in a browser, each is an interactive orbit-controllable 3D scene with its
own sliders/buttons):

| File | Shows |
|------|-------|
| `papers/src_paper14/hopfion_torus.html` | The $Q_H=2$ Hopfion torus (topological soliton of the condensate, $R_0=3$, $C^*=3.4318$) |
| `papers/src_paper15/hopfion_trefoil.html` | The $Q_H=3$ trefoil (quark sector — three quark colours from one knot, $T_{2,3}$, $R_0=3$, $r_0=0.874$) |
| `papers/src_paper11/phi_spiral_masses.html` | The $\varphi$-spiral fermion mass hierarchy ("The Golden Tower") |

---

### Using Claude Code to interact with this repository

This repo carries three Claude-facing instruction files, each for a
different mode of working with it. Point Claude Code at whichever matches
what you're trying to do — they're not interchangeable. Point your session
to one of these (or copy one of them to `CLAUDE.md` in your own clone so it
loads automatically).

- **`CLAUDE.md`** — for *asking questions* ("does this framework predict
  X?", "what's the accuracy of Y?", "is Z proved or still open?") without
  editing anything. If you're cloning or browsing this repository and just
  want a read-only Q&A assistant. It tells Claude
  to check `master_table.tex`, `paper_registry.yaml`, and `open_problems.tex`
  before answering, to respect each result's `\status{}` tag rather than
  presenting a draft or conjectural result as settled fact, and to flag any
  of its own speculation about connections the framework doesn't actually
  derive. This is the one that's active by default in this directory (Claude
  Code auto-loads a root `CLAUDE.md`).

- **`CLAUDE-MAINTAINER.md`** — for *maintaining* the repo: extracting results
  from a paper into the sector files (`SCAN PAPER`), adding a single new result,
  checking a label's content is consistent across every file it appears in,
  registering a new paper, or generating a derivative paper from existing
  labeled results. It's meant to be used when you're actively changing the repo,
  not just querying it.

- **`CLAUDE-RESEARCH-ASSISTANT.md`** — working principles for open-ended
  *research* sessions (chasing down a bug in a derivation, exploring whether
  a numeric coincidence means anything, running numerical experiments),
  distilled from an actual multi-session investigation that found and fixed
  real errors propagated across papers. It's not a set of commands like
  the other two files — it's guidance for Claude on how to stay honest 
  mid-investigation (verify computationally rather than trust hand algebra,
  distinguish a restated coincidence from independent evidence, track what's
  established vs. conjectural at every step). Worth reading before a long
  back-and-forth research session.

---

### Git Workflow

```bash
# Adding a new result after a derivation session:
git checkout -b new-result-name
# Edit the relevant .tex file, master_table.tex, open_problems.tex, paper_registry.yaml
git add relevant/file.tex master_table.tex open_problems.tex paper_registry.yaml
git commit -m "Brief description of what was derived and from what"
git checkout main && git merge new-result-name

# Generating a paper from repository sources:
# In your paper's .tex file:
\input{../hopfion-framework/foundations/pentagon}
\input{../hopfion-framework/quantum/schrodinger}
# No re-derivation. No duplication. Guaranteed consistency.
```

---

### Metadata Commands (defined in `preamble.tex`)

| Command | Usage |
|---------|-------|
| `\status{proved}` | Formally proved in a published paper |
| `\status{draft}` | From a draft paper — may change |
| `\status{open}` | Open problem or conjecture |
| `\status{prediction}` | Untested experimental prediction |
| `\source{PaperN, label}` | Where the proof lives |
| `\residual{0.013\%}` | Numerical accuracy vs. experiment |
| `\depends{file1, file2}` | Dependencies on other repository files |

---

### Current Status

- **87 numbered outputs** in `master_table.tex` (some split into lettered
  sub-rows to preserve the refinement history of a quantity across papers)
- **27 open problems** tracked in `open_problems.tex`, filed per-paper
  (e.g. `P16:OP1`, `P19:OP6`)
- **21 registry entries** (foundational paper + reader's guide + Papers
  I–XIX), all published except Paper XIX (draft)
- **580 individually labeled results** cross-referenced in
  `paper_registry.yaml`
- Label convention: every result in a sector file is prefixed by its source
  paper (`P5:thm:...`, `P19:prop:...`) to prevent collisions when a sector
  file aggregates results from multiple papers
