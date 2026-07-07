# Claude Code Instructions — Hopfion Framework (Public / Read-Only Mode)

## What this repository is

This repository is a structured index of every formal result (theorems,
propositions, corollaries, conjectures, constructions, remarks, definitions,
and labeled equations) from the Density-Feedback Faddeev–Niemi Hopfion paper
series — a physics framework that derives Standard Model and cosmological
parameters from a single measured input, the CMB temperature
`T_CMB = 2.7255 K`, plus geometry and group theory, with no fitted
parameters.

If you are a Claude Code session that has just been pointed at this
repository, this file tells you how to help a user explore it accurately.
See `README.md` first for the directory layout and paper list.

**This file governs read-only, exploratory use.** Treat the repository's
`.tex`, `.yaml` files as reference material to read, not to edit, unless the
user explicitly asks you to change something — and if they do, confirm what
they want changed before touching any file.

## The one rule that matters most

**Never present a result more confidently than its own status tag allows.**
Every extracted result in this repo carries a `\status{}` tag:

| Tag | Meaning |
|-----|---------|
| `\status{published}` | Formally proved in a **published** paper |
| `\status{draft}` | From a **draft, unpublished** paper — may still change |
| `\status{open}` | An open problem or conjecture, not resolved |
| `\status{prediction}` | An untested experimental prediction |

The text of a result is often written confidently (it's copied near-verbatim
from the source paper), but that confidence is about the paper's internal
argument, not about the result's standing. Check the tag — and the paper's
status in `paper_registry.yaml` — before repeating a claim as settled fact.
As of this writing, Paper XIX is a draft; anything sourced only to Paper XIX
should be flagged as such.

## How to answer a question about the framework

1. **Identify the sector.** Use the table below to map the topic to a
   subfolder, then read the relevant `.tex` file(s) directly
   — they are pre-digested extracts, not the raw papers, and are the
   fastest and most reliable source. Extract the macro from preamble.tex
   and replace the macro in the response with the expanded macro.

   | Topic | Subfolder |
   |-------|-----------|
   | Pentagon, suppression, α⁻¹ | `foundations/` |
   | Lepton masses, Weinberg angle, Yukawa couplings | `electroweak/` |
   | Colour, α_s, quark masses, strong CP, Q_H=3 confinement | `strong/` |
   | Gravity, cosmology, dark sector | `gravity/` |
   | Hilbert space, Schrödinger eq., Bell tests, chameleon screening, neutrinos, Q_H=1 | `quantum/` |
   | φ-spiral, chemistry, cooperativity | `chemistry/` |
   | Profile normalisation, virial theorem | `profile/` |
   | Higgs sector | `scalar/` |

2. **Use the ground-truth index files, not just grep:**
   - `paper_registry.yaml` — which paper a result comes from, that paper's
     publication status, and the DOI. This is the fastest way to check
     whether something is trustworthy-as-published or still a draft.
   - `master_table.tex` — the one-line summary of every numbered result:
     value, residual against experiment, status glyph
     (✓ = proved, ○ = conjectured, P = untested prediction), and source
     label. Best single place to answer "what's the accuracy of X?" —
     including its lettered sub-rows (e.g. `8a`–`8e`, `48d`–`48e`), which
     record how a quantity's derivation was refined across multiple papers.
     Don't just quote the first row you find for a quantity — check
     whether later sub-rows supersede it.
   - `open_problems.tex` — whether something is actually settled. Numbered
     per-paper (`P16:OP1`, `P19:OP6`, etc.). Check this before presenting a
     conjecture as resolved.
   - `preamble.tex` — resolves macros you don't recognize while reading a
     sector file (e.g. `\vp` = φ, `\QH` = $Q_H$, `\twoI`/`\twoT` = the binary
     icosahedral/tetrahedral groups).

3. **Fall back to the original papers only when the repo falls short.** The
   papers live in `papers/main_paperN.tex` (plus
   `main_paper_foundational.tex` and `main_reader_guide.tex`). Use them when:
   a result hasn't been filed in the repo yet (check `paper_registry.yaml`
   first — the last extraction pass may predate a paper's addition), the
   repo's terse extracted form lacks context the question needs (e.g. the
   actual proof, or surrounding discussion), or the question is about a
   paper's narrative rather than a specific formal result. When you do fall
   back, say so explicitly, e.g. "not yet indexed in the repo — this is
   drawn directly from Paper N."

4. **Cite what you used**, the way the repo cross-references itself: by
   `PN:label` and file path (e.g. `P19:thm:r0_derivation`,
   `strong/confinement_topology.tex`). If a quantity has a residual
   progression across several master-table sub-rows, mention the
   progression, not just the most recent number.

5. **If it genuinely isn't in the repo or the papers, say so.** This
   framework is explicit about its own gaps — `main_reader_guide.tex` has
   entire sections titled "What the framework does not predict" and "The
   genuine tensions and open questions"; read those before assuming a
   result exists.

## Speculating about missing connections

You are free to speculate about how the framework *might* connect to
established physics or open problems it doesn't currently address — e.g.
"this quark-mass mechanism resembles X in the SM literature," or "the same
$k+2=5$ pattern could plausibly extend to Y, though the repo doesn't derive
that." This is often exactly what a curious reader wants.

But you must **clearly and unmistakably flag it as your own speculation**,
separate from what the repo actually derives — never blend a guess into a
sentence that reads like an extracted result. Say "the repo does not derive
this, but one could speculate that..." or similar, every time, not just
once at the top of a long answer. If you're not sure whether something is
in the repo or you're extrapolating, check first (see steps 2–3 above)
before presenting it as established.

## Do not, in this mode

- Do not edit `.tex` files, `master_table.tex`, `open_problems.tex`, or
  `paper_registry.yaml` unless the user has explicitly asked you to change
  something in the repository — confirm what they want changed first.
- Do not treat a `\status{draft}` or `\status{open}` result as equivalent to
  a `\status{proved}` one, even if asked to summarize "what the framework
  predicts" in aggregate.
- Do not silently patch an inconsistency you notice (a stale cross-reference,
  a residual that doesn't match its source, a label that resolves to the
  wrong file). Report it to the user and ask before changing anything.
