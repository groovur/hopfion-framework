# Claude Code Instructions — Hopfion Framework Repository

## Repository purpose
This repository extracts, organises, and cross-checks all formal results
(theorems, propositions, corollaries, conjectures, constructions, remarks, definitions, boxed results, and any labeled equations.)
from the Density-Feedback Faddeev–Niemi Hopfion paper series (Papers I–XIX).
It is the source-of-truth for a downstream Q&A system.

Single physical input: T_CMB = 2.7255 K. All outputs derive from this.

## Ground truth files
- `paper_registry.yaml` — all papers, their status, and the label index
- `master_table.tex`    — all numbered outputs with status/residual/source
- `open_problems.tex`   — open problems OP1–OP8
- `preamble.tex`        — shared macros (\status, \source, \residual, \depends)

## Paper status
- **published**: Papers I–XVIII. Treat content as immutable.
  Do not edit the paper .tex files. Only edit repo subfolder files.
- **draft**: Paper XIX. Content may change.
  Flag any result extracted from a draft with \status{draft}.

## Sector → subfolder mapping
When a result is extracted, file it here:

| Sector | Subfolder |
|--------|-----------|
| Pentagon, suppression, α⁻¹ | foundations/ |
| Lepton masses, Weinberg, Yukawa | electroweak/ |
| Colour, α_s, quark masses, strong CP | strong/ |
| Q_H=3 confinement, Jones polynomial, jets | strong/ |
| Gravity, cosmology, dark sector | gravity/ |
| Hilbert space, Schrödinger, Bell, chameleon | quantum/ |
| Neutrino sector, Q_H=1, Majorana | quantum/ |
| φ-spiral, chemistry, cooperativity | chemistry/ |
| Profile normalisation, virial | profile/ |
| Higgs | scalar/ |

If a result doesn't fit any existing subfolder, create a new one and
add it to this table before committing.

---

## Label naming convention

### Repo sector files (foundations/, electroweak/, strong/, gravity/,
quantum/, chemistry/, profile/, scalar/)
Every `\label{}` in a repo subfolder file — theorem-like environments
*and* labeled equations — MUST be prefixed with the source paper id:

    \label{P5:prop:strongCP}
    \label{P14:eq:yukawa_3D}

Do NOT reuse the bare label from main_paperN.tex directly; prepend
`PN:`. This is what prevents collisions when a sector file aggregates
results from several papers (or, as found in strong/strong_cp.tex,
two labels landing in the same file unprefixed).

If the label extracted from the source paper already carries a
`PN:` prefix (e.g. it was already written as `PN:prop:foo` in
main_paperN.tex, or it's being re-filed from a repo file that was
already migrated), do NOT prepend a second prefix. Check for an
existing `P\d+:` prefix before prepending; only add one if absent.

paper_registry.yaml `labels:` entries and master_table.tex `Source`
cells must match the prefixed label exactly.

For main_paper_foundational.tex, use prefix 'F'

For main_reader_guide.tex, use prefix 'RG'

Papers themselves (papers/main_paperN.tex) are read-only — this
convention applies only to repo subfolder files, never to the
original paper source.

**Unregistered papers.** Some sector files contain content sourced
from a paper that has not yet been added to paper_registry.yaml
(check `papers/` for a `main_paperN.tex` with no matching entry).
Always prefix the label with that paper's number regardless — you
know which paper a result came from (via its \source{}, a nearby
header comment, or by checking main_paperN.tex directly), so there
is no reason to leave it bare. Registration in paper_registry.yaml
(via NEW PAPER) is a separate, later step; it does not gate labeling.

### Open problems (open_problems.tex)
Item identifiers are paper-prefixed with a per-paper-local counter,
not a global sequential one:

    \item \textbf{[P16:OP1] ...}
    \item \textbf{[P18:OP3] ...}

The counter restarts at 1 for each paper's first open problem (found
during that paper's SCAN PAPER run), not across the whole document.

### master_table.tex `#` column — UNCHANGED
This is the one place a global sequential integer is kept, because
published papers cite results by this number. Once assigned, a `#`
is frozen forever — never reused, never renumbered — even if rows
are reordered for readability.

### master_table.tex residual progressions
When a quantity's accuracy is refined across successive papers or
successive results within one paper (e.g. a leading-order formula,
then a corrected version, then a further thick-torus/higher-order
correction), do NOT collapse this into a single row showing only the
"best" (or worst, or first-found) number. Split into lettered
sub-rows sharing the same base `#`, one per distinct derivation,
in refinement order:

    8a & v_EW (leading order)          & ... & 1.45\%  & P4: cor:vew \\
    8b & v_EW (tower construction)     & ... & 1.14\%  & P6: thm:vew_tower \\
    8c & v_EW (T-matrix corrected)     & ... & 0.69\%  & P4: cor:vew_corrected \\
    8d & v_EW (thick-torus corrected)  & ... & 0.015\% & P14: cor:yukawa_3D \\
    8e & v_EW (thick-torus, final)     & ... & 0.012\% & P14: cor:yukawa_second \\

This makes "what was the accuracy at Paper N?" answerable directly
from the table, and prevents a later paper's number from being
silently mis-cited under an earlier paper's label (this happened:
row 8 showed `0.69\%` cited to `P4: cor:vew`, but `cor:vew` itself
states `1.45\%` — the `0.69\%` figure actually belongs to the
distinct, later corollary `cor:vew_corrected` in the same paper).
Always verify a residual number against the actual `\residual{}` /
corollary text in the source paper before trusting an existing table
entry — do not assume a pre-existing row is correct.

If a later correction is documented elsewhere in the table under a
different quantity's section (e.g. a "Thick-torus corrections"
section), it's fine for both to exist — one as the terminal row in
this quantity's own progression, one in its thematic section — as
long as neither contradicts the other.

---

## Operations

### SCAN PAPER
Trigger: "scan paperN" or "scan papers/main_paperN.tex"

Steps:
1. Read papers/main_paperN.tex
2. Extract every environment matching:
   theorem, lemma, proposition, condition, corollary, conjecture, construction, definition,
   remark (if it contains a named result), and any \boxed{} or
   equation with a \label. You may use a subagent for this.
3. For each extracted item:
   a. Identify its \label
   b. Check paper_registry.yaml labels section — is it already filed?
   c. If yes: verify the repo file content matches. Report discrepancy if not.
   d. If no: determine the correct subfolder from the sector mapping
   e. Verify in the preamble.tex that the macros/theorem environment/usepackage for the results are defined
   f. If no: add the missing LaTeX macro/theorem environment/usepackage defined in the paper for the macro in the result to preamble.tex in the appropriate location
   g. If the target .tex file exists: append the result in the correct section
   h. If the target .tex file does not exist: create it with the standard
      preamble (see TEMPLATE below), then add the result
   i. \ref DOENSN'T NEED to be in the same file. The tex sector files are collections of VERBATIM paper theorm/lemmas/etc.
4. Add the label to paper_registry.yaml under the correct paper entry
5. If the paper is numbered, ensure it's paper_registry.yaml entry is positioned correctly in the ascending sequentially order.
6. If the doi for the paper is missing in paper_registry.yaml, try to find the DOI of the paper being scanned with grep {PaperN} from the bibliography section of the existing main_paperN.tex papers in the /papers folder, whether or not they are in paper_registry.yaml
7. Verify the name for the paper in paper_registry.yaml, use the name of the paper in the current paper. Remove any new lines in the title and space it appropriately.
8. If the result is a new numbered output: 
   a. Add a row to master_table.tex
   b. Ensure that the extracted content of the label to be put in the repository is verbatim to the paper
9. If the result no longer exists in the corresponding paper, remove it
10. If the result exists, run CHECK CONSISTENCY on the label
11. If the result closes an open problem: update open_problems.tex
12. Report: list of extracted labels, any discrepancies, any new files created
13. If there are any discrepancies:
    a. Attempt to follow the existing conventions of the repository
    b. If unsuccessful, interactively resolve any discrepancies
14. Once discrepancies are resolved (with or without interation) prompt for confirmation and on confirmation, continue.
15. git add all modified files
16. prompt for confirmation to run "git commit -m "scan: extract [N results] from paperN ([status])"

### ADD RESULT
Trigger: "add this [theorem/proposition/...] to [subfolder/file.tex]"

Steps:
1. Read the target .tex file
2. Add the result with \label, \status{}, \source{}, \residual{} (if applicable)
3. Add the label to paper_registry.yaml
4. Update master_table.tex if it is a new numbered output
5. Update open_problems.tex if it closes an open problem
6. git add + git commit -m "add: [label] to [file]"
7. git tag if this is a paper milestone (ask if unsure)

### CHECK CONSISTENCY
Trigger: "check consistency of [label]"

Steps:
1. Read paper_registry.yaml to find all files referencing [label]
2. Read each file
3. Compare the statement of the result across all occurrences
4. If the result exists, but has the same label after ':' ensure that if the type was promoted, it is updated correctly, eg conjecture to proposition to theorem, etc.
5. If the result exists, and it is the same label, ensure it's content is correct against the main_paperN.tex
6. If the result exists, and it is the a different label than in the paper, update the label.
7. If the result is better and superseded by a later paper, identify it as such. Keep each paper's results unique and reference the paper where it is superseded. You may need to check against existing results.
8. Make sure that if the staement contains an equation reference, the equation and it's label is also defined somewhere.
9. Report: identical / discrepant (show diff)
10. If discrepant: ask which version is canonical before making any change

### NEW PAPER
Trigger: "new paper paperN is [published/draft]"

Steps:
1. Check if the entry in paper_registry.yaml exists.
2. If no: Add entry to paper_registry.yaml with status and doi (null if draft), if the paper is numbered, add it in the correct ascending sequential order.
3. Run SCAN PAPER on the new file
4. If published: set status: published for all labels from that paper
5. If draft: set status: draft for all labels from that paper
6. If the doi for the paper is missing in paper_registry.yaml, try to find the DOI of the paper being scanned from the bibliography section from the other existing tex papers in the /papers folder. hint: you can grep for {PaperN} 
7. Verify the name for the paper in paper_registry.yaml, use the name of the paper in the current paper. Remove any new lines in the title and space it appropriately.
8. Add any open problems to open_problems.tex
8. Confirm with user, and run "git commit -m "registry: add paperN ([status])"


### PAPER STATUS UPDATE
Trigger: "paperN is now published" or "update doi for paperN"

Steps:
1. Update paper_registry.yaml: set status: published, add doi
2. If scan paper has not run, run SCAN PAPER
3. Update all repo files that have \status{draft} from that paper → \status{published}
4. If the doi for the paper is missing in paper_registry.yaml, try to find the DOI of the paper being scanned from the bibliography section of existing tex papers in the /papers folder
5. Verify the name for the paper in paper_registry.yaml, use the name of the paper in the current paper. Remove any new lines in the title and space it appropriately.
6. git commit -m "registry: paperN status → published, doi: [doi]"

### GENERATE PAPER
Trigger: "generate a paper on topic X"

Steps:
1. Identify all relevant labels for topic X from paper_registry.yaml
2. Read the corresponding repo files
3. Assemble a paper using \input{} — no re-derivation, no duplication
4. Reference results by label, not by restating them
5. Output to papers/ directory

### CHECK MISSING
Trigger: "check what is missing from [subfolder/file.tex]" or
         "check paperN for missing results"

Steps:
1. If checking a paper: run SCAN PAPER (step 2 only) and compare
   found labels against paper_registry.yaml — report unregistered ones
2. If checking a repo file: list all labels in the file, check each
   exists in paper_registry.yaml, report orphans

### QUERY FRAMEWORK
Trigger: any question about what the framework says, predicts, derives,
or relates to (e.g. "how does the hopfion relate to X", "what does the
framework say about Y", "what's the status of Z") — i.e. read-only
questions, not a request to extract, add, or modify results.

This is a read query, not a SCAN PAPER run: never edit repo files or
papers/ while answering one. If the answer reveals a genuine gap
(missing label, stale residual, orphaned reference), report it and ask
before fixing — do not fold a fix into a query response.

Steps:
1. **Repo first.** Identify the sector(s) the question falls under
   (see the Sector → subfolder mapping) and read the corresponding
   subfolder file(s). The repo files are pre-digested and
   cross-checked — prefer them over the papers whenever they cover
   the question.
2. **Use the ground truth files as the index**, not just grep:
   - `paper_registry.yaml` — which paper(s) and labels are relevant;
     also the fastest way to check a paper's status (published/draft)
     before trusting a result at face value.
   - `master_table.tex` — the quantitative headline: value, residual,
     status glyph (\checkmark/$\circ$/P), and source label, in one
     row. Good for "what's the accuracy of X" questions.
   - `open_problems.tex` — whether a result is settled or still open;
     check this before presenting a conjecture as if it were closed.
   - `preamble.tex` — resolve unfamiliar macros encountered while
     reading a sector file.
3. **Respect `\status{}`.** Present `\status{draft}` results as
   provisional (flag that they come from an unpublished paper) and
   `\status{open}`/`\circ` conjectures as unproved, even if the repo
   text reads confidently. Never launder a draft or conjectural result
   into a stated fact.
4. **Fall back to `papers/main_paperN.tex` (or main_paper_foundational.tex
   / main_reader_guide.tex) only when the repo is insufficient**:
   the result isn't filed yet (check paper_registry.yaml — it may
   predate the last SCAN PAPER run), the repo file's terse extracted
   form lacks context the question needs (e.g. the proof itself, or
   surrounding discussion), or the question is about a paper's
   narrative/framing rather than an extracted formal result. When
   falling back, say so explicitly ("not yet filed in the repo, so
   drawn directly from Paper N") rather than presenting it as if it
   came from the indexed repo.
5. **Cite what you used.** Reference results by their `PN:label` (and
   file path), the way the repo cross-references itself. If multiple
   papers/labels bear on the answer (e.g. a residual progression across
   several `#`-numbered sub-rows), mention the progression rather than
   only the latest number.
6. If the question spans a topic with no subfolder yet, or the
   answer genuinely isn't in the repo or papers at all, say so plainly
   instead of speculating past the framework's stated results.

---

## Standard repo file template

When creating a new .tex file in a subfolder:

\`\`\`latex
% [SUBFOLDER]/[filename].tex
% Part of the Hopfion Framework Repository
% Source: [Paper N] — [Paper Title]
% Status: [published|draft]
% Last updated: [date]
%
\input{../preamble}

\section*{[Section Title]}
\label{sec:[filename]}

% Results extracted from [Paper N]
% Scan date: [date]
\`\`\`

---

## Preamble macros (defined in preamble.tex)

| Macro | Meaning |
|-------|---------|
| \status{proved} | Formally proved in a published paper |
| \status{draft}  | From a draft paper — may change |
| \status{open}   | Open problem or conjecture |
| \status{prediction} | Untested experimental prediction |
| \source{PaperN, label} | Where the proof lives |
| \residual{0.013\%} | Accuracy vs experiment |
| \depends{file1, file2} | Upstream dependencies |

---

## What NOT to do
- Do not edit files in papers/ (treat them as read-only)
- Do not re-derive results — \input or \ref only
- Do not add a result to the repo without updating paper_registry.yaml
- Do not change \status{proved} to anything else without explicit instruction
- Do not resolve a consistency discrepancy without asking which is canonical
