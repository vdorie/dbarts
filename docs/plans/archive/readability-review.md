# readability-review

agent: fable orchestrates; sonnet readers, opus for engine-header calls
rng: neutral (comments, docs, naming only)
budget: findings list first; fix diffs are comment/doc-only

## Goal

Every file the bartcore branch touched reads for a maintainer who never
saw an agent conversation. This is the retroactive pass; the standing
per-review check lives in docs/plans/README.md.

## Context

- Scope list: git diff main...bartcore --name-only (236 files; R/,
  src/bartcore/, src/R_interface*, inst/tinytest/, docs/design/, man/).
- Known hotspot: docs/design/core-generalization.md is a chronological
  landing-note accretion (768 lines) - correct as history, poor as a
  current-state reference. public-surface.md has stray artifacts (e.g.
  an orphaned sentence fragment near its section 6 RESOLVED note).
- The house style is already written down (CLAUDE.local.md: Doxygen
  LLVM style, ASCII, brief).

## Rubric (per file)

1. Do doc comments state contracts (inputs, invariants, units) rather
   than narrate implementation or history?
2. Any references to plans, conversations, requests, or "the classic/
   reference engine" where compatibility is not the actual constraint?
3. Naming: would the identifier mislead without the design docs open?
4. Comment density and idiom match the file's surroundings?
5. For design docs: can a newcomer find the CURRENT architecture
   without reading the landing-note history?

## Constraints

- Comment/doc/naming changes only; zero behavior diffs (renames limited
  to internal symbols, verified by the component tests compiling).
- Fixes batched per directory, one implementer each, reviewed under the
  standing process.

## Steps

1. Fan out readers (one per directory group: R/, src/bartcore/, bridge,
   tests, docs/design/, man/) with the rubric; each returns a findings
   list with file:line, quoted text, and the objection - no fixes yet.
2. Fable dedups and adjudicates into a fix list; VD skims it (taste
   calls are his).
3. Fix wave: mechanical comment/doc edits per directory. For
   core-generalization.md: add a short "current architecture" section
   up top (or split the history into an appendix) rather than rewriting
   the record.
4. Add the summary of what changed to this file; close.

## Verification

- R CMD INSTALL . + component tests + full tinytest (proves
  comment-only).
- git diff --stat shows only comments/docs/man (spot-check any .R/.hpp
  hunks).

## Status (2026-07-08)

CLOSED. Six readers (opus on the engine headers, sonnet elsewhere)
swept the R surface, engine, bridge + C API + touched support-lib
hunks, tests, design docs, and man/README/NEWS; benchmarks scripts, CI
workflows, and the build system were excluded and recorded as such.
Landings, in order: feda8fe (two REAL pre-existing R bugs the reader
surfaced: a live browser() in predict.rbart and predict.bart's error
message listing predict.rbart's type options), 3477c57 (engine header
contracts: deleted-symbol references, plan-step leaks, the stale
phase-1 scope line, a false "thinning-free" claim), c0d482d (bridge
bcfParams length-8 contract, stale step-4 refusal claim, tests'
scratchpad cites restated as inline derivations, gp-leaves.md status +
dead references), 0b66a15 (design docs: core-generalization.md gained
a verified "Current architecture (as shipped)" section and in-place
corrections incl. BartData -> ColumnStore; linear-leaves.md's Model
moved above its deltas; public-surface.md orphan fragment + mislabeled
section; weighted-logistic.md stray tool-call tag), 6bfca8a (VD-
directed sweep: 31 classic-engine attributions deleted, 3 reworded to
the real constraint, 2 kept as "the pre-1.0 engine", CGM literature
reference preserved; ClassicSampler renamed ConstantLeafSampler, the
only code-token change, equivalence exact 18/18), 7498154 (seven
code-verified Rd accuracy bugs, typos, xbart weights completeness,
summary.bart precision, NEWS.Rd codenames -> user-facing wording,
checkRd 14/14 clean). VD taste rulings: sweep aggressively (done),
restate scratchpad derivations inline (done), rename ClassicSampler
(done), leave the bartcoreSampler*/bartcore* R naming, no
core-generalization restructure beyond the new section, leave
mia-missingness ordering. Dropped mid-review by CI rather than a
reader: the warm-start installForests refusal-path leak (a91e5d4) and
the pinned-rbart-draws cross-platform tolerance + missing pkgdown
index entry (29b6243). The standing per-review check in
docs/plans/README.md continues; architecture-doc inherits the
current-architecture section as its seed.
