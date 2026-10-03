# k-scale-name: the value k is measured against is k.scale

Status: PLANNED 2026-10-02 under dec-B201 in [decisions.md](../decisions.md).
Follows [state-not-model.md](state-not-model.md), landed.

agent: sonnet implementer, one; opus reviewer.
rng: NEUTRAL. A rename: no draw, no stored state and no engine value changes.
window: pre-release, before the 1.0-0 merge.
budget: ~250 lines (R ~40, manual ~40, tests ~120, records ~50).

## Goal

The leaf-prior reader's field `anchor`, the value the engine's k is relative to so that the spread in force is
that value over `getK()`, is named `k.scale` everywhere a user meets it: [`getLeafPrior`](../../R/dbarts.R), the
leaf-prior list a fit carries, the manual and the messages. Its value and k's are unchanged (dec-B201).

## Constraints

- k stays the engine's k on every prior. Nothing about what the field holds moves.
- No engine or bridge change. Internal C++ and bridge names (a calibration map's `map.anchor` key among them)
  stay unless they reach a user.
- The model attribute `response.anchor`, the sampler's recorded response transform, is renamed
  `response.range`: it holds the response's (min, max), and the word goes from the surface with the field.
  No shipped object carries either name, so neither is read under its old name.
- Prose: where the manual or a message says "anchor" for this value, it says k.scale or, under a k-named prior,
  the data's scale; "relative to the data's anchor" becomes "relative to the data's scale".
- Present-facing docs that name the field (docs/architecture.md, the design notes for the leaf
  prior, the state and the fit) follow. Plans and the ledger keep their history and are not rewritten.
- No NEWS entry: the reader and the field are new in 1.0-0.

## Steps

1. R: the reader, the fit's leaf-prior list, every internal reader of the field (`extract`, the diagnostics),
   the docstrings and the messages; the model attribute rename.
2. Manual: `dbartsSampler-class.Rd`, `bart.Rd`, `bartBT.Rd`, `dbartsPriors.Rd` and any other page naming the
   value.
3. Tests: every case reading `$anchor` or the attribute; one case that the reader's names carry `k.scale` and
   no `anchor`, on a single-forest and a multi-forest sampler and on a fit.
4. Records: the present-facing docs above, the index row, the Landing note.

## Verification

Against a private library:

- `tinytest::test_package("dbarts")` passes with no new warning; `cd tests/cpp && make && ./test_bartcore`
  passes.
- `git grep -n -w anchor -- R man inst/tinytest` shows no user-facing use of the word for this value; each
  remaining hit is internal and listed in the Landing note.
- stan4bart's, bartCause's and treatSens's suites pass against a private-library chain built on this tip,
  and no consumer reads the field (`git -C <repo> grep -n anchor`).
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status.
