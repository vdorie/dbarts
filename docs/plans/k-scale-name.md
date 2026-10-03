# k-scale-name: the value k is measured against is k.scale

Status: LANDED 2026-10-02 (499cb3c3) under dec-B201 in [decisions.md](../decisions.md).
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

## Landing

Landed in 499cb3c3.

- Renamed: the reader's field `anchor` is `k.scale` in `getLeafPrior()`, in the fit's `leaf.prior`, in the
  readers of it (`extract`'s leaf.prior.sd, the diagnostics) and in the docstrings, the manual pages and the
  design notes that name it; the model attribute `response.anchor` is `response.range`, and the bridge's
  refusal message for a bad record names it. "Relative to the data's anchor" is "relative to the data's scale"
  in the messages and the manual. Values are unchanged.
- One test, test-k-scale-name.R: the names carry `k.scale` and no `anchor` on a single-forest sampler, a
  multi-forest sampler and a fit; the attribute is `response.range`.
- Hits of the word kept in R and the manual: internal names (`recordAnchor`, `applyAnchor`, the control's
  per-forest `anchor` record and the bridge's `map.anchor`, which no reader returns) and comments on them; "re-anchor"
  and "anchored", the verb for a response-transform change; internal comments on the model's `prior.scale`
  and the calibration map's leaf scale. In the C++ tests (tests/cpp) the engine's own names stay. In
  inst/tinytest the remaining hits are the control's `anchor` record and comments.
- Gates: tinytest 12341 tests, 0 failures, 0 new warnings; tests/cpp passes; lintr, air, check-rc-codoc,
  check-win-drift and check-doc-freshness pass. Against a chain built on this tip, stan4bart's suite at home
  passes (570 tests), and bartCause's and treatSens's testthat suites pass. No consumer reads the field or the attribute.
- Review corrections, in 56f3fa9e: the two vignettes and the backfit-exact gate read `k.scale`; the bridge's refusal
  of a non-positive leaf-prior sd no longer says anchor; the nbinom manual line and the prior-defaults note
  name the `prior.scale` slot where the slot is meant and qualify the k-named case. Afterwards tinytest passes,
  12341 tests, tests/cpp and the lint chain pass, every gate of exact-gates.yaml passes in quick mode
  (the two cross-host equivalence compares need the reference build and were not run), and the package builds
  with its vignettes, whose purled code runs.
