# restore-missing-direction: a state install drops a missing direction its column no longer routes

Status: LANDED 2026-10-05 (5cac562f, ee7ce3e5).

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. Every `setState`, `copy`, reload, `installTrees` and warm start accepted today installs no rule
with a stale missing direction, so its draws are bitwise unchanged; installs refused today now proceed.
window: pre-release. The same defect restore-empty-leaf closed, on a feature new in 1.0-0.
budget: ~300 lines (C++ ~60, tests/cpp ~80, tinytest ~120, manual and records ~40). Plans have run 1.5-2x low.

## Goal

Installing a saved state whose rules send a missing value to a side succeeds when the sampler's column no
longer holds a missing value: the direction is dropped, the handling a forced `setPredictor` already applies.
A sampler can always restore its own state and `copy()` itself after its predictors changed.

## Context

- A rule records which side a missing value goes. Whether a column routes one is derived from the current
  predictors and re-derived by `setPredictor` and `setData`.
- [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) refuses an ordinal rule that sends missing right, and a
  categorical rule whose missing bit is set, when the column has no missing value now. Every install path
  builds through it, so a state stored while a column had missing values is refused once
  `setPredictor(x, forceUpdate = TRUE)` fills them ("state is not consistent with this sampler"): by
  `setState`, by `copy()` and by a reload, 10 of 10 on an ordinal and on a categorical column. A state from
  another sampler over the same rows whose column had missing values is refused the same way.
- The live paths already drop such directions
  (retired: [`Chain::dropStaleMissingDirections`](../../src/bartcore/chain.hpp), gone since
  setstate-force-update Part A, from the forced update and `setData`): the bit routes nothing without a missing observation, so clearing it moves nothing.
- 0.9-34 had no missing predictors, so no released behaviour is involved.
- Since restore-empty-leaf, a restore merges leaves no row reaches. An outer step that stores, changes the
  predictor, and on rejection restores is therefore exact only in the order "re-set the predictor, then
  `setState`" (20 of 20); the other order merged in 4 of 20 runs. The manual does not say so.

## Constraints

- No accepted install changes: the seeded snapshot files, equivalence baselines and exact gates are untouched.
- Dropping a direction changes no routing, so it is not a merge and raises no signal.
- A state that is malformed for other reasons is still refused, a categorical rule whose mask names levels past
  the column's count included.
- Saved draws replay by value against the predictors they are given and are not touched.
- No NEWS: missing predictors are new in 1.0-0.
- Out of scope: states drawn on a different cut grid or standardization (state-frame-prior).

## Steps

1. Let every install path (`setState`, `copy`, a reload, `installTrees`, `warm.start`, and the validators that
   build scratch trees before a live tree is touched) accept a rule whose missing direction the column cannot
   route, and clear the direction in the installed tree, mean and variance forests alike. tests/cpp: an ordinal
   and a categorical rule, a pooled-mask column, a variance tree, and a multi-forest sampler.
2. tinytest: a sampler's own stale state restores through `setState`, `copy()` and a save-reload after
   `setPredictor` fills a column's missing values, ordinal and categorical, then runs, and its restored trees
   predict as the stored ones do on the filled rows; a state from another sampler over the same rows does the
   same; an install accepted before the change is bitwise identical in its next draws.
3. Mutation check: reinstating the refusal makes the new tests fail.
4. Manual, `setState`: to undo a predictor change exactly, re-set the predictor first and then restore; the
   other order restores against the changed predictor and merges what it leaves empty.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes; the
  four seeded snapshot files pass unchanged on a reference build.
- The equivalence compare in statistical (z) mode against the current baseline: every scenario is expected to
  report identical draws.
- ASan and UBSan on tests/cpp.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean.

## Landing note

Landed 2026-10-05 as 5cac562f (the install accepts and drops the direction, in the one build every install
path and validator shares) and ee7ce3e5 (review corrections: the manual's undo recipe, the pooled gate's test,
two-chain and sparse-backed cases). Review SOUND WITH CORRECTIONS, then SOUND. Full tinytest 13746 results, 0
failures; tests/cpp passes, clean under ASan and UBSan; the four seeded snapshot files pass on the reference
build; equivalence in statistical mode 55 of 55 identical, BCF 15 of 15 and multinomial 11 of 11; installs
accepted before the change draw identically across 24 sampler kinds.

Three things differ from the plan as written.
- A pooled categorical rule keeps its bit, as the forced update leaves it: the bit sits in the pool words,
  which rule equality compares, and clearing it at install alone would make a restore differ from the forced
  update. Such a rule routes nothing wrongly; tree moves leave its subtree alone until the rule is redrawn.
  Clearing it on both paths changes draws after a pooled column is filled and is its own item (TODO).
- A categorical mask that sends every reachable level one way with no missing flag, on a column that never
  held a missing value, was refused and now installs and is merged: a state does not record which columns held
  missing values, so the build cannot tell it from a rule that split the missing value from every level.
- "Predict as the stored ones do" holds for trees the filled predictors left no empty leaf in; the others are
  compared with the forced update.

To undo a predictor change exactly the old predictor goes back with `forceUpdate = TRUE` before `setState`; a
factor column does not take missing values back through a column update, so that case goes through `setData`.
