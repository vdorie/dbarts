# monotone-unforced-refusal: an unforced update that would put a monotone tree out of order is refused

Status: PLANNED (dec-B278).

agent: opus implementer, one (engine, tests); opus reviewer.
rng: by call sequence.
- POSTERIOR-CHANGING, as a correction, on one sequence: an unforced predictor update (whole matrix, by
  column, row by row) on a sampler with a monotone constraint that brings an unordered factor column its
  first missing value while some tree's leaf values would then be out of order. Today it is accepted and
  every leaf of such a tree is set to zero; afterwards it is refused and the sampler is as it was. The call
  itself draws what it draws today: nothing in the whole-matrix and column forms, one scan order in the row
  forms.
- NEUTRAL, bit for bit, for everything else: every sampler without a monotone constraint (the check
  compiles out), every forced update, `setData`, `setCutPoints`, `setState` and a warm start, every update
  refused for an empty leaf, and every unforced update on a monotone sampler that leaves each tree in order.
Proved by the bitwise gates below (no baseline scenario, exact gate or snapshot file pairs a monotone
constraint with a predictor update; checked by search at the tip), by a seeded digest of accepted and
refused updates on the base and slice builds, and by the new tests' identity with an untouched twin.
window: pre-release, before 1.0-0 (dec-B278). Serial with [state-zero-weight-rows.md](state-zero-weight-rows.md)
and [leaf-conversions.md](leaf-conversions.md): all three edit chain.hpp, sampler.hpp, the bridge file and
the sampler's manual page, in other functions and other items. Recommended order: after
state-zero-weight-rows, then before leaf-conversions or between its two pushes; see Calls.
budget: ~450 lines (C++ engine ~60, bridge ~6, R none, tests/cpp ~150, tinytest ~170, manual, design note,
comments and TODO ~65), upper figure 800. Plans have run 1.5-2x low; the ledger entry this comes from
estimated 200.

## Goal

On a sampler with a monotone constraint, an unforced predictor update never changes a fit it does not
refuse. One that would leave a tree's leaf values out of order along a constrained predictor returns
`FALSE`, or `FALSE` for each row that would do it, and the sampler is as it was, exactly as for an update
that would empty a leaf. A forced update, `setData`, `setCutPoints`, a warm start and a `setState` that has
to merge leaves keep completing and keep resetting such a tree to zero without a message. The manual says
both.

## Context

All numbers were run on the tip's build (shipped mode). "In order" is what the engine calls feasible: every
pair of leaves the constraint relates has its values the right way round.

- What the order reads. [`MonotoneLeafGeometry`](../../src/bartcore/model.hpp) gives each leaf, per split
  variable, a code interval and whether a missing value reaches it (numeric and ordered columns) or a set
  of levels with the missing position (unordered factors), from the tree's rules, the column's cut and
  level counts and its has-missing flag ([`hasMissing`](../../src/bartcore/data.hpp)).
  [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp) adds the leaf values and reads no partition, so
  where rows sit never matters. An unforced update changes no rule, level table or cut count (a refresh
  keeps the count the column holds; [`ColumnStore::cutsWouldRemainValid`](../../src/bartcore/data.hpp)
  stops one that cannot), so it can change the order only through a has-missing flag.
- Which flag change relates new leaves. A rule on a column that holds no missing value carries no direction
  for one: none is drawn, and [`Tree::dropStaleMissingDirections`](../../src/bartcore/tree.hpp) clears what
  a column's last missing value leaves behind. A first missing value therefore goes left at every rule on
  its column.
  - Numeric or ordered column: a leaf it reaches is left of every rule on that axis above it, so the leaf's
    interval starts at the lowest code, and two such leaves already share that code. Nothing is added.
  - Unordered factor: the missing position joins the left level set of every rule, and two left sets need
    not share a level (a and b left in one branch, c and d left in another). Two leaves that shared no
    level then share the missing position, and become related if they touch along a constrained predictor.
  - A flag that clears only removes relations, and an in-order tree stays in order.
- Measured through the engine's own geometry: 20000 random trees of 2 to 10 births over two constrained
  numeric columns and a free numeric, ordered, 4-level and 70-level column, 110810 related leaf pairs.
  Setting one flag adds 371 pairs in 217 trees for the 4-level factor and none for the other five, nor for
  the four numeric and ordered columns set together. Clearing one (stale directions dropped) removes 200,
  162, 398, 428, 3177 and 1 pairs and adds none. No pair reverses. The 70-level factor adds none because
  two random halves of 70 levels share a level; the argument does not exclude it.
- Measured from R: 3780 unforced trial updates of 35 kinds on copies of 27 running fits (a first missing
  value in a numeric, an ordered and each of two constrained columns, by column, whole and row by row; a
  factor row moved to a level with no rows; a numeric, an ordered and an unordered column losing its last
  missing value; a cut refresh on a free and a constrained column under both cut rules). 3620 accepted,
  160 refused for an empty leaf, no tree's leaf values changed.
- How often a fit meets it: of 3600 unforced whole-matrix updates each bringing a 4-level factor its first
  missing value, on fits whose slope in the constrained predictor flips with the level (1, 5 and 20 trees,
  both priors), one reset one tree. The tests install a hand-built tree.
- What the tip does with that tree (the constrained x1 at the root, the factor split two ways below it,
  leaf values in order until the factor has a missing value), each call silent unless an error is shown:

  | call | returns | the tree |
  |---|---|---|
  | whole matrix of codes, `forceUpdate = FALSE` | `TRUE` | every leaf 0 |
  | `updatePredictorPerObservationJointly`, codes | every row `TRUE` | every leaf 0 |
  | whole matrix, `forceUpdate = TRUE` or missing | `TRUE`, invisibly | every leaf 0 |
  | `setData` with the missing value | `NULL`, invisibly | every leaf 0 |
  | `installTrees` into a sampler whose factor has one | `NULL`, invisibly | every leaf 0 |
  | `setState` of that state into such a sampler | error, `leaf values violate this sampler's monotone constraint` | untouched |
  | by column or `"partial"`, the factor as a factor | error, `has missing values, which its training values do not` | untouched |

  With leaf values that stay in order the two unforced rows return the same and keep the values; two
  chains behave as one. Where a merge of emptied leaves puts a tree out of order, a forced update by
  column, `setData`, `installTrees` and a merging `setState` (`FALSE`, invisibly) reset it in silence too.
  `setCutPoints` completed in silence on five grids and reset nothing.
- Which routes reach it from R. By column and `"partial"`,
  [`codeCategoricalColumnUpdate`](../../R/bartcore.R) refuses a factor column's first missing value by name
  before the engine is called, forced or not, as the manual's `x` item says. The whole-matrix form takes a
  numeric matrix with the factor as codes, and is forced unless `forceUpdate = FALSE` is stated.
  [`updatePredictorPerObservationJointly`](../../R/updatePredictorPerObservationJointly.R) takes codes. The
  engine's column form and one-sampler row form are reached from tests/cpp alone; the flat C header has no
  entry for a predictor update.
- Which samplers carry the constraint. Built and run: gaussian, probit, logistic, Student-t,
  negative-binomial, ordinal and aft responses, weights, two chains, `interactions`, `blocks`, one declared
  forest. Refused by name: a second forest, a variance forest, linear and gp leaves, a multinomial
  response. So one forest per chain is all there is.
- The unforced update is two phases ([`Sampler::runPredictorTransaction`](../../src/bartcore/sampler.hpp),
  [`Sampler::revalidateAllChains`](../../src/bartcore/sampler.hpp)). The store first takes the new values
  ([`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp), [`SubsetUpdate`](../../src/bartcore/sampler.hpp)):
  codes, has-missing flags, the cut grid under `updateCutPoints`, sparse storage, raw copies. Phase one,
  [`Chain::revalidateTrees`](../../src/bartcore/chain.hpp), reads each tree's leaf values, re-routes its
  rows and asks whether every leaf is occupied; it writes no fit, leaf value or rule. A failure restores
  the store from its snapshot and re-routes ([`Chain::repartitionTrees`](../../src/bartcore/chain.hpp)).
  Phase two, [`Chain::rebuildFitsFromParameters`](../../src/bartcore/chain.hpp), is not undone: it drops
  stale missing directions, rewrites the row-to-leaf map and the cached fits, and holds today's call of
  [`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp). The order reads the same before
  and after that drop: a column with its flag clear has no missing position, whatever its rules say.
- The standard a refusal meets today. 36 updates refused for an empty leaf (by column, two columns, whole
  matrix; with a first missing value; with `updateCutPoints = TRUE`; on a column losing its last missing
  value; monotone and plain samplers, one and two chains): `FALSE`, no warning, and afterwards the stored
  state is byte for byte an untouched twin's, `data@x`, predictions and the next five draws identical. A
  refused whole-matrix update is rolled back as a refused column update is. The fixture's 7 rules that send
  a missing value right keep doing so. One thing differs, seen from tests/cpp: after the two re-routes a
  leaf's rows can be held in another order than the twin's. No state, prediction or draw above shows it.
- Row by row. [`UpdateSessionImpl`](../../src/bartcore/sampler.hpp) judges rows one at a time in a scan
  order drawn from the first chain's generator (the first sampler's, in the joint form), and an accepted
  row counts against the next: of a leaf's two rows both moved out, the first is kept under three of six
  seeds and the second under the other three, where the same column by column or whole is refused
  outright. The joint form declines a row in every sampler when one declines it. The scan order is drawn
  even when no row moves (the next three draws differ from a twin's by 0.49). A call that refuses rows ends
  bit for bit a twin given the column the call left, 4 of 4.
- The session has no rollback: [`ColumnStore::setCell`](../../src/bartcore/data.hpp) writes each accepted
  cell and sets the flag at the first missing one, and the closing rebuild is expected to hold. When it
  does not, the bridge raises `$setPredictor produced a tree with an empty leaf`
  ([`bartcore_updatePredictorPerObservation`](../../src/R_interface_bartcore.cpp),
  [`bartcore_updatePredictorPerObservationJointly`](../../src/R_interface_bartcore.cpp)) with the cells
  written. A row that would break the order has to be stopped at the row.
- R needs nothing. [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) puts `data@x` back or never installs
  it on `FALSE`, returns `FALSE` visibly, and leaves a refused row's old value in `data@x`. `FALSE` is tied
  to an empty leaf in words in three places: the two bridge messages and the manual.
- The placement, tried on a scratch copy of the tip (58 lines in chain.hpp and sampler.hpp): the tree is
  refused whole and by column with the flag clear again, and codes, trees, the row-to-leaf map, fits and
  three further sweeps are bit for bit an untouched twin's, one and two chains. Row by row the two rows
  bringing the missing value are refused, 398 install, and the sampler is a twin given the column without
  them. Forced updates and `setData` still reset. Each of five mutations fails it; with the row guard
  removed the session ends invalid, the error above.
- Pins that move: [`testMonotoneMissingArrives`](../../tests/cpp/test_monotone.cpp) asserts the acceptance
  and the reset on the whole-matrix and the row path. No tinytest pins the unforced reset. Pins that stay:
  ["a forced predictor update that flattens x2"](../../inst/tinytest/test-monotone.R), the hand-built
  collapse above it, and that file's warm-start blocks.
- Consumers (read only). bairrtt calls the joint form on a numeric latent column; stan4bart passes
  `monotone` through and calls no predictor update. Neither can reach the refusal.

## The rule

An unforced update on a sampler with a monotone constraint is judged on two things together, for every tree
of every chain it can move: every leaf keeps a row, and the leaf values are in order as the store stands
with the new values. Either failure refuses it.

1. Whole matrix and by column. The order is judged in phase one, with the empty leaf, before anything that
   is not rolled back. A refusal is the empty-leaf refusal: `rolledBack` from the engine, `FALSE` from R, no
   error, no warning, no leaf value touched, the store restored (has-missing flags and cut grid included),
   no direction dropped.
2. Row by row, single and joint. A row whose new value is missing, in a column that holds no missing value,
   is refused when a missing value in that column would leave any tree of any chain out of order; in the
   joint form, any tree of any sampler, and the row is declined in all of them. The answer does not depend
   on the row or on the rows installed before it, so it is found once per call: either every such row is
   refused and the column still holds no missing value, or the first one installs and the rest are judged
   on occupancy alone. Every other row is judged as today, and the scan order is drawn as today.
3. The calls that always complete are not touched: `setPredictor(forceUpdate = TRUE)` (the default when the
   whole matrix is replaced), `setData`, `setCutPoints`, a warm start and a `setState` that merges set every
   leaf of a tree left out of order to zero, drawing nothing and saying nothing.
4. The order is judged on every unforced update of a monotone sampler, not only when a flag rises. It is
   the test phase two runs today on the same trees, moved, so an accepted update costs what it costs now.

## Constraints

- Every update accepted today that leaves each tree in order is accepted, bit for bit as now.
- A refusal leaves what an empty-leaf refusal leaves, held to the same twin: state, `data@x`, predictions
  and next draws identical.
- [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp) gains no value; `rolledBack` covers both
  reasons. No facade virtual, state format, flat C entry or R code changes. `--preclean` on every install
  all the same: two engine headers change.
- R's refusal by name of a factor column's first missing value on a column update stands, forced or not.
- The whole-matrix default stays forced, the column default unforced.
- No message, warning or condition is added to any forced route (dec-B278).
- No NEWS entry: no monotone fit has been in a release.
- Out of scope: see the end.

## Steps

1. Engine, whole matrix and by column. [`Chain::revalidateTrees`](../../src/bartcore/chain.hpp) asks, for
   each tree it re-routes, occupancy and then the order at the values it has just read, through one const
   reader on the chain that wraps [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp) and is true off
   the monotone leaf. The call of
   [`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp) in
   [`Chain::rebuildFitsFromParameters`](../../src/bartcore/chain.hpp) goes, with its comment. Its four
   other callers keep theirs ([`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp),
   [`Chain::applyNewData`](../../src/bartcore/chain.hpp),
   [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp),
   [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp)). The comments on
   [`Sampler::setPredictor`](../../src/bartcore/sampler.hpp),
   [`Sampler::runPredictorTransaction`](../../src/bartcore/sampler.hpp) and
   [`Sampler::revalidateAllChains`](../../src/bartcore/sampler.hpp) name both reasons; their code stands.
2. Engine, row by row. [`UpdateSessionImpl`](../../src/bartcore/sampler.hpp): before a row is staged, when
   its new value is missing (the test [`ColumnStore::setCell`](../../src/bartcore/data.hpp) marks the column
   by) and the column's flag is clear, ask once per session whether every tree of every chain is in order
   with the flag set: set it, ask each chain over all its live trees at their current values, clear it,
   keep the answer. A "no" refuses the row with nothing staged. The comments on
   [`PredictorUpdateSession`](../../src/bartcore/sampler.hpp),
   [`Sampler::updatePredictorPerObservation`](../../src/bartcore/sampler.hpp) and
   [`updatePredictorPerObservationJointly`](../../src/bartcore/facade.hpp) say what a refused row is. The
   joint form's code stands: it installs a row only when every session takes it.
3. Bridge. The two messages for a session that ends invalid stop naming one cause: `$setPredictor left a
   tree invalid`, and the joint form's equivalent. No test pins either text.
4. tests/cpp. [`testMonotoneMissingArrives`](../../tests/cpp/test_monotone.cpp) turns over, on its own
   tree, one and two chains; "fails today" names what the tip does. The row given the missing value is one
   it moves to another leaf (x1 above the cut, level a or b), so a rebuilt fit shows.
   - Whole matrix and by column, unforced: `rolledBack`; the factor's flag clear; codes, flags, flat trees
     with their values, the row-to-leaf map and the cached fits identical to a twin built from the same
     seed, and identical again after three sweeps of each. Fails today: accepted, reset.
   - The same with the out-of-order tree in the second chain alone.
   - The same refused with `updateCutPoints` set: the cut grids are the twin's.
   - Row by row, a column with two missing values and one other changed row: both missing rows refused,
     the 398 others installed, the session valid, the flag clear, and trees, map, fits and three sweeps
     identical to a twin given the column without the two missing values. Fails today: 400 installed, reset.
   - Two samplers swept as [`updatePredictorPerObservationJointly`](../../src/bartcore/facade.hpp) sweeps
     them, the one in order listed first: the row is declined in both, both end valid, neither flag is set.
   - Leaf values that stay in order: accepted whole, by column and row by row, values kept, flag set. A
     column that already holds a missing value takes another. Both hold today.
   - Forced whole and by column, and [`Sampler::setData`](../../src/bartcore/sampler.hpp) with the missing
     value: taken, every leaf zero. Holds today.
   - The geometry claim, on random trees over a store with no missing value (the file's own tree builder):
     setting the flag of a numeric or an ordered column changes no relation, setting an unordered
     factor's adds some and removes none, and clearing one adds none.
5. tinytest, a new file `test-monotone-unforced.R`: y on x1 and a 4-level factor, x1 increasing, one tree,
   seeded, the tree above written into the stored state as test-monotone.R writes its hand-built tree (the
   level mask in the machine's byte order).
   - Whole matrix of codes, `forceUpdate = FALSE`: `FALSE`, visible, no warning (counted with
     `withCallingHandlers`); `data@x` identical to before, class included; the stored state identical to a
     twin's, predictions and the next five draws identical. Again with `updateCutPoints = TRUE`, the cut
     points identical, and on two chains with the tree out of order in the second alone. Fails today:
     `TRUE`, leaves zero.
   - The joint form with codes, two missing rows and a changed row: `FALSE` at exactly the two rows, no
     warning, `data@x` holding its old value at those rows, the sampler identical to a twin given the
     column without them. With a plain sampler listed first, both designs hold the same column. Fails today.
   - Refused again on a second call with other values and the same missing row; accepted once a state
     install has put the values in order. This is the manual's sentence on proposing again.
   - Values that stay in order: `TRUE` and every row `TRUE`, leaf values unchanged, the missing value in
     `data@x`. A sampler built with a missing value in the factor takes another. A first missing value in
     x1, by column, whole and `"partial"`: accepted, leaf values unchanged. All hold today.
   - The calls that complete: `forceUpdate = TRUE` and `forceUpdate` missing (invisible `TRUE`), `setData`,
     `installTrees` into a sampler whose factor has a missing value: no warning, every leaf zero, the fit
     monotone at every level after a sweep. On a numeric fixture whose forced update merges two leaves
     below a third: by column forced and a merging `setState` (`FALSE`, invisibly) reset in silence, and
     `setCutPoints` completes in silence with the fit in order. All hold today.
   - By column and `"partial"` with the factor given as a factor: the error by name, forced or not, the
     stored state unchanged. Holds today.
6. Mutations (Verification): apply, install with `--preclean`, run, report the failing counts, revert,
   `touch`.
7. Records.
   - Manual, [`dbartsSampler$setPredictor`](../../man/dbartsSampler-class.Rd), the `forceUpdate` item, after
     the sentence on what `"partial"` returns: "Under a `monotone` constraint a new predictor can also
     leave a tree's leaf values out of order along a constrained predictor, and the three values treat that
     as they treat an empty leaf: `FALSE` refuses the update and rolls it back, `"partial"` refuses the rows
     that would do it and installs the others, and `TRUE` completes the update, setting every leaf value of
     such a tree to zero without a message; the next iteration draws them again."
   - The Value paragraph's first three sentences become: "For `setPredictor`, `TRUE` if the new predictor
     was installed and `FALSE` if it was refused. An unforced update is refused, with no error and no
     warning, when it would leave a leaf empty in any tree of any forest of the sampler (see 'Multi-forest
     and heteroscedastic predictor mutation') or, under a `monotone` constraint, leave a tree's leaf values
     out of order. A refused update is rolled back, whether it named columns or replaced the whole matrix:
     the sampler and its predictors are as they were before the call. A forced update always installs and
     returns `TRUE` invisibly." The sentence on `"partial"` stands.
   - [`monotone`](../../man/monotone.Rd), a paragraph after the one that ends "nothing is claimed where the
     constrained predictor itself is missing": "New predictor values given to a sampler can leave a tree's
     leaf values out of order. The calls that always complete - `setPredictor` with `forceUpdate = TRUE`
     (its default when the whole matrix is replaced), `setData`, `setCutPoints`, `installTrees`, and a
     `setState` that has to merge leaves - set every leaf value of such a tree to zero, without a message,
     and the next iteration draws them again; that suits burn-in and is not a draw from the posterior. An
     unforced update - `setPredictor` with `forceUpdate = FALSE` (its default when columns are named) or
     `"partial"`, and `updatePredictorPerObservationJointly` - is refused instead, as one that would empty
     a leaf is: it returns `FALSE`, or `FALSE` for each row that would do it, and the sampler is left as it
     was. One unforced change can do this: the first missing value in an unordered factor column that has
     none, which gives the column a position that can relate leaves the order did not. While a tree is in
     that position every update that brings the missing value is refused, whatever its other values, so
     code that proposes again until an update is accepted should run the sampler, or leave the missing
     value out, before the next proposal."
   - [updatePredictorPerObservationJointly.Rd](../../man/updatePredictorPerObservationJointly.Rd), Details:
     "In a sampler with a `monotone` constraint a value is also declined, in every sampler, when installing
     it would leave a tree's leaf values out of order."
     [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd): the install mask is "`FALSE` where the new value
     would empty a leaf or, under a `monotone` constraint, leave a tree's leaf values out of order, and was
     rolled back".
   - Design record: [monotone.md](../design/monotone.md), a dated section on predictor changes under the
     constraint: the rule, the argument for which flag change relates new leaves with the measured counts,
     and the rate in fits; its Status line gains the revision.
     [Stage 3b: leaf geometry](monotone-exact-birth-death.md#stage-3b-leaf-geometry) is a landed record and
     is left as written.
   - TODO: the entry `monotone-unforced-refusal` names this plan.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; `tests/cpp` builds and passes, clean under ASan
  and UBSan; the full tinytest suite; the new file under ASan on the R-loaded path.
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares
  are bitwise, every scenario reporting identical draws, counted per scenario with no `max |z|` line: 55
  against `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded: the two monotone scenarios call no
  predictor update and the scenarios that call one carry no constraint. A scenario that is not identical is
  a finding, not a re-record.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged; none pairs the
  constraint with a predictor update. No new exact gate: a refusal has no posterior of its own to state,
  and is held to identity with an untouched twin, the standard the empty-leaf refusal meets.
- One script on the base and slice builds digesting seeded states and draws: plain and monotone samplers,
  one and two chains, through unforced updates accepted and refused for an empty leaf (whole, by column,
  two columns, row by row, with and without `updateCutPoints`, a first missing value in a numeric, an
  ordered and the constrained column, a last missing value lost), forced updates, `setData`,
  `setCutPoints`, a copy and a reload. Equal.
- Mutations, each expected to fail the named test:
  - the order is not asked in phase one and phase two resets as today: tests/cpp "whole matrix and by
    column, unforced" and tinytest "whole matrix of codes";
  - the order is asked after the fits are rebuilt, the update then reported refused: the same two, at the
    row-to-leaf map and the fits in tests/cpp and at the next draws in tinytest;
  - the has-missing flags are not put back, on the whole-matrix restore and on the column restore each:
    tests/cpp "the factor's flag clear" and the three sweeps; tinytest's next five draws;
  - the row guard is removed: tests/cpp "row by row" (the session ends invalid) and tinytest's joint form
    (the bridge raises);
  - the session leaves the flag set after asking, or asks without setting it: tests/cpp "row by row" and
    the joint case;
  - the session refuses every missing value, or phase one refuses whenever a flag rises with the values
    unread: tests/cpp and tinytest "values that stay in order" and the column that already holds one;
  - only the first chain is asked, in phase one and in the session each: tests/cpp "the second chain
    alone" and tinytest's two chains;
  - a forced update is refused when the order would break: tinytest "the calls that complete", tests/cpp
    "forced whole and by column", and the two forced blocks of test-monotone.R;
  - `setData` refuses, or leaves the values, when the order breaks: tinytest and tests/cpp `setData`.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; `R CMD check
  --as-cran` on a tarball from a clean copy (man/ changes).
- Not a hot-path change: nothing is added to a sweep. An unforced update on a monotone sampler runs the
  order test once per tree as it does today, in the other phase; a row-by-row update tests one more
  condition per row and runs the order test at most once per call.

## Out of scope, and where it goes

- A reason with the refusal. `FALSE` does not say whether a leaf would empty or the order would break, and
  the second is not cured by proposing other values. The manual says what to do; a result that names the
  cause would be an addition to [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp) and to the
  method's value. Not planned.
- A missing position on every factor axis, which would leave a first missing value nothing to add. Tried
  and rejected in [Stage 3b: leaf geometry](monotone-exact-birth-death.md#stage-3b-leaf-geometry): the
  order then depends on which side of a level split is called left.
- What R takes as a factor column's new values. By column and `"partial"` take labels and refuse a first
  missing value by name; the whole-matrix form and the joint form take codes, a missing one included; the
  joint form given labels reads them as codes; a data frame given as the whole matrix fails in a coercion.
  None is this slice's; the coordinator has the measurements and the text of a backlog entry.
- The order a leaf's rows are held in after a refused update. Every refusal leaves it, today's included.

## Calls made in planning

- The order is judged on every unforced update of a monotone sampler, with no test for a risen flag first.
  The argument in Context says only an unordered factor's first missing value can break it; the check does
  not lean on that, and costs nothing new. The alternative, a check only when a touched column's flag
  rises, saves the test on other updates and makes the argument load-bearing.
- The phase-two reset is removed, not kept as a backstop. Kept, it would hide a failure of phase one behind
  the silent reset the ruling ends for unforced updates.
- No new result value. dec-B278 asks for the empty-leaf refusal; R returns `FALSE` for both and the engine's
  `rolledBack` says what was done. The cost is the first item under Out of scope.
- Row by row, each row bringing the first missing value is refused, not the call: dec-B278 says "`FALSE`
  for that row", and the other rows' moves are valid on their own. The verdict is found once per call
  because it is the same for each such row.
- The session asks every live tree of every chain, not only the trees it caches for the column. A tree
  that does not split on the column cannot change its answer, and asking it costs one pass once per call.
- The row guard is where the session refuses, not its closing rebuild: a session writes cells as it goes
  and has nothing to undo, so a refusal found at the end is the bridge's error.
- The geometry claim gets a tests/cpp check of its own, since the manual states it to users.
- The tinytests get their own file; test-monotone.R keeps its forced blocks untouched as pins.
- The `rng:` line calls the changed sequence posterior-changing and adds no exact gate. Under an embedding
  sampler the reset is a move no acceptance rule accounts for and the refusal is a rejected proposal, so
  the sequence changes law; it is held to the twin identity. Classed as a changed result with neutral
  gates, the gates run are the same.
- Order with the other two plans. state-zero-weight-rows is in implementation and changes a facade virtual;
  this slice goes after it and rebases over nothing it edits: other functions of chain.hpp and sampler.hpp,
  another function of the bridge file, other items of the manual page, another tests/cpp file.
  leaf-conversions edits [`Chain::applyNewData`](../../src/bartcore/chain.hpp), whose reset call stands
  here, and is three times this size in two pushes: this slice before it or between its pushes, not beside
  it.
- The tip against the texts of dec-B278 and the TODO entry.
  - Both say a column's first missing value; it is an unordered factor column's. A numeric or an ordered
    column cannot do it (argument, 20000 trees, 3780 trial updates).
  - The TODO names by column and row by row beside the whole matrix. From R the column form and
    `"partial"` refuse a factor's first missing value with an error before the engine sees it; the
    whole-matrix form with codes and the joint form reach the refusal. The engine is changed and tested on
    all four.
  - dec-B278 lists `setData`, `setCutPoints` and a warm start as the calls that always complete. A
    `setState` that merges resets too and returns `FALSE` invisibly; one whose state is out of order as
    stored is refused with an error. `setCutPoints` completes; no grid tried made it reset a tree.
  - The manual says a refused update is rolled back "if only single columns were replaced". A refused
    whole-matrix update is rolled back as well (16 of 16).
  - "The sampler left as it was" holds for state, predictors, predictions and draws; the order of rows
    within a leaf can differ, for every refusal.
  Nothing in this slice was found done already.
