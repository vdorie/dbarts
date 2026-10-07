# repeated-cut-restore: every split stays on a grid position that holds its value

Status: PLANNED (dec-B285; dec-B283 for what an install refuses).

agent: opus implementer, one (engine, bridge, tests); opus reviewer, told to refute.
rng: by call sequence.
- NEUTRAL, bit for bit, for every fit and every call on a sampler none of whose cut grids repeats a
  value: runs, `setState`, `copy`, a reload, `setData`, `setPredictor`, a warm start.
- NEUTRAL on a grid that repeats a value for every fit that does not restore, replace its data or
  warm-start.
- SHIFTING on a grid that repeats a value, after `setState`, `copy`, a reload, `setData` and a warm
  start: the splits on a repeated value stay on the positions they held, where
  today they move to the first position holding the value, so the draws that follow differ from
  today's. The posterior does not: the prior and the moves are untouched.
- Not a generator matter but a changed result: `setState` returns `FALSE` where it had to choose a
  position; a state whose positions block is malformed is refused; a split the grid cannot hold inside
  its node's interval is merged where today it is left.
Proved by the argument under The rule that on a grid without repeats the resolver returns the position
both callers compute today; by the bitwise gates below (no baseline scenario, exact gate or snapshot
file builds, refreshes or installs a grid with equal neighbours: 0 reports over the 81 scenarios and
four files from a build that prints one line per such grid, run); and by a seeded digest of fits on
grids without repeats through every door, on the base and slice builds.
window: before the merge to main (VD 2026-10-07, dec-B285). After [leaf-conversions.md](leaf-conversions.md),
which is landing in chain.hpp, model.hpp, sampler.hpp and the bridge's state reader and writer; see
What waits on what.
budget: ~920 lines (tree.hpp ~100, chain.hpp ~10, bridge ~100, R none, tests/cpp ~260, tinytest ~380,
manual, design notes, comments and TODO ~70), upper figure 1500. Plans have run 1.5-2x low. The
prototype under Context came to about 80 lines in tree.hpp, 4 in chain.hpp and 105 in the bridge.

## Goal

A restore, a copy, a reload, a rebuild after `setData` or `setPredictor` and a warm start put every
split on a grid position that holds its value and lies inside the interval its node's ancestors leave,
and on the same position wherever the grid there is unchanged. A restore that could not know the
position says so through the value `setState` returns. A grid without repeated values, which is every
grid of an ordinary fit, sees no change of any kind.

## Context

Measured on the build of cc8381b3 ("tip") and on a prototype of this plan built from it as a reference
build ("prototype"); the critic's figures are from de66bd48 and were re-run here on both where marked
"both". A column's cut grid is the list of thresholds a split on it may use; a split's position is its
index on that grid.

- Where a grid repeats a value. The uniform rule builds `n.cuts` points whatever the range
  ([`ColumnStore::fillCutsOverRange`](../../src/bartcore/data.hpp)): 100 points, 1 distinct, over a
  constant column and 100, 5 over one a few ulp wide (run). The quantile rule repeats where two
  midpoints round to one double ([`fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp)): five
  adjacent doubles give 4 points, 3 distinct (run). A constant column given values by `setPredictor`
  with `updateCutPoints` left at its default keeps 20 points, 1 distinct (run). A state brings any
  non-decreasing grid ([`cutGridIsValid`](../../src/bartcore/data.hpp), non-strict form in
  [`Sampler::setState`](../../src/bartcore/sampler.hpp)); `setCutPoints` refuses a caller's repeat
  (dec-B285). The engine draws a rule uniformly over positions, so equal points are separate rules.
- Where a split's value becomes a position: two places, by search over the engine and the bridge for
  every comparison of a cut value (read). [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp), the
  build from a flat tree, takes the first position holding the value: `setState`, `copy`, a reload, a
  same-grid warm start, the first half of a warm start from another grid, the validity checks on
  scratch trees, the rollback of a failed warm start.
  [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp), the move onto another grid, takes the
  nearest value inside the node's interval and for an exact match the first position holding it:
  [`Chain::applyNewData`](../../src/bartcore/chain.hpp) (`setData`),
  [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp) (a warm start from another grid)
  and the variance forest's twins. `setCutPoints` and `setPredictor` keep each split's position and
  map no value.
- The restore, tip. Seven ways a grid repeats, three of them with missing values in the column, 20
  trees, 2 chains, 6 seeds, through `setState`, `copy` and a reload: the restored trees differ from the
  stored ones in position on 6 of 6 seeds in all 21 rows, and the next 30 draws differ from an
  untouched twin's by up to 2.8 on 1 to 6 seeds a row, against 1e-14 on a grid without repeats (run;
  the critic the same on its build).
- A tree the prior gives probability zero, tip. Beside missing values a parent and its child land on
  one position: 22 of 200 restores left a split outside its interval (run); 9 of 120 warm starts from
  kept draws on the same data (run); 9 of 96 warm starts from kept draws recorded before the donor's
  column was regridded (critic; run on both).
- The move, tip. After `setData` with the same data every tree is where it was on a grid without
  repeats, 6 of 6, and none on a grid with them, 0 of 6; a warm start across grids that differ on
  another column only leaves no forest's splits on the repeated column as the donor's, 0 of 6
  (critic; run on both). `setPredictor` with the same column and `updateCutPoints = TRUE` keeps every
  tree, 6 of 6 on both grids.
- The verdict, tip. `setState` of a state on a repeated grid returns `TRUE` with every forest moved, 6
  of 6 (critic; run on both). The flag for it exists: [restore-status.md](restore-status.md) names "a
  split moved onto another grid" as a cause to come.
- What a restore promises. Two samplers restored from one state are bit-identical over the next 30
  draws, and a restored sampler and its copy, 6 of 6 on a repeated grid (critic; run). Against the
  sampler that never stopped, the first draw already differs on any grid, 7e-15 here, because fits are
  summed again: dec-B32, and the manual's "to the last few digits, not bitwise". Under
  `storage = "single"` that gap is 5e-6, and a gp sampler's restore is not a continuation on any grid
  (critic; read).
- 0.9-34 stored a live tree as variable and position and a kept tree by value (read on main,
  `Node::serialize`, `SavedNode::serialize`; not run), over the same uniform grids. The one
  value-coded flat tree is this branch's (b5919a6d).
- The prototype of this plan (run). Restore: 27 rows, no tree moved, the state stored after a restore
  identical to the one installed. 25 kinds of sampler on a repeated grid through three doors (the
  families, DART, a variance forest, linear, gp and monotone leaves, two forests, factor and ordered
  columns beside it, a sparse design, kept trees, missing values on every column): no tree moved. No
  split outside its interval in 200 restores, 120 warm starts from kept draws, 96 from stale kept
  draws, and 200 installs of a state stored on a grid with runs onto the same grid with each value
  once (32 of 200 on the tip, and on a prototype with the interval rule but no merge). `setData` with
  the same data: every tree kept, 6 of 6. The warm start across grids: the repeated column's splits
  the donor's, 6 of 6. A state with its positions removed: `FALSE`, 6 of 6. 40 states with positions
  redrawn among the copies of each value: no split outside its interval, 3 returning `FALSE` (trusted
  as written, 3 of 40 hold a split outside; critic). A positions block that is short, a double, holds
  an `NA`, a negative number or an entry past the grid: refused, the sampler as it was.
- The prototype against the gates (run). Reference build, `EQUIVALENCE_CORES=2`, `--bitwise`: 55 of
  55, 15 of 15 and 11 of 11 identical with no z line, and the four snapshot files pass. tests/cpp 351
  ok. The tinytest suite of cc8381b3, 232 files, 16681 results: 228 files pass untouched; four stop
  or fail, each on a state whose tree blocks it rewrites by hand beside a positions block it leaves
  (["handState"](../../inst/tinytest/test-monotone-unforced.R),
  ["handTree"](../../inst/tinytest/test-monotone.R),
  ["replaceFirstVarianceTree"](../../inst/tinytest/test-heteroscedastic-mutation.R),
  ["splicedState"](../../inst/tinytest/test-heteroscedastic-warm-start.R)); with that block dropped,
  six lines over the four files, they pass, 390 results.
- Size of a state, 4 chains, 75 trees, 200 kept draws (critic; run): positions on the live blocks add
  0.14 percent in memory and 0.08 as an rds; on the kept blocks 28 and 5.
- `FlatNode` does not grow: the position shares the four bytes of `numMaskWords`, a declared member
  only a pooled mask reads; 24 bytes on 64-bit targets and 20 on i386, by static assertion on eight
  targets (critic; read). It is written field by field, never bytewise (read).
- Consumers (read). The flat C header has no entry that builds or takes a grid, a tree or a state
  ([`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h)). stan4bart splices whole chains of its
  samplers' states and names no block; bartCause stores no state; bairrtt sets strictly increasing
  grids.

## The rule

1. A flat tree's threshold split records its position on the grid it was flattened against, beside
   its value. A leaf and a categorical split record none.
2. One resolver turns a value into a position, for the build and the move alike. It is given the
   grid, the value, the position the split had (a flat node's record; before a move, the rule's own;
   possibly none), the node's interval and its missing direction, and returns:
   - the position it had, when that lies inside the interval and the grid holds the value there;
   - else, of the positions inside the interval that hold the value, the first when missing values go
     right and the last when they go left, saying so when there was more than one;
   - else none.
3. With none: the move does what it does today (the nearest value inside the interval; the subtree
   merged when the interval is empty). The build takes the first position holding the value, as
   today, marks the tree, and once the tree is partitioned the marked split is merged by the pass
   that merges an empty side and a split past a shortened grid. A value the grid does not hold fails
   the build, as today.
4. `setState` returns `FALSE` when, in any live tree, the resolver chose among several positions or a
   split was merged under 3. A recorded position that no longer fits and has one position to go to is
   the stored tree and leaves it `TRUE`.
5. A state carries positions beside its live tree blocks always, and beside its kept blocks when some
   threshold column's grid in it repeats a value.
6. A positions block that is not integer, is not one entry per node, holds an `NA`, a negative number
   or one past the largest count a grid can have, or is not 0 on a leaf or a categorical split, is
   refused; so is a live block with an entry past the grid its own state carries. A well-formed
   position is never trusted: rule 2 checks it.

On a grid without repeats one position holds each value, so rule 2 returns it by its first or second
clause, which is what `buildFromFlat` and the move's exact match return today; the move's inexact arm
and the failed build are today's code. That is the whole of the neutrality argument, and the gates
below test it.

## Constraints

- A fit whose grids repeat no value draws, stores, restores and reports exactly what it does now.
- No grid changes: no builder, `setCutPoints` (dec-B285 as landed) or state check is touched. The
  sampler's own grids keep their repeated points; whether they should is the maintainer's question,
  put by the coordinator. If the answer is one point per value, steps are added (What it leaves) and
  none here is removed.
- The tree prior, the moves, `getTrees`, `printTrees` and `predict` are not touched; a kept tree
  replays by value.
- State format: four optional blocks are added; [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)
  and the readable floor stay at 1 by the registry rule. A state without the blocks installs.
- [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move and no facade virtual
  changes; `--preclean` all the same (tree.hpp is a header every object reads).
- `FlatNode` keeps its size and every offset; a static assertion says so.
- Only the live build sets the verdict: the checks on scratch trees do not, and a refused install
  leaves the sampler, its stored state and its generators as they were.
- No warm start or `installTrees` gains a return value.
- No test compares a restored sampler with an untouched twin under a fixed tolerance.
- No R code changes. No NEWS entry (see NEWS).

## Steps

1. The position on a flat node. [`FlatNode`](../../src/bartcore/tree.hpp): a union over
   `numMaskWords` names the position, one more than the index, 0 for none, with a static assertion on
   the node's size. [`Tree::flatten`](../../src/bartcore/tree.hpp) writes it for a threshold split.
   The comment above `FlatNode` stops saying a value names the first index holding it.
2. The resolver and its two callers, in [`Tree`](../../src/bartcore/tree.hpp). A static routine as
   rule 2. [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp)'s threshold arm reads the node's
   interval ([`Tree::splitInterval`](../../src/bartcore/tree.hpp); the ancestors are built by then)
   and calls it with the flat node's position; its out-flag, the one a dropped missing direction
   sets, is also set when the routine chose. With none it places and marks as rule 3 says.
   [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)'s walk calls it with the rule's own
   position before its nearest-value search and keeps that search for none.
   [`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp)'s test gains "a marked tree's split
   outside its interval" beside [`Tree::ruleIsUnrepresentable`](../../src/bartcore/tree.hpp) and
   clears the mark. [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) and
   [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp) run that pass when the tree is
   marked as well as when a bottom is empty; one condition each, the flag already set there.
3. The bridge. [`storeFlatTrees`](../../src/R_interface_bartcore.cpp) takes a slot for positions and
   [`storeState`](../../src/R_interface_bartcore.cpp) appends four names to its two slot lists:
   `tree.positions` and `variance.positions`, always written; `saved.positions` and
   `variance.saved.positions`, written when a threshold column's grid in the state has equal
   neighbours. [`readFlatTrees`](../../src/R_interface_bartcore.cpp) takes the block and, for a live
   block, the state's own grid, refuses by rule 6 with `malformed tree positions in bartcore state`,
   and fills the nodes; its eight call sites in [`setState`](../../src/R_interface_bartcore.cpp) and
   [`readWarmStartState`](../../src/R_interface_bartcore.cpp) pass the block of their name. Both
   functions read the grid before the trees today; the live check depends on that order.
4. tests/cpp. In [test_tree.cpp](../../tests/cpp/test_tree.cpp), beside
   [`testMissingMechanics`](../../tests/cpp/test_tree.cpp):
   - `testResolveCutPosition`: on a hand grid with runs, each clause of rule 2 and each way out of
     it: the position it had taken; refused when outside the interval, when on another value, when
     past the grid; first against last by missing direction; one copy inside (no choice reported);
     none inside; an empty interval; the interval clipped to the grid. Proves the routine alone.
   - `testFlatPositions`: a tree built by hand over a column with missing values, a parent on the
     last copy of a value with missing values left and its left child on an earlier copy with them
     right, a categorical split and a pooled one beside them. `flatten` records each threshold
     split's position and 0 elsewhere; `buildFromFlat` gives every position back, flag clear;
     flattened again it is the same nodes. Fails today.
   - `testFlatPositionsStale`: the same nodes with positions zeroed, shifted, set on another value
     and swapped between parent and child: every split holds its value inside its interval, and the
     flag is set exactly when two copies were inside. A value absent from the grid still fails.
   In [test_state.cpp](../../tests/cpp/test_state.cpp), beside
   [`testDegenerateGridRestores`](../../tests/cpp/test_state.cpp),
   [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp) and
   [`testRestoreStatus`](../../tests/cpp/test_state.cpp):
   - `testRepeatedGridRestores`: a sampler of 10 trees and 2 chains over a constant column with
     missing values, a narrow one, one whose grid a state tripled, a factor and an ordered factor,
     run until splits sit on later copies. `getState`, `setState` into a second sampler and into
     itself: every tree's positions equal in every chain, the verdict exact, the state read back
     equal ([`statesAgree`](../../tests/cpp/common.cpp), which compares the shared field), and the
     two restored samplers' next draws bit-identical. The same with a variance forest. Fails today.
   - `testRepeatedGridVerdict`: that state with positions zeroed installs, every split inside its
     interval, verdict inexact; the same on a grid without repeats is exact.
   - `testSplitWithoutRoom`: a state whose tree stacks two splits on one value, onto a grid holding
     that value once: the child is merged, no split outside its interval, verdict inexact, fits
     finite; the mean forest and the variance forest. The existing checks that a hand-built
     redundant split merges on a grid without repeats stay as they are.
   - `testMoveKeepsPositions`: `setData` with the same data keeps every position on a repeated grid
     (fails today), and with another column's values changed keeps the repeated column's; a warm
     start from a donor whose grid differs on one column keeps the other columns' positions (fails
     today); a split whose value the new grid does not hold still goes to the nearest inside its
     interval; the variance forest through both.
   - `testWarmStartFromKeptDraws`: kept draws on a repeated grid with missing values, a warm start
     from a kept slot ([`testVarianceWarmStartSlot`](../../tests/cpp/test_state.cpp) is the model):
     positions the kept ones on the same grid; after the donor's column is regridded, no split
     outside its interval (fails today).
5. tinytest, a new file `test-repeated-cut-restore.R`. Helpers read the live positions from the stored
   state's `tree.positions` and check each split against its ancestors' interval from `tree.vars`,
   the positions and the sizes.
   - Every way a repeat arises (uniform over a narrow column and quantile over adjacent doubles,
     each with and without missing values in the column; uniform over a constant column with
     missing values; a constant column given values with the grid kept; a grid a state tripled),
     20 trees and 2 chains, through `setState` on itself, `setState` into a second sampler, `copy`
     and a reload: the positions after are the positions stored, `setState`
     returns `TRUE`, the state stored again is identical, and two restored samplers draw
     bit-identical next draws. Fails today on every row.
   - A grid without repeats through the same doors: the same checks, and the stored `tree.positions`
     name the one holder of each value. Holds today but for the block.
   - 25 rounds of sweeps and restores beside missing values: no split outside its interval. Fails
     today (22 of 200).
   - `setData` with the same data keeps every position on a repeated grid (fails today) and on one
     without. `setPredictor` on the repeated column with the same values keeps every position with
     and without `updateCutPoints`; with new values and `forceUpdate = TRUE` no split lies outside
     its interval.
   - Warm starts, `installTrees` and `bart`'s `warm.start`: from live trees on the same grid the
     positions are the donor's; from kept draws likewise, and after the donor's regrid no split lies
     outside its interval (fails today, 9 of 96); across grids that differ on another column the
     repeated column's positions are the donor's (fails today).
   - Kept draws: `saved.positions` is written on a repeated grid and absent without; `predict` from
     kept trees is identical before and after `copy` and a reload on both.
   - A variance forest and two forests on a repeated grid: the first check, on every forest's block.
   - Factor and ordered-factor columns beside a repeated one: a categorical split's entry is 0 and
     its rule restores; an ordered factor's positions restore.
   - The verdict: positions removed, on a repeated grid `FALSE` and every split inside its interval;
     on a grid without repeats `TRUE`. Positions redrawn among the copies: installs, no split outside.
     A state stored on a grid with runs installed on the same grid with each value once: `FALSE`, no
     split outside, draws finite.
   - Refusals, each with the message and the sampler's next draws a twin's: a block one short, a
     double, an `NA`, a negative entry, an entry past the grid, a nonzero entry on a leaf; for the
     live, kept and variance blocks, through `setState` and through `installTrees`.
   The four files under Context drop the positions block where they rewrite a tree block by hand.
   [test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R) is not touched.
6. Mutations (Verification): apply, install with `--preclean`, run, report the counts, revert, `touch`.
7. Records.
   - Manual, [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd), in the value: after the
     dropped side for missing values, "or a split on a cut point its column's grid holds more than
     once had to be placed without its recorded position, as for a state stored by an earlier build
     or one whose positions do not fit its grid". In the paragraph on restoring, one sentence: a
     state that rewrites a tree by hand must drop the positions stored beside it.
   - [Tree storage forms](../architecture.md#tree-storage-forms): the flat form carries a threshold
     split's position, and the resolver's two callers are named.
     [What a state carries](../design/state-not-model.md#what-a-state-carries): one row, the
     positions: state, installed with the trees and checked against the grid.
   - TODO: `repeated-cut-restore` names this plan; `cut-point-weights` gains the paragraph under
     What it leaves.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; `tests/cpp` passes, clean under ASan and
  UBSan; the full tinytest suite; the new test file and the four edited ones under ASan on the
  R-loaded path (the reader indexes a caller's block: length before index, position against the count
  before the grid is read).
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three
  compares are bitwise with `EQUIVALENCE_CORES=2`, every scenario reporting identical draws, counted
  per scenario with no `max |z|` line: 55, 15 and 11 at the last landing's counts. Nothing is
  re-recorded.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged: the state gains
  blocks, which is a change to what a fit carries; none of the gates restores a state (by search).
- One script on the base and slice builds digesting, on grids without repeats, the draws of a run, a
  restore into a second sampler with the value `setState` returns, a copy, a reload, a `setData`, a
  `setPredictor` with and without `updateCutPoints`, a warm start from live and from kept trees on
  the same and on another grid, for the kinds of sampler under Context, and the stored states less
  the new blocks. Equal.
- The same kinds on a repeated grid through every door, each in a process of its own: no tree moved,
  no split outside its interval, fits finite.
- Mutations, each expected to fail the named test:
  - `flatten` records no position: `testFlatPositions`, tinytest's first check;
  - `buildFromFlat` ignores the position: the same two;
  - the position is taken without the value check: `testFlatPositionsStale`, tinytest's redrawn and
    regrid checks; without the interval check: `testFlatPositionsStale`'s swapped pair;
  - the fallback ignores the interval, or always takes the first copy: `testFlatPositionsStale`,
    tinytest's 25 rounds with positions removed;
  - a split without room is left: `testSplitWithoutRoom`, tinytest's grid with each value once;
  - the choice does not report: `testRepeatedGridVerdict`, tinytest's verdict; it reports on a grid
    without repeats: the same;
  - the move does not ask the resolver: `testMoveKeepsPositions`, tinytest's `setData` and
    cross-grid checks;
  - the writer drops a live block, the variance block, or the kept block on a repeated grid; writes
    the kept block on every grid; a reader call site passes no block or another's: tinytest's block
    checks and `testWarmStartFromKeptDraws`;
  - each refusal of rule 6 removed in turn: tinytest's refusals.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; `R CMD build`
  with every vignette rebuilt and `R CMD check --as-cran` on the tarball from a clean copy (man/
  changes).
- Not a hot-path change: nothing is added to a sweep. A build from flat nodes and a move walk each
  node's ancestors once more.

## What it leaves, and where it goes

- TODO `cut-point-weights`, one paragraph added: "The restore is ready for a grid with repeats
  (docs/plans/repeated-cut-restore.md). Three things are not, measured on a grid with every point
  tripled: `setPredictor(updateCutPoints = TRUE)` rebuilds the grid without the repeats (60 points,
  60 distinct) and `setData` derives `n.cuts` points (20), so both drop the weights in silence; the
  undo through `setCutPoints` takes a repeated grid only while the column still holds it; and a
  point with 5 copies is drawn at the root 4.1 times as often as a single one, not 5, the empty-leaf
  rule taking the rest, so a repeat is not an exact weight." Nothing built here reads a repeat as
  more than one more position, so neither form of weights is prejudged.
- If the maintainer rules one point per value on the grids the sampler builds: added, nothing above
  removed. The two fill routines drop repeats; a regrid that cannot keep its count needs a rule; a
  column created constant and then refreshed needs one (it would keep one point where today and in
  0.9-34 it gets `n.cuts`); 8 checks of test-cut-points-undo.R and 2 of
  [`testDegenerateGridRestores`](../../tests/cpp/test_state.cpp) are rewritten; the class becomes
  posterior-changing for fits with a narrow column or a constant column with missing values, which
  no baseline holds. Positions stay: a state may still bring a repeated grid, and a caller's repeats
  come with weights.
- TODO `state-frame-prior`: when a state stops installing its grid, the build on the state's grid and
  the move serve as they do for a warm start today, and the move now keeps positions.
- Not this plan's: `copy` takes the last stored state, so with `updateState = FALSE` a copy made
  after `setCutPoints` holds the old grid (critic; read). A gp sampler's restore is not a
  continuation on any grid. Bit for bit against the sampler that never stopped is dec-B32's.

## What waits on what

- [leaf-conversions.md](leaf-conversions.md) first. It edits `Chain::applyNewData`,
  `Chain::convertStateUnits`, `Chain::rebuildLiveForestRemapped`, `Chain::getState`, the leaf models,
  `Sampler::installForests`, and in the bridge the leaf calibration blocks of the state's writer and
  both readers. This plan edits tree.hpp, which it does not touch; one condition each in
  `Chain::rebuildLiveForest` and `Chain::rebuildVarianceForest`; and in the bridge the tree blocks of
  the same three functions. The two meet in those three bridge functions, at different blocks, in
  `Chain::rebuildLiveForestRemapped`, which this plan reads and does not edit (it builds and then
  moves, and both halves change underneath it), and in tests/cpp's test_state.cpp and its list of
  tests. Branch after it lands; no shared line is expected.
- dec-B285's sorting and refusal in `setCutPoints` have landed and are not touched.
- The maintainer's answer on the sampler's own grids does not block this plan (Constraints).

## stan4bart

No edit. The flat C header has no entry for a grid, a tree or a state and does not change. stan4bart
stores its samplers' states whole on its fit and installs them with `setState`, ignoring the value: a
fit saved before this change restores as it does today, and one saved after carries the new blocks
inside the chains it already splices whole.

## NEWS

None owed. The fault is a regression against 0.9-34, which stored a live tree by position, but it was
introduced and is removed on this branch and never reached a release; NEWS compares releases. The
value `setState` returns is new in 1.0-0 and its entry, if it has one, lists no causes.

## Calls made in planning

The coordinator's eight, on the critique of the design (agent-made, reversible):

1. Where a recorded position does not hold its value, the fallback is not the first holder but a
   position inside the node's interval, the first when missing values go right and the last when
   they go left. 9 of 96 warm starts from stale kept draws left a split outside its interval under
   the first-holder fallback, 0 of 96 under this.
2. The move by value is brought under the same resolver as the build, and the claim is the Goal's
   sentence: every split on a position that holds its value inside its interval, and the same
   position where the grid there is unchanged.
3. A state without positions on a repeated grid does not report itself exact: the engine's flag, as
   the other inexact installs use it, and the manual says what a caller sees.
4. The claim that weights by repeats is then one refusal lifted is withdrawn; what else would change
   is the paragraph for TODO `cut-point-weights`.
5. The refusal of a malformed block is built and tested, and a well-formed block is checked against
   the values: a position whose grid value is not the stored value is stale, never trusted.
6. The bar: two restores of one state give bit-identical next draws and trees identical position
   for position; against the uninterrupted sampler, what the manual and dec-B32 promise. No test
   asserts a tolerance single storage or gp leaves cannot meet.
7. Positions on the live blocks always and on the kept blocks where a grid repeats (the planner's
   answer to the coordinator's question; reasons below).
8. Planned for the sampler's own grids keeping their repeated points; the question is the
   maintainer's and the other answer adds steps.

The planner's:

- A position is checked against the interval as well as the value. The critic's 40 states with
  positions redrawn among copies show why: each held its value, and 3 lay outside.
- A split the grid cannot hold inside its interval is merged. Call 1's rule alone has nowhere to put
  it: 32 of 200 installs of a state stored on a grid with runs, onto that grid with each value once,
  left a split outside. Refusing instead was prototyped and refuses hand-built states with a
  redundant split that install and merge today (five tests/cpp checks). The merge is the treatment a
  split past a shortened grid already gets, by the same pass.
- Live always, kept on condition. Always-live costs 0.14 percent of a state, makes every restore in
  the suite run the reader, and leaves no branch in the writer for the blocks a continuation depends
  on. Its price is that a state whose tree block is rewritten by hand must drop the block beside it:
  four test files, six lines. The kept blocks would cost 28 percent in memory on every state that
  keeps trees, so they keep the condition; their reader is the routine the live blocks always run.
  Not taken: no kept positions at all, two blocks fewer, because a warm start from kept draws on an
  unchanged grid would then land on another copy than the one drawn, against call 2.
- Not taken: the critic's alternative of a live block holding positions in place of values. Live and
  kept trees would stop sharing one record, and the value is what a stale position is detected by.
- The verdict stays `TRUE` for a stale position with one copy to go to: the tree is the stored one.
- A kept block is held to the largest count a grid can have, not to the grid in force: kept draws
  recorded before a regrid hold positions past a shorter grid honestly.
- A live entry past its own state's grid is refused in the bridge, where both are at hand, so the
  engine keeps one rule for a position that does not fit.
- The tip against the design. The design's first text called its suite run full at 221 of 232 files;
  the figures above are from the whole suite. Its "fixes every grid" and "one refusal lifted" were
  wrong and are corrected in it. Nothing in this slice was found done already.
