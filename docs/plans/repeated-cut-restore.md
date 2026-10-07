# repeated-cut-restore: no cut grid repeats a point

Status: PLANNED (dec-B285, dec-B297, dec-B298, dec-B299, dec-B300). Rewritten 2026-10-07 after the four
rulings; the plan this replaces kept repeated points and stored each split's position, and none of that
is built. The file keeps its name because the TODO item, four register entries and
[cut-points-undo.md](cut-points-undo.md) cite it.

agent: opus implementer, one (engine, bridge, tests, gate arms); opus reviewer, told to refute.
rng: by fit and by call.
- POSTERIOR-CHANGING at creation and `setData`, for a sampler with a numeric column whose grid today
  holds a point twice and on which a split can survive: under the uniform rule a column narrower than
  its `n.cuts` evenly spaced points resolve, and a column with one finite value beside missing values
  or `Inf`; under the quantile rule a column where two of the chosen midpoints round to one double.
  The tree prior over that column's thresholds changes.
- POSTERIOR-CHANGING after `setPredictor(updateCutPoints = TRUE)` wherever the refreshed grid differs
  from today's: the count a column held is no longer kept (a column created with few distinct values
  under the quantile rule; a grid set by `setCutPoints` or brought by a state at another length), a
  refresh onto too few distinct values shrinks where it failed or kept the old grid, and a refreshed
  grid drops repeats.
- NEUTRAL, bit for bit, for every other fit and call. That includes a column of few distinct values
  that are far apart (a 0/1 column under the uniform rule keeps its `n.cuts` distinct points), a
  constant column under the quantile rule (one point already), and a constant column with no missing
  value under the uniform rule, whose grid goes from `n.cuts` copies to one point and whose draws do
  not move (run, Context).
- Not a generator matter but a changed result: a state or a warm-start donor whose grid repeats a
  point is refused; a stored split outside its node's interval is merged and `setState` returns
  `FALSE`; an unforced refresh that would strand a split past a shorter grid returns `FALSE`.
Proved by the bitwise gates below (no baseline scenario and no snapshot file builds, refreshes or
installs a repeating grid), by two exact arms that fail on the base build, and by a seeded digest on
the base and slice builds.
window: before the merge to main (VD 2026-10-07, dec-B285). After [leaf-conversions.md](leaf-conversions.md)
lands; see What waits on what. Step 1 may land as a push of its own.
budget: ~1400 lines (data.hpp ~90, tree.hpp ~45, chain.hpp ~20, sampler.hpp ~25, bridge ~60, R ~5,
tests/cpp ~340, tinytest ~470, the two gate arms and the workflow ~180, design note, manual, NEWS,
comments and TODO ~165), upper figure 2600. It rests on the prototype's deduplication (about 45 engine
lines) and build-side interval check (about 80), and on [leaf-conversions.md](leaf-conversions.md),
planned at 1270 and built at about 2900, most of the excess in tests.

## Goal

Every cut grid a sampler builds, refreshes or accepts holds each point once, so a stored split's value
names one position and a restore by value is exact with nothing more stored. A refreshed grid is the
grid a creation would give for the new values. A state whose grid repeats a point is refused by name.
A split built from a stored tree lies inside the interval its ancestors leave, or is merged.

## Context

At a355852d; the engine and bridge files are those of [leaf-conversions.md](leaf-conversions.md)'s
branch, whose base differs from the tip in documents only and which does not touch the cut code read
here. "Run" is a run of mine on the user library (the tip, shipped build) unless it says "prototype":
the design's prototype of cc8381b3, a reference build with a switch that drops repeated points at
derivation and one that reports every repeating grid; benchmarks/ has not changed since cc8381b3.
"Read" is the source or another agent's recorded run.

- Where a grid repeats today (run; `n.cuts = 100`). Uniform rule
  ([`fillCutsOverRange`](../../src/bartcore/data.hpp)): five adjacent doubles give 100 points, 5
  distinct; a constant column 100, 1, with or without missing values; a 0/1 column 100, 100. Quantile
  rule ([`fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp)): five adjacent doubles give 4, 3; a
  constant column 1; a 0/1 column 1; six distinct values 5, 5.
- A refresh today (run; `n.cuts = 20`; [`refreshCutsForColumn`](../../src/bartcore/data.hpp) counts
  from the count the column holds). Uniform: a placeholder column of zeros is created with 20 points, 1
  distinct, keeps them when `updateCutPoints` is left `FALSE`, and gets 20 distinct once refreshed
  onto 200 values; refreshed onto five adjacent doubles, 20 points, 5 distinct; onto a constant the
  call returns `TRUE` and the old grid stays, in silence. Quantile: the placeholder has 1 point and
  keeps 1 after a refresh onto 200 values; a six-valued column keeps 5; a refresh onto fewer distinct
  values stops with "number of induced cut points in new predictor less than previous". Either rule:
  after `setCutPoints` with 60 points a refresh gives 60. `setData` derives afresh (20 points, 5
  distinct, onto the adjacent doubles; 5 onto six values under the quantile rule).
- A refresh keeps each split's position and maps no value (read; run for a monotone map under the
  quantile rule, 25 of 25 positions). Nothing restores a column's count when an unforced update rolls
  back ([`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp), [`SubsetUpdate`](../../src/bartcore/sampler.hpp)
  put back the points alone), and the unforced check asks only that every leaf is occupied
  ([`Chain::revalidateTrees`](../../src/bartcore/chain.hpp)): both lean on the count never changing.
  A forced update and `setCutPoints` already merge a split past a shorter grid
  ([`Tree::ruleIsUnrepresentable`](../../src/bartcore/tree.hpp)). No scratch buffer is sized from a
  count ahead of use (read).
- What is accepted today (run). `setCutPoints(c(0, 0.5, 0.5, 1))` is refused: "a cut point may appear
  only once in 'cuts'"; a grid out of order is sorted; the repeating grid a column holds is taken back.
  A state with a column's grid tripled installs and `setState` returns `TRUE`; `installTrees` takes
  such a state as a donor; a sampler on a narrow column copies and reloads with its 100 points, 5
  distinct.
- Which fits move (prototype, seeded, 20 trees, 2 chains, with and without the drop). Move: the narrow
  column under both rules, the constant column beside missing values under the uniform rule.
  Identical draws: a constant column with no missing value (100 points to 1), a 0/1 and a six-valued
  column, a column with missing values, and the constant column under the quantile rule. The first
  is identical because no split on such a column is ever accepted and the column stays available at
  every node either way; run under the default move mixture only.
- Where a grid is made or taken, recounted: thirteen, and two prechecks that make none. Derived (6):
  [`buildCutsForColumn`](../../src/bartcore/data.hpp)
  (creation, `setData`) through the two fill routines above, the first serving the dense and the
  sparse scan, and
  [`fillCutsAtLevelMidpoints`](../../src/bartcore/data.hpp), which cannot repeat;
  [`refreshCutsForColumn`](../../src/bartcore/data.hpp) and
  [`refreshCutsForCscColumn`](../../src/bartcore/data.hpp) with their prechecks
  [`cutsWouldRemainValid`](../../src/bartcore/data.hpp) and
  [`cutsWouldRemainValidCsc`](../../src/bartcore/data.hpp). Taken (7):
  [`setCutPointsForColumn`](../../src/bartcore/data.hpp) from
  [`Sampler::setCutPoints`](../../src/bartcore/sampler.hpp) and
  [`Sampler::setState`](../../src/bartcore/sampler.hpp);
  [`ScopedCutGrid`](../../src/bartcore/data.hpp), a donor's grid for the length of a warm start
  ([`Sampler::installForests`](../../src/bartcore/sampler.hpp)); a view's copy of its parent's
  ([`buildFromParent`](../../src/bartcore/data.hpp)); in the bridge
  [`bartcore_setCutPoints`](../../src/R_interface_bartcore.cpp) and the `cutPoints` attribute read
  by the state reader (["malformed cut points in bartcore state"](../../src/R_interface_bartcore.cpp))
  and by [`readWarmStartState`](../../src/R_interface_bartcore.cpp); in R
  [`bartcoreSamplerSetCutPoints`](../../R/bartcore.R). Both state checks use
  [`cutGridIsValid`](../../src/bartcore/data.hpp) in its form that allows equal neighbours. Edited:
  the two fill routines, the two refreshes and their prechecks, the validity check and its three
  callers, and the two bridge readers; the rest hold the rule because what feeds them does.
- Gates (run on the prototype with the report on, `quick`): of the gates
  `.github/workflows/exact-gates.yaml` lists, all but `bcf-latent-exact.R` and
  `monotone-exact-enumeration.R` were run (those two: read, their cut count is one less than their
  cell count). None refreshes a grid. Two build a repeating one, `heteroscedastic-exact.R` and
  `multinomial-exact.R`, each over a constant column with no missing value, and print identical
  output with and without the drop. So no exact gate covers a moved fit. The three equivalence
  harnesses (81 scenarios) and the four snapshot files build no repeating grid (read: the design's
  run with the report on) and refresh none (run: by search, the one `updateCutPoints` in them is
  `FALSE`).
- Tests that pin what changes (read, by search). The refusal of a coarser quantile refresh:
  ["induced cut points in new predictor less than previous"](../../inst/tinytest/test-quantile-grid.R),
  ["induced cut points"](../../inst/tinytest/test-bartcore.R),
  ["induced cut points"](../../inst/tinytest/test-mutate-sparse-valued.R),
  [`testQuantilePredictorUpdate`](../../tests/cpp/test_moves.cpp). The held count:
  ["a refresh spreads the count the column holds"](../../inst/tinytest/test-quantile-grid.R). The old
  grid kept over a constant column: [`testDegenerateReCutRoundTrips`](../../tests/cpp/test_moves.cpp).
  Repeating grids: [test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R) (8 checks
  failed under the prototype's drop) and [`testDegenerateGridRestores`](../../tests/cpp/test_state.cpp)
  (2). Eight tinytest files name `updateCutPoints`, seven with `TRUE`.
- A stored tree that stacks two splits on one value installs today on a grid holding the value once,
  and beside missing values, where both sides of the lower split stay occupied, it is left as it is:
  32 of 200 such installs kept a split outside its node's interval (read: the design's run).
- A state written by an earlier build of this branch carries the same format number and installs,
  unless a grid in it repeats: then it is refused, and that is every such state over a constant
  numeric column under the uniform rule. No build of this branch was released. A 0.9-34 state is
  another format: a state with no format number is refused at the floor
  ([`minReadableStateFormatVersion`](../../src/R_interface_bartcore.cpp); run with the number
  removed, "state encoding version 0 ... re-fit with this version"; a real 0.9-34 object not run).
- 0.9-34 (read on main, not run): the uniform rule placed `n.cuts` points whatever the range, so a
  narrow column repeated points there too; the quantile rule gave one point fewer than the distinct
  values; a quantile refresh onto fewer stopped with the message above and onto more warned and kept
  the count.

## The rule

1. dec-B297, the maintainer on keeping repeats: "I don't see the motivation for A. I doubt it matches
   a user's intent." The entry's rule: "a grid is a set of distinct thresholds, and every grid the
   sampler builds or accepts holds each value once." Testable: after any call, no column's grid has
   two equal neighbours.

   ```r
   ctl <- dbartsControl(n.cuts = 100L, n.chains = 1L)
   narrow <- 1 + .Machine$double.eps * (seq_len(n) %% 5L)
   s <- dbarts(cbind(z, narrow), y, control = ctl)
   length(attr(s$state, "cutPoints")[[2L]])    # 5; today 100, 5 distinct
   const <- replace(rep(3, n), 1:40, NA)
   s <- dbarts(cbind(z, const), y, control = ctl)
   attr(s$state, "cutPoints")[[2L]]            # 3; today 100 copies of 3
   ```

2. dec-B298: "Up to n.cuts, from the new values." The entry's rule: "a refresh of a column's grid
   derives it as a creation would for those values, up to the sampler's n.cuts, and the number of
   points a column held before does not enter." Testable: after `setPredictor(x, j, updateCutPoints =
   TRUE)` column j's grid is identical to the grid of a sampler created on those values with the
   same `n.cuts` and rule; `setData` already is.

   ```r
   s <- dbarts(cbind(z, latent = 0), y, control = dbartsControl(n.cuts = 20L, useQuantiles = TRUE))
   s$setPredictor(rnorm(n), 2L, forceUpdate = TRUE, updateCutPoints = TRUE)
   length(attr(s$state, "cutPoints")[[2L]])    # 20; today 1
   s$setCutPoints(seq(-2, 2, length.out = 60L), 2L)
   s$setData(dbartsData(cbind(z, latent = rnorm(n)), y))
   length(attr(s$state, "cutPoints")[[2L]])    # 20, as today
   ```

3. dec-B299: "Shrink: the grid gets the 5 distinct points." The entry's rule: "a refresh never fails
   for want of distinct points and never keeps the old grid; the grid is the distinct points the
   values supply, up to n.cuts, under the uniform rule and the quantile rule alike." Testable: an
   accepted refresh onto values that supply fewer points leaves that shorter grid, under both rules,
   with no error.

   ```r
   s$setPredictor(narrow, 2L, forceUpdate = TRUE, updateCutPoints = TRUE)
   length(attr(s$state, "cutPoints")[[2L]])    # 3 (quantile rule); today an error
   ```

4. dec-B300: "Weights later come as probabilities beside the grid, never as repeats" and "The state
   question answers itself: refuse". The entry: such a state "is refused by setState, a copy and a
   reload, naming the repeat." Testable: `setCutPoints` given a repeat stops as today; a state whose
   grid repeats a point stops, naming the column and the value, and leaves the sampler as it was.

   ```r
   s$setCutPoints(c(0, 0.5, 0.5, 1), 2L)       # error, as today
   st <- s$state
   attr(st, "cutPoints")[[2L]] <- rep(attr(st, "cutPoints")[[2L]], each = 3L)
   s$setState(st)                              # error naming column 2; today TRUE
   ```

5. Kept from the plan this replaces (dec-B297: "the part of the fix that keeps a split inside its
   node's interval after a grid changes stands"). A split built from a stored tree sits on the one
   position holding its value. Where that position is outside the interval the node's ancestors
   leave, the split is merged once the tree is partitioned, by the pass that merges an empty side,
   and the install reports itself altered. Only a tree a sampler did not write can hold one: a
   hand-edited state, or draws kept by an earlier build on a grid that repeated. The move by value
   at `setData` and a warm start from another grid already clamps to the interval and merges an
   empty one ([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)); it is tested here, not
   edited.

What the rulings leave open is under Calls made in planning; the two that shape the build are how a
rule's points are counted (call 2) and what a refresh that changes a column's count does to the splits
on it (call 3).

## Constraints

- A sampler none of whose grids repeats a point today, and which never refreshes one, draws, stores,
  restores and reports exactly what it does now.
- The tree prior, the moves, `getTrees`, `printTrees` and `predict` are not touched. A constant column
  keeps one point and stays available, as it is under the quantile rule today.
- Cut points are derived over every row in the design, a row at weight 0 or masked out included
  (dec-B296 leaves them out of its rule).
- State format: no block added, renamed or changed; [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)
  and the floor stay at 1.
- [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move; no facade virtual changes;
  `--preclean` on every engine commit all the same.
- A refusal leaves the sampler, its stored state and its generators as they were.
- The line the mutation battery anchors in [`runPredictorTransaction`](../../src/bartcore/sampler.hpp)
  stays as it is.
- No R argument is added. `setCutPoints` keeps dec-B285's sort and refusal.

## Steps

1. The interval at a build (NEUTRAL; may be pushed alone).
   [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp)'s threshold arm compares the position found
   with [`Tree::splitInterval`](../../src/bartcore/tree.hpp) (the ancestors are built by then) and
   marks the tree when it lies outside. [`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp)
   merges a marked tree's split outside its interval beside
   [`Tree::ruleIsUnrepresentable`](../../src/bartcore/tree.hpp), and clears the mark.
   [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) and
   [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp) run that pass for a marked tree as
   for an empty bottom and set their altered flag. The checks on scratch trees set no verdict.
2. One point per value where a grid is derived. A routine drops equal neighbours and sets the count;
   [`fillCutsOverRange`](../../src/bartcore/data.hpp) and
   [`fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp) end with it. The comments on
   [`ColumnStore`](../../src/bartcore/data.hpp) and [`requestedNumCuts`](../../src/bartcore/data.hpp)
   stop saying a count is fixed once built.
3. The refresh. [`refreshCutsForColumn`](../../src/bartcore/data.hpp) and
   [`refreshCutsForCscColumn`](../../src/bartcore/data.hpp) derive a numeric column's grid by the
   numeric arm of [`buildCutsForColumn`](../../src/bartcore/data.hpp), counting from
   [`requestedNumCuts`](../../src/bartcore/data.hpp); neither refuses.
   [`cutsWouldRemainValid`](../../src/bartcore/data.hpp) and
   [`cutsWouldRemainValidCsc`](../../src/bartcore/data.hpp) keep their factor arm and pass every
   numeric column; [`valuesAreDegenerate`](../../src/bartcore/data.hpp) and
   [`cscColumnIsDegenerate`](../../src/bartcore/data.hpp) go if nothing else reads them.
   [`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp) and
   [`SubsetUpdate`](../../src/bartcore/sampler.hpp) snapshot and put back the counts with the points.
   [`Chain::revalidateTrees`](../../src/bartcore/chain.hpp) and
   [`Chain::revalidateVarianceTrees`](../../src/bartcore/chain.hpp) also fail on a tree holding a rule
   past its column's grid, so an unforced update that would strand one rolls back; the forced path
   merges it as it does after `setCutPoints`
   ([`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp), unchanged). The message of
   [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp) and
   [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp) for
   [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp)'s refusal, which only a factor value
   off its level table can now raise, says that.
4. What is accepted. [`cutGridIsValid`](../../src/bartcore/data.hpp) loses its form that allows
   equal neighbours; [`Sampler::setState`](../../src/bartcore/sampler.hpp) and
   [`Sampler::installForests`](../../src/bartcore/sampler.hpp) call the strict one. The two bridge
   readers check each numeric column's grid where they read it and stop with
   `cut points of column <j> in bartcore state repeat a value (<v>): a cut grid holds each point once`
   (the donor reader in its own wording), before anything is installed.
   [`bartcore_setCutPoints`](../../src/R_interface_bartcore.cpp) drops the exception for the grid the
   column holds, which can no longer repeat; a held grid still passes as any distinct grid does.
5. tests/cpp. In [test_tree.cpp](../../tests/cpp/test_tree.cpp): `testBuildOutsideInterval`, a tree
   by hand over a column with missing values, a parent and its left child on one value with the
   directions that leave both sides of the child occupied: the build marks it, the merge removes the
   child, flattened again no split lies outside; the same nodes with distinct values in order build
   unmarked. In [test_data.cpp](../../tests/cpp/test_data.cpp), beside
   [`testRequestedCutCount`](../../tests/cpp/test_data.cpp) and
   [`testQuantileCutPoints`](../../tests/cpp/test_data.cpp): `testDistinctCutGrids`, every derivation
   (uniform and quantile, dense and CSC, creation, `setData`, refresh) over a narrow, a constant, a
   constant-with-missing, an all-missing, a 0/1 and an ordinary column: strictly increasing, the
   counts of Context, and the ordinary and 0/1 grids bit for bit today's; `testRefreshDerivesAsCreation`,
   a refreshed grid equals a fresh store's on the same values, from a shorter, a longer and a set
   grid. In [test_moves.cpp](../../tests/cpp/test_moves.cpp): `testRefreshCountChange`, a shrink
   forced (splits past the end merged, fits finite, state round trip) and unforced (rolled back,
   points, counts and codes as before; beside missing values, where both sides of a stranded split
   stay occupied, still rolled back), a growth from one point, the variance forest through both;
   [`testQuantilePredictorUpdate`](../../tests/cpp/test_moves.cpp),
   [`testDegenerateReCutRoundTrips`](../../tests/cpp/test_moves.cpp) and
   [`testSetPredictorTransaction`](../../tests/cpp/test_moves.cpp) rewritten to the new rule. In
   [test_state.cpp](../../tests/cpp/test_state.cpp): `testRepeatedGridRefused`, a state and a donor
   with one column's grid doubled are refused and the sampler's next draws are a twin's;
   `testStackedSplitsMerge`, a state whose tree stacks two splits on one value installs, the child
   merged, verdict altered, mean and variance forest;
   [`testDegenerateGridRestores`](../../tests/cpp/test_state.cpp) rewritten: the constant column's
   grid is one point and restores. [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp) and
   [`testMapOldCutPointsOntoNew`](../../tests/cpp/test_data.cpp) gain a recipient grid shorter than
   the donor's: no split outside its interval.
6. tinytest, a new file `test-cut-grid-distinct.R`.
   - The four examples under The rule, each rule, with the counts stated there.
   - No grid repeats after creation, `setData`, a refresh, `setCutPoints`, `setState`, `copy` and a
     reload, over the columns of step 5, dense and sparse.
   - A refresh equals a creation's grid: placeholder then values; six values then 200; 200 then six;
     after `setCutPoints` at 60 points; under both rules. Without `updateCutPoints` the grid is kept.
   - A shrink: forced, the sampler runs on and its state restores into a twin with `TRUE`; unforced,
     `FALSE` with data, grid and next draws a twin's. A refresh onto a constant leaves one point.
   - Restores are exact: on a narrow column, with and without missing values, and a constant column
     beside missing values, through `setState`, `copy` and a reload, the printed trees are the stored
     ones, `setState` returns `TRUE`, the state stored again is identical and two restored samplers
     draw identical next draws. Fails today on every row.
   - Refusals, each by its message with the sampler's next draws a twin's: a state with a doubled
     grid through `setState`, a copied sampler's state edited and installed, a reloaded object whose
     state was edited, `installTrees` and `bart(warm.start = )` from such a state; `setCutPoints`
     with a repeat, and with the held grid (taken).
   - A hand-edited state whose tree stacks two splits on one value: installs, `FALSE`, no split
     outside its interval, draws finite.
   Rewritten: the 8 checks of [test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R)
   and the checks named under Context; each of the eight files that name `updateCutPoints` is run
   and any other pin reported.
7. Exact arms, in [bd-balance.R](../../benchmarks/R/bd-balance.R), each its own step in
   `.github/workflows/exact-gates.yaml` as `zeroweight` is.
   - `narrow`: the gate's four cells placed at 1 plus 0, 1, 2 and 3 times `.Machine$double.eps`,
     `n.cuts = 5`. Today's grid is 1, 1 + eps and three copies of 1 + 2 eps (run); the rule gives the
     three points that separate the cells, so the target is the gate's own enumeration. The arm first
     asserts the grid.
   - `constmissing`: one column, a single value on some rows and missing on the rest, `n.cuts = 5`.
     Two trees have mass, the stump and present against missing; under the rule both children are
     unsplittable, under today's five copies most are not, which moves the odds by about a third at
     the default tree prior (computed, not run). The arm asserts the one-point grid and compares the
     odds, the missing direction scored as the engine's tree prior scores it.
   Both must fail on the base build and pass on the slice; report both.
8. Mutations (Verification): apply, install with `--preclean`, run, report the counts, revert, `touch`.
9. Records.
   - Manual, [`dbartsSampler$setPredictor`](../../man/dbartsSampler-class.Rd): `updateCutPoints`
     derives the grid a new sampler would have for the values, up to `n.cuts`, fewer where they
     supply fewer; splits keep their positions; one past a shorter grid is merged when forced and
     declines the call otherwise. `cuts`: the exception for the held grid and its sentence on equal
     neighbours go. `setState`: a state whose grid repeats a point is refused; `FALSE` also where a
     split outside its interval was merged. [`dbartsControl`](../../man/dbartsControl.Rd) and
     [`bartBT`](../../man/bartBT.Rd): a column whose range cannot hold `n.cuts` distinct points gets
     fewer; a constant column one. The method string of
     [`setCutPoints`](../../R/dbarts.R), for `check-rc-codoc`.
   - A design note, docs/design/cut-grid.md, with its row in docs/design/INDEX.md: the rule, the
     measurements, which fits move, the calls. [quantile-grid.md](../design/quantile-grid.md) points
     to it where it states the refresh; [ColumnStore](../architecture.md#columnstore) and
     [What a state carries](../design/state-not-model.md#what-a-state-carries) each gain a sentence.
   - NEWS (below). TODO: `repeated-cut-restore` goes; `cut-point-weights` stands as written.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; `tests/cpp` passes, clean under ASan and
  UBSan; the full tinytest suite; the new file and the rewritten ones under ASan on the R-loaded path
  (a rule past a grid is the read this slice could introduce).
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged and the three
  compares are bitwise with `EQUIVALENCE_CORES=2`, counted per scenario with no `max |z|` line: 55
  against `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded and no snapshot is replayed. A
  scenario that is not identical is a finding; if it stands, the file is recorded again on the
  reference build as benchmarks/README.md says, its MANIFEST row naming the two arms as the oracle.
- Every gate `.github/workflows/exact-gates.yaml` lists, `quick`: pass, unchanged;
  `heteroscedastic-exact.R` and `multinomial-exact.R` print what the base build prints. The two new
  arms pass, and fail on the base build.
- One script on the base and slice builds digesting seeded draws and stored states. Equal: ordinary,
  0/1 and six-valued columns under both rules, a constant column with no missing value under each
  move kind in turn, `setData`, a refresh that keeps its count, `setCutPoints`, a copy, a reload, a
  warm start on each grid, each leaf model and family once. Different, and reported: the three
  columns that move. If a move kind draws differently over the constant column, that fit is SHIFTING
  and is reported, not fixed.
- Mutations, each expected to fail the named test:
  - the drop removed from the uniform fill, or from the quantile fill: `testDistinctCutGrids`,
    tinytest's first check, the `narrow` arm;
  - the refresh counts from the held count: `testRefreshDerivesAsCreation`, tinytest's refresh
    checks; it refuses on too few points, or keeps the old grid over a constant: the shrink checks;
  - a rollback puts back the points without the count: `testRefreshCountChange` unforced;
  - the unforced check passes a rule past the grid: `testRefreshCountChange` beside missing values;
  - the state reader, the donor reader, or the engine's check takes a repeat:
    `testRepeatedGridRefused`, tinytest's refusals; the message names no column: the same;
  - `setCutPoints` takes a repeat: its existing refusal test;
  - the build does not compare with the interval; marks and does not merge; merges and does not
    report: `testBuildOutsideInterval`, `testStackedSplitsMerge`, tinytest's hand-edited state.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; inst/NEWS.Rd
  parses with a non-NULL result; `R CMD check --as-cran` on a tarball from a clean copy.
- As two pushes. Step 1 with its tests and mutations: NEUTRAL; the install, tests/cpp with
  sanitizers, the full suite, the reference-build files and compares. The rest: every gate above.
- Not a hot-path change: a derivation makes one more pass over at most `n.cuts` points; no sweep code
  is edited. No bench compare.

## Consumers

Read at each package's working branch; none edited, none run. The flat C header has no entry that
builds or takes a grid ([`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h)).
- stan4bart passes `n.cuts` to the data object and calls neither `setCutPoints` nor a refresh. Its
  fits move only where a tree column is narrow or constant beside missing values. It splices states
  whole: a fit saved on an earlier build of this branch with a constant tree column is refused at
  restore. Its suite and posterior baselines are run against the slice's library by the reviewer.
- bartCause names none of this. treatSens sets `n.cuts` to 100 and replaces test predictors only.
- bairrtt sets a strictly increasing grid of normal quantiles with `setCutPoints` and updates its
  latent column without `updateCutPoints`: unchanged.

## NEWS

Compared with 0.9-34. One item under the user-visible changes: a predictor's cut points are distinct,
so a column whose range cannot hold `n.cuts` distinct points gets fewer where 0.9-34 repeated them, a
constant column gets one, and fits on such a column change; `setPredictor(updateCutPoints = TRUE)`
derives up to `n.cuts` points from the new values whatever the column held, and no longer stops when
they induce fewer. The existing item on quantile cut points, which says a refresh spreads the cut
points a column has, is amended to match.

## Out of scope, and where it goes

- Weighting cut points: TODO `cut-point-weights`, an argument of probabilities beside the grid, after
  1.0-0 (dec-B300). Nothing here reads a weight.
- A split's position stored in a state, and the positions blocks of the plan this replaces: not built.
- Moving splits by value at a refresh, and a kept draw whose cut value its sampler's present grid no
  longer holds (it cannot seed a warm start, today and after): TODO `state-frame-prior`, with the
  install that moves splits by value.
- `copy` taking the last stored state, so that with `updateState = FALSE` a copy made after a grid
  change holds the old grid: not this plan's.

## What waits on what

- [leaf-conversions.md](leaf-conversions.md) first; this plan is written against its branch for
  data.hpp, chain.hpp, sampler.hpp, the bridge, tests/cpp, the manual, NEWS and the design index. It
  edits `Chain::applyNewData`, `Chain::rebuildLiveForestRemapped`, `Sampler::installForests`'
  conversion, the leaf calibration blocks of the bridge's readers and one function at the top of
  data.hpp. This plan edits the cut code of data.hpp, tree.hpp (which it does not touch), one
  condition in each of four other chain functions, the grid check and the two update records in
  sampler.hpp, and the `cutPoints` read in the same two bridge readers. Branch after it lands.
- [aft-reanchor-observed-times.md](aft-reanchor-observed-times.md) and
  [response-scale-rows.md](response-scale-rows.md) edit the response's side of the same files; no
  shared function. Serial all the same.

## Calls made in planning

Each is the planner's, reversible, and open for the maintainer's mark.

1. The file keeps its name; the alternative was a new name with the register's four record lines and
   the TODO repointed.
2. "The distinct points the values supply" is read as today's derivation with equal neighbours
   dropped: the uniform rule's `n.cuts` evenly spaced points, the quantile rule's chosen midpoints.
   So a column of few values far apart keeps its grid and its bits. The alternative, one point per
   gap between distinct values under the uniform rule too, changes every fit with a discrete column.
3. A refresh that changes a column's count keeps each split's position, as a refresh and
   `setCutPoints` do today. A split past a shorter grid is merged when the update is forced; when it
   is not, the call returns `FALSE` with column and grid as they were, as for any column the trees
   cannot hold. dec-B299's "never fails" and "never keeps the old grid" are read of an accepted
   refresh. The alternatives: move the column's splits by value when its count changes, which keeps
   more of them and needs the trees put back on a rollback; or rescale positions to the new count.
4. An all-missing column gets one point at 0, as the quantile rule gives it today.
5. A warm-start donor whose grid repeats is refused as a state is; dec-B300 names `setState`, a copy
   and a reload. The alternative was dropping the donor's repeats for the scratch build.
6. The refusal is raised in the bridge, where the column and value can be named, with the engine's
   check as backstop; the message text is step 4's.
7. dec-B285's exception for the held grid is removed as dead, with its sentence in the manual.
8. The interval check at a build merges and reports `FALSE`; it does not refuse. Refusing was
   prototyped for the plan this replaces and refused hand-built states that install and merge today.
9. `setPredictor`'s error for too few induced cut points goes; the bridge's message is reworded for
   the one case left.
10. Two exact arms, no new equivalence scenario, no baseline recorded again.
11. A state of an earlier build of this branch whose grid repeats is refused, not converted, and the
    format number does not move.
12. Two pushes are allowed, not required.
13. The design and its critique against the tip: the design's version of this option refused a
    refresh that would hold fewer points and left a column created constant at one point, both since
    ruled otherwise; it called the exact gates foregone, and two of them build a repeating grid; its
    fourteen places are thirteen here and two prechecks, the two uniform scans sharing one fill and
    the bridge's two state readers counted singly.
    The critique's findings on positions, their blocks and their verdict no longer apply.
