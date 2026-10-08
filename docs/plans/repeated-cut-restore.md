# repeated-cut-restore: no cut grid repeats a point

Status: LANDED 2026-10-08 (92f76b34 to 2feba5c5; dec-B285, dec-B297 to dec-B300, dec-B311, dec-B312,
with dec-A183 for the calls made in it). Rewritten 2026-10-07 after the four rulings; the plan this
replaces kept repeated points and stored each split's position, and none of that is built. The file keeps its name because the TODO item, four register entries and
[cut-points-undo.md](cut-points-undo.md) cite it.

agent: opus implementer, one (engine, bridge, tests, gate arms); opus reviewer, told to refute.
rng: by fit and by call.
- POSTERIOR-CHANGING at creation and `setData`, for a sampler with a numeric column whose grid today
  holds a point twice and on which a split can survive: under the uniform rule a column narrower than
  its `n.cuts` evenly spaced points resolve, and a column with one finite value beside missing values
  or `Inf`; under the quantile rule a column where two of the chosen midpoints round to one double.
  The tree prior over that column's thresholds changes.
- POSTERIOR-CHANGING after a refresh, `setPredictor(updateCutPoints = "position")` or the logical
  `TRUE` it replaces, wherever the refreshed grid differs from today's: the count a column held is no
  longer kept (a column created with few distinct values under the quantile rule; a grid set by
  `setCutPoints` or brought by a state at another length), a refresh onto too few distinct values
  shrinks where it failed or kept the old grid, and a refreshed grid drops repeats. Where the count
  changes each split's position is rescaled (rule 6). A refresh whose count does not change keeps
  every position and is NEUTRAL, bit for bit.
- POSTERIOR-CHANGING after `setCutPoints` under its default, `splits = "position"`, on a grid of
  another length than the column holds: each split's position is rescaled (rule 6), where 0.9-34 and
  the base build keep its index and merge a split past the end of a shorter grid. So that call draws
  differently from the release; ruled as built (dec-B311). On a grid of the held length it is NEUTRAL,
  bit for bit.
- New surface: `updateCutPoints = "value"` and `setCutPoints(splits = "value")` move each split to the
  point nearest its old threshold. No call of the base build reaches them.
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
- Tests that pin what changes (read, by search; the pins named by their text are gone since the
  build, rewritten to the new rule). The refusal of a coarser quantile refresh:
  retired: ["induced cut points in new predictor less than previous"](../../inst/tinytest/test-quantile-grid.R),
  retired: ["induced cut points"](../../inst/tinytest/test-bartcore.R),
  retired: ["induced cut points"](../../inst/tinytest/test-mutate-sparse-valued.R),
  [`testQuantilePredictorUpdate`](../../tests/cpp/test_moves.cpp). The held count:
  retired: ["a refresh spreads the count the column holds"](../../inst/tinytest/test-quantile-grid.R). The old
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

6. dec-B311 and, for the spelling, dec-B312. The maintainer on 2026-10-07, asked where the splits
   already in the trees go when a refresh re-derives a grid: "Why not let it be an option? Defaulting to
   value.", and a minute later: "Oh, wait a sec, keep the default at relative so it stays the same
   as 0.9-34." On the spelling: "given that if it were an argument, it would only make sense if
   `updateCutPoints` was `TRUE`, it should probably be the value of the existing argument", "Sure,
   \"none\", \"position\", and \"value\" works for me." and "If `updateCutPoints` is a logical, a
   one-time warning should be issued, scheduled for removal in 1.1-0." The rule, which replaces call 3:
   - `updateCutPoints` is one of `"none"` (the default: the grid is kept), `"position"` and `"value"`.
   - `"position"`: each split keeps its fraction of the grid. With the count unchanged the position
     is kept, bit for bit what `TRUE` did in 0.9-34 and does on the base build. With the count
     changed the position is rescaled to the new count; it is not kept as an index and nothing is
     dropped for being past the end.
   - `"value"`: each split on a refreshed column moves to the new point nearest its old threshold,
     inside the interval its ancestors leave, an empty interval merged: the move `setData` and a warm
     start make ([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)).
   - Under either word an unforced refresh the trees cannot hold after the move returns `FALSE` with
     the column, the grid and the trees as they were; a forced one merges what cannot stand.
   - A logical is taken until 1.1-0 with a warning once per session, `TRUE` as `"position"` and
     `FALSE` as `"none"`; `NA` and anything else is refused by name.
   - `setCutPoints`, asked of the maintainer the same day against always by value and always by
     position: "The same choice as the refresh, same words." It gains `splits`, after `updateState`
     so positional calls keep their meaning: `"position"` (the default; with the count unchanged
     today's behaviour bit for bit, with it changed the position rescaled, where today and in
     0.9-34 the index is kept and a split past the end merged) or `"value"`. Shown that the release
     keeps the index there, against the rescale as built, the maintainer on 2026-10-07: "Rescale, as
     built." (dec-B311). It stays forced: what cannot stand is merged, and there is no unforced form.

   ```r
   s$setPredictor(x, 2L, updateCutPoints = "position")   # TRUE's rule
   s$setPredictor(x, 2L, updateCutPoints = "value")      # splits follow their thresholds
   s$setPredictor(x, 2L, updateCutPoints = TRUE)         # "position", with a warning once
   s$setCutPoints(cuts, 2L, splits = "value")            # the default is "position"
   ```

What the rulings leave open is under Calls made in planning and Calls made in building; the one of
planning that shapes the build is how a rule's points are counted (call 2).

## Constraints

- A sampler none of whose grids repeats a point today, and which never refreshes one, draws, stores,
  restores and reports exactly what it does now.
- The tree prior, the moves, `getTrees`, `printTrees` and `predict` are not touched. A constant column
  keeps one point and stays available, as it is under the quantile rule today.
- Cut points are derived over every row in the design, a row at weight 0 or masked out included
  (dec-B296 leaves them out of its rule).
- State format: no block added, renamed or changed; [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)
  and the floor stay at 1.
- [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move. Three facade virtuals gain
  a trailing argument, the placement of the splits
  ([`SamplerBase::setPredictor`](../../src/bartcore/facade.hpp),
  [`SamplerBase::updatePredictor`](../../src/bartcore/facade.hpp),
  [`SamplerBase::setCutPoints`](../../src/bartcore/facade.hpp)): the bridge holds a sampler through
  the facade and has no other way to hand the rule down (Calls made in building, 4; ruled out here
  before rule 6). `--preclean` on every engine commit.
- A refusal leaves the sampler, its stored state and its generators as they were.
- The line the mutation battery anchors in [`runPredictorTransaction`](../../src/bartcore/sampler.hpp)
  stays as it is.
- `updateCutPoints` takes three words where it took a logical, and `setCutPoints` gains `splits`
  (rule 6); no other R argument is added. `setCutPoints` keeps dec-B285's sort and refusal.

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
   numeric column; the two degeneracy checks, nothing else reading them, are gone
   (retired: [`valuesAreDegenerate`](../../src/bartcore/data.hpp),
   retired: [`cscColumnIsDegenerate`](../../src/bartcore/data.hpp)).
   [`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp) and
   [`SubsetUpdate`](../../src/bartcore/sampler.hpp) snapshot and put back the counts with the points.
   The splits (rule 6, amended 2026-10-07; the check for a rule past a grid that this step first
   named is not built, no refresh leaving one). The engine's two predictor entries take the rule as
   a value in place of the flag: none, position or value.
   [`runPredictorTransaction`](../../src/bartcore/sampler.hpp) keeps the old grid of each column the
   move concerns - under position one whose count changed, under value one whose grid changed - and
   hands them to the chains; with none to hand, the path is today's.
   [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp) takes the rule and skips a column
   with no old grid; under position the target of position i of n on a grid of m is
   floor((2 i + 1) m / (2 n)), then placed inside the interval as the move by value places it. Forced:
   [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp) moves each tree before it re-routes
   and merges, the variance forest through its own remap. Unforced:
   [`Chain::revalidateTrees`](../../src/bartcore/chain.hpp) and
   [`Chain::revalidateVarianceTrees`](../../src/bartcore/chain.hpp) move each surviving tree with a
   routine that merges nothing, records each position it changes and fails on an empty interval;
   a failed transaction puts the recorded positions back before the grid. The bridge's two entries
   read the rule from the argument that carried the flag, so no entry gains an argument. In R
   [`resolveUpdateCutPoints`](../../R/tombstones.R) reads a logical as its word and warns once, in
   the idiom of [`noOpThreadMethod`](../../R/tombstones.R) and as a row of
   [`dbartsTombstones`](../../R/tombstones.R), and
   [`matchCutPointRule`](../../R/bartcore.R) matches the word and refuses the rest; the partial
   update refuses any word but `"none"`. [`Sampler::setCutPoints`](../../src/bartcore/sampler.hpp) takes the same rule
   and hands the forced refresh the old grids of the columns it names; its bridge entry gains the
   argument and [`dbartsSampler$setCutPoints`](../../man/dbartsSampler-class.Rd) gains `splits`,
   matched as `updateCutPoints`' words are. The message of
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
   grid. In [test_moves.cpp](../../tests/cpp/test_moves.cpp): `testRefreshCountChange`, under
   position: a shrink forced (every position the rescaled one or merged, none past the end, fits
   finite, state round trip) and unforced (accepted with the rescaled positions, or rolled back with
   points, counts, codes and trees as before), a growth from one point, an unchanged count leaving
   every position, the variance forest through both; `testRefreshByValue`: forced, every surviving
   split on the point nearest its old threshold inside its interval; unforced, accepted or rolled
   back with the trees as they were, on a tree made to fail; a refresh that leaves the grid as it
   was moves nothing; the variance forest;
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
     where it returns `FALSE`, data, grid, trees and next draws a twin's. A refresh onto a constant
     leaves one point.
   - The three words: `"none"` keeps the grid; `"position"` with the count unchanged draws what the
     base build's `TRUE` draws and `TRUE` draws what `"position"` draws; `"value"` puts each split on
     the point nearest its old threshold. A logical warns once and only once in a session, counted
     in process as the thread methods' warning is; `NA`, a number, two words and an unknown word are
     refused by name; an abbreviation is taken; the partial update refuses both words.
   - `setCutPoints(splits = )`: `"position"` on a grid of the held count draws what the base build
     draws; on a longer and a shorter grid every position is the rescaled one; `"value"` leaves each
     threshold on the given point nearest it; an unknown word is refused.
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
     and its three words; a refresh derives the grid a new sampler would have for the values, up to
     `n.cuts`, fewer where they supply fewer; where the splits go under each word; a logical and its
     removal in 1.1-0. `cuts`: the exception for the held grid and its sentence on equal
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
  - a rollback puts back the points without the count, or without the positions:
    `testRefreshCountChange` and `testRefreshByValue` unforced;
  - position rescales by another formula than the identity at an unchanged count, or keeps the
    index when the count changes: `testRefreshCountChange`; value keeps positions:
    `testRefreshByValue`, tinytest's words;
  - a logical warns on every call, or `TRUE` is read as `"value"`: tinytest's words;
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
they induce fewer; `updateCutPoints` takes `"none"`, `"position"` and `"value"`, a logical being
taken with a warning until 1.1-0; `setCutPoints` takes the same choice as `splits`, and by position
rescales where 0.9-34 kept the index. The existing item on quantile cut points, which named a refresh
onto more distinct values than the column's cut points could separate, is amended in the fix round:
such a refresh now derives up to `n.cuts` of them. The item that lists what expires in 1.1-0 names the
logical.

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
3. Replaced by rule 6 on 2026-10-07, for a refresh and for `setCutPoints`. As planned: A
   refresh that changes a column's count keeps each split's position, as a refresh and
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

## Recheck against the tip

By the implementer, at b2ddb2a4, before any code. The plan was reworked at d734af08; the second push of
the leaf conversions, a forest's defaults by kind and selection by label have landed since.

- Nothing the steps edit moved: data.hpp, tree.hpp, chain.hpp, sampler.hpp, the four tests/cpp files
  of step 5, test-cut-points-undo.R, test-quantile-grid.R, benchmarks/ and the exact-gates workflow
  are byte for byte those of d734af08. The bridge moved at the response setters and `setData` only (the
  gp refusals); its two `cutPoints` readers, `bartcore_setCutPoints` and the two predictor entries are
  as planned. Every cite of this plan resolves.
- The planning probes, run again on a shipped build of the tip in a library of its own, print what
  they printed: the counts of Context at creation, a refresh and `setData`, the accepted and refused
  grids, the twelve seeded digests, and the two arms' base grids (0, 1, 2, 2, 2 in units of eps above
  1; five copies of 1).
- Eight tinytest files name `updateCutPoints`, seven with `TRUE`, as counted. Of the tinytest files
  that landed since, one names a grid: test-leaf-conversions.R sets a recipient's grid with
  `setCutPoints`; it is run with the eight.
- No step's premise moved. Calls the implementer makes beyond the plan are listed under
  Calls made in building, added as they are made.

## Calls made in building

The implementer's unless marked, reversible, and open for the maintainer's mark.

1. Steps 2 and 3 land as one commit: with the drop alone a refresh would shorten a grid whose count
   no rollback puts back.
2. The rescaling under `"position"`: position i of n, counted from 0, goes to
   floor((2 i + 1) m / (2 n)) on a grid of m, the point under the centre of the split's share of the
   old grid. It is the identity at m = n, keeps order, and is the index the quantile rule uses to
   spread a count over midpoints. Two splits that land on one point are then placed as the move by
   value places them: inside the interval their ancestors leave, an empty interval merged.
3. An unforced refresh "the trees cannot hold" is one after whose move a split's interval is empty,
   a leaf is empty, or a monotone tree is out of order. A split moved to another point of its
   interval is held.
4. The engine takes the rule beside its flag for a refresh, not in its place as step 3 words it: the
   flag stays a logical and a trailing value says position or value, so no engine caller changes.
   Three facade virtuals gain that argument, which the Constraints had ruled out before rule 6.
   Under `"position"` a column whose count did not change is not touched, and under `"value"` a
   column whose refreshed grid is bit for bit the old one is not: such a refresh runs the base
   build's path.
5. A forced refresh and `setCutPoints` merge an emptied interval weighing each leaf by the rows it
   held before the change, as the merge of an empty leaf weighs them. `setData`'s move, the same
   routine, keeps weighing by the statistics the last sweep left, which is all it has; they differ
   between a sampler and its copy, and a copy must draw what its original draws after one call.
6. The step-3 check for a rule past a column's grid in the unforced validation is not built: under
   rule 6 no refresh and no `setCutPoints` leaves one.
7. The orchestrator's readings, not the maintainer's words, awaiting the maintainer's mark: a position is
   rescaled when a refresh changes a column's count, where the release has no such case, and by the
   formula of call 2; an explicit `FALSE` warns as `TRUE` does; the word is matched as `match.arg`
   matches, a unique abbreviation taken; `setCutPoints` stays forced. The rescale at `setCutPoints`
   on a grid of another length is not among them: there the release keeps the index, and the
   maintainer ruled "Rescale, as built." (dec-B311; rule 6).
8. The bridge reads the rule from the argument that carried the logical, as an integer code, so
   neither predictor entry changes its arity; the `setCutPoints` entry gains one argument. The flat C header carries no predictor update and does
   not move.

## Calls made in the fix round

After the review, by the implementer who took the slice over; reversible, open for the maintainer's mark.

1. The logical's reader moves to the tombstone file and does one thing, a logical to its word with the
   warning; the word is matched by the caller. Deleting the tombstones leaves a call that no longer
   resolves, as for the thread methods. Its row of
   [`dbartsTombstones`](../../R/tombstones.R) is ["logical updateCutPoints"](../../R/tombstones.R),
   kind `behaviour`, owner `dbartsSampler`. The tombstone test asks the news list for functions,
   methods and arguments only; the logical is named there all the same.
2. The news item on quantile cut points is amended, not the plan: against 0.9-34 a refresh onto more
   distinct values no longer keeps the lowest ones, and it no longer keeps the column's count either.
3. How a forced merge weighs leaves is pinned twice.
   [`testMapOldCutPointsLiveRowsMerge`](../../tests/cpp/test_data.cpp) merges two leaves of one row and
   five under node statistics set the other way round, with unit weights and with weights of 4 and
   2.5. From R,
   ["a forced merge weighs each leaf by the rows it holds"](../../inst/tinytest/test-cut-grid-distinct.R)
   gives a column of a weighted fit one point and works every tree out by hand from the trees as they
   were, for a sampler and for its copy. It holds all twenty trees to the shape worked out, so it
   also pins that the re-routing after the move merges nothing in that fixture.
4. ["an unforced refresh moves the splits of every forest"](../../inst/tinytest/test-cut-grid-distinct.R)
   sets a column of a two-forest sampler to every other point of its grid and refreshes it onto the
   values it holds: each split's position doubles and its threshold stays, so the trees are identical
   where every forest was moved.
5. Checks tightened beyond the lines the review named: at `setCutPoints` on a finer grid each split is
   held to twice its position, as the split-by-split check on a grid of the held count is; the
   unforced refresh by value in [`testRefreshByValue`](../../tests/cpp/test_moves.cpp) is held to the
   outcome it has, taken; the composition test gathers the column's splits over twenty single draws
   and holds them to the one point. The `expect_equal` on a rolled-back sampler's next draws stays.
6. The equivalence harness passes `"none"` where it passed `FALSE`; its compares stay bitwise.

Mutants run for these, each a build of the fixed tree: node statistics for rows at the merge, counts
for weights at the merge, and an unforced refresh that moves the first forest only. Each fails the
check named for it in 3 and 4.

## Landing note

Landed 2026-10-08 as 92f76b34 to 2feba5c5, after two looks by independent reviewers told to refute.
The first found the engine sound, two blocking gaps and three mutants no test killed: how a forced
merge weighs leaves, by node statistics or counts on a weighted fit, and an unforced refresh moving
the first forest only. The fix round pinned each, shown failing on its mutant, and rebased over the
small R batch. The second could not break it: 468 split placements over 13 sampler kinds, both grid
rules, both words and three grid lengths against its own reference, no mismatch; copies and reloads
identical; every unforced refusal restored exactly; sanitizers clean on the R-loaded path. It left
two lines of the R tests that cannot fail, a split being read back off its grid, and a build keeping
the index on a shorter grid at setCutPoints that passes every R test and fails the C++ tests. Gates
on the fixed branch: tests/cpp, under ASan and UBSan too; the tinytest suite at home, 19288 results,
none failed; the snapshot files and the three bitwise compares on the reference build; the exact
gates in quick mode, 32 of 32; lintr, air, rc-codoc, win-drift, doc-freshness, anchors, the news
parse; check with one note. On the landed build stan4bart's suite ran 582 results, bartCause's 1412
expectations and bairrtt's 207 results, none failed.

The push failed one older test off arm64: a copy and its original, given the same setCutPoints, drew
different values on ubuntu, windows and the sanitizer jobs. On x86 the two stored states are bit for
bit the same after the call and the draws part at the third sweep by 7e-16 and stay there; a copy
with nothing set parts at the first sweep on both platforms and both commits, its leaves holding
their rows in another order, and a twin matches for 50 sweeps; valgrind is clean and no SIMD level
changes it. More splits survive a shrinking grid now, so more leaves sum in another order. The test
now uses a twin, as the block after it did. Call 5's promise holds where a leaf's weights sum exactly;
with weights that do not, a copy's merged leaf can differ from its original's by one unit in the last
place (root TODO, copy-merge-weight-order).
