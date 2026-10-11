# warm-start-salvage: a warm start takes what fits of a donor

Status: PLANNED 2026-10-11, in two slices, A then B. Nothing is built. The blind critique is owed before
either slice goes to an implementer. No call is open. Design: [warm-start-salvage.md](../design/warm-start-salvage.md), whose six choices are
settled: (a) `forceUpdate = FALSE`, (b) FALSE for a split that has to move or merge, (c) the mean weighted by
rows, (d) the three errors stay, (e) the mean forest's rule for the variance forest (dec-B448), (f) one line.

agent: per slice, one opus implementer and one opus reviewer who runs that slice's mutants (Verification).
Tier "Changes draws" of [README.md](README.md): both slices edit the engine's install path.
rng: slice A is NEUTRAL, bit for bit, for every fit that makes no warm start and for a forced warm start
from a donor whose every live threshold is a point of the recipient's cut grid, which runs the code that runs
today. It is SHIFTING for a forced warm start across grids, `bart(warm.start = )` onto other rows included:
the starting trees are the across-grids routine's and not today's remap's, so the draws that follow differ
and the posterior does not. A forced start into a recipient whose column restriction or interaction limit the
donor breaks is an error today and has no earlier draws. Slice B is NEUTRAL for every donor of the
recipient's tree counts, leaf kind, DART setting and variance forest; everything else it installs is an
error today. No install draws beyond dec-B440's. No recorded scenario, snapshot file or exact gate calls
`installTrees` or `warm.start` (searched), so no baseline is re-recorded.
window: after state-grid-install lands (this plan's base is its tip, 9fb70c48) and before the merge to main
(TODO `forced-update-returns-null`). Engine slices stay serial; B follows A.
budget: slice A about 1370 lines (engine 300 added and about 290 removed, bridge and R 150, tests 750, help,
NEWS and notes 170), upper figure 2200. Slice B about 1520 (engine 380, bridge and R 210, tests 800, help,
NEWS and notes 130), upper figure 2500. Stops: Budget and stops.

## Goal

`installTrees(donor, samples = NULL, forceUpdate = FALSE)` asks first: a donor every tree of which stands in
the recipient as stored is installed and the call returns TRUE; one that would have to be changed is not
installed, the call returns FALSE and the sampler is as it was. With `forceUpdate = TRUE` the donor is
installed with the changes and the call returns NULL invisibly. `bart(warm.start = )` forces, and with
`verbose = TRUE` prints one line counting the trees taken, moved and merged. After slice B a donor of another
tree count, of constant against linear leaves, or with DART or a variance forest on one side only is taken
in part. Closes TODO `warm-start-salvage`, `warm-start-tree-count` and `forced-update-returns-null`.

## Context

Rulings built, in the register's words unless quoted as the maintainer's:

- dec-B295, the maintainer: "Warm starts should try to get as much information out of the fit it can. If
  trees can be collapsed to be made into valid ones, that would make sense. Install trees might have modes
  that work a force flag, where if it was absent it would accept/reject." Slice A.
- dec-B310, the maintainer: "Make them NULL. The return value when not forced is an indicator of whether or
  not a statistically valid update was made." Built for `setPredictor` and `setState`; slice A builds it for
  `installTrees`.
- dec-B398, "Repair, as we do elsewhere.": the cut-back at a forbidden split and the merged leaf at the mean
  weighted by rows. Built for `setState`; slice A reaches it from the warm start.
- dec-B326, "One line under verbose = TRUE, alongside the setup bart already prints." Slice A, with slice
  B's phrases.
- dec-B323, "Sure, option B.": "a warm start into fewer trees takes the donor's trees whose fitted values
  vary most over the recipient's rows, in the donor's order and with their values as they are; into more
  trees it takes them all and starts the rest as single leaves at 0". dec-B448, "Use the same rule.": the
  variance forest, "selection by variation of the tree's log scale over the recipient's rows and padding
  with leaves at 1". dec-B247 is the request both answer. Slice B.
- dec-B324, "Take the trees; a part only one side has does not transfer." dec-B325, "Take what translates."
  Slice B.

Rulings that stand as built and bound this plan: dec-B343 (sigma and k only where the recipient draws them),
dec-B418 (nothing converted between response scales), dec-B440 (the side for missing values drawn before the
trees are judged), dec-B328 and dec-B329 (a linear slope converted by the formula), dec-B449 and dec-B450
(the cut grid and its weights are the sampler's; a warm start never wrote them), dec-B305 and dec-B234 (what
clean means), dec-B278 and dec-B355 (the monotone reset; dec-B355 puts the warm start under dec-B295),
dec-A199 and dec-A201 (the build calls of the two state slices).

Out of scope, each with its TODO item: columns and factor levels matched by position
(`warm-start-by-position`); a linear donor on other covariates in the same number
(`warm-start-linear-covariates`); `bart`'s refusal of `warm.start` for nbinom, ordinal, multinomial and
hurdle fits (`warm-start-bart-families`); a donor that is another model (`warm-start-across-models`, after
1.0-0, dec-B248); a test of the undo of a failed install (`warm-start-undo-test`). A recipient of two or
more forests stays refused ([`refuseMultiForestWarmStart`](../../R/bartcore.R)); no item asks for it.

Released 0.9-34 has neither `installTrees` nor `warm.start` (run on scratch/libs/cran: the sampler's methods
matching install, warm or State are `setState` and `storeState`; `bart`, `bart2` and `rbart_vi` have no
`warm.start`). So no call's behavior changes from a release, and NEWS describes a new feature.

Code as it stands:

- [`Sampler::installForests`](../../src/bartcore/sampler.hpp) assembles each chain's trees from the donor
  sample it takes, restates linear coefficients
  ([`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp)), draws the missing-value sides
  ([`Chain::drawDonorMissingDirections`](../../src/bartcore/chain.hpp)), refuses on four scratch passes
  ([`Chain::interactionStateFeasible`](../../src/bartcore/chain.hpp),
  [`Chain::columnMaskStateFeasible`](../../src/bartcore/chain.hpp),
  [`Chain::forestsRebuildable`](../../src/bartcore/chain.hpp),
  [`Chain::varianceForestRebuildable`](../../src/bartcore/chain.hpp)) and installs by
  [`Chain::installForest`](../../src/bartcore/chain.hpp) and
  [`Chain::installVarianceForest`](../../src/bartcore/chain.hpp).
- A donor on the recipient's cut points goes through
  [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), the routine `setState` installs by: it stages
  each tree ([`Tree::stageFromFlat`](../../src/bartcore/tree.hpp)), which arms the forbidden-split collapse
  of [`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp). Only the refusals ahead of it keep a warm
  start from cutting a tree back today.
- A donor on other cut points goes through
  [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp): built on the donor's grid with no rows
  under it, moved, then partitioned, so a merge at the move takes the plain mean.
- [`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp) is the state grid slice's routine: the
  recorded grid goes in with re-quantized codes, the install step it is handed builds and partitions the
  trees there, the store's grid is put back, and
  [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp) moves the splits by value with the rows each
  leaf held. [`Chain::liveThresholdsAreOnGrid`](../../src/bartcore/chain.hpp) reads whether a state's live
  thresholds are points of a grid with no tree built.
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) holds the verdict's staging in its two tree loops.
  It cannot be called on a warm start's assembled chain, which carries no latents or kept draws.
- The bridge's [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads the donor's leaf blocks in
  the recipient's shape, and [`bartcore_installForests`](../../src/R_interface_bartcore.cpp) takes three
  arguments and no predictor matrix. [`installTrees`](../../R/dbarts.R) returns NULL invisibly or stops.
- The verbose setup lines are printed at creation by the bridge
  ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)); `bart` installs the donor after creation
  and before its first run.

### Measured on this tree

Run on scratch/libs/sgi (shipped build of 9fb70c48 with the design note on top), gaussian unless said, 300
rows, 4 numeric columns, 10 to 50 trees, one chain. Every error left the recipient's stored state identical.

| donor against recipient | sgi |
|---|---|
| same rows, same shape | taken, NULL invisibly; starting fit equals the donor's (difference 0) |
| recipient on 40 of the donor's rows, the donor's cut points set | taken; 46 nodes to 42. `setState` of the same trees: FALSE; forced, the same trees exactly |
| recipient on other rows, its own grid | taken; 0 of 18 donor splits on the recipient's grid; 46 nodes to 46; rmse to the donor's function 0.18 against a fit sd of 1.79. `setState`: FALSE, and forced NULL with the recipient's grid kept |
| a column holding one value | taken; the donor's 5 splits on it gone (108 nodes to 96) |
| another number of predictors | "malformed cut points in warm-start donor" |
| numeric against an unordered factor | the predictor structure error |
| numeric, or an ordered factor of 6 levels, into an ordered factor of 4 | taken by position; `setState` of the same state is an error in both forms |
| levels a, b, c against c, b, a | taken; correlation with the donor's fit -0.998 |
| another tree count, either way | the shape error |
| DART on one side only | the DART error; both DART: taken |
| a variance forest on one side only; both, with another count | the shape error; the variance count error |
| constant into linear or gp, linear against gp, linear with another number of covariates, gp on another number of rows | "malformed leaf parameters in bartcore state" |
| gp on other rows of the same number | taken, every fit at one value |
| linear into constant | taken, leaf values the donor's intercepts exactly |
| gp into constant | taken; starting range -1.95 to 4.99 against -2.15 to 5.43 |
| linear on x2 into linear on x3 | taken |
| `interactions(max.order = 1)`, `blocks`, a restricted variance forest | each its error. `setState` of the same trees: FALSE; forced, 108 nodes to 92 and to 58, the fit's mean over the rows unchanged |
| free donor into a monotone recipient | taken; 11 of 30 trees with every leaf 0. `setState`: FALSE |
| missing values on one side only | taken both ways. `setState` of a donor with sides onto complete data: FALSE |
| gaussian into binary, binary into gaussian | taken, fits finite after 3 sweeps (2 of the design's 49 pairs) |
| `bart(family = "nbinom", warm.start = )` | "does not support 'warm.start' or 'n.grow.sweeps'" |
| recipient of two forests | the multi-forest error |
| `bart(verbose = TRUE, warm.start = )` | 34 lines of output, none about the warm start |

What the state grid slice changed under the design's sections 1, 3 and 4:

- A donor's state on another grid is declined by an unforced `setState` and moved by a forced one, the
  recipient's grid kept. The design ran its table when `setState` took the state's grid.
- Today's warm start and the across-grids routine give other trees on a small recipient. A 50-tree donor of
  196 nodes onto new rows, 8 recipients a size, nodes kept and rmse to the donor's function:

  | rows | today's warm start | the routine (forced `setState`) |
  |---|---|---|
  | 30 | 184.0, 0.694 | 173.8, 0.675 |
  | 100 | 188.5, 0.392 | 188.8, 0.380 |
  | 300 | 191.8, 0.388 | 194.5, 0.388 |
  | 1000 | 196.2, 0.283 | 196.0, 0.283 |

  The routine merges a leaf no row reaches under the donor's own thresholds before any split moves; today's
  arm moves first, and a moved split always leaves a row on each side.
- The planned cross-grid arm can be written by hand on sgi: `setCutPoints` with the donor's grid, today's
  same-grid `installTrees`, then `setCutPoints(own, splits = "value")`. At 30 and 60 rows that gives exactly
  the routine's trees (76 and 72 nodes, leaf values equal, the grid back identical) where today's arm gives
  78 and 76.
- Section 4's staging of a donor on another grid is not needed: across grids the unforced answer is FALSE by
  choice (b) with no tree staged, and the forced install is the routine.

## The rule

The first row that applies. "Error" leaves the sampler untouched in both forms.

| The donor against the recipient | Unforced | Forced | Slice |
|---|---|---|---|
| not a donor; an older encoding; a malformed block or tree; a donor grid that repeats a point; a donor threshold that is no point of the donor's own grid; a bad `samples`; a recipient of two or more forests | error | error | as today |
| another number of predictors; a column an unordered factor on one side only; on a column the recipient holds as an ordered factor, another donor grid; a gp leaf against another kind; a gp donor on another number of rows | error | error | A keeps today's texts and adds the ordered factor; B writes the texts of (d) |
| every tree stands as stored | TRUE | NULL; today's install bit for bit | A |
| a split whose threshold is no point of the recipient's cut grid | FALSE | NULL; the across-grids move | A |
| a leaf no row reaches; a split outside its interval; a side for missing values on a column that has never held one here | FALSE | NULL; merged or dropped, as today | A |
| a split on a column its forest may not use, or past its `interactions` limit | FALSE | NULL; cut back at the first such split from the root | A |
| a monotone tree with leaf values out of order | FALSE | NULL; every leaf 0 | A |
| gp leaves on both sides, the donor recorded on other cut points | FALSE | NULL; trees taken, fits started at 0, as today | A |
| a recipient of more trees, mean or variance | TRUE | NULL; all taken in the donor's order, the rest single leaves at 0, or at 1 in a variance forest | B |
| a recipient of fewer trees | FALSE | NULL; the trees that vary most | B |
| constant donor into linear leaves | TRUE; slopes 0 | NULL | B |
| linear donor into constant leaves, or into linear leaves on another number of covariates | FALSE | NULL; each leaf at its stored intercept, slopes 0 | B |
| DART, or a variance forest, on one side only | not judged | the trees taken; the donor's part left behind, the recipient's started as a new fit's | B |

- Values are as `setState`'s ([install-surface.md](../design/install-surface.md)): unforced TRUE and FALSE
  are returned visibly; forced returns NULL invisibly whether or not anything was changed.
- FALSE leaves the trees, drawn values, generators, kept draws and the `state` field as they were. The
  sides drawn under dec-B440 are drawn on scratch first, and a decline puts the generators back.
- An install, forced or not, drops the recipient's kept draws and stores its state, as today.
- The missing-value draw of dec-B440 is made before the trees are judged. Such a start is clean unless a
  drawn side leaves a leaf no row reaches or a monotone tree out of order.
- Another response family, another response scale, and the donor's k and sigma are not judged (dec-B418,
  dec-B343).
- The forced move across grids is `setState`'s three passes: each tree built as stored and the recipient's
  rows routed through it by the donor's thresholds; each split on a column whose grids differ moved to the
  nearest point of the recipient's grid inside the interval the splits above it leave, a split with no point
  left becoming one leaf at the mean of the leaves beneath weighted by the rows each held; the rows routed
  again and a leaf left empty merged. The recipient's cut points, weights and codes are the same before and
  after, byte for byte.
- A merged leaf is the mean of the leaves it replaces weighted by the recipient's rows in each under its
  working weights, geometric in a variance forest, the plain mean where none holds a row (dec-B398).
- Several chains: one chain's tree that needs a change declines the whole donor.
- The line. `bart` prints it after its setup lines when `verbose = TRUE`, at every warm start, and nothing
  otherwise; no warning is raised. Trees are counted per chain and changes summed over chains.

      warm start: 75 of 75 donor trees taken
      warm start: 50 of 50 donor trees taken; 12 moved to this fit's cut points, 21 merged, over 4 chains
      warm start: 30 of 30 donor trees taken; 11 reset to 0
      warm start: 50 of 200 donor trees taken; 9 merged; variance forest started fresh

  "moved" counts trees that kept every split and hold one on another cut point; "merged" counts trees that
  lost a split; the two are disjoint. Slice B adds "150 trees started as single leaves", "slopes left
  behind", "slopes started at 0", "DART weights started fresh", "donor's DART weights left behind",
  "variance forest started fresh", "donor's variance forest left behind" and "variance trees: 4 of 6 taken".

## Constraints

- A forced install from a donor whose live thresholds are all points of the recipient's grid runs
  [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) and
  [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp) with the arguments they get today.
- The unforced answer is FALSE exactly when the forced install changes a tree. The verdict stages a tree by
  the routine the install stages it by.
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) returns what it returns today for every state; its
  tree loops move into helpers unchanged.
- [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp),
  [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp) and
  [`Sampler::setState`](../../src/bartcore/sampler.hpp) are called, not edited.
  [`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp) gains only the restore on a throw.
- No change to a stored state's format, and no flat entry of
  [dbarts.h](../../inst/include/dbarts/dbarts.h) installs a donor, so the C API hash does not move.
- `installTrees` stays refused on a borrowed view: the forced move re-quantizes.
- The facade virtual changes in both slices: install with `--preclean`.
- The three errors of (d) stay errors under force, and slice A lifts no error but the two containment ones.

## Slice A: the two forms

### Steps

1. Staging shared with the state install (about 110 lines, [chain.hpp](../../src/bartcore/chain.hpp)).
   - Move the body of each per-tree loop of `stateIsValid`, mean and variance, into a private member, and
     call them from `stateIsValid`.
   - Add `donorTreesStage(state, notClean, varianceFailed)`: the two helpers over every live tree of a warm
     start's assembled chain. False where a tree does not build, `varianceFailed` saying which forest, so the
     two existing messages stand. With `notClean` null the trees are only built.
   - Add `coldStartFunctionFits(state)`, compiled for function leaves: every live tree's fits set to zeros.
   - Add `countDonorTrees(taken, counts)`: flattens the live trees
     ([`Chain::getState`](../../src/bartcore/chain.hpp) without the kept draws) and compares each with the
     donor tree it was given: fewer nodes is merged; equal nodes with a threshold changed is moved; in a
     monotone forest, equal structure with every leaf 0 where the donor's were not is reset.
   - Delete `interactionStateFeasible`, `columnMaskStateFeasible`, `forestsRebuildable`,
     `varianceForestRebuildable` and `rebuildLiveForestRemapped`. `installForest` and `installVarianceForest`
     lose their grid and store parameters, the second its remap arm.
2. The install (about 190 lines, [sampler.hpp](../../src/bartcore/sampler.hpp),
   [facade.hpp](../../src/bartcore/facade.hpp)).
   [`Sampler::installForests`](../../src/bartcore/sampler.hpp) becomes
   `installForests(donor, sampleMap, currentPredictors, force = true, report = nullptr)`.
   `currentPredictors` has no default: a dense column re-quantizes from it and a null there is read.
   - The checks, the assembly and the direction draw stay as they are. The grid check also refuses a donor
     grid other than the store's on a column the store holds as an ordered factor.
   - Where the donor's recorded points differ from the store's and the leaf is a function leaf, the fits are
     cold-started on the scratch and the donor is not clean.
   - A live threshold that is no point of the donor's own grid is the rebuild error. `acrossGrids` is true
     where some chain holds a live threshold that is no point of the store's grid.
   - Across grids the trees are built under [`ScopedCutGrid`](../../src/bartcore/data.hpp) on the donor's
     points for validity only, and the donor is not clean. Otherwise `donorTreesStage` judges them on the
     store's grid, with the verdict where the call is unforced.
   - Not clean and unforced returns a new `WarmStartResult::declined` with nothing written; the generator
     guard puts back any draw.
   - The install step is today's two loops on their same-grid arms. Across grids it is handed to
     `installAcrossGrids`; otherwise it is called. A failure or a throw restores every chain touched, as
     today.
   - `report`, when given, is filled per chain from `countDonorTrees`: offered, taken, moved, merged and
     reset for the mean forest, and the same for the variance forest.
   - [`WarmStartResult`](../../src/bartcore/sampler.hpp) loses `interactionMismatch` and
     `columnMaskMismatch`.
   - `installAcrossGrids` puts the store's grid back when its install step throws, and rethrows.
   - Rewrite the comment block above `installForests`. The virtual and its override take the new signature.
3. Bridge and R (about 150 lines).
   - [`bartcore_installForests`](../../src/R_interface_bartcore.cpp) takes the predictor matrix and the
     force flag as a fourth and fifth argument (its row in [R_interface.cpp](../../src/R_interface.cpp), its
     declaration in the bridge header) and returns a list: `installed`, and one integer matrix of counts per
     forest kind, a row per chain.
   - In [R/dbarts.R](../../R/dbarts.R) an internal `installDonorTrees(sampler, donor, samples, force)` does
     what the method does today, passes [`rawPredictorMatrix`](../../R/dbarts.R) of the sampler's
     predictors, stores the state only when the donor went in, and returns the list.
   - [`installTrees`](../../R/dbarts.R) gains `forceUpdate = FALSE`, checked as `setState` checks its own
     ("'forceUpdate' must be TRUE or FALSE"), and returns `installed` visibly when unforced and NULL
     invisibly when forced.
   - `bart` calls `installDonorTrees` with force and, when the control is verbose, prints `warmStartLine` of
     the report with `cat`; `warmStartLine` lives in [R/bart.R](../../R/bart.R).
   - Docstring of `installTrees`.
4. Existing tests (about 300 lines changed). Counts by search:
   - `installTrees(` has 69 hits in 23 files of inst/tinytest: 67 calls and 2 comments. `warm.start =` has
     13 hits in 8 files. `C_dbarts_bartcore_installForests` is called raw in 3 files. `installForests(` has
     44 hits in 4 files of tests/cpp.
   - First run a tracer over the 23 files on the slice's library that prints each call's unforced answer.
     Add `forceUpdate = TRUE` where a test means the repair, and leave a call that returns TRUE alone.
   - By reading, the calls that change are in 9 or 10 files:
     ["column restriction"](../../inst/tinytest/test-blocks.R),
     ["interaction constraint"](../../inst/tinytest/test-interactions.R) (through `bart`, which now fits),
     and the refusals and repairs in [test-single-forest-vars.R](../../inst/tinytest/test-single-forest-vars.R),
     [test-heteroscedastic-warm-start.R](../../inst/tinytest/test-heteroscedastic-warm-start.R),
     [test-monotone.R](../../inst/tinytest/test-monotone.R),
     [test-monotone-unforced.R](../../inst/tinytest/test-monotone-unforced.R),
     [test-state-empty-leaf-merge.R](../../inst/tinytest/test-state-empty-leaf-merge.R),
     [test-state-missing-record.R](../../inst/tinytest/test-state-missing-record.R) and
     [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R).
   - The 3 raw calls gain two arguments.
   - In tests/cpp the 44 calls gain the predictor argument. The 2 `interactionMismatch` and 2
     `columnMaskMismatch` expectations become `declined` unforced and a checked cut-back forced.
     [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp) checks the merged leaf by rows.
     [`testWarmStartStandardization`](../../tests/cpp/test_state.cpp) keeps its gp sentence, "copied on the
     same grid, cold on another".
   - [composition-matrix.R](../../benchmarks/R/composition-matrix.R) has 2 calls to read.
5. New tests (about 450 lines; Tests).
6. Help, NEWS and notes (about 170 lines).
   - [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the usage line; `forceUpdate` gains "For
     `installTrees`, FALSE (the default) takes the donor only if every tree it gives can stand in this
     sampler as the donor holds it, and TRUE takes it changing what must be changed"; `donor` drops "tree
     count, and DART setting" only in slice B.
   - Warm starts, replacing the sentences on refusal and remapping: "`installTrees(donor)` asks first. A
     donor every tree of which fits this sampler as stored is installed and the call returns TRUE. One with
     a tree that would have to be changed is not installed: the call returns FALSE and the sampler is as it
     was, its trees, drawn values, generators, kept draws and `state` field. A tree has to be changed when a
     split's threshold is not one of this sampler's cut points, as for a donor fitted to other rows; when a
     leaf holds none of this sampler's rows; when a split is on a column, or in an interaction, its forest
     does not allow; when it records a side for missing values on a column that has never held one here; or
     when a monotone forest's leaf values are out of order. `installTrees(donor, forceUpdate = TRUE)`
     installs such a donor and returns NULL invisibly, whether or not anything was changed: a split off the
     cut points moves to the nearest one inside the range the splits above it leave or, with none left,
     becomes one leaf with everything beneath it; a split with a side no row reaches, and the first split
     from the root that its forest does not allow, likewise become one leaf, at the mean of the leaves it
     replaces weighted by the rows each held; a side is dropped; and a monotone tree out of order has every
     leaf set to zero. `bart`'s `warm.start` always forces."
   - [bart.Rd](../../man/bart.Rd), `warm.start`: "takes all of the earlier fit that this one can hold: a
     tree this fit cannot hold as it stands is cut back at the split that does not fit. With `verbose =
     TRUE` one line after the setup lines says how many trees were taken, moved to this fit's cut points and
     merged."
   - [interactions.Rd](../../man/interactions.Rd): "enforced on any warm start" becomes "a warm start cuts
     back a tree that breaks it".
   - [monotone.Rd](../../man/monotone.Rd): "`installTrees`" in the list of calls that always complete
     becomes "`installTrees` with `forceUpdate = TRUE`".
   - [NEWS.Rd](../../inst/NEWS.Rd), the Starting points item: "`bart`'s `warm.start` starts the chains from
     a previous fit's trees, taking all of the earlier fit that the new one can hold, and with `verbose =
     TRUE` says in one line how many trees were taken, moved and merged ... the sampler offers both as
     `$installTrees` and `$growFromRoot`; `installTrees(donor)` returns TRUE for a donor taken as stored and
     FALSE, the sampler untouched, for one that would have to be changed, which `forceUpdate = TRUE`
     installs." No sentence compares branch builds.
   - [install-surface.md](../design/install-surface.md): the sentences "A warm start is not a state install"
     and "A warm start still moves", and the closing section, are rewritten to the rule.
     [architecture.md](../architecture.md): the containment sentence.
   - Twelve documents cite the five deleted routines (17 cites, by search, some of them of
     `installVarianceForest`, which stays): a frozen record gets the `retired:` marker and a clause saying
     the routine is gone, a present-tense document a rewritten sentence.
   - TODO: drop `forced-update-returns-null`; `warm-start-salvage` keeps what slice B builds.

Commits, each on a tree that passes: (1) steps 1 to 4; (2) step 5; (3) step 6.

### Tests

tinytest in a new `test-warm-start-forms.R`, tests/cpp in `test_state.cpp`. Each row: what it pins; the
mutant it must kill.

| Test | Pins | Mutant |
|---|---|---|
| clean is TRUE | same rows: TRUE visibly, trees and fit the donor's, the `state` field stored again, kept draws dropped | the verdict always FALSE |
| default declines | donor on 40 of its rows and its cut points: `installTrees(donor)` is FALSE visibly | `forceUpdate = TRUE` as the default |
| forced returns NULL | the same donor forced: NULL invisibly, trees equal a forced `setState`'s; forced on a clean donor: NULL too | TRUE, or a visible value, when forced |
| declined leaves all | after FALSE: trees, sigma, k, the next 5 draws of a twin that never called, `predict` from kept draws, the `state` field | kept draws dropped or the state stored before the decline |
| off the grid declines | donor on other rows: FALSE; forced, the recipient's cut points, weights and codes identical before and after | the verdict not reading thresholds |
| equals the two calls | the same donor forced: trees and leaf values equal a twin's that took the donor's grid, installed and called `setCutPoints(own, splits = "value")`, and a forced `setState`'s on an unweighted pair | the remap arm kept |
| finer grid is clean | a recipient grid holding the donor's points and others: TRUE, splits as stored | the verdict keyed on the grids differing |
| merged leaf by rows | the hand tree of install-surface.md (leaves 0.10, -0.20, 0.05 holding 73, 38, 89 rows) as a donor state forced onto the 50th and 75th points: -0.0027 | the plain mean (-0.05); rows counted after the move (0.10) |
| weights | a case-weighted recipient: a merged leaf weighs by its working weights, by hand | unit weights |
| cut back, same grid | `interactions(max.order = 1)` and `blocks` recipients: FALSE; forced, no forbidden split left, the splits above stay, one tree's merged leaf by hand, the fit's mean over the rows unchanged | the error kept; the leaf at 0 |
| cut back, across grids | the same on other rows | the predicate unarmed on the across-grids arm |
| variance forest | a variance threshold off the grid alone: FALSE; a split a restricted variance forest may not use: FALSE, forced cut back at the geometric mean | variance trees not judged |
| monotone | free donor: FALSE; forced, out-of-order trees at 0 and the rest as stored | the monotone test left out of the verdict |
| missing sides | donor complete, recipient not: TRUE with sides drawn; a decline after a draw leaves the generators (twin draws equal); a donor side on a column never missing here: FALSE, forced dropped | the generator not put back; a dropped side not judged |
| second chain | `samples` giving chain 2 a kept draw that needs a change: FALSE, chain 1 untouched | the first chain only |
| verdict equals forced change | 40 seeded pairs over the kinds above: unforced FALSE exactly when the forced trees differ from the donor's | any gap between verdict and install |
| gp across cut points | same rows, other cut points: FALSE, forced fits 0; the same grid: TRUE, fits copied | fits copied across grids |
| ordered factor | a numeric donor column, and an ordered factor of another level count, into an ordered factor: error in both forms | the routine run on a factor column |
| argument | `1`, `NA`, `c(TRUE, TRUE)` and `"partial"` are errors | coercion |
| errors stay | another number of predictors, factor against numeric, tree count, DART and variance one side, leaf kinds, two forests, a malformed donor: error in both forms, stored state identical | an error lifted by force |
| bart forces and says | `bart(warm.start = )` into an interaction-limited fit runs; one line starting "warm start:" before "Running mcmc loop:"; the four texts of slice A; nothing without `verbose`; no warning, counted with `withCallingHandlers` | `bart` unforced; the line missing or after the loop |
| line counts | moved, merged and reset equal counts made by hand from `getTrees` | merged counted as moved |
| same grid exact | twins on the donor's data, one forced and one not: identical trees and next 5 draws | the same-grid donor sent through the across-grids routine |
| sparse and mixed stores | forced across grids on a sparse matrix and a mixed frame: grid and codes identical after | the predictor pointer read where it is null |

tests/cpp: `declined` leaves each chain's state and generator bytes identical; the report's counts on a
hand-built donor; the 4 rewritten expectations; the cross-grid merged leaf.

## Slice B: what does not match in kind

### Steps

1. Tree count, both forests (about 230 lines, sampler.hpp and chain.hpp).
   - The donor forest's own count sets the stride of a kept draw; the two count refusals go.
   - More trees: the assembled chain is the donor's trees in order, then single leaves at 0 (slopes 0; gp
     fits 0; 1 in the variance forest). Clean.
   - Fewer trees: not clean. Forced, the donor's trees are tried in the recipient's live slots, the
     recipient's count at a time, by the install of slice A on either arm, the last round filled with single
     leaves.
   - After each round `liveTreeVariation` reads each tree's variance over the rows under the working
     weights, from its fits as installed ([`muByTree`](../../src/bartcore/combiner.hpp) through the leaf
     map, or the fit slab); `liveVarianceTreeVariation` reads the same of the log of
     [`factorByTree`](../../src/bartcore/chain.hpp).
   - Each chain keeps the trees that vary most, a tie to the earlier tree, in the donor's order, and the
     chosen trees are installed once more. Any failure restores every chain. The mean and the variance
     forest choose separately; each chain chooses among the trees of the donor sample it takes.
   - The report gains the single-leaf count, and `offered` is the donor's count.
2. Parts one side has (about 70 lines). The DART refusal and the variance presence refusal go. A donor's
   DART state or variance trees with nowhere to go are not assembled. A DART recipient of a plain donor, and
   a variance forest whose donor has none, are put back to their creation values (uniform split weights, the
   starting concentration and a zero delay counter; root trees at the starting variance). The report names
   each.
3. Leaf kinds and the errors of (d) (about 290 lines, the bridge and chain.hpp).
   - [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads each donor forest's leaf blocks in the
     donor's own shape: no `tree.params` is constant; `leaf.lengthscales` is gp, its row count the block's
     length over the tree count; otherwise linear, its covariate count the block's length over the leaf
     count. The donor's shape goes to the engine beside the state; the state's format does not change.
   - In the engine: constant into linear gives every leaf slopes of 0; linear into constant drops the slopes;
     linear into another number of covariates keeps the intercepts with slopes of 0 and skips the
     restatement; the last two are not clean.
   - Errors in both forms: "warm-start donor has gp leaves and this sampler constant ones" for each pair of
     kinds with gp in it; "warm-start donor's gp leaves were fit to 300 rows and this sampler holds 200";
     "warm-start donor was fit to 4 predictors and this sampler holds 3"; and the factor error naming the
     column by position, and by name where the bridge holds one.
4. Existing tests, by search: ["shape-compatible"](../../inst/tinytest/test-warm-start.R) has 3 hits in 2
   files, the DART refusal 1 in test-warm-start.R, the variance count 1 and the leaf parameters 1; in
   tests/cpp 2 `shapeMismatch` and 1 `varianceMismatch` expectations are read. Each becomes a check of what
   is taken.
5. New tests (about 500 lines; Tests).
6. Help, NEWS and notes (about 130 lines). The two "must share ... tree count, and DART setting" sentences
   and "is refused rather than silently reshaped" give way to the rule's last five rows; NEWS gains "A
   previous fit with another number of trees, with constant where the new fit has linear leaves or the
   reverse, or with DART or a variance forest on one side only is taken in part."; the design note's Status;
   TODO drops `warm-start-salvage` and `warm-start-tree-count`.

Commits: (1) step 1 with its rewritten tests; (2) step 2; (3) step 3; (4) step 5; (5) step 6. Commit 1
stands alone and is the split point named under Budget and stops.

### Tests

tinytest in a new `test-warm-start-kinds.R`, tests/cpp in `test_state.cpp`.

| Test | Pins | Mutant |
|---|---|---|
| more trees equal fit | 50 into 200: TRUE; the starting fit is the donor's; 150 single leaves at 0 | leaves not 0 |
| fewer trees declines | 200 into 50: FALSE, untouched | TRUE for a donor with trees left out |
| most varying, by hand | a hand-built donor of 6 trees into 3: the three of largest variance over the rows, by hand | the first three |
| donor order kept | the chosen trees sit in the donor's order with their values | sorted by variance |
| ties to the earlier | two trees of equal variance, one slot | the later tree |
| weights in the ranking | case weights under which the ranking differs from the unweighted one | unit weights |
| ranked as installed | an interaction-limited recipient where the donor's most varying tree is cut back to one leaf: not chosen | ranked as stored |
| each chain its own | two chains given two kept draws: each chooses among its own draw's trees | chain 1's choice for both |
| kept draw of another count | `samples` naming a kept draw of a 200-tree donor: that draw's trees | the recipient's count as the stride |
| variance forest, log scale | hand trees with factors (0.25, 1) and (1, 2.2) on equal halves: the first is chosen | variance of the raw factor |
| variance padding | 4 variance trees into 6: the variance function is the donor's | leaves not 1 |
| parts | plain into DART and DART into plain: taken, weights uniform or left; variance on one side: taken, the forest at its creation values on a sampler that has run | the error kept; donor weights installed; the forest left as it was |
| constant into linear | TRUE; the starting fit is the donor's; slopes 0 | slopes not 0 |
| linear into constant | FALSE; forced, leaf values the donor's intercepts | TRUE unforced |
| linear, other covariates | one covariate into two: FALSE; forced, intercepts kept and slopes 0 | the error kept |
| gp stays an error | gp against constant and against linear, each way, and gp on another number of rows: error in both forms, by name | taken under force |
| predictors and kinds | the two texts of (d), sampler untouched | a write before the error |
| the phrases | each phrase of the line by its case | a phrase on the wrong case |

## Verification

- Each slice: `cd tests/cpp && make && ./test_bartcore`, and again under
  `OPT="-O2 -g -fsanitize=address,undefined"`; `R CMD INSTALL --preclean -l <lib> .` and the full tinytest
  suite with `R_LIBS=<lib>`; the R-loaded sanitizer run over the slice's new test file.
- The neutral claims, each with its oracle:
  - No warm start. On the reference build the four `test-reproducibility-*.R` files pass unedited, and the
    compares against equivalence-308c45be.rds (`--bitwise --strict-coverage`), bcf-equivalence-1b7d730c.rds
    and multinomial-equivalence-80b1c8d4.rds print 55, 15 and 11 "identical draws" lines and no "max |z|"
    line. The 42 exact gates of [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) run in quick,
    expected unchanged. None of these scripts calls `installTrees` or `warm.start`: of benchmarks/R only
    composition-matrix.R does, and it is no gate.
  - A same-grid forced install is today's. A script under scratch/wsp is run on scratch/libs/sgi with
    `installTrees(donor)` and on the slice's library with `forceUpdate = TRUE`, and the saved results are
    compared with `identical`: trees, stored state and the next 5 draws for a gaussian pair on the same
    rows, a 40-row recipient on the donor's cut points, a monotone recipient, a heteroscedastic pair, linear
    leaves under another standardization, a recipient with missing values, a kept draw by `samples`, and two
    chains. A recipient on a finer grid is run too and reported: equal is expected, and a difference there
    makes that sequence shifting, not a stop.
  - `stateIsValid` is unchanged: tests/cpp's state tests, the 51 tinytest files that call `setState`, and a
    `setState` battery (clean, declined, forced repair) in the same cross-build script.
  - Slice B, same shape: the same script on slice A's library and slice B's.
- The shifting sequence has two oracles in the suite: "equals the two calls" above, and the equality with a
  forced `setState` that ["equals the two calls"](../../inst/tinytest/test-state-grid-install.R) already
  holds to the by-hand sequence.
- Mutants, slice A: each row of its Tests table names one; the reviewer builds each on its own library with
  `--preclean` and names the failing check. Not killable by any test: the restore on a throw in
  `installAcrossGrids` (TODO `warm-start-undo-test`).
- Mutants, slice B: each row of its Tests table.
- Lint, air, the codoc, win-drift and doc-freshness checks; `R CMD check --as-cran` from a clean tarball;
  the NEWS parse.
- No bench-sampler run is owed: no sweep path is edited. `stateIsValid`'s loops move without change.

## Budget and stops

| | Slice A | Slice B |
|---|---|---|
| engine | about 300 added, 290 removed | about 380 |
| bridge and R | about 150 | about 210 |
| tests | about 750 | about 800 |
| help, NEWS, notes | about 170 | about 130 |

Plan estimates have run 1.5 to 2 times low in tests (the last state slice landed at 2400 lines against
1035), so the upper figures are 2200 and 2500. Stop and report, building nothing further:

- if a bitwise compare prints a "max |z|" line or a snapshot file moves;
- if the cross-build script differs on any same-grid fixture;
- if a forced install across grids and the by-hand twin disagree on any tree or leaf value;
- if the tracer shows a call returning FALSE in a file outside the ten named in slice A's step 4 (list it,
  do not edit it);
- if a test outside those named fails because it expected a refusal that is now a repair;
- slice A at 450 added engine lines or 2200 in all;
- slice B at 1100 engine and bridge lines or 2500 in all: commit 1 then lands alone as the tree-count
  slice and the rest follows as a third slice.

## Interactions

- state-grid-install ([state-grid-install.md](state-grid-install.md)): lands first; this plan calls its
  routine and retires the remap arm it left to this item.
- The forest slices (forest-sd-unit, forest-kind-by-class, forest-multiplier) and collapsed-scale-move:
  chain.hpp, the bridge, R/dbarts.R and the Rd, in other functions. Rebase only. A recipient of two or more
  forests stays refused, so nothing of theirs is installed here.
- The record slices of TODO `state-frame-prior`: with the grid on the data object the across-grids routine
  reads the donor's grid as it does now.
- Consumers, read with `git grep` and not run: stan4bart (branch bartcore) calls `setState(forceUpdate =
  TRUE)` and neither `installTrees` nor `warm.start`; bartCause and treatSens call none of the three;
  bairrtt calls `setPredictor` only. None sees either slice. No reverse dependency can call a warm start
  that 0.9-34 did not have.

## Calls made in planning

- A forced warm start across grids keeps the splits the across-grids routine keeps: a leaf that none of the
  recipient's rows reaches under the donor's own thresholds is merged before any split moves, where today's
  arm moves first. Measured above: 173.8 of 196 nodes kept at 30 rows against today's 184.0, the same within
  one node from 100 rows on, the starting fit slightly nearer the donor's function (0.675 against 0.694).
  One routine serves `setState` and the warm start, and the rows a merged leaf is weighed by are these. The
  alternative, a warm start step that keeps unreached leaves through the first pass, is about 90 lines and
  a second rule, and a leaf with no row has no value to recover from its fits.
- A recipient column held as an ordered factor against a donor column on another grid (a numeric column,
  or an ordered factor of another level count) is an error in both forms, where today it is taken by
  position (9 of 11 splits kept in the run): the routine cannot re-quantize a factor column, and `setState`
  refuses the same state. The alternative keeps the remap arm for such columns, about 90 lines, and the
  plain mean with it. One sequence that installs today becomes an error.
- Under `blocks`, where a tree's allowed columns depend on its slot, a donor of more trees is ranked as cut
  back for the slot it is tried in and installed in the slot its rank gives. The two can differ for a tree
  near the cut, and the installed forest is valid either way. The alternative refuses another tree count
  under `blocks`.
- Two slices, with a third planned for only at slice B's stop. One slice would be about 2900 planned lines
  under one review; three from the start would add a plan, a critique and a review to a serial queue.
- A forced install returns NULL invisibly whether or not anything changed: choice (b) read through
  install-surface.md, as asked.
- The verdict is `stateIsValid`'s staging, moved into helpers, and not a new routine. The alternative, a
  second staging for donors, is what the state slices removed for states.
- The cross-grid arm is the across-grids routine handed today's same-grid installs; the remap arm and the
  four scratch passes are deleted, about 290 lines.
- The predictor matrix is a required argument of the C++ `installForests`: a default of null would be read
  on a dense column. The 44 calls in tests/cpp gain it.
- `force` defaults to true in C++, as `setState`'s does, and the bridge states it on every call.
- gp fits start at 0 whenever the donor's recorded points differ from the recipient's, as today, and such a
  donor is not clean: the starting fit is not the donor's.
- The counts are made by comparing flat trees after the install, so no install routine is instrumented.
  Moved and merged are disjoint.
- The line's words. They can change after release.
- `bart` calls the internal installer and not the method, since the method returns no report.
- After a forced install the `state` field is the recipient's stored state, as today; a donor's state is not
  the recipient's, so `setState`'s rule of keeping the state as given has nothing to keep.
- Fewer trees are ranked by trial installs in the live slots and not by a scratch staging: the tree ranked
  is the tree installed by construction, at the cost of one install per round.
- A part only the recipient has is put back to its creation values, dec-B324's "as a new fit does" read for
  a sampler that has run.
- The donor's leaf kind travels beside the state in memory; the stored format is untouched.
- The texts of (d) come with slice B, where the bridge's reader is rewritten once.

## Not verified

- How many of the 67 calls return FALSE under the default: read, not run. No suite and no tracer was run.
- 2 of the 49 family pairs were run.
- That the same-grid forced arm is today's bit for bit: argued from the code path. Nothing was built.
- That a donor on a finer grid installs bit for bit through the same-grid arm.
- That a linear leaf's cached statistics are invalidated on the across-grids arm as the remap arm
  invalidates them; "equals the two calls" on linear leaves is the check.
- `installAcrossGrids` and a throw from its install step: read; no test can reach it.
- The cost of the trial installs.
- The consumers: read, not run.
