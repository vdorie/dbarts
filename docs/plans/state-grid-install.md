# state-grid-install: a state goes onto the sampler's cut grid and does not replace it

Status: PLANNED 2026-10-11 (dec-B449, dec-B450; dec-B398 for the merged leaf); revised from its blind
critique 2026-10-11 ([Revision from critique](#revision-from-critique)). Nothing is built; no call is open.
The rule is a change from released 0.9-34, which stores the cut grid in its state and installs it, as the
branch does today.

agent: opus implementer, one; one opus reviewer who runs the mutants under [Tests](#tests). The blind
critique ran (scratch/sgp-crit.md). Tier "Changes draws" ([Process by risk](README.md#process-by-risk)):
the diff is in the engine's install path.
rng: POSTERIOR-CHANGING on one call sequence: a state installed into a live sampler whose cut grid, or its
weights, differ from what the state records. Afterwards the sampler holds its own grid and weights where
today it holds the state's, and both are part of the tree prior. NEUTRAL, bit for bit, for every install
onto the grid and weights the state was stored under, for every copy and reload, and for every fit that
never calls `setState`. No recorded scenario, snapshot file or exact gate reaches the changed path
([Verification](#verification)).
window: before the merge to main; engine slices stay serial. After response-scale-rows (landed on this
base) and before the record slices of TODO `state-frame-prior`. Shared files: [Interactions](#interactions).
budget: planned ~1475 lines (engine ~250, bridge and R ~85, tests ~1000, help, NEWS and notes ~140); upper
figure 2200. Stops: [Budget and stops](#budget-and-stops).

## Goal

`setState` on a live sampler leaves the sampler's cut grid and the grid's weights alone and puts the state's
trees onto that grid. A state whose every split threshold is a point of the sampler's grid, and whose trees
fit as they must today, goes in as stored and the call returns TRUE. One with a threshold the grid does not
hold is declined unforced (FALSE, the sampler as it was) and, forced, has that split moved to the nearest
point its tree still allows, or merged with what hangs under it into one leaf at the mean of the stored
leaves, weighted by the rows each held (NULL). The state keeps recording the grid it was stored under; a
copy and a reload, which re-create the engine and have no other record of a grid, still take it from
there. What 0.9-34's restore did is reached by hand, `setCutPoints(attr(state, "cutPoints"))` and then
`setState(state)`, which the help states and a test pins. Closes TODO `state-grid-replaces-sampler` and the
grid half of TODO `state-frame-prior`.

## Context

Rulings, in the register's words unless quoted as the maintainer's:

- dec-B449, the maintainer on 2026-10-11, shown that 0.9-34 restores the grid as the branch does: "Change
  it, but investigate for me first if there is a way to accomplish the old behavior by manually sequencing
  some function calls." The rule: "a live setState never writes the sampler's cut points or their weights";
  "the help states the sequence".
- dec-B450, the maintainer the same day: "And regardless, defer adding weights to setCutPoints until after
  the main merge." So no `probs = ` here; "a restore onto a grid with the stored points and other weights
  goes through, the chain running under the sampler's weights". dec-A194: `setCutPoints` "keeps a column's
  weights only when handed the grid it already holds, byte for byte, and installs any other grid
  unweighted".
- dec-B231, the maintainer on 2026-10-05: "the predictors may change without changing the cut points. That
  actually makes them a form of 'data', albeit one that is initially derived from the user input." The
  entry concludes "They therefore belong on the data object, and a saved state does not install them"; the
  second half is the entry's, not a quote. dec-B195, the maintainer on 2026-10-02: "The way I see it,
  state (and scratch) describe the current literal state of an MCMC sampler as it walks the posterior, as
  defined by a model and data."
- dec-B234 gave `setState` its status and lists what makes an install inexact: "leaves no row reaches
  merged, splits moved to the nearest current split point, or a missing-value direction dropped, which
  happen when the sampler's predictors or split points changed after the state was stored". dec-B305, the
  maintainer: "Cleanly was defined as leaving a sampler in a coherent state with regards to its model";
  the entry's rule names "a split off the grid" as not clean. dec-B310: forced forms return NULL
  invisibly. dec-B305: "copy() and a reload force, having no choice."
- dec-B398, shown a forced install that "collapses the first offending split on each path from the root
  into one leaf at the weighted mean of the leaves beneath it", the maintainer on 2026-10-09: "Repair, as
  we do elsewhere." It rules the value of a leaf this plan's move makes.
- dec-B300: a state whose recorded grid repeats a point is refused, forced or not. Stays.
- dec-B311 and dec-B312 give a replaced grid's splits "position" or "value" at `setPredictor` and
  `setCutPoints`. They name those two calls, not `setState` ([Calls made in planning](#calls-made-in-planning)).
- dec-B418: the help already tells a caller to put a re-anchored response back before `setState`; this
  plan says the same of a grid. dec-B230: a state from a sampler on other data is "neither supported nor
  prevented".

Code as it stands:

- A split is stored by value. [`FlatNode`](../../src/bartcore/tree.hpp) holds an ordinal rule's cut point,
  and [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) puts it on the grid position holding exactly
  that value, failing where none does. So a state's trees resolve only on a grid that holds their
  thresholds.
- [`Sampler::setState`](../../src/bartcore/sampler.hpp) therefore installs the state's grid first
  ([`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp), re-quantizing the columns from the
  raw predictors R hands in, and the state's [`cutMass`](../../src/bartcore/data.hpp) even where the points
  are equal), judges and installs the trees against it, and puts the old grid back only on a refusal or a
  decline. [`SamplerStateData`](../../src/bartcore/sampler.hpp) carries `cutPoints` and `cutMass`; the
  bridge writes them as attributes of the state and requires the first
  ([`bartcore_setState`](../../src/R_interface_bartcore.cpp),
  [`readStateCutMass`](../../src/R_interface_bartcore.cpp)).
- Kept draws (`keepTrees`) are flat trees too, are copied by an install and never built against a grid:
  [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) holds them to form only
  ([`flatTreeIsWellFormed`](../../src/bartcore/tree.hpp)), and `predict` routes raw values through them.
- `setCutPoints(splits = "value")` moves a live tree's splits onto a new grid:
  [`Sampler::setCutPoints`](../../src/bartcore/sampler.hpp) installs the grid and ends in
  [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp), which moves each split by value
  ([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp),
  [`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp): the nearest point inside the interval the
  splits above leave, an exact match kept, a tie going to the upper point), merges a split with no point
  left at the mean of the leaves beneath weighted by the rows each holds (`LiveRows`), re-routes the rows
  and merges a leaf left empty; mean forests of every leaf kind, the variance forest and the monotone
  reseed included.
- The warm start moves a donor's trees too ([`Sampler::installForests`](../../src/bartcore/sampler.hpp),
  [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp)), but builds each on the donor's grid
  with no rows under it ([`ScopedCutGrid`](../../src/bartcore/data.hpp)), so its merge has no row sums and
  takes the plain mean.
- R re-creates an engine from `(control, model, data)` in three places and then installs a state
  ([`installStateOnto`](../../R/dbarts.R)): [`getPointer`](../../R/dbarts.R) after a reload, `copy`, and
  `setState` on a dead pointer ([`recreatePointer`](../../R/dbarts.R)). The data object holds no grid, so
  the fresh engine's grid is the one a creation derives from the predictors as they now are; the stored
  state is the only record of a grid set by `setCutPoints` or held while a column moved.
- Nothing documented reads the grid in force. Tests read `attr(sampler$state, "cutPoints")`.

### Measured: 0.9-34 beside the branch

Scripts in scratch/sgp (p1 to p5), run on scratch/libs/cran (0.9-34) and scratch/libs/rsrl (d5b637ee,
shipped build). 200 rows, two numeric columns, 20 trees, 100 points per column at creation; "3 points" is
`setCutPoints(c(0.25, 0.5, 0.75), 1)`. A 0.9-34 state was copied with `unserialize(serialize())` (below).

| Case | 0.9-34 | d5b637ee |
|---|---|---|
| what a state carries | attributes `numCuts` and `cutPoints`, one grid per column | attributes `cutPoints` and `cutMass` |
| grid unchanged, restore | NULL; splits as stored | TRUE; splits as stored |
| 3 points set, earlier state restored | NULL; 100 points again, the stored grid; splits as stored | TRUE; the same |
| column 1 replaced and its grid refreshed, restore without putting the column back | the stored grid back (0.023 to 0.98, from 0.030 to 2.93) | the same |
| kept draws: `predict` before and after 3 points, no restore | unchanged, with a warning | unchanged |
| kept draws: restore after 3 points | 100 again; predictions unchanged | the same |
| `copy()` and reload, state stored after the change | 100 values on column 1: the 3 set, then 97 read past the stored vector (`numCuts` is not rewritten by a store) | 3 and 100 |
| `copy()` and reload, state stored before the change | 100 and 100, the stored state's | the same |
| two chains; a sampler given 3 points takes another sampler's state stored on 100 | 100 again | 100 again |

Three more facts of 0.9-34. Its `restoreState` reads the trees by index and then installs the state's grid
with the engine's own `setCutPoints`, on every restore (source read on main, the 0.9-34 tree). Its
`storeState` writes into the R object the sampler already holds, so a state kept by plain assignment is
overwritten by the next store: with `updateState = TRUE` on the changing call the kept state held 3 points
and a restore brought nothing back; kept by assignment with `updateState = FALSE`, or copied, it held 100
and the restore brought them back (p5). Its `setState` returns NULL and leaves the `state` field NULL.

So the branch does what the release did. The rule below is a change from the release, made for dec-B231
and ruled as dec-B449 against leaving the restore as 0.9-34 has it; it repairs no regression. TODO
`state-frame-prior` still says "0.9-34, whose state carried no grid"; step 6 corrects it.

### Measured: what the rule would do

On d5b637ee. The forced move planned here is two calls that exist: the state installed on its own grid,
then the sampler's grid put back with `setCutPoints(grid, splits = "value")` (p7; the critique's c1.R).

- One tree, root on the 50th of 100 points, left child on the 30th, leaves 0.10, -0.20 and 0.05 holding 73,
  38 and 89 rows, onto the grid of the 50th and 75th points: the child has no point left and the two
  leaves become one at -0.0027, their mean weighted by 73 and 38. The warm start gives -0.05, the plain
  mean; moving the child to the nearest point of the whole column and merging by today's install gives
  0.10. The mirror (right child, leaves holding 21 and 68 rows): -0.129. Two deep (left child and its left
  child both without a point; leaves 0.3, 0.1, -0.2 holding 17, 56, 38): one leaf at 0.0279.
- A sampler given 3 points takes a 20-tree state stored on 100: it keeps 3 points; the state's 17 splits on
  column 1 become 15, on 0.25, 0.5 and 0.75, and 29 splits in all become 27.
- A sampler given 199 points holding the stored 100 takes the trees with every split as stored (p4).
- A column replaced with `updateCutPoints = "position"` and replaced back with its old values, then
  `setState(stored)`: the grid is the stored one bit for bit, TRUE, and the next five draws are those of a
  twin that only restored. With the old values put back and no refresh, none of the state's 19 splits on
  the column is a point of the grid in force; today `setState` returns TRUE and brings the old grid (p4).
- By hand, the orchestrator's runs for dec-B449: `setCutPoints(attr(state, "cutPoints"))` then
  `setState(state)` gives the stored grid back identically, TRUE, and the next 5 draws of a direct
  restore, on plain numeric columns, under case weights, beside a factor column and on few-valued columns.
  The whole-set call skips a factor's and an ordered factor's entry; column by column it errors on them.
- A weighted grid (p6: an integer column of 30 distinct values, 29 points, each with a weight). Half 1,
  `setCutPoints` with the stored points: the points are back and the column has no weights; today
  `setState` then returns TRUE, brings the state's weights and draws as a direct restore does. Half 2,
  `setPredictor(x[, j], j, updateCutPoints = "value")` with the column's own values: points and weights
  back identically, TRUE, the same draws. A weighted column handed back the points it holds keeps its
  weights. Today a state with the sampler's points and no `cutMass` returns TRUE and strips the weights.
- Two samplers made on the same data hold identical grids and weights, which is stan4bart's restore.
- Kinds (p4; the critique's c4b.R). A numeric column's state into a sampler whose column is an unordered
  factor: "state is not consistent with this sampler", forced or not. Into one whose column is an ordered
  factor: the same error today, but only because its thresholds are not level midpoints. An ordered
  factor's state into a numeric sampler: TRUE today, the numeric grid replaced by the 3 midpoints.
- The critique's tracer over the 51 installed test files that call `setState`: 748 live installs, 19 with
  a recorded grid or weights other than the sampler's, in six files.

## The rule

On a live sampler (an engine that existed before the call), the first row that applies:

| The state against the sampler | Unforced | Forced |
|---|---|---|
| malformed: a live threshold that is no point of the grid the state itself records; a recorded grid that repeats a point; a column an unordered factor on one side only; on a column the sampler holds as an ordered factor, a recorded grid other than the sampler's | error, sampler untouched | error, sampler untouched |
| its recorded points and weights equal the sampler's on every column | as today, bit for bit | as today |
| the points equal, the weights differ on a column | as today on the trees (TRUE where they fit); the sampler's weights untouched | NULL; the same |
| the points differ; every live threshold is a point of the sampler's grid | as today on the trees (TRUE where they fit, FALSE for a leaf no row reaches and the like); the sampler's grid and weights untouched | NULL; today's repairs |
| a live threshold is not a point of the sampler's grid | FALSE; nothing touched | NULL; the move below |

The forced move, in three passes over each live tree, mean and variance (dec-B398):

1. The tree is read as stored. The sampler's rows are routed through it by its stored thresholds, so each
   stored leaf holds the rows the stored tree sends it, counted with their case weights. For a state from
   another sampler's data these are the receiving sampler's rows. A stored leaf no row reaches is merged
   here by today's repair.
2. Each split on a column whose grid differs goes to the nearest point of the sampler's grid inside the
   interval the splits above it leave, a tie to the upper point. A split with no point left becomes, with
   everything beneath it, one leaf at the mean of the leaves beneath weighted by the rows each held in
   pass 1; the plain mean where none holds a row.
3. The rows are routed again by the moved thresholds, and a leaf left without a row is merged by today's
   repair, as after a forced `setPredictor`.

That is what installing the state on its own grid and then calling `setCutPoints(grid, splits = "value")`
gives today, with the sampler's weights kept.

- The weights of a grid ([`cutMass`](../../src/bartcore/data.hpp)) are part of the grid: no live install
  writes them, whatever the state records. Other weights are clean: the rule asks whether every live
  threshold is a point of the sampler's grid, weights place no threshold, no tree is changed, and the
  sampler is left "in a coherent state with regards to its model" (dec-B305), as it is on a finer grid.
  TRUE then says the trees are the stored ones, not that the chain continues as the stored chain would:
  it "runs under the sampler's weights" (dec-B450), and its next draws are other draws.
- A column's kind. The state records none. An unordered factor has an empty entry, which tells it from
  every other column. An ordered factor's grid is fixed by its level count, so where the sampler holds
  one, any other recorded grid is another column and is refused. The other way round cannot be told: an
  ordered factor's entry (0.5, 1.5, 2.5) is also what a numeric column of four whole numbers records, so
  such a state into a numeric column follows the numeric rule, FALSE unforced where a threshold is off
  the grid. No kind record is added (dec-B283: "let setState do what it can with no guarantees").
- A rule's missing-value direction rides the rule through a move. The draw for a column that can be
  missing now and could not when the state was stored
  ([`Chain::completeMissingDirections`](../../src/bartcore/chain.hpp)) is made first, as today.
- A gp leaf's per-row fits are kept through a move, as a forced `setState` and `setCutPoints` keep them.
- Kept draws go in as stored, never moved: they predict what they predicted. A kept draw can therefore
  hold a threshold the grid no longer has, as it can today after `setCutPoints`.
- Several chains share one grid; one chain's tree that needs a move declines the whole state.
- The `state` field. After a live install onto a grid or weights other than the state records, forced or
  not, the field holds what `storeState()` gives, so a `copy` and a reload reach the sampler's grid,
  weights and trees. After an install onto the recorded grid and weights it holds the state as given, as
  today.
- `installTrees` and `bart(warm.start = )` are not changed here. The warm start's plain mean is a
  difference from this rule that TODO `warm-start-salvage` removes by switching its cross-grid arm to the
  routine of step 2, which is built here, lands first, and then serves both doors.

At a re-creation (`copy`; a reload's first use) the new engine takes the grid and weights the state
records and then installs the state, as today. `setState(newState)` on a reloaded sampler first builds the
engine from the sampler's own saved state, where it has one that installs, and then applies the live rule
to `newState`: the call gives the same grid whether or not the engine was live. Only a sampler with no
saved state, or one that is refused, takes `newState`'s grid, as today. A sampler copied or saved after a
grid change with no store since is re-created on the stored state's grid, as today.

Undoing a grid change is the caller's, before `setState` (dec-B449), with the state held in a variable
(`setCutPoints` replaces the `state` field under `updateState`): `setCutPoints(attr(st, "cutPoints"))`,
the whole set in one call, and then `setState(st)`. Where no changed column's stored grid carries weights
that is the stored grid and the stored chain, bit for bit, what the branch's `setState` alone gave before
this rule. 0.9-34's own restore brought back the grid and the splits and not the leaf values (the review's
run: leaf values off by up to 0.13). A column whose grid
was changed and whose stored grid carries weights (`attr(st, "cutMass")` holds an entry for it) gets its
points and not its weights, so the restore is TRUE and the chain runs unweighted on that column;
`setPredictor(x[, j], j, updateCutPoints = "value")` gives back both while the column holds the values the
grid was derived from, and nothing does once it holds others (TODO `cut-point-weights`, after the merge).

## Where the line falls

In this slice: what a live `setState` does with the grid, its weights and the splits. Not in it, and left
to the record slices of TODO `state-frame-prior`: the `cutPoints`, `standardization` and `response.range`
records on dbartsData (dec-B231 to dec-B233), `setData`'s three arguments and its carried-record rule, a
creation that uses a carried record, a linear or gp leaf's standardization and slopes at an install, the
manual's table of what puts each derived value back, and the gp `setData` refusal.

The grid half stands alone on a live sampler: the trees carry their thresholds, the move is built for
`setCutPoints`, and nothing has to be recorded anywhere new. It does not stand alone at a re-creation. There the
state is the only record of a grid that a creation would not derive, so a copy and a reload keep taking
the grid from the state until dbartsData records it; taking that away now would make a copy of a sampler
with a set grid something other than a copy.

## The state's format and its consumers

No change. `cutPoints` stays required and `cutMass` optional, both written as today;
[`stateFormatVersion`](../../src/R_interface_bartcore.cpp) and the floor stay at 1. `cutPoints` is now the
record of the grid the thresholds were stored on: compared on a live install, adopted at a re-creation.
A state stored before this slice is the same object and follows the same rule. No flat entry of
[dbarts.h](../../inst/include/dbarts/dbarts.h) stores or installs a state, so the C API hash does not move.

stan4bart stores one state per chain and restores each into a sampler it has just made from the fit's own
`(control, model, data)`, with `setState(state, forceUpdate = TRUE)` (its `restoreBartSampler`). It never
sets cut points or refreshes a grid, and a sampler made on the same data derives the same grid and weights
(run above), so its installs are same-grid installs and see nothing. bartCause and treatSens call neither
`setState` nor `setCutPoints`; bairrtt calls `setCutPoints` and never `setState`. No vignette restores a
state after a grid change.

## Constraints

- An install onto the recorded points and weights is the code that runs today. The proof is the existing
  suite, about 730 such live restores in 51 test files, and the "same grid exact" test; the equivalence
  compares and the snapshot files never reach `setState` and are run as the class's battery only.
- [`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp), [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp)
  and the warm start's routines are called or left, not edited.
- The unforced answer is FALSE exactly when the forced install changes a tree: a threshold off the grid, or
  what today's verdict reports ([`Tree::stageFromFlat`](../../src/bartcore/tree.hpp)).
- A declined or refused live install leaves the grid, weights, codes, chains, generators and kept draws as
  they were. Judging a state on its recorded grid swaps the points over the store for the build only
  ([`ScopedCutGrid`](../../src/bartcore/data.hpp)), as the warm start's checks do.
- `adoptGrid` defaults to true in C++, in [`Sampler::setState`](../../src/bartcore/sampler.hpp) and the
  [`SamplerBase`](../../src/bartcore/facade.hpp) virtual: the 89 calls in tests/cpp keep today's path, and
  the bridge states the flag on every call.
- `setState` stays refused on a borrowed view: a forced move and a re-creation re-quantize.
- No new argument, no message, no warning (dec-B234: "Return a status, no warning").
- Out of scope: everything under [Where the line falls](#where-the-line-falls); a reader for the grid in
  force; weights on `setCutPoints`; `installTrees`'s `forceUpdate` (TODO `warm-start-salvage`).

## Steps

1. Thresholds on the grid (~60 lines, [tree.hpp](../../src/bartcore/tree.hpp),
   [chain.hpp](../../src/bartcore/chain.hpp)). Whether every ordinal threshold of a chain state's live mean
   and variance trees is a point of a given grid for its column: exact equality, by binary search over the
   flat nodes, no tree built. Kept draws are not read. It is asked of the recorded grid, where a miss is
   the refusal, and then of the sampler's.
2. The install (~190 lines, [sampler.hpp](../../src/bartcore/sampler.hpp),
   [facade.hpp](../../src/bartcore/facade.hpp)). [`Sampler::setState`](../../src/bartcore/sampler.hpp)
   gains `adoptGrid` and reports whether the grid or weights it went onto are other than the state's.
   - The ordered-factor check joins the checks at the top.
   - `adoptGrid` true: today's body.
   - False, equal points: today's body with no grid step and no write of the weights.
   - False, other points, step 1 true for every chain: the same, the state judged and installed on the
     sampler's grid, where every threshold resolves.
   - False, a threshold off the sampler's grid: every chain is judged on the recorded grid, swapped in
     for the builds; a state that does not build is the refusal. Unforced, the decline. Forced, one private
     routine (working name `installAcrossGrids`): snapshot the sampler's points, counts, weights and codes;
     install the recorded grid with re-quantized codes, as today's body does; run the install step it is
     handed, here today's chain installs, which build and partition each tree on that grid (pass 1); move
     the snapshot back, byte for byte; call each chain's
     [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp) with the recorded grids of the columns
     that differ and `byValue` (passes 2 and 3). The install step is a parameter so that the warm start
     can hand in its own.
   The virtual takes both: install with `--preclean`. Rewrite the comment blocks above `setState` and on
   [`SamplerStateData`](../../src/bartcore/sampler.hpp).
3. Bridge and R (~85 lines). [`bartcore_setState`](../../src/R_interface_bartcore.cpp) takes the flag as
   a sixth argument (its row in [R_interface.cpp](../../src/R_interface.cpp)), returns declined, installed,
   or installed onto another grid or weights, and has its comment on views corrected.
   [`installStateOnto`](../../R/dbarts.R) passes the flag, TRUE from `getPointer` and `copy`. `setState`
   on a dead pointer re-creates as `getPointer` does where the sampler has a saved state that installs,
   binding that engine, and then installs live; otherwise it re-creates bare and adopts, as today. After a
   live install onto another grid or weights it caches the engine's state, not the one given. Docstrings
   of `setState`, `copy` and `storeState`.
4. Tests ([Tests](#tests); ~1000 lines). Six tinytest cases pin today's rule and are rewritten, not
   deleted (the first three are named by wording the rewrite has since replaced):
   retired: ["a state stored on the longer grid installs"](../../inst/tinytest/test-cut-points-undo.R) and
   its neighbor, retired: ["a state without them installs the grids unweighted"](../../inst/tinytest/test-cut-grid-distinct.R),
   retired: ["the setState cut-point replay"](../../inst/tinytest/test-data-mixed-mutation.R) and its sparse twin
   in [test-mutate-sparse-valued.R](../../inst/tinytest/test-mutate-sparse-valued.R) (both move to a
   reload, which still re-quantizes),
   ["a state on another cut grid"](../../inst/tinytest/test-setstate-force-update.R) (its last line
   expects the state's grid) and
   ["sendsRight(other$state)"](../../inst/tinytest/test-state-missing-direction.R) (a state from data with
   missing values forced onto the filled data, most thresholds off its grid, then an equality of split
   columns that holds only where no move merges). In tests/cpp
   [`testRestoreStatus`](../../tests/cpp/test_state.cpp) and
   [`testWeightedGridState`](../../tests/cpp/test_state.cpp) stay as tests of the adopting path, and the
   live path gets its own.
5. Help, NEWS and notes (~140 lines). [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the
   `state` field names the `cutPoints` attribute beside `cutMass` (one entry per column, the points a
   state's thresholds were stored on) and says what the field holds after an install onto another grid;
   Saving (the list of what does not fit gains the threshold; the move and its merged leaf; "The cut grid
   ... come with the state as they are" becomes the rule; the reload sentence); `setCutPoints` and
   `updateCutPoints`. The sentences for `setState`, in Saving and in the method's docstring: a state does
   not bring its cut points; to restore them too, hold the state in a variable and call
   `setCutPoints(attr(st, "cutPoints"))`, the whole set, before `setState(st)`; a column whose grid was
   changed and that `attr(st, "cutMass")` lists then has its points and not its weights, the restore still
   returns TRUE and the chain runs with that column's points equally likely; `setPredictor(x[, j], j,
   updateCutPoints = "value")` restores both while the column holds the values its grid was derived from.
   NEWS: `setState` no longer brings back a cut grid, which 0.9-34 did, with the sequence.
   [install-surface.md](../design/install-surface.md) gains the section the posterior-changing class asks
   for and has its sentence on a copy after a repair kept true; one sentence each in
   [state-not-model.md](../design/state-not-model.md), [cut-grid.md](../design/cut-grid.md) and
   [architecture.md](../architecture.md). TODO: drop `state-grid-replaces-sampler`, cut the grid half from
   `state-frame-prior` and correct its sentence on 0.9-34; `cut-point-weights` stays as it is.

## Tests

tinytest in a new `test-state-grid-install.R`, tests/cpp in `test_state.cpp`. Each line: what it pins;
the mutant it must kill.

| Test | Pins | Mutant |
|---|---|---|
| live grid kept | 3 points set, an earlier state forced in: the sampler reports 3 points and its weights as before | `adoptGrid` true on a live call |
| same grid exact | store, 5 sweeps, restore: TRUE, printed trees as stored, next draws equal a twin's | the across-grids routine run on a same-grid state |
| finer grid is clean | 199 points holding the stored 100: TRUE, 199 kept, every split as stored | verdict keyed on the grids differing, not on a threshold |
| on the grid, a leaf empty | the finer grid and a column changed so a stored leaf holds no row: FALSE, as on the stored grid | step 1's TRUE returned as the verdict |
| off the grid declines | FALSE; grid, weights, codes, trees, `state` field and next draws equal a twin that never called | the verdict not reading step 1; a write before the decline |
| nearest inside the interval | hand-built one-tree states: 0.37 on {0.25, 0.5, 0.75} goes to 0.25; midway goes up; a child stays below its parent's new point | the move by position; the interval dropped |
| merged leaf, left, mirror, two deep | the three hand trees of the measured section: -0.0027, -0.129 and 0.0279, each the mean of the stored leaves by the rows each held | the plain mean (-0.05); rows counted after the move (0.10); a child left below an emptied node |
| equals the two calls | a 20-tree state forced across grids: trees, leaf values and the next 5 draws equal a twin's that installed it on the state's grid and then called `setCutPoints(grid, splits = "value")` | pass 1 skipped (trees moved with no rows under them) |
| another sampler's data | a state from a sampler on other rows: the merged leaf weighs by the receiving sampler's rows, by hand | weights read from the state |
| field after crossing | after a forced move, after TRUE on a finer grid and after the weighted by-hand half 1: `copy` and a reload report the sampler's grid, weights and trees; after a same-grid forced repair the field is the state as given | the state as given cached after a crossing |
| reload, first call | a sampler saved on 3 points, `setState(st100)` as the first call and after another call: the same grid and the same value; with no saved state, or a refused one, the handed state's grid | adopting whenever the engine was created in the call |
| weights are the sampler's | equal points: a state with no `cutMass` into a weighted sampler, and a weighted state into an unweighted one, each TRUE unforced and NULL forced, weights as they were, trees as stored | the weight write left on the same-points path; other weights made not clean |
| by hand, five fixtures | plain numeric columns, case weights, an unordered and an ordered factor column beside, few-valued unweighted columns: after 3 points, with the state in a variable, `setCutPoints(attr(st, "cutPoints"))` then `setState(st)` gives the stored grid identically, TRUE, and the next 5 draws of a twin that only restored | (pins the help's sequence); the same-grid path altered |
| by hand, weighted, half 1 | the 30-value integer column: `setCutPoints` with the stored points leaves no weights on the changed column and keeps them on an unchanged weighted one; `setState` is TRUE, the weights stay as they are, the trees are the stored ones | weights brought by the state; FALSE for other weights |
| by hand, weighted, half 2 | `setPredictor(x[, j], j, updateCutPoints = "value")` gives points and weights identically; `setState` TRUE; draws equal the twin's | (pins the help's weighted case) |
| re-creation adopts | stored after 3 points: `copy` and a reload report 3 points and the weights; mixed and sparse stores re-quantize | `adoptGrid` false in `getPointer` or `copy` |
| kept draws stay | `keepTrees`: forced across grids, `predict` is bit for bit what it was and kept thresholds are the old ones | kept draws moved |
| variance forest | a variance split off the grid alone gives FALSE; forced, it lies on the grid | variance trees not read by step 1 |
| kinds | a numeric column's state into an unordered-factor and into an ordered-factor sampler, and an ordered factor's of another level count: error in both forms, sampler untouched; an ordered factor's state into a numeric sampler: FALSE unforced; subset rules beside a moved column identical | the ordered-factor check dropped (FALSE and a move in place of the error) |
| neither grid | a threshold that is a point of neither grid: error under an equal and under another recorded grid | the decline reached before the build |
| direction rides the move | a moved rule keeps its missing side; a state from before a first missing value draws the sides a same-grid install draws | the draw made after the move |
| second chain | an off-grid split in chain 2 only: FALSE, chain 1 untouched | step 1 read for the first chain only |
| linear and gp leaves | forced across grids: unmerged linear leaves keep intercept and slopes; a gp forest's per-row fits are the stored ones | fits zeroed as a warm start zeros them |
| refresh away and back | a column replaced under "position" and replaced back, restore: grid identical, TRUE, draws equal a twin that only restored | (pins the help) |

## Verification

- `cd tests/cpp && make && ./test_bartcore`; again under `OPT="-O2 -g -fsanitize=address,undefined"`.
- `R CMD INSTALL --preclean -l <lib> .`, then the full tinytest suite with `R_LIBS=<lib>`; the R-loaded
  sanitizer run over the new test file.
- Reference build, the class's battery: the four `test-reproducibility-*.R` files pass unedited, and
  `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-308c45be.rds --bitwise
  --strict-coverage`, `bcf-equivalence.R compare .../bcf-equivalence-1b7d730c.rds --bitwise` and
  `multinomial-equivalence.R compare .../multinomial-equivalence-80b1c8d4.rds --bitwise` print 55, 15 and
  11 "identical draws" lines and no "max |z|" line. None of them calls `setState`, `setCutPoints` or
  `copy`, or reloads a sampler.
- The 42 exact gates in quick, as the class requires; expected unchanged: only negbin-mixing calls
  `setState`, on its own sampler's state.
- Lint, air, the codoc, win-drift and doc-freshness checks; `R CMD check --as-cran` from a clean tarball.
- stan4bart's suite against the slice's library, expected unchanged.

## Budget and stops

| Step | Lines |
|---|---|
| 1 tree.hpp, chain.hpp | ~60 |
| 2 sampler.hpp, facade.hpp | ~190 |
| 3 bridge and R | ~85 |
| 4 tests | ~1000 |
| 5 help, NEWS, notes, TODO | ~140 |

Stop and report, building nothing further: if a bitwise compare prints a "max |z|" line or a snapshot file
moves; if the forced move and the twin's two calls disagree on any tree or leaf value; if a tinytest file
outside step 4's six fails because it expected a state's grid or weights to arrive (list it, do not edit
it); at 375 engine lines or 2200 in all.

## Interactions

- response-scale-rows, landed on this base: no shared routine; its slice C shares the Saving section of
  the Rd (a sentence that `setState` does not undo a re-derivation).
- The forest slices (forest-sd-unit, forest-kind-by-class, forest-multiplier): chain.hpp, combiner.hpp,
  the bridge, [R/dbarts.R](../../R/dbarts.R) and the Rd, none in `Sampler::setState`. Rebase only.
- collapsed-scale-move (branch wt/collapsed-move-plan): chain.hpp and model.hpp; it names
  `Chain::setState` for the generator a state carries. This slice adds one routine to chain.hpp.
- TODO `warm-start-salvage` (design on branch wt/warm-start-design, which asks for one staging routine
  for both doors): after this slice. It hands step 2's routine the warm start's own install step, which
  gives the donor's trees rows before they move, and retires the plain mean.
- The record slices of TODO `state-frame-prior` follow: with the grid on dbartsData a re-creation reads
  it there and `adoptGrid` goes.

## Revision from critique

The blind critique (scratch/sgp-crit.md, its probes beside it) said revise. Item by item:

1. A copy or reload after a cross-grid live install got the state's grid, the field holding the state as
   given. Taken: the field holds the engine's state after such an install (The rule; step 3; "field after
   crossing").
2. The forced merge computed another value than the rule stated, and the plan misdescribed the warm
   start's. Taken, under dec-B398: the three passes of The rule, built as the across-grids routine of
   step 2 in place of the flat walk the first draft had; the measured values and the three merged-leaf
   tests; the warm start's plain mean named, and left to warm-start-salvage with this routine.
3. Two test files that pin today's rule were not listed, and `adoptGrid`'s default was not named. Taken:
   step 4 lists six; the default is true in C++ (Constraints), so tests/cpp's existing calls stay.
4. The first call after a reload gave another grid than a later one. Taken as the orchestrator settled it:
   a reloaded sampler is built from its own saved state first (The rule; step 3; "reload, first call").
5. "Numeric against factor: error" was false for an ordered factor. Taken: the refusal is exact where the
   sampler holds the ordered factor, and the other direction, which no state can tell, is stated (The
   rule, A column's kind; "kinds").
6. The equivalence compares and snapshots prove nothing of the same-grid path. Taken: Constraints names
   the proof.
7. Unpinned: below an emptied node, the mirror, gp leaves, an ordered factor in the by-hand fixtures.
   Taken: the two-deep and mirror trees, "linear and gp leaves", the fifth by-hand fixture. The refusal
   the critique foresaw below an emptied node belonged to the flat walk, which is gone.
8. Wording. Taken: TRUE holds only where the trees fit as today (Goal; rows 3 and 4); the malformed row
   is judged first and wins (a threshold on neither grid is an error whatever grid is recorded); the
   recipe is the whole-set call with the state in a variable; weights are lost only on a column whose
   grid changed; the Rd names the `cutPoints` attribute; the bridge's comment on views is corrected; the
   engine estimate is redone.

## Calls made in planning

- A copy and a reload keep adopting the grid the state records until the data object records it (settled
  by the orchestrator, 2026-10-11). Example: set 3 points, store, `copy()`: the copy has 3 points. The
  alternatives: the `cutPoints` record on dbartsData now, about 550 more lines and `setData`'s
  carried-record rule designed first; or deriving afresh, a copy on 100 points with the trees moved, not a
  copy. The cost: a state still installs a grid on these paths, and a copy made after a grid change with
  no store since has the older grid, as today.
- A forced install moves a split by value only, with no argument (settled by the orchestrator,
  2026-10-11). Example: a split at 0.333, the 33rd of 100 points, into {0.25, 0.5, 0.75}, goes to 0.25.
  The alternative, `splits = "value"` or `"position"` as at `setCutPoints`: by position a restore onto a
  refreshed grid of the same length changes every threshold and counts as exact, so a rollback that forgot
  the grid would return TRUE.
- A reloaded sampler is built from its own saved state before `setState` installs another (settled by the
  orchestrator, 2026-10-11). The alternative, adopting the handed state's grid whenever the engine was
  created in the call, gives 100 points and TRUE as a first call and 3 points and FALSE after any other.
- The `state` field holds the engine's state after an install onto another grid or weights (settled by the
  orchestrator, from the critique). The alternative, the state as given, makes a copy another sampler. The
  cost is a `storeState` on such an install.
- Other weights on equal points are clean. Derived, not ruled as such: dec-B450 has the restore go
  through, and nothing in a tree changes. The alternative, FALSE, would make the by-hand sequence return
  FALSE on every weighted column and leave the caller no unforced way in.
- A threshold on the grid is clean and one off it is not: read from dec-B305's rule ("a split off the
  grid") and dec-B234's list, both the entries' words. The alternative, any grid difference being not
  clean, would decline the finer-grid case, where no tree changes.
- The forced move is the state's own install followed by `setCutPoints`'s move, not a rewrite of the flat
  trees: the merged leaf then has the rows dec-B398's mean needs, every leaf kind and the variance forest
  are covered by routines that are tested, and the cost is two re-quantizations of the columns that
  differ on a forced install that moves a split, and none elsewhere.
- An ordered factor in the sampler refuses any other recorded grid; no column kind is recorded in a
  state. The alternative, an added record, tells the reverse direction too and is what dec-B283 removed
  for the family.
- The state keeps both grid attributes: a re-creation has no other grid, pass 1 needs the recorded
  points, and the tests and the by-hand sequence read them.

## Review and fix round

The review of the built slice (2026-10-11) found one defect and test gaps. What the fix round changed, and
the calls in it:

- The ordered-factor refusal belongs to an install that does not adopt (settled by the orchestrator,
  2026-10-11). As built it ran when adopting too, and `copy()` and a reload errored after `setData` had
  replaced an ordered factor with one declaring fewer levels: a live engine keeps the level table it was
  created with, its stored state records that grid (5 points for 6 levels), and a new engine reads the
  levels the data declares (3 points for 4). [`Sampler::setState`](../../src/bartcore/sampler.hpp) now
  makes the check only without `adoptGrid`, and a copy and a reload are again made as before the slice,
  draw for draw.
- An adopting install does not take the grid recorded on a factor column, before the slice or after: no
  install writes one ([`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp)). So in that
  sequence the copy and the reload hold 3 points where the original holds 5, their draws are not the
  original's, and their `setState` refuses the state the original stored, which records 5 points on a column
  whose live grid is 3. Before the slice that call returned TRUE with the grid left at 3. Left as the live
  rule gives it, the state being recorded on another grid than the live one. The copy and the reload are
  made only where no rule of the stored state lies above the declared levels; a rule there that sends
  missing values right keeps both children occupied through `setData`, and then neither can be made. With
  an unordered factor in the same sequence they are refused where a rule sends a dropped level right, and
  made where none does, as of a sampler never run. Both refusals are the base build's. Having a
  re-creation take the level table is scoped out of this slice (settled by the orchestrator, 2026-10-11).
- A sampler's own `state` field always installs on that sampler (settled by the orchestrator, 2026-10-11,
  from the second review). Where a re-creation leaves the engine on other points than the adopted state
  records, which only an ordered factor's column can, the field takes the engine's state
  ([`heldAfterAdoption`](../../R/dbarts.R)), as after a live install onto another grid. So
  `copy$setState(copy$state)` and a reload's `setState` of its own field, as a first call too, return what
  the base build returned, and no draw moves. The alternative, the field as given, had a copy refuse its
  own field until its next `storeState()`.
- A sampler whose engine is gone keeps, at its first `setState`, the kept-draw capacity of its own stored
  state, where the base build took the control's: after `setControl` to 3 draws on a sampler stored with
  5, a reload's first `setState` of another state keeps 5 draws, as a later call does. Base kept 3 as a
  first call and 5 after any other. The other way round too: stored with 4 and the control then naming 6,
  base took 6 as a first call and 4 after any other, and the slice keeps 4. Pinned beside "reload, first
  call" for the shorter control.
- A throw during pass 1 puts the store's grid, weights and codes back before it leaves
  ([`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp)); pinned in
  [`testLiveGridInstall`](../../tests/cpp/test_state.cpp) through the test peer, for a throw out of the
  install step. That the re-quantization before it sits inside the same guard is not pinned: only an
  allocation failure reaches it. A throw that follows a chain's install leaves that chain's trees on the
  recorded grid over a store back on its own, and the sampler is not to be read or run until a state is
  installed; only an allocation failure inside the install reaches that too.
- Pins added, each run against the mutant it answers: the weights through a forced move; a move on the
  second column alone; a threshold off the grid in the second of two mean forests alone; "equals the two
  calls" with a test set and on a sparse store; the engine a declined call leaves a sampler whose engine was
  gone; the refusal on a view in the unforced form; and "direction rides the move" now holds a rule the move
  merges, so the stored state differs where the sides are drawn after it.
- The sentence on 0.9-34's restore is corrected: it brought back the grid and the splits, not the leaf
  values.
- From the second review, two more pins, each run against its mutant: "equals the two calls" on a sparse
  store with a sparse test set, for the test rows' sparse codes put back after a move; and a state recording
  as many points as an ordered factor's grid, other ones, refused.
