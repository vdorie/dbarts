# Installing a stored state: the checked form and the forced form

Status: LANDED 2026-10-10: Part A (missingness first seen, 2269918e) and Part B (setState's two forms, 062ef59a).
Plan: [setstate-force-update.md](../plans/setstate-force-update.md). Rulings: dec-B305, dec-B310, dec-B318,
dec-B398 to dec-B402 and dec-B418 in docs/decisions.md. The cut grid a state goes onto
([state-grid-install.md](../plans/state-grid-install.md), landed 2026-10-11): dec-B449, dec-B450 and dec-A209.

A state is the chain, and installing one puts a chain into a sampler that has its own model and its own data.
Most states fit. Some hold a tree the sampler cannot hold as it is stored. This note says what `setState` does
with each, and is the home of those rules; what a state carries is
[state-not-model.md](state-not-model.md)'s, and which columns can hold a missing value is
[mia-missingness.md](mia-missingness.md)'s.

## The two forms

`setState(newState, forceUpdate = FALSE)` asks first. A state every tree of which fits is installed and the
call returns TRUE. A state with a tree that would have to be changed is not installed: the call returns FALSE
and the sampler is as it was, its trees, drawn values, latents, generators, cut grid, kept draws and cached
`state` field. A sampler whose engine was gone is first rebuilt from its own saved state
([The cut grid is the sampler's](#the-cut-grid-is-the-samplers)), and is then as that rebuild leaves it; one
with no saved state that installs keeps its dead pointer. Both values are returned visibly, as an unforced
`setPredictor`'s are.

`setState(newState, forceUpdate = TRUE)` installs, repairing what does not fit, and returns NULL invisibly
whether or not anything was repaired (dec-B310: the unforced value says whether a statistically valid update
was made, which has no meaning once the update is forced). `copy()` and a reload force and report nothing,
having no other state to fall back on; so does stan4bart's reload (dec-B402).

A state that is broken, or of a shape that does not fit, is an error in both forms: not a state, an older
encoding, another count of chains, forests or trees, another leaf model, a malformed block or tree, a cut grid
that repeats a point, a live threshold that is no point of the grid the state itself records, another grid
recorded on a column the sampler holds as an ordered factor, latents the family cannot hold, saved gp draws
under other lengthscales, saved-tree blocks that name different numbers of draws.

## What does not fit, and its repair

Only the live trees are judged, mean and variance, each built from its stored form and partitioned over the
sampler's rows. A kept draw is replayed over new rows and never held to this partition.

| a tree that | forced |
|---|---|
| holds a split whose threshold is not a point of the sampler's cut grid | the split moves to the nearest point inside the range the splits above it leave, or, with none left, becomes one leaf with everything beneath |
| has a leaf no row reaches | the split above it becomes one leaf with everything beneath |
| holds a split outside the range the splits above it leave | the same |
| records a side for missing values on a column, not pooled, that has never held one here | the side is dropped |
| splits on a column its forest may not use (`vars`, a `blocks` group, a moderator's columns, a restricted variance forest) | cut back at the first such split from the root |
| breaks its forest's `interactions` limit | cut back at the first such split from the root |
| is monotone with leaf values out of order as stored | every leaf set to 0 |

A merged leaf takes the mean of the leaves it replaces, weighted by their rows under the state's own working
weights, geometric in the variance forest. On each path from the root the first offending split goes and the
splits above it stay. A merge can leave a monotone tree out of order, and it is then reseeded as well. The
last three rows are repaired under dec-B398 ("Repair, as we do elsewhere."). A warm start judges and
repairs a donor's trees by the same table and the same staging
([`Chain::donorTreesStage`](../../src/bartcore/chain.hpp)): `installTrees(donor)` returns FALSE where a
row applies and `forceUpdate = TRUE` makes the change ([Warm starts](#warm-starts)).

The verdict and the install stage a tree by one routine,
[`Tree::stageFromFlat`](../../src/bartcore/tree.hpp), the verdict on a scratch tree
([`Chain::stateIsValid`](../../src/bartcore/chain.hpp)) and the install on the live one
([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp),
[`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp)), so the unforced call declines exactly the
states a forced install reports as altered. The repair is one pass of
[`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp), whose walk also takes a forbidden split
([`Tree::splitIsForbidden`](../../src/bartcore/tree.hpp)).

The verdict builds and partitions every tree once more than a forced install does, so an unforced restore
costs more. Measured with 200 trees on 10 columns, one chain and one thread (bench-sampler.R's `restore`
grid, 2026-10-10, on a machine under load, two runs within one percent), an unforced restore took 1.64 times
a forced one at n = 1000 and 1.72 times at n = 1e5, which is 0.9 and 1.3 sweeps more; a forced restore cost
what a restore cost before the check existed (1.000 and 1.001 times it). A loop that restores its own state
onto rows it has put back, where nothing can need repair, passes `forceUpdate = TRUE` and pays that.

## What is clean whatever it differs in

- The response mapping. Nothing is converted on install (dec-B418), so a state stored before a re-anchor or by
  a sampler on another response goes in as stored and is clean.
- The latents. Under other case weights or another censoring status they are redrawn after the install, from
  the installed generators.
- `k`. It is installed as recorded whatever leaf prior the recipient holds (dec-B401); the leaf sd is then
  `k.scale / k` against the recipient's `k.scale`. A recipient holding `k` fixed keeps its own.
- The kept draws (dec-B318: they are output, not model or chain). A state from a store the size of the
  sampler's has its ring copied slot for slot, so a copy continues its source's layout. From a store of any
  other size the newest recorded draws that fit go to the sampler's first slots, oldest first; `n.samples` and
  the store's size do not change, and a sampler keeping none takes none
  ([`Sampler::setState`](../../src/bartcore/sampler.hpp)). A re-created sampler still takes the state's
  capacity, having no other record of it.
- A value the sampler holds fixed, DART weights, and a generator of another kind, each left as the sampler's.

## A state from before a column could hold a missing value

A stored state records which columns could hold a missing value when it was stored, the store's flags and not
its content, as the attribute `missing.columns`
([`SamplerStateData`](../../src/bartcore/sampler.hpp)); a state without the attribute says nothing. Installed
into a sampler where a further column can, the state holds no direction on that column, and every rule on it,
in the live trees and in each recorded kept draw of the state's own ring, draws one with probability one half,
as the column's first missing value would have (dec-B435, extending dec-B400). The draw is made on a copy of
the state's chains before anything of the sampler is written
([`Chain::completeMissingDirections`](../../src/bartcore/chain.hpp)), so the verdict judges the completed
trees: such an install is clean, and TRUE, unless a drawn direction leaves a leaf no row reaches or a monotone
tree out of order. The coins are the state's own generators', which the install then leaves past them: a
state stored just before one whole-column or whole-matrix update raised its columns, and installed just after
it, holds the update's directions and its generators. Nothing else is claimed to reproduce. A declined or
refused state leaves each chain's generator where it was. The `state` field is the state as given, so a copy
or a reload of it draws the same directions; after an install onto another grid or weights the field is the
sampler's state past the draw, and a copy or a reload draws nothing.

A kept rule on a factor column that sends no level right and the missing value right is left as it is by this
draw and by the first-sight draw ([`drawFlatMissingDirections`](../../src/bartcore/tree.hpp)): sent left it
would send nothing right. Only kept draws taken from another sampler's state hold such a rule on a column
never flagged. A column the state records and the sampler has never seen missing is as above: its direction
is dropped, and a state never raises a flag (dec-B321).

A warm start draws the same for the trees it installs, from each receiving chain's own generator and before
the trees are judged, moved or merged ([`Sampler::installForests`](../../src/bartcore/sampler.hpp),
dec-B440); a donor state with no record draws nothing, and a start refused or declined puts the generators
back.

## The cut grid is the sampler's

A state's splits are stored by value, and a tree is built by putting each on the grid position that holds its
value. Until this rule a live `setState` therefore installed the grid the state records first, with its
weights, and judged the trees there, which is what released 0.9-34 did: three points set with `setCutPoints`
and an earlier state restored gave the earlier 100 back. Cut points are data, derived from the predictors and
then held while the predictors move (dec-B231), and a state is the chain (dec-B195), so a restore that
rewrites them changes the model the chain is read against without being asked. dec-B449 rules that a live
`setState` never writes the sampler's cut points or their weights. The grid and its weights are part of the
tree prior, so this changes the posterior a chain targets after one call sequence: a state installed into a
live sampler whose grid, or its weights, differ from what the state records. An install onto the recorded
grid and weights, every copy and reload, and every fit that never calls `setState` are bit for bit what they
were.

On a sampler whose engine existed before the call, the first row that applies:

| The state against the sampler | Unforced | Forced |
|---|---|---|
| a live threshold that is no point of the grid the state records; a recorded grid that repeats a point; a column an unordered factor on one side only; another recorded grid on a column the sampler holds as an ordered factor | error | error |
| recorded points and weights equal the sampler's | as before | as before |
| the points equal, the weights not | TRUE where the trees fit; the sampler's weights untouched | NULL; the same |
| other points, every live threshold a point of the sampler's grid | TRUE where the trees fit, FALSE for a leaf no row reaches and the like; grid and weights untouched | NULL; the repairs of the table above |
| a live threshold that is no point of the sampler's grid | FALSE, nothing touched | NULL; the move below |

Whether the thresholds are points of a grid is read from the flat trees with no tree built
([`Chain::liveThresholdsAreOnGrid`](../../src/bartcore/chain.hpp)), of the recorded grid first, where a miss
is the refusal, and then of the sampler's. A state with a threshold off the sampler's grid is judged on the
recorded one, swapped in for the builds ([`ScopedCutGrid`](../../src/bartcore/data.hpp)).

The forced move is three passes over each live tree, mean and variance
([`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp)):

1. The tree is built as stored and the sampler's rows are routed through it by its stored thresholds, the
   columns whose grids differ re-quantized on the recorded grid for the purpose. Each stored leaf then holds
   the rows the stored tree sends it, the receiving sampler's rows where the state is another sampler's. A
   leaf no row reaches is merged here, as on any forced install.
2. The sampler's grid, weights and codes are put back as they were, byte for byte, and each split on a column
   whose grid differs goes to the nearest point of the sampler's grid inside the interval the splits above it
   leave, a tie to the upper point ([`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp)). A split with no
   point left becomes, with everything beneath it, one leaf at the mean of the leaves beneath weighted by the
   rows each held in pass 1, the plain mean where none holds a row (dec-B398).
3. The rows are routed again by the moved thresholds and a leaf left without a row is merged, as after a
   forced `setPredictor`.

That is the state installed on its own grid followed by `setCutPoints(grid, splits = "value")`, the two calls
that existed, with the sampler's weights kept; passes 2 and 3 are that call's own routine
([`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp)), so linear and gp leaves, the variance forest
and the monotone reseed are handled as it handles them, and
["equals the two calls"](../../inst/tinytest/test-state-grid-install.R) holds the two to one result. Measured
on one tree of 200 rows, root on the 50th of 100 points, left child on the 30th, leaves 0.10, -0.20 and 0.05
holding 73, 38 and 89 rows, forced onto the 50th and 75th points: the child has no point left and its two
leaves become one at -0.0027. The plain mean is -0.05, and weighing by the rows after the child moved gives
0.10. A fitted state of 20 trees forced from 100 points onto {0.25, 0.5, 0.75} keeps the three points; its
17 splits on the column become 15 and 29 splits in all become 27.

- A move is by value only and takes no argument. By position a restore onto a refreshed grid of the same
  length would change every threshold and count as exact.
- The weights are part of the grid. Other weights on equal points place no threshold and change no tree, so
  such a state is clean (dec-B450): TRUE says the trees are the stored ones, not that the chain continues as
  the stored chain would, since it runs on under the sampler's weights.
- A state records no column kind (dec-B283). An unordered factor's entry is empty, which tells it from any
  other column. An ordered factor's grid is fixed by its level count, so on a column the sampler holds as one
  any other recorded grid is another column's and is refused. The reverse cannot be told: an ordered
  factor's entry is what a numeric column of the same whole numbers records, so such a state into a numeric
  column follows the numeric rule.
- A rule's side for missing values rides the rule through a move, and the sides a state was stored without
  are drawn before anything moves.
- Kept draws go in as stored and are never moved, so they predict what they predicted. A kept draw can hold
  a threshold the grid no longer has, as after `setCutPoints`.
- Several chains share one grid: one chain's tree that needs a move declines the whole state.
- The `state` field. After an install onto a grid or weights other than the state records, forced or not, the
  field holds what `storeState()` gives, so a copy and a reload reach the sampler's grid, weights and trees.
  After an install onto the recorded grid and weights it holds the state as given.

A copy and a reload create a new engine, whose grid is the one a creation derives from the predictors as they
now are. The stored state is the only record of a grid set with `setCutPoints` or held while a column moved,
so a re-creation still takes the grid and weights its state records and then installs the state; taking that
away before the data object records the grid (TODO state-frame-prior) would make a copy of a sampler with a
set grid something other than a copy. A factor column is outside this: no install writes its grid
([`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp)), which stays the one its level table
gives, so a re-creation does not hold its state to the grid recorded on an ordered factor, as a live install
does. A live engine keeps its level table through `setData`, and a new one reads the levels the data record
declares, which `setData` keeps the sampler's: it re-codes each factor row by label onto the sampler's table
([`alignDataFactorLevels`](../../R/bartcore.R), dec-B454), so a copy and a reload hold the original's grid
and take the state it stored. A record declaring other levels than its engine holds is reached by
assigning one by hand, by an object saved by an earlier build, and by a sampler made from a matrix of
codes, which records no table and takes the table of a data frame handed to `setData` (TODO
code-matrix-category-count); a copy and a reload of it hold the grid of the declared levels, are made only
where no rule of the stored state's trees lies above them or sends a dropped level right, and refuse a
state the original stored. Their `state` field is then the state their
engine holds, read once it is made ([`heldAfterAdoption`](../../R/dbarts.R)), so that a sampler's own field
is never refused for the cut points it records. `setState(newState)` on a sampler whose engine is gone first rebuilds
the engine from the sampler's own saved state, where it has one that installs, and then applies the rule
above, so the call gives one grid whether or not the engine was live. Only a sampler with no saved state, or
one that is refused, takes `newState`'s grid.

Putting a grid back is the caller's, before `setState`, with the state held in a variable:
`setCutPoints(attr(st, "cutPoints"))`, the whole set in one call, and then `setState(st)`. Where no changed
column's stored grid carries weights that is the stored grid and the stored chain, bit for bit, which is what
`setState` alone gave before this rule; 0.9-34's own restore brought back the grid and the splits and not
the leaf values. A column whose grid was changed and whose stored grid carries weights gets its points
and not its weights, `setCutPoints` installing any grid other than the held one unweighted: the restore is
TRUE and the chain runs with that column's points equally likely.
`setPredictor(x[, j], j, updateCutPoints = "value")` gives back both while the column holds the values the
grid was derived from, and nothing does once it holds others (TODO cut-point-weights).

## Warm starts

`installTrees(donor, samples = NULL, forceUpdate = FALSE)` has `setState`'s two forms and values
(dec-B295, dec-B310; [warm-start-salvage.md](warm-start-salvage.md)). Unforced it returns TRUE for a donor
installed as it holds its trees and FALSE, the sampler untouched, for one that would have to be changed;
forced it installs with the changes and returns NULL invisibly. A donor is changed where a row of the
table above applies to a tree a chain takes, where a live threshold is no point of the sampler's cut grid,
and where its leaves would not be the donor's: gp fits recorded on other cut points, which start at 0;
slopes or gp fits that a sampler of constant leaves does not read; and linear coefficients kept across a
covariate that holds one value on each side and not the same one. The sides drawn under dec-B440 and
coefficients restated in the sampler's standardization are how a donor is read and are not changes.

The forced move across grids is the state install's routine handed the warm start's own install step
([`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp)), so a merged leaf weighs by rows there
too and the sampler's cut points, weights and codes are the same before and after. What differs from a
state install: the donor's generators, latents and kept draws are not taken, the sampler's kept draws are
dropped, and its `state` field is stored again after every install. On a column the sampler holds as an
ordered factor a donor recorded on other points is an error in both forms, as for a state. A monotone tree
out of order is reset to 0, a tree count, DART setting or variance forest on one side only is an error, and
a leaf kind is read in the sampler's shape (TODO warm-start-salvage, warm-start-tree-count).
`bart(warm.start = )` forces and, under `verbose = TRUE`, prints one line saying what was taken and
changed, from counts read off the installed trees ([`Chain::countDonorTrees`](../../src/bartcore/chain.hpp)).
