# Warm starts that take part of a donor

Status: DECIDED 2026-10-11; slice A of the plan LANDED 2026-10-17 (d03760e3; calls dec-A211), the two forms,
the report and the line, and slice B, the monotone collapse, LANDED 2026-10-17 (ee1d7213; calls dec-A214);
slice C, another tree count, is built; slice D is planned, not built. The six choices of section 3 are
settled as recommended there: (a) A,
(b) A, (c) A, (d) the errors stay, (e) A (dec-B448), (f) one line at every warm start. Plan:
[warm-start-salvage.md](../plans/warm-start-salvage.md), which records what the state grid slice changed
under sections 1 and 4. TODO: warm-start-salvage; slice A closed forced-update-returns-null and slice C
warm-start-tree-count.
Rulings: dec-B247, dec-B295, dec-B310, dec-B323 to dec-B326, dec-B343, dec-B398, dec-B418, dec-B440 and
dec-B448 in docs/decisions.md.

A warm start copies an earlier fit's trees into a new sampler as its starting point: `bart(warm.start = )`, or
the sampler's `installTrees(donor, samples = NULL)`. The earlier fit is the donor and the new sampler the
recipient. This note says what a warm start does with a donor that does not fit the recipient as it stands,
and lists the choices still open. What `setState` does with a state that does not fit is
[install-surface.md](install-surface.md)'s, and a warm start is meant to read the same way.

## 1. Today

Run on 2026-10-11 on a plain build of the branch: gaussian unless said, 300 rows, 4 numeric columns, 10 or 50
trees. `installTrees` has one form. It installs and returns NULL invisibly, or stops with an error. Every
error below is raised before anything is written: the recipient's stored state was identical before and after
in each. An install drops the recipient's kept draws and stores its state again.

| donor against recipient | today |
|---|---|
| same rows, same shape | taken; the starting fit is the donor's |
| recipient on 40 of the donor's rows, same cut grid | taken; a split with a side no row reaches is merged, no word (50 nodes to 48). `setState` of the same trees returns FALSE |
| recipient on other rows, its own grid | taken; every split is moved to the nearest cut of the recipient's grid (0 of 58 donor splits sat on the grid of a new sample of 300; the starting fit differs from the donor's function there by rmse 0.43 against a fit sd of 2.56) and merged where a side empties (166 nodes to 162 on 30 rows) |
| recipient's column holds one value, so no cut | taken; the donor's 4 splits on it merged away (50 nodes to 42). Not refused |
| another number of predictors | error: "malformed cut points in warm-start donor" |
| a column numeric on one side, a factor on the other | error: "this sampler does not support the warm-start donor's predictor structure (a categorical/continuous column mismatch, or a malformed donor cut grid); its trees cannot be remapped onto this cut grid" |
| columns, or a factor's levels, in another order | taken by position (levels a, b, c against c, b, a: the starting fit correlates -0.9995 with the donor's) |
| another tree count, either way | error: "warm-start donor is not shape-compatible with this sampler (number of trees, forests, or predictors differ, or only one of the two carries a variance forest)" |
| DART on one side only | error: "DART state transfers only between two DART fits; the donor and destination disagree on dart". Both DART: taken |
| a variance forest on one side only | error, the shape message above. Both, with another variance tree count: "warm-start donor's variance forest has a different number of trees than this sampler's" |
| constant donor into linear or gp; linear into gp or gp into linear; linear with another number of covariates; gp on other rows | error: "malformed leaf parameters in bartcore state" |
| linear donor into constant | taken, slopes dropped, no word |
| gp donor into constant | taken, no word; the starting fit is not the donor's (range -0.7 to 7.3 against -2.2 to 8.5) |
| linear on x2 into linear on x3 | taken, no word; the x2 coefficients are read as x3's |
| a split on a column the recipient's forest may not use (`blocks`, a restricted variance forest) | error: "warm-start donor holds a tree that splits on a variable outside this forest's allowed column set; the donor's fit is incompatible with the column restriction (a forest's own column subset or a restricted variance forest) in force here" |
| a tree past the recipient's `interactions` limit | error: "warm-start donor holds a tree that violates this sampler's interaction constraint; the donor's fit is incompatible with the interactions() prior in force here" |
| unconstrained donor into a monotone recipient | taken; a tree out of order has every leaf set to 0 (3 of 30 leaves), no word |
| missing values in a column on one side only | taken; sides dropped, or drawn (dec-B440) |
| another response family, at the sampler (7 families, 49 pairs) | taken, fits finite after 3 sweeps, every pair |
| `bart(family = "nbinom", warm.start = )`, and ordinal, multinomial, hurdle | error: 'family = "nbinom" does not support 'warm.start' or 'n.grow.sweeps''; such a fit is not accepted as a donor either |
| recipient of two or more forests | error, by name |

The refusals sit in [`Sampler::installForests`](../../src/bartcore/sampler.hpp) (shape, grid, DART, variance
forest, and the two containment checks, removed in the build,
retired: [`Chain::interactionStateFeasible`](../../src/bartcore/chain.hpp) and
retired: [`Chain::columnMaskStateFeasible`](../../src/bartcore/chain.hpp)), in the bridge's
[`readWarmStartState`](../../src/R_interface_bartcore.cpp), which reads the donor's leaf blocks in the
recipient's shape, and in R ([`refuseMultiForestWarmStart`](../../R/bartcore.R),
[`checkFamilyUnsupportedArgs`](../../R/bart.R)). `bart` prints nothing about a warm start.

## 2. Ruled, and what follows without a further choice

- dec-B295, the maintainer: "Warm starts should try to get as much information out of the fit it can. If
  trees can be collapsed to be made into valid ones, that would make sense. Install trees might have modes
  that work a force flag, where if it was absent it would accept/reject." So `bart(warm.start = )` is best
  effort and always forces; `installTrees` has two forms. A tree the recipient cannot hold is cut back at the
  split that offends: that split and everything beneath it become one leaf, the splits above it stay.
- dec-B310: "Make them NULL. The return value when not forced is an indicator of whether or not a
  statistically valid update was made. That is irrelevant when forcing, because that kind of validity has no
  value." So the forced `installTrees` returns NULL invisibly, and the unforced one TRUE where it installed
  and FALSE where it did not, visibly, as `setState` does. FALSE leaves the sampler as it was: trees, drawn
  values, generators, kept draws, the `state` field. This is forced-update-returns-null.
- dec-B398, "Repair, as we do elsewhere.": a forced `setState` already cuts a tree back at a split on a
  column its forest may not use or past its interaction limit, and sets every leaf of a monotone tree out of
  order to 0. The routine is [`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp) with
  [`Tree::splitIsForbidden`](../../src/bartcore/tree.hpp) and
  [`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp). A forced warm start uses it, so
  the two errors on a column restriction and an interaction limit become a cut-back under force and FALSE
  without. A donor out of the monotone order, and one with a side no row reaches, are FALSE unforced and
  repaired forced, where today they are repaired with no word (dec-B278 for the first). A forced warm start
  collapses a monotone tree out of order where its order breaks, each merged leaf at the mean by rows
  ([`collapseMonotoneTree`](../../src/bartcore/model.hpp)), and does not reset it (dec-B355); the other
  forced calls keep the reset.
- dec-B323, "Sure, option B.": into fewer trees a warm start takes the donor trees whose fitted values vary
  most over the recipient's rows, in the donor's order, values as they are; into more it takes them all and
  starts the rest as single leaves at 0. Settled here, the ruling not naming them: "vary most" is the
  variance of a tree's fitted values over the recipient's rows under its working weights, of the tree as it
  will be installed (moved to the recipient's grid, cut back); each chain chooses among the trees of the
  donor sample it takes; ties go to the earlier tree. Section 3 (e) shows what this does to the sum.
- dec-B324, "Take the trees; a part only one side has does not transfer.": DART weights go only between two
  DART fits, a variance forest only where both have one; the recipient's own part starts as a new fit's.
- dec-B325, "Take what translates.": a constant donor seeds a linear recipient exactly, slopes 0. A linear
  donor seeds a constant recipient with each leaf's value at its covariates' center, which is the stored
  intercept, so today's install is that one; it is FALSE unforced. A gp leaf on either side against another
  kind is an error naming both kinds. Settled here: a linear donor with another number of covariates is
  reduced to its values at the center and then seeds slopes of 0, the two halves of the ruling in turn, FALSE
  unforced. The bridge reads the donor's kind from the donor's own blocks.
- dec-B326, "One line under verbose = TRUE, alongside the setup bart already prints." Its wording is section
  3 (f). Nothing is printed without `verbose`, and no warning is raised.
- dec-B343, dec-B418, dec-B440 stand as built: sigma and k only where the recipient draws them, nothing
  converted between response scales, and the side for missing values drawn on a column the donor could not
  have missing, before the trees are judged. An install that returns FALSE puts the generators back.
- Another response family is not a refusal and does not become one: under dec-B418 the trees are numbers on
  the internal scale, and every pair installs today, judged like any other donor.
- An error in both forms, as for `setState`: not a donor, an older encoding, a malformed block or tree, a bad
  `samples`, a recipient of two or more forests.

## 3. Open choices

Each stands alone. The running example: a donor of 50 trees fitted to 400 rows and 5 columns, and a recipient
on the same rows that allows no interactions (`interactions(max.order = 1)`), so that 21 of the donor's 50
trees hold a split it may not have. One of them, with leaf values on the internal scale:

    donor tree 2                              cut back
    x4 <= 0.228                               x4 <= 0.228
      yes: leaf -0.0075  (84 rows)              yes: leaf -0.0075
      no:  x3 <= 0.288                          no:  leaf  0.0264
             yes: leaf 0.0564  (92 rows)
             no:  leaf 0.0141  (224 rows)

### (a) The argument's name and default

Background. Three sampler methods change a live sampler and may have to alter its trees to do it.
`setPredictor(x, column, forceUpdate, ...)` writes no default: it forces a whole matrix and asks first for one
column. `setState(newState, forceUpdate = FALSE)` asks first. `setCutPoints` always forces and has no
argument. `installTrees` today always installs.

- A. `installTrees(donor, samples = NULL, forceUpdate = FALSE)`.
  `s$installTrees(donor)` prints `[1] FALSE` on the example and `s` is as it was;
  `s$installTrees(donor, forceUpdate = TRUE)` prints nothing and installs the 50 trees, 21 cut back;
  `if (!s$installTrees(donor)) s$installTrees(donor, forceUpdate = TRUE)` is the idiom.
- B. The same argument with `forceUpdate = TRUE` as the default: `s$installTrees(donor)` installs, as today,
  and returns NULL; a caller who wants the answer writes `forceUpdate = FALSE`.
- C. A as `force = FALSE`, base R's spelling (`unlink`, `detach`).

Cost. A: about 10 lines of R. A call that today installs a donor needing repair returns FALSE and installs
nothing until it says `forceUpdate = TRUE`: 69 calls in 23 test files to read, none in stan4bart, bartCause,
treatSens or bairrtt, and no release has had `installTrees`. A caller who drops the value runs on from a cold
start with no message, as with `setState`. B: no call changes, and the class's two installs default opposite
ways. C: a second name for one idea in one class. A default cannot change after release.

Recommended: A. It is what "if it was absent it would accept/reject" says, and what the maintainer ruled for
`setState` ("Yes, unforced should be the default.").

### (b) What the unforced form declines

Background. Unforced, `installTrees` returns TRUE or FALSE, and something must define TRUE. Proposed: TRUE
means the sampler now starts at the donor's function, every tree it took as the donor holds it and no donor
tree left out. That gives FALSE for a tree that must be cut back, merged or have a side dropped, for a
recipient of fewer trees, and for linear leaves reduced to constants; TRUE for a recipient of more trees (the
added leaves are 0) and for a constant donor into linear leaves. DART weights or a variance forest that only
one side has are not judged, as `setState` does not judge DART weights. One case needs a ruling: a donor
fitted to other rows. Its splits sit between its own data's values, and the recipient can only split at its
own cut points, so each split is moved to the nearest one.

- A. A moved split is a changed tree: FALSE. On a new sample of 300 rows from the donor's distribution, where
  none of 58 splits sat on the recipient's grid and none was lost, `s$installTrees(donor)` is FALSE.
- B. A moved split is not a change; only a lost split, tree or slope is. The same call is TRUE; on a sample
  of 30 rows, where two splits merged, it is FALSE.

Cost. A: nothing beyond the check, about one sweep as `setState`'s is. The unforced form then says TRUE only
for a donor on the recipient's own grid, and a caller on new rows must force and is not told whether anything
was lost. B: about 20 lines; TRUE then covers a starting fit that is not the donor's (rmse 0.43 against a fit sd
of 2.56 in the run), and `installTrees` and `setState` mean different things by FALSE once `setState` moves
splits between grids (state-grid-replaces-sampler), dec-B305's rule counting "a split off the grid" as not
clean. Either can be changed later only by changing what a released call returns.

Recommended: A, for one meaning of FALSE across the two installs.

### (c) The value of a leaf made by cutting a tree back

Background. In the example the right-hand leaf replaces two. Its value matters until the recipient's first
sweep reaches that tree and draws its leaves again.

- A. The mean of the leaves it replaces, weighted by the recipient's rows in each:
  (92 x 0.0564 + 224 x 0.0141) / 316 = 0.0264. The fit's mean over the rows is unchanged.
- B. Zero. The tree's contribution on those 316 rows is dropped.
- C. A fresh draw from the leaf prior.

Run: the example, 20 recipients an arm, rmse of the current fit to the truth (the donor's 0.84; sd of the
response 4.9). "cold" is no warm start.

| sweep | A, weighted mean | B, zero | C, draw | cold |
|---|---|---|---|---|
| 0 | 2.33 | 3.48 | 4.26 | 4.83 |
| 1 | 1.57 | 1.62 | 1.71 | 1.99 |
| 5 | 1.44 | 1.45 | 1.48 | 1.51 |
| 20 | 1.36 | 1.37 | 1.37 | 1.36 |

At sweep 0 A is 2.21 from the donor's fit with no shift of the mean, and B 3.34 with the mean 1.99 lower.

Cost. A: no new code; a forced `setState` does exactly this (its result on the example equals A to 4e-15),
geometric for a variance forest. B: about 15 lines and a second rule in one routine. C: about 40 lines, and a
forced warm start then moves the generator by the number of leaves merged.

Recommended: A. dec-B398 already answers it for `setState`; this confirms it for the warm start.

### (d) What stays an error under force

Background. Force repairs trees. It cannot repair a donor whose trees are not about the recipient's columns.
Proposed errors in both forms, beyond those of section 2:

1. Another number of predictors. Columns are matched by position, so there is nothing to take. New text:
   "warm-start donor was fit to 4 predictors and this sampler holds 3", for today's "malformed cut points".
2. A column numeric on one side and a factor on the other, naming the column.
3. A gp leaf against another kind (dec-B325): "warm-start donor has gp leaves and this sampler constant
   ones", and a gp donor on another number of rows, by name where it is "malformed leaf parameters" today.

Options. For 2: A, the error; B, under force every split on that column is cut back, the merged leaf at the
plain mean of the leaves beneath it, since the split cannot be evaluated on the recipient's rows to weigh
them. For 3: A, the error as ruled; B, under force a gp donor gives its splits with every leaf at 0, and a
constant or linear donor gives a gp recipient its splits with the fits started at 0, after the maintainer's
later "The goal of a warm start is to take tree structure (for the most part)" (dec-B418).

Cost. The errors: about 30 lines, mostly messages. 2 B: about 150 lines with tests and a second leaf-value
rule; such a donor almost always has another number of columns as well. 3 B: about 120 lines, and it sets
dec-B325's refusal aside. An error can be lifted after release without breaking a caller; a donor once taken
cannot be refused again without breaking one. Recommended: A for both.

### (e) A donor of another tree count: the variance forest

Background. dec-B323 rules the mean forest. Run on the example's data, the truth's rmse as in (c):

| donor into recipient | rule | starting fit against the donor's | rmse to truth at sweep 0, 1, 5, 20 |
|---|---|---|---|
| 200 into 50 | the 50 most varying (ruled) | rmse 1.43, sd 3.99 against 4.77 | 1.55, 1.11, 0.98, 0.93 |
| 200 into 50 | the first 50 | rmse 4.15, sd 1.36 | 4.19, 1.93, 1.36, 1.03 |
| 200 into 50 | cold | | 4.83, 2.12, 1.40, 1.11 |
| 50 into 200 | all, then leaves at 0 (ruled) | equal | 0.84, 0.86, 0.87, 0.88 |
| 50 into 200 | cold | | 4.83, 1.64, 1.09, 0.93 |

The 50 most varying of 200 trees carry 78 percent of the summed per-tree variance and the first 50 carry 24.
The sum of 50 trees is not the donor's prediction; it is 1.43 from it, and the sampler closes the rest. A
heteroscedastic fit has a second forest, for the variance, with its own tree count, which the ruling does not
name; today another count there is an error.

- A. The same rule: fewer trees take those whose log scale varies most over the recipient's rows; more
  trees take them all and add leaves at 1, which leaves the donor's variance function as it is.
- B. A donor's variance forest of another count is left behind and the recipient's starts as a new fit's.
- C. The error stays.

Cost. A: about 40 lines beside the mean forest's. B: about 10. C: none, and `bart(variance = , n.trees = ,
warm.start = )` keeps one refusal that the mean forest lost. Recommended: A.

### (f) What `bart` prints

Background. With `verbose = TRUE`, `bart` prints its setup (trees, chains, prior, data, cut points) and then
"Running mcmc loop:". dec-B326 puts one line there, its example "warm start: 20 of 50 donor trees taken; DART
weights started fresh". Open: the words for a cut-back, and whether the line appears when all was taken.

- A. A line at every warm start, trees counted per chain and repairs summed over chains:
  `warm start: 50 of 50 donor trees taken; 21 cut back to fit` on the example;
  `warm start: 75 of 75 donor trees taken` when all went in as they were;
  `warm start: 50 of 200 donor trees taken; 9 cut back to fit, over 4 chains; variance forest started fresh`.
- B. The line only when something was left, changed or started fresh; nothing for a whole donor.
- C. As A, counting splits: `50 of 50 donor trees taken; 27 of 74 splits removed`.

Further phrases under any option: "150 trees started as single leaves", "slopes
left behind", "slopes started at 0", "DART weights started fresh", "donor's DART weights left behind",
"donor's variance forest left behind".

Cost. About 40 lines of R and a count returned by the bridge. The words can change after release.

Recommended: A. A line that is always there says a warm start happened; its absence means none.

## 4. The build

One entry point carries both forms. [`Sampler::installForests`](../../src/bartcore/sampler.hpp) takes a force
flag and fills a report: per chain, trees taken of trees offered, trees cut back, parts left
behind or started fresh. It assembles each chain's trees as now, adapts the count and the leaf kind, draws the
missing-value sides, and then judges every tree on scratch as it will be installed, the staging
[`Chain::stateIsValid`](../../src/bartcore/chain.hpp) runs for a state
([`Tree::stageFromFlat`](../../src/bartcore/tree.hpp)) extended to a donor on another grid: built against
the donor's grid under [`ScopedCutGrid`](../../src/bartcore/data.hpp), moved
([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)), partitioned and marked. Unforced, anything
marked returns a new [`WarmStartResult`](../../src/bartcore/sampler.hpp) value, declined, with nothing
written. Forced, [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) and the remap arm, removed in
the build, retired: [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp), install, the second arming the forbidden-
split predicate as the first does; the two containment checks stop being refusals and feed the verdict. The
bridge's [`bartcore_installForests`](../../src/R_interface_bartcore.cpp) takes the flag and returns the
report; [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads the donor's leaf blocks in the
donor's shape. In R one internal function returns the report; [`installTrees`](../../R/dbarts.R) turns it
into TRUE, FALSE or NULL, and `bart` into its line.

Two slices, each on the engine tier (a plan with a blind critique, one review that runs mutants), serial:

1. The two forms: the verdict on both grid paths, force, the cut-back at a forbidden split on the path for
   another grid, the monotone reset and the empty side into the verdict, the return values, the report, the
   argument, `bart` forcing, the line for a cut-back, the help, and the 69 calls. About 950 lines planned,
   500 of them tests; at the 1.5 to 2 times such slices have run, 1400 to 1900 landed. A donor that installs
   whole today installs bit for bit; no stored draw moves.
2. What does not match in kind: the tree count for both forests, the parts one side has, the leaf kinds and
   the errors of (d), with their phrases in the line. About 1150 lines planned, 600 of them tests; 1700 to
   2300 landed. It needs the first slice's FALSE and its report.

Shared with other queued engine work: [`Sampler::installForests`](../../src/bartcore/sampler.hpp) and the
virtual of that name on the facade, the rebuild and validity routines in chain.hpp, the staging, moving and
collapse routines in tree.hpp, the bridge's state parsers, the `installTrees` method and the two help pages.
state-grid-replaces-sampler needs the same staging of a tree against a grid that is not the sampler's, with
the same verdict; whichever lands second reuses the first's routine and should not write another. The forest
slices touch the same files in other functions (the combiner, creation, the leaf prior); a recipient of two or
more forests stays refused, so nothing of theirs is installed here. Engine slices land one at a time.

## Not in this note

- Columns and a factor's levels are matched by position, with no check and no word (section 1). A donor given
  as a sampler or a fit carries its names and levels, a stored state does not.
- A linear donor on other covariates, in the same number, is read as on the recipient's; a state does not
  record which columns a linear leaf uses.
- `bart` refuses `warm.start` for ordinal, negative binomial, multinomial and hurdle fits, and refuses such a
  fit as a donor, where the sampler takes every pair of families; and a recipient of two or more forests.
