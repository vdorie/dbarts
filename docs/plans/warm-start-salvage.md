# warm-start-salvage: a warm start takes what fits of a donor

Status: slices A and B LANDED 2026-10-17 (A a5f3c9bc..d03760e3, calls dec-A211; B f5c16b4a..ee1d7213, calls dec-A214;
[Slice A as built](#slice-a-as-built), [Landing note](#landing-note), [Slice B as built](#slice-b-as-built),
[Landing note, slice B](#landing-note-slice-b)); slice C is built ([Slice C as built](#slice-c-as-built)); D is
planned, not built. Planned 2026-10-11, revised from its blind
critique, no call open. Four slices, A to D, serial. Design: [warm-start-salvage.md](../design/warm-start-salvage.md), whose six choices are settled: (a)
`forceUpdate = FALSE`, (b) FALSE for a split that has to move or merge, (c) the mean weighted by rows, (d) the
three errors stay, (e) the mean forest's rule for the variance forest (dec-B448), (f) one line. Where this plan
and the design note part, on the monotone tree and the number of slices, the plan is the later record
([Revision from critique](#revision-from-critique)).

agent: per slice, one opus implementer and one opus reviewer who runs that slice's mutants (Verification).
Tier "Changes draws" of [README.md](README.md): every slice edits the engine's install path.
rng, per slice:
- A is NEUTRAL, bit for bit, for every fit that makes no warm start, and gives today's result for a forced
  warm start from a donor whose every live threshold is a point of the recipient's cut grid. Where the
  donor's recorded points are the recipient's it runs the code that runs today. Where they differ and every
  threshold is still a point of the grid it runs the same-grid arm where today's build runs the remap arm:
  the two are identical on 8 fixtures (scratch/wsc/p6.R, run on sgi), and a difference on any of them is a
  stop. The one such donor that no longer installs is the ordered factor on another grid, an error with no
  draws. A is SHIFTING for a forced warm start across grids, `bart(warm.start = )` onto other rows included:
  the starting trees are the across-grids routine's and not today's remap's, so the draws that follow differ
  and the posterior does not.
- B is NEUTRAL except for a forced warm start into a monotone recipient that holds a donor tree out of
  order, which is SHIFTING: the tree starts collapsed where today every leaf is 0.
- C and D are NEUTRAL for every donor of the recipient's tree counts, leaf kind, DART setting and variance
  forest. What they newly install is an error today.
- No install draws beyond dec-B440's. No recorded scenario, snapshot file or exact gate calls `installTrees`
  or `warm.start` (searched), so no baseline is re-recorded.
window: after state-grid-install, which has landed (a15b05ca to 480d7d69), and before the merge to main
(TODO `forced-update-returns-null`). Engine slices stay serial: A, B, C, D.
budget: A about 1730 lines, B about 455, C about 820, D about 950; the table and the stops are under Budget
and stops.

## Goal

`installTrees(donor, samples = NULL, forceUpdate = FALSE)` asks first: a donor every tree of which stands in
the recipient as stored is installed and the call returns TRUE; one that would have to be changed is not
installed, the call returns FALSE and the sampler is as it was. With `forceUpdate = TRUE` the donor is
installed with the changes and the call returns NULL invisibly. `bart(warm.start = )` forces, and with
`verbose = TRUE` prints one line saying what was taken and what was changed. After slice B a monotone tree
out of order is collapsed where the order breaks and not set to 0. After slices C and D a donor of another
tree count, of constant against linear leaves, or with DART or a variance forest on one side only is taken
in part. Closes TODO `warm-start-salvage`, `warm-start-tree-count` and `forced-update-returns-null`.

## Context

Rulings built, in the register's words unless quoted as the maintainer's:

- dec-B295, the maintainer: "Warm starts should try to get as much information out of the fit it can. If
  trees can be collapsed to be made into valid ones, that would make sense. Install trees might have modes
  that work a force flag, where if it was absent it would accept/reject." Slices A and B.
- dec-B355: "A warm start is no longer under this call: dec-B295 collapses a donor tree at a split it cannot
  hold." The call it leaves is the monotone reset of dec-B278, which `setData`, `setCutPoints`, a forced
  `setPredictor` and `setState` keep. Slice B.
- dec-B310, the maintainer: "Make them NULL. The return value when not forced is an indicator of whether or
  not a statistically valid update was made." Built for `setPredictor` and `setState`; slice A builds it for
  `installTrees`.
- dec-B398, "Repair, as we do elsewhere.": the cut-back at a forbidden split and the merged leaf at the mean
  weighted by rows. Built for `setState`; slice A reaches it from the warm start, and slice B merges by it.
- dec-B326, "One line under verbose = TRUE, alongside the setup bart already prints." Slice A, with the
  later slices' phrases.
- dec-B323, "Sure, option B.": "a warm start into fewer trees takes the donor's trees whose fitted values
  vary most over the recipient's rows, in the donor's order and with their values as they are; into more
  trees it takes them all and starts the rest as single leaves at 0". dec-B448, "Use the same rule.": the
  variance forest, "selection by variation of the tree's log scale over the recipient's rows and padding
  with leaves at 1". dec-B247 is the request both answer. Slice C.
- dec-B324, "Take the trees; a part only one side has does not transfer.": "a DART recipient of a plain
  donor starting its weights as a new fit does and drawing them from the trees at its first sweep; a
  variance forest transfers only where both have one, a recipient's own starting as a new fit's and a
  donor's left behind". dec-B325, "Take what translates." Slice D.

Rulings that stand as built and bound this plan: dec-B343 (sigma and k only where the recipient draws them),
dec-B418 (nothing converted between response scales), dec-B440 (the side for missing values drawn before the
trees are judged), dec-B328 and dec-B329 (a linear slope converted by the formula), dec-B449 and dec-B450
(the cut grid and its weights are the sampler's; a warm start never wrote them), dec-B305 and dec-B234 (what
clean means), dec-B278 (the monotone reset on a forced update), dec-A199 and dec-A201 (the build calls of
the two state slices).

Out of scope, each with its TODO item: columns and factor levels matched by position
(`warm-start-by-position`); a linear donor on other covariates in the same number
(`warm-start-linear-covariates`); `bart`'s refusal of `warm.start` for nbinom, ordinal, multinomial and
hurdle fits (`warm-start-bart-families`); a donor that is another model (`warm-start-across-models`, after
1.0-0, dec-B248); a test of the undo of a failed install (`warm-start-undo-test`). A recipient of two or
more forests stays refused ([`refuseMultiForestWarmStart`](../../R/bartcore.R)); no item asks for it.

Released 0.9-34 has neither `installTrees` nor `warm.start` (run on scratch/libs/cran: the sampler's methods
matching install, warm or State are `setState` and `storeState`; `bart`, `bart2` and `rbart_vi` have no
`warm.start`). So no call's behavior changes from a release, and NEWS describes a new feature.

Code as it stood before slice A:

- [`Sampler::installForests`](../../src/bartcore/sampler.hpp) assembles each chain's trees from the donor
  sample it takes, restates linear coefficients
  ([`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp)), draws the missing-value sides
  ([`Chain::drawDonorMissingDirections`](../../src/bartcore/chain.hpp)), refuses on four scratch passes,
  which slice A removed (retired: [`Chain::interactionStateFeasible`](../../src/bartcore/chain.hpp);
  retired: [`Chain::columnMaskStateFeasible`](../../src/bartcore/chain.hpp);
  retired: [`Chain::forestsRebuildable`](../../src/bartcore/chain.hpp);
  retired: [`Chain::varianceForestRebuildable`](../../src/bartcore/chain.hpp)) and installs by
  [`Chain::installForest`](../../src/bartcore/chain.hpp) and
  [`Chain::installVarianceForest`](../../src/bartcore/chain.hpp). It reports nothing of what the install
  changed.
- [`Sampler::setState`](../../src/bartcore/sampler.hpp) reports `altered`: set by the installing build where
  a live tree is repaired or moved across grids, never by a scratch build, "so a forced install sets it
  exactly when the same call without force declines".
- A donor on the recipient's cut points goes through
  [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), the routine `setState` installs by: it stages
  each tree ([`Tree::stageFromFlat`](../../src/bartcore/tree.hpp)), which arms the forbidden-split collapse
  of [`Tree::collapseEmptyNodes`](../../src/bartcore/tree.hpp), copies a gp leaf's fits as stored, and sets
  every leaf of a monotone tree out of order to 0
  ([`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp)). Only the refusals ahead of it
  keep a warm start from cutting a tree back today.
- A donor on other cut points goes through the remap arm, which slice A removed,
  retired: [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp): built on the donor's grid with no rows
  under it, moved, then partitioned, so a merge at the move takes the plain mean, and gp fits start at 0.
- [`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp) is the state grid slice's routine: the
  recorded grid goes in with re-quantized codes, the install step it is handed builds and partitions the
  trees there, the store's grid is put back, and
  [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp) moves the splits by value with the rows each
  leaf held, and resets a monotone tree the move leaves out of order. A throw from its install step, or from
  the re-quantization before it, leaves with the store's grid, weights and codes back and nothing moved.
  [`Chain::liveThresholdsAreOnGrid`](../../src/bartcore/chain.hpp) reads whether a
  state's live thresholds are points of a grid with no tree built.
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) holds the verdict's staging in its two tree loops.
  It cannot be called on a warm start's assembled chain, which carries no latents or kept draws.
- The monotone order is [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp): two leaves are related
  along a monotone threshold column when their intervals of cut points touch there and some value reaches
  both on every other column the tree splits on ([`MonotoneLeafGeometry::relation`](../../src/bartcore/model.hpp)),
  and the lower one's value may not exceed the upper one's by more than 1e-9. It reads rules and values,
  never rows.
- The bridge's [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads the donor's leaf blocks in
  the recipient's shape: a recipient of constant leaves reads no `tree.params`, so a linear or gp donor
  installs there. [`bartcore_installForests`](../../src/R_interface_bartcore.cpp) takes three arguments and
  no predictor matrix. [`installTrees`](../../R/dbarts.R) returns NULL invisibly or stops.
- The verbose setup lines are printed at creation by the bridge
  ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)); `bart` installs the donor after creation
  and before its first run.
- DART at creation ([`DartPrior::initialize`](../../src/bartcore/model.hpp)): split weights equal over the
  columns the forest may use and 0 elsewhere, the concentration at the value the options give, the delay
  counter at 0. [`DartPrior::update`](../../src/bartcore/model.hpp) draws nothing until the counter reaches
  `updateDelay`, by default half the burn-in, and the weights stay equal until then. A variance forest at
  creation ([`Chain::buildVarianceForest`](../../src/bartcore/chain.hpp),
  [`VarianceForest::initialize`](../../src/bartcore/chain.hpp)): one root per tree, their product the
  starting variance on every row, sigma fixed at 1.

### Measured on this tree

Run on scratch/libs/sgi (shipped build of state-grid-install at acfd56ac, before its review's fix rounds, with the design note on top), gaussian unless said, 300
rows, 4 numeric columns, 10 to 50 trees, one chain. Every error left the recipient's stored state identical.
A row that names a script is a run of the critique's or of this revision's, on the same build, the script
under scratch/wsc.

| donor against recipient | sgi |
|---|---|
| same rows, same shape | taken, NULL invisibly; starting fit equals the donor's (difference 0) |
| recipient on 40 of the donor's rows, the donor's cut points set | taken; 46 nodes to 42. `setState` of the same trees: FALSE; forced, the same trees exactly |
| recipient on other rows, its own grid | taken; 0 of 18 donor splits on the recipient's grid; 46 nodes to 46; rmse to the donor's function 0.18 against a fit sd of 1.79. `setState`: FALSE, and forced NULL with the recipient's grid kept |
| a column holding one value | taken; the donor's 5 splits on it gone (108 nodes to 96) |
| another number of predictors | "malformed cut points in warm-start donor" |
| numeric against an unordered factor | the predictor structure error |
| numeric, or an ordered factor of 6 levels, into an ordered factor of 4 | taken by position (6 of 7 and 7 of 7 splits on the column kept, p2.R); `setState` of the same state is an error in both forms |
| an ordered factor of 4 into a numeric column | taken; 4 of 5 splits on the column kept, all at 0.99, the recipient's top cut point. `setState`: FALSE; forced, 0 of 5 kept (p2.R) |
| levels a, b, c against c, b, a | taken; correlation with the donor's fit -0.998 |
| another tree count, either way | the shape error |
| DART on one side only | the DART error; both DART: taken |
| a variance forest on one side only; both, with another count | the shape error; the variance count error |
| constant into linear or gp, linear against gp, linear with another number of covariates, gp on another number of rows | "malformed leaf parameters in bartcore state" |
| gp on other rows of the same number | taken, every fit at one value |
| linear into constant | taken, NULL invisibly, leaf values the donor's intercepts exactly |
| gp into constant | taken, NULL invisibly; starting range -1.95 to 4.99 against -2.15 to 5.43 |
| linear on x2 into linear on x3 | taken |
| linear, the covariate at one value on one side only | taken; the start is the donor's function at the recipient's rows to 4e-15 (f-linear-one-value.R) |
| linear, the covariate at 0 for the donor and at 3 for the recipient | taken, coefficients as the donor built them; the start differs from the donor's function at the recipient's rows by up to 9.1 (sd 2.55; f-linear-one-value.R) |
| `interactions(max.order = 1)`, `blocks`, a restricted variance forest | each its error. `setState` of the same trees: FALSE; forced, 108 nodes to 92 and to 58, the fit's mean over the rows unchanged |
| free donor into a monotone recipient | taken; 11 of 30 trees with every leaf 0 (14 of 30 in p2.R). `setState`: FALSE |
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
  same-grid `installTrees`, then `setCutPoints(own, splits = "value")`. On the critique's fixture (a donor
  of 182 nodes, scratch/wsc/p1.R) that gives exactly the routine's trees at 30, 60 and 300 rows (166, 180
  and 180 nodes, the grid back identical). Today's arm gives 174 at 30 rows and the same trees as the
  routine at 60 and 300, so only the 30-row fixture tells the two arms apart.
- With the recipient's grid written into the donor state, which reaches the same-grid arm on sgi, stored
  states and the next 5 draws are identical to today's remap arm on 8 fixtures (scratch/wsc/p6.R): constant
  leaves on a finer grid, 40 rows with empty leaves, monotone, sides drawn, sides dropped, linear on a finer
  grid, linear under another standardization, heteroscedastic. A gp pair differs: the same-grid arm copies
  the fits and the remap arm starts them at 0.
- Section 4's staging of a donor on another grid is not needed: across grids the unforced answer is FALSE by
  choice (b) with no tree staged, and the forced install is the routine.

## The rule

The first row that applies. "Error" leaves the sampler untouched in both forms.

| The donor against the recipient | Unforced | Forced | Slice |
|---|---|---|---|
| not a donor; an older encoding; a malformed block or tree; a donor grid that repeats a point; a donor threshold that is no point of the donor's own grid; a bad `samples`; a recipient of two or more forests | error | error | as today |
| another number of predictors; a column an unordered factor on one side only; a pair of leaf kinds that is an error today (constant into linear or gp, linear against gp, linear with another number of covariates, gp on another number of rows) | error | error | A keeps today's texts; D rules the kinds (below) and writes the texts of (d) |
| on a column the recipient holds as an ordered factor, a donor grid other than the recipient's: a numeric donor column, or an ordered factor of another level count | error | error | A; taken by position today. The predictor structure text in A, a text naming the column in D |
| every tree stands as stored; sides drawn under dec-B440 and linear coefficients restated under dec-B328 and dec-B329 are how the donor is read | TRUE | NULL; today's install bit for bit | A |
| a split whose threshold is no point of the recipient's cut grid, a donor column held as an ordered factor against a numeric one here included | FALSE | NULL; the across-grids move | A |
| a leaf no row reaches; a split outside its interval; a side for missing values on a column that has never held one here | FALSE | NULL; merged or dropped, as today | A |
| a split on a column its forest may not use, or past its `interactions` limit | FALSE | NULL; cut back at the first such split from the root | A |
| a monotone tree with leaf values out of order | FALSE | NULL; in A every leaf 0, as today; from B collapsed where the order breaks | A, B |
| gp leaves on both sides, the donor recorded on other cut points | FALSE | NULL; trees taken, fits started at 0, as today | A |
| leaf blocks the recipient does not read: a linear or gp donor into constant leaves | FALSE | NULL; as today, each leaf at the value the donor stores for it, a linear leaf's intercept | A; D makes the gp pair an error |
| linear leaves on both sides, a covariate that holds one value in the donor's data and another single value here, under a slope that is not 0 | FALSE | NULL; as today, the coefficients as the donor built them (dec-B329) | A |
| a recipient of more trees, mean or variance | TRUE | NULL; all taken in the donor's order, the rest single leaves at 0, or at 1 in a variance forest | C |
| a recipient of fewer trees | FALSE | NULL; the trees that vary most | C |
| constant donor into linear leaves | TRUE; slopes 0 | NULL | D |
| linear donor into linear leaves on another number of covariates | FALSE | NULL; each leaf at its stored intercept, slopes 0 | D |
| a gp leaf against another kind; a gp donor on another number of rows | error | error | D, by name |
| DART, or a variance forest, on one side only | not judged | not judged; the trees taken, the donor's part left behind, the recipient's started as a new fit's when the donor goes in | D |

- Values are as `setState`'s ([install-surface.md](../design/install-surface.md)): unforced TRUE and FALSE
  are returned visibly; forced returns NULL invisibly whether or not anything was changed.
- FALSE leaves the trees, drawn values, generators, kept draws, DART weights, variance forest and the `state`
  field as they were. The sides drawn under dec-B440 are drawn on scratch first, and a decline puts the
  generators back. After an install, taken or forced, each chain's generator stands past its coins.
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
- A donor column held as an ordered factor against a numeric column here is moved like any other: the
  routine re-quantizes a numeric column against any grid. The donor's thresholds there lie between level
  codes, so on a column whose values lie outside them every such split has a side no row reaches and merges
  (0 of 5 kept in the run), as a forced `setState` leaves it. The other direction is the error above: the
  level count fixes an ordered factor's grid, and the routine cannot re-quantize it.
- A merged leaf is the mean of the leaves it replaces weighted by the recipient's rows in each under its
  working weights, geometric in a variance forest, the plain mean where none holds a row (dec-B398).
- Several chains: one chain's tree that needs a change declines the whole donor.

### Altered

`installForests` reports `altered` as [`Sampler::setState`](../../src/bartcore/sampler.hpp) does: set by the
installing build, never by a scratch build, and false after an error or a decline. The unforced form returns
FALSE exactly when the forced form of the same call sets it.

Set, the donor not installed as it holds its trees:

- a tree changed: a leaf no row reaches, a split outside its interval or a forbidden split merged
  (`rebuildLiveForest`, `rebuildVarianceForest`); a split moved or merged across grids, which every install
  across grids does; a monotone tree out of order reset (A) or collapsed (B); from C, a donor tree left out;
- a side changed: a missing-value side dropped on a column that has never held a missing value here;
- a fit changed: gp fits started at 0 because the donor is recorded on other cut points;
- a leaf value not the donor's function: leaf blocks the recipient does not read; coefficients kept across a
  covariate that holds one value on each side and not the same one, under a slope that is not 0; from D, a
  linear leaf reduced to its intercept.

Not set, each being how a donor is read and not a change to it:

- a side drawn under dec-B440, though a merge or a collapse that a drawn side causes is set as above;
- linear coefficients restated in the recipient's standardization where the leaf stays the same function of
  the raw covariates, a covariate that holds one value on one side only included (dec-B328, dec-B329; the
  start equals the donor's function at the recipient's rows to 4e-15 in the run);
- sigma, k or a DART concentration the recipient holds fixed (dec-B343); another response family or scale
  (dec-B418); the recipient's kept draws dropped; a pooled categorical rule keeping a missing bit its column
  no longer routes, as for `setState`;
- from C, single leaves added to fill a larger forest; from D, slopes of 0 under a constant donor, and a
  part only one side has, which is not judged.

### The line

`bart` prints it after its setup lines when `verbose = TRUE`, at every warm start, and nothing otherwise; no
warning is raised. Trees are counted per chain and changes summed over chains. The words are the design
note's (f) where it gives them.

    warm start: 75 of 75 donor trees taken
    warm start: 50 of 50 donor trees taken; 12 moved to this fit's cut points, 21 cut back to fit, over 4 chains
    warm start: 30 of 30 donor trees taken; 14 reset to 0
    warm start: 30 of 30 donor trees taken; 3 with a side for missing values dropped; variance trees: 6 of 6 taken, 2 cut back to fit
    warm start: 8 of 8 donor trees taken; gp fits started at 0

- Each taken tree is counted once, under the first that holds: "reset to 0" (a monotone tree whose leaves
  are all 0 where the donor's were not; slice A only), "cut back to fit" (it lost a split), "moved to this
  fit's cut points" (it kept every split and holds one on another point), "with a side for missing values
  dropped". A tree cut back and then reset is one reset tree.
- The counts follow the first clause in the order moved, cut back, side dropped, reset, joined by ", ", and
  ", over K chains" closes them where there are several chains.
- A variance forest has its own clause, "variance trees: T of N taken" with its own counts, printed where
  one of its trees changed or, from C, was left out.
- Whole-donor phrases, each its own clause: "gp fits started at 0"; "slopes left behind" and "gp fits left
  behind" (leaf blocks the recipient does not read); "slopes kept across a covariate of one value". Slice C
  adds "150 trees started as single leaves"; slice D "slopes started at 0", "DART weights started fresh",
  "donor's DART weights left behind", "variance forest started fresh" and "donor's variance forest left
  behind".
- A donor that is not clean never prints the first clause alone. Slice B retires "reset to 0": a monotone
  tree collapsed is cut back.

## Constraints

- A forced install from a donor whose live thresholds are all points of the recipient's grid runs
  [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) and
  [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp) with the arguments they get today and
  `altered`, and from slice B the argument that names the monotone repair.
- The unforced answer is FALSE exactly when the forced install is altered ([Altered](#altered)). The verdict
  stages a tree by the routine the install stages it by. tests/cpp asserts declined equal to altered on
  every fixture; R compares trees only where no side is drawn and no slope restated.
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) returns what it returns today for every state; its
  tree loops move into helpers unchanged.
- [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp) and
  [`Sampler::setState`](../../src/bartcore/sampler.hpp) are called, not edited.
  [`Sampler::installAcrossGrids`](../../src/bartcore/sampler.hpp) already restores on a throw, and A calls
  it unedited; B gives it, with [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp), the argument
  that names the monotone repair, whose default is the reset every other caller keeps.
- `setState`, `setData`, `setCutPoints` and a forced `setPredictor` keep the monotone reset, bit for bit
  (dec-B278, dec-B355). Only the warm start collapses.
- No change to a stored state's format, and no flat entry of
  [dbarts.h](../../inst/include/dbarts/dbarts.h) installs a donor, so the C API hash does not move.
- `installTrees` stays refused on a borrowed view: the forced move re-quantizes.
- The facade virtual changes in A and D: install with `--preclean`, as on every engine commit.
- Errors. Slice A lifts two, the column restriction and the interaction limit, which become FALSE and a
  cut-back, and adds one, the ordered factor against another grid, which installs by position today. Slice
  B touches none. Slice C lifts the two tree-count errors. Slice D lifts those for DART or a variance forest
  on one side, for constant into linear leaves and for linear leaves on another number of covariates, and
  adds one: a gp donor into constant leaves, which A declines and force installs. Every other error of
  today stays an error in both forms.

## Slice A: the two forms

### Steps

1. Staging shared with the state install (about 150 lines, [chain.hpp](../../src/bartcore/chain.hpp)).
   - Move the body of each per-tree loop of `stateIsValid`, mean and variance, into a private member, and
     call them from `stateIsValid`.
   - Add `donorTreesStage(state, notClean, varianceFailed)`: the two helpers over every live tree of a warm
     start's assembled chain. False where a tree does not build, `varianceFailed` saying which forest, so the
     two existing messages stand. With `notClean` null the trees are only built.
   - Add `coldStartFunctionFits(state)`, compiled for function leaves: every live tree's fits set to zeros.
   - `convertDonorStandardization` gains an out flag, set where it leaves a column as built that holds one
     value on both sides, another center here, under a slope that is not 0.
   - Add `countDonorTrees(assembled, counts)`: flattens the live trees
     ([`Chain::getState`](../../src/bartcore/chain.hpp) without the kept draws) and compares each with the
     tree the install read, the assembled scratch after the direction draw and the restatement, so neither
     is counted. Per tree, the first that holds: in a monotone forest every leaf 0 where the assembled
     tree's were not is reset; fewer nodes is cut back; equal nodes with a threshold changed is moved; equal
     nodes and thresholds with a missing-value flag changed is a side dropped.
   - `installForest` and `installVarianceForest` take `altered` and hand it to `rebuildLiveForest` and
     `rebuildVarianceForest`; they lose their grid and store parameters, the second its remap arm.
   - Delete `interactionStateFeasible`, `columnMaskStateFeasible`, `forestsRebuildable`,
     `varianceForestRebuildable` and `rebuildLiveForestRemapped`.
2. The install (about 280 lines, [sampler.hpp](../../src/bartcore/sampler.hpp),
   [facade.hpp](../../src/bartcore/facade.hpp)).
   [`Sampler::installForests`](../../src/bartcore/sampler.hpp) becomes
   `installForests(donor, sampleMap, currentPredictors, force = true, altered = nullptr, report = nullptr,
   leafBlocksUnread = false)`. `currentPredictors` has no default: a dense column re-quantizes from it and a
   null there is read. `leafBlocksUnread` is the bridge's word that the donor carries leaf blocks this
   sampler's leaves do not read; slice D replaces it by the donor's shape.
   - The checks, the assembly and the direction draw stay as they are. The grid check also refuses a donor
     grid other than the store's on a column the store holds as an ordered factor
     (`WarmStartResult::gridMismatch`); a grid equal to the store's passes.
   - Where the donor's recorded points differ from the store's and the leaf is a function leaf, the fits are
     cold-started on the scratch.
   - A live threshold that is no point of the donor's own grid is the rebuild error. `acrossGrids` is true
     where some chain holds a live threshold that is no point of the store's grid.
   - Every chain's trees are built before any answer is given, in both forms: across grids under
     [`ScopedCutGrid`](../../src/bartcore/data.hpp) on the donor's points, for validity only; otherwise by
     `donorTreesStage` on the store's grid, with the verdict where the call is unforced. A donor that does
     not build is the rebuild error, never a decline.
   - Not clean is any of: `acrossGrids`; the verdict's flag; fits cold-started; `leafBlocksUnread`; the
     restatement's flag. Not clean and unforced returns a new `WarmStartResult::declined` with nothing
     written; the generator guard puts back any draw.
   - The install step is today's two loops on their same-grid arms, reading the slot of the donor that
     `sampleMap` names for each chain. Across grids it is handed to `installAcrossGrids`; otherwise it is
     called. Both arms install sigma, k and the scalars as today and drop the kept draws.
   - The undo. Every chain is snapshotted first, on the store's grid. Across grids every chain counts as
     touched from the start, and the restore runs only after `installAcrossGrids` has put the store's grid
     back, never inside the install step: the snapshots' thresholds resolve on the store's grid alone.
     `installAcrossGrids` puts the store's grid back when its install step throws, and rethrows.
   - `altered`, when given, is false on every exit but an install, and then true where the donor was not
     clean by any of the five above; the same-grid arm reads it from the rebuilds.
   - `report`, when given, is filled per chain after the install: offered and taken; moved, cut back, side
     dropped and reset from `countDonorTrees`, for the mean forest and for the variance forest; and the
     flags fits cold, leaf blocks unread and coefficients kept across one value.
   - [`WarmStartResult`](../../src/bartcore/sampler.hpp) loses `interactionMismatch` and
     `columnMaskMismatch`.
   - Rewrite the comment block above `installForests`. The virtual and its override take the new signature.
3. Bridge and R (about 170 lines).
   - [`bartcore_installForests`](../../src/R_interface_bartcore.cpp) takes the predictor matrix and the
     force flag as a fourth and fifth argument (its row in [R_interface.cpp](../../src/R_interface.cpp), its
     declaration in the bridge header) and returns a list: `installed`, `altered`, one integer matrix of
     counts per forest kind, a row per chain, and the flags.
   - [`readWarmStartState`](../../src/R_interface_bartcore.cpp) notes, for a recipient of constant leaves, a
     donor forest that carries `tree.params`: about 10 lines, and the flag above.
   - retired: [`columnMaskMismatchMessage`](../../src/R_interface_bartcore.cpp) and the two cases that read
     the deleted results go, and went.
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
   - First run a tracer over the 23 files on the slice's library that prints each executed call's unforced
     answer with its file and line. The rule is per call, not per line: a line that a loop takes through
     both answers ([test-state-missing-record.R](../../inst/tinytest/test-state-missing-record.R) and
     [test-monotone-unforced.R](../../inst/tinytest/test-monotone-unforced.R) each hold one) takes
     `forceUpdate` from what that pass means. Add `forceUpdate = TRUE` where a test means the repair, and
     leave a call that returns TRUE alone.
   - The critique's tracer on sgi (scratch/wsc/tracer.R, which approximates the verdict by comparing trees)
     finds 15 to 17 executed calls that turn FALSE, in 8 files:
     retired: ["column restriction"](../../inst/tinytest/test-blocks.R) 1, an expectation since rewritten,
     [test-single-forest-vars.R](../../inst/tinytest/test-single-forest-vars.R) 1,
     [test-heteroscedastic-warm-start.R](../../inst/tinytest/test-heteroscedastic-warm-start.R) 2,
     [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) 2 to 4,
     test-monotone-unforced.R 4, [test-monotone.R](../../inst/tinytest/test-monotone.R) 1,
     [test-state-empty-leaf-merge.R](../../inst/tinytest/test-state-empty-leaf-merge.R) 1 and
     test-state-missing-record.R 3 or more (the tracer cannot tell a dropped side from a drawn one); and
     retired: ["interaction constraint"](../../inst/tinytest/test-interactions.R), since rewritten, through
     `bart`, which now fits.
     These nine are the named files.
   - An expectation in a named file that pins today's cross-grid values is changed only to the by-hand
     twin's value. 6 executed calls of test-leaf-conversions.R give a donor on other recorded points with
     every threshold on the grid: they move from the remap arm to the same-grid arm and are expected
     unchanged.
   - The 3 raw calls gain two arguments.
   - In tests/cpp the 44 calls gain the predictor argument. The 2 `interactionMismatch` and 2
     `columnMaskMismatch` expectations become `declined` unforced and a checked cut-back forced.
     [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp) checks the merged leaf by rows.
     [`testWarmStartStandardization`](../../tests/cpp/test_state.cpp) keeps its gp sentence, "copied on the
     same grid, cold on another".
   - [composition-matrix.R](../../benchmarks/R/composition-matrix.R) has 2 calls to read.
5. New tests (about 650 lines; Tests).
6. Help, NEWS and notes (about 180 lines).
   - [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the usage line; `forceUpdate` gains "For
     `installTrees`, FALSE (the default) takes the donor only if every tree it gives can stand in this
     sampler as the donor holds it, and TRUE takes it changing what must be changed"; `donor` drops "tree
     count, and DART setting" only in slices C and D.
   - Warm starts, replacing the sentences on refusal and remapping: "`installTrees(donor)` asks first. A
     donor every tree of which fits this sampler as stored is installed and the call returns TRUE. One that
     would have to be changed is not installed: the call returns FALSE and the sampler is as it was, its
     trees, drawn values, generators, kept draws and `state` field. A tree has to be changed when a split's
     threshold is not one of this sampler's cut points, as for a donor fitted to other rows; when a leaf
     holds none of this sampler's rows; when a split is on a column, or in an interaction, its forest does
     not allow; when it records a side for missing values on a column that has never held one here; or when
     a monotone forest's leaf values are out of order. A donor is also declined when it holds slopes or
     Gaussian-process fits that a sampler of constant leaves does not read, when its Gaussian-process fits
     were recorded on other cut points, and when a linear leaf's covariate holds one value in both fits, a
     different one in each. `installTrees(donor, forceUpdate = TRUE)` installs such a donor and returns
     NULL invisibly, whether or not anything was changed: a split off the cut points moves to the nearest
     one inside the range the splits above it leave or, with none left, becomes one leaf with everything
     beneath it; a split with a side no row reaches, and the first split from the root that its forest does
     not allow, likewise become one leaf, at the mean of the leaves it replaces weighted by the rows each
     held; a side is dropped; a monotone tree out of order has every leaf set to zero; Gaussian-process
     fits start at zero; and leaf values the sampler does not read are left behind. `bart`'s `warm.start`
     always forces."
   - [bart.Rd](../../man/bart.Rd), `warm.start`: "takes all of the earlier fit that this one can hold: a
     split that is not on this fit's cut points is moved to the nearest one, a tree this fit cannot hold as
     it stands is cut back at the split that does not fit, and a tree that breaks a monotone constraint
     starts with every leaf at zero. With `verbose = TRUE` one line after the setup lines says how many
     trees were taken and how many were moved, cut back or reset."
   - [interactions.Rd](../../man/interactions.Rd): "enforced on any warm start" becomes "a warm start cuts
     back a tree that breaks it".
   - [monotone.Rd](../../man/monotone.Rd): "`installTrees`" in the list of calls that always complete
     becomes "`installTrees` with `forceUpdate = TRUE`".
   - [NEWS.Rd](../../inst/NEWS.Rd), the Starting points item: "`bart`'s `warm.start` starts the chains from
     a previous fit's trees, taking all of the earlier fit that the new one can hold, and with `verbose =
     TRUE` says in one line how many trees were taken, moved and cut back ... the sampler offers both as
     `$installTrees` and `$growFromRoot`; `installTrees(donor)` returns TRUE for a donor taken as stored and
     FALSE, the sampler untouched, for one that would have to be changed, which `forceUpdate = TRUE`
     installs." No sentence compares branch builds.
   - [install-surface.md](../design/install-surface.md): the sentences "A warm start is not a state install"
     and "A warm start still moves", and the closing section, are rewritten to the rule.
     [architecture.md](../architecture.md): the containment sentence.
   - The five deleted routines are cited 13 times in 10 other documents, 3 times in the design note and in
     this plan's Context (by search over the link form): a frozen record gets the `retired:` marker and a
     clause saying the routine is gone, a present-tense document a rewritten sentence.
   - TODO: drop `forced-update-returns-null`; `warm-start-salvage` keeps what slices B to D build.

Commits, each on a tree that passes: (1) the two forms: steps 1 to 4 less the report; (2) the report and the
line: `countDonorTrees`, `report`, the bridge's counts, `warmStartLine` and `bart`'s print; (3) step 5; (4)
step 6. Commit 2 is the split point named under Budget and stops.

### Tests

tinytest in a new `test-warm-start-forms.R`, tests/cpp in `test_state.cpp`. Each row: what it pins; the
mutant it must kill. "m" and a number is the critique's missing mutant of that number.

| Test | Pins | Mutant |
|---|---|---|
| clean is TRUE | same rows: TRUE visibly, trees and fit the donor's, the `state` field stored again, kept draws dropped | the verdict always FALSE |
| default declines | donor on 40 of its rows and its cut points: `installTrees(donor)` is FALSE visibly | `forceUpdate = TRUE` as the default |
| forced returns NULL | the same donor forced: NULL invisibly, trees equal a forced `setState`'s; forced on a clean donor: NULL too | TRUE, or a visible value, when forced |
| declined leaves all | a recipient run after its last `storeState`, so that its `state` field is stale; after FALSE: trees, sigma, k, the next 5 draws of a twin that never called, `predict` from kept draws, and the `state` field identical to the stale one | kept draws dropped; the state stored before the decline |
| off the grid declines | donor on other rows: FALSE; forced, the recipient's cut points, weights and codes identical before and after | the verdict not reading thresholds |
| equals the two calls | 30 rows, the critique's fixture, where sgi's remap arm gives 174 nodes and the routine 166, with constant and with linear leaves: forced, the trees, leaf values and next 5 draws equal a twin's that took the donor's grid, installed and called `setCutPoints(own, splits = "value")`; and a forced `setState`'s trees on an unweighted pair | the remap arm kept; a linear leaf's cached statistics left stale across grids |
| finer grid is clean | a recipient grid holding the donor's points and others: TRUE, splits as stored | the verdict keyed on the grids differing |
| merged leaf by rows | the hand tree of install-surface.md (leaves 0.10, -0.20, 0.05 holding 73, 38, 89 rows) as a donor state forced onto the 50th and 75th points: -0.0027 | the plain mean (-0.05); rows counted after the move (0.10) |
| weights | a case-weighted recipient: a merged leaf weighs by its working weights, by hand | unit weights |
| cut back, same grid | `interactions(max.order = 1)` and `blocks` recipients: FALSE; forced, no forbidden split left, the splits above stay, one tree's merged leaf by hand, the fit's mean over the rows unchanged | the error kept; the leaf at 0 |
| cut back, across grids | the same on other rows | the predicate unarmed on the across-grids arm |
| variance forest | a variance threshold off the grid alone: FALSE; a split a restricted variance forest may not use: FALSE, forced cut back at the geometric mean | variance trees not judged |
| monotone | free donor: FALSE; forced, the fit is in order along the monotone columns and the trees in order as stored are the donor's. What an out-of-order tree becomes is slice B's | the monotone test left out of the verdict |
| missing sides | donor complete, recipient not: TRUE with sides drawn, and each generator past its coins (its bytes in the stored state differ from those before the call); a decline after a draw leaves the generators (twin draws equal); a donor side on a column never missing here: FALSE, forced dropped | the generator not put back; a dropped side not judged; m5, the generator rewound after an install |
| verdict after the draw | tests/cpp: a hand tree whose lower split is on a column missing here, with one side that only rows missing it can reach; over seeds, the unforced answer is TRUE on the seeds where a forced twin of the same seed keeps that split and FALSE where it merges it, both seen | m1, the verdict taken before the direction draw |
| the slot `samples` picks | a donor whose live trees are clean and whose named kept draw needs a change: FALSE; the reverse: TRUE; and with two chains, chain 2 given the draw that needs a change: FALSE, chain 1 untouched | m2, the verdict on the donor's live trees; the first chain only |
| malformed is an error | tests/cpp: a donor with a threshold off the recipient's grid and a tree that does not build: the rebuild error in both forms, the sampler untouched | m3, the across-grids decline returned before the build |
| declined equals altered | tests/cpp, over a hand-built fixture for every row of the rule and 40 seeded pairs: unforced `declined` exactly when a forced twin reports `altered`, and `altered` exactly when the report holds a count or a flag. R, on the pairs with no side drawn and no slope restated: FALSE exactly when the forced trees differ from the donor's | any gap between verdict and install |
| gp across cut points | same rows, the recipient's grid the donor's points and the midpoints between them, so every threshold is on it: FALSE, forced fits 0; the same grid: TRUE, fits copied | fits copied across grids; m4, the cold start not marked |
| leaf blocks unread | linear into constant and gp into constant: FALSE; forced, the leaf values sgi gives | TRUE unforced |
| one value | linear, the covariate at 0 for the donor and 3 here: FALSE, forced the donor's coefficients; at 0 on both: TRUE; at one value on one side only, each way: TRUE, the start equal to the donor's function at these rows to 1e-12 | every restatement declined; the two-sided case judged clean |
| ordered factor | a numeric donor column, and an ordered factor of another level count, into an ordered factor: error in both forms, state identical; the same levels on both sides: TRUE; an ordered factor into a numeric column: FALSE, forced, the splits with a side no row reaches merged | the routine run on a factor column; m11, a recipient with equal grids refused |
| scalars across grids | forced on other rows, a recipient that draws sigma and k: both the donor's after | m6, sigma and k not installed on the across-grids arm |
| kept draws across grids | forced on other rows, a recipient holding kept draws: none recorded after | m12, kept draws not dropped on the across-grids arm |
| argument | `1`, `NA`, `c(TRUE, TRUE)` and `"partial"` are errors | coercion |
| errors stay | another number of predictors, factor against numeric, tree count, DART and variance one side, constant into linear, constant into gp, linear against gp each way, linear with another number of covariates, gp on another number of rows, two forests, a malformed donor: error in both forms, stored state identical | an error lifted by force |
| bart forces and says | `bart(warm.start = )` into an interaction-limited fit runs; one line starting "warm start:" before "Running mcmc loop:"; nothing without `verbose`; no line from `bart(verbose = TRUE)` with no warm start; no warning, counted with `withCallingHandlers` | `bart` unforced; the line missing or after the loop; m10, the line with no warm start |
| line counts | donor states built by hand with the counts known from how they are built: four trees, one as stored, one moved, one cut back, one with a side dropped; two chains given two kept draws by `samples`, with other counts in each; a monotone tree cut back and reset; a variance forest with one tree cut back; each line exact | cut back counted as moved; m7, taken summed over chains or changes from chain 1 only; m8, the tree counted twice; m9, counts against the live trees and not the named draw |
| a declined donor says why | every hand-built fixture that is not clean: its forced line holds a clause beyond the first | "N of N donor trees taken" alone for a donor that returns FALSE |
| sparse and mixed stores | forced across grids on a sparse matrix and a mixed frame: grid and codes identical after | the predictor pointer read where it is null |

tests/cpp besides: `declined` leaves each chain's state and generator bytes identical; the report's counts on
a hand-built donor; the 4 rewritten expectations; the cross-grid merged leaf.

Not in the table: a same-grid donor sent through the across-grids routine. Twins on one build take one path,
so no test in the suite sees it; the cross-build script of Verification is its only oracle.

## Slice B: a monotone tree out of order is collapsed

A forced warm start into a monotone recipient collapses a donor tree where its leaf values break the order,
until [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp) holds. Unforced such a donor is FALSE, as in
slice A. `setState`, `setData`, `setCutPoints` and a forced `setPredictor` keep the reset.

### The collapse

Two leaves are related only across the lowest split that separates them, and only where that split is on a
monotone column: the two sides of a split share no value of its column, so they cannot be related along any
other. So the order can break only at a split on a monotone column, between a leaf of its left side and a
leaf of its right.

The tree is walked from the leaves up, a split after both its sides. At a split on a monotone column, if a
leaf on its left and a related leaf on its right are out of order, the first of these that puts every such
pair in order is taken:

1. one side merged into one leaf, the side holding fewer nodes first and the left on a tie;
2. the other side merged;
3. both sides merged, the split kept between two leaves;
4. the split itself and everything beneath it merged.

Each merge is one leaf at the mean of the leaves it replaces weighted by the rows each held, under the
working weights, as every merge of the install is (dec-B398). A side that is one leaf already is skipped.

Tree 7 of the fixture below, monotone increasing in x1 and x2, values as `getTrees` gives them:

    as stored                                   collapsed
    x2 <= 0.671                                 x2 <= 0.671
      yes: x3 <= 0.762                            yes: leaf -0.0590  (198 rows)
             yes: leaf -0.0670  (150 rows)        no:  x2 <= 0.750
             no:  leaf -0.0339  (48 rows)                yes: leaf -0.0471
      no:  x2 <= 0.750                                   no:  leaf  0.0366
             yes: leaf -0.0471  (24 rows)
             no:  leaf  0.0366  (78 rows)

The split on x3 is on a free column and is not tested. The lower split on x2 is in order (-0.0471 below
0.0366). At the root the two left leaves touch the leaf at -0.0471 along x2, and -0.0339 above -0.0471 is
out of order. Both sides hold three nodes, so the left goes first:
(150 x -0.0670 + 48 x -0.0339) / 198 = -0.0590, below -0.0471, and the root stays. Merging the root itself
would leave one leaf at -0.0332. In tree 4 a split on x1 separates two leaves, -0.0207 on 15 rows below and
-0.0233 on 116 above: neither side can be merged, so the split becomes one leaf at
(15 x -0.0207 + 116 x -0.0233) / 131 = -0.0230.

What it guarantees:

- It ends: one pass, each split visited once, at most four candidates tested at it.
- The result is in order. A related pair of the result is separated by one split on a monotone column; that
  split was tested with both its sides final and kept only with every pair across it in order, and a merge
  above it removes the pair.
- It does not depend on the order in which sides are walked: whether the pairs across a split are in order
  depends on the splits above it, which stand while it does, and on what is beneath it.
- A tree in order as stored is left bit for bit. Nothing is changed outside a split that separates an
  out-of-order pair and what is beneath it.
- Each tree's mean over the rows is unchanged.
- It is not the least cut in every tree. A side holding two or more splits is merged at its top, where
  merging one of its lower splits alone might do. The least cut is a search over every pruning, 5, 26, 677
  and 458330 of them for full trees of depth 2 to 5, since a merged value changes which pairs are in order.

Measured by a prototype in R on `getTrees` output from sgi (scratch/wsc/f-monotone-collapse.R), whose
predicate fails exactly the trees the engine sets to 0. The critique's fixture: 30 trees, 134 nodes, 300
rows, monotone in x1 and x2; 14 trees out of order, holding 72 nodes.

| what a forced start does with the 14 | nodes kept with their values | trees holding a split with its values | rmse to the donor's fit (sd 0.190) | the fit's mean over the rows |
|---|---|---|---|---|
| every leaf 0 (today) | 0 of 72 | 0 of 14 | 0.151 | moves by -0.050 |
| step 4 alone, the split itself | 22 | 3 | 0.123 | unchanged |
| the collapse | 32 | 6 | 0.113 | unchanged |
| the least cut, by search | 32 | 6 | | unchanged |

The reset keeps every split and no value; the collapse keeps the values and the splits that can stand. Over
the whole forest 94 of 134 nodes carry donor values after the collapse, against 62 after the reset.
The collapse keeps as many nodes as the least cut on each of the 14 trees, and on each of the 9 and the 20
out-of-order trees of the same donor with x2 alone and with all four columns monotone (23 and 44 nodes
kept, against 13 and 34 for step 4 alone). Unit row weights, one donor.

### Steps

1. The merge as a routine (about 15 lines, [tree.hpp](../../src/bartcore/tree.hpp)): the body of
   `collapseEmptyNodesBelow` that turns a node and everything beneath it into one leaf at the weighted mean
   becomes a member that it calls, the sums in the order they run today.
2. The collapse (about 75 lines, [model.hpp](../../src/bartcore/model.hpp), beside
   `monotoneTreeIsFeasible`): the walk above, over a tree whose rows are routed. A candidate is tested
   without changing the tree, by [`MonotoneLeafGeometry`](../../src/bartcore/model.hpp) built over the
   nodes that would be its leaves, with the predicate's relation and its tolerance. Returns whether it
   merged anything.
3. Reaching it (about 40 lines). `rebuildLiveForest` and `forceRefreshTrees` take an argument naming what a
   monotone tree out of order gets, the reset by default and the collapse from the warm start;
   `installForest` passes the collapse, and `installAcrossGrids` hands the argument to the refresh, so a
   tree the move across grids leaves out of order is collapsed too. The collapse sets `altered`. The verdict
   does not change.
4. The report loses `reset` and the line "reset to 0"; such a tree is counted as cut back (about 15 lines,
   bridge and R).
5. Existing tests: the expectations of test-monotone.R and test-monotone-unforced.R that pin every leaf 0
   after a warm start (by search before editing; 12 of 20 trees in one of them) become the collapse's; those
   that pin the reset after `setState`, `setData`, `setCutPoints` or `setPredictor` do not move.
6. New tests (about 180 lines; Tests).
7. Help, NEWS and notes (about 50 lines): the two help sentences of slice A on a monotone tree become "a
   monotone tree whose leaf values are out of order is cut back where they are"; monotone.Rd says that a
   forced warm start collapses where `setData`, `setCutPoints` and a forced update reset; the design note's
   section 2 and its (f) phrase "3 reset to 0"; install-surface.md's row on the monotone tree for a warm
   start.

Commits: (1) steps 1 to 5; (2) step 6; (3) step 7.

### Tests

tests/cpp on hand trees in `test_state.cpp`, tinytest in `test-warm-start-forms.R`.

| Test | Pins | Mutant |
|---|---|---|
| the split itself | two leaves out of order across a monotone split, 15 and 116 rows: one leaf at -0.0230 | the plain mean (-0.0220); every leaf 0 |
| a side first | tree 7 above: the left side one leaf at -0.0590, the root and the right side as stored | step 4 with no side tried |
| the smaller side | a tree where merging either side restores the order and the sides differ in size: the smaller is merged | the larger side first |
| both sides | a tree where neither side alone restores the order and both do: the split kept between the two means | step 3 skipped |
| from the leaves up | a tree where a lower merge puts an upper split in order: the upper split stays | the root tested first |
| unrelated leaves | two leaves whose values are out of numeric order and that share no value on another split column: untouched | every left leaf tested against every right leaf |
| a free column | values out of order across a split on a free column: untouched | the column's direction not read |
| in order after | 200 seeded random trees and values: `monotoneTreeIsFeasible` holds after, each tree's weighted mean over the rows is unchanged to 1e-12, and a tree in order before is identical | a tolerance other than the predicate's; a merge not weighted by rows |
| weights | a case-weighted recipient: the merged leaf by hand | unit weights |
| the warm start | the critique's fixture: FALSE unforced; forced, no tree with every leaf 0, the trees in order as stored the donor's, each tree's mean over the rows the donor's, the fit in order; the line says "14 cut back to fit" and no "reset" | the reset kept on the warm start |
| across grids | the same donor forced onto other rows: no tree with every leaf 0, the fit in order | the refresh left on the reset |
| the others reset | `setState` forced, `setCutPoints` and `setData` on a tree out of order: every leaf 0, as the existing tests of test-monotone.R pin | the collapse as the default |
| declined equals altered | slice A's tests/cpp check with the monotone fixtures added | the collapse not reported |

On a monotone recipient the by-hand twin of "equals the two calls" ends in `setCutPoints`, which resets, so
that test holds no monotone recipient and "across grids" above stands alone.

## Slice C: another tree count

### Steps

1. Tree count, both forests (about 260 lines, sampler.hpp and chain.hpp).
   - The donor forest's own count sets the stride of a kept draw; the two count refusals go.
   - The sides of dec-B440 are drawn once, over every tree of the donor sample a chain takes, mean then
     variance, before anything else; every round and the final install read those sides. After a forced
     install each generator stands past one coin for every rule on a raised column in the whole sample,
     chosen or not; a decline puts it back.
   - More trees: the assembled chain is the donor's trees in order, then single leaves at 0 (slopes 0; gp
     fits 0; 1 in the variance forest). Clean.
   - Fewer trees: not clean. Forced, the donor's trees are tried in the recipient's live slots, the
     recipient's count at a time and in the donor's order, by the install of slices A and B on either arm.
     The last round holds the trees left over and then single leaves, written afresh.
   - After each round `liveTreeVariation` reads each tree's variance over the rows from its fits as
     installed ([`muByTree`](../../src/bartcore/combiner.hpp) through the leaf map, or the fit slab);
     `liveVarianceTreeVariation` reads the same of the log of
     [`factorByTree`](../../src/bartcore/chain.hpp). A row weighs by its case weight, and a row out of the
     likelihood, of weight 0 or switched off by its family's mask
     ([`isRowIn`](../../src/bartcore/model.hpp)), is left out. The working weights are not read: dec-B323
     says "over the recipient's rows", and a ranking does not depend on a latent draw.
   - Each chain keeps the trees that vary most, a tie to the earlier tree, in the donor's order, and the
     chosen trees are installed once more. The snapshot is taken once, before the first round, and any
     failure restores every chain. The mean and the variance forest choose separately; each chain chooses
     among the trees of the donor sample it takes.
   - `altered` is set for fewer trees. The report gains the single-leaf count, and `offered` is the donor's
     count.
2. Bridge and R (about 60 lines): the report's new counts, the phrases "150 trees started as single leaves"
   and "variance trees: 4 of 6 taken", and the two count messages removed.
3. Existing tests. The tracer of slice A's step 4 is run again over the files a count refusal is pinned in:
   retired: ["shape-compatible"](../../inst/tinytest/test-warm-start.R) has 3 hits in 2 files, gone from this
   one in slice C, the variance count 1;
   in tests/cpp 2 `shapeMismatch` and 1 `varianceMismatch` expectations are read. Each becomes a check of
   what is taken.
4. New tests (about 400 lines; Tests).
5. Help, NEWS and notes (about 50 lines): "tree count" leaves the two "must share" sentences; NEWS gains "A
   previous fit with another number of trees is taken in part: into more trees all of it, into fewer the
   trees whose fitted values vary most over the new fit's rows."; TODO drops `warm-start-tree-count`.

Commits: (1) steps 1 to 3; (2) step 4; (3) step 5.

### Tests

tinytest in a new `test-warm-start-kinds.R`, tests/cpp in `test_state.cpp`.

| Test | Pins | Mutant |
|---|---|---|
| more trees equal fit | 50 into 200: TRUE; the starting fit is the donor's; the donor's trees in slots 1 to 50 in its order, then 150 single leaves at 0 | leaves not 0; m18, the single leaves ahead of the donor's trees |
| fewer trees declines | 200 into 50: FALSE, untouched | TRUE for a donor with trees left out |
| most varying, by hand | a hand-built donor of 6 trees into 3, the leaves of unequal row counts: the three of largest variance over the rows, by hand | the first three; m14, the variance over leaves and not rows |
| variance, not size | a tree of two leaves at 5.0 and 5.1 beside trees of small values that vary more: not chosen | m15, the second moment |
| a padded round | 7 into 3, the two most varying trees in the second round and the seventh the third most: those three, each once, the forest equal to the one built by hand | m16, a stale tree left in the last round and chosen twice |
| donor order kept | the chosen trees sit in the donor's order with their values | sorted by variance |
| ties to the earlier | two trees of equal variance, one slot | the later tree |
| weights in the ranking | case weights under which the ranking differs from the unweighted one; a row of weight 0 that would change it left out; a logistic recipient: the choice computed by hand with unit weights, the same on two seeds | unit weights; the working weights |
| ranked as installed | an interaction-limited recipient where the donor's most varying tree is cut back to one leaf: not chosen | ranked as stored |
| blocks | a `blocks` recipient of fewer trees: no installed tree holds a split its slot forbids | a tree installed as its trial slot cut it |
| sides drawn once | tests/cpp: a recipient with missing values, fewer trees: each generator's bytes equal those of one pass of coins over the donor sample's trees, and the installed trees' sides are that pass's | m13, the coins drawn again at the final install |
| each chain its own | two chains given two kept draws: each chooses among its own draw's trees | chain 1's choice for both |
| kept draw of another count | `samples` naming a kept draw of a 200-tree donor: that draw's trees | the recipient's count as the stride |
| variance forest, log scale | hand trees with factors (0.25, 1) and (1, 2.2) on equal halves: the first is chosen | variance of the raw factor |
| each forest its own | a heteroscedastic pair, both counts smaller, where the mean and the variance rankings differ: each forest holds its own choice | m17, the variance forest taking the mean forest's indices |
| variance padding | 4 variance trees into 6: the variance function is the donor's | leaves not 1 |
| the phrases | the two phrases of slice C by their cases | a phrase on the wrong case |

"Ranked as installed" holds by construction except under `blocks`, where a tree's allowed columns depend on
its slot: the tree is ranked as cut back for the slot it is tried in and installed in the slot its rank
gives, cut back for that one. The row "blocks" pins what holds there.

## Slice D: parts and leaf kinds

### Steps

1. Parts one side has (about 70 lines, chain.hpp and sampler.hpp). The DART refusal and the variance
   presence refusal go. A donor's DART state or variance trees with nowhere to go are not assembled. A part
   only the recipient has is started as a new fit's when the donor goes in, by a forced install or an
   unforced one that is taken, and never on a decline or an error:
   - DART, following dec-B324's words: the split weights equal over the columns the forest may use, and the
     concentration at the value the options give, as `DartPrior::initialize` leaves them; and the delay
     counter at `updateDelay`, so that the first sweep's update draws the weights from the installed trees.
     A new fit's counter is 0 and holds the weights equal for the delay; the ruling's second clause sets
     that aside, the delay being there for a forest that knows nothing yet.
   - A variance forest: one root per tree at the starting variance, as `VarianceForest::initialize` leaves
     it, the forest's column restriction set again, and the variance of every row rebuilt.
   - Both sides DART, or both with a variance forest: as today, the donor's.
   The report names each.
2. Leaf kinds and the errors of (d) (about 350 lines, the bridge and chain.hpp).
   - [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads each donor forest's leaf blocks in the
     donor's own shape: no `tree.params` is constant; `leaf.lengthscales` is gp, its row count the block's
     length over the tree count; otherwise linear, its covariate count the block's length over the leaf
     count. The donor's shape goes to the engine beside the state, in place of slice A's flag; the state's
     format does not change.
   - The size checks the bridge makes today against the recipient's shape move to the engine, against the
     donor's own: a slope block is the donor's leaves times its covariates, a gp block the donor's rows a
     tree. A block of another size is the malformed-donor error.
   - The anchor ["m23"](../../benchmarks/R/mutation-battery.R) sits in the reader:
     ["verify-anchors"](../../.github/workflows/lint.yaml) runs on every push, so the reader keeps the
     anchor's two lines as they are or the anchor moves with them, and the check is run before the commit.
   - In the engine: constant into linear gives every leaf slopes of 0; linear into another number of
     covariates drops the donor's slopes first, then keeps the intercepts with slopes of 0, and makes no
     restatement, the donor's standardization being of other columns; linear into constant stays slice A's.
     The second is not clean.
   - Errors in both forms: "warm-start donor has gp leaves and this sampler constant ones" for each pair of
     kinds with gp in it, the pair slice A declined included; "warm-start donor's gp leaves were fit to 300
     rows and this sampler holds 200"; "warm-start donor was fit to 4 predictors and this sampler holds 3";
     and the factor errors naming the column by position, and by name where the bridge holds one.
3. Existing tests. The tracer is run again over the files these refusals are pinned in:
   [test-sampler-errors.R](../../inst/tinytest/test-sampler-errors.R),
   [test-warm-start.R](../../inst/tinytest/test-warm-start.R),
   [test-sampler-boundary-refusals.R](../../inst/tinytest/test-sampler-boundary-refusals.R) and
   [test-composition-sequences.R](../../inst/tinytest/test-composition-sequences.R). The DART refusal has 1
   hit and the leaf parameters 1; in tests/cpp the remaining `shapeMismatch` and `dartMismatch` expectations
   are read. Each becomes a check of what is taken, or of the new text.
4. New tests (about 350 lines; Tests).
5. Help, NEWS and notes (about 110 lines). The rest of the two "must share" sentences and "is refused rather
   than silently reshaped" give way to the rule's last rows; NEWS gains "A previous fit with constant where
   the new fit has linear leaves or the reverse, or with DART or a variance forest on one side only, is
   taken in part."; the design note's Status; TODO drops `warm-start-salvage`.

Commits: (1) step 1 with its rewritten tests; (2) steps 2 and 3; (3) step 4; (4) step 5.

### Tests

tinytest in `test-warm-start-kinds.R`, tests/cpp in `test_state.cpp`.

| Test | Pins | Mutant |
|---|---|---|
| plain into DART | taken on a sampler that has run: weights equal over the allowed columns, the concentration its creation value, and after one sweep the weights no longer equal, with `update.delay` above 0 | the error kept; the counter at 0 |
| DART into plain | taken; no weight installed | donor weights installed |
| both DART | the donor's weights, concentration and counter, as today | m20, the reset where both sides have the part |
| variance on one side | recipient only, on a sampler that has run: every variance tree one leaf and, in tests/cpp, every row's variance the starting variance; donor only: taken, the donor's variance trees nowhere | the forest left as it was; m21, the variance of the rows not rebuilt |
| both with a variance forest | the donor's variance trees, as today | m20 |
| a decline leaves the parts | a plain donor that needs a change, into a DART recipient with a variance forest that has run: FALSE; the weights, the concentration, the counter and the variance trees identical | m19, the reset ahead of the decline |
| parts are not judged | a clean plain donor into a DART recipient: TRUE | FALSE for a part on one side |
| constant into linear | TRUE; the starting fit is the donor's; slopes 0 | slopes not 0 |
| linear, other covariates | one covariate into two, the donor under another standardization: FALSE; forced, each intercept the donor's stored one exactly and slopes 0 | the error kept; m22, the restatement run before the slopes are dropped or with the recipient's covariate count |
| linear into constant | still FALSE; forced, each leaf the donor's stored intercept exactly | m22, the intercept restated |
| gp is an error | gp against constant and against linear, each way, and gp on another number of rows: error in both forms, by name, the sampler untouched | taken under force; slice A's decline kept |
| block sizes | tests/cpp: a slope block and a gp block of another size than the donor's own shape gives: the malformed-donor error | a block read past its end |
| predictors and factors | the texts of (d), sampler untouched | a write before the error |
| the phrases | each phrase of slice D by its case | a phrase on the wrong case |

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
    leaves under another standardization, a recipient with missing values, a kept draw by `samples`, two
    chains, linear into constant and gp into constant; and for the 8 fixtures of scratch/wsc/p6.R and a gp
    pair, where the recorded points differ and every threshold is on the grid. Any difference is a stop. It
    is the only oracle for a same-grid donor sent through the across-grids routine.
  - `stateIsValid` is unchanged: tests/cpp's state tests, the 52 tinytest files that call `setState`, and a
    `setState` battery (clean, declined, forced repair) in the same cross-build script.
  - Slices B to D: the same script on the slice before's library and the slice's own. In B the fixtures
    with a monotone recipient holding a tree out of order are expected to differ and every other to be
    identical; a `setState`, `setData` and `setCutPoints` battery on monotone samplers is identical.
- The shifting sequence of slice A has two oracles in the suite: "equals the two calls" above, with its next
  draws, and the equality with a forced `setState` that
  ["equals the two calls"](../../inst/tinytest/test-state-grid-install.R) already holds to the by-hand
  sequence. Slice B's has the hand trees and "in order after".
- Mutants: each row of a slice's Tests table names its own; the reviewer builds each on its own library with
  `--preclean` and names the failing check. Not killable by any test in the suite: a same-grid donor sent
  through the across-grids routine (the cross-build script); the restore on a throw in `installAcrossGrids`
  and the order of the undo across grids, which no test can make throw (TODO `warm-start-undo-test`).
- Lint, air, the codoc, win-drift and doc-freshness checks; `R CMD check --as-cran` from a clean tarball;
  the NEWS parse; in slice D `Rscript benchmarks/R/mutation-battery.R verify-anchors`.
- No bench-sampler run is owed: no sweep path is edited. `stateIsValid`'s loops move without change.

## Budget and stops

| | A | B | C | D |
|---|---|---|---|---|
| engine | about 430 added, 290 removed | about 130 | about 260 | about 220 |
| bridge and R | about 170 | about 15 | about 60 | about 200 |
| tests | about 950 | about 260 | about 450 | about 420 |
| help, NEWS, notes | about 180 | about 50 | about 50 | about 110 |
| planned | about 1730 | about 455 | about 820 | about 950 |
| stop, engine lines added | 600 | 220 | 400 | 350 |
| stop, all | 2700 | 800 | 1350 | 1550 |

A's engine figure is the critique's itemization, 340 to 380 lines, and about 60 for `altered`, the flags
that are not clean and the added counts. Plan estimates have run 1.5 to 2 times low in tests (the last state
slice landed at 2400 lines against 1035), and the stops follow. Stop and report, building nothing further:

- if a bitwise compare prints a "max |z|" line or a snapshot file moves;
- if the cross-build script differs on any same-grid fixture or any of the 8 of p6.R, or in slice B on a
  fixture that holds no monotone tree out of order;
- if a forced install across grids and the by-hand twin disagree on any tree, leaf value or next draw;
- if the tracer shows a call returning FALSE in a file outside the nine named in slice A's step 4, or in
  slices C and D outside the files their step 3 names (list it, do not edit it);
- if a test outside those named fails because it expected a refusal that is now a repair, or because a
  forced warm start across grids gives other trees than today's: list it with both values and do not
  record the new one;
- at a slice's stop in the table. Slice A then lands without its commit 2, the report and the line, which
  follows as a slice of its own: it reads the installed trees and changes no draw, so it takes the
  R-surface tier's review. `bart` forces from commit 1 either way.

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
  bairrtt calls `setPredictor` only. None sees any slice. No reverse dependency can call a warm start that
  0.9-34 did not have.

## Calls made in planning

- A forced warm start across grids keeps the splits the across-grids routine keeps: a leaf that none of the
  recipient's rows reaches under the donor's own thresholds is merged before any split moves, where today's
  arm moves first. That is dec-B295's collapse "as a forced setPredictor and a restore already collapse a
  split with an empty side", and the only order in which dec-B398's mean by rows has rows to weigh.
  Measured above: 173.8 of 196 nodes kept at 30 rows against today's 184.0 (166 against 174 of 182 on the
  critique's fixture), the same within one node from 100 rows on, the starting fit slightly nearer the
  donor's function (0.675 against 0.694). The alternative, a warm start step that keeps unreached leaves
  through the first pass, is about 90 lines and a second rule, and a leaf with no row has no value to
  recover from its fits.
- The monotone collapse tries a split's sides before the split. Step 4 alone is about 25 lines less and
  keeps 22 of the 72 nodes where the sides first keep 32, as many as the least cut on every tree measured. Not taken:
  the search for the least cut, unbounded in the tree's depth; and moving the values into order with every
  split kept, an isotonic fit over a partial order, which is a third rule for a leaf's value where dec-B398
  names the merge.
- A linear leaf whose covariate holds one value in both fits, a different one in each, is not clean: the
  coefficients stay as built (dec-B329), so the start is the donor's fit to its own rows and not its
  function at the recipient's value (9.1 apart against an sd of 2.55 in the run). Forced it installs as
  today.
- An ordered factor column of the donor against a numeric one here follows the routine and is no error; the
  other direction is the error, as for `setState`.
- gp fits start at 0 whenever the donor's recorded points differ from the recipient's, as today, even on
  the same rows, and such a donor is not clean. The fits are values by row, and equal recorded points are
  the one sign the engine holds that the rows are the donor's.
- A merged leaf weighs rows by the working weights, the one routine `setState` merges by; the ranking of
  slice C weighs by case weights.
- Under `blocks`, where a tree's allowed columns depend on its slot, a donor of more trees is ranked as cut
  back for the slot it is tried in and installed in the slot its rank gives. The two can differ for a tree
  near the cut, and the installed forest is valid either way. The alternative refuses another tree count
  under `blocks`. Superseded in slice C's fix round: each tree is tried in every block, which chooses its own
  ([Slice C as built](#slice-c-as-built)).
- Four slices. The monotone collapse is its own because slice A with it plans about 560 engine lines and
  2200 in all, which at the overrun such slices have run is past 3000, more than one review has held. The
  tree count is cut from the kinds and parts for the same reason: together they are three features and 18
  mutant builds.
- A forced install returns NULL invisibly whether or not anything changed: choice (b) read through
  install-surface.md, as asked.
- The verdict is `stateIsValid`'s staging, moved into helpers, and not a new routine. The alternative, a
  second staging for donors, is what the state slices removed for states.
- The cross-grid arm is the across-grids routine handed today's same-grid installs; the remap arm and the
  four scratch passes are deleted, about 290 lines.
- The predictor matrix is a required argument of the C++ `installForests`: a default of null would be read
  on a dense column. The 44 calls in tests/cpp gain it.
- `force` defaults to true in C++, as `setState`'s does, and the bridge states it on every call.
- The counts are made by comparing flat trees after the install with the trees the install read, so no
  install routine is instrumented; `altered` comes from the routines, and tests/cpp holds the two to each
  other.
- The line's words beyond the design note's: "moved to this fit's cut points", "with a side for missing
  values dropped", "gp fits started at 0", "gp fits left behind", "slopes kept across a covariate of one
  value" and the variance clause. They can change after release.
- `bart` calls the internal installer and not the method, since the method returns no report.
- After a forced install the `state` field is the recipient's stored state, as today; a donor's state is not
  the recipient's, so `setState`'s rule of keeping the state as given has nothing to keep.
- Fewer trees are ranked by trial installs in the live slots and not by a scratch staging, at the cost of
  one install per round: the tree ranked is the tree installed, except under `blocks` as above.
- The donor's leaf kind travels beside the state in memory; the stored format is untouched.
- The texts of (d) come with slice D, where the bridge's reader is rewritten once.

## Revision from critique

The blind critique (2026-10-11, of this plan as first written; its probes in scratch/wsc) and what the plan now says.
Slice B of the critique's text is slices C and D here.

- BLOCKING 1, the leaf-kind row: right, a recipient of constant leaves reads no `tree.params`. In A a linear
  or gp donor into constant leaves is not clean; the pairs that are errors today stay errors; "errors stay"
  lists those pairs only.
- BLOCKING 2, verdict against forced change: right. `installForests` reports `altered`; the constraint is
  "FALSE exactly when altered"; [Altered](#altered) says what counts; tests/cpp asserts it and R compares
  trees only where nothing is drawn or restated.
- 3, monotone: right, dec-B355 puts the warm start under dec-B295. Slice B collapses; the reset stays for
  the other calls; "reset to 0" lives in slice A's line only.
- 4, ordered factor: it stays an error, stated as the one error A adds; the Constraints no longer say A
  adds none. The reverse case has its row: moved by the routine, not an error.
- 5, report and line: counts for sides dropped, gp fits cold and variance trees; one count per tree; the
  test builds its donors by hand; "cut back to fit" and ", over K chains" are the design note's words.
- 6, "equals the two calls": pinned at 30 rows, where the arms differ on sgi, with linear leaves and the
  next draws.
- 7, "same grid exact": dropped as a row; the cross-build script is named as the mutant's only oracle.
- 8, restore order: step 2 says the restore runs after the grid is back.
- 9, fewer trees: the sides are drawn once over the whole donor sample, and where the generator stands is
  stated; the `blocks` exception is stated where "by construction" was claimed; an exact multiple and a
  padded round are both tested.
- 10, parts: a part is started as a new fit's only when the donor goes in; a decline is tested to leave
  it; the delay counter goes to `updateDelay`, by dec-B324's second clause.
- 11, restatement: wrong as put. With one value on one side only the start is the donor's function at the
  recipient's rows to 4e-15 (scratch/wsc/f-linear-one-value.R); what moves is the donor's fit to its own
  rows, as for any donor on other rows. The case that does differ is one value on both sides, another on
  each, which has its row and is not clean.
- 12, existing tests: the tracer's rule is per call; a stop covers a test that fails on shifted trees;
  slices C and D run the tracer over their four files.
- MINOR: the rng line no longer says today's code runs for a donor on other recorded points, and those 8
  fixtures are a stop; "declined leaves all" makes the `state` field stale first; "gp across cut points"
  uses a grid that holds every threshold; the anchor m23 is in slice D's step 2 and Verification;
  `columnMaskMismatchMessage` goes in A; the block size checks move to the engine in D; the bart.Rd text
  names moved and reset; A's engine budget is the critique's; the split is made from the start; the cites
  of the deleted routines are counted as the critique counts them, 13 in 10 documents; `setState` is in 52
  files.
- Missing mutants 1 to 12 are in slice A's table under "m" and their number, 13 to 18 in slice C's, 19 to
  22 in slice D's. None is unkillable.
- Settled without a ruling. The monotone reset and the ordered factor: 3 and 4 above. The reset of a run
  sampler's own part: 10. What the line counts and its words: 5. The restatement: 11. gp fits cold on the
  same rows, and the nodes lost at 30 rows against today's arm: kept, each with its reason under Calls made
  in planning, the second resting on dec-B295 and dec-B398. The ranking under working weights: changed to
  case weights.

## Not verified

- Nothing is built: every statement about new code is from reading.
- The monotone collapse is prototyped in R on one donor with unit row weights; its equality with the least
  cut is measured on 43 trees of at most 7 prunings each, and is not a property.
- How many of the 67 calls return FALSE under the default: the critique's tracer approximates the verdict;
  none ran on a build that has it.
- 2 of the 49 family pairs were run.
- That the same-grid forced arm is today's bit for bit: the 8 fixtures of p6.R on sgi, and argued from the
  code path otherwise. A gp donor on other recorded points with every threshold on the grid is argued only:
  sgi's same-grid arm copies the fits, so the planned cold start there has no build to compare on.
- That a linear leaf's cached statistics are invalidated on the across-grids arm as the remap arm
  invalidates them; "equals the two calls" on linear leaves with its next draws is the check.
- `installAcrossGrids` and a throw from its install step: read; no test can reach it.
- That the starting variance a run sampler's variance forest is put back to is the one its creation used
  after the response's range has been re-derived: read, not run.
- Whether expectations in the nine named files pin today's remap values.
- The cost of the trial installs, across grids above all; several chains across grids.
- The 44 tests/cpp calls and the family pairs beyond the two.
- The consumers: read, not run.

## Slice A as built

Built 2026-10-11 on 480d7d69 in four commits: the two forms, the report and the line, the tests, the help and
notes. What was measured on the base build (480d7d69, the state grid slice with its fix rounds) and where the
build parts from the steps above.

- Every number of [Measured on this tree](#measured-on-this-tree) that the build's oracles re-ran is the base
  build's: the cross-build script (scratch/wsa/cross.R) is identical on its 30 fixtures, the same-grid ones,
  the 8 of p6.R, the gp pair on a finer grid and a `setState` battery. So a gp donor on other recorded points
  with every threshold on the grid, argued only above, is the base build's remap arm bit for bit: trees, stored
  state and the next 5 draws.
- The tracer on the built verdict: 11 lines in eight of the nine named files return FALSE, the ninth,
  test-interactions.R, expecting the refusal through `bart`, and none elsewhere. One line of
  test-monotone-unforced.R is FALSE on each of its 55 passes, every one across grids, and one of
  test-state-missing-record.R takes both answers.
- The donor's trees are built on its recorded points wherever those are not the sampler's, and then on the
  sampler's unless the install goes across grids. A threshold that is no point of the donor's own grid is so
  refused by the build that names its forest, with today's two messages, and no separate check is made.
- "Equals the two calls" holds for the chains on 30 rows and for the next draws only under equal weights:
  cut points derived from 30 rows carry weights, the forced start keeps them and the twin's `setCutPoints`
  leaves its points unweighted, so the two samplers draw under other tree priors from identical trees. A
  forced `setPredictor` of each column as it stands, deriving its grid by value, puts the twin's weights back
  and moves no tree. The test holds the forced sampler with its weights to that twin, and a pair whose points
  were set to the donor's and back, which leaves them unweighted, to each other: each pair to one state and
  one stream.
- [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp) pinned the remap arm's 3 live leaves of 11; the
  routine leaves 1, every split on the donor's wider column having a side none of the destination's rows
  reaches. It now holds the trees to a forced state install's and a hand tree to its merge by rows.
- The line's unread leaf blocks are told apart by the bridge, slopes from gp fits by the lengthscales a gp
  forest records; the engine takes one flag. The codes of a store are compared before and after in tests/cpp
  only: R has no reader for them.
- Left to the reviewer's build: the mutant that keeps the remap arm, the routine being deleted.
- The variance clause also closes its counts with ", over K chains" where there are several chains, as the
  first clause does: its counts are summed over the chains too.

## Landing note

Slice A landed 2026-10-17 as a5f3c9bc..d03760e3 on f514c1d8, 11 commits: the plan and its revision from the blind
critique (8a8f84ee, 3b3b24d7, a546599e), the two forms (b9103182), the report and the line (7ab9d1a5), the tests
(94c551d8), the help, NEWS and notes (aeed60ce), two pins (7251c20e), and a fix round (175508a6, 1e728842,
d03760e3). The slice was built on 480d7d69 on 2026-10-11 and rebased onto f514c1d8, the records of the
forest-kind landing, in the fix round. One planner, a blind critique, one opus implementer, one opus review that
ran 23 mutants, and one fix round, which changed tests, the help and the plan and no engine line. Slices B, C
and D are not built.

Built as planned: `installTrees` unforced installs a donor every tree of which stands as stored and returns
TRUE, or returns FALSE with the sampler untouched, and forced installs it with the changes and returns NULL
invisibly; FALSE is exactly the forced install's `altered`; the across-grids routine of the state grid slice
replaces the remap arm and its four scratch passes; the install reports per chain what it moved, cut back,
left without a side for missing values or reset, and the leaf values that are not the donor's; `bart` forces
and prints one line. Departures, the calls of dec-A211: the by-hand twin of "equals the two calls" is the
weighted twin of [Slice A as built](#slice-a-as-built), and the variance clause closes with ", over K chains".

The review found no must-fix and no defect in the install. Eight of its 23 mutants survived the suite: six
are should-fix test gaps, and two (W18, W19) remove the undo after a failed install, which only an injected
throw reaches. The gaps: a linear leaf with a covariate at one value on both sides beside one
restated (W06); a donor of two chains with only the second off the grid (W14); a variance tree that splits
on more columns than the mean trees' interaction limit allows, which `setState` shares (W17); a clean
heteroscedastic donor's line (W13); slopes of 0 and a slope past the first leaf at one value on each side
(W04, W05). The fix round pinned the six in test-warm-start-forms.R, each
failing under its own mutant, held the forced sampler to the weighted twin, wrote "center" in the help, said in the help that a
donor whose splits are all the sampler's cut points needs no move, corrected line 119 of this plan (the grid
is put back after a throw), and recorded the variance clause's chain count. The undo is not pinned (TODO
`warm-start-undo-test`), and the restore of a chain no re-quantization touched is left as it is.

The same reviewer re-checked the fix round on the rebased tip 1e728842. It reran seven mutants (W04, W05, W06,
W13, W14, W17 and D1a), and each was killed by its pin. It confirmed that the rebase kept every one of the
slice's 47 files patch for patch and that bartcore's intervening changes are intact, and it upheld leaving the
undo unpinned. It found one untrue line of this plan, in Constraints, that `installAcrossGrids` "gains the
restore on a throw in A", which d03760e3 corrected. Verdict: land.

Gates, on the fix round's tip 1e728842 (d03760e3 adds a sentence of this plan): arm64, tests/cpp 400 ok plain
and under ASAN and UBSAN with no diagnostics; the full tinytest suite 26715 results, 0 failed; the reference
build's three bitwise compares 55 of 55, 15 of 15 and 11 of 11 and the four snapshot files 7, 3, 14 and 3
results, 0 failed; the 42 exact-gate commands at `quick`, all pass; the composition matrix no disagreement;
lintr, air, rc-codoc, win-drift and doc-freshness clean; `R CMD check --as-cran --no-manual` one NOTE (the Date
field). The reviewer's 23 mutants ran on 959de5c9 and were not rerun after the rebase, the engine being
unchanged. The x86-64 leg on 1e728842: tests/cpp 400 ok, plain and under ASAN with no diagnostics; the tinytest
suite 26748 results, 0 failed; the equivalence compare 55 of 55 statistically indistinguishable. The bcf and
multinomial bitwise compares differed in latent_logistic_weighted and k3countsswap, the two scenarios the
control build differs in as well: they compare bitwise against baselines recorded on arm64 for the shipped
build, so an x86 run is not expected to match.

Not run: Windows, the Linux leak sanitizer beyond the x86 leg's ASAN, the R-loaded sanitizer build on the fix
round's tip (the review ran six test files and its count probe under it), any consumer's suite (stan4bart's
only call, `setState(forceUpdate = TRUE)`, was read against the moved loops of `stateIsValid` and is unchanged;
bartCause and treatSens have no warm start), a gaussian-to-binary and the other family pairs beyond the two,
and the cost of an across-grids install over many chains.

Size against the plan: the plan estimated about 1730 lines for slice A and stopped at 2700. The range adds
3505 lines and removes 771 over 47 files: the engine 442 added and 479 removed (plan: about 430 and 290),
the bridge and R 220 (about 170), the tinytest files 1150 and tests/cpp 530 (about 950 in all), man and NEWS
27, docs 1123 (this plan among them). Without the docs it adds 2382, inside the stop; the tests ran past
the plan because the review and its fix round added the pins above, and the engine's deletions past it.

Left: slices B, C and D; the undo's test (`warm-start-undo-test`); `warm-start-by-position`,
`warm-start-linear-covariates` and `warm-start-bart-families`, which stay in TODO with the later slices.

## Slice B as built

Built 2026-10-17 on d03760e3 and rebased onto 92ae2714 in the fix round, in three commits: the collapse, its
tests, the help and notes.

- The merge is two members of [`Tree`](../../src/bartcore/tree.hpp): `subtreeMergedParameter`, the value a
  candidate is judged by, and `mergeSubtreeToLeaf`, which `collapseEmptyNodesBelow` calls, the sums in the
  order they ran. [`collapseMonotoneTree`](../../src/bartcore/model.hpp) tests only the pairs across a split,
  by `MonotoneLeafGeometry` over the would-be leaves. The reset stays
  [`Chain::reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp), unedited, behind
  `Chain::repairInfeasibleMonotoneTree`. `setData` reaches it through `applyNewData`, `setCutPoints` and a
  forced `setPredictor` through `forceRefreshTrees`'s default, `setState` through `rebuildLiveForest`'s.
- The hand trees are in tests/cpp's test_monotone.cpp, beside its tree builder; the verdict-against-install
  check is a test in test_model.cpp, 20 free donors, half on other rows, all with a tree out of order. "The
  others reset" adds no test: test-monotone.R and test-monotone-unforced.R pin the four calls' reset.
- "Across grids" is a hand tree the move puts out of order, not the critique's donor onto other rows: there
  the trees collapsed on the donor's points stay in order through the move, and the mutant leaving the
  refresh on the reset survived it.
- The critique's fixture: 14 trees collapsed, 94 of 134 nodes kept. Against the prototype's rule
  (scratch/wsc/f-monotone-collapse.R, collapseSides) every tree has the same structure and leaves within
  7e-18, for x1 and x2 increasing, x2 alone, all four, and x1 and x4 decreasing with x2 increasing.
- The cross-build script (scratch/wsb/cross.R, slice A's with a monotone battery of `setState`,
  `setPredictor`, `setCutPoints`, `setData` and a copy) is identical on 36 of 38 fixtures; the two monotone
  warm starts differ.
- The review (LAND AFTER FIXES) found the engine correct and four of its mutants surviving, besides one equivalent.
  The fix round holds the across-grids start to its order along x1 and, under case weights, to its merged leaf by
  hand; pins in tests/cpp the empty-leaf merge at a stride of two and the plain mean of rows that weigh 0; and pins
  that test-monotone.R's gap warm start is taken as stored. The 800-line stop was passed by those tests alone,
  under the coordinator's ruling; the engine stays inside its own stop.

## Landing note, slice B

Slice B landed 2026-10-17 as f5c16b4a..ee1d7213 on f5c16b4a, 5 commits: the collapse and the report (8bdb9504), its
tests (cee8031b), the help and notes with "Slice B as built" (750b2152), a fix round (7148b9f7) and three nits
(ee1d7213). The slice was built on d03760e3 and rebased onto 92ae2714 in the fix round and onto f5c16b4a, the records of
the setData landing, at the landing. One opus implementer, one opus review that ran ten mutants, and one fix round,
which changed tests, text and the plan and no engine line. Slices C and D are not built.

Built as planned: a forced `installTrees` collapses a monotone tree out of order where the order breaks, from the
leaves up and by the rows each leaf holds, on the same-grid arm and across grids; `setState`, `setData`,
`setCutPoints` and a forced `setPredictor` keep the reset; the report loses `reset` and such a tree counts as cut
back; the unforced answer is FALSE exactly when the forced install alters the donor. The engine matches the
critique's prototype tree for tree on four direction sets (7e-18). Departures, the calls of dec-A214: the hand
trees are in test_monotone.cpp and not test_state.cpp; "declined equals altered" is a test in test_model.cpp over
20 free donors, all with a tree out of order; "across grids" is a hand tree the move puts out of order; "the others
reset" adds no test, the four calls' reset being pinned by test-monotone.R and test-monotone-unforced.R.

The review found the engine correct. It proved the collapse always leaves the tree in order, and a stress probe on
36,000 random trees found 0 failures against the predicate, a point oracle and its own reimplementation. It found the
reset callers unchanged (the cross-build script identical on 36 of 38 fixtures; mutants R2 and R3 killed), the report's
readers updated, and consumers reading no counts. R-loaded ASAN over the three test files ran 738 tests with 0
diagnostics. Of its ten mutants R1, R4, R5 and R7 survived the suite, R9 survived as equivalent, and the rest were
killed. One must-fix: the across-grids test in test-warm-start-forms.R never checked the order (R1). Should-fix:
the gap warm start of test-monotone.R was described as merging though the unforced answer is TRUE and nothing merges
(pin it as taken as stored); this plan's account of `applyNewData` (`setData` only, `setCutPoints` and a forced
`setPredictor` going through `forceRefreshTrees`'s default); a weighted across-grids check (R7); the linear leaf's
merged block at a stride of two (R4, which predates the slice) and the plain mean of rows that weigh 0 (R5), both in
tests/cpp; and a line of the design note that still named reset as a repair of a warm start. Verdict: land after
fixes.

The fix round holds the across-grids start to its order along x1 and, under case weights, to its merged leaf by hand,
pins the empty-leaf merge at a stride of two and the plain mean of rows that weigh 0, pins that the gap warm start is
taken as stored, and corrects the text; each new test fails under its mutant. The same reviewer re-checked the round:
R1, R5 and R7 rerun and killed, the rebase lost nothing. Verdict: land. It left three nits, fixed in ee1d7213: a
comment of the monotone tests and a line of TODO the rebase left short, and a sentence of this plan saying that the
fifth surviving mutant was equivalent.

Gates, on the fix round's tip: arm64, tests/cpp 401 ok plain and under ASAN and UBSAN with no diagnostics; the suite
27084 results, 0 failed; the reference build's bitwise compares 55 of 55 against equivalence-32418b6e, 15 of 15
(bcf) and 11 of 11 (multinomial), the snapshot files 0 failed; the 42 exact-gate commands at `quick`, all pass;
`R CMD check` one NOTE (the Date field). On the tip rebased onto f5c16b4a: the suite 27154 results, 0 failed, and the
lint chain clean. The cross-build against slice A's tip: 36 of 38 fixtures identical, the two monotone warm starts
moving as planned. The x86-64 leg on the fix round's tip, whose engine is the landed one: tests/cpp 401 ok, plain and
under ASAN with no diagnostics; the tinytest suite 27117 results, 0 failed; the equivalence compare 55 of 55
statistically indistinguishable. The bcf and multinomial compares differed bitwise in latent_logistic_weighted and
k3countsswap, as the build before does on x86: they compare bitwise against baselines recorded on arm64.

Not run: Windows, and any consumer's suite (the review read the report's readers and found none that reads the
counts).

Size against the plan: the plan estimated about 455 lines and stopped at 800. The range adds 876 and removes 146 over
19 files (871 added at 7148b9f7, before the nits): the engine in src/bartcore 209 added, against about 130 and a stop
of 220, the bridge 15; tinytest 190 and tests/cpp 399 (plan: about 260 in all); man, R and docs 62. The builder stood
at 795 before the review; the fix round's tests carried the range past the stop, and landing went ahead on the
orchestrator's call (dec-A214). The engine stayed inside its own stop.

Left: slices C and D; `warm-start-undo-test`; `warm-start-by-position`, `warm-start-linear-covariates` and
`warm-start-bart-families`, which stay in TODO with the later slices.

## Slice C as built

Built 2026-10-17 on ee1d7213 and rebased onto 6fe5237a, the records of slice B's landing and the setData follow-up,
in four commits: the install, its tests, the help and notes, and a pin of the ranking across cut points. The
fix round, in four more, was rebased onto 738b1015, the records of forest-multiplier push 1.

- The donor forest's own count is the stride of a kept draw, and a donor forest of no trees is the shape
  refusal. A forest given fewer trees than it holds is filled on the assembled scratch after the sides are
  drawn ([`Chain::appendSingleLeaves`](../../src/bartcore/chain.hpp),
  [`Chain::appendSingleVarianceLeaves`](../../src/bartcore/chain.hpp)), so the verdict and the install read it
  as any donor. The scratch builds stage donor tree t in slot t modulo the count, the slot it is tried in
  ([`Chain::stageLiveTree`](../../src/bartcore/chain.hpp)); `stateIsValid` is unchanged, its t never reaching
  the count.
- The rounds are [`Sampler::chooseDonorTrees`](../../src/bartcore/sampler.hpp): one install of every chain a
  round, on the arm the donor takes, the ranking read after it, after the move across cut points included. A
  forest with no round left tries its first trees again and is not read. The ranking
  ([`Chain::liveTreeVariation`](../../src/bartcore/chain.hpp),
  [`Chain::liveVarianceTreeVariation`](../../src/bartcore/chain.hpp)) reads a new virtual,
  [`ResponseModel::caseWeight`](../../src/bartcore/model.hpp): the weight a row carries as given, 0 for a row
  out of the likelihood by weight, mask or membership. The gaussian and Student-t families return their
  weights under their mask, logistic its counts, the others their mask.
- The line, a call the plan leaves open: a padded variance forest is said in its clause, "variance trees: 4 of
  4 taken, 2 started as single leaves", and a variance forest that left trees out prints its clause with no
  change, "variance trees: 4 of 6 taken". The mean forest's "150 trees started as single leaves" follows the
  mean forest's counts. A donor of more trees with nothing else changed prints "50 of 200 donor trees taken"
  alone, its first clause saying why it was declined.
- Step 3. The tracer, now asking the built verdict on a copy before each call (scratch/wsc/tracer.R), over
  test-warm-start.R, test-sampler-boundary-refusals.R and test-sampler-errors.R: FALSE only at two `bart`
  calls of test-warm-start.R, which force; the 5-tree fit is now taken, and the variance count of
  test-sampler-boundary-refusals.R, 4 trees into 5, is TRUE. Both became checks of what is taken. The tests/cpp
  expectations step 3 names pin other refusals, a record of missable columns of another length and a saved
  variance stride, and stand: tests/cpp passed unedited. The "errors stay" row of test-warm-start-forms.R
  listed "tree count", written in slice A with the error this slice lifts (Constraints); it is dropped, the
  one edit outside the files step 3 names.
- Tests beyond the table: "across cut points", a donor onto other rows into fewer trees, held to the most
  varying of the trees a recipient of as many holds after the move; tests/cpp's verdict against the install
  over 24 seeded pairs of other counts. "Ranked as installed" cuts its tree back to the root, the interaction
  limit allowing a root split. Under `blocks` the fixture's four trees kept each landed in the other block's
  slots and were cut back to single leaves: the plan's rule, ranked for the trial slot and installed for the
  final one, at its worst. The fix round replaced the rule (below).
- The cross-build script of slice B (scratch/wsc/cross.R) on slice B's tip and this slice's: identical on all 38
  fixtures, none holding another tree count. Every mutant of the table and five of the builder's own, each alone,
  failed a test of test-warm-start-kinds.R or tests/cpp; the working weights (C11) fail on one of the six
  logistic seeds, so two seeds, as the table has it, could miss them (the fix round's kill takes both).
  "Ranked as stored" is built as a ranking of
  each flat tree routed with no constraint, and "a tree installed as its trial slot cut it" as the final install
  under the trial slot's mask, the containment backstop making the plain reading unreachable.

The review (LAND AFTER FIXES) found the engine right but under `blocks`, where a donor tree tried in a slot of
another block is cut to a leaf and never chosen: a donor of 100 and 100 trees into 25 and 25 held 25 of its 50
most varying trees in that band and kept none. It found nine of its eleven mutants surviving, and wording. The fix
round:

- Blocks, the orchestrator's call, which replaces the plan's "a tree near the cut" under Calls made in planning:
  no donor tree is ruled out by the slot it is tried in. Every donor tree is tried in the slots of every
  `blocks` row, a row's count at a time; each row in turn keeps its count of the trees that vary most as
  installed there among those no earlier row kept, a tie to the earlier tree, in the donor's order within the
  row. Without blocks nothing moves: the rounds and the choice are the old ones. The 100 and 100 donor into 25
  and 25 now keeps the 25 most varying of each half's trees by hand, 25 of them from the band the old rule never
  chose, and into 40 and 10 the 40 and the 10. Cost at 5000 rows and 10 columns, 200 trees into 25 and 25: 23 ms
  against 15 ms for 50 trees with no blocks on the same grid, 31 ms against 20 ms across grids; into 40 and 10,
  20 rounds, 43 ms (scratch/wsc/fix/time-blocks.R).
- A tree cut back to one leaf varied by about 1e-33, rounding in the variance, so a row that keeps more trees
  than vary there took its leaves in rounding order and not the earlier first. The variance is taken from the
  first row's value, and such trees tie at exactly 0.
- A live mean forest whose count does not divide the saved buffer, which a hand-built state can hold, is the
  shape refusal and no longer a slice across two kept draws, as the variance forest already refused it.
- The single leaves are said before the changes in the mean clause, as in the variance clause.
- Kills added to test-warm-start-kinds.R: an active-row mask, gaussian and probit; Student-t weights; logistic
  counts; a variance forest that needs rounds the mean forest does not; the variance ranking under case weights;
  linear and gp single leaves at 0; a pooled factor's masks with the trees chosen; the uneven kept draw. The
  review's deterministic kill for C11, a heteroscedastic recipient, does not kill it: a variance forest's
  precisions are not the response's working weights, so C11 reads the case weights there. It stays as a pin, and
  C11's kill is a logistic sampler whose response is 1 above x1's 50th point, its latent draws there a third to
  two thirds of the others on five seeds, which fails C11 on both seeds it runs.
- Wording: a smaller donor is taken whole and a larger one gives some trees, in NEWS, `bart.Rd`,
  `dbartsSampler-class.Rd` and install-surface.md; NEWS says it after the posterior sentence; "weighted"; the
  `installTrees` docstring; the design note's Status.

Each of the review's mutants K1 to K8 and K11 and the builder's C11, the pre-fix rule under blocks (BK1), a row
reading another row's slots (BK2), rows sharing a tree (BK3), the stride guard off (G1) and the clause order (V1),
each alone on an export, fails test-warm-start-kinds.R; K11 by an install error. C6 and C10, whose code the round
moved, fail it again. Size: the engine adds 417 lines against the stop of 400 and the range 1426 against 1350,
the fix round's tests and the blocks rule (133 added, 68 removed in src) carrying it past.
