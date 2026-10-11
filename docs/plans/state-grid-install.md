# state-grid-install: a state goes onto the sampler's cut grid and does not replace it

Status: FOR DECISION 2026-10-11. Nothing is built. Three calls are open ([Open calls](#open-calls)); the
first decides whether this is built at all, because a premise of the TODO entry failed its run: released
0.9-34 stores the cut grid in its state and installs it, as the branch does today.

agent: opus implementer, one; one opus reviewer who runs the mutants under [Tests](#tests); a blind critique
of this plan first. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)): the diff is in the
engine's install path.
rng: POSTERIOR-CHANGING on one call sequence: a state installed into a live sampler whose cut grid differs
from the grid the state records. Afterwards the sampler holds its own grid where today it holds the
state's, and the set of split points is part of the tree prior. NEUTRAL, bit for bit, for every install
onto the grid the state was stored under, for every copy and reload, and for every fit that never calls
`setState`. No recorded scenario, snapshot file or exact gate reaches the changed path
([Verification](#verification)).
window: before the merge to main; engine slices stay serial. After response-scale-rows (landed on this
base) and before the record slices of TODO `state-frame-prior`. Shared files: [Interactions](#interactions).
budget: planned ~1065 lines (engine ~245, bridge and R ~60, tests ~650, help, NEWS and notes ~110); upper
figure 1700. Stops: [Budget and stops](#budget-and-stops).

## Goal

`setState` on a live sampler leaves the sampler's cut grid and the grid's weights alone and puts the state's
trees onto that grid. A state whose every split threshold is a point of the sampler's grid goes in as
stored and the call returns TRUE. One with a threshold the grid does not hold is declined unforced (FALSE,
the sampler as it was) and, forced, has that split moved to the nearest point its tree still allows, or
merged where none is left (NULL). The state keeps recording the grid it was stored under; a copy and a
reload, which re-create the engine and have no other record of a grid, still take it from there. Closes
TODO `state-grid-replaces-sampler` and the grid half of TODO `state-frame-prior`.

## Context

Rulings, in the register's words unless quoted as the maintainer's:

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
- dec-B300: a state whose recorded grid repeats a point is refused, forced or not. Stays.
- dec-B311 and dec-B312 give a replaced grid's splits "position" or "value" at `setPredictor` and
  `setCutPoints`. They name those two calls, not `setState` ([Open calls](#open-calls), 3).
- dec-B418, the maintainer: "Does restoring across a re-anchor even make sense? The mental model of the
  sampler is a linear one through time, with each sample being a snapshot." The help already tells a
  caller to put a re-anchored response back before `setState`; this plan says the same of a grid.
- dec-B230: a state from a sampler on other data is "neither supported nor prevented".

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
- The warm start already puts another grid's trees onto the sampler's own:
  [`Sampler::installForests`](../../src/bartcore/sampler.hpp) builds each donor tree against the donor's
  grid ([`ScopedCutGrid`](../../src/bartcore/data.hpp)) and moves each split by value
  ([`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp),
  [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp),
  [`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp)): the nearest point inside the interval the
  splits above leave, an exact match kept, a tie going to the upper point, a subtree with no point left
  merged.
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
| 199 points holding the stored 100, restore | 100 again | 100 again |
| column 1 replaced and its grid refreshed, restore without putting the column back | the stored grid back (0.023 to 0.98, from 0.030 to 2.93) | the same |
| kept draws: `predict` before and after 3 points, no restore | unchanged, with a warning | unchanged |
| kept draws: restore after 3 points | 100 again; predictions unchanged | the same |
| `copy()` and reload, state stored after the change | 100 values on column 1: the 3 set, then 97 read past the stored vector (`numCuts` is not rewritten by a store) | 3 and 100 |
| `copy()` and reload, state stored before the change | 100 and 100, the stored state's | the same |
| two chains | 100 again | 100 again |
| a sampler given 3 points takes another sampler's state stored on 100 | 100 | 100 |

Three more facts of 0.9-34. Its `restoreState` reads the trees by index and then installs the state's grid
with the engine's own `setCutPoints`, on every restore (source read on main, the 0.9-34 tree). Its
`storeState` writes into the R object the sampler already holds, so a state kept by plain assignment is
overwritten by the next store: with `updateState = TRUE` on the changing call the kept state held 3 points
and a restore brought nothing back; kept by assignment with `updateState = FALSE`, or copied, it held 100
and the restore brought them back (p5). Its `setState` returns NULL and leaves the `state` field NULL.

So the TODO entry's "against ... 0.9-34, whose state carried no grid" is wrong: the branch does what the
release did. The rule below is a change from the release, made for dec-B231, not the repair of a
regression. That is open call 1.

### Measured: what the rule would do

On d5b637ee, with `installTrees`, which makes the move planned here (p4):

- A sampler given 3 points takes a donor's trees stored on 100: it keeps 3 points; the donor's 19 splits
  on column 1 become 18, on 0.25, 0.5 and 0.75 (5, 6 and 7 of them), and one is merged with what hangs
  under it. A split at 0.333 goes to 0.25 and one at 0.954 to 0.75; the largest move is 0.208.
- A sampler given 199 points that include the donor's 100 takes the donor's trees with every split as
  stored.
- Undo by hand: a column replaced with `updateCutPoints = "position"`, replaced again with the old values
  and the same word, then `setState(stored)`: the grid is the stored one bit for bit, TRUE, and the next
  five draws are identical to those of a twin that only restored its own state.
- The same with the old values put back and no refresh: none of the state's 19 splits on that column is a
  point of the grid derived from the new values. Today `setState` returns TRUE and brings the old grid.
- A weighted grid (a column of four distinct values: points 0.5, 1.5, 3.5 with weights) is put back by
  `setPredictor(values, column, updateCutPoints = "value")`, points and weights. `setCutPoints` with the
  old points puts back the points and not the weights; it takes no weights (TODO `cut-point-weights`).
- Two samplers made on the same data hold identical grids and weights, which is stan4bart's restore.
- A state whose column is numeric goes into a sampler whose column is a factor: "state is not consistent
  with this sampler", forced or not.

## The rule

On a live sampler (an engine that existed before the call):

| The state against the sampler | Unforced | Forced |
|---|---|---|
| its recorded grid equals the sampler's on every column | as today, bit for bit | as today |
| grids differ; every live split's threshold is a point of the sampler's grid | TRUE; trees as stored; the sampler's grid and weights untouched | NULL; the same |
| a live split's threshold is not a point of the sampler's grid | FALSE; nothing touched | NULL; the split moves to the nearest point inside the interval the splits above it leave (a tie to the upper point); with no point left it and everything under it become one leaf at the mean of the leaves replaced, weighted by rows |
| a column numeric on one side and a factor on the other; a recorded grid that repeats a point or is malformed | error, as today | error, as today |

- The weights of a grid ([`cutMass`](../../src/bartcore/data.hpp)) are part of the grid: no live install
  writes them, whatever the state records.
- Splits on a factor column carry no threshold and are not touched. An ordered factor's grid is its level
  midpoints and follows the same rule as a numeric column.
- A rule's missing-value direction rides the rule through a move. The draw for a column that can be
  missing now and could not when the state was stored
  ([`Chain::completeMissingDirections`](../../src/bartcore/chain.hpp)) is made first, as today.
- After a move the trees are judged as today: a leaf no row reaches, a forbidden split, a monotone tree
  out of order are each not clean, and repaired when forced.
- The variance forest's live trees follow the mean forests' rule.
- Kept draws go in as stored, never moved: they predict what they predicted. A kept draw can therefore
  hold a threshold the grid no longer has, as it can today after `setCutPoints`.
- Several chains share one grid; one chain's tree that needs a move declines the whole state.
- `installTrees` and `bart(warm.start = )` are not changed.

At a re-creation (the engine is created inside the call: a reload's first use, `copy`, `setState` on a
dead pointer) the new engine takes the grid and weights the state records and then installs the state, as
today. A sampler copied or saved after a grid change with no store since is re-created on the stored
state's grid, as today; `storeState()` after the change makes the copy match.

Undoing a grid change is the caller's, before `setState`: replace the column with its old values under the
same `updateCutPoints` word, or `setCutPoints` with the old points. Then the state is on its own grid and
the install is the stored chain.

## Where the line falls

In this slice: what a live `setState` does with the grid, its weights and the splits. Not in it, and left
to the record slices of TODO `state-frame-prior`: the `cutPoints`, `standardization` and `response.range`
records on dbartsData (dec-B231 to dec-B233), `setData`'s three arguments and its carried-record rule, a
creation that uses a carried record, a linear or gp leaf's standardization and slopes at an install, the
manual's table of what puts each derived value back, and the gp `setData` refusal.

The grid half stands alone on a live sampler: the trees carry their thresholds, the warm start's move is
built, and nothing has to be recorded anywhere new. It does not stand alone at a re-creation. There the
state is the only record of a grid that a creation would not derive, so a copy and a reload keep taking
the grid from the state until dbartsData records it; taking that away now would make a copy of a sampler
with a set grid something other than a copy. Open call 2.

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

- The same-grid install is the code that runs today, reached by one comparison of the grids; the three
  bitwise compares and the four snapshot files are the proof.
- The nearest-point rule is the one function the refresh, `setData` and the warm start use
  ([`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp), by value); its results there do not move.
- The verdict and the install judge the same trees by the same routine
  ([`Tree::stageFromFlat`](../../src/bartcore/tree.hpp)), as now, so the unforced answer is FALSE exactly
  when the forced install changes a tree.
- A declined or refused live install writes nothing: no grid, code, chain, generator or kept draw.
- No new argument, no message, no warning (dec-B234: "Return a status, no warning").
- Out of scope: everything under [Where the line falls](#where-the-line-falls); a reader for the grid in
  force; weights on `setCutPoints`; `installTrees`'s `forceUpdate` (TODO `warm-start-salvage`).

## Steps

1. Nearest point by value (~45 lines, [tree.hpp](../../src/bartcore/tree.hpp)). Split the by-value arm of
   [`Tree::mappedSplitIndex`](../../src/bartcore/tree.hpp) into a routine that takes a threshold, a
   column's points and the open interval of positions, and returns the position; `mappedSplitIndex` calls
   it. Verify: tests/cpp unchanged and green.
2. Flat trees onto a grid (~75 lines, tree.hpp). A walk over one flat tree in pre-order that carries each
   column's interval down, as [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp) does over a
   built tree, and rewrites each ordinal node's value to the point step 1 gives; where the interval holds
   no point it takes the nearest point of the whole column, which
   [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) then marks as outside its interval for the
   existing merge. Returns whether any value changed. Leaves, masks, flags and slopes are not touched; a
   malformed tree is left as given for the existing refusal. Verify: step 2's own tests/cpp cases.
3. The install (~125 lines, [sampler.hpp](../../src/bartcore/sampler.hpp),
   [chain.hpp](../../src/bartcore/chain.hpp), [facade.hpp](../../src/bartcore/facade.hpp)).
   [`Sampler::setState`](../../src/bartcore/sampler.hpp) gains `adoptGrid`. True: today's body. False: no
   grid snapshot, install or rollback; where the recorded grid differs from the sampler's, the live mean
   and variance trees of a copy of the chains (the copy the missing-direction draw already makes, or a
   new one) go through step 2; "a value changed" is folded into the verdict beside `needsRepair` and into
   `altered`; the rest is unchanged. The column-kind and grid-validity checks at the top stay. The
   virtual in [`SamplerBase`](../../src/bartcore/facade.hpp) takes the flag: install with `--preclean`.
   Rewrite the comment block above `setState` and the one on `SamplerStateData`.
4. Bridge and R (~60 lines). [`bartcore_setState`](../../src/R_interface_bartcore.cpp) takes the flag as
   a sixth argument (its row in [R_interface.cpp](../../src/R_interface.cpp)); the raw predictors are read
   only under it. [`installStateOnto`](../../R/dbarts.R) passes it: TRUE from `getPointer` and `copy`, and
   from `setState` exactly when it re-created the engine. Method docstrings of `setState`, `copy` and
   `storeState`.
5. Tests ([Tests](#tests); ~650 lines). The cases that pin today's rule are rewritten, not deleted:
   ["a state stored on the longer grid installs"](../../inst/tinytest/test-cut-points-undo.R) and its
   neighbor, ["a state without them installs the grids unweighted"](../../inst/tinytest/test-cut-grid-distinct.R),
   ["the setState cut-point replay"](../../inst/tinytest/test-data-mixed-mutation.R) and its sparse twin
   in [test-mutate-sparse-valued.R](../../inst/tinytest/test-mutate-sparse-valued.R) (both move to a
   reload, which still re-quantizes), and in tests/cpp
   [`testRestoreStatus`](../../tests/cpp/test_state.cpp)'s cross-grid block and
   [`testWeightedGridState`](../../tests/cpp/test_state.cpp).
6. Help, NEWS and notes (~110 lines). [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the
   `state` field, Saving (the list of what does not fit gains the threshold; the repair; "The cut grid ...
   come with the state as they are" becomes the rule and the undo recipe; the re-creation sentence),
   `setCutPoints` and `updateCutPoints`. NEWS: `setState` no longer brings back a cut grid, which 0.9-34
   did, with the two-line undo. [install-surface.md](../design/install-surface.md) gains the section the
   posterior-changing class asks for; one sentence each in
   [state-not-model.md](../design/state-not-model.md), [cut-grid.md](../design/cut-grid.md) and
   [architecture.md](../architecture.md). TODO: drop `state-grid-replaces-sampler`, cut the grid half from
   `state-frame-prior` and correct its sentence on 0.9-34.

## Tests

tinytest in a new `test-state-grid-install.R`, tests/cpp in `test_state.cpp`. Each line: what it pins;
the mutant it must kill.

| Test | Pins | Mutant |
|---|---|---|
| live grid kept | 3 points set, an earlier state forced in: the sampler reports 3 points and its weights as before | `adoptGrid` true on a live call |
| same grid exact | store, 5 sweeps, restore: TRUE, printed trees as stored, next draws equal a twin's | step 2 run on a same-grid state and reported as a change |
| finer grid is clean | 199 points holding the stored 100: TRUE, 199 kept, every split as stored | verdict keyed on the grids differing, not on a moved split |
| off the grid declines | FALSE; grid, trees, `state` field and next draws equal a twin that never called | "changed" not folded into the verdict; a write before the decline |
| nearest inside the interval | hand-built one-tree states: 0.37 on {0.25, 0.5, 0.75} goes to 0.25; midway goes up; a child stays below its parent's new point | nearest replaced by the point below; tie down; interval dropped |
| emptied interval merges | parent on 0.5, left child at 0.3, grid {0.5, 0.75}: one leaf at the row-weighted mean, by hand | plain mean; the child kept outside its interval |
| forced leaves the field | NULL; `state` holds the state as given; installing the field unforced is FALSE | the repaired state cached |
| weights are the sampler's | a state with no `cutMass`, and one with other weights, change no weight | the weight write left on the live path |
| re-creation adopts | stored after 3 points: `copy`, reload and a dead-pointer `setState` report 3 points and the weights; mixed and sparse stores re-quantize | `adoptGrid` false in `getPointer` or `copy` |
| kept draws stay | `keepTrees`: forced across grids, `predict` is bit for bit what it was and kept thresholds are the old ones | step 2 applied to saved trees |
| variance forest | a variance split off the grid alone gives FALSE; forced, it lies on the grid | variance trees skipped |
| factor column beside a changed grid | subset rules identical after a forced install; numeric against factor errors in both forms | a mask read as a value |
| direction rides the move | a moved rule keeps its missing side; a state from before a first missing value draws the sides a same-grid install draws | flags cleared by step 2 |
| second chain | an off-grid split in chain 2 only: FALSE, chain 1 untouched | verdict from the first chain |
| linear leaf | forced across grids: unmerged leaves keep intercept and slopes | slopes shifted by a move |
| undo recipe | refresh away and back, restore: grid identical, TRUE, draws equal a twin that only restored | (pins the help) |

## Verification

- `cd tests/cpp && make && ./test_bartcore`; again under `OPT="-O2 -g -fsanitize=address,undefined"`.
- `R CMD INSTALL --preclean -l <lib> .`, then the full tinytest suite with `R_LIBS=<lib>`; the R-loaded
  sanitizer run over the new test file.
- Reference build: the four `test-reproducibility-*.R` files pass unedited, and
  `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-308c45be.rds --bitwise
  --strict-coverage`, `bcf-equivalence.R compare .../bcf-equivalence-1b7d730c.rds --bitwise` and
  `multinomial-equivalence.R compare .../multinomial-equivalence-80b1c8d4.rds --bitwise` print 55, 15 and
  11 "identical draws" lines and no "max |z|" line. Read: no scenario of the three harnesses calls
  `setState` or `setCutPoints`; the snapshot files call neither.
- The 42 exact gates in quick, as the class requires; expected unchanged. Read: of the scripts under
  benchmarks/R only negbin-mixing, sbc, the frozen-ESS surface and bench-sampler call `setState`, each
  with a state its own sampler, or one made on the same data, stored.
- Lint, air, the codoc, win-drift and doc-freshness checks; `R CMD check --as-cran` from a clean tarball.
- stan4bart's suite against the slice's library, expected unchanged.

## Budget and stops

| Step | Lines |
|---|---|
| 1, 2 tree.hpp | ~120 |
| 3 sampler, chain, facade | ~125 |
| 4 bridge and R | ~60 |
| 5 tests | ~650 |
| 6 help, NEWS, notes, TODO | ~110 |

Stop and report, building nothing further: before open call 1 is ruled; if a bitwise compare prints a
"max |z|" line or a snapshot file moves (the same-grid path changed); if step 2's output is not accepted by
`buildFromFlat` for some grid (the walk and the build disagree on a point); if a tinytest file outside
step 5's list fails because it expected a state's grid to arrive (list it, do not edit it); at 370 engine
lines or 1700 in all.

## Interactions

- response-scale-rows, landed on this base: no shared routine. Its slice C adds a help sentence and a
  test that `setState` does not undo a re-derivation; it shares the Saving section of the Rd and nothing
  else.
- The forest slices (forest-sd-unit, forest-kind-by-class, forest-multiplier): chain.hpp, combiner.hpp,
  the bridge, [R/dbarts.R](../../R/dbarts.R) and the Rd, none in `Sampler::setState` or tree.hpp's flat
  routines. Rebase only; serial as engine slices.
- collapsed-scale-move (branch wt/collapsed-move-plan): chain.hpp and model.hpp, and it names
  `Chain::setState` for the generator a state carries. This slice adds one helper to chain.hpp at most.
- TODO `warm-start-salvage`: `installForests` and the remapped rebuild, which this slice does not edit;
  step 1 touches the function both call and must leave its results as they are.
- The record slices of TODO `state-frame-prior` follow: with the grid on dbartsData a re-creation reads
  it there and `adoptGrid` goes.

## Open calls

1. Build it, now that 0.9-34 is known to do what the branch does?
   Background: a sampler's cut points are the thresholds its trees may split on. A stored state holds the
   trees. In 0.9-34 and on the branch the state also holds the cut points and `setState` puts them back.
   dec-B231 classed cut points as data, and the register concluded a state should not install data; the
   TODO entry then called today's behavior a regression against 0.9-34. The run says it is 0.9-34's
   behavior.
   Example: 100 points per column; store the state; set 3 points on column 1; `setState(stored)`.
   - A, the rule of this plan. Unforced: FALSE; the sampler keeps its 3 points and the trees it had. Forced:
     3 points kept; the state's 19 splits on the column become 18 on 0.25, 0.5 and 0.75, a split at 0.333
     becoming 0.25, and one is merged. To get the stored chain back the caller first sets the 100 points
     again, or refreshes the column from its values, and then restores: TRUE, bit for bit. Cost: a change
     from 0.9-34. Code that undid a grid refresh with `setState` alone now gets FALSE and an untouched
     sampler, and code that never read the value (0.9-34 returned NULL) does not learn it; a weighted grid
     can only be put back by a refresh. ~1065 lines.
   - B, leave it. 100 points again, the 3 gone, TRUE, as in 0.9-34. Cost: `setState` changes a part of
     the sampler the caller did not name, so a sampler given a grid on purpose loses it to any older or
     foreign state (the run: a sampler on 3 points takes another's state and is on 100), and the grid stays
     the one derived value a state installs. No code; the TODO entries and the help are corrected.
   Recommendation: A. It is what dec-B231 and the division of data, model and state ask for, it is the
   rule the response range already follows (a state does not undo a re-anchor; the caller puts the
   response back first), and no consumer package or vignette restores across a grid change. What would
   change it: a known caller that rolls back a grid refresh through `setState` alone.
2. Under A, do a copy and a reload keep taking the grid from the stored state?
   Background: a copy or a reloaded sampler is rebuilt from the R objects, which do not record a grid; a
   fresh engine derives one from the predictors as they are.
   Example: set 3 points, store, `copy()`.
   - A, keep (this plan). The copy has 3 points. Cost: a state still installs a grid on these paths until
     the data object records it; a copy made after a grid change with no store since has the older grid,
     as today.
   - B, build the `cutPoints` record on dbartsData now. The copy has 3 points from the record, also with
     no store since the change. Cost: about 550 more lines in the store's creation path, and it needs
     `setData`'s rule for a carried record settled first (dec-B231's "re-derive if no values are present
     ... if values are supplied, use them"), which is the record slices' design.
   - C, derive afresh. The copy has 100 points and the state's trees moved onto them: not a copy.
   Recommendation: A.
3. Under A, where does a forced install put a split whose threshold the grid lacks: by value only, or by
   the caller's choice of "value" or "position" as at `setCutPoints`?
   Example: the state's split at 0.333, the 33rd of 100 points, into a sampler on {0.25, 0.5, 0.75}.
   - A, by value, no argument (this plan): 0.25, the nearest. Unforced, any threshold not on the grid is
     FALSE.
   - B, an argument `splits = "value"` or `"position"`. By position the split takes the same share of
     the new grid: 0.25 here, and on a grid of 100 points over another range the 33rd of those. Cost: an
     argument on `setState`; by position a restore onto a refreshed grid of the same length changes every
     threshold and still counts as exact, so a rollback that forgot to put the grid back returns TRUE.
   Recommendation: A. 0.9-34 kept every threshold (it brought the grid), the warm start moves by value,
   and dec-B234's entry already says "moved to the nearest current split point". An argument can be
   added later without breaking a call.

## Calls made in planning

- A threshold on the grid is clean and one off it is not: read from dec-B305's rule ("a split off the
  grid") and dec-B234's list, both the entries' words. The alternative, any grid difference being not
  clean, would decline the finer-grid case, where no tree changes.
- The state keeps both grid attributes. The alternative, dropping them and resolving thresholds on the
  sampler's grid alone, needs no record for a move by value but leaves a re-creation with no grid and
  the tests with no reader.
- A dead-pointer `setState` adopts the grid of the state it is given, as today, not the grid of the
  sampler's own cached state. The alternative re-creates from the cached state first, which a sampler
  whose cached state is refused could not get past.
- The forced merge weighs leaves by rows, as `setState`'s other repairs do, not by the last sweep's node
  statistics as the warm start's move does.
- The trees are moved in flat form ahead of the existing staging, not built against the state's grid and
  remapped as the warm start does: one staging routine then serves the verdict and the install, and
  nothing swaps a grid over the store.
