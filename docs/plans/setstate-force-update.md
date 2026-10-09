# setstate-force-update: setState's two forms, the kept store, the spread, and missingness first seen

Status: PLANNED (dec-B305, dec-B310, dec-B318, dec-B320, dec-B321, dec-B322, dec-B378, dec-B384). The
first slice of the install surface; the warm start (warm-start-salvage, warm-start-tree-count) is later.

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants
below, per part.
rng: Part A SHIFTING for a sampler whose column first gets a missing value after creation, or loses its
missing values (each chain's generator draws directions it did not draw before); NEUTRAL, bit for bit, for
every sampler whose columns' missingness never changes after creation. Part B NEUTRAL for a clean
install, a copy and a reload in the sampler's own store size and units; SHIFTING where a drawn k installs
across response units (dec-B384) or a store of another size installs (refused today).
window: before the merge to main; engine slices stay serial.
budget: ~1700 lines in two landings: Part A ~700 (engine ~200, bridge and R ~80, tests ~350, docs ~70),
Part B ~1000 (engine ~260, bridge and R ~90, tests ~520, help, NEWS and docs ~130).

## Goal

`setState(newState, forceUpdate = FALSE)` installs a state that goes in cleanly and returns TRUE, or
leaves the sampler exactly as it was and returns FALSE; forced, it installs, repairing what it must, and
returns NULL invisibly. A conversion of units, latents redrawn under other weights, a kept store of
another size and a missing direction on a column that has been missing are clean. Whether a column can
be missing is a property of the sampler, raised by the first missing value it is given and never
lowered, and that first value draws each existing split's direction at one half. A drawn k installs as
the spread in force in the state. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Inputs, untracked (dec-B389: near-term work): the design scratch/install-surface/design.md (2026-10-07)
and its critique scratch/install-surface-critique/critique.md (BUILD AFTER CORRECTIONS). The design
predates the rulings below; where they differ the ruling wins:

| Design said | Ruled since | This plan |
|---|---|---|
| forced FALSE when repaired | dec-B310: forced forms return NULL invisibly | NULL, invisibly |
| a dropped missing direction needs force | dec-B320: kept, clean; dec-B321: first seen, ever after | Part A |
| a larger kept store needs force; slice 2 | dec-B318: clean, newest kept, TRUE; "the install surface's first slice" | Part B, step B4 |
| (not covered) | dec-B384: k re-expressed against the recipient's k.scale | Part B, step B5 |
| (not covered) | dec-B378: a one-way categorical rule refused on a column never missing | Part A, step A4 |
| (not covered) | dec-B322: predict refuses a missing value only where the column could never be missing | Part A, step A5 |
| grid half (slice 4) | critique: out; TODO state-frame-prior | out: the grid arrives with the state, as built |

Bound by the critique's corrections and kept here: a forbidden column, an interaction limit and a
monotone order out of the cone stay errors in setState, as now (the warm start alone collapses them,
warm-start-salvage); one error text per cause and a new floor on leaf values and sigma are out (TODO
precision-floor-edge-values); an engine `force` argument defaults so tests/cpp call sites stay; the
verdict leaves a dead pointer dead; the help leads with the checked form; the verdict and the install
cannot disagree (one routine, step B2).

Code read on the tip 17df64ff:
- [`Sampler::setState`](../../src/bartcore/sampler.hpp) validates every chain before any is touched,
  except the state's cut grid, installed over the store before the containment and validity checks and
  put back by `restoreGrid` on a refusal (integers and vectors moved back, so exact). Its `altered`
  report is set by the units pass (`valuesMoved` from
  [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp)) and by
  [`rebuildLiveForest`](../../src/bartcore/chain.hpp) and
  [`rebuildVarianceForest`](../../src/bartcore/chain.hpp): a bottom no row reaches, a split outside its
  interval ([`holdsSplitOutsideInterval`](../../src/bartcore/tree.hpp)), or a direction dropped by
  [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp).
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) already builds every live tree on a scratch
  `Tree` but never partitions it; it refuses a saved block whose size is not capacity times tree count.
- [`Chain::setState`](../../src/bartcore/chain.hpp) installs the response's latents before the trees, so
  a merge weighs leaves by the state's working weights; it writes `forest.k = fs.k` where the forest
  draws k.
- `ColumnStore::hasMissing` ([data.hpp](../../src/bartcore/data.hpp)) is rewritten from content on every
  quantize ([`quantizeDenseObserved`](../../src/bartcore/data.hpp),
  [`quantizeDenseCodesInto`](../../src/bartcore/data.hpp),
  [`quantizeCscColumnInto`](../../src/bartcore/data.hpp)); [`setCell`](../../src/bartcore/data.hpp) only
  raises it. A raised flag makes the birth and change moves draw a direction at one half and halves a
  rule's prior on the column, one rule per direction;
  [`dropStaleMissingDirections`](../../src/bartcore/chain.hpp) clears directions on a column whose flag
  is down, from the data-mutation paths (applyNewData, forceRefreshTrees, rebuildFitsFromParameters and
  the variance forest's two) and from `buildFromFlat`.
- R: [`installStateOnto`](../../R/dbarts.R) is the install shared by
  [`dbartsSampler$setState`](../../R/dbarts.R), [`dbartsSampler$getPointer`](../../R/dbarts.R) and
  [`dbartsSampler$copy`](../../R/dbarts.R); the test-data refusal
  [`unroutableTestColumns`](../../R/data.R) reads whether the TRAINING matrix holds a missing value now.

Run on the tip (scratch/isplan/01-behaviours.R and 02-restore-cost.R, library scratch/libs/isplan):
setState returns TRUE invisibly on its own state; a state of 10 kept draws into a store of 4, 4 into 10
and 4 into none are refused ("state is not consistent with this sampler") and none into 4 installs, TRUE;
a state stored while a numeric column held missing values, restored after a forced setPredictor filled
them, gives FALSE; a first missing value by an unforced setPredictor is taken (TRUE) and predict then
answers for it; a state in other response units
gives FALSE; across units a drawn k installs as stored, k.scale 7.83 to 23.48 and the spread 2.47 to
7.40.

## Part A: missingness first seen (dec-B320, dec-B321, dec-B322, dec-B378)

Lands first: once a column's flag stays raised, a direction is never dropped on the sampler's own trees,
and Part B's verdict has one direction case left (B2, case c).

A1. The flag is sticky. Every train-side quantize ORs into `hasMissing[j]` instead of assigning it;
creation sets it from content (as today) or from the record (A3). The rollback snapshots of
`WholeMatrixUpdate` and `SubsetUpdate` ([sampler.hpp](../../src/bartcore/sampler.hpp)) already put the
flag back on a refusal. setData keeps today's rule, the flags reset and taken from the new data's
content before it re-quantizes, and its directions as built (dec-B321 does not cover setData; Open call
4). The test store tracks no flag.

A2. The first missing value draws. Where an accepted predictor change raises a column's flag from 0,
each chain draws, from its own generator, a direction at one half for every split already on that
column, in a fixed order: forests in order, trees in order, nodes in pre-order, then the variance
forest's trees, with `ext_rng_simulateBernoulli(rng, 0.5)` as the birth draw takes it. Sites:
- [`runPredictorTransaction`](../../src/bartcore/sampler.hpp), forced: after `applyForced`, before
  `forceRefreshTrees`.
- unforced: after `snapshotApply` and before [`revalidateAllChains`](../../src/bartcore/sampler.hpp),
  since occupancy and the monotone order are judged under the drawn directions. A rollback puts back the
  bits drawn and each chain's generator (serialized before the draw, as
  [`Chain::setState`](../../src/bartcore/chain.hpp) reads one back), so a refused update is untouched
  bit for bit.
- the per-observation session [`UpdateSessionImpl`](../../src/bartcore/sampler.hpp): in
  `observationWouldRemainValid`, the first missing row raises the flag provisionally and draws, then
  judges the row and the monotone order under the draws (replacing `orderHoldsWithMissing`'s temporary
  raise). The draw stays pending until `commitObservation`; a row not committed, because this session
  or, in [`updatePredictorPerObservationJointly`](../../src/bartcore/facade.hpp), another sampler's
  declined it, has the flag lowered and the bits and generators put back before the next row is judged
  or at `finalize`, so the next missing row draws the same directions.
A column already flagged draws nothing; a pooled column (more than 63 levels) keeps its bit in the pool
words, as built (dec-A165), and draws for its rules alike.

A3. The record survives a re-creation. A copy, a reload and a dead-pointer setState re-create the engine
from data@x, which may hold no missing value in a column the sampler has seen missing; without a record
the copy would hold other flags than its source and draw differently. The record is an attribute,
"missing.seen" (a logical per column, absent when none is raised beyond content), on data@x:
- written R-side after every accepted predictor change (every path of
  [`bartcoreSamplerSetPredictor`](../../R/bartcore.R), the joint row update) from a new bridge reader of
  the engine's flags, and carried wherever R rebuilds data@x on a mutation, beside the design attributes
  it carries today;
- read by the bridge at creation ([`bartcore_create`](../../src/R_interface_bartcore.cpp) and the
  handle path), ORed into the content flags;
- read by [`unroutableTestColumns`](../../R/data.R): a recorded column is routable (A5).
setData replaces data@x and so drops the record (A1). Placement is Open call 2.

A4. No direction is dropped afterwards, and dec-B378. With A1, `dropStaleMissingDirections` clears
nothing on a column that has been missing; its sites stay for setData and for a flat tree from another
sampler. In [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) the gauge on a column whose flag is down
keeps counting the missing position reachable (`missingReaches`), but a categorical rule that, with its
missing bit dropped, sends every reachable level one way is refused as malformed (dec-B378), where today
it is built and merged; on a flagged column it is built, and merged where no row is missing now (not
clean, B2). A numeric or inline categorical rule with a direction on a column whose flag is down loses
it, reported as today through `directionDropped` (Open call 3).

A5. predict (dec-B322). The R refusal reads the record beside the training content, so a column made
missable by setPredictor, then filled, is still routable; a column never missing is refused as built.
The engine's own check (store.hasMissing at the predict entries) follows A1 with no change.

A6. Docs: [mia-missingness.md](../design/mia-missingness.md) takes an amendment under its Status
paragraph and in [Bridge and R surface](../design/mia-missingness.md#bridge-and-r-surface) ("a column
that gains NAs mid-run routes them by ... (left)" and "cleared when a column loses its NAs" become the
first-seen rule).

## Part B: setState's two forms (dec-B305, dec-B310, dec-B318, dec-B384)

B1. The flag. `Sampler::setState` (and the facade virtual on
[`SamplerBase`](../../src/bartcore/facade.hpp)) take, after `adoptCapacity`, `bool force = true`, `bool
kFollowsUnits = false` and `bool* notClean = nullptr`, so the tests/cpp call sites (87 lines matching
`setState(`, chain-level calls among them, ran: grep) keep today's behaviour. `--preclean` on every
install: a virtual changes.

B2. Clean, in code. After every refusal (unchanged) and with only the state's grid written, an unforced
call asks each chain for a verdict and, if any chain is not clean, puts the grid back, sets `*notClean`
and returns false, the chains, store, generators, flags and latents untouched (and, in R, the mirrors,
B6). A chain is not clean when
any live mean tree or variance tree, built from its flat form and partitioned over the sampler's rows:
- (a) has a bottom node no row reaches;
- (b) holds a split outside the interval its ancestors leave;
- (c) carries a missing direction on a non-pooled column this sampler has never seen missing.

Everything else that differs from the state is clean: the units pass (`valuesMoved` stops setting
`altered`), latents redrawn under other weights or censoring (`reapplyWeights`, `reapplySurvivalStatus`,
after the install as today), a store of another size (B4), k re-expressed (B5), DART weights, a fixed
value or a generator of another kind left as the sampler's.

One routine: a new Tree-level staging call (build from flat, partition, report a, b and c) used by the
verdict on a scratch tree with one n-length index buffer, as `stateIsValid` builds today, and by
`rebuildLiveForest` and `rebuildVarianceForest` on the live tree, so the verdict and the install's report
cannot disagree; a tests/cpp check asserts it over the fuzz harness's states (unforced false exactly when
the forced install reports a repair). Forced, the verdict is skipped and the install is today's.

B3. Its cost, and the skip. Measured on the tip (ran: 02-restore-cost.R, one chain, one thread, 200
trees, 10 columns, machine under load 5 to 6): a clean restore of a fitted state costs 74.0 ms against a
41.0 ms sweep at n = 1e5 (1.81 sweeps) and 3.3 against 2.1 ms at n = 5000 (1.53); restoring single-leaf
trees costs 19.2 and 1.0 ms, so the per-tree build and partition is 74 and 70 percent of a restore, 1.34
and 1.07 sweeps. The scratch verdict repeats that work less the fits, so an unforced restore costs up to
about 1.1 to 1.3 sweeps more (bound, not run); the critique estimated 0.8 to 1.4 (its crit-04b).

The package's own loop (store, propose, put the data back, restore) pays it on every rejection, so the
verdict is skipped where it cannot find anything. storeState writes a new optional state attribute,
`rows.digest`: one 64-bit hash of the row count, the column count, each column's grid, flag and training
codes as stored (dense and CSC storage alike), and the structure of the state's live trees (variables,
flags, split values, masks, not leaf values). An unforced setState whose state carries it, equal to the
same hash taken of the sampler now (after the state's grid is in) and of the trees as read, skips B2: a
sampler's own live trees hold none of a, b or c over the rows they were written on (the empty-leaf veto,
the forced collapse, directions drawn only on a flagged column), and equal codes give equal partitions.
A leaf-value edit keeps the skip; a structure edit, a changed code, grid or flag loses it. The hash is
word-at-a-time, not weightsDigest's byte loop: about 2 MB of codes at n = 1e5, p = 10 (xint_t is 2 bytes),
estimated well under 1 ms (not run; measured in B7). Absent attribute (a state written before): full
verdict. The attribute is registered beside
[`stateFormatVersion`](../../src/R_interface_bartcore.cpp); the floor does not move. Open call 8.

B4. The kept store is clean (dec-B318). `stateIsValid` takes a saved block whose size is a multiple of
the tree count, every forest and the variance forest naming one state capacity S; mismatched blocks stay
refused (the stripped-variance case of [test-heteroscedastic.R](../../inst/tinytest/test-heteroscedastic.R)).
On a live install with the sampler's capacity C:
- S equal to C: the ring is copied as today, slot for slot, `currentSampleNum_` and `recordedDraws_`
  from the state, so a copy continues its source's layout bit for bit.
- S not C: with R = min(recordedDraws, S) and the state's draw d at slot (cur + S - R + d) mod S
  ([`savedSlotForDraw`](../../src/bartcore/sampler.hpp)'s rule), the newest K = min(R, C) draws go to
  slots 0 to K - 1, oldest first, with their params, masks and variance trees; `recordedDraws_` = K,
  `currentSampleNum_` = K mod C; n.samples and C do not change. C of 0 takes none.
- A re-creation (copy, reload) still adopts S ([`resizeSavedTrees`](../../src/bartcore/sampler.hpp)),
  as dec-B294 and dec-A188 have it.
`lengthscaleStateFeasible` is asked about saved gp draws only where K > 0.

B5. The spread (dec-B384). k.scale is the leaf scale times the transform's multiplier times sqrt(m)
(`priorScaleFactor`), so under a k-named prior (model@prior.scale NA) it moves with the response units
and under an sd-named one it is twice the sd whatever the units. The units pass gives each chain its
ratio r (state multiplier over the sampler's); `Chain::setState` installs `k = fs.k / r` where the forest
draws k, the state holds one and `kFollowsUnits`, and `fs.k` otherwise. R passes `kFollowsUnits =
is.na(model@prior.scale)`. A forest at a fixed or map-pinned k, or a state without a k block, is
untouched, as [`scaleDrawnK`](../../src/bartcore/chain.hpp) leaves them. The state's k.scale is read as
the recipient's prior under the state's units (Open call 7).

B6. Bridge and R.
- [`bartcore_setState`](../../src/R_interface_bartcore.cpp) takes `force` and `kFollowsUnits` (arity 6 in
  R_interface.cpp) and returns TRUE when it installed, FALSE when unforced and not clean; every refusal
  stays an error. reapplyWeights and reapplySurvivalStatus run only after an install.
- [`installStateOnto`](../../R/dbarts.R) takes `force`; on FALSE it returns FALSE before
  reapplyForestWeights, reapplyActiveRows and reissueNamedLeafSd.
- `setState(newState, forceUpdate = FALSE)`: `forceUpdate` a single TRUE or FALSE, else "'forceUpdate'
  must be TRUE or FALSE" ("partial" included); on FALSE neither the pointer (a dead one stays dead) nor
  `$state` is assigned and FALSE is returned visibly; on a clean unforced install TRUE visibly; forced,
  NULL invisibly (dec-B310). Visibility is Open call 5.
- getPointer and copy pass `force = TRUE`, silent (dec-B305, dec-B234).

B7. Speed. bench-sampler.R gains a restore scenario (store, forced setPredictor, run, put the predictors
back, unforced setState) at n = 1000 and 1e5, and 02-restore-cost.R is rerun against the base build on
one quiet machine: the skip path within 1.10 of the base restore; the full-verdict path recorded, not
gated.

## Refusals and messages

Unchanged, in both forms, the sampler untouched: not a bartcoreState, an older encoding, another chain
count, another tree, forest or variance-forest count, another leaf model, a malformed block or tree, a
repeated grid point (dec-B300), latents out of range (dec-A174's floor), units that cannot be
converted, saved gp draws under other lengthscales, a forbidden column, an interaction limit, a
monotone order out of the cone as stored. Messages as built, "state is not consistent with this
sampler" among them. New: "'forceUpdate' must be TRUE or FALSE"; a one-way categorical rule on a
column never missing (A4) reads as the other malformed trees. No longer refused: a kept store of another size.

## Callers

Package (ran: grep): copy and getPointer, forced; no other R caller. bart's warm.start goes through
installTrees, unchanged here. benchmarks/R: sbc.R (calibration), surfaces/C1-frozen-ess.R and
negbin-mixing.R, an exact gate's script, which the design missed: each hands a sampler its own state
with a scalar edited, so each is clean and returns TRUE, its draws unchanged; none reads the value.

stan4bart (bartcore 963956b, ran: grep): one call, `restoreBartSampler` in R/generics.R, reached from
getBartSampler and stan4bart_fit.R, value dropped. Its restore is clean under this plan (a conversion,
dec-A167; the design's probe 05 and the critique's crit-09 ran it TRUE on the prototype). Changed on its
bartcore branch to stop when the value is FALSE ("the saved BART state does not fit the rebuilt
sampler"), and its suite run against the build. B5 can move its drawn k on a restore across units.

bartCause (dbarts-1.0 e833be7), treatSens (dbarts-1.0 babfaa6) and bairrtt (main 3f57f61): no setState,
copy, storeState or installTrees (ran: grep). bairrtt calls setPredictor forced and the joint row update
on a latent trait with no missing values (read), so Part A reaches its code path and not its draws.

## Tests

tinytest, new:
- test-setstate-force-update.R. Clean, each unforced TRUE visibly and forced NULL invisibly, then the two
  twins' next 20 draws identical: the sampler's own state; one stored before a predictor change that
  empties no leaf; other units; other weights (Student-t) and censoring (aft); stores 10 into 4, 4 into
  10, 4 into none, none into 4; a filled column's directions. Not clean, each unforced FALSE visibly with
  getTrees, sigma, k, kept-draw predictions, cut points, `$state` and the next 20 draws identical to a
  twin that made no call, then forced NULL with today's repaired install (for the leaf no row reaches,
  the trees a forced setPredictor gives): a leaf no row reaches (mean, and variance forest); a split
  outside its interval; a direction on a column never missing, from another sampler's state. A dead
  pointer stays dead after a refused call. `forceUpdate` NA, "partial" and c(TRUE,
  FALSE) refused.
- Store: 10 into 4 predicts as the source's last 4 draws; a wrapped source ring keeps its newest in
  order; 4 into 10 then 3 sweeps holds 7; equal capacity, a copy run 5 sweeps and its source run 5 store
  the same state.
- Spread: across units under chi(1.5, 2), k.scale / getK() per chain equals the source's to 1e-12;
  under an sd-named drawn prior, under a fixed k and with the k block removed, getK unchanged.
- Skip (tests/cpp, through a count of skips the sampler keeps for tests): equal hash, the skip taken and
  the full verdict agreeing; a grid, a flag or one code changed, or a split value edited, the full verdict
  runs; a leaf value edited, the skip still taken.
- test-missingness-first-seen.R: forced, unforced, column and per-observation updates raise the flag and
  draw; over 40 seeds the share of directions drawn right is within a binomial 1e-3 band of one half; a
  second missing-bearing update draws nothing; an unforced first missing value refused (an emptied leaf)
  leaves the sampler bit for bit, the generator included; filling the column keeps every direction and
  predict accepts a missing value there; a column never missing still refuses (dec-B322); copy and
  reload of a filled sampler continue it bit for bit; dec-B378 refused on a never-missing column, built
  and merged on a flagged one.

tinytest, changed (the design's prototype failed the first five, 130 expectations, ran there):
[test-state-empty-leaf-merge.R](../../inst/tinytest/test-state-empty-leaf-merge.R),
[test-state-missing-direction.R](../../inst/tinytest/test-state-missing-direction.R) (TRUE where no leaf
merges), [test-monotone-unforced.R](../../inst/tinytest/test-monotone-unforced.R),
[test-heteroscedastic-mutation.R](../../inst/tinytest/test-heteroscedastic-mutation.R),
[test-state-not-model.R](../../inst/tinytest/test-state-not-model.R) (conversion TRUE), and those that
assert a merged or split-moved install's value:
[test-cut-grid-distinct.R](../../inst/tinytest/test-cut-grid-distinct.R),
[test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R). Each merged install takes
`forceUpdate = TRUE` and gains its unforced FALSE. The full suite finds the rest.

tests/cpp:
- the sticky flag and the first-sight draw in [test_data.cpp](../../tests/cpp/test_data.cpp) and
  [test_sampler.cpp](../../tests/cpp/test_sampler.cpp); the pins that read the flag going down change:
  "drop every missing value: hasMissing flips false" in
  [test_moves.cpp](../../tests/cpp/test_moves.cpp) and
  [`testPerObservationMissingCommit`](../../tests/cpp/test_moves.cpp)'s neighbours;
- in [test_state.cpp](../../tests/cpp/test_state.cpp): the unforced refusal leaves `getState` bytes
  and the generator equal; verdict and install agree over the fuzz harness's states; the store repack;
  k across units.

Mutants the reviewer runs, each failing a test: the flag assigned from content; no first-sight draw, or
one on every update; the generator not put back on a refused first missing value; the record not written
or not read at creation; dec-B378's rule built and merged; install then report FALSE; the verdict
blind to the variance forest; a conversion unclean; case c clean; forced TRUE or visible; the pointer or
`$state` bound on FALSE; the oldest draws kept, or the ring repacked at equal capacity; k divided under
an sd-named prior or where the state has no k; the hash without the flags or the grid.

## Equivalence and snapshots

Expected: nothing recorded moves. No equivalence scenario calls setState (ran: grep; bcf-equivalence.R
and multinomial-equivalence.R say a restore is deliberately not recorded), and the missing-value
scenarios (`missing`, `nafactor`, `testswap` in equivalence.R, the only places any of the three
harnesses writes NA, ran: grep) put their missing values in at creation and mutate no training predictor
(read). The four reproducibility files call neither setState nor setPredictor
(ran: grep). The exact gates' one caller, negbin-mixing.R, is a clean own-state install. So both parts
run against the MANIFEST's current files with every scenario "identical draws (same RNG stream)" (55,
15 and 11 at the last engine landing) and no |z| line; any move is a stop.

## Docs and NEWS

- [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the [`setState`](../../man/dbartsSampler-class.Rd)
  usage, the newState and forceUpdate items (forceUpdate is shared with setPredictor), the Saving
  subsection led by the checked form (`if (!sampler$setState(st)) ...`), TRUE stated first as not bit
  for bit (trees and leaf values as stored or converted, drawn scalars, the newest kept draws that fit;
  latents may be redrawn), copy and reload forced and silent, and Value. The docstrings of setState and
  copy in [dbarts.R](../../R/dbarts.R) (rc-codoc).
- [NEWS.Rd](../../inst/NEWS.Rd), changes from 0.9-34 only: setState gains forceUpdate, unforced by
  default, and returns TRUE or FALSE (0.9-34 installed and returned NULL); a state with more or fewer
  kept draws than the sampler keeps installs what fits (0.9-34 took a larger or smaller store and crashed
  R on a state with none, dec-B318); a first missing value given to a
  column has its splits' directions drawn (0.9-34 sent them down an arbitrary branch).
- [mia-missingness.md](../design/mia-missingness.md) (A6); [state-not-model.md](../design/state-not-model.md),
  the k row of its first table (installed as the spread). A short docs/design/install-surface.md, the
  rules this plan builds, with its INDEX row; the scratch design stays scratch.

## Gates

Per part, on its tip and its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting):
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`); stan4bart's and bairrtt's suites against the build.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) and the four arm64 snapshot files.
- Every exact gate in quick mode; `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc,
  win-drift, doc-freshness.
- Part B: B7 on a quiet machine; Part A: bench-sampler.R compare (the setPredictor scenarios).

## Budget and stops

Planned against landed on this surface ran 1.1 to 2.9 times (the critique's table of six plans); the
last engine slice ran 1.76 times ([small-rulings-1008-engine.md](small-rulings-1008-engine.md), 950 to
1675). So expect about 3000 to 3400 lines. Stop and report, without working around, when:
- Part A passes 1100 lines or Part B 1600;
- any equivalence scenario or snapshot moves;
- the verdict and the install disagree on a fuzz state and one routine cannot be made to serve both;
- the skip (B3) passes 200 lines, or the skip path misses 1.10: land B without it and file the skip;
- a reading needs a state-format floor change, or a rewritten test needs a call no ruling or this plan
  makes.

## Interactions

- The probit k move (probit-k-mixing, fix round in flight) writes `Forest::k` in chain.hpp; B5 writes it
  once at install. Serial: the second to land rebases.
- gp-copy-continuation (in flight) reorders a gp leaf's members on a restore; B4 and B2 do not touch the
  gp draw. Serial.
- state-frame-prior's grid half later moves the grid out of the state; B3's hash then covers the
  sampler's own grid and B2 gains "a split off the grid".
- warm-start-salvage builds installTrees's force argument, its NULL (dec-B310) and the collapse at a
  forbidden split, reusing B2's routine.

## Open calls

1. Two landings, Part A first. Recommend: yes; each reviewed and gated on its own, A making B's verdict
   simpler.
2. Where the missing-value record lives across a re-creation. Recommend: an attribute on data@x, read at
   creation and by the predict refusal, as data derived and held (dec-B231's class); alternatives: a
   state block (a state would then install data, which state-frame-prior removes) or a model attribute.
3. A direction on a column this sampler has never seen missing, in a state from another sampler.
   Recommend: dropped, not clean (a tree that would have to change, dec-B305); alternatives: dropped as
   clean (the bit routes nothing here), or the install raising the column's flag (a state setting a
   property of the sampler).
4. setData and missingness. Recommend: as built for now, the flag from the new data (dec-B321 leaves
   setData out); the alternative, the same first-seen rule, is a small follow-up if ruled.
5. Whether an unforced verdict prints. Recommend: visible, as setPredictor's unforced value is (dec-B305
   asked for setPredictor's idiom); today's setState is invisible and 0.9-34's returned NULL.
6. The first argument's name. Recommend: keep `newState`, 0.9-34's; dec-B305's entry writes `state` in
   its own words, not the maintainer's.
7. dec-B384's k.scale of the state. Recommend: the recipient's prior read in the state's units (B5), no
   format change; a state stored under another prior keeps its k. The alternative is a state block
   holding the spread, which carries model in the state.
8. The skip ships in this slice. Recommend: yes; without it every unforced restore pays up to about 1.3
   sweeps more on the package's hot path.
9. installTrees's NULL (forced-update-returns-null's other half). Recommend: with warm-start-salvage,
   where its force argument is built.
10. A constraint the state's tree breaks (column, interaction, monotone) stays an error in setState.
   Recommend: yes (the critique, against the design's needs-force reading).

## Evidence

Ran on 17df64ff (private library scratch/libs/isplan; probes and outputs in scratch/isplan/): the
behaviours in Context (01-behaviours.R); the restore costs in B3 (02-restore-cost.R, under load); the
greps for callers, test call sites, harness and snapshot calls. Read: every engine and bridge claim
cited by symbol, the TODO entries, and dec-B283 to dec-B397 and dec-A160 to dec-A189. Bounds not run:
the scratch verdict's added cost and the hash's cost.
