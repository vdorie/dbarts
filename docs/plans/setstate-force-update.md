# setstate-force-update: setState's two forms, the kept store, the spread, and missingness first seen

Status: PLANNED (dec-B305, dec-B310, dec-B318, dec-B320, dec-B321, dec-B322, dec-B378, dec-B384). The
first slice of the install surface; the warm start (warm-start-salvage, warm-start-tree-count) is later.
Revised from its blind critique (Revision from critique); five questions are held for the maintainer
(Held for the maintainer) and the plan builds either way.

agent: opus implementer, one; one opus reviewer who runs the mutants below, per part.
rng: Part A SHIFTING for a sampler whose column first gets a missing value after creation, or loses its
missing values (each chain's generator draws directions it did not draw before); NEUTRAL, bit for bit, for
every sampler whose columns' missingness never changes after creation. Part B NEUTRAL for a clean
install, a copy and a reload in the sampler's own store size and units; SHIFTING where a drawn k installs
across response units (dec-B384) or a store of another size installs (refused today). The held branches
add shifts where they say so.
window: before the merge to main; engine slices stay serial.
budget: planned ~1580 lines in two landings, Part A ~730 (engine ~200, bridge and R ~110, tests ~350,
docs ~70), Part B ~850 (engine ~170, bridge and R ~80, tests ~470, help, NEWS and docs ~130), before the
held branches (each states its own lines). Forecast and stops: Budget and stops.

## Goal

`setState(newState, forceUpdate = FALSE)` installs a state that goes in cleanly and returns TRUE, or
leaves the sampler as it was and returns FALSE; forced, it installs, repairing what it must, and returns
NULL invisibly. A conversion of units, latents redrawn under other weights, a kept store of another size
and a missing direction on a column that has been missing are clean. Whether a column can be missing is a
property of the sampler, raised by the first missing value it is given and never lowered by a predictor
change, and that first value draws each existing split's direction at one half. A drawn k installs as the
spread in force in the state. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Inputs, untracked (dec-B389: near-term work): the design scratch/install-surface/design.md (2026-10-07),
its critique scratch/install-surface-critique/critique.md, this plan's blind critique
scratch/iscrit/critique.md and the orchestrator's calls on it, scratch/iscrit/orchestrator-calls.md
(2026-10-09, recorded in the ledger at landing). The design predates the rulings below; where they
differ the ruling wins:

| Design said | Ruled since | This plan |
|---|---|---|
| forced FALSE when repaired | dec-B310: forced forms return NULL invisibly | NULL, invisibly |
| a dropped missing direction needs force | dec-B320: kept, clean; dec-B321: first seen, ever after | Part A |
| a larger kept store needs force | dec-B318: clean, newest kept, TRUE | Part B, step B4 |
| (not covered) | dec-B384: k re-expressed against the recipient's k.scale | Part B, step B5 |
| (not covered) | dec-B378: a one-way categorical rule refused on a column never missing | Part A, step A4 |
| (not covered) | dec-B322: predict refuses a missing value only where the column could never be missing | Part A, step A5 |
| grid half (slice 4) | the design's critique: out; TODO state-frame-prior | out: the grid arrives with the state, as built |

That these items go in "the install surface's first slice" is the register author's scheduling (dec-B318,
dec-B320, TODO missingness-first-seen and state-install-keeps-spread), not the maintainer's words; the
plan follows it because Part A is what dec-B320's clean rule needs (a sticky flag and a record that
survives a copy) and B4 and B5 are recorded for this slice.

Settled (orchestrator, 2026-10-09):
- Two landings, Part A first; each reviewed and gated on its own.
- The argument is named `forceUpdate`, setPredictor's, which settles dec-B305's "one name for both":
  warm-start-salvage gives installTrees `forceUpdate` too.
- The first argument keeps 0.9-34's name, `newState` (dec-B305's entry writes `state` in its author's
  words, not the maintainer's).
- installTrees's forced NULL (dec-B310, TODO forced-update-returns-null's other half) is built with
  warm-start-salvage, where its force argument is.
- The missing-value record is a dbartsData slot (A3).
- No verdict skip in this slice (B3).

Out, as the design's critique had it: one error text per cause; a new floor on leaf values and sigma (TODO
precision-floor-edge-values); the grid half. An engine `force` argument defaults so tests/cpp call sites
stay; the verdict leaves a dead pointer dead; the help leads with the checked form; the verdict and the
install share one routine (B2). Whether a constraint break is an error or a repair under force is held
(H1).

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

Run on the tip (scratch/isplan/01-behaviours.R and 02-restore-cost.R, library scratch/libs/isplan; the
critique reran both, rerun-01 identical, rerun-02 within a few percent): setState returns TRUE invisibly
on its own state; a state of 10 kept draws into a store of 4, 4 into 10 and 4 into none are refused
("state is not consistent with this sampler") and none into 4 installs, TRUE; a state stored while a
numeric column held missing values, restored after a forced setPredictor filled them, gives FALSE; a
first missing value by an unforced setPredictor is taken (TRUE) and predict then answers for it; a state
in other response units gives FALSE; across units a drawn k installs as stored, k.scale 7.83 to 23.48 and
the spread 2.47 to 7.40.

## Part A: missingness first seen (dec-B320, dec-B321, dec-B322, dec-B378)

Lands first: once a column's flag stays raised, a direction is never dropped on the sampler's own trees,
and Part B's verdict has one direction case left (B2, case c).

A1. The flag is sticky. Every train-side quantize ORs into `hasMissing[j]` instead of assigning it;
creation sets it from content (as today) or from the record (A3). The rollback snapshots of
`WholeMatrixUpdate` and `SubsetUpdate` ([sampler.hpp](../../src/bartcore/sampler.hpp)) already put the
flag back on a refusal. setData is held (H2): either it resets the flags from the new data's content
before it re-quantizes, as built, or it keeps them sticky. The test store tracks no flag.

A2. The first missing value draws. Where an accepted predictor change raises a column's flag from 0,
each chain draws, from its own generator, a direction at one half for every split already on that
column on its live trees, in a fixed order: forests in order, trees in order, nodes in pre-order, then the
variance forest's trees, with `ext_rng_simulateBernoulli(rng, 0.5)` as the birth draw takes it. Whether
the kept draws' trees draw too is held (H3). Sites:
- [`runPredictorTransaction`](../../src/bartcore/sampler.hpp), forced: after `applyForced`, before
  `forceRefreshTrees`.
- unforced: after `snapshotApply` and before [`revalidateAllChains`](../../src/bartcore/sampler.hpp),
  since occupancy and the monotone order are judged under the drawn directions. A rollback puts back the
  bits drawn and each chain's generator (serialized before the draw, as
  [`Chain::setState`](../../src/bartcore/chain.hpp) reads one back; each chain owns its generator, so
  the round trip is exact).
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
from data, which may hold no missing value in a column the sampler has seen missing; without a record the
copy would hold other flags than its source and draw differently. The record is a new dbartsData slot,
`missing.seen`, beside `rowNames` in [A_class.R](../../R/A_class.R): NULL, or a logical per predictor
column, TRUE where the sampler's flag is raised. dbartsData is where TODO state-frame-prior puts the other
records derived from data and held, and a slot survives what a matrix attribute does not (subsetting
drops an attribute, and on a dgCMatrix `[<-` does too; ran in the critique).
- Read through one helper that answers NULL for a data object saved before the slot existed (the
  [`dataRowNames`](../../R/data.R) pattern, `methods::.hasSlot`); validity: NULL or a logical of
  ncol(x) with no NA.
- Written R-side after every accepted predictor change (every path of
  [`bartcoreSamplerSetPredictor`](../../R/bartcore.R), the joint row update) from a new bridge reader of
  the engine's flags, and after setData by the rule H2 takes, so a data object handed to setData never
  keeps a record the call has reset (`d <- s$data; d@y <- y2; s$setData(d)`).
- Read by the bridge at creation ([`bartcore_create`](../../src/R_interface_bartcore.cpp) and the handle
  path), ORed into the content flags; and by the predict refusal (A5).
- A fit object changes: a fitted sampler's data carries the slot, so the exact gates run in quick mode at
  Part A's landing (CLAUDE.local.md's rule for a change to what a fit object carries), and the slot is
  described in [dbartsData.Rd](../../man/dbartsData.Rd)'s Value paragraph beside `rowNames`. A fit saved
  before has no slot and reads as none.

A4. No direction is dropped afterwards, and dec-B378. With A1, `dropStaleMissingDirections` clears
nothing on a column that has been missing; its sites stay for setData (H2) and for a flat tree from
another sampler. In [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) the gauge on a column whose flag
is down keeps counting the missing position reachable (`missingReaches`), but a categorical rule that,
with its missing bit dropped, sends every reachable level one way is refused as malformed (dec-B378),
where today it is built and merged; on a flagged column it is built, and merged where no row is missing
now (not clean, B2). A numeric or inline categorical rule with a direction on a column whose flag is down
loses it, reported as today through `directionDropped` (Open call 1). What this refusal does after a
setData is H2.

A5. predict (dec-B322). [`unroutableTestColumns`](../../R/data.R) and
[`refuseTestMissingness`](../../R/data.R) take the record beside the training matrix (their callers hold
the data object), so a column made missable by setPredictor, then filled, is still routable; a column
never missing is refused as built. The engine's own check (store.hasMissing at the predict entries)
follows A1 with no change.

A6. Docs: [mia-missingness.md](../design/mia-missingness.md) takes an amendment under its Status
paragraph, in [Representation](../design/mia-missingness.md#representation) ("cleared when a column
loses its NAs") and in [Bridge and R surface](../design/mia-missingness.md#bridge-and-r-surface) ("a
column that gains NAs mid-run routes them by ... (left)"), both becoming the first-seen rule.

## Part B: setState's two forms (dec-B305, dec-B310, dec-B318, dec-B384)

B1. The flag. `Sampler::setState` (and the facade virtual on
[`SamplerBase`](../../src/bartcore/facade.hpp)) take, after `adoptCapacity`, `bool force = true`, `bool
kFollowsUnits = false` and `bool* notClean = nullptr`, so the tests/cpp call sites (88 lines matching
`setState(`, chain-level calls among them; the critique's count) keep today's behaviour. `--preclean` on
every install: a virtual changes.

B2. Clean, in code. After every refusal (unchanged) and with only the state's grid written, an unforced
call asks each chain for a verdict and, if any chain is not clean, puts the grid back, sets `*notClean`
and returns false, the chains, store, generators, flags and latents untouched (and, in R, the mirrors,
B6). A chain is not clean when any live mean tree or variance tree, built from its flat form and
partitioned over the sampler's rows:
- (a) has a bottom node no row reaches;
- (b) holds a split outside the interval its ancestors leave;
- (c) carries a missing direction on a non-pooled column this sampler has never seen missing.

Everything else that differs from the state is clean: the units pass (`valuesMoved` stops setting
`altered`), latents redrawn under other weights or censoring (`reapplyWeights`, `reapplySurvivalStatus`,
after the install as today), a store of another size (B4), k re-expressed (B5), DART weights, a fixed
value or a generator of another kind left as the sampler's. Under H1's repair branch a constraint break
joins (a) to (c).

One routine: a new Tree-level staging call (build from flat, partition, report a, b and c) used by the
verdict on a scratch tree with one n-length index buffer, as `stateIsValid` builds today, and by
`rebuildLiveForest` and `rebuildVarianceForest` on the live tree, so the verdict and the install's report
cannot disagree; a tests/cpp check asserts it over the fuzz harness's states (unforced false exactly when
the forced install reports a repair). Forced, the verdict is skipped and the install is today's.

B3. Its cost. Measured on the tip (ran: 02-restore-cost.R, one chain, one thread, 200 trees, 10 columns,
machine under load 5 to 6; rerun by the critique within a few percent): a clean restore of a fitted state
costs 74.0 ms against a 41.0 ms sweep at n = 1e5 (1.81 sweeps) and 3.3 against 2.1 ms at n = 5000 (1.53);
restoring single-leaf trees costs 19.2 and 1.0 ms, so the per-tree build and partition is 74 and 70
percent of a restore, 1.34 and 1.07 sweeps. The scratch verdict repeats that work less the fits, so an
unforced restore costs up to about 1.1 to 1.3 sweeps more (bound, not run; B7 measures it).

Weighed and not taken:
- building each tree once into spare trees over a second index buffer and swapping them in on a clean
  or forced verdict: one partition, at a second n x m buffer of 4-byte indices per chain, 80 MB a chain
  at n = 1e5, m = 200, against the memory the large-n work holds down;
- a hash in the state that lets a restore onto unchanged rows skip the verdict (dropped by the
  orchestrator: it answers no ruling, its benefit is unmeasured, and storeState, which runs after every
  run() under updateState = TRUE, would pay about 70 percent more; the critique measured 0.33 ms on a 0.46
  ms storeState). A TODO note records it, to revisit only with a measured need.
A caller that restores its own state onto rows it has put back, where nothing can need repair, can pass
`forceUpdate = TRUE` and pay today's restore; the help says so.

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
untouched, as [`scaleDrawnK`](../../src/bartcore/chain.hpp) leaves them. A state stored under another
prior is H4.

B6. Bridge and R.
- [`bartcore_setState`](../../src/R_interface_bartcore.cpp) takes `force` and `kFollowsUnits` (arity 4
  to 6 in R_interface.cpp) and returns TRUE when it installed, FALSE when unforced and not clean; every
  refusal stays an error. reapplyWeights and reapplySurvivalStatus run only after an install.
- [`installStateOnto`](../../R/dbarts.R) takes `force`; on FALSE it returns FALSE before
  reapplyForestWeights, reapplyActiveRows and reissueNamedLeafSd.
- `setState(newState, forceUpdate = FALSE)`: `forceUpdate` a single TRUE or FALSE, else "'forceUpdate'
  must be TRUE or FALSE" ("partial" included); on FALSE neither the pointer (a dead one stays dead) nor
  `$state` is assigned and FALSE is returned; on a clean unforced install TRUE; forced, NULL invisibly
  (dec-B310). Whether the unforced value prints is Open call 2.
- getPointer and copy pass `force = TRUE`, silent (dec-B305, dec-B234).

B7. Speed, recorded. bench-sampler.R gains a restore scenario (store, forced setPredictor, run, put the
predictors back, unforced setState, and the same forced) at n = 1000 and 1e5, and 02-restore-cost.R is
rerun against the base build on one quiet machine. The unforced restore's ratio to the base restore and
the forced restore's (expected 1.00) go in the landing note; neither gates the landing. storeState and
run() do not change in this slice.

## Refusals and messages

Unchanged, in both forms, the sampler untouched: not a bartcoreState, an older encoding, another chain
count, another tree, forest or variance-forest count, another leaf model, a malformed block or tree, a
repeated grid point (dec-B300), latents out of range (dec-A174's floor), units that cannot be converted,
saved gp draws under other lengthscales; and, under H1's error branch, a forbidden column, an interaction
limit and a monotone order out of the cone as stored. Messages as built, "state is not consistent with
this sampler" among them. New: "'forceUpdate' must be TRUE or FALSE"; a one-way categorical rule on a
column never missing (A4) reads as the other malformed trees. No longer refused: a kept store of another
size.

## Callers

Package (ran: grep): copy and getPointer, forced; no other R caller. bart's warm.start goes through
installTrees, unchanged here. benchmarks/R: sbc.R (calibration) and negbin-mixing.R, an exact gate's
script, which the design missed, each hand a sampler its own state with a scalar edited;
surfaces/C1-frozen-ess.R transplants a recorded chain's state into a fresh sampler built on the same data
with levelGibbs. Each is clean (same rows, same units) and returns TRUE, its draws unchanged; none reads
the value.

stan4bart (bartcore 963956b, ran: grep): one call, `restoreBartSampler` in R/generics.R, reached from
getBartSampler after a reload and from stan4bart_fit.R at fit end, value dropped. Its restore is clean
under this plan (a conversion, dec-A167; the design's probe 05 and the first critique's crit-09 ran it
TRUE on the prototype). What it changes to is H5; its suite runs against the build either way. B5 can
move its drawn k on a restore across units.

bartCause (dbarts-1.0 e833be7), treatSens (dbarts-1.0 babfaa6) and bairrtt (main 3f57f61): no setState,
copy, storeState or installTrees (ran: grep). bairrtt calls setPredictor forced and the joint row update
on a latent trait with no missing values (read), so Part A reaches its code path and not its draws.

## Tests

tinytest, new:
- test-setstate-force-update.R. Clean, each unforced TRUE and forced NULL invisibly, then the two twins'
  next 20 draws identical: the sampler's own state; a recorded chain transplanted into a fresh sampler;
  one stored before a predictor change that empties no leaf; other units; other weights (Student-t) and
  censoring (aft); stores 10 into 4, 4 into 10, 4 into none, none into 4; a filled column's directions.
  Not clean, each unforced FALSE with getTrees, sigma, k, kept-draw predictions, cut points, `$state` and
  the next 20 draws identical to a twin that made no call, then forced NULL with today's repaired install
  (for the leaf no row reaches, the trees a forced setPredictor gives): a leaf no row reaches (mean, and
  variance forest); a split outside its interval; a direction on a column never missing, from another
  sampler's state. A dead pointer stays dead after a refused call. `forceUpdate` NA, "partial" and
  c(TRUE, FALSE) refused.
- Store: 10 into 4 predicts as the source's last 4 draws; a wrapped source ring keeps its newest in
  order; 4 into 10 then 3 sweeps holds 7; equal capacity, a copy run 5 sweeps and its source run 5 store
  the same state.
- Spread: across units under chi(1.5, 2), k.scale / getK() per chain equals the source's to 1e-12;
  under an sd-named drawn prior, under a fixed k and with the k block removed, getK unchanged.
- test-missingness-first-seen.R: forced, unforced, column and per-observation updates raise the flag and
  draw; over 40 seeds the share of directions drawn right is within a binomial 1e-3 band of one half; a
  second missing-bearing update draws nothing; filling the column keeps every direction and predict
  accepts a missing value there; a column never missing still refuses (dec-B322); copy and reload of a
  filled sampler continue it bit for bit; a data object saved before the slot existed reads as no
  record; dec-B378 refused on a never-missing column, built and merged on a flagged one.
- A refused first missing value (an emptied leaf) leaves each chain's generator state identical. On the
  column and per-observation paths the next 20 draws are identical to an untouched twin's; on the
  whole-matrix path they are equal within 1e-12, the rounding TODO refused-update-rounding records there
  (about 1.8e-15, a rollback's repartition), which this slice does not own.

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
one on every update; the generator not put back on a refused first missing value; the record not written,
not read at creation, or kept by setData against H2's branch; dec-B378's rule built and merged; install
then report FALSE; the verdict blind to the variance forest; a conversion unclean; case c clean; forced
TRUE or visible; the pointer or `$state` bound on FALSE; the oldest draws kept, or the ring repacked at
equal capacity; k divided under an sd-named prior or where the state has no k.

## Equivalence and snapshots

Expected: nothing recorded moves. No equivalence scenario calls setState (ran: grep; bcf-equivalence.R
and multinomial-equivalence.R say a restore is deliberately not recorded), and the missing-value
scenarios (`missing`, `nafactor`, `testswap` in equivalence.R, the only places any of the three
harnesses writes NA, ran: grep) put their missing values in at creation and mutate no training predictor
(read); the mutation and setData scenarios use complete data (the critique, read). The four
reproducibility files call neither setState nor setPredictor (ran: grep). The exact gates' one caller,
negbin-mixing.R, is a clean own-state install. So both parts run against the MANIFEST's current files
with every scenario "identical draws (same RNG stream)" (55, 15 and 11 at the last engine landing) and no
|z| line; any move is a stop. This holds on every held branch.

## Docs and NEWS

- [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the [`setState`](../../man/dbartsSampler-class.Rd)
  usage, the newState and forceUpdate items (forceUpdate is shared with setPredictor), the Saving
  subsection led by the checked form (`if (!sampler$setState(st)) ...`) and naming `forceUpdate = TRUE`
  for a loop that restores its own state onto rows it has put back, TRUE stated first as not bit for bit
  (trees and leaf values as stored or converted, drawn scalars, the newest kept draws that fit; latents
  may be redrawn), copy and reload forced and silent, and Value. The docstrings of setState and copy in
  [dbarts.R](../../R/dbarts.R) (rc-codoc). [dbartsData.Rd](../../man/dbartsData.Rd): the slot (A3).
- [NEWS.Rd](../../inst/NEWS.Rd), changes from 0.9-34 only: setState gains forceUpdate, unforced by
  default, and returns TRUE or FALSE (0.9-34 installed and returned NULL); a state with more or fewer
  kept draws than the sampler keeps installs what fits (0.9-34 took a larger or smaller store and crashed
  R on a state with none, dec-B318); a first missing value given to a column has its splits' directions
  drawn (0.9-34 sent them down an arbitrary branch).
- [mia-missingness.md](../design/mia-missingness.md) (A6); [state-not-model.md](../design/state-not-model.md),
  the k row of its first table (installed as the spread). A short docs/design/install-surface.md, the
  rules this plan builds and the held calls as ruled, with its INDEX row; the scratch design stays
  scratch.
- TODO at landing: setstate-force-update, state-store-size-install, state-install-keeps-spread and
  missingness-first-seen close; forced-update-returns-null keeps installTrees's half; a note on the
  verdict skip (B3).

## Gates

Per part, on its tip and its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting):
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`); stan4bart's and bairrtt's suites against the build.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) and the four arm64 snapshot files.
- Every exact gate in quick mode (on Part A also because a fit object gains a slot, A3); `R CMD check
  --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.
- Part A: bench-sampler.R compare (the setPredictor scenarios); Part B: B7, recorded.

## Budget and stops

Planned against landed on this surface ran 1.1 to 2.9 times (the first critique's table of six plans);
the last engine slice ran 1.76 times ([small-rulings-1008-engine.md](small-rulings-1008-engine.md), 950
to 1675). Forecast at 1.76, before the held branches: Part A ~1300 lines, Part B ~1500. The stops are
1.5 times the forecast: Part A ~1950 and Part B ~2250, each held branch taken adding 1.5 times its
forecast lines (below) to its part's stop. Stop and report, without working
around, when:
- a part passes its stop;
- any equivalence scenario or snapshot moves;
- the verdict and the install disagree on a fuzz state and one routine cannot be made to serve both;
- a reading needs a state-format floor change, or a rewritten test needs a call no ruling, settled call
  or held answer makes.

## Held for the maintainer

Each is open; the plan builds whichever way it is ruled. Lines are planned (forecast at 1.76 in
brackets). The first branch of each is the cheaper one, not a recommendation.

H1. A constraint the state's tree breaks, under force. A split on a column the forest may not use, a
split breaking an interactions() limit, or a monotone tree whose leaf values are out of the cone, in the
state as stored. Today each is an error in both forms. dec-B305's force is to "'cram' a state in there and
patch up whatever we could"; its entry's not-clean is "a tree that would have to be changed to fit";
dec-B295, of the warm start: "If trees can be collapsed to be made into valid ones, that would make
sense".
- Error, as now. Scope: nothing added; `columnMaskStateFeasible`, `interactionStateFeasible` and
  `monotoneStateFeasible` stay refusals. Tests: the existing expect_error pins (test-single-forest-vars,
  test-interactions, test-blocks, test-monotone) stand; one each that forced errors too. Lines: +20
  (+35).
- Repair. Scope: the three feed B2 as cases (d) to (f); forced, the first offending split on each path
  from the root becomes one leaf with everything beneath it, at the row-weighted mean of the leaves
  beneath (the merge `rebuildLiveForest` already does), and an out-of-cone tree is reseeded
  (`reseedInfeasibleMonotoneLeaves`, as a warm start does today); the collapse routine is warm-start-
  salvage's, pulled into Part B, which that item then reuses. Tests: the four files' errors become
  unforced FALSE and forced repaired; new: the collapse at the first offending split from the root, its
  value, the variance forest's column mask. Lines: +350 (+620), Part B.

H2. setData and missingness. Under A1 a predictor change never lowers a column's flag; setData is not
covered by dec-B321. Run by the critique (p1): a 4-level factor with 80 missing values, 200 sweeps, 11 of
96 factor rules in the kept draws one-way (every level right, missing left); setData to complete data;
setState of the stored state gives FALSE on the tip (merged).
- setData resets the flags from the new data, as built. Scope: setData writes the record from the
  engine after the call (cleared where content has none). Consequence: after setData the column is
  never missing, so A4 refuses a one-way categorical rule there as malformed, an error in both forms; a
  setState of a state stored before the setData errors where the tip merges, and so do a copy and a
  reload of a sampler whose stored state predates it (updateState = FALSE, the help's setting in a
  loop), at every use. Tests: that error pinned, forced and unforced, and on copy; the help names it.
  Lines: +40 (+70), Part A.
- setData keeps the flags sticky. Scope: setData ORs the old flags into the new data's; a column it first
  makes missing draws at first sight as A2 (its site in applyNewData, under each chain's generator);
  the record carries across. Consequence: no new error; a setData that removes a column's missing values
  keeps the directions and moves keep drawing them, a shift for such samplers (no recorded scenario
  reaches it). Tests: setData then setState of the earlier state clean or merged as B2 says, copy and
  reload of it; the first-sight draw through setData. Lines: +120 (+210), Part A.

H3. Kept draws recorded before a column's first missing value. A2 draws directions on live trees, as
dec-B321 names "each split already on the column"; the kept draws' trees keep every direction left, and
predict answers a missing row in that column from every kept draw (dec-B322). Run by the critique (p3): 20
kept draws, a first missing value, 5 more; in the 15 older draws every split on the column sends missing
left, and predict uses all 20, so 75 percent of that row's posterior mean comes from the fixed left.
- Live trees only, documented. Scope: the help says older kept draws route a missing value left. Tests:
  that share pinned at 0 in older draws. Lines: +20 (+35).
- Kept trees draw too. Scope: at A2's sites, after the live trees, each chain draws for every split on
  the column in its kept slots, oldest first, mean then variance, setting the flat flag; the rollback
  covers them. A shift for such samplers only. Tests: the share in older draws within the binomial band;
  the refused case leaves the store bitwise. Lines: +100 (+175), Part A. (The critique's third option,
  predict refusing or dropping older draws for such rows, needs a per-column record of when the column
  was first seen and a ragged predict; about +250, not planned.)

H4. The spread of a state stored under another prior. B5 takes k.scale_state as the recipient's prior
read in the state's units, which is exact when the state was stored under the prior in force. A state
stored before a setLeafPrior (dec-B393: "Changing the prior never moves the state"), or by a sampler with
another prior, keeps its k, so the spread jumps on install by the ratio of the two k.scale values, the
jump dec-B384 stops across units. The state holds no prior; the spread in force is chain state (dec-B356,
dec-B369, dec-B392 and dec-B393 each keep it across a prior change).
- As B5. Scope: none added; the help names the case. Tests: a setLeafPrior between store and restore
  keeps k, pinned. Lines: +15 (+25).
- The state carries the spread. Scope: an optional per-forest block, the spread in force (k.scale / k,
  response units) where the forest draws k, written by storeState and read by setState, a copy and a
  reload, which set k = k.scale_recipient / spread; absent (an older state) falls back to B5; registered
  beside [`stateFormatVersion`](../../src/R_interface_bartcore.cpp), floor unmoved; `kFollowsUnits`
  becomes unnecessary. Tests: the setLeafPrior case keeps the spread to 1e-12; an older state falls back.
  Lines: +90 (+160), Part B.

H5. stan4bart's restore. restoreBartSampler is stan4bart's own re-creation (a new dbartsSampler from
control, model and data, then the state), the analogue of dbarts's reload, which dec-B305 forces
("copy() and a reload force, having no choice").
- Force: `sampler$setState(state, forceUpdate = TRUE)`, as dbarts's reload does. A restore that is not
  clean is repaired in silence and the verdict's cost is not paid. Scope: one line on stan4bart's
  bartcore branch. Tests: its suite as is. Lines: +2.
- Stop on FALSE: `if (!sampler$setState(state)) stop("the saved BART state does not fit the rebuilt
  sampler")`. A fit whose restore is not clean cannot predict after a reload, where dbarts's own reload
  would repair; the restore pays the verdict. Scope: three lines there. Tests: its suite as is (no clean
  failure is known to provoke). Lines: +4.

## Interactions

- The probit k move (probit-k-mixing, fix round in flight) writes `Forest::k` in chain.hpp; B5 writes it
  once at install. Serial: the second to land rebases.
- gp-copy-continuation (in flight) reorders a gp leaf's members on a restore; B4 and B2 do not touch the
  gp draw. Serial.
- state-frame-prior's grid half later moves the grid out of the state; B2 then gains "a split off the
  grid", and the missing-value record sits beside the other derived records on dbartsData.
- warm-start-salvage builds installTrees's `forceUpdate`, its NULL (dec-B310) and, unless H1 pulls it
  forward, the collapse at a forbidden split, reusing B2's routine.

## Open calls

For the orchestrator; each has a recommendation.
1. A direction on a column this sampler has never seen missing, in a state from another sampler.
   Recommend: dropped, not clean (a tree that would have to change, dec-B305's entry); alternatives:
   dropped as clean (the bit routes nothing here), or the install raising the column's flag (a state
   setting a property of the sampler).
2. Whether an unforced verdict prints. Recommend: visible, as setPredictor's unforced value is (dec-B305
   asked for setPredictor's idiom); today's setState is invisible and 0.9-34's returned NULL.

## Revision from critique

Against aa79ac10, from scratch/iscrit/critique.md and the orchestrator's calls:
- The verdict skip (hash, state attribute, its tests, mutants and stop) is dropped; B3 states the
  measured cost, the spare-buffer swap weighed with its 80 MB a chain, and `forceUpdate = TRUE` as the
  loop's way out; B7 records the cost and gates nothing.
- Stops are 1.5 times the forecast per part (about 1950 and 2250), not below it.
- The refused first-missing-value test is bitwise on the column and per-observation paths and within
  1e-12 on the whole-matrix path, citing TODO refused-update-rounding.
- The record is a dbartsData slot, `missing.seen`, read through a `.hasSlot` helper; the fit-object
  consequence (exact gates at Part A's landing, the Rd) is stated, and setData rewrites the record.
- Settled and moved out of Open calls: two landings with Part A first, `newState`, installTrees's NULL
  with warm-start-salvage, the slot; added: `forceUpdate` is the one name with installTrees.
- Held for the maintainer, with both branches: constraint breaks under force (H1, with the quoted force),
  setData and dec-B378's refusal (H2), kept draws across a first missing value (H3), the spread under
  another prior (H4, its ground restated: the spread is chain state), stan4bart's restore (H5).
- Fixed: C1-frozen-ess.R transplants another sampler's state; mia-missingness.md's "cleared when a column
  loses its NAs" is under Representation; 88 tests/cpp setState lines; the "first slice" scheduling is
  marked as the register author's.

## Evidence

Ran on 17df64ff (private library scratch/libs/isplan; probes and outputs in scratch/isplan/): the
behaviours in Context (01-behaviours.R); the restore costs in B3 (02-restore-cost.R, under load); the
greps for callers, test call sites, harness and snapshot calls. From the critique's runs (scratch/iscrit/):
the store and hash costs, the one-way rules after setData (p1), the kept draws across a first missing
value (p3), attribute loss on subsetting. Read: every engine and bridge claim cited by symbol, the TODO
entries, and dec-B283 to dec-B397 and dec-A160 to dec-A189. Bounds not run: the scratch verdict's added
cost.
