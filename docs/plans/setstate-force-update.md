# setstate-force-update: setState's two forms, the kept store, constraint repair, and missingness first seen

Status: PART A LANDED 2026-10-10 (acf6a6fd..2269918e); Part B PLANNED (dec-B305, dec-B310, dec-B318, dec-B320, dec-B321, dec-B322, dec-B378, dec-B398 to
dec-B402, dec-B418). The first slice of the install surface; the warm start (warm-start-salvage,
warm-start-tree-count) is later. Revised from its blind critique (Revision from critique) and again for
the maintainer's rulings of 2026-10-09 (Revision for the 2026-10-09 rulings); nothing is held.

agent: opus implementer, one; one opus reviewer who runs the mutants below, per part.
rng: Part A SHIFTING for a sampler whose column first becomes able to hold a missing value after
creation (by setPredictor, the per-observation update or setData) or whose column loses its missing
values: each chain's generator draws directions, on live trees and kept draws, that it did not draw
before, and directions are kept where today they are dropped. NEUTRAL, bit for bit, for every sampler
whose columns' missingness never changes after creation. Part B NEUTRAL for every install, copy and
reload that goes in today; new draws only where an install is refused today (a kept store of another
size, a constraint break forced). What an install across a change of response mapping does is the
k-internal slice's (dec-B418), not this plan's.
window: before the merge to main; engine slices stay serial. Its place: after the k-internal slice
(dec-A191).
budget: planned ~2070 lines in two landings, Part A ~970 (engine ~290, bridge and R ~125, tests ~460,
help and docs ~95), Part B ~1100 (engine ~325, bridge and R ~75, tests ~560, help, NEWS and docs ~140).
Forecast and stops: Budget and stops.

## Goal

`setState(newState, forceUpdate = FALSE)` installs a state that goes in cleanly and returns TRUE, or
leaves the sampler as it was and returns FALSE; forced, it installs, repairing what it must, and returns
NULL invisibly. A kept store of another size, latents redrawn under other weights, a missing direction on
a column that has been missing and an install across a change of response mapping are clean. A tree that
breaks the recipient's column restriction, interaction limit or monotone cone is not clean and, forced,
is repaired: collapsed at the first offending split, or reseeded. Whether a column can be missing is a
property of the sampler, raised by the first missing value it is given, by any path, and never lowered;
that first value draws the direction of each split already on the column at one half, in the live trees
and in the kept draws. A state's k installs as recorded. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

Inputs, untracked (dec-B389: near-term work): the design scratch/install-surface/design.md (2026-10-07),
its critique scratch/install-surface-critique/critique.md, this plan's blind critique
scratch/iscrit/critique.md and the orchestrator's calls on it, scratch/iscrit/orchestrator-calls.md
(2026-10-09, recorded in the ledger at landing). The design predates the rulings below; where they
differ the ruling wins:

| Design or earlier plan said | Ruled since | This plan |
|---|---|---|
| forced FALSE when repaired | dec-B310: forced forms return NULL invisibly | NULL, invisibly |
| a dropped missing direction needs force | dec-B320: kept, clean; dec-B321: first seen, ever after | Part A |
| setData resets the flags (as built) | dec-B399: sticky | A1, A2 |
| kept draws keep their left directions | dec-B400: they draw too; the help says so | A2, A6 |
| a larger kept store needs force | dec-B318: clean, newest kept, TRUE | B4 |
| a conversion of units is clean; k re-expressed (dec-B384, dec-B407) | dec-B418: nothing converted on install, the install clean; reverses dec-B200's conversion, dec-B384 and dec-B407 | the k-internal slice's; nothing here |
| constraint breaks stay errors | dec-B398: repaired under force | B2 cases d to f, B8 |
| (not covered) | dec-B401: k transfers literally across a change of prior | nothing built; one help sentence |
| stan4bart stops on FALSE | dec-B402: forces | Callers |
| (not covered) | dec-B378: a one-way categorical rule refused on a column never missing | A4 |
| (not covered) | dec-B322: predict refuses a missing value only where the column could never be missing | A5 |
| grid half (slice 4) | the design's critique: out; TODO state-frame-prior | out: the grid arrives with the state, as built |

That these items go in "the install surface's first slice" is the register author's scheduling (dec-B318,
dec-B320, TODO missingness-first-seen), not the maintainer's words; the plan follows it because Part A is
what dec-B320's clean rule needs (a sticky flag and a record that survives a copy) and B4 is recorded for
this slice.

Built on the k-internal slice. dec-A191 orders the engine slices gp-copy (landed, dec-A192), the
width-weighted default cut rule, docs/plans/k-internal-parameter.md (branch wt/k-internal-plan, not on
this branch) carrying dec-B418, then this plan's two parts, then response-scale-rows. dec-B418, the
maintainer: "Does restoring across a re-anchor even make sense? The mental model of the sampler is a
linear one through time, with each sample being a snapshot." and "OK, let's do it."; the rule, in the
register's words: setState, copy() and a reload "install the chain's parameters ... as stored on the
internal scale, against the recipient's response mapping, and never convert them", such an install being
clean. So this plan
assumes, when Part B starts, a [`Sampler::setState`](../../src/bartcore/sampler.hpp) with no units pass
(today's retired: [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp) call and its `unitsRefused`
refusal, read on 7f948d02), no `valuesMoved` feeding `altered`, and k installed as stored. Read: the
k-internal plan at 13e4818b predates dec-B418 (its step 7 still divides k by the units ratio); its
revision for dec-B418 is the orchestrator's, and Part B checks the landed form before it starts (Budget
and stops). dec-B416 (the maintainer: "Shouldn't `k` be `k`, the data scale be fixed, and `sd` syntactic
sugar?") fixes k.scale by the data, so a literal k (dec-B401) keeps the spread across a change of prior
too; dec-B417 (a drawn sd lags at a re-anchor) touches nothing here.

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
- A direction on a column this sampler has never seen missing, in a state from another sampler, is
  dropped and the install is not clean (B2 case c): a tree that would have to change, dec-B305's entry.
- The unforced TRUE or FALSE is returned visibly, as setPredictor's unforced value is (dec-B305 asked for
  setPredictor's idiom).

Out, as the design's critique had it: one error text per cause; a new floor on leaf values and sigma (TODO
precision-floor-edge-values); the grid half. An engine `force` argument defaults so tests/cpp call sites
stay; the verdict leaves a dead pointer dead; the help leads with the checked form; the verdict and the
install share one routine (B2).

Code read on 7f948d02:
- [`Sampler::setState`](../../src/bartcore/sampler.hpp) validates every chain before any is touched,
  except the state's cut grid, installed over the store before the containment and validity checks and
  put back by `restoreGrid` on a refusal (integers and vectors moved back, so exact). Its `altered`
  report is set by the units pass (gone with the k-internal slice) and by
  [`rebuildLiveForest`](../../src/bartcore/chain.hpp) and
  [`rebuildVarianceForest`](../../src/bartcore/chain.hpp): a bottom no row reaches, a split outside its
  interval ([`holdsSplitOutsideInterval`](../../src/bartcore/tree.hpp)), or a direction dropped by
  [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp). It refuses a live split outside a forest's
  column mask ([`columnMaskStateFeasible`](../../src/bartcore/chain.hpp)), a live tree breaking an
  interaction limit ([`interactionStateFeasible`](../../src/bartcore/chain.hpp)) and a monotone tree
  out of its cone (retired: [`monotoneStateFeasible`](../../src/bartcore/chain.hpp), folded by Part B into
  [`Chain::stateIsValid`](../../src/bartcore/chain.hpp)'s verdict), each by its own message
  from [`bartcore_setState`](../../src/R_interface_bartcore.cpp)'s helper.
- [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) already builds every live tree on a scratch
  `Tree` but never partitions it; it refuses a saved block whose size is not capacity times tree count.
- [`Chain::setState`](../../src/bartcore/chain.hpp) installs the response's latents before the trees, so
  a merge weighs leaves by the state's working weights; it writes `forest.k = fs.k` where the forest
  draws k.
- The merge a repair needs exists: [`collapseEmptyNodes`](../../src/bartcore/tree.hpp) walks the tree in
  pre-order and, at the first node meeting its predicate on a path (an empty child, an unrepresentable
  rule, a split outside its interval), replaces the subtree by one leaf at the weighted mean of the
  leaves beneath (geometric for the variance forest). [`reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp)
  sets an out-of-cone tree's leaves to 0; `rebuildLiveForest` already calls it after a merge, and the
  warm start calls it on an unconstrained donor.
- `ColumnStore::hasMissing` ([data.hpp](../../src/bartcore/data.hpp)) is rewritten from content on every
  quantize ([`quantizeDenseObserved`](../../src/bartcore/data.hpp),
  [`quantizeDenseCodesInto`](../../src/bartcore/data.hpp),
  [`quantizeCscColumnInto`](../../src/bartcore/data.hpp)); [`setCell`](../../src/bartcore/data.hpp) only
  raises it. A raised flag makes the birth and change moves draw a direction at one half and halves a
  rule's prior on the column, one rule per direction;
  retired: [`dropStaleMissingDirections`](../../src/bartcore/chain.hpp) (removed by Part A) clears directions on a column whose flag
  is down, from the data-mutation paths ([`applyNewData`](../../src/bartcore/chain.hpp),
  `forceRefreshTrees`, `rebuildFitsFromParameters`, and the variance forest's
  `rebuildVarianceFactors` and `refreshVarianceForest`) and from `buildFromFlat`.
- R: [`installStateOnto`](../../R/dbarts.R) is the install shared by
  [`dbartsSampler$setState`](../../R/dbarts.R), [`dbartsSampler$getPointer`](../../R/dbarts.R) and
  [`dbartsSampler$copy`](../../R/dbarts.R); after the engine call it runs reapplyForestWeights,
  reapplyActiveRows and retired: [`reissueNamedLeafSd`](../../R/dbarts.R), the last of which the k-internal
  slice removes. The test-data refusal [`unroutableTestColumns`](../../R/data.R) reads whether the
  TRAINING matrix holds a missing value now.

Run on 17df64ff (scratch/isplan/01-behaviours.R and 02-restore-cost.R, library scratch/libs/isplan; the
critique reran both, rerun-01 identical, rerun-02 within a few percent): setState returns TRUE invisibly
on its own state; a state of 10 kept draws into a store of 4, 4 into 10 and 4 into none are refused
("state is not consistent with this sampler") and none into 4 installs, TRUE; a state stored while a
numeric column held missing values, restored after a forced setPredictor filled them, gives FALSE; a
first missing value by an unforced setPredictor is taken (TRUE) and predict then answers for it. (Two
behaviours across response units run then are the k-internal slice's under dec-B418.)

## Part A: missingness first seen (dec-B320, dec-B321, dec-B322, dec-B378, dec-B399, dec-B400)

Lands first: once a column's flag stays raised, a direction is never dropped on the sampler's own trees,
and Part B's verdict has one direction case left (B2, case c).

A1. The flag is sticky. Every train-side quantize ORs into `hasMissing[j]` instead of assigning it, setData
included (dec-B399): setData ORs the sampler's flags into the new data's content when it rebuilds the
store.
Creation sets the flag from content (as today) or from the record (A3). The rollback snapshots of
`WholeMatrixUpdate` and `SubsetUpdate` ([sampler.hpp](../../src/bartcore/sampler.hpp)) already put the
flag back on a refusal. The test store tracks no flag.

A2. The first missing value draws. Where an accepted predictor change or a setData raises a column's
flag from 0, each chain draws from its own generator, with `ext_rng_simulateBernoulli(rng, 0.5)` as the
birth draw takes it, a direction for every split already on that column: first on its live trees
(forests in order, trees in order, nodes in pre-order, then the variance forest's trees), then on its
kept draws (dec-B400), oldest draw first in [`savedSlotForDraw`](../../src/bartcore/sampler.hpp)'s
order, each draw's mean forests then its variance trees, setting the flat flag; a pooled column (more
than 63 levels) keeps its bit in the pool words, as built (dec-A165), live and kept alike. A refusal puts
back every bit drawn, live and kept, and each chain's generator (serialized before the draw, as
[`Chain::setState`](../../src/bartcore/chain.hpp) reads one back; each chain owns its generator, so the
round trip is exact). Sites:
- [`runPredictorTransaction`](../../src/bartcore/sampler.hpp), forced: after `applyForced`, before
  `forceRefreshTrees`.
- unforced: after `snapshotApply` and before [`revalidateAllChains`](../../src/bartcore/sampler.hpp),
  since occupancy and the monotone order are judged under the drawn directions.
- the per-observation session [`UpdateSessionImpl`](../../src/bartcore/sampler.hpp): in
  `observationWouldRemainValid`, the first missing row raises the flag provisionally and draws, live and
  kept, then judges the row and the monotone order under the draws (replacing `orderHoldsWithMissing`'s
  temporary raise). The draw stays pending until `commitObservation`; a row not committed, because this
  session or, in [`updatePredictorPerObservationJointly`](../../src/bartcore/facade.hpp), another
  sampler's declined it, has the flag lowered and the bits and generators put back before the next row
  is judged or at `finalize`, so the next missing row draws the same directions.
- setData (dec-B399): in [`applyNewData`](../../src/bartcore/chain.hpp), the per-chain phase that runs
  once the store holds the new data, where its `dropStaleMissingDirections` call sits today (read): before
  the trees are remapped and collapsed, so a collapse merges under the drawn directions.
A column already flagged draws nothing.

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
  [`bartcoreSamplerSetPredictor`](../../R/bartcore.R), the joint row update) and after setData, from a new
  bridge reader of the engine's flags. setData does not read the slot of the object it is handed: the
  sampler's flags are the record's source, and the slot is overwritten after the call, so
  `d <- s$data; d@y <- y2; s$setData(d)` leaves a record equal to the engine's.
- Read by the bridge at creation ([`bartcore_create`](../../src/R_interface_bartcore.cpp) and the handle
  path), ORed into the content flags; and by the predict refusal (A5).
- A fit object changes: a fitted sampler's data carries the slot, so the exact gates run in quick mode at
  Part A's landing (CLAUDE.local.md's rule for a change to what a fit object carries), and the slot is
  described in [dbartsData.Rd](../../man/dbartsData.Rd)'s Value paragraph beside `rowNames`. A fit saved
  before has no slot and reads as none.

A4. No direction is dropped afterwards, and dec-B378. With A1 every train-side path keeps the flag, so
`dropStaleMissingDirections` clears nothing on the sampler's own trees: its data-mutation sites go, and
`buildFromFlat`'s stays for a flat tree from another sampler. In [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp)
the gauge on a column whose flag is down keeps counting the missing position reachable
(`missingReaches`), but a categorical rule that, with its missing bit dropped, sends every reachable
level one way is refused as malformed (dec-B378), where today it is built and merged; on a flagged column
it is built, and merged where no row is missing now (not clean, B2). Under dec-B399 a column setData
filled stays flagged, so the critique's case (one-way rules stored before a setData to complete data) is
built and merged, not refused; dec-B378's refusal is reached only by a state from another sampler. A
numeric or inline categorical rule with a direction on a column whose flag is down loses it, reported as
today through `directionDropped`, and the install is not clean (Settled).

A5. predict (dec-B322). [`unroutableTestColumns`](../../R/data.R) and
[`refuseTestMissingness`](../../R/data.R) take the record beside the training matrix (their callers hold
the data object), so a column made missable by setPredictor or setData, then filled, is still routable; a
column never missing is refused as built. The engine's own check (store.hasMissing at the predict
entries) follows A1 with no change.

A6. Help and docs. The help wording dec-B400 requires ("we should make sure what happens is clear in the
documentation"; the rule: "the help says so where it describes how a column becomes missable and what
predict does with such a row"), as written for the implementer:
- [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd), a paragraph "Missing values in
  predictors" in Details: "A predictor column can hold missing values from the first time the sampler is
  given one in it - at creation, or later by \code{setPredictor}, \code{setData} or
  \code{\link{updatePredictorPerObservationJointly}} - and from then on for the life of the sampler,
  whatever values the column holds afterwards. Each split rule on such a column sends missing values down
  one side, and the sampler draws that side as part of the model. When a column first becomes able to
  hold missing values, every rule already on it, in the current trees and in every saved draw
  (\code{keepTrees}), has its side drawn at random, left or right with probability one half each, from
  each chain's own generator; a change that is refused undoes these draws with the rest of it. A saved
  draw was made when no value in the column was missing, so nothing it was drawn from bears on that side,
  and one half is that draw's posterior for it. \code{predict} then routes a missing value in the column
  through every saved draw, those made before its first missing value and those made after. A column
  that has never held a missing value has no side to route by, and \code{predict} refuses a missing value
  in it. \code{copy} and a reload keep which columns can hold missing values."
- [na.keepPredictors.Rd](../../man/na.keepPredictors.Rd), At Prediction: "routable when its column had
  missing values in training" becomes "routable when its column can hold missing values: it had them in
  training or, on a sampler, has been given them since (see \code{\link{dbartsSampler}})"; "a column that
  was complete in training" becomes "a column that has never held one".
- The setPredictor and setData items point to the paragraph; [dbartsData.Rd](../../man/dbartsData.Rd)
  the slot (A3).
- [mia-missingness.md](../design/mia-missingness.md) takes an amendment under its Status paragraph, in
  [Representation](../design/mia-missingness.md#representation) ("cleared when a column loses its NAs")
  and in [Bridge and R surface](../design/mia-missingness.md#bridge-and-r-surface) ("a column that gains
  NAs mid-run routes them by ... (left)"), both becoming the first-seen rule with kept draws drawing.

## Part B: setState's two forms (dec-B305, dec-B310, dec-B318, dec-B398, dec-B401)

B1. The flag. `Sampler::setState` (and the facade virtual on
[`SamplerBase`](../../src/bartcore/facade.hpp)) take, after `adoptCapacity`, `bool force = true` and
`bool* notClean = nullptr`, so the tests/cpp call sites (87 lines matching `setState(`, chain-level calls
among them; ran: grep on 7f948d02) keep today's behaviour. The `columnMaskRefused`, `interactionRefused`
and `monotoneRefused` out-parameters go: those cases are verdicts now (B2, B8). `--preclean` on every
install: a virtual changes.

B2. Clean, in code. After every refusal (unchanged) and with only the state's grid written, an unforced
call asks each chain for a verdict and, if any chain is not clean, puts the grid back, sets `*notClean`
and returns false, the chains, store, generators, flags and latents untouched (and, in R, the mirrors,
B6). A chain is not clean when any live mean tree or variance tree, built from its flat form and
partitioned over the sampler's rows:
- (a) has a bottom node no row reaches;
- (b) holds a split outside the interval its ancestors leave;
- (c) carries a missing direction on a non-pooled column this sampler has never seen missing;
- (d) splits on a column its forest's mask forbids (a moderator subset, a blocks() row, a restricted
  variance forest);
- (e) breaks its forest's interactions() limit;
- (f) is a monotone tree whose leaf values, as stored, leave the cone.

Everything else that differs from the state is clean: an install across a change of response mapping,
nothing converted (dec-B418: "any parameter values are a position the chain could hold, so such an
install is clean", the register's words); latents redrawn under other weights or censoring
(`reapplyWeights`, `reapplySurvivalStatus`, after the install as today); a store of another size (B4); a
k recorded under another prior, installed literally (B5); DART weights, a fixed value or a generator of
another kind left as the sampler's.

One routine: a new Tree-level staging call (build from flat, partition, report a to e) used by the
verdict on a scratch tree with one n-length index buffer, as `stateIsValid` builds today, and by
`rebuildLiveForest` and `rebuildVarianceForest` on the live tree, so the verdict and the install's report
cannot disagree; f is the monotone check both already run. A tests/cpp check asserts it over the fuzz
harness's states (unforced false exactly when the forced install reports a repair). Forced, the verdict is
skipped and the install repairs (B8).

B3. Its cost. Measured on 17df64ff (ran: 02-restore-cost.R, one chain, one thread, 200 trees, 10 columns,
machine under load 5 to 6; rerun by the critique within a few percent): a clean restore of a fitted state
costs 74.0 ms against a 41.0 ms sweep at n = 1e5 (1.81 sweeps) and 3.3 against 2.1 ms at n = 5000 (1.53);
restoring single-leaf trees costs 19.2 and 1.0 ms, so the per-tree build and partition is 74 and 70
percent of a restore, 1.34 and 1.07 sweeps. The scratch verdict repeats that work less the fits, so an
unforced restore costs up to about 1.1 to 1.3 sweeps more (bound, not run; B7 measures it). Cases d to f
add nothing measurable: the containment, interaction and monotone scratch builds run today on every
restore.

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
  (`savedSlotForDraw`'s rule), the newest K = min(R, C) draws go to slots 0 to K - 1, oldest first, with
  their params, masks and variance trees; `recordedDraws_` = K, `currentSampleNum_` = K mod C; n.samples
  and C do not change. C of 0 takes none.
- A re-creation (copy, reload) still adopts S ([`resizeSavedTrees`](../../src/bartcore/sampler.hpp)),
  as dec-B294 and dec-A188 have it.
`lengthscaleStateFeasible` is asked about saved gp draws only where K > 0.

B5. k (dec-B401). A state's k installs as recorded whatever prior the recipient holds, as
`Chain::setState` does today (`forest.k = fs.k` where the forest draws k, read); a recipient holding k
fixed keeps its own. Nothing is built. Under dec-B416 k.scale is the data's whatever the prior, so the
leaf sd k.scale / k is kept across a change of prior as well; the k-internal slice owns the test of that
row. dec-B401 asks that "the help says the sd then follows the recipient's k.scale": one sentence in
Saving (B9).

B6. Bridge and R.
- [`bartcore_setState`](../../src/R_interface_bartcore.cpp) takes `force` (arity 4 to 5 in
  R_interface.cpp) and returns TRUE when it installed, FALSE when unforced and not clean; every remaining
  refusal stays an error, and the three constraint messages (the column restriction, the interaction
  constraint, the monotone constraint) go. reapplyWeights and reapplySurvivalStatus run only after an
  install.
- [`installStateOnto`](../../R/dbarts.R) takes `force`; on FALSE it returns FALSE before
  reapplyForestWeights and reapplyActiveRows.
- `setState(newState, forceUpdate = FALSE)`: `forceUpdate` a single TRUE or FALSE, else "'forceUpdate'
  must be TRUE or FALSE" ("partial" included); on FALSE neither the pointer (a dead one stays dead) nor
  `$state` is assigned and FALSE is returned, visibly; on a clean unforced install TRUE, visibly; forced,
  NULL invisibly (dec-B310).
- getPointer and copy pass `force = TRUE`, silent (dec-B305, dec-B234), so a copy or reload of a state
  breaking a constraint is repaired where today it errors at every use.

B7. Speed, recorded. bench-sampler.R gains a restore scenario (store, forced setPredictor, run, put the
predictors back, unforced setState, and the same forced) at n = 1000 and 1e5, and 02-restore-cost.R is
rerun against the base build on one quiet machine. The unforced restore's ratio to the base restore and
the forced restore's (expected 1.00) go in the landing note; neither gates the landing. storeState and
run() do not change in this slice.

B8. Repair under force (dec-B398). Forced, cases d to f are installed as a warm start would take them:
- (d) and (e): the first offending split on each path from the root becomes one leaf with everything
  beneath it, at the row-weighted mean of the leaves beneath; built as a predicate on
  `collapseEmptyNodes`'s pre-order walk (a split on a masked column; a split whose variable, with its
  ancestors', breaks the forest's limit), so the merge, its weights and the variance forest's geometric
  merge are the routine already in use. The routine is warm-start-salvage's, pulled into this part; that
  item reuses it for installTrees's force.
- (f): the tree is reseeded, every leaf 0 ([`reseedInfeasibleMonotoneLeaves`](../../src/bartcore/chain.hpp)),
  as a warm start does today; `monotoneStateFeasible` stops being a refusal and feeds the verdict.
- A collapse can leave a monotone tree out of the cone; the reseed after it runs as built.
- The containment backstops in `rebuildLiveForest` (`interactionSubtreeIsValid`,
  `columnMaskSubtreeIsValid`) stay, now after the collapse, as the invariant's check.
- Saved draws are copied as stored: containment holds only for live trees, as built.

B9. Help. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd)'s Saving subsection is rewritten
led by the checked form (Docs and NEWS), and says: a state breaking the recipient's column restriction,
interaction limit or monotone cone returns FALSE, and forced is collapsed at the first offending split or
reseeded, where it was refused; a copy and a reload repair it silently; and (dec-B401) "a state's k is
installed as recorded, whatever leaf prior this sampler holds; the leaf sd is then k.scale / k against
this sampler's k.scale". What a restore across a change of response mapping means is the k-internal
slice's sentence (dec-B418), kept.

## Refusals and messages

Unchanged, in both forms, the sampler untouched: not a bartcoreState, an older encoding, another chain
count, another tree, forest or variance-forest count, another leaf model, a malformed block or tree, a
repeated grid point (dec-B300), latents out of range (dec-A174's floor), saved gp draws under other
lengthscales. Messages as built, "state is not consistent with this sampler" among them. New: "'forceUpdate'
must be TRUE or FALSE"; a one-way categorical rule on a column never missing (A4) reads as the other
malformed trees. No longer refused: a kept store of another size; a split on a forbidden column, an
interaction break and a monotone tree out of the cone (B8). Gone with the k-internal slice: units that
cannot be converted (dec-B418).

## Callers

Package (ran: grep on 7f948d02): copy and getPointer, forced; no other R caller. bart's warm.start goes
through installTrees, unchanged here. benchmarks/R: sbc.R (calibration) and negbin-mixing.R, an exact
gate's script, each hand a sampler its own state with a scalar edited; surfaces/C1-frozen-ess.R
transplants a recorded chain's state into a fresh sampler built on the same data with levelGibbs. Each is
clean (same rows, same mapping) and returns TRUE, its draws unchanged; none reads the value.

stan4bart (bartcore 963956b, ran: grep): one call, `sampler$setState(state)` in restoreBartSampler in
R/generics.R, reached from getBartSampler after a reload and from stan4bart_fit.R at fit end, value
dropped. dec-B402, the maintainer: "Force." It becomes `sampler$setState(state, forceUpdate = TRUE)`,
one line on stan4bart's bartcore branch, landing in lockstep with Part B (before it, the argument does
not exist); the orchestrator pushes it. Its suite runs against the build.

bartCause (dbarts-1.0 e833be7), treatSens (dbarts-1.0 babfaa6) and bairrtt (main 3f57f61): no setState,
copy, storeState or installTrees (ran: grep). bairrtt calls setPredictor forced and the joint row update
on a latent trait with no missing values (read), so Part A reaches its code path and not its draws.

## Tests

tinytest, new:
- test-setstate-force-update.R. Clean, each unforced TRUE (visible) and forced NULL invisibly, then the
  two twins' next 20 draws identical: the sampler's own state; a recorded chain transplanted into a fresh
  sampler; one stored before a predictor change that empties no leaf; other weights (Student-t) and
  censoring (aft); stores 10 into 4, 4 into 10, 4 into none, none into 4; a filled column's directions.
  Not clean, each unforced FALSE with getTrees, sigma, k, kept-draw predictions, cut points, `$state` and
  the next 20 draws identical to a twin that made no call, then forced NULL with the repaired install
  (for the leaf no row reaches, the trees a forced setPredictor gives): a leaf no row reaches (mean, and
  variance forest); a split outside its interval; a direction on a column never missing, from another
  sampler's state. A dead pointer stays dead after a refused call. `forceUpdate` NA, "partial" and
  c(TRUE, FALSE) refused.
- Repair (B8): a split on a forbidden column (a moderator forest, a blocks() row, a restricted variance
  forest), an interaction break and an out-of-cone monotone tree, each unforced FALSE with the sampler
  untouched and forced NULL with the collapse at the first offending split from the root (an allowed
  split above it kept, everything beneath gone), the collapsed leaf at the row-weighted mean, the
  variance forest's at the geometric mean, the monotone tree at 0; copy and reload of such a state repair
  silently.
- Store: 10 into 4 predicts as the source's last 4 draws; a wrapped source ring keeps its newest in
  order; 4 into 10 then 3 sweeps holds 7; equal capacity, a copy run 5 sweeps and its source run 5 store
  the same state.
- test-missingness-first-seen.R: forced, unforced, column, per-observation and setData updates raise the
  flag and draw; over 40 seeds the share of directions drawn right is within a binomial 1e-3 band of one
  half, in the live trees and in the kept draws recorded before (the critique's p3 shape: 20 kept draws,
  a first missing value, 5 more; today 0 of the older draws' splits send right); predict on a missing row
  averages all 20 draws; a second missing-bearing update draws nothing; filling the column (setPredictor
  or setData) keeps every direction and predict accepts a missing value there; setData to complete data,
  then setState of a state stored before it, merges the one-way rules (the critique's p1) and copy and
  reload of that sampler continue it bit for bit; a column never missing still refuses (dec-B322); a data
  object saved before the slot existed reads as no record; setData handed an object with a stale slot
  leaves the engine's record; dec-B378 refused on a never-missing column from another sampler's state,
  built and merged on a flagged one.
- A refused first missing value (an emptied leaf) leaves each chain's generator state and every kept
  draw identical. On the column and per-observation paths the next 20 draws are identical to an
  untouched twin's; on the whole-matrix path they are equal within 1e-12, the rounding TODO
  refused-update-rounding records there (about 1.8e-15, a rollback's repartition), which this slice does
  not own.

tinytest, changed (the design's prototype failed the first five, 130 expectations, ran there):
[test-state-empty-leaf-merge.R](../../inst/tinytest/test-state-empty-leaf-merge.R),
[test-state-missing-direction.R](../../inst/tinytest/test-state-missing-direction.R) (TRUE where no leaf
merges), [test-monotone-unforced.R](../../inst/tinytest/test-monotone-unforced.R),
[test-heteroscedastic-mutation.R](../../inst/tinytest/test-heteroscedastic-mutation.R),
[test-state-not-model.R](../../inst/tinytest/test-state-not-model.R), and those that assert a merged or
split-moved install's value: [test-cut-grid-distinct.R](../../inst/tinytest/test-cut-grid-distinct.R),
[test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R). Each merged install takes
`forceUpdate = TRUE` and gains its unforced FALSE. The setState refusal pins that become FALSE and
repaired (ran: grep on 7f948d02): ["allowed column set"](../../inst/tinytest/test-single-forest-vars.R)
(the warm-start pin beside it stays), retired: ["violates this sampler's interaction constraint"](../../inst/tinytest/test-interactions.R)
(two), retired: ["leaf values violate this sampler's monotone constraint"](../../inst/tinytest/test-monotone.R),
both messages gone with Part B.
test-blocks.R pins no setState refusal (ran: grep). The full suite finds the rest.

tests/cpp:
- the sticky flag (setData included) and the first-sight draw, live and kept, in
  [test_data.cpp](../../tests/cpp/test_data.cpp) and [test_sampler.cpp](../../tests/cpp/test_sampler.cpp);
  no data mutation lowers a flag; the pins that read the flag going down change: "drop every missing
  value: hasMissing flips false" in [test_moves.cpp](../../tests/cpp/test_moves.cpp) and
  [`testPerObservationMissingCommit`](../../tests/cpp/test_moves.cpp)'s neighbours;
- in [test_state.cpp](../../tests/cpp/test_state.cpp): the unforced refusal leaves `getState` bytes
  and the generator equal; verdict and install agree over the fuzz harness's states, constraint cases
  included; the store repack; the collapse at the first offending split.

Mutants the reviewer runs, each failing a test: the flag assigned from content, or reset by setData; no
first-sight draw, one on every update, or none on the kept draws or through setData; the generator or a
kept draw's bit not put back on a refused first missing value; the record not written, not read at
creation, or read by setData from its argument; dec-B378's rule built and merged on a never-missing
column; install then report FALSE; the verdict blind to the variance forest, or to case d, e or f; a
forced constraint break still refused, or collapsed below the first offending split; case c clean; forced
TRUE or visible; the pointer or `$state` bound on FALSE; the oldest draws kept, or the ring repacked at
equal capacity.

## Equivalence and snapshots

Expected: nothing recorded moves. No equivalence scenario calls setState (ran: grep; bcf-equivalence.R
and multinomial-equivalence.R say a restore is deliberately not recorded), and the missing-value
scenarios (`missing`, `nafactor`, `testswap` in equivalence.R, the only places any of the three
harnesses writes NA, ran: grep) put their missing values in at creation and mutate no training predictor
(read); the mutation and setData scenarios use complete data (the critique, read), so a sticky setData
flag changes nothing there. The four reproducibility files call neither setState nor setPredictor (ran:
grep on 7f948d02). The exact gates' one caller, negbin-mixing.R, is a clean own-state install. So both
parts run against the MANIFEST's current files with every scenario "identical draws (same RNG stream)"
and no |z| line; any move is a stop.

## Docs and NEWS

- [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the [`setState`](../../man/dbartsSampler-class.Rd)
  usage, the newState and forceUpdate items (forceUpdate is shared with setPredictor), the Saving
  subsection led by the checked form (`if (!sampler$setState(st)) ...`) and naming `forceUpdate = TRUE`
  for a loop that restores its own state onto rows it has put back, TRUE stated first as not bit for bit
  (trees and leaf values as stored, drawn scalars, the newest kept draws that fit; latents may be
  redrawn), the repair (B9), copy and reload forced and silent, and Value; Part A's paragraph (A6). The
  docstrings of setState and copy in [dbarts.R](../../R/dbarts.R) (rc-codoc).
  [dbartsData.Rd](../../man/dbartsData.Rd): the slot (A3). [na.keepPredictors.Rd](../../man/na.keepPredictors.Rd)
  (A6).
- [NEWS.Rd](../../inst/NEWS.Rd), changes from 0.9-34 only: setState gains forceUpdate, unforced by
  default, and returns TRUE or FALSE (0.9-34 installed and returned NULL); a state with more or fewer
  kept draws than the sampler keeps installs what fits (0.9-34 took a larger or smaller store and crashed
  R on a state with none, dec-B318); a first missing value given to a column has its splits' directions
  drawn, in the kept draws too (0.9-34 sent them down an arbitrary branch).
- [mia-missingness.md](../design/mia-missingness.md) (A6). A short docs/design/install-surface.md, the
  rules this plan builds and the held calls as ruled (dec-B398 to dec-B402), with its INDEX row; the
  scratch design stays scratch.
- TODO at landing: setstate-force-update, state-store-size-install and missingness-first-seen close;
  forced-update-returns-null keeps installTrees's half; warm-start-salvage notes that its collapse
  routine landed here; a note on the verdict skip (B3). state-install-keeps-spread is the k-internal
  slice's to close (dec-B418 reverses dec-B384).

## Gates

Per part, on its tip and its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting):
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`); stan4bart's (with its one-line change, Part B) and
  bairrtt's suites against the build.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) and the four arm64 snapshot files.
- Every exact gate in quick mode (on Part A also because a fit object gains a slot, A3); `R CMD check
  --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.
- Part A: bench-sampler.R compare (the setPredictor scenarios); Part B: B7, recorded.

## Budget and stops

Planned, from the critique-revised plan's ~730 and ~850, with the rulings' branches folded in (lines
planned):

| part | before | added | removed | planned |
|---|---|---|---|---|
| A | 730 | dec-B399 sticky setData +120 (engine 40, R 15, tests 60, help 5); dec-B400 kept draws +100 and its help +20 (engine 50, tests 50, help 20) | none | ~970 |
| B | 850 | dec-B398 repair +350 (engine 180, bridge and R 10, tests 150, help 10); dec-B401 help +5 | B5's re-expression, `kFollowsUnits`, its bridge and R argument and its tests, the units-clean case and its test, the state-not-model.md k row: about -100 (dec-B418, the k-internal slice) | ~1100 |

Planned against landed on this surface ran 1.1 to 2.9 times (the first critique's table of six plans);
the last engine slice measured ran 1.76 times ([small-rulings-1008-engine.md](small-rulings-1008-engine.md),
950 to 1675). Forecast at 1.5 to 2 times: Part A 1460 to 1940, Part B 1650 to 2200. The stops are 1.5
times the forecast's midpoint (1.75): Part A ~2550, Part B ~2900. stan4bart's one line is outside these
counts. Stop and report, without working around, when:
- a part passes its stop;
- any equivalence scenario or snapshot moves;
- Part B starts and the landed k-internal slice still converts a state's leaves or k on install, keeps
  the units refusal, or left `Sampler::setState`'s signature other than this plan assumes (Context): the
  plan is written against dec-B418's form;
- the verdict and the install disagree on a fuzz state and one routine cannot be made to serve both;
- the collapse at a forbidden or interaction-breaking split cannot be a predicate on
  `collapseEmptyNodes`'s walk and needs a second merge;
- the kept-draw first-sight draw cannot be undone with the rest on a refusal without copying the whole
  store;
- a reading needs a state-format floor change, or a rewritten test needs a call no ruling or settled
  call makes.

## Held for the maintainer: resolved

Each was put on 2026-10-09 with both branches costed; the branch taken is folded into the parts above and
its lines into Budget and stops.

H1. A constraint the state's tree breaks, under force (dec-B398). The maintainer: "Repair, as we do
elsewhere." Taken: repair. Not clean, cases d to f (B2); forced, the collapse at the first offending
split and the reseed (B8); the three refusals and their messages go (B6); the four pins become FALSE and
repaired, with new collapse tests (Tests); the collapse routine is warm-start-salvage's, built here. +350
planned, Part B.

H2. setData and missingness (dec-B399). The maintainer: "setData keeps the flags (sticky)." Taken:
sticky. setData ORs the flags (A1), draws at first sight live and kept (A2), and the record is written
from the engine, the argument's slot not read (A3); dec-B378's refusal after a setData does not arise
(A4). +120 planned, Part A.

H3. Kept draws recorded before a column's first missing value (dec-B400). The maintainer: "OK, you can
use your recommendation but we should make sure what happens is clear in the documentation." Taken: kept
draws draw too, at every site, undone with the rest on a refusal (A2); the help wording (A6); the share
in older draws within the binomial band and the refused case bitwise (Tests). +120 planned with the help,
Part A.

H4. The spread of a state stored under another prior (dec-B401). The maintainer: "The k itself should
literally transfer - a parameter is a parameter, regardless of the prior. It may be a bad fit, but so be
it." Taken: k literal, as built (B5); one help sentence (B9). Under dec-B416 k.scale no longer moves with
the prior, so the spread is kept too; its test is the k-internal slice's. +5 planned, Part B.

H5. stan4bart's restore (dec-B402). The maintainer: "Force." Taken: `forceUpdate = TRUE` in
restoreBartSampler, one line on stan4bart's bartcore branch, in lockstep with Part B (Callers); its suite
as is. Outside this plan's counts.

## Interactions

- The k-internal slice lands first (dec-A191) and removes the units pass, its refusal and
  `reissueNamedLeafSd` from the install; Part B's B1 and B6 are written against that (Budget and stops).
  Part A touches no leaf scale and does not depend on it.
- The width-weighted default cut rule lands before; it changes grids at creation and setData, not the
  install, so B2's partition reads whatever grid the state carries. Rebase only.
- gp-copy-continuation landed (dec-A192); B4 and B2 do not touch the gp draw.
- response-scale-rows follows this plan; its setData re-derivation meets A1's sticky OR only through
  setData's re-quantize. Rebase only.
- state-frame-prior's grid half later moves the grid out of the state; B2 then gains "a split off the
  grid", and the missing-value record sits beside the other derived records on dbartsData.
- warm-start-salvage builds installTrees's `forceUpdate`, its NULL (dec-B310) and reuses B8's collapse.

## Revision from critique

Against aa79ac10, from scratch/iscrit/critique.md and the orchestrator's calls:
- The verdict skip (hash, state attribute, its tests, mutants and stop) is dropped; B3 states the
  measured cost, the spare-buffer swap weighed with its 80 MB a chain, and `forceUpdate = TRUE` as the
  loop's way out; B7 records the cost and gates nothing.
- Stops are 1.5 times the forecast per part, not below it.
- The refused first-missing-value test is bitwise on the column and per-observation paths and within
  1e-12 on the whole-matrix path, citing TODO refused-update-rounding.
- The record is a dbartsData slot, `missing.seen`, read through a `.hasSlot` helper; the fit-object
  consequence (exact gates at Part A's landing, the Rd) is stated, and setData rewrites the record.
- Settled and moved out of Open calls: two landings with Part A first, `newState`, installTrees's NULL
  with warm-start-salvage, the slot; added: `forceUpdate` is the one name with installTrees.
- Five questions held for the maintainer with both branches (resolved since; see the next section).
- Fixed: C1-frozen-ess.R transplants another sampler's state; mia-missingness.md's "cleared when a column
  loses its NAs" is under Representation; the "first slice" scheduling is marked as the register
  author's.

## Revision for the 2026-10-09 rulings

Against dd0872da, rebased onto 7f948d02:
- H1 to H5 resolved as ruled (dec-B398 to dec-B402), each branch folded into the parts, tests, mutants
  and budget, the section rewritten as resolved with the maintainer's words.
- dec-B418 (with dec-B416 and dec-A191): B5's re-expression of k, `kFollowsUnits`, the bridge's sixth
  argument, the units-clean case, the spread tests and mutants, the "units that cannot be converted"
  refusal and the state-not-model.md k row are gone, the k-internal slice removing what is built; "clean"
  now counts an install across a change of mapping, nothing converted; Part B is written against that
  slice and stops if it landed otherwise. TODO state-install-keeps-spread leaves this plan's closing list.
- The help wording dec-B400 requires is written out (A6), na.keepPredictors's At Prediction included.
- The orchestrator's two calls after the first revision (case c dropped and not clean; the unforced value
  visible) moved from Open calls to Settled; Open calls is gone.
- Re-sized: planned ~970 and ~1100, forecast at 1.5 to 2 times, stops at ~2550 and ~2900, with three new
  stops (the k-internal slice's form, the collapse as a predicate, the kept-draw undo).
- Re-read on 7f948d02: the cited engine, bridge and R symbols; the dropStaleMissingDirections sites by
  name; the existing merge and reseed a repair reuses; tests/cpp `setState(` lines 87 (88 on 17df64ff);
  test-blocks.R pins no setState refusal (the earlier H1 named it); the callers and snapshot greps.

## Evidence

Ran on 17df64ff (private library scratch/libs/isplan; probes and outputs in scratch/isplan/): the
behaviours in Context (01-behaviours.R); the restore costs in B3 (02-restore-cost.R, under load). Ran on
7f948d02 (grep): callers, tests/cpp setState lines, the refusal pins, the snapshot files. From the
critique's runs (scratch/iscrit/): the store and hash costs, the one-way rules after setData (p1), the
kept draws across a first missing value (p3), attribute loss on subsetting. Read on 7f948d02: every engine,
bridge and R claim cited by symbol; the TODO entries; dec-B283 to dec-B418 and dec-A160 to dec-A192;
stan4bart's restoreBartSampler at 963956b; the k-internal plan at 13e4818b. Bounds not run: the scratch
verdict's added cost. The help wording in A6 is a draft for the implementer, not run through R CMD check.

## Landing note: Part A

Landed 2026-10-10 as acf6a6fd..2269918e after one opus review that ran the plan's mutants and 19 of its
own, one fix round (tests and wording) the same reviewer checked, and an x86 run of tests/cpp plain and
under the sanitizers and of tinytest on the fix round's tip. What was built is Part A, with these
departures and additions:

- A rule on an unordered factor draws its direction only where a missing value can reach it; a rule out of
  reach holds none, as on a sampler that saw the missing value at creation. Drawing everywhere breaks the
  sampler's reinstall of its own state. The help and the design note say so (dec-A198).
- On the per-observation paths the first-sight draw is made for the row being judged and taken back when
  that row is declined, which leaves rules, kept draws and generators as one draw at the first accepted
  row leaves them.
- The record of columns seen missing is the data object's `missing.seen`, NULL until a predictor change
  writes it; creation reads it, and so does the handle path that xbart's folds use.
- A factor column's first missing value by setPredictor stays refused at the R surface, the refusal lifted
  once the record holds the column; `dropStaleMissingDirections` is removed.
- The review found no wrong behavior and four untested ones (two forests, a declined row, an uncommitted
  draw at the end of a joint sweep, the handle path); the fix round's tests took the part to about 2,000
  lines outside this plan.

Open after this part, in the root TODO: a state stored before a column's first missing value and
installed after it (state-before-first-missing); the factor refusal against dec-B321's words
(factor-first-missing-setpredictor); three tests that hold by one seed (first-missing-test-seeds).
