# state-before-first-missing: a state records which columns could be missing, and an install draws what it lacks

Status: LANDED 2026-10-11 as a66b8f65..d180dcc2 (dec-B435, dec-B436, dec-B440, dec-A201). See the
[Landing note](#landing-note).

agent: opus implementer, one; one opus reviewer who runs the mutants below.
rng: SHIFTING on four paths: an install whose state's record lacks a column the sampler has flagged (coins from
the state's generator); a warm start whose donor's record lacks one (coins from each receiving chain's
generator, before its first sweep); a factor column's first missing value by an update refused today; a
first-sight draw over a kept rule of section 2's exempt form (one coin fewer). NEUTRAL, bit for bit, otherwise.
window: before the merge to main; engine slices stay serial (sparse-starting-sigma: section 8).
budget: planned ~1035 lines (step G ~105, A ~620, B ~185, C ~125); forecast and stops below.

## Goal

A stored state says which predictor columns could hold a missing value when it was stored. When `setState`, a
copy or a reload installs it into a sampler where a further column can hold one, the direction of every split on
that column in the state's live trees and kept draws is drawn at one half, before anything is judged (dec-B435).
A warm start draws the same for the trees it installs, from each receiving chain's generator (dec-B440). A
factor column takes its first missing value as a numeric column does (dec-B436). The tests the review of
setState's two forms left open are added. Tier: "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

dec-B435 (extends dec-B400), the maintainer on 2026-10-10: "Sure, draw it." The register's rule is the Goal's
first two sentences. dec-B436, the same day: "Lift the refusal." dec-B440: section 5. Rerun by me on the rebased
base 82678655, library scratch/libs/smr2: 01-review-case.R gives the ledger's numbers again (a forced
setPredictor's 3 missing values in x1 leave live trees 3 left and 5 right, kept draws 95 and 95; `setState` of
the older state then leaves 8 and 0, 190 and 0, and returns TRUE; a reload, a `copy()` and a twin created with
the values the same). Part B changed one thing: the TRUE is now visible.

Constraints: no change to the state-format version or floor, to [dbarts.h](../../inst/include/dbarts/dbarts.h)
or to stan4bart; nothing recorded moves (section 8); Part B's forms are used, not changed. Out of scope:
installTrees's `forceUpdate` (TODO forced-update-returns-null); the mirror of the exempt rule.

## 1. The record (step A)

- Engine: [`SamplerStateData`](../../src/bartcore/sampler.hpp) gains `missingColumns`, one byte per predictor,
  empty meaning no record; [`Sampler::getState`](../../src/bartcore/sampler.hpp) always fills it from the
  store's flags ([`hasMissing`](../../src/bartcore/data.hpp)), not the content. Of another length, `setState`
  returns false before anything is written.
- Stored object: a top-level attribute `missing.columns`, a logical per predictor with no NA, always written,
  all FALSE on a complete sampler (so not through [`bartcore_getMissingSeen`](../../src/R_interface_bartcore.cpp),
  NULL when no flag is up): absence means not known. The bridge's [`storeState`](../../src/R_interface_bartcore.cpp)
  alone writes it (no R code builds a state; ran: grep); read by name beside the digests, another type or an NA
  is an error in both forms, another length "state is not consistent with this sampler".
- A state WITHOUT the attribute (this unreleased branch; stan4bart fits saved on it) draws nothing and installs
  as landed: read as "no column could be missing" it would redraw directions a chain learned.
- Format: additive under [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)'s registry rule, version and
  floor staying 1; the flat C API has no state entry (ran: grep). The exact gates run quick: a fit's state grew.

## 2. The install rule (step A)

With S the state's record and F the sampler's flags, `raised[j] = F[j] && !S[j]` is drawn (dec-B435). A column
in both keeps its directions, clean (dec-B320); one in S and not in F is as landed, the record playing no part
(a direction dropped and not clean, a pooled rule's bit kept, a rule sending every level one way refused,
dec-B378); a state never raises a flag (dec-B321). Equal records, or none: the state installs as today.

The draw completes a copy of the state's chains before anything of the sampler is written: in
[`Sampler::setState`](../../src/bartcore/sampler.hpp) after the lengthscale check and before the grid snapshot
(section 6, steps 3 and 4). It reads column kinds and the state's own ring, never the grid, so a throw or a
later refusal has nothing to undo. With a column raised and the state's blocks naming one capacity (else step 5
refuses the state): each flat tree, live and kept, is checked with [`flatTreeIsWellFormed`](../../src/bartcore/tree.hpp)
and its mask channel's pairing (the draw recurses with no bounds; a failure leaves the refusal to step 5);
`state.chains` is copied; each chain draws on its part (a new `Chain::completeMissingDirections`, section 3)
with [`drawFlatMissingDirections`](../../src/bartcore/tree.hpp) in
[`Chain::drawMissingDirections`](../../src/bartcore/chain.hpp)'s order: live trees by forest, variance trees,
then the state's recorded kept draws oldest first by its own ring (capacity S, R = min(recordedDraws, S), slot
(currentSampleNum + S - R + d) mod S). Only steps 5 and 8 read the copy.
- The exempt rule (critique M1). A subset rule with no category on the right sends only the missing value right:
  its direction is its form. `drawFlatMissingDirections`, for both its callers (the first-sight draw's kept loop
  and this completion), leaves it as it is, in reach or out, takes no coin for it and reads its direction for
  the rules beneath; it is handed the mask channel, a pooled rule's words lying at its record's offset. A coin
  to the left would leave a rule sending nothing right, which the form check refuses. Live trees never flagged
  in the column cannot hold the form ([`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) asks that the
  categories alone be split; read), so `Tree::drawMissingDirections` needs none.
- The landed defect it closes (TODO first-sight-kept-one-way-rule), rerun by me on 82678655 (the critique's p2,
  p2b, p2c): with the other sampler's state FORCED in, a first missing value by setData flipped such kept rules
  in 11 of 11 seeds holding one (39 tried) and the sampler could not install its own state (4 of 4 controls
  could); unforced, as the probes were written, Part B declines that state and 11 of 11 reinstall. It narrowed
  the way in and fixed nothing.
- Unforced: TRUE for an install that drew (Calls made), unless the completed trees fail the verdict on another
  ground (a leaf no row reaches, a monotone tree out of its cone): FALSE, on every repeat. Forced (a copy, a
  reload, stan4bart's restore): completed, then repaired as landed; NULL.
- Declined or refused: step 6 as landed, the copy dropped; each chain's generator went back inside its own draw
  (section 3). A throw (the copy's allocation) finds nothing written. Completed, only a hand-edited state that
  contradicts its record is refused, an error in both forms.

## 3. The generator (step A)

A state carries each chain's generator, which [`Chain::setState`](../../src/bartcore/chain.hpp) installs last
when its bytes are of the chain's length; the coins come from it. Per chain, inside a scope guard that reads the
chain's own bytes back on every exit: the chain's generator is serialized, the state's bytes read into it when
they fit, the coins drawn and the advanced bytes written to the copy's `rngState` (a Mersenne Twister through
the bridge, an exact round trip; read). A hand-built state without the bytes draws from the chain's own position
and the copy carries the advance: step 8 installs it and a decline leaves the chain where it was.
- One state installed twice into one sampler draws the same directions, so the cached `$state`, left the object
  `setState` was given (section 6), still describes the sampler: a copy or reload of it holds its directions.
- Reproducible, in this scope only: a store; ONE whole-column or whole-matrix setPredictor raising its columns
  in one call, no generator use between; `setState(stored)`: the install draws the update's coins and leaves the
  generator where the update left it. Rerun by me on 82678655, output identical (the critique's p1 replay in R:
  unforced, forced, two columns in one call). Not claimed, not a stop: the per-observation and joint paths,
  columns raised by separate calls, a sweep between, setData.
- Nothing to draw: the routine returns before the walkability check, the copy and any generator read.

## 4. A factor column's first missing value (step B, dec-B436)

Its own step and commit, after step A, whose exempt rule a factor's first-sight draw can meet; no engine change.
Rerun by me on 82678655 with outputs identical to the base's: 04-factor-first-missing.R (the refusal bypassed by
hand: every path draws, reinstalls its own state, predicts and copies) and 05-test-column-na.R.
- The refusal "has missing values, which its training values do not" and the `missing.seen` argument that fed it
  go from [`codeCategoricalColumnUpdate`](../../R/bartcore.R), `codePredictorFrame` and their callers: a missing
  label is coded as a missing value; an unknown label and a number stay refused.
- The data object's `missing.seen` for a factor is written after every accepted change, as for a numeric column
  ([`recordMissingSeen`](../../R/data.R); ran), and no training update reads it. Row by row and jointly the
  engine's session ([`UpdateSessionImpl`](../../src/bartcore/sampler.hpp)) draws at the first missing row and
  takes the draw back if the row is declined.
- setTestPredictor is not named by dec-B436 and keeps refusing a missing value in a column that has never held
  one (dec-B322). Its column form borrowed the refusal that goes, for factors only: a numeric test column so
  updated is taken today and the next run fits the row, where the whole test set and predict refuse (rerun). So
  it calls [`refuseTestMissingness`](../../R/data.R).

## 5. The warm start (step C, dec-B440)

The maintainer on 2026-10-10: "Sure, draw it." Rerun by me on 82678655, as before Part B (02-other-routes.R): 7
of 7 splits on x1 send left before any sweep, 8 and 4 after 31. Built after step A, its own commit.
- The record: the bridge's [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads `missing.columns`
  into the donor's `missingColumns`; with D that record and F the recipient's flags,
  `raised[j] = F[j] && !D[j]`. [`warmStartState`](../../R/dbarts.R) stores a donor sampler's or fit's state at
  the call, so it carries one; a raw state stored before the record has none and installs as today, all left,
  with no coin and no generator read, as when nothing is raised. Of another length: `shapeMismatch`.
- Where: in [`Sampler::installForests`](../../src/bartcore/sampler.hpp), on `install[c]`, the scratch it
  reassembles per receiving chain from the donor's live trees or the slot asked for, after
  [`convertDonorStandardization`](../../src/bartcore/chain.hpp) and before `checkContainment`, whose checks then
  judge the completed trees. The draw so precedes both collapses of emptied splits
  ([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), `rebuildLiveForestRemapped`) and the monotone
  reseed. It reads no grid, and the scratch is its copy.
- Order, source and undo: chain c draws for its own trees, each forest's in order and then the variance trees,
  with `drawFlatMissingDirections` (exempt rule included: a slot-sourced donor draw can hold one), from its own
  generator in place, so two chains seeded from one donor draw take their own sides. Each tree is first checked
  with `flatTreeIsWellFormed` (failing: `rebuildFailed`, before any coin). A warm start refused after the draw
  (containment, a rebuild, a throw) leaves nothing touched, as today: each chain's generator bytes are saved
  before its coins and read back by a guard dismissed only on `ok`.
- What shifts: only a warm start with a raised column and a split on it in the installed trees: each receiving
  chain's generator stands its coins further on before the first sweep, so that run's draws differ from today's
  build's. No stored draw, value, message or refusal changes. Ran (grep on 82678655): no equivalence scenario,
  reproducibility file or exact-gate script warm-starts from a donor; benchmarks' one, composition-matrix.R,
  builds both on the same complete data and records no draw.

## 6. The landed install this slice builds on

Read on 82678655. `Sampler::setState` runs: (1) shape and grid checks; (2)
[`Chain::savedStateCapacity`](../../src/bartcore/chain.hpp), one capacity asked across chains, mean forests and
the variance forest (a disagreement is refused at 5); (3) the lengthscale check, AFTER WHICH THE DRAW GOES; (4)
the grid snapshotted and the state's written; (5) [`Chain::stateIsValid`](../../src/bartcore/chain.hpp): shape,
each kept tree's form, each live tree built and, unforced, staged by
[`Tree::stageFromFlat`](../../src/bartcore/tree.hpp), the install's own routine, with its cone checked: validity
and the verdict, monotone feasibility included, are one loop, and after the first tree that would change the
rest are only built; (6) refused or declined: the grid put back, nothing else written; (7) a re-creation's store
resize; (8) `Chain::setState`: latents, trees staged again and repaired, kept draws copied at equal capacity or
the newest that fit repacked to slots 0 on, the generator's bytes last; (9) the ring's position. The bridge
passes null for the engine's `altered` report and returns FALSE for a decline before
[`reapplyWeights`](../../src/bartcore/facade.hpp) and its censoring twin. R's `setState` asks a logical
`forceUpdate` and, on an install, assigns the field `state` the object given, after a forced repair too
(dec-A199). Of the critique's rechecks: the copy's readers are steps 5 and 8; `altered` has no R reader, and
tests/cpp's [`verdictMatchesInstall`](../../tests/cpp/common.hpp) holds, its two calls completing with the same
coins and its comparison across a decline covering the generator put back. The ring: step 8's repack reads the
state's own cursor and capacity, so the draw over the state's R recorded draws precedes it. `$state` stays the
object given: section 3's first bullet stands. On FALSE the generators are back before step 4.

## 7. Callers

- Package (read): `setState`, `getPointer` and `copy` share `installStateOnto`; a re-created engine takes its
  flags from the data object, so a reload or copy of an older state draws. benchmarks (ran: grep; read): sbc.R,
  negbin-mixing.R and surfaces/C1-frozen-ess.R install a sampler's own state or a same-data sampler's.
- stan4bart (bartcore 3d2a295, read; f0084aa..3d2a295 is the one line forcing the install): restoreBartSampler
  installs, into a sampler built per chain from the fit's control, model and data, the state that chain's
  sampler stored, kept whole, with no predictor call or warm start (ran: git grep): S equals F, and saved fits
  carry no record. bartCause (e833be7) and treatSens (babfaa6) make no state, predictor or warm-start call,
  bairrtt (3f57f61) no state or warm-start call (ran: git grep).

## 8. Recorded baselines, snapshots and existing tests

Expected: nothing recorded moves (greps rerun by me on 82678655). The equivalence trio's scripts call no
setState, copy, installTrees or warm.start and read no `$state` (no factor is given a missing label after
creation: the critique's reading, not redone); the four reproducibility files hold none of those names, nor
setPredictor, setData or readRDS; of the 30 exact-gate scripts only negbin-mixing.R installs a state, its own,
and none warm-starts from a donor. So against the MANIFEST's current equivalence-e4faed5c,
bcf-equivalence-1b7d730c and multinomial-equivalence-80b1c8d4 every scenario reads "identical draws (same RNG
stream)", no |z| line, and this slice re-records nothing. sparse-starting-sigma re-records the equivalence
baseline for its own shift. Landing after it: rebase, compare against the file it made current and expect the
same, first rerunning a scenario that differs on the rebased base without this slice's commits, to tell whose it
is. Landing before it: compare against e4faed5c; the other slice re-records on top and inherits no shift.

Existing tinytest, run by me on 82678655 (the suite probe redone for the landed `installStateOnto` and for warm
starts: 72 files, all passing, 689 installs, 72 warm starts): 2 installs put a state with no right-going
direction onto a flagged column (test-data-missing.R, test-monotone-unforced.R), each the sampler's own or a
same-data twin's, so S equals F (read); 1 warm start has a donor that never flagged a column its recipient has:
test-monotone-unforced.R's cross-grid one, rewritten in step C. tests/cpp: 6 files call `setState` and mention
missing values (grep); not probed. A test pinning a defect this slice fixes is rewritten and named in the
landing note; any other move is a stop.

## 9. Steps, tests and mutants

Build order, a commit or more each: G, A1 (engine), A2 (bridge), A3 (tests, mutants, docs), B, C. Step G (TODO
setstate-test-gaps, first-missing-test-seeds) passes on the base and after every later commit; each of its
mutants is run against the base's code:
- G1, a pooled factor's kept store: [`testStateStoreSizes`](../../tests/cpp/test_state.cpp) gains a 70-level
  column, 10 kept into 4 and 4 into 10, the kept draws' mask channels and predictions the source's newest; the
  same from R in test-setstate-force-update.R. Mutant: the repack's copy of `savedTreeMasks` dropped or cleared,
  mean or variance (it crashed R in getTrees for the review: a crash is a kill).
- G2, a decline reconciles nothing (test-setstate-force-update.R): a logistic sampler given other weights and a
  predictor change that empties a leaf, then its older state unforced: FALSE, its generators, stored state and
  next 20 draws a twin's; an aft twin. Mutant: the decline returned after a reconciliation.
- G3, test-state-missing-direction.R: the three expectations after
  ["status <- tryCatch(other$setState(stale)"](../../inst/tinytest/test-state-missing-direction.R) and the
  forced one in `checkRestores` for the sparse-backed column never run (ran by me: the first call is refused as
  malformed, the second is TRUE). Each arm gets a searched seed and an assertion that it ran. Mutants: a forced
  install keeping a direction its column cannot route; a forced value visible.
- G4, one capacity (`testStateStoreSizes`): a two-chain state whose second chain's blocks, and a two-forest
  state whose second forest's block, hold another whole number of draws: refused in both forms, `getState`
  unchanged. Mutants: the comparison across chains dropped in `Sampler::setState`; across forests in
  `savedStateCapacity`. G5, a decline across cut grids (`testRestoreStatus`, and a tinytest block): a not-clean
  state on another grid, unforced: FALSE with cut points, predictions and the next 20 draws a twin's. Mutant:
  the decline not putting the grid back.
- G6: negbin-mixing.R, sbc.R and surfaces/C1-frozen-ess.R (the third found by grep) stop on a setState that is
  not TRUE; each drops the value today. G7: the help sentence (Docs).
- G8, seeds. test-missingness-first-seen.R: the row declined for emptying a leaf searches for a seed whose coin
  goes right (declined, asserted) and one whose coin goes left (taken); the two-forest block searches for the
  first seed whose forests each hold a rule on x1, live and kept, before the update (asserted) and both sides
  after. [`testMonotoneMissingArrives`](../../tests/cpp/test_monotone.cpp): the every-row form asserts that the
  plain sampler alone finds the row judged last valid. Mutants, Part A's: no draw on a second forest; a declined
  row's draw not taken back; the draw left on at a joint sweep's end.

Step A, tinytest: new test-state-missing-record.R on test-missingness-first-seen.R's fixture, seeds as in G8:
- the record: all FALSE on a complete sampler; TRUE for x1 after a first missing value by each path and after
  the column is filled; kept by saveRDS; of another length or type, an error. Copy and reload (asserted first:
  the cached record FALSE for x1, the data object's TRUE) hold the source's directions and next 5 draws (1e-12).
- the review's case, the update UNFORCED and asserted TRUE, on one forest and on two: `setState(old)` is TRUE
  and visible, and directions on x1 (live, kept), generators and predict on a missing row are those between the
  update and the install; a second update draws nothing. The same on a wrapped ring (asserted: cursor not 0,
  rules on x1 on both sides of the seam); 20 kept into a store of 5 keeps the newest 5 as drawn.
- nothing to draw: a state stored after the value and 10 sweeps keeps its directions and generator; one with the
  attribute removed installs as landed; directions on x1 into a complete twin give FALSE.
- a factor: a complete twin's state into a sampler created with missing values holds both directions where one
  can reach, none on a hand-built rule out of reach; kept rules with no category on the right, forced in
  (asserted present), keep their direction through a first missing value and an install that draws, and the
  sampler reinstalls its own state.
- not clean: a stored state whose drawn directions break a monotone tree's order, and a hand state whose drawn
  direction empties a leaf (Part A's lone-left tree), each coin asserted: FALSE unforced with getTrees,
  generators, `$state` and the next 20 draws a twin's; NULL forced; else TRUE. Completes to invalid: a hand
  state, record FALSE, a live exempt rule beneath a rule on its column: on the coin that takes the missing value
  away (asserted) an error, the sampler a twin's.

Step A, tests/cpp, `testStateMissingRecord` in test_state.cpp, after
[`testMissingFirstSeen`](../../tests/cpp/test_sampler.cpp): `getState` writes the flags; state, one whole-column
first missing value, state again, install the first: trees, kept draws and generator bytes equal the second's
(two chains, two forests, a variance forest, a pooled column, a wrapped ring, a store of another size); exempt
rules, nested and out of reach, take no coin in either caller; cleared `rngState` draws from the chain's
position and leaves it on a decline; a decline and a refusal after the draw leave generators and grid; the fuzz
states with a column cleared from the record complete to valid ones. Mutants, each failing a test: the record
not written, written from content, or through getMissingSeen; the raised set reversed; no record read as none
missable; no draw on live trees, kept draws, the variance forest, a second forest or later chains; coins from
the chain's generator; the state's generator installed without the advance; kept draws in slot order, by the
recipient's ring, or after the repack; a draw at every factor rule; a coin for an exempt rule; step 5 or 8
reading the state as stored; a generator not put back; a flag raised.

Step B tests, 16 pins. test-monotone-unforced.R, 12: the section returns to what the engine does on the tree a
missing f breaks, each outcome checked against the side an unconstrained twin on the same seed draws (FALSE and
untouched, the row declined, NULL and reseeded when forced, TRUE where the order is kept).
test-joint-update-factor.R, 3: every sampler takes a missing label; only the third, complete in f, draws.
test-data-categorical-declared.R, 1, setTestPredictor's: still refused, in the test path's words, with a numeric
twin. test-missingness-first-seen.R runs its paths for f too. Mutants: the refusal left on one path; a missing
label coded as a level; the test setter taking one.

Step C tests. test-state-missing-record.R: the measured case holds both sides on x1 before any sweep, and over
40 seeds the share sent right is within the binomial 1e-3 band of one half; two chains given one donor draw
(`samples`) differ; a donor state with the attribute removed, and a donor created with the missing values, leave
the recipient's generators unmoved, all left and as learned. The draw before the collapse: a hand donor whose
one split on x1 lies above every observed value of the recipient, on its grid and on another: where the coin
goes right (seed searched, asserted) the split stands with the missing rows beneath it, where left it is merged.
test-monotone-unforced.R's block "a warm start from a donor on another cut grid": the recipient's seed is
searched, as the setData arrival beside it is, for the sides that leave `breaks` out of order (reseeded, as
pinned) and for a pair that keeps it (read, not run). tests/cpp, `testWarmStartMissingDirections` in
test_state.cpp: each chain's directions equal a replay of its own generator in tree order, variance trees after;
a slot-sourced exempt rule takes no coin; no record, and a refusal after the draw, leave trees and generator
bytes. Mutants: no draw; the draw after the collapse; coins from the donor state's generator, or chain 0's for
every chain; no record read as none missable; a column drawn where D and F agree; variance trees skipped; the
generator left advanced on a refusal.

Docs. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd). Step G, Saving: after a forced setState that
repaired a tree the `state` field holds the state as it was given, so installing the field again unforced
returns FALSE and a copy or a reload repairs it the same way; storeState() stores what the sampler now holds.
Step A: the same of an install that drew a side; the `state` field gains "and which predictor columns could hold
a missing value when it was stored"; Missing values in predictors says what the Goal's first two sentences say,
in the help's words (side, rule, saved draws), and that setState returns TRUE for such an install unless a drawn
side leaves a leaf with no row or a monotone tree out of order. Step B: the refusal of "a missing label in a
column that has never held a missing value" becomes setTestPredictor's alone; the paragraph beginning "A column
coded from a factor" goes. Step C, there and in installTrees's docstring: a warm start from a fit whose data
could not hold a missing value in a column draws the side of every rule on that column in the trees it installs.
[NEWS.Rd](../../inst/NEWS.Rd): its first-missing-value item gains one clause (0.9-34 had no warm start). The
design note is amended. TODO items closed at landing: state-before-first-missing,
factor-first-missing-setpredictor, first-sight-kept-one-way-rule, test-column-missing-not-refused,
setstate-test-gaps, first-missing-test-seeds.

## Verification

Independently of the implementer, `R_LIBS=$LIB` on every R call, one R process at a time
([RNG classes and their gates](README.md#rng-classes-and-their-gates): shifting, nothing to re-record expected):
- `R CMD INSTALL --preclean --library=$LIB .` per engine commit; `cd tests/cpp && make && ./test_bartcore`,
  plain and with `OPT="-O2 -g -fsanitize=address,undefined"`; `tinytest::run_test_file` on each touched file,
  also under R-loaded ASAN; `tinytest::test_package("dbarts", at_home = TRUE)`: no failure after any step.
- On a reference build: the trio's `compare <file> --bitwise` (equivalence.R also `--strict-coverage`) on
  section 8's files, in the landing order it names: 55, 15 and 11 scenarios, each "identical draws (same RNG
  stream)"; the four `test-reproducibility-*.R` files pass unchanged.
- Every script of `.github/workflows/exact-gates.yaml` with `quick`; `R CMD check --as-cran` from a clean
  tarball; lintr, `air format --check .`, and tools/check-rc-codoc.R, check-win-drift.R and
  check-doc-freshness.R; stan4bart's suite at 3d2a295 on the build, and a fit saved on the base reloaded.

## Calls made

- An install that drew is clean: unforced TRUE. dec-B435 names the draw, not the value; the grounds are
  dec-B305's clean ("a coherent state with regards to its model"), dec-B310's "statistically valid update",
  dec-B321 (a direction no row reached has its prior as its conditional) and the unforced setPredictor that
  draws and returns TRUE (ran). After such a TRUE the directions and generator are not the given state's, and
  `$state` is still the object given, as after a forced repair. Alternative: FALSE, drawing only when forced.
- The draw completes a copy before anything is judged, with coins from the state's generator (alternatives: draw
  after the install, judging under the left default; the sampler's generator, which the install overwrites).
- New, against the landed code. The draw sits before the grid is written and runs its own walkability pass,
  where the plan had it after a first validation: the landed validity and verdict are one loop behind the grid
  write (alternative: split `stateIsValid`). Only `state.chains` is copied. A state without generator bytes
  hands its advance to the copy, so a decline leaves the chain's generator (the plan advanced it in place). The
  warm start draws in the scratch `installForests` builds and undoes its coins on a refusal (alternative: draw
  after the checks, which would then judge other trees). Step G is built first, so its tests guard steps A to C.
- The exempt rule covers "no category on the right" only: its mirror, every category right and the missing value
  left, stays well formed either way (rerun). No record draws nothing; the attribute is `missing.columns`,
  always written; every recorded kept draw of the state draws, before the repack. Step B: the test column form
  takes the test path's refusal (dec-B322; alternative: keep the factor refusal there).

## Budget, stops and evidence

Planned ~1035. Step G ~105: tests 95 (G8 about 40), benchmarks 6, help 4. Step A ~620: engine 135, bridge 30,
tinytest 230, tests/cpp 180, help, NEWS and docs 45. Step B ~185: R 25, tests 150, help 10. Step C ~125: engine
35, bridge 15, tests 70, help 5. Forecast at 1.5 to 2 times, as this surface has run (Part B: 1.8): 1550 to
2070; stop at 1.5 times the midpoint, ~2720. Stop and report, without working around, when: the slice passes its
stop; anything recorded, or an existing test other than section 8's, moves; section 3's scoped loop does not
reproduce the update; a `getState` state completes to an invalid one; a step G test fails on the base (a defect
in landed code: report it); the draw cannot sit before the grid write without changing `stateIsValid`'s
signature; step B needs an engine change; step C's refusal needs more than the generators put back.

Evidence. Run by me on 82678655 (library scratch/libs/smr2; scripts in scratch/smr/ and its crit/): what the
sections mark rerun, the suite probe, the arms of test-state-missing-direction.R, the greps. Read there: every
claim cited by symbol; Part B's landing note; dec-A199, dec-B440; step G's two TODO items. Not run, no build of
the slice existing: the completion, the exempt rule, the warm start's draw, every mutant and new test; tests/cpp
was not built; the review's crash under the masks mutant and 06 (0.9-34) were not rerun.

## Landing note

Landed 2026-10-11 as a66b8f65..d180dcc2: one opus implementer, one opus review that ran 55 mutants of its own and found
nothing blocking, a fix round of tests, and an x86 run of tests/cpp plain and under the address, undefined-behavior
and leak sanitizers and of tinytest. Nothing recorded moved: the equivalence baselines (against equivalence-2b48939b,
which replaced the baseline section 8 names while the slice was in review), the snapshot files and the exact gates
are bit for bit; every install that draws nothing is what the base build did (253 of 253 in the review's comparison).
What was built is sections 1 to 5 and step G, with these departures and additions (dec-A201):

- An install that drew is a clean install: TRUE unforced, NULL forced. The sampler's cached state still shows the
  state as given.
- A warm start whose donor tree cannot be walked makes no draw and keeps the refusal it had; the plan returned a new
  one.
- A record of another length is refused by the engine; the bridge refuses a record that is not logical, holds NA or
  is empty.
- The plan's randomized test was unsound: a sampler's own state with a flagged column cleared from its record can
  complete to a state it would never produce (refused in 26 of 40 seeds). In its place an earlier state completed at
  an install must reinstall unforced and round-trip (826 of 826 at 600 seeds).
- The engine tests of the record and of the warm-start draw sit beside the first-seen fixture they reuse.
- The help sentence that a factor column does not take missing values back through a column update is removed.
- One test's lambda took a local as its default argument, which g++ refuses and clang accepts; found on the x86 run.
- Size: about 2,400 lines added beside this plan against about 1,035 planned; about 2,000 of them tests.

The live defect of section 6 is fixed: over 91 seeds, 14 held a kept one-way rule; the sampler reinstalled its own
state after a first missing value, unforced and forced, in 14 of 14, for 1 of 14 on the base build. Not verified:
Windows, more than two threads, bench-sampler's timing.
