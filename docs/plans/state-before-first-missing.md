# state-before-first-missing: a state records which columns could be missing, and an install draws what it lacks

Status: PLANNED (dec-B435; read with dec-B305, dec-B310, dec-B318, dec-B320 to dec-B322, dec-B378, dec-B398 to
dec-B402, dec-A198). One open call (O1). Written against Part B of [setstate-force-update.md](setstate-force-update.md)
AS PLANNED; What touches Part B lists the rechecks.

agent: opus implementer, one; one opus reviewer who runs the mutants below; a blind critique first.
rng: SHIFTING only for an install whose state's record lacks a column the receiving sampler has flagged: each
chain then draws coins no install draws today. NEUTRAL, bit for bit, for every other install, copy and reload, a
state with no record included: no copy is taken and no generator is read.
window: after Part B lands; before the merge to main; engine slices stay serial.
budget: planned ~510 lines (engine ~115, bridge ~30, tests ~320, help and docs ~45); forecast and stops below.

## Goal

A stored state says which predictor columns could hold a missing value when it was stored. When `setState`, a
copy or a reload installs it into a sampler where a further column can hold one, the direction of every split on
that column in the state's live trees and kept draws is drawn at one half, before anything is judged or
installed. Any other install is what it is today. Tier: "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

dec-B435 (extends dec-B400), the maintainer on 2026-10-10: "Sure, draw it." The register's rule: "a state
records the columns that could hold a missing value when it was stored; installing it into a sampler where a
further column can draws, with probability one half, the direction of every split on that column in the
installed trees and kept draws, as a first missing value does." Ran on c4ea26ab (scratch/smr/01-review-case.R:
150 rows, 2 chains, 15 trees, 20 kept draws): a state stored, then a forced setPredictor putting 3 missing
values in x1: live trees 3 left and 5 right, kept draws 95 and 95. `setState(old)`: live 8 and 0, kept 190 and
0, TRUE, the generators the old state's; predict answers a missing row; a second update with missing values
draws nothing. Two more routes give 8 and 0, 190 and 0: a reload and a `copy()` of a sampler whose `$state`
predates the value (its source holding 3 and 5, 95 and 95), and that state into a twin created with the values.

Read on c4ea26ab: the first-sight draw, [`Chain::drawMissingDirections`](../../src/bartcore/chain.hpp), takes
the live trees (forests in order, then the variance forest), then the kept draws oldest first. Kept draws are
flat trees, drawn by [`drawFlatMissingDirections`](../../src/bartcore/tree.hpp), which writes the direction to
the record's flag whatever the column's kind and applies dec-A198's reach rule. A state's live trees are flat
too, and [`Tree::buildFromFlat`](../../src/bartcore/tree.hpp) reads that flag for every kind of rule.

## 1. The record

- Engine: [`SamplerStateData`](../../src/bartcore/sampler.hpp) gains `missingSeen`, a `std::vector<std::uint8_t>`,
  one byte per predictor, empty meaning no record. [`Sampler::getState`](../../src/bartcore/sampler.hpp) copies
  the store's flags ([`hasMissing`](../../src/bartcore/data.hpp)) into it, not the content, so a column since
  filled stays recorded. Sampler-level: the chains and forests share the store.
- Stored object: a top-level attribute `missing.seen` on the bartcoreState, a logical with one element per
  predictor column and no NA, as the data object's record is ([`dataMissingSeen`](../../R/data.R)). The bridge's
  [`storeState`](../../src/R_interface_bartcore.cpp) writes it, and every R writer goes through that (the
  method, the `state` field's promise, `copy`; read). No R code builds or edits a state (ran: grep), so R does
  not change. It rides saveRDS, and a copy installs its source's object as is. One that is not such a logical,
  or of another length, is an error in both forms ("malformed missing-value record in bartcore state").
- A state WITHOUT the attribute (only on this unreleased branch and in stan4bart fits saved on it) has no
  record and draws nothing: it installs exactly as on the landed build. Safe: all that is lost is this slice's
  draw, which stan4bart cannot reach (section 4). Read as "no column could be missing" it would redraw directions
  a chain learned; refused, those states would be orphaned for nothing.
- Format: additive. [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)'s registry rule reads top-level
  attributes by name, defaults an absent one and holds the version and the floor at 1 until a release. The flat
  C API has no state entry (read; ran: grep of its entry table): no signature changes, no ABI event, and
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move (dec-B424 does not arise). A fit
  holding a state carries one attribute more, so the exact gates run in quick mode.

## 2. The install rule

With S the state's record and F the sampler's flags at the call:
- In F, not in S: drawn (dec-B435). `raised[j] = F[j] && !S[j]`.
- In S, not in F: unchanged; the record plays no part. As built and as Part B plans, a direction on a threshold
  or inline subset rule is dropped and the install not clean (Part B's case c), a pooled rule keeps its bit,
  which routes nothing, and a subset rule sending every level one way is refused in both forms (dec-B378). In
  both, a direction is kept, clean (dec-B320). A state never raises a flag: that takes data (dec-B321).
- Equal, or no record: nothing, and the caller's state is installed as today.

The draw completes a copy of the state before anything is judged. In
[`Sampler::setState`](../../src/bartcore/sampler.hpp), after the last check that reads no direction
([`Chain::stateIsValid`](../../src/bartcore/chain.hpp), which has found every flat tree well formed) and before
the first that does (the monotone cone, then Part B's verdict), when a column is raised: the state is copied
and each chain draws on its part, from its own generator (section 3), with `drawFlatMissingDirections` in the
first-sight order: each forest's live trees, the variance trees, then the state's recorded kept draws oldest
first by the state's own ring ([`savedSlotForDraw`](../../src/bartcore/sampler.hpp)'s rule on the state's
cursor, count and capacity), each draw's forests then its variance trees. The copy is validated again, which
refuses only a hand-edited state whose rules contradict its record, and stands in for the state from there on:
leaves are counted and merged, and the cone judged, under the directions the rows will be routed by.

- Unordered factors (dec-A198): the routine draws a subset rule only where a missing value reaches it and
  leaves it left elsewhere, the gauge `buildFromFlat` checks, so the sampler reinstalls its own state after.
  Ran (02-other-routes.R): a complete twin's state into a sampler created with 3 missing values in a factor of 4
  levels gives 10 and 0 live, 238 and 0 kept, TRUE today.
- Unforced: a draw does not make an install unclean. dec-B305, the maintainer: clean is "leaving a sampler in a
  coherent state with regards to its model", and its rule names as not clean "a tree that would have to be
  changed to fit". A drawn direction changes no tree to fit: no row had reached it, so its conditional was its
  prior (dec-B321's reasoning). dec-B435 put "reporting the install as not clean" to the maintainer as the
  alternative, and he took the draw; dec-B320 and dec-B318 ("Option A, clean.") make neither a direction nor
  the kept store a ground for FALSE. So TRUE, unless the completed trees fail the verdict on another ground (a
  leaf no row reaches, a monotone tree out of its cone): then FALSE on every repeat, the sampler untouched.
- Forced (so every copy and reload, dec-B305, and stan4bart's restore, dec-B402): completed, then installed
  with Part B's repairs under the drawn directions; NULL invisibly (dec-B310).
- Refused, not clean or thrown: nothing of the sampler is written but the grid, put back as today. The copy is
  taken first; each chain's generator is borrowed for its coins and put back, byte for byte, before the next.

## 3. The generator

A state carries each chain's generator, which [`Chain::setState`](../../src/bartcore/chain.hpp) installs last.
The coins come from the state's generator: per chain, the chain's generator is serialized, the state's bytes are
read into it, the coins are drawn, the advanced bytes go into the copy, and the chain's own bytes are read back
(a dedicated generator per chain and an exact round trip, as Part A relies on; read). The install then reads the
advanced bytes, and the bridge's redraw of latents under other weights follows as today. A state whose generator
is of another kind leaves the chain its own (as built): that one draws, and stays advanced only on an install.
- One state installs the same way every time. So the cached `$state`, left the object `setState` was given,
  still describes the sampler: a copy or reload of it holds the directions the live sampler holds.
- The loop store; a first missing value, no sweep between; `setState(stored)` is reproducible: the update drew
  from the generator the state recorded, over the same trees in the same order, so the install draws the same
  coins and leaves the generator where the update left it. Derived from the two orders (read); a test pins it.
- Nothing to draw: the routine returns before the copy and before any generator is serialized, so the install
  runs today's instructions on today's bytes.

## 4. Callers

- Package (read): `setState`, `getPointer` and `copy` share [`installStateOnto`](../../R/dbarts.R); a re-created
  engine takes its flags from the data object, so a reload or copy of a state stored before the value draws.
  benchmarks (ran: grep; read): sbc.R and negbin-mixing.R install a sampler's own state and C1-frozen-ess.R
  transplants into a fresh sampler on the same data, so S equals F in each.
- The warm start ([`Sampler::installForests`](../../src/bartcore/sampler.hpp)) installs trees, not a state: O1.
- stan4bart (bartcore f0084aa, read): restoreBartSampler builds a sampler per chain from the fit's control,
  model and data and installs the state that chain's sampler stored at the fit's end. Its R/ calls neither
  setPredictor nor setData (ran: grep) and dbarts.h has no predictor entry, so S equals F; its saved fits
  carry no record and install as today. No change there.
- bartCause (dbarts-1.0 e833be7), treatSens (dbarts-1.0 babfaa6): none of setState, storeState, copy,
  installTrees, setPredictor, setData. bairrtt (main 3f57f61): predictor updates, no state call (ran: git grep).

## 5. Equivalence, snapshots and gates

Expected: nothing recorded moves (ran: grep). The equivalence trio installs no state and reads no `$state`,
its missing values there at creation; under benchmarks/R only negbin-mixing.R, C1-frozen-ess.R and sbc.R call
`setState`. The four reproducibility files hold none of setState, setPredictor, setData, copy, readRDS,
installTrees, warm.start, `$state` or NA. Of the exact gates' scripts only negbin-mixing.R installs a state, its
own. Every scenario reports "identical draws (same RNG stream)"; a move is a stop. Existing tinytest (ran:
03-suite-probe.R, an R-level wrap of `installStateOnto` over the 61 files that install a state and either change
a predictor or write NA; 522 installs, 9521 expectations passing): 2 installs put a state with no right-going
direction onto a flagged column, in test-data-missing.R and test-monotone-unforced.R, each a sampler's own or a
same-data twin's state (read). None is expected to move. tests/cpp was not probed.

Gates ([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting, nothing to re-record
expected): tests/cpp plain and under the sanitizers; R-loaded ASAN over the touched test files; full tinytest; a
`--preclean` reference build for the trio's `compare --bitwise` and the four snapshot files; the exact gates,
quick; `R CMD check --as-cran`; lint, air, the doc checks; stan4bart's suite, and a stan4bart fit saved on the
Part B tip reloaded here with equal predictions.

## 6. Steps, tests, mutants and docs

Steps: (1) engine: the field, `getState`, the completion, the second validation (`--preclean` installs);
(2) bridge: the attribute written and read, the registry comment; (3) tests and mutants; (4) docs and records.

tinytest, new test-state-missing-record.R on Part A's fixture
([test-missingness-first-seen.R](../../inst/tinytest/test-missingness-first-seen.R)):
- the record: all FALSE on a complete sampler; TRUE for x1 after a first missing value by each path, and after
  the column is filled; kept by saveRDS.
- the review's case: store; forced setPredictor with 3 missing values; `setState(old)` is TRUE and the
  directions on x1, live and kept, the generators and predict on a missing row are identical to those taken
  between the update and the install; a second update with missing values draws nothing; a copy and a reload
  of a sampler whose `$state` predates the value hold its directions, their next 5 draws equal within 1e-12.
- nothing to draw: a state stored after the value and 10 sweeps keeps its learned directions and generator; one
  with the attribute removed installs as landed (all left, TRUE); a record of another length or type is an error.
- two forests and a variance forest: each forest's rules, live and kept, equal what the update drew; kept
  draws: a wrapped ring, and 20 kept into a store of 5 (Part B's B4) keeps the newest 5 as the update drew them.
- a factor: the complete twin's state into a sampler created with missing values holds both directions where
  one can reach and none on a hand-built rule out of reach; the sampler then reinstalls its own state, TRUE.
- not clean: a hand state whose drawn direction empties a leaf (Part A's lone-left tree) gives FALSE with
  getTrees, the generators, `$state` and the next 20 draws identical to a twin's; forced, NULL, the leaf merged.
- S beyond F: a state with directions on x1 into a complete twin is FALSE; the twin's record stays NULL.

tests/cpp, `testStateMissingRecord` in [test_state.cpp](../../tests/cpp/test_state.cpp), by the pattern of
[`testMissingFirstSeen`](../../tests/cpp/test_sampler.cpp): `getState` writes the flags; state, first missing
value, state again, install the first: trees, kept draws and generator bytes equal the second's, over two
chains, two forests, a variance forest, a pooled column and a wrapped ring; another generator kind draws from the
chain's own, put back on FALSE; the fuzz states with a column cleared from the record complete to valid states.

Mutants, each failing a test: the record not written, or written from content; the raised set reversed; no
record read as none missable; no draw on live trees, kept draws, the variance forest, a second forest or chains
past the first; coins from the chain's generator; the state's generator installed without the advance; kept
draws newest first, or drawn after the repack; a draw at every factor rule; the verdict or the cone judged on
the state as stored; a generator not put back on FALSE; the copy not validated; a flag raised from the record.

Docs:
- [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd), the `state` field: "and which predictor columns
  could hold a missing value when it was stored". Missing values in predictors gains: "A stored state records
  which columns could hold missing values when it was stored. Installing it, by \code{setState}, \code{copy} or
  a reload, into a sampler where a further column can hold them draws the side of every rule on that column,
  in the state's trees and saved draws, as a column's first missing value does, from the generator state the
  state carries, so the same state installs the same way each time. \code{setState} counts such an install as
  clean and returns \code{TRUE}, unless a drawn side leaves a leaf with no row." Saving, on what TRUE covers:
  "and sides drawn for a column that could not hold a missing value when the state was stored".
- [NEWS.Rd](../../inst/NEWS.Rd): the first-missing-value item gains "; a state stored before that value and
  installed after it has them drawn at the install". [mia-missingness.md](../design/mia-missingness.md), Flat
  format and serialization, is amended; Part B's design note gains a row; the TODO item closes.

## What touches Part B

Against [Part B: setState's two forms (dec-B305, dec-B310, dec-B318, dec-B398, dec-B401)](setstate-force-update.md#part-b-setstates-two-forms-dec-b305-dec-b310-dec-b318-dec-b398-dec-b401);
recheck each when it lands:
- B1, `Sampler::setState`'s body and signature: the completion sits inside it and reads neither `force` nor
  `notClean`. Recheck where the checks stand and that every reader after it takes the copy.
- B2, the verdict and its staging routine read the completed copy, as do case f and B8's repairs; case c stays
  keyed on the tree's bits and the sampler's flags. Recheck that the verdict is inside `Sampler::setState`.
- B4, the kept-store repack: the draw runs before it, over the state's ring. Recheck its capacity and slot rule.
- B6: the bridge reader gains the attribute beside `force`, with no arity change; R is not edited. B9 and the
  design note: this slice's sentences go into text Part B rewrites. Tests use Part B's forms and its fuzz check.

## Open call

O1. Should a warm start draw these directions too? A warm start (bart's `warm.start`, a sampler's
`installTrees`) copies a donor fit's trees into a new sampler as its starting point. If the donor's data were
complete in a column and the new sampler's are not, the copied trees hold no direction for missing values
there, and every split on the column starts out sending them left. Run on 2026-10-10 (150 rows, 3 of them
missing in x1, 15 trees): 7 of 7 splits on x1 send left at the start; 8 left and 1 right after 1 sweep, 8 and 4
after 31, so the start wears off slowly. dec-B435 names a state's install; a warm start takes trees and no
chain. With this slice the donor's state says which columns could be missing, so the same draw is available.
- (a) Draw: each receiving chain draws at one half, from its own generator, the direction of every split on a
  column the donor could not have missing and the recipient can. About 110 lines (engine 30, bridge 15, tests
  60, help 5). Such a warm start's draws shift; no equivalence scenario, snapshot or exact gate has one (grep).
- (b) Leave it: the splits start left and the sampler's moves redraw them over the run; about 3 lines of help.
Recommended: (a), so trees grown where a column could not be missing get their directions the same way by every
route. How realistic: a refit seeded from an earlier fit after rows with gaps arrive in a column that had none;
uncommon, and it moves a chain's starting point, not a stored posterior draw.

## Calls made

- The draw completes a copy of the state before the verdict and the install. Alternative: install, then draw
  and refresh as a forced update does, merging leaves and judging under the left default the draw replaces.
- Coins from the state's generator. Alternative: the sampler's current one, which the install overwrites; one
  state would install differently each time, and a reload of `$state` would not reproduce its source.
- No record draws nothing (section 1); a state never raises a flag (section 2); the alternatives are there. The
  record is a top-level logical attribute, `missing.seen` (alternatives: a block per chain; raw bytes).
- Every recorded kept draw of the state draws, before the repack, so the coins do not depend on the recipient's
  store (alternative: only the draws it keeps). `$state` after `setState` stays the object given, as today
  (alternative: store again after a draw). The copy is validated again, on the drawing path only (alternative:
  one check with the raised flags lowered).
- The engine entry keeps its const state and copies only when a column is raised. Alternative: a mutable
  state, changing the 87 tests/cpp call lines Part B counts.

## Budget, stops and evidence

Planned ~510: engine 115, bridge 30, tinytest 170, tests/cpp 150, help, NEWS and docs 45; O1 as (a) adds ~110.
Forecast at 1.5 to 2 times, as this surface has run: 770 to 1020; stop at ~1340 (~1630 with O1). Stop and
report, without working around, when: the slice passes its stop; a scenario, a snapshot or an existing test
moves; section 3's loop does not reproduce the update's directions and generator; a `getState` state completes
to an invalid one; Part B landed with the verdict outside `Sampler::setState` or another ring rule.

Evidence. Ran on c4ea26ab, library scratch/libs/smr, probes and outputs in scratch/smr/: the review's case and
its copy, reload and twin forms (01); the factor and warm-start forms (02); the install probe (03); the greps of
sections 1, 4 and 5. Read on c4ea26ab: every claim cited by symbol; the rulings in Status; Part B as planned;
stan4bart at f0084aa. Not run, there being no build of the slice: the loop of section 3, the factor gauge of a
completed state, and any tests/cpp pin installing one sampler's state into another that flags more columns.
