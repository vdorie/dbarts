# state-before-first-missing: a state records which columns could be missing, and an install draws what it lacks

Status: PLANNED (dec-B435 and dec-B436; read with dec-B305, dec-B310, dec-B318, dec-B320 to dec-B322, dec-B378,
dec-B398 to dec-B402, dec-A198). Revised from its blind critique (last paragraph). One open call (O1). Written
against Part B of [setstate-force-update.md](setstate-force-update.md) AS PLANNED.

agent: opus implementer, one; one opus reviewer who runs the mutants below.
rng: SHIFTING for an install whose state's record lacks a column the sampler has flagged, for a factor column's
first missing value by an update refused today, and for a first-sight draw over a kept rule of section 2's
exempt form. NEUTRAL, bit for bit, otherwise: an install with nothing to draw takes no copy, reads no generator.
window: after Part B lands; before the merge to main; engine slices stay serial.
budget: planned ~805 lines (step A ~620, step B ~185); forecast and stops below.

## Goal

A stored state says which predictor columns could hold a missing value when it was stored. When `setState`, a
copy or a reload installs it into a sampler where a further column can hold one, the direction of every split on
that column in the state's live trees and kept draws is drawn at one half, before anything is judged (dec-B435).
A factor column takes its first missing value as a numeric column does (dec-B436). Tier: "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

dec-B435 (extends dec-B400), the maintainer on 2026-10-10: "Sure, draw it." The register's rule: "a state
records the columns that could hold a missing value when it was stored; installing it into a sampler where a
further column can draws, with probability one half, the direction of every split on that column in the
installed trees and kept draws, as a first missing value does." dec-B436, the same day: "Lift the refusal."
Ran on c4ea26ab (01-review-case.R: 150 rows, 2 chains, 15 trees, 20 kept draws): a state stored, then a forced
setPredictor putting 3 missing values in x1: live trees 3 left and 5 right, kept draws 95 and 95; `setState(old)`:
live 8 and 0, kept 190 and 0, TRUE. A reload and a `copy()` of a sampler whose `$state` predates the value, and
that state into a twin created with the values, give the same.

Constraints: no change to the state-format version or floor, to [dbarts.h](../../inst/include/dbarts/dbarts.h)
or to stan4bart; nothing recorded may move (section 6); Part B's forms are used, not changed. Out of scope: the
warm start unless O1 is ruled (a); keeping, not redrawing, directions a kept draw already holds on a
never-flagged column (dec-B435's "every split on that column" stands); the mirror of section 2's exempt rule.

## 1. The record (step A)

- Engine: [`SamplerStateData`](../../src/bartcore/sampler.hpp) gains `missingColumns`, one byte per predictor,
  empty meaning no record; [`Sampler::getState`](../../src/bartcore/sampler.hpp) always fills it from the store's
  flags ([`hasMissing`](../../src/bartcore/data.hpp)), not the content. Of another length, `setState` returns false.
- Stored object: a top-level attribute `missing.columns` on the bartcoreState, a logical per predictor column
  with no NA; a malformed one is an error in both forms. Always written, all FALSE on a complete sampler, and
  never through [`bartcore_getMissingSeen`](../../src/R_interface_bartcore.cpp), which answers NULL when no flag
  is up; the name differs from the data object's `missing.seen`, where NULL means no column: absence here means
  not known. Every R writer goes through the bridge's [`storeState`](../../src/R_interface_bartcore.cpp) (read),
  and no R code builds or edits a state (ran: grep).
- A state WITHOUT the attribute (only on this unreleased branch and in stan4bart fits saved on it) has no
  record and draws nothing: it installs exactly as on the landed build, which stan4bart's use cannot fault
  (section 5). Read as "no column could be missing" it would redraw directions a chain learned.
- Format: additive under [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)'s registry rule (attributes
  read by name, an absent one defaulted): version and floor stay 1. The flat C API has no state entry (read;
  ran: grep): no signature change, no ABI event, no API-hash move. The exact gates run quick: a fit's state grew.

## 2. The install rule (step A)

With S the state's record and F the sampler's flags at the call:
- In F, not in S: drawn (dec-B435). `raised[j] = F[j] && !S[j]`.
- In S, not in F: unchanged; the record plays no part. As built and as Part B plans, a direction there is
  dropped and the install not clean (case c), a pooled rule keeps its bit, and a subset rule sending every
  level one way is refused (dec-B378). A state never raises a flag: that takes data (dec-B321).
- In both, a direction is kept, clean (dec-B320). Equal, or no record: nothing; the state installs as today.

The draw completes a copy of the state before anything is judged. Its place in
[`Sampler::setState`](../../src/bartcore/sampler.hpp) is set by what it needs: every flat tree, live and kept,
walkable (the routine recurses with no bounds), and nothing yet built, routed or judged. When a column is
raised: each flat tree is checked with [`flatTreeIsWellFormed`](../../src/bartcore/tree.hpp), the state is
copied, and each chain draws on its part from its own generator (section 3) with
[`drawFlatMissingDirections`](../../src/bartcore/tree.hpp) in the first-sight draw's order: each forest's live
trees, the variance trees, then the state's recorded kept draws oldest first by its own ring. The rest runs once,
on the copy: [`Chain::stateIsValid`](../../src/bartcore/chain.hpp), the cone, Part B's verdict, install, repairs.

- The exempt rule (critique M1). A subset rule with no category on the right sends only the missing value
  right: its direction is fixed by its form, not drawn. `drawFlatMissingDirections` leaves such a rule as it is,
  in reach or out of it, takes no coin for it and reads its direction for the rules beneath (it is handed the
  mask channel for a pooled rule), for both its callers: the first-sight draw's kept loop and this install.
  Otherwise a coin to the left leaves a rule sending nothing right, which the well-formedness check refuses.
  Ran (scratch/smr/crit/, the critique's p2 probes rerun, same output): a sampler complete in a 4-level factor
  that takes the state of one with 30 missing values in it holds 198 right-going directions on the factor in
  its kept draws, its record all FALSE, one on such a rule; on the landed build a later first missing value by
  setData flipped such rules in 11 of 11 seeds that held one, and the sampler could no longer install its own
  state. The mirror (every category right, the missing value left) keeps its coin and stays well formed (ran).
- Unordered factors: the routine's reach rule (dec-A198) is the gauge `buildFromFlat` checks on the copy.
- Unforced: TRUE for an install that drew (Calls made), unless the completed trees fail the verdict on another
  ground (a leaf no row reaches, a monotone tree out of its cone): FALSE, on every repeat. Forced (so a copy, a
  reload, dec-B305, and stan4bart's restore, dec-B402): completed, then repaired as Part B has it; NULL.
- Refused, not clean or thrown: nothing of the sampler is written but the grid, put back as today; each chain's
  generator is borrowed inside a scope guard and put back, byte for byte, on any exit. The validation after the
  draw refuses, as an error in both forms, only a hand-edited state that contradicts its record.

## 3. The generator (step A)

A state carries each chain's generator, which [`Chain::setState`](../../src/bartcore/chain.hpp) installs last.
The coins come from it: per chain, the chain's generator is serialized, the state's bytes are read into it, the
coins are drawn, the advanced bytes go into the copy and the chain's own bytes are read back (a dedicated
Mersenne Twister per chain, an exact round trip; read). A hand-built state with no generator bytes draws from
the chain's own, which stays advanced on an install.
- One state installed twice into one sampler draws the same directions, so the cached `$state`, left the object
  `setState` was given, still describes the sampler: a copy or reload of it holds the live sampler's directions.
- Reproducible, in this scope only: a store; ONE whole-column or whole-matrix setPredictor raising its columns
  in one call, with no generator use between; `setState(stored)`. The install then draws the update's coins and
  leaves the generator where the update left it. Ran by replay in R (the critique's p1 probes rerun; 2 chains,
  ring cursor 7 of 20, 134 and 55 coins): unforced, forced, and two columns in one call, all equal. NOT claimed
  and not a stop: the per-observation and joint paths (chain 1 first draws the scan order, ending 283 draws on,
  not 134); columns raised by separate calls; a sweep between; setData.
- Nothing to draw: the routine returns before the walkability check, the copy and any generator read.

## 4. A factor column's first missing value (step B, dec-B436)

Its own step and commit, after step A, whose exempt rule a factor's first-sight draw can meet. No engine change.
Ran (04-factor-first-missing.R, the refusal bypassed by hand, 3 missing labels in a 4-level factor): forced
column, unforced column and whole data frame, live 4 left and 6 right, kept 118 and 120; row by row and jointly,
150 of 150 rows, live 4 and 6, kept 115 and 123; each reinstalls its own state, predicts and copies.
- The refusal "has missing values, which its training values do not" and the `missing.seen` argument that fed
  it go from [`codeCategoricalColumnUpdate`](../../R/bartcore.R), [`codePredictorFrame`](../../R/bartcore.R) and
  their callers. A missing label is coded as a missing value; an unknown label and a number stay refused.
- The data object's `missing.seen` for a factor is written after every accepted change, as for a numeric column
  ([`recordMissingSeen`](../../R/data.R); ran: FALSE FALSE TRUE on each path), and no training update reads it.
- Row by row and jointly, the engine's session ([`UpdateSessionImpl`](../../src/bartcore/sampler.hpp)) draws at
  the first missing row and takes the draw back if the row is declined; only a sampler new to it draws.
- setTestPredictor is not named by dec-B436 and keeps refusing a missing value in a column that has never held
  one (dec-B322). Its column form borrowed the refusal that goes, for factors only: ran (05-test-column-na.R),
  a numeric test column so updated is taken today and the next run fits the row, where the whole test set and
  predict refuse. So the column form calls [`refuseTestMissingness`](../../R/data.R): one refusal for both.

## 5. Callers

- Package (read): `setState`, `getPointer` and `copy` share [`installStateOnto`](../../R/dbarts.R); a re-created
  engine takes its flags from the data object, so a reload or copy of a state stored before the value draws.
  benchmarks (ran: grep; read): sbc.R, negbin-mixing.R and C1-frozen-ess.R install a sampler's own state or a
  same-data sampler's. The warm start ([`Sampler::installForests`](../../src/bartcore/sampler.hpp)) is O1.
- stan4bart (bartcore f0084aa, read): restoreBartSampler builds a sampler per chain from the fit's control,
  model and data and installs the state that chain's sampler stored, kept whole; no predictor call in its R/ or
  in dbarts.h (ran: grep), so S equals F; saved fits carry no record and install as today. bartCause (e833be7),
  treatSens (babfaa6): no state or predictor call; bairrtt (3f57f61): no state call (ran: git grep).

## 6. Equivalence, snapshots and existing tests

Expected: nothing recorded moves (ran: grep). The equivalence trio installs no state, reads no `$state` and
gives no factor a missing label after creation; the four reproducibility files hold none of setState,
setPredictor, setData, copy, readRDS, installTrees, warm.start, `$state` or NA; of the exact gates' scripts only
negbin-mixing.R installs a state, its own. Existing tinytest (ran: 03-suite-probe.R, a wrap of
`installStateOnto` over the 61 files that install a state and change a predictor or write NA; 522 installs):
2 put a state with no right-going direction onto a flagged column, each the sampler's own or a same-data
twin's (read). tests/cpp: 7 files call `setState` and mention missing values; not probed. A test that pins the
defect itself (all left after an older state; a kept exempt rule flipped) is rewritten and named in the landing
note; any other move is a stop.

## 7. Steps, tests and mutants

Steps: (A1) engine: the field, `getState`, the exempt rule, the completion; (A2) bridge: the attribute written
and read, the registry comment; (A3) tests/cpp, tinytest, mutants, docs; (B) the refusal, its tests and help.

Step A, tinytest, new test-state-missing-record.R on the fixture of
[test-missingness-first-seen.R](../../inst/tinytest/test-missingness-first-seen.R), seed-held preconditions asserted:
- the record: present and all FALSE on a complete sampler; TRUE for x1 after a first missing value by each
  path, and after the column is filled; kept by saveRDS; of another length or type, an error.
- the review's case with the update UNFORCED and asserted TRUE (no tree merged), on one forest and on two:
  `setState(old)` is TRUE and the directions on x1, live and kept, the generators and predict on a missing row
  are identical to those taken between the update and the install; a second such update draws nothing.
- a wrapped ring: asserted a cursor other than 0 and rules on x1 in kept slots on both sides of the seam, then
  the same equality; 20 kept into a store of 5 (Part B's B4) keeps the newest 5 as the update drew them.
- copy and reload: asserted first that the cached state's record is FALSE for x1 and the data object's TRUE;
  each then holds the source's directions, their next 5 draws equal within 1e-12.
- nothing to draw: a state stored after the value and 10 sweeps keeps its learned directions and generator; one
  with the attribute removed installs as landed; directions on x1 into a complete twin give FALSE, record NULL.
- a factor: a complete twin's state into a sampler created with missing values holds both directions where one
  can reach and none on a hand-built rule out of reach. The critique's case: kept rules with no category on the
  right (asserted present) keep their direction through a first missing value and through an install that
  draws, and the sampler reinstalls, copies and reloads its own state.
- monotone, and not clean: a stored state whose drawn directions break the order of a monotone tree, and a hand
  state whose drawn direction empties a leaf (Part A's lone-left tree), each coin asserted: FALSE unforced with
  getTrees, generators, `$state` and the next 20 draws a twin's; NULL forced, reseeded or merged; else TRUE.
- completes to invalid: a hand state, record FALSE, with a live exempt rule beneath a rule on its column: on the
  coin that takes the missing value away (asserted), an error, the sampler a twin's; on the other, installed.

Step A, tests/cpp, `testStateMissingRecord` in [test_state.cpp](../../tests/cpp/test_state.cpp), by the pattern
of [`testMissingFirstSeen`](../../tests/cpp/test_sampler.cpp): `getState` writes the flags; state, one
whole-column first missing value, state again, install the first: trees, kept draws and generator bytes equal
the second's (two chains, two forests, a variance forest, a pooled column, a wrapped ring); exempt rules, nested
and out of reach, take no coin in either caller; cleared `rngState` draws from the chain's own; a throw leaves
the generators; the fuzz states with a column cleared from the record complete to valid ones.

Step A mutants, each failing a test: the record not written, written from content, or through getMissingSeen;
the raised set reversed; no record read as none missable; no draw on live trees, kept draws, the variance
forest, a second forest or later chains; coins from the chain's generator; the state's generator installed
without the advance; kept draws in slot order or after the repack; a draw at every factor rule; a coin for an
exempt rule; validation, verdict or cone run on the state as stored; a generator not put back; a flag raised.

Step B tests. [test-monotone-unforced.R](../../inst/tinytest/test-monotone-unforced.R), 12 pins: the section
returns to what the engine does on the tree a missing f breaks, each outcome checked against the side an
unconstrained twin on the same seed draws: FALSE and untouched, the row declined, NULL and reseeded when forced,
TRUE where the drawn side keeps the order. [test-joint-update-factor.R](../../inst/tinytest/test-joint-update-factor.R),
3 pins: every sampler takes a missing label; the third, complete in f, draws and the other two do not.
[test-data-categorical-declared.R](../../inst/tinytest/test-data-categorical-declared.R), 1 pin, setTestPredictor's:
still refused, in the test path's words, with a numeric twin. test-missingness-first-seen.R runs its paths for
f too. Mutants: the refusal left on one path; a missing label coded as a level; the test setter taking one.

Docs. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the `state` field gains "and which predictor
columns could hold a missing value when it was stored"; Missing values in predictors gains "Installing a stored
state, by \code{setState}, \code{copy} or a reload, into a sampler where a column can hold missing values that
could not when the state was stored draws the side of every rule on that column, in the state's trees and saved
draws, as a column's first missing value does. \code{setState} counts such an install as clean and returns
\code{TRUE}, unless a drawn side leaves a leaf with no row or takes a monotone tree out of order." Step B: the
refusal of "a missing label in a column that has never held a missing value" becomes setTestPredictor's alone,
and the paragraph beginning "A column coded from a factor" goes. [NEWS.Rd](../../inst/NEWS.Rd) states no refusal
(ran: grep); its first-missing-value item gains a clause. The design note is amended; both TODO items close.

## Verification

Independently of the implementer, `R_LIBS=$LIB` on every R call
([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting, nothing to re-record expected):
`R CMD INSTALL --preclean --library=$LIB .` per engine commit; `cd tests/cpp && make && ./test_bartcore`, plain
and with `OPT="-O2 -g -fsanitize=address,undefined"`; `tinytest::run_test_file` on each touched file (also under
R-loaded ASAN), then `tinytest::test_package("dbarts", at_home = TRUE)`; on a reference build, `Rscript
benchmarks/R/equivalence.R compare <MANIFEST file> --bitwise --strict-coverage`, the bcf and multinomial
compares and the four `test-reproducibility-*.R` files: "identical draws (same RNG stream)" throughout; the exact
gates with `quick`; `R CMD check --as-cran`; lint, air, the doc checks; stan4bart's suite and a saved fit reloaded.

## What touches Part B

Against [Part B: setState's two forms (dec-B305, dec-B310, dec-B318, dec-B398, dec-B401)](setstate-force-update.md#part-b-setstates-two-forms-dec-b305-dec-b310-dec-b318-dec-b398-dec-b401).
Recheck when it lands, most likely to be wrong first (the critique's order):
1. The seam. B2's staging routine may fold the verdict into `stateIsValid`'s loop. Section 2 asks only that
   the trees be walkable and nothing built, routed or judged before the draw; place it by that.
2. B1's signature: the refusal out-parameters go, `force` and `notClean` arrive, the monotone check feeds the
   verdict. Every reader after the draw takes the copy, `altered`'s successor too; recheck install-equals-state asserts.
3. B4's ring: the draw uses the state's capacity S, R = min(recorded, S) and slot (cursor + S - R + d) mod S,
   before the repack. Recheck that nothing reduces the cursor against the recipient's capacity first.
4. B6's R side: section 3 relies on `$state` being the object given after TRUE and after a forced install. If
   Part B stores again after a forced repair, that bullet and the copy-and-reload test change.
5. B2's "generators untouched on FALSE": the borrowed generator goes back on every exit (the scope guard).
6. Part B's fuzz check (unforced false exactly when forced reports a repair) runs on the completed copy.
Also: the bridge reader gains the attribute beside B6's `force`; B9's text takes this slice's sentences.

## Open call

O1. Should a warm start draw the directions for missing values?

A warm start (`bart(warm.start = fit)`, or `sampler$installTrees(fit)`) copies an earlier fit's trees into a
new sampler as its starting point; the new sampler then runs its own sweeps.

The case: the earlier fit's data had no missing value in a column, and the new data have some. The copied
trees were grown when nothing there could be missing, so their splits on that column say nothing about where a
missing value goes; today every one sends it left. 0.9-34 has no equivalent: no warm start (no `warm.start`
argument, no `installTrees`), and rows with a missing predictor value dropped (147 of 150 kept; run 2026-10-10).

Measured on this branch, 2026-10-10: 150 rows, 3 of them missing in x1, 2 chains of 15 trees, warm-started
from a fit to the complete data. Before any sweep, 7 of 7 splits on x1 send missing values left. After 1 sweep,
8 left and 1 right. After 31 sweeps, 8 left and 4 right.

- (a) Draw. Each chain of the new sampler draws each such split's direction at one half, from its own
  generator, before the warm start merges the splits the new rows leave empty. It rides this slice: about 110
  lines (engine 30, bridge 15, tests 60, help 5), same review. It reads the donor's record of which columns
  could be missing, which this slice adds; a donor state stored before it has none and installs as today, all
  left. Such a warm start's draws shift; no recorded baseline or gate warm-starts (grep).
- (b) Leave it. The splits start left and the sampler's moves change them over the run. About 3 lines of help.

Recommended: (a). By dec-B321's reasoning a direction no missing row has reached has its prior, one half, as
its conditional, so all left is not a draw from anything the earlier fit says. And the left start lasts: 8 of
12 splits still left after 31 sweeps. Neither option is wrong in the long run: a warm start is a starting
point. What would change it: a left start gone within a few sweeps, or a view that such warm starts will not
occur. How realistic: a refit seeded from an earlier fit after rows with gaps arrive; uncommon.

## Calls made

- An install that drew is clean: unforced TRUE. dec-B435 names the draw, not the return value. The grounds:
  dec-B305's clean is "leaving a sampler in a coherent state with regards to its model", its not clean "a tree
  that would have to be changed to fit"; dec-B310's unforced value marks "a statistically valid update"; a
  direction no row had reached has its prior as its conditional (dec-B321); an unforced setPredictor that draws
  returns TRUE (ran). What a user sees: after such a TRUE the directions and generator are not the given
  state's, and `$state` is still the object given. Alternative: FALSE, the draw under force.
- The draw completes a copy before validation, verdict and install (alternative: install, then draw and
  refresh, merging and judging under the left default). Coins from the state's generator (alternative: the
  sampler's current one, which the install overwrites; a reload of `$state` would not reproduce its source).
- The exempt rule covers "no category on the right" only (alternative: the mirror too, about 35 lines for the
  reachable categories in the flat walk, to keep a rule that stays well formed either way).
- No record draws nothing; a state never raises a flag; the attribute is `missing.columns`, always written
  (alternative: the data slot's name, NULL meaning two things); every recorded kept draw of the state draws,
  before the repack (alternative: only those the store keeps); `$state` after `setState` stays the object given.
- Step B: setTestPredictor's column form refuses through the test path's own check, so also the numeric test
  column it takes today (dec-B322; alternative: keep the factor refusal for that caller, a TODO for the rest).

## Budget, stops and evidence

Planned ~805. Step A ~620: engine 135, bridge 30, tinytest 230, tests/cpp 180, help, NEWS and docs 45. Step B
~185: R 25, tests 150 (the monotone file about 100), help 10. O1 as (a) adds ~110. Forecast at 1.5 to 2 times,
as this surface has run: 1210 to 1610; stop at ~2110 (~2400 with O1). Stop and report, without working around,
when: the slice passes its stop; anything recorded or an existing test other than section 6's exception moves;
section 3's scoped loop does not reproduce the update; a `getState` state completes to an invalid one; Part B
landed so that section 2's seam cannot be placed, or with another ring rule; step B needs an engine change.

Evidence. Ran on c4ea26ab, library scratch/libs/smr, in scratch/smr/: 01 (the review's case, copy, reload,
twin), 02 (warm start), 03 (install probe), 04 (factor paths), 05 (test setter), 06 (0.9-34, library
scratch/libs/cran), crit/ (the critique's p1 and p2 probes rerun, identical output; the mirror); the greps. Read:
every claim cited by symbol; Part B as planned. Not run, no build of the slice existing: the install path, the
exempt rule, every mutant; tests/cpp. Revised from the critique (7dcdb6f1, scratch/smrcrit/critique.md): M1 and
S2 (the exempt rule; one validation, on the copy), M2 (scope), M3 (tests, seeds), S1, S3, S4, S5, step B.
