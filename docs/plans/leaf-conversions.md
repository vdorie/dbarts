# leaf-conversions: kept draws and seeded forests keep their function when the units under them change

Status: PLANNED (dec-B200, dec-B231, dec-B233, dec-B237).

agent: opus implementer, one (engine, bridge, flat C entries); opus reviewer.
rng: three classes, by call.
- NEUTRAL for every fit and call not named below, and for every install between equal standardizations:
  draws bit for bit unchanged.
- SHIFTING for a linear-leaf sampler after `setData`, and for a linear-leaf warm start from a donor under
  another covariate standardization: the chain starts from converted coefficients, so its draws differ;
  the posterior does not.
- POSTERIOR-CHANGING on one sequence: a warm start from a donor on another cut grid into a linear or gp
  sampler whose leaf covariate was replaced after creation by `setPredictor`. Today that install re-derives
  the recipient's standardization (and default lengthscale) from its current predictors; afterwards the
  recipient keeps the one it had, which shapes its slope prior or kernel. And on one kind of fit: a linear or
  gp leaf covariate that holds a single value whose mean, taken through the sum, misses it by rounding (0.1
  over 150 rows; not 1000, 0 or 0.5). Today such a column gets that rounding as its scale (2.5e-16) and reads
  about 1 on every row; afterwards it is centred at its value, divided by 1 and reads zero, like any other
  constant column.
- Not a generator matter but a changed result: `predict` and every other replay of kept draws after a
  re-anchor or a `setData`, and five new refusals on a gp sampler that holds kept draws.
Proved by the bitwise gates below (no baseline scenario, exact gate or snapshot file pairs a leaf model
with `setData` or a warm start, or re-anchors with draws kept; checked by search at the tip), by a seeded
digest of unaffected fits on the base and slice builds, and by the new tests' function identities.
window: pre-release; any time. Serial with any other work in chain.hpp, model.hpp, sampler.hpp,
combiner.hpp or the bridge file. The `setState` return value ([restore-status.md](restore-status.md)) has
landed and is not touched. Steps 1 to 4 and step 5 may land as two pushes, each gated on its own
(Verification names each push's gates).
budget: ~1270 lines (C++ engine ~290, bridge and flat C entries ~95, R ~15, tests/cpp ~340, tinytest ~400,
manual, design note, NEWS and TODO ~130), upper figure 1900. Plans have run 1.5-2x low; the design this
comes from estimated 850.

## Goal

A draw the sampler has kept is the function that was drawn, and stays so: `predict` returns the same values
before and after `setData`, and before and after the response range is re-derived. A forest seeded from
another fit starts at that fit's function even when the two standardize a leaf covariate differently, and
the seeded sampler keeps its own standardization. A Gaussian-process sampler that holds kept draws, which
cannot be rewritten, refuses by name the calls that would silently change them.

## Context

Run on the tip's build; "kept draw" is what the messages and the manual call a saved draw. A linear leaf's
fit is an intercept plus slopes on covariates standardized by a centre and scale held by the leaf
([`LinearGaussianLeaf`](../../src/bartcore/model.hpp)); a gp leaf holds the same pair and a lengthscale
([`GPGaussianLeaf`](../../src/bartcore/model.hpp)). Every replay of a kept draw reads the leaf's current
values ([`Chain::addFlatPredictions`](../../src/bartcore/chain.hpp),
[`addFlatLinearPredictionsBelow`](../../src/bartcore/tree.hpp)).

Defects, measured (response sd about 2 throughout):
- Kept draws after `setData`. [`Chain::applyNewData`](../../src/bartcore/chain.hpp) re-derives the
  standardization ([`LinearGaussianLeaf::reinitialize`](../../src/bartcore/model.hpp)) and converts nothing.
  The same raw points, replayed before and after a `setData` whose covariate is 3 x + 5: constant leaf, no
  change; linear leaf, off by up to 7.9 (8.1 with two chains); gp leaf, off by up to 7.9 (8.0).
- The live linear chain after `setData`. With 40 rows appended inside every column's range (cut grid
  identical; centre 0.04 to 0.50, scale 0.99 to 1.23) the live fit on the old rows moves by up to 1.10;
  a constant leaf's moves by 4e-15. The gp leaf's fits restart at zero by design and are not part of this.
- Kept draws after a re-anchor. After `setResponse(3 y + 10, updateScale = TRUE)` kept draws come back as
  3 f + 10 (off by up to 20.7, equal to 3 f + 10 to 4e-15), on constant, linear and gp leaves alike, and
  the same through `setData` with that response. `setOffset(-5, updateScale = TRUE)`: off by 5. Student-t
  and monotone samplers as gaussian; aft by the log-time shift (0.5 for 0.5); nbinom by the log-mean shift
  (log 3 for tripled counts); probit not at all (it has no range). Under a variance forest the kept
  variance draws come back times 9 exactly. With `updateScale = FALSE` nothing moves.
- A warm start on the same grid reads the donor's slopes under the recipient's standardization. Donor centre
  and scale -0.019, 0.973, recipient 4.94, 2.92: the seeded fit differs from the donor's on the same rows
  by up to 6.30 (donor fit sd 2.07), from each of four kept draws by 6.2 to 6.4, and by 7e-15 when the
  standardizations are equal ([`Sampler::installForests`](../../src/bartcore/sampler.hpp) does not copy the
  donor's blocks; [`readWarmStartState`](../../src/R_interface_bartcore.cpp) does not read them).
- A warm start from another grid replaces the recipient's standardization:
  [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp) calls `reinitialize`. A recipient
  holding 4.94, 2.92 ends at -0.019, 0.973, on linear and gp leaves (the gp lengthscale is re-derived too).

What already works, and is kept:
- A state stored in other response units is converted on install, kept draws included
  ([`Chain::convertStateUnits`](../../src/bartcore/chain.hpp), dec-B200): intercepts and leaf values by
  ratio and shift, slopes and gp weights by ratio, variance factors by a power of the ratio. A gp leaf and
  forests with amplitudes refuse another shift. The re-anchor rewrite below is that arithmetic applied to
  the sampler's own store.
- After a re-anchor the live chain sits at the same internal values, so its fit in response units is
  3 f + 10 too. That is what re-anchoring a live chain means and it does not change.
- A gp warm start on the same grid copies the donor's per-row fits whatever the standardizations (5e-15).
- A covariate with no spread is standardized with a placeholder scale of 1: a constant column at 1000 gives
  centre 1000, scale 1, and slopes of sd 0.08 that no observation informs.
- Who reaches these paths: linear and gp leaves are taken by gaussian, aft, probit and nbinom samplers; all
  accept `setResponse(updateScale = TRUE)`; `setData` is refused on aft. A variance forest refuses a leaf
  model. A sampler with several forests refuses `updateScale = TRUE`, `setData` and a warm start, so every
  path here has one mean forest.
- On a gp sampler holding kept draws every call below is accepted today. `setResponse` to a response with
  the same midpoint and twice the spread is accepted and the kept draws come back doubled about it; onto
  the response already in force it changes nothing. With `updateScale = FALSE` (the default) `setResponse`
  and `setOffset` leave the kept draws identical.
- A gp sampler that holds no kept draw (trees not kept; kept but nothing run yet; kept draws dropped by a
  warm start) takes `updateScale = TRUE` on both calls: the live fit follows the new range (3 f + 10 to
  4e-15; an offset of -5 moves it by 5) and the next sweeps are finite. That stays exactly as it is.

## The rules

- Conversion of one leaf's coefficients from centre m, scale s to m', s', per covariate j:
  slope' = slope * s' / s, and the intercept gains slope * (m' - m) / s. The function is unchanged wherever
  the covariate is observed. A row missing covariate j is read at the centre in force, so it moves by that
  intercept term; this is stated in the manual, not hidden.
- No spread. A column that holds a single value is centred at that value exactly and divided by 1, so the
  leaf reads zero for it on every training row: its slope adds nothing to the fit and no observation
  informs it. Whether a column is one is read from what the leaf holds when asked (the scale is 1 and every
  gathered row is zero), never remembered, so a column given other values by `setPredictor`, or a state
  installed over other rows, is what it then is. Its scale is written to the state as `NA` beside a finite
  centre, and the engine still divides by 1. When a LIVE coefficient block is converted:
  - spread on neither side: the fit is the intercept before and after, so neither the intercept nor the
    slope moves;
  - spread on the new side only: the intercept stays (it is the fit) and the slope is set to zero: a slope no
    observation informed is not multiplied onto a real scale (0.13 drawn against the constant column above
    would become 40 on the scale of 307 that varying data gave it);
  - spread on the old side only: the intercept takes the function's value at the new centre, the one value
    the column now holds, and the slope is set to zero.
  KEPT blocks are converted by the formula with 1 for `NA` in every case: they are only replayed. A state
  written before this carries 1 and reads as 1.
- Equal standardizations, or a donor with no block: no arithmetic at all, so such installs stay bitwise.
- gp (dec-B237). A kept gp draw replays only under the centre, scale and lengthscale it was drawn with and
  has no mean term. While a gp sampler holds kept draws, `setData` is refused, and `setResponse` and
  `setOffset` with `updateScale = TRUE` are refused outright, whatever the new response or offset is and
  whatever the family, on the R methods and on both flat C entries. With `updateScale = FALSE` (the
  default) both are served as now. The refusal does not look at the values because one that did (refusing
  only when the new values move the range's midpoint) would let the same line of calling code pass on one
  sweep and stop on the next: a caller cannot predict it, and the case it would let through, a re-anchor
  that changes the spread and not the midpoint, is not one anybody relies on.

Messages, exact; `<caller>` is `$setResponse`, `$setOffset`, `dbarts_sampler_setResponse` or
`dbarts_sampler_setOffset`:

    $setData: 'newData' cannot replace the data of a sampler with gp leaves that holds saved draws: a saved gp draw replays only under the covariate standardization and response range it was drawn with; make a new sampler instead
    <caller>: 'updateScale' cannot be TRUE for a sampler with gp leaves that holds saved draws: a saved gp draw replays only under the response range it was drawn with; make a new sampler, or call without 'updateScale = TRUE'

## Constraints

- Every sampler on constant leaves, and every sampler that never calls `setData`, a re-anchor or a warm
  start, draws exactly what it draws now, but for the fit named under `rng:` whose leaf covariate is constant
  at a value its mean misses. On constant leaves the one thing that changes is what a replay of kept draws
  returns after a re-anchor.
- The recipient of a warm start keeps its centre, scale and lengthscale on either grid (dec-B231: derived
  data changes only by a call that names it). gp fits at a warm start stay as today: copied on the same
  grid, zero on another.
- `setState` is not changed: it still installs a state's standardization. Its move to comparing and
  converting, the state-side gp refusal and `setData`'s arguments are the later work (Out of scope).
- A refusal is raised before anything is touched: the sampler's state, its data object and its kept draws
  are identical after it.
- A gp sampler that holds no kept draw takes `setData` and `updateScale = TRUE` as it does today, and one
  that holds them takes `updateScale = FALSE` as it does today.
- facade.hpp is not edited: "gp leaves and a kept draw" is already on
  [`SamplerShape`](../../src/bartcore/facade.hpp) (`usesFunctionLeaves`, `numSavedDraws`), and no step adds
  a virtual. `--preclean` is still required on every engine commit, for the edited headers (data.hpp,
  model.hpp, chain.hpp, sampler.hpp), not for a changed virtual.
- The flat C header keeps every signature; its two re-anchoring entries gain a raised refusal, as the
  refusal on several forests is raised today. No consumer package uses gp leaves or re-anchors with draws
  kept (stan4bart re-anchors in warm-up with tree storage off; treatSens keeps no trees).
- State format: no new block and no version change. The scale block may now hold `NA`; a build before this
  refuses such a state, which the rule at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) allows
  before the first release.
- The two lines the mutation battery anchors in the warm-start reader
  (["malformed parameters in warm-start donor"](../../benchmarks/R/mutation-battery.R)) stay as they are.

## Steps

1. The conversion and the no-spread mark; no behaviour changes but one stored value.
   [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp) centres a constant column at its value;
   both leaf models report the mark from their gathered covariates
   ([`standardizedColumnHasSpread`](../../src/bartcore/model.hpp)), and
   [`LinearGaussianLeaf::restoreCalibration`](../../src/bartcore/model.hpp) (and the gp twin) reads `NA` as 1.
   One routine converts a coefficient block between two standardizations, with a switch for the live guard.
   [`Chain::getState`](../../src/bartcore/chain.hpp) writes a marked scale as NaN,
   [`Chain::leafCalibrationIsValid`](../../src/bartcore/chain.hpp) accepts it beside a finite centre, and the
   bridge writes and reads it as `NA`. tests/cpp, beside
   [`testLinearLeafFormats`](../../tests/cpp/test_model.cpp): there and back returns the block to rounding;
   function equality on rows with every covariate observed, two covariates; a row missing covariate j moves
   by slope_j (m'_j - m_j) / s_j and by nothing else; no spread on the old side, on the new side and on
   both, live (the three cases under The rules) and kept (formula with 1); a marked scale
   survives a state round trip and still divides by 1; a scale of 0, a negative one and an infinite one are
   still refused. tinytest: a constant leaf covariate's stored scale is `NA` (1 today), and the sampler
   restores, copies and reloads; ["leaf.covariate.scale <- 0"](../../inst/tinytest/test-mutate-then-serialize.R)
   still refuses.
2. `setData` converts, linear leaves. [`Chain::applyNewData`](../../src/bartcore/chain.hpp) reads the
   standardization before `reinitialize` and afterwards converts the live blocks (guard on) and every kept
   block in the forest's store. tests/cpp, beside
   [`testLinearLeafMutation`](../../tests/cpp/test_moves.cpp) (whose recovery check must still pass on its
   shifted draws): rows appended inside every column's range, so the grid and every route are unchanged:
   the live fit on the old rows is unchanged to 1e-12 (off by about 1 today), and kept draws replayed at
   fixed raw points are unchanged to 1e-12. tinytest, new file `test-leaf-conversions.R`: kept draws across
   a `setData` on a stretched covariate, one chain and two (off by about 8 today); a constant leaf beside
   it, identical. For the no-spread rule, in tests/cpp
   ([`testNoSpreadLiveConversion`](../../tests/cpp/test_state.cpp)) and in tinytest
   (["live coefficients across a covariate without spread"](../../inst/tinytest/test-leaf-conversions.R)),
   through `setData` and through a warm start: one constant replaced by another; a placeholder of zeros
   given values by `setPredictor`; a constant moved to another constant and then spread; values made
   constant; rows appended that give a constant column spread. Each holds the live fit on the old rows and
   on appended ones, and kept draws at fresh points, to 1e-10.
3. A re-anchor rewrites kept draws. The arithmetic of
   [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp) is factored so that it also runs over the
   chain's own store; `Chain::setResponse`, `Chain::setOffset` and
   [`Chain::applyNewData`](../../src/bartcore/chain.hpp) read the transform before and after and, when it
   moved, rewrite every kept mean draw and every kept variance factor. Live values are left. With nothing
   kept the call does what it does now. An install's own move of the transform (`setState`, a warm start,
   creation) is not such a site: the draws it brings are already converted, or there are none. On a gp leaf
   the shared routine carries the ratio and leaves the store alone when the shift moved, as the state
   conversion refuses it; step 5 refuses every such call from R and the flat entries, so only the engine's
   own tests reach it. tests/cpp, beside
   [`testVarianceSavedPredict`](../../tests/cpp/test_model.cpp): `predict` unchanged to 1e-12 across
   `setResponse(a y + b)`, `setOffset(c)` and `setData(x, a y + b)` on constant leaves (gaussian and
   nbinom), monotone and linear leaves, both channels of a variance forest, and, through the engine alone,
   a gp leaf at an equal shift (the shared routine's gp arm); out and back returns the store to rounding;
   the live fit after the re-anchor is a f + b, as today. tinytest: the same three calls on constant,
   monotone and linear leaves and a variance forest, aft and Student-t included, two chains once; no gp
   case here, step 5 refusing it. Each fails today by the numbers in Context. No existing test was found
   that reads a kept draw after a re-anchor (the four tinytest files that both keep trees and re-anchor read
   none afterwards); the suite decides, and one that does is rewritten and reported.
4. A warm start converts and stops re-deriving.
   [`readWarmStartState`](../../src/R_interface_bartcore.cpp) reads the donor's centre and scale blocks
   into [`ForestStateData`](../../src/bartcore/combiner.hpp);
   [`Sampler::installForests`](../../src/bartcore/sampler.hpp) converts the donor's coefficients, live or
   from a kept slot, into the recipient's standardization where the two differ (guard on);
   [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp) no longer calls `reinitialize` and
   keeps only the leaf's cache invalidation, the covariates no longer changing. tests/cpp, beside
   [`testCrossGridWarmStart`](../../tests/cpp/test_state.cpp): same grid, other standardization: the
   recipient's live fit equals the donor's function on the recipient's rows to 1e-12 (off by about 6
   today); a recipient grid that holds every donor split point and more: the same equality, and the
   recipient's centre, scale and lengthscale are bit for bit what they were (replaced today); equal
   standardizations: the installed slopes are identical to the donor's; gp: fits copied on the same grid,
   zero on another, kernel constants unchanged on both. tinytest: the three linear cases through
   `installTrees`, one through `bart(warm.start = )`, a donor's kept draw as the seed;
   ["warmLinear"](../../inst/tinytest/test-composition-sequences.R) stays identical.
5. The gp refusals. One guard beside
   [`refuseMultiForestResponseMutation`](../../src/R_interface_bartcore.cpp), called with it by
   [`bartcore_setResponse`](../../src/R_interface_bartcore.cpp),
   [`bartcore_setOffset`](../../src/R_interface_bartcore.cpp),
   [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) and
   [`dbarts_sampler_setOffset`](../../src/C_interface.cpp): gp leaves, a kept draw, and `updateScale` anything
   but an explicit FALSE raise the second message, before a value is read.
   [`bartcore_setData`](../../src/R_interface_bartcore.cpp) raises the first on gp leaves and a kept draw.
   Both read [`SamplerShape`](../../src/bartcore/facade.hpp) (`usesFunctionLeaves`, `numSavedDraws`); no
   file under src/bartcore changes in this step and tests/cpp gains nothing. tinytest:
   - refused, each by its message, with state, data object and kept draws identical after: `setData`;
     `setResponse` and `setOffset` with `updateScale = TRUE` onto another range, onto a response with the
     same midpoint and twice the spread (accepted today, kept draws doubled about it), and onto the response
     already in force (accepted today, a no-op); a probit gp sampler the same;
   - the two flat entries through ["capi_set_response"](../../inst/tinytest/test-capi.R) and its offset twin:
     an R error naming the entry when `updateScale` is true, the status 1 and identical kept draws when not;
   - still served, kept draws held, `updateScale = FALSE`: `setResponse(y + 1)`, `setOffset(rep(1, n))` and
     `setOffset(NULL)`, each leaving `predict` identical; `setPredictor`; `setCutPoints`;
   - still served, nothing kept, `updateScale = TRUE` and `setData`: trees not kept; trees kept and nothing
     run; kept draws dropped by a warm start. The live fit after `setResponse(3 y + 10, updateScale = TRUE)`
     is 3 f + 10, as today;
   - a linear sampler holding kept draws takes all three calls (steps 2 and 3).
6. Mutations (Verification): apply, install with `--preclean`, run, report the failing counts, revert, `touch`.
7. Records. [`dbartsSampler$setData`](../../man/dbartsSampler-class.Rd) and the `updateScale`,
   `installTrees` and Saving text: kept draws are rewritten across `setData` and a re-anchor; a row missing
   a leaf covariate is read at the new centre; a warm start converts and the sampler keeps its
   standardization; under `newData`, that a gp sampler holding saved draws refuses it; under `updateScale`:
   "`TRUE` is refused by a sampler with Gaussian-process leaves that holds saved draws, whatever the new
   response or offset: make a new sampler, or leave `updateScale` at `FALSE`." The two entries' comments in
   [dbarts.h](../../inst/include/dbarts/dbarts.h) (text only; the API hash does not move). A design note,
   docs/design/leaf-conversions.md, with its index row: the rules above, the measurements, why gp refuses;
   [linear-leaves.md](../design/linear-leaves.md), [gp-leaves.md](../design/gp-leaves.md) and
   [state-not-model.md](../design/state-not-model.md) point to it where they describe a data replacement or
   a warm start; check [feature-matrix.md](../design/feature-matrix.md) for a cell the refusals change.
   NEWS, user-visible changes, one item: saved draws keep their values when the response scale is
   re-derived. TODO: the two entries below.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library on every engine commit; `tests/cpp` builds and
  passes, clean under ASan and UBSan; the full tinytest suite; the new tinytest file and
  test-state-not-model.R under ASan on the R-loaded path.
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged; the three compares are
  bitwise, every scenario reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. That is the shifting class's statistical gate in its strongest
  form: nothing is re-recorded and no snapshot is regenerated. The two scenarios that call `setData` are on
  constant leaves with nothing kept. A scenario that is not identical is a finding, not a re-record.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged, `linear-exact.R`,
  `heteroscedastic-exact.R`, `aft-exact.R` and `negbin-exact.R` among them. None calls `setData` or a warm
  start, and each `setResponse` in them holds the range. No new exact gate: after the changed install the
  recipient is a linear or gp sampler under the standardization it reports, which those gates cover, and
  the conversions are held to an identity (same function) rather than to a distribution.
- One script on the base and slice builds digesting seeded draws: constant, linear and gp fits that call
  none of the three; a constant-leaf `setData`; a re-anchor with nothing kept; a constant-leaf warm start
  on each grid; a linear and a gp warm start between equal standardizations; a copy and a reload of each
  leaf model. Equal.
- Mutations, each expected to fail the named test:
  - the intercept's centre term is dropped (slopes scaled only): tests/cpp "function equality" of step 1,
    tinytest "kept draws across a `setData`";
  - kept blocks are not converted at `setData` (live only): tinytest "kept draws across a `setData`"; and
    the reverse, live blocks not converted: tests/cpp "the live fit on the old rows";
  - the live guard is removed: tests/cpp "no spread ... live";
  - the re-anchor rewrite drops the shift (ratio only): step 3's constant-leaf test, off by b; and skips
    the variance factors: step 3's variance channel, off by a squared;
  - the same-grid warm start skips the conversion: step 4's first test; the other-grid arm re-derives
    again: step 4's "bit for bit what they were";
  - the equal-standardization shortcut is removed: step 4's "identical to the donor's";
  - the guard ignores the kept-draw count: step 5's "nothing kept" calls; the guard also fires at
    `updateScale = FALSE`: step 5's "still served, kept draws held" calls; the guard is left off the offset
    entries: step 5's `setOffset` refusals, R and flat; `setData`'s check is removed: its refusal test.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")` is not needed (no new topic); inst/NEWS.Rd parses with a non-NULL result;
  `R CMD check --as-cran` on a tarball from a clean copy.
- As two pushes. Steps 1 to 4, with their mutations and their share of step 7: every gate above (the
  install, tests/cpp with sanitizers, the full suite, the reference-build files and compares, the exact
  gates, the digest, the lint chain, NEWS, `R CMD check`). Between the pushes a gp sampler holding kept
  draws is as it is today after `setData` or a re-anchor that moves the shift. Step 5, with its mutations,
  the manual and header text: NEUTRAL, no file under src/bartcore changing, so the `--preclean` install,
  tests/cpp unchanged and passing, the full tinytest suite with test-capi.R, the four reference-build files
  and the three bitwise compares (the bridge and the flat C file change), step 5's tests under ASan on the
  R-loaded path, the lint chain with `check-rc-codoc` and `check-doc-freshness`, and `R CMD check
  --as-cran`. The exact gates and the digest are not required for it.
- Cost: with nothing kept a re-anchor does one comparison more per chain. With draws kept it walks the
  store once per re-anchor; code that re-anchors every sweep and keeps trees pays that, and none of the
  consumer packages does. Not a sweep-path change; no bench compare.

## Out of scope, and where it goes

- TODO `state-frame-prior`, restated (the coordinator has the text). Left to it from here: `setData`'s
  arguments, after which the gp `setData` refusal narrows to the calls that do not hold both the
  standardization and the range; `setPredictor(updateStandardization = )` and its gp refusal; a state install
  that compares a stored standardization, converts linear slopes with step 1's routine and never installs
  it; the refusal of a state whose kept gp draws were made under another standardization; `setState`
  returning `FALSE` on such a difference; the record of all three on the data object.
- TODO `gp-kept-draw-kernel` (exists, dec-B237), two sentences added: "The refusals are step 5 of
  docs/plans/leaf-conversions.md. The kernel record lifts the standardization half of the `setData`
  refusal; the range half and the `updateScale = TRUE` refusals stand until a kept gp draw also records
  the response range it was drawn under."
- Not planned anywhere: converting the live chain at a re-anchor. The design keeps today's behaviour there.

## Calls made in planning

- The messages say "saved draws", the word the existing messages and the manual use, where the design says
  "kept".
- `updateScale = TRUE` is refused outright on a gp sampler that holds kept draws (the coordinator's
  decision, 2026-10-06), for the reason under The rules. Not taken: the design's "refuse only when the
  shift would move". It needs a preview of the transform a re-anchor would install (a virtual on the
  response models and the facade, about 80 lines with its tests) and makes the refusal depend on the data.
  The cost of the rule taken: a re-anchor that would have changed nothing a kept draw depends on, onto the
  response already in force or at an unchanged midpoint, is refused too.
- The refusal covers every family, probit included, where `updateScale = TRUE` re-anchors nothing (accepted
  today, kept draws identical). The refusal on several forests is family-blind in the same way, and a rule
  that names no family is the predictable one. The alternative is an exception for families with no range.
- Both surfaces share one message, so the flat entries' text says `updateScale = TRUE` of an int argument,
  as the existing several-forests message says `updateScale = FALSE` on both.
- gp fits at a warm start are left as the tip has them (copied on the same grid whatever the
  standardizations, zero on another grid); only the re-derivation of the recipient's constants stops. An
  earlier draft of the design zeroed them whenever the standardizations differ; its final text does not
  say, and per-row fits hold no kernel.
- The `rng:` line names a posterior-changing sequence the design filed under shifting: holding the
  recipient's standardization at a warm start from another grid changes the prior that sampler runs under
  when its predictors had been replaced. For a recipient unchanged since creation re-deriving gave the same
  numbers, so nothing moves there.
- `setData` on a gp sampler with kept draws is refused even when the new data would derive the same
  standardization and range (a fresh data object over the same rows moves kept draws by 0 today). The
  design's call, until `setData` can be told to hold both; the cost is that such a caller must make a new
  sampler or keep no draws.
- A marked scale is `NA` in the existing block, not a new block: one value says "no spread" to the reader
  that needs it. The cost is a state that an earlier development build refuses.
- NEWS carries one item, the re-anchor: released 0.9-34 replayed kept trees in the range in force (read
  from its source, not run). Leaf models are new in 1.0-0 and get none.
- The tip rechecked before building, 2026-10-07 (five slices had landed in the same files since the plan
  was written: cut-points-undo, state-zero-weight-rows, monotone-unforced-refusal,
  cross-family-state-install and the written-surface pushes). Every symbol the plan cites is still there
  under its name, and the six probes behind Context give the same output to the digit on a build of the
  tip. What moved is where things sit: `Chain::reapplyWeights` and the bridge's state reader and writer
  gained the zero-weight record, the state gained a `family` attribute beside the blocks step 1 touches,
  `Chain::stateIsValid` gained a check of the latents, and the predictor transaction no longer reseeds a
  monotone tree. None is in a function a step edits, and no rule, refusal or message changes.
- Rebased onto the tip after the review, 2026-10-07. Since the recheck above the state stopped recording a
  family, so the attribute named there is gone again; `setCutPoints` sorts a grid and refuses a repeated
  point; and the R front end evaluates a fit's rows once. None is in a function a step edits. One conflict, in
  [state-not-model.md](../design/state-not-model.md), where the warm start's sentence was added to a paragraph
  the tip had rewritten.
- Which kept draws a rewrite touches. The chain writes the ring of saved slots but only the sampler knows
  how many of them the runs have filled, so [`Chain::setResponse`](../../src/bartcore/chain.hpp),
  `Chain::setOffset` and `Chain::applyNewData` are handed the filled slots
  ([`SavedDrawSlots`](../../src/bartcore/chain.hpp)) and write only those. A slot nothing was recorded
  into stays the empty slot it was, which keeps "with nothing kept the call does what it does now" true of
  a sampler that keeps trees and has not run, down to the state it stores. The state conversion still
  walks every slot of a state, as before.
- The conversion is written as slope times (s' / s), and a column whose centre, scale and mark agree is
  skipped. At equal standardizations the arithmetic would therefore be exact anyway, apart from the sign
  of a zero intercept, so the mutation "the equal-standardization shortcut is removed" is caught only
  through a negative zero the engine tests plant (the conversion's own test and step 4's "identical to
  the donor's" in tests/cpp); the tinytest twin of that check cannot see it. The whole-forest comparison
  in the warm start is an early exit on top of the per-column one, and removing it alone changes nothing.
- Neither side with spread, at different centres (a constant leaf covariate replaced by another constant,
  through `setData` or a warm start): a live block does not move (the coordinator's ruling, 2026-10-07). The
  formula would add slope (m' - m) to an intercept that is the whole fit on both sides: with the constant
  going from 1000 to 2000 a live fit in (-6.1, 0.8) became (-1979, 8519). Kept draws keep the formula.
- The mark is read, not remembered. A flag set when the centre and scale were derived went stale when
  `setPredictor` gave a placeholder column of zeros its values: the guard then zeroed slopes the data had
  informed (the live fit on the old rows moved by 6.25 at `setData`, and a warm start from that sampler
  seeded a fit 6.25 off its own). The leaf now answers from its gathered covariates each time it is asked,
  which costs one pass over the rows per column still at the scale 1 and is not on a sweep's path. A column
  counts as without spread only where the scale is 1 as well as every row zero, so that a column held at its
  centre under a real scale keeps that scale in the state (an `NA` would restore as 1 and change what its
  kept draws predict off the centre) and converts by the formula.
- A constant column is centred at its value exactly. The mean taken through the sum misses 0.1 by rounding
  over 150 rows, which left such a column a scale of 2.5e-16 and a reading of about 1 on every row: on the
  base build its test rows at 0.11 are predicted at 1.8e14, and with the conversion built here the live fit
  after a `setData` that spreads the column ran to 2e16. This changes what such a fit draws without any of
  the three calls, the one departure from the first constraint; it is its own commit.
- [`SavedDrawSlots`](../../src/bartcore/chain.hpp) keeps its first slot. No run leaves a part-filled store
  whose draws start elsewhere, but a state whose write position is not its draw count installs one, the
  replay reads it from that slot on, and the rewrite must find the same slots;
  [`testNoSpreadLiveConversion`](../../tests/cpp/test_state.cpp) holds it.
- A donor's standardization record that does not fit the leaf (a length, a centre that is not finite, a
  scale that is neither positive and finite nor `NA`) is refused by the bridge with the existing
  "malformed parameters in warm-start donor"; the engine, reached without the bridge, answers with its
  shape refusal. A gp recipient does not read the record at all, so a gp warm start is exactly as before.
- A view built from a parent store carries the no-spread mark with the constants it inherits, so a fold
  sampler's leaf remembers it as a top-level one does.
- `bart(warm.start = )` runs a sweep before anything can be read, and one sweep at this size brings a
  cold start as close to the donor as a seeded one (measured), so "the first draw sits by the donor"
  cannot tell the builds apart. The test instead writes one donor forest under two standardizations, the
  second restated in R by the rule, and requires the two seeded fits to agree: to 1e-8 after `bart`'s
  sweeps and to 1e-12 through `installTrees` before any.
- Between the pushes the manual says what a gp sampler holding saved draws does today (its draws change
  at `setData` and at a re-anchor that moves the midpoint); step 5 replaces those sentences with the
  refusals. The two TODO entries are left for step 5 and the landing.
- The tip against the design. Its three probes were re-run on other fixtures: the facts hold, the sizes
  differ (7.9 where it reports 9.8 and 8.5; 6.30 where it reports 6.3). The design ran this work in series
  with the `setState` return value because both edit chain.hpp; that has landed, so the constraint is gone.
  Nothing in this slice was found done already.
