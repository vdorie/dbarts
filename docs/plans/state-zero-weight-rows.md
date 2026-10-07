# state-zero-weight-rows: a Student-t state names its zero-weight rows

Status: PLANNED (dec-B277).

agent: opus implementer, one (engine, facade, bridge); opus reviewer.
rng: by call sequence.
- SHIFTING for one install: a Student-t state that carries the new record, installed by `setState`, `copy`
  or a reload under weights other than the ones it was stored under. Each chain then draws one gamma variate
  per row that enters the likelihood where today it draws one per row in the likelihood, so every later draw
  differs; both are draws from the scales' own conditionals and the posterior of the run that follows is the
  same. Measured per chain, 81 rows, 20 at zero when stored: 61 today against 0 when the zero rows are the
  same, 61 against 20 when 20 rows enter and 20 leave, 81 against 20 under all-positive weights, 51 against
  0 when more rows go to zero.
- POSTERIOR-CHANGING, as a correction, for code that undoes a weight change by `setState` first and
  `setWeights` second, sweep after sweep: today the rows positive under both vectors are left, until the
  next sweep, with scales drawn under the proposed weights; afterwards they hold the stored ones.
- NEUTRAL, bit for bit, for everything else: every other family, an install under the state's own weights,
  a state without the record, `setWeights`, `setActiveRows`, and every sampler that installs no state.
Proved by the bitwise gates below (no baseline scenario, exact gate or snapshot file installs a Student-t
state under other weights; checked by search at the tip), by a seeded digest of unaffected fits and restores
on the base and slice builds, and by the new tests' identity with the same `setWeights` call.
window: pre-release, before 1.0-0 (dec-B277). Serial with [cut-points-undo.md](cut-points-undo.md) and
[leaf-conversions.md](leaf-conversions.md), which edit the same bridge file and manual page, the second
also chain.hpp, model.hpp, sampler.hpp and the state reader and writer. Recommended order: cut-points-undo
(landed), then this, then leaf-conversions rebased onto it; see Calls. This slice changes one facade
virtual and adds one, so every worktree that rebases over it reinstalls with `--preclean`.
budget: ~480 lines (C++ engine and facade ~55, bridge ~45, R none, tests/cpp ~110, tinytest ~220, manual,
design notes, comments and TODO ~50), upper figure 720. Plans have run 1.5-2x low; the design this comes
from estimated 150, tests not counted.

## Goal

A Student-t state records which rows were at weight zero when it was stored. Installed under other weights,
by `setState`, a copy or a reload, it has the scale of exactly the rows that enter the likelihood redrawn
(at zero then, positive now), as the same `setWeights` call would, and every other row keeps its stored
scale. Code that runs the sampler inside a larger sampler can then undo a weight change with `setState`
and `setWeights` in either order and be left with scales from the right conditionals.

## Context

All numbers were run on the tip's build (shipped mode): 81 rows, 20 trees, weights `rep(c(1, 4, 4))` and
`rep(c(2, 2, 5))` with blocks of 20 rows set to zero, one and two chains, the degrees of freedom fixed at 5
and estimated. The counts are per chain and were the same in all four arms.

- The scales and their record. [`TResponse`](../../src/bartcore/model.hpp) holds a scale per row and a mark
  per row saying whether the row was in the likelihood (weight times mask positive) when last marked.
  [`TResponse::redrawEntering`](../../src/bartcore/model.hpp) redraws, in row order from the chain's own
  generator, each row that is in now and was marked out, and marks every row afresh;
  [`TResponse::setWeights`](../../src/bartcore/model.hpp) and
  [`TResponse::setActiveRows`](../../src/bartcore/model.hpp) call it. A row at weight zero still gets a
  scale every sweep ([`TResponse::refreshLatents`](../../src/bartcore/model.hpp)), drawn without its
  residual. [`TResponse::restoreLatents`](../../src/bartcore/model.hpp) copies stored scales in and leaves
  the marks alone.
- What a state says about weights. The bridge's [`storeState`](../../src/R_interface_bartcore.cpp) writes
  ["weights.digest"](../../src/R_interface_bartcore.cpp), 8 raw bytes at the top level, from
  [`Chain::weightsDigest`](../../src/bartcore/chain.hpp); the weights are not stored. The top-level
  attributes today are `cutPoints`, `currentSampleNum`, `recordedDraws`, `formatVersion`, `weights.digest`,
  `survival.digest`, `packageVersion` and the class; the state of the fixture is 11672 bytes.
- The install today. The bridge's [`setState`](../../src/R_interface_bartcore.cpp) installs the state and
  then, when the stored digest differs from the destination's, calls
  [`Sampler::reapplyWeights`](../../src/bartcore/sampler.hpp), each chain in turn
  ([`Chain::reapplyWeights`](../../src/bartcore/chain.hpp)).
  [`TResponse::reapplyWeights`](../../src/bartcore/model.hpp) forgets every mark and so redraws every row
  in the likelihood. Measured against a twin that takes the same state under the stored weights (nothing
  drawn) and then the same `setWeights` call:

  | destination weights | rows in the likelihood | rows entering | restore redraws | `setWeights` twin redraws |
  |---|---|---|---|---|
  | the stored ones | 61 | 0 | 0 | 0 |
  | other positive values, same zero rows | 61 | 0 | 61 | 0 |
  | 20 rows enter, 20 other rows leave | 61 | 20 | 61 | 20 |
  | all positive | 81 | 20 | 81 | 20 |
  | 10 more rows at zero | 51 | 0 | 51 | 0 |

  The twin redraws exactly the entering rows in every arm; the restore leaves the stored scale at every row
  at zero now; `setState` returns `TRUE` throughout. Two restores of one state under the same other weights
  are the same sampler, scales and next draws identical.
- Under the state's own weights the install draws nothing: the scales are the stored ones bit for bit, the
  next five draws are identical to a restore of the same state with the digest removed, and they differ
  from the uninterrupted source's by 4.8e-15 (the fit) and 3.3e-16 (sigma), the size any restore leaves.
- A state with no digest (one written before the digest existed) redraws nothing under other weights
  (0 of 81), as [test-state-weight-pairing.R](../../inst/tinytest/test-state-weight-pairing.R) pins.
- Undoing a weight change, one chain: store, `setWeights(proposed)`, three sweeps, then go back.
  - Old weights first, then `setState`: 0 scales differ from the stored ones and the next five draws are
    identical to a sampler that never proposed. Both cases below.
  - `setState` first, same zero rows: after `setState` 61 scales differ, at the end 61, and the next draws
    differ from the never-proposed sampler's by up to 2.1. Under the rule planned here (run at the tip as
    the twin's calls): 0 differ and the next draws are identical.
  - `setState` first, 20 rows entering and 20 leaving under the proposal: at the end 81 differ today (the
    41 rows positive under both among them). Under the planned rule 40 differ: the 20 that entered and
    left again, which are out of the likelihood, and the 20 that left and came back, which the final
    `setWeights` redraws. The 41 hold the stored scales.
- The law of what is left, `setState` first. Transformed by the conditional at the stored fit and sigma
  (uniform for a draw from it): today the rows positive under both vectors are uniform against the
  conditional under the proposed weights (KS p 0.55 and 0.56) and not under the old ones when the change is
  large (proposed 20 times the old: mean 0.184, n 3660, p 0); for the (1, 4, 4) to (2, 2, 5) change the
  difference is at the edge of what 12200 scales show (mean 0.5003, p 0.035). Under the planned rule the
  rows that left and came back are uniform against the old weights' conditional (n 1600, mean 0.509,
  p 0.27).
- A copy and a reload made after the weights changed without a store (state stored, then
  `setWeights(w, updateState = FALSE)`): the live sampler redrew the 20 entering rows; the copy and the
  reload differ from the stored scales at 61 rows and from the live sampler at 53 (one chain) and 53 and 51
  (two), and their next three draws differ from the live sampler's by 1.4 and 1.9. The same calls on a
  fresh sampler that restores under the stored weights and then takes the `setWeights` agree with the live
  sampler to 1e-15 in the scales and 7e-15 in the next five draws, which is what the copy must do.
- Under a mask. A live `setState` under a mask and other weights redraws the rows positive now and active
  (49 of 49, no masked row). A copy and a reload get the mask back after the state
  ([`reapplyActiveRows`](../../R/dbarts.R), called by [`getPointer`](../../R/dbarts.R), `setState` and
  `copy`), so they redraw every row at positive weight, the 12 masked ones included. The same `setWeights`
  call under the mask redraws the entering rows that are active (14 of 20); the 6 masked ones are redrawn
  when the mask lifts, with every other masked row at positive weight (12 in all).
- The mask's own undo has the shape the weights will have: `setState` first and the old mask second
  redraws the 20 rows the proposal had masked; the old mask first, 0. `setState` returns `TRUE` in both.
- Refusals that exist. A state from a sampler with another row count: `state is not consistent with this
  sampler`. A gaussian state into a Student-t sampler and the reverse: the same message. A Student-t
  sampler refuses a variance forest. `setWeights(NULL)` is refused, so weights are never removed.
- An unweighted Student-t sampler's digest equals that of all-ones weights. Its state installed under
  weights with 20 zeros redraws 61 scales today although no row enters.
- The reader ignores a top-level attribute it does not know: a state given an extra one installs
  identically on the tip, scales and next draws, and `setState` keeps it on the `state` field.
- Logistic. A count swap redraws every active row
  ([`LogisticResponse::setWeights`](../../src/bartcore/model.hpp)), and so does a restore under other
  counts: with 10 of 81 counts changed both redraw 81 and the two samplers are identical. `setWeights`
  refuses a zero count. There is no difference between the two routes for a record to remove.
- Consumers (read only). stan4bart stores `sampler$state` on its fit and puts it back with `setState`; its
  samplers are gaussian or probit, which hold no such scales, and its tests check that the state is there,
  not its shape. bartCause, treatSens and bairrtt store no state. The flat C header has no
  entry for states or weights.

## The rule

What is stored. Every Student-t sampler writes, beside the digest, a top-level attribute `weights.zero`: a
raw vector with one byte per row, 1 where the weight in force is zero and 0 elsewhere, all 0 when the
sampler has no weights. No other family writes it. The mask is not in it.

The install, Z the recorded rows, w and a the destination's weights and mask in force at that moment:

1. Stored digest equal to the destination's: nothing is drawn and the marks are not touched. As today.
2. Digests differ and the state has the record: after the whole state is in, each chain in chain order
   draws, from its own generator as the state restored it and for rows in increasing order, a new scale for
   each row with Z, w positive and a active, from Gamma((nu + 1) / 2, rate (nu + w r^2 / sigma^2) / 2) at
   the installed fit, sigma and degrees of freedom. Every other row keeps the stored scale. Every row is
   then marked by whether it is in the likelihood here.
   - No row enters (the same zero rows, or more of them): no random number is used.
   - A row positive then and zero now keeps its stored scale; a later `setWeights` that brings it in
     redraws it.
   - An entering row that is masked is not drawn; `setActiveRows` redraws it when the mask brings it in.
   - On a copy, a reload, and `setState` on a sampler whose engine has to be made again, the mask goes
     back after the state, so no mask is in force at the install and an entering masked row is drawn then.
3. Digests differ and the state has no record: today's rule, every row in the likelihood here redrawn.
4. No digest: nothing is drawn, record or not. As today.
5. A record that is not a raw vector of 0s and 1s is refused with `malformed zero-weight rows in bartcore
   state`; one whose length is not the sampler's row count with today's `state is not consistent with this
   sampler`. Both before anything is touched.
6. A record put on another family's state is checked the same way and otherwise unused.

`setState` returns what it returns today: `TRUE` promises the trees, not the scales (dec-B234).

## Constraints

- `setWeights`, `setActiveRows` and `setData` draw exactly what they draw now.
- A state of any other family is byte for byte what it is now, so nothing changes for a consumer that
  stores gaussian or probit states.
- The mask stays out of the state (dec-A116).
- State format: an added top-level attribute. `formatVersion` and the oldest readable version stay 1,
  under the rule at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp).
- A state stored before this installs as it does today (rule 3 and rule 4), which the tests hold by
  removing the attribute. A state stored after this goes into an earlier build of this branch as a state
  without the record (measured: the attribute is ignored); nothing is promised about that before the first
  release, and nothing breaks.
- The flat C header is not touched and [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not
  move. facade.hpp changes ([`SamplerBase::reapplyWeights`](../../src/bartcore/facade.hpp) gains an
  argument, one reader is added): `--preclean` on every install, and the table of virtuals in
  [`FacadeVirtual`](../../tests/cpp/test_facade.cpp) gains its row.
- A warm start reads no scale and no digest and is not changed.
- No R code changes, and no NEWS entry: Student-t fits, the digest and the redraw are all new in 1.0-0.
- Out of scope: the logistic family; drawing nothing for a row that leaves the likelihood and comes back
  with no sweep between; putting the mask back before the state on a copy or a reload. See the end.

## Steps

1. Engine. [`ResponseModel::reapplyWeights`](../../src/bartcore/model.hpp) takes the record, a byte per
   row or null. [`TResponse::reapplyWeights`](../../src/bartcore/model.hpp): with null, as today; with a
   record, mark each row in exactly when the record has it at positive weight, then proceed as `setWeights`
   does, which draws a row only when it is in the likelihood now, mask included, and was marked out. A
   reader on the response, the chain, the sampler and [`SamplerBase`](../../src/bartcore/facade.hpp) fills
   a caller's buffer with the flags and returns whether the family keeps such scales (Student-t: yes,
   flags all 0 with no weights; every other: no).
   [`Chain::reapplyWeights`](../../src/bartcore/chain.hpp) and
   [`Sampler::reapplyWeights`](../../src/bartcore/sampler.hpp) pass the record through, defaulting to null
   so the engine's own callers compile unchanged. tests/cpp, in the Student-t block of
   [`testMembershipAcrossForests`](../../tests/cpp/test_sampler.cpp), which already holds a state stored
   under one zero pattern and a destination under another with a mask:
   - the source reports flags equal to its zero weights; an unweighted Student-t sampler all 0; a gaussian
     sampler none;
   - install, re-derive with the record: the scales are the stored ones except at rows zero then, positive
     now and active, which are the conditional's draws in row order from a generator cloned before the
     call; the generators agree after; the precisions the trees read are weight times scale times mask;
   - then the mask lifted redraws every masked row at positive weight, and weights that bring the
     remaining rows in redraw those (the existing two checks, after the new call);
   - a record under which no row enters: scales identical to the stored ones, generator untouched;
   - the existing check without a record stays: every scale in the likelihood redrawn.
   [`FacadeVirtual`](../../tests/cpp/test_facade.cpp): the reader's row (the boundary reports the
   implementation's flags) and the changed call.
2. Bridge. [`storeState`](../../src/R_interface_bartcore.cpp) allocates the raw vector, has the sampler
   fill it and attaches it when the sampler says it keeps one.
   [`setState`](../../src/R_interface_bartcore.cpp) reads the attribute by name where it reads the digest,
   applies rule 5, and hands the record's own bytes to `reapplyWeights` when the digests differ (null when
   absent). The bytes are read in place from the R object: no copy is made, so nothing is held across a
   raised R error. The comment at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) names the
   attribute with the two digests.
3. tinytest, a new file `test-state-zero-weight-rows.R` on the fixture of
   [test-state-weight-pairing.R](../../inst/tinytest/test-state-weight-pairing.R); "fails today" names what
   the tip does.
   - What a state carries: raw, length n, equal to `as.raw(w == 0)`; all 0 for an unweighted sampler and
     for all-positive weights; stored under a mask it still names the zero-weight rows only; absent on a
     gaussian, a logistic and a probit state; `formatVersion` 1. Fails today: absent.
   - The identity, for the five destinations of the table, one and two chains, df fixed and estimated: the
     restored sampler's scales are identical to those of a twin that took the state under the stored
     weights and then `setWeights(destination)`; rows that did not enter hold the stored scales bit for
     bit; rows that entered differ from them; `setState` is `TRUE`; the next three draws are identical to
     the twin's. Fails today in four of the five (61, 61, 81 and 51 redrawn against 0, 20, 20 and 0).
   - An unweighted sampler's state under weights with zeros: nothing redrawn. Fails today: 61.
   - Under a mask on a live sampler: identical to the twin under the same mask; an entering masked row
     holds its stored scale and is redrawn when the mask lifts. Fails today. A masked sampler copied after
     a swap between weights with the same zero rows, made without a store, holds the stored scales at
     every row. Fails today: every row at positive weight is redrawn.
   - A copy and a reload after `setWeights(w, updateState = FALSE)`: scales within 1e-13 of the live
     sampler's and the next three draws within 1e-13, one and two chains. Fails today: 53 of 81 scales
     differ and the draws by more than 1.
   - The undo. Old weights first: the stored scales and the never-proposed sampler's draws (holds today).
     `setState` first, same zero rows: the same (fails today: 61 differ). `setState` first, zero rows
     moved: rows positive under both hold the stored scales (fails today: 41 differ) and the rows that left
     and came back do not.
   - A state from before: with the attribute removed, every row at positive weight and active differs from
     its stored scale (today's pin, moved here from the Student-t block of the pairing file); with the
     digest removed too, none does.
   - Refusals: a logical vector, a raw vector holding a 2, and a raw vector of another length, each by its
     message, the sampler's scales identical after. A well-formed record put by hand on a gaussian state
     changes nothing.
   - Two restores of one state under the same other weights are identical samplers.
   - The law: forty rounds of a sampler that stores under weights with twenty rows at zero and a second
     sampler at all-positive weights that restores the state; the entering rows' scales, transformed by
     the conditional at the stored fit and sigma, have mean within 0.05 of 0.5 and a KS p above 1e-6, the
     bounds of ["The redrawn Student-t scale is a draw from its conditional"](../../inst/tinytest/test-active-rows-reactivation.R).
     It passes today, when those rows are redrawn with all the others; it is there for the mutation that
     leaves a stored scale at an entering row (that test's comment gives a mean of 0.65 for a kept scale).
   In the pairing file the block
   ["Student-t: the scales are re-derived too"](../../inst/tinytest/test-state-weight-pairing.R) keeps its
   matched half and loses the two "every scale differs" pins; Student-t joins gaussian and the variance
   forest in ["neutrality: the repair is a measured no-op elsewhere"](../../inst/tinytest/test-state-weight-pairing.R),
   whose two weight vectors are positive (holds after, fails today); the header comment drops "Student-t
   redraws more".
4. Mutations (Verification): apply, install with `--preclean`, run, report the failing counts, revert,
   `touch`.
5. Records. The manual, [`dbartsSampler$setWeights`](../../man/dbartsSampler-class.Rd)'s `weights` item:
   the sentences from "Student-t redraws the scale of every row at positive weight and active" to "change
   them with `setWeights`" become: "A Student-t state also records which rows were at weight zero when it
   was stored, and a restore under other weights redraws the scale of exactly the rows that enter the
   likelihood, at weight zero then and positive now, as the same `setWeights` call would: each chain from
   its own restored generator, in row order. Every other row keeps its stored scale, and a restore between
   weight vectors with the same zero rows draws nothing. An entering row that is inactive under `active`
   waits for `setActiveRows`; on a copy or a reload, where the mask is put back after the state, it is
   redrawn at the install. A state stored before the record existed is installed as it was then, every row
   at positive weight and active redrawn." The Saving paragraph's last two sentences become: "A weight
   change is undone by `setWeights` with the old weights and `setState` with the stored state, in either
   order. With the weights put back first the sampler is the stored chain bit for bit. With `setState`
   first a Student-t sampler differs from it only at rows the change had taken to weight zero: their
   scales are redrawn from their conditional at the stored fit as the old weights bring them back, which
   is a valid draw and not the stored value. `setState` returns `TRUE` in both. A sampler saved or copied
   after its weights changed without a store is re-created from the state stored under the earlier
   weights, moved to the current ones as `setWeights` moved the live sampler; call `storeState()` after
   the change to have it reload or copy as it stands." The Value text for `setState` stands.
   Design record: [The saved state (follow-on)](../design/weighted-logistic.md#the-saved-state-follow-on),
   amended and dated, with the table above and the rule; the sentence in
   [active-rows-mask.md](../design/active-rows-mask.md) that says a state does not name its rows out; the
   list of added attributes in [public-surface.md](../design/public-surface.md). The comments on
   [`TResponse::reapplyWeights`](../../src/bartcore/model.hpp),
   [`Chain::reapplyWeights`](../../src/bartcore/chain.hpp),
   [`SamplerBase::reapplyWeights`](../../src/bartcore/facade.hpp) and after the install in the bridge.
   TODO: the entry `state-zero-weight-rows` names this plan.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library (facade virtuals change); `tests/cpp` builds and
  passes, clean under ASan and UBSan; the full tinytest suite; the new file and
  test-state-weight-pairing.R under ASan on the R-loaded path (the bridge reads R's bytes in place).
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares
  are bitwise, every scenario reporting identical draws, counted per scenario with no `max |z|` line: 55
  against `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded and no snapshot regenerated: the gaussian
  harness installs no state, the other two do not record a restore, and the snapshot files hold no
  Student-t fit. A scenario that is not identical is a finding, not a re-record.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged, `t-exact.R` and
  `mask-redraw-exact.R` among them. None installs a Student-t state at all. No new exact gate: the
  changed install is held to an identity with `setWeights`, whose redraw those tests and
  [`testMembershipAcrossForests`](../../tests/cpp/test_sampler.cpp) already hold to its conditional.
- One script on the base and slice builds digesting seeded scales and draws: each family through a restore
  of its own state, a copy and a reload, with and without a mask, one and two chains; gaussian, logistic
  and a variance forest under other weights; Student-t under its own weights; a Student-t state with no
  digest under other weights; and a Student-t state under other weights, its record removed on the slice
  build, against the base build's plain install. Equal.
- Mutations, each expected to fail the named test:
  - the response ignores the record (forgets every mark): tests/cpp "re-derive with the record" and
    tinytest "the identity";
  - the record is read inverted (rows positive then are taken as entering): the same two, and tinytest
    "rows that did not enter hold the stored scales";
  - the marks are set from the record and nothing is drawn: tests/cpp "re-derive with the record",
    tinytest "rows that entered differ from them" and "the law";
  - every row is marked in after the install: tests/cpp "the mask lifted redraws";
  - the bridge passes null although the state has a record: tinytest "the identity", "the undo" and "a
    copy and a reload";
  - only the first chain is given the record: tinytest "the identity" on two chains;
  - an unweighted Student-t sampler writes no record: tinytest "what a state carries" and "an unweighted
    sampler's state";
  - the writer flags rows out by weight or by mask: tinytest "stored under a mask it still names the
    zero-weight rows only" and "a masked sampler copied after a swap between weights with the same zero
    rows";
  - a state without the record draws nothing: tinytest "a state from before";
  - the reader skips the value check: tinytest "refusals".
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; `R CMD check
  --as-cran` on a tarball from a clean copy (man/ changes).
- Not a hot-path change: nothing is added to a sweep. A Student-t `storeState` writes n more bytes.

## Out of scope, and where it goes

- Logistic. Its restore under other counts already lands where the same `setWeights` call does; both
  redraw every active row, the rows whose count did not change included, and an undo is right in law in
  either order. Drawing fewer there needs the stored counts, not a flag. No TODO is opened; the
  coordinator has the text of one should that change.
- A row that leaves the likelihood and comes back with no sweep between still holds a scale drawn in the
  likelihood, and `setWeights` and `setActiveRows` redraw it all the same. Not redrawing it would make the
  undo bit for bit in either order; it changes what those two setters draw on sequences that exist today
  and is not what dec-B277 names. Not planned.
- Putting the mask back before the state on a copy or a reload, so an entering masked row waits there too.
  It moves every family's re-creation; the difference is one draw order in a corner (Student-t, a mask, a
  state stored under other weights, a row both masked and entering), right in law either way. Not planned.

## Calls made in planning

- The record is a raw vector, a byte per row, beside the digest. The engine reads R's bytes in place, so
  the bridge copies nothing and holds nothing across an R error, and the check is a length and a scan. The
  alternatives: the indices of the zero rows (smaller when few are at zero; needs a range and order check
  and a buffer for the engine) and a logical vector (four bytes a row). The state is opaque to users
  either way.
- Every Student-t sampler writes it, weights or not, and no other family. dec-B277 says "a Student-t fit
  that has weights"; with that, the state of an unweighted Student-t sampler would carry no record and,
  installed under weights with zeros, would still have every row redrawn (61 today where no row enters).
  Written by every family it would widen the gaussian and probit states consumers store for a record
  nothing reads.
- It records zero weights, not rows out by the mask. The mask is kept out of the state (dec-A116); a copy
  and a reload put the mask back after the state, so a record of masked rows would make a plain copy of a
  masked sampler draw; and a masked row's scale is not the defect, being drawn every sweep at the row's
  own weight with its residual and redrawn by `setActiveRows` when the row comes back. The alternative was
  the response's own mark, weight times mask.
- An absent record means "not known" and never "no rows": a state with a digest and no record installs by
  today's rule. That is what keeps a state stored before this installing as it does.
- The reconciliation stays where it is, in the bridge after the engine's install, with the record handed
  to `reapplyWeights`. The alternative, carrying digest and record on the engine's state structure and
  comparing inside the engine's `setState`, is tidier and about twice the change.
- A record of another length is reported as a state of another row count, with the message that refusal
  has today, so a Student-t state from a sampler of another size is refused in the words it is now.
- The tests get their own file; the pairing file keeps the digest's tests and loses the two Student-t pins
  that reverse.
- `setState`'s value does not change. A restore that redraws entering rows is still `TRUE`, as a restore
  that redrew every row is today.
- The manual keeps a sentence on the order of an undo, reworded: either order is right, weights first is
  bit for bit. The sentence that named one order as the only exact one goes.
- The `rng:` line names a posterior-changing sequence, the undo with `setState` first. It is a correction
  and is held to the identity with `setWeights` and the uniformity test, not to a new exact gate; the
  coordinator may prefer an arm in `mask-redraw-exact.R` that proposes and rejects weights (about 150
  lines more).
- Order with the other two plans. cut-points-undo, since landed, touches another function of the bridge
  file and another item of the manual page: first. leaf-conversions edits chain.hpp, model.hpp,
  sampler.hpp, the state reader and writer and the Saving text, in other functions and paragraphs, and is
  three times this size in two pushes: this slice goes before it, or between its pushes if it has already
  branched, never beside it, because a changed facade virtual under a rebase needs `--preclean` and the
  two edit the same reader.
- The tip against the texts of dec-B277 and dec-A160.
  - Both say an install under other weights redraws "every row in the likelihood". That is so for
    `setState` on a live sampler; a copy and a reload also redraw the masked rows at positive weight (12
    of 12), the mask going back afterwards. The manual already says so.
  - dec-B277 says a state stored before this has "every scale in the likelihood redrawn". That holds for
    a state with a digest; one without a digest has none redrawn, today and after.
  - Neither says what the other order of an undo leaves after the change. It is exact in law in every
    case and bit for bit only when no row went to zero and back (40 of 81 scales differ in the measured
    case, 20 of them out of the likelihood).
  - dec-A160's "for one sweep" was checked as a law at the tip: clear when the proposed weights are far
    from the old ones, at the edge of detection for a modest change.
  Nothing in this slice was found done already.
