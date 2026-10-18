# response-scale-rows, slice C: setWeights and setActiveRows re-derive on request

Status: LANDED 2026-10-11 (46f2e2a2 to 6f5b9610; calls dec-A207; [Landing note](#landing-note)). Planned
against d5b637ee (bartcore with slice B, creation) and revised from its blind critique. Implements, for `setWeights` and `setActiveRows`: dec-B296 (the
setters that change the rows take the response setters' argument), dec-B327 (its effect on the live
chain), dec-B364 and dec-B405 (the keep, its status and its one warning) and dec-B363 (the by-hand call,
in dbartsSpec's help). Rests on dec-B362 (several forests refuse), dec-B330 and dec-B331 (kept draws),
dec-B314 (sigest is never re-estimated), dec-B381 (the argument's position), dec-B386 (a constant
response's window) and dec-B417 (a drawn k lags). Slice B is [response-scale-rows.md](response-scale-rows.md),
whose rule, constraints and keep this slice extends without restating; the standing statement is the
design note of the same name. Slice D is not planned ([Slice D](#slice-d)).

agent: one push. Opus implementer for the engine, the facade, the bridge and tests/cpp; sonnet for the R
methods, the tinytest file, the matrix runner, the help and the records; one opus reviewer told to refute,
who runs the mutants of [Verification](#verification) and reads the diff twice: once for bits that must
not move where the argument is absent, once for the order of refusal, install and re-derivation. Tier:
changes draws ([Process by risk](README.md#process-by-risk)), the engine and a facade virtual being
touched. Its blind critique ran on 2026-10-11 and is answered in the text.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them.
- NEUTRAL, bit for bit: every call that exists today. The two setters without the argument, FALSE by
  default; `dbarts` and `setData` on a count fit under 0/1 weights, whose follow-up call becomes the
  same engine operation through the new entry.
- NEW, and defined by an identity: `setWeights(w, updateScale = TRUE)` and
  `setActiveRows(a, updateScale = TRUE)` leave the sampler where today's `setWeights(w)` or
  `setActiveRows(a)`, then the bridge's offset entry called with the offset in force and `updateScale`
  true as retired: [`deriveCenterUnderMask`](../../R/bartcore.R) called it, then the record, leave it:
  bit for bit, the keep's status included. Not the R method `setOffset`, which also resets a test offset that
  follows the training one ([Context](#context)). No draw of an existing call changes.
- Moves no recorded scenario: none passes the argument, which does not exist, and none creates a count
  fit under 0/1 weights ([MANIFEST](../../benchmarks/baselines/MANIFEST), `equivalence-308c45be`'s row).
  `equivalence-308c45be` 55 of 55, `bcf-equivalence-1b7d730c` 15 of 15 and
  `multinomial-equivalence-80b1c8d4` 11 of 11 stay identical and the four snapshot files unchanged; no
  baseline is re-recorded and no exact arm is added ([Verification](#verification)).
window: after slice B has landed; before 1.0-0, and before the merge to main if it can be (TODO
`response-scale-rows`). Serial with any other work in the response and weight setters of
[`Sampler`](../../src/bartcore/sampler.hpp), in [`SamplerBase`](../../src/bartcore/facade.hpp)'s
virtuals, in [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) and
[`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp), or in the setter methods of
[dbarts.R](../../R/dbarts.R). Independent of [forest-sd-unit.md](forest-sd-unit.md): several forests
refuse the argument. TODO `sigma-prior-setter` comes with or after this slice.
budget: ~1120 lines planned: engine and facade ~40; bridge ~50; R ~75; tests/cpp ~200; tinytest ~510;
the matrix and its runner ~50; help and NEWS ~75; design note, indexes, amendments, TODO and ledger ~120.
Forecast 2000 to 2600; stop at ~2600, or at ~240 of package code ([Budget and stops](#budget-and-stops)).
No ABI event.

## Goal

A caller whose rows enter and leave the likelihood can ask, at the call that changes them, for the
response's scale and the count center to be taken again from the rows now in: `updateScale = TRUE` on
`setWeights` and `setActiveRows`, FALSE by default, with the response setters' effect on the live chain
and their keep. A count sampler made by hand from `dbartsSpec` becomes the one `dbarts()` makes by one
documented call, and a loop re-derives in burn-in and holds afterwards by one argument.

## Context

Marks: (ran) on 2026-10-11 on a shipped build of d5b637ee; (critic) run by the blind critique on that
build and not rerun; (0.9-34) on the released build; (read) in that source. Gaussian fixture unless a
line says otherwise: 200 rows, the last 100 out, their response 50 higher.

What the two setters do today.
- (read) [`Sampler::setWeights`](../../src/bartcore/sampler.hpp) and
  [`Sampler::setActiveRows`](../../src/bartcore/sampler.hpp) install in every chain and derive nothing.
  The bridge entries take two arguments and return NULL. The R methods are
  `setWeights(weights, updateState = NULL)` and `setActiveRows(active, updateState = NULL)`; on a
  probit, ordinal or count sampler `setWeights` forwards a 0/1 vector to `setActiveRows`
  ([`isMaskedWeightFamily`](../../R/spec.R)).
- (ran) Their latents, on a Student-t, a probit and a count sampler of 300 rows after 5 sweeps, by
  `setWeights` (0/1) and by `setActiveRows` alike: 100 rows leaving redraw nothing; 50 of them returning
  redraw those 50 and no other ([`TResponse::setWeights`](../../src/bartcore/model.hpp),
  [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp)); the same vector again, nothing; positive
  fractional weights on the Student-t sampler, no row entering or leaving, nothing. Then, on each,
  `setOffset(<offset in force>, updateScale = TRUE)` redraws 0 latents, and
  `setResponse(<response in force>, updateScale = TRUE)` redraws all 300 on Student-t and the 250 rows
  in on probit and count.
- (ran) Held: after `setWeights(w)` the range is creation's, shift 25.974 and scale 53.427. Positional
  `setWeights(w, TRUE)` stores the state and moves nothing.
- (0.9-34) `setWeights(weights, updateState = NA)` swapped the weights and held the scale: the sd of a
  prior draw of the fit was 14.3 at creation, 13.5 after `setWeights(w)` and 13.7 after
  `setOffset(0, updateScale = TRUE)` under w, where every row's range over 4 is 13.4 and the rows in
  give 0.95. `setWeights(w, updateScale = TRUE)` is an unused-argument error; `setActiveRows` reached no
  release.

What a re-derivation is today, by `setOffset(<offset in force>, updateScale = TRUE)` after the change.
- (ran) It reads the rows in: after `setWeights(w)` the call gives shift 1.153 and scale 3.785, the
  range [-0.739, 3.046] of rows 1 to 100 and the readers of the sampler created under w. Under weights
  and a mask together, gaussian and Student-t, it gives the literal over the rows both keep (center
  21.004, scale 2.949 on the critique's fixture).
- (ran, read) The live chain. The trees and the live leaves keep their internal numbers, so the live fit
  is read in the new transform until the next sweep (dec-B327). Kept draws of constant and linear leaves
  are rewritten and `predict` returns what it did (moved by 2.8e-14); a sampler with gp leaves that holds
  saved draws refuses (dec-B330; ran, the same sampler taking the call before a draw is saved). sigma
  holds its value in response units (moved by 0), its internal value and its prior restated; a variance
  forest is re-anchored ([`Chain::reanchorVarianceForest`](../../src/bartcore/chain.hpp)). A fixed k
  stays (2 before and after) while `k.scale` follows the range (26.71 to 1.89); on a second fixture a
  named sd of 0.7 has its k restated (24.20 to 2.64) and a drawn k stays (3.40 and 2.56 on two chains,
  unchanged; dec-B417). It draws nothing ([`Chain::setOffset`](../../src/bartcore/chain.hpp) hands the
  response model no generator).
- (ran) The record. The R method writes the range into the model (retired: [`recordAnchor`](../../R/dbarts.R), gone since state-frame-prior slice 1),
  so a copy and a reload hold it; `setState` with a state stored before the call returns TRUE and leaves
  the re-derived range in force (on the second fixture 33.886 before, 3.702 after the call and after the
  install).
- (ran) The test offset. With a test set and `offset = 0.7`, the R method `setOffset(data@offset,
  updateScale = TRUE)` sets `data@offset.test` to NULL and the test draws move by 0.7; the bridge's
  offset entry alone leaves both, with identical training draws and sigma.
- (ran) The keep. One row in: the call completes, the range in force stays, one warning of class
  `dbartsScaleFallbackWarning` ([`warnScaleKept`](../../R/bartcore.R)); rows back in, the next
  re-derivation reads them and warns nothing. From a range that is not every row's, shift 40.96 and
  scale 42.86 under a mask, one row in keeps both where every row gives 31.04 and 62.70. From creation's
  every-row range the keep and a read of every row give one number, so that start tells nothing.
- (critic, ran in part) The same range found again is not a no-op in bits.
  [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp) takes sigma and its prior through
  `* range_` then `/ range_`, which can move sigma by one ulp. After
  `setOffset(<offset in force>, updateScale = TRUE)` with no other change, 15 of 200 one-chain gaussian
  samplers then draw differently from an untouched twin (`setResponse` on the same 15 seeds); with two
  chains and an offset 7 of 25 gaussian, 11 of 25 Student-t and 4 of 25 survival; at most 2.5e-15 in a
  draw over 1000 sweeps (ran here: 2 of 40 after five sweeps, 0 of 40 fresh samplers). `updateScale =
  FALSE` with the same offset is inert, 200 of 200, and probit, logistic, ordinal and count samplers are
  inert at TRUE, 0 of 25 each (ran: probit, also after a mask that brings rows back).

Shapes that refuse.
- (ran) Several forests: a gaussian and a probit `forests =` sampler take `setWeights(w)` and
  `setActiveRows(a)` and refuse `setOffset` and `setResponse` with `updateScale = TRUE`
  (["does not support a sampler that carries forest amplitudes"](../../R/bartcore.R)).
- (ran) A multinomial sampler takes `setActiveRows(a)` and refuses `setOffset`, `setResponse` and
  `setWeights` outright ([`refuseCountsMutation`](../../R/bartcore.R)), and weights at creation, so only
  a later `setActiveRows` can bring it the argument. The matrix's multinomial unit-scale cell is R by
  that refusal and its runner expects the text
  (["is not available on a multinomial sampler"](../../benchmarks/R/composition-matrix.R)).
- (read) No flat entry sets weights or a mask: the two that re-anchor are
  [`dbarts_sampler_setResponse`](../../inst/include/dbarts/dbarts.h) and `dbarts_sampler_setOffset`. A
  hazard fit's sampler is a probit one on its person-period rows, and no call re-derives its periods.

A sampler made by hand (dec-B363), emulated with today's calls.
- (ran) A count list, 120 rows, half of them structural zeros left out by 0/1 weights, an offset: by
  hand the center is 0.2274 at creation and after the mask, and 0.9206 after the re-derivation, which
  is `dbarts()`'s and the log rate of the rows in. Then `getLeafPrior()`, `$model`, `$data`,
  `$activeRows`, the stored state and 20 sweeps of 2 chains are `identical()` to `dbarts()`'s.
- (ran) Where the mask keeps three rows of one count: `dbartsSpec` raises the creation warning, the
  re-derivation keeps and raises the keep's, and the sampler is `dbarts()`'s in every reader and draw
  (center 0.3542); a later keep is silent on it and warns once on `dbarts()`'s.
- (ran) A probit list with `bases =` and 0/1 weights has `active` not NULL and carries amplitudes:
  `setActiveRows(spec$active)` makes it `dbarts()`'s, draws identical, and `updateScale = TRUE` on its
  response setters is refused. A count fit refuses `forest()` terms.
- (ran) The loop, the mask redrawn each of 8 sweeps and re-derived in the first 5: the center is the
  log rate of the rows in at each of the 5 (0.483, 0.452, 0.442, -0.066, 0.720) and 0.720 for the last 3.
- (read) [`dbarts`](../../R/dbarts.R) and [`bartcoreSamplerSetData`](../../R/bartcore.R) follow a count
  fit's mask with retired: [`deriveCenterUnderMask`](../../R/bartcore.R), gone with step 3: the
  bridge's offset entry, called with the offset in force and `updateScale` true, then the record; no
  slot of the data object is touched.

## The rule

Slice B's rules 2 and 3 become, for the two setters:

1. The argument. `setWeights(weights, updateState = NULL, updateScale = FALSE)` and
   `setActiveRows(active, updateState = NULL, updateScale = FALSE)`: a single TRUE or FALSE
   ([`checkUpdateScale`](../../R/bartcore.R)), after `updateState`, so 0.9-34's positional
   `setWeights(w, TRUE)` still names `updateState` (dec-B327, dec-B381). FALSE holds the scale, as today.
2. What TRUE reads. After the weights or the mask are installed, the transform is derived again from the
   response and the offset in force over the rows then in, as creation takes it (slice B's rule 1):
   the range of the response less its offset for a gaussian or Student-t sampler, under any leaf model
   or a variance forest; the range of the observed log times less the offset for a survival sampler
   under a mask; the center for a count sampler. With every row in it is every row's, a constant
   response taking dec-B386's window. It is derived whether or not the call moved a row in or out: the
   caller asked. Nothing for one-forest probit, logistic, ordinal or hazard samplers, which take the
   argument and ignore it. sigest and the residual prior it placed are not touched (dec-B314); cut
   points, a hazard fit's periods, leaf covariate standardization and a test offset are not either.
3. The live chain. dec-B327's rule: the two calls "derive the response's scale from the rows in and
   treat the live chain and the kept draws as setResponse(updateScale = TRUE) does, the live leaves
   keeping their internal numbers and the kept draws converted". So the leaves, the kept draws, k, sigma
   and a variance forest are as Context's second list has them. The latents: the response has not
   changed, so nothing that depends on the response is redrawn; the setter's own effect on latents,
   which is to redraw those of the rows it brings in and no other (Context, ran), is what it is today
   and comes first, at the transform in force. A call on a probit sampler is the call without the
   argument. This is an orchestrator's reading and not the maintainer's words
   ([Calls made in planning](#calls-made-in-planning)).
4. Refusals of TRUE. Several forests, whatever the family, with the response setters' text (dec-B362:
   "A. Keep the refusal for 1.0-0."; its rule: "in 1.0-0 a sampler of several forests refuses
   updateScale = TRUE, the help saying the response and a re-derived basis scale may span different
   rows"); so the unit of a forest's sd is never re-derived. A multinomial sampler, through
   [`refuseCountsMutation`](../../R/bartcore.R), with the text a multinomial user already sees. A sampler
   with gp leaves that holds saved draws (dec-B330).
5. The keep. Where a row is out and the rows in hold fewer than two distinct values after the change
   (of the response less its offset, of the count), the weights or the mask are installed and the
   transform in force is kept, whatever rows it was read from: the call is the `updateScale = FALSE`
   call, bit for bit (dec-B364: "B. Keep the scale in force."; its rule names all four setters, which
   "complete, keep the prior in force and warn once per sampler"). The engine returns `scaleKept`
   ([`ResponseUpdateResult`](../../src/bartcore/sampler.hpp)) and the R method raises slice B's keep
   warning, same class and text, last. One flag for the four setters, the holder's
   `scaleKeptWarned`: dec-B364's "once per sampler" and dec-B405's "the R interface warns once per
   sampler" count samplers, not setters, so a keep through `setOffset` spends the warning a later keep
   through `setWeights` would raise. A copy or a reload is a new sampler; `setData` is not.
6. One call. The order is: every check, R's and the bridge's, the refusals of rule 4 among them; the
   install; the read of the rows in; the re-derivation or the keep; the record; the store; the warning.
   A call refused for its values (a mask that is not 0 and 1, a logistic count of 0, a survival
   sampler's weights, an `updateScale` that is not TRUE or FALSE) or by rule 4 installs nothing and
   derives nothing: the weights, the mask, their mirrors on the R object, the transform and the stored
   state are as they were, so a refused call loses its weights or mask change, as a refused response
   setter loses its response. `updateState` is the response setters': after the re-derivation the
   state is stored where it resolves to TRUE, so the stored state holds the new range; where it is
   FALSE the stored state is left, and a later `setState` of it does not put the old range back.
7. By hand (dec-B363: "We can do A', as long as the help is correctly updated when the feature is
   implemented."). dbartsSpec's help gives, for a count list whose `active` is not NULL,

       sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
       sampler$setActiveRows(spec$active, updateScale = TRUE)

   as the call that makes the sampler `dbarts()`'s, and the loop form, each pass calling
   `sampler$setActiveRows(mask, updateScale = i <= n.burn)`. A probit or ordinal list installs
   `setActiveRows(spec$active)` with no argument, having nothing to re-derive, and a several-forest one
   would be refused the argument; a gaussian or Student-t list has `active` NULL, its weights riding
   the data object into creation, and needs no call. "Equal" is `identical()`: `getLeafPrior()`,
   `$model`, `$data`, `$activeRows`, the stored state after `storeState()`, and every array `run()`
   returns from the same seed. For the count list it holds because `dbarts()` runs that call (step 3).
   Where the mask keeps fewer than two distinct counts the two are still equal in each of those; the
   by-hand route has then warned twice, `dbartsSpec` at creation and the call's keep, and has spent its
   once (rule 5), so a later keep is silent on it and audible on `dbarts()`'s. The help says so.

## Constraints

- dbarts.h, [C_interface.cpp](../../src/C_interface.cpp), the token and the C API test are not touched:
  no flat entry sets weights or a mask, and the two that re-anchor are slice B's. No ABI event.
- One new virtual on [`SamplerBase`](../../src/bartcore/facade.hpp), the re-derivation; no existing
  virtual changes its signature. `R CMD INSTALL --preclean` is mandatory on that commit.
- The re-derivation is [`Chain::setOffset`](../../src/bartcore/chain.hpp) at `updateScale` true with
  the offset in force: no response model, no chain setter and no tail gains a line or a branch. The
  keep and the family test are decided once, in the sampler, from chain 0, by `lacksDistinctRowsIn` as
  slice B's two setters decide it; on a keep, and on a family that carries no units, no chain is called.
- The new path touches no offset and no test offset, in the engine or on the data object.
- The two bridge entries gain one argument and a logical return; every other bridge entry keeps its
  arity. With the argument FALSE each makes the engine calls it makes today.
- No state format change, no new record, no new field on the reference class; the mask still does not
  ride the state. sigest, [`scaleRowsFallBack`](../../R/spec.R)'s creation rule, `xbart`, `rbart_vi`,
  the hazard grid and stan4bart are not touched.
- No test compares a call with TRUE against the call with FALSE bit for bit where a range is derived
  and found unchanged (gaussian, Student-t, survival): that pair differs in sigma's last bit at the
  rates Context gives. The pairs compared bit for bit are the identity's, the keep's, and the families
  with nothing to read.

## Open calls

None. The one this plan raised, a multinomial sampler given the argument, is settled as a refusal
(rule 4; [Calls made in planning](#calls-made-in-planning)).

## Steps

1. Engine. [`Sampler`](../../src/bartcore/sampler.hpp) gains `rederiveScale()`, returning
   `ResponseUpdateResult`: `accepted`, calling no chain, where the family carries no units
   ([`Chain::carriesUnits`](../../src/bartcore/chain.hpp)); `scaleKept`, calling no chain, where chain 0
   answers `lacksDistinctRowsIn` for the offset in force, as an incoming offset
   ([`IncomingVector`](../../src/bartcore/model.hpp)); otherwise every chain's
   `setOffset(<offset in force>, true, kept)` with the saved draw slots, then `recordAnchor`, and
   `accepted`. `Chain` gains a one-line read of the offset in force;
   [`SamplerBase`](../../src/bartcore/facade.hpp) gains the virtual and `SamplerFacade` forwards it.
2. Bridge. [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) and
   [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp) take `updateScale`. With it true each
   first refuses a sampler of several forests, a multinomial one among them
   ([`refuseMultiForestResponseMutation`](../../src/R_interface_bartcore.cpp), whose second arm today
   exempts the weights conduit), and a gp sampler holding saved draws
   ([`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp)), with the caller named; then
   validates and installs as today; then calls the re-derivation and returns
   [`firstScaleKept`](../../src/R_interface_bartcore.cpp)'s answer. With it false each returns FALSE
   after the install. Their registration and prototypes take the third argument.
3. R. The two methods of [dbarts.R](../../R/dbarts.R) take the argument (rule 1) and check it; with it
   TRUE they raise [`refuseAmplitudeMutation`](../../R/bartcore.R) and, on `setActiveRows`,
   [`refuseCountsMutation`](../../R/bartcore.R) named `$setActiveRows(updateScale = TRUE)`, before
   anything is installed; they pass it to the bridge, write the record when it is TRUE
   (retired: [`recordAnchor`](../../R/dbarts.R), gone since state-frame-prior slice 1), store by `updateState`, and raise
   [`warnScaleKept`](../../R/bartcore.R) last. `setWeights` forwards it with a 0/1 vector on a mask
   family and keeps its rollback of `data@weights` on a refused call. `reapplyActiveRows` passes FALSE:
   a re-creation restates the record and derives nothing. One helper in [bartcore.R](../../R/bartcore.R)
   holds the mask's install for the method, creation and `setData`;
   retired: [`deriveCenterUnderMask`](../../R/bartcore.R) goes, creation and `setData` passing `updateScale`
   true exactly where it ran (a count sampler, a mask installed, no amplitudes, not slice B's rule 4)
   and keeping its internal error where the bridge answers a keep. They do not store and do not warn,
   as today.
4. tests/cpp. `testRederiveScale` beside [`testScaleKeptStatus`](../../tests/cpp/test_sampler.cpp), two
   chains, its variants (plain, a variance forest, kept draws, a named sd, Student-t) and a survival and
   a count sampler:
   - after `setWeights` or `setActiveRows`, `rederiveScale()` leaves the sampler bit for bit where
     `setOffset(<offset in force>, true)` leaves a twin: the transform, the working response,
     `getAnchor`, the state and 20 sweeps of both chains;
   - the same twins where the install brings rows in, on Student-t by weights and on the survival and
     the count sampler by a mask: the latents of the entering rows moved at the install and no latent
     moved at the re-derivation;
   - the transform equals the literal over the rows in under the new weights, on a sampler whose
     previous weights left other rows out; under weights and a mask together, gaussian and Student-t,
     the literal over the rows both keep;
   - the keep, each arm starting from a range in force that is not every row's (a sampler created with
     rows out, or re-derived first): none in, one in, four in at one value give `scaleKept` and the
     sampler bit for bit a twin given the install alone, generators included; rows back in give
     `accepted` and the rows-in transform;
   - every row in: `accepted`, every row's transform; a constant response with every row in: the window
     and `accepted`; a probit sampler: `accepted` and bit for bit the install alone.
   [`SpySampler`](../../tests/cpp/test_facade.cpp) and the list of
   [`FacadeVirtual`](../../tests/cpp/test_facade.cpp) take the new virtual, which forwards the status.
5. tinytest, a new file test-update-scale-rows.R, fixtures as slice B's (rows out extreme, no missing
   predictor value) and, but for one arm, no test set. "Fails today" on every arm that passes the
   argument: it is unused.
   - The identity. For gaussian with a drawn and a fixed sigma, Student-t, a variance forest, linear
     leaves, gp leaves with no saved draw, kept draws, a named sd, a drawn k, a count sampler by
     `setActiveRows` and by 0/1 `setWeights`, and a survival sampler under a mask: the one call against
     a twin given the plain setter and then `setOffset(<offset in force>, updateScale = TRUE)`, which
     with no test set is the bridge's entry and the record: `identical()` in `getLeafPrior()`, `$model`,
     the stored state and 5 sweeps of 2 chains; with kept draws `predict` is what it was to 1e-10. Each
     family that holds latents once more with rows entering at the call.
   - The test offset. A test set and a scalar offset: the new call leaves `data@offset.test` and the
     test fits of the next sweeps those of a twin whose range was put there by the bridge's offset
     entry; the R method `setOffset` on a third sampler moves both (the pin that the two differ).
   - Literals. After each call the range, `k.scale` and the center equal values written from the rows
     in, and the readers of a sampler created under those rows; back to every row, every row's; weights
     and a mask together, on gaussian and Student-t, the rows both keep.
   - The argument. Without it, and with FALSE, the readers do not move; positional `setWeights(w, TRUE)`
     and `setActiveRows(a, TRUE)` store the state and re-derive nothing; `updateScale = NA` and `= 1`
     stop with nothing installed.
   - The keep, through each new call, on a gaussian sampler by weights and by a mask, a Student-t, a
     count and a survival sampler, each from a range in force that is not every row's, with none in,
     one in and several in at one value: the weights or the mask are installed (the mirrors, and NaN in
     the pointwise log-likelihood at the rows out); the readers and 5 sweeps equal an
     `updateScale = FALSE` twin's, and differ from every row's; exactly one warning, its class and text,
     counted with `withCallingHandlers`; none on a second keep, and none across setters in either order
     (`setOffset` then `setWeights`, `setActiveRows` then `setResponse`); a copy and a reload warn on
     their own first keep; rows back in re-derive and warn nothing; under `options(warn = 2)` the
     mirrors, the record and the stored state agree with the engine.
   - Refused, and nothing installed. Several forests on a gaussian and a probit sampler, a multinomial
     sampler, a gp sampler holding saved draws, a mask holding 0.5, a logistic count of 0, a survival
     sampler's weights: the message (the amplitude text, the multinomial text, the gp text); then
     `data@weights`, `$activeRows`, the readers and 3 sweeps equal an untouched twin's. Once more by a
     direct call of each bridge entry on a several-forest sampler, past R's refusal: the bridge's own,
     and nothing installed.
   - Nothing to read. One-forest probit, ordinal, logistic and hazard samplers: the call with TRUE is
     the call with FALSE, bit for bit, under a mask that brings rows back.
   - By hand (rule 7). A count list with an offset, 2 chains: the documented call gives `dbarts()`'s
     sampler by every comparison of rule 7. A probit and an ordinal list, and a probit `bases =` list:
     `setActiveRows(spec$active)` gives `dbarts()`'s, and on the `bases =` one the argument is refused
     with the mask not installed. A gaussian and a Student-t list with zero weights are `dbarts()`'s
     with no call. The count mask keeping three rows of one count: equal again, two warnings by hand,
     the creation text then the keep's, and a later keep silent by hand and warned on `dbarts()`'s. The
     loop form on a count sampler by `setActiveRows` and on a gaussian one by `setWeights`: in burn-in
     the center or the range is the literal over the rows in at each pass, afterwards the last burn-in
     pass's; and the loop equals a twin's written with the identity's calls.
   - The record. After a re-derivation a copy and a reload report its readers and the copy's 5 sweeps
     agree with the original's to 1e-10; `setState` of a state stored before it leaves the re-derived
     readers; `updateState = TRUE` stores a state whose range is the new one and FALSE leaves `$state`
     the object it was.
   Slice B's file is not edited: its held and creation-route arms pass as they are, through step 3.
6. The matrix. [feature-matrix.md](../design/feature-matrix.md): the sentence that `updateScale = TRUE`
   follows the unit-scale column names the four setters, and footnote f7 says a count fit's center, a
   shift, is re-derived though its cell reads "-". No cell changes and the rule needs no exception:
   the multinomial and bcf cells are R and both refuse. The runner, beside
   ["updateScale = function(s, v)"](../../benchmarks/R/composition-matrix.R): a prose probe anchored to
   that sentence runs `setActiveRows(mask, updateScale = TRUE)` on every variant and
   `setWeights(w, updateScale = TRUE)` where the case-weights cell is S, expecting what the unit-scale
   cell says: on multinomial the refusal text its `multinom|updateScale` pin already holds, reached now
   through `setActiveRows`; on bcf the amplitude text; a run elsewhere. Its pinned count of prose checks
   moves with it.
7. Help and NEWS. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the two usage lines; the
   `updateScale` item names the four setters, what the two new calls read (also when no row enters or
   leaves, and never a test offset), the order of rule 6, the keep (the weights or the mask still
   change), the refusals with dec-B362's sentence on the rows a basis scale and the response may span,
   and that a stored state does not undo a re-derivation; the weights and `active` items say "held
   unless the call carries `updateScale = TRUE`"; the multinomial paragraph's list of refused calls; the
   methods' own strings. The four enumerations of the calls that re-anchor the transform take the two
   new ones: in [dbarts.Rd](../../man/dbarts.Rd) the sentence that the sampler's `setResponse` and
   `setOffset` "only re-anchor it" and the leaf prior's list of channels a named `sd` is restated
   after; in dbartsSampler-class.Rd the state paragraph ("re-anchors it, and recorded on the") and the
   `getLeafPrior` paragraph ("whenever the response transform is re-anchored").
   [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the `active` value: rule 7's call for a count list and its
   loop, the plain install for a probit or ordinal list, and that the call under a mask keeping fewer
   than two distinct counts warns and spends the sampler's one warning, so a later keep is silent there
   and warned on `dbarts()`'s sampler; this replaces the sentence that a count sampler made by hand
   keeps every row's center. The "Whose prior is in force" paragraph of
   [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd) names the argument for a caller whose rows
   change. NEWS, new features, one item: `setWeights` takes `updateScale`.
8. Records. The design note [response-scale-rows.md](../design/response-scale-rows.md): rules 2 and 3,
   the keep's section, creation under a mask (the install with the argument) and the slices.
   [active-rows-mask.md](../design/active-rows-mask.md): a dated line where it says a mask leaves the
   transform as creation set it. TODO: `response-scale-rows` is cut to slice D;
   `rederive-same-range-ulp`, added with this plan, stays. The ledger entry for the calls made, the
   indexes, and the landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, plain and under
  ASan and UBSan; the full tinytest suite; the new file under ASan on the R-loaded path; tests/cpp under
  the sanitizers and tinytest once on x86-64 Linux, as the recent engine landings ran them.
- On a reference build, counted per scenario with no `max |z|` line: the main corpus against
  `equivalence-308c45be.rds` 55 identical; the two-forest corpus against
  `bcf-equivalence-1b7d730c.rds` 15; the multinomial corpus against
  `multinomial-equivalence-80b1c8d4.rds` 11; the four `test-reproducibility-*.R` files unchanged, 27
  results. Any other count is a stop.
- Count creation under a mask against the base build, which no corpus, snapshot or bitwise gate
  covers: [rows-pair.R](../../benchmarks/R/rows-pair.R) `record` on a library of d5b637ee and on the
  slice's, same seeds, and the two files `identical()` as R objects, every scenario; `nbinommask`,
  `countsmasked` and `countmover` are the ones that go through step 3. Not its `compare` mode, which
  holds slice B's two builds to slice B's groups.
- [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`: 42 of 42 as recorded. No arm
  is added: the new calls add no numerics, being two existing engine operations in one call, and the
  third arm of negbin-exact.R already fits through creation under a mask, which step 3 routes through
  the new entry, against an oracle whose center is over the rows in.
- [composition-matrix.R](../../benchmarks/R/composition-matrix.R) exits 0 with its new probes.
- Mutations, each to fail the named check:
  - the argument dropped in the bridge (never re-derives): step 5 "The identity", "Literals", "By hand";
  - the re-derivation before the install: step 4's literal under new weights and its entering-row
    twins; step 5 "Literals" and the identity with rows entering;
  - the keep re-deriving, over the rows in or over every row: steps 4 and 5's keep twins, which start
    from a range that is neither; the keep skipping the install: step 5's mirrors and NaN
    log-likelihood;
  - a second warned-once flag for the new setters: step 5's count across setters; the warning raised
    before the store: step 5's `options(warn = 2)` arm;
  - a refusal after the install, in R for several forests and a multinomial sampler and in the bridge
    for gp: step 5 "Refused"; in the bridge for several forests: its direct-call arm;
  - the multinomial refusal dropped in R: step 5's text (the bridge's wording differs); the runner;
  - the new path going through the R `setOffset` helper: step 5 "The test offset";
  - the saved draw slots not handed on: step 4's kept-draws twin; step 5's `predict`;
  - the latents redrawn or cold-started at the re-derivation: the entering-row twins of steps 4 and 5;
  - the engine's anchor not recorded: step 4's `getAnchor`; R's record not written: step 5's reload;
  - the argument placed before `updateState`, on either method: step 5 "The argument"; an `NA` let
    through: the same; `updateState` ignored: "The record";
  - `setWeights` on a mask family not forwarding it: step 5's count sampler by 0/1 `setWeights`;
  - creation's `updateScale` dropped in the helper: slice B's creation-route arm; the negbin arm; the
    rows-pair files.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; `R CMD check --as-cran` from a clean tarball.
- Cost: one logical read per call, and with the argument `setOffset(updateScale = TRUE)`'s passes. No
  sweep path; no bench compare.

## Budget and stops

By the header's parts: package code ~165, tests ~710, the rest ~245. On slice B's assumption (tests to
about 2.5 times their figure, the rest to 1.3) that is ~2300 lines, inside a forecast of 2000 to 2600.

Stop and report when:
- the diff passes ~2600 lines, or the package code passes ~240;
- a recorded scenario, a snapshot, an exact gate or a rows-pair file moves;
- a new call is not the identity's twin, or a keep not the `updateScale = FALSE` call, bit for bit;
- the re-derivation needs a line in a response model, a chain setter or a tail, a second facade
  virtual, a changed signature of an existing virtual, or anything in dbarts.h or the flat entries;
- a test of slice B's file fails, or `dbarts()`'s count sampler under a mask changes a reader or a draw.

## What a consumer sees

Read in each package's source on 2026-10-11; no consumer suite was run.

- stan4bart, on its branch whose range install lifts the weights (slice B's item 9): the install sets
  the weights to 1 and puts them back by `setWeights(..., updateState = FALSE)`, and a binary fit's mask
  goes on by `setActiveRows(active, updateState = FALSE)`. Each names `updateState` and passes no
  `updateScale`, so each is the call it is today and the installed range is held across them. It could
  not use the new argument for the install: its range is a supplied pair, not the rows in.
- treatSens (dbarts-1.0) makes the flat `setResponse`, `setOffset` and `setSigma` calls and none that
  sets weights or a mask. bartCause and bairrtt call neither setter. `rbart_vi` calls `$setOffset` alone.

## Slice D

Not planned and not unblocked by this slice. dec-B304 takes a linear or gp leaf's covariate center and
scale, and the default lengthscale, over the rows in at creation, held, re-derived by
`updateStandardization`. It waits on two things that are not built: the records of TODO
`state-frame-prior` (today a state install replaces the recipient's center and scale and a re-creation
derives them afresh, so nothing holds a number derived under other rows), not yet planned and without a
design in docs/design; and `updateStandardization` itself (dec-B233). One question stays open in it:
creation under a mask must derive the standardization under the mask too, and the leaf is built from
the column store, which holds no weights. This slice's argument moves the response's transform only.

## Out of scope, and where it goes

- `updateScale = TRUE` on a sampler of several forests: TODO `multiforest-scale-rederive`, after 1.0-0.
  A setter for the residual prior, how a caller moves what sigest set: TODO `sigma-prior-setter`.
- The same range found again made a no-op in bits: TODO `rederive-same-range-ulp`. It is a property of
  the two existing response setters, which this slice inherits and does not change; a fix shifts their
  draws. After 1.0-0 unless a slice touching that function takes it.

## Calls made in planning

Agent-made, each for the maintainer's later mark; the first two are the orchestrator's settlements of
the critique.

- A multinomial sampler refuses the argument on `setActiveRows`, through `refuseCountsMutation` with
  the text its other refused calls give. The feature matrix's multinomial unit-scale cell is R and its
  runner already expects that refusal; dec-B362's rule sentence is "a sampler of several forests refuses
  updateScale = TRUE"; and the package's refusals go by shape, as the gp refusal takes a probit sampler
  too. Not taken: accepting it as a call with nothing to re-derive, this plan's first recommendation,
  which made the matrix's rule need an exception and rested on a reader (`getLeafPrior()` reporting no
  `response.scale`) that a refusing probit `forests =` sampler shares.
- The latents follow `setOffset` and not `setResponse`: the response has not changed, so nothing that
  depends on it is redrawn, and the setter's own effect on latents is unchanged. The maintainer's words
  in dec-B327, "we should treat it the same as `setResponse`", were said of the scale and the leaves;
  read to the letter they would make the call redraw every latent in the likelihood (ran: 300 of 300
  on Student-t, 250 of 250 rows in on probit and count) and a call on a probit sampler a redraw where
  it is otherwise nothing. The install comes first, its latents drawn at the transform in force; the
  other order would draw them against a fit read in a transform it was not drawn in.
- The identity is stated against the bridge's offset entry, not the R method, so the new path leaves a
  test offset alone, as creation's follow-up does on purpose.
- dec-B362 and dec-B330 are read as covering the two new arguments: each rules on `updateScale = TRUE`
  for a kind of sampler, and the alternative to refusing is the work each declined. Every refusal comes
  before the install, so a refused call changes nothing.
- `setActiveRows` takes the argument after `updateState`, as `setWeights` must (dec-B327's record names
  both; only `setWeights` has a released positional form). Not taken: second, as on the response
  setters, which would make the two setters of rows differ from each other.
- One new facade virtual, called by the bridge after the install inside one entry. Not taken: the
  argument on the two existing virtuals, which changes two signatures.
- retired: [`deriveCenterUnderMask`](../../R/bartcore.R) goes and `dbarts()` runs the documented call's own path,
  so rule 7 holds by construction and not by a test alone.
- The help's by-hand call is for a count list, dec-B363's case; a probit or ordinal list keeps the plain
  install, which a several-forest one needs. The by-hand call in slice B's rule 4 case warns a second
  time and spends the sampler's once, as dec-B364 words a keep. The keep's warning keeps slice B's text
  for all four setters.
- The argument re-derives even when no row enters or leaves, so it is never the call without it in
  bits where a range is read; no fix of the last bit is planned here.
- No exact arm and no re-record: the class is neutral and the new calls are an identity with existing
  ones. A new tinytest file, not arms in slice B's, which is 1466 lines.

## Not verified

Nothing is built: no mutant and no new test was run, and the identity was run as today's calls. The
critique's rates for the same range found again, and its count of latents under `setResponse` without a
mask, were not rerun; its test-offset, weights-with-mask, keep-from-another-range and probit `bases =`
runs were. The 0.9-34 reading is from prior draws, not from a reader of the range. No reference build,
sanitizer, x86 run, `R CMD check` or consumer suite was run; that no recorded scenario moves is read
from the manifest, and rows-pair.R was read, not recorded. stan4bart and treatSens were read, not run.

## Landing note

Landed 2026-10-11 as 46f2e2a2 to 6f5b9610, in two halves: an opus implementer for the engine's
`Sampler::rederiveScale`, the facade virtual and the bridge, and a sonnet implementer for the R methods, the
test file, the matrix runner and the help. Neutral: the three bitwise compares, the snapshot files, the 42
exact gates and the recording of count creation under a mask against the base build did not move. One opus
review with 35 mutants found no defect in the package and three gaps, since closed: the position of
`setActiveRows`'s argument is pinned by a mask under which a re-derivation moves the readers, a cleared mask
with the argument is tested, and the help's fourth enumeration of the channels that re-anchor the transform
names the two setters. Two mutants are left alive on purpose: the bridge reading `NA` as TRUE on one forest,
which R's check keeps from the bridge, and the test for a family with no units, which is equivalent.

Departures from the steps: the sampler has no pointwise log-likelihood channel, so the keep's install is
pinned by the mirrors and by equality with the `updateScale = FALSE` twin's sweeps; the mask has its own
conduit in the bridge's refusal, so a multinomial sampler is not told to use the counts channel; the
multinomial refusal on `setActiveRows` names the call without the argument as the remedy. The calls made
building it are dec-A207.

