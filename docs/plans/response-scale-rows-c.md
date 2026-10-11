# response-scale-rows, slice C: setWeights and setActiveRows re-derive on request

Status: PLANNED 2026-10-11, against d5b637ee (bartcore with slice B, creation). One open call
([Open calls](#open-calls)). Implements, for `setWeights` and `setActiveRows`: dec-B296 (the setters
that change the rows take the response setters' argument), dec-B327 (its effect on the live chain),
dec-B364 and dec-B405 (the keep, its status and its one warning) and dec-B363 (the by-hand call, in
dbartsSpec's help). Rests on dec-B362 (several forests refuse), dec-B330 and dec-B331 (kept draws),
dec-B314 (sigest is never re-estimated), dec-B381 (the argument's position), dec-B386 (a constant
response's window) and dec-B417 (a drawn k lags). Slice B is [response-scale-rows.md](response-scale-rows.md),
whose rule, constraints and keep this slice extends without restating; the standing statement is the
design note of the same name. Slice D is not planned ([Slice D](#slice-d)).

agent: one push. Opus implementer for the engine, the facade, the bridge and tests/cpp; sonnet for the R
methods, the tinytest file, the matrix runner, the help and the records; one opus reviewer told to refute,
who runs the mutants of [Verification](#verification) and reads the diff twice: once for bits that must
not move where the argument is absent, once for the order of refusal, install and re-derivation. Tier:
changes draws ([Process by risk](README.md#process-by-risk)), the engine and a facade virtual being
touched; a blind critique of this plan comes first.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them.
- NEUTRAL, bit for bit: every call that exists today. The two setters without the argument, FALSE by
  default; `dbarts` and `setData` on a count fit under 0/1 weights, whose follow-up call becomes the
  same engine operation through the new entry.
- NEW, and defined by an identity: `setWeights(w, updateScale = TRUE)` and
  `setActiveRows(a, updateScale = TRUE)` leave the sampler where today's `setWeights(w)` or
  `setActiveRows(a)` followed by `setOffset(<the offset in force>, updateScale = TRUE)` leaves it, bit
  for bit, the keep and its warning included. No draw of either existing call changes.
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
budget: ~1000 lines planned: engine and facade ~40; bridge ~50; R ~70; tests/cpp ~180; tinytest ~430;
the matrix and its runner ~55; help and NEWS ~60; design note, indexes, amendments, TODO and ledger ~110.
Forecast 1800 to 2300; stop at ~2300, or at ~240 of package code ([Budget and stops](#budget-and-stops)).
No ABI event.

## Goal

A caller whose rows enter and leave the likelihood can ask, at the call that changes them, for the
response's scale and the count center to be taken again from the rows now in: `updateScale = TRUE` on
`setWeights` and `setActiveRows`, FALSE by default, with the response setters' effect on the live chain
and their keep. A sampler made by hand from `dbartsSpec` becomes the one `dbarts()` makes by one
documented call, and a loop re-derives in burn-in and holds afterwards by one argument.

## Context

Marks: (ran) on 2026-10-11 on a shipped build of d5b637ee; (0.9-34) on the released build; (read) in that
source. Gaussian fixture unless a line says otherwise: 200 rows, the last 100 out, their response 50 higher.

What the two setters do today.
- (read) [`Sampler::setWeights`](../../src/bartcore/sampler.hpp) and
  [`Sampler::setActiveRows`](../../src/bartcore/sampler.hpp) install in every chain and derive nothing;
  a chain that holds a latent for each row redraws that of every row the call brings into the
  likelihood, from its own generator, at the transform in force
  ([`TResponse::setWeights`](../../src/bartcore/model.hpp),
  [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp)). The bridge entries take two arguments
  and return NULL. The R methods are `setWeights(weights, updateState = NULL)` and
  `setActiveRows(active, updateState = NULL)`; on a probit, ordinal or count sampler `setWeights`
  forwards a 0/1 vector to `setActiveRows` ([`isMaskedWeightFamily`](../../R/spec.R)).
- (ran) Held: after `setWeights(w)` the range is creation's, shift 25.974 and scale 53.427.
- (0.9-34) `setWeights(weights, updateState = NA)` swapped the weights and held the scale: the sd of a
  prior draw of the fit was 14.3 at creation, 13.5 after `setWeights(w)` and 13.7 after
  `setOffset(0, updateScale = TRUE)` under w, where every row's range over 4 is 13.4 and the rows in
  give 0.95. `setWeights(w, updateScale = TRUE)` is an unused-argument error; `setActiveRows` reached no
  release.

What a re-derivation is today, through `setOffset(<offset in force>, updateScale = TRUE)` after the
change, which is the call slice C folds into the setters.
- (ran) It reads the rows in: after `setWeights(w)` the call gives shift 1.153 and scale 3.785, the
  range [-0.739, 3.046] of rows 1 to 100 and the readers of the sampler created under w.
- (ran, read) The live chain. The trees and the live leaves keep their internal numbers, so the live fit
  is read in the new transform until the next sweep (dec-B327). Kept draws of constant and linear leaves
  are rewritten and `predict` returns what it did (moved by 2.8e-14); a sampler with gp leaves that holds
  saved draws refuses (dec-B330; ran, the same sampler taking the call before a draw is saved). sigma
  holds its value in response units (moved by 0), its internal value and its prior restated; a variance
  forest is re-anchored ([`Chain::reanchorVarianceForest`](../../src/bartcore/chain.hpp)). A fixed k
  stays (2 before and after) while `k.scale` follows the range (26.71 to 1.89); on a second fixture a
  named sd of 0.7 has its k restated (24.20 to 2.64) and a drawn k stays (3.40 and 2.56 on two chains,
  unchanged; dec-B417). The re-derivation draws nothing: Student-t scales and survival times moved at 0
  rows ([`Chain::setOffset`](../../src/bartcore/chain.hpp) hands the response model no generator). That
  is `setOffset`'s effect and not `setResponse`'s, which cold-starts a Student-t chain's scales and df
  because its response changed.
- (ran) The record. The R method writes the range into the model ([`recordAnchor`](../../R/dbarts.R)),
  so a copy and a reload hold it; `setState` with a state stored before the call returns TRUE and leaves
  the re-derived range in force (on the second fixture 33.886 before, 3.702 after the call and after the
  install). Positional `setWeights(w, TRUE)` stores the state and moves nothing.
- (ran) The keep. One row in: the call completes, the range stays 53.427, one warning of class
  `dbartsScaleFallbackWarning` ([`warnScaleKept`](../../R/bartcore.R)); with rows back in the next
  re-derivation reads them and warns nothing.
- (ran) Families with nothing to read: on a probit sampler the call changes no bit of the next five
  sweeps, also after a mask that brings rows back; a logistic sampler takes it in silence.
- (ran) A re-derivation that finds the range in force is not exactly inert: of 40 gaussian samplers
  given `setOffset(NULL, updateScale = TRUE)` after five sweeps, 2 then draw differently from their
  twins, sigma having gone through the range and back; of 40 fresh samplers, none.

Shapes that refuse or have no say.
- (ran) Several forests: a gaussian and a probit `forests =` sampler take `setWeights(w)` and
  `setActiveRows(a)` and refuse `setOffset` and `setResponse` with `updateScale = TRUE`
  (["does not support a sampler that carries forest amplitudes"](../../R/bartcore.R)).
- (ran) A multinomial sampler takes `setActiveRows(a)`, refuses `setOffset` and `setWeights` for reasons
  of its own, and its `getLeafPrior()` reports no `response.scale`. A hazard fit's sampler is a probit
  one on its person-period rows; its periods are fixed by its rows and no call re-derives them.
- (read) No flat entry sets weights or a mask: the two that re-anchor are
  [`dbarts_sampler_setResponse`](../../inst/include/dbarts/dbarts.h) and `dbarts_sampler_setOffset`.

The by-hand count sampler (dec-B363), emulated with today's two calls.
- (ran) 120 rows, half of them structural zeros left out by 0/1 weights, an offset: by hand the center
  is 0.2274 at creation and after the mask, and 0.9206 after the re-derivation, which is `dbarts()`'s
  and the log rate of the rows in. Then `getLeafPrior()`, `$model`, `$data`, `$activeRows`, the stored
  state and 20 sweeps of 2 chains are `identical()` to `dbarts()`'s.
- (ran) Where the mask keeps three rows of one count: `dbartsSpec` raises the creation warning, the
  re-derivation keeps and raises the keep's, and the sampler is `dbarts()`'s in every reader and draw
  (center 0.3542); a later keep is silent on it and warns once on `dbarts()`'s.
- (ran) The loop, the mask redrawn each of 8 sweeps and re-derived in the first 5: the center is the
  log rate of the rows in at each of the 5 (0.483, 0.452, 0.442, -0.066, 0.720) and 0.720 for the last 3.
- (read) [`dbarts`](../../R/dbarts.R) and [`bartcoreSamplerSetData`](../../R/bartcore.R) follow a count
  fit's mask with [`deriveCenterUnderMask`](../../R/bartcore.R): the bridge's offset entry, called with
  the offset in force and `updateScale` true.

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
   response taking dec-B386's window. Nothing for probit, logistic, ordinal or hazard samplers, which
   take the argument and ignore it. sigest and the residual prior it placed are not touched (dec-B314);
   cut points, a hazard fit's periods and leaf covariate standardization are not either.
3. The live chain is treated as `setOffset(updateScale = TRUE)` treats it (dec-B327: "If they set it to
   true they're asking for the internal scale to be updated, so we should treat it the same as
   `setResponse`."): Context's second list, each line. The latents are the install's own, drawn before
   the transform moves; the re-derivation draws none.
4. Several forests refuse TRUE (dec-B362: "A. Keep the refusal for 1.0-0."; its rule: "in 1.0-0 a
   sampler of several forests refuses updateScale = TRUE, the help saying the response and a re-derived
   basis scale may span different rows"), whatever the family, with the response setters' text; so the
   unit of a forest's sd is never re-derived. A sampler with gp leaves that holds saved draws refuses
   TRUE (dec-B330). A multinomial sampler: [Open calls](#open-calls).
5. The keep. Where a row is out and the rows in hold fewer than two distinct values after the change
   (of the response less its offset, of the count), the weights or the mask are installed and the
   transform in force is kept: the call is the `updateScale = FALSE` call, bit for bit (dec-B364:
   "B. Keep the scale in force."; its rule names all four setters, which "complete, keep the prior in
   force and warn once per sampler"). The engine returns `scaleKept`
   ([`ResponseUpdateResult`](../../src/bartcore/sampler.hpp)) and the R method raises slice B's keep
   warning, same class and text, last. One flag for the four setters, the holder's
   `scaleKeptWarned`: dec-B364's "once per sampler" and dec-B405's "the R interface warns once per
   sampler" count samplers, not setters, so a keep through `setOffset` spends the warning a later keep
   through `setWeights` would raise. A copy or a reload is a new sampler; `setData` is not.
6. One call. The order is: every check, R's and the bridge's, the refusals of rule 4 among them; the
   install; the read of the rows in; the re-derivation or the keep; the record; the store; the warning.
   A call refused for its values (a mask that is not 0 and 1, a logistic count of 0, a survival
   sampler's weights) or by rule 4 installs nothing and derives nothing: the weights, the mask, their
   mirrors on the R object, the transform and the stored state are as they were. `updateState` is the
   response setters': after the re-derivation the state is stored where it resolves to TRUE, so the
   stored state holds the new range; where it is FALSE the stored state is left, and a later `setState`
   of it does not put the old range back (the help says so).
7. By hand (dec-B363: "We can do A', as long as the help is correctly updated when the feature is
   implemented."). dbartsSpec's help gives, for a list whose `active` is not NULL,

       sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
       sampler$setActiveRows(spec$active, updateScale = TRUE)

   as the call that makes the sampler `dbarts()`'s, and the loop form, each pass calling
   `sampler$setActiveRows(mask, updateScale = i <= n.burn)`. "Equal" is `identical()`: `getLeafPrior()`,
   `$model`, `$data`, `$activeRows`, the stored state after `storeState()`, and every array `run()`
   returns from the same seed. It holds because `dbarts()` runs that call (step 3). Where the mask keeps
   fewer than two distinct counts the two are still equal in each of those; the by-hand route has then
   warned twice, `dbartsSpec` at creation and the call's keep, and has spent its once (rule 5). A probit
   or ordinal list needs no re-derivation and the call is `dbarts()`'s too. A gaussian list has `active`
   NULL, its weights riding the data object into creation, and needs no call.

## Constraints

- dbarts.h, [C_interface.cpp](../../src/C_interface.cpp), the token and the C API test are not touched:
  no flat entry sets weights or a mask, and the two that re-anchor are slice B's. No ABI event.
- One new virtual on [`SamplerBase`](../../src/bartcore/facade.hpp), the re-derivation; no existing
  virtual changes its signature. `R CMD INSTALL --preclean` is mandatory on that commit.
- The re-derivation is [`Chain::setOffset`](../../src/bartcore/chain.hpp) at `updateScale` true with
  the offset in force: no response model, no chain setter and no tail gains a line or a branch. The
  keep and the family test are decided once, in the sampler, from chain 0, by `lacksDistinctRowsIn` as
  slice B's two setters decide it; on a keep, and on a family that carries no units, no chain is called.
- The two bridge entries gain one argument and a logical return; every other bridge entry keeps its
  arity. With the argument FALSE each makes the engine calls it makes today.
- No state format change, no new record, no new field on the reference class; the mask still does not
  ride the state. sigest, [`scaleRowsFallBack`](../../R/spec.R)'s creation rule, `xbart`, `rbart_vi`,
  the hazard grid and stan4bart are not touched.

## Open calls

One. Not settled here; written for a reader who has read nothing above.

1. Should a multinomial sampler accept `updateScale = TRUE` on `setActiveRows`, or refuse it?

   Background. Asking a sampler to update its scale means: read the response again and move the range
   the leaf prior is stated against. Some samplers have no such range. A probit, ordinal or logistic fit
   with one forest is one of them, and there the argument is accepted and does nothing (ran: a probit
   and a logistic sampler take `setOffset(NULL, updateScale = TRUE)` in silence). A sampler built with
   `forests =` is refused whatever its family (ran: a gaussian and a probit one both stop), which
   dec-B362 keeps for 1.0-0. A multinomial sampler is neither: it has one forest per category, takes
   nothing from the response, and has never met the argument, because its `setOffset` and `setWeights`
   are refused for reasons of their own. `setActiveRows` is the one setter of these it takes (ran), so
   this slice is where the argument first reaches it. The loop the help will show,
   `setActiveRows(mask, updateScale = i <= n.burn)`, is what would make the call.

   What a user sees, on a three-category fit of 150 rows with a mask keeping 100 of them:
   - A, accepted: the call returns in silence, the mask is installed, nothing else moves. The loop runs
     as it does on a probit fit.
   - B, refused: the first burn-in pass stops with an error and the mask is not installed. The caller
     writes `updateScale = FALSE` for this family.

   Cost. A: about three lines to let a multinomial sampler past the bridge's several-forests guard, and
   one test; a probit `forests =` sampler still refuses although it too has nothing to move, so the rule
   reads "refused where the response setters refuse" and not "refused above one forest". B: about eight
   lines for a message of its own, and one test; every caller of the loop form branches on the family.

   Recommended: A. The refusal dec-B362 keeps exists because each forest's size is stated against a
   response scale. A multinomial sampler has none (ran: its `getLeafPrior()` reports no
   `response.scale`), and the package's rule for a family with nothing to re-anchor is to ignore the
   argument, not to refuse it ([1. Structural signature](../design/feature-matrix.md#1-structural-signature),
   its footnote f7). The steps below are written to A; B changes step 2's guard, one test arm and one
   help sentence.

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
   first refuses several forests ([`refuseMultiForestResponseMutation`](../../src/R_interface_bartcore.cpp),
   whose second arm today exempts the weights conduit; a multinomial sampler passes under option A) and
   a gp sampler holding saved draws ([`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp)),
   with the caller named; then validates and installs as today; then calls the re-derivation and returns
   [`firstScaleKept`](../../src/R_interface_bartcore.cpp)'s answer. With it false each returns FALSE
   after the install. Their registration and prototypes take the third argument.
3. R. The two methods of [dbarts.R](../../R/dbarts.R) take the argument (rule 1), check it, raise
   [`refuseAmplitudeMutation`](../../R/bartcore.R) for several forests before anything is installed,
   pass it to the bridge, write the record when it is TRUE ([`recordAnchor`](../../R/dbarts.R)), store
   by `updateState`, and raise [`warnScaleKept`](../../R/bartcore.R) last. `setWeights` forwards it with
   a 0/1 vector on a mask family and keeps its rollback of `data@weights` on a refused call.
   `reapplyActiveRows` passes FALSE: a re-creation restates the record and derives nothing. One helper
   in [bartcore.R](../../R/bartcore.R) holds the mask's install for the method, creation and `setData`;
   [`deriveCenterUnderMask`](../../R/bartcore.R) goes, creation and `setData` passing `updateScale`
   true exactly where it ran (a count sampler, a mask installed, no amplitudes, not slice B's rule 4)
   and keeping its internal error where the bridge answers a keep. They do not store and do not warn,
   as today.
4. tests/cpp. `testRederiveScale` beside [`testScaleKeptStatus`](../../tests/cpp/test_sampler.cpp), two
   chains, its variants (plain, a variance forest, kept draws, a named sd, Student-t) and a survival and
   a count sampler:
   - after `setWeights` or `setActiveRows`, `rederiveScale()` leaves the sampler bit for bit where
     `setOffset(<offset in force>, true)` leaves a twin: the transform, the working response,
     `getAnchor`, the state and 20 sweeps of both chains;
   - the transform equals the literal over the rows in under the new weights, on a sampler whose
     previous weights left other rows out;
   - none in, one in, four in at one value: `scaleKept`, and the sampler bit for bit a twin given the
     install alone, generators included; rows back in: `accepted` and the rows-in transform;
   - every row in: `accepted`, every row's transform; a constant response with every row in: the window
     and `accepted`; a probit sampler: `accepted` and bit for bit the install alone.
   [`SpySampler`](../../tests/cpp/test_facade.cpp) and the list of
   [`FacadeVirtual`](../../tests/cpp/test_facade.cpp) take the new virtual, which forwards the status.
5. tinytest, a new file test-update-scale-rows.R, fixtures as slice B's (rows out extreme, no missing
   predictor value). "Fails today" on every arm that passes the argument: it is unused.
   - The two calls. For gaussian with a drawn and a fixed sigma, Student-t, a variance forest, linear
     leaves, gp leaves with no saved draw, kept draws, a named sd, a drawn k, a count sampler by
     `setActiveRows` and by 0/1 `setWeights`, and a survival sampler under a mask: the one call against
     a twin given the two (rng block), `identical()` in `getLeafPrior()`, `$model`, the stored state and
     5 sweeps of 2 chains; with kept draws `predict` is what it was to 1e-10.
   - Literals. After each call the range, `k.scale` and the center equal values written from the rows
     in, and the readers of a sampler created under those rows; back to every row, every row's.
   - The default. Without the argument, and with FALSE, the readers do not move; positional
     `setWeights(w, TRUE)` stores the state and re-derives nothing.
   - The keep, through each new call, on a gaussian sampler by weights and by a mask, a Student-t, a
     count and a survival sampler, with none in, one in and several in at one value: the weights or the
     mask are installed (the mirrors, and NaN in the pointwise log-likelihood at the rows out); the
     readers and 5 sweeps equal an `updateScale = FALSE` twin's; exactly one warning, its class and
     text, counted with `withCallingHandlers`; none on a second keep, and none across setters in either
     order (`setOffset` then `setWeights`, `setActiveRows` then `setResponse`); a copy and a reload warn
     on their own first keep; rows back in re-derive and warn nothing; under `options(warn = 2)` the
     mirrors, the record and the stored state agree with the engine.
   - Refused, and nothing installed. Several forests on a gaussian and a probit sampler, a gp sampler
     holding saved draws, a mask holding 0.5, a logistic count of 0, a survival sampler's weights: the
     message; then `data@weights`, `$activeRows`, the readers and 3 sweeps equal an untouched twin's.
   - Nothing to read. Probit, ordinal, logistic and hazard samplers, and a multinomial one under option
     A: the call with TRUE is the call with FALSE, bit for bit, under a mask that brings rows back.
   - By hand (rule 7). A count sampler with an offset, 2 chains: the documented call gives `dbarts()`'s
     sampler by every comparison of rule 7; the same for a probit and an ordinal list; a gaussian list
     with zero weights is `dbarts()`'s with no call. The mask keeping three rows of one count: equal
     again, with two warnings by hand, the creation text then the keep's. The loop form on a count
     sampler by `setActiveRows` and on a gaussian one by `setWeights`: in burn-in the center or the
     range is the literal over the rows in at each pass, afterwards the last burn-in pass's; and the
     loop equals a twin's written with the two calls.
   - The record. After a re-derivation a copy and a reload report its readers and the copy's 5 sweeps
     agree with the original's to 1e-10; `setState` of a state stored before it leaves the re-derived
     readers; `updateState = TRUE` stores a state whose range is the new one and FALSE leaves `$state`
     the object it was.
   Slice B's file is not edited: its held and creation-route arms pass as they are, through step 3.
6. The matrix. [feature-matrix.md](../design/feature-matrix.md): the sentence that `updateScale = TRUE`
   follows the unit-scale column names the four setters, and footnote f7 says a count fit's center, a
   shift, is re-derived though its cell reads "-". No cell changes.
   The runner, beside ["updateScale = function(s, v)"](../../benchmarks/R/composition-matrix.R): a prose
   probe anchored to that sentence runs `setActiveRows(mask, updateScale = TRUE)` on every variant and
   `setWeights(w, updateScale = TRUE)` where the case-weights cell is S, expecting what the unit-scale
   cell says (the amplitude refusal on bcf; a run elsewhere, multinomial included under option A); its
   pinned count of prose checks moves with it.
7. Help and NEWS. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the two usage lines; the
   `updateScale` item names the four setters, what the two new calls read, the order of rule 6, the
   keep (the weights or the mask still change), the two refusals with dec-B362's sentence on the rows a
   basis scale and the response may span, and that a stored state does not undo a re-derivation; the
   weights and `active` items say "held unless the call carries `updateScale = TRUE`"; the methods' own
   strings. [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the `active` value: rule 7's call and loop,
   replacing the sentence that a count sampler made by hand keeps every row's center. The "Whose prior
   is in force" paragraph of [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd) names the argument for
   a caller whose rows change. NEWS, new features, one item: `setWeights` takes `updateScale`.
8. Records. The design note [response-scale-rows.md](../design/response-scale-rows.md): rules 2 and 3,
   the keep's section, creation under a mask (the install with the argument) and the slices.
   [active-rows-mask.md](../design/active-rows-mask.md): a dated line where it says a mask leaves the
   transform as creation set it. TODO: `response-scale-rows` is cut to slice D. The ledger entry for the
   calls made, the indexes, and the landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, plain and under
  ASan and UBSan; the full tinytest suite; the new file under ASan on the R-loaded path; tests/cpp under
  the sanitizers and tinytest once on x86-64 Linux, as the recent engine landings ran them.
- On a reference build, counted per scenario with no `max |z|` line: the main corpus against
  `equivalence-308c45be.rds` 55 identical; the two-forest corpus against
  `bcf-equivalence-1b7d730c.rds` 15; the multinomial corpus against
  `multinomial-equivalence-80b1c8d4.rds` 11; the four `test-reproducibility-*.R` files unchanged, 27
  results. Any other count is a stop.
- [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`: 42 of 42 as recorded. No arm
  is added: the new calls add no numerics, being two existing engine operations in one call, and the
  third arm of negbin-exact.R already fits through creation under a mask, which step 3 routes through
  the new entry, against an oracle whose center is over the rows in.
- [composition-matrix.R](../../benchmarks/R/composition-matrix.R) exits 0 with its new probes.
- Mutations, each to fail the named check:
  - the argument dropped in the bridge (never re-derives): step 5 "The two calls", "Literals", "By hand";
  - the re-derivation before the install: step 4's literal under new weights; step 5 "Literals";
  - the keep re-deriving over every row: step 4's keep twins; step 5 "The keep"; the keep skipping the
    install: step 5's mirrors and NaN log-likelihood;
  - a second warned-once flag for the new setters: step 5's count across setters; the warning raised
    before the store: step 5's `options(warn = 2)` arm;
  - a refusal after the install, for several forests and for gp: step 5 "Refused";
  - the saved draw slots not handed on: step 4's kept-draws twin; step 5's `predict`;
  - the latents cold-started at the re-derivation: the Student-t, count and survival twins of steps 4, 5;
  - the engine's anchor not recorded: step 4's `getAnchor`; R's record not written: step 5's reload;
  - the argument placed before `updateState`: step 5 "The default"; `updateState` ignored: "The record";
  - `setWeights` on a mask family not forwarding it: step 5's count sampler by 0/1 `setWeights`;
  - creation's `updateScale` dropped in the helper: slice B's creation-route arm; the negbin arm.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; `R CMD check --as-cran` from a clean tarball.
- Cost: with the argument absent, one logical read per call; with it, the passes
  `setOffset(updateScale = TRUE)` makes. No sweep path; no bench compare.

## Budget and stops

By the header's parts: package code ~160, tests ~610, the rest ~225. On slice B's assumption (tests to
about 2.5 times their figure, the rest to 1.3) that is ~2000 lines, inside a forecast of 1800 to 2300.

Stop and report when:
- the diff passes ~2300 lines, or the package code passes ~240;
- a recorded scenario, a snapshot or an exact gate moves;
- a new call is not the two-call twin, or a keep not the `updateScale = FALSE` call, bit for bit;
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
- The composition matrix: no cell changes; one sentence and its probes are added (step 6).

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
- A setter for the residual prior, how a caller moves what sigest set: TODO `sigma-prior-setter`.
- A re-derivation that finds the range in force made exactly inert (Context, 2 of 40): it would change
  the bits of the two existing response setters, a draw-shifting change with a re-record; not taken.

## Calls made in planning

Agent-made, each for the maintainer's later mark.

- The new calls are the install followed by `setOffset`'s re-derivation with the offset in force, so
  the latents are treated as `setOffset` treats them and not cold-started as `setResponse` does for a
  Student-t chain: dec-B327 names the live leaves and the kept draws, and the call changes no response.
  The install comes first, its latents drawn at the transform in force; the other order would draw them
  against a fit read in a transform it was not drawn in.
- dec-B362 and dec-B330 are read as covering the two new arguments: each rules on `updateScale = TRUE`
  for a kind of sampler, and the alternative to refusing is the work each declined. Both refusals come
  before the install, so a refused call changes nothing.
- `setActiveRows` takes the argument after `updateState`, as `setWeights` must (dec-B327's record names
  both; only `setWeights` has a released positional form). Not taken: second, as on the response
  setters, which would make the two setters of rows differ from each other.
- One new facade virtual, called by the bridge after the install inside one entry. Not taken: the
  argument on the two existing virtuals, which changes two signatures.
- [`deriveCenterUnderMask`](../../R/bartcore.R) goes and `dbarts()` runs the documented call's own path,
  so rule 7 holds by construction and not by a test alone.
- The by-hand call in slice B's rule 4 case warns a second time and spends the sampler's once, as
  dec-B364 words a keep; the numbers are `dbarts()`'s. The help gives the call for a list whose `active`
  is not NULL. The keep's warning keeps slice B's text for all four setters.
- No exact arm, no pair script and no re-record: the class is neutral and the new calls are an identity
  with existing ones. A new tinytest file, not arms in slice B's, which is 1466 lines.

## Not verified

Nothing is built: no mutant and no new test was run, and the identity was run as today's two calls. The
0.9-34 reading is from prior draws, not from a reader of the range. No reference build, sanitizer, x86
run, `R CMD check` or consumer suite was run; that no recorded scenario moves is read from the manifest.
stan4bart and treatSens were read, not run, and so was the bridge's guard that option A lets a
multinomial sampler past.
