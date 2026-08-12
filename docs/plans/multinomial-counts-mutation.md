# multinomial-counts-mutation

agent: S1-S4 opus (each can move a posterior; S4 records one). S5 sonnet
  (R surface, messages, docs). Serialized: one implementer, each slice lands
  before the next starts.
rng: neutral for every path that does not call the new entries, EVERY slice.
  The counts channel adds no rng consumer; the offset path is reached only
  under a non-null offset, and the null path is today's statements verbatim
  (see "The offset indirection"). The trio bitwise is the gate at each tip. A
  trio deviation is a bug in the slice, NEVER a re-record - except at S4, whose
  whole content is one deliberate scenario addition.
window: after the cheap-uniformity S4 slice and the per-forest zero-weight
  exactness fix; BEFORE the pre-release dbarts.h reshape, so the C surface
  mirrors the final mutation shape. No `inst/include/dbarts/dbarts.h` edit in
  this arc (see "What the reshape must reserve").
budget: S1 ~210 engine/bridge + ~40 R + ~290 R test + ~110 tests/cpp;
  S2 ~320 engine/bridge + ~60 R + ~380 R/benchmark test + ~50 tests/cpp;
  S3 ~255 engine/bridge + ~60 R + ~260 test; S4 ~40 R + one baseline;
  S5 ~25 bridge + ~140 R + ~120 test + docs.

## Goal

The multinomial (softmax) sampler stops fixing its response at creation. A
counts mutation channel at fixed n and K, and an n x K category offset (train,
then test and predict), let a multinomial chain be a conditional inside a larger
Gibbs/MH sampler, as every other family already can - enabling discrete-time
competing-risks hazard, nominal IRT, discrete choice, and SUR-style
composition. Nothing is coerced and no reported channel silently omits the
offset: at every tip, a channel that cannot carry the offset refuses.

## Mechanism (the one fact the whole arc rests on)

The multinomial response is NOT the chain's `y`. It is the combiner's borrowed
`counts_` (an n x K category-major integer matrix, column k at `k*n`) plus the
derived `trials_` (`n_i = sum_k y_ik >= 1`), both set once from
`MultinomialSpec` and held in the bridge as `BartcoreHolder::ownedCounts` /
`ownedTrials`. `MultinomialForestCombiner::formForestResponse` (combiner.hpp)
names out and IGNORES both the `y` and `w` it is passed, so a `setResponse` that
installed a new `y` would write a vector nothing reads. Every response-side
entrance refuses cleanly today (audited entrance by entrance in the critique,
Part 0.1; no silently-wrong reachable path, so NO stop-loss slice is needed).
All combiner mutable state (`omega_`, `margins_`, `suffix_`, `prefix_`,
`lastF_`, `combined_`) is per-sweep scratch: the `f == 0` disjunct in
`drawForestGlue` rebuilds it at every sweep entry. Design artifacts are durable
at `.claude/multinomial-counts-mutation-design/` (`memo.md`, `critique.md`,
`synthesis.md`). **Read `synthesis.md` before starting** - it carries the B1
redesign this plan encodes and three extensions neither of the other two
documents has. `.claude/` is gitignored, so the directory does NOT arrive with
a `git worktree add`; the orchestrator copies it into the implementer's
worktree at spawn, and an implementer who cannot find it must ask rather than
proceed without it.

## Binding decisions inherited (do not reopen)

1. **Dedicated `setCounts` / `setCategoryOffset` entries, not a widened
   `setResponse`.** Multinomial handles are bare envs holding `$ptr`/`$x`
   (`bartcoreMultinomialSampler`, R/bartcore.R), no R5 method can hold one,
   `bartcore_setResponse` requires `Rf_isReal` of length n, and the shipped flat
   `dbarts_sampler_setResponse(sampler, const double* y)` cannot mirror a
   matrix weeks before the freeze.
2. **The channel carries its OWN capability predicate.**
   `Chain::supportsResponseMutation` is `combiner opt-in AND family ==
   gaussian`, and multinomial sets `family_ = logistic`, so flipping the
   combiner flag would not open `setResponse` anyway. ONE new field,
   `SamplerShape::supportsCountsMutation`; K comes off the existing
   `numReportedLocations` (see S1 item 1). Never route through
   `supportsResponseMutation`, whose contract (re-derive everything from the
   chain's `y` and `w`) multinomial cannot satisfy.
3. **A category offset enters the raw fits BEFORE the softmax, never after the
   blend.** The simplex precedent is binding: `docs/plans/multiforest-mutation-
   gaps.md` hole H1 measured reported "probabilities" summing to 10 because
   `storeSample` added a test offset after the K-blend. **The FLAT test-offset
   refusal stays in force forever** - and the helper on that entry is
   `refuseMultiForestMutation` (`C_interface.cpp:330`), NOT
   `refuseMultiForestTestOffset`, which is R-bridge-only
   (`R_interface_bartcore.cpp:2345`, called at :3959 and :3988). Both refuse
   the multinomial on the same predicate (`numForests >= 2`) and differ only in
   message, so the decision is unchanged; the earlier spelling of it simply
   named the wrong guard. The category offset must never reach
   `ResponseModel::offset()` (which `storeSample` adds post-blend to every
   location channel).
4. **Case weights stay refused**, on model grounds: a non-integer weight makes
   the augmentation `PG(w_i n_i, .)`, a real-shape draw with no exact sampler
   (the gap `docs/plans/negative-binomial.md` records), and an INTEGER weight is
   exactly row-wise count replication, already expressible. Improve the message
   to say the second half.
5. **`subset` stays refused** (Door 2 of `docs/design/model-space-survey.md`),
   as do `samplerOnly`, `warm.start` and `n.grow.sweeps`
   (`checkFamilyUnsupportedArgs`, R/bart.R).
6. **No row-centring at the entrance.** Only the row-centred part of the offset
   is identified, but `o'_ik = o_ik + a_i` gives `C'_ik = a_i + C_ik`, hence
   identical `psi`, identical working response and identical reported softmax in
   exact arithmetic. Silently rewriting a user's input buys nothing.
7. **No `dbarts.h` symbol now.** Multinomial has no flat creation path
   (`dbarts_sampler_create` routes to `createHolder` -> `createSampler` only),
   so `dbarts_sampler_setCounts` would be dead on arrival and would freeze a
   signature before use.

## Context (seams, first read at d3cb94b / engine 06c0254; REVALIDATED at
## f737702, 2026-08-12, after the multiforest-predictor-mutation arc landed
## SL+S0-S4 (33f6fdc..a825263, records through f737702) - corrections applied
## in place below and in the slices)

- Live bridge entries (`R_interface_bartcore.cpp` at f737702), since the
  multiforest arc moved every one of them: `setTreatment` :3295,
  `setForestWeights` :3319, `setOffset` :3782, `setResponse` :3800,
  `setSigma` :3825, `setData` :3832, `setTestOffset` :3951,
  `setTestPredictorAndOffset` :3974, `setWeights` :4037, `setControl` :4059,
  `setModel` :4089, `setPredictor` :4195, `predict` :4884 (whose offset guard
  now lives in the shared `predictFromSource` :4785, predicate unchanged:
  `offset != NULL && numLocations > 1`). Every multinomial combiner anchor is
  intact, line-shifted: `drawForestGlue` :800 (`lastF_ = f` :833),
  `formForestResponse` :854, `combinedFits` :870, `combinedTestFits` :885,
  `blendSoftmax` :995, `counts_` :1018, `trials_` :1019, the per-sweep scratch
  :1020-1025, `ForestCombiner::supportsResponseMutation` :427. The
  zero-weight arc's exact-zero snap landed in `BCFForestCombiner::
  formForestResponse` only; the multinomial combiner is untouched.
- **A mutation-time `totalFits` writer the original reading predates.**
  `Chain::revalidateTrees` / `rebuildFitsFromParameters` (chain.hpp:1572-1830,
  the multiforest arc's two-phase revalidation and per-observation session)
  write `forests_[f].totalFits` OUTSIDE any sweep and never touch `combiner_`.
  Multinomial reaches them: the equivalence harness's k3swap/k3txn/k3txncol/
  k3reject/k3perobs/k3perobspartial scenarios all run this path. The in-sweep
  statement below ("between `drawForestGlue(f-1)` and `drawForestGlue(f)` the
  only `totalFits` write is `finalizeTotalFits`") is still true WITHIN a sweep,
  but it is no longer the whole census, and S2's `raw_` must survive a rebuild
  that happened since the last sweep.

- `blendSoftmax` (combiner.hpp) is a template over a `fitsOf(k)` lambda SHARED
  by `combinedFits` (n, `totalFits`) and `combinedTestFits` (nTest,
  `totalTestFits`). It must stay shared and unchanged; an n x K slab cannot
  serve the test blend.
- `combinedFits()` is called THREE times per recorded sweep-equivalent, and
  still exactly three at f737702: chain.hpp:1125 post-forest-loop
  (pre-`afterCombine`), chain.hpp:4368 inside `storeSample` (which runs AFTER
  `afterCombine`), and chain.hpp:1443 in `growForestFromRoot`. Any per-column
  lazy refresh scheme leaves the last category stale at all three.
  `combinedTestFits` has one call site, chain.hpp:4405.
- `finalizeTotalFits` computes `total = forestY - resid + lastFit`, i.e. the sum
  of the forest's own tree fits. `forestY` will carry `- o_if`, so `totalFits`
  stays offset-FREE. That is the invariant `raw = totalFits + offset` needs.
- `drawForestGlue` sets `lastF_ = f` unconditionally, so a `lastF_ == f` assert
  in `formForestResponse` would be tautological. Add the invariant COMMENT; skip
  the assert.
- `predictFromSavedSampleMulti` and `predictFromCurrentTreesMulti`
  (chain.hpp:2282 / :2310) build their own raw slab and call
  `softmaxLocationMajor` directly - they never touch the combiner.
  `bartcore_predict` therefore is a THIRD reported channel. Both are now
  `template <typename Columns>` over a borrowed predictor view rather than a
  bare `const double* x_test` (the typed-ingestion re-sign), so S3's per-call
  offset argument threads through the templates, not through a pointer
  parameter.
- `MultinomialResponse` (model.hpp) reports `fitScale() == 1`,
  `fitShift() == 0`, `sigmaScale() == 1`, `initialSigma() == 1`, `offset() ==
  nullptr`, `workingWeights() == nullptr`; sigma is pinned; the leaf scale is
  the data-independent `pi*sqrt(3)/sqrt(2) / k`. **There is no response
  transform to pin**, which is why the create-vs-swap parity oracles are EXACT
  and unconditional rather than conditional on `resid.prior = fixed()` the way
  BCF's setWeights arm was.
- `PROT_COUNT` is a compile-time enum sized into every sampler's protection
  vector (R_interface_bartcore.cpp), so the counts are COPIED into the owned
  buffer and nothing is retained.
- Creation refuses non-integer TYPE via `Rf_isInteger`, which is stricter than
  a value test; `NA_INTEGER == INT_MIN` so NA is already caught by the
  nonnegativity check. Keep both.
- State: counts and offset are DATA, like predictors and `y`, and none of the
  three is serialized. `docs/design/multinomial.md` already records "State
  carries NO combiner wire blocks". **No wire-format change, no
  `stateFormatVersion` bump, no `restoreGlue` override** - the obligation is a
  contract comment: a restore reinstalls trees against WHATEVER counts the
  sampler currently holds, exactly as a single-forest restore does against the
  current `y`.
- **No comparator gates this channel.** Counts are not in `ChainStateData`, so
  `statesAgree` (tests/cpp/common.cpp) has nothing to widen; the gate must be
  behavioural. The one qualification: `trials_` drives the PG draw count, so a
  botched trials recompute desynchronizes `rngState` and would surface in a
  post-run state comparison - but no such test exists today.
- No R path reaches warm start, grow-from-root or prior-predictive for a
  multinomial handle (`installForests`/`growFromRoot`/`sampleTreesFromPrior` are
  `.Call`ed only from R5 methods in R/dbarts.R; `checkFamilyUnsupportedArgs`
  refuses the three bart2 arguments). Engine-side both sweep loops already fire
  `drawForestGlue` + `formForestResponse` per forest, so an eager swap at the
  entrance needs nothing further. An in-package direct `.Call` could still
  reach the bridge entries; the entrance validation is the guard.
- Design docs: `docs/design/multinomial.md` (the standing contract),
  `docs/plans/multiforest-mutation-gaps.md` (H1 and the capability-probe
  precedent), `docs/design/data-store.md` (read first - data-adjacent work).

## The offset indirection (decided; do not re-derive)

The semantics: the latent is `f_ik + o_ik`, so `C_ik = log sum_{j != k}
exp(f_ij + o_ij)`, `psi_ik = (f_ik + o_ik) - C_ik`, `omega_ik ~ PG(n_i,
psi_ik)`, forest f's working response is `(y_if - n_i/2)/omega_if + C_if -
o_if` under weight `omega_if`, and the reported channels are `softmax(f + o)`.
`afterCombine` is UNCHANGED: the offset is not a leaf value, the conditional for
the level shift `c` is derived from the leaf prior alone, and a common shift of
all K is a softmax null direction.

The implementation shape, which is what makes neutrality structural:

1. One private accessor is the single definition of "the fit the softmax sees":
   `rawFits(k, forests)` returns `forests[k].totalFits.data()` when
   `offset_ == nullptr`, and `raw_.data() + k*n` otherwise. Every train-side
   reader goes through it - the suffix fold, the prefix seed, the continuation
   fold, `fFits`, and `combinedFits`' `fitsOf`. On the null path it yields
   today's pointer, so every loop is byte-for-byte today's loop.
2. `raw_` is maintained ONLY under an offset, as
   `raw_[k*n+i] = totalFits_k[i] + offset_[k*n+i]`, and materialized in FULL at
   the `f == 0` fresh-sweep branch and at the TOP OF `combinedFits`. The
   in-order continuation may refresh only column f-1, because WITHIN A SWEEP,
   between `drawForestGlue(f-1)` and `drawForestGlue(f)`, the only `totalFits`
   write is `finalizeTotalFits(forests_[f-1])`. Every path that REPORTS
   rematerializes in full, so no reported fit is ever mixed-vintage.
   **The full-rematerialization pair is load-bearing for a second reason the
   original reading predates**: `revalidateTrees`/`rebuildFitsFromParameters`
   (chain.hpp:1572-1830) rewrite `totalFits` at PREDICTOR-MUTATION time,
   between sweeps, without entering the combiner - so a `raw_` maintained only
   by the in-sweep continuation would be stale from the mutation until the next
   sweep. The `f == 0` and `combinedFits`-top rematerializations both dominate
   that write, so the scheme is correct as specified; do NOT weaken either one
   into a lazy refresh, and do not add a "refresh at mutation" hook (the
   combiner is not on that path and must not be put on it).
3. `blendSoftmax` is UNCHANGED. `combinedTestFits` gets a symmetric
   `rawTestFits(k)` over its own nTest x K `rawTest_` slab, materialized in full
   inside `combinedTestFits` under a test offset and returning
   `totalTestFits.data()` otherwise.
4. `formForestResponse` SUBTRACTS `o_if` and so cannot route through `raw_`:
   a hoisted two-branch loop whose null branch is today's statement verbatim.

Cost: the offset path pays one extra O(nK) pass per `combinedFits`; the null
path pays zero - no inner-loop branch, no extra traffic.

## Constraints

- Refuse, never coerce. A refusal leaves the sampler BYTE-IDENTICAL: build and
  fully validate a `std::vector<int>` scratch, then `std::swap` it into
  `ownedCounts` (likewise trials), then call the engine setter with the new
  `.data()`. The combiner borrows `ownedCounts.data()`, so an in-place write
  IS the mutation and a mid-write error leaves half-new counts.
- **n and K are out of scope.** n is pinned by every combiner buffer and the
  forests' n-sized allocations; K is the forest count and cannot change on a
  live sampler. Refuse a length mismatch naming BOTH n and K.
- A one-ULP multinomial deviation on the null-offset path is a codegen
  falsifier, not a tolerance to widen - it means `rawFits` did not resolve to
  today's pointer.
- Refusal text must not promise a future. A temporary refusal says so in this
  plan, not in the message.
- Out of scope: n or K change; real-shape (non-integer) counts; case weights;
  per-forest row subsetting; a public multinomial sampler surface; the
  tests/cpp `ConfigSpec` multi-forest fuzz arm; any `dbarts.h` edit.

## S1. The counts channel

1. Engine: `ForestCombiner::setCounts` base virtual (inert) plus the
   multinomial override swapping `counts_` and `trials_`; a
   `supportsCountsMutation()` probe defaulting false on the base combiner,
   exactly as `supportsResponseMutation` does, so a future coupling stays
   refused until audited. `Chain` fan-in, `Sampler` fan-out (the
   `setTreatment` pattern), `SamplerBase` virtual and ONE new `SamplerShape`
   field. **No `numCategories` field**: `SamplerShape::numReportedLocations`
   already is K for multinomial (facade.hpp:42-44) and 1 for every additive
   model, so a second name for the same number would be a second thing to keep
   in sync. Read K off `numReportedLocations`, gated by
   `supportsCountsMutation`.
2. Bridge `bartcore_setCounts`: capability probe FIRST (never a `numForests`
   test - the `setTreatment`/`bcfGlue` precedent), then validate into scratch
   (integer type, `length == n*K` naming both, nonnegative, `sum_k y_ik >= 1`
   per row restated HERE not inherited from creation, `sum_k y_ik <= INT_MAX`),
   then swap, then call the setter with fresh pointers. Copy; retain nothing.
3. Add the invariant comment on `formForestResponse` naming the
   `drawForestGlue(f)`-immediately-precedes contract. No assert.
4. Document the PG cost cliff at the entrance: the sweep draws `n_i` PGs per
   observation per category, so a swap to large row sums multiplies sweep cost
   by `mean(n_i)`.
5. **Message repair (B3): the response conduit.** The guard is no longer an
   inline `numForests >= 2 && !supportsResponseMutation` test in
   `bartcore_setResponse`: the multiforest arc folded response, offset and
   weights into ONE conduit-parameterized helper,
   `refuseMultiForestResponseMutation(sampler, caller, ResponseConduit,
   updateScale)` (`R_interface_bartcore.cpp:2483`, declared in
   `R_interface_bartcore_common.hpp`), and the FLAT C entries call the same
   helper (`C_interface.cpp:206` / :226 / :238). So the repair goes inside the
   shared helper's `!supportsResponseMutation` branch, keyed on the NEW
   `supportsCountsMutation` probe, and the counts hint lands on BOTH surfaces.
   That is correct and deliberate: the hint is true on both, and a helper
   whose whole reason for existing is that "the two surfaces cannot state
   different rules" must not be made to state two. The generic text stays
   verbatim for a future non-opt-in coupling, and the `ResponseConduit`
   wording (`fixes its response at creation` / `carries no offset` / `fixes its
   case weights at creation`) is preserved per conduit. Re-pin in this slice:
   `inst/tinytest/test-multi-forest-seam.R` (the multinomial block, currently
   :231-282, which pins `setResponse` at both `updateScale`, `setOffset` at
   TRUE/FALSE/NA, `setWeights`, `setTreatment`, `setTestOffset`) AND the flat
   surface's expectations, `inst/tinytest/capi/consumer.c` with
   `inst/tinytest/test-capi.R` - the shared helper means an unpinned flat
   message drifts silently.

Tests: **O1** counts create-vs-swap parity - `create(A) -> setCounts(B) ->
run(burn, m)` bitwise identical on all five recorded channels to
`create(B) -> run(burn, m)` at the same seed, with the non-vacuity arm
`create(A) -> run` differing from `create(B) -> run`. **O5 BURNED-IN self-swap**
- `create(A) -> run(b, m1) -> setCounts(A) -> run(0, m2)` bitwise identical to
`create(A) -> run(b, m1) -> run(0, m2)`; the control is the SAME SPLIT without
the swap, never `run(b, m1+m2)`, so the oracle does not depend on split == single
(pin split == single separately if you like; do not build on it). **O11
same-vintage cross-check** at the null offset, as a pre-existing invariant:
for the last recorded sample of a single-chain run,
`softmax_k(bartcoreForestFits(bc, k)) == train[, k, last]` to 1e-12 (NOT
bitwise - an R-side softmax does not reproduce `softmaxLocationMajor`'s
reduction order). **O7 counts half** plus the multinomial `setData`/`setModel`/
`setSigma` refusal pins the existing battery lacks. **O8 tests/cpp** combiner
unit: build with A, `setCounts(B)`, glue + form, compare against a combiner
built with B from the same seeded rng, WITH a burned-in arm (several
glue/form cycles before the swap). **`tests/cpp/test_shape.cpp`**: the shape
POD is field-oracled there (`CHECK_SHAPE_FIELD` per member, with a multinomial
arm at :298-310), so the new `supportsCountsMutation` field is not optional
test work - extend the oracle in this slice, asserting true on the multinomial
arm and false on the single-forest and BCF arms.
Gate: trio bitwise (the a825263 baselines - see "Verification"); tests/cpp from
clean plus the ASAN/UBSAN leg (new engine code on a live path); full tinytest,
no snapshot regenerated; `air format --check .` and `lintr::lint_package()`;
rchk next scheduled run. `R CMD INSTALL --preclean` MANDATORY (facade.hpp
virtuals move - stale objects bus-error) and delete the `benchmarks/kernels`
binaries.

## S2. The n x K train offset, engine + internal creation, with the floors

1. Engine: the indirection of "The offset indirection" above - `raw_`,
   `rawFits`, the two materialization points, the offset in `drawForestGlue`
   and `formForestResponse`, `setCategoryOffset` on the combiner,
   `MultinomialSpec.offset`; chain/sampler/facade passthrough. NULL clears.
2. Bridge `bartcore_setCategoryOffset` plus the offset argument on BOTH internal
   creation entries: real matrix of exactly n x K, every entry FINITE (an `Inf`
   propagates through the log-sum-exp into a NaN margin and poisons every
   category), NULL clears. The category offset goes to the COMBINER; it must
   never reach `ResponseModel::offset()`.
3. **Mandatory refusal floors (B2), regardless of S3.** With a category offset
   installed: refuse test data, refuse installing an offset while test data is
   present, and refuse `bartcore_predict` on a multinomial carrying one. All
   three are TEMPORARY; S3 lifts them. Without these, `combinedTestFits`,
   `storeSample`'s test channel and the two predict replays each report
   `softmax(f_test)` - the H1 class exactly.
4. **Message repairs (B3): `bartcore_setOffset` and `bartcore_predict`**, both
   probe-conditioned, both re-pinned in the same slice. Note that predict's
   existing text is NARROW, not false: a flat per-row offset really is the
   softmax null direction and stays undefined; the repair points at the matrix
   form.
5. R: `bartcoreSetCategoryOffset` plus the `offset` argument on both internal
   creators. Validate R-side too (safe over fast in R).

Tests: **O2** offset create-vs-swap parity - `create(offset = O)` bitwise
identical to `create(offset = NULL) -> setCategoryOffset(O)` over the same run,
with a non-vacuity arm proving the offset moves the answer. **O3
null/zero-offset neutrality as a HARD gate**: `setCategoryOffset(NULL)` returns
to the exact null path bitwise, AND an all-zero matrix is bitwise identical to
NULL (settled: `logSumExp2` compares, `exp(-0.0) == 1`, and
`ext_rng_simulatePolyaGamma` opens with `0.5 * fabs(psi)`, so -0.0 is absorbed
at every consumer). **O4 EXACT-GATE offset arm, MANDATORY** - a new
intercept-only K = 3 arm in `benchmarks/R/multinomial-exact.R` with a fixed
nonzero offset, matching the posterior mean of the identified probabilities to
quadrature with the offset IN the likelihood. This is the only oracle that
catches a consistently-wrong PLACEMENT or SIGN (parity oracles agree happily if
the offset is added to the margin but not subtracted from the working response);
it is scoped to sign and placement only, since the common-component question is
an exact identity, not a measurement. **O11 with the offset** -
`softmax_k(forestFits[, k] + offset[, k]) == train[, k, last]` - the direct
falsifier for a mixed-vintage `raw_`. **O11-after-session**, its second arm and
the cheapest guard on the mutation-time rebuild named in the Context: run,
install a category offset, drive a per-observation predictor session
(`bartcoreUpdatePredictorPerObservation`, which the multinomial handle accepts
since the multiforest S2 widening) so `rebuildFitsFromParameters` rewrites
`totalFits` between sweeps, run one more sample, and re-check the same
identity. **O6** simplex invariant with an offset installed. **O7 offset half**
plus the three new floor refusals.
Gate: as S1, plus the FULL `benchmarks/R/multinomial-exact.R` (all arms, not
only the new one; note `.github/workflows/exact-gates.yaml` runs this harness
on every push to bartcore, so O4's new arm becomes a CI gate the moment it
lands). `--preclean`; delete the kernel binaries.

## S3. The nTest x K test offset and the predict replays

1. Engine: `rawTest_` and `rawTestFits` inside `combinedTestFits`; the
   nTest x K test offset storage and its setter; the two multi predict replays
   (`predictFromSavedSampleMulti` chain.hpp:2282, `predictFromCurrentTreesMulti`
   :2310 - both now `template <typename Columns>` over a borrowed predictor
   view, so the offset is a new function parameter threaded through the
   templates and their `bartcore_predict` / `predictFromSource` callers) take a
   per-call nNew x K offset added to their raw slab BEFORE
   `softmaxLocationMajor`.
2. Bridge: `bartcore_setCategoryTestOffset`; the creation-side test-offset lift
   for the MATRIX form only, on the R BRIDGE entries only. The FLAT entries stay
   refused, by `refuseMultiForestMutation` on `dbarts_sampler_setTestOffset`
   (`C_interface.cpp:330`) - the R-bridge-only `refuseMultiForestTestOffset`
   (`R_interface_bartcore.cpp:2345`) is the guard this slice conditions, not the
   flat one. Do not touch the flat guard. `bartcore_predict` gains an
   nNew x K matrix offset form.
   **Predict reads its ARGUMENT, never the resident test offset** - a row-count
   coincidence is not consent. If the sampler carries a category offset and no
   offset argument is supplied, REFUSE; an explicit all-zero matrix is how a
   caller asks for the offset-free surface.
3. Lift S2's three floors. Message repair: `bartcore_setTestOffset`'s text,
   probe-conditioned, re-pinned here.
4. R: `bartcoreSetCategoryTestOffset` and the predict passthrough.

Tests: test-channel create-vs-swap parity; predict-vs-run agreement (a
`keepTrees` predict on the resident test rows with the same offset matrix equals
the run's recorded test channel); the simplex guard with a test offset
installed; the flat `setTestOffset` refusal STILL firing; the
predict-without-offset-argument refusal and its explicit-zero escape; the three
lifted floors now accepting.
Gate: **as S2, plus a MANDATORY ASAN/UBSAN leg** - `rawTest_` is a new heap
buffer sized by `numTestObservations`, the exact class the variance arc's
mandatory S4 ASAN leg caught a live defect in. `--preclean` (combiner.hpp and
chain.hpp are headers); delete the kernel binaries.

## S4. Fixture scenario and baseline re-record (its own commit)

Not folded into S5: a baseline re-record IS a snapshot move, and S5's gate is
"no snapshot moves". Landing after S3 also means ONE recording covering train
and test offsets instead of two.

1. Add a `k3offset` scenario to `benchmarks/R/multinomial-equivalence.R`
   exercising a train AND a test category offset across all five recorded
   channels (`recordChannels`: train, test, per-category `forestFits`,
   per-category `varcount`, `runVarcount`). (Written for the S3-in-scope shape;
   if S3 is deferred - see "Open items" - the scenario carries the TRAIN offset
   only and its test channel is the S2 floor's refusal, not a recorded array.)
2. Re-record to `multinomial-equivalence-<hash>.rds`; mark
   `multinomial-equivalence-a825263.rds` historical in
   `benchmarks/baselines/MANIFEST` with the PARTITION statement the 2bd34db /
   8c2b5fc entries establish: inside a compare against a825263 the NINE
   existing scenarios (k3, k2, k3counts, k3swap, k3txn, k3txncol, k3reject,
   k3perobs, k3perobspartial) reproduce bitwise and only the new scenario is
   absent.
3. Re-pin the three CI legs in `.github/workflows/equivalence.yaml` - :61
   (gaussian, `equivalence-a825263.rds --strict-coverage`), :87 (bcf,
   `bcf-equivalence-a825263.rds`), :113 (multinomial) - so the multinomial line
   names the NEW hash and the other two stay where they are. The workflow does
   not fire from bartcore (its own header records this; it is manual-dispatch
   only until it lands on main), so the edit is a bookkeeping obligation the
   local run will not catch.
4. Every later gate (S5 and everything downstream) names the NEW hash,
   10 scenarios x their channels bitwise.

Gate: the partition statement above verified, not asserted; the new baseline
reproduces itself; `equivalence-a825263` (35 scenarios, strict) and
`bcf-equivalence-a825263` (11 scenarios) unmoved.

## S5. Public surface, remaining messages, docs

1. R/bart.R: LIFT the multinomial `offset` refusal for an n x K NUMERIC MATRIX
   only, threaded into `bart2Multinomial` / `bart2MultinomialCounts`. A flat
   vector stays refused, naming why (it points exactly along the softmax's null
   direction and is identically inert). Keep `parseMultinomialData`'s refusal of
   the HOST data object's flat offset - a different pointer, still inert.
2. The improved `weights` refusal message (binding decision 4).
3. Docs: `docs/design/multinomial.md` "The surface" gains the offset and the
   mutation contract, replacing the whole-data-mutation-is-refused bullet;
   the landing record in this file; `?bart2` Rd; `inst/NEWS.Rd`.
   **Record the reservation in `docs/plans/c-api-growth.md`**, not only here, or
   it is lost when the reshape is scoped (next section).

Gate: full tinytest with NO snapshot regenerated (S5 is RNG-neutral surface);
`R CMD check`; `air format --check .` and `lintr::lint_package()`; trio bitwise
against the S4 baseline (the new multinomial hash, plus `equivalence-a825263`
strict and `bcf-equivalence-a825263`).

## What the reshape must reserve (record in c-api-growth.md)

1. **Source-shaped response and offset parameters.** The reshape already
   re-signs predictor entries onto a borrowed `PredictorSource` view; do the
   same on the response side, replacing the bare `const double* y` /
   `const double* offset` with a tagged source able to express at minimum
   `{ double* vector }` and `{ int* counts, size_t numCategories }` for the
   response, and `{ double* vector }` / `{ double* matrix, size_t
   numCategories }` for the offset. A later flat multinomial creation entry then
   needs a new tag, not an ABI break.
2. **A size-first spec struct** for any future flat multinomial creation, per
   the `dbarts_results` `structSize` precedent this arc's sibling established.
3. **Whatever tagged struct the reshape adopts must PRESERVE the
   `refuseMultiForestMutation` guards** already on the flat response-side
   entries (the D4 precedent in `docs/plans/multiforest-mutation-gaps.md`). The
   reshape is exactly where they get dropped by accident.

## Falsifiers

- **F1 (every slice).** Any trio deviation is a bug in the slice, never a
  re-record. A one-ULP multinomial deviation on the null-offset path means
  `rawFits` did not resolve to today's pointer - abort, do not tune.
- **F2 (S2).** O4's offset arm must be written FAILING-FIRST against a
  deliberately sign-flipped build. An exact gate that has never been red is not
  a gate.
- **F3 (S1).** O5 burned-in must be shown to go red under two deliberate
  breaks: the trials recompute disabled, and `omega_` re-seeded at `setCounts`.
- **F4 (S1, S2).** Both parity oracles need their non-vacuity arms
  (`create(A)` differs from `create(B)`; the offset moves the answer). A parity
  test whose two arms cannot differ proves nothing.
- **F5 (S2).** O11 with an offset must be shown to go red when `combinedFits`'
  full rematerialization is reduced to the lazy per-column refresh - that is
  the B1 defect, and it is the one thing the null-path trio cannot catch.
- **F6 (S3).** Predict-vs-run agreement must be shown to go red when the offset
  is threaded into `combinedTestFits` but not into the two predict replays.

## Edge cases the tests must name

The burned-in swap has NO exact create-side twin (grown trees, `lastF_ == K-1`,
live `omega_`), so O1 is not coverage of it and O5 is the only discriminating
check - say so in the test file. A zero-row-sum row in a REPLACEMENT counts
matrix (`PG(0, .)` is a point mass at zero and the working response divides by
omega). `sum_k y_ik` overflowing `int`. A length that matches `n*K` numerically
but comes from a transposed matrix (refused on shape, naming n and K). NA in
counts (already caught by the nonnegativity test since `NA_INTEGER == INT_MIN`;
pin it anyway). A non-finite offset entry, each of NA/NaN/+Inf/-Inf. An
all-zero offset matrix (bitwise equal to NULL - hard gate). An offset with a
constant added to every entry of a row (identical reported probabilities, the
exact identity of binding decision 6). `setCategoryOffset(NULL)` after an
offset run, returning to the exact null path. `setCounts` and
`setCategoryOffset` on a gaussian / BCF / heteroscedastic sampler (refused by
capability probe, message naming the family situation, not the forest count).
`setState` after a `setCounts` (trees restore against the CURRENT counts).
Multi-chain: every chain sees the swap.

## Doors held open (recorded, not scheduled)

- **A public multinomial sampler surface**, the direct D3 `bcf-public-surface`
  analog. After this arc the channel is reachable only through
  `dbarts:::bartcore*` - sufficient for bairrtt, stan4bart and an in-repo
  competing-risks front end, NOT for an arbitrary CRAN user. `inst/NEWS.Rd`
  must not describe the channel in terms an external user reads as available.
- **The missing R wrapper for `installForests` on a multinomial handle.** The
  engine path is correct (the K check rides `src.forests.size() !=
  numForests()`), and a donor's trees against the destination's counts IS the
  right warm-start semantics; there is simply no wrapper.
- **`samplePriorPredictive` for multinomial.** If a low-level wrapper is ever
  added, `type = "ppd"` must refuse - there is no `n_i` to draw counts against
  (the S1 precedent from the variance arc).
- **`forceRefreshTrees`' unweighted collapse under multinomial.** It loops all
  of `forests_`, but `collapseEmptyNodes` is called with
  `response_->workingWeights()`, which `MultinomialResponse` returns as
  `nullptr` (unit weights), whereas the sweep merges against `omega`. The merged
  leaf value is the unweighted mean; it is redrawn on the next sweep, so the
  blast radius is the reported fits between the mutation and the next draw.
  Pre-existing, out of arc. `setCutPoints` carries the same exposure.
- **The tests/cpp `ConfigSpec` multi-forest fuzz arm.** `ConfigSpec` has never
  covered ANY multi-forest sampler, BCF included, so scope it as the general
  gap plus an `OP_SET_COUNTS` op, not as a multinomial follow-up.

## Verification (every slice)

- `R CMD INSTALL --preclean` for S1-S3 (facade.hpp virtuals at S1; combiner.hpp
  and chain.hpp at S2-S3). Delete the `benchmarks/kernels` binaries after any
  header edit - they have no header dependency tracking; `tests/cpp` tracks
  headers via `-MMD`.
- `cd tests/cpp && make clean && make && ./test_bartcore` - all pass. ASAN/UBSAN
  leg for S1, S2 and S3 (each makes new engine code reachable), MANDATORY at S3.
- Full `tinytest::test_package("dbarts")` from a preclean install. New tests
  ADD; no snapshot is regenerated at any slice. A forced snapshot is a signal
  the slice changed more than intended - stop and report.
- The trio, EVERY slice, expecting no deviation. The multiforest arc re-recorded
  all three at a825263, so these are the CURRENT names per
  `benchmarks/baselines/MANIFEST` (which is authoritative; re-read it at spawn
  rather than trusting this line):
  `benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-a825263.rds --strict-coverage`
  -> 35/35 "identical draws (same RNG stream)";
  `bcf-equivalence-a825263.rds` -> 11 scenarios x their channels bitwise;
  `multinomial-equivalence-a825263.rds` -> 9 scenarios x their channels bitwise
  (S1-S3; the new hash at 10 scenarios from S4 on). The multinomial nine are
  k3, k2, k3counts, k3swap, k3txn, k3txncol, k3reject, k3perobs,
  k3perobspartial - six of them predictor-mutation scenarios the multiforest
  arc added, so this leg now exercises `revalidateTrees` and the
  per-observation session on a multinomial chain and is a real guard on S2's
  `raw_`, not only on the untouched paths. No max-|z| line anywhere. THE TRIO IS
  NECESSARY, NOT SUFFICIENT: it proves only that nothing moved while the channel
  is unused. O1-O11 are the real oracle, and no state comparator can ever gate
  this channel.
- `air format --check .` and `lintr::lint_package()` on any slice touching R/
  (both are enforced by `.github/workflows/lint.yaml`). `R CMD check` on any
  slice touching R/ or Rd - S2 and S5 both do. rchk on the next scheduled run
  for S1-S3 and S5 (the bridge moves).
- `tests/cpp/test_shape.cpp` field-oracles `SamplerShape` member by member, so
  any slice adding a shape field extends it in the same slice (S1).
- Speed: no slice touches `run()`'s inner loop and no new entry is inside
  `bench-sampler.R`'s timed arms, so `bench-sampler-ab1dc52.csv` is a formality.
  The offset path's extra O(nK) pass is off every timed arm; confirm before
  skipping.

Stop conditions per `docs/plans/README.md`: a step fails twice, the diff exceeds
1.5x budget, or a needed change is out of scope - report and stop. Budgets here
are sized to the MANDATED ORACLE, not the engine delta, per the variance arc's
standing lesson; S2 and S3 are the test-heavy slices.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S1: a multinomial (softmax) sampler's count response can be replaced between
  runs, so it can serve as a conditional inside a larger sampler; the trees
  carry over, fitted to the previous counts.
- S2: a multinomial sampler accepts an n x K category offset, entering the
  linear predictor before the softmax. While one is installed, test data and
  out-of-sample prediction are refused rather than reported without it.
- S3: the category offset extends to test data and to out-of-sample
  prediction, which takes its own offset matrix; the S2 refusals are lifted.
- S5: `bart2(family = "multinomial")` accepts a matrix `offset`; a vector
  offset is refused by name, as it is inert under the softmax.

## Open items

- **Q7 (public reachability) is decided: `dbarts:::`-only, no public
  multinomial sampler surface scheduled alongside.** The door is recorded
  above. The ONE user-visible change in the arc is S5's `bart2` matrix
  `offset` - a creation-time argument on a one-shot fit, not a mutation
  surface. Flagged rather than buried.
- S3 is the droppable slice. Because B2's floors land in S2, the arc closes
  honestly at the S2 tip; if the queue tightens before the reshape, defer S3 -
  it is the only slice whose absence costs capability rather than correctness.
  This is the one scope question to settle with VD BEFORE S1, because it
  decides S4's content: with S3 in, the `k3offset` fixture carries a train AND
  a test category offset; with S3 out, it carries the train offset only and
  the test channel records S2's refusal. S4 is written above for S3-in.
  RESOLVED 2026-08-12 (VD): S3 stays IN SCOPE; S4 as written stands.
- Design artifacts are durable at
  `.claude/multinomial-counts-mutation-design/` (`memo.md`, `critique.md`,
  `synthesis.md`; the critique was an independent adversarial pass, read-only
  on source, re-anchored against engine tip 06c0254). The critique's verdict
  was STANDS WITH AMENDMENTS; all five blocking and all eight advisory
  amendments are ADOPTED, with NO overturns.
  Three are adopted with corrections recorded in the synthesis: B1's refresh
  redesign (a `rawFits` accessor that is today's pointer on the null path, full
  rematerialization at every reporting point, `blendSoftmax` untouched), B2's
  predict-offset-argument resolution, and B5's split-run control arm plus the
  new O11 same-vintage oracle.

## Landing notes

S1 LANDED 729bbdd, 2026-08-12. The counts channel: base-inert
ForestCombiner::setCounts + supportsCountsMutation with the
multinomial override; Chain/Sampler/facade fan-out; SamplerShape
gains the capability field with K read off numReportedLocations (no
numCategories field, per the revalidation); test_shape.cpp oracle
extended (multinomial true, single-forest and BCF false); bridge
bartcore_setCounts runs the capability probe first, whole-scratch
validation under unwindProtect, then swap-then-install - the
verifier probed 13 hostile inputs (wrong n, K mismatch, transposed,
no dim attr, 3-d array, real, logical, NA, negative, zero row,
overflow row, out-of-int-range double, valid self-swap): every one
refused and the sampler's subsequent 4-sample run stayed bitwise
identical to an untouched baseline; dbarts:::bartcoreSetCounts; the
B3 message repair in the shared conduit-parameterized
refuseMultiForestResponseMutation - the counts hint is emitted on
both surfaces by the helper, while the flat-surface re-pin is
documentation-only, verified sound: dbarts_sampler_create has no
multinomial branch, so the flat counts path is unreachable and
consumer.c's BCF legs pin the capability conditioning; NEWS bullet.
Deviations accepted at verification: the entrance checks
DIMENSIONS, not length (a transposed n x K matrix has n*K entries
and slipped a length test; pinned in the test file); O5 built over
grouped counts so a lost trials recompute shows in the self-swap;
O7 covers setData/setModel since multinomial setSigma was already
pinned; 807 changed lines vs the ~650 budget, overage all oracle
content. Falsifiers: F3a (trials recompute disabled) red, 14
failures; F3b (omega re-seed) PROVEN unfalsifiable by construction,
not a masked hole - omega_ has exactly four sites in combiner.hpp
and formForestResponse's only two engine call sites each sit
immediately after drawForestGlue for the same forest, so no
first-sweep, post-setState or post-swap path reads omega before
that sweep wrote it (verified independently by implementer and
verifier); replaced by the stronger F3b' (an rng draw inserted at
Chain::setCounts), red with 14 failures - O1's five channels plus
label-entry and two-chain arms plus O5's five, proving O5 live
(the implementer counted 12; the verifier measured 14 and derived
14 analytically; 14 is the record). The four oracle identifiers in
the test file's section comments were reworded to self-contained
rationale at landing per the comment-anchor discipline.
Message-quality carry to S2: R/bartcore.R coerces counts to integer
after the nonnegativity check with no range check, so a count above
INT_MAX errors with a misleading nonnegativity message (sampler
untouched either way) - range-check before coercion in S2's pass
through the same entry. Implementer carries to S2: O11 is pinned at
the null offset, so S2 only adds the offset term; the counts
refusal message needs its setCategoryOffset half in the same helper
branch (the seam file's setOffset pins move with it); raw_ cannot
lean on parity oracles - O11 is the oracle that can see a mixed
vintage. Gates double-run (implementer + independent verifier,
fresh privlibs): install --preclean, tests/cpp plain + ASAN/UBSAN
clean, tinytest 4047/0 (44 new) no snapshot regenerated, trio
35/35 strict / 11/11 / 9/9 vs a825263 with no re-record,
multinomial-exact all five arms at gaps 3e-4/8e-4/5e-4/5e-4/1e-4
vs tolerances 8e-3..1.5e-2, bcf oracles at the MANIFEST values,
air 0, lintr::lint_package 0, R CMD check clean-copy tarball
Status OK.

S2 (the n x K train offset) is next; S3 confirmed in scope (VD
2026-08-12); then S4 (fixture + re-record), S5.
