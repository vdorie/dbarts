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
(`memo.md`, `critique.md`,
`synthesis.md`). **Read `synthesis.md` before starting** - it carries the B1
redesign this plan encodes and three extensions neither of the other two
documents has. These are untracked session files, so they do NOT arrive with
a `git worktree add`; the orchestrator copies them into the implementer's
worktree at spawn, and an implementer who cannot find them must ask rather than
proceed without them.

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
   `refuseMultiForestMutation` (`[[C_interface.cpp:330@4c018187]]`), NOT
   `refuseMultiForestTestOffset`, which is R-bridge-only
   (`[[R_interface_bartcore.cpp:2345@4c018187]]`, called at [[R_interface_bartcore.cpp:3959@4c018187]] and [[R_interface_bartcore.cpp:3988@4c018187]]). Both refuse
   the multinomial on the same predicate (`numForests >= 2`) and differ only in
   message, so the decision is unchanged; the earlier spelling of it simply
   named the wrong guard. The category offset must never reach
   `ResponseModel::offset()` (which `storeSample` adds post-blend to every
   location channel).
4. **Case weights stay refused**, on model grounds: a non-integer weight makes
   the augmentation `PG(w_i n_i, .)`, a real-shape draw with no exact sampler
   (the gap `docs/plans/archive/negative-binomial.md` records), and an INTEGER weight is
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
  multiforest arc moved every one of them: `setTreatment` [[R_interface_bartcore.cpp:3295@4c018187]],
  `setForestWeights` [[R_interface_bartcore.cpp:3319@4c018187]], `setOffset` [[R_interface_bartcore.cpp:3782@4c018187]], `setResponse` [[R_interface_bartcore.cpp:3800@4c018187]],
  `setSigma` [[R_interface_bartcore.cpp:3825@4c018187]], `setData` [[R_interface_bartcore.cpp:3832@4c018187]], `setTestOffset` [[R_interface_bartcore.cpp:3951@4c018187]],
  `setTestPredictorAndOffset` [[R_interface_bartcore.cpp:3974@4c018187]], `setWeights` [[R_interface_bartcore.cpp:4037@4c018187]], `setControl` [[R_interface_bartcore.cpp:4059@4c018187]],
  `setModel` [[R_interface_bartcore.cpp:4089@4c018187]], `setPredictor` [[R_interface_bartcore.cpp:4195@4c018187]], `predict` [[R_interface_bartcore.cpp:4884@4c018187]] (whose offset guard
  now lives in the shared `predictFromSource` [[R_interface_bartcore.cpp:4785@4c018187]], predicate unchanged:
  `offset != NULL && numLocations > 1`). Every multinomial combiner anchor is
  intact, line-shifted: `drawForestGlue` [[R_interface_bartcore.cpp:800@4c018187]] (`lastF_ = f` [[R_interface_bartcore.cpp:833@4c018187]]),
  `formForestResponse` [[R_interface_bartcore.cpp:854@4c018187]], `combinedFits` [[R_interface_bartcore.cpp:870@4c018187]], `combinedTestFits` [[R_interface_bartcore.cpp:885@4c018187]],
  `blendSoftmax` [[R_interface_bartcore.cpp:995@4c018187]], `counts_` [[R_interface_bartcore.cpp:1018@4c018187]], `trials_` [[R_interface_bartcore.cpp:1019@4c018187]], the per-sweep scratch
  [[R_interface_bartcore.cpp:1020-1025@4c018187]], `ForestCombiner::supportsResponseMutation` [[R_interface_bartcore.cpp:427@4c018187]]. The
  zero-weight arc's exact-zero snap landed in `BCFForestCombiner::
  formForestResponse` only; the multinomial combiner is untouched.
- **A mutation-time `totalFits` writer the original reading predates.**
  `Chain::revalidateTrees` / `rebuildFitsFromParameters` ([[chain.hpp:1572-1830@4c018187]],
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
  still exactly three at f737702: [[chain.hpp:1125@4c018187]] post-forest-loop
  (pre-`afterCombine`), [[chain.hpp:4368@4c018187]] inside `storeSample` (which runs AFTER
  `afterCombine`), and [[chain.hpp:1443@4c018187]] in `growForestFromRoot`. Any per-column
  lazy refresh scheme leaves the last category stale at all three.
  `combinedTestFits` has one call site, [[chain.hpp:4405@4c018187]].
- `finalizeTotalFits` computes `total = forestY - resid + lastFit`, i.e. the sum
  of the forest's own tree fits. `forestY` will carry `- o_if`, so `totalFits`
  stays offset-FREE. That is the invariant `raw = totalFits + offset` needs.
- `drawForestGlue` sets `lastF_ = f` unconditionally, so a `lastF_ == f` assert
  in `formForestResponse` would be tautological. Add the invariant COMMENT; skip
  the assert.
- `predictFromSavedSampleMulti` and `predictFromCurrentTreesMulti`
  ([[chain.hpp:2282@4c018187]] / [[chain.hpp:2310@4c018187]]) build their own raw slab and call
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
  `docs/plans/archive/multiforest-mutation-gaps.md` (H1 and the capability-probe
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
   ([[chain.hpp:1572-1830@4c018187]]) rewrite `totalFits` at PREDICTOR-MUTATION time,
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
   already is K for multinomial ([[facade.hpp:42-44@4c018187]]) and 1 for every additive
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
   updateScale)` (`[[R_interface_bartcore.cpp:2483@4c018187]]`, declared in
   `R_interface_bartcore_common.hpp`), and the FLAT C entries call the same
   helper (`[[C_interface.cpp:206@4c018187]]` / [[C_interface.cpp:226@4c018187]] / [[C_interface.cpp:238@4c018187]]). So the repair goes inside the
   shared helper's `!supportsResponseMutation` branch, keyed on the NEW
   `supportsCountsMutation` probe, and the counts hint lands on BOTH surfaces.
   That is correct and deliberate: the hint is true on both, and a helper
   whose whole reason for existing is that "the two surfaces cannot state
   different rules" must not be made to state two. The generic text stays
   verbatim for a future non-opt-in coupling, and the `ResponseConduit`
   wording (`fixes its response at creation` / `carries no offset` / `fixes its
   case weights at creation`) is preserved per conduit. Re-pin in this slice:
   `inst/tinytest/test-multi-forest-seam.R` (the multinomial block, currently
   [[C_interface.cpp:231-282@4c018187]], which pins `setResponse` at both `updateScale`, `setOffset` at
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
arm at [[C_interface.cpp:298-310@4c018187]]), so the new `supportsCountsMutation` field is not optional
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
   (`predictFromSavedSampleMulti` [[chain.hpp:2282@4c018187]], `predictFromCurrentTreesMulti`
   [[chain.hpp:2310@4c018187]] - both now `template <typename Columns>` over a borrowed predictor
   view, so the offset is a new function parameter threaded through the
   templates and their `bartcore_predict` / `predictFromSource` callers) take a
   per-call nNew x K offset added to their raw slab BEFORE
   `softmaxLocationMajor`.
2. Bridge: `bartcore_setCategoryTestOffset`; the creation-side test-offset lift
   for the MATRIX form only, on the R BRIDGE entries only. The FLAT entries stay
   refused, by `refuseMultiForestMutation` on `dbarts_sampler_setTestOffset`
   (`[[C_interface.cpp:330@4c018187]]`) - the R-bridge-only `refuseMultiForestTestOffset`
   (`[[R_interface_bartcore.cpp:2345@4c018187]]`) is the guard this slice conditions, not the
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
3. Re-pin the three CI legs in `.github/workflows/equivalence.yaml` - [[R_interface_bartcore.cpp:61@4c018187]]
   (gaussian, `equivalence-a825263.rds --strict-coverage`), [[R_interface_bartcore.cpp:87@4c018187]] (bcf,
   `bcf-equivalence-a825263.rds`), [[R_interface_bartcore.cpp:113@4c018187]] (multinomial) - so the multinomial line
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
   **Record the reservation in `docs/plans/archive/c-api-growth.md`**, not only here, or
   it is lost when the reshape is scoped (next section).

Gate: full tinytest with NO snapshot regenerated (S5 is RNG-neutral surface);
`R CMD check`; `air format --check .` and `lintr::lint_package()`; trio bitwise
against the S4 baseline (`multinomial-equivalence-1027be5.rds`, 10 scenarios,
plus `equivalence-a825263` strict and `bcf-equivalence-a825263`).

## What the reshape must reserve (record in c-api-growth.md)

1. **DECLINED (dbarts-h-reshape, resolved question 3).** Source-shaped
   response and offset parameters. The reshape already re-signs predictor
   entries onto a borrowed `PredictorSource` view; do the same on the
   response side, replacing the bare `const double* y` / `const double*
   offset` with a tagged source able to express at minimum `{ double*
   vector }` and `{ int* counts, size_t numCategories }` for the response,
   and `{ double* vector }` / `{ double* matrix, size_t numCategories }` for
   the offset. A later flat multinomial creation entry then needs a new tag,
   not an ABI break.

   Declined: `dbarts_sampler_create` has no multinomial branch to reach it
   from, so every tag but the vector one would ship as unreachable dead
   surface; secondarily, a tagged struct would also re-sign `setResponse` a
   second time in the same window. Reserved in its place, meeting this
   item's own goal at lower cost: `dbarts_sampler_setCounts(dbarts_sampler*,
   const int* counts, size_t numCategories)` and
   `dbarts_sampler_setOffsetMatrix(dbarts_sampler*, const double* offset,
   size_t numCategories)`, both appends - recorded in
   `docs/plans/archive/c-api-growth.md`, "Reservations closed and opened at the
   reshape".
2. **ADOPTED.** A size-first spec struct for any future flat multinomial
   creation, per the `dbarts_results` `structSize` precedent this arc's
   sibling established.
3. **SUPERSEDED.** Whatever tagged struct the reshape adopts must PRESERVE
   the `refuseMultiForestMutation` guards already on the flat response-side
   entries (the D4 precedent in `docs/plans/archive/multiforest-mutation-gaps.md`).
   bcf-public-surface S3 item 2 RELAXED exactly those guards, by decision.

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
- **The tests/cpp `ConfigSpec` multi-forest fuzz arm.** DISCHARGED at
  `53525f4d` (see Landing notes). This bullet was stale on arrival:
  `8afd6eac` (2026-08-11) had already landed `fuzzRunBCF`/
  `fuzzRunMultinomial` and four multi-forest configs before this text was
  written on a descendant commit, never re-anchored - `ConfigSpec` has
  covered multi-forest samplers since `8afd6eac`. The genuine gap was
  narrower: `Sampler::setCounts`/`setCategoryOffset` had no tests/cpp
  coverage and no invariant tied a multi-forest sampler's combined output
  to its parts. `53525f4d` adds `OP_SET_COUNTS`/`OP_SET_CATEGORY_OFFSET`
  plus the two parts identities I1/I2.
- **A `k3countsswap` equivalence scenario, at the next multinomial
  re-record.** The counts-swap stream is transitively pinned today
  (test-multinomial-counts-mutation.R pins swap == rebuild bitwise on all five
  channels with a non-vacuity arm, and k3counts pins the rebuild against the
  frozen baseline), so a change moving swap draws is caught either by k3counts
  or by the tinytest identity. The single point of failure is a slice that
  DELIBERATELY revises swap semantics and regenerates that identity as
  expected - the transitive pin dissolves silently. A direct scenario costs a
  full baseline re-record whether recorded now or later, so it rides the next
  multinomial re-record rather than forcing one (S4 verifier adjudication,
  2026-08-12).

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
  (S1-S3; `multinomial-equivalence-1027be5.rds` at 10 scenarios from S4 on,
  k3offset the addition). The multinomial nine are
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
- Design artifacts are durable
  (`memo.md`, `critique.md`,
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

S2 LANDED d7d27a3, 2026-08-12. The category offset channel through
the offset indirection as decided: rawFits returns
totalFits.data() off an offset so the null path is today's
pointer; full rematerialization at the fresh-sweep branch (f == 0
or out-of-order) and at the top of combinedFits, per-column
refresh only in the in-order continuation; blendSoftmax and
combinedTestFits untouched; formForestResponse's null branch is
the old statement verbatim; the working response subtracts o[i]
with psi = rawFits - margin folded over offset fits. Internal
creation surface (matrix offset at creation). B2's floors sit at
bridge level via refuseCategoryOffsetTestSurface (anonymous
namespace) at FIVE call sites - creation, setCategoryOffset,
setTestPredictor, setTestPredictorAndOffset, predict - which S3
lifts by deletion (the implementer counted four; the verifier
found five; five is the record). The setCategoryOffset message
half sits in the shared conduit helper's supportsCountsMutation
branch with the seam-file pins edited in place. O11 carries the
offset term and the after-a-per-observation-session arm. The
INT_MAX range check precedes integer coercion on both counts
entrances, pinned. Deviations verified sound: ONE capability probe
(binding decision 2 fixes one field; no combiner can own an n x K
count response without the n x K predictor; the split note lives
on Chain::setCategoryOffset); predict's flat-offset predicate
numLocations > 1 is exactly supportsCountsMutation today;
creation's refusal reworded with the pre-existing pin edited in
place; multinomial handles carry $K, which cannot go stale because
setCounts refuses any K change; floors bridge-level, the flat C
API needing none. Falsifiers: F2 red - the sign-flipped working
response drives the new exact arm to gap 0.5648 vs tol 0.012, and
its FIRST run was a false green from a stale non-preclean install,
caught and re-run (the --preclean discipline is the recorded
lesson); F5 red 4/53 - the lazy per-column refresh fails both O11
same-vintage arms AND both all-zero-offset neutrality arms, since
a zero offset routes through raw_ and exposes the stale column,
proving raw_ rests on O11 exactly as the S1 carry predicted; F4
non-vacuity on both parities with the row-shift invariance stated
as softmax identifiability. The verifier's independent engine-risk
probe: 24 hostile inputs through both the R wrapper and the C
entry all refused with the sampler's subsequent run bitwise
identical to baseline; six set/clear cycles, set-then-setCounts,
three setState interleavings, a whole-matrix setPredictor, a
per-observation session and a two-chain run all hold
softmax(forestFits + offset) == train at gap 0; the offset is
never silently dropped. Budget: 1213 insertions vs ~810 (1.50x on
the insertions metric, 74 percent of the overage oracle content -
tests/cpp 187 vs 50, R/bench 582 vs 380) - accepted at
orchestrator discretion under the mandated-oracle rule. Note to
S5: bartcoreForestFits reports the offset-FREE totalFits, correct
and load-bearing - its surface docs must say so. Carries to S3:
lift the five floors by deleting the refuseCategoryOffsetTestSurface
call sites; the predict message repair is already the
per-category-matrix pointer S3 makes true; rawTestFits/rawTest_
mirror rawFits/materializeRawFits, both private, in place. Carry
to S4: k3offset can reuse recordChannels from
test-multinomial-category-offset.R. Gates double-run (implementer
+ independent verifier, fresh privlibs, --preclean): tests/cpp
plain + ASAN/UBSAN clean, tinytest 4100/0 (53 new) no snapshot
regenerated, trio 35/35 strict / 11/11 / 9/9 vs a825263 no
re-record, multinomial-exact all six arms (the five pre-existing
at exactly the S1-recorded gaps; arm6 1e-4 vs 8e-3), bcf oracles
at MANIFEST values, air 0, lintr::lint_package 0, R CMD check
clean-copy tarball Status OK.

S3 LANDED 1027be5, 2026-08-12. The test offset channel: creation
offset.test on both internal creators (the matrix form on the R
bridge entries, per item 2) + bartcoreSetCategoryTestOffset +
resident semantics; combinedTestFits carries the offset with
rawTest_/rawTestFits/materializeRawTestFits mirroring the train
trio, private, in place, blendSoftmax untouched, rematerialization
EAGER at every report (no install-time caching - the lazy variant
is a demonstrated-red falsifier); the two predict replays thread
the per-call offset through the template Columns forms with the
null path unchanged. All five refuseCategoryOffsetTestSurface call
sites AND the helper deleted, zero references repo-wide; the S2
floor pins inverted in place. dbarts.h and C_interface.cpp zero
diff - the flat refuseMultiForestMutation guard stands (binding
decision 3). No seventh exact arm, adjudicated correct: S3's gate
stanza is "as S2 plus ASAN", the new-arm mandate was S2's O4, and
a test offset enters no likelihood so there is no quadrature
target. The implementer-added refuseStaleCategoryTestOffset guard
was adjudicated STRONGER than a design choice - a memory-safety
requirement: rawTest_ is sized off numTestObservations and
materializeRawTestFits reads nTest*K doubles from the BORROWED
testOffset_ sized at install, so replacing test predictors at a
different nTest would overrun it; all four post-creation resize
sites are guarded and setData is independently refused. The
refusal asymmetry is the PLAN's own floors semantics: replacements
under a resident test offset refuse (clear first), zero-surface
extensions are allowed and PINNED AS INTENDED (all-zero == no
offset bitwise; the no-test-offset recorded channel == the
zero-offset predict), and predict re-imposes the refusal since it
receives both at once. Predict precedence pinned both ways: the
argument always wins and the resident offset is never read.
Deviations accepted: creation-side lift mandated by item 2 with
create-vs-swap parity bitwise; flat-offset-on-multinomial refusals
pinned R-side and by direct .Call; predict-vs-run asserted
BITWISE (expect_identical - the replay path has no RNG).
Falsifiers: F6 red on exactly the four replay-side arms with every
blend-side arm green (verifier re-ran it independently from a
scratch copy: 4 red of 59, same lines), proving the halves
independently gated; the lazy-rematerialization negative half red
on 3 tests/cpp checks; F4 non-vacuity in both new files. DEFECT
CARRIED TO S5 (verifier, low severity, non-blocking):
[[R_interface_bartcore.cpp:5239@4c018187]] - predict's missing-offset refusal
keys on the TRAIN offset only, so a sampler carrying a resident
TEST offset and no train offset accepts a no-offset predict and
returns the well-defined offset-free surface; the plan's wording
is the unqualified "carries a category offset"; fix is one ||
clause plus a pin, folded into S5's surface-message pass. Budget:
1068 insertions / 1254 changed lines vs ~575/~860 (1.86x on
insertions, 1.46x on changed lines) - the verifier mapped the
430-line test file section by section onto the plan's six S3 test
items and the full edge-case list, genuinely mandated; accepted at
orchestrator discretion under the mandated-oracle rule, as S2's
was. Engine-risk probes: every hostile input refused with the
subsequent run bitwise the untouched arm; set/clear cycles,
setCounts, setState round trips, whole-matrix setPredictor,
per-observation session all safe under a resident test offset
(nothing changes nTest; eager rematerialization dominates the
mutation-time totalFits writers). Carries to S4: k3offset takes
the test offset at creation or via the setter; reuse
recordChannels from test-multinomial-test-offset.R. Carries to S5:
the [[R_interface_bartcore.cpp:5239@4c018187]] predicate fix + pin; bartcoreForestFits reports
offset-FREE totalFits (docs must say so); bart2 does not yet
forward a test-side matrix (thread it if test.offset is wanted for
multinomial); no test-side analog of bartcoreForestFits. Gates
double-run (implementer + independent verifier, fresh privlibs,
--preclean): tests/cpp plain + ASAN/UBSAN clean, tinytest 4162/0
(62 new) no snapshot regenerated, trio 35/35 strict / 11/11 / 9/9
vs a825263 no re-record, multinomial-exact all six arms at exactly
the recorded gaps, bcf oracles at MANIFEST values, air 0,
lintr::lint_package 0, R CMD check clean-copy tarball Status OK.

S4 LANDED 2b96a9f, 2026-08-12. The k3offset scenario in
benchmarks/R/multinomial-equivalence.R: K = 3, an n x K train
offset and a 25 x K test offset taken at CREATION, a first run,
then BOTH replaced mid-chain via bartcoreSetCategoryOffset /
bartcoreSetCategoryTestOffset, a second run; five channels
recorded, no verdict channel (offsets install unconditionally);
seeds 6010/6110/7010 are literals outside the guarded seeds
vector, so settingsList() is byte-identical to the a825263
recording and the neutrality compare runs. ONLY the multinomial
baseline re-recorded, as multinomial-equivalence-1027be5.rds -
named by the ENGINE tip per the exact 33f6fdc precedent (HEAD at
recording was 89ddb0f, one docs-only commit later; meta.rev says
so, MANIFEST states it); equivalence-a825263 and
bcf-equivalence-a825263 re-verified at this tip and left current;
[[equivalence.yaml:113@4c018187]] re-pinned, [[equivalence.yaml:61@4c018187]]/[[equivalence.yaml:87@4c018187]] correctly untouched.
Neutrality: 35/35 strict / 11/11 / 9/9 with k3offset printing no
line (verified ran); self-reproduction 10/10; the six exact arms
at the recorded gaps. Non-vacuity, measured at recording and
reproduced independently by the verifier: with both offsets NULL
at the same seeds train moves 0.83 and test 0.77 max-abs and tree
structure moves (both varcount channels differ);
softmax(forestFits + offset) reproduces the final recorded train
sample at machine epsilon while softmax(forestFits) alone is 0.37
away - the offset-free forestFits report is load-bearing, not
redundant. The verifier additionally falsified the BCF harness
against the historical 99205ee baseline (5 MISMATCH in the exact
documented pattern) to prove the compare is live. No NEWS bullet
(the plan assigns S1/S2/S3/S5 only). COVERAGE ADJUDICATION
(verifier, accepted): no scenario drives bartcoreSetCounts
mid-chain; the swap stream is transitively pinned (tinytest swap
== rebuild bitwise + k3counts pins the rebuild), the residual
being only a deliberate contract revision that regenerates the
identity - recorded as the k3countsswap door in "Doors held open",
riding the next multinomial re-record rather than forcing one.
Gates double-run (implementer + independent verifier, fresh
privlibs, --preclean): tinytest 4162/0 unchanged, trio neutrality
+ self-reproduction as above, air 0, lintr 0 new (3 pre-existing
brace lints in the harness), diff confined to harness + baselines
+ MANIFEST + yaml (src/R/inst zero diff).

S5 LANDED 69b28a3e, 2026-08-12 (amended once at verification).
Public surface: bart2's multinomial offset refusal lifted for an
n x K numeric matrix only, threaded to the internal creators -
bitwise identical to the dbarts::: creator at the same seed and
offset, 0.807 from the no-offset fit (verifier functional probe);
offset.test NOT wired, per the plan's item 1 and the synthesis
For-VD note, both naming the train side only. The verifier's one
defect drove the amendment: a train-offset fit with test data
silently reported an offset-free yhat.test, against the arc's
no-channel-silently-omits-the-offset contract - the amendment
DOCUMENTS the asymmetry (bart.Rd: offset is train-side only,
yhat.test is always computed without any category offset, the
test twin is dbarts:::-only; multinomial.md carries the matching
sentence), adds an R-boundary refusal for an explicit offset.test
naming the internal channel (previously fell through to a
bridge-depth message several layers from the caller), and
re-qualifies parseMultinomialData's now-false unqualified "do not
support a test offset" to name the nTest x K channel; no prior
test pinned the old text, and both new messages are pinned.
predict's missing-offset predicate fixed (the S3 carry): the
refusal keys on EITHER resident offset; the verifier probed all
four configurations - neither predicts, train-only, test-only and
both refuse, and an explicit all-zero matrix escapes all three
refusing configurations. bartcoreForestFits documented offset-free
with the reconstruction identity (measured: gap 0 with the offset
added back, 0.451 without). Deviations accepted: the creation-time
weights message improved beyond the plan's literal singular (no
pinned text; accurate on both surfaces); two stale multinomial.md
surface bullets corrected on initiative (the transactional-refusal
and test-offset-at-creation claims, superseded by the multiforest
arc and S3 - verified live by the verifier). Wrong-shape public
offsets fail at the R boundary with named messages, never
bridge-deep. Door candidates recorded: no test-side
bartcoreForestFits analog; the k3countsswap door stands. Gates
double-run plus the amendment re-run: tests/cpp plain green,
tinytest 4174/0 (12 new), trio 35/35 strict / 11/11 / 10/10 vs
the current baselines, the six exact arms at the recorded gaps,
air 0, lintr::lint_package 0, R CMD check clean-copy tarball
Status OK on every round, no new exported topic.

THE ARC IS COMPLETE. Landed: S1 729bbdd (counts channel), S2
d7d27a3 (train category offset + floors), S3 1027be5 (test offset
+ predict replays, floors lifted), S4 2b96a9f (k3offset fixture +
multinomial-equivalence-1027be5), S5 69b28a3e (public surface +
messages + docs); records 99dfd10, dff706b, 89ddb0f, ed1ffb9, and
this commit. The multinomial family now serves as a conditional
inside a larger Gibbs/MH sampler on the same terms as every other
family: counts mutate at fixed n and K, train and test category
offsets install and clear mid-chain, predict replays honor a
per-call offset, no reported channel silently omits the offset,
and the one public surface is bart2's train-side matrix offset.

Fuzz-arm ops LANDED 53525f4d, 2026-08-24. Premise correction: the
"Doors held open" bullet and TODO's entry were half stale - 8afd6eac (2026-08-11)
landed fuzzRunBCF/fuzzRunMultinomial before this door text was
written on a descendant commit and never re-anchored. The real gap:
setCounts/setCategoryOffset had no tests/cpp callers, and no
invariant tied combined output to parts; both existed only as
fixed-scenario R pins. Landed: OP_SET_COUNTS (unit one-hot, grouped
multi-trial, self-swap asserted bitwise inert) and
OP_SET_CATEGORY_OFFSET (install or clear), guarded on
supportsCountsMutation; a multinomial-multichain (2-chain) config
drives the per-chain fan-out. Two parts identities run after every
op: I1, reported channel == softmax of the per-category totals plus
the installed offset, chain-strided, 1e-12 relative (mirrors
softmaxLocationMajor's max-subtraction, ulps apart; absolute 1e-12
vacuous at ~1e-16); I2, an amplitude coupling's fitsWithoutOffset
== the exact-z blend of forest totals, last-forest-down, via
forestCalibration's responseScale/responseShift, 1e-9 relative. I3
(finiteness), I4 (state round trip) pre-existed. Discrimination, by
planting: P-A drop raw_.resize -> SEGV under ASAN
(refreshRawColumn), I3 never runs; P-B conditional
materializeRawFits -> I1 fails, 18; P-C cache glue_.combined -> I2
fails, 24; P-D drop trials_=trials -> NOT caught here; caught by
the R identity pin (test-multinomial-counts-mutation.R) - honest
non-catch; P-E setCategoryOffset ignores a clear -> I1 fails, 12,
SEQUENCE-ONLY (install, clear, run) - the only tinytest offset
clear ([[test-multinomial-r5-surface.R:296@5a3bc276]]) never sweeps after.
Single-forest streams unchanged (opt-in DBARTS_FUZZ_TRACE FNV-1a
digest, 108 lines, 9 configs x 12 seeds, diffed pre/post). Gates:
tests/cpp full green 15.7s; fuzz 0.270s -> 0.281s (3 seeds), 0.604s
-> 0.703s (12 seeds, mostly the 2-chain config); ASAN/UBSAN clean;
R build byte-identical (tests/cpp outside the package), so the
equivalence trio is trivially bitwise, no re-record. Residual door:
bartcore_predictPerForest's forestReportingIsDefined gate still
refuses per-forest off-sample replay on multinomial - a separate
surface item, recorded in TODO.
