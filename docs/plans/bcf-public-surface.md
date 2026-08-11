# bcf-public-surface

agent: S0 sonnet (tests only). S1 opus (creation routing, an S4 slot, two live
  defects). S2 sonnet (R5 methods, mirroring). S3 opus (the flat C surface and
  the hash re-bake). S4 opus (an engine output channel). S5 sonnet (the fit
  function and its class). S6 sonnet (docs, records). Serialized: one
  implementer, each slice lands before the next starts.
rng: NEUTRAL on every existing path in every slice. The one path with no
  baseline is S1's public BCF creation, which is a NEW stream (it builds one
  engine where `bartcoreBCFSampler` builds and discards two, so under
  `rngSeed = NA` it consumes n.chains fewer R draws) - nothing to re-record.
  S4 fills output arrays and must be GATED neutral, not assumed. Every slice
  gates the trio bitwise; a divergence on the single-forest or multinomial
  baseline is a LEAK, never a re-record.
window: this arc lands INSIDE the pre-release breaking window and BEFORE the
  queued dbarts.h reshape (TODO). Everything is breakable; the four sister
  packages in /Users/vdorie/Repositories migrate in lockstep (VD 2026-08-10).
  Runs after zero-weight-exactness, which OWNS a bcf-equivalence re-record;
  this arc re-records NOTHING and gates against the then-current hash.
budget: S0 ~180 test; S1 ~130 R + ~30 S4 class + ~200 bridge + ~230 test;
  S2 ~180 R + ~230 test; S3 ~60 header + ~90 C_interface + ~30 bridge +
  ~120 capi test; S4 ~50 engine + ~120 bridge + ~30 R + ~190 test;
  S5 ~350 R + ~160 man + ~270 test; S6 ~60 docs. Total ~2350.

## Goal

BCF stops being reachable only through `dbarts:::bartcoreBCFSampler`. A user
writes `dbarts(x, y, treatment = z, ...)` and gets an ordinary `dbartsSampler`
whose R5 methods, `storeState`/`saveRDS` round trip, and re-creation all work;
a C consumer reaches the same sampler through `dbarts_sampler_create` and can
drive it; and `bcf()` fits one and reports per-draw mu, tau and glue. The
engine's BCF model does not change.

## Binding decisions inherited (do not reopen)

1. **Shape: BOTH, layered.** A creation selector on `dbarts()`/`dbartsSpec()`
   plus a `bcf()` front end on top. The two layers serve DISJOINT verified
   consumers - C and R5 mutation drivers below (stan4bart, treatSens, bairrtt),
   bartCause above - and that is the fact that decides the fork. Not
   `family = "bcf"`: `family` selects the RESPONSE model, `dbartsModel@family`'s
   validity enumerates response distributions, and BCF's response is gaussian.
   (State the exception explicitly wherever the rule is written down:
   `family = "multinomial"` selects a K-FOREST sampler, so "family never selects
   forest structure" is false on its face; a multinomial RESPONSE really is
   categorical, a BCF response is not.)
2. **Neither z nor per-forest weights ride the engine state.** No
   `stateFormatVersion` bump. Zero-weight-exactness deferred the per-forest
   weight half to this arc; it is decided here, the same way.
3. **Per-forest weights are NOT folded into a BCF creation spec** - membership
   is resampled per sweep (zero-weight-exactness, unchanged). They stay
   `dbarts:::`-only in this arc.
4. **`update.a`/`update.b` stay CREATION-TIME.** Freezing a glue block mid-chain
   is a model change, not conditioning, and no consumer asks. Door, not built.
5. **`bartcore_createBCF` and `bartcoreBCFSampler` are LEFT UNCHANGED.** The
   bcf-equivalence fixture and five test files drive them. Retiring them is a
   separate arc with its own re-record.
6. **This arc re-records no baseline.** Every gate is bitwise against an
   existing hash.

## Adjudication: corrected facts this plan is built on

The design memo and its blind refuting critique sit in
`.claude/bcf-public-surface-design/`, which is GIT-IGNORED. This file is the
only tracked record and carries the load-bearing facts rather than pointing at
them. Verdict STANDS WITH AMENDMENTS. Every anchor below was re-read against
the live tree at b374d9e; probes were re-run and extended.

- **The flat C mutation entries are guarded and become LIVE.**
  `dbarts_sampler_setResponse`, `setOffset`, `setWeights` and `setTestOffset`
  (C_interface.cpp) each call `refuseMultiForestMutation`, which errors on
  `numForests >= 2`. Their comments say the guard is unreachable "today"
  precisely because `dbarts_sampler_create` cannot reach a multi-forest model.
  S1 makes it reachable. PRECISION the plan carries: no EXISTING consumer code
  breaks - the guards fire only on a sampler the consumer chose to create as
  BCF - but the ADOPTION path is closed until S3 opens it, so the memo's
  "creation alone unblocks treatSens" is RETRACTED.
- **`dbarts_sampler_setResponse` hardcodes `updateScale = true`**, which is
  exactly what `bartcore_setResponse` refuses for a multi-forest sampler ("a
  multi-forest sampler supports a response swap only with updateScale = FALSE,
  which pins the response transform its per-forest leaf calibrations are stated
  against"). The semantics BCF requires are INEXPRESSIBLE through the shipped
  signature. So: no new dbarts.h symbol for CREATION, but the DRIVE loop needs
  new and changed entries, the hash DOES move, and the re-bake is forced.
- **Per-draw mu, tau and glue are obtainable TODAY.** Measured: `run(0, 1)` then
  `bartcoreForestFits(bc, 0L)`, `bartcoreForestFits(bc, 1L)`,
  `bartcoreBCFGlue(bc)` reconstructs the recorded `train` draw through
  `a*mu + b_z*tau` to 3.4e-14 at 1 chain and 2.9e-15 per chain at 4 chains, at
  1.04x the cost of a batched run. `bartcoreGetTrees(forest = 1L)` is a second
  channel. NEW MEASUREMENT, made here and not in either document: a per-sweep
  driver loop is BITWISE IDENTICAL to the batched run (`train` and `sigma`,
  2 chains, 10 samples, max diff exactly 0). So the two routes are
  interchangeable and S4 is ERGONOMICS AND C-REACHABILITY, not a prerequisite.
  "S4 without the reporting channel is impossible" is FALSE and is dropped.
- **The `setControl` attribute drop fails LOUDLY.** `setControl` copies only
  `binary` and `call` and replaces the field wholesale, so `bartcore.*` control
  attributes are dropped - reproduced. But `Chain::setState` refuses on a forest
  block-count mismatch, so re-creation ERRORS ("state is not consistent with
  this sampler") on heteroscedastic and grouped samplers, and without a stored
  state `getPointer()` refuses outright. Measured here for BCF specifically, in
  BOTH directions: a 2-forest state into a 1-forest sampler and a 1-forest state
  into a BCF handle both refuse. Severity is "an unnecessary loud failure after
  a legitimate `setControl`", not "silently changes the model".
- **The null-factory defect is not a leak.** `createBCFSampler` returns nullptr
  when `options.numVarianceTrees > 0`; `createBCFHolder` stores the result with
  no test. The rngs are MOVED INTO the holder and `~BartcoreHolder` destroys
  them, and `holderFinalizer` deletes the holder - so nothing leaks. What
  happens is that a live external pointer is returned wrapping a NULL
  `sampler` unique_ptr, and every entry dereferences `holder.sampler->...`
  unguarded. SEGFAULT, not a leak.
- **The engine state DOES carry data-derived content.** Measured:
  `attr(state, "cutPoints")` is the per-column split grid and FOLLOWS x (scaling
  a column by 10 and shifting by 3 moves that column's grid to [3.10, 12.45]
  and leaves the untouched column's grid identical); `fit.scale` is `range(y)`
  to the last digit. `setState` even takes a `currentPredictors` argument so a
  cross-grid restore can re-quantize. The memo's "the state carries no y, no
  weights, no offset, no x" is FALSE. The decision survives on the CORRECTED
  argument, recorded in "State carriage" below.
- **The `dbartsData`-slot rejection cited a mechanism that does not exist.**
  Unnamed `new()` arguments are matched by CLASS, not by slot POSITION -
  measured: the same four prior objects in reversed order error, and a plain
  two-slot basic-typed class refuses unnamed construction outright. Every
  `dbartsData` slot is basic-typed, so unnamed construction of one is impossible
  in the first place, and the cited receipt constructs a `dbartsModel` anyway.
  This is load-bearing: it is what reopens the slot route below.
- **bairrtt is a fourth consumer, present and causal.**
  `/Users/vdorie/Repositories/bairrtt` (`Imports: dbarts (>= 1.0-0)`), headline
  model `irt_causal_bart`, drives dbarts from R with `$run(0L, 1L)` per scan,
  `$setCutPoints`, `$sampleTreesFromPrior`,
  `dbarts::updatePredictorPerObservationJointly`, `$setPredictor(forceUpdate =
  TRUE)` and `$predict`. Under BCF the joint per-observation update is refused
  (`refuseMultiForestTransactionalUpdate`) and `$predict` is refused
  (`refuseBCFTestSurface`). So multiforest-predictor-mutation and the
  saved-tree-replay door HAVE a named consumer; they are not consumer-less.
- **Only one attribute-carried feature is cross-checked in both directions.**
  `resolveSamplerSpec` attaches FOUR bridge-read control attributes, not two:
  `bartcore.n.categories` (ordinal), `bartcore.dispersion` (nbinom),
  `bartcore.survival` (aft) and `bartcore.variance`. (`bartcore.hazard.periods`
  is attached there but the bridge never reads it; `bartcore.groups` is attached
  by `rbart_vi`, not by `resolveSamplerSpec` - there is no grouped precedent in
  that function to copy.) Of the four, only `applySurvivalAttribute` errors in
  BOTH directions - present without `family = "aft"`, and absent with it.
  Stripping `bartcore.variance` silently changes the model with no warning
  (measured). A BCF marker inherits `bartcore.variance`'s posture unless it is
  given a partner.

## State carriage - DECIDED, on the corrected argument

Neither z nor the per-forest weights enter the engine state.

The argument is NOT "the state carries chain state only" - it does not (see
above). It is: **the state carries DERIVED quantities the engine cannot
recompute from its inputs (the quantized cut grid, the pinned response
transform, tree structure, rng position), never RAW conditioning vectors.** y,
weights, offset and x are all absent, and z is exactly `weights` with a
different name. The `data@weights`/`setWeights` symmetry is untouched and
carries the decision on its own.

Consequences, binding:

- **z rides `dbartsData`, as a slot, not a control attribute.** This DEPARTS
  from the memo (see Departures). `$setTreatment(z)` writes the engine AND the
  slot, exactly as `$setWeights` writes the engine and `data@weights`; the
  mirror is then MECHANICAL rather than a discipline, `setControl` cannot drop
  it by construction, and it is what `dbarts_sampler_create` already carries
  (data is already a creation argument, so the ABI does not move).
- **The treatment-forest CONFIGURATION rides `attr(control, "bartcore.bcf")`** -
  tree count, base, power, moderator column mask, `sd.control`, `sd.moderate`,
  `b.prior.variance`, `update.a`, `update.b`, per-forest interactions and
  blocks. `bartcore.variance` is the exact precedent: a second forest's tree
  count, structure prior and column subset on a control attribute.
- **The bridge cross-checks the two in BOTH directions.** Config present with no
  treatment vector, or a treatment vector with no config, is an error naming the
  missing half. This is the mechanism F7 needs; the memo pre-registered the
  falsifier without one.
- **Per-forest weights stay OUT of any spec.** They gain an R5 field
  re-applied after re-creation, on the same rule. While the channel is
  `dbarts:::`-only against a raw handle the hazard is unreachable; the rule is
  recorded now and enforced when the channel goes public.

HONEST LIMITS, stated so a reader does not have to find them:

1. A state moved to a DIFFERENT sampler (`installForests`, cross-sampler
   `setState`) carries no z; the destination's z governs. Correct - z is data -
   but a warm start across two different treatment vectors is silent, not
   refused. It is not silently WRONG: forest counts must still match.
2. BCF state continuation is STRUCTURAL, not bitwise (the shipped `test-bcf.R`
   pins forest fits at `tolerance = 1e-5` after a round trip). No falsifier in
   this arc may assert bitwise continuation.
3. The control now carries the O(p) moderator mask and the config list; the data
   carries an O(n) vector it already carries three of. No new cost class.

## Verified context and seams (read at b374d9e)

- `createHolder` (R_interface_bartcore.cpp) already applies three optional
  control attributes - `applyGroupAttribute`, `applySurvivalAttribute`,
  `applyVarianceAttributes` - and its null-factory path is the teardown
  template: the holder does not yet exist, so the rngs are destroyed by hand
  before `Rf_error`. `createBCFHolder`'s is NOT that case (see the defect).
- `createBCFHolder` already parses everything S1 needs: `bcfParams` (length 8),
  the moderator mask, per-forest interactions and blocks. S1 FACTORS that parse
  into an `applyBCFAttributes` reading the control attribute, and calls it from
  `createHolder`. `createBCFHolder` keeps its own arguments and keeps working.
- `optionsFromParsed` fills `splitProbabilities`, `monotoneDirections`,
  `numLeafCovariates` and `gpLeaves`; the BCF chain constructor ignores all
  four, nulls `maxNumCutsPerVariable` and `forestColumns`, `createBCFSampler`
  instantiates only `SamplerFacade<ConstantGaussianLeaf>`, and `buildBCFForest`
  hardcodes `useDart = false`, `updateK = false`, `k = 1.0`. Every one of those
  is silently dropped today and becomes an explicit refusal in S1.
- The single-forest chain constructor wraps `GroupedResponse`; the BCF
  constructor does not. A `bartcore.groups` attribute reaching a shared
  `createHolder` would silently produce an UNGROUPED fit. S1 REFUSES the
  composition explicitly rather than relying on the attribute's current
  unreachability.
- `bartcore_setTreatment` is the validation and ownership template: capability
  probe FIRST (`Chain::bcfGlue`, never a forest count - a K-forest multinomial
  defeats a count test), then length, then copy into `holder.ownedTreatment`,
  then hand the engine the borrowed `.data()`. `PROT_COUNT` is a fixed enum;
  copy, retain nothing.
- `bartcore_setResponse` and `bartcore_setOffset` carry the identical two-part
  multi-forest condition (`supportsResponseMutation` AND
  `updateScale == FALSE`); `bartcore_setWeights` carries the first half only
  (there is no scale to pin). S3 factors that predicate into ONE helper both the
  bridge and the flat API call, so the two can never diverge.
- `createChainRngs`: with `control@rngSeed` set, a dedicated Mersenne generator
  hands each chain its seed and R's stream is untouched; without one, exactly
  one `unif_rand()` per chain. This is the premise of F1, both halves.
- `Results` (chain.hpp) has no per-forest fits channel and no glue channel;
  `varianceFits`/`varianceTestFits` are the precedent - a SEPARATELY-typed
  forest channel, explicitly "not a numReportedLocations widening", null unless
  a variance forest is present.
- `refuseBCFTestSurface` refuses `setTestPredictor`, `setTestOffset`,
  `setTestPredictorAndOffset` and `predict`; `refuseMultiForestMutation` refuses
  `setData` and `setModel`; `refuseMultiForestTransactionalUpdate` refuses the
  transactional and per-observation predictor sessions. `Chain::revalidateTrees`
  revalidates `forests_[0]` only, which is why.
- `DBARTS_C_API_LIST` (inst/include/dbarts/dbarts.h) is the single source; the
  prototypes, the dispatch table and the FNV-1a token are all expansions of it,
  and `static_assert(dbarts_fnv1a(DBARTS_C_API_DECLS) == DBARTS_C_API_HASH)` in
  C_interface.cpp fails dbarts's own compile until the literal is re-baked. The
  house return convention in that list is **1 = accepted, 0 = refused**
  (`dbarts_sampler_setPredictor`, `dbarts_sampler_getLatents`); the size_t
  probes (`numObservations`, `numChains`, `numTrees`) carry no error channel.
- Design docs: `docs/design/bcf.md`, `docs/design/forest-combiner.md`,
  `docs/design/model-space-survey.md` D3. Related plans:
  `docs/plans/forest-combiner.md` (the mandated BCF oracle set),
  `docs/plans/zero-weight-exactness.md`, `docs/plans/c-api-growth.md`.

## Open decision (VD) - naming only

The public spellings become the release compatibility contract and cannot be
renamed after it. Everything else in this plan is mechanically derived.

**What is being named.** (i) the creation selector and its companions on
`dbarts()`/`dbartsSpec()`; (ii) the `dbartsData` slot and the `dbartsData()`
argument that carry z; (iii) the DSL constructor for the treatment forest's
knobs; (iv) the `bartBCF` element names.

**The alternatives.**

- **A (recommended): `treatment =`, `moderators =`, `treatmentForest(...)`;
  slot `data@treatment`; elements `mu`, `tau`, `glue`, `mu.hat.obs`,
  `mu.hat.cf`, `sigma`, `varcount`.** Descriptive, collision-unlikely, and
  `treatment`/`moderators` are the words the causal literature and the `bcf`
  package use. Cost: `treatment` is a common column name, so a user with a
  `treatment` column in a data frame reads the argument twice; `moderators`
  presumes the reader knows tau's covariate restriction is what is meant.
- **B: `z =`, `tau.covariates =`, `tauForest(...)`.** Matches the engine's own
  vocabulary (`setTreatment(z)`, forest 1 = tau) and `bcf.md` throughout. Cost:
  `z` is opaque at the R surface and collides with nothing but reads as nothing;
  it is the worst of the three for a user who did not read the design doc.
- **C: `causal =`, `causalForest(...)`.** Reads best in a call. Cost: it
  violates the no-generic-names rule (`docs/design/public-surface.md`), and
  `causalForest` is taken in the ecosystem (grf).

Fourteen BCF knobs do NOT become fourteen `dbarts()` arguments under any
option: they ride ONE constructor, the `interactions()`/`blocks()` precedent, so
`dbarts()` grows exactly three arguments. `bcf()` itself takes the flat
arguments - its whole signature is BCF - and builds the spec internally, so
there is one resolver and one validator.

**Recommendation: A.** What would change it: a preference for engine-vocabulary
symmetry over reader-facing description, which argues B.

**RESOLVED: A, FINAL (VD 2026-08-10).** VD's ratification was conditional -
"I assume that your recommendation matches user expectations. If so, you can
proceed" - and the condition holds by construction: `treatment` and
`moderators` are the words the causal-inference literature and the reference
`bcf` package put in front of users, which is the ground the recommendation
was made on. The public spellings are the option-A list above.

**AMENDED (VD 2026-08-10, after S2 landed): the option-A creation surface is
PROVISIONAL.** VD: "I don't think `bcf` belongs in the `dbarts` function" -
the causal vocabulary (`treatment =`, `moderators =`, `treatmentForest()`)
is SCHEDULED FOR REMOVAL from `dbarts()`/`dbartsSpec()` once BCF finds its
proper home ("we can schedule it for removal when we find it a proper home";
the multiforest-extension-surface design arc owns the home question -
candidates include bartCause's reserved "bcf" slot and a declarative
engine-vocabulary K-forest surface). Removal sequences AFTER a replacement
creation route exists, since `treatment =` is today's only public R creation
route. Not a rollback: the S1/S2 mechanism (spec resolution, the
`bartcore.bcf` attributes, the bidirectional cross-check, the R5 methods)
survives re-skinning; only the public argument names move. The flat C names
S3 shipped (`setTreatment`, `bcfGlue`) are likewise renameable at the queued
dbarts.h reshape re-bake if the settled surface uses engine vocabulary. The
opening sentence of this section is thereby amended: the option-A spellings
are the surface UNTIL the home lands; what survives to release becomes the
contract.

## S0. Pin the current surface. No engine change.

tinytest asserting today's BCF behavior one by one, so S1's routing cannot
silently open one: `setData`, `setModel`, `predict`, `setTestPredictor`,
`setTestOffset`, transactional `setPredictor` and the per-observation session
all REFUSE with their current messages; forced whole-matrix `setPredictor`,
`setTreatment`, and the scale-pinned `setResponse`/`setOffset`/`setWeights` all
SUCCEED; `setResponse(updateScale = TRUE)` REFUSES. Plus the driver-loop
identity as a pinned fact: `a*mu + b_z*tau` reconstructs `train` per chain, and
a per-sweep loop is bitwise identical to the batched run.
rng: NEUTRAL (tests only). Gates: full tinytest.

## S1. Creation through the spec surface. THE architectural slice.

**R.** `dbartsData` gains one slot (prototype NULL) and `dbartsData()` one
trailing argument; validity accepts NULL or a length-n vector of 0/1.
`resolveSamplerSpec` gains the BCF branch: resolve the moderator columns (lift
the resolution from `bartcoreBCFSampler`), build `attr(control,
"bartcore.bcf")`, install z on the data object. Three new `dbarts()` and
`dbartsSpec()` arguments and the `treatmentForest()` constructor. Every
overridden option REFUSES rather than being silently dropped: DART,
`split.probs`, `monotone`, a linear or GP node prior, a k hyperprior, grouped
random effects, `variance =`, `storage = "single"`, per-column `n.cuts`, `k`,
`node.scale`.

**Bridge.** `applyBCFAttributes` - the `applyVarianceAttributes` twin - factored
out of `createBCFHolder`'s existing parse; `createHolder` reads the treatment
slot off the data object and routes to `createBCFSampler`; the BIDIRECTIONAL
cross-check (config without z, z without config); the missing null-factory test;
the explicit refusals mirroring the R layer (direct-API consumers get the same
answer). Correct the stale 5-argument `createBCFHolder` declaration in
R_interface_bartcore_common.hpp, whose Doxygen contract describes only
`bcfParams` and `z` while the definition takes ten - a declared-but-undefined
overload that becomes a LINK ERROR the moment `createHolder` tries to call it.

rng: NEUTRAL on every existing path (the new attribute readers return early when
absent). The new creation path is a new stream with no baseline.
Gates: `R CMD INSTALL --preclean` into a PRIVATE library; `tests/cpp` from
clean, plain AND ASAN; full tinytest with NO snapshot regenerated; the trio
BITWISE (`equivalence-c8f661a` 27/27, `multinomial-equivalence-ec2a3d0`, and
bcf-equivalence at whatever hash zero-weight-exactness left current);
`air format --check .`; `lintr::lint` on every touched R file; rchk next
scheduled run.
ABORT: any trio divergence.

## S2. The R5 BCF surface.

`$setTreatment(z)` writing the engine AND `data@treatment`; `$getForestFits`,
`$getBCFGlue`, `$getForestVariableCounts`; BCF-specific messages on the refused
methods; `getPointer()` re-creation works unchanged because the slot rides the
object. `$setControl` PRESERVES `bartcore.*` control attributes - SEPARATE
COMMIT, since it fixes ordinal, nbinom, aft, grouped and heteroscedastic too,
and its severity is "an unnecessary loud failure", not a silent model change.
rng: NEUTRAL. Gates: as S1 minus `tests/cpp`, plus a `saveRDS`/`readRDS`
round trip.

## S3. The flat C surface. ONE slice, ONE hash re-bake.

Every dbarts.h change this arc needs lands here, together, because each re-bake
forces a lockstep consumer recompile.

1. **Widen** `dbarts_sampler_setResponse` to
   `(dbarts_sampler*, const double* y, int updateScale)`, matching
   `dbarts_sampler_setOffset` and the R bridge exactly. Replacement, not a
   parallel name: pre-release reshape-by-replacement is the sanctioned mode
   (the `dbarts_results` precedent in `docs/plans/c-api-growth.md`, where a
   parallel-symbol proposal was superseded once it was clear nothing on CRAN
   links the unreleased header).
2. **Relax** the guards on `setResponse`, `setOffset` and `setWeights` from
   `refuseMultiForestMutation` to the bridge's own predicate, factored into ONE
   shared helper so the two surfaces cannot diverge. `setTestOffset` KEEPS its
   refusal - BCF refuses the whole test surface, so relaxing it would open a
   path the engine does not define.
2b. **Guard `dbarts_sampler_predict` and `dbarts_sampler_setTestPredictors`
   with `refuseBCFTestSurface`** (amendment applied 2026-08-10 from
   dbarts-h-reshape's adjudication, before this slice started). S1 makes BCF
   flat-creatable and thereby makes a SILENT WRONG ANSWER reachable:
   `Sampler::predict` routes to `predictColumns`, whose
   `Chain::predictFromSavedSample` / `predictFromCurrentTrees` both open
   `const Forest& forest = forests_[0]` and loop `forests_[0].numTrees`, so a
   flat BCF consumer receives mu(x) labelled as the fit. The R bridge already
   guards all four of its siblings (`bartcore_predict`,
   `bartcore_setTestPredictor`, `bartcore_setTestOffset`,
   `bartcore_setTestPredictorAndOffset`); of the flat siblings only
   `setTestOffset` is guarded, and only incidentally, by
   `refuseMultiForestMutation` (C_interface.cpp:306). This is NOT two lines.
   `refuseBCFTestSurface` is defined INSIDE the anonymous namespace of
   `R_interface_bartcore.cpp`: not declared in
   `R_interface_bartcore_common.hpp` and not reachable from `C_interface.cpp`.
   The one guard C_interface can already see, `refuseMultiForestMutation`, is
   the WRONG predicate - it fires on `numForests >= 2`, while
   `refuseBCFTestSurface` fires on `numForests >= 2 &&
   !shape.testFitsAreDefined`, deliberately gated on `testFitsAreDefined` so a
   future flat-creatable multinomial (whose test blend IS defined) passes
   through. Using it here would over-refuse exactly the case the engine
   comment says must be allowed. Mechanics, the five steps commit 7299b8b took
   for `refuseMultiForestTransactionalUpdate`: (i) move the definition down
   out of the anonymous namespace into the `bartcore_bridge` block, beside
   `refuseMultiForestMutation`; (ii) append one sentence to its comment -
   "External linkage: the flat C API reuses this guard on its own predict and
   test-predictor entries."; (iii) copy the doc comment as a declaration into
   `R_interface_bartcore_common.hpp` inside `namespace bartcore_bridge`;
   (iv) add `using bartcore_bridge::refuseBCFTestSurface;` at the top of BOTH
   `R_interface_bartcore.cpp` (so its existing unqualified call sites still
   resolve) and `C_interface.cpp`; (v) call it at
   `dbarts_sampler_setTestPredictors` and `dbarts_sampler_predict`, naming
   each entry. Budget: ~35 bridge/common-header + ~25 consumer.c. Extend the
   capi tests with two legs: on a flat-created BCF, `predict` and
   `setTestPredictors` each refuse with the BCF message. NEGATIVE HALF: swap
   the guard for `refuseMultiForestMutation` and a flat multinomial (when one
   exists) must go red - until then, assert the predicate directly in
   `tests/cpp` on a shape with `numForests >= 2 && testFitsAreDefined`.
3. **Append**, at the END of the X-list:
   `size_t dbarts_sampler_numForests(const dbarts_sampler*)` (no error channel,
   matching its `numChains`/`numTrees` siblings; always >= 1);
   `int dbarts_sampler_setTreatment(dbarts_sampler*, const double* z)`;
   `int dbarts_sampler_forestFits(const dbarts_sampler*, size_t forest,
   double* out)` (numObservations x numChains);
   `int dbarts_sampler_bcfGlue(const dbarts_sampler*, double* out)`
   (3 x numChains). The three `int` entries return **1 on success, 0 on
   refusal** - the shipped convention, not the inverse.
4. Re-bake `DBARTS_C_API_HASH` and LEAVE BOTH VERSION CONSTANTS ALONE (VD
   2026-08-10: "No need to increment versions - no bartcore version has been
   released"; the constants stay 1.0 until the first release, when whatever
   they read becomes the initial contract; the hash is the lockstep
   stale-library signal). Document the creation route in
   `dbarts_sampler_create`'s Doxygen.
5. `inst/tinytest/capi/consumer.c` exercises creation, the scale-pinned
   response swap, `numForests`, `setTreatment`, `forestFits`, `bcfGlue`, and
   the four surviving refusals.

**Coordination with the queued dbarts.h reshape.** The two arcs touch DISJOINT
entries: this one widens `setResponse` and appends four names; the reshape
replaces the PREDICTOR entries with source-shaped signatures. This arc lands
FIRST, so the reshape re-bakes over a list that already contains these entries
and adopts them verbatim. RESERVE for the reshape, do not build here: forest-
indexed variants of `setTreeStorage`, `getTrees`, `printTrees`, `predict` and
`numTrees`, all of which are forest-blind today and ambiguous on a BCF sampler.

rng: NEUTRAL. Gates: as S1, plus `test-capi.R`.
ABORT: any trio divergence; a re-signed or appended X-list that leaves
`DBARTS_C_API_HASH` unchanged; or any movement in `DBARTS_C_API_MAJOR` /
`DBARTS_C_API_MINOR` (VD 2026-08-10 - no increments pre-release).

## S4. Per-draw per-forest reporting. Ergonomics and C reach.

`Results` gains `forestFits` (numObservations x numForests x numSamples) and
`glue` (3 x numSamples), on the `varianceFits` precedent: null by default,
filled only for a combiner that defines them, so single-forest and multinomial
allocate and compute nothing. Bridge and R plumbing; the flat C reader appended
in S3 already exists, so nothing here moves the hash.

**This slice is NOT a prerequisite for S5** and must not be described as one.
Its real justification: one `.Call` per run instead of one per sweep; a channel
a C consumer can reach without owning an R sampling loop; and not forcing a
front end to own the loop. Measured cost of the alternative: 1.04x, bitwise
identical draws.
rng: NEUTRAL - filling an output array consumes no rng and touches no state -
and that is GATED, not assumed. Gates: as S1, plus the mandated BCF oracle set
(`bcf-exact.R` quick, `bcf-exact-restricted.R`) UNMOVED.

## S5. `bcf()` and the fit class.

`bcf(formula|x, y, treatment, ...)` returning `bartBCF` carrying per-draw mu,
tau, glue, sigma, per-forest varcount, and the two counterfactual surfaces
`a*mu + b_z*tau` and `a*mu + b_{1-z}*tau`. S3 methods:
fitted/extract/predict/print/residuals.

**The output contract is set by bartCause, and it is the arc's biggest
migration cost - an output cost, not an API one.** bartCause's estimator layer
is written against two counterfactual SURFACES, not against tau: it builds
`x.test` as the z-flipped training matrix, extracts train and test, assigns
`mu.hat.obs`/`mu.hat.cf`, and computes common support from the two surfaces'
posterior sds. A `bcf()` returning only "the treatment effect" is unusable to
it. Under BCF both surfaces follow from mu, tau and glue with NO test matrix at
all - the counterfactual is free.

`predict` on NEW rows returns COMPONENTS (mu(x), tau(x)), never a blended
surface: the test-fit refusal stands, and a component predict needs the
per-forest saved-tree replay, which is a DOOR.
rng: NEUTRAL. Gates: as S2 plus `R CMD check --as-cran`.

## S6. Documentation and records.

`docs/design/bcf.md` (status; the treatment slot and the control attribute; and
the CORRECTION that bartCause drives from R over the sampler mutation API - it
makes ZERO sampler-mutation calls anywhere in `R/`, verified by grep, and every
sibling in its response-method switch dispatches to a FIT FUNCTION);
`docs/design/model-space-survey.md` D3 closed and its stale anchors;
`docs/design/public-surface.md`; `docs/plans/c-api-growth.md` (the reshape
reservations from S3); `inst/NEWS.Rd`; the Landing note in this file.

## Pre-registered falsifiers

- **F1 (S1), the load-bearing one, BOTH halves.** With `control@rngSeed` SET,
  `createChainRngs` is deterministic and independent of R's stream, so the
  public path and `bartcoreBCFSampler` receive identical chain seeds from
  identical (control, model, data). Public BCF creation must therefore reproduce
  the internal path BITWISE on all six bcf-equivalence channels (both forests'
  fits, glue, sigma, train, varcount). A divergence means the spec resolution or
  the calibration input differs. NEGATIVE HALF, mandatory: with `rngSeed = NA`
  the two paths must DIFFER, because the internal path builds and discards a
  host engine that consumes n.chains `unif_rand()` draws first. Write both
  halves or an implementer chases the wrong bug.
- **F2 (every slice).** The trio is bitwise. The single-forest and multinomial
  baselines are LEAK detectors: any divergence aborts the slice.
- **F3 (S2), both halves.** A BCF `dbartsSampler`, after `$setTreatment(z2)`
  with z2 = 0 everywhere, `storeState()`, `saveRDS`/`readRDS`, `getPointer()`,
  reports a combined fit with no tau contribution. NEGATIVE HALF: with the
  mirror removed the re-created sampler uses the creation z - assert the
  difference is OBSERVABLE, so the test proves the mirror rather than the
  plumbing. Continuation is STRUCTURAL; do not assert bitwise.
- **F4 (S1).** Every overridden option refuses at creation rather than being
  ignored, one assertion each: DART, `split.probs`, `monotone`, a linear or GP
  node prior, a k hyperprior, grouped random effects, `variance =`,
  `storage = "single"`, per-column `n.cuts`, `k`, `node.scale`. A silently
  accepted one is the failure mode this arc exists to prevent.
- **F5 (S1), REWRITTEN.** The memo's version asserted "leaks no rng" under
  ASAN; there is no leak to find, so it passed vacuously and would keep passing
  after a bad fix. Correct form: a BCF creation the factory refuses must raise
  an R CONDITION, and the handle must never be usable. Drive it by requesting a
  variance forest alongside a treatment forest (the one live nullptr return) and
  assert `expect_error`; then assert that no external pointer escaped by
  checking the error is raised before any handle is bound.
- **F6 (S2), REWRITTEN.** The memo's version asserted a silent difference on aft
  and grouped samplers that does not occur. Correct form: `setControl` followed
  by `getPointer()` re-creation SUCCEEDS and yields the same z on a BCF sampler,
  and the same variance channel on a heteroscedastic one - i.e. the fix turns a
  loud failure into a success. Pin the PRE-fix behavior in the same file as the
  exact error message, so the commit shows the inversion. The aft leg is
  UNVERIFIED and stays a comment: `dbarts(Surv(t, s) ~ x)` is refused by the
  formula interface and the matrix path's two-column response is undocumented.
- **F7 (S1), the shape's own risk, NOW WITH A MECHANISM.** A control whose
  `bartcore.bcf` attribute has been stripped must produce a LOUD refusal, never
  a silent single-forest fit. Three legs: (i) config stripped, treatment slot
  intact -> the bridge cross-check errors naming the missing config; (ii)
  treatment slot cleared, config intact -> errors naming the missing treatment;
  (iii) `setControl` then re-creation -> succeeds after S2, and before S2 fails
  loudly on the state block mismatch (measured in both directions). The memo
  pre-registered this falsifier with no mechanism; the bidirectional
  cross-check is the mechanism.
- **F8 (S4).** `forestFits[, 0, s] * a_s + forestFits[, 1, s] * b_{z,s}` equals
  the recorded internal-scale train draw for every sample s, and matches what
  the driver loop reports for the same seed. If it does not, the channels are
  not the quantities the front end claims.
- **F9 (S5).** `bcf()` with `rngSeed` set reproduces the equivalent
  `dbarts(treatment = z, ...)` + `run()` sequence bitwise on mu, tau, glue and
  sigma - the front end adds no draw of its own.
- **F10 (S3).** On a BCF sampler created through `dbarts_sampler_create`:
  `setResponse(y, 0)` succeeds and `setResponse(y, 1)` refuses; `setOffset` with
  `updateScale = 0` succeeds and with `1` refuses; `setWeights` succeeds;
  `setTestOffset` refuses. Each assertion is run from `consumer.c`, not from R,
  because the point is that the flat API and the bridge now agree.

## Migration costs per consumer (enumerated, never constraining)

Verified at stan4bart bartcore@6ce0440, bartCause dbarts-1.0@695c603, treatSens
worktree dbarts-1.0@bb9a121, bairrtt main@6167423.

- **stan4bart.** Existing code: ZERO functional change; it must RECOMPILE and
  add `, 1` to its `dbarts_sampler_setResponse` calls if any (it creates via
  `dbarts_sampler_create` and builds the SEXP triple with `dbartsSpec()`). To
  ADOPT BCF: its recorded multilevel-BCF blocker ("it needs a new dbarts.h entry
  point") dissolves for CREATION and is ANSWERED by S3 for driving - it passes
  `updateScale = 1` during warmup today and must pass 0 on a BCF sampler.
  Residual, stated plainly rather than as a dissolved blocker: `setTreeStorage`,
  `getTrees`, `printTrees` and `predict` are forest-blind, and `predict` on a
  BCF sampler errors; grouped x BCF is refused, though stan4bart supplies random
  effects through its own parametric half via `setOffset`, so it does not need
  in-engine grouping.
- **treatSens.** Existing code: ZERO functional change plus the recompile. To
  ADOPT BCF: its grid mutates the RESPONSE at four sites through
  `dbarts_sampler_setResponse`, so it needs S3, not S1 alone - the memo's
  "creation alone unblocks treatSens" was wrong. Each site gains an explicit
  `updateScale` argument - 0 for a BCF, but its EXISTING gaussian grid needs
  **1** at all four sites to preserve behavior (the current entry hardcodes
  true), noting three of the four are post-burn-in calls where dbarts.h's own
  advice is `false`. It also hand-assembles the model with
  `dbarts:::parsePriors` and a positional `methods::new("dbartsModel", ...)`;
  a BCF spec reachable through `dbartsSpec()` lets it retire both `:::`
  reach-ins in the same lockstep. It inherits the `data@sigma` precalibration
  contract, which the BCF map shares.
- **bartCause.** Existing code: ZERO (pure R, no LinkingTo). To ADOPT BCF:
  uncomment the reserved `#bcf = redirectCall(...)` slot in its response-method
  switch and write a `bcf` response fitter against S5's output contract. Second
  cost: its propensity-score-as-covariate option becomes pihat as a MU-ONLY
  column, so it needs a moderator EXCLUSION, not just a covariate injection. No
  weights cost - its weighted-probit strip does not apply, BCF is gaussian.
- **bairrtt.** Existing code: ZERO. To ADOPT BCF for its outcome model: BLOCKED
  by `refuseMultiForestTransactionalUpdate` and `refuseBCFTestSurface`, i.e. by
  multiforest-predictor-mutation and the saved-tree-replay door, not by anything
  in this arc. This arc is what makes those doors consumable at all.

Four altitudes: 23 C entry points and mutation every sweep (stan4bart); 10 and
mutation every sweep (treatSens); ~7 R5 methods with per-observation predictor
mutation (bairrtt); zero and never (bartCause). Only the layered shape serves
all four.

## Defects owned here

1. **A handle wrapping a null sampler (S1).** `createBCFHolder` stores
   `createBCFSampler`'s return with no test; the factory returns nullptr when
   `options.numVarianceTrees > 0`. Nothing leaks - the rngs are moved into the
   holder and its destructor, reached through `holderFinalizer`, destroys them.
   What escapes is a live external pointer whose `sampler` unique_ptr is null,
   and every entry dereferences it unguarded: SEGFAULT on first use.
   Unreachable today only because `createBCFHolder` never calls
   `applyVarianceAttributes`; S1's shared routing makes it reachable, and S1
   closes it in the same commit. Severity: memory-unsafe, user-reachable after
   S1, must not land open.
2. **`setControl` drops `bartcore.*` control attributes (S2, separate commit).**
   Reproduced. Consequence measured: re-creation ERRORS on heteroscedastic and
   grouped samplers, and on a BCF sampler in both directions, because the state
   block counts no longer match; without a stored state `getPointer()` refuses
   outright. No silent path was constructible. Severity: an unnecessary loud
   failure after a legitimate `setControl` - user-reachable, not a release
   blocker, and not a silent model change.
3. **A stale 5-argument `createBCFHolder` declaration (S1).** The header
   declares five parameters with a Doxygen contract describing only `bcfParams`
   and `z`; the definition takes ten. Harmless today because nothing calls the
   5-argument form; a LINK ERROR the moment `createHolder` does. Fixed in the
   same commit that factors the parse.

## Doors held open (recorded, not scheduled)

- **Per-forest saved-tree replay** (out-of-sample mu(x) and tau(x)), needed by
  `predict.bartBCF` on new rows AND by bairrtt. The multinomial K-forest replay
  is the precedent. NOT consumer-less.
- **Per-observation and transactional predictor mutation on a multi-forest
  sampler** - multiforest-predictor-mutation, whose named consumer is bairrtt.
- **A test treatment vector**, which would make BCF's blended test fits defined
  and retire `refuseBCFTestSurface`. A modeling decision, not plumbing.
- **Per-forest `setModel`** (base/power per forest mid-chain).
- **Mutable `update.a`/`update.b`.** Two engine lines plus a bridge entry when a
  consumer pulls; it would immediately inherit the mirroring obligation.
- **Grouped x BCF** (stan4bart's multilevel BCF). The combiner composes with
  `GroupedResponse` by design; the BCF chain constructor does not build the
  decorator. S1 REFUSES the composition; building it is its own arc.
- **Retiring `bartcoreBCFSampler`** in favour of the public path, with the
  bcf-equivalence re-record it forces. Own arc.
- **Public `setForestWeights`**; inherits the mirroring rule on exposure.
- **Forest-indexed `setTreeStorage`/`getTrees`/`printTrees`/`predict`/
  `numTrees`** - reserved for the dbarts.h reshape, not built here.
- **Continuous treatment** (out of scope).
- **`bcf()` over a data handle / row-subset view.**

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S1: `dbarts()` and `dbartsSpec()` accept a treatment vector, which fits the
  two-forest Bayesian causal forest instead of a single ensemble. The sampler is
  an ordinary `dbartsSampler`: it stores, saves, and re-creates like any other.
- S2: BCF samplers gain treatment mutation and per-forest fit, glue and variable
  count accessors; `setControl` no longer discards a sampler's model
  configuration.
- S3: the C API gains a forest count, treatment mutation and per-forest fit and
  glue readers, and `dbarts_sampler_setResponse` takes an explicit
  `updateScale` argument matching `dbarts_sampler_setOffset`.
- S4: a run reports per-forest fits and the treatment glue for every draw.
- S5: `bcf()` fits a Bayesian causal forest and reports posterior draws of the
  prognostic and treatment surfaces and both counterfactual surfaces.

## TODO edits at landing

- `bcf-public-surface`: replace the entry with the landing record - plan doc
  path, artifacts under `.claude/bcf-public-surface-design/`, the slices landed,
  and the naming VD chose. Correct its own stale anchor `R/bartcore.R:536` to
  `bartcoreBCFSampler (R/bartcore.R)`, and drop the "Shape fork at design time"
  sentence, which is answered.
- `multiforest-predictor-mutation`: its two anchors are STALE -
  `R_interface_bartcore.cpp:1909-1923` is now
  `refuseMultiForestTransactionalUpdate`, and `chain.hpp:1484` is now
  `Chain::revalidateTrees`. Also strike "The BCF half is only consumable behind
  bcf-public-surface" once this arc lands.
- `docs/design/model-space-survey.md` D3: same stale `R/bartcore.R:536`; mark
  the prerequisite CLOSED.
- The dbarts.h reshape entry: record that this arc already widened
  `dbarts_sampler_setResponse` and appended four BCF entries, that the reshape
  adopts them verbatim, and that the forest-indexed predictor/tree entries are
  reserved for it.
- The release checklist: add "re-verify all four sister packages against the
  re-baked header" alongside the existing stan4bart line, and note the version
  constants stay 1.0 until the first release, when whatever they read becomes
  the initial contract (VD 2026-08-10).

## Departures from the memo and the critique (record)

1. **z rides `dbartsData`, not a control attribute** - AGAINST the memo. Its
   stated reason for rejecting the slot (positional S4 construction) was
   measured not to exist: unnamed `new()` arguments match by CLASS, every
   `dbartsData` slot is basic-typed so unnamed construction is impossible, and
   the cited receipt constructs a `dbartsModel`. The slot makes the memo's own
   `setWeights`/`data@weights` symmetry mechanical instead of disciplinary,
   removes z from `setControl`'s reach by construction, gives F7 a mechanism,
   and does not touch the ABI (data is already a creation argument). The
   treatment-forest CONFIG stays on the control, where `bartcore.variance` is
   the precedent - it is forest configuration, not data.
2. **The state-carriage argument is re-based, the conclusion kept.** The memo's
   "the state carries no y, no weights, no offset, no x" is refuted by
   measurement (`cutPoints` follows x; `fit.scale` is `range(y)` exactly). The
   correct form is derived-not-raw, plus the weights symmetry, which stands
   alone.
3. **The reporting channel is DEMOTED from prerequisite to ergonomics** - with
   the critique, and further than it went. New measurement made here: a
   per-sweep driver loop is BITWISE identical to a batched run, so `bcf()`
   built either way reports the same draws from the same seed. The arc may
   therefore be split at S4 with no correctness consequence.
4. **The critique's A6 is OVERTURNED.** It asked that
   `dbarts_sampler_numForests` gain an error channel because "0 is not a valid
   forest count either". The defect was the memo's SENTENCE ("all four return
   nonzero on refusal"), not the signature: `numForests` cannot fail, and its
   siblings `numObservations`/`numChains`/`numTrees` all return `size_t` with no
   error channel. Kept as `size_t`, and the sentence corrected in the other
   direction too - the house convention is 1 = accepted, 0 = refused, the
   inverse of what the memo wrote.
5. **The critique's A8 first item is OVERTURNED.** It reported citation drift on
   `combiner.hpp:636/445/250` -> "live 637/446/251". At b374d9e those three
   constructs are at 636, 445 and 250 exactly; the memo's numbers were right.
6. **B1's consequence is SHARPENED, not adopted verbatim.** The critique wrote
   that treatSens and stan4bart "would ERROR". No existing consumer code errors:
   the guards fire only on a sampler the consumer chose to create as BCF. What
   B1 correctly kills is the ADOPTION claim and the "enabling half" offer, and
   its conclusion - a new dbarts.h symbol IS needed and the hash DOES move - is
   adopted in full.
7. **B3's severity downgrade is adopted and EXTENDED to BCF**, which neither
   document could test: a 2-forest state into a 1-forest sampler and a 1-forest
   state into a BCF handle both refuse loudly, measured. So F7's `setControl`
   leg is loud for BCF too, and the genuinely silent window is fresh creation
   from a stripped control - which departure 1 closes.
8. **A4 is adopted and CORRECTED.** The critique fixed the memo's three
   attachment sites to two; the live count is FOUR - `bartcore.n.categories`,
   `bartcore.dispersion`, `bartcore.survival`, `bartcore.variance`. The
   precedent for "`resolveSamplerSpec` attaches a bridge-read control attribute"
   is therefore stronger than either document credited, and the F7 census is
   sharper: only `applySurvivalAttribute` cross-checks in both directions.
9. **A1, A2, A3, A5, A7, B2, B4, B5, B6, B7 adopted** as stated, with the
   amendments above folded into the slices, the falsifiers and the migration
   census.

## Landing notes

S3 LANDED 1622eb9 (2026-08-10). The flat C surface in one commit, one hash
re-bake: setResponse widened in place to take updateScale;
setResponse/setOffset/setWeights route through the NEW shared
bartcore_bridge::refuseMultiForestResponseMutation(conduit, updateScale) -
one predicate for both surfaces, refusal texts byte-identical (S0 pins pass
unchanged); refuseBCFTestSurface promoted to the bridge block per the
7299b8b five-step template and guarding dbarts_sampler_predict +
setTestPredictors (amendment 2b - the silent mu(x)-as-fit hole closed);
numForests/setTreatment/forestFits/bcfGlue appended at list end, int
entries 1 = accepted / 0 = refused; DBARTS_C_API_HASH 0xf760898d116cb3a3 ->
0x1a911c00bb26dcd7, both version constants untouched. The hash is NOT
script-generated (no generator exists; the C_interface static_assert is the
machinery) - the implementer validated its recomputation by reproducing the
old baked value from the pre-edit header first. F10 implemented as
capi_bcf_surface: 12 legs run under R_tryCatchError inside consumer.c,
accept/refuse outcome AND message text checked in C; the 2b negative half
asserted in tests/cpp/test_shape.cpp on both a BCF shape (condition fires)
and a multinomial shape (does not), since the guard cannot link into that
build. Unpredicted facts: dbartsSpec() leaves control@n.samples NA and
dbarts_sampler_create refuses that, so a flat consumer hand-building a spec
must set n.samples; dbarts.h has no hand-written stubs block - stubs and
registration are both DBARTS_C_API_LIST expansions, so the list edit
carried them. Budget 489+/90- (~1.36x; moved comment blocks + the 12-leg
harness). Gates double-run (implementer + orchestrator independently):
install --preclean fresh privlib, tests/cpp plain + ASAN/UBSAN from clean,
tinytest 3900/0, test-capi.R 70 -> 84, trio bitwise (27/27 identical draws
+ bcf 5x7 + multinomial 3x5), air 0, lintr clean on test-capi.R.

S2 LANDED 339aeb0 (2026-08-10; R and tests only, zero src delta, so the
tests/cpp leg was correctly skipped). The R5 surface: $setTreatment mirrors
data@treatment (F3 both halves: a getPointer() re-creation reconstructs
train with the MIRRORED z, not the creation-time z); $setControl now
PRESERVES bartcore.* control attributes (F6 both legs, pre-fix loud
failures quoted in test comments; the heteroscedastic leg keeps its
variance channel); BCF-specific messages on setData/setModel/transactional
setPredictor/setResponse+setOffset(updateScale=TRUE) - the last extended
beyond the plan's list by symmetry, accepted; accessors left un-wrapped
where the bridge message already names the capability, accepted; no S0 pin
inverted (the pins drive the low-level handle, a different path). Trio
bitwise, tinytest 3886/0, air 0, lintr 0, pkgdown clean (new aliases ride
the existing dbartsSampler-class page). Budget 0.71x R / 0.86x test.
Double battery (implementer + orchestrator independent).

S1 LANDED a1dbde7 (2026-08-10; includes the _pkgdown.yml reference entry the
orchestrator added at review - the new-exported-Rd-topic lesson's fourth
occurrence, caught pre-push this time; pkgdown::check_pkgdown clean). Public
`dbarts(treatment =, moderators =, treatmentForest =)` under the option-A
names; the decisive falsifier PASSED - public creation is BITWISE IDENTICAL
to `bartcoreBCFSampler` on all six channels (train, sigma, varcount, mu
fits, tau fits, glue) under a shared seed, and differs unseeded; 17 refusal
assertions fire (the pre-registered eleven plus proposal.probs, Student-t
residuals, and a test set - each silently dropped before; and non-gaussian
family); the null-factory segfault is fixed (factory result tested, clean
R error, no handle escapes) with the variance refusal firing first so the
nullptr is defensive; the stale 5-arg createBCFHolder declaration corrected.
Deviations accepted: resolveModerators EXTRACTED to R/model.R and called
from bartcoreBCFSampler (the lift the plan asked for; behavior proven
unchanged by the 5x7 bitwise gate); Rd docs added in-slice (an exported
function without Rd breaks R CMD check long before S6); the third public
argument is spelled treatmentForest= (argument name == constructor name, the
interactions/blocks precedent). Budget: R side ~1.85x (air one-arg-per-line
inflation + the early Rd work), bridge/test near budget. Gates double-run
(implementer + orchestrator independently): install --preclean, tests/cpp
plain + ASAN from clean, tinytest 3853/0 (S0 pins passing UNCHANGED), trio
bitwise (27/27 + 3x5 + 5x7), air 0, lintr 0 on all touched R files,
pkgdown check clean.

S0 LANDED 994c161 (2026-08-10; tests only, src proven untouched by git
status, so the trio was correctly not run). inst/tinytest/
test-bcf-mutation-pins.R, 19 assertions: all eight refusals pinned with
their current messages (setData, setModel, predict, setTestPredictor,
setTestOffset - the last never previously pinned on a BCF sampler -
transactional setPredictor, the per-observation session, and
setResponse(updateScale = TRUE)); all five allowed mutations pinned
(forced whole-matrix setPredictor, setTreatment, scale-pinned
setResponse/setOffset/setWeights); both halves of the driver-loop
identity pinned (scale*(a*mu + b_z*tau) + shift reconstructs train to
~9e-16 per chain at 2 chains; a per-sweep loop is bitwise identical to a
batched run). Gates: tinytest 3810/0 (implementer) and the new file
19/0 re-run independently by the orchestrator against a fresh install;
air 0; lintr 0. No deviations from the pin list.
