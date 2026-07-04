# Public surface for the major version

Status: reviewed 2026-07-03; decisions from that review are recorded inline
as DECIDED. Companion to core-generalization.md: the engine reached cutover
readiness (full R5 parity, statistical equivalence, and the zero-regression
speed bar all gated at 209c09b), so what remains is what the major version
exposes.

Standing decisions this builds on: the next release is a major version with
no backwards-compatibility requirement on internal surfaces; the classic
engine and the `R_C_interface.hpp` ABI work until cutover and are then
retired; weights with binary responses are rejected by design (the classic
weighted probit is incorrect); MATCH_BAYES_TREE is not preserved; RNG-shifting
changes regenerate test snapshots.

## 1. Engine selection and cutover sequencing

Every wrapper already routes through `dbartsControl` (bart, bart2, rbart_vi,
and pdbart construct one internally; xbart accepts one as an argument), so
engine selection needs no new public arguments: flip the `engine` default to
`"bartcore"` and everything follows.

Proposed sequence:

1. Port xbart's backend. Its C++ driver (`R_interface_crossvalidate.cpp`)
   uses `BARTFit` directly and is the one wrapper with classic-only plumbing.
   This wants the shared data handle (section 5) so folds bin against cuts
   computed once, rather than re-sorting per fold.
2. Flip the `dbartsControl` default. Regenerate the hardcoded test snapshots
   (established convention), re-record benchmark baselines.
3. Remove the classic engine, the `engine` argument, and the old ABI in the
   same release. Recommendation: no deprecation cycle for `engine="classic"`;
   the escape hatch is the previous CRAN release, and keeping both engines
   doubles maintenance for a path the equivalence harness already guarantees.
   `LinkingTo: dbarts` consumers (stan4bart) migrate to the new C API
   (section 6) in the same window.

Open: whether any external `LinkingTo` consumers beyond stan4bart exist and
need a transition release with both ABIs. RESOLVED 2026-07-03: the CRAN
reverse dependencies (EBcoBART, bartCause, bartMan, countSTAR, glossa,
riAFTBART, stan4bart, tidytreatment, voi) were checked; only stan4bart
links, the rest are R-level. No transition release needed beyond the
stan4bart port.

DECIDED: break freely now and patch at the end; xbart may be redesigned
around the new data layer rather than accommodated, or shelved and
completed later - it does not gate the rest.

Landed 2026-07-03 (step 3, C++ deletion): src/dbarts/ is gone, along with
R_interface_sampler.cpp, R_interface_common.cpp, R_C_interface.cpp, and
every classic header in inst/include/dbarts/ - dbarts.h is the only
shipped header left. The classic .Call and CCallable tables are gone
(only the dbarts_sampler_* registrations remain); the package unload
finalizer is a no-op (bartcore samplers free through per-pointer R
finalizers). Build system: SUBDIRS/PKG_LIBS drop dbarts.a on both
platforms, configure loses match-bayes-tree, thread-safe-unload, the
dbarts config headers, and the generated inst/include/dbarts/types.hpp;
--with-xint-size STAYS - misc's partition kernels (used by bartcore)
consume XINT_TYPE through src/include/misc/types.h. The cutover is
complete: stan4bart ports to dbarts.h in lockstep (Vincent publishes
both).

Landed 2026-07-03 (step 3, R surface): the `engine` argument is gone from
`dbartsControl` (with the slot, its validity checks, and the class
prototype), along with `rngKind`/`rngNormalKind` (explicit generator
kinds were a classic feature; `rngSeed` stays) and the classic
`dbartsState` class. Every R5 method runs its bartcore body
unconditionally; `startThreads`/`stopThreads` remain as no-ops for their
callers. Binary responses with weights are now refused at creation
everywhere (previously only the bartcore path refused).
test-sampler-customMCMC.R (classic-stream snapshots; behaviors covered on
bartcore in test-bartcore.R) is deleted, test-rng.R keeps its
engine-generic blocks, and the cross-engine comparisons in
test-bartcore.R (prior-sampling KS test, byte-identical verbose creation
summaries) became bartcore-only structural checks. The equivalence
harness still accepts `engine=new` and compares the installed engine
against saved classic baselines; re-recording under classic is no longer
possible, so the existing baseline rds files are the permanent
cross-engine reference. The classic C++ (src/dbarts/, R_C_interface,
call tables) is unreachable and deleted in the next step.

Landed 2026-07-03 (step 2, the flip): `dbartsControl` defaults to
`engine = "bartcore"` and `factors` defaults to `"categorical"` on
`dbarts`/`dbartsData`/`bart2`/`xbart` (which gained the argument); `bart`
stays the frozen shim (indicators, probit, but the new engine). The flip
flushed out and fixed real parity gaps the R5 harness missed: recorded
training fits omitted the offset (diverging rbart's ranef Gibbs),
`rngSeed` was ignored, test data could not be removed (bart2's burn-in
dance), `getLatents` refused preallocated results (rbart's in-place
contract), `setResponse`/`setOffset` skipped length validation
(segfault), `copy()` walked classic state slots, and the state was not
lazily materialized for saveRDS. Fixed residual priors now work on
bartcore with the documented variance semantics. rngKind/rngNormalKind
are refused (classic-only). RNG-locked snapshots regenerated; classic
mechanics tests pin `engine = "classic"` until step 3 removes them.

## 2. Data ingestion: factors as categorical columns

Current behavior: data.frame and formula inputs pass through
`makeModelMatrixFromDataFrame`, which dummy-expands factors; `dbartsData`
types every column ordinal. The engine's categorical rules (canonical-gauge
subset splits, up to 53 levels) are reachable only by flipping
`data@varTypes` by hand.

Proposed (landed 2026-07-03 behind `factors = "categorical"` on `dbarts` and
`dbartsData`; the default remains `"indicators"` until cutover flips it):

- Unordered factor columns ingest as categorical: codes 0..K-1 in the column,
  `varTypes` set, level names retained on the data object for reporting.
  Ordered factors ingest as ordinal on their integer codes (their order is
  meaningful; threshold splits are the right vocabulary).
- Matrix input stays all-ordinal - no change for the most common call.
- This changes fitted models for data.frame inputs with factors, which is the
  point: subset splits search 2^(K-1) - 1 partitions instead of K indicator
  thresholds. Given the major-version license, categorical is the default;
  an escape argument on `dbartsData` (and threaded through bart2/xbart
  formulas), tentatively `factors = c("categorical", "indicators")`, keeps
  the old expansion reachable for comparisons.
- Factors with more than 53 levels: error, naming the cap. No silent
  fallback to dummy expansion (that silently changes the model class).
  Lifting the cap needs pooled mask storage plus a wide-mask story for the
  flat tree format; see section 7.
- Missing values: still rejected at ingestion. The reserved-NA-code +
  missing-direction-bit design (MIA) remains future work and gets its own
  proposal when it lands. (LANDED 2026-07-04 as designed - see
  mia-missingness.md and section 7.)

DECIDED: users supply a data.frame and dbarts builds its own internal
representation (DMatrix-style); whether that representation is visible in
R is an implementation choice, so the handle starts internal. NAs stay
rejected in this pass, but nothing may make handling them hard later:
per-column cut counts cap one below the code type's maximum so a reserved
NA code always fits, and the missing-direction bit already has a home -
the 53-category cap leaves mask bits 53-63 unused and ordinal rules'
high union word is invariantly zero. The flat tree format may grow a
flags field when MIA lands; state objects are opaque, so that is not a
compatibility break.

Landed (2026-07-04): reporting format for categorical rules in
`getTrees`/`plotTree`. The flat format stores the direction mask as a
double and `getTrees` keeps the raw mask in the `value` column; when the
sampler has any categorical predictors the data.frame gains a `directions`
column decoding each categorical rule into one "L"/"R" character per level
in level order (bit set sends the level right, matching the engine), with
ordinal rules and leaves NA. `plotTree` labels a categorical rule with the
level names sent down the left branch ("g in {blue,red}"), reading as the
left-branch condition like the ordinal "<=". The decode lives R-side
(`decodeCategoricalSplits`); the bare `bartcoreGetTrees` helper stays raw.
This is the reporting decision pooled masks (section 7) waited on: a
wider-than-double mask changes only the flat format and the R decode,
not the reported vocabulary.

## 3. Response families

Current: binary responses silently mean probit; logistic (Polya-Gamma) is
internal-only via `dbarts:::bartcoreSampler(sampler, family = "logistic")`.

Proposed:

- `family = c("auto", "gaussian", "probit", "logistic")` on `dbarts()` and
  bart2 (bart keeps its BayesTree-compatible signature and stays
  gaussian/probit-by-response). `"auto"` preserves today's dispatch:
  continuous -> gaussian, binary -> probit. Probit stays the binary default;
  logistic is opt-in (heavier-tailed latents, roughly 1.6x the scale, and
  the established defaults/expectations are probit's).
- Node-scale default under logistic: match the probit convention on the
  latent scale, i.e. leaf prior spans the same +-3 sd after accounting for
  the logistic latent variance. Exact constant to be fixed against the
  exact-posterior gate when the argument lands.
- Weights plus a binary family error with a message explaining that the
  weighted-probit likelihood the old engine implemented was incorrect and
  what to use instead (replicate rows, or a gaussian model on latents).

Open: whether `family` also belongs on `xbart` (its losses already branch on
binary responses; logistic would want log-loss by default).

DECIDED: `family` as proposed; probit stays the binary default.
Landed 2026-07-03: `family = c("auto", "gaussian", "probit", "logistic")`
on `dbarts()` and bart2, resolved against the response into a new
`dbartsControl@family` slot (which also drives pointer re-creation after
save/load). `"gaussian"` on a 0/1 response fits a continuous model;
logistic requires the bartcore engine, with node.scale = pi * sqrt(3) -
probit's 3.0 widened by the logistic latent sd. The R5 surface reports
binary fits on the latent scale. The wrappers' probability transforms
(predict/extract/fitted/plot) went link-aware 2026-07-03: packaged fits
carry a `family` element and transform through it (fits saved without
one are probit). rbart stays probit-only - its ranef Gibbs step assumes
normal latents and `rbart_vi` rejects a `family` argument.

## 3a. Prior specification

The quoted-prior DSL (`tree.prior = cgm(power, base)`, `node.prior =
normal(chi(1.25))`) is re-evaluated in a namespace environment - which is
also why it never polluted the search path: `normal`, `chisq`, `chi`, and
`fixed` are too generic to export (rstanarm exports `normal`; masking bugs
by attach order become support requests).

DECIDED: evolve in place, keeping that no-pollution property.
Landed 2026-07-03: `dbartsPriors` exports the constructors (including
`dart`); the bare-name sugar evaluates against the vocabulary layered over
the caller's environment; specs resolve against data/control at fit time
(named split probabilities, the binary k default, the DART update delay).

- Prior specifications become first-class classed objects built by ordinary
  standard-evaluation constructors, validated at construction. `dart()`
  joins `cgm()`; hyperprior nesting (`normal(k = chi(1.25))`) is plain
  composition; family-aware defaults live on NULL arguments resolved at fit
  time.
- Exactly one new symbol is exported: a container (working name
  `dbartsPriors`) holding the constructors, e.g.
  `dbartsPriors$normal(k = dbartsPriors$chi(1.25))`. No generic names
  enter the search path.
- The bare-name sugar keeps working: fitting functions evaluate their
  prior arguments with the constructor vocabulary layered over the
  caller's environment, so `node.prior = normal(chi(1.25))` resolves
  against dbarts' vocabulary inside those arguments regardless of what is
  attached. A value that is already a prior object passes through, so the
  sugar and the programmatic path compose. NSE-only tricks
  (`split.probs = 1 / numvars`) get documented replacements (NULL means
  uniform over available variables).
- No new front end: `dbarts()` and `bart2` keep their names and signatures;
  `bart` stays the frozen BayesTree-compatibility shim.

## 4. DART exposure

The Dirichlet split-variable machinery is engine-complete (s-update, alpha
grid update, update delay) but only reachable through `splitprobs` overrides
or internal handles.

Proposed: a tree-prior spec alongside `cgm` in the existing quoted-prior DSL:

    dbarts(..., tree.prior = dart(power, base, a = 0.5, b = 1, rho = NULL,
                                  alpha = 1, update.alpha = TRUE))

with `update.delay` defaulting to half the burn-in (the BART-package
startdart convention; a cold forest under sampled alpha is bistable without
it). The fit exposes split-probability samples (`varprobs`, chains x samples
x predictors) next to `varcount`. bart2 gains a `dart = FALSE` convenience
flag mapping to the spec with defaults.

Landed 2026-07-03: `tree.prior = dart(...)` on `dbarts()` (bartcore engine
only; the classic engine and `setModel` refuse it; the sampled
probabilities appear in the stored state). Also landed same day:
per-sample probabilities come back as `varprobs` in run results (engine
Results member, filled only under DART) and on packaged bart2 fits next
to `varcount`; bart2 takes `dart = TRUE` (defaults, with its `power` and
`base`) or a full `dbartsPriors$dart(...)` spec, refusing `split.probs`
alongside it.

## 5. Standalone data handle and CV views

The design doc's remaining phase-4 item. Two consumers are concrete: xbart's
backend (folds sharing one cut grid, no per-fold re-sort) and multi-sampler
setups over shared predictors (BCF-style prognostic + treatment forests, the
IRT embedding).

Proposed:

- A built data object holding the two-layer store (cuts + codes), constructed
  once from a `dbartsData`; external-pointer handle, serializable. Working
  name `dbartsDataHandle`; constructing samplers from it skips cut
  construction and quantization.
- Row-subset views: a sampler over (handle, rows) gathers the subset's codes
  and copies the parent's cut points, so folds bin identically to the full
  data by construction. Views hold no raw predictor matrix, which naturally
  refuses the predictor-mutation surface (setPredictor and friends need raw
  values); response-side mutation stays available.
- Single-writer rule: samplers created over a shared handle treat predictors
  as read-only. Mutable-predictor workflows keep the current path (each
  sampler owns its store).

DECIDED: internal first; naming and serialization contracts are cheaper to
change before exposure. Exposure waits until multi-forest or user demand
needs it.

Landed (2026-07-04, internal): `ColumnStore::buildFromParent` +
`Sampler`/facade view construction; bridge entry points
`bartcore_createDataHandle` (external-pointer handle over a built store)
and `bartcore_createFromHandle` (view sampler; slices y/weights/offset by
trainRows C-side, test offset from offset[testRows]); R helpers
`bartcoreDataHandle`/`bartcoreSamplerFromHandle`. Views are self-contained
copies (no lifetime coupling to the handle) and refuse the raw-x mutation
surface, including setState (restoring cut points re-quantizes from raw
values); setTestPredictors and predict work off the copied grid. xbart
builds one handle per worker chunk, so folds bin on the full data's grid -
its hardcoded regression values were regenerated. Update 2026-07-04:
views compose with linear leaves - buildFromParent gathers the designated
columns' raw values with standardization constants from the parent's full
data (the same calibration inheritance as the copied cut grid), and xbart
gained a node.prior argument (linear-leaves.md). Still open, as decided:
serialization and any public exposure.

## 6. C API and callbacks

`inst/include/dbarts/R_C_interface.hpp` exposes C++ classes (Control, Data,
Model, BARTFit) across a compiled boundary - fragile against layout changes
and gone with the classic engine. bartcore is header-only C++20, which no
CRAN package should be forced to compile against.

Audit of stan4bart 0.0-13 (src/init.cpp, bart_util.{hpp,cpp}), 2026-07-03,
resolving the open question below:

- CCallables called (19): the lifecycle set (initializeControl,
  initializeData/invalidateData, initializeModel/invalidateModel,
  initializeFit/invalidateFit, setControl), runSamplerWithResults,
  sampleTreesFromPrior, the Gibbs conditioning set (setOffset with the
  update-scale flag every iteration, setSigma every iteration, storeLatents
  into its own buffer), predict, createStateExpression/initializeState,
  printInitialSummary, printTrees, getTrees. setResponse is looked up but
  never called (latents flow BART -> Stan, not back).
- ABI couplings beyond the function table: allocates BARTFit by sizeof;
  stack Control/Data/Model with direct field access (keepTrees toggled off
  for warmup and on for sampling, defaultNumSamples/defaultNumBurnIn sized
  to the run, responseIsBinary, verbose, numChains written on restored
  samplers, numTrees read); fit->currentNumSamples to size predict output;
  fit->sharedScratch.dataScale read (exported as range.bart to untransform
  predictions R-side) and written (forced to -0.5/0.5 on samplers restored
  for predict); subclasses Results and bumps its public pointers to land
  one draw at a time in caller-controlled storage; reads FlattenedTrees
  members and frees them field by field; model.kPrior->isFixed to size k
  samples.
- R-level internals (not C ABI, same migration window): dbarts:::parsePriors,
  new("dbartsModel") with node.scale, data@sigma/@weights/@n.cuts,
  control@binary.

Every coupling maps onto the bartcore facade, so the v1 header needs no
engine additions: caller-buffer run kills the Results subclass (bartcore
Results pointers are caller-owned with null-means-skip), setTreeStorage
replaces the keepTrees toggling, savedTreeCapacity replaces
currentNumSamples, kIsSampled/usesDart replace the prior pokes, and the
opaque complete state plus original-scale predict kills both the dataScale
write and range.bart. getTrees returns the data.frame directly - both
stan4bart entry points only marshal it into one. stan4bart inherits the
binary+weights refusal (it forwards user weights verbatim; by-design). The
create-time verbose summary replaces the printInitialSummary entry point.

DECIDED, v1 surface of `inst/include/dbarts/dbarts.h`:

- Flat C ABI: opaque `dbarts_sampler*`, `DBARTS_C_API_VERSION` plus
  `dbarts_apiVersion()`, everything reached via `R_GetCCallable`. Symbol
  names are `dbarts_sampler_*` (camelCase methods after the prefix):
  classic already exports `dbarts_setResponse` et al. until removal, and
  handle-prefixed names stay defensible after it.
- SEXPs at exactly two boundaries: creation takes the R spec objects
  (dbartsControl/dbartsModel/dbartsData - the prior/data language IS R),
  and state/trees are R-serializable/R-facing by purpose. All else is
  plain arrays.
- Entry points: create/destroy; run into a caller-owned `dbarts_results`
  struct (sigma, train, test, varcount, k, varprobs; NULL skips);
  sampleTreesFromPrior/sampleNodeParametersFromPrior; setResponse,
  setOffset(updateScale), setWeights, setSigma, getLatents; setPredictor/
  updatePredictor (transactional, rollback on failure); setTestPredictors
  (NULL removes)/setTestOffset; predict into a caller buffer;
  setTreeStorage, getTrees (data.frame), printTrees; storeState/setState;
  setNumThreads/setNumThin/setVerbose; queries (numObservations,
  numPredictors, numTestObservations, numChains, numTrees,
  numSavedSamples, kIsSampled, usesDart).
- Contracts documented in the header: errors longjmp via Rf_error; R RNG
  bracketing is internal (call from the main R thread, no caller
  bracketing); creation preserves the data object, raw-array setters
  borrow for the sampler's lifetime.
- Additive evolution is free (name lookup), so deferred without cost:
  per-observation predictor updates and the joint session, setCutPoints,
  setData, and the per-iteration observer callback - classic
  `Control::callback` had no reachable consumer, so its replacement waits
  for one.
- The C++ headers ship for `LinkingTo` use without ABI promises (header-only,
  recompiled per consumer), documented as unstable.

RESOLVED by the port (2026-07-04, stan4bart 8c96206 against dbarts
1.0-0): stan4bart needed two dbarts-side additions the audit missed, both
in state semantics rather than the entry-point list. First, the state now
carries the gaussian response transform (698f630): setOffset(updateScale)
moves the scale after creation and nothing else could reproduce it, which
is what actually kills the dataScale write and range.bart. Second,
multi-chain restores may gather single-chain states, with a generator-kind
mismatch leaving the destination stream unrestored (1eee83e): stan4bart
fits chains in separate single-chain samplers riding R's stream and
splices their states into one dedicated-generator sampler for prediction.
The port itself: per-iteration caller buffers replace the Results
subclass, setTreeStorage replaces the keepTrees control toggling, the
stored sampler is created from the spec objects with the control shaped to
the state, getTrees returns the data.frame directly, and the create-time
verbose summary replaces printInitialSummary. stan4bart pins
factors = "indicators" against the new categorical default, and its one
seed-pinned stochastic test was lengthened past trajectory luck. Suite:
206 passing, 0 failing.
before freezing the header. RESOLVED above, 2026-07-03.

Landed 2026-07-03: dbarts.h v1 with all entry points registered as
CCallables; the bridge's creation/state/tree cores extracted into
bartcore_bridge:: (R_interface_bartcore_common.hpp) so C_interface.cpp is
a thin adapter; an actual-consumer test (test-capi.R compiles
inst/tinytest/capi/consumer.c against the installed headers and drives
the stan4bart workout through R_GetCCallable, self-gating on toolchain
availability). Remaining for this section: the stan4bart port, observer
callbacks when a consumer exists, and retiring R_C_interface.hpp with the
classic engine.

## 7. Deferred, with their blockers

- Pooled category masks (K > 53): LANDED 2026-07-04 (design and landing
  notes in pooled-masks.md; the cap is now 65535, masks up to 63
  categories stay inline in the rule word, wider ones pool per tree, and
  the flattened format's masks past 53 categories ride a side channel
  with getTrees reporting value = NA plus the usual directions decode).
- Sparse columns: PROTOTYPED 2026-07-04 (docs/design/sparse-columns.md,
  benchmarks/kernels/sparse.c). The representation is settled - a
  rank-bitmap column layout at parity with the dense kernel where
  sparsity is real, in a tenth the memory - and order-preserving
  partitions are ruled out. Landing waits on the per-column data model
  (u8 widths) plus a dgCMatrix ingestion story.
- MIA missingness: LANDED 2026-07-04 (design and landing notes in
  mia-missingness.md; surface is missing = c("incorporate", "error"),
  incorporate the default).
- Wave-2 models (linear leaves, in-core grouped random effects retiring the
  rbart_vi R loop): engine work, independent of this document except that
  rbart_vi's public signature would eventually gain the in-core option.
  Linear leaves are LANDED in full (linear-leaves.md, 2026-07-04): a
  designated column set per leaf regression via node.prior =
  linear(columns, k) on dbarts() and xbart(), including data-handle views
  (section 5 update). Grouped random effects have a proposal
  (grouped-random-effects.md, 2026-07-04): a ResponseModel decorator
  running rbart_vi's Gibbs blocks in-core for the built-in tau priors,
  with the R loop retained for custom priors and callbacks.
