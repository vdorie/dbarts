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
need a transition release with both ABIs.

DECIDED: break freely now and patch at the end; xbart may be redesigned
around the new data layer rather than accommodated, or shelved and
completed later - it does not gate the rest.

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
  proposal when it lands.

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

Open: reporting format for categorical rules in `getTrees`/`plotTree`. The
flat format stores the direction mask as a double; proposal is to keep the
raw mask in the `value` column and add a decoded convenience (per-level
direction characters, using stored level names) at the R level.

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

## 6. C API and callbacks

`inst/include/dbarts/R_C_interface.hpp` exposes C++ classes (Control, Data,
Model, BARTFit) across a compiled boundary - fragile against layout changes
and gone with the classic engine. bartcore is header-only C++20, which no
CRAN package should be forced to compile against.

Proposed:

- A flat C ABI, `inst/include/dbarts/dbarts.h`: opaque `dbarts_sampler*`,
  `extern "C"` functions over plain arrays, an integer ABI version constant.
  Entry points sized by what the R bridge and stan4bart actually use:
  create/destroy, run (results into caller buffers), the mutation surface
  (response, offset, weights, sigma, predictors incl. transactional
  updates), latents, predict, state save/restore, prior sampling, and the
  control/model setters. Names mirror the R5 methods.
- Callbacks: a per-iteration observer registered with a `void*` closure,
  receiving (chain, iteration, is-burn-in, training fits, sigma). Same
  threading contract as the progress sink: with one thread it runs inline
  and may call into R; with worker threads it runs on the worker and must
  not touch R. This replaces classic `Control::callback`, which no R-level
  consumer could reach.
- The C++ headers ship for `LinkingTo` use without ABI promises (header-only,
  recompiled per consumer), documented as unstable.

Open: whether stan4bart needs anything beyond the list above (it currently
reaches into fit internals in places); worth enumerating against its source
before freezing the header.

## 7. Deferred, with their blockers

- Pooled category masks (K > 53): needs a flat-format extension (masks no
  longer fit a double exactly) and the section-2 reporting decision; state
  objects are opaque and engine-specific, so serialization is free to
  change.
- Sparse columns: kernel design plus an ingestion story; prototype after the
  handle exists.
- MIA missingness: reserved NA code + rule direction bit; interacts with
  ingestion defaults, so it follows section 2 landing.
- Wave-2 models (linear leaves, in-core grouped random effects retiring the
  rbart_vi R loop): engine work, independent of this document except that
  rbart_vi's public signature would eventually gain the in-core option.
