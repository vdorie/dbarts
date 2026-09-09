# Response-model feature matrix

Status: living reference, updated in place whenever a cell changes.

Amended by [pure-c-header](../plans/pure-c-header.md#pure-c-header): the flat C header creates no sampler and
no longer declares the predictor, test-data, weight, active-row, per-forest, state,
tree-extraction or augmentation entries - each is a method on the R sampler object the
handle is now read from. The `retired:` cites below name constructs that are gone; what
this record says about the R and engine sides still holds.

What each shipped response model can and cannot do, one row per model and one column per
capability that bears on scheduling; every SHIPPED and REFUSED cell carries a cite verified
against the live tree.

## Legend

| code | meaning |
|---|---|
| `S` | SHIPPED. Works today; the cite is the site that makes it work. |
| `R` | REFUSED on model or identification grounds - the refusal is part of the model, not a hole. The cite is the refusal site. |
| `M` | MISSING. Not built, no schedule. A cite, when given, is a guard that errors *because the thing is unbuilt* - a recorded open item, not a model refusal. |
| `-` | N/A. The concept does not apply to this row; the row footnote says why. |

No cell carries `P` (planned) or `?` (unverified), and `[fN by family]` flags a cell whose
code varies with the base family the row is built over.

Cites are markdown links whose text is the symbol and whose target is the file, and they are
existence-checked by `tools/check-doc-freshness.R`; a cell's VALUE is adjudicated separately
from its cite.

## Rows

Twelve rows. Six are response models proper, reached through the engine's own `ResponseFamily`
enum ([`ResponseFamily`](../../src/bartcore/model.hpp): gaussian, probit, logistic, aft, ordinal, nbinom); the other six
are reached some other way, so they need rows rather than an enum read. Leaf models (constant,
monotone, linear, GP) are an orthogonal axis, not rows; where one gates a capability the cell
or rule says so.

| key | model |
|---|---|
| gaussian | Gaussian (`ResponseFamily::gaussian`, [`GaussianResponse`](../../src/bartcore/model.hpp)) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, [`TResponse`](../../src/bartcore/model.hpp)) |
| probit | Binary probit ([`ProbitResponse`](../../src/bartcore/model.hpp)) |
| logistic | Binary logistic, weights = observation counts ([`LogisticResponse`](../../src/bartcore/model.hpp)) |
| ordinal | Ordered categorical, cumulative probit ([`OrdinalResponse`](../../src/bartcore/model.hpp)) |
| nbinom | Negative binomial, positive-integer dispersion ([`NBResponse`](../../src/bartcore/model.hpp)) |
| multinom | Multinomial softmax, K forests ([`MultinomialResponse`](../../src/bartcore/model.hpp) + combiner) |
| aft | AFT survival, log-normal ([`AFTResponse`](../../src/bartcore/model.hpp)) |
| hazard | Discrete-time hazard (person-period sugar, [`expandDiscreteTimeHazard`](../../R/dbarts.R)) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, [`bart2Hurdle`](../../R/bart.R)) |
| bcf | K-forest amplitude family, bcf's two forests being its K = 2 instance ([`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)) |
| hetero | Heteroscedastic variance forest ([`buildVarianceForest`](../../src/bartcore/chain.hpp)) |

## 1. Structural signature

Five facts a bridge or engine predicate decides for every row, whichever R entry point reaches
it: a live case-weight channel ([`familyCarriesNoWeights`](../../src/R_interface_bartcore.cpp),
[`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp)); a drawn rather than pinned sigma - pinned by a variance
forest owning the scale, or by the family's own definition ([`sigmaIsPinned`](../../src/R_interface_bartcore.cpp),
[`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp)); a persisted per-observation augmentation vector, a
`latents()` override ([`bartcore_getLatents`](../../src/R_interface_bartcore.cpp)); a non-trivial `fitScale`/`fitShift` at
creation, which makes `updateScale = TRUE` a re-anchor rather than a no-op ([f7]) or a refusal
([`refuseVarianceForestScaleUpdate`](../../src/R_interface_bartcore.cpp)); and a combined out-of-sample fit defined at all
([`refuseUndefinedTestFits`](../../src/R_interface_bartcore.cpp)).

| model | case weights | sigma | latents | unit-scale transform | test fits defined |
|---|---|---|---|---|---|
| gaussian | S [`GaussianResponse::setWeights`](../../src/bartcore/model.hpp) | S [`bartcore_setSigma`](../../src/R_interface_bartcore.cpp) | - [`bartcore_getLatents`](../../src/R_interface_bartcore.cpp) | S [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp) | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| student | S [`TResponse::setWeights`](../../src/bartcore/model.hpp) | S [`bartcore_setSigma`](../../src/R_interface_bartcore.cpp) | S [`TResponse::latents`](../../src/bartcore/model.hpp) | S [`TResponse::setOffset`](../../src/bartcore/model.hpp) | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| probit | R [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp) | R [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | S [`ProbitResponse::latents`](../../src/bartcore/model.hpp) | - [f7] | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| logistic | S [`LogisticResponse::setWeights`](../../src/bartcore/model.hpp) [f8] | R [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | S [`LogisticResponse::latents`](../../src/bartcore/model.hpp) | - [f7] | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| ordinal | R [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp) | R [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | S [`OrdinalResponse::latents`](../../src/bartcore/model.hpp) | - [f7] | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| nbinom | R [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp) | R [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | S [`NBResponse::latents`](../../src/bartcore/model.hpp) | - [f7] | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| multinom | R [`parseMultinomialData`](../../src/R_interface_bartcore.cpp) [f9] | R [`refuseCountsMutation`](../../R/bartcore.R) [f9] | R ["reports nothing, by a DECIDED decline"](multinomial.md) | R [`refuseCountsMutation`](../../R/bartcore.R) [f9] | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| aft | R [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp) | S, no variance forest [`bartcore_setSigma`](../../src/R_interface_bartcore.cpp) / R, variance forest [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | S [`AFTResponse::latents`](../../src/bartcore/model.hpp) | S, no variance forest [`AFTResponse::setOffset`](../../src/bartcore/model.hpp) / R, variance forest [`refuseVarianceForestScaleUpdate`](../../src/R_interface_bartcore.cpp) | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |
| hazard | as probit [f5] | as probit | as probit | as probit | as probit |
| hurdle | - [f10] | - | - | - | - |
| bcf | S, gaussian/logistic [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) [f11] | S, gaussian only [`bartcore_setSigma`](../../src/R_interface_bartcore.cpp) [f11] | S, probit/logistic only [`Chain::latents`](../../src/bartcore/chain.hpp) [f11] | R [`refuseAmplitudeMutation`](../../R/bartcore.R) | R [`refuseUndefinedTestFits`](../../src/R_interface_bartcore.cpp) |
| hetero | S, gaussian [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) / R, aft [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp) [f13] | R [`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) | -, gaussian [`bartcore_getLatents`](../../src/R_interface_bartcore.cpp) / S, aft [`AFTResponse::latents`](../../src/bartcore/model.hpp) [f13] | R [`refuseVarianceForestScaleUpdate`](../../src/R_interface_bartcore.cpp) | S [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) |

Mutation channels and row subsetting read off the table: `setWeights`/`setSigma`/`getLatents`
follow their columns; `updateScale = TRUE` follows the unit-scale column (a no-op where the
transform is fixed, refused where a re-anchor would break a coupling's or a variance forest's
calibration); test predictors, offsets and `predict()` follow the test-fits column;
`setResponse`/`setOffset`/`setPredictor` (+ per-observation) are open wherever a sampler
exists, except multinomial blocks the first two ([f9]) and hurdle has none. Three exceptions
follow from no column:

- Zero-weight subsetting follows the case-weights column except logistic (and bcf's logistic
  sub-case): its weights are positive-integer Polya-Gamma trial counts, so a zero is refused
  although the channel is open ([`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp)).
- Whole-data `setData` (n free) is narrower than `setPredictor`: it needs a sampler that owns
  its predictors ([`refusePredictorMutation`](../../src/R_interface_bartcore.cpp) refuses a data-handle view and a CSC-built
  design), and is refused above one forest ([`refuseMultiForestMutation`](../../src/R_interface_bartcore.cpp)) and, since a
  replacement may change n, for aft too
  (["fix the censoring structure at creation"](../../src/R_interface_bartcore.cpp) on that conduit; the status itself moves
  through `$setResponse(y, status = )`).
- Hurdle's `-` cells are not "ask the two components": `bart2()` refuses `weights`, `subset`
  and `offset`/`offset.test` on that family at its own entrance
  (["does not support 'weights'"](../../R/bart.R)).

Two capabilities are universal, not columns: the mid-chain active-rows mask
([`Chain::setActiveRows`](../../src/bartcore/chain.hpp), every family except hurdle, multinomial's being GLOBAL,
[The contract](active-rows-mask.md#the-contract)), and `extract(type = "loglik")` (every row,
including hurdle's composed density and bcf's combined-fit score; multinomial's engine-side
channel stays undefined, see Gaps). Named calibration (`$getCalibration`/`$setCalibration`, a
per-forest `prior.scale`, [2. The surface](nameable-calibration.md#2-the-surface)) is open on
every single-forest sampler, hetero included with its variance forest not counted
([`buildVarianceForest`](../../src/bartcore/chain.hpp)), and refused on both couplings, whose leaf scale comes from a
calibration map ([f11]).

## 2. Reach

`xbart()` and the flat C API, cited by the token or by the refusal/absence. The flat column
reads DRIVABILITY rather than construction: `dbarts.h` declares no creation entry, so every
sampler named here is built in R and the handle a flat entry takes is the address in that
object's external pointer.

| model | `xbart()` | flat C `dbarts.h` |
|---|---|---|
| gaussian | S [`xbart`](../../R/xbart.R), [`gaussian`](../../R/xbart.R) | S [`DBARTS_FAMILY_GAUSSIAN`](../../inst/include/dbarts/dbarts.h) |
| student | M [`xbart`](../../R/xbart.R) | S [`parseSamplerSpecification`](../../src/R_interface_bartcore.cpp), [`residualDf`](../../src/R_interface_bartcore.cpp) [f2] |
| probit | S [`xbart`](../../R/xbart.R), [`probit`](../../R/xbart.R) | S [`DBARTS_FAMILY_PROBIT`](../../inst/include/dbarts/dbarts.h) |
| logistic | S [`xbart`](../../R/xbart.R), [`logistic`](../../R/xbart.R) | S [`resolveFamily`](../../src/R_interface_bartcore.cpp), [`logistic`](../../src/R_interface_bartcore.cpp) |
| ordinal | R [`resolveClassificationFamily`](../../R/data.R) | S [`resolveFamily`](../../src/R_interface_bartcore.cpp), [`ordinal`](../../src/R_interface_bartcore.cpp) [f3] |
| nbinom | M [`xbart`](../../R/xbart.R) | S [`resolveFamily`](../../src/R_interface_bartcore.cpp), [`nbinom`](../../src/R_interface_bartcore.cpp) [f3] |
| multinom | R [`resolveClassificationFamily`](../../R/data.R) | M [f4] |
| aft | M [`xbart`](../../R/xbart.R) | S [`DBARTS_FAMILY_AFT`](../../inst/include/dbarts/dbarts.h) |
| hazard | M [`xbart`](../../R/xbart.R) | M [f5] |
| hurdle | M | M [f10] |
| bcf | M [`xbart`](../../R/xbart.R) | S [`dbartsSpec`](../../R/spec.R), [`forests`](../../R/spec.R) |
| hetero | M [`xbart`](../../R/xbart.R) | S [`applyVarianceAttributes`](../../src/R_interface_bartcore.cpp) [f3] |

Construction reach through `bart()`, `bartBT()` and `dbarts()` + R5: gaussian, student, probit,
logistic and aft reach all three ([`bart`](../../R/bart.R), [`bartBT`](../../R/bart.R), [`dbarts`](../../R/dbarts.R)).
Ordinal, nbinom, multinomial, hurdle, hazard, bcf and hetero are out of `bartBT()`'s reach
entirely: it carries no `family` formal, so the by-name refusals it once held are gone
(retired: [`refuseBartOwnClassFamily`](../../R/bart.R), [`refuseBartRedirectedFamily`](../../R/bart.R)), and its
one response refusal is the categorical one at [f1]. `bart()` ships all six: ordinal ([`bart2Ordinal`](../../R/bart.R)), nbinom
([`bart2Negbin`](../../R/bart.R)), multinomial ([`bart2Multinomial`](../../R/bart.R)), hurdle ([f10]), hazard as
person-period sugar ([f5]), bcf and hetero through `forests =` / `variance =` ([f6]).
`dbarts()` + R5 ships every one but hurdle, which it refuses
(["two independent samplers"](../../R/dbarts.R)): ordinal and nbinom
([`dbarts`](../../R/dbarts.R), [`ordinal`](../../R/dbarts.R), [`dbarts`](../../R/dbarts.R), [`nbinom`](../../R/dbarts.R)), multinomial on the matrix
interface only ([f4]), and hazard, bcf and hetero by the same routes.

## 3. How a fit is built

[`dbartsSpec`](../../R/spec.R) resolves the seven single-forest tokens - auto, gaussian, probit,
logistic, aft, ordinal, nbinom - takes `family = "multinomial"` directly
([`dbartsSpec`](../../R/spec.R), [`multinomial`](../../R/spec.R)), reaches the K-forest amplitude family through
`forests =` ([`dbartsSpec`](../../R/spec.R), [`forests`](../../R/spec.R), each declared `forest(basis = ...)`) and a
variance forest through `variance =` ([`dbartsSpec`](../../R/spec.R), [`variance`](../../R/spec.R)); only hazard and hurdle
stay out of its reach. A `forests =` fit resolves gaussian, probit or logistic only, aft,
ordinal and nbinom being refused by name at the R layer, the bridge and the factory alike
(["a treatment forest does not support family"](../../R/spec.R),
[`refusedAmplitudeFamilyReason`](../../src/R_interface_bartcore.cpp), [`createAmplitudeSampler`](../../src/bartcore/facade.hpp)); `bart2()` reaches the
same machinery through a `forest()` formula term, under an identical gate ([f6]).

## 4. Composition rules

A variance forest requires `family = "gaussian"` or `"aft"`
(["a variance forest requires family"](../../R/spec.R), [`varianceForestIsRefused`](../../src/bartcore/facade.hpp)) -
heteroscedastic IS that capability, so its row carries no variance-forest cell. The other four
single-forest families - probit, logistic, ordinal, nbinom - are refused because each already
routes precision through its own latent channel, the one a variance forest divides into; aft
has no such channel, so under aft each censored latent is instead redrawn at its own s(x_i)
([`AFTResponse::refreshLatents`](../../src/bartcore/model.hpp)). Even under either family it refuses Student-t residuals
(["does not support Student-t residuals"](../../R/spec.R),
[15. Post-landing: Student-t residuals refused (2026-08-17)](heteroscedastic.md#15-post-landing-student-t-residuals-refused-2026-08-17)) and monotone constraints
(["not supported with monotone constraints"](../../R/spec.R)). It never takes DART either
([`buildVarianceForest`](../../src/bartcore/chain.hpp) leaves `useDart` at its default false).

DART is refused for both multi-forest couplings by name - bcf at
["a DART tree prior"](../../R/spec.R), multinomial at ["'dart' or a DART 'tree.prior'"](../../R/bart.R) -
and each coupling's forest builder hard-sets `forest.useDart = false` whatever the route asked
for ([`buildSpecifiedForest`](../../src/bartcore/chain.hpp), [`buildMultinomialForest`](../../src/bartcore/chain.hpp)). Warm start and
grow-from-root are refused, as `M` rather than `R`, for the four alternate-family `bart2` arcs
and for the multi-forest donor warm start ([f12]); grow-from-root itself ships, covered at two
forests, and is gated by the LEAF model rather than the family - linear and GP leaves are
refused in [`growFromRoot`](../../R/dbarts.R) and a no-op in [`growForestFromRoot`](../../src/bartcore/chain.hpp), so every family
reads "constant leaf".

## 5. Combiners and couplings

Two rows sit on a COMBINER, an object holding K forests plus the rule that combines their fits
into one location per observation ([`ForestCombiner`](../../src/bartcore/combiner.hpp),
[The ForestCombiner<L> hierarchy](forest-combiner.md#the-forestcombinerl-hierarchy)): multinomial's is a softmax over K
category forests, bcf's the AMPLITUDE family, forest f's fit at row i scaled by
`dot(a_f, B_f(i, .))` ([The amplitude layout](multiplier-combiner.md#the-amplitude-layout)). A
K-forest chain takes its response model from `AmplitudeSpec::family`
([`AmplitudeSpec::family`](../../src/bartcore/combiner.hpp)) at the `switch (spec.family)` arm of its constructor
(["switch (spec.family)"](../../src/bartcore/chain.hpp)), which makes the bcf row family-dependent ([f11]).

## Evidence

The per-model equivalence baseline, SBC verdict and tinytest inventory are
[9. Per-model evidence](../plans/review-2026-08-24/gate-ledger.md#9-per-model-evidence). Three canonical baselines:
`equivalence-fbff1989.rds` (50 scenarios), `bcf-equivalence-fbff1989.rds` (12 scenarios) and
`multinomial-equivalence-fbff1989.rds` (11 scenarios), all in benchmarks/baselines/MANIFEST.

## Gaps

Candidate work items grouped by what would need to change, not by which model asks; scheduling
is VD's. REFUSED (`R`) cells are absent, being part of the models.

| work item | unblocks | pointer |
|---|---|---|
| `xbart()` family coverage ([`xbart`](../../R/xbart.R) admits only auto/gaussian/probit/logistic) | student, nbinom, aft, hazard, hurdle, bcf, hetero | ordinal/multinom redirect to `bart2()` instead |
| Flat C reach for the K-forest softmax family | multinomial | [f4] |
| Warm start / grow-from-root for the alternate-family `bart2` arcs | ordinal, nbinom, multinomial, hurdle | [`checkFamilyUnsupportedArgs`](../../R/bart.R), [f12] |
| Multi-forest donor warm start | bcf (multinomial hits the same guard independently) | [`refuseMultiForestWarmStart`](../../src/R_interface_bartcore.cpp), [f12] |
| Real-valued (continuous) dispersion | nbinom | TODO `negbin-real-dispersion` |
| SBC at full chain length (r/agg.psi ridge) | nbinom | docs/plans/sbc-family-tiers.md |
| SBC gamma3 re-run at full chain length | ordinal | docs/plans/sbc-family-tiers.md |
| The heteroscedastic arms behind the variance-forest prior draw | hetero, heteroscedastic aft | docs/design/aft-status-setter.md, slices 3-4 |
| Register the exact oracle in the baseline MANIFEST | aft | benchmarks/R/aft-exact.R |
| An engine per-observation log-likelihood channel | multinomial | [`multinomialLogLik`](../../R/generics.R) |
| Whole-data `setData` | bcf, multinomial | docs/design/model-space-survey.md, Doors 1 and 3 |
| Equivalence scenarios and active-rows-mask evidence for the latent sub-families | bcf | docs/plans/bcf-latent-evidence.md; the exact gate is recorded and SBC is measured, both latent arms a chain-length finding rather than a matrix member |
| SBC coverage, deferred not blocked; liftable via `setState` | hetero | docs/plans/sbc-family-tiers.md |

**Not gaps** - structurally impossible or settled by decision, not open work:

- Hazard's and hurdle's absence from the flat C API and `xbart()`: neither owns engine code to
  expose ([f5], [f10]).
- Their absence from the SBC matrix: both designs break exchangeability on `y0` and neither
  owns sampling code (docs/plans/sbc-family-tiers.md).
- BCF's missing `bcf()`/`bartBCF` verb: it ships in bartCause, the K-forest capability itself
  being reachable ([Public creation surface (2026-08-10 to 2026-08-11)](bcf.md#public-creation-surface-2026-08-10-to-2026-08-11), [f6]).
- Multinomial's per-forest active-rows mask: refused permanently on softmax log-sum-exp
  grounds, not unbuilt ([Per family](active-rows-mask.md#per-family)).
- Multinomial's `$getLatents()`: a decided decline - the augmentation is meaningless between
  sweeps (["reports nothing, by a DECIDED decline"](multinomial.md)).

## Footnotes

[f1] `bartBT()` carries 0.9-34's argument list and no `family` formal at all
([`bartBT`](../../R/bart.R)), so the by-name family refusals and the token tables behind them are
gone (retired: [`bartRedirectedFamilies`](../../R/bart.R), [`bartOwnClassFamilies`](../../R/bart.R)); its one
response refusal is a factor of three or more levels, whose message names both remedies
([`refuseLegacyFactorResponse`](../../R/bart.R)). `resid.dist` is the separate Student-t lever, at the
modern door alone. `"twopart"` is no longer an alias at either door: it is refused by name
([`refuseTwopartFamily`](../../R/tombstones.R)).

[f2] Student-t is no `family` token and not in `dbarts_sampler_create`'s admission list: a
finite `resid.df` on the model SEXP selects it ([`parseSamplerSpecification`](../../src/R_interface_bartcore.cpp), [`residualDf`](../../src/R_interface_bartcore.cpp),
gaussian-only, refused elsewhere by ["student residuals require a continuous"](../../R/spec.R)), and
the engine family stays `gaussian`; the header's [`DBARTS_FAMILY_STUDENT`](../../inst/include/dbarts/dbarts.h) serves the
augmentation entries alone.

[f3] Ordinal and nbinom each ship a `DBARTS_FAMILY_*` enumerator; heteroscedastic has none,
being a control-attribute decoration. The header's specification-attribute block
(retired: ["SPECIFICATION ATTRIBUTES"](../../inst/include/dbarts/dbarts.h)) documents all three selectors - `bartcore.n.categories`,
`bartcore.dispersion` ([`parseControl`](../../src/R_interface_bartcore.cpp)) and `bartcore.variance`
([`applyVarianceAttributes`](../../src/R_interface_bartcore.cpp)).

[f4] `dbarts(x, y, family = "multinomial")` (matrix interface only) takes a counts matrix or a
one-hot-expanded factor response ([`dbarts`](../../R/dbarts.R), [`multinomial`](../../R/dbarts.R),
[`resolveMultinomialCounts`](../../R/data.R)); there is no separate creation entry and no `dbarts.h`
one at all, retired: [`creationFamilyName`](../../src/C_interface.cpp) refusing the token.

[f5] The three `"hazard"` spellings are person-period ingestion sugar:
[`expandDiscreteTimeHazard`](../../R/dbarts.R) expands the design and remaps the token -
`"hazard"`/`"hazard.probit"` to `"probit"`, `"hazard.logistic"` to `"logistic"` - before any
model is built, adding no engine code. The row is therefore the probit row, or the logistic
one under that third spelling: case weights `S`, latents the Polya-Gamma omegas
([Discrete-time hazard (LANDED 2026-07-18, 4bcdccf)](survival.md#discrete-time-hazard-landed-2026-07-18-4bcdccf)).

[f6] `treatment` is not a `bart2()` formal; the K-forest amplitude capability comes from a
`forest()` formula term, rewritten into the same `forests =` channel `dbarts()`/`dbartsSpec()`
use ([`ingestFormulaTerms`](../../R/formulaTerms.R)).

[f7] `updateScale` re-derives the internal response transform. The latent families have
`fitScale() == 1` and `fitShift() == 0` by definition, so there is nothing to re-anchor and
the flag is ignored rather than refused.

[f8] Logistic weights are the counts its Polya-Gamma latents are built from, so a swap is a
model change: [`LogisticResponse::setWeights`](../../src/bartcore/model.hpp) redraws omega against the new counts, and
the creation-time positive-integer policy holds on every conduit.

[f9] Multinomial carries a counts response and no weight vector: creation refuses one
([`parseMultinomialData`](../../src/R_interface_bartcore.cpp)) and the shared counts guard, keyed on `samplerCarriesCounts`,
blocks every mutation conduit that would touch it ([`refuseCountsMutation`](../../R/bartcore.R); full
inventory [1. Problem statement and inventory](multinomial-mutation-arc.md#1-problem-statement-and-inventory)).

[f10] Hurdle has no sampler of its own: [`bart2Hurdle`](../../R/bart.R) composes two ordinary `bart2()`
fits - occupancy probit and lognormal positive part - glued at report time
([2. Decision (fork 1, the gating question) - COMPOSE IN R, do not build in the engine](hurdle.md#2-decision-fork-1-the-gating-question---compose-in-r-do-not-build-in-the-engine)).

[f11] Under a latent sub-family the amplitude combination is the index on the link's fixed
scale, so sigma is pinned and the transform the identity
([The model](multiplier-combiner.md#the-model)), while [`Chain::latents`](../../src/bartcore/chain.hpp) bare-delegates
to the sub-family's own model with no coupling gate. `prior.scale` is refused at creation and
mid-chain alike ([`ForestSpec::amplitudePriorScale`](../../src/bartcore/combiner.hpp); [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp)
returns false whenever a combiner is installed; R-side
[`refuseAmplitudeMutation`](../../R/bartcore.R)).

[f12] [`checkFamilyUnsupportedArgs`](../../R/bart.R) raises a bare
`does not support 'warm.start' or 'n.grow.sweeps'` for the ordinal/nbinom/multinomial/hurdle
arcs, with no model reason stated. Independently a multi-forest DONOR warm start is refused at
the forest count everywhere ([`refuseMultiForestWarmStart`](../../src/R_interface_bartcore.cpp)): the install takes a saved
slot's trees but the donor's LIVE amplitudes, and nothing tests the result above one forest
([The mutation-legality table](bart-as-a-component.md#the-mutation-legality-table),
[Mutation surface](bcf.md#mutation-surface)).

[f13] A variance forest is built over gaussian or aft ([`varianceForestIsRefused`](../../src/bartcore/facade.hpp)); its
row's channels otherwise follow the base family's own row, except sigma and `updateScale`,
which the variance forest pins or refuses for both regardless of family. `setData` stays
refused under aft as always, since a whole-data replacement may change n
(["fix the censoring structure at creation"](../../src/R_interface_bartcore.cpp) on that conduit); the status itself moves
through `$setResponse(y, status = )`.

