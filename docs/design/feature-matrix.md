# Response-model feature matrix

Status: living reference, updated in place whenever a cell changes.

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

Path aliases used in cites:

    RIB   src/R_interface_bartcore.cpp      CAPI  inst/include/dbarts/dbarts.h
    MOD   src/bartcore/model.hpp            CH    src/bartcore/chain.hpp
    FAC   src/bartcore/facade.hpp           COM   src/bartcore/combiner.hpp
    MOV   src/bartcore/moves.hpp            SAM   src/bartcore/sampler.hpp
    bart.R, dbarts.R, spec.R, xbart.R, data.R, generics.R,
    A_class.R, bartcore.R, formulaTerms.R    -> R/<name>
    C_interface.cpp -> src/C_interface.cpp (the flat C entry points)
    test-*.R -> inst/tinytest/test-*.R
    sampler.Rd -> man/dbartsSampler-class.Rd; every other *.Rd -> man/<name>

Cites are by symbol and are existence-checked by `tools/check-doc-freshness.R`; a cell's VALUE
is adjudicated separately from its cite.

## Rows

Twelve rows. Six are response models proper, reached through the engine's own `ResponseFamily`
enum ([[MOD#ResponseFamily]]: gaussian, probit, logistic, aft, ordinal, nbinom); the other six
are reached some other way, so they need rows rather than an enum read. Leaf models (constant,
monotone, linear, GP) are an orthogonal axis, not rows; where one gates a capability the cell
or rule says so.

| key | model |
|---|---|
| gaussian | Gaussian (`ResponseFamily::gaussian`, [[MOD#GaussianResponse]]) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, [[MOD#TResponse]]) |
| probit | Binary probit ([[MOD#ProbitResponse]]) |
| logistic | Binary logistic, weights = observation counts ([[MOD#LogisticResponse]]) |
| ordinal | Ordered categorical, cumulative probit ([[MOD#OrdinalResponse]]) |
| nbinom | Negative binomial, positive-integer dispersion ([[MOD#NBResponse]]) |
| multinom | Multinomial softmax, K forests ([[MOD#MultinomialResponse]] + combiner) |
| aft | AFT survival, log-normal ([[MOD#AFTResponse]]) |
| hazard | Discrete-time hazard (person-period sugar, [[dbarts.R#expandDiscreteTimeHazard]]) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, [[bart.R#bart2Hurdle]]) |
| bcf | K-forest amplitude family, bcf's two forests being its K = 2 instance ([[COM#AmplitudeForestCombiner]]) |
| hetero | Heteroscedastic variance forest ([[CH#buildVarianceForest]]) |

## 1. Structural signature

Five facts a bridge or engine predicate decides for every row, whichever R entry point reaches
it: a live case-weight channel ([[RIB#familyCarriesNoWeights]],
[[RIB#refuseBinaryWeightChange]]); a drawn rather than pinned sigma - pinned by a variance
forest owning the scale, or by the family's own definition ([[RIB#sigmaIsPinned]],
[[RIB#refusePinnedSigmaChange]]); a persisted per-observation augmentation vector, a
`latents()` override ([[RIB#bartcore_getLatents]]); a non-trivial `fitScale`/`fitShift` at
creation, which makes `updateScale = TRUE` a re-anchor rather than a no-op ([f7]) or a refusal
([[RIB#refuseVarianceForestScaleUpdate]]); and a combined out-of-sample fit defined at all
([[RIB#refuseUndefinedTestFits]]).

| model | case weights | sigma | latents | unit-scale transform | test fits defined |
|---|---|---|---|---|---|
| gaussian | S [[MOD#GaussianResponse::setWeights]] | S [[RIB#bartcore_setSigma]] | - [[RIB#bartcore_getLatents]] | S [[MOD#GaussianResponse::setOffset]] | S [[RIB#bartcore_setTestPredictor]] |
| student | S [[MOD#TResponse::setWeights]] | S [[RIB#bartcore_setSigma]] | S [[MOD#TResponse::latents]] | S [[MOD#TResponse::setOffset]] | S [[RIB#bartcore_setTestPredictor]] |
| probit | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[MOD#ProbitResponse::latents]] | - [f7] | S [[RIB#bartcore_setTestPredictor]] |
| logistic | S [[MOD#LogisticResponse::setWeights]] [f8] | R [[RIB#refusePinnedSigmaChange]] | S [[MOD#LogisticResponse::latents]] | - [f7] | S [[RIB#bartcore_setTestPredictor]] |
| ordinal | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[MOD#OrdinalResponse::latents]] | - [f7] | S [[RIB#bartcore_setTestPredictor]] |
| nbinom | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[MOD#NBResponse::latents]] | - [f7] | S [[RIB#bartcore_setTestPredictor]] |
| multinom | R [[RIB#parseMultinomialData]] [f9] | R [[bartcore.R#refuseCountsMutation]] [f9] | R [[docs/design/multinomial.md#"reports nothing, by a DECIDED decline"]] | R [[bartcore.R#refuseCountsMutation]] [f9] | S [[RIB#bartcore_setTestPredictor]] |
| aft | R [[RIB#refuseBinaryWeightChange]] | S [[RIB#bartcore_setSigma]] | S [[MOD#AFTResponse::latents]] | S [[MOD#AFTResponse::setOffset]] | S [[RIB#bartcore_setTestPredictor]] |
| hazard | as probit [f5] | as probit | as probit | as probit | as probit |
| hurdle | - [f10] | - | - | - | - |
| bcf | S, gaussian/logistic [[RIB#bartcore_setWeights]] [f11] | S, gaussian only [[RIB#bartcore_setSigma]] [f11] | S, probit/logistic only [[CH#Chain::latents]] [f11] | R [[bartcore.R#refuseAmplitudeMutation]] | R [[RIB#refuseUndefinedTestFits]] |
| hetero | S [[RIB#bartcore_setWeights]] | R [[RIB#refusePinnedSigmaChange]] | - [[RIB#bartcore_getLatents]] | R [[RIB#refuseVarianceForestScaleUpdate]] | S [[RIB#bartcore_setTestPredictor]] |

Mutation channels and row subsetting read off the table: `setWeights`/`setSigma`/`getLatents`
follow their columns; `updateScale = TRUE` follows the unit-scale column (a no-op where the
transform is fixed, refused where a re-anchor would break a coupling's or a variance forest's
calibration); test predictors, offsets and `predict()` follow the test-fits column;
`setResponse`/`setOffset`/`setPredictor` (+ per-observation) are open wherever a sampler
exists, except multinomial blocks the first two ([f9]) and hurdle has none. Three exceptions
follow from no column:

- Zero-weight subsetting follows the case-weights column except logistic (and bcf's logistic
  sub-case): its weights are positive-integer Polya-Gamma trial counts, so a zero is refused
  although the channel is open ([[RIB#enforceBinaryWeightPolicy]]).
- Whole-data `setData` (n free) is narrower than `setPredictor`: it needs a sampler that owns
  its predictors ([[RIB#refusePredictorMutation]] refuses a data-handle view and a CSC-built
  design), and is refused above one forest ([[RIB#refuseMultiForestMutation]]) and for aft
  ([[RIB#"fix the censoring structure at creation"]]).
- Hurdle's `-` cells are not "ask the two components": `bart2()` refuses `weights`, `subset`
  and `offset`/`offset.test` on that family at its own entrance
  ([[bart.R#"does not support 'weights'"]]).

Two capabilities are universal, not columns: the mid-chain active-rows mask
([[CH#Chain::setActiveRows]], every family except hurdle, multinomial's being GLOBAL,
[[docs/design/active-rows-mask.md#The contract]]), and `extract(type = "loglik")` (every row,
including hurdle's composed density and bcf's combined-fit score; multinomial's engine-side
channel stays undefined, see Gaps). Named calibration (`$getCalibration`/`$setCalibration`, a
per-forest `prior.scale`, [[docs/design/nameable-calibration.md#The surface]]) is open on
every single-forest sampler, hetero included with its variance forest not counted
([[CH#buildVarianceForest]]), and refused on both couplings, whose leaf scale comes from a
calibration map ([f11]).

## 2. Reach

`xbart()` and the flat C API, cited by the token or by the refusal/absence.

| model | `xbart()` | flat C `dbarts.h` |
|---|---|---|
| gaussian | S [[xbart.R#xbart, gaussian]] | S [[CAPI#DBARTS_FAMILY_GAUSSIAN]] |
| student | M [[xbart.R#xbart]] | S [[RIB#parseSamplerSpecification, residualDf]] [f2] |
| probit | S [[xbart.R#xbart, probit]] | S [[CAPI#DBARTS_FAMILY_PROBIT]] |
| logistic | S [[xbart.R#xbart, logistic]] | S [[RIB#resolveFamily, logistic]] |
| ordinal | R [[data.R#resolveClassificationFamily]] | S [[RIB#resolveFamily, ordinal]] [f3] |
| nbinom | M [[xbart.R#xbart]] | S [[RIB#resolveFamily, nbinom]] [f3] |
| multinom | R [[data.R#resolveClassificationFamily]] | M [f4] |
| aft | M [[xbart.R#xbart]] | S [[CAPI#DBARTS_FAMILY_AFT]] |
| hazard | M [[xbart.R#xbart]] | M [f5] |
| hurdle | M | M [f10] |
| bcf | M [[xbart.R#xbart]] | S [[CAPI#dbarts_sampler_create]] |
| hetero | M [[xbart.R#xbart]] | S [[RIB#applyVarianceAttributes]] [f3] |

Construction reach through `bart()`, `bart2()` and `dbarts()` + R5: gaussian, student, probit,
logistic and aft reach all three ([[bart.R#bart]], [[bart.R#bart2]], [[dbarts.R#dbarts]]).
Ordinal, nbinom, multinomial and hurdle are refused at `bart()` by name
([[bart.R#refuseBartOwnClassFamily]]), hazard likewise
([[bart.R#refuseBartRedirectedFamily]]), and bcf and hetero are no `family` tokens at all
(each refusal at [f1]). `bart2()` ships all six: ordinal ([[bart.R#bart2Ordinal]]), nbinom
([[bart.R#bart2Negbin]]), multinomial ([[bart.R#bart2Multinomial]]), hurdle ([f10]), hazard as
person-period sugar ([f5]), bcf and hetero through `forests =` / `variance =` ([f6]).
`dbarts()` + R5 ships every one but hurdle, which it refuses
([[dbarts.R#"is only available through bart2()"]]): ordinal and nbinom
([[dbarts.R#dbarts, ordinal]], [[dbarts.R#dbarts, nbinom]]), multinomial on the matrix
interface only ([f4]), and hazard, bcf and hetero by the same routes.

## 3. How a fit is built

[[spec.R#dbartsSpec]] resolves the seven single-forest tokens - auto, gaussian, probit,
logistic, aft, ordinal, nbinom - takes `family = "multinomial"` directly
([[spec.R#dbartsSpec, multinomial]]), reaches the K-forest amplitude family through
`forests =` ([[spec.R#dbartsSpec, forests]], each declared `forest(basis = ...)`) and a
variance forest through `variance =` ([[spec.R#dbartsSpec, variance]]); only hazard and hurdle
stay out of its reach. A `forests =` fit resolves gaussian, probit or logistic only, aft,
ordinal and nbinom being refused by name at the R layer, the bridge and the factory alike
([[spec.R#"a treatment forest does not support family"]],
[[RIB#refusedAmplitudeFamilyReason]], [[FAC#createAmplitudeSampler]]); `bart2()` reaches the
same machinery through a `forest()` formula term, under an identical gate ([f6]).

## 4. Composition rules

A variance forest requires `family = "gaussian"`
([[spec.R#"a variance forest requires family"]]) - heteroscedastic IS that capability, so its
row carries no variance-forest cell - and even under gaussian refuses Student-t residuals
([[spec.R#"does not support Student-t residuals"]],
[[docs/design/heteroscedastic.md#Student-t residuals refused]]) and monotone constraints
([[spec.R#"not supported with monotone constraints"]]). It never takes DART either
([[CH#buildVarianceForest]] leaves `useDart` at its default false).

DART is refused for both multi-forest couplings by name - bcf at
[[spec.R#"a DART tree prior"]], multinomial at [[bart.R#"'dart' or a DART 'tree.prior'"]] -
and each coupling's forest builder hard-sets `forest.useDart = false` whatever the route asked
for ([[CH#buildSpecifiedForest]], [[CH#buildMultinomialForest]]). Warm start and
grow-from-root are refused, as `M` rather than `R`, for the four alternate-family `bart2` arcs
and for the multi-forest donor warm start ([f12]); grow-from-root itself ships, covered at two
forests, and is gated by the LEAF model rather than the family - linear and GP leaves are
refused in [[dbarts.R#growFromRoot]] and a no-op in [[CH#growForestFromRoot]], so every family
reads "constant leaf".

## 5. Combiners and couplings

Two rows sit on a COMBINER, an object holding K forests plus the rule that combines their fits
into one location per observation ([[COM#ForestCombiner]],
[[docs/design/forest-combiner.md#The ForestCombiner]]): multinomial's is a softmax over K
category forests, bcf's the AMPLITUDE family, forest f's fit at row i scaled by
`dot(a_f, B_f(i, .))` ([[docs/design/multiplier-combiner.md#The amplitude layout]]). A
K-forest chain takes its response model from `AmplitudeSpec::family`
([[COM#AmplitudeSpec::family]]) at the `switch (spec.family)` arm of its constructor
([[CH#"switch (spec.family)"]]), which makes the bcf row family-dependent ([f11]).

## Evidence

The per-model equivalence baseline, SBC verdict and tinytest inventory are
[[docs/plans/review-2026-08-24/gate-ledger.md#Per-model evidence]]. Three canonical baselines:
`equivalence-1e5f80b2.rds` (50 scenarios), `bcf-equivalence-3c81d6df.rds` (12 scenarios) and
`multinomial-equivalence-4d9a3337.rds` (11 scenarios), all in benchmarks/baselines/MANIFEST.

## Gaps

Candidate work items grouped by what would need to change, not by which model asks; scheduling
is VD's. REFUSED (`R`) cells are absent, being part of the models.

| work item | unblocks | pointer |
|---|---|---|
| `xbart()` family coverage ([[xbart.R#xbart]] admits only auto/gaussian/probit/logistic) | student, nbinom, aft, hazard, hurdle, bcf, hetero | ordinal/multinom redirect to `bart2()` instead |
| Flat C reach for the K-forest softmax family | multinomial | [f4] |
| Warm start / grow-from-root for the alternate-family `bart2` arcs | ordinal, nbinom, multinomial, hurdle | [[bart.R#checkFamilyUnsupportedArgs]], [f12] |
| Multi-forest donor warm start | bcf (multinomial hits the same guard independently) | [[RIB#refuseMultiForestWarmStart]], [f12] |
| Real-valued (continuous) dispersion | nbinom | TODO `negbin-real-dispersion` |
| SBC at full chain length (r/agg.psi ridge) | nbinom | docs/plans/sbc-family-tiers.md |
| SBC gamma3 re-run at full chain length | ordinal | docs/plans/sbc-family-tiers.md |
| A censoring-status setter, the SBC-coverage enabler | aft | docs/plans/sbc-family-tiers.md |
| Register the exact oracle in the baseline MANIFEST | aft | benchmarks/R/aft-exact.R |
| An engine per-observation log-likelihood channel | multinomial | [[generics.R#multinomialLogLik]] |
| Whole-data `setData` | bcf, multinomial | docs/design/model-space-survey.md, Doors 1 and 3 |
| Equivalence, SBC and active-rows-mask evidence for the latent sub-families | bcf | gaussian sub-family only measured today |
| SBC coverage, deferred not blocked; liftable via `setState` | hetero | docs/plans/sbc-family-tiers.md |

**Not gaps** - structurally impossible or settled by decision, not open work:

- Hazard's and hurdle's absence from the flat C API and `xbart()`: neither owns engine code to
  expose ([f5], [f10]).
- Their absence from the SBC matrix: both designs break exchangeability on `y0` and neither
  owns sampling code (docs/plans/sbc-family-tiers.md).
- BCF's missing `bcf()`/`bartBCF` verb: it ships in bartCause, the K-forest capability itself
  being reachable ([[docs/design/bcf.md#Public creation surface]], [f6]).
- Multinomial's per-forest active-rows mask: refused permanently on softmax log-sum-exp
  grounds, not unbuilt ([[docs/design/active-rows-mask.md#Per family]]).
- Multinomial's `$getLatents()`: a decided decline - the augmentation is meaningless between
  sweeps ([[docs/design/multinomial.md#"reports nothing, by a DECIDED decline"]]).

## Footnotes

[f1] `bart()`'s `family` formal is the narrow `c("auto", "logistic", "aft")`
([[bart.R#bart]]); `resid.dist` is the separate Student-t lever. The ten other tokens of
[[bart.R#bartRedirectedFamilies]] are refused BY NAME ahead of `match.arg`: five own-class
ones - the four above plus the `"twopart"` alias - through [[bart.R#bartOwnClassFamilies]],
and `"gaussian"`, `"probit"` and the three `"hazard"` spellings with their own message.

[f2] Student-t is no `family` token and not in `dbarts_sampler_create`'s admission list: a
finite `resid.df` on the model SEXP selects it ([[RIB#parseSamplerSpecification, residualDf]],
gaussian-only, refused elsewhere by [[spec.R#"student residuals require a continuous"]]), and
the engine family stays `gaussian`; the header's [[CAPI#DBARTS_FAMILY_STUDENT]] serves the
augmentation entries alone.

[f3] Ordinal and nbinom each ship a `DBARTS_FAMILY_*` enumerator; heteroscedastic has none,
being a control-attribute decoration. The header's specification-attribute block
([[CAPI#"SPECIFICATION ATTRIBUTES"]]) documents all three selectors - `bartcore.n.categories`,
`bartcore.dispersion` ([[RIB#parseControl]]) and `bartcore.variance`
([[RIB#applyVarianceAttributes]]).

[f4] `dbarts(x, y, family = "multinomial")` (matrix interface only) takes a counts matrix or a
one-hot-expanded factor response ([[dbarts.R#dbarts, multinomial]],
[[data.R#resolveMultinomialCounts]]); there is no separate creation entry and no `dbarts.h`
one at all, [[C_interface.cpp#creationFamilyName]] refusing the token.

[f5] The three `"hazard"` spellings are person-period ingestion sugar:
[[dbarts.R#expandDiscreteTimeHazard]] expands the design and remaps the token -
`"hazard"`/`"hazard.probit"` to `"probit"`, `"hazard.logistic"` to `"logistic"` - before any
model is built, adding no engine code. The row is therefore the probit row, or the logistic
one under that third spelling: case weights `S`, latents the Polya-Gamma omegas
([[docs/design/survival.md#Discrete-time hazard]]).

[f6] `treatment` is not a `bart2()` formal; the K-forest amplitude capability comes from a
`forest()` formula term, rewritten into the same `forests =` channel `dbarts()`/`dbartsSpec()`
use ([[formulaTerms.R#ingestFormulaTerms]]).

[f7] `updateScale` re-derives the internal response transform. The latent families have
`fitScale() == 1` and `fitShift() == 0` by definition, so there is nothing to re-anchor and
the flag is ignored rather than refused.

[f8] Logistic weights are the counts its Polya-Gamma latents are built from, so a swap is a
model change: [[MOD#LogisticResponse::setWeights]] redraws omega against the new counts, and
the creation-time positive-integer policy holds on every conduit.

[f9] Multinomial carries a counts response and no weight vector: creation refuses one
([[RIB#parseMultinomialData]]) and the shared counts guard, keyed on `samplerCarriesCounts`,
blocks every mutation conduit that would touch it ([[bartcore.R#refuseCountsMutation]]; full
inventory [[docs/design/multinomial-mutation-arc.md#Problem statement]]).

[f10] Hurdle has no sampler of its own: [[bart.R#bart2Hurdle]] composes two ordinary `bart2()`
fits - occupancy probit and lognormal positive part - glued at report time
([[docs/design/hurdle.md#COMPOSE IN R, do not build in the engine]]).

[f11] Under a latent sub-family the amplitude combination is the index on the link's fixed
scale, so sigma is pinned and the transform the identity
([[docs/design/multiplier-combiner.md#The model]]), while [[CH#Chain::latents]] bare-delegates
to the sub-family's own model with no coupling gate. `prior.scale` is refused at creation and
mid-chain alike ([[COM#ForestSpec::amplitudePriorScale]]; [[CH#Chain::setForestPriorScale]]
returns false whenever a combiner is installed; R-side
[[bartcore.R#refuseAmplitudeMutation]]).

[f12] [[bart.R#checkFamilyUnsupportedArgs]] raises a bare
`does not support 'warm.start' or 'n.grow.sweeps'` for the ordinal/nbinom/multinomial/hurdle
arcs, with no model reason stated. Independently a multi-forest DONOR warm start is refused at
the forest count everywhere ([[RIB#refuseMultiForestWarmStart]]): the install takes a saved
slot's trees but the donor's LIVE amplitudes, and nothing tests the result above one forest
([[docs/design/bart-as-a-component.md#mutation-legality table]],
[[docs/design/bcf.md#Mutation surface]]).

