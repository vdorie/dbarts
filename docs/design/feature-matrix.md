# Response-model feature matrix

Status: living reference, updated in place whenever a cell changes.

What each shipped response model can and cannot do, one row per model and one
column per capability that bears on scheduling. Every SHIPPED and REFUSED cell
carries a cite verified against the live tree; a wrong cell misdirects
scheduling, so a cell that cannot be verified would be marked `?` rather
than guessed.

## Legend

| code | meaning |
|---|---|
| `S` | SHIPPED. Works today; the cite is the site that makes it work. |
| `R` | REFUSED on model or identification grounds - the refusal is part of the model, not a hole. The cite is the refusal site. |
| `P` | PLANNED. A named design arc covers it; the cell names which stage. No cell currently carries this code. |
| `M` | MISSING. Not built, no schedule. A cite, when given, is a guard that errors *because the thing is unbuilt* - a recorded open item, not a model refusal. |
| `-` | N/A. The concept does not apply to this row; the row footnote says why. |
| `?` | UNVERIFIED. Constructs today with no refusal site, but no test, doc or adjudication backs it. Do not schedule against the cell until it is settled; every `?` would be listed under "Gaps". No cell currently carries this code. |

Path aliases used in cites:

    RIB   src/R_interface_bartcore.cpp      CAPI  inst/include/dbarts/dbarts.h
    MOD   src/bartcore/model.hpp            CH    src/bartcore/chain.hpp
    FAC   src/bartcore/facade.hpp           COM   src/bartcore/combiner.hpp
    MOV   src/bartcore/moves.hpp            SAM   src/bartcore/sampler.hpp
    bart.R, dbarts.R, spec.R, rbart.R, xbart.R, data.R, generics.R,
    A_class.R, bartcore.R      -> R/<name>
    bartcoreHandle.R -> inst/common/bartcoreHandle.R (the unexported low-level
    bartcoreXxx test handles)
    sampler.Rd -> man/dbartsSampler-class.Rd     bart.Rd -> man/bart.Rd
    bart2.Rd -> man/bart2.Rd

Cites are by symbol and are existence-checked by `tools/check-doc-freshness.R`; a cell's VALUE is adjudicated separately from its cite.

## Rows

Thirteen rows. The first ten are response models proper; the last three are
couplings or decorations over a base response that a user selects the same way
and schedules against the same way, so they earn rows.

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
| grouped | Grouped random intercepts ([[MOD#GroupedResponse]]) |
| hetero | Heteroscedastic variance forest ([[CH#buildVarianceForest]]) |

The engine's `ResponseFamily` enum has only six tokens ([[MOD#ResponseFamily]]:
gaussian, probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf,
grouped and hetero are all reached some other way, which is exactly why they
need rows here rather than an enum read. The bcf row is family-dependent
below because a K-forest chain selects its response model from
`AmplitudeSpec::family` ([[COM#AmplitudeSpec::family]]) through the
`switch (spec.family)` arm of the K-forest chain constructor
([[CH#"switch (spec.family)"]]), and the K-forest `Sampler` constructor takes
`family_(spec.family)` ([[SAM#Sampler]]) like any other chain. Leaf
models (constant, monotone, linear, GP) are an orthogonal axis and are not rows;
where a leaf model gates a capability the cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `rbart_vi()` | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|---|
| gaussian | S [[bart.R#bart]] | S [[bart.R#bart2, gaussian]] | S [[dbarts.R#dbarts, gaussian]] | S [[rbart.R#rbart_vi, gaussian]] | S [[xbart.R#xbart, gaussian]] | S [[CAPI#DBARTS_FAMILY_GAUSSIAN]] |
| student | S [[bart.R#residDist]] | S [[bart.R#bart2, resid.dist]] | S [[dbarts.R#dbarts, resid.dist]] | M [[rbart.R#"family <- match.arg(family)"]] | M [[xbart.R#xbart]] | S [[RIB#parseSamplerSpecification, residualDf]] [f2] |
| probit | S [[bart.Rd#y.train]] | S [[bart.R#bart2, probit]] | S [[dbarts.R#dbarts, probit]] | S [[data.R#resolveClassificationFamily]] | S [[xbart.R#xbart, probit]] | S [[CAPI#DBARTS_FAMILY_PROBIT]] |
| logistic | S [[bart.R#bart, logistic]] [f1] | S [[bart.R#bart2, logistic]] | S [[dbarts.R#dbarts, logistic]] | R [[rbart.R#rbart_vi]] | S [[xbart.R#xbart, logistic]] | S [[RIB#resolveFamily, logistic]] |
| ordinal | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Ordinal]] | S [[dbarts.R#dbarts, ordinal]] | R [[data.R#resolveClassificationFamily]] | R [[data.R#resolveClassificationFamily]] | S [[RIB#resolveFamily, ordinal]] [f3] |
| nbinom | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Negbin]] | S [[dbarts.R#dbarts, nbinom]] | M [[rbart.R#rbart_vi]] | M [[xbart.R#xbart]] | S [[RIB#resolveFamily, nbinom]] [f3] |
| multinom | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Multinomial]] | S [[dbarts.R#dbarts, multinomial]] [f4] | R [[data.R#resolveClassificationFamily]] | R [[data.R#resolveClassificationFamily]] | M [f4] |
| aft | S [[bart.R#bart, aft]] [f1] [f5] | S [[bart.R#bart2, aft]] | S [[dbarts.R#dbarts, aft]] | S [[rbart.R#rbart_vi, aft]] | M [[xbart.R#xbart]] | S [[CAPI#DBARTS_FAMILY_AFT]] |
| hazard | R [[bart.R#refuseBartRedirectedFamily]] [f1] | S [[bart.R#bart2, hazard]] | S [[dbarts.R#dbarts, hazardTokens]] | M [[rbart.R#rbart_vi]] | M [[xbart.R#xbart]] | M [f6] |
| hurdle | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Hurdle]] | R [[dbarts.R#"is only available through bart2()"]] | M | M | M [f6] |
| bcf | - [f1] | S [[R/formulaTerms.R#ingestFormulaTerms]] [f7] | S [[dbarts.R#dbarts, forests]] | M [[rbart.R#"family <- match.arg(family)"]] | M | S [[CAPI#dbarts_sampler_create]] |
| grouped | - [f1] | M [f8] | M [f8] | S [[rbart.R#"bartcore.groups"]] | M | S [[RIB#applyGroupAttribute]] [f3] |
| hetero | - [f1] | S [[bart.R#bart2, variance]] | S [[dbarts.R#dbarts, variance]] | M [[rbart.R#"family <- match.arg(family)"]] | M | S [[RIB#applyVarianceAttributes]] [f3] |

Ten tokens outside `bart()`'s own `c("auto", "logistic", "aft")` hit a
named refusal ahead of `match.arg`
([[bart.R#refuseBartRedirectedFamily]]); hazard's `bart()` cell is `R`
(deliberately refused by name), not `-` (no mechanism reaches it). [f1] has
the full breakdown of which token gets which reason.

[[spec.R#dbartsSpec]] resolves the seven single-forest tokens - auto,
gaussian, probit, logistic, aft, ordinal, nbinom - plus BCF through its
`forests =` argument ([[spec.R#dbartsSpec, forests]], `forest(basis = ...)`
replacing the removed `treatment =`) and a
variance forest through `variance =` ([[spec.R#dbartsSpec, variance]]). Its
`family` formal also accepts `"multinomial"` directly
([[spec.R#dbartsSpec, multinomial]], with dedicated body logic in
[[spec.R#resolveSamplerSpec, unsupportedMultinomial]]); only hazard, hurdle and
grouped stay out of `dbartsSpec()`'s reach. A `forests =` fit resolves
**gaussian, probit or logistic**; aft, ordinal and nbinom are
refused there by name, each stating what it is missing
([[spec.R#"a treatment forest does not support family"]], which now also
carries an explicit `multinomial =` arm refusing a `forests =` declaration for
that family by name - "its forests are its categories... not an amplitude
coupling"), with the same three-family gate at the bridge
([[RIB#refusedAmplitudeFamilyReason]], called from both creation routes,
[[RIB#refuseUnsupportedAmplitudeComposition]] and [[RIB#createBCFHolder]]) and
at the factory ([[FAC#createAmplitudeSampler]]).

`bart2()`'s formula interface reaches the same `forests =` machinery through
a `forest()` term rather than a formal of its own; VD decided the term
route over a flat `forests =` formal on `bart2()`, which was considered and
declined (docs/plans/archive/bart2-argument-consolidation.md) - so this is
a THIRD construction route into the bcf row, alongside
`dbarts()`/`dbartsSpec()`'s `forests =` and bartCause's `bcf()`; the family gate
is identical (gaussian, probit or logistic only) since it resolves through the
same `forests =` machinery. See [f7].

## 2. Mutation channels on the R5 `dbartsSampler`

The channels that make dbarts a conditional model inside an outer sampler.
`updateScale` is broken out because it is refused independently of the setter
it rides on.

| model | `setResponse` | `setOffset` | `updateScale = TRUE` | `setPredictor` (+ per-obs) | `setWeights` | `setSigma` | test surface |
|---|---|---|---|---|---|---|---|
| gaussian | S [[MOD#GaussianResponse::setResponse]] | S [[MOD#GaussianResponse::setOffset]] | S [[MOD#GaussianResponse::setOffset]] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] | S [[MOD#GaussianResponse::setWeights]] | S [[RIB#bartcore_setSigma]] | S [[RIB#bartcore_setTestPredictor, bartcore_setTestOffset]] |
| student | S [[MOD#TResponse::setResponse]] | S [[MOD#TResponse::setOffset]] | S [[MOD#TResponse::setOffset]] | S [[RIB#bartcore_setPredictor]] | S [[MOD#TResponse::setWeights]] | S [[RIB#bartcore_setSigma]] | S [[RIB#bartcore_setTestPredictor]] |
| probit | S [[MOD#ProbitResponse::setResponse]] | S [[MOD#ProbitResponse::setOffset]] | - [f9] | S [[RIB#bartcore_setPredictor]] | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |
| logistic | S [[MOD#LogisticResponse::setResponse]] | S [[MOD#LogisticResponse::setOffset]] | - [f9] | S [[RIB#bartcore_setPredictor]] | S [[MOD#LogisticResponse::setWeights]] [f10] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |
| ordinal | S [[MOD#OrdinalResponse::setResponse]] | S [[MOD#OrdinalResponse::setOffset]] | - [f9] | S [[RIB#bartcore_setPredictor]] | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |
| nbinom | S [[MOD#NBResponse::setResponse]] | S [[MOD#NBResponse::setOffset]] | - [f9] | S [[RIB#bartcore_setPredictor]] | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |
| multinom | R [[bartcore.R#refuseCountsMutation]] [f11] | R [f11] | R [f11] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] [f11] | R [f11] | R [f11] | S [[RIB#bartcore_setTestPredictor]] [f11] |
| aft | S [[MOD#AFTResponse::setResponse]] | S [[MOD#AFTResponse::setOffset]] | S [[MOD#AFTResponse::setOffset]] | S [[RIB#bartcore_setPredictor]] | R [[RIB#refuseBinaryWeightChange]] | S [[RIB#bartcore_setSigma]] | S [[RIB#bartcore_setTestPredictor]] |
| hazard | S [[MOD#ProbitResponse::setResponse]] [f6] | S [[MOD#ProbitResponse::setOffset]] | - [f9] | S [[RIB#bartcore_setPredictor]] | R [[RIB#refuseBinaryWeightChange]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |
| hurdle | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] |
| bcf | S [[CH#Chain::setResponse]] [f48] | S [[CH#Chain::setOffset]] [f48] | R [[bartcore.R#refuseAmplitudeMutation]] [f48] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] | S [[RIB#bartcore_setWeights]] [f48] | S [[RIB#bartcore_setSigma]] [f48] | R [[RIB#refuseUndefinedTestFits]] [f49] |
| grouped | S [[MOD#GroupedResponse::setResponse]] [f13] | S [[MOD#GroupedResponse::setOffset]] | R [[RIB#refuseGroupedScaleUpdate]] [f13] | S [[RIB#bartcore_setPredictor]] | S [[MOD#GroupedResponse::setWeights]] [f14] | S [[RIB#bartcore_setSigma]] [f14] | S [[RIB#bartcore_setTestPredictor]] |
| hetero | S [[RIB#bartcore_setResponse]] | S [[RIB#bartcore_setOffset]] | R [[RIB#refuseVarianceForestScaleUpdate]] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] | S [[RIB#bartcore_setWeights]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |

`setData` (whole-data replacement, n free) is dense-store and single-forest
only ([[RIB#refusePredictorMutation, refuseMultiForestMutation]]) and is
refused for grouped ([[RIB#"grouped random effects fix the data at creation"]])
and aft ([[RIB#"fix the censoring structure at creation"]]);
BCF/multinomial whole-data `setData` stays undesigned by the model-space
survey's verdict (model-space-survey.md's open questions 1 and 3).

## 3. Row subsetting, latents, calibration

| model | zero-weight row subset | active-rows mask [f15] | `getLatents` | pointwise loglik | nameable calibration [f16] |
|---|---|---|---|---|---|
| gaussian | S [[sampler.Rd#weights]], [[MOD#GaussianResponse]] [f17] | S [[MOD#GaussianResponse::setActiveRows]] | - [[RIB#bartcore_getLatents]] [f18] | S [[generics.R#pointwiseLogLikelihood]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| student | S [[MOD#TResponse::refreshLatents]] [f17] | S [[MOD#TResponse::setActiveRows]] | S [[MOD#TResponse::latents]] | S [[generics.R#isStudent]] [f19] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| probit | R [[RIB#enforceBinaryWeightPolicy]] | S [[MOD#ProbitResponse::setActiveRows]] | S [[MOD#ProbitResponse::latents]] | S [[generics.R#pointwiseLogLikelihood]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| logistic | R [[RIB#enforceBinaryWeightPolicy]] [f20] | S [[MOD#LogisticResponse::setActiveRows]] | S [[MOD#LogisticResponse::latents]] | S [[generics.R#pointwiseLogLikelihood]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| ordinal | R [[RIB#parseSamplerSpecification]] | S [[MOD#OrdinalResponse::setActiveRows]] | S [[MOD#OrdinalResponse::latents]] | S [[generics.R#ordinalLogLik]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| nbinom | R [[RIB#parseSamplerSpecification]] | S [[MOD#NBResponse::setActiveRows]] | S [[MOD#NBResponse::latents]] | S [[generics.R#negbinLogLik]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| multinom | R [[RIB#parseMultinomialData]] | S [[COM#MultinomialForestCombiner::setActiveRows]] [f21] | M [[MOD#MultinomialResponse]] [f22] | S [[generics.R#multinomialLogLik]] | R [f23] |
| aft | R [[RIB#parseSamplerSpecification]] | S [[MOD#AFTResponse::setActiveRows]] | S [[MOD#AFTResponse::latents]] | S [[generics.R#pointwiseLogLikelihood]] | S [[dbarts.R#getCalibration, setCalibration]] [f16] |
| hazard | R [[RIB#enforceBinaryWeightPolicy]] [f6] | S [[MOD#ProbitResponse::setActiveRows]] [f6] | S [[MOD#ProbitResponse::latents]] | S [[generics.R#pointwiseLogLikelihood]] [f24] | S [[dbarts.R#getCalibration, setCalibration]] [f6] |
| hurdle | R [[bart.R#"does not support 'weights'"]] | - [f12] | - [f12] | S [[generics.R#hurdleLogLik]] [f25] | - [f12] |
| bcf | S [[COM#AmplitudeForestCombiner::formForestResponse]] [f17] [f48] | S [[MOD#GaussianResponse::setActiveRows]], [[CH#composeForestWeights]] [f26] | S [[CH#Chain::latents]] [f18] | M | R [f23] |
| grouped | S [[MOD#drawGroupEffects]] | S [[MOD#GroupedResponse::setActiveRows]] [f27] | S [[MOD#GroupedResponse::latents]] | S [[generics.R#pointwiseLogLikelihood]] | S [[MOD#GroupedResponse::fitScale, GroupedResponse::fitShift]] [f27] |
| hetero | S [[CH#sweepVarianceForest]], [[MOD#ConstantVarianceLeaf::accumulate]] | S [[CH#formMeanWeights]] [f27] | - [f18] | S [[generics.R#heteroscedasticScale]] [f28] | S [[test-calibration-prior-draws.R#"heteroscedastic = anchorSampler"]] [f29] |

## 4. Model composition

| model | variance forest | grouped ranef | DART | warm start | grow-from-root |
|---|---|---|---|---|---|
| gaussian | S [[FAC#createSampler]] | S [[CH#GroupedResponse]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| student | R [[spec.R#"does not support Student-t residuals"]] [f30] | S [[CH#GroupedResponse]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| probit | R [[spec.R#"a variance forest requires family"]] | S [[CH#GroupedResponse]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| logistic | R [[spec.R#"a variance forest requires family"]] | S [[CH#GroupedResponse]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| ordinal | R [[spec.R#"a variance forest requires family"]] | M [[RIB#"not supported for ordinal responses"]] [f31] | S [[CH#useDart]] | R [[bart.R#checkFamilyUnsupportedArgs]] | R [[bart.R#checkFamilyUnsupportedArgs]] |
| nbinom | R [[spec.R#"a variance forest requires family"]] | M [[RIB#"not supported for count (nbinom) responses"]] [f31] | S [[CH#useDart]] | R [[bart.R#checkFamilyUnsupportedArgs]] | R [[bart.R#checkFamilyUnsupportedArgs]] |
| multinom | R [[bart.R#unsupported]] | M [[RIB#applyGroupAttribute]] [f32] | R [[bart.R#"'dart' or a DART 'tree.prior'"]] [f33] | R [[bart.R#checkFamilyUnsupportedArgs]] | R [[bart.R#checkFamilyUnsupportedArgs]] |
| aft | R [[spec.R#"a variance forest requires family"]] | S [[CH#GroupedResponse]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| hazard | R [[spec.R#"a variance forest requires family"]] | M [[rbart.R#rbart_vi]] [f6] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| hurdle | R [[spec.R#"a variance forest requires family"]] [f34] | M | S [[bart.R#redirectCall]] [f35] | R [[bart.R#checkFamilyUnsupportedArgs]] | R [[bart.R#checkFamilyUnsupportedArgs]] |
| bcf | R [[FAC#createAmplitudeSampler]] [f48] | R [[RIB#refuseUnsupportedAmplitudeComposition]] | R [[spec.R#"a DART tree prior"]] | S [[SAM#Sampler::installForests]] [f36] | S [[CH#growForestFromRoot]] [f36] |
| grouped | R [[spec.R#"does not support grouped random effects"]] [f30] | - | S [[rbart.R#usesDart]] | M [[rbart.R#rbart_vi]] [f37] | M [[rbart.R#rbart_vi]] [f37] |
| hetero | - | R [[spec.R#"does not support grouped random effects"]] [f30] | S [[CH#useDart]] [f38] | S [[SAM#Sampler::installForests]] | S [[CH#growForestFromRoot]] |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused in [[dbarts.R#growFromRoot]] and a no-op in
[[CH#growForestFromRoot]], so every family above
reads "constant leaf" in that column.

## 5. Evidence

| model | equivalence baseline [f39] | SBC verdict [f40] | dedicated tinytest files |
|---|---|---|---|
| gaussian | 33 scenarios (friedman, weighted, splitprobs, chains, setdata, wtoffset, quants, categorical, missing, dart, linear, gp, zeroweights, sparse, wtgp, chik2, predswap, predcol, predpartial, predreject, predforce, ordfactor, nafactor, sparsefactor, testswap, leaffactor, leaffactormixed, factorpartial, xbartmixed, xbart, bart2gauss, bart2twoforest, mixedmatrix) | PASS 7/7 | ~20 (test-sampler-*.R, test-bart-bart2.R, test-zero-weights.R) |
| student | `student` | PASS 4/4 | test-robust-errors.R only |
| probit | `probit`, `chik` | PASS | test-binaryResponse-hyperprior.R, test-family.R, test-weighted-binary-ppd.R |
| logistic | `logistic`, `wtlogistic` | PASS 6/6 | test-weighted-logistic.R, test-family.R |
| ordinal | `ordinal` | 9/10 [f41] | test-ordinal.R only |
| nbinom | `nbinom` | 1/3 [f42] | test-nbinom.R, test-dispersion-channel.R |
| multinom | 11 scenarios, own harness | aggregate PASS, raw `f_ik` PASS [f43] | 6 (test-multinomial-*.R) |
| aft | `grouped_aft` only [f44] | OUT [f45] | test-aft.R, test-rbart-aft.R |
| hazard | `hazard` | OUT [f45] | test-hazard.R only |
| hurdle | `hurdle` | OUT [f45] | test-hurdle.R, test-hurdle-surface.R |
| bcf | 12 scenarios, own harness, gaussian only [f48] | PASS, gaussian only [f46] | 8 (test-bcf*.R) |
| grouped | `grouped`, `grouped_aft` | PASS (tier A) | 14 (test-rbart-*.R) |
| hetero | `hetforce`, `hetswap`, `hetpartial` | OUT [f47] | 4 (test-heteroscedastic*.R) |

The rows are keyed by response model; predictor SHAPE cuts across them.
Eight gaussian scenarios carry the factor shapes: `ordfactor` (an ordered
factor at K = 150, three interior levels absent from training),
`nafactor` (two NA-bearing factor columns, the MIA-on-a-factor-column
anchor), `sparsefactor` (a CSC-backed `sparseFactor`, its count from the
container's declared level table), and the five store-path scenarios
`testswap`, `leaffactor`, `leaffactormixed`, `factorpartial` and
`xbartmixed`. Four of them - `ordfactor`, `leaffactor`, `leaffactormixed`
and `xbartmixed` - carry an ordered-factor predictor, so the ordered-factor
cut grid rests on those plus `benchmarks/R/categorical-exact.R`'s
ordered-factor case.

`inst/tinytest/test-capi.R` drives `""`/`"probit"` (`ptr1`/`ptrBinary`),
logistic, ordinal, aft and nbinom (each through `dbarts_sampler_create` with
its family token, run and checked for finite, correctly-shaped output) - the
whole single-forest family list - plus grouped
([[test-capi.R#"bartcore.groups"]]) and heteroscedastic
([[test-capi.R#"bartcore.variance"]]) by control attribute, and BCF
([[test-capi.R#"zBCF"]]) through `forests = list(forest(basis = ...))`.
Multinomial has no flat-C creation path to test ([f4]).

## Footnotes

[f1] `bart()` carries a narrow, appended `family` formal, `c("auto",
"logistic", "aft")` ([[bart.R#bart]]'s formals), forwarded
to `dbarts()` verbatim; `"auto"` (the default) still infers gaussian/probit
from the response, and `resid.dist` remains the separate Student-t lever,
untouched by `family`. All ten tokens outside that vocabulary are refused BY
NAME ahead of `match.arg`'s generic message, by the shared
[[bart.R#refuseBartRedirectedFamily, bartRedirectedFamilies]], called from
`bart()`: five own-class tokens - `"multinomial"`, `"ordinal"`, `"nbinom"`,
`"hurdle.lognormal"` and `"twopart"` - redirect through
[[bart.R#refuseBartOwnClassFamily, bartOwnClassFamilies]] (ordinal's own
auto-detection off an ordered-factor response routes through the same
helper). The remaining five - `"gaussian"`, `"probit"`, `"hazard.probit"`,
`"hazard"`, `"hazard.logistic"` - get their own named message instead of
`match.arg`'s generic one: the first three because `"auto"` already gives no
extra capability, the last two because the discrete-time expansion needs
`breaks`/`max.rows`, which `bart()` has no formal for (both reasons stated in
man/bart.Rd's `family` item). `hazard`'s section-1 `bart()` cell accordingly
reads `R` (deliberately refused by name) rather than `-` (no mechanism
reaches it). `bcf`, `grouped` and `hetero` are not `family` tokens at all and
stay out of reach by signature.

[f2] Student-t is not a token anywhere. It is selected by a finite `resid.df`
attribute on the model SEXP ([[RIB#parseSamplerSpecification, residualDf]],
gaussian-only gate) and refused for every non-gaussian family R-side by
[[spec.R#"student residuals require a continuous gaussian response"]]. The
engine family stays `gaussian`, which is why the whole gaussian row applies to
it in tables 2-4.

[f3] Ordinal and nbinom are reachable through `dbarts_sampler_create` and
ship a `DBARTS_FAMILY_*` enumerator each; grouped and heteroscedastic are
decorations with no enumerator of their own, selected instead by a
control attribute. None of the four control attributes is documented in
the shipped header: ordinal needs `bartcore.n.categories`
([[RIB#parseControl, bartcore.n.categories]]), nbinom `bartcore.dispersion`
([[RIB#parseControl, bartcore.dispersion]]), grouped `bartcore.groups`
([[RIB#applyGroupAttribute]]), heteroscedastic `bartcore.variance`
([[RIB#applyVarianceAttributes]]). The K-forest paragraph beside
[[CAPI#dbarts_sampler_create]]'s own comment block names gaussian, probit
and logistic, having replaced its "Gaussian responses only"
([[CAPI#"on a single forest or on K of them"]]).

[f4] `dbarts(x, y, family = "multinomial")` (matrix interface only) accepts a
counts matrix or a factor/character/integer-code response, one-hot expanded
([[dbarts.R#dbarts, multinomial]]); `resolveMultinomialCounts` builds the
counts matrix - now defined in [[data.R#resolveMultinomialCounts]], called from
[[dbarts.R#resolveMultinomialCounts]]. Creation routes through the same public
dispatch every family uses, [[RIB#bartcore_create]], whose multinomial arm is
[[RIB#createMultinomialDataHolder]]; there is no separate multinomial
creation entry point. Still has no dbarts.h creation path.

[f5] `family = "aft"` is now an explicit, appended token on `bart()`
(`c("auto", "logistic", "aft")`), documented
in man/bart.Rd's `family` item, rather than an undocumented `Surv()`
auto-dispatch quirk - see [f1]. The
underlying `Surv()`/two-column-`y.train` detection is a separate mechanism and
resolves the response shape on its own.

[f6] `family = "hazard"` / `"hazard.probit"` / `"hazard.logistic"` is
person-period ingestion sugar: [[dbarts.R#expandDiscreteTimeHazard]] expands the
design and [[dbarts.R#"the remap: the engine-facing family is now an ordinary binary link"]] remaps the token to `"probit"` or
`"logistic"` before any model is built. The resulting sampler *is* an ordinary
binary one, so its whole row equals the probit (or logistic) row, and the fit
records `family = "probit"`. No engine code, hence no C-API token and no SBC
arm. Refusals inside the expander: no formula interface
([[dbarts.R#"discrete-time hazard fits currently use the matrix interface"]]),
no `subset` ([[dbarts.R#"survival responses do not support 'subset'"]]), no
`test` ([[dbarts.R#"discrete-time hazard fits do not take a 'test' set"]]).

[f7] `treatment` is still not a `bart2()` formal: its full formal list
([[bart.R#bart2]]) carries no `...` at all, so an unrecognized argument like
`treatment =` hits ordinary R argument-matching ("unused argument"), not a
dbarts-authored construct. The
general K-forest amplitude capability this row names IS reachable from
`bart2()`'s FORMULA interface: a `forest()` term - `z:forest(x1 + x2)`, or
the general `forest(x1 +
x2, basis = ~z)` form - rewrites the formula and feeds the same `forests =
list(forest(), forest(basis = ))` channel `dbarts()`/`dbartsSpec()` already used
([[R/formulaTerms.R#ingestFormulaTerms]], called from `dbarts()`'s formula path;
`bart2` dispatches through it for every family except multinomial, which refuses
a `forest()` term by name since it has no amplitude-coupled slot). A two-forest
fit built this way carries the same `bartcore.forests` control attribute as one
built through `forests =` ([[test-formula-terms.R#"Block B"]], byte-identical at
the same seed). What remains absent is the NAMED causal verb, not the general
capability: `bcf()` and `bartBCF` still ship in **bartCause**, not dbarts, a
decision VD resolved on 2026-08-11
([[docs/plans/archive/bcf-public-surface.md#S5. bcf() and the fit class]];
echoed [[docs/design/bcf.md#Public creation surface]]). `forests =` on
`bart2()` itself was considered and declined as a second spelling - the
formula route is the only one bart2 gets; forest 1's `sd`/`update.amplitude`
have no term spelling either, which is why `bart2()`'s own knobs still don't
reach a term's forest the way `forests =` on `dbarts()` would
(docs/plans/archive/bart2-argument-consolidation.md, section 5).

[f8] `bartcore.groups` is written at exactly one site in `R/`,
[[rbart.R#"bartcore.groups"]] (two tinytest files also set the attribute
directly for engine-level tests), and no other
entry point carries a `group.by` formal, so grouped random effects are an
`rbart_vi()`-only surface.

[f9] `updateScale` re-derives the internal response transform. The latent
families have `fitScale() == 1` and `fitShift() == 0` by definition, so there
is no transform to re-anchor and the flag is ignored rather than refused.

[f10] Logistic weights are the observation counts its Polya-Gamma latents are
built from, so a swap is a model change rather than a reweighting:
[[MOD#LogisticResponse::setWeights]] redraws omega against the new counts
before returning, which is what makes the
conduit coherent and is why it was the one weight refusal recorded here as
"unbuilt" rather than "incoherent". The positive-integer policy creation states
([[RIB#enforceBinaryWeightPolicy]]) holds on every mutation conduit too, and
`setData` hands the replacement counts through the same conduit, so a data swap
draws rather than cold-starts, and replacement data given without weights is
single-trial. Probit, ordinal, aft and nbinom stay refused by identification.
The saved state carries a `weights.digest` attribute (a byte hash of the
weights in force, top level, additive - no version bump) and `setState`
re-derives the latents through this same conduit when the destination's
weights differ from it, so a restore lands where a swap lands; on a match it
re-derives nothing and the round trip stays byte-identical.

[f11] `bart2(family = "multinomial")` builds its `dbartsSampler` directly
(multinomial-mutation-arc.md): `$fit` is the K-forest engine that ran, one
`bartcore_create`, no host shell and no `$bc`. `setResponse`, `setOffset`,
`setWeights`, `setSigma`, `setCalibration`, `setForestWeights`,
`setForestBasis`, `getFitsWithoutOffset`, `setModel` and `setData` - ten
methods in all - are refused by name, naming the capability and the channel
that serves the caller where one exists
([[bartcore.R#refuseCountsMutation]], the shared multinomial guard, called
from all ten). `setCounts`, `setCategoryOffset` and `setCategoryTestOffset`
are the three response channels, public R5 methods
([[dbarts.R#setCounts, setCategoryOffset, setCategoryTestOffset]]).
`setPredictor`, `setCutPoints`, `setTestPredictor` and the global
`setActiveRows` stay open, which is why section 2's `setPredictor` and test-
surface cells read `S`.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction
([[dbarts.R#"is only available through bart2()"]]) and [[bart.R#bart2Hurdle]]
composes two ordinary `bart2()` fits - an occupancy probit and a lognormal
positive part, both built through [[bart.R#redirectCall]] - glued at report
time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f13] A grouped sampler accepts a
same-length setResponse and setOffset at the fixed scale - faithful
delegation, [[MOD#GroupedResponse::setResponse]] - and
[[RIB#refuseGroupedScaleUpdate]] refuses
updateScale != FALSE only under a base family with a data-derived transform
(gaussian, which is Student-t's report, and aft): b and tau are held on the
base's internal scale and converted by nothing, so a re-anchoring swap would
silently restate both in response units. Grouped probit and logistic take
updateScale = TRUE as the no-op it always was. The flat C API guards through
the same call ([[C_interface.cpp#refuseGroupedScaleUpdate]]). setData stays
refused.

[f14] Reads off the BASE family: grouped gaussian takes `setWeights`
([[MOD#GroupedResponse::setWeights]]) and `setSigma`; grouped probit is refused
on both ([[RIB#refuseBinaryWeightChange]], [[RIB#refusePinnedSigmaChange]]);
grouped aft takes `setSigma` and refuses `setWeights`.

[f15] `setActiveRows` is a first-class 0/1 per-observation mask, ARC COMPLETE
(docs/plans/latent-subset-mask.md), that each family composes into its own
precision vector, with the latent draw skipped for an inactive row. The
engine channel - [[CH#Chain::setActiveRows]] (the single validating and
normalizing scan), [[SAM#Sampler::setActiveRows]], the facade's pure virtual
[[FAC#SamplerBase::setActiveRows]] and its shape probe
[[FAC#SamplerShape::supportsActiveRows]] - reaches every family through the
R5 `$setActiveRows` ([[dbarts.R#setActiveRows]]) and the bridge entry
([[RIB#bartcore_setActiveRows]]). Logistic
([[MOD#LogisticResponse::workingWeights]]) and nbinom
([[MOD#NBResponse::workingWeights]]) serve a SEPARATE a_i omega_i composite
rather than writing the zero into omega_ itself, since the working response
divides by it and 0 * inf in the node kernels is a NaN; nbinom's
[[MOD#NBResponse::setActiveRows]] additionally restricts the collapsed
statistic S the dispersion grid draw reads and rebuilds the count-histogram
kernel behind L_k ([[MOD#NBDispersionPrior::computeKernel]]) over the active
rows at every mask change. aft's [[MOD#AFTResponse::setActiveRows]] composes
into its contained Gaussian, inheriting the sigma degrees-of-freedom
recount, and skips the censored redraw at an inactive row. All three report
NaN pointwise log-likelihood at an inactive row.

Multinomial's mask is GLOBAL, landing on the softmax coupling rather than
the response, which holds no precisions of its own:
[[MOD#MultinomialResponse::setActiveRows]] is a pass-through that only
advertises the capability ([[MOD#MultinomialResponse::supportsActiveRows]]),
and `Chain::setActiveRows` forwards the mask to
[[COM#MultinomialForestCombiner::setActiveRows]] after the response's own
install ([[COM#ForestCombiner::setActiveRows]] is the inert default every
additive coupling relies on instead). An inactive row's K interleaved
Polya-Gamma draws are SKIPPED, not drawn and discarded, in
[[COM#drawForestGlue]], and its composed precision is zeroed in every
category in [[COM#MultinomialForestCombiner::formForestResponse]]; the row
keeps its leaf occupancy and its reported softmax probabilities, and omega
is never zeroed since the working response divides by it. PER-FOREST
masking is refused permanently on model grounds at the only reachable
per-forest, per-observation channel, [[RIB#bartcore_setForestWeights]] - see
[f21]. The bridge's active-row refusal
([[RIB#"active-row masking is not implemented for this response family"]])
no longer names multinomial by hand: the message is family-generic, reached
only by a future family that does not override the base refusal.

Oracles: per-family kernel comparisons against the compacted case, bitwise
in value and in RNG stream
([[tests/cpp/test_model.cpp#testActiveRowsLogisticKernel]],
[[tests/cpp/test_model.cpp#testActiveRowsNBKernels]],
[[tests/cpp/test_model.cpp#testActiveRowsAFTCensored]],
[[tests/cpp/test_sampler.cpp#testActiveRowsMultinomialKernel]] - each latent
being a rejection sampler means a discard-rather-than-skip at an inactive
row fails the check outright, and the shape probe flips in
[[tests/cpp/test_shape.cpp#testMultinomial]]); a sampler-level conditional
independence check under substituted inactive responses for every affected
family ([[test-active-rows-pins.R#"logistic, nbinom and aft"]],
[[test-active-rows-pins.R#"multinomial, GLOBAL only"]]: substituting
arbitrary in-support values at the inactive rows leaves every active row's
recorded draw bitwise, including both successes and trial counts for
multinomial's PG(n_i, .) draws); heteroscedastic and grouped tested bitwise
against `setWeights(w * a)` ([[test-active-rows-pins.R#"heteroSampler"]];
[[tests/cpp/test_sampler.cpp#testActiveRows]]); and, for multinomial, an
all-zeros mask run (every category forest at its prior, every row still
reporting a simplex) plus the `setForestWeights` model-grounds refusal
(test-forest-weights.R). The flat-C entry, `dbarts_sampler_setActiveRows`
([[C_interface.cpp#dbarts_sampler_setActiveRows]]), probes
`shape.supportsActiveRows` first and never switches on family, so a probit
sampler is reachable from C too, since every `ResponseModel` subclass
reports `supportsActiveRows`
([[test-capi.R#"capi_set_active_rows"]] pins a genuine mask, an all-ones
no-op, a fractional refusal, a NULL clear, and a probit mask moving draws).
ARC COMPLETE, R and flat C alike.

[f16] `prior.scale` names the per-forest prior ANCHOR (the forest-total
prior scale at k = 1, in response units) rather than an sd, with a
`$getCalibration`/`$setCalibration` pair that reads and writes it on every
chain of a single-forest sampler (1-based `forest` arg,
[[bartcore.R#resolveForestIndex]]). Refused under a `k` hyperprior (a
sampled `k` has no single value to divide by, and the same holds once the
chains' `k` have diverged) and for BCF/multinomial forests, both at
creation and again mid-chain (see [f23]); `prior.mean` is refused as not
writable, naming the `setOffset` recipe instead. Engine:
[[CH#Chain::resolvedNodeScale]] (creation-time conversion),
[[CH#Chain::forestCalibration, Chain::setForestPriorScale]] (the mid-chain
reader/writer, sharing one conversion so neither can drift from the
other); R5 [[dbarts.R#getCalibration, setCalibration]]; bridge
[[RIB#bartcore_getCalibration, bartcore_setCalibration]]; flat C
`dbarts_sampler_getForestCalibration`/`setForestPriorScale`
(`inst/include/dbarts/dbarts.h`). Tests: test-calibration-creation.R,
test-calibration-prior-draws.R, test-calibration-midchain.R,
test-bcf-family.R, test-forest-basis-r5.R, test-capi.R, and
`tests/cpp/test_sampler.cpp#testForestCalibration`.

[f17] Zero weights are accepted, not refused
([[A_class.R#"'weights' must all be non-negative"]] errors only below zero and
warns that zeros are ignored; bridge [[RIB#enforceBinaryWeightPolicy]]). The
conditionals are exact - leaf suffstats multiply by `w`
([[MOD#ConstantVarianceLeaf::accumulate]],
[[MOD#LinearGaussianLeaf::accumulateNodeStatistics]]), and the sigma posterior
counts only positive-weight rows ([[MOD#numPositiveWeights_]], recounted on
every install in [[MOD#GaussianResponse::installWeights]], consumed in
[[MOD#GaussianResponse::drawSigma]]). The one named
inexactness against a true subset fit is CLOSED (`empty-leaf-veto-fix`,
2026-08-12): the empty-leaf veto counts
POSITIVE-WEIGHT members, so a leaf held alive only by zeroed rows is empty and
its branch is vetoed, on the conjugate path (MOV) and the constrained-leaf path
(MOD) alike. Occupancy elsewhere - the birth scan's `count`,
`collapseEmptyNodes`' trigger, `stateIsValid` - still counts members
deliberately, so this does NOT make zero-weight occupancy match a compacted fit;
see [[docs/design/empty-leaf-veto.md#What counts as empty]]. The same fix covers
the Student-t row and a GAUSSIAN K-forest. It says nothing about a probit or
logistic one, where a zero weight cannot exist to begin with: probit refuses
weights entirely and logistic holds them to positive integer counts, both at
creation, so the cell is family-dependent ([f48]).

[f18] For gaussian and heteroscedastic no latent vector exists: both leave
[[MOD#ResponseModel::latents]] at its nullptr default, and the bridge returns
`R_NilValue`. A K-forest sampler is no longer one of them.
[[CH#Chain::latents]] is a bare delegation to `response_->latents()` carrying
no coupling gate and no family switch, and [[RIB#bartcore_getLatents]] gates
only on that pointer being null, so a probit K-forest reports its
truncated normals ([[MOD#ProbitResponse::latents]]) and a logistic one its
Polya-Gamma omegas ([[MOD#LogisticResponse::latents]]). A
GAUSSIAN K-forest still reports none, which is why this cell is
family-dependent rather than plainly `S`.

[f19] A Student-t fit records `family = "gaussian"`, so `extract(type =
"loglik")` takes the gaussian branch of
[[generics.R#pointwiseLogLikelihood]], which now distinguishes the two: an
[[generics.R#isStudent]] check on `resid.dist` scores the MARGINAL t density,
`dt((y - ev) / sd, df, log = TRUE) - log(sd)` ([[generics.R#dt]]), rather than
folding a Student-t fit into the gaussian `dnorm` call ([[generics.R#dnorm]])
the way it once did. Pinned by value:
[[test-pointwise-loglik.R#"ll.t"]] is checked against `dt(...)` directly at
three indices, so the cell is `S`.

[f20] Weights on logistic are PG copy counts and a zero count is refused by name
at creation ("drop zero-count rows", [[RIB#enforceBinaryWeightPolicy]]; R mirror
[[spec.R#enforceWeightPolicy]]), so zero-weight subsetting is foreclosed for
this family by the weight semantics themselves - it is exactly the hole the
mid-chain `setActiveRows` channel later filled, rather than by any change to
the zero-count creation refusal, which stands.

[f21] PER-FOREST masking stays REFUSED, permanently and on model
grounds: the softmax margin is a log-sum-exp over the other K-1 forests,
so a row absent from category k's forest is still in every other
category's likelihood, and "row i is out of category k only" restricts no
likelihood at all. The refusal lands at the only reachable per-forest,
per-observation channel, [[RIB#bartcore_setForestWeights]], naming the
model reason rather than "unbuilt". BCF's per-forest weight acceptance at
that same channel stands unaffected - a different (additive) coupling
where the per-forest mask is redundant with, not incoherent under, the
combined likelihood (see [f26]).

[f22] Multinomial's omegas live in the combiner, not the response model, and
[[MOD#MultinomialResponse]] does not override `latents()`, so
`getLatents` returns NULL. No accessor exposes them.

[f23] A named `prior.scale` is refused for BCF and multinomial forests both AT
CREATION and MID-CHAIN, by design - their per-forest leaf scales come from a
calibration map that owns them ([[COM#ForestSpec::amplitudePriorScale]]), so a
named value has nowhere to land. Three creation-time refusal sites: R-side
`dbartsSpec()`'s BCF composition
([[spec.R#"a named 'prior.scale'"]], an entry of the `unsupported` vector,
beside the non-default-`k` entry [[spec.R#"a non-default 'k'"]]), the engine's
own BCF-composition gate ([[RIB#refuseUnsupportedAmplitudeComposition]]), and
the multinomial forest builder ([[RIB#buildMultinomialSampler]]). A later
change adds the mid-chain refusals too, at TWO independent sites rather
than one shared gate: `$setCalibration`'s R5 method refuses BCF
through [[dbarts.R#refuseAmplitudeMutation]] (measured at
[[test-calibration-midchain.R#"multi-forest calibration map"]]) before ever
reaching the bridge; a multinomial fit's `$fit` refuses the same call the same
way, for the same softmax-calibration-map reason - not through a now-deleted
host-shell mechanism (multinomial-mutation-arc.md); underneath both, the
engine-level gate any DIRECT low-level call still hits -
[[CH#Chain::setForestPriorScale]] returning false whenever `combiner_ !=
nullptr`, surfaced as `Rf_error(...calibrationMapName...)` at the bridge
([[RIB#bartcore_setCalibration]]) - is what the unexported
`dbarts:::bartcoreSetForestPriorScale` hits on a multinomial forest's low-level
handle (MEASURED [[test-calibration-midchain.R#"softmax calibration map"]]); the
R5 layer never routes a BCF sampler there since `refuseAmplitudeMutation`
refuses first, so only the multinomial arm exercises the bridge gate directly.
These cells stay `R`.

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and now
supports `extract(type = "loglik")` directly: [[generics.R#hurdleLogLik]]
combines the occupancy's `log(1 - pi)` / `log(pi)` with the positive part's
lognormal density (a `-log(y)` Jacobian against the stored log-scale channel)
at every row, reached from [[generics.R#extract.bartHurdle]]'s
`type == "loglik"` branch. This is NOT the sum of the two components' own loglik
channels - the positive fit's own channel covers only its y > 0 rows and
carries no Jacobian - but each component fit (`$occupancy` probit, `$positive`
gaussian) still supports `extract(type = "loglik")` independently too.

[f26] SHIPPED, tested bitwise. Nothing on the path gates on the coupling: the
shape probe ([[FAC#SamplerShape::supportsActiveRows]]) reports
`supportsActiveRows`, the mask composes into whatever precision the
installed response owns, and then into the per-forest weights at
[[CH#composeForestWeights]]. A K-forest chain's response is not necessarily
a `GaussianResponse`: it builds a `ProbitResponse` or a `LogisticResponse`
just as readily ([[CH#"switch (spec.family)"]]), and each overrides
`setActiveRows` on its own terms ([[MOD#ProbitResponse::setActiveRows]],
[[MOD#LogisticResponse::setActiveRows]]) exactly as it does off a coupling,
so [f15]'s logistic and probit coverage carries the latent K-forest case
with no edit of its own. Gaussian's composition into the case weights,
which is what inherits the sigma df, is
[[MOD#GaussianResponse::setActiveRows]]. The measurements below are all
GAUSSIAN two-forest ones; no latent K-forest mask is measured anywhere.
Measured on a 200-row two-forest sampler: `$setActiveRows(a)` and the bridge
`bartcoreSetActiveRows` are both accepted; on a sampler carrying `w` the
mask is BITWISE `setWeights(w * a)` in `train` and in `sigma`; an all-zeros
mask runs finite; a fractional element is refused. Tested:
[[test-active-rows-pins.R#"masked.bcf"]] (bitwise vs `setWeights(w * a)` on
train and sigma) and the `bcf-equivalence.R` `masked` scenario, carried by
the current `bcf-equivalence-3c81d6df.rds`. A per-forest mask is refused as
REDUNDANT rather than unbuilt: [[RIB#bartcore_setForestWeights]] already
expresses it - though that channel is deliberately NOT row removal
([[CH#Chain::setForestWeights]]: it does not remove the row from occupancy,
the combination or the sigma df; it DOES reach that forest's empty-leaf
veto, which counts positive composed weights). It is a PUBLIC R5 method,
[[dbarts.R#setForestWeights]], 1-based via `resolveForestIndex` (a BCF
basis forest is `2L`) and mirrored across re-creation through a dedicated
[[dbarts.R#reapplyForestWeights]] method, called from `getPointer` and
`setState`; the unexported `bartcoreSetForestWeights`, now
[[bartcoreHandle.R#bartcoreSetForestWeights]] (moved out of R/bartcore.R
along with every other low-level test-handle wrapper), stays the 0-based
internal wrapper the R5 method does not call.

[f27] Delegating / decorating: neither row needed an engine edit of its own
for the active-rows column. `GroupedResponse` forwards `setActiveRows` to
its base ([[MOD#GroupedResponse::setActiveRows]]) exactly as it forwards
`setWeights` ([[MOD#GroupedResponse::setWeights]]), advertising the base's
capability ([[MOD#GroupedResponse::supportsActiveRows]]), and
[[MOD#drawGroupEffects]] already weights its per-group sums by
`workingWeights()`, so an inactive row leaves its group's mean and
precision and an all-inactive group falls back to its prior through the
same formula. The heteroscedastic [[CH#formMeanWeights]] reads
`response_->workingWeights()` - the COMPOSED `w * a` while a mask is
installed - and divides by `s^2(x_i)`, so a zero stays a zero. Both are
tested bitwise against a composed-weight sampler: grouped's effects,
training fits and sigma all agree ([[tests/cpp/test_sampler.cpp#testActiveRows]];
an entirely inactive group draws its effect from the prior, finite);
heteroscedastic likewise for train and varcount
([[test-active-rows-pins.R#"heteroSampler"]],
[[test-active-rows-pins.R#"draws.hetero$varcount"]]). The same delegation
carries both halves of the nameable-calibration column for grouped:
[[MOD#GroupedResponse::fitScale, GroupedResponse::fitShift]] forward to
their base, so [[CH#Chain::resolvedNodeScale]] at creation and
[[CH#Chain::forestCalibration, Chain::setForestPriorScale]] mid-chain all
convert a named `prior.scale` exactly as they do for the undecorated
family, with no edit of grouped's own. Creation-time, grouped is one of
nine family/decoration paths a shared test measures
([[test-calibration-prior-draws.R#"grouped = anchorSampler"]] - a recorded
run measured 0.74210 against a 0.75 target); the mid-chain half rides the
same generic mechanism every non-coupled family does, with no dedicated
grouped case of its own in test-calibration-midchain.R.

[f28] A heteroscedastic fit also records `family = "gaussian"`, but the same
gaussian branch of `pointwiseLogLikelihood` reads its `s.train` surface
first: [[generics.R#heteroscedasticScale]] supplies the per-observation scale
whenever the fit carries one, taking precedence over the scalar `object$sigma`
([[generics.R#chainFastest]]) - the surface is stored on the fit in
[[bart.R#packageBartResults]]. As with [f19], pinned by value:
[[test-heteroscedastic-channels.R#"loglik scores at s(x)"]] checks it against
`dnorm(...)` at tolerance 1e-12, so the cell is `S`.

[f29] The variance forest is a separate leaf model entirely, outside
`forests_`, and is not addressable by the mid-chain `setCalibration` -
[[CH#Chain::setForestPriorScale]]'s bounds check `f >= forests_.size()`
never sees the variance forest, so a heteroscedastic sampler's
`shape.numForests` is 1 and `forest = 2` is refused by the ORDINARY
out-of-range check every single-forest sampler hits
([[RIB#bartcore_getCalibration, bartcore_setCalibration]]), not by a
hetero-specific gate. The MEAN forest's own calibration - both halves - is
not gated by that check: [[CH#Chain::resolvedNodeScale]] runs at the
[[CH#"forest.leaf.scale = resolvedNodeScale"]] assignment before the
variance-forest branch ([[CH#buildVarianceForest]]), and the mid-chain
reader/writer
([[CH#Chain::forestCalibration, Chain::setForestPriorScale]]) read
the same `response_->fitScale()`/`fitShift()` - none of the three reads any
family flag. Creation-time it is tested rather than merely built: a
heteroscedastic case is the tenth of
`test-calibration-prior-draws.R`'s `anchorSamplers` sweep
([[test-calibration-prior-draws.R#"heteroscedastic = anchorSampler"]]),
measuring a named `scale = 1.5` against the shared 0.75 target inside the same
9% band ([[test-calibration-prior-draws.R#"familyBand <- 0.09"]]) every other
family in that loop is held to
([[test-calibration-prior-draws.R#"abs(measured / priorSd - 1) < familyBand"]]).
Mid-chain, the mean forest rides the same generic mechanism
as every other single-forest family, with no dedicated
heteroscedastic case of its own beyond that shared mechanism and the
`forest = 2` refusal above.

[f30] Split verdicts, measured by benchmarks/R/composition-matrix.R. `resid.dist
= student()` + `variance =` now REFUSES at
[[spec.R#"does not support Student-t residuals"]] ("a variance forest does not
support Student-t residuals: the two are not yet shown to compose") - a
validation error only; the formal stays, and adjudicating the composition of
two scale mixtures on the same precision channel would reopen it. grouped +
`variance =` and hetero + grouped ranef are the same
construction reached from either spelling ([[CH#GroupedResponse]] decorates
before [[CH#buildVarianceForest]] builds the variance forest), and now REFUSE at
[[spec.R#"does not support grouped random effects"]] ("a variance forest does
not support grouped random effects: the group effects draw at a scalar residual
scale, which the variance forest replaces row by row"), with a
[[RIB#createHolder]] backstop closing the entrances that never reach spec.R;
[[docs/design/heteroscedastic.md#refused with a variance forest]]
records what either adjudication would need.

[f31] Recorded but UNBUILT extensions, refused with that reason in the comment:
grouped ordinal because the cutpoint block and the group block are not yet shown
to interleave ([[RIB#"not supported for ordinal responses"]], ordinal.md
section 8), grouped nbinom the same for the dispersion block
([[RIB#"not supported for count (nbinom) responses"]], negative-binomial.md
section 7).

[f32] No surface at all: [[RIB#applyGroupAttribute]] is called from exactly
one site, in [[RIB#createHolder]], on the single-forest holder path, so
`bartcore.groups` is
never read for a multinomial sampler.

[f33] [[CH#buildMultinomialForest]] hard-sets
`forest.useDart = false`, and [[RIB#buildMultinomialSampler]] copies only
power/base/proposal-probability fields, so a DART tree prior built from either
the `dart` argument or a `tree.prior` object never reaches the K-forest engine.
`bart2` refuses both routes by name
([[bart.R#"'dart' or a DART 'tree.prior'"]]), matching BCF's own named refusal
([[spec.R#"a DART tree prior"]], [[CH#buildSpecifiedForest]]) before either
reaches the host sampler.

[f34] [[bart.R#bart2Hurdle]] builds both component calls with
[[bart.R#redirectCall]], so a user's `variance =` is forwarded to BOTH -
including the occupancy component, which then sets `family = "probit"` and hits
the non-gaussian variance refusal
([[spec.R#"a variance forest requires family"]]) before either component
fits.
That refusal is deliberate, not a bug: the positive fit is always
homoscedastic because the gate makes a heteroscedastic component
unreachable. No comment in R/ states that; the evidence is the live
[[bart.R#redirectCall]] /
[[spec.R#"a variance forest requires family"]] mechanism above and the Rd text,
[[bart2.Rd#variance]] and [[dbarts.Rd#variance]].

[f35] `dart` is forwarded to both components ([[bart.R#redirectCall]]), each of
which is an ordinary single-forest chain that takes it.

[f36] No family gate: [[SAM#Sampler::installForests]] checks shape, grid, DART,
and the variance forest's presence and saved slot, and matches donor forest
counts, and [[CH#growForestFromRoot]] loops every forest, with a
variance-forest pre-step above the loop. Neither is exercised by a BCF test,
and BCF has no `bart2()` surface, so both are reached only through the R5
[[dbarts.R#installTrees]] / [[dbarts.R#growFromRoot]].

[f37] `rbart_vi()` carries no `warm.start` or `n.grow.sweeps` formal
([[rbart.R#rbart_vi]]) and, since it also carries no `...` formal at all, an
unrecognized `warm.start =`/`n.grow.sweeps =` argument now hits ordinary R
argument-matching ("unused argument") rather than a dedicated
unknown-argument check - the same retirement as `bart2()`'s `treatment =`
case (see [f7]). The underlying R5 sampler carries no group gate on either
path, so this is a surface gap, not an engine one.

[f38] The MEAN forest keeps DART; the variance forest never takes it
([[CH#buildVarianceForest]] never sets `useDart`, whose default is false in
[[CH#SamplerOptions]]).

[f39] Current baselines: `equivalence-d4bca4ce.rds` (51 scenarios),
`bcf-equivalence-3c81d6df.rds` (12), `multinomial-equivalence-4d9a3337.rds` (11)
- benchmarks/baselines/MANIFEST. Scenario names are the keys in
[[benchmarks/R/equivalence.R#makeScenarios]].

[f40] docs/plans/sbc-family-tiers.md (status BUILT) plus
docs/plans/sbc-calibration.md (DONE). The A/B/C "tiers" in the latter are
FEATURE tiers (A baselines/DART/grouped/weighted/BCF, B linear leaf, C GP leaf),
not family tiers - there is no per-family tier ladder, only per-family recorded
verdicts, which is what this column carries.

[f41] gamma3 flagged in one stream and RESOLVED as the cutpoint-vs-mean-level
ridge mixing slowly (does not reproduce across streams; 0.31 of the band at 3x
the chain length). The cutpoint block, the latent eta and all K category
probabilities calibrate.

[f42] `avg.mu` - the identified mean - passes cleanly; `r` and `agg.psi` flag at
thin = 30 and cross into the band at 5x the spacing. Read as slow mixing
along the r-psi ridge, measured at two thinning settings rather than
three; a third, longer run is still owed.

[f43] Aggregate `p_k(x*)` and the three raw per-forest `f_ik` cells all pass.
[[COM#MultinomialForestCombiner::afterCombine]] draws the level from its exact
leaf-space conditional; the acceptance run at R = 200,
`Rscript benchmarks/R/sbc.R multinom 200 150 30`, scores every functional PASS
at band 0.1282, the three cells at 0.0688/0.0824/0.0675
([[docs/plans/archive/multinomial-level-centering.md#Landing]]). The function that
ranks those cells is [[benchmarks/R/sbc.R#cellNames]].

[f44] AFT is exercised only in combination, through the `grouped_aft` scenario
([[benchmarks/R/equivalence.R#grouped_aft]]). There is no standalone AFT
equivalence scenario; the
separate exact oracle benchmarks/R/aft-exact.R is not a MANIFEST entry.

[f45] Out of the SBC matrix by scope, each for its own recorded reason
([[docs/plans/sbc-family-tiers.md#Decision - scope]]): aft because its
censoring status is fixed at
creation, so a prior-draw replication cannot vary it (the enabler is a status
setter); hazard and hurdle because their person-period / two-part designs depend
on `y0`, which breaks exchangeability, and because neither owns any sampling
code.

[f46] Tier A PASS with the sigma channel resolved as slow mixing along the
(a, mu) ridge ([[docs/plans/sbc-calibration.md#Final summary]]). Explicitly out of the
family-tiers matrix. `runSbcBCF` once errored; it is now repaired
(docs/plans/archive/runsbcbcf-repair.md, acceptance PASS across
thin=30/90/120).

[f47] OUT but DEFERRED rather than blocked: prior draws never reach
`varianceForest_` today, and the capability is liftable R-side through `setState`
([[docs/plans/sbc-family-tiers.md#Decision - scope]]).

[f48] The K-forest coupling's family reach (docs/plans/archive/multiforest-extension-surface.md). gaussian, probit and logistic
build; aft, ordinal and nbinom are refused at all three creation routes, each
naming what it is missing
([[spec.R#"a treatment forest does not support family"]],
[[RIB#refusedAmplitudeFamilyReason]], [[FAC#createAmplitudeSampler]] - the last
sitting directly beside the variance-forest refusal in that same factory, which
is unchanged and family-independent). The calibration map's anchor is now
family-keyed, [[CH#latentScaleAnchor]]: sd(y) under gaussian, 1 under probit,
pi/sqrt(3) under logistic, and stated per unit of basis row norm
([[CH#basisRowNorm]]). Cell by cell in section 2: `setResponse`/`setOffset` are
OPEN under every family, because [[CH#Chain::setResponse]] now hands the
response `combinedFits()` rather than forest 0's bare totals, which is what let
the gaussian conjunct come off [[CH#Chain::supportsResponseMutation]]; the
combiner's own opt-in
([[COM#AmplitudeForestCombiner::supportsResponseMutation]]) is unchanged.
`setSigma` is REFUSED for probit and logistic; `setWeights` is refused for
probit and OPEN for logistic, whose counts are its Polya-Gamma shape - both
through the ORDINARY single-forest guards now that the sampler answers
`shape.family` for itself ([[RIB#refuseBinaryWeightChange]],
[[RIB#refusePinnedSigmaChange]]). The shared [[RIB#enforceBinaryWeightPolicy]]
refuses a probit weight outright and holds a logistic one to positive integer
counts, at creation and on every mutation conduit alike, which is what makes the
zero-weight-subset cell family-dependent too. `updateScale = TRUE` stays REFUSED
under EVERY family - NOT the latent convention `- [f9]` the original design
proposed: [[bartcore.R#refuseAmplitudeMutation]] keys on the sampler carrying
bases, never on the family, and the bridge's
[[RIB#refuseMultiForestResponseMutation]] keys on `numForests >= 2`, so a probit
K-forest is refused too, though its transform is the identity and the
re-anchoring the refusal guards against cannot occur. Tested, with the open
conduit and both refusals above, at [[test-bcf-family.R#"binaryFits"]]. That
file is the whole of the latent K-forest's evidence: no equivalence scenario and
no SBC coverage reaches one.

[f49] The test-surface cell stays `R`, and truthfully: what shipped is the
PER-FOREST replay ([[RIB#bartcore_predictPerForest]],
[[CH#Chain::predictPerForestFromSavedSample]]), not a test surface. The
resident test store, `run()$yhat.test` and the SAMPLER's own `predict()` remain
refused through [[RIB#refuseUndefinedTestFits]] because the blend
`sum_f dot(a_f, B_f(i, .)) f_f(x_i)` needs an off-sample basis the sampler does
not have; the entry sidesteps that by reporting the `f_f(x_i)` alone and
leaving the contraction to the caller, whose bases they are. What is no longer
refused is the FIT-level combined `predict()`: it performs that
contraction in R (`predictBlend`, R/generics.R), from the replay, the packaged
`glue` and either a `bases =` at the caller's own rows or a `forest()` term's
formula re-evaluated against `newdata` - which is precisely the off-sample
basis the sampler lacks and the caller holds. Evidence:
`inst/tinytest/test-predict-forest.R` (the replay-at-training-rows identity
against the in-sample channel at 1e-12, the recombination identity against
`yhat.train`, and the offset / no-reporting / multinomial refusals),
`inst/tinytest/test-predict-blend.R` (the blend-at-training-rows identity
against `yhat.train` at 1e-12 on gaussian, probit, logistic, two-chain
uncombined, q = 3 and both-forests-with-a-basis fixtures; bitwise identity with
the documented manual recombination on both counterfactual cases; the bart2
term route's auto-derivation against an explicit `bases =`; and the shape,
keepTrees, offset and weights refusals), plus `tests/cpp/test_sampler.cpp`'s
`testAmplitudePerForestReplay` (both replay routes against `forestTotalFits`,
and the raw-scale pin against `predict`).

The IN-SAMPLE per-forest OUTPUT channel (as opposed to construction, above)
shipped separately: `packageBartResults` packages `forestFits`/`glue` -
response-scale raw per-forest totals and the ragged multiplier channel -
onto any bcf-shaped fit regardless of how it was built, with `forest1..K`
names and a `forest.labels` attribute, and
`extract(type = "forest", forest =, contribution = )` serves it (raw slice
or the on-demand `basis %*% glue` contribution) under both `combineChains`
conventions; refused by name on a fit without forest reporting or with
`sample = "test"`. `inst/tinytest/test-formula-terms.R` adds its own
identity checks (a term-built two-forest fit is byte-identical to the
equivalent `forests =` one at the same seed) as further bitwise evidence for
the gaussian/probit/logistic construction route, though it is not an
equivalence-harness scenario and does not change the counts in table 5.
Out-of-sample per-forest replay is a separate, later addition: a new engine
entry replays every forest at new rows RAW - no `fitScale`, no `fitShift`, no
offset - and `predict(type = "forest", forest =)` and the R5
`$predictForests()` serve it, gated on `forestReportingIsDefined` like the
in-sample channel. It does not open a test-basis path: the resident test
surface, `run()$yhat.test` and `extract(type = "forest", sample = "test")` stay
refused for an amplitude coupling, since the RECOMBINATION off the training rows
needs bases only the caller has
([[docs/plans/archive/bart2-argument-consolidation.md#5.4 subset, weights, offset, test, missing]]).
The section-2 `bcf` test-surface cell is unchanged for that reason.

The recombination itself is in R rather than the engine, and does not touch
that cell: `predict(type = "ev"/"ppd"/"bart")` on an
amplitude-coupled fit blends the per-forest replay with the packaged `glue` and
the response shift, taking the bases at the predicted rows from a `bases =`
argument or, for a `bart2` `forest()` term, re-evaluating the declaring formula
against `newdata` under the fit's own factor levels (`basis.terms` on the fit).
The engine is untouched, so the sampler-level closure above still holds; what
changed is that the caller no longer writes the three-line contraction by hand.
Both amplitude routes now require `keeptrees`/`keepTrees`, which
`predict(type = "forest")` silently did without - it reported the current trees
once for every draw, a shape no amplitude pairing exists for.

## Gaps

Every MISSING (`M`) cell above, as a candidate work item, grouped by family.
Scheduling is VD's and the orchestrator's; nothing here carries a schedule,
and REFUSED cells are deliberately absent - they are part of the models.

**gaussian.** None; every column is S or an intentional `-`. The one standing
inexactness (zero-weight rows survived the empty-leaf veto) closed at
`empty-leaf-veto-fix` ([f17]).

**student (Gaussian + Student-t residuals).** No `rbart_vi()` surface
([[rbart.R#"family <- match.arg(family)"]]); no `xbart()` surface
([[xbart.R#xbart]]). Pointwise loglik is
pinned by value ([f19], now `S`), so it is not a gap. Composition with a
variance forest is refused by name ([f30]). Only one dedicated tinytest
file.

**probit.** None.

**logistic.** No `rbart_vi()` token ([[rbart.R#rbart_vi]]), so grouped logistic
is engine-reachable but not R-reachable.

**ordinal.** No `xbart()` or `rbart_vi()` reach (both refused in
[[data.R#resolveClassificationFamily]] as unsupported response shapes). Grouped
ordinal is a recorded unbuilt item
([[RIB#"not supported for ordinal responses"]]). `warm.start` and
`n.grow.sweeps` unbuilt for the arc ([[bart.R#checkFamilyUnsupportedArgs]]).
Its selecting control attribute is undocumented in the shipped header ([f3]).
One dedicated tinytest file. SBC gamma3 resolved but not re-run at full R.
Pointwise loglik ([[generics.R#ordinalLogLik]]) is no longer a gap - built,
though untested beyond the shape checks in test-ordinal.R.

**nbinom.** As ordinal, plus: no `bart()` token, no `xbart()` token, no
`rbart_vi()` token. Grouped nbinom a recorded unbuilt item
([[RIB#"not supported for count (nbinom) responses"]]). Header
attribute undocumented. Two dedicated tinytest files (test-nbinom.R,
test-dispersion-channel.R). SBC `r`/`agg.psi` flag
standing (read as slow mixing along the r-psi ridge; a third, longer run
is owed). Real-valued dispersion remains
a recorded open item (TODO's `negbin-real-dispersion` item). Pointwise loglik
([[generics.R#negbinLogLik]]) is no longer a gap.

**multinomial.** No flat-C creation path at all ([f4]); `dbarts()`,
`dbartsSpec()` and `bart2()` all build the R5 `dbartsSampler` directly
(multinomial-mutation-arc.md), one `bartcore_create`, `$fit` the
engine that ran - no host shell, no `$bc` ([f11]). Its three response
channels (`setCounts`, `setCategoryOffset`, `setCategoryTestOffset`) are
public R5 methods. No `getLatents` ([f22]). No grouped surface ([f32]). No
warm start / grow-from-root. DART, `split.probs`, `monotone`, and `variance`
are all refused by name rather than silently dropped ([f33]). Pointwise
loglik ([[generics.R#multinomialLogLik]]) is no longer a gap - it scores
the multinomial log-pmf on the REPORTED probabilities, distinct from the
engine's own per-observation channel, which stays undefined for this family.
SBC raw `f_ik` is no longer open ([f43]).

**aft.** No `xbart()` token. Pointwise loglik ships, but `setWeights` is refused
and the censoring status is fixed at creation, which is also what keeps AFT out
of the SBC matrix - a status setter is the named enabler ([f45]). No standalone
equivalence scenario ([f44]).

**hazard.** No flat-C token and no engine code of its own ([f6]); no
`rbart_vi()`, `xbart()` or `dbartsSpec()` reach. Out of SBC by design. One
dedicated tinytest file.

**hurdle.** No `dbarts()` sampler by construction ([f12]); no `rbart_vi()`,
`xbart()` or flat-C reach. No warm start / grow-from-root. A heteroscedastic
positive part is REFUSED via the occupancy component's own gate, deliberately,
not a partial feature ([f34]). Top-level pointwise loglik
([[generics.R#hurdleLogLik]]) is no longer a gap.

**bcf.** No `bart2()` surface for the NAMED `bcf()` causal verb, by VD's
resolution that puts it in bartCause; the general K-forest amplitude
capability itself IS reachable from `bart2()`'s formula interface via a
`forest()` term ([f7]). Warm start and grow-from-root are unrefused
and untested for two forests ([f36]). Whole-data `setData` stays undesigned
(open question 1 of the model-space survey, docs/design/model-space-survey.md).
`setForestWeights` now has a public R5
method (multiforest-extension-surface.md; [f26]) - no longer a gap.
The probit and logistic cases have no equivalence scenario, no SBC
coverage and no measured active-rows mask ([f48], [f26]); aft, ordinal and
nbinom are recorded open items, not gaps. Pointwise loglik on a bcf fit runs (a BCF-shaped fit's `family` is
gaussian, probit or logistic, so `extract(type = "loglik")` reaches the
shared [[generics.R#pointwiseLogLikelihood]] dispatcher unrefused), but has
never been checked to be the right quantity for a combined per-forest
location.

**grouped.** `rbart_vi()`-only surface ([f8]): no `bart2()`, `dbarts()`,
`xbart()` or `dbartsSpec()` reach, and no `warm.start` / `n.grow.sweeps`
formals ([f37]) though the engine paths carry no group gate. The `setResponse`
gap CLOSED at adoption-slate.md ([f13]; the section-2 cell is the record).
Composition with a variance forest now REFUSES at
[[spec.R#"does not support grouped random effects"]] ([f30]).

**heteroscedastic.** No `xbart()` reach. Pointwise loglik reads the
per-observation `s.train` surface it stores and is tested by value
([f28], now `S`). Selecting attribute undocumented in the
header ([f3]). Out of the SBC matrix, deferred not blocked ([f47]). The one
item its own design recorded as unbuilt - the `setState` variance column-mask
gap - is BUILT: [[CH#Chain::columnMaskStateFeasible]] carries a variance pass
of its own over the state's variance trees, [[CH#rebuildVarianceForest]] holds
every restored variance tree to the forest's mask, and both install entries
gate on the predicate - [[SAM#Sampler::setState]] and
[[SAM#Sampler::installForests]] - each surfacing the one refusal by name
([[RIB#columnMaskMismatchMessage]]).
