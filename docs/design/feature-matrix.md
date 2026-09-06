# Response-model feature matrix

Status: living reference, updated in place whenever a cell changes.

What each shipped response model can and cannot do, one row per model and one
column per capability that bears on scheduling. Every SHIPPED and REFUSED cell
carries a cite verified against the live tree; a wrong cell misdirects
scheduling.

## Legend

| code | meaning |
|---|---|
| `S` | SHIPPED. Works today; the cite is the site that makes it work. |
| `R` | REFUSED on model or identification grounds - the refusal is part of the model, not a hole. The cite is the refusal site. |
| `P` | PLANNED. A named design arc covers it; the cell names which stage. No cell currently carries this code. |
| `M` | MISSING. Not built, no schedule. A cite, when given, is a guard that errors *because the thing is unbuilt* - a recorded open item, not a model refusal. |
| `-` | N/A. The concept does not apply to this row; the row footnote says why. |
| `?` | UNVERIFIED. Constructs today with no refusal site, but no test, doc or adjudication backs it. Do not schedule against the cell until it is settled; every `?` is listed under "Gaps". No cell currently carries this code. |

A marker written `[fN by family]` flags a cell whose code varies with the base
family the row is built over; the footnote says which family takes which.

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

Cites are by symbol and are existence-checked by `tools/check-doc-freshness.R`; a cell's VALUE is adjudicated separately from its cite.

## Rows

Twelve rows. The first nine are response models proper; the last three are
compositions, couplings or decorations over a base response that a user
selects the same way and schedules against the same way, so they earn rows.
Hurdle is one of the four: it has no sampler of its own and is an R-side
composition of two ordinary fits ([f12]).

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

The engine's `ResponseFamily` enum has only six tokens ([[MOD#ResponseFamily]]:
gaussian, probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf
and hetero are all reached some other way, which is exactly why they need rows
here rather than an enum read.

Two rows sit on a COMBINER, equivalently a COUPLING: an object holding K
forests and the rule that combines their fits into one location per
observation, kept alongside the response model rather than inside it
([[COM#ForestCombiner]], [[docs/design/forest-combiner.md#The ForestCombiner]]).
Multinomial's coupling is a softmax over K category forests. The bcf row's is
the AMPLITUDE family: forest f's fit at row i is scaled by
`dot(a_f, B_f(i, .))`, the contraction of that forest's own BASIS row with its
own amplitude vector, and the per-forest amplitude block a fit packages for
later recombination is its GLUE
([[docs/design/multiplier-combiner.md#The amplitude layout]]). That row is
family-dependent below because a K-forest chain selects its response model
from `AmplitudeSpec::family` ([[COM#AmplitudeSpec::family]]) through the
`switch (spec.family)` arm of its constructor ([[CH#"switch (spec.family)"]]).
Leaf models (constant, monotone, linear, GP) are an orthogonal axis and are not
rows; where a leaf model gates a capability the cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|
| gaussian | S [[bart.R#bart]] | S [[bart.R#bart2, gaussian]] | S [[dbarts.R#dbarts, gaussian]] | S [[xbart.R#xbart, gaussian]] | S [[CAPI#DBARTS_FAMILY_GAUSSIAN]] |
| student | S [[bart.R#residDist]] | S [[bart.R#bart2, resid.dist]] | S [[dbarts.R#dbarts, resid.dist]] | M [[xbart.R#xbart]] | S [[RIB#parseSamplerSpecification, residualDf]] [f2] |
| probit | S [[bart.Rd#y.train]] | S [[bart.R#bart2, probit]] | S [[dbarts.R#dbarts, probit]] | S [[xbart.R#xbart, probit]] | S [[CAPI#DBARTS_FAMILY_PROBIT]] |
| logistic | S [[bart.R#bart, logistic]] [f1] | S [[bart.R#bart2, logistic]] | S [[dbarts.R#dbarts, logistic]] | S [[xbart.R#xbart, logistic]] | S [[RIB#resolveFamily, logistic]] |
| ordinal | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Ordinal]] | S [[dbarts.R#dbarts, ordinal]] | R [[data.R#resolveClassificationFamily]] | S [[RIB#resolveFamily, ordinal]] [f3] |
| nbinom | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Negbin]] | S [[dbarts.R#dbarts, nbinom]] | M [[xbart.R#xbart]] | S [[RIB#resolveFamily, nbinom]] [f3] |
| multinom | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Multinomial]] | S [[dbarts.R#dbarts, multinomial]] [f4] | R [[data.R#resolveClassificationFamily]] | M [f4] |
| aft | S [[bart.R#bart, aft]] [f1] [f5] | S [[bart.R#bart2, aft]] | S [[dbarts.R#dbarts, aft]] | M [[xbart.R#xbart]] | S [[CAPI#DBARTS_FAMILY_AFT]] |
| hazard | R [[bart.R#refuseBartRedirectedFamily]] [f1] | S [[bart.R#bart2, hazard]] | S [[dbarts.R#dbarts, hazardTokens]] | M [[xbart.R#xbart]] | M [f6] |
| hurdle | R [[bart.R#refuseBartOwnClassFamily]] [f1] | S [[bart.R#bart2Hurdle]] | R [[dbarts.R#"is only available through bart2()"]] | M | M [f6] |
| bcf | - [f1] | S [[formulaTerms.R#ingestFormulaTerms]] [f7] | S [[dbarts.R#dbarts, forests]] | M | S [[CAPI#dbarts_sampler_create]] |
| hetero | - [f1] | S [[bart.R#bart2, variance]] | S [[dbarts.R#dbarts, variance]] | M | S [[RIB#applyVarianceAttributes]] [f3] |

[[spec.R#dbartsSpec]] resolves the seven single-forest tokens - auto,
gaussian, probit, logistic, aft, ordinal, nbinom - plus the K-forest amplitude
family through its `forests =` argument ([[spec.R#dbartsSpec, forests]]), each
forest declared as `forest(basis = ...)`, a FOREST BASIS being the n x q matrix
whose row for observation i contracts with that forest's amplitude vector into
the scalar the forest's fit is multiplied by. A variance forest comes through
`variance =` ([[spec.R#dbartsSpec, variance]]), and the `family` formal also
accepts `"multinomial"` directly ([[spec.R#dbartsSpec, multinomial]], with
dedicated body logic in [[spec.R#resolveSamplerSpec, unsupportedMultinomial]]);
only hazard and hurdle stay out of `dbartsSpec()`'s reach. A
`forests =` fit resolves **gaussian, probit or logistic**; aft, ordinal and
nbinom are refused there by name, each stating what it is missing
([[spec.R#"a treatment forest does not support family"]], whose `multinomial =`
arm refuses a `forests =` declaration for that family too - "its forests are
its categories... not an amplitude coupling"), with the same three-family gate
at the bridge ([[RIB#refusedAmplitudeFamilyReason]], called from both creation
routes, [[RIB#refuseUnsupportedAmplitudeComposition]] and
[[RIB#createBCFHolder]]) and at the factory ([[FAC#createAmplitudeSampler]]).

`bart2()` reaches that same `forests =` machinery through a `forest()` term in
its formula rather than a formal of its own, under an identical family gate.
See [f7].

## 2. Mutation channels on the R5 `dbartsSampler`

The channels that make dbarts a conditional model inside an outer sampler; the
two families of mutation and their transaction semantics are
[[docs/architecture.md#The mutation surface]].
"Test surface" means the test-set predictor and offset conduits and the fits
they feed, not unit tests. `updateScale` is broken out because it is refused
independently of the setter it rides on. Probit, ordinal, aft and nbinom carry
no case weights at all, so `setWeights` is refused for them on every surface
([[RIB#refuseBinaryWeightChange]], [[RIB#familyCarriesNoWeights]]); logistic's
weights are its Polya-Gamma trial counts and stay open ([f10]).

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
| bcf | S [[CH#Chain::setResponse]] [f48] | S [[CH#Chain::setOffset]] [f48] | R [[bartcore.R#refuseAmplitudeMutation]] [f48] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] | S [[RIB#bartcore_setWeights]] [f48 by family] | S [[RIB#bartcore_setSigma]] [f48 by family] | R [[RIB#refuseUndefinedTestFits]] [f49] |
| hetero | S [[RIB#bartcore_setResponse]] | S [[RIB#bartcore_setOffset]] | R [[RIB#refuseVarianceForestScaleUpdate]] | S [[RIB#bartcore_setPredictor, bartcore_updatePredictorPerObservation]] | S [[RIB#bartcore_setWeights]] | R [[RIB#refusePinnedSigmaChange]] | S [[RIB#bartcore_setTestPredictor]] |

`setData` (whole-data replacement, n free) is dense-store and single-forest
only ([[RIB#refusePredictorMutation, refuseMultiForestMutation]]) and is
refused for aft ([[RIB#"fix the censoring structure at creation"]]);
BCF/multinomial whole-data `setData` is undesigned (model-space-survey.md,
open questions 1 and 3).

## 3. Row subsetting, latents, calibration

"Zero-weight row subset" asks whether a fit with zeroed rows equals a fit on
the retained rows; "active-rows mask" is the mid-chain per-observation channel
that answers the same question when weights cannot ([f15]); "nameable
calibration" is a named leaf-prior scale in response units ([f16]).

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
| bcf | S [[COM#AmplitudeForestCombiner::formForestResponse]] [f17] [f48 by family] | S [[MOD#GaussianResponse::setActiveRows]], [[CH#composeForestWeights]] [f26] | S [[CH#Chain::latents]] [f18 by family] | S [[generics.R#pointwiseLogLikelihood]], [[test-bcf-loglik.R#"the COMBINED one"]] | R [f23] |
| hetero | S [[CH#sweepVarianceForest]], [[MOD#ConstantVarianceLeaf::accumulate]] | S [[CH#formMeanWeights]] [f27] | - [f18] | S [[generics.R#heteroscedasticScale]] [f28] | S [[CH#Chain::resolvedNodeScale]], [[CH#Chain::forestCalibration, Chain::setForestPriorScale]] [f29] |

## 4. Model composition

| model | variance forest | DART | warm start | grow-from-root |
|---|---|---|---|---|
| gaussian | S [[FAC#createSampler]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| student | R [[spec.R#"does not support Student-t residuals"]] [f30] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| probit | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| logistic | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| ordinal | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] |
| nbinom | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] |
| multinom | R [[bart.R#unsupported]] | R [[bart.R#"'dart' or a DART 'tree.prior'"]] [f33] | M [[bart.R#checkFamilyUnsupportedArgs]], [[RIB#refuseMultiForestWarmStart]] [f50] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] |
| aft | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| hazard | R [[spec.R#"a variance forest requires family"]] | S [[CH#useDart]] | S [[bart.R#warm.start]] | S [[dbarts.R#growFromRoot]] |
| hurdle | R [[spec.R#"a variance forest requires family"]] [f34] | S [[bart.R#redirectCall]] [f35] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] | M [[bart.R#checkFamilyUnsupportedArgs]] [f50] |
| bcf | R [[FAC#createAmplitudeSampler]] [f48] | R [[spec.R#"a DART tree prior"]] | M [[RIB#refuseMultiForestWarmStart]] [f36] | S [[CH#growForestFromRoot]] [f36] |
| hetero | - | S [[CH#useDart]] [f38] | S [[SAM#Sampler::installForests]] | S [[CH#growForestFromRoot]] |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused in [[dbarts.R#growFromRoot]] and a no-op in
[[CH#growForestFromRoot]], so every family above
reads "constant leaf" in that column.

## 5. Evidence

| model | equivalence baseline [f39] | SBC verdict [f40] | dedicated tinytest files |
|---|---|---|---|
| gaussian | 33 scenarios (friedman, weighted, splitprobs, chains, setdata, wtoffset, quants, categorical, missing, dart, linear, gp, zeroweights, sparse, wtgp, chik2, predswap, predcol, predpartial, predreject, predforce, ordfactor, nafactor, sparsefactor, testswap, leaffactor, leaffactormixed, factorpartial, xbartmixed, xbart, bart2gauss, bart2twoforest, mixedmatrix) | PASS 7/7 | ~20 (test-sampler-*.R, test-bart-bart2.R, test-zero-weights.R) |
| student | `student` | PASS 4/4 | test-robust-errors.R only |
| probit | `probit`, `chik`, `maskprobit`, `bart2probit` | PASS | test-binaryResponse-hyperprior.R, test-family.R, test-weighted-binary-ppd.R |
| logistic | `logistic`, `wtlogistic` | PASS 6/6 | test-weighted-logistic.R, test-family.R |
| ordinal | `ordinal`, `maskordinal` | 9/10 [f41] | test-ordinal.R only |
| nbinom | `nbinom` | 1/3 [f42] | test-nbinom.R, test-dispersion-channel.R |
| multinom | 11 scenarios, own harness | aggregate PASS, raw `f_ik` PASS [f43] | 6 (test-multinomial-*.R) |
| aft | `aft` | OUT [f45] | test-aft.R only |
| hazard | `hazard` | OUT [f45] | test-hazard.R only |
| hurdle | `hurdle` | OUT [f45] | test-hurdle.R, test-hurdle-surface.R |
| bcf | 12 scenarios, own harness, gaussian only [f48] | PASS, gaussian only [f46] | 9 (test-bcf*.R) |
| hetero | `hetforce`, `hetswap`, `hetpartial` | OUT [f47] | 4 (test-heteroscedastic*.R) |

The rows are keyed by response model; predictor SHAPE cuts across them. Eight
gaussian scenarios carry the factor shapes - `ordfactor`, `nafactor` (the
anchor for MIA, missingness incorporated in attributes,
[[docs/design/mia-missingness.md#Model]]), `sparsefactor`, and the five
store-path scenarios `testswap`, `leaffactor`, `leaffactormixed`,
`factorpartial` and `xbartmixed` - four of them carrying an ordered-factor
predictor, so the ordered-factor cut grid rests on those plus
`benchmarks/R/categorical-exact.R`'s ordered-factor case.

`inst/tinytest/test-capi.R` drives the whole single-forest family list through
`dbarts_sampler_create` - `""`/`"probit"`, logistic, ordinal, aft and nbinom,
each run and checked for finite, correctly-shaped output - plus
heteroscedastic ([[test-capi.R#"bartcore.variance"]]) by control attribute,
and BCF
([[test-capi.R#"zBCF"]]) through `forests = list(forest(basis = ...))`.
Multinomial has no flat-C creation path to test ([f4]).

## Footnotes

[f1] `bart()`'s `family` formal is the narrow `c("auto", "logistic", "aft")`
([[bart.R#bart]]), forwarded to `dbarts()` verbatim; `"auto"` (the default)
infers gaussian or probit from the response, and `resid.dist` remains the
separate Student-t lever, untouched by `family`. All ten tokens outside that
vocabulary are refused BY NAME ahead of `match.arg`'s generic message, by the
shared [[bart.R#refuseBartRedirectedFamily, bartRedirectedFamilies]]: five
own-class tokens - `"multinomial"`, `"ordinal"`, `"nbinom"`,
`"hurdle.lognormal"` and `"twopart"` - redirect through
[[bart.R#refuseBartOwnClassFamily, bartOwnClassFamilies]] (ordinal's own
auto-detection off an ordered-factor response routes through the same helper),
and five more - `"gaussian"`, `"probit"`, `"hazard.probit"`, `"hazard"` and
`"hazard.logistic"` - get their own named message, the first three because
`"auto"` already gives no extra capability and the last two because the
discrete-time expansion needs `breaks`/`max.rows`, which `bart()` has no formal
for (both reasons stated in man/bart.Rd's `family` item). `hazard`'s section-1
`bart()` cell therefore reads `R`, deliberately refused by name. `bcf` and
`hetero` are not `family` tokens at all and stay out of reach by signature,
which is why their cells read `-`.

[f2] Student-t is not a `family` token on any R entry point, and it is not in
`dbarts_sampler_create`'s admission list: it is selected instead by a finite
`resid.df` attribute on the model SEXP
([[RIB#parseSamplerSpecification, residualDf]], a gaussian-only gate) and
refused for every non-gaussian family R-side by
[[spec.R#"student residuals require a continuous gaussian response"]]. The
shipped header does carry an enumerator, [[CAPI#DBARTS_FAMILY_STUDENT]], which
the per-observation augmentation entries read. The engine family stays
`gaussian`, which is why the whole gaussian row applies to student in tables
2-4.

[f3] Ordinal and nbinom are reachable through `dbarts_sampler_create` and ship
a `DBARTS_FAMILY_*` enumerator each; heteroscedastic is a decoration with no
enumerator of its own, selected instead by a control attribute. All three
selecting attributes are documented in the shipped header, in the
specification-attribute block on `dbarts_sampler_create`
([[CAPI#"SPECIFICATION ATTRIBUTES"]]): ordinal's `bartcore.n.categories`
([[RIB#parseControl, bartcore.n.categories]]), nbinom's `bartcore.dispersion`
([[RIB#parseControl, bartcore.dispersion]]) and heteroscedastic's
`bartcore.variance` ([[RIB#applyVarianceAttributes]]). The block states each
one's type and shape, its allowed values, the family or composition it selects
and what a malformed value does, and covers every other attribute creation
reads off the control or model beside them: `bartcore.survival`,
`bartcore.forests`, `resid.df`, the monotone directions, and the interaction
and block partitions.

[f4] `dbarts(x, y, family = "multinomial")` (matrix interface only) accepts a
counts matrix or a factor/character/integer-code response, one-hot expanded
([[dbarts.R#dbarts, multinomial]]), with the counts matrix built by
[[data.R#resolveMultinomialCounts]]. Creation routes through the same public
dispatch every family uses, [[RIB#bartcore_create]], whose multinomial arm is
[[RIB#createMultinomialDataHolder]]; there is no separate multinomial creation
entry point, and no dbarts.h one at all -
[[C_interface.cpp#creationFamilyName]] refuses the multinomial token as a
sampler that entry cannot build.

[f5] `family = "aft"` is an explicit token on `bart()`, documented in
man/bart.Rd's `family` item; the underlying `Surv()` / two-column-`y.train`
detection is a separate mechanism and resolves the response shape on its own.

[f6] `family = "hazard"` / `"hazard.probit"` / `"hazard.logistic"` is
person-period ingestion sugar: [[dbarts.R#expandDiscreteTimeHazard]] expands the
design and [[dbarts.R#"the remap: the engine-facing family is now an ordinary binary link"]]
remaps the token to `"probit"` or `"logistic"` before any model is built. The
resulting sampler *is* an ordinary binary one, so its whole row equals the
probit (or logistic) row and the fit records `family = "probit"`; there is no
engine code of its own, hence no C-API token and no SBC arm. Refusals inside
the expander: no formula interface
([[dbarts.R#"discrete-time hazard fits currently use the matrix interface"]]),
no `subset` ([[dbarts.R#"survival responses do not support 'subset'"]]), no
`test` ([[dbarts.R#"discrete-time hazard fits do not take a 'test' set"]]).

[f7] `treatment` is not a `bart2()` formal, and `bart2()`'s formal list
([[bart.R#bart2]]) carries no `...` at all, so `treatment =` hits ordinary R
argument matching ("unused argument") rather than a dbarts-authored construct.
The K-forest amplitude capability this row names is reached from `bart2()`'s
FORMULA interface instead: a `forest()` term - `z:forest(x1 + x2)`, or the
general `forest(x1 + x2, basis = ~z)` - rewrites the formula and feeds the same
`forests = list(forest(), forest(basis = ))` channel `dbarts()` and
`dbartsSpec()` use ([[formulaTerms.R#ingestFormulaTerms]]; `bart2` dispatches
through it for every family except multinomial, which refuses a `forest()`
term by name, having no amplitude-coupled slot). A two-forest fit built this
way carries the same `bartcore.forests` control attribute as one built through
`forests =`, byte-identical at the same seed
([[test-formula-terms.R#"Block B"]]). What is absent is the NAMED causal verb,
not the capability: `bcf()` and `bartBCF` ship in **bartCause**, not dbarts
([[docs/design/bcf.md#Public creation surface]]). `bart2()`'s own `sd` and
`update.amplitude` knobs have no term spelling, so they do not reach a term's
forest the way `forests =` on `dbarts()` does.

[f9] `updateScale` re-derives the internal response transform. The latent
families have `fitScale() == 1` and `fitShift() == 0` by definition, so there
is no transform to re-anchor and the flag is ignored rather than refused.

[f10] Logistic weights are the observation counts its Polya-Gamma latents are
built from, so a swap is a model change rather than a reweighting:
[[MOD#LogisticResponse::setWeights]] redraws omega against the new counts
before returning, which is what makes the conduit coherent. The
positive-integer policy stated at creation
([[RIB#enforceBinaryWeightPolicy]]) holds on every mutation conduit too, and
`setData` hands the replacement counts through the same conduit, so a data swap
draws rather than cold-starts and replacement data given without weights is
single-trial. `setState` re-derives the latents through this same conduit when
the destination's weights differ from the state's recorded `weights.digest`
attribute, so a restore lands where a swap lands; on a match it re-derives
nothing and the round trip stays byte-identical.

[f11] `bart2(family = "multinomial")` builds its `dbartsSampler` directly:
`$fit` is the K-forest engine that ran, from one `bartcore_create`, with no
wrapper sampler interposed and no `$bc` element on the fit. Ten methods -
`setResponse`, `setOffset`, `setWeights`, `setSigma`, `setCalibration`,
`setForestWeights`, `setForestBasis`, `getFitsWithoutOffset`, `setModel` and
`setData` - are refused by name through the shared counts guard
([[bartcore.R#refuseCountsMutation]]), each naming the capability and the
channel that serves the caller where one exists. `setCounts`,
`setCategoryOffset` and `setCategoryTestOffset` are the three response
channels, public R5 methods
([[dbarts.R#setCounts, setCategoryOffset, setCategoryTestOffset]]).
`setPredictor`, `setCutPoints`, `setTestPredictor` and the global
`setActiveRows` stay open, which is why section 2's `setPredictor` and
test-surface cells read `S`.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction
([[dbarts.R#"is only available through bart2()"]]) and [[bart.R#bart2Hurdle]]
composes two ordinary `bart2()` fits - an occupancy probit and a lognormal
positive part, both built through [[bart.R#redirectCall]] - glued at report
time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f15] `setActiveRows` is a first-class 0/1 per-observation mask
(docs/plans/latent-subset-mask.md): one validating and normalizing scan owns it
([[CH#Chain::setActiveRows]]), each family composes it into its own precision
vector, and a rejection-drawn latent is SKIPPED - not drawn and discarded - at
an inactive row. The 0/1 domain is a refusal, not a convention: a fractional
value is a weighted likelihood and belongs to `setWeights`
([[dbarts.h#dbarts_sampler_setActiveRows]], pinned on the flat surface
alongside the all-ones no-op and the NULL clear at
[[test-capi.R#"capi_set_active_rows"]]). Every shipped family implements the
channel, multinomial included, so the bridge's active-row refusal
([[RIB#"active-row masking is not implemented for this response family"]]) is
family-generic, reached only by a future family that does not override the base
refusal. Per-family composition, Student-t's drawn-and-annihilated exception,
multinomial's global mask, the NaN pointwise log-likelihood, the surfaces and
the evidence are in [[docs/design/active-rows-mask.md#The contract]].

[f16] `prior.scale` names the per-forest prior ANCHOR - the forest-total prior
scale at k = 1, in response units - rather than an sd, and
`$getCalibration`/`$setCalibration` read and write it on every chain of a
single-forest sampler (1-based `forest` argument,
[[bartcore.R#resolveForestIndex]]). It is refused under a `k` hyperprior, a
sampled `k` having no single value to divide by, and likewise once the chains'
`k` have diverged; refused for BCF and multinomial forests at creation and
mid-chain alike ([f23]); and `prior.mean` is refused as not writable, naming
the `setOffset` recipe instead. The engine sites are
[[CH#Chain::resolvedNodeScale]] at creation and
[[CH#Chain::forestCalibration, Chain::setForestPriorScale]] mid-chain, sharing
one conversion so neither can drift from the other; the surface, the per-leaf-
model meaning and the rationale for each refusal are
[[docs/design/nameable-calibration.md#2. The surface]] and
[[docs/design/nameable-calibration.md#5. Refusals, and why each is one]].

[f17] Zero weights are accepted, not refused
([[A_class.R#"'weights' must all be non-negative"]] errors only below zero and
warns that zeros are ignored; bridge [[RIB#enforceBinaryWeightPolicy]]). The
conditionals are exact - leaf suffstats multiply by `w`
([[MOD#ConstantVarianceLeaf::accumulate]],
[[MOD#LinearGaussianLeaf::accumulateNodeStatistics]]), and the sigma posterior
counts only positive-weight rows ([[MOD#numPositiveWeights_]], recounted on
every install in [[MOD#GaussianResponse::installWeights]], consumed in
[[MOD#GaussianResponse::drawSigma]]). The empty-leaf veto counts POSITIVE-WEIGHT
members, so a leaf held alive only by zeroed rows is empty and its branch is
vetoed. The veto is a RANK, taken once over the branch's leaves in
[[MOV#logLikelihoodForBranch]] for EVERY leaf model, the branch-owning
constrained ones included, so model.hpp keeps only that model's own
feasibility sentinel. Occupancy elsewhere
still counts members deliberately, so this does NOT make zero-weight occupancy
match a compacted fit
([[docs/design/empty-leaf-veto.md#What counts as empty]]). The same law covers
the Student-t row and a GAUSSIAN K-forest. It says nothing about a probit or
logistic one, where a zero weight cannot exist to begin with: probit refuses
weights entirely and logistic holds them to positive integer counts, both at
creation, so the cell is family-dependent ([f48]).

[f18] For gaussian and heteroscedastic no latent vector exists: both leave
[[MOD#ResponseModel::latents]] at its nullptr default, and the bridge returns
`R_NilValue`. A K-forest sampler is not necessarily one of them:
[[CH#Chain::latents]] is a bare delegation to `response_->latents()` carrying
no coupling gate and no family switch, and [[RIB#bartcore_getLatents]] gates
only on that pointer being null, so a probit K-forest reports its truncated
normals ([[MOD#ProbitResponse::latents]]) and a logistic one its Polya-Gamma
omegas ([[MOD#LogisticResponse::latents]]) while a GAUSSIAN K-forest reports
none - which is why the bcf cell is family-dependent rather than plainly `S`.

[f19] A Student-t fit records `family = "gaussian"`, so `extract(type =
"loglik")` takes the gaussian branch of
[[generics.R#pointwiseLogLikelihood]], where an [[generics.R#isStudent]] check
on `resid.dist` scores the MARGINAL t density,
`dt((y - ev) / sd, df, log = TRUE) - log(sd)` ([[generics.R#dt]]) - the
observation-level likelihood loo/waic are defined on - rather than folding the
fit into the gaussian `dnorm` call ([[generics.R#dnorm]]). Pinned by value
against `dt(...)` at three indices ([[test-pointwise-loglik.R#"ll.t"]]).

[f20] Weights on logistic are PG copy counts and a zero count is refused by name
at creation ("drop zero-count rows", [[RIB#enforceBinaryWeightPolicy]]; R mirror
[[spec.R#enforceWeightPolicy]]), so zero-weight subsetting is foreclosed for
this family by the weight semantics themselves. The mid-chain `setActiveRows`
channel is what serves that caller instead ([f15]).

[f21] PER-FOREST masking is REFUSED, permanently and on model
grounds: the softmax margin is a log-sum-exp over the other K-1 forests,
so a row absent from category k's forest is still in every other
category's likelihood, and "row i is out of category k only" restricts no
likelihood at all. The refusal lands at the only reachable per-forest,
per-observation channel, [[RIB#bartcore_setForestWeights]], naming the
model reason rather than "unbuilt". BCF accepts a per-forest weight at
that same channel - a different (additive) coupling
where the per-forest mask is redundant with, not incoherent under, the
combined likelihood (see [f26]).

[f22] Multinomial's omegas live in the combiner, not the response model, and
[[MOD#MultinomialResponse]] does not override `latents()`, so
`getLatents` returns NULL. No accessor exposes them.

[f23] A named `prior.scale` is refused for BCF and multinomial forests both at
creation and mid-chain, by design: their per-forest leaf scales come from a
calibration map that owns them ([[COM#ForestSpec::amplitudePriorScale]]), so a
named value has nowhere to land. Three creation-time sites: `dbartsSpec()`'s
BCF composition ([[spec.R#"a named 'prior.scale'"]], an entry of the
`unsupported` vector beside the non-default-`k` entry
[[spec.R#"a non-default 'k'"]]), the bridge's own BCF-composition gate
([[RIB#refuseUnsupportedAmplitudeComposition]]) and the multinomial forest
builder ([[RIB#buildMultinomialSampler]]). Mid-chain there are two independent
sites rather than one shared gate: `$setCalibration` refuses a BCF sampler
R-side through [[bartcore.R#refuseAmplitudeMutation]] before reaching the
bridge, and a multinomial fit's `$fit` refuses the same call the same way for
the same softmax-calibration-map reason; underneath both,
[[CH#Chain::setForestPriorScale]] returns false whenever `combiner_ !=
nullptr`, surfaced at [[RIB#bartcore_setCalibration]], which is the gate a
direct low-level call on a multinomial forest hits
([[test-calibration-midchain.R#"softmax calibration map"]]). The R5 layer never
routes a BCF sampler there, so only the multinomial arm exercises the bridge
gate directly.

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and supports
`extract(type = "loglik")` directly: [[generics.R#hurdleLogLik]], reached from
[[generics.R#extract.bartHurdle]]'s `type == "loglik"` branch, combines the
occupancy's `log(1 - pi)` / `log(pi)` with the positive part's lognormal
density (a `-log(y)` Jacobian against the stored log-scale channel) at every
row. This is NOT the sum of the two components' own loglik channels - the
positive fit's covers only its y > 0 rows and carries no Jacobian - though each
component fit (`$occupancy` probit, `$positive` gaussian) still supports
`extract(type = "loglik")` independently.

[f26] Nothing on the active-rows path gates on the coupling: the shape probe
([[FAC#SamplerShape::supportsActiveRows]]) reports `supportsActiveRows`, the
mask composes into whatever precision the installed response owns
([[MOD#GaussianResponse::setActiveRows]] being the gaussian composition into
the case weights, which is what inherits the sigma df), and then into the
per-forest weights at [[CH#composeForestWeights]]. A K-forest chain's response
is a `ProbitResponse` or a `LogisticResponse` just as readily as a
`GaussianResponse` ([[CH#"switch (spec.family)"]]), each overriding
`setActiveRows` exactly as it does off a coupling, so [f15]'s coverage carries
the latent case - but every MEASUREMENT is gaussian two-forest
([[test-active-rows-pins.R#"masked.bcf"]], bitwise against `setWeights(w * a)`
in `train` and in `sigma`, plus the `bcf-equivalence.R` `masked` scenario; an
all-zeros mask runs finite and a fractional element is refused), and
no latent K-forest mask is measured anywhere. A per-forest mask is refused as
REDUNDANT rather than unbuilt, since the public
[[dbarts.R#setForestWeights]] / [[RIB#bartcore_setForestWeights]] channel
already expresses it - though that channel is deliberately NOT row removal
([[CH#Chain::setForestWeights]]: the row stays in occupancy, in the combination
and in the sigma df, while a zeroed composed weight does reach that forest's
empty-leaf veto). See
[[docs/design/bcf.md#The multiplier snap and the per-forest weight]].

[f27] A delegating decoration: both columns fall out of the base family's
implementation. The heteroscedastic [[CH#formMeanWeights]] reads
`response_->workingWeights()` - the COMPOSED `w * a` while a mask is installed
- and divides by `s^2(x_i)`, so a zero stays a zero. It is tested bitwise
against a composed-weight sampler ([[tests/cpp/test_sampler.cpp#testActiveRows]],
[[test-active-rows-pins.R#"heteroSampler"]],
[[test-active-rows-pins.R#"draws.hetero$varcount"]]), and its creation-time
calibration is one of the family and decoration paths a shared test measures
([[test-calibration-prior-draws.R#"heteroscedastic = anchorSampler"]]).

[f28] A heteroscedastic fit also records `family = "gaussian"`, but the same
gaussian branch of `pointwiseLogLikelihood` reads its `s.train` surface first:
[[generics.R#heteroscedasticScale]] supplies the per-observation scale whenever
the fit carries one, taking precedence over the scalar `object$sigma`
([[generics.R#chainFastest]]); the surface is stored on the fit in
[[bart.R#packageBartResults]]. Pinned by value as in [f19]
([[test-heteroscedastic-channels.R#"loglik scores at s(x)"]], tolerance 1e-12).

[f29] The MEAN forest's calibration reads no family flag:
[[CH#Chain::resolvedNodeScale]] runs at the
[[CH#"forest.leaf.scale = resolvedNodeScale"]] assignment before the
variance-forest branch ([[CH#buildVarianceForest]]), and the mid-chain
reader/writer ([[CH#Chain::forestCalibration, Chain::setForestPriorScale]])
reads the same `response_->fitScale()`/`fitShift()`. The variance forest itself
is a separate leaf model entirely, outside `forests_`, and is not addressable
by `setCalibration`: a heteroscedastic sampler's `shape.numForests` is 1, so
`forest = 2` is refused by the ORDINARY out-of-range check every single-forest
sampler hits ([[RIB#bartcore_getCalibration, bartcore_setCalibration]]), not by
a hetero-specific gate. Creation-time it is measured, as the tenth case of
`test-calibration-prior-draws.R`'s `anchorSamplers` sweep
([[test-calibration-prior-draws.R#"heteroscedastic = anchorSampler"]]), inside
the same band every other family in that loop is held to
([[test-calibration-prior-draws.R#"familyBand <- 0.09"]]); mid-chain it rides
the same generic mechanism as every other single-forest family.

[f30] Split verdicts, measured by benchmarks/R/composition-matrix.R.
`resid.dist = student()` with `variance =` is refused at
[[spec.R#"does not support Student-t residuals"]] - "the two are not yet shown
to compose" - a validation error only; the formal stays, and adjudicating the
composition of two scale mixtures on the same precision channel would reopen
it. [[docs/design/heteroscedastic.md#Student-t residuals refused]] records
the refusal and what lifting it would take.

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
([[spec.R#"a variance forest requires family"]]) before either component fits.
That refusal is deliberate, not a bug: it makes a heteroscedastic component
unreachable, so the positive fit is always homoscedastic
([[bart2.Rd#variance]], [[dbarts.Rd#variance]]).

[f35] `dart` is forwarded to both components ([[bart.R#redirectCall]]), each of
which is an ordinary single-forest chain that takes it.

[f36] The two initializations part company at more than one forest.

WARM START is refused at the forest count.
`refuseMultiForestWarmStart` ([[RIB#refuseMultiForestWarmStart]]) errors on
`numForests >= 2` at [[RIB#bartcore_installForests]]; the R5 method raises the
same rule in its own wording ahead of it
([[bartcore.R#refuseMultiForestWarmStart]], called from
[[dbarts.R#installTrees]]); and `bart2`'s `warm.start` refuses there too
([[bart.R#refuseMultiForestWarmStart]]), which is the entrance a `forest()`
term or a bases-carrying data object arrives by ([f7]). The cell is `M` and not
`R` for [f50]'s reason: the guard records that the install is not validated
above one forest, not that it is incoherent. It would ANSWER rather than raise
- [[SAM#Sampler::installForests]] reassembles the trees from a saved slot but
takes `hasAmplitudes`, `amplitudes` and `amplitudeVariances` off the donor's
LIVE chain state, so a slot-sourced install pairs one draw's forests with
another's glue - and nothing covers the result: no R test, no component pin
above one forest, no equivalence scenario, no calibration evidence. The same
count catches the K-forest softmax.

GROW-FROM-ROOT is shipped and covered at two forests, so it carries no guard.
[[CH#growForestFromRoot]] composes through the combiner inside its own sweep -
`drawForestGlue`, `formForestResponse`, then `drawGlue` and `afterCombine`, the
same calls [[CH#Chain::run]] makes - with a variance-forest pre-step above the
loop. The two-forest branch is exercised from R by
[[test-prior-init-composed-law.R#"grow-from-root holds the same law"]], which
pins the per-forest no-empty-leaf conditioning against a composed per-forest
weight, and pinned in the component suite by
[[tests/cpp/test_sampler.cpp#testBCFGrowForestFromRoot]] (a hardcoded
characteristic value of `a mu + b_z tau`) with the K-forest twin
[[tests/cpp/test_sampler.cpp#testMultinomialGrowForestFromRoot]] beside it, and
walked by the mutation fuzzer. What it still lacks is an equivalence scenario
and calibration evidence, which is the standing multi-forest evidence gap
([f48]) rather than anything specific to the initializer.

[f38] The MEAN forest keeps DART; the variance forest never takes it
([[CH#buildVarianceForest]] never sets `useDart`, whose default is false in
[[CH#SamplerOptions]]).

[f39] The equivalence gate replays a fixed scenario set against a recorded
baseline and requires BITWISE-identical draws
([[docs/plans/README.md#RNG classes and their gates]]). Current baselines:
`equivalence-1e5f80b2.rds` (50 scenarios),
`bcf-equivalence-3c81d6df.rds` (12) and
`multinomial-equivalence-4d9a3337.rds` (11) - benchmarks/baselines/MANIFEST.
The names in this column are the keys in
[[benchmarks/R/equivalence.R#makeScenarios]]; each row lists only the scenarios
whose family it is.

[f40] SBC is simulation-based calibration: draw parameters from the prior,
simulate data, refit, and check that each true value's rank within its
posterior draws is uniform. A fraction here counts FUNCTIONALS - how many of
the family's checked functionals land inside the band - and the verdicts are
recorded in docs/plans/sbc-family-tiers.md (status BUILT) and
docs/plans/sbc-calibration.md (DONE). The A/B/C "tiers" in the latter are
FEATURE tiers (A baselines/DART/weighted/BCF, B linear leaf, C GP
leaf), not family tiers: there is no per-family tier ladder, only the
per-family verdicts this column carries.

[f41] gamma3 flagged in one stream and RESOLVED as the cutpoint-vs-mean-level
ridge mixing slowly: it does not reproduce across streams and sits at 0.31 of
the band at 3x the chain length. The cutpoint block, the latent eta and all K
category probabilities calibrate.

[f42] `avg.mu` - the identified mean - passes cleanly; `r` and `agg.psi` flag at
thin = 30 and cross into the band at 5x the spacing, read as slow mixing along
the r-psi ridge. Measured at two thinning settings rather than three; a third,
longer run is still owed.

[f43] Aggregate `p_k(x*)` and the three raw per-forest `f_ik` cells all pass,
the functionals at band 0.1282 and the three cells at 0.0688/0.0824/0.0675
(the acceptance run `Rscript benchmarks/R/sbc.R multinom 200 150 30`;
[[benchmarks/R/sbc.R#cellNames]] is the function that ranks the cells).
[[COM#MultinomialForestCombiner::afterCombine]] draws the level from its exact
leaf-space conditional.

[f45] Out of the SBC matrix by scope, each for its own recorded reason
([[docs/plans/sbc-family-tiers.md#Decision - scope]]): aft because its
censoring status is fixed at
creation, so a prior-draw replication cannot vary it (the enabler is a status
setter); hazard and hurdle because their person-period / two-part designs depend
on `y0`, which breaks exchangeability, and because neither owns any sampling
code.

[f46] Tier A PASS, with the sigma channel resolved as slow mixing along the
(a, mu) ridge ([[docs/plans/sbc-calibration.md#Final summary]]). Explicitly out
of the family-tiers matrix.

[f47] OUT but DEFERRED rather than blocked: prior draws never reach
`varianceForest_` today, and the capability is liftable R-side through `setState`
([[docs/plans/sbc-family-tiers.md#Decision - scope]]).

[f48] The amplitude coupling's family reach: gaussian, probit and logistic
build; aft, ordinal and nbinom are refused at all three creation routes, each
naming what it is missing ([[spec.R#"a treatment forest does not support family"]],
[[RIB#refusedAmplitudeFamilyReason]], [[FAC#createAmplitudeSampler]]). The
calibration map's anchor is family-keyed ([[CH#latentScaleAnchor]]), and the
model, the map and the reporting surfaces are
[[docs/design/multiplier-combiner.md#The calibration map, general in K]] and
[[docs/design/multiplier-combiner.md#Surfaces]].

Which cells vary by base family. `setResponse` and `setOffset` are OPEN under
every family, [[CH#Chain::setResponse]] handing the response `combinedFits()`
rather than forest 0's bare totals. `setWeights` is REFUSED under probit and
open under gaussian and logistic, whose counts are its Polya-Gamma shape;
`setSigma` is open under gaussian and REFUSED under probit and logistic - both
through the ORDINARY single-forest guards ([[RIB#refuseBinaryWeightChange]],
[[RIB#refusePinnedSigmaChange]]), the sampler answering `shape.family` for
itself. The same [[RIB#enforceBinaryWeightPolicy]] that refuses a probit weight
outright and holds a logistic one to positive integer counts makes the
zero-weight-subset cell family-dependent, and [f18] makes `getLatents` one.
`updateScale = TRUE` stays REFUSED under EVERY family:
[[bartcore.R#refuseAmplitudeMutation]] keys on the sampler carrying bases and
[[RIB#refuseMultiForestResponseMutation]] on `numForests >= 2`, never on
family, so a probit K-forest is refused too although its transform is the
identity. [[test-bcf-family.R#"binaryFits"]] is the whole of the latent
K-forest's evidence: no equivalence scenario and no SBC coverage reaches one.

[f49] The test-surface cell is `R` because the blend
`sum_f dot(a_f, B_f(i, .)) f_f(x_i)` needs an off-sample basis the sampler does
not have: the resident test store, `run()$yhat.test` and the SAMPLER's own
`predict()` all refuse through [[RIB#refuseUndefinedTestFits]], gated on
`testFitsAreDefined` rather than on a forest count
([[docs/design/multiplier-combiner.md#What this family does not do]]). What
ships instead reports the per-forest pieces and leaves the contraction to the
caller, whose bases they are: the engine replay
([[RIB#bartcore_predictPerForest]],
[[CH#Chain::predictPerForestFromSavedSample]]), which replays every forest at
the new rows RAW - on the forests' internal scale, with no `fitScale`, no
`fitShift` and no offset, because the combination carries all three - and, at
the FIT level, an R contraction ([[generics.R#predictBlend]]) taking the bases from a `bases =`
argument or, for a `bart2` `forest()` term, re-evaluating the declaring formula
against `newdata` under the fit's own factor levels. The in-sample channel has
the same shape: [[bart.R#packageBartResults]] packages `forestFits` and `glue`
onto any bcf-shaped fit and `extract(type = "forest", forest =, contribution = )`
serves it, refused by name on a fit without forest reporting or with
`sample = "test"` ([[docs/design/multiplier-combiner.md#Surfaces]]). Both
amplitude predict routes require `keeptrees`/`keepTrees`. Evidence:
inst/tinytest/test-predict-forest.R, inst/tinytest/test-predict-blend.R and
[[tests/cpp/test_sampler.cpp#testAmplitudePerForestReplay]].

[f50] [[bart.R#checkFamilyUnsupportedArgs]] raises a bare `family = "X" does
not support 'warm.start' or 'n.grow.sweeps'` for the four alternate-family
`bart2` arcs, with no model reason stated there or anywhere else: the guard
records that the arc was never built for either argument, not that either is
incoherent under the family, so these cells are `M` and not `R`. The two
columns are separate although one check refuses both today: a warm start
installs a previous fit's trees and family state, while grow-from-root
needs the family's working response inside the root-growth recursion, and
either could arrive for a family without the other. The warm-start column is
covered twice over for multinomial: its K forests also meet the forest-count
guard [f36] states, which closes the R5 sampler surface this `bart2` check
does not reach. The grow-from-root column is not - `$growFromRoot` runs on a
K-forest sampler, [f36]'s second half - so that cell records the `bart2` arc
alone.

## Gaps

Every MISSING (`M`) and UNVERIFIED (`?`) cell above, as a candidate work item,
grouped by family, together with the standing evidence gaps no cell records.
Scheduling is VD's and the orchestrator's; nothing here carries a schedule, and
REFUSED cells are deliberately absent - they are part of the models.

**gaussian.** None; every column is S or an intentional `-`.

**student (Gaussian + Student-t residuals).** No `xbart()` surface
([[xbart.R#xbart]]). Only one dedicated tinytest file.

**probit.** None.

**ordinal.** `warm.start` / `n.grow.sweeps` are unbuilt for the arc ([f50]).
One dedicated tinytest
file, in which pointwise loglik is exercised only for shape. SBC gamma3
resolved but not re-run at full R.

**nbinom.** No `xbart()` token, and warm start / grow-from-root are unbuilt
([f50]). SBC
`r`/`agg.psi` flag standing, read as slow mixing along the r-psi ridge, with a
third, longer run owed ([f42]). Real-valued dispersion remains a recorded open
item (TODO's `negbin-real-dispersion` item).

**multinomial.** No flat-C creation path at all ([f4]). No `getLatents`
([f22]), and no warm start or grow-from-root ([f50]).
The ENGINE's own per-observation log-likelihood channel stays undefined for
this family; `extract(type = "loglik")` scores the multinomial log-pmf on the
REPORTED probabilities instead ([[generics.R#multinomialLogLik]]).

**aft.** No `xbart()` token. Out of the SBC matrix because the censoring
status is fixed at creation; a status setter is the named enabler ([f45]). The
separate exact oracle benchmarks/R/aft-exact.R is not a MANIFEST entry.

**hazard.** No flat-C token and no engine code of its own ([f6]); no
`xbart()` reach. Out of SBC by design ([f45]). One dedicated tinytest file.

**hurdle.** No `xbart()` or flat-C reach, and no warm start or
grow-from-root ([f50]).

**bcf.** No `xbart()` reach. The
NAMED `bcf()` causal verb is bartCause's, not dbarts's; the K-forest amplitude
capability itself is reachable from `bart2()`'s formula interface ([f7]).
A donor warm start is refused at the forest count, untested above
one ([f36]); grow-from-root runs there and is covered, R-side and in the
component suite, but by no equivalence scenario ([f36]). Whole-data `setData` stays undesigned (open question 1 of the
model-space survey, docs/design/model-space-survey.md). The probit and logistic
cases have no equivalence scenario, no SBC coverage and no measured active-rows
mask ([f48], [f26]).

**heteroscedastic.** No `xbart()` reach. Out of the SBC matrix, deferred
not blocked ([f47]).
