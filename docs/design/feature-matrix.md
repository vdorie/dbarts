# Response-model feature matrix

Status: living reference, current at a262cd26 (2026-08-13). Carries no landing
date and is not a design proposal - the orchestrator updates it in place at
every landing that changes a cell, and VD uses it to schedule feature
completion.

What each shipped response model can and cannot do, one row per model and one
column per capability that bears on scheduling. Every SHIPPED and REFUSED cell
carries an anchor verified against the live tree; a wrong cell misdirects
scheduling, so a cell that cannot be verified is marked `?` rather than
guessed.

## Legend

| code | meaning |
|---|---|
| `S` | SHIPPED. Works today; the anchor is the site that makes it work. |
| `R` | REFUSED on model or identification grounds - the refusal is part of the model, not a hole. The anchor is the refusal site. |
| `P` | PLANNED. A named design arc covers it; the slice is named in the cell. |
| `M` | MISSING. Not built, no schedule. An anchor, when given, is a guard that errors *because the thing is unbuilt* - a recorded door, not a model refusal. |
| `-` | N/A. The concept does not apply to this row; the row footnote says why. |
| `?` | UNVERIFIED. Constructs today with no refusal site, but no test, doc or adjudication backs it. Do not schedule against the cell until it is settled; every `?` is listed under "Gaps". |

Path aliases used in anchors (line numbers are at a262cd26):

    RIB   src/R_interface_bartcore.cpp      CAPI  inst/include/dbarts/dbarts.h
    MOD   src/bartcore/model.hpp            CH    src/bartcore/chain.hpp
    FAC   src/bartcore/facade.hpp           COM   src/bartcore/combiner.hpp
    MOV   src/bartcore/moves.hpp            SAM   src/bartcore/sampler.hpp
    bart.R, dbarts.R, spec.R, rbart.R, xbart.R, data.R, generics.R,
    A_class.R, bartcore.R      -> R/<name>
    sampler.Rd -> man/dbartsSampler-class.Rd     bart.Rd -> man/bart.Rd

## Rows

Thirteen rows. The first ten are response models proper; the last three are
couplings or decorations over a base response that a user selects the same way
and schedules against the same way, so they earn rows.

| key | model |
|---|---|
| gaussian | Gaussian (`ResponseFamily::gaussian`, `GaussianResponse` MOD:2721) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, `TResponse` MOD:3960) |
| probit | Binary probit (`ProbitResponse` MOD:2997) |
| logistic | Binary logistic, weights = observation counts (`LogisticResponse` MOD:3456) |
| ordinal | Ordered categorical, cumulative probit (`OrdinalResponse` MOD:3136) |
| nbinom | Negative binomial, positive-integer dispersion (`NBResponse` MOD:4257) |
| multinom | Multinomial softmax, K forests (`MultinomialResponse` MOD:3635 + combiner) |
| aft | AFT survival, log-normal (`AFTResponse` MOD:3713) |
| hazard | Discrete-time hazard (person-period sugar, dbarts.R:437-486) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, bart.R:1958) |
| bcf | Bayesian causal forest, two forests (`BCFForestCombiner` COM:518) |
| grouped | Grouped random intercepts (`GroupedResponse` MOD:4624) |
| hetero | Heteroscedastic variance forest (CH:679) |

The engine's `ResponseFamily` enum has only six tokens (MOD:2524: gaussian,
probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf, grouped
and hetero are all reached some other way, which is exactly why they need rows
here rather than an enum read. Leaf models (constant, monotone, linear, GP) are
an orthogonal axis and are not rows; where a leaf model gates a capability the
cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `rbart_vi()` | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|---|
| gaussian | S bart.R:2398 | S bart.R:559 | S dbarts.R:359 | S rbart.R:49 | S xbart.R:27 | S CAPI:341 |
| student | S bart.R:2301 | S bart.R:573 | S dbarts.R:340 | M rbart.R:60 | M xbart.R:9-38 | S RIB:2464 [f2] |
| probit | S bart.Rd:151 | S bart.R:560 | S dbarts.R:360 | S data.R:523 | S xbart.R:27 | S CAPI:340 |
| logistic | - [f1] | S bart.R:561 | S dbarts.R:361 | R rbart.R:49 | S xbart.R:27 | S RIB:1534 |
| ordinal | R bart.R:2402 | S bart.R:564 | S dbarts.R:363 | R data.R:468 | R data.R:468 | S RIB:1542 [f3] |
| nbinom | - [f1] | S bart.R:565 | S dbarts.R:364 | M rbart.R:49 | M xbart.R:27 | S RIB:1549 [f3] |
| multinom | - [f1] | S bart.R:563 | R dbarts.R:357 | R data.R:468 | R data.R:468 | M [f4] |
| aft | ? bart.R:2398 [f5] | S bart.R:562 | S dbarts.R:362 | S rbart.R:49 | M xbart.R:27 | S CAPI:343 |
| hazard | - [f1] | S bart.R:566 | S dbarts.R:365 | M rbart.R:49 | M xbart.R:27 | M [f6] |
| hurdle | - [f1] | S bart.R:569 | R dbarts.R:395 | M | M | M [f6] |
| bcf | - [f1] | R bart.R:625 [f7] | S dbarts.R:350 | R rbart.R:60 | M | S CAPI:348 |
| grouped | - [f1] | M [f8] | M [f8] | S rbart.R:367 | M | S RIB:1911 [f3] |
| hetero | - [f1] | S bart.R:549 | S dbarts.R:346 | R rbart.R:60 | M | S RIB:2001 [f3] |

`dbartsSpec()` (spec.R:492-500) resolves the seven single-forest tokens - auto,
gaussian, probit, logistic, aft, ordinal, nbinom - plus BCF through its
`treatment =` argument (spec.R:487) and a variance forest through `variance =`
(spec.R:483); it does not reach multinomial, hazard, hurdle or grouped.

## 2. Mutation channels on the R5 `dbartsSampler`

The channels that make dbarts a conditional model inside an outer sampler.
`updateScale` is broken out because it is refused independently of the setter
it rides on.

| model | `setResponse` | `setOffset` | `updateScale = TRUE` | `setPredictor` (+ per-obs) | `setWeights` | `setSigma` | test surface |
|---|---|---|---|---|---|---|---|
| gaussian | S MOD:2754 | S MOD:2834 | S MOD:2834 | S RIB:4707, 4863 | S MOD:2803 | S RIB:2674 | S RIB:4402, 4455 |
| student | S MOD:4029 | S MOD:4036 | S MOD:4036 | S RIB:4707 | S MOD:4053 | S RIB:2674 | S RIB:4402 |
| probit | S MOD:3050 | S MOD:3056 | - [f9] | S RIB:4707 | R RIB:1610 | R RIB:2676 | S RIB:4402 |
| logistic | S MOD:3524 | S MOD:3530 | - [f9] | S RIB:4707 | R RIB:1603 [f10] | R RIB:2676 | S RIB:4402 |
| ordinal | S MOD:3187 | S MOD:3195 | - [f9] | S RIB:4707 | R RIB:1610 | R RIB:2676 | S RIB:4402 |
| nbinom | S MOD:4324 | S MOD:4331 | - [f9] | S RIB:4707 | R RIB:1610 | R RIB:2676 | S RIB:4402 |
| multinom | - dbarts.R:776 [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] |
| aft | S MOD:3788 | S MOD:3801 | S MOD:3801 | S RIB:4707 | R RIB:1610 | S RIB:2674 | S RIB:4402 |
| hazard | S MOD:3050 [f6] | S MOD:3056 | - [f9] | S RIB:4707 | R RIB:1610 | R RIB:2676 | S RIB:4402 |
| hurdle | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] |
| bcf | S COM:733 | S COM:733 | R bartcore.R:285 | S RIB:4707, 4863 | S RIB:4549 | S RIB:4325 | R RIB:2650 |
| grouped | R RIB:4309 [f13] | S MOD:4718 | S MOD:4718 | S RIB:4707 | S MOD:4732 [f14] | S RIB:2674 [f14] | S RIB:4402 |
| hetero | S RIB:2581 | S RIB:2581 | R RIB:2581 | S RIB:4707, 4863 | S RIB:2577 | R RIB:2672 | S RIB:4402 |

`setData` (whole-data replacement, n free) is single-forest and dense-store
only (RIB:3365) and is refused for grouped (RIB:4339) and aft (RIB:4341);
BCF/multinomial whole-data `setData` stays undesigned by the model-space
survey's verdict (model-space-survey.md doors 1 and 3).

## 3. Row subsetting, latents, calibration

| model | zero-weight row subset | active-rows mask [f15] | `getLatents` | pointwise loglik | nameable calibration [f16] |
|---|---|---|---|---|---|
| gaussian | S sampler.Rd:143, MOD:2729 [f17] | S MOD:2814 | - RIB:5599 [f18] | S generics.R:48 | S dbarts.R:1359, 1364 [f16] |
| student | S MOD:4009-4016 [f17] | S MOD:4064 | S MOD:4076 | ? generics.R:48 [f19] | S dbarts.R:1359, 1364 [f16] |
| probit | R RIB:1581 | S MOD:3040 | S MOD:3076 | S generics.R:59 | S dbarts.R:1359, 1364 [f16] |
| logistic | R RIB:1585 [f20] | S MOD:3476 | S MOD:3550 | S generics.R:59 | S dbarts.R:1359, 1364 [f16] |
| ordinal | R RIB:2455 | S MOD:3172 | S MOD:3217 | M generics.R:86 | S dbarts.R:1359, 1364 [f16] |
| nbinom | R RIB:2461 | S MOD:4304 | S MOD:4357 | M generics.R:86 | S dbarts.R:1359, 1364 [f16] |
| multinom | R RIB:2990 | S COM:936 [f21] | M MOD:3635 [f22] | M generics.R:86 | R [f23] |
| aft | R RIB:2451 | S MOD:3770 | S MOD:3831 | S generics.R:64 | S dbarts.R:1359, 1364 [f16] |
| hazard | R RIB:1581 [f6] | S MOD:3040 [f6] | S MOD:3076 | S generics.R:59 [f24] | S dbarts.R:1359, 1364 [f6] |
| hurdle | R bart.R:944 | - [f12] | - [f12] | M generics.R:86 [f25] | - [f12] |
| bcf | S COM:722-733 [f17] | S MOD:2814, CH:1160 [f26] | - [f18] | M generics.R:86 | R [f23] |
| grouped | S MOD:4587-4609 | S MOD:4743 [f27] | S MOD:4752 | S generics.R:1396 | S MOD:4774 [f27] |
| hetero | S CH:3615, MOD:306 | S CH:3611 [f27] | - [f18] | ? generics.R:48 [f28] | S test-calibration-prior-draws.R:251 [f29] |

## 4. Model composition

| model | variance forest | grouped ranef | DART | warm start | grow-from-root |
|---|---|---|---|---|---|
| gaussian | S FAC:680 | S CH:581 | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| student | ? [f30] | S CH:581 | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| probit | R spec.R:338 | S CH:581 | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| logistic | R spec.R:338 | S CH:581 | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| ordinal | R spec.R:338 | M RIB:2789 [f31] | S CH:525 | R bart.R:493 | R bart.R:493 |
| nbinom | R spec.R:338 | M RIB:2794 [f31] | S CH:525 | R bart.R:493 | R bart.R:493 |
| multinom | R FAC:771 | M RIB:1908 [f32] | M CH:4524 [f33] | R bart.R:493 | R bart.R:493 |
| aft | R spec.R:338 | S CH:581 | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| hazard | R spec.R:338 | M rbart.R:49 [f6] | S CH:525 | S bart.R:1014 | S dbarts.R:841 |
| hurdle | ? [f34] | M | S bart.R:922 [f35] | R bart.R:493 | R bart.R:493 |
| bcf | R FAC:758 | R RIB:2244 | R spec.R:393 | S SAM:717 [f36] | S CH:1589 [f36] |
| grouped | ? [f30] | - | S rbart.R:578 | M rbart.R:45-51 [f37] | M rbart.R:45-51 [f37] |
| hetero | - | ? [f30] | S CH:525 [f38] | S SAM:722 | S CH:1577 |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused at dbarts.R:844-853 and no-op at CH:1567, so every family above
reads "constant leaf" in that column.

## 5. Evidence

| model | equivalence baseline [f39] | SBC verdict [f40] | dedicated tinytest files |
|---|---|---|---|
| gaussian | 21 scenarios (friedman, weighted, splitprobs, chains, setdata, wtoffset, quants, categorical, missing, dart, linear, gp, zeroweights, sparse, wtgp, chik2, pred*) | PASS 7/7 | ~20 (test-sampler-*.R, test-bart-bart2.R, test-zero-weights.R) |
| student | `student` | PASS 4/4 | test-robust-errors.R only |
| probit | `probit`, `chik` | PASS | test-binaryResponse-hyperprior.R, test-family.R, test-weighted-binary-ppd.R |
| logistic | `logistic`, `wtlogistic` | PASS 6/6 | test-weighted-logistic.R, test-family.R |
| ordinal | `ordinal` | 9/10 [f41] | test-ordinal.R only |
| nbinom | `nbinom` | 1/3 [f42] | test-nbinom.R only |
| multinom | 10 scenarios, own harness | aggregate PASS, raw `f_ik` OPEN [f43] | 5 (test-multinomial-*.R) |
| aft | `grouped_aft` only [f44] | OUT [f45] | test-aft.R, test-rbart-aft.R |
| hazard | `hazard` | OUT [f45] | test-hazard.R only |
| hurdle | `hurdle` | OUT [f45] | test-hurdle.R, test-hurdle-surface.R |
| bcf | 11 scenarios, own harness | PASS [f46] | 6 (test-bcf*.R) |
| grouped | `grouped`, `grouped_aft` | PASS (tier A) | 14 (test-rbart-*.R) |
| hetero | `hetforce`, `hetswap`, `hetpartial` | OUT [f47] | test-heteroscedastic.R, -mutation.R |

Flat-C test coverage is thinner than family reach: inst/tinytest/test-capi.R
drives only the `""`/`"probit"` tokens plus grouped (:437), heteroscedastic
(:189) and BCF (:525) by control attribute. logistic, aft, ordinal and nbinom
are reachable through dbarts.h and untested there.

## Footnotes

[f1] `bart()` carries no `family` formal (bart.R:2268-2303); the family is
inferred from the response (numeric -> gaussian, 2-level -> probit) and
`resid.dist` is its only response-law lever. Families needing an explicit token
are out of reach by signature, not by refusal. Ordinal is the exception: it is
explicitly backstopped at bart.R:2402.

[f2] Student-t is not a token anywhere. It is selected by a finite `resid.df`
attribute on the model SEXP (RIB:2464-2475, gaussian-only gate) and refused for
every non-gaussian family R-side at spec.R:274-280. The engine family stays
`gaussian`, which is why the whole gaussian row applies to it in tables 2-4.

[f3] Reachable through `dbarts_sampler_create` but not discoverable from the
shipped header: ordinal needs the `bartcore.n.categories` control attribute
(RIB:441), nbinom `bartcore.dispersion` (RIB:454), grouped `bartcore.groups`
(RIB:1911), heteroscedastic `bartcore.variance` (RIB:2001). The header's
`family` documentation (CAPI:338-357) names only probit, logistic, gaussian,
aft and BCF.

[f4] Multinomial has no dbarts.h creation path at all; it is reached only by
the R-internal `.Call` entries `C_dbarts_bartcore_createMultinomial` /
`...Counts` (bartcore.R:798, 861).

[f5] A `Surv()` response passed to `bart()` auto-dispatches through `dbarts()`
and nothing refuses it (bart.R:2398), but man/bart.Rd:151 does not document it.
Reachable, undocumented, intent unadjudicated.

[f6] `family = "hazard"` / `"hazard.probit"` / `"hazard.logistic"` is
person-period ingestion sugar: dbarts.R:437-486 expands the design and remaps
the token to `"probit"` or `"logistic"` at dbarts.R:483 before any model is
built. The resulting sampler *is* an ordinary binary one, so its whole row
equals the probit (or logistic) row, and the fit records `family = "probit"`.
No engine code, hence no C-API token and no SBC arm. Refusals inside the
expander: no formula interface (:438), no `subset` (:453), no `test` (:456).

[f7] `treatment` is not a `bart2()` formal, so the unknown-argument check at
bart.R:625-635 kills it. BCF's public creation surface is `dbarts()` /
`dbartsSpec()`; `bcf()` and `bartBCF` ship in **bartCause**, not dbarts
(bcf-public-surface.md:483-496, fork 4 RESOLVED VD 2026-08-11; echoed
bcf.md:344-352). This row therefore tracks the engine capability, not a dbarts
user-facing verb.

[f8] `bartcore.groups` is written at exactly one site, rbart.R:367, and no other
entry point carries a `group.by` formal, so grouped random effects are an
`rbart_vi()`-only surface.

[f9] `updateScale` re-derives the internal response transform. The latent
families have `fitScale() == 1` and `fitShift() == 0` by definition, so there
is no transform to re-anchor and the flag is ignored rather than refused.

[f10] Logistic weights are the observation counts its Polya-Gamma latents were
built from, so they cannot change after creation (RIB:1603). This is the one
weight refusal that is "unbuilt" rather than "incoherent": creation accepts
positive-integer weights at RIB:1585.

[f11] `bart2(family = "multinomial")` returns a `dbartsSampler` that is only a
HOST SHELL - it carries the design and priors but not the model, which lives on
the separate bartcore handle `$bc` (bart.R:1258, 1336). Every R5 mutator is
refused by `refuseHostMutation` (dbarts.R:770-783). At the handle level:
`setResponse` R RIB:2547 and `setOffset` R RIB:2544 (both redirect to the counts
channel), `setWeights` R RIB:2554, `setSigma` R RIB:2676, `setTestOffset` R
RIB:2319 (a flat offset leaves the simplex), `setPredictor` S RIB:4707,
`setTestPredictor` S RIB:4402. Its own channels - `setCounts` RIB:3518,
`setCategoryOffset` RIB:3596, `setCategoryTestOffset` RIB:3635 - are S at the
bridge but unexported (`dbarts:::`) and absent from the R5 object.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction at
dbarts.R:395 and `bart2Hurdle` (bart.R:1958) composes two ordinary `bart2()`
fits - an occupancy probit (bart.R:1990) and a lognormal positive part
(bart.R:1999) - glued at report time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f13] Grouped samplers fix the response at creation (RIB:4309). The engine
method exists and delegates faithfully (MOD:4707), so this is an unbuilt door,
not a model refusal - but it errors, so callers see a refusal.

[f14] Reads off the BASE family: grouped gaussian takes `setWeights` (MOD:4732)
and `setSigma`; grouped probit is refused on both (RIB:1610, RIB:2676); grouped
aft takes `setSigma` and refuses `setWeights`.

[f15] Arc `latent-subset-mask` (TODO:155), design FINAL, ARC COMPLETE (S0
through S4 LANDED); artifacts .claude/latent-subset-mask-design/. A first-class
0/1 `setActiveRows` channel each family composes into its own precision vector,
with the latent draw skipped for inactive rows. Slices: **S0** pins (no engine
change); **S1** the channel plus gaussian, Student-t, probit, ordinal; **S2**
logistic, nbinom, aft; **S3** multinomial (global only); **S4** surface,
records, baselines. S0 landed at dc11a805 (the pins, now
inst/tinytest/test-active-rows-pins.R). S1 landed at 6db22aee: the engine
channel - `Chain::setActiveRows` CH:1341, which owns the single validating and
normalizing scan, `Sampler` SAM:1204, the facade's pure virtual FAC:298 and its
shape probe FAC:86 - plus gaussian, Student-t, probit and ordinal, the R5
`$setActiveRows` (dbarts.R:1109) and the bridge entry (RIB:3738). S2 landed at
87d370ea: logistic (`workingWeights()` MOD:3476) and nbinom
(`workingWeights()` MOD:4277) serve a SEPARATE a_i omega_i composite rather
than writing the zero into omega_ itself, since the working response divides
by it and 0 * inf in the node kernels is a NaN; nbinom's `setActiveRows`
(MOD:4304) additionally restricts the collapsed statistic S the dispersion
grid draw reads and REBUILDS the count-histogram kernel behind L_k
(`NBDispersionPrior::computeKernel` MOD:4199) over the active rows at every
mask change - the channel's one per-install cost; aft's `setActiveRows`
(MOD:3770) composes into its contained Gaussian, inheriting the sigma
degrees-of-freedom recount, and skips the censored redraw at an inactive row
(MOD:3754). All three report NaN pointwise log-likelihood at an inactive row.
Oracles: per-family kernel comparisons against the compacted arm, bitwise in
value and in RNG stream (`testActiveRowsLogisticKernel`
tests/cpp/test_model.cpp:5308, `testActiveRowsNBKernels` :5384,
`testActiveRowsAFTCensored` :5466 - each latent being a rejection sampler
means a discard-rather-than-skip at an inactive row fails the arm outright),
plus a sampler-level conditional independence oracle under substituted
inactive responses (inst/tinytest/test-active-rows-pins.R's S2 block:
substituting arbitrary in-support values at the inactive rows leaves every
active row's recorded draw bitwise). The heteroscedastic and grouped pins
move from FINITENESS to bitwise `setWeights(w * a)`
(test-active-rows-pins.R's hetero block; tests/cpp/test_sampler.cpp's grouped
block in `testActiveRows`), and ordinal - already S1 - gains a sampler-level
independence arm of its own beside the kernel-level coverage. S3 landed at
8b047f8b: multinomial's mask is GLOBAL and lands on the softmax coupling
rather than the response, which holds no precisions of its own -
`MultinomialResponse::setActiveRows` (MOD:3671) is a pass-through that only
advertises the capability (MOD:3670), and `Chain::setActiveRows` forwards the
mask to `MultinomialForestCombiner::setActiveRows` (COM:936) after the
response's own install (`ForestCombiner::setActiveRows` COM:486 is the inert
default every additive coupling relies on instead). An inactive row's K
interleaved Polya-Gamma draws are SKIPPED, not drawn and discarded, in
`drawForestGlue` (COM:1007), and its composed precision is zeroed in every
category in `formForestResponse` (COM:1064-1065); the row keeps its leaf
occupancy and its reported softmax probabilities, and omega is never zeroed
since the working response divides by it. PER-FOREST masking is refused
permanently on model grounds at the only reachable per-forest,
per-observation channel, `bartcore_setForestWeights` (RIB:3700-3704) - see
[f21] for the full statement. The bridge's active-row refusal (RIB:3748) no
longer names multinomial: the old per-family `activeRowsFamilyName` helper is
deleted, and the message is now family-generic, reached only by a future
family that does not override the base refusal. Oracles: the kernel-level
`testActiveRowsMultinomialKernel` (tests/cpp/test_sampler.cpp:4089) pins the
SKIP semantics and the zeroed composed precision bitwise against a compacted
combiner, in both value and Polya-Gamma stream; the shape probe flips in
`testMultinomial` (tests/cpp/test_shape.cpp:335); and
inst/tinytest/test-active-rows-pins.R's S3 block adds the same
substituted-response independence arm the S2 families got (moving successes
AND trial counts at the inactive rows, since PG(n_i, .) sums n_i variates),
plus an all-zeros mask run (every category forest at its prior, every row
still reporting a simplex) and the `setForestWeights` model-grounds refusal
(also pinned in test-forest-weights.R). S4 landed at 93afd635
(implemented as 76fd3ba6, amended during independent review): Rd
(man/dbartsSampler-class.Rd), NEWS (inst/NEWS.Rd), a named recipe
(man/bart.Rd), the dbarts.h reservation (docs/plans/c-api-growth.md),
two new equivalence.R scenarios (maskprobit, maskordinal) and one
bcf-equivalence.R scenario (masked, pinning BCF - see [f26]). ARC
COMPLETE; what remains is the flat-C entry, which rides the dbarts.h
reshape's S1.

[f16] Arc `nameable-calibration` (TODO:279), design AMENDED FINAL, LANDED
through S2; artifacts .claude/nameable-calibration-design/. Names the
per-forest prior ANCHOR (`prior.scale`, the forest-total prior scale at k = 1,
in response units) rather than an sd, with a `$getCalibration` /
`$setCalibration` pair. Slices: **S0** signature freeze, LANDED 4c866286;
**S1** creation half, LANDED c2a7e89b; **S2** mid-chain get/set, LANDED
d809b944 (+ a records correction, 7da36dc3); **S3** the flat-C half, executed
inside the dbarts.h reshape's S1. The R surface is now COMPLETE:
`$getCalibration`/`$setCalibration` read and write every chain of any
single-forest sampler, with a 1-based `forest` arg (`resolveForestIndex`,
bartcore.R:1051) mapped onto the engine's 0-based one. S1 names the model's
`prior.scale` slot (A_class.R:398), resolved from `node.prior`'s `scale =` /
`sd =` spelling at `dbartsSpec()` (spec.R:251, `resolvePriorScale` in
R/model.R), and converted against the response transform by a private
`Chain::resolvedNodeScale` helper (CH:3416) shared by the single-forest
constructor and every `setModel` reinstall, so a round trip through the model
SEXP no longer reverts a named calibration. S2 adds the reader,
`Chain::forestCalibration` (CH:980) - the AUTHORITATIVE report of what is in
force, independent of the model's recorded intent, so a `setResponse` /
`setOffset` at `updateScale = TRUE` or a `setData` shows up as a move rather
than staying silent - and the writer, `Chain::setForestPriorScale` (CH:1003),
sharing one `priorScaleFactor` conversion (CH:3426) with the reader so neither
direction can drift from the other; both are total over the four leaf models
and carry no family switch (facade FAC:284, 291; `Sampler` SAM:1182, 1191; R5
`dbartsSampler$getCalibration`/`$setCalibration` dbarts.R:1359, 1364; bridge
`bartcore_getCalibration`/`bartcore_setCalibration` RIB:3820, 3867). Refused
under a `k` hyperprior (the `sd` spelling only, since a sampled `k` has no
single value to divide by, or once the chains' `k` have diverged) and for
BCF/multinomial forests at creation and again mid-chain (see [f23]);
`prior.mean` is refused as not writable, naming the `setOffset` recipe. `NaN`
is refused as a malformed value rather than read as the unnamed spelling, both
at creation (R/model.R:1031) and mid-chain (`validateLiveScale`,
R/model.R:1040). Shipped tests: inst/tinytest/test-calibration-creation.R (two
composed probit arms at construction ranges 16x apart agree to 1e-12 under a
shared name, against 8.6 and 2.5 unnamed), inst/tinytest/test-calibration-prior-draws.R
(what the named quantity means per leaf model - exact for the constant leaf,
bounded inequalities for linear/gp/monotone - across ten family and decoration
paths, [f29] adding the tenth), inst/tinytest/test-calibration-midchain.R (a
read-then-write is BITWISE inert on the internal scale, not merely on the
report, since the setter skips a write reproducing either spelling of what is
in force; a write-then-read is pinned at ULP level rather than bitwise, since
`(P / f) * f != P` for about a tenth of positive pairs - 7da36dc3 corrected
which of the fixture's four `(m, P)` pairs actually rounds and added the
rounding cell, `m = 25` at `P = 0.25`, so the tolerance is exercised rather
than merely permitted; plus the refusal matrix, every mutation channel the
reported value must not surprise on, and the save/load round trip), and a
component test at the engine boundary (`testForestCalibration`,
tests/cpp/test_sampler.cpp:4300). What remains is **S3**: the flat C entries
(`inst/include/dbarts/dbarts.h`) - the R-user-facing capability above is
already `S`, since the flat-C gap does not gate a `dbartsSampler` caller -
riding the dbarts.h reshape's S1, with signatures frozen since S0.

[f17] Zero weights are accepted, not refused (A_class.R:572-576 errors only
below zero and warns that zeros are ignored; bridge RIB:4557). The conditionals
are exact - leaf suffstats multiply by `w` (MOD:314, 1122), and the sigma
posterior counts only positive-weight rows (`numPositiveWeights_` MOD:2729,
recounted on every install at MOD:2919, consumed MOD:2743-2747). The one named
inexactness against a true subset fit is CLOSED (`empty-leaf-veto-fix`,
2026-08-12): the empty-leaf veto counts
POSITIVE-WEIGHT members, so a leaf held alive only by zeroed rows is empty and
its branch is vetoed, on the conjugate path (MOV) and the constrained-leaf path
(MOD) alike. Occupancy elsewhere - the birth scan's `count`,
`collapseEmptyNodes`' trigger, `stateIsValid` - still counts members
deliberately, so this does NOT make zero-weight occupancy match a compacted fit;
see docs/design/empty-leaf-veto.md, "What counts as empty". The same fix covers
BCF and the Student-t row.

[f18] No latent vector exists: gaussian, BCF and heteroscedastic all leave
`ResponseModel::latents()` at its nullptr default (MOD:2636), and the bridge
returns `R_NilValue` (RIB:5599).

[f19] A Student-t fit records `family = "gaussian"`, so
`extract(type = "loglik")` takes the gaussian branch and evaluates `dnorm`
against the scalar `object$sigma` (generics.R:48-53). The t marginal is not
distinguished and the per-row `lambda_i` is not stored on the fit. No guard, no
doc, no TODO covers this - flagged, unadjudicated.

[f20] Weights on logistic are PG copy counts and a zero count is refused by
name at creation ("drop zero-count rows", RIB:1585-1590; R mirror
spec.R:131-140), so zero-weight subsetting is foreclosed for this family by the
weight semantics themselves - it is exactly the hole `latent-subset-mask` S2
filled (87d370ea), by the mid-chain `setActiveRows` channel rather than by any
change to the zero-count creation refusal, which stands.

[f21] The GLOBAL channel shipped at 8b047f8b, landing on the softmax coupling
rather than the response, which holds no precisions of its own to compose a
mask into: `MultinomialResponse::setActiveRows` (MOD:3671) is a pass-through
that only advertises the capability (`supportsActiveRows` MOD:3670), and
`Chain::setActiveRows` forwards the mask to
`MultinomialForestCombiner::setActiveRows` (COM:936) after the response's own
install. An inactive row's K interleaved Polya-Gamma draws are SKIPPED rather
than drawn and discarded, in `drawForestGlue` (COM:1007), and its composed
precision is zeroed in every category in `formForestResponse`
(COM:1064-1065); the row keeps its leaf occupancy and its reported softmax
probabilities, and omega is never zeroed since the working response divides
by it. PER-FOREST masking stays REFUSED, permanently and on model grounds:
the softmax margin is a log-sum-exp over the other K-1 forests, so a row
absent from category k's forest is still in every other category's
likelihood, and "row i is out of category k only" restricts no likelihood at
all. The refusal lands at the only reachable per-forest, per-observation
channel, `bartcore_setForestWeights` (RIB:3700-3704), naming the model reason
rather than "unbuilt". BCF's per-forest weight acceptance at that same
channel stands unaffected - a different (additive) coupling where the
per-forest mask is redundant with, not incoherent under, the combined
likelihood (see [f26]).

[f22] Multinomial's omegas live in the combiner, not the response model, and
`MultinomialResponse` does not override `latents()` (MOD:3635-3712), so
`getLatents` returns NULL. No accessor exposes them.

[f23] A named `prior.scale` is refused for BCF and multinomial forests both AT
CREATION and MID-CHAIN, by design - their per-forest leaf scales come from a
calibration map that owns them (map at COM:285-287), so a named value has
nowhere to land. Three creation-time refusal sites shipped at c2a7e89b: R-side
`dbartsSpec()`'s BCF composition (spec.R:408), the engine's own
BCF-composition gate (`refuseUnsupportedBCFComposition`, RIB:2238), and the
multinomial forest builder (`buildMultinomialSampler`, RIB:3041). S2
(d809b944) adds the mid-chain refusals, at TWO independent sites rather than
one shared gate: `$setCalibration`'s R5 method refuses BCF through
`refuseBCFMutation` (dbarts.R:1373, MEASURED "two-forest calibration map",
test-calibration-midchain.R:341-344) before ever reaching the bridge, and
refuses a multinomial fit's host shell through `refuseHostMutation`
(dbarts.R:1372, MEASURED "host sampler of a bart2", line 375-378); underneath
both, the engine-level gate any DIRECT low-level call still hits -
`Chain::setForestPriorScale` returning false whenever `combiner_ != nullptr`
(CH:1003), surfaced as `Rf_error(...calibrationMapName...)` at the bridge
(RIB:3880) - is what the unexported `dbarts:::bartcoreSetForestPriorScale`
hits on a multinomial forest's low-level handle (MEASURED "softmax calibration
map", line 356-359); the R5 layer never routes a BCF sampler there since
`refuseBCFMutation` refuses first, so only the multinomial arm exercises the
bridge gate directly. These cells stay `R`.

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and errors at
generics.R:86-91, but each component fit (`$occupancy` probit, `$positive`
gaussian) supports `extract(type = "loglik")` on its own.

[f26] SHIPPED, pinned bitwise at mask S4 (93afd635). A BCF sampler's response
IS a `GaussianResponse`, so nothing on the path gates on the coupling: the
shape probe (FAC:86) reports `supportsActiveRows`, the mask composes into the
case weights at MOD:2814 (so the sigma df is inherited) and then into the
per-forest weights at CH:1160. MEASURED at 6db22aee on a 200-row two-forest
sampler: `$setActiveRows(a)` and the bridge `bartcoreSetActiveRows` are both
accepted; on a sampler carrying `w` the mask is BITWISE `setWeights(w * a)` in
`train` and in `sigma`; an all-zeros mask runs finite; a fractional element is
refused. PINNED at mask S4: `inst/tinytest/test-active-rows-pins.R:84-112`
(masked-bcf, bitwise vs `setWeights(w * a)` on train and sigma) and the
`bcf-equivalence.R` `masked` scenario, recorded in
`bcf-equivalence-8b047f8b.rds`. A per-forest mask is refused as REDUNDANT rather
than unbuilt: `setForestWeights` (RIB:3688) already expresses it - though note
that channel is deliberately NOT row removal (CH:946-965: it does not remove
the row from occupancy, the combination or the sigma df; it DOES reach that
forest's empty-leaf veto, which counts positive composed weights), and it is
MISSING from the R5 object, reachable only through the unexported
`bartcoreSetForestWeights` (bartcore.R:999).

[f27] Delegating / decorating, and that is what 6db22aee landed for the
active-rows column: neither row took an edit of its own. `GroupedResponse`
forwards `setActiveRows` to its base (MOD:4743) exactly as it forwards
`setWeights` (MOD:4732), advertising the base's capability (MOD:4740), and
`drawGroupEffects` already weights its per-group sums by `workingWeights()`
(MOD:4668), so an inactive row leaves its group's mean and precision and an
all-inactive group falls back to its prior through the same formula. The
heteroscedastic `formMeanWeights` (CH:3611-3618) reads
`response_->workingWeights()` at CH:3614 - the COMPOSED `w * a` while a mask is
installed - and divides by `s^2(x_i)`, so a zero stays a zero. Both are
PINNED: grouped at tests/cpp/test_sampler.cpp:1591 (an entirely inactive group
draws its effect from the prior, finite), heteroscedastic at
inst/tinytest/test-active-rows-pins.R:212-237. Both pins were STRENGTHENED from
FINITENESS-only to bitwise masked-vs-`setWeights(w * a)` at mask S2
(87d370ea): grouped's effects, training fits and sigma all agree bitwise
against a composed-weight sampler (test_sampler.cpp:1638-1642); heteroscedastic
likewise for train and varcount (test-active-rows-pins.R:235-236). The same
delegation carries BOTH halves of the
nameable-calibration column for grouped: `GroupedResponse::fitScale()`/`fitShift()`
forward to their base (MOD:4774), so `Chain::resolvedNodeScale` (CH:3416) at
creation and `Chain::forestCalibration`/`setForestPriorScale` (CH:980, 1003)
mid-chain all convert a named `prior.scale` exactly as they do for the
undecorated family, with no edit of grouped's own. Creation-time, grouped is
one of the nine family/decoration paths c2a7e89b's own test measures
(inst/tinytest/test-calibration-prior-draws.R:221, MEASURED 0.74210 vs a 0.75
target); the mid-chain half rides the same generic mechanism every
non-coupled family does, with no dedicated grouped arm of its own in
test-calibration-midchain.R.

[f28] A heteroscedastic fit also records `family = "gaussian"` and takes the
same gaussian loglik branch with the SCALAR `object$sigma`, ignoring the
per-observation `s.train` surface the fit stores separately (bart.R:210-226).
Same standing as [f19]: no guard, no doc, unadjudicated.

[f29] The variance forest is a separate leaf model entirely, outside
`forests_`, and is not addressable by the mid-chain `setCalibration`: a
SHIPPED door (nameable-calibration synthesis 2.6 item 7), not an open one -
`Chain::setForestPriorScale`'s bounds check `f >= forests_.size()` (CH:1003)
never sees the variance forest, so a heteroscedastic sampler's
`shape.numForests` is 1 and `forest = 2` is refused by the ORDINARY
out-of-range check every single-forest sampler hits (`bartcore_getCalibration`
RIB:3825, `bartcore_setCalibration` RIB:3873), not by a hetero-specific gate.
The MEAN forest's own calibration - both halves - is not gated by that door:
`Chain::resolvedNodeScale` (CH:3416) runs at forest.leaf.scale assignment
(CH:588-590) before the variance-forest branch (CH:679), and the mid-chain
reader/writer (`forestCalibration`/`setForestPriorScale`, CH:980, 1003) read
the same `response_->fitScale()`/`fitShift()` - none of the three reads any
family flag. Creation-time it is now PINNED rather than merely constructing:
7da36dc3 added a heteroscedastic arm, the tenth of
`test-calibration-prior-draws.R`'s `anchorSamplers` sweep
(inst/tinytest/test-calibration-prior-draws.R:251-256), measuring a named
`scale = 1.5` against the shared 0.75 target inside the same 9% band
(`familyBand`, line 178) every other family in that loop is held to (line
258-264) - closing the earlier gap where the conversion ran unguarded but no
test exercised it. Mid-chain, the mean forest rides the same generic S2
mechanism as every other single-forest family, with no dedicated
heteroscedastic arm of its own beyond that shared mechanism and the
`forest = 2` refusal above.

[f30] Constructs today with no refusal site and no test. `spec.R:337-345`
refuses a variance forest only for a non-gaussian family or monotone
constraints, and a Student-t or grouped sampler still reports family
`gaussian`, so `resid.dist = student()` + `variance =` and grouped + `variance =`
both build (CH:581 decorates before CH:679 builds the variance forest). Whether
two scale mixtures on the same precision channel is a model anyone wants is not
adjudicated anywhere.

[f31] Recorded but UNBUILT doors, refused with that reason in the comment:
grouped ordinal because the cutpoint block and the group block are not yet shown
to interleave (RIB:2786-2790, ordinal.md section 8), grouped nbinom the same for
the dispersion block (RIB:2791-2795, negative-binomial.md section 7).

[f32] No surface at all: `applyGroupAttribute` (RIB:1908) is called from exactly
one site, RIB:2784, on the single-forest holder path, so `bartcore.groups` is
never read for a multinomial sampler.

[f33] DEFECT, not a refusal. `buildMultinomialForest` hard-sets
`forest.useDart = false` (CH:4524) while `bart2Multinomial` accepts and threads
`dart` through (bart.R:1219, 1235) and `buildMultinomialSampler` (RIB:3026-3068)
refuses nothing - unlike BCF, which names the option it drops (spec.R:393,
CH:4459). A user passing `dart = TRUE` gets a non-DART model silently.
multinomial.md does not mention DART.

[f34] `bart2Hurdle` builds both component calls with `redirectCall`
(bart.R:1987, 1995), so a user's `variance =` is forwarded to BOTH - including
the occupancy component, which then sets `family = "probit"` and would hit the
non-gaussian variance refusal at spec.R:338. bart.Rd:281 nonetheless describes
consuming the positive part's per-observation `sigma(x)`. The two readings are
in tension; unresolved.

[f35] `dart` is forwarded to both components (bart.R:922), each of which is an
ordinary single-forest chain that takes it.

[f36] No family gate: `installForests` checks shape, grid, DART and
variance-forest presence only (SAM:685-735), matching donor forest counts at
SAM:717, and `growForestFromRoot` loops every forest (CH:1589). Neither is
exercised by a BCF test, and BCF has no `bart2()` surface, so both are reached
only through the R5 `$installTrees` (dbarts.R:1489) / `$growFromRoot`
(dbarts.R:841).

[f37] `rbart_vi()` carries no `warm.start` or `n.grow.sweeps` formal
(rbart.R:45-51) and its unknown-argument check (rbart.R:60-69) rejects them. The
underlying R5 sampler carries no group gate on either path, so this is a surface
gap, not an engine one.

[f38] The MEAN forest keeps DART; the variance forest never takes it
(`buildVarianceForest` CH:3575 never sets `useDart`, default false at CH:125).

[f39] Current baselines: `equivalence-8b047f8b.rds` (37 scenarios),
`bcf-equivalence-8b047f8b.rds` (12), `multinomial-equivalence-1027be5.rds` (10)
- benchmarks/baselines/MANIFEST:16, 42, 48. Scenario names are the keys in
`makeScenarios()`, benchmarks/R/equivalence.R:60.

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
thin = 30 and cross into the band at 5x the spacing. Read as the r-vs-psi ridge
mixing slowly (H-MIX), on two ladder points rather than three; the recorded
full-R third point is still owed.

[f43] Aggregate `p_k(x*)` passes at both chain lengths; the three raw per-forest
`f_ik` cells carry a persistent U the ladder does not shrink. Pre-registered
suspect: `MultinomialForestCombiner::afterCombine`'s level-centering draw. OPEN.

[f44] AFT is exercised only in combination, through the `grouped_aft` scenario
(equivalence.R:436). There is no standalone AFT equivalence scenario; the
separate exact oracle benchmarks/R/aft-exact.R is not a MANIFEST entry.

[f45] Out of the SBC matrix by scope, each for its own recorded reason
(sbc-family-tiers.md:43-51): aft because its censoring status is fixed at
creation, so a prior-draw replication cannot vary it (the enabler is a status
setter); hazard and hurdle because their person-period / two-part designs depend
on `y0`, which breaks exchangeability, and because neither owns any sampling
code.

[f46] Tier A PASS with the sigma channel resolved as H-MIX on the (a, mu) ridge
(sbc-calibration.md:650-660). Explicitly out of the family-tiers matrix, and
`runSbcBCF` errors at that plan's HEAD (repair tracked separately in
docs/plans/runsbcbcf-repair.md).

[f47] OUT but DEFERRED rather than blocked: prior draws never reach
`varianceForest_` today, and the arm is liftable R-side through `setState`
(sbc-family-tiers.md:50-51).

## Gaps

Every MISSING (`M`) and UNVERIFIED (`?`) cell above, as candidate work items,
grouped by family. Scheduling is VD's and the orchestrator's; nothing here
carries a schedule, and REFUSED cells are deliberately absent - they are part of
the models.

**gaussian.** None; every column is S or an intentional `-`. The one standing
inexactness (zero-weight rows survived the empty-leaf veto) closed at
`empty-leaf-veto-fix` ([f17]).

**student (Gaussian + Student-t residuals).** No `rbart_vi()` surface
(rbart.R:60); no `xbart()` surface (xbart.R:9-38). Pointwise loglik evaluates
the Gaussian density, not the t marginal, with no guard ([f19]). Composition
with a variance forest constructs unrefused and untested ([f30]). Only one
dedicated tinytest file.

**probit.** None.

**logistic.** No `bart()` reach ([f1] - by signature). No `rbart_vi()` token
(rbart.R:49), so grouped logistic is engine-reachable but not R-reachable. No
flat-C test coverage.

**ordinal.** No `xbart()` or `rbart_vi()` reach (both refused at data.R:468 as
unsupported response shapes). Pointwise log-likelihood unbuilt
(generics.R:86-91). Grouped ordinal is a recorded unbuilt door (RIB:2789).
`warm.start` and `n.grow.sweeps` unbuilt for the arc (bart.R:493). Its selecting
control attribute is undocumented in the shipped header ([f3]). One dedicated
tinytest file. SBC gamma3 resolved but not re-run at full R.

**nbinom.** As ordinal, plus: no `bart()` token, no `xbart()` token, no
`rbart_vi()` token. Grouped nbinom a recorded unbuilt door (RIB:2794). Pointwise
loglik unbuilt. Header attribute undocumented. One dedicated tinytest file. SBC
`r`/`agg.psi` flag standing (H-MIX read, third ladder point owed). Real-valued
dispersion remains a recorded door (TODO `negbin-real-dispersion`).

**multinomial.** No flat-C creation path at all ([f4]). No mutable
`dbartsSampler` - only a host shell ([f11]) - and its three real channels
(`setCounts`, `setCategoryOffset`, `setCategoryTestOffset`) are unexported with
no R5 methods. No `getLatents` ([f22]). No pointwise loglik. No grouped surface
([f32]). No warm start / grow-from-root. **DART is accepted and silently
dropped** ([f33]) - the only cell in this matrix that is a defect rather than an
absence. SBC raw `f_ik` OPEN.

**aft.** No `xbart()` token. Pointwise loglik ships, but `setWeights` is refused
and the censoring status is fixed at creation, which is also what keeps AFT out
of the SBC matrix - a status setter is the named enabler ([f45]). No standalone
equivalence scenario ([f44]). `bart()` reach is undocumented ([f5]).

**hazard.** No flat-C token and no engine code of its own ([f6]); no
`rbart_vi()`, `xbart()` or `dbartsSpec()` reach. Out of SBC by design. One
dedicated tinytest file.

**hurdle.** No `dbarts()` sampler by construction ([f12]); no `rbart_vi()`,
`xbart()` or flat-C reach. No top-level pointwise loglik ([f25]). No warm start
/ grow-from-root. Composition with a variance forest is in tension between the
code path and the documentation ([f34]) - resolve before scheduling anything
that depends on it.

**bcf.** No `bart2()` surface, by the resolved fork that puts `bcf()` in
bartCause ([f7]). `setForestWeights` reaches the bridge but has no R5 method
([f26]). No pointwise loglik. Warm start and
grow-from-root are unrefused and untested for two forests ([f36]). Whole-data
`setData` stays undesigned (door 1 of the model-space survey).

**grouped.** `rbart_vi()`-only surface ([f8]): no `bart2()`, `dbarts()`,
`xbart()` or `dbartsSpec()` reach, and no `warm.start` / `n.grow.sweeps` formals
([f37]) though the engine paths carry no group gate. `setResponse` is an unbuilt
door (RIB:4309). Composition with a variance forest constructs unrefused and
untested ([f30]).

**heteroscedastic.** No `xbart()` reach. Pointwise loglik ignores the
per-observation `s(x)` surface it stores ([f28]). Selecting attribute
undocumented in the header ([f3]). Out of the SBC matrix, deferred not blocked
([f47]). One recorded unbuilt door from its own arc: the `setState` variance
column-mask gap (variance-forest-mutation-routing.md:499-500).

**Cross-cutting.** `nameable-calibration` is now R-SIDE COMPLETE, both halves
LANDED (S0 4c866286, S1 c2a7e89b, S2 d809b944 + 7da36dc3) - what remains is
S3, the flat-C entries, which rides the dbarts.h reshape's S1 and does not
gate any `dbartsSampler` caller. `latent-subset-mask` is ARC COMPLETE (S0
dc11a805, S1 6db22aee, S2 87d370ea, S3 8b047f8b, S4 93afd635) - the mask
covers every response family and every R-facing surface; what remains is the
flat-C entry, which rides the dbarts.h reshape's S1, same as
`nameable-calibration` S3; table 3 now carries no `P` cell.
