# Response-model feature matrix

Status: living reference, current at e218393d (2026-08-12). Carries no landing
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

Path aliases used in anchors (line numbers are at e218393d):

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
| gaussian | Gaussian (`ResponseFamily::gaussian`, `GaussianResponse` MOD:2683) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, `TResponse` MOD:3695) |
| probit | Binary probit (`ProbitResponse` MOD:2903) |
| logistic | Binary logistic, weights = observation counts (`LogisticResponse` MOD:3289) |
| ordinal | Ordered categorical, cumulative probit (`OrdinalResponse` MOD:3016) |
| nbinom | Negative binomial, positive-integer dispersion (`NBResponse` MOD:3956) |
| multinom | Multinomial softmax, K forests (`MultinomialResponse` MOD:3416 + combiner) |
| aft | AFT survival, log-normal (`AFTResponse` MOD:3481) |
| hazard | Discrete-time hazard (person-period sugar, dbarts.R:437-486) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, bart.R:1936) |
| bcf | Bayesian causal forest, two forests (`BCFForestCombiner` COM:698) |
| grouped | Grouped random intercepts (`GroupedResponse` MOD:4258) |
| hetero | Heteroscedastic variance forest (CH:650) |

The engine's `ResponseFamily` enum has only six tokens (MOD:2501: gaussian,
probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf, grouped
and hetero are all reached some other way, which is exactly why they need rows
here rather than an enum read. Leaf models (constant, monotone, linear, GP) are
an orthogonal axis and are not rows; where a leaf model gates a capability the
cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `rbart_vi()` | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|---|
| gaussian | S bart.R:2374 | S bart.R:550 | S dbarts.R:359 | S rbart.R:48 | S xbart.R:27 | S CAPI:341 |
| student | S bart.R:2278 | S bart.R:564 | S dbarts.R:340 | M rbart.R:59 | M xbart.R:9-38 | S RIB:2503 [f2] |
| probit | S bart.Rd:151 | S bart.R:551 | S dbarts.R:360 | S data.R:523 | S xbart.R:27 | S CAPI:340 |
| logistic | - [f1] | S bart.R:552 | S dbarts.R:361 | R rbart.R:48 | S xbart.R:27 | S RIB:1550 |
| ordinal | R bart.R:2378 | S bart.R:555 | S dbarts.R:363 | R data.R:468 | R data.R:468 | S RIB:1558 [f3] |
| nbinom | - [f1] | S bart.R:556 | S dbarts.R:364 | M rbart.R:48 | M xbart.R:27 | S RIB:1565 [f3] |
| multinom | - [f1] | S bart.R:554 | R dbarts.R:357 | R data.R:468 | R data.R:468 | M [f4] |
| aft | ? bart.R:2374 [f5] | S bart.R:553 | S dbarts.R:362 | S rbart.R:48 | M xbart.R:27 | S CAPI:343 |
| hazard | - [f1] | S bart.R:557 | S dbarts.R:365 | M rbart.R:48 | M xbart.R:27 | M [f6] |
| hurdle | - [f1] | S bart.R:560 | R dbarts.R:395 | M | M | M [f6] |
| bcf | - [f1] | R bart.R:616 [f7] | S dbarts.R:350 | R rbart.R:59 | M | S CAPI:348 |
| grouped | - [f1] | M [f8] | M [f8] | S rbart.R:365 | M | S RIB:1954 [f3] |
| hetero | - [f1] | S bart.R:540 | S dbarts.R:346 | R rbart.R:59 | M | S RIB:2044 [f3] |

`dbartsSpec()` (spec.R:484-492) resolves the seven single-forest tokens - auto,
gaussian, probit, logistic, aft, ordinal, nbinom - plus BCF through its
`treatment =` argument (spec.R:479) and a variance forest through `variance =`
(spec.R:475); it does not reach multinomial, hazard, hurdle or grouped.

## 2. Mutation channels on the R5 `dbartsSampler`

The channels that make dbarts a conditional model inside an outer sampler.
`updateScale` is broken out because it is refused independently of the setter
it rides on.

| model | `setResponse` | `setOffset` | `updateScale = TRUE` | `setPredictor` (+ per-obs) | `setWeights` | `setSigma` | test surface |
|---|---|---|---|---|---|---|---|
| gaussian | S MOD:2716 | S MOD:2771 | S MOD:2771 | S RIB:4533, 4689 | S MOD:2759 | S RIB:2713 | S RIB:4232, 4285 |
| student | S MOD:3758 | S MOD:3765 | S MOD:3765 | S RIB:4533 | S MOD:3781 | S RIB:2713 | S RIB:4232 |
| probit | S MOD:2937 | S MOD:2943 | - [f9] | S RIB:4533 | R RIB:1626 | R RIB:2715 | S RIB:4232 |
| logistic | S MOD:3329 | S MOD:3335 | - [f9] | S RIB:4533 | R RIB:1619 [f10] | R RIB:2715 | S RIB:4232 |
| ordinal | S MOD:3051 | S MOD:3059 | - [f9] | S RIB:4533 | R RIB:1626 | R RIB:2715 | S RIB:4232 |
| nbinom | S MOD:4005 | S MOD:4012 | - [f9] | S RIB:4533 | R RIB:1626 | R RIB:2715 | S RIB:4232 |
| multinom | - dbarts.R:776 [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] |
| aft | S MOD:3537 | S MOD:3550 | S MOD:3550 | S RIB:4533 | R RIB:1626 | S RIB:2713 | S RIB:4232 |
| hazard | S MOD:2937 [f6] | S MOD:2943 | - [f9] | S RIB:4533 | R RIB:1626 | R RIB:2715 | S RIB:4232 |
| hurdle | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] |
| bcf | S COM:716 | S COM:716 | R bartcore.R:285 | S RIB:4533, 4689 | S RIB:4379 | S RIB:4155 | R RIB:2689 |
| grouped | R RIB:4139 [f13] | S MOD:4352 | S MOD:4352 | S RIB:4533 | S MOD:4366 [f14] | S RIB:2713 [f14] | S RIB:4232 |
| hetero | S RIB:2620 | S RIB:2620 | R RIB:2620 | S RIB:4533, 4689 | S RIB:2616 | R RIB:2711 | S RIB:4232 |

`setData` (whole-data replacement, n free) is single-forest and dense-store
only (RIB:3326) and is refused for grouped (RIB:4169) and aft (RIB:4171);
BCF/multinomial whole-data `setData` stays undesigned by the model-space
survey's verdict (model-space-survey.md doors 1 and 3).

## 3. Row subsetting, latents, calibration

| model | zero-weight row subset | active-rows mask [f15] | `getLatents` | pointwise loglik | nameable calibration [f16] |
|---|---|---|---|---|---|
| gaussian | S sampler.Rd:137, MOD:2691 [f17] | P S1 | - RIB:5423 [f18] | S generics.R:48 | P S1 + S2 |
| student | S MOD:3740-3745 [f17] | P S1 | S MOD:3791 | ? generics.R:48 [f19] | P S1 + S2 |
| probit | R RIB:1597 | P S1 | S MOD:2962 | S generics.R:54 | P S1 + S2 |
| logistic | R RIB:1601 [f20] | P S2 | S MOD:3353 | S generics.R:54 | P S1 + S2 |
| ordinal | R RIB:2494 | P S1 | S MOD:3080 | M generics.R:81 | P S1 + S2 |
| nbinom | R RIB:2500 | P S2 | S MOD:4036 | M generics.R:81 | P S1 + S2 |
| multinom | R RIB:2959 | P S3 [f21] | M MOD:3416 [f22] | M generics.R:81 | R [f23] |
| aft | R RIB:2490 | P S2 | S MOD:3579 | S generics.R:59 | P S1 + S2 |
| hazard | R RIB:1597 [f6] | P S1 [f6] | S MOD:2962 | S generics.R:54 [f24] | P S1 + S2 |
| hurdle | R bart.R:931 | - [f12] | - [f12] | M generics.R:81 [f25] | - [f12] |
| bcf | S COM:705-716 [f17] | P S1 [f26] | - [f18] | M generics.R:81 | R [f23] |
| grouped | S MOD:4221-4243 | P [f27] | S MOD:4375 | S generics.R:1391 | P [f27] |
| hetero | S CH:3470, MOD:305 | P [f27] | - [f18] | ? generics.R:48 [f28] | P [f29] |

## 4. Model composition

| model | variance forest | grouped ranef | DART | warm start | grow-from-root |
|---|---|---|---|---|---|
| gaussian | S FAC:633 | S CH:551 | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| student | ? [f30] | S CH:551 | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| probit | R spec.R:334 | S CH:551 | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| logistic | R spec.R:334 | S CH:551 | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| ordinal | R spec.R:334 | M RIB:2758 [f31] | S CH:495 | R bart.R:485 | R bart.R:485 |
| nbinom | R spec.R:334 | M RIB:2763 [f31] | S CH:495 | R bart.R:485 | R bart.R:485 |
| multinom | R FAC:724 | M RIB:1951 [f32] | M CH:4379 [f33] | R bart.R:485 | R bart.R:485 |
| aft | R spec.R:334 | S CH:551 | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| hazard | R spec.R:334 | M rbart.R:48 [f6] | S CH:495 | S bart.R:1000 | S dbarts.R:841 |
| hurdle | ? [f34] | M | S bart.R:910 [f35] | R bart.R:485 | R bart.R:485 |
| bcf | R FAC:711 | R RIB:2283 | R spec.R:389 | S SAM:717 [f36] | S CH:1467 [f36] |
| grouped | ? [f30] | - | S rbart.R:576 | M rbart.R:44-50 [f37] | M rbart.R:44-50 [f37] |
| hetero | - | ? [f30] | S CH:495 [f38] | S SAM:722 | S CH:1455 |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused at dbarts.R:844-853 and no-op at CH:1445, so every family above
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

[f1] `bart()` carries no `family` formal (bart.R:2246-2280); the family is
inferred from the response (numeric -> gaussian, 2-level -> probit) and
`resid.dist` is its only response-law lever. Families needing an explicit token
are out of reach by signature, not by refusal. Ordinal is the exception: it is
explicitly backstopped at bart.R:2378.

[f2] Student-t is not a token anywhere. It is selected by a finite `resid.df`
attribute on the model SEXP (RIB:2503-2514, gaussian-only gate) and refused for
every non-gaussian family R-side at spec.R:270-276. The engine family stays
`gaussian`, which is why the whole gaussian row applies to it in tables 2-4.

[f3] Reachable through `dbarts_sampler_create` but not discoverable from the
shipped header: ordinal needs the `bartcore.n.categories` control attribute
(RIB:433), nbinom `bartcore.dispersion` (RIB:446), grouped `bartcore.groups`
(RIB:1954), heteroscedastic `bartcore.variance` (RIB:2044). The header's
`family` documentation (CAPI:338-357) names only probit, logistic, gaussian,
aft and BCF.

[f4] Multinomial has no dbarts.h creation path at all; it is reached only by
the R-internal `.Call` entries `C_dbarts_bartcore_createMultinomial` /
`...Counts` (bartcore.R:798, 861).

[f5] A `Surv()` response passed to `bart()` auto-dispatches through `dbarts()`
and nothing refuses it (bart.R:2374), but man/bart.Rd:151 does not document it.
Reachable, undocumented, intent unadjudicated.

[f6] `family = "hazard"` / `"hazard.probit"` / `"hazard.logistic"` is
person-period ingestion sugar: dbarts.R:437-486 expands the design and remaps
the token to `"probit"` or `"logistic"` at dbarts.R:483 before any model is
built. The resulting sampler *is* an ordinary binary one, so its whole row
equals the probit (or logistic) row, and the fit records `family = "probit"`.
No engine code, hence no C-API token and no SBC arm. Refusals inside the
expander: no formula interface (:438), no `subset` (:453), no `test` (:456).

[f7] `treatment` is not a `bart2()` formal, so the unknown-argument check at
bart.R:616-626 kills it. BCF's public creation surface is `dbarts()` /
`dbartsSpec()`; `bcf()` and `bartBCF` ship in **bartCause**, not dbarts
(bcf-public-surface.md:483-496, fork 4 RESOLVED VD 2026-08-11; echoed
bcf.md:344-352). This row therefore tracks the engine capability, not a dbarts
user-facing verb.

[f8] `bartcore.groups` is written at exactly one site, rbart.R:365, and no other
entry point carries a `group.by` formal, so grouped random effects are an
`rbart_vi()`-only surface.

[f9] `updateScale` re-derives the internal response transform. The latent
families have `fitScale() == 1` and `fitShift() == 0` by definition, so there
is no transform to re-anchor and the flag is ignored rather than refused.

[f10] Logistic weights are the observation counts its Polya-Gamma latents were
built from, so they cannot change after creation (RIB:1619). This is the one
weight refusal that is "unbuilt" rather than "incoherent": creation accepts
positive-integer weights at RIB:1601.

[f11] `bart2(family = "multinomial")` returns a `dbartsSampler` that is only a
HOST SHELL - it carries the design and priors but not the model, which lives on
the separate bartcore handle `$bc` (bart.R:1242, 1318). Every R5 mutator is
refused by `refuseHostMutation` (dbarts.R:770-783). At the handle level:
`setResponse` R RIB:2586 and `setOffset` R RIB:2583 (both redirect to the counts
channel), `setWeights` R RIB:2593, `setSigma` R RIB:2715, `setTestOffset` R
RIB:2358 (a flat offset leaves the simplex), `setPredictor` S RIB:4533,
`setTestPredictor` S RIB:4232. Its own channels - `setCounts` RIB:3479,
`setCategoryOffset` RIB:3557, `setCategoryTestOffset` RIB:3596 - are S at the
bridge but unexported (`dbarts:::`) and absent from the R5 object.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction at
dbarts.R:395 and `bart2Hurdle` (bart.R:1936) composes two ordinary `bart2()`
fits - an occupancy probit (bart.R:1968) and a lognormal positive part
(bart.R:1977) - glued at report time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f13] Grouped samplers fix the response at creation (RIB:4139). The engine
method exists and delegates faithfully (MOD:4341), so this is an unbuilt door,
not a model refusal - but it errors, so callers see a refusal.

[f14] Reads off the BASE family: grouped gaussian takes `setWeights` (MOD:4366)
and `setSigma`; grouped probit is refused on both (RIB:1626, RIB:2715); grouped
aft takes `setSigma` and refuses `setWeights`.

[f15] Arc `latent-subset-mask` (TODO:155), design FINAL and decision-gated on
its own scope; artifacts .claude/latent-subset-mask-design/. A first-class 0/1
`setActiveRows` channel each family composes into its own precision vector,
with the latent draw skipped for inactive rows. Slices: **S0** pins (no engine
change); **S1** the channel plus gaussian, Student-t, probit, ordinal; **S2**
logistic, nbinom, aft; **S3** multinomial (global only); **S4** surface,
records, baselines. Nothing is built today (no `setActiveRows` symbol exists in
src/, R/, inst/ or man/).

[f16] Arc `nameable-calibration` (TODO:279), design AMENDED FINAL, SCHEDULED
pre-release; artifacts .claude/nameable-calibration-design/. Names the
per-forest prior ANCHOR (`prior.scale`, the forest-total prior scale at k = 1,
in response units) rather than an sd, with a `$getCalibration` /
`$setCalibration` pair. Slices: **S0** signature freeze; **S1** creation half;
**S2** mid-chain get/set; **S3** the flat-C half, executed inside the dbarts.h
reshape's S1. Nothing is built today.

[f17] Zero weights are accepted, not refused (A_class.R:536-540 errors only
below zero and warns that zeros are ignored; bridge RIB:4387). The conditionals
are exact - leaf suffstats multiply by `w` (MOD:313, 1118), and the sigma
posterior counts only positive-weight rows (`numPositiveWeights_` MOD:2691,
consumed MOD:2705-2709). The one named inexactness against a true subset fit is
CLOSED (`empty-leaf-veto-fix`, 2026-08-12): the empty-leaf veto counts
POSITIVE-WEIGHT members, so a leaf held alive only by zeroed rows is empty and
its branch is vetoed, on the conjugate path (MOV) and the constrained-leaf path
(MOD) alike. Occupancy elsewhere - the birth scan's `count`,
`collapseEmptyNodes`' trigger, `stateIsValid` - still counts members
deliberately, so this does NOT make zero-weight occupancy match a compacted fit;
see docs/design/empty-leaf-veto.md, "What counts as empty". The same fix covers
BCF and the Student-t row.

[f18] No latent vector exists: gaussian, BCF and heteroscedastic all leave
`ResponseModel::latents()` at its nullptr default (MOD:2598), and the bridge
returns `R_NilValue` (RIB:5423).

[f19] A Student-t fit records `family = "gaussian"`, so
`extract(type = "loglik")` takes the gaussian branch and evaluates `dnorm`
against the scalar `object$sigma` (generics.R:48-53). The t marginal is not
distinguished and the per-row `lambda_i` is not stored on the fit. No guard, no
doc, no TODO covers this - flagged, unadjudicated.

[f20] Weights on logistic are PG copy counts and a zero count is refused by
name at creation ("drop zero-count rows", RIB:1601-1606; R mirror
spec.R:131-140), so zero-weight subsetting is foreclosed for this family by the
weight semantics themselves - it is exactly the hole `latent-subset-mask` S2
fills.

[f21] Per-forest masking is REFUSED on model grounds and lands with that reason
in S3: the softmax margin is a log-sum-exp over the other K-1 forests, so a row
absent from category k's forest is still in every other category's likelihood
and "row i is out of category k only" restricts no likelihood. Only a GLOBAL
mask is well-posed here.

[f22] Multinomial's omegas live in the combiner, not the response model, and
`MultinomialResponse` does not override `latents()` (MOD:3416-3480), so
`getLatents` returns NULL. No accessor exposes them.

[f23] `setCalibration` refuses BCF and multinomial forests by design - their
per-forest scales come from a calibration map that owns them
(nameable-calibration synthesis 2.6 item 3, map at COM:285-287).

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and errors at
generics.R:81-86, but each component fit (`$occupancy` probit, `$positive`
gaussian) supports `extract(type = "loglik")` on its own.

[f26] BCF composes the mask into its `GaussianResponse` (so the sigma df is
inherited) and then into the per-forest weights at CH:1086. A per-forest mask is
refused as REDUNDANT rather than unbuilt: `setForestWeights` (RIB:3649) already
expresses it - though note that channel is deliberately NOT row removal
(CH:915-931: it does not remove the row from occupancy, the combination or the
sigma df; it DOES reach that forest's empty-leaf veto, which counts positive
composed weights), and it is MISSING from the R5 object, reachable only through
the unexported `bartcoreSetForestWeights` (bartcore.R:999).

[f27] Delegating / decorating: `GroupedResponse` forwards `setActiveRows` to its
base exactly as it forwards `setWeights` (MOD:4366), and the heteroscedastic
`formMeanWeights` divides the composed weight by `s^2(x_i)` so a zero stays a
zero (CH:3466-3473). Both need a PIN, not an edit. The same delegation carries
nameable calibration for grouped.

[f28] A heteroscedastic fit also records `family = "gaussian"` and takes the
same gaussian loglik branch with the SCALAR `object$sigma`, ignoring the
per-observation `s.train` surface the fit stores separately (bart.R:210-226).
Same standing as [f19]: no guard, no doc, unadjudicated.

[f29] The variance forest is outside `forests_` and is not addressable by the
calibration setter; recorded as a door in the design (synthesis 2.6 item 7).
The MEAN forest of a heteroscedastic sampler is addressable as normal.

[f30] Constructs today with no refusal site and no test. `spec.R:333-341`
refuses a variance forest only for a non-gaussian family or monotone
constraints, and a Student-t or grouped sampler still reports family
`gaussian`, so `resid.dist = student()` + `variance =` and grouped + `variance =`
both build (CH:551 decorates before CH:650 builds the variance forest). Whether
two scale mixtures on the same precision channel is a model anyone wants is not
adjudicated anywhere.

[f31] Recorded but UNBUILT doors, refused with that reason in the comment:
grouped ordinal because the cutpoint block and the group block are not yet shown
to interleave (RIB:2755-2759, ordinal.md section 8), grouped nbinom the same for
the dispersion block (RIB:2760-2764, negative-binomial.md section 7).

[f32] No surface at all: `applyGroupAttribute` (RIB:1951) is called from exactly
one site, RIB:2753, on the single-forest holder path, so `bartcore.groups` is
never read for a multinomial sampler.

[f33] DEFECT, not a refusal. `buildMultinomialForest` hard-sets
`forest.useDart = false` (CH:4379) while `bart2Multinomial` accepts and threads
`dart` through (bart.R:1205, 1219) and `buildMultinomialSampler` (RIB:2995-3029)
refuses nothing - unlike BCF, which names the option it drops (spec.R:389,
CH:4314). A user passing `dart = TRUE` gets a non-DART model silently.
multinomial.md does not mention DART.

[f34] `bart2Hurdle` builds both component calls with `redirectCall`
(bart.R:1965, 1973), so a user's `variance =` is forwarded to BOTH - including
the occupancy component, which then sets `family = "probit"` and would hit the
non-gaussian variance refusal at spec.R:334. bart.Rd:278 nonetheless describes
consuming the positive part's per-observation `sigma(x)`. The two readings are
in tension; unresolved.

[f35] `dart` is forwarded to both components (bart.R:910), each of which is an
ordinary single-forest chain that takes it.

[f36] No family gate: `installForests` checks shape, grid, DART and
variance-forest presence only (SAM:685-735), matching donor forest counts at
SAM:717, and `growForestFromRoot` loops every forest (CH:1467). Neither is
exercised by a BCF test, and BCF has no `bart2()` surface, so both are reached
only through the R5 `$installTrees` (dbarts.R:1397) / `$growFromRoot`
(dbarts.R:841).

[f37] `rbart_vi()` carries no `warm.start` or `n.grow.sweeps` formal
(rbart.R:44-50) and its unknown-argument check (rbart.R:59-68) rejects them. The
underlying R5 sampler carries no group gate on either path, so this is a surface
gap, not an engine one.

[f38] The MEAN forest keeps DART; the variance forest never takes it
(`buildVarianceForest` CH:3430 never sets `useDart`, default false at CH:118).

[f39] Current baselines: `equivalence-21fc29c.rds` (35 scenarios),
`bcf-equivalence-a825263.rds` (11), `multinomial-equivalence-1027be5.rds` (10) -
benchmarks/baselines/MANIFEST:15, 40, 46. Scenario names are the keys in
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
(rbart.R:59); no `xbart()` surface (xbart.R:9-38). Pointwise loglik evaluates
the Gaussian density, not the t marginal, with no guard ([f19]). Composition
with a variance forest constructs unrefused and untested ([f30]). Only one
dedicated tinytest file.

**probit.** None.

**logistic.** No `bart()` reach ([f1] - by signature). No `rbart_vi()` token
(rbart.R:48), so grouped logistic is engine-reachable but not R-reachable. No
flat-C test coverage.

**ordinal.** No `xbart()` or `rbart_vi()` reach (both refused at data.R:468 as
unsupported response shapes). Pointwise log-likelihood unbuilt
(generics.R:81-86). Grouped ordinal is a recorded unbuilt door (RIB:2758).
`warm.start` and `n.grow.sweeps` unbuilt for the arc (bart.R:485). Its selecting
control attribute is undocumented in the shipped header ([f3]). One dedicated
tinytest file. SBC gamma3 resolved but not re-run at full R.

**nbinom.** As ordinal, plus: no `bart()` token, no `xbart()` token, no
`rbart_vi()` token. Grouped nbinom a recorded unbuilt door (RIB:2763). Pointwise
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
([f26]). No pointwise loglik. Warm start and grow-from-root are unrefused and
untested for two forests ([f36]). Whole-data `setData` stays undesigned (door 1
of the model-space survey).

**grouped.** `rbart_vi()`-only surface ([f8]): no `bart2()`, `dbarts()`,
`xbart()` or `dbartsSpec()` reach, and no `warm.start` / `n.grow.sweeps` formals
([f37]) though the engine paths carry no group gate. `setResponse` is an unbuilt
door (RIB:4139). Composition with a variance forest constructs unrefused and
untested ([f30]).

**heteroscedastic.** No `xbart()` reach. Pointwise loglik ignores the
per-observation `s(x)` surface it stores ([f28]). Selecting attribute
undocumented in the header ([f3]). Out of the SBC matrix, deferred not blocked
([f47]). One recorded unbuilt door from its own arc: the `setState` variance
column-mask gap (variance-forest-mutation-routing.md:499-500).

**Cross-cutting.** The `latent-subset-mask` and `nameable-calibration` arcs
account for every `P` cell in tables 3 and are the two largest scheduled
additions to this matrix; both are designed and unbuilt, and both are sequenced
before the dbarts.h reshape's S1.
