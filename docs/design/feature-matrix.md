# Response-model feature matrix

Status: living reference, current at c2a7e89b (2026-08-13). Carries no landing
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

Path aliases used in anchors (line numbers are at c2a7e89b):

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
| gaussian | Gaussian (`ResponseFamily::gaussian`, `GaussianResponse` MOD:2701) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, `TResponse` MOD:3842) |
| probit | Binary probit (`ProbitResponse` MOD:2977) |
| logistic | Binary logistic, weights = observation counts (`LogisticResponse` MOD:3436) |
| ordinal | Ordered categorical, cumulative probit (`OrdinalResponse` MOD:3116) |
| nbinom | Negative binomial, positive-integer dispersion (`NBResponse` MOD:4132) |
| multinom | Multinomial softmax, K forests (`MultinomialResponse` MOD:3563 + combiner) |
| aft | AFT survival, log-normal (`AFTResponse` MOD:3628) |
| hazard | Discrete-time hazard (person-period sugar, dbarts.R:437-486) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, bart.R:1958) |
| bcf | Bayesian causal forest, two forests (`BCFForestCombiner` COM:698) |
| grouped | Grouped random intercepts (`GroupedResponse` MOD:4434) |
| hetero | Heteroscedastic variance forest (CH:660) |

The engine's `ResponseFamily` enum has only six tokens (MOD:2504: gaussian,
probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf, grouped
and hetero are all reached some other way, which is exactly why they need rows
here rather than an enum read. Leaf models (constant, monotone, linear, GP) are
an orthogonal axis and are not rows; where a leaf model gates a capability the
cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `rbart_vi()` | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|---|
| gaussian | S bart.R:2398 | S bart.R:559 | S dbarts.R:359 | S rbart.R:49 | S xbart.R:27 | S CAPI:341 |
| student | S bart.R:2301 | S bart.R:573 | S dbarts.R:340 | M rbart.R:60 | M xbart.R:9-38 | S RIB:2525 [f2] |
| probit | S bart.Rd:151 | S bart.R:560 | S dbarts.R:360 | S data.R:523 | S xbart.R:27 | S CAPI:340 |
| logistic | - [f1] | S bart.R:561 | S dbarts.R:361 | R rbart.R:49 | S xbart.R:27 | S RIB:1566 |
| ordinal | R bart.R:2402 | S bart.R:564 | S dbarts.R:363 | R data.R:468 | R data.R:468 | S RIB:1574 [f3] |
| nbinom | - [f1] | S bart.R:565 | S dbarts.R:364 | M rbart.R:49 | M xbart.R:27 | S RIB:1581 [f3] |
| multinom | - [f1] | S bart.R:563 | R dbarts.R:357 | R data.R:468 | R data.R:468 | M [f4] |
| aft | ? bart.R:2398 [f5] | S bart.R:562 | S dbarts.R:362 | S rbart.R:49 | M xbart.R:27 | S CAPI:343 |
| hazard | - [f1] | S bart.R:566 | S dbarts.R:365 | M rbart.R:49 | M xbart.R:27 | M [f6] |
| hurdle | - [f1] | S bart.R:569 | R dbarts.R:395 | M | M | M [f6] |
| bcf | - [f1] | R bart.R:625 [f7] | S dbarts.R:350 | R rbart.R:60 | M | S CAPI:348 |
| grouped | - [f1] | M [f8] | M [f8] | S rbart.R:367 | M | S RIB:1972 [f3] |
| hetero | - [f1] | S bart.R:549 | S dbarts.R:346 | R rbart.R:60 | M | S RIB:2062 [f3] |

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
| gaussian | S MOD:2734 | S MOD:2814 | S MOD:2814 | S RIB:4612, 4768 | S MOD:2783 | S RIB:2735 | S RIB:4307, 4360 |
| student | S MOD:3911 | S MOD:3918 | S MOD:3918 | S RIB:4612 | S MOD:3935 | S RIB:2735 | S RIB:4307 |
| probit | S MOD:3030 | S MOD:3036 | - [f9] | S RIB:4612 | R RIB:1642 | R RIB:2737 | S RIB:4307 |
| logistic | S MOD:3476 | S MOD:3482 | - [f9] | S RIB:4612 | R RIB:1635 [f10] | R RIB:2737 | S RIB:4307 |
| ordinal | S MOD:3167 | S MOD:3175 | - [f9] | S RIB:4612 | R RIB:1642 | R RIB:2737 | S RIB:4307 |
| nbinom | S MOD:4181 | S MOD:4188 | - [f9] | S RIB:4612 | R RIB:1642 | R RIB:2737 | S RIB:4307 |
| multinom | - dbarts.R:776 [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] |
| aft | S MOD:3684 | S MOD:3697 | S MOD:3697 | S RIB:4612 | R RIB:1642 | S RIB:2735 | S RIB:4307 |
| hazard | S MOD:3030 [f6] | S MOD:3036 | - [f9] | S RIB:4612 | R RIB:1642 | R RIB:2737 | S RIB:4307 |
| hurdle | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] |
| bcf | S COM:716 | S COM:716 | R bartcore.R:285 | S RIB:4612, 4768 | S RIB:4454 | S RIB:4230 | R RIB:2711 |
| grouped | R RIB:4214 [f13] | S MOD:4528 | S MOD:4528 | S RIB:4612 | S MOD:4542 [f14] | S RIB:2735 [f14] | S RIB:4307 |
| hetero | S RIB:2642 | S RIB:2642 | R RIB:2642 | S RIB:4612, 4768 | S RIB:2638 | R RIB:2733 | S RIB:4307 |

`setData` (whole-data replacement, n free) is single-forest and dense-store
only (RIB:3356) and is refused for grouped (RIB:4244) and aft (RIB:4246);
BCF/multinomial whole-data `setData` stays undesigned by the model-space
survey's verdict (model-space-survey.md doors 1 and 3).

## 3. Row subsetting, latents, calibration

| model | zero-weight row subset | active-rows mask [f15] | `getLatents` | pointwise loglik | nameable calibration [f16] |
|---|---|---|---|---|---|
| gaussian | S sampler.Rd:137, MOD:2709 [f17] | S MOD:2794 | - RIB:5502 [f18] | S generics.R:48 | P S2 [f16] |
| student | S MOD:3891-3898 [f17] | S MOD:3946 | S MOD:3958 | ? generics.R:48 [f19] | P S2 [f16] |
| probit | R RIB:1613 | S MOD:3020 | S MOD:3056 | S generics.R:59 | P S2 [f16] |
| logistic | R RIB:1617 [f20] | P S2 | S MOD:3500 | S generics.R:59 | P S2 [f16] |
| ordinal | R RIB:2516 | S MOD:3152 | S MOD:3197 | M generics.R:86 | P S2 [f16] |
| nbinom | R RIB:2522 | P S2 | S MOD:4212 | M generics.R:86 | P S2 [f16] |
| multinom | R RIB:2981 | P S3 [f21] | M MOD:3563 [f22] | M generics.R:86 | R [f23] |
| aft | R RIB:2512 | P S2 | S MOD:3726 | S generics.R:64 | P S2 [f16] |
| hazard | R RIB:1613 [f6] | S MOD:3020 [f6] | S MOD:3056 | S generics.R:59 [f24] | P S2 [f16] |
| hurdle | R bart.R:944 | - [f12] | - [f12] | M generics.R:86 [f25] | - [f12] |
| bcf | S COM:705-716 [f17] | ? MOD:2794 [f26] | - [f18] | M generics.R:86 | R [f23] |
| grouped | S MOD:4397-4419 | S MOD:4553 [f27] | S MOD:4562 | S generics.R:1396 | P S2 [f27] |
| hetero | S CH:3539, MOD:305 | S CH:3535 [f27] | - [f18] | ? generics.R:48 [f28] | ? [f29] |

## 4. Model composition

| model | variance forest | grouped ranef | DART | warm start | grow-from-root |
|---|---|---|---|---|---|
| gaussian | S FAC:649 | S CH:562 | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| student | ? [f30] | S CH:562 | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| probit | R spec.R:338 | S CH:562 | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| logistic | R spec.R:338 | S CH:562 | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| ordinal | R spec.R:338 | M RIB:2780 [f31] | S CH:506 | R bart.R:493 | R bart.R:493 |
| nbinom | R spec.R:338 | M RIB:2785 [f31] | S CH:506 | R bart.R:493 | R bart.R:493 |
| multinom | R FAC:740 | M RIB:1969 [f32] | M CH:4448 [f33] | R bart.R:493 | R bart.R:493 |
| aft | R spec.R:338 | S CH:562 | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| hazard | R spec.R:338 | M rbart.R:49 [f6] | S CH:506 | S bart.R:1014 | S dbarts.R:841 |
| hurdle | ? [f34] | M | S bart.R:922 [f35] | R bart.R:493 | R bart.R:493 |
| bcf | R FAC:727 | R RIB:2305 | R spec.R:393 | S SAM:717 [f36] | S CH:1525 [f36] |
| grouped | ? [f30] | - | S rbart.R:578 | M rbart.R:45-51 [f37] | M rbart.R:45-51 [f37] |
| hetero | - | ? [f30] | S CH:506 [f38] | S SAM:722 | S CH:1513 |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused at dbarts.R:844-853 and no-op at CH:1503, so every family above
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
attribute on the model SEXP (RIB:2525-2536, gaussian-only gate) and refused for
every non-gaussian family R-side at spec.R:274-280. The engine family stays
`gaussian`, which is why the whole gaussian row applies to it in tables 2-4.

[f3] Reachable through `dbarts_sampler_create` but not discoverable from the
shipped header: ordinal needs the `bartcore.n.categories` control attribute
(RIB:438), nbinom `bartcore.dispersion` (RIB:451), grouped `bartcore.groups`
(RIB:1972), heteroscedastic `bartcore.variance` (RIB:2062). The header's
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
built from, so they cannot change after creation (RIB:1635). This is the one
weight refusal that is "unbuilt" rather than "incoherent": creation accepts
positive-integer weights at RIB:1617.

[f11] `bart2(family = "multinomial")` returns a `dbartsSampler` that is only a
HOST SHELL - it carries the design and priors but not the model, which lives on
the separate bartcore handle `$bc` (bart.R:1258, 1336). Every R5 mutator is
refused by `refuseHostMutation` (dbarts.R:770-783). At the handle level:
`setResponse` R RIB:2608 and `setOffset` R RIB:2605 (both redirect to the counts
channel), `setWeights` R RIB:2615, `setSigma` R RIB:2737, `setTestOffset` R
RIB:2380 (a flat offset leaves the simplex), `setPredictor` S RIB:4612,
`setTestPredictor` S RIB:4307. Its own channels - `setCounts` RIB:3509,
`setCategoryOffset` RIB:3587, `setCategoryTestOffset` RIB:3626 - are S at the
bridge but unexported (`dbarts:::`) and absent from the R5 object.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction at
dbarts.R:395 and `bart2Hurdle` (bart.R:1958) composes two ordinary `bart2()`
fits - an occupancy probit (bart.R:1990) and a lognormal positive part
(bart.R:1999) - glued at report time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f13] Grouped samplers fix the response at creation (RIB:4214). The engine
method exists and delegates faithfully (MOD:4517), so this is an unbuilt door,
not a model refusal - but it errors, so callers see a refusal.

[f14] Reads off the BASE family: grouped gaussian takes `setWeights` (MOD:4542)
and `setSigma`; grouped probit is refused on both (RIB:1642, RIB:2737); grouped
aft takes `setSigma` and refuses `setWeights`.

[f15] Arc `latent-subset-mask` (TODO:155), design FINAL, S0 and S1 LANDED;
artifacts .claude/latent-subset-mask-design/. A first-class 0/1
`setActiveRows` channel each family composes into its own precision vector,
with the latent draw skipped for inactive rows. Slices: **S0** pins (no engine
change); **S1** the channel plus gaussian, Student-t, probit, ordinal; **S2**
logistic, nbinom, aft; **S3** multinomial (global only); **S4** surface,
records, baselines. S0 landed at dc11a805 (the pins, now
inst/tinytest/test-active-rows-pins.R). S1 landed at 6db22aee: the engine
channel - `Chain::setActiveRows` CH:1282, which owns the single validating and
normalizing scan, `Sampler` SAM:1178, the facade's pure virtual FAC:278 and its
shape probe FAC:82 - plus gaussian, Student-t, probit and ordinal, the R5
`$setActiveRows` (dbarts.R:1109) and the bridge entry (RIB:3735). S2, S3 and S4
remain; a family S1 does not build is refused BY NAME at RIB:3735 (measured:
logistic, nbinom and aft each error there).

[f16] Arc `nameable-calibration` (TODO:279), design AMENDED FINAL, SCHEDULED
pre-release; artifacts .claude/nameable-calibration-design/. Names the
per-forest prior ANCHOR (`prior.scale`, the forest-total prior scale at k = 1,
in response units) rather than an sd, with a `$getCalibration` /
`$setCalibration` pair. Slices: **S0** signature freeze, LANDED 4c866286;
**S1** creation half, LANDED c2a7e89b; **S2** mid-chain get/set, remaining;
**S3** the flat-C half, executed inside the dbarts.h reshape's S1. S1 names the
model's `prior.scale` slot (A_class.R:398), resolved from `node.prior`'s
`scale =` / `sd =` spelling at `dbartsSpec()` (spec.R:251, `resolvePriorScale`
in R/model.R), and converted against the response transform by a private
`Chain::resolvedNodeScale` helper (CH:3350) shared by the single-forest
constructor and every `setModel` reinstall, so a round trip through the model
SEXP no longer reverts a named calibration; refused under a `k` hyperprior (the
`sd` spelling only, since a sampled `k` has no single value to divide by) and
for BCF/multinomial forests (three new refusal sites, see [f23]). Shipped
tests: inst/tinytest/test-calibration-creation.R (two composed probit arms at
construction ranges 16x apart agree to 1e-12 under a shared name, against 8.6
and 2.5 unnamed) and inst/tinytest/test-calibration-prior-draws.R (what the
named quantity means per leaf model - exact for the constant leaf, bounded
inequalities for linear/gp/monotone - across nine family and decoration
paths). What remains for a `P` cell here to flip to `S`: S2's mid-chain
`getCalibration` / `setCalibration` pair.

[f17] Zero weights are accepted, not refused (A_class.R:568-572 errors only
below zero and warns that zeros are ignored; bridge RIB:4462). The conditionals
are exact - leaf suffstats multiply by `w` (MOD:313, 1121), and the sigma
posterior counts only positive-weight rows (`numPositiveWeights_` MOD:2709,
recounted on every install at MOD:2899, consumed MOD:2723-2727). The one named
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
`ResponseModel::latents()` at its nullptr default (MOD:2616), and the bridge
returns `R_NilValue` (RIB:5502).

[f19] A Student-t fit records `family = "gaussian"`, so
`extract(type = "loglik")` takes the gaussian branch and evaluates `dnorm`
against the scalar `object$sigma` (generics.R:48-53). The t marginal is not
distinguished and the per-row `lambda_i` is not stored on the fit. No guard, no
doc, no TODO covers this - flagged, unadjudicated.

[f20] Weights on logistic are PG copy counts and a zero count is refused by
name at creation ("drop zero-count rows", RIB:1617-1622; R mirror
spec.R:131-140), so zero-weight subsetting is foreclosed for this family by the
weight semantics themselves - it is exactly the hole `latent-subset-mask` S2
fills.

[f21] Per-forest masking is REFUSED on model grounds and lands with that reason
in S3: the softmax margin is a log-sum-exp over the other K-1 forests, so a row
absent from category k's forest is still in every other category's likelihood
and "row i is out of category k only" restricts no likelihood. Only a GLOBAL
mask is well-posed here.

[f22] Multinomial's omegas live in the combiner, not the response model, and
`MultinomialResponse` does not override `latents()` (MOD:3563-3627), so
`getLatents` returns NULL. No accessor exposes them.

[f23] A named `prior.scale` is refused for BCF and multinomial forests AT
CREATION, by design - their per-forest leaf scales come from a calibration map
that owns them (map at COM:285-287), so a named value has nowhere to land.
Three refusal sites shipped at c2a7e89b: R-side `dbartsSpec()`'s BCF
composition (spec.R:408), the engine's own BCF-composition gate
(`refuseUnsupportedBCFComposition`, RIB:2299), and the multinomial forest
builder (`buildMultinomialSampler`, RIB:3032). The mid-chain `setCalibration`
(S2, not yet built) is designed to refuse the same way (nameable-calibration
synthesis 2.6 item 3), so these cells are expected to stay `R` rather than flip
to `S` when S2 ships.

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and errors at
generics.R:86-91, but each component fit (`$occupancy` probit, `$positive`
gaussian) supports `extract(type = "loglik")` on its own.

[f26] WORKS but is UNCOVERED, hence `?`. A BCF sampler's response IS a
`GaussianResponse`, so nothing on the path gates on the coupling: the shape
probe (FAC:82) reports `supportsActiveRows`, the mask composes into the case
weights at MOD:2794 (so the sigma df is inherited) and then into the per-forest
weights at CH:1101. MEASURED at 6db22aee on a 200-row two-forest sampler:
`$setActiveRows(a)` and the bridge `bartcoreSetActiveRows` are both accepted;
on a sampler carrying `w` the mask is BITWISE `setWeights(w * a)` in `train`
and in `sigma`; an all-zeros mask runs finite; a fractional element is refused.
No tinytest and no man/ entry exercises BCF under a mask, so no test, doc or
adjudication backs the cell. A per-forest mask is refused as REDUNDANT rather
than unbuilt: `setForestWeights` (RIB:3679) already expresses it - though note
that channel is deliberately NOT row removal (CH:927-946: it does not remove
the row from occupancy, the combination or the sigma df; it DOES reach that
forest's empty-leaf veto, which counts positive composed weights), and it is
MISSING from the R5 object, reachable only through the unexported
`bartcoreSetForestWeights` (bartcore.R:999).

[f27] Delegating / decorating, and that is what 6db22aee landed for the
active-rows column: neither row took an edit of its own. `GroupedResponse`
forwards `setActiveRows` to its base (MOD:4553) exactly as it forwards
`setWeights` (MOD:4542), advertising the base's capability (MOD:4550), and
`drawGroupEffects` already weights its per-group sums by `workingWeights()`
(MOD:4478), so an inactive row leaves its group's mean and precision and an
all-inactive group falls back to its prior through the same formula. The
heteroscedastic `formMeanWeights` (CH:3535-3542) reads
`response_->workingWeights()` at CH:3538 - the COMPOSED `w * a` while a mask is
installed - and divides by `s^2(x_i)`, so a zero stays a zero. Both are
PINNED: grouped at tests/cpp/test_sampler.cpp:1590 (an entirely inactive group
draws its effect from the prior, finite), heteroscedastic at
inst/tinytest/test-active-rows-pins.R:226-238. Both pins are FINITENESS-only; a
bitwise masked-vs-`setWeights(w * a)` strengthening for the two decorations is
slated for mask S2. The same delegation carries the nameable-calibration
CREATION half for grouped: `GroupedResponse::fitScale()` forwards to its base
(MOD:4584), so `Chain::resolvedNodeScale` (CH:3350) converts a named
`prior.scale` exactly as it does for the undecorated family, with no edit of
its own, and grouped is one of the nine family/decoration paths c2a7e89b's own
test measures (inst/tinytest/test-calibration-prior-draws.R:221, MEASURED
0.74210 vs a 0.75 target).

[f28] A heteroscedastic fit also records `family = "gaussian"` and takes the
same gaussian loglik branch with the SCALAR `object$sigma`, ignoring the
per-observation `s.train` surface the fit stores separately (bart.R:210-226).
Same standing as [f19]: no guard, no doc, unadjudicated.

[f29] The variance forest is a separate leaf model entirely, outside
`forests_`, and is not addressable by the mid-chain `setCalibration` (S2, not
yet built); recorded as a door in the design for when S2 ships (synthesis 2.6
item 7). The MEAN forest's own creation-time calibration is NOT gated by that
door - `Chain::resolvedNodeScale` (CH:3350) runs at forest.leaf.scale
assignment (CH:569-571) before the variance-forest branch (CH:660) and reads no
family flag - but it is UNCOVERED: c2a7e89b's nine-path test
(test-calibration-prior-draws.R) has no heteroscedastic arm. MEASURED here on a
`variance = ~x1` sampler: a named `scale = 1.5` prior-draw sd of 0.726 against
a 0.75 target (ratio 0.97, matching the covered families' band) versus 2.20
unnamed. No refusal site, no test, no doc backs the cell, hence `?` rather than
`P`.

[f30] Constructs today with no refusal site and no test. `spec.R:337-345`
refuses a variance forest only for a non-gaussian family or monotone
constraints, and a Student-t or grouped sampler still reports family
`gaussian`, so `resid.dist = student()` + `variance =` and grouped + `variance =`
both build (CH:562 decorates before CH:660 builds the variance forest). Whether
two scale mixtures on the same precision channel is a model anyone wants is not
adjudicated anywhere.

[f31] Recorded but UNBUILT doors, refused with that reason in the comment:
grouped ordinal because the cutpoint block and the group block are not yet shown
to interleave (RIB:2777-2781, ordinal.md section 8), grouped nbinom the same for
the dispersion block (RIB:2782-2786, negative-binomial.md section 7).

[f32] No surface at all: `applyGroupAttribute` (RIB:1969) is called from exactly
one site, RIB:2775, on the single-forest holder path, so `bartcore.groups` is
never read for a multinomial sampler.

[f33] DEFECT, not a refusal. `buildMultinomialForest` hard-sets
`forest.useDart = false` (CH:4448) while `bart2Multinomial` accepts and threads
`dart` through (bart.R:1219, 1235) and `buildMultinomialSampler` (RIB:3017-3059)
refuses nothing - unlike BCF, which names the option it drops (spec.R:393,
CH:4383). A user passing `dart = TRUE` gets a non-DART model silently.
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
SAM:717, and `growForestFromRoot` loops every forest (CH:1525). Neither is
exercised by a BCF test, and BCF has no `bart2()` surface, so both are reached
only through the R5 `$installTrees` (dbarts.R:1420) / `$growFromRoot`
(dbarts.R:841).

[f37] `rbart_vi()` carries no `warm.start` or `n.grow.sweeps` formal
(rbart.R:45-51) and its unknown-argument check (rbart.R:60-69) rejects them. The
underlying R5 sampler carries no group gate on either path, so this is a surface
gap, not an engine one.

[f38] The MEAN forest keeps DART; the variance forest never takes it
(`buildVarianceForest` CH:3499 never sets `useDart`, default false at CH:125).

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
(generics.R:86-91). Grouped ordinal is a recorded unbuilt door (RIB:2780).
`warm.start` and `n.grow.sweeps` unbuilt for the arc (bart.R:493). Its selecting
control attribute is undocumented in the shipped header ([f3]). One dedicated
tinytest file. SBC gamma3 resolved but not re-run at full R.

**nbinom.** As ordinal, plus: no `bart()` token, no `xbart()` token, no
`rbart_vi()` token. Grouped nbinom a recorded unbuilt door (RIB:2785). Pointwise
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
([f26]). The active-row mask is ACCEPTED on a BCF sampler and measures bitwise
`setWeights(w * a)`, but no test or doc covers BCF under a mask - the one `?`
the mask S1 landing leaves ([f26]). No pointwise loglik. Warm start and
grow-from-root are unrefused and untested for two forests ([f36]). Whole-data
`setData` stays undesigned (door 1 of the model-space survey).

**grouped.** `rbart_vi()`-only surface ([f8]): no `bart2()`, `dbarts()`,
`xbart()` or `dbartsSpec()` reach, and no `warm.start` / `n.grow.sweeps` formals
([f37]) though the engine paths carry no group gate. `setResponse` is an unbuilt
door (RIB:4214). Composition with a variance forest constructs unrefused and
untested ([f30]).

**heteroscedastic.** No `xbart()` reach. Pointwise loglik ignores the
per-observation `s(x)` surface it stores ([f28]). Selecting attribute
undocumented in the header ([f3]). Out of the SBC matrix, deferred not blocked
([f47]). One recorded unbuilt door from its own arc: the `setState` variance
column-mask gap (variance-forest-mutation-routing.md:499-500). The mean
forest's nameable-calibration creation half constructs and MEASURES correctly
but ships no test, doc, or man/ entry ([f29]) - the variance forest's own
calibration is a separate, permanent door (S2's `setCalibration` will not
reach it either).

**Cross-cutting.** The `latent-subset-mask` and `nameable-calibration` arcs
account for every `P` cell in table 3 and are the two largest scheduled
additions to this matrix; both are now half built - `latent-subset-mask` (S0
dc11a805, S1 6db22aee - what remains is S2, S3 and S4) and `nameable-calibration`
(S0 4c866286, S1 c2a7e89b - what remains is S2 and S3) - and both are sequenced
before the dbarts.h reshape's S1.
