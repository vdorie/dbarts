# Heteroscedastic AFT: lifting the variance-forest refusal

Status: LANDED, 2026-09-06 (2468e76f)

Admits `variance =` under `family = "aft"`, giving log-normal AFT a covariate-dependent
dispersion s(x) beside its mean surface. The refusal was R-side
(["a variance forest requires family"](../../R/spec.R)) and engine-side ([`varianceForestIsRefused`](../../src/bartcore/facade.hpp)
and the `family == ResponseFamily::gaussian` guard before `buildVarianceForest` in the single-forest
Chain constructor); the reason both give does not hold for aft.

## 1. Why the stated reason is false for aft

"Routes precision through its own latent channel instead" holds for probit, logistic, ordinal and
nbinom, whose `workingWeights()` carries their own latent precisions - the channel the variance
forest divides into ([5. Decision (fork 4) - per-observation residual variance channel](heteroscedastic.md#5-decision-fork-4---per-observation-residual-variance-channel)).
Not for aft: [`AFTResponse`](../../src/bartcore/model.hpp) builds its contained [`GaussianResponse`](../../src/bartcore/model.hpp) with a NULL weight
pointer, leaves `workingWeightsVaryPerSweep()` at the base false and delegates `drawSigma`, so
`workingWeights()` is null - or, under a mask, the composite a masked heteroscedastic gaussian
carries. The channel is free. The real collision is [`TResponse`](../../src/bartcore/model.hpp), whose Gaussian is built over
the composite mixing precisions and reports the varying flag; that refusal stands.

## 2. The model

    log T_i = f(x_i) + offset_i + s(x_i) eps_i,   eps_i ~ N(0, 1)

f is the constant-leaf mean forest, s^2 the multiplicative variance forest. An uncensored row is
gaussian data on log T; a censored one's latent log-time is redrawn each sweep from N(f(x_i) +
offset_i, s^2(x_i)) truncated below at its log censoring time - the shipped draw at that row's own
scale rather than the shared sigma.

Coherent as a Gibbs cycle over three blocks. Given (f, s) each censored latent is a conditional
normal draw at its own scale truncated at its own bound, the rows conditionally independent, so a
per-observation sd is the exact conditional. Given the latents the augmented log-times are complete
gaussian data and the variance forest backfits against e_i = y_i - f(x_i) as under gaussian
([`Chain::sweepVarianceForest`](../../src/bartcore/chain.hpp)), its conjugacy assuming nothing about the provenance of y. The
mean forest sees w_i / s^2(x_i) ([`Chain::formMeanWeights`](../../src/bartcore/chain.hpp)), aft's free channel supplying 1 or
the mask, and [`Chain::run`](../../src/bartcore/chain.hpp) already orders the blocks so: mean weights, mean forest,
`refreshLatents`, no sigma draw, variance forest.

## 3. Engine interface for the per-observation variance

[`Chain::buildVarianceForest`](../../src/bartcore/chain.hpp) pins the scalar sigma at 1 on the working scale, leaving
[`VarianceForest::combinedVariance`](../../src/bartcore/chain.hpp) the residual variance there. `AFTResponse` reads a scalar
sigma in three places - [`AFTResponse::refreshLatents`](../../src/bartcore/model.hpp),
[`AFTResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) and [`AFTResponse::setResponse`](../../src/bartcore/model.hpp), which redraws the
latents itself - all needing the vector.

**A, a parameter:** `const double* variance`, null when homoscedastic, on those three virtuals -
explicit and unstaleable, at three signatures across the base and eight response models plus some 30
tests/cpp call sites. **B, an installed pointer:** `ResponseModel::setVarianceSurface(const
double*)`, default no-op, overridden by `AFTResponse`, installed where `combinedVariance` is
allocated - one virtual, one override, two Chain call sites, no test edits.

Staleness decides, and is in hand: the `VarianceForest` is constructed once and never reset, swapped
or moved, chains are held by `unique_ptr`, every rollback writes `combinedVariance` elementwise, and
only `VarianceForest::initialize` and [`Chain::resizeVarianceStorage`](../../src/bartcore/chain.hpp) reallocate - the second
unreachable for aft, running only from the whole-data arm [`bartcore_setData`](../../src/R_interface_bartcore.cpp) refuses.
RECOMMEND B: one named Chain helper called from both allocation sites, with a component test
asserting the pointer equals `varianceFits()` after each. It owns the TRAIN vector only,
`combinedVarianceTest` reallocating in [`Chain::resizeTestStorage`](../../src/bartcore/chain.hpp) where no latent draw reads
it.

RESTORE CONTRACT, in the shape [`NBResponse`](../../src/bartcore/model.hpp)'s dispersion states: [`Chain::setState`](../../src/bartcore/chain.hpp)
restores the latents BEFORE rebuilding the surface, safe only because `AFTResponse::restoreLatents`
is a memcpy plus a working rebuild reading neither sigma nor surface. Keep `restoreLatents`
surface-free, or move the rebuild ahead of it.

Call sites: `Chain::run`'s `refreshLatents`, [`Chain::growForestFromRoot`](../../src/bartcore/chain.hpp)'s,
[`Chain::setResponse`](../../src/bartcore/chain.hpp) through `AFTResponse::setResponse`'s redraw, and
[`Chain::storeSample`](../../src/bartcore/chain.hpp)'s `computeLogLikelihood`; A must also edit
[`LogisticResponse::setWeights`](../../src/bartcore/model.hpp), a fifth `refreshLatents` caller outside any `setResponse`
body. LANDING NOTE: a new `ResponseModel` virtual is a full-recompile hazard - `--preclean`, or
stale objects bus-error.

FINDING, independent of this proposal and load-bearing for it.
[`GaussianResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) takes the pinned sigma of 1 and divides only by the
USER weights, so under a variance forest it reports a density at the response range, not at s(x_i)
times it. Reachable only through the flat C API's `logLikelihood` results member - R's `extract(type
= "loglik")` recomputes off `s.train` and is right - and one variance-aware log-likelihood repairs
gaussian and aft together.

## 4. Factory and chain gates

`varianceForestIsRefused`'s `family != gaussian` becomes `family != gaussian && family != aft`.
Every factory that CAN build one asks it ([`createSampler`](../../src/bartcore/facade.hpp), [`createSamplerOverStore`](../../src/bartcore/facade.hpp));
the amplitude and multinomial factories carry their own bare `numVarianceTrees > 0` refusals
([`createAmplitudeSampler`](../../src/bartcore/facade.hpp), [`createMultinomialSampler`](../../src/bartcore/facade.hpp)) and stay refusing. The Chain
constructor's `family == ResponseFamily::gaussian && options.numVarianceTrees > 0` admits aft too,
both arms inside the constant-leaf `if constexpr` guard.

Three texts change: spec.R's message; the flat header's `"bartcore.variance"` comment, "Gaussian
constant-leaf models only" (["bartcore.variance"](../../inst/include/dbarts/dbarts.h)); and the bridge's null-factory message
(["variance forest is combined with a family other than"](../../src/R_interface_bartcore.cpp)). Only the third reaches a
`LinkingTo` consumer, `applyVarianceAttributes` parsing the attribute with no family test, so the
header comment is load-bearing.

Refusals that stay: `sigmaIsPinned` is `hasVarianceForest || (family != gaussian && family != aft)`,
so a heteroscedastic aft loses `setSigma` as a heteroscedastic gaussian does,
[`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) naming the variance forest before the family;
[`refuseVarianceForestScaleUpdate`](../../src/R_interface_bartcore.cpp) is family-blind, so `setResponse`/`setOffset` are taken only
at `updateScale = FALSE`; `setWeights` stays refused ([`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp)), the
declined user channel being what frees the internal one; `setData` likewise.

Four channels change meaning without changing code. [`Chain::setModel`](../../src/bartcore/chain.hpp) under a variance forest
skips its `gaussian || aft` sigma clause and recalibrates the scale leaf: DECIDE that aft FOLLOWS
the heteroscedastic gaussian, one residual prior addressing one scale leaf, so a family split would
give one object two laws; it inherits, not widens, the gap that comment records. A warm start
([`Chain::installVarianceForest`](../../src/bartcore/chain.hpp), `installForest`) redraws no latents, so destination latents
stand under the donor's surface until the next sweep, and `growForestFromRoot` never sweeps the
variance forest, so its redraws run against the constant initial surface - the staleness
[`AFTResponse::setOffset`](../../src/bartcore/model.hpp) documents. The active-rows mask composes, `sweepVarianceForest`
handing the masked composite to the scale leaf, whose [`ConstantVarianceLeaf::accumulate`](../../src/bartcore/model.hpp) drops
non-positive weights from n and ssr alike; [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp) has no variance-forest
gate, so that composition arrives with the lift, its residue an inactive censored row's stale latent
entering `vf.meanResidual` before its zero weight annihilates it.

## 5. R surface

No packaging change: `s.train`/`s.test` are built above the packager's binary branch and attached in
the non-binary one, keyed on `control@binary`, which [`isBinaryFamily`](../../R/spec.R) sets for probit and
logistic alone. The fit's `$sigma` becomes the pinned constant carrying no posterior content;
[`resolveDrawsVars`](../../R/diagnostics.R) swaps `"sigma"` for `"mean.s"` on any fit carrying `s.train`,
but only on the draws-array path ([`presentDrawsVars`](../../R/diagnostics.R)), so
[`plotSigmaTrace`](../../R/plot.R), gated on `"sigma" %in% names(x)` alone, gives a heteroscedastic aft the
flat constant trace a heteroscedastic gaussian already gets - pre-existing.

Two functions read `object$sigma` where they must read the surface.
[`pointwiseLogLikelihood`](../../R/generics.R)'s aft branch must take [`heteroscedasticScale`](../../R/generics.R) of
`s.train` when present, with the gaussian branch's length check.
[`survivalProbabilitiesFromDraws`](../../R/bart.R) recycles one sigma per draw across observations where the
scale is per draw AND per observation; training rows always have `s.train`, and at `newdata` the
choice is to refuse, mirroring the ppd wording, or read the replay, [`predict.bart`](../../R/generics.R)
parking `sqrt(variance)` on `attr(result, "s")` for every type including the `"bart"` one
[`survivalProbabilities.bart`](../../R/bart.R) asks for. RECOMMEND the replay where available, the refusal
where not - it needs saved trees.

[`dbartsDrawLatents`](../../R/augmentation.R), the exported replay of `AFTResponse::refreshLatents`, takes
`sigma` through [`augScalar`](../../R/augmentation.R) as one positive scalar. DECIDE it stays scalar and
REFUSES a length-n sigma by a named message, because the flat header's [`dbarts_drawLatents`](../../inst/include/dbarts/dbarts.h)
takes `double sigma`: widening the R helper alone forks the two replays, widening both moves a
signature on the only shipped header - a priced door; per-row calls are the workaround.

Nothing else moves: the ppd branch already hands [`sampleFromPPD`](../../R/generics.R) the same
`heteroscedasticScale(s)` with no family test, `variance` is a formal of both [`bart2`](../../R/bart.R) and
[`dbarts`](../../R/dbarts.R) with `"aft"` on the ordinary single-forest route, and
`checkFamilyUnsupportedArgs` never gated aft.

## 6. Gates

**(a) Reduction, bitwise.** An all-uncensored heteroscedastic aft fit must be bit-identical to a
heteroscedastic gaussian fit on log T at the same seed, `AFTResponse::refreshLatents` returning
before drawing and every other hook delegating: [`testAFTReduction`](../../tests/cpp/test_model.cpp) with a
variance forest on both arms. It pins RNG-stream equality, reaching no truncated draw.

**(b) Per-observation redraw.** Extend [`testAFTCensoredMoments`](../../tests/cpp/test_model.cpp): one
`AFTResponse`, a fixed two-level variance vector whose censored rows differ by a large factor, many
redraws, each row's empirical mean and sd against the analytic lower-truncated normal at THAT row's
sd, and a poison arm substituting the mean of the two sds that must fail. It separates
per-observation from per-fit and catches the censored-index/row-index confusion. Its UNEXTENDED body
pins the implementation: the homoscedastic path must keep the literal `sigma * scale` expression
when the variance pointer is null, so no family or leaf gains or loses a draw and every RNG-locked
baseline stays valid.

**(c) Wiring, end to end.** The only gate testing surface-versus-pinned-1 and
working-versus-original scale Chain-wide: make the surface degenerate at a known scalar and reduce
to benchmarks/R/aft-exact.R's enumeration - one variance tree with [`chisq`](../../R/model.R) at a large
`df` and `data@sigma` anchored at the known sigma, so the scale leaf's posterior is prior-dominated.
NOT `resid.prior = fixed()`, which under a variance forest fixes nothing: the bridge's fixed arm
sets only `sigmaIsFixed` and `fixedSigmaSq`, leaving `sigmaDf`/`sigmaRawScale` at their defaults, so
`buildVarianceForest` still calibrates a df-3 leaf the data dominates. On that fixture the mean
censored latent moves from 0.221 at the correct internal sd to 0.843 at a pinned 1, against a
tolerance of 0.012; it discriminates no per-observation error.

**(d) Latent PIT, on the recovery leg.** Simulate a two-level true s(x) under heavy censoring;
assert the posterior surface separates the two levels and the mean surface stays unbiased, and on
the same fit, after burn-in, form for each censored row and recorded sweep

    v = (Phi((z - mu)/s) - Phi((b - mu)/s)) / (1 - Phi((b - mu)/s))

with z the drawn latent, b its bound, mu the current fit and s = sqrt(varianceFits()[i]) on the
internal scale the draw used. Pool v by x-cell: under the correct per-row scale each cell's v are
U(0,1), so a per-cell `ks.test(v, "punif")` passes - the idiom ["ks.test"](../../R/validateComposition.R)
already uses. Poison arm, the substitution (b) makes: drive the redraw at a scale constant across
rows and the low-s cell's v pile at 0 while the high-s cell's pile at 1, which the per-cell KS must
reject. It is the only gate sensitive to a wrong per-observation scale INSTALLED BY THE CHAIN, and
needs no oracle. HONEST GAP after all four: the joint calibration of (mean forest, variance forest,
censored latents) is still not SBC-tested, aft being out of that matrix until a censoring-status
setter lands ([Decision - scope](../plans/sbc-family-tiers.md#decision---scope)).

**(e) Matrix cells.** In [4. Composition rules](feature-matrix.md#4-composition-rules) the variance
forest's family rule names gaussian or aft, the four latent families refused for owning the weight
channel. In [1. Structural signature](feature-matrix.md#1-structural-signature) aft's sigma and unit-scale
cells become conditional on the variance forest, and the hetero row's case-weights and latents
cells become by-family (its footnote states the rule); the aft status-setter gap now names
heteroscedastic aft too. benchmarks/R/composition-matrix.R needs no code change: it probes only S
cells, derives the aft variance-forest probe from the matrix itself, and its base fixture already
threads an extra `variance =` into the aft recipe (["extra:variance"](../../benchmarks/R/composition-matrix.R)).

## 7. Consumers, and what stays out

Neither downstream composes both: stan4bart on bartcore maps an `"aft"` token in `getBARTFamily` but
fits only gaussian and probit, with no variance-forest reference in R or src, and bartCause on
dbarts-1.0 has neither. Nothing to migrate, no lockstep release constraint.

Out of scope: heteroscedastic probit, logistic, ordinal and nbinom, the latent-channel collision
being real for all four ([11. Out of scope, and the doors](heteroscedastic.md#11-out-of-scope-and-the-doors)); hazard,
probit over person-period rows, by inheritance; Student-t with a variance forest, `TResponse` owning
the weight channel; the censoring-status setter and with it aft's SBC arm; a vector `sigma` on the
two `drawLatents` surfaces; left and interval censoring and competing risks
([Out of scope (v1)](survival.md#out-of-scope-v1)); non-constant variance leaves; heteroscedastic BCF
and multinomial.

## 8. Alternatives and recommendation

**A, correct the message and keep the refusal:** ten lines in one file, one tinytest, one footnote;
no gates, no engine change, nothing a user can fit. **B, lift it:** engine, two predicates plus one
virtual with one override plus the variance-aware log-likelihood that also repairs the gaussian cell
section 3 records; R, two functions and one named `drawLatents` refusal; docs, 6(e)'s cells, one
header comment, one bridge message; gates (a) through (e), beside roughly two hundred and fifty
lines of engine and R.

RECOMMEND B. The value is nameable: covariate-dependent dispersion of log survival time - the spread
of the survival distribution varying with x, not only its location - and predictive survival curves
whose width is itself estimated. The reference suite ships AFT and heteroscedastic gaussian as
separate models and composes neither ([R packages: engines](bart-landscape.md#r-packages-engines)). Price
does not rank the item, and B is cheap: the weight channel is free, the scale leaf conjugate, the
only new code passing one vector where a scalar goes today.

## 9. Landing

Engine: `ResponseModel::setVarianceSurface` (default no-op) with `GaussianResponse` and
`AFTResponse` overrides; [`Chain::installVarianceSurface`](../../src/bartcore/chain.hpp) from both allocation points
([`Chain::buildVarianceForest`](../../src/bartcore/chain.hpp), [`Chain::resizeVarianceStorage`](../../src/bartcore/chain.hpp));
[`varianceForestIsRefused`](../../src/bartcore/facade.hpp) and the single-forest Chain constructor admit aft. The
homoscedastic arm of every reader keeps the literal `sigma * scale` expression.

Gates, this host, arm64/macOS.

(a) Reduction, bitwise: [`testAFTReduction`](../../tests/cpp/test_model.cpp), fits and variance surface
both, with the surface asserted non-constant. Also at the R level
(["an uncensored heteroscedastic aft IS the gaussian fit"](../../inst/tinytest/test-aft-heteroscedastic.R)).

(b) Per-observation redraw: [`testAFTCensoredMoments`](../../tests/cpp/test_model.cpp), two censored rows
at sds 0.25 and 1.0 under an identical bound-minus-mean gap, 60000 redraws, mean and sd each within
0.03 of the lower-truncated normal at that row's sd, with the pooled sd kept as a negative
expectation. Poison (row-constant scale, a temporary mutation of the redraw): both per-row
assertions fail, mean 0.922 against 1.498 and sd 0.106 against 0.551.

(c) Wiring: benchmarks/R/aft-exact.R, second arm. One variance tree, chisq at df 1e6, data sigma
anchored at 0.6. Quick: max gap 0.0020 against a tolerance of 0.020, surface 0.5995 against 0.6000.
Full: 0.0019 against 0.012, surface 0.5995.

(d) Latent PIT: benchmarks/R/aft-hetero-pit.R, per-cell `ks.test(v, "punif")` at alpha 1e-5. Quick:
p 0.922 over 9900 draws and 0.546 over 8280; surface 0.325 and 1.110 against a truth of 0.30 and
1.20, mean-surface bias -0.003. Full: p 0.804 over 33000 and 0.561 over 27600. Poison (the same
row-constant mutation): p 0.048 and 0, cell 2's mean PIT 0.279 against 0.5, and the recovery leg
fails alongside it at a surface ratio of 2.31 and a bias of -0.151. Added to the exact-gates
workflow; its quick mode is 0.4 s.

Suite: tests/cpp 277 checks, 0 failures, clean under -fsanitize=address,undefined; tinytest 7424
passes, 0 failures, from 7393 before. Equivalence bitwise identical on 50 gaussian, 12 bcf and 11
multinomial scenarios - no existing family or leaf gains or loses a draw.

