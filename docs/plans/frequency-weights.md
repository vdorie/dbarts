# frequency-weights: weights read as counts of identical rows beside the precision default, and the binomial response

Status: PLANNED 2026-11-03 on main (3e4c002d8); amended the same day after a blind critique (scratch/fwp-crit/report.md, probes
p1 to p7 there), verdict BUILD AFTER AMENDMENTS, rebased onto main at 5a99172e6 (forest-door landed at a6a16c31a), and onto
main at 7937d252f on 2026-11-05; amended 2026-11-05 for the work planned or built since
([Amendments of 2026-11-05](#amendments-of-2026-11-05)): latents-first, hurdle-completion, approximate-draws,
coefficient-prior-live and host-neutral-engine's landed S1. S0 landed on main (1a23f9a60). The feature is not built.
Rulings since, read whole for that pass: dec-B586 (a sweep draws its augmentation variables first; no setter draws them),
dec-B588 (a hurdle fit takes weights of 0 and 1 on both parts), dec-B589 (setData keeps the drawn family parameters), dec-B583
(an accept-or-keep latent is exact). Rulings, quoted where they settle a call: dec-B125 (frequency weights join precision weights by explicit choice; precision
stays the default), dec-A11 (sigma's count under precision), dec-B563 (the gaussian-type and Polya-Gamma families take
frequency weights; "Let's do (b) plus the Polya-Gamma families."), dec-B564 (weights.type, "I like (a)."), dec-B565 (a
frequency row's predictive draw is one observation), dec-B566 (a count family's training draw is one observation; dec-B470's
replace clause revised), dec-B567 (the binomial response; "Yes, have the binomial added. It can look like similar R
interfaces."), dec-B568 (real weights where the draw is exact), dec-B569 (membership as 0/1 weights on the data object),
dec-B570 (a zero weight alike on every family), dec-B580 (a 0 row's predictive draws NA on every family, correcting
dec-B570's premise; "0 weight means out of sample"), dec-B573 (a
proportion with the trials as weights refused; "Use (a)."), dec-B574 (probabilities for expected values, counts for draws;
"Use (a)."), with dec-B133 and dec-B177 (zero-trial rows), dec-B202 (the one-million cap), dec-B470 and dec-B472 (predict's
weights), dec-B494 (a feature merges once), dec-B85 (the header change frequency weights' field rides). The runs behind it:
the note's probes (scratch/fw-note, the first critique's p1 to p11 under scratch/fw-crit/probes), p12cost.R in scratch/fw-note
(a unit Polya-Gamma draw at 164 ns, psi 0, and 159 ns, psi 3, on the shipped build of c3d268c8f on the shared laptop), and the
plan critique's p1 to p7, rerun here on its library of main (scratch/libs/fwp-crit) and 0.9-34's (scratch/libs/cran), and the
amendment's probes in scratch/fw-amend/probes on a shipped build of main 7937d252f (scratch/libs/fw-amend). No prototype was
built here.

agent: S0 landed (one sonnet implementer and review, on main at 1a23f9a60). The feature: S1, S2, S3 and S4 one opus
implementer each (effort high), S1 and S3 the riskiest; S5 and S6 one sonnet implementer each (effort medium); S7 and S8
sonnet. One opus review per engine slice (S1 to S4) that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior; one sonnet review for S5, S6 and the gate slice S7
([Process by risk](README.md#process-by-risk)). Tier "Changes draws" for S1 to S4: new draws on new calls, engine read-offs
under frequency, and the creation and mask paths of every family rewritten, though no call taken today is meant to move;
tier "R surface" for S5 and S6, both pure R.
rng: neutral, bit for bit, for every engine draw of a call taken on the build it starts from, main after coefficient-prior-live
and latents-first (whose draws latents-first moved and re-recorded once): the equivalence trio, the four snapshot files,
tests/cpp and the exact gates do not move against that build's baselines. Reported values that move are listed under
[RNG class](#rng-class): a weighted logistic fit's training predictive draws (dec-B566), a 0 row's predictive draws on
probit, ordinal, hazard, hazard.probit and nbinom (NA, dec-B580), `$activeRows` after 0/1 weights (NULL), a logistic zero
weight (taken). New calls take new draws, each held to the expanded-data fit or, at predict, to its stated distribution.
window: the feature for 1.0-0 (dec-B567: "before 1.0-0"), one branch, frequency-weights, merged once (dec-B494). Build order,
the orchestrator's call as latents-first's plan states it: coefficient-prior-live, latents-first, this feature, then
host-neutral-engine S2a onward; hurdle-completion, R only, lands before this whenever green; see
[Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order).
budget: the feature ~4,300 lines, ~1,560 of them package code: engine ~570 (S1 ~285, S2 ~110, S3 ~125, S4 ~50), bridge ~200,
R ~765 (S4 ~80, S5 ~435, S6 ~250), standalone math ~25 written (S3; R's choose.c copied unedited, not counted); help ~190,
NEWS ~30; tests/cpp ~820, tinytest ~1,200, gate scripts ~450, records ~60. Agent time, estimated and not measured: S1 ~12.5, S2 ~7, S3 ~9.5, S4 ~6.5, S5 ~14, S6 ~9, S7
~6, S8 ~1, four engine reviews with mutants and fix rounds ~18, three R and gate reviews ~7.5, gates and sanitizers ~5,
rebases ~3: ~99. Plans here run 1.5 to 2 times low, so 150 to 200. The open call's recommendation is counted; the other
options change S5 by under an hour.

## Summary for the maintainer

You ruled that a weight can also mean a count: a row of weight 3 stands for three identical observations, and the fit is the
fit to the data with every row written out that many times. Precision stays the default, chosen with one new argument beside
the weights. This plan builds that, the membership and zero-weight rulings that came with it, and the binomial response.

What a user sees. bart, dbarts, dbartsSpec and xbart take weights.type. Left unset it means what weights mean today: a
precision on gaussian, Student-t and heteroscedastic fits, a count on logistic, negative binomial and multinomial fits.
"frequency" on the first three changes how many observations the noise level and the other data-derived defaults are read
from, and takes weights of any positive size. Negative binomial and multinomial fits take whole-number weights for the first
time. Probit and ordinal keep their 0/1 weights for membership, which now stay on the sampler's data object. A weight of 0
takes a row out of the likelihood on every family, logistic included, which refused it. A row's predictive draw is one new
observation like the row. A hurdle fit keeps the 0 and 1 you ruled for it, and treats weights.type as probit does.

The binomial response is written as R writes it, cbind(successes, failures), under family = binomial or "logistic", and
xbart takes it too. Its fitted values are probabilities and its predictive draws counts, as glm's are. A proportion with the
trial counts as weights is refused with a message showing the cbind form.

As you ruled after the critique corrected the zero-weight premise, a row of weight 0 has missing predictive draws on every
family; probit, ordinal, the probit hazard and negative binomial fits drew values there, and its fitted value is still
reported. You named predict(weights = 1) as the way to draw such a row as a new observation, but the ordinal, negative
binomial and multinomial fits' predict takes no weights today and refuses them by name. This plan adds the argument to all
three. 1 draws a new observation on each; on ordinal, whose weights mark rows in and out, 0 gives a missing draw; the
negative binomial reads a whole number as that many new observations summed and refuses 0, as logistic reads its weights as
trials and refuses 0. One choice is yours, below: what predict's weights mean where one draw is a category, on ordinal and
multinomial fits.

The xbart fix that was to land first has landed on main. Since this plan was written you also ruled that every sweep draws
its hidden values first, so no setter draws them; this plan is built after that change and is simpler for it: nothing here
redraws when weights, counts or rows change, the next sweep does. The approximate-draw work you opened is where the cost
below changes, and where fractional counts on logistic and binomial fits arrive; this plan keeps every draw exact and whole.

What it costs. Nothing per sweep on the gaussian-type families. On the count families a sweep grows with the total count, as
logistic's does now: each trial of each row is about 160 nanoseconds a sweep, so 1,000 rows of 1,000 trials take about 0.16
seconds a sweep.

What does not move. No posterior draw of any fit that runs on the build this starts from changes: the engine's random
numbers, the four recorded snapshot files and the equivalence baselines stay bit for bit, and each slice proves it. stan4bart,
bartCause and treatSens need no change; stan4bart can drop about 30 lines that become dead.

One open call, below. The feature is about 4,300 lines over eight slices, about 99 agent-hours as estimated, 150 to 200 by
the record of plans here.

## Goal

After it lands: a frequency-weighted fit on the gaussian-type families equals the fit to the expanded rows draw for draw; the
Polya-Gamma families take whole-number counts with shape times w, each drawn as its adjacent copies are; a binomial row is
one logistic row of shape n w drawn as the two-row expansion is; probit and ordinal membership rides the data object; a 0
means the same on every family; the R surface of predict, fitted, residuals, extract, the log-likelihood and xbart reads both
readings and the binomial as the rulings state, predict's weights reaching the ordinal, negative binomial and multinomial fits
as dec-B580's door asks; and no call taken on the build it starts from moves an engine draw.

## Context

- Design: docs/design/frequency-weights.md (its Calls made 1 to 14 are the design calls this plan builds; this plan's own are
  below), its section The two readings differ only in the scale's count, and its table What moves, place by place.
- Engine: [`GaussianResponse`](../../src/bartcore/model.hpp) (the positive-row count `numPositiveWeights_`, recounted by
  `installWeights`, and `computeLogLikelihood`),
  [`ChiSquaredScalePrior::drawSigmaSqFromPosterior`](../../src/bartcore/model.hpp),
  [`ConstantVarianceLeaf`](../../src/bartcore/model.hpp) (its `accumulate` counts positive rows),
  [`TResponse`](../../src/bartcore/model.hpp) and [`ResidualDfPrior::drawIndex`](../../src/bartcore/model.hpp),
  [`LogisticResponse`](../../src/bartcore/model.hpp) (shape from `lround(weights_[i])`, working response `w (y - 1/2) / omega -
  offset`, `computeLogLikelihood`), [`NBResponse`](../../src/bartcore/model.hpp) (its center from `computeShift` and
  `deriveShift`) and its `NBShapePrior`, the multinomial Polya-Gamma loop in
  [`MultinomialForestCombiner`](../../src/bartcore/combiner.hpp), [`ProbitResponse`](../../src/bartcore/model.hpp) and
  [`OrdinalResponse`](../../src/bartcore/model.hpp) (the active-row mask only), the leaf covariates' moments
  ([`standardizationMomentsForColumn`](../../src/bartcore/data.hpp)) and a forest's default sd
  ([`basisColumnScale`](../../src/bartcore/combiner.hpp)), with the setter doors that re-derive them (`updateStandardization`,
  `updateBasisScale`). The weights digest ([`SamplerBase::weightsDigest`](../../src/bartcore/facade.hpp)), the Student-t
  zero-weight record and setState's reconciliation retire under latents-first before this builds; so does every setter draw.
- The Polya-Gamma draw seam: [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp), the unit-draw loops in
  `LogisticResponse::refreshLatents`, `NBResponse::drawOmega` and the combiner, which approximate-draws' route replaces where
  selected.
- Math: [nmath.h](../../src/include/external/nmath.h), the six ext_ math names (Rmath on the R build, the dbarts_nmath_ copy
  under standalone/nmath otherwise), landed with host-neutral-engine's S1; R's `lchoose` (choose.c) is not among them.
- Bridge: [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp),
  [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp), [`refuseCountOverCap`](../../src/R_interface_bartcore.cpp),
  [`validateResponseSupport`](../../src/R_interface_bartcore.cpp), the creation refusals of ordinal, nbinom and aft weights
  beside them, the weights slot read in the data conversion.
- R: [`enforceWeightPolicy`](../../R/spec.R), [`isMaskedWeightFamily`](../../R/spec.R),
  [`refuseNonCountWeights`](../../R/spec.R), [`refuseNonBinaryResponse`](../../R/spec.R), [`dbartsSpec`](../../R/spec.R) and its
  `active` element, [`estimateStartingSigma`](../../R/spec.R) over [`residualStandardError`](../../R/utility.R) and
  [`sparseResidualStandardError`](../../R/utility.R), [`restateMatrixResponseError`](../../R/data.R) and the classed
  `dbartsMatrixResponseError`, [`sampleFromPPD`](../../R/generics.R), [`pointwiseLogLikelihood`](../../R/generics.R),
  [`predict.bart`](../../R/generics.R), [`predict.bartOrdinal`](../../R/generics.R),
  [`predict.bartNegbin`](../../R/generics.R) and [`predict.bartMultinomial`](../../R/generics.R) with the `weights` reason in
  `predictForeignReasons` that refuses them there, [`fitted.bart`](../../R/generics.R), [`residuals.bart`](../../R/generics.R),
  [`extract.bart`](../../R/generics.R), [`xbartLossFunction`](../../R/xbart.R), [`xbartSigestFunction`](../../R/xbart.R) and the
  0/1 refusal in [`xbart`](../../R/xbart.R), [`refuseCountsMutation`](../../R/bartcore.R), the dbartsModel class in
  [R/A_class.R](../../R/A_class.R), `$setModel`'s refusals in [R/fittingLevel.R](../../R/fittingLevel.R) (forest-door's, with
  coefficient-prior-live's coefficient refusal beside them), the hurdle arm of bart and its weights refusal (hurdle-completion's,
  on its branch), the
  sampler's `$setWeights`, `$setData`, `$setCounts`, `$setResponse` and `$setActiveRows` ([R/dbarts.R](../../R/dbarts.R)).
- Help: man/bart.Rd, man/dbarts.Rd, man/dbartsSpec.Rd, man/xbart.Rd, man/dbartsSampler-class.Rd (`setData`'s and
  `setWeights`' mask sentences), man/dbarts-embedding.Rd (its "survey weights" row), man/samplePriorPredictive.Rd (its
  weighted logistic sentence), man/predict.bart.Rd and the ordinal, negative binomial and multinomial predict pages. latents-first
  rewrites the sampler's setter sentences before this plan edits them.
- Tests that pin today's weight behavior: test-binary-weight-mask.R, test-weighted-logistic.R, test-weighted-binary-ppd.R,
  test-logistic-weight-swap.R, test-zero-weights.R, test-xbart-weight.R, test-pointwise-loglik.R, test-nbinom.R,
  test-xbart-loss.R; test-state-weight-pairing.R and test-state-zero-weight-rows.R as latents-first leaves them (the second
  retired, the first's digest arms gone).
- Plans it follows or meets: docs/plans/latents-first.md (Order and overlaps), docs/plans/hurdle-completion.md (Interaction with
  frequency-weights, and order), docs/plans/approximate-draws.md (Real shapes, Setters and marked rows),
  docs/plans/host-neutral-engine.md, each on its branch.
- TODO items: frequency weights (dec-B125), mask-on-own-data-object (closed by S4), hurdle-weight-extensions (its logistic zero
  part takes counts once this lands), approximate-distributions (the binomial's cost case, not this plan's); xbart-binary-rmse
  closed by S0.

## Constraints

- No draw of a call taken today moves: every engine path a precision fit, an unweighted fit or a binary logistic fit runs
  keeps its expression and its order of random numbers. Where a new term is zero under precision, it is not added, rather
  than added as 0.0.
- Every draw stays exact; real weights on the Polya-Gamma families are refused until a real-shape draw exists (dec-B568),
  which approximate-draws' S3 brings for logistic, binomial and hazard.logistic (below).
- Each new Polya-Gamma shape (nbinom's `w (y + r)`, the binomial's `w n`, multinomial's `w n_i`) is formed and handed to the
  unit-draw loop the family already has; no slice adds a loop of its own, so approximate-draws' route replaces one loop per
  family whichever lands first.
- No setter draws: under latents-first a weight, count or mask change stores what it is given and rebuilds what is
  deterministic, and the next sweep draws the augmentation variables first.
- No row out of the likelihood carries a zero precision or a non-finite working response: the node kernels read
  `omega` and `working / omega`, and 0 there is a NaN. latents-first's step 2 draws active rows only, so a row out of the
  likelihood keeps its creation start for the whole run.
- Out of scope: frequency weights on probit, ordinal, the probit hazard and aft (dec-B563); hurdle.lognormal's weights
  beyond 0 and 1 (dec-B588, hurdle-completion's; here only its weights.type row); the probit binomial (dec-B567); any
  approximate draw (TODO approximate-distributions); bartBT and rbart_vi, which take no weights.type; the shipped header (its
  field rides host-neutral-engine's S3, below).
- No new dependency; base R calls within DESCRIPTION's R floor.

## Design

### The reading and where it lives

- R: dbartsModel gains a slot `weights.type`, a character of length one, NA for the family's own reading (weights.type = NULL
  as given), else "precision" or "frequency". bart, dbarts, dbartsSpec and xbart take the argument and write the slot; the
  spec's resolution, where the family is known, writes the resolved value, so dbartsSpec's returned model carries it and
  creation reads it. A model made by hand with NA resolves at creation the same way.
- Bridge: reads the resolved slot into a sampler-level flag, `bartcore::WeightsReading { precision, frequency }`, passed to the
  response family's constructor beside the weights.
- Engine: the flag is held by the families that read it. A state carries no reading: after latents-first a state carries no
  weights digest and no zero-weight record, and the next sweep draws every augmentation variable first under the reading,
  weights and mask in force, so a state stored under one reading and installed into a sampler of the other continues under
  the destination's (design Calls made, 5, as amended). What the reading derives (the gaussian count, the nu constant) is a
  function of the destination's weights and mask, which a state does not carry, so an install recomputes none of it.

Resolution and what each value takes, per family (the refusal numbers are [Refusals](#refusals)'):

| family | NULL resolves to | "precision" | "frequency" | predict's weights and weights.test under the resolved reading |
|---|---|---|---|---|
| gaussian, Student-t, heteroscedastic | precision | taken, real w | taken, real w | precision: as today; frequency: 0 and 1 only (9) |
| logistic, hazard.logistic | frequency (counts) | refused (3) | taken, whole w (5) | trial counts (dec-B470), as today |
| nbinom | frequency (counts) | refused (3) | taken, whole w (6) | positive whole w, w new observations summed (dec-B566); new at predict |
| multinomial | frequency (counts) | refused (3) | taken, whole w (7) | 1 a new observation; the rest is [Open calls](#open-calls), 1; new at predict |
| probit, hazard, hazard.probit | frequency (membership) | refused (3) | 0 and 1 only (2) | as today: trial counts at predict (dec-B470); weights.test 0 and 1 |
| ordinal | frequency (membership) | refused (3) | 0 and 1 only (2) | 0 and 1 (dec-B580); above 1 [Open calls](#open-calls), 1; new at predict |
| hurdle.lognormal | frequency (membership) | refused (4) | 0 and 1 only (hurdle-completion's refusal) | 0 and 1, hurdle-completion's (dec-B588) |
| aft | precision (weights refused) | weights refused (4) | refused (4) | as today |

`$setModel` resolves the model it is given against the sampler's family before comparing, so an NA model on a gaussian
sampler created at NULL passes and one stating "frequency" is refused (10).

### Engine: the gaussian-type families and their read-offs (S1)

- `GaussianResponse` holds a per-row count source: null under precision, where the count is 1 per positive-weight row in, as
  now, and the user weights under frequency. The positive-row count becomes `rowCount_`, a double: under precision exactly
  `static_cast<double>(numPositiveWeights_)` as today, under frequency the sum of w over rows in the likelihood (weight
  positive, mask active), recomputed in `installWeights` and wherever the mask recomposes.
  `ChiSquaredScalePrior::drawSigmaSqFromPosterior` takes the count as a double. The SSR, the leaf statistics and the moves read
  the precisions w unchanged.
- Student-t: under frequency a row's composite precision is Lambda, one `Gamma(w (nu + 1) / 2, rate (nu + r^2 / sigma^2) / 2)`
  with the residual unweighted, in place of `w lambda` with `lambda ~ Gamma((nu + 1) / 2, rate (nu + w r^2 / sigma^2) / 2)`; the
  contained gaussian gets the user weights as its counts. The nu step reads the note's generalized statistics, `sum m_i`, `sum
  m_i log Lambda_i`, `sum Lambda_i`, plus the per-grid constant `sum_i [m_i lgamma((nu + 1) / 2) - lgamma(m_i (nu + 1) / 2)]`
  over rows with m_i > 0 only (at m_i = 0 the term is -lgamma(0), infinite), computed when the weights or the mask install
  (creation, `setWeights`, `setData`, a mask change); a state install changes neither and recomputes nothing. `setData` keeps
  the drawn nu (dec-B589) and rebuilds the constant from the new weights. Under precision `drawIndex` runs its shipped
  expression, the constant never built. A row at weight 0 keeps the lambda it draws, as now, out of the composite; under
  latents-first every lambda, such a row's included, is drawn at the sweep's start, and `setWeights` draws none.
- Heteroscedastic: `ConstantVarianceLeaf::accumulate` adds w for a row under frequency and 1 under precision, in its leaf draw
  and its integrated likelihood, which the variance tree moves read. The leaf learns the reading as a member set at
  calibration; any other site that counts a variance leaf's rows (the fused scan path, the combiner's variance glue) reads the
  same member, which the implementer lists in the slice's report.
- The engine's read-offs under frequency, at creation and at their setter doors: a linear or gp leaf covariate's center and
  scale (`standardizationMomentsForColumn`, weighted moments over the rows in) and its `updateStandardization` door; a forest's
  default sd (`basisColumnScale`, weighted) and its `updateBasisScale` door. The response range stays over the rows in, which
  copies do not move. The quantile cut rule is not weighted: its grid is built over unique values, which copies leave as they
  are (plan critique p5: identical grids, original against expanded, with and without useQuantiles).
- The log-likelihood channel: `GaussianResponse::computeLogLikelihood` and `TResponse::computeLogLikelihood` report
  `w log f(y_i | one observation at sigma)` under frequency, the unit-weight density times w; under precision as now.
- A frequency weight of 0 is a row out of the likelihood, as now; a positive real w is taken (dec-B568).

### Engine: the Polya-Gamma families (S2)

- logistic, hazard.logistic: their weights already are counts; the slice only checks the reading flag (precision refused
  above). The engine change is zero, the zero weight being S4's.
- nbinom: under counts a row draws `PG(w (y + r), psi)` and its working response uses `kappa = w (y - r) / 2`; the r grid step
  multiplies each row's log mass by w (`NBShapePrior::computeKernel`'s histogram adds w); the center (`computeShift`,
  `deriveShift`) is the weighted one. The shape `w (y + r)` is a whole number (r's grid and a fixed r hold whole numbers),
  drawn as that many unit draws. The cold start keeps today's `omega` for a row out of the likelihood, never `w` times it.
- multinomial: the combiner's per-category loop draws `PG(w n_i, .)` with kappa `w (y_ik - n_i / 2)`, the counts times w; its
  log-likelihood stays R's. A zero weight composes into the combiner's effective mask beside the zero-trial rows, with no
  zero-trial warning.
- The cap: today's caps stand at w = 1 (y and a fixed r up to one million on nbinom, row totals on multinomial); for w > 1 the
  product, `w (y + r)` with r at its grid maximum for a drawn r and `w n_i`, above dec-B202's one million is refused at creation
  and at every weight or count change. (A product cap at w = 1 with r at the grid maximum of 50 would refuse y from 999,950,
  accepted today.)
- A weight change on either family (`setWeights`, `setData`) stores the counts and rebuilds what is deterministic (nbinom's
  count kernel, center and composite, as latents-first's setActiveRows row has them); it draws nothing, the next sweep drawing
  omega at the new shapes first. `setData` keeps a drawn r (dec-B589).

### Engine: the binomial row (S3)

- `LogisticResponse` gains an optional `trials_` pointer, null for a binary response, where it reads as 1. With trials, `y_`
  holds the successes k (for a binary row k = y, n = 1, the same values). The shape is `lround(w n)` unit draws (w the
  frequency weight, 1 absent), the working response `w (k - n / 2) / omega - offset`, written so that at null trials the
  expression is today's `w (y - 0.5) / omega - offset`, operand for operand. The cold start is `omega = 0.25 max(n w, 1)` and
  the working response `4 (k / max(n, 1) - 1/2) - offset`, so a row of zero trials or zero weight, out of the likelihood under
  the mask, keeps a positive omega and a finite working response.
- The log-likelihood: `w [k log p + (n - k) log(1 - p)] + w c_i` with `c_i = log choose(n_i, k_i)`, per-row constants computed
  once per install of the counts with R's `lchoose`, so a binomial row reports `w log dbinom(k; n, p)` and a binary row today's
  value. `std::lgamma` is not used: at n = 1e6, k = n / 1000 its constant is 2.3e-10 off in relative terms against `dbinom`,
  R's route 2.3e-12 (plan critique p4, rerun here). The engine reaches R's math only through nmath.h's names (host-neutral
  S1, landed), which hold six routines and no `lchoose`; S3 adds a seventh, `ext_lchoose`, `Rf_lchoose` on the R build and
  `dbarts_nmath_lchoose` under `DBARTS_HOST_STANDALONE`, copying R's choose.c into standalone/nmath beside the lbeta.c,
  lgamma.c and lgammacor.c it calls, which the copy already holds, with its line in the standalone-math driver's argument
  check and the tests/cpp nmath suite (bitwise against Rmath on the reference host).
- The bridge hands successes `counts[, 1]` and trials `rowSums(counts)` from the data object's counts slot, the slot multinomial
  count rows use (design Calls made, 9); `validateResponseSupport` takes logistic successes above 1 where trials are given, up
  to the trials; `$setCounts` on a binomial sampler swaps both and draws nothing, as latents-first leaves a count swap: it
  stores the successes and trials, recomputes the log choose constants and rebuilds the working response from omega, the
  counts and the offset at every row, inactive rows included, as latents-first's logistic `setWeights` does; the next sweep
  draws omega at shape `w n` first. A row of 0 trials installs into the mask as dec-B133's zero-trial row does, with its
  warning.
- Draw for draw: a binomial row of k successes in n trials consumes the n unit draws at the same psi that two adjacent rows of
  weights k and n - k consume, and its kappa is their kappa's sum in exact arithmetic (`k - n/2` against `k/2 - (n - k)/2`,
  each exact for whole numbers below 2^52). The omega sums differ only in the association of n additions, so the gate is a
  tolerance of 1e-10 on the fits, not bitwise.
- The spec's form (for host-neutral-engine's S2b): `y` the successes and a `trials` field, per row, NULL for a binary
  response, under LOGISTIC only; the build refuses `trials` on any other family.

### Membership and zero weights (S4)

- `ProbitResponse` and `OrdinalResponse` take the caller's 0/1 weights at creation, `setWeights` and `setData` and compose them
  with the active-row mask into the effective mask, as `GaussianResponse::recomposeActiveRows` composes weights and mask: a row
  is in when its weight is 1 and the mask leaves it in. A row the composition reactivates draws nothing at the setter: under
  latents-first no mask change draws, and the next sweep draws its latent first. A weight other than 0 or 1 there keeps
  today's refusal; the bridge's creation refusals of probit and ordinal
  weights (`enforceBinaryWeightPolicy`, the ordinal refusal beside it) are lifted to admit 0/1 vectors.
- The R policy stops clearing the weights slot on probit and ordinal and stops returning a mask: `s$data@weights` holds the
  0/1 vector, `s$setData(s$data)` keeps the rows out, `$activeRows` is only what `setActiveRows` installed (NULL after
  creation from 0/1 weights, where today it shows the mask: plan critique p2, 5 zeros). `dbartsSpec()` drops its `active`
  element, and `new("dbartsSampler", ..., active = )` stays the creation door of a mask that composes with the weights.
- `LogisticResponse` takes a 0 as the mask: the row is skipped in the draws (consuming no variate, the property its mask skip
  protects), keeps its cold-start omega, and is NaN in the log-likelihood; `enforceBinaryWeightPolicy` and
  `refuseNonCountWeights` accept 0. dbartsData's warning that zero-weight rows among other weights stay in the fit stays, now
  true on logistic as on every family; what goes is the refusal that followed it there.
- nbinom's 0/1 weights are counts from S2, a 0 out of the likelihood under the same composition.
- xbart's refusal of 0/1 weights on probit goes (design Calls made, 14): each fold's data object carries its rows' weights, the
  fold's creation installs them, and the held-out loss weights a 0 row out.

### The R surface (S5, S6)

- weights.type (S5): validated, resolved and refused in [Refusals](#refusals)' words, per the resolution table; fixed at
  creation, so `$setModel` with another resolved reading is refused beside forest-door's `$setModel` refusals, and
  `setWeights`, `setData` and predict read new weights under the one in force.
- Read off the expanded rows under frequency in R (S5, design Calls made, 3): the starting sigma's regression with residual df
  `sum(w) - p`, by both of `estimateStartingSigma`'s routes, the dense `residualStandardError` and the sparse
  `sparseResidualStandardError`, and xbart's per-fold estimate through `xbartSigestFunction`; the stand-in, a weighted sd; the
  chi prior's anchor. A stated sigma is taken as stated. The engine's read-offs are S1's and S2's.
- Predictive draws (S5): `sampleFromPPD` stops reading the fit's weights for a frequency or count fit's training rows: one
  observation at sd sigma (dec-B565) or Bernoulli(p), one nbinom count, one multinomial category or count vector at the row's
  trials (dec-B566). A row of weight 0 is NA in the predictive draws on every family (dec-B570, dec-B580), drawn and then
  overwritten so the random stream is the one it is today; on probit, ordinal, hazard, hazard.probit and nbinom, which draw such
  rows today (plan critique p1: 0/1 and count draws, fitted(type = "ppd") of 0.25 and 1.05 at those rows), this moves a value by
  dec-B580; the row's fitted value is still reported, and predict(weights = 1) draws it as a new observation. predict(weights =
  ) and weights.test on a frequency gaussian-type fit take 1 and 0 only (dec-B565); on a count family they ask for that many
  trials, at training rows and at newdata (dec-B566 revising dec-B470).
- predict's weights on the ordinal, negative binomial and multinomial fits (S5). dec-B580 names predict(weights = 1) as the door
  to a 0 row's draw "on every family", but predict.bartOrdinal, predict.bartNegbin and predict.bartMultinomial have no weights
  formal and refuse the name with "this family's posterior-predictive draw takes no per-observation weight" (run on main
  7937d252f, scratch/fw-amend/probes/predict-weights.R); 0.9-34 had none of the three families. Each gains `weights = NULL`
  after `offset`, predict.bart's place (dec-B07) and the hurdle's, with dec-B470's rule as on every family: with no newdata
  given weights replace the fit's own for the training rows' draws, one value or one per training row, worked from the
  stored draws; the types that draw nothing read no weight; a missing value follows na.action (dec-B482). What each takes:
  - ordinal, whose weights are membership (dec-B569): 1 draws one new observation, a category, and 0 gives NA, dec-B580's
    and dec-B565's doors, as the hurdle's predict takes them (dec-B588, built on hurdle-completion); above 1 is
    [Open calls](#open-calls), 1.
  - nbinom, a count family (dec-B566): a positive whole w draws the sum of w new counts at the row's mean mu and shape r,
    one draw from NB(size w r, mean w mu) per posterior draw, as a logistic fit's w draws Binomial(w, p); 0 and a fraction
    refused in `refuseNonCountWeights`' words, as logistic's predict refuses them (dec-B470: "positive integer trial counts";
    run on main: "the posterior predictive 'weights' of a probit fit are trial counts and must be positive integers").
  - multinomial: 1 draws one new observation in the shape that call returns without weights (a category, or at a count
    fit's training rows a count vector at the row's trials); the rest is [Open calls](#open-calls), 1.
  - `foreignArgsFor` already drops a method's own formals from `predictForeignReasons`, so the three stop refusing the name
    once they take it; the list's `weights` entry, with no method left to refuse it once hurdle-completion has landed,
    goes.
- The hurdle (S5): hurdle-completion takes 0/1 weights on both parts (dec-B588). Its weights.type row is probit's: NULL
  resolves to frequency, "frequency" takes 0 and 1 with hurdle-completion's refusal of other values unchanged, "precision" is
  refused in the hurdle's words (Refusals, 4) before the zero part could refuse it in probit's; the hurdle arm hands the
  resolved value to both parts, where 0/1 weights draw the same under either reading (the positive part's count is the
  number of 1s either way), so every hurdle fit hurdle-completion takes draws as before, pinned by a before/after arm in
  test-hurdle.R. About 20 lines and that test.
- The log-likelihood (S5): `pointwiseLogLikelihood` reports `w log f(y_i | one observation)` under frequency (dec-A116,
  dec-B565), NaN at a 0 on every family (dec-B570); it and the engine's channel (S1) agree, which a test pins.
- The binomial (S6): the formula parse takes `cbind(k, f)` under family = binomial, stats::binomial() and "logistic", and the
  matrix interface takes a two-column `y` under those families; dbartsData takes the pair as `counts =`, and a data object
  with a two-column counts slot fits a binomial under family = binomial or "logistic" while "auto" on it fits multinomial, as
  it does today (plan critique p7); the "auto", probit-binomial, logistic counts-slot and xbart two-column refusals are
  reworded; the proportion spelling is refused (dec-B573). fit$y keeps the two-column response; fitted(type = "ev") is p and
  its residual `k / n - p`; fitted(type = "ppd") is `n p` and its residual `k - n p`; extract and predict at training rows
  draw Binomial(n, p) (dec-B574); predict at newdata draws Binomial(w, p) with w from its weights, one trial without (design
  Calls made, 10); the zero-trial row is dec-B133's (design Calls made, 13).
- xbart (S6, on S0's rule): the built-in losses on a binomial row read it as its trials (design Calls made, 12): S0's
  probability-scale rule with each row weighted by its trials times its weight and y the proportion k / n; a binomial xbart
  receives `cbind` as bart does; a user loss receives `y.test` as `k / n`, the latent draws, and `weights` as trials times
  weight.
- Help (S4 to S6): weights and weights.type on bart, dbarts, dbartsSpec and xbart, with the resolution table; the help's
  cautions (design Calls made, 8; dec-B568's rescaled and survey weights); the two scales of a binomial fit; the sampler's
  `setData` and `setWeights` mask sentences; dbarts-embedding.Rd's `$setWeights` row, which names survey weights as a
  precision; samplePriorPredictive.Rd's Binomial(w, p) sentence, which becomes one observation.
- NEWS (S5, S6), against 0.9-34 only: the features (weights.type and frequency weights, the binomial response, counts on nbinom
  and multinomial, zero weights on every family, probit and ordinal membership kept on the data object), and NA predictive draws
  at a 0 row on every family (dec-B580), where 0.9-34 gave a probit fit with any 0 weight NA at every row. Nothing about
  logistic's training draw or its zero weight, nor about predict's weights on ordinal, nbinom and multinomial fits: none of
  those families is in 0.9-34. S0's NEWS line is its own, on main.

## S0: the binary xbart rmse, alone on main (LANDED)

Landed on main at 1a23f9a60, before the rebase of 2026-11-05; the section stands as the record of what it built. A defect in
0.9-34 as on main (TODO xbart-binary-rmse, now closed): on a binary response `xbartLossFunction`'s rmse averages the latent
draws, before the link, against the 0/1 response. Plan critique p3, rerun here: 0.9-34's binary rmse 1.70 against a response
sd of 0.50; main's, on the first critique's probe, 0.965 against 0.501.

- Tier "R surface"; RNG class neutral (no draw moves; the loss value does). About 60 lines: R ~10, help ~6, NEWS ~4, tests ~40.
  Branch xbart-binary-rmse off main, one sonnet implementer, one sonnet review; lands on main when green.
- Change: in `xbartLossFunction`'s rmse, when `control@binary`, take `p.test <- rowMeans(probFromLatent(testSamples))` (the
  link per draw, averaged last, as the log and mcr losses already do) and return the root of the mean, or the weighted mean by
  `weights`, of `(y.test - p.test)^2`; a continuous response is unchanged. A user-supplied loss still receives the latent draws.
- Help (man/xbart.Rd, loss): rmse on a binary response compares the response with the posterior mean probability; a supplied
  function's `y.test.hat` is on the latent scale for a binary response.
- NEWS: on a binary response xbart's rmse loss compares the response with the predicted probability, where it compared it
  with the latent scale, as 0.9-34 did.
- Tests (test-xbart-loss.R): a binary arm, probit and logistic, weighted and unweighted, where the rmse equals the hand
  computation from `pnorm` or `plogis` of the same draws (one fit, its test draws read through a supplied loss that returns
  them, then the built-in on the same seed); test-reproducibility-xbart.R is gaussian and stays.
- Mutant: rmse on `rowMeans(testSamples)` - the new arm fails.
- Gates: the touched test file, lint_package, air, codoc, doc-freshness; CI on the branch.
- S6 builds on it: the binomial reading adds the trials to the weights and k / n to y, no second rewrite.

## Refusals

R-level words, refused in the R layer with the bridge keeping a C backstop worded the same where one exists today. `<family>`
is the family token as the caller wrote it.

1. `'weights.type' must be NULL, "precision" or "frequency"`.
2. On probit, ordinal, hazard, hazard.probit with weights not all 0 and 1 under "frequency": `weights.type = "frequency" is not
   taken under family = "<family>" for weights other than 0 and 1: a frequency weight there needs a hidden value for each copy,
   which the sampler does not draw; weights of 0 and 1 mark rows out of and in the likelihood under either reading`.
3. "precision" on logistic, hazard.logistic, nbinom, multinomial, probit, ordinal, hazard, hazard.probit, with or without
   weights: `weights.type = "precision" needs a residual scale, and family = "<family>" has none: its weights count identical
   observations or, at 0 and 1, mark rows in and out; leave weights.type NULL`.
4. aft keeps `aft (survival) models do not support case weights`; "frequency" there, with or without weights:
   `weights.type = "frequency" is not taken under family = "aft", which takes no weights`. hurdle.lognormal, which takes 0/1
   weights (dec-B588), keeps hurdle-completion's refusal of other values under NULL and "frequency"; "precision" there, with
   or without weights: `weights.type = "precision" is not taken under family = "hurdle.lognormal": its weights of 0 and 1
   mark rows in and out of both parts, and its zero part, a probit, has no residual scale; leave weights.type NULL`.
5. logistic and hazard.logistic, a weight negative, fractional or NA: `logistic weights are observation counts and must be
   non-negative whole numbers`.
6. nbinom: `nbinom weights are observation counts and must be non-negative whole numbers; an exposure belongs in the offset as
   a log-exposure term`.
7. multinomial: `multinomial weights are observation counts and must be non-negative whole numbers`.
8. The cap for w > 1, through `refuseCountOverCap`'s form: `<caller>: family "logistic" requires 'weights' times the trials no
   larger than 1000000; each unit of a count is a Polya-Gamma draw per sweep`, and the same with `"nbinom"` and `'weights'
   times (y + r)`, `"multinomial"` and `'weights' times the 'counts' row totals`. At w = 1 today's cap messages stand.
9. predict or weights.test on a frequency gaussian-type fit with weights other than 0 and 1: `on a fit with weights.type =
   "frequency" <what> takes weights of 1, a new observation, and 0, a row held out: a row of weight w stands for w observations,
   and its draw is one of them`, `<what>` being `predict` or `'weights.test'`.
10. `$setModel` with another resolved reading: `$setModel cannot change weights.type, which is fixed when the sampler is
    created; create a new sampler`.
11. The proportion spelling (dec-B573), under family = binomial, stats::binomial() and "logistic", a response strictly between
    0 and 1 that is not 0/1: `under family = binomial the response is 0/1 or cbind(successes, failures): a proportion with the
    trials as weights is written cbind(k, n - k), with k = round(p * n) from a proportion p and its trials n`.
12. The probit binomial: `family = binomial("probit") does not take a cbind(successes, failures) response: the probit binomial
    is not fitted; use family = binomial, the logit link, or a 0/1 response`.
13. A binomial count negative, fractional or missing outside na.action: `a binomial response's successes and failures must be
    non-negative whole numbers; row(s) <rows>`; failures written so that k exceeds n reach it as a negative failure count.
14. "auto" with a two-column response: today's message with the binomial named first: `'y' is an n x 2 matrix, which family =
    "auto" leaves ambiguous - successes and failures need family = binomial; per-category counts need family = "multinomial";
    a (time, status) pair needs a survival::Surv response (which "auto" reads as survival) or family = "aft" / "hazard"`.
15. `$setResponse` on a binomial sampler, through `refuseCountsMutation`: `$setResponse: its response is the n x 2 count pair of
    successes and failures, which a flat vector cannot express; replace it with $setCounts`.
16. dbartsData's two-column y keeps its classed refusal, reworded: `'y' is an n x 2 matrix; dbartsData() takes a single-column
    response - successes and failures and per-category counts go to dbartsData(counts = ), a (time, status) pair to
    dbarts()/bart() with family = "aft"/"hazard"`.
17. Logistic on a data object whose counts slot has K >= 3 columns, reworded from today's `the data carry an n x K count
    matrix, which only family "multinomial" fits; ...`: `the data carry an n x <K> count matrix: family "logistic" takes two
    columns, successes and failures, and family "multinomial" fits <K> categories`.
18. xbart with a matrix response other than a binomial pair, reworded from today's `xbart() takes a single-column response but
    'y' is an n x K matrix, ...`: `xbart() takes a single-column response or, under family = binomial, cbind(successes,
    failures); 'y' is an n x <K> matrix, and it cross-validates neither multinomial counts nor a survival response; use bart()
    instead - per-category counts need family = "multinomial"; a (time, status) pair needs a survival::Surv response (which
    "auto" reads as survival) or family = "aft" / "hazard"`.
19. predict on an nbinom fit with weights 0, fractional or negative, through `refuseNonCountWeights`' form: `the posterior
    predictive 'weights' of a nbinom fit count new observations, whose draws are summed, and must be positive integers`.
20. Under the open call's recommendation, predict on an ordinal fit, or on a multinomial fit whose rows are categories, with
    weights other than 0 and 1: `on a fit whose response is one category per row predict takes weights of 1, a new
    observation, and 0, a row held out: a draw over several trials would be a count per category, which only a fit to count
    rows predicts`; on a multinomial count fit, 0, fractional or negative weights: `the posterior predictive 'weights' of a
    multinomial count fit are trial counts and must be positive integers`.

Gone: xbart's `xbart does not accept weights of 0 and 1 under family ...`; logistic's refusal of a 0; the logistic refusal of
a two-column response under family = binomial and "logistic"; the multinomial refusal of weights other than all ones; the
bridge's creation refusal of 0/1 weights on probit and ordinal; predict's refusal of `weights` by name on ordinal, nbinom and
multinomial fits ("this family's posterior-predictive draw takes no per-observation weight").

## RNG class

Neutral on the engine for every call taken on the build this branch starts from, main after latents-first: latents-first moves
every augmenting family's draws once and re-records the equivalence trio and test-reproducibility-binaryResponse.R at the end
of its slice 2b, so "today" below means that build and those baselines, never main before it. Shown three ways:

- Bitwise harnesses, run against each engine slice's own `--preclean` library: `benchmarks/R/equivalence.R compare`,
  `bcf-equivalence.R compare` and `multinomial-equivalence.R compare`, every scenario "identical draws (same RNG stream)",
  counted per [Gate hygiene](README.md#gate-hygiene); the four `test-reproducibility-*.R` on the reference build, unchanged
  (two of them carry precision weights, which pins the precision path); tests/cpp unchanged in every existing test. If
  approximate-draws' S2 has landed first with a default that routes large shapes, the harness runs as built; every
  draw-for-draw gate here against adjacent copies runs at the exact sum's setting, since under the accept-or-keep route one
  row of shape `w n` and its w copies are different chains.
- Construction: under precision every changed expression is today's, operand for operand (the count is
  `static_cast<double>(numPositiveWeights_)`, the variance leaf adds 1.0, the nu constant is not built, the moments and the
  basis scale read unit weights, a null `trials_` reads as today's working response); each slice's reviewer reads the diff for
  it and runs the mutants that would break it.
- The membership paths (S4): a probit and an ordinal fit with 0/1 weights equal today's draw for draw, the mask installed by
  the engine at creation against today's install by R after creation, which draws no variate (under latents-first no setter
  draws, and creation's starts are deterministic); an nbinom fit with 0/1 weights equals the documented by-hand route,
  `new("dbartsSampler", ...)` then `$setActiveRows(w, updateScale = TRUE)`, on the same build, a live oracle; pinned by
  test-binary-weight-mask.R and test-nbinom.R unchanged in their draws.

Reported values that move, each re-pinned in the test that holds it:

- S0, landed on main: a binary xbart's rmse loss moved to the probability scale; its default, log, did not. A defect fix
  against 0.9-34, in NEWS.
- A weighted logistic or hazard.logistic fit's training predictive draws (extract, fitted(type = "ppd"), predict at training
  rows without weights): Bernoulli(p) where they drew Binomial(w, p) (dec-B566). R's random stream after such a draw differs
  too, since rbinom consumes by size; test-weighted-binary-ppd.R and test-weighted-logistic.R re-pin, replayed as whole files.
  No snapshot file holds a weighted logistic draw. Not in NEWS: logistic is not in 0.9-34.
- A 0 row's predictive draws on probit, ordinal, hazard, hazard.probit and nbinom: NA where they draw values today (plan
  critique p1), by dec-B580. Drawn then overwritten, so R's stream does not move; the fitted value stays reported. In NEWS
  against 0.9-34, whose probit fit with any 0 weight gave NA at every row.
- `$activeRows` and fit$active after creation from 0/1 weights on probit, ordinal and nbinom: NULL where they show the mask
  (plan critique p2); the weights carry it.
- A logistic zero weight, refused today, is taken: NaN in the log-likelihood, NA in the draws.

New calls take new draws, each held to an oracle: the expanded rows draw for draw (gaussian, heteroscedastic, nbinom at a
fixed and at a drawn r, multinomial, binomial), the expanded rows in distribution (Student-t only), the exact-posterior arms
of S7; predict's weights on nbinom and ordinal (and multinomial as the open call rules) against the same R draw written by
hand from the stored draws under one seed. A hurdle fit with 0/1 weights does not move (S5's before/after arm). The
exact-gates list grows by one script (S7); no existing gate's arm changes.

## Calls made

1. The engine form of the binomial is successes in `y` and a `trials` vector, not the proportion k / n with the trials as a
   count: kappa `w (k - n/2)` is then exact, which the draw-for-draw gate against the expansion needs, and the spec a C caller
   fills holds whole numbers. Alternative: y = k / n, where `n (k / n)` need not round back to k.
2. The log choose constants are per-row values computed at each install of the counts with R's `lchoose` route, not
   `std::lgamma`, and not evaluated in the log-likelihood loop. Alternative: three `std::lgamma` per row, 2.3e-10 off at the
   cap (plan critique p4).
3. The frequency count reaches the gaussian family as a second pointer (counts beside precisions), null under precision; a
   Student-t sampler hands its contained gaussian the user weights as counts and Lambda as precisions. Alternative: a count
   derived from the precisions, which Student-t's composite breaks.
4. Probit and ordinal compose the 0/1 weights with the mask in the family, as gaussian composes weights and mask, so a C caller
   gets membership from the weights with no R policy. Alternative: the bridge composes and installs the product through
   `setActiveRows`, which leaves the engine unaware of the weights and host-neutral's build without the rule.
5. `dbartsSpec()` drops `active`. A caller reading `spec$active` gets NULL and `setActiveRows(NULL)` installs nothing, so
   stan4bart's code stays correct unchanged and its plumbing is dead; dbartsSpec is not in 0.9-34. Alternative: keep `active` as
   a copy of the 0/1 weights, which a caller following today's help would install twice, harmless but misleading.
6. "precision" on a family with no residual scale is refused even without weights, since the reading is fixed at creation and
   a later `setWeights` would read under it. Alternative: refuse only when weights are given, which defers the refusal to a
   setter.
7. weights.type on the model record is NA for "the family's own" and the spec's resolution writes the resolved value into the
   model dbartsSpec returns; the bridge reads only the resolved two values. Alternative: a third engine value, which
   host-neutral's spec declines (its Calls made, 22).
8. The multinomial counts times w are formed in the combiner's loop, not as a bridge-owned product matrix, so `$setCounts`,
   `$setWeights` and a state install each have one place to read. Alternative: the product in the bridge, which needs its own
   lifecycle across three setters.
9. On dbartsData the binomial pair is `counts =` (the counts slot, design Calls made, 9) and the two-column `y` stays refused,
   reworded; on bart, dbarts and xbart the pair is the response under the logistic family; "auto" on a counts data object stays
   multinomial. Alternative: dbartsData taking a two-column y as counts, which it cannot tell from a (time, status) pair
   without the family.
10. The binary rmse fix is S0, alone on main now: it fixes a defect 0.9-34 has, needs nothing of the feature, and S6 then adds
    the trials. Alternative: inside S6, which ties a defect fix to a feature of 140 to 190 agent-hours.
11. The engine read-offs under frequency (leaf covariate moments, a forest's default sd, nbinom's center) and the engine's
    log-likelihood channel are built in S1 and S2 with the families they serve, so the R slices hold no engine code.
    Alternative: one read-off slice, which would review engine and R code at the R tier.
12. The cap applies to the product only for w > 1, today's caps standing at w = 1. Alternative: the product at every w, which
    refuses nbinom counts from 999,950 that are accepted today.
13. A state carries no weights reading (design Calls made, 5, as amended): under latents-first nothing in a state is read by
    the next sweep's draws but the trees, sigma, the family parameters and the generator, so the destination's reading
    governs. Alternative: storing the reading and refusing an install under the other, a refusal of an install that is
    valid.
14. nbinom's predict weights read as w new observations summed, one NB(size w r, mean w mu) draw, following dec-B566's "given
    weights ask for that many trials" on the count families and logistic's Binomial(w, p); 0 refused as on logistic.
    Alternative: w independent counts returned per row, a new shape no other count family returns.
15. The hurdle hands its resolved reading to both parts, where 0/1 weights draw as under either reading. Alternative: the
    zero part's reading alone, the positive part left at NULL, the same draws and two readings on one fit's record.
16. `lchoose` reaches the engine as a seventh nmath.h name with choose.c copied into standalone/nmath, as host-neutral S1 did
    for the other six. Alternative: a hand-written lchoose over the copied lbeta, a second implementation to keep equal to R's.
17. Each new Polya-Gamma shape goes through the family's existing unit-draw loop ([Constraints](#constraints)), so
    approximate-draws' route reaches it with no edit here. Alternative: a loop per new shape, which that plan would replace
    three more times.

## Open calls

One, put to the maintainer as written. Every other user- or consumer-visible choice the plan makes is settled by dec-B563 to
dec-B570, dec-B573, dec-B574, dec-B580, dec-B588 or the design note's calls made.

### 1. What predict's weights mean on a fit whose one draw is a category

On an ordinal or multinomial fit, what should predict(fit, type = "ppd", weights = 3) draw?

Background. predict's weights say what kind of new observation a row's draw is. On a logistic or probit fit a weight of 3
asks for three new trials and the draw is the number of successes, 0 to 3, a count on the scale of the 0/1 response (run
today on a shipped build of main: every draw between 0 and 3). This plan reads 3 on a negative binomial fit as three new
counts summed. An ordinal or multinomial draw is a category, so three trials make a count for each category, a different
shape from what predict returns now. Run the same day: on a multinomial fit whose rows are categories, extract and predict
return one category per row and draw (50 draws by 60 rows); on a fit whose rows are count vectors, extract at the training
rows returns a count per category at each row's trials (50 by 60 by 3) while predict at new rows returns one category per row
(50 by 5), as you ruled on 2026-10-01, a new row's number of trials being unknown. Today all three predicts refuse weights by
name. Settled already: the argument itself, weights = 1 drawing one new observation, and on ordinal weights = 0 giving a
missing draw (your zero-weight rulings). 0.9-34 had neither family. Base R gives no lead: predict.lm's weights scale its
prediction intervals, and glm's, MASS's polr and nnet's multinom predict take none.

- (a) Counts only where the fit's rows are counts. A multinomial fit to count rows takes any positive whole number as the
  row's trials and returns a count per category, at the training rows as now and at new rows, which gives a count fit the
  trials at new rows it lacked; 0 is refused there, as on logistic. An ordinal fit and a multinomial fit to categories take
  0 and 1 and refuse larger values by name. About 60 lines of R and 80 of tests. Cost: a multinomial fit reads predict's
  weights two ways by its response, and a count fit's predict at new rows returns categories without weights and counts with
  them.
- (b) 0 and 1 only, on both families, every draw a category. About 20 lines. Cost: a count fit still cannot draw counts at
  new rows, the one thing such weights would be for, and logistic and the negative binomial take 3 where these refuse it.
- (c) Counts on both families whenever predict is given weights, a single trial written as a count of 1 in its category.
  About 90 lines. Cost: predict(weights = 1), the call you named for drawing a held-out row, returns another shape than the
  same call without weights, and an ordinal fit's draws become counts.

Recommendation: (a). A draw keeps the shape of the response the fit models, your reason on 2026-10-01 for count vectors at a
count fit's training rows; the weights supply the trials that ruling lacked at new rows, as a binomial fit's weights do here;
a fit to one category per row keeps one category per draw, as the hurdle's predict keeps 0 and 1.

## Slices

Each slice is a commit series on frequency-weights, its own review and gates; the feature merges once after S8. S0 landed
apart (above). The branch is rebased onto main once latents-first has merged, before S1 starts, and every gate below runs
against that build.

| slice | tier | content | agent-hours (was) |
|---|---|---|---|
| S1 | changes draws | gaussian-type families, read-offs, log-likelihood channel; no digest | ~12.5 (13) |
| S2 | changes draws | nbinom and multinomial shapes times w, caps, setters storing only | ~7 (7) |
| S3 | changes draws | the binomial row, `ext_lchoose` with choose.c in the standalone copy, `$setCounts` storing only | ~9.5 (9) |
| S4 | changes draws | membership on the data object, logistic's 0; no reactivation redraw | ~6.5 (7) |
| S5 | R surface | weights.type, read-offs in R, predictive draws, predict's weights on ordinal, nbinom, multinomial, the hurdle row | ~14 (10) |
| S6 | R surface | the binomial surface and xbart | ~9 (9) |
| S7 | gate-only | frequency-exact.R, the SBC arm | ~6 (6) |
| S8 | records | records | ~1 (1) |
| reviews, gates, rebases | | four engine reviews ~18, three R and gate reviews ~7.5, gates ~5, rebases ~3 | ~33.5 (32) |

About 99 agent-hours (was 94; S0's 3 spent on main). The moves: S1 and S4 lose the digest and the reactivation redraw that
latents-first retires (about half an hour each); S3 gains the seventh math name (about an hour, less half an hour for a
`$setCounts` that only stores); S5 gains predict's weights on three methods and the hurdle's row (about 4 hours, with
review half an hour more); rebases gain an hour for latents-first's re-recorded baselines.

### S1: engine, the gaussian-type families and their read-offs (~285 code, ~350 tests/cpp, ~80 bridge and R plumbing)

Tier "Changes draws"; RNG class neutral on every call taken today, new draws under frequency.
Steps: `WeightsReading` and its bridge read from the model slot (the slot and its R plumbing land here, refusals in S5); the
gaussian count pointer and `rowCount_`; `drawSigmaSqFromPosterior` taking a double count; Student-t's frequency draw, its
counts, the nu constant over m_i > 0; the variance leaf's count and every site that counts its rows; the weighted leaf
covariate moments and default forest sd with their setter doors; the gaussian and Student-t log-likelihood channel. No digest
(latents-first retired it).
Gates: tests/cpp against w adjacent copies, draw for draw (sigma and the fits within 1e-12 over 300 sweeps on gaussian;
heteroscedastic, its variance leaves additive, the same bound; the critique's p6 gave 8e-15 at a fixed sigma); Student-t
against the copies in distribution and nu's grid probabilities on fixed residuals (the note's p8 arm, within Monte Carlo
error); a linear leaf and a two-forest fit against the copies draw for draw, covering the moments and the default sd; a 0
weight under frequency on each of the three families, finite throughout; the w = 1 identity (frequency at all-ones weights
equals precision bitwise, and at 0/1 weights on a gaussian fit with constant leaves); a frequency sampler masked by
`setActiveRows` against one created with those rows' weights at 0, draw for draw, on gaussian and Student-t (the count and
the nu constant recomputed at a mask change); a state, copy and reload round trip of a frequency sampler, and a Student-t
frequency state installed into a precision sampler of the same data consuming no variate at the install; the log-likelihood
channel against
`dnorm` and `dt` times w; the bitwise trio; the full tinytest suite; sanitizers on tests/cpp and on the R-loaded path;
`mutation-battery.R verify-anchors`.

### S2: engine, the Polya-Gamma families (~110 code, ~250 tests/cpp, ~40 bridge)

Tier "Changes draws"; neutral today (nbinom and multinomial weights were refused or masks).
Steps: nbinom's shape, kappa, r-step power, weighted center and cold start; the combiner's multinomial shape and kappa; the
caps at w > 1; the bridge stops refusing nbinom and multinomial weights; their weight changes store and rebuild, drawing
nothing.
Gates: tests/cpp: an nbinom row of weight w against w adjacent copies draw for draw, at a fixed r and at a drawn r (the r
draw reads the same histogram, the unit draws are consumed in row order), within 1e-10; a multinomial row of weight w against
w adjacent copies (unit draws in row order inside each category's loop) and against the counts times w, draw for draw; a 0
weight on each family, finite throughout; y at 999,999 with w = 1 taken, the product cap at w = 2 refused; a weight change
consuming no variate; a state round trip; the trio; the suite; sanitizers.

### S3: engine and bridge, the binomial row (~125 code, ~25 standalone math and choose.c copied, ~220 tests/cpp, ~50 bridge)

Tier "Changes draws"; neutral for a binary response.
Steps: `trials_`, the working response, the cold start, the log choose constants and the log-likelihood; the bridge reads the
counts slot under logistic and `validateResponseSupport` admits successes up to the trials; `$setCounts` reaches a binomial
sampler and stores, drawing nothing; the zero-trial row into the mask; `ext_lchoose` in nmath.h, choose.c in
standalone/nmath, its line in the standalone-math driver and the nmath suite.
Gates: tests/cpp: a binomial row against two adjacent rows of weights k and n - k for 0 < k < n, fits within 1e-10 over 300
sweeps, the same number of variates consumed; a binomial row at frequency weight w against the row of k w successes in n w
trials, bitwise; the log-likelihood against R's `dbinom(k, n, p, log = TRUE)` times w within 1e-10 relative, at n up to 1e6
with k = n / 2 and k = n / 1000 named; a zero-trial row and a zero-weight binomial row finite throughout; a state, copy and
reload round trip of a binomial sampler; `$setCounts` consuming no variate and leaving the working response equal to the one
rebuilt from omega, the counts and the offset at every row; `ext_lchoose` against Rmath bitwise on the reference host and the
`standalone-math` target building with it; a null-trials sampler bitwise today's; the trio; the suite. Today a weight-0 logistic
row still draws one variate, so k = 0 and k = n join the draw-for-draw test in S4, once the expansion's weight-0 row is masked
and draws none.

### S4: membership and zero weights (~50 engine, ~30 bridge, ~80 R, ~180 tinytest)

Tier "Changes draws" (the mask paths of four families change); neutral on every call taken today.
Steps: probit and ordinal take weights and compose; the bridge's creation refusals lifted for 0/1; logistic takes a 0;
nbinom's 0/1 through S2's counts; the R policy keeps the weights slot; `dbartsSpec()` drops `active`; xbart's refusal goes;
help for `$setWeights`, `$setData`, `$setActiveRows`, dbartsSpec.
Gates: the binomial draw-for-draw test of S3 at k = 0 and k = n; test-binary-weight-mask.R's draws unchanged (its
`s$data@weights`, `$activeRows` and `spec$active` assertions re-pinned to the ruling); nbinom with 0/1 weights bitwise against
`new("dbartsSampler", ...)` then `$setActiveRows(w, updateScale = TRUE)` on the same build; a new test that `s$setData(s$data)`
keeps masked rows out and that a mask set by `setActiveRows` composes with the weights in either order, a composition that
reactivates a row consuming no variate; a logistic fit with a
zero weight equal draw for draw to the fit with that row masked by `setActiveRows`; xbart on probit with 0/1 weights equal,
fold for fold under a fixed fold assignment, to a by-hand loop fitting each fold with those weights and scoring its held-out
rows of weight 1 (the test-xbart-fold-oracle.R pattern); the suite; the trio.

### S5: R surface, the reading (~435 R, ~135 help, ~600 tinytest)

Tier "R surface", pure R. Moves the weighted logistic training draw and, by dec-B580, a 0 row's draws on probit, ordinal,
hazard, hazard.probit and nbinom; re-pinned.
Steps: weights.type on bart, dbarts, dbartsSpec and xbart, its validation and the resolution table, the hurdle's row and its
routing to both parts; Refusals 1 to 10, 19 and 20; the starting sigma on both routes and xbart's per fold, the chi anchor;
`sampleFromPPD` (one observation, NA at a 0 on every family, dec-B580, the fitted value kept), predict's and weights.test's
doors, predict's `weights` on predict.bartOrdinal, predict.bartNegbin and predict.bartMultinomial and their help pages,
`pointwiseLogLikelihood`; help and NEWS.
Gates: the touched test files and the new test-frequency-weights.R (each refusal by its words, counted warnings; the
resolution table row by row, hurdle.lognormal's included; sigma's starting value and the chi anchor against the expanded data
on the dense and the sparse routes and in xbart's folds; predict's and weights.test's doors; NA draws at a 0 on every family
with the stream unmoved; predict(weights = 1) at a 0 training row of an ordinal, an nbinom and a multinomial fit drawing a
value, equal under one seed to the draw written by hand from the stored draws; nbinom predict(weights = 3) equal under one seed
to `rnbinom(size = 3 r, mu = 3 mu)` from the stored draws; the open call's arms as ruled; `pointwiseLogLikelihood` against
the engine's channel); test-hurdle.R's before/after arm (a hurdle fit with 0/1 weights at weights.type NULL and "frequency"
equal draw for draw to the build before S5); lint_package, air, codoc, win drift, doc-freshness; the snapshots on the
reference build.

### S6: R surface, the binomial and xbart (~250 R, ~55 help, ~400 tinytest)

Tier "R surface", pure R. Builds on S0, landed.
Steps: the parse and routing, the counts slot under the logistic family; Refusals 11 to 18; fit$y, fitted, residuals,
extract, predict's trials, the zero-trial row; the log-likelihood channel; xbart's binomial reading and its response; help and
NEWS.
Gates: test-binomial.R (against glm on the same data: fitted(type = "ev") near glm's fit on a large-n design, residuals'
definition exactly, the log-likelihood's constant against `dbinom`, predict's trials, the zero-trial row's 0 log-likelihood and
its -p residual, "auto" on a counts data object still multinomial, every refusal by its words); test-xbart-loss.R with an arm
showing a binomial xbart equals xbart on the trial-level expansion for the log and misclassification losses exactly; the
R-surface checks of S5.

### S7: gates (~450 in benchmarks/R, ~20 workflow)

Tier: gate-only, the R-surface review.
Steps: benchmarks/R/frequency-exact.R with arms at the R level: gaussian frequency against the expanded rows draw for draw
over 300 sweeps (the note's p2 emulation as the seed, sigma and the fits within 1e-12); the same with a linear leaf, with two
forests and with useQuantiles; heteroscedastic, nbinom (fixed and drawn r) and multinomial against adjacent copies, draw for
draw; binomial against its two-row expansion within 1e-10; an exact-posterior arm for sigma on a single-leaf tree under real
weights, against the analytic inverse-chi-squared posterior; Student-t's nu on fixed residuals against its exact grid
probabilities, and Student-t against the copies by Welch z on the fits, quick and long. It joins exact-gates.yaml's list. An
SBC arm in sbc.R for a frequency gaussian fit at real weights (dispatch, not the matrix, until its ladder earns admission).
Gates: each arm run green, and its discriminating mutant (below) run red.

### S8: records (~60)

The plan's Status and Landing note, the ledger entry for the calls made,
TODO (frequency weights, mask-on-own-data-object), the design Status line, the indexes.

## Tests and gates

Per engine slice, against its own library (`R CMD INSTALL --preclean -l <lib> .`, `buildInfo()$mode` checked): tests/cpp under
`OPT="-O2 -g -fsanitize=address,undefined"`; the full tinytest suite in one process at its default ncpu; the three
equivalence compares counted per scenario; the four snapshot files on a reference build; the R-loaded ASan run of
[Gate hygiene](README.md#gate-hygiene); `mutation-battery.R verify-anchors`; lint_package, air, codoc, win drift,
doc-freshness. Per R slice: the touched test files, the snapshots, the static checks; CI on the branch tip runs the rest. At
the feature's end: the full battery once on the merged tree, `R CMD check --as-cran` from a clean tarball, the exact gates
with `quick`, the x86 leak-sanitizer leg, and bench-sampler.R compare on a quiet machine for the hot paths S1 to S3 touch (the
sigma draw, the variance leaf, the Polya-Gamma loops), unweighted and precision-weighted within noise.

## Mutants

A reviewer applies each to a copy, rebuilds with `--preclean`, and runs the named killer; each must fail it.

1. S1: count each positive row as 1 under frequency (sigma's df) - the gaussian draw-for-draw test, frequency-exact.R.
2. S1: the frequency count under precision too (sigma's df the sum of the weights, not the positive rows) -
   test-reproducibility-continuousResponse-singleThreaded.R on the reference build (precision weights 1 and 2) and the
   equivalence harness's weighted scenarios, bitwise.
3. S1: drop the lgamma constant under frequency - tests/cpp nu grid test (a probability off by up to 0.963 in the note's p8).
4. S1: variance leaf adds 1 under frequency - the heteroscedastic draw-for-draw test.
5. S1: the frequency count and the nu constant left at their creation values after `setActiveRows` - the masked-by-setter
   against created-at-weight-0 draw-for-draw test. (The first plan's mutant 5, the digest ignoring the reading, goes with the
   digest under latents-first.)
6. S1: unweighted leaf covariate moments under frequency - the linear-leaf draw-for-draw test.
7. S1: the nu constant over every row, m_i = 0 included - the zero-weight Student-t test (an infinite constant).
8. S2: nbinom kappa without w - the nbinom draw-for-draw test.
9. S2: multinomial shape `n_i` not `w n_i` - the multinomial draw-for-draw test.
10. S2: the product cap at w = 1 - the y = 999,999 test.
11. S3: kappa `w (k - 1/2)`, the binary expression on the successes - the binomial draw-for-draw test.
12. S3: the log choose by `std::lgamma` - the log-likelihood test at n = 1e6, k = n / 1000.
13. S3: the shape `lround(n)` without the frequency weight - the weighted binomial row against k w in n w, bitwise.
14. S3: cold-start omega `n w / 4` - the zero-trial row test (a NaN in the fits).
15. S4: probit composition takes the mask alone - the setData test and the either-order composition test.
16. S4: logistic draws a 0-weight row - the zero-weight logistic against the masked fit, draw for draw.
17. S5: `sampleFromPPD` keeps Binomial(w, p) on a training row - test-weighted-binary-ppd.R re-pinned.
18. S5: no NA overwrite at a 0 on probit - test-frequency-weights.R's NA arm.
19. S5: predict on a frequency fit takes weights 2 - its refusal test.
20. S5: the sparse sigest route unweighted under frequency - its expanded-data test.
21. S6: the binomial residual `k - p` - test-binomial.R.
22. S3: `ext_lchoose` mapped to `std::lgamma` arithmetic under `DBARTS_HOST_STANDALONE` - the nmath suite's bitwise compare
    against Rmath on the reference host.
23. S3: `$setCounts` keeps the old working response - the tests/cpp working-response read after the swap.
24. S5: nbinom predict(weights = w) draws one count at mean mu - the `rnbinom(size = 3 r, mu = 3 mu)` test.
25. S5: the hurdle hands NULL to its zero part and "frequency" to its positive part with no check of "precision" - the
    resolution table's hurdle row ("precision" refused in the hurdle's words).

## Consumers

- stan4bart (main 27ef163): passes its weights to the BART part at weights.type NULL, a precision on its continuous fits as in
  its own Stan target; its probit rule passes 0/1 weights, which now ride `data.bart` into every sampler it builds, and its
  `active.bart` from `bart_spec$active` becomes NULL, so its `setActiveRows` installs are skipped and `restoreBartSampler`
  rebuilds from the data, which carries the weights (`data.bart.light` drops only x); the critique confirmed its draws do not
  move. No call changes; an optional cleanup of about 30 lines drops the dead `active` plumbing, in the consumer sweep, not
  this feature. Its `install_bart_range` reads `sampler$data@weights` on continuous fits only, unchanged.
- bartCause (1e4b882): passes population weights to bart's gaussian response fits at weights.type NULL, unchanged; drops
  weights from its propensity fits; its xbart call (bayesOpt) uses the default loss, which on a binary treatment is log and does
  not move under S0. No call changes.
- treatSens (268b2db): passes no weights to dbarts; no change.
- bairrtt: reads none.
- latents-first re-records stan4bart's, bartCause's and treatSens's moved draws before this lands (its slice 4); nothing here
  moves them again. None of the four calls predict on an ordinal, nbinom or multinomial fit or fits a hurdle (grep of their R
  trees for those family tokens and hurdle, hurdle-completion's Consumers).

## Interaction with the other pending plans, and order

- forest-door landed (a6a16c31a), with setmodel-parameters: its forest records carry shares, `$setModel` gained refusals in
  R/fittingLevel.R beside which Refusal 10 sits, and chain.hpp, combiner.hpp, model.hpp and R_interface_bartcore.cpp moved; this
  plan's symbols were re-read on 5a99172e6 and none it names changed role.
- Build order, the orchestrator's call as latents-first's plan states it: coefficient-prior-live, latents-first, this
  feature, host-neutral-engine S2a onward. Engine features merge serially; this branch rebases onto main after each of the
  first two merges and re-gates what the rebase touches.
- coefficient-prior-live (built, landing on main 2026-11-05; read at 05a12bb79 against main): touches no file's weight code and
  none of this plan's engine symbols. It does not touch model.hpp; its combiner.hpp and chain.hpp hunks are the live coefficient
  prior (`setForestAmplitudePrior`, `setForestCoefficientPrior`), not `basisColumnScale`, `updateBasisScale` or the
  Polya-Gamma loop. Two textual neighbors: R/fittingLevel.R's `setModelOnForests` gains `refuseCoefficientDifference`, beside
  which Refusal 10 sits; R/spec.R's `resolveSamplerSpec` gains `leafPriorGivesNoSize`, near where the spec resolves
  weights.type. A forest released to a drawn coefficient takes creation's sd, which under frequency is the weighted
  `basisColumnScale` S1 builds, so the two agree with no edit.
- latents-first (planned, builds before this; docs/plans/latents-first.md, Order and overlaps): every sweep draws its
  augmentation variables first; no setter draws them; the weights and survival digests, the Student-t zero-weight record and
  setState's reconciliation retire. What that does here, each rewritten above: design call 5 and the engine bullet of
  [The reading and where it lives](#the-reading-and-where-it-lives) (a state carries no reading; Calls made, 13); mutant 5 and
  the test-state-weight-pairing.R arm (replaced by the mask-change count test); `$setCounts` on a binomial sampler stores
  (S3); a reactivated row draws nothing at the setter (S4); the Student-t install list drops the state install (S1); the
  baselines this plan's neutrality is measured against are latents-first's (RNG class). It also shares `ProbitResponse`,
  `OrdinalResponse`, `LogisticResponse` and `NBResponse`'s setters and the setter help, which latents-first rewrites first.
- hurdle-completion (built, R only, landing before this whenever green; its Interaction with frequency-weights, and order):
  the hurdle takes 0/1 weights on both parts (dec-B588). Here: its resolution-table row is probit's, Refusal 4 refuses
  "precision" in its words, the arm hands the resolved reading to both parts (S5, Calls made, 15); textual overlaps in
  R/bart.R, R/generics.R, man/bart.Rd, inst/NEWS.Rd and test-zero-weights.R. If this feature were to land first instead,
  hurdle-completion's S3 would write its weights in weights.type's terms, as its plan says.
- host-neutral-engine: its S1 (build-time math and printing) landed on main (dec-A260). Its copy of R's math holds six
  routines and no `lchoose`, so S3 here adds the seventh name and copies choose.c (Calls made, 16). Its main feature rewrites
  every creation route's SEXP reads into a plain-C spec (S2b) that carries `weightsType` per dec-B564 and, under this plan,
  `trials` beside y under LOGISTIC. Under the build order this feature lands first: host-neutral's S2b translates the weights
  reading, the counts slot under logistic and the membership rule as code that exists, and its S3 adds the header field and
  the getLatents sentence on a frequency row's latent. Cost: host-neutral's S2a rebases over this feature's bridge changes
  (~200 lines in the creation and mutation conduits it rewrites; latents-first removed about 140 bridge lines there first),
  and S2a waits for this feature to merge. Were the order reversed, this feature would rebase its bridge reads onto S2b's
  spec form and add `trials` to the header itself, about the same lines in the other direction.
- approximate-draws (planned, its three open calls with the maintainer; docs/plans/approximate-draws.md at e0c62e566): the
  binomial at large n is its cost case. How this plan's shapes times w behave:
  - Before it, as this plan builds them: a row's shape (`w` on logistic, `w n` on a binomial row, `w (y + r)` on nbinom,
    `w n_i` per multinomial category) is a whole number drawn as that many unit Polya-Gamma draws, about 160 ns each (164 ns
    at psi 0, 159 ns at psi 3, p12cost.R on the shipped build), so 1,000 rows of 1,000 trials cost about 0.16 s a sweep and a
    1,000-sweep chain about 160 s, linear in the total; fractional weights refused on every Polya-Gamma family (Refusals 5 to
    7).
  - After its S2 (the exact accept-or-keep route at the sweep's start for logistic, binomial, hazard.logistic and nbinom):
    each row's shape is unchanged and the same rows are taken; which way a row is drawn is that plan's setting, and under its
    recommended default a row whose shape is above 64 takes the route, exact, at about 3 to 9 us a row whatever the shape, and
    the rest the sum as here. A row a setter moved (a weight, count or mask change here among them) is drawn by the sum at the
    next sweep. Multinomial's `w n_i` keeps the sum under every exact setting. Draws shift only where the route is selected.
  - After its S3: real weights on logistic, binomial and hazard.logistic are taken under the route at creation, lifting the
    fractional clause of Refusal 5 there and drawing the real shape where this plan rounds `w n` (its Real shapes); under an
    exact setting a between-sweep change that would mark a fractional-shape row is refused by that plan's name. nbinom's and
    multinomial's weights stay whole (Refusals 6 and 7).
  - Order: its S3 must land after this feature, since it reads weights.type, the binomial row and this plan's Refusal 5. Its
    S1 and S2 may land on either side: every new shape here goes through the family's existing loop (Constraints), which S2
    replaces once per family, and if S2 lands first every draw-for-draw gate here runs at the exact sum's setting (RNG class).
    Whichever lands second rebases over the other's `LogisticResponse` and `NBResponse` hunks.
- weighted-binary (parked, after 1.0-0): its integer-weight probit is dec-B563's truncated-normal case, not scheduled; its
  real-weight logistic is the real-shape draw above.
- Nothing else touching src/bartcore/model.hpp's response families merges while S1 to S4 are in review.

## Cost

S0, landed: about 60 lines and 3 agent-hours. The feature about 4,300 lines, 1,560 of them package code; about 99 agent-hours
estimated, 150 to 200 by the record of plans here ([Slices](#slices)). No per-sweep cost on the gaussian-type families at either reading (the
count is fixed between weight installs, the variance leaf adds w in the pass that counts rows now, the nu constant is per grid
point, the read-offs run at creation and their doors); the count families' cost is today's per unit draw, multiplied by the
counts the user asks for.

## Not verified

- The full list of sites that count a variance leaf's rows (the fused scan path, the combiner's variance glue); S1's
  implementer lists them and mutant 4 shows the list complete.
- That the multinomial combiner forms kappa in one place for both the train and a state install; S2 reads it.
- The Student-t nu step's generalized statistics were checked by the first critic's p8 within Monte Carlo error and
  re-derived by the plan critic, not derived anew here.
- The cost figures were measured once on a shared, loaded laptop; bench-sampler.R on a quiet machine is the gate for the hot
  paths, not these.
- stan4bart's, bartCause's and treatSens's call sites were read on main 27ef163, 1e4b882 and master 268b2db; their 1.0 compat
  branches may differ.
- latents-first is planned, not built: what it retires (the digests, the record, the setter draws) and the baselines it
  re-records are read from its plan at 83b23342c. Each rewritten clause here assumes that plan as written; S1's implementer
  confirms on the rebased tree that no `weightsDigest`, `reapplyWeights` or `redrawEntering` remains before building.
- hurdle-completion was read at 8cd1626d1 on its branch (its refusal of other weights and its predict door), not merged; the
  hurdle positive part's draws at 0/1 under "frequency" equalling precision's is argued from the count and constant leaves,
  and S5's before/after arm shows it.
- approximate-draws' route, its default and its real-weight surface wait on its open calls; the before and after above follow
  its plan's recommendations.
- The open call's line counts are estimates; no option was prototyped.
- That choose.c needs nothing beyond the copy's lbeta.c and lgamma.c (`lgammafn_sign`) and its prefix header; S3 confirms.
  Checked: libR on the laptop exports `Rf_lchoose` (nm), and Rmath.h declares it, so the R build's seventh name links as the
  six do.

## Amendments of 2026-11-05

Brought up to date with the work planned or built since 2026-11-03, on main 7937d252f; probes in scratch/fw-amend/probes on a
shipped build of that main.

| # | change | where it went |
|---|---|---|
| 1 | latents-first: a state carries no reading; design call 5 rewritten, the digest gone from the engine bullet, Context and S1 | The reading and where it lives; Calls made, 13; Context; S1; the design note's Calls made 5 and its What moves table (with the hurdle row) |
| 2 | latents-first: mutant 5 replaced by the count and nu constant at a mask change, with its draw-for-draw test; the state-pairing arm dropped | Mutants, 5; S1 gates |
| 3 | latents-first: `$setCounts` on a binomial sampler stores and rebuilds, drawing nothing | Engine: the binomial row (S3); S3 steps and gates; mutant 23 |
| 4 | latents-first: a row the weight composition reactivates draws nothing at the setter | Membership and zero weights (S4); RNG class; S4 gates |
| 5 | latents-first, the rest: the Student-t install list without a state install, nbinom and multinomial weight changes storing only, no setter draws as a constraint, a row out of the likelihood keeping its start, baselines measured against latents-first's re-recorded ones, the tests it retires | Engine (S1, S2); Constraints; RNG class; Context |
| 6 | hurdle-completion: the hurdle's resolution row is probit's, "precision" refused in its words, the resolved reading routed to both parts with a before/after arm | Resolution table; Refusals, 4; The R surface (S5); Calls made, 15; mutant 25; S5 |
| 7 | coefficient-prior-live: read against main; touches no weight code; two textual neighbors named; its "not verified" line closed | Interaction, coefficient-prior-live; Not verified |
| 8 | approximate-draws: shapes times w before its S2, after its S2 and after its S3; its S3 after this feature, its S1 and S2 either side; new shapes through the existing loops; gates at the sum's setting if S2 is first | Interaction, approximate-draws; Constraints; RNG class; Calls made, 17 |
| 9 | predict's weights on ordinal, nbinom and multinomial fits: the formal on all three, 1 a new observation, 0 NA on ordinal, nbinom's w as summed observations; refusals 19 and 20; the category case as an open call | The R surface (S5); Resolution table; Refusals; Calls made, 14; Open calls, 1; mutant 24; S5 |
| 10 | host-neutral S1 landed without `lchoose`: a seventh nmath.h name and choose.c in the standalone copy, in S3 | Engine: the binomial row (S3); Context; Calls made, 16; mutant 22; Not verified |
| 11 | S0 landed on main (1a23f9a60); the order is coefficient-prior-live, latents-first, this feature, host-neutral S2a | Status; front block; S0; Interaction |
| 12 | Summary rewritten for the above; slices re-estimated, ~99 agent-hours (150 to 200) | Summary; front block; Slices; Cost |
