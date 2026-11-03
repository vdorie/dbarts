# frequency-weights: weights read as counts of identical rows beside the precision default, and the binomial response

Status: PLANNED 2026-11-03 on main (3e4c002d8); amended the same day after a blind critique (scratch/fwp-crit/report.md, probes
p1 to p7 there), verdict BUILD AFTER AMENDMENTS, rebased onto main at 5a99172e6 (forest-door landed at a6a16c31a). Not built.
Rulings, quoted where they settle a call: dec-B125 (frequency weights join precision weights by explicit choice; precision
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
plan critique's p1 to p7, rerun here on its library of main (scratch/libs/fwp-crit) and 0.9-34's (scratch/libs/cran). No
prototype was built here.

agent: S0 one sonnet implementer (effort medium), its own small R-surface landing on main. The feature: S1, S2, S3 and S4 one
opus implementer each (effort high), S1 and S3 the riskiest; S5 and S6 one sonnet implementer each (effort medium); S7 and S8
sonnet. One opus review per engine slice (S1 to S4) that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior; one sonnet review for S0, S5, S6 and the gate slice S7
([Process by risk](README.md#process-by-risk)). Tier "Changes draws" for S1 to S4: new draws on new calls, engine read-offs
under frequency, and the creation and mask paths of every family rewritten, though no call taken today is meant to move;
tier "R surface" for S0, S5 and S6, S5 and S6 now pure R.
rng: neutral, bit for bit, for every engine draw of a call taken today: the equivalence trio, the four snapshot files,
tests/cpp and the exact gates do not move. Reported values that move are listed under [RNG class](#rng-class): the binary
xbart rmse (S0), a weighted logistic fit's training predictive draws (dec-B566), a 0 row's predictive draws on probit,
ordinal, hazard, hazard.probit and nbinom (NA, dec-B580), `$activeRows` after 0/1 weights (NULL), a logistic zero weight
(taken). New calls take new draws, each held to the expanded-data fit.
window: S0 now, on main. The feature for 1.0-0 (dec-B567: "before 1.0-0"), one branch, frequency-weights, merged once
(dec-B494), after coefficient-prior-live lands; its order against host-neutral-engine's main feature is a scheduling choice,
both orders stated under
[Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order).
budget: S0 ~60 lines. The feature ~4,000 lines, ~1,450 of them package code: engine ~575 (S1 ~285, S2 ~110, S3 ~120, S4 ~60),
bridge ~210, R ~650 (S4 ~80, S5 ~320, S6 ~250); help ~165, NEWS ~25; tests/cpp ~800, tinytest ~1,050, gate scripts ~450,
records ~60. Agent time, estimated and not measured: S0 ~3 with its review. The feature: S1 ~13, S2 ~7, S3 ~9, S4 ~7, S5 ~10,
S6 ~9, S7 ~6, S8 ~1, four engine reviews with mutants and fix rounds ~18, three R and gate reviews ~7, gates and sanitizers
~5, rebases ~2: ~94. Plans here run 1.5 to 2 times low, so 140 to 190.

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
observation like the row.

The binomial response is written as R writes it, cbind(successes, failures), under family = binomial or "logistic", and
xbart takes it too. Its fitted values are probabilities and its predictive draws counts, as glm's are. A proportion with the
trial counts as weights is refused with a message showing the cbind form.

One fix lands first and alone: xbart's rmse loss on a binary response compared the 0/1 response with the latent scale, in
0.9-34 as now, and will compare it with the probability instead.

As you ruled after the critique corrected the zero-weight premise, a row of weight 0 has missing predictive draws on every
family; probit, ordinal, the probit hazard and negative binomial fits drew values there, and its fitted value is still
reported.

What it costs. Nothing per sweep on the gaussian-type families. On the count families a sweep grows with the total count, as
logistic's does now: each trial of each row is about 160 nanoseconds a sweep, so 1,000 rows of 1,000 trials take about 0.16
seconds a sweep. The approximate-draw research you opened is where that changes; this plan keeps every draw exact.

What does not move. No posterior draw of any fit that runs today changes: the engine's random numbers, the four recorded
snapshot files and the equivalence baselines stay bit for bit, and each slice proves it. stan4bart, bartCause and treatSens
need no change; stan4bart can drop about 30 lines that become dead.

There are no open calls. The fix alone is about 60 lines; the feature about 4,000 lines over eight
slices, about 94 agent-hours as estimated, 140 to 190 by the record of plans here.

## Goal

After it lands: a frequency-weighted fit on the gaussian-type families equals the fit to the expanded rows draw for draw; the
Polya-Gamma families take whole-number counts with shape times w, each drawn as its adjacent copies are; a binomial row is
one logistic row of shape n w drawn as the two-row expansion is; probit and ordinal membership rides the data object; a 0
means the same on every family; the R surface of predict, fitted, residuals, extract, the log-likelihood and xbart reads both
readings and the binomial as the rulings state; and no call taken today moves an engine draw.

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
  `updateBasisScale`), the weights digest ([`SamplerBase::weightsDigest`](../../src/bartcore/facade.hpp)).
- Bridge: [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp),
  [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp), [`refuseCountOverCap`](../../src/R_interface_bartcore.cpp),
  [`validateResponseSupport`](../../src/R_interface_bartcore.cpp), the creation refusals of ordinal, nbinom and aft weights
  beside them, the weights slot read in the data conversion.
- R: [`enforceWeightPolicy`](../../R/spec.R), [`isMaskedWeightFamily`](../../R/spec.R),
  [`refuseNonCountWeights`](../../R/spec.R), [`refuseNonBinaryResponse`](../../R/spec.R), [`dbartsSpec`](../../R/spec.R) and its
  `active` element, [`estimateStartingSigma`](../../R/spec.R) over [`residualStandardError`](../../R/utility.R) and
  [`sparseResidualStandardError`](../../R/utility.R), [`restateMatrixResponseError`](../../R/data.R) and the classed
  `dbartsMatrixResponseError`, [`sampleFromPPD`](../../R/generics.R), [`pointwiseLogLikelihood`](../../R/generics.R),
  [`predict.bart`](../../R/generics.R), [`fitted.bart`](../../R/generics.R), [`residuals.bart`](../../R/generics.R),
  [`extract.bart`](../../R/generics.R), [`xbartLossFunction`](../../R/xbart.R), [`xbartSigestFunction`](../../R/xbart.R) and the
  0/1 refusal in [`xbart`](../../R/xbart.R), [`refuseCountsMutation`](../../R/bartcore.R), the dbartsModel class in
  [R/A_class.R](../../R/A_class.R), `$setModel`'s refusals in [R/fittingLevel.R](../../R/fittingLevel.R) (forest-door's), the
  sampler's `$setWeights`, `$setData`, `$setCounts`, `$setResponse` and `$setActiveRows` ([R/dbarts.R](../../R/dbarts.R)).
- Help: man/bart.Rd, man/dbarts.Rd, man/dbartsSpec.Rd, man/xbart.Rd, man/dbartsSampler-class.Rd (`setData`'s and
  `setWeights`' mask sentences), man/dbarts-embedding.Rd (its "survey weights" row), man/samplePriorPredictive.Rd (its
  weighted logistic sentence), man/predict.bart.Rd.
- Tests that pin today's weight behavior: test-binary-weight-mask.R, test-weighted-logistic.R, test-weighted-binary-ppd.R,
  test-logistic-weight-swap.R, test-zero-weights.R, test-xbart-weight.R, test-state-weight-pairing.R,
  test-state-zero-weight-rows.R, test-pointwise-loglik.R, test-nbinom.R, test-xbart-loss.R.
- TODO items: frequency weights (dec-B125), xbart-binary-rmse (S0), mask-on-own-data-object (closed by S4),
  approximate-distributions (the binomial's cost case, not this plan's).

## Constraints

- No draw of a call taken today moves: every engine path a precision fit, an unweighted fit or a binary logistic fit runs
  keeps its expression and its order of random numbers. Where a new term is zero under precision, it is not added, rather
  than added as 0.0.
- Every draw stays exact; real weights on the Polya-Gamma families are refused until a real-shape draw exists (dec-B568).
- No row out of the likelihood carries a zero precision or a non-finite working response: the node kernels read
  `omega` and `working / omega`, and 0 there is a NaN.
- Out of scope: frequency weights on probit, ordinal, the probit hazard and aft (dec-B563); hurdle.lognormal's weights
  (dec-A15); the probit binomial (dec-B567); any approximate draw (TODO approximate-distributions); bartBT and rbart_vi, which
  take no weights.type; the shipped header (its field rides host-neutral-engine's S3, below).
- No new dependency; base R calls within DESCRIPTION's R floor.

## Design

### The reading and where it lives

- R: dbartsModel gains a slot `weights.type`, a character of length one, NA for the family's own reading (weights.type = NULL
  as given), else "precision" or "frequency". bart, dbarts, dbartsSpec and xbart take the argument and write the slot; the
  spec's resolution, where the family is known, writes the resolved value, so dbartsSpec's returned model carries it and
  creation reads it. A model made by hand with NA resolves at creation the same way.
- Bridge: reads the resolved slot into a sampler-level flag, `bartcore::WeightsReading { precision, frequency }`, passed to the
  response family's constructor beside the weights.
- Engine: the flag is held by the families that read it and enters the weights digest, so a state stored under one reading and
  installed into a sampler of the other redraws as a weight change does (design Calls made, 5).

Resolution and what each value takes, per family (the refusal numbers are [Refusals](#refusals)'):

| family | NULL resolves to | "precision" | "frequency" | weights.test under the resolved reading |
|---|---|---|---|---|
| gaussian, Student-t, heteroscedastic | precision | taken, real w | taken, real w | precision: as today; frequency: 0 and 1 only (9) |
| logistic, hazard.logistic | frequency (counts) | refused (3) | taken, whole w (5) | trial counts (dec-B470), as today |
| nbinom, multinomial | frequency (counts) | refused (3) | taken, whole w (6, 7) | trial counts, as today |
| probit, ordinal, hazard, hazard.probit | frequency (membership) | refused (3) | 0 and 1 only (2) | 0 and 1, as today |
| aft, hurdle.lognormal | precision (weights refused) | weights refused (4) | refused (4) | as today |

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
  over rows with m_i > 0 only (at m_i = 0 the term is -lgamma(0), infinite), computed when the weights install (creation,
  `setWeights`, `setData`, a mask change, a state install). Under precision `drawIndex` runs its shipped expression, the
  constant never built. A row at weight 0 keeps the lambda it draws, as now, out of the composite.
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

### Engine: the binomial row (S3)

- `LogisticResponse` gains an optional `trials_` pointer, null for a binary response, where it reads as 1. With trials, `y_`
  holds the successes k (for a binary row k = y, n = 1, the same values). The shape is `lround(w n)` unit draws (w the
  frequency weight, 1 absent), the working response `w (k - n / 2) / omega - offset`, written so that at null trials the
  expression is today's `w (y - 0.5) / omega - offset`, operand for operand. The cold start is `omega = 0.25 max(n w, 1)` and
  the working response `4 (k / max(n, 1) - 1/2) - offset`, so a row of zero trials or zero weight, out of the likelihood under
  the mask, keeps a positive omega and a finite working response.
- The log-likelihood: `w [k log p + (n - k) log(1 - p)] + w c_i` with `c_i = log choose(n_i, k_i)`, per-row constants computed
  once per install of the counts with R's `lchoose` (`Rf_lchoose` on the R build, the corrected lbeta route), so a binomial row
  reports `w log dbinom(k; n, p)` and a binary row today's value. `std::lgamma` is not used: at n = 1e6, k = n / 1000 its
  constant is 2.3e-10 off in relative terms against `dbinom`, R's route 2.3e-12 (plan critique p4, rerun here).
- The bridge hands successes `counts[, 1]` and trials `rowSums(counts)` from the data object's counts slot, the slot multinomial
  count rows use (design Calls made, 9); `validateResponseSupport` takes logistic successes above 1 where trials are given, up
  to the trials; `$setCounts` on a binomial sampler swaps both and redraws as a count swap does now. A row of 0 trials installs
  into the mask as dec-B133's zero-trial row does, with its warning.
- Draw for draw: a binomial row of k successes in n trials consumes the n unit draws at the same psi that two adjacent rows of
  weights k and n - k consume, and its kappa is their kappa's sum in exact arithmetic (`k - n/2` against `k/2 - (n - k)/2`,
  each exact for whole numbers below 2^52). The omega sums differ only in the association of n additions, so the gate is a
  tolerance of 1e-10 on the fits, not bitwise.
- The spec's form (for host-neutral-engine's S2b): `y` the successes and a `trials` field, per row, NULL for a binary
  response, under LOGISTIC only; the build refuses `trials` on any other family.

### Membership and zero weights (S4)

- `ProbitResponse` and `OrdinalResponse` take the caller's 0/1 weights at creation, `setWeights` and `setData` and compose them
  with the active-row mask into the effective mask, as `GaussianResponse::recomposeActiveRows` composes weights and mask: a row
  is in when its weight is 1 and the mask leaves it in. A row the composition reactivates redraws its latent as a mask change
  does now. A weight other than 0 or 1 there keeps today's refusal; the bridge's creation refusals of probit and ordinal
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
  logistic's training draw or its zero weight: logistic is not in 0.9-34. S0's NEWS line is its own.

## S0: the binary xbart rmse, alone on main

A defect in 0.9-34 as on main (TODO xbart-binary-rmse): on a binary response `xbartLossFunction`'s rmse averages the latent
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
4. aft and hurdle.lognormal keep `aft (survival) models do not support case weights` and the hurdle refusal for weights;
   "frequency" there, with or without weights: `weights.type = "frequency" is not taken under family = "<family>", which takes
   no weights`.
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

Gone: xbart's `xbart does not accept weights of 0 and 1 under family ...`; logistic's refusal of a 0; the logistic refusal of
a two-column response under family = binomial and "logistic"; the multinomial refusal of weights other than all ones; the
bridge's creation refusal of 0/1 weights on probit and ordinal.

## RNG class

Neutral on the engine for every call taken today, shown three ways:

- Bitwise harnesses, run against each engine slice's own `--preclean` library: `benchmarks/R/equivalence.R compare`,
  `bcf-equivalence.R compare` and `multinomial-equivalence.R compare`, every scenario "identical draws (same RNG stream)",
  counted per [Gate hygiene](README.md#gate-hygiene); the four `test-reproducibility-*.R` on the reference build, unchanged
  (two of them carry precision weights, which pins the precision path); tests/cpp unchanged in every existing test.
- Construction: under precision every changed expression is today's, operand for operand (the count is
  `static_cast<double>(numPositiveWeights_)`, the variance leaf adds 1.0, the nu constant is not built, the moments and the
  basis scale read unit weights, a null `trials_` reads as today's working response); each slice's reviewer reads the diff for
  it and runs the mutants that would break it.
- The membership paths (S4): a probit and an ordinal fit with 0/1 weights equal today's draw for draw, the mask installed by
  the engine at creation against today's install by R after creation, which draws no variate (no row reactivates; the critique
  checked that `ProbitResponse` draws nothing at creation); an nbinom fit with 0/1 weights equals the documented by-hand route,
  `new("dbartsSampler", ...)` then `$setActiveRows(w, updateScale = TRUE)`, on the same build, a live oracle; pinned by
  test-binary-weight-mask.R and test-nbinom.R unchanged in their draws.

Reported values that move, each re-pinned in the test that holds it:

- S0: a binary xbart's rmse loss moves to the probability scale; its default, log, does not move. A defect fix against 0.9-34,
  not a ruling's move; in NEWS.
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
of S7. The exact-gates list grows by one script (S7); no existing gate's arm changes.

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

## Open calls

None. Every user- or consumer-visible choice the plan makes is settled by dec-B563 to dec-B570, dec-B573, dec-B574, dec-B580
or the design note's calls made; the plan's own calls above are internal or follow from a ruling (Calls made, 5, from
dec-B569). dec-B580 settled the one premise the critique found false.

## Slices

Each slice is a commit series on frequency-weights, its own review and gates; the feature merges once after S8. S0 is apart
(above) and lands first.

### S1: engine, the gaussian-type families and their read-offs (~285 code, ~350 tests/cpp, ~80 bridge and R plumbing)

Tier "Changes draws"; RNG class neutral on every call taken today, new draws under frequency.
Steps: `WeightsReading` and its bridge read from the model slot (the slot and its R plumbing land here, refusals in S5); the
gaussian count pointer and `rowCount_`; `drawSigmaSqFromPosterior` taking a double count; Student-t's frequency draw, its
counts, the nu constant over m_i > 0; the variance leaf's count and every site that counts its rows; the weighted leaf
covariate moments and default forest sd with their setter doors; the gaussian and Student-t log-likelihood channel; the digest
covers the reading.
Gates: tests/cpp against w adjacent copies, draw for draw (sigma and the fits within 1e-12 over 300 sweeps on gaussian;
heteroscedastic, its variance leaves additive, the same bound; the critique's p6 gave 8e-15 at a fixed sigma); Student-t
against the copies in distribution and nu's grid probabilities on fixed residuals (the note's p8 arm, within Monte Carlo
error); a linear leaf and a two-forest fit against the copies draw for draw, covering the moments and the default sd; a 0
weight under frequency on each of the three families, finite throughout; the w = 1 identity (frequency at all-ones weights
equals precision bitwise); a state, copy and reload round trip of a frequency sampler; the log-likelihood channel against
`dnorm` and `dt` times w; the bitwise trio; the full tinytest suite; sanitizers on tests/cpp and on the R-loaded path;
`mutation-battery.R verify-anchors`.

### S2: engine, the Polya-Gamma families (~110 code, ~250 tests/cpp, ~40 bridge)

Tier "Changes draws"; neutral today (nbinom and multinomial weights were refused or masks).
Steps: nbinom's shape, kappa, r-step power, weighted center and cold start; the combiner's multinomial shape and kappa; the
caps at w > 1; the bridge stops refusing nbinom and multinomial weights.
Gates: tests/cpp: an nbinom row of weight w against w adjacent copies draw for draw, at a fixed r and at a drawn r (the r
draw reads the same histogram, the unit draws are consumed in row order), within 1e-10; a multinomial row of weight w against
w adjacent copies (unit draws in row order inside each category's loop) and against the counts times w, draw for draw; a 0
weight on each family, finite throughout; y at 999,999 with w = 1 taken, the product cap at w = 2 refused; a state round trip;
the trio; the suite; sanitizers.

### S3: engine and bridge, the binomial row (~120 code, ~200 tests/cpp, ~60 bridge)

Tier "Changes draws"; neutral for a binary response.
Steps: `trials_`, the working response, the cold start, the log choose constants and the log-likelihood; the bridge reads the
counts slot under logistic and `validateResponseSupport` admits successes up to the trials; `$setCounts` reaches a binomial
sampler; the zero-trial row into the mask.
Gates: tests/cpp: a binomial row against two adjacent rows of weights k and n - k for 0 < k < n, fits within 1e-10 over 300
sweeps, the same number of variates consumed; a binomial row at frequency weight w against the row of k w successes in n w
trials, bitwise; the log-likelihood against R's `dbinom(k, n, p, log = TRUE)` times w within 1e-10 relative, at n up to 1e6
with k = n / 2 and k = n / 1000 named; a zero-trial row and a zero-weight binomial row finite throughout; a state, copy and
reload round trip of a binomial sampler; a null-trials sampler bitwise today's; the trio; the suite. Today a weight-0 logistic
row still draws one variate, so k = 0 and k = n join the draw-for-draw test in S4, once the expansion's weight-0 row is masked
and draws none.

### S4: membership and zero weights (~60 engine, ~30 bridge, ~80 R, ~200 tinytest)

Tier "Changes draws" (the mask paths of four families change); neutral on every call taken today.
Steps: probit and ordinal take weights and compose; the bridge's creation refusals lifted for 0/1; logistic takes a 0;
nbinom's 0/1 through S2's counts; the R policy keeps the weights slot; `dbartsSpec()` drops `active`; xbart's refusal goes;
help for `$setWeights`, `$setData`, `$setActiveRows`, dbartsSpec.
Gates: the binomial draw-for-draw test of S3 at k = 0 and k = n; test-binary-weight-mask.R's draws unchanged (its
`s$data@weights`, `$activeRows` and `spec$active` assertions re-pinned to the ruling); nbinom with 0/1 weights bitwise against
`new("dbartsSampler", ...)` then `$setActiveRows(w, updateScale = TRUE)` on the same build; a new test that `s$setData(s$data)`
keeps masked rows out and that a mask set by `setActiveRows` composes with the weights in either order; a logistic fit with a
zero weight equal draw for draw to the fit with that row masked by `setActiveRows`; xbart on probit with 0/1 weights equal,
fold for fold under a fixed fold assignment, to a by-hand loop fitting each fold with those weights and scoring its held-out
rows of weight 1 (the test-xbart-fold-oracle.R pattern); the suite; the trio.

### S5: R surface, the reading (~320 R, ~110 help, ~450 tinytest)

Tier "R surface", pure R. Moves the weighted logistic training draw and, by dec-B580, a 0 row's draws on probit, ordinal,
hazard, hazard.probit and nbinom; re-pinned.
Steps: weights.type on bart, dbarts, dbartsSpec and xbart, its validation and the resolution table; Refusals 1 to 10; the
starting sigma on both routes and xbart's per fold, the chi anchor; `sampleFromPPD` (one observation, NA at a 0 on every
family, dec-B580, the fitted value kept), predict's and weights.test's doors, `pointwiseLogLikelihood`; help and NEWS.
Gates: the touched test files and the new test-frequency-weights.R (each refusal by its words, counted warnings; the
resolution table row by row; sigma's starting value and the chi anchor against the expanded data on the dense and the sparse
routes and in xbart's folds; predict's and weights.test's doors; NA draws at a 0 on every family with the stream unmoved;
`pointwiseLogLikelihood` against the engine's channel); lint_package, air, codoc, win drift, doc-freshness; the snapshots on
the reference build.

### S6: R surface, the binomial and xbart (~250 R, ~55 help, ~400 tinytest)

Tier "R surface", pure R. Builds on S0.
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
5. S1: the digest ignores the reading - test-state-weight-pairing.R's new arm (a state under frequency into a precision sampler
   redraws).
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

## Interaction with the other pending plans, and order

- forest-door landed (a6a16c31a), with setmodel-parameters: its forest records carry shares, `$setModel` gained refusals in
  R/fittingLevel.R beside which Refusal 10 sits, and chain.hpp, combiner.hpp, model.hpp and R_interface_bartcore.cpp moved; this
  plan's symbols were re-read on 5a99172e6 and none it names changed role. coefficient-prior-live lands next (it reworks
  chain.hpp, combiner.hpp, model.hpp, R_interface_bartcore.cpp, R/dbarts.R and R/spec.R); engine features merge serially, so
  this feature branches from main after it or rebases onto it, re-gating what the rebase touches.
- S0 touches only R/xbart.R, its help, NEWS and one test file, and lands whenever green.
- host-neutral-engine: its S1 (build-time math, branch host-neutral-math) is independent, except that this plan's S3 calls R's
  `lchoose` (with `lbeta` and the `lgammacor` correction under it), which the standalone build's copied math set then needs;
  whichever lands second adds it. Its main feature rewrites every creation route's SEXP reads into a plain-C spec (S2b) that
  carries `weightsType` per dec-B564 and, under this plan, `trials` beside y under LOGISTIC. The order is a scheduling choice:
  - This feature first: host-neutral's S2b translates the weights reading, the counts slot under logistic and the membership
    rule as code that exists, and its S3 adds the header field. Cost: host-neutral's S2a rebases over this feature's bridge
    changes (~210 lines in the creation and mutation conduits it rewrites), and S2a waits for this feature to merge.
  - Host-neutral first: this feature rebases its bridge reads (the reading flag, the trials, the membership weights, the
    lifted refusals) onto S2b's spec form, writing them as spec fields and build checks, and adds `trials` to the header
    itself, a second hash move before the release, which rebuilds no one. Cost: roughly the same ~210 lines rewritten in the
    other direction, and this feature waits for host-neutral's merge.
- approximate-draws (being designed): the binomial at large n is its cost case. Today's exact cost, kept here: a row of n
  trials at frequency weight w takes `n w` unit Polya-Gamma draws a sweep, about 160 ns each (164 ns at psi 0, 159 ns at psi 3,
  p12cost.R on the shipped build), so 1,000 rows of 1,000 trials cost about 0.16 s a sweep and a 1,000-sweep chain about
  160 s, linear in the total trials; the tree pass adds what a row of trial 1 costs. The same holds for logistic counts,
  nbinom's `w (y + r)` and multinomial's `w n_i` per category. This plan keeps the shape a whole number and every draw exact;
  any approximate or real-shape draw replaces the unit-draw loop in `LogisticResponse::refreshLatents`, the nbinom omega draw
  and the combiner's loop, and is that arc's.
- weighted-binary (parked, after 1.0-0): its integer-weight probit is dec-B563's truncated-normal case, not scheduled; its
  real-weight logistic is the real-shape draw above.
- Nothing else touching src/bartcore/model.hpp's response families merges while S1 to S4 are in review.

## Cost

S0 about 60 lines and 3 agent-hours. The feature about 4,000 lines, 1,450 of them package code; about 94 agent-hours
estimated, 140 to 190 by the record of plans here. No per-sweep cost on the gaussian-type families at either reading (the
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
- How coefficient-prior-live leaves model.hpp's gaussian family and the combiner; its branch was read for the files it touches,
  not line by line.
- That `Rf_lchoose` is reachable from model.hpp under the engine's math layer as the other `Rf_` densities are; S3 confirms.
- stan4bart's, bartCause's and treatSens's call sites were read on main 27ef163, 1e4b882 and master 268b2db; their 1.0 compat
  branches may differ.
