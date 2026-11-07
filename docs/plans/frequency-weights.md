# frequency-weights: weights read as counts of identical rows beside the precision default, and the binomial response

Status: PLANNED 2026-11-03 on main (3e4c002d8); amended the same day after a blind critique (scratch/fwp-crit/report.md, probes
p1 to p7 there), verdict BUILD AFTER AMENDMENTS, rebased onto main at 5a99172e6 (forest-door landed at a6a16c31a), and onto
main at 7937d252f on 2026-11-05; amended 2026-11-05 for the work planned or built since
([Amendments of 2026-11-05](#amendments-of-2026-11-05)): latents-first, hurdle-completion, approximate-draws,
coefficient-prior-live and host-neutral-engine's landed S1; rebased onto main at 6f84b5a66 (latents-first, hurdle-completion
and coefficient-prior-live landed) and amended 2026-11-06 for dec-B590 to dec-B594
([Amendments of 2026-11-06](#amendments-of-2026-11-06)): its open call ruled, two left by that ruling put as Open calls 2
and 3. S0 landed on main (1a23f9a60). The feature is not built.
Rulings read for the 2026-11-06 pass: dec-B594 (predict's weights on ordinal and multinomial fits are trials, a draw with
them the count per category; dec-B177's predict clause revised; "Use your revised recommendation."), dec-B590 (a
Polya-Gamma row above a shape of 64 or at a fractional shape takes the accept-or-keep route once proven, the sum otherwise),
dec-B591 (that route proven by independent agents before it is built), dec-B592 (every sampler of that arc in 1.0-0, the
approximate ones opt-in, the default from a study before 1.0-0), dec-B593 (the sampler chosen by a method name, an
approximate one's accuracy by a separate tolerance). Rulings since the first plan, read whole for the 2026-11-05 pass: dec-B586 (a sweep draws its augmentation variables first; no setter draws them),
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
"Use (a)."), with dec-B133 and dec-B177 (zero-trial rows; its predict clause revised by dec-B594), dec-B202 (the one-million cap), dec-B470 and dec-B472 (predict's
weights), dec-B494 (a feature merges once), dec-B85 (the header change frequency weights' field rides). The runs behind it:
the note's probes (scratch/fw-note, the first critique's p1 to p11 under scratch/fw-crit/probes), p12cost.R in scratch/fw-note
(a unit Polya-Gamma draw at 164 ns, psi 0, and 159 ns, psi 3, on the shipped build of c3d268c8f on the shared laptop), and the
plan critique's p1 to p7, rerun here on its library of main (scratch/libs/fwp-crit) and 0.9-34's (scratch/libs/cran), and the
amendment's probes in scratch/fw-amend/probes on a shipped build of main 7937d252f (scratch/libs/fw-amend); the 2026-11-06 pass read the code on main 6f84b5a66 and ran
nothing. No prototype was built here.

agent: S0 landed (one sonnet implementer and review, on main at 1a23f9a60). The feature: S1, S2, S3 and S4 one opus
implementer each (effort high), S1 and S3 the riskiest; S5 and S6 one sonnet implementer each (effort medium); S7 and S8
sonnet. One opus review per engine slice (S1 to S4) that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior; one sonnet review for S5, S6 and the gate slice S7
([Process by risk](README.md#process-by-risk)). Tier "Changes draws" for S1 to S4: new draws on new calls, engine read-offs
under frequency, and the creation and mask paths of every family rewritten, though no call taken today is meant to move
but a weighted logistic fit with a linear leaf or a numeric basis once S5 resolves logistic to counts (Calls made, 21);
tier "R surface" for S5 and S6, both pure R.
rng: neutral, bit for bit, for every engine draw of a call taken on the build it starts from, main 6f84b5a66, after
coefficient-prior-live and latents-first (whose draws latents-first moved and re-recorded once): the equivalence trio
(equivalence-a6e83f38.rds, bcf-equivalence-a6e83f38.rds, multinomial-equivalence-80b1c8d4.rds), the four snapshot files,
tests/cpp and the exact gates do not move. Reported values that move are listed under [RNG class](#rng-class): a weighted
logistic fit's training predictive draws (dec-B566), a 0 row's predictive draws on probit, ordinal, hazard, hazard.probit
and nbinom (NA, dec-B580), `$activeRows` after 0/1 weights (NULL), a logistic zero weight (taken), and a multinomial count
fit's predictive draws at new rows, a one-trial count vector where they drew a category (dec-B594). New calls take new
draws, each held to the expanded-data fit or, at predict, to its stated distribution.
window: the feature for 1.0-0 (dec-B567: "before 1.0-0"), one branch, frequency-weights, merged once (dec-B494). Build order,
the orchestrator's call as latents-first's plan states it: coefficient-prior-live and latents-first (both landed), this
feature, then host-neutral-engine S2a onward; hurdle-completion landed before it; see
[Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order).
budget: the feature ~4,450 lines, ~1,610 of them package code: engine ~570 (S1 ~285, S2 ~110, S3 ~125, S4 ~50), bridge ~200,
R ~815 (S4 ~80, S5 ~485, S6 ~250), standalone math ~25 written (S3; R's choose.c copied unedited, not counted); help ~200,
NEWS ~30; tests/cpp ~820, tinytest ~1,280, gate scripts ~450, records ~60. Agent time, estimated and not measured: S1 ~12.5,
S2 ~7, S3 ~9.5, S4 ~6.5, S5 ~15.5, S6 ~9, S7 ~6, S8 ~1, four engine reviews with mutants and fix rounds ~18, three R and
gate reviews ~8, gates and sanitizers ~5, rebases ~3: ~101. Plans here run 1.5 to 2 times low, so 150 to 200. The
recommendations of Open calls 2 and 3 are counted; their other options change S5 by under half an hour.

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
three. The negative binomial reads a whole number as that many new observations summed and refuses 0, as logistic reads its
weights as trials and refuses 0. On ordinal and multinomial fits, as you ruled on 2026-11-05, predict without weights draws
in the shape of the fit's response (a category, or a count vector on a fit to count rows, one trial at a new row), and with
weights the count in each category over that many new trials; so weights = 1 on a fit to categories returns a one-trial
count vector where predict() returns the category. Two choices that ruling leaves are yours, below: what weights = 0 gives
there, and what a test sample's weights draw on those two families.

The xbart fix that was to land first has landed on main. Since this plan was written you also ruled that every sweep draws
its hidden values first, so no setter draws them; this plan is built after that change and is simpler for it: nothing here
redraws when weights, counts or rows change, the next sweep does. The approximate-draw work you ruled on 2026-11-05 is where
the cost below changes and where fractional counts on logistic and binomial fits arrive: once its exact route is proven and
built, a row whose count here is above 64 takes it. Until then, and in every gate here, this plan draws as today, exactly and
in whole numbers; nothing in it waits on that proof.

What it costs. Nothing per sweep on the gaussian-type families. On the count families a sweep grows with the total count, as
logistic's does now: each trial of each row is about 160 nanoseconds a sweep, so 1,000 rows of 1,000 trials take about 0.16
seconds a sweep.

What does not move. No posterior draw of any fit that runs on the build this starts from changes: the engine's random
numbers, the four recorded snapshot files and the equivalence baselines stay bit for bit, and each slice proves it. stan4bart,
bartCause and treatSens need no change; stan4bart can drop about 30 lines that become dead.

Three open calls, below, the third raised in S3. The feature is about 4,450 lines over eight slices, about 101 agent-hours as estimated, 150 to 200
by the record of plans here.

## Goal

After it lands: a frequency-weighted fit on the gaussian-type families equals the fit to the expanded rows draw for draw; the
Polya-Gamma families take whole-number counts with shape times w, each drawn as its adjacent copies are; a binomial row is
one logistic row of shape n w drawn as the two-row expansion is; probit and ordinal membership rides the data object; a 0
means the same on every family; the R surface of predict, fitted, residuals, extract, the log-likelihood and xbart reads both
readings and the binomial as the rulings state, predict's weights reaching the ordinal, negative binomial and multinomial fits
as dec-B580's door asks, as trials on ordinal and multinomial fits (dec-B594); and no call taken on the build it starts from moves an engine draw
but a weighted logistic fit with a linear leaf or a numeric basis, whose read-offs then read the counts (Calls made, 21).

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
  `updateBasisScale`). The weights and survival digests, the Student-t zero-weight record, setState's reconciliation and
  every setter draw retired with latents-first (dec-A263); none of `weightsDigest`, `reapplyWeights` or `redrawEntering`
  remains in src or R on main 6f84b5a66, and a sweep draws its augmentation variables first in `Chain::drawAugmentation`
  ([chain.hpp](../../src/bartcore/chain.hpp)).
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
  coefficient-prior-live's coefficient refusal beside them), the hurdle arm of bart and its weights refusal
  ([`refuseHurdleWeights`](../../R/bart.R), hurdle-completion's) and [`predict.bartHurdle`](../../R/generics.R)'s 0/1
  weights, [`multinomialPpdFromProbs`](../../R/generics.R) and [`multinomialCountPpdFromProbs`](../../R/generics.R) (the
  category draw and the count vector draw, shared by extract and predict), the
  sampler's `$setWeights`, `$setData`, `$setCounts`, `$setResponse` and `$setActiveRows` ([R/dbarts.R](../../R/dbarts.R)).
- Help: man/bart.Rd, man/dbarts.Rd, man/dbartsSpec.Rd, man/xbart.Rd, man/dbartsSampler-class.Rd (`setData`'s and
  `setWeights`' mask sentences), man/dbarts-embedding.Rd (its "survey weights" row), man/samplePriorPredictive.Rd (its
  weighted logistic sentence), man/bartBT.Rd (predict.bart's page) and man/bart.Rd's ordinal, negative binomial and
  multinomial generics paragraphs (their predict is documented there; the multinomial paragraph carries dec-B177's
  one-trial category sentence that dec-B594 revises). latents-first rewrote the sampler's setter sentences.
- Tests that pin today's weight behavior: test-binary-weight-mask.R, test-weighted-logistic.R, test-weighted-binary-ppd.R,
  test-logistic-weight-swap.R, test-zero-weights.R, test-xbart-weight.R, test-pointwise-loglik.R, test-nbinom.R,
  test-xbart-loss.R; test-state-weight-pairing.R as latents-first left it (its digest arms now check that no digest
  attribute is stored; test-state-zero-weight-rows.R deleted, its undo and copy arms moved there). Tests that pin what
  dec-B594 moves: test-predict-training-rows.R (predict's refusal of `weights` by name on nbinom and ordinal fits) and
  test-multinomial-surface.R's count-matrix arm (`extract(sample = "test")` equal to predict at newdata under one seed,
  which holds at the new shape while the draws themselves move).
- Plans it follows or meets: docs/plans/latents-first.md (Order and overlaps; landed), docs/plans/hurdle-completion.md
  (Interaction with frequency-weights, and order; landed), docs/plans/approximate-draws.md (Real shapes, Setters and marked
  rows; on its branch at e0c62e566, its open calls ruled as dec-B590, dec-B592 and dec-B593),
  docs/plans/host-neutral-engine.md (on its branch).
- TODO items: frequency weights (dec-B125), mask-on-own-data-object (closed by S4), hurdle-weight-extensions (its logistic zero
  part takes counts once this lands), approximate-distributions (the binomial's cost case, not this plan's); xbart-binary-rmse
  closed by S0.

## Constraints

- No draw of a call taken today moves: every engine path a precision fit, an unweighted fit or a binary logistic fit runs
  keeps its expression and its order of random numbers. The one exception is a weighted logistic or hazard.logistic fit
  with a linear or gp leaf or a numeric forest basis: its weights are counts (dec-B470), so once S5 resolves its NULL to
  frequency, its read-offs are the expanded rows' and its draws move (Calls made, 21; [RNG class](#rng-class)). Both
  features are new in 1.0-0; 0.9-34 had neither. Where a new term is zero under precision, it is not added, rather
  than added as 0.0.
- Every draw stays exact; real weights on the Polya-Gamma families are refused until a real-shape draw exists (dec-B568),
  which approximate-draws' S3 brings for logistic, binomial and hazard.logistic (below).
- Each new Polya-Gamma shape (nbinom's `w (y + r)`, the binomial's `w n`, multinomial's `w n_i`) is formed and handed to the
  unit-draw loop the family already has; no slice adds a loop of its own, so approximate-draws' route, chosen per row
  (dec-B590), meets one loop per family whichever lands first.
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
- Engine: the flag is a sampler option, `SamplerOptions::weightsReading`, fixed at creation. The gaussian and Student-t
  responses take it at construction and hold it; the read-offs (a leaf covariate's center and scale, a gp leaf's default
  lengthscale, a forest's default sd, the multi-forest anchor) read it from the sampler and chain options, since the sampler
  and chain compute them, for every family alike (Calls made, 22). A state carries no reading: after latents-first a state carries no
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
| multinomial | frequency (counts) | refused (3) | taken, whole w (7) | trials: whole w > 0 a count per category over w new trials (dec-B594), 0 [Open calls](#open-calls), 2; new at predict; weights.test [Open calls](#open-calls), 3 |
| probit, hazard, hazard.probit | frequency (membership) | refused (3) | 0 and 1 only (2) | as today: trial counts at predict (dec-B470); weights.test 0 and 1 |
| ordinal | frequency (membership) | refused (3) | 0 and 1 only (2) | trials: whole w > 0 a count per category over w new trials (dec-B594), 0 [Open calls](#open-calls), 2; new at predict; weights.test [Open calls](#open-calls), 3 |
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
  scale (`standardizationMomentsForColumn`, weighted moments over the rows in) and its `updateStandardization` door; a gp
  leaf's default lengthscale, the median pairwise distance over the expanded rows' pairs (Calls made, 20); a forest's
  default sd (`basisColumnScale`, weighted, with or without a mask) and its `updateBasisScale` door; the multi-forest
  anchor, the response's sd over the rows read. The anchor's sd is 0 exactly, which takes the no-spread anchor, where the
  rows read hold one value, as the expansion's, or where their counts total no more than 1 to the rounding of their sum
  (weights rescaled to sum to 1 sum to 1 plus a few units in the last place); the basis scale and the leaf moments take
  no scale on the same total. Where too few rows are in and a read-off falls back to every row, it reads
  every row at its count, as the expansion's every row (Calls made, 23). `gpMaxLeafSize` counts rows, the one stated
  exception (Calls made, 20). The response range stays over the rows in, which
  copies do not move. The quantile cut rule is not weighted: its grid is built over unique values, which copies leave as they
  are (plan critique p5: identical grids, original against expanded, with and without useQuantiles).
- The log-likelihood channel: `GaussianResponse::computeLogLikelihood` and `TResponse::computeLogLikelihood` report
  `w log f(y_i | one observation at sigma)` under frequency, the unit-weight density times w; under precision as now.
- A frequency weight of 0 is a row out of the likelihood, as now; a positive real w is taken (dec-B568).

### Engine: the Polya-Gamma families (S2)

- logistic, hazard.logistic: their weights already are counts; the slice only checks the reading flag (precision refused
  above). The engine change in the family is zero, the zero weight being S4's. The read-offs S1 weights under frequency
  reach logistic once S5 resolves its NULL to frequency, so a weighted logistic fit with a linear or gp leaf or a numeric
  basis moves ([RNG class](#rng-class); Calls made, 21).
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
  count kernel and composite, as latents-first's setActiveRows row has them, the center held unless the caller asks for a
  rescale, Calls made 24); it draws nothing, the next sweep drawing
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
  - nbinom, a count family (dec-B566): a positive whole w draws the sum of w new counts at the row's mean mu and shape r,
    one draw from NB(size w r, mean w mu) per posterior draw, as a logistic fit's w draws Binomial(w, p); 0 and a fraction
    refused in `refuseNonCountWeights`' words, as logistic's predict refuses them (dec-B470: "positive integer trial counts";
    run on main: "the posterior predictive 'weights' of a probit fit are trial counts and must be positive integers").
  - ordinal and multinomial, ruled by dec-B594: without weights predict(type = "ppd") "draws in the shape of the fit's
    response", and "with weights, a draw is the count per category over that many new trials on every one of these fits",
    each draw Multinomial(w, p(x)) at a posterior draw of the category probabilities, "so nothing well-defined is refused".
    The maintainer, 2026-11-05: "What's a correct draw from the posterior predictive distribution here? I always thought it
    was well-defined.", then "Use your revised recommendation." It revises dec-B177's predict clause. What predict returns
    (S draws, n rows, K categories; under combineChains = FALSE a chain margin leads, as today):

    | fit and rows | without weights | with weights w, a whole number above 0 |
    |---|---|---|
    | ordinal, training rows or newdata | one category per draw and row, S x n codes, as today | a count per category over w trials, S x n x K |
    | multinomial to one category per row, training rows or newdata | one category, S x n, as today | the same S x n x K counts |
    | multinomial to count rows, training rows (no newdata) | a count vector at the row's trials n_i, S x n x K (dec-B177), as today | counts over w trials, w in place of n_i (dec-B470's replace rule as dec-B566 reads it on a count family) |
    | multinomial to count rows, newdata | a count vector at one trial, S x n x K, where today it is one category (dec-B594) | counts over w trials |

    predict(weights = 1) on a fit to categories returns a one-trial count vector, each a unit vector of K, where predict()
    returns the category, "the same draw written two ways" (dec-B594); it is dec-B580's door to a 0 training row, whose
    draws without weights are NA. weights = 0 is [Open calls](#open-calls), 2 (recommended: NA draws); a fraction, a
    negative value or Inf is refused (Refusals, 20). The count vectors are drawn by `multinomialCountPpdFromProbs`, the
    categories by `multinomialPpdFromProbs` as today (Calls made, 18); ordinal's probabilities are its type = "ev" category
    probabilities, from the latent and the thresholds. A count fit's test sample draws as predict at those rows, a one-trial
    count vector without weights.test (Calls made, 19); what weights.test draws on either family is
    [Open calls](#open-calls), 3. 0.9-34 had neither family, so nothing here reaches NEWS.
  - `foreignArgsFor` already drops a method's own formals from `predictForeignReasons`, so the three stop refusing the name
    once they take it; the list's `weights` entry, with no method left to refuse it (predict.bartHurdle took `weights`
    with hurdle-completion), goes.
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
20. predict on an ordinal or multinomial fit with a weight fractional, negative or infinite (dec-B594; 0 is taken under
    Open call 2's recommendation): `the posterior predictive 'weights' of <a/an> <family> fit are trial counts, each draw
    the count per category over that many new trials, and must be non-negative whole numbers`. Under Open call 2's option
    (c), "non-negative" reads "positive" and 0 joins the refusal.

Gone: xbart's `xbart does not accept weights of 0 and 1 under family ...`; logistic's refusal of a 0; the logistic refusal of
a two-column response under family = binomial and "logistic"; the multinomial refusal of weights other than all ones; the
bridge's creation refusal of 0/1 weights on probit and ordinal; predict's refusal of `weights` by name on ordinal, nbinom and
multinomial fits ("this family's posterior-predictive draw takes no per-observation weight").

## RNG class

Neutral on the engine for every call taken on the build this branch starts from, main 6f84b5a66, after latents-first, which
moved every augmenting family's draws once and re-recorded equivalence-a6e83f38.rds, bcf-equivalence-a6e83f38.rds and
test-reproducibility-binaryResponse.R (5eb90d2a6), multinomial-equivalence-80b1c8d4.rds staying bitwise (benchmarks/baselines
MANIFEST, their current rows); "today" below means that build and those baselines, never main before it. Shown three ways:

- Bitwise harnesses, run against each engine slice's own `--preclean` library: `benchmarks/R/equivalence.R compare`,
  `bcf-equivalence.R compare` and `multinomial-equivalence.R compare`, every scenario "identical draws (same RNG stream)",
  counted per [Gate hygiene](README.md#gate-hygiene); the four `test-reproducibility-*.R` on the reference build, unchanged
  (two of them carry precision weights, which pins the precision path); tests/cpp unchanged in every existing test. Before
  approximate-draws' S2 lands every Polya-Gamma row here draws by the exact sum, today's sampler. If S2 lands first, its
  default (dec-B590) routes a row above a shape of 64 to the accept-or-keep route, and the harnesses run as that build
  records them; every draw-for-draw gate here against adjacent copies sets the method to the exact sum by name (dec-B593),
  since under the route one row of shape `w n` and its w copies are different chains. No gate here reads the route, so none
  waits on its proof (dec-B591).
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
- A weighted logistic or hazard.logistic fit with a linear or gp leaf or a numeric forest basis (S1's read-offs, reached
  when S5 resolves logistic's NULL to frequency): its leaf covariates' center and scale, a gp leaf's default lengthscale
  and a numeric basis's default sd are read off the expanded rows, so its draws move (Calls made, 21). Both features are
  new in 1.0-0, 0.9-34 having neither, so not in NEWS; S5 re-pins any test that holds such a fit's draws.
- A multinomial count fit's predictive draws at new rows, predict at newdata and its test sample (extract(sample = "test")):
  a one-trial count vector, S x n x K, where they drew one category, S x n (dec-B594 revising dec-B177). The vector is drawn
  by the sequential binomials (Calls made, 18), K - 1 `rbinom` draws per draw and row where `sample.int` took one, so R's
  stream after it moves too; test-multinomial-surface.R's count-matrix arm re-pins. Not in NEWS: multinomial is not in
  0.9-34.

New calls take new draws, each held to an oracle: the expanded rows draw for draw (gaussian, heteroscedastic, nbinom at a
fixed and at a drawn r, multinomial, binomial), the expanded rows in distribution (Student-t only), the exact-posterior arms
of S7; predict's weights on nbinom, ordinal and multinomial against the same R draw written by hand from the stored draws
under one seed (`rnbinom`, and the sequential binomials over the category probabilities). A hurdle fit with 0/1 weights does not move (S5's before/after arm). The
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
18. Every count vector an ordinal or multinomial predictive draw returns (dec-B594) is drawn by the one routine extract's
    count rows use today, `multinomialCountPpdFromProbs` (sequential binomials, one category at a time across every draw
    and row), ordinal's category probabilities included; every category by `multinomialPpdFromProbs`, as today. So on a
    count fit predict(newdata) and predict(newdata, weights = 1) agree draw for draw under one seed, and at the training
    rows of a fit with no zero-trial row predict(weights = rowSums(y)) equals predict() there. Cost: on a fit to categories, predict(weights = 1) and the
    category predict() returns agree in distribution, not draw for draw, and a count fit's draws at new rows move R's
    stream. Alternative: a one-trial row written as the unit vector of its category draw, which keeps that stream and draws
    one vector two ways.
19. A count fit's test sample, extract(sample = "test"), draws as predict at the same rows, a one-trial count vector, so
    main's identity extract(sample = "test") equal to predict(newdata = the test rows) under one seed holds at the new
    shape; man/bart.Rd's sentence that such a fit's test rows, "whose trial counts are unknown, draw one category" (dec-B177),
    follows its predict clause.
    Alternative: test rows keeping one category, which breaks that identity. A test row is a new row in dec-B594's sense,
    having no trial count; if the maintainer reads the ruling as naming predict alone, this is the call to put.
20. A gp leaf's default lengthscale under frequency is the median pairwise distance over the expanded rows' pairs: two
    distinct rows a pair of weight w_a w_b, a row's own copies w_a (w_a - 1) / 2 pairs at distance 0 where that is
    positive, and the 512-point subsample the expansion's own, its stride counted over the total count and each row
    weighted by its overlap with the sampled unit cells of the expanded line. At whole weights the value is the
    expansion's exactly; at unit weights and under precision it is the shipped heuristic, bitwise. A real count extends it
    as a length of that line: a row of count w covers w units and at stride 1 carries w itself. A row split into adjacent
    parts reads the same where every part carries at least 1 (their pairs at distance 0 sum to the same); the self-pair
    term v (v - 1) / 2 is held at 0 below 1, so a split leaving a part below 1 shifts the value, by 0.32 percent for 1 into
    0.5 + 0.5, 0.23 percent for 1.5 into 0.75 + 0.75 and 0.011 percent for 2.5 into 1.25 + 1.25 at stride 2 in the review's
    probe. Where counts below 1 put more than 512 rows into the subsample, a further stride thins them, which whole counts
    never reach. `gpMaxLeafSize`
    keeps counting rows, the one stated exception to the design's Calls made 3: it bounds the leaf's kernel matrix, which
    has one row per data row whatever its count. Alternatives: the heuristic unweighted and recorded as an exception, which
    left sigma 0.057 off the expansion's in the review's probe; and signed self-pair weights, v (v - 1) / 2 negative below
    1, which make every split exact and stay bitwise at whole counts, but drift the median upward where the total count is
    near 2.
21. Logistic's weights are counts (dec-B470), so under its reading every read-off is the expanded rows', the basis scale
    and a linear or gp leaf's standardization included; the basis scale reads the counts with or without a mask, where it
    read them only under one. A weighted logistic fit with a linear leaf or a numeric basis therefore moves once S5
    resolves logistic's NULL to frequency ([RNG class](#rng-class)); both features are new in 1.0-0. Alternative: the
    read-offs keyed on the family, logistic's left unweighted, a fit equal to its expansion in the likelihood and not in
    its scales.
22. The reading lives on the sampler, `SamplerOptions::weightsReading`, the gaussian-type responses holding it from
    construction, rather than in each family: the read-offs are computed by the sampler and the chain, beside the families
    (the leaf standardization, the gp heuristic, the basis scale, the anchor), and one flag read there serves every family.
    The review's logistic finding came from the basis scale reading the counts only under a mask, fixed, not from where the
    flag lives. Alternative: each family holding it and the read-offs asking the family, an interface on every response
    model for a value only the gaussian-type families and logistic read.
23. Where too few rows are in (fewer than two, or one response value over the rows in) and a read-off falls back to every
    row, under frequency it reads every row at its count, as the expansion's every row; a row of weight 0 still counts 0.
    Alternative: the fallback unweighted, as under precision, which no expansion matches.
24. nbinom's center is held at a weight change, as at a mask change, and re-derived over the counts only where the caller
    asks (`updateScale = TRUE`, which reads the weighted center over the rows in). Alternative: a re-derivation at every
    `setWeights`, which moves the transform unasked where gaussian's `setWeights` never does.
25. An nbinom row of count w > 0 cold-starts at PG(w (y + r), 0)'s mean, w (y + r) / 4, as logistic's w / 4, so its working
    response starts where its copies' do; a row of count 0 keeps (y + r) / 4. Alternative: (y + r) / 4 at every row. Under
    latents-first a row in draws omega before the trees read it, so neither moves a draw.
26. The cap at w > 1 lives in the bridge, in refusal 8's words, and holds at creation, `$setWeights`, and every change of
    what it multiplies: `$setResponse` and `$setData` (y), a `$setModel` shape install (r), and `$setCounts` (the row
    totals). It is tested in tinytest at the entries themselves, since the fitting functions refuse these weights until
    S5. Alternative: creation and `$setWeights` alone, which a later y or r change walks past.
27. The bridge reads the weights and the shape's mode in force from the engine, a `weights()` accessor on the sampler
    interface and `SamplerShape::drawsShape`. Alternative: a record on the bridge's holder of whether weights are
    installed, a second copy of engine state.
28. The total-count guard of the anchor, the basis scale and the leaf moments reads "not above 1" to the rounding of the
    total's sum (`countTotalIsAboveOne`): whole counts sum exactly and are unaffected, and weights rescaled to sum to 1,
    which sum to 1 plus 4.4e-16 in the test's probe, take the no-spread values. Alternative: the total strictly above 1,
    which leaves that probe an anchor of 1.2e7.
29. nbinom's and multinomial's engine read their weights as counts under either reading, the likelihood having no other
    one; "precision" on them is refused in R (Refusals, 3, S5). Alternative: the engine refusing them under precision, a
    rule the bridge's resolved reading already states.
30. A numeric basis that varies over rows whose counts total at most 1 is refused at creation, and its scale kept at a
    rescaling `$setForestBasis`, in words naming that total ("the rows its basis is read over carry a total count of at
    most 1"; the warning's "the rows in the likelihood carry a total count of at most 1"), a column constant there keeping
    "constant". Alternative: S1's "constant" for both, which is false of a column that varies.
31. A binomial row's read-offs (a numeric basis's scale, a linear or gp leaf's standardization, the case weight a tree's
    variation is read at) read it at its trials times its weight, 1 absent, under either reading, as the n binary rows it
    stands for; a row of 0 trials is out of the rows in, as a masked row, so a fallback counts it out. Alternative: the
    trials only under frequency, which reads a binomial row as one binary row under precision, a reading no expansion
    matches and S5 refuses on logistic.
32. The engine's names: `SamplerOptions::trials`, `ResponseModel::trials()`, `setBinomialCounts` on the response, chain,
    sampler and facade, and `SamplerShape::holdsTrials`, beside the multinomial `setCounts` rather than through it:
    `supportsCountsMutation` is what `dbarts_sampler_family` reads as multinomial. `Sampler::setData` returns false on a
    binomial response. Alternative: one `setCounts` virtual taking either layout, which mixes an int K-column matrix with
    two double vectors.
33. Under logistic, a data object's counts slot is read for the successes and trials and its y is not read, as
    multinomial does not read it. A slot of K >= 3 columns takes refusal 17; one of another malformed shape "a binomial
    response's counts must be an n x 2 integer matrix of successes and failures, one row per observation". Alternative:
    requiring y to equal the trials at the bridge, which S6's R constructor settles either way.
34. The trials cap at w = 1 is dec-B202's million in `refuseCountOverCap`'s form, "<caller>: family "logistic" requires
    the trials no larger than 1000000; each trial is a Polya-Gamma draw per sweep", at creation and `$setCounts`; for
    w > 1, refusal 8 at creation, `$setWeights` and `$setCounts`. Alternative: no cap at w = 1, which leaves one row's
    draws unbounded.
35. Refusal 13 at the bridge names the first five rows, then ", ...". Alternative: every row, which a large count
    matrix makes unreadable.
36. `$setData` on a binomial sampler is refused at the bridge (Open call 4): "$setData: a binomial sampler's rows and
    their trials are fixed at creation; replace its successes and failures with $setCounts, or make a new sampler".
37. The flat C API's `dbarts_sampler_setResponse` refuses a binomial sampler in refusal 15's words, as `$setResponse`
    does. Alternative: taking a 0/1 vector as the successes against the trials in force, which a row of 0 trials
    makes invalid.
38. A fold over a data handle (`bartcore_createFromHandle`, xbart's) carries the trials with the successes, so S6's
    binomial xbart needs no bridge change. Alternative: refusing a binomial fold until S6.
39. A row of 0 trials goes out of the likelihood in the engine and warns nowhere in S3: the zero-trial warning is the R
    layer's, raised after creation or `$setCounts` returns as multinomial's is, so S6 adds the binomial wording with its
    constructor and setter. Alternative: a warning from the bridge, which would warn at a call the R layer may still
    refuse.

## Open calls

Three, put to the maintainer as written: two left by the ruling of Open call 1, and one raised in S3. Every other user- or
consumer-visible choice the plan makes is settled by dec-B563 to dec-B570, dec-B573, dec-B574, dec-B580, dec-B588,
dec-B594 or the design note's calls made. S5 builds the recommendations; the other options change it by under half an
hour.

### 1. What predict's weights mean on a fit whose one draw is a category (RULED, dec-B594)

Put on 2026-11-05 as options (a) counts only on count fits, (b) 0 and 1 only, (c) counts whenever weights are given,
recommending (a). The maintainer: "What's a correct draw from the posterior predictive distribution here? I always thought
it was well-defined." Shown that it is, Multinomial(w, p(x)) at each posterior draw of the category probabilities, and the
recommendation revised to the rule below: "Use your revised recommendation." Recorded as dec-B594: without weights a draw
has the shape of the fit's response, a category on an ordinal fit and on a multinomial fit to categories, a count vector on
a fit to count rows, one trial at a new row; with weights, the count in each category over that many new trials on all
three; predict(weights = 1) on a fit to categories returns a one-trial count vector. It revises dec-B177's clause that a
count fit's predict at new rows draws one category. Built in S5 as [The R surface (S5, S6)](#the-r-surface-s5-s6) states
it.

### 2. What weights = 0 gives in predict on an ordinal or multinomial fit

On an ordinal or multinomial fit, what should predict(fit, type = "ppd", weights = 0) draw at a row?

Background. Under your ruling of 2026-11-05 predict's weights on these fits count new trials, and a draw with them is the
count in each category; the ruling does not say what 0 trials give. Elsewhere 0 reads three ways. On gaussian-type fits
predict(weights = 0) gives missing draws, a row held out (your ruling of 2026-10-23), and a training row of weight 0 has
missing draws on every family (2026-11-03: "0 weight means out of sample"). Logistic's and probit's predict refuse 0 as not a
positive trial count, and this plan refuses it on the negative binomial the same way. A multinomial count row with no trials
draws all zeros (2026-10-01). 0.9-34 had neither family.

- (a) Missing draws, K missing counts at the row, as a held-out row is on gaussian fits and a 0 training row is everywhere.
  About 5 lines. Cost: these two families read 0 unlike logistic, probit and the negative binomial, which refuse it.
- (b) All-zero counts, the draw of no trials, as a count row of no trials draws. About 3 lines. Cost: a draw that is zero at
  every posterior draw says nothing about the row, and differs from a 0 training row's missing draws.
- (c) Refused, as logistic, probit and the negative binomial refuse it. About 3 lines. Cost: the ruling's reason, that a
  well-defined draw is not refused, and a weight vector copied from an ordinal fit's own 0/1 weights fails on its 0s.

Recommendation: (a). A 0 means out of sample, your words of 2026-11-03, and under (a) predict at the training rows with the
fit's own weights gives missing draws at its 0 rows, as predict without weights does.

### 3. What a test sample's weights draw on an ordinal or multinomial fit

When an ordinal or multinomial fit is given test rows whose data carry the weights column, what do those weights draw?

Background. A fit given test rows through the formula reads the test rows' weights from the same column as the training
weights: 0 and 1 marking rows in and out on an ordinal fit, counts on a multinomial one. On logistic and the binomial the
test rows' weights are their trials, as predict's are, and probit's test weights are 0 and 1, where one trial and a 0/1
draw are the same thing. Under your ruling predict's weights on these two families are trials and give a count per
category, so read the same way an ordinal fit's 0/1 membership column would turn its test draws into one-trial count
vectors whenever the test data carry it. 0.9-34 had neither family.

- (a) Each family as its sibling: a multinomial fit's test weights are trials, as logistic's are, a count per category over
  that many trials; an ordinal fit's are 0 and 1, as probit's are, 1 drawing a category and 0 missing. About 15 lines.
  Cost: on an ordinal fit, a test row of weight 1 draws a category where predict(weights = 1) at the same row draws a
  one-trial count vector.
- (b) As predict's weights on both families: a count per category wherever test weights are given. About 10 lines. Cost:
  an ordinal fit's test draws change shape with whether its test data happen to carry the membership column.
- (c) The test sample reads only a 0, as missing draws, and otherwise draws as predict does without weights. About 10
  lines. Cost: a multinomial fit's test rows cannot draw counts at their trials, which logistic's can.

Recommendation: (a). An ordinal fit's weights mark rows in and out wherever they come from, so its draws stay categories,
and a multinomial fit's test rows take their trials as logistic's do.

### 4. What `$setData` does on a binomial sampler (raised in S3)

When a sampler fitted to successes and failures is handed a new data object, what happens?

Background. The design note's call 9 swaps a binomial response through `$setCounts` and refuses a flat vector at
`$setResponse`; it says nothing of `$setData`, which on every other family replaces the rows, the response and the
weights together and may change the number of rows. A multinomial sampler refuses `$setData` ("make a new sampler
instead"). S3's bridge refuses it on a binomial sampler, the engine's trials being sized to the rows it was built on;
nothing reaches it until S6 lets the fitting functions take a binomial response. 0.9-34 had no binomial response.

- (a) Refused, as on a multinomial sampler: `$setData: a binomial sampler's rows and their trials are fixed at creation;
  replace its successes and failures with $setCounts, or make a new sampler`. No further work. Cost: a binomial sampler
  inside a larger Gibbs sampler cannot swap its predictors and counts in one call, as a binary one can.
- (b) Taken: the new object's counts slot gives the successes and trials, n free, as its y does on a binary sampler.
  About 40 lines in the engine and bridge, and S6's R method. Cost: a fourth route by which the trials change, each with
  the caps and the zero-trial mask.

Recommendation: (a). It is multinomial's rule for the same reason, the count response is not a column of the data,
and `$setCounts` plus a new sampler cover the cases (b) would; (b) can follow if a consumer asks.

## Slices

Each slice is a commit series on frequency-weights, its own review and gates; the feature merges once after S8. S0 landed
apart (above). The branch is rebased onto main 6f84b5a66, after latents-first merged, and every gate below runs against that
build or a later main it rebases onto before S1 starts.

| slice | tier | content | agent-hours (2026-11-05) |
|---|---|---|---|
| S1 | changes draws | gaussian-type families, read-offs, log-likelihood channel; no digest | ~12.5 (12.5) |
| S2 | changes draws | nbinom and multinomial shapes times w, caps, setters storing only | ~7 (7) |
| S3 | changes draws | the binomial row, `ext_lchoose` with choose.c in the standalone copy, `$setCounts` storing only | ~9.5 (9.5) |
| S4 | changes draws | membership on the data object, logistic's 0; no reactivation redraw | ~6.5 (6.5) |
| S5 | R surface | weights.type, read-offs in R, predictive draws, predict's weights on ordinal, nbinom, multinomial (dec-B594's count vectors), the hurdle row | ~15.5 (14) |
| S6 | R surface | the binomial surface and xbart | ~9 (9) |
| S7 | gate-only | frequency-exact.R, the SBC arm | ~6 (6) |
| S8 | records | records | ~1 (1) |
| reviews, gates, rebases | | four engine reviews ~18, three R and gate reviews ~8, gates ~5, rebases ~3 | ~34 (33.5) |

About 101 agent-hours (99 on 2026-11-05, 94 when first planned; S0's 3 spent on main). The 2026-11-05 moves: S1 and S4 lost
the digest and the reactivation redraw that latents-first retired (about half an hour each); S3 gained the seventh math name
(about an hour, less half an hour for a `$setCounts` that only stores); S5 gained predict's weights on three methods and the
hurdle's row (about 4 hours, with review half an hour more); rebases gained an hour for latents-first's re-recorded
baselines. The 2026-11-06 move: dec-B594's rule against the old recommendation, count vectors with weights on all three
category fits, the count fit's one-trial vector at new rows and in its test sample, Open calls 2's and 3's recommendations,
and the re-pins in test-multinomial-surface.R and test-predict-training-rows.R (about 50 lines of R, 10 of help and 80 of
tests: S5 about 1.5 hours, its review half an hour).

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
`dnorm` and `dt` times w; the doors against the expansion's same door, draw for draw (`updateStandardization` after a mask
change, the heteroscedastic re-anchor by a rescaled response, `updateBasisScale`); a Student-t `setWeights` to a new pointer
against the sampler created with those weights; a row out, at weight 0 or masked, inert on gaussian and Student-t; the gp
lengthscale against the expansion's on both sides of the subsample; logistic's basis scale against the expansion's with and
without a mask; the fallback against the expansion's every row; the anchor at a total count of exactly 1 and at whole
weights rescaled to sum to 1 the no-spread anchor, the basis scale and the leaf moments there taking no scale; a tinytest Student-t fit at whole weights against its
expansion's df; the bitwise trio; the full tinytest suite; sanitizers on tests/cpp and on the R-loaded path;
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
consuming no variate; a state round trip; the trio; the suite; sanitizers. As built: the nbinom fits within 2.5e-14 of the
copies' at a fixed r and 3.1e-14 at a drawn r over the run, its log-likelihood against the sum over the copies; the
multinomial fits within 6.3e-15 of the copies' and bitwise against the counts and trials times w; a 0 weight on each
family against the row masked, bitwise at the rows in and within 1e-12 at the row out (its working response rides the
trees); `$setWeights` against the sampler created with those weights, consuming no variate; a state round trip within
1e-9; the caps in tinytest at the bridge's entries (test-frequency-counts-bridge.R: creation, `$setWeights`,
`$setResponse` and `$setCounts`, a drawn r capped at its grid's largest point, a fixed r at itself); the anchor's guard at
a total of exactly 1 and at rescaled weights (S1's follow-up, amendments row 7); `mutation-battery.R verify-anchors`.
The review's fix round added the cap's refusal at `$setData`, at a `$setModel` shape install (a drawn r and a fixed r of
600,000 refused, a fixed r of 10 taken) and at a multinomial `$setWeights` over a row total already installed; the
multinomial expansion and counts-times-w test repeated under a category offset (within 3.8e-15, and bitwise); and the
basis refusal and kept-scale warning at a total count of at most 1 by their words (Calls made, 30).

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
and draws none. As built: the binomial row within 5.5e-14 of its split rows over 300 sweeps, the same variates consumed;
at frequency weight w bitwise the row of k w successes in n w trials, and a one-trial row bitwise the binary response,
weighted and not, its log-likelihood included; the log-likelihood against `dbinom` times w over 252 rows, worst 1.5e-11
relative, at n = 1e6 2.3e-12 at k = n / 2 and 1.3e-13 at k = n / 1000, and after a successes swap the constants recomputed;
a zero-trial row a masked row of one trial, bitwise at the rows in, reporting a log-likelihood of 0, and a zero-weight
binomial row finite; the cold start by its expressions; a state, copy and reload round trip within 1e-9; `$setCounts`
consuming no variate, the latents unchanged and the working response the one rebuilt at every row; the read-offs (a basis
scale under both readings and its rescale, a linear leaf's moments, the fallback, after a count or weight change) against
the n binary rows within 1e-14 (Calls made, 31); `ext_lchoose` against Rmath bitwise over 4,115 (n, k) pairs, choose.c
needing nothing beyond the copy, `standalone-math` passing with no R symbol and `r-math-symbols` naming `Rf_lchoose`; the
bridge in tinytest (test-binomial-bridge.R, 20 results: creation, a fold, `$setCounts`, the caps at w = 1 and w > 1, the
refusals); the trio bitwise (55, 19 and 11 scenarios); the suite, 41,227 results; tests/cpp and the binomial and
logistic tinytest files under ASan and UBSan, clean. Mutant 22 is killed by the `standalone-math` driver (2 failures),
not by the nmath suite, which calls the copy's routine by its own name rather than through the macro.

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

### S5: R surface, the reading (~485 R, ~145 help, ~680 tinytest)

Tier "R surface", pure R. Moves the weighted logistic training draw, by dec-B580 a 0 row's draws on probit, ordinal, hazard,
hazard.probit and nbinom, and by dec-B594 a multinomial count fit's draws at new rows and in its test sample; re-pinned.
Steps: weights.type on bart, dbarts, dbartsSpec and xbart, its validation and the resolution table, the hurdle's row and its
routing to both parts; Refusals 1 to 10, 19 and 20; the starting sigma on both routes and xbart's per fold, the chi anchor;
`sampleFromPPD` (one observation, NA at a 0 on every family, dec-B580, the fitted value kept), predict's and weights.test's
doors, predict's `weights` on predict.bartOrdinal, predict.bartNegbin and predict.bartMultinomial and their paragraphs in
man/bart.Rd; dec-B594's count vectors through `multinomialCountPpdFromProbs` (Calls made, 18), the count fit's one-trial
vector at new rows and in its test sample (Calls made, 19), weights = 0 and the test sample's weights as Open calls 2 and 3
rule; the weights entry of `predictForeignReasons` removed; `pointwiseLogLikelihood`; help and NEWS.
Gates: the touched test files and the new test-frequency-weights.R (each refusal by its words, counted warnings; the
resolution table row by row, hurdle.lognormal's included; sigma's starting value and the chi anchor against the expanded data
on the dense and the sparse routes and in xbart's folds; predict's and weights.test's doors; NA draws at a 0 on every family
with the stream unmoved; nbinom predict(weights = 3) equal under one seed to `rnbinom(size = 3 r, mu = 3 mu)` from the stored
draws; `pointwiseLogLikelihood` against the engine's channel). dec-B594's cases, each on an ordinal fit, a multinomial fit to
categories and one to count rows, at training rows and at newdata, with combineChains TRUE and FALSE:
  - without weights: a category where the fit models categories, as today (the existing draws unchanged under one seed),
    a count vector at each row's trials at a count fit's training rows, as today, and at its new rows a one-trial count
    vector, S x n x K, each draw a unit vector;
  - with weights w (a vector with 1, 2 and 7 among its values): the S x n x K counts equal under one seed to the sequential
    binomials written by hand from the stored or replayed category probabilities at those trials, each row's counts
    summing to its w; on a count fit at the training rows, w in place of the row's trials;
  - weights = 1 on a fit to categories: a one-trial count vector, the unit vector of a category, whose category frequencies
    over 4,000 draws match predict()'s within Monte Carlo error; at a 0 training row it draws where predict() gives NA;
  - on a count fit, predict(newdata) and predict(newdata, weights = 1) identical under one seed, and extract(sample = "test")
    identical to predict at the test rows (test-multinomial-surface.R's identity, re-pinned at the new shape);
  - weights = 0 as Open call 2 rules (recommended: NA counts at the row, the stream unmoved), the test sample's weights as
    Open call 3 rules (recommended: an ordinal test row of weight 1 a category and of 0 NA; a multinomial test row's counts
    over its weight in trials);
  - type = "ev" and "class" unchanged by weights, reading none; Refusal 20 at 0.5, -1 and Inf; an NA weight under
    na.keepPredictors, na.pass and na.omit as at newdata (dec-B482); test-predict-training-rows.R's two refusals of
    `weights` by name on nbinom and ordinal replaced by these arms.
Also test-hurdle.R's before/after arm (a hurdle fit with 0/1 weights at weights.type NULL and "frequency"
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
doc-freshness. Per R slice: the touched test files, the snapshots, the static checks; CI on the branch tip runs the rest.
S5's predict arms (dec-B594's cases, listed under S5) are by-hand oracles under one seed, not snapshots, so they hold on
every build; no snapshot file holds an ordinal or multinomial predictive draw. At
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
10. S2: the product cap at w = 1 - the y = 999,999 test (tinytest, test-frequency-counts-bridge.R: the cap is the
    bridge's).
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
22. S3: `ext_lchoose` mapped to `std::lgamma` arithmetic under `DBARTS_HOST_STANDALONE` - the `standalone-math` driver's
    lchoose checks (the nmath suite calls the copy by its own name, past the macro).
23. S3: `$setCounts` keeps the old working response - the tests/cpp working-response read after the swap.
24. S5: nbinom predict(weights = w) draws one count at mean mu - the `rnbinom(size = 3 r, mu = 3 mu)` test.
25. S5: the hurdle hands NULL to its zero part and "frequency" to its positive part with no check of "precision" - the
    resolution table's hurdle row ("precision" refused in the hurdle's words).
26. S5: predict(weights = 1) on a fit to categories returns the category, the weights read as 0/1 membership - the
    one-trial count vector arm (shape S x n x K, unit vectors).
27. S5: a count fit's predict at newdata without weights keeps dec-B177's one category - the predict(newdata) against
    predict(newdata, weights = 1) identity.
28. S5: at a count fit's training rows predict(weights = w) draws at the row's own trials - the by-hand sequential binomials
    at w, each row's counts summing to w.
29. S5: ordinal's count vectors drawn from the cumulative probabilities rather than the category probabilities - the
    ordinal by-hand arm.

## Consumers

- stan4bart (main 27ef163): passes its weights to the BART part at weights.type NULL, a precision on its continuous fits as in
  its own Stan target; its probit rule passes 0/1 weights, which now ride `data.bart` into every sampler it builds, and its
  `active.bart` from `bart_spec$active` becomes NULL, so its `setActiveRows` installs are skipped and `restoreBartSampler`
  rebuilds from the data, which carries the weights (`data.bart.light` drops only x); the critique confirmed its draws do not
  move. No call changes; an optional cleanup of about 30 lines drops the dead `active` plumbing, in the consumer sweep, not
  this feature. Its `install_bart_range` reads `sampler$data@weights` on continuous fits only, unchanged.
- bartCause (main 11e3030, read at 1e4b882; only its regression test changed between): passes population weights to bart's gaussian response fits at weights.type NULL, unchanged; drops
  weights from its propensity fits; its xbart call (bayesOpt) uses the default loss, which on a binary treatment is log and does
  not move under S0. No call changes.
- treatSens (268b2db): passes no weights to dbarts; no change.
- bairrtt: reads none.
- latents-first's slice 4 found stan4bart, bairrtt and treatSens unchanged and restated bartCause's three propensity-score
  pins (11e3030); nothing here moves them again. None of the four calls predict on an ordinal, nbinom or multinomial fit or
  fits a hurdle (grep of their R trees on 2026-11-06 for those family tokens and hurdle: only stan4bart's Stan prior help
  and a comment on the probit and ordinal mask), so dec-B594 reaches none of them.

## Interaction with the other pending plans, and order

- forest-door landed (a6a16c31a), with setmodel-parameters: its forest records carry shares, `$setModel` gained refusals in
  R/fittingLevel.R beside which Refusal 10 sits, and chain.hpp, combiner.hpp, model.hpp and R_interface_bartcore.cpp moved; this
  plan's symbols were re-read on 5a99172e6 and none it names changed role.
- Build order, the orchestrator's call as latents-first's plan states it: coefficient-prior-live, latents-first, this
  feature, host-neutral-engine S2a onward. Engine features merge serially; the first two have merged and this branch is
  rebased onto main 6f84b5a66 after both.
- coefficient-prior-live (landed on main 2026-11-05; read at 05a12bb79 against main): touches no file's weight code and
  none of this plan's engine symbols. It does not touch model.hpp; its combiner.hpp and chain.hpp hunks are the live coefficient
  prior (`setForestAmplitudePrior`, `setForestCoefficientPrior`), not `basisColumnScale`, `updateBasisScale` or the
  Polya-Gamma loop. Two textual neighbors: R/fittingLevel.R's `setModelOnForests` gains `refuseCoefficientDifference`, beside
  which Refusal 10 sits; R/spec.R's `resolveSamplerSpec` gains `leafPriorGivesNoSize`, near where the spec resolves
  weights.type. A forest released to a drawn coefficient takes creation's sd, which under frequency is the weighted
  `basisColumnScale` S1 builds, so the two agree with no edit.
- latents-first (landed on main 2026-11-05, its calls in dec-A263; docs/plans/latents-first.md, Order and overlaps): every
  sweep draws its augmentation variables first; no setter draws them; the weights and survival digests, the Student-t
  zero-weight record and setState's reconciliation retired. What that does here, each rewritten above: design call 5 and the engine bullet of
  [The reading and where it lives](#the-reading-and-where-it-lives) (a state carries no reading; Calls made, 13); mutant 5 and
  the test-state-weight-pairing.R arm (replaced by the mask-change count test); `$setCounts` on a binomial sampler stores
  (S3); a reactivated row draws nothing at the setter (S4); the Student-t install list drops the state install (S1); the
  baselines this plan's neutrality is measured against are latents-first's (RNG class). It also shares `ProbitResponse`,
  `OrdinalResponse`, `LogisticResponse` and `NBResponse`'s setters and the setter help, which latents-first rewrote.
- hurdle-completion (landed on main 2026-11-05, R only, its calls in dec-A262; its Interaction with frequency-weights, and
  order): the hurdle takes 0/1 weights on both parts (dec-B588). Here: its resolution-table row is probit's, Refusal 4
  refuses "precision" in its words, the arm hands the resolved reading to both parts (S5, Calls made, 15); textual overlaps
  in R/bart.R, R/generics.R, man/bart.Rd, inst/NEWS.Rd and test-zero-weights.R.
- host-neutral-engine: its S1 (build-time math and printing) landed on main (dec-A260). Its copy of R's math holds six
  routines and no `lchoose`, so S3 here adds the seventh name and copies choose.c (Calls made, 16). Its main feature rewrites
  every creation route's SEXP reads into a plain-C spec (S2b) that carries `weightsType` per dec-B564 and, under this plan,
  `trials` beside y under LOGISTIC. Under the build order this feature lands first: host-neutral's S2b translates the weights
  reading, the counts slot under logistic and the membership rule as code that exists, and its S3 adds the header field and
  the getLatents sentence on a frequency row's latent. Cost: host-neutral's S2a rebases over this feature's bridge changes
  (~200 lines in the creation and mutation conduits it rewrites; latents-first removed about 140 bridge lines there first),
  and S2a waits for this feature to merge. Were the order reversed, this feature would rebase its bridge reads onto S2b's
  spec form and add `trials` to the header itself, about the same lines in the other direction.
- approximate-draws (planned, docs/plans/approximate-draws.md on its branch at e0c62e566; its three open calls ruled
  2026-11-05): until the defaults study a logistic, binomial, logistic-hazard or negative binomial row draws by the
  accept-or-keep route where its shape is above 64 or not a whole number and by the exact sum otherwise, the route taking
  that place "only once its exactness is proven (dec-B591) and its gates pass" (dec-B590; the maintainer: "And yes, you can
  use your recommendation."); the route is proven by independent agents before it is built (dec-B591); every sampler of
  that arc ships in 1.0-0, the approximate ones opt-in, and the default comes from the defaults study, run before 1.0-0
  right after the build, so dec-B590's rule holds only for the development builds between (dec-B592); the method is chosen
  by name, NULL for the automatic rule, an approximate one's accuracy by a separate tolerance (dec-B593). The binomial at
  large n is its cost case. How this plan's shapes times w behave, the rulings changing only which sampler draws a row:
  - Before it, as this plan builds them: a row's shape (`w` on logistic, `w n` on a binomial row, `w (y + r)` on nbinom,
    `w n_i` per multinomial category) is a whole number drawn as that many unit Polya-Gamma draws, about 160 ns each (164 ns
    at psi 0, 159 ns at psi 3, p12cost.R on the shipped build), so 1,000 rows of 1,000 trials cost about 0.16 s a sweep and a
    1,000-sweep chain about 160 s, linear in the total; fractional weights refused on every Polya-Gamma family (Refusals 5 to
    7).
  - After its S2 (the exact accept-or-keep route at the sweep's start for logistic, binomial, hazard.logistic and nbinom),
    once the route is proven and gated: each row's shape is unchanged and the same rows are taken; under dec-B590's rule a
    row whose shape here is above 64 (`w` on logistic, `w n` on a binomial row, `w (y + r)` on nbinom) takes the route,
    exact, at about 3 to 9 us a row whatever the shape, and the rest the sum as here, until the defaults study's default
    replaces that rule before 1.0-0 (dec-B592); the method setting forces either (dec-B593). A row a setter moved (a weight,
    count or mask change here among them) is drawn by the sum at the next sweep. Multinomial's `w n_i` keeps the sum under
    every exact method; its approximate series is opt-in (dec-B592). Draws shift only where the route is selected.
  - After its S3: real weights on logistic, binomial and hazard.logistic are taken, lifting the fractional clause of
    Refusal 5 there; a fractional shape `w` or `w n` is drawn by the route under dec-B590's rule (its Real shapes), while
    before S3 this plan refuses such weights, so no fractional shape reaches the sum; under an exact method a between-sweep change that would mark a
    fractional-shape row is refused by that plan's name. nbinom's and multinomial's weights stay whole (Refusals 6 and 7).
  - Nothing in this plan depends on the route or its proof, which is still under independent check (dec-B591): every
    slice here draws by the exact sum, today's sampler, and is gated there.
  - Order: its S3 must land after this feature, since it reads weights.type, the binomial row and this plan's Refusal 5. Its
    S1 and S2 may land on either side: every new shape here goes through the family's existing loop (Constraints), which S2
    meets once per family, and if S2 lands first every draw-for-draw gate here runs with the method set to the exact sum
    (RNG class). Whichever lands second rebases over the other's `LogisticResponse` and `NBResponse` hunks.
- weighted-binary (parked, after 1.0-0): its integer-weight probit is dec-B563's truncated-normal case, not scheduled; its
  real-weight logistic is the real-shape draw above.
- Nothing else touching src/bartcore/model.hpp's response families merges while S1 to S4 are in review.

## Cost

S0, landed: about 60 lines and 3 agent-hours. The feature about 4,450 lines, 1,610 of them package code; about 101 agent-hours
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
- stan4bart's, bartCause's and treatSens's call sites were read on main 27ef163, 1e4b882 and master 268b2db; bartCause's main
  has moved only its regression test since (11e3030); their 1.0 compat branches may differ.
- latents-first landed; checked on main 6f84b5a66 (2026-11-06) that no `weightsDigest`, `reapplyWeights` or `redrawEntering`
  remains in src or R, and the re-recorded baselines read from the MANIFEST's current rows. Its calls as built (dec-A263) were
  read for this plan's clauses; none changes one.
- hurdle-completion landed (dec-A262); its refusal of other weights (`refuseHurdleWeights`) and predict.bartHurdle's
  `weights` were read on main. The hurdle positive part's draws at 0/1 under "frequency" equalling precision's is argued
  from the count and constant leaves, and S5's before/after arm shows it.
- approximate-draws' open calls are ruled (dec-B590, dec-B592, dec-B593); its plan on its branch (e0c62e566) predates the
  rulings, and its route's proof is under independent check (dec-B591). The before and after above follow the rulings.
- dec-B594's cases were read from the code (`multinomialPpdFromProbs`, `multinomialCountPpdFromProbs`, predict.bartOrdinal
  drawing its categories through the first) and not run; the line counts of the rule and of Open calls 2 and 3 are
  estimates, no option prototyped.
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

## Amendments of 2026-11-06

For dec-B590 to dec-B594 and the rebase onto main 6f84b5a66; read, not run.

| # | change | where it went |
|---|---|---|
| 1 | dec-B594 rules Open call 1: predict's weights on ordinal and multinomial fits are trials, a draw with them the count per category; without them the response's shape, a count fit's new rows a one-trial count vector; dec-B177's predict clause revised | Open calls, 1 (ruled, quoting the maintainer); front block; Summary; Goal; The R surface (S5), with a table of what predict returns; Resolution table; Refusals, 20 rewritten; RNG class (a count fit's new-row draws move); S5 steps and gates; mutants 26 to 29; Consumers; the design note's Calls made 10; TODO; INDEX |
| 2 | Calls made 18 (one count-vector routine, `multinomialCountPpdFromProbs`, for every count draw) and 19 (a count fit's test sample draws as predict at those rows) | Calls made; The R surface (S5); RNG class |
| 3 | Two calls the ruling leaves, put as Open calls 2 (weights = 0: recommended NA) and 3 (a test sample's weights: recommended ordinal 0/1, multinomial trials) | Open calls, 2 and 3; Resolution table; Summary |
| 4 | dec-B590 to dec-B593: before approximate-draws lands every Polya-Gamma row here draws by the exact sum; after its S2, rows above 64 take the route under dec-B590 until the defaults study's default (dec-B592), the gates here setting the method to the sum by name (dec-B593); after its S3 a fractional shape takes the route; nothing here waits on the proof (dec-B591) | Interaction, approximate-draws; RNG class; Constraints; Summary; Not verified |
| 5 | Re-read against main 6f84b5a66: the retired `weightsDigest` link replaced by `Chain::drawAugmentation`; man/predict.bart.Rd, which does not exist, replaced by man/bartBT.Rd and man/bart.Rd's paragraphs; the baselines named from the MANIFEST (equivalence and bcf-equivalence a6e83f38, multinomial 80b1c8d4, not re-recorded); test-state-zero-weight-rows.R deleted; latents-first, hurdle-completion and coefficient-prior-live as landed; bartCause at 11e3030; the design note's link to the retired `TResponse::reapplyWeights` dropped | Context; front block; RNG class; Slices; Interaction; Consumers; Not verified; the design note's Calls made 5 |
| 6 | Slices re-estimated: S5 ~15.5 (14), R reviews ~8 (7.5); ~101 agent-hours (150 to 200), ~4,450 lines | front block; Slices; Cost; Summary |

## Amendments of 2026-11-06, the S1 review

For the S1 review's FIX FIRST and the calls on its two open items.

| # | change | where it went |
|---|---|---|
| 1 | A gp leaf's default lengthscale read off the expanded rows' pairs; `gpMaxLeafSize` counting rows, the one stated exception | Calls made, 20; Engine (S1) |
| 2 | Logistic's read-offs read its counts, the basis scale with or without a mask; a weighted logistic fit with a linear leaf or a numeric basis moves once S5 resolves it | Calls made, 21; front block; Goal; Constraints; Engine (S2); RNG class |
| 3 | The reading described where it is built, on the sampler | The reading and where it lives; Calls made, 22 |
| 4 | A read-off's fallback to every row reads each row at its count | Calls made, 23; Engine (S1) |
| 5 | Tests for the review's surviving mutants and for rows 1, 2 and 4 | S1 gates |
| 6 | Calls made 20's split sentence narrowed to parts of at least 1, the shift below it stated, signed self-pair weights recorded as the alternative; under frequency the anchor's sd over rows read holding one value is 0 exactly, as the expansion's | Calls made, 20; the anchor rule in Engine (S1) |
| 7 | The anchor's sd, the basis scale and the leaf moments take no scale where the counts over the rows read total no more than 1 to their sum's rounding, so weights rescaled to sum to 1 take the no-spread values | Engine (S1); S1 gates |
