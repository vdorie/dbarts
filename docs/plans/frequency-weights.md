# frequency-weights: weights read as counts of identical rows beside the precision default, and the binomial response

Status: PLANNED 2026-11-03 on main (3e4c002d8), from docs/design/frequency-weights.md. Not built. Rulings, quoted where they
settle a call: dec-B125 (frequency weights join precision weights by explicit choice; precision stays the default), dec-A11
(sigma's count under precision), dec-B563 (the gaussian-type and Polya-Gamma families take frequency weights; "Let's do (b)
plus the Polya-Gamma families."), dec-B564 (weights.type, "I like (a)."), dec-B565 (a frequency row's predictive draw is one
observation), dec-B566 (a count family's training draw is one observation; dec-B470's replace clause revised), dec-B567 (the
binomial response; "Yes, have the binomial added. It can look like similar R interfaces."), dec-B568 (real weights where the
draw is exact), dec-B569 (membership as 0/1 weights on the data object), dec-B570 (a zero weight alike on every family),
dec-B573 (a proportion with the trials as weights refused; "Use (a)."), dec-B574 (probabilities for expected values, counts
for draws; "Use (a)."), with dec-B133 and dec-B177 (zero-trial rows), dec-B202 (the one-million cap), dec-B470 and dec-B472
(predict's weights), dec-B494 (a feature merges once), dec-B85 (the header change frequency weights' field rides). The runs
behind it: the note's probes (scratch/fw-note, the critique's p1 to p11 under scratch/fw-crit/probes), and p12cost.R in
scratch/fw-note (a unit Polya-Gamma draw at 164 ns, psi 0, and 159 ns, psi 3, on the shipped build of c3d268c8f on the
shared laptop). No prototype was built here.

agent: S1, S2, S3 and S4 one opus implementer each (effort high), S1 and S3 the riskiest; S5 and S6 one sonnet implementer
each (effort medium); S7 sonnet; S8 sonnet. One opus review per engine slice (S1 to S4) that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior; one sonnet review per R slice (S5, S6) and for the
gate slice S7 ([Process by risk](README.md#process-by-risk)). Tier "Changes draws" for S1 to S4: new draws on new calls and the
creation and mask paths of every family rewritten, though no call taken today is meant to move; tier "R surface" for S5 and
S6, each of which moves one user-visible value class named under [RNG class](#rng-class).
rng: neutral, bit for bit, for every engine draw of a call taken today: the equivalence trio, the four snapshot files,
tests/cpp and the exact gates do not move. Two R-level values move by ruling and are re-pinned: a weighted logistic or
hazard.logistic fit's training predictive draws (dec-B566: Bernoulli(p), not Binomial(w, p)) and a zero weight's pointwise
log-likelihood on logistic (dec-B570: NaN, not 0). New calls (weights.type = "frequency" on the gaussian-type families,
weights on nbinom and multinomial, a binomial response, a 0 on logistic) take new draws, each held to the expanded-data fit.
See [RNG class](#rng-class).
window: for 1.0-0 (dec-B567: "before 1.0-0"). One branch, frequency-weights, merged once (dec-B494). After
coefficient-prior-live lands; ahead of host-neutral-engine's main feature (S2a), recommended under
[Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order).
budget: ~3,500 lines in the tree, ~1,330 of them package code: engine ~430 (S1 ~170, S2 ~90, S3 ~110, S4 ~60), bridge ~170, R
~730 (S4 ~80, S5 ~380, S6 ~270); help ~150, NEWS ~25; tests/cpp ~600, tinytest ~1,000, gate scripts ~350, records ~60. Agent
time, estimated and not measured: S1 ~9, S2 ~6, S3 ~8, S4 ~6, S5 ~12, S6 ~10, S7 ~5, S8 ~1, four engine reviews with mutants and
fix rounds ~16, two R reviews and S7's ~7, gates and sanitizers ~5, rebases ~2: ~87. Plans here run 1.5 to 2 times low, so 130
to 175.

## Summary for the maintainer

You ruled that a weight can also mean a count: a row of weight 3 stands for three identical observations, and the fit is the
fit to the data with every row written out that many times. Precision stays the default, chosen with one new argument beside
the weights. This plan builds that, the membership and zero-weight rulings that came with it, and the binomial response.

What a user sees. bart, dbarts, dbartsSpec and xbart take weights.type. Left unset it means what weights mean today: a
precision on gaussian, Student-t and heteroscedastic fits, a count on logistic, negative binomial and multinomial fits.
"frequency" on the first three changes one thing, how many observations the noise level is estimated from, and takes weights
of any positive size. Negative binomial and multinomial fits take whole-number weights for the first time. Probit and
ordinal keep their 0/1 weights for membership, which now stay on the sampler's data object. A weight of 0 takes a row out of
the likelihood on every family, logistic included, which refused it. A row's predictive draw is one new observation like the
row; on a weighted logistic fit that replaces today's count of successes among w copies, which the help and NEWS state.

The binomial response is written as R writes it, cbind(successes, failures), under family = binomial or "logistic", and
xbart takes it too. Its fitted values are probabilities and its predictive draws counts, as glm's are. A proportion with the
trial counts as weights is refused with a message showing the cbind form. xbart's losses on a binary or binomial fit compare
the observed proportion with the posterior mean probability; that also mends the binary rmse on main, which compared the 0/1
response with the latent scale.

What it costs. Nothing per sweep on the gaussian-type families. On the count families a sweep grows with the total count, as
logistic's does now: each trial of each row is about 160 nanoseconds a sweep, so 1,000 rows of 1,000 trials take about 0.16
seconds a sweep. The approximate-draw research you opened is where that changes; this plan keeps every draw exact.

What does not move. No posterior draw of any fit that runs today changes: the engine's random numbers, the four recorded
snapshot files and the equivalence baselines stay bit for bit, and each slice proves it. Three reported values move by your
rulings and are stated in NEWS: a weighted logistic fit's training predictive draws, a logistic zero weight's
log-likelihood (refused today) and a binary xbart's rmse. stan4bart, bartCause and treatSens need no
change; stan4bart can drop about 30 lines that become dead.

There are no open calls: the rulings settle every choice a user or a consumer sees. About 3,500 lines over eight slices,
about 87 agent-hours as estimated, 130 to 175 by the record of plans here. It should land after coefficient-prior-live and
before the host-neutral engine's main work, so that work translates the weights reading and the trial counts as they exist
rather than reserving room for them.

## Goal

After it lands: a frequency-weighted fit on the gaussian-type families equals the fit to the expanded rows draw for draw; the
Polya-Gamma families take whole-number counts with shape times w; a binomial row is one logistic row of shape n w drawn as the
two-row expansion is; probit and ordinal membership rides the data object; a 0 means the same on every family; the R surface
of predict, fitted, residuals, extract, the log-likelihood and xbart reads both readings and the binomial as the rulings
state; and no call taken today moves an engine draw.

## Context

- Design: docs/design/frequency-weights.md (its Calls made 1 to 14 are the design calls this plan builds; this plan's own are
  below), its section The two readings differ only in the scale's count, and its table What moves, place by place.
- Engine: [`GaussianResponse`](../../src/bartcore/model.hpp) (the positive-row count `numPositiveWeights_`, recounted by
  `installWeights`), [`ChiSquaredScalePrior::drawSigmaSqFromPosterior`](../../src/bartcore/model.hpp),
  [`ConstantVarianceLeaf`](../../src/bartcore/model.hpp) (its `accumulate` counts positive rows),
  [`TResponse`](../../src/bartcore/model.hpp) and [`ResidualDfPrior::drawIndex`](../../src/bartcore/model.hpp),
  [`LogisticResponse`](../../src/bartcore/model.hpp) (shape from `lround(weights_[i])`, working response
  `w (y - 1/2) / omega - offset`, `computeLogLikelihood`), [`NBResponse`](../../src/bartcore/model.hpp) and its
  `NBShapePrior`, the multinomial Polya-Gamma loop in [`MultinomialForestCombiner`](../../src/bartcore/combiner.hpp),
  [`ProbitResponse`](../../src/bartcore/model.hpp) and [`OrdinalResponse`](../../src/bartcore/model.hpp) (the active-row mask
  only), the weights digest ([`SamplerBase::weightsDigest`](../../src/bartcore/facade.hpp)).
- Bridge: [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp) and
  [`refuseBinaryWeightChange`](../../src/R_interface_bartcore.cpp), [`refuseCountOverCap`](../../src/R_interface_bartcore.cpp),
  the creation refusals of ordinal, nbinom and aft weights beside them, the weights slot read in the data conversion.
- R: [`enforceWeightPolicy`](../../R/spec.R), [`isMaskedWeightFamily`](../../R/spec.R),
  [`refuseNonCountWeights`](../../R/spec.R), [`refuseNonBinaryResponse`](../../R/spec.R), [`dbartsSpec`](../../R/spec.R) and its
  `active` element, [`restateMatrixResponseError`](../../R/data.R) and the classed `dbartsMatrixResponseError`,
  [`sampleFromPPD`](../../R/generics.R), [`pointwiseLogLikelihood`](../../R/generics.R), [`predict.bart`](../../R/generics.R),
  [`fitted.bart`](../../R/generics.R), [`residuals.bart`](../../R/generics.R), [`extract.bart`](../../R/generics.R),
  [`xbartLossFunction`](../../R/xbart.R) and the 0/1 refusal in [`xbart`](../../R/xbart.R), the dbartsModel class in
  [R/A_class.R](../../R/A_class.R), the sampler's `$setWeights`, `$setCounts`, `$setResponse` and `$setActiveRows`
  ([R/dbarts.R](../../R/dbarts.R)).
- Tests that pin today's weight behavior: test-binary-weight-mask.R, test-weighted-logistic.R, test-weighted-binary-ppd.R,
  test-logistic-weight-swap.R, test-zero-weights.R, test-xbart-weight.R, test-state-weight-pairing.R,
  test-state-zero-weight-rows.R, test-pointwise-loglik.R, test-nbinom.R, test-xbart-loss.R.
- TODO items: frequency weights (dec-B125), mask-on-own-data-object (closed by S4), xbart-binary-rmse (closed by S6),
  approximate-distributions (the binomial's cost case, not this plan's).

## Constraints

- No draw of a call taken today moves: every engine path a precision fit, an unweighted fit or a binary logistic fit runs
  keeps its expression and its order of random numbers. Where a new term is zero under precision, it is not added, rather
  than added as 0.0.
- Every draw stays exact; real weights on the Polya-Gamma families are refused until a real-shape draw exists (dec-B568).
- Out of scope: frequency weights on probit, ordinal, the probit hazard and aft (dec-B563); hurdle.lognormal's weights
  (dec-A15); the probit binomial (dec-B567); any approximate draw (TODO approximate-distributions); bartBT and rbart_vi, which
  take no weights.type; the shipped header (its field rides host-neutral-engine's S3, below).
- No new dependency; base R calls within DESCRIPTION's R floor.

## Design

### The reading and where it lives

- R: dbartsModel gains a slot `weights.type`, a character of length one, NA for the family's own reading (weights.type = NULL
  as given), else "precision" or "frequency". bart, dbarts, dbartsSpec and xbart take the argument and write the slot; the
  resolution against the family happens in the spec's resolution, where the family is known, so dbartsSpec's returned model
  carries the resolved value and creation reads it.
- Bridge: reads the resolved slot into a sampler-level flag, `bartcore::WeightsReading { precision, frequency }`, passed to the
  response family's constructor beside the weights. A family whose own reading is counts gets frequency.
- Engine: the flag is held by the families that read it (gaussian, Student-t, heteroscedastic through the variance leaf,
  logistic, nbinom, multinomial) and enters the weights digest, so a state stored under one reading and installed into a
  sampler of the other redraws as a weight change does (design Calls made, 5).

### Engine: the gaussian-type families (S1)

- `GaussianResponse` holds a per-row count source: null under precision, where the count is 1 per positive-weight row in, as
  now, and the user weights under frequency, where the count is the weight. The positive-row count becomes `rowCount_`, a
  double: under precision exactly `static_cast<double>(numPositiveWeights_)` as today, under frequency the sum of w over rows
  in the likelihood (weight positive, mask active), recomputed in `installWeights` and wherever the mask recomposes.
  `ChiSquaredScalePrior::drawSigmaSqFromPosterior` takes the count as a double; the precision path passes the same value it
  passes today. The SSR, the leaf statistics and the moves read the precisions w unchanged.
- Student-t: under frequency a row's composite precision is Lambda, drawn as one
  `Gamma(w (nu + 1) / 2, rate (nu + r^2 / sigma^2) / 2)` with the residual unweighted, in place of `w lambda` with
  `lambda ~ Gamma((nu + 1) / 2, rate (nu + w r^2 / sigma^2) / 2)`; the contained gaussian gets the user weights as its counts.
  The nu step reads the note's generalized statistics, `sum m_i`, `sum m_i log Lambda_i`, `sum Lambda_i`, plus the per-grid
  constant `sum_i [m_i lgamma((nu + 1) / 2) - lgamma(m_i (nu + 1) / 2)]`, computed when the weights install (at creation,
  `setWeights`, `setData`, a mask change, a state install). Under precision `drawIndex` runs its shipped expression, the
  constant array never built.
- Heteroscedastic: `ConstantVarianceLeaf::accumulate` adds w for a row under frequency and 1 under precision, in its leaf draw
  and its integrated likelihood, which the variance tree moves read. The leaf learns the reading as a member set at
  calibration; any other site that counts a variance leaf's rows (the fused scan path, the combiner's variance glue) reads the
  same member, which the implementer lists in the slice's report.
- A frequency weight of 0 is a row out of the likelihood, as now; a positive real w is taken (dec-B568).

### Engine: the Polya-Gamma families (S2)

- logistic, hazard.logistic: their weights already are counts; the slice only checks the reading flag (precision refused
  above). The engine change is zero, the zero weight being S4's.
- nbinom: under counts a row draws `PG(w (y + r), psi)` and its working response uses `kappa = w (y - r) / 2`; the r grid step
  multiplies each row's log mass by w, and `NBShapePrior::computeKernel` takes the weights. The shape `w (y + r)` is a whole
  number (r fixed is a whole number, its grid holds whole numbers, y and w are whole), drawn as that many unit draws.
- multinomial: the combiner's per-category loop draws `PG(w n_i, .)` with kappa `w (y_ik - n_i / 2)`, the counts times w; its
  log-likelihood stays R's. A zero weight composes into the combiner's effective mask beside the zero-trial rows, with no
  zero-trial warning.
- The cap: `w` times the row's shape-defining count (`w (y + r)` with r at its grid maximum for a drawn r, `w n_i`) above
  dec-B202's one million is refused at creation and at every weight or count change.

### Engine: the binomial row (S3)

- `LogisticResponse` gains an optional `trials_` pointer, null for a binary response, where it reads as 1. With trials, `y_`
  holds the successes k (for a binary row k = y, n = 1, the same values). The shape is `lround(w n)` unit draws (w the
  frequency weight, 1 absent), the working response `w (k - n / 2) / omega - offset`, written so that at null trials the
  expression is today's `w (y - 0.5) / omega - offset`, operand for operand. The cold start is `omega = n w / 4`.
- The log-likelihood: `w [k log p + (n - k) log(1 - p)] + w c_i` with `c_i = log choose(n_i, k_i)`, the per-row constants
  computed once per install of the counts from `std::lgamma`, so a binomial row reports `w log dbinom(k; n, p)` and a binary
  row today's value.
- The bridge hands successes `counts[, 1]` and trials `rowSums(counts)` from the data object's counts slot, the slot multinomial
  count rows use (design Calls made, 9); `$setCounts` on a binomial sampler swaps both and redraws as a count swap does now. A
  row of 0 trials installs into the mask as dec-B133's zero-trial row does, with its warning.
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
  does now. A weight other than 0 or 1 there keeps today's refusal.
- The R policy stops clearing the weights slot on probit and ordinal and stops returning a mask: `s$data@weights` holds the
  0/1 vector, `s$setData(s$data)` keeps the rows out, `$activeRows` is only what `setActiveRows` installed.
  `dbartsSpec()` drops its `active` element, and `new("dbartsSampler", ..., active = )` stays the creation door of a mask that
  composes with the weights.
- `LogisticResponse` takes a 0 as the mask: the row is skipped in the draws (consuming no variate, the property its mask skip
  protects) and is NaN in the log-likelihood; `enforceBinaryWeightPolicy` and `refuseNonCountWeights` accept 0. dbartsData's
  warning that zero-weight rows among other weights stay in the fit stays, now true on logistic as on every family; what goes is
  the refusal that followed it there.
- nbinom's 0/1 weights are counts from S2, a 0 out of the likelihood under the same composition.
- xbart's refusal of 0/1 weights on probit goes (design Calls made, 14): each fold's data object carries its rows' weights, the
  fold's creation installs them, and the held-out loss weights a 0 row out.

### The R surface (S5, S6)

- weights.type (S5): validated, resolved and refused as [Refusals](#refusals) words; fixed at creation, so `$setModel` with a
  model of another reading is refused, and `setWeights`, `setData` and predict read new weights under the one in force.
- Read off the expanded rows under frequency (S5, design Calls made, 3): the starting sigma's regression with residual df
  `sum(w) - p` and its stand-in, a weighted sd; the chi prior's anchor; a linear or gp leaf covariate's center and scale; a
  forest's default sd; the quantile cut rule's quantiles; nbinom's center. A stated sigma is taken as stated.
- Predictive draws (S5): `sampleFromPPD` stops reading the fit's weights for a frequency or count fit's training rows: one
  observation at sd sigma (dec-B565) or Bernoulli(p), one nbinom count, one multinomial category or count vector at the row's
  trials (dec-B566). predict(weights = ) on a frequency gaussian-type fit takes 1 and 0 only (dec-B565); on a count family it
  asks for that many trials, at training rows and at newdata (dec-B566 revising dec-B470).
- The log-likelihood (S5): `w log f(y_i | one observation)` under frequency (dec-A116, dec-B565), NaN at a 0 on every family
  (dec-B570); the engine's channel and `pointwiseLogLikelihood` agree, which a test pins.
- The binomial (S6): the formula parse takes `cbind(k, f)` under family = binomial, stats::binomial() and "logistic", and the
  matrix interface takes a two-column `y` under those families; dbartsData takes the pair as `counts =`; the "auto" refusal and
  the probit-binomial refusal are reworded; the proportion spelling is refused (dec-B573). fit$y keeps the two-column response;
  fitted(type = "ev") is p and its residual `k / n - p`; fitted(type = "ppd") is `n p` and its residual `k - n p`; extract and
  predict at training rows draw Binomial(n, p) (dec-B574); predict at newdata draws Binomial(w, p) with w from its weights, one
  trial without (design Calls made, 10); the zero-trial row is dec-B133's (design Calls made, 13).
- xbart (S6): the built-in losses on a logistic or probit fit follow design Calls made, 12, one rule on the probability scale
  (observed proportion y against the posterior mean probability p, the link applied per draw and averaged last, each row
  weighted by its trials times its weight); a binomial xbart receives `cbind` as bart does; a user loss receives `y.test` as
  `k / n`, the latent draws, and `weights` as trials times weight. On a binary fit this changes only rmse.
- Help and NEWS (S5, S6): weights and weights.type on bart, dbarts, dbartsSpec and xbart; the help's cautions (design Calls
  made, 8; dec-B568's rescaled and survey weights); the two scales of a binomial fit; NEWS states the moved training draw of a
  weighted logistic fit, the moved zero-weight log-likelihood, the binary xbart rmse, and the binomial.

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
4. aft and hurdle.lognormal keep `aft (survival) models do not support case weights` and the hurdle refusal; "frequency" there
   reaches the same words.
5. logistic and hazard.logistic, a weight negative, fractional or NA: `logistic weights are observation counts and must be
   non-negative whole numbers; use a gaussian model with weights.type = "frequency" for real-valued counts`.
6. nbinom: `nbinom weights are observation counts and must be non-negative whole numbers; an exposure belongs in the offset as
   a log-exposure term`.
7. multinomial: `multinomial weights are observation counts and must be non-negative whole numbers`.
8. The cap, through `refuseCountOverCap`'s form: `<caller>: family "logistic" requires 'weights' times the trials no larger than
   1000000; each unit of a count is a Polya-Gamma draw per sweep`, and the same with `"nbinom"` and `'weights' times (y + r)`,
   `"multinomial"` and `'weights' times the 'counts' row totals`.
9. predict on a frequency gaussian-type fit with weights other than 0 and 1: `on a fit with weights.type = "frequency" predict
   takes weights of 1, a new observation, and 0, a row held out: a row of weight w stands for w observations, and its draw is
   one of them`.
10. `$setModel` with another reading: `$setModel cannot change weights.type, which is fixed when the sampler is created;
    create a new sampler`.
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

Gone: xbart's `xbart does not accept weights of 0 and 1 under family ...`; logistic's refusal of a 0; the logistic refusal of
a two-column response under family = binomial and "logistic"; the multinomial refusal of weights other than all ones.

## RNG class

Neutral on the engine for every call taken today, shown three ways:

- Bitwise harnesses, run against each engine slice's own `--preclean` library: `benchmarks/R/equivalence.R compare`,
  `bcf-equivalence.R compare` and `multinomial-equivalence.R compare`, every scenario "identical draws (same RNG stream)",
  counted per [Gate hygiene](README.md#gate-hygiene); the four `test-reproducibility-*.R` on the reference build, unchanged
  (two of them carry precision weights, which pins the precision path); tests/cpp unchanged in every existing test.
- Construction: under precision every changed expression is today's, operand for operand (the count is
  `static_cast<double>(numPositiveWeights_)`, the variance leaf adds 1.0, the nu constant is not built, a null `trials_`
  reads as today's working response); each slice's reviewer reads the diff for it and runs the mutants that would break it.
- The membership paths (S4): a probit and an ordinal fit with 0/1 weights equal today's draw for draw, the mask installed by
  the engine at creation against today's install by R after creation, which draws no variate (no row reactivates); pinned by
  test-binary-weight-mask.R unchanged in its draws.

R-level values that move, by ruling, each re-pinned in the test that holds it and stated in NEWS:

- A weighted logistic or hazard.logistic fit's training predictive draws (extract, fitted(type = "ppd"), predict at training
  rows without weights): Bernoulli(p) where they drew Binomial(w, p) (dec-B566). R's random stream after such a draw differs
  too, since rbinom consumes by size; test-weighted-binary-ppd.R and test-weighted-logistic.R re-pin, replayed as whole files.
  No snapshot file holds a weighted logistic draw (checked: their weights are gaussian).
- A logistic zero weight's pointwise log-likelihood, refused today, is NaN (dec-B570).
- A binary xbart's rmse loss (S6, TODO xbart-binary-rmse): moves to the probability scale; its default, log, does not move.
  test-xbart-loss.R re-pins; test-reproducibility-xbart.R is gaussian and does not move.

New calls take new draws, each held to an oracle: the expanded rows draw for draw (gaussian, S1; binomial, S3), the expanded
rows in distribution (Student-t, heteroscedastic, nbinom, multinomial, S1 and S2), the exact-posterior arms of S7. The
exact-gates list grows by one script (S7); no existing gate's arm changes.

## Calls made

1. The engine form of the binomial is successes in `y` and a `trials` vector, not the proportion k / n with the trials as a
   count: kappa `w (k - n/2)` is then exact, which the draw-for-draw gate against the expansion needs, and the spec a C caller
   fills holds whole numbers. Alternative: y = k / n, where `n (k / n)` need not round back to k.
2. The log choose constants are per-row values computed at each install of the counts, not evaluated in the log-likelihood
   loop. Alternative: `std::lgamma` three times per row per stored draw.
3. The frequency count reaches the gaussian family as a second pointer (counts beside precisions), null under precision; a
   Student-t sampler hands its contained gaussian the user weights as counts and Lambda as precisions. Alternative: a count
   derived from the precisions, which Student-t's composite breaks.
4. Probit and ordinal compose the 0/1 weights with the mask in the family, as gaussian composes weights and mask, so a C caller
   gets membership from the weights with no R policy. Alternative: the bridge composes and installs the product through
   `setActiveRows`, which leaves the engine unaware of the weights and host-neutral's build without the rule.
5. `dbartsSpec()` drops `active`. A caller reading `spec$active` gets NULL and `setActiveRows(NULL)` installs nothing, so
   stan4bart's code stays correct unchanged and its plumbing is dead. Alternative: keep `active` as a copy of the 0/1 weights,
   which a caller following today's help would install twice, harmless but misleading.
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
   reworded; on bart, dbarts and xbart the pair is the response under the logistic family. Alternative: dbartsData taking a
   two-column y as counts, which it cannot tell from a (time, status) pair without the family.
10. The binary rmse fix lands inside this feature's S6, not alone on main: it shares S6's loss rewrite, and dec-B494 asks a
    feature to merge once. Alternative: a one-slice R-surface landing on main now, about 60 lines, which S6 then rewrites; the
    coordinator may prefer it to retire a defect sooner.

## Open calls

None. Every user- or consumer-visible choice the plan makes is settled by dec-B563 to dec-B570, dec-B573, dec-B574 or the
design note's calls made; the plan's own calls above are internal or follow from a ruling (Calls made, 5, from dec-B569).

## Slices

Each slice is a commit series on frequency-weights, its own review and gates; the feature merges once after S8.

### S1: engine, the gaussian-type families (~170 code, ~250 tests/cpp, ~80 bridge and R plumbing)

Tier "Changes draws"; RNG class neutral on every call taken today, new draws under frequency.
Steps: `WeightsReading` and its bridge read from the model slot (the slot and its R plumbing land here, refusals in S5); the
gaussian count pointer and `rowCount_`; `drawSigmaSqFromPosterior` taking a double count; Student-t's frequency draw, its
counts, the nu constant; the variance leaf's count and every site that counts its rows; the digest covers the reading.
Gates: tests/cpp new tests against the expanded data (sigma, the fits and nu's grid probabilities on fixed residuals, the
note's p8 arm, within Monte Carlo error) and the w = 1 identity (frequency at all-ones weights equals precision bitwise); the
bitwise trio; the full tinytest suite; sanitizers on tests/cpp and on the R-loaded path; `mutation-battery.R verify-anchors`.

### S2: engine, the Polya-Gamma families (~90 code, ~200 tests/cpp, ~40 bridge)

Tier "Changes draws"; neutral today (nbinom and multinomial weights were refused or masks).
Steps: nbinom's shape, kappa and r-step power; the combiner's multinomial shape and kappa; the cap on the products; the
bridge stops refusing nbinom and multinomial weights.
Gates: tests/cpp: an nbinom row of weight w against w copies in distribution (a fixed-r arm, where the copies' omega sum is
the row's in distribution, and a drawn-r arm); a multinomial row of weight w against the counts times w, draw for draw (the
same unit draws in the same order); the trio; the suite; sanitizers.

### S3: engine and bridge, the binomial row (~110 code, ~150 tests/cpp, ~50 bridge)

Tier "Changes draws"; neutral for a binary response.
Steps: `trials_`, the working response, the cold start, the log choose constants and the log-likelihood; the bridge reads the
counts slot under logistic; `$setCounts` reaches a binomial sampler; the zero-trial row into the mask.
Gates: tests/cpp: a binomial row against two adjacent rows of weights k and n - k for 0 < k < n, fits within 1e-10 over 300
sweeps, the same number of variates consumed; a binomial row at frequency weight w against the row of k w successes in n w
trials, bitwise; the log-likelihood against R's `dbinom(k, n, p, log = TRUE)` times w within 1e-10 relative; a null-trials
sampler bitwise today's; the trio; the suite. Today a weight-0 logistic row still draws one variate, so k = 0 and k = n join
the draw-for-draw test in S4, once the expansion's weight-0 row is masked and draws none.

### S4: membership and zero weights (~60 engine, ~80 R, ~150 tinytest)

Tier "Changes draws" (the mask paths of three families change); neutral on every call taken today.
Steps: probit and ordinal take weights and compose; logistic takes a 0; the R policy keeps the weights slot; `dbartsSpec()`
drops `active`; xbart's refusal goes; help for `$setWeights`, `$setActiveRows`, dbartsSpec.
Gates: the binomial draw-for-draw test of S3 at k = 0 and k = n; test-binary-weight-mask.R's draws unchanged (its
`s$data@weights` and `spec$active` assertions re-pinned to the ruling); a new test that `s$setData(s$data)` keeps masked rows
out and that a mask set by `setActiveRows` composes with the weights in either order; a logistic fit with a zero weight equal
draw for draw to the fit with that row masked by `setActiveRows`; xbart on probit with 0/1 weights equal, fold for fold under
a fixed fold assignment, to a by-hand loop fitting each fold with those weights and scoring its held-out rows of weight 1 (the
test-xbart-fold-oracle.R pattern); the suite; the trio.

### S5: R surface, the reading (~380 R, ~100 help, ~450 tinytest)

Tier "R surface". Moves the weighted logistic training draw and the logistic zero-weight log-likelihood, re-pinned.
Steps: weights.type on bart, dbarts, dbartsSpec and xbart, its validation and resolution; Refusals 1 to 10; the quantities
read off the expanded rows; `sampleFromPPD`, predict's doors, the log-likelihood; help and NEWS.
Gates: the touched test files and the new test-frequency-weights.R (each refusal by its words, counted warnings, the
expanded-data identity at the R level for sigma's starting value and the chi anchor, predict's doors, the log-likelihood
against the engine's channel); lint_package, air, codoc, win drift, doc-freshness; the snapshots on the reference build.

### S6: R surface, the binomial and xbart (~270 R, ~50 help, ~400 tinytest)

Tier "R surface". Moves the binary xbart rmse, re-pinned.
Steps: the parse and routing; Refusals 11 to 16; fit$y, fitted, residuals, extract, predict's trials, the zero-trial row; the
log-likelihood channel; xbart's losses and its binomial response; help and NEWS.
Gates: test-binomial.R (against glm on the same data: fitted(type = "ev") near glm's fit on a large-n design, residuals'
definition exactly, the log-likelihood's constant against `dbinom`, predict's trials, the zero-trial row's NaN-free 0 and its
-p residual, every refusal by its words); test-xbart-loss.R re-pinned with an arm showing a binomial xbart equals xbart on
the trial-level expansion for the log and misclassification losses exactly; the R-surface checks of S5.

### S7: gates (~350 in benchmarks/R, ~20 workflow)

Tier: gate-only, the R-surface review.
Steps: benchmarks/R/frequency-exact.R with arms: gaussian frequency against the expanded rows draw for draw over 300 sweeps
(the note's p2 emulation as the seed, sigma within 1e-12 and the fits within 1e-12); binomial against its two-row expansion,
the same bound; an exact-posterior arm for sigma on a single-leaf tree under real weights, against the analytic
inverse-chi-squared posterior; Student-t's nu on fixed residuals against its exact grid probabilities; heteroscedastic,
nbinom and multinomial against the expanded rows by Welch z on the fits, quick and long. It joins exact-gates.yaml's list.
An SBC arm in sbc.R for a frequency gaussian fit at real weights (dispatch, not the matrix, until its ladder earns admission).
Gates: each arm run green, and its discriminating mutant (below) run red.

### S8: records (~60)

The plan's Status and Landing note, the ledger entry for the calls made, TODO (frequency weights, mask-on-own-data-object,
xbart-binary-rmse), the design Status line, the indexes.

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

1. S1: count each positive row as 1 under frequency (sigma's df) - frequency-exact.R gaussian arm, tests/cpp expanded test.
2. S1: the frequency count under precision too (sigma's df the sum of the weights, not the positive rows) -
   test-reproducibility-continuousResponse-singleThreaded.R on the reference build (precision weights 1 and 2) and the
   equivalence harness's weighted scenarios, bitwise.
3. S1: drop the lgamma constant under frequency - tests/cpp nu grid test (a probability off by up to 0.963 in the note's p8).
4. S1: variance leaf adds 1 under frequency - the heteroscedastic expanded arm.
5. S1: the digest ignores the reading - test-state-weight-pairing.R's new arm (a state under frequency into a precision sampler
   redraws).
6. S2: nbinom kappa without w - the nbinom expanded arm.
7. S2: multinomial shape `n_i` not `w n_i` - the multinomial draw-for-draw test.
8. S3: kappa `w (k - 1/2)`, the binary expression on the successes - the binomial draw-for-draw test.
9. S3: omit the log choose constant - the log-likelihood against dbinom.
10. S3: the shape `lround(n)` without the frequency weight - the weighted binomial row against k w in n w, bitwise.
11. S4: probit composition takes the mask alone - the setData test and the either-order composition test.
12. S4: logistic draws a 0-weight row - the zero-weight logistic against the masked fit, draw for draw.
13. S5: `sampleFromPPD` keeps Binomial(w, p) on a training row - test-weighted-binary-ppd.R re-pinned.
14. S5: predict on a frequency fit takes weights 2 - its refusal test.
15. S6: rmse on the latent scale - test-xbart-loss.R's binary rmse arm.
16. S6: the binomial residual `k - p` - test-binomial.R.

## Consumers

- stan4bart (main 27ef163): passes its weights to the BART part at weights.type NULL, a precision on its continuous fits as in
  its own Stan target; its probit rule passes 0/1 weights, which now ride `data.bart` into every sampler it builds, and its
  `active.bart` from `bart_spec$active` becomes NULL, so its `setActiveRows` installs are skipped and `restoreBartSampler`
  rebuilds from the data, which carries the weights (`data.bart.light` drops only x). No call changes; an optional cleanup of
  about 30 lines drops the dead `active` plumbing, in the consumer sweep, not this feature. Its `install_bart_range` reads
  `sampler$data@weights` on continuous fits only, unchanged.
- bartCause (1e4b882): passes population weights to bart's gaussian response fits at weights.type NULL, unchanged; drops
  weights from its propensity fits; its xbart call (bayesOpt) uses the default loss, which on a binary treatment is log and does
  not move. No call changes.
- treatSens (268b2db): passes no weights to dbarts; no change.
- bairrtt: reads none.

## Interaction with the other pending plans, and order

- forest-door lands on main shortly, then coefficient-prior-live (its branch reworks chain.hpp, combiner.hpp, model.hpp,
  R_interface_bartcore.cpp, R/dbarts.R and R/spec.R, all of which this plan touches); engine features merge serially. This
  feature branches from main after coefficient-prior-live lands, or rebases onto it, re-gating what the rebase touches.
- host-neutral-engine: its S1 (build-time math, branch host-neutral-math) is independent of this plan and lands when green.
  Its main feature rewrites every creation route's SEXP reads into a plain-C spec (S2b), which carries `weightsType` per
  dec-B564 and, under this plan, `trials` beside y under LOGISTIC (Design, the binomial row). Recommended order: this feature
  merges before host-neutral's S2a starts, so S2b translates the weights reading, the counts slot under logistic and the
  membership rule as code that exists, and S3's header change carries the field once (dec-B564: "the shipped header gains
  one field in dec-B85's change"). Cost: S2a waits for this feature, about 130 to 175 agent-hours by the record, much of it
  parallel to coefficient-prior-live's tail if S1 and S2 start on a branch of forest-door's main. The other order makes this
  plan write its bridge reads into the spec form (S2b's) and add `trials` to the header itself, a second hash move. Nothing
  here edits dbarts.h.
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

About 3,500 lines, 1,330 of them package code; about 87 agent-hours estimated, 130 to 175 by the record of plans here. No
per-sweep cost on the gaussian-type families at either reading (the count is fixed between weight installs, the variance leaf
adds w in the pass that counts rows now, the nu constant is per grid point); the count families' cost is today's per unit
draw, multiplied by the counts the user asks for.

## Not verified

- The full list of sites that count a variance leaf's rows (the fused scan path, the combiner's variance glue); S1's
  implementer lists them and the mutant 4 run shows the list complete.
- That `std::lgamma`'s log choose agrees with R's `dbinom` constant to 1e-10 relative at n up to one million; S3's test is the
  check.
- That the multinomial combiner forms kappa in one place for both the train and a state install; S2 reads it.
- The Student-t nu step's generalized statistics were checked by the critic's p8 within Monte Carlo error, not derived anew
  here.
- The cost figures were measured once on a shared, loaded laptop; bench-sampler.R on a quiet machine is the gate for the hot
  paths, not these.
- How coefficient-prior-live leaves model.hpp's gaussian family and the combiner; its branch was read for the files it touches,
  not line by line.
- stan4bart's, bartCause's and treatSens's call sites were read on main 27ef163, 1e4b882 and master 268b2db; their 1.0 compat
  branches may differ.
