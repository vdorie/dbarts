# Frequency weights and the binomial response

Status: DECIDED 2026-11-03, design note before the plan; Calls made 5 and the state and hurdle rows of What moves
amended 2026-11-05 for dec-B586 and dec-B588, and Calls made 10 on 2026-11-06 for dec-B594. Rulings, in docs/decisions.md: dec-B125 (frequency
weights join precision weights by explicit choice; precision stays the default), dec-A11 (under precision weights
sigma's degrees of freedom count only positive-weight rows), dec-B120 (the comparison against 0.9-34 is a release
gate), and the note's own rulings dec-B563 to dec-B570, dec-B573, dec-B574 and dec-B580, each cited where it applies. The
plan is docs/plans/frequency-weights.md.

## Summary

A weight in dbarts meant one of three things, depending on the family: a precision on gaussian and Student-t
fits (the row's noise variance divided by the weight, as in R's lm), a count of identical observations on
logistic fits, and a 0 or 1 saying whether a row is in the data on probit, ordinal and negative binomial fits.
The three are really two. Both readings raise a row's pull on the fit by its weight; they differ in how many
observations the noise is estimated from, once per row under precision, as many times as the weight under
frequency. Families with no noise level to estimate have one reading only, and 0 and 1 mean the same under
both: out, or in once.

The maintainer ruled the following. A frequency-weighted fit is the fit to the data with every row written out
as many times as its weight. It is chosen with one argument beside the weights on the fitting functions, left
unset for the family's own reading. The gaussian, Student-t and heteroscedastic families take frequency weights
of any positive size, at no extra cost per sweep. Logistic, negative binomial and multinomial take whole
numbers, since their one hidden draw per row simply gets larger, at a cost that grows with the weight. Probit,
ordinal, the probit hazard family and survival fits keep what they take today, because each copy would need a
hidden value of its own. A frequency-weighted row's predictive draw is one new observation like the row, on
every family. Probit and ordinal keep 0/1 weights for membership, now kept on the sampler's data object. A zero
weight means the same everywhere, logistic included.

The maintainer also ruled that a binomial response joins the logistic family before 1.0-0: a row recording k
successes in n trials, written as R writes it, cbind(successes, failures). Its draw is the one logistic already
makes, so the work is mostly in the R surface. Wanting it with many trials is what opened the research into
approximate draws, which the package had deferred.

Two choices about the binomial response were then ruled. glm's second spelling, a proportion with the trial
count as weights, is refused with a message pointing to cbind, because on a 0/1 response it would collide with a
frequency weight. A binomial row's fitted values and residuals are on glm's scale, the probability of a success
and the observed proportion minus it, while its predictive draws stay counts.

## What a weight meant at the merge

| family | values taken | read as | a 0 |
|---|---|---|---|
| gaussian, Student-t, heteroscedastic | finite, non-negative | precision: noise variance sigma^2 / w | out of the likelihood, in the design |
| logistic, hazard.logistic | positive integers up to 1e6 | copies: Polya-Gamma shape w | refused, after a warning that it stays in the fit |
| probit, ordinal, nbinom, hazard and hazard.probit | 0 and 1 | the active-row mask | out of the likelihood, in the design |
| multinomial | all ones | refused otherwise | refused |
| aft, hurdle.lognormal | none | refused | refused |

The routing is [`enforceWeightPolicy`](../../R/spec.R), with [`isMaskedWeightFamily`](../../R/spec.R) naming the
mask families and [`refuseNonCountWeights`](../../R/spec.R) the logistic rule; hazard.logistic is remapped to
logistic before the policy, so it takes counts. Run on c3d268c8f: a logistic training predictive draw at weight 5
ranged over 0 to 5, Binomial(w, p); its pointwise log-likelihood was w times the Bernoulli log density (dec-A116);
a probit fit refused weights of 5 while its predict(type = "ppd", weights = 5) drew 0 to 5 (dec-B470); hazard and
hazard.probit refused weights of 5 and hazard.logistic took them (critique probe p10); and cbind(k, n - k) was
refused under family = binomial, "logistic" and "auto", each message pointing to the multinomial family, while a
proportion with weights = n was refused as not 0/1. rbart_vi, deprecated and removed at 1.1-0, keeps precision
weights and takes none of what follows.

0.9-34's probit weights were neither reading: its latent was a truncated normal of sd 1/sqrt(w) about the fit,
so it fitted P(y = 1) = Phi(sqrt(w) f(x)), a rescaled link per row (classic engine, bartFit.cpp at the 0.9-34
tag). dec-B13 removed that, and its 0.9-x error stands.

How others read weights, from installed packages and their sources, run where marked:

| fitter | reading | non-integer | run |
|---|---|---|---|
| lm, glm gaussian, mgcv gam, lme4 lmer | precision; glm: "each response y_i is the mean of w_i unit-weight observations" | taken | lm, glm, gam: sigma differs from the expanded-data fit |
| MASS rlm | wt.method = c("inv.var", "case") | taken | "case" equals the expanded-data fit |
| glm binomial, Poisson | counts; binomial as cbind(k, n - k) or as k / n with weights = n, one likelihood | warns "non-integer counts" | the two binomial spellings give one fit and one log-likelihood |
| survival survreg, coxph | case weights, "treated as replication weights" | taken; coxph then defaults to a robust variance | equal the expanded-data fit |
| MASS polr, nnet multinom | case weights | polr warns from its starting glm | polr equals the expanded-data fit |
| rstanarm 2.32.2 | log_lik times the weight; posterior_predict one observation per row | taken | a gaussian identity stan_glm with no offset ignores the weights entirely |
| brms | `weights()` multiplies the log-likelihood; `y | trials(n)` the binomial | positive values | not run |
| survey | design weights, rescaled to mean 1 inside svyglm | taken | |
| BART (wbart), bartMachine | `w` multiplies the noise sd; no weights | | |
| Stata, SAS, statsmodels | fweight and aweight, FREQ and WEIGHT, freq_weights and var_weights | Stata refuses non-integer fweights | not run, from their documentation |

## The two readings differ only in the scale's count

In glm's form a row's likelihood contribution is a kernel times a normalizer that may depend on a scale. A
precision weight multiplies the kernel's exponent: for a gaussian, exp(-w r^2 / 2 sigma^2) with one factor of
1/sigma per row of positive weight, dec-A11's count. A frequency weight raises the whole contribution to the
power w: the same kernel and w factors of 1/sigma. So on a family with no estimated scale the two readings are
one likelihood, and 0 and 1 are both readings everywhere. On a family with a scale, the mean forest's leaf
updates and moves are unchanged given the row precisions, the residual scale's posterior counts sum(w)
observations, and on a heteroscedastic fit the variance leaves count the same way, in their draws and in the
integrated likelihood their tree moves read ([`ConstantVarianceLeaf`](../../src/bartcore/model.hpp)).

The target is the fit on the expanded data, every row written out w times. The tree machinery does not separate
them: the empty-leaf rule counts rows reaching a leaf, the default cut rule places one point per gap between
distinct values, and the moves choose among the values a node holds. Every quantity the package reads off the
rows does, so those are read off the expanded rows (Calls made, 3). Emulated on c3d268c8f, a precision fit
whose chisq prior carries the frequency degrees of freedom and whose sigest is the expanded data's matched the
expanded-data fit draw for draw: sigma to 2.2e-16 and the fit to 2.9e-15 over 300 sweeps, 40 rows of weights 1
to 4 (critique probe p2, rerun). A gaussian sweep's random numbers do not depend on how many times a row
appears, which makes that a gate.

The same probe shows what frequency means for a continuous response: the expanded data's sigma settled at 0.05
against a true 0.4 (the precision fit 0.60), because w copies of y at one x are exact ties. The help says so
(Calls made, 8). A zero weight is the one place the target is not the expanded data: the row stays in the
design (dec-B238, dec-B251).

## Family by family

The sampler draws a hidden value per row instead of evaluating the likelihood, so the question for each family
is whether the likelihood to the power w keeps a one-draw form (dec-B563).

| group | families | w copies | draw | frequency weights |
|---|---|---|---|---|
| gaussian-type | gaussian, Student-t, heteroscedastic | N^w, t_nu^w | sigma's count sum(w); the copies' Student-t scales sum to one gamma draw; a variance leaf counts w | real w, exact, no cost per sweep |
| Polya-Gamma | logistic and hazard.logistic, nbinom, multinomial | the power multiplies the Polya-Gamma shape by w | one draw per row: shape w for a logistic row, w (y + r) for nbinom, w times each category's trials for multinomial; the shape is drawn as that many unit draws | whole numbers; a sweep grows with w |
| truncated-normal | probit, ordinal, the probit hazard, censored aft rows | Phi^w and its kin, with no one-draw form | one latent per copy | not taken (dec-B563): they keep what they take |

On the Polya-Gamma families nothing about the latent changes but its shape, exactly as logistic's counts
already work ([weighted-logistic.md](weighted-logistic.md)); nbinom's shape grid and multinomial's softmax
margins read the likelihood and take the power. At the cap of dec-B202 a sweep is long: about 160 ns per unit
draw, 0.16 s a sweep for one row of a million and 16 s for 100 such rows, and 5.2 s a sweep for 100
multinomial rows of 100,000 trials over three categories. Real w there waits on a real-shape draw, dated after
1.0-0 on 2026-08-28 and now inside the approximate-draw research (dec-B567). hurdle.lognormal's weights stay
with its feature completion (dec-A15), which takes 0 and 1 on both parts (dec-B588).

Student-t's derivation. For w copies the copies' mixing scales are iid Gamma((nu + 1) / 2, (nu + r^2 / sigma^2)
/ 2) given the residual, and the gaussian part reads only their sum Lambda_i, Gamma(w (nu + 1) / 2, same rate),
so the trees and sigma given nu are exact at any real w, sigma reading sum(w) observations and the sum of
Lambda r^2. With m_i the row's count (w_i under frequency, 1 under precision, rows out of the likelihood
skipped), the nu step on its grid reads sum(m_i) for the row count, sum(m_i log Lambda_i) and sum(Lambda_i) for
the lambda statistics, plus sum_i [m_i lgamma((nu + 1) / 2) - lgamma(m_i (nu + 1) / 2)], a constant of the
weights per grid point computed when they are installed. At m_i = 1 the term is zero and the draw is the shipped
one ([`ResidualDfPrior::drawIndex`](../../src/bartcore/model.hpp)). Critique probe p8 ran the step on fixed
residuals against the exact grid probabilities of nu, at integer and real counts: within 0.0026 (Monte Carlo sd
about 0.006); without the lgamma term a grid probability is off by up to 0.963.

## Settled by the rulings

### The families (dec-B563)

Frequency weights reach the gaussian-type and the Polya-Gamma families. The maintainer, 2026-11-03, on the question
whether every family is as simple as multiplying the log-likelihood: "Let's do (b) plus the Polya-Gamma families."
Probit, ordinal, the probit hazard (hazard and hazard.probit) and aft keep what they take today, and the ruling
does not schedule them. hazard.logistic, a logistic fit underneath, is a Polya-Gamma family and takes frequency
weights with logistic.

### The argument, weights.type (dec-B564)

bart (formerly bart2), dbarts, dbartsSpec and xbart gain weights.type, after MASS's rlm wt.method: NULL by
default, the family's own reading (precision where there is a residual scale, counts where there is none);
"frequency" taken on the families of dec-B563 and refused by name elsewhere, except 0/1 weights, which read the
same either way; "precision" refused by name on a family with no residual scale. It is fixed at creation and
held on the model record, so setWeights, setData and predict read new weights under it; dbartsData keeps one
weight vector; the shipped header's specification gains one field in dec-B85's change. bartBT and rbart_vi have
no such argument.

### A frequency row's predictive draw (dec-B565)

One new observation like the row, at sd sigma, as rstanarm and brms draw, so predict at a training row agrees
with predict at the same x as newdata. On a frequency fit predict keeps dec-B472's doors (weights = 1 a new
observation, weights = 0 a row held out as NA) and refuses other values by name. The log-likelihood stays w
times the one-observation value (dec-A116).

### The count families' training draw (dec-B566)

A logistic training row of count w draws Bernoulli(p), and the other count families one observation like the
row: an nbinom count, a multinomial category or, on a count row, a count vector at the row's trials (dec-B177).
predict(weights = w) keeps dec-B470's trial counts; the ruling revises dec-B470's clause that predict's weights
with no newdata "replace the fit's own" on the count families, where they now ask for that many trials.

### Real weights (dec-B568)

Under weights.type = "frequency" the gaussian-type families take real weights; the Polya-Gamma families take
whole numbers until a real-shape draw exists. The help cautions that weights rescaled to sum to 1 leave sigma one
observation's worth of data, and that survey or inverse-probability weights are not counts.

### Membership (dec-B569)

Probit and ordinal are the families left whose weights are membership only. Their 0/1 weights stay the spelling
and ride the data object as weights, the engine installing them as the mask, so `s$data@weights` is no longer
NULL there and `s$setData(s$data)` keeps the rows out; `setActiveRows` stays the door of the code that runs the
BART sampler inside a larger sampler. That resolves TODO mask-on-own-data-object. xbart's refusal of 0/1 weights
on those families goes (Calls made, 14).

### Zero weights (dec-B570, dec-B580)

A 0 is out of the likelihood and in the design on every family and under both readings, NaN in the pointwise
log-likelihood and NA in the predictive draws; logistic accepts it as the mask, and its contradictory warning
goes. The premise that zero weights already stood so everywhere held on gaussian, Student-t and aft only: a probit,
ordinal, hazard or nbinom fit drew predictive values at a 0 row. dec-B580 corrects it: those rows' draws are NA on every
family, drawn then overwritten so the random stream is unchanged, the fitted value still reported. The maintainer,
2026-11-03: "(a) is correct, right? 0 weight means out of sample, and we're doing what is standard for other R
functions?"

### A proportion with the trials as weights (dec-B573)

glm takes a binomial response as cbind(k, n - k) or as the proportion k / n with weights = n, one fit either way
(critique probe p9: one log-likelihood, -42.118), and glm reads binomial weights as trials throughout, its simulate
drawing proportions of 5 trials from a 0/1 response at weights 5. rstanarm 2.32.2 switches on the data: trials
when a proportion is fractional, Bernoulli draws at a 0/1 response with weights. In dbarts a logistic weight is a
frequency weight whose predictive draw is one observation (dec-B566), so the two readings share the posterior and
differ in the draws. Under family = binomial a response strictly between 0 and 1 that is not 0/1 is refused, the
message naming cbind(k, n - k) and how to form it from a proportion and its trials, round(p * n). Shown (a) that,
recommended; (b) rstanarm's reading, one call two models by the data, which a sampler swapping its response
between draws could change mid-chain; (c) glm's, weights on binomial always trials, which revises dec-B566. The
maintainer, 2026-11-03: "Use (a)."

### The scale of a binomial row's fitted values and residuals (dec-B574)

fitted(type = "ev") is the probability p and its residual k / n - p, glm's fitted values and type = "response"
residuals (glm's own default residual is the deviance type; dbarts's residuals are response residuals on every
family); the predictive draws are counts, Binomial(n, p), so fitted(type = "ppd") is their mean, n p, with residual
k - n p. The help states the two scales, as glm's does; rstanarm's fitted values are the probability, brms's
expected value the count. Shown (a) that, recommended, against (b) counts throughout. The maintainer, 2026-11-03:
"Use (a)."

## The binomial response (dec-B567)

A row recording k successes in n trials has the binomial likelihood, whose Polya-Gamma augmentation is the one
logistic draws: one draw per row of shape n, centered at k - n/2, the same posterior as two frequency rows of
weights k and n - k. The maintainer, 2026-11-03: "Yes, have the binomial added. It can look like similar R
interfaces." The probit link's binomial is a truncated-normal family and waits with them.

The ruled surface:

- The response is written cbind(successes, failures) on the left of a formula, or a two-column matrix on the
  x/y interface, under family = binomial, stats::binomial() or "logistic"; stats::binomial("probit") with such
  a response is refused by name, saying the probit binomial is not taken.
- family = "auto" still refuses a two-column response, which it cannot tell from a (time, status) pair, and its
  message names the binomial beside the multinomial and survival readings; the logistic refusal of a
  two-column response goes.
- One predictive draw like the row, Binomial(n, p), a count (dec-B566).
- A frequency weight multiplies n and k: shape w n, the row's log-likelihood times w, its predictive draw still
  Binomial(n, p).
- predict's trial counts through weights (dec-B470): at newdata a row draws Binomial(w, p).

What the rulings leave:

- dec-B573 and dec-B574 settle glm's proportion spelling and the reporting scale (Settled by the rulings).
- Calls made 9 to 13 settle the representation, predict's default trial count, the log-likelihood, xbart's
  losses and the refusals.

The engine change is small: [`LogisticResponse`](../../src/bartcore/model.hpp) already draws the shape from its
count and centers the working response at w (y - 1/2), so a row handed y = k / n and the count n (times its
frequency weight) is the binomial row; what changes is that y may lie between 0 and 1 where a count is given,
and the predictive and log-likelihood code reads n.

## Calls made

1. The engine carries one reading flag per sampler and derives a per-row count from it, 1 for a positive weight
   under precision and w under frequency, used by sigma's degrees of freedom
   ([`ChiSquaredScalePrior::drawSigmaSqFromPosterior`](../../src/bartcore/model.hpp), which takes the positive
   count today), a variance leaf's count, Student-t's shape and nu step, and the Polya-Gamma shapes.
   Alternative: a second per-row vector of counts beside the precisions, as statsmodels allows both at once; no
   user case and n more doubles.
2. Student-t's nu step is the shipped grid draw generalized as derived above. Alternative: the marginal-t grid
   draw, also exact and better mixing, but a second mechanism and a changed draw at w = 1.
3. Under frequency every quantity the package reads off the rows is read off the expanded rows: the starting
   sigma's regression with residual df sum(w) - p ([`residualStandardError`](../../R/utility.R) reads lm.wfit's
   positive-row df today), its stand-in a weighted sd, a linear or gp leaf covariate's center and scale, a
   forest's default sd and nbinom's center. The quantile cut rule is not among them: its grid is built over unique
   values, which copies leave as they are (plan critique p5). Where sum(w) - p is not positive the stand-in is
   used. The data object's sigma is derived under the model's reading; a stated one is taken as
   stated. Alternative: read them as precision does, which breaks the draw-for-draw gate.
4. The reading is model: held on the model record, fixed at creation, and a setModel that changes it is
   refused. A fit object saved before it reads the family's own reading. Alternative: on the data object, where
   setWeights or setData could change the model.
5. A saved state carries no reading. As first written, a state's weights digest covered the reading so that a state
   installed into a sampler of the other reading redrew the latents a weight change redrew (Student-t's
   reapplyWeights, since retired); under dec-B586 every sweep draws its augmentation variables first, the digest and
   those redraws retired with latents-first (dec-A263), and a state continues under
   the reading of the sampler it is installed into. Alternative: refuse such an install.
6. A row of frequency weight 0 reports NaN in the pointwise log-likelihood
   ([`pointwiseLogLikelihood`](../../R/generics.R)), as dec-B570 states. Alternative: 0, which loo would count as
   an observation.
7. xbart keeps a row's copies in one fold, weights the held-out loss by w as now
   ([`xbartLossFunction`](../../R/xbart.R)), and takes each fold's starting sigma under the fit's reading.
   Alternative: splitting copies across folds, which leaks identical observations into training.
8. The help states: copies of a continuous response are exact ties, read as nearly noiseless data; under
   Student-t the latents the sampler reports
   ([`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) and the R method) are the row's composite
   scale; a 0/1 response's weight is a frequency, so a binomial row is written with cbind. Alternative:
   documenting only the arithmetic.
9. A binomial response is held as two count columns in the data object's counts slot, the slot multinomial
   count rows use, and the fit keeps the two-column response as given. The sampler replaces it as its formula
   writes it (ruled, dec-B628): `$setResponse` takes the n x 2 successes and failures, or a vector of successes
   with each row's trials held, a count above its row's trials refused; `$setData` takes new data whose counts
   give the successes and failures, the number of rows free. A multinomial sampler's `$setResponse` takes its
   n x K count matrix and refuses a vector, which with the trials held cannot say how they split. Alternative:
   a trials slot beside y, a second representation of the same counts.
10. At newdata with no weights a binomial fit draws one trial, as multinomial count fits draw one trial at new
    rows (dec-B177; a one-trial count vector under dec-B594), and a test sample's weights give its trials. rstanarm and brms require the
    trials at newdata (critique probe p9: rstanarm stops on the missing columns); dbarts's newdata is a predictor
    matrix with the response's columns absent, so its trials come through weights as on the other count families.
    Alternative: refuse ppd at newdata without weights, as rstanarm and brms in effect do, which no other dbarts
    count family does.
11. The pointwise log-likelihood of a binomial row is log dbinom(k; n, p), its binomial coefficient included as in
    glm's logLik, times the row's frequency weight. Alternative: the kernel alone, the two-frequency-row value,
    which differs by log choose(n, k) and is not the row's density.
12. xbart's losses on a logistic or probit fit, binomial rows included, are one rule on the probability scale: a
    row's observed proportion y (0 or 1 on a binary row, k / n on a binomial one) against its posterior mean
    probability p, the link applied per draw and averaged last, with the row weighted by its trials times its
    weight. The log loss is then -[y log p + (1 - y) log(1 - p)], the misclassification rate the share of trials
    on the wrong side of 0.5 and rmse the root of the weighted mean of (y - p)^2. The log loss and the
    misclassification rate equal the expansion into single trials; rmse on the expansion would add
    (k / n)(1 - k / n) per trial, a constant of the data that does not change which setting xbart picks, and the
    row rule is kept. A user loss on a binomial fit receives y.test as k / n, the latent draws as on a binary fit,
    and weights as the trials times the row's weight. Today's binary rmse is not this rule: it averages the
    latent draws, unlinked, against the 0/1 response ([`xbartLossFunction`](../../R/xbart.R)); run on
    c3d268c8f, the draws ranged over -2.9 to 3.8 and the rmse was 0.965 against a response sd of 0.501. That is
    a defect on main, filed apart from this note; the binomial slice reads the corrected rule. Alternative: one
    loss per row, which weights a row of one trial like a row of a thousand.
13. Refused by name: successes or failures negative, not whole, or missing outside na.action's handling;
    failures written so that k exceeds n is a negative failure count; trials above dec-B202's cap of one
    million, frequency weight included. A row of 0 trials is taken as dec-B133 takes a multinomial zero-trial
    row: out of the likelihood under the row mask, with its warning; its pointwise log-likelihood 0 in every
    draw, not a zero weight's NaN; its predictive draw 0, Binomial(0, p), as a zero-trial multinomial row draws
    zeros (dec-B177); its residual -p, its observed proportion read as 0, glm's response residual at a
    zero-weight row; and weight 0 in xbart's losses. Alternative for the zero row: refusing it, as nnet's
    multinom does and dec-B133 declined.
14. xbart takes a probit fit's 0/1 weights. Its refusal's reason was that each fold's fit is built inside the
    C loop with no channel to install a mask; under dec-B569 the weights ride the fold's data object and
    creation installs them as the mask, so a 0 row is out of each fold's fit and, at weight 0, out of its
    held-out loss. The refusal in [`xbart`](../../R/xbart.R) goes. Alternative: keep refusing, which leaves a
    membership vector that works everywhere else stopping at cross-validation.

## What moves, place by place

On a fit that chooses frequency, or a count family whose own reading is counts; a precision fit moves nowhere.

| where a weight is read | precision (default) | frequency |
|---|---|---|
| sigma's posterior | df nu + positive rows, SSR sum(w r^2) | df nu + sum(w), same SSR |
| sigest and the chi prior's anchor | lm.wfit, df positive rows - p | df sum(w) - p |
| starting sigma's stand-in | unweighted sd over positive rows | weighted sd |
| mean forest: leaf updates, tree moves, k | precision w | unchanged |
| variance forest | leaf count: positive rows | leaf count sum(w), in its leaf draw and tree moves |
| Student-t | composite w lambda | composite Lambda, shape w (nu + 1) / 2, nu step generalized |
| predictive draw | sd sigma / sqrt(w); 0 gives NA | sd sigma; 0 gives NA; predict takes weights 0 and 1 only |
| pointwise log-likelihood | log N(y; f, sigma^2 / w); 0 gives NaN | w log N(y; f, sigma^2); 0 gives NaN |
| fitted, residuals | unweighted | unchanged |
| xbart | loss weighted by w | loss unchanged; fold sigest frequency |
| response range; the default cut rule's grid | positive rows; every row | unchanged |
| leaf covariate standardization, forest default sd | unweighted | weighted by w |
| the quantile cut rule's grid | unique values | unchanged |
| state, copy, reload | no digest (dec-B586) | no digest; the destination's reading |
| logistic, hazard.logistic | | counts as before; 0 taken (dec-B570); training draw Bernoulli(p) (dec-B566) |
| nbinom, multinomial | 0/1 mask; all ones | whole counts, shape times w; training draw one observation |
| probit, ordinal | 0/1 mask, cleared from the data object | 0/1 kept on the data object (dec-B569) |
| hazard and hazard.probit, aft | as before | as before |
| hurdle.lognormal | 0 and 1 on both parts (dec-B588) | the same; "precision" refused |

Consumers, read without editing: stan4bart passes its weights to the BART part and reads them as precisions in
its own Stan target, so it is unaffected at the default; its predict(newdata, type = "ppd") ignores weights
(stan4bart TODO predict-ppd-weights), which is dec-B565's frequency draw and a gap under precision; its probit
rule passes 0/1 weights, which dec-B569 keeps. bartCause passes population weights to bart and stan4bart's
response fits; its sampleFromPPD draws at sigma whatever the weights, the same gap; population weights are
neither reading's case (dec-B568's caution). treatSens's BART path passes no weights to dbarts; bairrtt reads
none. None fits a binomial response through dbarts.

## Cost

Measured for TODO frequency weights on a loaded machine, 20000 rows, 50 trees: the weighted path costs about 14
percent over the unweighted one, an all-ones mask nothing, a masked row as much per sweep as an active one.
Frequency adds nothing per sweep on the gaussian-type families: sigma's count is fixed between weight
installs, a variance leaf sums w in the pass that counts rows now, and Student-t's extra nu term is a per-grid
constant. On the Polya-Gamma families a sweep grows with the total shape, as logistic's does now (figures above);
a binomial row costs as two frequency rows of weights k and n - k. No change moves a draw of a fit that runs
today except the count families' training predictive draws (dec-B566).

## Slices a plan would take

1. Engine, gaussian-type: the reading flag and per-row count, sigma's df, the variance leaf, Student-t; C++ tests against the expanded data and the w = 1 identity. About 250 lines, 250 of tests.
2. Engine, Polya-Gamma: nbinom's and multinomial's shapes times w, the shape grid and softmax margins under the
   power, the cap on the product. About 150 lines, 200 of tests.
3. R and bridge: weights.type, the model record, the refusals, sigest and the read-off quantities, the
   log-likelihood, the predictive draws of dec-B565 and dec-B566, help and NEWS. About 500 lines, 400 of tests.
4. Binomial: the cbind and two-column parse, family routing and the "auto" message, the counts slot and
   `$setCounts`, the logistic response taking a proportion with its count, the predictive draw, predict's
   trials, fitted and residuals, the log-likelihood, xbart's losses, the refusals. Gate: a binomial row of
   0 < k < n against its expansion into two adjacent rows of frequency weights k and n - k, draw for draw, since
   the two consume the same n unit draws at the same linear predictor (at k = 0 or n the expansion's weight-0
   row still draws one variate, so those rows are checked in distribution). About 400 lines, 300 of tests.
5. Gates: the gaussian frequency fit against the expanded rows draw for draw (probe p2's emulation as the
   seed); Student-t, heteroscedastic, nbinom and multinomial against the expanded rows in distribution; an
   exact-posterior arm for sigma on a single tree; an SBC arm. About 300 lines.
6. Header: the field in dec-B85's change and the getLatents comment. About 40 lines and the ABI hash.
7. dec-B569 and dec-B570: 0/1 weights on the data object, xbart taking them, logistic's 0. About 230 lines with
   tests.

In all about 3,000 lines, code and tests about evenly split, the binomial about 700 of them. The approximate-
draw research that dec-B567 opens is the TODO approximate-distributions arc, not a slice here.
