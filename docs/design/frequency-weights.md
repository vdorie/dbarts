# Frequency weights

Status: PROPOSED 2026-11-03, design note before any plan; amended the same day after its blind critique.
Rulings it builds on, in docs/decisions.md: dec-B125 (frequency weights join precision weights by explicit
choice; precision stays the default), dec-A11 (sigma's degrees of freedom count only positive-weight rows under
precision weights), dec-B120 (the comparison against 0.9-34 is a release gate), and the weight rulings named
where they apply below.

## Summary

A weight in dbarts today means one of three things, depending on the family. On a gaussian or Student-t fit it
is a precision: the row's noise variance is sigma squared over the weight, as in R's lm. On a logistic fit it is
a count of identical observations. On probit, ordinal and negative binomial fits only 0 and 1 are taken, and they
say which rows are in the data. The maintainer ruled that frequency weights, where a row stands for that many
identical observations, join precision weights as an explicit choice, with precision staying the default.

The main finding is that the three meanings are really two, and they differ in one place only. Both readings
raise a row's pull on the fit by its weight. They differ in how many observations the noise is estimated from:
under precision each row with a positive weight counts once, and under frequency it counts as many times as its
weight. On a fit whose noise level varies with the predictors, that count also reaches the trees that model the
noise. Families with no noise level to estimate (probit, logistic, ordinal, negative binomial, multinomial) have
one reading only. A logistic count is a frequency weight, and 0 and 1 mean the same thing under every reading:
out, or in once. So "membership" is not a third meaning; it is the part of the frequency reading those families
can draw today. R reads weights the same way: lm and glm's gaussian family take them as precisions; glm's
binomial and Poisson families, survreg, coxph, polr and multinom take them as counts; and MASS's rlm lets the
user say which, through an argument whose two values are the variance reading and the count reading.

A frequency-weighted fit should be the fit to the data with every row written out as many times as its weight,
and the sampler can produce exactly that. On the gaussian family a frequency fit and the fit to the written-out
data then draw the same numbers, sweep for sweep, which a probe confirmed by emulation. The families whose draws
need a hidden value per observation (probit, ordinal, negative binomial, multinomial, and censored survival
times) take whole-number weights by drawing one hidden value per copy, which is what logistic already does.

One caution belongs in the help: copies of a continuous response are exact ties, so a frequency fit reads them
as nearly noiseless data. In one probe the written-out data gave a noise sd of 0.05 where the truth was 0.4.

Seven choices need the maintainer:

1. Which families take frequency weights. Recommended: every family, as whole numbers wherever each copy needs
   its own hidden value, so probit, ordinal, negative binomial, multinomial and survival fits all take counts.
2. How the choice is written. Recommended: one argument beside the weights on the fitting functions, as rlm has
   one, read the same on every family, with the words "precision" and "frequency".
3. What a frequency-weighted row's predictive draw is. Recommended: one new observation like the row, as
   rstanarm and brms draw.
4. Whether logistic's predictive draw at a row of count w stays the number of successes among w trials.
   Recommended: it becomes one new observation too, so that every count reads alike; asking predict for w trials
   still gives the count.
5. Whether frequency weights must be whole numbers. Recommended: yes; real-valued weights keep the precision
   reading.
6. Whether probit and ordinal keep 0/1 weights or take membership under its own name. Recommended: keep the
   weights, and let the sampler's data object carry them as it does on a gaussian fit.
7. Zero weights. Recommended: as now everywhere, and the logistic family accepts a weight of 0 as every other
   family does.

## What a weight means today

| family | values taken | read as | a 0 |
|---|---|---|---|
| gaussian, Student-t, heteroscedastic | finite, non-negative | precision: noise variance sigma^2 / w | out of the likelihood, in the design |
| logistic | positive integers up to 1e6 | copies: Polya-Gamma shape w | refused, after a warning that it stays in the fit |
| probit, ordinal, nbinom, hazard | 0 and 1 | the active-row mask | out of the likelihood, in the design |
| multinomial | all ones | refused otherwise | refused |
| aft, hurdle.lognormal | none | refused | refused |

The routing is [`enforceWeightPolicy`](../../R/spec.R), with [`isMaskedWeightFamily`](../../R/spec.R) naming the
mask families and [`refuseNonCountWeights`](../../R/spec.R) the logistic rule. Run on the main tip (scratch
library built from c3d268c8f): the refusals fire as stated; a logistic fit's training predictive draw at a row of
weight 5 ranges over 0 to 5, Binomial(w, p), and its pointwise log-likelihood is w times the Bernoulli log
density, the count rule of dec-A116; a probit fit refuses weights of 5 while its predict(type = "ppd", weights =
5) draws 0 to 5, the trial counts of dec-B470; logistic weights mixing 0 and 2 raise the zero-weight warning and
then the refusal. rbart_vi, deprecated and removed at 1.1-0, keeps precision weights and takes none of what
follows.

0.9-34's probit weights were neither reading: its latent draw was a truncated normal of sd 1/sqrt(w) about the
fit, so the model it fitted was P(y = 1) = Phi(sqrt(w) f(x)), a rescaled link per row and not a weighted
likelihood (classic engine, bartFit.cpp at the 0.9-34 tag). dec-B13 removed it: "The weighted probit path is gone
rather than carried into the new engine, because scaling the latent draws by one over the square root of the
weight is not a coherent weighted likelihood. Weights on a logistic binary fit are supported under the common
reading of a weight as a number of identical observations, restricted to positive integers. A 0.9-x weighted
probit call now errors, and the only migration is an integer-count logistic fit." It ruled out the latent
rescaling and put counts on logistic; it did not rule on counts on probit. The maintainer accepted the 0.9-x
error, so nothing below argues from it as a regression.

How others read the same argument, from installed packages and their sources, and run where marked:

| fitter | reading | non-integer | run |
|---|---|---|---|
| lm, glm gaussian, mgcv gam, lme4 lmer | precision; glm: "each response y_i is the mean of w_i unit-weight observations" | taken | lm, glm, gam: sigma differs from the expanded-data fit |
| MASS rlm | wt.method = c("inv.var", "case"): "a weight of 2 means there are two of these" or "this error is half as variable" | taken | "case" equals the expanded-data fit, "inv.var" does not |
| glm binomial, Poisson | copies (binomial: successes w y of w trials) | warns "non-integer #successes" on binomial | |
| survival survreg (lognormal AFT included), coxph | case weights, "treated as replication weights" | taken; coxph then defaults to a robust variance, treating them as sampling weights | scale and standard errors equal the expanded-data fit |
| MASS polr, nnet multinom | case weights | polr warns from its starting glm | polr: standard errors equal the expanded-data fit |
| rstanarm 2.32.2 | log_lik multiplies by the weight; posterior_predict draws one observation per row on every family | taken | a gaussian identity stan_glm with no offset ignores the weights entirely: its draws equal the unweighted fit's |
| brms | `weights()` multiplies the log-likelihood; `posterior_predict` ignores it; `trials()` is a separate term | positive numeric values | not run |
| survey | design weights, rescaled to mean 1 inside svyglm; inference from the design | taken | |
| BART (wbart) | `w` multiplies the noise sd | taken | |
| bartMachine | no weights | | |
| Stata, SAS, statsmodels | two kinds by name: fweight and aweight, FREQ and WEIGHT, freq_weights and var_weights | Stata refuses non-integer fweights; SAS truncates FREQ | not run, from their documentation |

## The two readings differ only in the scale's count

Write a row's likelihood contribution in glm's form, its kernel times a normalizer that may depend on a scale.
A precision weight multiplies the kernel's exponent: for a gaussian, exp(-w r^2 / 2 sigma^2) with one factor of
1/sigma per row of positive weight, which is dec-A11's count, "the likelihood's own". A frequency weight raises
the whole contribution to the power w: the same kernel and w factors of 1/sigma. So:

- On a family with no estimated scale the two readings are one likelihood. Logistic's counts are frequency
  weights already (the Polya-Gamma shape w is exactly w copies, [weighted-logistic.md](weighted-logistic.md)), a
  binomial row of k successes in n trials is two frequency rows, y = 1 at weight k and y = 0 at weight n - k, and
  0 and 1 are both readings on every family.
- On a family with a scale, the mean forest's leaf updates and moves are unchanged given the row precisions;
  the residual scale's posterior counts sum(w) observations instead of the positive rows. On a heteroscedastic
  fit the variance forest's leaves count the same way, and since that count enters the variance leaf's
  integrated likelihood ([`ConstantVarianceLeaf`](../../src/bartcore/model.hpp)), the variance trees' moves
  change too.
- A precision weight on a latent family is not a weighted likelihood at all (the 0.9-x probit above). Away from
  0 and 1, those families have only the frequency reading.

The target a frequency-weighted fit must hit is the fit on the expanded data, every row written out w times.
Nothing in the tree machinery separates the two: the empty-leaf rule counts rows reaching a leaf (one copy
reaches it if any does), the default cut rule places one point per gap between distinct values, and the moves
choose among variables and cuts by the values a node holds. What does separate them is every quantity the
package reads off the rows: the starting sigma, the quantile cut rule, a leaf covariate's center and scale, a
forest's default sd. Those are read off the expanded rows under frequency (Calls made, 3). Emulated on the main
tip, a precision fit whose chisq prior carries the frequency degrees of freedom and whose sigest is the
expanded data's matches the expanded-data fit draw for draw: sigma to 2.2e-16 and the fit to 2.9e-15 over 300
sweeps, 40 rows of weights 1 to 4 (critique probe p2, rerun). A gaussian sweep's random numbers do not depend on
how many times a row appears, so the gate is bitwise in all but summation order.

The same probe shows what frequency means for a continuous response: the expanded data's sigma settled at 0.05
with the truth at 0.4 (the precision fit at 0.60), because w copies of y at one x are exact ties and the trees
can fit each distinct x. That is the frequency model, not a defect; the help says it (Calls made, 8).

A zero weight is the one place the target is not the expanded data: a row of weight 0 stays in the design, as
it does under precision weights and the mask (dec-B238, dec-B251), so a frequency fit with zeros differs from
the expanded data as a precision fit with zeros differs from the subsetted one.

## Family by family

What w copies of a row mean for the likelihood, the latent draws and the residual prior, and whether the engine
draws it exactly:

| family | w copies | latent | residual prior | exact draw |
|---|---|---|---|---|
| gaussian | N(f, sigma^2)^w: leaf statistics as precision w | none | df nu + sum(w), SSR sum(w r^2) | yes |
| heteroscedastic | N(f, s^2(x))^w | none | a variance leaf's count is sum(w) over its rows, in its draw and its integrated likelihood | yes |
| Student-t | t_nu^w: the w copies' mixing scales sum to one Gamma(w (nu + 1) / 2, (nu + r^2 / sigma^2) / 2) | one composite scale per row | df nu + sum(w); the nu step reads sum(w log Lambda) and a per-grid constant of the weights | yes |
| logistic | Bernoulli^w | PG(w, psi), a sum of w PG(1) draws | none | yes, shipped |
| probit, hazard | Phi^w | the mean of w truncated normals | none | yes, w draws per row per sweep |
| ordinal | the category probability^w | as probit, on the category's interval; the cutpoint step reads the likelihood and takes the power | none | yes, w draws |
| nbinom | NB^w | PG(w (y + r), psi) | the shape grid reads the likelihood and takes the power | yes, w (y + r) draws |
| multinomial | each row's softmax likelihood^w | PG(w n) per category | none | yes, w n draws |
| aft | events N^w; censored S(c)^w | w censored log-times; sigma also reads their within-row sum of squares | df nu + sum(w) | yes, w draws on censored rows |
| hurdle.lognormal | both parts^w | the zero part's | the positive part's | follows its two parts; its weights are owed by feature completion (dec-A15) |

Every row is exact for whole-number w. The gaussian, heteroscedastic and Student-t rows are exact for real w as
well. No exact draw of a censored or binary row's w latents as one aggregate was found, so those families pay
one draw per copy, as logistic does (dec-B202 measured 160 ns per trial and capped a row at one million).

Student-t is the one family with a derivation to state. For w copies, the copies' mixing scales lambda_ij are
iid Gamma((nu + 1) / 2, (nu + r^2 / sigma^2) / 2) given the residual, and the gaussian part reads only their sum
Lambda_i = sum_j lambda_ij, which is Gamma(w (nu + 1) / 2, same rate). The trees and sigma given nu are then
exact, sigma's posterior reading sum(w) observations and the sum of Lambda r^2. With m_i the row's count (w_i
under frequency, 1 under precision, rows out of the likelihood skipped), the nu step on its grid reads sum(m_i)
where it reads the row count, sum(m_i log Lambda_i) and sum(Lambda_i) where it reads the lambda statistics, plus
sum_i [m_i lgamma((nu + 1) / 2) - lgamma(m_i (nu + 1) / 2)], which depends on the weights alone and is computed
per grid point when they are installed. At m_i = 1 that term is zero and the draw is the shipped one
([`ResidualDfPrior::drawIndex`](../../src/bartcore/model.hpp)). The critique's probe p8 ran this step on fixed residuals, at
integer and at real counts, against the exact grid probabilities of nu: within 0.0026 (Monte Carlo sd about
0.006); without the lgamma term a grid probability is off by up to 0.963.

## Open calls

### 1. Which families take frequency weights

Background: dec-B125 asks for frequency weights "for models where it makes sense". Every family above is a
coherent frequency model with an exact draw. A user with aggregated records means the same thing on any family,
and R's own fitters read it that way on binary, ordinal, multinomial and survival models (glm, polr, multinom,
survreg, coxph). Today the answer depends on the family: logistic takes counts, probit refuses them at the fit
and takes them at predict, ordinal and nbinom refuse them, aft takes no weights. The TODO stamp on
weighted-binary, "decision-gated, scheduled after 1.0-0 (VD 2026-08-28)", records the ruling of 2026-08-28,
"scheduled after 1.0-0, shares the decision above", where the decision above was the real-shape Polya-Gamma
draw for real-valued dispersion; that door's undecided part was real weights on the binary families. It
predates dec-B125, and it did not rule on whole-number counts.

- (a) Frequency weights on every family: real or whole (open call 5) on the gaussian, Student-t and
  heteroscedastic families and everything built on them, whole numbers on probit, hazard, ordinal, nbinom,
  multinomial and aft. Since the two readings coincide where there is no scale, a probit, ordinal, nbinom or
  multinomial fit takes counts under either value of the setting; aft takes both readings, precision drawing
  its censored latent at sd sigma / sqrt(w). One rule: a weight counts a row w times in the likelihood, and on a
  family with a residual scale it is a precision instead unless frequency is chosen. Cost: one latent draw per
  copy per sweep on those families, as logistic pays; about 200 to 300 lines per family with its tests and an
  exact-posterior arm, aft also handing sigma a within-row sum of squares; the count cap of dec-B202 on each.
- (b) The choice on the gaussian, Student-t and heteroscedastic families only; the others keep what they take
  today. The families then read one user statement four ways, and probit's fit and predict disagree on what a
  weight of 5 is. Cost: slices 1 to 4 below.
- (c) Gaussian alone. Student-t and the variance forest, which share the gaussian's residual prior home
  (dec-B123), refuse a choice their model takes. Cost: less than (b).

Recommendation: (a).

### 2. How the choice is written

Background: the reading is part of the model, the likelihood, not of the data; a sampler holding weights reads
them under the model it was created with. The one R fitter that offers both readings, MASS's rlm, puts the
choice in an argument of the fitting function beside the weights: wt.method = c("inv.var", "case"). Nothing in R
puts it on a family. Because the readings coincide on every family without a scale, one spelling written at the
top reads the same on every family, which a family setting cannot: logistic and probit would have no place to
write it. The words in use: "case weights" in survival, MASS and nnet ("treated as replication weights" in
coxph's help), "prior weights" in glm and lme4, "frequency weights" in dec-B125, this package's documentation,
Stata (fweight), SAS (FREQ) and statsmodels (freq_weights), against "inverse-variance", "analytic" or
"precision" for the other.

- (a) An argument of the fitting functions, weights.type = c("precision", "frequency"), on bart (formerly
  bart2), dbarts, dbartsSpec and xbart; the value is held on the model record and fixed at creation. The
  sampler's setWeights, setData and predict read the vector under it; dbartsData keeps one vector; the shipped
  header's specification gains one enumeration field on its family part in dec-B85's change, so consumers
  rebuild once. bartBT, the 0.9-x door, and the deprecated rbart_vi keep precision with no argument. Cost: a
  formal on four entry points and in a consumer that forwards it; on a family without a scale both values give
  one model, which the help says.
- (b) The same argument with rlm's words, wt.method = c("inv.var", "case"). Cost: the package's own help, the
  ruling and the other languages say precision and frequency; "case weights" also names sampling weights in
  survey practice.
- (c) The same argument named weights.type with the values "precision" and "case". Cost: as (b) for "case".
- (d) A family setting, gaussian(weights = "frequency"), beside the residual prior. Cost: no place to write it
  on the families without a scale, so the statement differs by family; no R precedent.
- (e) A marker on the value, weights = frequency(w), or a second vector, frequency.weights. Cost: the marker
  carries the model on the data, so setWeights could change it, and model.frame drops the class; the second
  vector doubles the slots and setters for a rule, both at once, that no user case asks for.

Recommendation: (a).

### 3. A frequency-weighted row's predictive draw

Background: under precision weights a row's draw is an observation of variance sigma^2 / w (and 0 gives NA by
every route, dec-B472). Under frequency a row stands for w observations, each of sd sigma. glm's help calls
"the mean of w_i unit-weight observations" the precision reading, so drawing that mean for a frequency row
gives it the precision row's draw. rstanarm and brms draw one observation per row whatever its weight
(rstanarm run: sd sigma at weights 1 and 5; brms read from source).

- (a) One new observation like the row: sd sigma. predict at a training row then agrees with predict at newdata
  holding the same x, and the draw is on the scale of y, so a check of y against its draws compares like with
  like. A weight given to predict on a frequency fit has no use for this draw and is refused by name. The
  log-likelihood stays w times the one-observation value (dec-A116), which is the row's whole contribution to
  the likelihood, not the density of what is drawn. Cost: a predictive draw that ignores the training weight,
  documented.
- (b) The mean of the row's w observations, sd sigma / sqrt(w), the precision formula at frequency's sigma.
  Totals over a weighted sum of rows come out with the right spread without the user reweighting. Cost: a
  frequency row gets the precision row's draw, and predict at a training row and at the same x as newdata
  disagree.

Recommendation: (a).

### 4. Logistic's predictive draw at a row of count w

Background: a logistic row of count w draws Binomial(w, p) today, the successes among its w copies (run: 0 to 5
at w = 5), and dec-B470 ruled that weights given to predict on a probit or logistic fit are "positive integer
trial counts". A count is a frequency weight, so under open call 3 (a) the gaussian frequency row draws one
observation while the logistic count row draws the total of its copies. Logistic is not in 0.9-34.

- (a) The training-row draw becomes one observation, Bernoulli(p), as on every frequency row and as rstanarm
  draws; predict(weights = w) keeps dec-B470's trial counts, an explicit request for the successes among w new
  trials, and with open call 1 (a) the same holds on probit and the other count families. Cost: extract and
  fitted(type = "ppd") on a weighted logistic fit move, and the help states the two cases.
- (b) Keep Binomial(w, p) at training rows. The count families draw a row's total and the scale families one
  observation, two conventions for one reading.

Recommendation: (a).

### 5. Whether frequency weights must be whole numbers

Background: dec-B125 left "the treatment of non-integer values" to this note. A frequency weight says how many
identical observations a row stands for; a weight of 2.5 has no such meaning and makes the likelihood a power,
which is what Stata calls an importance weight and refuses as an fweight. dec-B103 ("weights may be fractional
wherever the family supports real weights") separated counts from weights in general; it did not rule on this
reading. The gaussian-type draws are exact for real w; the latent families' are not. The use that brings real
weights is survey or inverse-probability weighting, where a likelihood weight gives the right point estimate and
uncertainty sized to sum(w) rather than to the design ([weighted-logistic.md](weighted-logistic.md)); coxph,
which takes case weights as replication, switches to a sampling-weight variance as soon as one is not whole.

- (a) Whole numbers only under frequency, on every family; real values keep the precision reading. One value
  domain for the reading across families, and the refusal says the weights are counts. Cost: a user with
  fractional counts (from a model, or rescaled) rounds them or uses precision.
- (b) Real values taken where the draw is exact (gaussian, Student-t, heteroscedastic), whole numbers elsewhere,
  as brms, rstanarm and survreg take them. Cost: the domain differs by family, and normalized weights give a
  sigma posterior from sum(w) observations, one observation's worth if they sum to 1.

Recommendation: (a).

### 6. Membership on probit and ordinal: weights, or its own name

Background (the maintainer's question of 2026-10-05): 0/1 weights on probit, ordinal and nbinom install the
active-row mask and are cleared from the data object, so `s$data` carries no weights and `s$setData(s$data)`
drops the rows out (TODO mask-on-own-data-object); xbart refuses 0/1 weights on those families, its folds having
no mask channel. Since 0 and 1 mean out and in once under every reading, membership is the part of the frequency
reading those families can draw now, and counts on them under open call 1 (a) extend it.

What the one argument costs today: one paragraph per family in every weights item (bart, dbarts, the sampler,
stan4bart's), the policy branches, the data object that drops the mask, and the xbart refusal. With the finding
above the weights help becomes one rule plus which values each family takes.

- (a) Keep 0/1 weights, and let them ride the data object as weights on every family, the engine installing
  them as the mask; `setActiveRows` stays the door of a sampler embedded in a larger sampler and composes with
  the weights. Cost: `s$data@weights` is no longer NULL on a mask family; resolves mask-on-own-data-object;
  about 150 lines with tests.
- (b) As now.
- (c) Membership under its own name on the fitting functions (such as active = ), weights refused on probit and
  ordinal. Cost: a new formal on bart, dbarts, dbartsData and xbart; stan4bart, whose probit rule passes 0/1
  weights, and any caller doing the same break; the name means what 0/1 weights mean anyway.

Recommendation: (a).

### 7. Zero weights

Background: a weight of 0 takes a row out of the likelihood and keeps it in the design (dec-B238, dec-B251), out
of sigma's count (dec-A11), NaN in the pointwise log-likelihood (dec-A116), NA in the predictive draws
(dec-B472). Under frequency a 0 is no copies, the same statement. Logistic alone refuses a 0, after warning that
such rows stay in the fit, though its engine takes the mask
([`LogisticResponse::workingWeights`](../../src/bartcore/model.hpp)).

- (a) As now on every family and under both readings, and logistic accepts 0 as the mask. Cost: about 80 lines
  with tests.
- (b) As now, logistic still refusing; its warning is then dropped where the refusal follows.

Recommendation: (a).

### Settled by rulings, not open

- pdbart averages unweighted unless given average.weights (dec-B214); the help says that average.weights = w
  gives the expanded data's average.
- Precision stays the default (dec-B125), so no existing fit and no scenario of dec-B120's comparison moves.

## Calls made

1. The engine carries one reading flag per sampler and derives a per-row count from it, 1 for a positive
   weight under precision and w under frequency, used by sigma's degrees of freedom
   ([`ChiSquaredScalePrior::drawSigmaSqFromPosterior`](../../src/bartcore/model.hpp), which takes the positive
   count today), a variance leaf's count in its draw and integrated likelihood, and Student-t's shape and nu
   step. Alternative: a second per-row vector of counts beside the precisions, which would allow both kinds at
   once as statsmodels does; no user case and n more doubles.
2. Student-t's nu step is the shipped grid draw generalized as derived above, identical to the shipped draw
   under precision. Alternative: the marginal-t grid draw, nu given sigma and the fit with the scales integrated
   out, also exact and better mixing, but a second mechanism and a changed draw at w = 1.
3. Under frequency every quantity the package reads off the rows is read off the expanded rows: the starting
   sigma's regression with residual df sum(w) - p ([`residualStandardError`](../../R/utility.R) reads lm.wfit's
   positive-row df today; at sum(w) - p it equals lm on the expanded rows), its stand-in a weighted sd, a linear
   or gp leaf covariate's center and scale, a forest's default sd, the quantile cut rule's quantiles and nbinom's
   center. Where sum(w) - p is not positive the stand-in is used, as where the regression is undefined today.
   The data object's sigma, the starting estimate, is derived under the model's reading wherever it is derived,
   and a stated one is taken as stated. Alternative: read them as precision does, over the positive rows
   unweighted; the fit would then differ from the expanded-data fit through its prior's anchors, and the
   draw-for-draw gate would not hold.
4. The reading is model: held on the model record, fixed at creation as the family is, and a setModel that
   changes it is refused. A fit object saved before it reads precision. Alternative: on the data object, where
   setWeights or setData could change the model.
5. A saved state's weights digest covers the reading, so a state installed into a sampler of the other reading
   redraws Student-t's scales as a weight change does
   ([`TResponse::reapplyWeights`](../../src/bartcore/model.hpp)), and the count families' latents as their
   weight change does; gaussian and heteroscedastic states carry no latent and install as they are.
   Alternative: refuse such an install, which no other weight change does.
6. A row of frequency weight 0 reports NaN in the pointwise log-likelihood, as a precision 0 and a masked row
   do ([`pointwiseLogLikelihood`](../../R/generics.R)). Alternative: 0, the product, which a tool such as loo
   would count as an observation.
7. xbart keeps a row's copies in one fold, weights the held-out loss by w as it does now
   ([`xbartLossFunction`](../../R/xbart.R)), and takes each fold's starting sigma under the fit's reading.
   Alternative: splitting a row's copies across folds, which leaks identical observations into training.
8. The help for frequency weights states: copies of a continuous response are exact ties, read by the model as
   nearly noiseless data where the trees can separate the distinct rows; weights describing a population or a
   design are not counts; under Student-t, the latents the sampler reports
   ([`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) and the R method) are the row's composite
   scale, the sum over its copies; and a binomial row of k successes in n trials is two rows. Alternative:
   documenting only the arithmetic.

## What moves, place by place

Under the recommendations, on a fit that chooses frequency; a precision fit moves nowhere except where open
calls 1 (a), 4, 6 and 7 change the count families.

| where a weight is read | precision (default, unchanged) | frequency |
|---|---|---|
| sigma's posterior | df nu + positive rows, SSR sum(w r^2) | df nu + sum(w), same SSR |
| sigest and the chi prior's anchor | lm.wfit, df positive rows - p | df sum(w) - p; the prior's scale moves with it |
| starting sigma's stand-in | unweighted sd over positive rows (TODO starting-sigma-fallback-zero-weights) | weighted sd |
| mean forest: leaf updates, tree moves, k | precision w | unchanged |
| variance forest | leaf count: positive rows | leaf count: sum(w), in its leaf draw and its tree moves; its leaf prior follows sigest |
| Student-t | composite w lambda, nu from lambda statistics | composite Lambda, shape w (nu + 1) / 2, nu step generalized |
| predictive draw | sd sigma / sqrt(w); 0 gives NA | sd sigma, one observation; 0 gives NA; predict's weights refused |
| pointwise log-likelihood | log N(y; f, sigma^2 / w); 0 gives NaN | w log N(y; f, sigma^2); 0 gives NaN |
| fitted, residuals | unweighted | unchanged |
| xbart | loss weighted by w; fold sigest precision | loss unchanged; fold sigest frequency |
| response range; the default cut rule's grid | range over positive rows; grid over every row | unchanged: neither depends on how many times a row appears |
| leaf covariate standardization, forest default sd, the quantile cut rule | unweighted | weighted by w |
| state, copy, reload | digest of weights | digest of weights and reading |
| logistic | counts; training draw Binomial(w, p) | the same model under either value; training draw Bernoulli(p) (open call 4), predict's weights trial counts |
| probit, hazard, ordinal, nbinom, multinomial, aft | 0/1 mask, or no weights | counts under either value (open call 1 (a)); aft precision or frequency |

Under the alternatives: 1(b) leaves the count families as they are and 1(c) also Student-t and the variance
forest. 2(b) and 2(c) change the words only; 2(d) and 2(e) move the same rows through a family or the data.
3(b) draws sd sigma / sqrt(w) and takes predict's weights as a count of observations whose mean is drawn. 4(b)
keeps logistic's Binomial(w, p). 5(b) admits real weights on the three gaussian-type families. 6(b) and 7(b)
move nothing; 6(c) moves every probit and ordinal call that passes 0/1 weights, stan4bart's included.

Consumers, read without editing: stan4bart passes its weights to the BART part and reads them as precisions in
its own Stan target ("the weighted normalizer's N-not-sum(w)" in parametric_model.hpp), so it is unaffected at the
default and would need its own normalizer to offer frequency; its predict(newdata, type = "ppd") ignores weights
(stan4bart TODO predict-ppd-weights), drawing one observation: open call 3 (a)'s frequency draw, and a gap under
precision. bartCause passes population weights to bart and stan4bart's response fits and uses them in its
estimators; its sampleFromPPD draws at sigma whatever the weights, the same gap. Population weights are neither
reading's case (open call 5), so neither consumer should move to frequency by default. treatSens's BART path
passes no weights to dbarts; bairrtt reads none.

## Cost

Measured for TODO frequency weights on a loaded machine, 20000 rows, 50 trees: the weighted path costs about 14
percent over the unweighted one, an all-ones mask nothing, and a masked row as much per sweep as an active one.
Frequency adds nothing per sweep on the gaussian-type families: sigma's count is fixed between weight installs,
a variance leaf sums w in the pass that counts rows now, and Student-t's extra nu term is a per-grid constant
computed at install. The count families of open call 1 (a) draw one latent per copy per sweep, as logistic does.

Engine changes beyond the gaussian df: the variance leaf's count, Student-t's shape, composite and nu step,
the state digest, and under open call 1 (a) the per-copy latent draws of probit, ordinal, nbinom, multinomial
and aft. None changes a posterior draw of a fit that runs today; open call 4 (a) changes logistic's
predictive draw at rows of count above 1.

## Slices a plan would take

1. Engine, gaussian-type families: the reading flag and per-row count, sigma's df, the variance leaf, Student-t,
   the digest; C++ tests against the expanded data and the w = 1 identity. About 250 lines, 250 of tests.
2. R and bridge: weights.type, the model record, the refusals, sigest and the read-off quantities, the
   log-likelihood and predictive draw, help and NEWS; the flag through the bridge. About 450 lines, 350 of tests.
3. Gates: the gaussian frequency fit against the expanded rows draw for draw (critique probe p2's emulation as the
   seed), Student-t and heteroscedastic in distribution, an exact-posterior arm for sigma on a single tree, an SBC
   arm. About 250 lines.
4. Header: the enumeration field, in dec-B85's change, and the getLatents comment. About 40 lines and the ABI
   hash.
5. Count families under open call 1 (a): probit and hazard, ordinal, nbinom, multinomial, aft, each its per-copy
   draw, exact-posterior arm and replication check. About 200 to 300 lines each, 1,000 to 1,500 in all.
6. Open calls 4, 6 and 7: logistic's training draw, 0/1 weights on the data object, logistic's 0. About 300 lines
   with tests.
