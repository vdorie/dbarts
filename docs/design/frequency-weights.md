# Frequency weights

Status: PROPOSED 2026-11-03, design note before any plan. Rulings it builds on, in docs/decisions.md:
dec-B125 (frequency weights join precision weights by explicit choice; precision stays the default), dec-A11
(sigma's degrees of freedom count only positive-weight rows under precision weights), dec-B120 (the comparison
against 0.9-34 is a release gate), and the weight rulings named where they apply below.

## Summary

A weight in dbarts today means one of three things, depending on the family. On a gaussian or Student-t fit it
is a precision: the row's noise variance is sigma squared over the weight, as in R's lm. On a logistic fit it is
a count of identical observations. On probit, ordinal and negative binomial fits only 0 and 1 are taken, and they
say which rows are in the data. The maintainer ruled that frequency weights, where a row stands for that many
identical observations, join precision weights as an explicit choice, with precision staying the default.

The main finding is that the three meanings are really two, and they differ in one place only. Both readings
raise a row's contribution to the fit by its weight. They differ only in how many observations the residual
scale is estimated from: under precision each row with a positive weight counts once, and under frequency it
counts as many times as its weight. Families with no residual scale to estimate (probit, logistic, ordinal,
negative binomial, multinomial) have one reading only, the frequency one. A logistic count is a frequency weight,
and 0 and 1 mean the same thing under every reading: out, or in once. So "membership" is not a third meaning;
it is all that the latent families can draw today. R itself reads weights this way: lm and glm's gaussian family
take them as precisions, while glm's binomial and Poisson families, survreg, coxph, polr and multinom take them
as counts. rstanarm and brms take them as counts on every family.

What a frequency-weighted fit should be is clear: the fit to the data with every row written out as many times
as its weight. On the gaussian, Student-t and heteroscedastic families, which are where the choice makes a
difference, the sampler can draw that exactly for any positive weight, whole or not, at no extra cost per sweep.
The families whose draws need a latent value per observation (probit, ordinal, negative binomial, multinomial,
and censored survival times) can be drawn exactly only for whole-number weights, by drawing a latent for each
copy, at a cost that grows with the total weight. Integer-weight probit is already scheduled after 1.0-0.

Five choices need the maintainer:

1. Which families take the choice in 1.0-0. Recommended: the gaussian, Student-t and heteroscedastic families,
   with every composition built on them; the rest keep what they take today, and the whole-number cases join
   after 1.0-0 with the integer-weight probit work.
2. Where the choice is written. Recommended: on the family, as gaussian(weights = "frequency"), beside the
   residual prior that the choice changes, so the fitting functions, the data object, the sampler and the
   shipped C header gain no new argument.
3. What a frequency-weighted row's predictive draw is. Recommended: what the row records, as every family
   does today, which for a gaussian row of weight w is the mean of its w observations; a single observation
   is one call away.
4. Whether probit and ordinal keep 0/1 weights or take membership under its own name. Recommended: keep the
   weights, and let the sampler's data object carry them as it does on a gaussian fit.
5. Zero weights. Recommended: as now everywhere, and the logistic family accepts a weight of 0 as every other
   family does.

Whether non-integer frequency weights are taken is settled by an earlier ruling: wherever the family's draw
is exact for them, which is the three recommended families.

## What a weight means today

| family | values taken | read as | a 0 |
|---|---|---|---|
| gaussian, Student-t, heteroscedastic | finite, non-negative | precision: noise variance sigma^2 / w | out of the likelihood, in the design |
| logistic | positive integers up to 1e6 | copies: Polya-Gamma shape w | refused |
| probit, ordinal, nbinom, hazard | 0 and 1 | the active-row mask | out of the likelihood, in the design |
| multinomial | all ones | refused otherwise | refused |
| aft, hurdle.lognormal | none | refused | refused |

The routing is [`enforceWeightPolicy`](../../R/spec.R), with [`isMaskedWeightFamily`](../../R/spec.R) naming the
mask families and [`refuseNonCountWeights`](../../R/spec.R) the logistic rule. Run on the main tip (scratch
library built from c3d268c8f): the table's refusals fire as stated, a logistic fit's training predictive draw at
a row of weight 5 ranges over 0 to 5, Binomial(w, p), and its pointwise log-likelihood is w times the Bernoulli
log density, the count rule of dec-A116.

How others read the same argument, read from installed packages and their sources, and run where marked:

| fitter | reading | non-integer | run |
|---|---|---|---|
| lm, glm gaussian, mgcv gam, lme4 lmer | precision; glm: "each response y_i is the mean of w_i unit-weight observations" | taken | lm, glm, gam: sigma differs from the expanded-data fit |
| glm binomial, Poisson | copies (binomial: successes w y of w trials) | warns "non-integer #successes" on binomial | |
| survival survreg (lognormal AFT included), coxph | case weights "treated as replication weights" | taken; coxph defaults to a robust variance | scale and standard errors equal the expanded-data fit |
| MASS polr, nnet multinom | case weights | polr warns from its starting glm | polr: standard errors equal the expanded-data fit |
| rstanarm | the log-likelihood times the weight on every family, though its help says "Same as glm" (continuous.stan's weighted branch, `.weighted`) | taken | |
| brms | `weights()` multiplies the log-likelihood; `posterior_predict` ignores it; `trials()` is a separate term | positive numeric values | |
| survey | design weights, rescaled to mean 1 inside svyglm; inference from the design | taken | |
| BART (wbart) | `w` multiplies the noise sd | taken | |
| bartMachine | no weights | | |
| Stata, SAS, statsmodels | two kinds by name: fweight and aweight, FREQ and WEIGHT, freq_weights and var_weights | Stata refuses non-integer fweights; SAS truncates FREQ | not run |

## The two readings differ only in the scale's count

Write a row's likelihood contribution in glm's form, its kernel times a normalizer that may depend on a scale.
A precision weight multiplies the kernel's exponent: for a gaussian, exp(-w r^2 / 2 sigma^2) with one factor of
1/sigma per row of positive weight, which is dec-A11's count, "the likelihood's own". A frequency weight raises
the whole contribution to the power w: the same kernel and w factors of 1/sigma. So:

- On a family with no estimated scale, the two readings are the same likelihood. Logistic's counts are
  frequency weights already (the Polya-Gamma shape w is exactly w copies,
  [weighted-logistic.md](weighted-logistic.md)), and 0 and 1 are both readings on every family.
- On a family with a scale, the leaf updates, the trees and every tree move are unchanged; only the residual
  scale's posterior counts differently: sum(w) observations instead of the positive rows.
- A precision weight on a latent family (0.9-x's probit latent scaled by 1/sqrt(w)) is not a weighted likelihood
  at all, which is why dec-B13 removed it. Away from 0 and 1, those families have only the frequency reading.

The target a frequency-weighted fit must hit is the fit on the expanded data, every row written out w times.
Nothing in the tree machinery separates the two: the empty-leaf rule counts rows reaching a leaf (one copy
reaches it if any does), the default cut rule places one point per gap between distinct values, and the moves
choose among variables and cuts by the values a node holds. What does separate them is every quantity the
package reads off the rows: the starting sigma, the quantile cut rule, a leaf covariate's center and scale, a
forest's default sd. Those are read off the expanded rows under frequency (Calls made, 3).

A zero weight is the one place the target is not the expanded data: a row of weight 0 stays in the design, as
it does under precision weights and the mask (dec-B238, dec-B251), so a frequency fit with zeros differs from
the expanded data as a precision fit with zeros differs from the subsetted one.

## Family by family

What w copies of a row mean for the likelihood, the latent draws and the residual prior, and whether the engine
draws it exactly:

| family | w copies | latent | residual prior | exact draw |
|---|---|---|---|---|
| gaussian | N(f, sigma^2)^w: leaf statistics as precision w | none | df nu + sum(w), SSR sum(w r^2) | yes, any real w |
| heteroscedastic | N(f, s^2(x))^w | none | a variance leaf's count is sum(w) over its rows, not its positive rows | yes, any real w |
| Student-t | t_nu^w: the w copies' mixing scales sum to one Gamma(w (nu + 1) / 2, (nu + r^2 / sigma^2) / 2) | one composite scale per row | df nu + sum(w); the nu draw takes sum(w log Lambda) and a per-grid constant of the weights | yes, any real w (below) |
| logistic | Bernoulli^w | PG(w, psi), a sum of w PG(1) draws | none | yes, integers; shipped |
| probit, hazard | Phi^w | the mean of w truncated normals | none | integers, O(w) draws per row per sweep |
| ordinal | the category probability^w | as probit, on the category's interval; the cutpoint step reads the likelihood and takes the power | none | integers, O(w) |
| nbinom | NB^w | PG(w (y + r), psi) | the shape grid reads the likelihood and takes the power | integers, O(w (y + r)) |
| multinomial | each row's softmax likelihood^w | PG(w n) per category | none | integers, O(w n) |
| aft | events N^w; censored S(c)^w | w censored log-times; sigma needs their within-row sum of squares | df nu + sum(w) | integers, O(w) on censored rows |
| hurdle.lognormal | both parts^w | the zero part's | the positive part's | follows its two parts; takes no weights today |

Student-t is the one family with a derivation to state. For w copies, the copies' mixing scales lambda_ij are
iid Gamma((nu + 1) / 2, (nu + r^2 / sigma^2) / 2) given the residual, and the gaussian part reads only their sum
Lambda_i = sum_j lambda_ij, which is Gamma(w (nu + 1) / 2, same rate). Taking that conditional for real w is
exact too: the joint of r and Lambda it implies, with a normalizing factor K(nu, w) free of r and sigma,
integrates over Lambda to t_nu(r / sigma)^w / sigma^w, so the trees and sigma are exact for any positive w given
nu, sigma's posterior reading sum(w) observations and the sum of Lambda r^2. With m_i the row's count
(w_i under frequency, 1 under precision, rows out of the likelihood skipped), the nu step on its grid reads
sum(m_i) where it reads the row count, sum(m_i log Lambda_i) and sum(Lambda_i) where it reads the lambda
statistics, plus sum_i [m_i lgamma((nu + 1) / 2) - lgamma(m_i (nu + 1) / 2)], which depends on the weights alone
and is computed per grid point when they are installed. At m_i = 1 that term is zero and the draw is the shipped
one ([`ResidualDfPrior::drawIndex`](../../src/bartcore/model.hpp)).

The families in the integer rows are coherent frequency models but need a latent per copy, so their cost grows
with the total weight, as logistic's already does (dec-B202 measured 160 ns per trial and capped a row at one
million). Real weights there need a real-shape Polya-Gamma or truncated-normal-sum sampler, the
approximate-distributions arc. Integer-weight probit is TODO weighted-binary, "decision-gated, scheduled after
1.0-0 (VD 2026-08-28)".

## Open calls

### 1. The family list for 1.0-0

Background: dec-B125 asks for frequency weights "for models where it makes sense". They make sense on every
family above; they differ from what the package does today only on the families with a residual scale, and
the exact draw is cheap only where no latent is needed per copy.

- (a) The choice on the gaussian, Student-t and heteroscedastic families, real-valued, with every composition
  on them (several forests, BCF, linear and gp leaves); logistic counts and the 0/1 weights of the latent families
  are frequency weights already and keep what they take; the integer cases (probit, ordinal, nbinom, multinomial,
  aft) join after 1.0-0 with weighted-binary. Cost: slices 1 to 4 below, about 1,700 lines with their tests.
- (b) (a) plus whole-number frequency weights on probit, ordinal, nbinom, multinomial and aft before 1.0-0.
  Cost: a per-copy latent draw for each family, its exact-posterior gate, and the replication equivalence;
  about 200 to 300 lines per family, aft also passing a within-row sum of squares to sigma; it pulls forward a
  scheduled item.
- (c) Gaussian alone, as the TODO's certain case. Cost: Student-t and the variance forest, which share the
  gaussian's residual prior home (dec-B123), would refuse a choice their model takes as easily.

Recommendation: (a).

### 2. The argument's name and shape

Background: the reading is part of the model, the likelihood, not of the data; a sampler holding weights reads
them under the model it was created with. Base R has no frequency spelling; glm lets the family decide what
prior weights do. statsmodels and Stata name the two kinds, and dbarts already keeps the residual scale's prior
on the family object (dec-B123). The shipped header has no weight entry today; the plain-C creation entry that
dec-B85 places before 1.0-0 will carry the weights in its specification.

- (a) A family setting, weights = "precision" (the default) or "frequency", on gaussian(), student() and, when
  they take weights, aft() and hurdle.lognormal(); the variance forest follows the family. bart (formerly
  bart2), dbarts, dbartsSpec and xbart gain nothing; dbartsData keeps one weight vector; the sampler's
  setWeights, setData and predict read the vector under the fit's reading; the header's specification gains
  one enumeration field on its family part in dec-B85's change, so consumers rebuild once. stats::gaussian()
  reads precision, as lm does; bartBT, the 0.9-x door, keeps precision and has no setting. Cost: a user looking
  only at the weights item must be pointed to the family, which the weights help does in its first sentence.
- (b) A fitting-function argument beside weights, such as weights.type = "frequency". Cost: a new formal on
  bart, dbarts, dbartsSpec and xbart and in every consumer that forwards arguments, and a default that has to
  be "precision" on some families and the only reading on others.
- (c) A marker on the value, weights = frequency(w), in the manner of chisq() inside its argument. Cost: the
  reading travels on the data, so setWeights could change the model and has to refuse it; model.frame and
  subsetting drop the class.
- (d) A second vector, frequency.weights, as statsmodels has freq_weights and var_weights, exclusive with
  weights. Cost: two slots, two setters and a rule for both at once, which no user case asks for.

Recommendation: (a).

### 3. A frequency-weighted row's predictive draw

Background: a row's predictive draw today is a draw of what the row records: a precision-weighted gaussian row
an observation of variance sigma^2 / w, a logistic row the successes among its w trials (and predict's binary
weights are trial counts, dec-B470), a multinomial count row a count vector (dec-B177). A row of weight 0 draws
NA by every route (dec-B472). The log-likelihood of a count row is w times the one-observation value (dec-A116).
The TODO's sketch drew a frequency row at sd sigma, one observation, as brms and rstanarm do.

- (a) What the row records: for a frequency gaussian row, the mean of its w observations, sd sigma / sqrt(w),
  the same formula as precision at frequency's sigma. One rule over every family, logistic unchanged, the
  predictive code unchanged; a weight given to predict on a frequency fit is the number of new observations
  whose mean is drawn, as on a binary fit it is the number of trials. A single new observation at a training
  row is predict(fit, type = "ppd", weights = 1), the door dec-B470 and dec-B472 opened.
- (b) One observation: sd sigma on a gaussian row. To hold one rule, logistic's training draw moves from
  Binomial(w, p) to Bernoulli(p) (logistic is not in 0.9-34, so no shipped behavior moves), against the
  premise of dec-B470 for binary training rows; predict's weights then have no use on a frequency gaussian
  fit and are refused there by name.
- (c) One observation on the scale families, logistic as now: the TODO's sketch. Two conventions for the one
  reading.

Recommendation: (a). The log-likelihood is w times the one-observation value under every option (dec-A116).

### 4. Membership on probit and ordinal: weights, or its own name

Background (the maintainer's question of 2026-10-05): 0/1 weights on probit, ordinal and nbinom install the
active-row mask and are cleared from the data object, so `s$data` carries no weights and `s$setData(s$data)`
drops the rows out (TODO mask-on-own-data-object); xbart refuses 0/1 weights on those families, its folds having
no mask channel; predict reads probit weights as trial counts. Since 0 and 1 mean out and in once under every
reading, membership is the part of the frequency reading those families can draw now, and whole-number weights
on them, when they come, extend it.

What the one argument costs today: one paragraph per family in every weights item (bart, dbarts, the sampler,
stan4bart's), the policy branches, the data object that drops the mask, and the xbart refusal. With the
finding above the weights help becomes one rule, "a weight counts a row w times in the likelihood; on a family
with a residual scale it divides the variance instead unless the family says frequency", plus which values each
family can draw.

- (a) Keep 0/1 weights, and let them ride the data object as weights on every family, the engine installing
  them as the mask; `setActiveRows` stays the door of a sampler embedded in a larger sampler and composes with
  the weights. Cost: `s$data@weights` is no longer NULL on a mask family; resolves mask-on-own-data-object;
  about 150 lines with tests.
- (b) As now.
- (c) Membership under its own name on the fitting functions (such as active = ), weights refused on probit and
  ordinal until whole numbers can be drawn. Cost: a new formal on bart, dbarts, dbartsData and xbart; stan4bart,
  whose probit rule passes 0/1 weights, and any caller doing the same break; the name means what 0/1 weights
  mean anyway.

Recommendation: (a).

### 5. Zero weights

Background: a weight of 0 takes a row out of the likelihood and keeps it in the design (dec-B238, dec-B251),
out of sigma's count (dec-A11), NaN in the pointwise log-likelihood (dec-A116), NA in the predictive draws
(dec-B472). Under frequency a 0 is no copies, the same statement. Logistic alone refuses a 0, telling the user to
drop the row, though its engine takes the mask
([`LogisticResponse::workingWeights`](../../src/bartcore/model.hpp)).

- (a) As now on every family and under both readings, and logistic accepts 0 as the mask. Cost: about 80 lines
  with tests.
- (b) As now, logistic still refusing.

Recommendation: (a).

### Settled by rulings, not open

- Non-integer frequency weights are taken wherever the family's draw is exact for them: dec-B103, "weights may
  be fractional wherever the family supports real weights"; the integer cases come from the binary-count and
  negative-binomial rulings. Under option 1(a) that is the three recommended families. brms, rstanarm and
  survey take them; Stata refuses non-integer fweights and offers iweights for the same power likelihood. The
  help says that a frequency weight below 1 makes the row count for less than an observation in sigma's
  posterior, and that weights normalized to sum to 1 leave sigma about one observation's worth of data.
- pdbart averages unweighted unless given average.weights (dec-B214); the help says that average.weights = w
  gives the expanded data's average.
- Precision stays the default (dec-B125), so no existing fit and no scenario of dec-B120's comparison moves.

## Calls made

1. The engine carries one reading flag per sampler and derives a per-row count from it, 1 for a positive
   weight under precision and w under frequency, used by sigma's degrees of freedom
   ([`ChiSquaredScalePrior::drawSigmaSqFromPosterior`](../../src/bartcore/model.hpp), which takes the positive
   count today), a variance leaf's count ([`ConstantVarianceLeaf`](../../src/bartcore/model.hpp)) and
   Student-t's shape and nu step. Alternative: a second per-row vector of counts beside the precisions, which
   would allow both kinds at once as statsmodels does; no user case and n more doubles.
2. Student-t's nu step is the shipped grid draw generalized as derived above, exact for real w and identical to
   the shipped draw under precision. Alternative: the marginal-t grid draw, nu given sigma and the fit with the
   scales integrated out, also exact and better mixing, but a second mechanism and a changed draw at w = 1.
3. Under frequency every quantity the package reads off the rows is read off the expanded rows: the starting
   sigma's regression with residual df sum(w) - p ([`residualStandardError`](../../R/utility.R) reads lm.wfit's
   positive-row df today), its stand-in a weighted sd, a linear or gp leaf covariate's center and scale, a
   forest's default sd, the quantile cut rule's quantiles and, when nbinom takes the reading, its center.
   Where sum(w) - p is not positive the stand-in is used, as where the regression is undefined today.
   Alternative: read them as precision does, over the positive rows unweighted; the fit would then differ from
   the expanded-data fit through its prior's anchors.
4. The reading is model: it lives on the model record, is fixed at creation as the family is, and a setModel
   that changes it is refused. A fit object saved before it reads precision. Alternative: on the data object,
   where setWeights or setData could change the model.
5. A saved state's weights digest covers the reading, so a state installed into a sampler of the other reading
   redraws Student-t's scales as a weight change does
   ([`TResponse::reapplyWeights`](../../src/bartcore/model.hpp)); gaussian and heteroscedastic states carry no
   latent and install as they are. Alternative: refuse such an install, which no other weight change does.
6. A row of frequency weight 0 reports NaN in the pointwise log-likelihood, as a precision 0 and a masked row
   do ([`pointwiseLogLikelihood`](../../R/generics.R)). Alternative: 0, the product, which a tool such as loo
   would count as an observation.
7. xbart keeps a row's copies in one fold, weights the held-out loss by w as it does now
   ([`xbartLossFunction`](../../R/xbart.R)), and takes each fold's starting sigma under the fit's reading.
   Alternative: splitting a row's copies across folds, defined only for whole numbers and leaking identical
   observations into training.

## What moves, place by place

Under the recommendations, on a fit that chooses frequency; a precision fit moves nowhere.

| where a weight is read | precision (default, unchanged) | frequency |
|---|---|---|
| sigma's posterior | df nu + positive rows, SSR sum(w r^2) | df nu + sum(w), same SSR |
| sigest and the chi prior's anchor | lm.wfit, df positive rows - p | df sum(w) - p; the prior's scale moves with it |
| starting sigma's stand-in | unweighted sd over positive rows (TODO starting-sigma-fallback-zero-weights) | weighted sd |
| leaf updates, tree moves, k | precision w | unchanged |
| variance forest | leaf count: positive rows | leaf count: sum(w); its leaf prior follows sigest |
| Student-t | composite w lambda, nu from lambda statistics | composite Lambda, shape w (nu + 1) / 2, nu step generalized |
| predictive draw | sd sigma / sqrt(w); 0 gives NA | the same formula at frequency's sigma (open call 3 (a)) |
| pointwise log-likelihood | log N(y; f, sigma^2 / w); 0 gives NaN | w log N(y; f, sigma^2); 0 gives NaN |
| fitted, residuals | unweighted | unchanged |
| xbart | loss weighted by w; fold sigest precision | loss unchanged; fold sigest frequency |
| response range; the default cut rule's grid | range over positive rows; grid over every row | unchanged: neither depends on how many times a row appears |
| leaf covariate standardization, forest default sd, the quantile cut rule | unweighted | weighted by w |
| state, copy, reload | digest of weights | digest of weights and reading |
| logistic, probit, ordinal, nbinom | unchanged by the choice | unchanged (open calls 4 and 5 move the 0/1 data object and logistic's 0) |

Under the alternatives: 1(b) adds the integer families' latent draws, the replication cost per sweep and an
exact-posterior arm each; 1(c) leaves Student-t and the variance forest at precision. 2(b) to 2(d) move the
same rows through another door, adding a formal or a slot to every entry point and consumer that forwards one.
3(b) moves the predictive row to sd sigma, logistic's training draw to Bernoulli(p), and refuses predict's
weights on a frequency gaussian fit; 3(c) moves only the first. 4(b) and 5(b) move nothing; 4(c) moves every
probit and ordinal call that passes 0/1 weights, stan4bart's included.

Consumers, read without editing: stan4bart passes its weights to the BART part and reads them as precisions in
its own Stan target ("the weighted normalizer's N-not-sum(w)" in parametric_model.hpp), so it is unaffected at the
default and would need its own normalizer to offer frequency; its predict(newdata, type = "ppd") ignores weights
(stan4bart TODO predict-ppd-weights), which draws one observation: option 3(b)'s draw, and a gap under
precision and under option 3(a). bartCause passes population weights to bart and stan4bart's response fits and
uses them in its estimators; its sampleFromPPD draws at sigma whatever the weights, the same gap. Population
weights are neither reading's case: design-based inference takes its variance from the design
([weighted-logistic.md](weighted-logistic.md) on inverse-probability weights), so neither consumer should move to
frequency by default. treatSens's BART path passes no weights to dbarts; bairrtt reads none.

## Cost

Measured for TODO frequency weights on a loaded machine, 20000 rows, 50 trees: the weighted path costs about 14
percent over the unweighted one, an all-ones mask nothing, and a masked row as much per sweep as an active one.
Frequency adds nothing per sweep on the recommended families: sigma's count is fixed between weight installs, a
variance leaf sums w in the pass that counts rows now, and Student-t's extra nu term is a per-grid constant
computed at install. The integer families of option 1(b) cost one latent draw per copy per sweep.

Engine changes beyond the gaussian df: the variance leaf's count, Student-t's shape, composite and nu step,
and the state digest. None changes a draw at the default.

## Slices a plan would take

1. Engine: the reading flag and per-row count, sigma's df, the variance leaf, Student-t, the digest; C++ tests
   against the expanded data and the w = 1 identity. About 250 lines, 250 of tests.
2. R and bridge: the family setting, the model record, the refusals, sigest and the read-off quantities,
   the log-likelihood, help and NEWS; the flag through the bridge. About 400 lines, 300 of tests.
3. Gates: frequency at integer w against the expanded rows in distribution on all three families, an
   exact-posterior arm for sigma on a single tree, an SBC arm. About 250 lines.
4. Header: the enumeration field, in dec-B85's change. About 40 lines and the ABI hash.
5. Open calls 4 and 5, if taken: 0/1 weights on the data object, logistic's 0. About 230 lines with tests.

The integer families are a plan of their own after 1.0-0, with weighted-binary.
