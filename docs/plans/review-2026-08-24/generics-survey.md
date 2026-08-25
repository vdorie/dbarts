# j2 survey: plot(), extract(type = "loglik") and as_draws_* for bartMultinomial,
# bartOrdinal, bartNegbin, bartHurdle

Every claim is tool-verified: dbarts by reading branch `bartcore` @ 0045507c and by fitting all
four families against a clean `git archive` build in a scratch library; external claims by reading
the cited source / URL or rendering the installed package's help.  "MEASURED" marks something I
ran.  Nothing from memory.

## 0. dbarts today - the baseline the four siblings must match

### 0.1 plot.bart (R/plot.R:47)
`plot.bart(x, plquants = c(0.05, 0.95), cols = c("blue", "black"), ...)`; errors by name when
`yhat.train` is absent.  Gaussian gets two panels.  LEFT (gated on `"sigma" %in% names(x)`):
`plotSigmaTrace(x$first.sigma, x$sigma)` sets `par(mfrow = c(1, 2))` and draws the residual-scale
trace - burn-in (`first.sigma`) RED, sampling run BLACK, bridged at a midpoint; multi-chain draws
one line per chain with `lty = i`, single-chain a point scatter. RIGHT: `x$y` on x against the
posterior MEDIAN of `yhat.train` on y, a vertical `plquants` interval per observation in
`cols[1]`, `ylab = "posterior interval for E(Y|x)"`, `abline(0, 1, lty = 2, col = cols[2])`.
Binary gets ONE full-device panel (no `sigma`, so no left panel): the same construction on
`probabilityFromLatents(...)` but plotted `qm` against `qm` - median p on BOTH axes, `xlab =
"median of p"`; a caterpillar, not an observed-vs-fitted panel.  `plot.rbart` is the same with
`yhat.train + ranef` substituted. `plot.bartMultinomial(x, cols = NULL, ...)` (R/plot.R:191)
already ships but is NOT a sibling: ONE full-device panel, a per-category trace of
`multinomialMeanProbArray(x)` (training-mean predicted probability per draw, chains pooled), `xlab
= "iteration"`, `ylab = "mean predicted probability"`, level legend.  MEASURED: `plot()` on a
bartOrdinal / bartNegbin / bartHurdle falls to `plot.default` and dies with `'x' is a list, but
does not have components 'x' and 'y'`.
### 0.2 extract(type = "loglik") (R/generics.R:37-160; man/bart.Rd:201)
Computed entirely in R by `pointwiseLogLikelihood(object, ev)`, keyed on `object$family`. The
engine's own channel (chain.hpp:5394) is reachable only through the flat C API
(`dbarts_results.logLikelihood`, dbarts.h:201), never from `bart2`.  Implemented: gaussian
`dnorm(y, ev, sigma/sqrt(w), log = TRUE)` ; probit/logistic `dbinom(y, 1, p, log = TRUE)` times
`w` for logistic (trial counts); aft, normal density for events and log upper tail for censored
rows.  Anything else errors by name - all four target families land there (feature-matrix.md cells
read `M generics.R:129`; [f25] records hurdle). Contract (man/bart.Rd:201,
inst/tinytest/test-pointwise-loglik.R): extract-only (`predict`/`fitted` error "type must be in");
`sample = "test"` refused by name; `dim(loglik) == dim(extract(type = "ev"))`; combined default =
samples-by-observations, "directly consumable by WAIC/PSIS-LOO ... such as those in the loo
package"; uncombined = chains-first `n.chains x n.samples x n.obs` with a documented `aperm(x,
c(2, 1, 3))` to reach loo's iteration x chain x N.  Offsets never appear in the formula - they are
already inside the stored `ev`.
### 0.3 as_draws_array / as_draws_df (R/diagnostics.R:177-185, R/hooks.R)
Registered dynamically into posterior's namespace on its onLoad (posterior is Suggests-only), for
`bart` and `rbart` ONLY, signature `(x, vars = c("sigma", "k", "tau"), ...)`.  `bartDrawsArray`
gathers named fields into one (iteration, chain, variable) array; scalar fields keep the bare
name, others become `field[inner]`.  Already written and reusable: `scalarFields` lists
`"dispersion"`; `ordinalCutpointsArray()` emits `cutpoint[1]..cutpoint[K-1]`;
`multinomialMeanProbArray()` emits the (iteration, chain, category) array
`summary.bartMultinomial` names `meanProb[<level>]`.  MEASURED: `posterior::as_draws_array()` on
any of the four fits today fails INSIDE posterior with "All list elements must be lists
themselves" - an ugly fall-through, not a refusal.

### 0.4 What the four fit objects carry (MEASURED, 2 chains x 6 draws, n = 40)
| | bartMultinomial | bartOrdinal | bartNegbin | bartHurdle |
|---|---|---|---|---|
| y | factor(n) OR n x K count matrix | ordered factor(n) | counts(n) | numeric(n), all rows |
| yhat.train | 12 x 40 x 3 softmax probs, levels on K | 12 x 40 x 4 probs | 12 x 40 mean counts mu | composed |
| latent | refused by name (non-identified) | 12 x 40 eta | 12 x 40 psi | via components |
| scalar params | NONE | cutpoints 12 x (K-1) | dispersion, length 12 | via components |
| burn-in channel | none | none | none | positive$first.sigma |
| components | - | - | - | $occupancy (class bart, probit on 1{y>0}, all n), $positive (class bart, gaussian on log y over the n+ rows, x.test = the FULL n rows) |

Constraining facts, MEASURED or read:
- `k` and `sigma` are NOT posterior quantities for ordinal / nbinom / multinomial: the engine run
  channel returns `k = NULL` and `sigma` identically 1.0 (fixed unit scale) and the packagers
  store neither.  The only scalar posterior parameters in the whole set are ordinal's K-1
  cutpoints and nbinom's r.  Multinomial has none.
- ordinal `cutpoints[, 1]` is pinned at exactly 0 (scheme A, ordinal.md section 2);
  `summary.bartOrdinal` reports it with rhat = NA.  nbinom `dispersion` draws are INTEGERS
  (MEASURED 3,4,4,4): r rides an integer grid.
- WEIGHTS ARE REFUSED for all four (R/spec.R:45-88 for ordinal and nbinom, each with an all-ones
  courtesy that NULLs the field; R/bart.R:831 multinomial; R/bart.R:1152 hurdle). No family's
  loglik needs a weight term and no fit object carries a `weights` field.
- OFFSETS are already inside every stored channel (multinomial's n x K category offset enters
  before the softmax; ordinal's eta and nbinom's psi are `f(x) + o`).  No explicit offset term is
  needed either.
- The ENGINE already implements ordinal (model.hpp:3304, the Phi difference) and nbinom
  (model.hpp:4480) pointwise densities, and deliberately does NOT for multinomial
  (`logLikelihoodIsDefined()` false - model.hpp:3710, chain.hpp:5396, multinomial.md:304).
  dbarts.h:174 still says the channel covers "gaussian, binary, and aft", understating it.
- The three K/count families' `extract` methods have NO `combineChains` formal; they return the
  fit's stored layout.  Only `extract.bartHurdle` has one.

## 1. Bayesian tree packages
BART 2.9.10 (installed; sources read).  **Zero plot methods** - `grep("^plot",
ls("package:BART"))` is empty (MEASURED); NAMESPACE registers 9 `predict` methods and nothing
else.  **No per-observation log-likelihood anywhere** (grep `loglik|logLik|LOO|waic` over all
eight fitters' deparsed sources: 0 hits, MEASURED).  The one LOO-ish quantity is `gbart`'s
undocumented scalar `LPML`: gbart.R:206-230 builds exactly the matrix dbarts already returns
- `dnorm(Y, res$yhat.train, SD, TRUE)`, or `dbinom(Y, 1, res$prob.train, TRUE)` - collapses
it to the naive harmonic-mean CPO and returns only `sum(log.CPO)`; the `res$CPO` line is commented
out.  Draws are bare draws x observations matrices, no dimnames, no chain index, no coda/posterior
support. mbart is a continuation-ratio model (K-1 forests on nested subsets) mbart2 is K
one-vs-all probit forests whose latents are exponentiated and normalised.  Both flatten n x K into
`ndpost x (n*K)`, observation-major `j = (i-1)*K + h`, with NO column names.  Neither is ordinal;
no cutpoints exist in BART. BayesTree 0.3-1.5: `plot.bart` as in 0.1; three registered methods,
all `plot`; no likelihood quantity crosses into R at all. bartMachine 1.4.2: no `plot` S3 method;
three exported ggplot2 functions. `plot_y_vs_yhat` puts ACTUAL on x and FITTED on y with interval
segments, a dashed identity line and empirical coverage in the title.  Posterior draws come back
OBSERVATIONS x DRAWS - transposed relative to every other package.  No loglik, no LOO/WAIC.
SoftBart 1.0.3 is the only tree package computing the right quantity: legacy `softbart()` returns
`loglik_train`, draws x observations, `= dnorm(Y, Yhat, sigma/sqrt(w), log = TRUE)`
(soft_bart.cpp:1257-1266) - undocumented, absent from its four modern fitters, on the NORMALISED y
scale, never wired to loo.  Net: dbarts's loglik channel is already the best in the tree family;
extend it rather than adopt anything from these packages.

## 2. Bayesian regression packages (the formulas worth copying)
brms (master, R/log_lik.R, R/distributions.R, inst/chunks/*.stan; not installed locally).
- Hurdle: every hurdle family routes through one helper `.dhurdle` (distributions.R:2580). `y == 0
  -> log(hu)`; `y > 0 -> log(1 - hu) + log g(y) - lccdf0`, with `lccdf0 = log(1 - g(0))` when
  `type == "int"` and `lccdf0 = 0` when `type == "real"`.  So the COUNT hurdles are zero-truncated
  and the CONTINUOUS ones (gamma, lognormal) are NOT - continuous positive support already puts
  zero mass at 0.  `hu` is P(y = 0), the opposite orientation to dbarts's pi = P(y > 0).
  `hurdle_lognormal` calls `dlnorm(x, meanlog, sdlog)`: the density on the NATURAL y scale,
  Jacobian included.
- Ordinal `log_lik_cumulative` (log_lik.R:880): `eta <- disc * (thres - mu)`; category 1 is `log
  F(eta_1)`, category K is `log(1 - F(eta_{K-1}))`, interior is `log_diff_exp(log F(eta_y), log
  F(eta_{y-1}))`.  Cutpoint-minus-linear-predictor is dbarts's own orientation.
- Categorical: `log softmax(eta)[y]`, reference category spliced in as a zero column. Multinomial
  `dmultinomial` (distributions.R:2836) is `lgamma(size+1) + rowSums(x*log_prob - lgamma(x+1))` -
  the multinomial coefficient IS included, `size = sum(x)` from the row.
- Negbinomial: `dnbinom(y, mu = mu, size = shape)`; brms's `shape` IS R's `size`, and rstanarm's
  `reciprocal_dispersion` is the same thing (log_lik.R:446).
- Weights: `log_lik_weight` is `x * weight`, applied at the end of every family (log_lik.R:1104);
  rstanarm's `.weighted` is identical.  Offsets never appear in log_lik - they are a summand of
  the linear predictor (predictor.R:18).
- Dimensions: roxygen at log_lik.R:22 promises "an S x N matrix"; rstanarm says the same; `?loo`
  (rendered from installed loo 2.10.0) says `loo(matrix)`: "An S by N matrix ... with all chains
  merged", `loo(array)`: "An I by C by N array".  dbarts already matches both.
- `plot.brmsfit` is `bayesplot::mcmc_combo(combo = c("hist","trace"))`, five variables per page;
  `plot.stanreg` defaults to `mcmc_intervals`.  Both parameter-space, both ggplot2.
rstanarm `.ll_polr_i` (log_lik.R:463) is brms's three-branch ordinal form computed on the
probability scale then logged - exactly what dbarts would do from its stored probabilities. No
hurdle, categorical or multinomial in rstanarm.  posterior 1.7.0 (installed, MEASURED):
`draws_array` dims are named iteration / chain / variable in that order; `draws_matrix` is draw x
variable, chains merged; `draws_df` adds `.chain`/`.iteration`/`.draw`.  Multi-dimensional
quantities are flat names with bracket indices, comma-separated, no space, first index fastest
(`theta[1]`, `p[1,2]`).  MEASURED: NON-NUMERIC indices are accepted - `variables()`,
`summarise_draws()` and `as_draws_df()` all round-trip `p[1,lo]` and `meanProb[lo]`; only
`extract_variable_array()` reads the label as a dimension level, harmlessly.  Level names inside
brackets are safe.

## 3. Classical packages - the logLik convention
`?logLik` documents the return as "a number with at least one attribute df ... [and] nobs"; every
stats/MASS/nnet/nlme method listed there returns a scalar.
- pscl::hurdle (hurdle.R:9-63, 149): two SUMS added.  The count part IS zero-truncated -
  `sum(w*dpois(y,mu,log=TRUE)) - sum(w*log(1 - exp(-mu)))` over positives only; `ppois(0,
  lower.tail = FALSE, log.p = TRUE)` is the same quantity, used in the gradient and fitted-value
  code.  `$loglik` is a SCALAR; no per-observation vector exists.  `zeroinfl` is the mixture form
  (`log(pi + (1-pi) f(0))` for zeros), likewise scalar.
- countreg::zerotrunc (R-Forge only) writes the truncation as `dpois(y, mu, log=TRUE) - ppois(0,
  mu, lower.tail=FALSE, log.p=TRUE)`; scalar.
- MASS::glm.nb: scalar `twologlik/2`; `plot()` inherits `stats:::plot.lm`'s four panels (Residuals
  vs Fitted, Q-Q, Scale-Location, Residuals vs Leverage).
- ordinal::clm: scalar `logLik`, but `fm$fitted.values` IS the per-observation P(Y_i = y_i^obs)
  (verified `sum(log(fitted)) == logLik`).  No `plot.clm`; the only plots are parameter-space
  likelihood slices/profiles.  nnet::multinom: `logLik = -deviance/2`, scalar, `nobs =
  sum(weights)`; no plot method.
- VGAM is the one real per-observation precedent: `logLik(fit, summation = FALSE)` returns an
  n-vector or n-row matrix with prior weights folded in, verified across poissonff, negbinomial,
  zipoisson, zapoisson, cumulative and multinomial.  The price: `logLik(fit)` itself returns a
  bare unattributed numeric.
None of pscl::hurdle, pscl::zeroinfl, ordinal::clm or nnet::multinom has a plot method; all four
die with the same `'x' is a list ...` error dbarts's three families die with today. There is no
classical convention to match.

## 4. RECOMMENDATIONS

Shared design rule, from plot.bart: LEFT panel = trace of the family's scalar posterior
parameter(s); RIGHT panel = an observed-vs-fitted posterior-interval panel with a reference line.
Keep the `plquants` / `cols` argument names and the by-name errors. Shared loglik rule:
extract-only, `sample = "test"` refused by name, no weight term (refused at fit time), no offset
term (already inside the stored channels), combined default = S x N, uncombined = chains-first
n.chains x n.samples x n.obs. Shared as_draws rule: expose the fit's SCALAR posterior parameters -
what `summary()` already tables - never the per-observation channels; that is the rule
`as_draws_array.bart` already follows, and it keeps as_draws and summary in agreement.

### 4.1 bartNegbin - the cleanest sibling
(a) plot: `par(mfrow = c(1, 2))`.  P1 the dispersion trace (`x$dispersion` is sigma-shaped, so
    `plotSigmaTrace`'s body applies), `ylab = "dispersion (r)"`.  P2 observed count `x$y` on x
    against the posterior median of `yhat.train` (= mu) with `plquants` bars and `abline(0, 1)`,
    `ylab = "posterior interval for E(Y|x)"` - plot.bart's gaussian panel verbatim on counts.
    REJECTED: a rootogram (the count community's real idiom, countreg/topmodels) - it needs a
    posterior-predictive binning choice and is not a sibling of plot.bart; record it as a door.
(b) loglik: `l[s,i] = dnbinom(y_i, size = r_s, mu = mu_si, log = TRUE)` with `r_s =
    dispersion[s]`, `mu_si = yhat.train[s,i]`.  Identical to the engine's form (model.hpp:4480)
    and to `dnbinom(y, size = r, prob = 1 - plogis(psi))`, both VERIFIED numerically on a fitted
    object, as was `mu == r*exp(psi)`.  Prefer the `mu =` spelling: R-core vocabulary, matches
    brms `shape` and rstanarm `reciprocal_dispersion`, and reads the two stored channels directly.
    `r_s` pairs with the fit draws by the same `chainFastest` recycling sigma uses (`"dispersion"`
    is already in `scalarFields`).  The log-exposure offset is inside psi and hence inside mu.
    Dims: exactly `dim(ev)`. Nothing new must be stored.
(c) as_draws: `as_draws_array.bartNegbin(x, vars = c("dispersion", "sigma", "k", "tau"))` =
    `as_draws_array.bart`, matching `summary.bartNegbin`'s defaults exactly.  Zero new machinery;
    `presentDrawsVars` already drops the absent names.  Variable: `dispersion`.

### 4.2 bartOrdinal
(a) plot: `par(mfrow = c(1, 2))`.  P1 the cutpoint traces, one line per FREE cutpoint
    (gamma_2..gamma_{K-1}); gamma_1 is pinned at 0 and drawing it is noise.  At K = 2 there is no
    free cutpoint (the fit is probit bit-for-bit), so the panel is skipped and the plot degrades
    to one full-device panel exactly as binary plot.bart does.  P2 the per-observation latent
    posterior interval: `latent.train` median + `plquants` bars, observations ordered by median
    eta on the x-axis, coloured by observed level, with horizontal dashed lines at the
    posterior-median cutpoints.  This is the ordinal-native observed-vs-fitted panel - it shows
    whether the latent separates the observed categories at the fitted thresholds - and needs no
    quantity the fit does not already store. REJECTED: median P(y=k|x) vs interval in K colours
    (the binary panel generalised) - it hides the ORDER, which is the whole point of the family,
    and duplicates what a multinomial plot shows.  REJECTED: observed category index vs the
    interval of `sum_k k*p_k` with an identity line - treats the ordered levels as numeric, which
    the model does not.
(b) loglik: `l[s,i] = log P(y_i = k | eta, gamma)`, `k = as.integer(y_i)`, `= log(Phi(gamma_{s,k}
    - eta_si) - Phi(gamma_{s,k-1} - eta_si))` with `gamma_0 = -Inf`, `gamma_K = +Inf`.  That IS
    `log(yhat.train[s, i, k])` - VERIFIED: the stored probability equals a hand-built Phi
    difference to machine precision.  One line, needs no cutpoints, works on a fit whose sampler
    was not kept.
On "already implemented and exported": the cumulative-probit density is implemented
TWICE - in the engine at model.hpp:3304 (`OrdinalResponse::computeLogLikelihood`,
reachable through `dbarts_results.logLikelihood`), and in R at R/bart.R:1735
`ordinalCategoryProbabilities(eta, cutpoints)`, which is INTERNAL, not in NAMESPACE.
So the C channel is exported, the R function is not; both agree with the
log(stored probability) form.
Orientation matches the ecosystem: brms and rstanarm both use `cutpoint - eta`.
Dims: `dim(ev)` WITHOUT the trailing K margin - i.e. `dim(latent.train)`, the same
shape `type = "ppd"` already returns for this family.
(c) as_draws: `vars = c("cutpoints", "sigma", "k", "tau")`, matching `summary.bartOrdinal`.  Zero
    new machinery: `bartDrawsArray` already special-cases `"cutpoints"` to
    `ordinalCutpointsArray`, which already emits `cutpoint[1]..cutpoint[K-1]` in posterior's
    bracket convention.

### 4.3 bartMultinomial
(a) plot: keep the existing per-category trace as P1 - it is this family's only convergence
    instrument, and `summary.bartMultinomial` already chose the same channel - and ADD P2 under
    `par(mfrow = c(1, 2))`.  P2: for every observation, the posterior median and `plquants`
    interval of the predicted probability of the OBSERVED category, plotted against that median
    with `abline(0, 1)` - the binary plot.bart panel exactly, `xlab = "median of p(observed
    category)"`.  For the COUNTS ingestion with multi-trial rows the observed proportion `y_ik /
    n_i` is a real number, so plot.bart's GAUSSIAN panel form applies over the n*K cells (observed
    proportion on x, interval on y, identity line) - richer, and the branch is `is.matrix(y) &&
    any(rowSums(y) > 1)`.  See ambiguity 2. Note this changes the current one-panel device
    behaviour; nothing is frozen pre-release.
(b) loglik: ONE formula covers both ingestions.  With `p[s,i,k] = yhat.train[s,i,k]` and `n_i =
    sum_k y_ik` (`= 1` for the label ingestion), `l[s,i] = lgamma(n_i + 1) - sum_k lgamma(y_ik +
    1) + sum_k y_ik * log p[s,i,k]` `= dmultinom(y[i,], prob = p[s,i,], log = TRUE)` - VERIFIED to
    machine precision against `dmultinom` on a fitted counts model.  The label case is its `n_i =
    1` reduction, where the coefficient is 0 and the value is `log p[s,i,y_i]`.  INCLUDE the
    coefficient: brms does (`dmultinomial`, distributions.R:2836), it makes the number a genuine
    log density, and elpd DIFFERENCES are unaffected either way.  The likelihood unit is the
    OBSERVATION (a row of counts, n_i trials) - not a trial and not a (row, category) cell - so
    LOO is leave-one-row-out; say so in the Rd.  Dims: `dim(ev)` without the trailing K margin,
    the `type = "ppd"` shape.  Nothing new must be stored.
Add a comment that this does NOT contradict the engine's NaN flag: the engine's
response model cannot see the K-blend, while R holds the softmax probabilities.  The
flag should stay false so nobody "fixes" it.
(c) as_draws: `vars = "meanProb"` by default, giving the K variables `meanProb[<level>]` that
    `summary.bartMultinomial` already builds via `multinomialMeanProbArray`.  This family has no
    scalar parameter at all, so the pooled per-category mean predicted probability is its only
    convergence instrument - and summary already chose it, so the two agree.  Do NOT default to
    `yhat.train` (n*K variables).  Level names inside brackets are safe (MEASURED, 0.3.4).
Rejected: BART's `(i-1)*K + h` flat, unnamed column layout - it loses the level names
and the observation/category structure posterior's bracket indices carry for free.

### 4.4 bartHurdle
(a) plot: `par(mfrow = c(2, 2))`, since there are two component fits and a composed model.  P1
    `plotSigmaTrace(x$positive$first.sigma, x$positive$sigma)` - verbatim reuse, red burn-in
    included, because the positive fit IS an ordinary bart gaussian fit carrying both channels
    (the helper's own `par(mfrow = c(1,2))` call must become optional).  P2 occupancy: median
    pi(x) vs its interval, `ylab = "posterior interval for P(Y > 0 | x)"` - the binary panel.  P3
    positive part: observed `log(y)` over the y > 0 rows on x against the interval of f(x) at
    those rows, identity line - the gaussian panel on the subset, on the scale the component
    actually fit.  P4 combined: observed y over ALL n (zeros included) against the interval of
    `E[y|x] = pi*exp(f + sigma^2/2)`, identity line - the only panel that shows the composed
    model, and the family's whole reason to exist. REJECTED: plotting each component separately
    (two devices, no combined view), and a 1x3 dropping P4.
(b) loglik: with pi_si = P(y_i > 0 | x_i) (occupancy, all n), f_si the positive part's log-scale
    linear predictor at ALL n rows (its x.test channel), sigma_s the positive fit's residual sd:
    y_i == 0 : l[s,i] = log1p(-pi_si)
    y_i  > 0 : l[s,i] = log(pi_si) + dlnorm(y_i, meanlog = f_si, sdlog = sigma_s, log = TRUE)
                      = log(pi_si) + dnorm(log y_i, f_si, sigma_s, log = TRUE) - log(y_i)
    NO TRUNCATION.  hurdle.md section 1 fixes the v1 positive part as LOGNORMAL, whose
    support is already (0, inf), so `1 - G(0) = 1`.  This is precisely brms's
    `type = "real"` branch, and precisely the contrast against the count hurdles, where
    brms, pscl and countreg all divide by `1 - g(0)`.  When `hurdle.nbinom` lands (a
    recorded door) THAT positive part must be truncated and this one still must not - state
    the rule in the Rd so the two are never conflated.
    ORIENTATION: brms's `hu` is P(y = 0); dbarts stores pi = P(y > 0), so the zero branch is
    `log1p(-pi)`, not `log(hu)`.  Note it, or users porting a brms formula will invert it.
    Reuse `hurdleParts()` verbatim: it already returns pi, f and sigma as flat draw-aligned
    vectors plus the uncombined shape, and already glues the two independently seeded
    components by draw index (valid because the parts are conditionally independent,
    hurdle.md section 0).  VERIFIED: `extract(fit$occupancy, "loglik")` equals
    `dbinom(1{y>0}, 1, pi, log = TRUE)` exactly, and the composed `ev` equals
    `pi*exp(f + sigma^2/2)` exactly.  The composed channel is NOT the sum of the two
    components' own loglik channels: the positive fit's channel covers only its n+ rows,
    sits on the log scale, and carries no Jacobian.  Dims: `dim(ev)` = S x n over ALL n.
(c) as_draws: implement, do not refuse.  `as_draws_array.bartHurdle(x, vars = c("sigma", "k",
    "tau"))` returning the union of the two components' present fields with `occupancy.` /
    `positive.` prefixes - the same two labelled blocks `print.summary.bartHurdle` already prints.
    Use a dot, not a bracket: posterior parses brackets as indices.  Refusing and telling the user
    to call as_draws on `$occupancy` / `$positive` is defensible (they ARE two independent fits)
    but is a worse answer than a two-line prefix, and one object per fit is posterior's whole
    point.

### 4.5 Conventions adopted vs rejected
ADOPTED for interoperability: loo/WAIC input shape (S x N combined, iteration x chain x N
uncombined - unanimous across brms, rstanarm and `?loo`, and already dbarts's); posterior variable
naming (flat bracket indices, comma-separated, first index fastest); brms's single
discrete/continuous truncation switch, adopted as the RULE dbarts states rather than as code;
brms's inclusion of the multinomial coefficient. REJECTED: brms's ordinal `disc` and rstanarm's
scobit `alpha` (dbarts's ordinal fixes the latent scale at 1 by construction, so there is no
parameter to report); weights multiplying the pointwise loglik (all four families refuse weights
at fit time - do not add a term "for symmetry"); BART's LPML (computes the right matrix, throws it
away, returns a scalar harmonic-mean CPO - copying it would foreclose PSIS-LOO and every
loo/bayesplot workflow); a `logLik()` method (base R's convention is a scalar total; VGAM's
`summation = FALSE` buys per-obs by breaking that contract); classical diagnostic plots (none of
pscl/ordinal/nnet has one at all; glm.nb only inherits plot.lm's normal-theory panels); brms's
`mcmc_combo` and rstanarm's `mcmc_intervals` defaults (both parameter-space ggplot2; keeping
plot.bart's half-parameter, half-observation base-graphics shape is what makes these four siblings
rather than a second idiom).

## 5. Ambiguities for the maintainer to decide
1. HURDLE SCALE (the one genuine ambiguity).  Natural-y density with the `-log y` Jacobian
   (`dlnorm`, brms's choice, comparable to any other model OF y - a hurdle-gamma, a tobit, a
   zero-inflated count) vs the log-y density without it (comparable only to other fits of log y).
   RECOMMEND natural. Consequence: it is not the sum of the components' own loglik channels, and
   differs from a log-y computation by the constant `-sum_{y>0} log y_i`, which shifts LPML/WAIC.
   Document that constant.
2. MULTINOMIAL PLOT P2.  One unconditional panel (median P(observed category) vs interval) or a
   branch to observed proportion `y_ik/n_i` on multi-trial count rows. RECOMMEND the branch; the
   unconditional version is simpler and defensible.
3. ORDINAL `cutpoint[1]`, pinned at 0, in as_draws and summary.  Keep and document the pin, or
   drop the degenerate column (it reports rhat = NA today).  RECOMMEND keep: the vector IS the
   model's gamma and a user reconstructing probabilities needs gamma_1.
4. ORDINAL loglik precision.  `log(stored Phi difference)` vs recomputing `log_diff_exp(log Phi,
   log Phi)` from eta and cutpoints (brms's numerically safer construction; rstanarm does it
   dbarts's way).  The difference is spent at packaging time, not loglik time, since dbarts stores
   the difference.  RECOMMEND the stored form for v1; record the tail-precision door.
5. NEGBIN plot P1 burn-in.  `bart2Negbin` drives one `run(n.burn, n.samples)`, so there is no
   `first.dispersion` and manufacturing one would restructure the run and move draws. RECOMMEND no
   red burn-in segment, noted in the Rd. Separately: r is integer on the current grid (MEASURED),
   so draw the trace with `type = "s"`.
6. `combineChains` on the three K/count families.  `extract.bartMultinomial`, `.bartOrdinal` and
   `.bartNegbin` have no such formal and return the fit's stored layout, so a fit packaged with
   `combineChains = FALSE` could never produce the S x N matrix loo wants.  RECOMMEND adding the
   formal (default TRUE) to all three and routing every channel through `reshapeChainedChannel` /
   `combineOrUncombineChains`, both of which already exist and already handle the K margin.
7. The K-margin dim contract.  man/bart.Rd currently promises `dim(loglik) == dim(extract(type =
   "ev"))`.  That is false by construction for multinomial and ordinal, whose loglik drops the
   trailing K margin.  Restate it as "the `ev` shape without the trailing category margin" - the
   `type = "ppd"` shape - and adjust test-pointwise-loglik.R's dim assertions.
