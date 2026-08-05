\name{bart}
\alias{bart}
\alias{bart2}
\alias{plot.bart}
\alias{extract}
\alias{predict.bart}
\alias{extract.bart}
\alias{fitted.bart}
\alias{residuals.bart}
\alias{extract.bartMultinomial}
\alias{fitted.bartMultinomial}
\alias{predict.bartMultinomial}
\alias{print.bartMultinomial}
\alias{residuals.bartMultinomial}
\alias{plot.bartMultinomial}
\alias{summary.bartMultinomial}
\alias{extract.bartOrdinal}
\alias{fitted.bartOrdinal}
\alias{predict.bartOrdinal}
\alias{print.bartOrdinal}
\alias{residuals.bartOrdinal}
\alias{extract.bartNegbin}
\alias{fitted.bartNegbin}
\alias{predict.bartNegbin}
\alias{print.bartNegbin}
\alias{residuals.bartNegbin}
\alias{extract.bartHurdle}
\alias{fitted.bartHurdle}
\alias{predict.bartHurdle}
\alias{print.bartHurdle}
\alias{residuals.bartHurdle}
\title{Bayesian Additive Regression Trees}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner.

\itemize{
  \item For numeric response \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary response \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi} denotes the standard normal cdf (probit link).
  \item \code{bart2} fits further response families through its \code{family} argument, described below: logistic, accelerated failure time (\code{family = "aft"}), discrete-time survival hazard (\code{family = "hazard"}, \code{"hazard.logistic"}), K-category multinomial (\code{family = "multinomial"}), ordered categorical (\code{family = "ordinal"}, a cumulative probit), negative-binomial counts (\code{family = "nbinom"}), and semicontinuous two-part responses (\code{family = "hurdle.lognormal"}).
}
}
\usage{
bart(
    x.train, y.train, x.test = matrix(0.0, 0, 0),
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2.0,
    power = 2.0, base = 0.95, splitprobs = 1 / numvars,
    binaryOffset = 0.0, weights = NULL,
    ntree = 200,
    ndpost = 1000, nskip = 100,
    printevery = 100, keepevery = 1, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100, printcutoffs = 0,
    verbose = TRUE, nchain = 1, nthread = 1, combinechains = TRUE,
    keeptrees = FALSE, keepcall = TRUE, sampleronly = FALSE,
    seed = NA_integer_,
    proposalprobs = NULL,
    keepsampler = keeptrees,
    resid.dist = gaussian)

bart2(
    formula, data, test, subset, weights, offset, offset.test = offset,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL,
    power = 2.0, base = 0.95, split.probs = 1 / num.vars,
    dart = FALSE,
    n.trees = 75L,
    n.samples = 500L, n.burn = 500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = TRUE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 1L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE, keepTrees = FALSE, 
    keepCall = TRUE, samplerOnly = FALSE,
    seed = NA_integer_,
    proposal.probs = NULL,
    monotone = NULL,
    interactions = NULL,
    blocks = NULL,
    variance = NULL,
    n.trees.variance = 40L,
    power.variance = NULL,
    base.variance = NULL,
    keepSampler = keepTrees,
    warm.start = NULL,
    n.grow.sweeps = 0L,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic", "aft", "multinomial",
               "ordinal", "nbinom", "hazard", "hazard.probit",
               "hazard.logistic", "hurdle.lognormal", "twopart"),
    missing = c("incorporate", "error"),
    resid.dist = gaussian,
    dispersion = NA_real_,
    breaks = NULL, max.rows = 1e7,
    \dots)

\method{plot}{bart}(
    x,
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    \dots)

\method{predict}{bart}(
    object, newdata, offset, weights,
    type = c("ev", "ppd", "bart"),
    combineChains = TRUE,
    n.threads,
    ci.level = NULL,
    \dots)

extract(object, \dots)
\method{extract}{bart}(
    object,
    type = c("ev", "ppd", "bart", "loglik", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE, \dots)

\method{fitted}{bart}(
    object,
    type = c("ev", "ppd", "bart"),
    sample = c("train", "test"),
    ci.level = NULL,
    \dots)

\method{residuals}{bart}(object, type = "ev", \dots)

\method{extract}{bartMultinomial}(
    object,
    type = c("ev", "ppd", "bart"),
    sample = c("train", "test"), \dots)

\method{fitted}{bartMultinomial}(object, type = c("ev", "class"), \dots)

\method{predict}{bartMultinomial}(
    object, newdata,
    type = c("ev", "ppd"),
    combineChains = TRUE, \dots)

\method{print}{bartMultinomial}(x, \dots)

\method{residuals}{bartMultinomial}(object, \dots)

\method{plot}{bartMultinomial}(x, cols = NULL, \dots)

\method{summary}{bartMultinomial}(object, \dots)
}
\arguments{
   \item{x.train}{
     Explanatory variables for training (in sample) data. May be a matrix or a data frame, with rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that \eqn{q} dummies are created if \eqn{q > 2} and one dummy is created if \eqn{q = 2}, where \eqn{q} is the number of levels of the factor. A sparse \code{Matrix::dgCMatrix} of ordinal predictors is also accepted, as is a data frame mixing ordinary columns with \code{Matrix::sparseVector} or \code{dgCMatrix} columns; see \code{\link{dbarts}}.
   }
   \item{y.train}{
     Dependent variable for training (in sample) data. If \code{y.train} is numeric a continuous response model is fit (normal errors). If \code{y.train} is a two-level factor (or logical, or two-level character) or has only values 0 and 1, then a binary response model with a probit link is fit. A factor with three or more levels is an error directing to \code{\link{bart2}}'s \code{family = "multinomial"} (unordered) or \code{family = "ordinal"} (ordered), neither of which \code{bart} fits.
   }
   \item{x.test}{
     Explanatory variables for test (out of sample) data. Should have same column structure as \code{x.train}. \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test}.
   }
   \item{sigest}{
     For continuous response models, an estimate of the residual standard deviation (residual standard error), \eqn{\sigma}, used to calibrate an inverse-chi-squared prior on the error variance. If not supplied, the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary. Same concept as \code{sigma} in \code{\link{dbarts}} and \code{\link{xbart}}.
   }
   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary. Default 3, Chipman, George, and McCulloch's calibration (see References); an aggressive choice rather than a derived one.
   }
   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate (\code{sigest}) is placed at. The closer the quantile is to 1, the more aggressive the fit will be as you are putting more prior weight on error standard deviations (\eqn{\sigma}) less than the rough estimate. Not applicable when \eqn{y} is binary. Default 0.90, from the same source as \code{sigdf}.
   }
   \item{resid.dist}{
     The residual error \emph{law} for a continuous (gaussian) response, passed through to \code{\link{dbarts}}. The default \code{gaussian()} gives normal errors; \code{student(df)} fits outlier-robust Student-t errors by a Gaussian scale-mixture augmentation, with \code{student(df = nu)} fixing the degrees of freedom and \code{student()} estimating them on a capped grid. Only continuous gaussian responses carry it (an error otherwise). Two caveats: with a flexible tree-forest mean, tail inference is partially confounded with fit flexibility, so an estimated \eqn{\nu} should be posterior-checked; and \eqn{\sigma} is then the \emph{conditional} scale, the marginal residual variance being \eqn{\sigma^2\,\nu/(\nu-2)}{sigma^2 nu / (nu - 2)}. See \code{\link{dbarts}} for the full description.
   }
   \item{dispersion}{
     The negative-binomial dispersion \eqn{r} (\code{family = "nbinom"} only; ignored otherwise), passed through to \code{\link{dbarts}}. \code{NA} (the default) estimates \eqn{r} on a capped positive-integer grid under a renormalized \eqn{\mathrm{gamma}(2, 0.1)}{gamma(2, 0.1)} prior; a supplied value fixes it and must be a positive integer (v1 ships the exact integer envelope, so a real fixed dispersion is refused). Larger \eqn{r} approaches the Poisson limit; smaller \eqn{r} is more overdispersed.
   }
   \item{breaks, max.rows}{
     Discrete-time hazard controls (\code{family = "hazard"} only; ignored otherwise), passed through to \code{\link{dbarts}}. \code{breaks} sets the period grid: \code{NULL} (the default) uses the sorted distinct observed times, a single positive integer bins at that many quantiles, and a numeric boundary vector gives explicit right-closed intervals. \code{max.rows} (default \eqn{10^7}) refuses an over-large person-period expansion, naming the coarsening levers.
   }
   \item{k}{
     For numeric \eqn{y}, \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}. The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}. For binary \eqn{y} fit with the default probit link, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from \eqn{\pm 3}{+/- 3}, the probit reference scale; \code{family = "logistic"} (\code{bart2} only) widens this to \eqn{\pm \pi \sqrt{3}}{+/- pi*sqrt(3)}, which is three standard deviations of the standard logistic latent variable (\eqn{\pi / \sqrt{3}}{pi/sqrt(3)} each), the same three-sd span probit's \eqn{\pm 3}{+/- 3} covers. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be. The value can be either a fixed number, or a \emph{hyperprior} of the form \code{chi(degreesOfFreedom = 1.5, scale = 2)}. For \code{bart2}, the default of \code{NULL} uses the value 2 for continuous responses and the \code{chi(1.5, 2)} hyperprior for binary ones, which centers the sampled \code{k} near the field-standard fixed value of 2 (prior median 1.9) while adapting to the data; pass \code{k = 2} for the fixed BART-package default, or \code{chi(1.5, Inf)} for the old improper prior. The default of 2 for continuous responses follows Chipman, George, and McCulloch's argument (see References) that with node prior standard deviation \eqn{\sigma_\mu = 0.5 / (k \sqrt{m})}{sigma_mu = 0.5 / (k * sqrt(m))} for \eqn{m} trees, \eqn{k} prior standard deviations of \eqn{f(x)} span the whole coded response range regardless of \eqn{m} -- so \eqn{k = 2} places that range at roughly a 95\% prior interval. The range scaling itself is an outlier-sensitive convention shared with BayesTree/bartMachine: extreme \eqn{y} values compress the effective prior on everything else, for which the published workaround is to log-transform or winsorize such values before fitting, or to let \code{k} adapt via the \code{chi} hyperprior.
   }
   \item{power}{
     Power parameter for tree prior. Default 2, Chipman, George, and McCulloch's empirical recommendation (see References) for the split-probability decay \eqn{base (1 + depth)^{-power}}{base * (1 + depth)^(-power)}, not derived from a formula.
   }
   \item{base}{
     Base parameter for tree prior. Default 0.95, from the same source as \code{power}.
   }
   \item{splitprobs, split.probs}{
     Prior and transition probabilities of variables used to generate splits. Can be missing/empty/\code{NULL} for equiprobability, a numeric vector of length equal to the number variables, or a named numeric vector with only a subset of the variables specified and a \code{.default} named value. Values given for factor variables are replicated for each resulting column in the generated model matrix. \code{numvars} and \code{num.vars} symbols will be rebound before execution to the number of columns in the model matrix.
   }
   \item{dart}{
     \code{TRUE} places a DART prior (a Dirichlet prior over the split-variable probabilities, inducing variable selection; Linero 2018) on the trees with default settings, using \code{power} and \code{base} for the growth probabilities; a prior created by \code{\link[=dbartsPriors]{dbartsPriors$dart}} supplies its own settings, overriding \code{power} and \code{base}. Cannot be combined with \code{split.probs}. Per-sample probabilities are returned in the \code{varprobs} element of the fit.
   }
   \item{binaryOffset}{
     Used for binary \eqn{y}. When present, the model is \eqn{P(Y = 1 \mid x) = \Phi(f(x) + \mathrm{binaryOffset})}{P(Y = 1 | x) = \Phi(f(x) + binaryOffset)}, allowing fits with probabilities shrunk towards values other than \eqn{0.5}. \eqn{\Phi} is the link for the default probit family; \code{bart2} fits with \code{family = "logistic"} use \code{plogis} instead, as noted in the value section below.
   }
   \item{weights}{
      An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. A probit fit (\code{bart}, or \code{bart2} with the default binary family) does not support weights, except that weights identically 1 are treated as absent; \code{bart2} with \code{family = "logistic"} treats them as observation counts and requires positive integers. For a weighted logistic fit, the \code{"ppd"} draw at an observation with weight \eqn{w} is the number of successes among \eqn{w} trials, \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} with \eqn{p} the fitted probability.
   }
   \item{ntree, n.trees}{
     The number of trees in the sum-of-trees formulation.
   }
   \item{ndpost, n.samples}{
     The number of posterior draws after burn in, \code{ndpost / keepevery} will actually be returned.
   }
   \item{nskip, n.burn}{
     Number of MCMC iterations to be treated as burn in.
   }
   \item{printevery, printEvery}{
     As the MCMC runs, a message is printed every \code{printevery} draws.
   }
   \item{keepevery, n.thin}{
     Every \code{keepevery} draw is kept to be returned to the user. Useful for \dQuote{thinning} samples.
   }
   \item{keeptrainfits, keepTrainingFits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of \code{x.train} are returned.
   }
   \item{usequants, useQuantiles}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived from the \code{x.train} variables. When \code{FALSE}, splits are determined using values equally spaced across the range of a variable. See details for more information.
   }
   \item{numcut, n.cuts}{
     The maximum number of possible values used in decision rules (see \code{usequants}, details). If a single number, it is recycled for all variables; otherwise must be a vector of length equal to \code{ncol(x.train)}. Fewer rules may be used if a covariate lacks enough unique values.
   }
   \item{printcutoffs, printCutoffs}{
     The number of cutoff rules to printed to screen before the MCMC is run. Given a single integer, the same value will be used for all variables. If 0, nothing is printed.
   }
   \item{verbose}{
     Logical; if \code{FALSE} suppress printing.
   }
   \item{nchain, n.chains}{
     Integer specifying how many independent tree sets and fits should be calculated. Default 1 for \code{bart}, BayesTree's historical single-chain convention; \code{bart2}'s (and \code{\link{rbart_vi}}'s) default is 4.
   }
   \item{nthread, n.threads}{
     Integer specifying how many threads to use. Depending on the CPU architecture, using more than the number of chains can degrade performance for small/medium data sets. As such some calculations may be executed single threaded regardless.
   }
   \item{combinechains, combineChains}{
     Logical; if \code{TRUE}, samples will be returned in arrays of dimensions equal to \code{nchain} \eqn{\times} \code{ndpost} \eqn{\times} number of observations. Default \code{TRUE} across \code{bart}, \code{bart2}, and \code{\link{rbart_vi}}; the number of chains is tracked regardless on the returned object's \code{n.chains} component (see \sQuote{Value}).
   }
   \item{keeptrees, keepTrees}{
     Logical; must be \code{TRUE} in order to use \code{predict} with the result of a \code{bart} fit. Note that for models with a large number of observations or a large number of trees, keeping the trees can be very memory intensive.
   }
   \item{keepcall, keepCall}{
     Logical; if \code{FALSE}, returned object will have \code{call} set to \code{call("NULL")}, otherwise the call used to instantiate BART.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. A \code{\link{set.seed}} beforehand suffices for reproducibility; supplying \code{seed} instead gives reproducible results without touching R's stream. See Reproducibility section below.
   }
   \item{proposalprobs, proposal.probs}{
     Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step. Defaults are 0.5, 0.1, 0.4, and 0.5 respectively.
   }
   \item{monotone}{
     Optional per-predictor monotonicity constraints, passed through to \code{\link{dbarts}}; see its \code{monotone} argument. A named vector selects predictors by model-matrix column name, each direction one of \code{"+"}/\code{"increasing"}/\code{1} or \code{"-"}/\code{"decreasing"}/\code{-1}; only numeric and ordered columns are eligible. A constraint forces birth/death-only proposals and a fixed \code{k = 2}. \code{NULL} (the default) fits the unconstrained model.
   }
   \item{interactions}{
     Optional interaction constraints, passed through to \code{\link{dbarts}}; see its \code{interactions} argument and \code{\link{interactions}}. Built with \code{interactions(max.order = , groups = , forbid = )}: cap the distinct predictors per root-to-leaf path, confine interactions to declared groups, and/or forbid named predictors from co-occurring. \code{NULL} (the default) fits the unconstrained model.
   }
   \item{blocks}{
     Optional block-additive constraint, passed through to \code{\link{dbarts}}; see its \code{blocks} argument and \code{\link{blocks}}. Built with \code{blocks(groups = , trees.per.group = )}: confine each whole tree to one declared group of predictors so the ensemble is exactly a sum of per-group functions. The groups must form a total, disjoint partition of the predictors, validated at fit time. \code{NULL} (the default) fits the ordinary model.
   }
   \item{variance, n.trees.variance, power.variance, base.variance}{
     Heteroscedastic BART (Pratola et al. 2020): passed through to \code{\link{dbarts}}. \code{variance} declares a second forest modeling the residual variance surface \eqn{s^2(x)} as a product of scaled-inverse-chi-squared leaves, coupled to the mean forest through the per-observation precision. It is a one-sided formula or column selector naming which predictors drive the variance (\code{~ x1 + x2}), \code{TRUE} for every predictor, or \code{NULL} (the default) for a homoscedastic fit. \code{n.trees.variance} is the number of variance trees; \code{power.variance}/\code{base.variance} their tree-structure prior (defaulting to the mean forest's). Gaussian responses only. The fit gains \code{s.train}/\code{s.test} (posterior draws of \eqn{s(x)}), and \code{predict} attaches an \code{"s"} attribute carrying \eqn{s(x)} for new data.
   }
   \item{keepsampler, keepSampler}{
     Logical that can be used to save the underlying \code{\link{dbartsSampler-class}} object even if \code{keepTrees} is false.
   }
   \item{warm.start}{
     \code{bart2} only. A previous fit whose forests seed the initial state instead of drawing trees from the prior: either a \code{bart} object fit with \code{keepSampler = TRUE} or a \code{\link{dbartsSampler-class}}. The donor must share this fit's predictors, tree count, and DART setting; a donor fit on a different cut grid is supported, its splits remapped onto this fit's grid (collapsing any the grid starves), the same way a data replacement remaps existing splits. Only its trees, \code{sigma}, and \code{k} carry over, and each chain starts from a different donor sample so multiple chains stay overdispersed. A warm start biases early draws toward the donor, so it shortens burn-in rather than removing it; keep a non-zero \code{n.burn}. See \code{\link{dbartsSampler-class}}'s \code{installTrees} method for finer control.
   }
   \item{n.grow.sweeps}{
     \code{bart2} only. A non-negative integer; when positive, the initial forest is built by \code{n.grow.sweeps} sweeps of XBART-style grow-from-root (He, Yalov and Hahn 2019) instead of a draw from the prior, an init-then-refine workflow that reaches a good fit in far fewer sweeps before the exact sampler takes over. The default \code{0L} draws trees from the prior, exactly as before. Like \code{warm.start}, a grow-from-root start biases the early draws toward the grown fit, so it shortens burn-in rather than removing it; keep a non-zero (if smaller) \code{n.burn}. The posterior the sampler targets is unchanged - only the starting point moves. Constant-leaf models only, and mutually exclusive with \code{warm.start} (both request an initialization). See \code{\link{dbartsSampler-class}}'s \code{growFromRoot} method.
   }
   \item{factors, family, missing}{
     As in \code{\link{dbarts}}; \code{bart2} only. \code{bart} keeps the historical behavior: factors always expand into indicator columns (\code{bart2}'s default, \code{"categorical"}, instead keeps each unordered factor as a single predictor - a different model-matrix representation, so a factor-predictor fit changes if moved from one interface to the other), binary responses are probit, and missing data are rejected.

     \code{family = "multinomial"} fits a K-category softmax classifier: K forests, one per category, coupled through an interleaved Polya-Gamma one-vs-rest augmentation (see \code{docs/design/multinomial.md}). \code{y.train} must be a factor (or character, coerced with \code{\link{factor}}) with at least two levels and no \code{NA}s or \code{NA} levels; \code{K = nlevels(y.train)} is inferred, never given explicitly, and \code{levels(y.train)} is captured at fit time and carried on every K-shaped output (\code{yhat.train}'s trailing dimension, \code{fitted}'s columns, the \code{fitted(type = "class")} factor). Alternatively, \code{y.train} may be an n x K numeric matrix of nonnegative integer trial counts, one row per observation and one column per category, with row sum (trials) \eqn{n_i \ge 1}{n_i >= 1}: category k's count is then a binomial(\eqn{n_i}{n_i}, p_k) draw under the same one-vs-rest augmentation, \code{K = ncol(y.train)} is inferred, and the levels carried on every K-shaped output are \code{colnames(y.train)} when present, otherwise \code{as.character(seq_len(K))}. A one-hot count matrix (every \eqn{n_i = 1}{n_i = 1}) is the same model as the factor response above, reproduced bit for bit. The formula interface is supported alongside the matrix one: \code{bart2(formula, data, family = "multinomial")} with a factor (or character) left-hand side routes to the factor response above, and a \code{cbind(c1, ..., cK) ~ x} left-hand side (the same idiom \code{glm}'s binomial family uses) routes to the count-matrix response above, \code{K} and the levels taken from the \code{cbind} column names; the right-hand side is coded exactly as it is for every other family's formula fit, including \code{predict} on a data frame \code{newdata}. Under the default \code{family = "auto"}, a factor (or character) response with three or more levels is detected and fit as multinomial, reporting the choice in a one-line message; a two-level factor (or logical) response resolves to probit, and a numeric response is unchanged. A count-matrix (\code{cbind(c1, ..., cK) ~ x}) response is only ever multinomial when \code{family = "multinomial"} is given explicitly - it is never inferred. \code{weights}, \code{offset}, \code{subset}, \code{samplerOnly}, \code{warm.start}, and \code{n.grow.sweeps} are all refused with an error naming the limitation. \code{test} is supported: an \code{x.test} of the same column structure as \code{x.train} reports the K-category softmax probabilities on the held-out rows as \code{yhat.test}, shaped and levels-named exactly like \code{yhat.train} (see \sQuote{Value}). \code{keepTrees} is supported too: it retains every one of the K forests' trees so \code{predict} can replay them at new predictors afterward, reproducing \code{yhat.test} bitwise when \code{newdata} matches the fit-time \code{test}; without \code{keepTrees}, \code{predict} errors. The per-forest leaf scale follows its own K-dependent calibration (the K = 2 anchor is the logistic scale \eqn{\pi\sqrt{3}}{pi*sqrt(3)} divided by \eqn{\sqrt{2}}{sqrt(2)}, for the identified pairwise log-odds); \code{k} is read from the usual node prior exactly as for any other family, but the node prior's \code{node.scale} itself is NOT consulted - the multinomial engine calibrates its own. The fit's class is \code{"bartMultinomial"}, not \code{"bart"}: see \sQuote{Value} below and the \code{extract}/\code{fitted}/\code{predict} methods for \code{bartMultinomial} objects.

    \code{family = "ordinal"} fits an ordered categorical response by a cumulative probit (a single forest, unlike multinomial's K): a latent \eqn{z = f(x) + \epsilon}{z = f(x) + eps} with \eqn{\epsilon \sim N(0, 1)}{eps ~ N(0, 1)} is cut at ordered thresholds \eqn{-\infty = \gamma_0 < \gamma_1 = 0 < \gamma_2 < \ldots < \gamma_{K-1} < \gamma_K = \infty}{-Inf = gamma_0 < gamma_1 = 0 < gamma_2 < ... < gamma_{K-1} < gamma_K = Inf}, so \eqn{P(Y = k \mid x) = \Phi(\gamma_k - f(x)) - \Phi(\gamma_{k-1} - f(x))}{P(Y = k | x) = Phi(gamma_k - f(x)) - Phi(gamma_{k-1} - f(x))}. The free cutpoints are sampled by a one-at-a-time marginal Metropolis update (latents integrated out) under a log-gap prior \eqn{N(0, 1.5^2)}{N(0, 1.5^2)}; at \eqn{K = 2} the model is exactly the probit family. \code{y.train} should be an ordered factor (\code{\link{is.ordered}}); its level order defines the category order, and \code{levels(y.train)} is captured at fit time and carried on every K-shaped output. Under the default \code{family = "auto"}, an ordered factor with three or more levels is detected and fit as ordinal, reporting the choice in a one-line message (an unordered factor stays multinomial - the two never overlap, the dispatch key is \code{is.ordered}); a two-level ordered factor is binary and resolves to probit. Given explicitly, \code{family = "ordinal"} also accepts an unordered factor or character response (with a message that the category order is taken from the level order - usually alphabetical, so order the levels yourself) and a numeric response, whose ordered levels are \code{sort(unique(y.train))}. \code{weights} are not supported (a weighted truncated-normal latent likelihood is not a coherent model); \code{samplerOnly}, \code{warm.start}, and \code{n.grow.sweeps} are refused with an error naming the limitation. \code{test} and \code{keepTrees} are supported as for multinomial, and \code{predict} replays the saved trees at new predictors, differencing the cumulative probit at the stored per-draw cutpoints. \code{rbart_vi} and \code{xbart} do not fit ordinal responses (grouped ordinal models and ordinal-scale cross-validation losses are recorded follow-ups). The fit's class is \code{"bartOrdinal"}, not \code{"bart"}: see \sQuote{Value} below.

    \code{family = "nbinom"} fits a non-negative integer (count) response by a negative-binomial model with the Polya-Gamma augmentation (a single forest, like logistic): the forest fits a log-odds latent \eqn{\psi = f(x) + o}{psi = f(x) + o}, the count law is \eqn{y_i \sim \mathrm{NB}(r, \mathrm{plogis}(\psi_i))}{y_i ~ NB(r, plogis(psi_i))} with dispersion \eqn{r}, and the mean is \eqn{E[y_i \mid x] = r e^{\psi_i} = r e^{f(x_i) + o_i}}{E[y_i | x] = r exp(psi_i) = r exp(f(x_i) + o_i)}. The offset \eqn{o_i} therefore enters the mean multiplicatively and is a log-exposure (\eqn{o_i = \log(\mathrm{exposure}_i)}{o_i = log(exposure_i)}). \code{y.train} must be a non-negative integer count, and \code{family = "nbinom"} is always explicit - a count carries no unambiguous class, so it is never inferred under \code{family = "auto"} (a numeric response there stays gaussian). The dispersion \code{r} is estimated by default on a capped positive-integer grid (\eqn{\{1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50\}}{{1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50}}) under a renormalized \eqn{\mathrm{gamma}(2, 0.1)}{gamma(2, 0.1)} prior by a closed-form discrete full conditional; passing \code{dispersion} as a positive integer fixes \eqn{r} instead (v1 ships the exact integer envelope, so a real fixed \code{dispersion} is refused - continuous dispersion is a recorded follow-up). Like probit, the latent scale is fixed at 1 and \code{weights} are not supported (exposure belongs in the offset, not in observation replication); \code{samplerOnly}, \code{warm.start}, and \code{n.grow.sweeps} are refused with an error naming the limitation. \code{test} and \code{keepTrees} are supported as for ordinal, and \code{predict} replays the saved trees at new predictors, reporting mean counts \eqn{r e^{\psi}}{r exp(psi)} at the stored per-draw dispersion (a log-exposure \code{offset.test} enters \eqn{\psi}{psi} additively). \code{rbart_vi} and \code{xbart} do not fit count responses (grouped negative-binomial models, real dispersion, and a Poisson family are recorded follow-ups). The fit's class is \code{"bartNegbin"}, not \code{"bart"}: see \sQuote{Value} below.

    \code{family = "hazard"} and \code{family = "hazard.logistic"} fit a discrete-time survival hazard model as \emph{ingestion sugar} over the binary families - the model adds no engine code (\code{"hazard.probit"} is an accepted alias for \code{"hazard"}). The response is the same \code{Surv} object or two-column \code{(time, status)} the \code{"aft"} family takes; each subject observed to time \eqn{t_i} is person-period-expanded into its at-risk rows (one per period \eqn{k = 1, \ldots, t_i}), each carrying the subject's covariates, an ordinal period column (appended last, so trees split on contiguous periods and borrow strength across time), and the binary indicator \eqn{y_{ik} = \mathrm{status}_i \cdot 1\{k = t_i\}}. The discrete hazard is \eqn{h(k \mid x) = g(f(x, k) + o)} for the chosen binary link \eqn{g} (probit for \code{"hazard"}, logistic for \code{"hazard.logistic"}), and the fit is EXACTLY an ordinary \code{"probit"}/\code{"logistic"} fit on the expanded rows - a hazard fit is byte-identical, draw for draw, to the binary fit on the hand-expanded design with the same seed. The time grid defaults to the sorted distinct observed times (the \pkg{BART} \code{surv.bart} convention); \code{breaks} coarsens it (a single integer bins at that many quantiles, a boundary vector gives explicit right-closed intervals), and \code{max.rows} guards the expansion size. Offsets are on the link scale and replicate per subject; \code{weights} replicate and then follow the chosen binary family's policy (probit refuses non-unit weights, logistic requires integer counts). A \code{Surv} response requires the family to be requested explicitly (\code{family = "auto"} selects \code{"aft"}); the matrix interface is required (no formula, \code{subset}, or \code{test} - expand held-out subjects with \code{survivalProbabilities(fit, times, newdata = )}, which needs \code{keepTrees}). The returned object is an ordinary \code{"bart"} fit whose \code{$family} records the binary link (so every prediction generic transforms correctly with no special case) and whose \code{$periods} element carries the grid; \code{predict}/\code{extract} return the per-(subject, period) hazards on the expanded rows, and \code{\link{survivalProbabilities}} produces survival-curve draws \eqn{S(t \mid x) = \prod_{k : \mathrm{periods}[k] \le t} (1 - h(k \mid x))}. \code{rbart_vi} and \code{xbart} do not fit hazard responses (grouped/frailty hazard and a \code{cloglog} link are recorded follow-ups).

    \code{family = "hurdle.lognormal"} (alias \code{"twopart"}, which resolves and prints as \code{"hurdle.lognormal"}) fits a semicontinuous two-part (hurdle) model for a non-negative response with exact zeros: an OCCUPANCY probit fit of \eqn{z = 1\{y > 0\}}{z = 1{y > 0}} over all n observations, glued at report time to a POSITIVE-PART gaussian fit of \eqn{\log y}{log y} over the subset \eqn{\{i : y_i > 0\}}{{i : y_i > 0}}. The two component fits share no parameters and are composed entirely R-side from two ordinary \code{bart2} fits at independently derived seeds - no engine code, and no coupling between the parts - so a shared variable-selection prior across the occupancy and positive parts is foreclosed by this composition (a recorded limitation, alongside a Duan-smearing retransformation and a heteroscedastic positive part, both follow-ups). \code{y.train} must be non-negative and finite and must carry at least one exact zero and one positive value (a response with no zeros, or none positive, is refused by name). By default \code{predict}/\code{fitted}/\code{extract} report the NATURAL (response) scale via posterior-predictive Monte Carlo, heteroscedasticity-aware: \eqn{E[y \mid x] = P(y > 0 \mid x)\,e^{f(x) + \sigma^2 / 2}}{E[y | x] = Pr(y > 0 | x) * exp(f(x) + sigma^2 / 2)} computed per posterior draw, consuming the positive part's per-observation \eqn{\sigma(x)} when a heteroscedastic \code{variance = ~x} surface is enabled on it (a single \eqn{\sigma} per draw otherwise). \code{type = "prob"} returns the occupancy probability \eqn{\pi(x)} through the probit link; \code{type = "link"}/\code{"log"} the positive part's log-scale linear predictor \eqn{f(x)}; \code{type = "ppd"} draws the proper BIMODAL predictive - a Bernoulli(\eqn{\pi}) spike at zero, else a lognormal draw - which the plain gaussian ppd path cannot produce. Fits currently use the matrix interface only (\code{bart2(x.train, y.train, family = "hurdle.lognormal")}); \code{weights}, \code{subset}, \code{offset}/\code{offset.test}, and \code{test} are all refused with an error naming the limitation - the positive-part fit is instead given the full training \code{x} as its own \code{x.test}, so its in-sample fitted values cover the zero rows it never trained on. \code{predict} requires \code{keepTrees = TRUE} and replays both saved forests at \code{newdata}. \code{dbarts()} does not fit this family - it composes two samplers, which only \code{bart2()} builds, so requesting it there is an error directing here. \code{rbart_vi} and \code{xbart} do not fit it either. The fit's class is \code{"bartHurdle"}, holding both component fits (\code{$occupancy}, a \code{"bart"} probit fit of the occupancy indicator; \code{$positive}, a \code{"bart"} gaussian fit of \eqn{\log y}{log y} on the positive subset) under their own \code{extract}/\code{fitted}/\code{predict}/\code{residuals}/\code{print} methods.
   }
   \item{formula}{
     The same as \code{x.train}, the name reflecting that a formula object can be used instead.
   }
   \item{data}{
     The same as \code{y.train}, the name reflecting that a data frame can be specified when a formula is given instead.
   }
   \item{test}{
     The same as \code{x.train}. Can be missing.
   }
   \item{subset}{
     A vector of logicals or indices used to subset of the data. Can be missing.
   }
   \item{offset}{
     The same as \code{binaryOffset}. Can be missing.
   }
   \item{offset.test}{
     A vector of offsets to be used with test data, in case it is different than the training offset. If \code{offset.test} is missing, defaults to \code{NULL}.
   }
   \item{object}{
     An object of class \code{bart}, returned from either the function \code{bart} or \code{bart2}.
   }
   \item{newdata}{
     Test data for prediction. Obeys all the same rules as \code{x.train} but cannot be missing.
   }
   \item{sampleronly, samplerOnly}{
     Builds the sampler from its arguments and returns it without running it. Useful to use the \code{bart2} interface in more complicated models.
   }
   \item{x}{
     Object of class \code{bart}, returned by function \code{bart}, which contains the information to be plotted.
   }
   \item{plquants}{
     In the plots, beliefs about \eqn{f(x)} are indicated by plotting the posterior median and a lower and upper quantile. \code{plquants} is a double vector of length two giving the lower and upper quantiles.
   }
   \item{cols}{
     Vector of two colors. First color is used to plot the median of \eqn{f(x)} and the second color is used to plot the lower and upper quantiles.
   }
   \item{type}{
     The quantity to be returned by generic functions. Options are \code{"ev"} - samples from the posterior of the individual level expected value, \code{"bart"} - the sum of trees component; same as \code{"ev"} for linear models but on the latent scale (probit or logistic, matching the fit's family) for binary ones, \code{"ppd"} - samples from the posterior predictive distribution, \code{"loglik"} - for \code{extract} only, the log-likelihood of each training observation at each posterior draw, and \code{"trees"} - a data frame with tree information for when model was fit with \code{keepTrees} equal to \code{TRUE}. For \code{"ppd"}, a weighted logistic fit draws the number of successes among the observation-count weight, \eqn{\mathrm{Binomial}(w_i, p_i)}{Binomial(w_i, p_i)} (see \code{weights}), and an aft (survival) fit draws on the log-time scale (the fit models \eqn{\log T}{log T}). For \code{"loglik"}, gaussian fits evaluate \eqn{y_i \mid x_i \sim N(\hat{f}(x_i), \sigma^2 / w_i)}{y_i | x_i ~ N(f(x_i), sigma^2 / w_i)} in logs at each draw of \eqn{f} and \eqn{\sigma}{sigma}, binary fits evaluate the Bernoulli log-likelihood of the fitted probability, multiplied for a weighted logistic fit by the observation-count weight \eqn{w_i}, and an aft fit contributes the log density for an event and the log survival tail \eqn{\log P(T > C)}{log P(T > C)} for a right-censored observation, both on the log-time scale. When chains are combined the result is a samples-by-observations matrix directly consumable by WAIC/PSIS-LOO implementations such as those in the \pkg{loo} package; the chains-first convention is kept, so a per-chain array (\code{combineChains = FALSE}, dimension chains-by-samples-by-observations) is reordered to the draws-by-chains-by-observations that \code{loo::relative_eff} expects with \code{aperm(x, c(2, 1, 3))}. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"ev"} and \code{"link"} can be used as a synonym for \code{"bart"}. For information on extracting trees, see the subsection below.
   }
   \item{sample}{
     Either \code{"train"} or \code{"test"}.
   }
   \item{ci.level}{
     For \code{fitted} and \code{predict}, an optional single number in \eqn{(0, 1)}. When \code{NULL} (the default) the posterior mean is returned as before. When given, the result is instead a matrix with one row per observation and columns \code{est} (the posterior mean), \code{ci.lower}, and \code{ci.upper} giving a symmetric credible band at that level, pooled over all samples and chains. The kind of interval follows \code{type}: \code{"ev"} yields a credible interval for the conditional mean (a probability for binary responses), \code{"ppd"} a prediction interval that additionally carries the residual noise (and so is wider), and \code{"bart"} a credible interval on the latent scale.
   }
   \item{\dots}{
     Additional arguments passed on to \code{plot}, \code{dbartsControl}, or \code{extract} when \code{type} is \code{"trees"}. Not used in \code{predict}.
   }
}
\details{
  BART is an Bayesian MCMC method. At each MCMC interation, we produce a draw from the joint posterior \eqn{(f, \sigma) \mid (x, y)}{(f, \sigma) | (x, y)} in the numeric \eqn{y} case and just \eqn{f} in the binary \eqn{y} case.

  Thus, unlike a lot of other modeling methods in R, \code{bart} does not produce a single model object from which fits and summaries may be extracted. The output consists of values \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{\sigma*} in the numeric case) where * denotes a particular draw. The \eqn{x} is either a row from the training data (\code{x.train}) or the test data (\code{x.test}).

  \subsection{Decision Rules}{
    Decision rules for any tree are of the form \eqn{x \le c}  vs. \eqn{x > c} for each \sQuote{\eqn{x}} corresponding to a column of \code{x.train}. \code{usequants} determines the means by which the set of possible \eqn{c} is determined. If \code{usequants} is \code{TRUE}, then the \eqn{c} are a subset of the values interpolated half-way between the unique, sorted values obtained from the corresponding column of \code{x.train}. If \code{usequants} is \code{FALSE}, the cutoffs are equally spaced across the range of values taken on by the corresponding column of \code{x.train}.

    The number of possible values of \eqn{c} is determined by \code{numcut}. If \code{usequants} is \code{FALSE}, \code{numcut} equally spaced cutoffs are used covering the range of values in the corresponding column of \code{x.train}. If \code{usequants} is \code{TRUE}, then for a variable the minimum of \code{numcut} and one less than the number of unique elements for that variable are used.
  }
  \subsection{End-node prior parameter \code{k}}{
    The amount of shrinkage of the node parameters is controlled by \code{k}. \code{k} can be given as either a fixed, positive number, or as any value that can be used to build a supported hyperprior. At present, only \eqn{\chi_\nu s} priors are supported, where \eqn{\nu} is a degrees of freedom and \eqn{s} is a scale. Both values must be positive, however the scale can be infinite which yields an improper prior, which is interpreted as just the polynomial part of the distribution. If \eqn{\nu} is 1 and \eqn{s} is \eqn{\infty}, the prior is \dQuote{flat}.

    For BART on binary outcomes, the degree of overfitting can be highly sensitive to \code{k} so it is encouraged to consider a number of values. The default hyperprior for binary \code{bart2}/\code{dbarts} fits, \code{chi(1.5, 2)}, centers the sampled \code{k} near the field-standard fixed value of 2 (prior median 1.9) while adapting to the data; pass \code{k = 2} for the fixed BART-package default. Crossvalidation may still be helpful, and running for a short time with a flat prior (\code{scale = Inf}) can show the range of \code{k} values the data are consistent with.
  }
  \subsection{Generics}{
    \code{bart} and \code{\link{rbart_vi}} support \code{\link{fitted}} to return the posterior mean of a predicted quantity, as well as \code{\link{predict}} to return a set of posterior samples for a different sample. In addition, the \code{extract} generic can be used to obtain the posterior samples for the training data or test data supplied during the initial fit.

    Using \code{predict} with a \code{bart} object requires that it be fitted with the option \code{keeptrees}/\code{keepTrees} as \code{TRUE}. Keeping the trees for a fit can require a sizeable amount of memory and is off by default.

    All generics return values on the scale of expected value of the response by default. This means that \code{predict}, \code{extract}, and \code{fitted} for binary outcomes return probabilities unless specifically the sum-of-trees component is requested (\code{type = "bart"}). This is in contrast to \code{yhat.train}/\code{yhat.test} that are returned with the fitted model.

    \code{residuals} returns \eqn{y} minus \code{fitted(object, type = type, sample = "train")}; the default \code{type = "ev"} gives ordinary response-scale residuals, so for a binary fit this is \eqn{y - P(Y = 1 \mid x)}{y - P(Y = 1 | x)}. Binary fits do not have \code{yhat.train.mean}/\code{yhat.test.mean} components (see \sQuote{Value}); use \code{fitted} to obtain the equivalent posterior-mean probabilities directly.
  }
  
  \subsection{Saving}{
    \code{\link{save}}ing and \code{\link{load}}ing a fitted BART object for use with \code{predict} requires that it be fit with \code{keeptrees}/\code{keepTrees} as \code{TRUE}, and that the sampler's tree state be captured before saving by calling \code{storeState()} on the sampler: \code{bartFit$fit$storeState()} (for \code{rbart_vi}, on each element of \code{fit}: \code{lapply(rbartFit$fit, function(f) f$storeState())}). The state is not captured automatically because it duplicates the trees on the R side and can exceed the fit's own sample blocks; it is materialized only on request. A fit saved without it reloads, but \code{predict} stops with an error naming \code{storeState()}.
  }
  
  \subsection{Reproducibility}{
    Every chain runs its own pseudo random number generator, so worker threads never touch R's generator and results do not depend on the thread count. Without a \code{seed}, chain generators are seeded from R's stream when the sampler is created, so calling \code{\link{set.seed}} beforehand makes results reproducible. Sampling itself never advances R's stream: after a fit, \code{\link{.Random.seed}} has moved only by the draws taken at creation.

    Setting \code{seed} makes results reproducible without involving R's stream at all: the seed drives a dedicated generator that hands each chain its own seed. A single-chain run with a given seed reproduces the first chain of a multi-chain run with the same seed.
  }

  \subsection{Extracting Trees}{
    When a model is fit with \code{keeptrees} (\code{bart}) or \code{keepTrees} (\code{bart2}) equal to \code{TRUE}, the generic \code{extract} can be used to retrieve a data frame containing the tree fit information. In this case, \code{extract} will accept the additional, optional arguments: \code{chainNums}, \code{sampleNums}, and \code{treeNums}. Each should be an integer vector detailing the desired trees to be returned. A further optional argument \code{newdata} routes a new set of predictors (in the same form accepted by \code{predict}) through the frozen trees so that the \code{n} column counts those observations instead of the training data.

    The result of \code{extract} will be a data frame with columns:
    \itemize{
      \item \code{chain}, \code{sample}, \code{tree} - index variables
      \item \code{n} - number of observations in node, drawn from the training data or from \code{newdata} when supplied
      \item \code{var} - either the index of the variable used for splitting or -1 if the node is a leaf
      \item \code{value} - either the value such that observations less than or equal to it are sent down the left path of the tree or the predicted value for a leaf node
    }
    The order of nodes in the result corresponds to a depth-first traversal, going down the left-side first. The names of variables used in splitting can be recovered by examining the column names of the \code{fit$data@x} element of a fitted \code{bart} or \code{bart2} model. See the package vignette \dQuote{Working with dbarts Saved Trees}.
  }
}
\value{
  \code{bart} and \code{bart2} return lists assigned the class \code{bart}. For applicable quantities, \code{ndpost / keepevery} samples are returned. In the numeric \eqn{y} case, the list has components:

  \item{\code{call}}{
    The matched call, or \code{call("NULL")} when the fit was made with \code{keepcall}/\code{keepCall} equal to \code{FALSE}.
  }
  \item{\code{yhat.train}}{
    A array/matrix of posterior samples. The \eqn{(i, j, k)} value is the \eqn{j}th draw of the posterior of \eqn{f} evaluated at the \eqn{k}th row of \code{x.train} (i.e. \eqn{f^*(x_k)}{f(x_k)}) corresponding to chain \eqn{i}. When \code{nchain} is one or \code{combinechains} is \code{TRUE}, the result is collapsed down to a matrix.
  }
  \item{\code{yhat.test}}{
    Same as \code{yhat.train} but now the \eqn{x}s are the rows of the test data.
  }
  \item{\code{yhat.train.mean}}{
    Vector of means of \code{yhat.train} across columns and chains, with length equal to the number of training observations.
  }
  \item{\code{yhat.test.mean}}{
    Vector of means of \code{yhat.test} across columns and chains.
  }
  \item{\code{sigma}}{
    Matrix of posterior samples of \code{sigma}, the residual/error standard deviation. Dimensions are equal to the number of chains times the number of samples unless \code{nchain} is one or \code{combinechains} is \code{TRUE}.
  }
  \item{\code{first.sigma}}{
    Burn-in draws of \code{sigma}.
  }
  \item{\code{varcount}}{
    A matrix with number of rows equal to the number of kept draws and each column corresponding to a training variable. Contains the total count of the number of times that variable is used in a tree decision rule (over all trees).
  }
  \item{\code{varprobs}}{
    Present only under a DART prior: the sampled split-variable probabilities, in the same layout as \code{varcount}. Each draw sums to one across variables.
  }
  \item{\code{sigest}}{
    The rough error standard deviation (\eqn{\sigma}) used in the prior.
  }
  \item{\code{y}}{
    The input dependent vector of values for the dependent variable. This is used in \code{plot.bart}.
  }
  \item{\code{fit}}{
    Optional sampler object which stores the values of the tree splits. Required for using \code{predict} and only stored if \code{keeptrees} or \code{keepsampler} is \code{TRUE}.
  }
  \item{\code{n.chains}}{
    Information that can be lost if \code{combinechains} is \code{TRUE} is tracked here.
  }
  \item{\code{k}}{
    Optional matrix of posterior samples of \code{k}. Only present when \code{k} is modeled, i.e. there is a hyperprior.
  }
  \item{\code{first.k}}{
    Burn-in draws of \code{k}, if modeled.
  }
  \item{\code{s.train}, \code{s.test}}{
    Present only under a heteroscedastic fit (\code{variance} supplied): posterior draws of the residual scale \eqn{s(x)}, laid out as \code{yhat.train}/\code{yhat.test}. See the \code{variance} argument.
  }
  \item{\code{binaryOffset}}{
    Present only for a binary fit: the offset value used.
  }
  \item{\code{periods}}{
    Present only for a discrete-time hazard fit (\code{family = "hazard"}/\code{"hazard.logistic"}): the ordered period grid the person-period expansion used. \code{\link{survivalProbabilities}} dispatches its survival-curve branch on this element.
  }
  \item{\code{family}}{
    The resolved response family (\code{"gaussian"}, \code{"probit"}, \code{"logistic"}, or \code{"aft"}). \code{predict}, \code{extract}, \code{fitted}, and \code{plot} use it to transform latent draws to probabilities. For an \code{"aft"} fit, predictions and fitted values are on the LOG-TIME scale: these functions return the linear predictor \eqn{E[\log T \mid x]}, exactly as for a gaussian fit of \eqn{\log T}, and never the time scale that \code{survreg} or \code{flexsurv} users might expect. The \code{type} aliases \code{"response"} and \code{"link"} are inert for \code{"aft"} - both return the same log-scale linear predictor; neither exponentiates. The posterior-median survival time is \code{exp} of the linear predictor, and \code{\link{survivalProbabilities}} gives survival-probability draws. A discrete-time hazard fit (\code{family = "hazard"}/\code{"hazard.logistic"}) records the binary link here (\code{"probit"}/\code{"logistic"}) - the fit IS that binary fit on the expanded rows - and carries a separate \code{$periods} element (the ordered period grid); \code{predict}/\code{extract} return the per-(subject, period) hazard through the link, and \code{\link{survivalProbabilities}} dispatches its survival-curve branch on the \code{$periods} marker.
  }

  In the binary \eqn{y} case, the returned list has the components \code{yhat.train}, \code{yhat.test}, and \code{varcount} as above, but not \code{yhat.train.mean}/\code{yhat.test.mean} - use \code{fitted} to get the posterior mean of \eqn{P(Y = 1 \mid x)}{P(Y = 1 | x)} instead. In addition the list has a \code{binaryOffset} component giving the value used.

  Note that in the binary \eqn{y}, case \code{yhat.train} and \code{yhat.test} are \eqn{f(x) + \mathrm{binaryOffset}}{f(x) + binaryOffset}. For draws of the probability \eqn{P(Y = 1 | x)}, apply the family's inverse link to these values: the normal cdf (\code{pnorm}) for \code{"probit"} fits, \code{plogis} for \code{"logistic"} ones.

  \code{bart2(family = "multinomial")} returns a \emph{different} list, of class \code{"bartMultinomial"} (never \code{"bart"} - its generics are separate, below, precisely so a K-category fit cannot silently flow through a method that assumes the single-category shape). Components: \code{call}, \code{family} (\code{"multinomial"}), \code{levels} (the resolved category names, length K: \code{levels(y.train)} for a factor response, or \code{colnames(y.train)} - falling back to \code{as.character(seq_len(K))} when unnamed - for a count-matrix response), \code{K}, \code{n.chains}, \code{n.trees}, \code{y} (the original response: the factor \code{y.train}, or the validated n x K count matrix when a count response was supplied), and \code{yhat.train} - the posterior draws of the K softmax probabilities, an array of dimension (\code{n.chains} \eqn{\times}, when \code{n.chains > 1} and \code{combineChains = FALSE}) \code{n.samples} \eqn{\times} number of training observations \eqn{\times} K, with the resolved category levels as the trailing dimension's names; \code{combineChains = TRUE} (the default) folds the chain margin into the samples margin as usual. \code{yhat.train} already holds PROBABILITIES, not a latent score - there is no \code{"bart"}-scale component to convert, unlike the binary families. When \code{test} was supplied, \code{yhat.test} is the same K softmax probabilities on the held-out rows, dimensioned like \code{yhat.train} with the training-observation margin replaced by the test one. \code{varcount} is the per-sample per-category split-usage channel: each category forest's per-draw variable counts, dimensioned like \code{yhat.train} but with the number of training predictors in place of the number of observations, and \code{colnames(x.train)} (when present) named on that margin. \code{bc} and \code{fit} are present only under \code{keepTrees} - the K-forest sampler handle and the host sampler \code{predict.bartMultinomial} replays through; absent \code{keepTrees}, neither survives the call and \code{predict} errors.

  Generics for a \code{"bartMultinomial"} fit: \code{fitted(object)} (\code{type = "ev"}, the default) returns the posterior-mean n \eqn{\times} K probability matrix, columns named by the fit's captured \code{levels}; \code{fitted(object, type = "class")} returns the argmax of that mean matrix as a factor over the original levels - a classification convenience. \code{extract(object, type = "ev", sample = "train")} returns the full \code{yhat.train} draws, and \code{sample = "test"} returns \code{yhat.test} (an error if the fit carries no test channel); \code{extract(object, type = "ppd")} draws one category per posterior draw from its probability vector, returned as an integer code (1-based, indexing the fit's captured \code{levels}) in an array shaped like \code{"ev"} minus the K margin. \code{extract}/\code{fitted} with \code{type = "bart"} error, naming the reason: the run records only the identified softmax probabilities (the raw per-category latent fits are non-identified, in the same sense as BCF's \eqn{a}, and are not recorded). \code{predict(object, newdata)} requires \code{object} fit with \code{keepTrees = TRUE} (error otherwise); it codes \code{newdata} to the training columns exactly as \code{predict.bart} does and returns a levels-named (\code{n.chains} \eqn{\times}) \code{n.samples} \eqn{\times} number of new observations \eqn{\times} K probability array, the \code{yhat.test}/\code{yhat.train} convention; \code{type = "ppd"} draws one category per posterior draw the same way \code{extract(object, type = "ppd")} does; \code{type = "bart"} is not offered, for the same non-identification reason. \code{residuals(object)} returns an n \eqn{\times} K matrix, the observed proportion (an indicator for a factor response, \code{y / rowSums(y)} for a count-matrix one) minus the fitted probability in \code{fitted(object)}, per category. \code{plot(object)} traces each category's training-mean predicted probability over the kept draws. \code{summary(object)} pools that same per-category mean-probability channel into posterior mean/sd/quantiles (R-hat/ESS when the \pkg{posterior} package is installed), the multinomial analog of \code{summary.bart}'s \eqn{\sigma}/k/\eqn{\tau} summary.

  \code{bart2(family = "ordinal")} likewise returns its own list, of class \code{"bartOrdinal"}. Components: \code{call}, \code{family} (\code{"ordinal"}), \code{levels} (the ordered category names, length K), \code{K}, \code{n.chains}, \code{n.trees}, \code{y} (the observed categories as an ordered factor over \code{levels}), \code{yhat.train} (and \code{yhat.test} when \code{test} was supplied) - the posterior draws of the K category probabilities, shaped and levels-named exactly as the multinomial arrays above - \code{latent.train} (and \code{latent.test}) - the corresponding draws of the latent \eqn{f(x)} on the unit-variance probit scale, shaped like a binary family's \code{yhat.train} - \code{cutpoints} - the posterior draws of the K - 1 finite thresholds \eqn{(\gamma_1 = 0, \gamma_2, \ldots)}{(gamma_1 = 0, gamma_2, ...)}, dimensioned draws \eqn{\times}{x} (K - 1), the ordinal analog of gaussian's \code{sigma}, from which probabilities at any latent value can be reconstructed - and \code{varcount}. \code{bc}, \code{fit}, and \code{cutpoints.raw} (the per-draw thresholds in the internal layout \code{predict} consumes) are present only under \code{keepTrees}.

  Generics for a \code{"bartOrdinal"} fit: \code{fitted(object)} returns the posterior-mean n \eqn{\times}{x} K probability matrix (columns named by \code{levels}); \code{fitted(object, type = "class")} the argmax category of that mean as an ordered factor; \code{fitted(object, type = "bart")} the posterior-mean latent. \code{extract(object, type = "ev")} returns the probability draws, \code{type = "bart"} (alias \code{"link"}) the latent draws, and \code{type = "ppd"} draws one category per posterior draw (1-based codes indexing \code{levels}). \code{predict(object, newdata)} requires \code{keepTrees = TRUE} and returns the probability draws at \code{newdata} (\code{type = "bart"} the replayed latent; \code{type = "ppd"} category draws); when \code{newdata} matches the fit-time \code{test} it reproduces \code{yhat.test}. \code{residuals(object)} returns the n \eqn{\times}{x} K matrix of observed-indicator minus fitted probability.

  \code{bart2(family = "nbinom")} likewise returns its own list, of class \code{"bartNegbin"}. Components: \code{call}, \code{family} (\code{"nbinom"}), \code{n.chains}, \code{n.trees}, \code{y} (the observed counts), \code{yhat.train} (and \code{yhat.test} when \code{test} was supplied) - the posterior draws of the mean counts \eqn{\mu = r e^{f(x) + o}}{mu = r exp(f(x) + o)}, shaped like a binary family's \code{yhat.train} - \code{latent.train} (and \code{latent.test}) - the corresponding draws of the log-odds latent \eqn{\psi = f(x) + o}{psi = f(x) + o} - \code{dispersion} - the per-draw dispersion \eqn{r}, the count analog of gaussian's \code{sigma} - and \code{varcount}. \code{bc}, \code{fit}, and \code{dispersion.raw} (the per-draw \eqn{r} in the internal layout \code{predict} consumes) are present only under \code{keepTrees}.

  Generics for a \code{"bartNegbin"} fit: \code{fitted(object)} returns the posterior-mean count per observation; \code{fitted(object, type = "bart")} the posterior-mean log-odds latent. \code{extract(object, type = "ev")} returns the mean-count draws, \code{type = "bart"} the latent draws, and \code{type = "ppd"} draws one count per posterior draw from \eqn{\mathrm{NB}(r, \mathrm{plogis}(\psi))}{NB(r, plogis(psi))}; \code{sample = "test"} selects the test channel. \code{predict(object, newdata)} requires \code{keepTrees = TRUE} and returns the mean-count draws at \code{newdata} (\code{type = "bart"} the replayed latent; \code{type = "ppd"} count draws), with an optional log-exposure \code{offset.test}; when \code{newdata} matches the fit-time \code{test} and no offset is given it reproduces \code{yhat.test}. \code{residuals(object)} returns the observed count minus the posterior-mean count per observation.

  For continuous response fits, the \code{plot} method sets \code{mfrow} to \code{c(1, 2)} and makes two plots. The first plot is the sequence of kept draws of \eqn{\sigma} including the burn-in draws. Initially these draws will decline as BART finds a good fit and then level off when the MCMC has burnt in. The second plot has \eqn{y} on the horizontal axis and posterior intervals for the corresponding \eqn{f(x)} on the vertical axis. For binary response fits, only this second kind of plot is drawn, with the posterior median of \eqn{P(Y = 1 \mid x)} on the horizontal axis and its posterior interval on the vertical axis.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2010)
   BART: Bayesian additive regression trees.
        \emph{The Annals of Applied Statistics}, \bold{4}(1), 266--298.
        \doi{10.1214/09-AOAS285}.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
   \url{https://www.rob-mcculloch.org}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.

He, J., Yalov, S., and Hahn, P.R. (2019)
   XBART: accelerated Bayesian additive regression trees.
        \emph{Proceedings of the 22nd International Conference on Artificial
        Intelligence and Statistics (AISTATS)}, PMLR \bold{89}, 1130--1138.

Linero, A.R. (2018)
   Bayesian regression trees for high-dimensional prediction and variable selection.
        \emph{Journal of the American Statistical Association}, \bold{113}(522), 626--636.

Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch, R.E. (2020)
   Heteroscedastic BART via multiplicative regression trees.
        \emph{Journal of Computational and Graphical Statistics}, \bold{29}(2), 405--417.
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com},
Robert McCulloch: \email{robert.mcculloch1@gmail.com},
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{pdbart}}
}
\examples{
## simulate data (example from Friedman MARS paper)
## y = f(x) + epsilon , epsilon ~ N(0, sigma)
## x consists of 10 variables, only first 5 matter

f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

## run BART
set.seed(99)
bartFit <- bart(x, y)

plot(bartFit)

## compare BART fit to linear matter and truth = Ey
lmFit <- lm(y ~ ., data.frame(x, y))

fitmat <- cbind(y, Ey, lmFit$fitted, bartFit$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))

## binary response with a logistic link, via bart2
set.seed(0)
n <- 60L
x.bin <- matrix(runif(n * 2), n, 2)
f.bin <- 3 * x.bin[,1] - 1.5
y.bin <- rbinom(n, 1L, plogis(f.bin))
fit.logit <- bart2(y.bin ~ x.bin, family = "logistic",
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L)

## fit with missing predictor values, using missing = "incorporate"
## (the default): every split rule learns a direction for NAs
set.seed(1)
x.na <- x.bin
x.na[1:5, 1] <- NA
y.na <- x.na[,2] + rnorm(n, 0, 0.3)
y.na[is.na(x.na[,1])] <- y.na[is.na(x.na[,1])] + 2
fit.na <- dbarts(y.na ~ x.na, control = dbartsControl(
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L))
samples.na <- fit.na$run()
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
