\name{bart2}
\alias{bart2}
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
\alias{summary.bartOrdinal}
\alias{plot.bartOrdinal}
\alias{extract.bartNegbin}
\alias{fitted.bartNegbin}
\alias{predict.bartNegbin}
\alias{print.bartNegbin}
\alias{residuals.bartNegbin}
\alias{summary.bartNegbin}
\alias{plot.bartNegbin}
\alias{extract.bartHurdle}
\alias{fitted.bartHurdle}
\alias{predict.bartHurdle}
\alias{print.bartHurdle}
\alias{residuals.bartHurdle}
\alias{summary.bartHurdle}
\alias{plot.bartHurdle}
\alias{print.summary.bartHurdle}
\title{Bayesian Additive Regression Trees with Extended Response Families}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner; see \code{\link{bart}} for the base gaussian/probit model definition. \code{bart2} is dbarts' own formula-first fitting interface: besides the default gaussian/probit response it reaches further families through \code{family} - logistic, accelerated failure time (\code{family = "aft"}), discrete-time survival hazard (\code{family = "hazard"}, \code{"hazard.logistic"}), K-category multinomial (\code{family = "multinomial"}), ordered categorical (\code{family = "ordinal"}, a cumulative probit), negative-binomial counts (\code{family = "nbinom"}), and semicontinuous two-part responses (\code{family = "hurdle.lognormal"}) - and, on its formula interface, an additional per-observation-modulated forest declared inline with a \code{forest()} term (see \sQuote{Formula Terms} below). These remaining own-class families are not reachable from \code{\link{bart}}: naming one of them to \code{bart}'s own \code{family} is refused, pointing back here.
}
\usage{
bart2(
    formula, data, test, subset, weights, offset, offset.test = offset,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL, prior.scale = NA_real_,
    power = 2.0, base = 0.95, split.probs = NULL,
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
    proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL,
    interactions = NULL,
    blocks = NULL,
    variance = NULL,
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
    tree.prior = NULL, node.prior = NULL, resid.prior = NULL,
    storage = c("double", "single"), updateState = TRUE)

\method{extract}{bartMultinomial}(
    object,
    type = c("ev", "ppd", "bart", "forest", "loglik"),
    sample = c("train", "test"),
    combineChains = TRUE, \dots)

\method{fitted}{bartMultinomial}(
    object, type = c("ev", "class", "bart"),
    ci.level = NULL, \dots)

\method{predict}{bartMultinomial}(
    object, newdata,
    type = c("ev", "ppd", "bart", "forest"),
    offset = NULL,
    combineChains = TRUE, ci.level = NULL, \dots)

\method{print}{bartMultinomial}(x, \dots)

\method{residuals}{bartMultinomial}(object, \dots)

\method{plot}{bartMultinomial}(x, cols = NULL, plquants = c(0.05, 0.95), \dots)

\method{plot}{bartOrdinal}(x, plquants = c(0.05, 0.95), cols = NULL, \dots)

\method{plot}{bartNegbin}(x, plquants = c(0.05, 0.95), cols = c("blue", "black"), \dots)

\method{plot}{bartHurdle}(x, plquants = c(0.05, 0.95), cols = c("blue", "black"), \dots)

\method{summary}{bartMultinomial}(object, \dots)

\method{summary}{bartOrdinal}(object, vars = c("cutpoints", "sigma", "k", "tau"), \dots)

\method{summary}{bartNegbin}(object, vars = c("dispersion", "sigma", "k", "tau"), \dots)

\method{summary}{bartHurdle}(object, vars = c("sigma", "k", "tau"), \dots)

\method{print}{summary.bartHurdle}(x, \dots)
}
\arguments{
   \item{formula}{
     An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}, optionally carrying one or more \code{forest()} terms that declare additional per-observation-modulated forests (see \sQuote{Formula Terms} below). For backward compatibility, can also be the \code{\link{bart}} matrix \code{x.train} (\code{data} then plays the role of \code{y.train}); see \code{\link{dbarts}}'s \code{formula} item for the full data-frame-to-predictor mapping, including sparse \code{Matrix::dgCMatrix} and \code{\link{sparseFactor}} columns.
   }
   \item{data}{
     An optional data frame, list, or environment containing predictors, used together with a formula \code{formula}. For backward compatibility, can also be the \code{\link{bart}} vector \code{y.train}.
   }
   \item{test}{
     An optional matrix or data frame with the same number of predictors as \code{data}, or \code{formula} in backward-compatibility mode; if column names are present, a matching algorithm is used. Refused, by name, when \code{formula} carries a \code{forest()} term (\sQuote{Formula Terms}): an amplitude-coupled fit has no per-observation test replay in this version.
   }
   \item{subset}{
     An optional vector of logicals or indices used to subset the data. Can be missing. A \code{forest()} term's basis is evaluated inside the same subsetted model frame as every other predictor (\sQuote{Formula Terms}); a basis given directly through \code{\link{dbarts}}'s \code{forests =} instead follows that argument's own subsetting rule (see \code{\link{forest}}'s \code{basis} item).
   }
   \item{weights}{
     An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. A probit fit (the default binary family here, and \code{\link{bart}}'s) does not support weights, except that weights identically 1 are treated as absent; a \code{family = "logistic"} fit treats them as observation counts and requires positive integers. For a weighted logistic fit, the \code{"ppd"} draw at an observation with weight \eqn{w} is the number of successes among \eqn{w} trials, \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} with \eqn{p} the fitted probability.
   }
   \item{offset}{
     An additive offset entering the family's mean or link, applied to training data; see \code{\link{dbarts}}'s \code{offset} item for the exact per-family formula. The same concept as \code{\link{bart}}'s \code{binaryOffset} (which only ever applies to that function's probit/logistic default family). Can be missing.

     For \code{predict} on a \code{bartMultinomial} fit, the same name carries a different shape: the per-category shift at the PREDICTED rows, an \code{nrow(newdata)} x K numeric matrix entering the K raw fits before the softmax, in the same layout as \code{bart2}'s own train-side \code{offset} above. It is never taken from the fit - these rows are not the fit's rows - so a fit trained with an \code{offset} REQUIRES one here, refusing by name rather than reporting the offset-free surface, which an explicit all-zero matrix asks for. \code{NULL} (the default) is the only value a fit trained without one takes by default, and passing the training offset back at the training rows reproduces \code{yhat.train}. A large negative entry is a shift like any other: the probabilities it drives are the exact softmax, so one far enough below its row's others underflows to zero (and its log to \code{-Inf}) as it would anywhere else.
   }
   \item{offset.test}{
     The equivalent of \code{offset} for test observations. Defaults to tracking \code{offset}: a scalar, or a vector already matching the test row count, is reused directly, and any other length is refused by name rather than silently recycled. See \code{\link{dbarts}}'s \sQuote{Test offset synchronization} details for how a live sampler keeps the two linked after fitting.
   }
   \item{sigest}{
     For continuous response models, an estimate of the residual standard deviation (residual standard error), \eqn{\sigma}, used to calibrate an inverse-chi-squared prior on the error variance. If not supplied, the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary. Same concept as \code{sigma} in \code{\link{dbarts}}; \code{\link{xbart}} spells it \code{sigest} too.
   }
   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary. Default 3, Chipman, George, and McCulloch's calibration (see References); an aggressive choice rather than a derived one.
   }
   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate (\code{sigest}) is placed at. The closer the quantile is to 1, the more aggressive the fit will be as you are putting more prior weight on error standard deviations (\eqn{\sigma}) less than the rough estimate. Not applicable when \eqn{y} is binary. Default 0.90, from the same source as \code{sigdf}.
   }
   \item{k}{
     For numeric \eqn{y}, \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}. The response is internally scaled to range from \eqn{-0.5} to \eqn{0.5}. For binary \eqn{y} fit with the default probit link, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from \eqn{\pm 3}{+/- 3}, the probit reference scale; \code{family = "logistic"} widens this to \eqn{\pm \pi \sqrt{3}}{+/- pi*sqrt(3)}, three standard deviations of the standard logistic latent variable. The value can be either a fixed number, or a \emph{hyperprior} of the form \code{chi(degreesOfFreedom = 1.5, scale = 2)}. The default, \code{NULL}, uses the value 2 for continuous responses and the \code{chi(1.5, 2)} hyperprior for binary ones, which centers the sampled \code{k} near the field-standard fixed value of 2 (prior median 1.9) while adapting to the data; pass \code{k = 2} for the fixed BART-package default, or \code{chi(1.5, Inf)} for the old improper prior. See \code{\link{bart}}'s \code{k} item, and the \sQuote{End-node prior parameter \code{k}} details there, for the full calibration argument and its outlier-sensitivity caveat.
   }
   \item{prior.scale}{
     Names the leaf calibration in response units instead of inheriting it from the range of the response: the prior standard deviation of the forest total \eqn{f} at \code{k = 1}, so the prior standard deviation in force is \code{prior.scale / k} and the prior mean is the response transform's shift. \code{NA} (the default) leaves the family's own calibration in place and nothing changes. Useful when the fit is one block of a larger sampler and the response it is handed varies in scale between sweeps. For the leaf-model qualifications on what the named quantity bounds, and for the \code{sd} spelling and its refusal under a sampled \code{k}, see \code{\link{dbarts}}'s Details and \code{\link{dbartsPriors}}.
   }
   \item{power}{
     Power parameter for tree prior. Default 2, Chipman, George, and McCulloch's empirical recommendation (see References) for the split-probability decay \eqn{base (1 + depth)^{-power}}{base * (1 + depth)^(-power)}, not derived from a formula.
   }
   \item{base}{
     Base parameter for tree prior. Default 0.95, from the same source as \code{power}.
   }
   \item{split.probs}{
     Prior and transition probabilities of variables used to generate splits. \code{NULL} (the default) requests equiprobability; can also be a numeric vector of length equal to the number of variables, or a named numeric vector with only a subset of the variables specified and a \code{.default} named value. Values given for factor variables are replicated for each resulting column in the generated model matrix. The symbol \code{num.vars} is rebound before execution to the number of columns in the model matrix. Collides with a supplied \code{tree.prior} (see below).
   }
   \item{dart}{
     \code{TRUE} places a DART prior (a Dirichlet prior over the split-variable probabilities, inducing variable selection; Linero 2018) on the trees with default settings, using \code{power} and \code{base} for the growth probabilities; a prior created by \code{\link[=dbartsPriors]{dbartsPriors$dart}} supplies its own settings, overriding \code{power} and \code{base}. Cannot be combined with \code{split.probs}. Per-sample probabilities are returned in the \code{varprobs} element of the fit.
   }
   \item{n.trees}{
     The number of trees in the sum-of-trees formulation.
   }
   \item{n.samples}{
     The number of posterior samples requested; \code{n.samples \%/\% n.thin} draws are actually kept and returned. This is a SWEEP budget in the BayesTree tradition - the same role \code{\link{bart}}'s \code{ndpost} plays (\code{ndpost / keepevery} are what is actually returned) and \code{\link{rbart_vi}}'s own \code{n.samples} plays too. It differs from \code{\link{dbartsControl}}'s (and the lower-level \code{\link{dbarts}}'s) \code{n.samples}, which already counts the draws one \code{run()} call returns and is NOT divided by thinning: that layer's sampler can be \code{run()} repeatedly with a different count each time, so its \code{n.samples} is a per-call return count rather than a one-shot total budget.
   }
   \item{n.burn}{
     Number of MCMC iterations to be treated as burn in.
   }
   \item{n.chains}{
     Integer specifying how many independent tree sets and fits should be calculated. Default 4; \code{\link{bart}}'s default is 1, BayesTree's historical single-chain convention.
   }
   \item{n.threads}{
     Integer specifying how many threads to use. Depending on the CPU architecture, using more than the number of chains can degrade performance for small/medium data sets. As such some calculations may be executed single threaded regardless.
   }
   \item{combineChains}{
     Logical; if \code{TRUE}, samples are returned collapsed across chains rather than in a chains-by-samples array. Default \code{TRUE}, as for \code{\link{bart}} and \code{\link{rbart_vi}}; the number of chains is tracked regardless on the returned object's \code{n.chains} component (see \sQuote{Value}).
   }
   \item{n.cuts}{
     The maximum number of possible values used in decision rules (see \code{useQuantiles}, and \code{\link{bart}}'s \sQuote{Decision Rules} details). If a single number, it is recycled for all variables; otherwise must be a vector of length equal to the number of predictor columns. Fewer rules may be used if a covariate lacks enough unique values.
   }
   \item{useQuantiles}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived from the predictor variables. When \code{FALSE}, splits are determined using values equally spaced across the range of a variable. See \code{\link{bart}}'s \sQuote{Decision Rules} details.
   }
   \item{n.thin}{
     Every \code{n.thin} draw is kept to be returned to the user. Useful for \dQuote{thinning} samples against serial correlation; see \code{n.samples} above for how the two combine.
   }
   \item{keepTrainingFits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of the training data are returned.
   }
   \item{printEvery}{
     As the MCMC runs, a message is printed every \code{printEvery} draws.
   }
   \item{printCutoffs}{
     The number of cutoff rules printed to screen before the MCMC is run. Given a single integer, the same value is used for all variables. If 0, nothing is printed.
   }
   \item{verbose}{
     Logical; if \code{FALSE} suppress printing.
   }
   \item{keepTrees}{
     Logical; must be \code{TRUE} in order to use \code{predict} with the result of a fit. Note that for models with a large number of observations or a large number of trees, keeping the trees can be very memory intensive.
   }
   \item{keepCall}{
     Logical; if \code{FALSE}, the returned object has \code{call} set to \code{call("NULL")}, otherwise the call used to instantiate BART.
   }
   \item{samplerOnly}{
     Builds the sampler from its arguments and returns it without running it. Useful to use the \code{bart2} interface in more complicated models.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. A \code{\link{set.seed}} beforehand suffices for reproducibility; supplying \code{seed} instead gives reproducible results without touching R's stream. See \code{\link{bart}}'s Reproducibility section.
   }
   \item{proposal.probs}{
     Named numeric vector, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step. The default is the named vector \code{c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)}, identical to \code{\link{dbarts}}'s own default.
   }
   \item{monotone}{
     Optional per-predictor monotonicity constraints, passed through to \code{\link{dbarts}}; see its \code{monotone} argument. A named vector selects predictors by model-matrix column name, each direction one of \code{"+"}/\code{"increasing"}/\code{1} or \code{"-"}/\code{"decreasing"}/\code{-1} - matching is case-insensitive, so \code{"Increasing"} is accepted; only numeric and ordered columns are eligible. A constraint forces birth/death-only proposals and a fixed \code{k = 2}. \code{NULL} (the default) fits the unconstrained model.
   }
   \item{interactions}{
     Optional interaction constraints, passed through to \code{\link{dbarts}}; see its \code{interactions} argument and \code{\link{interactions}}. Built with \code{interactions(max.order = , groups = , forbid = )}: cap the distinct predictors per root-to-leaf path, confine interactions to declared groups, and/or forbid named predictors from co-occurring. \code{NULL} (the default) fits the unconstrained model.
   }
   \item{blocks}{
     Optional block-additive constraint, passed through to \code{\link{dbarts}}; see its \code{blocks} argument and \code{\link{blocks}}. Built with \code{blocks(groups = , trees.per.group = )}: confine each whole tree to one declared group of predictors so the ensemble is exactly a sum of per-group functions. The groups must form a total, disjoint partition of the predictors, validated at fit time. \code{NULL} (the default) fits the ordinary model.
   }
   \item{variance}{
     Heteroscedastic BART (Pratola et al. 2020): passed through to \code{\link{dbarts}}. Declares a second forest modeling the residual variance surface \eqn{s^2(x)} as a product of scaled-inverse-chi-squared leaves, coupled to the mean forest through the per-observation precision. Either the plain selector naming which predictors drive the variance - a one-sided formula or column selector (\code{~ x1 + x2}), \code{TRUE} for every predictor, or \code{NULL} (the default) for a homoscedastic fit - or a \code{\link{varianceForest}(vars = , n.trees = , base = , power = )} object, which additionally sets the variance forest's own tree count and tree-structure prior (defaulting to \code{40} trees and the mean forest's \code{tree.prior}, exactly as the plain selector does). Gaussian responses only; \code{resid.dist = student()} residuals and a grouped (\code{\link{rbart_vi}}) fit are refused with it too, unadjudicated rather than unsupported by design (see \code{\link{dbarts}}). The fit gains \code{s.train}/\code{s.test} (posterior draws of \eqn{s(x)}), and \code{predict} attaches an \code{"s"} attribute carrying \eqn{s(x)} for new data. \eqn{s(x)} is the residual scale the reporting channels use: \code{extract(type = "loglik")}, the \code{type = "ppd"} draws, and \code{\link{summary.bart}}'s \code{mean.s} row all read it instead of the fixed \code{sigma} such a fit stores (see \code{\link{dbarts}}).
   }
   \item{keepSampler}{
     Logical that can be used to save the underlying \code{\link{dbartsSampler-class}} object even if \code{keepTrees} is false.
   }
   \item{warm.start}{
     A previous fit whose forests seed the initial state instead of drawing trees from the prior: either a \code{bart} object fit with \code{keepSampler = TRUE} or a \code{\link{dbartsSampler-class}}. The donor must share this fit's predictors, tree count, and DART setting; a donor fit on a different cut grid is supported, its splits remapped onto this fit's grid (collapsing any the grid starves), the same way a data replacement remaps existing splits. Only its trees, \code{sigma}, and \code{k} carry over, and each chain starts from a different donor sample so multiple chains stay overdispersed. A warm start biases early draws toward the donor, so it shortens burn-in rather than removing it; keep a non-zero \code{n.burn}. See \code{\link{dbartsSampler-class}}'s \code{installTrees} method for finer control.
   }
   \item{n.grow.sweeps}{
     A non-negative integer; when positive, the initial forest is built by \code{n.grow.sweeps} sweeps of XBART-style grow-from-root (He, Yalov and Hahn 2019) instead of a draw from the prior, an init-then-refine workflow that reaches a good fit in far fewer sweeps before the exact sampler takes over. The default \code{0L} draws trees from the prior, exactly as before. Like \code{warm.start}, a grow-from-root start biases the early draws toward the grown fit, so it shortens burn-in rather than removing it; keep a non-zero (if smaller) \code{n.burn}. The posterior the sampler targets is unchanged - only the starting point moves. Constant-leaf models only, and mutually exclusive with \code{warm.start} (both request an initialization). See \code{\link{dbartsSampler-class}}'s \code{growFromRoot} method.
   }
   \item{factors}{
     How factor columns in a data frame enter the model, exactly as in \code{\link{dbarts}}. The default \code{"categorical"} keeps each unordered factor as a single predictor whose splits send a subset of its levels down each branch - a different model-matrix representation than \code{\link{bart}}'s, which always expands every factor into indicator columns; \code{"indicators"} instead expands each factor into binary indicator columns, matching \code{bart}. A factor-predictor fit therefore changes if moved between the two interfaces.
   }
   \item{family}{
     The response family. \code{"auto"} (the default) fits gaussian models to continuous responses and probit models to those coded 0/1, exactly as \code{\link{dbarts}}; a two-level factor, logical, or two-level character response is also detected and fit as probit, an unordered factor (or character) response with three or more levels is detected and fit as multinomial, and an ordered factor with three or more levels is detected and fit as ordinal - each reporting the choice in a one-line message. \code{"gaussian"}, \code{"probit"}, and \code{"logistic"} force those fits directly; \code{"aft"}, \code{"hazard"}, \code{"hazard.probit"} (an accepted alias for \code{"hazard"}), \code{"hazard.logistic"}, \code{"multinomial"}, \code{"ordinal"}, \code{"nbinom"}, and \code{"hurdle.lognormal"} (alias \code{"twopart"}) reach the extended families described below. \code{\link{bart}} fits only a narrower, appended \code{c("auto", "logistic", "aft")} of this token set directly; naming one of the remaining tokens to \code{bart}'s \code{family} is refused, pointing here.

     \code{family = "multinomial"} fits a K-category softmax classifier: K forests, one per category, coupled through an interleaved Polya-Gamma one-vs-rest augmentation. When \code{formula} is a formula, \code{data} is required (not optional as for every other family): \code{bart2(y ~ x1 + x2, family = "multinomial", ...)} with no \code{data} errors by name rather than resolving \code{y}/\code{x1}/\code{x2} from the calling environment. \code{y.train} must be a factor (or character, coerced with \code{\link{factor}}) with at least two levels and no \code{NA}s or \code{NA} levels; \code{K = nlevels(y.train)} is inferred, never given explicitly, and \code{levels(y.train)} is captured at fit time and carried on every K-shaped output (\code{yhat.train}'s trailing dimension, \code{fitted}'s columns, the \code{fitted(type = "class")} factor). Alternatively, \code{y.train} may be an n x K numeric matrix of nonnegative integer trial counts, one row per observation and one column per category, with row sum (trials) \eqn{n_i \ge 1}{n_i >= 1}: category k's count is then a binomial(\eqn{n_i}{n_i}, p_k) draw under the same one-vs-rest augmentation, \code{K = ncol(y.train)} is inferred, and the levels carried on every K-shaped output are \code{colnames(y.train)} when present, otherwise \code{as.character(seq_len(K))}. A one-hot count matrix (every \eqn{n_i = 1}{n_i = 1}) is the same model as the factor response above, reproduced bit for bit. The formula interface is supported alongside the matrix one: \code{bart2(formula, data, family = "multinomial")} with a factor (or character) left-hand side routes to the factor response above, and a \code{cbind(c1, ..., cK) ~ x} left-hand side (the same idiom \code{glm}'s binomial family uses) routes to the count-matrix response above, \code{K} and the levels taken from the \code{cbind} column names; the right-hand side is coded exactly as it is for every other family's formula fit, including \code{predict} on a data frame \code{newdata}. Under the default \code{family = "auto"}, a factor (or character) response with three or more levels is detected and fit as multinomial, reporting the choice in a one-line message; a two-level factor (or logical) response resolves to probit, and a numeric response is unchanged. A count-matrix (\code{cbind(c1, ..., cK) ~ x}) response is only ever multinomial when \code{family = "multinomial"} is given explicitly - it is never inferred. \code{weights}, \code{subset}, \code{samplerOnly}, \code{warm.start}, \code{n.grow.sweeps}, \code{dart} (or a DART \code{tree.prior}), \code{split.probs}, \code{monotone}, and \code{variance} are all refused with an error naming the limitation (an integer weight is already expressible as row-wise count replication in the response, a non-integer one has no exact augmentation sampler, and the K-forest engine copies only power/base/proposal-probability fields from the host sampler it briefly builds, so DART, fixed split probabilities, monotone direction constraints, and a variance forest never reach it). \code{offset} is accepted only as an n x K numeric matrix, entering the K forests' raw fits before the softmax, one column per category, in the same layout as a count-matrix \code{y.train}; a flat (length-n) offset is refused by name, since a common per-observation shift is the softmax's own null direction and is identically inert. \code{offset} is a TRAIN-side argument only: \code{offset.test} is refused by name too, and \code{yhat.test} is always computed WITHOUT any category offset, even when \code{offset} was supplied for training - a caller comparing an offset-fitted \code{yhat.train} against \code{yhat.test} should keep this asymmetry in mind. A category test offset on the fit-time \code{test} rows is a sampler-level capability only (a \code{\link[=dbartsSampler-class]{dbartsSampler}}'s own \code{$setCategoryTestOffset} method, reached through \code{keepSampler = TRUE}, or the internal creators' own \code{offset.test} argument), not reachable from \code{bart2}; \code{predict}'s own \code{offset} is the supported route to an offset test surface, taking a matrix for the rows it is given. \code{test} is supported: an \code{x.test} of the same column structure as \code{x.train} reports the K-category softmax probabilities on the held-out rows as \code{yhat.test}, shaped and levels-named exactly like \code{yhat.train} (see \sQuote{Value}). \code{keepTrees} is supported too: it retains every one of the K forests' trees so \code{predict} can replay them at new predictors afterward, reproducing \code{yhat.test} bitwise when \code{newdata} matches the fit-time \code{test}; without \code{keepTrees}, \code{predict} errors. A fit trained with an \code{offset} replays as well, given \code{predict}'s \code{offset} at the new rows, and refuses by name without it: no resident offset describes rows the fit never saw. The per-forest leaf scale follows its own K-dependent calibration (the K = 2 anchor is the logistic scale \eqn{\pi\sqrt{3}}{pi*sqrt(3)} divided by \eqn{\sqrt{2}}{sqrt(2)}, for the identified pairwise log-odds); \code{k} is read from the usual node prior exactly as for any other family, but the node prior's \code{node.scale} itself is NOT consulted - the multinomial engine calibrates its own. The fit's class is \code{"bartMultinomial"}, not \code{"bart"}: see \sQuote{Value} below and the \code{extract}/\code{fitted}/\code{predict} methods for \code{bartMultinomial} objects.

     \code{family = "ordinal"} fits an ordered categorical response by a cumulative probit (a single forest, unlike multinomial's K): a latent \eqn{z = f(x) + \epsilon}{z = f(x) + eps} with \eqn{\epsilon \sim N(0, 1)}{eps ~ N(0, 1)} is cut at ordered thresholds \eqn{-\infty = \gamma_0 < \gamma_1 = 0 < \gamma_2 < \ldots < \gamma_{K-1} < \gamma_K = \infty}{-Inf = gamma_0 < gamma_1 = 0 < gamma_2 < ... < gamma_{K-1} < gamma_K = Inf}, so \eqn{P(Y = k \mid x) = \Phi(\gamma_k - f(x)) - \Phi(\gamma_{k-1} - f(x))}{P(Y = k | x) = Phi(gamma_k - f(x)) - Phi(gamma_{k-1} - f(x))}. The free cutpoints are sampled by a one-at-a-time marginal Metropolis update (latents integrated out) under a log-gap prior \eqn{N(0, 1.5^2)}{N(0, 1.5^2)}; at \eqn{K = 2} the model is exactly the probit family. \code{y.train} should be an ordered factor (\code{\link{is.ordered}}); its level order defines the category order, and \code{levels(y.train)} is captured at fit time and carried on every K-shaped output. Under the default \code{family = "auto"}, an ordered factor with three or more levels is detected and fit as ordinal, reporting the choice in a one-line message (an unordered factor stays multinomial - the two never overlap, the dispatch key is \code{is.ordered}); a two-level ordered factor is binary and resolves to probit. Given explicitly, \code{family = "ordinal"} also accepts an unordered factor or character response (with a message that the category order is taken from the level order - usually alphabetical, so order the levels yourself) and a numeric response, whose ordered levels are \code{sort(unique(y.train))}. \code{weights} are not supported (a weighted truncated-normal latent likelihood is not a coherent model); \code{warm.start} and \code{n.grow.sweeps} are refused with an error naming the limitation. \code{samplerOnly} is supported, unlike multinomial: it returns the fitted sampler unrun, the same object \code{fit} is under \code{keepSampler}. \code{test} and \code{keepTrees} are supported as for multinomial, and \code{predict} replays the saved trees at new predictors, differencing the cumulative probit at the stored per-draw cutpoints. \code{rbart_vi} and \code{xbart} do not fit ordinal responses (grouped ordinal models and ordinal-scale cross-validation losses are recorded follow-ups). The fit's class is \code{"bartOrdinal"}, not \code{"bart"}: see \sQuote{Value} below.

     \code{family = "nbinom"} fits a non-negative integer (count) response by a negative-binomial model with the Polya-Gamma augmentation (a single forest, like logistic): the forest fits a log-odds latent \eqn{\psi = f(x) + o}{psi = f(x) + o}, the count law is \eqn{y_i \sim \mathrm{NB}(r, \mathrm{plogis}(\psi_i))}{y_i ~ NB(r, plogis(psi_i))} with dispersion \eqn{r}, and the mean is \eqn{E[y_i \mid x] = r e^{\psi_i} = r e^{f(x_i) + o_i}}{E[y_i | x] = r exp(psi_i) = r exp(f(x_i) + o_i)}. The offset \eqn{o_i} therefore enters the mean multiplicatively and is a log-exposure (\eqn{o_i = \log(\mathrm{exposure}_i)}{o_i = log(exposure_i)}). \code{y.train} must be a non-negative integer count, and \code{family = "nbinom"} is always explicit - a count carries no unambiguous class, so it is never inferred under \code{family = "auto"} (a numeric response there stays gaussian). The dispersion \code{r} is estimated by default on a capped positive-integer grid (\eqn{\{1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50\}}{{1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50}}) under a renormalized \eqn{\mathrm{gamma}(2, 0.1)}{gamma(2, 0.1)} prior by a closed-form discrete full conditional; passing \code{dispersion} as a positive integer fixes \eqn{r} instead (v1 ships the exact integer envelope, so a real fixed dispersion is refused - continuous dispersion is a recorded follow-up). Like probit, the latent scale is fixed at 1 and \code{weights} are not supported (exposure belongs in the offset, not in observation replication); \code{warm.start} and \code{n.grow.sweeps} are refused with an error naming the limitation. \code{samplerOnly} is supported, unlike multinomial: it returns the fitted sampler unrun, the same object \code{fit} is under \code{keepSampler}. \code{test} and \code{keepTrees} are supported as for ordinal, and \code{predict} replays the saved trees at new predictors, reporting mean counts \eqn{r e^{\psi}}{r exp(psi)} at the stored per-draw dispersion (a log-exposure \code{offset.test} enters \eqn{\psi}{psi} additively). \code{rbart_vi} and \code{xbart} do not fit count responses (grouped negative-binomial models, real dispersion, and a Poisson family are recorded follow-ups). \code{y.train} is additionally capped at \eqn{10^6}{1e6}: the dispersion grid's count histogram is sized from the largest count, so a larger one allocates without bound, and a count above it is refused at creation and at every response swap alike. The fit's class is \code{"bartNegbin"}, not \code{"bart"}: see \sQuote{Value} below.

     \code{family = "hazard"} and \code{family = "hazard.logistic"} fit a discrete-time survival hazard model as \emph{ingestion sugar} over the binary families - the model adds no engine code (\code{"hazard.probit"} is an accepted alias for \code{"hazard"}). The response is the same \code{Surv} object or two-column \code{(time, status)} the \code{"aft"} family takes; each subject observed to time \eqn{t_i} is person-period-expanded into its at-risk rows (one per period \eqn{k = 1, \ldots, t_i}), each carrying the subject's covariates, an ordinal period column (appended last, so trees split on contiguous periods and borrow strength across time), and the binary indicator \eqn{y_{ik} = \mathrm{status}_i \cdot 1\{k = t_i\}}. The discrete hazard is \eqn{h(k \mid x) = g(f(x, k) + o)} for the chosen binary link \eqn{g} (probit for \code{"hazard"}, logistic for \code{"hazard.logistic"}), and the fit is EXACTLY an ordinary \code{"probit"}/\code{"logistic"} fit on the expanded rows - a hazard fit is byte-identical, draw for draw, to the binary fit on the hand-expanded design with the same seed. The time grid defaults to the sorted distinct observed times (the \pkg{BART} \code{surv.bart} convention); \code{breaks} coarsens it (a single integer bins at that many quantiles, a boundary vector gives explicit right-closed intervals), and \code{max.rows} guards the expansion size. Offsets are on the link scale and replicate per subject; \code{weights} replicate and then follow the chosen binary family's policy (probit refuses non-unit weights, logistic requires integer counts). A \code{Surv} response requires the family to be requested explicitly (\code{family = "auto"} selects \code{"aft"}); the matrix interface is required (no formula, \code{subset}, or \code{test} - expand held-out subjects with \code{survivalProbabilities(fit, times, newdata = )}, which needs \code{keepTrees}). The returned object is an ordinary \code{"bart"} fit whose \code{$family} records the binary link (so every prediction generic transforms correctly with no special case) and whose \code{$periods} element carries the grid; \code{predict}/\code{extract} return the per-(subject, period) hazards on the expanded rows, and \code{\link{survivalProbabilities}} produces survival-curve draws \eqn{S(t \mid x) = \prod_{k : \mathrm{periods}[k] \le t} (1 - h(k \mid x))}. \code{rbart_vi} and \code{xbart} do not fit hazard responses (grouped/frailty hazard and a \code{cloglog} link are recorded follow-ups).

     \code{family = "hurdle.lognormal"} (alias \code{"twopart"}, which resolves and prints as \code{"hurdle.lognormal"}) fits a semicontinuous two-part (hurdle) model for a non-negative response with exact zeros: an OCCUPANCY probit fit of \eqn{z = 1\{y > 0\}}{z = 1{y > 0}} over all n observations, glued at report time to a POSITIVE-PART gaussian fit of \eqn{\log y}{log y} over the subset \eqn{\{i : y_i > 0\}}{{i : y_i > 0}}. The two component fits share no parameters and are composed entirely R-side from two ordinary \code{bart2} fits at independently derived seeds - no engine code, and no coupling between the parts - so a shared variable-selection prior across the occupancy and positive parts is foreclosed by this composition (a recorded limitation, alongside a Duan-smearing retransformation and a heteroscedastic positive part, both follow-ups). \code{y.train} must be non-negative and finite and must carry at least one exact zero and one positive value (a response with no zeros, or none positive, is refused by name). By default \code{predict}/\code{fitted}/\code{extract} report the NATURAL (response) scale via posterior-predictive Monte Carlo: \eqn{E[y \mid x] = P(y > 0 \mid x)\,e^{f(x) + \sigma^2 / 2}}{E[y | x] = Pr(y > 0 | x) * exp(f(x) + sigma^2 / 2)} computed per posterior draw from the positive part's single \eqn{\sigma} per draw, recycled across observations (the positive part is always homoscedastic; a heteroscedastic positive part is not reachable, see the recorded limitation above). \code{type = "prob"} returns the occupancy probability \eqn{\pi(x)} through the probit link; \code{type = "link"}/\code{"log"} the positive part's log-scale linear predictor \eqn{f(x)}; \code{type = "ppd"} draws the proper BIMODAL predictive - a Bernoulli(\eqn{\pi}) spike at zero, else a lognormal draw - which the plain gaussian ppd path cannot produce. Fits currently use the matrix interface only (\code{bart2(x.train, y.train, family = "hurdle.lognormal")}); \code{weights}, \code{subset}, \code{offset}/\code{offset.test}, and \code{test} are all refused with an error naming the limitation - the positive-part fit is instead given the full training \code{x} as its own \code{x.test}, so its in-sample fitted values cover the zero rows it never trained on. \code{predict} requires \code{keepTrees = TRUE} and replays both saved forests at \code{newdata}. \code{dbarts()} does not fit this family - it composes two samplers, which only \code{bart2()} builds, so requesting it there is an error directing here. \code{rbart_vi} and \code{xbart} do not fit it either. The fit's class is \code{"bartHurdle"}, holding both component fits (\code{$occupancy}, a \code{"bart"} probit fit of the occupancy indicator; \code{$positive}, a \code{"bart"} gaussian fit of \eqn{\log y}{log y} on the positive subset) under their own \code{extract}/\code{fitted}/\code{predict}/\code{residuals}/\code{print} methods.
   }
   \item{missing}{
     How missing predictor values enter the model, exactly as in \code{\link{dbarts}}. The default \code{"incorporate"} keeps them, every split rule learning a direction for \code{NA}; \code{"error"} rejects predictors containing \code{NA} - the historical behavior \code{\link{bart}} always uses.
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
   \item{tree.prior, node.prior, resid.prior}{
     The full \code{\link{dbarts}} prior objects, forwarded unevaluated so a bare \code{\link{dbartsPriors}} vocabulary name inside them (\code{cgm}, \code{dart}, \code{normal}, \code{linear}, \code{gp}, \code{chisq}, \code{fixed}) resolves exactly as it does for \code{dbarts}: \code{node.prior = linear(columns = 1:3)}, \code{gp(...)}, and \code{resid.prior = fixed(1)} are reachable this way. \code{NULL} (the default for all three) instead builds the tree/node/residual priors from \code{power}/\code{base}/\code{split.probs}/\code{dart}, \code{k}/\code{prior.scale}, and \code{sigdf}/\code{sigquant} respectively, exactly as before this argument existed. Supplying an object alongside a shorthand that would otherwise help build the same prior is an error naming both: \code{tree.prior} collides with any of \code{power}/\code{base}/\code{split.probs}/\code{dart}; \code{node.prior} with \code{k}/\code{prior.scale}; \code{resid.prior} with \code{sigdf}/\code{sigquant}/\code{sigest}. A supplied \code{resid.prior} under a fixed-unit-scale family (\code{"probit"}, \code{"logistic"}, \code{"ordinal"}, \code{"nbinom"}, the hazard families, \code{"multinomial"}) is diagnosed the same way \code{sigest}/\code{sigdf}/\code{sigquant} already are - a \code{dbartsFamilyGatedWarning}, since the family overwrites it with \code{fixed(1)} regardless. \code{tree.prior}/\code{node.prior} are honored on every family, including both component fits of \code{family = "hurdle.lognormal"}.
   }
   \item{storage}{
     Passed through to \code{\link{dbartsControl}}; a formal on both \code{bart} and \code{bart2}. A character string selecting the precision of the internal running residual; see \code{dbartsControl}'s \code{storage} for the full description. The default \code{"double"} reproduces existing draws bitwise; \code{"single"} changes the numbers a fit returns - the sampled values shift slightly - in exchange for a bandwidth-bound speedup, and is currently supported only for continuous (gaussian) responses with constant leaves.
   }
   \item{updateState}{
     Passed through to \code{\link{dbartsControl}}. A logical setting the default behavior for many \code{\link[=dbartsSampler-class]{sampler}} methods with regards to the immediate updating of the cached state of the object; see \code{dbartsControl}'s \code{updateState}.
   }
   \item{\dots}{
     Unused; present on the \code{bartMultinomial}/\code{bartOrdinal}/\code{bartNegbin}/\code{bartHurdle} methods below only for S3 generic compatibility.
   }
   \item{object}{
     An object of class \code{bartMultinomial}, \code{bartOrdinal}, \code{bartNegbin}, or \code{bartHurdle}, returned by \code{bart2} under the corresponding \code{family}.
   }
   \item{newdata}{
     Test data for prediction. Obeys the same rules as \code{data}/\code{test} but cannot be missing.
   }
   \item{type}{
     The quantity to return; see \sQuote{Value} below for what each family's generics accept and produce. As on a \code{"bart"} fit (\code{\link{bart}}'s \code{type} item), \code{"response"} is a synonym for \code{"ev"} and \code{"link"} for \code{"bart"}, on every method here; a value a method does not offer is an error naming the ones it does, and one it names only to refuse is an error naming the reason.
   }
   \item{sample}{
     Either \code{"train"} or \code{"test"}.
   }
   \item{vars}{
     Character vector of fields to gather, as \code{\link{summary.bart}}'s own \code{vars}; requested fields absent from \code{object} are silently dropped. \code{summary.bartOrdinal}'s \code{"cutpoints"} contributes one draws variable per threshold \eqn{\gamma_1, \ldots, \gamma_{K-1}}{gamma_1, ..., gamma_{K-1}}, the ordinal analog of \code{sigma}; \code{summary.bartNegbin}'s \code{"dispersion"} contributes the per-draw dispersion \eqn{r}, the count analog of \code{sigma}. \code{summary.bartHurdle} applies \code{vars} separately to its \code{$occupancy} and \code{$positive} component fits. \code{summary.bartMultinomial} has no \code{vars} formal at all: it always pools the per-category mean-probability channel (its only scalar convergence instrument), so a supplied \code{vars} is refused by name rather than silently ignored.
   }
   \item{x}{
     An object of one of the four own-class families (\code{bartMultinomial}, \code{bartOrdinal}, \code{bartNegbin}, \code{bartHurdle}), as returned by \code{bart2} under the matching \code{family}, for \code{print}/\code{plot}; an object of class \code{summary.bartHurdle}, as returned by \code{summary} on a \code{bartHurdle} fit, for \code{print.summary.bartHurdle}.
   }
   \item{cols}{
     For \code{plot.bartMultinomial}/\code{plot.bartOrdinal}: a vector of colors, one per category, tracing each category's training-mean predicted probability (\code{bartOrdinal}'s cutpoint panel, when drawn) or each category's own P2 panel point; \code{NULL} (the default) chooses colors automatically. For \code{plot.bartNegbin}/\code{plot.bartHurdle}: two colors, as \code{\link{bart}}'s own \code{cols} item - the first for a posterior interval's bars, the second for the reference line.
   }
   \item{plquants}{
     For \code{plot.bartOrdinal}/\code{plot.bartNegbin}/\code{plot.bartHurdle}/\code{plot.bartMultinomial}, as \code{\link{bart}}'s own \code{plquants} item: the lower and upper quantiles bounding each posterior-interval panel's bars.
   }
   \item{ci.level}{
     For \code{fitted}/\code{predict} on a \code{bartMultinomial}, \code{bartOrdinal}, or \code{bartNegbin} fit, as \code{\link{bart}}'s own \code{ci.level} item. On a K-margined \code{type} (\code{"ev"}/\code{"ppd"} for the first two), the result keeps the category margin: an (observation \eqn{\times} K \eqn{\times} 3) array of \code{est}/\code{ci.lower}/\code{ci.upper} rather than a 3-column matrix, since categories are not poolable. On a non-K-margined \code{type} (\code{"bart"}; any type for \code{bartNegbin}) the result is the ordinary 3-column matrix. \code{bartHurdle}'s \code{ci.level} is documented under \code{\link{bart}} directly, since it never carries a K margin.
   }
}
\section{Formula Terms}{
A \code{forest()} call written inside \code{bart2}'s \code{formula} - bare, or as one operand of \code{:} - declares a second forest whose amplitudes multiply a basis, without a separate \code{\link{dbarts}} \code{forests =} list (the two routes are mutually exclusive: a term-bearing formula given together with \code{forests =} is refused by name). The canonical spelling uses colon sugar:

\preformatted{
y ~ x1 + x2 + z:forest(x1 + x2)             # numeric z: one amplitude
y ~ x1 + x2 + factor(z):forest(x1 + x2)     # factor z: one amplitude per level
y ~ x1 + x2 + (a + b):forest(x1 + x2)       # two-column basis, one forest
}

Each desugars to the general named form, \code{forest(x1 + x2, basis = ~ z)} (respectively \code{~ factor(z)}, \code{~ cbind(a, b)}) - textually identical to a \code{\link{forest}} object built directly for \code{\link{dbarts}}'s \code{forests =}, so the term route and the constructor route share one vocabulary and one validator; write \code{basis =} explicitly for anything the colon sugar cannot say (an arbitrary formula, or an already-evaluated vector/matrix). \code{z * forest(...)} is refused, naming both explicit spellings (\code{z:forest(...)}, or \code{z + z:forest(...)} to also keep \code{z} as a main effect) instead of guessing which one was meant.

In this position \code{forest()}'s unnamed first slot is a SYMBOLIC predictor set, not \code{basis}: \code{forest(x1 + x2)} means the variable set \code{c("x1", "x2")}, read from the call and never evaluated - a bare symbol, or symbols joined by \code{+}, and nothing else; a transformation, a removal, \code{.}, or \code{:}/\code{*} inside the slot is refused by name, pointing at \code{vars = } for what the grammar cannot say. A factor name in that slot expands to every column the factor produces (all of its indicator columns under \code{factors = "indicators"}). Every NAMED argument inside a term-context \code{forest()} - including \code{basis} itself, for the general form \code{forest(x1 + x2, basis = ~ z)} - evaluates normally, exactly as it does for \code{\link{dbarts}}'s own \code{forests =}.

A term's basis is evaluated inside the fit's own model frame, so \code{subset}, missing-data handling, and \code{data} scoping apply to it exactly as they do to any predictor - unlike a basis given directly through \code{forests =}, which is evaluated against the raw \code{data} and then restricted to the subsetted rows (see \code{\link{forest}}'s \code{basis} item for that route's rule; the two reach the same rows whenever \code{subset} keeps every row unaltered).

A handful of shapes are refused outright rather than silently misinterpreted, each error naming the offending expression: a \code{forest()} term anywhere other than a top-level additive term or a \code{:} operand (inside \code{I()}, a removal, or the formula's left-hand side); a compound \code{:} operand containing anything other than bare numeric/logical symbols joined by \code{+} (a factor or character member, a transformation, a multi-way \code{:} chain, or a second \code{forest()} on the other side); \code{test} given together with a term (an amplitude-coupled fit has no per-observation test replay in this version); a family the multiplier model does not support (\code{"aft"}, \code{"ordinal"}, \code{"nbinom"}, the hazard families, \code{"multinomial"}, and \code{"hurdle.lognormal"}/\code{"twopart"} - gaussian, probit, and logistic are the only families a term can join); and a formula whose only right-hand-side content is the term, which would leave nothing to split on. See \sQuote{Value} below, and \code{\link{extract}}'s \code{type = "forest"}, for reading the resulting per-forest fits back out.
}
\value{
Under the default and \code{\link{bart}}-compatible families (\code{"gaussian"}, \code{"probit"}, \code{"logistic"}, \code{"aft"}, \code{"hazard"} and its aliases), \code{bart2} returns a list of class \code{"bart"}, documented in full under \code{\link{bart}}'s \sQuote{Value} section - the same \code{yhat.train}/\code{sigma}/\code{varcount}/\dots{} components, including \code{forestFits}, \code{glue}, \code{bases}, \code{n.forests}, and a \code{"forest.labels"} attribute when the formula carries a \code{forest()} term (\sQuote{Formula Terms} above): that term route is the only way \code{bart2} itself reaches the per-forest channel, since it has no \code{forests =} formal of its own. A heteroscedastic fit (\code{variance} supplied) additionally carries \code{s.train}/\code{s.test}, posterior draws of the residual scale \eqn{s(x)} laid out as \code{yhat.train}/\code{yhat.test}; see the \code{variance} argument above. A discrete-time hazard fit (\code{family = "hazard"}/\code{"hazard.logistic"}) additionally carries \code{periods}, the ordered period grid the person-period expansion used; \code{\link{survivalProbabilities}} dispatches its survival-curve branch on it.

The remaining families - \code{"multinomial"}, \code{"ordinal"}, \code{"nbinom"}, and \code{"hurdle.lognormal"}/\code{"twopart"} - return their own list class instead, never \code{"bart"}, each with its own generics, described below: precisely so a K-category or two-part fit cannot silently flow through a method that assumes the single-forest shape.

\code{bart2(family = "multinomial")} returns a \emph{different} list, of class \code{"bartMultinomial"}. Components: \code{call}, \code{family} (\code{"multinomial"}), \code{levels} (the resolved category names, length K: \code{levels(y.train)} for a factor response, or \code{colnames(y.train)} - falling back to \code{as.character(seq_len(K))} when unnamed - for a count-matrix response), \code{K}, \code{n.chains}, \code{n.trees}, \code{y} (the original response: the factor \code{y.train}, or the validated n x K count matrix when a count response was supplied), and \code{yhat.train} - the posterior draws of the K softmax probabilities, an array of dimension (\code{n.chains} \eqn{\times}, when \code{n.chains > 1} and \code{combineChains = FALSE}) \code{n.samples} \eqn{\times} number of training observations \eqn{\times} K, with the resolved category levels as the trailing dimension's names; \code{combineChains = TRUE} (the default) folds the chain margin into the samples margin as usual. \code{yhat.train} already holds PROBABILITIES, not a latent score - there is no \code{"bart"}-scale component to convert, unlike the binary families. When \code{test} was supplied, \code{yhat.test} is the same K softmax probabilities on the held-out rows, dimensioned like \code{yhat.train} with the training-observation margin replaced by the test one. \code{varcount} is the per-sample per-category split-usage channel: each category forest's per-draw variable counts, dimensioned like \code{yhat.train} but with the number of training predictors in place of the number of observations, and \code{colnames(x.train)} (when present) named on that margin. \code{fit} is present whenever \code{keepTrees} is \code{TRUE} \emph{or} \code{keepSampler} is set, independent of \code{keepTrees}: it is the \code{dbartsSampler} whose K-forest engine actually ran (one \code{bartcore_create}, not a discarded host), fully mutable and readable on the channels the softmax gives meaning to - \code{$setCounts}, \code{$setCategoryOffset}, \code{$setPredictor}, and the rest - and refused by name on the ones it does not (\code{$setResponse}, \code{$setOffset}, \code{$setSigma}, \code{$setCalibration}, \code{$setForestWeights}); \code{fit$storeState()} followed by \code{\link{save}}/\code{\link{load}} restores a sampler \code{predict.bartMultinomial} can replay through. \code{predict} still requires \code{keepTrees = TRUE} - a kept \code{fit} alone carries no saved trees to replay - and errors otherwise.

Generics for a \code{"bartMultinomial"} fit: \code{fitted(object)} (\code{type = "ev"}, the default) returns the posterior-mean n \eqn{\times} K probability matrix, columns named by the fit's captured \code{levels}; \code{fitted(object, type = "class")} returns the argmax of that mean matrix as a factor over the original levels - a classification convenience. \code{extract(object, type = "ev", sample = "train")} returns the full \code{yhat.train} draws, and \code{sample = "test"} returns \code{yhat.test} (an error if the fit carries no test channel); \code{extract(object, type = "ppd")} draws one category per posterior draw from its probability vector, returned as an integer code (1-based, indexing the fit's captured \code{levels}) in an array shaped like \code{"ev"} minus the K margin. \code{extract}/\code{fitted}/\code{predict} with \code{type = "bart"}, and \code{extract}/\code{predict} with \code{type = "forest"}, error naming the reason: the run records only the identified softmax probabilities, and a category's forest is a latent whose level is reproducibly structured yet not identified (in the same sense as BCF's \eqn{a}), so a raw replay would read as signal - the identified content is the log-ratio, which the logs of the reported probabilities carry. \code{predict(object, newdata)} requires \code{object} fit with \code{keepTrees = TRUE} (error otherwise); it codes \code{newdata} to the training columns exactly as \code{predict.bart} does and returns a levels-named (\code{n.chains} \eqn{\times}) \code{n.samples} \eqn{\times} number of new observations \eqn{\times} K probability array, the \code{yhat.test}/\code{yhat.train} convention; \code{type = "ppd"} draws one category per posterior draw the same way \code{extract(object, type = "ppd")} does; \code{type = "bart"} and \code{type = "forest"} are named only to be refused, for the non-identification reason above. \code{offset} supplies the per-category shift at the new rows, required on a fit trained with an \code{offset} and reproducing \code{yhat.train} when the training offset is passed back at the training rows. \code{residuals(object)} returns an n \eqn{\times} K matrix, the observed proportion (an indicator for a factor response, \code{y / rowSums(y)} for a count-matrix one) minus the fitted probability in \code{fitted(object)}, per category. \code{extract}'s \code{combineChains} (default \code{TRUE}) reshapes \code{"ev"}/\code{"ppd"}/\code{"loglik"} to a combined or per-chain layout, exactly as \code{extract.bart}'s does. \code{fitted}/\code{predict}'s \code{ci.level} returns an (observation \eqn{\times} K \eqn{\times} 3) array of \code{est}/\code{ci.lower}/\code{ci.upper} instead of a 3-column matrix, since a category margin survives alongside the observation one. \code{extract(object, type = "loglik")} is the multinomial log-density of the observed row (its coefficient included, so it reduces to \eqn{\log p_{i, y_i}}{log p[i, y_i]} for a one-trial row), extract-only, \code{sample = "test"} refused by name, shaped like \code{"ppd"} (the K margin dropped) - loo/WAIC on it is leave-one-\emph{row}-out, since the likelihood unit is the whole count row. \code{forest}/\code{contribution} on \code{extract}/\code{predict}, \code{sample} on \code{fitted}, \code{type} on \code{residuals}, and \code{vars} on \code{summary} are all refused by name rather than silently ignored: none is meaningful on this K-widened, non-identified-latent shape. \code{plot(object)} adds a second panel to the per-category trace: the posterior median and interval of the predicted probability of each observation's OWN observed category against that median (a count response with multi-trial rows instead plots the observed proportion \eqn{y_{ik} / n_i}{y_ik / n_i} against the interval of the drawn \eqn{p_{ik}}{p_ik}, one point per row-category cell). \code{plotTree} and \code{\link{survivalProbabilities}} are refused by name (a multinomial fit's trees live on its sampler; it has no hazard channel). \code{as_draws_array}/\code{as_draws_df} (\pkg{posterior}, if installed) expose the same per-category mean-probability channel \code{summary(object)} pools, named \code{meanProb[<level>]} - this family has no other scalar posterior parameter. \code{summary(object)} pools that same per-category mean-probability channel into posterior mean/sd/quantiles (R-hat/ESS when the \pkg{posterior} package is installed), the multinomial analog of \code{summary.bart}'s \eqn{\sigma}/k/\eqn{\tau} summary.

\code{bart2(family = "ordinal")} likewise returns its own list, of class \code{"bartOrdinal"}. Components: \code{call}, \code{family} (\code{"ordinal"}), \code{levels} (the ordered category names, length K), \code{K}, \code{n.chains}, \code{n.trees}, \code{y} (the observed categories as an ordered factor over \code{levels}), \code{yhat.train} (and \code{yhat.test} when \code{test} was supplied) - the posterior draws of the K category probabilities, shaped and levels-named exactly as the multinomial arrays above - \code{latent.train} (and \code{latent.test}) - the corresponding draws of the latent \eqn{f(x)} on the unit-variance probit scale, shaped like a binary family's \code{yhat.train} - \code{cutpoints} - the posterior draws of the K - 1 finite thresholds \eqn{(\gamma_1 = 0, \gamma_2, \ldots)}{(gamma_1 = 0, gamma_2, ...)}, dimensioned draws \eqn{\times}{x} (K - 1), the ordinal analog of gaussian's \code{sigma}, from which probabilities at any latent value can be reconstructed - and \code{varcount}. \code{cutpoints.raw} (the per-draw thresholds in the internal layout \code{predict} consumes) is present only under \code{keepTrees}. \code{fit} is present whenever \code{keepTrees} is \code{TRUE} \emph{or} \code{keepSampler} is set, independent of \code{keepTrees}: it is the \code{dbartsSampler} whose engine actually ran, fully mutable and readable like any sampler \code{\link{dbarts}} returns, and \code{fit$storeState()} followed by \code{\link{save}}/\code{\link{load}} restores a sampler \code{predict.bartOrdinal} can replay through. \code{summary(object)} reports the cutpoints alongside whatever mean-function scale \code{vars} finds present, pooled into posterior mean/sd/quantiles (R-hat/ESS when \pkg{posterior} is installed), the ordinal analog of \code{summary.bart}'s \eqn{\sigma}/k/\eqn{\tau} summary.

Generics for a \code{"bartOrdinal"} fit: \code{fitted(object)} returns the posterior-mean n \eqn{\times}{x} K probability matrix (columns named by \code{levels}); \code{fitted(object, type = "class")} the argmax category of that mean as an ordered factor; \code{fitted(object, type = "bart")} the posterior-mean latent. \code{extract(object, type = "ev")} returns the probability draws, \code{type = "bart"} (alias \code{"link"}) the latent draws, and \code{type = "ppd"} draws one category per posterior draw (1-based codes indexing \code{levels}). \code{predict(object, newdata)} requires \code{keepTrees = TRUE} and returns the probability draws at \code{newdata} (\code{type = "bart"} the replayed latent; \code{type = "ppd"} category draws); when \code{newdata} matches the fit-time \code{test} it reproduces \code{yhat.test}. \code{residuals(object)} returns the n \eqn{\times}{x} K matrix of observed-indicator minus fitted probability. \code{extract}'s \code{combineChains} (default \code{TRUE}) and \code{fitted}/\code{predict}'s \code{ci.level} work as for \code{bartMultinomial} above, including the (observation \eqn{\times} K \eqn{\times} 3) interval shape for a K-margined \code{type}. \code{extract(object, type = "loglik")} is \eqn{\log P(y_i \mid \eta, \gamma)}{log P(y_i | eta, gamma)} - exactly the stored category probability at the observed level, since \code{yhat.train} already holds the cumulative-probit difference - extract-only, \code{sample = "test"} refused by name, shaped like \code{type = "ppd"}. \code{forest}/\code{contribution} on \code{extract}/\code{predict} and \code{sample} on \code{fitted} are refused by name (this family has a single forest and \code{fitted} is always the training rows). \code{plot(object)} traces the free cutpoints \eqn{\gamma_2, \ldots, \gamma_{K-1}}{gamma_2, ..., gamma_{K-1}} (\eqn{\gamma_1 = 0}{gamma_1 = 0} is pinned and carries no information; at K = 2 there is none, and the plot degrades to one full-device panel) alongside the per-observation latent interval, observations ordered by median \eqn{\eta}{eta}, coloured by observed level, with dashed reference lines at the posterior-median cutpoints. \code{plotTree} and \code{\link{survivalProbabilities}} are refused by name. \code{as_draws_array}/\code{as_draws_df} default to \code{vars = c("cutpoints", "sigma", "k", "tau")}, matching \code{summary}'s own default; \code{cutpoints} contributes \code{cutpoint[1]} (the pinned 0) through \code{cutpoint[K - 1]}.

\code{bart2(family = "nbinom")} likewise returns its own list, of class \code{"bartNegbin"}. Components: \code{call}, \code{family} (\code{"nbinom"}), \code{n.chains}, \code{n.trees}, \code{y} (the observed counts), \code{yhat.train} (and \code{yhat.test} when \code{test} was supplied) - the posterior draws of the mean counts \eqn{\mu = r e^{f(x) + o}}{mu = r exp(f(x) + o)}, shaped like a binary family's \code{yhat.train} - \code{latent.train} (and \code{latent.test}) - the corresponding draws of the log-odds latent \eqn{\psi = f(x) + o}{psi = f(x) + o} - \code{dispersion} - the per-draw dispersion \eqn{r}, the count analog of gaussian's \code{sigma} - and \code{varcount}. \code{dispersion.raw} (the per-draw \eqn{r} in the internal layout \code{predict} consumes) is present only under \code{keepTrees}. \code{fit} is present whenever \code{keepTrees} is \code{TRUE} \emph{or} \code{keepSampler} is set, independent of \code{keepTrees}: it is the \code{dbartsSampler} whose engine actually ran, fully mutable and readable - \code{$getDispersion()} on it answers with the fit's own \eqn{r} - and \code{fit$storeState()} followed by \code{\link{save}}/\code{\link{load}} restores a sampler \code{predict.bartNegbin} can replay through. \code{summary(object)} reports the dispersion draws alongside whatever mean-function scale \code{vars} finds present, pooled into posterior mean/sd/quantiles (R-hat/ESS when \pkg{posterior} is installed), the count analog of \code{summary.bart}'s \eqn{\sigma}/k/\eqn{\tau} summary.

Generics for a \code{"bartNegbin"} fit: \code{fitted(object)} returns the posterior-mean count per observation; \code{fitted(object, type = "bart")} the posterior-mean log-odds latent. \code{extract(object, type = "ev")} returns the mean-count draws, \code{type = "bart"} the latent draws, and \code{type = "ppd"} draws one count per posterior draw from \eqn{\mathrm{NB}(r, \mathrm{plogis}(\psi))}{NB(r, plogis(psi))}; \code{sample = "test"} selects the test channel. \code{predict(object, newdata)} requires \code{keepTrees = TRUE} and returns the mean-count draws at \code{newdata} (\code{type = "bart"} the replayed latent; \code{type = "ppd"} count draws), with an optional log-exposure \code{offset.test}; when \code{newdata} matches the fit-time \code{test} and no offset is given it reproduces \code{yhat.test}. \code{residuals(object)} returns the observed count minus the posterior-mean count per observation. \code{extract}'s \code{combineChains} and \code{fitted}/\code{predict}'s \code{ci.level} (a plain 3-column matrix here - this family has no K margin) work as for \code{\link{bart}}. \code{extract(object, type = "loglik")} is \eqn{\mathrm{dnbinom}(y_i, \mathrm{size} = r_s, \mathrm{mu} = \mu_{s,i})}{dnbinom(y_i, size = r_s, mu = mu_[s,i])}, extract-only, \code{sample = "test"} refused by name. \code{forest}/\code{contribution} on \code{extract}/\code{predict} and \code{sample} on \code{fitted} are refused by name. \code{plot(object)} traces the integer-gridded dispersion \eqn{r} (a step plot; there is no burn-in channel to bridge from, since \code{bart2} drives one \code{n.burn}/\code{n.samples} sweep for this family) alongside \code{\link{bart}}'s own gaussian observed-vs-fitted panel, applied to the counts. \code{plotTree} and \code{\link{survivalProbabilities}} are refused by name. \code{as_draws_array}/\code{as_draws_df} default to \code{vars = c("dispersion", "sigma", "k", "tau")}, matching \code{summary}'s own default.

\code{bart2(family = "hurdle.lognormal")} returns its own list, of class \code{"bartHurdle"}, holding both component fits (\code{$occupancy}, a \code{"bart"} probit fit of the occupancy indicator; \code{$positive}, a \code{"bart"} gaussian fit of \eqn{\log y}{log y} on the positive subset) plus \code{call} and \code{family} (\code{"hurdle.lognormal"}). \code{extract}/\code{fitted}/\code{predict}/\code{residuals}/\code{print} methods for \code{"bartHurdle"} compose the two components: \code{type = "ev"} (the default) is the natural-scale posterior mean or draws described under the \code{family} argument above, \code{type = "prob"} the occupancy probability, \code{type = "link"}/\code{"log"} the positive part's log-scale linear predictor, and \code{type = "ppd"} the bimodal predictive draw. \code{ci.level} on \code{fitted}/\code{predict} works as for \code{\link{bart}}. \code{residuals(object)} returns the observed response minus \code{fitted(object)} on the natural scale. \code{extract(object, type = "loglik")} is the density of \eqn{y_i} at its own natural scale (a \eqn{-\log y_i}{-log y_i} Jacobian against the components' log-scale channels, so it is comparable to any other model OF y, and differs from a log-y computation by the constant \eqn{-\sum_{y_i > 0} \log y_i}{-sum(log y_i), y_i > 0}): \eqn{\log(1 - \pi_{s,i})}{log(1 - pi_[s,i])} where \eqn{y_i = 0}, else \eqn{\log \pi_{s,i} + \mathrm{dnorm}(\log y_i, f_{s,i}, \sigma_s, \log = \mathrm{TRUE}) - \log y_i}{log pi_[s,i] + dnorm(log y_i, f_[s,i], sigma_s, log = TRUE) - log y_i}, with NO truncation (the lognormal positive part's support is already \eqn{(0, \infty)}{(0, Inf)}, unlike a count hurdle's positive part, which would have to be truncated at zero - the two rules must not be conflated); it is NOT the sum of the two components' own \code{"loglik"} channels, extract-only, \code{sample = "test"} refused by name. \eqn{\pi}{pi} here is \eqn{P(Y > 0 \mid x)}{P(Y > 0 | x)}, the complement of the zero probability other implementations report (\pkg{brms}'s \code{hu}), so the zero branch is \eqn{\log(1 - \pi)}{log(1 - pi)} rather than \eqn{\log \mathrm{hu}}{log hu}; a formula ported from one of those must be inverted. \code{forest}/\code{contribution} on \code{extract}/\code{predict} are refused by name (each component is a single forest). \code{plot(object)} draws four panels: the positive component's sigma trace; the occupancy probability \eqn{\pi(x)}{pi(x)}; the positive part on the scale it fit (\eqn{\log y}{log y} over the \eqn{y > 0} rows); and the composed natural-scale mean \eqn{E[y \mid x] = \pi e^{f + \sigma^2/2}}{E[y|x] = pi*exp(f + sigma^2/2)} over all n rows, the only panel showing the model this family exists for. \code{plotTree} and \code{\link{survivalProbabilities}} are refused by name, naming \code{object$occupancy$fit}/\code{object$positive$fit} as the route to the trees. \code{as_draws_array}/\code{as_draws_df} return the union of both components' present scalar fields, labelled \code{occupancy.<field>}/\code{positive.<field>} (a dot, not a bracket, since \pkg{posterior} parses a bracket as an index) - the same two blocks \code{summary(object)} prints under. \code{summary(object)} reports both components - \code{$occupancy}'s and \code{$positive}'s own \code{summary.bart} tables, under their own headers - rather than pooling them into one, since the two component fits share no parameters.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2010)
   BART: Bayesian additive regression trees.
        \emph{The Annals of Applied Statistics}, \bold{4}(1), 266--298.
        \doi{10.1214/09-AOAS285}.

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
\code{\link{bart}} for the BayesTree-compatible interface and the shared \code{"bart"}-class Value/Details (Decision Rules, end-node \code{k}, Generics, Saving, Reproducibility, Extracting Trees).

\code{\link{dbarts}} for the mutable sampler bart2 builds on, \code{\link{forest}} for the constructor \code{forest()} formula terms desugar to, and \code{\link{rbart_vi}}/\code{\link{xbart}} for the other formula-first entry points.
}
\examples{
## binary response with a logistic link
set.seed(0)
n <- 60L
x.bin <- matrix(runif(n * 2), n, 2)
f.bin <- 3 * x.bin[,1] - 1.5
y.bin <- rbinom(n, 1L, plogis(f.bin))
fit.logit <- bart2(y.bin ~ x.bin, family = "logistic",
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L)

## a forest() formula term: a second forest whose amplitude is modulated
## by a binary z, alongside the fit's own (prognostic) forest - the
## Bayesian causal forest shape, declared inline rather than through
## dbarts()'s forests =
set.seed(1)
n <- 60L
x1 <- runif(n); x2 <- runif(n)
z  <- rbinom(n, 1L, 0.5)
y  <- x1 + z * (1 + x2) + rnorm(n, 0, 0.2)
fit.bcf <- bart2(y ~ x1 + x2 + z:forest(x1 + x2),
                 n.samples = 20L, n.burn = 20L, n.chains = 1L,
                 n.trees = 10L, n.threads = 1L)
fit.bcf$n.forests
modulatingFit <- extract(fit.bcf, type = "forest", forest = 2L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
