\name{bart}
\alias{bart}
\alias{plot.bart}
\alias{extract}
\alias{predict.bart}
\alias{extract.bart}
\alias{fitted.bart}
\alias{residuals.bart}
\title{Bayesian Additive Regression Trees}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner.

\itemize{
  \item For numeric response \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary response \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi} denotes the standard normal cdf (probit link).
  \item \code{bart} is the BayesTree-compatible surface: besides the default gaussian/probit response it reaches \code{family = "logistic"} and \code{family = "aft"} directly (below), each packaging as an ordinary \code{"bart"} object. \code{\link{bart2}} reaches further response families - logistic, accelerated failure time, discrete-time survival hazard, K-category multinomial, ordered categorical, negative-binomial counts, and semicontinuous two-part responses - plus its own formula-term multi-forest interface; naming one of \code{bart2}'s remaining own-class families to \code{bart}'s \code{family} is refused, pointing there.
}
}
\usage{
bart(
    x.train, y.train, x.test = matrix(0.0, 0, 0),
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2.0, prior.scale = NA_real_,
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
    resid.dist = gaussian,
    subset = NULL,
    storage = c("double", "single"),
    family = c("auto", "logistic", "aft"))

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
    type = c("ev", "ppd", "bart", "loglik", "trees", "forest"),
    sample = c("train", "test"),
    combineChains = TRUE,
    forest = NULL,
    contribution = FALSE, \dots)

\method{fitted}{bart}(
    object,
    type = c("ev", "ppd", "bart"),
    sample = c("train", "test"),
    ci.level = NULL,
    \dots)

\method{residuals}{bart}(object, type = "ev", \dots)
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
     For continuous response models, an estimate of the residual standard deviation (residual standard error), \eqn{\sigma}, used to calibrate an inverse-chi-squared prior on the error variance. If not supplied, the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary. Same concept as \code{sigma} in \code{\link{dbarts}}; \code{\link{xbart}} spells it \code{sigest} too.
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
   \item{k}{
     For numeric \eqn{y}, \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}. The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}. For binary \eqn{y} fit with the default probit link, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from \eqn{\pm 3}{+/- 3}, the probit reference scale; \code{family = "logistic"} widens this to \eqn{\pm \pi \sqrt{3}}{+/- pi*sqrt(3)}, which is three standard deviations of the standard logistic latent variable (\eqn{\pi / \sqrt{3}}{pi/sqrt(3)} each), the same three-sd span probit's \eqn{\pm 3}{+/- 3} covers. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be. The value can be either a fixed number, or a \emph{hyperprior} of the form \code{chi(degreesOfFreedom = 1.5, scale = 2)}. \code{bart}'s own default is the fixed value 2; see \code{\link{bart2}}'s \code{k} item for its \code{NULL} default, a \code{chi(1.5, 2)} hyperprior on binary responses. The default of 2 for continuous responses follows Chipman, George, and McCulloch's argument (see References) that with node prior standard deviation \eqn{\sigma_\mu = 0.5 / (k \sqrt{m})}{sigma_mu = 0.5 / (k * sqrt(m))} for \eqn{m} trees, \eqn{k} prior standard deviations of \eqn{f(x)} span the whole coded response range regardless of \eqn{m} -- so \eqn{k = 2} places that range at roughly a 95\% prior interval. The range scaling itself is an outlier-sensitive convention shared with BayesTree/bartMachine: extreme \eqn{y} values compress the effective prior on everything else, for which the published workaround is to log-transform or winsorize such values before fitting, or to let \code{k} adapt via the \code{chi} hyperprior.
   }
   \item{prior.scale}{
     Names the leaf calibration in response units instead of inheriting it from the range of \code{y.train}: the prior standard deviation of the forest total \eqn{f} at \code{k = 1}, so the prior standard deviation in force is \code{prior.scale / k} and the prior mean is the response transform's shift. \code{NA} (the default) leaves the family's own calibration in place and nothing changes. Useful when the fit is one block of a larger sampler and the response it is handed varies in scale between sweeps. For the leaf-model qualifications on what the named quantity bounds, and for the \code{sd} spelling and its refusal under a sampled \code{k}, see \code{\link{dbarts}}'s Details and \code{\link{dbartsPriors}}.
   }
   \item{power}{
     Power parameter for tree prior. Default 2, Chipman, George, and McCulloch's empirical recommendation (see References) for the split-probability decay \eqn{base (1 + depth)^{-power}}{base * (1 + depth)^(-power)}, not derived from a formula.
   }
   \item{base}{
     Base parameter for tree prior. Default 0.95, from the same source as \code{power}.
   }
   \item{splitprobs}{
     Prior and transition probabilities of variables used to generate splits. Can be missing/empty/\code{NULL} for equiprobability, a numeric vector of length equal to the number variables, or a named numeric vector with only a subset of the variables specified and a \code{.default} named value. Values given for factor variables are replicated for each resulting column in the generated model matrix. The symbol \code{numvars} is rebound before execution to the number of columns in the model matrix.
   }
   \item{binaryOffset}{
     Used for binary \eqn{y}. When present, the model is \eqn{P(Y = 1 \mid x) = \Phi(f(x) + \mathrm{binaryOffset})}{P(Y = 1 | x) = \Phi(f(x) + binaryOffset)}, allowing fits with probabilities shrunk towards values other than \eqn{0.5}. \eqn{\Phi} is the link for the default probit family; a \code{family = "logistic"} fit uses \code{plogis} instead, as noted in the value section below.
   }
   \item{weights}{
      An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. A probit fit (\code{bart}, or \code{\link{bart2}} with the default binary family) does not support weights, except that weights identically 1 are treated as absent; a \code{family = "logistic"} fit treats them as observation counts and requires positive integers. For a weighted logistic fit, the \code{"ppd"} draw at an observation with weight \eqn{w} is the number of successes among \eqn{w} trials, \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} with \eqn{p} the fitted probability.
   }
   \item{ntree}{
     The number of trees in the sum-of-trees formulation.
   }
   \item{ndpost}{
     The number of posterior draws after burn in, \code{ndpost / keepevery} will actually be returned.
   }
   \item{nskip}{
     Number of MCMC iterations to be treated as burn in.
   }
   \item{printevery}{
     As the MCMC runs, a message is printed every \code{printevery} draws.
   }
   \item{keepevery}{
     Every \code{keepevery} draw is kept to be returned to the user. Useful for \dQuote{thinning} samples.
   }
   \item{keeptrainfits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of \code{x.train} are returned.
   }
   \item{usequants}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived from the \code{x.train} variables. When \code{FALSE}, splits are determined using values equally spaced across the range of a variable. See details for more information.
   }
   \item{numcut}{
     The maximum number of possible values used in decision rules (see \code{usequants}, details). If a single number, it is recycled for all variables; otherwise must be a vector of length equal to \code{ncol(x.train)}. Fewer rules may be used if a covariate lacks enough unique values.
   }
   \item{printcutoffs}{
     The number of cutoff rules to printed to screen before the MCMC is run. Given a single integer, the same value will be used for all variables. If 0, nothing is printed.
   }
   \item{verbose}{
     Logical; if \code{FALSE} suppress printing.
   }
   \item{nchain}{
     Integer specifying how many independent tree sets and fits should be calculated. Default 1 for \code{bart}, BayesTree's historical single-chain convention; \code{\link{bart2}}'s (and \code{\link{rbart_vi}}'s) default is 4.
   }
   \item{nthread, n.threads}{
     Integer specifying how many threads to use. Depending on the CPU architecture, using more than the number of chains can degrade performance for small/medium data sets. As such some calculations may be executed single threaded regardless.
   }
   \item{combinechains, combineChains}{
     Logical; if \code{TRUE}, samples will be returned in arrays of dimensions equal to \code{nchain} \eqn{\times} \code{ndpost} \eqn{\times} number of observations. Default \code{TRUE} across \code{bart}, \code{\link{bart2}}, and \code{\link{rbart_vi}}; the number of chains is tracked regardless on the returned object's \code{n.chains} component (see \sQuote{Value}).
   }
   \item{keeptrees}{
     Logical; must be \code{TRUE} in order to use \code{predict} with the result of a \code{bart} fit. Note that for models with a large number of observations or a large number of trees, keeping the trees can be very memory intensive.
   }
   \item{keepcall}{
     Logical; if \code{FALSE}, returned object will have \code{call} set to \code{call("NULL")}, otherwise the call used to instantiate BART.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. A \code{\link{set.seed}} beforehand suffices for reproducibility; supplying \code{seed} instead gives reproducible results without touching R's stream. See Reproducibility section below.
   }
   \item{proposalprobs}{
     Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step. Defaults are 0.5, 0.1, 0.4, and 0.5 respectively.
   }
   \item{keepsampler}{
     Logical that can be used to save the underlying \code{\link{dbartsSampler-class}} object even if \code{keeptrees} is false.
   }
   \item{storage}{
     Passed through to \code{\link{dbartsControl}}; a formal on both \code{bart} and \code{\link{bart2}}. A character string selecting the precision of the internal running residual; see \code{dbartsControl}'s \code{storage} for the full description. The default \code{"double"} reproduces existing draws bitwise; \code{"single"} changes the numbers a fit returns - the sampled values shift slightly - in exchange for a bandwidth-bound speedup, and is currently supported only for continuous (gaussian) responses with constant leaves.
   }
   \item{family}{
     \code{bart}'s own \code{family} is a narrower, appended formal, \code{c("auto", "logistic", "aft")}: the first-token default \code{"auto"} reproduces today's behavior exactly (a binary \code{y.train} still resolves to probit), and \code{"logistic"}/\code{"aft"} are newly reachable directly, both packaging as an ordinary \code{"bart"} object (\code{"aft"} needs a \code{\link[survival]{Surv}} or two-column \code{(time, status)} \code{y.train}, as below). The other ten \code{\link{bart2}} tokens are refused BY NAME rather than falling through to \code{\link{match.arg}}'s generic \dQuote{should be one of} message, which names neither the token nor \code{bart2}: \code{"multinomial"}, \code{"ordinal"}, \code{"nbinom"}, and \code{"hurdle.lognormal"} package as their own S3 class or compose two samplers, so \code{bart2(x.train, y.train, family = ...)} is where they are fit; \code{"gaussian"}, \code{"probit"}, and \code{"hazard.probit"} are already reachable through \code{family = "auto"} plus the response, so they add no capability as separate \code{bart} tokens; and \code{"hazard"}/\code{"hazard.logistic"} need \code{breaks}/\code{max.rows}, which \code{bart} does not have.
   }
   \item{subset}{
     A vector of logicals or indices used to subset the data; an appended formal defaulting to \code{NULL} (today's behavior, unchanged), forwarded to \code{\link{dbarts}} the same way. See \code{\link{bart2}}'s \code{subset} item for how a \code{forest()} term composes with it on that interface.
   }
   \item{object}{
     An object of class \code{bart}, returned from either the function \code{bart} or \code{bart2}.
   }
   \item{newdata}{
     Test data for prediction. Obeys all the same rules as \code{x.train} but cannot be missing.
   }
   \item{offset}{
     For \code{predict}: an optional offset for \code{newdata}, applied the same way \code{\link{dbarts}}'s \code{offset.test} is at fit time. \code{NULL} (the default, when missing) applies none.
   }
   \item{sampleronly}{
     Builds the sampler from its arguments and returns it without running it. Useful to use the \code{bart} interface in more complicated models.
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
     The quantity to be returned by generic functions. Options are \code{"ev"} - samples from the posterior of the individual level expected value, \code{"bart"} - the sum of trees component; same as \code{"ev"} for linear models but on the latent scale (probit or logistic, matching the fit's family) for binary ones, \code{"ppd"} - samples from the posterior predictive distribution, \code{"loglik"} - for \code{extract} only, the log-likelihood of each training observation at each posterior draw, \code{"trees"} - a data frame with tree information for when model was fit with \code{keepTrees} equal to \code{TRUE}, and \code{"forest"} - for \code{extract} only, the per-forest channels of an amplitude-coupled multi-forest fit (see \code{forest}, \code{contribution}, and the \code{forestFits}/\code{glue}/\code{bases} components under \sQuote{Value}); an error naming the reason on any other fit. For \code{"ppd"}, a weighted logistic fit draws the number of successes among the observation-count weight, \eqn{\mathrm{Binomial}(w_i, p_i)}{Binomial(w_i, p_i)} (see \code{weights}), and an aft (survival) fit draws on the log-time scale (the fit models \eqn{\log T}{log T}). For \code{"loglik"}, gaussian fits evaluate \eqn{y_i \mid x_i \sim N(\hat{f}(x_i), \sigma^2 / w_i)}{y_i | x_i ~ N(f(x_i), sigma^2 / w_i)} in logs at each draw of \eqn{f} and \eqn{\sigma}{sigma}, a \code{resid.dist = student()} fit the corresponding marginal \eqn{t_\nu}{t_nu} density at that draw's \eqn{\nu}{nu} (the fit's \code{$resid.df}) rather than the normal one, binary fits evaluate the Bernoulli log-likelihood of the fitted probability, multiplied for a weighted logistic fit by the observation-count weight \eqn{w_i}, and an aft fit contributes the log density for an event and the log survival tail \eqn{\log P(T > C)}{log P(T > C)} for a right-censored observation, both on the log-time scale. When chains are combined the result is a samples-by-observations matrix directly consumable by WAIC/PSIS-LOO implementations such as those in the \pkg{loo} package; the chains-first convention is kept, so a per-chain array (\code{combineChains = FALSE}, dimension chains-by-samples-by-observations) is reordered to the draws-by-chains-by-observations that \code{loo::relative_eff} expects with \code{aperm(x, c(2, 1, 3))}. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"ev"} and \code{"link"} can be used as a synonym for \code{"bart"}. For information on extracting trees, see the subsection below.
   }
   \item{sample}{
     Either \code{"train"} or \code{"test"}. \code{"test"} is refused, by name, when \code{type} is \code{"forest"}: an amplitude-coupled fit has no test fits (see \code{\link{dbarts}}'s \code{forests =}).
   }
   \item{forest}{
     For \code{extract(type = "forest")} only: which forest(s) to return, by 1-based index or by margin name (\code{"forest1"}, \code{"forest2"}, \dots); \code{NULL} (the default) returns every forest. The returned array always keeps the trailing forest margin, subset to the requested forests, even when only one is selected.
   }
   \item{contribution}{
     For \code{extract(type = "forest")} only: \code{FALSE} (the default) returns the raw per-forest total \code{forestFits} already carries (see \sQuote{Value}). \code{TRUE} instead returns each selected forest's per-observation CONTRIBUTION to the fit, \eqn{(\mathrm{basis}_k \, \mathrm{glue}_k) \times \mathrm{forestFits}_k}{(basis_k \%*\% glue_k) * forestFits_k}, computed on demand from the stored \code{forestFits}, \code{glue}, and \code{bases} rather than stored itself.
   }
   \item{ci.level}{
     For \code{fitted} and \code{predict}, an optional single number in \eqn{(0, 1)}. When \code{NULL} (the default) the posterior mean is returned as before. When given, the result is instead a matrix with one row per observation and columns \code{est} (the posterior mean), \code{ci.lower}, and \code{ci.upper} giving a symmetric credible band at that level, pooled over all samples and chains. The kind of interval follows \code{type}: \code{"ev"} yields a credible interval for the conditional mean (a probability for binary responses), \code{"ppd"} a prediction interval that additionally carries the residual noise (and so is wider), and \code{"bart"} a credible interval on the latent scale.
   }
   \item{\dots}{
     Additional arguments passed on to \code{plot} (via \code{plot.bart}), or to the sampler's \code{getTrees} method via \code{extract} when \code{type} is \code{"trees"} (\code{chainNums}, \code{sampleNums}, \code{treeNums}, \code{newdata}; see \sQuote{Extracting Trees} below). Not used in \code{predict}.
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
  \code{bart} and \code{\link{bart2}} return lists assigned the class \code{bart} for the shared gaussian/probit/logistic/aft/hazard* families (\code{bart2} reaches further response families under their own list classes and generics; see \code{\link{bart2}}'s Value section). For applicable quantities, \code{ndpost / keepevery} samples are returned. In the numeric \eqn{y} case, the list has components:

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
  \item{\code{forestFits}, \code{glue}, \code{bases}, \code{n.forests}}{
    Present only for an amplitude-coupled multi-forest fit: built through \code{\link{dbarts}}'s \code{forests =}, or through \code{\link{bart2}}'s \code{forest()} formula term (see its \sQuote{Formula Terms} section); not reachable from \code{bart} directly, which has no multi-forest front door. \code{n.forests} is the forest count K. \code{forestFits} is a (\code{n.chains} \eqn{\times}, when uncombined) \code{n.samples} \eqn{\times} number of training observations \eqn{\times} K array: forest k's own RESPONSE-scale raw total at each draw, i.e. \eqn{\mathrm{response.scale} \times f_k(x)}{response.scale * f_k(x)}, with NO amplitude (\code{glue}) folded in - the same quantity the sampler's \code{$getForestFits} reports, up to that one scalar. The trailing margin is named \code{forest1}, ..., \code{forestK}; a declaration's own names (\code{names(forests)} on the \code{dbarts()} route, or a term's own text on the \code{bart2} route) ride separately as a \code{"forest.labels"} attribute on the fit object when given. \code{glue} is the ragged (\code{n.chains} \eqn{\times}, when uncombined) \code{n.samples} \eqn{\times} sum(q_k) matrix of amplitudes multiplying each forest's basis columns, forest-major (the layout \code{$getForestAmplitudes} documents); its \code{"forest"} attribute names each column's forest. \code{bases} is a length-K list of each forest's expanded basis matrix (no draw axis; \code{NULL} for a forest with no declared basis, whose amplitude then scales an implicit all-ones column). Together they reconstruct the training fit exactly: \eqn{\hat y = \mathrm{response.shift} + \sum_k (\mathrm{bases}_k \, \mathrm{glue}_k) \times \mathrm{forestFits}_k}{yhat = response.shift + sum_k (bases_k \%*\% glue_k) * forestFits_k}, with \code{response.scale}/\code{response.shift} read from the sampler's \code{$getCalibration}. See \code{\link{extract}}'s \code{type = "forest"}.
  }
  \item{\code{sigest}}{
    The rough error standard deviation (\eqn{\sigma}) used in the prior.
  }
  \item{\code{resid.dist}, \code{resid.df}}{
    The residual error law the fit was made under, \code{"gaussian"} or \code{"student"} (see \code{\link{dbarts}}'s \code{resid.dist}). \code{resid.df} is present only for a \code{student()} fit: the degrees of freedom \eqn{\nu} each draw was conditioned on, in \code{sigma}'s layout - a fixed \eqn{\nu} repeats its value, an estimated one gives that draw's grid value. \code{extract(type = "loglik")} reads both.
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
  \item{\code{binaryOffset}}{
    Present only for a binary fit: the offset value used.
  }
  \item{\code{family}}{
    The resolved response family (\code{"gaussian"}, \code{"probit"}, \code{"logistic"}, or \code{"aft"}). \code{predict}, \code{extract}, \code{fitted}, and \code{plot} use it to transform latent draws to probabilities. For an \code{"aft"} fit, predictions and fitted values are on the LOG-TIME scale: these functions return the linear predictor \eqn{E[\log T \mid x]}, exactly as for a gaussian fit of \eqn{\log T}, and never the time scale that \code{survreg} or \code{flexsurv} users might expect. The \code{type} aliases \code{"response"} and \code{"link"} are inert for \code{"aft"} - both return the same log-scale linear predictor; neither exponentiates. The posterior-median survival time is \code{exp} of the linear predictor, and \code{\link{survivalProbabilities}} gives survival-probability draws. A discrete-time hazard fit (\code{family = "hazard"}/\code{"hazard.logistic"}, reachable only through \code{\link{bart2}}) records the binary link here (\code{"probit"}/\code{"logistic"}) - the fit IS that binary fit on the expanded rows - and carries a separate \code{$periods} element (the ordered period grid); \code{predict}/\code{extract} return the per-(subject, period) hazard through the link, and \code{\link{survivalProbabilities}} dispatches its survival-curve branch on the \code{$periods} marker.
  }

  In the binary \eqn{y} case, the returned list has the components \code{yhat.train}, \code{yhat.test}, and \code{varcount} as above, but not \code{yhat.train.mean}/\code{yhat.test.mean} - use \code{fitted} to get the posterior mean of \eqn{P(Y = 1 \mid x)}{P(Y = 1 | x)} instead. In addition the list has a \code{binaryOffset} component giving the value used.

  Note that in the binary \eqn{y}, case \code{yhat.train} and \code{yhat.test} are \eqn{f(x) + \mathrm{binaryOffset}}{f(x) + binaryOffset}. For draws of the probability \eqn{P(Y = 1 | x)}, apply the family's inverse link to these values: the normal cdf (\code{pnorm}) for \code{"probit"} fits, \code{plogis} for \code{"logistic"} ones.

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
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com},
Robert McCulloch: \email{robert.mcculloch1@gmail.com},
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{bart2}} for the formula-first interface, further response families, and formula-term multi-forest models.

\code{\link{pdbart}}

\code{\link{dbarts-embedding}} and \code{vignette("dbarts-as-a-component", package = "dbarts")} for
using BART as one block of a sampler you write, rather than as the whole fit: these functions own
their own MCMC loop, while \code{\link{dbarts}} returns a
\code{\linkS4class{dbartsSampler}} that runs one sweep at a time.
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

## fit with missing predictor values, using missing = "incorporate"
## (the default): every split rule learns a direction for NAs
set.seed(1)
n <- 60L
x.na <- matrix(runif(n * 2), n, 2)
x.na[1:5, 1] <- NA
y.na <- x.na[,2] + rnorm(n, 0, 0.3)
y.na[is.na(x.na[,1])] <- y.na[is.na(x.na[,1])] + 2
fit.na <- dbarts(y.na ~ x.na, control = dbartsControl(
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L))
samples.na <- fit.na$run()

## recipe: between-draws row subsetting on a latent family (e.g. principal
## stratification or mediation, where an outer Gibbs step redraws which rows
## belong to a stratum every iteration). $setActiveRows changes which rows
## are "in" the model from one run() to the next, without rebuilding the
## sampler and without losing any row's fitted value; see the "active" entry
## of \link{dbartsSampler-class} for the full semantics. Unlike zero case
## weights (gaussian only), this reaches probit, ordinal, logistic, nbinom
## and aft samplers - and an emptied stratum is a legal all-zeros mask
## rather than a case the caller must special-case and skip.
set.seed(3)
n <- 120L
x.strat <- matrix(runif(n * 2), n, 2)
y.strat <- rbinom(n, 1L, plogis(2 * x.strat[,1] - 1))
sampler.strat <- dbarts(x.strat, y.strat, family = "probit",
                        control = dbartsControl(n.chains = 1L, n.threads = 1L,
                          n.trees = 20L, updateState = FALSE))

stratum <- rbinom(n, 1L, 0.6)  # this iteration's stratum membership
sampler.strat$setActiveRows(stratum)
draws.strat <- sampler.strat$run(20L, 20L)

## the next outer iteration redraws the stratum; even a stratum that empties
## out is accepted and runs, rather than needing to be skipped
stratum <- rep(0, n)
sampler.strat$setActiveRows(stratum)
draws.empty <- sampler.strat$run(0L, 5L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
