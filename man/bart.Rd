\name{bart}
\alias{bart}
\alias{bart2}
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
    keepsampler = keeptrees)

bart2(
    formula, data, test, subset, weights, offset, offset.test = offset,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL,
    power = 2.0, base = 0.95, split.probs = 1 / num.vars,
    dart = FALSE,
    n.trees = 75L,
    n.samples = 500L, n.burn = 500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = FALSE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 1L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE, keepTrees = FALSE, 
    keepCall = TRUE, samplerOnly = FALSE,
    seed = NA_integer_,
    proposal.probs = NULL,
    keepSampler = keepTrees,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    missing = c("incorporate", "error"),
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
    type = c("ev", "ppd", "bart", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE, \dots)

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
     Dependent variable for training (in sample) data. If \code{y.train} is numeric a continous response model is fit (normal errors). If \code{y.train} is a binary factor or has only values 0 and 1, then a binary response model with a probit link is fit.
   }
   \item{x.test}{
     Explanatory variables for test (out of sample) data. Should have same column structure as \code{x.train}. \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test}.
   }
   \item{sigest}{
     For continuous response models, an estimate of the error variance, \eqn{\sigma^2}, used to calibrate an inverse-chi-squared prior used on that parameter. If not supplied, the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary.
   }
   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary. Default 3, Chipman, George, and McCulloch's calibration (see References); an aggressive choice rather than a derived one.
   }
   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate (\code{sigest}) is placed at. The closer the quantile is to 1, the more aggresive the fit will be as you are putting more prior weight on error standard deviations (\eqn{\sigma}) less than the rough estimate. Not applicable when \eqn{y} is binary. Default 0.90, from the same source as \code{sigdf}.
   }
   \item{k}{
     For numeric \eqn{y}, \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}. The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}. For binary \eqn{y} fit with the default probit link, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from \eqn{\pm 3}{+/- 3}, the probit reference scale; \code{family = "logistic"} (\code{bart2} only) widens this to \eqn{\pm \pi \sqrt{3}}{+/- pi*sqrt(3)}, the standard deviation of the standard logistic latent variable. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be. The value can be either a fixed number, or the a \emph{hyperprior} of the form \code{chi(degreesOfFreedom = 1.25, scale = Inf)}. For \code{bart2}, the default of \code{NULL} uses the value 2 for continuous reponses and a \code{chi} hyperprior for binary ones. The default \code{chi} hyperprior is improper, and slightly penalizes small values of \code{k}. The default of 2 for continuous responses follows Chipman, George, and McCulloch's argument (see References) that with node prior standard deviation \eqn{\sigma_\mu = 0.5 / (k \sqrt{m})}{sigma_mu = 0.5 / (k * sqrt(m))} for \eqn{m} trees, \eqn{k} prior standard deviations of \eqn{f(x)} span the whole coded response range regardless of \eqn{m} -- so \eqn{k = 2} places that range at roughly a 95\% prior interval. The range scaling itself is an outlier-sensitive convention shared with BayesTree/bartMachine: extreme \eqn{y} values compress the effective prior on everything else, for which the published workaround is to log-transform or winsorize such values before fitting, or to let \code{k} adapt via the \code{chi} hyperprior.
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
      An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. A probit fit (\code{bart}, or \code{bart2} with the default binary family) does not support weights, except that weights identically 1 are treated as absent; \code{bart2} with \code{family = "logistic"} treats them as observation counts and requires positive integers.
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
     Logical; if \code{FALSE} supress printing.
   }
   \item{nchain, n.chains}{
     Integer specifying how many independent tree sets and fits should be calculated.
   }
   \item{nthread, n.threads}{
     Integer specifying how many threads to use. Depending on the CPU architecture, using more than the number of chains can degrade performance for small/medium data sets. As such some calculations may be executed single threaded regardless.
   }
   \item{combinechains, combineChains}{
     Logical; if \code{TRUE}, samples will be returned in arrays of dimensions equal to \code{nchain} \eqn{\times} \code{ndpost} \eqn{\times} number of observations.
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
   \item{keepsampler, keepSampler}{
     Logical that can be used to save the underlying \code{\link{dbartsSampler-class}} object even if \code{keepTrees} is false.
   }
   \item{factors, family, missing}{
     As in \code{\link{dbarts}}; \code{bart2} only. \code{bart} keeps the historical behavior: factors expand into indicator columns, binary responses are probit, and missing data are rejected.
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
     A vector of logicals or indicies used to subset of the data. Can be missing.
   }
   \item{offset}{
     The same as \code{binaryOffset}. Can be missing.
   }
   \item{offset.test}{
     A vector of offsets to be used with test data, in case it is different than the training offset. If \code{offest} is missing, defaults to \code{NULL}.
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
     The quantity to be returned by generic functions. Options are \code{"ev"} - samples from the posterior of the individual level expected value, \code{"bart"} - the sum of trees component; same as \code{"ev"} for linear models but on the latent scale (probit or logistic, matching the fit's family) for binary ones, \code{"ppd"} - samples from the posterior predictive distribution, and \code{"trees"} - a data frame with tree information for when model was fit with \code{keepTrees} equal to \code{TRUE}. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"ev"} and \code{"link"} can be used as a synonym for \code{"bart"}. For information on extracting trees, see the subsection below.
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
    The amount of shrinkage of the node parameters is controlled by \code{k}. \code{k} can be given as either a fixed, positive number, or as any value that can be used to build a supported hyperprior. At present, only \eqn{\chi_\nu s} priors are supported, where \eqn{\nu} is a degrees of freedom and \eqn{s} is a scale. Both values must be positive, however the scale can be infinite which yields an improper prior, which is interpretted as just the polynomial part of the distribution. If \eqn{nu} is 1 and \eqn{s} is \eqn{\infty}, the prior is \dQuote{flat}.

    For BART on binary outcomes, the degree of overfitting can be highly sensitive to \code{k} so it is encouraged to consider a number of values. The default hyperprior for binary BART, \code{chi(1.25, Inf)}, has been shown to work well in a large number of datasets, however crossvalidation may be helpful. Running for a short time with a flat prior may be helpful to see the range of values of \code{k} that are consistent with the data.
  }
  \subsection{Generics}{
    \code{bart} and \code{\link{rbart_vi}} support \code{\link{fitted}} to return the posterior mean of a predicted quantity, as well as \code{\link{predict}} to return a set of posterior samples for a different sample. In addition, the \code{extract} generic can be used to obtain the posterior samples for the training data or test data supplied during the initial fit.

    Using \code{predict} with a \code{bart} object requires that it be fitted with the option \code{keeptrees}/\code{keepTrees} as \code{TRUE}. Keeping the trees for a fit can require a sizeable amount of memory and is off by default.

    All generics return values on the scale of expected value of the response by default. This means that \code{predict}, \code{extract}, and \code{fitted} for binary outcomes return probabilities unless specifically the sum-of-trees component is requested (\code{type = "bart"}). This is in contrast to \code{yhat.train}/\code{yhat.test} that are returned with the fitted model.
  }
  
  \subsection{Saving}{
    \code{\link{save}}ing and \code{\link{load}}ing fitted BART objects for use with \code{predict} requires that R's serialization mechanism be able to access the underlying trees, in addition to being fit with \code{keeptrees}/\code{keepTrees} as \code{TRUE}. For memory purposes, the trees are not stored as R objects unless specifically requested. To do this, one must \dQuote{touch} the sampler's state object before saving, e.g. for a fitted object \code{bartFit}, execute \code{invisible(bartFit$fit$state)}.
  }
  
  \subsection{Reproducibility}{
    Every chain runs its own pseudo random number generator, so worker threads never touch R's generator and results do not depend on the thread count. Without a \code{seed}, chain generators are seeded from R's stream when the sampler is created, so calling \code{\link{set.seed}} beforehand makes results reproducible. Sampling itself never advances R's stream: after a fit, \code{\link{.Random.seed}} has moved only by the draws taken at creation.

    Setting \code{seed} makes results reproducible without involving R's stream at all: the seed drives a dedicated generator that hands each chain its own seed. A single-chain run with a given seed reproduces the first chain of a multi-chain run with the same seed.
  }

  \subsection{Extracting Trees}{
    When a model is fit with \code{keeptrees} (\code{bart}) or \code{keepTrees} (\code{bart2}) equal to \code{TRUE}, the generic \code{extract} can be used to retrieve a data frame containing the tree fit information. In this case, \code{extract} will accept the additional, optional arguments: \code{chainNums}, \code{sampleNums}, and \code{treeNums}. Each should be an integer vector detailing the desired trees to be returned. A further optional argument \code{newdata} routes a new set of predictors (in the same form accepted by \code{predict}) through the frozen trees so that the \code{n} column counts those observations instead of the training data.

    The result of \code{extract} will be a data frame with columns:
    \itemize{
      \item \code{sample}, \code{chain}, \code{tree} - index variables
      \item \code{n} - number of observations in node, drawn from the training data or from \code{newdata} when supplied
      \item \code{var} - either the index of the variable used for splitting or -1 if the node is a leaf
      \item \code{value} - either the value such that observations less than or equal to it are sent down the left path of the tree or the predicted value for a leaf node
    }
    The order of nodes in the result corresponds to a depth-first traversal, going down the left-side first. The names of variables used in splitting can be recovered by examining the column names of the \code{fit$data@x} element of a fitted \code{bart} or \code{bart2} model. See the package vignette \dQuote{Working with dbarts Saved Trees}.
  }
}
\value{
  \code{bart} and \code{bart2} return lists assigned the class \code{bart}. For applicable quantities, \code{ndpost / keepevery} samples are returned. In the numeric \eqn{y} case, the list has components:

  \item{\code{yhat.train}}{
    A array/matrix of posterior samples. The \eqn{(i, j, k)} value is the \eqn{j}th draw of the posterior of \eqn{f} evaluated at the \eqn{k}th row of \code{x.train} (i.e. \eqn{f^*(x_k)}{f(x_k)}) corresponding to chain \eqn{i}. When \code{nchain} is one or \code{combinechains} is \code{TRUE}, the result is a collapsed down to a matrix.
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
    Matrix of posterior samples of \code{sigma}, the residual/error standard deviation. Dimensions are equal to the number of chains times the numbers of samples unless \code{nchain} is one or \code{combinechains} is \code{TRUE}.
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
  \item{\code{family}}{
    The resolved response family (\code{"gaussian"}, \code{"probit"}, or \code{"logistic"}). \code{predict}, \code{extract}, \code{fitted}, and \code{plot} use it to transform latent draws to probabilities.
  }

  In the binary \eqn{y} case, the returned list has the components \code{yhat.train}, \code{yhat.test}, and \code{varcount} as above.  In addition the list has a \code{binaryOffset} component giving the value used.

  Note that in the binary \eqn{y}, case \code{yhat.train} and \code{yhat.test} are \eqn{f(x) + \mathrm{binaryOffset}}{f(x) + binaryOffset}. For draws of the probability \eqn{P(Y = 1 | x)}, apply the family's inverse link to these values: the normal cdf (\code{pnorm}) for \code{"probit"} fits, \code{plogis} for \code{"logistic"} ones.

  The \code{plot} method sets \code{mfrow} to \code{c(1, 2)} and makes two plots. The first plot is the sequence of kept draws of \eqn{\sigma} including the burn-in draws. Initially these draws will decline as BART finds a good fit and then level off when the MCMC has burnt in. The second plot has \eqn{y} on the horizontal axis and posterior intervals for the corresponding \eqn{f(x)} on the vertical axis.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2009)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning. 
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

both of the above at:
\url{https://www.rob-mcculloch.org}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.

Linero, A.R. (2018)
   Bayesian regression trees for high-dimensional prediction and variable selection.
        \emph{Journal of the American Statistical Association}, \bold{113}(522), 626--636.
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
