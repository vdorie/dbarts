\name{dbarts}
\title{Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbarts}
\description{
  Creates a sampler object for a given problem which fits a Bayesian Additive Regreesion Trees model. Internally stores state in such a way as to be mutable.
}
\usage{
dbarts(
    formula, data, test, subset, weights, offset, offset.test = offset,
    verbose = FALSE, n.samples = 800L,
    tree.prior = cgm, node.prior = normal, resid.prior = chisq,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    control = dbarts::dbartsControl(), sigma = NA_real_, seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    missing = c("incorporate", "error"))}
\arguments{
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bart}} matrix \code{x.train}, including a sparse \code{Matrix::dgCMatrix}: its columns enter as ordinal predictors, sufficiently sparse columns are stored in a compact rank-bitmap layout instead of being expanded, and the predictor-mutation surface (\code{setPredictor} and relatives, \code{setData}) is fixed at creation. Sparse inputs are not supported by \code{\link{rbart_vi}} or the \code{node.prior = linear()} and \code{gp()} leaf models; test matrices remain dense. A data frame may mix ordinary columns with \code{Matrix::sparseVector} or \code{dgCMatrix} columns (assign them into the frame; they do not survive \code{data.frame(...)} or \code{I()}): sparse columns behave as in an all-sparse design, while dense columns keep categorical splits and linear-leaf designation.
  }
  \item{data}{
    An optional data frame, list, or environment containing predictors to be used with the model. For backwards compatibility, can also be the \code{\link{bart}} vector \code{y.train}.
  }
  \item{test}{
    An optional matrix or data frame with the same number of predictors as \code{data}, or \code{formula} in backwards compatibility mode. If column names are present, a matching algorithm is used.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
    An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. Binary responses differ: a \code{"probit"} model does not support weights (a weighted probit has no tractable latent-variable form), except that weights identically 1 are treated as absent; a \code{"logistic"} model treats them as observation counts and so requires positive integers (its Polya-Gamma latent for a count \eqn{w} is a sum of \eqn{w} unit draws).
  }
  \item{offset}{
    An optional vector specifying an offset from 0 for the relationship between the underyling function, \eqn{f(x)}, and the response \eqn{y}. Only is useful for binary responses, in which case the model fit is to assume \eqn{P(Y = 1 \mid X = x) = \Phi(f(x) + \mathrm{offset})}{P(Y = 1 | X = x) = \Phi(f(x) + offset)}, where \eqn{\Phi} is the standard normal cumulative distribution function.
  }
  \item{offset.test}{
    The equivalent of \code{offset} for test observations. Will attempt to use \code{offset} when applicable.
  }
  \item{verbose}{
    A logical determining if additional output is printed to the console. See \code{\link{dbartsControl}}.
  }
  \item{n.samples}{
    A positive integer setting the default number of posterior samples to be returned for each run of the sampler. Can be overriden at run-time. See \code{\link{dbartsControl}}.
  }
  \item{tree.prior}{
    An expression of the form \code{cgm} or \code{cgm(power, base, split.probs)} setting the tree prior used in fitting, or a prior object built with \code{\link{dbartsPriors}}.
  }
  \item{node.prior}{
    An expression of the form \code{normal} or \code{normal(k)} that sets the prior used on the averages within nodes, or a prior object built with \code{\link{dbartsPriors}}. \code{linear(columns, k)} instead fits each leaf with an intercept plus a linear term in the designated continuous predictor columns (character names or numeric indices into the model matrix), standardized internally, all coefficients sharing the \code{normal(k)} prior; factor columns cannot be designated. \code{gp(columns, k, lengthscale, max.leaf.size)} fits each leaf with a smooth Gaussian-process function of the designated columns under a squared-exponential kernel; leaves larger than \code{max.leaf.size} fall back to constant fits. \code{\link{xbart}} accepts the same specifications through its own \code{node.prior} argument.
  }
  \item{resid.prior}{
    An expression of the form \code{chisq} or \code{chisq(df, quant)} that sets the prior used on the residual/error variance, or a prior object built with \code{\link{dbartsPriors}}.
  }
  \item{proposal.probs}{
    Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step.
  }
  \item{control}{
    An object inheriting from \code{dbartsControl}, created by the \code{\link{dbartsControl}} function.} \item{sigma}{A positive numeric estimate of the residual standard deviation. If \code{NA}, a linear model is used with all of the predictors to obtain one.
  }
  \item{seed}{
    Optional integer seed for the random number generator, a convenience mirror of \code{dbartsControl(rngSeed = )}. When not \code{NA} it overrides the seed in \code{control}; the fitting-function wrappers (\code{\link{bart2}}, \code{\link{xbart}}) accept the same argument.
  }
  \item{factors}{
    How factor columns in a data frame enter the model. The default \code{"categorical"} keeps each unordered factor as a single predictor whose splits send a subset of its levels (at most 65535) down each branch, and codes ordered factors as ordinal; level tables are retained so that test data are coded identically. \code{"indicators"} expands each factor into binary indicator columns, as previous versions always did and as \code{\link{bart}} still does.
  }
  \item{family}{
    The response model. \code{"auto"} fits gaussian models to continuous responses and probit models to those coded 0/1, as always. \code{"gaussian"} forces a continuous fit even for a 0/1 response; \code{"probit"} and \code{"logistic"} require a 0/1 response and fit latent-variable models, with fits and predictions on the latent scale. \code{"logistic"} uses Polya-Gamma augmentation. See \code{weights} for how the families differ in their support for weights.
  }
  \item{missing}{
    How missing values in the predictors enter the model. The default \code{"incorporate"} keeps them: every split rule learns a direction for missing values on its variable, so an observation whose split value is \code{NA} follows that rule's chosen branch (\dQuote{Missingness Incorporated in Attributes}, Twala et al. 2008), in training and test data alike. \code{"error"} rejects predictors containing \code{NA}. The response, \code{weights}, and \code{offset} must always be complete; note that previous versions silently dropped incomplete rows for formula inputs.
  }
}
\details{
  \dQuote{Discrete sampler} refers to that \code{dbarts} is implemented using \link{ReferenceClasses}, so that there exists a mutable object constructed in C++ that is largely obscured from R. The \code{dbarts} function is the primary way of creating a \code{\link[=dbartsSampler-class]{dbartsSampler}}, for which a variety of methods exist.

  \subsection{Test offset synchronization}{
    When \code{offset.test} is left at its default of tracking \code{offset}, the sampler links the two: the sampler's \code{setOffset} method re-derives \code{offset.test} from each new \code{offset} it is given. Calling the sampler's \code{setTestOffset} or \code{setTestPredictorAndOffset} methods breaks this link -- \code{offset.test} is set independently from then on, and the link never re-forms, even if the two are later set to equal values.
  }
}
\value{
  A reference object of \code{\link[=dbartsSampler-class]{dbartsSampler}}.
}
\references{
Twala, B.E.T.H., Jones, M.C., and Hand, D.J. (2008)
   Good methods for coping with missing data in decision trees.
        \emph{Pattern Recognition Letters}, \bold{29}(7), 950--956.
}
