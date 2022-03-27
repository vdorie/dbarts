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
    control = dbartsControl(), sigma = NA_real_)}
\arguments{
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bart}} matrix \code{x.train}.
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
    An optional vector of weights to be used in the fitting process. When present, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function.
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
    An expression of the form \code{cgm} or \code{cgm(power, base, split.probs)} setting the tree prior used in fitting.
  }
  \item{node.prior}{
    An expression of the form \code{normal} or \code{normal(k)} that sets the prior used on the averages within nodes.
  }
  \item{resid.prior}{
    An expression of the form \code{chisq} or \code{chisq(df, quant)} that sets the prior used on the residual/error variance.
  }
  \item{proposal.probs}{
    Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step.
  }
  \item{control}{
    An object inheriting from \code{dbartsControl}, created by the \code{\link{dbartsControl}} function.} \item{sigma}{A positive numeric estimate of the residual standard deviation. If \code{NA}, a linear model is used with all of the predictors to obtain one.
  }
}
\details{
  \dQuote{Discrete sampler} refers to that \code{dbarts} is implemented using \link{ReferenceClasses}, so that there exists a mutable object constructed in C++ that is largely obscured from R. The \code{dbarts} function is the primary way of creating a \code{\link[=dbartsSampler-class]{dbartsSampler}}, for which a variety of methods exist.
}
\value{
  A reference object of \code{\link[=dbartsSampler-class]{dbartsSampler}}.
}
