\docType{class}
\name{dbartsSampler-class}
\title{Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbartsSampler-class}
\alias{\S4method{run}{dbartSampler}}
\description{
  A reference class object that contains a Bayesian Additive Regression Trees sampler in such a
  way that it can be modified, stopped, and started all while maintaining its own state.
}
\usage{
\S4method{run}{dbartSampler}(numBurnIn, numSamples, updateState = NA)
\S4method{setControl}{dbartSampler}(control)
\S4method{setResponse}{dbartSampler}(y, updateState = NA)
\S4method{setOffset}{dbartSampler}(offset, updateState = NA)
\S4method{setPredictor}{dbartSampler}(x, column, updateState = NA)
}
\arguments{
  \item{numBurnIn}{A non-negative integer determining how many iterations the sampler should
    skip before storing results. If missing or \code{NA}, the default is filled in from the sampler's
    \code{\link[dbartsControl]{control}} object.}
  \item{numSamples}{A positive integer determining how many posterior samples should be
  	returned. If missing or \code{NA}, the default is also filled in from the control object.}
  \item{updateState}{A logical determining if the local cache of the sampler's state
  	should be updated after the completion of the run. If \code{NA}, the default is also
  	filled in from the control object.}
  \item{y}{A numeric vector of length equal to that with which the sampler was created.}
  \item{offset}{A numeric vector of length equal to that with which the sampler was created, or \code{NULL}.}
  \item{control}{An object inheriting from \code{\link{dbartsControl}}.}
}
\details{
  Huzzah!
}
