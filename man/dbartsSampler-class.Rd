\docType{class}
\name{dbartsSampler-class}
\title{Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbartsSampler}
\alias{dbartsSampler-class}
\alias{\S4method{run}{dbartSampler}}
\alias{\S4method{setControl}{dbartSampler}}
\alias{\S4method{setResponse}{dbartSampler}}
\alias{\S4method{setOffset}{dbartSampler}}
\alias{\S4method{setTestOffset}{dbartSampler}}
\alias{\S4method{setPredictor}{dbartSampler}}
\alias{\S4method{setTestPredictor}{dbartSampler}}
\alias{\S4method{setTestPredictors}{dbartSampler}}
\description{
  A reference class object that contains a Bayesian Additive Regression Trees sampler in such a
  way that it can be modified, stopped, and started all while maintaining its own state.
}
\usage{
\S4method{run}{dbartSampler}(numBurnIn, numSamples, updateState = NA)
\S4method{setControl}{dbartSampler}(control)
\S4method{setResponse}{dbartSampler}(y, updateState = NA)
\S4method{setOffset}{dbartSampler}(offset, updateState = NA)
\S4method{setTestOffset}{dbartSampler}(offset.test, updateState = NA)
\S4method{setPredictor}{dbartSampler}(x, column, updateState = NA)
\S4method{setTestPredictor}{dbartSampler}(x, column, updateState = NA)
\S4method{setTestPredictors}{dbartSampler}(x, offset.test, updateState = NA)
}
\arguments{
  \item{numBurnIn}{A non-negative integer determining how many iterations the sampler should
    skip before storing results. If missing or \code{NA}, the default is filled in from the sampler's
    \code{\link[=dbartsControl]{control}} object.}
  \item{numSamples}{A positive integer determining how many posterior samples should be
  	returned. If missing or \code{NA}, the default is also filled in from the control object.}
  \item{updateState}{A logical determining if the local cache of the sampler's state
  	should be updated after the completion of the run. If \code{NA}, the default is also
  	filled in from the control object.}
  \item{y}{A numeric response vector of length equal to that with which the sampler was created.}
  \item{x}{A numeric predictor vector of length equal to that with which the sampler was created. Can be
  	an entirely matrix of new number of rows for \code{setTestPredictor}.}
  \item{offset}{A numeric vector of length equal to that with which the sampler was created, or \code{NULL}.
  	If \code{offset.test} was set from \code{offset}, will attempt to update that as well.}
  \item{offset.test}{A numeric vector of length equal to that of the test matrix, or \code{NULL}. Can be missing
  	for \code{setTestPredictors}.}
  \item{control}{An object inheriting from \code{\link{dbartsControl}}.}
  \item{column}{An integer or character string specifying which column of the predictor matrix is to be replaced.}
}
\details{
  A \code{dbartsSampler} is a \dQuote{mutable} object which contains information pertaining to
  fitting a Bayesian additive regression tree model. The sampler is first created and then,
  in a separate instruction, run or modified. In this way, MCMC samplers can be constructed
  with BART components filling arbitrary roles.
}
