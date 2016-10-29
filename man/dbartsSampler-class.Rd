\docType{class}
\name{dbartsSampler-class}
\title{Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbartsSampler}
\alias{dbartsSampler-class}
\alias{\S4method{run}{dbartsSampler}}
\alias{\S4method{copy}{dbartsSampler}}
\alias{\S4method{show}{dbartsSampler}}
\alias{\S4method{setControl}{dbartsSampler}}
\alias{\S4method{setModel}{dbartsSampler}}
\alias{\S4method{setData}{dbartsSampler}}
\alias{\S4method{setResponse}{dbartsSampler}}
\alias{\S4method{setOffset}{dbartsSampler}}
\alias{\S4method{setPredictor}{dbartsSampler}}
\alias{\S4method{setTestPredictor}{dbartsSampler}}
\alias{\S4method{setTestPredictorAndOffset}{dbartsSampler}}
\alias{\S4method{setTestOffset}{dbartsSampler}}
\alias{\S4method{printTrees}{dbartsSampler}}
\alias{\S4method{plotTree}{dbartsSampler}}
\description{
  A reference class object that contains a Bayesian Additive Regression Trees sampler in such a
  way that it can be modified, stopped, and started all while maintaining its own state.
}
\usage{
\S4method{run}{dbartsSampler}(numBurnIn, numSamples, updateState = NA)
\S4method{copy}{dbartsSampler}(shallow = FALSE)
\S4method{show}{dbartsSampler}()
\S4method{setControl}{dbartsSampler}(control)
\S4method{setModel}{dbartsSampler}(model)
\S4method{setData}{dbartsSampler}(data)
\S4method{setResponse}{dbartsSampler}(y, updateState = NA)
\S4method{setOffset}{dbartsSampler}(offset, updateState = NA)
\S4method{setPredictor}{dbartsSampler}(x, column, updateState = NA)
\S4method{setTestPredictor}{dbartsSampler}(x.test, column, updateState = NA)
\S4method{setTestPredictorAndOffset}{dbartsSampler}(x.test, offset.test, updateState = NA)
\S4method{setTestOffset}{dbartsSampler}(offset.test, updateState = NA)
\S4method{printTrees}{dbartsSampler}(treeNums)
\S4method{plotTree}{dbartsSampler}(treeNum, treePlotPars = list(nodeHeight = 12, nodeWidth = 40, nodeGap = 8), ...)
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
  \item{shallow}{A logical determining if the copy should retain the underlying data of the sampler
    (\code{TRUE}) or have its own copies (\code{FALSE}).}
  \item{control}{An object inheriting from \code{\link{dbartsControl}}.}
  \item{model}{An object inheriting from \code{dbartsModel}.}
  \item{data}{An object inheriting from \code{\link{dbartsData}}.}
  \item{y}{A numeric response vector of length equal to that with which the sampler was created.}
  \item{x}{A numeric predictor vector of length equal to that with which the sampler was created. Can be
  	an entirely matrix of new number of rows for \code{setTestPredictor}.}
  \item{x.test}{A new matrix of test predictors, of the number of columns equal to that in the current model.}
  \item{offset}{A numeric vector of length equal to that with which the sampler was created, or \code{NULL}.
  	If \code{offset.test} was set from \code{offset}, will attempt to update that as well.}
  \item{offset.test}{A numeric vector of length equal to that of the test matrix, or \code{NULL}. Can be missing
  	for \code{setTestPredictors}.}
  \item{column}{An integer or character string vector specifying which column/columns of the predictor matrix is
  	to be replaced. If missing, the entire matrix is substitude.}
  \item{treeNums}{An integer vector listing the indices of the trees to print.}
  \item{treeNum}{An integer listing the indices of the tree to plot.}
  \item{treePlotPars}{A list containing the number quantities \code{nodeHeight}, \code{nodeWidth}, and
    \code{nodeGap}, all of which control aspects of the resulting plot.}
  \item{...}{Extra arguments to \code{\link{plot}}.}
}
\details{
  A \code{dbartsSampler} is a \dQuote{mutable} object which contains information pertaining to
  fitting a Bayesian additive regression tree model. The sampler is first created and then,
  in a separate instruction, run or modified. In this way, MCMC samplers can be constructed
  with BART components filling arbitrary roles.
}
\value{
  For \code{run}, a named-list with contents \code{sigma}, \code{train}, \code{test}, and \code{varcount}.
  
  For \code{setPredictor}, \code{TRUE}/\code{FALSE} depending on whether or not the operation was successful.
  The operation can fail if the new predictor results in a tree with an empty leaf-node. If only single columns
  were replaced, on the update is rolled-back so that the sampler remains in a valid state.
}
