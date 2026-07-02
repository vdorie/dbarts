\docType{class}
\name{dbartsSampler-class}
\title{Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbartsSampler}
\alias{dbartsSampler-class}
\alias{\S4method{run}{dbartsSampler}}
\alias{\S4method{sampleTreesFromPrior}{dbartsSampler}}
\alias{\S4method{sampleNodeParametersFromPrior}{dbartsSampler}}
\alias{\S4method{copy}{dbartsSampler}}
\alias{\S4method{show}{dbartsSampler}}
\alias{\S4method{predict}{dbartsSampler}}
\alias{\S4method{setControl}{dbartsSampler}}
\alias{\S4method{setModel}{dbartsSampler}}
\alias{\S4method{setData}{dbartsSampler}}
\alias{\S4method{setResponse}{dbartsSampler}}
\alias{\S4method{setOffset}{dbartsSampler}}
\alias{\S4method{setSigma}{dbartsSampler}}
\alias{\S4method{setPredictor}{dbartsSampler}}
\alias{\S4method{setCutPoints}{dbartsSampler}}
\alias{\S4method{setTestPredictor}{dbartsSampler}}
\alias{\S4method{setTestPredictorAndOffset}{dbartsSampler}}
\alias{\S4method{setTestOffset}{dbartsSampler}}
\alias{\S4method{printTrees}{dbartsSampler}}
\alias{\S4method{getTrees}{dbartsSampler}}
\alias{\S4method{plotTree}{dbartsSampler}}
\alias{\S4method{startThreads}{dbartsSampler}}
\alias{\S4method{stopThreads}{dbartsSampler}}
\description{
  A reference class object that contains a Bayesian Additive Regression Trees sampler in such a way that it can be modified, stopped, and started all while maintaining its own state.
}
\usage{
\S4method{run}{dbartsSampler}(
  numBurnIn, numSamples, updateState = NA, n.threads = control@n.threads
)
\S4method{sampleTreesFromPrior}{dbartsSampler}(updateState = NA)
\S4method{sampleNodeParametersFromPrior}{dbartsSampler}(updateState = NA)
\S4method{copy}{dbartsSampler}(shallow = FALSE)
\S4method{show}{dbartsSampler}()
\S4method{predict}{dbartsSampler}(x.test, offset.test, n.threads = control@n.threads)
\S4method{setControl}{dbartsSampler}(control)
\S4method{setModel}{dbartsSampler}(model)
\S4method{setData}{dbartsSampler}(data)
\S4method{setResponse}{dbartsSampler}(y, updateState = NA)
\S4method{setOffset}{dbartsSampler}(offset, updateScale = FALSE, updateState = NA)
\S4method{setSigma}{dbartsSampler}(sigma, updateState = NA)
\S4method{setPredictor}{dbartsSampler}(x, column, forceUpdate, updateCutPoints = FALSE, updateState = NA)
\S4method{setCutPoints}{dbartsSampler}(cuts, column, updateState = NA)
\S4method{setTestPredictor}{dbartsSampler}(x.test, column, updateState = NA)
\S4method{setTestPredictorAndOffset}{dbartsSampler}(x.test, offset.test, updateState = NA)
\S4method{setTestOffset}{dbartsSampler}(offset.test, updateState = NA)
\S4method{printTrees}{dbartsSampler}(treeNums)
\S4method{getTrees}{dbartsSampler}(
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL
)
\S4method{plotTree}{dbartsSampler}(
  treeNum, treePlotPars = c(
    nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
  ...
)
\S4method{startThreads}{dbartsSampler}(n.threads = control@n.threads)
\S4method{stopThreads}{dbartsSampler}()
}
\arguments{
  \item{numBurnIn}{
    A non-negative integer determining how many iterations the sampler should skip before storing results. If missing or \code{NA}, the default is filled in from the sampler's \code{\link[=dbartsControl]{control}} object.
  }
  \item{numSamples}{
    A positive integer determining how many posterior samples should be returned. If missing or \code{NA}, the default is also filled in from the control object.
  }
  \item{updateState}{
    A logical determining if the local cache of the sampler's state should be updated after the completion of the run. If \code{NA}, the default is also filled in from the control object.
  }
  \item{shallow}{
    A logical determining if the copy should retain the underlying data of the sampler (\code{TRUE}) or have its own copies (\code{FALSE}).
  }
  \item{control}{
    An object inheriting from \code{\link{dbartsControl}}.
  }
  \item{model}{
    An object inheriting from \code{dbartsModel}.
  }
  \item{data}{
    An object inheriting from \code{\link{dbartsData}}.
  }
  \item{y}{
    A numeric response vector of length equal to that with which the sampler was created.
  }
  \item{x}{
    A numeric predictor vector of length equal to that with which the sampler was created. Can be of a distinct number of rows for \code{setTestPredictor}.
  }
  \item{x.test}{
    A new matrix of test predictors, of the number of columns equal to that in the current model.
  }
  \item{offset}{
    A numeric vector of length equal to that with which the sampler was created, or \code{NULL}. If \code{offset.test} was set from \code{offset}, will attempt to update that as well.
  }
  \item{updateScale}{
    Logical indicating whether BART's internal scale should update with the new offset. Should only be \code{TRUE} during burn-in.
  }
  \item{offset.test}{
    A numeric vector of length equal to that of the test matrix, or \code{NULL}. Can be missing for \code{setTestPredictors}.
  }
  \item{n.threads}{
    If greater than one, chain predictions will take place in parallel.
  }
  \item{sigma}{
    Numeric vector of residual standard deviations, one for each chain.
  }
  \item{cuts}{
    Vector of cut points for use with \code{setCutPoints}.
  }
  \item{column}{
    An integer or character string vector specifying which column/columns of
    the predictor matrix is to be replaced or for which cuts should be used.
    When updating predictors, it can be missing and the entire matrix will be
    substituted.
  }
  \item{forceUpdate}{
    For \code{setPredictor}, controls how the trees respond to a new predictor
    that would leave a tree with an empty leaf-node. \code{FALSE} rolls the whole
    update back and leaves the sampler unchanged (rejection sampling); \code{TRUE}
    forces the change, collapsing any empty nodes; the string \code{"partial"}
    installs each observation's new value individually and rolls back only those
    observations whose value would empty a leaf, returning a per-observation
    logical of what was installed. \code{"partial"} requires a single \code{column}
    and cannot be combined with \code{updateCutPoints}. When missing, defaults to
    \code{TRUE} if the replacement predictor is a matrix and \code{FALSE}
    otherwise.
  }
  \item{updateCutPoints}{
    For \code{setPredictor}, a logical; when \code{TRUE} the cut points (split
    candidate locations) for the replaced column(s) are recomputed from the new
    values, otherwise the existing cut points are kept. Defaults to \code{FALSE}
    and cannot be combined with \code{forceUpdate = "partial"}.
  }
  \item{treeNums}{
    An integer vector listing the indices of the trees to print or return.
  }
  \item{chainNums}{
    An integer vector listing the chains to return from \code{getTrees}. When
    missing, all chains are returned.
  }
  \item{sampleNums}{
    An integer vector listing the saved samples to return from \code{getTrees}.
    Applies only when \code{keepTrees} is \code{TRUE} and \code{current} is
    \code{FALSE}; otherwise the live working trees have no sample dimension and
    it is ignored.
  }
  \item{current}{
    For \code{getTrees}, a logical; when \code{TRUE} the live working trees (the
    current sampler position) are returned even if \code{keepTrees} is
    \code{TRUE}, rather than the saved samples.
  }
  \item{newdata}{
    For \code{getTrees}, an optional matrix of predictors (in the form accepted
    by \code{predict}). When supplied, its observations are routed through each
    tree so the \code{n} column counts them instead of the training data; the
    tree structure, splits, and leaf values are unchanged.
  }
  \item{treeNum}{
    An integer listing the indices of the tree to plot.
  }
  \item{treePlotPars}{
    A named numeric vector containing the quantities \code{nodeHeight}, \code{nodeWidth}, and \code{nodeGap}, all of which control aspects of the resulting plot.
  }
  \item{...}{
    Extra arguments to \code{\link{plot}}.
  }
}
\details{
  A \code{dbartsSampler} is a mutable object which contains information pertaining to fitting a Bayesian additive regression tree model. The sampler is first created and then, in a separate instruction, run or modified. In this way, MCMC samplers can be constructed with BART components filling arbitrary roles.

  \subsection{Saving}{
    \code{\link{save}}-ing and \code{\link{load}}ing a \code{dbarts} sampler for future use requires that R's serialization mechanism be able to access the state of the sampler which, for memory purposes, is only made available to R on request. To do this, one must \dQuote{touch} the sampler's state object before saving, e.g. for the object \code{sampler}, execute \code{invisible(sampler$state)}. This is in addition to guaranteeing that the \code{state} object is not \code{NULL}, which can be done by setting the sampler's control to an object with \code{updateState} as \code{TRUE} or passing \code{TRUE} as the \code{updateState} argument to any of the sampler's applicable methods.
  }
}
\value{
  For \code{run}, a named-list with contents \code{sigma}, \code{train}, \code{test}, and \code{varcount}.

  For \code{setPredictor}, \code{TRUE}/\code{FALSE} depending on whether or not the operation was successful. The operation can fail if the new predictor results in a tree with an empty leaf-node. If only single columns were replaced, on the update is rolled-back so that the sampler remains in a valid state. When \code{forceUpdate} is \code{"partial"}, instead returns a logical vector of length equal to the number of observations, \code{TRUE} where that observation's new value was installed and \code{FALSE} where it was rolled back to its previous value to keep every tree valid.

  \code{predict} keeps the current test matrix in place and uses the current set of tree splits. This function has two use cases. The first is when \code{keepTrees} of \code{\link{dbartsControl}} is \code{TRUE}, in which case the sampler should be run to completion and the function can be used to interrogate the existing fit. When \code{keepTrees} is \code{FALSE}, the function can be used to obtain the likelihood as part of a proposed new set of covariates in a Metropolis-Hastings step in a full-Bayes sampler. This would typically be followed by a call to \code{setPredictor} if the step is accepted.

  For \code{getTrees}, a \code{data.frame} with one row per tree node in depth-first, left-hand-side pre-order, with columns \code{chain}, \code{sample} (present only for saved samples), \code{tree}, \code{n} (the number of observations in the node), \code{var} (the splitting variable, or -1 at a leaf), and \code{value} (the split value, or the leaf prediction).
}
\seealso{
  \code{\link{updatePredictorPerObservationJointly}} for applying a shared per-observation predictor update across several samplers at once.
}
