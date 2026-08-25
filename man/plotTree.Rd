\name{plotTree}
\alias{plotTree}
\alias{plotTree.bart}
\alias{plotTree.rbart}
\alias{plotTree.dbartsSampler}
\alias{plotTree.bartMultinomial}
\alias{plotTree.bartOrdinal}
\alias{plotTree.bartNegbin}
\alias{plotTree.bartHurdle}
\title{Plot a Single Tree From a Fitted BART Model}
\description{
  Minimalist visualization of the branching and leaf contents of one tree
  in a fitted \code{\link{bart}}/\code{\link{bart2}} or
  \code{\link{rbart_vi}} model, or in a \code{\linkS4class{dbartsSampler}}.
  A fit-level convenience wrapper around the sampler's \code{plotTree}
  method, so the trees can be plotted without reaching into the sampler
  stored on the fit.

  \code{\link{bart2}}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}) refuse
  by name instead: their trees live on the sampler(s) those fits carry
  (\code{object$fit} for the first three, \code{object$occupancy$fit}/
  \code{object$positive$fit} for the hurdle), which \code{plotTree} can
  be called on directly.
}
\usage{
plotTree(object, \dots)

\method{plotTree}{bart}(object, treeNum = 1L, chainNum, sampleNum, \dots)

\method{plotTree}{rbart}(object, treeNum = 1L, chainNum = 1L, sampleNum, \dots)

\method{plotTree}{dbartsSampler}(object, \dots)

\method{plotTree}{bartMultinomial}(object, \dots)

\method{plotTree}{bartOrdinal}(object, \dots)

\method{plotTree}{bartNegbin}(object, \dots)

\method{plotTree}{bartHurdle}(object, \dots)
}
\arguments{
  \item{object}{
    A fitted model of class \code{bart} (from \code{\link{bart}} or
    \code{\link{bart2}}) or \code{rbart} (from \code{\link{rbart_vi}}), or a
    \code{\linkS4class{dbartsSampler}}. Fits must have been made with the
    trees kept (\code{keeptrees}/\code{keepTrees} equal to \code{TRUE}).
  }
  \item{treeNum}{
    An integer, the index of the tree to plot.
  }
  \item{chainNum}{
    An integer, the index of the chain to plot from. For a \code{bart} fit
    with a single chain it may be omitted; for \code{rbart} it selects the
    chain and defaults to the first.
  }
  \item{sampleNum}{
    An integer, the index of the saved sample to plot from. When omitted,
    the last sample is used (or the current working trees when the trees
    were not kept as samples).
  }
  \item{\dots}{
    Additional arguments, including \code{treePlotPars} (a named numeric
    vector controlling \code{nodeHeight}, \code{nodeWidth}, and
    \code{nodeGap}) and further arguments passed on to \code{\link{plot}}.
  }
}
\value{
  Invisibly returns \code{NULL}; called for the side effect of drawing.
}
\details{
  The rendering handles constant, linear, and Gaussian-process leaves. For a
  multiple-chain fit, \code{chainNum} selects the chain. Equivalent to
  calling the \code{plotTree} method on the underlying sampler, i.e.
  \code{object$fit$plotTree(\dots)} for a \code{bart} fit.
}
\seealso{
  \code{\link{bart}}, \code{\link{rbart_vi}}, \code{\linkS4class{dbartsSampler}}
}
\examples{
f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

fit <- bart2(x, y, n.samples = 40L, n.burn = 10L, n.trees = 25L,
             n.chains = 1L, keepTrees = TRUE, verbose = FALSE)

plotTree(fit, treeNum = 1L, sampleNum = 40L)
}
\keyword{hplot}
