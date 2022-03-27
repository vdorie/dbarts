\name{dbartsData}
\alias{dbartsData}
\title{Discrete Bayesian Additive Regression Trees Sampler Data}
\description{
  Convenience function to create a data object for use with a \code{\link{dbarts}} sampler.
}
\usage{
dbartsData(
    formula, data, test, subset, weights,
    offset, offset.test = offset)
}
\arguments{
   \item{formula, data, test, subset, weights, offset, offset.test}{
     As in \code{\link{dbarts}}. Retains backwards compatibility with \code{\link{bart}}, so that \code{formula}/\code{data} can be a \code{\link{formula}}/\code{\link{data.frame}} pair, or a pair of \code{x.train}/\code{y.train} matrices/vector.
   }
}
\value{
  An object of class \code{dbartData}.
}
\seealso{
  \code{\link{dbarts}}
}
