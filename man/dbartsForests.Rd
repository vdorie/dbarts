\name{dbartsForests}
\alias{dbartsForests}
\docType{data}
\title{Forest Specification Constructors}
\description{
  A list of the constructor functions building the forest specifications that
  the \code{monotone}, \code{interactions}, \code{blocks}, \code{variance} and
  \code{forests} arguments of \code{\link{dbarts}}, \code{\link{dbartsSpec}} and
  \code{\link{bart}} accept. The constructors are not exported, so that names
  like \code{forest} and \code{blocks} stay out of the search path, where
  another package could mask them or be masked depending on attach order.

  Inside those arguments, and inside a formula's \code{forest()} term, the
  constructors resolve by bare name: a call is always the constructor, and a
  bare name the caller has assigned a value to is that value. Anywhere else -
  a specification built ahead of the call, or an argument a wrapper forces
  before passing it on - write \code{dbartsForests$interactions(...)}.
  \code{forest()}'s first argument, the predictors a forest splits on, is
  read as terms when it names a predictor of the fit and otherwise takes the
  value it has where and when the constructor is called; see
  \code{\link{forest}}.
}
\format{
  A list of functions:
  \describe{
    \item{\code{interactions}}{
      Interaction constraints; see \code{\link{interactions}}.
    }
    \item{\code{blocks}}{
      Block-additive constraints; see \code{\link{blocks}}.
    }
    \item{\code{monotone}}{
      Monotone constraints and their prior; see \code{\link{monotone}}.
    }
    \item{\code{forest}}{
      One forest of a model, in a \code{forests} list or as a term of a
      formula; see \code{\link{forest}}.
    }
    \item{\code{varianceForest}}{
      The heteroscedastic variance forest; see \code{\link{varianceForest}}.
    }
  }
}
\seealso{
  \code{\link{dbarts}}, \code{\link{bart}}, \code{\link{dbartsPriors}},
  \code{\link{dbartsFamilies}}
}
\examples{
## built ahead of the call, then passed as a value
constraint <- dbartsForests$interactions(max.order = 1)

set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] + ifelse(x[, 2] > 0.5, 1, -1) + rnorm(n, 0, 0.2)

fit <- bart(x, y,
             interactions = constraint,
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## the same constraint, written by bare name inside the argument
fit.bare <- bart(x, y,
                  interactions = interactions(max.order = 1),
                  n.trees = 25L, n.samples = 20L, n.burn = 20L,
                  n.chains = 1L, verbose = FALSE)
}
\keyword{datasets}
