\name{monotone}
\alias{monotone}
\title{Monotone Constraints for BART}
\description{
  Build a monotone-constraint specification making a BART fit monotone
  increasing or decreasing in chosen predictors (monotone BART; Chipman,
  George, McCulloch, and Shively 2022), and name the prior it is read under.
  Pass the result as the \code{monotone} argument of \code{\link{dbarts}},
  \code{\link{bart}} or \code{\link{dbartsSpec}}; a plain direction vector
  passed there is shorthand for \code{monotone(directions)} at the default
  prior.

  The constructor is not exported: it resolves by bare name inside the
  arguments that take it, and elsewhere is written
  \code{dbartsForests$monotone(...)}; see \code{\link{dbartsForests}}.
}
\usage{
monotone(directions, prior = c("leaf", "joint"))
}
\arguments{
  \item{directions}{
    The per-predictor directions. A named vector selects predictors by
    model-matrix column name, each element \code{"increasing"} or \code{1},
    \code{"decreasing"} or \code{-1}, or \code{0} for unconstrained. An
    unnamed vector as long as the model matrix is wide assigns the directions
    by position. Matching is case-sensitive and takes no abbreviations. A
    predictor named \code{prior} is written like any other,
    \code{monotone(c(prior = "increasing"))}.
  }
  \item{prior}{
    The prior the constraint is read under. \code{"leaf"} restricts the
    leaf-value prior to the leaf values that are monotone given the tree,
    leaving the tree prior unchanged. \code{"joint"} conditions the tree
    structure and its leaf values on monotonicity together; it is not
    available yet, and a fit naming it is an error.
  }
}
\details{
  Only numeric and ordered columns are eligible: a direction on a categorical
  (unordered factor) predictor is an error. Names, positions and every
  direction are validated against the model matrix at fit time; \code{prior}
  is matched when the specification is built. An all-zero specification fits
  the unconstrained model.

  A constraint forces birth/death-only tree proposals (a \code{control} naming
  a non-default \code{proposal.probs} is then an error) and a fixed
  \code{k = 2} (an explicit \code{k} hyperprior is an error); linear and
  Gaussian-process leaves, variance forests and multi-forest models are not
  supported under it. It may be combined with \code{\link{interactions}} and
  \code{\link{blocks}}.
}
\value{
  A \code{dbartsMonotone} specification object, a list of \code{directions}
  and \code{prior}, resolved when a sampler is built.
}
\references{
  Chipman, H. A., George, E. I., McCulloch, R. E., and Shively, T. S. (2022)
  mBART: multidimensional monotone BART. \emph{Bayesian Analysis}
  \bold{17}(2), 515--544.
}
\seealso{
  \code{\link{dbartsForests}},
  \code{\link{dbarts}}, \code{\link{bart}}, \code{\link{interactions}},
  \code{\link{blocks}}
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] - x[, 2] + sin(2 * pi * x[, 3]) + rnorm(n, 0, 0.2)
df <- data.frame(y, x)

## increasing in x1, decreasing in x2, free in x3
fit <- bart(y ~ x1 + x2 + x3, df,
             monotone = monotone(c(x1 = "increasing", x2 = "decreasing"),
                                 prior = "leaf"),
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## the plain vector, positionally, at the default prior
fit.plain <- bart(y ~ x1 + x2 + x3, df,
                   monotone = c(1, -1, 0),
                   n.trees = 25L, n.samples = 20L, n.burn = 20L,
                   n.chains = 1L, verbose = FALSE)
}
\keyword{models}
