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
monotone(directions, prior = c("joint", "leaf"))
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
    How the constraint enters the prior, \code{"joint"} (the default) or
    \code{"leaf"}; see Details.
  }
}
\details{
  Each tree is held monotone in the constrained predictors, so the sum of
  trees is too. The two priors differ in where the constraint acts.
  \describe{
    \item{\code{"joint"}}{
      The ordinary BART prior on the trees and their leaf values, conditioned
      on every tree being monotone. A tree whose leaf values would seldom be
      in order under the ordinary prior becomes less likely, so the
      constrained predictors get somewhat fewer splits than in an
      unconstrained fit.
    }
    \item{\code{"leaf"}}{
      The ordinary tree prior, unchanged; given a tree, its leaf values are
      drawn from the ordinary prior restricted to values in order. Every tree
      move then counts the orderings of the tree's leaves that the constraint
      allows.
    }
  }
  \code{"joint"} is the default because in a simulation study the two
  priors fit alike, and \code{"joint"} never counts leaf orderings, so it is
  never slow.

  Under \code{"leaf"} the count is negligible at 10 or more trees, but with
  one to five trees a tree can grow large enough that a count takes seconds
  to minutes per move, and its memory grows with it, to a gigabyte or more
  for the largest counts. There is no limit: the run finishes on the stated
  model, and afterwards a warning of class \code{"dbartsSlowCountWarning"}
  reports counts that took more than about a second and names the remedies,
  more trees or \code{prior = "joint"}. A count can be interrupted, and a
  count that runs out of memory is an ordinary error; either way the
  sampler is left in a valid state.

  The tree prior is the fit's \code{tree.prior}, as without a constraint:
  \code{cgm()} at its defaults, \code{power = 2} and \code{base = 0.95}
  (see \code{\link{dbartsPriors}}), unless another is given. The mBART paper used \code{cgm(power = 0.8, base =
  0.25)}, which gives smaller trees; in a simulation study under
  \code{"joint"} at 200 trees it gave sharper and more accurate fits, but
  narrower intervals that under-covered with plentiful, low-noise data. Pass
  it as \code{tree.prior} to opt in (see Examples).

  The fit is monotone in a constrained predictor along every line on which
  that predictor is observed, the other predictors at any values, missing
  ones and every factor level included; nothing is claimed where the
  constrained predictor itself is missing.

  Only numeric and ordered columns are eligible: a direction on a categorical
  (unordered factor) predictor is an error. Names, positions and every
  direction are validated against the model matrix at fit time; \code{prior}
  is matched when the specification is built. An all-zero specification fits
  the unconstrained model.

  A constraint forces birth/death-only tree proposals (a \code{control}
  naming a non-default \code{proposal.probs} is then an error). The leaf
  scale is fixed: \code{k} is 2 unless another number is given, or a fixed
  \code{sd} may be named, and a hyperprior on either is an error. Linear and
  Gaussian-process leaves, variance forests, multi-forest models and
  \code{family = "multinomial"} are not supported under it. It may be
  combined with \code{\link{interactions}} and \code{\link{blocks}}.
}
\value{
  A \code{dbartsMonotone} specification object, a list of \code{directions}
  and \code{prior}, resolved when a sampler is built. A fit records the prior
  as \code{monotone.prior}, which \code{print} and \code{summary} show.
}
\references{
  Chipman, H. A., George, E. I., McCulloch, R. E., and Shively, T. S. (2022)
  mBART: multidimensional monotone BART. \emph{Bayesian Analysis}
  \bold{17}(2), 515--544.
}
\seealso{
  \code{\link{dbartsForests}}, \code{\link{dbartsPriors}},
  \code{\link{dbarts}}, \code{\link{bart}}, \code{\link{interactions}},
  \code{\link{blocks}}
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] - x[, 2] + sin(2 * pi * x[, 3]) + rnorm(n, 0, 0.2)
df <- data.frame(y, x)

## increasing in x1, decreasing in x2, free in x3, under the default prior
fit <- bart(y ~ x1 + x2 + x3, df,
            monotone = monotone(c(x1 = "increasing", x2 = "decreasing")),
            n.trees = 25L, n.samples = 20L, n.burn = 20L,
            n.chains = 1L, verbose = FALSE)
fit$monotone.prior

## the plain vector, positionally, is the same specification
fit.plain <- bart(y ~ x1 + x2 + x3, df,
                  monotone = c(1, -1, 0),
                  n.trees = 25L, n.samples = 20L, n.burn = 20L,
                  n.chains = 1L, verbose = FALSE)

## the "leaf" prior
fit.leaf <- bart(y ~ x1 + x2 + x3, df,
                 monotone = monotone(c(x1 = "increasing", x2 = "decreasing"),
                                     prior = "leaf"),
                 n.trees = 25L, n.samples = 20L, n.burn = 20L,
                 n.chains = 1L, verbose = FALSE)

## the mBART paper's tree prior, as an opt-in
fit.mbart <- bart(y ~ x1 + x2 + x3, df,
                  monotone = c(x1 = "increasing", x2 = "decreasing"),
                  tree.prior = cgm(power = 0.8, base = 0.25),
                  n.trees = 25L, n.samples = 20L, n.burn = 20L,
                  n.chains = 1L, verbose = FALSE)
}
\keyword{models}
