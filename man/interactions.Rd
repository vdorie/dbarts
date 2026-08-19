\name{interactions}
\alias{interactions}
\title{Interaction Constraints for BART}
\description{
  Build an interaction-constraint specification restricting which predictors
  may jointly shape a BART fit. Pass the result as the \code{interactions}
  argument of \code{\link{dbarts}} or \code{\link{bart2}}. The constraint is
  applied per forest, so in a multi-forest model one forest can be held
  additive-or-low-order while another stays free.
}
\usage{
interactions(max.order = NULL, groups = NULL, forbid = NULL)
}
\arguments{
  \item{max.order}{
    Optional cap on the number of \emph{distinct} predictors that may appear
    together on any single root-to-leaf path. \code{max.order = 1} makes the
    forest additive (each tree splits on one predictor); \code{max.order = 2}
    allows bounded pairwise interactions; and so on. A single integer of at
    least 1, or \code{NULL} for no cap.
  }
  \item{groups}{
    Optional co-occurrence allow-list: a list of predictor sets. Two named
    predictors may share a root-to-leaf path only if some group contains both;
    every other pair of named predictors is forbidden. Predictors named in no
    group are unconstrained. Each set is a character vector of model-matrix
    column names (a bare factor term name expands to its indicator columns) or
    a numeric vector of column indices. A single vector (rather than a list)
    is taken as one group.
  }
  \item{forbid}{
    Optional co-occurrence deny-list: a list whose entries each name two or
    more predictors barred from ever sharing a root-to-leaf path (an entry of
    more than two forbids every pair within it). Columns are named or indexed
    as for \code{groups}. A single vector (rather than a list) is taken as one
    forbidden set.
  }
}
\details{
  All three restrictions lower to a single per-node admissibility rule over the
  set of predictors already used by a node's ancestors: \code{max.order} bars a
  new predictor once the path's distinct count reaches the cap, while
  \code{groups} and \code{forbid} bar a predictor whose forbidden partner sits
  above it. The constraint is hard (an availability ban), applied inside the
  tree proposals and enforced on any warm start, so a constrained fit never
  visits a path that violates it; an unconstrained fit is unaffected.

  Names and indices resolve against the model matrix at fit time, where every
  value is validated: an unrecognized name, an empty group, a \code{max.order}
  below 1, and a reference to a column dropped from the design are each an
  error there. Supplying none of the three arguments is instead an error
  immediately, at construction: \code{interactions()} with no arguments never
  produces an object to pass to \code{\link{dbarts}}/\code{\link{bart2}}.

  Monotone constraints (see \code{\link{dbarts}}'s \code{monotone}) and
  interaction constraints are orthogonal and may be combined.

  A \code{groups} partition that is both \emph{disjoint} and \emph{total}
  (every predictor named in exactly one group) already makes the fit
  block-additive: the root split fixes each tree's group, and the
  co-occurrence rule confines every descendant to that same group, so the
  ensemble is exactly a sum of per-group functions
  \eqn{f = \sum_G f_G}{f = sum_G f_G} - a functional-ANOVA / grouped-GAMI
  decomposition that holds at every MCMC state. Under this idiom each tree's
  group is chosen adaptively (by its root draw each sweep), so the number of
  trees devoted to a group floats with the data. \code{\link{blocks}} gives
  the same decomposition with a \emph{fixed} per-group tree count and a
  static per-tree mask instead of an adaptive per-path allow-list; use it
  when a group needs guaranteed representation or a deterministic layout.
}
\value{
  A \code{dbartsInteractions} specification object, resolved when a sampler is
  built.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{bart2}}, \code{\link{blocks}}
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] + ifelse(x[, 2] > 0.5, 1, -1) + rnorm(n, 0, 0.2)
df <- data.frame(y, x)

## additive fit: no tree splits on more than one predictor on any path
fit.additive <- bart2(y ~ x1 + x2 + x3, df,
                      interactions = interactions(max.order = 1),
                      n.trees = 25L, n.samples = 20L, n.burn = 20L,
                      n.chains = 1L, verbose = FALSE)

## forbid x1 and x2 from ever interacting, but allow order-2 interactions
## among the rest
fit.forbid <- bart2(y ~ x1 + x2 + x3, df,
                    interactions = interactions(max.order = 2,
                                                forbid = list(c("x1", "x2"))),
                    n.trees = 25L, n.samples = 20L, n.burn = 20L,
                    n.chains = 1L, verbose = FALSE)
}
\keyword{models}
