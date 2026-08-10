\name{treatmentForest}
\alias{treatmentForest}
\title{Treatment Forest Specification for Bayesian Causal Forests}
\description{
  Build the specification of the treatment forest of a Bayesian causal forest.
  Pass the result as the \code{treatmentForest} argument of
  \code{\link{dbarts}} or \code{\link{dbartsSpec}}, alongside a
  \code{treatment} vector. One constructor carries every knob the second forest
  and its glue own, so the fitting functions grow three arguments rather than
  fourteen.
}
\usage{
treatmentForest(
    n.trees = 50L, base = 0.25, power = 3,
    sd.control = 2, sd.moderate = 1, b.prior.variance = 0.5,
    update.a = TRUE, update.b = TRUE,
    interactions = NULL, blocks = NULL)
}
\arguments{
  \item{n.trees}{
    Number of trees in the treatment forest. A single integer of at least 1.
    Fewer trees than the prognostic forest is the usual choice: the treatment
    surface is normally the smoother of the two.
  }
  \item{base, power}{
    The treatment forest's tree-structure prior, as \code{\link{dbarts}}'s
    \code{tree.prior} is for the prognostic forest. The defaults are shallower
    than the prognostic forest's, for the same reason.
  }
  \item{sd.control, sd.moderate}{
    Prior scales, in standard deviations of the response. \code{sd.control} is
    the half-Cauchy median of the prognostic scalar \eqn{a}{a}, placing the
    prognostic total at \code{sd.control} response standard deviations;
    \code{sd.moderate} places the treatment effect
    \eqn{(b_1 - b_0)\tau}{(b1 - b0) tau} at \code{sd.moderate} of them.
  }
  \item{b.prior.variance}{
    Prior variance of the glue coefficients \eqn{b_0}{b0} and \eqn{b_1}{b1},
    which scale the treatment forest for control and treated observations.
  }
  \item{update.a, update.b}{
    Whether the corresponding glue block is redrawn each sweep. \code{FALSE}
    fixes it at its prior center for the sampler's life; both are fixed at
    creation and cannot be toggled afterwards.
  }
  \item{interactions, blocks}{
    Optional \code{\link{interactions}} and \code{\link{blocks}} constraints on
    the treatment forest. The arguments of the same names on
    \code{\link{dbarts}} constrain the prognostic forest, so the two can be
    held to different structures - the calibrated-additivity idiom is an
    additive or low-order treatment forest beside a free prognostic one. A
    \code{blocks} partition here covers the columns the forest may split on,
    i.e. the \code{moderators} subset when one is given.
  }
}
\details{
  A Bayesian causal forest fits
  \eqn{y = a\,\mu(x) + b_z\,\tau(x) + \epsilon}{y = a mu(x) + b_z tau(x) +
  eps}: a prognostic forest \eqn{\mu}{mu} over every predictor, and a treatment
  forest \eqn{\tau}{tau} over the moderators, joined by a three-parameter glue
  \eqn{(a, b_0, b_1)}{(a, b0, b1)}. This constructor describes the second
  forest and the glue; the prognostic forest takes the fitting function's own
  \code{control@n.trees} and \code{tree.prior}.

  Both forests' leaf scales come from the model's own calibration map rather
  than from the node prior, which is why a \code{k} hyperprior, a non-default
  \code{k}, and a linear or Gaussian-process node prior are refused when a
  treatment vector is supplied. Every value here is validated at fit time.
}
\value{
  A \code{dbartsTreatmentForest} specification object, resolved when a sampler
  is built.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{dbartsSpec}}, \code{\link{interactions}},
  \code{\link{blocks}}
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
z <- rbinom(n, 1L, 0.5)
y <- 2 * x[, 1] + z * (1 + 2 * x[, 3]) + rnorm(n, 0, 0.2)

sampler <- dbarts(x, y, treatment = z, moderators = c("x1", "x3"),
                  treatmentForest = treatmentForest(n.trees = 25L,
                                                    sd.moderate = 1.5),
                  control = dbartsControl(n.chains = 1L, n.trees = 25L,
                                          n.samples = 20L, n.burn = 20L))
samples <- sampler$run(20L, 20L)
}
\keyword{models}
