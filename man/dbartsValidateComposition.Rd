\name{dbartsValidateComposition}
\alias{dbartsValidateComposition}
\alias{print.dbartsCompositionValidation}
\title{Validate A Composed Sampler By Simulation-Based Calibration}
\description{
  Simulation-based calibration (Talts et al. 2018) over a one-sweep step \emph{you} supply, for a
  model \emph{you} state as a prior draw and a simulation. For each replication it draws
  \eqn{\theta_0}{theta0} from the prior, simulates data at it, runs the step, and ranks each scalar
  functional's value at \eqn{\theta_0}{theta0} among its values over the retained draws. Under a
  kernel that targets the posterior the prior and simulation imply, those ranks are uniform; a
  systematic departure means the composition does not sample what it claims to - the defect class in
  which a posterior \emph{mean} is reported where a draw belongs. It is a diagnostic, not a fix.
}
\usage{
dbartsValidateComposition(drawPrior, simulate, init, step, functionals,
                          n.replications = 200L, n.draws = 200L,
                          n.thin = 30L, n.burn = 200L, alpha = 0.05,
                          seed = NULL)
}
\arguments{
  \item{drawPrior}{
    A function of no arguments returning one draw \eqn{\theta_0}{theta0} from the prior, in whatever
    shape \code{simulate} and \code{init} expect. Nothing is read out of it by name.
  }
  \item{simulate}{
    A function of \eqn{\theta_0}{theta0} returning one simulated data set from the likelihood the
    composition assumes.
  }
  \item{init}{
    A function of \eqn{\theta_0}{theta0} and the simulated data returning the chain's initial state.
    It \strong{must} start the chain \emph{at} \eqn{\theta_0}{theta0}: the rank is taken of
    \code{functionals(init(theta0, data))}, so a state built anywhere else ranks the wrong quantity.
    The state carries whatever \code{step} and \code{functionals} need, the data included.
  }
  \item{step}{
    A function of the state returning the state after \strong{one} sweep of the composed kernel -
    the whole Gibbs/Metropolis-Hastings scan, not one block of it.
  }
  \item{functionals}{
    A function of the state returning a numeric vector of scalar summaries to rank, the same ones in
    the same order at every state. They may be \emph{derived} - a mean fit, a probability, an
    indicator - and need not be parameters; names, if given, label the report.
  }
  \item{n.replications}{
    Number of prior draw / simulation / refit replications. The rank resolution and the band both
    come from this.
  }
  \item{n.draws}{
    Number of draws retained per replication, \eqn{L}{L}; ranks run over \eqn{0}{0} to \eqn{L}{L}.
  }
  \item{n.thin}{
    Sweeps per retained draw. SBC assumes the retained draws are independent, so this should exceed
    the composition's autocorrelation time; \code{1} is right only for a step that draws the
    posterior outright.
  }
  \item{n.burn}{
    Sweeps discarded before the first retained draw. The default is not \code{0}: an \code{init} that
    is only approximately at \eqn{\theta_0}{theta0} would otherwise be flagged for its early draws.
    An \code{init} that is exact can use \code{0}.
  }
  \item{alpha}{
    Level of the uniformity band, Bonferroni-corrected within the call over the number of functionals
    ranked, so that a whole call passes with probability about \eqn{1 - \alpha}{1 - alpha}.
  }
  \item{seed}{
    Seeds R's random number stream for the whole run, making it reproducible; \code{NULL} uses the
    stream as it stands. The band's own stream is fixed internally and restored on exit either way,
    so a call never leaves that fix behind.
  }
}
\details{
  \strong{The verdict.} The primary statistic is the maximum absolute difference between the ranks'
  empirical cdf and the uniform one, against a simultaneous band simulated from the null (Talts et
  al. 2018, fig. 1), which is already corrected for the multiple looks across the rank grid. A
  chi-square goodness of fit and a jittered Kolmogorov-Smirnov test are reported beside it as
  secondary signals; the verdict does not hinge on them.

  \strong{Ties.} A rank is the number of draws below \eqn{\theta_0}{theta0}, which is uniform only
  for an atomless law. Every functional is ranked with the tie-break instead: draws and
  \eqn{\theta_0}{theta0} alike carry an independent uniform tag, so an atom - an indicator, a
  discrete parameter, a probability that underflows to zero - contributes a uniform increment rather
  than piling every replication onto one rank. With no ties it is the plain count, exactly.

  \strong{The blind spot.} The validator compares the composition against the model \code{drawPrior}
  and \code{simulate} state \emph{between them}. It cannot detect a host whose prior draw and
  simulation disagree with each other: a step that correctly targets the posterior of the wrong
  model passes, because that model is the one being validated against. Nor does it localize a
  failure to a block of the sweep; it reports which functionals are miscalibrated, not why.
}
\value{
  An object of class \code{dbartsCompositionValidation}, a list with components \code{ranks} (the
  \code{n.replications} by number-of-functionals integer matrix of ranks, columns named for the
  functionals), \code{L}, \code{n.replications}, \code{n.thin}, \code{n.burn}, \code{alpha},
  \code{alpha.adjusted} (the Bonferroni-corrected level actually used), \code{mean.rank.target}
  (\eqn{L/2}{L/2}), \code{pass}, and \code{verdicts}, a data frame of one row per functional carrying
  its mean rank, ecdf-difference statistic, band, chi-square and KS p-values, and its \code{"PASS"}
  or \code{"FLAG"} verdict. \code{print} reports the table.
}
\references{
  Talts, S., Betancourt, M., Simpson, D., Vehtari, A., and Gelman, A. (2018) Validating Bayesian
  inference algorithms with simulation-based calibration. \emph{arXiv preprint} arXiv:1804.06788.
}
\seealso{
  \code{\link[=dbartsSampler-class]{dbartsSampler}}, the sampler a composition drives, and
  \code{\link{dbartsAugmentation}}, the augmentation draws a composed step is built from.
}
\examples{
# a conjugate normal-normal one-sweep step: theta ~ N(0, 1) and
# y_i | theta ~ N(theta, 1) give theta | y ~ N(v * sum(y), v), v = 1 / (1 + n)
n <- 5L
v <- 1 / (1 + n)

drawPrior <- function() c(theta = rnorm(1L))
simulateData <- function(theta0) rnorm(n, theta0[["theta"]], 1)
initState <- function(theta0, data)
  list(theta = theta0[["theta"]], mean = v * sum(data), sd = sqrt(v))
functionals <- function(state)
  c(theta = state$theta, positive = as.numeric(state$theta > 0))

# the step draws the posterior outright, so its states are independent and its
# init is exact at theta0: this call needs neither thinning nor burn-in
stepExact <- function(state) {
  state$theta <- rnorm(1L, state$mean, state$sd)
  state
}
dbartsValidateComposition(drawPrior, simulateData, initState, stepExact,
                          functionals, n.replications = 50L, n.draws = 100L,
                          n.thin = 1L, n.burn = 0L, seed = 8L)

# the defect class this exists to catch: the posterior MEAN reported where a
# draw belongs. 'theta' flags; 'positive' is coarse enough to survive it, power
# being per functional.
stepMean <- function(state) {
  state$theta <- state$mean
  state
}
dbartsValidateComposition(drawPrior, simulateData, initState, stepMean,
                          functionals, n.replications = 50L, n.draws = 100L,
                          n.thin = 1L, n.burn = 0L, seed = 8L)
}
\keyword{models}
