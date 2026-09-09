\name{dbartsFamilies}
\alias{dbartsFamilies}
\docType{data}
\title{Response Family Constructors}
\description{
  A list of constructor functions building the response families that the
  \code{family} argument of \code{\link{bart}}, \code{\link{dbarts}},
  \code{\link{dbartsSpec}} and \code{\link{xbart}} accepts. Bundling them
  keeps generic names like \code{gaussian} and \code{probit} out of the
  search path, where another package - \code{stats::gaussian} first among
  them - could mask them or be masked depending on attach order.

  \code{family} takes either a token string (\code{family = "gaussian"}),
  meaning that family at its defaults, or one of these calls
  (\code{family = student(df = 4)}), meaning that family with settings.
  Inside the \code{family} argument the constructors resolve by bare name
  whatever the caller has attached, exactly as the prior vocabulary does
  inside \code{tree.prior}; elsewhere reach them through this list.

  Every setting that only one family reads rides its family object rather
  than a formal of the fitting function. The retired spellings
  \code{resid.dist}, \code{dispersion}, \code{breaks} and \code{max.rows}
  are accepted for one release with a once-per-session warning and are
  removed in dbarts 1.1-0.
}
\format{
  A list of functions:
  \describe{
    \item{\code{gaussian()}}{
      A continuous response with normal errors; the default for a numeric
      response under \code{family = "auto"}.
    }
    \item{\code{student(df = NULL)}}{
      A continuous response with outlier-robust Student-t errors, drawn by
      the Gaussian scale-mixture augmentation. \code{df = NULL} estimates
      the degrees of freedom on a capped grid; a positive number fixes
      them. The former \code{resid.dist = student(df)} spelling.
    }
    \item{\code{probit()}, \code{logistic()}}{
      Binary latent-variable models for a 0/1 response.
    }
    \item{\code{multinomial()}, \code{ordinal()}}{
      Unordered K-category softmax and ordered cumulative-probit models.
    }
    \item{\code{nbinom(dispersion = NA)}}{
      Negative-binomial counts. \code{NA} estimates the dispersion
      \eqn{r} on a capped positive-integer grid; a positive integer fixes
      it. The former \code{dispersion} argument.
    }
    \item{\code{aft()}}{
      Accelerated failure time (log-normal) survival.
    }
    \item{\code{hazard(breaks = NULL, max.rows = 1e7,
                       link = c("probit", "logistic"))}}{
      Discrete-time survival hazard by person-period expansion.
      \code{breaks} sets the period grid: \code{NULL} uses the sorted
      distinct observed times, a single positive integer bins at that many
      quantiles, and a numeric boundary vector gives explicit right-closed
      intervals. \code{max.rows} refuses an over-large expansion, naming
      the coarsening levers. \code{link} selects the binary link the
      expanded rows are fit under, and is what the tokens
      \code{"hazard.probit"} and \code{"hazard.logistic"} name. The former
      \code{breaks} and \code{max.rows} arguments.
    }
    \item{\code{hurdle.lognormal()}}{
      A semicontinuous two-part response: an occupancy probit glued to a
      lognormal positive part. Composed from two samplers, so only
      \code{\link{bart}} fits it.
    }
  }
}
\details{
  Which entry points accept which family:

  \tabular{lllll}{
    \strong{family} \tab \strong{bart} \tab \strong{dbarts} \tab
      \strong{dbartsSpec} \tab \strong{xbart} \cr
    \code{gaussian} \tab yes \tab yes \tab yes \tab yes \cr
    \code{student} \tab yes \tab yes \tab yes \tab no \cr
    \code{probit} \tab yes \tab yes \tab yes \tab yes \cr
    \code{logistic} \tab yes \tab yes \tab yes \tab yes \cr
    \code{aft} \tab yes \tab yes \tab yes \tab no \cr
    \code{multinomial} \tab yes \tab yes \tab yes \tab no \cr
    \code{ordinal} \tab yes \tab yes \tab yes \tab no \cr
    \code{nbinom} \tab yes \tab yes \tab yes \tab no \cr
    \code{hazard} \tab yes \tab yes \tab no \tab no \cr
    \code{hurdle.lognormal} \tab yes \tab no \tab no \tab no
  }

  A family the entry point does not fit is refused by name. See
  \code{\link{dbarts}}'s \code{family} argument for what each family
  models and what it supports.
}
\seealso{
  \code{\link{bart}}, \code{\link{dbarts}}, \code{\link{dbartsPriors}}
}
\examples{
f <- dbartsFamilies$student(df = 4)
f

set.seed(99)
n <- 100L
x <- matrix(rnorm(n * 2L), n, 2L)
y <- x[, 1L] + rt(n, 4)

fit <- bart(
  x, y,
  family = student(df = 4),
  n.samples = 50L, n.burn = 50L, n.chains = 1L, verbose = FALSE
)
}
\keyword{datasets}
