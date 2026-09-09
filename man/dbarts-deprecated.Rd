\name{dbarts-deprecated}
\alias{dbarts-deprecated}
\alias{bart2}
\alias{rbart_vi}
\alias{predict.rbart}
\alias{extract.rbart}
\alias{fitted.rbart}
\alias{residuals.rbart}
\title{Deprecated and Removed Functions in dbarts}
\description{
  Names kept reachable for one release past their removal or rename, each pointing at its successor. They are removed in dbarts 1.1-0.

  \code{bart2} is the former name of the modern front door, which is now \code{\link{bart}}; a \code{bart2} call is forwarded unchanged, with a warning the first time in a session.

  \code{rbart_vi} and its methods fit a grouped random-intercept BART model. That model now lives in the \pkg{stan4bart} package (\code{stan4bart::stan4bart}), whose prior on the group spread differs from the one dbarts drew from, so a refit moves the results rather than reproducing them. These stubs error.
}
\details{
  The rest of the registry is names still reachable on entry points that keep their own page, rather than functions of their own: each is accepted for one release, mapped onto its successor with a once-per-session warning (or refused by name, where noted), and removed in dbarts 1.1-0.

  \describe{
    \item{\code{rngSeed}}{
      On \code{\link{bart}} and \code{\link{dbartsControl}}: the engine seed, now \code{seed}.
    }
    \item{\code{sigma}}{
      On \code{\link{dbarts}} and \code{\link{dbartsSpec}}: the creation-time residual estimate, now \code{sigest}.
    }
    \item{\code{resid.dist}, \code{dispersion}, \code{breaks}, \code{max.rows}}{
      On \code{\link{bart}}, \code{\link{dbarts}}, and (the first two) \code{\link{dbartsSpec}}: family-only settings, now arguments of the \code{\link{dbartsFamilies}} constructor that owns each (\code{family = student(df)}, \code{family = nbinom(dispersion)}, \code{family = hazard(breaks)}, \code{family = hazard(max.rows)}).
    }
    \item{\code{dart}, \code{levelGibbs}}{
      On \code{\link{bart}} (\code{dart} also on \code{\link{xbart}}): tree-prior-only settings, now arguments of \code{\link{dbartsPriors}}'s tree priors (\code{tree.prior = dart()}, \code{tree.prior = cgm(levelGibbs = )}).
    }
    \item{\code{control} on \code{\link{xbart}}}{
      Refused outright, naming the flat arguments its settings became (\code{n.cuts}, \code{useQuantiles}, \code{n.thin}, \code{storage}): \code{xbart} builds its own control, so there is nothing left for one to set.
    }
    \item{A three-element \code{n.burn} on \code{\link{xbart}}}{
      Refused outright, naming \code{n.burn = c(fresh, warm)}: 0.9-x's third element was a per-replication burn-in, and a chain is never carried between replications now.
    }
    \item{\code{family = "twopart"}}{
      Refused outright on every entry point that takes \code{family}, naming \code{family = "hurdle.lognormal"}.
    }
    \item{\code{dbartsSampler}'s \code{$startThreads}/\code{$stopThreads}}{
      No-op methods (\code{\linkS4class{dbartsSampler}}): the hierarchical thread manager they drove is gone, threads now owned per run, but a 0.9-x Gibbs loop that brackets its sweeps with them still runs.
    }
    \item{A BayesTree-spelled \code{\link{bart}} call}{
      Forwarded to \code{\link{bartBT}} after a once-per-session warning, when the call carries a name from \code{bartBT}'s own formals that \code{bart} does not take.
    }
    \item{A fourth or later positional \code{\link{bart}} argument}{
      Refused, naming both doors: 0.9-x's \code{bart} read a fourth positional argument as \code{sigest}, where \code{bart}'s own fourth formal is \code{subset}.
    }
    \item{The front-door startup message}{
      Attaching the package prints one \code{packageStartupMessage} naming both doors, once per session (\code{suppressPackageStartupMessages} silences it); no successor, since nothing replaces it.
    }
  }

  Every entry above, and \code{bart2} and the \code{rbart_vi} stubs, is one row of \code{dbarts:::dbartsTombstones}, the registry this page documents by hand.
}
\usage{
rbart_vi(\dots)

\method{predict}{rbart}(object, \dots)
\method{extract}{rbart}(object, \dots)
\method{fitted}{rbart}(object, \dots)
\method{residuals}{rbart}(object, \dots)
}
\arguments{
  \item{object}{
    An \code{rbart} fit object saved by dbarts 0.9-x.
  }
  \item{\dots}{
    Ignored; the stubs do not run. \code{bart2}'s own arguments are \code{\link{bart}}'s, copied rather than restated - \code{formals(bart2)} is \code{formals(bart)} - so a caller reading the defaults programmatically still finds them.
  }
}
\value{
  \code{bart2} returns whatever \code{\link{bart}} returns. The \code{rbart_vi} stubs do not return.
}
\seealso{
  \code{\link{bart}}, \code{\link{bartBT}}
}
\keyword{internal}
