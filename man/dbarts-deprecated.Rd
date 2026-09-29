\name{dbarts-deprecated}
\alias{dbarts-deprecated}
\alias{bart2}
\title{Deprecated and Removed Functions in dbarts}
\description{
  Names kept reachable for one release past their removal or rename, each pointing at its successor. They are removed in dbarts 1.1-0.

  \code{bart2} is the former name of the modern front door, which is now \code{\link{bart}}; a \code{bart2} call is forwarded unchanged, with a warning the first time in a session.

  \code{rbart_vi} and its methods are deprecated, and documented on \code{\link{rbart_vi}}.
}
\details{
  The rest of the registry is names still reachable on entry points that keep their own page, rather than functions of their own: each is accepted for one release, mapped onto its successor with a once-per-session warning (or refused by name, where noted), and removed in dbarts 1.1-0.

  \describe{
    \item{\code{rngSeed}}{
      On \code{\link{bart}} and \code{\link{dbartsControl}}: the engine seed, now \code{seed}.
    }
    \item{\code{sigma}}{
      On \code{\link{dbarts}} and \code{\link{xbart}}: the creation-time residual estimate, now \code{sigest}; like \code{sigest}, it defaults to \code{NULL}.
    }
    \item{\code{node.prior}}{
      On \code{\link{dbarts}}: the leaf-value prior, now \code{leaf.prior}.
    }
    \item{\code{degreesOfFreedom}}{
      On \code{\link{dbartsPriors}}'s \code{chi()}: the degrees of freedom of the hyperprior on \code{k}, now \code{df}, as \code{invchi()} and base R's \code{dchisq} spell it. Positional calls, \code{chi(1.5, 2)}, are unchanged.
    }
    \item{\code{power}, \code{base}}{
      On \code{\link{bart}}: tree-prior-only settings, now arguments of \code{\link{dbartsPriors}}'s \code{cgm()} or \code{dart()} (\code{tree.prior = cgm(power, base)} or \code{tree.prior = dart(power, base)}).
    }
    \item{\code{split.probs}}{
      On \code{\link{bart}}: a tree-prior-only setting, now an argument of \code{\link{dbartsPriors}}'s \code{cgm()} only (\code{tree.prior = cgm(split.probs = )}); a DART tree prior draws its own split probabilities, so combining it with \code{tree.prior = dart()} is refused.
    }
    \item{\code{resid.prior}, \code{sigdf}, \code{sigquant}}{
      The residual scale's prior is a setting of the families that draw one, so it rides the family object (\code{family = gaussian(sigma = chisq(df, quant))}, \code{family = gaussian(sigma = fixed(value))}); likewise \code{student(df, sigma = )}, \code{aft(sigma = )} and \code{hurdle.lognormal(sigma = )}. \code{sigdf} and \code{sigquant} are retired on \code{\link{bart}}; \code{resid.prior} is retired the same way on \code{\link{bart}}, \code{\link{dbarts}} and \code{\link{xbart}}, and \code{dbartsSpec} refuses it, so the prior has one home on every door. A call that writes it both ways - a retired spelling and a family whose call named \code{sigma} - is refused where the two disagree. \code{sigest}, the residual-scale estimate supplied at creation, is unaffected.
    }
    \item{\code{proposal.probs}}{
      On \code{\link{bart}} and \code{\link{dbarts}}: the tree-move mixture, now a \code{\link{dbartsControl}} slot and formal (\code{control = dbartsControl(proposal.probs = )}).
    }
    \item{A three-element \code{n.burn} on \code{\link{xbart}}}{
      Refused outright, naming \code{n.burn = c(fresh, warm)}: 0.9-x's third element was a per-replication burn-in, and a chain is never carried between replications now.
    }
    \item{\code{dbartsSampler}'s \code{$startThreads}/\code{$stopThreads}}{
      No-op methods (\code{\linkS4class{dbartsSampler}}): the hierarchical thread manager they drove is gone, threads now owned per run, but a 0.9-x Gibbs loop that brackets its sweeps with them still runs.
    }
    \item{\code{dbartsSampler}'s \code{$sampleNodeParametersFromPrior}}{
      Forwarded to \code{$sampleLeafParametersFromPrior} (\code{\linkS4class{dbartsSampler}}) after a once-per-session warning.
    }
    \item{A thread count on \code{dbartsSampler}'s \code{$run}}{
      Ignored after a once-per-session warning, whether passed as \code{n.threads}, as 0.9-x's formal \code{numThreads}, or positionally after \code{updateState}: a run uses the sampler's own count, which \code{$setControl} changes, and its draws do not depend on it.
    }
    \item{\code{NA} where \code{NULL} now means not given}{
      \code{NULL}, not \code{NA}, is the one spelling of an absent value. Where 0.9-34's own default was \code{NA}, an explicit \code{NA} is still read as \code{NULL} after a once-per-session warning: \code{sigest} on \code{\link{bart}}; \code{seed} on \code{bart}, \code{\link{bartBT}} and \code{\link{xbart}}; \code{n.samples} on \code{\link{dbartsControl}}; \code{updateState} on the \code{\linkS4class{dbartsSampler}} methods and \code{\link{updatePredictorPerObservationJointly}}; and \code{numBurnIn} and \code{numSamples} on \code{$run}. An argument new in 1.0-0 refuses an \code{NA}, naming \code{NULL}: \code{sigest} and \code{seed} on \code{\link{dbarts}}, \code{sigest} on \code{xbart}, \code{seed} on \code{dbartsControl}, and both on \code{dbartsSpec}, which is new, and \code{seed} on \code{\link{dbartsValidateComposition}}.
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
    \item{The front-door defaults message}{
      A \code{\link{bart}} call that names no \code{bartBT} argument and whose first argument is not a formula - a 0.9-x call such as \code{bart(x.train, y.train, x.test)} - binds the same way under both doors, so nothing forwards it and it runs under \code{bart}'s own defaults. The first such call in a session prints a message of class \code{dbartsFrontDoorMessage} naming the change and \code{\link{bartBT}}; \code{suppressMessages} silences it, and a call made from package code prints nothing.
    }
  }

  \code{control} on \code{\link{xbart}} is no longer a tombstone: both \code{bart} and \code{xbart} now take a \code{control = \link{dbartsControl}()} formal, reaching every setting neither door spells flatly, so the earlier retirement (which refused the argument outright) is reversed.

  Every entry above, \code{bart2} and \code{rbart_vi}, is one row of \code{dbarts:::dbartsTombstones}, the registry this page documents by hand.
}
\value{
  \code{bart2} returns whatever \code{\link{bart}} returns.
}
\seealso{
  \code{\link{bart}}, \code{\link{bartBT}}
}
\keyword{internal}
