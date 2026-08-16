\name{dbarts-embedding}
\alias{dbarts-embedding}
\title{Embedding dbarts In A Larger Sampler}
\description{
  Reference for using BART as a conditional model inside a Gibbs or Metropolis-Hastings scheme the
  caller writes: which channel an outer block writes, which fit it reads back, whose prior is in
  force after a swap, and how a composition is checked. The runnable set is
  \code{vignette("dbarts-as-a-component", package = "dbarts")}, six named recipes with their
  output; this page is the rule each of them applies.
}
\details{
  A \code{\linkS4class{dbartsSampler}} is created once by \code{\link{dbarts}} and then run one
  sweep at a time with \code{run(0L, 1L)}, its inputs replaced between sweeps. The fitting
  functions \code{\link{bart}}, \code{\link{bart2}} and \code{\link{rbart_vi}} own their own loop
  and are not the surface for this; they return fit objects, and the R modelling accessors
  (\code{\link{fitted}}, \code{\link{predict}}, \code{\link{extract}}) apply to those.

  \subsection{The channels an outer block writes}{
    Each method replaces one input in place, on every chain, taking effect at the next sweep. All
    of them are documented in \code{\link[=dbartsSampler-class]{dbartsSampler}}.

    \tabular{lll}{
      \strong{method} \tab \strong{what the outer block owns} \tab \strong{typical use} \cr
      \code{$setResponse(y)} \tab the regression target \tab a partial residual, a working response, an imputed outcome \cr
      \code{$setOffset(o)} \tab an additive term BART does not fit \tab another block's linear predictor, a fixed exposure \cr
      \code{$setWeights(w)} \tab per-observation precision \tab a Polya-Gamma precision, a variance mixture, survey weights \cr
      \code{$setSigma(s)} \tab the residual standard deviation \tab a scale drawn elsewhere; requires \code{resid.prior = fixed()} \cr
      \code{$setPredictor(x, column)} \tab one predictor column \tab a latent covariate, an imputed predictor \cr
      \code{$setActiveRows(a)} \tab which rows enter the likelihood \tab a subset indicator drawn by another block \cr
      \code{$setCalibration(prior.scale = )} \tab the leaf prior in force \tab restating the prior after a response swap
    }

    A predictor moved per observation rather than wholesale takes
    \code{$setPredictor(x, column, forceUpdate = "partial")}, or
    \code{\link{updatePredictorPerObservationJointly}} when several samplers share the column.
    Both return a per-observation install mask, \code{FALSE} where the new value would empty a leaf
    and was rolled back. \strong{That mask is part of a Metropolis accept decision, not a
    diagnostic}: a declined observation holds its old value, so the likelihood read after the call
    is its old likelihood, and treating the resulting ratio as an acceptance leaves the host's copy
    of the latent disagreeing with what the samplers hold.

    Not every channel is open on every model. A refusal is part of the model rather than a
    limitation - \code{setSigma} on a family that fixes the residual scale by definition, a
    response swap that would re-anchor a calibration stated at creation - and each names its
    reason.
  }

  \subsection{Which fit to read back}{
    Five methods report a fitted quantity and no two report the same one; the table in
    \sQuote{Reading the fit}, \code{\link[=dbartsSampler-class]{dbartsSampler}}, is authoritative.
    The rule an outer block needs: it conditions on \eqn{f(x_i)}{f(x_i)}, so it reads
    \code{$getFitsWithoutOffset()} and adds back whatever offset it installed.
    \code{run()$train} reports the same location with the offset folded in, and
    \code{$getLatents()} reports the family's augmentation variable, which carries no offset and
    for \code{"logistic"}, \code{"nbinom"} and a Student-t residual distribution is a precision
    rather than a location at all. Differencing the two is short by exactly the installed offset,
    which reaches the outer block as a coefficient shrunk toward zero rather than as an error.
  }

  \subsection{Whose prior is in force}{
    The leaf prior is calibrated once, at creation, from the response the sampler was built on. A
    per-sweep \code{setResponse} deliberately does \emph{not} re-anchor it (\code{updateScale =
    FALSE}, the default), so that the sweeps are comparable - which also means a sampler built on a
    cold-start vector keeps that vector's range as its prior scale for the whole run. State the
    prior instead of inheriting it: \code{node.prior = normal(scale = )} at creation (see
    \code{\link{dbartsPriors}}), or \code{$setCalibration(prior.scale = )} afterward, with
    \code{$getCalibration()} reporting what is actually in force.

    A decomposition across \eqn{K}{K} samplers needs the same care in the other direction. Each
    carries a prior sized for the whole of \eqn{f}{f}, so their sum has \eqn{\sqrt{K}}{sqrt(K)}
    times the prior standard deviation one sampler would have; dividing each
    \code{prior.scale} by \eqn{\sqrt{K}}{sqrt(K)} restates them as shares of one budget.
  }

  \subsection{The seeding contract}{
    The value in force during a sweep is the one installed before it. Create the sampler with the
    outer block's initial value already in place, and write the channel before \code{run}, never
    after: a sweep run against a stale offset, weight or response has drawn from the wrong
    conditional and nothing later reports it. The mutators store no state unless
    \code{updateState = TRUE} is passed explicitly, so a per-sweep write costs nothing beyond the
    write itself.
  }

  \subsection{Augmentation written in R}{
    The per-observation augmentation a non-gaussian sampler runs internally is callable directly:
    \code{\link{dbartsDrawLatents}} draws the latent at a supplied fit and
    \code{\link[=dbartsDrawLatents]{dbartsWorkingResponse}} turns it into the quantity a host
    regresses on. A host that needs the latent for its own blocks runs the step itself over a
    gaussian sampler with \code{resid.prior = fixed(1)}, rather than reimplementing the family.
    Their \code{fit} argument is offset-free, as \code{$getFitsWithoutOffset()} reports it.
  }

  \subsection{Checking a composition}{
    Recovering a known truth on one simulated data set is weak evidence.
    \code{\link{dbartsValidateComposition}} runs simulation-based calibration over the composed
    one-sweep step and reports, per scalar functional, whether the kernel targets the posterior the
    stated prior and likelihood imply. It catches the defect class in which a posterior mean is
    reported where a draw belongs.
  }
}
\seealso{
  \code{\link[=dbartsSampler-class]{dbartsSampler}} for every method named here,
  \code{\link{dbarts}} for creating one, \code{\link{dbartsAugmentation}} for the augmentation
  helpers, \code{\link{dbartsValidateComposition}} for the calibration check,
  \code{\link{updatePredictorPerObservationJointly}} for the shared-predictor sweep, and
  \code{\link{samplePriorPredictive}} for inspecting a prior before fitting.
}
\keyword{documentation}
