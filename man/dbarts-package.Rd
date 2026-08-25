\name{dbarts-package}
\alias{dbarts-package}
\docType{package}
\title{Discrete Bayesian Additive Regression Trees Sampler}
\description{
  Fits Bayesian additive regression trees (BART; Chipman, George, and McCulloch 2010), a Bayesian \dQuote{sum-of-trees} model in which each tree is held to be a weak learner by its prior. The package offers two ways in, described below, and serves as a drop-in replacement for package \pkg{BayesTree}.

  What distinguishes \pkg{dbarts} from other BART implementations is that its sampler is \emph{mutable}: predictors, response, offset, and weights can be swapped \emph{between} MCMC draws, so that BART can be a conditional model inside a larger Gibbs or Metropolis-Hastings sampler rather than only a standalone fit.
}
\section{Fitting a model}{
  For an ordinary fit, use one of the fitting functions. Each runs the sampler to completion and returns posterior draws.

  \itemize{
    \item \code{\link{bart}} - the \pkg{BayesTree}-compatible interface, kept at its historical defaults (200 trees, one chain, factors expanded to indicator columns, binary responses probit, missing data rejected).
    \item \code{\link{bart2}} - the current interface, taking a formula or matrices and reaching the full feature set, including every response \code{family}. This is the one to reach for first.
    \item \code{\link{rbart_vi}} - adds a random intercept per group to the model of \code{bart2}.
    \item \code{\link{xbart}} - crossvalidates over \code{k}, \code{power}, \code{base}, and the tree count.
    \item \code{\link{pdbart}} - partial dependence plots for one or two variables.
  }

  Fits support \code{\link[=predict.bart]{predict}}, \code{\link[=fitted.bart]{fitted}}, \code{\link[=residuals.bart]{residuals}}, \code{\link{extract}}, \code{\link[=plot.bart]{plot}}, and \code{\link[=summary.bart]{summary}}; see \code{\link{bart}}'s \sQuote{Generics} section for what each returns and on which scale. \code{\link{survivalProbabilities}} turns a survival fit into survival-curve draws.
}
\section{Driving the sampler directly}{
  \code{\link{dbarts}} builds a \code{\link[=dbartsSampler-class]{dbartsSampler}} without running it. The sampler is a \link[=ReferenceClasses]{reference class} whose methods - \code{run}, \code{setResponse}, \code{setPredictor}, \code{setOffset}, \code{setWeights}, and the rest - map onto the underlying C++ engine, so an outer sampler can alternate its own draws with BART's. \code{\link{dbartsData}} and \code{\link{dbartsControl}} build its data and control objects, \code{\link{dbartsPriors}} its priors, and \code{\link{samplePriorPredictive}} draws from the prior for calibration before any fitting. \code{\link{dbartsSpec}} resolves the same specification without constructing a sampler, for packages that embed dbarts and supply their own design matrix.
}
\section{Model features}{
  Reached through \code{\link{bart2}} and \code{\link{dbarts}} unless noted.

  \itemize{
    \item Response families (\code{family}): \code{"gaussian"}, \code{"probit"}, \code{"logistic"}, accelerated failure time (\code{"aft"}), discrete-time survival hazard (\code{"hazard"}, \code{"hazard.logistic"}), multinomial (\code{"multinomial"}, \code{bart2} only), ordered categorical (\code{"ordinal"}), negative-binomial counts (\code{"nbinom"}), and semicontinuous two-part (\code{"hurdle.lognormal"}, \code{bart2} only). \code{"auto"}, the default, resolves the family from the response.
    \item Outlier-robust Student-t errors (\code{resid.dist = student(...)}) and a heteroscedastic variance forest (\code{variance}).
    \item Structural constraints: monotonicity (\code{monotone}), interaction limits (\code{\link{interactions}}), and block-additivity (\code{\link{blocks}}).
    \item Leaf models: constant, linear, or Gaussian-process (\code{node.prior}; see \code{\link{dbartsPriors}}).
    \item Variable selection by the DART prior (\code{dart}).
    \item Predictors: missing values incorporated in place (\code{missing}), unordered factors split on level subsets (\code{factors}), and sparse or mixed dense/sparse input (\code{Matrix::dgCMatrix}, \code{\link{sparseFactor}}).
    \item Warm starts from a previous fit (\code{warm.start}) or by XBART-style grow-from-root (\code{n.grow.sweeps}).
  }
}
\section{Reproducibility}{
  Every chain runs its own random number generator, so results never depend on the thread count and sampling does not advance R's stream. Calling \code{\link{set.seed}} beforehand suffices; passing \code{seed} instead makes a fit reproducible without touching R's stream at all. See \code{\link{bart}}'s \sQuote{Reproducibility} section.
}
\section{Using dbarts from another package}{
  Calling the R interface needs no special linkage. A package that drives the sampler from C should declare \code{LinkingTo: dbarts} and use the flat C API declared in the installed header \file{dbarts/dbarts.h}, whose entry points are reached through \code{R_GetCCallable} and versioned by the two-component handshake \code{DBARTS_C_API_MAJOR}/\code{DBARTS_C_API_MINOR}, with \code{dbarts_apiHash()} as an opt-in exact-ABI check, which moves on additive releases as well. \code{\link{dbartsSpec}} produces the specification objects that API's sampler constructor takes. The C++ interface of releases before 1.0-0, \file{dbarts/R_C_interface.hpp}, has been removed.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2010)
   BART: Bayesian additive regression trees.
        \emph{The Annals of Applied Statistics}, \bold{4}(1), 266--298.
        \doi{10.1214/09-AOAS285}.
}
\author{
  Vincent Dorie \email{vdorie@gmail.com}, with Hugh Chipman and Robert McCulloch. See \code{citation("dbarts")} for how to cite the package.
}
\seealso{
  \code{\link{bart2}} to fit a model, \code{\link{dbarts}} to build a sampler, and the package vignettes: \code{vignette("gibbs_sampler_mixture_model", package = "dbarts")} and \code{vignette("working_with_saved_trees", package = "dbarts")}.
}
\keyword{package}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
