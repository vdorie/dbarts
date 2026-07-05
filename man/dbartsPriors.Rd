\name{dbartsPriors}
\alias{dbartsPriors}
\docType{data}
\title{Prior Specification Constructors}
\description{
  A list of constructor functions building the prior specifications that the
  \code{tree.prior}, \code{node.prior}, and \code{resid.prior} arguments of
  \code{\link{dbarts}} and the fitting functions accept. Bundling them keeps
  generic names like \code{normal} out of the search path, where another
  package could mask them or be masked depending on attach order.
}
\format{
  A list of functions:
  \describe{
    \item{\code{cgm(power = 2, base = 0.95, split.probs = NULL)}}{
      The Chipman, George, and McCulloch tree prior. \code{split.probs}
      weights the choice of split variable: \code{NULL} is uniform, a named
      vector assigns by column or term name (with an optional
      \code{".default"} element), an unnamed vector assigns by position.
      Named or positional probabilities are matched to the data when a
      sampler is built.
    }
    \item{\code{dart(power = 2, base = 0.95, a = 0.5, b = 1, rho = NULL,
                     alpha = 1, update.alpha = TRUE, update.delay = NULL)}}{
      The CGM structure prior with DART (Linero 2018): a Dirichlet prior
      over the split-variable probabilities inducing variable selection.
      \code{alpha} is the concentration, optionally sampled
      (\code{update.alpha}) on a grid with a Beta(\code{a}, \code{b}) prior
      on \code{alpha / (alpha + rho)}; \code{rho} defaults to the number of
      predictors. Updates hold until \code{update.delay} iterations have
      passed (default: half the control's burn-in), so the forest is
      likelihood-informed when counts first enter the Dirichlet.
    }
    \item{\code{normal(k = NULL)}}{
      Normal prior on the node means. \code{k} scales the standard
      deviation and can be a positive scalar, a hyperprior built with
      \code{chi}, or \code{NULL} for the default: 2 for continuous
      responses, \code{chi(1.25, Inf)} for binary ones.
    }
    \item{\code{linear(columns, k = NULL)}}{
      Each leaf fits an intercept plus a linear term in the designated
      continuous predictor columns instead of a constant, so the forest
      models smoothly-varying coefficients. \code{columns} names model
      matrix columns (character) or indexes them (numeric) and is matched
      to the data when a sampler is built; factor columns cannot be
      designated. The covariates are standardized internally and every
      coefficient shares the \code{normal(k)} prior. Reported leaf values
      keep the intercept; \code{getTrees} adds one \code{beta.<column>}
      column per covariate. \code{\link{xbart}} accepts the same
      specification through its own \code{node.prior} argument.
    }
    \item{\code{gp(columns, k = NULL, lengthscale = NULL,
                   max.leaf.size = 256L)}}{
      Each leaf fits a smooth Gaussian-process function of the designated
      continuous predictor columns, drawn under a squared-exponential
      kernel whose prior scale ties to \code{k} exactly as the other leaf
      priors; \code{columns} resolves as for \code{linear}. \code{k} must
      be fixed (a positive scalar; binary responses therefore need it
      supplied). \code{lengthscale} fixes the kernel lengthscales on the
      standardized covariate scale, one value per column or one recycled
      scalar; \code{NULL} uses the median pairwise-distance heuristic per
      column. Leaves holding more than \code{max.leaf.size} observations
      fall back to constant fits, confining the cubic kernel cost.
      Function-valued fits ride prediction only: \code{getTrees} reports
      \code{NA} leaf values, and \code{keepTrees} storage grows with the
      leaf sizes. \code{\link{xbart}} accepts the same specification
      through its own \code{node.prior} argument.
    }
    \item{\code{chi(degreesOfFreedom = 1.25, scale = Inf)}}{
      Chi hyperprior over \code{k}, sampled along with the rest of the
      model.
    }
    \item{\code{chisq(df = 3, quant = 0.9)}}{
      Chi-squared prior on the residual variance, scaled so that an
      estimate of the residual standard deviation falls at the given
      quantile.
    }
    \item{\code{fixed(value = 1)}}{
      Fixed residual variance.
    }
  }
}
\details{
  Inside the prior arguments of the fitting functions the same constructors
  are available by bare name, so \code{dbarts(..., node.prior =
  normal(chi(1.25)))} works regardless of what packages are attached: those
  arguments are evaluated with this vocabulary layered over the calling
  environment, along with \code{num.vars}, the number of predictor columns.
  A prior object built ahead of time with \code{dbartsPriors$...} can be
  passed to the same arguments.
}
\examples{
prior <- dbartsPriors$normal(dbartsPriors$chi(1.25))

x <- matrix(runif(200), ncol = 2)
y <- x[, 1] + rnorm(100, 0, 0.5)
sampler <- dbarts(y ~ x, node.prior = prior,
                  tree.prior = dbartsPriors$cgm(power = 1.5))
}
\seealso{
  \code{\link{dbarts}}, \code{\link{bart2}}
}
\keyword{datasets}
