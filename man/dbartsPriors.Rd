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
      responses, \code{chi(1.5, 2)} for binary ones. The continuous
      default follows Chipman, George, and McCulloch's argument that with
      node standard deviation \code{sigma_mu = 0.5 / (k * sqrt(m))} for
      \code{m} trees, \code{k} prior standard deviations of \eqn{f(x)}
      span the whole response range regardless of \code{m}; see
      \code{\link{dbarts}}'s Details for the response-scaling caveat this
      relies on.
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
      priors; \code{columns} resolves as for \code{linear}, and \code{k}
      may be fixed or given the \code{chi} hyperprior, sampled from the
      drawn functions' standardized magnitudes.
      \code{lengthscale} fixes the kernel lengthscales on the
      standardized covariate scale, one value per column or one recycled
      scalar; \code{NULL} uses the median pairwise-distance heuristic per
      column. Leaves holding more than \code{max.leaf.size} observations
      fall back to constant fits, confining the cubic kernel cost.
      Function-valued fits ride prediction only: \code{getTrees} reports
      \code{NA} leaf values, and \code{keepTrees} storage grows with the
      leaf sizes. \code{\link{xbart}} accepts the same specification
      through its own \code{node.prior} argument.

      \code{predict} at a training row re-krigs the jitter-free
      posterior mean from the cached kernel and drawn training values,
      while the fit recorded during sampling includes a small
      conditioning nugget on the kernel diagonal; the two differ by
      roughly 2e-3 to 3e-3 at training locations (never during MCMC,
      which reads the recorded fit directly - only \code{predict} and
      test-data prediction re-krig).
    }
    \item{\code{chi(degreesOfFreedom = 1.5, scale = 2)}}{
      Chi hyperprior over \code{k}, sampled along with the rest of the
      model: \code{k} is given a chi distribution with the stated degrees
      of freedom and scale. The default, \code{chi(1.5, 2)}, centers the
      sampled \code{k} near the field-standard fixed value of 2 (prior
      median 1.9) while letting it adapt to the data. \code{scale = Inf}
      remains accepted for the old improper prior; the sampled \code{k} is
      capped at 1e6, a backstop that only engages at that improper scale.
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
  normal(chi(1.5)))} works regardless of what packages are attached: those
  arguments are evaluated with this vocabulary layered over the calling
  environment, along with \code{num.vars}, the number of predictor columns.
  A prior object built ahead of time with \code{dbartsPriors$...} can be
  passed to the same arguments.
}
\examples{
prior <- dbartsPriors$normal(dbartsPriors$chi(1.5))

x <- matrix(runif(200), ncol = 2)
y <- x[, 1] + rnorm(100, 0, 0.5)
sampler <- dbarts(y ~ x, node.prior = prior,
                  tree.prior = dbartsPriors$cgm(power = 1.5))

## DART: a Dirichlet prior over split-variable probabilities, useful
## when only a handful of predictors out of many are relevant
set.seed(0)
n <- 60L
x.dart <- matrix(runif(n * 8), n)
y.dart <- 4 * x.dart[, 1] + rnorm(n, 0, 0.2)
fit.dart <- dbarts(y.dart ~ x.dart, tree.prior = dart(),
                   control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                            n.threads = 1L))
samples.dart <- fit.dart$run(20L, 20L)

## linear leaves: each leaf fits an intercept plus a slope in x1,
## useful for a function that is smooth and roughly linear in a subset
## of the predictors
set.seed(1)
x1 <- runif(n)
x2 <- runif(n)
y.lin <- 3 * x1 + rnorm(n, 0, 0.2)
df.lin <- data.frame(x1, x2, y.lin)
fit.lin <- dbarts(y.lin ~ x1 + x2, df.lin, node.prior = linear("x1"),
                  control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                           n.threads = 1L))
samples.lin <- fit.lin$run(20L, 20L)

## gp leaves: each leaf fits a smooth Gaussian-process function of x1,
## useful for a smoothly-varying, nonlinear function
set.seed(2)
y.gp <- sin(2 * pi * x1) + rnorm(n, 0, 0.2)
df.gp <- data.frame(x1, x2, y.gp)
fit.gp <- dbarts(y.gp ~ x1 + x2, df.gp,
                 node.prior = gp("x1", max.leaf.size = 30L),
                 control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                          n.threads = 1L))
samples.gp <- fit.gp$run(20L, 20L)
}
\references{
Chipman, H., George, E., and McCulloch, R. (2010)
   BART: Bayesian additive regression trees.
        \emph{The Annals of Applied Statistics}, \bold{4}(1), 266--298.
        \doi{10.1214/09-AOAS285}.

Linero, A.R. (2018)
   Bayesian regression trees for high-dimensional prediction and variable selection.
        \emph{Journal of the American Statistical Association}, \bold{113}(522), 626--636.

Gramacy, R.B. and Lee, H.K.H. (2008)
   Bayesian Treed Gaussian Process Models With an Application to Computer Modeling.
        \emph{Journal of the American Statistical Association}, \bold{103}(483), 1119--1130.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{bart2}}
}
\keyword{datasets}
