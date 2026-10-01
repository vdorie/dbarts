\name{dbartsPriors}
\alias{dbartsPriors}
\docType{data}
\title{Prior Specification Constructors}
\description{
  A list of constructor functions building the prior specifications that the
  \code{tree.prior} and \code{leaf.prior} arguments of \code{\link{dbarts}}
  and the fitting functions accept. Bundling them keeps generic names like
  \code{normal} out of the search path, where another package could mask them
  or be masked depending on attach order.

  The residual priors \code{chisq} and \code{fixed} are what a family
  object's \code{sigma} setting takes - \code{family = gaussian(sigma =
  chisq(3, 0.9))} - which is the one place every entry point reaches the
  residual prior; this vocabulary resolves by bare name inside that argument
  too.
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
      predictors. \code{NA} is a missing value, not a spelling of the default,
      and is refused for \code{rho} and \code{update.delay} alike. Updates hold
      until \code{update.delay} iterations have
      passed (default: half the control's burn-in), so the forest is
      likelihood-informed when counts first enter the Dirichlet.
    }
    \item{\code{normal(k = NULL, sd = NULL)}}{
      Normal prior on the leaf values, named one of two ways, never both.
      Either names the spread of the whole forest - the sum of trees - not
      of one tree.

      \code{k} is relative to the anchor the data fixes: a positive scalar,
      a hyperprior on \code{k} built with \code{chi}, or \code{NULL} for
      the default, 2 for continuous responses and \code{chi(1.5, 2)} for
      binary ones. The continuous default follows Chipman, George, and
      McCulloch's argument that with leaf standard deviation
      \code{sigma_mu = 0.5 / (k * sqrt(m))} for \code{m} trees, \code{k}
      prior standard deviations of \eqn{f(x)} span the whole response
      range regardless of \code{m}; see \code{\link{dbarts}}'s Details for
      the response-scaling caveat this relies on. A string such as
      \code{"chi(1.5)"} is kept for 0.9-x.

      \code{sd} is the standard deviation of the normal prior on the leaf
      model's own parameter, for the forest total, on the scale the
      family's forest fits: a positive number, or a hyperprior on the sd
      itself built with \code{invchi}. It takes no string form. The scale
      and the anchor \code{k} is relative to, per family:
      \tabular{lll}{
        family \tab \code{sd} is stated in \tab anchor at \code{k = 1} \cr
        gaussian, student \tab response units \tab half the training response range \cr
        aft \tab log survival time \tab half the observed log-time range \cr
        probit, ordinal \tab probit latent \tab 3 \cr
        logistic, nbinom \tab log-odds latent \tab \eqn{\pi\sqrt{3}}{pi * sqrt(3)} \cr
        hazard \tab its link's latent \tab 3 or \eqn{\pi\sqrt{3}}{pi * sqrt(3)}
      }
      so a fixed \code{k} is the \code{sd} anchor / \code{k}. For the
      constant leaf \code{sd} is exactly the prior standard deviation of
      \eqn{f(x)} at every \eqn{x}; for the other leaf models it is the sd of
      their parameter, and the prior spread of \eqn{f(x)} it implies is a
      consequence stated under each. Under a \code{monotone} constraint
      (see \code{\link{dbarts}}) \code{sd} is the prior sd of a leaf value
      that no ordered neighbor bounds; a leaf a neighbor bounds is drawn
      from a normal \eqn{\sqrt{\pi / (\pi - 1)}}{sqrt(pi / (pi - 1))} times
      as wide, truncated to the ordering, which matches its marginal prior
      variance to \code{sd^2}, and \code{sd} is a lower bound on the prior
      sd of \eqn{f(x)} in the interior. A named \code{sd} is absolute: a
      sampler restates it after every channel that re-anchors the response
      transform (see \code{\link{dbartsSampler-class}}), where a \code{k}
      moves with the data. It is refused on multinomial and multi-forest
      models (whose spreads come from their calibration maps, the latter
      stated through \code{\link{forest}(sd = )}, a different quantity),
      and on a hurdle fit, whose two parts are on different scales.
      \code{NULL} leaves either unnamed; \code{NA} is a missing value and is
      refused, here and in \code{linear} and \code{gp}.
    }
    \item{\code{linear(columns, k = NULL, sd = NULL)}}{
      Each leaf fits an intercept plus a linear term in the designated
      continuous predictor columns instead of a constant, so the forest
      models smoothly-varying coefficients. \code{columns} names model
      matrix columns (character) or indexes them (numeric) and is matched
      to the data when a sampler is built; unordered factor columns cannot
      be designated. The covariates are standardized internally and every
      coefficient shares the \code{normal(k)} prior. Reported leaf values
      keep the intercept; \code{getTrees} adds one \code{beta.<column>}
      column per covariate. \code{\link{xbart}} accepts the same
      specification through its own \code{leaf.prior} argument.
      \code{k} and \code{sd} name the spread as they do for
      \code{normal}, \code{sd} per standardized covariate: every
      coefficient, the intercept included, has prior sd \code{sd} for the
      forest total, so the prior sd of \eqn{f(x)} is
      \code{sd * sqrt(1 + ||z||^2)} at standardized covariates \eqn{z}, and
      \code{sd} is its lower bound, attained at the covariate means.
    }
    \item{\code{gp(columns, k = NULL, lengthscale = NULL,
                   max.leaf.size = 256L, sd = NULL)}}{
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
      Fewer trees are needed than for constant leaves: 10 to 25 is a
      reasonable range, since each leaf already fits a smooth function, while
      past about 50 the per-tree shrinkage leaves each function too small to
      improve on constant fits and the cost, linear in the number of trees and
      cubic in the leaf size, keeps rising. The two settings interact: a tree
      typically holds only two or three leaves, so a leaf carries a sizable
      fraction of the sample and \code{max.leaf.size} must be comparable to
      that fraction for the kernel to be used at all; lowering \code{n.trees}
      grows deeper trees and is the effective way to bring leaves under the
      cap, whereas lowering \code{max.leaf.size} only sends more leaves to the
      constant fallback.
      Because that fallback is silent - the fit stays coherent, it is simply
      not a Gaussian process over the leaves that took it - every leaf
      evaluation is counted. The counts ride the fit as \code{gp.fallback}
      (a named pair, \code{evaluations} and \code{fallbacks}), and a run in
      which more than a quarter of evaluations fell back warns once (class \code{dbartsGPFallbackWarning}), giving
      the share: at the default cap on a design of any size that is the
      normal outcome, not an error, and the warning is the notice that the
      cap rather than the kernel is doing the modeling. The remedies are the
      two the paragraph above names, in the order it names them.
      Function-valued fits ride prediction only: \code{getTrees} reports
      \code{NA} leaf values, and \code{keepTrees} storage grows with the
      leaf sizes. \code{\link{xbart}} accepts the same specification
      through its own \code{leaf.prior} argument. \code{k} and \code{sd}
      name the spread as they do for \code{normal}, \code{sd} as the
      amplitude of the forest total's Gaussian process: the prior sd of
      \eqn{f(x)} is \code{sd} at a leaf's own rows and decays away from
      them, so \code{sd} is its upper bound.

      \code{predict} at a training row re-krigs the jitter-free
      posterior mean from the cached kernel and drawn training values,
      while the fit recorded during sampling includes a small
      conditioning nugget on the kernel diagonal; the two differ by
      roughly 2e-3 to 3e-3 at training locations (never during MCMC,
      which reads the recorded fit directly - only \code{predict} and
      test-data prediction re-krig).
    }
    \item{\code{chi(df = 1.5, scale = 2)}}{
      Chi hyperprior over \code{k}, sampled along with the rest of the
      model: \code{k} is given a chi distribution with \code{df} degrees
      of freedom and scale \code{scale}. The default, \code{chi(1.5, 2)},
      centers the sampled \code{k} near the field-standard fixed value of
      2 (prior median 1.9) while letting it adapt to the data. It is the
      same prior as \code{sd = invchi(df, anchor / scale)}. \code{scale =
      Inf} remains accepted, but the posterior it gives \code{k} is
      improper: with little signal, few trees or many degrees of freedom,
      \code{k} can drift to infinity, where every leaf is zero and the
      trees add nothing to the fit. The first argument was
      \code{degreesOfFreedom} in 0.9-x; that name is accepted, with a
      warning, until dbarts 1.1-0.
    }
    \item{\code{invchi(df = 1.5, scale)}}{
      Scaled inverse chi hyperprior over a leaf prior's \code{sd}, on the
      sd's own scale: \eqn{sd = scale / \chi_{df}}{sd = scale / chi_df},
      equivalently \eqn{sd^2}{sd^2} is scaled inverse chi-square with
      \code{df} degrees of freedom. \code{scale} has no default, since a
      spread on the family's scale has no data-free one; \code{scale = 0}
      is the improper \eqn{sd^{-(df + 1)}}{sd^-(df + 1)} limit. It is the
      same prior as \code{k = chi(df, anchor / scale)}, and
      \code{chi(df, Inf)} is \code{invchi(df, 0)}; the binary defaults are
      \code{invchi(1.5, 1.5)} on probit and
      \code{invchi(1.5, pi * sqrt(3) / 2)} on logistic. In Stan's terms it
      is \code{scaled_inv_chi_square(df, scale / sqrt(df))} on the
      variance, and an inverse gamma with shape \code{df / 2} and scale
      \code{scale^2 / 2}. No other law on the sd is offered: this is the
      one conjugate to the normal leaves, which the sampler draws exactly;
      a half-t or half-Cauchy would need a non-conjugate update. Refused
      under a monotone constraint, as a \code{chi} on \code{k} is.
    }
    \item{\code{chisq(df = 3, quant = 0.9)}}{
      Chi-squared prior on the residual variance. Its scale is set from
      \code{sigest}, the residual-scale estimate supplied (or derived) at
      creation: \code{quant} is the prior probability that the residual
      standard deviation is below that estimate, so \code{sigest = 1} gives
      the prior under which \eqn{\sigma} is below 1 with probability
      \code{quant}.
    }
    \item{\code{fixed(value = 1)}}{
      Fixed residual variance, in squared response units. It is the residual
      scale rather than a law over one, so it ignores \code{sigest} - which
      is refused beside it rather than accepted and overwritten - and
      suppresses the sampler's own \eqn{\sigma} draw.
    }
  }
}
\details{
  Inside the prior arguments of the fitting functions - and inside a family's
  \code{sigma} - the same constructors are available by bare name, so \code{dbarts(..., leaf.prior =
  normal(chi(1.5)))} works regardless of what packages are attached: those
  arguments are evaluated with this vocabulary layered over the calling
  environment (for an argument a wrapper forwards through its \code{\dots},
  the one it was written in), along with \code{num.vars}, the number of predictor columns.
  A prior object built ahead of time with \code{dbartsPriors$...} can be
  passed to the same arguments. A wrapper's named formal that forces a bare
  constructor call before passing it on reaches R's own could-not-find-function
  error instead, extended with a hint naming \code{dbartsPriors$cgm(...)} and
  its siblings.
}
\examples{
prior <- dbartsPriors$normal(dbartsPriors$chi(1.5))

x <- matrix(runif(200), ncol = 2)
y <- x[, 1] + rnorm(100, 0, 0.5)
sampler <- dbarts(y ~ x, leaf.prior = prior,
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
fit.lin <- dbarts(y.lin ~ x1 + x2, df.lin, leaf.prior = linear("x1"),
                  control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                           n.threads = 1L))
samples.lin <- fit.lin$run(20L, 20L)

## gp leaves: each leaf fits a smooth Gaussian-process function of x1,
## useful for a smoothly-varying, nonlinear function
set.seed(2)
y.gp <- sin(2 * pi * x1) + rnorm(n, 0, 0.2)
df.gp <- data.frame(x1, x2, y.gp)
fit.gp <- dbarts(y.gp ~ x1 + x2, df.gp,
                 leaf.prior = gp("x1", max.leaf.size = 52L),
                 control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                          n.threads = 1L))
samples.gp <- fit.gp$run(20L, 20L)
## the cap of 52 observations per leaf sends some leaves to the constant
## fallback (about a tenth of the evaluations here), counted in gp.fallback;
## a warning comes only when more than a quarter fall back
print(unlist(attr(samples.gp, "gp.fallback")))
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
  \code{\link{dbarts}}, \code{\link{bart}}
}
\keyword{datasets}
