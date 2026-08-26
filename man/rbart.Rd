\name{rbart}
\alias{rbart_vi}
\alias{plot.rbart}
\alias{fitted.rbart}
\alias{extract.rbart}
\alias{predict.rbart}
\alias{residuals.rbart}
\title{Bayesian Additive Regression Trees with Random Effects}
\description{
  Fits a varying intercept/random effect BART model.
}
\usage{
rbart_vi(
    formula, data, test, subset, weights, offset, offset.test = offset,
    group.by, group.by.test, prior = cauchy,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL, prior.scale = NA_real_,
    power = 2.0, base = 0.95,
    split.probs = NULL,
    dart = FALSE,
    n.trees = 75L,
    n.samples = 1500L, n.burn = 1500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = TRUE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 5L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE,
    keepTrees = TRUE, keepCall = TRUE,
    seed = NA_integer_,
    keepSampler = keepTrees,
    keepTestFits = TRUE,
    callback = NULL,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "aft"),
    missing = c("incorporate", "error"),
    storage = c("double", "single"), updateState = TRUE)

\S3method{plot}{rbart}(
    x, plquants = c(0.05, 0.95), cols = c('blue', 'black'), \dots)

\S3method{fitted}{rbart}(
    object,
    type = c("ev", "ppd", "bart", "ranef"),
    sample = c("train", "test"),
    ci.level = NULL,
    \dots)

\S3method{extract}{rbart}(
    object,
    type = c("ev", "ppd", "bart", "loglik", "ranef", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE,
    \dots)

\S3method{predict}{rbart}(
    object, newdata,
    type = c("ev", "ppd", "bart", "ranef"),
    offset = NULL, weights = NULL,
    combineChains = TRUE,
    ci.level = NULL,
    n.threads,
    \dots, group.by)

\S3method{residuals}{rbart}(object, type = "ev", \dots)
}
\arguments{
  \item{group.by}{
    Grouping factor. Can be an integer vector/factor, or a reference to such in \code{data}. For \code{predict} and \code{\link{survivalProbabilities}}, supplied by name only - it follows \code{\dots} in the signature, so it is never matched positionally; a missing one is refused, naming itself.
  }
  \item{group.by.test}{
    Grouping factor for test data, of the same type as \code{group.by}. Can be missing.
  }
  \item{prior}{
    A function or symbolic reference to built-in priors. Determines the prior over the standard deviation of the random effects. Supplied functions take two arguments, \code{x} - the standard deviation, and \code{rel.scale} - the standard deviation of the response variable before random effects are fit. Built in priors are \code{cauchy} with a scale of 2.5 times the relative scale and \code{gamma} with a shape of 2.5 and scale of 2.5 times the relative scale. With a built-in prior and no \code{callback}, the whole Gibbs sampler runs inside the engine on one multi-chain sampler; a custom prior function runs the random effect updates in R instead.
  }
  \item{n.thin}{
    The number of tree jumps taken for every stored sample, but also the number of samples from the posterior of the standard deviation of the random effects before one is kept.
  }
  \item{keepTestFits}{
    Logical where, if false, test fits are obtained while running but not returned. Useful with \code{callback}.
  }
  \item{callback}{
    Optional function of \code{trainFits}, \code{testFits}, \code{ranef}, \code{sigma}, and \code{tau}. Called after every post-burn-in iteration and the results of which are collected and stored in the final object.
  }
  \item{formula, data, test, subset, weights, offset, offset.test, sigest, sigdf, sigquant, power, base, split.probs, dart, n.trees, n.samples, n.burn, n.chains, n.threads, combineChains, n.cuts, useQuantiles, keepTrainingFits, printEvery, printCutoffs, verbose, keepTrees, keepCall, seed, keepSampler, factors, missing, storage, updateState}{
    Same as in \code{\link{bart2}}, with one default difference: \code{keepTrees} defaults to \code{TRUE} here (\code{bart2}'s default is \code{FALSE}). With \code{dart}, split probability samples appear as \code{varprobs} on the fit. Unlike \code{bart2}, \code{rbart_vi} offers a reduced \code{family} set (see below) with no logistic option, and does not accept sparse (\code{Matrix::dgCMatrix}) predictors.
  }
  \item{k}{
    As in \code{\link{bart2}}: \code{NULL} (the default) uses the value 2 for continuous responses and the \code{chi(1.5, 2)} hyperprior for binary ones; a supplied fixed value or hyperprior overrides that default.
  }
  \item{prior.scale}{
    As in \code{\link{bart2}}: names the forest's leaf calibration in response units (the prior standard deviation of \eqn{f} at \code{k = 1}) instead of inheriting it from the response range. The group intercepts' own prior is a separate object and is unaffected. \code{NA} (the default) changes nothing.
  }
  \item{family}{
    One of \code{"auto"}, \code{"gaussian"}, or \code{"aft"}. With \code{"auto"} (the default) the response family is resolved as gaussian or probit from the response, exactly as before this argument existed; a two-level factor, logical, or two-level character response is also detected and fit as probit (reporting the choice in a one-line message), while a factor response with three or more levels is an error, since \code{rbart_vi}'s random-effects model does not fit the multinomial family. \code{"gaussian"} forces the continuous model. \code{"aft"} fits a grouped accelerated failure time (AFT) log-normal survival model, \eqn{\log T_i = f(x_i) + \alpha_{g[i]} + \sigma \epsilon_i} with \eqn{\epsilon_i \sim N(0, 1)}, adding random intercepts to the model of \code{\link{bart2}}'s \code{family = "aft"}. The survival response enters through the formula's left-hand side as a \code{survival::Surv} object (right-censoring only) or a two-column \code{cbind(time, status)} of the event/censoring time and the 0/1 event indicator; a \code{Surv} response selects \code{"aft"} on its own from \code{"auto"}, while a bare two-column response needs \code{family = "aft"} explicitly. As in \code{\link{bart2}}, \code{predict}, \code{extract}, and \code{fitted} then return the linear predictor \eqn{E[\log T \mid x, g]} on the LOG-TIME scale (never the time scale), and \code{\link{survivalProbabilities}} gives survival-probability draws that include the drawn intercepts. Survival fits do not support \code{weights} or \code{subset} in this version, and enter only through the formula interface. This vocabulary is narrower than \code{\link{bart2}}'s by design - a grouped random-effects model has no multinomial/ordinal/nbinom/hurdle counterpart in this version; the wider family set lives on \code{bart2}.
  }
  \item{object}{
    A fitted \code{rbart} model.
  }
  \item{newdata}{
    Same as \code{test}, but named to match \code{\link{predict}} generic.
  }
  \item{type}{
    One of \code{"ev"}, \code{"ppd"}, \code{"bart"}, \code{"loglik"}, \code{"ranef"}, or \code{"trees"} for the posterior of the expected value, posterior predictive distribution, non-parametric/BART component, training log-likelihood (\code{extract} only), random effect, or saved trees respectively. The expected value is the sum of the BART component and the random effects, while the posterior predictive distribution is a response sampled with that mean. \code{"loglik"} evaluates the log-likelihood of each training observation at each posterior draw, conditioning on the drawn random intercepts as well as \eqn{\sigma}{sigma} (gaussian) or the fitted probability (binary); when chains are combined the result is a samples-by-observations matrix directly consumable by WAIC/PSIS-LOO implementations such as those in the \pkg{loo} package; the chains-first convention is kept, so a per-chain array (\code{combineChains = FALSE}, dimension chains-by-samples-by-observations) is reordered to the draws-by-chains-by-observations that \code{loo::relative_eff} expects with \code{aperm(x, c(2, 1, 3))}. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"ev"} and \code{"link"} can be used as a synonym for \code{"bart"}. For additional details on tree extraction, see the corresponding subsection in \code{\link{bart}}.
  }
  \item{sample}{
    One of \code{"train"} or \code{"test"}, referring to the training or test samples respectively.
  }
  \item{ci.level}{
    For \code{fitted} and \code{predict}, an optional single number in \eqn{(0, 1)}. As in \code{\link{bart}}: \code{NULL} (the default) returns the posterior mean, while a level returns a matrix of \code{est}, \code{ci.lower}, and \code{ci.upper}, with the interval kind following \code{type}.
  }
  \item{x, plquants, cols}{
    Same as in \code{\link{plot.bart}}.
  }
  \item{\dots}{
    Unused; present on \code{plot}/\code{fitted}/\code{extract}/\code{predict}/\code{residuals} only for S3 generic compatibility.
  }
}
\details{
  Fits a BART model with additive random intercepts, one for each factor level of \code{group.by}. For continuous responses:

\itemize{
  \item \eqn{y_i \sim N(f(x_i) + \alpha_{g[i]}, \sigma^2)}{y_i ~ N(f(x_i) + \alpha_{g[i]}, \sigma^2)}
  \item \eqn{\alpha_j \sim N(0, \tau^2)}{\alpha_j ~ N(0, \tau^2)}.
}

For binary outcomes the response model is changed to \eqn{P(Y_i = 1) = \Phi(f(x_i) + \alpha_{g[i]})}. \eqn{i} indexes observations, \eqn{g[i]} is the group index of observation \eqn{i}, \eqn{f(x)} and \eqn{\sigma_y} come from a BART model, and \eqn{\alpha_j} are the independent and identically distributed random intercepts. Draws from the posterior of \eqn{\tau} use an exact conjugate update for the built-in \code{cauchy} prior (a Makalic-Schmidt inverse-gamma scale mixture, two fixed inverse-gamma draws per sweep) and a slice sampler otherwise; the in-engine slice sampler for the built-in \code{gamma} prior steps out with a width tied to the prior's scale, while the R implementation used with custom priors determines its width from the curvature of the posterior at its mode.

\subsection{Out Of Sample Groups}{
  Predicting random effects for groups not in the training sample is supported by sampling from their posterior predictive distribution, that is a draw is taken from \eqn{p(\alpha \mid y) = \int p(\alpha \mid \tau)p(\tau \mid y)d\tau}. For out-of-sample groups in the test data, these random effect draws can be kept with the saved object. For those supplied to \code{predict}, they cannot and may change for subsequent calls.
}

\subsection{Generics}{
  See the generics section of \code{\link{bart}}.
}
}
\value{
  An object of class \code{rbart}. Contains all of the same elements of an object of class \code{\link{bart}}, as well as the elements:

  \item{ranef}{
    Samples from the posterior of the random effects. A array/matrix of posterior samples. The \eqn{(k, l, j)} value is the \eqn{l}th draw of the posterior of the random effect for group \eqn{j} (i.e. \eqn{\alpha^*_j}{\alpha*_j}) corresponding to chain \eqn{k}. When \code{n.chains} is one or \code{combineChains} is \code{TRUE}, the result is collapsed down to a matrix: with \eqn{m = }\code{n.samples \%/\% n.thin} draws kept per chain, row \eqn{r} is chain \eqn{1 + (r - 1) \%/\% m}'s draw \eqn{1 + (r - 1) \%\% m} - chain 1's whole run first, then chain 2's (chain-major), not interleaved draw by draw. Every collapsed field the fit returns (\code{sigma}, \code{tau}, \code{yhat.train}, ...) shares this row order, so row \eqn{r} of one pairs with row \eqn{r} of another.
  }
  \item{ranef.mean}{
    Posterior mean of random effects, derived by taking mean across group index of samples.
  }
  \item{tau}{
    Matrix of posterior samples of \code{tau}, the standard deviation of the random effects. Dimensions are equal to the number of chains times the number of samples unless \code{n.chains} is one or \code{combineChains} is \code{TRUE}, in which case it collapses to a vector in \code{ranef}'s chain-major order (see above).
  }
  \item{\code{first.tau}}{
    Burn-in draws of \code{tau}.
  }
  \item{\code{callback}}{
    Optional results of \code{callback} function.
  }
}
\references{
Makalic, E. and Schmidt, D.F. (2016)
   A simple sampler for the horseshoe estimator.
        \emph{IEEE Signal Processing Letters}, \bold{23}(1), 179--182.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}
}
\seealso{
  \code{\link{bart}}, \code{\link{dbarts}}
}
\examples{
f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

n.g <- 10
g <- sample(n.g, length(y), replace = TRUE)
sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

y <- y + b[g]

df <- as.data.frame(x)
colnames(df) <- paste0("x_", seq_len(ncol(x)))
df$y <- y
df$g <- g

## low numbers to reduce run time
rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                     n.samples = 40L, n.burn = 10L, n.thin = 2L,
                     n.chains = 1L,
                     n.trees = 25L, n.threads = 1L)

## with dart = TRUE, split-variable probabilities are sampled under a
## Dirichlet prior, inducing variable selection
rbartFit.dart <- rbart_vi(y ~ . - g, df, group.by = g, dart = TRUE,
                          n.samples = 40L, n.burn = 10L, n.thin = 2L,
                          n.chains = 1L,
                          n.trees = 25L, n.threads = 1L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{randomeffects}

