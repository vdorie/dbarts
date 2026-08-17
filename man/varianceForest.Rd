\name{varianceForest}
\alias{varianceForest}
\alias{print.dbartsVarianceForest}
\alias{format.dbartsVarianceForest}
\title{Variance Forest Specification for Heteroscedastic BART}
\description{
  Build the specification of the heteroscedastic variance forest (HBART;
  Pratola, Chipman, George, and McCulloch 2020). Pass the result as the
  \code{variance} argument of \code{\link{dbarts}}, \code{\link{dbartsSpec}},
  or \code{\link{bart2}} - the same argument that already accepts the plain
  selector (\code{NULL}/\code{FALSE} for none, \code{TRUE}/a one-sided
  formula/a character or integer column selector for the predictors driving
  the variance). This constructor is the other accepted type of that one
  argument, adding the variance forest's own tree count and structure prior
  alongside the column selection.
}
\usage{
varianceForest(vars = NULL, n.trees = NULL, base = NULL, power = NULL)
}
\arguments{
  \item{vars}{
    Which predictors drive the variance surface: a one-sided formula
    (\code{~ x1 + x2}), or a character/integer column selector, resolved the
    same way the plain \code{variance} selector is (a bare factor term name
    expands to its indicator columns). \code{NULL} (the default) reads
    \emph{every} predictor - the object's own \dQuote{unrestricted} reading,
    distinct from \code{variance = NULL}, which declares no variance forest
    at all.
  }
  \item{n.trees}{
    The number of variance trees. \code{NULL} (the default) takes the
    variance forest's own default of \code{40}, the same value \code{variance
    = TRUE} and the other plain-selector forms use.
  }
  \item{base, power}{
    The variance forest's tree-structure prior. \code{NULL} (the default)
    inherits the mean forest's \code{tree.prior} \code{base}/\code{power}.
  }
}
\details{
  A second ensemble models the residual variance surface \eqn{s^2(x)} as a
  product of scaled-inverse-chi-squared (multiplicative) leaves, so
  \eqn{y_i = f(x_i) + s(x_i)\epsilon_i}, coupled to the mean forest through
  the per-observation precision. Gaussian responses and constant leaves only;
  monotone constraints are not supported. Every value is validated when a
  sampler is built.
}
\value{
  A \code{dbartsVarianceForest} specification object, resolved when a
  sampler is built.
}
\references{
Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch, R.E. (2020)
   Heteroscedastic BART via multiplicative regression trees.
        \emph{Journal of Computational and Graphical Statistics}, \bold{29}(2), 405--417.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{dbartsSpec}}, \code{\link{bart2}}
}
\examples{
set.seed(0)
n <- 200L
x <- matrix(runif(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
s <- ifelse(x[, 1] > 0.5, 0.3, 1.2)
y <- 2 * x[, 1] + s * rnorm(n)

## every predictor drives the variance, 20 variance trees
fit <- bart2(x, y,
             variance = varianceForest(n.trees = 20L),
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## restricted to x1, with its own tree-structure prior
fit.restricted <- bart2(x, y,
                        variance = varianceForest(vars = ~x1, n.trees = 20L,
                                                  base = 0.9, power = 1.5),
                        n.trees = 25L, n.samples = 20L, n.burn = 20L,
                        n.chains = 1L, verbose = FALSE)

varianceForest(vars = ~x1, n.trees = 20L, base = 0.9, power = 1.5)
}
\keyword{models}
