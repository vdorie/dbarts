\name{rhat}
\alias{rhat}
\title{R-hat of One Quantity's Draws}
\description{
  The rank-normalized split R-hat (Vehtari, Gelman, Simpson, Carpenter and
  Buerkner, 2021) of the draws of one quantity from several chains: the larger
  of the bulk and the folded (tail) values. It is the statistic
  \code{\link{summary.bart}} reports, and from four iterations on it returns
  the value \code{posterior::rhat} does for the inputs that function takes
  (see Details for the inputs where the two differ).
}
\usage{
rhat(x)
}
\arguments{
  \item{x}{
    A numeric (or logical) matrix of the draws of one quantity, with a row for
    each iteration and a column for each chain. A vector is taken as one chain.
  }
}
\details{
  Each chain is cut in two, the middle draw of an odd length dropped, and the
  draws are rank-normalized as one pool; the bulk value compares the
  half-chains' means and variances, and the tail value does the same on the
  draws folded about their median. R-hat compares each chain's spread with the
  pooled spread, so a value above 1 says the chains disagree.

  The result is \code{NA}, as for \code{posterior::rhat}, when \code{x} has an
  \code{NA}, \code{NaN} or constant draws, or has a single iteration. A
  non-finite draw other than \code{NA} and \code{NaN} is ranked, as
  \pkg{posterior} ranks it. An argument that is not numeric or logical, is an
  array of more than two dimensions, or has no iteration or no chain is an
  error. \pkg{posterior} refuses a data frame and a character vector too,
  but takes a complex vector (and returns an infinite value with a warning
  for a matrix with no chain); dbarts refuses those, as they are not draws.

  With two or three iterations and several chains, \pkg{posterior} (1.7.0)
  loses the dimension of its one-draw half-chains and treats each iteration as
  a chain; dbarts returns \code{NA}, as a half-chain of one draw has no
  variance. At four or more iterations the two agree to rounding.

  R-hat should be read on the estimand reported: at the defaults BART's chains
  routinely disagree on \code{sigma}, \code{k} and single predictions, while
  averaged estimands, such as a subgroup's mean prediction or a contrast
  between subgroups, mix far better. See \code{\link{summary.bart}}.
}
\section{Calling \code{rhat}}{
  Call it as \code{dbarts::rhat}. Other packages export a function named
  \code{rhat}: \pkg{posterior} and \pkg{bayesplot} (and \pkg{brms}
  re-exports \pkg{posterior}'s), and \pkg{rstan} exports \code{Rhat}.
  Whichever is attached last masks the others, and \pkg{bayesplot}'s takes a
  fitted Stan model, not a matrix. The \pkg{stan4bart} package attaches
  \pkg{dbarts}, so it brings this name in.
}
\value{
  A single number, \code{NA_real_} where the draws carry no information.
}
\references{
  Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Buerkner, P.-C.
  (2021) Rank-normalization, folding, and localization: an improved
  \eqn{\hat{R}}{Rhat} for assessing convergence of MCMC. \emph{Bayesian
  Analysis} 16(2), 667-718.
}
\seealso{
\code{\link{summary.bart}}, \code{\link{extract}}
}
\examples{
set.seed(1)
draws <- matrix(rnorm(400), 100L, 4L)
dbarts::rhat(draws)

# four chains that disagree
dbarts::rhat(draws + rep(0:3, each = 100L))
}
