\name{summary.bart}
\alias{summary.bart}
\alias{print.summary.bart}
\title{Convergence Diagnostics for BART Fits}
\description{
  Reports a per-variable posterior summary of the scalar parameters
  (\code{sigma} and the leaf prior's \code{k} or \code{leaf.prior.sd}) of a \code{bartBT}/\code{bart} fit, along
  with split-\eqn{\hat{R}}{Rhat} and bulk/tail effective sample size,
  computed by dbarts itself (no \pkg{posterior} dependency). A
  heteroscedastic fit (\code{\link{dbarts}}'s \code{variance}) has no
  scalar residual scale and reports \code{mean.s} in place of
  \code{sigma}; see \code{vars}.

  \code{bart}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}; see
  \code{\link{bart}}) summarize through this same method, each exposing
  the scalar posterior parameters that family carries rather than its
  per-observation channels - never \code{yhat.train} itself. A chain-
  separated (rather than pooled) view of these same draws is
  \code{\link{extract}}, whose \code{"sigma"}, \code{"k"}, \code{"leaf.prior.sd"}, and other
  scalar types return a vector or, with \code{combineChains = FALSE}, a
  chains-by-samples matrix. The table holds the parameters the fit
  sampled; one line under it names each requested parameter the fit held
  fixed, with its value, in the manner of \code{summary.glm}'s dispersion
  line. The summary of a fit of several forests adds the full block of what
  each forest resolved to, after the table (see the \code{sd} argument of
  \code{\link{forest}}), and the line under the table then leaves out each
  forest's \code{leaf.prior.sd}, which the block gives. Under the block, each
  forest with a basis has a line for the effect of one unit of its basis (for
  a factor, of its second level against its first): its \emph{average} over
  the training rows in the likelihood, and its \emph{sd across rows}, the
  standard deviation of the effect across those rows, each with its posterior median and
  the 5\% to 95\% interval \code{stats} uses. Both are computed for each draw
  and the median taken last, on the scale the forests add on (the response's
  units, or the probit or logit scale of a binary response), from
  each forest's own fits times the coefficient on its basis, which is what
  \code{predict} gives with the basis one unit higher less the basis as it
  stands; a fit that kept no per-forest fits (\code{keepFits = FALSE}) says so
  in one line instead. A change in probability is \code{predict}'s job (\code{predict} with
  the basis one unit up less \code{predict} as it stands). Where the weights
  are counts of identical rows (the count families' default, and gaussian with
  \code{weights.type = "frequency"}), the average and the sd across rows, as
  the forest line's column mean and sd, count each row at its count, as the
  data written out row by row would.
}
\usage{
\method{summary}{bart}(object, vars = c("sigma", "k", "leaf.prior.sd", "resid.df"), \dots)
\method{print}{summary.bart}(x, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{bart}, as returned by
    \code{\link{bartBT}} or \code{\link{bart}}; for the four own-class
    fits, a \code{bart} fit of the matching class.
  }
  \item{vars}{
    Character vector of fields to summarize. By default the leaf-scale
    field follows how the leaf prior was named: \code{k} for a fit that
    named it by \code{k}, \code{leaf.prior.sd} (the \code{k.scale} over \code{k};
    see \code{\link{bartBT}}, where on a fit with several forests it is each
    forest's \code{sd} in the response's units, printed in the block after the table) for one that named it by \code{sd}, and
    either can be requested by name, as can \code{first.k}, the burn-in
    draws. A requested parameter the fit held fixed (\code{sigma} under
    \code{fixed()}, a fixed \code{k} or \code{leaf.prior.sd}, a fixed
    \code{shape} or \code{resid.df}) is not tabulated but named on the
    line under the table; with the default \code{vars} that includes a
    fixed Student-t \code{resid.df}, and a sampled one is tabulated; other requested fields absent from
    \code{object} are silently dropped. \code{sigma}, \code{k} and
    \code{leaf.prior.sd} contribute one summarized variable each; any other field (\code{varcount}, \code{varprobs},
    \code{yhat.train}, \code{yhat.test}) contributes one variable per
    column, named \code{"field[column]"}.

    On a heteroscedastic fit the \code{sigma} token resolves to
    \code{mean.s}: the mean over the training observations of that
    draw's variance surface \eqn{s(x)}, one value per draw, since
    such a fit has no scalar \code{sigma}; \code{s.train} itself is
    reachable by name, contributing one variable per observation.

    For the four own-class fits, \code{vars} is scoped to that
    family's own vocabulary (see \code{\link{bart}}): a
    \code{"bartOrdinal"} fit's \code{"thresholds"} contributes
    \code{threshold[2]} through \code{threshold[K - 1]}, the first being pinned at 0 and named under the table;
    a \code{"bartNegbin"} fit's \code{"shape"} contributes the
    per-draw shape \eqn{r}; a \code{"bartMultinomial"} fit has no
    \code{vars} argument at all - its only scalar posterior parameter
    is the per-category mean predicted probability, reported as
    \code{prob[<level>]}; a \code{"bartHurdle"} fit applies \code{vars}
    to both components separately.
  }
  \item{\dots}{Unused.}
}
\value{
  An object of class \code{summary.bart} with elements \code{call},
  \code{stats} (a data frame with one row per requested sampled scalar, or
  \code{NULL} if none of \code{vars} are present on the fit), \code{vars},
  and \code{fixed}, a named list of the requested parameters the fit held
  fixed and their values (on a fit with several forests a vector with one
  entry for each forest, \code{k} and \code{leaf.prior.sd} named by the
  forests' labels), and, on a fit of
  several forests that kept its fits, \code{forests}, the lines of the full block
  that the print method gives after the table, and \code{effects}, where a
  forest has a basis, a data frame with a row for each effect (columns
  \code{forest}, \code{label}, \code{effect}, \code{average},
  \code{average.q5}, \code{average.q95}, \code{spread}, \code{spread.q5} and
  \code{spread.q95}, with the number of rows averaged over and the scale as
  attributes \code{"rows"} and \code{"scale"}) or the string \code{"dropped"}
  for a fit that kept no per-forest fits. \code{stats} has columns \code{variable}, \code{mean}, \code{median},
  \code{sd}, \code{mad}, \code{q5}, \code{q95}, \code{rhat},
  \code{ess_bulk}, and \code{ess_tail} - the same columns
  \code{posterior::summarise_draws} reports for a plain array, computed
  without that package. Its print method shows each R-hat as a value, with
  no bar and no note, and dbarts does not refuse to summarize a fit whose
  chains disagree.
}
\section{Reading R-hat}{
  R-hat compares each chain's spread with the spread of the chains pooled, so
  a value above 1 says the chains disagree; pooling them does not cure it. At
  the defaults BART's chains routinely disagree on \code{sigma}, on \code{k}
  and on a single row's prediction, so values above 1 there are expected, and
  no single bar says whether a fit is usable. Averaged estimands, such as the
  mean prediction over a subgroup or a contrast between two subgroups, mix far
  better, so read R-hat on the estimand you report. \code{\link{extract}} with
  \code{combineChains = FALSE} gives each chain's draws; the example computes
  R-hat of a contrast this way, with \code{posterior::rhat} when
  \pkg{posterior} is installed.
}
\seealso{
\code{\link{bart}}, \code{\link{bartBT}}, \code{\link{extract}}
}
\examples{
\donttest{
fit <- bart(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)
}

\donttest{
# R-hat of the estimand reported: the difference between the mean prediction
# for two subgroups, one value per draw and chain
set.seed(1)
n <- 200L
x <- matrix(runif(n * 5L), n, 5L)
y <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
  10 * x[, 4L] + 5 * x[, 5L] + rnorm(n)
fit <- bart(x, y, n.chains = 4L, n.samples = 500L, n.burn = 100L,
            n.trees = 50L, n.threads = 1L, verbose = FALSE)
# chains by samples by rows
yhat <- extract(fit, "ev", combineChains = FALSE)
group <- x[, 1L] > 0.5
contrast <- apply(yhat[, , group, drop = FALSE], 1:2, mean) -
  apply(yhat[, , !group, drop = FALSE], 1:2, mean)
if (requireNamespace("posterior", quietly = TRUE))
  posterior::rhat(t(contrast)) # rhat() takes samples by chains
}
}
\keyword{htest}
