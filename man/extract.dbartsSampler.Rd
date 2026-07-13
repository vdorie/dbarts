\name{extract.dbartsSampler}
\alias{extract.dbartsSampler}
\title{Extract the Predictor Matrix From a Sampler}
\description{
  Returns the predictor matrix as the sampler currently sees it - one numeric column per predictor, with factor columns as their integer codes. This is the canonical source for tree replay: it is the same matrix \code{getTrees} uses to count training observations through saved trees.
}
\usage{
\method{extract}{dbartsSampler}(object, type = "predictors", \dots)
}
\arguments{
  \item{object}{
    A \code{\linkS4class{dbartsSampler}}.
  }
  \item{type}{
    A character string; only \code{"predictors"} is supported.
  }
  \item{\dots}{
    Unused; included for compatibility with the \code{\link{extract}} generic.
  }
}
\details{
  Data frame input is stored columnar in the sampler's data object, so its \code{x} slot need not be a plain matrix. \code{extract} materializes the numeric coding regardless of the storage form: numeric and logical columns as doubles, and factor columns - ordered or not - as their integer codes, \code{0} to one less than the number of levels, in level order. This is the same coding a hand-built predictor matrix supplies, so the result can be fed back to \code{predict} or used to construct an equivalent sampler.

  After predictor mutation (\code{setPredictor} and relatives) the result reflects the sampler's current predictor values, not the originally ingested data. Each mutating call collects the accepted change into the sampler's data object directly; nothing is maintained per iteration, so the returned codes always reflect current values, materialized fresh on every call.
}
\value{
  A numeric matrix with one row per training observation and one column per predictor, carrying the predictors' column names.
}
\seealso{
  \code{\linkS4class{dbartsSampler}}, \code{\link{dbartsData}}, \code{\link{extract}}
}
\examples{
df <- data.frame(
  x = runif(50L),
  g = factor(sample(c("a", "b", "c"), 50L, replace = TRUE))
)
df$y <- df$x + 0.5 * (df$g == "b") + rnorm(50L, 0, 0.1)

sampler <- dbarts(y ~ x + g, df, control = dbartsControl(
  n.chains = 1L, n.threads = 1L, updateState = FALSE))
predictors <- extract(sampler, "predictors")
str(predictors)
stopifnot(all(predictors[, "g"] == as.integer(df$g) - 1L))
}
\keyword{models}
