\name{na.keepPredictors}
\alias{na.keepPredictors}
\title{Drop Rows With A Missing Response, Keep Missing Predictors}
\description{
  The \code{na.action} that \code{\link{bart}}, \code{\link{dbarts}} and
  \code{\link{dbartsData}} default to. It drops the rows whose
  \emph{response} is missing and keeps the rows whose predictors are.
}
\usage{
na.keepPredictors(object, \dots)
}
\arguments{
  \item{object}{
    A model frame, as \code{\link{model.frame}} hands its \code{na.action};
    the response is the column the frame's \code{terms} attribute names,
    and a frame with no response loses no rows.
  }
  \item{\dots}{
    Ignored, for compatibility with the other \code{na.action} functions.
  }
}
\details{
  BART routes a missing predictor value down a learned side of each split
  rule (\dQuote{Missingness Incorporated in Attributes}), so an incomplete
  row is data rather than a hole and dropping it discards information the
  model can use. A missing response, on the other hand, is nothing to fit.
  This function is that distinction: \code{\link{na.omit}} would drop both.
  \pkg{rpart}'s \code{na.rpart} is the behavioural precedent.

  The dropped rows are recorded exactly as \code{\link{na.exclude}} records
  them - a named integer vector of class \code{"exclude"} - so the training
  fits pad back to the caller's own row count through
  \code{\link{naresid}}: \code{fitted} and \code{residuals} on the fit are
  as long as the data, \code{NA} where a row was dropped.

  The base functions keep their usual meaning on these entry points:
  \code{\link{na.omit}} drops any row with a missing value anywhere and
  pads nothing, \code{\link{na.exclude}} drops the same rows and pads,
  \code{\link{na.fail}} refuses any missing value, and
  \code{\link{na.pass}} keeps every row - after which a missing response is
  an error, since there is nothing to fit it against.

  On the matrix (\code{x.train}, \code{y.train}) interface, where there is
  no model frame, the \code{na.action} is applied to a frame standing for
  the (response, predictors) pair, so every one of these functions means on
  that interface what it means on the formula one.
}
\value{
  \code{object} with the missing-response rows removed, carrying an
  \code{"na.action"} attribute of class \code{"exclude"} naming them; or
  \code{object} unchanged when no response value is missing.
}
\seealso{
  \code{\link{na.omit}}, \code{\link{na.exclude}}, \code{\link{naresid}},
  \code{\link{bart}}, \code{\link{dbartsData}}
}
\examples{
df <- data.frame(y = c(1, NA, 3, 4), x = c(1, 2, NA, 4))

# the response row goes, the predictor row stays
model.frame(y ~ x, df, na.action = na.keepPredictors)

# and the fit's training values pad back to the data's own length
fit <- bart(
  y ~ x, df,
  n.samples = 20L, n.burn = 20L, n.chains = 1L, verbose = FALSE
)
length(fitted(fit)) == nrow(df)
}
\keyword{NA}
