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
\section{At Prediction}{
  Every \code{predict} method on a fit (\code{\link{predict.bart}} and
  \code{\link{bart}}'s four own classes) and
  \code{\link{survivalProbabilities}} take an \code{na.action} for
  \code{newdata}, defaulting to this function; \code{NULL} means the
  default, and \code{getOption("na.action")} is never consulted.

  A missing value in \code{newdata} is \emph{routable} when its column can
  hold missing values: it had them in training or, on a sampler, has been
  given them since (see \code{\link{dbartsSampler}}). Every split rule on
  that column then holds a side for them, and the row predicts like any
  other. On a column that has never held one no rule holds a side, so the
  value is \emph{unroutable}. Only the default and \code{\link{na.pass}} look at
  this:
  \itemize{
    \item \code{na.keepPredictors} predicts every routable row and refuses
    an unroutable one, naming its columns;
    \item \code{na.pass} predicts every routable row and returns \code{NA}
    for an unroutable one;
    \item \code{\link{na.omit}} drops every row with a missing value;
    \item \code{\link{na.exclude}} drops the same rows and pads them back
    as \code{NA}, so the result stays aligned with \code{newdata}, as the
    fit-time record does (\code{\link{predict.lm}} instead treats it as
    \code{na.omit} on new data);
    \item \code{\link{na.fail}} refuses any missing value, naming its
    columns.
  }
  Any other function is applied to a one-column frame that is \code{NA} on
  the rows with a missing value, on every interface, and the rows it keeps
  then meet the default's refusal. The result carries \code{newdata}'s row
  names and no \code{"na.action"} attribute. When no row survives, or
  \code{newdata} has none, the result has the usual draw dimensions and
  zero rows (\code{NA} rows under \code{na.exclude} and \code{na.pass}),
  and R's random number stream is left as it was.

  A per-row \code{offset}, \code{weights} or \code{bases} must match
  \code{newdata}'s rows and loses the rows \code{na.action} drops; a single
  value passes through unaffected. A missing value in a per-row
  \code{offset} or \code{weights} (an \code{n x K} category shift counts as
  one, on a multinomial fit) marks its row incomplete exactly as an
  unroutable predictor does, and every \code{na.action} above governs it the
  same way: the default and \code{na.fail} refuse it by name, \code{na.pass}
  returns \code{NA} for it, and \code{na.omit}/\code{na.exclude} drop it
  (\code{na.exclude} padding it back). A missing value in a single
  (length-one) \code{offset} or \code{weights} is always an error, whatever
  \code{na.action} is, since there is then no row for it to mark; the same
  goes for any missing value in \code{bases}. The sampler's own
  \code{predict}, \code{predictForests} and \code{getTrees} methods take no
  \code{na.action}, keep the default, and check no per-row
  \code{offset.test} for a missing value at all: an \code{NA} there reaches
  the engine unrouted.
}
\value{
  \code{object} with the missing-response rows removed, carrying an
  \code{"na.action"} attribute of class \code{"exclude"} naming them; or
  \code{object} unchanged when no response value is missing.
}
\seealso{
  \code{\link{na.omit}}, \code{\link{na.exclude}}, \code{\link{naresid}},
  \code{\link{bart}}, \code{\link{dbartsData}},
  \code{\link{survivalProbabilities}}
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
