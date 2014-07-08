\name{makeModelMatrixFromDataFrame}
\title{Make Model Matrix From Data Frame}
\alias{makeModelMatrixFromDataFrame}
\alias{makeind}
\description{
  Converts a data frame with numeric and factor contents into a matrix,
  suitable for use with \code{\link{bart}}. Unlike in linear regression,
  factors containing more than two levels result in dummy variables being
  created for each level.
}
\usage{
  makeModelMatrixFromDataFrame(x)
  makeind(x, all = TRUE)
}
\arguments{
  \item{x}{Data frame of explanatory variables.}

  \item{all}{Not currently implemented.
   %If all=TRUE, a factor with p levels will be replaced by all p dummies.
   %If all=FALSE, the pth dummy is dropped.
   }
}
\details{
  Note that if you have train and test data frames, it may be best
  to \code{\link{rbind}} the two together, apply \code{makeModelMatrixFromDataFrame}
  to the result, and then pull them back apart.
}
\value{
   A matrix with columns corresponding to the elements of the data frame.
}
\author{
  Vincent Dorie: \email{vjd4@nyu.edu}.
}
\examples{
x1 <- 1:10
x2 <- as.factor(c(rep(1, 5), rep(2, 5)))
x3 <- as.factor(rep(c(1, 2, 4), c(3, 3, 4)))
levels(x3) <- c('alice', 'bob', 'charlie')
x <- data.frame(x1, x2, x3)

ignored <- makeModelMatrixFromDataFrame(x)
}
\keyword{factor}
