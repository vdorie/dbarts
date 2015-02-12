\name{makeModelMatrixFromDataFrame}
\title{Make Model Matrix from Data Frame}
\alias{makeModelMatrixFromDataFrame}
\alias{makeind}
\description{
  Converts a data frame with numeric and factor contents into a matrix,
  suitable for use with \code{\link{bart}}. Unlike in linear regression,
  factors containing more than two levels result in dummy variables being
  created for each level.
}
\usage{
  makeModelMatrixFromDataFrame(x, drop = TRUE)
  makeind(x, all = TRUE)
}
\arguments{
  \item{x}{Data frame of explanatory variables.}
  \item{drop}{Logical or list controling whether or not columns that are constants or factor
    levels with no instances are omitted from the result. When a list, must be of
    length equal to \code{x}. Elements correspond to \code{x} according to:
    \itemize{
      \item vector - single logical
      \item matrix - vector of logicals, one per column
      \item factor - table of factor levels to be referenced; levels with counts of 0 are to be dropped
    }}
  \item{all}{Not currently implemented.
   }
}
\details{
  Note that if you have train and test data frames, it may be best
  to \code{\link{rbind}} the two together, apply \code{makeModelMatrixFromDataFrame}
  to the result, and then pull them back apart. Alternatively, save the drop
  attribute used in creating the train data and use it when creating a matrix
  from the test data. Example given below.
}
\value{
  A matrix with columns corresponding to the elements of the data frame. If \code{drop = TRUE}
  or is a list, the attribute \code{drop} on the result is set to the list used when creating
  the matrix.
}
\author{
  Vincent Dorie: \email{vjd4@nyu.edu}.
}
\examples{
iv <- 1:10
rv <- runif(10)
f <- factor(rep(seq.int(3), c(4L, 4L, 2L)),
            labels = c("alice", "bob", "charlie"))
df <- data.frame(iv, rv, f)

mm <- makeModelMatrixFromDataFrame(df)

## create test and train matrices with disjoint factor levels
train.df <- df[1:8,]
test.df <- df[9:10,]
train.mm <- makeModelMatrixFromDataFrame(train.df)
test.mm <- makeModelMatrixFromDataFrame(test.df, attr(train.mm, "drop"))
}
\keyword{factor}
