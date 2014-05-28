\name{makeind}
\title{Build x matrix from x data frame (convert factors to dummies)}
\alias{makeind}
\description{
Converts factors to dummies.\cr
Note that with more than one level, BART needs a dummy for each level of a factor
(unlike in linear regression where one of the dummies is dropped).
}
\usage{
makeind(
   x,
   all=TRUE)
}
\arguments{
   \item{x}{
   Data frame of explanatory variables.
   }

   \item{all}{
   If all=TRUE, a factor with p levels will be replaced by all p dummies.
   If all=FALSE, the pth dummy is dropped.
   }
}
\details{
   Uses function class.ind from the nnet library.
   Note that if you have train and test data frames, it may be best
   to rbind the two together, apply makeind to the result, and then pull them back apart.
}
\value{
   A matrix.\cr
   Numerical variables come first, and then the appended dummies.
}
\author{
Hugh Chipman: \email{hugh.chipman@acadiau.ca}\cr
Robert McCulloch: \email{robert.mcculloch@chicagogsb.edu}.
}
\seealso{
\code{\link{bart}}
}
\examples{
x1 = 1:10
x2 = as.factor(c(rep(1,5),rep(2,5)))
x3 = as.factor(c(1,1,1,2,2,2,4,4,4,4))
levels(x3) = c('rob','hugh','ed')
x = data.frame(x1,x2,x3)

junk = makeind(x)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
