\name{guessNumCores}
\title{Guess Number of Cores}
\alias{guessNumCores}
\alias{defaultNumThreads}
\description{
  Attempts to guess the number of CPU \sQuote{cores}, both physical and logical, and gives the thread count the fitting functions default to.
}
\usage{
guessNumCores(logical = FALSE)

defaultNumThreads(n.tasks = Inf)
}
\arguments{
  \item{logical}{
    A logical value. When \code{FALSE}, an estimate of the number of physical cores is returned. When \code{TRUE}, so-called \dQuote{logical} cores are also included.
  }
  \item{n.tasks}{
    A positive number, the units of work that can run at once (for \code{\link{bart}}, the number of chains); \code{Inf} for no limit, which is also what a value that is no task count is read as, so that the function whose default this is can refuse it itself.
  }
}
\details{
  Because of different definitions of cores used by different manufacturers, the distinction between logical and physical cores is not universally recognized. This function will attempt to use operating system definitions when available, which should usually match the CPU itself.

  \code{defaultNumThreads} returns \code{guessNumCores()}, or \code{n.tasks} if that is smaller. When the environment variable \env{_R_CHECK_LIMIT_CORES_} is set to anything but \code{""} or \code{"false"}, as \command{R CMD check} sets it, the result is at most two. This cap applies to the default alone: an \code{n.threads} or \code{mc.cores} that is stated is never changed. \code{guessNumCores} is unaffected by the variable.
}
\value{
  An integer, or \code{NA} if no clear answer was obtained. The fitting functions take their default \code{n.threads} from this value, and where it is \code{NA} the default is one thread, with no message; an \code{n.threads} that is stated as \code{NA} is refused. \code{defaultNumThreads} returns an integer, \code{NA} where \code{guessNumCores()} is.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\examples{
## physical cores, the default used by dbartsControl(n.threads = )
guessNumCores()

## including hyperthreads
guessNumCores(logical = TRUE)

## the default n.threads of bart() with four chains
defaultNumThreads(4L)
}
\keyword{parallel}
