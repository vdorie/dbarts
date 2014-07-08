\name{pdbart}
\title{Partial Dependence Plots for BART}
\alias{pdbart}
\alias{plot.pdbart}
\alias{pd2bart}
\alias{plot.pd2bart}
\description{
   Run \code{\link{bart}} at test observations constructed so that
   a  plot can be created
   displaying the effect of
   a single variable (\code{pdbart}) or pair of variables (\code{pd2bart}).
   Note that if \eqn{y} is a binary with \eqn{P(Y=1 | x) = F(f(x))}, \eqn{F} the standard
   normal cdf, then the plots are all on the \eqn{f} scale.
}
\usage{
   pdbart(
      x.train, y.train,
      xind = seq_len(ncol(x.train)),
      levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
      pl = TRUE, plquants = c(0.05, 0.95),
      \dots)
      
   \method{plot}{pdbart}(
      x,
      xind = seq_len(x$fd),
      plquants = c(0.05, 0.95), cols = c('black', 'blue'),
      \dots)
      
   pd2bart(
      x.train, y.train,
      xind = c(1, 2),
      levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
      pl = TRUE, plquants = c(0.05, 0.95),
      \dots)
      
   \method{plot}{pd2bart}(
      x,
      plquants = c(0.05, 0.95), contour.color = 'white',
      justmedian = TRUE,
      \dots)
}
\arguments{
  \item{x.train}{
    Explanatory variables for training (in sample) data. Must be a matrix of numeric type
    with rows corresponding to observations and columns to variables. Categorical
    variables/factors need to be converted to dummies, with a full set of columns present
    if there are more than two levels.
   }

  \item{y.train}{
    Dependent variable for training (in sample) data. Must be a numeric vector with length
    equal to the number of rows in \code{x.train}.
  }

   \item{xind}{
     Integer vector indicating which variables are to be plotted.
     In \code{pdbart}, corresponds to the variables (columns of \code{x.train}) for which a plot is to be constructed.
     In \code{plot.pdbart}, corresponds to the indices in list returned by \code{pdbart} for which plot is to be constructed.
     In \code{pd2bart}, the indicies of a pair of variables (columns of \code{x.train}) to plot.
   }

   \item{levs}{
     Gives the values of a variable at which the plot is to be constructed.
     Must be a list, where the \eqn{i}th component gives the values for the \eqn{i}th variable.
     In \code{pdbart}, it should have same length as \code{xind}. In \code{pd2bart}, it should have length 2.
     See also argument \code{levquants}.
   }

  \item{levquants}{
    If \code{levs} in \code{NULL}, the values of each variable used in the plot is
    set to the quantiles (in \code{x.train}) indicated by levquants.
    Must be a vector of numeric type.
  }

  \item{pl}{
    For \code{pdbart} and \code{pd2bart}, if \code{TRUE}, plot is subsequently made (by calling \code{plot.*}).
  }

  \item{plquants}{
    In the plots, beliefs about \eqn{f(x)} are indicated by plotting the
    posterior median and a lower and upper quantile.
    \code{plquants} is a double vector of length two giving the lower and upper quantiles.
  }

  \item{\dots}{
    Additional arguments.
    In \code{pdbart} and \code{pd2bart}, arguments are passed on to \code{\link{bart}}.
    In \code{plot.pdbart}, they are passed on to \code{\link{plot}}.
    In \code{plot.pd2bart}, they are passed on to \code{\link{image}}.
  }

  \item{x}{
    For \code{plot.*}, object returned from \code{pdbart} or \code{pd2bart}.
  }

  \item{cols}{
    Vector of two colors. The first color is for the median of \eqn{f}, while the second color
    is for the upper and lower quantiles.
  }

  \item{contour.color}{
    Color for contours plotted on top of the image.
  }

  \item{justmedian}{
    A logical where if \code{TRUE} just one plot is created for
    the median of \eqn{f(x)} draws. If \code{FALSE}, three plots are created
    one for the median and two additional ones for the lower and upper quantiles.
    In this case, \code{\link[par]{mfrow}} is set to \code{c(1,3)}.
  }
}
\details{
   We divide the predictor vector \eqn{x} into a subgroup of interest,
   \eqn{x_s} and the complement \eqn{x_c = x \setminus x_s}{x_c = x - x_s}.
   A prediction \eqn{f(x)} can
   then be written as \eqn{f(x_s, x_c)}. To estimate the effect of \eqn{x_s}
   on the prediction, Friedman suggests the partial dependence
   function
   \deqn{
   f_s(x_s) = \frac{1}{n}\sum_{i=1}^n f(x_s,x_{ic})
   }{
   f_s(x_s) = (1/n) \sum_{i=1}\^n f(x_s,x_{ic})
   }
   where \eqn{x_{ic}} is the \eqn{i}th observation of \eqn{x_c} in the data. Note
   that \eqn{(x_s, x_{ic})} will generally not be one of the observed data
   points. Using BART it is straightforward to then estimate and even
   obtain uncertainty bounds for \eqn{f_s(x_s)}. A draw of \eqn{f^*_s(x_s)}{f*_s(x_s)}
   from the induced BART posterior on \eqn{f_s(x_s)} is obtained by
   simply computing \eqn{f^*_s(x_s)}{f*_s(x_s)} as a byproduct of each MCMC draw
   \eqn{f^*}{f*}. The median (or average)
   of these MCMC draws \eqn{f^*_s(x_s)}{f*_s(x_s)} then yields an
   estimate of \eqn{f_s(x_s)}, and lower and upper quantiles can be used
   to obtain intervals for \eqn{f_s(x_s)}.

   In \code{pdbart} \eqn{x_s} consists of a single variable in \eqn{x} and in
   \code{pd2bart} it is a pair of variables.

   This is a computationally intensive procedure.
   For example, in \code{pdbart}, to compute the partial dependence plot
   for 5 \eqn{x_s} values, we need
   to compute \eqn{f(x_s, x_c)} for all possible \eqn{(x_s, x_{ic})} and there
   would be \eqn{5n} of these where \eqn{n} is the sample size.
   All of that computation would be done for each kept BART draw.
   For this reason running BART with \code{keepevery} larger than 1 (eg. 10)
   makes the procedure much faster.
}
\value{
  The plot methods produce the plots and don't return anything.

  \code{pdbart} and \code{pd2bart} return lists with components
  given below. The list returned by \code{pdbart} is assigned class
  \code{pdbart} and the list returned by \code{pd2bart} is assigned
  class \code{pd2bart}.

  \item{fd}{
    A matrix whose \eqn{(i, j)} value is the \eqn{i}th
    draw of \eqn{f_s(x_s)} for the \eqn{j}th value of \eqn{x_s}.
    \dQuote{fd} is for \dQuote{function draws}.

    For \code{pdbart} \code{fd} is actually a list whose
    \eqn{k}th component is the matrix described above
    corresponding to the \eqn{k}th variable chosen by argument xind.\
    The number of columns in each matrix will equal the number of values
    given in the corresponding component of argument \code{levs} (or number of values in \code{levquants}).

    For \code{pd2bart}, \code{fd} is a single matrix.
    The columns correspond to all possible pairs of values for the pair
    of variables indicated by \code{xind}.
    That is, all possible \eqn{(x_i, x_j)} where \eqn{x_i} is a value in
    the levs component corresponding to the first \eqn{x} and
    \eqn{x_j} is a value in the levs components corresponding to the second one.
    The first \eqn{x} changes first.
  }

  \item{levs}{
    The list of levels used, each component corresponding to a variable.
    If argument \code{levs} was supplied it is unchanged.
    Otherwise, the levels in \code{levs} are as constructed using argument \code{levquants}.
  }

  \item{xlbs}{
    A vector of character strings which are the plotting labels used for the variables.
  }
   
  The remaining components returned in the list are the same as in the value of \code{\link{bart}}.
  They are simply passed on from the BART run used to create the partial dependence plot.
  The function \code{\link{plot.bart}} can be applied to the object returned by \code{pdbart} or
  \code{pd2bart} to examine the BART run.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2006)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.

both of the above at:
\url{http://faculty.chicagogsb.edu/robert.mcculloch/research/rob-mcculloch-cv.html}

Friedman, J.H. (2001)
   Greedy function approximation: A gradient boosting machine.
     \emph{The Annals of Statistics}, \bold{29}, 1189--1232.
}
\author{
Hugh Chipman: \email{hugh.chipman@acadiau.ca}.\cr
Robert McCulloch: \email{robert.mcculloch@chicagogsb.edu}.
}
\examples{
## simulate data 
f <- function(x) { return(0.5 * x[,1] + 2 * x[,2] * x[,3]) }

sigma <- 0.2
n     <- 100

set.seed(27)
x <- matrix(2 * runif(n * 3) -1, ncol = 3);
colnames(x) <- c('rob', 'hugh', 'ed')

Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

## first two plot regions are for pdbart, third for pd2bart
par(mfrow = c(1, 3))

## pdbart: one dimensional partial dependence plot
set.seed(99)
pdb1 <-
  pdbart(x, y, xind = c(1, 2),
         levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
         pl = FALSE, keepevery = 10, ntree = 100)
plot(pdb1, ylim = c(-0.6,.6))

## pd2bart: two dimensional partial dependence plot
set.seed(99)
pdb2 <-
  pd2bart(x, y, xind = c(2, 3),
          levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
          pl = FALSE, ntree = 100, keepevery = 10, verbose = FALSE)
plot(pdb2)

## compare BART fit to linear model and truth = Ey
lmFit  <- lm(y ~., data.frame(x, y))
fitmat <- cbind(y, Ey, lmFit$fitted, pdb1$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
\keyword{dplot}
