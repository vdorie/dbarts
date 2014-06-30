\name{bart}
\alias{bart}
\alias{plot.bart}
\title{Bayesian Additive Regression Trees}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner.

\itemize{
  \item For numeric response \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon},
where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary response \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi}
denotes the standard normal cdf (probit link).
}
}
\usage{
bart(x.train, y.train, x.test = matrix(0.0, 0, 0),
     sigest = NA, sigdf = 3, sigquant = 0.90,
     k = 2.0,
     power = 2.0, base = 0.95,
     binaryOffset = 0.0,
     ntree = 200,
     ndpost = 1000, nskip = 100,
     printevery = 100, keepevery = 1, keeptrainfits = TRUE,
     usequants = FALSE, numcut = 100, printcutoffs = 0,
     verbose = TRUE)
\method{plot}{bart}(x,
     plquants = c(0.05, 0.95), cols = c('blue', 'black'),
     \dots)
}
\arguments{
   \item{x.train}{
     Explanatory variables for training (in sample) data. May be a matrix or a data frame,
     with rows corresponding to observations and columns to variables.
     If a variable is a factor in a data frame, it is replaced with dummies.
     Note that \eqn{q} dummies are created if \eqn{q > 2} and
     one dummy is created if \eqn{q = 2}, where \eqn{q} is the number of levels of the factor.
   }

   \item{y.train}{
     Dependent variable for training (in sample) data. If \code{y.train} is numeric a continous response
     model is fit (normal errors). If \code{y.train} is a binary factor or has only values 0 and 1, then
     a binary response model with a probit link is fit.
   }

   \item{x.test}{
     Explanatory variables for test (out of sample) data. Should have same column structure as
     \code{x.train}. \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a
     row of \code{x.test}.
   }

   \item{sigest}{
   	 For continuous response models, an estimate of the error variance, \eqn{\sigma^2},
   	 used to calibrate an inverse-chi-squared prior used on that parameter. If not supplied,
   	 the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary.
   }

   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary.
   }

   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate
     (\code{sigest}) is placed at. The closer the quantile is to 1,
     the more aggresive the fit will be as you are putting more prior weight
     on error standard deviations (\eqn{\sigma}) less than the rough estimate.
     Not applicable when \eqn{y} is binary.
   }

   \item{k}{
     For numeric \eqn{y},
     \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}.
     The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}.
     For binary \eqn{y}, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from
     \eqn{\pm 3}{+/- 3}. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be.
   }

   \item{power}{
     Power parameter for tree prior.
   }

   \item{base}{
     Base parameter for tree prior.
   }

   \item{binaryOffset}{
     Used for binary \eqn{y}. When present, the model is
       \eqn{P(Y = 1 \mid x) = \Phi(f(x) + \mathrm{binaryOffset})}{P(Y = 1 | x) = \Phi(f(x) + binaryOffset)},
     allowing fits with probabilities shrunk towards values other than \eqn{0.5}.
   }

   \item{ntree}{
     The number of trees in the sum-of-trees formulation.
   }

   \item{ndpost}{
     The number of posterior draws after burn in, \code{ndpost / keepevery} will actually be returned.
   }

   \item{nskip}{
     Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
     As the MCMC runs, a message is printed every \code{printevery} draws.
   }

   \item{keepevery}{
     Every \code{keepevery} draw is kept to be returned to the user. Useful for \dQuote{thinning}
     samples.
   }

   \item{keeptrainfits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of \code{x.train} are returned.
   }

   \item{usequants}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived
     from the \code{x.train} variables. When \code{FALSE}, splits are determined using
     values equally spaced across the range of a variable.
     %Decision rules in the tree are of the form
     %\eqn{x \le c}{x \le c}  vs. \eqn{x > c}
     %for each variable corresponding to a column of \code{x.train}.
     %\code{usequants} determines how the set of possible \eqn{c} is determined.
     %If \code{usequants} is \code{TRUE}, then the \eqn{c} are
     %a subset of the values \code{(xs[i] + xs[i + 1]) / 2} where \code{xs} is
     %unique sorted values obtained from
     %the corresponding column of \code{x.train}.
     %If \code{usequants} is \code{FALSE}, the cutoffs are equally spaced across the
     %range of values taken on by the corresponding
     %column of \code{x.train}.
   }
   
   \item{numcut}{
   	 The number of possible values used in decision rules (see \code{usequants}). If a single
   	 number, it is recycled for all variables; otherwise must be a vector of length equal to
   	 \code{ncol(x.train)}. Can be truncated if more values are requested than present in the
   	 training data.
   	 
%     The number of possible values of c (see usequants).
%   If a single number if given, this is used for all variables.
%   Otherwise a vector with length equal to ncol(x.train) is required,
%   where the \eqn{i^{th}}{i^th} element gives the number of c used for
%   the \eqn{i^{th}}{i^th} variable in x.train.
%   If usequants is false, numcut equally spaced cutoffs
%   are used covering the range of values in the corresponding
%   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
%   corresponding columns of x.train - 1) c values are used.
   }

   \item{printcutoffs}{
     The number of cutoff rules to printed to screen before the MCMC is run.
     Given a single integer, the same value will be used for all variables.
     If 0, nothing is printed.
   }

   \item{verbose}{
     Logical; if \code{FALSE} supress printing.
   }

   \item{x}{
     Object of class \code{bart}, returned by function \code{bart}, which contains the information to be plotted.
   }

   \item{plquants}{
     In the plots, beliefs about \eqn{f(x)} are indicated by plotting the
     posterior median and a lower and upper quantile. \code{plquants} is a double vector
     of length two giving the lower and upper quantiles.
   }

   \item{cols}{
     Vector of two colors. First color is used to plot the median of \eqn{f(x)} and
     the second color is used to plot the lower and upper quantiles.
   }

   \item{\dots}{
     Additional arguments passed on to plot.
   }
}
\details{
   BART is an Bayesian MCMC method. 
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   including the burn-in draws.  Initially these draws will decline as BART finds fit
   and then level off when the MCMC has burnt in.\cr
   The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   the corresponding \eqn{f(x)} on the vertical axis. 
   

   \code{bart} returns a list assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }
   
   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }
 
   \item{y}{
   The input dependent vector of values for the dependent variable.\cr
   This is used in plot.bart.
   }

   In the binary \eqn{y} case, the returned list has the components
   yhat.train, yhat.test, and varcount as above.  In addition the list 
   has a binaryOffset component giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2009)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning. 
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

both of the above at:
\url{http://www.rob-mcculloch.org}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\note{
\emph{There was a bug in BayesTree\_0.1-0} (now fixed of course).\cr
If the number of test observations was less than the number of trees
(200 is the default), the yhat.test and yhat.test.mean components were suspect.
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com}\cr
Robert McCulloch: \email{robert.mcculloch1@gmail.com}.
}
\seealso{
\code{\link{pdbart}}
}
\examples{
##simulate data (example from Friedman MARS paper)
f = function(x){
10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y)
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
