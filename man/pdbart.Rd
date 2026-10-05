\name{pdbart}
\title{Partial Dependence Plots for BART}
\alias{pdbart}
\alias{plot.pdbart}
\alias{pd2bart}
\alias{plot.pd2bart}
\description{
   Fit a model with \code{\link{bart}}, or take a fit or sampler, and average its predictions over the fit's rows with one variable (\code{pdbart}) or a pair of variables (\code{pd2bart}) set to each of a grid of values, so that a plot can display that variable's effect. The averages are on the scale of the model's linear predictor: if \eqn{y} is binary with \eqn{P(Y=1 | x) = F(f(x))}, \eqn{F} the standard normal cdf, then the plots are all on the \eqn{f} scale.
}
\usage{
pdbart(
    formula, data,
    xind = NULL,
    levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
    pl = TRUE, plquants = c(0.05, 0.95),
    \dots)

\method{plot}{pdbart}(
    x,
    xind = seq_along(x$fd),
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    \dots)

pd2bart(
    formula, data,
    xind = NULL,
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
  \item{formula}{
    The data to fit, as \code{\link{bart}}'s own \code{formula}: a matrix, data frame or formula. Also accepted are a fitted \code{bart} model kept with \code{keepTrees = TRUE}, or a \code{\link[=dbartsSampler-class]{dbartsSampler}}, passed first. A \code{bart} model kept without its trees or sampler is refit from its stored call by the function that made it, with a warning of class \code{dbartsFallbackWarning}, and is refused if no call was kept; a sampler without \code{keepTrees} is run, generating fresh samples and changing its state, with the same warning.
  }
  \item{data}{
    As \code{\link{bart}}'s \code{data}: the response when \code{formula} is a matrix, or the data frame in which to evaluate a formula. Not used with a fit or sampler.
  }
  \item{xind}{
    Integer, character vector, or the right-hand side of a formula indicating which variables are to be plotted. In \code{pdbart}, corresponds to the variables (columns of the predictor matrix) for which a plot is to be constructed. In \code{plot.pdbart}, corresponds to the indices in list returned by \code{pdbart} for which plot is to be constructed. In \code{pd2bart}, the indices of a pair of variables (columns of the predictor matrix) to plot. If \code{NULL} a default of all columns is used for \code{pdbart} and the first two columns is used for \code{pd2bart}.
  }
   \item{levs}{
     Gives the values of a variable at which the plot is to be constructed. Must be a list, where the \eqn{i}th component gives the values for the \eqn{i}th variable. In \code{pdbart}, it should have same length as \code{xind}. In \code{pd2bart}, it should have length 2. The values of a factor predictor are level names. See also argument \code{levquants}.
   }
  \item{levquants}{
    If \code{levs} is \code{NULL}, the values of each variable used in the plot are set to the quantiles (in the training data, missing values left out) indicated by levquants, and a factor predictor takes every level. Must be a vector of numeric type.
  }
  \item{pl}{
    For \code{pdbart} and \code{pd2bart}, if \code{TRUE}, plot is subsequently made (by calling \code{plot.*}).
  }
  \item{plquants}{
    In the plots, beliefs about \eqn{f(x)} are indicated by plotting the posterior median and a lower and upper quantile. \code{plquants} is a double vector of length two giving the lower and upper quantiles.
  }
  \item{\dots}{
    In \code{pdbart} and \code{pd2bart}, arguments of \code{\link{bart}}, under its names and with its defaults; an argument \code{bart} does not take is an error. \code{pdbart} sets \code{keepTrees} and \code{samplerOnly} itself and has no test rows, so \code{keepTrees = FALSE}, \code{samplerOnly}, \code{test} and \code{offset.test} are refused; \code{keepSampler = FALSE} drops the sampler from the result. With a fit or sampler passed in, only \code{keepSampler} is taken. BayesTree's spellings, as \code{\link{bartBT}} takes them (\code{x.train}, \code{y.train}, \code{ntree}, \code{ndpost}, \code{nskip}, \code{keepevery}, \code{binaryOffset} and the rest), are translated to \code{bart}'s with a once-per-session warning until dbarts 1.1-0, when they are refused; \code{power}, \code{base} and \code{splitprobs} become \code{tree.prior = cgm()}, and \code{proposalprobs} the control's \code{proposal.probs}. A setting given under both spellings is refused.

    In \code{plot.pdbart}, arguments are passed on to \code{\link{plot}}, except \code{type}, which is passed to the lines drawn, and \code{xlab} and \code{ylab}, which replace the method's own labels. In \code{plot.pd2bart}, they are passed on to \code{\link{image}}.
  }
  \item{x}{
    For \code{plot.*}, object returned from \code{pdbart} or \code{pd2bart}.
  }
  \item{cols}{
    Vector of two colors. The first color is for the median of \eqn{f}, while the second color is for the upper and lower quantiles.
  }
  \item{contour.color}{
    Color for contours plotted on top of the image.
  }
  \item{justmedian}{
    A logical where if \code{TRUE} just one plot is created for the median of \eqn{f(x)} draws. If \code{FALSE}, three plots are created one for the median and two additional ones for the lower and upper quantiles. In this case, \code{\link[=par]{mfrow}} is set to \code{c(1,3)}.
  }
}
\details{
   We divide the predictor vector \eqn{x} into a subgroup of interest, \eqn{x_s} and the complement \eqn{x_c = x \setminus x_s}{x_c = x - x_s}. A prediction \eqn{f(x)} can then be written as \eqn{f(x_s, x_c)}. To estimate the effect of \eqn{x_s} on the prediction, Friedman suggests the partial dependence function \deqn{f_s(x_s) = \frac{1}{n}\sum_{i=1}^n f(x_s,x_{ic})}{f_s(x_s) = (1/n) \sum_{i=1}\^n f(x_s,x_{ic})} where \eqn{x_{ic}} is the \eqn{i}th observation of \eqn{x_c} in the data. Note that \eqn{(x_s, x_{ic})} will generally not be one of the observed data points. Using BART it is straightforward to then estimate and even obtain uncertainty bounds for \eqn{f_s(x_s)}. A draw of \eqn{f^*_s(x_s)}{f*_s(x_s)} from the induced BART posterior on \eqn{f_s(x_s)} is obtained by simply computing \eqn{f^*_s(x_s)}{f*_s(x_s)} as a byproduct of each MCMC draw \eqn{f^*}{f*}. The median (or average) of these MCMC draws \eqn{f^*_s(x_s)}{f*_s(x_s)} then yields an estimate of \eqn{f_s(x_s)}, and lower and upper quantiles can be used to obtain intervals for \eqn{f_s(x_s)}.

   In \code{pdbart} \eqn{x_s} consists of a single variable in \eqn{x} and in \code{pd2bart} it is a pair of variables.

   The rows averaged over are the fit's own: those it was fit to, less any it gives a weight of 0, which for a binary response masks a row out. Each row's offset, \code{bart}'s \code{offset}, is added to its prediction before averaging, so a value equals the row mean of \code{\link[=predict.bart]{predict}(fit, newdata, type = "bart")} with the variable set.

   A data call fits through \code{\link{bart}} at its defaults, keeping trees and sampler, and predicts each grid value from the saved trees, so \code{pdbart(x, y, seed = s)} and \code{pdbart(bart(x, y, seed = s, keepTrees = TRUE))} agree. \code{pdbart(bartBT(x, y, keeptrees = TRUE))} gives the model dbarts 0.9-34 fit. The first data call in a session prints a message saying so; a call from package code does not.

   \code{pdbart} and \code{pd2bart} serve gaussian, Student-t, probit and logistic fits. Multinomial and ordinal fits, whose prediction is a probability per category, are refused; \code{predict} on rows with the variable set gives each category's. Negative binomial, hurdle, accelerated failure time and hazard fits are refused in this version. Each refusal comes before anything is fit, whether the family is named, resolved by \code{family = "auto"} from the response, or that of a fit or sampler passed in.

   This is a computationally intensive procedure. For example, in \code{pdbart}, to compute the partial dependence plot for 5 \eqn{x_s} values, we need to compute \eqn{f(x_s, x_c)} for all possible \eqn{(x_s, x_{ic})} and there would be \eqn{5n} of these where \eqn{n} is the sample size. All of that computation would be done for each kept BART draw. For this reason thinning the draws with \code{n.thin} larger than 1 (eg. 10) makes the procedure much faster.
}
\value{
  The plot methods produce the plots and don't return anything.

  \code{pdbart} and \code{pd2bart} return lists with components given below. The list returned by \code{pdbart} is assigned class \code{pdbart} and the list returned by \code{pd2bart} is assigned class \code{pd2bart}. The class carries a \code{plot} method only; it is not a fit, so \code{predict}, \code{extract}, \code{fitted}, and \code{residuals} are not defined for it (\code{fitted} and \code{residuals} fall through to \pkg{stats}' defaults and return \code{NULL}).

  \item{fd}{
    A matrix whose \eqn{(i, j)} value is the \eqn{i}th draw of \eqn{f_s(x_s)} for the \eqn{j}th value of \eqn{x_s}. \dQuote{fd} is for \dQuote{function draws}.

    For \code{pdbart} \code{fd} is actually a list whose \eqn{k}th component is the matrix described above corresponding to the \eqn{k}th variable chosen by argument \code{xind}. The number of columns in each matrix will equal the number of values given in the corresponding component of argument \code{levs} (or number of values in \code{levquants}).

    For \code{pd2bart}, \code{fd} is a single matrix. The columns correspond to all possible pairs of values for the pair of variables indicated by \code{xind}. That is, all possible \eqn{(x_i, x_j)} where \eqn{x_i} is a value in the levs component corresponding to the first \eqn{x} and \eqn{x_j} is a value in the levs components corresponding to the second one. The first \eqn{x} changes first.

    The draws of every chain are merged, each chain's in turn. When the fit keeps its chains apart (\code{combineChains = FALSE} with more than one chain), each matrix gains a leading chain margin, chains x draws x values, as the fit's own components do.
  }
  \item{levs}{
    The list of levels used, each component corresponding to a variable. If argument \code{levs} was supplied it is unchanged. Otherwise, the levels in \code{levs} are as constructed using argument \code{levquants}, or a factor predictor's level names. A factor predictor is plotted with one point and interval per level, and without contours in \code{plot.pd2bart}.
  }
  \item{xlbs}{
    A vector of character strings which are the plotting labels used for the variables.
  }
  \item{n.chains}{
    The number of chains.
  }
   
  The remaining components are passed on from the fit, under \code{\link{bart}}'s own names (see its \sQuote{Value} section) with one rename and several omissions - the fit's call is under \code{bartcall} rather than \code{call}, and not every \code{bart} component is carried over. For a continuous response the remaining components are \code{bartcall}, \code{yhat.train}, \code{first.sigma}, \code{sigma}, \code{yhat.train.mean}, \code{sigest}, \code{y}, and \code{fit}, the sampler, which \code{keepSampler = FALSE} leaves out. For a binary response only \code{bartcall}, \code{yhat.train}, \code{y}, and \code{fit} are present. \code{pd2bart} follows the same rules. The function \code{\link{plot.bart}} can be applied to the object returned by \code{pdbart} or \code{pd2bart} to examine the BART run.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2010)
   BART: Bayesian additive regression trees.
        \emph{The Annals of Applied Statistics}, \bold{4}(1), 266--298.
        \doi{10.1214/09-AOAS285}.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
   \url{https://www.rob-mcculloch.org/}

Friedman, J.H. (2001)
   Greedy function approximation: A gradient boosting machine.
     \emph{The Annals of Statistics}, \bold{29}, 1189--1232.
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com}.\cr
Robert McCulloch: \email{robert.mcculloch1@gmail.com}.
}
\examples{
## simulate data
f <- function(x)
    return(0.5 * x[,1] + 2 * x[,2] * x[,3])

sigma <- 0.2
n     <- 100

set.seed(27)
x <- matrix(2 * runif(n * 3) - 1, ncol = 3)
colnames(x) <- c('rob', 'hugh', 'ed')

Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

\donttest{
## pdbart: one dimensional partial dependence plot
pdb1 <- pdbart(
    x, y, xind = c(1, 2),
    levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
    pl = FALSE, n.thin = 10, n.trees = 100, n.chains = 2, n.threads = 1,
    seed = 99, verbose = FALSE
)

## pd2bart: two dimensional partial dependence plot
pdb2 <- pd2bart(
    x, y, xind = c(2, 3),
    levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
    pl = FALSE, n.thin = 10, n.trees = 100, n.chains = 2, n.threads = 1,
    seed = 99, verbose = FALSE)

## the plot methods draw into the current device, so set up the layout and
## put the caller's graphical parameters back afterwards. The first two plot
## regions are for pdbart, the third for pd2bart.
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(pdb1, ylim = c(-0.6, 0.6))
plot(pdb2)
par(oldpar)

## compare BART fit to linear model and truth = Ey
lmFit  <- lm(y ~ ., data.frame(x, y))
fitmat <- cbind(y, Ey, lmFit$fitted, pdb1$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
}

\donttest{
## example showing the use of a pre-fitted model
df <- data.frame(y, x)
bartFit <- bart(
    y ~ rob + hugh + ed, df,
    n.thin = 10, n.trees = 100, n.chains = 2, n.threads = 1, seed = 99,
    keepTrees = TRUE, verbose = FALSE)
pdb3 <- pdbart(bartFit, xind = rob + ed, pl = FALSE)
}}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
\keyword{dplot}
