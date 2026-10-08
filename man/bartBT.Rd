\name{bartBT}
\alias{bartBT}
\alias{plot.bart}
\alias{print.bart}
\alias{extract}
\alias{predict.bart}
\alias{extract.bart}
\alias{fitted.bart}
\alias{residuals.bart}
\alias{family.bart}
\title{Bayesian Additive Regression Trees}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner.

\itemize{
  \item For numeric response \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary response \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi} denotes the standard normal cdf (probit link).
  \item \code{bartBT} is the BayesTree-compatible surface: 0.9-34's argument list exactly, no new features, and no \code{family} argument at all - every fit is the default gaussian/probit pair above. \code{\link{bart}} reaches every other response family - logistic, accelerated failure time, discrete-time survival hazard, K-category multinomial, ordered categorical, negative-binomial counts, and semicontinuous two-part responses - plus its own formula-term multi-forest interface.
}
}
\usage{
bartBT(
    x.train, y.train, x.test = matrix(0.0, 0, 0),
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2.0,
    power = 2.0, base = 0.95, splitprobs = 1 / numvars,
    binaryOffset = 0.0, weights = NULL,
    ntree = 200,
    ndpost = 1000, nskip = 100,
    printevery = 100, keepevery = 1, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100, printcutoffs = 0,
    verbose = TRUE, nchain = 1, nthread = 1, combinechains = TRUE,
    keeptrees = FALSE, keepcall = TRUE, sampleronly = FALSE,
    seed = NULL,
    proposalprobs = NULL,
    keepsampler = keeptrees)

\method{plot}{bart}(
    x,
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    \dots)

\method{predict}{bart}(
    object, newdata,
    type = c("ev", "ppd", "bart", "forest", "sigma"),
    offset = NULL, weights = NULL,
    combineChains = TRUE,
    ci.level = NULL,
    forest = NULL,
    bases = NULL,
    na.action = dbarts::na.keepPredictors,
    n.threads,
    \dots)

\method{print}{bart}(x, \dots)

extract(object, \dots)
\method{extract}{bart}(
    object,
    type = c("ev", "ppd", "bart", "loglik", "trees", "forest",
             "sigma", "k", "leaf.prior.sd", "varcount"),
    sample = c("train", "test"),
    combineChains = TRUE,
    forest = NULL,
    contribution = FALSE, \dots)

\method{fitted}{bart}(
    object,
    type = c("ev", "ppd", "bart", "forest"),
    sample = c("train", "test"),
    ci.level = NULL,
    \dots)

\method{residuals}{bart}(object, type = "ev", \dots)

\method{family}{bart}(object, \dots)
}
\arguments{
   \item{x.train}{
     Explanatory variables for training (in sample) data. May be a matrix or a data frame, with rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that \eqn{q} dummies are created if \eqn{q > 2} and one dummy is created if \eqn{q = 2}, where \eqn{q} is the number of levels of the factor. A sparse \code{Matrix::dgCMatrix} of ordinal predictors is also accepted, as is a data frame mixing ordinary columns with \code{Matrix::sparseVector} or \code{dgCMatrix} columns; see \code{\link{dbarts}}.
   }
   \item{y.train}{
     Dependent variable for training (in sample) data. If \code{y.train} is numeric a continuous response model is fit (normal errors). If \code{y.train} is a two-level factor (or logical, or two-level character) or has only values 0 and 1, then a binary response model with a probit link is fit. A factor with three or more levels is an error, whose message names both remedies: \code{\link{bart}(x, y, family = "multinomial")} for a multinomial fit, or \code{as.integer(y) - 1L} to keep the numbers dbarts 0.9-x fit from the level codes.
   }
   \item{x.test}{
     Explanatory variables for test (out of sample) data. Should have same column structure as \code{x.train}. \code{bartBT} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test}.
   }
   \item{sigest}{
     For continuous response models, an estimate of the residual standard deviation (residual standard error), \eqn{\sigma}, used to calibrate an inverse-chi-squared prior on the error variance. If not supplied, the least-squares estimate is derived instead. That estimate falls back to the marginal standard deviation of the response when the linear model's residual standard error comes out non-finite, warning as it does so (class \code{dbartsSigmaFallbackWarning}); a design with sparse-backed predictor columns skips the linear model altogether and falls back the same way (class \code{dbartsSparseSigmaFallbackWarning}, a \code{dbartsSigmaFallbackWarning}). See \code{sigquant} for more information. Not applicable when \eqn{y} is binary. Same concept as \code{sigest} in \code{\link{dbarts}} and \code{\link{xbart}}.
   }
   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary. Default 3, Chipman, George, and McCulloch's calibration (see References); an aggressive choice rather than a derived one.
   }
   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate (\code{sigest}) is placed at. The closer the quantile is to 1, the more aggressive the fit will be as you are putting more prior weight on error standard deviations (\eqn{\sigma}) less than the rough estimate. Not applicable when \eqn{y} is binary. Default 0.90, from the same source as \code{sigdf}.
   }
   \item{k}{
     For numeric \eqn{y}, \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}. The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}. For binary \eqn{y} fit with the default probit link, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from \eqn{\pm 3}{+/- 3}, the probit reference scale; \code{family = "logistic"} widens this to \eqn{\pm \pi \sqrt{3}}{+/- pi*sqrt(3)}, which is three standard deviations of the standard logistic latent variable (\eqn{\pi / \sqrt{3}}{pi/sqrt(3)} each), the same three-sd span probit's \eqn{\pm 3}{+/- 3} covers. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be. The value can be either a fixed number, or a \emph{hyperprior} of the form \code{chi(df = 1.5, scale = 2)}. \code{bartBT}'s own default is the fixed value 2; see \code{\link{bart}}'s \code{k} item for its \code{NULL} default, a \code{chi(1.5, 2)} hyperprior on binary responses. The default of 2 for continuous responses follows Chipman, George, and McCulloch's argument (see References) that with leaf prior standard deviation \eqn{\sigma_\mu = 0.5 / (k \sqrt{m})}{sigma_mu = 0.5 / (k * sqrt(m))} for \eqn{m} trees, \eqn{k} prior standard deviations of \eqn{f(x)} span the whole coded response range regardless of \eqn{m} -- so \eqn{k = 2} places that range at roughly a 95\% prior interval. The range scaling itself is an outlier-sensitive convention shared with BayesTree/bartMachine: extreme \eqn{y} values compress the effective prior on everything else, for which the published workaround is to log-transform or winsorize such values before fitting, or to let \code{k} adapt via the \code{chi} hyperprior.
   }
   \item{power}{
     Power parameter for tree prior. Default 2, Chipman, George, and McCulloch's empirical recommendation (see References) for the split-probability decay \eqn{base (1 + depth)^{-power}}{base * (1 + depth)^(-power)}, not derived from a formula.
   }
   \item{base}{
     Base parameter for tree prior. Default 0.95, from the same source as \code{power}.
   }
   \item{splitprobs}{
     Prior and transition probabilities of variables used to generate splits. Can be missing/empty/\code{NULL} for equiprobability, a numeric vector of length equal to the number variables, or a named numeric vector with only a subset of the variables specified and a \code{.default} named value. Values given for factor variables are replicated for each resulting column in the generated model matrix. The symbol \code{numvars} is rebound before execution to the number of columns in the model matrix.
   }
   \item{binaryOffset}{
     Used for binary \eqn{y}. When present, the model is \eqn{P(Y = 1 \mid x) = \Phi(f(x) + \mathrm{binaryOffset})}{P(Y = 1 | x) = \Phi(f(x) + binaryOffset)}, allowing fits with probabilities shrunk towards values other than \eqn{0.5}. \eqn{\Phi} is the link for the default probit family; a \code{family = "logistic"} fit uses \code{plogis} instead, as noted in the value section below.
   }
   \item{weights}{
      An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. A probit fit (\code{bartBT}, or \code{\link{bart}} with the default binary family) does not support weights, except that weights identically 1 are treated as absent, and weights that are all 0 and 1 name the rows in the data set and install as the fit's active-row mask, a 0 row leaving the likelihood but keeping its leaf occupancy, its latent and its fitted value (the same holds for \code{family = "ordinal"}); a \code{family = "logistic"} fit treats them as observation counts and requires positive integers no larger than \eqn{10^6}{1e6}, as a sweep's time grows with the total count. For a weighted logistic fit, the \code{"ppd"} draw at an observation with weight \eqn{w} is the number of successes among \eqn{w} trials, \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} with \eqn{p} the fitted probability.
   }
   \item{ntree}{
     The number of trees in the sum-of-trees formulation. Where the data object carries \code{bases}, it is the count of the forest with no basis, a forest with one taking 50 trees, and it is refused, when given, where every forest has a basis; see \code{\link{forest}}. Like every count on this page (\code{nskip}, \code{ndpost}, \code{numcut}, \code{keepevery}, \code{printevery}, \code{nthread}, \code{nchain}, \code{printcutoffs}), a fractional value is refused, naming the argument, rather than silently truncated; see \code{\link{dbartsControl}} for the shared whole-number rule.
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
     Every \code{keepevery} draw is kept to be returned to the user. Useful for \dQuote{thinning} samples.
   }
   \item{keeptrainfits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of \code{x.train} are returned.
   }
   \item{usequants}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived from the \code{x.train} variables. When \code{FALSE}, splits are determined using values equally spaced across the range of a variable. See details for more information.
   }
   \item{numcut}{
     The maximum number of possible values used in decision rules (see \code{usequants}, details). If a single number, it is recycled for all variables; otherwise must be a vector of length equal to \code{ncol(x.train)}. Fewer rules may be used if a covariate lacks enough unique values.
   }
   \item{printcutoffs}{
     The number of cutoff rules to printed to screen before the MCMC is run. Given a single integer, the same value will be used for all variables. If 0, nothing is printed.
   }
   \item{verbose}{
     Logical; if \code{FALSE} suppress printing.
   }
   \item{nchain}{
     Integer specifying how many independent tree sets and fits should be calculated. Default 1 for \code{bartBT}, BayesTree's historical single-chain convention; \code{\link{bart}}'s default is 4.
   }
   \item{nthread, n.threads}{
     A positive integer giving a total thread budget, distinct from \code{nchain}: chains run in parallel across it, but tree sampling itself uses at most one thread per chain, so a budget above \code{nchain} buys the sampler's own sweep nothing. The surplus still feeds the test-fit pool and \code{predict}'s own fan-out, and a fit warns once, naming both counts, when an explicit budget exceeds the chain count. \code{bartBT}'s default is 1.

     On \code{predict}, \code{n.threads} is a per-call worker count for the saved-tree replay, defaulting to the fit's own: the replay is partitioned by (chain, posterior draw), each partition writing its own rows and nothing being reduced across workers, so the answer is identical bit for bit at every value and only the time taken changes. It is the last positional argument of every \code{predict} method, so a caller that passes arguments by position is unaffected by its presence. A value that is not a single positive number is an error rather than a silent fallback.
   }
   \item{combinechains, combineChains}{
     Logical; if \code{TRUE}, samples from every chain are collapsed into a single matrix/vector of \code{nchain} \eqn{\times} \code{ndpost} rows (see \sQuote{Value} for the row order); if \code{FALSE}, they are kept in arrays of dimensions equal to \code{nchain} \eqn{\times} \code{ndpost} \eqn{\times} number of observations. Default \code{TRUE} across \code{bartBT} and \code{\link{bart}}; the number of chains is tracked regardless on the returned object's \code{n.chains} component (see \sQuote{Value}). On \code{extract}, \code{predict} and \code{\link{survivalProbabilities}}, \code{FALSE} keeps that leading \code{nchain} margin - length 1 on a one-chain fit - rather than dropping it there; a scalar type (\code{sigma}, \code{k}, a component's \code{shape}) is then a 1 x samples matrix instead of a plain vector. A fit's own STORED \code{combinechains = FALSE} fields keep 0.9-x's one-chain collapse instead, so \code{$yhat.train} and its siblings never gain that margin on their own.
   }
   \item{keeptrees}{
     Logical; must be \code{TRUE} in order to use \code{predict} with the result of a \code{bart} fit. Note that for models with a large number of observations or a large number of trees, keeping the trees can be very memory intensive.
   }
   \item{keepcall}{
     Logical; if \code{FALSE}, returned object will have \code{call} set to \code{NULL}, otherwise the call used to instantiate BART.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. \code{NULL} (the default) leaves the fit unseeded, and an explicit \code{NA} is a missing value, as in \code{\link{set.seed}}: it reads as \code{NULL} for one release, with a once-per-session warning; a \code{\link{set.seed}} beforehand then suffices for reproducibility. Supplying \code{seed} instead gives reproducible results without touching R's stream. See the Reproducibility section of \link{bart}.
   }
   \item{proposalprobs}{
     Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"swap"}, \code{"change"}, \code{"perturb"} and \code{"rule_gibbs"} to control tree structure proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step. \code{NULL} (the default), or a vector naming none of \code{"birth_death"}, \code{"swap"} and \code{"change"}, takes BayesTree's mixture for those three, \code{c(birth_death = 0.5, swap = 0.1, change = 0.4)}, as 0.9-x did (\code{birth} defaults to 0.5 either way); any other vector is completed as \code{\link{dbartsControl}} completes one. \code{\link{bart}}'s default differs: \code{c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0, rule_gibbs = 0, birth = 0.5)}. All five structural probabilities zero is the frozen mixture: no structural proposal is made, the tree structures stand where they are, and only the leaf values, \code{sigma} and the family's latents keep being drawn, which is how a fitted forest is re-sampled as a fixed basis. Under \code{\link{dbartsControl}}'s default \code{treeShift = "auto"} a frozen forest additionally takes the level-shifting Gibbs step each iteration, the leaf values then being the only thing left to move. \code{bartBT} has no \code{treeShift} formal of its own and always runs at this default; \code{\link{bart}} and a \code{\link{dbarts}} sampler's control object can set it to \code{"always"} or \code{"never"} instead, which changes the sampled values without changing the posterior. A \code{"swap"} element exchanges a parent's split rule with a child's; \code{bart} ships it at zero because at production forest sizes it measures as a no-op, and \code{bartBT} keeps it at 0.1 to match BayesTree. With \code{n.trees = 1} it is the only move that rotates a rule up the tree. A \code{"perturb"} element displaces one node's split point by a single cut position while keeping its variable and the tree's shape; it defaults to zero, and only ordinal (numeric) columns can be perturbed. A \code{"rule_gibbs"} element replaces one nog node's rule - a node whose two children are both leaves - with a draw from that rule's own full conditional over the available ordinal variables and their admissible cuts, so its acceptance is one; it defaults to zero, it acts only where the node's own rule is ordinal, and it is inert on an all-categorical design.
   }
   \item{keepsampler}{
     Logical that can be used to save the underlying \code{\link{dbartsSampler-class}} object even if \code{keeptrees} is false.
   }
   \item{object}{
     An object of class \code{bart}, returned from either the function \code{bartBT} or \code{bart}.
   }
   \item{newdata}{
     Test data for prediction. Obeys all the same rules as \code{x.train}. An unnamed matrix given to a fit whose predictors are named is matched by position, with a warning. Left out (or \code{NULL}), \code{predict} returns each \code{type} at the training rows from the draws the fit stored, as \code{extract} with \code{sample = "train"} does, and needs no saved trees; an \code{offset} given then replaces the fit's own at the training rows, one value or one per training row, while \code{weights} and \code{bases} are refused there.
   }
   \item{offset}{
     For \code{predict}: an optional offset for \code{newdata}, applied the same way \code{\link{dbarts}}'s \code{offset.test} is at fit time. It shifts the PREDICTED rows, not the training ones - a different channel from \code{offset.test}, which is why it carries its own name rather than that one's. \code{NULL} (the default, when missing) applies none. Refused, by name, with \code{type = "forest"}: an offset shifts the recombination of the forests, never any one forest's own total, so it belongs to the arithmetic under \sQuote{Value} rather than to what comes back. On an amplitude-coupled fit's combined arms (\code{type = "ev"}, \code{"ppd"}, \code{"bart"}) it is that recombination, and enters the sum of trees before the family's link, exactly as it does for a single forest.
   }
   \item{sampleronly}{
     Builds the sampler from its arguments and returns it without running it. Useful to use the \code{bartBT} interface in more complicated models.
   }
   \item{x}{
     Object of class \code{bart}, returned by function \code{bartBT}, which contains the information to be plotted.
   }
   \item{plquants}{
     In the plots, beliefs about \eqn{f(x)} are indicated by plotting the posterior median and a lower and upper quantile. \code{plquants} is a double vector of length two giving the lower and upper quantiles.
   }
   \item{cols}{
     Vector of two colors. First color is used to plot the median of \eqn{f(x)} and the second color is used to plot the lower and upper quantiles.
   }
   \item{type}{
     The quantity to be returned by generic functions. Options are \code{"ev"} - samples from the posterior of the individual level expected value, \code{"bart"} - the sum of trees component; same as \code{"ev"} for linear models but on the latent scale (probit or logistic, matching the fit's family) for binary ones, \code{"ppd"} - samples from the posterior predictive distribution, \code{"loglik"} - for \code{extract} only, the log-likelihood of each training observation at each posterior draw (a row whose weight is a count has its value multiplied by the count; a row with a precision weight of zero, or masked out of the data, is not in the likelihood and is \code{NaN}, so drop those columns before handing the matrix to a tool such as \pkg{loo}), \code{"trees"} - a data frame with tree information for when model was fit with \code{keepTrees} equal to \code{TRUE}, and \code{"forest"} - for \code{extract}, \code{predict} and \code{fitted}, the per-forest channels of an amplitude-coupled multi-forest fit (see \code{forest}, \code{contribution}, and the \code{forestFits}/\code{glue}/\code{bases} components under \sQuote{Value}); an error naming the reason on any other fit, a \code{\link{bart}} multinomial one included, whose K forests are per-category latents rather than additive components of one location. \code{extract(type = "forest")} reads the stored in-sample channel, while \code{predict(type = "forest")} replays each forest at \code{newdata}; both report the raw per-forest total, leaving the recombination to the caller, where the combined arms (\code{"ev"}, \code{"ppd"}, \code{"bart"}) perform it given the bases at those rows (see \code{bases}). \code{fitted(type = "forest")} is the posterior mean of what \code{extract} reports, a rows-by-forests matrix with \code{extract}'s dimnames, pooled over all samples and chains and refused wherever \code{extract} refuses it; it is each forest's own total BEFORE its basis and amplitude multiply it, so it is not that forest's contribution to the fitted value, and \code{residuals} refuses it. For \code{"ppd"}, a weighted logistic fit draws the number of successes among the observation-count weight, \eqn{\mathrm{Binomial}(w_i, p_i)}{Binomial(w_i, p_i)} (see \code{weights}), and \code{weights} given to \code{predict} on a probit or logistic fit are likewise trial counts, positive integers, each new row drawing \eqn{\mathrm{Binomial}(w_i, p_i)}{Binomial(w_i, p_i)}, an aft (survival) fit draws on the log-time scale (the fit models \eqn{\log T}{log T}), a \code{family = student()} fit draws \eqn{f + \sigma t_\nu / \sqrt{w_i}}{f + sigma t_nu / sqrt(w_i)} at each draw's own \eqn{\sigma}{sigma} and \eqn{\nu}{nu} (the fit's \code{$resid.df}), and a heteroscedastic fit (\code{\link{dbarts}}'s \code{variance}) draws its noise at that observation's own \eqn{s(x_i)/\sqrt{w_i}}{s(x_i)/sqrt(w_i)}. For \code{"loglik"}, gaussian fits evaluate \eqn{y_i \mid x_i \sim N(\hat{f}(x_i), \sigma^2 / w_i)}{y_i | x_i ~ N(f(x_i), sigma^2 / w_i)} in logs at each draw of \eqn{f} and \eqn{\sigma}{sigma} - or, for a heteroscedastic fit, \eqn{y_i \mid x_i \sim N(\hat{f}(x_i), s^2(x_i) / w_i)}{y_i | x_i ~ N(f(x_i), s^2(x_i) / w_i)} at each draw of \eqn{f} and \eqn{s}, a \code{family = student()} fit the corresponding marginal \eqn{t_\nu}{t_nu} density at that draw's \eqn{\nu}{nu} (the fit's \code{$resid.df}) rather than the normal one, binary fits evaluate the Bernoulli log-likelihood of the fitted probability, multiplied for a weighted logistic fit by the observation-count weight \eqn{w_i}, and an aft fit contributes the log density for an event and the log survival tail \eqn{\log P(T > C)}{log P(T > C)} for a right-censored observation, both on the log-time scale. When chains are combined the result is a samples-by-observations matrix directly consumable by WAIC/PSIS-LOO implementations such as those in the \pkg{loo} package; the chains-first convention is kept, so a per-chain array (\code{combineChains = FALSE}, dimension chains-by-samples-by-observations) is reordered to the draws-by-chains-by-observations that \code{loo::relative_eff} expects with \code{aperm(x, c(2, 1, 3))}. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"ev"} and \code{"link"} can be used as a synonym for \code{"bart"}; the \code{\link{bart}} extended-family methods take the same two synonyms, each against its own set of types. \code{extract} also takes \code{"sigma"}, \code{"k"} and \code{"leaf.prior.sd"}, the residual scale, the leaf prior's \code{k} and its standard deviation, and \code{"varcount"}, the per-predictor split counts shaped as \code{"ev"} is with predictors in place of observations - none of them per-observation, so a supplied \code{sample} is refused by name on all of them. A parameter the fit sampled comes back as its draws, a chain-combined vector or (\code{combineChains = FALSE}) a chains-by-samples matrix; one the fit held fixed comes back as that one number, under either \code{combineChains}, and \code{forest} selects among the forests of a fit with several (see \code{forest} below and the \sQuote{Leaf prior parameter \code{k}} subsection). \code{sigma} on a weighted fit is the scale at weight 1 (row i's is \eqn{\sigma / \sqrt{w_i}}{sigma / sqrt(w_i)}); on a \code{family = student()} fit it is the t scale, not the standard deviation. A fit whose family has no residual scale to estimate - probit, logistic, hazard, and the ordinal, negative-binomial and multinomial fits of \code{\link{bart}} - answers \code{"sigma"} with 1, as \code{\link[stats]{sigma}} does on a binomial \code{glm}. A heteroscedastic fit's \code{"sigma"} is its per-observation scale \eqn{s(x)}, laid out as \code{"ev"} is and taking \code{sample}: \code{extract} returns the stored \code{s.train} or \code{s.test}, and \code{predict}, which takes \code{"sigma"} on such a fit only, evaluates the variance forest at \code{newdata} (refusing \code{weights}, since a case weight does not change \eqn{s(x)}). For information on extracting trees, see the subsection below.
   }
   \item{sample}{
     Either \code{"train"} or \code{"test"}. \code{"test"} is refused, by name, when \code{type} is \code{"forest"}: an amplitude-coupled fit has no test fits (see \code{\link{dbarts}}'s \code{forests =}). It is \code{extract}'s and \code{fitted}'s own argument - \code{fitted} is always the training rows - and is refused by name on \code{predict}, whose stored train and test channels are \code{extract}'s \code{sample} instead; also refused by name on \code{extract(type = "sigma")}/\code{"k"}/\code{"leaf.prior.sd"}/\code{"varcount"}, none of them per-observation (\code{"sigma"} excepted on a heteroscedastic fit, where it is).
   }
   \item{forest}{
     For \code{extract(type = "forest")}, \code{extract(type = "trees")} and \code{predict(type = "forest")}: which forest(s) to return, by 1-based position, as base R reads a list element by number, or by label (a string; see \code{\link{forest}}), \code{"forest1"}, \code{"forest2"}, \dots being the position names every per-forest margin carries; \code{NULL} (the default) returns every forest. For \code{type = "forest"} the returned array always keeps the trailing forest margin, subset to the requested forests, even when only one is selected; \code{type = "trees"} returns a data frame and keeps no such margin. For \code{extract(type = "k")} and \code{extract(type = "leaf.prior.sd")} on a fit with several forests, it selects among the named per-forest numbers (one number for a single forest). Selecting a forest outside those is refused by name: every other arm has already recombined the forests into the location it reports, except \code{"sigma"} (a model parameter with no forest axis) and \code{"varcount"} (which already keeps every forest on its own trailing margin), each naming its own reason instead.
   }
   \item{bases}{
     For \code{predict} on an amplitude-coupled multi-forest fit: the bases each forest's amplitudes multiply AT the predicted rows, which off the training rows only the caller knows. A length-K list, in the forests' order, entry k either \code{NULL} (for a forest that declared no basis, whose amplitude multiplies an implicit all-ones column) or anything \code{\link{forest}}'s \code{basis} accepts as a value - a numeric vector or matrix, a factor, a character or a logical vector - expanded by the same rule and required to have \code{nrow(newdata)} rows and the width that forest's amplitudes take. A name on an entry is refused unless it is that forest's label or \code{forest<i>}, since the list is read by position. When exactly one forest carries a basis the bare value may be given on its own, without the list: for a Bayesian causal forest, \code{bases = cbind(1 - zstar, zstar)} is the counterfactual arm being predicted under. Unlike at fit time, a column of all zeros and a single-level factor are accepted, since a constant arm is the point of a counterfactual. \code{NULL} (the default) builds each basis again from \code{newdata} when the fit's forest wrote it as code, in a \code{\link{bart}} or \code{\link{dbarts}} formula's \code{forest()} term or in a \code{forests} list (it then needs the columns of the fit's data that basis names, coded against the fit's own factor levels; see \code{\link{forest}}'s \code{basis} item), and is an error naming the forest when the basis was handed over as a value or needs a variable that had a value for every fitted row and is no column of \code{newdata}. Refused on a single-forest fit, and with \code{type = "forest"}, which reports each forest's total before any basis.
   }
   \item{contribution}{
     For \code{extract(type = "forest")} only: \code{FALSE} (the default) returns the raw per-forest total \code{forestFits} already carries (see \sQuote{Value}). \code{TRUE} instead returns each selected forest's per-observation CONTRIBUTION to the fit, \eqn{(\mathrm{basis}_k \, \mathrm{glue}_k) \times \mathrm{forestFits}_k}{(basis_k \%*\% glue_k) * forestFits_k}, computed on demand from the stored \code{forestFits}, \code{glue}, and \code{bases} rather than stored itself.
   }
   \item{na.action}{
     For \code{predict}: what to do with the rows of \code{newdata} that have missing values; see \sQuote{At Prediction} in \code{\link{na.keepPredictors}}.
   }
   \item{ci.level}{
     For \code{fitted} and \code{predict}, an optional single number in \eqn{(0, 1)}. When \code{NULL} (the default), each generic returns what it always has: \code{fitted} the posterior mean, \code{predict} the full array of posterior samples (see \code{type} above), unsummarized. When given, the result is instead a matrix with one row per observation and columns \code{est} (the posterior mean), \code{ci.lower}, and \code{ci.upper} giving a symmetric credible band at that level, pooled over all samples and chains. The kind of interval follows \code{type}: \code{"ev"} yields a credible interval for the conditional mean (a probability for binary responses), \code{"ppd"} a prediction interval that additionally carries the residual noise (and so is wider), and \code{"bart"} a credible interval on the latent scale. For \code{fitted(type = "forest")} it returns the rows-by-forests-by-3 array of \code{est}, \code{ci.lower} and \code{ci.upper}, as for a multinomial fit; \code{predict(type = "forest")} refuses it by name, since that arm reports each forest's own total before any basis. It is also refused on \code{residuals}, where it would band a constant minus an increasing triple and so mislabel two of its three columns; take \code{y - fitted(object, ci.level = )} instead, where the reversal is visible at the call site.
   }
   \item{\dots}{
     Additional arguments passed on to \code{plot} (via \code{plot.bart}), or to the sampler's \code{getTrees} method via \code{extract} when \code{type} is \code{"trees"} (\code{chainNums}, \code{sampleNums}, \code{treeNums}, \code{newdata}, \code{forest}; see \sQuote{Extracting Trees} below) - the one place a name belonging to a sibling method is still accepted rather than refused. Elsewhere, a name that is a formal on \code{predict}, \code{extract}, \code{fitted}, or \code{residuals} but foreign to the method actually called is refused by name rather than silently discarded. Any other unrecognized name warns instead - a warning, not an error, so that a subclass method forwarding its own formals through \code{\link{NextMethod}} is never refused for an argument its caller legitimately supplied.
   }
}
\details{
  BART is an Bayesian MCMC method. At each MCMC interation, we produce a draw from the joint posterior \eqn{(f, \sigma) \mid (x, y)}{(f, \sigma) | (x, y)} in the numeric \eqn{y} case and just \eqn{f} in the binary \eqn{y} case.

  Thus, unlike a lot of other modeling methods in R, \code{bartBT} does not produce a single model object from which fits and summaries may be extracted. The output consists of values \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{\sigma*} in the numeric case) where * denotes a particular draw. The \eqn{x} is either a row from the training data (\code{x.train}) or the test data (\code{x.test}).

  \subsection{Decision Rules}{
    Decision rules for any tree are of the form \eqn{x \le c}  vs. \eqn{x > c} for each \sQuote{\eqn{x}} corresponding to a column of \code{x.train}. \code{usequants} determines the means by which the set of possible \eqn{c} is determined. If \code{usequants} is \code{TRUE}, then the \eqn{c} are a subset of the values interpolated half-way between the unique, sorted, finite values obtained from the corresponding column of \code{x.train}. If \code{usequants} is \code{FALSE}, the cutoffs are equally spaced across the range of finite values taken on by the corresponding column of \code{x.train}; an infinite value falls past an end cutoff, as it does on the quantile grid.

    The number of possible values of \eqn{c} is determined by \code{numcut}. If \code{usequants} is \code{FALSE}, \code{numcut} equally spaced cutoffs are used covering the range of values in the corresponding column of \code{x.train}; a cutoff is held once, so a column whose range is too narrow to hold that many distinct ones gets fewer, and a constant column gets one, where dbarts 0.9-34 repeated them. If \code{usequants} is \code{TRUE}, then for a variable the minimum of \code{numcut} and one less than the number of unique elements for that variable are used. When the half-way values outnumber \code{numcut}, the ones used are spread evenly over all of them: with \eqn{M} half-way values and \eqn{m} = \code{numcut}, cutoff \eqn{k}, counting from one, is the half-way value at position \eqn{\lfloor (2k - 1) M / (2m) \rfloor + 1}{floor((2k - 1) M / (2m)) + 1}, so the cutoffs reach both ends of the variable. dbarts 0.9-34 and BayesTree stepped through them from the bottom at a fixed stride and could stop well short of the top, leaving that part of the variable without a cutoff; fits with \code{usequants = TRUE} on such a variable therefore generally differ from theirs.
  }
  \subsection{Leaf prior parameter \code{k}}{
    The amount of shrinkage of the leaf parameters is controlled by \code{k}. \code{k} can be given as either a fixed, positive number, or as any value that can be used to build a supported hyperprior. At present, only \eqn{\chi_\nu s} priors are supported, where \eqn{\nu} is a degrees of freedom and \eqn{s} is a scale. Both values must be positive, however the scale can be infinite which yields an improper prior, which is interpreted as just the polynomial part of the distribution. If \eqn{\nu} is 1 and \eqn{s} is \eqn{\infty}, the prior is \dQuote{flat}.

    For BART on binary outcomes, the degree of overfitting can be highly sensitive to \code{k} so it is encouraged to consider a number of values. The default hyperprior for binary \code{bart}/\code{dbarts} fits, \code{chi(1.5, 2)}, centers the sampled \code{k} near the field-standard fixed value of 2 (prior median 1.9) while adapting to the data; pass \code{k = 2} for the fixed BART-package default. Crossvalidation may still be helpful, and running for a short time with a flat prior (\code{scale = Inf}) can show the range of \code{k} values the data are consistent with.

    A fit stores the \code{k} its sampler recorded, whatever terms the leaf prior was named in (\code{k}, or \code{sd} through \code{leaf.prior}; see \code{\link{dbartsPriors}}); \code{extract(type = "k")} returns it on every fit. Under a \code{k}-named prior it is relative to the data's scale, and under an \code{sd}-named one relative to the \code{k.scale} the fit stores (\code{leaf.prior$k.scale}), not to the data. \code{extract(type = "leaf.prior.sd")} returns the \code{k.scale} over \code{k}: the prior standard deviation of the forest's total, which is the spread the leaf prior states, in the units the forest fits - the response's on a gaussian or Student-t fit, log time on an accelerated failure time fit, the latent scale of the link on a probit, logistic or ordinal fit, and the log mean on a negative-binomial fit. What the standard deviation is the standard deviation of depends on the leaf model: the leaf value for a constant leaf, each coefficient for a linear leaf and the amplitude for a Gaussian-process leaf. Under a monotone constraint the total's prior spread runs a few percent wider than this number. On a fit with several forests it is each forest's own total before its amplitude multiplies it, one named number per forest (\code{forest1}, \code{forest2}, ...), or the one that \code{forest} selects; \code{k} is pinned at 1 on each. A hurdle fit returns a list of its zero and positive parts, each in its own units.
  }
  \subsection{Generics}{
    \code{bartBT} supports \code{\link{fitted}} to return the posterior mean of a predicted quantity, as well as \code{\link{predict}} to return a set of posterior samples for a different sample. In addition, the \code{extract} generic can be used to obtain the posterior samples for the training data or test data supplied during the initial fit.

    Using \code{predict} with a \code{bart} object requires that it be fitted with the option \code{keeptrees}/\code{keepTrees} as \code{TRUE}. Keeping the trees for a fit can require a sizeable amount of memory and is off by default.

    All generics return values on the scale of expected value of the response by default. This means that \code{predict}, \code{extract}, and \code{fitted} for binary outcomes return probabilities unless specifically the sum-of-trees component is requested (\code{type = "bart"}). This is in contrast to \code{yhat.train}/\code{yhat.test} that are returned with the fitted model.

    \code{residuals} returns \eqn{y} minus \code{fitted(object, type = type, sample = "train")}; the default \code{type = "ev"} gives ordinary response-scale residuals, so for a binary fit this is \eqn{y - P(Y = 1 \mid x)}{y - P(Y = 1 | x)}. Binary fits do not have \code{yhat.train.mean}/\code{yhat.test.mean} components (see \sQuote{Value}); use \code{fitted} to obtain the equivalent posterior-mean probabilities directly. \code{residuals} is always against the training response, so \code{sample} is refused by name; \code{ci.level} is refused too (see its own entry above).
  }

  \subsection{Saving}{
    \code{\link{save}}ing and \code{\link{load}}ing a fitted BART object for use with \code{predict} requires that it be fit with \code{keeptrees}/\code{keepTrees} as \code{TRUE}, and that the sampler's tree and generator state be captured before saving by calling \code{storeState()} on the sampler: \code{bartFit$fit$storeState()}. The state is not captured automatically because it duplicates the trees on the R side and can exceed the fit's own sample blocks; it is materialized only on request. A fit saved without it reloads, but \code{predict}, \code{extract(type = "trees")}, and \code{plotTree} all stop identically, each with an error naming \code{storeState()}. The same convention covers \code{\link{bart}}'s own-class families - \code{family = "multinomial"}, \code{"ordinal"}, and \code{"nbinom"} - whose \code{$fit} is the sampler (K-forest or single-forest) that actually ran; see \code{\link{bart}}'s Value section for each. A \code{family = "hurdle.lognormal"} fit has no \code{$fit} of its own: it is the pair \code{$zero} and \code{$positive}, each a \code{bart} fit in its own right, so both of their samplers need the call before saving: \code{fH$zero$fit$storeState(); fH$positive$fit$storeState()}.
  }

  \subsection{Reproducibility}{
    See the Reproducibility section of \link{bart}: it covers \code{seed} versus \code{\link{set.seed}}, chains and thread count, and continuing, copying, or restoring a sampler, in full for both \code{bartBT} and \code{\link{bart}}.
  }

  \subsection{Extracting Trees}{
    When a model is fit with \code{keeptrees} (\code{bartBT}) or \code{keepTrees} (\code{bart}) equal to \code{TRUE}, the generic \code{extract} can be used to retrieve a data frame containing the tree fit information. In this case, \code{extract} will accept the additional, optional arguments: \code{chainNums}, \code{sampleNums}, and \code{treeNums}. Each should be an integer vector detailing the desired trees to be returned. A further optional argument \code{newdata} routes a new set of predictors (in the same form accepted by \code{predict}) through the frozen trees so that the \code{n} column counts those observations instead of the training data. A last optional argument \code{forest} - a single position or label, a vector of them, or the default \code{NULL} for every forest (a label is a string, as \code{extract(type = "trees", forest = "dose")} reads it; see \code{\link{forest}}) - selects among a fit's forests (only an amplitude-coupled multi-forest fit - one more forest per additional \code{forest()} term or \code{forests =} entry - or a \code{"bartMultinomial"} fit's K categories, carry more than one); see \code{\linkS4class{dbartsSampler}}'s \code{forest} item.

    The result of \code{extract} will be a data frame with columns:
    \itemize{
      \item \code{forest} - which forest the row belongs to, from 1; present only on a fit with several forests (a \code{"bartMultinomial"} fit, or one declared with \code{forests =}, whatever the count); a single-forest fit has no such column, and still accepts \code{forest = 1}
      \item \code{chain}, \code{sample}, \code{tree} - index variables; \code{chain} is always present, 1 on a single-chain fit, since tree draws cannot be combined across chains (the sampler's own \code{getTrees} omits it at one chain)
      \item \code{n} - number of observations in node, drawn from the training data or from \code{newdata} when supplied
      \item \code{var} - either the index of the variable used for splitting or -1 if the node is a leaf
      \item \code{value} - either the value such that observations less than or equal to it are sent down the left path of the tree or the predicted value for a leaf node
    }
    The order of nodes in the result corresponds to a depth-first traversal, going down the left-side first. The names of variables used in splitting can be recovered by examining the column names of the \code{fit$data@x} element of a fitted \code{bartBT} or \code{bart} model. See the package vignette \dQuote{Working with dbarts Saved Trees}.

    A fit made with \code{keeptrees}/\code{keepTrees} \code{FALSE} but a sampler kept anyway (\code{keepSampler = TRUE}, or a \code{bart} fit with \code{fit} present per \sQuote{Value} below) still answers \code{extract(fit, "trees")}, the same fallback \code{\link{plotTree}} uses: the frame holds the sampler's CURRENT trees rather than a saved history, so it carries no \code{sample} column and reflects whatever predictors and structure the sampler holds at the time of the call, not a fixed draw.
  }
}
\value{
  \code{bartBT} and \code{\link{bart}} return lists assigned the class \code{bart} for the shared gaussian/probit/logistic/aft/hazard* families (\code{bart} reaches further response families under their own list classes and generics; see \code{\link{bart}}'s Value section). For applicable quantities, \code{ndpost / keepevery} samples are returned. In the numeric \eqn{y} case, the list has components:

  \item{\code{call}}{
    The matched call, or \code{NULL} when the fit was made with \code{keepcall}/\code{keepCall} equal to \code{FALSE}.
  }
  \item{\code{yhat.train}}{
    A array/matrix of posterior samples. The \eqn{(i, j, k)} value is the \eqn{j}th draw of the posterior of \eqn{f} evaluated at the \eqn{k}th row of \code{x.train} (i.e. \eqn{f^*(x_k)}{f(x_k)}) corresponding to chain \eqn{i}. When \code{nchain} is one or \code{combinechains} is \code{TRUE}, the result is collapsed down to a matrix: with \eqn{m = }\code{ndpost \%/\% keepevery} draws kept per chain, row \eqn{r} is chain \eqn{1 + (r - 1) \%/\% m}'s draw \eqn{1 + (r - 1) \%\% m} - chain 1's whole run first, then chain 2's, and so on (chain-major), not interleaved draw by draw. Every collapsed field a fit returns (\code{sigma}, \code{k}, \code{varcount}, ...) uses this same row order, so row \eqn{r} of one pairs with row \eqn{r} of another - e.g. \code{sigma[r]} is the residual scale that produced \code{yhat.train[r, ]}.
  }
  \item{\code{yhat.test}}{
    Same as \code{yhat.train} but now the \eqn{x}s are the rows of the test data.
  }
  \item{\code{yhat.train.mean}}{
    Vector of means of \code{yhat.train} across columns and chains, with length equal to the number of training observations.
  }
  \item{\code{yhat.test.mean}}{
    Vector of means of \code{yhat.test} across columns and chains.
  }
  \item{\code{sigma}}{
    Matrix of posterior samples of \code{sigma}, the residual/error standard deviation; absent on a heteroscedastic fit, which has no scalar residual scale and whose scale draws are \code{s.train}/\code{s.test}. Dimensions are equal to the number of chains times the number of samples unless \code{nchain} is one or \code{combinechains} is \code{TRUE}, in which case it collapses to a vector in the same chain-major order as \code{yhat.train} (see above), so combined \code{sigma[r]} pairs with combined \code{yhat.train[r, ]}.
  }
  \item{\code{first.sigma}}{
    Burn-in draws of \code{sigma}; absent where \code{sigma} is.
  }
  \item{\code{varcount}}{
    A matrix with number of rows equal to the number of kept draws and each column corresponding to a training variable. Contains the total count of the number of times that variable is used in a tree decision rule (over all trees). Row order (when combined) is \code{yhat.train}'s chain-major order (see above).
  }
  \item{\code{varprobs}}{
    Present only under a DART prior: the sampled split-variable probabilities, in the same layout as \code{varcount}. Each draw sums to one across variables.
  }
  \item{\code{forestFits}, \code{glue}, \code{bases}}{
    Present only for an amplitude-coupled multi-forest fit: built through \code{\link{dbarts}}'s \code{forests =}, or through \code{\link{bart}}'s \code{forest()} formula term (see its \sQuote{Formula Terms} section); not reachable from \code{bartBT} directly, which has no multi-forest front door. \code{forestFits} is a (\code{n.chains} \eqn{\times}, when uncombined) \code{n.samples} \eqn{\times} number of training observations \eqn{\times} K array: forest k's own RESPONSE-scale raw total at each draw, i.e. \eqn{\mathrm{response.scale} \times f_k(x)}{response.scale * f_k(x)}, with NO amplitude (\code{glue}) folded in - the same quantity the sampler's \code{$getForestFits} reports, up to that one scalar. The trailing margin is named \code{forest1}, ..., \code{forestK}; each forest's label (a name in \code{forests}, else the text of the forest's \code{basis}, else \code{forest<i>}; see \code{\link{dbarts}}'s \code{forests} item) rides separately as a \code{"forest.labels"} attribute on the fit object. \code{glue} is the ragged (\code{n.chains} \eqn{\times}, when uncombined) \code{n.samples} \eqn{\times} sum(q_k) matrix of amplitudes multiplying each forest's basis columns, forest-major (the layout \code{$getForestAmplitudes} documents); its \code{"forest"} attribute names each column's forest. \code{bases} is a length-K list of each forest's expanded basis matrix (no draw axis; \code{NULL} for a forest with no declared basis, whose amplitude then scales an implicit all-ones column). Together they reconstruct the training fit exactly: \eqn{\hat y = \mathrm{response.shift} + \sum_k (\mathrm{bases}_k \, \mathrm{glue}_k) \times \mathrm{forestFits}_k}{yhat = response.shift + sum_k (bases_k \%*\% glue_k) * forestFits_k}, with \code{response.scale}/\code{response.shift} read from the sampler's \code{$getLeafPrior}. See \code{\link{extract}}'s \code{type = "forest"}. \code{predict(object, newdata, type = "forest")} reports the same quantity at NEW rows - an (\code{n.chains} \eqn{\times}, when uncombined) \code{n.samples} \eqn{\times} \code{nrow(newdata)} \eqn{\times} K array, same trailing margin, same \code{forest} selection - replayed from the saved trees, and so requiring \code{keeptrees}/\code{keepTrees}. It carries no \code{glue}, no \code{response.shift} and no offset either, because off the training rows the bases are the caller's, so an offset is refused there rather than folded in and there is no \code{contribution} argument for the same reason. The identity itself is what \code{predict(object, newdata, type = "ev")} (or \code{"ppd"}, or \code{"bart"}) performs at those rows in one call, taking each \code{bases} entry's replacement - the basis that forest's amplitudes multiply at the new rows, for a Bayesian causal forest the \eqn{(1 - z^*, z^*)}{(1 - z*, z*)} indicator pair of whichever assignment is being predicted under - through its own \code{bases} argument, or built again from \code{newdata} for a basis that was written as code. Written out by hand it remains the fallback for a caller holding only these components and no sampler. A fit one of whose bases was written as code, in a \code{forest()} term of a formula or in a \code{forests} list, additionally carries \code{basis.terms}, a length-K list holding for each such basis what builds it again at new rows (R's \code{\link[stats]{terms}} object for it, the levels of its factors, its text and the columns of the fit's data it names) and \code{NULL} for a forest with no basis or with one given as a value; it is absent when no basis was written as code.
  }
  \item{\code{n.forests}}{
    The forest count K, on every fit: 1 for a single forest.
  }
  \item{\code{row.names.train}, \code{row.names.test}}{
    The names of the training and test rows, which every per-observation output carries; absent when the rows are unnamed. See the row names paragraph of \code{\link{bart}}'s \sQuote{Value} section.
  }
  \item{\code{sigest}}{
    The rough error standard deviation (\eqn{\sigma}) used in the prior.
  }
  \item{\code{resid.scale}, \code{resid.df}}{
    \code{resid.scale} is present only on a family with a residual law (\code{"gaussian"}, \code{"student"}, \code{"aft"}): the law's scale model, \code{"constant"} or \code{"forest"} for a heteroscedastic fit (\code{variance} on \code{\link{bart}}), whose draws are \code{s.train}/\code{s.test}. The law's shape is the family's own (see \code{\link{dbarts}}'s \code{family = student()}). \code{resid.df} is present only for a \code{student()} fit: the degrees of freedom \eqn{\nu} each draw was conditioned on, in \code{sigma}'s layout - a fixed \eqn{\nu} repeats its value, an estimated one gives that draw's grid value. \code{extract(type = "loglik")} reads both.
  }
  \item{\code{y}}{
    The input dependent vector of values for the dependent variable. This is used in \code{plot.bart}.
  }
  \item{\code{fit}}{
    Optional sampler object which stores the values of the tree splits. Required for using \code{predict} and only stored if \code{keeptrees} or \code{keepsampler} is \code{TRUE}.
  }
  \item{\code{n.chains}}{
    Information that can be lost if \code{combinechains} is \code{TRUE} is tracked here.
  }
  \item{\code{k}}{
    Optional matrix of posterior samples of \code{k}, as the sampler recorded them whatever terms the leaf prior was named in. Only present when \code{k} is modeled, i.e. there is a hyperprior. Collapses (when combined) in \code{sigma}'s chain-major order.
  }
  \item{\code{first.k}}{
    Burn-in draws of \code{k}, if modeled.
  }
  \item{\code{leaf.prior}}{
    The leaf prior the fit ran under, on every fit: the sampler's own reading as the run ended (\code{\link{dbartsSampler}}'s \code{getLeafPrior}), a list holding the prior as it was named, what its standard deviation is the standard deviation of, and \code{k.scale}, the value \code{k} is relative to. A list of such lists, named \code{forest1}, \code{forest2}, ..., on a fit with several forests.
  }
  \item{\code{fixed}}{
    The scalars the sampler held fixed, a named list on every fit, empty when there are none: among \code{sigma}, \code{k} and \code{resid.df}, and \code{shape} on a negative-binomial fit, each read from the sampler. A fixed \code{sigma} is the sampler's, the square root of the variance \code{fixed()} names; a fixed \code{k} on a fit with several forests is a vector named by forest. Their draw channels, where there are any, repeat the value; \code{extract} and \code{summary} read this list in their place.
  }
  \item{\code{binaryOffset}}{
    Present only for a binary fit: the offset value used.
  }
  \item{\code{family}}{
    The response family as specified, once \code{"auto"} has resolved: \code{"gaussian"}, \code{"student"}, \code{"probit"}, \code{"logistic"}, \code{"aft"}, \code{"hazard.probit"}, or \code{"hazard.logistic"} (the bare token \code{"hazard"} is recorded as \code{"hazard.probit"}). \code{predict}, \code{extract}, \code{fitted}, and \code{plot} look its link up from it to transform latent draws to probabilities, and stop on a token they do not know. For an \code{"aft"} fit, predictions and fitted values are on the LOG-TIME scale: these functions return the linear predictor \eqn{E[\log T \mid x]}, exactly as for a gaussian fit of \eqn{\log T}, and never the time scale that \code{survreg} or \code{flexsurv} users might expect. The \code{type} aliases \code{"response"} and \code{"link"} are inert for \code{"aft"} - both return the same log-scale linear predictor; neither exponentiates. The posterior-median survival time is \code{exp} of the linear predictor, and \code{\link{survivalProbabilities}} gives survival-probability draws. A discrete-time hazard fit (\code{family = "hazard"}/\code{"hazard.logistic"}, reachable only through \code{\link{bart}}) records the binary link here (\code{"probit"}/\code{"logistic"}) - the fit IS that binary fit on the expanded rows - and carries a separate \code{$periods} element (the ordered period grid); \code{predict}/\code{extract} return the per-(subject, period) hazard through the link, and \code{\link{survivalProbabilities}} dispatches its survival-curve branch on the \code{$periods} marker.
  }

  \item{\code{family.spec}}{
    The same family as a family object carrying its settings (a Student-t \code{df}, a hazard link and grid, a residual prior), such as \code{student(df = 3)} or \code{hazard(link = "probit")}. \code{family(object)} returns it.
  }

  In the binary \eqn{y} case, the returned list has the components \code{yhat.train}, \code{yhat.test}, and \code{varcount} as above, but not \code{yhat.train.mean}/\code{yhat.test.mean} or \code{resid.scale} - use \code{fitted} to get the posterior mean of \eqn{P(Y = 1 \mid x)}{P(Y = 1 | x)} instead. In addition the list has a \code{binaryOffset} component giving the value used.

  Note that in the binary \eqn{y}, case \code{yhat.train} and \code{yhat.test} are \eqn{f(x) + \mathrm{binaryOffset}}{f(x) + binaryOffset}. For draws of the probability \eqn{P(Y = 1 | x)}, apply the family's inverse link to these values: the normal cdf (\code{pnorm}) for \code{"probit"} fits, \code{plogis} for \code{"logistic"} ones.

  For continuous response fits, the \code{plot} method sets \code{mfrow} to \code{c(1, 2)} and makes two plots, unless \code{sigma} was held fixed, when it draws only the second. The first plot is the sequence of kept draws of \eqn{\sigma} including the burn-in draws. Initially these draws will decline as BART finds a good fit and then level off when the MCMC has burnt in. The second plot has \eqn{y} on the horizontal axis and posterior intervals for the corresponding \eqn{f(x)} on the vertical axis. For binary response fits, only this second kind of plot is drawn, with the posterior median of \eqn{P(Y = 1 \mid x)} on the horizontal axis and its posterior interval on the vertical axis.

  For \code{print.bart}, the fit itself (\code{x}), returned invisibly; the call and a short synopsis (family, chain/tree/burn counts, kept-draw count) print to the console as the side effect.
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
   \url{https://www.rob-mcculloch.org}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com},
Robert McCulloch: \email{robert.mcculloch1@gmail.com},
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{bart}} for the formula-first interface, further response families, and formula-term multi-forest models.

\code{\linkS4class{dbartsSampler}}'s \code{$predictForests} for the same per-forest replay from a mutable sampler, which is what \code{predict(type = "forest")} reads and what \code{extract(type = "forest")} reports in sample.

\code{\link{pdbart}}

\code{\link{dbarts-embedding}} and \code{vignette("dbarts-as-a-component", package = "dbarts")} for
using BART as one block of a sampler you write, rather than as the whole fit: these functions own
their own MCMC loop, while \code{\link{dbarts}} returns a
\code{\linkS4class{dbartsSampler}} that runs one sweep at a time.
}
\examples{
## simulate data (example from Friedman MARS paper)
## y = f(x) + epsilon , epsilon ~ N(0, sigma)
## x consists of 10 variables, only first 5 matter

f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

## run BART
set.seed(99)
bartFit <- bartBT(x, y)

plot(bartFit)

## compare BART fit to linear matter and truth = Ey
lmFit <- lm(y ~ ., data.frame(x, y))

fitmat <- cbind(y, Ey, lmFit$fitted, bartFit$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))

## fit with missing predictor values: every split rule learns a
## direction for NAs (MIA, always on; no argument selects it)
set.seed(1)
n <- 60L
x.na <- matrix(runif(n * 2), n, 2)
x.na[1:5, 1] <- NA
y.na <- x.na[,2] + rnorm(n, 0, 0.3)
y.na[is.na(x.na[,1])] <- y.na[is.na(x.na[,1])] + 2
fit.na <- dbarts(y.na ~ x.na, control = dbartsControl(
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L))
samples.na <- fit.na$run()

## recipe: between-draws row subsetting on a latent family (e.g. principal
## stratification or mediation, where an outer Gibbs step redraws which rows
## belong to a stratum every iteration). $setActiveRows changes which rows
## are "in" the model from one run() to the next, without rebuilding the
## sampler and without losing any row's fitted value; see the "active" entry
## of \link{dbartsSampler-class} for the full semantics. Unlike zero case
## weights (gaussian, probit, ordinal and nbinom only), this also reaches
## logistic and aft samplers - and an emptied stratum is a legal all-zeros mask
## rather than a case the caller must special-case and skip.
set.seed(3)
n <- 120L
x.strat <- matrix(runif(n * 2), n, 2)
y.strat <- rbinom(n, 1L, plogis(2 * x.strat[,1] - 1))
sampler.strat <- dbarts(x.strat, y.strat, family = "probit",
                        control = dbartsControl(n.chains = 1L, n.threads = 1L,
                          n.trees = 20L, updateState = FALSE))

stratum <- rbinom(n, 1L, 0.6)  # this iteration's stratum membership
sampler.strat$setActiveRows(stratum)
draws.strat <- sampler.strat$run(20L, 20L)

## the next outer iteration redraws the stratum; even a stratum that empties
## out is accepted and runs, rather than needing to be skipped
stratum <- rep(0, n)
sampler.strat$setActiveRows(stratum)
draws.empty <- sampler.strat$run(0L, 5L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
