\name{dbarts}
\title{Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbarts}
\description{
  Creates a sampler object for a given problem which fits a Bayesian Additive Regression Trees model. Internally stores state in such a way as to be mutable.
}
\usage{
dbarts(
    formula, data, test, subset, weights, offset, offset.test = offset,
    verbose = FALSE, n.samples = 800L,
    tree.prior = cgm, node.prior = normal, resid.prior = chisq,
    resid.dist = gaussian,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL,
    variance = NULL, n.trees.variance = 40L,
    power.variance = NULL, base.variance = NULL,
    control = dbarts::dbartsControl(), sigma = NA_real_, seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic", "aft", "ordinal",
               "nbinom", "hazard", "hazard.probit", "hazard.logistic"),
    missing = c("incorporate", "error"), dispersion = NA_real_,
    breaks = NULL, max.rows = 1e7)}
\arguments{
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bart}} matrix \code{x.train}, including a sparse \code{Matrix::dgCMatrix}: its columns enter as ordinal predictors, sufficiently sparse columns are stored in a compact rank-bitmap layout instead of being expanded, and the predictor-mutation surface (\code{setPredictor} and relatives, \code{setData}) is fixed at creation. Sparse inputs are not supported by \code{\link{rbart_vi}} or the \code{node.prior = linear()} and \code{gp()} leaf models (a sparse-backed test column cannot serve a designated leaf covariate either). A sparse or data-frame test set stays resident (rank-bitmap or dense per column, the same rule as training) rather than densifying at ingestion; only \code{predict} and \code{getTrees(newdata = )} materialize it to a dense matrix, since those entry points evaluate already-built trees on resolved values. A data frame may mix ordinary columns with \code{Matrix::sparseVector} or \code{dgCMatrix} columns (assign them into the frame; they do not survive \code{data.frame(...)} or \code{I()}): sparse columns behave as in an all-sparse design, while dense columns keep categorical splits and linear-leaf designation. A \code{\link{sparseFactor}} column is accepted here (the x/y, data-frame path) exactly as the other sparse columns are; a bare \code{formula} object cannot carry it through \code{model.frame} and refuses it, with a message to use the x/y interface instead. See \code{\link{dbartsData}} for the full data-frame-to-predictor mapping.
  }
  \item{data}{
    An optional data frame, list, or environment containing predictors to be used with the model. For backwards compatibility, can also be the \code{\link{bart}} vector \code{y.train}.
  }
  \item{test}{
    An optional matrix or data frame with the same number of predictors as \code{data}, or \code{formula} in backwards compatibility mode. If column names are present, a matching algorithm is used.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
    An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. Binary responses differ: a \code{"probit"} model does not support weights (a weighted probit has no tractable latent-variable form), except that weights identically 1 are treated as absent; a \code{"logistic"} model treats them as observation counts and so requires positive integers (its Polya-Gamma latent for a count \eqn{w} is a sum of \eqn{w} unit draws).
  }
  \item{offset}{
    An optional vector specifying an offset from 0 for the relationship between the underlying function, \eqn{f(x)}, and the response \eqn{y}. Useful for both response families, though the two apply it differently. For a gaussian response, \eqn{y = f(x) + \mathrm{offset} + \epsilon}{y = f(x) + offset + epsilon}: the offset is netted out of \eqn{y} before BART's internal range-scaling (see \sQuote{Details}), so it acts as a fixed component of the mean rather than a free parameter. For binary responses it enters the link directly: a \code{"probit"} fit assumes \eqn{P(Y = 1 \mid X = x) = \Phi(f(x) + \mathrm{offset})}{P(Y = 1 | X = x) = \Phi(f(x) + offset)}, where \eqn{\Phi} is the standard normal cumulative distribution function, and a \code{"logistic"} fit the analogous model on the logistic link.
  }
  \item{offset.test}{
    The equivalent of \code{offset} for test observations. Will attempt to use \code{offset} when applicable.
  }
  \item{verbose}{
    A logical determining if additional output is printed to the console. See \code{\link{dbartsControl}}.
  }
  \item{n.samples}{
    A positive integer setting the default number of posterior samples to be returned for each run of the sampler. Can be overriden at run-time. See \code{\link{dbartsControl}}.
  }
  \item{tree.prior}{
    An expression of the form \code{cgm} or \code{cgm(power, base, split.probs)} setting the tree prior used in fitting, or a prior object built with \code{\link{dbartsPriors}}.
  }
  \item{node.prior}{
    An expression of the form \code{normal} or \code{normal(k)} that sets the prior used on the averages within nodes, or a prior object built with \code{\link{dbartsPriors}}. \code{linear(columns, k)} instead fits each leaf with an intercept plus a linear term in the designated continuous predictor columns (character names or numeric indices into the model matrix), standardized internally, all coefficients sharing the \code{normal(k)} prior; factor columns cannot be designated. \code{gp(columns, k, lengthscale, max.leaf.size)} fits each leaf with a smooth Gaussian-process function of the designated columns under a squared-exponential kernel; leaves larger than \code{max.leaf.size} fall back to constant fits. \code{\link{xbart}} accepts the same specifications through its own \code{node.prior} argument. See \dQuote{Response scaling} below for how \code{k} interacts with the response's internal scaling.
  }
  \item{resid.prior}{
    An expression of the form \code{chisq} or \code{chisq(df, quant)} that sets the prior used on the residual/error variance, or a prior object built with \code{\link{dbartsPriors}}.
  }
  \item{resid.dist}{
    The residual error \emph{law} for a continuous (gaussian) response, orthogonal to \code{resid.prior} (which is the prior on the error variance). The default \code{gaussian()} gives the usual normal errors. \code{student(df)} instead fits outlier-robust Student-t errors by the classic Gaussian scale-mixture augmentation: each observation carries a latent precision \eqn{\lambda_i \sim \mathrm{Gamma}(\nu/2, \nu/2)}{lambda_i ~ Gamma(nu/2, nu/2)} so that \eqn{\sqrt{w_i}\,\epsilon_i / \sigma \sim t_\nu}{sqrt(w_i) eps_i / sigma ~ t_nu}, and a gross outlier draws a small \eqn{\lambda_i}, downweighting its leverage on the fit and on \eqn{\sigma}. \code{student(df = nu)} fixes the degrees of freedom at \code{nu} (a positive number; \code{nu = 4} is a conventional robust default), while \code{student()} (equivalently \code{student(df = NULL)}) estimates them on a capped grid under a proper tail-bounding prior. Only a continuous gaussian response carries this: \code{student()} with a \code{"probit"}, \code{"logistic"}, or \code{"aft"} family is an error, as those have their own fixed latent scale. Two caveats: (a) with a flexible tree-forest mean, tail inference is partially confounded with fit flexibility, so an estimated \eqn{\nu} should be posterior-checked rather than read as a pure noise property; (b) \eqn{\sigma} is now the \emph{conditional} scale, so the marginal residual variance is \eqn{\sigma^2\,\nu/(\nu-2)}{sigma^2 nu / (nu - 2)}, not \eqn{\sigma^2}{sigma^2}. See \dQuote{Response scaling} below - robust errors mitigate an outlier's leverage but not the range compression that outlier causes.
  }
  \item{proposal.probs}{
    Named numeric vector or \code{NULL}, optionally specifying the proposal rules and their probabilities. Elements should be \code{"birth_death"}, \code{"change"}, and \code{"swap"} to control tree change proposals, and \code{"birth"} to give the relative frequency of birth/death in the \code{"birth_death"} step.
  }
  \item{monotone}{
    Optional per-predictor monotonicity constraints (monotone BART; Chipman, George, McCulloch, and Shively 2022). The forest is constrained so that the fitted function is monotone increasing or decreasing in each named predictor. Supply a named vector selecting predictors by model-matrix column name, each element one of \code{"+"}/\code{"increasing"}/\code{1} (increasing) or \code{"-"}/\code{"decreasing"}/\code{-1} (decreasing); an unnamed integer or character vector of length equal to the number of columns assigns directions positionally, with \code{0} for unconstrained. Only numeric and ordered columns are eligible - a direction on a categorical (unordered factor) predictor is an error. A constraint forces birth/death-only tree proposals (an explicit non-default \code{proposal.probs} is then an error) and a fixed \code{k = 2} (an explicit \code{k} hyperprior is an error); linear and Gaussian-process leaves are not supported under the constraint. \code{NULL} (the default) or an all-zero vector fits the ordinary unconstrained model.
  }
  \item{variance}{
    Optional heteroscedastic variance forest (HBART; Pratola, Chipman, George, and McCulloch 2020). When supplied, a second ensemble models the residual variance surface \eqn{s^2(x)} as a product of scaled-inverse-chi-squared (multiplicative) leaves, so \eqn{y_i = f(x_i) + s(x_i)\epsilon_i}, coupled to the mean forest through the per-observation precision. \code{variance} selects which predictors drive the variance: a one-sided formula or character/integer column selector (\code{~ x1 + x2}), \code{TRUE} for every predictor, or \code{NULL} (the default) for the ordinary homoscedastic fit. Gaussian responses and constant leaves only; monotone constraints and the latent families are not supported. The per-tree leaf prior is calibrated from the residual (\code{resid.prior}) hyperparameters so that a constant variance surface reproduces the homoscedastic \code{sigma} posterior. The fit gains posterior draws \code{s.train}/\code{s.test} of \eqn{s(x)}, and \code{predict} attaches an \code{"s"} attribute with \eqn{s(x)} at new predictors (requires \code{keepTrees}).
  }
  \item{n.trees.variance, power.variance, base.variance}{
    The number of variance trees (default \code{40}) and their tree-structure prior. \code{power.variance} and \code{base.variance} default to the mean forest's \code{tree.prior} values. Read only when \code{variance} is supplied.
  }
  \item{control}{
    An object inheriting from \code{dbartsControl}, created by the \code{\link{dbartsControl}} function.} \item{sigma}{A positive numeric estimate of the residual standard deviation. If \code{NA}, a linear model is used with all of the predictors to obtain one. Same concept as \code{sigest} in \code{\link{bart}}/\code{bart2} and \code{\link{rbart_vi}}; \code{\link{xbart}} uses this function's \code{sigma} name directly.
  }
  \item{seed}{
    Optional integer seed for the random number generator, a convenience mirror of \code{dbartsControl(rngSeed = )}. When not \code{NA} it overrides the seed in \code{control}; the fitting-function wrappers (\code{\link{bart2}}, \code{\link{xbart}}) accept the same argument.
  }
  \item{factors}{
    How factor columns in a data frame enter the model. The default \code{"categorical"} keeps each unordered factor as a single predictor whose splits send a subset of its levels (at most 65535) down each branch, and codes ordered factors as ordinal; level tables are retained so that test data are coded identically. \code{"indicators"} expands each factor into binary indicator columns, as previous versions always did and as \code{\link{bart}} still does.
  }
  \item{family}{
    The response model. \code{"auto"} fits gaussian models to continuous responses and probit models to those coded 0/1, as always; a two-level factor, logical, or two-level character response is also detected and fit as probit, reporting the choice in a one-line message, while a factor (or character) response with three or more levels is an error directing to \code{\link{bart2}}'s \code{family = "multinomial"} (which \code{dbarts} does not fit). An explicit family that a factor response cannot support (e.g. \code{"gaussian"}) is also an error rather than a silent fit of the integer level codes. \code{"gaussian"} forces a continuous fit even for a 0/1 numeric response; \code{"probit"} and \code{"logistic"} require a 0/1 response and fit latent-variable models, with fits and predictions on the latent scale. \code{"logistic"} uses Polya-Gamma augmentation. \code{"aft"} fits an accelerated failure time (log-normal) survival model: the response is a \code{Surv} object (from the \pkg{survival} package) or a two-column \code{(time, status)} matrix or data frame (status 1 an event, 0 right-censored; logical status also works, factor status does not), and \code{\link{survivalProbabilities}} produces survival-probability draws. Fits and predictions are on the log-time scale (\eqn{E[\log T \mid x]}), not the time scale. A \code{Surv} response selects \code{"aft"} automatically when \code{family} is \code{"auto"}; an explicitly conflicting family (e.g. \code{"gaussian"}) with a \code{Surv} response is an error. A \code{Surv}-like object carrying no \code{type} attribute is treated as right-censored; \code{type}s other than \code{"right"} are rejected. Survival fits currently use the matrix (\code{x.train}/\code{y.train}) interface and do not support \code{subset} or case weights. See \code{weights} for how the families differ in their support for weights. \code{"ordinal"} fits an ordered categorical response by a cumulative probit: a latent \eqn{z = f(x) + \epsilon}{z = f(x) + eps}, \eqn{\epsilon \sim N(0, 1)}{eps ~ N(0, 1)}, is cut at ordered thresholds \eqn{\gamma_1 = 0 < \gamma_2 < \ldots < \gamma_{K-1}}{gamma_1 = 0 < gamma_2 < ... < gamma_{K-1}} into the \eqn{K} ordered categories, the free cutpoints sampled by a marginal Metropolis update with the latents integrated out. The response should be an ordered factor (\code{\link{is.ordered}}), whose level order defines the category order; under \code{family = "auto"} an ordered factor with three or more levels is detected and fit as ordinal, reporting the choice in a one-line message (a two-level ordered factor is binary and resolves to probit). \code{family = "ordinal"} given explicitly also accepts an unordered factor or character response - with a message that the category order is taken from the level order - and a numeric response, whose ordered levels are \code{sort(unique(y))}. Like probit, the latent scale is fixed at 1, fits are on the latent scale, and weights are not supported. \code{bart2} reports ordinal fits as \eqn{n \times K}{n x K} category probabilities; see \code{\link{bart2}}. \code{"nbinom"} fits a non-negative integer (count) response by a negative-binomial model with the Polya-Gamma augmentation: the forest fits a log-odds latent \eqn{\psi = f(x) + o}{psi = f(x) + o} and \eqn{y \sim \mathrm{NB}(r, \mathrm{plogis}(\psi))}{y ~ NB(r, plogis(psi))}, with mean \eqn{E[y \mid x] = r e^{\psi}}{E[y | x] = r exp(psi)}, so the offset enters multiplicatively as a log-exposure. The response must be a non-negative integer, and \code{"nbinom"} is never inferred - a count carries no unambiguous class, so it must be requested explicitly. The dispersion \code{r} is estimated by default (see \code{dispersion}); like probit, the latent scale is fixed at 1, fits are on the latent (log-odds) scale, and weights are not supported (exposure belongs in the offset). \code{bart2} reports mean counts; see \code{\link{bart2}}. \code{"hazard"} and \code{"hazard.logistic"} fit a discrete-time survival hazard model by person-period expansion (\code{"hazard.probit"} is an accepted alias for \code{"hazard"}). The response is the same \code{Surv} object or two-column \code{(time, status)} the \code{"aft"} family takes, but each subject observed to time \eqn{t_i} is expanded into its at-risk rows (one per period \eqn{k = 1, \ldots, t_i}), each carrying the subject's covariates, an ordinal period column appended last, and the binary indicator \eqn{y_{ik} = \mathrm{status}_i \cdot 1\{k = t_i\}}; the discrete hazard is then \eqn{h(k \mid x) = g(f(x, k) + o)} for the chosen binary link \eqn{g} (probit for \code{"hazard"}, logistic for \code{"hazard.logistic"}), and the fit IS an ordinary \code{"probit"}/\code{"logistic"} fit on the expanded rows - it adds no engine code. The time grid is the sorted distinct observed times by default; \code{breaks} coarsens it. The offset is on the link scale and replicates per subject; weights replicate and follow the chosen binary family's policy. A \code{Surv} response requires the family to be requested explicitly (\code{"auto"} selects \code{"aft"}). Survival curves come from \code{\link{survivalProbabilities}}, which requires \code{keepTrees}; \code{bart2} reports the fit with its \code{$periods} grid, and its \code{$family} records the binary link. Like \code{"aft"}, hazard fits use the matrix interface and do not support \code{subset} or a \code{test} set (expand test subjects with \code{survivalProbabilities(..., newdata = )}).
  }
  \item{dispersion}{
    The negative-binomial dispersion \eqn{r} (family \code{"nbinom"} only; ignored otherwise). \code{NA} (the default) estimates \eqn{r} on a capped positive-integer grid under a renormalized \eqn{\mathrm{gamma}(2, 0.1)}{gamma(2, 0.1)} prior. A supplied value fixes \eqn{r} and must be a positive integer: v1 ships the exact integer envelope, so a real fixed dispersion is refused. Larger \eqn{r} approaches the Poisson limit (variance \eqn{\to}{->} mean); smaller \eqn{r} is more overdispersed.
  }
  \item{missing}{
    How missing values in the predictors enter the model. The default \code{"incorporate"} keeps them: every split rule learns a direction for missing values on its variable, so an observation whose split value is \code{NA} follows that rule's chosen branch (\dQuote{Missingness Incorporated in Attributes}, Twala et al. 2008), in training and test data alike. \code{"error"} rejects predictors containing \code{NA}. The response, \code{weights}, and \code{offset} must always be complete; note that previous versions silently dropped incomplete rows for formula inputs.
  }
  \item{breaks}{
    The discrete-time grid for a \code{"hazard"} fit (ignored otherwise). \code{NULL} (the default) uses the sorted distinct observed times, one period per distinct time (the \pkg{BART} \code{surv.bart} convention). A single positive integer bins the times at the \eqn{(1{:}K)/K} quantiles, giving \eqn{K} periods. A numeric vector of length two or more gives explicit strictly-increasing interval boundaries \eqn{b_0 < \ldots < b_K}, defining right-closed periods \eqn{(b_{k-1}, b_k]}; every time must lie in \eqn{(b_1, b_K]}.
  }
  \item{max.rows}{
    A guard on the person-period expansion for a \code{"hazard"} fit (ignored otherwise). If the expanded design would exceed \code{max.rows} rows (\eqn{\sum_i t_i}), the fit is refused with a message naming the coarsening levers. The default \eqn{10^7} catches an over-fine grid on heavily continuous times; coarsen with \code{breaks} or raise \code{max.rows}.
  }
}
\details{
  \dQuote{Discrete sampler} refers to that \code{dbarts} is implemented using \link{ReferenceClasses}, so that there exists a mutable object constructed in C++ that is largely obscured from R. The \code{dbarts} function is the primary way of creating a \code{\link[=dbartsSampler-class]{dbartsSampler}}, for which a variety of methods exist.

  \subsection{Data frame ingestion}{
    Data frame input is stored columnar: no \eqn{n \times p}{n x p} predictor matrix is retained. Code that previously reached into the sampler's stored data expecting a plain matrix should call \code{\link{extract}(sampler, "predictors")} - or \code{as.matrix} on it - to obtain the numeric predictor matrix, factor columns as their integer codes. See \code{\link{dbartsData}} for the column-type mapping (numeric, factor, sparse, and \code{\link{sparseFactor}}).
  }

  \subsection{Test offset synchronization}{
    When \code{offset.test} is left at its default of tracking \code{offset}, the sampler links the two: the sampler's \code{setOffset} method re-derives \code{offset.test} from each new \code{offset} it is given. Calling the sampler's \code{setTestOffset} or \code{setTestPredictorAndOffset} methods breaks this link -- \code{offset.test} is set independently from then on, and the link never re-forms, even if the two are later set to equal values.
  }

  \subsection{Response scaling}{
    Continuous responses are range-scaled internally: \code{y} (net of \code{offset}) is mapped to \eqn{[-0.5, 0.5]} by its observed minimum and maximum, the convention of the entire BART software lineage (BayesTree, BART, bartMachine), which is what lets \code{k} (see \code{node.prior} above) and its defaults transfer across packages and papers. The known caveat is outlier sensitivity: extreme \code{y} values stretch the range and compress the effective leaf prior for everything else; the published workaround is to log-transform or winsorize such values before fitting, or to let \code{k} adapt via the \code{chi} hyperprior (\code{\link{dbartsPriors}}). The scale is fixed when the sampler is created; the sampler's \code{setResponse} and \code{setOffset} methods only re-anchor it when called with \code{updateScale = TRUE}, intended for burn-in only, since re-anchoring mid-run makes fits across iterations no longer comparable (see \link{dbartsSampler-class}).
  }
}
\value{
  A reference object of \code{\link[=dbartsSampler-class]{dbartsSampler}}.
}
\seealso{
  \code{\link{dbartsData}}, \code{\link{extract.dbartsSampler}}, \code{\link{sparseFactor}}
}
\references{
Twala, B.E.T.H., Jones, M.C., and Hand, D.J. (2008)
   Good methods for coping with missing data in decision trees.
        \emph{Pattern Recognition Letters}, \bold{29}(7), 950--956.
}
