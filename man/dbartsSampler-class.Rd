\docType{class}
\name{dbartsSampler-class}
\title{Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler}
\alias{dbartsSampler}
\alias{dbartsSampler-class}
\alias{\S4method{run}{dbartsSampler}}
\alias{\S4method{sampleTreesFromPrior}{dbartsSampler}}
\alias{\S4method{sampleNodeParametersFromPrior}{dbartsSampler}}
\alias{\S4method{growFromRoot}{dbartsSampler}}
\alias{\S4method{copy}{dbartsSampler}}
\alias{\S4method{show}{dbartsSampler}}
\alias{\S4method{predict}{dbartsSampler}}
\alias{\S4method{setControl}{dbartsSampler}}
\alias{\S4method{setModel}{dbartsSampler}}
\alias{\S4method{setData}{dbartsSampler}}
\alias{\S4method{setResponse}{dbartsSampler}}
\alias{\S4method{setOffset}{dbartsSampler}}
\alias{\S4method{setWeights}{dbartsSampler}}
\alias{\S4method{setSigma}{dbartsSampler}}
\alias{\S4method{setPredictor}{dbartsSampler}}
\alias{\S4method{setCutPoints}{dbartsSampler}}
\alias{\S4method{setTestPredictor}{dbartsSampler}}
\alias{\S4method{setTestPredictorAndOffset}{dbartsSampler}}
\alias{\S4method{setTestOffset}{dbartsSampler}}
\alias{\S4method{printTrees}{dbartsSampler}}
\alias{\S4method{getTrees}{dbartsSampler}}
\alias{\S4method{getSigmas}{dbartsSampler}}
\alias{\S4method{getLatents}{dbartsSampler}}
\alias{\S4method{getSumsOfSquaredResiduals}{dbartsSampler}}
\alias{\S4method{installTrees}{dbartsSampler}}
\alias{\S4method{storeState}{dbartsSampler}}
\alias{\S4method{setState}{dbartsSampler}}
\alias{\S4method{plotTree}{dbartsSampler}}
\description{
  A reference class object that contains a Bayesian Additive Regression Trees sampler in such a way that it can be modified, stopped, and started all while maintaining its own state.
}
\note{
  \code{dbartsSampler} is a \link[=ReferenceClasses]{reference class}: its methods are not called as free functions but as \code{$}-dispatched calls on a sampler instance, e.g. \code{sampler$run(100, 100)} or \code{sampler$setResponse(newY)}. The S4-method notation in \sQuote{Usage} below is an artifact of how reference-class methods are documented and does not reflect the calling syntax; see \sQuote{Examples}.
}
\usage{
\S4method{run}{dbartsSampler}(
  numBurnIn, numSamples, updateState = NA, n.threads = control@n.threads
)
\S4method{sampleTreesFromPrior}{dbartsSampler}(updateState = NA)
\S4method{sampleNodeParametersFromPrior}{dbartsSampler}(updateState = NA)
\S4method{growFromRoot}{dbartsSampler}(n.sweeps = 2L, updateState = NA)
\S4method{copy}{dbartsSampler}(shallow = FALSE)
\S4method{show}{dbartsSampler}()
\S4method{predict}{dbartsSampler}(x.test, offset.test, n.threads = control@n.threads)
\S4method{setControl}{dbartsSampler}(control)
\S4method{setModel}{dbartsSampler}(model)
\S4method{setData}{dbartsSampler}(data)
\S4method{setResponse}{dbartsSampler}(y, updateScale = FALSE, updateState = NA)
\S4method{setOffset}{dbartsSampler}(offset, updateScale = FALSE, updateState = NA)
\S4method{setWeights}{dbartsSampler}(weights, updateState = NA)
\S4method{setSigma}{dbartsSampler}(sigma, updateState = NA)
\S4method{setPredictor}{dbartsSampler}(x, column, forceUpdate, updateCutPoints = FALSE, updateState = NA)
\S4method{setCutPoints}{dbartsSampler}(cuts, column, updateState = NA)
\S4method{setTestPredictor}{dbartsSampler}(x.test, column)
\S4method{setTestPredictorAndOffset}{dbartsSampler}(x.test, offset.test)
\S4method{setTestOffset}{dbartsSampler}(offset.test)
\S4method{printTrees}{dbartsSampler}(treeNums, chainNums, sampleNums)
\S4method{getTrees}{dbartsSampler}(
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL
)
\S4method{getSigmas}{dbartsSampler}(result)
\S4method{getLatents}{dbartsSampler}(result)
\S4method{getSumsOfSquaredResiduals}{dbartsSampler}(result)
\S4method{installTrees}{dbartsSampler}(donor, samples = NULL)
\S4method{storeState}{dbartsSampler}()
\S4method{setState}{dbartsSampler}(newState)
\S4method{plotTree}{dbartsSampler}(
  treeNum, chainNum, sampleNum, treePlotPars = c(
    nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
  ...
)
}
\arguments{
  \item{numBurnIn}{
    A non-negative integer determining how many iterations the sampler should skip before storing results. If missing or \code{NA}, the default is filled in from the sampler's \code{\link[=dbartsControl]{control}} object.
  }
  \item{numSamples}{
    A positive integer determining how many posterior samples should be returned. If missing or \code{NA}, the default is also filled in from the control object.
  }
  \item{updateState}{
    A logical determining if the local cache of the sampler's state should be updated after the call completes. Two conventions apply, by method: for \code{run}, \code{sampleTreesFromPrior}, and \code{sampleNodeParametersFromPrior}, \code{NA} (the default) fills in the sampler's \code{\link[=dbartsControl]{control}} object's \code{updateState}, and explicit \code{TRUE}/\code{FALSE} override it. For the mutators - \code{setData}, \code{setResponse}, \code{setOffset}, \code{setWeights}, \code{setSigma}, \code{setPredictor}, and \code{setCutPoints} - the state is stored only on explicit \code{TRUE}; \code{NA} (the default) and \code{FALSE} both store nothing, regardless of \code{control@updateState}. These are typically called once per sweep inside a larger Gibbs/MH loop (as \code{dbartsSampler} is designed for), where storing state on every mutation would be wasted work whenever the loop only reads \code{state} occasionally (or never); an unforced \code{state} promise materializes the sampler's \emph{current} state on first access regardless, so a mutate-then-first-read sequence needs no explicit store. Pass \code{TRUE} explicitly when \code{state} was already forced (read or saved) earlier and a later mutation must be reflected in the next save.
  }
  \item{shallow}{
    A logical determining if the copy should retain the underlying data of the sampler (\code{TRUE}) or have its own copies (\code{FALSE}).
  }
  \item{control}{
    An object inheriting from \code{\link{dbartsControl}}. When passed to \code{setControl}, it cannot change \code{n.trees}, \code{n.chains}, \code{useQuantiles}, or \code{rngSeed} from the values the sampler was created with, and cannot set \code{keepTrees = TRUE} without also giving \code{n.samples}; either is an error.
  }
  \item{model}{
    An object inheriting from \code{dbartsModel}. When passed to \code{setModel}, it cannot switch a DART tree prior on or off relative to the sampler's creation-time model; recreate the sampler instead.
  }
  \item{data}{
    An object inheriting from \code{\link{dbartsData}}.
  }
  \item{y}{
    A numeric response vector of length equal to that with which the sampler was created.
  }
  \item{x}{
    A numeric predictor vector of length equal to that with which the sampler was created. Can be of a distinct number of rows for \code{setTestPredictor}.
  }
  \item{x.test}{
    A new matrix, data frame, or sparse-bearing test set (the column types accepted by \code{\link{dbartsData}}), of the number of columns equal to that in the current model. A data-frame/sparse \code{x.test} is coded against the training levels; \code{setTestPredictor} and \code{setTestPredictorAndOffset} then install it as a resident, whole-object replacement of the test set (see \code{column}), while \code{predict} materializes it to a numeric matrix before evaluating.
  }
  \item{offset}{
    A numeric vector of length equal to that with which the sampler was created, or \code{NULL}. If \code{offset.test} was set from \code{offset}, will attempt to update that as well.
  }
  \item{updateScale}{
    Logical indicating whether BART's internal scale should re-anchor to the new offset (\code{setOffset}) or response (\code{setResponse}). Defaults to \code{FALSE}, locking the scale set at creation; should only be \code{TRUE} during burn-in, as re-anchoring mid-run makes the fits across iterations no longer comparable.
  }
  \item{offset.test}{
    A numeric vector of length equal to that of the test matrix, or \code{NULL}. Can be missing for \code{setTestPredictorAndOffset}.
  }
  \item{n.threads}{
    Currently has no effect: \code{run} and \code{predict} both execute
    serially regardless of the value passed here. The sampler's own thread
    count is fixed when it is created, from the \code{n.threads} of its
    \code{\link[=dbartsControl]{control}} object; this argument is reserved
    for a future per-call override.
  }
  \item{sigma}{
    Numeric vector of residual standard deviations, one for each chain.
  }
  \item{weights}{
    A numeric vector of non-negative case weights, of length equal to that with which the sampler was created. Weights of zero exclude an observation from the likelihood while keeping its fitted values. \code{setWeights} only applies to a gaussian-family sampler; calling it on a probit- or logistic-family sampler is refused, since a weighted probit has no tractable latent-variable form and logistic weights (observation counts) are fixed when the sampler is created. See \code{\link{dbarts}} for the family-specific weight rules that apply at creation time.
  }
  \item{cuts}{
    Vector of cut points for use with \code{setCutPoints}.
  }
  \item{column}{
    An integer or character string vector specifying which column/columns of
    the predictor matrix is to be replaced or for which cuts should be used.
    When updating predictors, it can be missing and the entire matrix will be
    substituted, which requires a dense design: a sparse or mixed
    dense/sparse design is mutated one column at a time, naming \code{column}
    (dense-backed and sparse-backed columns alike, and a single call may name
    both kinds). For \code{setTestPredictor}, a per-column update is refused
    when the current test set is a data frame or sparse container rather
    than a plain matrix; replace the whole test set instead.
  }
  \item{forceUpdate}{
    For \code{setPredictor}, controls how the trees respond to a new predictor
    that would leave a tree with an empty leaf-node. \code{FALSE} rolls the whole
    update back and leaves the sampler unchanged (rejection sampling); \code{TRUE}
    forces the change, collapsing any empty nodes; the string \code{"partial"}
    installs each observation's new value individually and rolls back only those
    observations whose value would empty a leaf, returning a per-observation
    logical of what was installed. \code{"partial"} requires a single \code{column}
    and cannot be combined with \code{updateCutPoints}; it also requires a
    dense-backed column, since a sparse column's storage fixes its nonzero
    pattern per cell - replace such a column whole instead. When missing, defaults to
    \code{TRUE} when the whole predictor matrix is being replaced (\code{column}
    is missing) and \code{FALSE} when a single column is being replaced.
  }
  \item{updateCutPoints}{
    For \code{setPredictor}, a logical; when \code{TRUE} the cut points (split
    candidate locations) for the replaced column(s) are recomputed from the new
    values, otherwise the existing cut points are kept. Defaults to \code{FALSE}
    and cannot be combined with \code{forceUpdate = "partial"}.
  }
  \item{treeNums}{
    An integer vector listing the indices of the trees to print or return.
  }
  \item{chainNums}{
    An integer vector listing the chains to return from \code{getTrees}. When
    missing, all chains are returned.
  }
  \item{sampleNums}{
    An integer vector listing the saved samples to return from \code{getTrees}.
    Applies only when \code{keepTrees} is \code{TRUE} and \code{current} is
    \code{FALSE}; otherwise the live working trees have no sample dimension and
    it is ignored.
  }
  \item{current}{
    For \code{getTrees}, a logical; when \code{TRUE} the live working trees (the
    current sampler position) are returned even if \code{keepTrees} is
    \code{TRUE}, rather than the saved samples.
  }
  \item{newdata}{
    For \code{getTrees}, an optional matrix of predictors (in the form accepted
    by \code{predict}). When supplied, its observations are routed through each
    tree so the \code{n} column counts them instead of the training data; the
    tree structure, splits, and leaf values are unchanged.
  }
  \item{result}{
    For \code{getLatents}, an optional pre-allocated numeric vector (or,
    with multiple chains, matrix) that it fills in place and returns instead
    of allocating a new one; its length must equal the number of
    observations times the number of chains. Omit it to let \code{getLatents}
    allocate. Accepted but not used by \code{getSigmas} and
    \code{getSumsOfSquaredResiduals}.
  }
  \item{treeNum}{
    An integer listing the indices of the tree to plot.
  }
  \item{chainNum}{
    For \code{plotTree}, an integer giving the chain to plot from. Required
    when the sampler has more than one chain; defaults to \code{1} when
    there is only one.
  }
  \item{sampleNum}{
    For \code{plotTree}, an integer giving the saved sample to plot.
    Applies only when \code{keepTrees} is \code{TRUE}, in which case it
    defaults to the most recently drawn sample; ignored (with the live
    working trees plotted instead) when \code{keepTrees} is \code{FALSE}.
  }
  \item{treePlotPars}{
    A named numeric vector containing the quantities \code{nodeHeight}, \code{nodeWidth}, and \code{nodeGap}, all of which control aspects of the resulting plot.
  }
  \item{donor}{
    For \code{installTrees}, the fit whose forests seed this sampler: another \code{dbartsSampler}, or a \code{bart} object fit with \code{keepSampler = TRUE}. The donor must share this sampler's predictors, tree count, and DART setting; a different cut grid is allowed, its splits remapped onto this sampler's grid.
  }
  \item{samples}{
    For \code{installTrees}, an optional integer vector with one entry per chain, each a 1-based index into the donor's pool of samples (its saved trees when it kept them, else its final trees per chain). \code{NULL} spreads the chains evenly across the pool.
  }
  \item{n.sweeps}{
    For \code{growFromRoot}, a single positive integer giving the number of grow-from-root sweeps run in place (default \code{2L}). Each sweep rebuilds every tree in the forest against the current residual, so a small handful reaches a good fit.
  }
  \item{newState}{
    For \code{setState}, a state object previously produced by this sampler (its \code{state} field, or the return of \code{storeState}) over the same model. Must inherit from \code{bartcoreState}.
  }
  \item{...}{
    Extra arguments to \code{\link{plot}}.
  }
}
\section{Fields}{
  Read with \code{$}, as the methods are called. Treat them as read-only: assigning to a field changes only the R-side copy and leaves the underlying engine untouched, so the sampler's behavior does not follow. Route changes through the \code{set*} methods instead.

  \describe{
    \item{\code{control}}{
      The sampler's \code{\link{dbartsControl}}, carrying \code{n.chains}, \code{n.threads}, \code{n.trees}, \code{keepTrees}, \code{updateState}, and the rest. Replace with \code{setControl}.
    }
    \item{\code{model}}{
      A \code{dbartsModel} holding the parsed tree, node, and residual priors, together with any monotonicity, interaction, and block constraints. Replace with \code{setModel}.
    }
    \item{\code{data}}{
      A \code{\link{dbartsData}} holding the current predictors, response, test data, weights, and offset. Its \code{x} slot is stored columnar for data-frame input and so need not be a plain matrix; use \code{\link[=extract.dbartsSampler]{extract}(sampler, "predictors")} for the numeric predictor matrix, which is what tree replay and \code{getTrees} count against. Replace with \code{setData}, or a piece at a time with the other \code{set*} methods.
    }
    \item{\code{state}}{
      The cached, serializable engine state, or \code{NULL} until one is materialized. Reading it forces the sampler's \emph{current} state; \code{storeState} refreshes it on demand and \code{updateState} governs when the methods do so themselves. It is the only field \code{\link{save}} needs, and restoring one requires \code{setState} - see \sQuote{Saving}.
    }
    \item{\code{pointer}}{
      The external pointer to the C++ engine. Internal; never manipulate it directly.
    }
  }
}
\details{
  A \code{dbartsSampler} is a mutable object which contains information pertaining to fitting a Bayesian additive regression tree model. The sampler is first created and then, in a separate instruction, run or modified. In this way, MCMC samplers can be constructed with BART components filling arbitrary roles.

  \subsection{Saving}{
    \code{\link{save}}-ing and \code{\link{load}}ing a \code{dbarts} sampler for future use requires that R's serialization mechanism be able to access the state of the sampler which, for memory purposes, is only made available to R on request. To make it available, call \code{storeState()} on the sampler before saving, e.g. for the object \code{sampler}, execute \code{sampler$storeState()}; this captures the sampler's current state into the serializable \code{state} field regardless of the \code{updateState} setting.

    The state object is opaque and engine-specific, and carries a format version and the package version that wrote it. A state loads only within the format version of the \code{dbarts} that wrote it; loading it with an incompatible version refuses cleanly, naming both versions, rather than risk a silent misread. There is no cross-version migration: re-fit the model, or restore the state with the \code{dbarts} release that wrote it.

    To restore a saved state into a sampler, call \code{setState(newState)}: it validates that \code{newState} inherits from \code{bartcoreState}, re-creates the underlying engine if needed, pushes the state into it, and caches it on the \code{state} field. Assigning the field directly (\code{sampler$state <- newState}) does \emph{not} restore the sampler - it only overwrites the R-side cache, leaving the engine untouched, so the next run continues from the engine's own state rather than the assigned one. Always route a restore through \code{setState}.
  }

  \subsection{Mutation cost}{
    Each accepted predictor mutation does two things: the engine updates its internal representation, and the R layer collects the accepted values into the sampler's data object so that replay, saving, and re-quantization see current data. The collection is by reference for data-frame-built samplers (only the affected column changes hands) but requires copying the full predictor matrix when the sampler was built from a matrix, so tight loops that mutate predictors every iteration are better served by data-frame input. The remaining per-call overhead is a few tens of microseconds of R method dispatch and bookkeeping - the price of R-level consistency. Clients for which that matters can drive the sampler through the C interface (see the \code{dbarts.h} header installed with the package), which invokes the engine directly and performs no R-side collection; such clients supply their own current predictor matrix when replaying saved trees.
  }

  \subsection{Warm starts}{
    \code{installTrees} seeds the sampler's forests from a \code{donor} instead of drawing trees from the prior, for scaling to more chains or embedding a fit in a larger sampler. Only the donor's trees, \code{sigma}, and \code{k} transfer; each chain keeps its own random-number stream and redraws everything else, so several chains seeded from one donor stay overdispersed. A donor with a different tree count or DART setting is refused rather than silently reshaped; a donor fit on a different cut grid is instead remapped onto this sampler's grid, collapsing any splits the grid starves, the same way a data replacement remaps existing splits. A warm start biases the early draws toward the donor, so it shortens burn-in rather than removing it; keep drawing a non-zero number of burn-in samples before treating the chain as converged.

    \code{growFromRoot} instead builds the sampler's initial forest by XBART-style root-down stochastic tree construction (He, Yalov and Hahn 2019): each of \code{n.sweeps} sweeps rebuilds every tree from the root, sampling each split from the integrated-likelihood weight of every candidate cut under the tree prior. This reaches a good fit in far fewer sweeps than the exact sampler, so it is a fast starting point rather than a posterior sampler - the exact MCMC sweeps own stationarity once \code{run} begins, and the posterior is unchanged. It is available for the constant-leaf model only; calling it on a \code{linear} or \code{gp} node prior is an error, not a silent fall-back, so initialize those forests with \code{sampleTreesFromPrior} instead. As with \code{installTrees}, the grown forest biases the early draws toward its fit, so shorten burn-in rather than skipping it. Each chain grows on its own random-number stream, so the result is independent of \code{n.threads}. A composable cross-sampler workflow follows for free: \code{donor$growFromRoot(k)} then \code{target$installTrees(donor)}.
  }
}
\value{
  For \code{run}, a named-list with contents \code{sigma}, \code{train}, \code{test}, and \code{varcount} (plus \code{k}, \code{varprobs}, \code{tau}, and \code{ranef} when applicable). \code{train} is an array of dimension n.obs x n.samples x n.chains, and likewise \code{test} (or \code{NULL} if the sampler has no test data) and \code{varcount} are n.predictors x n.samples x n.chains; \code{sigma} is n.samples x n.chains. When \code{n.chains} is \code{1} the trailing chain dimension is dropped, so \code{train} is a plain n.obs x n.samples matrix and \code{sigma} a plain vector of length n.samples. A run can be interrupted with \code{Ctrl-C}: it stops between iterations - joining any worker threads first - and signals an error, returning no samples from the interrupted run. The sampler's chains are left at the iteration they reached, which is a valid state to run again from.

  For \code{setPredictor}, \code{TRUE}/\code{FALSE} depending on whether or not the operation was successful. The operation can fail if the new predictor results in a tree with an empty leaf-node. If only single columns were replaced, the update is rolled back so that the sampler remains in a valid state. When \code{forceUpdate} is \code{"partial"}, instead returns a logical vector of length equal to the number of observations, \code{TRUE} where that observation's new value was installed and \code{FALSE} where it was rolled back to its previous value to keep every tree valid.

  \code{predict} keeps the current test matrix in place and uses the current set of tree splits. This function has two use cases. The first is when \code{keepTrees} of \code{\link{dbartsControl}} is \code{TRUE}, in which case the sampler should be run to completion and the function can be used to interrogate the existing fit. When \code{keepTrees} is \code{FALSE}, the function can be used to obtain the likelihood as part of a proposed new set of covariates in a Metropolis-Hastings step in a full-Bayes sampler. This would typically be followed by a call to \code{setPredictor} if the step is accepted.

  For \code{getTrees}, a \code{data.frame} with one row per tree node in depth-first, left-hand-side pre-order, with columns \code{chain}, \code{sample} (present only for saved samples), \code{tree}, \code{n} (the number of observations in the node), \code{var} (the splitting variable, or -1 at a leaf), and \code{value} (the split value, or the leaf prediction). An ordinal rule's value is its cut point and observations with values less than or equal to it go left; a categorical rule carries no data value (its \code{value} is \code{NA}) and its split is reported in the \code{directions} column instead. When the sampler has any categorical predictors the result gains a \code{directions} column decoding each categorical rule into one \code{"L"}/\code{"R"} character per level, in level order (level \code{k} goes right when its character is \code{"R"}); ordinal rules and leaves are \code{NA}. When any predictor contains missing values the result gains a \code{missing} column giving the branch (\code{"L"}/\code{"R"}) each rule sends missing values down; rules on complete columns and leaves are \code{NA}. Under a \code{linear} node prior each leaf's \code{value} is its intercept and the result gains one \code{beta.<column>} column per designated covariate holding that leaf's slope on the internal standardized scale; internal nodes are \code{NA}.

  For \code{getSigmas}, a numeric vector of length equal to the number of chains, giving each chain's current residual standard deviation on the original response scale.

  For \code{getSumsOfSquaredResiduals}, a numeric vector of length equal to the number of chains, giving each chain's residual sum of squares \eqn{\sum (y - \hat{y})^2}{sum((y - yhat)^2)} on the original response scale; a binary-response sampler reports on the latent scale instead.

  For \code{getLatents}, \code{NULL} when the model has no latent-variable representation (e.g. a gaussian response); otherwise the sampler's current latent values - a plain vector of length equal to the number of observations when there is a single chain, or an observations-by-chains matrix otherwise - written into \code{result} when one was supplied.

  For \code{storeState}, \code{NULL} invisibly; it is called for its side effect of capturing the sampler's current engine state into the serializable \code{state} field, which \code{\link{save}} then writes out. See \sQuote{Saving}.
}
\seealso{
  \code{\link{updatePredictorPerObservationJointly}} for applying a shared per-observation predictor update across several samplers at once.

  \code{\link{samplePriorPredictive}} for repeated \code{sampleTreesFromPrior}/\code{sampleNodeParametersFromPrior}/\code{predict} draws on a private sampler, for calibrating priors before fitting.
}
\examples{
# the embedded-Gibbs pattern: BART as a conditional model inside a larger
# sampler, alternating its own draws with updates to its response between
# calls to run()
n <- 100
x <- matrix(runif(n * 2), n, 2)
y <- x[, 1] - x[, 2] + rnorm(n, 0, 0.1)

sampler <- dbarts(
  y ~ x,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 0L,
    n.samples = 1L,
    updateState = TRUE
  )
)

## first draw: response as given at creation
samples <- sampler$run()
str(samples$train) # a plain n.obs x n.samples matrix (n.chains == 1)

## an outer Gibbs step changes the response (e.g. after updating some
## other part of a joint model); the sampler picks the new target up on
## the next run() without being recreated
newY <- y + rnorm(n, 0, 0.05)
sampler$setResponse(newY)
samples <- sampler$run()
str(samples$train)
}
