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
\alias{\S4method{predictForests}{dbartsSampler}}
\alias{\S4method{setControl}{dbartsSampler}}
\alias{\S4method{setModel}{dbartsSampler}}
\alias{\S4method{setData}{dbartsSampler}}
\alias{\S4method{setResponse}{dbartsSampler}}
\alias{\S4method{setOffset}{dbartsSampler}}
\alias{\S4method{setWeights}{dbartsSampler}}
\alias{\S4method{setCounts}{dbartsSampler}}
\alias{\S4method{setCategoryOffset}{dbartsSampler}}
\alias{\S4method{setCategoryTestOffset}{dbartsSampler}}
\alias{\S4method{setActiveRows}{dbartsSampler}}
\alias{\S4method{setForestWeights}{dbartsSampler}}
\alias{\S4method{setForestBasis}{dbartsSampler}}
\alias{\S4method{setSigma}{dbartsSampler}}
\alias{\S4method{setPredictor}{dbartsSampler}}
\alias{\S4method{setCutPoints}{dbartsSampler}}
\alias{\S4method{setTestPredictor}{dbartsSampler}}
\alias{\S4method{setTestPredictorAndOffset}{dbartsSampler}}
\alias{\S4method{setTestOffset}{dbartsSampler}}
\alias{\S4method{printTrees}{dbartsSampler}}
\alias{\S4method{getTrees}{dbartsSampler}}
\alias{\S4method{getSigmas}{dbartsSampler}}
\alias{\S4method{getDispersion}{dbartsSampler}}
\alias{\S4method{getLatents}{dbartsSampler}}
\alias{\S4method{getSumsOfSquaredResiduals}{dbartsSampler}}
\alias{\S4method{getFitsWithoutOffset}{dbartsSampler}}
\alias{\S4method{getForestFits}{dbartsSampler}}
\alias{\S4method{getForestAmplitudes}{dbartsSampler}}
\alias{\S4method{getForestVariableCounts}{dbartsSampler}}
\alias{\S4method{getCalibration}{dbartsSampler}}
\alias{\S4method{setCalibration}{dbartsSampler}}
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
\S4method{predictForests}{dbartsSampler}(
  x.test, offset.test, n.threads = control@n.threads
)
\S4method{setControl}{dbartsSampler}(newControl)
\S4method{setModel}{dbartsSampler}(newModel)
\S4method{setData}{dbartsSampler}(newData, updateState = NA)
\S4method{setResponse}{dbartsSampler}(y, updateScale = FALSE, updateState = NA)
\S4method{setOffset}{dbartsSampler}(offset, updateScale = FALSE, updateState = NA)
\S4method{setWeights}{dbartsSampler}(weights, updateState = NA)
\S4method{setCounts}{dbartsSampler}(counts, updateState = NA)
\S4method{setCategoryOffset}{dbartsSampler}(offset, updateState = NA)
\S4method{setCategoryTestOffset}{dbartsSampler}(offset.test, updateState = NA)
\S4method{setActiveRows}{dbartsSampler}(active, updateState = NA)
\S4method{setForestWeights}{dbartsSampler}(forest, weights, updateState = NA)
\S4method{setForestBasis}{dbartsSampler}(forest, basis, updateState = NA)
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
\S4method{getDispersion}{dbartsSampler}()
\S4method{getLatents}{dbartsSampler}(result)
\S4method{getSumsOfSquaredResiduals}{dbartsSampler}(result)
\S4method{getFitsWithoutOffset}{dbartsSampler}()
\S4method{getForestFits}{dbartsSampler}(forest = NULL)
\S4method{getForestAmplitudes}{dbartsSampler}(forest = NULL)
\S4method{getForestVariableCounts}{dbartsSampler}(forest = NULL)
\S4method{getCalibration}{dbartsSampler}(forest = NULL)
\S4method{setCalibration}{dbartsSampler}(
  prior.scale, prior.sd, prior.mean, forest = 1L, updateState = NA
)
\S4method{installTrees}{dbartsSampler}(donor, samples = NULL)
\S4method{storeState}{dbartsSampler}(ptr = getPointer())
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
    A logical determining if the local cache of the sampler's state should be updated after the call completes. Two conventions apply, by method: for \code{run}, \code{sampleTreesFromPrior}, and \code{sampleNodeParametersFromPrior}, \code{NA} (the default) fills in the sampler's \code{\link[=dbartsControl]{control}} object's \code{updateState}, and explicit \code{TRUE}/\code{FALSE} override it. For the mutators - \code{setData}, \code{setResponse}, \code{setOffset}, \code{setWeights}, \code{setForestBasis}, \code{setSigma}, \code{setCalibration}, \code{setPredictor}, and \code{setCutPoints} - the state is stored only on explicit \code{TRUE}; \code{NA} (the default) and \code{FALSE} both store nothing, regardless of \code{control@updateState}. These are typically called once per sweep inside a larger Gibbs/MH loop (as \code{dbartsSampler} is designed for), where storing state on every mutation would be wasted work whenever the loop only reads \code{state} occasionally (or never); an unforced \code{state} promise materializes the sampler's \emph{current} state on first access regardless, so a mutate-then-first-read sequence needs no explicit store. Pass \code{TRUE} explicitly when \code{state} was already forced (read or saved) earlier and a later mutation must be reflected in the next save.
  }
  \item{shallow}{
    A logical determining if the copy should retain the underlying data of the sampler (\code{TRUE}) or have its own copies (\code{FALSE}).
  }
  \item{newControl}{
    An object inheriting from \code{\link{dbartsControl}}. When passed to \code{setControl}, it cannot change \code{n.trees}, \code{n.chains}, \code{useQuantiles}, or \code{seed} from the values the sampler was created with, and cannot set \code{keepTrees = TRUE} without also giving \code{n.samples}; either is an error.
  }
  \item{newModel}{
    An object inheriting from \code{dbartsModel}. When passed to \code{setModel}, it cannot switch a DART tree prior on or off relative to the sampler's creation-time model; recreate the sampler instead.
  }
  \item{newData}{
    An object inheriting from \code{\link{dbartsData}}.
  }
  \item{y}{
    A numeric response vector of length equal to that with which the sampler was created, and lying in the response family's support: 0/1 for \code{probit} and \code{logistic}, an integer category index in \eqn{[1, K]}{[1, K]} for \code{ordinal}, and a finite non-negative integer count no larger than \eqn{10^6}{1e6} for \code{nbinom} (the dispersion grid's count histogram is sized from the largest count, so a larger one allocates without bound). Values off the support are refused, as they are at creation; \code{gaussian} and \code{aft} (log survival times) constrain nothing.
  }
  \item{x}{
    A numeric predictor vector of length equal to that with which the sampler was created. Can be of a distinct number of rows for \code{setTestPredictor}.
  }
  \item{x.test}{
    A new matrix, data frame, or sparse-bearing test set (the column types accepted by \code{\link{dbartsData}}), of the number of columns equal to that in the current model. A data-frame/sparse \code{x.test} is coded against the training levels; \code{setTestPredictor} and \code{setTestPredictorAndOffset} then install it as a resident, whole-object replacement of the test set (see \code{column}), while \code{predict} and \code{predictForests} route its rows through the trees off that same storage, retaining nothing.
  }
  \item{offset}{
    A numeric vector of length equal to that with which the sampler was created, or \code{NULL}. If \code{offset.test} was set from \code{offset}, will attempt to update that as well.

    For \code{setCategoryOffset}, a multinomial sampler's per-category shift instead: an \eqn{n \times K}{n x K} numeric matrix, or \code{NULL} to clear one. The latent becomes \eqn{f_{ik} + o_{ik}}{f_ik + o_ik}, so the shift enters the log-sum-exp margins, every category's working response and the reported probabilities, and never a leaf value. This is the response-side counterpart of \code{setCounts} rather than of \code{setOffset}, whose flat shift is added AFTER the categories are blended - the wrong side of the nonlinearity - and is the softmax's own null direction besides, which is why \code{setOffset} itself is refused on such a sampler. Only the row-centred part is identified: adding a constant to a whole row leaves every reported probability unchanged, and the entrance leaves the matrix as given rather than re-centring it. It shifts the TRAINING latent only. Mirrored into \code{data@offset.category}, so a re-created sampler carries it.
  }
  \item{updateScale}{
    Logical indicating whether BART's internal scale should re-anchor to the new offset (\code{setOffset}) or response (\code{setResponse}). Defaults to \code{FALSE}, locking the scale set at creation; should only be \code{TRUE} during burn-in, as re-anchoring mid-run makes the fits across iterations no longer comparable. \code{TRUE} is refused by samplers whose leaf calibrations are stated against the transform fixed at creation and are never restated: a multi-forest (\code{forests}) sampler - refused whatever its response family, the calibration map being pinned at creation whether or not there is a transform to re-anchor - and a heteroscedastic (\code{variance}) one, whose variance forest would otherwise keep reporting \eqn{s^2(x)}{s^2(x)} on the abandoned scale. A sampler carrying grouped random effects (see \code{\link{rbart_vi}}) is refused on the same grounds, but only where there is a data-derived transform to abandon: under a gaussian, Student-t (\code{resid.dist = student}) or \code{"aft"} response the random intercepts \eqn{b}{b} and their scale \eqn{\tau}{tau} are held against the transform fixed at creation and nothing converts them, so \code{TRUE} would silently restate both in response units while \code{sigma} moved with the scale; under \code{"probit"} or \code{"logistic"} the transform is the link's own and \code{TRUE} is accepted as the no-op it already is. On \code{setResponse}, supplying this argument positionally (\code{setResponse(y, TRUE)}) rather than by name warns once per session - it is the second argument, where a caller porting code written before this order had \code{updateState} there instead.
  }
  \item{offset.test}{
    A numeric vector of length equal to that of the test matrix, or \code{NULL}. Can be missing for \code{setTestPredictorAndOffset}. Refused, by name, for \code{predictForests}: an offset shifts the combination of the forests, exactly as the response transform's shift does, and neither is any one forest's own total.

    For \code{setCategoryTestOffset}, a multinomial sampler's \eqn{n_{\mathrm{test}} \times K}{nTest x K} per-category test shift, or \code{NULL} to clear one: the recorded test channel becomes \eqn{\mathrm{softmax}(f_{\mathrm{test}} + o_{\mathrm{test}})}{softmax(f_test + o_test)}, formed where the training blend forms its own. The test fits enter no likelihood, so this moves the reported test probabilities and nothing else. Its rows are the CURRENT test rows, so replacing those rows while it is installed is refused rather than silently reinterpreted - clear it first - and out-of-sample \code{predict} does not read it at all, taking its own matrix for the rows it is given. Mirrored into \code{data@offset.category.test}. For \code{predict} on a multinomial sampler, likewise a per-category matrix, one row per PREDICTED row: a flat vector is refused there for the same reason, and a sampler holding either resident category offset refuses an unnamed call rather than reporting the offset-free surface, which an explicit all-zero matrix asks for on purpose.
  }
  \item{n.threads}{
    For \code{predict} and \code{predictForests}, a per-call worker count that
    does not persist, defaulting to the sampler's own (the \code{n.threads} of
    its \code{\link[=dbartsControl]{control}} object, fixed at creation). The
    replay is partitioned by (chain, saved draw), each partition writing its
    own rows and nothing being reduced across workers, so the result is
    identical bit for bit at every value and only the time taken changes. A
    replay small enough that the partition could not pay for itself runs
    inline whatever is passed.

    These are native threads, so \code{_R_CHECK_LIMIT_CORES_} does not govern
    them: that variable is read only by \code{parallel:::.check_ncores}, which
    counts simultaneously spawned processes and never sees a thread.

    For \code{run} the argument still has no effect: the sampler's own count
    governs, and this one is reserved for a future per-call override.
  }
  \item{sigma}{
    A single positive residual standard deviation, on the original response scale; it is assigned to every chain. \code{setSigma} only applies to a sampler whose residual standard deviation is a free parameter: it is refused on a probit-, logistic-, multinomial-, ordinal-, or count-family sampler, which fix it at 1 by the model definition, and on a heteroscedastic (\code{variance}) sampler, whose variance forest owns the residual scale. A gaussian sampler with \code{resid.prior = fixed()} is not fixed in this sense and is still accepted: suppressing the internal draw so an outer sampler owns \code{sigma} is what that prior is for.
  }
  \item{counts}{
    For \code{setCounts}, the replacement response of a multinomial (softmax) sampler: an \eqn{n \times K}{n x K} matrix of non-negative integer counts, column \eqn{k}{k} holding category \eqn{k}'s successes, with trials \eqn{n_i = \sum_k}{n_i = sum_k} \code{counts[i, ]} at least 1. Both \eqn{n}{n} and \eqn{K}{K} are fixed at creation - every combiner buffer is sized by \eqn{n}{n}, and \eqn{K}{K} is the forest count - so only the values may change; a matrix of the wrong shape is refused naming the count in force. The trees carry over, fitted to the previous counts exactly as \code{setResponse} leaves a single-forest sampler's, and the next \code{run} forms every category's working response against the new matrix. The matrix is written to both the engine and \code{data@counts} (its row sums to \code{data@y}), so \code{getPointer}'s transparent re-creation after save and load carries the current response rather than the one the sampler was created with. Cost, not a defect: the sweep draws \eqn{n_i}{n_i} Polya-Gamma variates per observation per category, so replacing single-trial labels with grouped counts multiplies sweep cost by \code{mean(n_i)}. Refused, naming the reason, on any sampler that carries no count response.
  }
  \item{weights}{
    A numeric vector of case weights, of length equal to that with which the sampler was created. What a weight MEANS is family-specific: for a gaussian-family sampler it is a precision, non-negative, and a weight of zero excludes an observation from the likelihood while keeping its fitted values; for a \code{logistic}-family sampler it is an observation COUNT and must be a positive integer, a zero count being a dropped row rather than a down-weighted one (\code{active} is the supported way to take a row out of the data set for a sweep). Installed BETWEEN samples on a forest that has already grown, a vector of zeros can leave leaves that no positive-weight row reaches - the trees were drawn before it existed, and weights are not part of the saved \code{state} - which is a legal, transient state rather than an error: such a tree keeps moving under the tree prior at a constant likelihood and returns to leaves the likelihood reaches. \code{setWeights} applies to gaussian- and \code{logistic}-family samplers. A logistic swap is a model change with a defined meaning rather than a reweighting: the counts are the shape of the Polya-Gamma latents, which are redrawn against the new counts before the call returns - from the sampler's own generators, not R's stream - so an outer sampler can vary exposure between runs. \code{probit}, \code{ordinal}, \code{aft} and \code{nbinom} refuse it by identification: a weighted probit has no tractable latent-variable form, \code{ordinal} inherits that, \code{aft} fixes its censoring structure at creation, and \code{nbinom}'s Polya-Gamma shape is \eqn{y_i + r}{y_i + r} with no weight slot. \code{setData} carries the same rule on the whole-data conduit, and redraws the latents against the counts the replacement data carries; replacement data given without weights is single-trial, as at creation, so a logistic sampler built with counts and handed weightless data becomes an unweighted one. Under an installed mask a swap redraws only the ACTIVE rows - an inactive row consumes no random numbers and returns to its deterministic cold start against the new count. The weights themselves are not part of the saved \code{state}, but a digest of the ones in force when it was stored is: \code{setState} compares it against the destination's own weights and, where the two disagree, re-derives the weight-dependent latents against the DESTINATION's - so a restore lands where the same \code{setWeights} call would rather than pairing one vector's latents with another's counts, silently and off the restored generators. A matched round trip re-derives nothing and installs the stored state unchanged. Only a family whose augmentation is stated against the weights moves under it (\code{logistic}); for gaussian, Student-t and every weight-refusing family it is a no-op. See \code{\link{dbarts}} for the family-specific weight rules that apply at creation time.

    For \code{setForestWeights}, a distinct per-FOREST case weight on a Bayesian causal forest (built with \code{forests = }, see \code{\link{dbarts}}), refused off one: it composes with \code{weights} and \code{active} as \eqn{(w_i a_i) m_f^2 s_i}{(w_i * a_i) * m_f^2 * s_i}, where \eqn{m_f}{m_f} is the named forest's own multiplier, rather than widening either channel, so \code{s_i = 0} excludes row \eqn{i} from \emph{that forest's own leaf conditionals only} - its occupancy, its place in the combination, and the residual degrees of freedom are unaffected, and the reported location still carries the forest's full contribution even when every row is excluded. The third leg of a three-way degenerate-value contrast: \code{weights} installs and is measurably distinct from carrying none, \code{active} at all-ones installs nothing and clears any mask in force, and \code{setForestWeights} at all-ones INSTALLS - a round trip reports it - but is bitwise IDENTICAL to carrying no per-forest weight at all, the null gate its multiplicative composition (\eqn{m_f \times 1 = m_f}{m_f * 1 = m_f}) is built to guarantee. The weight does not ride the sampler's saved \code{state} - a sampler rebuilt with \code{setState} from a stored state silently drops it while \code{statesAgree} still reports agreement - so it is additionally mirrored on an R5 field that \code{getPointer}, \code{setState}, and \code{copy} all reinstall on every re-creation; a caller never reinstalls it by hand. \code{setData} needs no clearing rule here the way it does for \code{active}: it is refused outright on any multi-forest sampler (\dQuote{a multi-forest sampler fixes its data at creation}), so the two channels never interact.
  }
  \item{active}{
    A numeric vector of per-observation 0/1 membership indicators - "row \eqn{i} is not in the data set for this sampler this sweep" - of length equal to that with which the sampler was created, or \code{NULL}. \code{setActiveRows} reaches every family, the four that refuse case weights entirely included: gaussian, Student-t, \code{probit}, \code{ordinal}, \code{logistic}, \code{nbinom}, and \code{aft} samplers all accept it, including a multi-forest (\code{forests}) sampler, whose forests share one response of whichever family it was built under. It is absolute and independent of \code{weights}, composing with it in either call order (\eqn{w_i \times a_i}{w_i * a_i} in effect). \code{NULL} clears the mask (every row active). Length and \code{NA} are validated as for \code{weights}; the only legal values are exactly \code{0} and \code{1}, and a fractional element refuses the WHOLE call and installs nothing. An all-ones vector reports success but installs NOTHING, clearing any previously installed mask - the OPPOSITE of \code{weights}, where an all-ones vector installs and is measurably distinct from carrying no weights at all. One channel states membership, the other precision, and the two deliberately carry opposite degenerate-value policies.

    An inactive row (\code{active[i] == 0}) contributes nothing to any leaf sufficient statistic, branch log-likelihood, birth-scan weight total, leaf parameter draw, or family-level parameter update that sums over rows (residual degrees of freedom, a dispersion statistic, a group's per-group sums, and so on). It still occupies its leaf for COUNT-based accounting - \code{numObservations}, the birth scan's own member count, and leaf collapsing, which triggers on member count regardless of weights - and it still receives a fitted value: \code{run()$train}, \code{getForestFits}, \code{getFitsWithoutOffset}, and \code{predict} report \eqn{f(x_i)}{f(x_i)} at an inactive row exactly as at an active one, which is what makes this channel worth more than physically dropping the row. The empty-leaf VETO is the one exception, and it is conditional: it counts POSITIVE-WEIGHT members rather than members, so it degenerates to the member count exactly when no weight vector is installed - but a mask IS a weight vector, so with one installed a leaf all of whose rows are inactive is vetoed rather than counted as occupied. A mask installed mid-run on a grown forest can therefore strand a leaf the trees were drawn without it: the affected tree keeps moving under the tree prior at a constant likelihood and is absorbed back into the vetoed-free set as those leaves clear. With every row inactive that is the whole forest, which is the exact sense in which it \dQuote{sits at its prior} - a distribution over structures, not a frozen one.

    An inactive row's own latent draw is SKIPPED - no random numbers are consumed for it - for every family except Student-t, whose per-row \eqn{\lambda_i}{lambda_i} is drawn unconditionally at every row regardless of the mask; there the mask instead ANNIHILATES its contribution through the composed weight (\eqn{w_i \lambda_i a_i}{w_i * lambda_i * a_i}). For a skipping family, \code{getLatents} at an inactive row returns its LAST DRAWN value, which is stale; the correct read at an inactive row is its fit, not its latent. Reactivating a row (a mask change from inactive to active) is a one-sweep MODEL hazard, not only a read hazard: that sweep's tree moves and leaf draws run against the row's stale latent (or, for Student-t, a \eqn{\lambda_i}{lambda_i} drawn while the row was out) before that sweep's own latent refresh updates it; this does not disturb the posterior while a mask is held fixed, but a caller that moves the mask every sweep should expect it.

    An inactive row's pointwise log-likelihood - the value the flat C API's \code{logLikelihood} results field reports (see \sQuote{Mutation cost} above for the C interface) - is \code{NaN} rather than a finite value or \code{-Inf}: a row that is not in the model has no likelihood to report.

    Three quantities do NOT follow the mask, by design, so a caller wanting them recomputed on the active subset must re-create the sampler instead: the split-candidate cut grid stays the FULL-data grid, holding the tree prior fixed while the row set moves; the response transform and residual prior scale stay the FULL-data calibration (a masked gaussian - and so \code{aft}, and any \code{variance} or grouped decorator - keeps the range and residual-prior scale of every row); and an ordinal sampler keeps its FULL-data number of categories, free thresholds, and log-gap prior even when a mask empties a boundary category - only \code{setData} changes them. \code{setData} also CLEARS an installed mask - the mask is a per-row vector and \code{setData} may change the number of rows - so a caller replacing the data must reinstall any mask it wants kept.

    The mask does not ride a sampler's saved \code{state}: a sampler rebuilt with \code{setState} from a stored state silently drops the mask and computes different draws while \code{statesAgree} still reports agreement. It is DROPPED rather than reconciled, which is where it parts from the case weights: a destination that never had a mask has no other mask for the restore to re-derive against, only an absent channel, so a caller that wants one reinstalls it by hand. Conversely, the latents DO ride the state, and under an installed mask they are stale at every inactive row for a skipping family, so two samplers at the same posterior state but different mask histories carry different stored latents and \code{statesAgree} correctly reports DISAGREEMENT; a Student-t sampler is affected the same way, since an inactive row's stored \eqn{\lambda_i}{lambda_i} is drawn from a conditional its own data no longer informs.

    On a sampler with grouped random effects, a group all of whose rows are inactive draws its effect from the PRIOR through the group-effect update's own formula - the coherent answer, and NOT what a compacted sampler with that group physically deleted would do.
  }
  \item{basis}{
    The data the named forest's amplitudes multiply, with as many rows as the sampler was created with and one amplitude per column: a factor (or a one-sided formula naming one) expands to its level indicators with no reference level dropped - a two-level factor giving the pair whose amplitudes are \eqn{(b_0, b_1)}{(b0, b1)} - and a numeric vector or matrix is already those columns, which is the same expansion \code{\link{forest}} applies at creation. Any forest takes a basis of any width. \code{setForestBasis} is the \emph{sole} route by which a basis changes after creation, and it applies only to a sampler whose forests carry amplitudes, built with \code{forests = } (see \code{\link{dbarts}}). The amplitudes are preserved and remapped: a width-preserving install leaves every one of them bitwise, and a width change carries each forest's block to its new offset and enters the added coordinates at \code{1}. A basis that is not one of the canonical shapes (a dense all-ones column, or a complementary 0/1 pair) moves that forest onto the general amplitude conditional, which is a model fact and not revertible by a later data swap. The matrix is written to both the engine and \code{data@bases}, so \code{getPointer}'s transparent re-creation after save/load carries the current assignment rather than the one the sampler was created with; because the bases ride creation, a widening applied after a state restore preserves the \emph{restored} amplitudes.
  }
  \item{forest}{
    A single positive integer indexing the forests from \code{1}, or - for \code{getForestFits}, \code{getForestAmplitudes}, \code{getForestVariableCounts}, and \code{getCalibration} - \code{NULL} (the default on all four) for every forest, stacked (see \sQuote{Value} for each reader's stacked shape). A single-forest sampler's \code{NULL} read is bitwise its \code{forest = 1} read, so the default costs nothing on an ordinary sampler. \code{setForestWeights} and \code{setForestBasis} have no default: a writer names one target rather than reading every one. \code{setCalibration} keeps \code{forest = 1L}: it is refused on every multi-forest sampler, since a calibration map owns those forests' scales, so \code{1} is the only value that can succeed and \code{NULL} would name nothing writable. A Bayesian causal forest's prognostic forest is \code{1} and its basis forest \code{2}; \code{setForestWeights} and \code{setForestBasis} are both refused on a sampler whose forests carry no amplitudes, but not with the same message - \code{setForestWeights} names the missing capability, while \code{setForestBasis} raises \code{"forest index out of range"} - and \code{setForestBasis} accepts any forest of one that does. \code{getForestFits} and \code{getForestVariableCounts} accept \code{forest = 1} on any sampler - it selects the only forest - and refuse only an out-of-range index. \code{getCalibration} is likewise served on every forest of a multi-forest sampler, and its calibration-map columns are that forest's own.
  }
  \item{prior.scale}{
    For \code{setCalibration}, the prior standard deviation of the forest total at \code{k = 1}, in response units (the family's latent units where the response is not rescaled). This is the identified quantity: only the ratio of the leaf scale to \code{k} enters a draw law, so under a fixed \code{k} the pair has one degree of freedom, and \code{prior.scale} is the half of it a hyperprior on \code{k} leaves alone. Must be a single positive finite number; exactly one of \code{prior.scale} and \code{prior.sd} is given.
  }
  \item{prior.sd}{
    For \code{setCalibration}, the same statement at the \code{k} currently in force, so \code{prior.scale = prior.sd * k}. Refused when \code{k} is drawn from a hyperprior - it would name a different prior every sweep - and refused when the chains' \code{k} have diverged, since one number would then mean a different scale on each; \code{prior.scale} serves in both cases. Note the binary families default to a sampled \code{k}.
  }
  \item{prior.mean}{
    Not writable, and present only to say so with its remedy: the leaf values it would shift are already drawn and stored. The prior mean of the forest total is the response transform's shift, and the lever that moves the modelled quantity is the offset channel, \code{setOffset(rep_len(-getCalibration()[1, "prior.mean"], n))}.
  }
  \item{cuts}{
    Vector of cut points for use with \code{setCutPoints}, which applies to numeric predictors only: a factor column of either kind is refused, naming the kind, its grid following its level table rather than any externally chosen one.
  }
  \item{column}{
    An integer or character string vector specifying which column/columns of
    the predictor matrix is to be replaced or for which cuts should be used.
    When updating predictors, it can be missing and the entire matrix will be
    substituted; a sparse or mixed dense/sparse design accepts that too, and
    keeps its storage (the replacement is spliced into the container rather
    than replacing it with a plain matrix). Such a design may equally be
    mutated one column at a time by naming \code{column} (dense-backed and
    sparse-backed columns alike, and a single call may name both kinds), which
    is the cheaper spelling when only a few columns move. Note that replacing
    a sparse-backed column densifies its storage permanently: every row then
    differs from the column's implicit value, so it stores \code{n} entries
    from that point on. For \code{setTestPredictor}, a per-column update
    reaches a data-frame or sparse test container too: the replacement is
    spliced into its storage rather than densifying the whole test set.
  }
  \item{forceUpdate}{
    For \code{setPredictor}, controls how the trees respond to a new predictor
    that would leave a leaf empty in any tree of any forest of any chain - the
    ordinary mean forest, and, on a sampler with more than one forest, all of
    them together: a Bayesian causal forest's second (treatment) forest, a
    multinomial sampler's per-category forests, and a heteroscedastic
    (\code{variance}) sampler's variance forest are checked on exactly the same
    terms as the mean forest, not independently (see \sQuote{Multi-forest and
    heteroscedastic predictor mutation} below). \code{FALSE} rolls the whole
    update back and leaves the sampler unchanged (rejection sampling); \code{TRUE}
    forces the change, collapsing any empty nodes; the string \code{"partial"}
    installs each observation's new value individually and rolls back only those
    observations whose value would empty a leaf in any of those trees, returning
    a per-observation logical of what was installed. \code{"partial"} requires a
    single \code{column} and cannot be combined with \code{updateCutPoints}; it
    also requires a dense-backed column, since a sparse column's storage fixes
    its nonzero pattern per cell - replace such a column whole instead. When
    missing, defaults to \code{TRUE} when the whole predictor matrix is being
    replaced (\code{column} is missing) and \code{FALSE} when a single column is
    being replaced.
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
    An integer vector listing the saved samples to return from \code{getTrees}
    or print from \code{printTrees}. Applies only when \code{keepTrees} is
    \code{TRUE} and \code{current} is \code{FALSE}; otherwise the live working
    trees have no sample dimension and it is ignored. Sample numbers address
    recorded DRAWS on the oldest-first axis \code{predict} reports - \code{1} to
    the number of draws recorded, at most \code{control@n.samples} - and not
    the store's internal slots.
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
    allocate. Refused by name by \code{getSigmas} and
    \code{getSumsOfSquaredResiduals}: both readers allocate their own vector,
    and filling a caller's buffer in place is \code{getLatents}'s contract
    alone.
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
  \item{ptr}{
    For \code{storeState}, the sampler's engine pointer to capture from; defaults to \code{getPointer()} and is rarely passed explicitly.
  }
  \item{newState}{
    For \code{setState}, a state object previously produced by this sampler over the same model: its \code{state} field, calling \code{storeState()} first to refresh it if a later mutation may have gone unstored (\code{storeState} itself returns \code{NULL} invisibly, so its return value cannot be passed here). Must inherit from \code{bartcoreState}.
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
      The cached, serializable engine state, or \code{NULL} until one is materialized. The saved-tree store's write position and the number of draws it has recorded both ride it, so \code{predict} after a \code{setState} reports the same draws, in the same order, as before the store. Reading it forces the sampler's \emph{current} state only the first time, before any value has been materialized; once set, it is a cached snapshot that a later mutation does not refresh automatically - call \code{storeState} again, or pass \code{updateState = TRUE} to the mutating call (see \code{updateState} above), to bring it forward. It is the only field \code{\link{save}} needs, and restoring one requires \code{setState} - see \sQuote{Saving}.
    }
    \item{\code{pointer}}{
      The external pointer to the C++ engine. Internal; never manipulate it directly.
    }
  }
}
\details{
  A \code{dbartsSampler} is a mutable object which contains information pertaining to fitting a Bayesian additive regression tree model. The sampler is first created and then, in a separate instruction, run or modified. In this way, MCMC samplers can be constructed with BART components filling arbitrary roles.

  \subsection{Reading the fit}{
    Six methods report a fitted quantity and no two report the same one. They differ along three axes - the SCALE the values live on, whether the installed OFFSET is in them, and whether they answer from the CURRENT state or from stored draws:

    \tabular{llll}{
      \strong{method} \tab \strong{scale} \tab \strong{offset} \tab \strong{current or stored} \cr
      \code{run()$train} \tab response (latent under a binary family) \tab included \tab one slab per recorded draw \cr
      \code{$getFitsWithoutOffset()} \tab response (latent under a binary family) \tab excluded \tab current \cr
      \code{$getForestFits(f)} \tab internal, and one forest only \tab excluded \tab current \cr
      \code{$predict(x)} \tab response (latent under a binary family) \tab whatever you pass it \tab the saved samples under \code{keepTrees}, else the current trees \cr
      \code{$predictForests(x)} \tab internal, one channel per forest \tab excluded, and refused if passed \tab the saved samples under \code{keepTrees}, else the current trees \cr
      \code{$getLatents()} \tab the family's own augmentation variable, which is a location for some families and a precision for others \tab not a fit; see \sQuote{Value} \tab current
    }

    A multinomial (softmax) sampler reads and writes a strict subset of this surface, by the model rather than by omission. Its response channels are \code{$setCounts}, \code{$setCategoryOffset} and \code{$setCategoryTestOffset}; the predictor family (\code{$setPredictor} in all four shapes, \code{$setCutPoints}, \code{$setTestPredictor}) and the global \code{$setActiveRows} are open, and \code{$predict} answers with the \eqn{K}{K} category probabilities and takes the matching per-category offset. Everything else on the mutation surface is refused, each by a message naming the capability and, where one exists, the channel that serves the caller instead: \code{$setResponse} and \code{$setOffset} (a flat vector cannot express the matrix, and a flat shift is the softmax's null direction), \code{$setWeights} (an integer case weight is already row-wise replication in the count response, and a non-integer one has no exact augmentation sampler), \code{$setSigma} (no residual scale), \code{$setData} and \code{$setModel} (the K category forests fix their data and their calibration at creation), \code{$setCalibration} (the softmax map owns every category forest's leaf scale), \code{$setForestWeights} and \code{$setForestBasis} (its forests are its categories, which carry no amplitudes), and \code{$getFitsWithoutOffset} (its reported channels are probabilities, not one additive location). Per-forest masking is refused permanently rather than pending: a category's margin is a log-sum-exp over the other \eqn{K-1}{K - 1}, so a mask on one forest alone has no conditional to be a mask of.

    The composition rule follows from the table: an outer block conditions on \eqn{f(x_i)}{f(x_i)}, so it reads \code{$getFitsWithoutOffset()} and adds back whatever offset it installed, rather than differencing \code{$getLatents()} against \code{run()$train}, which mixes an offset-bearing quantity with one that does not carry the offset.

    Where the boundary runs. These are the methods of the mutable sampler, and their audience is a model writer composing BART into a larger scheme: they return engine primitives on the engine's own terms. The R modelling conventions - \code{\link{fitted}}, \code{\link{predict}}, \code{\link{extract}} - apply to FIT objects instead, the results of \code{\link{bart}}, \code{\link{bart2}} and \code{\link{rbart_vi}}. Those accessors read the stored \code{yhat.train}/\code{yhat.test} channels the engine already wrote the offset into, and no \code{type} arm removes it: \code{type = "bart"} returns those draws as they stand, \code{"ev"} maps them through the response transform, \code{"ppd"} samples from them, and \code{"loglik"} evaluates against them - all offset-inclusive. (\code{extract} on a \code{dbartsSampler} itself takes \code{type = "predictors"} only and returns the design matrix; it is not a fitted-value accessor.) So an offset-free fit is available from the sampler surface and nowhere on the fit-object surface, which is the division rather than an omission.
  }

  \subsection{Saving}{
    \code{\link{save}}-ing and \code{\link{load}}ing a \code{dbarts} sampler for future use requires that R's serialization mechanism be able to access the state of the sampler which, for memory purposes, is only made available to R on request. To make it available, call \code{storeState()} on the sampler before saving, e.g. for the object \code{sampler}, execute \code{sampler$storeState()}; this captures the sampler's current state into the serializable \code{state} field regardless of the \code{updateState} setting.

    The state object is opaque and engine-specific, and carries a format version and the package version that wrote it. A state loads only within the format version of the \code{dbarts} that wrote it; loading it with an incompatible version refuses cleanly, naming both versions, rather than risk a silent misread. There is no cross-version migration: re-fit the model, or restore the state with the \code{dbarts} release that wrote it.

    To restore a saved state into a sampler, call \code{setState(newState)}: it validates that \code{newState} inherits from \code{bartcoreState}, re-creates the underlying engine if needed, pushes the state into it, and caches it on the \code{state} field. Validation covers every forest the state carries: a state whose trees split outside the recipient forest's allowed columns - a \code{blocks}-constrained or moderator-restricted mean forest, or a \code{variance = ~ x1 + x2} variance forest - is refused with the message \code{installTrees} gives the same donor, the two entries sharing one rule so neither admits what the other refuses. Every check runs before any live state is touched, so a refused restore leaves the sampler exactly as it was. Assigning the field directly (\code{sampler$state <- newState}) does \emph{not} restore the sampler - it only overwrites the R-side cache, leaving the engine untouched, so the next run continues from the engine's own state rather than the assigned one. Always route a restore through \code{setState}. The case weights are not in the state, so a restore is reconciled against the DESTINATION's own rather than governed by the source's: where they differ from the weights the state was stored under, the weight-dependent latents are re-derived against the destination's before \code{setState} returns (see \code{weights} above).
  }

  \subsection{Mutation cost}{
    Each accepted predictor mutation does two things: the engine updates its internal representation, and the R layer collects the accepted values into the sampler's data object so that replay, saving, and re-quantization see current data. The collection is by reference for data-frame-built samplers (only the affected column changes hands) but requires copying the full predictor matrix when the sampler was built from a matrix, so tight loops that mutate predictors every iteration are better served by data-frame input. The remaining per-call overhead is a few tens of microseconds of R method dispatch and bookkeeping - the price of R-level consistency. Clients for which that matters can drive the sampler through the C interface (see the \code{dbarts.h} header installed with the package), which invokes the engine directly and performs no R-side collection; such clients supply their own current predictor matrix when replaying saved trees.
  }

  \subsection{Infrastructure methods}{
    \code{getPointer()} returns the sampler's underlying engine pointer, checking it for validity first. If the engine the cached pointer refers to is no longer live (typically after a \code{\link{load}}), it transparently re-creates one from \code{control}, \code{model}, \code{data} and the stored \code{state}, reinstalls every recorded forest weight (see \code{reapplyForestWeights} below, since forest weights do not ride the saved state) so the replacement matches a freshly built engine, and only then replaces the cached pointer; a re-creation that fails leaves the sampler exactly as it was. \code{adoptPointer(ptr)} is the write side of the same relationship: it rebinds the sampler onto \code{ptr}, an \code{externalptr} some caller built directly, outside this object, from this sampler's own \code{(control, model, data)} triple - typically by driving the C interface described under \sQuote{Mutation cost} - in place of the engine this object created at construction; the abandoned engine becomes unreachable and is released once by its own finalizer. \code{adoptPointer} trusts the caller that \code{ptr} was built from this object's own triple and performs no independent check. \code{reapplyForestWeights(ptr)} re-installs every per-forest weight mirrored on the sampler onto a freshly (re-)created \code{ptr}; it is called only from \code{getPointer} and \code{setState}, never directly.

    These three are the only R-level primitives for exchanging engine pointers with code that drives the engine directly through \code{dbarts.h}. There is no other R-level path to that boundary: the low-level \code{.Call} wrappers this class's methods delegate to are internal to the package and unexported, so a consumer below this reference class has \code{dbarts.h} and nothing else.
  }

  \subsection{Multi-forest and heteroscedastic predictor mutation}{
    \code{setPredictor} - whole-matrix, single- or multi-column (via \code{column}), and per-observation (\code{forceUpdate = "partial"}) - and \code{\link{updatePredictorPerObservationJointly}} accept a Bayesian causal forest, a multinomial sampler, and a heteroscedastic (\code{variance}) sampler on the same terms as an ordinary single-forest sampler. The acceptance criterion is a single conjunction over every tree of every forest of every chain of the sampler: a treatment forest, a multinomial sampler's per-category forests, and a variance forest are all checked together with the mean forest, never independently, so a change that a mean-forest tree would tolerate can still be rejected because it would empty a leaf in the treatment or variance forest, and vice versa. A rejected transactional whole-matrix or column update rolls \emph{every} forest back to its pre-call state; under \code{forceUpdate = "partial"} the per-observation install mask reflects that same joint criterion, so an observation installs only where it empties no leaf anywhere in the sampler, and a rejected observation is rolled back in every forest, leaving the rest of the run untouched. The one exemption is a forest whose trees never split on the column being changed - a treatment or moderator forest restricted to a column subset, or a \code{variance} forest declared over its own predictor subset - which structurally cannot veto a change to a column outside its own reach; nothing needs to be done to obtain that exemption; it falls out of which columns the forest's trees actually use. \code{updatePredictorPerObservationJointly} applies the identical per-sampler conjunction across every sampler passed to it as well: an observation installs only where it would remain valid in the union of every sampler's forests, and it is declined in all of them otherwise.

    Column names after a whole-matrix replacement: a transactional or forced whole-matrix \code{setPredictor(x, ...)} (\code{column} missing) carries column names onto the installed \code{sampler$data@x}. When \code{x} itself has \code{dimnames}, those names are installed; when it does not (a bare matrix, or a vector reshaped through \code{dim}), the sampler's existing column names are kept instead. Either way, a later \emph{by-name} column reference on that sampler - \code{setPredictor}'s own \code{column} argument given as a character string, or \code{\link{updatePredictorPerObservationJointly}}'s shared-column match - continues to resolve against \code{colnames(sampler$data@x)} as it did before the replacement.
  }

  \subsection{Grouped random effects}{
    A sampler carrying grouped random intercepts - the one \code{\link{rbart_vi}} returns in \code{$fit} when its prior is built in - is mutable on the response side. \code{setResponse} and \code{setOffset} accept a same-length replacement at the pinned scale (\code{updateScale = FALSE}, the default): the group effects and their scale are a Gibbs block the swap deliberately carries across, exactly as they are carried across a tree sweep, and the group indices are per-observation and unchanged, so every row stays in its group. \code{updateScale = TRUE} is refused under a re-anchoring family (see \code{updateScale} above). \code{setData} stays refused outright: it may change the number of rows, and the grouping is fixed at creation, so there is no coherent reading of the new rows' group membership - re-create the sampler instead.
  }

  \subsection{Warm starts}{
    \code{installTrees} seeds the sampler's forests from a \code{donor} instead of drawing trees from the prior, for scaling to more chains or embedding a fit in a larger sampler. Only the donor's trees, \code{sigma}, and \code{k} transfer; each chain keeps its own random-number stream and redraws everything else, so several chains seeded from one donor stay overdispersed. A donor with a different tree count or DART setting is refused rather than silently reshaped; a donor fit on a different cut grid is instead remapped onto this sampler's grid, collapsing any splits the grid starves, the same way a data replacement remaps existing splits. A warm start biases the early draws toward the donor, so it shortens burn-in rather than removing it; keep drawing a non-zero number of burn-in samples before treating the chain as converged.

    \code{growFromRoot} instead builds the sampler's initial forest by XBART-style root-down stochastic tree construction (He, Yalov and Hahn 2019): each of \code{n.sweeps} sweeps rebuilds every tree from the root, sampling each split from the integrated-likelihood weight of every candidate cut under the tree prior. This reaches a good fit in far fewer sweeps than the exact sampler, so it is a fast starting point rather than a posterior sampler - the exact MCMC sweeps own stationarity once \code{run} begins, and the posterior is unchanged. It is available for the constant-leaf model only; calling it on a \code{linear} or \code{gp} node prior is an error, not a silent fall-back, so initialize those forests with \code{sampleTreesFromPrior} instead. As with \code{installTrees}, the grown forest biases the early draws toward its fit, so shorten burn-in rather than skipping it. Each chain grows on its own random-number stream, so the result is independent of \code{n.threads}. A composable cross-sampler workflow follows for free: \code{donor$growFromRoot(k)} then \code{target$installTrees(donor)}.
  }
}
\value{
  For \code{run}, a named-list with contents \code{sigma}, \code{train}, \code{test}, and \code{varcount} (plus \code{k}, \code{varprobs}, \code{tau}, and \code{ranef} when applicable). \code{train} is an array of dimension n.obs x n.samples x n.chains, and likewise \code{test} (or \code{NULL} if the sampler has no test data) and \code{varcount} are n.predictors x n.samples x n.chains; \code{sigma} is n.samples x n.chains. When \code{n.chains} is \code{1} the trailing chain dimension is dropped, so \code{train} is a plain n.obs x n.samples matrix and \code{sigma} a plain vector of length n.samples. On a multi-forest (\code{forests}) sampler \code{varcount} gains a forest axis between the predictors and the samples - n.predictors x n.forests x n.samples x n.chains, forest-major within a draw and the prognostic forest first - so each forest's own per-draw split counts arrive from one call rather than only the reported forest's; a single-forest sampler's array keeps exactly its n.predictors x n.samples x n.chains shape, and the same widening is what a \code{family = "multinomial"} sampler's per-category counts ride. This is the RAW run shape; the packaged fit \code{\link[=bart]{bart2}} builds reshapes it draws-first with the forest names on the trailing margin, as it does for the multinomial channel. \code{$getForestVariableCounts} reads the same quantity for the CURRENT state, one forest at a time. A Bayesian causal forest adds two more: \code{forestFits}, an n.obs x n.forests x n.samples x n.chains array of each forest's fitted values on the internal scale (the prognostic \eqn{\mu}{mu} first, the treatment \eqn{\tau}{tau} second), and \code{glue}, a sum(q) x n.samples x n.chains array of the amplitudes each draw combines them through, stacked forest-major and as wide as each forest's own basis (a Bayesian causal forest's \eqn{(a, b_0, b_1)}{(a, b0, b1)}, three rows), so both surfaces and their recombination come from a single call; \code{train} carries the combination \eqn{a \mu(x_i) + b_{z_i} \tau(x_i)}{a mu(x_i) + b_z_i tau(x_i)} on the response scale, or on the latent scale when the sampler was built under \code{"probit"} or \code{"logistic"}, and \code{test} is filled with \code{NaN} (there is no test treatment vector to combine off-sample). No other model reports either element. A \code{"nbinom"} sampler adds one of its own: \code{dispersion}, the negative-binomial \eqn{r}{r} each draw is conditioned on, shaped exactly as \code{sigma} (a length n.samples vector at one chain, an n.samples x n.chains matrix otherwise) because it is the count analog of it - fixed at the value the sampler was created with under a fixed \code{dispersion}, and that sweep's grid draw otherwise. It is written from the same state \code{storeState} serializes and consumes no random numbers, so reading it costs a run nothing. No other family carries the element at all: it is absent from the list, not \code{NULL} within it, so \code{run()$dispersion} is \code{NULL} on every non-\code{"nbinom"} sampler and a test of the channel must be \code{!is.null(...)} rather than a comparison, which \code{NULL} would satisfy vacuously. A sampler built with \code{resid.dist = student()} adds \code{resid.df} on exactly the same terms - the degrees of freedom \eqn{\nu}{nu} each draw is conditioned on, shaped as \code{sigma}, written from settled state and consuming no random numbers, absent from the list under any other error law - fixed at the value supplied to \code{student(df = )}, and that sweep's grid draw when the degrees of freedom are estimated. A run can be interrupted with \code{Ctrl-C}: it stops between iterations - joining any worker threads first - and signals an error, returning no samples from the interrupted run. The sampler's chains are left at the iteration they reached, which is a valid state to run again from.

  For \code{setPredictor}, \code{TRUE}/\code{FALSE} depending on whether or not the operation was successful. The operation can fail if the new predictor results in an empty leaf-node in any tree of any forest of the sampler (see \sQuote{Multi-forest and heteroscedastic predictor mutation}). If only single columns were replaced, the update is rolled back so that the sampler remains in a valid state. When \code{forceUpdate} is \code{"partial"}, instead returns a logical vector of length equal to the number of observations, \code{TRUE} where that observation's new value was installed and \code{FALSE} where it was rolled back to its previous value to keep every tree of every forest valid.

  Under \code{keepTrees} the draws \code{predict}, \code{predictForests}, \code{getTrees} and \code{printTrees} report come out OLDEST FIRST: the store keeps the most recent \code{control@n.samples} recorded draws, so a sampler driven through several recorded \code{run} calls replays them in the order they were drawn and lines up draw for draw with the \code{train}/\code{sigma} channels those runs returned. A sampler that has recorded fewer draws than that reports only the draws it holds, and one that has recorded none - \code{keepTrees} set but nothing run since creation, a \code{setControl} that resized the store, or an \code{installTrees} warm start - is refused rather than answered from slots nothing has written.

  \code{predict} keeps the current test matrix in place and uses the current set of tree splits. This function has two use cases. The first is when \code{keepTrees} of \code{\link{dbartsControl}} is \code{TRUE}, in which case the sampler should be run to completion and the function can be used to interrogate the existing fit. When \code{keepTrees} is \code{FALSE}, the function can be used to obtain the likelihood as part of a proposed new set of covariates in a Metropolis-Hastings step in a full-Bayes sampler. This would typically be followed by a call to \code{setPredictor} if the step is accepted.

  Its return shape mirrors \code{run()}'s own \code{train}/\code{test} channel for a single-forest sampler, with two exceptions. An amplitude-coupled (multi-forest, including a Bayesian causal forest) sampler refuses \code{predict} outright rather than answering here - see \code{getFitsWithoutOffset} below, which already records the same refusal - so everything that follows describes a single-forest or multinomial sampler only. For a single-forest sampler: with saved trees to replay (\code{keepTrees} \code{TRUE} and at least one recorded draw), the result is an n.test x n.samples x n.chains array, the trailing chain dimension dropped to a plain n.test x n.samples matrix when \code{n.chains} is \code{1}; without saved trees to replay, the CURRENT (live) trees are replayed once per chain instead, with no samples axis at all - a plain vector of length n.test at one chain, an n.test x n.chains matrix otherwise. A \code{family = "multinomial"} sampler's surface reports the K category probabilities instead: an n.test x K x n.samples (x n.chains) array, the K axis inserted between the rows and the samples. Unlike the single-forest case, this samples axis is never dropped: it is sized \code{1} rather than omitted when there are no saved trees to replay. A heteroscedastic (variance forest) single-forest sampler with saved trees to replay returns not a bare array but a named list \code{list(mean = ..., variance = ...)}, both members the single-forest shape above; without saved trees there is nothing for the variance forest to replay, so the bare mean array returns alone, exactly as a non-heteroscedastic sampler's would.

  For \code{getTrees}, a \code{data.frame} with one row per tree node in depth-first, left-hand-side pre-order, with columns \code{chain} (present only when \code{n.chains > 1}), \code{sample} (present only for saved samples, and reporting the draw number asked for), \code{tree}, \code{n} (the number of observations in the node), \code{var} (the splitting variable, or -1 at a leaf), and \code{value} (the split value, or the leaf prediction). An ordinal rule's value is its cut point and observations with values less than or equal to it go left; on an ordered-factor column that cut point is an exact half-integer at a level boundary, so \code{levels(x)[floor(value) + 1]} is the last level going left (the split value is a 0-based code, R's \code{levels} is 1-based); a categorical rule carries no data value (its \code{value} is \code{NA}) and its split is reported in the \code{directions} column instead. When the sampler has any categorical predictors the result gains a \code{directions} column decoding each categorical rule into one \code{"L"}/\code{"R"} character per level, in level order (level \code{k} goes right when its character is \code{"R"}); ordinal rules and leaves are \code{NA}. When any predictor contains missing values the result gains a \code{missing} column giving the branch (\code{"L"}/\code{"R"}) each rule sends missing values down; rules on complete columns and leaves are \code{NA}. Under a \code{linear} node prior each leaf's \code{value} is its intercept and the result gains one \code{beta.<column>} column per designated covariate holding that leaf's slope on the internal standardized scale; internal nodes are \code{NA}.

  For \code{getSigmas}, a numeric vector of length equal to the number of chains, giving each chain's current residual standard deviation on the original response scale.

  For \code{getDispersion}, the negative-binomial dispersion \eqn{r}{r} currently in force, a numeric vector of length equal to the number of chains, and \code{NULL} on every other family - the count analog of \code{getSigmas}. It is the same scalar \code{run()$dispersion} records once per kept draw, read mid-sweep and without serializing state, so a host driving the sampler one sweep at a time reads it here rather than through \code{storeState()} and \code{state[[chain]]$dispersion}. Because the refusal for a family carrying no dispersion is a \code{NULL} and not an error, a caller distinguishing \dQuote{no dispersion} from a value must test \code{!is.null(...)}; comparing to an expected number would pass vacuously.

  For \code{getSumsOfSquaredResiduals}, a numeric vector of length equal to the number of chains, giving each chain's residual sum of squares \eqn{\sum (y - \hat{y})^2}{sum((y - yhat)^2)} on the original response scale; a binary-response sampler reports on the latent scale instead. On a multi-forest (\code{forests}) sampler \eqn{\hat{y}}{yhat} is the COMBINED location \eqn{\sum_f m_f(x_i) f_f(x_i)}{sum_f m_f(x_i) f_f(x_i)} - the same quantity \code{run()$train} records - and not any one forest's own total. Two qualifications, both keyed to the response family rather than to the forest count. A \code{logistic} sampler's working response is \eqn{w_i(y_i - 1/2)/\omega_i - o_i}{w_i (y_i - 1/2) / omega_i - o_i} with \eqn{\omega_i}{omega_i} a Polya-Gamma variate redrawn every sweep, so what it reports is a function of that sweep's auxiliary variables and is not a residual sum of squares in any scale. A multinomial fit reports \code{NaN}: its coupling reports a category-by-observation probability slab rather than a location, so no residual sum of squares is defined - the \code{NaN} is deliberate, in place of the sum of squared category probabilities a bare substitution would have given.

  For \code{getLatents}, \code{NULL} when the model augments nothing - a plain gaussian response - and, deliberately, on a multinomial one, which augments but reports nothing: its \eqn{K \times n}{K x n} Polya-Gamma \eqn{\omega}{omega} matrix is drawn every sweep, but each \eqn{\omega_{ik}}{omega_ik} is an interleaved ONE-VS-REST augmentation variable refreshed against a margin that the other \eqn{K-1}{K - 1} forests move within the same sweep, so a host reading it between sweeps reads a quantity whose conditioning it cannot reproduce, and the composition recipes reach the same posterior through the fits and the category offset without touching it. That decline is scoped to the augmentation families whose latent is a PRECISION refreshed against a moving margin; the location-valued latents below - probit's and ordinal's \eqn{z_i}{z_i}, and aft's imputed log survival time - are reported as documented. Otherwise the sampler's current draw of the augmentation variable, a plain vector of length equal to the number of observations when there is a single chain, or an observations-by-chains matrix otherwise, written into \code{result} when one was supplied. \strong{What that variable IS depends on the family and is not uniform.} It is a LOCATION, on the sampler's own latent scale, for \code{"probit"} (the truncated normal \eqn{z_i}{z_i}), \code{"ordinal"} (the same \eqn{z_i}{z_i} under the cut points) and \code{"aft"} (the imputed log survival time), all of which a host regresses on directly. It is a PRECISION, one per observation, for \code{"logistic"} and \code{"nbinom"} (the Polya-Gamma \eqn{\omega_i}{omega_i}) and for a Student-t residual distribution (the scale-mixing \eqn{\lambda_i}{lambda_i}); these WEIGHT a working response and are not on the response scale at all, so differencing them against a fit is meaningless. Note the last case: a sampler whose \code{family} is \code{"gaussian"} but whose \code{resid.dist} is \code{student()} (see \code{\link{dbarts}}) does report latents, and they are precisions.

  For \code{getFitsWithoutOffset}, the sampler's combined per-observation location on the response scale and \emph{without} the installed offset, an n.observations x n.chains matrix. \code{run()$train} reports the same quantity with the offset folded in, so this value plus the installed offset is that one; there is deliberately no with-offset twin, since the offset is the caller's own. It is the read a host driving the sampler one sweep at a time needs: an incremental scheme updates its outer block against \eqn{f(x_i)}{f(x_i)}, and taking \code{getLatents()} minus \code{run()$train} instead is BIASED, because the training channel carries the offset the latent scale does not. On a multi-forest (\code{forests}) sampler this is the only route, short of running the sampler, from R to the combined fit - \code{getForestFits} gives one forest's internal-scale totals and \code{predict} is refused. It is refused, naming the reason, on a multinomial sampler, whose reported channels are per-category softmax probabilities rather than one additive location; \code{predict(x)} serves that read, but reports the SAVED samples rather than the current state when the sampler was built with \code{keepTrees}.

  For \code{getForestFits}, a multi-forest sampler's requested forest's current internal-scale fitted values, an n.observations x n.chains matrix, or, at the default \code{forest = NULL}, every forest stacked with the forest margin between the observations and the chains, n.observations x n.forests x n.chains (bitwise the \code{forest = 1} read on a single-forest sampler). For \code{getForestAmplitudes}, the named forest's amplitudes - one per basis column - as a q x n.chains matrix, or, at the default \code{forest = NULL}, every forest's stacked forest-major into a sum(q) x n.chains matrix. The vector is \emph{ragged}: a Bayesian causal forest's forest \code{1} carries the single \eqn{a}{a} on its implicit intercept and its forest \code{2} the pair \eqn{(b_0, b_1)}{(b0, b1)}, so the stacked read is the \eqn{(a, b_0, b_1)}{(a, b0, b1)} of \eqn{y = a \mu(x) + b_z \tau(x) + \epsilon}{y = a mu(x) + b_z tau(x) + eps}. For \code{getForestVariableCounts}, the requested forest's current per-predictor split counts, an n.predictors x n.chains integer matrix, or, at the default \code{forest = NULL}, every forest stacked the same way as \code{getForestFits}, n.predictors x n.forests x n.chains; its rows are named by the training matrix's column names when that matrix carries them (on margin 1 in both shapes), and left unnamed otherwise.

  For \code{predictForests}, the off-sample twin of \code{getForestFits}: every forest replayed separately at \code{x.test}, as an n.test x n.forests x n.samples (x n.chains) array on the forests' internal scale, forest-major within a draw. It reports the saved samples when the sampler was built with \code{keepTrees}, and the current trees otherwise, exactly as \code{predict} does. Nothing but the forests' own totals is in it: no amplitude glue, no response transform and no offset, the last refused rather than ignored. That is the point of the channel rather than an omission - an amplitude coupling's location is \code{response.shift} + \eqn{\sum_f (B_f a_f) \times (\mathrm{response.scale} \times f_f)}{sum_f (B_f a_f) * (response.scale * f_f)}, and off the training rows only the caller knows the bases \eqn{B_f}{B_f}, so the whole recombination is theirs, with \code{$getForestAmplitudes()} and \code{$getCalibration()} supplying the rest. At the fit level it is packaged: \code{predict} on a \code{\link{bart}} object with \code{type = "ev"}, \code{"ppd"} or \code{"bart"} performs the same recombination off this replay, taking the bases at the predicted rows through its own \code{bases} argument. Only a sampler that composes its forests through scalar amplitude glue reports per-forest fits; every other one is refused by name, a multinomial sampler included, whose raw per-category totals are perfectly well defined but whose reported quantity is a softmax probability that \code{predict} serves. The refusal is off-sample only: \code{getForestFits} is not gated the same way, so those same per-category totals stay readable AT the training rows, for the current trees and one draw, and it is the off-sample replay - which would have to hand back a level nothing identifies at rows the sampler never saw - that has no counterpart.

  For \code{getCalibration}, the leaf-prior calibration a forest currently runs under, as a numeric matrix with one row per chain and the columns \code{prior.scale} (the prior standard deviation of the forest total at \code{k = 1}, in response units), \code{prior.sd} (\code{prior.scale / k}), \code{prior.mean}, \code{k}, \code{k.has.hyperprior} (1 when this forest's \code{k} is drawn every sweep, in which case \code{prior.sd} moves every sweep while \code{prior.scale} does not), \code{response.scale}, and \code{response.shift}. The leaf model rides as a \code{"leaf.model"} attribute, one of \code{"constant"}, \code{"monotone"}, \code{"linear"}, or \code{"gp"}, and qualifies what \code{prior.sd} means. At the default \code{forest = NULL}, every forest's calibration is stacked with the forest margin appended LAST, n.chains x 12 x n.forests (bitwise the \code{forest = 1} read on a single-forest sampler) - the row margin here is the chain axis, so there is no row margin for the forest axis to sit behind, unlike the other three readers. The column dimnames and the \code{"leaf.model"} attribute are carried from the first forest, since both are properties of the sampler rather than of any one forest.

  Five further columns report the multi-forest CALIBRATION MAP that fixed \code{prior.scale} on a sampler built with \code{forests =} or \code{dbartsData(bases = )} (see \code{\link{forest}}), and are \code{NaN} on every forest whose scale that map does not own - any single-forest sampler, and a multinomial one, whose scale is not map-derived. \code{amplitude.prior.variance} and \code{amplitude.prior.scale} are EXCLUSIVE per forest: a forest whose amplitudes carry a fixed prior variance reports that variance and a \code{NaN} scale, and one whose amplitude carries the half-Cauchy scale mixture (a forest declaring no basis) reports its median and a \code{NaN} variance. Each is a prior the caller may set; neither moves with the scale mixture's own variance auxiliary, which is a drawn quantity rather than a prior. \code{node.scale.factor} and \code{node.scale.divisor} are the map's two factors and \code{basis.row.norm} the median nonzero row norm of the forest's basis IN FORCE, which \code{setForestBasis} re-derives.

  Together they decompose the reported scale as \code{prior.scale = node.scale.factor * s / (node.scale.divisor * basis.row.norm)}, so the family's own latent anchor \eqn{s}{s} - the only quantity of the map with no column of its own, and data-dependent under a gaussian response - is recovered as \code{prior.scale * node.scale.divisor * basis.row.norm / node.scale.factor}. That identity holds whenever \code{node.scale.factor} is not \code{NaN}. It becomes \code{NaN}, on both \code{node.scale} columns, when \code{setState} or \code{installTrees} installs a leaf scale differing from the one in force: the donor's trees arrive with the donor's calibration, and the stored factor and divisor no longer decompose it. The other three columns are unaffected - the amplitude prior FOLLOWS the installed state, and the bases are not state - and a \code{setForestBasis} on that forest re-imposes the map and restores both columns. Restoring a sampler's own state installs a bitwise-identical scale and so changes nothing.

  \code{prior.scale} and \code{prior.sd} describe the LEAF-PARAMETER scale of the forest total. They equal the prior standard deviation of \eqn{f(x)}{f(x)} at every \eqn{x}{x} for the constant leaf only; for the other three the prior of \eqn{f(x)}{f(x)} is x-dependent and \code{prior.sd} bounds it in a leaf-specific direction. Under \code{"linear"} it is a LOWER bound attained at the standardized covariate origin, with \eqn{sd(f(x)) = }{sd(f(x)) = } \code{prior.sd} \eqn{\sqrt{1 + \|z(x)\|^2}}{* sqrt(1 + ||z(x)||^2)} in the internally standardized leaf covariates (a missing value maps to \eqn{z_j = 0}{z_j = 0}); \code{prior.mean} is exact. Under \code{"gp"} it is an UPPER bound over \eqn{x}{x}, attained at rows that reproduce a leaf member and on over-cap leaves, and elsewhere decaying to zero as \eqn{x}{x} leaves the leaf's data cloud, where every prior draw equals \code{prior.mean} exactly. Under \code{"monotone"} it is a LOWER bound in the interior - the realized standard deviation runs a few per cent to about twenty per cent above it - and \code{prior.mean} is NOT the prior mean of \eqn{f(x)}{f(x)} under an active constraint: that marginal is skew, with an x-dependent mean tracking the constraint's direction across several \code{prior.sd} (see \code{monotone} in \code{\link{dbarts}}).

  This is the authoritative reader of the calibration in force. A model's \code{prior.scale} slot records the named intent and is never rewritten by the engine, so a channel that re-anchors the response transform - \code{setResponse} or \code{setOffset} at \code{updateScale = TRUE}, or \code{setData} - moves what is in force while leaving the intent alone; \code{getCalibration} shows the move, and \code{setCalibration} or \code{setModel(sampler$model)} re-issues the intent.

  For \code{setCalibration}, \code{NULL} invisibly. The write lands on every chain and takes effect on the next sweep, reinterpreting no leaf value already drawn; a write that reproduces what is already in force is skipped bitwise, so a read followed by a write cannot perturb a draw. It is total over the four leaf models, each of which carries the one scale it writes. It is refused on a Bayesian causal forest and on a multinomial sampler, whose per-forest leaf scales come from their own calibration maps; a value that is not a single positive finite number is an error. Nothing else moves - not \code{k}, not the response transform, not \code{sigma}, not the tree prior, not a DART split prior - which is the difference from \code{setModel}, and the reason a DART sampler is served here and refused there. A heteroscedastic sampler's variance forest is a separate leaf model and is not addressable.

  For \code{storeState}, \code{NULL} invisibly; it is called for its side effect of capturing the sampler's current engine state into the serializable \code{state} field, which \code{\link{save}} then writes out. See \sQuote{Saving}.
}
\seealso{
  \code{\link{updatePredictorPerObservationJointly}} for applying a shared per-observation predictor update across several samplers at once.

  \code{\link{samplePriorPredictive}} for repeated \code{sampleTreesFromPrior}/\code{sampleNodeParametersFromPrior}/\code{predict} draws on a private sampler, for calibrating priors before fitting.

  \code{\link{dbarts-embedding}} for which channel an outer block writes, which fit it reads back, and whose prior is in force after a swap, and \code{vignette("dbarts-as-a-component", package = "dbarts")} for six worked compositions built from these methods.
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
