## Resolves a sampler specification - the (control, model, data) triple and the
## family token - from an already-materialized response. Everything upstream of
## this (formula/matrix dispatch, survival ingestion, the dbartsData build) is
## the entry point's business; everything here is shared, so dbarts() and
## dbartsSpec() can never resolve a family two ways.
##
## The prior parse stays call-shaped: the prior vocabulary is NSE (a bare
## normal(chi(1.5)) must resolve in dbarts's vocabulary no matter what the
## caller has attached), so the entry point hands over its own match.call() and
## formals, and parsePriors evaluates the argument expressions in evalEnv.
## `matchedCall` must therefore be the CALLER's, not this function's.
##
## Binary latent-variable families: probit and logistic draw a fixed-unit-
## scale latent index rather than an ordinary continuous response, which is
## what control@binary, the weight policy, and the resid.prior override below
## all key off. Shared so no entry point's own family gate can drift from
## this one - xbart once carried the weaker `family != "gaussian"`.
isBinaryFamily <- function(family) {
  family %in% c("probit", "logistic")
}

## Wraps estimateSigmaFromLinearModel so every caller needing a starting sigma
## estimate raises the same failure message instead of a bare lm() error.
estimateStartingSigma <- function(data) {
  tryResult <- tryCatch(
    estimateSigmaFromLinearModel(data),
    error = function(e) e
  )
  if (inherits(tryResult, "error")) {
    stop("unable to obtain a starting estimate of sigma; provide one instead")
  }
  tryResult
}

## The binary/ordinal/nbinom weight policy, shared by every entry point that
## can reach these families: a probit has no tractable weighted latent-
## variable form and is refused, except that weights identically 1 are the
## unweighted likelihood and are treated as absent (SuperLearner-style
## callers pass obsWeights = rep(1, n) unconditionally); a logistic model
## treats weights as observation counts (its Polya-Gamma latent is a sum of
## per-copy draws), so they must be positive integers; ordinal refuses for
## probit's reason, with the same all-ones courtesy; nbinom refuses outright
## (exposure belongs in the offset, not observation replication). Gaussian
## weights are unrestricted and reach here as a no-op.
enforceWeightPolicy <- function(data, family) {
  if (is.null(data@weights)) {
    return(data)
  }
  if (family == "probit") {
    if (all(data@weights == 1)) {
      data@weights <- NULL
    } else {
      stop(
        "probit models do not support weights; fit integer count weights ",
        "with family = \"logistic\", or model continuous weights' latents ",
        "directly"
      )
    }
  } else if (family == "logistic") {
    w <- data@weights
    if (anyNA(w) || any(w <= 0) || any(w != round(w))) {
      stop(
        "logistic weights are observation counts and must be positive ",
        "integers; drop zero-count rows, and use a gaussian model for ",
        "continuous weights"
      )
    }
  } else if (family == "ordinal") {
    if (all(data@weights == 1)) {
      data@weights <- NULL
    } else {
      stop(
        "ordinal models do not support weights: a weighted truncated-normal ",
        "latent likelihood is not a coherent model"
      )
    }
  } else if (family == "nbinom") {
    if (all(data@weights == 1)) {
      data@weights <- NULL
    } else {
      stop(
        "nbinom (count) models do not support weights: exposure belongs in ",
        "the offset as a log-exposure term"
      )
    }
  } else if (family == "multinomial") {
    if (all(data@weights == 1)) {
      data@weights <- NULL
    } else {
      stop(
        "multinomial (softmax) models do not support weights: an integer ",
        "weight is already row-wise replication in the count response, and a ",
        "non-integer one has no exact augmentation sampler"
      )
    }
  }
  data
}

## survivalStatus and hazardPeriods carry the two survival markers the entry
## point resolved from the raw response; both are NULL for every other family.
resolveSamplerSpec <- function(
  matchedCall,
  callFormals,
  control,
  data,
  family,
  requestedFamily,
  dispersion,
  proposal.probs,
  monotone,
  interactions,
  blocks,
  variance,
  survivalStatus,
  hazardPeriods,
  bases,
  forests,
  evalEnv
) {
  # a factor/logical/character response declares a classification model. The
  # single-forest engine here fits only the 2-level (probit) case; 3+ levels
  # are multinomial, which only bart2(family = "multinomial") implements. A
  # numeric response takes the historic 0/1-vs-continuous path unchanged.
  # dbarts() is also reached anonymously through bart(), which never sets an
  # explicit family; see resolveClassificationFamily's doc comment for why
  # its auto-branch message lists every single-forest entry point instead of
  # naming itself. probit/logistic on a 2-level categorical response proceed
  # as binary.
  family <- resolveClassificationFamily(
    data,
    family,
    "dbarts()/rbart_vi()/bart()/xbart",
    c("gaussian", "aft", "nbinom"),
    splitMultinomialMessage = TRUE,
    allowOrdinal = TRUE
  )
  # multinomial (K-forest softmax): the response is
  # the n x K count matrix on the data object, so the family is DECLARED by the
  # slot and not inferred from any response shape - data@y is that matrix's
  # trials vector. A counts-carrying object resolves to it from "auto" without
  # an announcement (unlike the ordered-factor auto-dispatch, there is nothing
  # ambiguous to report: the slot is itself the declaration), and every other
  # explicit family is refused rather than silently fitting the trials.
  counts <- dataCounts(data)
  if (!is.null(counts) && identical(family, "auto")) {
    family <- "multinomial"
  }
  if (identical(family, "multinomial") && is.null(counts)) {
    stop(
      "family \"multinomial\" needs an n x K count-matrix response; build ",
      "the data with dbartsData(counts = ), or pass a factor or count ",
      "matrix as the response to dbarts()"
    )
  }
  if (!identical(family, "multinomial") && !is.null(counts)) {
    stop(
      "the data carry an n x K count matrix, which only family ",
      "\"multinomial\" fits; drop 'counts' to fit family \"",
      family,
      "\" against their trials"
    )
  }
  if (identical(family, "multinomial")) {
    # the counts ARE the response and the bridge reads K off their column
    # count, so nothing is recoded and no control attribute carries the
    # category count; data@y is already the trials
    NULL
  } else if (identical(family, "ordinal")) {
    # ordinal (cumulative probit): a single-forest
    # fixed-unit-scale model like probit, but K-level. Recode the response to
    # the 1-based category codes the engine reads, and attach K on the control
    # attribute the bridge reads to select OrdinalResponse (the
    # bartcore.survival precedent below). The resolved ordered levels ride the
    # data object for the round-trip.
    ordinal <- resolveOrdinalResponse(data)
    data@y <- ordinal$y
    data@response.levels <- ordinal$levels
    attr(control, "bartcore.n.categories") <- ordinal$K
  } else if (identical(family, "nbinom")) {
    # negative-binomial counts: the
    # count response has no unambiguous class, so "nbinom" is never auto - only
    # explicit. y must be a non-negative integer count (the NB pmf has zero mass
    # off the integers and the grid kernel's count histogram presumes integer
    # y), validated here beside the binary 0/1 test.
    y <- data@y
    if (anyNA(y) || any(y < 0) || any(y != round(y))) {
      stop("family \"nbinom\" requires a non-negative integer (count) response")
    }
    # the dispersion r: NA (the default) estimates it on the capped integer grid;
    # a supplied value FIXES it and must be a positive integer (v1 ships the
    # exact integer envelope, section 2). The C bridge reads the resolved spec
    # off the control attribute the sampler build attaches below: a positive
    # value fixes r, a non-positive value estimates it on the grid.
    dispersionSpec <- resolveDispersion(dispersion)
    attr(control, "bartcore.dispersion") <- dispersionSpec
  } else if (data@response.type == "numeric") {
    uniqueResponses <- unique(data@y)
    responseIsBinary <- length(uniqueResponses) == 2 &&
      all(sort(uniqueResponses) == c(0, 1))
    if (family == "auto") {
      family <- if (responseIsBinary) "probit" else "gaussian"
    } else if (family != "gaussian" && family != "aft" && !responseIsBinary) {
      # gaussian on a 0/1 response is a legitimate request; the binary
      # families need latent-variable coding. aft fits continuous log-times.
      stop("family \"", family, "\" requires a response coded 0/1")
    }
  }
  # aft draws sigma and rescales like gaussian; only the binary families are
  # latent-variable models on a fixed unit scale
  control@binary <- isBinaryFamily(family)
  # ordinal (cumulative probit) shares probit's fixed unit latent scale - sigma
  # fixed at 1, resid.prior fixed(1), no sigma estimate, node.scale 3.0 - but is
  # NOT binary: the bridge selects it by the bartcore.n.categories attribute
  # (not control@binary), and it reports K category levels. nbinom (counts) is
  # likewise a fixed-unit-scale family (sigma fixed at 1, the log-odds latent
  # entering kappa directly), selected by the bartcore.dispersion attribute.
  # fixedUnitScale covers all
  # three families wherever the unit-scale handling matters.
  # multinomial (softmax) is the fourth: its K category forests take their leaf
  # scale from the softmax calibration map's own anchor, there is no residual
  # scale to draw, and a single-trial count response has an identically
  # constant trials vector that estimateStartingSigma would fit a degenerate
  # sigma against
  fixedUnitScale <- control@binary ||
    identical(family, "ordinal") ||
    identical(family, "nbinom") ||
    identical(family, "multinomial")

  # binary/ordinal/nbinom weight policy, enforced here in the R layer (the
  # bridge keeps the same checks as a backstop for direct-API consumers) -
  # see enforceWeightPolicy's own doc comment for the rule each family
  # follows. Gaussian weights, including a gaussian fit of a 0/1 response,
  # are unrestricted and pass through untouched.
  data <- enforceWeightPolicy(data, family)

  if (is.na(data@sigma) && !fixedUnitScale) {
    data@sigma <- estimateStartingSigma(data)
  }

  # bart will passthrough with offset == something no matter what, which we
  # can NULL out; a latent-scale fixed-unit family (binary, ordinal) keeps its
  # zero offset as the meaningful reference, as probit always has
  if (!fixedUnitScale && !is.null(data@offset) && all(data@offset == 0.0)) {
    data@offset <- NULL
  }
  if (
    !fixedUnitScale &&
      !is.null(data@offset.test) &&
      all(data@offset.test == 0.0)
  ) {
    data@offset.test <- NULL
  }
  # the softmax is invariant to a common per-observation shift, so a flat
  # offset points exactly along its null direction: an all-zero passthrough
  # names nothing and is dropped, and anything else is refused by the channel
  # that does mean something. The C bridge keeps the same refusal as a backstop.
  if (identical(family, "multinomial")) {
    if (!is.null(data@offset) && all(data@offset == 0.0)) {
      data@offset <- NULL
    }
    if (!is.null(data@offset.test) && all(data@offset.test == 0.0)) {
      data@offset.test <- NULL
    }
    if (!is.null(data@offset)) {
      refuseFlatOffsetOnMultinomial(data@offset)
    }
    if (!is.null(data@offset.test)) {
      refuseFlatOffsetOnMultinomial(data@offset.test, "offset.test")
    }
  }

  # the multi-forest declaration: forests = list(forest(),
  # forest(basis = ~ factor(z))) names the ensembles the mean is a weighted sum
  # of, and every knob is per forest. The FIRST forest's structural knobs are
  # this fit's own - its tree count is control@n.trees and its structure prior
  # is tree.prior - so they are applied to those here, before anything reads
  # them; its interactions/blocks are the top-level arguments of those names,
  # which is why declaring both refuses. The rest ride the forests control
  # attribute below.
  # PER FOREST, not one flag for the model: forest f is excused from declaring
  # a basis only when one reaches it some other way, the dbartsData(bases = )
  # route being the supported one. On the fitting path the declarations have
  # already been expanded onto the data object, so this reads them back and a
  # forest that declared none is refused by name.
  declaredBases <- if (is.null(bases)) data@bases else bases
  forestSpec <- resolveForests(
    forests,
    interactions,
    blocks,
    hasBasis = if (is.null(declaredBases)) {
      logical(0L)
    } else {
      !vapply(declaredBases, is.null, logical(1L))
    }
  )
  firstForest <- if (is.null(forestSpec)) NULL else forestSpec[[1L]]
  if (!is.null(firstForest$n.trees)) {
    control@n.trees <- firstForest$n.trees
  }
  if (!is.null(firstForest$interactions)) {
    interactions <- firstForest$interactions
  }
  if (!is.null(firstForest$blocks)) {
    blocks <- firstForest$blocks
  }

  # resolve the monotone spec here, where its argument is a forced value (a
  # wrapper forwarding it through ... would otherwise reach parsePriors as an
  # unevaluated ...-reference); the resolved direction vector is injected below
  monotoneDirections <- resolveMonotone(monotone, data)

  parsePriorsCall <- redirectCall(matchedCall, quoteInNamespace(parsePriors))
  parsePriorsCall <- setDefaultsFromFormals(
    parsePriorsCall,
    callFormals,
    "tree.prior",
    "node.prior",
    "resid.prior",
    "resid.dist"
  )
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$monotone <- monotoneDirections
  # a multi-forest fit is the one whose forests carry amplitude bases; the
  # calibration map pins every forest's k, so the binary chi-k default is
  # redirected to the fixed 2 rather than refused below
  parsePriorsCall$multiForest <- !is.null(declaredBases)
  parsePriorsCall$parentEnv <- evalEnv

  if (fixedUnitScale) {
    if (any(names(parsePriorsCall) == "resid.prior")) {
      parsePriorsCall[[which(
        names(parsePriorsCall) == "resid.prior"
      )]] <- quote(fixed(1))
    } else {
      parsePriorsCall[[length(parsePriorsCall) + 1L]] <- quote(fixed(1))
      names(parsePriorsCall)[length(parsePriorsCall)] <- "resid.prior"
    }
  }
  priors <- eval(parsePriorsCall)

  # the first forest's structure prior is this fit's tree.prior, so a knob
  # declared on it restates that prior's half rather than adding a second one
  if (!is.null(firstForest$base)) {
    priors$tree.prior@base <- firstForest$base
  }
  if (!is.null(firstForest$power)) {
    priors$tree.prior@power <- firstForest$power
  }

  # A monotone constraint restricts the forest to birth/death proposals (change
  # and swap would need a > 2-D constrained integral, monotone.md section 5): a
  # defaulted proposal.probs is forced to birth/death-only, an explicit
  # non-default one conflicts and errors.
  if (!is.null(monotoneDirections)) {
    defaultProbs <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
    if (
      !is.null(proposal.probs) &&
        !isTRUE(all.equal(proposal.probs[names(defaultProbs)], defaultProbs))
    ) {
      stop(
        "'monotone' forces birth/death-only proposals; a non-default ",
        "'proposal.probs' cannot be honored under the constraint"
      )
    }
    proposal.probs <- c(birth_death = 1, swap = 0, change = 0, birth = 0.5)
  }

  model <- newValidated(
    "dbartsModel",
    priors$tree.prior,
    priors$node.prior,
    priors$node.hyperprior,
    priors$resid.prior,
    proposal.probs = proposal.probs,
    family = family,
    # a named calibration (node.prior's scale = / sd =) overrides the
    # family default below in the engine, which converts it out of response
    # units against the transform; NA leaves that default in force
    prior.scale = resolvePriorScale(priors$node.prior, priors$node.hyperprior),
    node.scale = defaultNodeScale(family)
  )

  # Student-t residuals: only a continuous
  # gaussian response carries them (the binary families and aft have their own
  # latent scale), refused here R-side to match the C bridge's backstop. The
  # resolved degrees of freedom ride the model's resid.df attribute the bridge
  # reads - the bartcore.survival precedent above - absent for the Gaussian law.
  if (is(priors$resid.dist, "dbartsStudentDist") && family != "gaussian") {
    stop(
      "student residuals require a continuous gaussian response; family \"",
      requestedFamily,
      "\" has its own fixed error scale"
    )
  }
  residDf <- residDistDf(priors$resid.dist)
  if (!is.null(residDf)) {
    attr(model, "resid.df") <- residDf
  }

  # the resolved per-column monotone directions ride the model attribute the C
  # bridge reads into SamplerOptions.monotoneDirections (the resid.df precedent)
  if (!is.null(monotoneDirections)) {
    attr(model, "monotone") <- monotoneDirections
  }

  # the resolved per-forest interaction constraint (max-order cap + forbidden
  # co-occurrence pairs) rides two model attributes the C bridge reads into
  # SamplerOptions (the monotone precedent).
  # Absent when no interactions() prior is supplied, so the availability path is
  # byte-for-byte unchanged.
  interactionSpec <- resolveInteractions(interactions, data)
  if (!is.null(interactionSpec)) {
    attr(model, "interaction.max.order") <- interactionSpec$max.order
    attr(model, "interaction.forbidden") <- interactionSpec$forbidden
  }

  # the resolved per-forest block-additive constraint (variant A): each whole
  # tree is confined to one declared group, so the ensemble is exactly
  # f = sum_G f_G. Rides two model attributes the C bridge reads (the
  # interactions precedent); absent when no blocks() prior is supplied, so the
  # path is byte-for-byte unchanged. The partition covers the full design.
  blockSpec <- resolveBlocks(blocks, data, control@n.trees)
  if (!is.null(blockSpec)) {
    attr(model, "block.of.column") <- blockSpec$block.of.column
    attr(model, "block.tree.counts") <- blockSpec$block.tree.counts
  }

  # The K category forests are built from the softmax calibration map and the
  # CONSTANT-leaf instantiation only: a monotone
  # constraint or a non-constant leaf selects an instantiation the multinomial
  # factory does not build, a DART prior and a drawn k are unadjudicated
  # against the map's fixed anchor, the map owns every leaf scale so a named
  # prior.scale has nowhere to land, and the grouped decorator wraps ONE
  # response model, which a K-forest blend is not. Every one of these would
  # otherwise be dropped in silence, changing the fitted model without a word;
  # name each one instead. The bridge keeps its own backstops for the callers
  # that reach it without this layer.
  if (identical(family, "multinomial")) {
    unsupportedMultinomial <- c(
      "a DART tree prior" = is(priors$tree.prior, "dbartsDartPrior"),
      "'split.probs'" = length(priors$tree.prior@splitProbabilities) > 0L,
      "'monotone'" = !is.null(monotoneDirections),
      "a linear node prior" = is(priors$node.prior, "dbartsLinearPrior"),
      "a Gaussian-process node prior" = is(priors$node.prior, "dbartsGPPrior"),
      "a 'k' hyperprior" = is(priors$node.hyperprior, "dbartsChiHyperprior"),
      "a named 'prior.scale'" = !is.na(model@prior.scale),
      "grouped random effects" = !is.null(attr(control, "bartcore.groups")),
      "storage = \"single\"" = identical(control@storage, "single")
    )
    if (any(unsupportedMultinomial)) {
      stop(
        "a multinomial (softmax) model does not support ",
        paste0(
          names(unsupportedMultinomial)[unsupportedMultinomial],
          collapse = ", "
        ),
        "; drop it or fit a single-forest model"
      )
    }
  }

  # the AFT survival family reads its per-observation status off this control
  # attribute (the bartcore.groups precedent); the C bridge validates it
  if (!is.null(survivalStatus)) {
    if (length(survivalStatus) != length(data@y)) {
      stop("survival status must have length ", length(data@y))
    }
    attr(control, "bartcore.survival") <- survivalStatus
  }
  # the discrete-time hazard marker: the
  # period grid, parked here for packageBartResults to read into $periods. The
  # C bridge never reads this attribute (unlike bartcore.survival), so a hazard
  # fit's draw stream is byte-identical to the by-hand binary fit's - the
  # reduction gate (benchmarks/R/hazard-reduction.R).
  if (!is.null(hazardPeriods)) {
    attr(control, "bartcore.hazard.periods") <- hazardPeriods
  }

  # the heteroscedastic variance forest: a
  # `variance` selector installs a second forest modeling s^2(x); its config
  # rides the control attribute the C bridge reads. Gaussian + constant leaf
  # only (the C factory refuses otherwise; a friendly R check for the family).
  # `variance` accepts either the plain shorthand (NULL/FALSE/TRUE/formula/
  # character/index) or a varianceForest() object: its `vars`
  # slot routes through the SAME resolveVarianceColumns the shorthand uses -
  # one selector vocabulary - and its n.trees/base/power knobs land exactly
  # where the removed flat n.trees.variance/power.variance/base.variance
  # formals' values used to, byte-identically.
  varianceSpec <- if (inherits(variance, "dbartsVarianceForest")) {
    variance
  } else {
    NULL
  }
  # vars = NULL on the OBJECT means every column (resolveVarianceColumns'
  # TRUE reading), unlike variance = NULL's "no variance forest" - the two
  # NULLs mean opposite things, so the object's own NULL is translated to
  # TRUE before it reaches the shared resolver rather than passed through
  varianceSelector <- if (is.null(varianceSpec)) {
    variance
  } else if (is.null(varianceSpec$vars)) {
    TRUE
  } else {
    varianceSpec$vars
  }
  varianceColumns <- resolveVarianceColumns(varianceSelector, data)
  if (!is.null(varianceColumns)) {
    if (family != "gaussian") {
      stop(
        "a variance forest requires family = \"gaussian\"; family \"",
        requestedFamily,
        "\" routes precision through its own latent channel instead"
      )
    }
    # the scale-mixture reweighting and the variance forest's weight-channel
    # routing are unadjudicated together; refuse rather than silently fit
    # an uncomposed model
    if (!is.null(residDf)) {
      stop(
        "a variance forest does not support Student-t residuals: the two ",
        "are not yet shown to compose"
      )
    }
    if (!is.null(monotoneDirections)) {
      stop("a variance forest is not supported with monotone constraints")
    }
    allColumns <- setequal(varianceColumns, seq_len(ncol(data@x)))
    varianceNTrees <- if (is.null(varianceSpec)) NULL else varianceSpec$n.trees
    varianceBase <- if (is.null(varianceSpec)) NULL else varianceSpec$base
    variancePower <- if (is.null(varianceSpec)) NULL else varianceSpec$power
    attr(control, "bartcore.variance") <- list(
      n.trees = coerceOrError(
        if (is.null(varianceNTrees)) 40L else varianceNTrees,
        "integer"
      ),
      base = if (is.null(varianceBase)) {
        model@tree.prior@base
      } else {
        as.double(varianceBase)
      },
      power = if (is.null(variancePower)) {
        model@tree.prior@power
      } else {
        as.double(variancePower)
      },
      columns = if (allColumns) NULL else as.integer(varianceColumns)
    )
  }

  # the Bayesian causal forest: a second forest with a
  # two-level factor basis selects the model y = a mu(x) + b_z tau(x) + eps.
  # The 0/1 column that basis expands to is conditioning DATA and rides the
  # data object beside the weights it mirrors; the second forest's
  # CONFIGURATION - its tree count and structure prior, its column mask, the
  # amplitude scales and the per-forest constraints - rides the control
  # attribute the C bridge reads, exactly as the variance forest's does above.
  # The bridge cross-checks the two halves in both directions, so a stripped
  # attribute is a loud error, never a silent single-forest fit.
  if (!is.null(bases)) {
    # the only caller supplies these as forests' bases, so a refusal names that
    data@bases <- validateForestBases(
      bases,
      length(data@y),
      argument = "basis"
    )
  }
  if (!is.null(data@bases)) {
    # the RESOLVED forest count, which is the data object's own bases or the
    # declaration that replaced them just above, and which the refusal below
    # names the source of: a length-1 declaration over a data object already
    # carrying two bases resolves to one, and telling that caller they wrote
    # one basis would be false
    numForests <- length(data@bases)
    # K = 1 is not a shipped configuration. Both creation routes reach here -
    # the dbartsData(bases = ) one and the forests = one, whose declarations
    # forestBasisDeclarations now carries down at any length - so this is the
    # single site the refusal is owed at. What a lone amplitude forest is
    # missing is the second ensemble its amplitudes would distinguish it from;
    # reaching it as a configuration is VCBART's shape (the general per-forest
    # multiplier door), which owes acceptance evidence of its own before it
    # could be opened.
    if (numForests < 2L) {
      stop(
        "a multi-forest model needs at least two forests, and ",
        if (is.null(bases)) {
          "the data object carries "
        } else {
          "this call's 'basis' declarations resolve to "
        },
        numForests,
        "; the amplitudes multiplying a forest's basis are what distinguish ",
        "it from another forest, so one forest alone has nothing to be ",
        "distinguished from. Declare a second - forests = list(forest(), ",
        "forest(basis = ...)) is the two-forest spelling - or drop the basis ",
        "for a single-forest model"
      )
    }
    # the families the calibration map has a latent scale to state its node
    # scales against, and whose own parameter block is shown to interleave with
    # the amplitude block. A fixed error scale is what makes the binary
    # families work here rather than a reason they cannot: the combined index
    # is stated in the link's own units and sigma is pinned there.
    if (family %not_in% c("gaussian", "probit", "logistic")) {
      stop(
        "a treatment forest does not support family \"",
        family,
        "\": ",
        switch(
          family,
          aft = paste0(
            "it draws sigma, which the calibration map pins, and its ",
            "censoring status reaches no multi-forest creation path"
          ),
          ordinal = paste0(
            "its cutpoint block is not shown to interleave with the ",
            "amplitude block"
          ),
          nbinom = paste0(
            "its dispersion block is not shown to interleave with the ",
            "amplitude block"
          ),
          multinomial = paste0(
            "its forests are its categories, and the softmax blend that ",
            "combines them is not an amplitude coupling"
          ),
          "the calibration map states no scale for it"
        )
      )
    }
    # The amplitude chain builds every forest from its own calibration map (fixed
    # k = 1, leaf scales from the family's own latent scale) and reads neither
    # the DART machinery, the split probabilities, the monotone directions, a
    # non-constant leaf, the grouped decorator, a variance forest, an fp32
    # residual, a per-column cut cap, nor the Student-t error law. Every one of
    # those would otherwise be dropped in silence, changing the fitted model
    # without a word; name each one instead.
    defaultProbs <- defaultProposalProbs
    unsupported <- c(
      "a DART tree prior" = is(priors$tree.prior, "dbartsDartPrior"),
      "'split.probs'" = length(priors$tree.prior@splitProbabilities) > 0L,
      "'monotone'" = !is.null(monotoneDirections),
      "a linear node prior" = is(priors$node.prior, "dbartsLinearPrior"),
      "a Gaussian-process node prior" = is(priors$node.prior, "dbartsGPPrior"),
      "a 'k' hyperprior" = is(priors$node.hyperprior, "dbartsChiHyperprior"),
      "a non-default 'k'" = is(
        priors$node.hyperprior,
        "dbartsFixedHyperprior"
      ) &&
        priors$node.hyperprior@k != 2.0,
      # "differs from the family default", which is what the gaussian-only
      # literal 0.5 this used to compare against always meant
      "a non-default 'node.scale'" = model@node.scale !=
        defaultNodeScale(family),
      # the calibration map fixes every forest's leaf scale from the family's
      # own latent scale, so a named prior.scale has nowhere to land and the
      # node.scale gate above does not fire on it
      "a named 'prior.scale'" = !is.na(model@prior.scale),
      # a monotone constraint rewrites proposal.probs above, so only an
      # unconstrained fit's is the caller's own
      "a non-default 'proposal.probs'" = is.null(monotoneDirections) &&
        !isTRUE(all.equal(proposal.probs[names(defaultProbs)], defaultProbs)),
      "Student-t residuals" = !is.null(residDf),
      "grouped random effects" = !is.null(attr(control, "bartcore.groups")),
      "'variance'" = !is.null(varianceColumns),
      "storage = \"single\"" = identical(control@storage, "single"),
      "per-column 'n.cuts'" = length(unique(data@n.cuts)) > 1L,
      "test predictors" = !is.null(data@x.test)
    )
    if (any(unsupported)) {
      stop(
        "a treatment forest does not support ",
        paste0(names(unsupported)[unsupported], collapse = ", "),
        "; drop it or fit a single-forest model"
      )
    }
    # every forest resolves its own knobs; a data object carrying bases with no
    # forests = declaration at all resolves to the same defaults, which is what
    # keeps the dbartsData(bases = ) route a supported one
    specs <- if (is.null(forestSpec)) {
      rep(list(NULL), numForests)
    } else {
      forestSpec
    }
    if (length(specs) > numForests) {
      stop(
        "'forests' declares ",
        length(specs),
        " forests but the data carry ",
        numForests,
        " bases"
      )
    }
    # a declaration reaching only the first forests leaves the rest at the
    # engine's defaults, which is what keeps a bases-only data object
    # configurable one forest at a time
    if (length(specs) < numForests) {
      specs <- c(specs, rep(list(NULL), numForests - length(specs)))
    }
    hasBasis <- !vapply(data@bases, is.null, logical(1L))
    treeCounts <- vapply(
      seq_len(numForests),
      function(index) {
        if (index == 1L) {
          control@n.trees
        } else if (is.null(specs[[index]]$n.trees)) {
          50L
        } else {
          as.integer(specs[[index]]$n.trees)
        }
      },
      integer(1L)
    )
    forestColumns <- lapply(
      seq_len(numForests),
      function(index) resolveModerators(specs[[index]]$vars, data, "vars")
    )
    # a declaration's own per-forest names (forests = list(prognostic = ...,
    # treatment = ...)), distinct from the packaged channels' own
    # forest1..forestK vocabulary; NULL when the caller named none, "" for a
    # forest left unnamed among named siblings or reached only through
    # data@bases
    forestLabels <- names(forests)
    if (!is.null(forestLabels)) {
      length(forestLabels) <- numForests
      forestLabels[is.na(forestLabels)] <- ""
    }
    attr(control, "bartcore.forests") <- list(
      # one length-8 numeric per forest; the family selects the basis-free
      # channel's default median and the count the K-aware node scale factor
      params = forestParams(specs, hasBasis, family),
      # resolved 1-based column indices per forest, or NULL for unrestricted
      vars = forestColumns,
      # the declaration's own names, or NULL; see forestLabels above
      labels = forestLabels,
      # the first forest takes the fit's own interactions()/blocks() arguments,
      # already resolved above; the rest take their own, resolved against the
      # columns they may split on
      interactions = c(
        list(interactionSpec),
        lapply(
          specs[-1L],
          function(spec) resolveInteractions(spec$interactions, data)
        )
      ),
      blocks = c(
        list(blockSpec),
        lapply(
          seq_len(numForests)[-1L],
          function(index) {
            resolveBlocks(
              specs[[index]]$blocks,
              data,
              treeCounts[index],
              availableColumns = forestColumns[[index]]
            )
          }
        )
      )
    )
  }

  namedList(control, model, data, family)
}

## The exported consumer surface: resolves
## a specification without constructing a sampler, for a LinkingTo: dbarts
## consumer that holds its sampler C-side through dbarts.h and supplies its own
## design matrix. dbartsData() (already exported) builds the response half; this
## builds the rest, so no consumer has to reach into an internal to reach a
## feature that is otherwise complete.
dbartsSpec <- function(
  data,
  control = dbarts::dbartsControl(),
  tree.prior = cgm,
  node.prior = normal,
  resid.prior = chisq,
  resid.dist = gaussian,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  forests = NULL,
  sigma = NA_real_,
  seed = NA_integer_,
  family = c(
    "auto",
    "gaussian",
    "probit",
    "logistic",
    "aft",
    "multinomial",
    "ordinal",
    "nbinom"
  ),
  dispersion = NA_real_,
  survival = NULL,
  parentEnv = parent.frame()
) {
  matchedCall <- match.call()

  if (!inherits(data, "dbartsData")) {
    stop("'data' must be a dbartsData object; see ?dbartsData")
  }
  if (!inherits(control, "dbartsControl")) {
    stop("'control' must be a dbartsControl object; see ?dbartsControl")
  }
  family <- match.arg(family)

  # the survival status is the one piece of response ingestion this surface
  # cannot do for the caller: dbarts() reads it off a Surv or two-column
  # response, and a consumer supplying log-times directly supplies it here
  if (!is.null(survival)) {
    if (!identical(family, "aft")) {
      stop("'survival' status is only used by family \"aft\"")
    }
    survival <- as.double(survival)
    if (anyNA(survival) || any(survival != 0.0 & survival != 1.0)) {
      stop("survival status must be 0 (censored) or 1 (event)")
    }
  } else if (identical(family, "aft")) {
    stop(
      "family \"aft\" needs a 'survival' status vector (1 = event, ",
      "0 = right-censored) alongside a log-time response"
    )
  }

  seed <- coerceOrError(seed, "integer")
  if (!is.na(seed)) {
    control@seed <- seed
  }

  # the control owns the cut-point count, as it does inside dbarts(), but a data
  # object already carrying resolved per-column counts keeps them - a consumer
  # that set them deliberately is not silently overridden
  if (length(data@n.cuts) != ncol(data@x) || anyNA(data@n.cuts)) {
    data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))
  }
  # an explicit sigma overrides whatever the data carries; the default leaves it
  # alone, so a consumer's own starting estimate survives (an NA is estimated
  # during resolution, exactly as for dbarts())
  if (!missing(sigma)) {
    data@sigma <- coerceOrError(sigma, "numeric")
  }

  # this surface does no data ingestion of its own, so a declared basis is
  # evaluated here and reaches data@bases through the same validation
  # dbartsData() applies on the fitting path; the caller's data object has
  # already had its own 'subset' applied
  basisDeclarations <- forestBasisDeclarations(forests)
  basis <- if (is.null(basisDeclarations)) {
    NULL
  } else {
    expanded <- lapply(
      basisDeclarations,
      function(declaration) {
        expandForestBasis(evaluateForestBasis(declaration))
      }
    )
    # as in dbarts(): a list in which no forest declares a basis names no
    # multi-forest model, and falls through to resolveForests' refusal
    if (any(!vapply(expanded, is.null, logical(1L)))) expanded else NULL
  }

  resolveSamplerSpec(
    matchedCall,
    formals(dbartsSpec),
    control,
    data,
    family,
    requestedFamily = family,
    dispersion = dispersion,
    proposal.probs = proposal.probs,
    monotone = monotone,
    interactions = interactions,
    blocks = blocks,
    variance = variance,
    survivalStatus = survival,
    hazardPeriods = NULL,
    bases = basis,
    forests = forests,
    evalEnv = parentEnv
  )
}
