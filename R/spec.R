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
## survivalStatus and hazardPeriods carry the two survival markers the entry
## point resolved from the raw response; both are NULL for every other family.
resolveSamplerSpec <- function(
  matchedCall,
  callFormals,
  control,
  data,
  family,
  dispersion,
  proposal.probs,
  monotone,
  interactions,
  blocks,
  variance,
  n.trees.variance,
  power.variance,
  base.variance,
  survivalStatus,
  hazardPeriods,
  treatment,
  moderators,
  treatmentForest,
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
  # ordinal (cumulative probit, docs/design/ordinal.md): a single-forest fixed-
  # unit-scale model like probit, but K-level. Recode the response to the
  # 1-based category codes the engine reads, and attach K on the control
  # attribute the bridge (C2) reads to select OrdinalResponse (the
  # bartcore.survival precedent below). The resolved ordered levels ride the
  # data object for the round-trip.
  if (identical(family, "ordinal")) {
    ordinal <- resolveOrdinalResponse(data)
    data@y <- ordinal$y
    data@response.levels <- ordinal$levels
    attr(control, "bartcore.n.categories") <- ordinal$K
  } else if (identical(family, "nbinom")) {
    # negative-binomial counts (docs/design/negative-binomial.md section 4): the
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
  control@binary <- family %in% c("probit", "logistic")
  # ordinal (cumulative probit) shares probit's fixed unit latent scale - sigma
  # fixed at 1, resid.prior fixed(1), no sigma estimate, node.scale 3.0 - but is
  # NOT binary: the bridge selects it by the bartcore.n.categories attribute
  # (not control@binary), and it reports K category levels. nbinom (counts) is
  # likewise a fixed-unit-scale family (sigma fixed at 1, the log-odds latent
  # entering kappa directly, docs/design/negative-binomial.md section 1),
  # selected by the bartcore.dispersion attribute. fixedUnitScale covers all
  # three families wherever the unit-scale handling matters.
  fixedUnitScale <- control@binary ||
    identical(family, "ordinal") ||
    identical(family, "nbinom")

  # binary weight policy, enforced here in the R layer (the bridge keeps the
  # same checks as a backstop for direct-API consumers): a probit has no
  # tractable weighted latent-variable form and is refused, except that
  # weights identically 1 are the unweighted likelihood and are treated as
  # absent (SuperLearner-style callers pass obsWeights = rep(1, n)
  # unconditionally); a logistic model treats weights as observation counts
  # (its Polya-Gamma latent is a sum of per-copy draws), so they must be
  # positive integers. Gaussian weights, including a gaussian fit of a 0/1
  # response, are unrestricted.
  if (!is.null(data@weights)) {
    if (family == "probit") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "probit models do not support weights: a weighted probit has no ",
          "tractable latent-variable form. Integer count weights can be ",
          "fit exactly with family = \"logistic\"; for continuous weights, ",
          "model the latents directly."
        )
      }
    }
    if (family == "logistic") {
      w <- data@weights
      if (anyNA(w) || any(w <= 0) || any(w != round(w))) {
        stop(
          "logistic weights are observation counts and must be positive ",
          "integers; drop zero-count rows, and use a gaussian model for ",
          "continuous weights."
        )
      }
    }
    # ordinal refuses weights for probit's reason (a weighted truncated-normal
    # latent likelihood is not a coherent model), with the same all-ones-are-
    # absent courtesy the probit path extends to SuperLearner-style callers
    if (family == "ordinal") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "ordinal models do not support weights: a weighted truncated-normal ",
          "latent likelihood is not a coherent model."
        )
      }
    }
    # nbinom refuses weights in v1 (docs/design/negative-binomial.md section 4):
    # the usual count "weight" is exposure, which belongs in the offset as a
    # log-exposure term, not in observation replication. The R refusal mirrors
    # the bridge's backstop, with the probit all-ones-are-absent courtesy.
    if (family == "nbinom") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "nbinom (count) models do not support weights: exposure belongs in ",
          "the offset as a log-exposure term."
        )
      }
    }
  }

  if (is.na(data@sigma) && !fixedUnitScale) {
    tryResult <- tryCatch(
      data@sigma <- estimateSigmaFromLinearModel(data),
      error = function(e) e
    )
    if (inherits(tryResult, "error")) {
      stop("unable to obtain a starting estimate of sigma; provide one instead")
    }
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

  # A monotone constraint restricts the forest to birth/death proposals (change
  # and swap would need a > 2-D constrained integral, monotone.md section 5): a
  # defaulted proposal.probs is forced to birth/death-only, an explicit
  # non-default one conflicts and errors.
  if (!is.null(monotoneDirections)) {
    defaultProbs <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
    if (!isTRUE(all.equal(proposal.probs[names(defaultProbs)], defaultProbs))) {
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
    # a named calibration (node.prior's scale = / sd =) overrides the switch
    # below in the engine, which converts it out of response units against the
    # transform; NA leaves the family-keyed internal default in force
    prior.scale = resolvePriorScale(priors$node.prior, priors$node.hyperprior),
    # the logistic scale is probit's default widened by the logistic
    # latent's standard deviation, pi / sqrt(3)
    node.scale = switch(
      family,
      gaussian = 0.5,
      aft = 0.5,
      probit = 3.0,
      # ordinal reuses probit's latent scale (docs/design/ordinal.md section 2,
      # scheme C: the K = 2 anchor is probit exactly)
      ordinal = 3.0,
      # nbinom's psi is a log-odds, so it reuses logistic's node.scale
      # (docs/design/negative-binomial.md section 1)
      nbinom = pi * sqrt(3.0),
      logistic = pi * sqrt(3.0)
    )
  )

  # Student-t residuals (docs/design/robust-errors.md): only a continuous
  # gaussian response carries them (the binary families and aft have their own
  # latent scale), refused here R-side to match the C bridge's backstop. The
  # resolved degrees of freedom ride the model's resid.df attribute the bridge
  # reads - the bartcore.survival precedent above - absent for the Gaussian law.
  if (is(priors$resid.dist, "dbartsStudentDist") && family != "gaussian") {
    stop(
      "student residuals require a continuous gaussian response; family \"",
      family,
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
  # co-occurrence pairs, docs/design/interaction-constraints.md) rides two model
  # attributes the C bridge reads into SamplerOptions (the monotone precedent).
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

  # the AFT survival family reads its per-observation status off this control
  # attribute (the bartcore.groups precedent); the C bridge validates it
  if (!is.null(survivalStatus)) {
    if (length(survivalStatus) != length(data@y)) {
      stop("survival status length does not match the response")
    }
    attr(control, "bartcore.survival") <- survivalStatus
  }
  # the discrete-time hazard marker (docs/design/survival.md section 4): the
  # period grid, parked here for packageBartResults to read into $periods. The
  # C bridge never reads this attribute (unlike bartcore.survival), so a hazard
  # fit's draw stream is byte-identical to the by-hand binary fit's - the
  # reduction gate (benchmarks/R/hazard-reduction.R).
  if (!is.null(hazardPeriods)) {
    attr(control, "bartcore.hazard.periods") <- hazardPeriods
  }

  # the heteroscedastic variance forest (docs/design/heteroscedastic.md): a
  # `variance` selector installs a second forest modeling s^2(x); its config
  # rides the control attribute the C bridge reads. Gaussian + constant leaf
  # only (the C factory refuses otherwise; a friendly R check for the family).
  varianceColumns <- resolveVarianceColumns(variance, data)
  if (!is.null(varianceColumns)) {
    if (family != "gaussian") {
      stop(
        "a variance forest requires family = \"gaussian\"; the latent families ",
        "own the precision channel it routes through"
      )
    }
    if (!is.null(monotoneDirections)) {
      stop("a variance forest is not supported with monotone constraints")
    }
    allColumns <- setequal(varianceColumns, seq_len(ncol(data@x)))
    attr(control, "bartcore.variance") <- list(
      n.trees = coerceOrError(n.trees.variance, "integer"),
      base = if (is.null(base.variance)) {
        model@tree.prior@base
      } else {
        as.double(base.variance)
      },
      power = if (is.null(power.variance)) {
        model@tree.prior@power
      } else {
        as.double(power.variance)
      },
      columns = if (allColumns) NULL else as.integer(varianceColumns)
    )
  }

  # the Bayesian causal forest (docs/design/bcf.md): a 0/1 `treatment` vector
  # selects the two-forest model y = a mu(x) + b_z tau(x) + eps. The vector is
  # conditioning DATA and rides the data object beside the weights it mirrors;
  # the treatment forest's CONFIGURATION - its tree count and structure prior,
  # the moderator column mask, the glue scales and the per-forest constraints -
  # rides the control attribute the C bridge reads, exactly as the variance
  # forest's does above. The bridge cross-checks the two halves in both
  # directions, so a stripped attribute is a loud error, never a silent
  # single-forest fit.
  if (!is.null(treatment)) {
    data@treatment <- validateTreatment(treatment, length(data@y))
  }
  if (!is.null(data@treatment)) {
    if (family != "gaussian") {
      stop(
        "a treatment forest requires a continuous (gaussian) response; ",
        "family \"",
        family,
        "\" has its own fixed error scale"
      )
    }
    # The BCF chain builds both forests from its own calibration map (fixed
    # k = 1, leaf scales from the response sd) and reads neither the DART
    # machinery, the split probabilities, the monotone directions, a
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
      "a non-default 'node.scale'" = model@node.scale != 0.5,
      # the calibration map fixes both forests' leaf scales from the response
      # sd, so a named prior.scale has nowhere to land and the node.scale gate
      # above does not fire on it
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
    tauSpec <- resolveTreatmentForest(treatmentForest)
    moderatorColumns <- resolveModerators(moderators, data)
    attr(control, "bartcore.bcf") <- list(
      params = as.double(c(
        tauSpec$n.trees,
        tauSpec$base,
        tauSpec$power,
        tauSpec$sd.control,
        tauSpec$sd.moderate,
        tauSpec$b.prior.variance,
        tauSpec$update.a,
        tauSpec$update.b
      )),
      # resolved 1-based moderator indices, or NULL for an unrestricted forest
      moderators = moderatorColumns,
      # the prognostic forest takes the fit's own interactions()/blocks()
      # arguments, already resolved above; the treatment forest takes the
      # constructor's, resolved against the columns it may split on
      mu.interactions = interactionSpec,
      tau.interactions = resolveInteractions(tauSpec$interactions, data),
      mu.blocks = blockSpec,
      tau.blocks = resolveBlocks(
        tauSpec$blocks,
        data,
        tauSpec$n.trees,
        availableColumns = moderatorColumns
      )
    )
  } else if (!is.null(moderators) || !is.null(treatmentForest)) {
    stop(
      "'moderators' and 'treatmentForest' configure the treatment forest a ",
      "'treatment' vector selects; supply one or drop them"
    )
  }

  namedList(control, model, data, family)
}

## The exported consumer surface (docs/design/consumer-spec-surface.md): resolves
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
  n.trees.variance = 40L,
  power.variance = NULL,
  base.variance = NULL,
  treatment = NULL,
  moderators = NULL,
  treatmentForest = NULL,
  sigma = NA_real_,
  seed = NA_integer_,
  family = c(
    "auto",
    "gaussian",
    "probit",
    "logistic",
    "aft",
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
    control@rngSeed <- seed
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

  resolveSamplerSpec(
    matchedCall,
    formals(dbartsSpec),
    control,
    data,
    family,
    dispersion = dispersion,
    proposal.probs = proposal.probs,
    monotone = monotone,
    interactions = interactions,
    blocks = blocks,
    variance = variance,
    n.trees.variance = n.trees.variance,
    power.variance = power.variance,
    base.variance = base.variance,
    survivalStatus = survival,
    hazardPeriods = NULL,
    treatment = treatment,
    moderators = moderators,
    treatmentForest = treatmentForest,
    evalEnv = parentEnv
  )
}
