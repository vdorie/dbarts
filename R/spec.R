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
## this one.
## The discrete-time hazard tokens, each remapped to its binary link before
## the engine sees it.
hazardFamilyTokens <- c("hazard", "hazard.probit", "hazard.logistic")

## Refuses a response a binary family cannot fit, saying what is wrong with
## it: a 0/1 response with one class, or one not coded 0/1. A hazard fit
## ('hazard', the token the caller gave) is a binary fit on person-period rows
## the caller never wrote, so its refusal speaks of subjects and events.
refuseNonBinaryResponse <- function(uniqueResponses, family, hazard = NULL) {
  singleClass <- length(uniqueResponses) == 1L &&
    uniqueResponses %in% c(0, 1)
  if (!is.null(hazard) && singleClass) {
    stop(
      "family \"",
      hazard,
      "\" needs ",
      if (uniqueResponses == 0) {
        "an event; every subject is censored"
      } else {
        paste0(
          "a period at risk without an event; every subject has its event ",
          "in the first period"
        )
      },
      call. = FALSE
    )
  }
  if (singleClass) {
    stop(
      "family \"",
      family,
      "\" requires a response with both classes; the response has a single ",
      "class",
      call. = FALSE
    )
  }
  stop(
    "family \"",
    family,
    "\" requires a response coded 0/1",
    if (family == "logistic") {
      " (family = binomial is the logit link, a logistic fit)"
    },
    call. = FALSE
  )
}

## Whether a coded response holds one class of a binary one: every value 0,
## or every value 1, whatever encoding (numeric, logical, factor, character)
## it was coded from.
responseHasSingleClass <- function(y) {
  values <- unique(y[!is.na(y)])
  length(values) == 1L && values %in% c(0, 1)
}

## A single-class response that a binary family will refuse: an explicit
## binary family, or "auto" on a categorical encoding, which resolves to one.
refusesSingleClass <- function(data, family) {
  responseHasSingleClass(data@y) &&
    (family %in%
      c("probit", "logistic", hazardFamilyTokens) ||
      (identical(family, "auto") && data@response.type != "numeric"))
}

## dbartsData warns of a response whose values are indistinguishable at
## double precision before any family is known, and a constant response is
## one. Where that response is a single class a binary family refuses, the
## refusal names the actual problem, so the warning, whose remedy is to
## rescale, is held back there and raised everywhere else. A formula hazard
## fit's response at this point is the log time standing in for the binary
## rows it expands to, which the warning does not describe, so it is held back
## there too.
withBinaryResponsePrecision <- function(family, expr) {
  held <- NULL
  data <- withCallingHandlers(
    expr,
    warning = function(w) {
      if (startsWith(conditionMessage(w), responsePrecisionWarningStem)) {
        held <<- w
        invokeRestart("muffleWarning")
      }
    }
  )
  if (
    !is.null(held) &&
      family %not_in% hazardFamilyTokens &&
      !(is(data, "dbartsData") && refusesSingleClass(data, family))
  ) {
    warning(held)
  }
  data
}

isBinaryFamily <- function(family) {
  family %in% c("probit", "logistic")
}

## The families whose leaf-scale k defaults to the chi(1.5, 2) hyperprior: the
## binary links, and nbinom, whose log-mean leaf prior has a fixed anchor and
## no residual scale to calibrate against (dec-B183). Every other family's k
## defaults to a fixed 2.
drawsLeafKByDefault <- function(family) {
  isBinaryFamily(family) || identical(family, "nbinom")
}

## The latent-variable families that carry no case weight at all but do
## implement the active-row mask. A weight vector of 0s and 1s there is
## membership rather than precision - the row leaves the likelihood, keeps its
## leaf occupancy, its latent and its fitted value - which is exactly the mask,
## so such a vector installs as one instead of being refused. Any other value
## is a weighted latent likelihood, which these families have no coherent form
## for, and stays refused.
isMaskedWeightFamily <- function(family) {
  family %in% c("probit", "ordinal", "nbinom")
}

## Wraps estimateSigmaFromLinearModel so every caller needing a starting sigma
## estimate raises the same failure message instead of a bare lm() error.
estimateStartingSigma <- function(data) {
  tryResult <- tryCatch(
    estimateSigmaFromLinearModel(data),
    error = function(e) e
  )
  if (inherits(tryResult, "error")) {
    nonFinite <- nonFinitePredictorNames(data@x)
    if (length(nonFinite) > 0L) {
      stop(
        "unable to obtain a starting estimate of sigma: predictor ",
        paste0("'", nonFinite, "'", collapse = ", "),
        " has infinite values; remove them or provide 'sigest'"
      )
    }
    stop("unable to obtain a starting estimate of sigma; provide one instead")
  }
  tryResult
}

## The names (or 1-based positions) of a dense predictor source's columns
## holding an infinite value, which the linear fit behind the starting sigma
## estimate cannot take.
nonFinitePredictorNames <- function(x) {
  if (predictorSourceIsSparse(x)) {
    return(character())
  }
  x <- as.matrix(x)
  bad <- which(colSums(is.infinite(x)) > 0L)
  if (length(bad) == 0L) {
    return(character())
  }
  if (is.null(colnames(x))) as.character(bad) else colnames(x)[bad]
}

## The binary/ordinal/nbinom weight policy, shared by every entry point that
## can reach these families: a probit has no tractable weighted latent-
## variable form and is refused, except that weights identically 1 are the
## unweighted likelihood and are treated as absent (SuperLearner-style
## callers pass obsWeights = rep(1, n) unconditionally), and weights that are
## all 0 or 1 name the rows in the data set and resolve to the active-row mask
## (isMaskedWeightFamily above); a logistic model treats weights as
## observation counts (its Polya-Gamma latent is a sum of per-copy draws), so
## they must be positive integers; ordinal and nbinom follow probit, mask
## included (an nbinom exposure belongs in the offset, not in a weight).
## Gaussian weights are unrestricted and reach here as a no-op.
##
## Returns the data with the policy applied and the mask the weights resolved
## to, NULL where they resolved to none: the weights slot is cleared in both
## the all-ones and the mask case, since neither family carries a weight
## channel, so the caller must install the mask on the sampler it builds.
## A logistic fit's weights are observation counts, and the weights of a
## binary fit's posterior predictive draw at new rows are its trial counts;
## either way, positive integers. 'what' names whose weights they are.
refuseNonCountWeights <- function(
  w,
  remedy = "",
  what = "logistic weights are observation counts"
) {
  if (anyNA(w) || any(w <= 0) || !all(is.finite(w)) || any(w != round(w))) {
    stop(what, " and must be positive integers", remedy, call. = FALSE)
  }
  invisible(NULL)
}

enforceWeightPolicy <- function(data, family) {
  if (is.null(data@weights)) {
    return(list(data = data, active = NULL))
  }
  # a data object's slot can be edited past its validity check, and a gaussian
  # fit would otherwise fail on the starting sigma without naming the weights;
  # NA is left to each family's own rule
  observed <- data@weights[!is.na(data@weights)]
  if (any(observed < 0)) {
    stop("'weights' must all be non-negative")
  }
  if (!all(is.finite(observed))) {
    stop("'weights' must all be finite")
  }
  active <- NULL
  if (isMaskedWeightFamily(family)) {
    w <- data@weights
    if (!anyNA(w) && all(w == 1)) {
      data@weights <- NULL
    } else if (!anyNA(w) && all(w == 0 | w == 1)) {
      active <- w
      data@weights <- NULL
    } else if (family == "probit") {
      stop(
        "probit models do not support weights other than 0 and 1, which mark ",
        "rows in and out of the likelihood as the sampler's $setActiveRows ",
        "does; fit integer count weights with family = \"logistic\", or ",
        "model continuous weights' latents directly"
      )
    } else if (family == "nbinom") {
      stop(
        "nbinom (count) models do not support weights other than 0 and 1, ",
        "which mark rows in and out of the likelihood as the sampler's ",
        "$setActiveRows does: exposure belongs in the offset as a ",
        "log-exposure term"
      )
    } else {
      stop(
        "ordinal models do not support weights other than 0 and 1, which ",
        "mark rows in and out of the likelihood as the sampler's ",
        "$setActiveRows does: a weighted truncated-normal latent likelihood ",
        "is not a coherent model"
      )
    }
  } else if (family == "logistic") {
    refuseNonCountWeights(
      data@weights,
      "; drop zero-count rows, and use a gaussian model for continuous weights"
    )
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
  list(data = data, active = active)
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
  shape,
  residDf,
  proposal.probs,
  monotone,
  interactions,
  blocks,
  variance,
  survivalStatus,
  hazardPeriods,
  bases,
  forests,
  evalEnv,
  residPrior = NULL,
  familySpec = NULL,
  basisRecords = NULL,
  written = names(matchedCall)
) {
  # a control taken from a fit whose first forest has a basis holds that
  # forest's tree count, where the bridge reads it. The count this call is to
  # inherit rides the forests' record and goes back on the slot before the
  # record is cleared with the rest: a multiplied forest's count is never the
  # next fit's
  carriedTreeCount <- attr(
    control,
    "bartcore.forests",
    exact = TRUE
  )$control.n.trees
  if (!is.null(carriedTreeCount)) {
    control@n.trees <- carriedTreeCount
  }
  # a caller-supplied control may have been taken from another fit, and the
  # bartcore.* attributes are that fit's model configuration; this call
  # attaches its own, so none of the incoming ones is honored
  for (attrName in grep(
    "^bartcore\\.",
    names(attributes(control)),
    value = TRUE
  )) {
    attr(control, attrName) <- NULL
  }
  # a forests = declaration puts a forest column on getTrees' table, whatever
  # its count, and nothing else of the control tells a one-forest declaration
  # from none
  attr(control, "bartcore.forestsDeclared") <- if (length(forests) > 0L) {
    TRUE
  }
  # what the caller stated of the arguments that are the forest with no
  # basis's own: bart()'s record of its own caller where it built this
  # control, and otherwise the control and the names 'written' in this call
  stated <- attr(control, plainStatedAttr, exact = TRUE)
  if (is.null(stated)) {
    stated <- plainForestStated(control, written)
  }
  attr(control, plainStatedAttr) <- NULL

  # a factor/logical/character response declares a classification model. The
  # single-forest engine here fits only the 2-level (probit) case; 3+ levels
  # are multinomial, which only bart(family = "multinomial") implements. A
  # numeric response takes the 0/1-vs-continuous path.
  # dbarts() is also reached anonymously through bartBT(), which has no family
  # formal at all; see resolveClassificationFamily's doc comment for why
  # its auto-branch message lists every single-forest entry point instead of
  # naming itself. probit/logistic on a 2-level categorical response proceed
  # as binary.
  autoDescription <- if (identical(requestedFamily, "auto")) {
    describeAutoResponse(data, survivalStatus)
  }
  family <- resolveClassificationFamily(
    data,
    family,
    "dbarts()/bartBT()/xbart",
    c("gaussian", "aft", "nbinom"),
    splitMultinomialMessage = TRUE,
    allowOrdinal = TRUE
  )
  # multinomial (K-forest softmax): the response is
  # the n x K count matrix on the data object, so the family is DECLARED by the
  # slot and not inferred from any response shape - data@y is that matrix's
  # trials vector. A counts-carrying object resolves to it from "auto" (and
  # is announced like any other resolution, below), and every other
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
    # the shape r: NA (the default) estimates it on the capped integer grid;
    # a supplied value FIXES it and must be a positive integer (v1 ships the
    # exact integer envelope, section 2). The C bridge reads the resolved spec
    # off the control attribute the sampler build attaches below: a positive
    # value fixes r, a non-positive value estimates it on the grid.
    shapeSpec <- resolveShape(shape)
    attr(control, "bartcore.shape") <- shapeSpec
  } else if (data@response.type == "numeric") {
    uniqueResponses <- unique(data@y)
    responseIsBinary <- length(uniqueResponses) == 2 &&
      all(sort(uniqueResponses) == c(0, 1))
    if (family == "auto") {
      family <- if (responseIsBinary) "probit" else "gaussian"
    } else if (family != "gaussian" && family != "aft" && !responseIsBinary) {
      # gaussian on a 0/1 response is a legitimate request; the binary
      # families need latent-variable coding. aft fits continuous log-times.
      refuseNonBinaryResponse(
        uniqueResponses,
        family,
        if (!is.null(hazardPeriods)) requestedFamily
      )
    }
  }
  # a factor, logical or character response of one class codes to a single
  # 0/1 value, which the numeric check above never sees
  if (isBinaryFamily(family) && responseHasSingleClass(data@y)) {
    refuseNonBinaryResponse(
      unique(data@y[!is.na(data@y)]),
      family,
      if (!is.null(hazardPeriods)) requestedFamily
    )
  }
  # aft draws sigma and rescales like gaussian; only the binary families are
  # latent-variable models on a fixed unit scale
  control@binary <- isBinaryFamily(family)
  # ordinal (cumulative probit) shares probit's fixed unit latent scale - sigma
  # fixed at 1, resid.prior fixed(1), no sigma estimate, leaf.scale 3.0 - but is
  # NOT binary: the bridge selects it by the bartcore.n.categories attribute
  # (not control@binary), and it reports K category levels. nbinom (counts) is
  # likewise a fixed-unit-scale family (sigma fixed at 1, the counts entering
  # kappa directly), selected by the bartcore.shape attribute.
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
  # bridge keeps the same checks as a backstop) - see enforceWeightPolicy's
  # own doc comment for the rule each family follows. Gaussian weights,
  # including a gaussian fit of a 0/1 response, are unrestricted and pass
  # through untouched. A latent family's 0/1 weights resolve to an active-row
  # mask the caller installs on the sampler it builds from this spec; the
  # weights slot is cleared either way, since the bridge takes none there.
  weightPolicy <- enforceWeightPolicy(data, family)
  data <- weightPolicy$data
  active <- weightPolicy$active

  if (is.na(data@sigma) && !fixedUnitScale) {
    data@sigma <- estimateStartingSigma(data)
  }

  # bart passes offset == something through no matter what; a latent-scale
  # fixed-unit family (binary, ordinal) keeps its zero offset as the
  # meaningful reference, as probit always has
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
  # forest(basis = factor(z))) names the ensembles the mean is a weighted sum
  # of, and every knob is per forest. A forest's defaults go by its KIND, read
  # here from the bases the model ends with. The forest with no basis, the
  # plain one, takes the fitting function's tree count, tree prior,
  # 'interactions' and 'blocks' wherever it stands, and a forest with a basis
  # the multiplied forest's defaults and its own constraints. The bridge reads
  # the FIRST forest's tree count from control@n.trees and its structure prior
  # from the model, so that forest's are written there, here, before anything
  # reads them. The rest ride the forests control attribute below.
  # PER FOREST, not one flag for the model: forest f is excused from declaring
  # a basis only when one reaches it some other way, the dbartsData(bases = )
  # route being the supported one. On the fitting path the declarations have
  # already been expanded onto the data object, so this reads them back.
  declaredBases <- if (is.null(bases)) data@bases else bases
  declaredHasBasis <- if (is.null(declaredBases)) {
    logical(0L)
  } else {
    !vapply(declaredBases, is.null, logical(1L))
  }
  forestSpec <- resolveForests(forests, interactions, blocks, declaredHasBasis)
  firstForest <- if (is.null(forestSpec)) NULL else forestSpec[[1L]]
  # a model with no bases has one forest, which is the plain one
  plain <- if (is.null(declaredBases)) 1L else plainForest(declaredHasBasis)
  # ahead of every other reading of these arguments, so that each is told
  # where it belongs and not what a model of several forests cannot take
  if (plain == 0L && length(declaredBases) > 1L) {
    refuseStatedWithNoPlainForest(stated, interactions, blocks)
  }
  plainSpec <- if (plain > 0L && length(forestSpec) >= plain) {
    forestSpec[[plain]]
  }
  # bart()'s own 'n.trees' and the plain forest's are one count; the control's
  # beside the forest's own leaves the forest's to govern
  if (
    identical(unname(stated["n.trees"]), "n.trees") &&
      !is.null(plainSpec$n.trees)
  ) {
    stop(
      "'n.trees' is given to the fitting function and to the forest with no ",
      "basis, which are the same count; give one",
      call. = FALSE
    )
  }
  # the fitting function's own, read once before any forest's statement is
  # written over the slot
  fitTreeCount <- control@n.trees
  fitInteractions <- interactions
  fitBlocks <- blocks
  firstIsPlain <- plain == 1L
  if (!is.null(firstForest$n.trees)) {
    control@n.trees <- firstForest$n.trees
  } else if (!firstIsPlain) {
    control@n.trees <- multipliedForestDefaults$n.trees
  }
  # the model carries the first forest's constraints: its own, and the fitting
  # function's where it is the plain forest
  if (!firstIsPlain) {
    interactions <- NULL
    blocks <- NULL
  }
  if (!is.null(firstForest$interactions)) {
    interactions <- firstForest$interactions
  }
  if (!is.null(firstForest$blocks)) {
    blocks <- firstForest$blocks
  }

  # the monotone spec arrives resolved by name at the door; its direction
  # vector and its prior ride the model below
  monotoneResolved <- resolveMonotone(monotone, data)
  monotoneDirections <- monotoneResolved$directions

  parsePriorsCall <- redirectCall(
    matchedCall,
    quoteInNamespace(parsePriors),
    callFormals = callFormals
  )
  parsePriorsCall <- setDefaultsFromFormals(
    parsePriorsCall,
    callFormals,
    "tree.prior",
    "leaf.prior"
  )
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$monotone <- monotoneDirections
  # a multi-forest fit is the one whose forests carry amplitude bases; the
  # calibration map pins every forest's k, so the binary chi-k default is
  # redirected to the fixed 2 rather than refused below
  parsePriorsCall$multiForest <- !is.null(declaredBases)
  parsePriorsCall$kHyperprior <- drawsLeafKByDefault(family)
  parsePriorsCall$parentEnv <- evalEnv

  # The residual prior has one home, the family object it rides, so it
  # arrives here already resolved (the entry point reads it off that object,
  # or off the retired flat spelling it still accepts) and NULL is the
  # package default.
  parsePriorsCall <- setCallArgument(
    parsePriorsCall,
    "resid.prior",
    if (is.null(residPrior)) quote(chisq) else residPrior
  )
  if (fixedUnitScale) {
    parsePriorsCall <- setCallArgument(
      parsePriorsCall,
      "resid.prior",
      quote(fixed(1))
    )
  }
  priors <- eval(parsePriorsCall)

  # the plain forest's count and tree prior: the fitting function's, read
  # before the first forest's are written over the model's
  fitTree <- list(
    n.trees = fitTreeCount,
    base = priors$tree.prior@base,
    power = priors$tree.prior@power
  )
  # the model's tree prior is the first forest's: a knob declared on that
  # forest restates its half, and one left out takes the default of the
  # forest's kind
  if (!is.null(firstForest$base)) {
    priors$tree.prior@base <- firstForest$base
  } else if (!firstIsPlain) {
    priors$tree.prior@base <- multipliedForestDefaults$base
  }
  if (!is.null(firstForest$power)) {
    priors$tree.prior@power <- firstForest$power
  } else if (!firstIsPlain) {
    priors$tree.prior@power <- multipliedForestDefaults$power
  }

  # The tree-move mixture rides the control. A caller that named it flat -
  # dbartsSpec's own argument, or the retired spelling on an entry point that
  # shed it - wins over the control's slot; NULL leaves the slot standing.
  if (!is.null(proposal.probs)) {
    control@proposal.probs <- resolveProposalProbs(proposal.probs)
  }

  # A monotone constraint restricts the forest to birth/death proposals: a
  # defaulted proposal.probs is forced to birth/death-only, an explicit
  # non-default one conflicts and errors.
  if (!is.null(monotoneDirections)) {
    control@proposal.probs <- monotoneProposalProbs(control@proposal.probs)
  }
  validObject(control)

  model <- newValidated(
    "dbartsModel",
    priors$tree.prior,
    priors$leaf.prior,
    priors$leaf.hyperprior,
    priors$resid.prior,
    family = family,
    # a named leaf-prior sd, translated to its anchor, overrides the family
    # default below in the engine, which converts it out of response units
    # against the transform; NA leaves that default in force
    prior.scale = priors$prior.scale,
    leaf.scale = defaultLeafScale(family)
  )

  # Student-t residuals: only a continuous
  # gaussian response carries them (the binary families and aft have their own
  # latent scale), refused here R-side to match the C bridge's backstop. The
  # resolved degrees of freedom ride the model's resid.df attribute the bridge
  # reads - the bartcore.survival precedent above - absent for the Gaussian law.
  # NA in the settings means estimate, which the bridge spells 0.
  if (!is.null(residDf) && family != "gaussian") {
    stop(
      "student residuals require a continuous gaussian response; family \"",
      requestedFamily,
      "\" has its own fixed error scale"
    )
  }
  if (!is.null(residDf)) {
    attr(model, "resid.df") <- if (is.na(residDf)) 0.0 else as.double(residDf)
  }
  # the family as the caller specified it, which a packaged fit carries for
  # family(); model@family is the engine's token
  if (!is.null(familySpec)) {
    attr(model, "family.spec") <- specifiedFamily(familySpec, family)
  }

  # the resolved per-column monotone directions and the prior ride two model
  # attributes the C bridge reads into SamplerOptions (the resid.df
  # precedent); a copy or reload rebuilds from the model, so both persist
  if (!is.null(monotoneDirections)) {
    attr(model, "monotone") <- monotoneDirections
    attr(model, "monotone.prior") <- monotoneResolved$prior
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
  # path is byte-for-byte unchanged. The partition covers the columns the
  # first forest may split on.
  #
  # Those columns are the first forest's 'vars', resolved once, here. On a
  # single forest a hazard fit's own period column, which rides last and which
  # the caller did not supply, stays allowed whatever 'vars' names.
  firstColumns <- resolveForestVars(firstForest$vars, data)
  singleForest <- is.null(declaredBases)
  if (singleForest && !is.null(firstColumns) && !is.null(hazardPeriods)) {
    firstColumns <- union(firstColumns, ncol(data@x))
  }
  blockSpec <- resolveBlocks(
    blocks,
    data,
    control@n.trees,
    availableColumns = firstColumns
  )
  if (!is.null(blockSpec)) {
    attr(model, "block.of.column") <- blockSpec$block.of.column
    attr(model, "block.tree.counts") <- blockSpec$block.tree.counts
  }

  # The K category forests are built from the softmax calibration map and the
  # CONSTANT-leaf instantiation only: a monotone
  # constraint or a non-constant leaf selects an instantiation the multinomial
  # factory does not build, a DART prior and a drawn k are unadjudicated
  # against the map's fixed anchor, the map owns every leaf scale so a named
  # leaf-prior sd has nowhere to land. Every one of these would
  # otherwise be dropped in silence, changing the fitted model without a word;
  # name each one instead. The bridge keeps its own backstops for the callers
  # that reach it without this layer. interactions() and blocks() are not here:
  # the bridge installs both on every category forest.
  if (identical(family, "multinomial")) {
    unsupportedMultinomial <- c(
      "a DART tree prior" = is(priors$tree.prior, "dbartsDartPrior"),
      "'split.probs'" = length(priors$tree.prior@splitProbabilities) > 0L,
      "'monotone'" = !is.null(monotoneDirections),
      "a linear leaf prior" = is(priors$leaf.prior, "dbartsLinearPrior"),
      "a Gaussian-process leaf prior" = is(priors$leaf.prior, "dbartsGPPrior"),
      "a 'k' hyperprior" = is(priors$leaf.prior@k, "dbartsLeafHyperprior"),
      "a named leaf-prior 'sd'" = !is.null(priors$leaf.prior@prior.sd),
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

  # a single forest's column restriction is a model fact and rides a model
  # attribute the C bridge reads, beside the two constraints above, so a copy
  # or a reload rebuilds it. Naming every column restricts nothing and stores
  # nothing; a multi-forest fit carries every forest's columns on the forests
  # control attribute instead. After the multinomial refusals above, so a
  # prior a category forest cannot take is named before its entries are read.
  if (
    singleForest &&
      !is.null(firstColumns) &&
      length(firstColumns) < ncol(data@x)
  ) {
    refuseNoSplittableColumn(
      priors$tree.prior@splitProbabilities,
      firstColumns
    )
    attr(model, "forest.columns") <- firstColumns
  }

  # the AFT survival family reads its per-observation status off this control
  # attribute; the C bridge validates it
  if (!is.null(survivalStatus)) {
    if (length(survivalStatus) != length(data@y)) {
      stop("survival status must have length ", length(data@y))
    }
    attr(control, "bartcore.survival") <- survivalStatus
  }
  # the discrete-time hazard marker: the
  # period grid, parked here for packageBartResults to read into $periods. The
  # C bridge never reads this attribute (unlike bartcore.survival), so a hazard
  # fit's draw stream is byte-identical to the by-hand binary fit's.
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
  # one selector vocabulary - and its n.trees/base/power knobs land on the
  # same control attribute the C factory reads.
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
  # the numeric branch's fractional-index refusal names the spelling the
  # caller actually wrote: 'vars' when the selector arrived through a
  # varianceForest() object, 'variance' for the plain shorthand
  varianceArgument <- if (is.null(varianceSpec)) "variance" else "vars"
  varianceColumns <- resolveVarianceColumns(
    varianceSelector,
    data,
    varianceArgument
  )
  if (!is.null(varianceColumns)) {
    if (!family %in% c("gaussian", "aft")) {
      # a hazard fit is a binary fit underneath; the caller named the hazard
      stop(
        "a variance forest requires family = \"gaussian\" or \"aft\"; ",
        "family \"",
        if (!is.null(hazardPeriods)) requestedFamily else family,
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
    if (
      is(priors$leaf.prior, "dbartsLinearPrior") ||
        is(priors$leaf.prior, "dbartsGPPrior")
    ) {
      stop(
        "a variance forest is not supported with a ",
        if (is(priors$leaf.prior, "dbartsLinearPrior")) {
          "linear"
        } else {
          "Gaussian-process"
        },
        " leaf prior; it takes constant leaves only"
      )
    }
    allColumns <- setequal(varianceColumns, seq_len(ncol(data@x)))
    varianceNTrees <- if (is.null(varianceSpec)) NULL else varianceSpec$n.trees
    varianceBase <- if (is.null(varianceSpec)) NULL else varianceSpec$base
    variancePower <- if (is.null(varianceSpec)) NULL else varianceSpec$power
    n.trees <- if (is.null(varianceNTrees)) 40L else varianceNTrees
    attr(control, "bartcore.variance") <- list(
      n.trees = coerceOrError(n.trees, "integer"),
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
    # K = 1 is not a shipped configuration (dec-A109). Both creation routes
    # reach here - the dbartsData(bases = ) one and the forests = one, whose
    # declarations forestBasisDeclarations carries down at any length - so this
    # is the single site the refusal is owed at. A varying-coefficient model
    # is one forest per coefficient function, intercept included.
    if (numForests < 2L) {
      fromData <- is.null(bases) &&
        !any(lengths(forestBasisDeclarations(forests)) > 0L)
      stop(
        "a multi-forest model needs at least two forests, and ",
        if (fromData) {
          paste0(
            "the data object carries ",
            numForests,
            "; for varying coefficients declare an intercept forest plus ",
            "one basis forest per covariate - forests = list(forest(), ",
            "forest(basis = z1), ...) on dbarts() or dbartsSpec(), or a ",
            "data object with bases = list(NULL, z1, ...) - or use a single ",
            "forest with linear() leaves; otherwise drop the basis"
          )
        } else {
          paste0(
            "this call's 'basis' declarations resolve to ",
            numForests,
            ": a forest with a 'basis' stands beside another forest. Write ",
            "the forest with no multiplier too, as y ~ forest(x1 + x2) + ",
            "forest(x1 + x2, basis = z1) or forests = list(forest(), ",
            "forest(basis = z1)), or use a single forest with linear() ",
            "leaves; otherwise drop the basis"
          )
        }
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
            "its threshold block is not shown to interleave with the ",
            "amplitude block"
          ),
          nbinom = paste0(
            "its shape block is not shown to interleave with the ",
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
    # non-constant leaf, a variance forest, an fp32
    # residual, a per-column cut cap, nor the Student-t error law. Every one of
    # those would otherwise be dropped in silence, changing the fitted model
    # without a word; name each one instead. 'proposal.probs' is NOT among them:
    # the chain carries the control's mixture onto every forest it builds, so a
    # coupling honors it rather than dropping it.
    unsupported <- c(
      "a DART tree prior" = is(priors$tree.prior, "dbartsDartPrior"),
      "'split.probs'" = length(priors$tree.prior@splitProbabilities) > 0L,
      "'monotone'" = !is.null(monotoneDirections),
      "a linear leaf prior" = is(priors$leaf.prior, "dbartsLinearPrior"),
      "a Gaussian-process leaf prior" = is(priors$leaf.prior, "dbartsGPPrior"),
      "a 'k' hyperprior" = is(priors$leaf.hyperprior, "dbartsChiHyperprior"),
      "a non-default 'k'" = is(
        priors$leaf.hyperprior,
        "dbartsFixedHyperprior"
      ) &&
        priors$leaf.hyperprior@k != 2.0,
      # "differs from the family default": defaultLeafScale(family), not a
      # gaussian-only literal
      "a non-default 'leaf.scale'" = model@leaf.scale !=
        defaultLeafScale(family),
      # the calibration map fixes every forest's leaf scale from the family's
      # own latent scale, so a named anchor has nowhere to land and the
      # leaf.scale gate above does not fire on it; parsePriors refuses a named
      # sd first, so this backstops a model built by hand
      "a named leaf-prior 'sd'" = !is.na(model@prior.scale),
      "Student-t residuals" = !is.null(residDf),
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
    # a coefficient is held only where the engine holds it at the value the
    # help states, which goes by the width of the basis and the position
    for (index in seq_len(numForests)) {
      if (identical(specs[[index]]$amplitude, "fixed")) {
        refuseHeldShape(
          index,
          if (hasBasis[index]) NCOL(data@bases[[index]]) else 0L
        )
      }
    }
    params <- forestParams(specs, hasBasis, family, fitTree)
    treeCounts <- vapply(params, function(forest) as.integer(forest[1L]), 0L)
    forestColumns <- lapply(
      seq_len(numForests),
      function(index) {
        if (index == 1L) {
          firstColumns
        } else {
          resolveForestVars(specs[[index]]$vars, data)
        }
      }
    )
    # every forest's label, fixed here for the life of the sampler: a list's
    # own name, else the text of a basis written as code, else its position
    labels <- forestLabels(
      names(forests),
      lapply(basisRecords, function(record) record$label),
      numForests
    )
    attr(control, "bartcore.forests") <- list(
      # one length-8 numeric per forest; the family selects the basis-free
      # channel's default median and the count the K-aware leaf scale factor
      params = params,
      # resolved 1-based column indices per forest, or NULL for unrestricted
      vars = forestColumns,
      # one label for each forest; see forestLabels()
      labels = labels,
      # the first forest's constraints are the model's, already resolved
      # above; the rest take their own, and the plain forest the fitting
      # function's where it states none, each resolved against the columns
      # and the tree count of the forest it lands on
      interactions = c(
        list(interactionSpec),
        lapply(
          seq_len(numForests)[-1L],
          function(index) {
            own <- specs[[index]]$interactions
            resolveInteractions(
              if (is.null(own) && index == plain) fitInteractions else own,
              data
            )
          }
        )
      ),
      blocks = c(
        list(blockSpec),
        lapply(
          seq_len(numForests)[-1L],
          function(index) {
            own <- specs[[index]]$blocks
            resolveBlocks(
              if (is.null(own) && index == plain) fitBlocks else own,
              data,
              treeCounts[index],
              availableColumns = forestColumns[[index]]
            )
          }
        )
      )
    )
    # what builds a basis written as code again at new rows, positional
    # against the forests as data@bases is and NULL where a forest's basis
    # is a value, which has no code to build from. Inert to the run: nothing
    # the engine reads
    if (!is.null(basisRecords)) {
      forestInfo <- attr(control, "bartcore.forests", exact = TRUE)
      forestInfo$basisTerms <- basisRecords
      attr(control, "bartcore.forests") <- forestInfo
    }
    # where the first forest has a basis the control's slot holds a multiplied
    # forest's count, so the count a later fit given this control inherits is
    # kept beside it: the plain forest's as it runs, or the fitting function's
    # own where no forest is plain. Inert to the run
    if (!firstIsPlain) {
      forestInfo <- attr(control, "bartcore.forests", exact = TRUE)
      forestInfo$control.n.trees <- if (plain > 0L) {
        treeCounts[plain]
      } else {
        fitTreeCount
      }
      attr(control, "bartcore.forests") <- forestInfo
    }
  }

  # every resolution of "auto" is announced once, here, after the refusals
  # above: bart() reaches this through dbarts(), and dbartsSpec() reads its
  # control's verbose
  if (!is.null(autoDescription)) {
    announceAutoFamily(control@verbose, family, autoDescription)
  }

  namedList(control, model, data, family, active)
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
  leaf.prior = normal,
  proposal.probs = c(
    birth_death = 0.6,
    swap = 0,
    change = 0.4,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  ),
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  forests = NULL,
  sigest = NULL,
  seed = NULL,
  family = c(
    "auto",
    "gaussian",
    "student",
    "probit",
    "logistic",
    "aft",
    "multinomial",
    "ordinal",
    "nbinom"
  ),
  survival = NULL,
  parentEnv = parent.frame(),
  ...
) {
  matchedCall <- match.call()

  # '...' exists so a retired spelling is refused by name, with its successor
  supplied <- dotNames(...)
  refuseForeignFrontDoorArgs(
    supplied,
    "dbartsSpec",
    names(formals(dbarts::dbartsSpec))
  )
  sigest <- resolveSigestArg(sigest, "dbartsSpec", "refuse")

  if (!inherits(data, "dbartsData")) {
    stop("'data' must be a dbartsData object; see ?dbartsData")
  }
  if (!inherits(control, "dbartsControl")) {
    stop("'control' must be a dbartsControl object; see ?dbartsControl")
  }
  familySpec <- resolveFamily(
    matchedCall$family,
    eval(formals(dbarts::dbartsSpec)$family),
    "dbartsSpec",
    parentEnv
  )
  family <- familySpec@token
  shape <- familySetting(familySpec, "shape", NA_real_)
  residPrior <- familySetting(familySpec, "sigma", NULL)
  refuseSigestUnderFixedPrior(residPrior, sigest, "sigest")
  # Student-t is a gaussian response carrying a degrees-of-freedom attribute
  # on this side of the bridge; the remap happens once, here
  residDf <- NULL
  if (identical(family, "student")) {
    residDf <- familySetting(familySpec, "df", NA_real_)
    family <- "gaussian"
  }

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

  seed <- resolveSeedArg(seed, "dbartsSpec", refuse = TRUE)
  if (!is.na(seed)) {
    control@seed <- seed
  }

  # the control owns the cut-point count, as it does inside dbarts(), but a data
  # object already carrying resolved per-column counts keeps them - a consumer
  # that set them deliberately is not silently overridden
  if (length(data@n.cuts) != ncol(data@x) || anyNA(data@n.cuts)) {
    data@n.cuts <- recycleNumCuts(control@n.cuts, ncol(data@x))
  }
  # an explicit sigest overrides whatever the data carries; NULL leaves it
  # alone, so a consumer's own starting estimate survives (an unset one is
  # estimated during resolution, exactly as for dbarts())
  if (!is.na(sigest)) {
    data@sigma <- validateSigest(sigest, "dbartsSpec")
  }

  # as on dbarts(): the forest constructors resolve by bare name inside their
  # arguments, here in parentEnv
  forestArguments <- resolveForestArguments(matchedCall, parentEnv)
  forests <- forestArguments$forests
  interactions <- forestArguments$interactions
  blocks <- forestArguments$blocks
  monotone <- forestArguments$monotone
  variance <- forestArguments$variance

  # this surface does no data ingestion of its own, so a declared basis is
  # read here, as dbarts() reads one, and reaches data@bases through the same
  # validation dbartsData() applies on the fitting path. There is no data
  # frame for a basis's code to name a column of, so code is the value it had
  # where forest() was called; the caller's data object has already had its
  # own 'subset' applied, so a basis covers the rows it holds
  basis <- NULL
  basisRecords <- NULL
  if (!is.null(forestBasisDeclarations(forests))) {
    declared <- readDeclaredBases(forests, NULL, length(data@y))
    forests <- declared$forests
    # the data object holds the rows it holds: every row of a code basis
    built <- buildFitBases(
      declared$reads,
      lapply(declared$reads, function(read) {
        if (is.null(read)) {
          NULL
        } else if (is.null(read$frame)) {
          expandValueBasis(read$value)
        } else {
          basisRowNumbers(length(data@y))
        }
      }),
      length(data@y)
    )
    # as in dbarts(): a list in which no forest declares a basis names no
    # multi-forest model, and falls through to resolveForests' refusal
    if (any(!vapply(built$bases, is.null, logical(1L)))) {
      basis <- built$bases
      if (any(!vapply(built$records, is.null, logical(1L)))) {
        basisRecords <- built$records
      }
    }
  }

  resolveSamplerSpec(
    matchedCall,
    formals(dbartsSpec),
    control,
    data,
    family,
    requestedFamily = familySpec@token,
    shape = shape,
    residDf = residDf,
    # the flat argument wins where the caller named it, and leaves the
    # control's own slot standing where they did not
    proposal.probs = if ("proposal.probs" %in% names(matchedCall)) {
      proposal.probs
    } else {
      NULL
    },
    monotone = monotone,
    interactions = interactions,
    blocks = blocks,
    variance = variance,
    survivalStatus = survival,
    hazardPeriods = NULL,
    bases = basis,
    forests = forests,
    evalEnv = parentEnv,
    residPrior = residPrior,
    familySpec = familySpec,
    basisRecords = basisRecords
  )
}
