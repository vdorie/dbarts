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

## Refuses a response a binary family cannot fit: one not coded 0/1. A
## response of one class is fitted instead, with one warning (dec-B333); see
## warnSingleClass.
refuseNonBinaryResponse <- function(family, y = NULL) {
  # glm's other spelling of a binomial, a proportion with the trials as
  # weights, which a logistic weight cannot carry: it counts copies of a row
  observed <- y[!is.na(y)]
  if (
    identical(family, "logistic") &&
      length(observed) > 0L &&
      all(observed >= 0 & observed <= 1)
  ) {
    stop(
      "under family = binomial the response is 0/1 or cbind(successes, ",
      "failures): a proportion with the trials as weights is written ",
      "cbind(k, n - k), with k = round(p * n) from a proportion p and its ",
      "trials n",
      call. = FALSE
    )
  }
  stop(
    "family \"",
    family,
    "\" requires a response coded 0/1",
    if (family == "logistic") {
      paste0(
        " or cbind(successes, failures) (family = binomial is the logit ",
        "link, a logistic fit)"
      )
    },
    call. = FALSE
  )
}

## The probit binomial, refused wherever a probit family meets successes and
## failures.
probitBinomialMessage <- paste0(
  "family = binomial(\"probit\") does not take a cbind(successes, failures) ",
  "response: the probit binomial is not fitted; use family = binomial, the ",
  "logit link, or a 0/1 response"
)

## A logistic family on a count matrix that is not two columns. 'what' says
## where the matrix is: "the data carry" or "the response is".
binomialCountColumnsMessage <- function(K, what) {
  paste0(
    what,
    " an n x ",
    K,
    " count matrix: family \"logistic\" takes two columns, successes and ",
    "failures, and family \"multinomial\" fits ",
    K,
    " categories"
  )
}

## The columns a binomial response expression reads, checked before cbind
## turns a factor to its codes or a character column to text: either is
## refused in the words of a response that is not numeric, whatever the
## expression does with it. 'data' and 'env' are where the formula's
## variables are found; 'prefix' leads the message where a door names itself.
refuseNonNumericBinomialColumns <- function(
  expression,
  data,
  env,
  prefix = ""
) {
  for (name in all.vars(expression)) {
    value <- tryCatch(
      eval(as.name(name), data, env),
      error = function(e) NULL
    )
    if (is.factor(value) || is.character(value)) {
      stop(
        prefix,
        "a binomial response's successes and failures must be numeric",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

## A binomial response checked: an n x 2 matrix of non-negative whole numbers,
## a row with a missing cell left to the caller's na.action. 'prefix' leads
## the message where a door names itself. Returns the counts as integers.
validateBinomialCounts <- function(counts, prefix = "") {
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.numeric(counts) && !is.logical(counts)) {
    stop(
      prefix,
      "a binomial response's successes and failures must be numeric",
      call. = FALSE
    )
  }
  # a count past the integer range, an infinite one included, is past the
  # trials cap, and is refused here rather than read as NA by the integer
  # coercion below
  bad <- which(
    rowSums(is.na(counts)) == 0L &
      (rowSums(counts < 0) > 0L |
        rowSums(counts != round(counts)) > 0L |
        rowSums(counts > .Machine$integer.max) > 0L)
  )
  if (length(bad) > 0L) {
    stop(
      prefix,
      "a binomial response's successes and failures must be non-negative ",
      "whole numbers; row(s) ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ...",
      call. = FALSE
    )
  }
  storage <- counts
  storage.mode(storage) <- "integer"
  storage
}

## The one warning a binary fit gives when its response holds one class: the
## model is defined, and the prior alone pushes the probability toward 0 or 1.
## A hazard fit ('hazard', the token the caller gave) is a binary fit on
## person-period rows the caller never wrote, so its warning speaks of
## subjects and events; 'y' is the class the rows hold.
warnSingleClass <- function(family, hazard = NULL, y = NULL) {
  if (is.null(hazard)) {
    warning(
      "family \"",
      family,
      "\" is fitted to a response with a single class; the prior alone drives ",
      "the fitted probabilities toward 0 or 1",
      call. = FALSE
    )
  } else {
    warning(
      "family \"",
      hazard,
      "\" is fitted with ",
      if (y == 0) {
        "no event: every subject is censored; the prior alone drives the "
      } else {
        paste0(
          "no period at risk without an event: every subject has its event ",
          "in the first period; the prior alone drives the "
        )
      },
      "hazards toward ",
      if (y == 0) "0" else "1",
      ", so the survival curves sit where the prior puts them",
      call. = FALSE
    )
  }
}

## Whether a coded response holds one class of a binary one: every value 0,
## or every value 1, whatever encoding (numeric, logical, factor, character)
## it was coded from.
responseHasSingleClass <- function(y) {
  values <- unique(y[!is.na(y)])
  length(values) == 1L && values %in% c(0, 1)
}

## The one warning an ordinal or multinomial fit gives when its response shows
## one category: the model is defined, its priors proper, and the prior alone
## places the probabilities, as with a binary fit to one class.
warnSingleCategory <- function(family) {
  warning(
    "family \"",
    family,
    "\" is fitted to a response with a single observed category; the prior ",
    "alone drives the fitted probabilities",
    call. = FALSE
  )
}

## The families whose response is a count or a category, where dbartsData's
## precision warning, whose remedy is to rescale, does not apply. A formula
## hazard fit's response at that point is the log time standing in for the
## binary rows it expands to, which the warning does not describe either. A
## binary family is not listed: a near-constant response that is not one
## class keeps the warning beside the family's refusal.
discreteResponseFamilies <- c(
  "ordinal",
  "nbinom",
  "multinomial",
  hazardFamilyTokens
)

## Evaluates 'expr' with dbartsData's precision warning held back: the value,
## and the held warning or NULL.
holdingResponsePrecision <- function(expr) {
  held <- NULL
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      if (startsWith(conditionMessage(w), responsePrecisionWarningStem)) {
        held <<- w
        invokeRestart("muffleWarning")
      }
    }
  )
  list(value = value, held = held)
}

## dbartsData warns of a response whose values are indistinguishable at
## double precision before any family is known, and a constant response is
## one. A count or a category has nothing to rescale: a family that fits one,
## a binary family fit to a single class, or a response with a categorical
## encoding, which "auto" resolves to a binary, ordinal or multinomial fit,
## gets its own family's words where the response is degenerate (a binary
## fit's single class, an ordinal or multinomial fit's single category) or its
## refusal, so the warning is held back there and raised everywhere else.
withBinaryResponsePrecision <- function(family, expr) {
  result <- holdingResponsePrecision(expr)
  data <- result$value
  if (
    !is.null(result$held) &&
      family %not_in% discreteResponseFamilies &&
      !(is(data, "dbartsData") &&
        (data@response.type != "numeric" ||
          (family %in%
            c("probit", "logistic") &&
            responseHasSingleClass(data@y))))
  ) {
    warning(result$held)
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

## The families whose 0/1 weights are membership rather than precision - the
## row leaves the likelihood, keeps its leaf occupancy, its latent and its
## fitted value - which the engine composes with the active-row mask, and
## which a fit records as its rows in. Any other value is a weighted latent
## likelihood, which probit and ordinal have no coherent form for, and stays
## refused; on nbinom it is a count of identical observations.
isMaskedWeightFamily <- function(family) {
  family %in% c("probit", "ordinal", "nbinom")
}

## weights.type as the caller gave it: NA_character_ for NULL, the family's
## own reading, else "precision" or "frequency".
validateWeightsType <- function(weights.type) {
  if (is.null(weights.type)) {
    return(NA_character_)
  }
  if (
    !is.character(weights.type) ||
      length(weights.type) != 1L ||
      is.na(weights.type) ||
      !(weights.type %in% c("precision", "frequency"))
  ) {
    stop(
      "'weights.type' must be NULL, \"precision\" or \"frequency\"",
      call. = FALSE
    )
  }
  weights.type
}

## The reading 'stated' resolves to on the engine family 'family', refusing
## nothing: NA is the family's own.
familyWeightsType <- function(stated, family) {
  if (!is.na(stated)) {
    return(stated)
  }
  if (family %in% c("gaussian", "aft")) "precision" else "frequency"
}

## The reading a sampler of the engine family 'family' is created under, from
## the stated one (validateWeightsType's value): NA resolves to precision where
## there is a residual scale (gaussian, its student and heteroscedastic forms)
## and to frequency, counts or 0/1 membership, everywhere else. aft takes no
## weights and reads precision. 'shown' is the family as the caller wrote it,
## for the refusals; 'weights' the data's, for the 0/1 rule of the membership
## families under a stated "frequency".
resolveWeightsType <- function(stated, family, shown = family, weights = NULL) {
  if (identical(family, "aft")) {
    if (identical(stated, "frequency")) {
      stop(
        "weights.type = \"frequency\" is not taken under family = \"aft\", ",
        "which takes no weights",
        call. = FALSE
      )
    }
    return("precision")
  }
  if (identical(family, "gaussian")) {
    return(if (is.na(stated)) "precision" else stated)
  }
  if (identical(stated, "precision")) {
    stop(
      "weights.type = \"precision\" needs a residual scale, and family = \"",
      shown,
      "\" has none: its weights count identical observations or, at 0 and ",
      "1, mark rows in and out; leave weights.type NULL",
      call. = FALSE
    )
  }
  if (
    identical(stated, "frequency") &&
      family %in% c("probit", "ordinal") &&
      !is.null(weights) &&
      any(!is.na(weights) & weights != 0 & weights != 1)
  ) {
    stop(
      "weights.type = \"frequency\" is not taken under family = \"",
      shown,
      "\" for weights other than 0 and 1: a frequency weight there needs a ",
      "hidden value for each copy, which the sampler does not draw; weights of ",
      "0 and 1 mark rows out of and in the likelihood under either reading",
      call. = FALSE
    )
  }
  "frequency"
}

## The rows in the likelihood: weight positive and no mask switching them off.
## NA counts as out, as in sigest's regression.
rowsInLikelihood <- function(
  weights,
  active = NULL,
  n = max(length(weights), length(active))
) {
  rows <- rep_len(TRUE, n)
  if (!is.null(weights)) {
    rows <- rows & !is.na(weights) & weights > 0
  }
  if (!is.null(active)) {
    rows <- rows & active != 0
  }
  rows
}

## Whether a row is out of the likelihood and the rows in hold fewer than two
## distinct values of 'value'. Rule 4 of response-scale-rows: nothing is out,
## or two values are in, is not that case.
rowsInLackTwoValues <- function(value, rowsIn) {
  !all(rowsIn) && length(unique(value[rowsIn & !is.na(value)])) < 2L
}

## Whether a data object created with rows out has too little in the
## likelihood to read the response's scale or the prior's center from, for the
## families that take something from the response and can be created with a
## row out (gaussian with any residual distribution, and nbinom; no other: an
## aft sampler takes no weights, and a mask reaches it only after creation,
## where the engine answers). The scale and the center are then taken over
## every row, and so is the starting sigma: the caller hands the estimate the
## data without its weights. 'active' is a mask composed with the weights,
## NULL for none. Raises the creation warning unless 'warn' is FALSE.
scaleRowsFallBack <- function(data, family, active = NULL, warn = TRUE) {
  if (!(family %in% c("gaussian", "nbinom"))) {
    return(FALSE)
  }
  n <- length(data@y)
  rowsIn <- rowsInLikelihood(data@weights, active, n)
  value <- if (identical(family, "nbinom") || is.null(data@offset)) {
    data@y
  } else {
    data@y - data@offset
  }
  if (!rowsInLackTwoValues(value, rowsIn)) {
    return(FALSE)
  }
  if (warn) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      "fewer than two distinct values of the response are in the ",
      "likelihood (",
      countInLikelihood(sum(rowsIn), n, "observation"),
      "), so the response's scale and the prior's center are taken from ",
      takenFromAll(n)
    )
  }
  TRUE
}

## The starting sigma for a data object, by the rule 'sigest' named: the
## residual sd of the regression of the response on the predictors and on
## each forest's basis (sigmaBasisColumns), which the data object holds as
## the model ends with; under 'frequency' that of the regression on the rows
## written out as many times as their weights count. Wraps estimateSigmaFromLinearModel so every caller
## raises the same failure message instead of a bare lm() error. An infinite
## predictor value is refused first, for every design and before any routine
## runs; a basis is finite by its own check. Under 'verbose' one line says
## what is about to be computed, which cannot be interrupted once begun, and
## another where no regression was defined and the sd of the response stands
## in.
estimateStartingSigma <- function(
  data,
  rule = "auto",
  verbose = FALSE,
  frequency = FALSE
) {
  nonFinite <- nonFinitePredictorNames(data@x)
  if (length(nonFinite) > 0L) {
    stop(
      "unable to obtain a starting estimate of sigma: predictor ",
      paste0("'", nonFinite, "'", collapse = ", "),
      " has infinite values; remove them or provide 'sigest'"
    )
  }
  # the rows with a response, an offset and a positive weight
  offset <- if (is.null(data@offset)) 0 else data@offset
  weights <- if (is.null(data@weights)) 1 else data@weights
  size <- sum(!is.na(data@y + offset) & !is.na(weights) & weights > 0)
  fellBack <- FALSE
  bases <- sigmaBasisColumns(data@bases, dataBasisLevels(data))
  tryResult <- tryCatch(
    withCallingHandlers(
      {
        design <- startingSigmaDesign(data@x, rule, bases)
        sparse <- inherits(design, "dbartsSparseSigmaDesign")
        size[2L] <- if (sparse) ncol(design$X) else ncol(design)
        announceStartingSigma(
          verbose,
          "estimating the starting sigma by a ",
          if (sparse) "sparse" else "dense",
          " linear regression on ",
          size[1L],
          " rows and ",
          size[2L],
          " columns; see 'sigest'"
        )
        estimateSigmaFromLinearModel(data, rule, design, frequency)
      },
      dbartsSigmaFallback = function(condition) {
        fellBack <<- sigmaFallbackRecord(condition)
      }
    ),
    error = function(e) e
  )
  if (inherits(tryResult, "error")) {
    stop(
      "unable to obtain a starting estimate of sigma; provide one instead: ",
      conditionMessage(tryResult)
    )
  }
  if (isTRUE(fellBack)) {
    announceStartingSigma(
      verbose,
      "starting sigma is the ",
      if (isTRUE(attr(fellBack, "weighted"))) "weighted ",
      "sd of the response (",
      format(tryResult, digits = 4L),
      "): the linear model on ",
      size[2L],
      if (is.null(bases)) {
        " columns (factors as indicators)"
      } else {
        " columns (factors as indicators, each forest's basis among them)"
      },
      " leaves no residual degrees of freedom in ",
      size[1L],
      " rows; supply 'sigest' to set it"
    )
  }
  tryResult
}

## The data object with its weights slot emptied, for an estimate made over
## every row.
withoutWeights <- function(data) {
  data@weights <- NULL
  data
}

## One line about the starting sigma under 'verbose', in announceAutoFamily's
## form and classed, so a door that prints its own way can take the text.
announceStartingSigma <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    class <- c("dbartsStartingSigmaMessage", "dbartsMessage", "message")
    message(structure(
      class = c(class, "condition"),
      list(message = paste0(..., "\n"), call = NULL)
    ))
  }
}

## The names (or 1-based positions) of a predictor source's columns holding an
## infinite value, which the regression behind the starting sigma estimate
## cannot take. A sparse-stored column is read from its stored entries.
nonFinitePredictorNames <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    bad <- unique(sparseEntryColumns(x)[is.infinite(x@x)])
    names <- x@Dimnames[[2L]]
  } else {
    bad <- which(vapply(
      seq_len(ncol(x)),
      function(j) {
        column <- predictorSourceColumn(x, j, ncol(x), nrow(x))
        any(is.infinite(if (is.list(column)) column$x else column))
      },
      FALSE
    ))
    names <- colnames(x)
  }
  if (is.null(names)) as.character(bad) else names[bad]
}

## The weights of a count family's posterior predictive draw are the trial
## counts or the observations of the draw, non-negative integers, a 0 drawing
## a count of 0. 'what' names whose weights they are.
refuseNonCountWeights <- function(w, remedy = "", what) {
  if (anyNA(w) || any(w < 0) || !all(is.finite(w)) || any(w != round(w))) {
    stop(what, " and must be non-negative integers", remedy, call. = FALSE)
  }
  invisible(NULL)
}

## A logistic fit's weights: observation counts, a 0 a row out of the
## likelihood.
refuseLogisticWeights <- function(w) {
  if (anyNA(w) || any(w < 0) || !all(is.finite(w)) || any(w != round(w))) {
    stop(
      "logistic weights are observation counts and must be non-negative ",
      "whole numbers",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## Whole non-negative counts, the weights of the count families, refused in
## 'family''s words otherwise.
refuseNonCountFamilyWeights <- function(w, family) {
  if (anyNA(w) || any(w < 0) || !all(is.finite(w)) || any(w != round(w))) {
    stop(
      family,
      " weights are observation counts and must be non-negative whole numbers",
      if (identical(family, "nbinom")) {
        "; an exposure belongs in the offset as a log-exposure term"
      },
      call. = FALSE
    )
  }
  invisible(NULL)
}

## The binary/ordinal/nbinom weight policy, shared by every entry point that
## can reach these families: a probit has no tractable weighted latent-
## variable form and is refused, except that weights identically 1 are the
## unweighted likelihood and are treated as absent (SuperLearner-style
## callers pass obsWeights = rep(1, n) unconditionally), and weights that are
## all 0 or 1 name the rows in the data set (isMaskedWeightFamily above),
## which the engine composes with the active-row mask; ordinal follows probit.
## The count families, logistic, nbinom and multinomial, read weights as
## counts of identical observations (each copy's Polya-Gamma draw summed), so
## they must be whole numbers, a 0 a row out; an nbinom vector of 0s and 1s is
## kept as probit's is. Gaussian weights are unrestricted and reach here as a
## no-op.
##
## Returns the data with the policy applied: the weights slot cleared where
## the weights are all ones on a family whose weights are membership, and
## kept otherwise.
enforceWeightPolicy <- function(data, family) {
  if (!is.null(data@weights)) {
    refuseFamilyWeightValues(data@weights, family, "'weights'")
    w <- data@weights
    if (
      (isMaskedWeightFamily(family) || family == "multinomial") &&
        !anyNA(w) &&
        all(w == 1)
    ) {
      data@weights <- NULL
    }
  }
  # a test row is drawn as a training row of its weight is, so its weight is
  # held to the same rule, and a missing one is refused as a training row's is
  # (predict's na.action is where a row is left to NA draws)
  test <- data@weights.test
  refuseMissingTestWeights(test)
  if (!is.null(test) && any(!is.na(test))) {
    refuseFamilyWeightValues(test[!is.na(test)], family, "'weights.test'")
  }
  data
}

## A test frame's missing weight, refused naming the rows and the remedies.
refuseMissingTestWeights <- function(weights) {
  if (anyNA(weights)) {
    rows <- which(is.na(weights))
    stop(
      "the weights of the test rows ",
      toString(utils::head(rows, 5L)),
      if (length(rows) > 5L) paste0(" and ", length(rows) - 5L, " more"),
      " are missing; fill them, drop those rows from 'test', or fit without ",
      "them and give predict(fit, newdata, na.action = na.pass) the rows to ",
      "leave at NA draws",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## The values 'family' takes as weights, 'what' naming whose they are.
refuseFamilyWeightValues <- function(w, family, what) {
  # a data object's slot can be edited past its validity check, and a gaussian
  # fit would otherwise fail on the starting sigma without naming the weights;
  # NA is left to each family's own rule
  observed <- w[!is.na(w)]
  if (any(observed < 0)) {
    stop(what, " must all be non-negative")
  }
  if (!all(is.finite(observed))) {
    stop(what, " must all be finite")
  }
  if (isMaskedWeightFamily(family)) {
    if (!anyNA(w) && all(w == 0 | w == 1)) {
      # membership: kept on the data object, the engine composing it with
      # the active-row mask
    } else if (family == "nbinom") {
      # counts of identical observations, a 0 a row out
      refuseNonCountFamilyWeights(w, family)
    } else if (family == "probit") {
      stop(
        "probit models do not support weights other than 0 and 1, which mark ",
        "rows in and out of the likelihood as the sampler's $setActiveRows ",
        "does; fit integer count weights with family = \"logistic\", or ",
        "model continuous weights' latents directly"
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
    refuseLogisticWeights(w)
  } else if (family == "multinomial") {
    refuseNonCountFamilyWeights(w, family)
  }
  invisible(NULL)
}

## Whether a tree-prior record's split probabilities are other than equal:
## DART's, which it draws, or a vector that is not uniform.
varySplitProbabilities <- function(treePrior) {
  is(treePrior, "dbartsDartPrior") ||
    length(treePrior@splitProbabilities) > 0L
}

## survivalStatus and hazardPeriods carry the two survival markers the entry
## point resolved from the raw response; both are NULL for every other family.
## `treePrior` is a record a retired tree.prior handed over whole, which the
## tree prior is built from in place of 'split.probs'.
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
  forests,
  evalEnv,
  residPrior = NULL,
  familySpec = NULL,
  basisRecords = NULL,
  sigestRule = "auto",
  treePrior = NULL,
  weightsType = NA_character_
) {
  # a data object handed in whole is held to its own record here, at
  # dbarts() and at dbartsSpec() alike
  checkBasisRecord(data)
  # what the door set aside of a control taken from a fit (a carried value
  # states nothing), for the not-used warning below
  carried <- attr(control, controlCarriedAttr, exact = TRUE)
  attr(control, controlCarriedAttr) <- NULL
  # what the caller stated of the arguments that are the forest with no
  # basis's own, recorded at the door (readDoorCopies)
  stated <- attr(control, plainStatedAttr, exact = TRUE)
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
    allowOrdinal = TRUE,
    shownFamily = requestedFamily
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
  # a logistic family on two count columns is the binomial: successes and
  # failures, data@y their trials
  binomial <- identical(family, "logistic") && !is.null(counts)
  if (binomial && ncol(counts) != 2L) {
    stop(binomialCountColumnsMessage(ncol(counts), "the data carry"))
  }
  if (identical(family, "probit") && NCOL(counts) == 2L) {
    stop(probitBinomialMessage)
  }
  if (!identical(family, "multinomial") && !binomial && !is.null(counts)) {
    stop(
      "the data carry an n x K count matrix, which only family ",
      "\"multinomial\" fits; drop 'counts' to fit family \"",
      family,
      "\" against their trials"
    )
  }
  # the ordinal category count and the negative-binomial shape setting the
  # model records, NA on every other family
  numCategories <- NA_integer_
  shapeSpec <- NA_real_
  if (identical(family, "multinomial")) {
    # the counts ARE the response and the bridge reads K off their column
    # count, so nothing is recoded and the model records no category count;
    # data@y is already the trials. Two columns took a flat offset, which
    # might have been a binomial's
    if (!is.null(data@offset)) {
      refuseFlatOffsetOnMultinomial(data@offset)
    }
    if (length(data@offset.test) > 0L) {
      refuseFlatOffsetOnMultinomial(data@offset.test, "offset.test")
    }
    if (sum(colSums(counts, na.rm = TRUE) > 0) == 1L) {
      warnSingleCategory(family)
    }
  } else if (identical(family, "ordinal")) {
    # ordinal (cumulative probit): a single-forest
    # fixed-unit-scale model like probit, but K-level. Recode the response to
    # the 1-based category codes the engine reads, and record K on the model,
    # where the bridge reads it to select OrdinalResponse. The resolved ordered
    # levels ride the data object for the round-trip.
    ordinal <- resolveOrdinalResponse(data)
    data@y <- ordinal$y
    data@response.levels <- ordinal$levels
    numCategories <- as.integer(ordinal$K)
    if (length(unique(ordinal$y[!is.na(ordinal$y)])) == 1L) {
      warnSingleCategory(family)
    }
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
    # exact integer envelope, section 2). The model records the resolved spec,
    # where the bridge reads it: a positive value fixes r, a non-positive value
    # estimates it on the grid.
    shapeSpec <- resolveShape(shape)
  } else if (binomial) {
    # the counts are the response, checked where they were read; data@y is
    # their trials, which no 0/1 rule reads. Its offset is flat, as every
    # logistic fit's
    if (!is.null(data@offset.category)) {
      refuseMatrixOffset(data@offset.category, "offset")
    }
    if (
      all(counts[, 1L] == 0L, na.rm = TRUE) ||
        all(counts[, 2L] == 0L, na.rm = TRUE)
    ) {
      warnSingleClass(family)
    }
  } else if (data@response.type == "numeric") {
    uniqueResponses <- unique(data@y)
    responseIsBinary <- length(uniqueResponses) == 2 &&
      all(sort(uniqueResponses) == c(0, 1))
    if (family == "auto") {
      family <- if (responseIsBinary) "probit" else "gaussian"
    } else if (
      family != "gaussian" &&
        family != "aft" &&
        !responseIsBinary &&
        !responseHasSingleClass(data@y)
    ) {
      # gaussian on a 0/1 response is a legitimate request; the binary
      # families need latent-variable coding. aft fits continuous log-times.
      refuseNonBinaryResponse(family, data@y)
    }
  }
  # a factor, logical or character response of one class codes to a single
  # 0/1 value, which the numeric check above never sees
  if (isBinaryFamily(family) && !binomial && responseHasSingleClass(data@y)) {
    warnSingleClass(
      family,
      if (!is.null(hazardPeriods)) requestedFamily,
      unique(data@y[!is.na(data@y)])
    )
  }
  # aft draws sigma and rescales like gaussian; only the binary families are
  # latent-variable models on a fixed unit scale
  control@binary <- isBinaryFamily(family)
  # ordinal (cumulative probit) shares probit's fixed unit latent scale - sigma
  # fixed at 1, resid.prior fixed(1), no sigma estimate, leaf.scale 3.0 - but is
  # NOT binary: the bridge selects it by the model's category count, and it
  # reports K category levels. nbinom (counts) is likewise a fixed-unit-scale
  # family (sigma fixed at 1, the counts entering kappa directly), selected by
  # the model's shape setting.
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
  # through untouched. A latent family's 0/1 weights stay on the data object,
  # which the engine composes with any mask installed later.
  # the reading is resolved first, so a stated one is refused in its own
  # words before a family's rule on the values could refuse the weights
  weightsType <- resolveWeightsType(
    weightsType,
    family,
    if (identical(requestedFamily, "auto")) family else requestedFamily,
    data@weights
  )
  data <- enforceWeightPolicy(data, family)

  # the rows in the likelihood leave nothing to read a scale from: every
  # number is taken over every row, the starting sigma's regression included
  allRows <- scaleRowsFallBack(data, family)

  if (is.na(data@sigma) && !fixedUnitScale) {
    data@sigma <- if (is(residPrior, "dbartsFixedPrior")) {
      # a fixed residual standard deviation is the scale itself: no estimate
      # has anything to calibrate, so none is made and the slot holds it
      sqrt(residPrior@value)
    } else {
      estimateStartingSigma(
        if (allRows) withoutWeights(data) else data,
        sigestRule,
        control@verbose,
        frequency = identical(weightsType, "frequency")
      )
    }
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
  # the multiplied forest's defaults and its own constraints. The FIRST
  # forest's tree count is written on the control and its constraints on the
  # model here, before the priors read them; every forest's statement, the
  # first's among them, goes into the model's forest records below, which the
  # bridge reads.
  # PER FOREST, not one flag for the model: the declarations have already
  # been built onto the data object, so this reads them back.
  declaredBases <- data@bases
  declaredHasBasis <- if (is.null(declaredBases)) {
    logical(0L)
  } else {
    !vapply(declaredBases, is.null, logical(1L))
  }
  basisTexts <- lapply(basisRecords, function(record) record$label)
  forestSpec <- resolveForests(
    forests,
    interactions,
    blocks,
    declaredHasBasis,
    basisTexts
  )
  # how a refusal names a forest: by its label, as the fit's results do
  knobLabels <- if (!is.null(forestSpec)) {
    forestLabels(names(forests), basisTexts, length(forests))
  }
  firstForest <- if (is.null(forestSpec)) NULL else forestSpec[[1L]]
  # a model with no bases has one forest, which is the plain one
  plain <- if (is.null(declaredBases)) 1L else plainForest(declaredHasBasis)
  # ahead of every other reading of these arguments, so that each is told
  # where it belongs and not what a model of several forests cannot take
  if (length(plain) == 0L && length(declaredBases) > 1L) {
    refuseStatedWithNoPlainForest(
      stated,
      interactions,
      blocks,
      function(name) {
        plainTreePriorValue(
          matchedCall,
          name,
          evalEnv,
          multipliedForestDefaults[[name]]
        )
      }
    )
  }
  # a partition or a constraint of one forest's predictors cannot be each
  # forest's: beside two or more forests with no basis it is stated on one
  if (length(plain) > 1L) {
    refuseConstraintOnSeveralPlain(interactions, blocks, length(plain))
  }
  plainSpecs <- lapply(plain[plain <= length(forestSpec)], function(index) {
    forestSpec[[index]]
  })
  # the fitting function's own 'n.trees', 'power' and 'base' and the plain
  # forests' are one setting each, reaching the forests that state none, the
  # count in either spelling
  if (length(plainSpecs) == length(plain) && length(plain) > 0L) {
    refuseStatedOnPlainForest(
      stated,
      plainSpecs,
      forestKnobLabel(knobLabels, plain[1L])
    )
  }
  # the forests that take the fitting function's leaf prior as one forest's:
  # the share of a total, or the size of a drawn coefficient's trees
  shares <- !is.null(declaredBases) &&
    length(leafPriorForests(forestSpec, plain)) > 0L
  # the fitting function's own, read once before any forest's statement is
  # written over the slot
  fitTreeCount <- control@n.trees
  fitInteractions <- interactions
  fitBlocks <- blocks
  firstIsPlain <- 1L %in% plain
  firstSplit <- if (firstIsPlain && length(plain) > 1L) {
    allSpecs <- rep(list(NULL), length(declaredHasBasis))
    allSpecs[seq_along(forestSpec)] <- forestSpec
    plainTreeCounts(allSpecs, declaredHasBasis, fitTreeCount)[[1L]]
  }
  if (!is.null(firstForest$n.trees)) {
    control@n.trees <- firstForest$n.trees
  } else if (!firstIsPlain) {
    control@n.trees <- multipliedForestDefaults$n.trees
  } else if (!is.null(firstSplit)) {
    control@n.trees <- as.integer(firstSplit)
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
    "leaf.prior"
  )
  # a retired tree.prior's whole record, where it handed one over
  if (!is.null(treePrior)) {
    parsePriorsCall$treePrior <- treePrior
  }
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$monotone <- monotoneDirections
  # a multi-forest fit is the one whose forests carry amplitude bases; the
  # calibration map pins every forest's k, so the binary chi-k default is
  # redirected to the fixed 2 rather than refused below
  parsePriorsCall$multiForest <- !is.null(declaredBases) && !shares
  parsePriorsCall$kHyperprior <- drawsLeafKByDefault(family)
  parsePriorsCall$parentEnv <- evalEnv
  # the one forest of a single-forest model states its size as the fitting
  # function's leaf prior does: its 'sd' is that prior's, in the same units,
  # and the two together state one thing twice
  if (is.null(declaredBases) && !is.null(firstForest$sd)) {
    if ("leaf.prior" %in% names(stated)) {
      stop(
        "this model has one forest and states its size twice: 'sd' on the ",
        "forest and '",
        stated[["leaf.prior"]],
        "' given to the fitting function are one statement; give one",
        call. = FALSE
      )
    }
    parsePriorsCall <- setCallArgument(
      parsePriorsCall,
      "leaf.prior",
      bquote(normal(sd = .(firstForest$sd)))
    )
  }

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
  # before the first forest's are written over the model's. Where the first
  # forest is the plain one and states a knob, the fitting function stated
  # none, which was refused above
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
    for (index in seq_along(forestSpec)) {
      refuseMonotoneForestMixture(
        forestSpec[[index]]$proposal.probs,
        forestKnobLabel(knobLabels, index)
      )
    }
  }
  validObject(control)

  model <- newValidated(
    "dbartsModel",
    priors$tree.prior,
    priors$leaf.prior,
    priors$leaf.hyperprior,
    priors$resid.prior,
    family = family,
    # a named leaf-prior sd rides here as the model's encoding, twice the sd
    # beside a reference k of 2; the engine divides it back to the sd and
    # states it against the data's scale. NA leaves the default below in force
    prior.scale = priors$prior.scale,
    leaf.scale = defaultLeafScale(family)
  )
  # the reading as resolved, which the bridge reads and creation fixes
  model@weights.type <- weightsType

  # Student-t residuals: only a continuous
  # gaussian response carries them (the binary families and aft have their own
  # latent scale), refused here R-side to match the C bridge's backstop. The
  # resolved degrees of freedom ride the model's resid.df attribute the bridge
  # reads, absent for Gaussian errors.
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
  # Those columns are the first forest's 'vars', resolved once, here. A
  # hazard fit's own period column, which rides last and which the caller did
  # not supply, stays allowed in every forest whatever its 'vars' names.
  withPeriod <- function(columns) {
    if (is.null(columns) || is.null(hazardPeriods)) {
      columns
    } else {
      union(columns, ncol(data@x))
    }
  }
  firstColumns <- withPeriod(resolveForestVars(firstForest$vars, data))
  singleForest <- is.null(declaredBases)
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
    variedSplitProbs <- varySplitProbabilities(priors$tree.prior)
    unsupportedMultinomial <- c(
      "'split.probs' (a vector that is not uniform, or dart())" = variedSplitProbs,
      "'monotone'" = !is.null(monotoneDirections),
      "a linear leaf prior" = is(priors$leaf.prior, "dbartsLinearPrior"),
      "a Gaussian-process leaf prior" = is(priors$leaf.prior, "dbartsGPPrior"),
      "a 'k' hyperprior" = is(priors$leaf.prior@k, "dbartsLeafHyperprior"),
      "a named leaf-prior 'sd'" = !is.null(priors$leaf.prior@prior.sd),
      "storage = \"single\"" = identical(control@storage, "single")
    )
    if (any(unsupportedMultinomial)) {
      refuseMultinomialStatement(
        names(unsupportedMultinomial)[unsupportedMultinomial]
      )
    }
    if (any(nzchar(names(forests)))) {
      stop(
        "a multinomial model's one forest entry states every category forest ",
        "and takes no name",
        call. = FALSE
      )
    }
  }

  # a single forest's column restriction is a model fact and rides a model
  # attribute the C bridge reads, beside the two constraints above, so a copy
  # or a reload rebuilds it. Naming every column restricts nothing and stores
  # nothing; a multi-forest fit carries every forest's columns in the model's
  # forest records instead. After the multinomial refusals above, so a
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

  # the AFT survival family reads its per-observation status off the data
  # object; the C bridge validates it
  if (!is.null(survivalStatus)) {
    if (length(survivalStatus) != length(data@y)) {
      stop("survival status must have length ", length(data@y))
    }
    data@status <- survivalStatus
  }
  # the discrete-time hazard marker: the period grid, which packageBartResults
  # reads into $periods. The C bridge never reads it, so a hazard fit's draw
  # stream is byte-identical to the by-hand binary fit's.
  model@hazard.periods <- hazardPeriods
  model@n.categories <- numCategories
  model@shape <- shapeSpec

  # the heteroscedastic variance forest: a
  # `variance` selector installs a second forest modeling s^2(x); its config
  # rides the model's variance record, which the C bridge reads. Gaussian + constant leaf
  # only (the C factory refuses otherwise; a friendly R check for the family).
  # `variance` accepts either the plain shorthand (NULL/FALSE/TRUE/formula/
  # character/index) or a varianceForest() object: its `vars`
  # slot routes through the SAME resolveVarianceColumns the shorthand uses -
  # one selector vocabulary - and its n.trees/base/power knobs land on the
  # same record the C factory reads.
  varianceSpec <- if (inherits(variance, "dbartsVarianceForest")) {
    resolveForestMoves(variance, "varianceForest()")
  } else {
    NULL
  }
  refuseVarianceRuleGibbs(varianceSpec$proposal.probs)
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
    varianceRecord <- new("dbartsVarianceRecord")
    varianceRecord@n.trees <- coerceOrError(n.trees, "integer")
    refuseVarianceCalibration(model@resid.prior, varianceRecord@n.trees)
    # what the caller stated, so that $setModel installs a stated depth and
    # leaves one copied from the mean forest
    varianceRecord@stated <- c("n.trees", "power", "base")[c(
      !is.null(varianceNTrees),
      !is.null(variancePower),
      !is.null(varianceBase)
    )]
    # a stated power or base is held to what setVarianceForest holds it to
    varianceRecord@base <- if (is.null(varianceBase)) {
      model@tree.prior@base
    } else {
      validateForestTreePrior(varianceBase, "base")
    }
    varianceRecord@power <- if (is.null(variancePower)) {
      model@tree.prior@power
    } else {
      validateForestTreePrior(variancePower, "power")
    }
    varianceRecord@columns <- if (allColumns) {
      NULL
    } else {
      as.integer(varianceColumns)
    }
    model@variance <- varianceRecord
  }

  # the Bayesian causal forest: a second forest with a
  # two-level factor basis selects the model y = a mu(x) + b_z tau(x) + eps.
  # The block that basis expands to (a factor's level indicators) and the
  # factor's levels are conditioning DATA and ride the data object beside the
  # weights they mirror; a block dbartsSpec() declares was finished there, on
  # the rows kept, and installed above, and a data object handed in whole was
  # checked on entry; the second forest's
  # CONFIGURATION - its tree count and structure prior, its column mask, the
  # amplitude scales and the per-forest constraints - rides the model's forest
  # records the C bridge reads, as the variance forest's does above. The bridge
  # cross-checks the two halves in both directions, so a stripped record is a
  # loud error, never a silent single-forest fit.
  if (!is.null(data@bases)) {
    # the RESOLVED forest count: the bases the call's declarations built,
    # which replace any the data object held
    numForests <- length(data@bases)
    # K = 1 is not a shipped configuration (dec-A109). Every creation route
    # reaches here - a forest() term and a forests = list, whose declarations
    # forestBasisDeclarations carries down at any length - so this is the
    # single site the refusal is owed at. A varying-coefficient model is one
    # forest per coefficient function, intercept included.
    if (numForests < 2L) {
      stop(
        "a multi-forest model needs at least two forests, and this call's ",
        "'basis' declarations resolve to ",
        numForests,
        ": a forest with a 'basis' stands beside another forest. Write ",
        "the forest with no multiplier too, as y ~ forest(x1 + x2) + ",
        "forest(x1 + x2, basis = z1) or forests = list(forest(), ",
        "forest(basis = z1)), or use a single forest with linear() ",
        "leaves; otherwise drop the basis"
      )
    }
    # the families the calibration map has a latent scale to state its node
    # scales against, and whose own parameter block is shown to interleave with
    # the amplitude block. A fixed error scale is what makes the binary
    # families work here rather than a reason they cannot: the combined index
    # is stated in the link's own units and sigma is pinned there.
    if (family %not_in% c("gaussian", "probit", "logistic")) {
      stop(
        "a model of several forests does not support family \"",
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
    # non-constant leaf, a variance forest, an fp32 residual, a per-column cut
    # cap, nor the Student-t error distribution. Every one of those would
    # otherwise be dropped in silence, changing the fitted model without a
    # word; name each one instead. 'proposal.probs' is NOT among them: the
    # chain carries the control's mixture onto every forest it builds, so a
    # coupling honors it rather than dropping it.
    # the split probabilities are the model's, so a model of several forests
    # takes them where they state the default, equal probabilities
    if (varySplitProbabilities(priors$tree.prior)) {
      stop(
        "a model of several forests does not support 'split.probs' (a ",
        "vector that is not uniform, or dart()); drop it or fit a ",
        "single-forest model",
        call. = FALSE
      )
    }
    # every forest with no basis states its size, so none has this k to take
    if (
      !is.null(declaredBases) &&
        !shares &&
        (is(priors$leaf.hyperprior, "dbartsChiHyperprior") ||
          (is(priors$leaf.hyperprior, "dbartsFixedHyperprior") &&
            priors$leaf.hyperprior@k != 2.0))
    ) {
      stop(
        "no forest of this model takes a leaf prior's 'k': each forest with ",
        "no basis states its own size, as forest(sd = ), in the response's ",
        "units; drop 'k'",
        call. = FALSE
      )
    }
    unsupported <- c(
      "'monotone'" = !is.null(monotoneDirections),
      "a linear leaf prior" = is(priors$leaf.prior, "dbartsLinearPrior"),
      "a Gaussian-process leaf prior" = is(priors$leaf.prior, "dbartsGPPrior"),
      # a forest with no basis stating no sd takes the fitting function's
      # leaf prior as one forest does; with none it has no forest
      "a 'k' hyperprior" = !shares &&
        is(priors$leaf.hyperprior, "dbartsChiHyperprior"),
      "a non-default 'k'" = !shares &&
        is(priors$leaf.hyperprior, "dbartsFixedHyperprior") &&
        priors$leaf.hyperprior@k != 2.0,
      # "differs from the family default": defaultLeafScale(family), not a
      # gaussian-only literal
      "a non-default 'leaf.scale'" = model@leaf.scale !=
        defaultLeafScale(family),
      # the calibration map fixes every forest's leaf scale from the family's
      # own latent scale, so a named anchor has nowhere to land and the
      # leaf.scale gate above does not fire on it; parsePriors refuses a named
      # sd first, so this backstops a model built by hand
      "a named leaf-prior 'sd'" = !shares && !is.na(model@prior.scale),
      "Student-t residuals" = !is.null(residDf),
      "'variance'" = !is.null(varianceColumns),
      "storage = \"single\"" = identical(control@storage, "single"),
      "per-column 'n.cuts'" = length(unique(data@n.cuts)) > 1L,
      "test predictors" = !is.null(data@x.test)
    )
    if (any(unsupported)) {
      stop(
        "a model of several forests does not support ",
        paste0(names(unsupported)[unsupported], collapse = ", "),
        "; drop it or fit a single-forest model"
      )
    }
    # a drawn coefficient stating no 'sd' takes its trees' size from the leaf
    # prior, which a flat prior on k or an infinite k does not, nor a
    # chi prior whose median the engine's quantile rounds to zero (it tests
    # the size for positive and finite, the quotient of the prior's scale by
    # that median)
    hyperprior <- priors$leaf.hyperprior
    if (leafPriorGivesNoSize(hyperprior)) {
      for (index in seq_along(forestSpec)) {
        spec <- forestSpec[[index]]
        if (
          statesCoefficientPrior(spec$amplitude) &&
            is.null(spec$sd) &&
            index %in% plain
        ) {
          refuseSizelessCoefficient(
            forestKnobLabel(knobLabels, index),
            spec$amplitude,
            hyperprior
          )
        }
      }
    }
    # every forest resolves its own knobs; a data object carrying bases with no
    # forests = declaration at all resolves to the same defaults
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
    params <- forestParams(specs, hasBasis, fitTree)
    statements <- forestStatements(specs, hasBasis)
    treeCounts <- vapply(params, function(forest) as.integer(forest[1L]), 0L)
    forestColumns <- lapply(
      seq_len(numForests),
      function(index) {
        if (index == 1L) {
          firstColumns
        } else {
          withPeriod(resolveForestVars(specs[[index]]$vars, data))
        }
      }
    )
    # every forest's label, fixed here for the life of the sampler: a list's
    # own name, else the text of a basis written as code, else its position
    labels <- forestLabels(names(forests), basisTexts, numForests)
    refuseConstantBasis(data, allRows, specs, labels)
    # the first forest's constraints are the model's, already resolved
    # above; the rest take their own, and the plain forest the fitting
    # function's where it states none, each resolved against the columns
    # and the tree count of the forest it lands on
    forestInteractions <- c(
      list(interactionSpec),
      lapply(
        seq_len(numForests)[-1L],
        function(index) {
          own <- specs[[index]]$interactions
          resolveInteractions(
            if (is.null(own) && identical(index, plain)) {
              fitInteractions
            } else {
              own
            },
            data
          )
        }
      )
    )
    forestBlocks <- c(
      list(blockSpec),
      lapply(
        seq_len(numForests)[-1L],
        function(index) {
          own <- specs[[index]]$blocks
          resolveBlocks(
            if (is.null(own) && identical(index, plain)) fitBlocks else own,
            data,
            treeCounts[index],
            availableColumns = forestColumns[[index]]
          )
        }
      )
    )
    # one record per forest: what it runs under (see forestParams()), what it
    # states of its size and its coefficient's prior as written, NA where it
    # states none (the engine holds the calibration map and its defaults, and
    # reads what each forest multiplies from the data object), its resolved
    # 1-based columns or NULL for unrestricted, its label (see
    # forestLabels()), whether its coefficient is held (amplitude = fixed(),
    # which a fit records in fit$fixed), its constraints, and what builds a
    # basis written as code again at new rows (NULL where a forest's basis is
    # a value, which has no code to build from; inert to the run)
    # each sharer's part of the total, equal, as the engine derives it, NA on
    # every other forest
    shares <- sharerShares(hasBasis, statements)
    records <- lapply(seq_len(numForests), function(index) {
      record <- new("dbartsForestRecord")
      record@n.trees <- treeCounts[index]
      record@stated <- c(
        forestCountStated(specs[[index]], index %in% plain, stated),
        forestDepthStated(specs[[index]])
      )
      record@base <- params[[index]][2L]
      record@power <- params[[index]][3L]
      record@held.amplitude <- identical(specs[[index]]$amplitude, "fixed")
      record@sd <- statements$sd[index]
      if (!is.null(statements$amplitude.df)) {
        record@amplitude.df <- statements$amplitude.df[index]
      }
      record@share <- shares[index]
      record@vars <- forestColumns[[index]]
      record@interactions <- forestInteractions[[index]]
      record@blocks <- forestBlocks[[index]]
      record@label <- labels[index]
      # a held basis kept by its text has the text and no terms
      if (
        index <= length(basisRecords) && !is.null(basisRecords[[index]]$terms)
      ) {
        record@basis.term <- basisRecords[[index]]
      }
      record
    })
  } else {
    records <- list(plainForestRecord(
      control@n.trees,
      forestCountStated(firstForest, TRUE, stated)
    ))
  }

  # each forest's own move settings, where any forest states one; the
  # control's slots keep the control's, which every other forest runs
  moves <- forestMoveRecord(
    forestSpec,
    length(records),
    varianceSpec
  )
  if (!is.null(moves)) {
    for (index in seq_along(records)) {
      records[[index]]@moves <- moves$forests[[index]]
    }
    if (!is.null(moves$variance)) {
      model@variance@moves <- moves$variance
    }
  }
  model@forests <- records
  control <- mirrorControl(control, model, data)
  warnNotUsed(
    carried$values,
    heldCopies(model, data),
    notUsedWords(
      carried$door,
      model,
      data,
      !isTRUE(attr(stated, "noForestRemedy", exact = TRUE))
    )
  )

  # every resolution of "auto" is announced once, here, after the refusals
  # above: bart() reaches this through dbarts(), and dbartsSpec() reads its
  # control's verbose
  if (!is.null(autoDescription)) {
    announceAutoFamily(control@verbose, family, autoDescription)
  }

  namedList(control, model, data, family)
}

## What a multinomial model's creation and $setForest refuse, by its name.
refuseMultinomialStatement <- function(named) {
  stop(
    "a multinomial (softmax) model does not support ",
    paste0(named, collapse = ", "),
    "; drop it or fit a single-forest model",
    call. = FALSE
  )
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
  power = 2,
  base = 0.95,
  split.probs = NULL,
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
  n.trees = NULL,
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
  weights.type = NULL,
  survival = NULL,
  tree.prior = NULL,
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
  # the retired 'tree.prior' is read as the plain arguments it states, which
  # are what this call states of them from here on
  retiredTreePrior <- readRetiredTreePrior(
    matchedCall,
    parentEnv,
    "dbartsSpec"
  )
  matchedCall <- withRetiredTreePrior(matchedCall, retiredTreePrior)
  # judged without forcing: power may hold a retired tree prior in its
  # vocabulary
  for (name in nullWritten(matchedCall, c("power", "base"), parentEnv)) {
    refuseNullArgument(NULL, name)
  }
  sigest <- resolveSigestArg(sigest, "dbartsSpec", "refuse")
  resolvedSigest <- resolveSigestRule(sigest, "dbartsSpec")
  sigest <- resolvedSigest$value
  sigestRule <- resolvedSigest$rule

  if (!inherits(data, "dbartsData")) {
    stop("'data' must be a dbartsData object; see ?dbartsData")
  }
  # as handed over, which a forest's basis is read against
  handedData <- data
  data <- resolveDataObject(data)
  if (!inherits(control, "dbartsControl")) {
    stop("'control' must be a dbartsControl object; see ?dbartsControl")
  }
  # the tree-move mixture is the control's, as at dbarts() and bart(): the
  # retired flat spelling is used with a warning, and refused beside a control
  # whose own call named it
  if ("proposal.probs" %in% names(matchedCall)) {
    warnConsolidated("proposal.probs", "dbartsSpec")
    refuseCollidingMixture(control)
  }
  # the count as this call states it, flat or written on the control; the cut
  # settings have no flat spelling here
  copies <- readDoorCopies(
    releaseCarriedControl(control, "dbartsSpec"),
    list(n.trees = n.trees),
    names(matchedCall)
  )
  control <- copies$control
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

  # read as dbarts() reads them beside the same data object
  data <- applyCutSettings(data, copies)
  # an explicit sigest overrides whatever the data carries, a number by
  # replacing it and "dense" or "sparse" by estimating over it; NULL and
  # "auto" leave it alone, so a consumer's own starting estimate survives (an
  # unset one is estimated during resolution, exactly as for dbarts())
  if (!is.na(sigest)) {
    data@sigma <- validateSigest(sigest, "dbartsSpec")
  } else if (sigestRule != "auto") {
    data@sigma <- NA_real_
  }

  # as on dbarts(): the forest constructors resolve by bare name inside their
  # arguments, here in parentEnv
  forestArguments <- resolveForestArguments(matchedCall, parentEnv)
  forests <- forestArguments$forests
  interactions <- forestArguments$interactions
  blocks <- forestArguments$blocks
  monotone <- forestArguments$monotone
  variance <- forestArguments$variance

  # a declared basis is read against the data object as dbarts() reads one
  # beside it (readFitBasis) and built on the rows the object keeps. A forests
  # list states the forests in full: what the object held is replaced, a held
  # basis a forest names kept, and a list of one forest with none, or none,
  # gives one forest
  basisReads <- NULL
  if (!is.null(forestBasisDeclarations(forests))) {
    basisSource <- objectBasisSource(handedData, forests)
    declared <- readDeclaredBases(forests, basisSource)
    forests <- declared$forests
    basisReads <- declared$reads
    if (all(vapply(basisReads, is.null, NA)) && length(basisReads) < 2L) {
      basisReads <- NULL
    }
  }
  built <- if (!is.null(basisReads)) {
    buildFitBases(basisReads, basisSource, data@rows, length(data@y), forests)
  }
  data <- installFitBases(data, built)
  basisRecords <- fitBasisRecords(built)

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
    forests = forests,
    evalEnv = parentEnv,
    residPrior = residPrior,
    familySpec = familySpec,
    basisRecords = basisRecords,
    sigestRule = sigestRule,
    treePrior = retiredTreePrior$record,
    weightsType = validateWeightsType(weights.type)
  )
}
