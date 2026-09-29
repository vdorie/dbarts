ORDINAL_VARIABLE <- 0L
CATEGORICAL_VARIABLE <- 1L
ORDERED_FACTOR_VARIABLE <- 2L

## The package's own na.action, and the default of every fitting function
## that takes one. It drops the rows whose RESPONSE is missing and keeps the
## rows whose predictors are: BART routes a missing predictor value down a
## learned side of each rule, so an incomplete row is data rather than a
## hole, while a missing response is nothing to fit. rpart's na.rpart is the
## behavioural precedent. The dropped rows are recorded exactly as
## stats::na.exclude records them - a named integer vector of class
## "exclude" - so training fits pad back to the caller's own row count
## through stats::naresid.
##
## `object` is a model frame; the response is the column the frame's terms
## attribute names, and a frame with no response loses no rows at all.
na.keepPredictors <- function(object, ...) {
  response <- attr(attr(object, "terms"), "response")
  if (is.null(response) || length(response) != 1L || response == 0L) {
    return(object)
  }
  ## is.na() on a multi-column response is a matrix, except where a class
  ## defines its own (survival::Surv reduces to one value per row); the
  ## reduction keys off what came back, not off the response's own shape
  dropped <- is.na(object[[response]])
  if (!is.null(dim(dropped))) {
    dropped <- apply(dropped, 1L, any)
  }
  if (!any(dropped)) {
    return(object)
  }
  omit <- seq_along(dropped)[dropped]
  names(omit) <- rownames(object)[dropped]
  class(omit) <- "exclude"
  kept <- object[!dropped, , drop = FALSE]
  attr(kept, "na.action") <- omit
  kept
}

## Before R 4.3, terms.formula() runs as.data.frame() over a plain list
## 'data', which refuses a classed matrix that no as.data.frame method claims
## (a Surv-like response built without survival). Such a list is returned as
## a data frame holding the column as is, which terms() then leaves alone;
## any other list is returned unchanged.
asDataFrameableList <- function(data) {
  lengths <- vapply(data, NROW, 1L)
  bare <- vapply(
    data,
    function(column) {
      is.matrix(column) &&
        is.object(column) &&
        !any(vapply(
          class(column),
          function(k) !is.null(utils::getS3method("as.data.frame", k, TRUE)),
          NA
        ))
    },
    NA
  )
  if (!any(bare) || length(unique(lengths)) != 1L) {
    return(data)
  }
  structure(
    data,
    class = "data.frame",
    row.names = .set_row_names(lengths[1L])
  )
}

## Which rows of a predictor container hold a missing value. Written off the
## stored entries for the sparse flavors: an implicit zero is an observed
## value, so densifying to look for NAs would be both wasteful and wrong.
rowsWithMissingPredictors <- function(x) {
  n <- NROW(x)
  if (is.matrix(x) || is.data.frame(x)) {
    if (!anyNA(x)) {
      return(rep_len(FALSE, n))
    }
    return(rowSums(is.na(x)) > 0L)
  }
  if (inherits(x, "dbartsMixedMatrix")) {
    rows <- rep_len(FALSE, n)
    for (column in if (is.null(x$dense)) list() else x$dense) {
      if (anyNA(column)) {
        rows <- rows | is.na(column)
      }
    }
    if (!is.null(x$sparse)) {
      rows <- rows | rowsWithMissingPredictors(x$sparse)
    }
    return(rows)
  }
  if (inherits(x, "dgCMatrix")) {
    rows <- rep_len(FALSE, n)
    missingEntries <- is.na(x@x)
    if (any(missingEntries)) {
      rows[x@i[missingEntries] + 1L] <- TRUE
    }
    return(rows)
  }
  if (!anyNA(x)) {
    return(rep_len(FALSE, n))
  }
  is.na(as.vector(x))
}

## Applies an 'na.action' to the (y, x) pair the matrix interface supplies.
## Those functions take a MODEL FRAME, so the pair is presented as one: the
## response, and a single predictor column that is NA exactly on the rows
## where some predictor is. na.omit, na.exclude, na.fail and na.pass then
## mean on this interface what they mean on the formula one. Returns the
## rows to keep and the record of what was dropped, or NULL when nothing is
## missing at all and no na.action can have anything to say.
applyNaActionToXY <- function(na.action, y, x) {
  predictorNA <- rowsWithMissingPredictors(x)
  responseNA <- is.na(y)
  if (!any(predictorNA) && !any(responseNA)) {
    return(NULL)
  }
  frame <- data.frame(
    response = ifelse(responseNA, NA_real_, 0.0),
    predictors = ifelse(predictorNA, NA_real_, 0.0)
  )
  kept <- stats::model.frame(
    response ~ predictors,
    frame,
    na.action = na.action
  )
  omit <- attr(kept, "na.action")
  keep <- rep_len(TRUE, length(responseNA))
  if (!is.null(omit)) {
    keep[unclass(omit)] <- FALSE
  }
  list(keep = keep, na.action = omit)
}

## A multinomial fit's true response is 'counts' (an n x K matrix), not the
## 'y' the matrix-interface branches otherwise derive - a placeholder trials
## vector with nothing of its own to be missing. Folding a row's missingness
## in there routes it through applyNaActionToXY exactly as a single-column
## response's own NA would: a count-matrix row is missing when any of its
## cells is NA (resolveMultinomialCounts marks a missing label the same way).
multinomialResponseNA <- function(counts, y) {
  if (is.null(counts) || NROW(counts) != length(y)) {
    return(y)
  }
  ifelse(rowSums(is.na(counts)) > 0L, NA_real_, y)
}

## applyNaActionToXY's record names its rows by position in a synthetic
## frame; relabel them with the caller's own row names when there are any, as
## a model frame's record would be, so padding fills those names back in.
nameOmittedRows <- function(omit, rowNames) {
  if (!is.null(omit) && !is.null(rowNames)) {
    names(omit) <- rowNames[unclass(omit)]
  }
  omit
}

## The rows a formula fit kept, in the caller's own row numbering, once the
## model frame's na.action has taken its share: 'subset' chose them and the
## na.action then dropped some by position within that choice. Anything the
## caller supplied at the full, pre-'subset' shape - a forest's amplitude
## basis - is restricted through this.
alignSubsetRowsToFrame <- function(subsetRows, naOmitted, keptRows) {
  if (is.null(naOmitted)) {
    return(subsetRows)
  }
  if (is.null(subsetRows)) {
    # no 'subset' at all, so the frame's rows ARE the data's - less whatever
    # the na.action dropped, which is still a restriction a full-data basis
    # has to follow
    full <- keptRows + length(naOmitted)
    return(list(full = full, index = seq_len(full)[-unclass(naOmitted)]))
  }
  subsetRows$index <- subsetRows$index[-unclass(naOmitted)]
  subsetRows
}

## The same restriction on the matrix interface, where the bases were
## already validated and subset against the caller's own row count before
## the na.action ran: forest f's basis loses exactly the rows y did.
restrictBasesToRows <- function(bases, keep) {
  if (is.null(bases) || all(keep)) {
    return(bases)
  }
  lapply(bases, function(basis) {
    if (is.null(basis)) NULL else basis[keep, , drop = FALSE]
  })
}

## An out-of-range 'subset' is silent in base R: row indexing pads an
## unmatched row with NA rather than erroring, and the na.action would then
## drop exactly those rows and fit fewer observations than the caller asked
## for without a word. Named here instead, ahead of the model frame. `n` of
## NA means the row count is not known this early, which leaves the check
## unrun rather than guessed at.
refuseOutOfRangeSubset <- function(index, n, rowNames = NULL) {
  if (is.null(index) || is.na(n)) {
    return(invisible(NULL))
  }
  outOfRange <- if (is.logical(index)) {
    length(index) > n
  } else if (is.numeric(index)) {
    anyNA(index) || any(abs(index) > n)
  } else if (is.character(index)) {
    anyNA(index) || any(index %not_in% rowNames)
  } else {
    FALSE
  }
  if (outOfRange) {
    stop(
      "'subset' selects rows outside the data, which has ",
      n,
      " rows; check that 'subset' selects rows within range"
    )
  }
  invisible(NULL)
}

## Pads a training-side quantity back to the caller's own row count. A
## na.action of class "exclude" records the rows it dropped and naresid puts
## NA back in their places; "omit" records them and pads nothing, which is
## exactly na.omit's contract. The observation margin is the LAST one in
## every draws array this package reports, and naresid pads the first, so
## only vectors and observation-by-column matrices go through here.
padOmittedRows <- function(naOmitted, x) {
  if (is.null(naOmitted) || is.null(x)) {
    return(x)
  }
  stats::naresid(naOmitted, x)
}

## The row names of a predictor set as given, before any coding: a data
## frame's automatic names count as "1".."n", and a matrix, sparse matrix or
## vector without them gives NULL, as lm.fit does.
observationRowNames <- function(x) {
  if (is.data.frame(x) || is.matrix(x) || inherits(x, "dgCMatrix")) {
    rownames(x)
  } else {
    NULL
  }
}

## One channel ("train" or "test") of a data object's row-name record. A data
## object saved before the slot existed has none.
dataRowNames <- function(data, channel) {
  if (!methods::.hasSlot(data, "rowNames") || is.null(data@rowNames)) {
    return(NULL)
  }
  data@rowNames[[channel]]
}

## Stores one channel's names on a data object, keeping the slot NULL when
## neither channel is named.
setDataRowNames <- function(data, channel, names) {
  rowNames <- if (methods::.hasSlot(data, "rowNames")) data@rowNames
  if (is.null(rowNames)) {
    rowNames <- list(train = NULL, test = NULL)
  }
  rowNames[channel] <- list(names)
  data@rowNames <- if (is.null(rowNames$train) && is.null(rowNames$test)) {
    NULL
  } else {
    rowNames
  }
  data
}

## Names the observation margin of an output: a vector's elements, or
## otherwise the margin 'trailing' places before the last (0 when the
## observations are last, 1 when a category, forest or interval axis trails).
## Set in place, so an unshared array is not copied.
nameObservationMargin <- function(x, names, trailing = 0L) {
  if (is.null(names) || is.null(x)) {
    return(x)
  }
  d <- dim(x)
  if (is.null(d)) {
    if (length(x) != length(names)) {
      stop("internal error: row names do not match the observations")
    }
    names(x) <- names
    return(x)
  }
  margin <- length(d) - trailing
  if (d[margin] != length(names)) {
    stop("internal error: row names do not match the observations")
  }
  dn <- dimnames(x)
  if (is.null(dn)) {
    dn <- vector("list", length(d))
  }
  dn[margin] <- list(names)
  dimnames(x) <- dn
  x
}

# The multinomial capability probe: a data object carrying the n x K count
# response is a multinomial one, on both the fitting and the mutation surfaces.
# A capability test rather than a forest count, for the reason
# samplerCarriesAmplitudes is one.
dataCounts <- function(data) {
  data@counts
}

methods::setMethod(
  "initialize",
  "dbartsData",
  function(.Object, modelMatrices, n.cuts = 100L, sigma = NA_real_) {
    if (!missing(modelMatrices)) {
      .Object@y <- modelMatrices$y
      .Object@x <- modelMatrices$x
      # makeCategoricalModelMatrix types its columns; everything else is ordinal
      .Object@varTypes <- if (!is.null(attr(.Object@x, "varTypes"))) {
        as.integer(attr(.Object@x, "varTypes"))
      } else {
        rep.int(ORDINAL_VARIABLE, ncol(.Object@x))
      }
      .Object@x.test <- modelMatrices$x.test
      .Object@weights <- modelMatrices$weights
      .Object@weights.test <- modelMatrices$weights.test
      .Object@offset <- modelMatrices$offset
      .Object@offset.test <- modelMatrices$offset.test
      .Object@bases <- modelMatrices$bases
      .Object@counts <- modelMatrices$counts
      .Object@offset.category <- modelMatrices$offset.category
      .Object@offset.category.test <- modelMatrices$offset.category.test

      .Object@testUsesRegularOffset <- modelMatrices$testUsesRegularOffset
    }

    .Object@n.cuts <- rep_len(as.integer(n.cuts), ncol(.Object@x))
    .Object@sigma <- sigma

    validObject(.Object)
    .Object
  }
)

makeTestModelMatrix <- function(data, newdata) {
  validateXTest(newdata, data@x)
}

## A split rule learns a route for NA only on a column whose TRAINING values
## carried one: the missing direction is drawn only there and cannot be
## restored onto an NA-free column, so on a training-complete column every
## rule sends NA down one fixed branch. Refuse rather than answer from a
## route the model never learned.
sourceColumnHasNA <- function(source, j, numColumns, numObservations) {
  column <- predictorSourceColumn(source, j, numColumns, numObservations)
  if (is.list(column)) {
    return(anyNA(column$x) || is.na(column$implicit))
  }
  anyNA(column)
}

sourceHasNA <- function(source) {
  if (inherits(source, "dbartsMixedMatrix")) {
    return(
      any(vapply(source$dense, anyNA, logical(1L))) ||
        (!is.null(source$sparse) && anyNA(source$sparse@x))
    )
  }
  if (inherits(source, "dgCMatrix")) {
    return(anyNA(source@x))
  }
  anyNA(source)
}

## The predictor columns a test set leaves without a route: missing there, and
## complete in training.
unroutableTestColumns <- function(x.test, x.train) {
  # the whole-object probe short-circuits, so complete test data - the usual
  # case - pays one scan and never touches the training side
  if (!sourceHasNA(x.test)) {
    return(integer(0L))
  }
  # validateXTest has already matched the two sides' column counts by here,
  # so one count serves both
  numColumns <- NCOL(x.test)
  numTest <- NROW(x.test)
  numTrain <- NROW(x.train)
  offending <- integer(0L)
  for (j in seq_len(numColumns)) {
    if (!sourceColumnHasNA(x.test, j, numColumns, numTest)) {
      next
    }
    if (sourceColumnHasNA(x.train, j, numColumns, numTrain)) {
      next
    }
    offending <- c(offending, j)
  }
  offending
}

## The test rows missing a value in any of 'columns', read off the stored
## entries so that a sparse container is never densified.
testRowsMissingIn <- function(x.test, columns) {
  numTest <- NROW(x.test)
  numColumns <- NCOL(x.test)
  rows <- rep_len(FALSE, numTest)
  for (j in columns) {
    column <- predictorSourceColumn(x.test, j, numColumns, numTest)
    if (!is.list(column)) {
      rows <- rows | is.na(column)
    } else if (is.na(column$implicit)) {
      missingRows <- rep_len(TRUE, numTest)
      missingRows[column$i + 1L] <- is.na(column$x)
      rows <- rows | missingRows
    } else {
      rows[column$i[is.na(column$x)] + 1L] <- TRUE
    }
  }
  rows
}

unroutableTestRows <- function(x.test, x.train) {
  testRowsMissingIn(x.test, unroutableTestColumns(x.test, x.train))
}

## Labels predictor columns by name, or by position when training had none.
testColumnLabels <- function(columns, x.train) {
  predictorNames <- colnames(x.train)
  labels <- if (is.null(predictorNames)) {
    paste0("column ", columns)
  } else {
    paste0("'", predictorNames[columns], "'")
  }
  shown <- labels[seq_len(min(5L, length(labels)))]
  paste0(
    toString(shown),
    if (length(labels) > 5L) {
      paste0(" and ", length(labels) - 5L, " more column(s)")
    }
  )
}

## 'naActionHint' is set by the predict methods, whose callers can choose
## another na.action; the sampler's methods keep the default and say nothing.
refuseTestMissingness <- function(x.test, x.train, naActionHint = FALSE) {
  offending <- unroutableTestColumns(x.test, x.train)
  # every NA sits in a column that carried training NAs: every one has a
  # learned route, and nothing is refused
  if (length(offending) == 0L) {
    return(invisible(NULL))
  }
  stop(
    "test predictors have missing values in ",
    testColumnLabels(offending, x.train),
    ", which carried none in training: a split rule learns a route for NA ",
    "only on a column that had missing values when the trees were grown, ",
    "so these rows have no route to take",
    if (naActionHint) {
      paste0(
        "; use na.action = na.pass to return NA for them, or na.omit to ",
        "drop them"
      )
    }
  )
}

## The predict methods' 'na.action' (dec-B34): a function or its name, with
## NULL meaning the default, as predict.lm reads it - the option is never
## consulted.
resolvePredictNaAction <- function(na.action) {
  if (is.null(na.action)) {
    return(na.keepPredictors)
  }
  if (is.character(na.action) && length(na.action) == 1L) {
    return(match.fun(na.action))
  }
  if (!is.function(na.action)) {
    stop("'na.action' must be a function, the name of one, or NULL")
  }
  na.action
}

## A per-row 'offset' or 'weights' channel given to predict has no route an
## na.action can assign it but its own (dec-A89): a length-one value recycles
## and is out of scope, unconditionally refused if NA since there is then no
## row for na.action to drop; anything else must match newdata's rows, and a
## row with a missing cell anywhere in it (a flat NA, or any NA in a category
## offset's row) joins the incomplete-row set resolvePredictRows resolves.
## Returns a logical vector, or NULL when the channel is absent or complete.
predictChannelIncomplete <- function(value, n, argument) {
  if (is.null(value)) {
    return(NULL)
  }
  if (is.null(dim(value)) && length(value) == 1L) {
    if (anyNA(value)) {
      stop("'", argument, "' has missing values")
    }
    return(NULL)
  }
  if (NROW(value) != n) {
    stop("'", argument, "' must have the same number of rows as 'newdata'")
  }
  incomplete <- if (is.null(dim(value))) {
    is.na(value)
  } else {
    rowSums(is.na(value)) > 0L
  }
  if (!any(incomplete)) NULL else incomplete
}

## The refusal a 'channels' entry earns when na.action keeps every row
## despite it (dec-A89): unlike a predictor, an offset or weight is never
## routable, so a survivor is refused outright, with the same hint
## refuseTestMissingness gives a predictor's.
refuseExtraMissingness <- function(extra, naActionHint) {
  for (name in names(extra)) {
    if (any(extra[[name]])) {
      stop(
        "'",
        name,
        "' has missing values",
        if (naActionHint) {
          paste0(
            "; use na.action = na.pass to return NA for them, or na.omit ",
            "to drop them"
          )
        }
      )
    }
  }
  invisible(NULL)
}

## Joins quoted names in prose: one alone, or every name but the last
## comma-separated with "and" before it.
quotedNameList <- function(names) {
  quoted <- paste0("'", names, "'")
  if (length(quoted) <= 1L) {
    return(quoted)
  }
  paste0(
    paste(quoted[-length(quoted)], collapse = ", "),
    " and ",
    quoted[length(quoted)]
  )
}

## na.fail's refusal, naming whichever of the predictor columns and the
## 'extra' channels (dec-A89) carried a missing value.
refuseNaFail <- function(x.test, x.train, extra) {
  columns <- Filter(
    function(j) any(testRowsMissingIn(x.test, j)),
    seq_len(NCOL(x.test))
  )
  channels <- Filter(function(name) any(extra[[name]]), names(extra))
  parts <- character(0L)
  if (length(columns) > 0L) {
    parts <- c(
      parts,
      paste0(
        "test predictors have missing values in ",
        testColumnLabels(columns, x.train)
      )
    )
  }
  if (length(channels) > 0L) {
    parts <- c(
      parts,
      paste0(
        quotedNameList(channels),
        if (length(channels) > 1L) {
          " have missing values"
        } else {
          " has missing values"
        }
      )
    )
  }
  stop(paste(parts, collapse = "; "), ", which na.action = na.fail refuses")
}

## Which rows of a coded test set a predict method answers (dec-B34). The
## incomplete-row set is the union of an unroutable predictor and a missing
## 'extra' channel (offset, weights; dec-A89), named by 'channels' in
## preparePredictRows. na.pass drops only those rows and pads them back;
## na.fail refuses any of them, by name; every other function is applied to a
## one-column frame that is NA on the incomplete rows, as applyNaActionToXY
## does for the matrix interface, and whatever it keeps then meets the
## default's refusal. Returns the logical 'keep' and whether the dropped rows
## 'pad' back as NA, or NULL when every row is answered as it stands.
resolvePredictRows <- function(na.action, x.test, x.train, extra = list()) {
  numTest <- NROW(x.test)
  extraDropped <- Reduce(`|`, extra, rep_len(FALSE, numTest))
  if (identical(na.action, stats::na.pass)) {
    dropped <- unroutableTestRows(x.test, x.train) | extraDropped
    if (!any(dropped) && numTest > 0L) {
      return(NULL)
    }
    return(list(keep = !dropped, pad = TRUE))
  }
  hasNA <- sourceHasNA(x.test) || any(extraDropped)
  if (identical(na.action, stats::na.fail)) {
    if (hasNA) {
      refuseNaFail(x.test, x.train, extra)
    }
    return(if (numTest == 0L) list(keep = logical(0L), pad = FALSE))
  }
  if (!hasNA) {
    # a function passes a complete frame through unchanged
    return(if (numTest == 0L) list(keep = logical(0L), pad = FALSE))
  }
  incomplete <- testRowsMissingIn(x.test, seq_len(NCOL(x.test))) | extraDropped
  frame <- data.frame(predictors = ifelse(incomplete, NA_real_, 0.0))
  kept <- stats::model.frame(~predictors, frame, na.action = na.action)
  omit <- attr(kept, "na.action")
  # read off the kept rows' names, so a function that drops rows without
  # recording them is still followed
  keep <- rownames(frame) %in% rownames(kept)
  # the default's refusal, on the rows kept
  if (all(keep)) {
    refuseTestMissingness(x.test, x.train, naActionHint = TRUE)
    refuseExtraMissingness(extra, naActionHint = TRUE)
    return(NULL)
  }
  refuseTestMissingness(
    x.test[keep, , drop = FALSE],
    x.train,
    naActionHint = TRUE
  )
  refuseExtraMissingness(
    lapply(extra, `[`, keep),
    naActionHint = TRUE
  )
  list(keep = keep, pad = inherits(omit, "exclude"))
}

## The record padPredictedRows reads, from resolvePredictRows's answer over
## 'names', the rows' names.
keptRowsRecord <- function(resolved, names) {
  keep <- resolved$keep
  list(
    n = length(keep),
    keep = keep,
    pad = resolved$pad,
    names = names,
    keptNames = if (any(keep)) names[keep],
    placeholder = !any(keep)
  )
}

## The shared front of every predict method: codes newdata once, resolves its
## rows under 'na.action' - which a missing 'channels' entry (offset,
## weights) marks incomplete exactly as an unroutable predictor does
## (dec-A89) - and returns the rows to predict. When none survives, or
## newdata has none, the fit's first training row stands in and the caller
## slices its answer away, so the result keeps the draw dimensions (dec-B34).
preparePredictRows <- function(newdata, x.train, na.action, channels = NULL) {
  if (missing(newdata) || is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
  rowNames <- observationRowNames(newdata)
  x <- validateXTest(newdata, x.train, refuseMissing = FALSE)
  extra <- list()
  for (name in names(channels)) {
    incomplete <- predictChannelIncomplete(channels[[name]], NROW(x), name)
    if (!is.null(incomplete)) {
      extra[[name]] <- incomplete
    }
  }
  resolved <- resolvePredictRows(
    resolvePredictNaAction(na.action),
    x,
    x.train,
    extra
  )
  if (is.null(resolved)) {
    return(list(x = x, newdata = newdata, n = NROW(x), keptNames = rowNames))
  }
  result <- keptRowsRecord(resolved, rowNames)
  if (result$placeholder) {
    result$x <- suppressPositionalWarnings(validateXTest(
      x.train[1L, , drop = FALSE],
      x.train,
      refuseMissing = FALSE
    ))
  } else {
    result$x <- x[result$keep, , drop = FALSE]
    result$newdata <- if (is.null(dim(newdata))) {
      newdata
    } else {
      newdata[result$keep, , drop = FALSE]
    }
  }
  result
}

suppressPositionalWarnings <- function(expr) {
  withCallingHandlers(
    expr,
    dbartsPositionalArgsWarning = function(w) invokeRestart("muffleWarning")
  )
}

## A per-row input at the rows preparePredictRows kept: a length-one value
## recycles and passes through, anything else must match newdata's rows. A
## placeholder row takes 'stub'. Missingness in a per-row value is resolved
## upstream, in preparePredictRows's 'channels' (dec-A89): a row that reaches
## here survived na.action, so it carries none of its own.
subsetPredictInput <- function(value, rows, argument, stub = NULL) {
  if (is.null(rows$keep)) {
    return(value)
  }
  if (rows$placeholder) {
    return(stub)
  }
  if (is.null(value)) {
    return(NULL)
  }
  if (is.null(dim(value)) && length(value) == 1L) {
    return(value)
  }
  if (NROW(value) != rows$n) {
    stop(
      "'",
      argument,
      "' must have the same number of rows as 'newdata'"
    )
  }
  if (is.null(dim(value))) {
    value[rows$keep]
  } else {
    value[rows$keep, , drop = FALSE]
  }
}

## The same for 'bases', a bare value or one entry per forest.
subsetPredictBases <- function(bases, rows, stub) {
  if (!is.list(bases) || is.data.frame(bases) || !isFALSE(rows$placeholder)) {
    return(subsetPredictInput(bases, rows, "bases", stub))
  }
  lapply(bases, subsetPredictInput, rows = rows, argument = "bases")
}

## Puts a predicted result back on newdata's rows: under na.exclude and
## na.pass the dropped rows return as NA in place, under anything else they
## stay dropped, and a placeholder's answer is cut to zero width or all NA.
## The row margin is the first ('first') or the one 'trailing' places before
## the last. The "s" attribute, laid out as draws by rows, goes the same way.
padPredictedRows <- function(x, rows, trailing = 0L, first = FALSE) {
  if (is.null(rows$keep) || (!rows$pad && !rows$placeholder)) {
    return(x)
  }
  index <- if (rows$pad) {
    index <- rep_len(NA_integer_, rows$n)
    if (!rows$placeholder) {
      index[rows$keep] <- seq_len(sum(rows$keep))
    }
    index
  } else {
    integer(0L)
  }
  names <- if (rows$pad) rows$names else rows$names[rows$keep]
  s <- attr(x, "s")
  d <- dim(x)
  if (is.null(d)) {
    x <- x[index]
    names(x) <- names
  } else {
    margin <- if (first) 1L else length(d) - trailing
    indices <- rep(list(quote(expr = )), length(d))
    indices[[margin]] <- index
    x <- do.call(`[`, c(list(x), indices, list(drop = FALSE)))
    dn <- dimnames(x)
    if (!is.null(dn) || !is.null(names)) {
      if (is.null(dn)) {
        dn <- vector("list", length(d))
      }
      dn[margin] <- list(names)
      dimnames(x) <- dn
    }
  }
  if (!is.null(s)) {
    attr(x, "s") <- padPredictedRows(s, rows)
  }
  x
}

## Leaves R's generator exactly as it was found, removing '.Random.seed' again
## when it was absent: a placeholder row's answer is discarded, and so is any
## draw it made. Returns the restoring function for on.exit.
protectRandomSeed <- function() {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  function() {
    if (had) {
      assign(".Random.seed", seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
}

validateXTest <- function(x.test, x.train, refuseMissing = TRUE) {
  termLabels <- attr(x.train, "term.labels")
  numPredictors <- ncol(x.train)
  predictorNames <- colnames(x.train)
  drop <- attr(x.train, "drop")
  factorLevels <- attr(x.train, "factor.levels")
  varTypes <- attr(x.train, "varTypes")

  if (is.null(x.test)) {
    return(x.test)
  }
  if (is.numeric(x.test) && is.null(dim(x.test)) && length(x.test) > 0L) {
    x.test <- matrix(x.test, ncol = length(x.test))
  }
  if (is.numeric(x.test) && NCOL(x.test) == 0L) {
    return(NULL)
  }
  testFactorLevels <- NULL
  if (is.data.frame(x.test)) {
    # captured before any re-expansion below: on the indicators route
    # (factorLevels NULL, e.g. bart()'s x/y interface, which stores no
    # level table) a test factor with different levels re-expands to a
    # different column count than training's, and the mismatch is
    # otherwise reported as a bare column-count error naming neither the
    # factor nor its levels
    isFactorCol <- vapply(x.test, is.factor, FALSE)
    if (any(isFactorCol)) {
      testFactorLevels <- lapply(x.test[isFactorCol], levels)
    }
    isSparseColumn <- vapply(x.test, isSparseDataFrameColumn, FALSE)
    # the term replay: names an absent variable up front and runs
    # model.frame over the given term labels
    replayTerms <- function(data, labels) {
      testFormula <- as.formula(paste("~", paste(labels, collapse = " + ")))
      # model.frame resolves an absent term in the enclosing scope, so a
      # predictor missing from newdata that shares a name with a base object
      # (e.g. 'c') silently binds to it and fails with an opaque
      # "invalid type (builtin)"; name the missing variables up front instead
      neededVars <- all.vars(testFormula)
      missingVars <- neededVars[neededVars %not_in% names(data)]
      if (length(missingVars) > 0L) {
        stop(
          "'test' data is missing ",
          if (length(missingVars) > 1L) "variables" else "variable",
          " required by the model: '",
          toString(missingVars),
          "'"
        )
      }
      model.frame(
        formula = testFormula,
        data = data,
        na.action = stats::na.pass
      )
    }
    if (any(isSparseColumn)) {
      # sparse columns ride to the engine unexpanded, coded over the training
      # level table; the resulting container is preserved below. The model
      # frame replay takes no S4 columns, so they are lifted out, the dense
      # remainder replayed, and the columns the model uses re-attached.
      if (is.null(factorLevels)) {
        stop(
          "sparse test predictor columns require a categorical training ",
          "design; supply 'x' through the x/y interface"
        )
      }
      if (!is.null(termLabels)) {
        bareLabels <- sub("^`(.*)`$", "\\1", termLabels)
        sparseNames <- intersect(bareLabels, names(x.test)[isSparseColumn])
        denseLabels <- termLabels[bareLabels %not_in% sparseNames]
        sparseColumns <- x.test[sparseNames]
        dense <- x.test[!isSparseColumn]
        if (length(denseLabels) > 0L) {
          x.test <- replayTerms(dense, denseLabels)
          for (name in sparseNames) {
            x.test[[name]] <- sparseColumns[[name]]
          }
        } else {
          x.test <- sparseColumns
        }
      }
      x.test <- mapFactorColumnsToTrainingLevels(
        x.test,
        predictorNames,
        factorLevels
      )
      x.test <- makeCategoricalModelMatrix(x.test)
    } else {
      if (!is.null(termLabels)) {
        x.test <- replayTerms(x.test, termLabels)
      }
      if (!is.null(factorLevels)) {
        # trained with factors unexpanded: code against the training levels
        x.test <- mapFactorColumnsToTrainingLevels(
          x.test,
          predictorNames,
          factorLevels
        )
        x.test <- makeCategoricalModelMatrix(x.test)
      } else {
        if (is.list(drop) && length(drop) == length(x.test)) {
          refuseWiderTestColumns(x.test, drop)
        }
        x.test <- makeModelMatrixFromDataFrame(
          x.test,
          if (!is.null(drop)) drop else TRUE
        )
      }
    }
  }
  # a bare dgCMatrix test set takes the same resident path a mixed-container
  # test set's sparse columns already do (below), symmetric with a bare
  # dgCMatrix train set: wrap it as an all-sparse mixed container rather than
  # densifying. A bare numeric sparse matrix carries no factor levels, so it
  # cannot supply a categorical training column's values - refuse informatively
  # rather than let the bridge reject a malformed container.
  # any other sparse Matrix class takes the same path as its dgCMatrix
  x.test <- asDgCMatrix(x.test)
  if (inherits(x.test, "dgCMatrix")) {
    if (!is.null(factorLevels)) {
      stop(
        "a sparse matrix test set cannot supply values for a ",
        "categorical training column; supply 'test' as a dense matrix or ",
        "data frame instead"
      )
    }
    x.test <- wrapSparseTestMatrix(x.test)
  }
  # a container assembled elsewhere carries its own level order, so its factor
  # columns - CSC-backed and dense-backed alike - are re-coded against the
  # training tables before anything reads their codes. One branch at this
  # funnel aligns every entrance (creation, setTestPredictor, predict,
  # getTrees) by construction; a container this call built from a data frame
  # was already coded against those tables, so it passes through untouched.
  if (inherits(x.test, "dbartsMixedMatrix") && !is.null(factorLevels)) {
    x.test <- alignContainerFactorLevels(x.test, predictorNames, factorLevels)
  }
  # a sparse column's declared reference level means one thing to this
  # function's own densification (as.matrix.dbartsMixedMatrix, gated only
  # on is.na) and another to the engine (referenceCodeOf ignored for a
  # non-categorical column); refuse the mismatch here rather than let the two
  # disagree silently
  if (inherits(x.test, "dbartsMixedMatrix")) {
    refuseSparseTestReferenceAgainstTrainTypes(
      x.test,
      predictorNames,
      if (is.null(varTypes)) {
        rep.int(ORDINAL_VARIABLE, numPredictors)
      } else {
        varTypes
      }
    )
  }
  # a sparse-backed container stays resident (the engine codes it against the
  # training cuts); everything else densifies
  xTestIsSparseContainer <-
    inherits(x.test, "dbartsMixedMatrix") && predictorSourceIsSparse(x.test)
  if (!is.matrix(x.test) && !xTestIsSparseContainer) {
    x.test <- as.matrix(x.test)
  }

  if (!xTestIsSparseContainer) {
    if (!is.numeric(x.test)) {
      stop("test matrix must be numeric")
    }

    if (is.integer(x.test)) {
      # storage.mode<- keeps the dimnames that matrix() would drop
      storage.mode(x.test) <- "double"
    }
  }

  if (!identical(NCOL(x.test), numPredictors)) {
    if (is.null(factorLevels) && !is.null(testFactorLevels)) {
      mismatched <- Filter(
        function(nm) {
          !identical(
            length(resolveTermColumns(nm, predictorNames, termLabels)),
            length(resolveTermColumns(nm, colnames(x.test), termLabels))
          )
        },
        names(testFactorLevels)
      )
      if (length(mismatched) > 0L) {
        stop(
          "'test' factor '",
          mismatched[1L],
          "' does not match training's indicator columns ('test' levels: ",
          toString(testFactorLevels[[mismatched[1L]]]),
          "); use bart() or dbarts(), which track levels across predict ",
          "by default"
        )
      }
    }
    stop("number of columns in 'test' must be equal to that of 'x'")
  }
  if (numPredictors > 1) {
    xIsNamed <- !is.null(predictorNames)
    testIsNamed <- !is.null(colnames(x.test))

    columnIndices <- seq.int(numPredictors)
    if (xIsNamed && !testIsNamed) {
      # named fit, positional (unnamed) test: swapped columns would match
      # silently and return badly wrong numbers, so spell out the mapping the
      # positional match assumes rather than only noting that it happened
      shown <- min(numPredictors, 3L)
      mapping <- paste0(
        "column ",
        seq_len(shown),
        " = '",
        predictorNames[seq_len(shown)],
        "'",
        collapse = ", "
      )
      # shares dbartsPositionalArgsWarning with $setResponse's positional
      # warning and the massign position-only site: all three report the
      # same condition, columns/arguments matched by position rather than
      # by name
      warning(warningCondition(
        paste0(
          "'test' is unnamed but 'x' had named predictors, matched to 'x' by ",
          "position (",
          mapping,
          if (numPredictors > shown) ", ..." else "",
          "); supply 'test' with column names to match by name instead"
        ),
        class = c("dbartsPositionalArgsWarning", "dbartsWarning")
      ))
    } else if (
      (!xIsNamed && testIsNamed) ||
        length(unique(predictorNames)) != length(predictorNames)
    ) {
      warning(warningCondition(
        "'x' and 'test' are not both named; columns of 'test' will be matched by position",
        class = c("dbartsPositionalArgsWarning", "dbartsWarning")
      ))
    } else if (xIsNamed && testIsNamed) {
      matchIndices <- match(predictorNames, colnames(x.test))
      if (any(is.na(matchIndices))) {
        stop(
          "column names of 'test' do not match those of 'x': '",
          toString(predictorNames[is.na(matchIndices)]),
          "' present in 'x' but not in 'test' (whose columns are '",
          toString(colnames(x.test)),
          "')"
        )
      } else {
        columnIndices <- matchIndices
      }
    }

    if (xTestIsSparseContainer) {
      # reorder columns without densifying: the map re-points each predictor
      # position at its source, the names follow
      x.test$map <- x.test$map[columnIndices]
      if (!is.null(x.test$columnNames)) {
        x.test$columnNames <- x.test$columnNames[columnIndices]
      }
      if (xIsNamed) x.test$columnNames <- predictorNames
    } else {
      x.test <- x.test[, columnIndices, drop = FALSE]
      if (xIsNamed) colnames(x.test) <- predictorNames
    }
  }

  if (refuseMissing) {
    refuseTestMissingness(x.test, x.train)
  }

  x.test
}

findTermInFormulaData <- function(formula, data, term) {
  dataIsMissing <- missing(data)
  matchedCall <- match.call()

  if (is.numeric(matchedCall$term)) {
    return(term)
  }

  if (!dataIsMissing) {
    if (is.symbol(matchedCall$term)) {
      if (any(names(data) == as.character(matchedCall$term))) {
        return(data[[as.character(matchedCall$term)]])
      }
    } else if (is.language(matchedCall$term)) {
      tryResult <- with(
        data,
        tryCatch(eval(matchedCall$term), error = function(e) e)
      )
      if (!inherits(tryResult, "error")) return(tryResult)
    }
  }
  if (is.symbol(matchedCall$term)) {
    if (any(ls(environment(formula)) == as.character(matchedCall$term))) {
      return(get(as.character(matchedCall$term), envir = environment(formula)))
    }
    tryResult <- tryCatch(
      get(as.character(matchedCall$term)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error") && !is.null(tryResult)) return(tryResult)
  } else if (is.language(matchedCall$term)) {
    tryResult <- tryCatch(
      eval(matchedCall$term, environment(formula)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error")) {
      return(tryResult)
    }
    tryResult <- tryCatch(eval(matchedCall$term), error = function(e) e)
    if (!inherits(tryResult, "error")) return(tryResult)
  }

  NULL
}

## A block of code rather than a function: evaluating a function this way in
## the caller's frame triggers an R CMD check warning.
getTestOffset <- quote({
  if (is.numeric(matchedCall$offset.test)) {
    return(list(offset.test = offset.test, testUsesRegularOffset = FALSE))
  }
  if (is.null(matchedCall$offset.test)) {
    return(list(offset.test = NULL, testUsesRegularOffset = FALSE))
  }

  if (is.symbol(matchedCall$offset.test)) {
    testOffsetName <- as.character(matchedCall$offset.test)

    if (identical(testOffsetName, "offset") && !is.null(offset)) {
      return(list(
        offset.test = if (offsetGivenAsScalar == TRUE) offset[1] else offset,
        testUsesRegularOffset = TRUE
      ))
    }

    if (is.formula(formula)) {
      if (!dataIsMissing && any(names(data) == testOffsetName)) {
        return(list(
          offset.test = data[[testOffsetName]],
          testUsesRegularOffset = FALSE
        ))
      }
      if (any(ls(environment(formula)) == testOffsetName)) {
        return(list(
          offset.test = get(testOffsetName, environment(formula)),
          testUsesRegularOffset = FALSE
        ))
      }
    }
    tryResult <- tryCatch(get(testOffsetName), error = function(e) e)
    if (!inherits(tryResult, "error") && !is.null(tryResult)) {
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }

    stop("cannot find test offset '", testOffsetName, "'")
  } else if (is.language(matchedCall$offset.test)) {
    ## offset.test could have been something like (offset + 0.5), or (offset + variable)
    baseOffset <- if (is.null(offset)) {
      NA_real_
    } else {
      if (offsetGivenAsScalar == TRUE) offset[1] else offset
    }

    if (identical(matchedCall$offset.test, quote(offset))) {
      return(list(offset.test = baseOffset, testUsesRegularOffset = TRUE))
    }

    testOffset <- subTermInLanguage(
      matchedCall$offset.test,
      quote(offset),
      baseOffset
    )

    if (is.formula(formula)) {
      if (!dataIsMissing) {
        tryResult <- with(
          data,
          tryCatch(eval(testOffset), error = function(e) e)
        )
        if (!inherits(tryResult, "error")) {
          return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
        }
      }
      tryResult <- tryCatch(
        eval(testOffset, environment(formula)),
        error = function(e) e
      )
      if (!inherits(tryResult, "error")) {
        return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
    }
    tryResult <- tryCatch(
      eval(testOffset, parent.frame(3L)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error")) {
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }
  }

  stop("cannot construct test offset")
})

# Classify a raw response for family routing. A factor/ordered/logical/
# character response declares a classification model; numeric passes through.
# The level count drives the 2 -> probit vs 3+ -> multinomial (bart2) / refusal
# (single-forest) split downstream. Shared by dbartsData's response coding and
# bart2's family = "auto" peek.
classifyResponse <- function(y) {
  if (is.factor(y)) {
    list(
      type = if (is.ordered(y)) "ordered factor" else "factor",
      n.levels = nlevels(y)
    )
  } else if (is.logical(y) && !is.matrix(y)) {
    list(type = "logical", n.levels = 2L)
  } else if (is.character(y) && !is.matrix(y)) {
    list(type = "character", n.levels = length(unique(y[!is.na(y)])))
  } else {
    list(type = "numeric", n.levels = NA_integer_)
  }
}

# Code a raw response to the doubles the engine reads and report its original
# type. A factor (or a character coerced with factor(), or a logical) becomes
# 0-based codes; a numeric response passes through as.double. The original
# levels are returned alongside so an ordinal fit can round-trip them; they
# are NULL for a numeric response (a numeric ordinal derives sort(unique(y))
# itself).
# codeResponse flattens any matrix response with as.double(), column-major,
# to length ncol(y) * nrow(y); called on a Surv object or an n x K matrix
# without this guard first, it silently produces a length mismatch against
# 'x' rather than naming the response shape dbartsData() cannot ingest.
refuseMultiColumnResponse <- function(y) {
  if (inherits(y, "Surv")) {
    stop(
      "'y' is a survival response (Surv); dbartsData() takes a single-",
      "column response - fit through dbarts()/bart() with family = ",
      "\"aft\" or \"hazard\", which extract time and status first"
    )
  }
  if (is.matrix(y) && ncol(y) > 1L) {
    # classed so that dbarts() and xbart() can restate the refusal in terms
    # of the entry the caller used, which dbartsData() cannot know
    text <- if (ncol(y) == 2L) {
      paste0(
        "'y' is an n x 2 matrix; dbartsData() takes a single-column ",
        "response - a (time, status) pair goes to dbarts()/bart() with ",
        "family = \"aft\"/\"hazard\", per-category counts to ",
        "dbartsData(counts = )"
      )
    } else {
      paste0(
        "'y' is an n x ",
        ncol(y),
        " matrix; dbartsData() takes a single-column response - pass per-",
        "category counts as dbartsData(counts = ) and fit with family = ",
        "\"multinomial\""
      )
    }
    stop(structure(
      class = c("dbartsMatrixResponseError", "error", "condition"),
      list(message = text, call = NULL, y = y)
    ))
  }
  invisible(NULL)
}

# TRUE for an n x K (K >= 3) numeric matrix whose complete rows hold only
# non-negative whole numbers: what family = "auto" reads as multinomial
# counts. A row with a missing cell is a missing response, left to
# na.action exactly as the explicit multinomial path leaves it.
isAutoCountMatrix <- function(y) {
  if (!is.matrix(y) || !is.numeric(y) || inherits(y, "Surv") || ncol(y) < 3L) {
    return(FALSE)
  }
  observed <- y[rowSums(is.na(y)) == 0L, , drop = FALSE]
  all(observed >= 0) && all(observed == round(observed))
}

# Restates a matrix-response refusal from dbartsData() for the entry the
# caller used: what to write there, under the family they wrote. 'entry' is
# "bart()/dbarts()" or "xbart()"; 'family' is the family token as requested.
restateMatrixResponseError <- function(e, entry, family) {
  y <- e$y
  K <- ncol(y)
  readings <- paste0(
    "per-category counts need family = \"multinomial\"; a (time, status) ",
    "pair needs a survival::Surv response (which \"auto\" reads as ",
    "survival) or family = \"aft\" / \"hazard\""
  )
  text <- if (identical(entry, "xbart()")) {
    paste0(
      "xbart() takes a single-column response but 'y' is an n x ",
      K,
      " matrix, and it cross-validates neither multinomial counts nor a ",
      "survival response; fit those with bart() instead - ",
      readings
    )
  } else if (family != "auto") {
    paste0(
      "family = \"",
      family,
      "\" takes a single-column response but 'y' is an n x ",
      K,
      " matrix - ",
      readings
    )
  } else if (K == 2L) {
    paste0(
      "'y' is an n x 2 matrix, which family = \"auto\" leaves ambiguous - ",
      readings
    )
  } else if (isAutoCountMatrix(y)) {
    paste0(
      "dbarts() fits an n x ",
      K,
      " count matrix through the matrix interface only: ",
      "dbarts(x, counts, family = \"multinomial\"); bart() takes the formula ",
      "form, bart(cbind(c1, c2, c3) ~ x, data)"
    )
  } else {
    paste0(
      "'y' is an n x ",
      K,
      " matrix with negative, fractional or non-numeric entries, which ",
      "family = \"auto\" does not read as multinomial counts; per-category ",
      "counts must be non-negative whole numbers (NA marks a missing row)"
    )
  }
  stop(text, call. = FALSE)
}

withMatrixResponseRestated <- function(entry, family, expr) {
  tryCatch(
    expr,
    dbartsMatrixResponseError = function(e) {
      restateMatrixResponseError(e, entry, family)
    }
  )
}

codeResponse <- function(y) {
  info <- classifyResponse(y)
  coded <- if (info$type == "numeric" || info$type == "logical") {
    as.double(y)
  } else {
    as.double(as.integer(if (is.character(y)) factor(y) else y) - 1L)
  }
  levels <- if (is.factor(y)) {
    levels(y)
  } else if (info$type == "character") {
    levels(factor(y))
  } else {
    NULL
  }
  list(y = coded, type = info$type, n.levels = info$n.levels, levels = levels)
}

# Resolve family for a categorical (factor/logical/character) response ahead
# of a single-forest fit: a 2-level response is a binary classification
# (family = "auto" fits probit); 3+ levels is multinomial, which none of
# dbarts()/xbart() implement (only bart2(family = "multinomial")
# does). A numeric response is returned unchanged - a caller that also
# resolves the 0/1-vs-continuous ambiguity for numeric responses does so
# itself afterward (dbarts()/xbart() do).
#
# `caller` names the entry point for the 2-level conflict message and the
# non-split multinomial message ("CALLER does not fit a K-level ..."), used
# as-is by xbart(): it is reached only directly, so the
# auto/explicit distinction adds nothing at K >= 3 (every family choice is
# equally invalid). dbarts() passes splitMultinomialMessage = TRUE instead,
# because it is also reached anonymously through bart() (which never sets an
# explicit family): its auto-branch message cannot name a single caller and
# instead lists every single-forest entry point (passed via `caller`), while
# its explicit-family branch echoes the conflicting family like the 2-level
# message does.
resolveClassificationFamily <- function(
  data,
  family,
  caller,
  incompatibleFamilies,
  splitMultinomialMessage = FALSE,
  allowOrdinal = FALSE
) {
  responseType <- data@response.type
  K <- data@response.n.levels
  # ordinal (cumulative probit): only
  # dbarts()/bart() pass allowOrdinal = TRUE. family = "ordinal" is the
  # explicit primitive and forces the model on any response (numeric levels are
  # sort(unique(y)), resolved later); family = "auto" auto-dispatches an ORDERED
  # factor to it, announced. The other single-forest entries leave allowOrdinal
  # FALSE and fall through to their K >= 3 refusals below. is.ordered() is the
  # disjoint key: an unordered K >= 3 factor stays multinomial.
  if (allowOrdinal) {
    if (identical(family, "ordinal")) {
      return("ordinal")
    }
    # a 2-level ordered factor is binary (probit); only a 3+-level ordered
    # factor is a genuine ordinal scale worth auto-dispatching
    if (family == "auto" && responseType == "ordered factor" && K >= 3L) {
      return("ordinal")
    }
  }
  if (responseType == "numeric") {
    return(family)
  }
  if (K >= 3L) {
    # a 3+-level ORDERED factor is ordinal (reached here only from an entry that
    # cannot fit it - xbart; dbarts/bart route ordinal above); every
    # other 3+-level factor/character is unordered multinomial
    isOrdered <- identical(responseType, "ordered factor")
    model <- if (isOrdered) "ordinal" else "multinomial"
    suggestion <- if (isOrdered) {
      "bart(x, y, family = \"ordinal\")"
    } else {
      "bart(x, y, family = \"multinomial\")"
    }
    if (!splitMultinomialMessage) {
      stop(
        caller,
        " does not fit a ",
        K,
        "-level ",
        responseType,
        " response; ",
        model,
        " classification requires ",
        suggestion
      )
    }
    if (family == "auto") {
      stop(
        "a ",
        K,
        "-level ",
        responseType,
        " response is ",
        model,
        "; fit it with ",
        suggestion,
        " - ",
        caller,
        " fit only binary and continuous responses"
      )
    }
    stop(
      "family \"",
      family,
      "\" cannot fit a ",
      K,
      "-level ",
      responseType,
      " response; a 3+-level ",
      responseType,
      " is ",
      model,
      " (",
      suggestion,
      ")"
    )
  }
  if (family == "auto") {
    family <- "probit"
  } else if (family %in% incompatibleFamilies) {
    stop(
      "family \"",
      family,
      "\" cannot fit a ",
      responseType,
      " response; a 2-level factor is a binary classification ",
      "(family = \"auto\" fits probit)"
    )
  }
  family
}

# Resolve the ordered category structure for an ordinal (cumulative-probit) fit
# from a dbartsData whose family has resolved to "ordinal". Returns the
# ONE-based category codes the engine reads
# (its y_ holds 1..K), the count K, and the ordered level labels for the
# round-trip. A factor/character response takes its stored level order - an
# UNORDERED one with an informational note, since the factor default
# (alphabetical) is rarely the intended scale; a numeric/integer/logical
# response uses sort(unique(y)) as the ordered levels. Called once, where the
# recoding happens (dbarts()), so the note fires once.
resolveOrdinalResponse <- function(data) {
  labels <- data@response.levels
  if (is.null(labels)) {
    # numeric/integer/logical: the distinct sorted values are the ordered
    # levels, so a continuous response would silently become one category per
    # distinct value; only a whole-number response names a plausible set of
    # ordered categories
    if (any(data@y != round(data@y))) {
      stop(
        "family \"ordinal\" requires a factor or an integer-valued ",
        "response, not continuous values"
      )
    }
    levels <- sort(unique(data@y))
    codes <- match(data@y, levels)
    labels <- as.character(levels)
  } else {
    # data@y holds 0-based factor codes (codeResponse); the engine wants 1..K
    codes <- as.integer(data@y) + 1L
    if (data@response.type != "ordered factor") {
      message(
        "family = \"ordinal\": the ",
        data@response.type,
        " response is unordered; its category order is taken from the level ",
        "order (",
        paste(labels, collapse = " < "),
        ")"
      )
    }
  }
  K <- length(labels)
  if (K < 2L) {
    stop("family = \"ordinal\" requires a response with at least 2 categories")
  }
  list(y = as.double(codes), K = as.integer(K), levels = labels)
}

# Resolve the negative-binomial dispersion argument to the length-1 real the C
# bridge reads off the control's bartcore.dispersion attribute. NA (the
# default) estimates r on the capped integer grid, encoded as a non-positive
# spec; a supplied value FIXES r as a single positive integer - the exact
# integer envelope v1 ships, real dispersion not yet supported.
resolveDispersion <- function(dispersion) {
  if (length(dispersion) != 1L) {
    stop("'dispersion' must be a single value")
  }
  if (is.na(dispersion)) {
    return(-1) # a non-positive spec estimates r on the grid
  }
  if (!is.numeric(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a positive number")
  }
  if (dispersion != round(dispersion)) {
    stop(
      "family \"nbinom\" fits an integer dispersion; real dispersion is not ",
      "yet supported"
    )
  }
  as.double(dispersion)
}

# Validate and subset a user-supplied weights vector for the x/y interfaces
# (the sparse dgCMatrix branch and the numeric/data.frame/factor branch carry
# byte-identical validation): a NULL passes through, otherwise the vector is
# numeric-checked, length-1 recycled to the observation count, length-validated
# against 'y', and finally restricted to 'subset'. The missing(weights) guard
# that produces the NULL stays inline at each call site, since missing() must
# reference the dbartsData formal.
validateXYWeights <- function(weights, initialNumObservations, subset) {
  if (is.null(weights)) {
    return(NULL)
  }
  if (!is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  weights <- as.double(weights)
  if (length(weights) == 1L) {
    weights <- rep_len(weights, initialNumObservations)
  }
  if (length(weights) != initialNumObservations) {
    stop("'weights' must have the same length as 'y'")
  }
  weights[subset]
}

# Validate and subset a user-supplied list of per-forest amplitude bases, the
# matrices a multi-forest fit combines its forests through. A NULL passes
# through; otherwise every non-null element
# is coerced to a numeric matrix. 'subsetRows' (resolveFormulaBasisSubset's
# result, R/model.R), when given, resolves each element the way a forests =
# declaration's basis already does - alignForestBasisToSubset, the same
# function - so the formula path and the x/y path converge on one rule: a
# basis at the pre-'subset' data's row count is restricted to the rows a
# formula's 'subset' kept, and a basis matching the kept-row count but not the
# full data's is refused by name rather than guessed at. With no 'subsetRows'
# (the x/y interface, or a formula with no 'subset'), the row-index
# contract applies directly: 'initialNumObservations' is the row count to
# match and 'subset', if given, is applied as-is. 'argument' names the surface
# the value came from, so a caller who wrote a forest's basis is refused in
# those terms rather than in the data object's, and 'rows' names what the row
# count is checked against - the response at fit time, the predicted rows at
# predict time.
validateForestBases <- function(
  bases,
  initialNumObservations,
  subset = NULL,
  argument = "bases",
  rows = "'y'",
  subsetRows = NULL
) {
  if (is.null(bases)) {
    return(NULL)
  }
  if (!is.list(bases)) {
    stop("'", argument, "' must be a list of per-forest basis matrices")
  }
  lapply(seq_along(bases), function(i) {
    basis <- bases[[i]]
    if (is.null(basis)) {
      return(NULL)
    }
    if (!is.numeric(basis) && !is.logical(basis)) {
      stop("'", argument, "' must be numeric or logical")
    }
    basis <- as.matrix(basis)
    storage.mode(basis) <- "double"
    if (!is.null(subsetRows)) {
      basis <- alignForestBasisToSubset(basis, i, subsetRows)
    } else {
      if (nrow(basis) != initialNumObservations) {
        stop("'", argument, "' must have the same length as ", rows)
      }
      if (!is.null(subset)) {
        basis <- basis[subset, , drop = FALSE]
      }
    }
    if (anyNA(basis) || !all(is.finite(basis))) {
      stop("'", argument, "' values must all be finite")
    }
    basis
  })
}

# Validate and subset a user-supplied offset vector for the x/y interfaces
# (shared by both x/y branches): a NULL passes through unchanged, otherwise
# the vector is numeric-checked, length-1 recycled (which
# sets offsetGivenAsScalar = TRUE, a longer vector FALSE), length-validated
# against 'y', and restricted to 'subset'. offsetGivenAsScalar is threaded in
# and back out so a NULL offset leaves it untouched. The offsetIsMissing guard
# that produces the NULL stays inline at each call site.
validateXYOffset <- function(
  offset,
  initialNumObservations,
  subset,
  offsetGivenAsScalar
) {
  if (!is.null(offset)) {
    if (!is.numeric(offset)) {
      stop("'offset' must be numeric")
    }
    if (length(offset) == 1L) {
      offset <- rep_len(offset, initialNumObservations)
      offsetGivenAsScalar <- TRUE
    } else {
      offsetGivenAsScalar <- FALSE
    }
    if (length(offset) != initialNumObservations) {
      stop("'offset' must have the same length as 'y'")
    }
    offset <- offset[subset]
  }
  list(offset = offset, offsetGivenAsScalar = offsetGivenAsScalar)
}

# Coerce a raw multinomial response to the n x K count matrix the softmax
# engine's response is. Two shapes arrive: an
# n x K matrix (or data frame) of non-negative integer counts, whose trials are
# its row sums, and a length-n vector of labels - a factor, character, logical
# or integer code - which is the single-trial special case, one-hot expanded
# with every trial 1. The category labels ride the result's COLUMN NAMES, the
# carrier that survives both serialization and the engine's re-creation; the
# engine reads neither.
# A missing label one-hot-expands to an all-NA row rather than erroring: a
# single-trial spike has no cell to put a missing category in, so the row
# itself is what is missing, exactly as a count-matrix row is missing when
# any of its cells is NA (validateMultinomialCounts). na.action reads that
# row downstream, as it does for every other family's response.
resolveMultinomialCounts <- function(y) {
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
  if (is.matrix(y)) {
    if (!is.numeric(y)) {
      stop(
        "a multinomial count-matrix response must be numeric, not ",
        typeof(y)
      )
    }
    return(y)
  }
  if (is.character(y) || is.logical(y)) {
    y <- factor(y)
  }
  levels <- if (is.factor(y)) {
    levels(y)
  } else {
    if (!is.numeric(y)) {
      stop(
        "a multinomial response must be a factor, a character vector, or an ",
        "n x K count matrix"
      )
    }
    observed <- y[!is.na(y)]
    if (any(observed != round(observed)) || any(observed < 0)) {
      stop("multinomial category codes must be non-negative whole numbers")
    }
    as.character(seq.int(0L, if (length(observed) > 0L) max(observed) else 0L))
  }
  codes <- if (is.factor(y)) as.integer(y) else as.integer(y) + 1L
  K <- length(levels)
  counts <- matrix(0L, length(codes), K, dimnames = list(NULL, levels))
  complete <- !is.na(codes)
  counts[cbind(seq_along(codes)[complete], codes[complete])] <- 1L
  counts[!complete, ] <- NA_integer_
  counts
}

# Validate and subset a multinomial count response for dbartsData: an n x K
# matrix of non-negative whole numbers with at least two categories. A row with
# no trial is accepted, entering no likelihood, and this stays silent about it:
# the warning is warnZeroTrials', raised where a sampler takes the matrix. The
# engine re-derives the trials and re-checks every invariant; this is the R
# layer's own (safe over fast) refusal, and the one that names the argument
# the caller wrote. 'subset' already reflects na.action's own row selection
# when allowMissing is TRUE (dbartsData's own ingestion): a row dropped there
# (every na.action but na.pass) is gone after the subset below, and na.pass
# keeps it, all-NA, for the generic missing-response check downstream to
# name. allowMissing is FALSE for every
# other caller ($setCounts mutates a live sampler's response outright, with
# no na.action of its own to defer to), where any NA is refused here by name.
# The whole-number check runs on the ORIGINAL doubles, before asCountMatrix's
# integer coercion would silently floor a fraction into passing.
validateMultinomialCounts <- function(
  counts,
  initialNumObservations,
  subset,
  allowMissing = FALSE
) {
  if (is.null(counts)) {
    return(NULL)
  }
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("'counts' must be a numeric matrix")
  }
  if (nrow(counts) != initialNumObservations) {
    stop("'counts' must have the same number of rows as 'y' has elements")
  }
  if (ncol(counts) < 2L) {
    stop("'counts' must have at least two categories")
  }
  missingRow <- rowSums(is.na(counts)) > 0L
  if (!allowMissing && any(missingRow)) {
    stop("'counts' cannot be NA")
  }
  observed <- counts[!missingRow, , drop = FALSE]
  if (any(observed < 0)) {
    stop("'counts' must all be non-negative")
  }
  if (any(observed != round(observed))) {
    stop("'counts' must all be whole numbers")
  }
  asCountMatrix(counts)[subset, , drop = FALSE]
}

# Evaluates expr, a re-creation of a sampler that already exists, with the
# zero-trial warning held off: re-creating introduces no rows, so it neither
# warns nor spends the key.
withoutZeroTrialsWarning <- function(expr) {
  spent <- onceWarnState[["multinomialZeroTrials"]]
  onceWarnState[["multinomialZeroTrials"]] <- TRUE
  on.exit(onceWarnState[["multinomialZeroTrials"]] <- spent)
  expr
}

# Warns, once per session, that a count matrix a sampler has just taken has
# rows with no trial. Called only after the creation or $setCounts .Call has
# returned, so a refused call never spends the key; NULL (a sampler that
# carries no counts) and a matrix without such a row are silent.
warnZeroTrials <- function(counts) {
  if (is.null(counts)) {
    return(invisible(NULL))
  }
  numEmpty <- sum(rowSums(counts) == 0L, na.rm = TRUE)
  if (numEmpty > 0L) {
    warnOnce(
      "multinomialZeroTrials",
      warningCondition(
        sprintf(
          paste0(
            "multinomial count rows with zero trials (%d of %d) contribute ",
            "nothing to the likelihood and still receive fitted probabilities"
          ),
          numEmpty,
          nrow(counts)
        ),
        class = c("dbartsZeroTrialsWarning", "dbartsWarning")
      )
    )
  }
  invisible(NULL)
}

# A matrix 'offset'/'offset.test' is only ever meaningful on a counts-
# carrying (multinomial) data object, where it is the n x K category shift;
# reached whenever one arrives without 'counts' to pair it with. 'argument'
# is the spelling the caller used, which is the only one they can act on.
refuseMatrixOffset <- function(offset, argument) {
  stop(
    "'",
    argument,
    "' must be a numeric vector of length n or a single number; a ",
    nrow(offset),
    " x ",
    ncol(offset),
    " matrix was supplied, which only family = \"multinomial\" accepts"
  )
}

# A flat (per-observation) offset is the softmax's own null direction - a
# constant added to every category of a row leaves every reported
# probability unchanged - so it carries no information a multinomial fit can
# use; the meaningful shift is the n x K category offset instead.
refuseFlatOffsetOnMultinomial <- function(offset, argument = "offset") {
  stop(
    "family = \"multinomial\" requires an n x K matrix \"",
    argument,
    "\", one column per category; a length-",
    length(offset),
    " vector was supplied, and a common per-observation shift is the ",
    "softmax's own null direction - it cancels"
  )
}

# Validate (and, for the training-row twin, subset) a category offset ARRIVING
# ON A DATA OBJECT: the n x K matrix added to the raw per-category fits before
# the softmax. 'rows' is NULL for the test twin, whose row count belongs to the
# test store and is pinned by the validity method instead. Distinct from
# validateCategoryOffset, which checks a matrix arriving at a bartcore HANDLE
# against that handle's own already-fixed n and K.
validateDataCategoryOffset <- function(
  offset,
  rows,
  subset,
  argument,
  categories = NULL
) {
  if (is.null(offset)) {
    return(NULL)
  }
  if (is.data.frame(offset)) {
    offset <- as.matrix(offset)
  }
  if (!is.matrix(offset) || !is.numeric(offset)) {
    stop("'", argument, "' must be a numeric matrix")
  }
  # checked here rather than left to the class's validity method, which can
  # only name the slot the value lands in and not the argument it arrived as
  if (!is.null(categories) && ncol(offset) != categories) {
    stop(
      "'",
      argument,
      "' must have one column per category (K = ",
      categories,
      "); a ",
      nrow(offset),
      " x ",
      ncol(offset),
      " matrix was supplied"
    )
  }
  if (anyNA(offset) || !all(is.finite(offset))) {
    stop("'", argument, "' values must all be finite")
  }
  storage.mode(offset) <- "double"
  if (is.null(rows)) {
    return(offset)
  }
  if (nrow(offset) != rows) {
    stop(
      "'",
      argument,
      "' must have the same number of rows as 'y' has elements"
    )
  }
  offset[subset, , drop = FALSE]
}

dbartsData <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  factors = c("categorical", "indicators"),
  na.action = dbarts::na.keepPredictors,
  bases = NULL,
  counts = NULL
) {
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  offsetIsMissing <- missing(offset)
  testOffsetIsMissing <- missing(offset.test)
  basesIsMissing <- missing(bases)
  countsIsMissing <- missing(counts)
  matchedCall <- match.call()
  # a matrix-shaped 'offset'/'offset.test' declares a multinomial category
  # shift (one column per category), never a flat per-row one; the matrix
  # interface branches below set these aside as they resolve the ordinary
  # flat offset, and validateDataCategoryOffset installs them once 'counts'
  # is known
  categoryOffset <- NULL
  categoryTestOffset <- NULL

  # "indicators" dummy-expands factor columns as always; "categorical" keeps
  # them as single columns split by category subset, which only the bartcore
  # engine runs
  factors <- match.arg(factors)
  makeModelMatrix <- if (factors == "categorical") {
    makeCategoricalModelMatrix
  } else {
    makeModelMatrixFromDataFrame
  }
  # the rows an incomplete case costs, and the record of them the training
  # fits pad through; NULL until an na.action drops something
  naOmitted <- NULL
  # the kept rows' names, captured from the raw inputs as they are subset
  # (dataRowNames); NULL for unnamed rows
  trainRowNames <- NULL

  # a Surv formula response's raw event/censoring time and 0/1 status,
  # parked as attributes on the returned object (below) rather than a slot -
  # dbartsData() has no family vocabulary to log-transform or person-period
  # expand it itself, so dbarts() decodes them once the caller's requested
  # family is known. NULL for every response that is not Surv.
  survivalStatus <- NULL
  survivalTime <- NULL
  survivalTimeOmitted <- NULL

  offsetGivenAsScalar <- NA
  testUsesRegularOffset <- NA
  # the response's original type, recorded on the result so the fitters can
  # route family = "auto" and reject a categorical response an unsupported
  # family cannot fit; each y-producing branch below refreshes it
  responseInfo <- list(type = "numeric", n.levels = NA_integer_, levels = NULL)

  if (missing(formula)) {
    stop("first argument to dbartsData - 'formula'/'x.train' - must be present")
  }

  if (inherits(formula, "dbartsData")) {
    if (
      !dataIsMissing ||
        !testIsMissing ||
        !offsetIsMissing ||
        !testOffsetIsMissing ||
        !basesIsMissing ||
        !countsIsMissing
    ) {
      warning(warningCondition(
        "if data supplied as dbartsData, remaining arguments are ignored",
        class = c("dbartsIgnoredArgWarning", "dbartsWarning")
      ))
    }
    return(formula)
  }

  if (is.formula(formula)) {
    if (
      !dataIsMissing &&
        !is.data.frame(data) &&
        !is.list(data) &&
        !is.environment(data)
    ) {
      stop(
        "for formula/data specification, data must be a data frame, list, or environment"
      )
    }

    modelFrameArgs <- c("formula", "data", "subset", "weights", "offset")

    ## extract offset prematurely, if necessary
    if (offsetIsMissing) {
      offset <- NULL
      modelFrameArgs <- c("formula", "data", "subset", "weights")
    } else {
      offsetCall <- matchedCall
      offsetCall <- offsetCall[c(
        1L,
        match(c("formula", "data", "offset"), names(offsetCall), nomatch = 0L)
      )]
      names(offsetCall)[which(names(offsetCall) == "offset")] <- "term"
      offsetCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
      offset <- eval(offsetCall, parent.frame())

      # a matrix-shaped offset declares a multinomial category shift, one
      # column per category, never a flat per-row one; set it aside before the
      # model frame is built (which would read it as a per-row term) and
      # refuse it where there is no 'counts' for it to belong to
      if (!is.null(offset) && (is.matrix(offset) || is.data.frame(offset))) {
        if (is.null(counts)) {
          refuseMatrixOffset(offset, "offset")
        }
        categoryOffset <- if (is.data.frame(offset)) {
          as.matrix(offset)
        } else {
          offset
        }
        offset <- NULL
        modelFrameArgs <- c("formula", "data", "subset", "weights")
      }

      if (!is.null(offset)) {
        offsetGivenAsScalar <- length(offset) == 1
        if (offsetGivenAsScalar) {
          modelFrameArgs <- c("formula", "data", "subset", "weights")
        }
      }
    }

    # pre-validate lengths against a known y/data length so a mismatch reads
    # as our own message rather than model.frame's "variable lengths differ
    # (found for '(weights)')"; a data.frame is the only case where the
    # eventual length is known this early without duplicating model.frame's
    # own work
    dataLength <- if (!dataIsMissing && is.data.frame(data)) {
      nrow(data)
    } else {
      NA_integer_
    }
    if (!is.na(dataLength)) {
      if (
        !is.null(offset) &&
          !isTRUE(offsetGivenAsScalar) &&
          length(offset) != dataLength
      ) {
        stop("'offset' must have the same length as 'y'")
      }
      if (!missing(weights)) {
        weightsCall <- matchedCall[c(
          1L,
          match(
            c("formula", "data", "weights"),
            names(matchedCall),
            nomatch = 0L
          )
        )]
        names(weightsCall)[names(weightsCall) == "weights"] <- "term"
        weightsCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
        weightsValue <- tryCatch(
          eval(weightsCall, parent.frame()),
          error = function(e) NULL
        )
        if (!is.null(weightsValue) && length(weightsValue) != dataLength) {
          stop("'weights' must have the same length as 'y'")
        }
      }
    }

    modelFrameCall <- matchedCall
    modelFrameCall <- modelFrameCall[c(
      1L,
      match(modelFrameArgs, names(modelFrameCall), nomatch = 0L)
    )]
    modelFrameCall$drop.unused.levels <- FALSE
    # the one site the caller's own na.action governs; every other model
    # frame this function builds is a re-read of columns these rows already
    # settled, so those keep na.pass and align to whatever this frame kept
    modelFrameCall$na.action <- na.action
    modelFrameCall[[1L]] <- quote(stats::model.frame)
    ## this allows subset to be applied to offset, even if offset was a language construct (e.g. off + 0.1)
    if (identical(offsetGivenAsScalar, FALSE)) {
      modelFrameCall$offset <- offset
    }

    # a sparseVector/dgCMatrix/sparseFactor column would die inside
    # model.frame with a bare S4 type error; a data-frame 'data' has row
    # names to re-attach a pulled-out column by (below), so its sparse
    # columns are lifted out here rather than refused. A plain list or an
    # environment has no such row identity to align by, so those keep the
    # plain refusal.
    usedSparseNames <- character(0)
    sparseColumns <- list()
    if (!dataIsMissing && is.data.frame(data)) {
      pulledOut <- pullOutSparseFormulaColumns(data)
      if (length(pulledOut$sparseColumns) > 0L) {
        denseData <- pulledOut$denseData
        sparseColumns <- pulledOut$sparseColumns
        # '.' has to be able to reach a sparse name too, so it is expanded
        # by hand against a placeholder frame naming every column of the
        # ORIGINAL 'data' - dense ones densely, sparse ones as zero-length
        # stand-ins terms() never reads the values of, only the names -
        # before 'formula' is rewritten to name only its dense terms
        # explicitly
        placeholderFrame <- denseData[0L, , drop = FALSE]
        for (sparseName in names(sparseColumns)) {
          placeholderFrame[[sparseName]] <- numeric(0)
        }
        expandedTerms <- terms(formula, data = placeholderFrame)
        denseTermLabels <- character(0)
        for (label in attr(expandedTerms, "term.labels")) {
          bareLabel <- sub("^`(.*)`$", "\\1", label)
          sparseHits <- intersect(
            all.vars(str2lang(label)),
            names(sparseColumns)
          )
          if (length(sparseHits) == 0L) {
            denseTermLabels <- c(denseTermLabels, label)
          } else if (bareLabel %in% names(sparseColumns)) {
            usedSparseNames <- c(usedSparseNames, bareLabel)
          } else {
            stop(
              "sparse predictor '",
              sparseHits[1L],
              "' cannot appear inside '",
              label,
              "'; a sparse column must be its own term, not wrapped in ",
              "poly(), ns(), log(), offset(), or a ':'/'*' interaction"
            )
          }
        }
        # offset() terms never reach term.labels (they carry their own
        # "offset" attribute instead, a 1-based position over the response
        # plus predictors that "variables" also carries - but "variables"
        # is the UNEVALUATED call list(response, ...), so its own [[1]] is
        # the "list" symbol and a given position's element is at [[position
        # + 1]]). A sparse name inside one is checked separately, and a
        # dense one is carried into the rewritten formula's term list
        # explicitly - reformulate has no separate 'offset' argument, but
        # the deparsed offset(...) call is valid syntax as a term label
        # like any other
        offsetPositions <- attr(expandedTerms, "offset")
        if (!is.null(offsetPositions)) {
          variables <- attr(expandedTerms, "variables")
          for (position in offsetPositions) {
            offsetTerm <- variables[[position + 1L]]
            offsetHits <- intersect(
              all.vars(offsetTerm),
              names(sparseColumns)
            )
            if (length(offsetHits) > 0L) {
              stop(
                "sparse predictor '",
                offsetHits[1L],
                "' cannot appear inside 'offset()'; a sparse column must ",
                "be its own term"
              )
            }
            denseTermLabels <- c(denseTermLabels, deparse(offsetTerm))
          }
        }
        # reformulate refuses an empty term list before R 4.3
        formula <- if (length(denseTermLabels) > 0L) {
          stats::reformulate(
            denseTermLabels,
            response = formula[[2L]],
            intercept = attr(expandedTerms, "intercept"),
            env = environment(formula)
          )
        } else {
          # what current R builds for an empty list: y ~ 1, or y ~ -1
          stats::reformulate(
            if (isTRUE(attr(expandedTerms, "intercept") == 1L)) "1" else "-1",
            response = formula[[2L]],
            env = environment(formula)
          )
        }
        usedSparseNames <- unique(usedSparseNames)
        sparseColumns <- sparseColumns[usedSparseNames]
        data <- denseData
        modelFrameCall$formula <- formula
        modelFrameCall$data <- data
      }
    } else if (!dataIsMissing && (is.list(data) || is.environment(data))) {
      refuseSparseFormulaColumns(formula, data)
    }

    # an out-of-range 'subset' would otherwise reach the na.action as a set
    # of all-NA rows and be dropped in silence
    if (!is.null(matchedCall$subset) && !dataIsMissing && is.data.frame(data)) {
      subsetIndex <- tryCatch(
        eval(matchedCall$subset, data, environment(formula)),
        error = function(e) NULL
      )
      refuseOutOfRangeSubset(subsetIndex, nrow(data), rownames(data))
    }

    # terms.formula() before R 4.3 refuses a plain list holding a bare
    # classed matrix
    if (!dataIsMissing && is.list(data) && !is.data.frame(data)) {
      data <- asDataFrameableList(data)
      modelFrameCall$data <- data
    }
    # a lifted sparse column never reaches model.frame, so the rows it holds
    # a missing value at ride along as a numeric column of NA and 0, which
    # the caller's na.action and 'subset' then treat as any other predictor
    if (length(sparseColumns) > 0L) {
      sparseMissingRows <- unique(unlist(
        lapply(sparseColumns, sparseColumnMissingRows),
        use.names = FALSE
      ))
      if (length(sparseMissingRows) > 0L) {
        sparseMissing <- numeric(nrow(data))
        sparseMissing[sparseMissingRows] <- NA_real_
        modelFrameCall$dbartsSparseMissing <- sparseMissing
      }
    }
    modelFrame <- eval(modelFrameCall, parent.frame())
    # the test frame is built from this call again, against rows this vector
    # does not cover
    modelFrameCall$dbartsSparseMissing <- NULL
    naOmitted <- attr(modelFrame, "na.action")
    # a model frame always names its rows, "1".."n" when the data has none,
    # as lm does
    trainRowNames <- rownames(modelFrame)
    # the test frame built from this call below re-reads the test data, whose
    # rows this na.action never saw
    modelFrameCall$na.action <- stats::na.pass
    if (NROW(modelFrame) == 0) {
      if (!is.null(matchedCall$subset)) {
        stop("empty 'subset' specified")
      }
      stop("cannot construct model matrices from formula")
    }

    ## pull out y - NO type coercion, so a factor response keeps its levels
    ## (model.response(., "numeric") would leave it a factor and warn, then
    ## trip "range not meaningful for factors" downstream); codeResponse then
    ## routes it exactly as the x/y path does
    y <- model.response(modelFrame)
    # a Surv response short-circuits codeResponse, which has no vocabulary
    # for it: the working response becomes the log event/censoring time (the
    # aft transform, extractSurvivalResponse's own), and the raw time/status
    # ride as attributes on the returned object (below) for whichever family
    # dbarts() resolves this into - aft reads status directly, a hazard
    # token additionally re-expands x/y by the raw time (R/dbarts.R). Read
    # off THIS already-model.frame-subsetted response - subset and na.action
    # both already applied - never re-evaluated.
    if (inherits(y, "Surv")) {
      survival <- extractSurvivalTimes(y)
      survivalStatus <- survival$status
      survivalTime <- survival$time
      # a hazard fit expands only the kept subjects, and its na.action record
      # has to name the person-period rows the dropped ones would have had,
      # which their own times decide
      if (!is.null(naOmitted)) {
        survivalTimeOmitted <- as.double(unclass(model.response(
          eval(modelFrameCall, parent.frame())
        ))[unclass(naOmitted), 1L])
      }
      y <- log(survival$time)
      responseInfo <- list(
        type = "numeric",
        n.levels = NA_integer_,
        levels = NULL
      )
    } else {
      if (is.null(y)) {
        y <- rep(0, NROW(modelFrame))
      }
      # the same by-kind guard the matrix branches carry: codeResponse
      # flattens a multi-column response column-major, and the row count
      # then reads as a mismatch against 'x' rather than naming the shape. A
      # Surv left-hand side takes its own branch just above.
      refuseMultiColumnResponse(y)
      coded <- codeResponse(y)
      y <- coded$y
      responseInfo <- coded[c("type", "n.levels", "levels")]
    }
    numObservations <- NROW(y)
    # a 'bases' entry is resolved the same way a forests = declaration's
    # basis is (validateForestBases's 'subsetRows' branch): checked against
    # the FULL pre-'subset' data and aligned to the rows the model frame
    # kept. A basis already at the model frame's own row count - the common
    # case of no 'subset' at all - passes through unchanged, since 'subset'
    # is then just the identity
    subsetRows <- if (is.null(bases)) {
      NULL
    } else {
      alignSubsetRowsToFrame(
        resolveFormulaBasisSubset(
          formula,
          if (dataIsMissing) NULL else data,
          matchedCall$subset
        ),
        naOmitted,
        numObservations
      )
    }
    bases <- validateForestBases(
      bases,
      numObservations,
      subsetRows = subsetRows
    )
    # the count response is the model's response, and a formula already names
    # one on its left-hand side; taking both would silently discard whichever
    # lost. dbarts() refuses the same combination at its own entry, this being
    # the layer that owns response ingestion. 'subset' is why the refusal
    # matters twice over: the model frame has already applied it here, so a
    # count matrix at the caller's full row count could not be aligned to the
    # rows it kept the way 'weights' is.
    if (!is.null(counts)) {
      stop(
        "'counts' is a response and cannot be given with a formula, which ",
        "names one; use the matrix interface - dbartsData(x.train, ",
        "counts = )"
      )
    }
    countsRows <- numObservations
    countsSubset <- seq_len(numObservations)

    ## weights
    weights <- as.vector(model.weights(modelFrame))
    if (!is.null(weights)) {
      if (!is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
      }
      weights <- as.double(weights)
    }

    ## offset, when in data frame
    if (identical(offsetGivenAsScalar, FALSE)) {
      offset <- as.vector(model.offset(modelFrame))
    } else if (identical(offsetGivenAsScalar, TRUE)) {
      offset <- rep_len(offset, numObservations)
    }

    ## predictors
    modelTerms <- terms(modelFrame)
    # a formula naming only sparse predictors rewrites to a dense-only RHS
    # with no terms of its own (above) - empty in exactly the sense
    # is.empty.model checks for, but not actually empty of predictors
    if (is.empty.model(modelTerms) && length(sparseColumns) == 0L) {
      stop("predictors must be specified for regression tree analysis")
    }

    termLabels <- attr(modelTerms, "term.labels")
    badLabels <- grepl("`.* .*`", termLabels)
    if (sum(badLabels) > 0) {
      termLabels[badLabels] <- gsub("^`(.*)`$", "\\1", termLabels[badLabels])
    }

    # a ':'/'*' term expands to a label like "x1:x2" in term.labels, a
    # column the model frame never carries (only the individual predictors
    # it names); passed through unchecked, makeModelMatrix dies on
    # "undefined columns selected", naming neither the term nor the
    # unsupported syntax. A label entirely wrapped in one pair of backticks
    # is a single non-syntactic name, not an interaction, and is exempt.
    interactionLabels <- termLabels[
      !grepl("^`.*`$", termLabels) & grepl(":", termLabels, fixed = TRUE)
    ]
    if (length(interactionLabels) > 0L) {
      stop(
        "':' and '*' terms are not supported in 'formula'; write each ",
        "predictor as its own term - poly(), ns(), log(), and offset() are ",
        "supported"
      )
    }

    predictorFrame <- modelFrame[termLabels]
    if (length(sparseColumns) > 0L) {
      # rownames(modelFrame) is character; a sparse column carries no row
      # names of its own, so its rows are resolved by matching the model
      # frame's back into the (already sparse-column-pulled) 'data' this
      # sparse column itself still indexes by - a match that aligns under
      # 'subset' and na.action together, since both already shaped
      # modelFrame's own rows by the time this runs
      pos <- match(rownames(modelFrame), rownames(data))
      for (sparseName in names(sparseColumns)) {
        predictorFrame[[sparseName]] <-
          subsetSparseColumn(sparseColumns[[sparseName]], pos)
      }
    }
    x <- makeModelMatrix(predictorFrame)

    if (!testIsMissing) {
      testCall <- matchedCall
      testCall <- testCall[c(
        1L,
        match(c("formula", "data", "test"), names(testCall), nomatch = 0L)
      )]
      names(testCall)[which(names(testCall) == "test")] <- "term"
      testCall[[1L]] <- quoteInNamespace(findTermInFormulaData)

      temp <- eval(testCall, parent.frame())
      if (!is.null(temp)) test <- temp
    }
  } else if (inherits(formula, "dgCMatrix")) {
    ## sparse designs enter through the x/y interface only; columns are all
    ## ordinal and missing values are stored NaN entries (the Matrix
    ## convention), so no complete-case filtering applies
    if (dataIsMissing || is.null(data)) {
      data <- rep(0, nrow(formula))
    }
    if (
      !is.numeric(data) &&
        !is.factor(data) &&
        !is.logical(data) &&
        !is.character(data)
    ) {
      stop(
        "when 'formula' is a sparse matrix, 'data' must be numeric, a ",
        "factor, logical, or character"
      )
    }

    refuseMultiColumnResponse(data)
    coded <- codeResponse(data)
    y <- coded$y
    responseInfo <- coded[c("type", "n.levels", "levels")]
    if (nrow(formula) != NROW(y)) {
      stop("'x' must have the same number of observations as 'y'")
    }
    initialNumObservations <- NROW(y)
    # an empty training set has no cut grid to quantize and would fault deeper
    # (a subscript-out-of-bounds subsetting the zero rows below); name it here,
    # as the formula path already rejects an empty model frame
    if (initialNumObservations == 0L) {
      stop("data has zero rows")
    }

    if (missing(subset) || is.null(subset)) {
      subset <- seq.int(length(y))
    }
    refuseOutOfRangeSubset(subset, initialNumObservations)
    y <- y[subset]
    x <- formula[subset, , drop = FALSE]
    trainRowNames <- rownames(x)
    bases <- validateForestBases(bases, initialNumObservations, subset)
    countsRows <- initialNumObservations
    countsSubset <- subset

    if (missing(weights)) {
      weights <- NULL
    }
    weights <- validateXYWeights(weights, initialNumObservations, subset)

    if (offsetIsMissing) {
      offset <- NULL
    } else if (is.matrix(offset) || is.data.frame(offset)) {
      if (is.null(counts)) {
        refuseMatrixOffset(offset, "offset")
      }
      categoryOffset <- if (is.data.frame(offset)) as.matrix(offset) else offset
      offset <- NULL
    }
    offsetResult <- validateXYOffset(
      offset,
      initialNumObservations,
      subset,
      offsetGivenAsScalar
    )
    offset <- offsetResult$offset
    offsetGivenAsScalar <- offsetResult$offsetGivenAsScalar

    # the same (y, x) row rule the dense branch applies; a sparse container
    # holds its missing values among the STORED entries, which is where
    # rowsWithMissingPredictors looks
    naResult <- applyNaActionToXY(
      na.action,
      multinomialResponseNA(counts, y),
      x
    )
    if (!is.null(naResult)) {
      naOmitted <- nameOmittedRows(naResult$na.action, trainRowNames)
      if (!all(naResult$keep)) {
        keep <- naResult$keep
        trainRowNames <- trainRowNames[keep]
        y <- y[keep]
        x <- x[keep, , drop = FALSE]
        bases <- restrictBasesToRows(bases, keep)
        if (!is.null(weights)) {
          weights <- weights[keep]
        }
        if (!is.null(offset)) {
          offset <- offset[keep]
        }
        # 'counts' (family = "multinomial") is validated and row-selected
        # against this same subset, below, once it is known; it loses
        # exactly the rows y and x did
        countsSubset <- countsSubset[keep]
      }
    }
  } else if (
    is.numeric(formula) || is.data.frame(formula) || is.factor(formula)
  ) {
    ## backwards compatibility of bart(x.train, y.train, x.test)
    if (dataIsMissing || is.null(data)) {
      data <- rep(0, NROW(formula))
    }
    if (
      !is.numeric(data) &&
        !is.data.frame(data) &&
        !is.factor(data) &&
        !is.logical(data) &&
        !is.character(data)
    ) {
      stop(
        "when 'formula' is numeric, 'data' must be numeric, a factor, ",
        "logical, or character"
      )
    }

    refuseMultiColumnResponse(data)
    coded <- codeResponse(data)
    y <- coded$y
    responseInfo <- coded[c("type", "n.levels", "levels")]
    if (NROW(formula) != NROW(y)) {
      stop("'x' must have the same number of observations as 'y'")
    }
    initialNumObservations <- NROW(y)
    # an empty training set has no cut grid to quantize and would fault deeper
    # (a subscript-out-of-bounds subsetting the zero rows below); name it here,
    # as the formula path already rejects an empty model frame
    if (initialNumObservations == 0L) {
      stop("data has zero rows")
    }

    # a data frame's automatic names stay a deferred string until indexed,
    # so they are indexed only when rows are dropped
    trainRowNames <- observationRowNames(formula)
    if (missing(subset) || is.null(subset)) {
      subset <- seq.int(length(y))
    } else if (!is.null(trainRowNames)) {
      trainRowNames <- trainRowNames[subset]
    }
    refuseOutOfRangeSubset(subset, initialNumObservations)
    y <- y[subset]

    if (is.data.frame(formula)) {
      formula <- makeModelMatrix(formula)
    }
    xIsMixed <- inherits(formula, "dbartsMixedMatrix")
    x <- if (is.matrix(formula) || xIsMixed) {
      formula[subset, , drop = FALSE]
    } else {
      formula[subset]
    }
    bases <- validateForestBases(bases, initialNumObservations, subset)
    countsRows <- initialNumObservations
    countsSubset <- subset

    if (missing(weights)) {
      weights <- NULL
    }
    weights <- validateXYWeights(weights, initialNumObservations, subset)

    if (offsetIsMissing) {
      offset <- NULL
    } else if (is.matrix(offset) || is.data.frame(offset)) {
      if (is.null(counts)) {
        refuseMatrixOffset(offset, "offset")
      }
      categoryOffset <- if (is.data.frame(offset)) as.matrix(offset) else offset
      offset <- NULL
    }
    offsetResult <- validateXYOffset(
      offset,
      initialNumObservations,
      subset,
      offsetGivenAsScalar
    )
    offset <- offsetResult$offset
    offsetGivenAsScalar <- offsetResult$offsetGivenAsScalar

    # the (y, x) pair the matrix interface supplies is the model frame the
    # na.action reads; a mixed container keeps its own attributes across the
    # row selection, so it takes the shared subsetting below rather than this
    # branch's attribute-preserving one
    naResult <- applyNaActionToXY(
      na.action,
      multinomialResponseNA(counts, y),
      x
    )
    if (!is.null(naResult)) {
      naOmitted <- nameOmittedRows(naResult$na.action, trainRowNames)
      if (!is.null(trainRowNames) && !all(naResult$keep)) {
        trainRowNames <- trainRowNames[naResult$keep]
      }
    }
    if (!xIsMixed) {
      # NULL means nothing was missing anywhere (applyNaActionToXY's own
      # contract), and an all-TRUE mask selects every row: either way the
      # row selection below would copy y, x, the bases, the weights and the
      # offset to no effect. x is n*p doubles, the largest transient
      # ingestion has, so the common case skips it rather than paying it.
      completeCases <- if (is.null(naResult)) NULL else naResult$keep
      if (!is.null(completeCases) && !all(completeCases)) {
        y <- y[completeCases]
        x <- if (!is.matrix(x)) {
          x[completeCases]
        } else {
          x[completeCases, , drop = FALSE]
        }
        bases <- restrictBasesToRows(bases, completeCases)
        if (!is.null(weights)) {
          weights <- weights[completeCases]
        }
        if (!is.null(offset)) {
          offset <- offset[completeCases]
        }
        # 'counts' (family = "multinomial") is validated and row-selected
        # against this same subset, below, once it is known; it loses
        # exactly the rows y and x did
        countsSubset <- countsSubset[completeCases]
      }
      if (length(attributes(formula)) > 0L) {
        for (attributeName in names(attributes(formula))) {
          if (attributeName == "dim") {
            next
          }
          if (attributeName == "dimnames" && !identical(dim(formula), dim(x))) {
            next
          }
          attr(x, attributeName) <- attr(formula, attributeName)
        }
      }
    } else if (!is.null(naResult) && !all(naResult$keep)) {
      keep <- naResult$keep
      y <- y[keep]
      x <- x[keep, , drop = FALSE]
      bases <- restrictBasesToRows(bases, keep)
      if (!is.null(weights)) {
        weights <- weights[keep]
      }
      if (!is.null(offset)) {
        offset <- offset[keep]
      }
      countsSubset <- countsSubset[keep]
    }
  } else {
    stop(
      "unrecognized 'formula' type; must be coercible to numeric or a valid formula object"
    )
  }

  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  if (is.data.frame(x)) {
    x <- makeModelMatrix(x)
  }

  x.test <- NULL
  testRowNames <- NULL
  if (!testIsMissing && !is.null(test)) {
    testRowNames <- observationRowNames(test)
    x.test <- validateXTest(test, x)
  }

  if (!is.null(x.test)) {
    if (testOffsetIsMissing) {
      ## default is offset.test = offset
      if (identical(offsetGivenAsScalar, TRUE)) {
        offset.test <- rep_len(offset[1L], nrow(x.test))
        testUsesRegularOffset <- TRUE
      } else if (identical(offsetGivenAsScalar, FALSE)) {
        if (nrow(x.test) != length(y)) {
          stop(
            "vectored 'offset' cannot be directly applied to test data of unequal length"
          )
        }
        offset.test <- offset
        testUsesRegularOffset <- TRUE
      }
    } else {
      testOffsetInfo <- eval(getTestOffset)

      offset.test <- testOffsetInfo$offset.test
      testUsesRegularOffset <- testOffsetInfo$testUsesRegularOffset

      if (
        !is.null(offset.test) &&
          !is.matrix(offset.test) &&
          !is.data.frame(offset.test)
      ) {
        offset.test <- rep_len(offset.test, nrow(x.test))
      }
    }
  } else {
    if (testOffsetIsMissing) offset.test <- NULL
  }
  if (
    !is.null(offset.test) &&
      (is.matrix(offset.test) || is.data.frame(offset.test))
  ) {
    if (is.null(counts)) {
      refuseMatrixOffset(offset.test, "offset.test")
    }
    categoryTestOffset <- if (is.data.frame(offset.test)) {
      as.matrix(offset.test)
    } else {
      offset.test
    }
    offset.test <- NULL
  }

  weights.test <- NULL
  if (!is.null(x.test) && !is.null(matchedCall$weights)) {
    if (!is.formula(formula)) {
      warning(warningCondition(
        "'weights' are ignored for test data when model is not specified as a formula; this only impacts extracting samples from the posterior predictive distribution of the test data",
        class = c("dbartsIgnoredArgWarning", "dbartsWarning")
      ))
    } else {
      testFormula <- formula
      lhs <- testFormula[[2L]]
      remainder <- testFormula
      remainder[[2L]] <- NULL
      testFormula <- as.formula(paste0(deparse(remainder), " - ", deparse(lhs)))
      environment(testFormula) <- environment(formula)
      modelFrameCall$formula <- testFormula
      modelFrameCall$data <- test
      tryResult <- tryCatch(
        testFrame <- eval(modelFrameCall, parent.frame()),
        error = function(e) e
      )
      if (inherits(tryResult, "error")) {
        warning(warningCondition(
          "weights specified but not found in test data - ignoring",
          class = c("dbartsIgnoredArgWarning", "dbartsWarning")
        ))
      } else {
        weights.test <- testFrame[["(weights)"]]
      }
    }
  }

  # missingness is a predictor-only feature: rules route NAs in x, but the
  # response side must be complete. In a sparse x, NAs live only among the
  # stored entries and implicit zeros are observed values, so the checks
  # work off the slots without densifying. An out-of-range 'subset' is one
  # way to reach this silently: data.frame row-indexing pads unmatched rows
  # with NA rather than erroring, so a NA response here can be that instead
  # of a genuinely incomplete row - naming 'subset' when it was given is a
  # cheap, cause-agnostic hint rather than a claim of certainty.
  # The multinomial response: the n x K count
  # matrix IS the response, and 'y' is its trials vector n_i = sum_k
  # counts[i, k], which is what keeps every length(data@y) reader meaningful on
  # such an object. DERIVED rather than taken, so the two cannot disagree and a
  # caller supplying 'counts' supplies no separate response.
  counts <- validateMultinomialCounts(
    counts,
    countsRows,
    countsSubset,
    allowMissing = TRUE
  )
  if (!is.null(counts)) {
    y <- as.double(rowSums(counts))
    # a flat offset was already refused at the point it was set aside if
    # 'counts' was absent there; this is the twin case, a flat offset
    # alongside counts that WAS supplied
    if (!is.null(offset)) {
      refuseFlatOffsetOnMultinomial(offset)
    }
    if (!is.null(offset.test)) {
      refuseFlatOffsetOnMultinomial(offset.test, "offset.test")
    }
  }
  offset.category <- validateDataCategoryOffset(
    categoryOffset,
    countsRows,
    countsSubset,
    "offset",
    if (is.null(counts)) NULL else ncol(counts)
  )
  # the test twin's rows are the TEST rows, so it is not subset with the
  # training ones and its row count is pinned by the validity method
  offset.category.test <- validateDataCategoryOffset(
    categoryTestOffset,
    NULL,
    NULL,
    "offset.test",
    if (is.null(counts)) NULL else ncol(counts)
  )
  if (!is.null(offset.category.test) && is.null(x.test)) {
    stop("'offset.test' must be null when 'test' is null")
  }

  if (anyNA(y)) {
    if (!is.null(matchedCall$subset)) {
      stop(
        "response contains missing values; check that 'subset' selects ",
        "rows within range"
      )
    }
    stop("response contains missing values")
  }
  # Inf/-Inf survive the anyNA check (NaN does not), then poison the
  # precision-degeneracy ratio below into an NA condition; reject them here
  # with a named error instead.
  if (any(is.infinite(y))) {
    stop("response contains non-finite values")
  }

  # Precision-degenerate response: a large magnitude but tiny spread
  # quantizes to (near-)identical doubles before the engine ever sees it,
  # e.g. y in [1e15, 1e15 + 1e-3] rounds to one representable double. Doubles
  # near magnitude s are spaced ~2.22e-16 * s apart, so the 1e-10 threshold is
  # ~1e6x the ulp spacing, clear of legitimately discrete (binary, counts,
  # ordinal) responses. max(abs(y)) == 0 is guarded to avoid a 0/0. A
  # multinomial 'y' is the trials vector, not a modelled response, so the
  # check is skipped whenever 'counts' is carried, keyed on the data object
  # rather than on a family dbartsData has no formal for.
  yRange <- diff(range(y))
  yScale <- max(abs(y))
  if (is.null(counts) && yScale > 0 && yRange / yScale < 1e-10) {
    warning(warningCondition(
      paste0(
        "response values are indistinguishable, or nearly so, at double ",
        "precision (",
        length(unique(y)),
        " distinct value(s) among ",
        length(y),
        " observations); center and/or rescale the response before fitting"
      ),
      class = c("dbartsDegenerateResponseWarning", "dbartsWarning")
    ))
  }

  sparseAllMissingCheck <- function(x.sparse) {
    columnNnz <- diff(x.sparse@p)
    columnNumNA <- vapply(
      seq_len(ncol(x.sparse)),
      function(j) {
        sum(is.na(x.sparse@x[seq.int(
          x.sparse@p[j] + 1L,
          length.out = columnNnz[j]
        )]))
      },
      0L
    )
    if (any(columnNnz == nrow(x.sparse) & columnNumNA == nrow(x.sparse))) {
      stop("predictor columns cannot be entirely missing")
    }
  }
  if (is.matrix(x)) {
    xHasNA <- anyNA(x)
    if (xHasNA && any(colSums(!is.na(x)) == 0L)) {
      stop("predictor columns cannot be entirely missing")
    }
  } else if (inherits(x, "dbartsMixedMatrix")) {
    # both flavors hold dense columns as a per-column list; the mixed flavor
    # adds a sparse part
    denseColumns <- if (is.null(x$dense)) list() else x$dense
    denseHasNA <- any(vapply(denseColumns, anyNA, FALSE))
    sparseHasNA <- !is.null(x$sparse) && anyNA(x$sparse@x)
    xHasNA <- denseHasNA || sparseHasNA
    if (xHasNA) {
      denseAllMissing <- vapply(
        denseColumns,
        function(column) all(is.na(column)),
        FALSE
      )
      if (any(denseAllMissing)) {
        stop("predictor columns cannot be entirely missing")
      }
      if (!is.null(x$sparse)) {
        sparseAllMissingCheck(x$sparse)
      }
    }
  } else {
    xHasNA <- anyNA(x@x)
    if (xHasNA) sparseAllMissingCheck(x)
  }
  if (!is.null(offset) && anyNA(offset)) {
    stop("'offset' contains missing values")
  }
  if (!is.null(offset.test) && anyNA(offset.test)) {
    stop("'offset.test' contains missing values")
  }

  result <- newValidated(
    "dbartsData",
    modelMatrices = namedList(
      y,
      x,
      x.test,
      weights,
      weights.test,
      offset,
      offset.test,
      bases,
      counts,
      offset.category,
      offset.category.test,
      testUsesRegularOffset
    ),
    n.cuts = NA_integer_,
    sigma = NA_real_
  )
  result@na.action <- naOmitted
  if (!is.null(trainRowNames) || !is.null(testRowNames)) {
    result@rowNames <- list(train = trainRowNames, test = testRowNames)
  }
  result@response.type <- responseInfo$type
  result@response.n.levels <- as.integer(responseInfo$n.levels)
  result@response.levels <- responseInfo$levels
  # a Surv formula response's raw time/status, for dbarts() to decode
  # (above) - not a slot, since every other creation route (the matrix
  # interface, a pre-built dbartsData) has no use for it and carries it
  # through its own channel instead (control@bartcore.survival)
  if (!is.null(survivalStatus)) {
    attr(result, "survivalStatus") <- survivalStatus
    attr(result, "survivalTime") <- survivalTime
    attr(result, "survivalTimeOmitted") <- survivalTimeOmitted
  }
  result
}
