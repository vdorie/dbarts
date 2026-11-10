ORDINAL_VARIABLE <- 0L
CATEGORICAL_VARIABLE <- 1L
ORDERED_FACTOR_VARIABLE <- 2L

## The package's own na.action, and the default of every fitting function
## that takes one. It drops the rows whose RESPONSE is missing and keeps the
## rows whose predictors are: BART routes a missing predictor value down a
## learned side of each rule, so an incomplete row is data rather than a
## hole, while a missing response is nothing to fit. rpart's na.rpart is the
## behavioral precedent. The dropped rows are recorded exactly as
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
## The opening of dbartsData's precision-degenerate response warning, which
## withBinaryResponsePrecision recognizes it by.
responsePrecisionWarningStem <- "response values are indistinguishable"

applyNaActionToXY <- function(na.action, y, x) {
  predictorNA <- rowsWithMissingPredictors(x)
  responseNA <- is.na(y)
  if (!any(predictorNA) && !any(responseNA)) {
    return(NULL)
  }
  # unnamed: the caller's row names can hold a missing one, which
  # data.frame refuses; nameOmittedRows relabels the record from them
  frame <- data.frame(
    response = unname(ifelse(responseNA, NA_real_, 0.0)),
    predictors = unname(ifelse(predictorNA, NA_real_, 0.0))
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
  # a vector is not a count matrix: left for validateMultinomialCounts to
  # refuse by name, as a row-count mismatch is
  if (
    is.null(counts) || length(dim(counts)) != 2L || NROW(counts) != length(y)
  ) {
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

## The rows `rows` of a basis as it was handed over, a vector or a matrix.
cutBasisRows <- function(basis, rows) {
  if (is.null(basis)) {
    NULL
  } else if (length(dim(basis)) < 2L) {
    basis[rows]
  } else {
    basis[rows, , drop = FALSE]
  }
}

## An index past the data's rows is silent in base R: row indexing pads it
## with NA rather than erroring, and the na.action would then drop exactly
## that row and fit fewer observations than the caller asked for without a
## word. Named here instead, ahead of the model frame. A missing entry of a
## logical or a numeric subset is read as lm reads it, the row it pads left
## to the na.action; a missing name matches no row and is refused with the
## names the data lacks. `n` of NA means the row count is not known this
## early, which leaves the check unrun rather than guessed at.
refuseOutOfRangeSubset <- function(index, n, rowNames = NULL) {
  if (is.null(index) || is.na(n)) {
    return(invisible(NULL))
  }
  outOfRange <- if (is.logical(index)) {
    length(index) > n
  } else if (is.numeric(index)) {
    any(abs(index[!is.na(index)]) > n)
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

## The positions among n rows that a 'subset' given as positions, a logical
## or row names selects, so that a later cut by the na.action reads them.
subsetPositions <- function(index, n, rowNames = NULL) {
  if (is.character(index)) {
    return(match(index, rowNames))
  }
  seq_len(n)[index]
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

## The predictor columns a data object's sampler has been given a missing
## value in: NULL, or a logical per column. A data object saved before the
## slot existed has none.
dataMissingSeen <- function(data) {
  if (!methods::.hasSlot(data, "missing.seen")) {
    return(NULL)
  }
  data@missing.seen
}

## The levels a data object records for its forests' bases: NULL, or a list
## with one entry per forest. A data object saved before the slot existed has
## none, and its forests read as multiplying numbers.
dataBasisLevels <- function(data) {
  if (!methods::.hasSlot(data, "basis.levels")) {
    return(NULL)
  }
  data@basis.levels
}

## The scale a data object records for each numeric forest's column: NULL, or
## a list with one entry per forest. A data object saved before the slot
## existed has none.
dataBasisScale <- function(data) {
  if (!methods::.hasSlot(data, "basis.scale")) {
    return(NULL)
  }
  data@basis.scale
}

## What a data object's record must hold for as long as a sampler has it: a
## forest's recorded levels are as many as its block has columns, and its
## basis is of a size this version takes (basisSizeFault), read from the
## levels and the block's width. Nothing here looks at a row, a swap being
## free to leave a level or a column with none. Returns the fault as text
## where `refuse` is FALSE, for the class's validity.
checkBasisRecord <- function(data, refuse = TRUE) {
  levels <- dataBasisLevels(data)
  for (index in seq_len(max(length(levels), length(data@bases)))) {
    block <- if (index <= length(data@bases)) data@bases[[index]]
    recorded <- if (index <= length(levels)) levels[[index]]
    width <- if (is.null(block)) 0L else NCOL(block)
    fault <- if (!is.null(recorded) && length(recorded) != width) {
      sprintf(
        "'basis.levels' gives forest %d %d %s and its block in 'bases' has %d %s; dbartsData() and dbarts() write the two together",
        index,
        length(recorded),
        if (length(recorded) == 1L) "level" else "levels",
        width,
        if (width == 1L) "column" else "columns"
      )
    } else {
      basisSizeFault(block, recorded, paste("forest", index))
    }
    if (!is.null(fault)) {
      if (refuse) {
        stop(fault, call. = FALSE)
      }
      return(fault)
    }
  }
  NULL
}

## Brings a sampler's record up to its engine's flags, which only the engine
## raises; called after every call that can raise one. Written only where it
## moved, a write of the data field costing more than the read.
recordMissingSeen <- function(sampler, ptr) {
  seen <- .Call(C_dbarts_bartcore_getMissingSeen, ptr)
  if (!identical(seen, dataMissingSeen(sampler$data))) {
    sampler$dataHeld@missing.seen <- seen
  }
  invisible(NULL)
}

## The derived values a data object records for its sampler: NULL where it
## records none, as one saved before the slots existed does.
dataCutPoints <- function(data) {
  if (!methods::.hasSlot(data, "cutPoints")) {
    return(NULL)
  }
  data@cutPoints
}
dataResponseRange <- function(data) {
  if (!methods::.hasSlot(data, "response.range")) {
    return(NULL)
  }
  data@response.range
}
dataForestCalibration <- function(data) {
  if (!methods::.hasSlot(data, "forest.calibration")) {
    return(NULL)
  }
  data@forest.calibration
}
dataStandardization <- function(data) {
  if (!methods::.hasSlot(data, "standardization")) {
    return(NULL)
  }
  data@standardization
}

## The names of the columns a model's linear or gp leaves take, NA for one
## whose name cannot key the record: missing, empty, or not that of exactly
## one column, the test setData's alignment matches names by; NULL where the
## leaves take none or the columns have no names. The record is keyed by the
## names where none is NA and otherwise by the columns' positions, so no name
## key can match another column's entry.
leafCovariateNames <- function(model, data) {
  prior <- model@leaf.prior
  if (is(prior, "dbartsLinearPrior") || is(prior, "dbartsGPPrior")) {
    names <- colnames(data@x)
    if (!is.null(names)) {
      once <- !is.na(names) &
        nzchar(names) &
        !(duplicated(names) | duplicated(names, fromLast = TRUE))
      names[!once] <- NA_character_
      names[prior@columns]
    }
  }
}
## Whether names can key a record: present, non-empty and distinct.
usableKeys <- function(names) {
  !is.null(names) &&
    !anyNA(names) &&
    all(nzchar(names)) &&
    !anyDuplicated(names)
}

## Holds the derived values a data object records to the data before a
## sampler is created from it, and returns it with each grid in order, a
## count family's range as the window (c, c + 1) of its first number and the
## standardization matched to the model's leaves, by their names in
## covariates. Only the record's form, the columns' kinds and the leaf
## covariates are read, never a row: a record is the values to use, whatever
## data they came from. Each refusal names the slot and says that NULL
## derives the value instead.
checkDerivedRecord <- function(
  data,
  model,
  covariates = leafCovariateNames(model, data)
) {
  refuse <- refuseDerivedRecord
  data <- checkCutPointsRecord(data)
  # a count family's record is its center c, held as (c, c + 1), which is
  # never an equal pair; a pair of another width, as a range recorded under
  # another family is, is read by its first number and recorded in that form
  # at every door, so the sampler's record is the window it holds
  range <- dataResponseRange(data)
  count <- identical(model@family, "nbinom")
  if (
    !is.null(range) &&
      (!is.numeric(range) ||
        length(range) != 2L ||
        !all(is.finite(range)) ||
        range[[2L]] < range[[1L]] ||
        (count && range[[2L]] == range[[1L]]))
  ) {
    refuse(
      "response.range",
      "must be two finite numbers, the second ",
      if (count) "above" else "not below",
      " the first"
    )
  }
  if (count && !is.null(range)) {
    data@response.range <- as.double(range[[1L]]) + c(0, 1)
  }
  calibration <- dataForestCalibration(data)
  if (
    !is.null(calibration) &&
      (!is.numeric(calibration) ||
        length(calibration) != 2L ||
        !all(is.finite(calibration)) ||
        any(calibration <= 0) ||
        !(is.null(names(calibration)) ||
          identical(names(calibration), c("anchor", "unit"))))
  ) {
    refuse(
      "forest.calibration",
      "must be two positive finite numbers, c(anchor = , unit = )"
    )
  }
  standardization <- dataStandardization(data)
  if (!is.null(standardization)) {
    data@standardization <- checkStandardizationRecord(
      standardization,
      model,
      covariates,
      refuse
    )
  }
  data
}

refuseDerivedRecord <- function(slot, ...) {
  stop(
    "the data's '",
    slot,
    "' record ",
    ...,
    "; set it to NULL to derive it from the data",
    call. = FALSE
  )
}

## A data object's grid record held to its columns' kinds and put in the order a
## grid takes, as a sampler or xbart()'s store is built on it.
checkCutPointsRecord <- function(data) {
  refuse <- refuseDerivedRecord
  grids <- dataCutPoints(data)
  if (!is.null(grids)) {
    varTypes <- data@varTypes
    if (!is.list(grids) || length(grids) != length(varTypes)) {
      refuse(
        "cutPoints",
        "must be a list with one entry per predictor column, ",
        length(varTypes),
        " here"
      )
    }
    for (j in seq_along(grids)) {
      grid <- grids[[j]]
      if (varTypes[[j]] == CATEGORICAL_VARIABLE) {
        if (!is.null(grid)) {
          refuse(
            "cutPoints",
            "gives predictor ",
            j,
            ", an unordered factor, cut points; its entry must be NULL"
          )
        }
        next
      }
      if (!is.numeric(grid) || length(grid) == 0L || !all(is.finite(grid))) {
        refuse(
          "cutPoints",
          "for predictor ",
          j,
          " must be one or more finite numbers"
        )
      }
      mass <- attr(grid, "mass", exact = TRUE)
      grid <- as.double(grid)
      if (
        !is.null(mass) &&
          (!is.numeric(mass) ||
            length(mass) != length(grid) + 1L ||
            !all(is.finite(mass)) ||
            any(diff(mass) <= 0))
      ) {
        refuse(
          "cutPoints",
          "for predictor ",
          j,
          " carries weights, its 'mass' attribute, that are not ",
          length(grid) + 1L,
          " finite numbers strictly increasing"
        )
      }
      # a grid is a set of thresholds, so one out of order is sorted, as
      # setCutPoints sorts one, each point keeping its weight: the mass is
      # cumulative, point k weighing mass[k + 1] - mass[k]
      if (is.unsorted(grid)) {
        ordering <- order(grid)
        grid <- grid[ordering]
        if (!is.null(mass)) {
          mass <- mass[[1L]] + c(0, cumsum(diff(mass)[ordering]))
        }
      }
      if (anyDuplicated(grid) > 0L) {
        refuse(
          "cutPoints",
          "for predictor ",
          j,
          " repeats a point (",
          grid[[anyDuplicated(grid)]],
          "): a cut grid holds each point once"
        )
      }
      if (!is.null(mass)) {
        attr(grid, "mass") <- as.double(mass)
      }
      grids[[j]] <- grid
    }
    data@cutPoints <- grids
  }
  data
}

## The standardization record's part of checkDerivedRecord: a list of one
## entry per forest, each NULL or list(center = , scale = ) with, for gp
## leaves, lengthscale beside them, as many of each: a finite center, a
## positive or NA scale, a logical NA read as one, and a positive lengthscale,
## keyed by the names of center, present, non-empty and distinct, or by
## columns, their positions. Any other form is refused. Which covariates the
## leaves take is the model's, so the rest is matched and never refused: only
## the first forest's leaves can take any, and its entry goes to them by key,
## names to leaf columns whose names can be keys, as covariates gives them,
## the others deriving, and positions where any leaf column's name cannot be
## a key (leafCovariateNames); an entry with neither, as only a hand can
## write, goes by position where the counts agree. An entry for a covariate
## the leaves do not take is ignored, and a covariate no entry matches gets an
## NA center, which derives it. Returns the matched record, NULL where nothing
## matches.
checkStandardizationRecord <- function(record, model, covariates, refuse) {
  if (is.null(record)) {
    return(NULL)
  }
  given <- if (is.list(record)) !vapply(record, is.null, NA)
  if (
    !is.list(record) ||
      length(record) == 0L ||
      !all(vapply(record[given], is.list, NA))
  ) {
    refuse(
      "standardization",
      "must be a list with one entry per forest, each NULL or a list of ",
      "'center' and 'scale'"
    )
  }
  numbersOf <- function(value) {
    if (is.logical(value) && length(value) > 0L && all(is.na(value))) {
      as.double(value)
    } else {
      value
    }
  }
  for (entry in record[given]) {
    center <- entry$center
    scale <- numbersOf(entry$scale)
    lengthscale <- entry$lengthscale
    positions <- entry$columns
    keys <- names(center)
    q <- length(center)
    if (
      !is.numeric(center) ||
        q == 0L ||
        !all(is.finite(center)) ||
        !is.numeric(scale) ||
        length(scale) != q ||
        !all(is.na(scale) | (is.finite(scale) & scale > 0)) ||
        !(is.null(lengthscale) ||
          (is.numeric(lengthscale) &&
            length(lengthscale) == q &&
            all(is.finite(lengthscale) & lengthscale > 0))) ||
        !(is.null(positions) ||
          (is.numeric(positions) &&
            length(positions) == q &&
            all(is.finite(positions) & positions >= 1) &&
            all(positions == round(positions)))) ||
        !(is.null(keys) || usableKeys(keys))
    ) {
      refuse(
        "standardization",
        "must give a finite 'center', a positive or NA 'scale' and, for gp ",
        "leaves, a positive 'lengthscale', one value each per leaf covariate, ",
        "keyed by distinct non-empty names or by 'columns', their positions"
      )
    }
  }
  entry <- record[[1L]]
  prior <- model@leaf.prior
  if (
    is.null(entry) ||
      !(is(prior, "dbartsLinearPrior") || is(prior, "dbartsGPPrior"))
  ) {
    return(NULL)
  }
  q <- length(prior@columns)
  keys <- names(entry$center)
  none <- rep(NA_integer_, q)
  matched <- if (!is.null(keys)) {
    if (!is.null(covariates)) match(covariates, keys) else none
  } else if (!is.null(entry$columns)) {
    if (is.null(covariates) || anyNA(covariates)) {
      match(prior@columns, entry$columns)
    } else {
      none
    }
  } else if (length(entry$center) == q) {
    seq_len(q)
  } else {
    none
  }
  if (all(is.na(matched))) {
    return(NULL)
  }
  held <- list(
    center = as.double(entry$center)[matched],
    scale = as.double(numbersOf(entry$scale))[matched]
  )
  if (!is.null(entry$lengthscale)) {
    held$lengthscale <- as.double(entry$lengthscale)[matched]
  }
  list(held)
}

## Refuses data whose record gives an ordered factor other cut points than the
## ones the engine holds for it, those of its levels, which the record cannot
## move; records holds the engine's (C_dbarts_bartcore_derivedRecords), or the
## sampler's record of them at setData.
refuseCarriedOrderedGrid <- function(data, records) {
  grids <- dataCutPoints(data)
  for (j in which(data@varTypes == ORDERED_FACTOR_VARIABLE)) {
    if (j > length(grids) || is.null(grids[[j]])) {
      next
    }
    if (!identical(grids[[j]], records$cutPoints[[j]])) {
      stop(
        "the data's 'cutPoints' record gives predictor ",
        j,
        ", an ordered factor, other cut points than those of its levels; ",
        "set it to NULL to derive it from the data",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

## What an ordered factor's cut points are, from its levels alone, for a caller
## with no sampler to ask: one between each pair, by level code.
orderedFactorGrids <- function(data) {
  grids <- vector("list", length(data@varTypes))
  levelTables <- attr(data@x, "factor.levels")
  for (j in which(data@varTypes == ORDERED_FACTOR_VARIABLE)) {
    numLevels <- length(levelTables[[j]])
    grids[[j]] <- if (numLevels >= 2L) {
      seq_len(numLevels - 1L) - 0.5
    } else {
      dataCutPoints(data)[[j]]
    }
  }
  grids
}

## Brings a sampler's record of what it derived from its data up to its
## engine; called after creation and after every call that can change a grid,
## recordResponseRange after one that can move only the range.
## The object a caller handed in is never written, the record being on the
## sampler's own copy.
recordDerived <- function(sampler, ptr, records = NULL) {
  if (is.null(records)) {
    records <- .Call(C_dbarts_bartcore_derivedRecords, ptr, TRUE)
  }
  data <- sampler$data
  data@cutPoints <- records$cutPoints
  data@response.range <- records$response.range
  data@forest.calibration <- records$forest.calibration
  data@standardization <- keyStandardization(
    sampler$model,
    data,
    records$standardization
  )
  sampler$dataHeld <- data
  invisible(NULL)
}

## The engine's standardization as a data object records it: keyed by the
## leaf covariates' names, or by their positions where any of those names
## cannot be a key.
keyStandardization <- function(model, data, standardization) {
  if (is.null(standardization)) {
    return(NULL)
  }
  covariates <- leafCovariateNames(model, data)
  standardization[[1L]] <- if (!is.null(covariates) && !anyNA(covariates)) {
    lapply(standardization[[1L]], `names<-`, covariates)
  } else {
    c(
      standardization[[1L]],
      list(columns = as.integer(model@leaf.prior@columns))
    )
  }
  standardization
}

## Brings a sampler's record of the response range, and with standardization
## its leaf covariates' standardization, up to its engine after a call that
## can move them and no grid: a response, offset, weights or active-row change
## with updateScale or updateStandardization, which can run every iteration of
## a larger sampler. A calibration map's anchor and unit are fixed at
## construction, so none of those calls moves them, and the grids, whose read
## grows with the predictors, are not read.
recordResponseRange <- function(
  sampler,
  ptr,
  range = TRUE,
  standardization = FALSE
) {
  records <- .Call(C_dbarts_bartcore_derivedRecords, ptr, FALSE)
  if (
    range && !identical(records$response.range, dataResponseRange(sampler$data))
  ) {
    sampler$dataHeld@response.range <- records$response.range
  }
  if (standardization) {
    keyed <- keyStandardization(
      sampler$model,
      sampler$data,
      records$standardization
    )
    if (!identical(keyed, sampler$data@standardization)) {
      sampler$dataHeld@standardization <- keyed
    }
  }
  invisible(NULL)
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

## Whether a multinomial data object's response was written as a factor (TRUE)
## or as a count matrix (FALSE); NA where it carries no counts, and for an
## object saved before the slot existed.
dataCountsAsFactor <- function(data) {
  if (!methods::.hasSlot(data, "counts.as.factor")) {
    return(NA)
  }
  data@counts.as.factor
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
      .Object@offset <- integerToDouble(modelMatrices$offset)
      .Object@offset.test <- integerToDouble(modelMatrices$offset.test)
      .Object@bases <- modelMatrices$bases
      .Object@basis.levels <- modelMatrices$basisLevels
      .Object@counts <- modelMatrices$counts
      # a count matrix as far as this constructor can tell; dbarts() states a
      # factor it one-hot expanded
      .Object@counts.as.factor <- if (is.null(modelMatrices$counts)) {
        NA
      } else {
        FALSE
      }
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
  if (dataKeepsFrame(data)) {
    coding <- data@build$coding
    return(validateXTest(newdata, coding$x, missing.seen = coding$missing.seen))
  }
  validateXTest(newdata, data@x, missing.seen = dataMissingSeen(data))
}

## A split rule holds a route for NA only on a column that can hold one: its
## training values carry one or, on a sampler, it has been given one since
## ('missing.seen', the data object's record). On a column that never has,
## every rule sends NA down one fixed branch. Refuse rather than answer from a
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

## The predictor columns a test set leaves without a route: missing there,
## complete in training and never seen missing.
unroutableTestColumns <- function(x.test, x.train, missing.seen = NULL) {
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
    if (
      isTRUE(missing.seen[j]) ||
        sourceColumnHasNA(x.train, j, numColumns, numTrain)
    ) {
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

unroutableTestRows <- function(x.test, x.train, missing.seen = NULL) {
  testRowsMissingIn(
    x.test,
    unroutableTestColumns(x.test, x.train, missing.seen)
  )
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
refuseTestMissingness <- function(
  x.test,
  x.train,
  naActionHint = FALSE,
  missing.seen = NULL
) {
  offending <- unroutableTestColumns(x.test, x.train, missing.seen)
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
  incomplete <- if (length(dim(value)) < 2L) {
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
refuseExtraMissingness <- function(extra, naActionHint, labels = NULL) {
  for (name in names(extra)) {
    if (any(extra[[name]])) {
      stop(
        extraMissingSubject(name, labels),
        " missing values",
        if (naActionHint) {
          paste0(
            "; use na.action = na.pass to return NA for them, or na.omit ",
            "to drop them"
          )
        },
        labels[[name]]$remedy
      )
    }
  }
  invisible(NULL)
}

## How a refusal names a missing 'extra' channel: "'name' has", or the
## subject and verb 'labels' gives it, with the remedy that closes it (a
## binomial's response columns, dec-B641).
extraMissingSubject <- function(name, labels = NULL) {
  subject <- labels[[name]]$subject
  if (is.null(subject)) paste0("'", name, "' has") else subject
}

## Joins quoted names in prose: one alone, or every name but the last
## comma-separated with "and" before it.
quotedNameList <- function(names) {
  joinNames(paste0("'", names, "'"))
}

joinNames <- function(quoted) {
  if (length(quoted) <= 1L) {
    return(quoted)
  }
  paste0(
    paste(quoted[-length(quoted)], collapse = ", "),
    " and ",
    quoted[length(quoted)]
  )
}

## The same, for a list that can be long: the first 'max' names, then a count
## of the rest. A message that names every level of a factor with thousands
## of them runs past R's limit on one and is cut mid-list. 'comma' keeps a
## list that fits as a plain comma-separated one.
quotedLevelList <- function(names, max = 5L, width = 40L, comma = FALSE) {
  names <- unique(names)
  shown <- names[seq_len(min(length(names), max))]
  # a label past 'width' characters is cut, so five of them cannot run the
  # message past the limit either; a missing label is shown as R prints one
  size <- nchar(shown, "chars", allowNA = TRUE)
  long <- !is.na(size) & size > width
  shown[long] <- paste0(substr(shown[long], 1L, width - 3L), "...")
  quoted <- ifelse(is.na(shown), "<NA>", paste0("'", shown, "'"))
  if (length(names) <= max) {
    return(if (comma) paste(quoted, collapse = ", ") else joinNames(quoted))
  }
  paste0(
    paste(quoted, collapse = ", "),
    " and ",
    length(names) - max,
    " more"
  )
}

## na.fail's refusal, naming whichever of the predictor columns and the
## 'extra' channels (dec-A89) carried a missing value.
refuseNaFail <- function(x.test, x.train, extra, labels = NULL) {
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
  labeled <- channels[channels %in% names(labels)]
  channels <- setdiff(channels, labeled)
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
  for (name in labeled) {
    parts <- c(
      parts,
      paste0(extraMissingSubject(name, labels), " missing values")
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
## default's refusal; 'labels' gives a channel its own words there
## (extraMissingSubject). Returns the logical 'keep' and whether the dropped
## rows 'pad' back as NA, or NULL when every row is answered as it stands.
resolvePredictRows <- function(
  na.action,
  x.test,
  x.train,
  extra = list(),
  missing.seen = NULL,
  labels = NULL
) {
  numTest <- NROW(x.test)
  extraDropped <- Reduce(`|`, extra, rep_len(FALSE, numTest))
  if (identical(na.action, stats::na.pass)) {
    dropped <- unroutableTestRows(x.test, x.train, missing.seen) | extraDropped
    if (!any(dropped) && numTest > 0L) {
      return(NULL)
    }
    return(list(keep = !dropped, pad = TRUE))
  }
  hasNA <- sourceHasNA(x.test) || any(extraDropped)
  if (identical(na.action, stats::na.fail)) {
    if (hasNA) {
      refuseNaFail(x.test, x.train, extra, labels)
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
    refuseTestMissingness(x.test, x.train, naActionHint = TRUE, missing.seen)
    refuseExtraMissingness(extra, naActionHint = TRUE, labels)
    return(NULL)
  }
  refuseTestMissingness(
    x.test[keep, , drop = FALSE],
    x.train,
    naActionHint = TRUE,
    missing.seen
  )
  refuseExtraMissingness(
    lapply(extra, `[`, keep),
    naActionHint = TRUE,
    labels
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
preparePredictRows <- function(
  newdata,
  x.train,
  na.action,
  channels = NULL,
  missing.seen = NULL,
  labels = NULL
) {
  if (missing(newdata) || is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
  rowNames <- observationRowNames(newdata)
  refuseMatrixOfCodes(newdata, x.train, "'newdata'")
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
    extra,
    missing.seen,
    labels
  )
  if (is.null(resolved)) {
    return(list(x = x, newdata = newdata, n = NROW(x), keptNames = rowNames))
  }
  result <- keptRowsRecord(resolved, rowNames)
  if (result$placeholder) {
    result$x <- suppressTestMatchWarnings(validateXTest(
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

suppressTestMatchWarnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) invokeRestart("muffleWarning")
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
  if (length(dim(value)) < 2L) {
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
  # an "s" attribute shares the value's layout: draws, or a band per row
  if (!is.null(s)) {
    attr(x, "s") <- padPredictedRows(s, rows, first = first)
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

## The names an expression reads as data: its symbols, less the function in
## each call, the pieces of a pkg::name reference and what follows a '$' or
## '@'.
expressionDataNames <- function(expr) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    return(if (nzchar(name)) name)
  }
  if (!is.call(expr)) {
    return(NULL)
  }
  head <- expr[[1L]]
  if (is.symbol(head) && as.character(head) %in% c("::", ":::", "function")) {
    return(NULL)
  }
  arguments <- as.list(expr)[-1L]
  if (is.symbol(head) && as.character(head) %in% c("$", "@")) {
    arguments <- arguments[1L]
  }
  unique(c(
    if (is.call(head)) expressionDataNames(head),
    unlist(lapply(arguments, expressionDataNames), use.names = FALSE)
  ))
}

## The function a call names, as a stored expression resolves it with no
## environment: a base function by its own name, a package function as
## pkg::name (or the function itself, when it sits in its namespace under
## another name). A function from no package is refused by name, since nothing
## would carry it to another session.
qualifyFunctionName <- function(head, env) {
  name <- as.character(head)
  fn <- get0(name, envir = env, mode = "function")
  if (is.null(fn) || is.primitive(fn)) {
    return(head)
  }
  fnEnv <- environment(fn)
  if (identical(fnEnv, baseenv()) || identical(fnEnv, .BaseNamespaceEnv)) {
    return(head)
  }
  if (!isNamespace(fnEnv)) {
    stop(
      "the formula calls '",
      name,
      "', a function from no package; a fit stores its formula without ",
      "the environment that function lives in, so compute the column in ",
      "'data' instead"
    )
  }
  package <- getNamespaceName(fnEnv)
  if (identical(package, "base")) {
    return(head)
  }
  if (!identical(get0(name, envir = fnEnv, inherits = FALSE), fn)) {
    return(fn)
  }
  call(
    if (name %in% getNamespaceExports(fnEnv)) "::" else ":::",
    as.name(package),
    as.name(name)
  )
}

## An expression stored on a fit, made independent of the environment it was
## written in: a data column stays a name, read from new rows alone; every
## other name is replaced by its fit-time value and every function call is
## qualified (qualifyFunctionName). A name is a data column when 'data'
## carries it ('dataNames'), or, when there is no 'data' (NULL), when its
## value has one row per observation ('numRows').
freezeExpression <- function(expr, dataNames, env, numRows) {
  freeze <- function(e) {
    if (is.symbol(e)) {
      name <- as.character(e)
      if (!nzchar(name) || name %in% dataNames) {
        return(e)
      }
      value <- get0(name, envir = env, inherits = TRUE)
      if (is.null(value) && !exists(name, envir = env)) {
        return(e)
      }
      if (is.function(value)) {
        return(qualifyFunctionName(e, env))
      }
      # without 'data', a value with one row per observation is a column
      if (
        is.null(dataNames) &&
          (!is.null(dim(value)) ||
            length(value) > 1L && (is.na(numRows) || NROW(value) == numRows))
      ) {
        return(e)
      }
      return(value)
    }
    if (!is.call(e)) {
      return(e)
    }
    head <- e[[1L]]
    if (is.symbol(head) && as.character(head) %in% c("::", ":::", "function")) {
      return(e)
    }
    positions <- seq_along(e)[-1L]
    if (is.symbol(head) && as.character(head) %in% c("$", "@")) {
      positions <- 2L
    }
    if (is.symbol(head)) {
      e[[1L]] <- qualifyFunctionName(head, env)
    } else {
      e[[1L]] <- freeze(head)
    }
    for (i in positions) {
      if (!identical(e[[i]], quote(expr = ))) {
        e[i] <- list(freeze(e[[i]]))
      }
    }
    e
  }
  freeze(expr)
}

## A formula fit's predictor terms, kept on the design matrix as its "terms"
## attribute: the training terms without the response, carrying the predvars
## model.frame stamped on them, so that a data-dependent basis (poly(), ns(),
## scale()) is rebuilt on the test set and on predict's newdata from the
## training values, as predict.lm rebuilds it, and the offset() terms can be
## evaluated there. The predvars are frozen (freezeExpression) and the terms
## keep no environment but base, so a fit carries neither its caller's frame
## nor the global environment and predicts from newdata alone; the variables
## stay as written, since they name the columns.
predictorTerms <- function(modelTerms, dataNames, numRows) {
  terms <- stats::delete.response(modelTerms)
  env <- environment(modelTerms)
  predvars <- attr(terms, "predvars")
  if (is.null(predvars)) {
    predvars <- attr(terms, "variables")
  }
  attr(terms, "predvars") <- freezeExpression(
    predvars,
    dataNames,
    env,
    numRows
  )
  environment(terms) <- baseenv()
  terms
}

## The names a formula's 'data' carries, and the rows it has (NA when that is
## not known), against which freezeExpression tells data columns apart.
formulaDataNames <- function(data) {
  if (is.environment(data)) {
    return(ls(data, all.names = TRUE))
  }
  names(data)
}

## The same terms without their offset() terms, for rebuilding the predictors
## alone; model.frame reads only the variables and predvars.
dropOffsetTerms <- function(terms) {
  offsets <- attr(terms, "offset")
  if (is.null(offsets)) {
    return(terms)
  }
  keep <- -(offsets + 1L)
  attr(terms, "variables") <- attr(terms, "variables")[keep]
  if (!is.null(attr(terms, "predvars"))) {
    attr(terms, "predvars") <- attr(terms, "predvars")[keep]
  }
  if (!is.null(attr(terms, "dataClasses"))) {
    attr(terms, "dataClasses") <- attr(terms, "dataClasses")[-offsets]
  }
  attr(terms, "offset") <- NULL
  terms
}

## Names the variables an expression list needs that 'data' does not carry,
## rather than letting model.frame or eval find a same-named object in an
## enclosing scope (a predictor 'c' binding to base::c).
refuseMissingVariables <- function(neededVars, data, argument, what) {
  missingVars <- neededVars[neededVars %not_in% names(data)]
  if (length(missingVars) > 0L) {
    stop(
      "'",
      argument,
      "' data is missing ",
      if (length(missingVars) > 1L) "variables" else "variable",
      " required by ",
      what,
      ": '",
      toString(missingVars),
      "'"
    )
  }
  invisible(NULL)
}

## The fit's offset() terms evaluated on 'newdata', as predict.lm evaluates
## them, or NULL when the fit's formula has none.
formulaTermOffset <- function(x.train, newdata, argument) {
  terms <- attr(x.train, "terms")
  offsets <- attr(terms, "offset")
  if (is.null(offsets)) {
    return(NULL)
  }
  if (!is.data.frame(newdata) && !(is.list(newdata) && !is.object(newdata))) {
    stop(
      "the fit's formula has an offset() term, which is evaluated on '",
      argument,
      "'; supply '",
      argument,
      "' as a data frame"
    )
  }
  # the frozen expressions: data columns come from the new rows alone
  predvars <- attr(terms, "predvars")
  expressions <- lapply(offsets, function(i) predvars[[i + 1L]])
  refuseMissingVariables(
    unique(unlist(lapply(expressions, expressionDataNames))),
    newdata,
    argument,
    "the formula's offset() term"
  )
  n <- if (is.data.frame(newdata)) nrow(newdata) else NROW(newdata[[1L]])
  result <- numeric(n)
  for (expression in expressions) {
    value <- as.double(eval(expression, newdata, environment(terms)))
    if (length(value) != n && length(value) != 1L) {
      stop(
        "the formula's offset() term '",
        deparse1(expression),
        "' has ",
        length(value),
        " values for the ",
        n,
        " rows of '",
        argument,
        "'"
      )
    }
    result <- result + value
  }
  result
}

## The 'offset' argument as written, kept on the data object as its
## "offset.argument" attribute (not on the design matrix, which is often the
## caller's own and would be copied): a one-sided formula carrying the
## expression and
## the environment it is evaluated in, so that predict and the test-set default
## re-evaluate it on new rows, as predict.lm re-evaluates lm's call$offset. An
## argument forwarded through a wrapper's dots is recovered as written.
offsetArgumentFormula <- function(expr, env, dataNames, numRows) {
  if (isDotsReference(expr)) {
    written <- recoverForwardedArgument(expr, env)
    expr <- written$expr
    env <- written$env
  }
  # a function from no package is refused, as in a formula
  frozen <- tryCatch(
    freezeExpression(expr, dataNames, env, numRows),
    error = function(e) {
      stop(
        sub(
          "^the formula calls",
          "the 'offset' expression calls",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  result <- stats::as.formula(call("~", frozen), env = baseenv())
  # the expression as written, for messages
  attr(result, "written") <- if (is.language(expr)) deparse1(expr)
  result
}

## The 'weights' argument of a formula fit as written when it names columns of
## the fit's data, kept on the data object as its "weights.argument" attribute:
## a one-sided formula whose expression reads those columns from new rows alone
## (freezeExpression), so that predict finds them in newdata and never in the
## caller's frame. NULL where the weights name no column: a vector, a constant,
## or a fit with no data to name columns.
weightsArgumentFormula <- function(expr, env, dataNames) {
  # d$w reads a vector from an object, not a column of the data
  if (
    is.null(expr) ||
      is.null(dataNames) ||
      isDotsReference(expr) ||
      any(c("$", "@") %in% all.names(expr))
  ) {
    return(NULL)
  }
  frozen <- tryCatch(
    freezeExpression(expr, dataNames, env, NA_integer_),
    error = function(e) NULL
  )
  if (is.null(frozen) || length(all.vars(frozen)) == 0L) {
    return(NULL)
  }
  stats::as.formula(call("~", frozen), env = baseenv())
}

## A formula fit's 'offset' expression evaluated as model.frame evaluates
## one: in 'data', with 'env' for every name it does not carry. Returns the
## error rather than signaling it.
evaluateOffsetExpression <- function(expr, data, env) {
  tryCatch(
    if (is.environment(data)) {
      eval(expr, data)
    } else {
      eval(expr, data, env)
    },
    error = function(e) e
  )
}

## How an offset argument reads in a message: its expression, or a plain
## vector by its kind.
describeOffsetArgument <- function(argument) {
  text <- attr(argument, "written")
  if (is.null(text)) {
    return("a numeric vector")
  }
  paste0(
    "'",
    if (nchar(text) > 60L) paste0(substr(text, 1L, 57L), "...") else text,
    "'"
  )
}

## The fit's 'offset' argument evaluated on new rows: NULL when the fit had
## none, FALSE when it cannot be (its names found neither in the new rows nor
## in its environment, or a value of the wrong length - a plain vector given
## for the training rows), and otherwise one value or one per row.
evaluateOffsetArgument <- function(argument, newdata) {
  if (is.null(argument)) {
    return(NULL)
  }
  data <- if (
    is.data.frame(newdata) || (is.list(newdata) && !is.object(newdata))
  ) {
    newdata
  } else if (is.matrix(newdata) && !is.null(colnames(newdata))) {
    as.data.frame(newdata)
  }
  n <- NROW(newdata)
  value <- tryCatch(
    as.double(eval(argument[[2L]], data, environment(argument))),
    error = function(e) NULL
  )
  if (is.null(value) || (length(value) != 1L && length(value) != n)) {
    return(FALSE)
  }
  value
}

## Two shares of an offset at new rows summed, each a single value or one per
## row; either alone when the other is absent.
addOffsetShares <- function(share, offset, argument, rows) {
  if (is.null(share) || is.null(offset)) {
    return(if (is.null(share)) offset else share)
  }
  if (
    !is.null(dim(offset)) ||
      (length(offset) != 1L &&
        length(share) != 1L &&
        length(offset) != length(share))
  ) {
    stop("'", argument, "' must have the same number of rows as '", rows, "'")
  }
  share + offset
}

## The offset a formula fit's test set or predict's newdata carries: its
## offset() terms evaluated there plus the offset the caller gives for those
## rows, a single number or one per row; either alone when the other is
## absent.
addFormulaTermOffset <- function(x.train, newdata, offset, argument, rows) {
  termOffset <- formulaTermOffset(x.train, newdata, rows)
  if (is.null(termOffset) || is.null(offset)) {
    return(if (is.null(termOffset)) offset else termOffset)
  }
  if (
    !is.null(dim(offset)) ||
      (length(offset) != 1L && length(offset) != length(termOffset))
  ) {
    stop(
      "'",
      argument,
      "' must have the same number of rows as '",
      rows,
      "'"
    )
  }
  termOffset + offset
}

## Predictors given as a matrix or vector take the types a data frame column
## does: integer and logical values are numbers, any sparse Matrix class (a
## sparseVector as one column) is the dgCMatrix the engine ingests, and a
## dense Matrix class is a plain matrix. A factor is left as is.
asNumericPredictors <- function(x) {
  if (isS4(x) && methods::is(x, "sparseVector")) {
    x <- methods::as(x, "CsparseMatrix")
  } else if (
    isS4(x) && methods::is(x, "Matrix") && !methods::is(x, "sparseMatrix")
  ) {
    x <- as.matrix(x)
  }
  x <- asDgCMatrix(x)
  if (
    !is.factor(x) &&
      (is.matrix(x) || is.null(dim(x))) &&
      (is.integer(x) || is.logical(x))
  ) {
    # storage.mode<- keeps the dimnames that as.double() would drop
    storage.mode(x) <- "double"
  }
  x
}

## The data columns a fit's predictors are read from: its stored terms'
## variables, else its term labels' and column names.
predictorDataNames <- function(x.train) {
  trainTerms <- attr(x.train, "terms")
  if (!is.null(trainTerms)) {
    return(all.vars(attr(dropOffsetTerms(trainTerms), "variables")))
  }
  labels <- sub("^`(.*)`$", "\\1", attr(x.train, "term.labels"))
  unique(c(labels, colnames(x.train)))
}

## A sparse column of a data-frame test set made dense where the fit cannot
## take it resident: on a formula fit whose own column of that name is dense,
## the model frame replay reads its values (a data-dependent term such as
## poly() included), and on a fit with no level table the indicator expansion
## does. A sparseFactor becomes the factor over its levels, a sparseVector or
## a one-column sparse matrix a vector, a wider sparse matrix a matrix. A
## column the training design itself holds sparse, and every sparse column of
## a categorical fit made from a data frame, stays as given.
densifyUnheldSparseColumns <- function(x.test, x.train) {
  isSparseColumn <- vapply(x.test, isSparseDataFrameColumn, FALSE)
  if (
    !any(isSparseColumn) ||
      (is.null(attr(x.train, "terms")) &&
        !is.null(attr(x.train, "factor.levels")))
  ) {
    return(x.test)
  }
  heldSparse <- if (predictorSourceIsSparse(x.train)) {
    if (inherits(x.train, "dgCMatrix")) {
      colnames(x.train)
    } else {
      x.train$columnNames[x.train$map < 0L]
    }
  } else {
    character()
  }
  for (name in names(x.test)[isSparseColumn]) {
    if (
      name %in% heldSparse || any(startsWith(heldSparse, paste0(name, ".")))
    ) {
      next
    }
    column <- x.test[[name]]
    x.test[[name]] <- if (methods::is(column, "sparseFactor")) {
      sparseFactorToFactor(column)
    } else if (methods::is(column, "sparseVector")) {
      as.vector(column)
    } else if (ncol(column) == 1L) {
      as.vector(as.matrix(column))
    } else {
      as.matrix(column)
    }
  }
  x.test
}

validateXTest <- function(
  x.test,
  x.train,
  refuseMissing = TRUE,
  missing.seen = NULL
) {
  termLabels <- attr(x.train, "term.labels")
  numPredictors <- ncol(x.train)
  predictorNames <- colnames(x.train)
  drop <- attr(x.train, "drop")
  factorLevels <- attr(x.train, "factor.levels")
  varTypes <- attr(x.train, "varTypes")

  if (is.null(x.test)) {
    return(x.test)
  }
  # a plain list of columns is the data frame it holds, as predict.lm reads
  # one
  if (is.list(x.test) && !is.object(x.test)) {
    if (length(x.test) == 0L || length(unique(lengths(x.test))) != 1L) {
      stop("'test' must be a data frame, a matrix or a vector")
    }
    x.test <- as.data.frame(x.test, optional = TRUE)
  }
  x.test <- asNumericPredictors(x.test)
  # a vector is one row, its names those of the columns
  if (is.numeric(x.test) && is.null(dim(x.test)) && length(x.test) > 0L) {
    x.test <- matrix(
      x.test,
      nrow = 1L,
      dimnames = list(NULL, names(x.test))
    )
  }
  if (is.numeric(x.test) && NCOL(x.test) == 0L) {
    return(NULL)
  }
  testFactorLevels <- NULL
  if (is.data.frame(x.test)) {
    # a fit made from a matrix is given a frame as a frame fit is: every
    # training name among its columns is taken by name and the rest are
    # ignored; a missing name falls through to the refusals below
    if (
      is.null(termLabels) &&
        !is.null(predictorNames) &&
        anyDuplicated(predictorNames) == 0L &&
        all(predictorNames %in% names(x.test))
    ) {
      x.test <- x.test[predictorNames]
    }
    refusePOSIXltColumns(x.test, predictorDataNames(x.train))
    x.test <- densifyUnheldSparseColumns(x.test, x.train)
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
    # the term replay: names an absent variable up front and rebuilds the
    # given term labels on the test rows. A formula fit's training terms
    # rebuild them with the training predvars, as predict.lm does, so a
    # data-dependent basis means at a test row what it meant in training,
    # whatever other rows the test set holds; a label they do not cover (a
    # hazard fit's appended period) is read by its own name
    trainTerms <- attr(x.train, "terms")
    replayTerms <- function(data, labels) {
      replayLabels <- labels
      frame <- NULL
      if (!is.null(trainTerms)) {
        predictorOnly <- dropOffsetTerms(trainTerms)
        refuseMissingVariables(
          expressionDataNames(attr(predictorOnly, "predvars")),
          data,
          "test",
          "the model"
        )
        frame <- model.frame(
          predictorOnly,
          data = data,
          na.action = stats::na.pass
        )
        replayLabels <- labels[labels %not_in% names(frame)]
        frame <- frame[labels[labels %in% names(frame)]]
      }
      if (length(replayLabels) > 0L) {
        testFormula <- as.formula(paste(
          "~",
          paste(replayLabels, collapse = " + ")
        ))
        refuseMissingVariables(all.vars(testFormula), data, "test", "the model")
        replayed <- model.frame(
          formula = testFormula,
          data = data,
          na.action = stats::na.pass
        )
        if (is.null(frame)) {
          return(replayed)
        }
        for (label in replayLabels) {
          frame[[label]] <- replayed[[label]]
        }
        frame <- frame[labels]
      }
      frame
    }
    if (any(isSparseColumn)) {
      # sparse columns ride to the engine unexpanded, coded over the training
      # level table; the resulting container is preserved below. The model
      # frame replay takes no S4 columns, so they are lifted out, the dense
      # remainder replayed, and the columns the model uses re-attached.
      if (is.null(factorLevels)) {
        stop(
          "a data-frame 'test' with sparse columns is read against a fit ",
          "made from a data frame or a formula; give 'test' as a dgCMatrix, ",
          "or its columns dense"
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
        indicatorLevels <- attr(x.train, "indicator.levels")
        if (!is.null(indicatorLevels)) {
          x.test <- mapFactorColumnsToIndicatorLevels(
            x.test,
            indicatorLevels,
            drop
          )
        }
        if (is.list(drop) && length(drop) == length(x.test)) {
          refuseWiderTestColumns(x.test, drop)
        }
        # a fit with no drop record was made from a matrix, which drops no
        # column, so no column of the test frame is dropped either: one
        # constant on the test rows, as every column of a one-row frame is,
        # is a value to predict at
        x.test <- makeIndicatorModelMatrix(
          x.test,
          if (!is.null(drop)) drop else FALSE,
          storage = "dense"
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
    if (!is.numeric(x.test) && !is.logical(x.test)) {
      stop("test matrix must be numeric")
    }

    if (is.integer(x.test) || is.logical(x.test)) {
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
          quotedLevelList(testFactorLevels[[mismatched[1L]]]),
          "); give the test factor the training levels, factor(x, levels = ",
          "...), or fit with factors = \"categorical\""
        )
      }
    }
    stop("number of columns in 'test' must be equal to that of 'x'")
  }
  xIsNamed <- !is.null(predictorNames)
  testIsNamed <- !is.null(colnames(x.test))
  repeatedNames <- anyDuplicated(predictorNames) > 0L
  # both named, the fit's names distinct: 'test' is matched by name at any
  # predictor count, and a column that is not the fit's is refused, as at
  # setData and setPredictor
  if (xIsNamed && testIsNamed && !repeatedNames) {
    matchIndices <- match(predictorNames, colnames(x.test))
    if (any(is.na(matchIndices))) {
      stop(
        "column names of 'test' do not match those of 'x': '",
        toString(predictorNames[is.na(matchIndices)]),
        "' present in 'x' but not in 'test' (whose columns are '",
        toString(colnames(x.test)),
        "')",
        formulaPrefixHint(
          predictorNames[is.na(matchIndices)],
          colnames(x.test)
        )
      )
    }
  }
  if (numPredictors > 1) {
    columnIndices <- seq.int(numPredictors)
    if (repeatedNames) {
      warning(
        "'x' has repeated column names, so 'test' is matched by position",
        call. = FALSE
      )
    } else if (!xIsNamed && testIsNamed) {
      warning(
        "'x' and 'test' are not both named; columns of 'test' will be matched by position",
        call. = FALSE
      )
    } else if (xIsNamed && testIsNamed) {
      columnIndices <- matchIndices
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
    refuseTestMissingness(x.test, x.train, missing.seen = missing.seen)
  }

  x.test
}

## A formula fit's 'weights' or 'test' argument as model.frame reads its
## extras: the expression evaluated in 'data' with the formula's environment
## for every name 'data' does not carry, so a name is never looked up in this
## package's own frames, whatever the caller has defined globally. NULL when
## it cannot be evaluated there, which leaves the caller its plain argument.
findTermInFormulaData <- function(formula, data, term) {
  expr <- substitute(term)
  if (is.numeric(expr)) {
    return(term)
  }
  # any other constant is model.frame's to judge
  if (!is.language(expr)) {
    return(NULL)
  }
  enclos <- environment(formula)
  value <- tryCatch(
    if (missing(data) || is.null(data)) {
      eval(expr, enclos)
    } else if (is.environment(data)) {
      eval(expr, data)
    } else {
      eval(expr, data, enclos)
    },
    error = function(e) NULL
  )
  # a name that reached a function (base's t, say) found no value of the
  # caller's
  if (is.function(value)) NULL else value
}

## The nearest binding of 'name' along env's chain that holds what an offset
## can be - a number, a list or NULL - as list(value =), and NULL when there is
## none. get0() cannot serve: it reports a name bound to NULL as one that is
## not bound at all.
findOffsetBinding <- function(name, env) {
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
      if (mode(value) %in% c("numeric", "list", "NULL")) {
        return(list(value = value))
      }
    }
    env <- parent.env(env)
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

  # 'offset' named in offset.test is the argument: beside a formula's offset()
  # terms, which the training offset already includes and which are evaluated
  # on 'test' separately, it is the argument's own share - evaluated on 'test'
  # as the default is, or failing that its training value - and 0 when not
  # given
  argumentShare <- if (hasOffsetTerm) {
    testShare <- evaluateOffsetArgument(offsetArgument, test)
    if (is.null(argumentOffset)) {
      0
    } else if (!isFALSE(testShare)) {
      testShare
    } else if (isTRUE(offsetGivenAsScalar)) {
      argumentOffset[1L]
    } else {
      argumentOffset
    }
  }

  if (is.symbol(matchedCall$offset.test)) {
    testOffsetName <- as.character(matchedCall$offset.test)

    if (identical(testOffsetName, "offset") && hasOffsetTerm) {
      return(list(offset.test = argumentShare, testUsesRegularOffset = FALSE))
    }
    if (identical(testOffsetName, "offset") && !is.null(offset)) {
      return(list(
        offset.test = if (isTRUE(offsetGivenAsScalar)) offset[1] else offset,
        testUsesRegularOffset = TRUE
      ))
    }

    # a test-set column comes first: the test offset belongs to the test
    # rows, as predict.lm evaluates an offset in newdata
    if (
      !testIsMissing &&
        (is.data.frame(test) || (is.list(test) && !is.object(test))) &&
        any(names(test) == testOffsetName)
    ) {
      return(list(
        offset.test = test[[testOffsetName]],
        testUsesRegularOffset = FALSE
      ))
    }

    if (
      is.formula(formula) &&
        !dataIsMissing &&
        any(names(data) == testOffsetName)
    ) {
      return(list(
        offset.test = data[[testOffsetName]],
        testUsesRegularOffset = FALSE
      ))
    }
    # then the enclosure model.frame reads it in, and failing that the
    # caller's frame; never this function's own, which would answer a name
    # such as 'x' with the predictors. An offset is a value, so a function
    # met first along either chain - base's t, stats' df or sigma - is passed
    # over for the value it masks, and a name bound to NULL is a NULL offset,
    # as a literal NULL is
    for (env in unique(list(testOffsetEnclos, callerEnv))) {
      found <- findOffsetBinding(testOffsetName, env)
      if (!is.null(found)) {
        return(list(offset.test = found$value, testUsesRegularOffset = FALSE))
      }
    }

    stop("cannot find test offset '", testOffsetName, "'")
  } else if (is.language(matchedCall$offset.test)) {
    ## offset.test could have been something like (offset + 0.5), or (offset + variable)
    baseOffset <- if (hasOffsetTerm) {
      argumentShare
    } else if (is.null(offset)) {
      NA_real_
    } else {
      if (isTRUE(offsetGivenAsScalar)) offset[1] else offset
    }

    if (identical(matchedCall$offset.test, quote(offset))) {
      return(list(
        offset.test = baseOffset,
        testUsesRegularOffset = !hasOffsetTerm
      ))
    }

    testOffset <- subTermInLanguage(
      matchedCall$offset.test,
      quote(offset),
      baseOffset
    )

    # in 'test', then 'data', each with the enclosure model.frame reads an
    # expression in for every name it does not carry; then that enclosure and
    # the caller's frame alone
    evaluateTestOffset <- function(where) {
      tryCatch(
        if (is.null(where)) {
          eval(testOffset, testOffsetEnclos)
        } else if (is.environment(where)) {
          eval(testOffset, where)
        } else {
          eval(testOffset, where, testOffsetEnclos)
        },
        error = function(e) e
      )
    }
    if (
      !testIsMissing &&
        (is.data.frame(test) || (is.list(test) && !is.object(test))) &&
        any(all.vars(testOffset) %in% names(test))
    ) {
      tryResult <- evaluateTestOffset(test)
      if (!inherits(tryResult, "error")) {
        return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
    }
    if (is.formula(formula) && !dataIsMissing) {
      tryResult <- evaluateTestOffset(data)
      if (!inherits(tryResult, "error")) {
        return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
    }
    tryResult <- evaluateTestOffset(NULL)
    if (!inherits(tryResult, "error")) {
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }
    if (!identical(testOffsetEnclos, callerEnv)) {
      tryResult <- tryCatch(eval(testOffset, callerEnv), error = function(e) e)
      if (!inherits(tryResult, "error")) {
        return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
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
## The attribute dbarts() sets on a formula it hands the data object under a
## logistic family.
binomialResponseAttr <- "dbartsBinomialResponse"

## A data frame as the matrix interface's response, which no family reads:
## the counts a binomial or multinomial response gives are a matrix, as glm
## takes them.
refuseDataFrameResponse <- function(y) {
  if (is.data.frame(y)) {
    stop(
      "'y' is a data frame; give the response as a vector, or counts such ",
      "as a binomial's successes and failures as a matrix, ",
      "cbind(successes, failures) or as.matrix(y)",
      call. = FALSE
    )
  }
  invisible(NULL)
}

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
        "response - successes and failures and per-category counts go to ",
        "dbartsData(counts = ), a (time, status) pair to dbarts()/bart() ",
        "with family = \"aft\"/\"hazard\""
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
  text <- if (
    family == "probit" && K == 2L && entry %in% c("bart()/dbarts()", "xbart()")
  ) {
    probitBinomialMessage
  } else if (identical(entry, "xbart()")) {
    paste0(
      "xbart() takes a single-column response or, under family = binomial, ",
      "cbind(successes, failures); 'y' is an n x ",
      K,
      " matrix, and it cross-validates neither multinomial counts nor a ",
      "survival response; use bart() instead - ",
      readings
    )
  } else if (entry %in% c("bartBT()", "rbart_vi()")) {
    paste0(
      entry,
      " takes a single-column response but 'y' is an n x ",
      K,
      " matrix, and it ",
      if (identical(entry, "xbart()")) "cross-validates" else "fits",
      " neither multinomial counts nor a survival response; use bart() ",
      "instead - ",
      readings
    )
  } else if (
    family %in% c("aft", "hazard", "hazard.probit", "hazard.logistic")
  ) {
    # the caller already named the survival family; what is wrong is the
    # spelling of the response, which a formula writes as Surv()
    paste0(
      "family = \"",
      family,
      "\" takes a formula response written survival::Surv(time, status); ",
      "'y' is an n x ",
      K,
      " matrix, and cbind() on the left of a formula is not read as a ",
      "time and a status"
    )
  } else if (family == "logistic") {
    binomialCountColumnsMessage(K, "the response is")
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
      "successes and failures need family = binomial; ",
      readings
    )
  } else if (isAutoCountMatrix(y)) {
    paste0(
      "an n x ",
      K,
      " count matrix as a formula response is fit by bart(cbind(c1, c2, ",
      "c3) ~ x, data), which needs 'data'; dbarts() takes counts on the ",
      "matrix interface, dbarts(x, counts, family = \"multinomial\")"
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
  # classed, and carrying the response, so a door that wraps this one can
  # restate it once more in its own words
  stop(structure(
    class = c("dbartsMatrixResponseRestated", "error", "condition"),
    list(message = text, call = NULL, y = y)
  ))
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
# message does. `shownFamily` is the family those two messages name: the
# caller's own token where the entry point has remapped it ("student" fits
# as "gaussian").
resolveClassificationFamily <- function(
  data,
  family,
  caller,
  incompatibleFamilies,
  splitMultinomialMessage = FALSE,
  allowOrdinal = FALSE,
  shownFamily = family
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
        if (isOrdered) {
          paste0(" - ", caller, " fit only binary and continuous responses")
        } else {
          paste0(
            " or dbarts(x, y, family = \"multinomial\") - bartBT() and ",
            "xbart() fit only binary and continuous responses"
          )
        }
      )
    }
    stop(
      "family \"",
      shownFamily,
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
      shownFamily,
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
    stop(
      "family = \"ordinal\" requires a response with at least 2 categories; ",
      "to declare categories that no row takes, give the response as ",
      "factor(y, levels = ...)"
    )
  }
  list(y = as.double(codes), K = as.integer(K), levels = labels)
}

# Resolve the negative-binomial shape argument to the length-1 real the C
# bridge reads off the model's shape slot. NA (the
# default) estimates r on the capped integer grid, encoded as a non-positive
# spec; a supplied value FIXES r as a single positive integer - the exact
# integer envelope v1 ships, real shape not yet supported.
resolveShape <- function(shape) {
  if (length(shape) != 1L) {
    stop("'shape' must be a single value")
  }
  if (is.na(shape)) {
    return(-1) # a non-positive spec estimates r on the grid
  }
  if (!is.numeric(shape) || shape <= 0) {
    stop("'shape' must be a positive number")
  }
  if (shape != round(shape)) {
    stop(
      "family \"nbinom\" fits an integer shape; real shape is not ",
      "yet supported"
    )
  }
  as.double(shape)
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
# through; otherwise every non-null element is read by class (basisClass): one
# of a class a basis does not take is refused, numbers are coerced to a
# numeric matrix and a factor, a character vector or a logical vector stays as
# it is. Only the length is checked, against 'numRows', and 'subset', if
# given, is applied as it is; what the kept rows hold, a missing value among
# it, is judged where they are final (finishBasisRows). 'argument' names the
# surface the value came from and 'rows' what the row count is checked
# against.
validateForestBases <- function(
  bases,
  numRows,
  subset = NULL,
  argument = "bases",
  rows = "'y'"
) {
  if (is.null(bases)) {
    return(NULL)
  }
  if (!is.list(bases)) {
    stop("'", argument, "' must be a list of per-forest basis matrices")
  }
  labels <- attr(bases, "forest.labels", exact = TRUE)
  lapply(seq_along(bases), function(i) {
    basis <- readBasisValue(bases[[i]], argument, forestName(i, labels))
    if (is.null(basis)) {
      return(NULL)
    }
    if (NROW(basis) != numRows) {
      stop("'", argument, "' must have the same length as ", rows)
    }
    if (!is.null(subset)) {
      basis <- cutBasisRows(basis, subset)
    }
    basis
  })
}

# A basis handed over as a value, read by class (basisClass): one of a class a
# basis does not take is refused, naming 'forest', and numbers become a
# numeric matrix.
readBasisValue <- function(basis, argument, forest) {
  if (is.null(basis)) {
    return(NULL)
  }
  read <- basisClass(basis)
  if (is.na(read$kind)) {
    refuseBasisClass(read, argument, forest, NCOL(basis))
  }
  if (read$kind == "numeric") {
    basis <- as.matrix(basis)
    storage.mode(basis) <- "double"
  }
  basis
}

# An integer means what the same number as a double does; a factor is not one.
integerToDouble <- function(x) {
  if (is.integer(x)) {
    storage.mode(x) <- "double"
  }
  x
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

# The refusal for a multinomial response that is neither labels nor counts,
# said by bart() and dbarts() alike: a bare number is not a category label.
multinomialResponseMessage <- paste0(
  "family = \"multinomial\" requires a factor (or character) response, ",
  "or an n x K count matrix"
)

# The shapes a multinomial response is written in, read alike by bart() and
# dbarts(): a data frame is the numeric matrix of its columns (one that is
# not numeric is refused, a column of labels included), and a one-column character matrix is its vector of
# labels, as a one-column matrix is read where a vector is wanted. Returns a
# numeric matrix, or a vector for the caller's label path.
readMultinomialShape <- function(y) {
  refuse <- function(type) {
    stop(
      "a multinomial count-matrix response must be numeric, not ",
      type,
      "; give labels as a vector, or as a one-column matrix of labels",
      call. = FALSE
    )
  }
  if (is.data.frame(y)) {
    if (!all(vapply(y, is.numeric, NA))) {
      refuse(typeof(as.matrix(y)))
    }
    y <- as.matrix(y)
  }
  if (is.matrix(y) && !is.numeric(y)) {
    if (is.character(y) && ncol(y) == 1L) {
      return(as.vector(y))
    }
    refuse(typeof(y))
  }
  y
}

# Coerce a raw multinomial response to the n x K count matrix the softmax
# engine's response is. Two shapes arrive: an
# n x K matrix (or data frame) of non-negative integer counts, whose trials are
# its row sums, and a length-n vector of labels - a factor or character - which
# is the single-trial special case, one-hot expanded with every trial 1. A
# number or a logical is refused: neither names a category, and reading whole
# numbers as codes would build an empty category for y in 1, 2, 3. The
# category labels ride the result's COLUMN NAMES, the carrier that survives
# both serialization and the engine's re-creation; the engine reads neither.
# A missing label one-hot-expands to an all-NA row rather than erroring: a
# single-trial spike has no cell to put a missing category in, so the row
# itself is what is missing, exactly as a count-matrix row is missing when
# any of its cells is NA (validateMultinomialCounts). na.action reads that
# row downstream, as it does for every other family's response.
resolveMultinomialCounts <- function(y) {
  y <- readMultinomialShape(y)
  if (is.matrix(y)) {
    return(y)
  }
  if (is.character(y)) {
    y <- factor(y)
  }
  if (!is.factor(y)) {
    stop(multinomialResponseMessage)
  }
  codes <- as.integer(y)
  levels <- levels(y)
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
# other caller ($setResponse mutates a live sampler's response outright, with
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
# rows with no trial. Called only after the creation or $setResponse .Call has
# returned, so a refused call never spends the key; NULL (a sampler that
# carries no counts) and a matrix without such a row are silent.
warnZeroTrials <- function(counts, binomial = FALSE) {
  if (is.null(counts)) {
    return(invisible(NULL))
  }
  numEmpty <- sum(rowSums(counts) == 0L, na.rm = TRUE)
  if (numEmpty > 0L) {
    warnOnce(
      "multinomialZeroTrials",
      sprintf(
        paste0(
          if (binomial) "binomial rows" else "multinomial count rows",
          " with zero trials (%d of %d) contribute nothing to the ",
          "likelihood and still receive fitted probabilities"
        ),
        numEmpty,
        nrow(counts)
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
  counts = NULL,
  n.cuts = NULL,
  useQuantiles = NULL,
  keep.data = TRUE
) {
  stated <- checkFlatCopies(n.cuts = n.cuts, useQuantiles = useQuantiles)
  # set by dbarts() on a formula under a logistic family, whose two-column
  # response is then a binomial's successes and failures
  binomialResponse <- isTRUE(attr(formula, binomialResponseAttr, exact = TRUE))
  formulaCounts <- NULL
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  offsetIsMissing <- missing(offset)
  testOffsetIsMissing <- missing(offset.test)
  matchedCall <- match.call()
  # where the caller's own expressions are evaluated when no formula supplies
  # an environment for them, as model.frame falls back; never this function's
  # frame, whose locals ('x', 'data', 'offset') would shadow the caller's names
  # nolint next: object_usage_linter. getTestOffset reads it, through eval.
  callerEnv <- parent.frame()
  # an offset.test forwarded through a wrapper's dots arrives as ..N: the
  # expression it was written as stands in its place, and the names in it are
  # read where it was written, as the training offset's are
  if (isDotsReference(matchedCall$offset.test)) {
    written <- recoverForwardedArgument(matchedCall$offset.test, callerEnv)
    matchedCall["offset.test"] <- list(written$expr)
    callerEnv <- written$env
  }
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
    makeIndicatorModelMatrix
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
  # a formula's offset() terms, and the 'offset' argument's own share of the
  # training offset beside them (set in the formula branch)
  hasOffsetTerm <- FALSE
  # nolint next: object_usage_linter. getTestOffset reads it, through eval.
  argumentOffset <- NULL
  # the 'offset' argument as written, re-evaluated on new rows
  # (offsetArgumentFormula)
  offsetArgument <- NULL
  # the 'weights' argument as written when it names data columns, re-read from
  # new rows (weightsArgumentFormula)
  weightsArgument <- NULL
  # the response's original type, recorded on the result so the fitters can
  # route family = "auto" and reject a categorical response an unsupported
  # family cannot fit; each y-producing branch below refreshes it
  responseInfo <- list(type = "numeric", n.levels = NA_integer_, levels = NULL)
  # the frame handed over and what reads it again (dataFrame.R); the positions
  # of the rows kept, and the names of the data's columns
  frameData <- NULL
  rowsKept <- NULL
  rowsTotal <- NA_integer_
  frameColumns <- NULL
  offsetExpression <- NULL
  # whether the build drew from R's generator where a frame is read again
  drawsAtBuild <- FALSE

  if (missing(formula)) {
    stop("first argument to dbartsData - 'formula'/'x.train' - must be present")
  }

  if (inherits(formula, "dbartsData")) {
    ignored <- c(
      "data",
      "test",
      "subset",
      "weights",
      "offset",
      "offset.test",
      "factors",
      "na.action",
      "counts"
    )
    ignored <- c(
      ignored[ignored %in% names(matchedCall)],
      intersect(c("n.cuts", "useQuantiles"), names(stated))
    )
    if (length(ignored) > 0L) {
      warning(
        "if data supplied as dbartsData, remaining arguments are ignored: ",
        paste0("'", ignored, "'", collapse = ", "),
        call. = FALSE
      )
    }
    return(formula)
  }
  if (!is.formula(formula)) {
    formula <- asNumericPredictors(formula)
  }

  if (is.formula(formula)) {
    # a hurdle fit reads its zero part's offset rows through a formula that
    # names no predictor, marked so; that data object never reaches a sampler
    readsRowsOnly <- isTRUE(attr(formula, "readsRowsOnly", exact = TRUE))
    # the fitting functions hand over a formula whose forest() terms they have
    # already read; one arriving here was given to a door that reads none
    refuseForestTerm(formula, "dbartsData")
    formulaGiven <- formula
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

    # 'data' is evaluated once, when this function first reads it, and that
    # value is handed on (handOn) in the caller's expression's place in every
    # call built from the matched one below: the model frame and the readings
    # of 'weights' and 'test'. Evaluated again by one of them, a 'data' that
    # draws its rows would give that reader other rows than the ones 'subset'
    # was read against. What else is put in the model frame's call below is
    # handed on the same way, so that the call names no row.
    if (!dataIsMissing) {
      matchedCall <- handOn(matchedCall, list(data = data))
      frameData <- data
      frameColumns <- if (is.environment(data)) {
        ls(data, all.names = TRUE)
      } else {
        names(data)
      }
    }
    # 'weights' as the length check below read it, where it did
    weightsValue <- NULL
    modelFrameArgs <- c("formula", "data", "subset", "weights", "offset")

    # whether the expressions a frame-keeping object reads again draw from R's
    # generator (the formula, the weights and the offset; not 'data', which
    # is kept, nor 'subset', whose rows are): a frame that cannot give the
    # same values twice is not kept
    seedBefore <- readGlobalSeed()

    ## extract offset prematurely, if necessary
    if (offsetIsMissing) {
      offset <- NULL
      modelFrameArgs <- c("formula", "data", "subset", "weights")
    } else {
      # as model.frame evaluates an offset: data columns first, every other
      # name from the formula's environment, or failing that the caller's;
      # a name found in none of them is an error, never a dropped offset
      written <- if (isDotsReference(matchedCall$offset)) {
        recoverForwardedArgument(matchedCall$offset, parent.frame())
      } else {
        list(expr = matchedCall$offset, env = environment(formula))
      }
      offsetEnv <- written$env
      offsetExpression <- written$expr
      offset <- evaluateOffsetExpression(
        written$expr,
        if (dataIsMissing) NULL else data,
        offsetEnv
      )
      if (inherits(offset, "error") && !identical(offsetEnv, parent.frame())) {
        offsetEnv <- parent.frame()
        offset <- evaluateOffsetExpression(
          written$expr,
          if (dataIsMissing) NULL else data,
          offsetEnv
        )
      }
      if (inherits(offset, "error")) {
        stop(
          "'offset' cannot be evaluated: ",
          conditionMessage(offset),
          call. = FALSE
        )
      }

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
        offsetArgument <- offsetArgumentFormula(
          written$expr,
          offsetEnv,
          if (!dataIsMissing) formulaDataNames(data),
          if (!dataIsMissing && is.data.frame(data)) {
            nrow(data)
          } else {
            NA_integer_
          }
        )
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

    if (!missing(weights) && !dataIsMissing) {
      weightsArgument <- weightsArgumentFormula(
        matchedCall$weights,
        environment(formula),
        formulaDataNames(data)
      )
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
      modelFrameCall <- handOn(modelFrameCall, list(offset = offset))
    }
    # 'weights' is evaluated once: what the length check read is what the
    # model frame is given
    if (!is.null(weightsValue)) {
      modelFrameCall <- handOn(modelFrameCall, list(weights = weightsValue))
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
        modelFrameCall <- handOn(modelFrameCall, list(data = data))
      }
    } else if (!dataIsMissing && (is.list(data) || is.environment(data))) {
      refuseSparseFormulaColumns(formula, data)
    }

    # 'subset' is evaluated once, here, as the model frame would evaluate it,
    # in the one value of 'data' and then the formula's environment, and every
    # reader below is given that value: the model frame, the range check and
    # the rows a basis is cut to. A subset that draws its rows, or has a side
    # effect, is then one draw for all of them. One that cannot be evaluated
    # here is left for the model frame to refuse in R's own words.
    subsetValue <- NULL
    drawsAtBuild <- !identical(readGlobalSeed(), seedBefore)
    if (!is.null(matchedCall$subset)) {
      taken <- tryCatch(
        list(
          subset = if (dataIsMissing) {
            eval(matchedCall$subset, environment(formula))
          } else {
            eval(matchedCall$subset, data, environment(formula))
          }
        ),
        error = function(e) NULL
      )
      if (!is.null(taken)) {
        subsetValue <- taken[[1L]]
        modelFrameCall <- handOn(modelFrameCall, taken["subset"])
      }
    }
    # an out-of-range 'subset' would otherwise reach the na.action as a set
    # of all-NA rows and be dropped in silence
    if (!is.null(subsetValue) && !dataIsMissing && is.data.frame(data)) {
      refuseOutOfRangeSubset(subsetValue, nrow(data), rownames(data))
    }

    # terms.formula() before R 4.3 refuses a plain list holding a bare
    # classed matrix
    if (!dataIsMissing && is.list(data) && !is.data.frame(data)) {
      data <- asDataFrameableList(data)
      modelFrameCall <- handOn(modelFrameCall, list(data = data))
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
        modelFrameCall <- handOn(
          modelFrameCall,
          list(dbartsSparseMissing = sparseMissing)
        )
      }
    }
    # the source row of each frame row, which 'subset' and the na.action shape
    # as they shape the frame, so that the rows kept are the frame's whatever
    # the na.action records: the frame's row names cannot serve, since a
    # repeated row is renamed ("1.1"). Data that are not a data frame are
    # counted where the frame reads them, after its variables (rowIndexColumn)
    rowCount <- new.env(parent = emptyenv())
    rowCount$total <- NA_integer_
    if (!dataIsMissing && is.data.frame(data)) {
      rowCount$total <- nrow(data)
      modelFrameCall <- handOn(
        modelFrameCall,
        list(dbartsRowIndex = seq_len(nrow(data)))
      )
    } else {
      modelFrameCall$dbartsRowIndex <- rowIndexColumn(formula, rowCount)
    }
    if (!dataIsMissing) {
      formulaVars <- all.vars(formula)
      refusePOSIXltColumns(
        data,
        if ("." %in% formulaVars) NULL else formulaVars
      )
    }
    seedBefore <- readGlobalSeed()
    modelFrame <- eval(modelFrameCall, parent.frame())
    drawsAtBuild <- drawsAtBuild || !identical(readGlobalSeed(), seedBefore)
    # the test frame is built from this call again, against rows these vectors
    # do not cover
    modelFrameCall$dbartsSparseMissing <- NULL
    modelFrameCall$dbartsRowIndex <- NULL
    naOmitted <- attr(modelFrame, "na.action")
    rowsTotal <- rowCount$total
    framePositions <- modelFrame[["(dbartsRowIndex)"]]
    modelFrame[["(dbartsRowIndex)"]] <- NULL
    frameTerms <- attr(modelFrame, "terms")
    attr(frameTerms, "dataClasses") <- attr(frameTerms, "dataClasses")[
      names(attr(frameTerms, "dataClasses")) != "(dbartsRowIndex)"
    ]
    attr(modelFrame, "terms") <- frameTerms
    rowsKept <- if (
      !is.null(framePositions) &&
        !identical(framePositions, seq_len(rowsTotal))
    ) {
      framePositions
    }
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
      # a status missing beside a time is a missing response, which only
      # na.pass lets reach here
      y <- log(survival$time)
      y[is.na(survival$status)] <- NA_real_
      responseInfo <- list(
        type = "numeric",
        n.levels = NA_integer_,
        levels = NULL
      )
    } else if (binomialResponse && is.matrix(y) && ncol(y) == 2L) {
      # a binomial's successes and failures, the rows the model frame kept;
      # its trials are the response the data object holds. A factor or
      # character column the expression reads is no count, whatever cbind
      # made of it
      refuseNonNumericBinomialColumns(
        formula[[2L]],
        if (dataIsMissing) NULL else data,
        environment(formula)
      )
      formulaCounts <- validateBinomialCounts(y)
      y <- as.double(rowSums(formulaCounts))
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
    counts <- formulaCounts
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

    ## offset: every offset() term plus the 'offset' argument, as in lm. The
    ## argument's own share is kept apart, since its default test offset
    ## follows its own rule (below), where the terms are evaluated on 'test'
    hasOffsetTerm <- !is.null(attr(terms(modelFrame), "offset"))
    if (identical(offsetGivenAsScalar, FALSE)) {
      argumentOffset <- as.vector(modelFrame[["(offset)"]])
      offset <- as.vector(model.offset(modelFrame))
    } else if (identical(offsetGivenAsScalar, TRUE)) {
      argumentOffset <- offset
      offset <- rep_len(offset, numObservations)
      if (hasOffsetTerm) {
        offset <- offset + as.vector(model.offset(modelFrame))
      }
    } else if (hasOffsetTerm) {
      offset <- as.vector(model.offset(modelFrame))
    }

    ## predictors
    modelTerms <- terms(modelFrame)
    # no term label is no predictor, an intercept or an offset() term
    # notwithstanding (y ~ 1, y ~ offset(o), y ~ . less every column); a
    # formula naming only sparse predictors rewrites to a dense-only RHS with
    # no terms of its own (above), empty here but not of predictors
    if (
      !readsRowsOnly &&
        length(attr(modelTerms, "term.labels")) == 0L &&
        length(sparseColumns) == 0L
    ) {
      stop("predictors must be specified for regression tree analysis")
    }

    termLabels <- attr(modelTerms, "term.labels")
    # before the backticks of a non-syntactic name are stripped below, since
    # the scan reads them as the quoting they are
    refuseInteractionTerms(termLabels)
    badLabels <- grepl("^`[^`]*`$", termLabels)
    if (sum(badLabels) > 0) {
      termLabels[badLabels] <- gsub("^`(.*)`$", "\\1", termLabels[badLabels])
    }

    predictorFrame <- modelFrame[termLabels]
    if (length(sparseColumns) > 0L) {
      # a sparse column is indexed by the source rows the frame carries,
      # which 'subset' and na.action have already shaped
      for (sparseName in names(sparseColumns)) {
        predictorFrame[[sparseName]] <-
          subsetSparseColumn(sparseColumns[[sparseName]], framePositions)
      }
    }
    x <- makeModelMatrix(predictorFrame)
    attr(x, "terms") <- predictorTerms(
      modelTerms,
      if (!dataIsMissing) formulaDataNames(data),
      if (!dataIsMissing && is.data.frame(data)) {
        nrow(data)
      } else if (is.null(matchedCall$subset)) {
        NROW(modelFrame) + length(naOmitted)
      } else {
        NA_integer_
      }
    )

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
    if (NCOL(formula) == 0L) {
      stop("predictors must be specified for regression tree analysis")
    }

    if (missing(subset) || is.null(subset)) {
      subset <- seq.int(length(y))
    }
    refuseOutOfRangeSubset(subset, initialNumObservations)
    y <- y[subset]
    x <- formula[subset, , drop = FALSE]
    trainRowNames <- rownames(x)
    countsRows <- initialNumObservations
    countsSubset <- subsetPositions(subset, countsRows, rownames(formula))

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
    if (!is.null(offset)) {
      offsetArgument <- offsetArgumentFormula(
        matchedCall$offset,
        parent.frame(),
        as.character(colnames(formula)),
        initialNumObservations
      )
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

    refuseDataFrameResponse(data)
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
    if (NCOL(formula) == 0L) {
      stop("predictors must be specified for regression tree analysis")
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
      # a formula fit's predictor frame handed on by bart()'s multinomial
      # door carries its terms (extractMultinomialFormulaData)
      frameTerms <- attr(formula, "dbartsTerms")
      formula <- makeModelMatrix(formula)
      attr(formula, "terms") <- frameTerms
    }
    xIsMixed <- inherits(formula, "dbartsMixedMatrix")
    x <- if (is.matrix(formula) || xIsMixed) {
      formula[subset, , drop = FALSE]
    } else {
      formula[subset]
    }
    countsRows <- initialNumObservations
    countsSubset <- subsetPositions(subset, countsRows, rownames(formula))

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
    if (!is.null(offset)) {
      offsetArgument <- offsetArgumentFormula(
        matchedCall$offset,
        parent.frame(),
        as.character(colnames(formula)),
        initialNumObservations
      )
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
    refuseMatrixOfCodes(test, x, "'test'")
    x.test <- validateXTest(test, x)
  }

  if (!is.null(x.test)) {
    if (testOffsetIsMissing && (hasOffsetTerm || !is.null(offsetArgument))) {
      ## the 'offset' argument re-evaluated on 'test', as predict re-evaluates
      ## it on newdata, plus the offset() terms evaluated there; it stays
      ## linked to the training offset only where it is that offset
      argumentTest <- evaluateOffsetArgument(offsetArgument, test)
      if (isFALSE(argumentTest)) {
        stop(
          "'offset' was given as ",
          describeOffsetArgument(offsetArgument),
          ", which cannot be evaluated on the rows of 'test'; give ",
          "'offset.test' for them"
        )
      }
      offset.test <- rep_len(
        addFormulaTermOffset(x, test, argumentTest, "offset", "test"),
        nrow(x.test)
      )
      testUsesRegularOffset <- !hasOffsetTerm &&
        if (isTRUE(offsetGivenAsScalar)) {
          length(argumentTest) == 1L &&
            identical(argumentTest, as.double(offset[1L]))
        } else {
          identical(argumentTest, as.double(offset))
        }
    } else if (testOffsetIsMissing) {
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
      # nolint next: object_usage_linter. getTestOffset reads it, through eval.
      testOffsetEnclos <- if (is.formula(formula)) {
        environment(formula)
      } else {
        callerEnv
      }
      testOffsetInfo <- eval(getTestOffset)

      offset.test <- testOffsetInfo$offset.test
      testUsesRegularOffset <- testOffsetInfo$testUsesRegularOffset

      if (
        !is.null(offset.test) &&
          !is.matrix(offset.test) &&
          !is.data.frame(offset.test)
      ) {
        # one number recycles; anything else is one per test row, refused
        # rather than recycled or cut, as predict refuses its 'offset'
        if (length(offset.test) != 1L && length(offset.test) != nrow(x.test)) {
          stop("'offset.test' must have the same number of rows as 'test'")
        }
        offset.test <- rep_len(offset.test, nrow(x.test))
      }
      if (
        hasOffsetTerm && !is.matrix(offset.test) && !is.data.frame(offset.test)
      ) {
        offset.test <- addFormulaTermOffset(
          x,
          test,
          offset.test,
          "offset.test",
          "test"
        )
        testUsesRegularOffset <- FALSE
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
  # 'weights' bound to NULL gave the training rows none, so there are none to
  # look for in the test rows
  # a multinomial fit's test weights are read where its formula is, before
  # the count response reaches here; a binomial's are read here, as a binary
  # fit's are
  if (
    !is.null(x.test) &&
      !is.null(matchedCall$weights) &&
      !is.null(weights) &&
      (is.null(counts) || binomialResponse)
  ) {
    if (!is.formula(formula)) {
      warning(
        "'weights' are ignored for test data when model is not specified as a formula; this only impacts extracting samples from the posterior predictive distribution of the test data",
        call. = FALSE
      )
    } else {
      # the weights alone are read, in the test rows as model.frame reads
      # them: 'test' first, then the formula's environment. A frame built
      # from the predictors would need columns 'test' need not hold (the
      # response), and the training call's offset and 'subset' belong to the
      # training rows
      testFrameCall <- handOn(
        as.call(list(
          quote(stats::model.frame),
          stats::reformulate("1", env = environment(formula)),
          weights = matchedCall$weights,
          na.action = stats::na.pass
        )),
        list(data = test)
      )
      tryResult <- tryCatch(
        eval(testFrameCall, parent.frame()),
        error = function(e) e
      )
      testWeights <- if (inherits(tryResult, "error")) {
        NULL
      } else {
        tryResult[["(weights)"]]
      }
      # as many weights as test rows, as model.frame requires of a column
      if (
        inherits(tryResult, "error") ||
          (!is.null(testWeights) && length(testWeights) != nrow(x.test))
      ) {
        warning(
          "weights specified but not found in test data - ignoring",
          call. = FALSE
        )
      } else {
        # the training weights' own checks and coding
        if (!is.null(testWeights)) {
          testWeights <- as.vector(testWeights)
          if (!is.numeric(testWeights)) {
            stop("'weights' must be a numeric vector")
          }
          testWeights <- as.double(testWeights)
        }
        weights.test <- testWeights
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
    # alongside counts that WAS supplied. Two columns may be a binomial's,
    # whose offset is flat; the family, resolved later, refuses it on a
    # multinomial
    if (!is.null(offset) && ncol(counts) != 2L) {
      refuseFlatOffsetOnMultinomial(offset)
    }
    if (!is.null(offset.test) && ncol(counts) != 2L) {
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
  # matched to the categories by name, as the sampler's setters and predict
  # match them
  if (!is.null(counts)) {
    offset.category <- alignSamplerCategoryColumns(
      offset.category,
      counts,
      "offset"
    )
    offset.category.test <- alignSamplerCategoryColumns(
      offset.category.test,
      counts,
      "offset.test"
    )
  }
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
  # an infinite offset poisons sigma and every fit the same way, and the
  # sampler's offset swap refuses one too
  if (!is.null(offset) && any(is.infinite(offset))) {
    stop("'offset' contains non-finite values")
  }

  # Precision-degenerate response: a large magnitude but tiny spread
  # quantizes to (near-)identical doubles before the engine ever sees it,
  # e.g. y in [1e15, 1e15 + 1e-3] rounds to one representable double. Doubles
  # near magnitude s are spaced ~2.22e-16 * s apart, so the 1e-10 threshold is
  # ~1e6x the ulp spacing, clear of legitimately discrete (binary, counts,
  # ordinal) responses. max(abs(y)) == 0 is guarded to avoid a 0/0. A
  # multinomial 'y' is the trials vector, not a modeled response, so the
  # check is skipped whenever 'counts' is carried, keyed on the data object
  # rather than on a family dbartsData has no formal for.
  yRange <- diff(range(y))
  yScale <- max(abs(y))
  if (is.null(counts) && yScale > 0 && yRange / yScale < 1e-10) {
    warning(
      paste0(
        responsePrecisionWarningStem,
        ", or nearly so, at double precision (",
        length(unique(y)),
        " distinct value(s) among ",
        length(y),
        " observations); center and/or rescale the response before fitting"
      ),
      call. = FALSE
    )
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
      counts,
      offset.category,
      offset.category.test,
      testUsesRegularOffset
    ),
    n.cuts = NA_integer_,
    sigma = NA_real_
  )
  if (!is.null(stated$n.cuts)) {
    counts <- stated$n.cuts
    # a count over every column is read here. One a hazard fit's sampler could
    # read, over the columns before the period it appends, which dbartsData()
    # cannot know, is held as written for the sampler that reads it, and a
    # sampler of any other kind refuses it there. One no sampler could read is
    # refused here, with the message of the plain reading.
    result@n.cuts <- tryCatch(
      readCutCounts(counts, result@x),
      error = function(plain) {
        tryCatch(
          readCutCounts(counts, result@x, hazard = TRUE),
          error = function(hazard) stop(plain)
        )
      }
    )
    attr(result, statedCountAttr) <- counts
  }
  if (!is.null(stated$useQuantiles)) {
    result@useQuantiles <- stated$useQuantiles
  }
  result@na.action <- naOmitted
  if (is.formula(formula) && length(formula) == 3L) {
    attr(result, "response.label") <- paste(
      deparse(formula[[2L]]),
      collapse = " "
    )
  }
  if (!is.formula(formula)) {
    rowsKept <- if (identical(as.integer(countsSubset), seq_len(countsRows))) {
      NULL
    } else {
      as.integer(countsSubset)
    }
  }
  result@rows <- rowsKept
  if (!is.null(rowsKept)) {
    result@rows.total <- as.integer(
      if (is.formula(formula)) rowsTotal else countsRows
    )
  }
  result@frame.columns <- frameColumns
  # whose names those are: a basis named later is refused naming the kind of
  # data that would be kept
  if (is.environment(frameData)) {
    attr(result, "frame.environment") <- TRUE
  }
  if (!is.null(offsetArgument)) {
    attr(result, "offset.argument") <- offsetArgument
  }
  if (!is.null(weightsArgument)) {
    attr(result, "weights.argument") <- weightsArgument
  }
  if (!is.null(trainRowNames) || !is.null(testRowNames)) {
    result@rowNames <- list(train = trainRowNames, test = testRowNames)
  }
  result@response.type <- responseInfo$type
  result@response.n.levels <- as.integer(responseInfo$n.levels)
  result@response.levels <- responseInfo$levels
  # a Surv formula response's raw time/status, for dbarts() to decode
  # (above) - not a slot, since every other creation route (the matrix
  # interface, a pre-built dbartsData) has no use for it and carries it
  # through its own channel instead (the data object's status)
  if (!is.null(survivalStatus)) {
    attr(result, "survivalStatus") <- survivalStatus
    attr(result, "survivalTime") <- survivalTime
    attr(result, "survivalTimeOmitted") <- survivalTimeOmitted
  }
  if (
    keep.data &&
      (is.data.frame(frameData) || is.list(frameData)) &&
      !drawsAtBuild
  ) {
    weightsExpression <- matchedCall$weights
    environment(formulaGiven) <- frozenBuildEnvironment(
      formulaGiven,
      list(weightsExpression, offsetExpression),
      frameColumns,
      callerEnv
    )
    result <- keepFrameInObject(
      result,
      frameData,
      list(
        formula = formulaGiven,
        weights = weightsExpression,
        offset = offsetExpression,
        factors = factors
      )
    )
  }
  result
}

## A ':'/'*' term expands to a label like "x1:x2" in term.labels, a column
## the model frame never carries (only the individual predictors it names);
## passed through unchecked, makeModelMatrix dies on "undefined columns
## selected", naming neither the term nor the unsupported syntax. A ':'
## inside backticks (one non-syntactic name), brackets or a call's
## parentheses (x[, 1:3], I(a:b)) is part of one term, not an interaction
## between terms; so is the '::' or ':::' of a package-qualified call
## (splines::ns(w, 3)), which is one operator, not two ':'.
refuseInteractionTerms <- function(termLabels) {
  if (any(vapply(termLabels, hasTopLevelColon, NA))) {
    stop(
      "':' and '*' terms are not supported in 'formula'; write each ",
      "predictor as its own term - poly(), ns(), log(), and offset() are ",
      "supported",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## Whether a formula term label has a ':' outside every bracket and outside
## backticks and quotes, which is what makes a label an interaction of terms.
hasTopLevelColon <- function(label) {
  characters <- strsplit(label, "", fixed = TRUE)[[1L]]
  depth <- 0L
  quote <- ""
  isColon <- characters == ":"
  # a ':' next to another is the namespace operator, not an interaction
  inNamespace <- isColon &
    (c(FALSE, isColon[-length(isColon)]) | c(isColon[-1L], FALSE))
  for (position in seq_along(characters)) {
    character <- characters[[position]]
    if (nzchar(quote)) {
      if (character == quote) {
        quote <- ""
      }
    } else if (character %in% c("`", "\"", "'")) {
      quote <- character
    } else if (character %in% c("(", "[", "{")) {
      depth <- depth + 1L
    } else if (character %in% c(")", "]", "}")) {
      depth <- depth - 1L
    } else if (character == ":" && depth == 0L && !inNamespace[[position]]) {
      return(TRUE)
    }
  }
  FALSE
}

## A formula that names a matrix prefixes its columns with the term (X gives
## X.a, X.b), so a matrix of the unprefixed columns fails to match by name.
## When every missing name is a test column's name behind one common prefix,
## says so and says the two ways out; otherwise nothing is added.
formulaPrefixHint <- function(missingNames, testNames) {
  prefixes <- vapply(
    missingNames,
    function(name) {
      stripped <- testNames[
        endsWith(name, testNames) & nchar(name) > nchar(testNames)
      ]
      if (length(stripped) != 1L) {
        return(NA_character_)
      }
      substr(name, 1L, nchar(name) - nchar(stripped))
    },
    ""
  )
  if (anyNA(prefixes) || length(unique(prefixes)) != 1L) {
    return("")
  }
  paste0(
    "; a formula prefixes the columns of a matrix with its term, here '",
    prefixes[[1L]],
    "', so name the columns of 'test' as the fit's (",
    toString(missingNames),
    "), or pass 'test' as an unnamed matrix to match by position"
  )
}
