## A data object built from a formula and a data frame (or list) keeps the
## frame, the positions of the rows it kept and what it takes to read the
## model's variables from them again; its x, y, weights and offset slots hold
## a stand-in that stops when read. dataKeepsFrame() tells the two kinds of
## object apart, and resolveDataObject() makes a block-holding object (the
## kind the sampler's data and keep.data = FALSE give) from a frame-keeping one
## by running the build's own formula branch on the stored frame.

dataKeepsFrame <- function(object) {
  methods::.hasSlot(object, "frame") && !is.null(object@frame)
}

## The slot a stand-in names is the whole of it: dim(), subsetting, arithmetic
## and the rest stop, so nothing reads a slot that holds no values as if it did.
frameStandIn <- function(slot) {
  structure(slot, class = "dbartsStandIn")
}

refuseStandInRead <- function(...) {
  stop("the data are read with model.frame(o)", call. = FALSE)
}

dim.dbartsStandIn <- dimnames.dbartsStandIn <- length.dbartsStandIn <-
  `[.dbartsStandIn` <- `[[.dbartsStandIn` <- `$.dbartsStandIn` <-
    as.matrix.dbartsStandIn <- as.vector.dbartsStandIn <- mtfrm.dbartsStandIn <-
      Ops.dbartsStandIn <- Math.dbartsStandIn <- Summary.dbartsStandIn <-
        is.na.dbartsStandIn <- anyNA.dbartsStandIn <- is.finite.dbartsStandIn <-
          is.infinite.dbartsStandIn <- is.nan.dbartsStandIn <-
            refuseStandInRead

format.dbartsStandIn <- function(x, ...) {
  paste0("<'", unclass(x), "' is read with model.frame(o)>")
}

print.dbartsStandIn <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

str.dbartsStandIn <- function(object, ...) {
  cat(" ", format(object), "\n", sep = "")
  invisible(NULL)
}

## The positions, within the frame handed over, of the rows a build kept: those
## the one value of 'subset' selects (as model.frame reads it: a logical, a
## number or a row name), less the rows the na.action then dropped, which it
## numbers within that selection. NULL when every row is kept in order.
keptRowPositions <- function(subsetValue, total, rowNames, naOmitted) {
  if (is.na(total)) {
    return(NULL)
  }
  positions <- seq_len(total)
  if (!is.null(subsetValue)) {
    positions <- if (is.character(subsetValue)) {
      match(subsetValue, rowNames)
    } else {
      positions[subsetValue]
    }
  }
  if (length(naOmitted) > 0L) {
    positions <- positions[-as.integer(unclass(naOmitted))]
  }
  if (identical(positions, seq_len(total))) {
    return(NULL)
  }
  as.integer(positions)
}

frameNumRows <- function(data) {
  if (is.data.frame(data)) nrow(data) else NROW(data[[1L]])
}

## The kept rows' names as model.frame() gives them, which the rows alone
## determine unless a repeated row is one the na.action then drops.
derivedRowNames <- function(frame, rows) {
  names <- if (is.data.frame(frame)) {
    rownames(frame)
  } else {
    as.character(seq_len(frameNumRows(frame)))
  }
  make.unique(if (is.null(rows)) names else names[rows])
}

## An environment holding the caller's variables a build's expressions name that
## the data do not, as they were at the build, so that reading the frame again
## finds the values the build found (dbartsData() keeps a caller's value, not its
## name), and the functions defined locally among them. Its parent is the
## formula's environment where that is the global, a namespace or a package
## environment, which supply the functions, and otherwise the global one.
frozenBuildEnvironment <- function(formula, expressions, dataNames, fallback) {
  names <- unique(unlist(
    lapply(c(list(formula), expressions), all.names),
    use.names = FALSE
  ))
  home <- environment(formula)
  # a local environment is not kept whole, for it holds whatever else its
  # function did; what the expressions name from it is held instead
  keeps <- identical(home, globalenv()) ||
    isNamespace(home) ||
    startsWith(environmentName(home), "package:")
  held <- new.env(parent = if (keeps) home else globalenv())
  for (name in setdiff(names, c(dataNames, "."))) {
    value <- get0(name, home)
    if (is.null(value)) {
      value <- get0(name, fallback)
    }
    # a function the global environment would find as this one is not held
    found <- is.null(value) ||
      (is.function(value) &&
        identical(value, get0(name, globalenv(), mode = "function")))
    if (!found) {
      assign(name, value, envir = held)
    }
  }
  held
}

## What the object keeps beside its frame in place of the block: the stand-ins,
## the frame, the record that reads it again, and the row names the frame does
## not give.
keepFrameInObject <- function(result, frame, build) {
  result@frame <- frame
  # the coding record: the block with no rows, and whether each column held a
  # missing value, which is all makeTestModelMatrix() reads of the block
  x <- result@x
  seen <- dataMissingSeen(result)
  flags <- vapply(
    seq_len(ncol(x)),
    function(j) sourceColumnHasNA(x, j, ncol(x), nrow(x)),
    NA
  )
  build$coding <- list(
    x = x[0L, , drop = FALSE],
    missing.seen = if (is.null(seen)) flags else flags | seen
  )
  result@build <- build
  result@x <- frameStandIn("x")
  result@y <- frameStandIn("y")
  if (!is.null(result@weights)) {
    result@weights <- frameStandIn("weights")
  }
  if (!is.null(result@offset)) {
    result@offset <- frameStandIn("offset")
  }
  rowNames <- result@rowNames
  if (identical(rowNames$train, derivedRowNames(frame, result@rows))) {
    rowNames$train <- NULL
  }
  result@rowNames <- if (is.null(rowNames$train) && is.null(rowNames$test)) {
    NULL
  } else {
    rowNames
  }
  # the status and time are read again from the response
  attr(result, "survivalStatus") <- NULL
  attr(result, "survivalTime") <- NULL
  result
}

## What a block-holding object's slots hold, which the slot classes, widened to
## admit a stand-in, no longer say.
validateBlockSlotTypes <- function(object) {
  if (!is.numeric(object@y)) {
    return("'y' must be numeric")
  }
  for (name in c("weights", "offset")) {
    value <- methods::slot(object, name)
    if (!is.null(value) && !is.numeric(value)) {
      return(paste0("'", name, "' must be null or numeric"))
    }
  }
  NULL
}

validateFrameKeepingData <- function(object) {
  for (name in c("x", "y")) {
    if (!inherits(methods::slot(object, name), "dbartsStandIn")) {
      return(paste0(
        "'",
        name,
        "' must be a stand-in on an object that keeps a frame"
      ))
    }
  }
  for (name in c("weights", "offset")) {
    value <- methods::slot(object, name)
    if (!is.null(value) && !inherits(value, "dbartsStandIn")) {
      return(paste0(
        "'",
        name,
        "' must be null or a stand-in on an object that keeps a frame"
      ))
    }
  }
  frame <- object@frame
  if (!is.list(frame)) {
    return("'frame' must be a data frame or a list")
  }
  total <- frameNumRows(frame)
  rows <- object@rows
  if (
    !is.null(rows) && (!is.integer(rows) || any(rows > total, na.rm = TRUE))
  ) {
    return("'rows' must be positions within the rows of 'frame'")
  }
  if (!is.list(object@build)) {
    return("'build' must be a list")
  }
  TRUE
}

## The block-holding object a frame-keeping one stands for: dbartsData()'s own
## formula branch run again on the stored frame, with the frozen variables, the
## kept positions as 'subset' and every row kept as it stands (the na.action
## ran at the build and its record is put back), so that the block, the response,
## the weights, the offset and a Surv response's attributes come from the build's
## own code in the build's order. The build's warnings were given once, at the
## build, and are not given again. Any other object passes as it is.
resolveDataObject <- function(object) {
  if (!inherits(object, "dbartsData") || !dataKeepsFrame(object)) {
    return(object)
  }
  build <- object@build
  held <- new.env(parent = environment(build$formula))
  held$.frame <- object@frame
  held$.formula <- build$formula
  arguments <- list(
    formula = quote(.formula),
    data = quote(.frame),
    na.action = quote(stats::na.pass),
    factors = build$factors,
    keep.data = FALSE
  )
  if (!is.null(object@rows)) {
    arguments$subset <- object@rows
  }
  arguments$weights <- build$weights
  arguments$offset <- build$offset
  rebuilt <- withCallingHandlers(
    eval(as.call(c(list(quoteInNamespace(dbartsData)), arguments)), held),
    warning = function(w) invokeRestart("muffleWarning")
  )
  for (name in c("x", "y", "weights", "offset", "varTypes")) {
    methods::slot(object, name) <- methods::slot(rebuilt, name)
  }
  rowNames <- object@rowNames
  train <- if (is.null(rowNames$train)) {
    rebuilt@rowNames$train
  } else {
    rowNames$train
  }
  object@rowNames <- if (is.null(train) && is.null(rowNames$test)) {
    NULL
  } else {
    list(train = train, test = rowNames$test)
  }
  object@frame <- NULL
  object@build <- NULL
  attr(object, "survivalStatus") <- attr(rebuilt, "survivalStatus")
  attr(object, "survivalTime") <- attr(rebuilt, "survivalTime")
  object
}

## What a data object's rows number: every row in, whatever its weight.
nobs.dbartsData <- function(object, ...) {
  if (dataKeepsFrame(object)) {
    return(
      if (is.null(object@rows)) {
        frameNumRows(object@frame)
      } else {
        length(object@rows)
      }
    )
  }
  NROW(object@y)
}

## The model's variables over the rows in, as model.frame() of an lm fit has
## them: the response first, the predictors, '(weights)' and the offsets. From
## the frame where one is kept; otherwise decoded from the block where that is
## exact (the response to its recorded type and levels, a factor column as a
## factor), and the rest as held.
model.frame.dbartsData <- function(formula, ...) {
  object <- formula
  if (!dataKeepsFrame(object)) {
    return(decodeBlockFrame(object))
  }
  build <- object@build
  held <- new.env(parent = environment(build$formula))
  held$.frame <- object@frame
  held$.formula <- build$formula
  arguments <- list(
    quote(.formula),
    data = quote(.frame),
    na.action = quote(stats::na.pass),
    drop.unused.levels = FALSE
  )
  if (!is.null(object@rows)) {
    arguments$subset <- object@rows
  }
  arguments$weights <- build$weights
  arguments$offset <- build$offset
  frame <- tryCatch(
    eval(as.call(c(list(quote(stats::model.frame)), arguments)), held),
    error = function(e) {
      # a sparse column does not survive model.frame; the block the build
      # makes holds it, so the decoded frame answers as a sampler's data do
      if (any(vapply(object@frame, isS4, NA))) {
        NULL
      } else {
        stop(
          "model.frame() cannot read this object's frame: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    }
  )
  if (is.null(frame)) {
    return(decodeBlockFrame(resolveDataObject(object)))
  }
  attr(frame, "na.action") <- object@na.action
  frame
}

decodeBlockFrame <- function(object) {
  y <- object@y
  response <- switch(
    object@response.type,
    factor = ,
    `ordered factor` = factor(
      object@response.levels[y + 1L],
      levels = object@response.levels,
      ordered = object@response.type == "ordered factor"
    ),
    logical = y == 1,
    character = object@response.levels[y + 1L],
    y
  )
  if (!is.null(object@counts)) {
    response <- object@counts
  }
  if (!is.null(attr(object, "survivalStatus", exact = TRUE))) {
    response <- survival::Surv(
      attr(object, "survivalTime", exact = TRUE),
      attr(object, "survivalStatus", exact = TRUE)
    )
  }
  predictors <- as.data.frame(object@x)
  label <- attr(object, "response.label", exact = TRUE)
  frame <- data.frame(row.names = seq_len(NROW(y)))
  frame[[if (is.null(label)) "y" else label]] <- response
  # a matrix term is one column again, named as its columns were by the term
  joined <- character()
  for (term in attr(attr(object@x, "terms"), "term.labels")) {
    parts <- names(predictors)[startsWith(names(predictors), paste0(term, "."))]
    if (length(parts) > 1L && !term %in% names(predictors)) {
      columns <- as.matrix(predictors[parts])
      colnames(columns) <- substring(parts, nchar(term) + 2L)
      frame[[term]] <- columns
      joined <- c(joined, parts)
    }
  }
  for (name in setdiff(names(predictors), joined)) {
    frame[[name]] <- predictors[[name]]
  }
  if (!is.null(object@weights)) {
    frame[["(weights)"]] <- object@weights
  }
  if (!is.null(object@offset)) {
    frame[["(offset)"]] <- object@offset
  }
  attr(frame, "na.action") <- object@na.action
  frame
}
