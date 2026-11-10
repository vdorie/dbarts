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

## The slot a stand-in names is the whole of it: dim(), subsetting, arithmetic,
## mean(), c(), unique(), cbind(), sd() and the rest stop, so nothing reads a
## slot that holds no values as if it did. var() and paste() coerce in R's own
## code, which no method reaches.
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
            mean.dbartsStandIn <- c.dbartsStandIn <- rep.dbartsStandIn <-
              unique.dbartsStandIn <- as.double.dbartsStandIn <-
                cbind.dbartsStandIn <- rbind.dbartsStandIn <-
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

frameNumRows <- function(data) {
  if (is.data.frame(data)) nrow(data) else NROW(data[[1L]])
}

## The number of rows a formula fit was handed: a data frame's, and otherwise
## the length of what formulaRowsRead() reads; NA where that cannot be read.
## `data` is NULL where none was given.
formulaRowsTotal <- function(formula, data) {
  if (is.data.frame(data)) {
    return(nrow(data))
  }
  read <- formulaRowsRead(formula)
  if (is.null(read)) {
    return(NA_integer_)
  }
  tryCatch(
    NROW(eval(read, data, environment(formula))),
    error = function(e) NA_integer_
  )
}

## What gives the number of rows of data that are not a data frame: the
## response, read where the formula reads it, or the first variable a
## one-sided formula names (a list's first element need not be a column); NULL
## where a formula names none.
formulaRowsRead <- function(formula) {
  if (length(formula) == 3L) {
    return(formula[[2L]])
  }
  vars <- setdiff(all.vars(formula), ".")
  if (length(vars) > 0L) as.name(vars[1L])
}

## The row-number column of a model frame read from data that are not a data
## frame, as a call the frame evaluates after its variables, with the data it
## reads them from: it counts the rows where formulaRowsTotal() would, records
## the count in `rowCount$total` and gives the positions, or gives no column
## where the count cannot be read.
rowIndexColumn <- function(formula, rowCount) {
  read <- formulaRowsRead(formula)
  if (is.null(read)) {
    return(NULL)
  }
  count <- function(value) {
    rowCount$total <- tryCatch(NROW(value), error = function(e) NA_integer_)
    if (!is.na(rowCount$total)) seq_len(rowCount$total)
  }
  as.call(list(count, read))
}

## The number of rows a frame-keeping object was handed, as its build counted
## them: where every row is kept in order, the rows the build kept, recorded so
## that no reader evaluates the formula again.
frameRowsTotal <- function(object) {
  if (is.null(object@rows)) object@build$nobs else dataRowsTotal(object)
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
## name). A function from no package is held by value too (heldFunction), and a
## package's function by name. Its parent is the formula's environment where
## that is the global, a namespace or a package environment, which supply the
## package functions, and otherwise the global one.
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
  copies <- new.env(parent = emptyenv())
  copies$from <- list()
  copies$to <- list()
  for (name in setdiff(names, c(dataNames, "."))) {
    value <- get0(name, home)
    if (is.null(value)) {
      value <- get0(name, fallback)
    }
    holdValue(held, name, value, copies)
  }
  held
}

## Whether a function is a package's: a primitive, or one whose environment is
## a namespace or the base environment. A package function is read by name.
isPackageFunction <- function(fn) {
  env <- environment(fn)
  is.primitive(fn) ||
    identical(env, baseenv()) ||
    isNamespace(env) ||
    startsWith(environmentName(env), "package:")
}

## Holds `value` under `name` in `env`, as heldFunction holds a function; a
## package function the global environment finds as this one, and a name read
## nowhere, are not held.
holdValue <- function(env, name, value, copies) {
  if (is.null(value)) {
    return(invisible(NULL))
  }
  if (is.function(value)) {
    if (
      isPackageFunction(value) &&
        identical(value, get0(name, globalenv(), mode = "function"))
    ) {
      return(invisible(NULL))
    }
    value <- heldFunction(value, copies)
  }
  assign(name, value, envir = env)
  invisible(NULL)
}

## A function from no package, held by value: one defined at the top level or in
## a caller's function is copied onto an environment holding only the names it
## reads, each held as the build found it, so that neither a later change to
## the workspace nor the rest of a local frame is carried; one a package's code
## made (its environments reach a namespace before the global one) is held as
## it stands. `copies` pairs each function copied with its copy, so a function
## that reads itself, or two that read each other, are copied once.
heldFunction <- function(fn, copies) {
  if (isPackageFunction(fn)) {
    return(fn)
  }
  env <- environment(fn)
  scope <- env
  while (
    !identical(scope, globalenv()) &&
      !identical(scope, emptyenv()) &&
      !isNamespace(scope)
  ) {
    scope <- parent.env(scope)
  }
  if (!identical(scope, globalenv())) {
    return(fn)
  }
  for (index in seq_along(copies$from)) {
    if (identical(copies$from[[index]], fn)) {
      return(copies$to[[index]])
    }
  }
  scoped <- new.env(parent = globalenv())
  copy <- fn
  environment(copy) <- scoped
  copies$from <- c(copies$from, list(fn))
  copies$to <- c(copies$to, list(copy))
  reads <- setdiff(
    unique(c(all.names(body(fn)), unlist(lapply(formals(fn), all.names)))),
    names(formals(fn))
  )
  for (name in reads) {
    holdValue(scoped, name, get0(name, env), copies)
  }
  copy
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
  build$nobs <- NROW(result@y)
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
  if (!is.list(object@build)) {
    return("'build' must be a list")
  }
  rows <- object@rows
  if (!is.null(rows) && any(rows > dataRowsTotal(object), na.rm = TRUE)) {
    return("'rows' must be positions within the rows of 'frame'")
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
  if (dataKeepsFrame(object)) object@build$nobs else NROW(object@y)
}

## The model's variables over the rows in, as model.frame() of an lm fit has
## them: the response first, the predictors, '(weights)' and the offsets, with
## the terms record model.response() reads. From the frame where one is kept,
## its sparse columns as stored; otherwise decoded from the block where that is
## exact (the response to its recorded type and levels, a factor column as a
## factor), and the rest as held.
model.frame.dbartsData <- function(formula, ...) {
  object <- formula
  if (!dataKeepsFrame(object)) {
    return(decodeBlockFrame(object))
  }
  build <- object@build
  # a sparse column does not survive model.frame: it is read as its row
  # numbers, which 'subset' shapes as it shapes the frame, and taken as stored
  # at the rows they then name, never made dense
  stored <- object@frame
  sparse <- names(stored)[vapply(stored, isSparseDataFrameColumn, NA)]
  if (length(sparse) > 0L) {
    for (name in sparse) {
      stored[[name]] <- seq_len(frameNumRows(object@frame))
    }
  }
  held <- new.env(parent = environment(build$formula))
  held$.frame <- stored
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
      stop(
        "model.frame() cannot read this object's frame: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  terms <- attr(frame, "terms")
  classes <- attr(terms, "dataClasses")
  for (name in intersect(sparse, names(frame))) {
    frame[[name]] <- subsetSparseColumn(object@frame[[name]], frame[[name]])
    classes[[name]] <- "other"
  }
  attr(terms, "dataClasses") <- classes
  attr(frame, "terms") <- terms
  attr(frame, "na.action") <- object@na.action
  frame
}

## The frame of a data object that holds its block, in the formula's order:
## the response named as the formula wrote it ("y" where none did), each term
## as one column, a matrix term joined again under its own name, the columns
## no term names, '(weights)' and '(offset)'. A column
## the block holds sparse is decoded sparse. The terms record is the held one
## with the response put back, so that model.response() reads it, and is not
## given where no formula was or a column the formula does not name was added
## (a hazard fit's period).
decodeBlockFrame <- function(object) {
  y <- object@y
  response <- decodedResponse(object)
  x <- object@x
  predictors <- decodedPredictors(x)
  terms <- attr(x, "terms")
  classes <- attr(terms, "dataClasses")
  label <- dataResponseLabel(object)
  frame <- data.frame(row.names = seq_len(NROW(y)))
  frame[[if (is.null(label)) "y" else label]] <- response
  labels <- attr(x, "term.labels", exact = TRUE)
  if (is.null(labels)) {
    labels <- attr(terms, "term.labels")
  }
  labels <- sub("^`(.*)`$", "\\1", labels)
  formulaLabels <- sub("^`(.*)`$", "\\1", attr(terms, "term.labels"))
  placed <- character()
  # the columns the formula's terms name, which its record describes; a
  # sparse column, which the record does not hold, rides beside it
  described <- names(predictors)[vapply(predictors, isS4, NA)]
  for (term in labels) {
    parts <- if (term %in% names(predictors)) {
      frame[[term]] <- predictors[[term]]
      term
    } else {
      joinedTermParts(predictors, term, placed, classes[term])
    }
    if (length(parts) > 0L && !term %in% names(predictors)) {
      frame[[term]] <- joinedTermColumns(x, predictors, parts, term)
    }
    placed <- c(placed, parts)
    if (term %in% formulaLabels) {
      described <- c(described, parts)
    }
  }
  for (name in setdiff(names(predictors), placed)) {
    frame[[name]] <- predictors[[name]]
  }
  if (!is.null(object@weights)) {
    frame[["(weights)"]] <- object@weights
  }
  if (!is.null(object@offset)) {
    frame[["(offset)"]] <- object@offset
  }
  if (
    !is.null(terms) &&
      !is.null(label) &&
      all(names(predictors) %in% described)
  ) {
    attr(frame, "terms") <- termsWithResponse(terms, label)
  }
  attr(frame, "na.action") <- object@na.action
  frame
}

## The response of a block-holding object in its recorded type and levels.
decodedResponse <- function(object) {
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
  response
}

## The predictor columns of a block by name, a factor column as a factor over
## its declared levels and a column held sparse as sparse: a sparse vector, or
## a sparse factor where the column is a factor's.
decodedPredictors <- function(x) {
  if (!inherits(x, "dbartsMixedMatrix") || all(x$map > 0L)) {
    return(as.list(as.data.frame(x)))
  }
  declaredLevels <- attr(x, "factor.levels", exact = TRUE)
  numObservations <- x$numObservations
  columns <- lapply(seq_along(x$map), function(j) {
    levelsOfColumn <- if (j <= length(declaredLevels)) declaredLevels[[j]]
    if (x$map[j] > 0L) {
      column <- x$dense[[x$map[j]]]
      if (!is.null(levelsOfColumn) && !is.factor(column)) {
        column <- factor(levelsOfColumn[column + 1L], levels = levelsOfColumn)
      }
      return(column)
    }
    k <- -x$map[j]
    entries <- seq.int(x$sparse@p[k] + 1L, length.out = diff(x$sparse@p)[k])
    rows <- x$sparse@i[entries] + 1L
    values <- x$sparse@x[entries]
    if (is.null(levelsOfColumn)) {
      return(methods::new(
        "dsparseVector",
        i = rows,
        x = values,
        length = numObservations
      ))
    }
    reference <- x$sparseReference[k]
    sparseFactor(
      as.integer(values) + 1L,
      levels = levelsOfColumn,
      reference = if (!is.na(reference)) levelsOfColumn[reference + 1L],
      i = rows,
      length = numObservations
    )
  })
  names(columns) <- x$columnNames
  columns
}

## The block columns a matrix term was split into, named "term.suffix": two or
## more, or one where the term made a one-column matrix or the block holds it
## sparse; none where the term is no matrix.
joinedTermParts <- function(predictors, term, placed, class) {
  parts <- setdiff(
    names(predictors)[startsWith(names(predictors), paste0(term, "."))],
    placed
  )
  if (
    length(parts) > 1L ||
      (length(parts) == 1L &&
        (isS4(predictors[[parts]]) || identical(unname(class), "nmatrix.1")))
  ) {
    parts
  } else {
    character()
  }
}

## A matrix term's columns joined into one, named as the term named them: a
## sparse matrix where the block holds every part sparse.
joinedTermColumns <- function(x, predictors, parts, term) {
  suffixes <- substring(parts, nchar(term) + 2L)
  if (all(vapply(predictors[parts], isS4, NA))) {
    columns <- x$sparse[, -x$map[match(parts, x$columnNames)], drop = FALSE]
  } else {
    columns <- do.call(cbind, lapply(unname(predictors[parts]), as.vector))
  }
  dimnames(columns) <- list(NULL, suffixes)
  columns
}

## A one-sided terms record with `label`, the response as model.frame names it,
## put back in front of its variables, as the record of the two-sided formula.
termsWithResponse <- function(terms, label) {
  response <- str2lang(label)
  result <- call("~", response, terms[[2L]])
  attributes(result) <- attributes(terms)
  prepend <- function(variables) {
    as.call(c(as.list(variables)[1L], list(response), as.list(variables)[-1L]))
  }
  attr(result, "variables") <- prepend(attr(terms, "variables"))
  if (!is.null(attr(terms, "predvars"))) {
    attr(result, "predvars") <- prepend(attr(terms, "predvars"))
  }
  factors <- attr(terms, "factors")
  if (length(factors) > 0L) {
    factors <- rbind(0L, factors)
    rownames(factors)[1L] <- label
    attr(result, "factors") <- factors
  }
  ## The block holds every offset summed in one "(offset)" column, which
  ## model.offset() reads; an offset() term named here would count twice.
  attr(result, "offset") <- NULL
  attr(result, "response") <- 1L
  result
}

## The response's name as model.frame() gave it at the build, read from the
## held terms record, whose variables no longer hold it but whose data classes
## still name it first; NULL where no formula wrote one.
dataResponseLabel <- function(object) {
  terms <- attr(object@x, "terms")
  classes <- attr(terms, "dataClasses")
  if (length(classes) == 0L) {
    return(NULL)
  }
  first <- names(classes)[1L]
  variables <- vapply(
    as.list(attr(terms, "variables"))[-1L],
    function(variable) {
      paste(
        deparse(
          variable,
          width.cutoff = 500L,
          backtick = !is.symbol(variable) && is.language(variable)
        ),
        collapse = " "
      )
    },
    ""
  )
  if (first %in% variables || startsWith(first, "(")) NULL else first
}
