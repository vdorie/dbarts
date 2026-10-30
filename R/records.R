## The records a model and a data object carry of what a sampler runs: each
## mean forest's (dbartsForestRecord), the variance forest's
## (dbartsVarianceRecord), and the values the control mirrors from them.

## The record of a model of one forest: its count, what its caller stated of
## it, and its move settings.
plainForestRecord <- function(n.trees, stated = character(0L), moves = NULL) {
  record <- new("dbartsForestRecord")
  record@n.trees <- coerceOrError(n.trees, "integer")
  record@stated <- stated
  record@moves <- moves
  record
}

## The forests' configuration of a model of several forests as parallel lists,
## one entry per forest: `params`, the count, base, power and drawn flag of
## each (forestParams' layout), `sd`, `vars`, `labels`, `heldAmplitude`,
## `interactions` and `blocks`, with `amplitude.df` where some forest states
## one and `basisTerms` where some basis was written as code. NULL on a model
## of fewer than two forests.
forestsRecord <- function(model) {
  records <- model@forests
  if (length(records) < 2L) {
    return(NULL)
  }
  field <- function(name) lapply(records, methods::slot, name)
  result <- list(
    params = lapply(records, function(record) {
      as.double(c(
        record@n.trees,
        record@base,
        record@power,
        if (record@held.amplitude) 0 else 1
      ))
    }),
    sd = vapply(records, methods::slot, 0, "sd"),
    vars = field("vars"),
    labels = vapply(records, methods::slot, "", "label"),
    heldAmplitude = vapply(records, methods::slot, NA, "held.amplitude"),
    interactions = field("interactions"),
    blocks = field("blocks")
  )
  df <- vapply(records, methods::slot, 0, "amplitude.df")
  if (!all(is.na(df))) {
    result$amplitude.df <- df
  }
  terms <- field("basis.term")
  if (!all(vapply(terms, is.null, NA))) {
    result$basisTerms <- terms
  }
  result
}

## The move settings a model's forests state as one record: `forests`, the
## entry of each mean forest, and `variance`, the variance forest's or NULL;
## NULL where no forest states any.
modelMoveRecord <- function(model) {
  forests <- lapply(model@forests, methods::slot, "moves")
  variance <- if (!is.null(model@variance)) model@variance@moves
  if (all(vapply(forests, is.null, NA)) && is.null(variance)) {
    return(NULL)
  }
  list(forests = forests, variance = variance)
}

## Whether a model carries the forests of an amplitude coupling, two or more.
modelCarriesForests <- function(model) {
  length(model@forests) >= 2L
}

## The labels of a model's forests, NULL on a model of fewer than two.
modelForestLabels <- function(model) {
  forestsRecord(model)$labels
}

## The variance forest of a model as a list of its count, base, power and
## columns, NULL where it has none.
varianceRecordList <- function(model) {
  record <- model@variance
  if (is.null(record)) {
    return(NULL)
  }
  list(
    n.trees = record@n.trees,
    base = record@base,
    power = record@power,
    columns = record@columns
  )
}

## The attribute a control carries beside the slots the sampler writes into
## it from the model and data: the values written, the tree count a fit given
## the control inherits where its first forest has a basis, and `fit.state`,
## TRUE where the model or data hold a fact bart() and xbart() cannot take a
## control of.
controlMirrorAttr <- "dbarts.mirror"

controlMirror <- function(control) {
  attr(control, controlMirrorAttr, exact = TRUE)
}

## Writes the control's n.trees, useQuantiles and binary from the model and
## data a sampler takes: the first forest's count, the data object's rule, and
## whether the family is binary, a model whose family is "auto" leaving the
## control's own. `carried` is the count a later fit given the control
## inherits, NULL where it is the slot's; it defaults to the one the control
## already records.
mirrorControl <- function(
  control,
  model,
  data,
  carried = controlMirror(control)$control.n.trees
) {
  if (length(model@forests) > 0L && !is.null(model@forests[[1L]]@n.trees)) {
    control@n.trees <- model@forests[[1L]]@n.trees
  }
  control@useQuantiles <- isTRUE(data@useQuantiles)
  if (!identical(model@family, "auto")) {
    control@binary <- isBinaryFamily(model@family)
  }
  attr(control, controlMirrorAttr) <- list(
    n.trees = control@n.trees,
    useQuantiles = control@useQuantiles,
    binary = control@binary,
    control.n.trees = carried,
    fit.state = if (holdsFitState(model, data)) TRUE
  )
  control
}

## Whether a model and data hold a fact beyond a plain forest's count: several
## forests, a variance forest, a category count, a shape setting, hazard
## periods, a censoring status, a forests = declaration or a forest's own move
## settings.
holdsFitState <- function(model, data) {
  modelCarriesForests(model) ||
    !is.null(model@variance) ||
    !is.na(model@n.categories) ||
    !is.na(model@shape) ||
    !is.null(model@hazard.periods) ||
    !is.null(data@status) ||
    isTRUE(attr(model, "forests.declared", exact = TRUE)) ||
    any(vapply(model@forests, function(record) !is.null(record@moves), NA))
}

## The model and data a sampler takes from a control, a model and a data
## object handed to it: the control's tree count is the first forest's and
## its quantile rule the data object's, and a model that records no forest
## takes a record of its count.
takeControlReadings <- function(control, model, data) {
  records <- model@forests
  if (length(records) == 0L) {
    records <- list(new("dbartsForestRecord"))
  }
  records[[1L]]@n.trees <- control@n.trees
  model@forests <- records
  data@useQuantiles <- control@useQuantiles
  list(model = model, data = data)
}

## The family the bridge is handed for a model: its own, a model whose family
## is "auto" taking the control's binary flag as probit's.
bridgeFamily <- function(control, model) {
  if (!identical(model@family, "auto")) {
    model@family
  } else if (isTRUE(control@binary)) {
    "probit"
  } else {
    ""
  }
}

## The model facts a sampler holds for its life, written onto a model given
## to $setModel from the one it replaces: the forests, the variance forest,
## the category count, the shape setting and the hazard periods.
keepModelRecords <- function(newModel, model) {
  newModel@forests <- model@forests
  newModel@variance <- model@variance
  newModel@n.categories <- model@n.categories
  newModel@shape <- model@shape
  newModel@hazard.periods <- model@hazard.periods
  attr(newModel, "forests.declared") <- attr(
    model,
    "forests.declared",
    exact = TRUE
  )
  newModel
}
