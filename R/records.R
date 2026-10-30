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
## it from the model and data: the values written.
controlMirrorAttr <- "dbarts.mirror"

controlMirror <- function(control) {
  attr(control, controlMirrorAttr, exact = TRUE)
}

## The control's copies of what the model and data hold. A copy never states
## anything; the caller's own value in one of these slots does.
controlCopyNames <- c("n.trees", "useQuantiles", "n.cuts", "binary")

## How a control holds one of controlCopyNames: "written", the caller's
## dbartsControl() call named it or its slot differs from the value the sampler
## last mirrored there, or, on a control with no mirror, from a fresh control's;
## "carried", the slot equals the mirrored value; "silent", a fresh control's
## default. A written value equal to the mirror cannot be told from no edit.
controlValueKind <- function(control, name) {
  if (name %in% controlSuppliedSlots(control)) {
    return("written")
  }
  value <- methods::slot(control, name)
  mirror <- controlMirror(control)
  if (!is.null(mirror)) {
    return(if (identical(value, mirror[[name]])) "carried" else "written")
  }
  if (identical(value, freshControlValue(name))) "silent" else "written"
}

## A fresh control's value of one of controlCopyNames. No caller's call names
## binary, which only a slot edit writes.
freshControlValue <- function(name) {
  if (identical(name, "binary")) {
    FALSE
  } else {
    formals(dbarts::dbartsControl)[[name]]
  }
}

## The tree count at the level of the fitting function: the summed counts of
## the forests with no basis, the one forest's on a model of one, or the first
## forest's where every forest has a basis. NULL on a model recording none.
fittingLevelCount <- function(model, data) {
  counts <- lapply(model@forests, methods::slot, "n.trees")
  if (length(counts) == 0L || any(vapply(counts, is.null, NA))) {
    return(NULL)
  }
  counts <- unlist(counts)
  plain <- forestsWithNoBasis(model, data)
  if (length(plain) > 0L) sum(counts[plain]) else counts[1L]
}

## Which of a model's forests have no basis, by position.
forestsWithNoBasis <- function(model, data) {
  numForests <- length(model@forests)
  hasBasis <- rep_len(FALSE, numForests)
  bases <- data@bases
  if (!is.null(bases)) {
    hasBasis[seq_along(bases)] <- !vapply(bases, is.null, NA)
  }
  which(!hasBasis[seq_len(numForests)])
}

## Writes the control's n.trees, useQuantiles, n.cuts and binary from the model
## and data a sampler takes: the count at the fitting function's level, the
## data object's rule and cut count, and whether the family is binary, a model
## whose family is "auto" leaving the control's own. The mirror is written
## beside them and the slots stop reading as named by the caller.
mirrorControl <- function(control, model, data) {
  count <- fittingLevelCount(model, data)
  if (!is.null(count)) {
    control@n.trees <- count
  }
  control@useQuantiles <- isTRUE(data@useQuantiles)
  if (!identical(model@family, "auto")) {
    control@binary <- isBinaryFamily(model@family)
  }
  if (length(data@n.cuts) > 0L && !anyNA(data@n.cuts)) {
    cuts <- unique(data@n.cuts)
    control@n.cuts <- if (length(cuts) == 1L) cuts else data@n.cuts
  }
  named <- attr(control, controlSuppliedAttr, exact = TRUE)
  if (!is.null(named)) {
    attr(control, controlSuppliedAttr) <- setdiff(named, controlCopyNames)
  }
  attr(control, controlMirrorAttr) <- list(
    n.trees = control@n.trees,
    useQuantiles = control@useQuantiles,
    n.cuts = control@n.cuts,
    binary = control@binary
  )
  control
}

## The attribute a fitting function leaves on a control it took from a fit:
## the door and the copies it set aside, for the not-used warning once the
## call knows what it takes.
controlCarriedAttr <- "dbarts.carried"

## A control a door may read: what a fit carried states nothing, so each
## carried value goes back to a fresh control's and is kept aside. The mirror
## comes off. A control with no mirror comes back as it is.
releaseCarriedControl <- function(control, door) {
  if (is.null(controlMirror(control))) {
    return(control)
  }
  values <- list()
  for (name in setdiff(controlCopyNames, "binary")) {
    if (identical(controlValueKind(control, name), "carried")) {
      values[[name]] <- methods::slot(control, name)
      methods::slot(control, name, check = FALSE) <- freshControlValue(name)
    }
  }
  attr(control, controlMirrorAttr) <- NULL
  attr(control, controlCarriedAttr) <- list(door = door, values = values)
  control
}

## The not-used warning: on every call, where a control's value differs from
## what the call takes.
warnControlNotUsed <- function(name, value, what) {
  warnClassed(
    "dbartsControlNotUsedWarning",
    "the control's ",
    name,
    " (",
    paste(unique(value), collapse = ", "),
    ") is not used: ",
    what
  )
}

## Whether a control's value, recycled, is the one per entry of `held`: a cut
## count against a data object's per column, a count against xbart's grid.
sameRecycled <- function(value, held) {
  identical(
    rep_len(as.integer(value), length(held)),
    as.integer(held)
  )
}

## The not-used warning for the values a fitting function set aside
## (releaseCarriedControl), against the model and data it ended with.
warnCarriedNotUsed <- function(carried, model, data) {
  values <- carried$values
  door <- carried$door
  count <- fittingLevelCount(model, data)
  if ("n.trees" %in% names(values) && !identical(values$n.trees, count)) {
    plain <- length(forestsWithNoBasis(model, data))
    numForests <- length(model@forests)
    warnControlNotUsed(
      "n.trees",
      values$n.trees,
      if (numForests == 1L) {
        paste0(
          "this fit has ",
          count,
          " trees; state the count as ",
          door,
          "(n.trees = ) or forest(n.trees = )"
        )
      } else if (plain > 0L) {
        paste0(
          "this fit has ",
          count,
          " trees in its forest",
          if (plain > 1L) "s",
          " with no basis; state the count as ",
          door,
          "(n.trees = ) or forest(n.trees = )"
        )
      } else {
        paste0(
          "this fit's first forest has ",
          count,
          " trees; state the count as forest(n.trees = )"
        )
      }
    )
  }
  flat <- if (identical(door, "dbartsSpec")) "" else paste0(door, "(%s = ) or ")
  remedy <- function(name, what) {
    paste0(
      "; state the ",
      what,
      " as ",
      sprintf(flat, name),
      "dbartsData(",
      name,
      " = )"
    )
  }
  rule <- isTRUE(data@useQuantiles)
  if (
    "useQuantiles" %in% names(values) && !identical(values$useQuantiles, rule)
  ) {
    warnControlNotUsed(
      "useQuantiles",
      values$useQuantiles,
      paste0("this fit takes ", rule, remedy("useQuantiles", "rule"))
    )
  }
  if (
    "n.cuts" %in% names(values) && !sameRecycled(values$n.cuts, data@n.cuts)
  ) {
    warnControlNotUsed(
      "n.cuts",
      values$n.cuts,
      paste0(
        "this fit takes ",
        paste(unique(data@n.cuts), collapse = ", "),
        remedy("n.cuts", "count")
      )
    )
  }
  invisible(NULL)
}

## The model and data a sampler takes from a control, a model and a data
## object handed to it. A count or rule the caller wrote on the control is
## read, and a model that records no count takes the control's; anything else
## the model and data hold stands. A carried value that differs draws the
## not-used warning, after the model is known to be taken.
takeControlReadings <- function(control, model, data) {
  records <- model@forests
  if (length(records) == 0L) {
    records <- list(new("dbartsForestRecord"))
  }
  if (
    is.null(records[[1L]]@n.trees) ||
      identical(controlValueKind(control, "n.trees"), "written")
  ) {
    records[[1L]]@n.trees <- control@n.trees
  }
  model@forests <- records
  if (identical(controlValueKind(control, "useQuantiles"), "written")) {
    data@useQuantiles <- control@useQuantiles
  }
  list(model = model, data = data)
}

## The not-used warning at new("dbartsSampler", control, model, data): what the
## control carries beside the model and data taken, and any binary flag beside
## a model that names its family.
warnDirectDoorNotUsed <- function(control, model, data) {
  kind <- function(name) controlValueKind(control, name)
  count <- fittingLevelCount(model, data)
  if (
    identical(kind("n.trees"), "carried") && !identical(control@n.trees, count)
  ) {
    warnControlNotUsed(
      "n.trees",
      control@n.trees,
      paste0("a sampler assembled by hand takes the model's count, ", count)
    )
  }
  rule <- isTRUE(data@useQuantiles)
  if (
    identical(kind("useQuantiles"), "carried") &&
      !identical(control@useQuantiles, rule)
  ) {
    warnControlNotUsed(
      "useQuantiles",
      control@useQuantiles,
      paste0("a sampler assembled by hand takes the data object's rule, ", rule)
    )
  }
  if (
    identical(kind("n.cuts"), "carried") &&
      !anyNA(data@n.cuts) &&
      !sameRecycled(control@n.cuts, data@n.cuts)
  ) {
    warnControlNotUsed(
      "n.cuts",
      control@n.cuts,
      paste0(
        "a sampler assembled by hand takes the data object's count, ",
        paste(unique(data@n.cuts), collapse = ", ")
      )
    )
  }
  if (
    !identical(kind("binary"), "silent") &&
      !identical(model@family, "auto") &&
      !identical(control@binary, isBinaryFamily(model@family))
  ) {
    warnControlNotUsed(
      "binary",
      control@binary,
      paste0(
        "a sampler assembled by hand takes the model's family, \"",
        model@family,
        "\""
      )
    )
  }
  invisible(NULL)
}

## $setControl's reading of a control's count, rule or cut count against the
## sampler's: nothing where they agree, a refusal where the caller wrote
## another, and otherwise the name, the sampler's being kept. A silent value
## keeps it with no message.
setControlCopyReading <- function(newControl, control, name) {
  given <- methods::slot(newControl, name)
  own <- methods::slot(control, name)
  if (
    identical(given, own) ||
      (identical(name, "n.cuts") &&
        sameRecycled(given, own) &&
        sameRecycled(own, given))
  ) {
    return(NULL)
  }
  kind <- controlValueKind(newControl, name)
  if (identical(kind, "written")) {
    stop(
      "changing '",
      name,
      "' is not available on an existing sampler",
      call. = FALSE
    )
  }
  if (identical(kind, "carried")) name
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

## The moved structure of a model given to $setModel against the sampler's, as
## the text of a refusal or NULL: each forest's count, columns, interactions,
## blocks, basis and label, the variance forest's count and columns, the
## category count and the hazard periods. A model recording no forest, one made
## by hand, and a record stating no count state no change. What else a record
## holds is not structure: a parameter stays the sampler's and a move setting
## is installed (keepModelRecords).
movedStructureRefusal <- function(newModel, model) {
  shown <- function(value) {
    if (is.null(value)) "none" else paste(value, collapse = ", ")
  }
  differs <- function(what, own, given) {
    paste0(
      "$setModel cannot change ",
      what,
      ": this sampler's is ",
      shown(own),
      " and the model's is ",
      shown(given),
      "; make a new sampler instead"
    )
  }
  moved <- function(what) {
    paste0("$setModel cannot change ", what, "; make a new sampler instead")
  }
  if (length(newModel@forests) == 0L) {
    return(NULL)
  }
  if (length(newModel@forests) != length(model@forests)) {
    return(differs(
      "the number of forests",
      length(model@forests),
      length(newModel@forests)
    ))
  }
  for (index in seq_along(model@forests)) {
    own <- model@forests[[index]]
    given <- newModel@forests[[index]]
    forest <- paste0("forest ", index, "'s ")
    if (!is.null(given@n.trees) && !identical(given@n.trees, own@n.trees)) {
      return(differs(paste0(forest, "tree count"), own@n.trees, given@n.trees))
    }
    if (!identical(given@vars, own@vars)) {
      return(differs(paste0(forest, "'vars'"), own@vars, given@vars))
    }
    if (!identical(given@interactions, own@interactions)) {
      return(moved(paste0(forest, "'interactions'")))
    }
    if (!identical(given@blocks, own@blocks)) {
      return(moved(paste0(forest, "'blocks'")))
    }
    if (!identical(given@basis.term$label, own@basis.term$label)) {
      return(differs(
        paste0(forest, "'basis'"),
        own@basis.term$label,
        given@basis.term$label
      ))
    }
    if (!identical(given@label, own@label)) {
      return(differs(paste0(forest, "label"), own@label, given@label))
    }
  }
  ownVariance <- model@variance
  givenVariance <- newModel@variance
  if (is.null(ownVariance) != is.null(givenVariance)) {
    return(moved("whether there is a variance forest"))
  }
  if (!is.null(ownVariance)) {
    if (!identical(givenVariance@n.trees, ownVariance@n.trees)) {
      return(differs(
        "the variance forest's tree count",
        ownVariance@n.trees,
        givenVariance@n.trees
      ))
    }
    if (!identical(givenVariance@columns, ownVariance@columns)) {
      return(differs(
        "the columns of the variance forest",
        ownVariance@columns,
        givenVariance@columns
      ))
    }
  }
  if (!identical(newModel@n.categories, model@n.categories)) {
    return(differs(
      "the number of categories",
      model@n.categories,
      newModel@n.categories
    ))
  }
  if (!identical(newModel@hazard.periods, model@hazard.periods)) {
    return(moved("the hazard periods"))
  }
  NULL
}

## The model facts a sampler holds for its life, written onto a model given
## to $setModel from the one it replaces: the forests, the variance forest,
## the category count, the shape setting and the hazard periods. The move
## settings the given model's records state are the exception: they are not
## model, and the install takes them.
keepModelRecords <- function(newModel, model) {
  given <- newModel@forests
  variance <- newModel@variance
  newModel@forests <- model@forests
  newModel@variance <- model@variance
  newModel@n.categories <- model@n.categories
  newModel@shape <- model@shape
  newModel@hazard.periods <- model@hazard.periods
  if (length(given) > 0L) {
    for (index in seq_along(given)) {
      newModel@forests[[index]]@moves <- given[[index]]@moves
    }
    if (!is.null(variance)) {
      newModel@variance@moves <- variance@moves
    }
  }
  newModel
}
