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
controlCopyNames <- c("n.trees", "useQuantiles", "binary")

## The names a caller can write on a control for the model and data to hold,
## read for one release: the count and the rule on their slots, and the cut
## count on an attribute, which a control keeps only until a sampler stores it.
controlWrittenNames <- c("n.trees", "n.cuts", "useQuantiles")
controlCutsAttr <- "n.cuts"

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

## A fresh control's value of one of controlCopyNames.
freshControlValue <- function(name) {
  list(n.trees = 75L, useQuantiles = FALSE, binary = FALSE)[[name]]
}

## What the caller wrote on a control of controlWrittenNames, by name, judged as
## the values it stands for: a slot takes any value of its class.
writtenControlCopies <- function(control) {
  written <- list()
  for (name in c("n.trees", "useQuantiles")) {
    if (identical(controlValueKind(control, name), "written")) {
      written[[name]] <- methods::slot(control, name)
    }
  }
  written[["n.cuts"]] <- attr(control, controlCutsAttr, exact = TRUE)
  checked <- checkFlatCopies(
    written$n.trees,
    written$n.cuts,
    written$useQuantiles
  )
  checked[names(written)]
}

## The deprecation warning of a value written on a control, once a session per
## name wherever it is written: the constructor, or a slot edited on the way to
## a door.
warnWrittenCopies <- function(names) {
  for (name in names) {
    warnOnce(
      paste0("tombstone.control.", name),
      "dbartsControl(",
      name,
      " = ) is now ",
      if (identical(name, "n.trees")) {
        "dbarts(n.trees = ), dbartsSpec(n.trees = ) or forest(n.trees = )"
      } else {
        paste0("dbarts(", name, " = ) or dbartsData(", name, " = )")
      },
      ": the value was read, and on a control it is refused from dbarts ",
      tombstoneExpiry,
      ".",
      class = "dbartsDeprecatedWarning"
    )
  }
}

## The count, cut count and rule a caller gave as values: each coerced and
## judged, NULL standing for not given and left out of the list returned.
checkFlatCopies <- function(
  n.trees = NULL,
  n.cuts = NULL,
  useQuantiles = NULL
) {
  result <- list()
  if (!is.null(n.trees)) {
    n.trees <- coerceOrError(n.trees, "integer")
    if (length(n.trees) != 1L) {
      stop("'n.trees' must be of length 1", call. = FALSE)
    }
    if (is.na(n.trees) || n.trees <= 0L) {
      stop("'n.trees' must be a positive integer", call. = FALSE)
    }
    result$n.trees <- n.trees
  }
  if (!is.null(n.cuts)) {
    n.cuts <- coerceOrError(n.cuts, "integer")
    if (length(n.cuts) == 0L || anyNA(n.cuts) || any(n.cuts <= 0L)) {
      stop("'n.cuts' must contain only positive integers", call. = FALSE)
    }
    result$n.cuts <- n.cuts
  }
  if (!is.null(useQuantiles)) {
    useQuantiles <- coerceOrError(useQuantiles, "logical")
    if (length(useQuantiles) != 1L) {
      stop("'useQuantiles' must be of length 1", call. = FALSE)
    }
    if (is.na(useQuantiles)) {
      stop("'useQuantiles' must be TRUE/FALSE", call. = FALSE)
    }
    result$useQuantiles <- useQuantiles
  }
  result
}

## A door's reading of what a control has written beside the door's own values
## (checkFlatCopies): a written value stands for one not given, and one given
## beats it, which then waits with the carried ones for the not-used warning.
## Returns the values and that record.
readWrittenCopies <- function(control, flat) {
  written <- writtenControlCopies(control)
  warnWrittenCopies(names(written))
  carried <- attr(control, controlCarriedAttr, exact = TRUE)
  for (name in names(written)) {
    if (is.null(flat[[name]])) {
      flat[[name]] <- written[[name]]
    } else {
      carried$values[[name]] <- written[[name]]
    }
  }
  list(flat = flat, carried = carried)
}

## The control a door goes on with, holding nothing of controlWrittenNames.
clearWrittenCopies <- function(control) {
  control@n.trees <- freshControlValue("n.trees")
  control@useQuantiles <- freshControlValue("useQuantiles")
  attr(control, controlCutsAttr) <- NULL
  named <- attr(control, controlSuppliedAttr, exact = TRUE)
  if (!is.null(named)) {
    attr(control, controlSuppliedAttr) <- setdiff(named, controlWrittenNames)
  }
  control
}

## The attribute bart() leaves on the control it builds with the count and the
## cut settings it passes flat, for the call that builds the sampler to read.
wrapperCopiesAttr <- "dbarts.flat"

## What dbarts() and dbartsSpec() take of the count and the cut settings: the
## control they go on with, its slot holding the count the forest with no basis
## takes, its plainStatedAttr record what the call stated, and `n.cuts`,
## `useQuantiles` and `stated` for applyCutSettings. A control that holds a
## record comes from a wrapper, which has read its caller's control and passes
## its values flat. `written` holds the names the caller wrote.
readDoorCopies <- function(control, flat, written) {
  flat <- checkFlatCopies(flat$n.trees, flat$n.cuts, flat$useQuantiles)
  stated <- attr(control, plainStatedAttr, exact = TRUE)
  if (is.null(stated)) {
    # a written NULL is the name left out
    stated <- plainForestStated(
      control,
      setdiff(written, setdiff(controlWrittenNames, names(flat)))
    )
    read <- readWrittenCopies(control, flat)
    flat <- read$flat
    control <- clearWrittenCopies(control)
    attr(control, plainStatedAttr) <- stated
    attr(control, controlCarriedAttr) <- read$carried
  }
  if (!is.null(flat$n.trees)) {
    control@n.trees <- flat$n.trees
  }
  list(
    control = control,
    n.cuts = flat$n.cuts,
    useQuantiles = flat$useQuantiles,
    stated = stated
  )
}

## The cut settings a data object goes on with: a stated count or rule is used,
## and beside a grid the data object carries, the grid is derived again at it. A
## count the data object holds otherwise stands, recycled over any column a
## hazard fit has appended. With none, an unstated count is `copies$n.cuts`
## where a wrapper passed its own, and 100; an unstated rule stays NA.
applyCutSettings <- function(data, copies) {
  numColumns <- ncol(data@x)
  cutsStated <- "n.cuts" %in% names(copies$stated)
  ruleStated <- "useQuantiles" %in% names(copies$stated)
  if (cutsStated) {
    data@n.cuts <- recycleNumCuts(copies$n.cuts, numColumns)
  } else if (length(data@n.cuts) == 0L || anyNA(data@n.cuts)) {
    data@n.cuts <- recycleNumCuts(
      if (is.null(copies$n.cuts)) 100L else copies$n.cuts,
      numColumns
    )
  } else if (length(data@n.cuts) != numColumns) {
    data@n.cuts <- rep_len(data@n.cuts, numColumns)
  }
  if (ruleStated) {
    data@useQuantiles <- copies$useQuantiles
  }
  if (cutsStated || ruleStated) {
    data@cutPoints <- NULL
  }
  data
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

## Writes the control's n.trees, useQuantiles and binary from the model and data
## a sampler takes: the count at the fitting function's level, the data object's
## rule, and whether the family is binary, a model whose family is "auto"
## leaving the control's own. The mirror is written beside them, the slots stop
## reading as named by the caller and a written cut count comes off.
mirrorControl <- function(control, model, data) {
  count <- fittingLevelCount(model, data)
  if (!is.null(count)) {
    control@n.trees <- count
  }
  control@useQuantiles <- isTRUE(data@useQuantiles)
  if (!identical(model@family, "auto")) {
    control@binary <- isBinaryFamily(model@family)
  }
  attr(control, controlCutsAttr) <- NULL
  named <- attr(control, controlSuppliedAttr, exact = TRUE)
  if (!is.null(named)) {
    attr(control, controlSuppliedAttr) <- setdiff(
      named,
      c(controlCopyNames, controlWrittenNames)
    )
  }
  attr(control, controlMirrorAttr) <- list(
    n.trees = control@n.trees,
    useQuantiles = control@useQuantiles,
    binary = control@binary,
    n.cuts = data@n.cuts
  )
  control
}

## The attribute a fitting function leaves on a control it took from a fit:
## the door and the copies it set aside, for the not-used warning once the
## call knows what it takes.
controlCarriedAttr <- "dbarts.carried"

## A control a door may read: what a fit carried states nothing, so each
## carried value goes back to a fresh control's and is kept aside, and an edit
## of a fit's value goes on as named by the caller, whatever it equals. The
## mirror comes off. A control handed on by a wrapper keeps the wrapper's record.
releaseCarriedControl <- function(control, door) {
  carried <- attr(control, controlCarriedAttr, exact = TRUE)
  if (is.null(carried)) {
    carried <- list(door = door, values = list())
  }
  if (!is.null(controlMirror(control))) {
    named <- controlSuppliedSlots(control)
    for (name in setdiff(controlCopyNames, "binary")) {
      kind <- controlValueKind(control, name)
      if (identical(kind, "carried")) {
        carried$values[[name]] <- methods::slot(control, name)
        methods::slot(control, name) <- freshControlValue(name)
      } else if (identical(kind, "written")) {
        named <- union(named, name)
      }
    }
    carried$values$n.cuts <- controlMirror(control)$n.cuts
    attr(control, controlSuppliedAttr) <- named
    attr(control, controlMirrorAttr) <- NULL
  }
  attr(control, controlCarriedAttr) <- carried
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
  remedy <- function(name, what) {
    paste0(
      "; state the ",
      what,
      " as ",
      if (!identical(door, "dbartsSpec")) paste0(door, "(", name, " = ) or "),
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
## object handed to it. A model that records no count takes the control's. A
## count or rule written on the control is used where the model or the data
## object states none, the grid then derived again under a rule, and refused by
## name where it states another; a count has no forest to go to where every
## forest has a basis. A data object taken from a survival fit keeps its status
## only under a survival model. Anything else the model and data hold stands.
takeControlReadings <- function(control, model, data) {
  if (!identical(model@family, "aft")) {
    data@status <- NULL
  }
  written <- writtenControlCopies(control)
  if (length(model@forests) == 0L) {
    model@forests <- list(new("dbartsForestRecord"))
  }
  records <- model@forests
  refuse <- function(name, what) {
    stop(
      "the control's '",
      name,
      "' (",
      methods::slot(control, name),
      ") ",
      what,
      "; leave it out of dbartsControl()",
      call. = FALSE
    )
  }
  if (!is.null(written$n.trees)) {
    plain <- forestsWithNoBasis(model, data)
    if (length(plain) == 0L) {
      refuse(
        "n.trees",
        "is the tree count of the forest with no basis, and every forest of this model has a basis"
      )
    }
    own <- records[[plain[1L]]]
    if (!("n.trees" %in% own@stated)) {
      own@n.trees <- written$n.trees
      own@stated <- c(own@stated, "n.trees")
      records[[plain[1L]]] <- own
    } else if (!identical(own@n.trees, written$n.trees)) {
      refuse(
        "n.trees",
        paste0(
          "differs from the count the model states (",
          own@n.trees,
          "), which a sampler assembled by hand takes"
        )
      )
    }
  } else if (is.null(records[[1L]]@n.trees)) {
    records[[1L]]@n.trees <- control@n.trees
  }
  model@forests <- records
  if (!is.null(written$useQuantiles)) {
    if (is.na(data@useQuantiles)) {
      data@useQuantiles <- written$useQuantiles
      data@cutPoints <- NULL
    } else if (!identical(data@useQuantiles, written$useQuantiles)) {
      refuse(
        "useQuantiles",
        paste0(
          "differs from the rule the data object states (",
          data@useQuantiles,
          "), which a sampler assembled by hand takes"
        )
      )
    }
  }
  list(model = model, data = data)
}

## The warnings at new("dbartsSampler", control, model, data): the deprecation
## of what the caller wrote, and the not-used warning for what the control
## carries beside the model and data taken, a written cut count, which the data
## object's governs there, and any binary flag beside a model that names its
## family.
warnDirectDoorNotUsed <- function(control, model, data) {
  written <- writtenControlCopies(control)
  warnWrittenCopies(names(written))
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
  cuts <- if (is.null(written$n.cuts)) {
    controlMirror(control)$n.cuts
  } else {
    written$n.cuts
  }
  if (
    !is.null(cuts) && !anyNA(data@n.cuts) && !sameRecycled(cuts, data@n.cuts)
  ) {
    warnControlNotUsed(
      "n.cuts",
      cuts,
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

## $setControl's reading of a control's count or rule against the sampler's:
## nothing where they agree, a refusal where the caller wrote another, and
## otherwise the name, the sampler's being kept. A silent value keeps it with
## no message.
setControlCopyReading <- function(newControl, control, name) {
  given <- methods::slot(newControl, name)
  own <- methods::slot(control, name)
  if (identical(given, own)) {
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
