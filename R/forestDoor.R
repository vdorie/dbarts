## The door through which a live sampler changes one forest: $setForest and
## $getForest, $setVarianceForest and $getVarianceForest. A call changes only
## what it names, and the reader's function gives what a name is compared to, so
## a value equal to the one in force writes nothing and states nothing. Every
## check comes before the first write; the engine sees numbers, and the model's
## records and tree and leaf priors are written as a creation with the statement
## writes them, so a copy and a reload run what the live sampler runs.

## The names of a forest that say what the model is, each with what it is.
forestStructureNames <- c(
  vars = "the columns a forest may split on are",
  n.trees = "a forest's tree count is",
  interactions = "a forest's interaction limits are",
  blocks = "a forest's blocks are"
)

## What a call passed through its '...' is refused for, each name by what it
## says: structure whatever its value, a basis naming its own door, anything else
## as R refuses an argument. `taken` names what the method does take, `structure`
## what of forestStructureNames the forest has and `basis` whether it has one.
refuseForestDots <- function(supplied, method, taken, structure, basis = TRUE) {
  for (name in supplied) {
    if (name %in% structure) {
      stop(
        "'",
        name,
        "' is fixed at creation: ",
        forestStructureNames[[name]],
        " part of the model; make a new sampler",
        call. = FALSE
      )
    }
    if (basis && identical(name, "basis")) {
      stop(
        "'basis' is data, not a setting of the forest; change it with ",
        "$setForestBasis",
        call. = FALSE
      )
    }
    stop(
      "unused argument",
      if (nzchar(name)) paste0(" '", name, "'"),
      ": ",
      method,
      " takes ",
      paste(taken, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

## The names $setForest takes after the forest, and $setVarianceForest.
forestDoorNames <- c(
  "sd",
  "power",
  "base",
  "amplitude",
  "proposal.probs",
  "n.perturb.cuts",
  "treeShift"
)
varianceDoorNames <- c("power", "base", "proposal.probs", "n.perturb.cuts")

## The 1-based forest a door's 'forest' argument names: left out, the one forest
## of a sampler of one and, on a multinomial model, its one statement.
doorForestIndex <- function(sampler, forest, ptr) {
  numForests <- bartcoreNumForests(ptr)
  if (is.null(forest)) {
    if (numForests > 1L && !samplerCarriesCounts(sampler)) {
      stop(
        "'forest' must be given: this sampler has ",
        numForests,
        " forests; name one by its label or position",
        call. = FALSE
      )
    }
    return(1L)
  }
  index <- samplerForestIndex(forest, sampler$model, sampler$data, ptr) + 1L
  if (index > numForests) {
    stop("forest index out of range", call. = FALSE)
  }
  index
}

## A forest's depth as the bridge reads it: forest 1's from the model's tree
## prior, on every sampler, and the others' from their records. Either has the
## slots `base` and `power`.
forestDepthHolder <- function(model, index) {
  if (index == 1L) model@tree.prior else model@forests[[index]]
}

## A level-step slot's value in the words forest() takes it in.
treeShiftWords <- function(level) {
  if (is.na(level)) {
    "auto"
  } else if (level) {
    "always"
  } else {
    "never"
  }
}

## The move settings in force under a forest's own `moves`: each its own, else
## the control's.
movesInForce <- function(moves, control, mode = TRUE) {
  orControl <- function(own, slot) {
    if (is.null(own)) methods::slot(control, slot) else own
  }
  values <- list(
    proposal.probs = orControl(moves$proposal.probs, "proposal.probs"),
    n.perturb.cuts = orControl(moves$n.perturb.cuts, "n.perturb.cuts")
  )
  if (mode) {
    values$treeShift <- treeShiftWords(orControl(
      moves$levelGibbs,
      "levelGibbs"
    ))
  }
  values
}

## The refusal of a control a caller wrote a tree-move setting on that differs
## from the one in force when no forest would take it: every forest, the
## variance forest included, states its own, so $setControl would store a
## value nothing reads. NULL where some forest follows the control, or the
## setting is not the caller's (not named when the control was built, equal to
## a fresh control's) or equals the one in force.
setControlMovesRefusal <- function(newControl, control, model) {
  records <- model@forests
  if (!is.null(model@variance)) {
    records <- c(records, list(model@variance))
  }
  if (length(records) == 0L) {
    return(NULL)
  }
  fresh <- dbarts::dbartsControl()
  for (name in c("proposal.probs", "n.perturb.cuts")) {
    value <- methods::slot(newControl, name)
    written <- name %in%
      controlSuppliedSlots(newControl) ||
      !identical(value, methods::slot(fresh, name))
    ownEverywhere <- all(vapply(
      records,
      function(record) !is.null(record@moves[[name]]),
      NA
    ))
    if (
      written &&
        ownEverywhere &&
        !identical(value, methods::slot(control, name))
    ) {
      return(paste0(
        "'",
        name,
        "' on a control changes no forest: every forest states its own; ",
        "change one forest's with $setForest(forest, ",
        name,
        " = ) or the variance forest's with $setVarianceForest(",
        name,
        " = )"
      ))
    }
  }
  NULL
}

## The sd in force on a forest of a model that has one, NULL where there is no
## one positive finite number to give back: a drawn k, an sd under invchi(), a
## k of infinity and a unit that is not positive. A numeric forest's carries no
## column name here: it is a value for forest(sd = ).
forestSdInForce <- function(sampler, index, raw) {
  sd <- reportLeafPrior(sampler, raw, index - 1L)$sd
  if (is.na(sd) || !is.finite(sd) || sd <= 0) NULL else unname(sd)
}

## The coefficient's prior in force, from the engine's row for the forest. A
## forest with a basis stating none draws under its default, which is normal().
forestCoefficientInForce <- function(raw) {
  if (is.nan(raw[[1L, "map.anchor"]])) {
    return(fixed())
  }
  held <- raw[[1L, "amplitude.held"]] != 0 || !is.nan(raw[[1L, "share"]])
  prior <- statedCoefficientPrior(raw[[1L, "amplitude.df"]], held)
  if (is.null(prior)) normal() else prior
}

## Every value in force on forest `index`, in the names $setForest takes on its
## kind and in that order: what $getForest returns and step 2 of a call compares
## to. A multinomial model's forests take one statement, read from its one
## record.
forestInForce <- function(sampler, index, ptr) {
  force(index)
  model <- sampler$model
  multinomial <- samplerCarriesCounts(sampler)
  holder <- forestDepthHolder(model, if (multinomial) 1L else index)
  values <- list()
  raw <- NULL
  if (!multinomial) {
    raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, index - 1L)
    values["sd"] <- list(forestSdInForce(sampler, index, raw))
  }
  values$power <- holder@power
  values$base <- holder@base
  if (!multinomial) {
    values$amplitude <- forestCoefficientInForce(raw)
  }
  c(
    values,
    movesInForce(
      model@forests[[if (multinomial) 1L else index]]@moves,
      sampler$control
    )
  )
}

## Every value in force on the variance forest, NULL where the sampler has none.
varianceForestInForce <- function(sampler) {
  record <- sampler$model@variance
  if (is.null(record)) {
    return(NULL)
  }
  c(
    list(power = record@power, base = record@base),
    movesInForce(record@moves, sampler$control, FALSE)
  )
}

## The names of `given`, a door's arguments with those NULL gone, each validated
## as creation validates it, and with those equal to the value in force dropped
## too; a name with no value in force is not equal to one. Each is compared with
## the value in force before the call, whatever else the call changes. treeShift's
## validated value is the slot's tri-state, which forest() takes in words.
changedForestSettings <- function(given, inForce, label) {
  given <- given[!vapply(given, is.null, NA)]
  if (length(given) == 0L) {
    return(given)
  }
  given <- resolveForestMoves(validateForestKnobs(given, label), label)
  amplitudeChanges <- !is.null(given$amplitude) &&
    !sameCoefficientPrior(
      given$amplitude,
      validateForestAmplitude(inForce$amplitude, label)
    )
  changed <- vapply(
    names(given),
    function(name) {
      value <- given[[name]]
      current <- inForce[[name]]
      switch(
        name,
        amplitude = amplitudeChanges,
        treeShift = !identical(value, resolveTreeShift(current)),
        !identical(value, current)
      )
    },
    NA
  )
  given[changed]
}

## Whether two coefficient priors, each as validateForestAmplitude returns it,
## are one: a df is compared as a number.
sameCoefficientPrior <- function(a, b) {
  if (is(a, "dbartsStudentPrior") && is(b, "dbartsStudentPrior")) {
    return(isTRUE(as.double(a@df) == as.double(b@df)))
  }
  is(a, "dbartsNormalPrior") && is(b, "dbartsNormalPrior") || identical(a, b)
}

## How a door's refusal names forest `index`: by its label, as the fit's
## results do, and by position where it has none of its own.
doorForestLabel <- function(model, index) {
  forestKnobLabel(vapply(model@forests, methods::slot, "", "label"), index)
}

## Whether forest `index` takes its size from the fitting function's leaf
## prior: one with no basis that states no sd, its record read as
## resolveForests reads a call's statement.
takesLeafPriorSize <- function(sampler, index) {
  is.null(sampler$data@bases[[index]]) &&
    is.na(sampler$model@forests[[index]]@sd)
}

## The checks creation applies to a forest's statement that read more than the
## one name, with `changed` merged into what forest `index` states: a drawn
## coefficient where creation refuses one, on a model of one forest, on a
## forest sharing the prior total with others and stating no 'sd', and on a
## forest sized by the fitting function's leaf prior where that gives no size;
## the mixture a monotone sampler takes. Each forest's place is read from the
## engine's row and the records, never from a count of forests. Writes nothing.
checkForestStatement <- function(sampler, index, changed, label, ptr) {
  model <- sampler$model
  numForests <- bartcoreNumForests(ptr)
  if (samplerCarriesCounts(sampler)) {
    if (!is.null(changed$sd)) {
      refuseMultinomialStatement("a named leaf-prior 'sd'")
    }
  }
  amplitude <- changed$amplitude
  if (!is.null(amplitude)) {
    sharers <- if (numForests > 1L && is.null(changed$sd)) {
      samplerSharers(sampler, ptr)
    }
    refuseCoefficientPriorWithNoForest(
      amplitude,
      label,
      lone = numForests == 1L,
      sharing = length(sharers) > 1L && index %in% sharers
    )
    if (
      statesCoefficientPrior(amplitude) &&
        is.null(changed$sd) &&
        takesLeafPriorSize(sampler, index) &&
        leafPriorGivesNoSize(model@leaf.hyperprior)
    ) {
      refuseSizelessCoefficient(label, amplitude, model@leaf.hyperprior)
    }
  }
  if (!is.null(attr(model, "monotone"))) {
    refuseMonotoneForestMixture(changed$proposal.probs, label)
  }
  invisible(NULL)
}

## A coefficient's prior as validateForestAmplitude returns it, as a call writes
## it, for a message.
coefficientInWords <- function(amplitude) {
  if (identical(amplitude, "fixed")) {
    "fixed()"
  } else if (is(amplitude, "dbartsStudentPrior")) {
    paste0("student(", format(amplitude@df), ")")
  } else {
    "normal()"
  }
}

## Which forests a depth, by a door, is installed on: a multinomial model's are
## one statement.
depthInstallForests <- function(sampler, index, ptr) {
  if (samplerCarriesCounts(sampler)) {
    seq_len(bartcoreNumForests(ptr))
  } else {
    index
  }
}

## The names of a statement that are move settings.
forestMoveNames <- c("proposal.probs", "n.perturb.cuts", "treeShift")

## The move settings `moves`, a forest's own list, with `changed` written over
## them, as forest()'s names, the mode as its slot; a spec for forestMoveRecord.
moveSpecWith <- function(moves, changed = NULL) {
  spec <- list(
    proposal.probs = moves$proposal.probs,
    n.perturb.cuts = moves$n.perturb.cuts,
    treeShift = moves$levelGibbs
  )
  for (name in intersect(names(changed), forestMoveNames)) {
    spec[name] <- list(changed[[name]])
  }
  spec
}

## A model with `changed` written on the move record of forest `index`, or with
## `index` NULL on the variance forest's, through the function creation writes
## the record with.
modelWithMoves <- function(model, index, changed) {
  records <- model@forests
  specs <- lapply(seq_along(records), function(position) {
    moveSpecWith(
      records[[position]]@moves,
      if (identical(position, index)) changed
    )
  })
  variance <- moveSpecWith(
    if (!is.null(model@variance)) model@variance@moves,
    if (is.null(index)) changed
  )
  moves <- forestMoveRecord(specs, length(records), variance)
  for (position in seq_along(records)) {
    records[[position]]@moves <- moves$forests[[position]]
  }
  model@forests <- records
  if (!is.null(model@variance)) {
    model@variance@moves <- moves$variance
  }
  model
}

## A model with a changed depth written where creation writes it: forest 1's in
## the tree prior on every sampler and, on a model of several, every forest's in
## its record, with `stated` gaining what was written.
modelWithDepth <- function(model, index, changed) {
  for (name in intersect(c("power", "base"), names(changed))) {
    if (index == 1L) {
      methods::slot(model@tree.prior, name) <- changed[[name]]
    }
    if (length(model@forests) > 1L) {
      methods::slot(model@forests[[index]], name) <- changed[[name]]
      model@forests[[index]]@stated <- withDepthStated(
        model@forests[[index]]@stated,
        name
      )
    }
  }
  model
}

## What a forest's record says it stated, with `name` among it, in creation's
## order.
withDepthStated <- function(stated, name) {
  order <- c("n.trees", "power", "base")
  all <- union(stated, name)
  c(intersect(order, all), setdiff(all, order))
}

## A model with forest `index`'s sd written as creation writes it: on a model of
## several its record, its share of the prior total gone, so that a sharer
## leaves the group and the others keep their parts, and where the last forest
## taking the fitting function's leaf prior states its own, that prior is left
## with no forest. On a model of
## one, the leaf prior itself, in whichever leaf model the sampler has, as
## $setLeafPrior writes it.
modelWithForestSd <- function(sampler, model, index, sd) {
  if (length(model@forests) == 1L) {
    spec <- model@leaf.prior
    spec@k <- NULL
    spec@prior.sd <- sd
    written <- restateLeafPrior(sampler, spec, NULL)
    model@leaf.prior <- written@leaf.prior
    model@leaf.hyperprior <- written@leaf.hyperprior
    model@prior.scale <- written@prior.scale
    return(model)
  }
  model@forests[[index]]@sd <- sd
  model@forests[[index]]@share <- NA_real_
  hasBasis <- !vapply(sampler$data@bases, is.null, NA)
  recordedSd <- vapply(model@forests, methods::slot, 0, "sd")
  if (!any(!hasBasis & is.na(recordedSd))) {
    model <- modelWithLeafPriorUnclaimed(model)
  }
  model
}

## A model with forest `index`'s coefficient `amplitude`, as
## validateForestAmplitude returns it, written as creation writes it
## (forestStatements), with `sd`, NULL where none, written first as
## modelWithForestSd writes it. A forest the fitting function's leaf prior
## sizes stops sharing the total when drawn, and shares it alone when held;
## the last left of a group, drawn, states an sd the engine gives it, which
## modelWithStatedSd records once the engine holds it.
modelWithCoefficient <- function(sampler, model, index, amplitude, sd) {
  if (!is.null(sd)) {
    model <- modelWithForestSd(sampler, model, index, sd)
  }
  hasBasis <- !is.null(sampler$data@bases[[index]])
  record <- model@forests[[index]]
  record@held.amplitude <- identical(amplitude, "fixed")
  record@amplitude.df <- if (is(amplitude, "dbartsStudentPrior")) {
    as.double(amplitude@df)
  } else if (is(amplitude, "dbartsNormalPrior") && !hasBasis) {
    Inf
  } else {
    NA_real_
  }
  if (!hasBasis && is.na(record@sd)) {
    record@share <- if (record@held.amplitude) 1 else NA_real_
  }
  model@forests[[index]] <- record
  model
}

## A model with the sd the engine states on forest `index` recorded where its
## record states none, as the last left of a group released states its part:
## read back from the engine's row, so no number of R's is compared with one
## of the engine's. Where no forest then takes the fitting function's leaf
## prior, that prior is left with no forest, as modelWithForestSd leaves it.
modelWithStatedSd <- function(sampler, model, index, ptr) {
  if (
    !takesLeafPriorSize(sampler, index) ||
      model@forests[[index]]@held.amplitude
  ) {
    return(model)
  }
  raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, index - 1L)
  stated <- raw[[1L, "named.sd"]]
  if (is.na(stated) || !is.na(model@forests[[index]]@sd)) {
    return(model)
  }
  model@forests[[index]]@sd <- stated
  model@forests[[index]]@share <- NA_real_
  hasBasis <- !vapply(sampler$data@bases, is.null, NA)
  recordedSd <- vapply(model@forests, methods::slot, 0, "sd")
  if (!any(!hasBasis & is.na(recordedSd))) {
    model <- modelWithLeafPriorUnclaimed(model)
  }
  model
}

## Installs a forest's depth on the forests given, 0-based, -1L the variance
## forest, from `holder`'s slots.
installForestDepth <- function(ptr, forests, holder) {
  for (forest in forests) {
    .Call(
      C_dbarts_bartcore_setForestTreePrior,
      ptr,
      forest,
      holder@base,
      holder@power
    )
  }
  invisible(NULL)
}

## $setForest. The checks, then the engine in the order that leaves only the sd
## or the coefficient able to refuse after a write, each undone if a later one
## refuses; the last, the coefficient with any sd beside it or else the sd, has
## nothing to undo. The records are written once the engine holds the values.
writeForestDoor <- function(sampler, forest, given, supplied, updateState) {
  updateState <- checkUpdateState(updateState)
  multinomial <- samplerCarriesCounts(sampler)
  refuseForestDots(
    supplied,
    "$setForest",
    c("forest", forestDoorNames, "updateState"),
    names(forestStructureNames)
  )
  if (multinomial && !is.null(forest)) {
    stop(
      "a multinomial model's forests take one statement, so this would change ",
      "every category; write s$setForest(power = 4)",
      call. = FALSE
    )
  }
  ptr <- sampler$getPointer()
  index <- doorForestIndex(sampler, forest, ptr)
  model <- sampler$model
  label <- if (multinomial) "forest 1" else doorForestLabel(model, index)
  if (multinomial && !is.null(given$amplitude)) {
    # no category has a coefficient to draw: fixed() states what is so
    refuseCoefficientPriorWithNoForest(
      validateForestAmplitude(given$amplitude, label),
      label,
      lone = TRUE,
      sharing = FALSE
    )
    given$amplitude <- NULL
  }
  inForce <- forestInForce(sampler, index, ptr)
  changed <- changedForestSettings(given, inForce, label)
  if (length(changed) == 0L) {
    if (resolveUpdateState(updateState, sampler$control)) {
      sampler$storeState(ptr)
    }
    return(invisible(NULL))
  }
  checkForestStatement(sampler, index, changed, label, ptr)

  written <- modelWithDepth(model, index, changed)
  movesChanged <- any(forestMoveNames %in% names(changed))
  if (movesChanged) {
    written <- modelWithMoves(written, index, changed)
  }
  coefficientChanged <- !is.null(changed$amplitude)
  if (coefficientChanged) {
    written <- modelWithCoefficient(
      sampler,
      written,
      index,
      changed$amplitude,
      changed$sd
    )
  } else if (!is.null(changed$sd)) {
    written <- modelWithForestSd(sampler, written, index, changed$sd)
  }

  depthChanged <- any(c("power", "base") %in% names(changed))
  depthForests <- depthInstallForests(sampler, index, ptr) - 1L
  undo <- function() {
    if (depthChanged) {
      installForestDepth(
        ptr,
        depthForests,
        forestDepthHolder(model, if (multinomial) 1L else index)
      )
    }
    if (movesChanged) {
      bartcoreSetMoveSettings(ptr, sampler$control, model)
    }
  }
  tryCatch(
    {
      if (movesChanged) {
        bartcoreSetMoveSettings(ptr, sampler$control, written)
      }
      if (depthChanged) {
        installForestDepth(
          ptr,
          depthForests,
          forestDepthHolder(written, if (multinomial) 1L else index)
        )
      }
      if (coefficientChanged) {
        # one engine call for the coefficient and an sd beside it, so every
        # chain is checked for both before either is written
        record <- written@forests[[index]]
        .Call(
          C_dbarts_bartcore_setForestCoefficient,
          ptr,
          index - 1L,
          record@held.amplitude,
          record@amplitude.df,
          if (is.null(changed$sd)) NA_real_ else as.double(changed$sd)
        )
      } else if (!is.null(changed$sd)) {
        if (length(model@forests) == 1L) {
          writeLeafPrior(sampler, ptr, written)
        } else {
          .Call(C_dbarts_bartcore_setForestSd, ptr, index - 1L, changed$sd)
        }
      }
    },
    error = function(e) {
      undo()
      stop(e)
    }
  )
  if (coefficientChanged) {
    written <- modelWithStatedSd(sampler, written, index, ptr)
  }
  sampler$model <- written
  if (resolveUpdateState(updateState, sampler$control)) {
    sampler$storeState(ptr)
  }
  invisible(NULL)
}

## A model with the variance forest's depth in `changed` written on its record,
## each name written stated.
modelWithVarianceDepth <- function(model, changed) {
  for (name in intersect(c("power", "base"), names(changed))) {
    methods::slot(model@variance, name) <- changed[[name]]
    model@variance@stated <- withDepthStated(model@variance@stated, name)
  }
  model
}

## The variance forest's depth a model given to $setModel states and that
## differs from the sampler's, validated as $setVarianceForest validates it;
## an empty list where none does. A depth the record does not state, as one
## copied from the mean forest's at creation, is left.
statedVarianceDepth <- function(sampler, given) {
  record <- given@variance
  if (is.null(record) || is.null(sampler$model@variance)) {
    return(list())
  }
  names <- intersect(c("power", "base"), record@stated)
  stated <- stats::setNames(
    lapply(names, function(name) methods::slot(record, name)),
    names
  )
  changedForestSettings(
    stated,
    varianceForestInForce(sampler),
    "varianceForest()"
  )
}

## $setVarianceForest, on the same steps; it has no sd to refuse after a write.
writeVarianceForestDoor <- function(sampler, given, supplied, updateState) {
  updateState <- checkUpdateState(updateState)
  refuseForestDots(
    supplied,
    "$setVarianceForest",
    c(varianceDoorNames, "updateState"),
    c("vars", "n.trees"),
    basis = FALSE
  )
  inForce <- varianceForestInForce(sampler)
  if (is.null(inForce)) {
    stop(
      "this sampler has no variance forest: whether there is one is part of ",
      "the model; make a new sampler",
      call. = FALSE
    )
  }
  changed <- changedForestSettings(given, inForce, "varianceForest()")
  ptr <- sampler$getPointer()
  if (length(changed) == 0L) {
    if (resolveUpdateState(updateState, sampler$control)) {
      sampler$storeState(ptr)
    }
    return(invisible(NULL))
  }
  refuseVarianceRuleGibbs(changed$proposal.probs)
  written <- modelWithVarianceDepth(sampler$model, changed)
  movesChanged <- any(forestMoveNames %in% names(changed))
  if (movesChanged) {
    written <- modelWithMoves(written, NULL, changed)
    bartcoreSetMoveSettings(ptr, sampler$control, written)
  }
  if (any(c("power", "base") %in% names(changed))) {
    installForestDepth(ptr, -1L, written@variance)
  }
  sampler$model <- written
  if (resolveUpdateState(updateState, sampler$control)) {
    sampler$storeState(ptr)
  }
  invisible(NULL)
}
