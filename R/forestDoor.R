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
## too. treeShift's validated value is the slot's tri-state, which forest()
## takes in words.
changedForestSettings <- function(given, inForce, label) {
  given <- given[!vapply(given, is.null, NA)]
  if (length(given) == 0L) {
    return(given)
  }
  given <- resolveForestMoves(validateForestKnobs(given), label)
  changed <- vapply(
    names(given),
    function(name) {
      value <- given[[name]]
      current <- inForce[[name]]
      switch(
        name,
        amplitude = !sameCoefficientPrior(
          value,
          validateForestAmplitude(current)
        ),
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

## How a door's refusal names forest `index`: by its label, as a creation names
## one by its list's name, and by position where the label is its position's.
doorForestLabel <- function(model, index) {
  labels <- vapply(model@forests, methods::slot, "", "label")
  labels[labels == paste0("forest", seq_along(labels))] <- ""
  forestKnobLabel(stats::setNames(as.list(labels), labels), index)
}

## The forests with no basis that state no sd, by position, a model's records
## read as resolveForests reads a call's statements.
unsizedForests <- function(sampler) {
  records <- sampler$model@forests
  bases <- sampler$data@bases
  which(vapply(
    seq_along(records),
    function(index) {
      is.null(bases[[index]]) && is.na(records[[index]]@sd)
    },
    NA
  ))
}

## The checks creation applies to a forest's statement that read more than the
## one name, with `changed` merged into what forest `index` states: a coefficient
## where creation refuses one and, here where a live sampler keeps its coefficient
## prior, any other; a size where the forest's is another forest's; the mixture
## a monotone sampler takes. Writes nothing.
checkForestStatement <- function(sampler, index, changed, label, ptr) {
  model <- sampler$model
  numForests <- bartcoreNumForests(ptr)
  if (samplerCarriesCounts(sampler)) {
    if (!is.null(changed$sd)) {
      refuseMultinomialStatement("a named leaf-prior 'sd'")
    }
    if (!is.null(changed$amplitude)) {
      stop(
        "a multinomial model's forests have no coefficient: its forests take ",
        "no 'amplitude'",
        call. = FALSE
      )
    }
  }
  amplitude <- changed$amplitude
  if (!is.null(amplitude)) {
    unsized <- if (numForests > 1L && is.null(changed$sd)) {
      unsizedForests(sampler)
    }
    refuseCoefficientPriorWithNoForest(
      amplitude,
      label,
      lone = numForests == 1L,
      sharing = length(unsized) > 1L && index %in% unsized
    )
    stop(
      "a forest's coefficient prior cannot be changed live yet: ",
      label,
      " keeps ",
      coefficientInWords(validatedCoefficient(sampler, index, ptr)),
      "; make a new sampler",
      call. = FALSE
    )
  }
  if (!is.null(changed$sd) && numForests > 1L) {
    sharers <- samplerSharers(sampler, ptr)
    if (length(sharers) > 1L && index %in% sharers) {
      stop(
        "this forest shares one prior total with forests ",
        paste(sharers, collapse = ", "),
        "; a forest's own size cannot be changed live yet: make a new sampler",
        call. = FALSE
      )
    }
  }
  if (!is.null(attr(model, "monotone"))) {
    refuseMonotoneForestMixture(changed$proposal.probs, label)
  }
  invisible(NULL)
}

## The coefficient's prior in force on a forest, as validateForestAmplitude
## returns one.
validatedCoefficient <- function(sampler, index, ptr) {
  validateForestAmplitude(forestCoefficientInForce(
    .Call(C_dbarts_bartcore_getLeafPrior, ptr, index - 1L)
  ))
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
## several its record, and where the last forest taking the fitting function's
## leaf prior states its own, that prior is left with no forest. On a model of
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
## able to refuse after a write, each undone if a later one refuses; the sd,
## being last, has nothing to undo. The records are written once the engine
## holds the values.
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
  inForce <- forestInForce(sampler, index, ptr)
  changed <- changedForestSettings(given, inForce, label)
  if (length(changed) == 0L) {
    return(invisible(NULL))
  }
  checkForestStatement(sampler, index, changed, label, ptr)

  written <- modelWithDepth(model, index, changed)
  movesChanged <- any(forestMoveNames %in% names(changed))
  if (movesChanged) {
    written <- modelWithMoves(written, index, changed)
  }
  if (!is.null(changed$sd)) {
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
      if (!is.null(changed$sd)) {
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
  sampler$model <- written
  if (resolveUpdateState(updateState, sampler$control)) {
    sampler$storeState(ptr)
  }
  invisible(NULL)
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
  if (length(changed) == 0L) {
    return(invisible(NULL))
  }
  refuseVarianceRuleGibbs(changed$proposal.probs)
  model <- sampler$model
  ptr <- sampler$getPointer()
  written <- model
  for (name in intersect(c("power", "base"), names(changed))) {
    methods::slot(written@variance, name) <- changed[[name]]
  }
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
