## The fitting function's level on a live sampler of several forests or a
## multinomial one: the depth and the leaf prior that $setTreePrior,
## $setLeafPrior and $setModel write. Nothing records the level apart from the
## forests that run it, so its value in force is that of the first forest it
## reaches. A door's candidates are the forests with no basis: the one
## candidate, or a multinomial model's one entry, is written whether it states
## its own or not; of several, those stating none of their own are, and where
## all state their own the door is refused, naming $setForest. A value equal to
## the one in force writes nothing, and a value written makes no forest state
## it. The variance forest follows nothing and stays as it is.

## The 1-based forests a fitting-function door can reach: those with no basis,
## and a multinomial model's one entry.
levelCandidates <- function(sampler) {
  if (samplerCarriesCounts(sampler)) {
    return(1L)
  }
  which(vapply(sampler$data@bases, is.null, NA))
}

## The forests of `candidates` a depth name reaches on the sampler's model:
## the one candidate whatever it states, else those not stating `name`.
levelReach <- function(model, candidates, name) {
  if (length(candidates) <= 1L) {
    return(candidates)
  }
  stated <- vapply(
    model@forests[candidates],
    function(record) name %in% record@stated,
    NA
  )
  candidates[!stated]
}

## Forest `index`'s depth value `name` in `model`: forest 1's from the tree
## prior, the others' from their records. A record past what `model` holds, as
## on a model made by hand, gives `fallback`'s.
forestDepthValue <- function(model, index, name, fallback = model) {
  if (index == 1L) {
    return(methods::slot(model@tree.prior, name))
  }
  if (index > length(model@forests)) {
    return(forestDepthValue(fallback, index, name))
  }
  methods::slot(model@forests[[index]], name)
}

## The call a refusal names for forest `index`: $setForest with its label, or
## its position where the label is its position's.
setForestCall <- function(model, index, name, value) {
  labels <- vapply(model@forests, methods::slot, "", "label")
  label <- labels[index]
  forest <- if (
    length(label) == 1L &&
      nzchar(label) &&
      !identical(label, paste0("forest", index))
  ) {
    paste0("\"", label, "\"")
  } else {
    index
  }
  paste0("$setForest(", forest, ", ", name, " = ", format(value), ")")
}

## The refusal of a depth `name` given to `door` that reaches no forest, naming
## the forests that state it and the call to write.
refuseLevelReachesNone <- function(model, candidates, name, value, door) {
  if (length(candidates) == 0L) {
    stop(
      "'",
      name,
      "' given to ",
      door,
      " reaches no forest: every forest of this sampler has a basis; change ",
      "one with ",
      setForestCall(model, 1L, name, value),
      call. = FALSE
    )
  }
  labels <- vapply(candidates, doorForestLabel, "", model = model)
  stop(
    "'",
    name,
    "' given to ",
    door,
    " reaches no forest: ",
    paste(labels, collapse = " and "),
    if (length(candidates) == 1L) {
      ", the forest with no basis, states its own"
    } else {
      ", the forests with no basis, state their own"
    },
    "; change ",
    if (length(candidates) == 1L) "it" else "one",
    " with ",
    setForestCall(model, candidates[[1L]], name, value),
    call. = FALSE
  )
}

## The refusal of a difference in forest `index`'s `name`, which its door
## changes.
refuseOwnDifference <- function(model, index, name, value) {
  label <- doorForestLabel(model, index)
  record <- model@forests[[index]]
  states <- if (name == "sd") {
    !is.na(record@sd)
  } else {
    name %in% record@stated
  }
  stop(
    "$setModel cannot change ",
    label,
    "'s '",
    name,
    "'",
    if (states) ", which the forest states as its own",
    ": change it with ",
    setForestCall(model, index, name, value),
    call. = FALSE
  )
}

## The depths a model given to $setModel states for a sampler of several
## forests or a multinomial one, as a list of the model to store and the bases
## and powers to install, one per forest. Per name, the level first, as
## $setTreePrior writes it: its value is the given one of the first forest it
## reaches that the given model does not state. Then each forest whose given
## record states the name at another value than it holds once the level is
## written takes it, as $setForest after $setTreePrior writes it. A number the
## given record does not state is left.
levelDepths <- function(sampler, given, stored) {
  model <- sampler$model
  numForests <- bartcoreNumForests(sampler$getPointer())
  multinomial <- samplerCarriesCounts(sampler)
  candidates <- levelCandidates(sampler)
  givenRecords <- length(given@forests) == length(model@forests)
  for (name in c("power", "base")) {
    if (multinomial) {
      methods::slot(stored@tree.prior, name) <- methods::slot(
        given@tree.prior,
        name
      )
      next
    }
    # forest 1's depth, in the tree prior on every sampler, starts at its own
    methods::slot(stored@tree.prior, name) <- methods::slot(
      model@tree.prior,
      name
    )
    reach <- levelReach(model, candidates, name)
    own <- vapply(
      seq_len(numForests),
      forestDepthValue,
      0,
      model = model,
      name = name
    )
    stated <- vapply(
      seq_len(numForests),
      forestDepthValue,
      0,
      model = given,
      name = name,
      fallback = model
    )
    givenStates <- vapply(
      seq_len(numForests),
      function(index) givenRecords && name %in% given@forests[[index]]@stated,
      NA
    )
    sources <- reach[!givenStates[reach]]
    if (length(sources) > 0L) {
      level <- stated[[sources[[1L]]]]
      for (index in sources[-1L]) {
        if (
          !identical(stated[[index]], own[[index]]) &&
            !identical(stated[[index]], level)
        ) {
          refuseOwnDifference(model, index, name, stated[[index]])
        }
      }
      if (!identical(level, own[[sources[[1L]]]])) {
        stored <- withLevelDepth(stored, reach, name, level)
      }
    }
    for (index in which(givenStates)) {
      value <- validateForestTreePrior(stated[[index]], name)
      if (!identical(value, forestDepthValue(stored, index, name))) {
        stored <- modelWithDepth(
          stored,
          index,
          stats::setNames(list(value), name)
        )
      }
    }
  }
  depth <- function(name) {
    if (multinomial) {
      return(rep(methods::slot(stored@tree.prior, name), numForests))
    }
    vapply(
      seq_len(numForests),
      forestDepthValue,
      0,
      model = stored,
      name = name
    )
  }
  list(model = stored, bases = depth("base"), powers = depth("power"))
}

## `model` with depth `name` at `value` on the forests `reach`, where creation
## writes it: forest 1's in the tree prior and, on a model of several, every
## forest's in its record, what each states left as it is.
withLevelDepth <- function(model, reach, name, value) {
  for (index in reach) {
    if (index == 1L) {
      methods::slot(model@tree.prior, name) <- value
    }
    if (length(model@forests) > 1L) {
      methods::slot(model@forests[[index]], name) <- value
    }
  }
  model
}

## Whether two models state one fitting-function leaf prior: its resolved
## hyperprior and named scale.
sameLeafPrior <- function(a, b) {
  identical(a@leaf.hyperprior, b@leaf.hyperprior) &&
    identical(a@prior.scale, b@prior.scale)
}

## The model to store and whether to install the fitting function's leaf prior
## `given` states, on a sampler of several forests or a multinomial one: the
## one writer behind $setLeafPrior and $setModel. `stored` is the model being
## written, carrying the sampler's records; `givenSds` the sd each forest's
## record in the given model states, NA where none, which matters here only on
## the one candidate: the leaf prior reaches it whatever it states, unless the
## given model states its sd. Where no forest takes the leaf prior, anything
## but the one in force is refused.
levelLeafPrior <- function(sampler, given, stored, givenSds, door) {
  model <- sampler$model
  if (samplerCarriesCounts(sampler)) {
    if (sameLeafPrior(given, model)) {
      return(list(model = stored, install = FALSE))
    }
    refuseMultinomialLeafPrior(given@leaf.prior)
    if (is(given@leaf.hyperprior, "dbartsChiHyperprior")) {
      stop(
        multinomialLeafPriorMessage,
        "a 'k' hyperprior is not supported on a multinomial sampler, at ",
        "creation or after",
        call. = FALSE
      )
    }
    return(list(model = withLeafPrior(stored, given), install = TRUE))
  }
  candidates <- levelCandidates(sampler)
  ownSds <- vapply(model@forests, methods::slot, 0, "sd")
  changes <- !sameLeafPrior(given, model)
  if (length(candidates) == 1L) {
    claims <- is.na(givenSds[[candidates]]) && !is.na(ownSds[[candidates]])
    if (!changes && !claims) {
      return(list(model = stored, install = FALSE))
    }
    stored <- withLeafPrior(stored, given)
    stored@forests[[candidates]]@sd <- NA_real_
    # it then shares the total alone, unless it draws its own coefficient
    if (is.na(stored@forests[[candidates]]@amplitude.df)) {
      stored@forests[[candidates]]@share <- 1
    }
    return(list(model = stored, install = TRUE))
  }
  if (!changes) {
    return(list(model = stored, install = FALSE))
  }
  if (length(samplerSharers(sampler)) == 0L) {
    named <- if (length(candidates) > 0L) candidates[[1L]] else 1L
    statesDf <- any(
      !is.na(vapply(
        model@forests[candidates],
        methods::slot,
        0,
        "amplitude.df"
      ))
    )
    stop(
      door,
      ": no forest of this sampler takes the fitting function's leaf prior, ",
      if (length(candidates) == 0L) {
        "every forest having a basis"
      } else if (statesDf) {
        "every forest with no basis stating its own sd or a coefficient's prior"
      } else {
        "every forest with no basis stating its own sd"
      },
      ", so it is the default in force; state a forest's size with ",
      setForestCall(model, named, "sd", "v"),
      call. = FALSE
    )
  }
  list(model = withLeafPrior(stored, given), install = TRUE)
}

## The refusal of a forest's share of the prior total that differs in `given`
## from `planned`'s, the model $setModel would store: the parts are set once, at
## creation, and a forest stating its sd leaves the group. The one candidate is
## not read, the leaf prior reaching it whatever it states.
refuseShareDifference <- function(sampler, given, planned) {
  model <- planned
  candidates <- levelCandidates(sampler)
  if (
    length(given@forests) != length(model@forests) ||
      samplerCarriesCounts(sampler)
  ) {
    return(invisible(NULL))
  }
  for (index in seq_along(model@forests)) {
    if (length(candidates) == 1L && index %in% candidates) {
      next
    }
    if (
      !identical(given@forests[[index]]@share, model@forests[[index]]@share)
    ) {
      stop(
        "$setModel cannot give ",
        doorForestLabel(model, index),
        " another part of the prior total (",
        partFraction(model@forests[[index]]@share),
        " here, ",
        partFraction(given@forests[[index]]@share),
        " in the model given): the parts are set when a sampler is created, ",
        "from which forests share it; make a new sampler instead",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

## The refusal of a forest's coefficient that differs in `given` from the
## sampler's, held or drawn or under another prior: it is the forest's own,
## which $setForest changes, named in the value the given model states.
refuseCoefficientDifference <- function(sampler, given) {
  model <- sampler$model
  if (
    length(given@forests) != length(model@forests) ||
      samplerCarriesCounts(sampler)
  ) {
    return(invisible(NULL))
  }
  for (index in seq_along(model@forests)) {
    own <- model@forests[[index]]
    stated <- given@forests[[index]]
    if (
      identical(stated@held.amplitude, own@held.amplitude) &&
        identical(stated@amplitude.df, own@amplitude.df)
    ) {
      next
    }
    prior <- statedCoefficientPrior(stated@amplitude.df, stated@held.amplitude)
    if (is.null(prior)) {
      prior <- if (is.null(sampler$data@bases[[index]])) fixed() else normal()
    }
    refuseOwnDifference(
      model,
      index,
      "amplitude",
      coefficientInWords(validateForestAmplitude(prior))
    )
  }
  invisible(NULL)
}

## `model` with the fitting function's leaf prior `given` states.
withLeafPrior <- function(model, given) {
  model@leaf.prior <- given@leaf.prior
  model@leaf.hyperprior <- given@leaf.hyperprior
  model@prior.scale <- given@prior.scale
  model
}

## The refusal of split probabilities other than equal ones given to `door`,
## in creation's words.
refuseLevelSplitProbabilities <- function(sampler, door) {
  stop(
    door,
    ": a ",
    if (samplerCarriesCounts(sampler)) {
      "multinomial (softmax) model"
    } else {
      "model of several forests"
    },
    " does not support 'split.probs' (a vector that is not uniform, or ",
    "dart())",
    call. = FALSE
  )
}

## $setModel on a sampler of several forests or a multinomial one: the level's
## depth and leaf prior, the residual prior and the move settings, through one
## bridge entry that never takes the single forest's install. `newModel` has
## passed the structural checks and carries the sampler's records
## (keepModelRecords); `given` is the model as the caller wrote it. Returns the
## model to store.
setModelOnForests <- function(sampler, newModel, given) {
  model <- sampler$model
  if (varySplitProbabilities(given@tree.prior)) {
    refuseLevelSplitProbabilities(sampler, "$setModel")
  }
  if (!identical(given@leaf.scale, model@leaf.scale)) {
    stop(
      "$setModel cannot change the leaf scale of a sampler of several ",
      "forests: the calibration map sets every forest's",
      call. = FALSE
    )
  }
  depths <- levelDepths(sampler, given, newModel)
  # a model recording no forest, as one made by hand, states no sd
  givenRecords <- length(given@forests) == length(model@forests)
  ownSds <- vapply(model@forests, methods::slot, 0, "sd")
  statedSds <- if (givenRecords) {
    vapply(given@forests, methods::slot, 0, "sd")
  } else {
    rep(NA_real_, length(ownSds))
  }
  givenSds <- if (givenRecords) statedSds else ownSds
  leaf <- levelLeafPrior(sampler, given, depths$model, givenSds, "$setModel")
  refuseCoefficientDifference(sampler, given)
  # the forests whose size the level's leaf prior writes: the one candidate
  # whatever it states, and a forest it sizes
  candidates <- levelCandidates(sampler)
  levelSized <- vapply(
    seq_along(ownSds),
    function(index) {
      leaf$install &&
        ((length(candidates) == 1L && index %in% candidates) ||
          takesLeafPriorSize(sampler, index))
    },
    NA
  )
  sds <- statedForestSds(sampler, leaf$model, statedSds, ownSds, levelSized)
  stored <- sds$model
  refuseShareDifference(sampler, given, stored)
  stored@tree.prior@splitProbabilities <- numeric()
  ptr <- sampler$getPointer()
  withSizeRefusal(
    sampler,
    "$setModel",
    .Call(
      C_dbarts_bartcore_setModelForests,
      ptr,
      stored,
      sampler$data,
      sampler$control,
      depths$bases,
      depths$powers,
      leaf$install
    )
  )
  # checked above, so none of these refuses
  for (index in names(sds$installs)) {
    .Call(
      C_dbarts_bartcore_setForestSd,
      ptr,
      as.integer(index) - 1L,
      sds$installs[[index]]
    )
  }
  stored
}

## The sds `givenSds` a model given to $setModel states, each written on
## `stored`, the model after the level's writes, as $setForest(f, sd = ) writes
## it, as a list of the model and the sds to install by forest position. One is
## compared, as the door compares, with the sd in force once the level is
## written: where the level's leaf prior does not size the forest
## (`levelSized`), the sd in force now; where it does, the size it gives, which
## a statement the sampler does not already make is taken to differ from. A
## forest whose given record states no sd, or the sd its sampler's record
## states, is left to the level.
statedForestSds <- function(sampler, stored, givenSds, ownSds, levelSized) {
  installs <- list()
  if (samplerCarriesCounts(sampler)) {
    return(list(model = stored, installs = installs))
  }
  ptr <- sampler$getPointer()
  for (index in which(!is.na(givenSds))) {
    sd <- validateForestSd(givenSds[[index]])
    if (identical(sd, unname(ownSds[[index]]))) {
      next
    }
    if (!levelSized[[index]]) {
      raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, index - 1L)
      if (identical(sd, forestSdInForce(sampler, index, raw))) {
        next
      }
    }
    stored <- modelWithForestSd(sampler, stored, index, sd)
    installs[[as.character(index)]] <- sd
  }
  list(model = stored, installs = installs)
}

## Evaluates `expr`, an install of the fitting function's leaf prior. Where one
## forest has no basis the engine finds no forest to take it only for a size
## the prior cannot give a forest stating a coefficient's prior, which is said
## by name.
withSizeRefusal <- function(sampler, door, expr) {
  candidates <- levelCandidates(sampler)
  tryCatch(expr, error = function(e) {
    taken <- grepl(
      "no forest of this sampler takes the fitting function's leaf prior",
      conditionMessage(e),
      fixed = TRUE
    )
    if (!taken || length(candidates) != 1L || samplerCarriesCounts(sampler)) {
      stop(e)
    }
    stop(
      door,
      ": ",
      doorForestLabel(sampler$model, candidates),
      " states a coefficient's prior and takes its size from the fitting ",
      "function's leaf prior, which here gives no positive finite size; state ",
      "the forest's size with ",
      setForestCall(sampler$model, candidates, "sd", "v"),
      call. = FALSE
    )
  })
}

## $setLeafPrior on a sampler of several forests or a multinomial one: the
## restated leaf prior through levelLeafPrior, the writer $setModel shares,
## the forests' own sd unchanged. A bare normal() where no forest takes the
## leaf prior is the default restated. Returns the model to store, or NULL
## where nothing changes.
setLeafPriorOnForests <- function(sampler, spec, expr) {
  model <- sampler$model
  if (
    !samplerCarriesCounts(sampler) &&
      length(levelCandidates(sampler)) != 1L &&
      length(samplerSharers(sampler)) == 0L &&
      is(spec, "dbartsNormalPrior") &&
      is.null(spec@k) &&
      is.null(spec@prior.sd)
  ) {
    return(NULL)
  }
  if (samplerCarriesCounts(sampler)) {
    refuseMultinomialLeafPrior(spec)
  }
  given <- restateLeafPrior(sampler, spec, expr)
  # the forests' own sd as they stand, the one candidate's excepted: it takes
  # the leaf prior whatever it states
  givenSds <- vapply(model@forests, methods::slot, 0, "sd")
  candidates <- levelCandidates(sampler)
  if (!samplerCarriesCounts(sampler) && length(candidates) == 1L) {
    givenSds[candidates] <- NA_real_
  }
  leaf <- levelLeafPrior(sampler, given, model, givenSds, "$setLeafPrior")
  if (!leaf$install) {
    return(NULL)
  }
  withSizeRefusal(
    sampler,
    "$setLeafPrior",
    .Call(
      C_dbarts_bartcore_setFittingLeafPrior,
      sampler$getPointer(),
      leaf$model
    )
  )
  leaf$model
}

## $setTreePrior: the stored model with the level changed, through $setModel's
## route. `given` holds power and base, NULL where left out; `splitExpr` the
## split.probs expression, NULL where left out, evaluated in `env`, a NULL
## value left out too. A DART setting left out takes creation's default, the
## update delay half the burn-in.
writeTreePriorLevel <- function(sampler, given, splitExpr, env) {
  given <- given[!vapply(given, is.null, NA)]
  for (name in names(given)) {
    value <- given[[name]]
    valid <- is.numeric(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      value > 0 &&
      (name == "power" || value < 1)
    if (!valid) {
      stop(
        "'",
        name,
        "' must be a single ",
        if (name == "power") "positive number" else "number in (0, 1)",
        call. = FALSE
      )
    }
    given[[name]] <- as.double(value)
  }
  model <- sampler$model
  ptr <- sampler$getPointer()
  numForests <- bartcoreNumForests(ptr)
  written <- model
  record <- if (!is.null(splitExpr)) {
    splitProbsRecord(splitExpr, ncol(sampler$data@x), env, ifNull = NULL)
  }
  if (!is.null(record)) {
    record@power <- model@tree.prior@power
    record@base <- model@tree.prior@base
    if (!is(record, "dbartsDartPrior")) {
      record <- resolveSplitProbabilities(record, sampler$data)
    } else if (is.na(record@update.delay)) {
      record@update.delay <- as.numeric(sampler$control@n.burn %/% 2L)
    }
    if (numForests > 1L && varySplitProbabilities(record)) {
      refuseLevelSplitProbabilities(sampler, "$setTreePrior")
    }
    ownDart <- is(model@tree.prior, "dbartsDartPrior")
    if (is(record, "dbartsDartPrior") != ownDart) {
      stop(
        if (ownDart) {
          paste0(
            "a vector of split probabilities makes an ordinary sampler of a ",
            "DART one"
          )
        } else {
          "split.probs = dart() makes a DART sampler of an ordinary one"
        },
        ", which changes the model's structure; make a new sampler instead",
        call. = FALSE
      )
    }
    written@tree.prior <- record
  }
  multinomial <- samplerCarriesCounts(sampler)
  candidates <- levelCandidates(sampler)
  for (name in names(given)) {
    value <- given[[name]]
    if (numForests == 1L || multinomial) {
      methods::slot(written@tree.prior, name) <- value
      next
    }
    reach <- levelReach(model, candidates, name)
    if (length(reach) == 0L) {
      refuseLevelReachesNone(model, candidates, name, value, "$setTreePrior")
    }
    if (identical(value, forestDepthValue(model, reach[[1L]], name))) {
      next
    }
    written <- withLevelDepth(written, reach, name, value)
  }
  if (!identical(written, model)) {
    # a refusal on $setModel's route names the door called
    tryCatch(sampler$setModel(written), error = function(e) {
      message <- conditionMessage(e)
      if (!startsWith(message, "$setModel")) {
        stop(e)
      }
      stop("$setTreePrior", substring(message, 10L), call. = FALSE)
    })
  }
  invisible(NULL)
}
