## The block of a model of several forests: one line a forest, saying what it
## multiplies, what its column looks like now and the sd it resolved to, in the
## units the reader reports it in. The short block is what a fit's print and a
## sampler's show give: the label, the kind, a numeric column's mean and sd and
## the forest's sd, a numeric forest's per sd of its column as it stands (per
## unit where the column is constant). The full block, which a fit's summary
## and a verbose creation give, has a numeric forest's sd per unit and per sd
## and adds the column's range, where the sd came from and whether the
## coefficient is held, under a header that names the unit.

forestPrintState <- new.env(parent = emptyenv())

## Evaluates expr, a re-creation of a sampler that already exists, with the
## block held off: a copy is not a creation.
withoutForestBlock <- function(expr) {
  held <- forestPrintState[["held"]]
  forestPrintState[["held"]] <- TRUE
  on.exit(forestPrintState[["held"]] <- held)
  expr
}

## A number as the block prints it: format(x, digits = 4), so at least four
## significant digits, and every digit of an integer part.
printedNumber <- function(value) {
  format(unname(value), digits = 4L)
}

## A numeric forest's column as it stands: its mean, standard deviation and
## range over the rows in the likelihood, every row where none is in. The mean
## is zapped against the sd, so a column centered to rounding prints 0. A
## column is constant where its sd is under 1e-11 of its largest value, the
## engine's bound.
printedColumn <- function(column, rowsIn) {
  values <- as.double(column)
  if (any(rowsIn)) {
    values <- values[rowsIn]
  }
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(list(mean = NaN, sd = 0, range = rep(NA_real_, 2L), constant = TRUE))
  }
  spread <- if (length(values) > 1L) stats::sd(values) else 0
  list(
    mean = zapsmall(c(mean(values), spread))[1L],
    sd = spread,
    range = range(values),
    constant = spread <= 1e-11 * max(abs(values))
  )
}

## What a fit keeps of its forests' bases when it drops them (keepFits =
## FALSE): for a numeric basis printedColumn's numbers over the rows in the
## likelihood, for a factor its levels (a factor of no rows), NULL for a forest
## with none. printedForest reads either in place of the basis.
storedBasisColumns <- function(bases, weights = NULL) {
  numRows <- max(0L, vapply(bases, NROW, 0L))
  rows <- rowsInLikelihood(if (length(weights) > 0L) weights, NULL, numRows)
  lapply(bases, function(basis) {
    if (is.null(basis)) {
      NULL
    } else if (is.factor(basis)) {
      basis[0L]
    } else {
      printedColumn(basis, rows)
    }
  })
}

## One forest's line, after its label. `entry` is the reader's entry for it,
## `basis` the forest's stored basis (NULL, a factor or a column), `rowsIn`
## the rows in the likelihood, `full` whether the line is the full block's,
## `label` the forest's, which names the column in a short line, and `family`
## the response's.
printedForest <- function(entry, basis, rowsIn, full, label, family) {
  how <- if (isTRUE(entry$sd.stated)) "stated" else "default"
  size <- printedNumber(entry$sd)
  if (identical(entry$prior.sd.of, "forest, shared")) {
    # a share of the total: its sd where k is held, nothing to print where it
    # is drawn, and never a coefficient of its own
    return(paste0(
      "no basis; ",
      if (is.na(entry$sd)) {
        "a share of the total"
      } else {
        paste0("sd ", size, ", a share of the total")
      }
    ))
  }
  if (is(entry$leaf.prior, "dbartsLeafPrior")) {
    # the fitting function's leaf prior: its size where one number states it,
    # else that prior, and never a coefficient of its own
    if (is.na(entry$sd)) {
      return(paste0("no basis; ", leafPriorText(entry$leaf.prior)))
    }
    return(paste0(
      "no basis; sd ",
      size,
      if (full) {
        paste0(
          " (",
          # a held k of 2 is the default only where the family does not draw k
          if (
            identical(entry$leaf.prior@k, 2) && !drawsLeafKByDefault(family)
          ) {
            "default"
          } else {
            "leaf prior"
          },
          ")"
        )
      }
    ))
  }
  amplitude <- entry[["amplitude"]]
  held <- if (!full || is.null(amplitude)) {
    ""
  } else if (is(amplitude, "dbartsFixedPrior")) {
    "; coefficient held"
  } else if (is(amplitude, "dbartsNormalPrior")) {
    "; coefficient normal"
  } else {
    paste0("; coefficient student(", printedNumber(amplitude@df), ")")
  }
  source <- function(origin = "") {
    if (full) paste0(" (", how, origin, ")") else ""
  }
  switch(
    entry$multiplier,
    # with no sd stated, a drawn forest with no basis has its size from the
    # leaf prior, whichever that is
    none = paste0(
      "no basis; sd ",
      size,
      if (full && how == "default") " (leaf prior)" else source(),
      held
    ),
    factor = paste0(
      "factor basis, ",
      nlevels(basis),
      " levels; sd ",
      size,
      " between the two levels",
      source(),
      held
    ),
    numeric = {
      column <- if (is.list(basis)) basis else printedColumn(basis, rowsIn)
      recorded <- entry[["basis.scale"]]
      range <- if (full) {
        paste0(
          ", range ",
          printedNumber(column$range[1L]),
          " to ",
          printedNumber(column$range[2L])
        )
      } else {
        ""
      }
      perSd <- printedNumber(entry$sd * column$sd)
      origin <- if (
        how == "default" &&
          !is.null(recorded) &&
          abs(recorded - column$sd) > 1e-8 * abs(recorded)
      ) {
        paste0(", taken from a column of sd ", printedNumber(recorded))
      } else {
        ""
      }
      sizes <- if (full) {
        paste0(
          "; sd ",
          size,
          " per unit",
          if (!column$constant) paste0(", ", perSd, " per sd"),
          source(origin)
        )
      } else if (column$constant) {
        paste0("; sd ", size, " per unit")
      } else {
        paste0("; sd ", perSd, " per sd of ", label)
      }
      paste0(
        "numeric basis, mean ",
        printedNumber(column$mean),
        ", sd ",
        printedNumber(column$sd),
        range,
        sizes,
        held
      )
    }
  )
}

## The block, as lines, from the reader's entries of a model of several
## forests and the data they stand on; character(0) for anything else. The
## short block's header is plain and its last line points to `pointer`, where
## the full block is to be had; the full block's header names `unit`, the
## number the sds are stated in, the engine's: sd(y) as the help's `sd` item
## takes it for a continuous response, the link's error sd for a latent one.
forestPriorLines <- function(
  priors,
  labels,
  bases,
  weights = NULL,
  active = NULL,
  family = "gaussian",
  unit,
  full = FALSE,
  pointer = NULL,
  treeCounts = NULL,
  shares = NULL
) {
  if (
    !is.list(priors) ||
      length(priors) < 2L ||
      !is.null(priors[["leaf.prior"]]) ||
      length(bases) != length(priors) ||
      !all(vapply(priors, function(entry) !is.null(entry$multiplier), NA))
  ) {
    return(character(0L))
  }
  # stored numbers (a list) have no rows to count
  numRows <- max(
    0L,
    vapply(bases, function(basis) if (is.list(basis)) 0L else NROW(basis), 0L)
  )
  rows <- rowsInLikelihood(
    if (length(weights) > 0L) weights,
    if (length(active) > 0L) active,
    numRows
  )
  if (length(labels) != length(priors)) {
    labels <- names(priors)
  }
  latent <- family %in% c("probit", "logistic")
  header <- if (full) {
    paste0(
      "forests, sd ",
      if (latent) "on the latent scale" else "in the response's units",
      " (",
      if (latent) "error sd " else "sd(y) ",
      printedNumber(unit),
      "):"
    )
  } else {
    "forests:"
  }
  lines <- vapply(
    seq_along(priors),
    function(index) {
      paste0(
        "  ",
        formatC(index, width = nchar(length(priors))),
        " ",
        labels[index],
        ": ",
        printedForest(
          priors[[index]],
          bases[[index]],
          rows,
          full,
          labels[index],
          family
        )
      )
    },
    ""
  )
  lines <- c(lines, forestTotalLines(priors, treeCounts, shares))
  footer <- if (!full && !is.null(pointer)) {
    paste0("(each forest's prior in full: ", pointer, ")")
  }
  c(header, lines, footer)
}

## The lines after the forests: each forest's tree count and their total,
## where the counts are known, and, where forests share one prior total,
## which do and the group's k, a number where it is held and its prior where
## it is drawn, or the total's named sd.
forestTotalLines <- function(priors, treeCounts, shares = NULL) {
  sharers <- which(vapply(
    priors,
    function(entry) identical(entry$prior.sd.of, "forest, shared"),
    NA
  ))
  counts <- if (length(treeCounts) == length(priors)) {
    paste0(
      "  trees: ",
      paste(treeCounts, collapse = ", "),
      "; ",
      sum(treeCounts),
      " in all"
    )
  }
  # one left of a group, the others having stated their own sd, keeps its
  # part still, printed as 1/n of the total, as a verbose creation prints it
  share <- if (length(sharers) > 0L) {
    spec <- priors[[sharers[1L]]]$leaf.prior
    paste0(
      if (length(sharers) > 1L) "  forests " else "  forest ",
      paste(sharers, collapse = ", "),
      if (length(sharers) > 1L) {
        " share one prior total, "
      } else {
        paste0(
          " keeps its part",
          partText(shares[sharers]),
          " of the prior total, "
        )
      },
      leafPriorText(spec)
    )
  }
  c(counts, share)
}

## A forest's part of the prior total as the block and a verbose creation
## print it, 1/n; "none" where it has no part (NA), "unknown" for anything else.
partFraction <- function(share) {
  if (length(share) != 1L || is.nan(share) || is.na(share)) {
    return("none")
  }
  if (!is.finite(share) || share <= 0) {
    return("unknown")
  }
  if (share == 1) "1" else paste0("1/", round(1 / share))
}

## The part a lone sharer keeps, set off by commas for the line it is in; empty
## where the shares are not known.
partText <- function(share) {
  fraction <- partFraction(share)
  if (fraction %in% c("none", "unknown")) "" else paste0(", ", fraction, ",")
}

## A leaf prior as the block prints it: its k or its named sd.
leafPriorText <- function(spec) {
  if (!is.null(spec@prior.sd)) {
    paste0("sd ", deparsedPrior(spec@prior.sd))
  } else {
    paste0("k ", deparsedPrior(spec@k))
  }
}

## A prior object as the block prints it, chi(1.5, 2) for a chi prior on k
## and invchi(df, scale) for one on the sd; a number as printedNumber gives it.
deparsedPrior <- function(prior) {
  if (is(prior, "dbartsChiHyperprior")) {
    paste0(
      "chi(",
      printedNumber(prior@degreesOfFreedom),
      ", ",
      printedNumber(prior@scale),
      ")"
    )
  } else if (is(prior, "dbartsSdHyperprior")) {
    paste0(
      "invchi(",
      printedNumber(prior@df),
      ", ",
      printedNumber(prior@scale),
      ")"
    )
  } else {
    printedNumber(prior)
  }
}

## The block of a fit, from what the fit stores: the short one of its synopsis,
## the full one of its summary.
fitForestLines <- function(fit, full = FALSE) {
  forestPriorLines(
    fit[["leaf.prior"]],
    attr(fit, "forest.labels", exact = TRUE),
    if (is.null(fit[["bases"]])) fit[["basis.columns"]] else fit[["bases"]],
    fit[["weights"]],
    fit[["active"]],
    fit[["family"]],
    attr(fit, "forest.unit", exact = TRUE),
    full,
    "summary()",
    attr(fit, "forest.trees", exact = TRUE)
  )
}

## What one forest does per unit of its basis, each draw's effect at each row:
## the forest's own fits times the amplitude on the basis column, which is
## what predict() at the column one unit higher gives beyond the column as it
## stands, whatever the other columns hold. A factor's effect is a level's
## against the first, the difference of the two levels' amplitudes. Returns a
## list of draws-by-rows matrices, one for each column the forest multiplies
## (a factor's second and later levels), named by what the effect is of.
forestEffectDraws <- function(fits, glue, multiplier, basis) {
  width <- ncol(glue)
  if (multiplier == "factor") {
    levels <- levels(basis)
    if (length(levels) != width) {
      return(list())
    }
    steps <- sweep(glue[, -1L, drop = FALSE], 1L, glue[, 1L])
    names <- paste0("level ", levels[-1L], " against ", levels[1L])
  } else {
    steps <- glue
    names <- if (width == 1L) {
      "per unit"
    } else if (!is.null(colnames(basis))) {
      colnames(basis)
    } else {
      paste0("column ", seq_len(width))
    }
  }
  effects <- lapply(seq_len(ncol(steps)), function(column) {
    steps[, column] * fits
  })
  names(effects) <- names
  effects
}

## One effect's summary: the average over the rows in the likelihood and its
## spread across them (sd), each taken per draw, so each is a posterior of the
## quantity itself, then the median and the 5% and 95% quantiles summary()
## tabulates elsewhere.
effectSummary <- function(draws, rowsIn) {
  draws <- draws[, rowsIn, drop = FALSE]
  average <- rowMeans(draws)
  spread <- if (ncol(draws) > 1L) {
    sqrt(rowSums((draws - average)^2) / (ncol(draws) - 1L))
  } else {
    rep(NA_real_, nrow(draws))
  }
  quantiles <- function(x) {
    if (anyNA(x)) {
      rep(NA_real_, 3L)
    } else {
      stats::quantile(x, c(0.5, 0.05, 0.95), names = FALSE)
    }
  }
  data.frame(
    average = quantiles(average)[1L],
    average.q5 = quantiles(average)[2L],
    average.q95 = quantiles(average)[3L],
    spread = quantiles(spread)[1L],
    spread.q5 = quantiles(spread)[2L],
    spread.q95 = quantiles(spread)[3L]
  )
}

## The scale the forests add on, as a phrase for the line that reports an
## effect: the response's units unless the family's forests live on a link or
## a log scale.
effectScaleText <- function(family) {
  switch(
    as.character(family)[1L],
    probit = "on the probit scale",
    ordinal = "on the probit scale",
    logistic = "on the logit scale",
    aft = "on the log time scale",
    nbinom = "on the log mean scale",
    hazard.probit = "on the probit scale",
    "in the response's units"
  )
}

## A fit's forests' effects per unit of their bases, from its per-forest fits
## and amplitudes: a data frame of one row for each column of each forest with
## a basis (forest, label, effect, and effectSummary's columns), with the rows
## it averages over and the scale as attributes; NULL for a fit with no forest
## that has a basis, and, for one that has but kept no per-forest fits, the
## string "dropped".
fitForestEffects <- function(fit) {
  priors <- fit[["leaf.prior"]]
  multipliers <- vapply(
    priors,
    function(entry) as.character(entry$multiplier)[1L],
    ""
  )
  has <- which(multipliers %in% c("factor", "numeric"))
  if (length(has) == 0L) {
    return(NULL)
  }
  if (is.null(fit[["forestFits"]]) || is.null(fit[["glue"]])) {
    return("dropped")
  }
  chains <- fitNChains(fit)
  fits <- reshapeChainedChannel(fit$forestFits, chains, TRUE, 2L)
  glue <- reshapeChainedChannel(fit$glue, chains, TRUE, 1L)
  glueForest <- attr(fit$glue, "forest")
  forestNames <- dimnames(fits)[[3L]]
  labels <- attr(fit, "forest.labels", exact = TRUE)
  if (length(labels) != length(priors)) {
    labels <- names(priors)
  }
  rowsIn <- rowsInLikelihood(
    if (length(fit[["weights"]]) > 0L) fit[["weights"]],
    if (length(fit[["active"]]) > 0L) fit[["active"]],
    dim(fits)[2L]
  )
  rows <- lapply(has, function(k) {
    effects <- forestEffectDraws(
      fits[,, k],
      glue[, glueForest == forestNames[k], drop = FALSE],
      multipliers[k],
      fit$bases[[k]]
    )
    do.call(
      rbind,
      Map(
        function(draws, what) {
          cbind(
            data.frame(
              forest = k,
              label = labels[k],
              effect = what,
              stringsAsFactors = FALSE
            ),
            effectSummary(draws, rowsIn)
          )
        },
        effects,
        names(effects)
      )
    )
  })
  result <- do.call(rbind, rows)
  if (is.null(result)) {
    return(NULL)
  }
  rownames(result) <- NULL
  attr(result, "rows") <- sum(rowsIn)
  attr(result, "scale") <- effectScaleText(fit[["family"]])
  result
}

## The lines under the block that report fitForestEffects: what each number is
## once, then a line for each effect.
forestEffectLines <- function(effects) {
  if (is.null(effects)) {
    return(character(0L))
  }
  if (identical(effects, "dropped")) {
    return("  effects of the bases not shown: the fit kept no per-forest fits")
  }
  shown <- function(median, lower, upper) {
    paste0(
      printedNumber(median),
      " (",
      printedNumber(lower),
      " to ",
      printedNumber(upper),
      ")"
    )
  }
  lines <- vapply(
    seq_len(nrow(effects)),
    function(row) {
      with(
        effects[row, ],
        paste0(
          "  ",
          forest,
          " ",
          label,
          ", ",
          effect,
          ": average ",
          shown(average, average.q5, average.q95),
          ", spread across rows ",
          shown(spread, spread.q5, spread.q95)
        )
      )
    },
    ""
  )
  c(
    paste0(
      "effect of one unit of each basis (a factor's: one level against the ",
      "first), ",
      attr(effects, "scale"),
      ", over ",
      attr(effects, "rows"),
      " training rows (posterior median, 5% to 95%):"
    ),
    lines
  )
}

## The block of a sampler, from its data as it stands: the short one of its
## show, the full one of a verbose creation; character(0) where it has one
## forest, or its pointer has not survived a reload, which show does not
## re-create.
samplerForestLines <- function(sampler, full = FALSE) {
  records <- sampler$model@forests
  if (
    length(records) < 2L ||
      !.Call(C_dbarts_bartcore_isValidPointer, sampler$pointer) ||
      bartcoreNumForests(sampler$pointer) < 2L
  ) {
    return(character(0L))
  }
  forestPriorLines(
    sampler$getLeafPrior(),
    vapply(records, methods::slot, "", "label"),
    storedForestBases(sampler$data),
    sampler$data@weights,
    sampler$activeRows,
    sampler$model@family,
    dataForestCalibration(sampler$data)[["unit"]],
    full,
    "$getLeafPrior()",
    vapply(records, function(record) as.integer(record@n.trees), 0L),
    vapply(records, methods::slot, 0, "share")
  )
}
