## The block of a model of several forests: one line a forest, saying what it
## multiplies, what its column looks like now and the sd it resolved to, in the
## units the reader reports it in. The short block is what a fit's print and a
## sampler's show give: the label, the kind, a numeric column's mean and sd and
## the sd. The full block, which a fit's summary and a verbose creation give,
## adds the column's range, the sd per sd of the column, where the sd came
## from and whether the coefficient is held, under a header that names the unit.

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

## One forest's line, after its label. `entry` is the reader's entry for it,
## `basis` the forest's stored basis (NULL, a factor or a column), `rowsIn`
## the rows in the likelihood and `full` whether the line is the full block's.
printedForest <- function(entry, basis, rowsIn, full) {
  how <- if (isTRUE(entry$sd.stated)) "stated" else "default"
  size <- printedNumber(entry$sd)
  held <- if (full && !is.null(entry[["amplitude"]])) {
    "; coefficient held"
  } else {
    ""
  }
  source <- function(origin = "") {
    if (full) paste0(" (", how, origin, ")") else ""
  }
  switch(
    entry$multiplier,
    none = paste0("no basis; sd ", size, source(), held),
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
      column <- printedColumn(basis, rowsIn)
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
      perSd <- if (!full || column$constant) {
        ""
      } else {
        paste0(", ", printedNumber(entry$sd * column$sd), " per sd")
      }
      origin <- if (
        how == "default" &&
          !is.null(recorded) &&
          abs(recorded - column$sd) > 1e-8 * abs(recorded)
      ) {
        paste0(", taken from a column of sd ", printedNumber(recorded))
      } else {
        ""
      }
      paste0(
        "numeric basis, mean ",
        printedNumber(column$mean),
        ", sd ",
        printedNumber(column$sd),
        range,
        "; sd ",
        size,
        " per unit",
        perSd,
        source(origin),
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
  pointer = NULL
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
  numRows <- max(0L, vapply(bases, NROW, 0L))
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
        printedForest(priors[[index]], bases[[index]], rows, full)
      )
    },
    ""
  )
  footer <- if (!full && !is.null(pointer)) {
    paste0("(each forest's prior in full: ", pointer, ")")
  }
  c(header, lines, footer)
}

## The block of a fit, from what the fit stores: the short one of its synopsis,
## the full one of its summary.
fitForestLines <- function(fit, full = FALSE) {
  forestPriorLines(
    fit[["leaf.prior"]],
    attr(fit, "forest.labels", exact = TRUE),
    fit[["bases"]],
    fit[["weights"]],
    fit[["active"]],
    fit[["family"]],
    attr(fit, "forest.unit", exact = TRUE),
    full,
    "summary()"
  )
}

## The block of a sampler, from its data as it stands: the short one of its
## show, the full one of a verbose creation; character(0) where it has one
## forest, or its pointer has not survived a reload, which show does not
## re-create.
samplerForestLines <- function(sampler, full = FALSE) {
  info <- attr(sampler$control, "bartcore.forests", exact = TRUE)
  if (
    is.null(info) ||
      !.Call(C_dbarts_bartcore_isValidPointer, sampler$pointer) ||
      bartcoreNumForests(sampler$pointer) < 2L
  ) {
    return(character(0L))
  }
  forestPriorLines(
    sampler$getLeafPrior(),
    info$labels,
    storedForestBases(sampler$data),
    sampler$data@weights,
    sampler$activeRows,
    sampler$model@family,
    info$unit,
    full,
    "$getLeafPrior()"
  )
}
