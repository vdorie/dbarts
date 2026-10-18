## The block a fit's print, a sampler's show and a verbose creation give a
## model of several forests: one line a forest, saying what it multiplies,
## what its columns look like now and the sd it resolved to, in the units the
## reader reports it in.

forestPrintState <- new.env(parent = emptyenv())

## Evaluates expr, a re-creation of a sampler that already exists, with the
## block held off: a copy is not a creation.
withoutForestBlock <- function(expr) {
  held <- forestPrintState[["held"]]
  forestPrintState[["held"]] <- TRUE
  on.exit(forestPrintState[["held"]] <- held)
  expr
}

## A number as the block prints it: four significant digits.
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
## `basis` the forest's stored basis (NULL, a factor or a column) and `rowsIn`
## the rows in the likelihood.
printedForest <- function(entry, basis, rowsIn) {
  how <- if (isTRUE(entry$sd.stated)) "stated" else "default"
  size <- printedNumber(entry$sd)
  held <- if (!is.null(entry[["amplitude"]])) "; coefficient held" else ""
  switch(
    entry$multiplier,
    none = paste0("no basis; sd ", size, " (", how, ")", held),
    factor = paste0(
      "factor basis, ",
      nlevels(basis),
      " levels; sd ",
      size,
      " between the two levels (",
      how,
      ")",
      held
    ),
    numeric = {
      column <- printedColumn(basis, rowsIn)
      recorded <- entry[["basis.scale"]]
      perSd <- if (column$constant) {
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
        ", range ",
        printedNumber(column$range[1L]),
        " to ",
        printedNumber(column$range[2L]),
        "; sd ",
        size,
        " per unit",
        perSd,
        " (",
        how,
        origin,
        ")",
        held
      )
    }
  )
}

## The block, as lines, from the reader's entries of a model of several
## forests and the data they stand on; character(0) for anything else. `unit`
## is the number the sds are stated in, the response's sd or the link's error
## sd, and is left out of the header where it is not known.
forestPriorLines <- function(
  priors,
  labels,
  bases,
  weights = NULL,
  active = NULL,
  family = "gaussian",
  unit = NULL
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
  header <- paste0(
    "forests, sd ",
    if (latent) "on the latent scale" else "in the response's units",
    if (!is.null(unit)) {
      paste0(
        " (",
        if (latent) "error sd " else "its sd ",
        printedNumber(unit),
        ")"
      )
    },
    ":"
  )
  lines <- vapply(
    seq_along(priors),
    function(index) {
      paste0(
        "  ",
        formatC(index, width = nchar(length(priors))),
        " ",
        labels[index],
        ": ",
        printedForest(priors[[index]], bases[[index]], rows)
      )
    },
    ""
  )
  c(header, lines)
}

## The block of a fit's synopsis, from what the fit stores.
fitForestLines <- function(fit) {
  forestPriorLines(
    fit[["leaf.prior"]],
    attr(fit, "forest.labels", exact = TRUE),
    fit[["bases"]],
    fit[["weights"]],
    fit[["active"]],
    fit[["family"]],
    attr(fit, "forest.unit", exact = TRUE)
  )
}

## The block of a sampler, from its data as it stands; character(0) where it
## has one forest, or its pointer has not survived a reload, which show does
## not re-create.
samplerForestLines <- function(sampler) {
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
    info$unit
  )
}
