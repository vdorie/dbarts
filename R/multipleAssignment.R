## allows multiple assignment
massign <- structure(NA, class = "lval")

"[<-.lval" <- function(x, ..., value) {
  callingEnv <- parent.frame(1L)

  args <- as.list(match.call())
  args <- args[-c(1L, 2L, length(args))]
  argNames <- names(args)

  valueNames <- names(value)
  if (isS4(value) && !is.list(value)) {
    value <- list(value)
  }
  if (length(value) < length(args)) {
    value <- rep_len(value, length(args))
    if (!is.null(valueNames)) names(value) <- rep_len(valueNames, length(args))
  }

  varNames <- as.character(args)
  if (any(argNames != "")) {
    varNames[argNames != ""] <- argNames[argNames != ""]
  }

  varNamesNoSkip <- varNames[varNames != ""]

  duplicateVarNames <- duplicated(varNamesNoSkip)
  if (any(duplicateVarNames)) {
    warning(
      "names on left-hand-side of assignment appear more than once: ",
      paste0(varNamesNoSkip[duplicateVarNames], collapse = ", "),
      "; result undefined",
      sep = ""
    )
  }

  ## for unnamed rhs, we go entirely by position
  if (is.null(valueNames)) {
    if (any(argNames != "")) {
      warning("right-hand-side of assignment is unnamed; using position only")
    }

    for (i in seq_along(varNames)) {
      var <- args[[i]]

      if (!missing(var)) {
        assign(varNames[i], value[[1L]], envir = callingEnv)
      }

      value <- value[-1L]
    }
    return(x)
  }

  ## go through named args first
  for (i in seq_along(varNames)) {
    if (argNames[i] == "") {
      next
    }

    varName <- varNames[i]
    valueName <- as.character(args[[i]])

    sel <- valueNames == valueName
    numMatches <- sum(sel)

    if (numMatches == 0L) {
      stop(
        "'",
        valueName,
        "' not present in right-hand-side of assignment",
        sep = ""
      )
    }

    if (numMatches > 1L) {
      warning(
        "'",
        valueName,
        "' present multiple times in right-hand-side of assignment; only first will be used",
        sep = ""
      )
      selectionIndex <- which.max(sel)
      sel <- logical(length(value))
      sel[selectionIndex] <- TRUE
    }

    assign(varName, value[sel][[1L]], envir = callingEnv)
    ## check to see if the value is named later, if not pop it off
    if (
      i == length(argNames) ||
        !any(valueName %in% args[seq.int(i + 1L, length(args))])
    ) {
      value <- value[!sel]
      valueNames <- valueNames[!sel]
    }
  }

  ## now for unnamed args
  for (i in seq_along(args)) {
    if (argNames[i] != "") {
      next
    }

    var <- args[[i]]

    if (!missing(var)) {
      assign(as.character(var), value[[1L]], envir = callingEnv)
    }

    ## pop selected values
    value <- value[-1L]
    valueNames <- valueNames[-1L]
  }
  x
}
