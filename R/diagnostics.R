# convergence diagnostics: posterior-package draws accessors and a
# summary() method for the scalar parameters of bart/bart2/rbart fits.
# posterior is Suggests-only; as_draws_array/as_draws_df are registered
# for it conditionally in zzz.R. summary() itself is a base generic and
# degrades to a plain quantile table when posterior is unavailable.

# indirection point so tests can simulate posterior's absence without
# uninstalling it
posteriorAvailable <- function() requireNamespace("posterior", quietly = TRUE)

# n.chains survives on the object whether or not the sampler was kept (see
# packageBartResults/packageRbartResults); fit is a single dbartsSampler for
# bart/bart2 and a list of them (length n.chains, or 1 for the in-core
# multi-chain rbart path, which always keeps n.chains on the object too) for
# rbart
fitNChains <- function(object) {
  if (!is.null(object[["n.chains"]])) {
    return(object[["n.chains"]])
  }
  fit <- object[["fit"]]
  if (inherits(fit, "dbartsSampler")) fit$control@n.chains else length(fit)
}

# maps one field's native bart-convention samples - (n.chains, n.samples[,
# n.vars]), collapsed to drop the chain dimension when n.chains == 1, or
# flattened to a vector/matrix when combineChains was requested at fit time
# - to posterior's (iteration, chain, variable) array. uncombineChains
# already knows how to invert the flattening; only the trailing transpose
# is new here.
toDrawsArray <- function(x, n.chains) {
  d <- dim(x)
  if (n.chains <= 1L) {
    if (is.null(d)) {
      arr <- array(x, c(length(x), 1L, 1L))
      varNames <- NULL
    } else {
      arr <- array(x, c(d[1L], 1L, d[2L]))
      varNames <- dimnames(x)[[2L]]
    }
  } else if (is.null(d)) {
    mat <- uncombineChains(x, n.chains) # n.chains x n.samples
    arr <- array(t(mat), c(ncol(mat), n.chains, 1L))
    varNames <- NULL
  } else if (length(d) == 2L) {
    arr <- array(t(x), c(d[2L], d[1L], 1L))
    varNames <- NULL
  } else {
    arr <- aperm(x, c(2L, 1L, 3L))
    varNames <- dimnames(x)[[3L]]
  }
  dimnames(arr) <- list(
    NULL,
    NULL,
    if (is.null(varNames)) as.character(seq_len(dim(arr)[3L])) else varNames
  )
  arr
}

# fields with no per-variable axis; every other requested field (varcount,
# varprobs, yhat.train, yhat.test, ranef, ...) contributes one draws
# variable per column/observation, named "field[inner]"
scalarFields <- c("sigma", "k", "tau")

# gathers one or more chain-dimensioned fields off a bart/bart2/rbart fit
# into a single (iteration, chain, variable) base array
bartDrawsArray <- function(object, vars) {
  n.chains <- fitNChains(object)
  present <- vars[!vapply(vars, function(v) is.null(object[[v]]), logical(1L))]
  if (length(present) == 0L) {
    stop(
      "none of 'vars' (",
      paste0(vars, collapse = ", "),
      ") are present on this fit",
      call. = FALSE
    )
  }
  pieces <- lapply(present, function(v) {
    piece <- toDrawsArray(object[[v]], n.chains)
    dimnames(piece)[[3L]] <- if (v %in% scalarFields) {
      v
    } else {
      paste0(v, "[", dimnames(piece)[[3L]], "]")
    }
    piece
  })
  varNames <- unlist(lapply(pieces, function(p) dimnames(p)[[3L]]))
  array(
    unlist(pieces),
    dim = c(dim(pieces[[1L]])[1L], n.chains, length(varNames)),
    dimnames = list(NULL, NULL, varNames)
  )
}

as_draws_array.bart <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_array(bartDrawsArray(x, vars))
}
as_draws_array.rbart <- as_draws_array.bart

as_draws_df.bart <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_df(bartDrawsArray(x, vars))
}
as_draws_df.rbart <- as_draws_df.bart

# plain per-variable mean/sd/quantiles, pooling samples and chains; the
# degrade path when posterior is not installed
quantileSummary <- function(arr) {
  pooled <- matrix(arr, prod(dim(arr)[1L:2L]), dim(arr)[3L])
  data.frame(
    variable = dimnames(arr)[[3L]],
    mean = colMeans(pooled),
    sd = apply(pooled, 2L, sd),
    q2.5 = apply(pooled, 2L, quantile, probs = 0.025, names = FALSE),
    median = apply(pooled, 2L, quantile, probs = 0.5, names = FALSE),
    q97.5 = apply(pooled, 2L, quantile, probs = 0.975, names = FALSE),
    row.names = NULL
  )
}

# rhat > 1.01 is noted in the printed summary, not enforced: dbarts does not
# refuse to summarize a non-converged fit
summary.bart <- function(object, vars = c("sigma", "k", "tau"), ...) {
  present <- vars[!vapply(vars, function(v) is.null(object[[v]]), logical(1L))]
  havePosterior <- length(present) > 0L && posteriorAvailable()
  stats <- if (length(present) == 0L) {
    NULL
  } else if (havePosterior) {
    posterior::summarise_draws(posterior::as_draws_array(bartDrawsArray(
      object,
      present
    )))
  } else {
    quantileSummary(bartDrawsArray(object, present))
  }
  structure(
    list(call = object[["call"]], stats = stats, posterior = havePosterior),
    class = "summary.bart"
  )
}
summary.rbart <- summary.bart

print.summary.bart <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  if (is.null(x$stats)) {
    cat("No scalar parameters (sigma, k, tau) to summarize.\n")
    return(invisible(x))
  }
  print(x$stats, ...)
  if (!x$posterior) {
    cat("\nNote: install 'posterior' for R-hat/ESS convergence diagnostics.\n")
  } else if (any(x$stats$rhat > 1.01, na.rm = TRUE)) {
    cat(
      "\nNote: some R-hat values exceed 1.01; chains may not have converged.\n"
    )
  }
  invisible(x)
}
