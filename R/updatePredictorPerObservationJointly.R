# Joint per-observation update of a single shared column across several samplers.
# This is deliberately a package-level function rather than a dbartsSampler method:
# it acts on a *set* of peer samplers, so no single sampler is a privileged receiver.
updatePredictorPerObservationJointly <- function(samplers, x, column, updateState = NA) {
  if (is(samplers, "dbartsSampler")) samplers <- list(samplers)
  if (!is.list(samplers) || length(samplers) == 0L)
    stop("'samplers' must be a dbartsSampler or non-empty list of dbartsSampler objects")
  for (f in samplers)
    if (!is(f, "dbartsSampler")) stop("'samplers' must be a dbartsSampler or list of dbartsSampler objects")

  if (missing(column)) stop("joint updates require a single 'column' to be specified")
  if (length(column) != 1L) stop("joint updates can only be applied to a single column")

  # resolve the (single, shared) column to a per-sampler index by NAME, since the same variable
  # can sit at different positions in each sampler's design matrix
  columnName <- if (is.character(column)) {
    column
  } else {
    nm <- colnames(samplers[[1L]]$data@x)
    if (is.null(nm))
      stop("column names not specified at initialization, so a numeric column cannot be matched across samplers")
    column <- as.integer(column)
    if (is.na(column) || column < 1L || column > length(nm)) stop("column is out of range")
    nm[column]
  }

  cols <- integer(length(samplers))
  for (i in seq_along(samplers)) {
    nm <- colnames(samplers[[i]]$data@x)
    if (is.null(nm)) stop("column names not specified at initialization, so cannot be replaced by name")
    idx <- match(columnName, nm)
    if (is.na(idx)) stop(gettextf("column name '%s' not found in names of X for sampler %d", columnName, i))
    cols[i] <- idx
  }

  x <- as.double(x)

  ptrs <- lapply(samplers, function(f) f$getPointer())
  installed <- .Call(C_dbarts_updatePredictorPerObservationJointly, ptrs, x, as.integer(cols))

  for (f in samplers)
    if ((is.na(updateState) && f$control@updateState == TRUE) || identical(updateState, TRUE))
      f$storeState()

  installed
}
