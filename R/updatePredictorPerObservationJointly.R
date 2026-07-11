# Joint per-observation update of a single shared column across several
# samplers. Following dbartsSampler convention, NA updateState pulls from the
# sampler's control; if it's a logical scalar, that value is used.
updatePredictorPerObservationJointly <- function(
  samplers,
  x,
  column,
  updateState = NA
) {
  if (inherits(samplers, "dbartsSampler")) {
    samplers <- list(samplers)
  }
  if (!is.list(samplers) || length(samplers) == 0L) {
    stop(
      "'samplers' must be a dbartsSampler or non-empty list of dbartsSampler ",
      "objects"
    )
  }
  for (sampler in samplers) {
    if (!inherits(sampler, "dbartsSampler")) {
      stop(
        "'samplers' must be a dbartsSampler or list of dbartsSampler objects"
      )
    }
  }

  if (missing(column)) {
    stop("joint updates require a single 'column' to be specified")
  }
  if (length(column) != 1L) {
    stop("joint updates can only be applied to a single column")
  }

  # resolve the (single, shared) column to a per-sampler index by NAME, since
  # the same variable can sit at different positions in each sampler's design
  # matrix
  if (is.character(column)) {
    columnName <- column
  } else {
    columNames <- colnames(samplers[[1L]]$data@x)
    if (is.null(columNames)) {
      stop(
        "column names not specified at initialization, so a numeric column ",
        "cannot be matched across samplers"
      )
    }
    column <- as.integer(column)
    if (is.na(column) || column < 1L || column > length(columNames)) {
      stop("column is out of range")
    }
    columnName <- columNames[column]
  }

  columnIndices <- integer(length(samplers))
  for (i in seq_along(samplers)) {
    columNames <- colnames(samplers[[i]]$data@x)
    if (is.null(columNames)) {
      stop(
        "column names not specified at initialization, so cannot be replaced ",
        "by name"
      )
    }
    columnIndex <- match(columnName, columNames)
    if (is.na(columnIndex)) {
      stop(sprintf(
        "column name '%s' not found in names of X for sampler %d",
        columnName,
        i
      ))
    }
    columnIndices[i] <- columnIndex
  }

  x <- as.double(x)

  ptrs <- lapply(samplers, function(sampler) sampler$getPointer())
  installed <- .Call(
    C_dbarts_bartcore_updatePredictorPerObservationJointly,
    ptrs,
    x,
    as.integer(columnIndices)
  )

  # the engine keeps no predictor matrix, so maintain each sampler's data@x
  # R-side for the observations the shared scan installed (design plan-1 interim)
  for (i in seq_along(samplers)) {
    samplers[[i]]$data@x[installed, columnIndices[i]] <- x[installed]
  }

  for (sampler in samplers) {
    if (
      (is.na(updateState) && sampler$control@updateState == TRUE) ||
        identical(updateState, TRUE)
    ) {
      sampler$storeState()
    }
  }

  installed
}
