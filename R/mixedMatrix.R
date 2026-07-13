## Internal containers for predictor sets that stay columnar instead of
## materializing one cbound double matrix. Two flavors share the
## dbartsMixedMatrix class, told apart by the type of their "dense" field:
##
## - the mixed flavor (dense = matrix or NULL): built when a data frame
##   holds Matrix::sparseVector or dgCMatrix columns alongside ordinary
##   ones - the dense columns as one matrix, the sparse columns cbound into
##   one dgCMatrix, and a map from overall column to its source, positive k
##   for dense column k, negative -k for sparse column k. The dense block
##   stays resident: the engine retains per-column slices of it.
## - the dense flavor (dense = list): built for every other data frame -
##   one element per predictor column, factors keeping their integer codes
##   and everything else held as doubles; sparse is NULL and the map all
##   positive. The contiguous coded block exists only transiently,
##   assembled by the bridge inside create/setData calls and by as.matrix
##   for R-side consumers.
##
## The model-matrix builders attach "varTypes"/"factor.levels"/
## "term.labels"/"drop" exactly as they do on matrices. A plain R list
## underneath, so serialization and sampler re-creation work unchanged.

isSparseDataFrameColumn <- function(column) {
  isS4(column) &&
    (methods::is(column, "sparseVector") ||
      methods::is(column, "sparseMatrix") ||
      methods::is(column, "sparseFactor"))
}

## The 0-based row indices and values of each predictor column a sparse
## data-frame column contributes - a sparseVector contributes one, a
## dgCMatrix its columns. Missing entries are stored NaNs, the Matrix
## convention.
sparseColumnSlices <- function(column, name, numObservations) {
  # S-CAT: the wrapper is recognized, but the CSC-categorical engine path
  # is plan 5's; refuse before anything reaches the bridge
  if (methods::is(column, "sparseFactor")) {
    stop(
      "sparse categorical predictors are not yet supported; column '",
      name,
      "' is a sparseFactor"
    )
  }
  if (methods::is(column, "sparseVector")) {
    if (length(column) != numObservations) {
      stop(
        "sparse column '",
        name,
        "' must have length equal to the ",
        "number of observations"
      )
    }
    values <- if (methods::.hasSlot(column, "x")) {
      as.double(column@x)
    } else {
      rep_len(1.0, length(column@i))
    }
    return(list(
      i = list(as.integer(column@i) - 1L),
      x = list(values),
      names = name
    ))
  }
  if (!inherits(column, "dgCMatrix")) {
    stop("sparse matrix column '", name, "' must be a Matrix::dgCMatrix")
  }
  if (nrow(column) != numObservations) {
    stop(
      "sparse column '",
      name,
      "' must have rows equal to the number ",
      "of observations"
    )
  }
  columnNnz <- diff(column@p)
  i <- vector("list", ncol(column))
  values <- vector("list", ncol(column))
  for (k in seq_len(ncol(column))) {
    range <- seq.int(column@p[k] + 1L, length.out = columnNnz[k])
    i[[k]] <- column@i[range]
    values[[k]] <- column@x[range]
  }
  columnNames <- paste(
    name,
    if (!is.null(colnames(column))) colnames(column) else seq_len(ncol(column)),
    sep = "."
  )
  list(i = i, x = values, names = columnNames)
}

## Assemble expanded per-input-column blocks - dense numeric matrices or
## sparseColumnSlices results, flagged by columnIsSparse - into the mixed
## container, preserving the input order. blockNames supplies the overall
## column names per block.
assembleMixedMatrix <- function(
  columns,
  columnIsSparse,
  blockNames,
  numObservations
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("sparse predictor columns require the Matrix package")
  }

  widths <- integer(length(columns))
  for (j in seq_along(columns)) {
    widths[j] <- if (columnIsSparse[j]) {
      length(columns[[j]]$i)
    } else {
      ncol(columns[[j]])
    }
  }

  map <- integer(sum(widths))
  numDense <- 0L
  numSparse <- 0L
  offset <- 0L
  for (j in seq_along(columns)) {
    if (widths[j] == 0L) {
      next
    }
    map[offset + seq_len(widths[j])] <- if (columnIsSparse[j]) {
      -(numSparse + seq_len(widths[j]))
    } else {
      numDense + seq_len(widths[j])
    }
    if (columnIsSparse[j]) {
      numSparse <- numSparse + widths[j]
    } else {
      numDense <- numDense + widths[j]
    }
    offset <- offset + widths[j]
  }

  dense <- if (numDense > 0L) {
    do.call(cbind, columns[!columnIsSparse & widths > 0L])
  } else {
    NULL
  }

  rowIndices <- unlist(lapply(columns[columnIsSparse], `[[`, "i"))
  if (is.null(rowIndices)) {
    rowIndices <- integer(0)
  }
  values <- unlist(lapply(columns[columnIsSparse], `[[`, "x"))
  if (is.null(values)) {
    values <- double(0)
  }
  entryCounts <- unlist(lapply(columns[columnIsSparse], function(slices) {
    lengths(slices$i)
  }))
  pointers <- c(0L, cumsum(as.integer(entryCounts)))
  sparse <- methods::new(
    "dgCMatrix",
    i = as.integer(rowIndices),
    p = as.integer(pointers),
    x = as.double(values),
    Dim = c(as.integer(numObservations), numSparse)
  )

  result <- list(
    dense = dense,
    sparse = sparse,
    map = map,
    numObservations = as.integer(numObservations),
    columnNames = unlist(blockNames)
  )
  class(result) <- "dbartsMixedMatrix"
  result
}

## Assemble per-input-column blocks - factors, double vectors, or double
## matrices (spliced per column) - into the dense columnar flavor.
assembleDenseColumnMatrix <- function(columns, blockNames, numObservations) {
  widths <- vapply(
    columns,
    function(block) if (is.matrix(block)) ncol(block) else 1L,
    0L
  )
  dense <- vector("list", sum(widths))
  offset <- 0L
  for (j in seq_along(columns)) {
    block <- columns[[j]]
    if (is.matrix(block)) {
      for (k in seq_len(ncol(block))) {
        dense[[offset + k]] <- block[, k]
      }
    } else {
      dense[[offset + 1L]] <- block
    }
    offset <- offset + widths[j]
  }
  result <- list(
    dense = dense,
    sparse = NULL,
    map = seq_along(dense),
    numObservations = as.integer(numObservations),
    columnNames = unlist(blockNames)
  )
  class(result) <- "dbartsMixedMatrix"
  result
}

dim.dbartsMixedMatrix <- function(x) {
  n <- if (!is.null(x$numObservations)) {
    x$numObservations
  } else {
    nrow(x$sparse)
  }
  c(n, length(x$map))
}

dimnames.dbartsMixedMatrix <- function(x) {
  if (is.null(x$columnNames)) NULL else list(NULL, x$columnNames)
}

## Row subsetting keeps the container (and its builder attributes); any
## column selection densifies, with matrix drop semantics.
`[.dbartsMixedMatrix` <- function(x, i, j, drop = TRUE) {
  if (!missing(j)) {
    return(as.matrix(x)[if (missing(i)) TRUE else i, j, drop = drop])
  }
  if (missing(i)) {
    return(x)
  }
  result <- x
  if (is.list(result$dense)) {
    result$dense <- lapply(result$dense, function(column) column[i])
    result$numObservations <- length(result$dense[[1L]])
  } else {
    if (!is.null(result$dense)) {
      result$dense <- result$dense[i, , drop = FALSE]
    }
    result$sparse <- result$sparse[i, , drop = FALSE]
    result$numObservations <- nrow(result$sparse)
  }
  result
}

as.matrix.dbartsMixedMatrix <- function(x, ...) {
  result <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  if (is.list(x$dense)) {
    # the dense flavor: every column vector-backed, factors materializing
    # as their zero-based codes
    for (j in seq_along(x$map)) {
      column <- x$dense[[x$map[j]]]
      result[, j] <- if (is.factor(column)) {
        as.double(as.integer(column) - 1L)
      } else {
        column
      }
    }
    return(result)
  }
  denseColumns <- x$map > 0L
  if (any(denseColumns)) {
    result[, denseColumns] <- x$dense[, x$map[denseColumns], drop = FALSE]
  }
  if (any(!denseColumns)) {
    result[, !denseColumns] <-
      as.matrix(x$sparse[, -x$map[!denseColumns], drop = FALSE])
  }
  result
}

## Copy-modify a block of predictor columns in a sampler's stored source -
## the plan-1 interim write-back that keeps data@x current under mutation,
## extended to the dense columnar flavor: only the addressed columns'
## vectors are replaced (a factor column decays to its code vector; the
## engine validated the new values as existing codes), the rest stay
## shared. rows = NULL replaces the columns whole. Mutation is refused for
## sparse sources before this runs, so a container here is the dense
## flavor.
assignIntoPredictorSource <- function(x, rows, columns, values) {
  if (is.matrix(x)) {
    if (is.null(rows)) {
      x[, columns] <- values
    } else {
      x[rows, columns] <- values
    }
    return(x)
  }
  values <- matrix(as.double(values), ncol = length(columns))
  for (k in seq_along(columns)) {
    sourceIndex <- x$map[columns[k]]
    column <- x$dense[[sourceIndex]]
    if (is.factor(column)) {
      column <- as.double(as.integer(column) - 1L)
    }
    if (is.null(rows)) {
      column[] <- values[, k]
    } else {
      column[rows] <- values[, k]
    }
    x$dense[[sourceIndex]] <- column
  }
  x
}

## Install a block of predictor columns into a sampler's stored source BY
## REFERENCE (design/data-ownership.md plan 3) - keeps data@x current under
## mutation without adding an R-side copy on top of R's own copy-on-write.
## rows = NULL replaces columns whole: slot j is repointed straight at the
## supplied vector (or, for a multi-column call, the caller's per-column
## slice of it), so an unmutated column stays shared with the prior
## container and only the O(p) list spine is duplicated; a factor source
## column decays to a plain code vector automatically, since the old column
## is discarded rather than mutated in place. rows given (a partial merge)
## starts from the old column and overwrites only the addressed rows - the
## merge is inherent for a freshly supplied vector (O(n) for one column, no
## worse than before); the O(spine) win there needs the caller's own
## already-merged vector reinstalled in place. A matrix source has no
## columnar spine to share, so it keeps copy-modify. Mutation is refused for
## sparse sources before this runs, so a container here is the dense
## flavor.
installPredictorColumns <- function(x, rows, columns, values) {
  if (is.matrix(x)) {
    if (is.null(rows)) {
      x[, columns] <- values
    } else {
      x[rows, columns] <- values
    }
    return(x)
  }
  if (is.null(rows)) {
    values <- as.double(values)
    n <- x$numObservations
    for (k in seq_along(columns)) {
      x$dense[[x$map[columns[k]]]] <- if (length(columns) == 1L) {
        values
      } else {
        values[seq.int((k - 1L) * n + 1L, length.out = n)]
      }
    }
    return(x)
  }
  values <- matrix(as.double(values), ncol = length(columns))
  for (k in seq_along(columns)) {
    sourceIndex <- x$map[columns[k]]
    column <- x$dense[[sourceIndex]]
    if (is.factor(column)) {
      column <- as.double(as.integer(column) - 1L)
    }
    column[rows] <- values[, k]
    x$dense[[sourceIndex]] <- column
  }
  x
}
