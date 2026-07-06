## Internal container for predictor sets mixing dense and sparse sources,
## built when a data frame holds Matrix::sparseVector or dgCMatrix columns
## alongside ordinary ones: the dense columns as one matrix, the sparse
## columns cbound into one dgCMatrix, and a map from overall column to its
## source - positive k for dense column k, negative -k for sparse column k.
## The model-matrix builders attach "varTypes"/"factor.levels"/
## "term.labels"/"drop" exactly as they do on matrices. A plain R list
## underneath, so serialization and sampler re-creation work unchanged.

isSparseDataFrameColumn <- function(column) {
  isS4(column) &&
    (methods::is(column, "sparseVector") ||
      methods::is(column, "sparseMatrix"))
}

## The 0-based row indices and values of each predictor column a sparse
## data-frame column contributes - a sparseVector contributes one, a
## dgCMatrix its columns. Missing entries are stored NaNs, the Matrix
## convention.
sparseColumnSlices <- function(column, name, numObservations) {
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
    columnNames = unlist(blockNames)
  )
  class(result) <- "dbartsMixedMatrix"
  result
}

dim.dbartsMixedMatrix <- function(x) {
  c(nrow(x$sparse), length(x$map))
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
  if (!is.null(result$dense)) {
    result$dense <- result$dense[i, , drop = FALSE]
  }
  result$sparse <- result$sparse[i, , drop = FALSE]
  result
}

as.matrix.dbartsMixedMatrix <- function(x, ...) {
  result <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
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
