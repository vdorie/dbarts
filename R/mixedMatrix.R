## Internal containers for predictor sets that stay columnar instead of
## materializing one cbound double matrix. Both flavors hold their dense
## columns as a per-column list (one element per predictor column, factors
## keeping their integer codes and everything else held as doubles), told
## apart by the "sparse" field:
##
## - the dense flavor (sparse = NULL): built for a data frame of ordinary
##   columns; the map is all positive, so column j reads dense[[map[j]]].
## - the mixed flavor (sparse = a dgCMatrix): built when a data frame holds
##   Matrix::sparseVector, dgCMatrix, or sparseFactor columns alongside
##   ordinary ones. The map carries positive k for dense column k and
##   negative -k for sparse column k.
##
## Either way the contiguous coded block the engine quantizes is assembled
## transiently - by the bridge inside create/setData calls, and by as.matrix
## for R-side consumers - so no dense block stays resident in the container.
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
## data-frame column contributes - a sparseVector or sparseFactor
## contributes one, a dgCMatrix its columns. Missing entries are stored
## NaNs, the Matrix convention. Per contributed column the result also
## carries "reference" and "K": for a sparseFactor the 0-based level-order
## code of the implicit level and the level count, and NA/0 for an ordinal
## sparse column (whose implicit rows are numeric zero).
sparseColumnSlices <- function(column, name, numObservations) {
  if (methods::is(column, "sparseFactor")) {
    if (length(column) != numObservations) {
      stop(
        "sparse column '",
        name,
        "' must have length equal to the ",
        "number of observations"
      )
    }
    # column@i is already 0-based (unlike a sparseVector's), and values are
    # 1-based level codes; the engine reads 0-based codes, so subtract one.
    # The reference level's own level-order code is the implicit rows' code.
    return(list(
      i = list(as.integer(column@i)),
      x = list(as.double(column@values) - 1),
      names = name,
      reference = match(column@reference, column@levels) - 1L,
      K = length(column@levels)
    ))
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
      names = name,
      reference = NA_integer_,
      K = 0L
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
  list(
    i = i,
    x = values,
    names = columnNames,
    reference = rep(NA_integer_, ncol(column)),
    K = rep(0L, ncol(column))
  )
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
    } else if (is.matrix(columns[[j]])) {
      ncol(columns[[j]])
    } else {
      1L
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

  # the dense columns as a per-column list (matrices spliced per column,
  # factors and vectors kept whole), aligned with the map's positive values
  dense <- if (numDense > 0L) {
    denseList <- vector("list", numDense)
    denseOffset <- 0L
    for (j in seq_along(columns)) {
      if (columnIsSparse[j] || widths[j] == 0L) {
        next
      }
      block <- columns[[j]]
      if (is.matrix(block)) {
        for (k in seq_len(ncol(block))) {
          denseList[[denseOffset + k]] <- block[, k]
        }
      } else {
        denseList[[denseOffset + 1L]] <- block
      }
      denseOffset <- denseOffset + widths[j]
    }
    denseList
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

  # per sparse column (aligned with the assembled dgCMatrix columns), the
  # reference level's 0-based code and the level count K; both carry the
  # ordinal sparse columns' NA/0, and the bridge reads them only for the
  # CSC-backed columns the varTypes flag categorical
  sparseReference <- unlist(lapply(columns[columnIsSparse], `[[`, "reference"))
  if (is.null(sparseReference)) {
    sparseReference <- integer(0)
  }
  sparseCategoryCount <- unlist(lapply(columns[columnIsSparse], `[[`, "K"))
  if (is.null(sparseCategoryCount)) {
    sparseCategoryCount <- integer(0)
  }

  result <- list(
    dense = dense,
    sparse = sparse,
    map = map,
    sparseReference = as.integer(sparseReference),
    sparseCategoryCount = as.integer(sparseCategoryCount),
    numObservations = as.integer(numObservations),
    columnNames = unlist(blockNames)
  )
  class(result) <- "dbartsMixedMatrix"
  result
}

## Wrap a bare Matrix::dgCMatrix 'test' set as an all-sparse mixed container,
## so it takes exactly the resident path a mixed-container test set's sparse
## columns already do (validateXTest's xTestIsSparseContainer branch, the
## bridge's parseTestContainer) instead of densifying. All columns map
## negative (CSC-sourced); a bare numeric sparse matrix carries no factor
## levels, so it is all-ordinal by construction - validateXTest refuses it
## against a categorical training design before this runs. Column names come
## from x.test itself (NULL when unnamed, matching a plain matrix's
## colnames()), not synthesized, so the unnamed-test position-matching
## warning still fires when appropriate.
wrapSparseTestMatrix <- function(x.test) {
  p <- ncol(x.test)
  result <- list(
    dense = NULL,
    sparse = x.test,
    map = -seq_len(p),
    sparseReference = rep(NA_integer_, p),
    sparseCategoryCount = rep(0L, p),
    numObservations = as.integer(nrow(x.test)),
    columnNames = colnames(x.test)
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
  }
  if (!is.null(result$sparse)) {
    result$sparse <- result$sparse[i, , drop = FALSE]
  }
  result$numObservations <- if (!is.null(result$sparse)) {
    nrow(result$sparse)
  } else {
    length(result$dense[[1L]])
  }
  result
}

as.matrix.dbartsMixedMatrix <- function(x, ...) {
  result <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  denseColumns <- x$map > 0L
  # every dense column is vector-backed (a factor materializing as its
  # zero-based codes); the map's positive value indexes the dense list
  if (any(denseColumns)) {
    for (j in which(denseColumns)) {
      column <- x$dense[[x$map[j]]]
      result[, j] <- if (is.factor(column)) {
        as.double(as.integer(column) - 1L)
      } else {
        column
      }
    }
  }
  # a sparse column's implicit rows are numeric zero, except a
  # sparse-categorical column's, whose implicit rows are the reference
  # code - possibly itself zero, so the fill has to happen before the
  # explicit entries scatter over it, not after (an explicit entry can
  # legitimately be code 0 too)
  if (any(!denseColumns)) {
    p <- x$sparse@p
    rowIndex <- x$sparse@i
    values <- x$sparse@x
    for (j in which(!denseColumns)) {
      k <- -x$map[j]
      implicitValue <- if (!is.na(x$sparseReference[k])) {
        x$sparseReference[k]
      } else {
        0
      }
      result[, j] <- implicitValue
      entries <- seq.int(p[k] + 1L, length.out = p[k + 1L] - p[k])
      if (length(entries) > 0L) {
        result[rowIndex[entries] + 1L, j] <- values[entries]
      }
    }
  }
  result
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
