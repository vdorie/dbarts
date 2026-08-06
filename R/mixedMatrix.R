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

## The k-th column of a column-major replacement block of numColumns columns
## over numObservations rows. A single-column call hands back the supplied
## vector itself, so the reference install below stays a pointer swap.
predictorColumnSlice <- function(values, k, numColumns, numObservations) {
  if (numColumns == 1L) {
    values
  } else {
    values[seq.int(
      (k - 1L) * numObservations + 1L,
      length.out = numObservations
    )]
  }
}

## Replace one column of a dgCMatrix in place by DIRECT SLOT SURGERY on
## @i/@p/@x: splice the column's entries out and the replacement's in, shifting
## the pointers of every column after it. Matrix's `[<-` is disqualified here -
## it runs drop0 MATRIX-WIDE, so writing one column strips the explicit zeros
## every OTHER column holds, and a categorical column whose reference level is
## not levels[1] stores code 0 explicitly. rank is the 1-based column index
## within the sparse block; implicit is the value the column's unstored rows
## read (numeric zero for an ordinal column, the reference level's code for a
## categorical one), so the stored entries are the cells that differ from it,
## missing values included - the engine's own pattern rule
## (mutateCscColumnFromDense, src/bartcore/data.hpp), mirrored.
replaceSparseColumn <- function(sparse, rank, implicit, values) {
  stored <- is.na(values) | values != implicit
  newRows <- as.integer(which(stored) - 1L)
  newValues <- as.double(values[stored])

  numEntries <- length(sparse@i)
  start <- sparse@p[rank]
  end <- sparse@p[rank + 1L]
  prefix <- seq_len(start)
  suffix <- if (end < numEntries) seq.int(end + 1L, numEntries) else integer(0)

  sparse@i <- c(sparse@i[prefix], newRows, sparse@i[suffix])
  sparse@x <- c(sparse@x[prefix], newValues, sparse@x[suffix])
  shifted <- seq.int(rank + 1L, length(sparse@p))
  sparse@p[shifted] <-
    sparse@p[shifted] + (length(newRows) - (end - start))
  # any cached factorization describes the pre-change values
  sparse@factors <- list()
  sparse
}

## Replace SEVERAL columns of a dgCMatrix in ONE pass over @i/@p/@x, the plural
## form of replaceSparseColumn: splicing them one at a time would rebuild the
## whole entry vector per column, so a whole-matrix replacement would cost
## O(p * nnz) rather than O(nnz). ranks are the 1-based column indices within
## the sparse block, implicits the value each one's unstored rows read, and
## values a list of replacement vectors parallel to ranks. Matrix's `[<-` is
## disqualified here for the same matrix-wide drop0 reason the singular form
## records.
replaceSparseColumns <- function(sparse, ranks, implicits, values) {
  numColumns <- length(sparse@p) - 1L
  entryRows <- vector("list", numColumns)
  entryValues <- vector("list", numColumns)
  counts <- integer(numColumns)
  for (k in seq_len(numColumns)) {
    index <- match(k, ranks)
    if (is.na(index)) {
      # untouched: carry the column's stored entries across verbatim
      start <- sparse@p[k]
      end <- sparse@p[k + 1L]
      keep <- if (end > start) seq.int(start + 1L, end) else integer(0)
      entryRows[[k]] <- sparse@i[keep]
      entryValues[[k]] <- sparse@x[keep]
    } else {
      column <- values[[index]]
      stored <- is.na(column) | column != implicits[index]
      entryRows[[k]] <- as.integer(which(stored) - 1L)
      entryValues[[k]] <- as.double(column[stored])
    }
    counts[k] <- length(entryRows[[k]])
  }
  sparse@i <- as.integer(unlist(entryRows, use.names = FALSE))
  sparse@x <- as.double(unlist(entryValues, use.names = FALSE))
  sparse@p <- as.integer(c(0L, cumsum(counts)))
  # any cached factorization describes the pre-change values
  sparse@factors <- list()
  sparse
}

## Install a block of predictor columns into a sampler's stored source BY
## REFERENCE (the reference-install mutation semantics,
## docs/design/data-ownership.md, "Mutation: reference-install") - keeps
## data@x current under
## mutation without adding an R-side copy on top of R's own copy-on-write.
## rows = NULL replaces columns whole: a dense-backed slot j is repointed
## straight at the supplied vector (or, for a multi-column call, the caller's
## per-column slice of it), so an unmutated column stays shared with the prior
## container and only the O(p) list spine is duplicated; a factor source
## column decays to a plain code vector automatically, since the old column
## is discarded rather than mutated in place. A CSC-backed column instead has
## its entries spliced into the container's sparse block (replaceSparseColumn),
## which is O(nnz); a mixed call naming both kinds dispatches per column on the
## sign of the map. rows given (a partial merge)
## starts from the old column and overwrites only the addressed rows - the
## merge is inherent for a freshly supplied vector (O(n) for one column, no
## worse than before); the O(spine) win there needs the caller's own
## already-merged vector reinstalled in place. A matrix source has no
## columnar spine to share, so it keeps copy-modify. Per-observation mutation
## of a CSC-backed column is refused upstream, so the partial merge below only
## ever addresses dense-backed columns. Naming several sparse-backed columns -
## a whole-matrix replacement of a sparse design does - splices them in one
## pass (replaceSparseColumns).
installPredictorColumns <- function(x, rows, columns, values) {
  if (is.matrix(x)) {
    if (is.null(rows)) {
      x[, columns] <- values
    } else {
      x[rows, columns] <- values
    }
    return(x)
  }
  values <- as.double(values)
  # a pure dgCMatrix design (the sparse-column in-place mutation extension,
  # docs/design/sparse-columns.md): every column is CSC-backed and ordinal, so
  # the implicit rows read numeric zero
  if (inherits(x, "dgCMatrix")) {
    numObservations <- nrow(x)
    if (length(columns) == 1L) {
      return(replaceSparseColumn(x, columns, 0, values))
    }
    columnValues <- lapply(
      seq_along(columns),
      function(k) {
        predictorColumnSlice(values, k, length(columns), numObservations)
      }
    )
    return(replaceSparseColumns(
      x,
      columns,
      rep(0, length(columns)),
      columnValues
    ))
  }
  numObservations <- x$numObservations
  if (is.null(rows)) {
    sparseRanks <- integer(0)
    sparseImplicits <- numeric(0)
    sparseValues <- list()
    for (k in seq_along(columns)) {
      column <-
        predictorColumnSlice(values, k, length(columns), numObservations)
      source <- x$map[columns[k]]
      if (source > 0L) {
        x$dense[[source]] <- column
      } else {
        rank <- -source
        sparseRanks <- c(sparseRanks, rank)
        sparseImplicits <- c(
          sparseImplicits,
          if (!is.na(x$sparseReference[rank])) x$sparseReference[rank] else 0
        )
        sparseValues[[length(sparseValues) + 1L]] <- column
      }
    }
    if (length(sparseRanks) == 1L) {
      x$sparse <- replaceSparseColumn(
        x$sparse,
        sparseRanks,
        sparseImplicits,
        sparseValues[[1L]]
      )
    } else if (length(sparseRanks) > 1L) {
      x$sparse <- replaceSparseColumns(
        x$sparse,
        sparseRanks,
        sparseImplicits,
        sparseValues
      )
    }
    return(x)
  }
  values <- matrix(values, ncol = length(columns))
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
