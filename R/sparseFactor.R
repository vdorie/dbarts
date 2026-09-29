# Constructor and methods for the sparseFactor class (defined in
# R/A_class.R): an unordered factor stored sparsely, entries at the given
# positions carrying x's levels and every other row the implicit reference
# level. Accepted as a predictor through the x/y interface only - a bare S4
# column cannot survive model.frame, so the formula path refuses it - and it
# rides to the engine unexpanded as one categorical column over its level
# table, binning bitwise-identically to a dense factor of the same values.
#
# The 'length' argument shadows base::length inside the constructor, so
# every internal length is written base::length (the Matrix::sparseVector
# argument convention).
sparseFactor <- function(x, levels, reference, i, length) {
  # resolve the supplied entries to 1-based codes over the level table
  if (is.factor(x)) {
    if (missing(levels)) {
      levels <- base::levels(x)
      codes <- as.integer(x)
    } else {
      levels <- as.character(levels)
      codes <- match(as.character(x), levels)
    }
  } else if (is.character(x)) {
    levels <- if (missing(levels)) {
      sort(unique(x)) # the factor() convention
    } else {
      as.character(levels)
    }
    codes <- match(x, levels)
  } else if (is.numeric(x)) {
    if (missing(levels)) {
      stop("integer 'x' requires explicit 'levels'")
    }
    levels <- as.character(levels)
    codes <- as.integer(x)
    if (
      any(
        !is.na(x) &
          (x != codes | codes < 1L | codes > base::length(levels))
      )
    ) {
      stop("'x' must hold integer level codes in [1, length(levels)]")
    }
  } else {
    stop("'x' must be a factor, character, or integer level codes")
  }
  if (anyNA(codes)) {
    if (anyNA(x)) {
      stop("missing values are not supported in a sparseFactor")
    }
    stop("'x' contains values absent from 'levels'")
  }

  reference <- if (missing(reference)) {
    levels[1L] # the baseline-contrast convention
  } else {
    as.character(reference)
  }
  if (
    base::length(reference) != 1L ||
      is.na(reference) ||
      reference %not_in% levels
  ) {
    stop("'reference' must be a single element of 'levels'")
  }
  referenceCode <- match(reference, levels)

  if (missing(i)) {
    # x is the dense vector; rows off the reference level become the
    # stored entries
    n <- if (missing(length)) base::length(codes) else as.integer(length)
    if (base::length(n) != 1L || is.na(n) || n != base::length(codes)) {
      stop("'length' must match the length of a dense 'x'")
    }
    keep <- which(codes != referenceCode)
    rows <- keep - 1L
    storedValues <- codes[keep]
  } else {
    if (missing(length)) {
      stop("'length' is required when 'i' is supplied")
    }
    n <- as.integer(length)
    if (base::length(n) != 1L || is.na(n) || n < 0L) {
      stop("'length' must be a single non-negative integer")
    }
    rows <- as.integer(i)
    if (anyNA(rows) || any(i != rows)) {
      stop("'i' must hold integer positions")
    }
    if (base::length(rows) != base::length(codes)) {
      stop("'i' must pair one position with each entry of 'x'")
    }
    if (any(rows < 1L | rows > n)) {
      stop("'i' must hold 1-based positions in [1, length]")
    }
    if (anyDuplicated(rows) > 0L) {
      stop("'i' cannot contain duplicated positions")
    }
    # store ascending 0-based rows; explicit reference-coded entries are
    # the implicit value, so they drop in canonicalization
    ordering <- order(rows)
    rows <- rows[ordering] - 1L
    storedValues <- codes[ordering]
    keep <- storedValues != referenceCode
    rows <- rows[keep]
    storedValues <- storedValues[keep]
  }

  newValidated(
    "sparseFactor",
    i = rows,
    values = as.integer(storedValues),
    levels = levels,
    reference = reference,
    length = n
  )
}

methods::setMethod("show", "sparseFactor", function(object) {
  numStored <- base::length(object@i)
  cat(
    "sparseFactor of length ",
    object@length,
    ", ",
    numStored,
    " stored entr",
    if (numStored == 1L) "y" else "ies",
    "\n",
    "  levels: ",
    paste0(object@levels, collapse = ", "),
    "\n",
    "  reference (implicit): ",
    object@reference,
    "\n",
    sep = ""
  )
  invisible(NULL)
})

# data.frame insertion sizes columns through NROW/length, so the class
# needs its observation count here to ride in a frame at all
methods::setMethod("length", "sparseFactor", function(x) x@length)

# The 1-based level code of every position.
sparseFactorCodes <- function(x) {
  codes <- rep.int(match(x@reference, x@levels), x@length)
  codes[x@i + 1L] <- x@values
  codes
}

# The level label of every position, as a character vector.
methods::setMethod("as.character", "sparseFactor", function(x, ...) {
  x@levels[sparseFactorCodes(x)]
})

methods::setMethod("levels", "sparseFactor", function(x) x@levels)

# no missing values are held, so none are reported
methods::setMethod("is.na", "sparseFactor", function(x) {
  logical(x@length)
})

# as.factor keeps every level, as it does for a factor
methods::setMethod("as.factor", "sparseFactor", function(x) {
  structure(sparseFactorCodes(x), levels = x@levels, class = "factor")
})

methods::setMethod("as.vector", "sparseFactor", function(x, mode = "any") {
  as.vector(as.character(x), mode)
})

# order, sort and rank read the level codes, as for a factor
methods::setMethod("xtfrm", "sparseFactor", function(x) sparseFactorCodes(x))

# format is what print, head and str of a data frame call per column.
format.sparseFactor <- function(x, ...) {
  format(as.character(x), ...)
}

# str shows the same line a factor does
str.sparseFactor <- function(object, ...) {
  shown <- object[seq_len(min(length(object), 100L))]
  utils::str(
    factor(as.character(shown), levels = object@levels),
    ...
  )
}

# The stored positions a row index selects, as for a factor: positive,
# negative, zero, logical (recycled) and repeated indices.
sparseFactorPositions <- function(x, i, extend = FALSE) {
  if (is.character(i) || !is.null(dim(i)) || is.factor(i)) {
    stop("a sparseFactor can be indexed by position only")
  }
  positions <- if (extend && is.numeric(i) && all(i >= 0, na.rm = TRUE)) {
    i[i > 0]
  } else {
    seq_len(x@length)[i]
  }
  if (anyNA(positions)) {
    stop(
      "a sparseFactor cannot hold NA, so an NA or out-of-range row index ",
      "is refused"
    )
  }
  as.integer(positions)
}

# Drops the levels no row takes. Storage is sparse in the reference level, so
# when that one is unused the most common remaining level takes over.
dropSparseFactorLevels <- function(x) {
  used <- tabulate(sparseFactorCodes(x), length(x@levels))
  if (!any(used > 0L)) {
    return(x)
  }
  reference <- match(x@reference, x@levels)
  if (used[reference] == 0L) {
    reference <- which.max(used)
  }
  levels <- x@levels[used > 0L]
  sparseFactor(
    as.character(x),
    levels = levels,
    reference = x@levels[reference]
  )
}

# Row subset over the same levels and reference, mapping the stored positions
# rather than densifying; drop = TRUE drops unused levels as a factor does.
methods::setMethod("[", "sparseFactor", function(x, i, j, ..., drop = FALSE) {
  if (!missing(j) || ...length() > 0L) {
    stop("incorrect number of dimensions")
  }
  result <- if (missing(i)) {
    x
  } else {
    subsetSparseFactorRows(x, sparseFactorPositions(x, i))
  }
  if (isTRUE(drop)) dropSparseFactorLevels(result) else result
})

methods::setMethod("[[", "sparseFactor", function(x, i, j, ...) {
  if (length(i) != 1L) {
    stop("attempt to select more or less than one element")
  }
  x[i]
})

# Assignment by position: the value's labels must be levels of x, since a
# sparseFactor cannot hold NA, and an index past the end extends the vector
# provided nothing is left unassigned.
methods::setReplaceMethod(
  "[",
  "sparseFactor",
  function(x, i, j, ..., value) {
    if (!missing(j) || ...length() > 0L) {
      stop("incorrect number of dimensions")
    }
    positions <- if (missing(i)) {
      seq_len(x@length)
    } else {
      sparseFactorPositions(x, i, extend = TRUE)
    }
    if (length(positions) == 0L) {
      return(x)
    }
    labels <- as.character(value)
    if (length(labels) == 0L) {
      stop("replacement has length zero")
    }
    labels <- rep_len(labels, length(positions))
    codes <- match(labels, x@levels)
    if (anyNA(codes)) {
      stop("a sparseFactor cannot hold NA or a level it does not have")
    }
    newLength <- max(x@length, positions)
    # a later assignment to a position wins
    last <- !duplicated(positions, fromLast = TRUE)
    positions <- positions[last]
    codes <- codes[last]
    if (newLength > x@length) {
      unassigned <- setdiff(seq.int(x@length + 1L, newLength), positions)
      if (length(unassigned) > 0L) {
        stop(
          "a sparseFactor cannot hold NA, so it cannot be extended past a gap"
        )
      }
    }
    stored <- x@i + 1L
    keep <- stored %not_in% positions
    rows <- c(stored[keep], positions)
    values <- c(x@values[keep], codes)
    sparseFactor(
      values,
      levels = x@levels,
      reference = x@reference,
      i = rows,
      length = newLength
    )
  }
)

methods::setReplaceMethod("length", "sparseFactor", function(x, value) {
  stop("the length of a sparseFactor cannot be set")
})

# Combining takes the union of the levels in order of appearance, as c does
# for factors, and stays sparse over the first argument's reference.
methods::setMethod("c", "sparseFactor", function(x, ...) {
  parts <- c(list(x), list(...))
  isFactorLike <- vapply(
    parts,
    function(part) is.factor(part) || methods::is(part, "sparseFactor"),
    FALSE
  )
  if (!all(isFactorLike)) {
    stop("a sparseFactor can be combined with factors and sparseFactors only")
  }
  levels <- unique(unlist(lapply(parts, levels), use.names = FALSE))
  offset <- 0L
  rows <- integer(0L)
  labels <- character(0L)
  for (part in parts) {
    n <- length(part)
    if (methods::is(part, "sparseFactor") && part@reference == x@reference) {
      rows <- c(rows, part@i + 1L + offset)
      labels <- c(labels, part@levels[part@values])
    } else {
      rows <- c(rows, seq_len(n) + offset)
      labels <- c(labels, as.character(part))
    }
    offset <- offset + n
  }
  sparseFactor(
    labels,
    levels = levels,
    reference = x@reference,
    i = rows,
    length = offset
  )
})

# unique and duplicated compare labels, as they do for a factor
unique.sparseFactor <- function(x, incomparables = FALSE, ...) {
  x[!duplicated(sparseFactorCodes(x))]
}
duplicated.sparseFactor <- function(x, incomparables = FALSE, ...) {
  duplicated(sparseFactorCodes(x), ...)
}

# Only == and != mean anything for unordered factors; the rest give NA with a
# warning, as they do there.
sparseFactorOps <- function(e1, e2) {
  generic <- .Generic # nolint: object_usage_linter.
  if (generic %not_in% c("==", "!=")) {
    warning(gettextf("%s not meaningful for factors", sQuote(generic)))
    return(rep.int(NA, max(length(e1), length(e2))))
  }
  labels <- function(e) {
    if (methods::is(e, "sparseFactor") || is.factor(e)) as.character(e) else e
  }
  get(generic, mode = "function")(labels(e1), labels(e2))
}
methods::setMethod("Ops", signature("sparseFactor", "ANY"), sparseFactorOps)
methods::setMethod("Ops", signature("ANY", "sparseFactor"), sparseFactorOps)
methods::setMethod(
  "Ops",
  signature("sparseFactor", "sparseFactor"),
  sparseFactorOps
)

# data.frame(sf = x) and as.data.frame(x) need this to take the column in.
as.data.frame.sparseFactor <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  nm = deparse1(substitute(x))
) {
  force(nm)
  as.data.frame.vector(x, row.names, optional, ..., nm = nm)
}

# counts per level, as for a factor
summary.sparseFactor <- function(object, ...) {
  summary(as.factor(object), ...)
}
