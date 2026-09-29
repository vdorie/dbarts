# Constructor and methods for the sparseFactor class (defined in
# R/A_class.R): an unordered factor stored sparsely, entries at the given
# positions carrying x's levels and every other row the implicit reference
# level. A missing value is an explicit stored entry whose value is NA, never
# the reference; a vector whose every row is missing may have no levels, as a
# factor may. It rides to the engine unexpanded as one categorical column
# over its level table, binning bitwise-identically to a dense factor of the
# same values.
#
# The 'length' argument shadows base::length inside the constructor, so
# every internal length is written base::length (the Matrix::sparseVector
# argument convention).
sparseFactor <- function(x, levels, reference, i, length) {
  # a lone NA is logical; it means missing characters, as for factor()
  if (is.logical(x) && all(is.na(x))) {
    x <- as.character(x)
  }
  # resolve the supplied entries to 1-based codes over the level table
  if (is.factor(x)) {
    if (anyNA(base::levels(x))) {
      stop(
        "'x' has an NA level; a sparseFactor stores a missing value, not a ",
        "missing level"
      )
    }
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
  if (anyNA(levels)) {
    stop("'levels' cannot contain NA")
  }
  # a missing value is stored as an NA code; only a value outside 'levels'
  # is an error
  if (any(is.na(codes) & !is.na(x))) {
    stop("'x' contains values absent from 'levels'")
  }

  reference <- if (missing(reference)) {
    levels[1L] # the baseline-contrast convention; NA when there are none
  } else {
    as.character(reference)
  }
  if (
    base::length(reference) != 1L ||
      (base::length(levels) > 0L && is.na(reference)) ||
      (!is.na(reference) && reference %not_in% levels) ||
      (base::length(levels) == 0L && !is.na(reference))
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
    keep <- which(is.na(codes) | codes != referenceCode)
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
    keep <- is.na(storedValues) | storedValues != referenceCode
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
  numMissing <- sum(is.na(object@values))
  cat(
    "sparseFactor of length ",
    object@length,
    ", ",
    numStored,
    " stored entr",
    if (numStored == 1L) "y" else "ies",
    if (numMissing > 0L) paste0(" (", numMissing, " missing)"),
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

# a missing value is always a stored entry
methods::setMethod("is.na", "sparseFactor", function(x) {
  result <- logical(x@length)
  result[x@i + 1L] <- is.na(x@values)
  result
})

methods::setMethod("anyNA", "sparseFactor", function(x, recursive = FALSE) {
  anyNA(x@values)
})

# A factor over every level the sparseFactor declares.
sparseFactorToFactor <- function(x) {
  structure(sparseFactorCodes(x), levels = x@levels, class = "factor")
}

methods::setMethod("as.integer", "sparseFactor", function(x, ...) {
  sparseFactorCodes(x)
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
  utils::str(sparseFactorToFactor(shown), ...)
}

# The positions a row index selects, as for a factor: positive, negative,
# zero, logical (recycled) and repeated indices, with NA for an NA index and
# for a position past the end.
sparseFactorPositions <- function(x, i) {
  if (is.character(i) || !is.null(dim(i)) || is.factor(i)) {
    stop("a sparseFactor can be indexed by position only")
  }
  as.integer(seq_len(x@length)[i])
}

# The positions an assignment index names, as for a factor: NAs are kept
# (the caller drops or refuses them), and an index past the end, or a logical
# one longer than the vector, extends it.
sparseFactorAssignPositions <- function(x, i) {
  if (is.character(i) || !is.null(dim(i)) || is.factor(i)) {
    stop("a sparseFactor can be indexed by position only")
  }
  if (is.logical(i)) {
    newLength <- max(x@length, base::length(i))
    positions <- if (base::length(i) == 0L) {
      integer(0L)
    } else {
      seq_len(newLength)[rep_len(i, newLength)]
    }
  } else if (is.numeric(i) && all(i >= 0, na.rm = TRUE)) {
    positions <- as.integer(i[is.na(i) | i >= 1])
    newLength <- max(x@length, positions, na.rm = TRUE)
  } else {
    positions <- as.integer(seq_len(x@length)[i])
    newLength <- x@length
  }
  list(positions = positions, newLength = as.integer(newLength))
}

# Drops the levels no row takes. Storage is sparse in the reference level, so
# when that one is unused the most common remaining level takes over.
dropSparseFactorLevels <- function(x) {
  used <- tabulate(sparseFactorCodes(x), length(x@levels))
  reference <- match(x@reference, x@levels)
  if (!any(used > 0L)) {
    # every row is missing, so no level is used, as for a factor
    return(sparseFactor(
      rep.int(NA_character_, x@length),
      levels = character(0L)
    ))
  }
  if (used[reference] == 0L) {
    reference <- which.max(used)
  }
  sparseFactor(
    as.character(x),
    levels = x@levels[used > 0L],
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
  if (is.character(i) || !is.null(dim(i)) || is.factor(i)) {
    stop("a sparseFactor can be indexed by position only")
  }
  if (is.logical(i)) {
    i <- as.integer(i)
  }
  # base R's own wording, with the index type it names
  kind <- if (is.integer(i)) "<integer>" else "<real>"
  if (!is.na(i) && i == 0) {
    stop("attempt to select less than one element in get1index ", kind)
  }
  if (!is.na(i) && i < 0) {
    stop("invalid negative subscript in get1index ", kind)
  }
  position <- seq_len(x@length)[i]
  if (length(position) != 1L || is.na(position)) {
    stop("subscript out of bounds")
  }
  x[position]
})

# Assignment by position, as for a factor: a label that is not a level is
# stored as NA with a warning, an NA index is dropped for a value of length
# one and refused otherwise, and an index past the end extends the vector
# with NA.
methods::setReplaceMethod(
  "[",
  "sparseFactor",
  function(x, i, j, ..., value) {
    if (!missing(j) || ...length() > 0L) {
      stop("incorrect number of dimensions")
    }
    if (missing(i)) {
      positions <- seq_len(x@length)
      newLength <- x@length
    } else {
      target <- sparseFactorAssignPositions(x, i)
      positions <- target$positions
      newLength <- target$newLength
    }
    labels <- as.character(value)
    if (anyNA(positions)) {
      if (length(labels) != 1L) {
        stop("NAs are not allowed in subscripted assignments")
      }
      positions <- positions[!is.na(positions)]
    }
    if (length(positions) == 0L && newLength == x@length) {
      return(x)
    }
    codes <- integer(0L)
    if (length(positions) > 0L) {
      if (length(labels) == 0L) {
        stop("replacement has length zero")
      }
      multiple <- length(positions) %% length(labels) == 0L
      if (any(labels %not_in% x@levels & !is.na(labels))) {
        warning("invalid factor level, NA generated", call. = FALSE)
      }
      labels <- rep_len(labels, length(positions))
      codes <- match(labels, x@levels)
      if (!multiple) {
        warning(
          "number of items to replace is not a multiple of replacement length",
          call. = FALSE
        )
      }
    }
    # a later assignment to a position wins
    last <- !duplicated(positions, fromLast = TRUE)
    positions <- positions[last]
    codes <- codes[last]
    stored <- x@i + 1L
    keep <- stored %not_in% positions
    rows <- c(stored[keep], positions)
    values <- c(x@values[keep], codes)
    if (newLength > x@length) {
      gap <- setdiff(seq.int(x@length + 1L, newLength), rows)
      rows <- c(rows, gap)
      values <- c(values, rep.int(NA_integer_, length(gap)))
    }
    sparseFactor(
      values,
      levels = x@levels,
      reference = x@reference,
      i = rows,
      length = newLength
    )
  }
)

methods::setReplaceMethod("[[", "sparseFactor", function(x, i, j, ..., value) {
  if (length(i) != 1L) {
    stop("attempt to select more or less than one element")
  }
  if (is.character(i) || !is.null(dim(i)) || is.factor(i)) {
    stop("a sparseFactor can be indexed by position only")
  }
  if (is.logical(i)) {
    i <- as.integer(i)
  }
  if (is.na(i)) {
    stop("attempt to select more than one element in integerOneIndex")
  }
  kind <- if (is.integer(i)) "<integer>" else "<real>"
  if (i == 0) {
    stop("attempt to select less than one element in OneIndex ", kind)
  }
  if (i < 0) {
    stop("attempt to select more than one element in OneIndex ", kind)
  }
  x[as.integer(i)] <- value
  x
})

# Renaming as for a factor: values that repeat merge their levels, and an NA
# drops its level, turning the entries at it NA. When the reference is
# dropped its implicit rows become stored NA entries and the most common
# remaining level takes over.
methods::setReplaceMethod("levels", "sparseFactor", function(x, value) {
  if (is.list(value)) {
    stop("a sparseFactor's levels can be set from a character vector only")
  }
  value <- as.character(value)
  if (length(value) < length(x@levels)) {
    stop("number of levels differs")
  }
  levels <- unique(value[!is.na(value)])
  recode <- match(value, levels)
  reference <- recode[match(x@reference, x@levels)]
  values <- recode[x@values]
  rows <- x@i
  if (is.na(reference)) {
    implicit <- setdiff(seq_len(x@length) - 1L, rows)
    rows <- c(rows, implicit)
    values <- c(values, rep.int(NA_integer_, length(implicit)))
    ordering <- order(rows)
    rows <- rows[ordering]
    values <- values[ordering]
    used <- tabulate(values, length(levels))
    reference <- if (length(levels) == 0L) {
      NA_integer_
    } else if (any(used > 0L)) {
      which.max(used)
    } else {
      1L
    }
  }
  keep <- is.na(values) | values != reference
  newValidated(
    "sparseFactor",
    i = as.integer(rows[keep]),
    values = as.integer(values[keep]),
    levels = levels,
    reference = levels[reference],
    length = x@length
  )
})

# rep, as for a factor: a row subset by the repeated positions
methods::setMethod("rep", "sparseFactor", function(x, ...) {
  x[rep(seq_len(x@length), ...)]
})

droplevels.sparseFactor <- function(x, ...) {
  dropSparseFactorLevels(x)
}

# truncates, or pads with missing values, as for a factor
methods::setReplaceMethod("length", "sparseFactor", function(x, value) {
  if (length(value) != 1L || is.na(value) || value < 0) {
    stop("invalid value")
  }
  x[seq_len(as.integer(value))]
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
  # the first reference any sparseFactor argument holds; a missing one
  # belongs to a vector with no unstored rows
  references <- unlist(
    lapply(parts, function(part) {
      if (methods::is(part, "sparseFactor")) part@reference
    }),
    use.names = FALSE
  )
  references <- references[!is.na(references)]
  reference <- if (length(references) > 0L) references[1L] else levels[1L]
  offset <- 0L
  rows <- integer(0L)
  labels <- character(0L)
  for (part in parts) {
    n <- length(part)
    if (
      methods::is(part, "sparseFactor") &&
        (is.na(part@reference) || part@reference == reference)
    ) {
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
    reference = reference,
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
anyDuplicated.sparseFactor <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(sparseFactorCodes(x), ...)
}

# as na.omit.default does for a vector: the missing entries go and the
# positions dropped ride along as an "omit" na.action
na.omit.sparseFactor <- function(object, ...) {
  omit <- which(is.na(object))
  if (length(omit) == 0L) {
    return(object)
  }
  object <- object[-omit]
  attr(omit, "class") <- "omit"
  attr(object, "na.action") <- omit
  object
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
  isFactorLike <- function(e) methods::is(e, "sparseFactor") || is.factor(e)
  if (
    isFactorLike(e1) &&
      isFactorLike(e2) &&
      !setequal(levels(e1), levels(e2))
  ) {
    stop("level sets of factors are different")
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
  if (is.null(row.names)) {
    return(as.data.frame.vector(x, NULL, optional, ..., nm = nm))
  }
  # as.data.frame.vector takes integer row names only from R 4.3
  if (
    !(is.character(row.names) || is.integer(row.names)) ||
      length(row.names) != length(x)
  ) {
    stop(
      "'row.names' is not a character or integer vector of length ",
      length(x)
    )
  }
  value <- list(x)
  if (!optional) {
    names(value) <- nm
  }
  structure(value, row.names = row.names, class = "data.frame")
}

# counts per level, as for a factor
summary.sparseFactor <- function(object, ...) {
  summary(sparseFactorToFactor(object), ...)
}
