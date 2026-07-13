# Constructor and methods for the sparseFactor class (defined in
# R/A_class.R): an unordered factor stored sparsely, entries at the given
# positions carrying x's levels and every other row the implicit reference
# level. Per decision S-CAT the class is recognized by ingestion but
# refused at data construction until the CSC-categorical engine path
# lands; it exists so data sets can be assembled ahead of engine support.
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
