# The exported augmentation helpers: the per-observation draws a response
# family's engine site runs every sweep, made callable from R against R's own
# random number stream. 'fit' is the location
# WITHOUT the offset - what $getFitsWithoutOffset() reports - and the bridge
# forms the linear predictor as fit + offset.

augFamilies <- c("probit", "logistic", "ordinal", "aft", "nbinom", "student")

# A finite numeric vector, optionally of a required length ('reference' names
# the argument that length came from); counts = TRUE additionally requires
# positive whole numbers, the rule the logistic weights follow, each being a
# sum of that many independent PG(1, psi) draws.
augVector <- function(x, name, n = NULL, reference = "fit", counts = FALSE) {
  if (!is.numeric(x) || length(x) == 0L || !all(is.finite(x))) {
    stop(sprintf("'%s' must be a finite, non-empty numeric vector", name))
  }
  if (!is.null(n) && length(x) != n) {
    stop(sprintf("'%s' must have length %d, that of '%s'", name, n, reference))
  }
  if (counts && (any(x <= 0) || any(x != round(x)))) {
    stop(sprintf("'%s' are observation counts: positive whole numbers", name))
  }
  as.double(x)
}

# A single positive number; whole = TRUE additionally requires an integer, the
# rule the engine's integer-shape Polya-Gamma reduction rounds to anyway.
augScalar <- function(x, name, whole = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a single positive finite number", name))
  }
  if (whole && x != round(x)) {
    stop(sprintf("'%s' must be a whole number", name))
  }
  as.double(x)
}

# Each optional argument belongs to a named family: it is refused BY NAME
# elsewhere rather than ignored, a caller passing 'weights' to a probit draw
# having the wrong model in mind and not a spare argument, and its own family
# REQUIRES it, the draw law carrying no default to fall back on.
augRestrict <- function(supplied, name, family, families) {
  if (supplied && !family %in% families) {
    stop(sprintf(
      "'%s' applies only to family %s, not \"%s\"",
      name,
      paste0("\"", families, "\"", collapse = " and "),
      family
    ))
  }
  if (!supplied && family %in% families) {
    stop(sprintf("family \"%s\" requires '%s'", family, name))
  }
  invisible(NULL)
}

naIfNull <- function(x) if (is.null(x)) NA_real_ else x

dbartsDrawLatents <- function(
  family,
  fit,
  y,
  weights = NULL,
  offset = NULL,
  sigma = NULL,
  dispersion = NULL,
  cutpoints = NULL,
  df = NULL
) {
  family <- match.arg(family, augFamilies)
  fit <- augVector(fit, "fit")
  n <- length(fit)
  y <- augVector(y, "y", n)

  # 'weights' and 'sigma' are the two that are optional WITHIN their families
  if (!is.null(weights)) {
    augRestrict(TRUE, "weights", family, "logistic")
  }
  if (!is.null(sigma)) {
    augRestrict(TRUE, "sigma", family, c("aft", "student"))
  }
  augRestrict(!is.null(dispersion), "dispersion", family, "nbinom")
  augRestrict(!is.null(cutpoints), "cutpoints", family, "ordinal")
  augRestrict(!is.null(df), "df", family, "student")

  if (!is.null(weights)) {
    weights <- augVector(weights, "weights", n, counts = TRUE)
  }
  if (!is.null(offset)) {
    offset <- augVector(offset, "offset", n)
  }
  if (!is.null(dispersion)) {
    dispersion <- augScalar(dispersion, "dispersion", whole = TRUE)
  }
  if (!is.null(df)) {
    df <- augScalar(df, "df")
  }
  if (!is.null(sigma)) {
    sigma <- augScalar(sigma, "sigma")
  }
  if (!is.null(cutpoints)) {
    cutpoints <- augVector(cutpoints, "cutpoints")
    if (is.unsorted(cutpoints, strictly = TRUE)) {
      stop("'cutpoints' must be strictly increasing")
    }
  }

  .Call(
    C_dbarts_bartcore_drawLatents,
    family,
    fit,
    y,
    weights,
    offset,
    naIfNull(sigma),
    naIfNull(dispersion),
    cutpoints,
    naIfNull(df)
  )
}

dbartsWorkingResponse <- function(
  family,
  latent,
  y,
  weights = NULL,
  offset = NULL,
  dispersion = NULL
) {
  family <- match.arg(family, augFamilies)
  latent <- augVector(latent, "latent")
  n <- length(latent)
  y <- augVector(y, "y", n, "latent")

  if (!is.null(weights)) {
    augRestrict(TRUE, "weights", family, "logistic")
  }
  augRestrict(!is.null(dispersion), "dispersion", family, "nbinom")

  if (family %in% c("logistic", "nbinom") && any(latent <= 0)) {
    stop(
      "'latent' is a precision the working response divides by for family ",
      "\"",
      family,
      "\": it must be positive"
    )
  }
  if (!is.null(weights)) {
    weights <- augVector(weights, "weights", n, "latent", counts = TRUE)
  }
  if (!is.null(offset)) {
    offset <- augVector(offset, "offset", n, "latent")
  }
  if (!is.null(dispersion)) {
    dispersion <- augScalar(dispersion, "dispersion", whole = TRUE)
  }

  .Call(
    C_dbarts_bartcore_workingResponse,
    family,
    latent,
    y,
    weights,
    offset,
    naIfNull(dispersion)
  )
}
