setMethod(
  "initialize",
  "dbartsModel",
  function(
    .Object,
    tree.prior,
    node.prior,
    node.hyperprior,
    resid.prior,
    proposal.probs = c(
      birth_death = 0.5,
      swap = 0.1,
      change = 0.4,
      birth = 0.5
    ),
    node.scale = 0.5,
    family = "auto"
  ) {
    if (
      !missing(tree.prior) &&
        is(tree.prior, "dbartsCGMPrior") &&
        !is.null(tree.prior@splitProbabilitiesSpec)
    ) {
      stop(
        "tree prior split probabilities must be resolved against data; ",
        "pass the prior to a fitting function instead"
      )
    }
    if (
      !missing(node.prior) &&
        (is(node.prior, "dbartsLinearPrior") ||
          is(node.prior, "dbartsGPPrior")) &&
        !is.integer(node.prior@columns)
    ) {
      stop(
        "node prior columns must be resolved against data; ",
        "pass the prior to a fitting function instead"
      )
    }
    if (!missing(tree.prior)) {
      .Object@tree.prior <- tree.prior
    }
    if (!missing(node.prior)) {
      .Object@node.prior <- node.prior
    }
    if (!missing(node.hyperprior)) {
      .Object@node.hyperprior <- node.hyperprior
    }
    if (!missing(resid.prior)) {
      .Object@resid.prior <- resid.prior
    }

    if (is.null(proposal.probs)) {
      proposal.probs <- c(
        birth_death = 0.5,
        swap = 0.1,
        change = 0.4,
        birth = 0.5
      )
    }

    probs <- proposal.probs[c("birth_death", "swap", "change")]
    if (sum(is.na(probs)) == 1L) {
      probs[is.na(probs)] <- 1 - sum(probs[!is.na(probs)])
      names(probs) <- c("birth_death", "swap", "change")
    } else if (all(is.na(probs))) {
      probs <- c(birth_death = 0.5, swap = 0.1, change = 0.4)
    }

    .Object@p.birth_death <- probs[["birth_death"]]
    .Object@p.swap <- probs[["swap"]]
    .Object@p.change <- probs[["change"]]

    probs <- proposal.probs["birth"]
    if (is.na(probs)) {
      probs <- c(birth = 0.5)
    }

    .Object@p.birth <- probs[["birth"]]

    .Object@node.scale <- node.scale
    .Object@family <- family

    validObject(.Object)
    .Object
  }
)

parsePriors <- function(
  control,
  data,
  tree.prior,
  node.prior,
  resid.prior,
  resid.dist,
  monotone = NULL,
  parentEnv
) {
  matchedCall <- match.call()

  # the prior vocabulary shadows the caller's environment inside these
  # arguments only: bare names like normal(chi(1.5)) resolve here no matter
  # what packages are attached, and nothing is exported under generic names
  evalEnv <- new.env(parent = parentEnv)
  for (name in names(dbartsPriors)) {
    assign(name, dbartsPriors[[name]], envir = evalEnv)
  }
  # the residual-distribution vocabulary rides the same environment, so a bare
  # gaussian()/student() in resid.dist resolves without shadowing stats::gaussian
  evalEnv$gaussian <- gaussian
  evalEnv$student <- student
  evalEnv$control <- control
  evalEnv$data <- data
  evalEnv$num.vars <- evalEnv$numvars <- ncol(data@x)

  # a bare constructor name (tree.prior = cgm) means its defaults; a value
  # that is already a prior object passes through
  resolveSpec <- function(expr) {
    result <- eval(expr, evalEnv)
    if (is.function(result)) {
      result <- result()
    }
    result
  }

  tree.prior <- resolveSpec(matchedCall$tree.prior)
  if (!is(tree.prior, "dbartsTreePrior")) {
    stop("'tree.prior' must be a tree prior specification; see ?dbartsPriors")
  }
  tree.prior <- resolveSplitProbabilities(tree.prior, data)
  # BART package startdart convention: hold the Dirichlet updates until the
  # forest is likelihood-informed
  if (is(tree.prior, "dbartsDartPrior") && is.na(tree.prior@update.delay)) {
    tree.prior@update.delay <- as.numeric(control@n.burn %/% 2L)
  }

  resid.prior <- resolveSpec(matchedCall$resid.prior)
  if (!is(resid.prior, "dbartsResidPrior")) {
    stop(
      "'resid.prior' must be a residual prior specification; see ?dbartsPriors"
    )
  }

  resid.dist <- resolveSpec(matchedCall$resid.dist)
  if (!is(resid.dist, "dbartsResidDist")) {
    stop(
      "'resid.dist' must be a residual distribution: gaussian() or student()"
    )
  }

  node.prior <- resolveSpec(matchedCall$node.prior)
  if (!is(node.prior, "dbartsNodePrior")) {
    stop("'node.prior' must be a node prior specification; see ?dbartsPriors")
  }
  if (is(node.prior, "dbartsLinearPrior") || is(node.prior, "dbartsGPPrior")) {
    node.prior <- resolveLeafCovariates(node.prior, data)
  }

  # `monotone` arrives already resolved to a per-column direction vector (or
  # NULL); it only gates the leaf-model refusal and the fixed-k rule here
  if (
    !is.null(monotone) &&
      (is(node.prior, "dbartsLinearPrior") || is(node.prior, "dbartsGPPrior"))
  ) {
    stop(
      "monotone constraints require the constant leaf; they are not ",
      "supported with linear or gp node priors"
    )
  }
  node.hyperprior <- resolveNodeHyperprior(
    node.prior@k,
    control@binary,
    monotone = !is.null(monotone)
  )

  namedList(tree.prior, resid.prior, resid.dist, node.prior, node.hyperprior)
}

## The residual error degrees of freedom the C bridge reads off the model's
## resid.df attribute: NULL for gaussian errors (no attribute; the Gaussian
## law), 0 for an estimated Student-t df, and the fixed value otherwise.
residDistDf <- function(resid.dist) {
  if (!is(resid.dist, "dbartsStudentDist")) {
    return(NULL)
  }
  if (is.na(resid.dist@df)) 0.0 else resid.dist@df
}

## Turn a linear or gp node prior's raw columns specification into 1-based
## model matrix column indices: names match columns exactly, numbers pass
## through as indices. Categorical columns are rejected - their codes are
## unordered, so a linear term or a distance is meaningless; interact
## through splits instead. (Under factors = "indicators" the dummy columns
## are ordinary numeric columns and legal.) A gp prior's lengthscale also
## resolves here: NULL passes through (the median-distance heuristic), a
## scalar recycles per column.
resolveLeafCovariates <- function(prior, data) {
  label <- if (is(prior, "dbartsGPPrior")) "gp" else "linear"
  columns <- prior@columns
  if (is.null(columns) || length(columns) == 0L) {
    stop(label, " node prior requires at least one covariate column")
  }

  # the engine reads raw covariate values from contiguous dense columns; a
  # mixed container serves them for its dense-backed columns only
  if (!is.matrix(data@x) && !inherits(data@x, "dbartsMixedMatrix")) {
    stop(
      label,
      " node priors are not supported with sparse predictor ",
      "matrices"
    )
  }

  columnNames <- colnames(data@x)
  if (is.character(columns)) {
    if (is.null(columnNames)) {
      stop("cannot assign leaf covariates: model matrix has no column names")
    }
    columnIndices <- match(columns, columnNames)
    if (anyNA(columnIndices)) {
      stop(
        "cannot assign leaf covariates: unrecognized column name(s) ",
        paste0("'", columns[is.na(columnIndices)], "'", collapse = ", ")
      )
    }
  } else if (is.numeric(columns)) {
    columnIndices <- as.integer(columns)
    if (
      anyNA(columnIndices) ||
        any(columnIndices < 1L) ||
        any(columnIndices > ncol(data@x))
    ) {
      stop("cannot assign leaf covariates: column indices out of range")
    }
  } else {
    stop(label, " node prior 'columns' must be a character or numeric vector")
  }
  if (anyDuplicated(columnIndices) > 0L) {
    stop("cannot assign leaf covariates: duplicate columns")
  }

  if (any(data@varTypes[columnIndices] == CATEGORICAL_VARIABLE)) {
    stop(
      "leaf covariates must be continuous columns; interact with factors ",
      "through splits instead"
    )
  }
  if (
    inherits(data@x, "dbartsMixedMatrix") &&
      any(data@x$map[columnIndices] < 0L)
  ) {
    stop(
      "leaf covariates must be dense columns; sparse-backed columns hold ",
      "no raw values"
    )
  }

  # the engine's cap; blocks are solved on the stack
  if (length(columnIndices) > 8L) {
    stop("at most 8 leaf covariates are supported")
  }

  prior@columns <- columnIndices
  if (is(prior, "dbartsGPPrior") && !is.null(prior@lengthscale)) {
    lengthscale <- prior@lengthscale
    if (length(lengthscale) == 1L) {
      lengthscale <- rep_len(lengthscale, length(columnIndices))
    }
    if (length(lengthscale) != length(columnIndices)) {
      stop(
        "gp node prior 'lengthscale' must have length 1 or match the ",
        "number of columns"
      )
    }
    prior@lengthscale <- as.double(lengthscale)
  }
  prior
}

## Turn a cgm-family prior's raw split.probs specification into normalized
## per-column probabilities: NULL or a scalar is uniform, a named vector
## assigns by column or term name with an optional ".default", an unnamed
## vector assigns by position.
resolveSplitProbabilities <- function(prior, data) {
  split.probs <- prior@splitProbabilitiesSpec
  if (is.null(split.probs)) {
    return(prior)
  }

  if (length(split.probs) == 1L) {
    # if length 1, we can ignore it
    split.probs <- numeric()
  } else if (!is.null(names(split.probs))) {
    default <- NA_real_
    split.names <- names(split.probs)
    defaultMatch <- split.names %in% ".default"
    if (sum(defaultMatch) > 1L) {
      stop(
        "cannot assign split probabilities: default specified multiple times"
      )
    }
    if (sum(defaultMatch) == 1L) {
      default <- split.probs[[which(defaultMatch)]]
      split.probs <- split.probs[!defaultMatch]
      split.names <- names(split.probs)
    }

    result <- rep(default, ncol(data@x))
    names(result) <- colnames(data@x)

    if (is.null(names(result)) && length(split.names) > 0L) {
      stop(
        "cannot assign split probabilities: model matrix has no column names"
      )
    }

    namesMatch <- match(split.names, names(result))
    result[namesMatch[!is.na(namesMatch)]] <- split.probs[!is.na(namesMatch)]

    split.probs <- split.probs[is.na(namesMatch)]
    split.names <- names(split.probs)

    for (i in seq_along(split.probs)) {
      if (split.names[i] %not_in% attr(data@x, "term.labels")) {
        stop(
          "cannot assign split probabilities: unrecognized variable name '",
          split.names[i],
          "'"
        )
      }
      factorMatch <- which(startsWith(
        names(result),
        paste0(split.names[i], ".")
      ))
      result[factorMatch] <- split.probs[i]
    }

    split.probs <- result
  } else {
    if (length(split.probs) != ncol(data@x)) {
      stop(
        "cannot assign split probabilities: length of input (",
        length(split.probs),
        ") does not equal number of columns in model matrix (",
        ncol(data@x),
        ")"
      )
    }
  }

  if (length(split.probs) > 0L) {
    if (anyNA(split.probs)) {
      if (!is.null(names(split.probs))) {
        stop(
          "cannot assign split probabilities: missing values for columns ",
          paste0(
            paste0("'", names(split.probs)[is.na(split.probs)], "'"),
            collapse = ", "
          )
        )
      } else {
        stop(
          "cannot assign split probabilities: missing values for columns ",
          paste0(which(is.na(split.probs)), collapse = ", ")
        )
      }
    }

    split.probs <- split.probs / sum(split.probs)
    if (all(split.probs == split.probs[1L])) {
      split.probs <- numeric()
    }
  }

  prior@splitProbabilities <- split.probs
  prior@splitProbabilitiesSpec <- NULL
  validateObject(prior)
  prior
}

## Turn a normal prior's raw k into the model's node hyperprior: NULL is the
## family default (2 for continuous responses, chi(1.5, 2) for binary),
## a positive scalar is fixed, and a hyperprior object passes through. Under a
## monotone constraint k is fixed for both families (an unsupplied k resolves
## to 2, the truncated leaf law having no clean chi-k update, monotone.md
## section 6) and a chi hyperprior is refused.
resolveNodeHyperprior <- function(k, binary, monotone = FALSE) {
  if (is.null(k)) {
    k <- if (monotone || !binary) 2.0 else chi(1.5, 2.0)
  } else if (monotone && !is.numeric(k)) {
    stop(
      "a 'k' hyperprior is not supported under a monotone constraint; ",
      "supply a fixed numeric k (the truncated leaf law has no chi-k update)"
    )
  }
  if (is.numeric(k)) {
    return(newValidated("dbartsFixedHyperprior", k = k))
  }
  if (is(k, "dbartsNodeHyperprior")) {
    return(k)
  }
  stop("'k' must be a positive scalar or a hyperprior specification")
}

## Sign of a single monotone direction element: the sign glyphs, the words, or
## an integer in {-1, 0, 1}.
parseMonotoneSign <- function(value) {
  if (is.character(value)) {
    switch(
      tolower(value),
      "+" = 1L,
      "increasing" = 1L,
      "inc" = 1L,
      "-" = -1L,
      "decreasing" = -1L,
      "dec" = -1L,
      "0" = 0L,
      stop(
        "invalid monotone direction '",
        value,
        "'; use '+'/'-', 'increasing'/'decreasing', or +1/-1"
      )
    )
  } else {
    direction <- as.integer(round(as.numeric(value)))
    if (is.na(direction) || direction < -1L || direction > 1L) {
      stop("invalid monotone direction '", value, "'; use -1, 0, or +1")
    }
    direction
  }
}

## Resolve the 'monotone' argument into a per-column direction vector in
## {-1, 0, +1} of length ncol(data@x): a named vector assigns by column or
## expanded-term name, an unnamed vector of length p assigns by position.
## Categorical predictors refuse the constraint (their codes are unordered);
## numeric and ordinal columns are eligible. NULL when no constraint is active.
resolveMonotone <- function(monotone, data) {
  if (is.null(monotone) || length(monotone) == 0L) {
    return(NULL)
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  result <- integer(numColumns)

  monotoneNames <- names(monotone)
  if (!is.null(monotoneNames) && any(nzchar(monotoneNames))) {
    if (is.null(columnNames)) {
      stop(
        "cannot assign monotone constraints: model matrix has no column names"
      )
    }
    for (i in seq_along(monotone)) {
      direction <- parseMonotoneSign(monotone[[i]])
      name <- monotoneNames[i]
      index <- match(name, columnNames)
      if (!is.na(index)) {
        result[index] <- direction
      } else if (name %in% attr(data@x, "term.labels")) {
        expanded <- which(startsWith(columnNames, paste0(name, ".")))
        result[expanded] <- direction
      } else {
        stop(
          "cannot assign monotone constraints: unrecognized variable name '",
          name,
          "'"
        )
      }
    }
  } else {
    if (length(monotone) != numColumns) {
      stop(
        "unnamed 'monotone' must have length ",
        numColumns,
        " (the number of model matrix columns)"
      )
    }
    for (i in seq_len(numColumns)) {
      result[i] <- parseMonotoneSign(monotone[[i]])
    }
  }

  categorical <- which(result != 0L & data@varTypes == CATEGORICAL_VARIABLE)
  if (length(categorical) > 0L) {
    stop(
      "monotone constraints are undefined for categorical predictors: ",
      paste0("'", columnNames[categorical], "'", collapse = ", "),
      "; only numeric and ordered columns are eligible"
    )
  }

  if (all(result == 0L)) {
    return(NULL)
  }
  result
}

# Resolve the `variance` heteroscedastic selector to 1-based model-matrix
# column indices for the variance forest, or NULL for a homoscedastic fit.
# NULL/FALSE -> no variance forest; TRUE or ~. -> every column; a one-sided
# formula, character names, or numeric indices -> that subset (factor terms
# expand to their indicator columns, the resolveMonotone precedent). Returns
# an integer vector; a full-set selection returns every index (the caller may
# elide the mask, which is equivalent).
resolveVarianceColumns <- function(variance, data) {
  if (is.null(variance) || isFALSE(variance)) {
    return(NULL)
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  if (isTRUE(variance)) {
    return(seq_len(numColumns))
  }
  names <- if (inherits(variance, "formula")) {
    all.vars(variance)
  } else if (is.character(variance)) {
    variance
  } else {
    NULL
  }
  if (!is.null(names)) {
    if (length(names) == 0L) {
      return(seq_len(numColumns))
    }
    result <- integer(0)
    for (name in names) {
      index <- match(name, columnNames)
      if (!is.na(index)) {
        result <- c(result, index)
      } else if (name %in% attr(data@x, "term.labels")) {
        result <- c(result, which(startsWith(columnNames, paste0(name, "."))))
      } else {
        stop(
          "cannot resolve variance predictor: unrecognized variable name '",
          name,
          "'"
        )
      }
    }
    return(sort(unique(result)))
  }
  # numeric indices
  index <- as.integer(variance)
  if (anyNA(index) || any(index < 1L) || any(index > numColumns)) {
    stop("variance column indices must be in [1, number of columns]")
  }
  sort(unique(index))
}

num.vars <- numvars <- NULL # R CMD check
cgm <- function(power = 2, base = 0.95, split.probs = NULL) {
  result <- newValidated(
    "dbartsCGMPrior",
    power = power,
    base = base,
    splitProbabilities = numeric(),
    splitProbabilitiesSpec = NULL
  )
  if (length(split.probs) > 0L && !is.numeric(split.probs)) {
    stop("'split.probs' must be numeric")
  }
  result@splitProbabilitiesSpec <- split.probs
  result
}

linear <- function(columns, k = NULL) {
  if (missing(columns)) {
    stop("linear node prior requires 'columns' naming the leaf covariates")
  }
  if (!is.character(columns) && !is.numeric(columns)) {
    stop("linear node prior 'columns' must be a character or numeric vector")
  }
  normalPrior <- normal(k) # reuses its k validation and coercions
  new("dbartsLinearPrior", k = normalPrior@k, columns = columns)
}

gp <- function(columns, k = NULL, lengthscale = NULL, max.leaf.size = 256L) {
  if (missing(columns)) {
    stop("gp node prior requires 'columns' naming the leaf covariates")
  }
  if (!is.character(columns) && !is.numeric(columns)) {
    stop("gp node prior 'columns' must be a character or numeric vector")
  }
  if (
    !is.null(lengthscale) &&
      (!is.numeric(lengthscale) ||
        length(lengthscale) == 0L ||
        anyNA(lengthscale) ||
        any(lengthscale <= 0))
  ) {
    stop("gp node prior 'lengthscale' must be positive")
  }
  max.leaf.size <- coerceOrError(max.leaf.size, "integer")
  if (
    length(max.leaf.size) != 1L || is.na(max.leaf.size) || max.leaf.size < 1L
  ) {
    stop("gp node prior 'max.leaf.size' must be a positive integer")
  }
  normalPrior <- normal(k) # reuses its k validation and coercions
  new(
    "dbartsGPPrior",
    k = normalPrior@k,
    columns = columns,
    lengthscale = if (is.null(lengthscale)) NULL else as.double(lengthscale),
    max.leaf.size = max.leaf.size
  )
}

normal <- function(k = NULL) {
  if (is.character(k)) {
    # compatibility with string specifications like "chi(1.5)" or "2"
    if (startsWith(k, "chi")) {
      kExpr <- parse(text = k)[[1L]]
      if (!is.call(kExpr)) {
        kExpr <- call(as.character(kExpr))
      }
      k <- eval(kExpr, list2env(list(chi = chi), parent = baseenv()))
    } else {
      k <- coerceOrError(k, "numeric")
    }
  }
  if (is.function(k)) {
    k <- k()
  } # normal(chi)
  if (
    !is.null(k) &&
      !is(k, "dbartsNodeHyperprior") &&
      (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0.0)
  ) {
    stop("'k' must be a positive scalar or a hyperprior specification")
  }
  new("dbartsNormalPrior", k = k)
}

chisq <- function(df = 3, quant = 0.9) {
  newValidated("dbartsChiSqPrior", df = df, quantile = quant)
}

fixed <- function(value = 1.0) {
  newValidated("dbartsFixedPrior", value = value)
}

chi <- function(degreesOfFreedom = 1.5, scale = 2.0) {
  newValidated(
    "dbartsChiHyperprior",
    degreesOfFreedom = degreesOfFreedom,
    scale = scale
  )
}

## Residual-distribution constructors (docs/design/robust-errors.md), passed
## as resid.dist to dbarts()/bart()/bart2(). gaussian() is the default error
## law; student(df) selects outlier-robust Student-t errors. Not exported:
## resolved by bare name inside the resid.dist argument, like the priors, so
## nothing shadows stats::gaussian.
gaussian <- function() {
  new("dbartsGaussianDist")
}

student <- function(df = NULL) {
  if (!is.null(df)) {
    if (
      !is.numeric(df) ||
        length(df) != 1L ||
        is.na(df) ||
        !is.finite(df) ||
        df <= 0.0
    ) {
      stop(
        "student residual 'df' must be NULL (estimate the degrees of ",
        "freedom) or a single positive finite number"
      )
    }
  }
  newValidated(
    "dbartsStudentDist",
    df = if (is.null(df)) NA_real_ else as.double(df)
  )
}

dart <- function(
  power = 2,
  base = 0.95,
  a = 0.5,
  b = 1,
  rho = NULL,
  alpha = 1,
  update.alpha = TRUE,
  update.delay = NULL
) {
  newValidated(
    "dbartsDartPrior",
    power = power,
    base = base,
    splitProbabilities = numeric(),
    splitProbabilitiesSpec = NULL,
    a = a,
    b = b,
    rho = if (is.null(rho)) NA_real_ else rho,
    alpha = alpha,
    update.alpha = update.alpha,
    update.delay = if (is.null(update.delay)) {
      NA_real_
    } else {
      as.numeric(update.delay)
    }
  )
}

## The exported face of the prior constructors: one object, so that no
## generic name (normal, chisq, fixed, chi) enters the search path to be
## masked by or to mask another package by attach order. Inside the
## tree.prior/node.prior/resid.prior arguments of the fitting functions the
## same constructors are available by bare name.
dbartsPriors <- list(
  cgm = cgm,
  dart = dart,
  normal = normal,
  linear = linear,
  gp = gp,
  chisq = chisq,
  fixed = fixed,
  chi = chi
)
