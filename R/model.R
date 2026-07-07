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
  parentEnv
) {
  matchedCall <- match.call()

  # the prior vocabulary shadows the caller's environment inside these
  # arguments only: bare names like normal(chi(1.25)) resolve here no matter
  # what packages are attached, and nothing is exported under generic names
  evalEnv <- new.env(parent = parentEnv)
  for (name in names(dbartsPriors)) {
    assign(name, dbartsPriors[[name]], envir = evalEnv)
  }
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

  node.prior <- resolveSpec(matchedCall$node.prior)
  if (!is(node.prior, "dbartsNodePrior")) {
    stop("'node.prior' must be a node prior specification; see ?dbartsPriors")
  }
  if (is(node.prior, "dbartsLinearPrior") || is(node.prior, "dbartsGPPrior")) {
    node.prior <- resolveLeafCovariates(node.prior, data)
  }
  node.hyperprior <- resolveNodeHyperprior(node.prior@k, control@binary)

  namedList(tree.prior, resid.prior, node.prior, node.hyperprior)
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
  validObject(prior)
  prior
}

## Turn a normal prior's raw k into the model's node hyperprior: NULL is the
## family default (2 for continuous responses, chi(1.25, Inf) for binary),
## a positive scalar is fixed, and a hyperprior object passes through.
resolveNodeHyperprior <- function(k, binary) {
  if (is.null(k)) {
    k <- if (binary) chi(1.25, Inf) else 2.0
  }
  if (is.numeric(k)) {
    return(new("dbartsFixedHyperprior", k = k))
  }
  if (is(k, "dbartsNodeHyperprior")) {
    return(k)
  }
  stop("'k' must be a positive scalar or a hyperprior specification")
}

num.vars <- numvars <- NULL # R CMD check
cgm <- function(power = 2, base = 0.95, split.probs = NULL) {
  result <- new(
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
    # compatibility with string specifications like "chi(1.25)" or "2"
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
  new("dbartsChiSqPrior", df = df, quantile = quant)
}

fixed <- function(value = 1.0) {
  new("dbartsFixedPrior", value = value)
}

chi <- function(degreesOfFreedom = 1.25, scale = Inf) {
  new("dbartsChiHyperprior", degreesOfFreedom = degreesOfFreedom, scale = scale)
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
  new(
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
