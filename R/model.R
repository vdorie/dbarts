## The default tree-proposal mixture: P(birth/death), P(swap), P(change) select
## the structure move, and P(birth) splits birth vs. death within a birth/death
## move. One source for the formal default, the is.null reset, and the all-NA
## fallbacks below.
defaultProposalProbs <- c(
  birth_death = 0.5,
  swap = 0.1,
  change = 0.4,
  birth = 0.5
)

setMethod(
  "initialize",
  "dbartsModel",
  function(
    .Object,
    tree.prior,
    node.prior,
    node.hyperprior,
    resid.prior,
    proposal.probs = defaultProposalProbs,
    node.scale = 0.5,
    prior.scale = NA_real_,
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
      proposal.probs <- defaultProposalProbs
    }

    probs <- proposal.probs[c("birth_death", "swap", "change")]
    if (sum(is.na(probs)) == 1L) {
      probs[is.na(probs)] <- 1 - sum(probs[!is.na(probs)])
      names(probs) <- c("birth_death", "swap", "change")
    } else if (all(is.na(probs))) {
      probs <- defaultProposalProbs[c("birth_death", "swap", "change")]
    }

    .Object@p.birth_death <- probs[["birth_death"]]
    .Object@p.swap <- probs[["swap"]]
    .Object@p.change <- probs[["change"]]

    probs <- proposal.probs["birth"]
    if (is.na(probs)) {
      probs <- defaultProposalProbs["birth"]
    }

    .Object@p.birth <- probs[["birth"]]

    .Object@node.scale <- node.scale
    .Object@prior.scale <- as.double(prior.scale)
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
  multiForest = FALSE,
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
  for (name in names(dbartsResidDists)) {
    assign(name, dbartsResidDists[[name]], envir = evalEnv)
  }
  # both spellings are exposed for the split.probs vocabulary: num.vars is the
  # current name, numvars the backward-compatible alias, so a bare 1 / num.vars
  # or 1 / numvars in a split.probs expression resolves either way
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
  # NULL); it only gates the leaf-model refusal and the fixed-k rule here, as
  # `multiForest` - whether this fit declares per-forest bases - gates the
  # second half of that same rule
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
    monotone = !is.null(monotone),
    multiForest = isTRUE(multiForest)
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
    # a length-1 spec is a uniform scalar: drop it, uniform is the default
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

## The node scale a family that names no calibration takes, in the units its
## own latent scale is stated in: gaussian and aft the response range's 0.5,
## the probit-scale families 3, and the log-odds families that same 3 widened
## by the logistic latent's standard deviation pi / sqrt(3). One function
## rather than a switch per call site, since the multi-forest guard reads it
## too (a "non-default node scale" has always meant "differs from the family
## default"); the C bridge carries the twin it backstops direct-API consumers
## with. Ordinal reuses probit's latent scale (docs/design/ordinal.md section 2,
## scheme C: the K = 2 anchor is probit exactly) and nbinom's psi is a log-odds,
## so it reuses logistic's (docs/design/negative-binomial.md section 1).
defaultNodeScale <- function(family) {
  switch(
    family,
    gaussian = 0.5,
    aft = 0.5,
    probit = 3.0,
    ordinal = 3.0,
    nbinom = pi * sqrt(3.0),
    logistic = pi * sqrt(3.0)
  )
}

## The half-Cauchy median a forest carrying NO basis takes when its `sd` is not
## declared - the magnitude of its scalar amplitude, in units of the same latent
## scale defaultNodeScale states above. FAMILY-AWARE, because the unit is not
## the same quantity in the two cases. Under gaussian and aft the anchor is the
## RESPONSE's own sd, which already contains the signal, and sigma is DRAWN, so
## a prognostic total at a median 2 of them is both a statement about total
## variation and one the sampler can correct; that is bcf's use_muscale
## convention (docs/design/bcf.md). Under the latent families the anchor is the
## LINK's own error sd and sigma is PINNED, so 2 asserts the median prognostic
## signal is twice the noise before any other forest is added, and nothing
## absorbs it: at 1 - one anchor unit, which is also what the Rd can state - the
## induced index prior reproduces the shipped single-forest binary default's
## tail mass to within 3.4 percent (P(p < 0.01 or p > 0.99) 0.2468 against
## 0.2387, where 2 gave 0.3764).
##
## Total over the package's family vocabulary rather than over the three the
## multi-forest path builds, so an unknown family ERRORS instead of returning
## switch()'s invisible NULL: the internal bartcoreBCFSampler route has no
## R-side family gate to borrow (its refusal is the bridge's), and a NULL there
## collapses a length-8 transport vector to length 7. The families that route
## refuses still resolve, so the refusal stays where it is written.
##
## No C twin, unlike defaultNodeScale: applyBCFSpec always receives explicit
## per-forest parameter vectors, so there is no route on which this default
## could be silently taken engine-side and nothing to backstop.
defaultAmplitudePriorScale <- function(family) {
  switch(
    family,
    gaussian = 2.0,
    aft = 2.0,
    probit = 1.0,
    ordinal = 1.0,
    nbinom = 1.0,
    logistic = 1.0,
    stop("no amplitude prior scale is defined for family \"", family, "\"")
  )
}

## Turn a normal prior's raw k into the model's node hyperprior: NULL is the
## family default (2 for continuous responses, chi(1.5, 2) for binary),
## a positive scalar is fixed, and a hyperprior object passes through. Under a
## monotone constraint k is fixed for both families (an unsupplied k resolves
## to 2, the truncated leaf law having no clean chi-k update; see
## docs/design/monotone.md) and a chi hyperprior is refused.
##
## A multi-forest fit fixes it on the same terms and for the same reason: the
## calibration map pins every forest's k at 1 and never updates it, so the
## binary default's chi-k draw would be a hyperprior on a quantity no forest
## reads. Redirecting the DEFAULT rather than refusing it downstream is what
## keeps a plain binary two-forest call silent while an explicitly named chi()
## still refuses by name, and the forced value changes no fitted model.
resolveNodeHyperprior <- function(
  k,
  binary,
  monotone = FALSE,
  multiForest = FALSE
) {
  if (is.null(k)) {
    k <- if (monotone || multiForest || !binary) 2.0 else chi(1.5, 2.0)
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

## Resolve one predictor selector name to its 1-based model-matrix column
## indices: an exact column-name match returns that single index; otherwise a
## bare term label expands to its indicator columns
## (startsWith(columnNames, "<name>.")). Returns NULL when the name is neither a
## column nor a term, so each caller can raise its own diagnostic; a recognized
## term with no indicator columns yields integer(0), not NULL. columnNames may
## be NULL (no match, and expansion yields none).
resolveTermColumns <- function(name, columnNames, termLabels) {
  index <- match(name, columnNames)
  if (!is.na(index)) {
    return(index)
  }
  if (name %in% termLabels) {
    return(which(startsWith(columnNames, paste0(name, "."))))
  }
  NULL
}

## Resolve an interactions()/blocks() column selector -- a character vector of
## names (each an exact column or a bare term expanding to its indicator
## columns, via resolveTermColumns) or a numeric index vector -- to 1-based
## model-matrix indices. `what` names the selector in error messages.
resolveColumnVector <- function(
  cols,
  what,
  columnNames,
  termLabels,
  numColumns
) {
  if (is.character(cols)) {
    if (is.null(columnNames)) {
      stop("cannot resolve ", what, ": model matrix has no column names")
    }
    result <- integer(0)
    for (name in cols) {
      columns <- resolveTermColumns(name, columnNames, termLabels)
      if (is.null(columns)) {
        stop(
          "cannot resolve ",
          what,
          ": unrecognized variable name '",
          name,
          "'"
        )
      }
      result <- c(result, columns)
    }
    result
  } else if (is.numeric(cols)) {
    index <- as.integer(cols)
    if (anyNA(index) || any(index < 1L) || any(index > numColumns)) {
      stop("cannot resolve ", what, ": column indices out of range")
    }
    index
  } else {
    stop("cannot resolve ", what, ": expected column names or indices")
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
      columns <- resolveTermColumns(
        name,
        columnNames,
        attr(data@x, "term.labels")
      )
      if (is.null(columns)) {
        stop(
          "cannot assign monotone constraints: unrecognized variable name '",
          name,
          "'"
        )
      }
      result[columns] <- direction
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
  requestedNames <- if (inherits(variance, "formula")) {
    all.vars(variance)
  } else if (is.character(variance)) {
    variance
  } else {
    NULL
  }
  if (!is.null(requestedNames)) {
    if (length(requestedNames) == 0L) {
      return(seq_len(numColumns))
    }
    result <- integer(0)
    for (name in requestedNames) {
      columns <- resolveTermColumns(
        name,
        columnNames,
        attr(data@x, "term.labels")
      )
      if (is.null(columns)) {
        stop(
          "cannot resolve variance predictor: unrecognized variable name '",
          name,
          "'"
        )
      }
      result <- c(result, columns)
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

## Resolve a column selector - the columns one forest of a multi-forest model
## may split on (docs/design/bcf.md) - to sorted 1-based model-matrix column
## indices, or NULL for an unrestricted forest. Column names resolve against
## colnames(data@x), indices are range-checked; an explicitly empty selection
## is an error rather than a silent full forest. `argument` names the spelling
## the caller used, which is `vars` on a forest() specification and
## `moderators` on the internal two-forest constructor.
resolveModerators <- function(moderators, data, argument = "moderators") {
  if (is.null(moderators)) {
    return(NULL)
  }
  if (length(moderators) == 0L) {
    stop("'", argument, "' is empty; omit it to leave the forest unrestricted")
  }
  if (is.character(moderators)) {
    columnNames <- colnames(data@x)
    if (is.null(columnNames)) {
      stop("'", argument, "' given by name but the design has no column names")
    }
    moderators <- match(moderators, columnNames)
    if (anyNA(moderators)) {
      stop("'", argument, "' name not found in the design's column names")
    }
  } else {
    moderators <- as.integer(moderators)
    if (any(moderators < 1L | moderators > ncol(data@x))) {
      stop("'", argument, "' column index out of range")
    }
  }
  sort(unique(as.integer(moderators)))
}

## The `basis` declarations of a `forests` list, one element per forest and
## NULL where a forest declares none, or NULL when there is no usable list at
## all. Read before the structural validation resolveForests does, so anything
## it cannot make sense of falls through as NULL to be refused there, by name.
##
## The floor is ONE forest, not two: a length-1 list declaring a basis has to
## reach data@bases like every other one, or the declaration is dropped in
## silence and an ordinary single-forest model is fit instead of the model the
## caller wrote. What refuses it is the designed one-forest refusal in
## resolveSamplerSpec, the site the dbartsData(bases = ) route reaches too. A
## length-1 list declaring NO basis still expands to list(NULL), still fails the
## any-non-null gate at both call sites, and still falls through to
## resolveForests' own refusal by name.
forestBasisDeclarations <- function(forests) {
  if (!is.list(forests) || length(forests) < 1L) {
    return(NULL)
  }
  if (!all(vapply(forests, inherits, logical(1L), "dbartsForest"))) {
    return(NULL)
  }
  lapply(forests, function(spec) spec$basis)
}

## Evaluate a `basis` declaration to the vector it names. A one-sided formula
## is evaluated against the data the fit was given and then in its own
## environment, as a model formula's terms are; anything else is already a
## value. `data` is the fitting function's data argument when that is a frame,
## list or environment, and NULL otherwise (the x/y interface, dbartsSpec).
evaluateForestBasis <- function(basis, data = NULL) {
  if (!inherits(basis, "formula")) {
    return(basis)
  }
  if (length(basis) != 2L) {
    stop("a 'basis' formula must be one-sided, as ~ factor(z)")
  }
  if (is.data.frame(data) || is.list(data) || is.environment(data)) {
    eval(basis[[2L]], data, environment(basis))
  } else {
    eval(basis[[2L]], environment(basis))
  }
}

## Expand an evaluated basis to the matrix of columns a forest's amplitudes
## multiply, one amplitude per column. The rule is R's own model-matrix rule -
## a factor expands to its level indicators, one amplitude per level, with no
## reference level dropped, since the forest carries no intercept of its own -
## and a numeric vector or matrix is already those columns. Level ORDER is
## therefore load-bearing: amplitude j scales level j. A two-level factor
## expands to the (1 - z, z) pair whose amplitudes are exactly bcf's (b0, b1).
expandForestBasis <- function(basis) {
  if (is.null(basis)) {
    return(NULL)
  }
  if (is.character(basis)) {
    basis <- factor(basis)
  }
  if (is.logical(basis) && is.null(dim(basis))) {
    basis <- factor(basis, levels = c(FALSE, TRUE))
  }
  if (anyNA(basis)) {
    stop("a 'basis' cannot be NA")
  }
  if (is.factor(basis)) {
    if (nlevels(basis) < 2L) {
      stop("a 'basis' factor must have at least two levels")
    }
    codes <- as.integer(basis)
    expanded <- matrix(0, length(codes), nlevels(basis))
    expanded[cbind(seq_along(codes), codes)] <- 1
    return(expanded)
  }
  if (!is.numeric(basis)) {
    stop("a 'basis' must be numeric, a factor, or a character vector")
  }
  basis <- as.matrix(basis)
  storage.mode(basis) <- "double"
  if (ncol(basis) < 1L) {
    stop("a 'basis' must have at least one column")
  }
  if (!all(is.finite(basis))) {
    stop("a 'basis' must be finite")
  }
  basis
}

## Validate the knobs one forest() declares. Only a declared (non-NULL) knob is
## checked; an omitted one keeps its NULL, so the caller can tell "not
## declared" from "declared at the default" - which is what makes the
## top-level-versus-forest-0 ambiguity below detectable.
validateForestKnobs <- function(spec) {
  if (!is.null(spec$n.trees)) {
    numTrees <- suppressWarnings(as.integer(spec$n.trees))
    if (length(numTrees) != 1L || is.na(numTrees) || numTrees < 1L) {
      stop("forest 'n.trees' must be a single integer >= 1")
    }
    spec$n.trees <- numTrees
  }
  for (name in c("power", "sd", "amplitude.prior.variance")) {
    if (!is.null(spec[[name]])) {
      value <- suppressWarnings(as.double(spec[[name]]))
      if (length(value) != 1L || is.na(value) || value <= 0.0) {
        stop("forest '", name, "' must be a single positive number")
      }
      spec[[name]] <- value
    }
  }
  if (!is.null(spec$base)) {
    base <- suppressWarnings(as.double(spec$base))
    if (length(base) != 1L || is.na(base) || base <= 0.0 || base >= 1.0) {
      stop("forest 'base' must be a single number in (0, 1)")
    }
    spec$base <- base
  }
  if (!is.null(spec$update.amplitude)) {
    flag <- suppressWarnings(as.logical(spec$update.amplitude))
    if (length(flag) != 1L || is.na(flag)) {
      stop("forest 'update.amplitude' must be TRUE or FALSE")
    }
    spec$update.amplitude <- flag
  }
  spec
}

## Resolve a `forests` declaration into the per-forest knobs a sampler
## specification carries, or NULL for the single-forest path. Returns a LIST of
## K validated knob lists, one per forest, so nothing downstream is keyed on
## two. Everything the engine cannot honour refuses here, by name, rather than
## being dropped: an amplitude prior only where a basis is, a basis somewhere
## on every forest past the first, and the amplitude knobs only where an
## amplitude exists at all. `interactions` and `blocks` are the fit's own
## top-level arguments, which address the FIRST forest under the same spelling
## a forest() uses, so supplying both is ambiguous rather than layered.
## `hasBasis` is a per-forest logical recording whether a basis reached that
## forest some other way, the dbartsData(bases = ) route being a supported one.
resolveForests <- function(forests, interactions, blocks, hasBasis) {
  if (is.null(forests)) {
    return(NULL)
  }
  if (
    !is.list(forests) ||
      !all(vapply(forests, inherits, logical(1L), "dbartsForest"))
  ) {
    stop("'forests' must be a list of forest() specifications")
  }
  if (length(forests) == 0L) {
    stop("'forests' is empty; omit it to fit a single forest")
  }
  resolved <- lapply(forests, validateForestKnobs)
  numForests <- length(resolved)

  for (index in seq_len(numForests)) {
    spec <- resolved[[index]]
    excused <- length(hasBasis) >= index && hasBasis[index]
    if (
      is.null(spec$basis) && !excused && !is.null(spec$amplitude.prior.variance)
    ) {
      stop(
        "'amplitude.prior.variance' is the prior on a basis forest's ",
        "amplitudes, and forest ",
        index,
        " has no 'basis'"
      )
    }
    if (index >= 2L && is.null(spec$basis) && !excused) {
      stop(
        "forest ",
        index,
        " needs a 'basis': the amplitudes multiplying it are what ",
        "distinguishes it from the first"
      )
    }
  }

  first <- resolved[[1L]]
  if (numForests == 1L) {
    amplitudeKnobs <- if (any(hasBasis)) {
      character(0L)
    } else {
      c("sd", "update.amplitude")
    }
    for (name in amplitudeKnobs) {
      if (!is.null(first[[name]])) {
        stop(
          "'",
          name,
          "' configures the amplitudes a model combines its forests with; a ",
          "single-forest 'forests' has none"
        )
      }
    }
  }
  if (!is.null(first$interactions) && !is.null(interactions)) {
    stop(
      "'interactions' is declared both at the top level and on the first ",
      "forest, which are the same constraint; give one"
    )
  }
  if (!is.null(first$blocks) && !is.null(blocks)) {
    stop(
      "'blocks' is declared both at the top level and on the first forest, ",
      "which are the same constraint; give one"
    )
  }
  resolved
}

## The eight doubles attr(control, "bartcore.bcf")$params carries FOR EACH
## FOREST, in the order the C bridge reads them: the forest's tree count and
## structure prior, the node-scale factor and divisor the calibration map
## reads, the amplitude prior's variance and half-Cauchy scale, and the
## amplitude update flag. Forest 1 takes its tree count and structure prior
## from the fit's own control/tree.prior instead, so its first three are
## carried but unread.
##
## Which of the two magnitude channels a forest's `sd` reaches is decided by
## whether it carries a BASIS. A forest WITHOUT one has a plain scalar
## amplitude under a half-Cauchy scale mixture, so `sd` is that mixture's
## median and the node scale stays at the calibration map's anchor - bcf's
## prognostic forest. A forest WITH one has a fixed-variance amplitude block,
## so `sd` multiplies the node scale, divided through the half-normal median
## 0.674 so the amplitude spread times that scale sits at `sd` anchors - bcf's
## basis forest. The ANCHOR is the family's own latent scale: sd(y) under
## gaussian, 1 under probit and pi/sqrt(3) under logistic, per unit of basis
## row norm.
## Both defaults are stated here rather than engine-side, so a `forests`
## declaring nothing beyond its basis and a data object carrying bases with no
## declaration at all resolve to the same numbers.
##
## Each answers a question the caller did not. The basis-free channel's median
## is the FAMILY's (defaultAmplitudePriorScale above). The fixed-variance
## channel's node scale factor is K-AWARE, sqrt(2/K): the MAP carries no per-K
## renormalization - ForestSpec is per forest and the calibration loop reads no
## count - so at the literal 1 the all-basis index prior grew as
## 1.04912 sqrt(K) s without bound, making the prior on the combined location
## depend on how the caller DECOMPOSED the mean rather than on what they said
## about it, the same model written as two forests or as four differing by
## sqrt(2). sqrt(2/K) is that invariance, at K = 2 the identity, so bcf's own
## convention is preserved rather than judged. A DEFAULT and not map algebra: a
## declared `sd` keeps its per-forest reading at every K, and the basis-free
## channel is untouched, a Cauchy having no finite variance to enter the budget
## with - so the all-basis index prior is 1.4837 s at every K while the shipped
## shape's fixed-variance channel is only CAPPED by it, rising from 0.699x of
## the classic k = 2 budget toward 0.989x rather than passing 2x at K = 10.
forestParams <- function(specs, hasBasis, family) {
  declared <- function(value, default) {
    if (is.null(value)) default else value
  }
  nodeScaleDefault <- sqrt(2 / length(specs))
  amplitudeScaleDefault <- defaultAmplitudePriorScale(family)
  lapply(seq_along(specs), function(index) {
    spec <- specs[[index]]
    withBasis <- hasBasis[index]
    as.double(c(
      declared(spec$n.trees, 50L),
      declared(spec$base, 0.25),
      declared(spec$power, 3),
      if (withBasis) declared(spec$sd, nodeScaleDefault) else 1,
      if (withBasis) 0.674 else 1,
      if (withBasis) declared(spec$amplitude.prior.variance, 0.5) else 1,
      if (withBasis) 0 else declared(spec$sd, amplitudeScaleDefault),
      declared(spec$update.amplitude, TRUE)
    ))
  })
}

## Resolve an interactions() specification against the fitted model matrix into
## the engine's per-forest constraint (docs/design/interaction-constraints.md):
## a max-order cap (0 = uncapped) and a de-duplicated 2 x k integer matrix of
## 0-based forbidden co-occurrence pairs. Every validation happens here, at fit
## time, where the column set is known: unknown names, empty groups, a
## max.order below 1, and a forbid/group entry naming a dropped column all
## error. Returns NULL when nothing constrains anything (or no spec is given),
## leaving the availability path byte-for-byte unchanged.
resolveInteractions <- function(interactions, data) {
  if (is.null(interactions)) {
    return(NULL)
  }
  if (!inherits(interactions, "dbartsInteractions")) {
    stop("'interactions' must be an interactions() specification")
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  termLabels <- attr(data@x, "term.labels")

  maxOrder <- 0L
  if (!is.null(interactions$max.order)) {
    order <- suppressWarnings(as.integer(interactions$max.order))
    if (length(order) != 1L || is.na(order) || order < 1L) {
      stop("interactions 'max.order' must be a single integer >= 1")
    }
    maxOrder <- order
  }

  # forbidden pairs accumulate (1-based, unordered) from both forbid and groups
  pairs <- list()

  # forbid: each entry names >= 2 columns barred from sharing any path; a
  # >2-column entry forbids every pair within it
  if (!is.null(interactions$forbid)) {
    forbidList <- interactions$forbid
    if (!is.list(forbidList)) {
      forbidList <- list(forbidList) # a single vector shorthand
    }
    for (entry in forbidList) {
      cols <- unique(resolveColumnVector(
        entry,
        "forbidden interaction",
        columnNames,
        termLabels,
        numColumns
      ))
      if (length(cols) < 2L) {
        stop("each interactions 'forbid' entry must name two or more columns")
      }
      for (i in seq_along(cols)) {
        for (j in seq_len(i - 1L)) {
          pairs[[length(pairs) + 1L]] <- c(cols[i], cols[j])
        }
      }
    }
  }

  # groups: an allow-list. Two NAMED columns may co-occur on a path only if some
  # group holds both; every other pair of named columns is forbidden. Columns
  # named in no group are unconstrained.
  if (!is.null(interactions$groups)) {
    groupList <- interactions$groups
    if (!is.list(groupList)) {
      groupList <- list(groupList)
    }
    resolved <- lapply(groupList, function(group) {
      cols <- unique(resolveColumnVector(
        group,
        "interaction group",
        columnNames,
        termLabels,
        numColumns
      ))
      if (length(cols) == 0L) {
        stop("interactions 'groups' entries must each name at least one column")
      }
      cols
    })
    named <- sort(unique(unlist(resolved)))
    for (i in seq_along(named)) {
      for (j in seq_len(i - 1L)) {
        a <- named[i]
        b <- named[j]
        shareGroup <- any(vapply(
          resolved,
          function(group) a %in% group && b %in% group,
          logical(1)
        ))
        if (!shareGroup) {
          pairs[[length(pairs) + 1L]] <- c(a, b)
        }
      }
    }
  }

  if (maxOrder == 0L && length(pairs) == 0L) {
    return(NULL) # nothing constrains anything
  }

  # 0-based, low-index-first, de-duplicated 2 x k matrix (column-major = the
  # flat pair stream the C bridge reads)
  forbidden <- matrix(0L, nrow = 2L, ncol = 0L)
  if (length(pairs) > 0L) {
    columns <- lapply(pairs, function(pair) as.integer(sort(pair) - 1L))
    forbidden <- do.call(cbind, columns)
    forbidden <- forbidden[, !duplicated(t(forbidden)), drop = FALSE]
    storage.mode(forbidden) <- "integer"
  }

  list(max.order = maxOrder, forbidden = forbidden)
}

## Resolve a blocks() specification against the fitted model matrix into the
## engine's per-tree block-additive constraint (variant A,
## docs/design/interaction-constraints.md): each whole tree is confined to one
## declared group of predictors, so the ensemble is exactly f = sum_G f_G
## (functional ANOVA / grouped GAMI). Unlike interactions(groups=)'s per-path
## allow-list, blocks() lowers to a STATIC per-tree column MASK, so a predictor
## named in no block would be masked out of every tree and go dead; the
## partition must therefore be TOTAL and DISJOINT over the forest's available
## columns, validated here at fit time. Returns a list carrying a 0-based group
## index per predictor (block.of.column; -1 for a column in no block, only when
## availableColumns restricts the forest) and the deterministic per-group tree
## capacity (block.tree.counts, summing to nTrees). NULL when no spec is given.
## availableColumns (1-based) restricts the partitioned set for a moderator or
## variance forest; NULL partitions the full design.
resolveBlocks <- function(blocks, data, nTrees, availableColumns = NULL) {
  if (is.null(blocks)) {
    return(NULL)
  }
  if (!inherits(blocks, "dbartsBlocks")) {
    stop("'blocks' must be a blocks() specification")
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  termLabels <- attr(data@x, "term.labels")

  available <- if (is.null(availableColumns)) {
    seq_len(numColumns)
  } else {
    sort(unique(as.integer(availableColumns)))
  }

  groupList <- blocks$groups
  if (!is.list(groupList)) {
    groupList <- list(groupList) # a single vector shorthand: one block
  }
  if (length(groupList) == 0L) {
    stop("blocks() 'groups' must name at least one group")
  }
  numGroups <- length(groupList)

  # assign each declared column to its group, detecting overlap and columns
  # named outside the forest's available set as we go
  blockOfColumn <- rep(NA_integer_, numColumns)
  for (g in seq_len(numGroups)) {
    cols <- unique(resolveColumnVector(
      groupList[[g]],
      "block group",
      columnNames,
      termLabels,
      numColumns
    ))
    if (length(cols) == 0L) {
      stop("blocks() 'groups' entries must each name at least one column")
    }
    outside <- cols[cols %not_in% available]
    if (length(outside) > 0L) {
      stop(
        "blocks() group ",
        g,
        " names column(s) not among the forest's available predictors: ",
        paste(columnLabels(columnNames, outside), collapse = ", ")
      )
    }
    overlap <- cols[!is.na(blockOfColumn[cols])]
    if (length(overlap) > 0L) {
      stop(
        "blocks() groups must be disjoint; column(s) named in more than one ",
        "group: ",
        paste(columnLabels(columnNames, overlap), collapse = ", ")
      )
    }
    blockOfColumn[cols] <- g
  }

  # totality: every available predictor must be named, else it would be masked
  # out of every tree and go dead
  unnamed <- available[is.na(blockOfColumn[available])]
  if (length(unnamed) > 0L) {
    stop(
      "blocks() must name every predictor exactly once; unassigned column(s): ",
      paste(columnLabels(columnNames, unnamed), collapse = ", ")
    )
  }
  # columns outside the available set carry the -1 (no-block) sentinel
  blockOfColumn[is.na(blockOfColumn)] <- 0L # placeholder; shift below
  blockOfColumn <- blockOfColumn - 1L # 0-based; unavailable columns become -1

  # deterministic per-group tree capacity (consumes NO rng): explicit
  # trees.per.group, or an even split with the first (nTrees mod G) groups
  # getting one extra
  treesPerGroup <- blocks$trees.per.group
  if (is.null(treesPerGroup)) {
    if (nTrees < numGroups) {
      stop(
        "blocks(): the forest has ",
        nTrees,
        " tree(s) but ",
        numGroups,
        " group(s); every block needs at least one tree - use fewer groups, ",
        "more trees, or an explicit 'trees.per.group'"
      )
    }
    baseCount <- nTrees %/% numGroups
    remainder <- nTrees %% numGroups
    counts <- rep.int(baseCount, numGroups)
    if (remainder > 0L) {
      counts[seq_len(remainder)] <- baseCount + 1L
    }
  } else {
    counts <- suppressWarnings(as.integer(treesPerGroup))
    if (length(counts) != numGroups || anyNA(counts)) {
      stop(
        "blocks() 'trees.per.group' must be an integer vector with one entry ",
        "per group (",
        numGroups,
        ")"
      )
    }
    if (any(counts < 1L)) {
      stop("blocks() 'trees.per.group' entries must be positive")
    }
    if (sum(counts) != nTrees) {
      stop(
        "blocks() 'trees.per.group' must sum to the forest's tree count (",
        nTrees,
        "); got ",
        sum(counts)
      )
    }
  }

  list(
    block.of.column = as.integer(blockOfColumn),
    block.tree.counts = as.integer(counts)
  )
}

# name-or-index a column set for an error message
columnLabels <- function(columnNames, cols) {
  if (is.null(columnNames)) {
    return(as.character(cols))
  }
  columnNames[cols]
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

linear <- function(columns, k = NULL, sd = NULL, scale = NULL) {
  if (missing(columns)) {
    stop("linear node prior requires 'columns' naming the leaf covariates")
  }
  if (!is.character(columns) && !is.numeric(columns)) {
    stop("linear node prior 'columns' must be a character or numeric vector")
  }
  # reuses normal()'s k validation and coercions, and its named-scale rules
  normalPrior <- normal(k, sd, scale)
  new(
    "dbartsLinearPrior",
    k = normalPrior@k,
    columns = columns,
    prior.scale = normalPrior@prior.scale,
    prior.sd = normalPrior@prior.sd
  )
}

gp <- function(
  columns,
  k = NULL,
  lengthscale = NULL,
  max.leaf.size = 256L,
  sd = NULL,
  scale = NULL
) {
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
  # reuses normal()'s k validation and coercions, and its named-scale rules
  normalPrior <- normal(k, sd, scale)
  new(
    "dbartsGPPrior",
    k = normalPrior@k,
    columns = columns,
    lengthscale = if (is.null(lengthscale)) NULL else as.double(lengthscale),
    max.leaf.size = max.leaf.size,
    prior.scale = normalPrior@prior.scale,
    prior.sd = normalPrior@prior.sd
  )
}

## One named-calibration argument, wherever it is spelled: NULL or NA leaves it
## unnamed, anything else must be a single positive finite number. NaN is NOT
## the unnamed spelling even though is.na() says so - it carries no intent and
## cannot serve as a divisor - so it is refused here rather than surviving to
## the bridge's own last-line check.
validateNamedScale <- function(value, name) {
  if (is.null(value)) {
    return(NA_real_)
  }
  value <- coerceOrError(value, "numeric")
  if (length(value) != 1L) {
    stop("'", name, "' must be a single number")
  }
  if (is.nan(value) || (!is.na(value) && (!is.finite(value) || value <= 0.0))) {
    stop("'", name, "' must be positive")
  }
  value
}

## The mid-chain setter's value: the same rules, plus a refusal of the NA that
## spells "unnamed" at creation. There is no family default to fall back on
## once a sampler exists, so an absent value is a malformed one.
validateLiveScale <- function(value, name) {
  value <- validateNamedScale(value, name)
  if (is.na(value)) {
    stop("'", name, "' must be a positive finite number")
  }
  value
}

## The two spellings of the named leaf calibration, shared by normal(),
## linear() and gp(): 'scale' is the forest total's prior sd at k = 1, 'sd' the
## same at the resolved k. Both are response units and at most one may be
## given; the sd-to-scale conversion waits for k (resolvePriorScale).
resolveNamedScaleArgs <- function(sd, scale) {
  sd <- validateNamedScale(sd, "sd")
  scale <- validateNamedScale(scale, "scale")
  if (!is.na(sd) && !is.na(scale)) {
    stop("give at most one of 'sd' and 'scale' to a node prior")
  }
  list(prior.sd = sd, prior.scale = scale)
}

## The model's response-unit prior.scale, resolved from a node prior's named
## calibration against the k that will actually be in force. A sampled k has no
## single value to multiply an sd by and drifts every sweep, so the sd spelling
## is refused there rather than honored at the current draw.
resolvePriorScale <- function(node.prior, node.hyperprior) {
  if (is.na(node.prior@prior.sd)) {
    return(node.prior@prior.scale)
  }
  if (!is(node.hyperprior, "dbartsFixedHyperprior")) {
    stop(
      "'sd' names a prior sd at the current 'k', but 'k' is drawn every ",
      "sweep under a hyperprior and the named sd would drift: name ",
      "'scale' instead (the prior sd at k = 1), or fix 'k' at a number"
    )
  }
  node.prior@prior.sd * node.hyperprior@k
}

normal <- function(k = NULL, sd = NULL, scale = NULL) {
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
  named <- resolveNamedScaleArgs(sd, scale)
  new(
    "dbartsNormalPrior",
    k = k,
    prior.scale = named$prior.scale,
    prior.sd = named$prior.sd
  )
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

## Per-forest interaction constraint (docs/design/interaction-constraints.md),
## passed as interactions = to dbarts()/bart2() (and mu.interactions /
## tau.interactions to bcf()). Packages the raw specification; groups and forbid
## resolve against the model matrix, and every value is validated, at fit time
## in resolveInteractions. max.order caps the number of DISTINCT split variables
## on any root-to-leaf path; groups is a co-occurrence allow-list (named columns
## may share a path only with group-mates); forbid names column sets barred from
## sharing a path. Exported (a distinctive top-level name), unlike the priors.
interactions <- function(max.order = NULL, groups = NULL, forbid = NULL) {
  if (is.null(max.order) && is.null(groups) && is.null(forbid)) {
    stop("interactions() needs at least one of 'max.order', 'groups', 'forbid'")
  }
  structure(
    list(max.order = max.order, groups = groups, forbid = forbid),
    class = "dbartsInteractions"
  )
}

## Per-forest block-additive constraint (variant A,
## docs/design/interaction-constraints.md), passed as blocks = to dbarts() /
## bart2() (and mu.blocks / tau.blocks to bcf()). Confines each WHOLE tree to one
## declared group of predictors, so the ensemble is exactly f = sum_G f_G (a
## clean functional-ANOVA / grouped-GAMI decomposition). groups is a list
## partitioning the forest's predictors into disjoint blocks (by model-matrix
## column name - a bare factor term name expands to its indicator columns - or
## index); the partition must be TOTAL, so every predictor is named exactly once
## (a predictor named in no block would be masked out of every tree and go dead).
## trees.per.group optionally fixes how many of the n.trees trees each block
## gets; NULL distributes them as evenly as possible. Everything is validated at
## fit time in resolveBlocks. Exported (a distinctive top-level name), like
## interactions().
blocks <- function(groups, trees.per.group = NULL) {
  if (missing(groups) || is.null(groups)) {
    stop("blocks() requires 'groups', a list partitioning the predictors")
  }
  structure(
    list(groups = groups, trees.per.group = trees.per.group),
    class = "dbartsBlocks"
  )
}

## One forest of a multi-forest model (docs/design/bcf.md), passed inside the
## forests = list(forest(), forest(...)) argument of dbarts()/dbartsSpec().
## Every knob is per forest, so the fitting functions grow exactly one argument
## however many forests a model has. basis is the data the forest's amplitudes
## multiply, a one-sided formula or a vector, expanded by R's own model-matrix
## rule; vars restricts the columns the forest may split on; n.trees, base and
## power are its tree-structure prior; sd is its total's prior scale in units
## of the family's latent scale (sd(y) under gaussian, 1 under probit,
## pi/sqrt(3) under logistic) per unit of basis row norm;
## amplitude.prior.variance is the N(0, .) variance of the amplitudes on
## its basis; update.amplitude fixes those amplitudes at their prior center
## when FALSE; interactions and blocks are this forest's own constraints - the
## arguments of the same names on the fitting function are the FIRST forest's.
## Every knob defaults to NULL, "not declared", which is what lets a
## declaration that collides with one of those arguments refuse rather than
## silently win. Validated at fit time, in resolveForests. Exported, like
## interactions() and blocks().
forest <- function(
  basis = NULL,
  vars = NULL,
  n.trees = NULL,
  base = NULL,
  power = NULL,
  sd = NULL,
  interactions = NULL,
  blocks = NULL,
  amplitude.prior.variance = NULL,
  update.amplitude = NULL
) {
  structure(
    list(
      basis = basis,
      vars = vars,
      n.trees = n.trees,
      base = base,
      power = power,
      sd = sd,
      interactions = interactions,
      blocks = blocks,
      amplitude.prior.variance = amplitude.prior.variance,
      update.amplitude = update.amplitude
    ),
    class = "dbartsForest"
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

## The residual-distribution vocabulary, resolved by bare name inside the
## resid.dist argument the same way dbartsPriors handles the prior names.
dbartsResidDists <- list(
  gaussian = gaussian,
  student = student
)
