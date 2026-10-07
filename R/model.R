## The default tree-proposal mixture: P(birth/death), P(swap), P(change),
## P(perturb) and P(rule_gibbs) select the structure move, and P(birth) splits
## birth vs. death within a birth/death move. Swap ships at zero because at
## production forest sizes it is nearly all no-op, but the move stays
## reachable: it alone rotates a child's rule up the tree, which is how a
## single-tree fit crosses between rootings. Perturb and rule_gibbs ship at
## zero pending their benefit measurements. One source for the formal default,
## the is.null reset, and the all-NA fallbacks below.
defaultProposalProbs <- c(
  birth_death = 0.6,
  swap = 0,
  change = 0.4,
  perturb = 0,
  rule_gibbs = 0,
  birth = 0.5
)

## A monotone forest proposes birth and death only (change and swap would need
## a constrained integral over more than two leaves). Creation and $setControl
## share the rule: a defaulted mixture is rewritten to birth/death-only, and
## any other is refused, except that $setControl lets stand a mixture already
## free of the other moves (the sampler's own, or the frozen all-zero one).
monotoneProposalProbs <- function(proposal.probs, allowBirthDeath = FALSE) {
  if (
    isTRUE(all.equal(
      proposal.probs[names(defaultProposalProbs)],
      defaultProposalProbs
    ))
  ) {
    return(c(
      birth_death = 1,
      swap = 0,
      change = 0,
      perturb = 0,
      rule_gibbs = 0,
      birth = 0.5
    ))
  }
  otherMoves <- proposal.probs[c("swap", "change", "perturb", "rule_gibbs")]
  if (allowBirthDeath && all(otherMoves == 0)) {
    return(proposal.probs)
  }
  stop(
    "'monotone' forces birth/death-only proposals; a non-default ",
    "'proposal.probs' cannot be honored under the constraint"
  )
}

## The moves whose default is a NUMBER rather than a share of what is left:
## they resolve ahead of the fill below and never enter it.
zeroDefaultProposalNames <- c("perturb", "rule_gibbs")

## Resolves a caller's `proposal.probs` into the canonical six-name mixture the
## control slot carries and the bridge reads. NULL is the shipped default; a
## partial vector is filled by the rules below. No branch may leave an NA.
resolveProposalProbs <- function(proposal.probs) {
  if (is.null(proposal.probs)) {
    proposal.probs <- defaultProposalProbs
  }
  ## every entry names its move; an unnamed entry or an unknown name would
  ## otherwise be dropped by the lookups below and the default mix run
  validNames <- names(defaultProposalProbs)
  entryNames <- names(proposal.probs)
  if (length(proposal.probs) > 0L) {
    if (is.null(entryNames) || any(is.na(entryNames) | entryNames == "")) {
      stop(
        "'proposal.probs' must name each of its entries, from ",
        quotedNameList(validNames)
      )
    }
    duplicated <- unique(entryNames[duplicated(entryNames)])
    if (length(duplicated) > 0L) {
      stop(
        "'proposal.probs' names ",
        quotedNameList(duplicated),
        " more than once; the moves are ",
        quotedNameList(validNames)
      )
    }
    unknown <- unique(entryNames[entryNames %not_in% validNames])
    if (length(unknown) > 0L) {
      stop(
        "'proposal.probs' has unknown ",
        if (length(unknown) > 1L) "names " else "name ",
        quotedNameList(unknown),
        "; the moves are ",
        quotedNameList(validNames)
      )
    }
  }
  ## Perturb and rule_gibbs resolve AHEAD of the fill and never enter it.
  ## Their defaults are numbers rather than shares, so an unnamed one is
  ## zero and the residual below is taken against 1 minus their sum.
  ## Widening the fill's name set instead would silently re-resolve every
  ## vector leaving two of the widened set unnamed:
  ## c(birth_death = 0.5, change = 0.4) resolves swap to 0.1 and would
  ## resolve the zero-default moves to it as well.
  zeroDefaults <- vapply(
    zeroDefaultProposalNames,
    function(name) {
      value <- proposal.probs[name]
      if (is.na(value)) {
        value <- defaultProposalProbs[name]
      }
      value[[1L]]
    },
    0.0
  )

  ## The fill, over the three structural names. One unnamed element takes
  ## the residual. Two unnamed, one of them swap, resolve as well: swap is
  ## the one whose default is a number rather than a share, so it takes its
  ## zero and the other takes the residual. Naming only the moves that
  ## default to a number leaves the split between birth/death and change
  ## undetermined and is an error rather than a silent choice, and naming
  ## none of them is the default.
  ##
  ## The residual is a share of structural mass and an all-zero mixture has
  ## none: with birth/death and change both named zero and the zero-default
  ## moves zero there is nothing to distribute, so an unnamed swap keeps its
  ## zero rather than taking the whole of it, and the mixture stays frozen -
  ## no structural proposal is made at all. An unnamed birth/death or change
  ## still takes the residual, so c(change = 0) is birth/death 1 as before.
  zeroDefaultTotal <- sum(zeroDefaults)
  probs <- proposal.probs[c("birth_death", "swap", "change")]
  names(probs) <- c("birth_death", "swap", "change")
  unnamed <- is.na(probs)
  if (all(unnamed) && zeroDefaultTotal == 0) {
    probs <- defaultProposalProbs[c("birth_death", "swap", "change")]
  } else if (unnamed[["birth_death"]] && unnamed[["change"]]) {
    stop(
      "'proposal.probs' names only the zero-default moves 'swap', ",
      "'perturb' and 'rule_gibbs'; name at least one of 'birth_death' ",
      "and 'change'"
    )
  } else {
    if (sum(unnamed) == 2L) {
      probs[["swap"]] <- 0
    }
    named <- probs[!is.na(probs)]
    frozen <- zeroDefaultTotal == 0 &&
      !unnamed[["birth_death"]] &&
      !unnamed[["change"]] &&
      all(named == 0)
    probs[is.na(probs)] <- if (frozen) {
      0
    } else {
      1 - (zeroDefaultTotal + sum(named))
    }
  }

  birth <- proposal.probs["birth"]
  if (is.na(birth)) {
    birth <- defaultProposalProbs["birth"]
  }

  c(
    birth_death = probs[["birth_death"]],
    swap = probs[["swap"]],
    change = probs[["change"]],
    perturb = zeroDefaults[["perturb"]],
    rule_gibbs = zeroDefaults[["rule_gibbs"]],
    birth = birth[[1L]]
  )
}

setMethod(
  "initialize",
  "dbartsModel",
  function(
    .Object,
    tree.prior,
    leaf.prior,
    leaf.hyperprior,
    resid.prior,
    leaf.scale = 0.5,
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
      !missing(leaf.prior) &&
        (is(leaf.prior, "dbartsLinearPrior") ||
          is(leaf.prior, "dbartsGPPrior")) &&
        !is.integer(leaf.prior@columns)
    ) {
      stop(
        "leaf prior columns must be resolved against data; ",
        "pass the prior to a fitting function instead"
      )
    }
    if (!missing(tree.prior)) {
      .Object@tree.prior <- tree.prior
    }
    if (!missing(leaf.prior)) {
      .Object@leaf.prior <- leaf.prior
    }
    if (!missing(leaf.hyperprior)) {
      .Object@leaf.hyperprior <- leaf.hyperprior
    }
    if (!missing(resid.prior)) {
      .Object@resid.prior <- resid.prior
    }

    .Object@leaf.scale <- leaf.scale
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
  leaf.prior,
  resid.prior,
  monotone = NULL,
  multiForest = FALSE,
  kHyperprior = control@binary,
  parentEnv
) {
  matchedCall <- match.call()

  # the prior vocabulary shadows the caller's environment inside these
  # arguments only: bare names like normal(chi(1.5)) resolve here no matter
  # what packages are attached, and nothing is exported under generic names.
  # Both spellings are exposed for the split.probs vocabulary: num.vars is the
  # current name, numvars the backward-compatible alias, so a bare 1 / num.vars
  # or 1 / numvars in a split.probs expression resolves either way
  vocabulary <- c(
    dbartsPriors,
    list(num.vars = ncol(data@x), numvars = ncol(data@x))
  )
  # a bare constructor name (tree.prior = cgm) means its defaults; a value
  # that is already a prior object passes through
  resolveSpec <- function(expr, name, class, label) {
    evalInVocabulary(
      expr,
      vocabulary,
      parentEnv,
      resolvedAs(name, class, paste(label, "specification"))
    )
  }

  tree.prior <- resolveSpec(
    matchedCall$tree.prior,
    "tree.prior",
    "dbartsTreePrior",
    "tree prior"
  )
  tree.prior <- resolveSplitProbabilities(tree.prior, data)
  # BART package startdart convention: hold the Dirichlet updates until the
  # forest is likelihood-informed
  if (is(tree.prior, "dbartsDartPrior") && is.na(tree.prior@update.delay)) {
    tree.prior@update.delay <- as.numeric(control@n.burn %/% 2L)
  }

  resid.prior <- resolveSpec(
    matchedCall$resid.prior,
    "resid.prior",
    "dbartsResidPrior",
    "residual prior"
  )
  leaf.prior <- resolveSpec(
    matchedCall$leaf.prior,
    "leaf.prior",
    "dbartsLeafPrior",
    "leaf prior"
  )
  if (is(leaf.prior, "dbartsLinearPrior") || is(leaf.prior, "dbartsGPPrior")) {
    leaf.prior <- resolveLeafCovariates(leaf.prior, data)
  }

  # `monotone` arrives already resolved to a per-column direction vector (or
  # NULL); it only gates the leaf-model refusal and the fixed-k rule here, as
  # `multiForest` - whether this fit declares per-forest bases - gates the
  # second half of that same rule
  if (
    !is.null(monotone) &&
      (is(leaf.prior, "dbartsLinearPrior") || is(leaf.prior, "dbartsGPPrior"))
  ) {
    stop(
      "monotone constraints require the constant leaf; they are not ",
      "supported with linear or gp leaf priors"
    )
  }
  resolved <- resolveLeafPrior(
    leaf.prior,
    kHyperprior,
    monotone = !is.null(monotone),
    multiForest = isTRUE(multiForest)
  )
  leaf.hyperprior <- resolved$leaf.hyperprior
  prior.scale <- resolved$prior.scale

  namedList(tree.prior, resid.prior, leaf.prior, leaf.hyperprior, prior.scale)
}

## Turn a linear or gp leaf prior's raw columns specification into 1-based
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
    stop(label, " leaf prior requires at least one covariate column")
  }

  # the engine reads raw covariate values from contiguous dense columns; a
  # mixed container serves them for its dense-backed columns only
  if (!is.matrix(data@x) && !inherits(data@x, "dbartsMixedMatrix")) {
    stop(
      label,
      " leaf priors are not supported with sparse predictor ",
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
    columnIndices <- coerceOrError(columns, "integer")
    if (
      anyNA(columnIndices) ||
        any(columnIndices < 1L) ||
        any(columnIndices > ncol(data@x))
    ) {
      stop("cannot assign leaf covariates: column indices out of range")
    }
  } else {
    stop(label, " leaf prior 'columns' must be a character or numeric vector")
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
        "gp leaf prior 'lengthscale' must have length 1 or match the ",
        "number of columns"
      )
    }
    prior@lengthscale <- as.double(lengthscale)
  }
  prior
}

## The entry rules a split.probs specification is held to, as sample() holds
## its 'prob': no negative and no non-finite entry, and at least one positive
## one. A missing entry of a vector is refused ahead of this, by column.
refuseInvalidSplitProbabilities <- function(split.probs) {
  if (any(!is.finite(split.probs) | split.probs < 0)) {
    stop("'split.probs' must be non-negative and finite")
  }
  if (!any(split.probs > 0)) {
    stop("'split.probs' must give at least one column a positive probability")
  }
  invisible(NULL)
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
    # a length-1 spec is a uniform scalar: held to the entry rules like any
    # other, then dropped, uniform being the default
    refuseInvalidSplitProbabilities(split.probs)
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

    # ahead of the normalization, which a negative sum turns into
    # probabilities and a zero or infinite one into NaN
    refuseInvalidSplitProbabilities(split.probs)
    split.probs <- split.probs / sum(split.probs)
    if (all(split.probs == split.probs[1L])) {
      split.probs <- numeric()
    }
  }

  prior@splitProbabilities <- split.probs
  prior@splitProbabilitiesSpec <- NULL
  tryCatch(validObject(prior), error = rethrowValidityError)
  prior
}

## The leaf scale a family that names no calibration takes, in the units its
## own latent scale is stated in: gaussian and aft the response range's 0.5,
## the probit-scale families 3, logistic that same 3 widened by the logistic
## latent's standard deviation pi / sqrt(3), and nbinom 3 on the log mean,
## the anchor its probe preferred with k drawn (dec-B183). One function
## rather than a switch per call site, since the multi-forest guard reads it
## too (a "non-default leaf scale" has always meant "differs from the family
## default"); the C bridge carries the twin it backstops direct-API consumers
## with. Ordinal reuses probit's latent scale (scheme C: the K = 2 anchor is
## probit exactly).
defaultLeafScale <- function(family) {
  switch(
    family,
    gaussian = 0.5,
    aft = 0.5,
    probit = 3.0,
    ordinal = 3.0,
    nbinom = 3.0,
    logistic = pi * sqrt(3.0),
    # the K = 2 pairwise-log-odds anchor: the softmax calibration map owns
    # every category forest's leaf scale, so the engine never reads this value
    # for a multinomial sampler and it is recorded here only so the model
    # object states the anchor the map applies
    multinomial = pi * sqrt(3.0) / sqrt(2.0),
    stop("no leaf scale is defined for family \"", family, "\"")
  )
}

## The half-Cauchy median a forest carrying NO basis takes when its `sd` is
## not declared, in units of the latent scale defaultLeafScale states above.
## FAMILY-AWARE: under gaussian and aft the anchor is the RESPONSE's own sd
## (sigma is DRAWN); under the latent families the anchor is the LINK's own
## error sd (sigma is PINNED), so the anchor unit differs 2:1 between them to
## keep the induced index prior matched to the shipped single-forest binary
## default.
##
## Total over the package's family vocabulary rather than over the three the
## multi-forest path builds: an unknown family ERRORS instead of switch()'s
## invisible NULL, since the multi-forest creation path (dbarts(forests = ))
## borrows this vocabulary rather than declaring a family gate of its own.
##
## No C twin, unlike defaultLeafScale: applyAmplitudeSpec always receives
## explicit per-forest parameter vectors, so there is nothing to backstop.
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

## Refuses a prior object supplied together with a shorthand argument that
## would have built the same prior, naming both. The collision set is data
## because it differs by entry point: bart2 collides tree.prior with
## power/base, xbart deliberately does not, since there they are grid axes
## that override any supplied object every cell (man/xbart.Rd). Presence in
## the matched call, not an explicit-NULL value, is what collides.
## 'supplied' names the shorthands that reached the caller through '...'
## rather than as formals: a retired name is cleared from the matched call
## before anything is forwarded, so the call alone no longer sees it.
refuseColliding <- function(
  matchedCall,
  objectName,
  shorthands,
  supplied = character()
) {
  hit <- shorthands[shorthands %in% c(names(matchedCall), supplied)]
  if (length(hit) > 0L) {
    stop(
      "'",
      objectName,
      "' cannot be combined with '",
      hit[1L],
      "': supply the prior either as an object or through its shorthand ",
      "arguments, not both"
    )
  }
  invisible(NULL)
}

## Turn a normal prior's raw k into the model's leaf hyperprior: NULL is the
## family default (2 for continuous responses, chi(1.5, 2) for binary and
## nbinom, drawsLeafKByDefault; 'binary' carries that flag),
## a positive scalar is fixed, and a hyperprior object passes through. Under a
## monotone constraint k is fixed for both families (an unsupplied k resolves
## to 2, the truncated leaf law having no clean chi-k update) and a chi
## hyperprior is refused.
##
## A multi-forest fit fixes it on the same terms and for the same reason: the
## calibration map pins every forest's k at 1 and never updates it, so the
## binary default's chi-k draw would be a hyperprior on a quantity no forest
## reads. Redirecting the DEFAULT rather than refusing it downstream is what
## keeps a plain binary two-forest call silent while an explicitly named chi()
## still refuses by name, and the forced value changes no fitted model.
resolveLeafHyperprior <- function(
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
  if (is(k, "dbartsLeafHyperprior")) {
    return(k)
  }
  stop("'k' must be a positive scalar or a hyperprior specification")
}

## The monotone direction vocabulary, matched exactly and case-sensitively:
## the words, and the codes as the numbers 1, -1 and 0 or as the strings c()
## makes of them when words and codes share a vector.
MONOTONE_DIRECTION_CODES <- c(
  increasing = 1L,
  decreasing = -1L,
  "1" = 1L,
  "-1" = -1L,
  "0" = 0L
)

## Code of a single monotone direction element, in {-1, 0, 1}.
parseMonotoneSign <- function(value) {
  direction <- NA_integer_
  if (length(value) == 1L && is.character(value)) {
    direction <- unname(MONOTONE_DIRECTION_CODES[value])
  } else if (length(value) == 1L && is.numeric(value)) {
    direction <- match(value, c(-1, 0, 1)) - 2L
  }
  if (is.na(direction)) {
    stop(
      "monotone directions must be one of \"increasing\", \"decreasing\", ",
      "1, -1, 0; got '",
      toString(value),
      "'"
    )
  }
  direction
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
  numColumns,
  argument = what
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
    index <- coerceOrError(cols, "integer", name = argument)
    if (anyNA(index) || any(index < 1L) || any(index > numColumns)) {
      stop("cannot resolve ", what, ": column indices out of range")
    }
    index
  } else {
    stop("cannot resolve ", what, ": expected column names or indices")
  }
}

## Resolve the 'monotone' argument, a monotone() specification or the plain
## direction vector that is shorthand for monotone(directions) at the default
## prior, into list(directions = , prior = ): directions is a per-column
## vector in {-1, 0, +1} of length ncol(data@x), a named vector assigning by
## column or expanded-term name and an unnamed one of length p by position.
## Categorical predictors refuse the constraint (their codes are unordered);
## numeric and ordinal columns are eligible. NULL when no constraint is active.
resolveMonotone <- function(spec, data) {
  if (is.null(spec)) {
    return(NULL)
  }
  if (!inherits(spec, "dbartsMonotone")) {
    spec <- monotone(spec)
  }
  directions <- spec$directions
  if (length(directions) == 0L) {
    return(NULL)
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  result <- integer(numColumns)
  assigned <- logical(numColumns)

  monotoneNames <- names(directions)
  if (!is.null(monotoneNames) && any(nzchar(monotoneNames))) {
    if (!all(nzchar(monotoneNames))) {
      stop("'monotone' directions must be all named or all unnamed")
    }
    if (anyDuplicated(monotoneNames)) {
      stop(
        "'monotone' names a predictor more than once: '",
        monotoneNames[anyDuplicated(monotoneNames)],
        "'"
      )
    }
    if (is.null(columnNames)) {
      stop(
        "cannot assign monotone constraints: model matrix has no column names"
      )
    }
    for (i in seq_along(directions)) {
      direction <- parseMonotoneSign(directions[[i]])
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
      if (any(assigned[columns])) {
        stop(
          "'monotone' names a predictor more than once: '",
          name,
          "' overlaps a column named earlier"
        )
      }
      assigned[columns] <- TRUE
      result[columns] <- direction
    }
  } else {
    if (length(directions) != numColumns) {
      stop(
        "unnamed 'monotone' must have length ",
        numColumns,
        " (the number of model matrix columns)"
      )
    }
    for (i in seq_len(numColumns)) {
      result[i] <- parseMonotoneSign(directions[[i]])
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
  list(directions = result, prior = spec$prior)
}

# Resolve the `variance` heteroscedastic selector to 1-based model-matrix
# column indices for the variance forest, or NULL for a homoscedastic fit.
# NULL/FALSE -> no variance forest; TRUE or ~. -> every column; a one-sided
# formula, character names, or numeric indices -> that subset (factor terms
# expand to their indicator columns, the resolveMonotone precedent). Returns
# an integer vector; a full-set selection returns every index (the caller may
# elide the mask, which is equivalent).
resolveVarianceColumns <- function(variance, data, argument = "variance") {
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
  index <- coerceOrError(variance, "integer", name = argument)
  if (anyNA(index) || any(index < 1L) || any(index > numColumns)) {
    stop("variance column indices must be in [1, number of columns]")
  }
  sort(unique(index))
}

## Resolve a column selector - the columns one forest of a multi-forest model
## may split on - to sorted 1-based model-matrix column
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
  # an atomic vector only: anything else is left to the coercion below
  if (is.atomic(moderators) && anyNA(moderators)) {
    stop("'", argument, "' contains missing values")
  }
  if (is.character(moderators)) {
    columnNames <- colnames(data@x)
    if (is.null(columnNames)) {
      stop("'", argument, "' given by name but the design has no column names")
    }
    found <- match(moderators, columnNames)
    if (anyNA(found)) {
      stop(
        "'",
        argument,
        "' name not found in the design's column names: ",
        paste0("'", unique(moderators[is.na(found)]), "'", collapse = ", ")
      )
    }
    moderators <- found
  } else {
    moderators <- coerceOrError(moderators, "integer", name = argument)
    if (any(moderators < 1L | moderators > ncol(data@x))) {
      stop("'", argument, "' column index out of range")
    }
  }
  sort(unique(as.integer(moderators)))
}

## A term label or a column name as a model frame's column is named.
stripBackticks <- function(labels) {
  sub("^`(.*)`$", "\\1", labels)
}

## A name in forest()'s first argument that holds a formula: the terms are
## written in place, or the names handed over.
refuseHeldFormula <- function(expr) {
  stop(
    "forest()'s first argument, '",
    paste(deparse(expr, width.cutoff = 500L), collapse = " "),
    "', holds a formula; write its terms in place, as forest(x1 + x2), or ",
    "give the predictors' names as a character vector",
    call. = FALSE
  )
}

## The terms forest()'s first argument states in a 'forests' list, read as
## the right-hand side of a model formula over the fit's predictors: '.' is
## every predictor, '-' removes one. '.' is written out as the fit's own
## terms before R reads the formula, so that a removal names a term as the
## fit's formula does, log(x1) among them; a removal that names no predictor
## is refused, where a model formula would pass it by. `termCode` is the
## fit's terms as code, `predictors` their labels, and `refuseUnknown` the
## caller's refusal of a label that is none of them.
readSelectionTerms <- function(expr, termCode, predictors, refuseUnknown) {
  shown <- paste(deparse(expr, width.cutoff = 500L), collapse = " ")
  if (writesInterceptTerm(expr)) {
    stop(
      "'forest(",
      shown,
      ")': an intercept term (1, 0 or - 1) is the fit's and no forest's; ",
      "write it beside the forests",
      call. = FALSE
    )
  }
  everyPredictor <- if (length(termCode) > 0L) {
    call("(", Reduce(function(left, right) call("+", left, right), termCode))
  }
  # '.' stands for the predictors where a model formula expands it: at the
  # top of the '+' and '-' chain
  writeOut <- function(e) {
    if (identical(e, as.name(".")) && !is.null(everyPredictor)) {
      return(everyPredictor)
    }
    if (
      is.call(e) &&
        is.name(e[[1L]]) &&
        as.character(e[[1L]]) %in% c("+", "-", "(")
    ) {
      for (i in seq_along(e)[-1L]) {
        e[[i]] <- writeOut(e[[i]])
      }
    }
    e
  }
  labelsOf <- function(e) {
    rhs <- call("~", e)
    class(rhs) <- "formula"
    environment(rhs) <- baseenv()
    read <- tryCatch(stats::terms(rhs), error = function(error) {
      stop(
        "forest()'s first argument, '",
        shown,
        "': ",
        conditionMessage(error),
        call. = FALSE
      )
    })
    if (length(attr(read, "offset")) > 0L) {
      stop(
        "forest()'s first argument, '",
        shown,
        "': an offset() is the fit's and no forest's",
        call. = FALSE
      )
    }
    stripBackticks(attr(read, "term.labels"))
  }
  removed <- function(e) {
    if (isBinaryCall(e, "-")) {
      return(c(removed(e[[2L]]), labelsOf(e[[3L]])))
    }
    if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 2L) {
      return(labelsOf(e[[2L]]))
    }
    if (isBinaryCall(e, "+")) {
      return(c(removed(e[[2L]]), removed(e[[3L]])))
    }
    if (is.call(e) && identical(e[[1L]], as.name("(")) && length(e) == 2L) {
      return(removed(e[[2L]]))
    }
    character(0L)
  }
  written <- writeOut(expr)
  for (label in removed(written)) {
    if (label %not_in% predictors) {
      refuseUnknown(label)
    }
  }
  labelsOf(written)
}

## A forest's predictors, as forest()'s first argument states them, resolved
## to sorted 1-based design columns, or NULL for an unrestricted forest.
##
## Code is the terms of a model formula over the fit's predictors when any
## name in it is one of them, so that a predictor hides a variable of the
## caller's with its name (readSelectionTerms). A formula's forest() arrives
## with its labels already read against the data, and with the names its
## value gave that are columns of the design alone (readForestTerms). Any
## other code is the value it had where and when forest() was called, which
## forest() kept: nothing is looked up again here, so what the caller's
## variables hold by now changes nothing. Code that could not be evaluated
## at the call and names no predictor is refused. A value is names or
## positions; a repeat in it is refused, since a vector with a value for
## every row, a multiplier written where the predictors go, would otherwise
## be read as a few columns. `allIsNull` makes a selection of every column
## NULL, the unrestricted forest it is.
resolveForestVars <- function(vars, data, allIsNull = FALSE) {
  if (is.null(vars)) {
    return(NULL)
  }
  columnNames <- colnames(data@x)
  termLabels <- stripBackticks(attr(data@x, "term.labels"))
  predictors <- if (length(termLabels) > 0L) termLabels else columnNames
  refuseUnknown <- function(name) {
    stop(
      "'",
      name,
      "' is not a predictor of this fit (",
      paste(predictors, collapse = ", "),
      "); here a forest()'s first argument selects among them",
      call. = FALSE
    )
  }

  labels <- NULL
  named <- NULL
  if (inherits(vars, "dbartsForestTerms")) {
    expr <- vars$expr
    shown <- paste(deparse(expr, width.cutoff = 500L), collapse = " ")
    labels <- vars$labels
    named <- vars$columns
    if (is.null(labels)) {
      # the fit's terms as code; a term such as log(x1) is a predictor whose
      # name is not its variable's
      termCode <- if (length(termLabels) > 0L) {
        lapply(attr(data@x, "term.labels"), function(label) {
          tryCatch(str2lang(label), error = function(e) {
            as.name(stripBackticks(label))
          })
        })
      } else {
        lapply(columnNames, as.name)
      }
      known <- unique(c(
        ".",
        predictors,
        columnNames,
        unlist(lapply(termCode, all.vars))
      ))
      if (any(all.vars(expr) %in% known)) {
        labels <- readSelectionTerms(expr, termCode, predictors, refuseUnknown)
      } else if (isTRUE(vars$evaluated)) {
        vars <- vars$value
        if (inherits(vars, "formula")) {
          refuseHeldFormula(expr)
        }
        if (is.function(vars)) {
          refuseUnknown(shown)
        }
        if (is.null(vars)) {
          return(NULL)
        }
      } else if (length(vars$unbound) > 0L) {
        refuseUnknown(vars$unbound[1L])
      } else {
        stop(
          "forest()'s first argument, '",
          shown,
          "': ",
          if (is.null(vars$error)) "it is no selection" else vars$error,
          call. = FALSE
        )
      }
    }
  }

  if (!is.null(labels)) {
    columns <- integer(0L)
    for (label in c(stripBackticks(labels), named)) {
      found <- resolveTermColumns(label, columnNames, termLabels)
      if (is.null(found)) {
        refuseUnknown(label)
      }
      columns <- c(columns, found)
    }
    columns <- sort(unique(columns))
    if (length(columns) == 0L) {
      stop(
        "forest()'s first argument, '",
        shown,
        "', leaves the forest no predictor to split on",
        call. = FALSE
      )
    }
  } else {
    if (is.atomic(vars) && !anyNA(vars) && anyDuplicated(vars) > 0L) {
      repeated <- vars[[anyDuplicated(vars)]]
      if (!is.character(repeated)) {
        position <- suppressWarnings(as.integer(repeated))
        if (
          !is.na(position) &&
            position >= 1L &&
            position <= length(columnNames)
        ) {
          repeated <- columnNames[position]
        }
      }
      stop(
        "forest()'s first argument selects predictors of the fit and names '",
        format(repeated),
        "' more than once; name each once. A multiplier is given as 'basis ='",
        call. = FALSE
      )
    }
    columns <- resolveModerators(vars, data, "vars")
  }
  if (allIsNull && length(columns) == ncol(data@x)) {
    return(NULL)
  }
  columns
}

## A forest restricted to `columns` draws each split variable among them by
## their relative split probabilities, which a vector giving none of them a
## positive probability does not state: the engine would split on the first
## column available. Refused where a model reaches a restricted single forest,
## at creation and in setModel. `columns` NULL is no restriction, and an empty
## `splitProbabilities` is the uniform default.
refuseNoSplittableColumn <- function(splitProbabilities, columns) {
  if (
    !is.null(columns) &&
      length(splitProbabilities) > 0L &&
      !isTRUE(any(splitProbabilities[columns] > 0))
  ) {
    stop(
      "'split.probs' gives no positive probability to any column the ",
      "forest's 'vars' allows; give one of them a positive probability"
    )
  }
  invisible(NULL)
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
## 'atPrediction' drops the two rules that are about CONDITIONING DATA rather
## than about the expansion - a factor needing two levels, and a numeric column
## needing a nonzero entry - because at prediction a constant arm is the whole
## point: everyone under z = 1 gives an all-zero control column, and the width
## it must match is the fit's, checked against it there.
expandForestBasis <- function(
  basis,
  atPrediction = FALSE,
  allowEmptyLevels = atPrediction
) {
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
    if (!atPrediction && nlevels(basis) < 2L) {
      stop("a 'basis' factor must have at least two levels")
    }
    codes <- as.integer(basis)
    # a level no row takes expands to an all-zero column, refused below for a
    # numeric basis for the same reason; a mutation may keep a declared level
    # the data currently leave empty, so the amplitudes keep their columns
    empty <- tabulate(codes, nlevels(basis)) == 0L
    if (!allowEmptyLevels && any(empty)) {
      stop(
        "a 'basis' factor level with no observations contributes nothing to ",
        "a forest: '",
        levels(basis)[which(empty)[1L]],
        "'; drop it with droplevels()"
      )
    }
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
  # the calibration divides by the median nonzero row norm, so a norm that
  # overflows, or underflows to zero on a nonzero row, poisons the sampler
  rowNorms <- sqrt(rowSums(basis * basis))
  if (any(!is.finite(rowNorms) | (rowNorms == 0 & rowSums(basis != 0) > 0L))) {
    stop(
      "a 'basis' row's norm is not representable; rescale the basis to ",
      "moderate values"
    )
  }
  # an all-zero column has no observation for its amplitude to multiply, so it
  # is not a degenerate prior but a missing predictor wearing one
  if (!atPrediction && any(colSums(basis != 0) == 0L)) {
    stop("a 'basis' column of all zeros contributes nothing to a forest")
  }
  basis
}

## Re-evaluate one forest's stored basis TERM at new rows: the one-sided
## formula the forest was declared with, replayed through the same model frame
## the fit built it from, with the fit-time factor levels imposed so a level
## set that differs in ORDER cannot silently misalign amplitude j with a
## different level, and the fit-time levels of a CATEGORICAL basis re-imposed
## on the value itself (an expression such as ~ factor(z) derives its levels
## from the data it sees, so newdata alone would set the width). model.frame
## resolves an absent variable in the formula's own scope, which for a
## predicted row is a silent wrong answer rather than a missing predictor, so
## the variables are named up front the way validateXTest names its own. The
## value is built by the stored call when there is one, which carries the
## training rows' centre, scale and knots; a term stored without one is
## evaluated on the new rows.
replayForestBasis <- function(term, newdata, index) {
  vars <- all.vars(term$formula[[2L]])
  missingVars <- vars[vars %not_in% names(newdata)]
  if (length(missingVars) > 0L) {
    stop(
      "'newdata' is missing ",
      if (length(missingVars) > 1L) "variables" else "variable",
      " '",
      toString(missingVars),
      "', required by forest ",
      index,
      "'s basis (",
      deparse(term$formula),
      "); supply ",
      if (length(missingVars) > 1L) "them" else "it",
      ", or give that basis at the new rows with 'bases ='"
    )
  }
  basisFormula <- stats::reformulate(vars)
  environment(basisFormula) <- environment(term$formula)
  frame <- stats::model.frame(
    formula = basisFormula,
    data = newdata,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE,
    xlev = term$xlev
  )
  value <- if (is.null(term$predcall)) {
    evaluateForestBasis(term$formula, frame)
  } else {
    eval(term$predcall, frame, environment(term$formula))
  }
  if (!is.null(term$levels)) {
    value <- factor(as.character(value), levels = term$levels)
  }
  expandForestBasis(value, atPrediction = TRUE)
}

## The full (pre-'subset') row count of a formula fit's data, and the exact
## rows 'subset' keeps in it: the length of any one variable the formula
## names ('.' names no single column, so it is skipped for a real one), and
## that same expression read as an ordinary vector subscript into the full
## row sequence - the reading stats::model.frame() itself gives 'subset'.
## NULL when there is no rule to apply: not a formula, or no 'subset' given.
resolveFormulaBasisSubset <- function(formula, data, subsetExpr) {
  if (!is.formula(formula) || is.null(subsetExpr)) {
    return(NULL)
  }
  vars <- setdiff(all.vars(formula), ".")
  if (length(vars) == 0L) {
    return(NULL)
  }
  env <- environment(formula)
  hasData <- is.data.frame(data) || is.list(data) || is.environment(data)
  evalHere <- if (hasData) {
    function(expr) eval(expr, data, env)
  } else {
    function(expr) eval(expr, env)
  }
  full <- NROW(evalHere(as.name(vars[1L])))
  list(
    full = full,
    index = seq_len(full)[evalHere(subsetExpr)],
    kept = "'subset'"
  )
}

## Align one forest() declaration's evaluated basis to 'subsetRows'
## (resolveFormulaBasisSubset's result). A NULL basis, or a NULL 'subsetRows'
## (no rule to apply - the x/y interface, or no 'subset'), passes through
## unchanged: the basis was already validated against the un-subset data at
## the count dbartsData() checks it against. A basis at the FULL data's row
## count is restricted to the same rows the model frame keeps, the alignment
## 'weights' and every predictor column already get. A basis at the SUBSET's
## row count instead - matching the kept-row count but not the full data's -
## is an ambiguous shape, refused by name, naming the forest and both
## counts, rather than guessed at.
alignForestBasisToSubset <- function(basis, forestIndex, subsetRows) {
  if (is.null(basis) || is.null(subsetRows)) {
    return(basis)
  }
  n <- NROW(basis)
  if (n == subsetRows$full) {
    return(basis[subsetRows$index, , drop = FALSE])
  }
  if (n == length(subsetRows$index)) {
    stop(
      "forest ",
      forestIndex,
      "'s 'basis' has ",
      n,
      " rows, matching ",
      if (is.null(subsetRows$kept)) "'subset'" else subsetRows$kept,
      " (",
      length(subsetRows$index),
      ") but not the full data (",
      subsetRows$full,
      " rows); a 'basis' must cover the full data and is ",
      "subset with it"
    )
  }
  basis
}

## Validate the knobs one forest() declares. Only a declared (non-NULL) knob is
## checked; an omitted one keeps its NULL, so the caller can tell "not
## declared" from "declared at the default" - which is what makes the
## top-level-versus-forest-0 ambiguity below detectable.
validateForestKnobs <- function(spec) {
  if (!is.null(spec$n.trees)) {
    n.trees <- spec$n.trees
    numTrees <- coerceOrError(n.trees, "integer")
    if (length(numTrees) != 1L || is.na(numTrees) || numTrees < 1L) {
      stop("forest 'n.trees' must be a single integer >= 1")
    }
    spec$n.trees <- numTrees
  }
  for (name in c("power", "amplitude.prior.variance")) {
    if (!is.null(spec[[name]])) {
      value <- suppressWarnings(as.double(spec[[name]]))
      if (length(value) != 1L || is.na(value) || value <= 0.0) {
        stop("forest '", name, "' must be a single positive number")
      }
      spec[[name]] <- value
    }
  }
  if (!is.null(spec$sd)) {
    spec$sd <- validateForestSd(spec$sd)
  }
  if (!is.null(spec$base)) {
    base <- suppressWarnings(as.double(spec$base))
    if (length(base) != 1L || is.na(base) || base <= 0.0 || base >= 1.0) {
      stop("forest 'base' must be a single number in (0, 1)")
    }
    spec$base <- base
  }
  if (!is.null(spec$amplitude)) {
    spec$amplitude <- validateForestAmplitude(spec$amplitude)
  }
  spec
}

## The refusal of a held coefficient on a basis of one numeric column, which
## the engine holds at zero, so that the forest would drop out of the model.
refuseHeldOneColumn <- function(index) {
  stop(
    "forest ",
    index,
    ": amplitude = fixed() on a basis of one numeric column is not ",
    "supported yet; it would hold the forest at zero. Let the coefficient be ",
    "drawn; a column of two values can be held if it is written as a factor",
    call. = FALSE
  )
}

## A forest's 'amplitude': NULL, which draws its coefficient, or fixed(), which
## holds it at the value its forest's shape gives it. Returns "fixed" for the
## hold. A bare fixed is fixed(), as a constructor given where a value is.
validateForestAmplitude <- function(amplitude) {
  if (identical(amplitude, fixed)) {
    amplitude <- fixed()
  }
  if (!is(amplitude, "dbartsFixedPrior")) {
    stop(
      "a forest's 'amplitude' must be fixed(), which holds its coefficient, ",
      "or left out, which draws it"
    )
  }
  value <- amplitude@value
  if (!is.numeric(value) || length(value) != 1L || value != 1) {
    stop(
      "'amplitude = fixed(",
      paste(deparse(value), collapse = ""),
      ")': a held coefficient is 1 for a forest with no basis, and 0 for the ",
      "first level of a factor and 1 for the others; fixed() takes no other ",
      "value here. Write fixed(), and state the forest's size with 'sd'"
    )
  }
  "fixed"
}

## The kind of value a stated 'sd' is when a number is not read from it, for
## the message, and NULL for a number: nothing is coerced into one.
sdKindRefused <- function(sd) {
  if (is.character(sd)) {
    "a string"
  } else if (is.factor(sd)) {
    "a factor"
  } else if (inherits(sd, "Date")) {
    "a Date"
  } else if (inherits(sd, c("POSIXct", "POSIXlt"))) {
    "a date-time"
  } else if (inherits(sd, "difftime")) {
    "a time difference"
  } else if (is.list(sd)) {
    "a list"
  } else if (is.matrix(sd) || is.array(sd)) {
    "a matrix"
  } else if (is.logical(sd)) {
    "a logical"
  } else if (!is.numeric(sd)) {
    kind <- class(sd)[1L]
    paste0(if (grepl("^[aeiouAEIOU]", kind)) "an " else "a ", kind)
  }
}

## A forest's 'sd', at creation and on $setLeafPrior: one unnamed number, finite
## and positive. Infinity states no prior the map can scale, so it is refused
## rather than carried into a leaf scale or a half-Cauchy median. A number is
## a numeric or integer, whatever class it carries.
validateForestSd <- function(sd) {
  if (isSingleNA(sd) || (is.numeric(sd) && length(sd) == 1L && is.nan(sd))) {
    stop("forest 'sd' must not be NA; leave it out for the default")
  }
  if (is(sd, "dbartsSdHyperprior")) {
    stop(
      "forest 'sd' must be a number, not invchi(): a law on a forest's sd is ",
      "not supported yet"
    )
  }
  kind <- sdKindRefused(sd)
  if (!is.null(kind)) {
    stop("forest 'sd' must be a number, not ", kind)
  }
  if (length(sd) != 1L) {
    stop(
      "forest 'sd' must be a single number, not a vector of length ",
      length(sd),
      ": a forest states one sd, for every column of its basis; to size the ",
      "columns differently, rescale them in 'basis', as I(dose / 30)"
    )
  }
  if (!is.null(names(sd))) {
    stop(
      "forest 'sd' must not be named (\"",
      names(sd),
      "\"): it is one number, for every column of a basis; drop the name ",
      "with unname()"
    )
  }
  if (is.na(sd) || !is.finite(sd) || sd <= 0) {
    stop("forest 'sd' must be positive and finite")
  }
  as.double(sd)
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
    stop(
      "'forests' must be a list of forest() specifications; see ?dbartsForests"
    )
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
    if (!is.null(first$sd) && !any(hasBasis)) {
      stop(
        "this model has one forest, so its size is the fitting function's ",
        "leaf.prior = normal(sd = ), not forest(sd = )"
      )
    }
    if (!is.null(first$amplitude) && !any(hasBasis)) {
      stop(
        "this model has one forest, which has no coefficient to hold; ",
        "'amplitude' needs a model of several forests"
      )
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

## The eight doubles attr(control, "bartcore.forests")$params carries FOR EACH
## FOREST, in the order the C bridge reads them: the forest's tree count and
## structure prior, the leaf-scale factor and divisor the calibration map
## reads, the amplitude prior's variance and half-Cauchy scale, and the
## amplitude update flag. Forest 1 takes its tree count and structure prior
## from the fit's own control/tree.prior instead, so its first three are
## carried but unread.
##
## Which of the two magnitude channels a forest's `sd` reaches is decided by
## whether it carries a BASIS. A forest WITHOUT one has a plain scalar
## amplitude under a half-Cauchy scale mixture, so `sd` is that mixture's
## median and the leaf scale stays at the calibration map's anchor. A forest
## WITH one has a fixed-variance amplitude block, so `sd` multiplies the node
## scale, divided through the half-normal median 0.674. The ANCHOR is the
## family's own latent scale: sd(y) under gaussian, 1 under probit and
## pi/sqrt(3) under logistic, per unit of basis row norm.
##
## The fixed-variance channel's leaf scale factor is K-AWARE, sqrt(2/K): that
## keeps the prior on the combined location invariant to how the caller
## decomposed the mean across forests, the identity at K = 2.
forestParams <- function(specs, hasBasis, family) {
  declared <- function(value, default) {
    if (is.null(value)) default else value
  }
  leafScaleDefault <- sqrt(2 / length(specs))
  amplitudeScaleDefault <- defaultAmplitudePriorScale(family)
  lapply(seq_along(specs), function(index) {
    spec <- specs[[index]]
    withBasis <- hasBasis[index]
    as.double(c(
      declared(spec$n.trees, 50L),
      declared(spec$base, 0.25),
      declared(spec$power, 3),
      if (withBasis) declared(spec$sd, leafScaleDefault) else 1,
      if (withBasis) 0.674 else 1,
      if (withBasis) declared(spec$amplitude.prior.variance, 0.5) else 1,
      if (withBasis) 0 else declared(spec$sd, amplitudeScaleDefault),
      if (identical(spec$amplitude, "fixed")) 0 else 1
    ))
  })
}

## Resolve an interactions() specification against the fitted model matrix into
## the engine's per-forest constraint: a max-order cap (0 = uncapped) and a
## de-duplicated 2 x k integer matrix of
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
    stop(
      "'interactions' must be an interactions() specification; see ",
      "?dbartsForests"
    )
  }
  numColumns <- ncol(data@x)
  columnNames <- colnames(data@x)
  termLabels <- attr(data@x, "term.labels")

  maxOrder <- 0L
  if (!is.null(interactions$max.order)) {
    max.order <- interactions$max.order
    order <- coerceOrError(max.order, "integer")
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
        numColumns,
        "forbid"
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
        numColumns,
        "groups"
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
## engine's per-tree block-additive constraint (variant A): each whole tree
## is confined to one declared group of predictors, so the ensemble is
## exactly f = sum_G f_G
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
    stop("'blocks' must be a blocks() specification; see ?dbartsForests")
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
      numColumns,
      "groups"
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
    trees.per.group <- treesPerGroup
    counts <- coerceOrError(trees.per.group, "integer")
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

linear <- function(columns, k = NULL, sd = NULL) {
  if (missing(columns)) {
    stop("linear leaf prior requires 'columns' naming the leaf covariates")
  }
  if (!is.character(columns) && !is.numeric(columns)) {
    stop("linear leaf prior 'columns' must be a character or numeric vector")
  }
  # reuses normal()'s k and sd validation and coercions
  normalPrior <- normal(k, sd)
  new(
    "dbartsLinearPrior",
    k = normalPrior@k,
    columns = columns,
    prior.sd = normalPrior@prior.sd
  )
}

gp <- function(
  columns,
  k = NULL,
  lengthscale = NULL,
  max.leaf.size = 256L,
  sd = NULL
) {
  if (missing(columns)) {
    stop("gp leaf prior requires 'columns' naming the leaf covariates")
  }
  if (!is.character(columns) && !is.numeric(columns)) {
    stop("gp leaf prior 'columns' must be a character or numeric vector")
  }
  if (
    !is.null(lengthscale) &&
      (!is.numeric(lengthscale) ||
        length(lengthscale) == 0L ||
        anyNA(lengthscale) ||
        any(lengthscale <= 0))
  ) {
    stop("gp leaf prior 'lengthscale' must be positive")
  }
  max.leaf.size <- coerceOrError(max.leaf.size, "integer")
  if (
    length(max.leaf.size) != 1L || is.na(max.leaf.size) || max.leaf.size < 1L
  ) {
    stop("gp leaf prior 'max.leaf.size' must be a positive integer")
  }
  # reuses normal()'s k and sd validation and coercions
  normalPrior <- normal(k, sd)
  new(
    "dbartsGPPrior",
    k = normalPrior@k,
    columns = columns,
    lengthscale = if (is.null(lengthscale)) NULL else as.double(lengthscale),
    max.leaf.size = max.leaf.size,
    prior.sd = normalPrior@prior.sd
  )
}

## A named spread, wherever it is spelled: NULL leaves it unnamed, anything
## else must be a single positive finite number. NA is a missing value, not
## the unnamed spelling, and NaN carries no intent and cannot serve as a
## divisor, so both are refused here rather than surviving to the bridge's own
## last-line check. The result is NA_real_ for unnamed.
## A specification built outside the constructors, such as one read back from
## $getLeafPrior, is held to their rules.
refuseInvalidLeafPrior <- function(leaf.prior) {
  valid <- methods::validObject(leaf.prior, test = TRUE)
  if (!isTRUE(valid)) {
    stop(valid, call. = FALSE)
  }
  invisible(NULL)
}

validateNamedScale <- function(value, name) {
  if (is.null(value)) {
    return(NA_real_)
  }
  if (isSingleNA(value)) {
    stop(
      "'",
      name,
      "' must not be NA: NA is a missing value; NULL leaves it unnamed",
      call. = FALSE
    )
  }
  value <- coerceOrError(value, "numeric", name)
  if (length(value) != 1L) {
    stop("'", name, "' must be a single number")
  }
  if (is.na(value) || !is.finite(value) || value <= 0.0) {
    stop("'", name, "' must be positive")
  }
  value
}

## The mid-chain setter's value: the same rules, plus a refusal of the NA that
## spells "unnamed" at creation. There is no family default to fall back on
## once a sampler exists, so an absent value is a malformed one.
validateLiveScale <- function(value, name) {
  if (is.null(value) || isSingleNA(value)) {
    stop("'", name, "' must be a positive finite number")
  }
  validateNamedScale(value, name)
}

## A leaf prior's 'sd': NULL, a positive number, or an invchi() law on it. A
## law on k and the string forms k keeps for 0.9-x are refused by name.
validateLeafSd <- function(sd) {
  if (is.null(sd) || is(sd, "dbartsSdHyperprior")) {
    return(sd)
  }
  if (!is(sd, "dbartsLeafHyperprior") && !is.character(sd) && !isSingleNA(sd)) {
    kind <- sdKindRefused(sd)
    if (!is.null(kind)) {
      stop("'sd' must be a number or invchi(), not ", kind)
    }
  }
  if (is(sd, "dbartsLeafHyperprior")) {
    stop(
      "'sd' takes a number or invchi(), a law on the sd itself; a law on k ",
      "is spelled k = chi(), and the same prior on the sd is ",
      "sd = invchi(df, k.scale / scale)"
    )
  }
  if (is.character(sd)) {
    stop(
      "'sd' must be a number or invchi(); unlike 'k' it takes no string form"
    )
  }
  validateNamedScale(sd, "sd")
}

## The engine's inputs for a leaf prior, the model's prior.scale anchor and its
## leaf hyperprior, translated from whichever of 'k' and 'sd' it names. The
## engine's k is relative to its anchor, and only their ratio enters a draw, so
## a named sd rides a reference k of 2: a fixed sd x is anchor 2x with k fixed
## at 2, and invchi(df, c) is anchor 2c with k ~ chi(df, 2). The bridge starts a
## drawn k at 2, so the chain starts at the named spread, and the binary
## default and every k spelling at the defaults keep bitwise engine inputs.
## invchi(df, 0) is the improper limit, which no anchor can state and
## chi(df, Inf) is. A multi-forest fit refuses a named sd before any of this:
## its calibration map pins every forest's scale.
resolveLeafPrior <- function(
  leaf.prior,
  binary,
  monotone = FALSE,
  multiForest = FALSE
) {
  refuseInvalidLeafPrior(leaf.prior)
  sd <- leaf.prior@prior.sd
  if (is.null(sd)) {
    return(list(
      prior.scale = NA_real_,
      leaf.hyperprior = resolveLeafHyperprior(
        leaf.prior@k,
        binary,
        monotone = monotone,
        multiForest = multiForest
      )
    ))
  }
  if (multiForest) {
    stop(
      "a multi-forest model does not support a named leaf-prior 'sd': the ",
      "leaf prior's 'sd' is not a forest's 'sd': each forest's scale is set ",
      "by forest(sd = ), which states that forest's share of the combined ",
      "location's prior (see ?forest)"
    )
  }
  if (is.numeric(sd)) {
    return(list(
      prior.scale = 2.0 * sd,
      leaf.hyperprior = newValidated("dbartsFixedHyperprior", k = 2.0)
    ))
  }
  if (monotone) {
    stop(
      "an 'sd' hyperprior is not supported under a monotone constraint; ",
      "supply a fixed numeric sd (the truncated leaf law has no chi-k update)"
    )
  }
  if (sd@scale == 0.0) {
    return(list(
      prior.scale = NA_real_,
      leaf.hyperprior = chi(sd@df, Inf)
    ))
  }
  list(
    prior.scale = 2.0 * sd@scale,
    leaf.hyperprior = chi(sd@df, 2.0)
  )
}

normal <- function(k = NULL, sd = NULL) {
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
      !is(k, "dbartsLeafHyperprior") &&
      (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0.0)
  ) {
    if (is(k, "dbartsSdHyperprior")) {
      stop(
        "'k' takes a number or chi(), a law on k; invchi() is a law on the ",
        "sd and is spelled sd = invchi()"
      )
    }
    stop("'k' must be a positive scalar or a hyperprior specification")
  }
  if (is.function(sd)) {
    sd <- sd()
  }
  sd <- validateLeafSd(sd)
  if (!is.null(k) && !is.null(sd)) {
    stop(
      "give either 'k' (relative to the data's scale) or 'sd' (on the ",
      "family's scale) to a leaf prior, not both"
    )
  }
  new("dbartsNormalPrior", k = k, prior.sd = sd)
}

chisq <- function(df = 3, quant = 0.9) {
  newValidated("dbartsChiSqPrior", df = df, quantile = quant)
}

fixed <- function(value = 1.0) {
  newValidated("dbartsFixedPrior", value = value)
}

chi <- function(df = 1.5, scale = 2.0, degreesOfFreedom) {
  if (!missing(degreesOfFreedom)) {
    if (!missing(df)) {
      stop(
        "'degreesOfFreedom' and 'df' name the same value on chi(); supply one"
      )
    }
    warnOnce(
      "tombstone.degreesOfFreedom.chi",
      "chi()'s 'degreesOfFreedom' is now 'df'; the value was used. The old ",
      "name is removed in dbarts ",
      tombstoneExpiry,
      ".",
      class = "dbartsDeprecatedWarning"
    )
    df <- degreesOfFreedom
  }
  newValidated("dbartsChiHyperprior", degreesOfFreedom = df, scale = scale)
}

invchi <- function(df = 1.5, scale) {
  if (missing(scale)) {
    stop(
      "invchi() requires 'scale': an sd on the family's scale has no ",
      "data-free default"
    )
  }
  for (value in list(df, scale)) {
    if (!is.numeric(value) || length(value) != 1L) {
      stop("invchi() 'df' and 'scale' must be single numbers")
    }
  }
  newValidated(
    "dbartsSdHyperprior",
    df = as.double(df),
    scale = as.double(scale)
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
  refuseNaN(rho, "rho")
  refuseNaN(update.delay, "update.delay")
  if (isSingleNA(rho)) {
    refuseNAForNull("rho", "dart", "the default, the number of predictors")
  }
  if (isSingleNA(update.delay)) {
    refuseNAForNull("update.delay", "dart", "the default, half the burn-in")
  }
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

## Per-forest interaction constraint, passed as interactions = to
## dbarts()/bart2() (and mu.interactions /
## tau.interactions to bcf()). Packages the raw specification; groups and forbid
## resolve against the model matrix, and every value is validated, at fit time
## in resolveInteractions. max.order caps the number of DISTINCT split variables
## on any root-to-leaf path; groups is a co-occurrence allow-list (named columns
## may share a path only with group-mates); forbid names column sets barred from
## sharing a path. Not exported, like the priors: it resolves by bare name
## inside the arguments that take it (evalInForestVocabulary), and
## dbartsForests is its exported face.
interactions <- function(max.order = NULL, groups = NULL, forbid = NULL) {
  if (is.null(max.order) && is.null(groups) && is.null(forbid)) {
    stop("interactions() needs at least one of 'max.order', 'groups', 'forbid'")
  }
  structure(
    list(max.order = max.order, groups = groups, forbid = forbid),
    class = "dbartsInteractions"
  )
}

## The monotone priors, the default first: the one place the default is set.
## It is monotone()'s prior formal, so match.arg takes its first element when
## prior is not given, as the plain-vector shorthand does.
MONOTONE_PRIORS <- c("joint", "leaf")

## Per-predictor monotone constraint and the prior it is read under, passed as
## monotone = to dbarts()/bart()/dbartsSpec(). directions is the vector the
## plain shorthand takes, named or positional, so a predictor named "prior"
## needs nothing special; it resolves against the model matrix, and every
## element is validated, at fit time in resolveMonotone. prior is matched here.
## Not exported, like interactions().
monotone <- function(directions, prior) {
  if (missing(directions) || is.null(directions)) {
    stop("monotone() requires 'directions', the per-predictor directions")
  }
  if (
    !missing(prior) &&
      (!is.character(prior) || length(prior) != 1L || is.na(prior))
  ) {
    stop(
      "monotone() 'prior' must be one of ",
      paste0("\"", MONOTONE_PRIORS, "\"", collapse = ", ")
    )
  }
  prior <- match.arg(prior)
  structure(
    list(directions = directions, prior = prior),
    class = "dbartsMonotone"
  )
}
formals(monotone)$prior <- MONOTONE_PRIORS

## Per-forest block-additive constraint (variant A), passed as blocks = to
## dbarts() / bart2() (and mu.blocks / tau.blocks to bcf()). Confines each
## WHOLE tree to one
## declared group of predictors, so the ensemble is exactly f = sum_G f_G (a
## clean functional-ANOVA / grouped-GAMI decomposition). groups is a list
## partitioning the forest's predictors into disjoint blocks (by model-matrix
## column name - a bare factor term name expands to its indicator columns - or
## index); the partition must be TOTAL, so every predictor is named exactly once
## (a predictor named in no block would be masked out of every tree and go dead).
## trees.per.group optionally fixes how many of the n.trees trees each block
## gets; NULL distributes them as evenly as possible. Everything is validated at
## fit time in resolveBlocks. Not exported, like interactions().
blocks <- function(groups, trees.per.group = NULL) {
  if (missing(groups) || is.null(groups)) {
    stop("blocks() requires 'groups', a list partitioning the predictors")
  }
  structure(
    list(groups = groups, trees.per.group = trees.per.group),
    class = "dbartsBlocks"
  )
}

## One forest of a model, written as a term of a fitting function's formula
## or inside its forests = list(forest(), forest(...)) argument. Every knob is
## per forest, so the fitting functions grow exactly one argument however many
## forests a model has. vars is the predictors the forest splits on and the
## one argument that may be given unnamed, so that any later formal is an
## addition: it is kept as code, with the place it was written, and read as
## the right-hand side of a model formula once the fit's predictors are known
## (resolveForestVars); a value handed over is names or positions. basis is
## the data the forest's amplitudes multiply, a one-sided formula or a vector,
## expanded by R's own model-matrix rule; n.trees, base and power are its
## tree-structure prior; sd is its total's prior scale in units of the
## family's latent scale (sd(y) under gaussian, 1 under probit, pi/sqrt(3)
## under logistic) per unit of basis row norm; amplitude.prior.variance is the
## N(0, .) variance of the amplitudes on its basis; amplitude = fixed() holds
## those amplitudes at the value their shape gives them, and left out draws
## them; interactions and blocks are this forest's own constraints - the
## arguments of the same names on the fitting function are the FIRST forest's.
## Every knob defaults to NULL, "not declared", which is what lets a
## declaration that collides with one of those arguments refuse rather than
## silently win. Validated at fit time, in resolveForests. Not exported, like
## interactions() and blocks(); a formula's forest() term is recognized by
## name.
forest <- function(
  vars = NULL,
  basis = NULL,
  sd = NULL,
  n.trees = NULL,
  base = NULL,
  power = NULL,
  amplitude = NULL,
  interactions = NULL,
  blocks = NULL,
  amplitude.prior.variance = NULL
) {
  # the arguments as the caller gave them, a forwarded '...' spelled out
  given <- match.call(function(...) NULL, sys.call(), expand.dots = FALSE)$...
  if (numUnnamed(given) > 1L) {
    refuseSecondUnnamed()
  }
  written <- captureForestVars(substitute(vars), callingPlace(parent.frame()))
  if (inherits(written, "dbartsForestTerms")) {
    # the argument's value here and now, taken once: whatever the caller's
    # variables hold later, a forest built in a loop or by lapply() keeps the
    # selection it was given. Code that cannot be evaluated here, terms over
    # predictors among it, is kept as code alone
    taken <- tryCatch(
      withCallingHandlers(list(vars), warning = function(w) {
        invokeRestart("muffleWarning")
      }),
      error = function(e) e
    )
    if (inherits(taken, "error")) {
      symbols <- all.vars(written$expr)
      written$unbound <- symbols[
        !vapply(symbols, exists, NA, envir = written$env)
      ]
      written$error <- conditionMessage(taken)
    } else {
      written["value"] <- taken
      written$evaluated <- TRUE
    }
  }
  structure(
    list(
      basis = basis,
      vars = written,
      n.trees = n.trees,
      base = base,
      power = power,
      sd = sd,
      interactions = interactions,
      blocks = blocks,
      amplitude.prior.variance = amplitude.prior.variance,
      amplitude = amplitude
    ),
    class = "dbartsForest"
  )
}

## How many of a call's arguments are given without a name.
numUnnamed <- function(arguments) {
  given <- names(arguments)
  if (is.null(given)) length(arguments) else sum(!nzchar(given))
}

## The one-unnamed-argument rule, which keeps every later formal of forest()
## an addition: the text for a second unnamed argument, at either door.
refuseSecondUnnamed <- function() {
  stop(
    "forest() takes one unnamed argument, the predictors the forest splits ",
    "on, joined by '+' as forest(x1 + x2); every other argument is given by ",
    "name: a multiplier is 'basis =' and a size is 'sd ='",
    call. = FALSE
  )
}

## forest()'s first argument as written. A value (names, positions, NULL) is
## kept as it is; code is kept unevaluated with the environment it was
## written in, and forest() puts beside it the value it has at the call. A
## tilde there is a multiplier written where the predictors go, and is
## refused rather than read as predictors.
captureForestVars <- function(expr, env) {
  if (!is.language(expr)) {
    return(expr)
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("~"))) {
    stop(
      "forest()'s first argument is the predictors the forest splits on, ",
      "written without '~', as forest(x1 + x2); a multiplier is 'basis ='",
      call. = FALSE
    )
  }
  structure(list(expr = expr, env = env), class = "dbartsForestTerms")
}

## The heteroscedastic variance forest's own specification, passed as the
## SAME variance = argument of dbarts()/dbartsSpec()/bart2() that already
## takes the plain selector (NULL/FALSE for none, TRUE/a one-sided formula/
## character or index vector for the column subset the variance forest
## reads). vars resolves through the identical resolveVarianceColumns the
## shorthand uses - one selector vocabulary - with vars = NULL on THIS
## object meaning every column, distinct from variance = NULL's "no
## variance forest". n.trees/base/power default to NULL, "not declared",
## matching forest(); resolveSamplerSpec falls each back to a default when
## NULL (40 trees; the mean forest's tree.prior base/power). Validated
## at fit time. Not exported, like forest().
varianceForest <- function(
  vars = NULL,
  n.trees = NULL,
  base = NULL,
  power = NULL
) {
  structure(
    list(vars = vars, n.trees = n.trees, base = base, power = power),
    class = "dbartsVarianceForest"
  )
}

format.dbartsVarianceForest <- function(x, ...) {
  describe <- function(value, default) {
    if (is.null(value)) default else paste(deparse(value), collapse = " ")
  }
  c(
    paste0("vars    = ", describe(x$vars, "<all columns>")),
    paste0("n.trees = ", describe(x$n.trees, "<default 40>")),
    paste0("base    = ", describe(x$base, "<mean forest's>")),
    paste0("power   = ", describe(x$power, "<mean forest's>"))
  )
}

print.dbartsVarianceForest <- function(x, ...) {
  cat("dbarts variance forest specification\n")
  cat(paste0("  ", format(x)), sep = "\n")
  invisible(x)
}

## The exported face of the prior constructors: one object, so that no
## generic name (normal, chisq, fixed, chi) enters the search path to be
## masked by or to mask another package by attach order. Inside the
## tree.prior and leaf.prior arguments of the fitting functions - and inside
## the 'sigma' argument of the families that draw a residual scale - the
## same constructors are available by bare name.
dbartsPriors <- list(
  cgm = cgm,
  dart = dart,
  normal = normal,
  linear = linear,
  gp = gp,
  chisq = chisq,
  fixed = fixed,
  chi = chi,
  invchi = invchi
)

## The exported face of the forest constructors, for the reason dbartsPriors
## exists: 'forest' and 'blocks' are names other packages attach. Inside the
## arguments that take them they resolve by bare name (resolveForestArguments).
dbartsForests <- list(
  interactions = interactions,
  blocks = blocks,
  monotone = monotone,
  forest = forest,
  varianceForest = varianceForest
)

## What the arguments that take a forest constructor resolve by bare name:
## dbartsForests and the one prior constructor a forest's 'amplitude' takes,
## which stays under dbartsPriors for the caller who names it outside.
forestConstructors <- c(dbartsForests, dbartsPriors["fixed"])

## Each door argument taking a forest constructor, and the vocabulary it
## resolves over. In the order the doors forced them before they resolved by
## name, so that an error names the same argument.
FOREST_ARGUMENT_VOCABULARIES <- list(
  forests = c("forest", "interactions", "blocks", "fixed"),
  interactions = "interactions",
  blocks = "blocks",
  monotone = "monotone",
  variance = "varianceForest"
)

## The door arguments that take a forest constructor, resolved from the
## caller's own unevaluated arguments in 'matchedCall'. An absent argument is
## NULL, every door's default.
resolveForestArguments <- function(
  matchedCall,
  evalEnv,
  arguments = names(FOREST_ARGUMENT_VOCABULARIES)
) {
  resolved <- list()
  for (name in arguments) {
    expr <- matchedCall[[name]]
    resolved[name] <- list(
      if (is.null(expr)) {
        NULL
      } else {
        evalInForestVocabulary(
          expr,
          forestConstructors[FOREST_ARGUMENT_VOCABULARIES[[name]]],
          evalEnv
        )
      }
    )
  }
  resolved
}
