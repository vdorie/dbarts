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
  power = 2,
  base = 0.95,
  split.probs = NULL,
  leaf.prior,
  resid.prior,
  monotone = NULL,
  multiForest = FALSE,
  kHyperprior = control@binary,
  parentEnv,
  treePrior = NULL
) {
  matchedCall <- match.call()

  # the prior vocabulary shadows the caller's environment inside these
  # arguments only: bare names like normal(chi(1.5)) resolve here no matter
  # what packages are attached, and nothing is exported under generic names
  vocabulary <- c(
    dbartsPriors,
    list(num.vars = ncol(data@x), numvars = ncol(data@x))
  )
  # a bare constructor name (leaf.prior = normal) means its defaults; a value
  # that is already a prior object passes through
  resolveSpec <- function(expr, name, class, label) {
    evalInVocabulary(
      expr,
      vocabulary,
      parentEnv,
      resolvedAs(name, class, paste(label, "specification"))
    )
  }

  tree.prior <- resolveTreePrior(
    matchedCall,
    ncol(data@x),
    parentEnv,
    treePrior
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

## The vocabulary a split.probs expression is written in, over the caller's
## frame: dart() and the number of predictors, under its current name num.vars
## and its backward-compatible alias numvars, so a bare 1 / num.vars or
## 1 / numvars resolves either way.
splitProbsVocabulary <- function(numVars) {
  list(dart = dart, num.vars = numVars, numvars = numVars)
}

## Which of the constructors cgm() and dart() a call's head names, in any
## spelling a caller writes: bare, from dbartsPriors (qualified or not), or
## from the package's namespace. NULL for any other expression.
treePriorConstructorHead <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }
  head <- expr[[1L]]
  priorsList <- list(
    quote(dbartsPriors),
    quote(dbarts::dbartsPriors),
    call(":::", quote(dbarts), quote(dbartsPriors))
  )
  name <- if (is.symbol(head)) {
    as.character(head)
  } else if (
    is.call(head) &&
      identical(head[[1L]], as.name("$")) &&
      any(vapply(priorsList, identical, NA, head[[2L]])) &&
      is.symbol(head[[3L]])
  ) {
    as.character(head[[3L]])
  } else if (
    is.call(head) &&
      length(head) == 3L &&
      (identical(head[[1L]], as.name("::")) ||
        identical(head[[1L]], as.name(":::"))) &&
      identical(head[[2L]], as.name("dbarts"))
  ) {
    as.character(head[[3L]])
  }
  if (length(name) == 1L && name %in% c("cgm", "dart")) name
}

refuseCgmAsSplitProbs <- function() {
  stop(
    "'split.probs' takes a vector of split probabilities or dart(); cgm() is ",
    "the retired tree.prior's value: write power, base and split.probs",
    call. = FALSE
  )
}

## A fitting function's 'power' or 'base', as `matchedCall` names it and
## `parentEnv` evaluates it: the default where it is not given or is NULL,
## and otherwise a single positive number, a base below 1, kept as given.
plainTreePriorValue <- function(matchedCall, name, parentEnv, default) {
  expr <- matchedCall[[name]]
  value <- if (!is.null(expr)) eval(expr, parentEnv)
  if (is.null(value)) {
    return(default)
  }
  valid <- is.numeric(value) &&
    length(value) == 1L &&
    !is.na(value) &&
    value > 0 &&
    (name == "power" || value < 1)
  if (!valid) {
    stop(
      "'",
      name,
      "' must be a single ",
      if (name == "power") "positive number" else "number in (0, 1)",
      call. = FALSE
    )
  }
  as.double(value)
}

## The record a split.probs expression builds, resolved in its vocabulary over
## `parentEnv`: dart() or a bare dart a dbartsDartPrior, a vector of
## probabilities or NULL a dbartsCGMPrior carrying it as its specification,
## each at the default depth.
splitProbsRecord <- function(expr, numVars, parentEnv) {
  if (is.null(expr)) {
    return(cgm())
  }
  written <- recoverForwardedArgument(expr, parentEnv)$expr
  if (identical(treePriorConstructorHead(written), "cgm")) {
    refuseCgmAsSplitProbs()
  }
  value <- evalInVocabulary(expr, splitProbsVocabulary(numVars), parentEnv)
  if (identical(value, dart)) {
    value <- dart()
  }
  if (is(value, "dbartsDartPrior")) {
    return(value)
  }
  if (is(value, "dbartsCGMPrior") || identical(value, cgm)) {
    refuseCgmAsSplitProbs()
  }
  if (length(value) > 0L && !is.numeric(value)) {
    stop("'split.probs' must be a numeric vector or dart()", call. = FALSE)
  }
  cgm(split.probs = value)
}

## The model's tree-prior record: the split probabilities' record, or
## `record`, a whole record a retired tree.prior handed over, with the
## fitting function's power and base written on it. `matchedCall` is
## parsePriors's own call, whose arguments are written in `parentEnv`.
resolveTreePrior <- function(matchedCall, numVars, parentEnv, record = NULL) {
  power <- plainTreePriorValue(matchedCall, "power", parentEnv, 2)
  base <- plainTreePriorValue(matchedCall, "base", parentEnv, 0.95)
  if (is.null(record)) {
    record <- splitProbsRecord(
      matchedCall[["split.probs"]],
      numVars,
      parentEnv
    )
  }
  record@power <- power
  record@base <- base
  record
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

## Refuses a prior object supplied together with a shorthand argument that
## would have built the same prior, naming both: leaf.prior with k, the
## retired resid.prior with sigdf and sigquant. Presence in the matched call,
## not an explicit-NULL value, is what collides.
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

## The refusal of a language object given as k, which 0.9-34's bart2 took
## (k = quote(chi) evaluated to the chi hyperprior there). Each form named
## fits: chi() written in the call, and, for arguments held in a list
## (do.call, bartCause's args.rsp, stan4bart's bart_args), the object from
## dbartsPriors or its string.
refuseQuotedK <- function(k) {
  if (is.language(k)) {
    stop(
      "'k' is a quoted expression, which is not taken as a prior; write ",
      "chi(1.25, Inf) directly in the call, or, where the arguments are held ",
      "in a list, give dbartsPriors$chi(1.25, Inf) or the string ",
      "\"chi(1.25, Inf)\"",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## Turn a normal prior's raw k into the model's leaf hyperprior: NULL is the
## family default (2 for continuous responses, chi(1.5, 2) for binary and
## nbinom, drawsLeafKByDefault; 'binary' carries that flag),
## a positive scalar is fixed, and a hyperprior object passes through. Under a
## monotone constraint k is fixed for both families (an unsupplied k resolves
## to 2, the truncated leaf distribution having no clean chi-k update) and a
## chi hyperprior is refused.
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
      "supply a fixed numeric k (the truncated leaf prior has no chi-k update)"
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

## The monotone direction vocabulary: the words, matched as base R matches a
## choice (case-sensitive, a unique abbreviation taken), and the codes as the
## numbers 1, -1 and 0 or as the strings c() makes of them when words and codes
## share a vector.
MONOTONE_DIRECTION_WORDS <- c("increasing", "decreasing")
MONOTONE_DIRECTION_CODES <- c("1" = 1L, "-1" = -1L, "0" = 0L)

## Code of a single monotone direction element, in {-1, 0, 1}.
parseMonotoneSign <- function(value) {
  direction <- NA_integer_
  if (length(value) == 1L && is.character(value)) {
    direction <- unname(MONOTONE_DIRECTION_CODES[value])
    if (is.na(direction)) {
      direction <- c(1L, -1L)[
        pmatch(value, MONOTONE_DIRECTION_WORDS, duplicates.ok = TRUE)
      ]
    }
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

## The columns each term of a design produced, one integer vector per term in
## `termLabels`' order, read from the design's own layout: its columns stand
## in its terms' order, each term's as one run, named for the term itself or,
## where it gave several (a factor's indicators, a matrix's columns), for the
## term and a suffix. A run is told from its neighbor by position, never by
## a name's prefix alone: where the NEXT term's label looks like a column of
## this one ('g' beside 'g.total', or beside a numeric predictor 'g.u'), the
## run ends where that term's own columns begin. NULL when the names are not
## such a layout, a test matrix's columns in another order among them.
termColumnBlocks <- function(columnNames, termLabels) {
  numColumns <- length(columnNames)
  numTerms <- length(termLabels)
  if (numTerms == 0L || numColumns == 0L) {
    return(NULL)
  }
  blocks <- vector("list", numTerms)
  at <- 1L
  for (j in seq_len(numTerms)) {
    term <- termLabels[j]
    if (at <= numColumns && columnNames[at] == term) {
      blocks[[j]] <- at
      at <- at + 1L
      next
    }
    prefix <- paste0(term, ".")
    end <- at
    while (end <= numColumns && startsWith(columnNames[end], prefix)) {
      end <- end + 1L
    }
    if (j < numTerms && end > at && startsWith(termLabels[j + 1L], prefix)) {
      # the next term's own columns close the run: its one column, the last
      # that bears its name, or its several, from the first that does
      following <- termLabels[j + 1L]
      run <- seq.int(at, end - 1L)
      own <- run[columnNames[run] == following]
      if (length(own) > 0L) {
        end <- own[length(own)]
      } else {
        own <- run[startsWith(columnNames[run], paste0(following, "."))]
        if (length(own) > 0L) {
          end <- own[1L]
        }
      }
    }
    blocks[[j]] <- seq_len(end - at) + (at - 1L)
    at <- end
  }
  if (at != numColumns + 1L) {
    return(NULL)
  }
  blocks
}

## Resolve one predictor selector name to its 1-based model-matrix column
## indices: a term label is every column the term produced, by the design's
## layout (termColumnBlocks), and any other name the column that bears it.
## Returns NULL when the name is neither a term nor a column, so each caller
## can raise its own diagnostic; a recognized term with no columns yields
## integer(0), not NULL. columnNames may be NULL. Names that are no design in
## its terms' order have no layout to read, and there a term's columns are
## those named for it.
##
## `asCode` is for a name a caller wrote as code: an empty one names nothing,
## and one that is no label and no column as it stands is also matched as R
## writes the same code out. `refuseShared` refuses a name that two columns
## carry, which cannot tell them apart; without it the first is taken.
resolveTermColumns <- function(
  name,
  columnNames,
  termLabels,
  asCode = FALSE,
  refuseShared = FALSE
) {
  if (asCode && !nzchar(name)) {
    return(NULL)
  }
  written <- function(names) {
    vapply(
      names,
      function(name) {
        tryCatch(
          paste(deparse(str2lang(name), width.cutoff = 500L), collapse = " "),
          error = function(e) name
        )
      },
      "",
      USE.NAMES = FALSE
    )
  }
  term <- match(name, termLabels)
  if (asCode && is.na(term) && !any(columnNames == name)) {
    term <- match(name, written(termLabels))
  }
  if (!is.na(term)) {
    label <- termLabels[term]
    exact <- which(columnNames == label)
    # one column of that name, and none that could be another of the term's
    if (
      length(exact) == 1L && !any(startsWith(columnNames, paste0(label, ".")))
    ) {
      return(exact)
    }
    blocks <- termColumnBlocks(columnNames, termLabels)
    if (!is.null(blocks)) {
      return(blocks[[match(label, termLabels)]])
    }
    if (length(exact) > 0L) {
      return(exact[1L])
    }
    return(which(startsWith(columnNames, paste0(label, "."))))
  }
  columns <- which(columnNames == name)
  if (asCode && length(columns) == 0L) {
    columns <- which(nzchar(columnNames) & written(columnNames) == name)
  }
  if (length(columns) == 0L) {
    return(NULL)
  }
  if (length(columns) > 1L) {
    if (refuseShared) {
      refuseSharedName(name, length(columns))
    }
    return(columns[1L])
  }
  columns
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
  if (inherits(variance, "dbartsForestCode")) {
    # varianceForest()'s predictors, read as forest()'s are; a switch held in
    # a variable is the switch
    if (isTRUE(variance$evaluated) && is.logical(variance$value)) {
      return(resolveVarianceColumns(variance$value, data, argument))
    }
    columns <- resolveForestVars(
      variance,
      data,
      owner = "varianceForest()'s first argument"
    )
    return(if (is.null(columns)) seq_len(numColumns) else columns)
  }
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

## A name in forest()'s first argument that holds a formula with a left-hand
## side: only the right-hand side of a one-sided formula names predictors.
refuseHeldFormula <- function(expr, owner = "forest()'s first argument") {
  stop(
    owner,
    ", '",
    paste(deparse(expr, width.cutoff = 500L), collapse = " "),
    "', holds a formula with a left-hand side; write the predictors as the ",
    "right-hand side alone, as ~ x1 + x2",
    call. = FALSE
  )
}

## A name that several columns carry selects none of them.
refuseSharedName <- function(name, count) {
  stop(
    "'",
    name,
    "' is the name of ",
    count,
    " predictors of this fit; a name selects one, so select these by ",
    "position",
    call. = FALSE
  )
}

## The columns forest()'s first argument selects in a 'forests' list when it
## is read as terms over the fit's predictors: '+' adds a term, '-' removes
## one, parentheses group, and '.' is every column of the design, in the
## order written, as a model formula reads them. A term names a predictor by
## its label as the fit holds it, the term label on a formula fit and the
## column name on a matrix fit, written as code or as a backticked name
## (resolveTermColumns); one that names none is refused, added or removed alike,
## where a model formula would pass a removal by. The arithmetic is on
## columns, not on R's reading of the code, so that a column whose name is
## itself a call or an operator is told from that call, and a column with no
## name or with another's is still one of '.'.
readSelectionTerms <- function(expr, data, refuseUnknown) {
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
  columnNames <- designColumnNames(data)
  termLabels <- stripBackticks(attr(data@x, "term.labels"))
  columnsOf <- function(e) {
    if (identical(e, as.name("."))) {
      return(seq_len(ncol(data@x)))
    }
    if (is.numeric(e)) {
      # alone a number is a position; among terms it would be a column's name
      stop(
        "'",
        format(e),
        "' is not a predictor; give positions as the whole argument, as ",
        "c(2, 3)",
        call. = FALSE
      )
    }
    if (is.call(e) && identical(e[[1L]], as.name("offset"))) {
      stop(
        "forest()'s first argument, '",
        shown,
        "': an offset() is the fit's and no forest's",
        call. = FALSE
      )
    }
    label <- if (is.name(e)) {
      as.character(e)
    } else {
      paste(deparse(e, width.cutoff = 500L), collapse = " ")
    }
    columns <- resolveTermColumns(
      label,
      columnNames,
      termLabels,
      asCode = TRUE,
      refuseShared = TRUE
    )
    if (is.null(columns)) {
      refuseUnknown(label)
    }
    columns
  }
  # the columns `e` leaves of `selected`, its terms taken left to right
  take <- function(e, selected, removing = FALSE) {
    if (isBinaryCall(e, c("+", "-"))) {
      selected <- take(e[[2L]], selected, removing)
      flips <- identical(e[[1L]], as.name("-"))
      return(take(e[[3L]], selected, xor(removing, flips)))
    }
    if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 2L) {
      return(take(e[[2L]], selected, !removing))
    }
    if (is.call(e) && identical(e[[1L]], as.name("(")) && length(e) == 2L) {
      # a group is read on its own and then added or removed whole
      columns <- take(e[[2L]], integer(0L))
    } else {
      columns <- columnsOf(e)
    }
    if (removing) setdiff(selected, columns) else union(selected, columns)
  }
  sort(take(expr, integer(0L)))
}

## A design's column names, one for each column: a column with no name, an
## NA or an empty one, is "", which no name written by a caller matches.
designColumnNames <- function(data) {
  columnNames <- colnames(data@x)
  if (is.null(columnNames)) {
    return(character(ncol(data@x)))
  }
  columnNames[is.na(columnNames)] <- ""
  columnNames
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
## variables hold by now changes nothing, and what that evaluation warned of
## is raised now, the value being used. Code that could not be evaluated at
## the call and names no predictor is refused. A value is names or positions:
## a logical and a factor are refused, not read through their codes; a repeat
## is refused, since a vector with a value for every row, a multiplier
## written where the predictors go, would otherwise be read as a few columns;
## and so is a name that two columns carry. `allIsNull` makes a selection of
## every column NULL, the unrestricted forest it is.
resolveForestVars <- function(
  vars,
  data,
  allIsNull = FALSE,
  owner = "forest()'s first argument"
) {
  if (is.null(vars)) {
    return(NULL)
  }
  columnNames <- if (is.null(colnames(data@x))) {
    NULL
  } else {
    designColumnNames(data)
  }
  termLabels <- stripBackticks(attr(data@x, "term.labels"))
  predictors <- if (length(termLabels) > 0L) termLabels else columnNames
  refuseUnknown <- function(name) {
    stop(
      "'",
      name,
      "' is not a predictor of this fit (",
      paste(predictors, collapse = ", "),
      "); here a ",
      owner,
      " selects among them",
      call. = FALSE
    )
  }

  columns <- NULL
  if (inherits(vars, "dbartsForestTerms")) {
    expr <- vars$expr
    shown <- paste(deparse(expr, width.cutoff = 500L), collapse = " ")
    if (inherits(vars, "dbartsFormulaTerms")) {
      # a formula's forest(): term labels, and names of design columns
      columns <- integer(0L)
      for (label in c(stripBackticks(vars$labels), vars$columns)) {
        found <- resolveTermColumns(label, columnNames, termLabels)
        if (is.null(found)) {
          refuseUnknown(label)
        }
        columns <- c(columns, found)
      }
    } else {
      # a predictor by its own name, a column's, or a variable of its code:
      # log(x1) is a predictor whose name is not its variable's
      known <- unique(c(
        ".",
        predictors,
        unlist(lapply(
          c(attr(data@x, "term.labels"), columnNames),
          function(name) {
            tryCatch(all.vars(str2lang(name)), error = function(e) NULL)
          }
        ))
      ))
      if (any(all.vars(expr) %in% known)) {
        columns <- readSelectionTerms(expr, data, refuseUnknown)
      } else if (isTRUE(vars$evaluated)) {
        # the value is the selection, so what evaluating it warned of is the
        # caller's to hear
        for (kept in vars$warnings) {
          warning(kept)
        }
        vars <- vars$value
        if (inherits(vars, "formula")) {
          # a formula held in a variable: its right-hand side, as written
          if (length(vars) != 2L) {
            refuseHeldFormula(expr, owner)
          }
          columns <- readSelectionTerms(vars[[2L]], data, refuseUnknown)
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
          owner,
          ", '",
          shown,
          "': ",
          if (is.null(vars$error)) "it is no selection" else vars$error,
          call. = FALSE
        )
      }
    }
  }

  if (!is.null(columns)) {
    columns <- sort(unique(columns))
    if (length(columns) == 0L) {
      stop(
        owner,
        ", '",
        shown,
        "', leaves the forest no predictor to split on",
        call. = FALSE
      )
    }
  } else {
    if ((is.logical(vars) || is.factor(vars)) && !anyNA(vars)) {
      stop(
        owner,
        " is a ",
        if (is.factor(vars)) "factor" else "logical",
        "; a selection is names or positions. A multiplier is given as ",
        "'basis ='",
        call. = FALSE
      )
    }
    if (is.atomic(vars) && !anyNA(vars) && anyDuplicated(vars) > 0L) {
      repeated <- vars[[anyDuplicated(vars)]]
      if (!is.character(repeated)) {
        position <- suppressWarnings(as.integer(repeated))
        if (
          !is.na(position) &&
            position >= 1L &&
            position <= length(columnNames) &&
            nzchar(columnNames[position])
        ) {
          repeated <- columnNames[position]
        }
      }
      stop(
        owner,
        " selects predictors of the fit and names '",
        format(repeated),
        "' more than once; name each once. A multiplier is given as 'basis ='",
        call. = FALSE
      )
    }
    if (is.character(vars) && !anyNA(vars)) {
      if (!all(nzchar(vars))) {
        stop(
          owner,
          " has an empty name; a predictor with no ",
          "name is selected by position",
          call. = FALSE
        )
      }
      for (name in vars) {
        if (sum(columnNames == name) > 1L) {
          refuseSharedName(name, sum(columnNames == name))
        }
      }
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
##
## A declaration is a value handed over, kept as it is, or code
## (captureForestBasis), which the fit reads against its data
## (readForestBasis); the two are told apart by class and never by content.
forestBasisDeclarations <- function(forests) {
  if (!is.list(forests) || length(forests) < 1L) {
    return(NULL)
  }
  if (!all(vapply(forests, inherits, logical(1L), "dbartsForest"))) {
    return(NULL)
  }
  lapply(forests, function(spec) spec$basis)
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
## it must match is the fit's, checked against it there. The class is read by
## basisClass, and a class a basis does not take is the caller's to refuse
## before it comes here; a level no row takes expands to an all-zero column.
## The texts name the basis by `argument`, the caller's word for it.
expandForestBasis <- function(
  basis,
  atPrediction = FALSE,
  allowEmptyLevels = atPrediction,
  argument = "basis"
) {
  if (is.null(basis)) {
    return(NULL)
  }
  read <- basisClass(basis)
  if (anyNA(basis)) {
    stop("a '", argument, "' cannot be NA")
  }
  if (!is.null(read$factor)) {
    basis <- read$factor
    if (!atPrediction && nlevels(basis) < 2L) {
      stop("a '", argument, "' factor must have at least two levels")
    }
    codes <- as.integer(basis)
    expanded <- matrix(0, length(codes), nlevels(basis))
    expanded[cbind(seq_along(codes), codes)] <- 1
    return(expanded)
  }
  if (is.na(read$kind)) {
    stop("a '", argument, "' cannot be ", read$what)
  }
  basis <- as.matrix(basis)
  storage.mode(basis) <- "double"
  if (ncol(basis) < 1L) {
    stop("a '", argument, "' must have at least one column")
  }
  if (!all(is.finite(basis))) {
    stop("a '", argument, "' must be finite")
  }
  # the calibration divides by the median nonzero row norm, so a norm that
  # overflows, or underflows to zero on a nonzero row, poisons the sampler
  rowNorms <- sqrt(rowSums(basis * basis))
  if (any(!is.finite(rowNorms) | (rowNorms == 0 & rowSums(basis != 0) > 0L))) {
    stop(
      "a '",
      argument,
      "' row's norm is not representable; rescale the basis to moderate values"
    )
  }
  # an all-zero column has no observation for its amplitude to multiply, so it
  # is not a degenerate prior but a missing predictor wearing one
  if (!atPrediction && !allowEmptyLevels && any(colSums(basis != 0) == 0L)) {
    stop(
      "a '",
      argument,
      "' column of all zeros contributes nothing to a forest"
    )
  }
  basis
}

## One forest's basis at new rows, from the record its fit kept of a basis
## written as code (buildCodeBasis): R's terms rebuild every term from what
## the fitted rows gave it, a center, a scale or knots, through model.frame(),
## the fit-time levels of every factor are imposed, so that a level set in
## another order cannot move a coefficient to another level, and the columns
## are built as at the fit (basisColumns). model.frame() finds a variable the
## new rows lack in the terms' own environment, which for a column of the
## data would be the fitted rows' values taken for the new rows, so each such
## variable is looked at first: a column of the fit's data and any other
## variable with a value for every fitted row are refused by name, and
## anything else, a number written beside a column among it, is used as it is
## found, as lm() uses it.
replayForestBasis <- function(term, newdata, index) {
  terms <- term$terms
  env <- environment(terms)
  vars <- basisVariables(terms[[2L]])
  absent <- vars[vars %not_in% names(newdata)]
  # a column of the fit's data is a column of the new rows; any other variable
  # with a value for every fitted row has none for the new ones
  perRow <- vapply(
    absent,
    function(name) {
      if (name %in% term$columns) {
        return(TRUE)
      }
      found <- tryCatch(eval(as.name(name), env), error = function(e) e)
      inherits(found, "error") || (term$rows > 1L && NROW(found) == term$rows)
    },
    NA
  )
  missingVars <- absent[perRow]
  if (length(missingVars) > 0L) {
    stop(
      "'newdata' is missing ",
      if (length(missingVars) > 1L) "variables" else "variable",
      " '",
      toString(missingVars),
      "', required by forest ",
      index,
      "'s basis (",
      term$label,
      "); supply ",
      if (length(missingVars) > 1L) "them" else "it",
      ", or give that basis at the new rows with 'bases ='"
    )
  }
  frame <- stats::model.frame(
    terms,
    newdata,
    na.action = stats::na.pass,
    xlev = term$xlev
  )
  basisColumns(
    terms,
    frame,
    term$label,
    levels = term$levels,
    atPrediction = TRUE
  )$basis
}

## The full (pre-'subset') row count of a formula fit's data, and the exact
## rows 'subset' keeps in it, in the order it keeps them: the rows of a data
## frame, and without one the length of the first variable the formula names
## ('.' names no single column, so it is skipped for a real one), as for a
## basis (fitBasisRows); and `subset`, the VALUE the data object took of
## the caller's expression, used as an ordinary subscript into the full row
## sequence, row names among it - the reading stats::model.frame() itself
## gives 'subset'. The caller's expression is evaluated once, by the data
## object, and never here: evaluated again it could select other rows than
## the fit's. NULL when there is no rule to apply: not a formula, or no
## 'subset' given.
resolveFormulaBasisSubset <- function(formula, data, subset) {
  if (!is.formula(formula) || is.null(subset)) {
    return(NULL)
  }
  # a variable the formula names need not have a value for each row, the k
  # of I(k * y) ~ x, so a data frame's own rows are counted where there is one
  if (is.data.frame(data)) {
    full <- nrow(data)
  } else {
    vars <- setdiff(all.vars(formula), ".")
    if (length(vars) == 0L) {
      return(NULL)
    }
    env <- environment(formula)
    first <- if (is.list(data) || is.environment(data)) {
      eval(as.name(vars[1L]), data, env)
    } else {
      eval(as.name(vars[1L]), env)
    }
    full <- NROW(first)
  }
  rows <- seq_len(full)
  if (is.character(subset) && is.data.frame(data)) {
    names(rows) <- row.names(data)
  }
  list(
    full = full,
    index = unname(rows[subset]),
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
    return(cutBasisRows(basis, subsetRows$index))
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
  if (!is.null(spec$power)) {
    value <- suppressWarnings(as.double(spec$power))
    if (length(value) != 1L || is.na(value) || value <= 0.0) {
      stop("forest 'power' must be a single positive number")
    }
    spec$power <- value
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

## How a refusal names forest `index` of a `forests` list: by the list's name
## where it gives one, by position otherwise.
forestKnobLabel <- function(forests, index) {
  name <- names(forests)[index]
  if (length(name) == 1L && !is.na(name) && nzchar(name)) {
    paste0("forest \"", name, "\"")
  } else {
    paste0("forest ", index)
  }
}

## A forest's own move settings, each resolved as the control's slot is and
## left NULL where the forest states none: a mixture filled against the
## package defaults (not the control's) and checked whole, a window kept as a
## double, and treeShift as the slot's tri-state logical. `label` names the
## forest in a refusal, as forestKnobLabel's or "varianceForest()".
resolveForestMoves <- function(spec, label) {
  refuse <- function(knob, message) {
    quoted <- paste0("'", knob, "'")
    stop(
      label,
      "'s ",
      if (startsWith(message, quoted)) {
        message
      } else {
        paste0(quoted, ": ", message)
      },
      call. = FALSE
    )
  }
  if (!is.null(spec$proposal.probs)) {
    probs <- tryCatch(
      resolveProposalProbs(spec$proposal.probs),
      error = function(e) refuse("proposal.probs", conditionMessage(e))
    )
    refusal <- proposalProbsRefusal(probs)
    if (!is.null(refusal)) {
      refuse("proposal.probs", refusal)
    }
    spec$proposal.probs <- probs
  }
  if (!is.null(spec$n.perturb.cuts)) {
    width <- spec$n.perturb.cuts
    refusal <- if (!is.numeric(width) || is.object(width)) {
      "'n.perturb.cuts' must be a single positive whole number or Inf"
    } else {
      perturbCutsRefusal(width)
    }
    if (!is.null(refusal)) {
      refuse("n.perturb.cuts", refusal)
    }
    spec$n.perturb.cuts <- as.double(width)
  }
  if (!is.null(spec$treeShift)) {
    spec$treeShift <- tryCatch(
      resolveTreeShift(spec$treeShift),
      error = function(e) refuse("treeShift", conditionMessage(e))
    )
  }
  spec
}

## The record of what each forest states of its move settings, as
## attr(control, "bartcore.moves") carries it to the bridge: `forests`, one
## entry for each of the numForests mean forests (a multinomial model's one
## entry setting every category forest), and `variance`, the variance
## forest's or NULL. An entry holds proposal.probs, n.perturb.cuts and, on a
## mean forest, levelGibbs, each NULL where the forest states none. NULL
## where no forest states any, so a fit stating none carries no record.
forestMoveRecord <- function(specs, numForests, variance) {
  entry <- function(spec, mode = TRUE) {
    out <- list(
      proposal.probs = spec$proposal.probs,
      n.perturb.cuts = spec$n.perturb.cuts
    )
    if (mode) {
      out["levelGibbs"] <- list(spec$treeShift)
    }
    out
  }
  states <- function(e) !all(vapply(e, is.null, NA))
  forests <- lapply(seq_len(numForests), function(index) {
    entry(if (index <= length(specs)) specs[[index]])
  })
  varianceEntry <- entry(variance, FALSE)
  if (!any(vapply(forests, states, NA)) && !states(varianceEntry)) {
    return(NULL)
  }
  list(
    forests = forests,
    variance = if (states(varianceEntry)) varianceEntry
  )
}

## A forest's 'amplitude': NULL, which draws its coefficient under its kind's
## default prior; fixed(), which holds it at 1, the basis carrying the coding;
## or student(df), which draws it under a Student-t prior at the default's
## median. Returns "fixed" for the hold and the student() object as it is. A
## bare fixed or student is the constructor called, as a constructor given
## where a value is.
validateForestAmplitude <- function(amplitude) {
  if (identical(amplitude, fixed)) {
    amplitude <- fixed()
  } else if (identical(amplitude, studentPrior)) {
    amplitude <- studentPrior()
  }
  if (is(amplitude, "dbartsStudentPrior")) {
    return(amplitude)
  }
  if (!is(amplitude, "dbartsFixedPrior")) {
    stop(
      "a forest's 'amplitude' must be fixed(), which holds its coefficient, ",
      "student(df), which draws it under a Student-t prior, or left out, ",
      "which draws it under its default prior",
      if (is(amplitude, "dbartsFamily")) {
        "; dbartsFamilies$student is the response's family"
      }
    )
  }
  value <- amplitude@value
  if (!is.numeric(value) || length(value) != 1L || value != 1) {
    stop(
      "'amplitude = fixed(",
      paste(deparse(value), collapse = ""),
      ")': a held coefficient is 1, and what it multiplies is the basis: ",
      "nothing, a numeric column as given, or 0 and 1 for the two levels of ",
      "a factor; fixed() takes no other value. Write fixed(), state the ",
      "forest's size with 'sd', and write another coding as the basis, as ",
      "basis = z - 0.5"
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

## A forest's 'sd', at creation and on $setLeafPrior: one number, finite and
## positive, a name on it dropped as base R drops one. Infinity states no
## prior the map can scale, so it is refused rather than carried into a leaf
## scale or a half-Cauchy median. A number is a numeric or integer, whatever
## class it carries.
validateForestSd <- function(sd) {
  if (isSingleNA(sd) || (is.numeric(sd) && length(sd) == 1L && is.nan(sd))) {
    stop("forest 'sd' must not be NA; leave it out for the default")
  }
  if (is(sd, "dbartsSdHyperprior")) {
    stop(
      "forest 'sd' must be a number, not invchi(): a hyperprior on a forest's ",
      "sd is not supported yet"
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
      ": a forest states one sd; to size a column differently, rescale it ",
      "in 'basis', as I(dose / 30)"
    )
  }
  if (is.na(sd) || !is.finite(sd) || sd <= 0) {
    stop("forest 'sd' must be positive and finite")
  }
  as.double(sd)
}

## The position of the forest with no basis among a model's forests, or 0
## where every forest has one. `hasBasis` says of each forest whether the
## model ends with a basis for it. This is the one place a forest's kind is
## asked: its tree count and tree prior default by it and the fitting
## function's 'interactions' and 'blocks' go to it, whatever its position. A
## model has at most one forest with no basis.
plainForest <- function(hasBasis) {
  bare <- which(!hasBasis)
  if (length(bare) > 1L) {
    stop(
      "forests ",
      paste(bare[-length(bare)], collapse = ", "),
      " and ",
      bare[length(bare)],
      " have no 'basis': a model has one forest with no multiplier, and ",
      "every other forest states a 'basis'",
      call. = FALSE
    )
  }
  if (length(bare) == 1L) bare else 0L
}

## The tree count and tree prior of a forest with a basis that states none.
multipliedForestDefaults <- list(n.trees = 50L, base = 0.25, power = 3)

## The attribute bart() and bartBT() stamp on the control they build with
## what their own caller stated (plainForestStated). That control always names
## a tree count and bartBT() always hands dbarts() a power and a base, so
## dbarts() cannot read either from what it is given. Not a bartcore.*
## attribute, which is fit state: the sampler specification takes it off the
## control it resolves.
plainStatedAttr <- "dbarts.stated"

## Whether a control states a tree count: its constructor's call named one,
## or the slot is not that constructor's default. A slot edited back to the
## default and a control made by new() read as not stated.
controlStatesTreeCount <- function(control) {
  "n.trees" %in%
    controlSuppliedSlots(control) ||
    !identical(control@n.trees, formals(dbarts::dbartsControl)$n.trees)
}

## The names under which a fitting function takes the arguments that are the
## forest with no basis's own, each in the order it is looked for. The split
## probabilities are not among them: they are the model's, judged by their
## value where no forest can take them.
plainForestArgNames <- list(
  n.trees = "n.trees",
  power = "power",
  base = "base",
  leaf.prior = c("leaf.prior", "k", "node.prior")
)

## What the caller of a fitting function stated of the arguments that are the
## forest with no basis's own: the tree count, the tree prior's power and base,
## and the leaf prior. Stated means named, the default's own value included.
## `written` holds the argument names the caller wrote, `control` the control
## the caller gave, NULL for none, and `argNames` that function's names for
## them. Returns, for each of "n.trees", "power", "base" and "leaf.prior" that
## is stated, the name it was stated under, "control" for a count the control
## states. `remedy` says whether the fitting function can be told to state
## these on a forest; where it cannot (bartBT() has no `forests` and refuses
## forest() terms, and bart() given a data object has its bases from the
## object) the value carries a "noForestRemedy" attribute, which the refusal
## reads.
plainForestStated <- function(
  control,
  written,
  argNames = plainForestArgNames,
  remedy = TRUE
) {
  firstWritten <- function(names) {
    hits <- names[names %in% written]
    if (length(hits) > 0L) hits[[1L]]
  }
  countName <- firstWritten(argNames$n.trees)
  stated <- c(
    n.trees = if (!is.null(countName)) {
      countName
    } else if (!is.null(control) && controlStatesTreeCount(control)) {
      "control"
    },
    power = firstWritten(argNames$power),
    base = firstWritten(argNames$base),
    leaf.prior = firstWritten(argNames$leaf.prior)
  )
  if (is.null(stated)) {
    stated <- character(0L)
  }
  if (!remedy) {
    attr(stated, "noForestRemedy") <- TRUE
  }
  stated
}

## The refusal of what the fitting function states for the forest with no
## basis in a model that has none: each of the tree count, the tree prior's
## power and base, the leaf prior, 'interactions' and 'blocks' that is stated
## has no forest to belong to. In that order, the first one stated. `stated`
## is plainForestStated's value, `treeValue` a function of "power" or "base"
## returning the value the caller gave, which the remedy repeats.
refuseStatedWithNoPlainForest <- function(
  stated,
  interactions,
  blocks,
  treeValue = function(name) multipliedForestDefaults[[name]]
) {
  # where the fitting function cannot take an argument on a forest, the caller
  # is told what it can do: leave the argument out, or fit with dbarts()
  noRemedy <- isTRUE(attr(stated, "noForestRemedy", exact = TRUE))
  refuse <- function(subject, what, remedy, leaveOut = "leave it out") {
    if (noRemedy) {
      remedy <- paste0(leaveOut, ", or fit with dbarts()")
    }
    stop(
      subject,
      " is ",
      what,
      " the forest with no basis, and every forest of this model has a ",
      "basis; ",
      remedy,
      call. = FALSE
    )
  }
  given <- function(argument) {
    paste0("'", argument, "' given to the fitting function")
  }
  if ("n.trees" %in% names(stated)) {
    fromControl <- identical(stated[["n.trees"]], "control")
    refuse(
      if (fromControl) {
        "the control's 'n.trees'"
      } else {
        given(stated[["n.trees"]])
      },
      "the tree count of",
      paste0(
        "state a count on a forest, as forest(x1, basis = a, n.trees = 100)",
        if (fromControl) ", and leave it out of dbartsControl()"
      ),
      if (fromControl) "leave it out of dbartsControl()" else "leave it out"
    )
  }
  for (name in c("power", "base")) {
    if (name %in% names(stated)) {
      refuse(
        given(stated[[name]]),
        "the tree prior of",
        paste0(
          "state it on a forest, as forest(x1, basis = a, ",
          name,
          " = ",
          treeValue(name),
          ")"
        )
      )
    }
  }
  if ("leaf.prior" %in% names(stated)) {
    refuse(
      given(stated[["leaf.prior"]]),
      "the leaf prior of",
      "state a forest's size on the forest, as forest(x1, basis = a, sd = 2)"
    )
  }
  if (!is.null(interactions)) {
    refuse(
      given("interactions"),
      "a constraint on",
      paste0(
        "state it on a forest, as forest(x1, basis = a, interactions = ",
        "interactions(max.order = 1))"
      )
    )
  }
  if (!is.null(blocks)) {
    refuse(
      given("blocks"),
      "a constraint on",
      paste0(
        "state it on a forest, as forest(x1 + x2, basis = a, blocks = ",
        "blocks(list(\"x1\", \"x2\")))"
      )
    )
  }
  invisible(NULL)
}

## The fitting function's settings that are the forest with no basis's own and
## that the forest states too: its 'n.trees', 'power' and 'base', each one
## setting written twice and refused by name, the forest's statement of one
## leaving the others the fitting function's. A count the control states is
## not the fitting function's to give twice, and is left to the forest's.
## `stated` is plainForestStated's value, `plainSpec` that forest's resolved
## knobs and `where` how the refusal shows it.
refuseStatedOnPlainForest <- function(stated, plainSpec, where) {
  for (name in c("n.trees", "power", "base")) {
    isStated <- name %in% names(stated)
    if (
      isStated &&
        !identical(stated[[name]], "control") &&
        !is.null(plainSpec[[name]])
    ) {
      refuseGivenTwice(name, where)
    }
  }
  invisible(NULL)
}

## The refusal of one of 'n.trees', 'power' and 'base' given to the fitting
## function and to the forest with no basis, shown as `where`.
refuseGivenTwice <- function(name, where) {
  stop(
    "'",
    name,
    "' is given to the fitting function and to the forest with no basis (",
    where,
    "), which are the same ",
    if (name == "n.trees") "count" else "setting",
    "; give one",
    call. = FALSE
  )
}

## Resolve a `forests` declaration into the per-forest knobs a sampler
## specification carries, or NULL for the single-forest path. Returns a LIST of
## K validated knob lists, one per forest, so nothing downstream is keyed on
## two. Everything the engine cannot honor refuses here, by name, rather than
## being dropped: one forest with no basis at most, and the amplitude knobs
## only where an amplitude exists at all. `interactions` and `blocks` are the fit's own top-level arguments,
## which address the forest with NO BASIS under the same spelling a forest()
## uses, so supplying both is ambiguous rather than layered.
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
  resolved <- lapply(seq_along(forests), function(index) {
    resolveForestMoves(
      validateForestKnobs(forests[[index]]),
      forestKnobLabel(forests, index)
    )
  })
  numForests <- length(resolved)

  withBasis <- vapply(
    seq_len(numForests),
    function(index) {
      !is.null(resolved[[index]]$basis) ||
        (length(hasBasis) >= index && hasBasis[index])
    },
    NA
  )
  plain <- plainForest(withBasis)

  first <- resolved[[1L]]
  if (numForests == 1L) {
    if (is(first$amplitude, "dbartsStudentPrior") && !any(hasBasis)) {
      stop(
        "this model has one forest, which has no coefficient to give a ",
        "prior; a lone forest draws its size through the fitting ",
        "function's leaf prior, as leaf.prior = normal(k = chi(1.5, 2)) or ",
        "normal(sd = invchi(1.5, s)) with s a size in the response's units"
      )
    }
    if (!is.null(first$amplitude) && !any(hasBasis)) {
      stop(
        "this model has one forest, which has no coefficient to hold; ",
        "'amplitude' needs a model of several forests"
      )
    }
  }
  for (name in c("interactions", "blocks")) {
    top <- if (name == "interactions") interactions else blocks
    if (plain > 0L && !is.null(resolved[[plain]][[name]]) && !is.null(top)) {
      stop(
        "'",
        name,
        "' is given to the fitting function and to the forest with no basis, ",
        "which are the same constraint; give one",
        call. = FALSE
      )
    }
  }
  resolved
}

## The four doubles attr(control, "bartcore.forests")$params carries FOR EACH
## FOREST, in the order the C bridge reads them: the forest's tree count and
## structure prior, and the amplitude update flag, 0 where the coefficient is
## held. The first three are what the forest runs under, at every position:
## its own statement, else the default of its kind, `plain` holding the
## fitting function's count, base and power for the forest with no basis. The
## bridge reads forest 1's three from the control and the model's tree prior,
## which hold the same numbers.
##
## Nothing of a forest's size is here: its `sd` goes to the engine as written
## (forestStatements), and the engine holds the calibration map and every
## default.
forestParams <- function(specs, hasBasis, plain) {
  declared <- function(value, default) {
    if (is.null(value)) default else value
  }
  lapply(seq_along(specs), function(index) {
    spec <- specs[[index]]
    kind <- if (hasBasis[index]) multipliedForestDefaults else plain
    as.double(c(
      declared(spec$n.trees, kind$n.trees),
      declared(spec$base, kind$base),
      declared(spec$power, kind$power),
      if (identical(spec$amplitude, "fixed")) 0 else 1
    ))
  })
}

## What each forest states of its size and its coefficient's prior, as the
## forests' record carries it to the engine: `sd`, in the response's units, as
## written and NA where none is stated, and `amplitude.df`, a student()
## coefficient's df and NA elsewhere, present only where some forest states
## one, so that a record made without it is as it was.
forestStatements <- function(specs) {
  statements <- list(
    sd = vapply(
      specs,
      function(spec) if (is.null(spec$sd)) NA_real_ else spec$sd,
      0
    )
  )
  df <- vapply(
    specs,
    function(spec) {
      if (is(spec$amplitude, "dbartsStudentPrior")) {
        spec$amplitude@df
      } else {
        NA_real_
      }
    },
    0
  )
  if (!all(is.na(df))) {
    statements$amplitude.df <- df
  }
  statements
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
## The value of the retired tree.prior argument of dbarts() and dbartsSpec(),
## which leaves with it; its class stays as the model's tree-prior record.
cgm <- function(power = 2, base = 0.95, split.probs = NULL) {
  # NULL is absent, so it is the default
  if (is.null(power)) {
    power <- 2
  }
  if (is.null(base)) {
    base <- 0.95
  }
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

## A leaf prior's 'sd': NULL, a positive number, or an invchi() distribution
## on it. A distribution on k and the string forms k keeps for 0.9-x are
## refused by name.
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
      "'sd' takes a number or invchi(), a hyperprior on the sd itself; a ",
      "hyperprior on k is spelled k = chi(), and the same prior on the sd is ",
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

## The model's encoding of a leaf prior, its prior.scale and its leaf
## hyperprior, from whichever of 'k' and 'sd' it names. A named sd rides a
## reference k of 2: a fixed sd x is prior.scale 2x beside k fixed at 2, and
## invchi(df, c) is prior.scale 2c beside chi(df, 2). The pair is an encoding,
## not what the engine runs: the bridge divides it back to the sd, exactly,
## and the engine states that sd as k against the data's scale - a fixed k of
## k.scale / x, or k ~ chi(df, k.scale / c) started there, so the chain starts
## at the named spread - and restates it whenever the response mapping moves.
## Every k spelling leaves prior.scale NA and reaches the engine as written.
## invchi(df, 0) is the improper limit, which no scale can state and
## chi(df, Inf) is. A multi-forest fit refuses a named sd before any of this:
## each of its forests states its size on the forest.
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
      "a model of several forests states each forest's size on the forest: ",
      "write the size of the forest with no basis as forest(x1 + x2, sd = 2), ",
      "which is in the response's units as leaf.prior = normal(sd = ) is, and ",
      "leave 'sd' out of 'leaf.prior'"
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
      "supply a fixed numeric sd (the truncated leaf prior has no chi-k update)"
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
    refuseQuotedK(k)
    if (is(k, "dbartsSdHyperprior")) {
      stop(
        "'k' takes a number or chi(), a hyperprior on k; invchi() is a ",
        "hyperprior on the sd and is spelled sd = invchi()"
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

## A coefficient's Student-t prior, student(df) at a forest's 'amplitude' and
## dbartsPriors$student outside it; its own name here, the response family
## being student. student() is the half-Cauchy.
studentPrior <- function(df = 1) {
  refusal <- studentPriorRefusal(df)
  if (!is.null(refusal)) {
    stop(refusal, call. = FALSE)
  }
  new("dbartsStudentPrior", df = as.double(df))
}

## Why df is no coefficient's Student-t df, NULL where it is one.
studentPriorRefusal <- function(df) {
  if (
    !is.numeric(df) ||
      is.object(df) ||
      length(df) != 1L ||
      !is.finite(df) ||
      df < 1
  ) {
    paste0(
      "a coefficient's student(df) takes one finite df of 1 or more; ",
      "student() is student(1), the half-Cauchy"
    )
  }
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

## DART's Dirichlet prior on the split probabilities, the value of
## split.probs that selects it. Its record is the tree prior's, whose power
## and base, held here at the defaults, the fitting function writes.
dart <- function(
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
    power = 2,
    base = 0.95,
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
## tree-structure prior, and on the forest with no basis they are the fitting
## function's arguments of the same names, so stating one at both levels is
## refused, each name on its own; sd is its size in the response's units (the latent
## index's under probit and logistic), which the engine divides by the unit
## it derives at creation; amplitude = fixed() holds its coefficient at 1, the
## basis carrying the coding, and left out draws it; interactions and blocks are this forest's own constraints - the
## arguments of the same names on the fitting function are the FIRST forest's;
## proposal.probs, n.perturb.cuts and treeShift are its own move settings,
## each the control's where not stated (resolveForestMoves).
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
  proposal.probs = NULL,
  n.perturb.cuts = NULL,
  treeShift = NULL
) {
  # the arguments as the caller gave them, a forwarded '...' spelled out
  given <- match.call(function(...) NULL, sys.call(), expand.dots = FALSE)$...
  if (numUnnamed(given) > 1L) {
    refuseSecondUnnamed()
  }
  place <- callingPlace(parent.frame())
  matched <- match.call()
  # An argument that a function of base R wrote for the caller - lapply()'s
  # X[[i]], Map()'s dots[[2L]][[1L]] - is no code of the caller's: what it
  # holds is a value handed over, and neither its text nor its names mean
  # anything. One forwarded through dots was written by the caller of the
  # function that forwards it.
  written <- substitute(vars)
  if (
    is.language(written) &&
      !isDotsReference(matched[["vars"]]) &&
      isLoopMachinery(place, written)
  ) {
    written <- vars
  }
  written <- captureForestVars(written)
  if (inherits(written, "dbartsForestCode")) {
    # the argument's value here and now, taken once: whatever the caller's
    # variables hold later, a forest built in a loop or by lapply() keeps the
    # selection it was given. Code that cannot be evaluated here, terms over
    # predictors among it, is kept as code alone
    written <- takeAtCall(written, function() vars, place)
  }
  # the basis as its caller wrote it and where: one forwarded through a
  # wrapper's dots is read where the wrapper's caller wrote it
  declared <- matched[["basis"]]
  if (is.language(declared)) {
    recovered <- recoverForwardedArgument(declared, parent.frame())
    if (isDotsReference(recovered$expr)) {
      declared <- forwardedBasis(substitute(basis), function() basis, place)
    } else if (isLoopMachinery(callingPlace(recovered$env), recovered$expr)) {
      # such a function indexes its variables and never passes one whole, so
      # a bare name was handed to it, and stands for a variable of its own
      if (is.name(recovered$expr)) {
        stop(
          "forest()'s 'basis' is the name '",
          recovered$expr,
          "', which the function of base R that called forest() was handed ",
          "and has a variable of its own by, so the name cannot be read as ",
          "the caller's; hand the column itself over, or write the call of ",
          "forest() in a function of your own",
          call. = FALSE
        )
      }
      declared <- basis
    } else {
      declared <- recovered$expr
      place <- callingPlace(recovered$env)
    }
  }
  declared <- captureForestBasis(declared, place)
  structure(
    list(
      basis = declared,
      vars = written,
      n.trees = n.trees,
      base = base,
      power = power,
      sd = sd,
      interactions = interactions,
      blocks = blocks,
      amplitude = amplitude,
      proposal.probs = proposal.probs,
      n.perturb.cuts = n.perturb.cuts,
      treeShift = treeShift
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

## An argument of forest() that is read as code, its predictors and its
## basis, is held in one of three states, which every reader tells apart by
## class:
##
## - a value handed over (names, positions, a column, NULL): the value
##   itself, with no wrapper;
## - code a caller wrote in a call of forest() (forestCode): the code and,
##   once takeAtCall() has run, either the value it had at the call with the
##   warnings that raised, or the reason it had none;
## - the terms of a formula's forest() (forestFormulaTerms): the code and the
##   term labels read against the data.
##
## The last two share the class "dbartsForestTerms". The predictors' code
## keeps no environment: in a call its value is taken there and nothing is
## looked up again, and a formula's terms are labels. So a forest carries no
## frame of the function that built it. A basis's code (basisCode) keeps the
## environment its names are looked up in.

## Code written in a call of forest(), not yet evaluated.
forestCode <- function(expr) {
  structure(
    list(expr = expr),
    class = c("dbartsForestCode", "dbartsForestTerms")
  )
}

## The value code has where and when forest() is called, taken once and
## quietly by `evaluate`, and kept beside the code: `value`, `evaluated` and
## the `warnings` it raised, which belong to the fit that uses the value. Code
## that cannot be evaluated there keeps the `error` and the names it uses
## that nothing binds in `env`, the place of the call, as `unbound`; that is
## no error of forest()'s.
takeAtCall <- function(code, evaluate, env) {
  warned <- list()
  taken <- tryCatch(
    withCallingHandlers(list(evaluate()), warning = function(w) {
      warned[[length(warned) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  if (inherits(taken, "error")) {
    symbols <- all.vars(code$expr)
    code$unbound <- symbols[!vapply(symbols, exists, NA, envir = env)]
    code$error <- conditionMessage(taken)
  } else {
    code["value"] <- taken
    code$evaluated <- TRUE
    # raised when the value is what a fit uses, and never when the code
    # turns out to be terms
    code$warnings <- warned
  }
  code
}

## The terms of a formula's forest(): `labels`, the term labels its code
## gives against the data, and `columns`, names of columns the built design
## is to resolve. `named` marks terms that came from names given by value,
## whose labels are then the code to write for them.
forestFormulaTerms <- function(expr, labels, columns = NULL, named = NULL) {
  terms <- list(expr = expr, labels = labels)
  terms$columns <- columns
  terms$named <- named
  structure(terms, class = c("dbartsFormulaTerms", "dbartsForestTerms"))
}

## forest()'s first argument as written. A value (names, positions, NULL) is
## kept as it is; code is kept unevaluated, and forest() puts beside it the
## value it has at the call. A one-sided formula written in place is its
## right-hand side, the tilde accepted and never required; one with a
## left-hand side is no selection and is refused.
captureForestVars <- function(expr, owner = "forest()'s first argument") {
  if (!is.language(expr)) {
    return(expr)
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("~"))) {
    if (length(expr) != 2L) {
      stop(
        owner,
        " is a formula with a left-hand side; write the predictors as the ",
        "right-hand side alone, as x1 + x2 or ~ x1 + x2",
        call. = FALSE
      )
    }
    expr <- expr[[2L]]
  }
  forestCode(expr)
}

## Whether `expr`, an argument of forest() written in the frame `place`, is
## one a function of base R wrote: lapply(), sapply(), vapply(), Map(),
## mapply() and their like call forest() for the caller with an argument of
## their own making, which indexes a variable of their own frame, X[[i]] or
## dots[[2L]][[1L]]. A frame of the caller's, of a package's function or of a
## wrapper is not such a frame. And a name or a call that such a function was
## handed and passes on as it is, as mapply() passes on its MoreArgs, names
## no variable of that frame: it is the caller's code, as the same name or
## call handed over through do.call() is.
isLoopMachinery <- function(place, expr) {
  identical(topenv(place), .BaseNamespaceEnv) &&
    any(vapply(all.names(expr), exists, NA, envir = place, inherits = FALSE))
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
## NULL (40 trees; the base and power the mean forest resolved to, the fitting
## function's or a lone forest's own), and
## proposal.probs/n.perturb.cuts to the control's. It takes no treeShift: the
## variance forest takes no level step. Validated at fit time. Not exported,
## like forest().
varianceForest <- function(
  vars = NULL,
  n.trees = NULL,
  base = NULL,
  power = NULL,
  proposal.probs = NULL,
  n.perturb.cuts = NULL
) {
  # the predictors as forest()'s are read: code is kept with the value it has
  # here, and read as terms over the fit's predictors once they are known
  written <- captureForestVars(
    substitute(vars),
    "varianceForest()'s first argument"
  )
  if (inherits(written, "dbartsForestCode")) {
    written <- takeAtCall(
      written,
      function() vars,
      callingPlace(parent.frame())
    )
    vars <- written
  }
  structure(
    list(
      vars = vars,
      n.trees = n.trees,
      base = base,
      power = power,
      proposal.probs = proposal.probs,
      n.perturb.cuts = n.perturb.cuts
    ),
    class = "dbartsVarianceForest"
  )
}

## One knob of a forest specification as its statement reads: code as the
## caller wrote it, a constructor's object by its name, a value deparsed on
## one line, and `default` where the knob is not stated.
describeKnob <- function(value, default) {
  if (is.null(value)) {
    return(default)
  }
  if (is.list(value) && !is.null(value$expr)) {
    value <- value$expr
  } else if (
    isS4(value) || (is.list(value) && !is.null(attr(value, "class")))
  ) {
    name <- sub("Prior$", "", sub("^dbarts", "", class(value)[[1L]]))
    return(paste0(tolower(substring(name, 1L, 1L)), substring(name, 2L), "()"))
  }
  text <- paste(deparse(value, width.cutoff = 500L), collapse = " ")
  if (nchar(text) > 60L) paste0(substr(text, 1L, 57L), "...") else text
}

## The knobs of a specification, one line each.
formatKnobs <- function(x, defaults) {
  paste0(
    format(names(defaults)),
    " = ",
    vapply(
      names(defaults),
      function(knob) describeKnob(x[[knob]], defaults[[knob]]),
      ""
    )
  )
}

format.dbartsVarianceForest <- function(x, ...) {
  formatKnobs(
    x,
    c(
      vars = "<all columns>",
      n.trees = "<default 40>",
      base = "<mean forest's>",
      power = "<mean forest's>",
      proposal.probs = "<control's>",
      n.perturb.cuts = "<control's>"
    )
  )
}

print.dbartsVarianceForest <- function(x, ...) {
  cat("dbarts variance forest specification\n")
  cat(paste0("  ", format(x)), sep = "\n")
  invisible(x)
}

## A forest's statement: the knobs it states, one line each, named in the
## order forest() takes them.
format.dbartsForest <- function(x, ...) {
  knobs <- c(
    "vars",
    "basis",
    "sd",
    "n.trees",
    "base",
    "power",
    "amplitude",
    "interactions",
    "blocks",
    "proposal.probs",
    "n.perturb.cuts",
    "treeShift"
  )
  stated <- knobs[!vapply(knobs, function(knob) is.null(x[[knob]]), NA)]
  if (length(stated) == 0L) {
    return("every knob at its default")
  }
  formatKnobs(x, stats::setNames(rep("", length(stated)), stated))
}

print.dbartsForest <- function(x, ...) {
  cat("dbarts forest specification\n")
  cat(paste0("  ", format(x)), sep = "\n")
  invisible(x)
}

## The exported face of the prior constructors: one object, so that no
## generic name (normal, chisq, fixed, chi) enters the search path to be
## masked by or to mask another package by attach order. Inside the
## leaf.prior argument of the fitting functions - and inside the 'sigma'
## argument of the families that draw a residual scale - the same
## constructors are available by bare name, as dart() is inside split.probs.
## cgm() stays here as the value of the retired tree.prior until it leaves.
dbartsPriors <- list(
  cgm = cgm,
  dart = dart,
  normal = normal,
  linear = linear,
  gp = gp,
  chisq = chisq,
  fixed = fixed,
  chi = chi,
  invchi = invchi,
  student = studentPrior
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
## dbartsForests and the prior constructors a forest's 'amplitude' takes,
## which stay under dbartsPriors for the caller who names them outside.
forestConstructors <- c(
  dbartsForests,
  dbartsPriors[c("fixed", "student", "invchi")]
)

## Each door argument taking a forest constructor, and the vocabulary it
## resolves over. In the order the doors forced them before they resolved by
## name, so that an error names the same argument.
FOREST_ARGUMENT_VOCABULARIES <- list(
  forests = c(
    "forest",
    "interactions",
    "blocks",
    "fixed",
    "student",
    "invchi"
  ),
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
