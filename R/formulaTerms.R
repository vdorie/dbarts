## A forest() call written at the top of a formula's '+' chain declares one
## forest of the model, the same declaration a forests = list(forest(...))
## entry makes, and shares its ingestion between dbarts() and bart() (which
## reaches it only by forwarding its own formula, unchanged, into dbarts()).
## A forest's first argument is the right-hand side of a model formula of its
## own, the predictors it splits on; the fit's predictors are every forest's.
## Everything decidable from the unevaluated formula refuses without touching
## data; a basis's degenerate columns need real values and are checked once
## the model frame exists.

## The families a multiplier model's calibration map has no per-forest
## amplitude block for (the engine's aft/ordinal/nbinom refusal), plus the
## families whose OWN ingestion the amplitude coupling is incompatible with
## (hazard's person-period row expansion, multinomial's K-forest softmax,
## hurdle.lognormal's two-sampler composition) - closed here so a forest()
## term is refused by name instead of surfacing downstream as a row-count or
## type mismatch that names none of family, term, or reason.
TERM_UNSUPPORTED_FAMILIES <- c(
  "hazard",
  "hazard.probit",
  "hazard.logistic",
  "multinomial",
  "hurdle.lognormal",
  "aft",
  "ordinal",
  "nbinom"
)

## The term grammar names forest only: forest() is not exported, so no other
## spelling reaches it.
isForestCall <- function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("forest"))
}

## The arguments of a call that are written; x[, 1] has one that is not, and
## reading that one as a value is an error.
writtenArguments <- function(expr) {
  args <- as.list(expr)[-1L]
  written <- vapply(
    seq_along(args),
    function(i) !identical(args[[i]], quote(expr = )),
    NA
  )
  args[written]
}

containsForestCall <- function(expr) {
  if (isForestCall(expr)) {
    return(TRUE)
  }
  if (is.call(expr)) {
    for (a in writtenArguments(expr)) {
      if (containsForestCall(a)) return(TRUE)
    }
  }
  FALSE
}

## The cheap guard bart2's multinomial and hurdle.lognormal branches use
## before their own dispatch: neither reaches the shared dbarts() ingestion
## below, so each must refuse a term itself rather than silently walking past
## one into a formula it never expected to contain a forest() call.
formulaHasForestTerm <- function(formula) {
  is.formula(formula) && containsForestCall(formula)
}


## A piece of the caller's code as one line, for a message.
shownCode <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

## The first forest() call of an expression, reading left to right.
firstForestCall <- function(expr) {
  if (isForestCall(expr)) {
    return(expr)
  }
  if (is.call(expr)) {
    for (a in writtenArguments(expr)) {
      found <- firstForestCall(a)
      if (!is.null(found)) return(found)
    }
  }
  NULL
}

## A door that fits no model of several forests refuses a forest() term by
## name, ahead of its own reading of the formula, which would stop on the
## term with a message naming neither.
refuseForestTerm <- function(formula, caller) {
  if (!formulaHasForestTerm(formula)) {
    return(invisible(NULL))
  }
  stop(
    caller,
    "() does not take a forest() term ('",
    shownCode(firstForestCall(formula)),
    "'): the forests of a model are written in the formula of bart() or ",
    "dbarts(), or in a 'forests' list",
    call. = FALSE
  )
}

## A forest() call's arguments under forest()'s own names, unevaluated. An
## argument forest() does not have is R's own error.
forestCallArguments <- function(call) {
  as.list(match.call(forest, call))[-1L]
}

## The basis a crossed operand stands for: a name, factor() of a name, any
## other call on columns, or a parenthesised sum of names, which is those
## columns side by side. NULL where the operand is itself crossed or holds a
## forest(), for which no one forest() can be written.
crossedOperandBasis <- function(expr) {
  if (containsForestCall(expr)) {
    return(NULL)
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L
  ) {
    inner <- expr[[2L]]
    members <- all.vars(inner)
    if (
      length(members) >= 2L &&
        "." %not_in% members &&
        all(setdiff(all.names(inner), members) == "+")
    ) {
      return(as.call(c(quote(cbind), lapply(members, as.name))))
    }
    return(crossedOperandBasis(inner))
  }
  if (
    is.call(expr) &&
      is.name(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% c(":", "*")
  ) {
    return(NULL)
  }
  if (identical(expr, as.name("."))) {
    return(NULL)
  }
  expr
}

## The refusal of a forest() crossed with another term by ':' or '*'.
## `forestSide` is the forest() operand and `otherSide` what it is crossed
## with, both NULL where the forest() lies deeper in the crossing.
refuseCrossedForest <- function(expr, forestSide = NULL, otherSide = NULL) {
  lead <- paste0(
    "'",
    shownCode(expr),
    "': a forest() is not crossed with another term; "
  )
  basis <- if (is.null(forestSide)) NULL else crossedOperandBasis(otherSide)
  if (is.null(basis)) {
    stop(
      lead,
      "a forest's multiplier is its 'basis' argument, as ",
      "forest(x1 + x2, basis = ~ z)",
      call. = FALSE
    )
  }
  written <- as.list(forestSide)[-1L]
  given <- tryCatch(
    names(forestCallArguments(forestSide)),
    error = function(e) {
      names(written)
    }
  )
  if ("basis" %in% given) {
    stop(
      lead,
      "its multiplier is its 'basis' argument, which this one already states",
      call. = FALSE
    )
  }
  argumentNames <- names(written)
  if (is.null(argumentNames)) {
    argumentNames <- rep("", length(written))
  }
  arguments <- paste0(
    ifelse(nzchar(argumentNames), paste0(argumentNames, " = "), ""),
    vapply(written, shownCode, "")
  )
  stop(
    lead,
    "a forest's multiplier is its 'basis' argument: write forest(",
    paste(
      c(arguments, paste0("basis = ~ ", shownCode(basis))),
      collapse = ", "
    ),
    ")",
    call. = FALSE
  )
}

## Walks 'formula' for forest() terms. A forest() is a term of the right-hand
## side's top-level '+' chain and nothing else: crossed by ':' or '*', at any
## depth, it is refused with the form to write; on the left-hand side, in a
## removal, or anywhere it could only be reached by evaluating an expression,
## it is refused by name. Returns NULL when no forest() is present anywhere,
## leaving the caller's formula untouched; otherwise the forest() calls in
## the order written and `plain`, the right-hand side without them, NULL when
## nothing is left.
walkFormulaTerms <- function(formula) {
  if (!containsForestCall(formula)) {
    return(NULL)
  }
  hits <- list()

  isOperator <- function(expr, operators) {
    is.call(expr) &&
      length(expr) == 3L &&
      is.name(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% operators
  }
  refuseInside <- function(expr) {
    stop(
      "a forest() term must appear as a top-level additive term, not inside '",
      shownCode(expr),
      "'",
      call. = FALSE
    )
  }
  refuseCrossed <- function(expr) {
    if (!is.call(expr)) {
      return(invisible(NULL))
    }
    if (isOperator(expr, c(":", "*")) && containsForestCall(expr)) {
      left <- expr[[2L]]
      right <- expr[[3L]]
      if (isForestCall(right)) {
        refuseCrossedForest(expr, right, left)
      }
      if (isForestCall(left)) {
        refuseCrossedForest(expr, left, right)
      }
      refuseCrossedForest(expr)
    }
    for (a in writtenArguments(expr)) {
      refuseCrossed(a)
    }
    invisible(NULL)
  }
  walk <- function(expr) {
    if (isForestCall(expr)) {
      hits[[length(hits) + 1L]] <<- expr
      return(NULL)
    }
    if (isOperator(expr, "+")) {
      left <- walk(expr[[2L]])
      right <- walk(expr[[3L]])
      if (is.null(left)) {
        return(right)
      }
      if (is.null(right)) {
        return(left)
      }
      return(call("+", left, right))
    }
    # a trailing removal is the fit's, as - 1 is; a forest() is not removed
    if (isOperator(expr, "-") && !containsForestCall(expr[[3L]])) {
      left <- walk(expr[[2L]])
      if (is.null(left)) {
        return(call("-", expr[[3L]]))
      }
      return(call("-", left, expr[[3L]]))
    }
    if (containsForestCall(expr)) {
      refuseInside(expr)
    }
    expr
  }

  hasResponse <- length(formula) == 3L
  lhs <- if (hasResponse) formula[[2L]] else NULL
  rhs <- if (hasResponse) formula[[3L]] else formula[[2L]]
  if (!is.null(lhs) && containsForestCall(lhs)) {
    stop(
      "a forest() term cannot appear on the left-hand side of a formula: '",
      shownCode(lhs),
      "'",
      call. = FALSE
    )
  }
  refuseCrossed(rhs)
  plain <- walk(rhs)
  list(hits = hits, plain = plain)
}

## The terms of a right-hand side as R's own terms() reads them, '.' expanded
## over `data` and '-' applied: the predictor terms' labels, the offset()
## calls, and whether the intercept was removed. `response` is the fit's
## left-hand side, which '.' leaves out. `data` is read for '.' alone.
readRhsTerms <- function(expr, response, data, env) {
  rhs <- if (is.null(response)) call("~", expr) else call("~", response, expr)
  class(rhs) <- "formula"
  environment(rhs) <- env
  read <- if ("." %in% all.names(expr) && !is.null(data)) {
    stats::terms(rhs, data = data)
  } else {
    stats::terms(rhs)
  }
  offsets <- attr(read, "offset")
  list(
    labels = attr(read, "term.labels"),
    offsets = if (length(offsets) > 0L) {
      as.list(attr(read, "variables"))[-1L][offsets]
    },
    noIntercept = attr(read, "intercept") == 0L
  )
}

## Whether a right-hand side writes an intercept term, 1 or 0, at the top of
## its '+' and '-' chain; terms() records 0 and - 1 and not a written 1.
writesInterceptTerm <- function(expr) {
  if (is.numeric(expr)) {
    return(length(expr) == 1L && expr %in% c(0, 1))
  }
  if (
    is.call(expr) &&
      is.name(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% c("+", "-", "(")
  ) {
    return(any(vapply(writtenArguments(expr), writesInterceptTerm, NA)))
  }
  FALSE
}

## One forest() term of a formula as a forest() specification: its
## predictors kept as code with the formula's environment, its basis and
## every other argument evaluated there, as a basis formula written by hand
## already resolves, with the constraint constructors and fixed resolved by
## bare name. `call` is the term as written, for messages.
processHit <- function(hit, env) {
  if (numUnnamed(as.list(hit)[-1L]) > 1L) {
    refuseSecondUnnamed()
  }
  arguments <- forestCallArguments(hit)
  basis <- if ("basis" %in% names(arguments)) {
    eval(arguments[["basis"]], env)
  }
  spec <- do.call(
    forest,
    lapply(
      arguments[names(arguments) %not_in% c("vars", "basis")],
      evalInForestVocabulary,
      vocabulary = forestConstructors[c("interactions", "blocks", "fixed")],
      evalEnv = env
    )
  )
  spec["vars"] <- list(captureForestVars(arguments[["vars"]], env))
  spec["basis"] <- list(basis)
  list(call = hit, spec = spec)
}

## A forest() term's predictors read as the right-hand side of a model
## formula of its own, the labels kept on the specification for the fit's
## formula and for the forest's columns. A value written in the term is names
## of the fit's predictors and is left for the built design to resolve; a
## position is refused, the columns having no order a formula states.
readForestTerms <- function(entry, response, data, env) {
  vars <- entry$spec$vars
  if (is.null(vars)) {
    return(entry)
  }
  shown <- shownCode(entry$call)
  refusePosition <- function() {
    stop(
      "'",
      shown,
      "': in a formula a forest's predictors are written as terms, as ",
      "forest(x1 + x2), or named, as forest(c(\"x1\", \"x2\")); a position ",
      "is for a 'forests' list",
      call. = FALSE
    )
  }
  if (!inherits(vars, "dbartsForestTerms")) {
    if (!is.character(vars)) {
      refusePosition()
    }
    return(entry)
  }
  expr <- vars$expr
  if (length(all.vars(expr)) == 0L) {
    # no name in it: a selection written as a value, c("x1", "x2")
    value <- tryCatch(eval(expr, env), error = function(e) NULL)
    if (!is.character(value)) {
      refusePosition()
    }
    entry$spec$vars <- value
    return(entry)
  }
  if (is.name(expr) && !identical(expr, as.name("."))) {
    name <- as.character(expr)
    inData <- !is.null(data) && name %in% names(data)
    if (!inData && inherits(get0(name, envir = env), "formula")) {
      refuseHeldFormula(expr)
    }
  }
  if (writesInterceptTerm(expr)) {
    stop(
      "'",
      shown,
      "': an intercept term (1, 0 or - 1) is the fit's and no forest's; ",
      "write it beside the forests",
      call. = FALSE
    )
  }
  read <- tryCatch(
    readRhsTerms(expr, response, data, env),
    error = function(e) {
      stop("'", shown, "': ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.null(read$offsets)) {
    stop(
      "'",
      shown,
      "': an offset() is the fit's and no forest's; write it beside the ",
      "forests, as ",
      if (is.null(response)) "" else paste0(shownCode(response), " "),
      "~ ",
      shownCode(read$offsets[[1L]]),
      " + forest(...)",
      call. = FALSE
    )
  }
  if (read$noIntercept) {
    stop(
      "'",
      shown,
      "': an intercept term (1, 0 or - 1) is the fit's and no forest's; ",
      "write it beside the forests",
      call. = FALSE
    )
  }
  if (length(read$labels) == 0L) {
    stop(
      "'",
      shown,
      "': its terms leave the forest no predictor to split on",
      call. = FALSE
    )
  }
  entry$spec$vars$labels <- read$labels
  entry
}

## Phase 1, from the formula and its environment: walk, refuse, read every
## forest's terms, rebuild the fit's own formula from them, and evaluate
## every basis against a model frame built from the SAME data and subset the
## fit itself uses (post-subset - the ambiguity a basis evaluated against
## raw, pre-subset data would otherwise carry). NULL when 'formula' has no
## forest() term, leaving the caller's formula handling untouched. 'family'
## is checked against the multiplier-incompatible set here, at the point each
## entry point has just resolved its own requested token, before any
## family-specific dispatch (a hazard/aft remap, a diversion to bart's own
## multinomial/ordinal/nbinom/hurdle arcs) can make that token unrecoverable
## or unreachable.
##
## The forest with no basis is the forest with no multiplier: the plain
## terms, or the one forest() written without a basis, and forest 1 wherever
## it is written. Where every forest has a basis they keep the order written.
## `forests`, `bases` and `basisTerms` are positional against that order.
ingestFormulaTerms <- function(
  formula,
  family,
  data,
  subsetMissing,
  subsetExpr,
  evalEnv
) {
  if (!is.formula(formula)) {
    return(NULL)
  }
  walked <- walkFormulaTerms(formula)
  if (is.null(walked)) {
    return(NULL)
  }

  if (family %in% TERM_UNSUPPORTED_FAMILIES) {
    stop(
      "family \"",
      family,
      "\" does not support a forest() formula term ('",
      deparse(walked$hits[[1L]]),
      "'): an amplitude-coupled fit is ",
      "not defined for it"
    )
  }

  # a term's own arguments evaluate in the formula's environment, as a
  # hand-written basis = formula already would ("in its own environment");
  # 'subset' below instead resolves in the caller's frame, matching every
  # other model-frame special
  formulaEnv <- environment(formula)
  response <- if (length(formula) == 3L) formula[[2L]] else NULL
  # what '.' expands over, as the fit's own model frame expands it
  termData <- if (is.data.frame(data) || is.list(data)) data else NULL
  written <- lapply(walked$hits, processHit, env = formulaEnv)
  multiplied <- !vapply(written, function(entry) is.null(entry$spec$basis), NA)

  # the plain terms: predictors are the forest with no multiplier; an
  # offset() and an intercept term are the fit's and no forest's
  plain <- if (is.null(walked$plain)) {
    NULL
  } else {
    readRhsTerms(walked$plain, response, termData, formulaEnv)
  }
  hasPlainForest <- length(plain$labels) > 0L
  if (hasPlainForest && !all(multiplied)) {
    stop(
      "the formula has plain terms (",
      paste(plain$labels, collapse = " + "),
      ") and a forest() with no basis ('",
      shownCode(written[[which(!multiplied)[1L]]]$call),
      "'): each is the forest with no multiplier, and a model has one. ",
      "Write the plain terms inside that forest(), or give it a 'basis'",
      call. = FALSE
    )
  }
  if (sum(!multiplied) > 1L) {
    stop(
      "the formula has ",
      sum(!multiplied),
      " forest() terms with no basis; a model has one forest with no ",
      "multiplier, and every other forest states a 'basis'",
      call. = FALSE
    )
  }
  written <- lapply(
    written,
    readForestTerms,
    response = response,
    data = termData,
    env = formulaEnv
  )

  entries <- if (hasPlainForest) {
    first <- forest()
    first$vars <- structure(
      list(expr = walked$plain, env = formulaEnv, labels = plain$labels),
      class = "dbartsForestTerms"
    )
    c(list(list(call = NULL, spec = first)), written)
  } else {
    c(written[!multiplied], written[multiplied])
  }
  forests <- lapply(entries, function(entry) entry$spec)

  # the fit's predictors: every term any forest names, the first forest's
  # first and then the others' in the order written, each once; then the
  # fit's own offset() terms and intercept removal
  labels <- unique(unlist(lapply(forests, function(spec) {
    if (inherits(spec$vars, "dbartsForestTerms")) spec$vars$labels
  })))
  if (length(labels) == 0L) {
    stop(
      "the formula names no predictors: write them as plain terms or inside ",
      "a forest(), as forest(x1 + x2)",
      call. = FALSE
    )
  }
  rhs <- Reduce(
    function(left, right) call("+", left, right),
    c(lapply(labels, str2lang), plain$offsets)
  )
  if (isTRUE(plain$noIntercept)) {
    rhs <- call("-", rhs, 1)
  }
  rewritten <- formula
  rewritten[[length(rewritten)]] <- rhs

  declared <- lapply(forests, function(spec) spec$basis)
  basisVars <- unique(unlist(lapply(declared, function(basis) {
    if (inherits(basis, "formula")) {
      all.vars(basis[[2L]])
    } else {
      NULL
    }
  })))
  basisFrame <- NULL
  if (length(basisVars) > 0L) {
    basisFormula <- stats::reformulate(basisVars)
    environment(basisFormula) <- formulaEnv
    mfCall <- quote(stats::model.frame(
      formula = NULL,
      data = NULL,
      na.action = stats::na.pass,
      drop.unused.levels = FALSE
    ))
    mfCall$formula <- basisFormula
    mfCall$data <- data
    if (!subsetMissing) {
      mfCall$subset <- subsetExpr
    }
    basisFrame <- eval(mfCall, evalEnv)
  }

  evaluated <- lapply(declared, evaluateForestBasis, data = basisFrame)
  bases <- lapply(evaluated, expandForestBasis)

  # what a blend at NEW rows needs to rebuild the same basis, which the
  # expanded matrix alone cannot supply: the declaring formula, the fit-time
  # levels of every factor it reads (so a replay refuses a new level and keeps
  # the order amplitude j is stated against), and the levels of the evaluated
  # value itself when that is categorical, since an expression such as
  # ~ factor(z) derives its own from whatever data it sees and would otherwise
  # set the width from newdata. The call that rebuilds the value from the
  # training rows' centre, scale and knots is stored too. A basis given as a
  # value rather than a formula has no expression to replay and stores none.
  basisTerms <- lapply(seq_along(declared), function(i) {
    basis <- declared[[i]]
    if (!inherits(basis, "formula")) {
      return(NULL)
    }
    factorVars <- Filter(
      is.factor,
      basisFrame[intersect(all.vars(basis[[2L]]), names(basisFrame))]
    )
    value <- evaluated[[i]]
    # an operand is evaluated a second time to find its call; a draw it makes
    # must not move R's generator, which seeds the chains
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed())
    list(
      formula = basis,
      # the expression with what scale(), poly(), ns() and the like computed
      # from the training rows written into the call, as model.frame does for a
      # model formula's terms (stats::makepredictcall); an expression with no
      # such method comes back unchanged and is evaluated on whatever rows
      # predict is given, as lm evaluates it
      predcall = forestBasisPredictCall(
        basis[[2L]],
        basisFrame,
        environment(basis),
        value
      ),
      xlev = if (length(factorVars) > 0L) lapply(factorVars, levels) else NULL,
      levels = if (is.factor(value)) {
        levels(value)
      } else if (is.character(value)) {
        levels(factor(value))
      } else {
        NULL
      }
    )
  })

  list(
    formula = rewritten,
    forests = forests,
    bases = bases,
    basisTerms = basisTerms
  )
}

## The call that rebuilds a basis expression at new rows. model.frame gives
## every variable of a model formula its own stats::makepredictcall; a basis is
## one expression, so the same is done for the operands of the arithmetic
## operators, parentheses and cbind() at its top, and for the call at the top
## when it is none of those. Anything else, I() and indexing among it, comes
## back as written. 'value' is the already evaluated expression, which saves
## evaluating a call at the top again. An operand under the descent is
## evaluated a second time, so a function with a side effect runs twice at fit;
## that is safe for the fit because the caller restores R's generator around
## it and the first evaluation has already raised any warning, so the second
## is silenced.
forestBasisPredictCall <- function(expr, frame, env, value = NULL) {
  if (!is.call(expr)) {
    return(expr)
  }
  if (
    is.name(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% c("+", "-", "*", "/", "^", "(", "cbind")
  ) {
    for (i in seq_along(expr)[-1L]) {
      expr[[i]] <- forestBasisPredictCall(expr[[i]], frame, env)
    }
    return(expr)
  }
  if (is.null(value)) {
    value <- suppressWarnings(eval(expr, frame, env))
  }
  stats::makepredictcall(value, expr)
}

## Phase 2: each forest's predictors become design columns, which exist only
## once dbartsData() has returned. A forest naming every column is the
## unrestricted forest. A formula whose one forest() has no basis and states
## nothing else is the single-forest fit its terms written plainly are, and
## declares no forests.
finalizeTermForests <- function(forests, data) {
  for (index in seq_along(forests)) {
    forests[[index]]["vars"] <- list(resolveForestVars(
      forests[[index]]$vars,
      data,
      allIsNull = TRUE
    ))
  }
  if (length(forests) == 1L && all(vapply(forests[[1L]], is.null, NA))) {
    return(NULL)
  }
  forests
}

## The fitting function's 'n.trees' is the tree count of the forest with no
## basis, so a formula that writes that forest with a count of its own states
## one count twice. Read from the formula as written, nothing evaluated.
refuseTreeCountGivenTwice <- function(formula) {
  walked <- walkFormulaTerms(formula)
  for (hit in walked$hits) {
    given <- tryCatch(names(forestCallArguments(hit)), error = function(e) NULL)
    if ("n.trees" %in% given && "basis" %not_in% given) {
      stop(
        "'n.trees' is given to the fitting function and to the forest with ",
        "no basis ('",
        shownCode(hit),
        "'), which are the same count; give one",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
