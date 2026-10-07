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

## The term grammar names forest: by bare name, or as the constructor is
## written outside the arguments that resolve it, dbartsForests$forest and
## dbarts:::forest, with or without the package named. forest() is not
## exported, so dbarts::forest is no spelling of it.
isForestCall <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  head <- expr[[1L]]
  if (identical(head, as.name("forest"))) {
    return(TRUE)
  }
  if (!is.call(head) || length(head) != 3L) {
    return(FALSE)
  }
  inPackage <- function(owner, name) {
    identical(owner, as.name(name)) ||
      (is.call(owner) &&
        length(owner) == 3L &&
        (identical(owner[[1L]], as.name("::")) ||
          identical(owner[[1L]], as.name(":::"))) &&
        identical(owner[[2L]], as.name("dbarts")) &&
        identical(owner[[3L]], as.name(name)))
  }
  if (identical(head[[1L]], as.name("$"))) {
    return(
      identical(head[[3L]], as.name("forest")) &&
        inPackage(head[[2L]], "dbartsForests")
    )
  }
  identical(head[[1L]], as.name(":::")) && inPackage(head, "forest")
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

## The terms a '+' chain joins, in the order written.
plusTerms <- function(expr) {
  if (isBinaryCall(expr, "+")) {
    return(c(plusTerms(expr[[2L]]), plusTerms(expr[[3L]])))
  }
  list(expr)
}

## The basis a crossed operand stands for: a name, factor() of a name, any
## other call on columns, or a parenthesised sum, which is its members side
## by side, one column each: cbind() of them, a basis formula being evaluated
## as R code. NULL where no one forest() can be written for it: the operand is
## itself crossed or holds a forest(), or the sum has a member cbind() would
## not keep as a column of numbers, a factor() call or a column that
## `columnOf`, where the data are at hand, finds to be neither numeric nor
## logical.
crossedOperandBasis <- function(expr, columnOf = NULL) {
  if (containsForestCall(expr)) {
    return(NULL)
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L
  ) {
    members <- plusTerms(expr[[2L]])
    if (length(members) < 2L) {
      return(crossedOperandBasis(expr[[2L]], columnOf))
    }
    isColumn <- function(member) {
      if (is.name(member)) {
        if (identical(member, as.name("."))) {
          return(FALSE)
        }
        column <- if (!is.null(columnOf)) columnOf(as.character(member))
        return(is.null(column) || is.numeric(column) || is.logical(column))
      }
      is.call(member) &&
        is.name(member[[1L]]) &&
        as.character(member[[1L]]) %not_in%
          c("factor", "as.factor", "ordered", "as.character", "-", ":", "*")
    }
    if (!all(vapply(members, isColumn, NA))) {
      return(NULL)
    }
    return(as.call(c(quote(cbind), members)))
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
refuseCrossedForest <- function(
  expr,
  forestSide = NULL,
  otherSide = NULL,
  columnOf = NULL
) {
  lead <- paste0(
    "'",
    shownCode(expr),
    "': a forest() is not crossed with another term; "
  )
  basis <- if (!is.null(forestSide)) crossedOperandBasis(otherSide, columnOf)
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

## A call of one of `operators` on two operands.
isBinaryCall <- function(expr, operators) {
  is.call(expr) &&
    length(expr) == 3L &&
    is.name(expr[[1L]]) &&
    as.character(expr[[1L]]) %in% operators
}

## The top of a right-hand side with each forest() term replaced by what
## `replace` gives for it, in the order written, or dropped where that is
## NULL. Everything else stays as written, a removal after the forests among
## it. NULL when nothing is left.
replaceForestTerms <- function(expr, replace) {
  if (isForestCall(expr)) {
    return(replace(expr))
  }
  if (isBinaryCall(expr, "+")) {
    left <- replaceForestTerms(expr[[2L]], replace)
    right <- replaceForestTerms(expr[[3L]], replace)
    if (is.null(left)) {
      return(right)
    }
    if (is.null(right)) {
      return(left)
    }
    return(call("+", left, right))
  }
  if (isBinaryCall(expr, "-") && !containsForestCall(expr[[3L]])) {
    left <- replaceForestTerms(expr[[2L]], replace)
    if (is.null(left)) {
      return(call("-", expr[[3L]]))
    }
    return(call("-", left, expr[[3L]]))
  }
  expr
}

## Walks 'formula' for forest() terms. A forest() is a term of the right-hand
## side's top-level '+' chain and nothing else: crossed by ':' or '*', at any
## depth, it is refused with the form to write; on the left-hand side, in a
## removal, or anywhere it could only be reached by evaluating an expression,
## it is refused by name. Returns NULL when no forest() is present anywhere,
## leaving the caller's formula untouched; otherwise the forest() calls in
## the order written, `rhs`, the right-hand side as written, and `plain`, the
## right-hand side without them, NULL when nothing is left. What is written
## at the top beside the forests, a removal among it, stays in `plain`.
## `columnOf` gives a column by its name where the data are at hand, for the
## refusal of a crossed forest.
walkFormulaTerms <- function(formula, columnOf = NULL) {
  if (!containsForestCall(formula)) {
    return(NULL)
  }
  hits <- list()

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
    if (isBinaryCall(expr, c(":", "*")) && containsForestCall(expr)) {
      left <- expr[[2L]]
      right <- expr[[3L]]
      if (isForestCall(right)) {
        refuseCrossedForest(expr, right, left, columnOf)
      }
      if (isForestCall(left)) {
        refuseCrossedForest(expr, left, right, columnOf)
      }
      refuseCrossedForest(expr)
    }
    for (a in writtenArguments(expr)) {
      refuseCrossed(a)
    }
    invisible(NULL)
  }
  # every forest() that is not a term of the top '+' chain
  refuseBuried <- function(expr) {
    if (isForestCall(expr)) {
      return(invisible(NULL))
    }
    if (isBinaryCall(expr, "+")) {
      refuseBuried(expr[[2L]])
      refuseBuried(expr[[3L]])
    } else if (isBinaryCall(expr, "-") && !containsForestCall(expr[[3L]])) {
      refuseBuried(expr[[2L]])
    } else if (containsForestCall(expr)) {
      refuseInside(expr)
    }
    invisible(NULL)
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
  refuseBuried(rhs)
  plain <- replaceForestTerms(rhs, function(hit) {
    hits[[length(hits) + 1L]] <<- hit
    NULL
  })
  list(hits = hits, rhs = rhs, plain = plain)
}

## A right-hand side with what its top removes left out: the terms it
## writes. NULL when it writes none.
withoutRemovals <- function(expr) {
  if (is.null(expr)) {
    return(NULL)
  }
  if (isBinaryCall(expr, "-")) {
    return(withoutRemovals(expr[[2L]]))
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L
  ) {
    return(NULL)
  }
  if (isBinaryCall(expr, "+")) {
    left <- withoutRemovals(expr[[2L]])
    right <- withoutRemovals(expr[[3L]])
    if (is.null(left)) {
      return(right)
    }
    if (is.null(right)) {
      return(left)
    }
    return(call("+", left, right))
  }
  expr
}

## Whether the top of a right-hand side removes a term: a '-' on anything but
## the intercept's 1 or 0.
removesTerms <- function(expr) {
  if (is.null(expr) || !is.call(expr)) {
    return(FALSE)
  }
  isTerm <- function(removed) {
    !(is.numeric(removed) && length(removed) == 1L && removed %in% c(0, 1))
  }
  if (isBinaryCall(expr, "-")) {
    return(isTerm(expr[[3L]]) || removesTerms(expr[[2L]]))
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L
  ) {
    return(isTerm(expr[[2L]]))
  }
  if (isBinaryCall(expr, "+")) {
    return(removesTerms(expr[[2L]]) || removesTerms(expr[[3L]]))
  }
  FALSE
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
## formula and for the forest's columns. A value written in the term is
## names: each that is a column of the data, or with no data a variable where
## the formula was written, is a term as the name written out is, and brings
## its column to the fit; any other is left in `columns` for the built design
## to resolve, a column a term of several gives among them. A position is
## refused, the columns having no order a formula states.
readForestTerms <- function(entry, response, data, env) {
  vars <- entry$spec$vars
  if (is.null(vars)) {
    return(entry)
  }
  shown <- shownCode(entry$call)
  namedTerms <- function(value, expr) {
    # an empty, missing or repeated name is the built design's to refuse
    if (length(value) == 0L || anyNA(value) || anyDuplicated(value) > 0L) {
      return(value)
    }
    isTerm <- if (is.null(data)) {
      vapply(
        value,
        function(name) {
          found <- get0(name, envir = env)
          !is.null(found) && !is.function(found)
        },
        NA
      )
    } else {
      value %in% names(data)
    }
    structure(
      list(
        expr = expr,
        env = env,
        labels = vapply(
          value[isTerm],
          function(name) deparse(as.name(name), backtick = TRUE),
          "",
          USE.NAMES = FALSE
        ),
        columns = value[!isTerm],
        named = TRUE
      ),
      class = "dbartsForestTerms"
    )
  }
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
    entry$spec$vars <- namedTerms(vars, vars)
    return(entry)
  }
  expr <- vars$expr
  if (length(all.vars(expr)) == 0L) {
    # no name in it: a selection written as a value, c("x1", "x2")
    value <- tryCatch(eval(expr, env), error = function(e) NULL)
    if (!is.character(value)) {
      refusePosition()
    }
    entry$spec$vars <- namedTerms(value, expr)
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
## is checked against the multiplier-incompatible set here, for a formula
## with a multiplied forest, at the point each entry point has just resolved
## its own requested token, before any family-specific dispatch (a hazard/aft
## remap, a diversion to bart's own multinomial/ordinal/nbinom/hurdle arcs)
## can make that token unrecoverable or unreachable.
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
  walked <- walkFormulaTerms(formula, function(name) {
    if ((is.data.frame(data) || is.list(data)) && name %in% names(data)) {
      data[[name]]
    } else {
      get0(name, envir = environment(formula))
    }
  })
  if (is.null(walked)) {
    return(NULL)
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
  # a formula whose one forest() has no basis is a single-forest fit, in every
  # family; a multiplier is what these families have no fit for
  if (any(multiplied) && family %in% TERM_UNSUPPORTED_FAMILIES) {
    stop(
      "family \"",
      family,
      "\" does not support a forest() formula term ('",
      deparse(walked$hits[[which(multiplied)[1L]]]),
      "'): an amplitude-coupled fit is ",
      "not defined for it"
    )
  }

  termLabels <- function(vars) {
    if (inherits(vars, "dbartsForestTerms")) vars$labels
  }
  joined <- function(terms) {
    Reduce(function(left, right) call("+", left, right), terms)
  }

  # the plain terms the formula writes, whatever it then removes: predictor
  # terms beside a forest() with no basis are two forests with no multiplier
  written <- lapply(
    written,
    readForestTerms,
    response = response,
    data = termData,
    env = formulaEnv
  )
  plainWritten <- withoutRemovals(walked$plain)
  plainWritten <- if (!is.null(plainWritten)) {
    readRhsTerms(plainWritten, response, termData, formulaEnv)$labels
  }
  if (length(plainWritten) > 0L && !all(multiplied)) {
    stop(
      "the formula has plain terms (",
      paste(plainWritten, collapse = " + "),
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

  # The forest with no basis. Its terms are those of ONE formula, read by
  # terms(): the right-hand side with that forest() replaced, where it
  # stands, by its contents as a group of their own (the call tree keeps
  # them one operand, as parentheses would) and every forest with a basis
  # deleted. With the forest left as plain terms that formula is the plain
  # terms themselves, so the two spellings of a model are one reading and
  # agree whatever is written beside them: a removal before, between or after
  # the forests takes its term from this forest and from no forest with a
  # basis, whose terms are its own. The fit's offset() terms and intercept
  # term are read here too.
  #
  # Written with no first argument that forest is every predictor of the
  # fit, which are the terms of the forests with a basis: those are its
  # contents, so that a removal beside it takes its term from it as from a
  # forest that names them. Its predictors given by value may name columns
  # the fit builds, an indicator of a factor among them, which no formula can
  # be read against before the design exists: beside a removal that is
  # refused, not left with a column the removal was meant to take.
  bare <- which(!multiplied)
  if (length(bare) == 1L) {
    bareVars <- written[[bare]]$spec$vars
    if (is.null(bareVars)) {
      others <- unique(unlist(lapply(written[multiplied], function(entry) {
        termLabels(entry$spec$vars)
      })))
      if (length(others) > 0L) {
        written[[bare]]$spec$vars <- structure(
          list(
            expr = quote(.),
            env = formulaEnv,
            labels = others,
            named = TRUE
          ),
          class = "dbartsForestTerms"
        )
      }
    } else if (
      inherits(bareVars, "dbartsForestTerms") &&
        length(bareVars$columns) > 0L &&
        removesTerms(walked$plain)
    ) {
      stop(
        "'",
        shownCode(written[[bare]]$call),
        "': '",
        bareVars$columns[1L],
        "' is a column the fit builds, not a column of the data, and what ",
        "the formula removes beside the forest cannot be read against it; ",
        "write the forest's predictors as terms, or leave the removal out",
        call. = FALSE
      )
    }
  }
  position <- 0L
  reduced <- replaceForestTerms(walked$rhs, function(hit) {
    position <<- position + 1L
    vars <- written[[position]]$spec$vars
    if (multiplied[position] || length(termLabels(vars)) == 0L) {
      NULL
    } else if (isTRUE(vars$named)) {
      joined(lapply(vars$labels, str2lang))
    } else {
      vars$expr
    }
  })
  first <- if (!is.null(reduced)) {
    readRhsTerms(reduced, response, termData, formulaEnv)
  }
  bareTerms <- if (length(bare) == 1L) termLabels(written[[bare]]$spec$vars)
  if (length(first$labels) == 0L) {
    if (length(bareTerms) > 0L) {
      stop(
        "'",
        shownCode(written[[bare]]$call),
        "': what the formula removes beside it leaves the forest no ",
        "predictor to split on",
        call. = FALSE
      )
    }
    if (length(plainWritten) > 0L) {
      stop(
        "the formula removes every plain term it writes (",
        paste(plainWritten, collapse = " + "),
        "), which leaves the forest with no multiplier no predictor to ",
        "split on",
        call. = FALSE
      )
    }
  }

  entries <- if (length(bare) == 1L) {
    if (length(bareTerms) > 0L) {
      written[[bare]]$spec$vars$labels <- first$labels
    }
    c(written[bare], written[-bare])
  } else if (length(first$labels) > 0L) {
    plainForest <- forest()
    plainForest$vars <- structure(
      list(expr = reduced, env = formulaEnv, labels = first$labels),
      class = "dbartsForestTerms"
    )
    c(list(list(call = NULL, spec = plainForest)), written)
  } else {
    # every forest has a basis: the order written
    written
  }
  forests <- lapply(entries, function(entry) entry$spec)

  # the fit's predictors in the order of the design: the first forest's terms
  # and then the others' in the order written, each once
  labels <- unique(unlist(lapply(forests, function(spec) {
    termLabels(spec$vars)
  })))
  if (length(labels) == 0L) {
    stop(
      "the formula names no predictors: write them as plain terms or inside ",
      "a forest(), as forest(x1 + x2)",
      call. = FALSE
    )
  }
  # the fit's formula: the formula the first forest was read from, as it is
  # written, and after it the terms the forests with a basis add. A fit all
  # of whose forests' terms are among its plain terms so stores the terms of
  # the same formula with no forest() written, and predict asks of new rows
  # what it asks then
  rhs <- joined(c(
    if (!is.null(reduced)) list(reduced),
    lapply(setdiff(labels, first$labels), str2lang)
  ))
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

## A formula whose one forest() term states its predictors, as terms or as
## names, and nothing else, with those predictors written in its place: the
## single-forest formula it is, for a family that reads its formula before
## any forest could be declared. NULL for every other formula, one with a
## second forest() or with any other argument on the one among them.
loneForestFormula <- function(formula) {
  walked <- tryCatch(walkFormulaTerms(formula), error = function(e) NULL)
  if (length(walked$hits) != 1L) {
    return(NULL)
  }
  hit <- walked$hits[[1L]]
  arguments <- tryCatch(forestCallArguments(hit), error = function(e) NULL)
  if (
    numUnnamed(as.list(hit)[-1L]) > 1L ||
      !identical(names(arguments), "vars")
  ) {
    return(NULL)
  }
  contents <- arguments[["vars"]]
  if (length(all.vars(contents)) == 0L) {
    # names by value, each the term the name written out is
    named <- tryCatch(
      eval(contents, environment(formula)),
      error = function(e) NULL
    )
    if (
      !is.character(named) ||
        length(named) == 0L ||
        anyNA(named) ||
        !all(nzchar(named))
    ) {
      return(NULL)
    }
    contents <- Reduce(
      function(left, right) call("+", left, right),
      lapply(named, as.name)
    )
  } else if (writesInterceptTerm(contents)) {
    return(NULL)
  }
  formula[[length(formula)]] <- replaceForestTerms(walked$rhs, function(hit) {
    contents
  })
  formula
}

## The fitting function's 'n.trees' is the tree count of the forest with no
## basis, so a formula that writes that forest with a count of its own states
## one count twice. Read from the formula as written, nothing evaluated.
refuseTreeCountGivenTwice <- function(formula) {
  # what the formula itself refuses is the fitting function's to say, with
  # the data at hand
  walked <- tryCatch(walkFormulaTerms(formula), error = function(e) NULL)
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
