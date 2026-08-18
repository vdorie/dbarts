## A forest() call written inside a formula (bare, or a ':'/'*' operand)
## declares an ADDITIONAL amplitude-coupled forest: sugar over the existing
## forests = list(forest(), forest(basis = )) declaration
## (docs/design/bcf.md), sharing its ingestion between dbarts() and bart2()
## (bart2 reaches it only by forwarding its own formula, unchanged, into
## dbarts()). Everything decidable from the unevaluated formula - grammar,
## ancestor position, the two collisions below - refuses without touching
## data; a basis's degenerate columns and a compound operand's member types
## are the two checks that need real values and happen once the model frame
## exists.

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

isForestCall <- function(expr) {
  is.call(expr) &&
    (identical(expr[[1L]], as.name("forest")) ||
      identical(expr[[1L]], quote(dbarts::forest)))
}

containsForestCall <- function(expr) {
  if (isForestCall(expr)) {
    return(TRUE)
  }
  if (is.call(expr)) {
    for (a in as.list(expr)[-1L]) {
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

## A bare symbol, or a '+' chain of bare symbols, read from the AST without
## evaluating it - the closed grammar the term's unnamed (vars) slot and a
## compound ':' operand's member list both use. '.' parses as a symbol but is
## never an accepted predictor name here.
flattenPlusSymbols <- function(expr) {
  if (is.name(expr)) {
    name <- as.character(expr)
    return(if (name == ".") NULL else name)
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L
  ) {
    left <- flattenPlusSymbols(expr[[2L]])
    right <- flattenPlusSymbols(expr[[3L]])
    if (is.null(left) || is.null(right)) {
      return(NULL)
    }
    return(c(left, right))
  }
  NULL
}

## The desugar target for a ':'/'*' operand's non-forest side: a bare symbol
## or a factor() call pass through as the basis expression itself; a
## '('-wrapped '+' chain of two or more bare symbols becomes a cbind() of
## them. 'members' carries the compound member names back to the caller for
## the type check that must wait until those columns exist (anything else is
## refused here, quoting the operand, since it is not an expression the
## general basis = channel needs a term for at all).
desugarBasisOperand <- function(expr) {
  if (is.name(expr) && as.character(expr) != ".") {
    return(list(expr = expr, members = NULL))
  }
  if (
    is.call(expr) &&
      identical(expr[[1L]], as.name("factor")) &&
      length(expr) == 2L &&
      is.name(expr[[2L]])
  ) {
    return(list(expr = expr, members = NULL))
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L
  ) {
    members <- flattenPlusSymbols(expr[[2L]])
    if (!is.null(members) && length(members) >= 2L) {
      return(list(
        expr = as.call(c(quote(cbind), lapply(members, as.name))),
        members = members
      ))
    }
  }
  if (
    is.call(expr) && identical(expr[[1L]], as.name(":")) && length(expr) == 3L
  ) {
    stop(
      "'",
      deparse(expr),
      "' modulates a forest by more than one operand; ",
      "modulate by one, or combine variables via forest(basis = )"
    )
  }
  stop(
    "'",
    deparse(expr),
    "' is not a supported forest() modulator; use forest(X, basis = ~ ",
    deparse(expr),
    ") directly for anything else"
  )
}

## Walks both sides of 'formula' for forest() terms, refusing three shapes the
## moment they are found: a forest() reachable only by evaluating an
## expression rather than through the additive/':'-operand grammar, a '*'
## with a forest() operand, and a forest() call reachable from both sides of
## one ':'/'*' node (directly or nested inside either operand). A hit is a
## bare top-level forest() call or a ':'/'*' node with one forest() operand;
## legality is
## tracked over the WHOLE ancestor chain back to '~', not just the immediate
## parent, because only +, :, * and ~ preserve the reading that a forest()
## call declares a term rather than being an ordinary expression someone
## just happens to have written; anything else in the chain (I(), ^, a
## removal, any other call, the left-hand side) means the call would only be
## reached by EVALUATING it, which is refused rather than walked past.
## Returns NULL when no forest() term is present anywhere (the walked
## formula is returned byte-for-byte, so a caller with no term is unaffected
## down to object identity), otherwise the rewritten formula (AST surgery:
## each hit's node is deleted from its additive chain, everything else kept)
## and the raw hits, one per term, for the caller's grammar/collision/basis
## processing.
walkFormulaTerms <- function(formula) {
  hits <- list()

  refuseChain <- function(expr) {
    stop(
      "a forest() term must appear as a top-level additive term or a ':' ",
      "operand, not inside '",
      deparse(expr),
      "'"
    )
  }

  walk <- function(expr, legal) {
    if (isForestCall(expr)) {
      if (!legal) {
        refuseChain(expr)
      }
      hits[[length(hits) + 1L]] <<- list(op = NULL, forest = expr, other = NULL)
      return(NULL)
    }
    if (
      is.call(expr) && length(expr) == 3L && identical(expr[[1L]], as.name("+"))
    ) {
      left <- walk(expr[[2L]], legal)
      right <- walk(expr[[3L]], legal)
      if (is.null(left)) {
        return(right)
      }
      if (is.null(right)) {
        return(left)
      }
      if (identical(left, expr[[2L]]) && identical(right, expr[[3L]])) {
        return(expr)
      }
      return(call("+", left, right))
    }
    if (
      is.call(expr) &&
        length(expr) == 3L &&
        (identical(expr[[1L]], as.name(":")) ||
          identical(expr[[1L]], as.name("*")))
    ) {
      isColon <- identical(expr[[1L]], as.name(":"))
      left <- expr[[2L]]
      right <- expr[[3L]]
      leftIsForest <- isForestCall(left)
      rightIsForest <- isForestCall(right)
      if (leftIsForest || rightIsForest) {
        if (!legal) {
          refuseChain(expr)
        }
        forestSide <- if (rightIsForest) right else left
        otherSide <- if (rightIsForest) left else right
        if (!isColon) {
          stop(
            "'",
            deparse(expr),
            "' is not supported: write ",
            deparse(call(":", otherSide, forestSide)),
            " to modulate the forest by ",
            deparse(otherSide),
            ", or ",
            deparse(call("+", otherSide, call(":", otherSide, forestSide))),
            " to also include it as a main effect"
          )
        }
        if (containsForestCall(otherSide)) {
          stop(
            "'",
            deparse(expr),
            "' has a forest() call on both sides of ':';",
            " modulate by one and combine variables via forest(basis = )"
          )
        }
        hits[[length(hits) + 1L]] <<- list(
          op = ":",
          forest = forestSide,
          other = otherSide
        )
        return(NULL)
      }
      if (legal && (containsForestCall(left) || containsForestCall(right))) {
        # neither operand is itself a forest() call, but one holds a nested
        # ':'/'*' chain that buries one (both associativity directions -
        # forest(x1):z:a and a:(b:forest(x1)) - reach here): the same
        # one-operand rule as the flat colon-chain case above applies, so it
        # is named the same way rather than left to fall through unstripped
        stop(
          "'",
          deparse(expr),
          "' modulates a forest by more than one operand; ",
          "modulate by one, or combine variables via forest(basis = )"
        )
      }
      if (!legal) {
        if (containsForestCall(left)) {
          refuseChain(left)
        }
        if (containsForestCall(right)) refuseChain(right)
      }
      return(expr)
    }
    if (is.call(expr)) {
      for (a in as.list(expr)[-1L]) {
        if (containsForestCall(a)) walk(a, FALSE)
      }
      return(expr)
    }
    expr
  }

  hasResponse <- length(formula) == 3L
  lhs <- if (hasResponse) formula[[2L]] else NULL
  rhs <- if (hasResponse) formula[[3L]] else formula[[2L]]
  if (!is.null(lhs) && containsForestCall(lhs)) {
    stop(
      "a forest() term cannot appear on the left-hand side of a formula: '",
      deparse(lhs),
      "'"
    )
  }

  rewrittenRhs <- walk(rhs, TRUE)
  if (length(hits) == 0L) {
    return(NULL)
  }
  if (is.null(rewrittenRhs)) {
    rewrittenRhs <- 1
  }
  if (length(all.vars(rewrittenRhs)) == 0L) {
    stop(
      "a formula whose only right-hand-side content is a forest() term has ",
      "no predictors left for the base forest after rewriting"
    )
  }

  rewritten <- formula
  rewritten[[length(rewritten)]] <- rewrittenRhs
  list(hits = hits, formula = rewritten)
}

## Parses one raw hit into knob arguments and an effective basis value,
## refusing an unsupported symbolic-vars expression, an unsupported left-
## operand expression, and the two ways one knob can be named twice on one
## term - sugar together with an explicit basis =, or the unnamed slot
## together with an explicit vars = - wherever the closed unnamed-slot or
## left-operand grammar is violated. Every named
## argument evaluates in 'env' - the formula's own environment, since a term
## lives inside the formula literal, matching how a hand-written basis =
## formula already resolves "in its own environment" - except vars, which
## stays symbolic (a character vector of names) when it came from the
## unnamed slot, deferred to resolveTermVars once the design exists.
processHit <- function(hit, env) {
  args <- as.list(hit$forest)[-1L]
  argNames <- names(args)
  if (is.null(argNames)) {
    argNames <- rep("", length(args))
  }
  unnamedIndex <- which(argNames == "")
  if (length(unnamedIndex) > 1L) {
    stop(
      "'",
      deparse(hit$forest),
      "' gives more than one unnamed argument; a ",
      "forest() term's unnamed slot is vars, symbolically, and takes at ",
      "most one"
    )
  }
  named <- args[argNames != ""]
  namedNames <- argNames[argNames != ""]

  varsSymbols <- NULL
  if (length(unnamedIndex) == 1L) {
    if ("vars" %in% namedNames) {
      stop(
        "'",
        deparse(hit$forest),
        "' gives vars both positionally and by ",
        "name; the unnamed slot IS vars here - give one"
      )
    }
    varsSymbols <- flattenPlusSymbols(args[[unnamedIndex]])
    if (is.null(varsSymbols)) {
      stop(
        "'",
        deparse(args[[unnamedIndex]]),
        "' is not a supported forest() vars expression; use vars = for ",
        "anything else"
      )
    }
  }

  compoundMembers <- NULL
  basisValue <- NULL
  if (!is.null(hit$op)) {
    if ("basis" %in% namedNames) {
      stop(
        "'",
        deparse(hit$forest),
        "' declares a basis both through '",
        hit$op,
        "' and 'basis ='; give one"
      )
    }
    desugared <- desugarBasisOperand(hit$other)
    compoundMembers <- desugared$members
    basisValue <- eval(call("~", desugared$expr), env)
  } else if ("basis" %in% namedNames) {
    basisValue <- eval(named[[which(namedNames == "basis")]], env)
  }

  list(
    knobArgs = lapply(
      named[namedNames %not_in% c("vars", "basis")],
      eval,
      envir = env
    ),
    varsArg = if ("vars" %in% namedNames) {
      eval(named[[which(namedNames == "vars")]], env)
    } else {
      NULL
    },
    varsSymbols = varsSymbols,
    basisValue = basisValue,
    compoundMembers = compoundMembers,
    other = hit$other
  )
}

## Phase 1, AST- and formula-environment-only: walk, refuse, desugar, and
## evaluate every basis against a model frame built from the SAME data and
## subset the fit itself uses (post-subset, post-NA - the ambiguity a basis
## evaluated against raw, pre-subset data would otherwise carry). NULL when
## 'formula' has no forest() term, leaving the caller's existing formula
## handling byte-for-byte untouched. 'family' is checked against the
## multiplier-incompatible set here, at the point each entry point has just
## resolved its own requested token, before any family-specific dispatch
## (a hazard/aft remap, a diversion to bart2's own multinomial/ordinal/
## nbinom/hurdle arcs) can make that token unrecoverable or unreachable.
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
      deparse(walked$hits[[1L]]$forest),
      "'): an amplitude-coupled fit is ",
      "not defined for it"
    )
  }

  # a term's own arguments evaluate in the formula's environment, as a
  # hand-written basis = formula already would ("in its own environment");
  # 'subset' below instead resolves in the caller's frame, matching every
  # other model-frame special
  formulaEnv <- environment(formula)
  pending <- lapply(walked$hits, processHit, env = formulaEnv)

  basisVars <- unique(unlist(lapply(pending, function(p) {
    if (inherits(p$basisValue, "formula")) {
      all.vars(p$basisValue[[2L]])
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

  for (p in pending) {
    for (member in p$compoundMembers) {
      column <- basisFrame[[member]]
      if (!is.numeric(column) && !is.logical(column)) {
        stop(
          "a compound modulating operand takes numeric or logical ",
          "variables only; '",
          member,
          "' in '",
          deparse(p$other),
          "' is ",
          "not - give it its own term, or use forest(basis = ) with an ",
          "explicit basis"
        )
      }
    }
  }

  bases <- lapply(pending, function(p) {
    if (is.null(p$basisValue)) {
      return(NULL)
    }
    expandForestBasis(evaluateForestBasis(p$basisValue, basisFrame))
  })

  list(formula = walked$formula, pending = pending, bases = bases)
}

## Phase 2: resolve each term's vars against the built design (needs
## data@x, so this runs once dbartsData() has returned) and construct the
## final forests = declaration - the fit's own (unrestricted) forest 1, plus
## one forest() per term in term order. A symbolic name expands to every
## indicator column derived from it (resolveTermColumns' term-label
## expansion, the variance selector's own precedent), unlike the bare
## match() an explicit vars = performs; a name resolving to nothing is
## refused, quoting it and the rewritten right-hand side's own terms.
## Unrecognized knob names refuse by construction, through forest() itself.
finalizeTermForests <- function(pending, data) {
  colNames <- colnames(data@x)
  termLabels <- attr(data@x, "term.labels")
  forests <- vector("list", length(pending) + 1L)
  forests[[1L]] <- dbarts::forest()
  for (i in seq_along(pending)) {
    p <- pending[[i]]
    vars <- p$varsArg
    if (!is.null(p$varsSymbols)) {
      resolved <- character(0)
      for (name in p$varsSymbols) {
        columns <- resolveTermColumns(name, colNames, termLabels)
        if (is.null(columns)) {
          stop(
            "'",
            name,
            "' is not a term of the rewritten formula's right-",
            "hand side (",
            paste(termLabels, collapse = ", "),
            "); a ",
            "forest() term's vars names a subset of it"
          )
        }
        resolved <- c(resolved, colNames[columns])
      }
      vars <- unique(resolved)
    }
    args <- p$knobArgs
    args$vars <- vars
    forests[[i + 1L]] <- do.call(dbarts::forest, args)
  }
  forests
}
