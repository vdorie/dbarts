## A forest's basis written as code: the right-hand side of a model formula
## with no tilde, read as lm() reads one. '+' separates columns, I() holds
## arithmetic, a factor gives a column for each level, and scale() and poly()
## are rebuilt at new rows from the fitted rows. Both doors, a forest() term
## of a formula and a forest() in a 'forests' list, reduce a basis to the same
## code with an environment and build it with the same functions, so that one
## text is one model at either.

## A call at the top of a basis term that would later state a prior on the
## term is refused now, so that such a statement can be added without changing
## what any accepted call means.
BASIS_RESERVED_CALLS <- c(
  "normal",
  "fixed",
  "student",
  "cauchy",
  "linear",
  "gp",
  "cgm",
  "dart",
  "chisq",
  "chi",
  "invchi",
  "forest",
  "varianceForest"
)

## The operators whose operands are arithmetic and not terms of a model
## formula: a member of cbind() that is one of these calls is written inside
## I() when the refusal of cbind() says what to write.
BASIS_ARITHMETIC_CALLS <- c(
  "+",
  "-",
  "*",
  "/",
  "^",
  "%%",
  "%/%",
  "%*%",
  "%o%",
  "%in%",
  "<",
  ">",
  "<=",
  ">=",
  "==",
  "!=",
  "!",
  "&",
  "|",
  "&&",
  "||",
  ":",
  "("
)

## Checks the top of a basis against its grammar before anything is
## evaluated: '+' separates terms, parentheses group, ':' between terms is
## their product as in a model formula, 1 asks for a constant column and 0 or
## '- 1' for none. Every other operator a model formula gives a meaning of its
## own, or that a size per column could later be written with, is refused by
## name with the form to write; anything else at the top is a term, evaluated
## as written.
parseBasisGrammar <- function(expr) {
  isNumber <- function(e) is.numeric(e) && length(e) == 1L
  refuseOperator <- function(e, operator, means) {
    stop(
      "'basis' does not take '",
      operator,
      "' between its terms ('",
      shownCode(e),
      "'): in a model formula it ",
      means,
      "; write the arithmetic inside I(), as I(",
      shownCode(e),
      ")",
      call. = FALSE
    )
  }
  walk <- function(e, inProduct = FALSE) {
    if (isNumber(e)) {
      if (inProduct || (e != 0 && e != 1)) {
        stop(
          "'basis' has the number ",
          shownCode(e),
          " as a term; only 1, a constant column, and 0, none, are terms. ",
          "Write arithmetic on a column inside I()",
          call. = FALSE
        )
      }
      return(invisible(NULL))
    }
    if (is.name(e)) {
      if (identical(e, as.name("."))) {
        stop(
          "'basis' does not take '.': name the columns the forest is ",
          "multiplied by",
          call. = FALSE
        )
      }
      return(invisible(NULL))
    }
    if (!is.call(e)) {
      # a logical, a string, NULL or any other constant written as a term
      stop(
        "'basis' has the constant ",
        shownCode(e),
        " as a term; a term is a column, and only 1, a constant column, and ",
        "0, none, are written as numbers",
        call. = FALSE
      )
    }
    if (!is.name(e[[1L]])) {
      return(invisible(NULL))
    }
    operator <- as.character(e[[1L]])
    arity <- length(e) - 1L
    if (operator == "~") {
      stop("a 'basis' formula must be one-sided, as ~ dose", call. = FALSE)
    }
    if (operator %in% c("(", "+") && arity == 1L) {
      return(walk(e[[2L]], inProduct))
    }
    if (operator == "+" && arity == 2L) {
      walk(e[[2L]], inProduct)
      return(walk(e[[3L]], inProduct))
    }
    if (operator %in% c("/", "*") && arity == 2L) {
      scaled <- (operator == "*" && isNumber(e[[2L]])) || isNumber(e[[3L]])
      if (scaled) {
        stop(
          "'basis' term '",
          shownCode(e),
          "' ",
          if (operator == "/") "divides" else "multiplies",
          " a column by a number, which a model formula does not take; ",
          "write I(",
          shownCode(e),
          ") for the rescaled column",
          call. = FALSE
        )
      }
    }
    if (operator == "*" && arity == 2L) {
      stop(
        "'basis' does not take '*' between its terms ('",
        shownCode(e),
        "'): in a model formula it is both columns and their product. Write ",
        "I(",
        shownCode(e),
        ") for the product alone, or ",
        shownCode(e[[2L]]),
        " + ",
        shownCode(e[[3L]]),
        " + I(",
        shownCode(e),
        ") for all three",
        call. = FALSE
      )
    }
    if (operator == ":" && arity == 2L) {
      walk(e[[2L]], TRUE)
      return(walk(e[[3L]], TRUE))
    }
    if (operator == "-") {
      removed <- e[[arity + 1L]]
      if (isNumber(removed) && removed == 1 && !inProduct) {
        if (arity == 2L) {
          walk(e[[2L]], inProduct)
        }
        return(invisible(NULL))
      }
      refuseOperator(e, "-", "removes a term and subtracts nothing")
    }
    if (operator == "/") {
      refuseOperator(e, "/", "nests one term in another and divides nothing")
    }
    if (operator == "^") {
      refuseOperator(e, "^", "crosses terms and raises nothing to a power")
    }
    if (operator == "%in%") {
      stop(
        "'basis' does not take '%in%' ('",
        shownCode(e),
        "'): in a model formula it nests one term in another",
        call. = FALSE
      )
    }
    if (operator %in% c("|", "||")) {
      stop(
        "'basis' does not take '",
        operator,
        "' ('",
        shownCode(e),
        "')",
        call. = FALSE
      )
    }
    if (operator == "offset") {
      stop("'basis' does not take an offset() term", call. = FALSE)
    }
    if (operator == "cbind") {
      refuseBoundColumns(e)
    }
    if (operator %in% BASIS_RESERVED_CALLS) {
      stop(
        "'basis' term '",
        shownCode(e),
        "' calls ",
        operator,
        "(), which is not a column: a prior is not stated on a basis term. ",
        "Give the forest's size as 'sd' and its coefficient's law as ",
        "'amplitude'",
        call. = FALSE
      )
    }
    invisible(NULL)
  }
  walk(expr)
}

## The refusal of cbind() at the top of a basis, with the '+' form to write
## where there is one that names the same columns: each member as a term,
## arithmetic inside I(). The form is given only when the grammar takes it.
refuseBoundColumns <- function(expr) {
  members <- as.list(expr)[-1L]
  asTerm <- function(member) {
    if (
      is.name(member) ||
        (is.numeric(member) && length(member) == 1L && member == 1)
    ) {
      return(shownCode(member))
    }
    if (!is.call(member)) {
      return(NA_character_)
    }
    arithmetic <- is.name(member[[1L]]) &&
      as.character(member[[1L]]) %in% BASIS_ARITHMETIC_CALLS
    if (arithmetic) paste0("I(", shownCode(member), ")") else shownCode(member)
  }
  written <- vapply(members, asTerm, "")
  # a model formula keeps one of two terms alike and puts its constant column
  # first, so neither has a '+' form with cbind()'s columns in cbind()'s order
  sameColumns <- length(written) > 0L &&
    !anyNA(written) &&
    anyDuplicated(written) == 0L &&
    all(which(written == "1") == 1L)
  rewrite <- if (sameColumns) {
    text <- paste(written, collapse = " + ")
    taken <- tryCatch(
      {
        parseBasisGrammar(str2lang(text))
        TRUE
      },
      error = function(e) FALSE
    )
    if (taken) text
  }
  stop(
    "'basis' term '",
    shownCode(expr),
    "': the columns of a basis are separated by '+'",
    if (!is.null(rewrite)) paste0("; write ", rewrite),
    call. = FALSE
  )
}

## A basis as code: the right-hand side `expr`, any tilde taken off, with the
## environment its names are looked up in after the data's columns, and its
## text, which is the forest's label when no other is given. A formula brings
## its own environment. The grammar is checked here, and a basis that could
## have no column or the constant column alone is refused here, nothing being
## evaluated for either. As it leaves here the code is read when the model is
## built; `atCall` marks code that captureForestBasis() bound to a call of
## forest() instead.
basisCode <- function(expr, env) {
  isFormula <- is.call(expr) && identical(expr[[1L]], as.name("~"))
  if (isFormula) {
    if (length(expr) != 2L) {
      stop("a 'basis' formula must be one-sided, as ~ dose", call. = FALSE)
    }
    if (!is.null(environment(expr))) {
      env <- environment(expr)
    }
    expr <- expr[[2L]]
  }
  label <- shownCode(expr)
  parseBasisGrammar(expr)
  read <- tryCatch(
    stats::terms(basisFormula(expr, env)),
    error = function(e) {
      stop("'basis' (", label, "): ", conditionMessage(e), call. = FALSE)
    }
  )
  if (length(attr(read, "term.labels")) == 0L) {
    if (attr(read, "intercept") != 0L) {
      stop(
        "'basis' (",
        label,
        ") is a constant column and nothing else, which multiplies the ",
        "forest by a constant: leave 'basis' out for the forest with no ",
        "multiplier",
        call. = FALSE
      )
    }
    stop(
      "'basis' (",
      label,
      ") names no column; leave 'basis' out for the forest with no ",
      "multiplier",
      call. = FALSE
    )
  }
  structure(
    list(expr = expr, env = env, label = label, formula = isFormula),
    class = "dbartsForestBasis"
  )
}

## The model formula a basis is: no intercept unless it writes one, and the
## basis as one operand, so that it is grouped as parentheses would group it.
basisFormula <- function(expr, env) {
  formula <- call("~", call("+", 0, expr))
  class(formula) <- "formula"
  environment(formula) <- env
  formula
}

## Whether `name` is a column of the data a fit was given, a data frame, a
## list or an environment; with no such data nothing is a column.
isDataColumn <- function(name, data) {
  if (is.environment(data)) {
    return(vapply(name, exists, NA, envir = data, inherits = FALSE))
  }
  if (is.list(data)) {
    return(name %in% names(data))
  }
  rep(FALSE, length(name))
}

## What a basis's code gives against `data` and then `env`: its model frame
## over every row, with R's own record of how each term is rebuilt at new
## rows, as list(frame = ); or, where the code holds a one-sided formula or
## NULL and not columns, list(held = ) of that. A name that is no column and
## holds a formula stands for the formula, and code that gives NULL states no
## basis; code that is one call is evaluated a second time, plainly and
## quietly, only when no frame could be built from it, to see whether it gives
## one of the two.
evaluateBasisCode <- function(expr, env, data = NULL) {
  heldBy <- function(value) {
    if (is.null(value) || inherits(value, "formula")) list(held = value)
  }
  if (is.name(expr) && !isDataColumn(as.character(expr), data)) {
    found <- tryCatch(eval(expr, env), error = function(e) e)
    if (!inherits(found, "error") && !is.null(held <- heldBy(found))) {
      return(held)
    }
  }
  frameCall <- quote(stats::model.frame(
    formula = NULL,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  ))
  frameCall$formula <- basisFormula(expr, env)
  if (!is.null(data)) {
    frameCall$data <- data
  }
  frame <- tryCatch(eval(frameCall), error = function(e) e)
  if (!inherits(frame, "error")) {
    return(list(frame = frame))
  }
  isOneCall <- is.call(expr) &&
    !(is.name(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% c("+", "-", ":", "("))
  if (isOneCall) {
    again <- tryCatch(
      suppressWarnings(
        if (is.null(data)) eval(expr, env) else eval(expr, data, env)
      ),
      error = function(e) e
    )
    if (!inherits(again, "error") && !is.null(held <- heldBy(again))) {
      return(held)
    }
  }
  stop(frame)
}

## The names a basis's code reads as variables. What stands to the right of
## '$' or '@' is a component's name and no variable, and a function written
## in place has variables of its own.
basisVariables <- function(expr) {
  if (is.name(expr)) {
    name <- as.character(expr)
    return(if (nzchar(name)) name else character())
  }
  if (!is.call(expr)) {
    return(character())
  }
  head <- expr[[1L]]
  parts <- as.list(expr)[-1L]
  if (is.name(head)) {
    operator <- as.character(head)
    if (operator %in% c("$", "@")) {
      parts <- parts[1L]
    } else if (operator %in% c("::", ":::", "function")) {
      parts <- list()
    }
  } else {
    parts <- c(list(head), parts)
  }
  unique(unlist(lapply(parts, basisVariables), use.names = FALSE))
}

## What a caller binds where a basis is written in a call of forest(): every
## name the code uses that is bound in the frames from the call's own up to
## the workspace or the package the call was made in, copied into `env`,
## which then stands for the place of the call, so that what the caller's
## variables hold later changes nothing, a number written beside a column, a
## function of the caller's and a loop's own variable included. Names found
## further out, a package's functions and attached data, are left to be
## looked up. `unbound` is the variables nothing binds at the call and
## `unboundFunctions` the functions, none of which is looked up later; a name
## whose value cannot be had there keeps the reason in `broken`. Nothing is
## signalled from here: what taking a value warns of is returned in
## `warnings`.
bindBasisAtCall <- function(expr, place) {
  top <- topenv(place)
  bound <- new.env(parent = top)
  broken <- list()
  frameOf <- function(name) {
    frame <- place
    repeat {
      if (exists(name, envir = frame, inherits = FALSE)) {
        return(frame)
      }
      if (identical(frame, top) || identical(frame, emptyenv())) {
        return(NULL)
      }
      frame <- parent.env(frame)
    }
  }
  used <- all.names(expr, unique = TRUE)
  warned <- list()
  for (name in used) {
    frame <- frameOf(name)
    if (
      is.null(frame) ||
        (identical(frame, top) && !identical(top, globalenv()))
    ) {
      next
    }
    # a name may be an argument not yet evaluated, which this evaluates
    value <- tryCatch(
      withCallingHandlers(
        list(get(name, envir = frame, inherits = FALSE)),
        warning = function(w) {
          warned[[length(warned) + 1L]] <<- w
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    if (inherits(value, "error")) {
      broken[[name]] <- conditionMessage(value)
    } else {
      assign(name, value[[1L]], envir = bound)
    }
  }
  variables <- basisVariables(expr)
  functions <- setdiff(used, all.vars(expr))
  list(
    env = bound,
    unbound = variables[!vapply(variables, exists, NA, envir = place)],
    unboundFunctions = functions[
      !vapply(functions, exists, NA, envir = place, mode = "function")
    ],
    broken = broken,
    warnings = warned
  )
}

## forest()'s 'basis' as written, by the caller in a call of forest(). A
## value (a column, a factor, a matrix, NULL) is kept as it is. Code is kept
## unevaluated: a formula with its environment, to be read when the model is
## built as lm() reads one, and any other code bound to the call
## (bindBasisAtCall) with the value it has there taken once (takeAtCall),
## what it warns of while it is bound kept with that value.
captureForestBasis <- function(expr, env) {
  if (!is.language(expr)) {
    return(expr)
  }
  code <- basisCode(expr, env)
  if (code$formula) {
    return(code)
  }
  binding <- bindBasisAtCall(code$expr, env)
  code <- takeAtCall(
    code,
    function() evaluateBasisCode(code$expr, binding$env),
    env
  )
  if (isTRUE(code$evaluated)) {
    code$warnings <- c(binding$warnings, code$warnings)
  }
  code$atCall <- TRUE
  code$env <- binding$env
  code$unbound <- binding$unbound
  code$unboundFunctions <- binding$unboundFunctions
  code$broken <- binding$broken
  code
}

## A single string or number is no basis: a basis has a value for every
## observation. The string is the commonest slip of a caller who builds
## forests by program, so its text gives the form to write. `holder` is the
## code that holds the value, where it was not written in place.
refuseSingleBasisValue <- function(value, holder = NULL) {
  if (!is.atomic(value) || length(value) != 1L || !is.null(dim(value))) {
    return(invisible(NULL))
  }
  what <- if (is.null(holder)) {
    "is"
  } else {
    paste0("is '", holder, "', which holds")
  }
  if (is.character(value)) {
    stop(
      "'basis' ",
      what,
      " the string \"",
      value,
      "\": a basis is a column, not its name. Write basis = ",
      value,
      ", or from a program do.call(forest, list(basis = as.name(\"",
      value,
      "\"))), or hand over the column itself",
      call. = FALSE
    )
  }
  stop(
    "'basis' ",
    what,
    " the single value ",
    format(value),
    ": a basis has a value for every observation",
    call. = FALSE
  )
}

## A declared basis read for a fit of `numRows` observations, every row of
## the data the fit was given. NULL when the forest has none; list(value = )
## for a value handed over; and for code list(frame = , label = ,
## columns = ), its model frame over those rows, its text and the columns of
## the data it names.
##
## Code a formula holds, a term of the fit's formula included, is read
## against the data and then the formula's environment, when the model is
## built. Code written in a call of forest() is read against the data when it
## names a column of it, a column hiding a variable of the caller's with its
## name, and then against what the call bound, nothing being looked up again;
## code that names no column is the value taken at the call, whose warnings
## are raised here, the value being used. `data` is NULL where the fit was
## given no data frame, list or environment.
readForestBasis <- function(basis, data, numRows) {
  if (is.null(basis)) {
    return(NULL)
  }
  if (!inherits(basis, "dbartsForestBasis")) {
    refuseSingleBasisValue(basis)
    return(list(value = basis))
  }
  code <- basis
  refuse <- function(...) {
    stop("'basis' (", code$label, "): ", ..., call. = FALSE)
  }
  where <- if (is.null(data)) {
    " where forest() was called"
  } else {
    " in 'data' or where forest() was called"
  }
  # R's own message for a name found nowhere says where it was looked for
  refuseEvaluation <- function(message) {
    for (name in basisVariables(code$expr)) {
      notFound <- gettextf("object '%s' not found", name, domain = "R")
      if (identical(message, notFound)) {
        refuse("object '", name, "' not found", where)
      }
    }
    refuse(message)
  }
  # a formula may be held by a name a formula holds
  read <- NULL
  for (depth in seq_len(8L)) {
    variables <- basisVariables(code$expr)
    atColumns <- variables[isDataColumn(variables, data)]
    if (!isTRUE(code$atCall) || length(atColumns) > 0L) {
      # what nothing bound at the call is not looked up later than the call
      late <- setdiff(code$unbound, atColumns)
      late <- late[vapply(late, exists, NA, envir = code$env)]
      if (length(late) > 0L) {
        refuse("object '", late[1L], "' not found", where)
      }
      late <- code$unboundFunctions
      late <- late[vapply(
        late,
        exists,
        NA,
        envir = code$env,
        mode = "function"
      )]
      if (length(late) > 0L) {
        refuse(
          "could not find function \"",
          late[1L],
          "\" where forest() was called"
        )
      }
      for (name in setdiff(names(code$broken), atColumns)) {
        refuse(code$broken[[name]])
      }
      read <- tryCatch(
        evaluateBasisCode(code$expr, code$env, data),
        error = function(e) refuseEvaluation(conditionMessage(e))
      )
    } else if (isTRUE(code$evaluated)) {
      for (kept in code$warnings) {
        warning(kept)
      }
      read <- code$value
    } else {
      # an argument that could not be evaluated says why before the code that
      # used it does
      for (name in names(code$broken)) {
        refuse(code$broken[[name]])
      }
      refuseEvaluation(code$error)
    }
    if (!is.null(read$frame)) {
      break
    }
    if (is.null(read$held)) {
      return(NULL)
    }
    code <- basisCode(read$held, environment(read$held))
    read <- NULL
  }
  if (is.null(read)) {
    refuse("a formula holds a formula more deeply than is read")
  }
  frame <- read$frame
  if (nrow(frame) != numRows) {
    if (length(frame) == 1L && nrow(frame) == 1L) {
      refuseSingleBasisValue(frame[[1L]], code$label)
    }
    stop(
      "'basis' (",
      code$label,
      ") must have the same length as the data: it has ",
      nrow(frame),
      " rows and the data ",
      numRows,
      "; a basis covers every row of the data and is cut by 'subset' and ",
      "the na.action with it",
      call. = FALSE
    )
  }
  list(frame = frame, label = code$label, columns = atColumns)
}

## The bases a 'forests' list declares, read for a fit (readForestBasis), one
## element per forest and NULL where a forest has none. A forest whose code
## gives NULL has no basis and is returned saying so, for the checks made of
## the list.
readDeclaredBases <- function(forests, data, numRows) {
  reads <- vector("list", length(forests))
  for (index in seq_along(forests)) {
    read <- readForestBasis(forests[[index]]$basis, data, numRows)
    if (is.null(read)) {
      forests[[index]]["basis"] <- list(NULL)
    } else {
      reads[[index]] <- read
    }
  }
  list(forests = forests, reads = reads)
}

## The column names of a basis handed over as a value: kept when every column
## has one and they differ; none when any column lacks one, the columns then
## being known by position; two alike refused.
valueBasisNames <- function(names) {
  if (is.null(names) || anyNA(names) || !all(nzchar(names))) {
    return(NULL)
  }
  refuseSharedBasisNames(names)
  names
}

## The columns of a basis are told apart by name, so no two share one.
refuseSharedBasisNames <- function(names, label = NULL) {
  if (anyDuplicated(names) > 0L) {
    stop(
      "'basis' ",
      if (!is.null(label)) paste0("(", label, ") "),
      "has two columns named \"",
      names[duplicated(names)][[1L]],
      "\"; the columns of a basis are told apart by name",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## A value handed over as a basis, expanded to its columns
## (expandForestBasis) and named by the rule for a value.
expandValueBasis <- function(value, ...) {
  expanded <- expandForestBasis(value, ...)
  if (is.null(expanded)) {
    return(NULL)
  }
  given <- valueBasisNames(colnames(expanded))
  dimnames(expanded) <- if (is.null(given)) NULL else list(NULL, given)
  expanded
}

## The columns of a basis from the model frame of its code, on the rows the
## frame holds. A basis is one factor, a character vector or a logical
## vector, by class, with a column for each level that has a row, built by
## expandForestBasis and by no contrast; or numeric terms, with the columns
## model.matrix() gives them. The columns are named as model.matrix() names
## them. With `atPrediction` the rows are new rows of a fit whose basis had
## the levels `levels`, NULL for a numeric one: the same columns are built, a
## level the new rows lack included, and a level the fit did not have and a
## basis of the other kind are refused. Returns list(basis = , levels = ).
basisColumns <- function(
  terms,
  frame,
  label,
  levels = NULL,
  atPrediction = FALSE
) {
  if (any(vapply(frame, anyNA, NA))) {
    stop("a 'basis' cannot be NA")
  }
  termLabels <- attr(terms, "term.labels")
  isLevels <- vapply(
    frame,
    function(variable) {
      is.factor(variable) ||
        is.character(variable) ||
        (is.logical(variable) && is.null(dim(variable)))
    },
    NA
  )
  if (any(isLevels)) {
    if (
      length(frame) != 1L ||
        length(termLabels) != 1L ||
        attr(terms, "intercept") != 0L
    ) {
      stop(
        "'basis' (",
        label,
        ") mixes a factor with other terms: a basis is one factor, a ",
        "character or a logical vector, with one coefficient per level, or ",
        "numeric columns, with one each. Give the other terms a forest() of ",
        "their own, or write the factor as numbers",
        call. = FALSE
      )
    }
    value <- frame[[1L]]
    if (atPrediction && is.null(levels)) {
      stop(
        "'basis' (",
        label,
        ") is a factor, a character or a logical vector at the new rows and ",
        "was numeric at the fit",
        call. = FALSE
      )
    }
    if (!atPrediction) {
      if (is.logical(value)) {
        if (all(value) || !any(value)) {
          stop(
            "'basis' (",
            label,
            ") is ",
            if (all(value)) "TRUE" else "FALSE",
            " on every row the fit keeps, so its other level has no ",
            "observations and the forest would be multiplied by a constant",
            call. = FALSE
          )
        }
        levels <- c("FALSE", "TRUE")
      } else {
        # a level the kept rows leave empty is no column, whatever emptied it
        asFactor <- if (is.factor(value)) value else factor(value)
        levels <- levels(asFactor)[
          tabulate(as.integer(asFactor), nlevels(asFactor)) > 0L
        ]
      }
    }
    leveled <- factor(as.character(value), levels = levels)
    if (anyNA(leveled)) {
      stop(
        "'basis' (",
        label,
        ") has the level '",
        as.character(value)[is.na(leveled)][1L],
        "' at a new row, which no row of the fit had",
        call. = FALSE
      )
    }
    basis <- expandForestBasis(leveled, atPrediction = atPrediction)
    colnames(basis) <- paste0(termLabels, levels)
    return(list(basis = basis, levels = levels))
  }
  if (atPrediction && !is.null(levels)) {
    stop(
      "'basis' (",
      label,
      ") is numeric at the new rows and had levels at the fit",
      call. = FALSE
    )
  }
  for (variable in frame) {
    if (!is.numeric(variable)) {
      kind <- if (is.logical(variable)) {
        "a logical matrix"
      } else {
        paste0(
          if (grepl("^[aeiouAEIOU]", class(variable)[1L])) "an " else "a ",
          class(variable)[1L]
        )
      }
      stop(
        "'basis' (",
        label,
        ") has a term that is ",
        kind,
        ": a basis is a factor, a character or logical vector, or numeric",
        call. = FALSE
      )
    }
  }
  design <- stats::model.matrix(terms, frame)
  given <- colnames(design)
  refuseSharedBasisNames(given, label)
  basis <- matrix(as.double(design), nrow(design), ncol(design))
  basis <- expandForestBasis(basis, atPrediction = atPrediction)
  colnames(basis) <- given
  list(basis = basis, levels = NULL)
}

## A model frame's rows `rows`, with the terms it was built from.
cutBasisFrame <- function(frame, rows) {
  cut <- frame[rows, , drop = FALSE]
  attr(cut, "terms") <- attr(frame, "terms")
  cut
}

## A basis written as code, built on the rows a fit keeps: `read` is its
## reading over every row of the data (readForestBasis) and `rows` the rows
## 'subset' and the na.action leave, NULL for all. Whether a column is a
## factor's, which levels have a row and whether a value is missing are all
## decided on the rows kept, at either door, while what scale() and poly()
## compute comes from every row, as in lm(). Returns list(basis = ,
## record = ), the record being what builds the same columns at new rows: R's
## terms with their environment, the levels of the frame's factors, the levels
## kept, the basis's text, the columns of the data it names and the number of
## rows it was read over.
buildCodeBasis <- function(read, rows = NULL) {
  frame <- read$frame
  terms <- attr(frame, "terms")
  kept <- if (is.null(rows)) frame else cutBasisFrame(frame, rows)
  built <- basisColumns(terms, kept, read$label)
  list(
    basis = built$basis,
    record = list(
      terms = terms,
      xlev = stats::.getXlevels(terms, frame),
      levels = built$levels,
      label = read$label,
      columns = read$columns,
      rows = nrow(frame)
    )
  )
}

## The label of every forest of a model of several, fixed when the sampler is
## created. A list name is the label; without one it is the text of the
## forest's basis; without that, a forest with no basis or with one handed
## over as a value or carried by the data object, it is forest<i>, i the
## forest's position. Labels are unique: two list names alike are refused, as
## is a list name of the form forest<digits> that is not its own position's,
## and a basis text that repeats an earlier label takes make.unique()'s
## suffix, as the names of a data frame do.
forestLabels <- function(listNames, basisTexts, numForests) {
  given <- rep("", numForests)
  if (!is.null(listNames)) {
    count <- min(length(listNames), numForests)
    given[seq_len(count)] <- listNames[seq_len(count)]
    given[is.na(given)] <- ""
  }
  defaults <- paste0("forest", seq_len(numForests))
  named <- nzchar(given)
  if (anyDuplicated(given[named]) > 0L) {
    stop(
      "'forests' names two forests \"",
      given[named][duplicated(given[named])][[1L]],
      "\"; a label names one forest",
      call. = FALSE
    )
  }
  reserved <- named & grepl("^forest[0-9]+$", given) & given != defaults
  if (any(reserved)) {
    index <- which(reserved)[[1L]]
    stop(
      "'forests' names forest ",
      index,
      " \"",
      given[[index]],
      "\", which is the label an unnamed forest ",
      sub("^forest", "", given[[index]]),
      " has; choose another name",
      call. = FALSE
    )
  }
  texts <- rep("", numForests)
  for (index in seq_len(min(length(basisTexts), numForests))) {
    if (!is.null(basisTexts[[index]])) {
      texts[index] <- basisTexts[[index]]
    }
  }
  derived <- !named & nzchar(texts)
  plain <- !named & !derived
  labels <- given
  labels[plain] <- defaults[plain]
  taken <- labels[named | plain]
  for (index in which(derived)) {
    labels[index] <- utils::tail(make.unique(c(taken, texts[index])), 1L)
    taken <- c(taken, labels[index])
  }
  labels
}

## $setForestBasis takes columns by position, as a matrix is taken in base R,
## and the names recorded at creation stay. The one case in which names and
## positions visibly disagree is refused: a replacement whose column names are
## the recorded ones in another order.
refuseReorderedBasisNames <- function(value, current, forestIndex) {
  given <- colnames(value)
  recorded <- colnames(current)
  if (
    !is.null(given) &&
      !is.null(recorded) &&
      !anyNA(given) &&
      length(given) == length(recorded) &&
      setequal(given, recorded) &&
      !identical(given, recorded)
  ) {
    stop(
      "'basis' has the columns of forest ",
      forestIndex,
      "'s basis in another order (",
      paste(given, collapse = ", "),
      "; the forest has ",
      paste(recorded, collapse = ", "),
      "): columns are taken by position, so give them in the forest's order",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## A basis that reached forest() through dots whose writer is no frame on the
## stack, so that its code cannot be read where it was written: `value`, what
## the code gives there, is handed on as a value, and code whose top is a
## model formula's own, which evaluated plainly would be arithmetic, is
## refused.
forwardedBasisValue <- function(code, value) {
  if (
    is.call(code) &&
      is.name(code[[1L]]) &&
      as.character(code[[1L]]) %in% c("+", "-", ":", "(")
  ) {
    stop(
      "forest()'s 'basis' (",
      shownCode(code),
      ") was passed on through '...' from a call that is no longer running, ",
      "so it cannot be read where it was written; hand the code over, as ",
      "do.call(forest, list(basis = quote(",
      shownCode(code),
      ")))",
      call. = FALSE
    )
  }
  value
}

## The data a fit's bases are read against and the number of its rows, every
## row the fit was given, before 'subset' and the na.action: a formula's
## data frame, list or environment, and with the matrix interface no data and
## the rows of the predictors. Without a data frame the count is the length
## of the first variable the formula names, as it is for 'subset'
## (resolveFormulaBasisSubset).
fitBasisRows <- function(formula, data) {
  if (!is.formula(formula)) {
    return(list(data = NULL, full = NROW(formula)))
  }
  if (!is.list(data) && !is.environment(data)) {
    data <- NULL
  }
  if (is.data.frame(data)) {
    return(list(data = data, full = nrow(data)))
  }
  vars <- setdiff(all.vars(formula), ".")
  first <- if (length(vars) == 0L) {
    NULL
  } else if (is.null(data)) {
    eval(as.name(vars[1L]), environment(formula))
  } else {
    eval(as.name(vars[1L]), data, environment(formula))
  }
  list(data = data, full = NROW(first))
}

## The rows of a formula fit's data that its 'subset' keeps, of `full` rows,
## in the order it keeps them: the expression evaluated in the data and then
## in the formula's environment and used to subscript the rows, their names
## among it, which is the reading the fit's own model frame gives it. NULL
## for every row.
formulaSubsetRows <- function(formula, data, subsetExpr, full) {
  if (is.null(subsetExpr)) {
    return(NULL)
  }
  subset <- if (is.null(data)) {
    eval(subsetExpr, environment(formula))
  } else {
    eval(subsetExpr, data, environment(formula))
  }
  rows <- seq_len(full)
  if (is.character(subset)) {
    names(rows) <- if (is.data.frame(data)) row.names(data) else rows
  }
  unname(rows[subset])
}

## The bases of a fit's forests from their readings: each written as code is
## built on `rows` of its frame (buildCodeBasis), NULL for all, and takes its
## place in `bases`, which holds the values already restricted to those rows
## or is NULL. Returns `bases` and, beside them, `records`, the record of
## each code basis. A basis has one row for each observation the fit keeps.
buildFitBases <- function(reads, bases, rows, numObservations) {
  numForests <- length(reads)
  if (is.null(bases)) {
    bases <- vector("list", numForests)
  }
  records <- vector("list", numForests)
  for (index in seq_len(numForests)) {
    read <- reads[[index]]
    if (is.null(read$frame)) {
      next
    }
    built <- buildCodeBasis(read, rows)
    if (nrow(built$basis) != numObservations) {
      stop(
        "forest ",
        index,
        "'s 'basis' (",
        read$label,
        ") has ",
        nrow(built$basis),
        " rows where the fit has ",
        numObservations,
        " observations; a basis has one row for each",
        call. = FALSE
      )
    }
    bases[index] <- list(built$basis)
    records[index] <- list(built$record)
  }
  list(bases = bases, records = records)
}

## The basis $setForestBasis installs, expanded to its columns over the
## sampler's `numRows` observations. A one-sided formula is read as the basis
## of a forest() is, with no data, so its names are found where it was
## written; anything else is a value. A level that no row has keeps its
## column, a replacement being free to leave one unobserved for a while, so
## a factor's columns here are every level's.
replacementForestBasis <- function(basis, numRows) {
  if (!inherits(basis, "formula")) {
    refuseSingleBasisValue(basis)
    return(expandValueBasis(basis, allowEmptyLevels = TRUE))
  }
  read <- readForestBasis(basisCode(basis, environment(basis)), NULL, numRows)
  if (is.null(read)) {
    stop("'basis' cannot be NULL")
  }
  frame <- read$frame
  value <- if (length(frame) == 1L) frame[[1L]]
  if (
    is.factor(value) ||
      is.character(value) ||
      (is.logical(value) && is.null(dim(value)))
  ) {
    # every level of the factor as written, one with no row among them
    expanded <- expandForestBasis(value, allowEmptyLevels = TRUE)
    colnames(expanded) <- paste0(
      attr(attr(frame, "terms"), "term.labels"),
      if (is.logical(value)) {
        c("FALSE", "TRUE")
      } else if (is.factor(value)) {
        levels(value)
      } else {
        levels(factor(value))
      }
    )
    return(expanded)
  }
  buildCodeBasis(read)$basis
}
