## Response families as objects. A fitting function takes one 'family'
## argument; every setting that only one family reads - a Student-t degrees
## of freedom, a count shape, a discrete-time hazard's period grid -
## rides the object rather than a formal that is inert for every other
## family. glm()'s idiom, with the same two spellings: a token string
## ("gaussian") means the family at its defaults, and a call
## (student(df = 3)) means the family with settings.
##
## The constructors are NOT exported, for the reason the prior constructors
## are not: gaussian, student, probit and logistic are generic enough names
## that exporting them would mask - and be masked by - other attached
## packages, stats::gaussian first among them. They resolve by bare name
## inside the 'family' argument of every entry point that takes one
## (resolveFamily below), and dbartsFamilies is their exported face.

## An environment in which a package vocabulary resolves by bare name and
## everything else falls through to the caller's own frame: the prior
## constructors inside the prior arguments (parsePriors), the family
## constructors inside 'family'. Marked, so that code a constructor keeps
## unevaluated is kept with the caller's frame and not with this one
## (callingPlace).
vocabularyEnv <- function(vocabulary, evalEnv) {
  env <- new.env(parent = evalEnv)
  for (name in names(vocabulary)) {
    assign(name, vocabulary[[name]], envir = env)
  }
  attr(env, "dbarts.vocabulary") <- TRUE
  env
}

## Where a constructor's caller wrote the call: `env`, or the frame under it
## when `env` is one a vocabulary was layered over.
callingPlace <- function(env) {
  while (isTRUE(attr(env, "dbarts.vocabulary", exact = TRUE))) {
    env <- parent.env(env)
  }
  env
}

## An argument's expression as written, and the environment it was written
## in. match.call() records an argument forwarded through a wrapper's dots as
## ..N, and evaluating that forces the wrapper caller's promise outside any
## vocabulary. Each forwarding frame's own call names the Nth dots element,
## so the walk back reaches the original expression; it stops where that is
## not possible and leaves the reference to evaluate as an ordinary one.
recoverForwardedArgument <- function(expr, env) {
  while (isDotsReference(expr)) {
    recovered <- tryCatch(
      {
        ## a closure can reference dots its enclosing function owns
        while (!exists("...", envir = env, inherits = FALSE)) {
          env <- parent.env(env)
        }
        ## the frame's first place on the stack is its own call; later ones
        ## are evaluations in it (eval(call, env)), whose caller is not the
        ## one that wrote the dots
        frame <- Position(function(f) identical(f, env), sys.frames())
        parent <- sys.parents()[frame]
        ## NextMethod(name = value) replaces the method's dots, but the frame
        ## still records the generic's call
        if (frame > 1L && identical(sys.function(frame - 1L), NextMethod)) {
          stop("dots replaced by NextMethod")
        }
        ## a caller that is no frame on the stack, do.call(envir = ), numbers
        ## as the frame itself and cannot be named
        if (parent >= frame) {
          stop("unknown caller")
        }
        caller <- sys.frame(parent)
        dots <- match.call(
          sys.function(frame),
          sys.call(frame),
          expand.dots = FALSE,
          envir = caller
        )$...
        list(dots[[as.integer(substring(expr, 3L))]], caller)
      },
      error = function(e) NULL
    )
    if (is.null(recovered)) {
      break
    }
    expr <- recovered[[1L]]
    env <- recovered[[2L]]
  }
  list(expr = expr, env = env)
}

## A matched call as its caller wrote it, for storing on a fit: each argument
## forwarded through a wrapper's dots (..N) is replaced by the expression it
## was written as, so update() re-evaluates that rather than a reference to
## dots that no longer exist. A reference the walk cannot resolve is left as
## it stands. Must run while the forwarding frames are still on the stack.
expandForwardedCall <- function(call, env) {
  for (i in seq_along(call)[-1L]) {
    if (isDotsReference(call[[i]])) {
      call[i] <- list(recoverForwardedArgument(call[[i]], env)$expr)
    }
  }
  call
}

isDotsReference <- function(expr) {
  is.symbol(expr) && grepl("^\\.\\.[1-9][0-9]*$", expr)
}

## An entry point's unevaluated argument evaluated with a package vocabulary
## layered over evalEnv, then passed through 'resolve', the site's own
## normalization, which signals an error for a value the site refuses. A
## forwarded reference (..N) is evaluated as it stands first; the expression
## it was written as is recovered, and evaluated where it was written, only
## when that fails. Recovery can turn a failure into a value but never changes
## a value; the error raised is the recovered expression's own, or the first
## one when there is nothing to recover. The cost is that an expression that failed part way is evaluated again, side
## effects included, by every site that reads it; R's warning on forcing the
## failed promise again is expected here and muffled. A constructor call
## forced outside the argument that takes it - a wrapper's named formal
## passed on unevaluated, then forced where it was written - fails with R's
## own message, extended by forceCallerCode's hint, drawn only from
## 'hintLists', the vocabulary list(s) this site actually reads constructors
## from (every prior site but 'family' takes only dbartsPriors).
evalInVocabulary <- function(
  expr,
  vocabulary,
  evalEnv,
  resolve = identity,
  hintLists = "dbartsPriors"
) {
  evalIn <- function(expr, env) {
    resolve(forceCallerCode(
      eval(expr, vocabularyEnv(vocabulary, env)),
      hintLists = hintLists
    ))
  }
  if (!isDotsReference(expr)) {
    return(evalIn(expr, evalEnv))
  }
  restarted <- gettext(
    "restarting interrupted promise evaluation",
    domain = "R"
  )
  first <- function() {
    withCallingHandlers(evalIn(expr, evalEnv), warning = function(w) {
      if (identical(conditionMessage(w), restarted)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  tryCatch(first(), error = function(original) {
    written <- recoverForwardedArgument(expr, evalEnv)
    if (isDotsReference(written$expr)) {
      stop(original)
    }
    evalIn(written$expr, written$env)
  })
}

## An argument taking a forest constructor, evaluated over the site's
## vocabulary. Unlike evalInVocabulary, a vocabulary name in value position is
## the caller's when the caller binds it (callerBinding); every other
## vocabulary name, and any in call position, is the constructor, so neither an
## attached mask nor a caller's variable of the same name changes what a call
## builds. A top-level value that is a constructor is called for its defaults.
## A ..N anywhere in the expression is evaluated as it stands and, on failure,
## recovered and resolved under the same rule where it was written. There is
## no other recovery: a constructor evaluated outside the argument fails with
## R's own message, extended by a hint where dbarts forces the code.
evalInForestVocabulary <- function(expr, vocabulary, evalEnv) {
  expr <- inlineForwardedArguments(expr, vocabulary, evalEnv)
  env <- vocabularyEnv(vocabulary, evalEnv)
  claimed <- character()
  for (name in valueSymbols(expr, names(vocabulary))) {
    binding <- callerBinding(name, evalEnv)
    if (!is.null(binding)) {
      claimed <- c(claimed, name)
      assign(name, binding$value, envir = env)
    }
  }
  ## a claimed name's call position is still the constructor, which 'env' no
  ## longer binds
  expr <- inlineConstructorCalls(expr, vocabulary[claimed])
  value <- forceCallerCode(eval(expr, env))
  if (is.function(value) && any(vapply(vocabulary, identical, NA, value))) {
    value <- value()
  }
  value
}

## Forces caller code: R's warning on forcing a promise that failed before
## is expected and muffled, and a failure to find a constructor gains the
## spelling that works outside its argument - the exported list its name is
## reached through - but only among 'hintLists', the vocabulary list(s) the
## calling site actually reads constructors from. A name from an unrelated
## vocabulary (a family constructor forced at a prior argument, say) is left
## as R's own message: a name that argument could never have meant is not a
## hint, it is noise.
forceCallerCode <- function(
  value,
  hintLists = c("dbartsForests", "forestPriors")
) {
  restarted <- gettext(
    "restarting interrupted promise evaluation",
    domain = "R"
  )
  withCallingHandlers(
    tryCatch(value, error = function(e) {
      sources <- list(
        dbartsForests = dbartsForests,
        dbartsPriors = dbartsPriors,
        dbartsFamilies = dbartsFamilies,
        forestPriors = dbartsPriors["fixed"]
      )[hintLists]
      name <- unlist(lapply(sources, names), use.names = FALSE)
      topic <- sub(
        "^forestPriors$",
        "dbartsPriors",
        rep(names(sources), lengths(sources))
      )
      missingFunction <- gettextf(
        "could not find function \"%s\"",
        name,
        domain = "R"
      )
      i <- match(conditionMessage(e), missingFunction)
      if (!is.na(i)) {
        e$message <- paste0(
          missingFunction[i],
          "; outside the argument that takes it, write ",
          topic[i],
          "$",
          name[i],
          "(...)"
        )
      }
      stop(e)
    }),
    warning = function(w) {
      if (identical(conditionMessage(w), restarted)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

## The caller's value for a vocabulary name in value position, as list(value
## = ), or NULL when the caller does not claim it. The lookup stops at
## topenv(env), so neither the search path nor a package wrapper's user
## globals count. An unsupplied formal with an empty default means NULL, the
## door default; any other binding is forced, and a function is no site's
## value, so it leaves the constructor. That also discounts the dbarts
## namespace's own binding.
callerBinding <- function(name, env) {
  top <- topenv(env)
  while (!exists(name, envir = env, inherits = FALSE)) {
    if (identical(env, top) || identical(env, emptyenv())) {
      return(NULL)
    }
    env <- parent.env(env)
  }
  frame <- Position(function(f) identical(f, env), sys.frames())
  if (
    !is.na(frame) &&
      identical(formals(sys.function(frame))[[name]], quote(expr = )) &&
      eval(call("missing", as.name(name)), env)
  ) {
    return(list(value = NULL))
  }
  value <- forceCallerCode(get(name, envir = env, inherits = FALSE))
  if (is.function(value)) NULL else list(value = value)
}

## A formula, quote() or bquote() holds language rather than values, so the
## walks below leave it as written.
isQuotingCall <- function(expr) {
  is.call(expr) &&
    is.symbol(expr[[1L]]) &&
    as.character(expr[[1L]]) %in% c("~", "quote", "bquote")
}

## The vocabulary names an expression uses in value position.
valueSymbols <- function(expr, vocabularyNames) {
  if (is.symbol(expr)) {
    return(intersect(as.character(expr), vocabularyNames))
  }
  if (!is.call(expr) || isQuotingCall(expr)) {
    return(character())
  }
  parts <- as.list(expr)
  if (is.symbol(parts[[1L]])) {
    parts <- parts[-1L]
  }
  unique(unlist(lapply(parts, valueSymbols, vocabularyNames)))
}

inlineConstructorCalls <- function(expr, constructors) {
  if (length(constructors) == 0L || !is.call(expr) || isQuotingCall(expr)) {
    return(expr)
  }
  head <- expr[[1L]]
  if (is.symbol(head) && as.character(head) %in% names(constructors)) {
    expr[[1L]] <- constructors[[as.character(head)]]
  }
  as.call(lapply(as.list(expr), inlineConstructorCalls, constructors))
}

## Replaces each ..N in an expression with its resolved value, quoted so that
## a language value (a formula) is not evaluated again. A function literal's
## ..N are its own dots.
inlineForwardedArguments <- function(expr, vocabulary, env) {
  if (isDotsReference(expr)) {
    value <- tryCatch(
      forceCallerCode(eval(expr, env)),
      error = function(original) {
        written <- recoverForwardedArgument(expr, env)
        if (isDotsReference(written$expr)) {
          stop(original)
        }
        evalInForestVocabulary(written$expr, vocabulary, written$env)
      }
    )
    return(call("quote", value))
  }
  if (
    !is.call(expr) ||
      isQuotingCall(expr) ||
      identical(expr[[1L]], as.name("function"))
  ) {
    return(expr)
  }
  as.call(lapply(as.list(expr), inlineForwardedArguments, vocabulary, env))
}

## A site's 'resolve' for evalInVocabulary: a bare constructor name means its
## defaults, and a value of none of 'classes' is refused by name.
resolvedAs <- function(name, classes, what, topic = "dbartsPriors") {
  function(value) {
    if (is.function(value)) {
      value <- value()
    }
    if (!any(vapply(classes, is, NA, object = value))) {
      stop("'", name, "' must be a ", what, "; see ?", topic, call. = FALSE)
    }
    value
  }
}

## The residual scale's own prior, carried by every family that draws one.
## NULL leaves the package default (chisq(3, 0.9)) standing; a bare
## constructor name means its defaults; anything else must be a residual
## prior object. The name is the PARAMETER, as 'shape' and 'df' are on
## their families: 'sigest' is the estimate supplied at creation and
## $setSigma writes the parameter itself, so 'sigma' here is the law that
## parameter is drawn under.
validateFamilySigma <- function(sigma, caller) {
  if (is.null(sigma)) {
    return(NULL)
  }
  if (is.function(sigma)) {
    sigma <- sigma()
  }
  if (!is(sigma, "dbartsResidPrior")) {
    stop(
      caller,
      " 'sigma' must be a residual prior - chisq(df, quant) or ",
      "fixed(value); see ?dbartsPriors"
    )
  }
  sigma
}

## The settings list is sparse, so an unsupplied residual prior adds no entry
## at all and a family built from a token stays indistinguishable from one
## built by its constructor.
familySigmaSetting <- function(sigma, caller) {
  sigma <- validateFamilySigma(sigma, caller)
  if (is.null(sigma)) list() else list(sigma = sigma)
}

## The Gaussian (continuous) response, the package default for a numeric
## response under family = "auto". 'sigma' is the residual scale's prior.
## dbarts fits only the identity link, so this takes none; an R family object
## such as stats::gaussian() is mapped where family objects are resolved.
gaussian <- function(sigma = NULL, ...) {
  # a link however spelled - named or partially named, a string, or a bare
  # link name read before it would be evaluated - gets the link message
  sigmaExpr <- substitute(sigma)
  dotNames <- as.character(...names())
  dotStrings <- vapply(
    as.list(substitute(list(...)))[-1L],
    is.character,
    logical(1L)
  )
  if (
    (is.name(sigmaExpr) && as.character(sigmaExpr) %in% statsLinkNames) ||
      any(nzchar(dotNames) & startsWith("link", dotNames)) ||
      any(dotStrings) ||
      is.character(sigma)
  ) {
    stop(
      "gaussian() takes no 'link': dbarts fits the identity link only; drop ",
      "it, or pass R's stats::gaussian()",
      call. = FALSE
    )
  }
  if (...length() > 0L) {
    stop("gaussian() takes only 'sigma'", call. = FALSE)
  }
  newValidated(
    "dbartsFamily",
    token = "gaussian",
    settings = familySigmaSetting(sigma, "gaussian")
  )
}

## Outlier-robust Student-t errors, drawn by the Gaussian scale-mixture
## augmentation. df = NULL estimates the degrees of freedom on a capped
## grid; a finite positive value fixes them. A continuous response only:
## every other family has its own fixed latent scale.
student <- function(df = NULL, sigma = NULL) {
  if (!is.null(df)) {
    if (
      !is.numeric(df) ||
        length(df) != 1L ||
        is.na(df) ||
        !is.finite(df) ||
        df <= 0.0
    ) {
      stop(
        "student 'df' must be NULL (estimate the degrees of freedom) or a ",
        "single positive finite number"
      )
    }
  }
  newValidated(
    "dbartsFamily",
    token = "student",
    settings = c(
      list(df = if (is.null(df)) NA_real_ else as.double(df)),
      familySigmaSetting(sigma, "student")
    )
  )
}

probit <- function() {
  newValidated("dbartsFamily", token = "probit")
}

logistic <- function() {
  newValidated("dbartsFamily", token = "logistic")
}

multinomial <- function() {
  newValidated("dbartsFamily", token = "multinomial")
}

ordinal <- function() {
  newValidated("dbartsFamily", token = "ordinal")
}

## Negative-binomial counts. shape = NULL estimates the shape r;
## a positive value fixes it. The setting stores NA_real_ for "estimate".
nbinom <- function(shape = NULL) {
  if (is.null(shape)) {
    shape <- NA_real_
  } else if (
    length(shape) != 1L ||
      !is.numeric(shape) ||
      !is.finite(shape) ||
      shape <= 0.0 ||
      shape != round(shape)
  ) {
    stop(
      "nbinom 'shape' must be NULL (estimate it) or a single positive ",
      "whole number"
    )
  }
  newValidated(
    "dbartsFamily",
    token = "nbinom",
    settings = list(shape = as.double(shape))
  )
}

## Accelerated failure time (log-normal) survival. The log-time residual
## scale is drawn as a gaussian one is, so it carries the same prior.
aft <- function(sigma = NULL) {
  newValidated(
    "dbartsFamily",
    token = "aft",
    settings = familySigmaSetting(sigma, "aft")
  )
}

## Discrete-time hazard: (time, status) are person-period expanded onto a
## grid of period breaks and fit as an ordinary binary model under 'link'.
## breaks = NULL takes the grid from the observed event times.
##
## max.rows caps the expansion, whose row count grows with the grid's
## resolution. Time is not what the cap protects: the expansion is linear and
## cheap even well past the cap. Memory is - about 190 MB per million rows at
## ten predictor columns, and it scales with the column count, so the default
## cap bounds the expanded design's memory before the sampler has allocated
## anything of its own. Ten million rows is where that bound stops being
## something a modest machine absorbs; above it the refusal names both
## levers, coarsening the grid and raising the cap, because which one is
## right depends on the column count and the host.
hazard <- function(
  breaks = NULL,
  max.rows = 1e7,
  link = c("probit", "logistic")
) {
  link <- match.arg(link)
  if (!is.null(breaks) && !is.numeric(breaks)) {
    stop("hazard 'breaks' must be NULL or numeric")
  }
  if (
    length(max.rows) != 1L ||
      !is.numeric(max.rows) ||
      is.na(max.rows) ||
      max.rows <= 0.0
  ) {
    stop("hazard 'max.rows' must be a single positive number")
  }
  newValidated(
    "dbartsFamily",
    token = paste0("hazard.", link),
    settings = list(breaks = breaks, max.rows = as.double(max.rows))
  )
}

## The semicontinuous two-part model: a zero-part probit on 1{y > 0} and a
## gaussian on log(y) over the positive part. A composition of two samplers,
## so only bart() fits it.
## 'sigma' is the positive part's residual prior; the zero-part probit has a
## fixed unit latent scale and takes none.
hurdle.lognormal <- function(sigma = NULL) {
  newValidated(
    "dbartsFamily",
    token = "hurdle.lognormal",
    settings = familySigmaSetting(sigma, "hurdle.lognormal")
  )
}

## The exported face of the family constructors: one object, so no generic
## name enters the search path. Inside the 'family' argument of the fitting
## functions the same constructors resolve by bare name.
dbartsFamilies <- list(
  gaussian = gaussian,
  student = student,
  probit = probit,
  logistic = logistic,
  multinomial = multinomial,
  ordinal = ordinal,
  nbinom = nbinom,
  aft = aft,
  hazard = hazard,
  hurdle.lognormal = hurdle.lognormal
)

## A residual prior reads back as the constructor call that built it; the
## default deparse of an S4 object would show its class and slots instead.
formatResidPrior <- function(prior) {
  if (is(prior, "dbartsChiSqPrior")) {
    paste0("chisq(", prior@df, ", ", prior@quantile, ")")
  } else if (is(prior, "dbartsFixedPrior")) {
    paste0("fixed(", prior@value, ")")
  } else {
    class(prior)[1L]
  }
}

## The residual prior a door resolved from a retired flat spelling, stamped
## onto the family object that door forwards: the prior has one home, so a
## door that still reads an old spelling has to put it there. NULL leaves
## the family untouched, as does a family with no free residual scale, whose
## constructor takes no 'sigma'.
withResidPrior <- function(family, residPrior) {
  if (!is.null(residPrior) && familyTakesSetting(family@token, "sigma")) {
    family@settings$sigma <- residPrior
  }
  family
}

## 'sigest' is the residual-scale estimate a chisq prior calibrates against:
## its quantile says where the estimate falls, so the two belong together.
## A fixed residual scale has nothing to calibrate - the engine overwrites
## the estimate with the square root of the fixed variance - so a 'sigest'
## beside one is refused when it differs from that square root. One that
## equals it is accepted with a message, once per session, and is an error
## from the tombstone expiry whatever its value. 0.9-34 documented no fixed
## residual prior; reached internally, it stopped sigma being drawn, so the
## chain sat at its starting value, 'sigma', and loops that wrote the two
## equal ran at the intended scale.
## "Equals" is within 4 machine epsilons, relative to sqrt(fixed variance),
## not exact: sqrt is correctly rounded, so a caller who writes sqrt(v)
## matches exactly, but one who writes 0.3 beside fixed(0.09) is off by an
## ulp or two from the rounding of 0.3^2; anything looser would accept a
## different scale.
## The estimate itself is the test, not its name in the call: every entry
## point forwards 'sigest' to the one below it, defaulted to NA, so a name
## is no evidence a caller wrote one. 'sigestName' is the spelling the
## caller actually wrote - dbarts() and xbart() still read the retired
## 'sigma =' for one release - and defaults to 'sigest' at doors that carry
## no other spelling for it.
refuseSigestUnderFixedPrior <- function(
  residPrior,
  sigest,
  sigestName = "sigest"
) {
  if (
    is(residPrior, "dbartsFixedPrior") &&
      length(sigest) == 1L &&
      !is.na(sigest)
  ) {
    fixedSd <- sqrt(residPrior@value)
    if (
      is.numeric(sigest) &&
        abs(sigest - fixedSd) <= 4 * .Machine$double.eps * fixedSd
    ) {
      if (!isTRUE(onceWarnState[["tombstone.sigest.fixed"]])) {
        onceWarnState[["tombstone.sigest.fixed"]] <- TRUE
        message(
          "'",
          sigestName,
          "' has no effect under a fixed residual scale (it agrees with ",
          formatResidPrior(residPrior),
          ") and is an error from dbarts ",
          tombstoneExpiry,
          "; drop it."
        )
      }
      return(invisible(NULL))
    }
    stop(
      "'",
      sigestName,
      "' has no effect under a fixed residual scale: sigma = ",
      formatResidPrior(residPrior),
      " IS the residual scale, and the estimate is overwritten with its ",
      "square root. Drop '",
      sigestName,
      "', or give a chisq prior for it to calibrate.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

formatFamilyCall <- function(token, settings) {
  ## hazard's link is folded into the token, so it is restated as the
  ## argument the caller would write rather than dropped from the display
  if (startsWith(token, "hazard.")) {
    settings <- c(settings, list(link = sub("^hazard\\.", "", token)))
    token <- "hazard"
  }
  described <- vapply(
    names(settings),
    function(name) {
      value <- settings[[name]]
      shown <- if (is(value, "dbartsResidPrior")) {
        formatResidPrior(value)
      } else if (name %in% c("shape", "df") && identical(value, NA_real_)) {
        # stored as NA_real_, spelled NULL: the constructors refuse the NA
        "NULL"
      } else {
        paste0(deparse(value), collapse = " ")
      }
      paste0(name, " = ", shown)
    },
    character(1L)
  )
  paste0(token, "(", paste0(described, collapse = ", "), ")")
}

methods::setMethod("show", "dbartsFamily", function(object) {
  cat(
    "dbarts response family: ",
    formatFamilyCall(object@token, object@settings),
    "\n",
    sep = ""
  )
  invisible(object)
})

## Resolves the 'family' argument of an entry point into one family object.
## `expr` is the caller's own unevaluated argument (matchedCall$family), so a
## bare constructor call resolves in the family vocabulary no matter what the
## caller has attached - the rule the prior vocabulary already follows - and
## an ordinary variable holding a token or an object still resolves in the
## caller's frame. Both hold for an argument forwarded through a wrapper's
## dots, which resolves where it was written. `tokens` is the entry point's own admissible list, its
## first element the default.
resolveFamily <- function(
  expr,
  tokens,
  caller,
  evalEnv,
  refused = character()
) {
  if (is.null(expr)) {
    return(newValidated("dbartsFamily", token = tokens[1L]))
  }

  ## the prior constructors join the family vocabulary: the residual prior
  ## rides the family object, so gaussian(sigma = chisq(3, 0.9)) has to
  ## resolve 'chisq' here the same way parsePriors resolves it inside
  ## 'resid.prior'. The two name sets are disjoint.
  value <- evalInVocabulary(
    expr,
    c(dbartsFamilies, dbartsPriors),
    evalEnv,
    resolvedFamily,
    hintLists = c("dbartsFamilies", "dbartsPriors")
  )

  if (is.character(value)) {
    if (length(value) == 0L) {
      stop("'family' must name one response family")
    }
    ## a wrapper forwarding its own unevaluated formal hands over the whole
    ## default vector, whose first element is the default match.arg reads;
    ## anything else is one token, matched partially as match.arg does
    token <- if (length(value) > 1L) {
      value[1L]
    } else if (value %in% statsFamilyNames) {
      resolvedFamily(get(value, envir = asNamespace("stats")))@token
    } else {
      ## a token the caller refuses by name is matched, so its own message
      ## can say why, but never listed among the choices of a bad token
      hit <- pmatch(value, c(tokens, refused))
      if (!is.na(hit) && hit > length(tokens)) {
        refused[hit - length(tokens)]
      } else {
        hit <- pmatch(value, tokens)
        if (is.na(hit)) {
          stop(
            "'family' should be one of ",
            paste0("\"", tokens, "\"", collapse = ", "),
            call. = FALSE
          )
        }
        tokens[hit]
      }
    }
    value <- newValidated("dbartsFamily", token = token)
  }
  if (value@token %in% refused) {
    return(value)
  }
  refuseUnsupportedFamily(value@token, tokens, caller)
  # a family object given whole is held to what its constructor takes, as one
  # built by that constructor already is
  refusedSettings <- familyRefusedSettings(value@token, names(value@settings))
  if (!is.null(refusedSettings)) {
    stop(refusedSettings, call. = FALSE)
  }
  value
}

## The link names stats::make.link knows, which a bare symbol may spell.
statsLinkNames <- c(
  "logit",
  "probit",
  "cauchit",
  "cloglog",
  "identity",
  "log",
  "sqrt",
  "1/mu^2",
  "inverse"
)

refusedLink <- function(name, link) {
  paste0(
    name,
    "(link = \"",
    link,
    "\") is not supported; dbarts fits ",
    if (identical(name, "gaussian")) {
      "the identity link"
    } else {
      "the probit and logit links"
    },
    "; see ?dbartsFamilies"
  )
}

## A family named by a string: the stats family functions map as the object
## does ("binomial" is binomial(), the logit link), and the rest of them are
## refused by name; any other string is the entry point's own token.
statsFamilyNames <- c(
  "binomial",
  "poisson",
  "Gamma",
  "inverse.gaussian",
  "quasi",
  "quasibinomial",
  "quasipoisson"
)

## A value of 'family' as an object: a function is called, as glm() does, and
## one of base R's family objects maps to the dbarts family with the same
## likelihood and link, or is refused by name. binomial's default link is
## logit, as in glm(), where a 0/1 response under "auto" is probit.
resolvedFamily <- function(value) {
  if (is.function(value)) {
    value <- value()
  }
  if (is.character(value) || is(value, "dbartsFamily")) {
    return(value)
  }
  if (!inherits(value, "family")) {
    stop(
      "'family' must be a family name or a family object; see ?dbartsFamilies",
      call. = FALSE
    )
  }
  name <- as.character(value$family)[1L]
  link <- as.character(value$link)[1L]
  if (identical(name, "gaussian") && identical(link, "identity")) {
    return(dbartsFamilies$gaussian())
  }
  if (identical(name, "binomial") && link %in% c("probit", "logit")) {
    return(
      if (identical(link, "probit")) {
        dbartsFamilies$probit()
      } else {
        dbartsFamilies$logistic()
      }
    )
  }
  if (name %in% c("gaussian", "binomial")) {
    stop(refusedLink(name, link), call. = FALSE)
  }
  stop(
    "family \"",
    name,
    "\" is not supported; see ?dbartsFamilies for the families dbarts fits",
    call. = FALSE
  )
}

## The per-entry-point family list, refused by name rather than by
## match.arg's generic message, since a family object's token never passed
## through match.arg in the first place.
refuseUnsupportedFamily <- function(token, tokens, caller) {
  if (token %in% tokens) {
    return(invisible(NULL))
  }
  stop(
    "'",
    caller,
    "' does not fit family \"",
    token,
    "\"; it takes ",
    paste0("\"", tokens, "\"", collapse = ", "),
    call. = FALSE
  )
}

## One family setting, or its documented absence. The settings list is
## sparse: a family object built from a token carries none at all.
familySetting <- function(family, name, default) {
  value <- family@settings[[name]]
  if (is.null(value)) default else value
}

## The settings a family token takes: its constructor's arguments, a link
## being folded into the token itself; NULL for a token no constructor builds,
## which an entry point refuses by name later. "auto" carries the residual
## prior a door without a family of its own stamps on before the response
## settles the family.
familyTokenSettings <- function(token) {
  if (identical(token, "auto")) {
    return("sigma")
  }
  constructor <- dbartsFamilies[[sub("^hazard\\..*$", "hazard", token)]]
  if (is.null(constructor)) {
    return(NULL)
  }
  # a constructor with no arguments takes no setting at all, character(0)
  # rather than the NULL that means unjudged
  setdiff(as.character(names(formals(constructor))), c("link", "..."))
}

## Whether a family token takes a setting; a token outside the constructors is
## not judged here.
familyTakesSetting <- function(token, name) {
  allowed <- familyTokenSettings(token)
  is.null(allowed) || name %in% allowed
}

## The validity message for settings a family token does not take, so a
## family object built by hand, or read back from a fit that predates this
## check, is refused by name rather than carried inertly; NULL when all fit.
familyRefusedSettings <- function(token, names) {
  extra <- names[!vapply(names, familyTakesSetting, logical(1L), token = token)]
  if (length(extra) == 0L) {
    return(NULL)
  }
  allowed <- familyTokenSettings(token)
  paste0(
    "family \"",
    token,
    "\" takes no setting ",
    paste0("'", extra, "'", collapse = ", "),
    if (length(allowed) > 0L) {
      paste0("; it takes ", paste0("'", allowed, "'", collapse = ", "))
    }
  )
}

## The family a fit was specified with once "auto" has resolved, with its
## settings: what family() on the fit returns, and whose token is the fit's
## $family. The engine family the link and likelihood follow, which a
## Student-t fit ("gaussian") and a hazard fit ("probit" or "logistic")
## remap, is looked up from it (fitEngineFamily). 'resolved' is the engine
## token "auto" resolved to.
specifiedFamily <- function(familySpec, resolved) {
  token <- familySpec@token
  if (identical(token, "auto")) {
    token <- resolved
  } else if (identical(token, "hazard")) {
    token <- "hazard.probit"
  }
  # a residual prior stamped on "auto" means nothing once it has resolved to
  # a family with a fixed latent scale, whose constructor takes none
  settings <- familySpec@settings
  settings <- settings[vapply(
    names(settings),
    familyTakesSetting,
    logical(1L),
    token = token
  )]
  # an emptied list keeps no names, as a constructor's own empty one has none
  if (length(settings) == 0L) {
    settings <- list()
  }
  newValidated("dbartsFamily", token = token, settings = settings)
}
