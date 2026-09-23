## Response families as objects. A fitting function takes one 'family'
## argument; every setting that only one family reads - a Student-t degrees
## of freedom, a count dispersion, a discrete-time hazard's period grid -
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
## constructors inside 'family'.
vocabularyEnv <- function(vocabulary, evalEnv) {
  env <- new.env(parent = evalEnv)
  for (name in names(vocabulary)) {
    assign(name, vocabulary[[name]], envir = env)
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
  while (is.symbol(expr) && grepl("^\\.\\.[1-9][0-9]*$", expr)) {
    recovered <- tryCatch(
      {
        ## a closure can reference dots its enclosing function owns
        while (!exists("...", envir = env, inherits = FALSE)) {
          env <- parent.env(env)
        }
        ## the frame's first place on the stack is its own call; later ones
        ## are evaluations in it (eval(call, env)), whose caller is not the
        ## one that wrote the dots
        onStack <- vapply(sys.frames(), identical, NA, env)
        frame <- which(onStack)[1L]
        parent <- sys.parents()[frame]
        ## a caller that is no frame on the stack, do.call(envir = ), numbers
        ## as the frame itself. parent.frame() reads the most recent context
        ## on the frame, which is its own call only when nothing has
        ## evaluated in it since; otherwise the caller is unknown and the
        ## reference is left as it is rather than guessed
        caller <- if (parent < frame) {
          sys.frame(parent)
        } else if (sum(onStack) == 1L) {
          do.call(parent.frame, list(), envir = env)
        } else {
          stop("unknown caller")
        }
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

## An entry point's unevaluated argument evaluated with a package vocabulary
## shadowing the environment the argument was written in.
evalInVocabulary <- function(expr, vocabulary, evalEnv) {
  written <- recoverForwardedArgument(expr, evalEnv)
  eval(written$expr, vocabularyEnv(vocabulary, written$env))
}

## The residual scale's own prior, carried by every family that draws one.
## NULL leaves the package default (chisq(3, 0.9)) standing; a bare
## constructor name means its defaults; anything else must be a residual
## prior object. The name is the PARAMETER, as 'dispersion' and 'df' are on
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
gaussian <- function(sigma = NULL) {
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

## Negative-binomial counts. dispersion = NA estimates the dispersion r;
## a positive value fixes it.
nbinom <- function(dispersion = NA) {
  # the default is the logical NA, "estimate it", which is not numeric
  if (
    length(dispersion) != 1L ||
      (!is.na(dispersion) &&
        (!is.numeric(dispersion) ||
          !is.finite(dispersion) ||
          dispersion <= 0.0))
  ) {
    stop(
      "nbinom 'dispersion' must be NA (estimate it) or a single positive ",
      "finite number"
    )
  }
  newValidated(
    "dbartsFamily",
    token = "nbinom",
    settings = list(dispersion = as.double(dispersion))
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
## runs in 0.57 seconds at the cap and 1.6 seconds at three times it. Memory
## is - about 190 MB per million rows at ten predictor columns, and it scales
## with the column count, so the default cap holds 1.9 GB of expanded design
## before the sampler has allocated anything of its own, and a caller who
## raises it to 3e7 holds 5.7 GB. Ten million rows is where the expansion
## stops being something a 16 GB machine absorbs; above it the refusal names
## both levers, coarsening the grid and raising the cap, because which one is
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

## The semicontinuous two-part model: an occupancy probit on 1{y > 0} and a
## gaussian on log(y) over the positive part. A composition of two samplers,
## so only bart() fits it.
## 'sigma' is the positive part's residual prior; the occupancy probit has a
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
## the family untouched. A family with no free residual scale carries the
## setting inertly - the fixed-unit-scale rule overwrites it downstream.
withResidPrior <- function(family, residPrior) {
  if (!is.null(residPrior)) {
    family@settings$sigma <- residPrior
  }
  family
}

## 'sigest' is the residual-scale estimate a chisq prior calibrates against:
## its quantile says where the estimate falls, so the two belong together.
## A fixed residual scale has nothing to calibrate - the engine overwrites
## the estimate with the square root of the fixed variance - so the pair is
## refused rather than accepted and ignored.
## The estimate itself is the test, not its name in the call: every entry
## point forwards 'sigest' to the one below it, defaulted to NA, so a name
## is no evidence a caller wrote one. 'sigestName' is the spelling the
## caller actually wrote - dbarts()/dbartsSpec() still read the retired
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
resolveFamily <- function(expr, tokens, caller, evalEnv) {
  if (is.null(expr)) {
    return(newValidated("dbartsFamily", token = tokens[1L]))
  }

  ## the prior constructors join the family vocabulary: the residual prior
  ## rides the family object, so gaussian(sigma = chisq(3, 0.9)) has to
  ## resolve 'chisq' here the same way parsePriors resolves it inside
  ## 'resid.prior'. The two name sets are disjoint.
  value <- evalInVocabulary(expr, c(dbartsFamilies, dbartsPriors), evalEnv)
  ## a bare constructor name (family = probit) means its defaults
  if (is.function(value)) {
    value <- value()
  }

  if (is.character(value)) {
    if (length(value) == 0L) {
      stop("'family' must name one response family")
    }
    ## one model, one token: "twopart" is a retired spelling, refused by name
    ## ahead of the generic list rather than folded in silence
    if (identical(value, "twopart")) {
      refuseTwopartFamily(caller)
    }
    ## a wrapper forwarding its own unevaluated formal hands over the whole
    ## default vector, whose first element is the default match.arg reads;
    ## anything else is one token, matched partially as match.arg does
    token <- if (length(value) > 1L) value[1L] else match.arg(value, tokens)
    value <- newValidated("dbartsFamily", token = token)
  }

  if (!is(value, "dbartsFamily")) {
    stop(
      "'family' must be a family name or a family object; see ?dbartsFamilies"
    )
  }
  refuseUnsupportedFamily(value@token, tokens, caller)
  value
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
