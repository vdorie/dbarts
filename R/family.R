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
## everything else falls through to the caller's own frame. The prior
## constructors already do this inside the prior arguments (parsePriors);
## the family constructors need it inside 'family'.
vocabularyEnv <- function(vocabulary, evalEnv) {
  env <- new.env(parent = evalEnv)
  for (name in names(vocabulary)) {
    assign(name, vocabulary[[name]], envir = env)
  }
  env
}

## The Gaussian (continuous) response, the package default for a numeric
## response under family = "auto".
gaussian <- function() {
  newValidated("dbartsFamily", token = "gaussian")
}

## Outlier-robust Student-t errors, drawn by the Gaussian scale-mixture
## augmentation. df = NULL estimates the degrees of freedom on a capped
## grid; a finite positive value fixes them. A continuous response only:
## every other family has its own fixed latent scale.
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
        "student 'df' must be NULL (estimate the degrees of freedom) or a ",
        "single positive finite number"
      )
    }
  }
  newValidated(
    "dbartsFamily",
    token = "student",
    settings = list(df = if (is.null(df)) NA_real_ else as.double(df))
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

## Accelerated failure time (log-normal) survival.
aft <- function() {
  newValidated("dbartsFamily", token = "aft")
}

## Discrete-time hazard: (time, status) are person-period expanded onto a
## grid of period breaks and fit as an ordinary binary model under 'link'.
## breaks = NULL takes the grid from the observed event times. max.rows caps
## the expansion, which is quadratic in the grid's resolution - ten million
## rows is roughly a gigabyte of design at double precision, so a grid that
## would exceed it is a specification mistake rather than a long run.
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
hurdle.lognormal <- function() {
  newValidated("dbartsFamily", token = "hurdle.lognormal")
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
      paste0(name, " = ", paste0(deparse(value), collapse = " "))
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
## caller's frame. `tokens` is the entry point's own admissible list, its
## first element the default.
resolveFamily <- function(expr, tokens, caller, evalEnv) {
  if (is.null(expr)) {
    return(newValidated("dbartsFamily", token = tokens[1L]))
  }

  value <- eval(expr, vocabularyEnv(dbartsFamilies, evalEnv))
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
