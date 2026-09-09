## Tombstones: every part of the 0.9-x public surface that is gone or
## renamed but still reachable for one release. Each entry names its
## successor and the version it expires at, and the whole set expires
## together, so the release that drops them deletes this file and nothing
## survives its expiry by accident. A tombstone never adds a capability:
## it errors, or it forwards to the successor after saying so once.

tombstoneExpiry <- "1.1-0"

## name: the spelling a 0.9-x caller writes.
## kind: "function" (exported stub or alias), "method" (registered S3
##   method), "rcMethod" (dbartsSampler reference-class method),
##   "argument" (a renamed formal accepted on 'owner'), "family" (a family
##   token), "behaviour" (a call shape handled for the release).
## successor: what to write instead; NA_character_ when there is none.
## The registry is asserted against NAMESPACE and the news file by
## inst/tinytest/test-tombstones.R.
dbartsTombstones <- list(
  list(
    name = "bart2",
    kind = "function",
    owner = NA_character_,
    successor = "bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "rbart_vi",
    kind = "function",
    owner = NA_character_,
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "predict.rbart",
    kind = "method",
    owner = "rbart",
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "extract.rbart",
    kind = "method",
    owner = "rbart",
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "fitted.rbart",
    kind = "method",
    owner = "rbart",
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "residuals.rbart",
    kind = "method",
    owner = "rbart",
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "startThreads",
    kind = "rcMethod",
    owner = "dbartsSampler",
    successor = NA_character_,
    expires = tombstoneExpiry
  ),
  list(
    name = "stopThreads",
    kind = "rcMethod",
    owner = "dbartsSampler",
    successor = NA_character_,
    expires = tombstoneExpiry
  ),
  list(
    name = "rngSeed",
    kind = "argument",
    owner = "bart",
    successor = "seed",
    expires = tombstoneExpiry
  ),
  list(
    name = "rngSeed",
    kind = "argument",
    owner = "dbartsControl",
    successor = "seed",
    expires = tombstoneExpiry
  ),
  list(
    name = "sigma",
    kind = "argument",
    owner = "dbarts",
    successor = "sigest",
    expires = tombstoneExpiry
  ),
  list(
    name = "sigma",
    kind = "argument",
    owner = "dbartsSpec",
    successor = "sigest",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.dist",
    kind = "argument",
    owner = "bart",
    successor = "family = student(df)",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.dist",
    kind = "argument",
    owner = "dbarts",
    successor = "family = student(df)",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.dist",
    kind = "argument",
    owner = "dbartsSpec",
    successor = "family = student(df)",
    expires = tombstoneExpiry
  ),
  list(
    name = "dispersion",
    kind = "argument",
    owner = "bart",
    successor = "family = nbinom(dispersion)",
    expires = tombstoneExpiry
  ),
  list(
    name = "dispersion",
    kind = "argument",
    owner = "dbarts",
    successor = "family = nbinom(dispersion)",
    expires = tombstoneExpiry
  ),
  list(
    name = "dispersion",
    kind = "argument",
    owner = "dbartsSpec",
    successor = "family = nbinom(dispersion)",
    expires = tombstoneExpiry
  ),
  list(
    name = "breaks",
    kind = "argument",
    owner = "bart",
    successor = "family = hazard(breaks)",
    expires = tombstoneExpiry
  ),
  list(
    name = "breaks",
    kind = "argument",
    owner = "dbarts",
    successor = "family = hazard(breaks)",
    expires = tombstoneExpiry
  ),
  list(
    name = "max.rows",
    kind = "argument",
    owner = "bart",
    successor = "family = hazard(max.rows)",
    expires = tombstoneExpiry
  ),
  list(
    name = "max.rows",
    kind = "argument",
    owner = "dbarts",
    successor = "family = hazard(max.rows)",
    expires = tombstoneExpiry
  ),
  list(
    name = "dart",
    kind = "argument",
    owner = "bart",
    successor = "tree.prior = dart()",
    expires = tombstoneExpiry
  ),
  list(
    name = "dart",
    kind = "argument",
    owner = "xbart",
    successor = "tree.prior = dart()",
    expires = tombstoneExpiry
  ),
  list(
    name = "levelGibbs",
    kind = "argument",
    owner = "bart",
    successor = "tree.prior = cgm(levelGibbs)",
    expires = tombstoneExpiry
  ),
  list(
    name = "twopart",
    kind = "family",
    owner = NA_character_,
    successor = "hurdle.lognormal",
    expires = tombstoneExpiry
  ),
  list(
    name = "BayesTree-spelled bart call",
    kind = "behaviour",
    owner = "bart",
    successor = "bartBT",
    expires = tombstoneExpiry
  ),
  list(
    name = "fourth positional bart argument",
    kind = "behaviour",
    owner = "bart",
    successor = "bartBT",
    expires = tombstoneExpiry
  ),
  list(
    name = "front-door startup message",
    kind = "behaviour",
    owner = ".onAttach",
    successor = NA_character_,
    expires = tombstoneExpiry
  )
)

## ------------------------------------------------------------------
## bart2, the one-release alias for the modern front door
## ------------------------------------------------------------------

## bart2's formals are bart's, copied rather than restated: a consumer
## reads formals(dbarts::bart2) for its defaults (bartCause does), so the
## alias has to be a real closure over the real list and not function(...).
## The matched call forwards with only the function replaced, so every
## argument is still promised in - and evaluated in - the caller's frame,
## and an unsupplied argument stays unsupplied.
bart2 <- function() {
  matchedCall <- match.call()
  matchedCall[[1L]] <- quote(dbarts::bart)
  warnOnce(
    "tombstone.bart2",
    "'bart2' is now 'bart'; this call was forwarded. The alias is removed ",
    "in dbarts ",
    tombstoneExpiry
  )
  eval(matchedCall, parent.frame())
}
formals(bart2) <- formals(bart)

## ------------------------------------------------------------------
## rbart_vi and its four generics
## ------------------------------------------------------------------

refuseGroupedRandomEffects <- function(what) {
  stop(
    what,
    " was removed in dbarts 1.0-0; grouped random effects live in the ",
    "stan4bart package (stan4bart::stan4bart). Its group-spread prior is ",
    "not the one dbarts drew from, so results move rather than reproduce. ",
    "This stub is removed in dbarts ",
    tombstoneExpiry,
    ".",
    call. = FALSE
  )
}

rbart_vi <- function(...) {
  refuseGroupedRandomEffects("'rbart_vi'")
}

predict.rbart <- function(object, ...) {
  refuseGroupedRandomEffects("'predict' on an rbart fit")
}

extract.rbart <- function(object, ...) {
  refuseGroupedRandomEffects("'extract' on an rbart fit")
}

fitted.rbart <- function(object, ...) {
  refuseGroupedRandomEffects("'fitted' on an rbart fit")
}

residuals.rbart <- function(object, ...) {
  refuseGroupedRandomEffects("'residuals' on an rbart fit")
}

## ------------------------------------------------------------------
## The BayesTree-spelled bart() call
## ------------------------------------------------------------------

## The names that select the legacy door: a formal of bartBT that bart
## does not itself take. Derived, so a name added to either signature
## cannot leave the shim behind; the list is pinned by
## inst/tinytest/test-tombstones.R.
bartBTOnlyFormals <- setdiff(
  names(formals(bartBT)),
  names(formals(bart))
)

## Fires on an exact BayesTree spelling only. A modern call that merely
## looks 0.9-x - three positional arguments, or nothing but shared names
## like 'k' and 'sigest' - is fit by the modern door with modern defaults,
## which is what dec-era callers of bart2 already get.
## The call forwarded is the one the caller WROTE, not bart's matched
## version of it: matching renames every positional argument to bart's own
## formals, which the legacy door does not have.
forwardToLegacyDoor <- function(suppliedCall, supplied, callingEnv) {
  legacy <- intersect(supplied, bartBTOnlyFormals)
  if (length(legacy) == 0L) {
    return(NULL)
  }
  warnOnce(
    "tombstone.bartShim",
    "'",
    legacy[1L],
    "' is dbarts 0.9-x's BayesTree-style 'bart' argument; that function is ",
    "now 'bartBT' and this call was forwarded to it. 'bart' is the modern ",
    "front door and takes different names and defaults. Forwarding is ",
    "removed in dbarts ",
    tombstoneExpiry,
    "."
  )
  suppliedCall[[1L]] <- quote(dbarts::bartBT)
  list(value = eval(suppliedCall, envir = callingEnv))
}

## bart's fourth formal is 'subset' where 0.9-x's fourth was 'sigest', so a
## positional 0.9-x call that carried no BayesTree-spelled name would bind a
## number to a row selector and fit a different model in silence. Refused
## for the transition release; after it, the ordinary argument rules apply.
refuseLegacyPositionalCall <- function(suppliedCall) {
  argNames <- names(suppliedCall)
  positional <- if (is.null(argNames)) {
    length(suppliedCall) - 1L
  } else {
    sum(!nzchar(argNames[-1L]))
  }
  if (positional >= 4L) {
    stop(
      "'bart' takes at most three positional arguments (formula, data, ",
      "test); ",
      positional,
      " were supplied. dbarts 0.9-x's 'bart' read the fourth as 'sigest', ",
      "where this one binds it to 'subset' - name the arguments, or call ",
      "'bartBT' for the BayesTree-style door. This refusal is removed in ",
      "dbarts ",
      tombstoneExpiry,
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## ------------------------------------------------------------------
## Names carried through '...' on the three entry points that grew one
## ------------------------------------------------------------------

## bart, dbartsControl and xbart take '...' only so that a name from this
## registry reaches a message instead of R's own "unused argument", which
## fires before any body runs. The per-door lists are keyed by name so
## foreignArgsFor can drop any entry the door has since taken as a real
## formal - a promoted name then stops being a dots name with no second
## edit here.
seedRenameReason <- paste0(
  "the engine seed is spelled 'seed'; 'rngSeed' is removed in dbarts ",
  tombstoneExpiry
)

## The family-only and feature-only formals dec-B98's consolidation moved
## onto the objects that own them. Each is still accepted for one release,
## carried on '...' and mapped onto its object after saying so once.
consolidatedArgReasons <- list(
  resid.dist = paste0(
    "the residual law is a family: write family = student(df) or ",
    "family = gaussian(); 'resid.dist' is removed in dbarts ",
    tombstoneExpiry
  ),
  dispersion = paste0(
    "the count dispersion rides its family: write ",
    "family = nbinom(dispersion = ); 'dispersion' is removed in dbarts ",
    tombstoneExpiry
  ),
  breaks = paste0(
    "the hazard period grid rides its family: write ",
    "family = hazard(breaks = ); 'breaks' is removed in dbarts ",
    tombstoneExpiry
  ),
  max.rows = paste0(
    "the hazard expansion cap rides its family: write ",
    "family = hazard(max.rows = ); 'max.rows' is removed in dbarts ",
    tombstoneExpiry
  ),
  dart = paste0(
    "a DART prior is a tree prior: write tree.prior = dart(); 'dart' is ",
    "removed in dbarts ",
    tombstoneExpiry
  ),
  levelGibbs = paste0(
    "the categorical-split level Gibbs step is declared on the tree prior: ",
    "write tree.prior = cgm(levelGibbs = ); 'levelGibbs' is removed in ",
    "dbarts ",
    tombstoneExpiry
  )
)

## Which of them each entry point used to carry.
consolidatedArgsFor <- list(
  bart = c(
    "resid.dist",
    "dispersion",
    "breaks",
    "max.rows",
    "dart",
    "levelGibbs"
  ),
  dbarts = c("resid.dist", "dispersion", "breaks", "max.rows"),
  dbartsSpec = c("resid.dist", "dispersion"),
  xbart = "dart"
)

tombstoneDotsReasons <- list(
  bart = c(
    list(rngSeed = seedRenameReason),
    consolidatedArgReasons[consolidatedArgsFor$bart]
  ),
  dbarts = consolidatedArgReasons[consolidatedArgsFor$dbarts],
  dbartsSpec = consolidatedArgReasons[consolidatedArgsFor$dbartsSpec],
  dbartsControl = list(rngSeed = seedRenameReason),
  xbart = consolidatedArgReasons[consolidatedArgsFor$xbart]
)

## Reads the consolidated names out of an entry point's '...', warning once
## per name and per entry point. The values come back under their old
## spellings; the caller maps each onto the object that now owns it, so the
## fit is the one the old spelling asked for.
resolveConsolidatedArgs <- function(matchedCall, supplied, caller, evalEnv) {
  names <- intersect(consolidatedArgsFor[[caller]], supplied)
  values <- list()
  for (name in names) {
    warnOnce(
      paste0("tombstone.consolidated.", name, ".", caller),
      "'",
      name,
      "' has left '",
      caller,
      "': ",
      consolidatedArgReasons[[name]],
      ". The value was used."
    )
    # each old spelling keeps the vocabulary it was written in: resid.dist
    # took the residual-law constructors, now the family ones, and dart took
    # the prior constructors, neither of which is exported
    env <- switch(
      name,
      resid.dist = vocabularyEnv(dbartsFamilies, evalEnv),
      dart = vocabularyEnv(dbartsPriors, evalEnv),
      evalEnv
    )
    values[name] <- list(eval(matchedCall[[name]], env))
  }
  values
}

## Maps the family-only names onto the family object they now ride. A
## 'resid.dist' of student() is the family itself, so it can only be
## reconciled with an explicit family that is already gaussian.
applyConsolidatedFamilyArgs <- function(family, consolidated) {
  residDist <- consolidated[["resid.dist"]]
  if (!is.null(residDist)) {
    if (is.function(residDist)) {
      residDist <- residDist()
    }
    if (
      !is(residDist, "dbartsFamily") ||
        residDist@token %not_in% c("gaussian", "student")
    ) {
      stop(
        "'resid.dist' takes gaussian() or student(df), and is now spelled ",
        "family = gaussian() or family = student(df)",
        call. = FALSE
      )
    }
    if (identical(residDist@token, "student")) {
      if (family@token %not_in% c("auto", "gaussian", "student")) {
        stop(
          "student residuals require a continuous gaussian response; family ",
          "\"",
          family@token,
          "\" has its own fixed error scale",
          call. = FALSE
        )
      }
      family@settings <- c(
        family@settings,
        residDist@settings[setdiff(
          names(residDist@settings),
          names(family@settings)
        )]
      )
      family@token <- "student"
    }
  }
  for (name in c("dispersion", "breaks", "max.rows")) {
    if (name %in% names(consolidated)) {
      family@settings[[name]] <- consolidated[[name]]
    }
  }
  family
}

## The names in a '...', without forcing one of them: a retired argument may
## be spelled in a vocabulary that only this package holds (resid.dist =
## student()), so its promise must not be evaluated in the caller's frame.
dotNames <- function(...) {
  count <- ...length()
  if (count == 0L) {
    return(character(0L))
  }
  supplied <- ...names()
  if (is.null(supplied)) rep_len("", count) else supplied
}

## Anything else in '...' is a caller mistake: refused by name, naming the
## entry point, rather than dropped without a word. 'supplied' is the dots
## names, "" for an unnamed one.
refuseForeignFrontDoorArgs <- function(supplied, caller, own) {
  accepted <- names(foreignArgsFor(tombstoneDotsReasons[[caller]], own))
  if (identical(caller, "bart")) {
    accepted <- c(accepted, bartBTOnlyFormals)
  }
  foreign <- setdiff(supplied[nzchar(supplied)], accepted)
  if (length(foreign) > 0L) {
    stop(
      "unused argument",
      if (length(foreign) > 1L) "s" else "",
      " ",
      paste0("'", foreign, "'", collapse = ", "),
      " passed to '",
      caller,
      "'",
      call. = FALSE
    )
  }
  unnamed <- which(!nzchar(supplied))
  if (length(unnamed) > 0L) {
    stop(
      "'",
      caller,
      "' does not take unnamed extra arguments: ",
      length(unnamed),
      " supplied",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## 'rngSeed' arrives through '...'; it is the same value under the old
## name, so it is accepted rather than refused, once per session.
resolveRenamedSeed <- function(rngSeed, caller, seed) {
  if (is.null(rngSeed)) {
    return(seed)
  }
  warnOnce(
    paste0("tombstone.rngSeed.", caller),
    "'rngSeed' is now 'seed'; the value was used. The old name is removed ",
    "in dbarts ",
    tombstoneExpiry,
    "."
  )
  rngSeed
}

## ------------------------------------------------------------------
## dbarts(sigma = ) and dbartsSpec(sigma = ), now 'sigest'
## ------------------------------------------------------------------

## The estimate supplied at creation is 'sigest' everywhere; the sampler's
## own setSigma, which sets the parameter rather than an estimate of it,
## keeps its name. Both entry points keep the old formal for the release,
## so no '...' is needed to reach this message.
resolveRenamedSigma <- function(
  sigmaIsMissing,
  sigestIsMissing,
  sigma,
  sigest,
  caller
) {
  if (sigmaIsMissing) {
    return(sigest)
  }
  if (!sigestIsMissing) {
    stop(
      "'sigma' and 'sigest' name the same estimate on '",
      caller,
      "'; supply one",
      call. = FALSE
    )
  }
  warnOnce(
    paste0("tombstone.sigma.", caller),
    "'sigma' is now 'sigest' on '",
    caller,
    "'; the value was used. The old name is removed in dbarts ",
    tombstoneExpiry,
    "."
  )
  sigma
}

## ------------------------------------------------------------------
## family = "twopart"
## ------------------------------------------------------------------

## One model, one token: the alias is gone rather than folded, so the
## message names the surviving spelling instead of quietly fitting it.
refuseTwopartFamily <- function(caller) {
  stop(
    "family = \"twopart\" is now family = \"hurdle.lognormal\" on '",
    caller,
    "'. The old token is removed in dbarts ",
    tombstoneExpiry,
    ".",
    call. = FALSE
  )
}

## ------------------------------------------------------------------
## dbartsSampler$startThreads / $stopThreads
## ------------------------------------------------------------------

## The hierarchical thread manager these drove is gone; threads are owned
## per run. Kept as no-ops so a 0.9-x Gibbs loop that brackets its sweeps
## with them still runs.
noOpThreadMethod <- function(name) {
  warnOnce(
    paste0("tombstone.", name),
    "'$",
    name,
    "' does nothing in dbarts 1.0-0: threads are owned by each run. The ",
    "method is removed in dbarts ",
    tombstoneExpiry,
    "."
  )
  invisible(NULL)
}

## ------------------------------------------------------------------
## A fit saved by dbarts 0.9-x
## ------------------------------------------------------------------

## 0.9-x wrote no state format field, so its absence is the version test:
## every state this engine writes carries one, and the bridge's own floor
## reads a missing field as encoding 0 and would blame the encoding rather
## than the release. Named here instead, at both R-side restore points, so
## a saved 0.9-x fit reaching predict through object$fit says what it is.
## 0.9-x held each chain's state as an S4 dbartsState carrying these fields;
## a class definition that no longer exists still leaves the class attribute
## on a loaded object, so either the class or the field names identifies the
## shape. Recognizing the OLD SHAPE rather than merely the missing attribute
## is what keeps an object that is no state at all - a bare list, a string -
## on the ordinary type refusal instead of being called a 0.9-x fit.
legacyStateFields <- c(
  "fit.tree",
  "fit.total",
  "sigma",
  "runningTime",
  "trees",
  "treeFits",
  "savedTrees"
)

isLegacyState <- function(state) {
  looksLegacy <- function(block) {
    inherits(block, "dbartsState") ||
      sum(legacyStateFields %in% names(block)) >= 3L
  }
  if (looksLegacy(state)) {
    return(TRUE)
  }
  if (is.list(state)) {
    for (block in state) {
      if (looksLegacy(block)) {
        return(TRUE)
      }
    }
  }
  FALSE
}

refuseLegacyState <- function(state) {
  if (!is.null(attr(state, "formatVersion")) || !isLegacyState(state)) {
    return(invisible(NULL))
  }
  stop(
    "this fit was saved by dbarts 0.9-x (no state format field); dbarts ",
    "1.0-0 cannot read it; refit with this version",
    call. = FALSE
  )
}
