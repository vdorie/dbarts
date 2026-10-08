## Tombstones: every part of the 0.9-x public surface that is gone or
## renamed but still reachable for one release. Each entry names its
## successor and the version it expires at, and the whole set expires
## together, so the release that drops them deletes this file and nothing
## survives its expiry by accident. A tombstone never adds a capability:
## it errors, or it forwards to the successor after saying so once.
## rbart_vi and its methods are the exception: they run the 0.9-x
## implementation, kept in R/rbart.R, after warning once that it is
## deprecated. That file is deleted with this one at expiry.

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
    name = "plot.rbart",
    kind = "method",
    owner = "rbart",
    successor = "stan4bart::stan4bart",
    expires = tombstoneExpiry
  ),
  list(
    name = "print.rbart",
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
    name = "sampleNodeParametersFromPrior",
    kind = "rcMethod",
    owner = "dbartsSampler",
    successor = "sampleLeafParametersFromPrior",
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
    name = "node.prior",
    kind = "argument",
    owner = "dbarts",
    successor = "leaf.prior",
    expires = tombstoneExpiry
  ),
  list(
    name = "degreesOfFreedom",
    kind = "argument",
    owner = "chi",
    successor = "df",
    expires = tombstoneExpiry
  ),
  list(
    name = "power",
    kind = "argument",
    owner = "bart",
    successor = "tree.prior = cgm(power)",
    expires = tombstoneExpiry
  ),
  list(
    name = "base",
    kind = "argument",
    owner = "bart",
    successor = "tree.prior = cgm(base)",
    expires = tombstoneExpiry
  ),
  list(
    name = "split.probs",
    kind = "argument",
    owner = "bart",
    successor = "tree.prior = cgm(split.probs)",
    expires = tombstoneExpiry
  ),
  list(
    name = "sigdf",
    kind = "argument",
    owner = "bart",
    successor = "family = gaussian(sigma = chisq(df))",
    expires = tombstoneExpiry
  ),
  list(
    name = "sigquant",
    kind = "argument",
    owner = "bart",
    successor = "family = gaussian(sigma = chisq(quant))",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.prior",
    kind = "argument",
    owner = "bart",
    successor = "family = gaussian(sigma = )",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.prior",
    kind = "argument",
    owner = "dbarts",
    successor = "family = gaussian(sigma = )",
    expires = tombstoneExpiry
  ),
  list(
    name = "resid.prior",
    kind = "argument",
    owner = "xbart",
    successor = "family = gaussian(sigma = )",
    expires = tombstoneExpiry
  ),
  list(
    name = "proposal.probs",
    kind = "argument",
    owner = "bart",
    successor = "control = dbartsControl(proposal.probs)",
    expires = tombstoneExpiry
  ),
  list(
    name = "proposal.probs",
    kind = "argument",
    owner = "dbarts",
    successor = "control = dbartsControl(proposal.probs)",
    expires = tombstoneExpiry
  ),
  list(
    name = "three-element n.burn",
    kind = "behaviour",
    owner = "xbart",
    successor = "n.burn = c(fresh, warm)",
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
    name = "thread count on $run",
    kind = "behaviour",
    owner = "dbartsSampler",
    successor = "setControl",
    expires = tombstoneExpiry
  ),
  list(
    name = "logical updateCutPoints",
    kind = "behaviour",
    owner = "dbartsSampler",
    successor = "updateCutPoints = \"position\" or \"none\"",
    expires = tombstoneExpiry
  ),
  list(
    name = "front-door startup message",
    kind = "behaviour",
    owner = ".onAttach",
    successor = NA_character_,
    expires = tombstoneExpiry
  ),
  list(
    name = "front-door defaults message",
    kind = "behaviour",
    owner = "bart",
    successor = "bartBT",
    expires = tombstoneExpiry
  ),
  list(
    name = "BayesTree-spelled pdbart call",
    kind = "behaviour",
    owner = "pdbart",
    successor = "bart's argument names",
    expires = tombstoneExpiry
  ),
  list(
    name = "BayesTree-spelled pd2bart call",
    kind = "behaviour",
    owner = "pd2bart",
    successor = "bart's argument names",
    expires = tombstoneExpiry
  ),
  list(
    name = "pdbart defaults message",
    kind = "behaviour",
    owner = "pdbart",
    successor = "bartBT",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA sigest",
    kind = "behaviour",
    owner = "bart",
    successor = "sigest = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA seed",
    kind = "behaviour",
    owner = "bart",
    successor = "seed = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA seed",
    kind = "behaviour",
    owner = "bartBT",
    successor = "seed = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA seed",
    kind = "behaviour",
    owner = "xbart",
    successor = "seed = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA updateState",
    kind = "behaviour",
    owner = "dbartsSampler",
    successor = "updateState = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA numBurnIn and numSamples",
    kind = "behaviour",
    owner = "dbartsSampler",
    successor = "$run(numBurnIn = NULL, numSamples = NULL)",
    expires = tombstoneExpiry
  ),
  list(
    name = "NA n.samples",
    kind = "behaviour",
    owner = "dbartsControl",
    successor = "n.samples = NULL",
    expires = tombstoneExpiry
  ),
  list(
    name = "sigma",
    kind = "argument",
    owner = "xbart",
    successor = "sigest",
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
    tombstoneExpiry,
    class = "dbartsDeprecatedWarning"
  )
  onceWarnState[[frontDoorDefaultsKey]] <- TRUE
  eval(matchedCall, parent.frame())
}
formals(bart2) <- formals(bart)

## ------------------------------------------------------------------
## The BayesTree-spelled bart() call
## ------------------------------------------------------------------

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
    ".",
    class = "dbartsDeprecatedWarning"
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

## A 0.9-x 'bart' call that names no BayesTree argument - bart(x, y, x.test)
## - binds the same way under both doors, so nothing forwards it and it runs
## under this door's defaults. The first such call in a session says so, as
## an informational message rather than a warning: the call is valid and its
## arguments land where they were meant. Only a call whose first argument is
## not a formula or data object qualifies, since 0.9-x's 'bart' took neither.
## A call from package code is exempt: its user never wrote 'bart', and the
## package's author is the one to act. bart2 sets the key before forwarding,
## since its own warning already says what changed.
frontDoorDefaultsKey <- "tombstone.bartDefaultsMessage"

noteFrontDoorDefaults <- function(formula, callingEnv) {
  if (is.formula(formula) || inherits(formula, "dbartsData")) {
    return(invisible(NULL))
  }
  if (isNamespace(topenv(callingEnv))) {
    return(invisible(NULL))
  }
  if (isTRUE(onceWarnState[[frontDoorDefaultsKey]])) {
    return(invisible(NULL))
  }
  onceWarnState[[frontDoorDefaultsKey]] <- TRUE
  # built by hand: messageCondition() is newer than the R this package
  # supports
  message(structure(
    class = c(
      "dbartsFrontDoorMessage",
      "dbartsMessage",
      "message",
      "condition"
    ),
    list(
      message = paste0(
        "dbarts: 'bart' is the function 0.9-x called 'bart2', with its own ",
        "defaults (75 trees; four chains, their draws merged) rather than ",
        "those of 0.9-x's 'bart' (200 trees, one chain). Call 'bartBT' for ",
        "the BayesTree-style fit and its defaults. Shown once per session ",
        "until dbarts ",
        tombstoneExpiry,
        ".\n"
      ),
      call = NULL
    )
  ))
  invisible(NULL)
}

## ------------------------------------------------------------------
## pdbart and pd2bart calls in BayesTree's spellings
## ------------------------------------------------------------------

## pdbart and pd2bart fit through bart. Every name bartBT takes and bart does
## not is translated to the bart argument it lands on, after a warning per
## name, except the two pdbart sets itself (x.test, sampleronly), which are
## refused with bart's own spellings of them. The table is pinned by
## inst/tinytest/test-pdbart.R.
pdbartBayesTreeNames <- c(
  x.train = "formula",
  y.train = "data",
  sigdf = "sigdf",
  sigquant = "sigquant",
  power = "tree.prior",
  base = "tree.prior",
  splitprobs = "tree.prior",
  binaryOffset = "offset",
  ntree = "n.trees",
  ndpost = "n.samples",
  nskip = "n.burn",
  printevery = "printEvery",
  keepevery = "n.thin",
  keeptrainfits = "keepTrainingFits",
  usequants = "useQuantiles",
  numcut = "n.cuts",
  printcutoffs = "printCutoffs",
  nchain = "n.chains",
  nthread = "n.threads",
  combinechains = "combineChains",
  keeptrees = "keepTrees",
  keepcall = "keepCall",
  proposalprobs = "control",
  keepsampler = "keepSampler"
)

## The spelling a translation warning names, where it is not the bare name.
pdbartBayesTreeSpelling <- c(
  sigdf = "family = gaussian(sigma = chisq(df = ))",
  sigquant = "family = gaussian(sigma = chisq(quant = ))",
  power = "tree.prior = cgm(power = )",
  base = "tree.prior = cgm(base = )",
  splitprobs = "tree.prior = cgm(split.probs = )",
  proposalprobs = "control = dbartsControl(proposal.probs = )"
)

## Rewrites a pdbart call's BayesTree spellings into bart's. Warns once per
## name and function, package callers included, since their authors are the
## ones to change the call. sigdf and sigquant stay under their own names,
## which bart still reads: the residual prior rides a family that cannot be
## built before the response is known. A setting given in both spellings is
## refused naming both. 'fits' is FALSE where nothing is fit, a fit or sampler
## having been passed in, and the warning then says nothing of defaults.
translatePdbartCall <- function(call, legacy, callingEnv, caller, fits = TRUE) {
  argNames <- names(call)[-1L]
  for (old in legacy) {
    new <- pdbartBayesTreeNames[[old]]
    if (new != old && new != "control" && new %in% argNames) {
      stop(
        "'",
        old,
        "' and '",
        new,
        "' both set one setting of '",
        caller,
        "'; give '",
        new,
        "' alone",
        call. = FALSE
      )
    }
  }
  for (old in legacy) {
    spelling <- if (old %in% names(pdbartBayesTreeSpelling)) {
      pdbartBayesTreeSpelling[[old]]
    } else {
      pdbartBayesTreeNames[[old]]
    }
    warnOnce(
      paste0("tombstone.", caller, ".", old),
      "'",
      old,
      "' is BayesTree's spelling; '",
      caller,
      if (fits) {
        "' now fits through 'bart', which takes '"
      } else {
        "' takes bart's spelling, '"
      },
      spelling,
      "', and the value was used. ",
      if (fits) {
        paste0(
          "Settings not named take bart's defaults, so the model is not the ",
          "one dbarts 0.9-34 fit. "
        )
      },
      "BayesTree spellings are refused from dbarts ",
      tombstoneExpiry,
      ".",
      class = "dbartsDeprecatedWarning"
    )
  }

  priorNames <- c(power = "power", base = "base", splitprobs = "split.probs")
  priorArgs <- list()
  for (old in intersect(names(priorNames), legacy)) {
    priorArgs[priorNames[[old]]] <- list(call[[old]])
    call[[old]] <- NULL
  }
  if (length(priorArgs) > 0L) {
    call$tree.prior <- as.call(c(list(quote(cgm)), priorArgs))
  }

  # the mixture is set on a copy of the caller's control, so the control's
  # other settings stand as bart would read them
  if ("proposalprobs" %in% legacy) {
    probs <- eval(call$proposalprobs, callingEnv)
    call$proposalprobs <- NULL
    if ("control" %in% argNames) {
      control <- eval(call$control, callingEnv)
      control@proposal.probs <- dbartsControl(
        proposal.probs = probs
      )@proposal.probs
      validObject(control)
      call$control <- control
    } else {
      controlCall <- quote(dbarts::dbartsControl(proposal.probs = NULL))
      controlCall["proposal.probs"] <- list(probs)
      call$control <- controlCall
    }
  }

  renamed <- setdiff(legacy, c(names(priorNames), "proposalprobs"))
  callNames <- names(call)
  callNames[callNames %in% renamed] <- pdbartBayesTreeNames[
    callNames[callNames %in% renamed]
  ]
  names(call) <- callNames
  call
}

## A pdbart or pd2bart call that passes data and names no BayesTree argument
## fits under bart's defaults where 0.9-34 fit under BayesTree's. The first
## such call in a session says so; a call from package code is exempt. One
## key serves both functions.
pdbartDefaultsKey <- "tombstone.pdbartDefaultsMessage"

notePdbartDefaults <- function(callingEnv) {
  if (isNamespace(topenv(callingEnv))) {
    return(invisible(NULL))
  }
  if (isTRUE(onceWarnState[[pdbartDefaultsKey]])) {
    return(invisible(NULL))
  }
  onceWarnState[[pdbartDefaultsKey]] <- TRUE
  # built by hand: messageCondition() is newer than the R this package
  # supports
  message(structure(
    class = c(
      "dbartsFrontDoorMessage",
      "dbartsMessage",
      "message",
      "condition"
    ),
    list(
      message = paste0(
        "dbarts: 'pdbart' and 'pd2bart' fit through 'bart', with its ",
        "defaults (75 trees; four chains, their draws merged) rather than ",
        "those of 0.9-x (200 trees, one chain). ",
        "pdbart(bartBT(x, y, keeptrees = TRUE)) gives the 0.9-34 model. ",
        "Shown once per session until dbarts ",
        tombstoneExpiry,
        ".\n"
      ),
      call = NULL
    )
  ))
  invisible(NULL)
}

## The once-per-session notices bart shows that pdbart replaces with its
## own: the defaults message and the warnings for sigdf and sigquant, which
## pdbart translates. Held back for the duration of 'expr', each key
## restored as it was found; every other retired name still warns through
## bart.
holdingBartNotices <- function(expr) {
  keys <- c(
    frontDoorDefaultsKey,
    paste0("tombstone.consolidated.", c("sigdf", "sigquant"), ".bart")
  )
  saved <- mget(
    keys,
    envir = onceWarnState,
    ifnotfound = rep_len(list(NULL), length(keys))
  )
  on.exit(
    for (key in keys) {
      onceWarnState[[key]] <- saved[[key]]
    }
  )
  for (key in keys) {
    onceWarnState[[key]] <- TRUE
  }
  expr
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
  power = paste0(
    "the branching decay is a tree prior: write tree.prior = cgm(power = ) ",
    "or tree.prior = dart(power = ); 'power' is removed in dbarts ",
    tombstoneExpiry
  ),
  base = paste0(
    "the branching base is a tree prior: write tree.prior = cgm(base = ) or ",
    "tree.prior = dart(base = ); 'base' is removed in dbarts ",
    tombstoneExpiry
  ),
  split.probs = paste0(
    "the per-predictor split probabilities are a tree prior: write ",
    "tree.prior = cgm(split.probs = ); 'split.probs' is removed in dbarts ",
    tombstoneExpiry
  ),
  resid.prior = paste0(
    "the residual scale's prior rides its family: write ",
    "family = gaussian(sigma = chisq(df, quant)) or ",
    "family = gaussian(sigma = fixed(value)); 'resid.prior' is removed in ",
    "dbarts ",
    tombstoneExpiry
  ),
  sigdf = paste0(
    "the residual prior's degrees of freedom rides its family: write ",
    "family = gaussian(sigma = chisq(df = )); 'sigdf' is removed in dbarts ",
    tombstoneExpiry
  ),
  sigquant = paste0(
    "the residual prior's quantile rides its family: write ",
    "family = gaussian(sigma = chisq(quant = )); 'sigquant' is removed in ",
    "dbarts ",
    tombstoneExpiry
  ),
  proposal.probs = paste0(
    "the tree-move mixture is a sampler setting: write ",
    "control = dbartsControl(proposal.probs = ); 'proposal.probs' is ",
    "removed in dbarts ",
    tombstoneExpiry
  )
)

## The prior scalars dec-B116 moved onto the prior objects and the control.
## Named as a set because bart's hurdle arc has to hand them back to each
## component call, where the family-only names never travel.
consolidatedPriorScalars <- c(
  "power",
  "base",
  "split.probs",
  "resid.prior",
  "sigdf",
  "sigquant",
  "proposal.probs"
)

## The consolidated names whose value must reach its object UNEVALUATED: the
## split probabilities are written in a vocabulary only the prior resolver
## holds (num.vars, numvars), which does not resolve in the caller's frame.
unevaluatedConsolidatedArgs <- "split.probs"

## Which of them each entry point used to carry.
consolidatedArgsFor <- list(
  bart = c(
    "power",
    "base",
    "split.probs",
    "resid.prior",
    "sigdf",
    "sigquant",
    "proposal.probs"
  ),
  dbarts = c("resid.prior", "proposal.probs"),
  xbart = c("resid.prior")
)

tombstoneDotsReasons <- list(
  bart = c(
    list(rngSeed = seedRenameReason),
    consolidatedArgReasons[consolidatedArgsFor$bart]
  ),
  dbarts = consolidatedArgReasons[consolidatedArgsFor$dbarts],
  dbartsControl = list(rngSeed = seedRenameReason),
  xbart = consolidatedArgReasons[consolidatedArgsFor$xbart]
)

## 0.9-x's xbart read a third n.burn element as a per-replication burn-in.
## Chains are never carried between replications now, so the element names
## nothing; refused by name rather than dropped in silence.
refuseThreeElementBurn <- function(n.burn) {
  if (length(n.burn) > 2L) {
    stop(
      "'n.burn' must be of length 1 or 2: the burn-in of a freshly started ",
      "chain and the burn-in of a warm start onto the same data split. ",
      "dbarts 0.9-x read a third element as a per-replication burn-in, and ",
      "a chain is never carried between replications now. The three-element ",
      "form is removed in dbarts ",
      tombstoneExpiry,
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## The names that select the legacy door: a formal of bartBT that bart
## neither takes as a formal nor carries on its own '...'. Derived, so a name
## added to either signature cannot leave the shim behind; the list is pinned
## by inst/tinytest/test-tombstones.R.
bartBTOnlyFormals <- setdiff(
  names(formals(bartBT)),
  c(names(formals(bart)), consolidatedArgsFor$bart)
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
      ". The value was used.",
      class = "dbartsDeprecatedWarning"
    )
    # resid.prior alone is still written in a vocabulary this package holds:
    # the prior constructors, which are not exported
    values[name] <- list(
      if (name %in% unevaluatedConsolidatedArgs) {
        forwardedSplitProbs(matchedCall[[name]], evalEnv)
      } else if (name == "resid.prior") {
        evalInVocabulary(matchedCall[[name]], dbartsPriors, evalEnv)
      } else {
        eval(matchedCall[[name]], evalEnv)
      }
    )
  }
  values
}

## The retired 'split.probs', left for the prior resolver, whose vocabulary
## (num.vars, numvars) it may be written in. A direct expression stays as
## written. A forwarded reference (..N) is forced first and, when that
## succeeds, stands as its value; when it fails, the recovered expression
## goes on as a call that evaluates it in the environment it was written in,
## layered under the vocabulary the resolver supplies where the call runs.
forwardedSplitProbs <- function(expr, evalEnv) {
  if (!isDotsReference(expr)) {
    return(expr)
  }
  # an expression that fails when forced runs a second time below, the cost
  # evalInVocabulary also pays
  forced <- tryCatch(list(eval(expr, evalEnv)), error = function(e) NULL)
  if (!is.null(forced)) {
    return(forced[[1L]])
  }
  written <- recoverForwardedArgument(expr, evalEnv)
  if (isDotsReference(written$expr)) {
    return(expr)
  }
  as.call(list(
    evalInVocabularyOver,
    call("quote", written$expr),
    written$env,
    quote(environment())
  ))
}

## Evaluates 'expr' in 'writtenEnv' with the bindings of the vocabulary
## environment 'vocabEnv' laid over it.
evalInVocabularyOver <- function(expr, writtenEnv, vocabEnv) {
  eval(expr, vocabularyEnv(as.list(vocabEnv, all.names = TRUE), writtenEnv))
}

## The retired flat 'resid.prior', resolved to an object: a bare constructor
## name means its defaults, and anything that is not a residual prior is
## refused here, where the spelling the caller wrote is still known. An
## explicit NULL is a supplied value, not an absent one, so it is refused
## rather than read as the default.
consolidatedResidPrior <- function(consolidated) {
  if ("resid.prior" %not_in% names(consolidated)) {
    return(NULL)
  }
  value <- consolidated[["resid.prior"]]
  if (is.function(value)) {
    value <- value()
  }
  if (!is(value, "dbartsResidPrior")) {
    stop(
      "'resid.prior' must be a residual prior specification; see ?dbartsPriors",
      call. = FALSE
    )
  }
  value
}

## The residual prior written twice: a retired flat spelling and a family
## object whose own call named 'sigma'. Agreeing spellings are one statement
## said twice and stand; disagreeing ones are a conflict no precedence rule
## can settle, so the call is refused naming both and saying which to delete.
## Priors agree when they read back as the same constructor call, which is
## also how the message renders them, so the test and the message cannot
## disagree. 'flatName' is every retired spelling the caller actually wrote,
## named together in one message: bart's retired sigdf/sigquant build the
## same prior under two names at once, and naming only one would have the
## caller delete it, rerun, and hit the same refusal on the other.
reconcileResidPrior <- function(flat, flatName, family) {
  familySigma <- familySetting(family, "sigma", NULL)
  if (is.null(flat) || is.null(familySigma)) {
    if (is.null(flat)) familySigma else flat
  } else if (identical(formatResidPrior(flat), formatResidPrior(familySigma))) {
    flat
  } else {
    names <- paste0("'", flatName, "'", collapse = " and ")
    plural <- length(flatName) > 1L
    stop(
      names,
      " and the family's own 'sigma' set different residual priors: ",
      names,
      if (plural) " together say " else " says ",
      formatResidPrior(flat),
      ", family = ",
      family@token,
      "(sigma = ) says ",
      formatResidPrior(familySigma),
      ". Delete ",
      names,
      if (plural) {
        ", both removed in dbarts "
      } else {
        ", which is removed in dbarts "
      },
      tombstoneExpiry,
      ", and keep the prior on the family.",
      call. = FALSE
    )
  }
}

## The names in a '...', without forcing one of them: a retired argument may
## be spelled in a vocabulary that only this package holds (resid.prior =
## chisq(3, 0.9)), so its promise must not be evaluated in the caller's frame.
dotNames <- function(...) {
  count <- ...length()
  if (count == 0L) {
    return(character(0L))
  }
  supplied <- ...names()
  if (is.null(supplied)) rep_len("", count) else supplied
}

## A name refused from a front door's '...' whose successor is not the
## obvious next guess gets a hint appended to the "unused argument" message.
## This is a refusal, not a tombstone: the name never worked on this door,
## so there is no registry entry, no warning and no NEWS text. A hint keyed
## by door wins over the shared one.
foreignArgHints <- list(
  node.prior = "the leaf prior is 'leaf.prior'"
)
foreignArgHintsByDoor <- list(
  dbartsSpec = list(
    sigma = "the starting estimate of sigma is 'sigest'",
    resid.prior = paste0(
      "the residual prior rides its family: write ",
      "family = gaussian(sigma = )"
    )
  )
)

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
    known <- c(foreignArgHintsByDoor[[caller]], foreignArgHints)
    hints <- unlist(
      known[unique(foreign[foreign %in% names(known)])],
      use.names = FALSE
    )
    stop(
      "unused argument",
      if (length(foreign) > 1L) "s" else "",
      " ",
      paste0("'", foreign, "'", collapse = ", "),
      " passed to '",
      caller,
      "'",
      if (length(hints) > 0L) {
        paste0("; ", paste(hints, collapse = "; "))
      } else {
        ""
      },
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
    ".",
    class = "dbartsDeprecatedWarning"
  )
  # 0.9-x's rngSeed = NA meant no seed; the warning above covers the name
  refuseNaN(rngSeed, "rngSeed")
  if (isSingleNA(rngSeed)) NULL else rngSeed
}

## ------------------------------------------------------------------
## dbarts(sigma = ), now 'sigest'
## ------------------------------------------------------------------

## The estimate supplied at creation is 'sigest' everywhere; the sampler's
## own setSigma, which sets the parameter rather than an estimate of it,
## keeps its name. The entry points keep the old formal for the release,
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
    ".",
    class = "dbartsDeprecatedWarning"
  )
  sigma
}

## ------------------------------------------------------------------
## NA where NULL now means "not given"
## ------------------------------------------------------------------

## A 0.9-x entry point that took NA for an absent value keeps reading it that
## way for the release, after saying so once per entry point and argument;
## the one message serves every such site, so the wording cannot drift.
warnNAForNull <- function(argument, caller) {
  warnOnce(
    paste0("tombstone.NA.", argument, ".", caller),
    "'",
    argument,
    " = NA' is now '",
    argument,
    " = NULL' on '",
    caller,
    "'; the NA was read as NULL. NA is a missing value, and reading it as ",
    "absent is removed in dbarts ",
    tombstoneExpiry,
    ".",
    class = "dbartsDeprecatedWarning"
  )
}

## An argument that never shipped with NA for absent refuses it now.
refuseNAForNull <- function(argument, caller, meaning = NULL) {
  stop(
    "'",
    argument,
    "' must not be NA on '",
    caller,
    "': NA is a missing value; use NULL instead",
    if (is.null(meaning)) "" else paste0(", which gives ", meaning),
    call. = FALSE
  )
}

## NaN is never "absent": it is an invalid value, refused by name.
refuseNaN <- function(x, argument) {
  if (is.numeric(x) && length(x) == 1L && is.nan(x)) {
    stop(
      "'",
      argument,
      "' must be a number or NULL; NaN is neither",
      call. = FALSE
    )
  }
  invisible(x)
}

## TRUE for a single NA of any type, the value a caller writes for "missing".
## NaN is not one: it carries no such intent, and the argument's own checks
## refuse it.
isSingleNA <- function(x) {
  is.atomic(x) && length(x) == 1L && is.na(x) && !is.nan(x)
}

## The residual-scale estimate supplied at creation: NULL means estimate it
## by least squares, and NA_real_ is the internal spelling of the same thing
## that the data object's slot keeps. 'onNA' says what an explicit NA is: a
## "warn" where 0.9-34 took it, a "refuse" where the argument never shipped,
## and "silent" for a value that arrived under a retired name whose own
## warning has already spoken. NaN is refused in every case.
resolveSigestArg <- function(
  sigest,
  caller,
  onNA = c("warn", "refuse", "silent"),
  name = "sigest"
) {
  onNA <- match.arg(onNA)
  if (is.null(sigest)) {
    return(NA_real_)
  }
  refuseNaN(sigest, name)
  if (isSingleNA(sigest)) {
    switch(
      onNA,
      warn = warnNAForNull(name, caller),
      refuse = refuseNAForNull(name, caller),
      silent = NULL
    )
    return(NA_real_)
  }
  sigest
}

## ------------------------------------------------------------------
## dbarts(node.prior = ), now 'leaf.prior'
## ------------------------------------------------------------------

## The prior vocabulary is NSE, so this never forces the argument: only
## presence is tested, via the caller's own missing() on both formals, read
## before either is assigned. dbarts keeps the old formal for the release,
## so no '...' is needed to reach this message.
resolveRenamedLeafPrior <- function(
  matchedCall,
  nodePriorSupplied,
  leafPriorSupplied,
  caller
) {
  if (!nodePriorSupplied) {
    return(matchedCall)
  }
  if (leafPriorSupplied) {
    stop(
      "'node.prior' and 'leaf.prior' name the same prior on '",
      caller,
      "'; supply one",
      call. = FALSE
    )
  }
  warnOnce(
    paste0("tombstone.node.prior.", caller),
    "'node.prior' is now 'leaf.prior' on '",
    caller,
    "'; the value was used. The old name is removed in dbarts ",
    tombstoneExpiry,
    ".",
    class = "dbartsDeprecatedWarning"
  )
  # a plain $<- assignment of NULL deletes the element instead of setting
  # it, and node.prior = NULL is a supplied value here (nodePriorSupplied
  # is TRUE), not an absent one - the wrapping list() keeps it
  matchedCall["leaf.prior"] <- list(matchedCall[["node.prior"]])
  matchedCall$node.prior <- NULL
  matchedCall
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
    "' does nothing in dbarts 1.0-0: threads are owned by each run, which ",
    "uses control@n.threads; $setControl changes it. The method is removed ",
    "in dbarts ",
    tombstoneExpiry,
    ".",
    class = "dbartsDeprecatedWarning"
  )
  invisible(NULL)
}

## ------------------------------------------------------------------
## dbartsSampler$run's thread count
## ------------------------------------------------------------------

## 0.9-x's run took a per-call thread count as its fourth argument, a formal
## spelled numThreads and documented as n.threads. A run now uses the
## sampler's own count, which setControl changes, and the draws do not depend
## on it, so a value passed either way - by name, or as the one unnamed
## argument after updateState - is ignored after saying so once. Anything
## else in the dots is refused, as an unused argument would be.
ignoreRunThreadCount <- function(...) {
  supplied <- dotNames(...)
  if (length(supplied) == 0L) {
    return(invisible(NULL))
  }
  legacy <- supplied %in%
    c("n.threads", "numThreads") |
    (!nzchar(supplied) & seq_along(supplied) == 1L)
  if (!all(legacy)) {
    foreign <- supplied[!legacy]
    foreign[!nzchar(foreign)] <- "<unnamed>"
    stop(
      "unused argument",
      if (length(foreign) > 1L) "s" else "",
      " ",
      paste0("'", foreign, "'", collapse = ", "),
      " passed to '$run'",
      call. = FALSE
    )
  }
  warnOnce(
    "tombstone.run.n.threads",
    "'$run' no longer takes a thread count; the value was ignored. A run ",
    "uses the sampler's own, control@n.threads, which $setControl changes, ",
    "and the draws do not depend on it. The argument is removed in dbarts ",
    tombstoneExpiry,
    ".",
    class = "dbartsDeprecatedWarning"
  )
  invisible(NULL)
}

## ------------------------------------------------------------------
## dbartsSampler$setPredictor's logical updateCutPoints
## ------------------------------------------------------------------

## 0.9-x's updateCutPoints was a logical, TRUE re-deriving the grid with
## every split left on its position. It is now one of three words; a logical
## is read as the word it stood for after saying so once, and anything else
## is returned as it came, for the caller to match or refuse.
resolveUpdateCutPoints <- function(updateCutPoints) {
  if (
    !is.logical(updateCutPoints) ||
      length(updateCutPoints) != 1L ||
      is.na(updateCutPoints)
  ) {
    return(updateCutPoints)
  }
  warnOnce(
    "tombstone.updateCutPoints.logical",
    "'updateCutPoints' is now one of \"none\", \"position\" or \"value\"; ",
    "TRUE was taken as \"position\" and FALSE as \"none\". A logical is no ",
    "longer taken in dbarts ",
    tombstoneExpiry,
    ".",
    class = "dbartsDeprecatedWarning"
  )
  if (updateCutPoints) "position" else "none"
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
