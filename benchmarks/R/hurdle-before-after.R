#!/usr/bin/env Rscript

# The hurdle before/after gate: every hurdle.lognormal call a build takes
# today must give the same draws and the same output from every method after
# the hurdle completion work as before it. Record on each build, then compare:
#
#   R_LIBS=<absolute library path> Rscript benchmarks/R/hurdle-before-after.R \
#     record <out.rds>
#   Rscript benchmarks/R/hurdle-before-after.R compare <a.rds> <b.rds>
#
# Both records are made on one machine (two threads at most), one from a
# library built from the base commit and one from the branch, both in the
# shipped mode (buildInfo()$mode is checked and recorded). compare reports each
# element by identical() and exits with status 1 on any difference.
#
# A scenario is one fit call; it records the fit less its samplers (each part's
# trees are read through extract(type = "trees") instead, and a part's own call,
# which holds the data it was handed, is left out) and the output of every
# method on it. An error or warning is an element too, so a refusal's words
# are held.

args <- commandArgs(trailingOnly = TRUE)

# --- data --------------------------------------------------------------------

semicontinuous <- function(n, seed, p = 3L) {
  set.seed(seed)
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, letters[seq_len(p)]))
  occupied <- rbinom(n, 1L, pnorm(0.8 * x[, 1L] - 0.4)) == 1L
  y <- numeric(n)
  y[occupied] <- exp(rnorm(sum(occupied), 1 + 0.5 * x[, 2L][occupied], 0.4))
  # at least one zero and one positive, whatever the draw gave
  y[1L] <- 0
  y[2L] <- 1.5
  list(x = x, y = y)
}

withFactors <- function(n, seed) {
  d <- semicontinuous(n, seed, 2L)
  x <- as.data.frame(d$x)
  x$f <- factor(sample(c("u", "v", "w"), n, TRUE))
  x$g <- factor(sample(c("p", "q"), n, TRUE), levels = c("p", "q", "r"))
  list(x = x, y = d$y)
}

# a level that occurs only on zero rows, and one on no row at all
zeroRowsLevel <- function(n, seed) {
  d <- withFactors(n, seed)
  x <- d$x
  x$f <- as.character(x$f)
  x$f[d$y == 0][1:3] <- "z"
  x$f[d$y > 0][x$f[d$y > 0] == "z"] <- "u"
  x$f <- factor(x$f)
  list(x = x, y = d$y)
}

missingAt <- function(d, rows, columns) {
  for (j in columns) {
    d$x[rows[seq_len(min(4L, length(rows)))], j] <- NA
  }
  d
}

zeroOnly <- function(d) which(d$y == 0)[-1L]
positiveOnly <- function(d) which(d$y > 0)[-1L]

# --- the scenarios -----------------------------------------------------------

n <- 70L
base <- list(
  n.trees = 8L,
  n.burn = 10L,
  n.samples = 12L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 20261103L
)

scenarioList <- function() {
  num <- semicontinuous(n, 1L)
  noNames <- num
  colnames(noNames$x) <- NULL
  named <- num
  rownames(named$x) <- paste0("r", seq_len(n))
  fac <- withFactors(n, 2L)
  framed <- fac
  rownames(framed$x) <- paste0("s", seq_len(n))
  zl <- zeroRowsLevel(n, 3L)
  mixed <- fac
  mixed$x$h <- sample(c(TRUE, FALSE), n, TRUE)
  mixed$x$s <- sample(c("m", "n", "o"), n, TRUE)
  naZero <- missingAt(num, zeroOnly(num), 1L)
  naPositive <- missingAt(num, positiveOnly(num), 2L)
  naBoth <- missingAt(
    missingAt(num, zeroOnly(num), 3L),
    positiveOnly(num),
    3L
  )
  naFactor <- missingAt(fac, positiveOnly(fac), "f")
  naY <- num
  naY$y[c(5L, 9L)] <- NA
  naYFrame <- fac
  naYFrame$y[c(5L, 9L)] <- NA

  with <- function(d, ...) c(list(x = d$x, y = d$y), list(...))
  list(
    numeric = with(num),
    numeric_nocolnames = with(noNames),
    numeric_rownames = with(named),
    frame_factors = with(fac),
    frame_rownames = with(framed),
    frame_zero_rows_level = with(zl),
    frame_empty_level = with(fac),
    na_omit_zero_rows = with(naZero, na.action = "na.omit"),
    na_exclude_zero_rows = with(naZero, na.action = "na.exclude"),
    na_default_zero_rows = with(naZero),
    na_pass_zero_rows = with(naZero, na.action = "na.pass"),
    na_omit_positive_rows = with(naPositive, na.action = "na.omit"),
    na_exclude_positive_rows = with(naPositive, na.action = "na.exclude"),
    na_default_positive_rows = with(naPositive),
    na_both = with(naBoth),
    na_omit_both = with(naBoth, na.action = "na.omit"),
    na_exclude_factor = with(naFactor, na.action = "na.exclude"),
    na_default_factor = with(naFactor),
    na_y_omit = with(naY, na.action = "na.omit"),
    na_y_exclude = with(naY, na.action = "na.exclude"),
    na_y_default = with(naY),
    na_y_exclude_frame = with(naYFrame, na.action = "na.exclude"),
    chains3_threads2 = with(num, n.chains = 3L, n.threads = 2L),
    threads2 = with(fac, n.threads = 2L),
    unseeded = with(num, seed = NA_integer_, unseeded = TRUE),
    sigdf = with(num, sigdf = 5),
    sigquant = with(num, sigquant = 0.95),
    resid_prior = with(num, resid.prior = quote(chisq(4, 0.8))),
    proposal_probs = with(
      num,
      proposal.probs = c(birth_death = 0.6, swap = 0.1, change = 0.3)
    ),
    sigma_chisq = with(
      num,
      family = quote(hurdle.lognormal(sigma = chisq(4, 0.8)))
    ),
    sigma_fixed = with(
      num,
      family = quote(hurdle.lognormal(sigma = fixed(0.5)))
    ),
    keep_trees_false = with(num, keepTrees = FALSE),
    leaf_prior_k = with(num, k = 1.5),
    power_base = with(num, power = 1.5, base = 0.9),
    n_cuts = with(num, n.cuts = 12L, useQuantiles = TRUE),
    factors_indicators = with(fac, factors = "indicators"),
    keep_call_false = with(num, keepCall = FALSE),
    storage_single = with(num, storage = "single"),
    thinned = with(num, n.thin = 2L),
    mixed_columns = with(mixed),
    one_column = with(list(x = num$x[, 1L, drop = FALSE], y = num$y)),
    negative_y = with(list(x = num$x, y = replace(num$y, 3L, -1))),
    infinite_y = with(list(x = num$x, y = replace(num$y, 3L, Inf))),
    all_positive = with(list(x = num$x, y = abs(num$y) + 1)),
    all_zero = with(list(x = num$x, y = num$y * 0)),
    short_y = with(list(x = num$x, y = num$y[-1L])),
    variance_refused = with(num, variance = "hetero"),
    keep_fits_refused = with(num, keepFits = FALSE),
    sampler_only_refused = with(num, samplerOnly = TRUE),
    leaf_prior_sd_refused = with(num, leaf.prior = quote(normal(sd = 0.5)))
  )
}

# --- recording ---------------------------------------------------------------

# a value, the warnings it raised, or the error that stopped it
guarded <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "recordedError"
      )
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

# the same value, without anything that cannot be compared across processes
plain <- function(value) {
  if (is.language(value)) {
    return(paste(deparse(value), collapse = "\n"))
  }
  value
}

stripFit <- function(fit) {
  fit$call <- plain(fit$call)
  for (part in c("zero", "positive")) {
    fit[[part]][["call"]] <- NULL
    fit[[part]][["fit"]] <- NULL
  }
  fit
}

seeded <- function(seed, expr) {
  set.seed(seed)
  guarded(expr)
}

recordFit <- function(fit, x, y, name) {
  isFrame <- is.data.frame(x)
  out <- list()
  put <- function(label, value) out[[label]] <<- value

  put("fit", stripFit(fit))
  for (part in c("zero", "positive")) {
    put(
      paste0("trees.", part),
      guarded(extract(fit[[part]], type = "trees"))
    )
  }
  put("family", guarded(plain(family(fit))))

  types <- c("ev", "ppd", "prob", "bart")
  for (type in types) {
    put(
      paste0("extract.", type),
      seeded(11L, extract(fit, type = type))
    )
    put(
      paste0("extract.unchained.", type),
      seeded(12L, extract(fit, type = type, combineChains = FALSE))
    )
    put(
      paste0("fitted.", type),
      seeded(13L, fitted(fit, type = type))
    )
    put(
      paste0("fitted.ci.", type),
      seeded(14L, fitted(fit, type = type, ci.level = 0.8))
    )
  }
  for (type in c("loglik", "sigma", "k", "leaf.prior.sd", "varcount")) {
    put(paste0("extract.", type), guarded(extract(fit, type = type)))
  }
  put(
    "extract.loglik.unchained",
    guarded(
      extract(fit, type = "loglik", combineChains = FALSE)
    )
  )
  put("residuals", guarded(residuals(fit)))
  put("residuals.prob", guarded(residuals(fit, type = "prob")))

  complete <- which(stats::complete.cases(x))
  rows <- complete[unique(pmin(c(1L, 2L, 5L, 11L, 30L, 60L), length(complete)))]
  newdata <- x[rows, , drop = FALSE]
  withGaps <- x[
    unique(c(rows, which(!stats::complete.cases(x))[1:2])),
    ,
    drop = FALSE
  ]
  for (type in types) {
    put(
      paste0("predict.", type),
      seeded(15L, predict(fit, newdata, type = type))
    )
    put(
      paste0("predict.ci.", type),
      seeded(16L, predict(fit, newdata, type = type, ci.level = 0.9))
    )
    put(
      paste0("predict.train.", type),
      seeded(17L, predict(fit, type = type))
    )
  }
  put(
    "predict.unchained",
    seeded(
      18L,
      predict(fit, newdata, type = "ev", combineChains = FALSE)
    )
  )
  put(
    "predict.gaps.pass",
    seeded(
      20L,
      predict(fit, withGaps, type = "ev", na.action = na.pass)
    )
  )
  put(
    "predict.gaps.default",
    seeded(
      20L,
      predict(fit, withGaps, type = "ev")
    )
  )
  put(
    "predict.na.omit",
    seeded(
      19L,
      predict(fit, newdata, type = "ev", na.action = na.omit)
    )
  )
  put(
    "predict.na.exclude",
    seeded(
      19L,
      predict(fit, newdata, type = "prob", na.action = na.exclude)
    )
  )

  put("print", guarded(utils::capture.output(print(fit))))
  put("summary", guarded(dropSummaryCalls(summary(fit))))
  put("summary.print", guarded(utils::capture.output(print(summary(fit)))))

  numericColumns <- if (isFrame) {
    which(vapply(x, is.numeric, NA))
  } else {
    seq_len(ncol(x))
  }
  if (length(numericColumns) >= 2L) {
    first <- numericColumns[1L]
    second <- numericColumns[2L]
    one <- guarded(dbarts::pdbart(
      fit,
      xind = first,
      levs = list(c(0.3, 0.7)),
      pl = FALSE
    ))
    put("pdbart", guarded(dropPd(one$value)))
    two <- guarded(dbarts::pd2bart(
      fit,
      xind = c(first, second),
      levs = list(c(0.3, 0.7), c(0.4, 0.6)),
      pl = FALSE
    ))
    put("pd2bart", guarded(dropPd(two$value)))
    put("pdbart.warnings", one$warnings)
  }

  put(
    "plot",
    guarded({
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off())
      plot(fit)
      TRUE
    })
  )
  out
}

# a part's own call holds the data the part was handed
dropSummaryCalls <- function(summary) {
  for (part in c("zero", "positive")) {
    summary[[part]][["call"]] <- NULL
  }
  summary
}

dropPd <- function(pd) {
  if (inherits(pd, "recordedError")) {
    return(pd)
  }
  pd$times <- NULL
  pd$fit <- NULL
  pd$bartcall <- plain(pd$bartcall)
  pd
}

fitOne <- function(spec) {
  unseeded <- isTRUE(spec$unseeded)
  spec$unseeded <- NULL
  x <- spec$x
  y <- spec$y
  spec$x <- NULL
  spec$y <- NULL
  callArgs <- base
  callArgs[names(spec)] <- spec
  if (unseeded) {
    callArgs$seed <- NULL
  }
  if (!"family" %in% names(callArgs)) {
    callArgs$family <- "hurdle.lognormal"
  }
  if (unseeded) {
    set.seed(4242L)
  }
  made <- guarded(do.call("bart", c(list(x, y), callArgs)))
  list(x = x, y = y, made = made)
}

recordAll <- function() {
  info <- dbarts:::buildInfo()
  if (!identical(info$mode, "shipped")) {
    stop("the library is not a shipped build: mode ", info$mode)
  }
  scenarios <- scenarioList()
  records <- list()
  for (name in names(scenarios)) {
    one <- fitOne(scenarios[[name]])
    if (inherits(one$made$value, "recordedError")) {
      records[[name]] <- list(fit = one$made)
      next
    }
    fit <- one$made$value
    records[[name]] <- c(
      list(warnings = one$made$warnings),
      recordFit(fit, one$x, one$y, name)
    )
    cat("recorded", name, "\n")
  }
  list(
    mode = info$mode,
    version = as.character(utils::packageVersion("dbarts")),
    library = dirname(system.file(package = "dbarts")),
    records = records
  )
}

# --- comparing ---------------------------------------------------------------

compareAll <- function(a, b) {
  flat <- function(record) {
    out <- list()
    for (scenario in names(record$records)) {
      for (element in names(record$records[[scenario]])) {
        out[[paste(scenario, element, sep = "/")]] <-
          record$records[[scenario]][[element]]
      }
    }
    out
  }
  left <- flat(a)
  right <- flat(b)
  labels <- union(names(left), names(right))
  different <- character()
  for (label in labels) {
    if (!label %in% names(left) || !label %in% names(right)) {
      different <- c(different, label)
    } else if (!identical(left[[label]], right[[label]])) {
      different <- c(different, label)
    }
  }
  cat(
    "scenarios: ",
    length(a$records),
    " / ",
    length(b$records),
    "\n",
    "elements compared: ",
    length(labels),
    "\n",
    "identical: ",
    length(labels) - length(different),
    "\n",
    "different: ",
    length(different),
    "\n",
    sep = ""
  )
  if (length(different) > 0L) {
    cat(paste0("  ", different, "\n"), sep = "")
    return(1L)
  }
  0L
}

if (length(args) >= 2L && args[1L] == "record") {
  suppressPackageStartupMessages(library(dbarts))
  saveRDS(recordAll(), args[2L])
} else if (length(args) >= 3L && args[1L] == "compare") {
  status <- compareAll(readRDS(args[2L]), readRDS(args[3L]))
  quit(status = status)
} else {
  stop("usage: record <out.rds> | compare <a.rds> <b.rds>")
}
