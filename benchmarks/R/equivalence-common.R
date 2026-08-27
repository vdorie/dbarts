# Shared helpers for bcf-equivalence.R and multinomial-equivalence.R: the
# draws-axis reductions, the cross-host tolerance/ESS math, and the two-tier
# cross-host verdict reporter. Each caller supplies its own channel taxonomy
# (statChannels/snapshotChannels), tolerance, and report column width; every
# function here reads only the two arms' recorded channels, so a defect
# degrades both arms of a comparison symmetrically rather than laundering a
# bug into the verdict.

# A channel's draws laid out one row per cell, columns the trailing (draws)
# axis.
drawMatrix <- function(a) {
  d <- dim(a)
  matrix(a, ncol = if (is.null(d)) length(a) else d[length(d)])
}

# Posterior mean/var over a channel's trailing (draws) axis, one cell per
# leading-index combination.
drawSummary <- function(a) {
  m <- drawMatrix(a)
  list(mean = rowMeans(m), var = apply(m, 1L, var), n = ncol(m))
}

# Per-cell lag-1 autocorrelation over the draws axis, 0 where a cell is
# constant.
lag1Acf <- function(m) {
  n <- ncol(m)
  if (n < 3L) {
    return(rep(0, nrow(m)))
  }
  cx <- m[, -n, drop = FALSE]
  cy <- m[, -1L, drop = FALSE]
  cx <- cx - rowMeans(cx)
  cy <- cy - rowMeans(cy)
  den <- sqrt(rowSums(cx^2) * rowSums(cy^2))
  ifelse(den > 0, rowSums(cx * cy) / den, 0)
}

# Tier 1 on a continuous channel: per-cell |a - b| against atol + rtol * |a|,
# atol taken from the channel's own scale in this scenario. Returns the ratio
# per cell, so <= 1 everywhere is the pass. A cell that agrees exactly scores 0
# even where the tolerance is 0; a non-finite difference scores Inf so it can
# never be read as a pass.
toleranceRatio <- function(x, y, rtol) {
  d <- abs(x - y)
  ratio <- ifelse(d == 0, 0, d / (rtol * max(abs(x)) + rtol * abs(x)))
  ratio[is.na(ratio)] <- Inf
  ratio
}

# Welch z per cell with an ESS-aware denominator. The draws axis here is ONE
# autocorrelated chain, not independent seeds, so the nominal n understates
# the Monte Carlo error; the effective size comes from the pooled lag-1
# autocorrelation, floored at 2 and capped at n so the denominator can only
# grow relative to the nominal statistic.
essCompare <- function(x, y) {
  mx <- drawMatrix(x)
  my <- drawMatrix(y)
  r <- 0.5 * (lag1Acf(mx) + lag1Acf(my))
  effective <- function(m) pmin(pmax(ncol(m) * (1 - r) / (1 + r), 2), ncol(m))
  ex <- effective(mx)
  ey <- effective(my)
  list(
    z = (rowMeans(mx) - rowMeans(my)) /
      sqrt(apply(mx, 1L, var) / ex + apply(my, 1L, var) / ey),
    ess = pmin(ex, ey)
  )
}

# Per-channel non-finite counts, zeros included so the line names every gated
# channel and a reader can see which one carries them.
nonFiniteParts <- function(gated, nf) {
  vapply(
    seq_along(gated),
    function(i) {
      n <- sum(nf[, i])
      if (n == 0) {
        return(sprintf("%s 0", gated[i]))
      }
      sprintf(
        "%s %.0f cells (%s)",
        gated[i],
        n,
        if (nf[1L, i] > 0 && nf[2L, i] > 0) {
          "baseline and this run"
        } else if (nf[1L, i] > 0) {
          "baseline"
        } else {
          "this run"
        }
      )
    },
    ""
  )
}

# One scenario's cross-host verdict, three name-prefixed lines: the exempt
# roster with a NON-GATING deviation diagnostic, the non-finite count, and the
# tier verdict. Tier 1 is the gate; tier 2 runs only after tier 1 fails and
# only adjudicates - did the posterior move, or did the stream merely decouple.
# Its bar is weak by construction (the line reports how weak), so a tier-2
# pass is evidence the failure is not gross, never evidence the builds agree.
# Combinatorial channels are integer split counts: tier 1 compares them
# exactly, and tier 2 reports them diagnostically rather than through a z,
# because a Welch z on counts is a bitwise test in disguise and a
# structurally unsplit cell makes it 0/0. `snapshotChannels` and
# `crossHostRtol` are the caller's own taxonomy and tolerance; `width` is the
# caller's report column width.
compareCrossHost <- function(
  name,
  a,
  b,
  snapshotChannels,
  crossHostRtol,
  width
) {
  exempt <- intersect(snapshotChannels, names(a))
  gated <- setdiff(names(a), c("summaries", exempt))
  for (ch in gated) {
    stopifnot(
      identical(dim(a[[ch]]), dim(b[[ch]])),
      length(a[[ch]]) == length(b[[ch]])
    )
  }

  n.differ <- 0L
  rel.dev <- 0
  for (ch in exempt) {
    if (!identical(a[[ch]], b[[ch]])) {
      n.differ <- n.differ + 1L
    }
    av <- unlist(a[[ch]])
    bv <- unlist(b[[ch]])
    if (is.double(av) && length(av) == length(bv) && max(abs(av)) > 0) {
      rel.dev <- max(rel.dev, max(abs(av - bv)) / max(abs(av)))
    }
  }
  cat(sprintf(
    "%-*s exempt (cross-host): %s - %d snapshot channels skipped by design [%d of %d differ, max rel dev %.1e]\n",
    width,
    name,
    paste(exempt, collapse = ", "),
    length(exempt),
    n.differ,
    length(exempt),
    rel.dev
  ))

  # Exempting the point-in-time channels deletes the only guard that saw a
  # NaN/Inf draw: every z built from one is NaN, and a max over NaNs with
  # na.rm reports -Inf and passes. Count them directly instead.
  nf <- vapply(
    gated,
    function(ch) c(sum(!is.finite(a[[ch]])), sum(!is.finite(b[[ch]]))),
    numeric(2L)
  )
  failed <- FALSE
  if (sum(nf) == 0) {
    cat(sprintf("%-*s NON-FINITE: none\n", width, name))
  } else {
    failed <- TRUE
    cat(sprintf(
      "%-*s NON-FINITE: %s <- FAIL\n",
      width,
      name,
      paste(nonFiniteParts(gated, nf), collapse = ", ")
    ))
  }

  is.continuous <- vapply(gated, function(ch) is.double(a[[ch]]), logical(1L))
  cont <- gated[is.continuous]
  comb <- gated[!is.continuous]
  ratios <- lapply(cont, function(ch) {
    toleranceRatio(a[[ch]], b[[ch]], crossHostRtol)
  })
  names(ratios) <- cont
  worst <- vapply(ratios, max, numeric(1L))
  cell <- vapply(ratios, which.max, integer(1L))
  comb.differ <- vapply(comb, function(ch) sum(a[[ch]] != b[[ch]]), numeric(1L))
  comb.n <- vapply(comb, function(ch) length(a[[ch]]), numeric(1L))

  bad.cont <- worst > 1
  bad.comb <- comb.differ > 0
  if (!any(bad.cont) && !any(bad.comb)) {
    cat(sprintf(
      "%-*s tier 1 PASS: max dev ratio %.1e (%s)%s%s\n",
      width,
      name,
      max(c(0, worst)),
      paste(sprintf("%s %.1e", cont, worst), collapse = ", "),
      if (length(comb) > 0L) {
        paste0(", ", paste(sprintf("%s identical", comb), collapse = ", "))
      } else {
        ""
      },
      if (all(worst == 0)) {
        sprintf(" - all %d GATED channels bitwise identical", length(gated))
      } else {
        ""
      }
    ))
    return(list(fail = failed, decoupled = FALSE))
  }

  cat(sprintf(
    "%-*s tier 1 FAIL: %s\n",
    width,
    name,
    paste(
      c(
        sprintf(
          "%s ratio %.1e (cell %d)",
          cont[bad.cont],
          worst[bad.cont],
          cell[bad.cont]
        ),
        sprintf(
          "%s %.0f of %.0f cells differ",
          comb[bad.comb],
          comb.differ[bad.comb],
          comb.n[bad.comb]
        )
      ),
      collapse = ", "
    )
  ))

  tier2 <- lapply(cont, function(ch) essCompare(a[[ch]], b[[ch]]))
  names(tier2) <- cont
  z <- unlist(lapply(tier2, function(s) s$z))
  ess <- unlist(lapply(tier2, function(s) s$ess))
  n.warn <- sum(abs(z) > 3, na.rm = TRUE)
  n.fail <- sum(abs(z) > 4, na.rm = TRUE)
  cat(sprintf(
    "%-*s decoupled: statistical - %d summaries, ESS-adjusted max |z| = %.2f (weak bar: tolerates %.2f posterior sd)%s%s%s\n",
    width,
    name,
    length(z),
    max(abs(z), na.rm = TRUE),
    4 * sqrt(2 / median(ess, na.rm = TRUE)),
    paste0(
      c(
        "",
        sprintf(
          "%s diagnostic-only, %.0f of %.0f cells differ, worst %.0f",
          comb[bad.comb],
          comb.differ[bad.comb],
          comb.n[bad.comb],
          vapply(
            comb[bad.comb],
            function(ch) max(abs(as.numeric(a[[ch]]) - as.numeric(b[[ch]]))),
            numeric(1L)
          )
        )
      ),
      collapse = "; "
    ),
    if (n.warn > 0L) sprintf(", %d with |z| > 3", n.warn) else "",
    if (n.fail > 0L) sprintf(", %d with |z| > 4 <- FAIL", n.fail) else ""
  ))
  if (n.fail > 0L) {
    failed <- TRUE
    offenders <- which(abs(z) > 4)
    cat(
      "  worst offenders:",
      paste0(
        names(z)[offenders],
        " (z=",
        round(z[offenders], 2L),
        ")",
        collapse = ", "
      ),
      "\n"
    )
  }
  list(fail = failed, decoupled = TRUE)
}
