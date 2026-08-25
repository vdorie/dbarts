#!/usr/bin/env Rscript
# Post-hoc band computation over every finished arm: the per-functional 5%
# ecdf band (the level every non-matrix recorded result was read at) and the
# band Bonferroni'd over THIS matrix's total functional count (the family-tiers
# admission convention). Prints one row per functional and the flag detail.
src <- Sys.getenv("SBC_SRC")
out <- Sys.getenv("SBC_OUT")
source(file.path(src, "benchmarks", "R", "sbc.R"))
files <- sort(Sys.glob(file.path(out, "fit-*.rds")))
fits <- lapply(files, readRDS)
names(fits) <- sub("^fit-(.*)\\.rds$", "\\1", basename(files))
M <- sum(vapply(fits, function(f) ncol(f$ranks), integer(1)))
alphaB <- 0.05 / M
cat(sprintf("arms %d  total functionals %d  Bonferroni alpha %.6f\n\n", length(fits), M, alphaB))
cat(sprintf("%-16s %-12s %8s %8s %9s %8s %8s %7s %6s\n", "arm", "functional",
            "chisqP", "ksP", "ecdfDiff", "band05", "bandBonf", "meanRk", "verd"))
flagged <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  for (fn in colnames(f$ranks)) {
    rk <- f$ranks[, fn]
    u5 <- rankUniformity(rk, f$L, alpha = 0.05)
    ub <- rankUniformity(rk, f$L, alpha = alphaB)
    v <- if (u5$ecdfDiff > ub$ecdfBand) "FLAG" else if (u5$ecdfDiff > u5$ecdfBand) "NOTE" else "PASS"
    cat(sprintf("%-16s %-12s %8.3f %8.3f %9.4f %8.4f %8.4f %7.1f %6s\n",
                nm, fn, u5$chisqP, u5$ksP, u5$ecdfDiff, u5$ecdfBand, ub$ecdfBand, u5$mean, v))
    if (v != "PASS") flagged[[paste(nm, fn)]] <- list(nm = nm, fn = fn, u = u5, ub = ub, L = f$L, rk = rk)
  }
}
cat("\n== flag detail ==\n")
for (k in names(flagged)) {
  z <- flagged[[k]]
  cat(sprintf("\n[%s] ecdfDiff %.4f  band05 %.4f  bandBonf %.4f  chisqP %.3f  mean rank %.1f (target %.1f)\n",
              k, z$u$ecdfDiff, z$u$ecdfBand, z$ub$ecdfBand, z$u$chisqP, z$u$mean, z$L / 2))
  cat("bin counts (20 bins, low rank -> high):\n  ", paste(z$u$counts, collapse = " "), "\n")
}
if (length(flagged) == 0L) cat("(none)\n")
