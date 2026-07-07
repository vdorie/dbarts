# End-to-end sweep-time cross-check for benchmarks/kernels/linear_leaf.cpp: the
# per-sweep wall time of a linear-leaves dbarts fit, to confirm the C++ shim's
# sweep denominators. Runs against the installed package (R CMD INSTALL . first).

suppressMessages(library(dbarts))
set.seed(1)

timePerSweep <- function(n, q, ntree = 50L, nburn = 20L, ndraw = 40L) {
  X <- matrix(runif(n * q), n, q)
  colnames(X) <- paste0("x", seq_len(q))
  f <- 2 * X[, 1] + 1.5 * X[, 2] + ifelse(X[, 3] > 0.5, 1, -1)
  y <- f + 0.5 * (runif(n) - 0.5)
  ctrl <- dbartsControl(
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE,
    verbose = FALSE
  )
  s <- dbarts(X, y, control = ctrl, node.prior = linear(seq_len(q), k = 2))
  s$run(nburn, 0L) # warm up
  el <- system.time(s$run(0L, ndraw))[["elapsed"]]
  el / ndraw * 1e6
}

for (cfg in list(c(1e4, 8), c(1e5, 4), c(1e5, 8))) {
  n <- as.integer(cfg[1])
  q <- as.integer(cfg[2])
  cat(sprintf("n=%d q=%d: %.0f us/sweep\n", n, q, timePerSweep(n, q)))
}
