# Cutpoint-sampler mixing study for cumulative-probit BART (design evidence).
#
# NO TREES. The mean is a linear intercept+slope drawn by a conjugate Gibbs
# step, so this measures the cutpoint update in isolation. A flexible BART mean
# competes with the cutpoints for the location far more aggressively than a
# 2-parameter linear mean, so these ESS numbers are an OPTIMISTIC bound on the
# in-sampler mixing, not a BART prediction. The RANKING across samplers is the
# transferable finding; the absolute ESS is not.
#
# Identification: sigma = 1 fixed, gamma_1 = 0 fixed (scheme A). The mean carries
# an intercept, which is what makes the location float and the AC pathology bite.
#
# Three cutpoint updates compared:
#   AC     Albert-Chib (1993) Gibbs: gamma_k ~ Uniform(order-stat window).
#   COWLES Cowles (1996) one-at-a-time truncated-normal RW-MH, latents
#          marginalized in the acceptance ratio, then z redrawn.
#   LOGGAP joint RW-MH on the unconstrained log-gaps delta = log(diff(gamma)),
#          marginal acceptance + Jacobian, then z redrawn.

set.seed(20260718L)

## ---- effective sample size: N / (1 + 2 sum rho_t), initial-positive-sequence
ess <- function(x) {
  x <- x[is.finite(x)]; x <- x - mean(x); n <- length(x)
  v0 <- sum(x * x) / n
  if (!is.finite(v0) || v0 <= 0) return(NA_real_)
  s <- 0
  for (t in 1:(n - 1L)) {
    rho <- sum(x[1:(n - t)] * x[(1 + t):n]) / (n * v0)
    if (rho <= 0) break
    s <- s + rho
  }
  n / (1 + 2 * s)
}

## ---- vectorized doubly-truncated normal via inverse CDF (unit sd). A pure
## inverse-CDF draw underflows for intervals far in the tail (Phi(lo)==Phi(hi));
## those rows fall back to the interval midpoint. A production engine would use
## Robert (1995) rejection there (a component to add C-side, section 5).
rtnorm <- function(mu, lo, hi) {
  plo <- pnorm(lo - mu); phi <- pnorm(hi - mu)
  u <- runif(length(mu))
  z <- mu + qnorm(plo + u * (phi - plo))
  flo <- ifelse(is.finite(lo), lo, mu - 8); fhi <- ifelse(is.finite(hi), hi, mu + 8)
  bad <- !is.finite(z); z[bad] <- 0.5 * (flo[bad] + fhi[bad])
  pmin(pmax(z, flo), fhi)
}

## ---- log marginal category likelihood sum over a subset of rows
loglik_rows <- function(idx, y, mu, g) {   # g = full cutpoints incl -Inf,0,..,Inf
  sum(log(pmax(pnorm(g[y[idx] + 1L] - mu[idx]) -
               pnorm(g[y[idx]] - mu[idx]), 1e-300)))
}

## ---- one data set: ordered probit, sigma = 1, true gamma, intercept + slope
simulate <- function(n, K, seed) {
  set.seed(seed)
  x <- rnorm(n)
  b0 <- 0.3; b1 <- 1.0
  gtrue <- c(-Inf, 0, cumsum(rep(1.1, K - 2L)), Inf)  # gamma_1 = 0, gaps 1.1
  z <- b0 + b1 * x + rnorm(n)
  y <- findInterval(z, gtrue[2:K]) + 1L               # 1..K
  list(x = x, y = y, K = K)
}

## ---- shared Gibbs shell; `cut_update` swaps the cutpoint sampler
run <- function(dat, method, iters = 6000L, burn = 1500L) {
  x <- dat$x; y <- dat$y; K <- dat$K; n <- length(y)
  X <- cbind(1, x); XtX <- crossprod(X); prior <- diag(c(1e-2, 1e-2))
  Vbeta <- solve(XtX + prior); Lb <- chol(Vbeta)
  g <- c(-Inf, 0, seq_len(K - 2L), Inf)               # start gaps = 1
  beta <- c(0, 1)
  free <- 3:K                                         # indices into g (gamma_2..)
  tau <- rep(0.3, K - 2L)                             # RW scales (adapted in burn)
  acc <- numeric(K - 2L)                              # rolling accept count (adaptation)
  accPost <- numeric(K - 2L)                          # post-burn accept count (report)
  keep <- matrix(NA_real_, iters - burn, K - 2L)
  z <- rtnorm(as.numeric(X %*% beta), g[y], g[y + 1L])
  for (it in seq_len(iters)) {
    mu <- as.numeric(X %*% beta)
    ## latents (all methods draw these; COWLES/LOGGAP redraw after cutpoints too)
    z <- rtnorm(mu, g[y], g[y + 1L])
    ## beta | z  (conjugate)
    bm <- Vbeta %*% crossprod(X, z)
    beta <- as.numeric(bm + t(Lb) %*% rnorm(2))
    mu <- as.numeric(X %*% beta)
    ## cutpoints
    if (method == "AC") {
      for (k in free) {
        lo <- max(g[k - 1L], if (any(y == k - 1L)) max(z[y == k - 1L]) else -Inf)
        hi <- min(g[k + 1L], if (any(y == k))     min(z[y == k])     else  Inf)
        if (is.finite(lo) && is.finite(hi) && hi > lo) # sparse boundary cell:
          g[k] <- runif(1, lo, hi)                     # one-sided window -> skip
      }
    } else if (method == "COWLES") {
      for (j in seq_along(free)) {
        k <- free[j]
        prop <- g[k] + tau[j] * rnorm(1)
        if (prop <= g[k - 1L] || prop >= g[k + 1L]) next
        idx <- which(y == k - 1L | y == k)
        gp <- g; gp[k] <- prop
        la <- loglik_rows(idx, y, mu, gp) - loglik_rows(idx, y, mu, g)
        if (is.finite(la) && log(runif(1)) < la) {
          g[k] <- prop; acc[j] <- acc[j] + 1; if (it > burn) accPost[j] <- accPost[j] + 1
        }
      }
      z <- rtnorm(mu, g[y], g[y + 1L])
    } else if (method == "LOGGAP") {
      d <- diff(g[2:K]); delta <- log(d)               # K-2 gaps (gamma_1..gamma_{K-1})
      step <- tau * rnorm(K - 2L)
      dp <- exp(delta + step); gp <- g; gp[2:K] <- c(0, cumsum(dp))
      la <- loglik_rows(seq_len(n), y, mu, gp) - loglik_rows(seq_len(n), y, mu, g) +
            sum(step)                                   # Jacobian: sum log(d'/d)
      if (is.finite(la) && log(runif(1)) < la) {
        g <- gp; acc <- acc + 1; if (it > burn) accPost <- accPost + 1
      }
      z <- rtnorm(mu, g[y], g[y + 1L])
    }
    ## adapt RW scale in burn-in toward ~30% acceptance (COWLES/LOGGAP)
    if (it <= burn && method != "AC" && it %% 100 == 0) {
      rate <- if (method == "COWLES") acc / 100 else rep(mean(acc) / 100, K - 2L)
      tau <- tau * exp(pmax(pmin(rate - 0.3, 0.5), -0.5))
      acc <- numeric(K - 2L)
    }
    if (it > burn) keep[it - burn, ] <- g[free]
  }
  list(keep = keep, ess = apply(keep, 2, ess),
       accept = if (method == "AC") NA else accPost / (iters - burn))
}

## ---- grid
scen <- expand.grid(n = c(200L, 2000L), K = c(3L, 5L))
methods <- c("AC", "COWLES", "LOGGAP")
cat(sprintf("%-6s %-4s %-7s %-8s %-8s %-8s\n",
            "n", "K", "method", "essMin", "essMean", "accept"))
for (r in seq_len(nrow(scen))) {
  dat <- simulate(scen$n[r], scen$K[r], seed = 1000L + r)
  for (m in methods) {
    t0 <- proc.time()[3]
    out <- run(dat, m)
    el <- proc.time()[3] - t0
    cat(sprintf("%-6d %-4d %-7s %-8.0f %-8.0f %-8s  (%.1fs)\n",
                scen$n[r], scen$K[r], m,
                min(out$ess), mean(out$ess),
                if (all(is.na(out$accept))) "-" else
                  paste0(sprintf("%.2f", mean(out$accept))), el))
  }
}
